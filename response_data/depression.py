"""离散栅格洼地分析：最小溢流高程、填洼体积与规范合并树。

口径（见 docs/data-contract.md）：
- 有限矩形网格，四邻接（上下左右）；NoData(None) 是不可穿越屏障，不是零高程。
- 出口仅由配置给出且必须是有效格；边界格不自动成为出口。
- 对可排格 c，溢流高程 H(c) = min_{c→出口路径 p} max_{格 x∈p} z(x)，
  即格 c 越过地形连到任一真实出口所需抬到的最低高程；填洼深度 H(c)-z(c)。
- 体积按每格独立正面积加权：area(c)*max(0, H(c)-z(c))；不可排格不赋 H，不参与体积。
- 合并树按原高程分层批处理：同高全部格同批激活，批内连接统一合并，
  每个含新格的批后分量产生一个多叉节点；不产生逐格顺序导致的等高二叉链。
本模块只做离散地形连通分析，不模拟降雨，不把填洼体积当作洪水蓄量预测。
"""
from fractions import Fraction
import heapq


def _grid(v):
    """解析并校验网格，返回 (z, area, outlet_idx, height, width)，值为精确 Fraction。"""
    z = v['elevation']
    area = v['area']
    if not z or not z[0] or any(len(r) != len(z[0]) for r in z):
        raise ValueError('网格形状无效')
    if len(area) != len(z) or any(len(r) != len(z[0]) for r in area):
        raise ValueError('面积形状不符')
    height, width = len(z), len(z[0])
    max_cells = v.get('max_cells')
    if max_cells is not None:
        if isinstance(max_cells, bool) or not isinstance(max_cells, int) or max_cells <= 0:
            raise ValueError('max_cells 应为正整数')
        if height * width > max_cells:
            raise ValueError('网格格数超过 max_cells')
    zz, aa = [], []
    for zrow, arow in zip(z, area):
        zr, ar = [], []
        for x, y in zip(zrow, arow):
            if isinstance(y, bool) or not isinstance(y, (str, int)):
                raise ValueError('格面积应正')
            a = Fraction(y)
            if a <= 0:
                raise ValueError('格面积应正')
            ar.append(a)
            if x is None:
                zr.append(None)
            elif isinstance(x, bool) or not isinstance(x, (str, int)):
                raise ValueError('精确高程使用整数或分数字符串')
            else:
                zr.append(Fraction(x))
        zz.append(tuple(zr))
        aa.append(tuple(ar))
    outlet_idx, seen = [], set()
    for pair in v['outlets']:
        i, j = pair
        if not 0 <= i < height or not 0 <= j < width or zz[i][j] is None:
            raise ValueError('出口不是有效格')
        idx = i * width + j
        if idx in seen:
            raise ValueError('出口重复')
        seen.add(idx)
        outlet_idx.append(idx)
    if not any(x is not None for row in zz for x in row):
        raise ValueError('网格没有有效格')
    return tuple(zz), tuple(aa), outlet_idx, height, width


def _neighbors(idx, width, height):
    r, c = divmod(idx, width)
    if r > 0:
        yield idx - width
    if r + 1 < height:
        yield idx + width
    if c > 0:
        yield idx - 1
    if c + 1 < width:
        yield idx + 1


def _node_id(layer, new_cells, width):
    """规范稳定 ID：层高规范分数 + 排序的新格行列坐标。"""
    coord = ','.join(f'{r}.{c}' for r, c in (divmod(i, width) for i in new_cells))
    return f'z={str(layer)}|{coord}'


def build_merge_tree(v, cell_order=None, edge_order=None):
    """按原高程从低到高批处理激活，构造规范多叉合并树并给出每格 H。

    cell_order / edge_order 仅用于测试顺序无关性；生产口径恒用规范序
    （同层格按索引升序、邻接按固定方向）。
    返回内部结构（Fraction 精确值），不可排格 H 为 None。
    """
    z, area, outlets, height, width = _grid(v)
    n = height * width
    zflat = tuple(z[i][j] for i in range(height) for j in range(width))
    is_outlet = [False] * n
    for idx in outlets:
        is_outlet[idx] = True
    valid = [idx for idx in range(n) if zflat[idx] is not None]

    if cell_order is None:
        groups = {}
        for idx in valid:
            groups.setdefault(zflat[idx], []).append(idx)
    else:
        groups = {}
        for idx in cell_order:  # 测试用：桶内顺序由参数决定，层高仍按数值排序
            if zflat[idx] is not None:
                groups.setdefault(zflat[idx], []).append(idx)
    layers = sorted(groups)

    parent = list(range(n))
    rank = [0] * n

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb:
            return ra
        if rank[ra] < rank[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        if rank[ra] == rank[rb]:
            rank[ra] += 1
        return ra

    root_node = {}        # DSU 根 -> 已归属的树节点 id（仅批前分量根持有）
    nodes = []            # 节点按创建（层高升序）排列
    node_of_cell = [None] * n
    id_to_node_pending = {}  # 创建期即可按 id 查节点（孩子总是更早创建）

    for layer in layers:
        born = groups[layer]
        # 阶段一：本批全部同高格同时激活，各自成 DSU 分量。
        scratch = {}      # 根 -> {'new': set(本批新格), 'old': set(吸收的批前节点 id)}
        for idx in born:
            scratch[idx] = {'new': {idx}, 'old': set()}
        # 阶段二：收集并统一处理全部相关连接（新-新平台边、新-旧跨层边），无逐格子批。
        edges = []
        for idx in born:
            for nb in _neighbors(idx, width, height):
                if zflat[nb] is None:
                    continue  # NoData 屏障：不可穿越，绝不进入连接
                if zflat[nb] > layer:
                    continue  # 更高格尚未激活
                edges.append((idx, nb))
        if edge_order is not None:
            edges = edge_order(edges)
        for a, b in edges:
            ra, rb = find(a), find(b)
            if ra == rb:
                continue
            sa, sb = scratch.get(ra), scratch.get(rb)
            na, nb_old = root_node.get(ra), root_node.get(rb)
            rnew = union(ra, rb)
            bag = {'new': set(), 'old': set()}
            for s in (sa, sb):
                if s is not None:
                    bag['new'] |= s['new']
                    bag['old'] |= s['old']
            for old in (na, nb_old):
                if old is not None:
                    bag['old'].add(old)
            scratch.pop(ra, None)
            scratch.pop(rb, None)
            root_node.pop(ra, None)
            root_node.pop(rb, None)
            scratch[rnew] = bag
        # 阶段三：每个含新格的批后分量恰生成一个多叉节点；
        # 同层多个节点按规范 node_id（即排序成员坐标）排序，与 DSU 根顺序无关。
        pending_nodes = []
        roots_here = sorted({find(idx) for idx in born})
        for r in roots_here:
            bag = scratch.get(r)
            members = sorted(bag['new']) if bag else []
            child_ids = sorted(bag['old']) if bag else []
            nid = _node_id(layer, members, width)
            if nid in id_to_node_pending:
                raise RuntimeError('节点 ID 不唯一')
            # 含出口：本批新格自带出口，或吸收的任一较低分量已含出口。
            has_outlet = (any(is_outlet[i] for i in members)
                          or any(id_to_node_pending[cid]['has_outlet']
                                 for cid in child_ids))
            pending_nodes.append((nid, r, {
                'node_id': nid,
                'layer': layer,
                'new_cells': members,
                'child_ids': child_ids,
                'has_outlet': has_outlet,
                'parent_id': None,
                'subtree_cell_count': len(members),
                'subtree_min_elevation': layer,
                'first_outlet_layer': None,
            }))
        for nid, r, node in sorted(pending_nodes, key=lambda t: t[0]):
            nodes.append(node)
            id_to_node_pending[nid] = node
            for i in node['new_cells']:
                node_of_cell[i] = nid
            root_node[r] = nid

    id_to_node = id_to_node_pending
    # 回填父子关系与子树统计（节点按层高升序，逆序即后序）。
    for nd in reversed(nodes):
        for cid in nd['child_ids']:
            child = id_to_node[cid]
            child['parent_id'] = nd['node_id']
            nd['subtree_cell_count'] += child['subtree_cell_count']
            nd['subtree_min_elevation'] = min(
                nd['subtree_min_elevation'], child['subtree_min_elevation'])
    # 节点溢流层：本节点分量在本层已含出口取本层；否则继承父节点。
    # 父节点总是更晚创建，故按创建顺序逆序（根→叶）传递。
    # 格 c 的 H 即其出生节点的溢流层，也就是“含 c 的分量首次含出口”的层高；
    # 链上无含出口祖先则不可排。
    for nd in reversed(nodes):
        if nd['has_outlet']:
            nd['first_outlet_layer'] = nd['layer']
        elif nd['parent_id'] is None:
            nd['first_outlet_layer'] = None
        else:
            nd['first_outlet_layer'] = id_to_node[nd['parent_id']]['first_outlet_layer']
    # 每格 H：出生节点首次含出口层高；不可排连通分量为 None（不伪装成数值）。
    spill = tuple(
        None if zflat[i] is None
        else id_to_node[node_of_cell[i]]['first_outlet_layer']
        for i in range(n))
    return {
        'z': zflat,
        'area': tuple(area[i][j] for i in range(height) for j in range(width)),
        'is_outlet': tuple(is_outlet),
        'outlets': tuple(outlets),
        'height': height,
        'width': width,
        'nodes': nodes,
        'node_of_cell': tuple(node_of_cell),
        'spill': spill,
    }


def minimax_heap(v):
    """有界优先队列求 H：标签 = 路径上最大原高程的最小可行值。

    max 半群上的 Dijkstra：严格更优才松弛，节点按 (标签, 索引) 最终化，
    先定型的同标签父被保留，因此父森林无环、确定且不依赖输入顺序。
    返回 (label, parent)。
    """
    z, _area, outlets, height, width = _grid(v)
    n = height * width
    zflat = tuple(z[i][j] for i in range(height) for j in range(width))
    label = [None] * n
    parent = [None] * n
    finalized = [False] * n
    heap = []
    for idx in outlets:
        label[idx] = zflat[idx]
        parent[idx] = idx  # 出口以自身为父，标识父森林的根
        heapq.heappush(heap, (zflat[idx], idx))
    while heap:
        d, u = heapq.heappop(heap)
        if finalized[u]:
            continue
        finalized[u] = True
        for w in _neighbors(u, width, height):
            if zflat[w] is None or finalized[w]:
                continue
            nd = d if d >= zflat[w] else zflat[w]  # max(d, z(w))
            old = label[w]
            # 严格更优才松弛；同标签先定型者保留，父结构确定且无环。
            if old is None or nd < old:
                label[w] = nd
                parent[w] = u
                heapq.heappush(heap, (nd, w))
    return label, parent


def witness_path(idx, parent, is_outlet, check_acyclic=True):
    """沿父结构还原见证路径；调用方另行校验路径峰值恰为 H(idx)。

    check_acyclic=False 时信任父森林（其无环性由测试与最终化规则保证），
    供大图批量重建使用。
    """
    if parent[idx] is None:
        return None
    path, cur = [], idx
    seen = set() if check_acyclic else None
    while True:
        if check_acyclic:
            if cur in seen:
                raise RuntimeError('父结构出现环')
            seen.add(cur)
        path.append(cur)
        if is_outlet[cur]:
            return path
        nxt = parent[cur]
        if nxt is None or nxt == cur:
            raise RuntimeError('父结构未到达真实出口')
        cur = nxt


def analyze(v):
    """完整离散洼地分析；DSU 批处理与堆式最劣路径独立互校后返回精确结果。"""
    tree = build_merge_tree(v)
    label, parent = minimax_heap(v)
    n = tree['height'] * tree['width']
    zflat, area, is_outlet = tree['z'], tree['area'], tree['is_outlet']
    # 互校一：两套独立算法的 H 逐格一致（含不可排判定）。
    for i in range(n):
        if tree['spill'][i] != label[i]:
            raise RuntimeError(f'溢流高程不一致 @ {i}')
    id_to_node = {nd['node_id']: nd for nd in tree['nodes']}
    # 互校二：每格出生节点的首次含出口层高与 H 相同。
    for i in range(n):
        if zflat[i] is None:
            continue
        fol = id_to_node[tree['node_of_cell'][i]]['first_outlet_layer']
        if fol != tree['spill'][i]:
            raise RuntimeError('合并树首接出口层高与 H 不一致')
    cells = []
    total_volume = Fraction(0)
    drainable_area = Fraction(0)
    for i in range(n):
        if zflat[i] is None:
            continue
        h = tree['spill'][i]
        depth = None if h is None else h - zflat[i]
        vol = Fraction(0) if h is None else area[i] * (h - zflat[i])
        if h is not None:
            total_volume += vol
            drainable_area += area[i]
        cells.append({
            'index': i,
            'row': i // tree['width'],
            'col': i % tree['width'],
            'elevation': zflat[i],
            'is_outlet': bool(is_outlet[i]),
            'drainable': h is not None,
            'spill_elevation': h,
            'fill_depth': depth,
            'fill_volume': vol,
            'node_id': tree['node_of_cell'][i],
            'parent_index': parent[i] if h is not None else None,
        })
    # 以最终树根节点标识物理连通分量；不可排分量显式标注，不赋任何数值高程。
    root_of_node = {}
    for nd in tree['nodes']:
        cur = nd['node_id']
        while id_to_node[cur]['parent_id'] is not None:
            cur = id_to_node[cur]['parent_id']
        root_of_node[nd['node_id']] = cur
    components = {}
    for rec in cells:
        rid = root_of_node[rec['node_id']]
        c = components.setdefault(rid, {
            'root_node_id': rid,
            'cell_count': 0,
            'min_elevation': rec['elevation'],
            'drainable': rec['drainable'],
            'contains_outlet': False,
        })
        c['cell_count'] += 1
        c['min_elevation'] = min(c['min_elevation'], rec['elevation'])
        c['contains_outlet'] = c['contains_outlet'] or rec['is_outlet']
    return {
        'height': tree['height'],
        'width': tree['width'],
        'cells': cells,
        'nodes': tree['nodes'],
        'components': sorted(components.values(), key=lambda c: c['root_node_id']),
        'total_fill_volume': total_volume,
        'drainable_cell_count': sum(1 for r in cells if r['drainable']),
        'undrainable_cell_count': sum(1 for r in cells if not r['drainable']),
        'drainable_area': drainable_area,
        '_parent': parent,
        '_zflat': zflat,
        '_is_outlet': is_outlet,
    }


def to_jsonable(result, with_paths=False):
    """把精确 Fraction 结果转为 JSON 安全结构（分数保留规范字符串）。"""
    parent = result['_parent']
    zflat = result['_zflat']
    is_outlet = result['_is_outlet']

    def f(x):
        return None if x is None else str(x)

    cells = []
    for rec in result['cells']:
        out = {k: rec[k] for k in ('index', 'row', 'col', 'is_outlet', 'drainable',
                                   'node_id', 'parent_index')}
        out['elevation'] = f(rec['elevation'])
        out['spill_elevation'] = f(rec['spill_elevation'])
        out['fill_depth'] = f(rec['fill_depth'])
        out['fill_volume'] = f(rec['fill_volume'])
        if with_paths and rec['drainable']:
            # 批量模式信任最终化规则保证的无环父森林；逐格峰值与终点仍逐一校验。
            path = witness_path(rec['index'], parent, is_outlet,
                                check_acyclic=False)
            peak = max(zflat[k] for k in path)
            if peak != rec['spill_elevation']:  # 见证：路径最大原高程恰为 H
                raise RuntimeError('见证路径峰值与 H 不符')
            if not is_outlet[path[-1]]:
                raise RuntimeError('见证路径终点不是真实出口')
            out['witness_path'] = path
        cells.append(out)
    nodes = [{
        'node_id': nd['node_id'],
        'layer': f(nd['layer']),
        'new_cells': list(nd['new_cells']),
        'child_ids': list(nd['child_ids']),
        'parent_id': nd['parent_id'],
        'has_outlet': nd['has_outlet'],
        'first_outlet_layer': f(nd['first_outlet_layer']),
        'subtree_cell_count': nd['subtree_cell_count'],
        'subtree_min_elevation': f(nd['subtree_min_elevation']),
    } for nd in result['nodes']]
    comps = [{**c, 'min_elevation': f(c['min_elevation'])}
             for c in result['components']]
    return {
        'height': result['height'],
        'width': result['width'],
        'method_note': ('离散地形连通分析：H 为最小最大路径高程；不模拟降雨，'
                        '填洼体积不等于洪水蓄量预测'),
        'cells': cells,
        'merge_tree': {'node_count': len(nodes), 'nodes': nodes},
        'components': comps,
        'totals': {
            'total_fill_volume': f(result['total_fill_volume']),
            'drainable_cell_count': result['drainable_cell_count'],
            'undrainable_cell_count': result['undrainable_cell_count'],
            'drainable_area': f(result['drainable_area']),
        },
    }
