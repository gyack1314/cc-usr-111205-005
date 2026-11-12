"""离散洼地分析：独立基准互校、批处理规范树与边界口径。"""
import json
import random
import unittest
from fractions import Fraction

from response_data.io import load
from response_data.depression import (
    analyze, build_merge_tree, minimax_heap, witness_path, to_jsonable, _grid,
    _neighbors)


def simple_path_oracle(v):
    """小图独立基准：DFS 穷举每条到出口的简单路径，取路径峰值的最小值。"""
    z, _area, outlets, height, width = _grid(v)
    n = height * width
    zf = [z[i][j] for i in range(height) for j in range(width)]
    adj = [[w for w in _neighbors(i, width, height) if zf[w] is not None]
           for i in range(n)]
    expected = [None] * n
    for start in range(n):
        if zf[start] is None:
            continue
        best = None
        stack = [(start, frozenset([start]), zf[start])]
        while stack:
            u, used, peak = stack.pop()
            if u in outlets:
                best = peak if best is None or peak < best else best
                continue
            for w in adj[u]:
                if w not in used:
                    stack.append((w, used | frozenset([w]), max(peak, zf[w])))
        expected[start] = best
    return expected, zf


def threshold_bfs_oracle(v):
    """阈值可达基准：对每个层高做 BFS，求每格可连通任一出口的最小阈值。"""
    z, _area, outlets, height, width = _grid(v)
    n = height * width
    zf = [z[i][j] for i in range(height) for j in range(width)]
    levels = sorted({x for x in zf if x is not None})

    def reachable(limit):
        seen = {o for o in outlets if zf[o] <= limit}
        stack = list(seen)
        while stack:
            u = stack.pop()
            for w in _neighbors(u, width, height):
                if w in seen or zf[w] is None or zf[w] > limit:
                    continue
                seen.add(w)
                stack.append(w)
        return seen

    out = [None] * n
    for level in levels:
        for idx in reachable(level):
            if out[idx] is None:
                out[idx] = level
    return out, zf


class DepressionTests(unittest.TestCase):
    def test_normal_hand_computed(self):
        v, _ = load('data/normal.json')
        r = analyze(v)
        by = {c['index']: c for c in r['cells']}
        # 两个低洼盆地：左盆地 H=3（经 (1,2) 鞍部），右盆地 H=2（出口同层）。
        self.assertEqual(by[6]['spill_elevation'], Fraction(3))
        self.assertEqual(by[11]['spill_elevation'], Fraction(3))
        self.assertEqual(by[8]['spill_elevation'], Fraction(2))
        self.assertEqual(by[13]['spill_elevation'], Fraction(2))
        self.assertEqual(by[14]['spill_elevation'], Fraction(2))  # 出口自身
        self.assertEqual(by[0]['spill_elevation'], Fraction(5))
        self.assertEqual(by[2]['spill_elevation'], Fraction(4))
        # 屏障格不出现在结果中。
        self.assertNotIn(12, by)
        # 填洼深度与独立正面积加权体积。
        self.assertEqual(by[6]['fill_depth'], Fraction(2))
        self.assertEqual(by[8]['fill_depth'], Fraction(1))
        self.assertEqual(r['total_fill_volume'], Fraction(6))
        self.assertEqual(r['undrainable_cell_count'], 0)
        self.assertEqual(r['drainable_cell_count'], 24)

    def test_simple_path_oracle_on_samples(self):
        for name in ('normal', 'no-outlet', 'single-plateau',
                     'enclosed-undrainable', 'multi-outlet-plateau'):
            with self.subTest(file=name):
                v, _ = load(f'data/{name}.json')
                expected, zf = simple_path_oracle(v)
                r = analyze(v)
                got = {c['index']: c['spill_elevation'] for c in r['cells']}
                for i, h in enumerate(expected):
                    if zf[i] is None:
                        self.assertNotIn(i, got)
                    else:
                        self.assertEqual(got[i], h)

    def test_threshold_bfs_oracle_random_grids(self):
        rng = random.Random(20261112)
        for trial in range(40):
            h = rng.randint(1, 4)
            w = rng.randint(1, 4)
            levels = [Fraction(1), Fraction('3/2'), Fraction(2), Fraction(4)]
            elev, outlets = [], []
            for i in range(h):
                row = []
                for j in range(w):
                    if rng.random() < 0.18:
                        row.append(None)
                    else:
                        row.append(rng.choice(levels))
                        if rng.random() < 0.2:
                            outlets.append([i, j])
                elev.append(row)
            if not any(x is not None for row in elev for x in row) or not outlets:
                continue
            # 保证出口落在有效格。
            outlets = [p for p in outlets if elev[p[0]][p[1]] is not None]
            if not outlets:
                continue
            v = {'elevation': [[None if x is None else str(x) for x in row]
                               for row in elev],
                 'area': [[str(rng.choice([Fraction(1), Fraction('1/3'),
                                           Fraction('5/2')])) for _ in row]
                          for row in elev],
                 'outlets': outlets}
            with self.subTest(trial=trial):
                expected, zf = threshold_bfs_oracle(v)
                r = analyze(v)
                got = {c['index']: c['spill_elevation'] for c in r['cells']}
                for i, hh in enumerate(expected):
                    if zf[i] is None:
                        continue
                    self.assertEqual(got[i], hh)
                # 体积恒等式：总体积 == 逐格面积*(H-z)，可排格无负深度。
                total = Fraction(0)
                for c in r['cells']:
                    if c['drainable']:
                        self.assertGreaterEqual(c['fill_depth'], 0)
                        total += c['fill_volume']
                    else:
                        self.assertIsNone(c['fill_depth'])
                        self.assertIsNone(c['spill_elevation'])
                self.assertEqual(total, r['total_fill_volume'])

    def test_no_outlet_all_undrainable_without_fake_numbers(self):
        v, _ = load('data/no-outlet.json')
        j = to_jsonable(analyze(v))
        self.assertTrue(all(c['spill_elevation'] is None for c in j['cells']))
        self.assertTrue(all(c['fill_depth'] is None for c in j['cells']))
        self.assertEqual(j['totals']['undrainable_cell_count'], 24)
        self.assertEqual(j['totals']['total_fill_volume'], '0')
        # 合并树仍完整记录层级，只是没有任何含出口节点。
        self.assertTrue(all(not nd['has_outlet'] for nd in j['merge_tree']['nodes']))

    def test_enclosed_component_split(self):
        v, _ = load('data/enclosed-undrainable.json')
        r = analyze(v)
        comps = r['components']
        self.assertEqual(len(comps), 2)  # 外环可排分量 + 屏障封死的内区
        drain = [c for c in comps if c['drainable']]
        blocked = [c for c in comps if not c['drainable']]
        self.assertEqual(len(drain), 1)
        self.assertEqual(len(blocked), 1)
        self.assertEqual(blocked[0]['cell_count'], 1)
        self.assertEqual(blocked[0]['min_elevation'], Fraction(1))
        self.assertFalse(blocked[0]['contains_outlet'])
        # 不可排格没有数值 H，也不贡献体积。
        for c in r['cells']:
            if not c['drainable']:
                self.assertIsNone(c['spill_elevation'])
        self.assertEqual(r['total_fill_volume'], Fraction(0))
        self.assertEqual(r['undrainable_cell_count'], 1)
        self.assertEqual(r['drainable_cell_count'], 16)

    def test_single_plateau_is_one_node(self):
        v, _ = load('data/single-plateau.json')
        j = to_jsonable(analyze(v))
        self.assertEqual(j['merge_tree']['node_count'], 1)
        node = j['merge_tree']['nodes'][0]
        self.assertEqual(node['layer'], '2')
        self.assertEqual(node['subtree_cell_count'], 25)
        self.assertEqual(node['child_ids'], [])  # 无更低孩子 → 新叶，无二叉链
        self.assertTrue(node['has_outlet'])

    def test_multi_outlet_unique_heights_nonunique_parents(self):
        v, _ = load('data/multi-outlet-plateau.json')
        r = analyze(v)
        # 多出口：每个有效格 H 仍唯一确定。
        for c in r['cells']:
            self.assertTrue(c['drainable'])
            self.assertEqual(c['spill_elevation'],
                             Fraction(1) if c['is_outlet'] else Fraction(2))
        # 分数面积：平台层 z=2 格深度为 0；两个非出口的 z=1 角格深 1、各贡献 1。
        self.assertEqual(r['total_fill_volume'], Fraction(2))
        by = {c['index']: c for c in r['cells']}
        self.assertEqual(by[4]['fill_volume'], Fraction(1))
        self.assertEqual(by[20]['fill_volume'], Fraction(1))
        # 半面积格恒零体积，验证按每格独立面积加权。
        self.assertTrue(all(by[i]['fill_volume'] == 0 for i in (6, 7, 8, 11, 12)))
        nodes = r['nodes']
        # 四个层1出口格各自成叶（对角的两个不相邻，互不连通），平台层多叉合并。
        self.assertEqual(len(nodes), 5)
        plateau = max(nodes, key=lambda nd: nd['layer'])
        self.assertEqual(plateau['layer'], Fraction(2))
        self.assertEqual(len(plateau['child_ids']), 4)
        # 父结构合法但不唯一：所有见证路径无环且终点为真实出口。
        j = to_jsonable(r, with_paths=True)
        for c in j['cells']:
            path = c['witness_path']
            self.assertEqual(len(path), len(set(path)))  # 无环
            self.assertTrue(r['_is_outlet'][path[-1]])

    def test_witness_path_peak_equals_H(self):
        v, _ = load('data/normal.json')
        r = analyze(v)
        zf, parent, outlet = r['_zflat'], r['_parent'], r['_is_outlet']
        for c in r['cells']:
            path = witness_path(c['index'], parent, outlet)
            self.assertEqual(max(zf[k] for k in path), c['spill_elevation'])
            self.assertTrue(outlet[path[-1]])
            # 路径相邻格必须四邻接（同行差 1 或同列差 width）。
            width = r['width']
            for a, b in zip(path, path[1:]):
                adjacent = abs(a - b) == 1 and min(a, b) // width == max(a, b) // width
                adjacent = adjacent or abs(a - b) == width
                self.assertTrue(adjacent)

    def test_batch_tree_no_plateau_binary_chain(self):
        """同高激活的平台必须单批合并：规范树不依赖逐格处理顺序。"""
        v, _ = load('data/single-plateau.json')
        ref = build_merge_tree(v)
        self.assertEqual(len(ref['nodes']), 1)
        rng = random.Random(7)
        valid = [i for i in range(25)]
        for _ in range(8):
            order = valid[:]
            rng.shuffle(order)
            t = build_merge_tree(v, cell_order=order,
                                 edge_order=lambda es: rng.sample(es, len(es)))
            self.assertEqual([nd['node_id'] for nd in t['nodes']],
                             [nd['node_id'] for nd in ref['nodes']])
            self.assertEqual(t['spill'], ref['spill'])
            self.assertEqual(len(t['nodes'][0]['new_cells']), 25)

    def test_order_invariance_random_grids(self):
        rng = random.Random(99)
        for trial in range(25):
            h, w = rng.randint(2, 4), rng.randint(2, 4)
            elev, outlets = [], []
            for i in range(h):
                row = []
                for j in range(w):
                    x = None if rng.random() < 0.15 else rng.randint(1, 3)
                    row.append(x)
                    if x is not None and rng.random() < 0.25:
                        outlets.append([i, j])
                elev.append(row)
            outlets = [p for p in outlets if elev[p[0]][p[1]] is not None]
            if not outlets:
                continue
            v = {'elevation': [[None if x is None else str(x) for x in row]
                               for row in elev],
                 'area': [['1'] * w for _ in range(h)],
                 'outlets': outlets}
            ref = build_merge_tree(v)
            ref_sig = [(nd['node_id'], nd['layer'], tuple(nd['child_ids']),
                        nd['has_outlet'], nd['first_outlet_layer'])
                       for nd in ref['nodes']]
            order = [i for i in range(h * w)
                     if ref['z'][i] is not None]
            rng.shuffle(order)
            t = build_merge_tree(v, cell_order=order,
                                 edge_order=lambda es: rng.sample(es, len(es)))
            sig = [(nd['node_id'], nd['layer'], tuple(nd['child_ids']),
                    nd['has_outlet'], nd['first_outlet_layer']) for nd in t['nodes']]
            self.assertEqual(sig, ref_sig)
            self.assertEqual(t['spill'], ref['spill'])

    def test_node_ids_stable_and_compact(self):
        v, _ = load('data/normal.json')
        j = to_jsonable(analyze(v))
        ids = [nd['node_id'] for nd in j['merge_tree']['nodes']]
        self.assertEqual(len(ids), len(set(ids)))  # 公开 ID 唯一
        # ID 仅由层高与本层新格坐标决定；节点只列本层新格与孩子引用，不复制全体成员。
        for nd in j['merge_tree']['nodes']:
            self.assertTrue(nd['node_id'].startswith(f"z={nd['layer']}|"))
            self.assertLessEqual(len(nd['new_cells']), nd['subtree_cell_count'])

    def test_merge_tree_saddle_layers(self):
        """两个层 1 叶洼地分别在层 2、层 3 鞍部高度接通出口侧分量。"""
        v, _ = load('data/normal.json')
        j = to_jsonable(analyze(v))
        nodes = {nd['node_id']: nd for nd in j['merge_tree']['nodes']}
        leaves = [nd for nd in nodes.values() if not nd['child_ids']]
        self.assertEqual(sorted(nd['layer'] for nd in leaves), ['1', '1'])
        saddle2 = next(nd for nd in nodes.values() if nd['layer'] == '2')
        saddle3 = next(nd for nd in nodes.values() if nd['layer'] == '3')
        self.assertTrue(saddle2['has_outlet'])
        self.assertTrue(saddle3['has_outlet'])
        # 右洼地首接出口在 2，左洼地在 3。
        self.assertEqual(nodes[saddle2['child_ids'][0]]['first_outlet_layer'], '2')
        self.assertEqual(nodes[saddle3['child_ids'][0]]['first_outlet_layer'], '3')

    def test_nodata_is_barrier_not_zero(self):
        # NoData 格不进入图：25 格中 1 个屏障，结果只有 24 格。
        v, _ = load('data/normal.json')
        r = analyze(v)
        self.assertEqual(len(r['cells']), 24)
        # 把屏障换成有效零高程后，原本被屏障分开的两块低洼地应在同一分量内。
        filled = json.loads(json.dumps(v))
        filled['elevation'][2][2] = '0'
        rf = analyze(filled)
        self.assertEqual(len(rf['cells']), 25)
        self.assertEqual(len(rf['components']), 1)
        # 面积必须为正：零面积被拒。
        bad = json.loads(json.dumps(v))
        bad['area'][0][0] = '0'
        with self.assertRaises(ValueError):
            analyze(bad)

    def test_invalid_inputs(self):
        v, _ = load('data/invalid-domain.json')
        with self.assertRaises(ValueError):
            analyze(v)
        # 出口落在 NoData、出口重复、max_cells 超限。
        base, _ = load('data/normal.json')
        bad = json.loads(json.dumps(base))
        bad['outlets'] = [[2, 2]]
        with self.assertRaises(ValueError):
            analyze(bad)
        bad = json.loads(json.dumps(base))
        bad['outlets'] = [[2, 4], [2, 4]]
        with self.assertRaises(ValueError):
            analyze(bad)
        bad = json.loads(json.dumps(base))
        bad['max_cells'] = 10
        with self.assertRaises(ValueError):
            analyze(bad)

    def test_json_output_is_finite_and_fraction_strings(self):
        v, _ = load('data/multi-outlet-plateau.json')
        text = json.dumps(to_jsonable(analyze(v), with_paths=True),
                          allow_nan=False)
        self.assertNotIn('Infinity', text)
        self.assertNotIn('NaN', text)
        j = json.loads(text)
        self.assertIsInstance(j['cells'][0]['elevation'], str)

    def test_heap_and_dsu_agree_everywhere(self):
        v, _ = load('data/normal.json')
        t = build_merge_tree(v)
        label, _ = minimax_heap(v)
        self.assertEqual(list(t['spill']), label)

    def test_tiny_and_boundary_grids(self):
        # 1x1 出口：H 即自身高程，零填洼，单叶节点。
        v = {'elevation': [['3/2']], 'area': [['7']], 'outlets': [[0, 0]]}
        j = to_jsonable(analyze(v))
        self.assertEqual(j['cells'][0]['spill_elevation'], '3/2')
        self.assertEqual(j['totals']['total_fill_volume'], '0')
        self.assertEqual(j['merge_tree']['node_count'], 1)
        # 边界格不自动成为出口：只有中间格是出口，两侧 H 被自身高程顶到 5。
        v = {'elevation': [['5', '1', '5']], 'area': [['1', '1', '1']],
             'outlets': [[0, 1]]}
        r = analyze(v)
        self.assertEqual([c['spill_elevation'] for c in r['cells']],
                         [Fraction(5), Fraction(1), Fraction(5)])
        # 全屏障与出口越界均拒绝。
        with self.assertRaises(ValueError):
            analyze({'elevation': [[None, None], [None, None]],
                     'area': [['1', '1'], ['1', '1']], 'outlets': []})
        with self.assertRaises(ValueError):
            analyze({'elevation': [['1', '2']], 'area': [['1', '1']],
                     'outlets': [[0, 9]]})


if __name__ == '__main__':
    unittest.main()
