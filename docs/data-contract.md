# 栅格连通与层级

行列坐标为零基(r,c)，只用上下左右四邻，NoData=null是不进入图的屏障。elevation和area为精确分数字符串，area对所有格为正，体积只累计可排有效格的(area*(H-z))。无法接出口的格与屏障分开，不赋数值H。

合并树按原z的不同层高批处理。每一层同时激活该高程的所有格，将新格的四邻连接及其相接的全部较低连通分量统一处理。每个含新格的批后分量构成一个节点，孩子是它吸收的批前分量的当前树根；没有孩子则叶，单孩子可保留为带平台的延伸节点，但规范输出须一致保留，不能随输入顺序压缩。节点成员只列该层新激活格及子节点引用，避免复制全体成员。node_id由层高规范分数与排序的新格坐标构造；多出口的首次含出口层高仍等于H。

## 离散洼地分析（response_data.depression）

### 定义

- 出口仅取`outlets`配置的有效格；边界格不自动成为出口，出口重复或落在NoData上为非法输入。
- 可排格c的溢流高程 `H(c) = min over 路径p(c→任一出口) of max over 格x∈p of z(x)`，即越过地形接到真实出口所需的最低连续水位；不可排连通分量的H为null，不用无穷数伪装。
- 填洼深度 `H(c)-z(c)` 恒非负；单格填洼体积 `area(c)*(H(c)-z(c))`，总体积只对可排格求和。面积按每格独立给出，支持分数权重。
- 该结果是离散地形连通量，不模拟降雨过程，填洼体积不等于洪水蓄量预测。

### 两套独立实现互校

1. `build_merge_tree`：并查集按层高批处理，一次扫描同时得到合并树与每格H（格出生节点沿父链上溯首个含出口祖先的层高）。
2. `minimax_heap`：有界优先队列上的最劣路径（max半群Dijkstra），标签为路径峰值；严格更优才松弛，最终化集合保证父森林无环且确定。

`analyze`逐格断言两者的H一致，并核对每格出生节点的首次含出口层高。小图另由测试中的简单路径DFS穷举与阈值BFS两个独立基准复核；大图复杂度近线性（DSU的O(nα(n))，堆式O(n log n)），不穷举所有路径。

### 批处理与规范ID

- 同高全部格在同一批次激活，批内所有边统一合并；同高激活的新平台在没有更低孩子时形成单个新叶，逐格顺序不会凭空产生等高二叉洼地链。
- 同层多个节点按node_id规范排序；node_id = `z=<规范分数>|<r.c>,...`，只含本层新激活格坐标，配合child_ids引用即可还原全体成员而无需逐层复制。
- 节点输出：layer、new_cells（紧凑行主索引）、child_ids、parent_id、has_outlet、first_outlet_layer、subtree_cell_count、subtree_min_elevation。
- 同高格的输入顺序、边的处理顺序均不改变规范树（测试以随机置换验证）。

### 见证与父结构

- 每格输出`parent_index`：父森林中通向出口的下一邻格，出口以自身为父；不可排格为null。
- `--with-paths`时附`witness_path`，无环、相邻格四邻接、终点为真实出口，且路径最大原高程恰好等于H。
- 平台与多出口场景下满足H的父结构不唯一（等高程邻居任选其一都保持同一峰值），实现用堆的最终化次序给出确定的一个；H本身唯一。

### 输出字段

每格：index/row/col、elevation（原高程）、is_outlet、drainable、spill_elevation（H，不可排为null）、fill_depth、fill_volume、node_id、parent_index，可选witness_path。
分量：root_node_id、cell_count、min_elevation（该分量最低原高程）、drainable、contains_outlet。
合计：total_fill_volume、drainable_cell_count、undrainable_cell_count、drainable_area，并保留原始资料SHA-256与`inspect`摘要。
