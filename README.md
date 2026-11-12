# 地形洼地的最小溢流高程与层级合并

这是独立的防汛演练资料Python项目。`data/`保存虚构人数、观测、时序或方案，`response_data.io`提供严格数值读取与来源字节SHA-256，`response_data.domain`给出已有原始摘要，`response_data.cli`输出JSON。字段、量纲与分析口径见[资料约定](docs/data-contract.md)。这些数据用于离线复盘，不连接真实人员、传感器或指挥系统。

Python 3.11及以上即可直接使用，无第三方Python依赖。在项目根目录执行：

```sh
python -m unittest discover -s tests -v
python -m compileall -q response_data tests
python -m response_data.cli inspect data/normal.json
python -m response_data.cli depressions data/normal.json
python -m response_data.cli depressions data/normal.json --with-paths
```

`inspect`显示原始资料概况，不覆盖源文件。`depressions`执行离散栅格洼地分析：按四邻接求每格的最小溢流高程H（到任一配置出口的路径中“路径最大原高程”的最小值）、填洼深度与按格面积加权的填洼体积，并给出按原高程分层批处理的规范多叉合并树（同高格同批激活，不产生逐格等高二叉链）。NoData是不可穿越屏障而非零高程；无出口可达的连通分量标为不可排，H留null而不写无穷数。`--with-paths`为每个可排格附无环见证路径，其最大原高程恰为H。该分析只做离散地形连通核算，不模拟降雨，填洼体积不等同洪水蓄量预测。字段、节点ID与层级口径见[资料约定](docs/data-contract.md)。

`invalid-`前缀资料为不合法输入，其余变体使用相同字段保留不同业务边界（`enclosed-undrainable`含屏障封死的不可排分量，`multi-outlet-plateau`含多出口与分数面积，`single-plateau`验证同高平台只生成一个叶节点）。数值null须按字段语义区别于0，不在读取时自动填补。

Dockerfile固定基础镜像，可用`docker build -t response-data .`及`docker run --rm response-data`读取资料；`docker run --rm response-data python -m unittest discover -s tests -v`执行已有检查。运行所需资料全部在本地，无需登录或实时来源。
