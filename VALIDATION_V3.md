# v3 开发与验证记录

日期：2026-11-08。版本：0.3.0。所有新示例均为明确标记的合成参数，未进行实际 DFT/NEB 或实验标定。

## 1. 改造内容

新增 `dftkmc/`，保留 `atomistic/`、`vcm3d/` 及原 C 文件。研究范围为固定晶格、固定电荷分支、等温、均匀给定场下的空位守恒扩散与聚集。

| 文件 | 职责 |
|---|---|
| `dftkmc/geometry.py` | 单斜晶格、完整周期镜像、稳定通道编号、独立迁移图和相互作用图 |
| `dftkmc/catalog.py` | 结构指纹、参数与来源、完整性审计、显式通道排除和未知值拒绝机制 |
| `dftkmc/model.py` | 同一能量模型与共享鞍点、KMC 事件、标签和未折回位移守恒 |
| `dftkmc/analysis.py` | 整条独立轨迹的统计、有限时间示踪估计、截断处理与组间比较 |
| `dftkmc/workflow.py` | 可追溯输出、机制消融研究和 NEB 端点准备 |
| `dftkmc/cli.py` | `template / audit / run / study / export-neb` 命令 |

本版将构型相互作用和过渡态环境修正分开，支持三个默认消融组：`noninteracting`、`pair-only`、`pair-environment`。固定 q 之间单独比较动力学指标，不比较独立零点下的绝对能量。

## 2. 最终测试结果

完整命令 `python -m unittest discover -s tests -v`：**68 项全部通过**，测试框架报告 7.261 秒。

| 测试组 | 数量 | 主要证据 |
|---|---:|---|
| 原有 atomistic / vcm3d | 24 | 旧接口、守恒、场、热、电阻网络与等待时间基线 |
| 新目录审计 | 22 | 缺项、非有限数、单位与来源、结构指纹、共享鞍点、反向重复、显式排除、pending 拒绝 |
| 新独立物理验证 | 12 | 总能差、详细平衡、NEB 重复计能回归、独立邻接、非正交周期位移、守恒、场方向、两态解析过程 |
| 新工作流与统计 | 10 | 未完成轨迹处理、replica 分组、退化不确定性、记录重放、NEB 原子与镜像、覆盖范围传播 |

运行环境为 Python 3.12.14、NumPy 2.3.5、SciPy 1.17.0、Matplotlib 3.10.8。完整日志和环境分别位于 `examples/dftkmc/validation/unittest.txt` 与 `test_environment.json`。软件与数学测试通过不等于参数或材料预测完成验证。

## 3. v2 回归证据

源压缩包 `kmc_for_rram_atomistic_v2.zip` 的 SHA-256：

```
59705368eac7b9c70392ee9985ace7bd1bb6d7272d22b9db5154e726917b373a
```

`atomistic/`、`vcm3d/` 源码与恢复的 v2 相同。重复原默认算例：3×3×5、600 K、空位比例 0.03、seed=17、+1 V / −1 V 各 1 μs、最大 1000 事件。

- 完成整个 2 μs 波形，122 个事件：120 次 hop、1 次 oxygen_out、1 次 oxygen_in。
- 最终 11 个空位，电导约 `5.4867400010193706e-12 S`。
- `summary.json`、`trace.csv`、`events.csv`、`initial.extxyz`、`final.extxyz` 与原包逐字节一致。
- `states.npz` 的全部数组逐项相同；不要求 ZIP 容器时间戳等元数据相同。

机器可读比对结果见 `examples/dftkmc/validation/artifact_checks.json`。

## 4. v3 已实际运行的演示

### 单条轨迹

`examples/dftkmc/single_run/`：3×3×3、q=2、800 K、8 个空位、零场、seed=17、时长 1 μs、5000 事件上限。实际完成 **1060 次事件**并达到终止物理时间，未截断。

提供完整参数、逐事件来源和有效势垒、轨迹、标签、未折回位移及图。图上明确标记不是 DFT 标定参数。

### 机制对比

`examples/dftkmc/mechanism_study/`：q=0、q=2，各 3 个消融组，每组 8 条独立轨迹，共 **48 条**；每条都是 3×3×3、800 K、8 个空位、零场、0.2 μs，base seed=20261108。

48 条全部完成共同物理时长，无事件上限截断；总计 **14616 次事件**，每条 166–441 次。置信区间的重复单位是整条轨迹。D 的区间按 replica bootstrap 计算；团簇标量均值采用 Student t 区间。所有结果都是单一观察时间下的演示，不代表长时间扩散收敛。

| 固定 q | 环境组 | 平均事件数 | 最终平均空位对数 | 最终平均最大簇 |
|---:|---|---:|---:|---:|
| 0 | noninteracting | 410.750 | 1.000 | 1.875 |
| 0 | pair-only | 385.875 | 2.250 | 2.625 |
| 0 | pair-environment | 323.500 | 2.500 | 2.750 |
| 2 | noninteracting | 249.375 | 1.750 | 2.375 |
| 2 | pair-only | 249.875 | 2.875 | 3.125 |
| 2 | pair-environment | 207.625 | 1.500 | 2.125 |

表中是实际输出的合成示例统计，不是 HfO₂ 的材料结论。不同 q 的任意演示鞍点由明确公式生成，不能将这些数值趋势称为已发现的真实电荷效应。各组初态都按独立种子均匀随机放置空位；尚未证明达到热平衡，单个终态簇大小不等于团簇寿命或保持稳定性。

`study.json` 和 `groups.csv` 包含区间、扩散张量与全部对照；`replicas.csv` 保留原始重复结果；每条 replica 子目录保留事件和状态。`mechanism_comparison.png` 由保存的数据生成，重绘脚本为 `examples/dftkmc/plot_study.py`。

## 5. DFT 工作清单及导出核对

`catalog.illustrative.json` 审计为内部完整；必须使用 `--allow-illustrative` 才能由 CLI 运行。`catalog.pending.json` 有 **158 个未填数值**，审计为未就绪，不能运行。

当前几何目录有 **34 条候选通道**、q=0/2 两个固定分支。`neb_worklist/` 包含 **68 组待弛豫端点**，每端点 108 Hf + 215 O = 323 原子。对全部 68 组重新核对了元素数、原子顺序、仅指定氧原子移动、移动向量与实际周期镜像一致。

每组含 `POSCAR.initial`、`POSCAR.final` 和 `mapping.json`；总表为 `tasks.csv`。没有计算势垒、分发 POTCAR 或提交作业。导出器重新枚举几何候选；源目录的任何排除理由作为来源保留，而不自动作用于新 pending 工作清单。

本版未建立全空间群下的最小不等价路径集合，没有完成超胞尺寸、泛函或 NEB 图像数收敛。真实数据接入步骤见 `DFT_INPUT_PLAN.md`。

## 6. 复现命令

从解压后的项目根目录执行。已安装 NumPy、SciPy、Matplotlib 时可直接使用 `python -m`；也可 `python -m pip install -e .`。输出目录必须新建或为空，避免覆盖已有结果。

```bash
python -m unittest discover -s tests -v

python -m dftkmc audit examples/dftkmc/catalog.illustrative.json

python -m dftkmc run \
  --catalog examples/dftkmc/catalog.illustrative.json --allow-illustrative \
  --repeats 3 3 3 --charge 2 --temperature 800 --nvac 8 \
  --duration 1e-6 --max-events 5000 --sample-every 10 --seed 17 \
  --plot --output results/v3_single

python -m dftkmc study \
  --catalog examples/dftkmc/catalog.illustrative.json --allow-illustrative \
  --repeats 3 3 3 --charges 0 2 --temperature 800 --nvac 8 \
  --duration 2e-7 --replicas 8 --max-events 5000 --sample-every 20 \
  --bootstrap-samples 2000 --seed 20261108 --output results/v3_study

python examples/dftkmc/plot_study.py results/v3_study/study.json

python -m dftkmc export-neb \
  --catalog examples/dftkmc/catalog.pending.json --repeats 3 3 3 \
  --charges 0 2 --output results/v3_neb
```

生成新模板：

```bash
python -m dftkmc template --status pending --charges 0 2 \
  --output results/new_catalog.pending.json
```

`audit` 对未就绪目录返回 1；参数错误返回 2；轨迹/研究被事件上限截断返回 3。保留这些状态是结果解释的一部分。

## 7. 仍待完成的科学验证

真实位点、迁移及局域环境数据；长程或更细致的相互作用模型；时间、尺寸和样本数收敛；具体电荷态的稳定性与转换机制。当前 `symmetric_pair` 是可逆的声明近似，未实现逐环境 NEB 查表。几何簇不等于导电细丝，固定电荷分支不等于动态电荷转变，有限时间示踪估计不等于实际器件行为。
