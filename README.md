# kmc_for_rram — 第一性原理参数接口与氧缺陷 KMC

本分支发布 v4.0 源码、回归测试、可复现输入和精简结果。六组主算例及对照的记录见
[`protocol.json`](examples/coupled_device_v4/protocol.json) 和
[`demo_report.json`](examples/coupled_device_v4/demo_report.json)。
完整事件日志、状态快照、NEB 端点输出和交互图等生成物不随 Git 分支分发。
下文所述随原始项目包提供的生成物，可按相应命令复现。

## v4.0：有限氧库、界面交换与局域电势反馈

新增 `device-run` 器件引擎：界面氧移出/填回、逐原子编号与氧/反应电子账本、独立有效电阻图、Kirchhoff 电势、实际限流电压及逐事件反馈。每次事件后更新电子网络；下一次事件速率使用冻结的当前网络电势。原有 `run` / `study` 的固定空位数模型及全部旧模块保留。

```bash
python examples/coupled_device_v4/run_demo.py \
  --allow-illustrative --workers 3 --output results/coupled_device
```

这会执行预先声明的 3 个偏压种子和零偏压、关闭交换、恒定电导 3 个对照。所有算例从 0 空位开始，波形与判据在运行前写入 `protocol.json`，不筛选成功种子。查看各组的 `coupled_overview.png`，或在浏览器打开自包含的 `viewer.html` 播放已保存的实际状态。完整说明与验证见 [COUPLED_DEVICE_V4.md](COUPLED_DEVICE_V4.md)。

**当前新增界面化学与电子电导参数均为合成数据。** 这是可以运行的有效模型形成过程演示，尚不是已标定的真实电极堆栈预测。Kirchhoff 电压作为离子驱动电势是明确的近似；未求解空间电荷 Poisson 方程、动态电荷态或热反馈。有限偏压时，只验证同一冻结电势下的局部详细平衡，不宣称更新场后的全局平衡自由能。几何贯通、电导阈值与持续时间分别记录。

## v3.2：电极—氧化层—电极边界

新增 `dftkmc/device.py`：上下理想电极、可配置两端电势、氧化物接触区域、侧向周期和法向氧阻挡边界。薄膜厚度使用单斜超胞的法向投影；电场由两端电压和厚度决定，并沿用共享过渡态的带符号能量修正。接触区标签不把氧化物变成金属，也不产生电子电流。完整说明见 [DEVICE_BOUNDARY_V3_2.md](DEVICE_BOUNDARY_V3_2.md)。

```bash
python -m dftkmc device-template --top-voltage 1 --output results/device.json
python -m dftkmc run --catalog examples/virtual_data_v3_1/catalog.virtual.json \
  --allow-illustrative --device results/device.json --repeats 3 3 8 \
  --charge 2 --temperature 800 --nvac 8 --duration 1e-6 \
  --max-events 20000 --sample-every 10 --plot --output results/device_run
```

`--device` 自动选择 `slab`，与显式 `--field` 或 `--boundary bulk` 冲突。输出增加两端边界、逐位点电势与接触标签、接触区空位数、外场能和器件边界图。当前仍为固定电荷、等温、固定空位数和线性电势近似；界面氧交换、电子输运、空间电荷与热反馈尚未接入。随包的 `examples/device_boundary_v3_2/` 提供配置、实际轨迹和验证记录。

## v3.1：已接通虚拟计算数据与批量验证

2026-11-11：新增结构化的虚拟端点/正反势垒输入、严格导入、单空位周期网络的确定性参考，以及温度、观察时间、环境消融和外场响应测试。标准批次已执行 24 组条件、480 条轨迹；结果位于 `examples/virtual_data_v3_1/`，方法与解读见 [VIRTUAL_DATA_TESTS.md](VIRTUAL_DATA_TESTS.md)。所有新参数均明确标为合成数据。

```bash
python -m dftkmc.virtual_suite --allow-illustrative --workers 4 \
  --output results/virtual_data
```

`--profile smoke` 可快速检查导入与输出流程。正式流程从原始虚拟数据生成目录，再与相同物理时刻的解析矩参考比较；不将有限时间结果直接当成长期扩散极限，也不通过更换种子强制置信区间覆盖参考。已有 v3 示例、旧模型与其验证记录保留。

## 当前开发入口：dftkmc/

v3 新增固定电荷态、局域环境相关的氧空位扩散与聚集模型，v3.2 增加理想电极边界下的有限薄膜模式。通过可追溯的迁移通道目录、共享过渡态能量、独立相互作用图和多随机种子消融实验，研究电荷分支与局域环境对动力学的影响。

**本版已实现可运行的第一性原理参数接口；提供的演示参数是合成数据，尚未完成真实 DFT/NEB 标定。** 固定电荷分支不等于动态电子俘获。新模块的几何聚集指标不等于导电细丝或器件开关。

完整方程、数据规范、命令和局限见 [DFT_KMC_V3.md](DFT_KMC_V3.md)。运行 `python -m dftkmc --help` 查看 `template`、`device-template`、`coupled-template`、`audit`、`run`、`device-run`、`study`、`export-neb`。随包的 `examples/dftkmc/` 包含明确标记的演示数据、待填的 DFT 目录和实际运行结果；其复现与验证记录见 [VALIDATION_V3.md](VALIDATION_V3.md)。

第一批真实计算的结构、电荷、能量映射、三类试算路径及局域环境数据计划见 [DFT_INPUT_PLAN.md](DFT_INPUT_PLAN.md)。

原有 `atomistic/`、`vcm3d/` 和 C 示例保留原样。它们是不同层次的模型，使用各自的运行入口；下面保留旧版使用说明。`device-run` 是 v4 新增的可变空位数入口，不能用它的动态人口直接套用 v3 的固定 cohort 扩散统计。

---

## 原有三维氧缺陷 VCM 参考模型

本仓库保留原始 `kmc-ecm.c`（二维金属离子 ECM 示例），新增独立的 **三维氧空位型 VCM 动力学蒙特卡洛**实现 `vcm3d/`。这不是简单把二维数组扩展成三维，也不将金属沉积机制冒充氧空位机制。

**定位：可运行、可审查、可扩展的粗粒化研究原型；尚不是经过实验标定的 HfO₂ 器件预测工具。** 示例参数是明确标记的占位值，500 K 示例温度用于在短计算内观察演变，不能当成室温预测。详见 [MODEL.md](MODEL.md)。

## 主要改进

- 真正的 `(x,y,z)` 三维体网格，六邻居迁移，顶部/底部电极，侧面封闭。
- 独立的氧空位与氧间隙原子占据数组；成对生成/复合和两类迁移，严格保持封闭体系 `Nv = Ni`。
- 所有事件率组成统一目录，Gillespie 直接法选择事件并推进连续物理时间。
- 带符号的电势能差、非负迁移势垒、可配置的空位近邻结合能。
- 每次事件重新解三维电阻网络与稳态热传导，电流合规先改变器件电压，再计算发热和速率。
- 定电压脉冲边界处理、随机种子、事件日志、三维状态快照、数值守恒检查。
- 可重复单元测试。保留原始 C 文件不变，原始文件不属于本实现的运行依赖。

## 运行

需要 Python 3.10+、NumPy、SciPy、Matplotlib。在已经有这些库的环境中无需安装本项目：

```bash
python -m unittest discover -s tests -v
python -m vcm3d --shape 5 5 8 --initial random \
  --pulse 2 0.000001 --pulse -2 0.000001 \
  --max-events 2000 --sample-every 5 --output results/random_pulses
python -m vcm3d --shape 5 5 8 --initial filament \
  --pulse -2 0.000001 --max-events 2000 --sample-every 5 \
  --output results/seeded_relaxation
```

也可以安装：`python -m pip install -e .`，然后使用 `vcm3d ...` 命令。模拟代码不联网。

修改参数：`python -m vcm3d --config examples/illustrative.json ...`。JSON 只需提供覆盖字段，其余使用 `Config` 默认值。所有输入/输出单位见 `MODEL.md`。

### 输出

- `summary.json`：完整参数、初态、脉冲、种子、事件统计、实际结束时间和结束原因。
- `trace.csv`：实际器件电压、电流、冻结态电导/电阻、0.1 V 假想读电流、温度、缺陷数、连通状态、守恒残差。
- `events.csv`：每个事件的时间、类型、位置、能量差和当时速率。
- `states.npz`：真实三维空位/间隙原子、电势和温度序列。
- `overview.png`：初末态三维图及时间曲线。

若 `status` 是 `event_limit_before_waveform_end`，只完成了部分脉冲，不能宣称整段波形完成。提高 `--max-events` 或缩短脉冲。`--no-plot` 可跳过绘图。

## 已运行示例（随代码提供）

`examples/results/` 是 seed=7、5×5×8 网格的真实计算结果，可用上述命令重现：

- `random_pulses`：±2 V、各 1 μs，完成 799 个事件。最终 46 个空位/46 个间隙原子；保存轨迹中有 16 个采样点出现模型内几何贯通，随后断开，末态不贯通。这是简化通路演变演示，尚未验证稳定细丝或实际 SET。
- `seeded_relaxation`：人为初始化一条贯通空位链，−2 V、1 μs，完成 44 个事件，链断开，冻结态电导由约 6.25×10⁻⁸ S 降至 1.18×10⁻¹¹ S。初态间隙原子随机分布，这说明模型能破坏初始导电通道，**不能据此认定已经验证实际器件的 RESET 机制或双极性开关**。

这些例子验证程序流程和数值一致性，不是实验拟合。切勿通过强制连通或修改图像来伪造 forming/SET/RESET。

## 下一步：成为具体器件模型

1. 明确材料、电极、几何尺寸、是否存在氧库/吸氧层以及可测的 I–V/脉冲数据。
2. 用 DFT/实验约束缺陷生成、迁移、复合、结合和界面交换参数及不确定性。
3. 引入有守恒记账的电极氧交换、氧库化学势、不同界面势垒，才能研究可信的非对称双极开关。
4. 根据器件导电机制替换二值局域电导：陷阱辅助隧穿、热激活跃迁、金属/亚氧化物通道等。
5. 多随机种子统计、网格/尺寸收敛、极性测试、温度/电流合规扫描，分别验证 SET、RESET、保持和循环分布。

## 文件

`vcm3d/model.py` 为物理/数值核心，`vcm3d/cli.py` 为运行及输出接口，`tests/` 为无外部服务的测试，`MODEL.md` 为完整假设和方程。

## 新增：真实单晶结构上的氧位点 KMC

`atomistic/` 基于实验单斜 HfO₂ 的 P2₁/c 原子坐标，不再使用立方六邻居格点。包含三/四配位氧、真实距离候选路径、顶部有限氧库、严格氧守恒和可旋转的离线 HTML 查看器。

运行：`python3 -m atomistic --output results/atomistic_demo`，然后打开 `results/atomistic_demo/viewer.html`。完整来源、方法、运行命令与限制见 [ATOMISTIC.md](ATOMISTIC.md)。**真实的是晶体几何，势垒/电导仍未标定；它不是 MD 或经过 DFT 验证的器件预测。**
