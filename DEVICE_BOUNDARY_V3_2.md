# v3.2：电极—氧化物—电极边界使用指南

本版在 `dftkmc/` 中加入**两端给定电势、氧阻挡、有限厚度**的器件边界。氧化物仍是固定晶格、固定空位数、固定电荷分支、等温 KMC；电极用等势平面表示。它可以研究给定电压下的空位再分布、接触区富集及几何团簇变化。

这里的“器件边界”尚不包含电子输运、界面氧交换或自洽空间电荷。`atomistic/`、`vcm3d/` 和历史 C 程序继续保留原有行为。**本次已执行 126 项测试和 3 组虚拟参数器件轨迹；结果与适用范围见第 7 节。**

## 1. 几何、电势与电场

设超胞三个行向量为 \(\mathbf A,\mathbf B,\mathbf C\)。上下电极平行于 \(\mathbf A,\mathbf B\) 张成的平面，底电极经过原点。定义

\[
\mathbf n=\frac{\mathbf A\times\mathbf B}{|\mathbf A\times\mathbf B|},\qquad
L_{\mathrm{Å}}=\mathbf C\cdot\mathbf n>0,\qquad
s_i=\mathbf r_i\cdot\mathbf n.
\]

底、顶电极分别位于 \(s=0\)、\(s=L\)。**斜晶胞厚度取法向投影 \(\mathbf C\cdot\mathbf n\)，不能用 \(|\mathbf C|\)**；也不以最外层氧原子的间距代替电极间距。沿 \(\mathbf A,\mathbf B\) 方向周期，沿第三个晶胞方向有限；这对应代码中的 `boundary="slab"`。对于旋转或倾斜晶胞，法向不必等于全局笛卡尔 z 方向。

当前唯一电势模型为 `linear_dielectric`。在均匀介质、背景自由体电荷忽略、端点电势给定的条件下，采用 Laplace 方程的解析解：

\[
\phi(s)=V_b+(V_t-V_b)\frac{s}{L},\qquad
\mathbf E\,[\mathrm{V/nm}]=-10\frac{V_t-V_b}{L_{\mathrm{Å}}}\mathbf n.
\]

换算中的 10 来自 \(1\,\mathrm{nm}=10\,\mathrm{Å}\)。当顶电极电势高于底电极时，电场指向底电极；正电荷空位在统计上偏向低电势方向迁移，但单次随机跳迁不必沿电场方向。

这一近似的依据是均匀介质中的无自由体电荷边值问题；固定端点电势时，均匀介电常数不改变该电势解和电场形状。参见 [MIT 电磁学在线教材，第 6.6 节：Piece-Wise Uniform Electrically Linear Dielectrics](https://web.mit.edu/6.013_book/www/chapter6/6.6.html)。本版并未据此计算电容、电极表面电荷或电场储能。

**空位作为给定背景电场中的测试电荷处理。** 空位移动后不重新求解 Poisson 方程，其带电分布不反过来改变 \(\phi\)。因此“无体电荷的 Laplace 背景”是场近似，不是宣称含带电空位的真实器件没有空间电荷。

## 2. 配置与接触区的含义

`device.json` 的完整结构如下。`DeviceConfig` 和 `ElectrodeConfig` 是不可变配置对象；JSON 加载要求字段完整，拒绝缺失或未知字段。

```json
{
  "schema_version": 1,
  "potential_model": "linear_dielectric",
  "bottom": {
    "potential_V": 0.0,
    "contact_depth_A": 2.0,
    "oxygen_boundary": "blocking"
  },
  "top": {
    "potential_V": 1.0,
    "contact_depth_A": 2.0,
    "oxygen_boundary": "blocking"
  }
}
```

| 字段 | 含义 |
|---|---|
| `potential_V` | 电极平面的电势，单位 V；作用于运动的是两端电势差 |
| `contact_depth_A` | 从相应电极面向氧化物内部量取的标记深度，单位 Å |
| `oxygen_boundary` | 当前只支持 `blocking`：氧不进入或离开器件，空位数守恒 |
| `potential_model` | 当前只支持 `linear_dielectric`，即上述给定线性电势 |

底接触区满足 \(s_i\le d_b\)，顶接触区满足 \(L-s_i\le d_t\)，几何判定容差为 \(10^{-9}+10^{-12}L_{\mathrm{Å}}\) Å。两个深度须为正，\(d_b+d_t<L\)，每个接触区都须包含氧化物位点且彼此不重叠；不满足时应报错，不自动扩层。

输出中的三类标签为 `bottom_contact`、`top_contact` 和 `interior`。接触区只是**氧化物位点标签**。修改其深度只改变分区统计，不把这些位点转成金属、不把整层电势钳成端点值，也不缩短电压降落的厚度 \(L\)。接触区内仍使用线性 \(\phi(s_i)\) 和声明的氧化物跃迁参数；未增加接触区专属的形成能、迁移势垒或界面反应。

阻挡边界意味着无跨越上下电极面的跳迁、无氧库交换事件。标为接触区的空位仍可在合法氧化物位点之间运动，接触区并非吸收区或冻结区。

## 3. 电压如何进入 KMC

这里的 \(q\) 是以基本电荷为单位的整数固定电荷分支，例如 `charge=2`。对位移为 \(\mathbf d_{ij}\,[\mathrm{Å}]\) 的空位跳迁，静电势能变化为

\[
\Delta U_{\mathrm{el},ij}\,[\mathrm{eV}]
=q\,[\phi_j-\phi_i]
=-\frac{q}{10}\,\mathbf E\,[\mathrm{V/nm}]\cdot\mathbf d_{ij}\,[\mathrm{Å}].
\]

该项通过已有事件字段 `field_work_eV` **只计入一次**。沿用的字段名在这里表示静电势能变化；电场对粒子所做的机械功符号相反。不得再把 \(q\phi_i\) 加进 catalogue 的 site energy 后重复施加同一场修正。

若某条 canonical 路径的过渡态比例为 \(f=\texttt{ts_fraction}\)，它在正向势垒中增加 \(f\Delta U_{\mathrm{el}}\)，在反向势垒中增加 \(-(1-f)\Delta U_{\mathrm{el}}\)。因此正反速率仍使用同一个过渡态，满足

\[
\frac{k_{i\to j}}{k_{j\to i}}
=\exp\!\left[-\frac{\Delta E_{\mathrm{chem},ij}+q\Delta\phi_{ij}}{k_B T}\right],
\]

其中化学项沿用 catalogue 的 site/pair 模型；声明的对称环境修正仍须单独保留其近似性质。若电压使势垒进入当前激活事件模型不支持的范围，应报告模型域错误，不靠任意截断势垒制造可运行结果。

### 电势零点与能量输出

同时给两端电势增加常数 \(V_0\) 时，\(\Delta V\)、电场、每个跳迁的 \(q\Delta\phi\) 和速率均不变。固定空位数 \(N_v\)、固定 \(q\) 时，绝对 \(q\sum_i v_i\phi_i\) 只增加常数 \(N_v qV_0\)。为便于比较，本版输出以底电极为参考的静电势能：

\[
U_{\mathrm{el,rel}}=q\sum_i v_i(\phi_i-V_b),\qquad
U_{\mathrm{total,rel}}=E_{\mathrm{chem}}+U_{\mathrm{el,rel}}.
\]

因此 `electrostatic_energy_relative_bottom_eV` 与 `total_energy_relative_bottom_eV` 在共同平移电势零点后保持不变。它们是缺陷在外加势中的能量，不乘自能表达式中的 1/2 因子，也不是器件电场储能或热量；也不能用不同固定电荷分支的独立能量零点比较电荷态稳定性。

## 4. CLI：生成配置并运行

以下示例使用合成参数，需明确写 `--allow-illustrative`。真实 DFT 参数应来自已完成、可追溯的计算；pending catalogue 仍不能运行。

```bash
python -m dftkmc template \
  --status illustrative \
  --output illustrative_device_catalog.json

python -m dftkmc device-template \
  --output device.json \
  --top-voltage 1 \
  --bottom-voltage 0 \
  --bottom-contact-depth 2 \
  --top-contact-depth 2

python -m dftkmc run \
  --catalog illustrative_device_catalog.json \
  --allow-illustrative \
  --device device.json \
  --repeats 3 3 8 \
  --charge 2 \
  --temperature 800 \
  --nvac 8 \
  --duration 1e-6 \
  --max-events 20000 \
  --sample-every 10 \
  --plot \
  --output results/device_example
```

`--device` 自动采用 slab 边界。显式指定 `--boundary bulk` 或同时给出 `--field` 会被拒绝，避免两套边界或两次施场；即使显式 `--field 0 0 0` 也不应混用。输出目录需为新目录或空目录，配置文件使用未被占用的路径。

`duration` 是目标物理时间；`max-events` 是计算上限。若先碰到事件上限，结果属于截断轨迹，不能把其终点解释为指定时刻的响应或物理稳态。跨电压、厚度或 seed 的比较仍应使用相同观察时间并报告未完成 replica。

器件 slab 存在有限边界，通常还施加非零场，因此不能直接套用 bulk 零场的 tracer diffusion 或八态无限周期解析参考。设备运行中的位移、接触区占据和团簇统计不等同于电流或 I–V 曲线。

### Python 配置

```python
from dftkmc.device import DeviceConfig, ElectrodeConfig
from dftkmc.model import SimulationConfig

device = DeviceConfig(
    bottom=ElectrodeConfig(potential_V=0.0, contact_depth_A=2.0),
    top=ElectrodeConfig(potential_V=1.0, contact_depth_A=2.0),
)

config = SimulationConfig(
    repeats=(3, 3, 8),
    boundary="slab",
    device=device,
    charge=2,
    temperature_K=800.0,
    seed=17,
)
```

Python API 必须显式使用 `boundary="slab"`。配置还可通过 `DeviceConfig.load(path)` / `.write(path)` 或 `.from_dict(data)` / `.to_dict()` 读写；不要直接修改 frozen 配置对象。

## 5. 输出怎么读

| 输出 | 新增内容与用途 |
|---|---|
| `summary.json` | `device_boundary`、`effective_field_V_nm`，用于核对实际采用的设备配置与解析电场；同时保留参数校准状态、catalogue hash 和截断状态 |
| `device_boundary.json` | 与 `summary.device_boundary` 一致；包含 `normal`、`thickness_A`、`planes`、`contacts`、`applied_voltage_V`、`field_V_nm` 及模型限制 |
| `device_sites.csv` | 每个氧化物位点的 `site_id, kind, x_A, y_A, z_A, depth_A, potential_V, region`，另含 `potential_relative_bottom_V` 和两个 contact 布尔标签 |
| `states.npz` | 在原有时序空位状态之外保存静态 `potential_V`、`potential_relative_bottom_V`、`device_depth_A`、`bottom_contact`、`top_contact` 数组；位点轴与 lattice site ID 对齐 |
| `trace.csv` | 增加 `electrostatic_energy_relative_bottom_eV`、`total_energy_relative_bottom_eV`、`bottom_contact_vacancies`、`top_contact_vacancies`、`interior_vacancies` |
| `events.csv` | 增加 `source_potential_V`、`target_potential_V`、`source_region`、`target_region`、`applied_voltage_V`；`field_V_nm` 应为设备导出的实际场 |
| `device_boundary.png` | 使用 `--plot` 时生成，展示电势/接触区、空位分布和分区占据随时间的变化；合成参数明确标注 `ILLUSTRATIVE PARAMETERS / NOT DFT CALIBRATED`，同时保留原有 `overview.png` |

实际电场应以 `effective_field_V_nm` 或事件中的设备场为准，不把 `config.field_V_nm` 的默认零向量误读为零偏压。应满足三类空位数之和始终等于 \(N_v\)。某接触区空位增多表示邻电极氧化物区域富集；几何团簇增大只表示相应 interaction graph 上的关联，不能据此认定导电细丝或阻变开关已经形成。固定电压、闭合阻挡 slab 的离子势是保守的，接触区富集也不表示持续穿过上下电极的离子电流。

## 6. 当前边界的物理范围

| 已有表示 | 尚需独立建模或标定的内容 |
|---|---|
| 两个给定电势的平面 | 显式金属原子、真实金属/氧化物结构弛豫及界面电子结构 |
| 固定线性背景电势、带电空位测试粒子 | 随缺陷/载流子分布变化的自洽 Poisson 场、空间电荷和屏蔽 |
| 氧阻挡与空位数守恒 | 氧吸收/释放、氧化还原、氧库化学势及电极氧交换势垒 |
| 固定 charge 分支的离子跳迁 | 电子注入、捕获/发射、动态电荷态占据与费米能级求解 |
| 给定温度下的几何轨迹 | 电子电流、电导、I–V、Joule 热与热反馈 |

`bottom` / `top` 是几何与电势标签，没有指定 Pt、TiN 或其他真实电极材料。配置中的电势不是功函数或 Schottky 势垒；本版没有凭空赋予这些量。接触区仍使用所声明的氧化物参数，因此运行成功也不构成表面化学或界面 DFT 验证。

## 7. 本次验证与示例结果

2026-11-11，版本 **0.3.2**。完整工程测试已实际执行，**126 项全部通过**，测试框架报告 13.813 秒。其中保留 v3.1 的 102 项，新增几何与配置验收 14 项、器件物理与工作流集成验收 10 项。日志与环境在 `examples/device_boundary_v3_2/validation/`。

```bash
python -m unittest discover -s tests -v
```

| 验收项目 | 核对的性质 | 本次记录 |
|---|---|---|
| 斜晶胞几何 | 法向与厚度的独立检查、刚体旋转协变、横向周期像不变 | 通过 |
| 电势与符号 | 端点电势、线性分布、数值梯度、场方向及单位换算 | 通过 |
| 跳迁能量 | `field_work_eV = q * (target_potential_V - source_potential_V)`；总能差无重复计能 | 通过 |
| 详细平衡与零点 | 正反速率比、闭合路径零做功、共同电势平移后的种子轨迹及相对能量一致 | 通过 |
| 中性对照 | `q=0` 在 −1/0/+1 V 下保持相同种子轨迹 | 通过 |
| 阻挡与分区 | 空位数和标签守恒，无上下周期连接，三类分区数之和恒定 | 通过 |
| 输入冲突与几何域 | bulk/field 冲突、非有限值、空/重叠接触层、未知与重复 JSON 字段明确拒绝 | 通过 |
| 输出一致性 | 真实工作流的配置、有效电场、位点电势/分区、事件日志及 NPZ 对齐 | 通过 |
| 统计范围 | 有限 slab 即使零偏压也不输出 bulk diffusion 推断 | 通过 |
| 初始化前提 | 未初始化的器件状态查询明确报错，无空数组统计警告 | 通过 |

### 7.1 随包的三组实际轨迹

使用 v3.1 已导入的同一套合成数据，不修改势垒以追求特定轨迹。配置为 `3×3×8` 超胞、576 个氧晶格位点、8 个初始空位、固定 `q=2`、800 K、seed=17。法向厚度为 **4.1774306956 nm**，每侧接触层深度 2 Å，对应每侧 27 个氧晶格位点。三个条件的底电极均为 0 V，顶电极为 −1、0、+1 V。

三个条件使用相同初始化规则和 seed，初始空位位点均为 `[61, 92, 212, 262, 320, 421, 437, 481]`，平均初始深度 **1.8297146447 nm**。它们是配对的边界演示，不能当作独立随机重复或用于估计形成概率。

| 顶电极电势 | 法向电场（V/nm） | 实际事件数 | 末态平均深度（nm） | 末态底接触/内部/顶接触空位数 |
|---:|---:|---:|---:|---:|
| −1 V | +0.23938159 | 984 | 4.07743345 | 0 / 1 / 7 |
| 0 V | 0 | 2118 | 1.48834023 | 2 / 6 / 0 |
| +1 V | −0.23938159 | 856 | 0.14474797 | 6 / 2 / 0 |

三条轨迹全部达到预先设定的 **1 μs**，共 **3958 次事件**，均未触及 20000 事件上限。保存状态中的空位总数始终为 8。正偏压下空位向底电极富集，负偏压下向顶电极富集，符合当前模型的正测试电荷受力方向。零偏压的单条有限样本轨迹允许有随机再分布，不要求末态平均位置恰好等于初态。

原始结果在 `examples/device_boundary_v3_2/negative/`、`zero/`、`positive/`，汇总在 `demo_report.json`。每组保存全部事件、采样状态、配置、位点标签和图。`positive/device_boundary.png` 已进行图面检查，其图数据来自已保存的实际轨迹。

复现相同三组演示：

```bash
python examples/device_boundary_v3_2/run_demo.py \
  --allow-illustrative --output results/device_boundary
```

单独复现 +1 V 案例：

```bash
python -m dftkmc run \
  --catalog examples/device_boundary_v3_2/catalog.virtual.json \
  --allow-illustrative \
  --device examples/device_boundary_v3_2/device_positive.json \
  --repeats 3 3 8 --charge 2 --temperature 800 --seed 17 --nvac 8 \
  --duration 1e-6 --max-events 20000 --sample-every 10 --plot \
  --output results/device_positive
```

这些结果验证声明的理想电极边界与离子 KMC 工作流。参数仍为 `illustrative`，未标定特定金属/氧化物界面，未计算电流或证明导电细丝形成。此次同时纠正了旧 README/MODEL 对随机粗粒化示例的描述：旧轨迹曾有短暂几何贯通，最终断开；旧代码及其结果数据没有修改。
