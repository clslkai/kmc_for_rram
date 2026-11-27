# 虚拟数据测试实测报告

**所有输入能量与势垒均为合成数据；本次没有进行真实 DFT/NEB 计算。**

- 条件数：24；轨迹执行数：480。
- 不重复的初始化/动力学子流对：448；中性场控制明确使用配对流。
- 实际事件总数：63831；全部完成指定物理时长：True。
- 中性空位配对外场控制保持相同路径：True。
- 同种子完整事件轨迹重放核对：True。

## 单空位扩散：与相同时间的确定性参考比较

物理时长在抽样前按平衡平均跳率确定，所有重复运行到固定终止时间。初始化随机数与动力学随机数使用独立子流。
矩方程参考给出同一时刻的 D(t)；长时间 D∞ 只用于考察有限时偏差。区间是否覆盖参考值仅作诊断。

| 条件 | KMC / 精确 D(t) | 95% 区间覆盖参考 | 精确 D(t)/D∞ |
|---|---:|---|---:|
| diffusion_q0_T700_h50 | 1.0370 | True | 1.00057 |
| diffusion_q0_T700_h200 | 1.0212 | True | 1.00014 |
| diffusion_q0_T800_h50 | 1.0596 | True | 1.00056 |
| diffusion_q0_T800_h200 | 1.2129 | True | 1.00014 |
| diffusion_q0_T900_h50 | 1.2073 | True | 1.00056 |
| diffusion_q0_T900_h200 | 1.2084 | True | 1.00014 |
| diffusion_q2_T700_h50 | 0.7754 | False | 1.00053 |
| diffusion_q2_T700_h200 | 0.9406 | True | 1.00013 |
| diffusion_q2_T800_h50 | 0.8546 | True | 1.00054 |
| diffusion_q2_T800_h200 | 1.0267 | True | 1.00014 |
| diffusion_q2_T900_h50 | 0.8655 | True | 1.00055 |
| diffusion_q2_T900_h200 | 0.9084 | True | 1.00014 |

## 环境消融和外场响应

| 条件 | 平均事件数 | 平均空位对数 | 平均最大簇 | 平均 z 位移 / Å |
|---|---:|---:|---:|---:|
| mechanism_q0_noninteracting | 58.81 | 0.8125 | 1.688 | -0.07642 |
| mechanism_q0_pair-only | 63.38 | 1.875 | 2.25 | 0.4723 |
| mechanism_q0_pair-environment | 54.94 | 1.75 | 2.25 | 0.1961 |
| mechanism_q2_noninteracting | 312.4 | 1.375 | 2.125 | -0.7759 |
| mechanism_q2_pair-only | 292.5 | 2.562 | 3.125 | 0.7352 |
| mechanism_q2_pair-environment | 246.5 | 2.562 | 2.75 | 0.4378 |
| field_q0_F-0.15 | 42 | 0 | 1 | -0.5962 |
| field_q0_F+0.00 | 42 | 0 | 1 | -0.5962 |
| field_q0_F+0.15 | 42 | 0 | 1 | -0.5962 |
| field_q2_F-0.15 | 201.4 | 0 | 1 | -120.9 |
| field_q2_F+0.00 | 194.9 | 0 | 1 | -4.189 |
| field_q2_F+0.15 | 206.1 | 0 | 1 | 124.4 |

环境鞍点项在两个方向对称加入，能改变有限时动力学，但不改变同一构型能量模型的平衡分布。不同模型的瞬时簇大小排序不应预设。

## 温度扫描的表观斜率

下面是三个温度点的描述性斜率，未作为独立 NEB 势垒或实际材料激活能：

- q=0，exact_long_time_reference：0.81262 eV。
- q=0，finite_time_KMC_protocol：0.86004 eV。
- q=2，exact_long_time_reference：0.70670 eV。
- q=2，finite_time_KMC_protocol：0.69967 eV。

## 文件与复现

suite_plan.json 固化参数、时长和所有种子；replicas.jsonl 保存每条轨迹初态、终态、累计位移及通道计数；groups.csv 与 suite_summary.json 保存统计；sample_trajectory/ 保存一个实际条件的完整事件重放。
原始虚拟计算输入为 mock_dft_input.json，导入后目录为 catalog.virtual.json。导入将正反势垒转换到同一鞍点，并检查闭合关系。

![虚拟数据测试总览](overview.png)

## 解释范围

- Synthetic input only; no DFT/NEB calculations or material claims.
- The reference is exact for the declared single-vacancy periodic rate network, not for real hafnia.
- The mechanism study starts from uniform random multi-vacancy configurations; final cluster metrics are not equilibrium lifetimes.
- Neutral field groups use paired streams; every group has independent replicas internally, but those three groups are not independent of one another.
- Fields +/-0.15 V/nm test driven response, not a claimed weak-field Einstein relation.
