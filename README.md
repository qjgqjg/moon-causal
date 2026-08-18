# moon-causal

moon-causal 是面向 MoonBit 的实用型因果推断工具箱，服务于教育评估、产品实验、公共政策分析和观测数据研究。它把数据校验、倾向得分、加权、匹配、平衡诊断、异质性、敏感性、不确定性和模拟基准组织成可组合的纯 MoonBit API。

本仓库是 2026 年 8 月 MoonBit 黑客松项目的结项版本，项目方向与开源大赛申报不同。原始申报书保留在仓库外层工作目录中，仓库本身包含可运行源码、测试、CLI 示例和 CI。

## 能力范围

- 数据集契约：行列校验、缺失值处理、筛选、汇总、标准化和特征变换。
- 建模：带标准化和 L2 正则的逻辑倾向模型、线性结果模型、设计矩阵与交互项。
- 因果估计：IPW/稳定化 IPW、ATT/ATC/ATO、g-computation、双重稳健、分层和最近邻匹配。
- 诊断与稳健性：SMD/方差比、共同支持、positivity、校准、有效样本量、Rosenbaum 敏感性分析。
- 拓展分析：多臂处理、连续剂量反应、DID、事件研究、合成控制、生存曲线、中介和策略价值。
- 不确定性：bootstrap、置换检验、jackknife、区间估计、误差指标和交叉验证。
- 工程化：确定性观测数据模拟、可重复 benchmark CLI、边界测试和 wasm/native CI。

## 快速开始

运行 moon check --deny-warn、moon test 和 moon run cmd/main。

CLI 会运行固定种子、200 条观测、4 个协变量、真实 ATE 为 1.75 的观测数据 benchmark。当前本地工具链输出为：

    sample_size=200
    true_ate=1.75
    estimated_ate=1.695501947510638
    absolute_error=0.0544980524893619
    standard_error=0.15562533888456542
    effective_sample_size=177.61530018366125
    propensity_auc=0.6822682268226823

该数据由 SyntheticConfig 和固定 UInt64 seed 生成；benchmark 会同时报告估计误差、标准误、有效样本量和 propensity AUC。测试会检查相同配置下数据和报告的关键字段保持确定性。

## 规模与质量证据

结项本地检查结果：

- MoonBit 源码 8,073 行，其中生产代码 5,350 行、测试代码 2,723 行，73 个 .mbt 文件。
- moon check --deny-warn 通过。
- 181 个测试全部通过，覆盖正常路径、空输入、单例、维度不一致、缺失/非有限值、极端概率、重复生存时间和确定性模拟。
- CLI benchmark 可重复运行；默认种子为 20260818UL。
- LICENSE 为 Apache-2.0；模块元数据、仓库地址和 README 已与 GitHub 项目一致。

## CI

.github/workflows/ci.yml 使用 MoonBit 官方 stable 安装脚本，并执行 moon fmt --check、moon check --deny-warn、moon info、生成接口差异检查、wasm-gc 测试、覆盖率摘要以及 native 检查/测试。这样可避免依赖过时的固定工具链版本，并覆盖组委会重点关注的 CI 缺失问题。

## 项目地址

- GitHub：<https://github.com/qjgqjg/moon-causal>
- Gitlink：<https://gitlink.org.cn/qjgqjg68/moon-causal>
- MoonBit：<https://www.moonbitlang.com/>

## 许可与原创说明

本项目使用 Apache License 2.0。源码为本项目围绕 MoonBit 因果推断场景编写的实现；未将第三方仓库源码复制进本仓库。若未来接入外部数据集或算法实现，应在对应文件和文档中单独标注来源、许可证及可复现实验方式。
