# moon-causal

moon-causal 是面向 MoonBit 的实用型因果推断工具箱，覆盖数据契约、倾向得分、IPW/ATT/ATC/ATO、匹配、平衡诊断、异质性、敏感性、不确定性、DID、生存和确定性模拟 benchmark。

这是 2026 年 8 月 MoonBit 黑客松项目的结项版本，模块元数据使用 qjgqjg/moon-causal，许可证为 Apache-2.0。

## 本地验证

运行 moon fmt --check、moon check --deny-warn、moon test 和 moon run cmd/main。

当前本地结果为 8,073 行 MoonBit 源码、181 个测试全部通过。CLI 使用固定种子生成 200 条观测、4 个协变量、真实 ATE 1.75 的 benchmark，并输出估计值、绝对误差、标准误、有效样本量和 propensity AUC。

## 目录与工程化

- 根包：核心数据结构、数值线性代数、模型、估计器和诊断。
- 测试：边界、性质、回归、集成、压力、覆盖矩阵和 benchmark 测试。
- cmd/main：可直接运行的确定性 benchmark。
- .github/workflows/ci.yml：官方 stable 工具链、格式、类型、接口、wasm/native 构建与测试。

项目地址：<https://github.com/qjgqjg/moon-causal>

```mbt check
///|
test "README example" {
  inspect(mean([1.0, 3.0, 5.0]), content="3")
}
```
