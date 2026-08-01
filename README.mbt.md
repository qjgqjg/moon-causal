# moon-causal

`moon-causal` 是一款专为 MoonBit 生态打造的因果推断基础数学/数值计算工具箱，作为底层支撑构建复杂因果图和数值模型。

## 核心特性
- **倾向得分估计 (Propensity Score)**
- **逆概率加权 (IPW)**
- **分层估计 (Stratification)**
- **标准化差异 (SMD)**
- **双重稳健估计 (Doubly Robust Estimation)**
- **敏感性分析 (Sensitivity Analysis)**

## 快速开始

```mbt check
///|
test "example" {
  inspect(1 + 1, content="2")
}
```

## 黑客松比赛声明

- **来源说明**：本项目为全新原创 MoonBit 模块，专为本次 8 月份的黑客松大赛构建。
- **技术定位**：致力于填补 MoonBit 生态在统计与因果推断领域的空白，避免与已有组件重复，提供高内聚的因果推断算子。
- **许可协议**：Apache-2.0