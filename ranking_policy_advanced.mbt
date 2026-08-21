///|
/// Policy ranking entry with estimated incremental value.
pub struct PolicyRankEntry {
  index : Int
  score : Double
  treatment_effect : Double
  selected : Bool
  cumulative_value : Double
}

///|
/// Policy ranking summary.
pub struct PolicyRanking {
  entries : Array[PolicyRankEntry]
  selected_count : Int
  selected_fraction : Double
  total_value : Double
  average_value : Double
  passes : Bool
}

///|
/// Ranks units by an estimated individualized effect.
pub fn rank_by_effect(
  effects : Array[Double],
  costs? : Array[Double] = [],
  budget? : Double = 1.0e300,
) -> PolicyRanking {
  let order : Array[Int] = Array::new(capacity=effects.length())
  for i in 0..<effects.length() {
    order.push(i)
  }
  for i in 1..<order.length() {
    let value = order[i]
    let mut cursor = i
    while cursor > 0 && effects[order[cursor - 1]] < effects[value] {
      order[cursor] = order[cursor - 1]
      cursor -= 1
    }
    order[cursor] = value
  }
  let entries : Array[PolicyRankEntry] = Array::new(capacity=effects.length())
  let mut cumulative = 0.0
  let mut selected = 0
  for position in 0..<order.length() {
    let index = order[position]
    let cost = if index < costs.length() { costs[index] } else { 1.0 }
    let selected_now = cumulative + cost <= budget
    if selected_now {
      cumulative += effects[index] - cost
      selected += 1
    }
    entries.push({
      index,
      score: if index < effects.length() {
        effects[index]
      } else {
        0.0
      },
      treatment_effect: effects[index],
      selected: selected_now,
      cumulative_value: cumulative,
    })
  }
  {
    entries,
    selected_count: selected,
    selected_fraction: if effects.length() == 0 {
      0.0
    } else {
      selected.to_double() / effects.length().to_double()
    },
    total_value: cumulative,
    average_value: if selected == 0 {
      0.0
    } else {
      cumulative / selected.to_double()
    },
    passes: selected == 0 || is_finite(cumulative),
  }
}

///|
/// Computes policy value at each top-k prefix.
pub fn policy_prefix_curve(ranking : PolicyRanking) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(
    capacity=ranking.entries.length(),
  )
  for i in 0..<ranking.entries.length() {
    let entry = ranking.entries[i]
    result.push([
      (i + 1).to_double(),
      entry.cumulative_value,
      entry.treatment_effect,
    ])
  }
  result
}

///|
/// Computes regret against an oracle effect ranking.
pub fn policy_regret(effects : Array[Double], selected : Array[Bool]) -> Double {
  let n = effects.length().min(selected.length())
  let selected_values = Array::new()
  let oracle = effects[:n].to_owned()
  for i in 0..<n {
    if selected[i] {
      selected_values.push(effects[i])
    }
  }
  for i in 1..<oracle.length() {
    let value = oracle[i]
    let mut cursor = i
    while cursor > 0 && oracle[cursor - 1] < value {
      oracle[cursor] = oracle[cursor - 1]
      cursor -= 1
    }
    oracle[cursor] = value
  }
  let oracle_total = sum(oracle[:selected_values.length()].to_owned())
  oracle_total - sum(selected_values)
}

///|
/// Applies a minimum-effect policy constraint.
pub fn select_effect_threshold(
  effects : Array[Double],
  threshold : Double,
) -> Array[Bool] {
  let result = Array::new(capacity=effects.length())
  for effect in effects {
    result.push(effect >= threshold)
  }
  result
}

///|
/// Computes treatment-rate disparity across integer groups.
pub fn policy_group_rates(
  groups : Array[Int],
  selected : Array[Bool],
) -> Array[Array[Double]] {
  let n = groups.length().min(selected.length())
  let unique : Array[Int] = Array::new()
  for group in groups[:n] {
    if !unique.contains(group) {
      unique.push(group)
    }
  }
  let result : Array[Array[Double]] = Array::new(capacity=unique.length())
  for group in unique {
    let mut count = 0
    let mut treated = 0
    for i in 0..<n {
      if groups[i] == group {
        count += 1
        if selected[i] {
          treated += 1
        }
      }
    }
    result.push([
      group.to_double(),
      treated.to_double(),
      count.to_double(),
      if count == 0 {
        0.0
      } else {
        treated.to_double() / count.to_double()
      },
    ])
  }
  result
}

///|
/// Returns a policy that equalizes group selection rates up to a tolerance.
pub fn constrained_policy(
  effects : Array[Double],
  groups : Array[Int],
  target_rate : Double,
) -> Array[Bool] {
  let n = effects.length().min(groups.length())
  let result = Array::make(n, false)
  let rates = policy_group_rates(groups[:n].to_owned(), Array::make(n, true))
  for rate in rates {
    let group = rate[0].to_int()
    let candidates : Array[Int] = Array::new()
    for i in 0..<n {
      if groups[i] == group {
        candidates.push(i)
      }
    }
    for i in 1..<candidates.length() {
      let value = candidates[i]
      let mut cursor = i
      while cursor > 0 && effects[candidates[cursor - 1]] < effects[value] {
        candidates[cursor] = candidates[cursor - 1]
        cursor -= 1
      }
      candidates[cursor] = value
    }
    let take = (candidates.length().to_double() * clamp(target_rate, 0.0, 1.0))
      .round()
      .to_int()
    for i in 0..<take.min(candidates.length()) {
      result[candidates[i]] = true
    }
  }
  result
}

///|
/// Returns a compact ranking summary vector.
pub fn policy_ranking_summary(ranking : PolicyRanking) -> Array[Double] {
  [
    ranking.entries.length().to_double(),
    ranking.selected_count.to_double(),
    ranking.selected_fraction,
    ranking.total_value,
    ranking.average_value,
    if ranking.passes {
      1.0
    } else {
      0.0
    },
  ]
}
