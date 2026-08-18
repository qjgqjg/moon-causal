///|
pub struct SubgroupEffect {
  name : String
  sample_size : Int
  treated_count : Int
  control_count : Int
  estimate : Double
  standard_error : Double
}

///|
fn subgroup_effect(
  outcomes : Array[Double],
  treatment : Array[Bool],
  indices : Array[Int],
  name : String,
) -> SubgroupEffect {
  let subgroup_outcomes = Array::new(capacity=indices.length())
  let subgroup_treatment = Array::new(capacity=indices.length())
  for index in indices {
    subgroup_outcomes.push(outcomes[index])
    subgroup_treatment.push(treatment[index])
  }
  let treated = group_values(subgroup_outcomes, subgroup_treatment, true)
  let control = group_values(subgroup_outcomes, subgroup_treatment, false)
  let estimate = mean(treated) - mean(control)
  let standard_error = if treated.length() == 0 || control.length() == 0 {
    0.0
  } else {
    (variance(treated) / treated.length().to_double() +
    variance(control) / control.length().to_double()).sqrt()
  }
  {
    name,
    sample_size: indices.length(),
    treated_count: treated.length(),
    control_count: control.length(),
    estimate,
    standard_error,
  }
}

///|
pub fn threshold_subgroup_effects(
  covariate : Array[Double],
  outcomes : Array[Double],
  treatment : Array[Bool],
  threshold : Double,
) -> Array[SubgroupEffect] {
  let below = Array::new()
  let above = Array::new()
  for i in 0..<covariate.length() {
    if covariate[i] < threshold {
      below.push(i)
    } else {
      above.push(i)
    }
  }
  [
    subgroup_effect(outcomes, treatment, below, "below-threshold"),
    subgroup_effect(outcomes, treatment, above, "at-or-above-threshold"),
  ]
}

///|
pub fn quantile_subgroup_effects(
  covariate : Array[Double],
  outcomes : Array[Double],
  treatment : Array[Bool],
  groups : Int,
) -> Array[SubgroupEffect] {
  let actual_groups = if groups < 1 { 1 } else { groups }
  let result : Array[SubgroupEffect] = Array::new(capacity=actual_groups)
  for group in 0..<actual_groups {
    let lower = quantile(
      covariate,
      group.to_double() / actual_groups.to_double(),
    )
    let upper = quantile(
      covariate,
      (group + 1).to_double() / actual_groups.to_double(),
    )
    let indices = Array::new()
    for i in 0..<covariate.length() {
      let is_last = group == actual_groups - 1
      if covariate[i] >= lower && (covariate[i] < upper || is_last) {
        indices.push(i)
      }
    }
    result.push(
      subgroup_effect(outcomes, treatment, indices, "quantile-{group + 1}"),
    )
  }
  result
}

///|
pub fn conditional_average_treatment_effect(
  predicted_treated : Array[Double],
  predicted_control : Array[Double],
) -> Array[Double] {
  let n = if predicted_treated.length() < predicted_control.length() {
    predicted_treated.length()
  } else {
    predicted_control.length()
  }
  let result = Array::new(capacity=n)
  for i in 0..<n {
    result.push(predicted_treated[i] - predicted_control[i])
  }
  result
}

///|
pub fn policy_value(
  outcomes : Array[Double],
  treatment : Array[Bool],
  individual_effects : Array[Double],
) -> Double {
  let n = if outcomes.length() < individual_effects.length() {
    outcomes.length()
  } else {
    individual_effects.length()
  }
  if n == 0 {
    return 0.0
  }
  let mut total = 0.0
  for i in 0..<n {
    if individual_effects[i] >= 0.0 && treatment[i] {
      total += outcomes[i]
    }
    if individual_effects[i] < 0.0 && !treatment[i] {
      total += outcomes[i]
    }
  }
  total / n.to_double()
}

///|
pub fn uplift_curve(
  outcomes : Array[Double],
  treatment : Array[Bool],
  scores : Array[Double],
  bins : Int,
) -> Array[Double] {
  let order = Array::new(capacity=scores.length())
  for i in 0..<scores.length() {
    order.push(i)
  }
  for i in 1..<order.length() {
    let key = order[i]
    let mut j = i
    while j > 0 && scores[order[j - 1]] < scores[key] {
      order[j] = order[j - 1]
      j -= 1
    }
    order[j] = key
  }
  let result = Array::new()
  let actual_bins = if bins < 1 { 1 } else { bins }
  for bin in 1..<=actual_bins {
    let count = scores.length() * bin / actual_bins
    let selected_y = Array::new(capacity=count)
    let selected_t = Array::new(capacity=count)
    for position in 0..<count {
      selected_y.push(outcomes[order[position]])
      selected_t.push(treatment[order[position]])
    }
    result.push(difference_in_means(selected_y, selected_t))
  }
  result
}
