///|
fn group_values(
  values : Array[Double],
  treatment : Array[Bool],
  treated : Bool,
) -> Array[Double] {
  let result = Array::new()
  let n = if values.length() < treatment.length() {
    values.length()
  } else {
    treatment.length()
  }
  for i in 0..<n {
    if treatment[i] == treated {
      result.push(values[i])
    }
  }
  result
}

///|
fn group_weights(
  weights : Array[Double],
  treatment : Array[Bool],
  treated : Bool,
) -> Array[Double] {
  let result = Array::new()
  let n = if weights.length() < treatment.length() {
    weights.length()
  } else {
    treatment.length()
  }
  for i in 0..<n {
    if treatment[i] == treated {
      result.push(weights[i])
    }
  }
  result
}

///|
/// Computes stabilized inverse-probability weights for a binary treatment.
pub fn stabilized_ipw_weights(
  treatment : Array[Bool],
  propensity_scores : Array[Double],
) -> Array[Double] {
  let n = if treatment.length() < propensity_scores.length() {
    treatment.length()
  } else {
    propensity_scores.length()
  }
  let mut treated_count = 0
  for value in treatment {
    if value {
      treated_count += 1
    }
  }
  let treated_rate = if n == 0 {
    0.5
  } else {
    treated_count.to_double() / n.to_double()
  }
  let result : Array[Double] = Array::new(capacity=n)
  for i in 0..<n {
    let p = safe_probability(propensity_scores[i])
    if treatment[i] {
      result.push(treated_rate / p)
    } else {
      result.push((1.0 - treated_rate) / (1.0 - p))
    }
  }
  result
}

///|
/// Clips propensity scores to a common support interval.
pub fn trim_propensity_scores(
  propensity_scores : Array[Double],
  lower? : Double = 0.01,
  upper? : Double = 0.99,
) -> Array[Double] {
  let result = Array::new(capacity=propensity_scores.length())
  for score in propensity_scores {
    result.push(clamp(score, lower, upper))
  }
  result
}

///|
fn weighted_group_mean(
  values : Array[Double],
  treatment : Array[Bool],
  weights : Array[Double],
  treated : Bool,
) -> Double {
  let selected_values = Array::new()
  let selected_weights = Array::new()
  let n = values.length()
  if treatment.length() < n {
    return 0.0
  }
  if weights.length() < n {
    return mean_or(group_values(values, treatment, treated), 0.0)
  }
  for i in 0..<n {
    if treatment[i] == treated {
      selected_values.push(values[i])
      selected_weights.push(weights[i])
    }
  }
  weighted_mean(selected_values, selected_weights)
}

///|
fn weighted_difference_se(
  outcomes : Array[Double],
  treatment : Array[Bool],
  weights : Array[Double],
) -> Double {
  let treated_values = group_values(outcomes, treatment, true)
  let control_values = group_values(outcomes, treatment, false)
  let treated_weights = group_weights(weights, treatment, true)
  let control_weights = group_weights(weights, treatment, false)
  let treated_variance = weighted_variance(treated_values, treated_weights)
  let control_variance = weighted_variance(control_values, control_weights)
  let treated_ess = effective_sample_size(treated_weights)
  let control_ess = effective_sample_size(control_weights)
  let treated_term = if treated_ess == 0.0 {
    0.0
  } else {
    treated_variance / treated_ess
  }
  let control_term = if control_ess == 0.0 {
    0.0
  } else {
    control_variance / control_ess
  }
  (treated_term + control_term).sqrt()
}

///|
/// Estimates the ATE with inverse-probability weighting and a normal interval.
pub fn estimate_ipw_ate(
  outcomes : Array[Double],
  treatment : Array[Bool],
  propensity_scores : Array[Double],
  stabilized? : Bool = true,
) -> Estimate {
  let weights = if stabilized {
    stabilized_ipw_weights(treatment, propensity_scores)
  } else {
    ipw_weights_ate(treatment, propensity_scores)
  }
  let treated_mean = weighted_group_mean(outcomes, treatment, weights, true)
  let control_mean = weighted_group_mean(outcomes, treatment, weights, false)
  let se = weighted_difference_se(outcomes, treatment, weights)
  Estimate::from_standard_error(
    treated_mean - control_mean,
    se,
    outcomes.length(),
    effective_sample_size(weights),
    "ATE (IPW)",
  )
}

///|
/// Estimates ATT using inverse-probability weighting.
pub fn estimate_ipw_att(
  outcomes : Array[Double],
  treatment : Array[Bool],
  propensity_scores : Array[Double],
) -> Estimate {
  let weights = ipw_weights_att(treatment, propensity_scores)
  let treated_mean = weighted_group_mean(outcomes, treatment, weights, true)
  let control_mean = weighted_group_mean(outcomes, treatment, weights, false)
  let se = weighted_difference_se(outcomes, treatment, weights)
  Estimate::from_standard_error(
    treated_mean - control_mean,
    se,
    outcomes.length(),
    effective_sample_size(weights),
    "ATT (IPW)",
  )
}

///|
/// Estimates the ATE from two potential-outcome predictions (g-computation).
pub fn g_computation_ate(
  predicted_treated : Array[Double],
  predicted_control : Array[Double],
) -> Estimate {
  let n = if predicted_treated.length() < predicted_control.length() {
    predicted_treated.length()
  } else {
    predicted_control.length()
  }
  let differences = Array::new(capacity=n)
  for i in 0..<n {
    differences.push(predicted_treated[i] - predicted_control[i])
  }
  let denominator = if n == 0 { 1.0 } else { n.to_double().sqrt() }
  let se = std_dev(differences, sample=true) / denominator
  Estimate::from_standard_error(
    mean(differences),
    se,
    n,
    n.to_double(),
    "ATE (g-computation)",
  )
}

///|
/// Computes a weighted mean potential outcome for a treatment arm.
pub fn weighted_potential_mean(
  outcomes : Array[Double],
  treatment : Array[Bool],
  weights : Array[Double],
  treated : Bool,
) -> Double {
  weighted_group_mean(outcomes, treatment, weights, treated)
}

///|
/// Returns a design-effect diagnostic for a set of weights.
pub fn weight_design_effect(weights : Array[Double]) -> Double {
  let n = weights.length().to_double()
  let ess = effective_sample_size(weights)
  if ess == 0.0 {
    0.0
  } else {
    n / ess
  }
}
