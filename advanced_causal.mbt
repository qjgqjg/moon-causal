///|
pub fn estimate_atc(
  outcomes : Array[Double],
  treatment : Array[Bool],
  propensity_scores : Array[Double],
) -> Estimate {
  let n = if outcomes.length() < treatment.length() {
    outcomes.length()
  } else {
    treatment.length()
  }
  let weights = Array::new(capacity=n)
  for i in 0..<n {
    let p = safe_probability(propensity_scores[i])
    weights.push(if treatment[i] { (1.0 - p) / p } else { 1.0 })
  }
  let effect = weighted_treatment_effect(outcomes, treatment, weights)
  let standard_error = weighted_difference_se(outcomes, treatment, weights)
  Estimate::from_standard_error(
    effect,
    standard_error,
    n,
    effective_sample_size(weights),
    "ATC (IPW)",
  )
}

///|
pub fn estimate_overlap_ate(
  outcomes : Array[Double],
  treatment : Array[Bool],
  propensity_scores : Array[Double],
) -> Estimate {
  let weights = overlap_weights(treatment, propensity_scores)
  let effect = weighted_treatment_effect(outcomes, treatment, weights)
  let standard_error = weighted_difference_se(outcomes, treatment, weights)
  Estimate::from_standard_error(
    effect,
    standard_error,
    outcomes.length(),
    effective_sample_size(weights),
    "ATO (overlap weighting)",
  )
}

///|
pub fn estimate_from_strata(
  outcomes : Array[Double],
  treatment : Array[Bool],
  propensity_scores : Array[Double],
  strata_count : Int,
) -> Estimate {
  let strata = stratify_by_score(propensity_scores, strata_count)
  aggregate_stratified_effects(stratified_effects(outcomes, treatment, strata))
}

///|
pub fn effect_curve_by_trim(
  outcomes : Array[Double],
  treatment : Array[Bool],
  propensity_scores : Array[Double],
  cuts : Array[Double],
) -> Array[Estimate] {
  let result : Array[Estimate] = Array::new(capacity=cuts.length())
  for cut in cuts {
    let lower = clamp(cut, 0.0, 0.5)
    let indices = Array::new()
    for i in 0..<propensity_scores.length() {
      if propensity_scores[i] >= lower && propensity_scores[i] <= 1.0 - lower {
        indices.push(i)
      }
    }
    let filtered_y = Array::new()
    let filtered_t = Array::new()
    let filtered_p = Array::new()
    for index in indices {
      filtered_y.push(outcomes[index])
      filtered_t.push(treatment[index])
      filtered_p.push(propensity_scores[index])
    }
    result.push(estimate_ipw_ate(filtered_y, filtered_t, filtered_p))
  }
  result
}

///|
pub struct SmdProfile {
  before : Array[BalanceMetric]
  maximum_before : Double
  balanced_before : Int
}

///|
pub fn smd_profile(
  covariates : Array[Array[Double]],
  treatment : Array[Bool],
  names : Array[String],
  threshold : Double,
) -> SmdProfile {
  let before = balance_table(covariates, treatment, names, threshold~)
  {
    before,
    maximum_before: maximum_absolute_smd(before),
    balanced_before: balanced_covariate_count(before),
  }
}

///|
pub fn treatment_effect_by_threshold(
  covariate : Array[Double],
  outcomes : Array[Double],
  treatment : Array[Bool],
  thresholds : Array[Double],
) -> Array[SubgroupEffect] {
  let result : Array[SubgroupEffect] = Array::new(
    capacity=thresholds.length() * 2,
  )
  for threshold in thresholds {
    for
      effect in threshold_subgroup_effects(
        covariate, outcomes, treatment, threshold,
      ) {
      result.push(effect)
    }
  }
  result
}

///|
pub fn policy_treatment_rate(
  individual_effects : Array[Double],
  threshold : Double,
) -> Double {
  if individual_effects.length() == 0 {
    return 0.0
  }
  let mut selected = 0
  for effect in individual_effects {
    if effect >= threshold {
      selected += 1
    }
  }
  selected.to_double() / individual_effects.length().to_double()
}

///|
pub fn uplift_auc(
  outcomes : Array[Double],
  treatment : Array[Bool],
  scores : Array[Double],
) -> Double {
  let curve = uplift_curve(outcomes, treatment, scores, scores.length())
  if curve.length() < 2 {
    return 0.0
  }
  let mut area = 0.0
  for i in 1..<curve.length() {
    area += (curve[i] + curve[i - 1]) / 2.0
  }
  area / curve.length().to_double()
}

///|
pub fn estimate_difference_in_means(
  outcomes : Array[Double],
  treatment : Array[Bool],
) -> Estimate {
  let treated = group_values(outcomes, treatment, true)
  let control = group_values(outcomes, treatment, false)
  let standard_error = if treated.length() == 0 || control.length() == 0 {
    0.0
  } else {
    (variance(treated) / treated.length().to_double() +
    variance(control) / control.length().to_double()).sqrt()
  }
  Estimate::from_standard_error(
    mean(treated) - mean(control),
    standard_error,
    outcomes.length(),
    outcomes.length().to_double(),
    "ATE (difference in means)",
  )
}
