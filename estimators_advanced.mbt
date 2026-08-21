///|
/// Generalized effect estimate with diagnostics and uncertainty.
pub struct AdvancedEffect {
  estimate : Double
  standard_error : Double
  lower : Double
  upper : Double
  effective_sample_size : Double
  estimand : String
  passes : Bool
}

///|
/// Augmented inverse-probability estimate.
pub struct AugmentedEffect {
  estimate : Double
  standard_error : Double
  influence : Array[Double]
  propensity_minimum : Double
  propensity_maximum : Double
  overlap_fraction : Double
  passes : Bool
}

///|
/// Instrumental-variable Wald estimate.
pub struct InstrumentalVariableResult {
  estimate : Double
  standard_error : Double
  first_stage_difference : Double
  reduced_form_difference : Double
  compliance_rate : Double
  weak_instrument : Bool
  passes : Bool
}

///|
/// Local regression discontinuity estimate.
pub struct RegressionDiscontinuityResult {
  estimate : Double
  standard_error : Double
  bandwidth : Double
  left_sample_size : Int
  right_sample_size : Int
  left_slope : Double
  right_slope : Double
  passes : Bool
}

///|
/// Cross-fitted nuisance prediction summary.
pub struct CrossFitSummary {
  folds : Int
  propensity_predictions : Array[Double]
  outcome_treated_predictions : Array[Double]
  outcome_control_predictions : Array[Double]
  fold_sizes : Array[Int]
  prediction_rmse : Double
}

///|
/// Effect curve point indexed by a dose or score threshold.
pub struct EffectCurvePoint {
  threshold : Double
  treated_count : Int
  control_count : Int
  effect : Double
  standard_error : Double
  policy_value : Double
}

///|
/// Bounds under an unmeasured confounding multiplier.
pub struct ConfoundingBounds {
  strength : Double
  lower : Double
  upper : Double
  center : Double
  width : Double
}

///|
fn ea_target(value : Bool) -> Double {
  if value {
    1.0
  } else {
    0.0
  }
}

///|
fn ea_normal_quantile(probability : Double) -> Double {
  let p = clamp(probability, 1.0e-7, 1.0 - 1.0e-7)
  let mut low = -9.0
  let mut high = 9.0
  for _ in 0..<70 {
    let middle = (low + high) / 2.0
    let sign = if middle < 0.0 { -1.0 } else { 1.0 }
    let x = middle.abs()
    let t = 1.0 / (1.0 + 0.2316419 * x)
    let polynomial = t *
      (
        0.319381530 +
        t *
        (
          -0.356563782 +
          t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))
        )
      )
    let cdf = 0.5 +
      sign * (0.5 - @math.exp(-0.5 * x * x) / 2.5066282746310002 * polynomial)
    if cdf < p {
      low = middle
    } else {
      high = middle
    }
  }
  (low + high) / 2.0
}

///|
/// Calculates overlap weights for a binary treatment.
pub fn advanced_overlap_weights(
  treatment : Array[Bool],
  propensity : Array[Double],
) -> Array[Double] {
  let n = treatment.length().min(propensity.length())
  let result = Array::new(capacity=n)
  for i in 0..<n {
    let probability = safe_probability(propensity[i])
    result.push(if treatment[i] { 1.0 - probability } else { probability })
  }
  result
}

///|
/// Calculates stabilized inverse-probability weights with marginal treatment probability.
pub fn advanced_stabilized_weights(
  treatment : Array[Bool],
  propensity : Array[Double],
) -> Array[Double] {
  let n = treatment.length().min(propensity.length())
  let mut treated = 0
  for i in 0..<n {
    if treatment[i] {
      treated += 1
    }
  }
  let marginal = if n == 0 { 0.5 } else { treated.to_double() / n.to_double() }
  let result = Array::new(capacity=n)
  for i in 0..<n {
    let p = safe_probability(propensity[i])
    result.push(
      if treatment[i] {
        marginal / p
      } else {
        (1.0 - marginal) / (1.0 - p)
      },
    )
  }
  result
}

///|
/// Returns a normalized weighted mean for one treatment arm.
pub fn advanced_weighted_arm_mean(
  outcomes : Array[Double],
  treatment : Array[Bool],
  weights : Array[Double],
  arm : Bool,
) -> Double {
  let n = outcomes.length().min(treatment.length()).min(weights.length())
  let mut numerator = 0.0
  let mut denominator = 0.0
  for i in 0..<n {
    if treatment[i] == arm {
      numerator += weights[i] * outcomes[i]
      denominator += weights[i]
    }
  }
  if denominator == 0.0 {
    0.0
  } else {
    numerator / denominator
  }
}

///|
/// Computes a weighted ATE and sandwich-like standard error.
pub fn advanced_weighted_ate(
  outcomes : Array[Double],
  treatment : Array[Bool],
  weights : Array[Double],
) -> AdvancedEffect {
  let n = outcomes.length().min(treatment.length()).min(weights.length())
  let treated_mean = advanced_weighted_arm_mean(
    outcomes, treatment, weights, true,
  )
  let control_mean = advanced_weighted_arm_mean(
    outcomes, treatment, weights, false,
  )
  let estimate = treated_mean - control_mean
  let treated_residuals = Array::new()
  let control_residuals = Array::new()
  let treated_weights = Array::new()
  let control_weights = Array::new()
  for i in 0..<n {
    if treatment[i] {
      treated_residuals.push(outcomes[i] - treated_mean)
      treated_weights.push(weights[i])
    } else {
      control_residuals.push(outcomes[i] - control_mean)
      control_weights.push(weights[i])
    }
  }
  let treated_variance = weighted_variance(treated_residuals, treated_weights)
  let control_variance = weighted_variance(control_residuals, control_weights)
  let treated_ess = effective_sample_size(treated_weights)
  let control_ess = effective_sample_size(control_weights)
  let treated_component = if treated_ess == 0.0 {
    0.0
  } else {
    treated_variance / treated_ess
  }
  let control_component = if control_ess == 0.0 {
    0.0
  } else {
    control_variance / control_ess
  }
  let standard_error = (treated_component + control_component).sqrt()
  let critical = ea_normal_quantile(0.975)
  {
    estimate,
    standard_error,
    lower: estimate - critical * standard_error,
    upper: estimate + critical * standard_error,
    effective_sample_size: treated_ess + control_ess,
    estimand: "weighted ATE",
    passes: n > 1 && treated_ess > 2.0 && control_ess > 2.0,
  }
}

///|
/// Computes an augmented inverse-probability ATE from nuisance predictions.
pub fn augmented_inverse_probability_effect(
  outcomes : Array[Double],
  treatment : Array[Bool],
  propensity : Array[Double],
  predicted_treated : Array[Double],
  predicted_control : Array[Double],
) -> AugmentedEffect {
  let n = outcomes
    .length()
    .min(treatment.length())
    .min(propensity.length())
    .min(predicted_treated.length())
    .min(predicted_control.length())
  if n == 0 {
    return {
      estimate: 0.0,
      standard_error: 0.0,
      influence: [],
      propensity_minimum: 0.0,
      propensity_maximum: 0.0,
      overlap_fraction: 0.0,
      passes: false,
    }
  }
  let influence = Array::new(capacity=n)
  let mut total = 0.0
  let mut minimum = 1.0
  let mut maximum = 0.0
  let mut overlap = 0
  for i in 0..<n {
    let p = safe_probability(propensity[i])
    let target = ea_target(treatment[i])
    let value = predicted_treated[i] -
      predicted_control[i] +
      target * (outcomes[i] - predicted_treated[i]) / p -
      (1.0 - target) * (outcomes[i] - predicted_control[i]) / (1.0 - p)
    influence.push(value)
    total += value
    if p < minimum {
      minimum = p
    }
    if p > maximum {
      maximum = p
    }
    if p > 0.05 && p < 0.95 {
      overlap += 1
    }
  }
  let estimate = total / n.to_double()
  let centered = Array::new(capacity=n)
  for value in influence {
    centered.push(value - estimate)
  }
  let standard_error = std_dev(centered) / n.to_double().sqrt()
  {
    estimate,
    standard_error,
    influence: centered,
    propensity_minimum: minimum,
    propensity_maximum: maximum,
    overlap_fraction: overlap.to_double() / n.to_double(),
    passes: overlap > 0 && is_finite(estimate) && is_finite(standard_error),
  }
}

///|
/// Computes the average treatment effect on the treated with augmentation.
pub fn augmented_att(
  outcomes : Array[Double],
  treatment : Array[Bool],
  propensity : Array[Double],
  predicted_control : Array[Double],
) -> AdvancedEffect {
  let n = outcomes
    .length()
    .min(treatment.length())
    .min(propensity.length())
    .min(predicted_control.length())
  let treated_outcomes = Array::new()
  let treated_residuals = Array::new()
  for i in 0..<n {
    if treatment[i] {
      let p = safe_probability(propensity[i])
      treated_outcomes.push(outcomes[i] - predicted_control[i])
      treated_residuals.push((outcomes[i] - predicted_control[i]) / p)
    }
  }
  let estimate = mean_or(treated_outcomes, 0.0)
  let standard_error = if treated_residuals.length() < 2 {
    0.0
  } else {
    std_dev(treated_residuals) / treated_residuals.length().to_double().sqrt()
  }
  let critical = ea_normal_quantile(0.975)
  {
    estimate,
    standard_error,
    lower: estimate - critical * standard_error,
    upper: estimate + critical * standard_error,
    effective_sample_size: treated_residuals.length().to_double(),
    estimand: "augmented ATT",
    passes: treated_residuals.length() > 2,
  }
}

///|
/// Computes a Wald instrumental-variable estimate from binary instrument and treatment.
pub fn instrumental_variable_wald(
  outcome : Array[Double],
  treatment : Array[Bool],
  instrument : Array[Bool],
) -> InstrumentalVariableResult {
  let n = outcome.length().min(treatment.length()).min(instrument.length())
  let treated_instrument = Array::new()
  let control_instrument = Array::new()
  let treated_outcome = Array::new()
  let control_outcome = Array::new()
  for i in 0..<n {
    if instrument[i] {
      treated_instrument.push(if treatment[i] { 1.0 } else { 0.0 })
      treated_outcome.push(outcome[i])
    } else {
      control_instrument.push(if treatment[i] { 1.0 } else { 0.0 })
      control_outcome.push(outcome[i])
    }
  }
  let first_stage = mean_or(treated_instrument, 0.0) -
    mean_or(control_instrument, 0.0)
  let reduced_form = mean_or(treated_outcome, 0.0) -
    mean_or(control_outcome, 0.0)
  let estimate = if first_stage.abs() < 1.0e-8 {
    0.0
  } else {
    reduced_form / first_stage
  }
  let weak = first_stage.abs() < 0.1
  let standard_error = if n == 0 || weak {
    0.0
  } else {
    std_dev(outcome) / n.to_double().sqrt() / first_stage.abs()
  }
  {
    estimate,
    standard_error,
    first_stage_difference: first_stage,
    reduced_form_difference: reduced_form,
    compliance_rate: first_stage,
    weak_instrument: weak,
    passes: n > 10 && !weak,
  }
}

///|
/// Computes a two-stage linear IV estimate using a scalar endogenous regressor.
pub fn two_stage_least_squares(
  outcome : Array[Double],
  treatment : Array[Double],
  instrument : Array[Double],
) -> InstrumentalVariableResult {
  let n = outcome.length().min(treatment.length()).min(instrument.length())
  if n < 2 {
    return {
      estimate: 0.0,
      standard_error: 0.0,
      first_stage_difference: 0.0,
      reduced_form_difference: 0.0,
      compliance_rate: 0.0,
      weak_instrument: true,
      passes: false,
    }
  }
  let first_slope = covariance(
      instrument[:n].to_owned(),
      treatment[:n].to_owned(),
    ) /
    variance(instrument[:n].to_owned())
  let reduced_slope = covariance(
      instrument[:n].to_owned(),
      outcome[:n].to_owned(),
    ) /
    variance(instrument[:n].to_owned())
  let estimate = if first_slope.abs() < 1.0e-8 {
    0.0
  } else {
    reduced_slope / first_slope
  }
  let fitted_treatment = Array::new(capacity=n)
  for value in instrument[:n] {
    fitted_treatment.push(
      mean(treatment[:n].to_owned()) +
      first_slope * (value - mean(instrument[:n].to_owned())),
    )
  }
  let residuals = Array::new(capacity=n)
  for i in 0..<n {
    residuals.push(outcome[i] - estimate * fitted_treatment[i])
  }
  let standard_error = std_dev(residuals) /
    n.to_double().sqrt() /
    first_slope.abs().max(1.0e-8)
  {
    estimate,
    standard_error,
    first_stage_difference: first_slope,
    reduced_form_difference: reduced_slope,
    compliance_rate: first_slope,
    weak_instrument: first_slope.abs() < 0.1,
    passes: first_slope.abs() >= 0.1,
  }
}

///|
/// Fits a local linear regression discontinuity estimate using a rectangular kernel.
pub fn regression_discontinuity(
  running : Array[Double],
  outcome : Array[Double],
  cutoff : Double,
  bandwidth : Double,
) -> RegressionDiscontinuityResult {
  let n = running.length().min(outcome.length())
  let width = bandwidth.abs()
  let left_x = Array::new()
  let left_y = Array::new()
  let right_x = Array::new()
  let right_y = Array::new()
  for i in 0..<n {
    let distance = running[i] - cutoff
    if distance < 0.0 && distance.abs() <= width {
      left_x.push(distance)
      left_y.push(outcome[i])
    } else if distance >= 0.0 && distance <= width {
      right_x.push(distance)
      right_y.push(outcome[i])
    }
  }
  let left_slope = covariance(left_x, left_y) / variance(left_x)
  let right_slope = covariance(right_x, right_y) / variance(right_x)
  let left_intercept = mean_or(left_y, 0.0) -
    (if is_finite(left_slope) { left_slope } else { 0.0 }) *
    mean_or(left_x, 0.0)
  let right_intercept = mean_or(right_y, 0.0) -
    (if is_finite(right_slope) { right_slope } else { 0.0 }) *
    mean_or(right_x, 0.0)
  let estimate = right_intercept - left_intercept
  let residuals = Array::new()
  for i in 0..<left_x.length() {
    residuals.push(left_y[i] - left_intercept - left_slope * left_x[i])
  }
  for i in 0..<right_x.length() {
    residuals.push(right_y[i] - right_intercept - right_slope * right_x[i])
  }
  let standard_error = if residuals.length() < 2 {
    0.0
  } else {
    std_dev(residuals) *
    (1.0 / left_x.length().max(1).to_double() +
    1.0 / right_x.length().max(1).to_double()).sqrt()
  }
  {
    estimate,
    standard_error,
    bandwidth: width,
    left_sample_size: left_x.length(),
    right_sample_size: right_x.length(),
    left_slope: if is_finite(left_slope) {
      left_slope
    } else {
      0.0
    },
    right_slope: if is_finite(right_slope) {
      right_slope
    } else {
      0.0
    },
    passes: left_x.length() >= 10 && right_x.length() >= 10,
  }
}

///|
/// Evaluates a regression-discontinuity estimate over a bandwidth grid.
pub fn regression_discontinuity_sensitivity(
  running : Array[Double],
  outcome : Array[Double],
  cutoff : Double,
  bandwidths : Array[Double],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=bandwidths.length())
  for bandwidth in bandwidths {
    let fit = regression_discontinuity(running, outcome, cutoff, bandwidth)
    result.push([
      bandwidth,
      fit.estimate,
      fit.standard_error,
      fit.left_sample_size.to_double(),
      fit.right_sample_size.to_double(),
    ])
  }
  result
}

///|
/// Computes an individualized policy value from outcome and treatment predictions.
pub fn advanced_policy_value(
  treatment : Array[Bool],
  outcomes : Array[Double],
  propensity : Array[Double],
  policy_treatment : Array[Bool],
) -> AdvancedEffect {
  let n = treatment
    .length()
    .min(outcomes.length())
    .min(propensity.length())
    .min(policy_treatment.length())
  let contributions = Array::new(capacity=n)
  for i in 0..<n {
    let p = safe_probability(propensity[i])
    let denominator = if treatment[i] { p } else { 1.0 - p }
    let observed = if treatment[i] == policy_treatment[i] {
      outcomes[i] / denominator
    } else {
      0.0
    }
    contributions.push(observed)
  }
  let estimate = mean_or(contributions, 0.0)
  let standard_error = if contributions.length() < 2 {
    0.0
  } else {
    std_dev(contributions) / contributions.length().to_double().sqrt()
  }
  let interval = normal_confidence_interval(estimate, standard_error)
  {
    estimate,
    standard_error,
    lower: interval[0],
    upper: interval[1],
    effective_sample_size: contributions.length().to_double(),
    estimand: "policy value",
    passes: contributions.length() > 2,
  }
}

///|
/// Computes policy value at a score threshold.
pub fn threshold_policy_curve(
  treatment : Array[Bool],
  outcomes : Array[Double],
  propensity : Array[Double],
  scores : Array[Double],
  thresholds : Array[Double],
) -> Array[EffectCurvePoint] {
  let n = treatment
    .length()
    .min(outcomes.length())
    .min(propensity.length())
    .min(scores.length())
  let result : Array[EffectCurvePoint] = Array::new(
    capacity=thresholds.length(),
  )
  for threshold in thresholds {
    let treated_outcomes = Array::new()
    let control_outcomes = Array::new()
    let mut policy_total = 0.0
    for i in 0..<n {
      if scores[i] >= threshold {
        treated_outcomes.push(outcomes[i])
        policy_total += outcomes[i]
      } else {
        control_outcomes.push(outcomes[i])
        policy_total += outcomes[i]
      }
    }
    let effect = mean_or(treated_outcomes, 0.0) - mean_or(control_outcomes, 0.0)
    let standard_error = if treated_outcomes.length() == 0 ||
      control_outcomes.length() == 0 {
      0.0
    } else {
      (variance(treated_outcomes) / treated_outcomes.length().to_double() +
      variance(control_outcomes) / control_outcomes.length().to_double()).sqrt()
    }
    result.push({
      threshold,
      treated_count: treated_outcomes.length(),
      control_count: control_outcomes.length(),
      effect,
      standard_error,
      policy_value: if n == 0 {
        0.0
      } else {
        policy_total / n.to_double()
      },
    })
  }
  result
}

///|
/// Computes sensitivity bounds by scaling an observed standard error.
pub fn confounding_bounds(
  estimate : Double,
  standard_error : Double,
  strengths : Array[Double],
) -> Array[ConfoundingBounds] {
  let result : Array[ConfoundingBounds] = Array::new(
    capacity=strengths.length(),
  )
  for strength in strengths {
    let multiplier = strength.abs().max(1.0)
    let width = multiplier * standard_error.abs()
    result.push({
      strength,
      lower: estimate - width,
      upper: estimate + width,
      center: estimate,
      width: 2.0 * width,
    })
  }
  result
}

///|
/// Computes a crude partial-R2 robustness value for a standardized effect.
pub fn robustness_value(
  estimate : Double,
  standard_error : Double,
  sample_size : Int,
) -> Double {
  if sample_size <= 2 || standard_error == 0.0 {
    return 0.0
  }
  let t = estimate / standard_error
  let numerator = t * t - 1.0
  if numerator <= 0.0 {
    0.0
  } else {
    (numerator / sample_size.to_double()).sqrt()
  }
}

///|
/// Builds a cross-fitting fold assignment with deterministic rotation.
pub fn cross_fit_folds(
  sample_size : Int,
  folds : Int,
  seed? : UInt64 = 20260819,
) -> Array[Int] {
  let n = if sample_size > 0 { sample_size } else { 0 }
  let count = if folds > 1 { folds } else { 2 }
  let order = shuffled_indices(n, seed)
  let result = Array::make(n, 0)
  for position in 0..<n {
    result[order[position]] = position % count
  }
  result
}

///|
/// Aggregates out-of-fold nuisance predictions and validation error.
pub fn cross_fit_summary(
  folds : Array[Int],
  propensity : Array[Double],
  outcome_treated : Array[Double],
  outcome_control : Array[Double],
  treatment : Array[Bool],
  outcomes : Array[Double],
) -> CrossFitSummary {
  let n = folds
    .length()
    .min(propensity.length())
    .min(outcome_treated.length())
    .min(outcome_control.length())
    .min(treatment.length())
    .min(outcomes.length())
  let fold_sizes = Array::new()
  let mut maximum_fold = -1
  for value in folds[:n] {
    if value > maximum_fold {
      maximum_fold = value
    }
  }
  for fold in 0..<=maximum_fold {
    let mut count = 0
    for value in folds[:n] {
      if value == fold {
        count += 1
      }
    }
    fold_sizes.push(count)
  }
  let mut squared_error = 0.0
  for i in 0..<n {
    let prediction = if treatment[i] {
      outcome_treated[i]
    } else {
      outcome_control[i]
    }
    let error = prediction - outcomes[i]
    squared_error += error * error
  }
  {
    folds: fold_sizes.length(),
    propensity_predictions: propensity[:n].to_owned(),
    outcome_treated_predictions: outcome_treated[:n].to_owned(),
    outcome_control_predictions: outcome_control[:n].to_owned(),
    fold_sizes,
    prediction_rmse: if n == 0 {
      0.0
    } else {
      (squared_error / n.to_double()).sqrt()
    },
  }
}

///|
/// Returns a compact advanced-estimator summary vector.
pub fn advanced_effect_summary(effect : AdvancedEffect) -> Array[Double] {
  [
    effect.estimate,
    effect.standard_error,
    effect.lower,
    effect.upper,
    effect.effective_sample_size,
    if effect.passes {
      1.0
    } else {
      0.0
    },
  ]
}
