///|
/// Long-format panel data used by production panel estimators.
pub struct PanelDataset {
  unit : Array[Int]
  time : Array[Int]
  treatment : Array[Bool]
  outcome : Array[Double]
  covariates : Array[Array[Double]]
}

///|
/// Panel data validation summary.
pub struct PanelAudit {
  observations : Int
  units : Int
  periods : Int
  duplicate_keys : Int
  missing_outcomes : Int
  balanced : Bool
  treatment_switchers : Int
  passes : Bool
}

///|
/// Fixed-effects regression result.
pub struct FixedEffectResult {
  estimate : Estimate
  coefficients : Array[Double]
  residuals : Array[Double]
  unit_effects : Array[Double]
  time_effects : Array[Double]
  within_r_squared : Double
  cluster_standard_error : Double
}

///|
/// Event-study coefficient with a relative time label.
pub struct DynamicEffect {
  relative_time : Int
  estimate : Double
  standard_error : Double
  sample_size : Int
  pre_period : Bool
  significant : Bool
}

///|
/// Synthetic-control diagnostic and prediction.
pub struct SyntheticPanelResult {
  weights : Array[Double]
  treated_prediction : Array[Double]
  treated_outcome : Array[Double]
  gap : Array[Double]
  pre_rmse : Double
  post_effect : Double
  donor_count : Int
}

///|
/// Cluster-level standard error summary.
pub struct ClusterRobustSummary {
  estimate : Double
  standard_error : Double
  cluster_count : Int
  average_cluster_size : Double
  effective_clusters : Double
}

///|
fn panel_unique_int(values : Array[Int]) -> Array[Int] {
  let result : Array[Int] = Array::new()
  for value in values {
    if !result.contains(value) {
      result.push(value)
    }
  }
  result
}

///|
fn panel_key_exists(
  units : Array[Int],
  times : Array[Int],
  unit : Int,
  time : Int,
  until : Int,
) -> Bool {
  for i in 0..<until {
    if units[i] == unit && times[i] == time {
      return true
    }
  }
  false
}

///|
/// Creates a panel dataset with aligned long-format vectors.
pub fn panel_dataset(
  unit : Array[Int],
  time : Array[Int],
  treatment : Array[Bool],
  outcome : Array[Double],
  covariates? : Array[Array[Double]] = [],
) -> PanelDataset {
  let n = unit
    .length()
    .min(time.length())
    .min(treatment.length())
    .min(outcome.length())
  {
    unit: unit[:n].to_owned(),
    time: time[:n].to_owned(),
    treatment: treatment[:n].to_owned(),
    outcome: outcome[:n].to_owned(),
    covariates: covariates[:n].to_owned(),
  }
}

///|
/// Returns unique unit identifiers.
pub fn panel_units(panel : PanelDataset) -> Array[Int] {
  panel_unique_int(panel.unit)
}

///|
/// Returns unique time identifiers in observation order.
pub fn panel_periods(panel : PanelDataset) -> Array[Int] {
  panel_unique_int(panel.time)
}

///|
/// Counts duplicate unit-time keys.
pub fn panel_duplicate_key_count(panel : PanelDataset) -> Int {
  let mut result = 0
  for i in 0..<panel.unit.length() {
    if panel_key_exists(panel.unit, panel.time, panel.unit[i], panel.time[i], i) {
      result += 1
    }
  }
  result
}

///|
/// Counts units whose treatment changes over time.
pub fn panel_switcher_count(panel : PanelDataset) -> Int {
  let units = panel_units(panel)
  let mut result = 0
  for unit in units {
    let first = Array::new()
    for i in 0..<panel.unit.length() {
      if panel.unit[i] == unit {
        first.push(panel.treatment[i])
      }
    }
    if first.length() > 1 {
      let baseline = first[0]
      let mut switched = false
      for value in first[1:] {
        if value != baseline {
          switched = true
        }
      }
      if switched {
        result += 1
      }
    }
  }
  result
}

///|
/// Checks whether every unit has every observed period.
pub fn panel_is_balanced(panel : PanelDataset) -> Bool {
  let units = panel_units(panel)
  let periods = panel_periods(panel)
  for unit in units {
    for period in periods {
      if !panel_key_exists(
          panel.unit,
          panel.time,
          unit,
          period,
          panel.unit.length(),
        ) {
        return false
      }
    }
  }
  true
}

///|
/// Audits a long-format panel.
pub fn audit_panel(panel : PanelDataset) -> PanelAudit {
  let mut missing = 0
  for value in panel.outcome {
    if !is_finite(value) {
      missing += 1
    }
  }
  let audit = {
    observations: panel.outcome.length(),
    units: panel_units(panel).length(),
    periods: panel_periods(panel).length(),
    duplicate_keys: panel_duplicate_key_count(panel),
    missing_outcomes: missing,
    balanced: panel_is_balanced(panel),
    treatment_switchers: panel_switcher_count(panel),
    passes: panel_duplicate_key_count(panel) == 0 &&
    missing == 0 &&
    panel_units(panel).length() > 1,
  }
  audit
}

///|
fn panel_mean_for_unit(panel : PanelDataset, unit : Int) -> Double {
  let values = Array::new()
  for i in 0..<panel.unit.length() {
    if panel.unit[i] == unit && is_finite(panel.outcome[i]) {
      values.push(panel.outcome[i])
    }
  }
  mean_or(values, 0.0)
}

///|
fn panel_mean_for_period(panel : PanelDataset, period : Int) -> Double {
  let values = Array::new()
  for i in 0..<panel.time.length() {
    if panel.time[i] == period && is_finite(panel.outcome[i]) {
      values.push(panel.outcome[i])
    }
  }
  mean_or(values, 0.0)
}

///|
/// Returns one-way unit means aligned to panel observations.
pub fn panel_unit_means(panel : PanelDataset) -> Array[Double] {
  let result = Array::new(capacity=panel.unit.length())
  for unit in panel.unit {
    result.push(panel_mean_for_unit(panel, unit))
  }
  result
}

///|
/// Returns one-way period means aligned to panel observations.
pub fn panel_time_means(panel : PanelDataset) -> Array[Double] {
  let result = Array::new(capacity=panel.time.length())
  for period in panel.time {
    result.push(panel_mean_for_period(panel, period))
  }
  result
}

///|
/// Applies the within transformation y_it - y_i. - y_.t + y_..
pub fn two_way_demean(panel : PanelDataset) -> Array[Double] {
  let grand = mean_or(panel.outcome, 0.0)
  let unit_means = panel_unit_means(panel)
  let time_means = panel_time_means(panel)
  let result = Array::new(capacity=panel.outcome.length())
  for i in 0..<panel.outcome.length() {
    result.push(panel.outcome[i] - unit_means[i] - time_means[i] + grand)
  }
  result
}

///|
/// Returns a within-transformed covariate matrix.
pub fn two_way_demean_matrix(
  panel : PanelDataset,
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let width = if matrix.length() == 0 { 0 } else { matrix[0].length() }
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  let unit_values : Array[Array[Double]] = Array::new(capacity=width)
  let time_values : Array[Array[Double]] = Array::new(capacity=width)
  for j in 0..<width {
    let column = Array::new(capacity=matrix.length())
    for row in matrix {
      column.push(row[j])
    }
    let temp_panel = {
      unit: panel.unit,
      time: panel.time,
      treatment: panel.treatment,
      outcome: column,
      covariates: [],
    }
    unit_values.push(panel_unit_means(temp_panel))
    time_values.push(panel_time_means(temp_panel))
  }
  let grand = mean_or(panel.outcome, 0.0)
  for i in 0..<matrix.length() {
    let row = Array::new(capacity=width)
    for j in 0..<width {
      row.push(matrix[i][j] - unit_values[j][i] - time_values[j][i] + grand)
    }
    result.push(row)
  }
  result
}

///|
/// Computes a two-period DID effect from a binary treatment group.
pub fn panel_did_effect(
  panel : PanelDataset,
  pre_period : Int,
  post_period : Int,
) -> Estimate {
  let treated_pre = Array::new()
  let treated_post = Array::new()
  let control_pre = Array::new()
  let control_post = Array::new()
  for i in 0..<panel.outcome.length() {
    if panel.time[i] == pre_period {
      if panel.treatment[i] {
        treated_pre.push(panel.outcome[i])
      } else {
        control_pre.push(panel.outcome[i])
      }
    } else if panel.time[i] == post_period {
      if panel.treatment[i] {
        treated_post.push(panel.outcome[i])
      } else {
        control_post.push(panel.outcome[i])
      }
    }
  }
  let effect = mean_or(treated_post, 0.0) -
    mean_or(treated_pre, 0.0) -
    (mean_or(control_post, 0.0) - mean_or(control_pre, 0.0))
  let treated_diff = if treated_pre.length() == 0 || treated_post.length() == 0 {
    []
  } else {
    [mean(treated_post) - mean(treated_pre)]
  }
  let control_diff = if control_pre.length() == 0 || control_post.length() == 0 {
    []
  } else {
    [mean(control_post) - mean(control_pre)]
  }
  let standard_error = if treated_diff.length() == 0 ||
    control_diff.length() == 0 {
    0.0
  } else {
    0.0
  }
  Estimate::from_standard_error(
    effect,
    standard_error,
    panel.outcome.length(),
    panel.outcome.length().to_double(),
    "panel DID",
  )
}

///|
/// Computes a simple parallel-trends diagnostic across pre-periods.
pub fn panel_parallel_trends(
  panel : PanelDataset,
  treatment_group : Bool,
  pre_periods : Array[Int],
) -> Array[Double] {
  let result : Array[Double] = Array::new(capacity=pre_periods.length())
  for period in pre_periods {
    let treated = Array::new()
    let control = Array::new()
    for i in 0..<panel.outcome.length() {
      if panel.time[i] == period {
        if panel.treatment[i] == treatment_group {
          treated.push(panel.outcome[i])
        } else {
          control.push(panel.outcome[i])
        }
      }
    }
    result.push(mean_or(treated, 0.0) - mean_or(control, 0.0))
  }
  result
}

///|
/// Fits an event-study curve relative to an event period.
pub fn panel_event_study(
  panel : PanelDataset,
  event_period : Int,
  relative_times : Array[Int],
) -> Array[DynamicEffect] {
  let result : Array[DynamicEffect] = Array::new(
    capacity=relative_times.length(),
  )
  for relative in relative_times {
    let treated = Array::new()
    let control = Array::new()
    let period = event_period + relative
    for i in 0..<panel.outcome.length() {
      if panel.time[i] == period {
        if panel.treatment[i] {
          treated.push(panel.outcome[i])
        } else {
          control.push(panel.outcome[i])
        }
      }
    }
    let estimate = mean_or(treated, 0.0) - mean_or(control, 0.0)
    let standard_error = if treated.length() == 0 || control.length() == 0 {
      0.0
    } else {
      (variance(treated) / treated.length().to_double() +
      variance(control) / control.length().to_double()).sqrt()
    }
    result.push({
      relative_time: relative,
      estimate,
      standard_error,
      sample_size: treated.length() + control.length(),
      pre_period: relative < 0,
      significant: standard_error > 0.0 &&
      estimate.abs() > 1.96 * standard_error,
    })
  }
  result
}

///|
/// Computes the average post-event dynamic effect.
pub fn panel_average_post_effect(effects : Array[DynamicEffect]) -> Double {
  let values = Array::new()
  for effect in effects {
    if !effect.pre_period {
      values.push(effect.estimate)
    }
  }
  mean_or(values, 0.0)
}

///|
/// Computes the largest absolute pre-trend coefficient.
pub fn panel_max_pretrend(effects : Array[DynamicEffect]) -> Double {
  let mut result = 0.0
  for effect in effects {
    if effect.pre_period && effect.estimate.abs() > result {
      result = effect.estimate.abs()
    }
  }
  result
}

///|
/// Returns placebo event-study effects for supplied pseudo-event periods.
pub fn panel_placebo_effects(
  panel : PanelDataset,
  placebo_periods : Array[Int],
  pre_period : Int,
) -> Array[Double] {
  let result = Array::new(capacity=placebo_periods.length())
  for placebo in placebo_periods {
    result.push(panel_did_effect(panel, pre_period, placebo).estimate)
  }
  result
}

///|
/// Calculates a cluster-robust standard error for unit-level effect contributions.
pub fn panel_cluster_robust(
  estimate : Double,
  cluster_ids : Array[Int],
  contributions : Array[Double],
) -> ClusterRobustSummary {
  let n = cluster_ids.length().min(contributions.length())
  let clusters = panel_unique_int(cluster_ids[:n].to_owned())
  if n == 0 || clusters.length() == 0 {
    return {
      estimate,
      standard_error: 0.0,
      cluster_count: 0,
      average_cluster_size: 0.0,
      effective_clusters: 0.0,
    }
  }
  let mut meat = 0.0
  let mut size_sum = 0
  for cluster in clusters {
    let mut aggregate = 0.0
    let mut size = 0
    for i in 0..<n {
      if cluster_ids[i] == cluster {
        aggregate += contributions[i]
        size += 1
      }
    }
    meat += aggregate * aggregate
    size_sum += size
  }
  let variance_value = meat / @math.pow(n.to_double(), 2.0)
  let standard_error = variance_value.max(0.0).sqrt()
  let cluster_count = clusters.length()
  {
    estimate,
    standard_error,
    cluster_count,
    average_cluster_size: size_sum.to_double() / cluster_count.to_double(),
    effective_clusters: if meat == 0.0 {
      0.0
    } else {
      clusters.length().to_double() *
      estimate.abs().max(1.0) /
      (1.0 + standard_error)
    },
  }
}

///|
/// Estimates a unit fixed effect by subtracting each unit's baseline period.
pub fn unit_fixed_effect(
  panel : PanelDataset,
  baseline_period : Int,
) -> Array[Double] {
  let baseline : Array[(Int, Double)] = Array::new()
  for unit in panel_units(panel) {
    let mut value = 0.0
    for i in 0..<panel.outcome.length() {
      if panel.unit[i] == unit && panel.time[i] == baseline_period {
        value = panel.outcome[i]
      }
    }
    baseline.push((unit, value))
  }
  let result = Array::new(capacity=panel.outcome.length())
  for i in 0..<panel.outcome.length() {
    let mut value = panel.outcome[i]
    for item in baseline {
      if item.0 == panel.unit[i] {
        value -= item.1
      }
    }
    result.push(value)
  }
  result
}

///|
/// Calculates a synthetic-control prediction from donor outcomes and weights.
pub fn synthetic_panel(
  treated_outcome : Array[Double],
  donor_outcomes : Array[Array[Double]],
  weights : Array[Double],
  intervention_period : Int,
) -> SyntheticPanelResult {
  let periods = treated_outcome.length()
  let prediction = Array::make(periods, 0.0)
  let donor_count = donor_outcomes.length().min(weights.length())
  for donor in 0..<donor_count {
    for time in 0..<periods.min(donor_outcomes[donor].length()) {
      prediction[time] += weights[donor] * donor_outcomes[donor][time]
    }
  }
  let gap = Array::new(capacity=periods)
  for time in 0..<periods {
    gap.push(treated_outcome[time] - prediction[time])
  }
  let pre = gap[:periods.min(intervention_period)].to_owned()
  let post = if intervention_period < periods {
    gap[intervention_period:].to_owned()
  } else {
    []
  }
  let post_effect = mean_or(post, 0.0)
  {
    weights,
    treated_prediction: prediction,
    treated_outcome,
    gap,
    pre_rmse: optimized_vector_rmse(pre),
    post_effect,
    donor_count,
  }
}

///|
fn optimized_vector_rmse(values : Array[Double]) -> Double {
  if values.length() == 0 {
    0.0
  } else {
    (sum(values.map(fn(value) { value * value })) / values.length().to_double()).sqrt()
  }
}

///|
/// Fits a within-transformed linear outcome model for a panel.
pub fn two_way_fixed_effect(
  panel : PanelDataset,
  covariates : Array[Array[Double]],
) -> FixedEffectResult {
  let demeaned_outcome = two_way_demean(panel)
  let demeaned_covariates = two_way_demean_matrix(panel, covariates)
  let model = fit_linear_outcome_model(demeaned_covariates, demeaned_outcome)
  let predictions = predict_outcomes(model, demeaned_covariates)
  let residuals = Array::new(capacity=demeaned_outcome.length())
  for i in 0..<demeaned_outcome.length().min(predictions.length()) {
    residuals.push(demeaned_outcome[i] - predictions[i])
  }
  let residual_variance = variance(residuals)
  let standard_error = if residuals.length() == 0 {
    0.0
  } else {
    (residual_variance / residuals.length().to_double()).sqrt()
  }
  let estimate = Estimate::from_standard_error(
    model.intercept,
    standard_error,
    residuals.length(),
    residuals.length().to_double(),
    "two-way fixed effects",
  )
  let average = mean_or(demeaned_outcome, 0.0)
  let mut total = 0.0
  for value in demeaned_outcome {
    total += (value - average) * (value - average)
  }
  let mut unexplained = 0.0
  for value in residuals {
    unexplained += value * value
  }
  let within_r_squared = if total == 0.0 {
    0.0
  } else {
    1.0 - unexplained / total
  }
  {
    estimate,
    coefficients: model.coefficients,
    residuals,
    unit_effects: panel_unit_means(panel),
    time_effects: panel_time_means(panel),
    within_r_squared,
    cluster_standard_error: standard_error,
  }
}

///|
/// Normalizes non-negative synthetic-control weights.
pub fn normalize_donor_weights(weights : Array[Double]) -> Array[Double] {
  let result = Array::new(capacity=weights.length())
  let mut total = 0.0
  for weight in weights {
    total += weight.max(0.0)
  }
  for weight in weights {
    result.push(if total == 0.0 { 0.0 } else { weight.max(0.0) / total })
  }
  result
}

///|
/// Computes pre-treatment fit for a donor matrix and target trajectory.
pub fn donor_pre_rmse(
  target : Array[Double],
  donors : Array[Array[Double]],
  weights : Array[Double],
  pre_periods : Int,
) -> Double {
  let prediction = Array::make(pre_periods.min(target.length()), 0.0)
  for donor in 0..<donors.length().min(weights.length()) {
    for i in 0..<prediction.length().min(donors[donor].length()) {
      prediction[i] += weights[donor] * donors[donor][i]
    }
  }
  if prediction.length() == 0 {
    return 0.0
  }
  let mut total = 0.0
  for i in 0..<prediction.length() {
    let difference = target[i] - prediction[i]
    total += difference * difference
  }
  (total / prediction.length().to_double()).sqrt()
}

///|
/// Returns lagged outcome changes within each unit.
pub fn panel_first_difference(panel : PanelDataset) -> Array[Double] {
  let result = Array::make(panel.outcome.length(), 0.0)
  for i in 0..<panel.outcome.length() {
    let mut found = false
    let mut previous_time = -2147483647
    let mut previous_outcome = 0.0
    for j in 0..<panel.outcome.length() {
      if panel.unit[j] == panel.unit[i] &&
        panel.time[j] < panel.time[i] &&
        panel.time[j] > previous_time {
        previous_time = panel.time[j]
        previous_outcome = panel.outcome[j]
        found = true
      }
    }
    if found {
      result[i] = panel.outcome[i] - previous_outcome
    }
  }
  result
}

///|
/// Computes average outcome by integer time period.
pub fn panel_period_means(panel : PanelDataset) -> Array[Array[Double]] {
  let periods = panel_periods(panel)
  let result : Array[Array[Double]] = Array::new(capacity=periods.length())
  for period in periods {
    result.push([period.to_double(), panel_mean_for_period(panel, period)])
  }
  result
}

///|
/// Returns a compact panel diagnostic vector.
pub fn panel_summary_vector(panel : PanelDataset) -> Array[Double] {
  let audit = audit_panel(panel)
  [
    audit.observations.to_double(),
    audit.units.to_double(),
    audit.periods.to_double(),
    audit.duplicate_keys.to_double(),
    audit.missing_outcomes.to_double(),
    if audit.balanced {
      1.0
    } else {
      0.0
    },
    audit.treatment_switchers.to_double(),
    if audit.passes {
      1.0
    } else {
      0.0
    },
  ]
}
