///|
/// A regularly ordered numeric time series.
pub struct CausalTimeSeries {
  time : Array[Int]
  values : Array[Double]
  frequency : Int
  name : String
}

///|
/// Time-series diagnostics used before an interrupted-series analysis.
pub struct TimeSeriesProfile {
  observations : Int
  missing : Int
  mean : Double
  standard_deviation : Double
  first_difference_mean : Double
  autocorrelation_lag1 : Double
  trend_slope : Double
  seasonal_period : Int
  passes : Bool
}

///|
/// Forecast result with residual scale and horizon labels.
pub struct ForecastResult {
  time : Array[Int]
  forecast : Array[Double]
  lower : Array[Double]
  upper : Array[Double]
  residual_standard_error : Double
  model_order : Int
}

///|
/// Interrupted time-series estimate.
pub struct InterruptedEffect {
  level_change : Double
  slope_change : Double
  standard_error : Double
  pre_slope : Double
  post_slope : Double
  pre_observations : Int
  post_observations : Int
  passes : Bool
}

///|
/// Candidate change point diagnostic.
pub struct ChangePoint {
  index : Int
  statistic : Double
  left_mean : Double
  right_mean : Double
  detected : Bool
}

///|
/// Newey-West long-run variance summary.
pub struct LongRunVariance {
  variance : Double
  standard_error : Double
  bandwidth : Int
  autocovariances : Array[Double]
}

///|
fn ts_mean(values : Array[Double]) -> Double {
  mean_or(values, 0.0)
}

///|
fn ts_time_mean(time : Array[Int]) -> Double {
  if time.length() == 0 {
    0.0
  } else {
    let mut total = 0.0
    for value in time {
      total += value.to_double()
    }
    total / time.length().to_double()
  }
}

///|
/// Creates a validated time series.
pub fn causal_time_series(
  time : Array[Int],
  values : Array[Double],
  frequency? : Int = 1,
  name? : String = "outcome",
) -> CausalTimeSeries {
  let n = time.length().min(values.length())
  {
    time: time[:n].to_owned(),
    values: values[:n].to_owned(),
    frequency: if frequency > 0 {
      frequency
    } else {
      1
    },
    name,
  }
}

///|
/// Returns a profile for a time series.
pub fn time_series_profile(
  series : CausalTimeSeries,
  seasonal_period? : Int = 1,
) -> TimeSeriesProfile {
  let missing = missing_count(series.values)
  let observed = series.values.filter(fn(value) { is_finite(value) })
  let differences = first_difference(series.values)
  let slope = linear_trend_slope(series.time, series.values)
  let autocorrelation = autocorrelation_at(series.values, 1)
  {
    observations: series.values.length(),
    missing,
    mean: ts_mean(observed),
    standard_deviation: std_dev(observed),
    first_difference_mean: ts_mean(differences),
    autocorrelation_lag1: autocorrelation,
    trend_slope: slope,
    seasonal_period: if seasonal_period > 0 {
      seasonal_period
    } else {
      1
    },
    passes: missing == 0 && series.values.length() >= 4,
  }
}

///|
/// Returns a lagged vector aligned with the original index.
pub fn lag_series(
  values : Array[Double],
  lag : Int,
  fill? : Double = 0.0,
) -> Array[Double] {
  let result = Array::make(values.length(), fill)
  if lag <= 0 {
    return values.copy()
  }
  for i in lag..<values.length() {
    result[i] = values[i - lag]
  }
  result
}

///|
/// Returns a leading vector aligned with the original index.
pub fn lead_series(
  values : Array[Double],
  lead : Int,
  fill? : Double = 0.0,
) -> Array[Double] {
  let result = Array::make(values.length(), fill)
  if lead <= 0 {
    return values.copy()
  }
  if lead < values.length() {
    for i in 0..<(values.length() - lead) {
      result[i] = values[i + lead]
    }
  }
  result
}

///|
/// Computes first differences with a configurable initial value.
pub fn first_difference(
  values : Array[Double],
  initial? : Double = 0.0,
) -> Array[Double] {
  let result = Array::new(capacity=values.length())
  if values.length() == 0 {
    return result
  }
  result.push(values[0] - initial)
  for i in 1..<values.length() {
    result.push(values[i] - values[i - 1])
  }
  result
}

///|
/// Computes a k-step difference.
pub fn k_difference(values : Array[Double], lag : Int) -> Array[Double] {
  if lag <= 0 {
    return values.copy()
  }
  let result = Array::new(capacity=values.length())
  for i in lag..<values.length() {
    result.push(values[i] - values[i - lag])
  }
  result
}

///|
/// Computes a trailing moving average with partial windows at the boundary.
pub fn moving_average(values : Array[Double], window : Int) -> Array[Double] {
  let result = Array::new(capacity=values.length())
  if values.length() == 0 {
    return result
  }
  let width = if window > 0 { window } else { 1 }
  for i in 0..<values.length() {
    let start = if i + 1 > width { i + 1 - width } else { 0 }
    let end = i + 1
    let slice = values[start:end].to_owned()
    result.push(ts_mean(slice))
  }
  result
}

///|
/// Computes a trailing moving sum.
pub fn moving_sum(values : Array[Double], window : Int) -> Array[Double] {
  let result = Array::new(capacity=values.length())
  let width = if window > 0 { window } else { 1 }
  let mut running = 0.0
  for i in 0..<values.length() {
    running += values[i]
    if i >= width {
      running -= values[i - width]
    }
    result.push(running)
  }
  result
}

///|
/// Computes exponentially weighted moving averages.
pub fn exponential_moving_average(
  values : Array[Double],
  alpha : Double,
) -> Array[Double] {
  let result = Array::new(capacity=values.length())
  if values.length() == 0 {
    return result
  }
  let weight = clamp(alpha, 1.0e-6, 1.0)
  let mut level = values[0]
  result.push(level)
  for value in values[1:] {
    level = weight * value + (1.0 - weight) * level
    result.push(level)
  }
  result
}

///|
/// Computes rolling standard deviations.
pub fn rolling_standard_deviation(
  values : Array[Double],
  window : Int,
) -> Array[Double] {
  let result = Array::new(capacity=values.length())
  let width = if window > 0 { window } else { 1 }
  for i in 0..<values.length() {
    let start = if i + 1 > width { i + 1 - width } else { 0 }
    let end = i + 1
    result.push(std_dev(values[start:end].to_owned()))
  }
  result
}

///|
/// Computes the lag-k sample autocorrelation.
pub fn autocorrelation_at(values : Array[Double], lag : Int) -> Double {
  if lag <= 0 || values.length() <= lag {
    return 0.0
  }
  let average = mean(values)
  let mut numerator = 0.0
  let mut denominator = 0.0
  for value in values {
    let difference = value - average
    denominator += difference * difference
  }
  for i in lag..<values.length() {
    numerator += (values[i] - average) * (values[i - lag] - average)
  }
  if denominator == 0.0 {
    0.0
  } else {
    numerator / denominator
  }
}

///|
/// Computes autocorrelations from lag one through a maximum lag.
pub fn autocorrelation_curve(
  values : Array[Double],
  maximum_lag : Int,
) -> Array[Double] {
  let result = Array::new(
    capacity=if maximum_lag > 0 { maximum_lag } else { 0 },
  )
  for lag in 1..<=maximum_lag {
    result.push(autocorrelation_at(values, lag))
  }
  result
}

///|
/// Computes a least-squares linear trend slope.
pub fn linear_trend_slope(time : Array[Int], values : Array[Double]) -> Double {
  let n = time.length().min(values.length())
  if n < 2 {
    return 0.0
  }
  let x_mean = ts_time_mean(time[:n].to_owned())
  let y_mean = ts_mean(values[:n].to_owned())
  let mut numerator = 0.0
  let mut denominator = 0.0
  for i in 0..<n {
    let x = time[i].to_double() - x_mean
    numerator += x * (values[i] - y_mean)
    denominator += x * x
  }
  if denominator == 0.0 {
    0.0
  } else {
    numerator / denominator
  }
}

///|
/// Computes the intercept of a linear trend.
pub fn linear_trend_intercept(
  time : Array[Int],
  values : Array[Double],
) -> Double {
  mean_or(values, 0.0) - linear_trend_slope(time, values) * ts_time_mean(time)
}

///|
/// Fits an autoregressive model of order one by ordinary least squares.
pub fn fit_ar1(values : Array[Double]) -> Array[Double] {
  if values.length() < 2 {
    return [0.0, ts_mean(values)]
  }
  let lagged = values[:values.length() - 1].to_owned()
  let current = values[1:].to_owned()
  let slope = covariance(lagged, current) / variance(lagged)
  let intercept = mean(current) - slope * mean(lagged)
  [intercept, if is_finite(slope) { slope } else { 0.0 }]
}

///|
/// Fits an autoregressive model with several lags using a regularized normal equation.
pub fn fit_ar(
  values : Array[Double],
  order : Int,
  ridge? : Double = 1.0e-6,
) -> Array[Double] {
  let width = if order > 0 { order } else { 1 }
  if values.length() <= width {
    return Array::make(width + 1, 0.0)
  }
  let rows : Array[Array[Double]] = Array::new(capacity=values.length() - width)
  let target = Array::new(capacity=values.length() - width)
  for i in width..<values.length() {
    let row = Array::new(capacity=width + 1)
    row.push(1.0)
    for lag in 1..<=width {
      row.push(values[i - lag])
    }
    rows.push(row)
    target.push(values[i])
  }
  let transposed = matrix_transpose(rows)
  let normal = add_ridge(matrix_multiply(transposed, rows), ridge)
  let rhs = matrix_vector_multiply(transposed, target)
  solve_linear_system(normal, rhs)
}

///|
/// Forecasts h steps from an autoregressive coefficient vector.
pub fn forecast_ar(
  values : Array[Double],
  coefficients : Array[Double],
  horizon : Int,
  start_time? : Int = 0,
) -> ForecastResult {
  let steps = if horizon > 0 { horizon } else { 0 }
  let history = values.copy()
  let forecasts : Array[Double] = Array::new(capacity=steps)
  let order = if coefficients.length() > 1 {
    coefficients.length() - 1
  } else {
    1
  }
  for _ in 0..<steps {
    let mut prediction = if coefficients.length() > 0 {
      coefficients[0]
    } else {
      0.0
    }
    for lag in 1..<=order {
      let index = history.length() - lag
      if index >= 0 && lag < coefficients.length() {
        prediction += coefficients[lag] * history[index]
      }
    }
    forecasts.push(prediction)
    history.push(prediction)
  }
  let residual_scale = if values.length() > order {
    std_dev(k_difference(values, 1))
  } else {
    0.0
  }
  let time = Array::new(capacity=steps)
  let lower = Array::new(capacity=steps)
  let upper = Array::new(capacity=steps)
  for i in 0..<steps {
    time.push(start_time + i + 1)
    let radius = 1.96 * residual_scale * (i + 1).to_double().sqrt()
    lower.push(forecasts[i] - radius)
    upper.push(forecasts[i] + radius)
  }
  {
    time,
    forecast: forecasts,
    lower,
    upper,
    residual_standard_error: residual_scale,
    model_order: order,
  }
}

///|
/// Interpolates finite endpoint values across a missing run.
pub fn interpolate_missing(values : Array[Double]) -> Array[Double] {
  let result = values.copy()
  let n = result.length()
  for i in 0..<n {
    if !is_finite(result[i]) {
      let mut left = i - 1
      while left >= 0 && !is_finite(result[left]) {
        left -= 1
      }
      let mut right = i + 1
      while right < n && !is_finite(result[right]) {
        right += 1
      }
      if left >= 0 && right < n {
        result[i] = result[left] +
          (result[right] - result[left]) *
          (i - left).to_double() /
          (right - left).to_double()
      } else if left >= 0 {
        result[i] = result[left]
      } else if right < n {
        result[i] = result[right]
      } else {
        result[i] = 0.0
      }
    }
  }
  result
}

///|
/// Computes seasonal means by position within a cycle.
pub fn seasonal_indices(values : Array[Double], period : Int) -> Array[Double] {
  let width = if period > 0 { period } else { 1 }
  let result = Array::make(width, 0.0)
  let counts = Array::make(width, 0)
  for i in 0..<values.length() {
    let index = i % width
    result[index] += values[i]
    counts[index] += 1
  }
  let average = mean_or(values, 0.0)
  for i in 0..<width {
    result[i] = if counts[i] == 0 || average == 0.0 {
      1.0
    } else {
      result[i] / counts[i].to_double() / average
    }
  }
  result
}

///|
/// Removes a seasonal component by dividing by seasonal indices.
pub fn deseasonalize(
  values : Array[Double],
  indices : Array[Double],
) -> Array[Double] {
  let result = Array::new(capacity=values.length())
  if indices.length() == 0 {
    return values.copy()
  }
  for i in 0..<values.length() {
    let index = i % indices.length()
    result.push(
      if indices[index] == 0.0 {
        values[i]
      } else {
        values[i] / indices[index]
      },
    )
  }
  result
}

///|
/// Adds a seasonal component to a level series.
pub fn seasonalize(
  values : Array[Double],
  indices : Array[Double],
) -> Array[Double] {
  let result = Array::new(capacity=values.length())
  if indices.length() == 0 {
    return values.copy()
  }
  for i in 0..<values.length() {
    result.push(values[i] * indices[i % indices.length()])
  }
  result
}

///|
/// Estimates an interrupted time-series level and slope change.
pub fn interrupted_time_series(
  time : Array[Int],
  values : Array[Double],
  intervention_time : Int,
) -> InterruptedEffect {
  let pre_time = Array::new()
  let pre_values = Array::new()
  let post_time = Array::new()
  let post_values = Array::new()
  for i in 0..<time.length().min(values.length()) {
    if time[i] < intervention_time {
      pre_time.push(time[i])
      pre_values.push(values[i])
    } else {
      post_time.push(time[i])
      post_values.push(values[i])
    }
  }
  let pre_slope = linear_trend_slope(pre_time, pre_values)
  let post_slope = linear_trend_slope(post_time, post_values)
  let pre_intercept = linear_trend_intercept(pre_time, pre_values)
  let post_intercept = linear_trend_intercept(post_time, post_values)
  let level = post_intercept +
    post_slope * intervention_time.to_double() -
    (pre_intercept + pre_slope * intervention_time.to_double())
  let differences = Array::new()
  for i in 0..<post_values.length() {
    differences.push(
      post_values[i] - (pre_intercept + pre_slope * post_time[i].to_double()),
    )
  }
  let standard_error = if differences.length() < 2 {
    0.0
  } else {
    std_dev(differences) / differences.length().to_double().sqrt()
  }
  {
    level_change: level,
    slope_change: post_slope - pre_slope,
    standard_error,
    pre_slope,
    post_slope,
    pre_observations: pre_values.length(),
    post_observations: post_values.length(),
    passes: pre_values.length() >= 3 && post_values.length() >= 3,
  }
}

///|
/// Scans a series for the largest mean-shift statistic.
pub fn detect_change_point(
  values : Array[Double],
  minimum_segment? : Int = 3,
) -> ChangePoint {
  let n = values.length()
  let minimum = if minimum_segment > 0 { minimum_segment } else { 3 }
  if n < 2 * minimum {
    return {
      index: 0,
      statistic: 0.0,
      left_mean: 0.0,
      right_mean: 0.0,
      detected: false,
    }
  }
  let mut best_index = minimum
  let mut best_statistic = 0.0
  let mut best_left = 0.0
  let mut best_right = 0.0
  for index in minimum..<(n - minimum + 1) {
    let left = values[:index].to_owned()
    let right = values[index:].to_owned()
    let difference = (mean(left) - mean(right)).abs()
    let pooled = (variance(left) / index.to_double() +
    variance(right) / (n - index).to_double()).sqrt()
    let statistic = if pooled == 0.0 { difference } else { difference / pooled }
    if statistic > best_statistic {
      best_statistic = statistic
      best_index = index
      best_left = mean(left)
      best_right = mean(right)
    }
  }
  {
    index: best_index,
    statistic: best_statistic,
    left_mean: best_left,
    right_mean: best_right,
    detected: best_statistic > 2.0,
  }
}

///|
/// Computes autocovariances through a Newey-West bandwidth.
pub fn newey_west_variance(
  values : Array[Double],
  bandwidth : Int,
) -> LongRunVariance {
  let n = values.length()
  if n < 2 {
    return {
      variance: 0.0,
      standard_error: 0.0,
      bandwidth: 0,
      autocovariances: [],
    }
  }
  let width = if bandwidth > 0 { bandwidth.min(n - 1) } else { 0 }
  let center = mean(values)
  let autocovariances : Array[Double] = Array::new(capacity=width + 1)
  for lag in 0..<=width {
    let mut covariance_value = 0.0
    for i in lag..<n {
      covariance_value += (values[i] - center) * (values[i - lag] - center)
    }
    autocovariances.push(covariance_value / n.to_double())
  }
  let mut long_run = autocovariances[0]
  for lag in 1..<=width {
    let weight = 1.0 - lag.to_double() / (width + 1).to_double()
    long_run += 2.0 * weight * autocovariances[lag]
  }
  let adjusted = long_run.max(0.0)
  {
    variance: adjusted,
    standard_error: (adjusted / n.to_double()).sqrt(),
    bandwidth: width,
    autocovariances,
  }
}

///|
/// Computes a Bartlett bandwidth from sample size.
pub fn automatic_bandwidth(sample_size : Int) -> Int {
  if sample_size <= 2 {
    0
  } else {
    (4.0 * @math.pow(sample_size.to_double() / 100.0, 0.25))
    .round()
    .to_int()
    .min(sample_size - 1)
  }
}

///|
/// Computes a block bootstrap index sequence for dependent observations.
pub fn block_bootstrap_indices(
  sample_size : Int,
  block_size : Int,
  seed : UInt64,
) -> Array[Int] {
  let n = if sample_size > 0 { sample_size } else { 0 }
  let width = if block_size > 0 { block_size.min(n.max(1)) } else { 1 }
  let result = Array::new(capacity=n)
  let rng = RandomState::new(seed)
  while result.length() < n {
    let start = if n == 0 {
      0
    } else {
      (rng.uniform() * n.to_double()).to_int()
    }
    for offset in 0..<width {
      if result.length() < n {
        result.push((start + offset) % n)
      }
    }
  }
  result
}

///|
/// Computes a moving z-score for anomaly monitoring.
pub fn rolling_z_scores(values : Array[Double], window : Int) -> Array[Double] {
  let means = moving_average(values, window)
  let deviations = rolling_standard_deviation(values, window)
  let result = Array::new(capacity=values.length())
  for i in 0..<values.length() {
    result.push(
      if deviations[i] == 0.0 {
        0.0
      } else {
        (values[i] - means[i]) / deviations[i]
      },
    )
  }
  result
}

///|
/// Detects rolling anomalies at an absolute z-score threshold.
pub fn anomaly_flags(
  values : Array[Double],
  window : Int,
  threshold? : Double = 3.0,
) -> Array[Bool] {
  let scores = rolling_z_scores(values, window)
  let result = Array::new(capacity=values.length())
  for score in scores {
    result.push(!is_finite(score) || score.abs() > threshold)
  }
  result
}

///|
/// Computes a seasonal naive forecast.
pub fn seasonal_naive_forecast(
  values : Array[Double],
  period : Int,
  horizon : Int,
  start_time? : Int = 0,
) -> ForecastResult {
  let width = if period > 0 { period } else { 1 }
  let steps = if horizon > 0 { horizon } else { 0 }
  let forecast = Array::new(capacity=steps)
  for i in 0..<steps {
    let index = values.length() - width + i % width
    forecast.push(
      if index >= 0 && index < values.length() {
        values[index]
      } else {
        ts_mean(values)
      },
    )
  }
  let error = if values.length() > width {
    std_dev(k_difference(values, width))
  } else {
    0.0
  }
  let time = Array::new(capacity=steps)
  let lower = Array::new(capacity=steps)
  let upper = Array::new(capacity=steps)
  for i in 0..<steps {
    time.push(start_time + i + 1)
    lower.push(forecast[i] - 1.96 * error)
    upper.push(forecast[i] + 1.96 * error)
  }
  {
    time,
    forecast,
    lower,
    upper,
    residual_standard_error: error,
    model_order: width,
  }
}

///|
/// Computes a compact time-series summary vector.
pub fn time_series_summary(series : CausalTimeSeries) -> Array[Double] {
  let profile = time_series_profile(series, seasonal_period=series.frequency)
  [
    profile.observations.to_double(),
    profile.missing.to_double(),
    profile.mean,
    profile.standard_deviation,
    profile.first_difference_mean,
    profile.autocorrelation_lag1,
    profile.trend_slope,
    profile.seasonal_period.to_double(),
    if profile.passes {
      1.0
    } else {
      0.0
    },
  ]
}
