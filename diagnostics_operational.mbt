///|
/// A named operational metric suitable for a release dashboard.
pub struct OperationalMetric {
  name : String
  value : Double
  target : Double
  tolerance : Double
  status : String
  weight : Double
}

///|
/// A point-in-time collection of production diagnostics.
pub struct OperationalSnapshot {
  run_id : String
  metrics : Array[OperationalMetric]
  score : Double
  passes : Bool
  fingerprint : UInt64
}

///|
/// Comparison of two operational snapshots.
pub struct OperationalComparison {
  changed : Int
  improved : Int
  degraded : Int
  largest_change : Double
  score_delta : Double
  passes : Bool
}

///|
/// Creates a metric with a two-sided tolerance band.
pub fn operational_metric(
  name : String,
  value : Double,
  target : Double,
  tolerance? : Double = 0.0,
  weight? : Double = 1.0,
) -> OperationalMetric {
  let safe_tolerance = tolerance.max(0.0)
  let safe_weight = if is_finite(weight) && weight > 0.0 { weight } else { 1.0 }
  let finite = is_finite(value) && is_finite(target)
  let status = if !finite {
    "invalid"
  } else if (value - target).abs() <= safe_tolerance {
    "passed"
  } else if value >= target {
    "above-target"
  } else {
    "below-target"
  }
  {
    name,
    value,
    target,
    tolerance: safe_tolerance,
    status,
    weight: safe_weight,
  }
}

///|
/// Returns whether a metric is within its release band.
pub fn operational_metric_passes(metric : OperationalMetric) -> Bool {
  is_finite(metric.value) &&
  is_finite(metric.target) &&
  (metric.value - metric.target).abs() <= metric.tolerance
}

///|
/// Returns a normalized metric score in [0, 1].
pub fn operational_metric_score(metric : OperationalMetric) -> Double {
  if !is_finite(metric.value) || !is_finite(metric.target) {
    0.0
  } else if metric.tolerance == 0.0 {
    if metric.value == metric.target {
      1.0
    } else {
      0.0
    }
  } else {
    clamp(
      1.0 - (metric.value - metric.target).abs() / metric.tolerance,
      0.0,
      1.0,
    )
  }
}

///|
/// Builds a snapshot and computes its weighted health score.
pub fn operational_snapshot(
  run_id : String,
  metrics : Array[OperationalMetric],
) -> OperationalSnapshot {
  let mut weighted = 0.0
  let mut total_weight = 0.0
  let mut passes = metrics.length() > 0
  let rows : Array[Array[Double]] = Array::new(capacity=metrics.length())
  for metric in metrics {
    let score = operational_metric_score(metric)
    weighted += score * metric.weight
    total_weight += metric.weight
    if !operational_metric_passes(metric) {
      passes = false
    }
    rows.push([metric.value, metric.target, metric.tolerance, metric.weight])
  }
  let score = if total_weight == 0.0 { 0.0 } else { weighted / total_weight }
  {
    run_id,
    metrics: metrics.copy(),
    score,
    passes,
    fingerprint: matrix_checksum(rows),
  }
}

///|
/// Appends a metric to an existing snapshot.
pub fn operational_snapshot_with(
  snapshot : OperationalSnapshot,
  metric : OperationalMetric,
) -> OperationalSnapshot {
  let metrics = snapshot.metrics.copy()
  metrics.push(metric)
  operational_snapshot(snapshot.run_id, metrics)
}

///|
/// Returns metric names in snapshot order.
pub fn operational_metric_names(
  snapshot : OperationalSnapshot,
) -> Array[String] {
  let result : Array[String] = Array::new(capacity=snapshot.metrics.length())
  for metric in snapshot.metrics {
    result.push(metric.name)
  }
  result
}

///|
/// Returns metric values in snapshot order.
pub fn operational_metric_values(
  snapshot : OperationalSnapshot,
) -> Array[Double] {
  let result : Array[Double] = Array::new(capacity=snapshot.metrics.length())
  for metric in snapshot.metrics {
    result.push(metric.value)
  }
  result
}

///|
/// Returns metrics that do not meet their target bands.
pub fn operational_failures(
  snapshot : OperationalSnapshot,
) -> Array[OperationalMetric] {
  let result : Array[OperationalMetric] = Array::new()
  for metric in snapshot.metrics {
    if !operational_metric_passes(metric) {
      result.push(metric)
    }
  }
  result
}

///|
/// Returns a seven-element snapshot summary.
pub fn operational_snapshot_summary(
  snapshot : OperationalSnapshot,
) -> Array[Double] {
  [
    snapshot.metrics.length().to_double(),
    operational_failures(snapshot).length().to_double(),
    snapshot.score,
    if snapshot.passes {
      1.0
    } else {
      0.0
    },
    snapshot.fingerprint.to_double(),
    mean_or(operational_metric_values(snapshot), 0.0),
    std_dev(operational_metric_values(snapshot)),
  ]
}

///|
/// Compares matching metrics from two snapshots.
pub fn compare_operational_snapshots(
  baseline : OperationalSnapshot,
  current : OperationalSnapshot,
  tolerance? : Double = 1.0e-12,
) -> OperationalComparison {
  let mut changed = 0
  let mut improved = 0
  let mut degraded = 0
  let mut largest = 0.0
  for current_metric in current.metrics {
    for baseline_metric in baseline.metrics {
      if current_metric.name == baseline_metric.name {
        let delta = current_metric.value - baseline_metric.value
        let magnitude = delta.abs()
        if magnitude > tolerance {
          changed += 1
          if operational_metric_score(current_metric) >
            operational_metric_score(baseline_metric) {
            improved += 1
          } else {
            degraded += 1
          }
        }
        if magnitude > largest {
          largest = magnitude
        }
        break
      }
    }
  }
  {
    changed,
    improved,
    degraded,
    largest_change: largest,
    score_delta: current.score - baseline.score,
    passes: degraded == 0 && current.passes,
  }
}

///|
/// Computes the average absolute standardized difference.
pub fn operational_mean_smd(metrics : Array[BalanceMetric]) -> Double {
  if metrics.length() == 0 {
    0.0
  } else {
    let mut total = 0.0
    for metric in metrics {
      total += metric.standardized_difference.abs()
    }
    total / metrics.length().to_double()
  }
}

///|
/// Returns the maximum absolute standardized difference.
pub fn operational_max_smd(metrics : Array[BalanceMetric]) -> Double {
  let mut result = 0.0
  for metric in metrics {
    if metric.standardized_difference.abs() > result {
      result = metric.standardized_difference.abs()
    }
  }
  result
}

///|
/// Builds balance diagnostics as operational metrics.
pub fn operational_balance_metrics(
  metrics : Array[BalanceMetric],
  maximum_smd? : Double = 0.1,
) -> Array[OperationalMetric] {
  let result : Array[OperationalMetric] = Array::new()
  let threshold = maximum_smd.max(1.0e-12)
  for metric in metrics {
    let value = metric.standardized_difference.abs()
    result.push(
      operational_metric(
        "smd:{metric.name}",
        value,
        0.0,
        tolerance=threshold,
        weight=1.0,
      ),
    )
  }
  result
}

///|
/// Checks whether all balance metrics meet an SMD threshold.
pub fn operational_balance_passes(
  metrics : Array[BalanceMetric],
  maximum_smd? : Double = 0.1,
) -> Bool {
  let threshold = maximum_smd.max(0.0)
  for metric in metrics {
    if metric.standardized_difference.abs() > threshold {
      return false
    }
  }
  metrics.length() > 0
}

///|
/// Returns a metric for an estimate's absolute standard error.
pub fn operational_estimate_precision(
  estimate : Estimate,
  maximum_standard_error : Double,
) -> OperationalMetric {
  operational_metric(
    "standard-error",
    estimate.standard_error.abs(),
    0.0,
    tolerance=maximum_standard_error.max(1.0e-12),
    weight=2.0,
  )
}

///|
/// Returns the width of an estimate's confidence interval.
pub fn operational_interval_width(estimate : Estimate) -> Double {
  (estimate.upper - estimate.lower).abs()
}

///|
/// Creates interval-width and effective-sample-size metrics.
pub fn operational_estimate_metrics(
  estimate : Estimate,
  maximum_interval_width : Double,
  minimum_effective_sample_size : Double,
) -> Array[OperationalMetric] {
  [
    operational_metric(
      "interval-width",
      operational_interval_width(estimate),
      0.0,
      tolerance=maximum_interval_width.max(1.0e-12),
      weight=1.0,
    ),
    operational_metric(
      "effective-sample-size",
      estimate.effective_sample_size,
      minimum_effective_sample_size.max(1.0),
      tolerance=0.0,
      weight=2.0,
    ),
  ]
}

///|
/// Creates metrics from a completed high-level pipeline.
pub fn operational_pipeline_snapshot(
  run_id : String,
  result : PipelineResult,
  minimum_quality_score? : Double = 0.8,
  minimum_effective_sample_size? : Double = 10.0,
) -> OperationalSnapshot {
  let metrics : Array[OperationalMetric] = Array::new()
  metrics.push(
    operational_metric(
      "quality-score",
      result.quality.score,
      clamp(minimum_quality_score, 0.0, 1.0),
      tolerance=0.0,
      weight=2.0,
    ),
  )
  metrics.push(
    operational_metric(
      "pipeline-score",
      result.score,
      clamp(minimum_quality_score, 0.0, 1.0),
      tolerance=0.0,
      weight=2.0,
    ),
  )
  metrics.push(
    operational_metric(
      "effective-sample-size",
      result.estimate.effective_sample_size,
      minimum_effective_sample_size.max(1.0),
      tolerance=0.0,
      weight=2.0,
    ),
  )
  metrics.push(
    operational_metric(
      "overlap-minimum",
      result.positivity.minimum_score,
      0.05,
      tolerance=0.0,
      weight=1.0,
    ),
  )
  operational_snapshot(run_id, metrics)
}

///|
/// Computes moving averages for operational time series.
pub fn operational_moving_average(
  values : Array[Double],
  window : Int,
) -> Array[Double] {
  let width = window.max(1)
  let result : Array[Double] = Array::new(capacity=values.length())
  for i in 0..<values.length() {
    let start = (i + 1 - width).max(0)
    let mut total = 0.0
    let mut count = 0
    for j in start..<(i + 1) {
      if is_finite(values[j]) {
        total += values[j]
        count += 1
      }
    }
    result.push(if count == 0 { 0.0 } else { total / count.to_double() })
  }
  result
}

///|
/// Computes first differences for a monitored series.
pub fn operational_differences(values : Array[Double]) -> Array[Double] {
  if values.length() < 2 {
    return []
  }
  let result : Array[Double] = Array::new(capacity=values.length() - 1)
  for i in 1..<values.length() {
    result.push(values[i] - values[i - 1])
  }
  result
}

///|
/// Computes a robust median absolute deviation.
pub fn operational_mad(values : Array[Double]) -> Double {
  let finite = causal_finite_values(values)
  if finite.length() == 0 {
    return 0.0
  }
  let center = quantile(finite, 0.5)
  let deviations = finite.map(fn(value) { (value - center).abs() })
  quantile(deviations, 0.5)
}

///|
/// Returns a robust z-score for each finite value.
pub fn operational_robust_zscores(values : Array[Double]) -> Array[Double] {
  let finite = causal_finite_values(values)
  let center = quantile(finite, 0.5)
  let deviation = operational_mad(finite).max(1.0e-12)
  values.map(fn(value) {
    if is_finite(value) {
      (value - center) / (1.4826 * deviation)
    } else {
      0.0
    }
  })
}

///|
/// Flags robust outliers beyond an absolute z-score threshold.
pub fn operational_outlier_flags(
  values : Array[Double],
  threshold? : Double = 3.5,
) -> Array[Bool] {
  let limit = threshold.max(0.0)
  operational_robust_zscores(values).map(fn(value) { value.abs() > limit })
}

///|
/// Computes drift metrics between two finite vectors.
pub fn operational_vector_drift(
  baseline : Array[Double],
  current : Array[Double],
) -> Array[OperationalMetric] {
  let base = causal_finite_values(baseline)
  let now = causal_finite_values(current)
  let base_mean = mean_or(base, 0.0)
  let now_mean = mean_or(now, 0.0)
  let base_scale = std_dev(base).max(1.0e-12)
  let mean_delta = (now_mean - base_mean).abs()
  let scale_ratio = std_dev(now) / base_scale
  [
    operational_metric(
      "mean-drift",
      mean_delta,
      0.0,
      tolerance=base_scale * 0.1,
    ),
    operational_metric("scale-ratio", scale_ratio, 1.0, tolerance=0.2),
    operational_metric(
      "sample-size",
      now.length().to_double(),
      base.length().to_double(),
      tolerance=base.length().to_double() * 0.2,
    ),
  ]
}

///|
/// Computes a trend slope using centered least squares.
pub fn operational_trend_slope(values : Array[Double]) -> Double {
  let n = values.length()
  if n < 2 {
    0.0
  } else {
    let x = Array::new(capacity=n)
    for i in 0..<n {
      x.push(i.to_double())
    }
    covariance(x, values) / variance(x).max(1.0e-12)
  }
}

///|
/// Returns a metric describing a bounded trend.
pub fn operational_trend_metric(
  values : Array[Double],
  maximum_absolute_slope : Double,
) -> OperationalMetric {
  operational_metric(
    "trend-slope",
    operational_trend_slope(values).abs(),
    0.0,
    tolerance=maximum_absolute_slope.max(1.0e-12),
  )
}

///|
/// Returns a compact trend summary.
pub fn operational_trend_summary(values : Array[Double]) -> Array[Double] {
  let differences = operational_differences(values)
  [
    values.length().to_double(),
    mean_or(values, 0.0),
    std_dev(values),
    operational_trend_slope(values),
    mean_or(differences, 0.0),
    operational_mad(values),
  ]
}

///|
/// Serializes an operational snapshot for a text artifact.
pub fn operational_snapshot_text(snapshot : OperationalSnapshot) -> String {
  let builder = StringBuilder::new()
  builder.write_string("run_id=")
  builder.write_string(snapshot.run_id)
  builder.write_string("\nscore=")
  builder.write_string(snapshot.score.to_string())
  builder.write_string("\npasses=")
  builder.write_string(snapshot.passes.to_string())
  builder.write_string("\n")
  for metric in snapshot.metrics {
    builder.write_string(metric.name)
    builder.write_string("=")
    builder.write_string(metric.value.to_string())
    builder.write_string(";target=")
    builder.write_string(metric.target.to_string())
    builder.write_string(";status=")
    builder.write_string(metric.status)
    builder.write_string("\n")
  }
  builder.write_string("fingerprint=")
  builder.write_string(snapshot.fingerprint.to_string())
  builder.to_string()
}

///|
/// Converts a monitoring snapshot into a report section.
pub fn operational_report_section(
  snapshot : OperationalSnapshot,
) -> ReportSection {
  let lines : Array[String] = Array::new()
  lines.push("run_id=\{snapshot.run_id}")
  lines.push("score=\{snapshot.score}")
  lines.push("passes=\{snapshot.passes}")
  for metric in snapshot.metrics {
    lines.push("\{metric.name}=\{metric.value};status=\{metric.status}")
  }
  {
    title: "operational-diagnostics",
    status: if snapshot.passes {
      "passed"
    } else {
      "warning"
    },
    lines,
    fingerprint: snapshot.fingerprint,
  }
}

///|
/// Computes a dashboard from quality, overlap, and effect diagnostics.
pub fn operational_causal_snapshot(
  run_id : String,
  quality : DatasetQuality,
  positivity : PositivityProfile,
  effect : AdvancedEffect,
) -> OperationalSnapshot {
  let metrics : Array[OperationalMetric] = Array::new()
  metrics.push(
    operational_metric("quality", quality.score, 0.8, tolerance=0.0, weight=2.0),
  )
  metrics.push(
    operational_metric(
      "overlap",
      positivity.effective_sample_size,
      10.0,
      tolerance=0.0,
      weight=2.0,
    ),
  )
  metrics.push(
    operational_metric(
      "precision",
      effect.standard_error.abs(),
      0.0,
      tolerance=1.0,
      weight=1.0,
    ),
  )
  metrics.push(
    operational_metric(
      "estimate-finite",
      if is_finite(effect.estimate) {
        1.0
      } else {
        0.0
      },
      1.0,
      tolerance=0.0,
      weight=2.0,
    ),
  )
  operational_snapshot(run_id, metrics)
}

///|
/// Calculates the fraction of valid metrics in a collection.
pub fn operational_valid_fraction(metrics : Array[OperationalMetric]) -> Double {
  if metrics.length() == 0 {
    0.0
  } else {
    let valid = metrics.fold(init=0, fn(total, metric) {
      if is_finite(metric.value) {
        total + 1
      } else {
        total
      }
    })
    valid.to_double() / metrics.length().to_double()
  }
}

///|
/// Returns whether a snapshot has no invalid metric values.
pub fn operational_snapshot_is_finite(snapshot : OperationalSnapshot) -> Bool {
  operational_valid_fraction(snapshot.metrics) == 1.0 &&
  is_finite(snapshot.score)
}
