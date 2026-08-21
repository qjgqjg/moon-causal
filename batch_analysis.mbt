///|
/// One completed dataset analysis in a batch run.
pub struct BatchRun {
  index : Int
  result : PipelineResult
  fingerprint : UInt64
  succeeded : Bool
}

///|
/// Aggregate operating characteristics from a batch.
pub struct BatchSummary {
  runs : Int
  successes : Int
  mean_estimate : Double
  standard_deviation : Double
  minimum_estimate : Double
  maximum_estimate : Double
  success_rate : Double
  passes : Bool
}

///|
/// Runs one analysis plan over aligned datasets.
pub fn run_batch(
  datasets : Array[CausalDataset],
  plan : AnalysisPlan,
) -> Array[BatchRun] {
  let result : Array[BatchRun] = Array::new(capacity=datasets.length())
  for i in 0..<datasets.length() {
    let analysis = run_causal_pipeline(datasets[i], plan)
    result.push({
      index: i,
      result: analysis,
      fingerprint: pipeline_fingerprint(analysis),
      succeeded: analysis.passes,
    })
  }
  result
}

///|
/// Summarizes batch estimates.
pub fn summarize_batch(runs : Array[BatchRun]) -> BatchSummary {
  let estimates : Array[Double] = Array::new()
  let mut successes = 0
  for run in runs {
    estimates.push(run.result.estimate.estimate)
    if run.succeeded {
      successes += 1
    }
  }
  let mut minimum = 0.0
  let mut maximum = 0.0
  if estimates.length() > 0 {
    minimum = estimates[0]
    maximum = estimates[0]
    for value in estimates {
      if value < minimum {
        minimum = value
      }
      if value > maximum {
        maximum = value
      }
    }
  }
  {
    runs: runs.length(),
    successes,
    mean_estimate: mean_or(estimates, 0.0),
    standard_deviation: std_dev(estimates),
    minimum_estimate: minimum,
    maximum_estimate: maximum,
    success_rate: if runs.length() == 0 {
      0.0
    } else {
      successes.to_double() / runs.length().to_double()
    },
    passes: runs.length() > 0 && successes > 0,
  }
}

///|
/// Returns only successful runs while preserving original order.
pub fn batch_successful_runs(runs : Array[BatchRun]) -> Array[BatchRun] {
  let result : Array[BatchRun] = Array::new()
  for run in runs {
    if run.succeeded {
      result.push(run)
    }
  }
  result
}

///|
/// Returns indices of failed runs.
pub fn batch_failed_indices(runs : Array[BatchRun]) -> Array[Int] {
  let result : Array[Int] = Array::new()
  for run in runs {
    if !run.succeeded {
      result.push(run.index)
    }
  }
  result
}

///|
/// Extracts point estimates from completed runs.
pub fn batch_estimates(runs : Array[BatchRun]) -> Array[Double] {
  let result : Array[Double] = Array::new(capacity=runs.length())
  for run in runs {
    result.push(run.result.estimate.estimate)
  }
  result
}

///|
/// Extracts standard errors from completed runs.
pub fn batch_standard_errors(runs : Array[BatchRun]) -> Array[Double] {
  let result : Array[Double] = Array::new(capacity=runs.length())
  for run in runs {
    result.push(run.result.estimate.standard_error)
  }
  result
}

///|
/// Extracts data quality scores from completed runs.
pub fn batch_quality_scores(runs : Array[BatchRun]) -> Array[Double] {
  let result : Array[Double] = Array::new(capacity=runs.length())
  for run in runs {
    result.push(run.result.quality.score)
  }
  result
}

///|
/// Extracts effective sample sizes from completed runs.
pub fn batch_effective_sample_sizes(runs : Array[BatchRun]) -> Array[Double] {
  let result : Array[Double] = Array::new(capacity=runs.length())
  for run in runs {
    result.push(run.result.estimate.effective_sample_size)
  }
  result
}

///|
/// Extracts interval widths from completed runs.
pub fn batch_interval_widths(runs : Array[BatchRun]) -> Array[Double] {
  let result : Array[Double] = Array::new(capacity=runs.length())
  for run in runs {
    result.push((run.result.estimate.upper - run.result.estimate.lower).abs())
  }
  result
}

///|
/// Returns a selected batch quantile of point estimates.
pub fn batch_estimate_quantile(
  runs : Array[BatchRun],
  probability : Double,
) -> Double {
  quantile(batch_estimates(runs), clamp(probability, 0.0, 1.0))
}

///|
/// Returns the point-estimate range in a batch.
pub fn batch_estimate_range(runs : Array[BatchRun]) -> Array[Double] {
  let estimates = batch_estimates(runs)
  if estimates.length() == 0 {
    [0.0, 0.0]
  } else {
    [causal_minimum(estimates), causal_maximum(estimates)]
  }
}

///|
/// Computes the weighted mean of successful estimates.
pub fn batch_successful_weighted_mean(runs : Array[BatchRun]) -> Double {
  let successful = batch_successful_runs(runs)
  let estimates = batch_estimates(successful)
  let weights : Array[Double] = Array::new(capacity=successful.length())
  for run in successful {
    weights.push(1.0 / run.result.estimate.standard_error.abs().max(1.0e-12))
  }
  weighted_mean(estimates, weights)
}

///|
/// Computes the bias of the batch mean against a known truth.
pub fn batch_bias(runs : Array[BatchRun], truth : Double) -> Double {
  mean_or(batch_estimates(runs), truth) - truth
}

///|
/// Computes the RMSE of batch estimates against a known truth.
pub fn batch_rmse(runs : Array[BatchRun], truth : Double) -> Double {
  let estimates = batch_estimates(runs)
  if estimates.length() == 0 {
    0.0
  } else {
    let mut total = 0.0
    for estimate in estimates {
      let error = estimate - truth
      total += error * error
    }
    (total / estimates.length().to_double()).sqrt()
  }
}

///|
/// Computes the fraction of runs with stable point estimates.
pub fn batch_reproducibility_rate(
  runs : Array[BatchRun],
  tolerance : Double,
) -> Double {
  if runs.length() < 2 {
    return if runs.length() == 1 { 1.0 } else { 0.0 }
  }
  let estimates = batch_estimates(runs)
  let center = mean(estimates)
  let threshold = tolerance.max(0.0)
  let mut stable = 0
  for estimate in estimates {
    if (estimate - center).abs() <= threshold {
      stable += 1
    }
  }
  stable.to_double() / estimates.length().to_double()
}

///|
/// Compares each batch estimate to a reference value.
pub fn batch_compare_to_reference(
  runs : Array[BatchRun],
  reference : Double,
) -> Array[Double] {
  let result : Array[Double] = Array::new(capacity=runs.length())
  for estimate in batch_estimates(runs) {
    result.push(estimate - reference)
  }
  result
}

///|
/// Returns the first successful run with the smallest absolute error.
pub fn batch_best_run(runs : Array[BatchRun], reference : Double) -> BatchRun? {
  let mut best : BatchRun? = None
  let mut best_error = 1.0e300
  for run in runs {
    if run.succeeded {
      let error = (run.result.estimate.estimate - reference).abs()
      if error < best_error {
        best_error = error
        best = Some(run)
      }
    }
  }
  best
}

///|
/// Returns the run with the largest effective sample size.
pub fn batch_best_overlap_run(runs : Array[BatchRun]) -> BatchRun? {
  let mut best : BatchRun? = None
  let mut best_size = -1.0
  for run in runs {
    let size = run.result.estimate.effective_sample_size
    if run.succeeded && size > best_size {
      best_size = size
      best = Some(run)
    }
  }
  best
}

///|
/// Counts pipeline stages across all runs.
pub fn batch_stage_counts(runs : Array[BatchRun]) -> Array[Int] {
  let mut maximum = 0
  for run in runs {
    if run.result.stages.length() > maximum {
      maximum = run.result.stages.length()
    }
  }
  let counts = Array::make(maximum, 0)
  for run in runs {
    for i in 0..<run.result.stages.length() {
      if run.result.stages[i].status == "passed" {
        counts[i] += 1
      }
    }
  }
  counts
}

///|
/// Returns the share of runs passing each pipeline stage.
pub fn batch_stage_pass_rates(runs : Array[BatchRun]) -> Array[Double] {
  let counts = batch_stage_counts(runs)
  counts.map(fn(count) {
    if runs.length() == 0 {
      0.0
    } else {
      count.to_double() / runs.length().to_double()
    }
  })
}

///|
/// Computes pairwise differences between adjacent batch estimates.
pub fn batch_adjacent_differences(runs : Array[BatchRun]) -> Array[Double] {
  let estimates = batch_estimates(runs)
  if estimates.length() < 2 {
    return []
  }
  let result : Array[Double] = Array::new(capacity=estimates.length() - 1)
  for i in 1..<estimates.length() {
    result.push(estimates[i] - estimates[i - 1])
  }
  result
}

///|
/// Computes a running mean of batch estimates.
pub fn batch_running_mean(runs : Array[BatchRun]) -> Array[Double] {
  let result : Array[Double] = Array::new(capacity=runs.length())
  let mut total = 0.0
  for i in 0..<runs.length() {
    total += runs[i].result.estimate.estimate
    result.push(total / (i + 1).to_double())
  }
  result
}

///|
/// Computes a running standard error of the batch mean.
pub fn batch_running_standard_error(runs : Array[BatchRun]) -> Array[Double] {
  let result : Array[Double] = Array::new(capacity=runs.length())
  for i in 0..<runs.length() {
    let prefix : Array[Double] = Array::new(capacity=i + 1)
    for j in 0..<(i + 1) {
      prefix.push(runs[j].result.estimate.estimate)
    }
    result.push(inference_monte_carlo_se(prefix))
  }
  result
}

///|
/// Returns a stable fingerprint for all batch runs.
pub fn batch_fingerprint(runs : Array[BatchRun]) -> UInt64 {
  let rows : Array[Array[Double]] = Array::new(capacity=runs.length())
  for run in runs {
    rows.push([
      run.index.to_double(),
      run.result.estimate.estimate,
      run.result.estimate.standard_error,
      run.result.estimate.effective_sample_size,
      if run.succeeded {
        1.0
      } else {
        0.0
      },
    ])
  }
  matrix_checksum(rows)
}

///|
/// Returns a vector used by benchmark result tables.
pub fn batch_summary_vector(runs : Array[BatchRun]) -> Array[Double] {
  let summary = summarize_batch(runs)
  [
    summary.runs.to_double(),
    summary.successes.to_double(),
    summary.mean_estimate,
    summary.standard_deviation,
    summary.minimum_estimate,
    summary.maximum_estimate,
    summary.success_rate,
    if summary.passes {
      1.0
    } else {
      0.0
    },
  ]
}

///|
/// Serializes a batch summary with failure indices.
pub fn batch_summary_text(runs : Array[BatchRun]) -> String {
  let summary = summarize_batch(runs)
  let builder = StringBuilder::new()
  builder.write_string("runs=")
  builder.write_string(summary.runs.to_string())
  builder.write_string(";successes=")
  builder.write_string(summary.successes.to_string())
  builder.write_string(";mean=")
  builder.write_string(summary.mean_estimate.to_string())
  builder.write_string(";sd=")
  builder.write_string(summary.standard_deviation.to_string())
  builder.write_string(";failures=")
  builder.write_string(batch_failed_indices(runs).length().to_string())
  builder.to_string()
}

///|
/// Checks that a batch meets minimum success and precision requirements.
pub fn batch_release_check(
  runs : Array[BatchRun],
  minimum_success_rate : Double,
  maximum_standard_error : Double,
) -> Bool {
  let summary = summarize_batch(runs)
  let errors = batch_standard_errors(runs)
  let mut precise = true
  for error in errors {
    if !is_finite(error) || error.abs() > maximum_standard_error {
      precise = false
    }
  }
  summary.success_rate >= clamp(minimum_success_rate, 0.0, 1.0) &&
  precise &&
  summary.passes
}

///|
/// Returns the minimum quality score observed in a batch.
pub fn batch_minimum_quality(runs : Array[BatchRun]) -> Double {
  let values = batch_quality_scores(runs)
  causal_minimum(values, fallback=0.0)
}

///|
/// Returns the minimum effective sample size observed in a batch.
pub fn batch_minimum_effective_sample_size(runs : Array[BatchRun]) -> Double {
  let values = batch_effective_sample_sizes(runs)
  causal_minimum(values, fallback=0.0)
}

///|
/// Returns whether all successful runs are finite and non-degenerate.
pub fn batch_all_successful_finite(runs : Array[BatchRun]) -> Bool {
  let successful = batch_successful_runs(runs)
  if successful.length() == 0 {
    return false
  }
  for run in successful {
    if !is_finite(run.result.estimate.estimate) ||
      !is_finite(run.result.estimate.standard_error) ||
      run.result.estimate.effective_sample_size <= 0.0 {
      return false
    }
  }
  true
}

///|
/// Returns the observed success rate without exposing batch internals.
pub fn batch_success_rate(runs : Array[BatchRun]) -> Double {
  summarize_batch(runs).success_rate
}

///|
/// Returns whether a batch has at least one usable successful result.
pub fn batch_has_success(runs : Array[BatchRun]) -> Bool {
  summarize_batch(runs).successes > 0
}

///|
/// Returns the number of successful batch runs.
pub fn batch_success_count(runs : Array[BatchRun]) -> Int {
  summarize_batch(runs).successes
}

///|
/// Returns the number of completed batch runs.
pub fn batch_run_count(runs : Array[BatchRun]) -> Int {
  summarize_batch(runs).runs
}

///|
/// Returns the number of failed batch runs.
pub fn batch_failure_count(runs : Array[BatchRun]) -> Int {
  batch_failed_indices(runs).length()
}

///|
/// Returns whether any batch run failed.
pub fn batch_has_failures(runs : Array[BatchRun]) -> Bool {
  batch_failure_count(runs) > 0
}
