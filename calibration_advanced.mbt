///|
/// Reliability bin for probability calibration.
pub struct ReliabilityBin {
  lower : Double
  upper : Double
  count : Int
  mean_prediction : Double
  observed_rate : Double
  absolute_error : Double
}

///|
/// Calibration model summary.
pub struct CalibrationModel {
  intercept : Double
  slope : Double
  log_loss : Double
  brier : Double
  ece : Double
  mce : Double
  passes : Bool
}

///|
/// Creates equal-width reliability bins.
pub fn reliability_bins(
  probabilities : Array[Double],
  outcomes : Array[Bool],
  bins : Int,
) -> Array[ReliabilityBin] {
  let width = if bins > 0 { bins } else { 1 }
  let counts = Array::make(width, 0)
  let predictions = Array::make(width, 0.0)
  let observations = Array::make(width, 0.0)
  let n = probabilities.length().min(outcomes.length())
  for i in 0..<n {
    let p = clamp(probabilities[i], 0.0, 0.999999)
    let index = (p * width.to_double()).to_int()
    counts[index] += 1
    predictions[index] += p
    if outcomes[i] {
      observations[index] += 1.0
    }
  }
  let result : Array[ReliabilityBin] = Array::new(capacity=width)
  for i in 0..<width {
    let count = counts[i].to_double()
    let predicted = if count == 0.0 { 0.0 } else { predictions[i] / count }
    let observed = if count == 0.0 { 0.0 } else { observations[i] / count }
    result.push({
      lower: i.to_double() / width.to_double(),
      upper: (i + 1).to_double() / width.to_double(),
      count: counts[i],
      mean_prediction: predicted,
      observed_rate: observed,
      absolute_error: (predicted - observed).abs(),
    })
  }
  result
}

///|
/// Computes expected calibration error.
pub fn calibration_ece(bins : Array[ReliabilityBin]) -> Double {
  let total = bins
    .fold(init=0, fn(sum_count, bin) { sum_count + bin.count })
    .to_double()
    .max(1.0)
  let mut result = 0.0
  for bin in bins {
    result += bin.count.to_double() / total * bin.absolute_error
  }
  result
}

///|
/// Computes maximum calibration error over non-empty bins.
pub fn calibration_mce(bins : Array[ReliabilityBin]) -> Double {
  let mut result = 0.0
  for bin in bins {
    if bin.count > 0 && bin.absolute_error > result {
      result = bin.absolute_error
    }
  }
  result
}

///|
/// Fits a calibration intercept and slope by weighted least squares on logits.
pub fn calibration_logit_fit(
  probabilities : Array[Double],
  outcomes : Array[Bool],
) -> CalibrationModel {
  let n = probabilities.length().min(outcomes.length())
  let logits = Array::new(capacity=n)
  let targets = Array::new(capacity=n)
  for i in 0..<n {
    let p = safe_probability(probabilities[i])
    logits.push(safe_logit(p))
    targets.push(if outcomes[i] { 1.0 } else { 0.0 })
  }
  let slope = if variance(logits) == 0.0 {
    0.0
  } else {
    covariance(logits, targets) / variance(logits)
  }
  let intercept = mean_or(targets, 0.0) - slope * mean_or(logits, 0.0)
  let brier = brier_score(probabilities[:n].to_owned(), outcomes[:n].to_owned())
  let loss = log_loss(probabilities[:n].to_owned(), outcomes[:n].to_owned())
  let bins = reliability_bins(probabilities, outcomes, 10)
  let ece = calibration_ece(bins)
  let mce = calibration_mce(bins)
  {
    intercept,
    slope,
    log_loss: loss,
    brier,
    ece,
    mce,
    passes: ece <= 0.1 && is_finite(loss),
  }
}

///|
/// Applies an intercept and slope calibration model to scores.
pub fn apply_logit_calibration(
  probabilities : Array[Double],
  intercept : Double,
  slope : Double,
) -> Array[Double] {
  let result = Array::new(capacity=probabilities.length())
  for probability in probabilities {
    result.push(
      probability_from_odds(
        @math.exp(intercept + slope * safe_logit(probability)),
      ),
    )
  }
  result
}

///|
/// Performs pool-adjacent-violators isotonic calibration.
pub fn isotonic_calibration(
  probabilities : Array[Double],
  outcomes : Array[Bool],
) -> Array[Double] {
  let n = probabilities.length().min(outcomes.length())
  let order : Array[Int] = Array::new(capacity=n)
  for i in 0..<n {
    order.push(i)
  }
  for i in 1..<n {
    let value = order[i]
    let mut cursor = i
    while cursor > 0 && probabilities[order[cursor - 1]] > probabilities[value] {
      order[cursor] = order[cursor - 1]
      cursor -= 1
    }
    order[cursor] = value
  }
  let fitted = Array::make(n, 0.0)
  let mut block_start = 0
  while block_start < n {
    let mut block_end = block_start + 1
    let mut sum_value = if outcomes[order[block_start]] { 1.0 } else { 0.0 }
    while block_end < n {
      let current_mean = sum_value / (block_end - block_start).to_double()
      let next_value = if outcomes[order[block_end]] { 1.0 } else { 0.0 }
      if next_value >= current_mean {
        sum_value += next_value
        block_end += 1
      } else {
        break
      }
    }
    let mean_value = sum_value / (block_end - block_start).to_double()
    for position in block_start..<block_end {
      fitted[order[position]] = mean_value
    }
    block_start = block_end
  }
  fitted
}

///|
/// Computes calibration drift between two score samples.
pub fn calibration_drift(
  reference : Array[Double],
  current : Array[Double],
) -> Array[Double] {
  let reference_mean = mean_or(reference, 0.0)
  let current_mean = mean_or(current, 0.0)
  let reference_scale = std_dev(reference)
  let current_scale = std_dev(current)
  [
    current_mean - reference_mean,
    current_scale - reference_scale,
    population_stability_index(reference, current),
    kolmogorov_smirnov_distance(reference, current),
  ]
}

///|
/// Returns a compact calibration summary vector.
pub fn advanced_calibration_summary(model : CalibrationModel) -> Array[Double] {
  [
    model.intercept,
    model.slope,
    model.log_loss,
    model.brier,
    model.ece,
    model.mce,
    if model.passes {
      1.0
    } else {
      0.0
    },
  ]
}
