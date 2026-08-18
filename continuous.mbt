///|
pub struct DoseResponsePoint {
  dose : Double
  count : Int
  mean_outcome : Double
  standard_error : Double
}

///|
pub fn dose_response_curve(
  dose : Array[Double],
  outcomes : Array[Double],
  bins : Int,
) -> Array[DoseResponsePoint] {
  let count_bins = if bins < 1 { 1 } else { bins }
  let result : Array[DoseResponsePoint] = Array::new(capacity=count_bins)
  let minimum = quantile(dose, 0.0)
  let maximum = quantile(dose, 1.0)
  let width = if maximum == minimum {
    1.0
  } else {
    (maximum - minimum) / count_bins.to_double()
  }
  for bin in 0..<count_bins {
    let values = Array::new()
    let doses = Array::new()
    let lower = minimum + bin.to_double() * width
    let upper = if bin == count_bins - 1 {
      maximum + 1.0e-12
    } else {
      lower + width
    }
    for i in 0..<dose.length() {
      if dose[i] >= lower && dose[i] < upper {
        values.push(outcomes[i])
        doses.push(dose[i])
      }
    }
    result.push({
      dose: mean(doses),
      count: values.length(),
      mean_outcome: mean(values),
      standard_error: if values.length() <= 1 {
        0.0
      } else {
        std_dev(values) / values.length().to_double().sqrt()
      },
    })
  }
  result
}

///|
pub fn generalized_propensity_scores(
  dose : Array[Double],
  covariates : Array[Array[Double]],
) -> Array[Double] {
  let model = fit_linear_outcome_model(covariates, dose, ridge=1.0e-6)
  let predictions = predict_outcomes(model, covariates)
  let residuals = Array::new(capacity=dose.length())
  for i in 0..<dose.length() {
    residuals.push(dose[i] - predictions[i])
  }
  let scale = std_dev(residuals, sample=false)
  let result = Array::new(capacity=dose.length())
  for residual in residuals {
    result.push(
      if scale == 0.0 {
        0.5
      } else {
        @math.exp(-0.5 * residual * residual / (scale * scale))
      },
    )
  }
  result
}

///|
pub fn overlap_weights(
  treatment : Array[Bool],
  propensity_scores : Array[Double],
) -> Array[Double] {
  let n = if treatment.length() < propensity_scores.length() {
    treatment.length()
  } else {
    propensity_scores.length()
  }
  let result = Array::new(capacity=n)
  for i in 0..<n {
    let p = safe_probability(propensity_scores[i])
    result.push(if treatment[i] { 1.0 - p } else { p })
  }
  result
}

///|
pub fn dose_weighted_mean(
  dose : Array[Double],
  outcomes : Array[Double],
  weights : Array[Double],
) -> Double {
  weighted_mean(outcomes, weights) + weighted_mean(dose, weights) * 0.0
}
