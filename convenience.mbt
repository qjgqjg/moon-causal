///|
pub fn feature_means(covariates : Array[Array[Double]]) -> Array[Double] {
  summarize_covariates(covariates).means
}

///|
pub fn feature_standard_deviations(
  covariates : Array[Array[Double]],
) -> Array[Double] {
  summarize_covariates(covariates).standard_deviations
}

///|
pub fn feature_minima(covariates : Array[Array[Double]]) -> Array[Double] {
  summarize_covariates(covariates).minima
}

///|
pub fn feature_maxima(covariates : Array[Array[Double]]) -> Array[Double] {
  summarize_covariates(covariates).maxima
}

///|
pub fn pairwise_difference(values : Array[Double]) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=values.length())
  for first in values {
    let row = Array::new(capacity=values.length())
    for second in values {
      row.push(first - second)
    }
    result.push(row)
  }
  result
}

///|
pub fn absolute_pairwise_difference(
  values : Array[Double],
) -> Array[Array[Double]] {
  let differences = pairwise_difference(values)
  for row in differences {
    for i in 0..<row.length() {
      row[i] = row[i].abs()
    }
  }
  differences
}

///|
pub fn binned_counts(values : Array[Double], bins : Int) -> Array[Int] {
  let count = if bins < 1 { 1 } else { bins }
  let result = Array::make(count, 0)
  let minimum = quantile(values, 0.0)
  let maximum = quantile(values, 1.0)
  for value in values {
    let mut index = if maximum == minimum {
      0
    } else {
      ((value - minimum) / (maximum - minimum) * count.to_double()).to_int()
    }
    if index < 0 {
      index = 0
    }
    if index >= count {
      index = count - 1
    }
    result[index] += 1
  }
  result
}

///|
pub fn treatment_outcome_covariance(
  outcomes : Array[Double],
  treatment : Array[Bool],
) -> Double {
  let treatment_numeric = Array::new(capacity=treatment.length())
  for value in treatment {
    treatment_numeric.push(if value { 1.0 } else { 0.0 })
  }
  covariance(outcomes, treatment_numeric)
}

///|
pub fn crude_effect_size(
  outcomes : Array[Double],
  treatment : Array[Bool],
) -> Double {
  standardized_effect_size(
    estimate_difference_in_means(outcomes, treatment).estimate,
    outcomes,
  )
}

///|
pub fn overlap_fraction(
  propensity_scores : Array[Double],
  treatment : Array[Bool],
) -> Double {
  overlap_report(propensity_scores, treatment).overlap_fraction
}

///|
pub fn common_support_width(
  propensity_scores : Array[Double],
  treatment : Array[Bool],
) -> Double {
  let report = overlap_report(propensity_scores, treatment)
  report.common_support_upper - report.common_support_lower
}

///|
pub fn estimate_quality_flags(estimate : Estimate) -> Array[Bool] {
  [
    estimate.sample_size > 0,
    estimate.effective_sample_size > 1.0,
    is_finite(estimate.estimate),
    is_finite(estimate.standard_error),
    estimate.lower <= estimate.upper,
  ]
}

///|
pub fn model_quality_flags(model : ModelFit) -> Array[Bool] {
  [
    model.converged,
    model.iterations > 0,
    is_finite(model.loss),
    finite_array(model.coefficients),
  ]
}
