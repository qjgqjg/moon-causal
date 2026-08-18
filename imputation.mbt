///|
pub struct ImputationSummary {
  imputed_values : Int
  columns : Int
  rows : Int
  column_means : Array[Double]
}

///|
pub fn impute_non_finite(
  covariates : Array[Array[Double]],
) -> (Array[Array[Double]], ImputationSummary) {
  if covariates.length() == 0 {
    return ([], { imputed_values: 0, columns: 0, rows: 0, column_means: [] })
  }
  let columns = covariates[0].length()
  let means = Array::make(columns, 0.0)
  let counts = Array::make(columns, 0)
  for row in covariates {
    for j in 0..<columns {
      if is_finite(row[j]) {
        means[j] += row[j]
        counts[j] += 1
      }
    }
  }
  for j in 0..<columns {
    if counts[j] > 0 {
      means[j] /= counts[j].to_double()
    }
  }
  let result : Array[Array[Double]] = Array::new(capacity=covariates.length())
  let mut imputed = 0
  for row in covariates {
    let copy = row.copy()
    for j in 0..<columns {
      if !is_finite(copy[j]) {
        copy[j] = means[j]
        imputed += 1
      }
    }
    result.push(copy)
  }
  (
    result,
    {
      imputed_values: imputed,
      columns,
      rows: covariates.length(),
      column_means: means,
    },
  )
}

///|
pub fn complete_case_indices(covariates : Array[Array[Double]]) -> Array[Int] {
  let result = Array::new()
  for i in 0..<covariates.length() {
    if finite_array(covariates[i]) {
      result.push(i)
    }
  }
  result
}

///|
pub fn complete_case_dataset(dataset : CausalDataset) -> CausalDataset {
  let indices = complete_case_indices(dataset.covariates)
  dataset.select(indices)
}

///|
pub fn missingness_indicator(
  covariates : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=covariates.length())
  for row in covariates {
    let indicators = Array::new(capacity=row.length())
    for value in row {
      indicators.push(if is_finite(value) { 0.0 } else { 1.0 })
    }
    result.push(indicators)
  }
  result
}

///|
pub fn missingness_rate(covariates : Array[Array[Double]]) -> Double {
  let mut missing = 0
  let mut total = 0
  for row in covariates {
    for value in row {
      total += 1
      if !is_finite(value) {
        missing += 1
      }
    }
  }
  if total == 0 {
    0.0
  } else {
    missing.to_double() / total.to_double()
  }
}

///|
pub fn mean_impute_vector(values : Array[Double]) -> Array[Double] {
  let finite = Array::new()
  for value in values {
    if is_finite(value) {
      finite.push(value)
    }
  }
  let center = mean(finite)
  let result = values.copy()
  for i in 0..<result.length() {
    if !is_finite(result[i]) {
      result[i] = center
    }
  }
  result
}
