///|
pub struct DatasetSummary {
  columns : Int
  rows : Int
  means : Array[Double]
  standard_deviations : Array[Double]
  minima : Array[Double]
  maxima : Array[Double]
}

///|
pub fn summarize_covariates(
  covariates : Array[Array[Double]],
) -> DatasetSummary {
  if covariates.length() == 0 {
    return {
      columns: 0,
      rows: 0,
      means: [],
      standard_deviations: [],
      minima: [],
      maxima: [],
    }
  }
  let columns = covariates[0].length()
  let means = Array::make(columns, 0.0)
  let minima = Array::make(columns, 1.7976931348623157e308)
  let maxima = Array::make(columns, -1.7976931348623157e308)
  for row in covariates {
    for j in 0..<columns {
      means[j] += row[j]
      if row[j] < minima[j] {
        minima[j] = row[j]
      }
      if row[j] > maxima[j] {
        maxima[j] = row[j]
      }
    }
  }
  for j in 0..<columns {
    means[j] /= covariates.length().to_double()
  }
  let standard_deviations = Array::make(columns, 0.0)
  for row in covariates {
    for j in 0..<columns {
      let difference = row[j] - means[j]
      standard_deviations[j] += difference * difference
    }
  }
  for j in 0..<columns {
    standard_deviations[j] = (standard_deviations[j] /
    covariates.length().to_double()).sqrt()
  }
  {
    columns,
    rows: covariates.length(),
    means,
    standard_deviations,
    minima,
    maxima,
  }
}

///|
pub fn column_values(
  covariates : Array[Array[Double]],
  column : Int,
) -> Array[Double] {
  let result = Array::new(capacity=covariates.length())
  for row in covariates {
    if column >= 0 && column < row.length() {
      result.push(row[column])
    }
  }
  result
}

///|
pub fn center_covariates(
  covariates : Array[Array[Double]],
) -> Array[Array[Double]] {
  let summary = summarize_covariates(covariates)
  let result : Array[Array[Double]] = Array::new(capacity=covariates.length())
  for row in covariates {
    let centered = Array::new(capacity=summary.columns)
    for j in 0..<summary.columns {
      centered.push(row[j] - summary.means[j])
    }
    result.push(centered)
  }
  result
}

///|
pub fn scale_covariates(
  covariates : Array[Array[Double]],
) -> Array[Array[Double]] {
  let summary = summarize_covariates(covariates)
  let result : Array[Array[Double]] = Array::new(capacity=covariates.length())
  for row in covariates {
    let scaled = Array::new(capacity=summary.columns)
    for j in 0..<summary.columns {
      let scale = if summary.standard_deviations[j] == 0.0 {
        1.0
      } else {
        summary.standard_deviations[j]
      }
      scaled.push((row[j] - summary.means[j]) / scale)
    }
    result.push(scaled)
  }
  result
}

///|
pub fn winsorize(
  values : Array[Double],
  lower_probability : Double,
  upper_probability : Double,
) -> Array[Double] {
  let lower = quantile(values, lower_probability)
  let upper = quantile(values, upper_probability)
  let result = Array::new(capacity=values.length())
  for value in values {
    result.push(clamp(value, lower, upper))
  }
  result
}

///|
pub fn rank_transform(values : Array[Double]) -> Array[Double] {
  let result = Array::make(values.length(), 0.0)
  for i in 0..<values.length() {
    let mut less = 0
    let mut equal = 0
    for other in values {
      if other < values[i] {
        less += 1
      }
      if other == values[i] {
        equal += 1
      }
    }
    let denominator = if values.length() == 0 {
      1.0
    } else {
      values.length().to_double()
    }
    result[i] = (less.to_double() + (equal - 1).to_double() / 2.0) / denominator
  }
  result
}

///|
pub fn polynomial_features(
  covariates : Array[Array[Double]],
  degree : Int,
) -> Array[Array[Double]] {
  if covariates.length() == 0 || degree <= 0 {
    return []
  }
  let result : Array[Array[Double]] = Array::new(capacity=covariates.length())
  for row in covariates {
    let expanded = Array::new()
    for value in row {
      expanded.push(value)
    }
    if degree >= 2 {
      for d in 2..<=degree {
        for value in row {
          expanded.push(@math.pow(value, d.to_double()))
        }
      }
    }
    result.push(expanded)
  }
  result
}

///|
pub fn interaction_features(
  covariates : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=covariates.length())
  for row in covariates {
    let expanded = row.copy()
    for i in 0..<row.length() {
      for j in (i + 1)..<row.length() {
        expanded.push(row[i] * row[j])
      }
    }
    result.push(expanded)
  }
  result
}

///|
pub struct DataSplit {
  training : CausalDataset
  validation : CausalDataset
}

///|
pub fn train_validation_split(
  dataset : CausalDataset,
  validation_fraction : Double,
  seed : UInt64,
) -> DataSplit {
  let fraction = clamp(validation_fraction, 0.0, 0.9)
  let order = Array::new(capacity=dataset.n())
  for i in 0..<dataset.n() {
    order.push(i)
  }
  let rng = RandomState::new(seed)
  for i in 0..<order.length() {
    let j = (rng.uniform() * order.length().to_double()).to_int()
    let temporary = order[i]
    order[i] = order[j]
    order[j] = temporary
  }
  let validation_size = (dataset.n().to_double() * fraction).to_int()
  let validation_indices = Array::new(capacity=validation_size)
  let training_indices = Array::new(capacity=dataset.n() - validation_size)
  for i in 0..<order.length() {
    if i < validation_size {
      validation_indices.push(order[i])
    } else {
      training_indices.push(order[i])
    }
  }
  {
    training: dataset.select(training_indices),
    validation: dataset.select(validation_indices),
  }
}
