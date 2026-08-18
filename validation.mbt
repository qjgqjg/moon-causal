///|
/// Returns whether a number is neither NaN nor infinite.
pub fn is_finite(value : Double) -> Bool {
  value == value && value.abs() < 1.7976931348623157e308
}

///|
pub fn clamp(value : Double, lower : Double, upper : Double) -> Double {
  if value < lower {
    lower
  } else if value > upper {
    upper
  } else {
    value
  }
}

///|
pub fn safe_probability(value : Double, epsilon? : Double = 1.0e-8) -> Double {
  clamp(value, epsilon, 1.0 - epsilon)
}

///|
pub fn same_length(a : Array[Double], b : Array[Double]) -> Bool {
  a.length() == b.length()
}

///|
pub fn same_length_bool(a : Array[Bool], b : Array[Double]) -> Bool {
  a.length() == b.length()
}

///|
pub fn has_variation(values : Array[Double]) -> Bool {
  if values.length() < 2 {
    return false
  }
  let first = values[0]
  for value in values[1:] {
    if value != first {
      return true
    }
  }
  false
}

///|
pub fn finite_array(values : Array[Double]) -> Bool {
  for value in values {
    if !is_finite(value) {
      return false
    }
  }
  true
}

///|
pub fn finite_matrix(values : Array[Array[Double]]) -> Bool {
  for row in values {
    if !finite_array(row) {
      return false
    }
  }
  true
}

///|
pub fn mean_or(values : Array[Double], fallback : Double) -> Double {
  if values.length() == 0 {
    fallback
  } else {
    mean(values)
  }
}

///|
pub fn weighted_mean(values : Array[Double], weights : Array[Double]) -> Double {
  if values.length() == 0 || values.length() != weights.length() {
    return 0.0
  }
  let mut numerator = 0.0
  let mut denominator = 0.0
  for i in 0..<values.length() {
    numerator += values[i] * weights[i]
    denominator += weights[i]
  }
  if denominator == 0.0 {
    0.0
  } else {
    numerator / denominator
  }
}

///|
pub fn weighted_variance(
  values : Array[Double],
  weights : Array[Double],
) -> Double {
  if values.length() == 0 || values.length() != weights.length() {
    return 0.0
  }
  let center = weighted_mean(values, weights)
  let mut numerator = 0.0
  let mut denominator = 0.0
  for i in 0..<values.length() {
    let difference = values[i] - center
    numerator += weights[i] * difference * difference
    denominator += weights[i]
  }
  if denominator == 0.0 {
    0.0
  } else {
    numerator / denominator
  }
}

///|
pub fn effective_sample_size(weights : Array[Double]) -> Double {
  if weights.length() == 0 {
    return 0.0
  }
  let mut sum_weights = 0.0
  let mut sum_squares = 0.0
  for weight in weights {
    sum_weights += weight
    sum_squares += weight * weight
  }
  if sum_squares == 0.0 {
    0.0
  } else {
    sum_weights * sum_weights / sum_squares
  }
}

///|
pub fn quantile(values : Array[Double], probability : Double) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let sorted = values.copy()
  for i in 1..<sorted.length() {
    let value = sorted[i]
    let mut j = i
    while j > 0 && sorted[j - 1] > value {
      sorted[j] = sorted[j - 1]
      j -= 1
    }
    sorted[j] = value
  }
  let p = clamp(probability, 0.0, 1.0)
  let position = p * (sorted.length() - 1).to_double()
  let lower = position.to_int()
  let upper = if lower + 1 < sorted.length() { lower + 1 } else { lower }
  sorted[lower] +
  (sorted[upper] - sorted[lower]) * (position - lower.to_double())
}

///|
pub fn standardize(values : Array[Double]) -> Array[Double] {
  let center = mean(values)
  let scale = std_dev(values, sample=false)
  let result = Array::new(capacity=values.length())
  for value in values {
    result.push(if scale == 0.0 { 0.0 } else { (value - center) / scale })
  }
  result
}
