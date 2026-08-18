///|
pub fn median(values : Array[Double]) -> Double {
  quantile(values, 0.5)
}

///|
pub fn median_absolute_deviation(values : Array[Double]) -> Double {
  let center = median(values)
  let deviations = Array::new(capacity=values.length())
  for value in values {
    deviations.push((value - center).abs())
  }
  median(deviations)
}

///|
pub fn trimmed_mean(values : Array[Double], trim_fraction : Double) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let sorted = values.copy()
  for i in 1..<sorted.length() {
    let key = sorted[i]
    let mut j = i
    while j > 0 && sorted[j - 1] > key {
      sorted[j] = sorted[j - 1]
      j -= 1
    }
    sorted[j] = key
  }
  let fraction = clamp(trim_fraction, 0.0, 0.49)
  let cut = (sorted.length().to_double() * fraction).to_int()
  let remaining = Array::new()
  for i in cut..<(sorted.length() - cut) {
    remaining.push(sorted[i])
  }
  mean(remaining)
}

///|
pub fn winsorized_mean(
  values : Array[Double],
  trim_fraction : Double,
) -> Double {
  mean(winsorize(values, trim_fraction, 1.0 - clamp(trim_fraction, 0.0, 0.49)))
}

///|
pub fn weighted_quantile(
  values : Array[Double],
  weights : Array[Double],
  probability : Double,
) -> Double {
  if values.length() == 0 || values.length() != weights.length() {
    return 0.0
  }
  let order = Array::new(capacity=values.length())
  for i in 0..<values.length() {
    order.push(i)
  }
  for i in 1..<order.length() {
    let key = order[i]
    let mut j = i
    while j > 0 && values[order[j - 1]] > values[key] {
      order[j] = order[j - 1]
      j -= 1
    }
    order[j] = key
  }
  let target = clamp(probability, 0.0, 1.0) * sum(weights)
  let mut cumulative = 0.0
  for index in order {
    cumulative += weights[index]
    if cumulative >= target {
      return values[index]
    }
  }
  values[order[order.length() - 1]]
}

///|
pub fn empirical_cdf(
  values : Array[Double],
  points : Array[Double],
) -> Array[Double] {
  let result = Array::new(capacity=points.length())
  for point in points {
    let mut count = 0
    for value in values {
      if value <= point {
        count += 1
      }
    }
    result.push(
      if values.length() == 0 {
        0.0
      } else {
        count.to_double() / values.length().to_double()
      },
    )
  }
  result
}

///|
pub fn standardized_mean_difference_from_groups(
  treated : Array[Double],
  control : Array[Double],
) -> Double {
  let pooled = (variance(treated) + variance(control)) / 2.0
  if pooled <= 0.0 {
    0.0
  } else {
    (mean(treated) - mean(control)) / pooled.sqrt()
  }
}

///|
pub fn robust_z_scores(values : Array[Double]) -> Array[Double] {
  let center = median(values)
  let scale = median_absolute_deviation(values) * 1.4826
  let result = Array::new(capacity=values.length())
  for value in values {
    result.push(if scale == 0.0 { 0.0 } else { (value - center) / scale })
  }
  result
}

///|
pub fn outlier_indices(
  values : Array[Double],
  threshold : Double,
) -> Array[Int] {
  let z = robust_z_scores(values)
  let result = Array::new()
  for i in 0..<z.length() {
    if z[i].abs() > threshold {
      result.push(i)
    }
  }
  result
}
