///|
pub fn shuffled_indices(size : Int, seed : UInt64) -> Array[Int] {
  let n = if size < 0 { 0 } else { size }
  let result = Array::new(capacity=n)
  for i in 0..<n {
    result.push(i)
  }
  let rng = RandomState::new(seed)
  for i in 0..<n {
    let j = if n == 0 { 0 } else { (rng.uniform() * n.to_double()).to_int() }
    if n > 0 {
      let temporary = result[i]
      result[i] = result[j]
      result[j] = temporary
    }
  }
  result
}

///|
pub fn systematic_sample(
  size : Int,
  sample_size : Int,
  seed : UInt64,
) -> Array[Int] {
  let n = if size < 0 { 0 } else { size }
  let target = if sample_size < 0 {
    0
  } else if sample_size > n {
    n
  } else {
    sample_size
  }
  if target == 0 {
    return []
  }
  let step = n.to_double() / target.to_double()
  let rng = RandomState::new(seed)
  let start = rng.uniform() * step
  let result = Array::new(capacity=target)
  for i in 0..<target {
    result.push((start + i.to_double() * step).to_int())
  }
  result
}

///|
pub fn bootstrap_counts(
  size : Int,
  replicates : Int,
  seed : UInt64,
) -> Array[Array[Int]] {
  let n = if size < 0 { 0 } else { size }
  let result : Array[Array[Int]] = Array::new(
    capacity=if replicates > 0 { replicates } else { 0 },
  )
  let rng = RandomState::new(seed)
  for _ in 0..<replicates {
    let counts = Array::make(n, 0)
    for _ in 0..<n {
      if n > 0 {
        let index = (rng.uniform() * n.to_double()).to_int()
        counts[index] += 1
      }
    }
    result.push(counts)
  }
  result
}

///|
pub fn sample_variance(
  values : Array[Double],
  probabilities : Array[Double],
) -> Double {
  if values.length() == 0 || values.length() != probabilities.length() {
    return 0.0
  }
  let center = weighted_mean(values, probabilities)
  let mut result = 0.0
  let mut total = 0.0
  for i in 0..<values.length() {
    let difference = values[i] - center
    result += probabilities[i] * difference * difference
    total += probabilities[i]
  }
  if total == 0.0 {
    0.0
  } else {
    result / total
  }
}

///|
pub fn influence_function(values : Array[Double]) -> Array[Double] {
  let center = mean(values)
  let result = Array::new(capacity=values.length())
  for value in values {
    result.push(value - center)
  }
  result
}

///|
pub fn delete_one_mean(values : Array[Double]) -> Array[Double] {
  let result = Array::new(capacity=values.length())
  for excluded in 0..<values.length() {
    let remaining = Array::new(capacity=values.length() - 1)
    for i in 0..<values.length() {
      if i != excluded {
        remaining.push(values[i])
      }
    }
    result.push(mean(remaining))
  }
  result
}

///|
pub fn resample_by_weights(
  values : Array[Double],
  weights : Array[Double],
  draws : Int,
  seed : UInt64,
) -> Array[Double] {
  let result = Array::new(capacity=if draws > 0 { draws } else { 0 })
  if values.length() == 0 || values.length() != weights.length() {
    return result
  }
  let total = sum(weights)
  let rng = RandomState::new(seed)
  for _ in 0..<draws {
    let target = rng.uniform() * total
    let mut cumulative = 0.0
    let mut selected = values.length() - 1
    for i in 0..<values.length() {
      cumulative += weights[i]
      if cumulative >= target {
        selected = i
        break
      }
    }
    result.push(values[selected])
  }
  result
}
