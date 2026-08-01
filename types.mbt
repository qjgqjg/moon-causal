///|
/// Calculates the sum of an array of doubles.
pub fn sum(arr : Array[Double]) -> Double {
  let mut s = 0.0
  for val in arr {
    s += val
  }
  s
}

///|
/// Calculates the mean of an array of doubles.
pub fn mean(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  sum(arr) / arr.length().to_double()
}

///|
/// Calculates the variance of an array of doubles.
pub fn variance(arr : Array[Double], sample? : Bool = true) -> Double {
  let n = arr.length()
  if n <= 1 {
    return 0.0
  }
  let m = mean(arr)
  let mut sum_sq = 0.0
  for val in arr {
    let diff = val - m
    sum_sq += diff * diff
  }
  let denominator = if sample { n - 1 } else { n }
  sum_sq / denominator.to_double()
}

///|
/// Calculates the standard deviation of an array of doubles.
pub fn std_dev(arr : Array[Double], sample? : Bool = true) -> Double {
  let v = variance(arr, sample~)
  v.sqrt()
}
