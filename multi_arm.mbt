///|
pub struct ArmEstimate {
  arm : Int
  count : Int
  mean_outcome : Double
  standard_error : Double
  difference_from_reference : Double
}

///|
fn arm_values(
  values : Array[Double],
  arms : Array[Int],
  arm : Int,
) -> Array[Double] {
  let result = Array::new()
  let n = if values.length() < arms.length() {
    values.length()
  } else {
    arms.length()
  }
  for i in 0..<n {
    if arms[i] == arm {
      result.push(values[i])
    }
  }
  result
}

///|
pub fn arm_levels(arms : Array[Int]) -> Array[Int] {
  let result = Array::new()
  for arm in arms {
    if !result.contains(arm) {
      result.push(arm)
    }
  }
  result
}

///|
pub fn arm_estimates(
  outcomes : Array[Double],
  arms : Array[Int],
  reference_arm : Int,
) -> Array[ArmEstimate] {
  let result : Array[ArmEstimate] = Array::new()
  let reference = mean(arm_values(outcomes, arms, reference_arm))
  for arm in arm_levels(arms) {
    let values = arm_values(outcomes, arms, arm)
    let standard_error = if values.length() == 0 {
      0.0
    } else {
      std_dev(values) / values.length().to_double().sqrt()
    }
    result.push({
      arm,
      count: values.length(),
      mean_outcome: mean(values),
      standard_error,
      difference_from_reference: mean(values) - reference,
    })
  }
  result
}

///|
pub fn pairwise_arm_difference(
  outcomes : Array[Double],
  arms : Array[Int],
  first_arm : Int,
  second_arm : Int,
) -> Estimate {
  let first = arm_values(outcomes, arms, first_arm)
  let second = arm_values(outcomes, arms, second_arm)
  let estimate = mean(first) - mean(second)
  let standard_error = if first.length() == 0 || second.length() == 0 {
    0.0
  } else {
    (variance(first) / first.length().to_double() +
    variance(second) / second.length().to_double()).sqrt()
  }
  Estimate::from_standard_error(
    estimate,
    standard_error,
    first.length() + second.length(),
    (first.length() + second.length()).to_double(),
    "pairwise arm difference",
  )
}

///|
pub fn arm_sample_counts(arms : Array[Int]) -> Array[(Int, Int)] {
  let result : Array[(Int, Int)] = Array::new()
  for arm in arm_levels(arms) {
    let mut count = 0
    for value in arms {
      if value == arm {
        count += 1
      }
    }
    result.push((arm, count))
  }
  result
}

///|
pub fn weighted_arm_means(
  outcomes : Array[Double],
  arms : Array[Int],
  weights : Array[Double],
) -> Array[ArmEstimate] {
  let result : Array[ArmEstimate] = Array::new()
  for arm in arm_levels(arms) {
    let values = Array::new()
    let selected_weights = Array::new()
    for i in 0..<outcomes.length() {
      if arms[i] == arm {
        values.push(outcomes[i])
        selected_weights.push(weights[i])
      }
    }
    result.push({
      arm,
      count: values.length(),
      mean_outcome: weighted_mean(values, selected_weights),
      standard_error: 0.0,
      difference_from_reference: 0.0,
    })
  }
  result
}
