///|
fn resample_indices(rng : RandomState, n : Int) -> Array[Int] {
  let result = Array::new(capacity=n)
  if n == 0 {
    return result
  }
  for _ in 0..<n {
    result.push((rng.uniform() * n.to_double()).to_int())
  }
  result
}

///|
fn difference_in_means(
  outcomes : Array[Double],
  treatment : Array[Bool],
) -> Double {
  mean(group_values(outcomes, treatment, true)) -
  mean(group_values(outcomes, treatment, false))
}

///|
fn resampled_values(
  values : Array[Double],
  indices : Array[Int],
) -> Array[Double] {
  let result = Array::new(capacity=indices.length())
  for index in indices {
    result.push(values[index])
  }
  result
}

///|
fn resampled_bool(values : Array[Bool], indices : Array[Int]) -> Array[Bool] {
  let result = Array::new(capacity=indices.length())
  for index in indices {
    result.push(values[index])
  }
  result
}

///|
/// Builds a reproducible nonparametric bootstrap distribution for the difference in means.
pub fn bootstrap_difference_in_means(
  outcomes : Array[Double],
  treatment : Array[Bool],
  replicates? : Int = 1000,
  seed? : UInt64 = 20260818UL,
) -> BootstrapSummary {
  let point = difference_in_means(outcomes, treatment)
  let rng = RandomState::new(seed)
  let samples : Array[Double] = Array::new(capacity=replicates)
  let mut successful = 0
  for _ in 0..<replicates {
    let indices = resample_indices(rng, outcomes.length())
    let sample_outcomes = resampled_values(outcomes, indices)
    let sample_treatment = resampled_bool(treatment, indices)
    let mut treated = 0
    let mut control = 0
    for value in sample_treatment {
      if value {
        treated += 1
      } else {
        control += 1
      }
    }
    if treated > 0 && control > 0 {
      samples.push(difference_in_means(sample_outcomes, sample_treatment))
      successful += 1
    }
  }
  let standard_error = std_dev(samples, sample=true)
  {
    point_estimate: point,
    standard_error,
    lower: quantile(samples, 0.025),
    upper: quantile(samples, 0.975),
    replicates,
    successful_replicates: successful,
  }
}

///|
/// Calculates a randomization p-value under the sharp null using treatment permutations.
pub fn permutation_p_value(
  outcomes : Array[Double],
  treatment : Array[Bool],
  permutations? : Int = 1000,
  seed? : UInt64 = 20260818UL,
) -> Double {
  let observed = difference_in_means(outcomes, treatment).abs()
  let rng = RandomState::new(seed)
  let mut exceedances = 0
  for _ in 0..<permutations {
    let shuffled = treatment.copy()
    for i in 0..<shuffled.length() {
      let j = (rng.uniform() * shuffled.length().to_double()).to_int()
      let temporary = shuffled[i]
      shuffled[i] = shuffled[j]
      shuffled[j] = temporary
    }
    if difference_in_means(outcomes, shuffled).abs() >= observed {
      exceedances += 1
    }
  }
  (exceedances + 1).to_double() / (permutations + 1).to_double()
}

///|
/// Computes the leave-one-out influence values for the difference in means.
pub fn jackknife_influence(
  outcomes : Array[Double],
  treatment : Array[Bool],
) -> Array[Double] {
  let result = Array::new(capacity=outcomes.length())
  let full = difference_in_means(outcomes, treatment)
  for excluded in 0..<outcomes.length() {
    let y = Array::new(capacity=outcomes.length() - 1)
    let t = Array::new(capacity=outcomes.length() - 1)
    for i in 0..<outcomes.length() {
      if i != excluded {
        y.push(outcomes[i])
        t.push(treatment[i])
      }
    }
    result.push(full - difference_in_means(y, t))
  }
  result
}
