///|
/// Distribution summary for a finite sample.
pub struct DistributionSummary {
  sample_size : Int
  mean : Double
  variance : Double
  standard_deviation : Double
  skewness : Double
  kurtosis : Double
  q01 : Double
  q25 : Double
  median : Double
  q75 : Double
  q99 : Double
}

///|
/// Probability interval for a scalar probability.
pub struct ProbabilityInterval {
  estimate : Double
  lower : Double
  upper : Double
  standard_error : Double
  sample_size : Int
}

///|
/// Computes a normal density.
pub fn normal_density(
  value : Double,
  mean_value? : Double = 0.0,
  standard_deviation? : Double = 1.0,
) -> Double {
  let scale = standard_deviation.abs().max(1.0e-12)
  let z = (value - mean_value) / scale
  @math.exp(-0.5 * z * z) / (2.5066282746310002 * scale)
}

///|
/// Computes a normal CDF with a rational approximation.
pub fn normal_distribution_cdf(
  value : Double,
  mean_value? : Double = 0.0,
  standard_deviation? : Double = 1.0,
) -> Double {
  let scale = standard_deviation.abs().max(1.0e-12)
  let z = (value - mean_value) / scale
  let sign = if z < 0.0 { -1.0 } else { 1.0 }
  let x = z.abs()
  let t = 1.0 / (1.0 + 0.2316419 * x)
  let polynomial = t *
    (
      0.319381530 +
      t *
      (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429)))
    )
  0.5 + sign * (0.5 - normal_density(x) * polynomial)
}

///|
/// Computes a Bernoulli log likelihood.
pub fn bernoulli_log_likelihood(
  outcomes : Array[Bool],
  probability : Double,
) -> Double {
  let p = safe_probability(probability)
  let mut result = 0.0
  for outcome in outcomes {
    result += if outcome { @math.ln(p) } else { @math.ln(1.0 - p) }
  }
  result
}

///|
/// Computes a Poisson log PMF for a non-negative count.
pub fn poisson_log_probability(count : Int, rate : Double) -> Double {
  if count < 0 || rate <= 0.0 {
    return -1.0e300
  }
  let mut factorial_log = 0.0
  for i in 2..<=count {
    factorial_log += @math.ln(i.to_double())
  }
  count.to_double() * @math.ln(rate) - rate - factorial_log
}

///|
/// Computes a Poisson PMF.
pub fn poisson_probability(count : Int, rate : Double) -> Double {
  @math.exp(poisson_log_probability(count, rate))
}

///|
/// Computes a binomial log PMF.
pub fn binomial_log_probability(
  successes : Int,
  trials : Int,
  probability : Double,
) -> Double {
  if successes < 0 || trials < 0 || successes > trials {
    return -1.0e300
  }
  let p = safe_probability(probability)
  let mut choose = 0.0
  for i in 1..<=successes {
    choose += @math.ln((trials - successes + i).to_double()) -
      @math.ln(i.to_double())
  }
  choose +
  successes.to_double() * @math.ln(p) +
  (trials - successes).to_double() * @math.ln(1.0 - p)
}

///|
/// Computes a binomial PMF.
pub fn binomial_probability(
  successes : Int,
  trials : Int,
  probability : Double,
) -> Double {
  @math.exp(binomial_log_probability(successes, trials, probability))
}

///|
/// Draws a Bernoulli sample with a reproducible RNG.
pub fn sample_bernoulli(
  probability : Double,
  size : Int,
  seed : UInt64,
) -> Array[Bool] {
  let rng = RandomState::new(seed)
  let result = Array::new(capacity=if size > 0 { size } else { 0 })
  for _ in 0..<size {
    result.push(rng.uniform() < clamp(probability, 0.0, 1.0))
  }
  result
}

///|
/// Draws Poisson counts using Knuth's algorithm for moderate rates.
pub fn sample_poisson(rate : Double, size : Int, seed : UInt64) -> Array[Int] {
  let rng = RandomState::new(seed)
  let result = Array::new(capacity=if size > 0 { size } else { 0 })
  let expected = rate.max(0.0)
  for _ in 0..<size {
    if expected == 0.0 {
      result.push(0)
    } else {
      let threshold = @math.exp(-expected)
      let mut product = 1.0
      let mut count = 0
      while product > threshold && count < 100000 {
        product *= rng.uniform().max(1.0e-12)
        count += 1
      }
      result.push(count - 1)
    }
  }
  result
}

///|
/// Summarizes finite sample moments and quantiles.
pub fn distribution_summary(values : Array[Double]) -> DistributionSummary {
  let observed = values.filter(fn(value) { is_finite(value) })
  let center = mean_or(observed, 0.0)
  let deviation = std_dev(observed)
  let mut third = 0.0
  let mut fourth = 0.0
  for value in observed {
    let z = if deviation == 0.0 { 0.0 } else { (value - center) / deviation }
    third += z * z * z
    fourth += z * z * z * z
  }
  let n = observed.length().to_double().max(1.0)
  {
    sample_size: observed.length(),
    mean: center,
    variance: variance(observed),
    standard_deviation: deviation,
    skewness: third / n,
    kurtosis: fourth / n - 3.0,
    q01: quantile(observed, 0.01),
    q25: quantile(observed, 0.25),
    median: quantile(observed, 0.5),
    q75: quantile(observed, 0.75),
    q99: quantile(observed, 0.99),
  }
}

///|
/// Computes a Wilson confidence interval for a proportion.
pub fn wilson_interval(
  successes : Int,
  trials : Int,
  confidence_level? : Double = 0.95,
) -> ProbabilityInterval {
  if trials <= 0 {
    return {
      estimate: 0.0,
      lower: 0.0,
      upper: 0.0,
      standard_error: 0.0,
      sample_size: 0,
    }
  }
  let estimate = clamp(successes.to_double() / trials.to_double(), 0.0, 1.0)
  let critical = if confidence_level >= 0.99 {
    2.5758293035489004
  } else if confidence_level >= 0.9 {
    1.959963984540054
  } else {
    1.6448536269514722
  }
  let denominator = 1.0 + critical * critical / trials.to_double()
  let center = (estimate + critical * critical / (2.0 * trials.to_double())) /
    denominator
  let radius = critical /
    denominator *
    (estimate * (1.0 - estimate) / trials.to_double() +
    critical * critical / (4.0 * trials.to_double() * trials.to_double())).sqrt()
  {
    estimate,
    lower: clamp(center - radius, 0.0, 1.0),
    upper: clamp(center + radius, 0.0, 1.0),
    standard_error: (estimate * (1.0 - estimate) / trials.to_double()).sqrt(),
    sample_size: trials,
  }
}

///|
/// Converts an integer count array to a normalized PMF.
pub fn normalize_counts(counts : Array[Int]) -> Array[Double] {
  let total = sum(counts.map(fn(value) { value.to_double() }))
  let result = Array::new(capacity=counts.length())
  for count in counts {
    result.push(if total == 0.0 { 0.0 } else { count.to_double() / total })
  }
  result
}

///|
/// Computes entropy of a probability vector.
pub fn probability_entropy(probabilities : Array[Double]) -> Double {
  let mut result = 0.0
  for probability in probabilities {
    let p = clamp(probability, 1.0e-12, 1.0)
    result -= p * @math.ln(p)
  }
  result
}

///|
/// Computes the Kullback-Leibler divergence between discrete PMFs.
pub fn discrete_kl_divergence(
  reference : Array[Double],
  current : Array[Double],
) -> Double {
  let n = reference.length().min(current.length())
  let mut result = 0.0
  for i in 0..<n {
    let p = reference[i].max(1.0e-12)
    let q = current[i].max(1.0e-12)
    result += p * @math.ln(p / q)
  }
  result
}

///|
/// Returns a compact distribution summary vector.
pub fn distribution_summary_vector(
  summary : DistributionSummary,
) -> Array[Double] {
  [
    summary.sample_size.to_double(),
    summary.mean,
    summary.variance,
    summary.standard_deviation,
    summary.skewness,
    summary.kurtosis,
    summary.q01,
    summary.q25,
    summary.median,
    summary.q75,
    summary.q99,
  ]
}
