///|
/// Minimum-cell and privacy policy for released aggregates.
pub struct PrivacyPolicy {
  minimum_group_size : Int
  maximum_group_count : Int
  noise_scale : Double
  suppress_small_cells : Bool
  seed : UInt64
}

///|
/// Privacy-safe aggregate output.
pub struct SafeAggregate {
  key : Int
  count : Int
  released_value : Double
  suppressed : Bool
  noisy : Bool
}

///|
/// K-anonymity and cell-size audit.
pub struct KAnonymityAudit {
  groups : Int
  small_groups : Int
  minimum_count : Int
  suppressed_fraction : Double
  passes : Bool
}

///|
/// Creates a privacy policy with validated limits.
pub fn privacy_policy(
  minimum_group_size? : Int = 5,
  maximum_group_count? : Int = 10000,
  noise_scale? : Double = 0.0,
  suppress_small_cells? : Bool = true,
  seed? : UInt64 = 20260819,
) -> PrivacyPolicy {
  {
    minimum_group_size: minimum_group_size.max(1),
    maximum_group_count: maximum_group_count.max(1),
    noise_scale: noise_scale.max(0.0),
    suppress_small_cells,
    seed,
  }
}

///|
/// Applies count suppression to grouped numeric values.
pub fn suppress_small_groups(
  keys : Array[Int],
  counts : Array[Int],
  values : Array[Double],
  policy : PrivacyPolicy,
) -> Array[SafeAggregate] {
  let n = keys.length().min(counts.length()).min(values.length())
  let result : Array[SafeAggregate] = Array::new(
    capacity=n.min(policy.maximum_group_count),
  )
  let rng = RandomState::new(policy.seed)
  for i in 0..<n.min(policy.maximum_group_count) {
    let suppressed = policy.suppress_small_cells &&
      counts[i] < policy.minimum_group_size
    let noise = if policy.noise_scale == 0.0 {
      0.0
    } else {
      policy.noise_scale * privacy_laplace(rng.uniform())
    }
    result.push({
      key: keys[i],
      count: counts[i],
      released_value: if suppressed {
        0.0
      } else {
        values[i] + noise
      },
      suppressed,
      noisy: !suppressed && policy.noise_scale > 0.0,
    })
  }
  result
}

///|
fn privacy_laplace(uniform : Double) -> Double {
  let centered = clamp(uniform, 1.0e-12, 1.0 - 1.0e-12) - 0.5
  if centered >= 0.0 {
    -@math.ln(1.0 - 2.0 * centered)
  } else {
    @math.ln(1.0 + 2.0 * centered)
  }
}

///|
/// Audits group sizes against a k-anonymity threshold.
pub fn audit_k_anonymity(
  counts : Array[Int],
  policy : PrivacyPolicy,
) -> KAnonymityAudit {
  let mut small = 0
  let mut minimum = if counts.length() == 0 { 0 } else { counts[0] }
  for count in counts {
    if count < policy.minimum_group_size {
      small += 1
    }
    if count < minimum {
      minimum = count
    }
  }
  {
    groups: counts.length(),
    small_groups: small,
    minimum_count: minimum,
    suppressed_fraction: if counts.length() == 0 {
      0.0
    } else {
      small.to_double() / counts.length().to_double()
    },
    passes: small == 0 && counts.length() <= policy.maximum_group_count,
  }
}

///|
/// Caps an aggregate contribution before release.
pub fn cap_contribution(
  value : Double,
  lower : Double,
  upper : Double,
) -> Double {
  clamp(value, lower.min(upper), upper.max(lower))
}

///|
/// Applies contribution caps to a numeric vector.
pub fn cap_contributions(
  values : Array[Double],
  lower : Double,
  upper : Double,
) -> Array[Double] {
  let result = Array::new(capacity=values.length())
  for value in values {
    result.push(cap_contribution(value, lower, upper))
  }
  result
}

///|
/// Computes a bounded mean for privacy-safe release.
pub fn bounded_mean(
  values : Array[Double],
  lower : Double,
  upper : Double,
) -> Double {
  mean_or(cap_contributions(values, lower, upper), 0.0)
}

///|
/// Computes the maximum single-row influence on a bounded mean.
pub fn bounded_mean_sensitivity(
  sample_size : Int,
  lower : Double,
  upper : Double,
) -> Double {
  if sample_size <= 0 {
    0.0
  } else {
    (upper - lower).abs() / sample_size.to_double()
  }
}

///|
/// Returns a compact privacy audit vector.
pub fn privacy_summary(
  audit : KAnonymityAudit,
  policy : PrivacyPolicy,
) -> Array[Double] {
  [
    audit.groups.to_double(),
    audit.small_groups.to_double(),
    audit.minimum_count.to_double(),
    audit.suppressed_fraction,
    policy.minimum_group_size.to_double(),
    policy.noise_scale,
    if audit.passes {
      1.0
    } else {
      0.0
    },
  ]
}
