///|
/// Dose-response curve point with uncertainty.
pub struct AdvancedDosePoint {
  dose : Double
  count : Int
  mean_outcome : Double
  standard_error : Double
  lower : Double
  upper : Double
  weight_sum : Double
}

///|
/// Dose-response curve diagnostic.
pub struct DoseResponseAudit {
  points : Array[AdvancedDosePoint]
  monotonicity_score : Double
  maximum_gap : Double
  effective_sample_size : Double
  passes : Bool
}

///|
/// Builds weighted dose-response points over a grid.
pub fn advanced_dose_response(
  dose : Array[Double],
  outcome : Array[Double],
  weights : Array[Double],
  grid : Array[Double],
  bandwidth : Double,
) -> Array[AdvancedDosePoint] {
  let n = dose.length().min(outcome.length()).min(weights.length())
  let result : Array[AdvancedDosePoint] = Array::new(capacity=grid.length())
  for center in grid {
    let selected = Array::new()
    let selected_weights = Array::new()
    for i in 0..<n {
      let distance = (dose[i] - center).abs()
      if distance <= bandwidth.abs() {
        selected.push(outcome[i])
        selected_weights.push(
          weights[i] * (1.0 - distance / bandwidth.abs().max(1.0e-12)),
        )
      }
    }
    let estimate = weighted_mean(selected, selected_weights)
    let effective = effective_sample_size(selected_weights)
    let error = if effective == 0.0 {
      0.0
    } else {
      (weighted_variance(selected, selected_weights) / effective).sqrt()
    }
    result.push({
      dose: center,
      count: selected.length(),
      mean_outcome: estimate,
      standard_error: error,
      lower: estimate - 1.96 * error,
      upper: estimate + 1.96 * error,
      weight_sum: sum(selected_weights),
    })
  }
  result
}

///|
/// Scores local monotonicity of a dose-response curve.
pub fn dose_monotonicity(points : Array[AdvancedDosePoint]) -> Double {
  if points.length() < 2 {
    return 1.0
  }
  let mut concordant = 0
  let mut total = 0
  for i in 1..<points.length() {
    total += 1
    if points[i].mean_outcome >= points[i - 1].mean_outcome {
      concordant += 1
    }
  }
  concordant.to_double() / total.to_double()
}

///|
/// Computes the maximum adjacent uncertainty gap.
pub fn dose_maximum_gap(points : Array[AdvancedDosePoint]) -> Double {
  let mut result = 0.0
  for i in 1..<points.length() {
    let gap = (points[i].mean_outcome - points[i - 1].mean_outcome).abs()
    if gap > result {
      result = gap
    }
  }
  result
}

///|
/// Audits dose-response support and smoothness.
pub fn audit_dose_response(
  points : Array[AdvancedDosePoint],
  minimum_effective_sample_size? : Double = 5.0,
) -> DoseResponseAudit {
  let mut effective = 0.0
  for point in points {
    effective += point.weight_sum
  }
  let monotonicity = dose_monotonicity(points)
  let gap = dose_maximum_gap(points)
  {
    points,
    monotonicity_score: monotonicity,
    maximum_gap: gap,
    effective_sample_size: effective,
    passes: points.length() > 1 && effective >= minimum_effective_sample_size,
  }
}

///|
/// Integrates a dose-response curve by the trapezoidal rule.
pub fn integrate_dose_response(points : Array[AdvancedDosePoint]) -> Double {
  let mut result = 0.0
  for i in 1..<points.length() {
    result += (points[i].dose - points[i - 1].dose) *
      (points[i].mean_outcome + points[i - 1].mean_outcome) /
      2.0
  }
  result
}

///|
/// Finds the dose maximizing estimated outcome.
pub fn optimal_dose(points : Array[AdvancedDosePoint]) -> Double {
  if points.length() == 0 {
    return 0.0
  }
  let mut best = 0
  for i in 1..<points.length() {
    if points[i].mean_outcome > points[best].mean_outcome {
      best = i
    }
  }
  points[best].dose
}

///|
/// Computes a finite-difference marginal effect along a dose grid.
pub fn dose_marginal_effect(
  points : Array[AdvancedDosePoint],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  for i in 1..<points.length() {
    let distance = points[i].dose - points[i - 1].dose
    result.push([
      points[i].dose,
      if distance == 0.0 {
        0.0
      } else {
        (points[i].mean_outcome - points[i - 1].mean_outcome) / distance
      },
    ])
  }
  result
}

///|
/// Returns a compact dose-response summary vector.
pub fn dose_response_summary(audit : DoseResponseAudit) -> Array[Double] {
  [
    audit.points.length().to_double(),
    audit.monotonicity_score,
    audit.maximum_gap,
    audit.effective_sample_size,
    if audit.passes {
      1.0
    } else {
      0.0
    },
  ]
}
