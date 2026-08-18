///|
pub struct MediationResult {
  total_effect : Double
  direct_effect : Double
  indirect_effect : Double
  mediated_fraction : Double
}

///|
/// Decomposes an effect using model-based counterfactual predictions.
pub fn mediation_decomposition(
  outcome_treated_mediator_treated : Array[Double],
  outcome_treated_mediator_control : Array[Double],
  outcome_control_mediator_treated : Array[Double],
  outcome_control_mediator_control : Array[Double],
) -> MediationResult {
  let n = outcome_treated_mediator_treated.length()
  let n1 = if n < outcome_treated_mediator_control.length() {
    n
  } else {
    outcome_treated_mediator_control.length()
  }
  let n2 = if n1 < outcome_control_mediator_treated.length() {
    n1
  } else {
    outcome_control_mediator_treated.length()
  }
  let count = if n2 < outcome_control_mediator_control.length() {
    n2
  } else {
    outcome_control_mediator_control.length()
  }
  let mut total = 0.0
  let mut direct = 0.0
  let mut indirect = 0.0
  for i in 0..<count {
    total += outcome_treated_mediator_treated[i] -
      outcome_control_mediator_control[i]
    direct += outcome_treated_mediator_control[i] -
      outcome_control_mediator_control[i]
    indirect += outcome_treated_mediator_treated[i] -
      outcome_treated_mediator_control[i]
  }
  let denominator = if count == 0 { 1.0 } else { count.to_double() }
  let total_effect = total / denominator
  let direct_effect = direct / denominator
  let indirect_effect = indirect / denominator
  {
    total_effect,
    direct_effect,
    indirect_effect,
    mediated_fraction: if total_effect == 0.0 {
      0.0
    } else {
      indirect_effect / total_effect
    },
  }
}

///|
pub fn mediation_sensitivity(
  result : MediationResult,
  confounding_strength : Double,
) -> MediationResult {
  let strength = clamp(confounding_strength, 0.0, 1.0)
  let indirect_effect = result.indirect_effect * (1.0 - strength)
  let direct_effect = result.total_effect - indirect_effect
  {
    total_effect: result.total_effect,
    direct_effect,
    indirect_effect,
    mediated_fraction: if result.total_effect == 0.0 {
      0.0
    } else {
      indirect_effect / result.total_effect
    },
  }
}

///|
pub fn mediation_balance_check(
  result : MediationResult,
  tolerance : Double,
) -> Bool {
  (result.total_effect - result.direct_effect - result.indirect_effect).abs() <=
  tolerance
}
