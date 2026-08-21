///|
/// Utility rule for turning causal estimates into actions.
pub struct DecisionRule {
  treatment_cost : Double
  treatment_benefit : Double
  harm_cost : Double
  minimum_net_benefit : Double
  threshold : Double
}

///|
/// Decision audit for a treatment policy.
pub struct DecisionAudit {
  selected : Int
  total : Int
  selected_fraction : Double
  expected_net_benefit : Double
  false_positive_cost : Double
  false_negative_cost : Double
  passes : Bool
}

///|
/// Policy decision for one unit.
pub struct UnitDecision {
  index : Int
  effect : Double
  net_benefit : Double
  treat : Bool
  confidence : Double
}

///|
/// Creates a decision rule from utility parameters.
pub fn decision_rule(
  treatment_cost : Double,
  treatment_benefit : Double,
  harm_cost? : Double = 0.0,
  minimum_net_benefit? : Double = 0.0,
) -> DecisionRule {
  let threshold = if treatment_benefit + harm_cost == 0.0 {
    1.0e300
  } else {
    (treatment_cost + harm_cost) / (treatment_benefit + harm_cost)
  }
  {
    treatment_cost,
    treatment_benefit,
    harm_cost,
    minimum_net_benefit,
    threshold,
  }
}

///|
/// Converts an effect estimate into expected net benefit.
pub fn expected_net_benefit(effect : Double, rule : DecisionRule) -> Double {
  let harm = if effect < 0.0 { rule.harm_cost * effect.abs() } else { 0.0 }
  effect * rule.treatment_benefit - rule.treatment_cost - harm
}

///|
/// Decides treatment from an estimated individual effect.
pub fn decide_unit(
  index : Int,
  effect : Double,
  standard_error : Double,
  rule : DecisionRule,
) -> UnitDecision {
  let net_benefit = expected_net_benefit(effect, rule)
  let treat = net_benefit >= rule.minimum_net_benefit
  let confidence = if standard_error == 0.0 {
    1.0
  } else {
    clamp(effect.abs() / standard_error.abs() / 3.0, 0.0, 1.0)
  }
  { index, effect, net_benefit, treat, confidence }
}

///|
/// Applies a decision rule to a vector of effects.
pub fn decide_policy(
  effects : Array[Double],
  standard_errors : Array[Double],
  rule : DecisionRule,
) -> Array[UnitDecision] {
  let result : Array[UnitDecision] = Array::new(capacity=effects.length())
  for i in 0..<effects.length() {
    result.push(
      decide_unit(
        i,
        effects[i],
        if i < standard_errors.length() {
          standard_errors[i]
        } else {
          0.0
        },
        rule,
      ),
    )
  }
  result
}

///|
/// Audits a policy against observed outcomes and treatment assignment.
pub fn audit_decision_policy(
  decisions : Array[UnitDecision],
  treatment : Array[Bool],
  outcomes : Array[Double],
  rule : DecisionRule,
) -> DecisionAudit {
  let n = decisions.length().min(treatment.length()).min(outcomes.length())
  let mut selected = 0
  let mut net_benefit = 0.0
  let mut false_positive = 0.0
  let mut false_negative = 0.0
  for i in 0..<n {
    if decisions[i].treat {
      selected += 1
      net_benefit += expected_net_benefit(outcomes[i], rule)
    }
    if decisions[i].treat && outcomes[i] < 0.0 {
      false_positive += outcomes[i].abs()
    }
    if !decisions[i].treat && outcomes[i] > 0.0 {
      false_negative += outcomes[i]
    }
  }
  {
    selected,
    total: n,
    selected_fraction: if n == 0 {
      0.0
    } else {
      selected.to_double() / n.to_double()
    },
    expected_net_benefit: if n == 0 {
      0.0
    } else {
      net_benefit / n.to_double()
    },
    false_positive_cost: false_positive,
    false_negative_cost: false_negative,
    passes: n > 0 &&
    net_benefit >= rule.minimum_net_benefit * selected.to_double(),
  }
}

///|
/// Computes net benefit for a range of treatment thresholds.
pub fn decision_curve(
  effects : Array[Double],
  rule : DecisionRule,
  thresholds : Array[Double],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=thresholds.length())
  for threshold in thresholds {
    let mut selected = 0
    let mut value = 0.0
    for effect in effects {
      if effect >= threshold {
        selected += 1
        value += expected_net_benefit(effect, rule)
      }
    }
    result.push([
      threshold,
      selected.to_double(),
      if effects.length() == 0 {
        0.0
      } else {
        value / effects.length().to_double()
      },
    ])
  }
  result
}

///|
/// Returns a compact decision summary vector.
pub fn decision_summary(audit : DecisionAudit) -> Array[Double] {
  [
    audit.selected.to_double(),
    audit.total.to_double(),
    audit.selected_fraction,
    audit.expected_net_benefit,
    audit.false_positive_cost,
    audit.false_negative_cost,
    if audit.passes {
      1.0
    } else {
      0.0
    },
  ]
}
