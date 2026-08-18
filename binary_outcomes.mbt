///|
pub struct RiskMeasure {
  treated_risk : Double
  control_risk : Double
  risk_difference : Double
  risk_ratio : Double
  odds_ratio : Double
  number_needed_to_treat : Double
}

///|
fn binary_group_rate(
  outcomes : Array[Bool],
  treatment : Array[Bool],
  treated : Bool,
) -> Double {
  let mut successes = 0
  let mut total = 0
  let n = if outcomes.length() < treatment.length() {
    outcomes.length()
  } else {
    treatment.length()
  }
  for i in 0..<n {
    if treatment[i] == treated {
      total += 1
      if outcomes[i] {
        successes += 1
      }
    }
  }
  if total == 0 {
    0.0
  } else {
    successes.to_double() / total.to_double()
  }
}

///|
pub fn risk_measure(
  outcomes : Array[Bool],
  treatment : Array[Bool],
) -> RiskMeasure {
  let treated_risk = binary_group_rate(outcomes, treatment, true)
  let control_risk = binary_group_rate(outcomes, treatment, false)
  let difference = treated_risk - control_risk
  let ratio = if control_risk == 0.0 {
    1.0e300
  } else {
    treated_risk / control_risk
  }
  let treated_odds = if treated_risk == 1.0 {
    1.0e300
  } else {
    treated_risk / (1.0 - treated_risk)
  }
  let control_odds = if control_risk == 1.0 {
    1.0e300
  } else {
    control_risk / (1.0 - control_risk)
  }
  let odds_ratio = if control_odds == 0.0 {
    1.0e300
  } else {
    treated_odds / control_odds
  }
  {
    treated_risk,
    control_risk,
    risk_difference: difference,
    risk_ratio: ratio,
    odds_ratio,
    number_needed_to_treat: if difference == 0.0 {
      1.0e300
    } else {
      1.0 / difference.abs()
    },
  }
}

///|
pub fn binary_outcome_to_double(outcomes : Array[Bool]) -> Array[Double] {
  let result = Array::new(capacity=outcomes.length())
  for outcome in outcomes {
    result.push(if outcome { 1.0 } else { 0.0 })
  }
  result
}

///|
pub fn threshold_outcome(
  outcomes : Array[Double],
  threshold : Double,
) -> Array[Bool] {
  let result = Array::new(capacity=outcomes.length())
  for outcome in outcomes {
    result.push(outcome >= threshold)
  }
  result
}

///|
pub fn log_risk_ratio(measure : RiskMeasure) -> Double {
  if measure.risk_ratio <= 0.0 {
    0.0
  } else {
    @math.ln(measure.risk_ratio)
  }
}

///|
pub fn relative_risk_reduction(measure : RiskMeasure) -> Double {
  if measure.control_risk == 0.0 {
    0.0
  } else {
    1.0 - measure.treated_risk / measure.control_risk
  }
}
