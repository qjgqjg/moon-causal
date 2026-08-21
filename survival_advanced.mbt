///|
/// Discrete-time hazard table.
pub struct HazardTable {
  times : Array[Double]
  risk_set : Array[Int]
  events : Array[Int]
  hazard : Array[Double]
  survival : Array[Double]
}

///|
/// Competing-risk cumulative incidence curve.
pub struct CumulativeIncidence {
  times : Array[Double]
  incidence : Array[Double]
  event_type : Int
  final_incidence : Double
}

///|
/// Survival-weighted treatment effect.
pub struct SurvivalEffect {
  estimate : Double
  standard_error : Double
  restricted_mean_treated : Double
  restricted_mean_control : Double
  horizon : Double
  passes : Bool
}

///|
/// Time-varying hazard ratio summary.
pub struct HazardRatioPoint {
  time : Double
  hazard_treated : Double
  hazard_control : Double
  hazard_ratio : Double
  risk_set : Int
}

///|
/// Computes a discrete hazard table from observed event times and event indicators.
pub fn discrete_hazard_table(
  times : Array[Double],
  events : Array[Bool],
) -> HazardTable {
  let n = times.length().min(events.length())
  let unique : Array[Double] = Array::new()
  for i in 0..<n {
    if !unique.contains(times[i]) {
      unique.push(times[i])
    }
  }
  for i in 1..<unique.length() {
    let value = unique[i]
    let mut cursor = i
    while cursor > 0 && unique[cursor - 1] > value {
      unique[cursor] = unique[cursor - 1]
      cursor -= 1
    }
    unique[cursor] = value
  }
  let risk = Array::new(capacity=unique.length())
  let event_counts = Array::new(capacity=unique.length())
  let hazards = Array::new(capacity=unique.length())
  let survival = Array::new(capacity=unique.length())
  let mut current_survival = 1.0
  for time in unique {
    let mut at_risk = 0
    let mut event_count = 0
    for i in 0..<n {
      if times[i] >= time {
        at_risk += 1
      }
      if times[i] == time && events[i] {
        event_count += 1
      }
    }
    let hazard = if at_risk == 0 {
      0.0
    } else {
      event_count.to_double() / at_risk.to_double()
    }
    current_survival *= 1.0 - hazard
    risk.push(at_risk)
    event_counts.push(event_count)
    hazards.push(hazard)
    survival.push(current_survival)
  }
  {
    times: unique,
    risk_set: risk,
    events: event_counts,
    hazard: hazards,
    survival,
  }
}

///|
/// Computes a cumulative incidence curve for one competing event type.
pub fn competing_risk_incidence(
  times : Array[Double],
  event_types : Array[Int],
  target_event : Int,
) -> CumulativeIncidence {
  let n = times.length().min(event_types.length())
  let unique : Array[Double] = Array::new()
  for i in 0..<n {
    if !unique.contains(times[i]) {
      unique.push(times[i])
    }
  }
  for i in 1..<unique.length() {
    let value = unique[i]
    let mut cursor = i
    while cursor > 0 && unique[cursor - 1] > value {
      unique[cursor] = unique[cursor - 1]
      cursor -= 1
    }
    unique[cursor] = value
  }
  let incidence = Array::new(capacity=unique.length())
  let mut cumulative = 0.0
  let mut survival = 1.0
  for time in unique {
    let mut risk = 0
    let mut target = 0
    let mut any_event = 0
    for i in 0..<n {
      if times[i] >= time {
        risk += 1
      }
      if times[i] == time && event_types[i] != 0 {
        any_event += 1
      }
      if times[i] == time && event_types[i] == target_event {
        target += 1
      }
    }
    let target_hazard = if risk == 0 {
      0.0
    } else {
      target.to_double() / risk.to_double()
    }
    cumulative += survival * target_hazard
    let all_hazard = if risk == 0 {
      0.0
    } else {
      any_event.to_double() / risk.to_double()
    }
    survival *= 1.0 - all_hazard
    incidence.push(cumulative)
  }
  {
    times: unique,
    incidence,
    event_type: target_event,
    final_incidence: cumulative,
  }
}

///|
/// Computes restricted mean survival to a finite horizon.
pub fn restricted_mean_from_hazard(
  table : HazardTable,
  horizon : Double,
) -> Double {
  if table.times.length() == 0 {
    return 0.0
  }
  let mut area = 0.0
  let mut previous = 0.0
  let mut survival = 1.0
  for i in 0..<table.times.length() {
    let time = table.times[i]
    if time >= horizon {
      break
    }
    area += (time - previous) * survival
    survival = table.survival[i]
    previous = time
  }
  area + (horizon - previous).max(0.0) * survival
}

///|
/// Computes a restricted mean survival difference between treatment arms.
pub fn restricted_mean_survival_effect(
  times : Array[Double],
  events : Array[Bool],
  treatment : Array[Bool],
  horizon : Double,
) -> SurvivalEffect {
  let n = times.length().min(events.length()).min(treatment.length())
  let treated_times = Array::new()
  let treated_events = Array::new()
  let control_times = Array::new()
  let control_events = Array::new()
  for i in 0..<n {
    if treatment[i] {
      treated_times.push(times[i].min(horizon))
      treated_events.push(events[i] && times[i] <= horizon)
    } else {
      control_times.push(times[i].min(horizon))
      control_events.push(events[i] && times[i] <= horizon)
    }
  }
  let treated_table = discrete_hazard_table(treated_times, treated_events)
  let control_table = discrete_hazard_table(control_times, control_events)
  let treated_mean = restricted_mean_from_hazard(treated_table, horizon)
  let control_mean = restricted_mean_from_hazard(control_table, horizon)
  let estimate = treated_mean - control_mean
  let standard_error = if n == 0 {
    0.0
  } else {
    (treated_mean.abs() + control_mean.abs() + 1.0) / n.to_double().sqrt()
  }
  {
    estimate,
    standard_error,
    restricted_mean_treated: treated_mean,
    restricted_mean_control: control_mean,
    horizon,
    passes: treated_times.length() > 2 && control_times.length() > 2,
  }
}

///|
/// Computes a hazard ratio at each distinct time for a binary treatment.
pub fn time_varying_hazard_ratio(
  times : Array[Double],
  events : Array[Bool],
  treatment : Array[Bool],
) -> Array[HazardRatioPoint] {
  let n = times.length().min(events.length()).min(treatment.length())
  let unique : Array[Double] = Array::new()
  for i in 0..<n {
    if !unique.contains(times[i]) {
      unique.push(times[i])
    }
  }
  let result : Array[HazardRatioPoint] = Array::new(capacity=unique.length())
  for time in unique {
    let mut treated_risk = 0
    let mut control_risk = 0
    let mut treated_event = 0
    let mut control_event = 0
    for i in 0..<n {
      if times[i] >= time {
        if treatment[i] {
          treated_risk += 1
        } else {
          control_risk += 1
        }
      }
      if times[i] == time && events[i] {
        if treatment[i] {
          treated_event += 1
        } else {
          control_event += 1
        }
      }
    }
    let treated_hazard = if treated_risk == 0 {
      0.0
    } else {
      treated_event.to_double() / treated_risk.to_double()
    }
    let control_hazard = if control_risk == 0 {
      0.0
    } else {
      control_event.to_double() / control_risk.to_double()
    }
    result.push({
      time,
      hazard_treated: treated_hazard,
      hazard_control: control_hazard,
      hazard_ratio: if control_hazard == 0.0 {
        0.0
      } else {
        treated_hazard / control_hazard
      },
      risk_set: treated_risk + control_risk,
    })
  }
  result
}

///|
/// Computes inverse-probability-of-censoring weights.
pub fn censoring_weights(censoring_survival : Array[Double]) -> Array[Double] {
  let result = Array::new(capacity=censoring_survival.length())
  for survival in censoring_survival {
    result.push(1.0 / safe_probability(survival))
  }
  result
}

///|
/// Computes a survival-weighted restricted mean.
pub fn weighted_restricted_mean(
  times : Array[Double],
  survival : Array[Double],
  weights : Array[Double],
  horizon : Double,
) -> Double {
  let n = times.length().min(survival.length()).min(weights.length())
  let mut area = 0.0
  let mut previous = 0.0
  for i in 0..<n {
    if times[i] >= horizon {
      break
    }
    area += (times[i] - previous) * survival[i] * weights[i]
    previous = times[i]
  }
  let final_level = if n == 0 { 1.0 } else { survival[n - 1] * weights[n - 1] }
  area + (horizon - previous).max(0.0) * final_level
}

///|
/// Returns a compact survival effect summary vector.
pub fn survival_effect_summary(effect : SurvivalEffect) -> Array[Double] {
  [
    effect.estimate,
    effect.standard_error,
    effect.restricted_mean_treated,
    effect.restricted_mean_control,
    effect.horizon,
    if effect.passes {
      1.0
    } else {
      0.0
    },
  ]
}
