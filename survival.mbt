///|
pub struct SurvivalPoint {
  time : Double
  at_risk : Int
  events : Int
  survival : Double
  standard_error : Double
}

///|
pub fn kaplan_meier(
  times : Array[Double],
  events : Array[Bool],
) -> Array[SurvivalPoint] {
  let order = Array::new(capacity=times.length())
  for i in 0..<times.length() {
    order.push(i)
  }
  for i in 1..<order.length() {
    let key = order[i]
    let mut j = i
    while j > 0 && times[order[j - 1]] > times[key] {
      order[j] = order[j - 1]
      j -= 1
    }
    order[j] = key
  }
  let unique_times = Array::new()
  for index in order {
    if unique_times.length() == 0 ||
      unique_times[unique_times.length() - 1] != times[index] {
      unique_times.push(times[index])
    }
  }
  let result : Array[SurvivalPoint] = Array::new()
  let mut survival = 1.0
  let mut variance_sum = 0.0
  for time in unique_times {
    let mut at_risk = 0
    let mut event_count = 0
    for i in 0..<times.length() {
      if times[i] >= time {
        at_risk += 1
      }
      if times[i] == time && events[i] {
        event_count += 1
      }
    }
    if event_count > 0 && at_risk > 0 {
      survival *= (at_risk - event_count).to_double() / at_risk.to_double()
      if at_risk > event_count {
        variance_sum += event_count.to_double() /
          (at_risk * (at_risk - event_count)).to_double()
      }
      result.push({
        time,
        at_risk,
        events: event_count,
        survival,
        standard_error: (survival * survival * variance_sum).sqrt(),
      })
    }
  }
  result
}

///|
pub fn restricted_mean_survival(
  curve : Array[SurvivalPoint],
  horizon : Double,
) -> Double {
  if curve.length() == 0 {
    return horizon
  }
  let mut last_time = 0.0
  let mut last_survival = 1.0
  let mut total = 0.0
  for point in curve {
    if point.time >= horizon {
      break
    }
    total += (point.time - last_time) * last_survival
    last_time = point.time
    last_survival = point.survival
  }
  total + (horizon - last_time).max(0.0) * last_survival
}

///|
pub fn survival_at(curve : Array[SurvivalPoint], time : Double) -> Double {
  let mut result = 1.0
  for point in curve {
    if point.time <= time {
      result = point.survival
    }
  }
  result
}

///|
pub fn hazard_ratio_from_curves(
  treated : Array[SurvivalPoint],
  control : Array[SurvivalPoint],
  horizon : Double,
) -> Double {
  let treated_survival = survival_at(treated, horizon)
  let control_survival = survival_at(control, horizon)
  if treated_survival <= 0.0 || control_survival <= 0.0 {
    1.0
  } else {
    @math.ln(treated_survival) / @math.ln(control_survival)
  }
}
