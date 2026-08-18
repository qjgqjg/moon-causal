///|
pub struct DifferenceInDifferencesResult {
  treated_pre_mean : Double
  treated_post_mean : Double
  control_pre_mean : Double
  control_post_mean : Double
  treatment_change : Double
  control_change : Double
  estimate : Double
  standard_error : Double
}

///|
fn select_panel(
  values : Array[Double],
  group : Array[Bool],
  period : Array[Bool],
  selected_group : Bool,
  selected_period : Bool,
) -> Array[Double] {
  let result = Array::new()
  for i in 0..<values.length() {
    if group[i] == selected_group && period[i] == selected_period {
      result.push(values[i])
    }
  }
  result
}

///|
/// Estimates a two-period, two-group difference-in-differences effect.
pub fn difference_in_differences(
  outcomes : Array[Double],
  treated_group : Array[Bool],
  post_period : Array[Bool],
) -> DifferenceInDifferencesResult {
  let treated_pre = select_panel(
    outcomes, treated_group, post_period, true, false,
  )
  let treated_post = select_panel(
    outcomes, treated_group, post_period, true, true,
  )
  let control_pre = select_panel(
    outcomes, treated_group, post_period, false, false,
  )
  let control_post = select_panel(
    outcomes, treated_group, post_period, false, true,
  )
  let treatment_change = mean(treated_post) - mean(treated_pre)
  let control_change = mean(control_post) - mean(control_pre)
  let estimate = treatment_change - control_change
  let treated_variance = if treated_pre.length() == 0 ||
    treated_post.length() == 0 {
    0.0
  } else {
    variance(treated_pre) / treated_pre.length().to_double() +
    variance(treated_post) / treated_post.length().to_double()
  }
  let control_variance = if control_pre.length() == 0 ||
    control_post.length() == 0 {
    0.0
  } else {
    variance(control_pre) / control_pre.length().to_double() +
    variance(control_post) / control_post.length().to_double()
  }
  {
    treated_pre_mean: mean(treated_pre),
    treated_post_mean: mean(treated_post),
    control_pre_mean: mean(control_pre),
    control_post_mean: mean(control_post),
    treatment_change,
    control_change,
    estimate,
    standard_error: (treated_variance + control_variance).sqrt(),
  }
}

///|
pub struct EventStudyPoint {
  relative_time : Int
  treated_mean : Double
  control_mean : Double
  effect : Double
  count : Int
}

///|
pub fn event_study(
  outcomes : Array[Double],
  treated_group : Array[Bool],
  relative_time : Array[Int],
) -> Array[EventStudyPoint] {
  let times = Array::new()
  for time in relative_time {
    if !times.contains(time) {
      times.push(time)
    }
  }
  let result : Array[EventStudyPoint] = Array::new(capacity=times.length())
  for time in times {
    let treated = Array::new()
    let control = Array::new()
    for i in 0..<outcomes.length() {
      if relative_time[i] == time {
        if treated_group[i] {
          treated.push(outcomes[i])
        } else {
          control.push(outcomes[i])
        }
      }
    }
    result.push({
      relative_time: time,
      treated_mean: mean(treated),
      control_mean: mean(control),
      effect: mean(treated) - mean(control),
      count: treated.length() + control.length(),
    })
  }
  result
}

///|
pub fn pretrend_difference(
  outcomes : Array[Double],
  treated_group : Array[Bool],
  relative_time : Array[Int],
) -> Double {
  let pre_treated = Array::new()
  let pre_control = Array::new()
  for i in 0..<outcomes.length() {
    if relative_time[i] < 0 {
      if treated_group[i] {
        pre_treated.push(outcomes[i])
      } else {
        pre_control.push(outcomes[i])
      }
    }
  }
  mean(pre_treated) - mean(pre_control)
}

///|
pub fn synthetic_control_weights(
  treated_profile : Array[Double],
  donor_profiles : Array[Array[Double]],
) -> Array[Double] {
  let weights = Array::make(donor_profiles.length(), 1.0)
  for i in 0..<donor_profiles.length() {
    let distance = l2_norm(subtract_vectors(treated_profile, donor_profiles[i]))
    weights[i] = 1.0 / (1.0 + distance)
  }
  let total = sum(weights)
  if total > 0.0 {
    for i in 0..<weights.length() {
      weights[i] /= total
    }
  }
  weights
}

///|
pub fn weighted_donor_outcome(
  outcomes : Array[Double],
  weights : Array[Double],
) -> Double {
  weighted_mean(outcomes, weights)
}
