///|
/// Effect and prediction diagnostics for one subgroup.
pub struct FairSubgroupEffect {
  group : Int
  count : Int
  treated_count : Int
  control_count : Int
  effect : Double
  standard_error : Double
  positive_rate : Double
}

///|
/// Fairness report across subgroups.
pub struct FairnessReport {
  subgroups : Array[FairSubgroupEffect]
  maximum_effect_gap : Double
  maximum_positive_rate_gap : Double
  equal_opportunity_gap : Double
  passes : Bool
}

///|
/// Computes subgroup treatment effects and positive prediction rates.
pub fn subgroup_effects(
  groups : Array[Int],
  treatment : Array[Bool],
  outcomes : Array[Double],
  predictions : Array[Bool],
) -> Array[FairSubgroupEffect] {
  let n = groups
    .length()
    .min(treatment.length())
    .min(outcomes.length())
    .min(predictions.length())
  let unique : Array[Int] = Array::new()
  for group in groups[:n] {
    if !unique.contains(group) {
      unique.push(group)
    }
  }
  let result : Array[FairSubgroupEffect] = Array::new(capacity=unique.length())
  for group in unique {
    let treated = Array::new()
    let control = Array::new()
    let mut positive = 0
    for i in 0..<n {
      if groups[i] == group {
        if treatment[i] {
          treated.push(outcomes[i])
        } else {
          control.push(outcomes[i])
        }
        if predictions[i] {
          positive += 1
        }
      }
    }
    let effect = mean_or(treated, 0.0) - mean_or(control, 0.0)
    let standard_error = if treated.length() == 0 || control.length() == 0 {
      0.0
    } else {
      (variance(treated) / treated.length().to_double() +
      variance(control) / control.length().to_double()).sqrt()
    }
    let count = treated.length() + control.length()
    result.push({
      group,
      count,
      treated_count: treated.length(),
      control_count: control.length(),
      effect,
      standard_error,
      positive_rate: if count == 0 {
        0.0
      } else {
        positive.to_double() / count.to_double()
      },
    })
  }
  result
}

///|
/// Computes true-positive rates by subgroup.
pub fn fair_true_positive_rates(
  groups : Array[Int],
  actual : Array[Bool],
  predicted : Array[Bool],
) -> Array[Array[Double]] {
  let n = groups.length().min(actual.length()).min(predicted.length())
  let unique : Array[Int] = Array::new()
  for group in groups[:n] {
    if !unique.contains(group) {
      unique.push(group)
    }
  }
  let result : Array[Array[Double]] = Array::new(capacity=unique.length())
  for group in unique {
    let mut positives = 0
    let mut true_positives = 0
    for i in 0..<n {
      if groups[i] == group && actual[i] {
        positives += 1
        if predicted[i] {
          true_positives += 1
        }
      }
    }
    result.push([
      group.to_double(),
      if positives == 0 {
        0.0
      } else {
        true_positives.to_double() / positives.to_double()
      },
      positives.to_double(),
    ])
  }
  result
}

///|
/// Computes false-positive rates by subgroup.
pub fn fair_false_positive_rates(
  groups : Array[Int],
  actual : Array[Bool],
  predicted : Array[Bool],
) -> Array[Array[Double]] {
  let n = groups.length().min(actual.length()).min(predicted.length())
  let unique : Array[Int] = Array::new()
  for group in groups[:n] {
    if !unique.contains(group) {
      unique.push(group)
    }
  }
  let result : Array[Array[Double]] = Array::new(capacity=unique.length())
  for group in unique {
    let mut negatives = 0
    let mut false_positives = 0
    for i in 0..<n {
      if groups[i] == group && !actual[i] {
        negatives += 1
        if predicted[i] {
          false_positives += 1
        }
      }
    }
    result.push([
      group.to_double(),
      if negatives == 0 {
        0.0
      } else {
        false_positives.to_double() / negatives.to_double()
      },
      negatives.to_double(),
    ])
  }
  result
}

///|
/// Produces a fairness report against maximum allowed gaps.
pub fn fairness_report(
  groups : Array[Int],
  treatment : Array[Bool],
  outcomes : Array[Double],
  predictions : Array[Bool],
  maximum_effect_gap? : Double = 0.2,
  maximum_rate_gap? : Double = 0.1,
) -> FairnessReport {
  let effects = subgroup_effects(groups, treatment, outcomes, predictions)
  let mut maximum_effect = 0.0
  let mut maximum_rate = 0.0
  for first in effects {
    for second in effects {
      let effect_gap = (first.effect - second.effect).abs()
      let rate_gap = (first.positive_rate - second.positive_rate).abs()
      if effect_gap > maximum_effect {
        maximum_effect = effect_gap
      }
      if rate_gap > maximum_rate {
        maximum_rate = rate_gap
      }
    }
  }
  let tpr = fair_true_positive_rates(groups, treatment, predictions)
  let mut opportunity_gap = 0.0
  for first in tpr {
    for second in tpr {
      let gap = (first[1] - second[1]).abs()
      if gap > opportunity_gap {
        opportunity_gap = gap
      }
    }
  }
  {
    subgroups: effects,
    maximum_effect_gap: maximum_effect,
    maximum_positive_rate_gap: maximum_rate,
    equal_opportunity_gap: opportunity_gap,
    passes: maximum_effect <= maximum_effect_gap &&
    maximum_rate <= maximum_rate_gap,
  }
}

///|
/// Computes demographic parity ratio from subgroup positive rates.
pub fn demographic_parity_ratio(report : FairnessReport) -> Double {
  if report.subgroups.length() == 0 {
    return 1.0
  }
  let mut minimum = report.subgroups[0].positive_rate
  let mut maximum = minimum
  for subgroup in report.subgroups {
    if subgroup.positive_rate < minimum {
      minimum = subgroup.positive_rate
    }
    if subgroup.positive_rate > maximum {
      maximum = subgroup.positive_rate
    }
  }
  if maximum == 0.0 {
    1.0
  } else {
    minimum / maximum
  }
}

///|
/// Returns a compact fairness summary vector.
pub fn fairness_summary(report : FairnessReport) -> Array[Double] {
  [
    report.subgroups.length().to_double(),
    report.maximum_effect_gap,
    report.maximum_positive_rate_gap,
    report.equal_opportunity_gap,
    demographic_parity_ratio(report),
    if report.passes {
      1.0
    } else {
      0.0
    },
  ]
}
