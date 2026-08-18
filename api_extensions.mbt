///|
pub fn CausalDataset::column(
  self : CausalDataset,
  index : Int,
) -> Array[Double] {
  column_values(self.covariates, index)
}

///|
pub fn CausalDataset::treatment_rate(self : CausalDataset) -> Double {
  if self.n() == 0 {
    0.0
  } else {
    self.treated_count().to_double() / self.n().to_double()
  }
}

///|
pub fn CausalDataset::outcome_mean(
  self : CausalDataset,
  treated : Bool,
) -> Double {
  mean(group_values(self.outcome, self.treatment, treated))
}

///|
pub fn CausalDataset::outcome_range(self : CausalDataset) -> (Double, Double) {
  (quantile(self.outcome, 0.0), quantile(self.outcome, 1.0))
}

///|
pub fn CausalDataset::with_feature_names(
  self : CausalDataset,
  names : Array[String],
) -> CausalDataset {
  {
    covariates: self.covariates,
    treatment: self.treatment,
    outcome: self.outcome,
    feature_names: names,
  }
}

///|
pub fn CausalDataset::with_outcome(
  self : CausalDataset,
  outcome : Array[Double],
) -> CausalDataset {
  {
    covariates: self.covariates,
    treatment: self.treatment,
    outcome,
    feature_names: self.feature_names,
  }
}

///|
pub fn CausalDataset::filter_by_score(
  self : CausalDataset,
  scores : Array[Double],
  lower : Double,
  upper : Double,
) -> CausalDataset {
  let indices = Array::new()
  let n = if self.n() < scores.length() { self.n() } else { scores.length() }
  for i in 0..<n {
    if scores[i] >= lower && scores[i] <= upper {
      indices.push(i)
    }
  }
  self.select(indices)
}

///|
pub fn CausalDataset::validate_alignment(self : CausalDataset) -> Bool {
  self.is_valid() && self.feature_names.length() == self.p()
}

///|
pub fn Estimate::is_statistically_significant(
  self : Estimate,
  alpha? : Double = 0.05,
) -> Bool {
  let critical = if alpha <= 0.01 {
    2.575829
  } else if alpha <= 0.05 {
    1.959964
  } else {
    1.644854
  }
  treatment_t_statistic(self.estimate, self.standard_error).abs() >= critical
}

///|
pub fn Estimate::relative_precision(self : Estimate) -> Double {
  if self.estimate == 0.0 {
    0.0
  } else {
    self.standard_error / self.estimate.abs()
  }
}

///|
pub fn Estimate::as_vector(self : Estimate) -> Array[Double] {
  [
    self.estimate,
    self.standard_error,
    self.lower,
    self.upper,
    self.effective_sample_size,
  ]
}

///|
pub fn ModelFit::coefficient_l2(self : ModelFit) -> Double {
  l2_norm(self.coefficients)
}

///|
pub fn ModelFit::has_converged(self : ModelFit) -> Bool {
  self.converged
}

///|
pub fn propensity_histogram(scores : Array[Double], bins : Int) -> Array[Int] {
  let actual = if bins < 1 { 1 } else { bins }
  let result = Array::make(actual, 0)
  for score in scores {
    let mut index = (clamp(score, 0.0, 1.0) * actual.to_double()).to_int()
    if index == actual {
      index -= 1
    }
    result[index] += 1
  }
  result
}

///|
pub fn stratum_counts(
  strata : Array[Int],
  number_of_strata : Int,
) -> Array[Int] {
  let result = Array::make(
    if number_of_strata > 0 {
      number_of_strata
    } else {
      0
    },
    0,
  )
  for stratum in strata {
    if stratum >= 0 && stratum < result.length() {
      result[stratum] += 1
    }
  }
  result
}

///|
pub fn group_mean(
  values : Array[Double],
  group : Array[Bool],
  selected : Bool,
) -> Double {
  mean(group_values(values, group, selected))
}

///|
pub fn group_standard_error(
  values : Array[Double],
  group : Array[Bool],
  selected : Bool,
) -> Double {
  let selected_values = group_values(values, group, selected)
  if selected_values.length() == 0 {
    0.0
  } else {
    std_dev(selected_values) / selected_values.length().to_double().sqrt()
  }
}

///|
pub fn common_support_indices(
  scores : Array[Double],
  lower : Double,
  upper : Double,
) -> Array[Int] {
  let result = Array::new()
  for i in 0..<scores.length() {
    if scores[i] >= lower && scores[i] <= upper {
      result.push(i)
    }
  }
  result
}

///|
pub fn safe_logit(probability : Double) -> Double {
  let p = safe_probability(probability)
  @math.ln(p / (1.0 - p))
}

///|
pub fn odds_from_probability(probability : Double) -> Double {
  let p = safe_probability(probability)
  p / (1.0 - p)
}

///|
pub fn probability_from_odds(odds : Double) -> Double {
  if odds <= 0.0 {
    0.0
  } else {
    odds / (1.0 + odds)
  }
}

///|
/// Exposes the deterministic benchmark fields to CLI and downstream users.
pub fn BenchmarkReport::sample_size(self : BenchmarkReport) -> Int {
  self.sample_size
}

///|
pub fn BenchmarkReport::true_ate(self : BenchmarkReport) -> Double {
  self.true_ate
}

///|
pub fn BenchmarkReport::estimated_ate(self : BenchmarkReport) -> Double {
  self.estimated_ate
}

///|
pub fn BenchmarkReport::absolute_error(self : BenchmarkReport) -> Double {
  self.absolute_error
}

///|
pub fn BenchmarkReport::standard_error(self : BenchmarkReport) -> Double {
  self.standard_error
}

///|
pub fn BenchmarkReport::effective_sample_size(self : BenchmarkReport) -> Double {
  self.effective_sample_size
}

///|
pub fn BenchmarkReport::propensity_auc(self : BenchmarkReport) -> Double {
  self.propensity_auc
}
