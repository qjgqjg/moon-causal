///|
pub struct ObservationalAnalysis {
  dataset : CausalDataset
  propensity_model : ModelFit
  propensity_scores : Array[Double]
  ipw_estimate : Estimate
  balance : Array[BalanceMetric]
  overlap : OverlapReport
}

///|
pub fn run_observational_analysis(
  dataset : CausalDataset,
) -> ObservationalAnalysis {
  let model = fit_logistic_regression(dataset.covariates, dataset.treatment)
  let propensity = predict_propensity(model, dataset.covariates)
  let estimate = estimate_ipw_ate(
    dataset.outcome,
    dataset.treatment,
    propensity,
  )
  {
    dataset,
    propensity_model: model,
    propensity_scores: propensity,
    ipw_estimate: estimate,
    balance: balance_table(
      dataset.covariates,
      dataset.treatment,
      dataset.feature_names,
    ),
    overlap: overlap_report(propensity, dataset.treatment),
  }
}

///|
pub fn analysis_quality_score(analysis : ObservationalAnalysis) -> Double {
  let overlap_score = analysis.overlap.overlap_fraction
  let balance_score = 1.0 -
    clamp(maximum_absolute_smd(analysis.balance), 0.0, 1.0)
  let convergence_score = if analysis.propensity_model.converged {
    1.0
  } else {
    0.5
  }
  (overlap_score + balance_score + convergence_score) / 3.0
}

///|
pub fn doubly_robust_from_dataset(analysis : ObservationalAnalysis) -> Estimate {
  let treated_model = fit_linear_outcome_model(
    analysis.dataset.covariates,
    analysis.dataset.outcome,
  )
  let control_model = fit_linear_outcome_model(
    analysis.dataset.covariates,
    analysis.dataset.outcome,
  )
  let predicted_treated = predict_outcomes(
    treated_model,
    analysis.dataset.covariates,
  )
  let predicted_control = predict_outcomes(
    control_model,
    analysis.dataset.covariates,
  )
  let estimate = doubly_robust_ate(
    analysis.dataset.outcome,
    analysis.dataset.treatment,
    analysis.propensity_scores,
    predicted_treated,
    predicted_control,
  )
  Estimate::from_standard_error(
    estimate,
    analysis.ipw_estimate.standard_error,
    analysis.dataset.n(),
    analysis.ipw_estimate.effective_sample_size,
    "ATE (doubly robust)",
  )
}

///|
pub fn analysis_summary_vector(
  analysis : ObservationalAnalysis,
) -> Array[Double] {
  [
    analysis.ipw_estimate.estimate,
    analysis.ipw_estimate.standard_error,
    analysis.overlap.overlap_fraction,
    maximum_absolute_smd(analysis.balance),
    analysis.quality_score(),
  ]
}

///|
pub fn ObservationalAnalysis::quality_score(
  self : ObservationalAnalysis,
) -> Double {
  analysis_quality_score(self)
}
