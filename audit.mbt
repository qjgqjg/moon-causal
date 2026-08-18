///|
pub struct AuditReport {
  valid_dataset : Bool
  adequate_sample : Bool
  has_treatment_variation : Bool
  has_outcome_variation : Bool
  propensity_overlap : Bool
  no_positivity_violation : Bool
  quality_score : Double
}

///|
fn audit_dataset_impl(
  dataset : CausalDataset,
  overlap : Bool,
  positivity : Bool,
) -> AuditReport {
  let valid = dataset.is_valid()
  let sample = dataset.n() >= 2
  let treatment_variation = dataset.treated_count() > 0 &&
    dataset.control_count() > 0
  let outcome_variation = has_variation(dataset.outcome)
  let score = if valid && sample && treatment_variation && outcome_variation {
    1.0
  } else {
    0.0
  }
  {
    valid_dataset: valid,
    adequate_sample: sample,
    has_treatment_variation: treatment_variation,
    has_outcome_variation: outcome_variation,
    propensity_overlap: overlap,
    no_positivity_violation: positivity,
    quality_score: score,
  }
}

///|
pub fn audit_dataset(dataset : CausalDataset) -> AuditReport {
  audit_dataset_impl(dataset, false, false)
}

///|
pub fn audit_estimate(estimate : Estimate) -> Array[Bool] {
  estimate_quality_flags(estimate)
}

///|
pub fn audit_model(model : ModelFit) -> Array[Bool] {
  model_quality_flags(model)
}

///|
pub fn audit_analysis(analysis : ObservationalAnalysis) -> AuditReport {
  audit_dataset_impl(
    analysis.dataset,
    analysis.overlap.overlap_fraction > 0.0,
    positivity_violations(analysis.propensity_scores, 0.01) == 0,
  )
}

///|
pub fn audit_passes(report : AuditReport) -> Bool {
  report.valid_dataset &&
  report.adequate_sample &&
  report.has_treatment_variation &&
  report.has_outcome_variation &&
  report.quality_score > 0.0
}
