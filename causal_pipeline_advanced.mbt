///|
/// Declarative analysis plan for a reproducible causal run.
pub struct AnalysisPlan {
  estimand : String
  trim_lower : Double
  trim_upper : Double
  bootstrap_replicates : Int
  cross_fit_folds : Int
  seed : UInt64
  require_overlap : Bool
  minimum_quality_score : Double
}

///|
/// Stage record emitted by the high-level causal pipeline.
pub struct PipelineStage {
  name : String
  status : String
  value : Double
  message : String
}

///|
/// High-level analysis result with production diagnostics.
pub struct PipelineResult {
  estimate : AdvancedEffect
  propensity : Array[Double]
  weights : Array[Double]
  quality : DatasetQuality
  positivity : PositivityProfile
  stages : Array[PipelineStage]
  score : Double
  passes : Bool
}

///|
/// Creates a validated causal analysis plan.
pub fn analysis_plan(
  estimand? : String = "ATE",
  trim_lower? : Double = 0.02,
  trim_upper? : Double = 0.98,
  bootstrap_replicates? : Int = 300,
  cross_fit_folds? : Int = 5,
  seed? : UInt64 = 20260819,
  require_overlap? : Bool = true,
  minimum_quality_score? : Double = 0.8,
) -> AnalysisPlan {
  {
    estimand,
    trim_lower: clamp(trim_lower, 0.0, 0.5),
    trim_upper: clamp(trim_upper, 0.5, 1.0),
    bootstrap_replicates: bootstrap_replicates.max(0),
    cross_fit_folds: cross_fit_folds.max(2),
    seed,
    require_overlap,
    minimum_quality_score: clamp(minimum_quality_score, 0.0, 1.0),
  }
}

///|
/// Validates an analysis plan's ordering and operational constraints.
pub fn validate_analysis_plan(plan : AnalysisPlan) -> Bool {
  plan.trim_lower < plan.trim_upper &&
  plan.cross_fit_folds >= 2 &&
  plan.minimum_quality_score >= 0.0 &&
  plan.minimum_quality_score <= 1.0 &&
  plan.estimand != ""
}

///|
fn pipeline_rows_from_columns(
  columns : Array[Array[Double]],
  rows : Int,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=rows)
  for i in 0..<rows {
    let row = Array::new(capacity=columns.length())
    for column in columns {
      if i < column.length() {
        row.push(column[i])
      } else {
        row.push(0.0)
      }
    }
    result.push(row)
  }
  result
}

///|
fn pipeline_stage(
  name : String,
  status : String,
  value : Double,
  message : String,
) -> PipelineStage {
  { name, status, value, message }
}

///|
/// Runs a complete observational estimate from a causal dataset.
pub fn run_causal_pipeline(
  dataset : CausalDataset,
  plan : AnalysisPlan,
) -> PipelineResult {
  let empty_effect : AdvancedEffect = {
    estimate: 0.0,
    standard_error: 0.0,
    lower: 0.0,
    upper: 0.0,
    effective_sample_size: 0.0,
    estimand: plan.estimand,
    passes: false,
  }
  let empty_quality : DatasetQuality = assess_matrix([])
  if !validate_analysis_plan(plan) || !dataset.is_valid() {
    return {
      estimate: empty_effect,
      propensity: [],
      weights: [],
      quality: empty_quality,
      positivity: {
        treated: 0,
        control: 0,
        minimum_score: 0.0,
        maximum_score: 0.0,
        near_zero: 0,
        near_one: 0,
        effective_sample_size: 0.0,
        passes: false,
      },
      stages: [
        pipeline_stage("validation", "failed", 0.0, "invalid plan or dataset"),
      ],
      score: 0.0,
      passes: false,
    }
  }
  let rows = pipeline_rows_from_columns(dataset.covariates, dataset.n())
  let quality = assess_matrix(rows)
  let stages : Array[PipelineStage] = Array::new()
  stages.push(
    pipeline_stage(
      "validation",
      "passed",
      quality.score,
      "dataset quality assessed",
    ),
  )
  let propensity_model = fit_logistic_regression(
    rows,
    dataset.treatment,
    max_iterations=500,
  )
  let propensity = predict_propensity(propensity_model, rows)
  let trimmed = trim_propensity_scores(
    propensity,
    lower=plan.trim_lower,
    upper=plan.trim_upper,
  )
  let positivity = positivity_profile(dataset.treatment, trimmed)
  stages.push(
    pipeline_stage(
      "propensity",
      if propensity_model.converged {
        "passed"
      } else {
        "warning"
      },
      propensity_model.loss,
      "propensity model fitted",
    ),
  )
  stages.push(
    pipeline_stage(
      "overlap",
      if positivity.passes {
        "passed"
      } else {
        "failed"
      },
      positivity.effective_sample_size,
      "positivity diagnostics",
    ),
  )
  let treated_rows : Array[Array[Double]] = Array::new()
  let treated_outcomes : Array[Double] = Array::new()
  let control_rows : Array[Array[Double]] = Array::new()
  let control_outcomes : Array[Double] = Array::new()
  for i in 0..<dataset.n() {
    if dataset.treatment[i] {
      treated_rows.push(rows[i])
      treated_outcomes.push(dataset.outcome[i])
    } else {
      control_rows.push(rows[i])
      control_outcomes.push(dataset.outcome[i])
    }
  }
  let treated_model = fit_linear_outcome_model(treated_rows, treated_outcomes)
  let control_model = fit_linear_outcome_model(control_rows, control_outcomes)
  let predicted_treated = predict_outcomes(treated_model, rows)
  let predicted_control = predict_outcomes(control_model, rows)
  let effect = augmented_inverse_probability_effect(
    dataset.outcome,
    dataset.treatment,
    trimmed,
    predicted_treated,
    predicted_control,
  )
  let advanced = {
    estimate: effect.estimate,
    standard_error: effect.standard_error,
    lower: effect.estimate - 1.96 * effect.standard_error,
    upper: effect.estimate + 1.96 * effect.standard_error,
    effective_sample_size: positivity.effective_sample_size,
    estimand: plan.estimand,
    passes: effect.passes &&
    quality_gate(quality, minimum_score=plan.minimum_quality_score) &&
    (!plan.require_overlap || positivity.passes),
  }
  stages.push(
    pipeline_stage(
      "outcome",
      if treated_model.converged && control_model.converged {
        "passed"
      } else {
        "warning"
      },
      treated_model.loss + control_model.loss,
      "potential-outcome models fitted",
    ),
  )
  stages.push(
    pipeline_stage(
      "estimate",
      if advanced.passes {
        "passed"
      } else {
        "failed"
      },
      advanced.estimate,
      "augmented estimate computed",
    ),
  )
  let weights = advanced_stabilized_weights(dataset.treatment, trimmed)
  let overlap_score = if positivity.passes { 1.0 } else { 0.0 }
  let effect_score = if effect.passes { 1.0 } else { 0.0 }
  let score = clamp(
    0.4 * quality.score + 0.3 * overlap_score + 0.3 * effect_score,
    0.0,
    1.0,
  )
  {
    estimate: advanced,
    propensity: trimmed,
    weights,
    quality,
    positivity,
    stages,
    score,
    passes: advanced.passes && score >= plan.minimum_quality_score,
  }
}

///|
/// Builds a stage audit vector for monitoring.
pub fn pipeline_stage_vector(stages : Array[PipelineStage]) -> Array[Double] {
  let result = Array::new(capacity=stages.length())
  for stage in stages {
    result.push(stage.value)
  }
  result
}

///|
/// Serializes stage statuses for human-readable run logs.
pub fn pipeline_stage_text(stages : Array[PipelineStage]) -> String {
  let builder = StringBuilder::new()
  for stage in stages {
    builder.write_string(stage.name)
    builder.write_string("=")
    builder.write_string(stage.status)
    builder.write_string("|")
    builder.write_string(stage.message)
    builder.write_string("\n")
  }
  builder.to_string()
}

///|
/// Computes a deterministic pipeline fingerprint.
pub fn pipeline_fingerprint(result : PipelineResult) -> UInt64 {
  let rows : Array[Array[Double]] = [
    advanced_effect_summary(result.estimate),
    pipeline_stage_vector(result.stages),
  ]
  matrix_checksum(rows)
}

///|
/// Returns a compact pipeline summary vector.
pub fn pipeline_summary(result : PipelineResult) -> Array[Double] {
  [
    result.estimate.estimate,
    result.estimate.standard_error,
    result.estimate.lower,
    result.estimate.upper,
    result.quality.score,
    result.positivity.effective_sample_size,
    result.score,
    if result.passes {
      1.0
    } else {
      0.0
    },
  ]
}
