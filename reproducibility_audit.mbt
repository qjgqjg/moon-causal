///|
/// Reproducible run manifest.
pub struct RunManifest {
  project : String
  version : String
  seed : UInt64
  data_fingerprint : UInt64
  plan_fingerprint : UInt64
  source_fingerprint : UInt64
  parameters : Array[Double]
}

///|
/// Comparison between two run manifests.
pub struct ManifestComparison {
  same_data : Bool
  same_plan : Bool
  same_source : Bool
  same_seed : Bool
  parameter_distance : Double
  reproducible : Bool
}

///|
/// Immutable audit checkpoint.
pub struct AuditCheckpoint {
  name : String
  timestamp : Int
  value : Double
  fingerprint : UInt64
  passed : Bool
}

///|
/// Creates a run manifest from core fingerprints.
pub fn run_manifest(
  project : String,
  version : String,
  seed : UInt64,
  data_fingerprint : UInt64,
  plan_fingerprint : UInt64,
  source_fingerprint : UInt64,
  parameters : Array[Double],
) -> RunManifest {
  {
    project,
    version,
    seed,
    data_fingerprint,
    plan_fingerprint,
    source_fingerprint,
    parameters: parameters.copy(),
  }
}

///|
/// Combines fingerprints using a stable order-sensitive hash.
pub fn combine_fingerprints(values : Array[UInt64]) -> UInt64 {
  let mut result : UInt64 = 1469598103934665603
  for value in values {
    result = (result ^ value) * 1099511628211
  }
  result
}

///|
/// Compares two manifests for exact reproducibility.
pub fn compare_manifests(
  first : RunManifest,
  second : RunManifest,
  tolerance? : Double = 1.0e-12,
) -> ManifestComparison {
  let n = first.parameters.length().min(second.parameters.length())
  let mut distance = 0.0
  for i in 0..<n {
    distance += (first.parameters[i] - second.parameters[i]).abs()
  }
  let same_data = first.data_fingerprint == second.data_fingerprint
  let same_plan = first.plan_fingerprint == second.plan_fingerprint
  let same_source = first.source_fingerprint == second.source_fingerprint
  let same_seed = first.seed == second.seed
  {
    same_data,
    same_plan,
    same_source,
    same_seed,
    parameter_distance: distance,
    reproducible: same_data &&
    same_plan &&
    same_source &&
    same_seed &&
    distance <= tolerance,
  }
}

///|
/// Computes a fingerprint for a causal dataset including treatment and outcome.
pub fn causal_dataset_fingerprint(dataset : CausalDataset) -> UInt64 {
  let rows : Array[Array[Double]] = Array::new(capacity=dataset.n())
  for i in 0..<dataset.n() {
    let row = Array::new(capacity=dataset.p() + 2)
    row.push(if dataset.treatment[i] { 1.0 } else { 0.0 })
    row.push(dataset.outcome[i])
    for column in dataset.covariates {
      row.push(column[i])
    }
    rows.push(row)
  }
  matrix_checksum(rows)
}

///|
/// Computes a plan fingerprint from numerical configuration fields.
pub fn analysis_plan_fingerprint(plan : AnalysisPlan) -> UInt64 {
  matrix_checksum([
    [
      plan.trim_lower,
      plan.trim_upper,
      plan.bootstrap_replicates.to_double(),
      plan.cross_fit_folds.to_double(),
      plan.seed.to_double(),
      if plan.require_overlap {
        1.0
      } else {
        0.0
      },
      plan.minimum_quality_score,
    ],
  ])
}

///|
/// Creates a checkpoint for an analysis stage.
pub fn audit_checkpoint(
  name : String,
  timestamp : Int,
  value : Double,
  fingerprint : UInt64,
  passed : Bool,
) -> AuditCheckpoint {
  { name, timestamp, value, fingerprint, passed }
}

///|
/// Computes a checkpoint chain fingerprint.
pub fn checkpoint_chain(checkpoints : Array[AuditCheckpoint]) -> UInt64 {
  let mut result : UInt64 = 2166136261
  for checkpoint in checkpoints {
    result = (result ^ checkpoint.fingerprint) * 16777619
    result = (result ^ checkpoint.timestamp.to_uint64()) * 16777619
  }
  result
}

///|
/// Returns whether all checkpoints passed.
pub fn checkpoints_pass(checkpoints : Array[AuditCheckpoint]) -> Bool {
  for checkpoint in checkpoints {
    if !checkpoint.passed {
      return false
    }
  }
  checkpoints.length() > 0
}

///|
/// Computes a data-contract fingerprint from table shape and values.
pub fn table_run_fingerprint(
  table : Array[Array[Double]],
  schema : TableSchema,
) -> UInt64 {
  combine_fingerprints([
    matrix_checksum(table),
    table_schema_fingerprint(schema),
    table.length().to_uint64(),
    (if table.length() == 0 { 0 } else { table[0].length() }).to_uint64(),
  ])
}

///|
/// Builds a manifest from a high-level pipeline result.
pub fn pipeline_manifest(
  project : String,
  version : String,
  dataset : CausalDataset,
  plan : AnalysisPlan,
  result : PipelineResult,
  source_fingerprint : UInt64,
) -> RunManifest {
  run_manifest(
    project,
    version,
    plan.seed,
    causal_dataset_fingerprint(dataset),
    analysis_plan_fingerprint(plan),
    source_fingerprint,
    pipeline_summary(result),
  )
}

///|
/// Returns a compact manifest summary vector.
pub fn manifest_summary(manifest : RunManifest) -> Array[Double] {
  [
    manifest.seed.to_double(),
    manifest.data_fingerprint.to_double(),
    manifest.plan_fingerprint.to_double(),
    manifest.source_fingerprint.to_double(),
    manifest.parameters.length().to_double(),
  ]
}
