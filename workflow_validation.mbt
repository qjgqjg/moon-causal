///|
/// A single preflight check emitted before a causal analysis is released.
pub struct WorkflowCheck {
  name : String
  passed : Bool
  score : Double
  severity : Int
  detail : String
}

///|
/// Aggregated validation result for a data-to-report workflow.
pub struct WorkflowReport {
  checks : Array[WorkflowCheck]
  passed : Bool
  score : Double
  critical_failures : Int
  warnings : Int
}

///|
/// Creates a normalized workflow check.
pub fn workflow_check(
  name : String,
  passed : Bool,
  score : Double,
  severity? : Int = 1,
  detail? : String = "",
) -> WorkflowCheck {
  {
    name,
    passed,
    score: clamp(score, 0.0, 1.0),
    severity: severity.max(1).min(3),
    detail,
  }
}

///|
/// Aggregates checks using severity-weighted scores.
pub fn workflow_report(checks : Array[WorkflowCheck]) -> WorkflowReport {
  let mut weighted = 0.0
  let mut total_weight = 0.0
  let mut critical = 0
  let mut warnings = 0
  for check in checks {
    let weight = check.severity.to_double()
    weighted += weight * check.score
    total_weight += weight
    if !check.passed && check.severity >= 3 {
      critical += 1
    }
    if !check.passed && check.severity == 2 {
      warnings += 1
    }
  }
  let score = if total_weight == 0.0 { 0.0 } else { weighted / total_weight }
  {
    checks: checks.copy(),
    passed: checks.length() > 0 &&
    critical == 0 &&
    checks.fold(init=true, fn(ok, check) { ok && check.passed }),
    score,
    critical_failures: critical,
    warnings,
  }
}

///|
/// Adds one check without mutating an existing report.
pub fn workflow_report_with(
  report : WorkflowReport,
  check : WorkflowCheck,
) -> WorkflowReport {
  let checks = report.checks.copy()
  checks.push(check)
  workflow_report(checks)
}

///|
/// Returns the number of checks that passed.
pub fn workflow_pass_count(report : WorkflowReport) -> Int {
  report.checks.fold(init=0, fn(total, check) {
    if check.passed {
      total + 1
    } else {
      total
    }
  })
}

///|
/// Returns the pass rate for a workflow report.
pub fn workflow_pass_rate(report : WorkflowReport) -> Double {
  if report.checks.length() == 0 {
    0.0
  } else {
    workflow_pass_count(report).to_double() / report.checks.length().to_double()
  }
}

///|
/// Returns failed check names in declaration order.
pub fn workflow_failed_names(report : WorkflowReport) -> Array[String] {
  let result : Array[String] = Array::new()
  for check in report.checks {
    if !check.passed {
      result.push(check.name)
    }
  }
  result
}

///|
/// Returns the first failed check, if any.
pub fn workflow_first_failure(report : WorkflowReport) -> WorkflowCheck? {
  for check in report.checks {
    if !check.passed {
      return Some(check)
    }
  }
  None
}

///|
/// Counts failed checks at or above a severity threshold.
pub fn workflow_failure_count(
  report : WorkflowReport,
  minimum_severity : Int,
) -> Int {
  report.checks.fold(init=0, fn(total, check) {
    if !check.passed && check.severity >= minimum_severity {
      total + 1
    } else {
      total
    }
  })
}

///|
/// Extracts check scores for dashboards and monitoring.
pub fn workflow_score_vector(report : WorkflowReport) -> Array[Double] {
  let result : Array[Double] = Array::new(capacity=report.checks.length())
  for check in report.checks {
    result.push(check.score)
  }
  result
}

///|
/// Serializes a workflow report as stable line-oriented text.
pub fn workflow_report_text(report : WorkflowReport) -> String {
  let builder = StringBuilder::new()
  builder.write_string("workflow_score=")
  builder.write_string(report.score.to_string())
  builder.write_string("\nworkflow_passed=")
  builder.write_string(report.passed.to_string())
  builder.write_string("\n")
  for check in report.checks {
    builder.write_string(check.name)
    builder.write_string("=")
    builder.write_string(if check.passed { "passed" } else { "failed" })
    builder.write_string(";score=")
    builder.write_string(check.score.to_string())
    builder.write_string(";severity=")
    builder.write_string(check.severity.to_string())
    if check.detail != "" {
      builder.write_string(";")
      builder.write_string(check.detail)
    }
    builder.write_string("\n")
  }
  builder.to_string()
}

///|
/// Checks that the dataset has aligned rectangular dimensions.
pub fn check_dataset_shape(dataset : CausalDataset) -> WorkflowCheck {
  let valid = dataset.is_valid()
  workflow_check(
    "dataset-shape",
    valid,
    if valid {
      1.0
    } else {
      0.0
    },
    severity=3,
    detail=if valid { "rectangular and aligned" } else { "invalid dimensions" },
  )
}

///|
/// Checks that both treatment arms are represented.
pub fn check_dataset_support(dataset : CausalDataset) -> WorkflowCheck {
  let treated = dataset.treated_count()
  let control = dataset.control_count()
  let total = dataset.n()
  let passed = treated > 0 && control > 0
  let score = if total == 0 {
    0.0
  } else {
    treated.min(control).to_double() / (total.to_double() / 2.0)
  }
  workflow_check(
    "treatment-support",
    passed,
    clamp(score, 0.0, 1.0),
    severity=3,
    detail="treated=\{treated};control=\{control}",
  )
}

///|
/// Checks that outcomes and covariates contain only finite values.
pub fn check_dataset_finite(dataset : CausalDataset) -> WorkflowCheck {
  let mut finite = true
  for value in dataset.outcome {
    if !is_finite(value) {
      finite = false
      break
    }
  }
  if finite {
    for row in dataset.covariates {
      for value in row {
        if !is_finite(value) {
          finite = false
          break
        }
      }
      if !finite {
        break
      }
    }
  }
  workflow_check(
    "finite-values",
    finite,
    if finite {
      1.0
    } else {
      0.0
    },
    severity=3,
    detail=if finite {
      "all numeric values finite"
    } else {
      "non-finite value found"
    },
  )
}

///|
/// Checks minimum observations per treatment arm.
pub fn check_dataset_minimum_arm(
  dataset : CausalDataset,
  minimum_per_arm : Int,
) -> WorkflowCheck {
  let minimum = minimum_per_arm.max(1)
  let observed = dataset.treated_count().min(dataset.control_count())
  let score = clamp(observed.to_double() / minimum.to_double(), 0.0, 1.0)
  workflow_check(
    "minimum-arm-size",
    observed >= minimum,
    score,
    severity=2,
    detail="minimum=\{minimum};observed=\{observed}",
  )
}

///|
/// Checks variation in each covariate column.
pub fn check_dataset_variation(dataset : CausalDataset) -> WorkflowCheck {
  let profiles = profile_matrix(dataset.covariates)
  let mut varying = 0
  for profile in profiles {
    if profile.std_dev > 0.0 && profile.count > 1 {
      varying += 1
    }
  }
  let columns = dataset.p()
  let score = if columns == 0 {
    0.0
  } else {
    varying.to_double() / columns.to_double()
  }
  workflow_check(
    "covariate-variation",
    columns == 0 || varying == columns,
    score,
    severity=2,
    detail="varying=\{varying};columns=\{columns}",
  )
}

///|
/// Checks a plan's bounds and operational settings.
pub fn check_plan_configuration(plan : AnalysisPlan) -> WorkflowCheck {
  let passed = validate_analysis_plan(plan) && plan.bootstrap_replicates >= 0
  let score = if !passed {
    0.0
  } else if plan.bootstrap_replicates == 0 {
    0.8
  } else {
    1.0
  }
  workflow_check(
    "analysis-plan",
    passed,
    score,
    severity=3,
    detail="estimand=\{plan.estimand};folds=\{plan.cross_fit_folds}",
  )
}

///|
/// Checks that a graph is acyclic and has a useful ordering.
pub fn check_graph_acyclic(graph : CausalGraph) -> WorkflowCheck {
  let order = graph.topological_order()
  let passed = order.length() > 0 || graph.nodes.length() == 0
  let score = if graph.nodes.length() == 0 {
    0.0
  } else {
    order.length().to_double() / graph.nodes.length().to_double()
  }
  workflow_check(
    "causal-graph",
    passed,
    clamp(score, 0.0, 1.0),
    severity=3,
    detail="nodes=\{graph.nodes.length()};edges=\{graph.edge_count()}",
  )
}

///|
/// Checks that treatment and outcome are declared in a graph.
pub fn check_graph_endpoints(
  graph : CausalGraph,
  treatment : String,
  outcome : String,
) -> WorkflowCheck {
  let found_treatment = graph.nodes.contains(treatment)
  let found_outcome = graph.nodes.contains(outcome)
  let passed = found_treatment && found_outcome && treatment != outcome
  workflow_check(
    "graph-endpoints",
    passed,
    if passed {
      1.0
    } else {
      0.0
    },
    severity=3,
    detail="treatment=\{found_treatment};outcome=\{found_outcome}",
  )
}

///|
/// Checks positivity diagnostics against an operational threshold.
pub fn check_positivity(
  profile : PositivityProfile,
  minimum_effective_sample_size : Double,
) -> WorkflowCheck {
  let threshold = minimum_effective_sample_size.max(1.0)
  let passed = profile.passes && profile.effective_sample_size >= threshold
  let score = clamp(profile.effective_sample_size / threshold, 0.0, 1.0)
  workflow_check(
    "positivity",
    passed,
    score,
    severity=3,
    detail="ess=\{profile.effective_sample_size};threshold=\{threshold}",
  )
}

///|
/// Checks the dataset quality score against a release threshold.
pub fn check_quality(
  quality : DatasetQuality,
  minimum_score : Double,
) -> WorkflowCheck {
  let threshold = clamp(minimum_score, 0.0, 1.0)
  let score = if threshold == 0.0 { 1.0 } else { quality.score / threshold }
  workflow_check(
    "dataset-quality",
    quality.score >= threshold && quality_gate(quality, minimum_score=threshold),
    clamp(score, 0.0, 1.0),
    severity=3,
    detail="score=\{quality.score};threshold=\{threshold}",
  )
}

///|
/// Checks an effect estimate for finite uncertainty and a positive sample size.
pub fn check_effect(
  effect : AdvancedEffect,
  minimum_effective_sample_size : Double,
) -> WorkflowCheck {
  let finite = is_finite(effect.estimate) &&
    is_finite(effect.standard_error) &&
    is_finite(effect.lower) &&
    is_finite(effect.upper)
  let threshold = minimum_effective_sample_size.max(1.0)
  let passed = finite &&
    effect.passes &&
    effect.effective_sample_size >= threshold
  let score = if !finite {
    0.0
  } else {
    clamp(effect.effective_sample_size / threshold, 0.0, 1.0)
  }
  workflow_check(
    "effect-estimate",
    passed,
    score,
    severity=3,
    detail="estimand=\{effect.estimand};estimate=\{effect.estimate}",
  )
}

///|
/// Checks all stages and the release score of a pipeline result.
pub fn check_pipeline_result(
  result : PipelineResult,
  minimum_score? : Double = 0.8,
) -> WorkflowCheck {
  let threshold = clamp(minimum_score, 0.0, 1.0)
  let mut stage_failures = 0
  for stage in result.stages {
    if stage.status == "failed" {
      stage_failures += 1
    }
  }
  let passed = result.passes && result.score >= threshold && stage_failures == 0
  let score = if threshold == 0.0 { 1.0 } else { result.score / threshold }
  workflow_check(
    "causal-pipeline",
    passed,
    clamp(score, 0.0, 1.0),
    severity=3,
    detail="stage_failures=\{stage_failures};score=\{result.score}",
  )
}

///|
/// Checks a named estimand registry.
pub fn check_estimand_registry(registry : EstimandRegistry) -> WorkflowCheck {
  let audit = audit_estimand_registry(registry)
  workflow_check(
    "estimand-registry",
    audit.passes,
    if audit.passes {
      1.0
    } else {
      0.0
    },
    severity=2,
    detail="registered=\{audit.registered};missing=\{audit.missing_fields}",
  )
}

///|
/// Checks report sections for non-empty titles and stable fingerprints.
pub fn check_report_sections(sections : Array[ReportSection]) -> WorkflowCheck {
  let mut valid = sections.length() > 0
  let mut nonempty = 0
  for section in sections {
    if section.title == "" || section.fingerprint == 0UL {
      valid = false
    }
    if section.lines.length() > 0 {
      nonempty += 1
    }
  }
  let score = if sections.length() == 0 {
    0.0
  } else {
    nonempty.to_double() / sections.length().to_double()
  }
  workflow_check(
    "report-sections",
    valid,
    score,
    severity=2,
    detail="sections=\{sections.length()};nonempty=\{nonempty}",
  )
}

///|
/// Runs the dataset, plan, and graph checks used at analysis ingress.
pub fn workflow_preflight(
  dataset : CausalDataset,
  plan : AnalysisPlan,
  graph : CausalGraph,
  treatment : String,
  outcome : String,
  minimum_arm? : Int = 10,
) -> WorkflowReport {
  workflow_report([
    check_dataset_shape(dataset),
    check_dataset_finite(dataset),
    check_dataset_support(dataset),
    check_dataset_minimum_arm(dataset, minimum_arm),
    check_dataset_variation(dataset),
    check_plan_configuration(plan),
    check_graph_acyclic(graph),
    check_graph_endpoints(graph, treatment, outcome),
  ])
}

///|
/// Audits a completed pipeline together with its registry and report.
pub fn workflow_release_audit(
  result : PipelineResult,
  registry : EstimandRegistry,
  sections : Array[ReportSection],
  minimum_score? : Double = 0.8,
  minimum_effective_sample_size? : Double = 10.0,
) -> WorkflowReport {
  workflow_report([
    check_quality(result.quality, minimum_score),
    check_positivity(result.positivity, minimum_effective_sample_size),
    check_effect(result.estimate, minimum_effective_sample_size),
    check_pipeline_result(result, minimum_score~),
    check_estimand_registry(registry),
    check_report_sections(sections),
  ])
}

///|
/// Returns a compact vector for quality gates and dashboards.
pub fn workflow_summary(report : WorkflowReport) -> Array[Double] {
  [
    report.checks.length().to_double(),
    workflow_pass_count(report).to_double(),
    workflow_pass_rate(report),
    report.score,
    report.critical_failures.to_double(),
    report.warnings.to_double(),
    if report.passed {
      1.0
    } else {
      0.0
    },
  ]
}

///|
/// Computes a stable fingerprint for a workflow report.
pub fn workflow_fingerprint(report : WorkflowReport) -> UInt64 {
  let rows : Array[Array[Double]] = Array::new(capacity=report.checks.length())
  for check in report.checks {
    rows.push([
      check.score,
      check.severity.to_double(),
      if check.passed {
        1.0
      } else {
        0.0
      },
    ])
  }
  matrix_checksum(rows)
}

///|
/// Returns checks that need human review, including warnings.
pub fn workflow_review_queue(report : WorkflowReport) -> Array[String] {
  let result : Array[String] = Array::new()
  for check in report.checks {
    if !check.passed || (check.severity <= 2 && check.score < 1.0) {
      result.push(check.name)
    }
  }
  result
}

///|
/// Returns whether the report can be promoted to a release artifact.
pub fn workflow_ready_for_release(
  report : WorkflowReport,
  minimum_score? : Double = 0.9,
) -> Bool {
  report.passed &&
  report.critical_failures == 0 &&
  report.score >= clamp(minimum_score, 0.0, 1.0)
}

///|
/// Compares two workflow reports by check name and score.
pub fn workflow_score_delta(
  baseline : WorkflowReport,
  current : WorkflowReport,
) -> Array[Double] {
  let result : Array[Double] = Array::new(capacity=current.checks.length())
  for current_check in current.checks {
    let mut baseline_score = 0.0
    for baseline_check in baseline.checks {
      if baseline_check.name == current_check.name {
        baseline_score = baseline_check.score
        break
      }
    }
    result.push(current_check.score - baseline_score)
  }
  result
}

///|
/// Checks for a material degradation in a monitored workflow.
pub fn workflow_degraded(
  baseline : WorkflowReport,
  current : WorkflowReport,
  tolerance? : Double = 0.05,
) -> Bool {
  let threshold = tolerance.max(0.0)
  current.score + threshold < baseline.score ||
  current.critical_failures > baseline.critical_failures
}
