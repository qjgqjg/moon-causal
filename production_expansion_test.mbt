///|
test "production transform boundaries" {
  let empty = causal_standardize([])
  assert_eq(empty.values, [])
  assert_eq(empty.missing, 0)
  assert_true(empty.passes)
  let clipped = causal_clip([-10.0, 0.0, 10.0], -1.0, 1.0)
  assert_eq(clipped.values, [-1.0, 0.0, 1.0])
  assert_eq(clipped.clipped, 2)
  let logits = causal_logit_transform([0.0, 0.5, 1.0], epsilon=0.01)
  assert_eq(logits.values.length(), 3)
  for value in logits.values {
    assert_true(is_finite(value))
  }
  assert_eq(causal_transpose([[], []]), [])
  assert_eq(causal_from_columns([], rows=2), [[], []])
  assert_eq(causal_matrix_missing([[1.0, 2.0], [3.0]]), 0)
}

///|
test "production data view split boundaries" {
  let data = CausalDataset::new(
    [[0.0], [1.0], [2.0], [3.0], [4.0], [5.0]],
    [false, true, false, true, false, true],
    [0.0, 2.0, 1.0, 3.0, 2.0, 4.0],
  )
  let empty_split = causal_train_validation_split(0)
  assert_eq(empty_split.train.length(), 0)
  assert_eq(empty_split.validation.length(), 0)
  let split = causal_stratified_split(
    data.treatment,
    train_fraction=0.5,
    seed=7UL,
  )
  assert_eq(split.train.length() + split.validation.length(), data.n())
  assert_true(split.passes)
  let assignment = causal_fold_assignment(data.n(), 3, seed=7UL)
  assert_eq(assignment.folds.length(), data.n())
  assert_eq(assignment.sizes.length(), 3)
  for fold in causal_all_folds(data.n(), 3, seed=7UL) {
    assert_true(fold.passes)
  }
  let bootstrap = causal_bootstrap_indices(data.n(), 4, seed=7UL)
  assert_eq(bootstrap.length(), 4)
  for indices in bootstrap {
    assert_eq(indices.length(), data.n())
  }
  assert_eq(causal_complete_rows(data).n(), data.n())
  assert_eq(causal_data_view_summary(data).rows, 6)
}

///|
test "production strata and design matrix" {
  let data = CausalDataset::new(
    [[0.0, 1.0], [1.0, 2.0], [2.0, 3.0], [3.0, 4.0]],
    [false, true, false, true],
    [1.0, 3.0, 2.0, 4.0],
  )
  let design = causal_design_matrix(data, include_intercept=true, interaction_pairs=[
    (1, 2),
  ])
  assert_eq(design.length(), 4)
  assert_eq(design[0].length(), 4)
  let strata = causal_score_strata([0.0, 0.25, 0.75, 1.0], 2)
  assert_eq(strata, [0, 0, 1, 1])
  let summaries = causal_stratum_summaries(
    strata,
    data.treatment,
    data.outcome,
    2,
  )
  assert_eq(summaries.length(), 2)
  assert_eq(summaries[0].count + summaries[1].count, 4)
  assert_eq(causal_select_columns(data.covariates, [1]).length(), 4)
  assert_eq(causal_add_intercept(data.covariates)[0].length(), 3)
}

///|
test "production inference boundaries" {
  assert_close(inference_normal_cdf(0.0), 0.5, 1.0e-4)
  assert_true(inference_normal_cdf(-2.0) < inference_normal_cdf(2.0))
  let interval = inference_interval(2.0, 0.5)
  assert_true(interval.lower < interval.estimate)
  assert_true(interval.upper > interval.estimate)
  let contrast = inference_contrast("ate", 2.0, 0.5)
  assert_true(is_finite(contrast.p_value))
  let summary = inference_summary([contrast])
  assert_eq(summary.contrasts, 1)
  assert_eq(inference_summary([]).contrasts, 0)
  let curve = inference_power_curve([10, 50, 100], 1.0, 1.0)
  assert_eq(curve.length(), 3)
  assert_true(curve[2].power >= curve[0].power)
  assert_true(inference_required_sample_size(1.0, 0.5) > 0)
  assert_eq(inference_containing_count([interval], 2.0), 1)
  assert_true(inference_rmse(2.0, [1.0, 2.0, 3.0]) > 0.0)
}

///|
test "production operational diagnostics" {
  let metric = operational_metric("quality", 0.95, 0.9, tolerance=0.1)
  assert_true(operational_metric_passes(metric))
  let snapshot = operational_snapshot("run-1", [metric])
  assert_true(snapshot.passes)
  assert_eq(operational_failures(snapshot).length(), 0)
  assert_true(operational_snapshot_is_finite(snapshot))
  let degraded = operational_snapshot("run-2", [
    operational_metric("quality", 0.5, 0.9, tolerance=0.0),
  ])
  let comparison = compare_operational_snapshots(snapshot, degraded)
  assert_true(comparison.degraded > 0)
  assert_true(
    workflow_degraded(
      workflow_report([workflow_check("quality", true, 1.0)]),
      workflow_report([workflow_check("quality", false, 0.2, severity=3)]),
    ),
  )
  assert_true(operational_mad([1.0, 1.0, 2.0, 3.0]) >= 0.0)
}

///|
test "production workflow preflight" {
  let data = CausalDataset::new(
    [[-2.0], [-1.0], [0.0], [1.0], [2.0], [3.0], [4.0], [5.0]],
    [false, true, false, true, false, true, false, true],
    [0.0, 2.0, 1.0, 3.0, 2.0, 4.0, 3.0, 5.0],
  )
  let graph = CausalGraph::new(["x", "t", "y"])
    .add_edge("x", "t")
    .add_edge("t", "y")
  let plan = analysis_plan(
    estimand="ATE",
    bootstrap_replicates=20,
    cross_fit_folds=2,
    minimum_quality_score=0.0,
  )
  let report = workflow_preflight(data, plan, graph, "t", "y", minimum_arm=2)
  assert_eq(report.checks.length(), 8)
  assert_true(report.score >= 0.0 && report.score <= 1.0)
  assert_eq(workflow_summary(report).length(), 7)
  assert_true(workflow_fingerprint(report) != 0UL)
  assert_eq(
    workflow_failed_names(report).length(),
    workflow_failure_count(report, 1),
  )
}

///|
test "production deterministic benchmark truth" {
  let data = CausalDataset::new(
    [[-1.0], [0.0], [1.0], [2.0], [3.0], [4.0], [5.0], [6.0]],
    [false, true, false, true, false, true, false, true],
    [0.0, 2.0, 1.0, 3.0, 2.0, 4.0, 3.0, 5.0],
  )
  let estimate = estimate_difference_in_means(data.outcome, data.treatment)
  assert_close(estimate.estimate, 2.0, 1.0e-12)
  assert_eq(estimate.sample_size, 8)
  assert_true(estimate.lower <= estimate.upper)
  let reports : Array[BenchmarkReport] = Array::new()
  for seed in [11UL, 12UL, 13UL, 14UL] {
    reports.push(
      run_synthetic_benchmark(SyntheticConfig::new(120, 3, 1.5, seed)),
    )
  }
  let operating = simulation_operating_point(reports)
  assert_eq(operating.scenarios, 4)
  assert_true(is_finite(operating.mean_estimate))
  assert_true(is_finite(operating.rmse))
  assert_true(operating.mean_effective_sample_size > 0.0)
  assert_true(operating.coverage >= 0.0 && operating.coverage <= 1.0)
}

///|
test "production batch empty and finite boundaries" {
  let empty : Array[BatchRun] = Array::new()
  let summary = summarize_batch(empty)
  assert_eq(summary.runs, 0)
  assert_eq(batch_estimates(empty), [])
  assert_eq(batch_estimate_range(empty), [0.0, 0.0])
  assert_eq(batch_failed_indices(empty), [])
  assert_true(!batch_release_check(empty, 0.8, 1.0))
  assert_true(!batch_all_successful_finite(empty))
  assert_eq(batch_summary_vector(empty).length(), 8)
  assert_eq(batch_running_mean(empty), [])
}
