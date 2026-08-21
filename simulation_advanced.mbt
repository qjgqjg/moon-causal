///|
/// Scenario specification for stress-testing causal estimators.
pub struct SimulationScenario {
  sample_size : Int
  covariates : Int
  treatment_effect : Double
  confounding_strength : Double
  nonlinear : Bool
  noise_scale : Double
  seed : UInt64
}

///|
/// Potential-outcome simulation output.
pub struct PotentialOutcomeSample {
  covariates : Array[Array[Double]]
  treatment : Array[Bool]
  outcome_treated : Array[Double]
  outcome_control : Array[Double]
  outcome_observed : Array[Double]
  propensity : Array[Double]
  true_ate : Double
}

///|
/// Repeated simulation operating characteristics.
pub struct SimulationOperatingPoint {
  scenarios : Int
  mean_estimate : Double
  bias : Double
  monte_carlo_standard_error : Double
  rmse : Double
  coverage : Double
  mean_effective_sample_size : Double
  passes : Bool
}

///|
/// Creates a simulation scenario with safe bounds.
pub fn simulation_scenario(
  sample_size? : Int = 500,
  covariates? : Int = 4,
  treatment_effect? : Double = 1.0,
  confounding_strength? : Double = 0.5,
  nonlinear? : Bool = false,
  noise_scale? : Double = 1.0,
  seed? : UInt64 = 20260819,
) -> SimulationScenario {
  {
    sample_size: sample_size.max(0),
    covariates: covariates.max(1),
    treatment_effect,
    confounding_strength,
    nonlinear,
    noise_scale: noise_scale.max(0.0),
    seed,
  }
}

///|
/// Simulates binary-treatment potential outcomes with known truth.
pub fn simulate_potential_outcomes(
  scenario : SimulationScenario,
) -> PotentialOutcomeSample {
  let rng = RandomState::new(scenario.seed)
  let covariates : Array[Array[Double]] = Array::new(
    capacity=scenario.sample_size,
  )
  let treatment : Array[Bool] = Array::new(capacity=scenario.sample_size)
  let treated_outcome : Array[Double] = Array::new(
    capacity=scenario.sample_size,
  )
  let control_outcome : Array[Double] = Array::new(
    capacity=scenario.sample_size,
  )
  let observed : Array[Double] = Array::new(capacity=scenario.sample_size)
  let propensity : Array[Double] = Array::new(capacity=scenario.sample_size)
  for _ in 0..<scenario.sample_size {
    let row = Array::new(capacity=scenario.covariates)
    for _ in 0..<scenario.covariates {
      row.push(rng.normal())
    }
    let nonlinear_term = if scenario.nonlinear {
      0.5 * row[0] * row[0] - 0.25 * row[0].abs()
    } else {
      0.0
    }
    let confounding = scenario.confounding_strength * row[0]
    let probability = safe_probability(sigmoid(confounding), epsilon=0.02)
    let assigned = rng.uniform() < probability
    let second_covariate = if scenario.covariates > 1 {
      0.5 * row[1]
    } else {
      0.0
    }
    let baseline = nonlinear_term +
      row[0] +
      second_covariate +
      scenario.noise_scale * rng.normal()
    let heterogeneity = if scenario.nonlinear && scenario.covariates > 1 {
      0.2 * row[1]
    } else {
      0.0
    }
    let control_value = baseline
    let treated_value = baseline + scenario.treatment_effect + heterogeneity
    covariates.push(row)
    treatment.push(assigned)
    treated_outcome.push(treated_value)
    control_outcome.push(control_value)
    observed.push(if assigned { treated_value } else { control_value })
    propensity.push(probability)
  }
  {
    covariates,
    treatment,
    outcome_treated: treated_outcome,
    outcome_control: control_outcome,
    outcome_observed: observed,
    propensity,
    true_ate: scenario.treatment_effect,
  }
}

///|
/// Computes the realized individual treatment effects.
pub fn realized_effects(sample : PotentialOutcomeSample) -> Array[Double] {
  let n = sample.outcome_treated.length().min(sample.outcome_control.length())
  let result = Array::new(capacity=n)
  for i in 0..<n {
    result.push(sample.outcome_treated[i] - sample.outcome_control[i])
  }
  result
}

///|
/// Computes an operating point from repeated benchmark reports.
pub fn simulation_operating_point(
  reports : Array[BenchmarkReport],
  confidence_level? : Double = 0.95,
) -> SimulationOperatingPoint {
  if reports.length() == 0 {
    return {
      scenarios: 0,
      mean_estimate: 0.0,
      bias: 0.0,
      monte_carlo_standard_error: 0.0,
      rmse: 0.0,
      coverage: 0.0,
      mean_effective_sample_size: 0.0,
      passes: false,
    }
  }
  let estimates = Array::new(capacity=reports.length())
  let errors = Array::new(capacity=reports.length())
  let mut covered = 0
  let critical_value = if confidence_level >= 0.99 {
    2.5758
  } else if confidence_level >= 0.9 {
    1.96
  } else {
    1.645
  }
  let mut ess = 0.0
  for report in reports {
    estimates.push(report.estimated_ate())
    errors.push(report.absolute_error())
    ess += report.effective_sample_size()
    let margin = critical_value * report.standard_error()
    if report.true_ate() >= report.estimated_ate() - margin &&
      report.true_ate() <= report.estimated_ate() + margin {
      covered += 1
    }
  }
  let average = mean(estimates)
  let mut squared = 0.0
  for estimate in estimates {
    let difference = estimate - reports[0].true_ate()
    squared += difference * difference
  }
  let coverage = covered.to_double() / reports.length().to_double()
  {
    scenarios: reports.length(),
    mean_estimate: average,
    bias: average - reports[0].true_ate(),
    monte_carlo_standard_error: std_dev(estimates) /
    reports.length().to_double().sqrt(),
    rmse: (squared / reports.length().to_double()).sqrt(),
    coverage,
    mean_effective_sample_size: ess / reports.length().to_double(),
    passes: coverage >= 0.8 && is_finite(average),
  }
}

///|
/// Runs a deterministic repeated synthetic benchmark.
pub fn repeated_synthetic_benchmark(
  scenario : SimulationScenario,
  repetitions : Int,
) -> SimulationOperatingPoint {
  let reports : Array[BenchmarkReport] = Array::new(
    capacity=if repetitions > 0 { repetitions } else { 0 },
  )
  for i in 0..<repetitions {
    let run_scenario = {
      sample_size: scenario.sample_size,
      covariates: scenario.covariates,
      treatment_effect: scenario.treatment_effect,
      confounding_strength: scenario.confounding_strength,
      nonlinear: scenario.nonlinear,
      noise_scale: scenario.noise_scale,
      seed: scenario.seed + i.to_uint64(),
    }
    let simulated = simulate_potential_outcomes(run_scenario)
    let dataset = CausalDataset::new(
      simulated.covariates,
      simulated.treatment,
      simulated.outcome_observed,
    )
    let model = fit_logistic_regression(dataset.covariates, dataset.treatment)
    let propensity = predict_propensity(model, dataset.covariates)
    let estimate = estimate_ipw_ate(
      dataset.outcome,
      dataset.treatment,
      propensity,
    )
    reports.push({
      sample_size: dataset.n(),
      true_ate: simulated.true_ate,
      estimated_ate: estimate.estimate,
      absolute_error: (estimate.estimate - simulated.true_ate).abs(),
      standard_error: estimate.standard_error,
      effective_sample_size: effective_sample_size(
        stabilized_ipw_weights(dataset.treatment, propensity),
      ),
      propensity_auc: auc_score(propensity, dataset.treatment),
    })
  }
  simulation_operating_point(reports)
}

///|
/// Computes empirical coverage of normal intervals.
pub fn empirical_interval_coverage(
  truth : Double,
  estimates : Array[Double],
  standard_errors : Array[Double],
  critical_value? : Double = 1.96,
) -> Double {
  let n = estimates.length().min(standard_errors.length())
  if n == 0 {
    return 0.0
  }
  let mut covered = 0
  for i in 0..<n {
    if truth >= estimates[i] - critical_value * standard_errors[i] &&
      truth <= estimates[i] + critical_value * standard_errors[i] {
      covered += 1
    }
  }
  covered.to_double() / n.to_double()
}

///|
/// Computes Monte Carlo bias and RMSE against a known truth.
pub fn simulation_bias_rmse(
  truth : Double,
  estimates : Array[Double],
) -> Array[Double] {
  if estimates.length() == 0 {
    return [0.0, 0.0]
  }
  let average = mean(estimates)
  let mut squared = 0.0
  for estimate in estimates {
    squared += (estimate - truth) * (estimate - truth)
  }
  [average - truth, (squared / estimates.length().to_double()).sqrt()]
}

///|
/// Returns a compact simulation summary vector.
pub fn simulation_summary(point : SimulationOperatingPoint) -> Array[Double] {
  [
    point.scenarios.to_double(),
    point.mean_estimate,
    point.bias,
    point.monte_carlo_standard_error,
    point.rmse,
    point.coverage,
    point.mean_effective_sample_size,
    if point.passes {
      1.0
    } else {
      0.0
    },
  ]
}
