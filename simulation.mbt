///|
pub struct SyntheticConfig {
  sample_size : Int
  covariate_count : Int
  treatment_effect : Double
  seed : UInt64
}

///|
pub fn SyntheticConfig::new(
  sample_size : Int,
  covariate_count : Int,
  treatment_effect : Double,
  seed : UInt64,
) -> SyntheticConfig {
  { sample_size, covariate_count, treatment_effect, seed }
}

///|
pub struct SyntheticData {
  dataset : CausalDataset
  true_ate : Double
  true_propensity : Array[Double]
}

///|
pub struct BenchmarkReport {
  sample_size : Int
  true_ate : Double
  estimated_ate : Double
  absolute_error : Double
  standard_error : Double
  effective_sample_size : Double
  propensity_auc : Double
}

///|
pub fn simulate_observational_data(config : SyntheticConfig) -> SyntheticData {
  let n = if config.sample_size < 0 { 0 } else { config.sample_size }
  let p = if config.covariate_count < 1 { 1 } else { config.covariate_count }
  let rng = RandomState::new(config.seed)
  let covariates : Array[Array[Double]] = Array::new(capacity=n)
  let treatment : Array[Bool] = Array::new(capacity=n)
  let outcomes : Array[Double] = Array::new(capacity=n)
  let propensities : Array[Double] = Array::new(capacity=n)
  for _ in 0..<n {
    let row = Array::new(capacity=p)
    for _ in 0..<p {
      row.push(rng.normal())
    }
    let second_covariate = if p > 1 { 0.2 * row[1] } else { 0.0 }
    let linear = 0.35 * row[0] + second_covariate
    let propensity = safe_probability(sigmoid(linear), epsilon=0.02)
    let assigned = rng.uniform() < propensity
    let second_outcome = if p > 1 { -0.5 * row[1] } else { 0.0 }
    let baseline = 1.0 + 0.8 * row[0] + second_outcome
    let treatment_contribution = if assigned {
      config.treatment_effect
    } else {
      0.0
    }
    let outcome = baseline + treatment_contribution + 0.4 * rng.normal()
    covariates.push(row)
    treatment.push(assigned)
    outcomes.push(outcome)
    propensities.push(propensity)
  }
  {
    dataset: CausalDataset::new(covariates, treatment, outcomes),
    true_ate: config.treatment_effect,
    true_propensity: propensities,
  }
}

///|
pub fn run_synthetic_benchmark(config : SyntheticConfig) -> BenchmarkReport {
  let simulated = simulate_observational_data(config)
  let model = fit_logistic_regression(
    simulated.dataset.covariates,
    simulated.dataset.treatment,
  )
  let estimated_propensity = predict_propensity(
    model,
    simulated.dataset.covariates,
  )
  let estimate = estimate_ipw_ate(
    simulated.dataset.outcome,
    simulated.dataset.treatment,
    estimated_propensity,
  )
  {
    sample_size: simulated.dataset.n(),
    true_ate: simulated.true_ate,
    estimated_ate: estimate.estimate,
    absolute_error: (estimate.estimate - simulated.true_ate).abs(),
    standard_error: estimate.standard_error,
    effective_sample_size: estimate.effective_sample_size,
    propensity_auc: auc_score(estimated_propensity, simulated.dataset.treatment),
  }
}
