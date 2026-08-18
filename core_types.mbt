///|
/// A rectangular observational data set used by the estimators in this module.
/// Rows are observations, columns are covariates, and treatment is binary.
pub struct CausalDataset {
  covariates : Array[Array[Double]]
  treatment : Array[Bool]
  outcome : Array[Double]
  feature_names : Array[String]
}

///|
/// Creates a data set with generated feature names.
pub fn CausalDataset::new(
  covariates : Array[Array[Double]],
  treatment : Array[Bool],
  outcome : Array[Double],
) -> CausalDataset {
  let width = if covariates.length() == 0 { 0 } else { covariates[0].length() }
  let names : Array[String] = Array::new(capacity=width)
  for j in 0..<width {
    names.push("x\{j + 1}")
  }
  { covariates, treatment, outcome, feature_names: names }
}

///|
/// Returns the number of observations.
pub fn CausalDataset::n(self : CausalDataset) -> Int {
  self.outcome.length()
}

///|
/// Returns the number of covariate columns.
pub fn CausalDataset::p(self : CausalDataset) -> Int {
  if self.covariates.length() == 0 {
    0
  } else {
    self.covariates[0].length()
  }
}

///|
/// Returns the number of treated observations.
pub fn CausalDataset::treated_count(self : CausalDataset) -> Int {
  let mut result = 0
  for value in self.treatment {
    if value {
      result += 1
    }
  }
  result
}

///|
/// Returns the number of control observations.
pub fn CausalDataset::control_count(self : CausalDataset) -> Int {
  self.n() - self.treated_count()
}

///|
/// Checks dimensions and basic finite-value constraints.
pub fn CausalDataset::is_valid(self : CausalDataset) -> Bool {
  if self.covariates.length() != self.n() {
    return false
  }
  if self.treatment.length() != self.n() {
    return false
  }
  if self.n() == 0 {
    return true
  }
  let width = self.p()
  if self.feature_names.length() != width {
    return false
  }
  for row in self.covariates {
    if row.length() != width {
      return false
    }
    for value in row {
      if !is_finite(value) {
        return false
      }
    }
  }
  for value in self.outcome {
    if !is_finite(value) {
      return false
    }
  }
  true
}

///|
/// Returns a row-filtered copy of a data set.
pub fn CausalDataset::select(
  self : CausalDataset,
  indices : Array[Int],
) -> CausalDataset {
  let x = Array::new(capacity=indices.length())
  let t = Array::new(capacity=indices.length())
  let y = Array::new(capacity=indices.length())
  for index in indices {
    if index >= 0 && index < self.n() {
      x.push(self.covariates[index])
      t.push(self.treatment[index])
      y.push(self.outcome[index])
    }
  }
  { covariates: x, treatment: t, outcome: y, feature_names: self.feature_names }
}

///|
/// A named point estimate with an uncertainty interval and diagnostic metadata.
pub struct Estimate {
  estimate : Double
  standard_error : Double
  lower : Double
  upper : Double
  sample_size : Int
  effective_sample_size : Double
  estimand : String
}

///|
/// Creates an estimate using a normal approximation interval.
pub fn Estimate::from_standard_error(
  estimate : Double,
  standard_error : Double,
  sample_size : Int,
  effective_sample_size : Double,
  estimand : String,
) -> Estimate {
  let margin = 1.959963984540054 * standard_error
  {
    estimate,
    standard_error,
    lower: estimate - margin,
    upper: estimate + margin,
    sample_size,
    effective_sample_size,
    estimand,
  }
}

///|
/// A coefficient vector and training diagnostics for a generalized linear model.
pub struct ModelFit {
  coefficients : Array[Double]
  intercept : Double
  iterations : Int
  converged : Bool
  loss : Double
  feature_means : Array[Double]
  feature_scales : Array[Double]
}

///|
/// A single matched treated-control pair.
pub struct MatchPair {
  treated_index : Int
  control_index : Int
  distance : Double
  weight : Double
}

///|
/// A balance diagnostic for one covariate.
pub struct BalanceMetric {
  name : String
  treated_mean : Double
  control_mean : Double
  standardized_difference : Double
  variance_ratio : Double
  absolute_difference : Double
  balanced : Bool
}

///|
/// Summary of a bootstrap distribution.
pub struct BootstrapSummary {
  point_estimate : Double
  standard_error : Double
  lower : Double
  upper : Double
  replicates : Int
  successful_replicates : Int
}

///|
/// A deterministic pseudo-random generator for reproducible simulations.
pub struct RandomState {
  mut state : UInt64
}

///|
pub fn RandomState::new(seed : UInt64) -> RandomState {
  { state: if seed == 0UL { 88172645463393265UL } else { seed } }
}

///|
pub fn RandomState::next_u64(self : RandomState) -> UInt64 {
  self.state = self.state * 2862933555777941757UL + 3037000493UL
  self.state
}

///|
pub fn RandomState::uniform(self : RandomState) -> Double {
  (self.next_u64() >> 11).to_double() / 9007199254740992.0
}

///|
pub fn RandomState::normal(self : RandomState) -> Double {
  let u1 = if self.uniform() < 1.0e-12 { 1.0e-12 } else { self.uniform() }
  let u2 = self.uniform()
  (-2.0 * @math.ln(u1)).sqrt() * @math.cos(2.0 * 3.141592653589793 * u2)
}
