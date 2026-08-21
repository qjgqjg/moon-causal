///|
/// Missingness pattern for one row.
pub struct MissingPattern {
  code : String
  missing_count : Int
  observed_count : Int
  frequency : Int
  complete : Bool
}

///|
/// Result of deterministic multiple imputation.
pub struct ImputationResult {
  datasets : Array[Array[Array[Double]]]
  imputations : Int
  changed_cells : Int
  seed : UInt64
  passes : Bool
}

///|
/// Rubin-style combination of scalar estimates.
pub struct RubinCombination {
  estimate : Double
  within_variance : Double
  between_variance : Double
  total_variance : Double
  standard_error : Double
  degrees_of_freedom : Double
}

///|
/// Missingness model diagnostics.
pub struct MissingnessAudit {
  rows : Int
  columns : Int
  missing_cells : Int
  complete_rows : Int
  pattern_count : Int
  maximum_pattern_frequency : Int
  passes : Bool
}

///|
fn md_pattern_code(row : Array[Double]) -> String {
  let builder = StringBuilder::new()
  for value in row {
    builder.write_string(if is_finite(value) { "1" } else { "0" })
  }
  builder.to_string()
}

///|
/// Returns binary missingness indicators for a matrix.
pub fn missingness_matrix(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    let pattern = Array::new(capacity=row.length())
    for value in row {
      pattern.push(if is_finite(value) { 0.0 } else { 1.0 })
    }
    result.push(pattern)
  }
  result
}

///|
/// Counts unique missingness patterns.
pub fn missingness_patterns(
  matrix : Array[Array[Double]],
) -> Array[MissingPattern] {
  let codes : Array[String] = Array::new()
  let counts : Array[Int] = Array::new()
  let missing_counts : Array[Int] = Array::new()
  let observed_counts : Array[Int] = Array::new()
  for row in matrix {
    let code = md_pattern_code(row)
    let mut missing = 0
    for value in row {
      if !is_finite(value) {
        missing += 1
      }
    }
    let mut index = -1
    for i in 0..<codes.length() {
      if codes[i] == code {
        index = i
        break
      }
    }
    if index < 0 {
      codes.push(code)
      counts.push(1)
      missing_counts.push(missing)
      observed_counts.push(row.length() - missing)
    } else {
      counts[index] += 1
    }
  }
  let result : Array[MissingPattern] = Array::new(capacity=codes.length())
  for i in 0..<codes.length() {
    let missing = missing_counts[i]
    let observed = observed_counts[i]
    result.push({
      code: codes[i],
      missing_count: missing,
      observed_count: observed,
      frequency: counts[i],
      complete: missing == 0,
    })
  }
  result
}

///|
/// Audits missing cells and missingness pattern concentration.
pub fn audit_missingness(matrix : Array[Array[Double]]) -> MissingnessAudit {
  let patterns = missingness_patterns(matrix)
  let mut missing = 0
  let mut complete = 0
  for row in matrix {
    for value in row {
      if !is_finite(value) {
        missing += 1
      }
    }
    if row.all(fn(value) { is_finite(value) }) {
      complete += 1
    }
  }
  let mut maximum = 0
  for pattern in patterns {
    if pattern.frequency > maximum {
      maximum = pattern.frequency
    }
  }
  {
    rows: matrix.length(),
    columns: if matrix.length() == 0 {
      0
    } else {
      matrix[0].length()
    },
    missing_cells: missing,
    complete_rows: complete,
    pattern_count: patterns.length(),
    maximum_pattern_frequency: maximum,
    passes: matrix.length() == 0 ||
    missing <
    matrix.length() *
    (if matrix.length() == 0 { 0 } else { matrix[0].length() }),
  }
}

///|
/// Returns a row-level complete-case mask.
pub fn complete_case_mask(matrix : Array[Array[Double]]) -> Array[Bool] {
  let result = Array::new(capacity=matrix.length())
  for row in matrix {
    let mut complete = true
    for value in row {
      if !is_finite(value) {
        complete = false
      }
    }
    result.push(complete)
  }
  result
}

///|
/// Imputes a matrix with column means and returns changed-cell count.
pub fn mean_impute_matrix(matrix : Array[Array[Double]]) -> ImputationResult {
  let profiles = profile_matrix(matrix)
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  let mut changed = 0
  for row in matrix {
    let imputed = row.copy()
    for j in 0..<row.length() {
      if !is_finite(imputed[j]) {
        imputed[j] = if j < profiles.length() { profiles[j].mean } else { 0.0 }
        changed += 1
      }
    }
    result.push(imputed)
  }
  {
    datasets: [result],
    imputations: 1,
    changed_cells: changed,
    seed: 0,
    passes: true,
  }
}

///|
/// Creates multiple deterministic imputations with small model-based jitter.
pub fn multiple_mean_imputation(
  matrix : Array[Array[Double]],
  imputations : Int,
  seed : UInt64,
) -> ImputationResult {
  let count = if imputations > 0 { imputations } else { 1 }
  let profiles = profile_matrix(matrix)
  let rng = RandomState::new(seed)
  let datasets : Array[Array[Array[Double]]] = Array::new(capacity=count)
  let mut changed = 0
  for _ in 0..<count {
    let dataset : Array[Array[Double]] = Array::new(capacity=matrix.length())
    for row in matrix {
      let imputed = row.copy()
      for j in 0..<row.length() {
        if !is_finite(imputed[j]) {
          let center = if j < profiles.length() {
            profiles[j].mean
          } else {
            0.0
          }
          let scale = if j < profiles.length() {
            profiles[j].std_dev
          } else {
            0.0
          }
          imputed[j] = center + scale * 0.1 * rng.normal()
          changed += 1
        }
      }
      dataset.push(imputed)
    }
    datasets.push(dataset)
  }
  { datasets, imputations: count, changed_cells: changed, seed, passes: true }
}

///|
/// Iteratively imputes with current column means until convergence.
pub fn iterative_mean_imputation(
  matrix : Array[Array[Double]],
  iterations : Int,
) -> Array[Array[Double]] {
  let result = mean_impute_matrix(matrix).datasets[0]
  let rounds = if iterations > 0 { iterations } else { 1 }
  for _ in 0..<rounds {
    let profiles = profile_matrix(result)
    for row in result {
      for j in 0..<row.length() {
        if !is_finite(row[j]) {
          row[j] = profiles[j].mean
        }
      }
    }
  }
  result
}

///|
/// Combines estimates and variances using Rubin's rules.
pub fn rubin_combine(
  estimates : Array[Double],
  variances : Array[Double],
) -> RubinCombination {
  let n = estimates.length().min(variances.length())
  if n == 0 {
    return {
      estimate: 0.0,
      within_variance: 0.0,
      between_variance: 0.0,
      total_variance: 0.0,
      standard_error: 0.0,
      degrees_of_freedom: 0.0,
    }
  }
  let estimate = mean(estimates[:n].to_owned())
  let within = mean(variances[:n].to_owned())
  let between = variance(estimates[:n].to_owned())
  let total = within + (1.0 + 1.0 / n.to_double()) * between
  let standard_error = total.max(0.0).sqrt()
  let degrees = if between == 0.0 {
    1.0e9
  } else {
    (n - 1).to_double() *
    @math.pow(1.0 + within / ((1.0 + 1.0 / n.to_double()) * between), 2.0)
  }
  {
    estimate,
    within_variance: within,
    between_variance: between,
    total_variance: total,
    standard_error,
    degrees_of_freedom: degrees,
  }
}

///|
/// Computes inverse-probability weights for a missing outcome indicator.
pub fn missingness_weights(
  observed : Array[Bool],
  observation_probability : Array[Double],
) -> Array[Double] {
  let n = observed.length().min(observation_probability.length())
  let result = Array::new(capacity=n)
  for i in 0..<n {
    let probability = safe_probability(observation_probability[i])
    result.push(if observed[i] { 1.0 / probability } else { 0.0 })
  }
  result
}

///|
/// Computes a missing-outcome weighted mean.
pub fn missingness_weighted_mean(
  outcomes : Array[Double],
  observed : Array[Bool],
  weights : Array[Double],
) -> Double {
  let n = outcomes.length().min(observed.length()).min(weights.length())
  let mut numerator = 0.0
  let mut denominator = 0.0
  for i in 0..<n {
    if observed[i] {
      numerator += outcomes[i] * weights[i]
      denominator += weights[i]
    }
  }
  if denominator == 0.0 {
    0.0
  } else {
    numerator / denominator
  }
}

///|
/// Creates missingness indicators for a causal dataset outcome and covariates.
pub fn causal_missingness_indicators(
  dataset : CausalDataset,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=dataset.n())
  for i in 0..<dataset.n() {
    let row = Array::new(capacity=dataset.p() + 1)
    row.push(if is_finite(dataset.outcome[i]) { 0.0 } else { 1.0 })
    for column in dataset.covariates {
      row.push(if is_finite(column[i]) { 0.0 } else { 1.0 })
    }
    result.push(row)
  }
  result
}

///|
/// Returns a complete causal dataset after dropping incomplete rows.
pub fn advanced_complete_case_dataset(dataset : CausalDataset) -> CausalDataset {
  let indices : Array[Int] = Array::new()
  for i in 0..<dataset.n() {
    let mut complete = is_finite(dataset.outcome[i])
    for column in dataset.covariates {
      if !is_finite(column[i]) {
        complete = false
      }
    }
    if complete {
      indices.push(i)
    }
  }
  dataset.select(indices)
}

///|
/// Returns a compact missingness summary vector.
pub fn missingness_summary(audit : MissingnessAudit) -> Array[Double] {
  [
    audit.rows.to_double(),
    audit.columns.to_double(),
    audit.missing_cells.to_double(),
    audit.complete_rows.to_double(),
    audit.pattern_count.to_double(),
    audit.maximum_pattern_frequency.to_double(),
    if audit.passes {
      1.0
    } else {
      0.0
    },
  ]
}
