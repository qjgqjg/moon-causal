///|
/// Data-quality issue emitted by a production analysis gate.
pub struct QualityIssue {
  code : String
  severity : Int
  score : Double
  message : String
}

///|
/// Descriptive profile for one numeric covariate.
pub struct ColumnProfile {
  index : Int
  count : Int
  missing : Int
  mean : Double
  std_dev : Double
  minimum : Double
  maximum : Double
  unique : Int
  zero_fraction : Double
  outlier_count : Int
}

///|
/// Dataset-level quality summary.
pub struct DatasetQuality {
  rows : Int
  columns : Int
  valid_rows : Int
  duplicate_rows : Int
  missing_cells : Int
  columns_profile : Array[ColumnProfile]
  issues : Array[QualityIssue]
  score : Double
}

///|
/// Distribution drift summary between a reference and current sample.
pub struct DriftReport {
  statistic : Double
  reference_size : Int
  current_size : Int
  threshold : Double
  drifted : Bool
  interpretation : String
}

///|
/// Positivity and overlap summary for a binary treatment.
pub struct PositivityProfile {
  treated : Int
  control : Int
  minimum_score : Double
  maximum_score : Double
  near_zero : Int
  near_one : Int
  effective_sample_size : Double
  passes : Bool
}

///|
fn finite_or(value : Double, fallback : Double) -> Double {
  if is_finite(value) {
    value
  } else {
    fallback
  }
}

///|
fn numeric_values(column : Array[Double]) -> Array[Double] {
  let result = Array::new(capacity=column.length())
  for value in column {
    if is_finite(value) {
      result.push(value)
    }
  }
  result
}

///|
fn min_value(values : Array[Double]) -> Double {
  if values.length() == 0 {
    0.0
  } else {
    let mut result = values[0]
    for value in values[1:] {
      if value < result {
        result = value
      }
    }
    result
  }
}

///|
fn max_value(values : Array[Double]) -> Double {
  if values.length() == 0 {
    0.0
  } else {
    let mut result = values[0]
    for value in values[1:] {
      if value > result {
        result = value
      }
    }
    result
  }
}

///|
fn sorted_copy(values : Array[Double]) -> Array[Double] {
  let result = values.copy()
  for i in 1..<result.length() {
    let value = result[i]
    let mut j = i
    while j > 0 && result[j - 1] > value {
      result[j] = result[j - 1]
      j -= 1
    }
    result[j] = value
  }
  result
}

///|
fn unique_count(values : Array[Double]) -> Int {
  if values.length() == 0 {
    return 0
  }
  let sorted = sorted_copy(values)
  let mut count = 1
  for i in 1..<sorted.length() {
    if sorted[i] != sorted[i - 1] {
      count += 1
    }
  }
  count
}

///|
fn median_value(values : Array[Double]) -> Double {
  let sorted = sorted_copy(values)
  let n = sorted.length()
  if n == 0 {
    0.0
  } else if n % 2 == 1 {
    sorted[n / 2]
  } else {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  }
}

///|
/// Returns a profile for one matrix column. Missing values are non-finite.
pub fn profile_column(column : Array[Double], index : Int) -> ColumnProfile {
  let values = numeric_values(column)
  let missing = column.length() - values.length()
  let average = mean_or(values, 0.0)
  let deviation = if values.length() > 1 { std_dev(values) } else { 0.0 }
  let lower = if values.length() > 3 {
    quantile(values, 0.25)
  } else {
    min_value(values)
  }
  let upper = if values.length() > 3 {
    quantile(values, 0.75)
  } else {
    max_value(values)
  }
  let iqr = upper - lower
  let low_fence = lower - 1.5 * iqr
  let high_fence = upper + 1.5 * iqr
  let mut outliers = 0
  let mut zeros = 0
  for value in values {
    if value == 0.0 {
      zeros += 1
    }
    if value < low_fence || value > high_fence {
      outliers += 1
    }
  }
  {
    index,
    count: values.length(),
    missing,
    mean: average,
    std_dev: finite_or(deviation, 0.0),
    minimum: min_value(values),
    maximum: max_value(values),
    unique: unique_count(values),
    zero_fraction: if values.length() == 0 {
      0.0
    } else {
      zeros.to_double() / values.length().to_double()
    },
    outlier_count: outliers,
  }
}

///|
/// Profiles each column in a row-major matrix.
pub fn profile_matrix(matrix : Array[Array[Double]]) -> Array[ColumnProfile] {
  if matrix.length() == 0 {
    return []
  }
  let columns = matrix[0].length()
  let result : Array[ColumnProfile] = Array::new(capacity=columns)
  for column_index in 0..<columns {
    let column = Array::new(capacity=matrix.length())
    for row in matrix {
      if column_index < row.length() {
        column.push(row[column_index])
      }
    }
    result.push(profile_column(column, column_index))
  }
  result
}

///|
/// Counts non-finite entries in a column.
pub fn missing_count(column : Array[Double]) -> Int {
  let mut result = 0
  for value in column {
    if !is_finite(value) {
      result += 1
    }
  }
  result
}

///|
/// Returns the number of missing cells in each row.
pub fn row_missing_counts(matrix : Array[Array[Double]]) -> Array[Int] {
  let result : Array[Int] = Array::new(capacity=matrix.length())
  for row in matrix {
    let mut count = 0
    for value in row {
      if !is_finite(value) {
        count += 1
      }
    }
    result.push(count)
  }
  result
}

///|
/// Returns a binary missingness indicator for a column.
pub fn missing_indicator(column : Array[Double]) -> Array[Double] {
  let result = Array::new(capacity=column.length())
  for value in column {
    result.push(if is_finite(value) { 0.0 } else { 1.0 })
  }
  result
}

///|
/// Imputes non-finite values with the observed mean.
pub fn impute_mean(column : Array[Double]) -> Array[Double] {
  let observed = numeric_values(column)
  let replacement = mean_or(observed, 0.0)
  let result = column.copy()
  for i in 0..<result.length() {
    if !is_finite(result[i]) {
      result[i] = replacement
    }
  }
  result
}

///|
/// Imputes non-finite values with the observed median.
pub fn impute_median(column : Array[Double]) -> Array[Double] {
  let replacement = median_value(numeric_values(column))
  let result = column.copy()
  for i in 0..<result.length() {
    if !is_finite(result[i]) {
      result[i] = replacement
    }
  }
  result
}

///|
/// Imputes non-finite values with a supplied constant.
pub fn impute_constant(
  column : Array[Double],
  replacement : Double,
) -> Array[Double] {
  let result = column.copy()
  for i in 0..<result.length() {
    if !is_finite(result[i]) {
      result[i] = replacement
    }
  }
  result
}

///|
/// Marks observations outside Tukey's 1.5-IQR fences.
pub fn iqr_outlier_flags(
  column : Array[Double],
  multiplier? : Double = 1.5,
) -> Array[Bool] {
  let values = numeric_values(column)
  if values.length() < 2 {
    return Array::make(column.length(), false)
  }
  let lower = quantile(values, 0.25)
  let upper = quantile(values, 0.75)
  let spread = upper - lower
  let low = lower - multiplier * spread
  let high = upper + multiplier * spread
  let result = Array::new(capacity=column.length())
  for value in column {
    result.push(!is_finite(value) || value < low || value > high)
  }
  result
}

///|
/// Marks observations farther than a supplied number of standard deviations.
pub fn zscore_outlier_flags(
  column : Array[Double],
  threshold? : Double = 3.0,
) -> Array[Bool] {
  let values = numeric_values(column)
  let average = mean_or(values, 0.0)
  let deviation = std_dev(values)
  let result = Array::new(capacity=column.length())
  for value in column {
    result.push(
      !is_finite(value) ||
      deviation == 0.0 ||
      (value - average).abs() > threshold * deviation,
    )
  }
  result
}

///|
/// Replaces outliers with the nearest Tukey fence.
pub fn winsorize_iqr(
  column : Array[Double],
  multiplier? : Double = 1.5,
) -> Array[Double] {
  let values = numeric_values(column)
  if values.length() < 2 {
    return impute_median(column)
  }
  let lower = quantile(values, 0.25)
  let upper = quantile(values, 0.75)
  let spread = upper - lower
  let low = lower - multiplier * spread
  let high = upper + multiplier * spread
  let result = impute_median(column)
  for i in 0..<result.length() {
    result[i] = clamp(result[i], low, high)
  }
  result
}

///|
/// Computes a stable row hash for audit and duplicate detection.
pub fn row_fingerprint(
  row : Array[Double],
  precision? : Double = 1.0e6,
) -> UInt64 {
  let mut hash : UInt64 = 1469598103934665603
  for value in row {
    let scaled = if is_finite(value) {
      (value * precision).round()
    } else {
      -922337203685477.0
    }
    let integer = scaled.to_int().to_uint64()
    hash = (hash ^ integer) * 1099511628211
    hash = (hash ^ 0xff) * 1099511628211
  }
  hash
}

///|
/// Returns a fingerprint for every row.
pub fn row_fingerprints(matrix : Array[Array[Double]]) -> Array[UInt64] {
  let result = Array::new(capacity=matrix.length())
  for row in matrix {
    result.push(row_fingerprint(row))
  }
  result
}

///|
/// Flags rows whose rounded values have appeared earlier in the sample.
pub fn duplicate_row_flags(
  matrix : Array[Array[Double]],
  precision? : Double = 1.0e6,
) -> Array[Bool] {
  let seen : Array[UInt64] = Array::new()
  let result = Array::new(capacity=matrix.length())
  for row in matrix {
    let fingerprint = row_fingerprint(row, precision~)
    let duplicate = seen.contains(fingerprint)
    result.push(duplicate)
    if !duplicate {
      seen.push(fingerprint)
    }
  }
  result
}

///|
/// Computes a deterministic order-sensitive matrix checksum.
pub fn matrix_checksum(matrix : Array[Array[Double]]) -> UInt64 {
  let mut checksum : UInt64 = 2166136261
  for row in matrix {
    checksum = (checksum ^ row_fingerprint(row)) * 16777619
  }
  checksum
}

///|
fn histogram(
  values : Array[Double],
  bins : Int,
  lower : Double,
  upper : Double,
) -> Array[Int] {
  let result = Array::make(bins, 0)
  if bins <= 0 || upper <= lower {
    return result
  }
  let width = (upper - lower) / bins.to_double()
  for value in values {
    if is_finite(value) {
      let raw = ((value - lower) / width).to_int()
      let index = if raw < 0 { 0 } else if raw >= bins { bins - 1 } else { raw }
      result[index] += 1
    }
  }
  result
}

///|
/// Computes population stability index using fixed reference quantile bins.
pub fn population_stability_index(
  reference : Array[Double],
  current : Array[Double],
  bins? : Int = 10,
) -> Double {
  let reference_observed = numeric_values(reference)
  let current_observed = numeric_values(current)
  if reference_observed.length() == 0 ||
    current_observed.length() == 0 ||
    bins <= 0 {
    return 0.0
  }
  let lower = min_value(reference_observed)
  let upper = max_value(reference_observed)
  if upper <= lower {
    return 0.0
  }
  let expected = histogram(reference_observed, bins, lower, upper)
  let actual = histogram(current_observed, bins, lower, upper)
  let expected_total = reference_observed.length().to_double()
  let actual_total = current_observed.length().to_double()
  let mut psi = 0.0
  for i in 0..<bins {
    let p = (expected[i].to_double() + 0.5) /
      (expected_total + bins.to_double() * 0.5)
    let q = (actual[i].to_double() + 0.5) /
      (actual_total + bins.to_double() * 0.5)
    psi += (q - p) * @math.ln(q / p)
  }
  psi
}

///|
/// Computes the two-sample Kolmogorov-Smirnov distance.
pub fn kolmogorov_smirnov_distance(
  reference : Array[Double],
  current : Array[Double],
) -> Double {
  let first = sorted_copy(numeric_values(reference))
  let second = sorted_copy(numeric_values(current))
  if first.length() == 0 || second.length() == 0 {
    return 0.0
  }
  let mut i = 0
  let mut j = 0
  let mut distance = 0.0
  while i < first.length() && j < second.length() {
    if first[i] <= second[j] {
      i += 1
    } else {
      j += 1
    }
    let left = i.to_double() / first.length().to_double()
    let right = j.to_double() / second.length().to_double()
    let gap = (left - right).abs()
    if gap > distance {
      distance = gap
    }
  }
  distance
}

///|
/// Builds a drift report from the PSI statistic.
pub fn psi_drift_report(
  reference : Array[Double],
  current : Array[Double],
  threshold? : Double = 0.2,
) -> DriftReport {
  let statistic = population_stability_index(reference, current)
  {
    statistic,
    reference_size: reference.length(),
    current_size: current.length(),
    threshold,
    drifted: statistic > threshold,
    interpretation: if statistic > threshold {
      "material distribution drift"
    } else {
      "no material distribution drift"
    },
  }
}

///|
/// Builds a KS drift report.
pub fn ks_drift_report(
  reference : Array[Double],
  current : Array[Double],
  threshold? : Double = 0.1,
) -> DriftReport {
  let statistic = kolmogorov_smirnov_distance(reference, current)
  {
    statistic,
    reference_size: reference.length(),
    current_size: current.length(),
    threshold,
    drifted: statistic > threshold,
    interpretation: if statistic > threshold {
      "material empirical CDF drift"
    } else {
      "no material empirical CDF drift"
    },
  }
}

///|
/// Scores a matrix with operational quality rules.
pub fn assess_matrix(
  matrix : Array[Array[Double]],
  missing_limit? : Double = 0.2,
  outlier_limit? : Double = 0.05,
) -> DatasetQuality {
  let rows = matrix.length()
  let columns = if rows == 0 { 0 } else { matrix[0].length() }
  let profiles = profile_matrix(matrix)
  let mut missing_cells = 0
  let mut valid_rows = 0
  for row in matrix {
    let mut row_missing = 0
    for value in row {
      if !is_finite(value) {
        row_missing += 1
        missing_cells += 1
      }
    }
    if row_missing == 0 {
      valid_rows += 1
    }
  }
  let duplicate_flags = duplicate_row_flags(matrix)
  let mut duplicates = 0
  for duplicate in duplicate_flags {
    if duplicate {
      duplicates += 1
    }
  }
  let issues : Array[QualityIssue] = Array::new()
  for profile in profiles {
    let missing_fraction = if rows == 0 {
      0.0
    } else {
      profile.missing.to_double() / rows.to_double()
    }
    let outlier_fraction = if profile.count == 0 {
      0.0
    } else {
      profile.outlier_count.to_double() / profile.count.to_double()
    }
    if missing_fraction > missing_limit {
      issues.push({
        code: "MISSINGNESS",
        severity: 2,
        score: missing_fraction,
        message: "column exceeds missingness limit",
      })
    }
    if outlier_fraction > outlier_limit {
      issues.push({
        code: "OUTLIERS",
        severity: 1,
        score: outlier_fraction,
        message: "column contains many Tukey outliers",
      })
    }
    if profile.unique <= 1 && profile.count > 0 {
      issues.push({
        code: "CONSTANT",
        severity: 1,
        score: 1.0,
        message: "column has no observed variation",
      })
    }
  }
  if duplicates > 0 {
    issues.push({
      code: "DUPLICATES",
      severity: 1,
      score: duplicates.to_double() / rows.max(1).to_double(),
      message: "duplicate rows were detected",
    })
  }
  let mut penalty = 0.0
  for issue in issues {
    penalty += issue.score * issue.severity.to_double()
  }
  let denominator = if columns == 0 { 1.0 } else { columns.to_double() }
  let score = clamp(1.0 - penalty / denominator, 0.0, 1.0)
  {
    rows,
    columns,
    valid_rows,
    duplicate_rows: duplicates,
    missing_cells,
    columns_profile: profiles,
    issues,
    score,
  }
}

///|
/// Returns whether a quality result meets minimum operational requirements.
pub fn quality_gate(
  quality : DatasetQuality,
  minimum_score? : Double = 0.8,
  maximum_issues? : Int = 3,
) -> Bool {
  quality.score >= minimum_score && quality.issues.length() <= maximum_issues
}

///|
/// Profiles treatment positivity from estimated propensity scores.
pub fn positivity_profile(
  treatment : Array[Bool],
  propensity : Array[Double],
  cutoff? : Double = 0.02,
) -> PositivityProfile {
  let n = if treatment.length() < propensity.length() {
    treatment.length()
  } else {
    propensity.length()
  }
  if n == 0 {
    return {
      treated: 0,
      control: 0,
      minimum_score: 0.0,
      maximum_score: 0.0,
      near_zero: 0,
      near_one: 0,
      effective_sample_size: 0.0,
      passes: false,
    }
  }
  let mut treated = 0
  let mut control = 0
  let mut near_zero = 0
  let mut near_one = 0
  let mut minimum = 1.0
  let mut maximum = 0.0
  let weights = Array::new(capacity=n)
  for i in 0..<n {
    let score = safe_probability(propensity[i])
    if treatment[i] {
      treated += 1
    } else {
      control += 1
    }
    if score < cutoff {
      near_zero += 1
    }
    if score > 1.0 - cutoff {
      near_one += 1
    }
    if score < minimum {
      minimum = score
    }
    if score > maximum {
      maximum = score
    }
    weights.push(if treatment[i] { 1.0 / score } else { 1.0 / (1.0 - score) })
  }
  let ess = effective_sample_size(weights)
  {
    treated,
    control,
    minimum_score: minimum,
    maximum_score: maximum,
    near_zero,
    near_one,
    effective_sample_size: ess,
    passes: near_zero == 0 && near_one == 0 && treated > 0 && control > 0,
  }
}

///|
/// Detects rows with a missing treatment or outcome.
pub fn valid_causal_rows(dataset : CausalDataset) -> Array[Bool] {
  let result = Array::new(capacity=dataset.n())
  for i in 0..<dataset.n() {
    let mut valid = is_finite(dataset.outcome[i])
    for column in dataset.covariates {
      if !is_finite(column[i]) {
        valid = false
      }
    }
    result.push(valid)
  }
  result
}

///|
/// Counts row-wise complete cases for a causal dataset.
pub fn complete_case_count(dataset : CausalDataset) -> Int {
  let flags = valid_causal_rows(dataset)
  let mut count = 0
  for flag in flags {
    if flag {
      count += 1
    }
  }
  count
}

///|
/// Creates a compact quality summary vector for monitoring.
pub fn quality_summary_vector(quality : DatasetQuality) -> Array[Double] {
  [
    quality.rows.to_double(),
    quality.columns.to_double(),
    quality.valid_rows.to_double(),
    quality.duplicate_rows.to_double(),
    quality.missing_cells.to_double(),
    quality.score,
    quality.issues.length().to_double(),
  ]
}

///|
/// Returns columns with at least the requested observed variation.
pub fn variable_columns(
  profiles : Array[ColumnProfile],
  minimum_unique? : Int = 2,
) -> Array[Int] {
  let result : Array[Int] = Array::new()
  for profile in profiles {
    if profile.unique >= minimum_unique && profile.count > 0 {
      result.push(profile.index)
    }
  }
  result
}

///|
/// Calculates a covariance-like missingness association between two indicators.
pub fn missingness_association(
  first : Array[Double],
  second : Array[Double],
) -> Double {
  if first.length() != second.length() || first.length() < 2 {
    return 0.0
  }
  correlation(missing_indicator(first), missing_indicator(second))
}

///|
/// Produces a reproducibility-friendly quality fingerprint.
pub fn quality_fingerprint(quality : DatasetQuality) -> UInt64 {
  let vector = quality_summary_vector(quality)
  let rows : Array[Array[Double]] = [vector]
  matrix_checksum(rows)
}
