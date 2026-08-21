///|
/// Feature matrix with reproducible column metadata.
pub struct FeatureSet {
  values : Array[Array[Double]]
  names : Array[String]
  row_count : Int
  column_count : Int
  fingerprint : UInt64
}

///|
/// Integer-category encoding map.
pub struct CategoryEncoding {
  categories : Array[Int]
  values : Array[Double]
  unknown_value : Double
  ordered : Bool
}

///|
/// Feature screening result.
pub struct FeatureScreening {
  selected : Array[Int]
  dropped : Array[Int]
  variance : Array[Double]
  outcome_correlation : Array[Double]
  reasons : Array[String]
}

///|
/// Train-validation split with explicit row indices.
pub struct FeatureSplit {
  train_indices : Array[Int]
  validation_indices : Array[Int]
  train_values : Array[Array[Double]]
  validation_values : Array[Array[Double]]
}

///|
/// Monotonic transformation specification.
pub struct TransformSpec {
  center : Double
  scale : Double
  power : Double
  lower : Double
  upper : Double
  missing_value : Double
}

///|
fn fe_unique_int(values : Array[Int]) -> Array[Int] {
  let result : Array[Int] = Array::new()
  for value in values {
    if !result.contains(value) {
      result.push(value)
    }
  }
  result
}

///|
fn fe_sorted_int(values : Array[Int]) -> Array[Int] {
  let result = values.copy()
  for i in 1..<result.length() {
    let value = result[i]
    let mut cursor = i
    while cursor > 0 && result[cursor - 1] > value {
      result[cursor] = result[cursor - 1]
      cursor -= 1
    }
    result[cursor] = value
  }
  result
}

///|
fn fe_column(matrix : Array[Array[Double]], index : Int) -> Array[Double] {
  let result = Array::new(capacity=matrix.length())
  for row in matrix {
    if index < row.length() {
      result.push(row[index])
    }
  }
  result
}

///|
/// Creates a feature set and computes a stable fingerprint.
pub fn feature_set(
  values : Array[Array[Double]],
  names? : Array[String] = [],
) -> FeatureSet {
  let row_count = values.length()
  let column_count = if row_count == 0 { 0 } else { values[0].length() }
  let actual_names : Array[String] = Array::new(capacity=column_count)
  if names.length() == column_count {
    for name in names {
      actual_names.push(name)
    }
  } else {
    for i in 0..<column_count {
      actual_names.push("x\{i + 1}")
    }
  }
  {
    values: values.copy(),
    names: actual_names,
    row_count,
    column_count,
    fingerprint: matrix_checksum(values),
  }
}

///|
/// Adds deterministic missingness indicators for every numeric column.
pub fn add_missing_features(matrix : Array[Array[Double]]) -> FeatureSet {
  if matrix.length() == 0 {
    return feature_set([])
  }
  let width = matrix[0].length()
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    let expanded = row.copy()
    for value in row {
      expanded.push(if is_finite(value) { 0.0 } else { 1.0 })
    }
    result.push(expanded)
  }
  let names : Array[String] = Array::new(capacity=width * 2)
  for i in 0..<width {
    names.push("x\{i + 1}")
  }
  for i in 0..<width {
    names.push("x\{i + 1}_missing")
  }
  feature_set(result, names~)
}

///|
/// Encodes integer categories as one-hot columns.
pub fn one_hot_encode(
  categories : Array[Int],
  levels? : Array[Int] = [],
) -> FeatureSet {
  let actual_levels = if levels.length() > 0 {
    fe_sorted_int(fe_unique_int(levels))
  } else {
    fe_sorted_int(fe_unique_int(categories))
  }
  let result : Array[Array[Double]] = Array::new(capacity=categories.length())
  for category in categories {
    let row = Array::make(actual_levels.length(), 0.0)
    for j in 0..<actual_levels.length() {
      if actual_levels[j] == category {
        row[j] = 1.0
      }
    }
    result.push(row)
  }
  let names : Array[String] = Array::new(capacity=actual_levels.length())
  for level in actual_levels {
    names.push("category_\{level}")
  }
  feature_set(result, names~)
}

///|
/// Encodes integer categories using observed frequency.
pub fn frequency_encode(
  categories : Array[Int],
  unknown_value? : Double = 0.0,
) -> Array[Double] {
  let levels = fe_unique_int(categories)
  let counts = Array::make(levels.length(), 0)
  for category in categories {
    for j in 0..<levels.length() {
      if levels[j] == category {
        counts[j] += 1
      }
    }
  }
  let result = Array::new(capacity=categories.length())
  for category in categories {
    let mut encoded = unknown_value
    for j in 0..<levels.length() {
      if levels[j] == category {
        encoded = counts[j].to_double() / categories.length().max(1).to_double()
      }
    }
    result.push(encoded)
  }
  result
}

///|
/// Builds a category-to-value encoding map.
pub fn category_encoding(
  categories : Array[Int],
  values : Array[Double],
  unknown_value? : Double = 0.0,
  ordered? : Bool = false,
) -> CategoryEncoding {
  let n = categories.length().min(values.length())
  let levels = fe_unique_int(categories[:n].to_owned())
  let encoded = Array::make(levels.length(), unknown_value)
  let counts = Array::make(levels.length(), 0)
  for i in 0..<n {
    for j in 0..<levels.length() {
      if levels[j] == categories[i] {
        encoded[j] += values[i]
        counts[j] += 1
      }
    }
  }
  for j in 0..<levels.length() {
    if counts[j] > 0 {
      encoded[j] /= counts[j].to_double()
    }
  }
  { categories: levels, values: encoded, unknown_value, ordered }
}

///|
/// Applies a category encoding map to new category values.
pub fn apply_category_encoding(
  encoding : CategoryEncoding,
  categories : Array[Int],
) -> Array[Double] {
  let result = Array::new(capacity=categories.length())
  for category in categories {
    let mut encoded = encoding.unknown_value
    for j in 0..<encoding.categories.length() {
      if encoding.categories[j] == category {
        encoded = encoding.values[j]
      }
    }
    result.push(encoded)
  }
  result
}

///|
/// Computes smoothed target encoding for integer categories.
pub fn smoothed_target_encoding(
  categories : Array[Int],
  outcome : Array[Double],
  smoothing? : Double = 10.0,
) -> Array[Double] {
  let n = categories.length().min(outcome.length())
  let global = mean_or(outcome[:n].to_owned(), 0.0)
  let levels = fe_unique_int(categories[:n].to_owned())
  let sums = Array::make(levels.length(), 0.0)
  let counts = Array::make(levels.length(), 0)
  for i in 0..<n {
    for j in 0..<levels.length() {
      if levels[j] == categories[i] {
        sums[j] += outcome[i]
        counts[j] += 1
      }
    }
  }
  let encoded = Array::new(capacity=n)
  for i in 0..<n {
    let mut value = global
    for j in 0..<levels.length() {
      if levels[j] == categories[i] {
        value = (sums[j] + smoothing * global) /
          (counts[j].to_double() + smoothing)
      }
    }
    encoded.push(value)
  }
  encoded
}

///|
/// Maps a continuous variable into quantile bins.
pub fn quantile_bin(values : Array[Double], bins : Int) -> Array[Int] {
  let width = if bins > 0 { bins } else { 1 }
  let observed = values.filter(fn(value) { is_finite(value) })
  let sorted = observed.copy()
  for i in 1..<sorted.length() {
    let value = sorted[i]
    let mut cursor = i
    while cursor > 0 && sorted[cursor - 1] > value {
      sorted[cursor] = sorted[cursor - 1]
      cursor -= 1
    }
    sorted[cursor] = value
  }
  let cutpoints = Array::new(capacity=width - 1)
  for bin in 1..<width {
    cutpoints.push(
      if sorted.length() == 0 {
        0.0
      } else {
        sorted[(bin * sorted.length() / width).min(sorted.length() - 1)]
      },
    )
  }
  let result = Array::new(capacity=values.length())
  for value in values {
    if !is_finite(value) {
      result.push(-1)
    } else {
      let mut bucket = 0
      for cutpoint in cutpoints {
        if value >= cutpoint {
          bucket += 1
        }
      }
      result.push(bucket.min(width - 1))
    }
  }
  result
}

///|
/// Computes a piecewise-linear spline basis with knots.
pub fn linear_spline_basis(
  values : Array[Double],
  knots : Array[Double],
) -> FeatureSet {
  let result : Array[Array[Double]] = Array::new(capacity=values.length())
  for value in values {
    let row = Array::new(capacity=knots.length() + 2)
    row.push(value)
    row.push(value * value)
    for knot in knots {
      row.push(if value > knot { value - knot } else { 0.0 })
    }
    result.push(row)
  }
  let names : Array[String] = ["linear", "quadratic"]
  for i in 0..<knots.length() {
    names.push("hinge_\{i + 1}")
  }
  feature_set(result, names~)
}

///|
/// Computes a polynomial basis up to a validated degree.
pub fn power_basis(values : Array[Double], degree : Int) -> FeatureSet {
  let order = if degree > 0 { degree } else { 1 }
  let result : Array[Array[Double]] = Array::new(capacity=values.length())
  for value in values {
    let row = Array::new(capacity=order)
    let mut power = 1.0
    for _ in 1..<=order {
      power *= value
      row.push(power)
    }
    result.push(row)
  }
  let names : Array[String] = Array::new(capacity=order)
  for i in 0..<order {
    names.push("power_\{i + 1}")
  }
  feature_set(result, names~)
}

///|
/// Computes pairwise products for selected numeric features.
pub fn pairwise_interactions(
  matrix : Array[Array[Double]],
  include_original? : Bool = true,
) -> FeatureSet {
  if matrix.length() == 0 {
    return feature_set([])
  }
  let width = matrix[0].length()
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    let expanded : Array[Double] = Array::new(
      capacity=if include_original {
        width + width * (width - 1) / 2
      } else {
        width * (width - 1) / 2
      },
    )
    if include_original {
      for value in row {
        expanded.push(value)
      }
    }
    for i in 0..<width {
      for j in (i + 1)..<width {
        expanded.push(row[i] * row[j])
      }
    }
    result.push(expanded)
  }
  let names : Array[String] = Array::new()
  if include_original {
    for i in 0..<width {
      names.push("x\{i + 1}")
    }
  }
  for i in 0..<width {
    for j in (i + 1)..<width {
      names.push("x\{i + 1}_x_\{j + 1}")
    }
  }
  feature_set(result, names~)
}

///|
/// Applies a signed log1p transform.
pub fn signed_log_transform(values : Array[Double]) -> Array[Double] {
  let result = Array::new(capacity=values.length())
  for value in values {
    if value >= 0.0 {
      result.push(@math.ln(1.0 + value))
    } else {
      result.push(-@math.ln(1.0 + value.abs()))
    }
  }
  result
}

///|
/// Applies a signed square-root transform.
pub fn signed_sqrt_transform(values : Array[Double]) -> Array[Double] {
  let result = Array::new(capacity=values.length())
  for value in values {
    result.push(if value >= 0.0 { value.sqrt() } else { -value.abs().sqrt() })
  }
  result
}

///|
/// Fits a robust center-scale-power transformation specification.
pub fn transform_spec(
  values : Array[Double],
  power? : Double = 1.0,
  lower? : Double = -1.0e12,
  upper? : Double = 1.0e12,
  missing_value? : Double = 0.0,
) -> TransformSpec {
  let observed = values.filter(fn(value) { is_finite(value) })
  let center = mean_or(observed, 0.0)
  let scale = std_dev(observed).max(1.0e-12)
  { center, scale, power, lower, upper, missing_value }
}

///|
/// Applies a configured center-scale-power transformation.
pub fn apply_transform_spec(
  spec : TransformSpec,
  values : Array[Double],
) -> Array[Double] {
  let result = Array::new(capacity=values.length())
  for value in values {
    if !is_finite(value) {
      result.push(spec.missing_value)
    } else {
      let bounded = clamp(value, spec.lower, spec.upper)
      let standardized = (bounded - spec.center) / spec.scale
      result.push(
        if spec.power == 1.0 {
          standardized
        } else if standardized >= 0.0 {
          @math.pow(standardized, spec.power)
        } else {
          -@math.pow(standardized.abs(), spec.power)
        },
      )
    }
  }
  result
}

///|
/// Clips every matrix column to supplied quantile limits.
pub fn clip_matrix(
  matrix : Array[Array[Double]],
  lower_probability? : Double = 0.01,
  upper_probability? : Double = 0.99,
) -> Array[Array[Double]] {
  if matrix.length() == 0 {
    return []
  }
  let width = matrix[0].length()
  let lower = Array::new(capacity=width)
  let upper = Array::new(capacity=width)
  for j in 0..<width {
    let column = fe_column(matrix, j)
    lower.push(quantile(column, lower_probability))
    upper.push(quantile(column, upper_probability))
  }
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    let clipped = Array::new(capacity=width)
    for j in 0..<width {
      clipped.push(clamp(row[j], lower[j], upper[j]))
    }
    result.push(clipped)
  }
  result
}

///|
/// Computes robust median and MAD scales for each feature.
pub fn robust_feature_scales(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  if matrix.length() == 0 {
    return []
  }
  let width = matrix[0].length()
  let result : Array[Array[Double]] = Array::new(capacity=width)
  for j in 0..<width {
    let column = fe_column(matrix, j)
    let center = quantile(column, 0.5)
    let deviations = Array::new(capacity=column.length())
    for value in column {
      deviations.push((value - center).abs())
    }
    result.push([center, 1.4826 * quantile(deviations, 0.5).max(1.0e-12)])
  }
  result
}

///|
/// Applies robust median/MAD scaling to a feature matrix.
pub fn robust_scale_matrix(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  if matrix.length() == 0 {
    return []
  }
  let scales = robust_feature_scales(matrix)
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    let scaled = Array::new(capacity=row.length())
    for j in 0..<row.length() {
      scaled.push((row[j] - scales[j][0]) / scales[j][1])
    }
    result.push(scaled)
  }
  result
}

///|
/// Selects features by variance and absolute outcome correlation.
pub fn screen_features(
  matrix : Array[Array[Double]],
  outcome : Array[Double],
  minimum_variance? : Double = 1.0e-8,
  maximum_correlation? : Double = 0.995,
) -> FeatureScreening {
  let width = if matrix.length() == 0 { 0 } else { matrix[0].length() }
  let variances = feature_variances(matrix)
  let correlations = Array::new(capacity=width)
  let selected : Array[Int] = Array::new()
  let dropped : Array[Int] = Array::new()
  let reasons : Array[String] = Array::new()
  for j in 0..<width {
    let correlation_value = correlation(fe_column(matrix, j), outcome)
    correlations.push(correlation_value)
    if variances[j] > minimum_variance &&
      correlation_value.abs() < maximum_correlation {
      selected.push(j)
    } else {
      dropped.push(j)
      reasons.push(
        if variances[j] <= minimum_variance {
          "low variance"
        } else {
          "near-perfect outcome correlation"
        },
      )
    }
  }
  {
    selected,
    dropped,
    variance: variances,
    outcome_correlation: correlations,
    reasons,
  }
}

///|
/// Computes a variance inflation proxy from pairwise correlations.
pub fn variance_inflation_proxy(matrix : Array[Array[Double]]) -> Array[Double] {
  if matrix.length() == 0 {
    return []
  }
  let width = matrix[0].length()
  let result = Array::make(width, 1.0)
  for i in 0..<width {
    let mut total = 0.0
    for j in 0..<width {
      if i != j {
        let value = correlation(fe_column(matrix, i), fe_column(matrix, j))
        total += value * value
      }
    }
    result[i] = 1.0 / (1.0 - total.min(0.999999))
  }
  result
}

///|
/// Returns columns below a VIF threshold.
pub fn select_low_collinearity(
  matrix : Array[Array[Double]],
  threshold? : Double = 10.0,
) -> Array[Int] {
  let vif = variance_inflation_proxy(matrix)
  let result : Array[Int] = Array::new()
  for i in 0..<vif.length() {
    if vif[i] <= threshold {
      result.push(i)
    }
  }
  result
}

///|
/// Splits rows into train and validation sets using a deterministic seed.
pub fn feature_train_validation_split(
  matrix : Array[Array[Double]],
  validation_fraction? : Double = 0.2,
  seed? : UInt64 = 20260819,
) -> FeatureSplit {
  let n = matrix.length()
  let validation_size = (n.to_double() * clamp(validation_fraction, 0.0, 0.9))
    .round()
    .to_int()
  let order = shuffled_indices(n, seed)
  let validation_indices = order[:validation_size].to_owned()
  let train_indices = order[validation_size:].to_owned()
  let train_values = Array::new(capacity=train_indices.length())
  let validation_values = Array::new(capacity=validation_indices.length())
  for index in train_indices {
    train_values.push(matrix[index])
  }
  for index in validation_indices {
    validation_values.push(matrix[index])
  }
  { train_indices, validation_indices, train_values, validation_values }
}

///|
/// Creates deterministic feature names for a matrix width.
pub fn generated_feature_names(
  width : Int,
  prefix? : String = "x",
) -> Array[String] {
  let result : Array[String] = Array::new(
    capacity=if width > 0 { width } else { 0 },
  )
  for i in 0..<width {
    result.push("\{prefix}\{i + 1}")
  }
  result
}

///|
/// Aligns a matrix to a requested set of columns.
pub fn select_feature_columns(
  matrix : Array[Array[Double]],
  indexes : Array[Int],
) -> FeatureSet {
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    let selected = Array::new(capacity=indexes.length())
    for index in indexes {
      if index >= 0 && index < row.length() {
        selected.push(row[index])
      }
    }
    result.push(selected)
  }
  let names = generated_feature_names(indexes.length())
  feature_set(result, names~)
}

///|
/// Concatenates two feature matrices column-wise.
pub fn concatenate_features(
  first : Array[Array[Double]],
  second : Array[Array[Double]],
) -> FeatureSet {
  let n = first.length().min(second.length())
  let result : Array[Array[Double]] = Array::new(capacity=n)
  for i in 0..<n {
    let row = first[i].copy()
    for value in second[i] {
      row.push(value)
    }
    result.push(row)
  }
  let width = if n == 0 { 0 } else { result[0].length() }
  feature_set(result, names=generated_feature_names(width))
}

///|
/// Centers a matrix using externally supplied feature means.
pub fn center_with_means(
  matrix : Array[Array[Double]],
  means : Array[Double],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    let centered = Array::new(capacity=row.length())
    for j in 0..<row.length() {
      let mean_value = if j < means.length() { means[j] } else { 0.0 }
      centered.push(row[j] - mean_value)
    }
    result.push(centered)
  }
  result
}

///|
/// Scales a matrix using externally supplied feature scales.
pub fn scale_with_scales(
  matrix : Array[Array[Double]],
  scales : Array[Double],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    let scaled = Array::new(capacity=row.length())
    for j in 0..<row.length() {
      let scale = if j < scales.length() {
        scales[j].abs().max(1.0e-12)
      } else {
        1.0
      }
      scaled.push(row[j] / scale)
    }
    result.push(scaled)
  }
  result
}

///|
/// Generates a deterministic hash feature without external dictionaries.
pub fn integer_hash_feature(
  values : Array[Int],
  buckets : Int,
) -> Array[Double] {
  let width = if buckets > 0 { buckets } else { 1 }
  let result = Array::make(width, 0.0)
  for value in values {
    let index = value.abs() % width
    result[index] += 1.0
  }
  result
}

///|
/// Computes row norms for regularization diagnostics.
pub fn feature_row_norms(matrix : Array[Array[Double]]) -> Array[Double] {
  let result = Array::new(capacity=matrix.length())
  for row in matrix {
    let mut norm = 0.0
    for value in row {
      norm += value * value
    }
    result.push(norm.sqrt())
  }
  result
}

///|
/// Computes column L2 norms.
pub fn feature_column_norms(matrix : Array[Array[Double]]) -> Array[Double] {
  if matrix.length() == 0 {
    return []
  }
  let width = matrix[0].length()
  let result = Array::make(width, 0.0)
  for row in matrix {
    for j in 0..<width {
      result[j] += row[j] * row[j]
    }
  }
  for j in 0..<width {
    result[j] = result[j].sqrt()
  }
  result
}

///|
/// Normalizes columns to unit L2 norm.
pub fn normalize_feature_columns(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let norms = feature_column_norms(matrix)
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    let normalized = Array::new(capacity=row.length())
    for j in 0..<row.length() {
      normalized.push(if norms[j] == 0.0 { 0.0 } else { row[j] / norms[j] })
    }
    result.push(normalized)
  }
  result
}

///|
/// Computes the sparse density of a feature matrix.
pub fn feature_density(
  matrix : Array[Array[Double]],
  tolerance? : Double = 1.0e-12,
) -> Double {
  let mut nonzero = 0
  let mut total = 0
  for row in matrix {
    for value in row {
      total += 1
      if value.abs() > tolerance {
        nonzero += 1
      }
    }
  }
  if total == 0 {
    0.0
  } else {
    nonzero.to_double() / total.to_double()
  }
}

///|
/// Returns a compact feature engineering summary vector.
pub fn feature_summary(features : FeatureSet) -> Array[Double] {
  [
    features.row_count.to_double(),
    features.column_count.to_double(),
    feature_density(features.values),
    features.fingerprint.to_double(),
  ]
}
