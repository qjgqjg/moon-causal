///|
/// Configuration for deterministic numeric feature transformations.
pub struct CausalTransformSpec {
  center : Double
  scale : Double
  lower : Double
  upper : Double
  log_shift : Double
  missing_value : Double
  fitted : Bool
}

///|
/// Diagnostics returned with a transformed vector.
pub struct CausalTransformResult {
  values : Array[Double]
  missing : Int
  clipped : Int
  mean : Double
  scale : Double
  minimum : Double
  maximum : Double
  passes : Bool
}

///|
/// Matrix-level transformation diagnostics.
pub struct CausalMatrixTransformSummary {
  rows : Int
  columns : Int
  missing_before : Int
  missing_after : Int
  clipped : Int
  constant_columns : Int
  rank : Int
  passes : Bool
}

///|
/// Creates a bounded transform specification.
pub fn causal_transform_spec(
  center? : Double = 0.0,
  scale? : Double = 1.0,
  lower? : Double = -1.0e300,
  upper? : Double = 1.0e300,
  log_shift? : Double = 1.0,
  missing_value? : Double = 0.0,
  fitted? : Bool = false,
) -> CausalTransformSpec {
  let safe_scale = if is_finite(scale) && scale.abs() > 1.0e-12 {
    scale.abs()
  } else {
    1.0
  }
  {
    center: if is_finite(center) {
      center
    } else {
      0.0
    },
    scale: safe_scale,
    lower: lower.min(upper),
    upper: upper.max(lower),
    log_shift: log_shift.max(1.0e-12),
    missing_value: if is_finite(missing_value) {
      missing_value
    } else {
      0.0
    },
    fitted,
  }
}

///|
/// Returns the finite values in a vector.
pub fn causal_finite_values(values : Array[Double]) -> Array[Double] {
  let result : Array[Double] = Array::new(capacity=values.length())
  for value in values {
    if is_finite(value) {
      result.push(value)
    }
  }
  result
}

///|
/// Fits location, scale, and winsorization limits from a reference vector.
pub fn fit_causal_transform(
  values : Array[Double],
  lower_quantile? : Double = 0.01,
  upper_quantile? : Double = 0.99,
  missing_value? : Double = 0.0,
) -> CausalTransformSpec {
  let finite = causal_finite_values(values)
  let center = mean_or(finite, 0.0)
  let deviation = if finite.length() > 1 { std_dev(finite) } else { 1.0 }
  let lower = if finite.length() == 0 {
    -1.0e300
  } else {
    quantile(finite, clamp(lower_quantile, 0.0, 1.0))
  }
  let upper = if finite.length() == 0 {
    1.0e300
  } else {
    quantile(finite, clamp(upper_quantile, 0.0, 1.0))
  }
  causal_transform_spec(
    center~,
    scale=deviation,
    lower~,
    upper~,
    missing_value~,
    fitted=true,
  )
}

///|
/// Returns the minimum of finite values with a fallback.
pub fn causal_minimum(
  values : Array[Double],
  fallback? : Double = 0.0,
) -> Double {
  let finite = causal_finite_values(values)
  if finite.length() == 0 {
    fallback
  } else {
    let mut result = finite[0]
    for value in finite[1:] {
      if value < result {
        result = value
      }
    }
    result
  }
}

///|
/// Returns the maximum of finite values with a fallback.
pub fn causal_maximum(
  values : Array[Double],
  fallback? : Double = 0.0,
) -> Double {
  let finite = causal_finite_values(values)
  if finite.length() == 0 {
    fallback
  } else {
    let mut result = finite[0]
    for value in finite[1:] {
      if value > result {
        result = value
      }
    }
    result
  }
}

///|
/// Applies clipping and standardization using a fitted specification.
pub fn apply_causal_transform(
  values : Array[Double],
  spec : CausalTransformSpec,
) -> CausalTransformResult {
  let result : Array[Double] = Array::new(capacity=values.length())
  let mut missing = 0
  let mut clipped = 0
  for value in values {
    let present = is_finite(value)
    let source = if present {
      value
    } else {
      missing += 1
      spec.missing_value
    }
    let bounded = clamp(source, spec.lower, spec.upper)
    if bounded != source {
      clipped += 1
    }
    result.push((bounded - spec.center) / spec.scale)
  }
  let minimum = causal_minimum(result)
  let maximum = causal_maximum(result)
  {
    values: result,
    missing,
    clipped,
    mean: mean_or(result, 0.0),
    scale: if result.length() > 1 {
      std_dev(result)
    } else {
      0.0
    },
    minimum,
    maximum,
    passes: finite_array(result),
  }
}

///|
/// Applies a fitted specification without clipping the reference limits.
pub fn apply_causal_standardization(
  values : Array[Double],
  spec : CausalTransformSpec,
) -> Array[Double] {
  let result : Array[Double] = Array::new(capacity=values.length())
  for value in values {
    let source = if is_finite(value) { value } else { spec.missing_value }
    result.push((source - spec.center) / spec.scale)
  }
  result
}

///|
/// Reverses a location-scale transformation.
pub fn invert_causal_standardization(
  values : Array[Double],
  spec : CausalTransformSpec,
) -> Array[Double] {
  let result : Array[Double] = Array::new(capacity=values.length())
  for value in values {
    result.push(value * spec.scale + spec.center)
  }
  result
}

///|
/// Standardizes a vector using its own finite observations.
pub fn causal_standardize(values : Array[Double]) -> CausalTransformResult {
  apply_causal_transform(values, fit_causal_transform(values))
}

///|
/// Clips a vector to explicit limits and records boundary hits.
pub fn causal_clip(
  values : Array[Double],
  lower : Double,
  upper : Double,
) -> CausalTransformResult {
  let spec = causal_transform_spec(lower~, upper~)
  let result = apply_causal_transform(values, spec)
  {
    values: result.values,
    missing: result.missing,
    clipped: result.clipped,
    mean: result.mean,
    scale: result.scale,
    minimum: result.minimum,
    maximum: result.maximum,
    passes: result.passes && lower <= upper,
  }
}

///|
/// Applies a numerically safe log1p transform.
pub fn causal_log1p(
  values : Array[Double],
  shift? : Double = 1.0,
) -> CausalTransformResult {
  let offset = shift.max(1.0e-12)
  let result : Array[Double] = Array::new(capacity=values.length())
  let mut missing = 0
  let clipped = 0
  for value in values {
    if !is_finite(value) || value + offset <= 0.0 {
      result.push(0.0)
      missing += 1
    } else {
      result.push(@math.ln(1.0 + value + offset))
    }
  }
  {
    values: result,
    missing,
    clipped,
    mean: mean_or(result, 0.0),
    scale: if result.length() > 1 {
      std_dev(result)
    } else {
      0.0
    },
    minimum: causal_minimum(result),
    maximum: causal_maximum(result),
    passes: finite_array(result),
  }
}

///|
/// Applies a safe logit transform to probabilities.
pub fn causal_logit_transform(
  probabilities : Array[Double],
  epsilon? : Double = 1.0e-6,
) -> CausalTransformResult {
  let result : Array[Double] = Array::new(capacity=probabilities.length())
  let mut clipped = 0
  for probability in probabilities {
    let safe = safe_probability(probability, epsilon~)
    if safe != probability {
      clipped += 1
    }
    result.push(@math.ln(safe / (1.0 - safe)))
  }
  {
    values: result,
    missing: 0,
    clipped,
    mean: mean_or(result, 0.0),
    scale: if result.length() > 1 {
      std_dev(result)
    } else {
      0.0
    },
    minimum: causal_minimum(result),
    maximum: causal_maximum(result),
    passes: finite_array(result),
  }
}

///|
/// Converts a score vector to average ranks with stable ties.
pub fn causal_average_ranks(values : Array[Double]) -> Array[Double] {
  let result = Array::make(values.length(), 0.0)
  for i in 0..<values.length() {
    if !is_finite(values[i]) {
      result[i] = 0.0
    } else {
      let lower = values.fold(init=0, fn(total, value) {
        if is_finite(value) && value < values[i] {
          total + 1
        } else {
          total
        }
      })
      let ties = values.fold(init=0, fn(total, value) {
        if is_finite(value) && value == values[i] {
          total + 1
        } else {
          total
        }
      })
      result[i] = lower.to_double() + (ties.to_double() + 1.0) / 2.0
    }
  }
  result
}

///|
/// Converts average ranks to percentile scores in [0, 1].
pub fn causal_percentile_ranks(values : Array[Double]) -> Array[Double] {
  let ranks = causal_average_ranks(values)
  let denominator = values.length().to_double().max(1.0)
  ranks.map(fn(rank) { (rank - 1.0) / (denominator - 1.0).max(1.0) })
}

///|
/// Encodes a binary treatment vector as a numeric column.
pub fn causal_treatment_code(treatment : Array[Bool]) -> Array[Double] {
  treatment.map(fn(value) { if value { 1.0 } else { 0.0 } })
}

///|
/// Encodes a numeric threshold into a binary treatment vector.
pub fn causal_threshold_treatment(
  scores : Array[Double],
  threshold : Double,
) -> Array[Bool] {
  scores.map(fn(value) { is_finite(value) && value >= threshold })
}

///|
/// Computes a centered outcome vector and returns its mean.
pub fn causal_center_outcome(outcome : Array[Double]) -> CausalTransformResult {
  let spec = fit_causal_transform(
    outcome,
    lower_quantile=0.0,
    upper_quantile=1.0,
  )
  let result = apply_causal_standardization(
    outcome,
    causal_transform_spec(center=spec.center),
  )
  {
    values: result,
    missing: missing_count(outcome),
    clipped: 0,
    mean: spec.center,
    scale: spec.scale,
    minimum: causal_minimum(result),
    maximum: causal_maximum(result),
    passes: finite_array(result),
  }
}

///|
/// Returns a matrix transpose, preserving short rows with a fill value.
pub fn causal_transpose(
  matrix : Array[Array[Double]],
  columns? : Int = -1,
  fill? : Double = 0.0,
) -> Array[Array[Double]] {
  let width = if columns >= 0 {
    columns
  } else {
    let mut maximum = 0
    for row in matrix {
      if row.length() > maximum {
        maximum = row.length()
      }
    }
    maximum
  }
  let result : Array[Array[Double]] = Array::new(capacity=width)
  for column in 0..<width {
    let values : Array[Double] = Array::new(capacity=matrix.length())
    for row in matrix {
      values.push(if column < row.length() { row[column] } else { fill })
    }
    result.push(values)
  }
  result
}

///|
/// Reconstructs a rectangular matrix from column vectors.
pub fn causal_from_columns(
  columns : Array[Array[Double]],
  rows? : Int = -1,
  fill? : Double = 0.0,
) -> Array[Array[Double]] {
  let height = if rows >= 0 {
    rows
  } else {
    let mut maximum = 0
    for column in columns {
      if column.length() > maximum {
        maximum = column.length()
      }
    }
    maximum
  }
  let result : Array[Array[Double]] = Array::new(capacity=height)
  for row_index in 0..<height {
    let row : Array[Double] = Array::new(capacity=columns.length())
    for column in columns {
      row.push(
        if row_index < column.length() {
          column[row_index]
        } else {
          fill
        },
      )
    }
    result.push(row)
  }
  result
}

///|
/// Selects columns by index, filling invalid indices with zeros.
pub fn causal_select_columns(
  matrix : Array[Array[Double]],
  indices : Array[Int],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    let selected : Array[Double] = Array::new(capacity=indices.length())
    for index in indices {
      selected.push(
        if index >= 0 && index < row.length() {
          row[index]
        } else {
          0.0
        },
      )
    }
    result.push(selected)
  }
  result
}

///|
/// Adds an intercept column to a design matrix.
pub fn causal_add_intercept(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    let next : Array[Double] = Array::new(capacity=row.length() + 1)
    next.push(1.0)
    for value in row {
      next.push(value)
    }
    result.push(next)
  }
  result
}

///|
/// Adds one pairwise interaction to a matrix.
pub fn causal_add_interaction(
  matrix : Array[Array[Double]],
  first : Int,
  second : Int,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    let next = row.copy()
    let left = if first >= 0 && first < row.length() { row[first] } else { 0.0 }
    let right = if second >= 0 && second < row.length() {
      row[second]
    } else {
      0.0
    }
    next.push(left * right)
    result.push(next)
  }
  result
}

///|
/// Adds polynomial powers of one column to a matrix.
pub fn causal_add_powers(
  matrix : Array[Array[Double]],
  column : Int,
  degree : Int,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  let highest = degree.max(0)
  for row in matrix {
    let next = row.copy()
    let base = if column >= 0 && column < row.length() {
      row[column]
    } else {
      0.0
    }
    let mut power = base
    for _ in 2..<(highest + 1) {
      power *= base
      next.push(power)
    }
    result.push(next)
  }
  result
}

///|
/// Standardizes every matrix column independently.
pub fn causal_standardize_matrix(
  matrix : Array[Array[Double]],
) -> (Array[Array[Double]], Array[CausalTransformSpec]) {
  let columns = causal_transpose(matrix)
  let transformed : Array[Array[Double]] = Array::new(capacity=columns.length())
  let specs : Array[CausalTransformSpec] = Array::new(capacity=columns.length())
  for column in columns {
    let spec = fit_causal_transform(column)
    specs.push(spec)
    transformed.push(apply_causal_standardization(column, spec))
  }
  (causal_from_columns(transformed, rows=matrix.length()), specs)
}

///|
/// Clips every matrix column to its own explicit bounds.
pub fn causal_clip_matrix(
  matrix : Array[Array[Double]],
  lower : Double,
  upper : Double,
) -> (Array[Array[Double]], Int) {
  let columns = causal_transpose(matrix)
  let transformed : Array[Array[Double]] = Array::new(capacity=columns.length())
  let mut clipped = 0
  for column in columns {
    let result = causal_clip(column, lower, upper)
    clipped += result.clipped
    transformed.push(result.values)
  }
  (causal_from_columns(transformed, rows=matrix.length()), clipped)
}

///|
/// Counts non-finite matrix entries.
pub fn causal_matrix_missing(matrix : Array[Array[Double]]) -> Int {
  let mut result = 0
  for row in matrix {
    for value in row {
      if !is_finite(value) {
        result += 1
      }
    }
  }
  result
}

///|
/// Counts columns with no finite variation.
pub fn causal_constant_columns(matrix : Array[Array[Double]]) -> Int {
  let mut result = 0
  for column in causal_transpose(matrix) {
    if !has_variation(causal_finite_values(column)) {
      result += 1
    }
  }
  result
}

///|
/// Produces matrix transform diagnostics after standardization.
pub fn causal_matrix_transform_summary(
  before : Array[Array[Double]],
  after : Array[Array[Double]],
  clipped : Int,
) -> CausalMatrixTransformSummary {
  let rows = before.length()
  let columns = if before.length() == 0 { 0 } else { before[0].length() }
  let missing_before = causal_matrix_missing(before)
  let missing_after = causal_matrix_missing(after)
  let constants = causal_constant_columns(before)
  let rank_value = causal_matrix_rank_proxy(after)
  {
    rows,
    columns,
    missing_before,
    missing_after,
    clipped,
    constant_columns: constants,
    rank: rank_value,
    passes: missing_after == 0 &&
    after.length() == rows &&
    (after.length() == 0 || after[0].length() == columns),
  }
}

///|
/// Estimates matrix rank by counting columns with finite variation.
pub fn causal_matrix_rank_proxy(matrix : Array[Array[Double]]) -> Int {
  let mut rank = 0
  for column in causal_transpose(matrix) {
    let finite = causal_finite_values(column)
    if finite.length() > 0 && has_variation(finite) {
      rank += 1
    }
  }
  rank
}

///|
/// Builds a treatment-design matrix with optional intercept and interactions.
pub fn causal_design_matrix(
  dataset : CausalDataset,
  include_intercept? : Bool = true,
  interaction_pairs? : Array[(Int, Int)] = [],
) -> Array[Array[Double]] {
  let base = if include_intercept {
    causal_add_intercept(dataset.covariates)
  } else {
    dataset.covariates
  }
  let mut result = base
  for pair in interaction_pairs {
    result = causal_add_interaction(result, pair.0, pair.1)
  }
  result
}

///|
/// Returns a dataset with transformed covariates and preserved treatment/outcome.
pub fn causal_dataset_transform(
  dataset : CausalDataset,
  specs : Array[CausalTransformSpec],
) -> CausalDataset {
  let columns = causal_transpose(dataset.covariates)
  let transformed : Array[Array[Double]] = Array::new(capacity=columns.length())
  for i in 0..<columns.length() {
    let spec = if i < specs.length() {
      specs[i]
    } else {
      fit_causal_transform(columns[i])
    }
    transformed.push(apply_causal_standardization(columns[i], spec))
  }
  CausalDataset::new(
    causal_from_columns(transformed, rows=dataset.n()),
    dataset.treatment,
    dataset.outcome,
  )
}

///|
/// Returns a dataset with a centered or transformed outcome.
pub fn causal_dataset_outcome_transform(
  dataset : CausalDataset,
  transformed_outcome : Array[Double],
) -> CausalDataset {
  let outcome : Array[Double] = Array::new(capacity=dataset.n())
  for i in 0..<dataset.n() {
    outcome.push(
      if i < transformed_outcome.length() {
        transformed_outcome[i]
      } else {
        0.0
      },
    )
  }
  CausalDataset::new(dataset.covariates, dataset.treatment, outcome)
}

///|
/// Selects observations above a finite score threshold.
pub fn causal_dataset_filter_score(
  dataset : CausalDataset,
  scores : Array[Double],
  lower : Double,
  upper : Double,
) -> CausalDataset {
  let indices : Array[Int] = Array::new()
  let n = dataset.n().min(scores.length())
  for i in 0..<n {
    if is_finite(scores[i]) && scores[i] >= lower && scores[i] <= upper {
      indices.push(i)
    }
  }
  dataset.select(indices)
}

///|
/// Returns a causal dataset summary vector for audit logs.
pub fn causal_transform_dataset_summary(
  dataset : CausalDataset,
) -> Array[Double] {
  let matrix_summary = causal_matrix_transform_summary(
    dataset.covariates,
    dataset.covariates,
    0,
  )
  [
    dataset.n().to_double(),
    dataset.p().to_double(),
    dataset.treated_count().to_double(),
    dataset.control_count().to_double(),
    matrix_summary.missing_before.to_double(),
    matrix_summary.constant_columns.to_double(),
    mean_or(dataset.outcome, 0.0),
    std_dev(dataset.outcome),
  ]
}

///|
/// Computes a deterministic fingerprint for a transform specification.
pub fn causal_transform_fingerprint(spec : CausalTransformSpec) -> UInt64 {
  matrix_checksum([
    [spec.center, spec.scale, spec.lower, spec.upper, spec.log_shift],
  ])
}

///|
/// Computes a combined fingerprint for a list of specifications.
pub fn causal_transform_specs_fingerprint(
  specs : Array[CausalTransformSpec],
) -> UInt64 {
  let rows : Array[Array[Double]] = Array::new(capacity=specs.length())
  for spec in specs {
    rows.push([
      spec.center,
      spec.scale,
      spec.lower,
      spec.upper,
      spec.log_shift,
      if spec.fitted {
        1.0
      } else {
        0.0
      },
    ])
  }
  matrix_checksum(rows)
}
