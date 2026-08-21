///|
/// Matrix-level numerical diagnostics.
pub struct MatrixDiagnostics {
  rows : Int
  columns : Int
  rank_proxy : Int
  condition_proxy : Double
  density : Double
  finite : Bool
  fingerprint : UInt64
}

///|
/// Computes row sums.
pub fn advanced_matrix_row_sums(matrix : Array[Array[Double]]) -> Array[Double] {
  let result = Array::new(capacity=matrix.length())
  for row in matrix {
    result.push(sum(row))
  }
  result
}

///|
/// Computes row means.
pub fn advanced_matrix_row_means(
  matrix : Array[Array[Double]],
) -> Array[Double] {
  let result = Array::new(capacity=matrix.length())
  for row in matrix {
    result.push(mean_or(row, 0.0))
  }
  result
}

///|
/// Computes row variances.
pub fn advanced_matrix_row_variances(
  matrix : Array[Array[Double]],
) -> Array[Double] {
  let result = Array::new(capacity=matrix.length())
  for row in matrix {
    result.push(variance(row))
  }
  result
}

///|
/// Computes column sums.
pub fn advanced_matrix_column_sums(
  matrix : Array[Array[Double]],
) -> Array[Double] {
  if matrix.length() == 0 {
    return []
  }
  let width = matrix[0].length()
  let result = Array::make(width, 0.0)
  for row in matrix {
    for j in 0..<width {
      result[j] += row[j]
    }
  }
  result
}

///|
/// Computes column means.
pub fn advanced_matrix_column_means(
  matrix : Array[Array[Double]],
) -> Array[Double] {
  let sums = advanced_matrix_column_sums(matrix)
  let result = Array::new(capacity=sums.length())
  for value in sums {
    result.push(
      if matrix.length() == 0 {
        0.0
      } else {
        value / matrix.length().to_double()
      },
    )
  }
  result
}

///|
/// Computes column variances.
pub fn advanced_matrix_column_variances(
  matrix : Array[Array[Double]],
) -> Array[Double] {
  if matrix.length() == 0 {
    return []
  }
  let means = advanced_matrix_column_means(matrix)
  let result = Array::make(means.length(), 0.0)
  for row in matrix {
    for j in 0..<means.length() {
      let difference = row[j] - means[j]
      result[j] += difference * difference
    }
  }
  for j in 0..<result.length() {
    result[j] /= matrix.length().to_double()
  }
  result
}

///|
/// Computes a covariance matrix.
pub fn advanced_covariance_matrix(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  if matrix.length() == 0 {
    return []
  }
  let means = advanced_matrix_column_means(matrix)
  let width = means.length()
  let result = Array::make(width, Array::make(width, 0.0))
  for row in matrix {
    for j in 0..<width {
      for k in 0..<width {
        result[j][k] += (row[j] - means[j]) * (row[k] - means[k])
      }
    }
  }
  for j in 0..<width {
    for k in 0..<width {
      result[j][k] /= matrix.length().to_double().max(1.0)
    }
  }
  result
}

///|
/// Computes a correlation matrix.
pub fn advanced_correlation_matrix(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let covariance_matrix = advanced_covariance_matrix(matrix)
  let variances = advanced_matrix_column_variances(matrix)
  let result = Array::make(
    covariance_matrix.length(),
    Array::make(covariance_matrix.length(), 0.0),
  )
  for j in 0..<covariance_matrix.length() {
    for k in 0..<covariance_matrix.length() {
      let scale = (variances[j] * variances[k]).sqrt()
      result[j][k] = if scale == 0.0 {
        0.0
      } else {
        covariance_matrix[j][k] / scale
      }
    }
  }
  result
}

///|
/// Computes Euclidean distance between two vectors.
pub fn advanced_euclidean_distance(
  first : Array[Double],
  second : Array[Double],
) -> Double {
  let n = first.length().min(second.length())
  let mut total = 0.0
  for i in 0..<n {
    let difference = first[i] - second[i]
    total += difference * difference
  }
  total.sqrt()
}

///|
/// Computes Manhattan distance between two vectors.
pub fn advanced_manhattan_distance(
  first : Array[Double],
  second : Array[Double],
) -> Double {
  let n = first.length().min(second.length())
  let mut total = 0.0
  for i in 0..<n {
    total += (first[i] - second[i]).abs()
  }
  total
}

///|
/// Computes a pairwise Euclidean distance matrix.
pub fn advanced_distance_matrix(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result = Array::make(matrix.length(), Array::make(matrix.length(), 0.0))
  for i in 0..<matrix.length() {
    for j in (i + 1)..<matrix.length() {
      let distance = advanced_euclidean_distance(matrix[i], matrix[j])
      result[i][j] = distance
      result[j][i] = distance
    }
  }
  result
}

///|
/// Returns a row-normalized matrix.
pub fn advanced_row_normalize(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    let total = sum(row)
    let normalized = Array::new(capacity=row.length())
    for value in row {
      normalized.push(if total == 0.0 { 0.0 } else { value / total })
    }
    result.push(normalized)
  }
  result
}

///|
/// Returns a column-normalized matrix.
pub fn advanced_column_normalize(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let sums = advanced_matrix_column_sums(matrix)
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    let normalized = Array::new(capacity=row.length())
    for j in 0..<row.length() {
      normalized.push(if sums[j] == 0.0 { 0.0 } else { row[j] / sums[j] })
    }
    result.push(normalized)
  }
  result
}

///|
/// Computes an element-wise absolute matrix.
pub fn advanced_matrix_abs(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    result.push(row.map(fn(value) { value.abs() }))
  }
  result
}

///|
/// Computes an element-wise square matrix.
pub fn advanced_matrix_square(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    result.push(row.map(fn(value) { value * value }))
  }
  result
}

///|
/// Computes a matrix trace.
pub fn advanced_matrix_trace(matrix : Array[Array[Double]]) -> Double {
  let mut result = 0.0
  for i in 0..<matrix.length() {
    if i < matrix[i].length() {
      result += matrix[i][i]
    }
  }
  result
}

///|
/// Computes a diagonal vector.
pub fn advanced_matrix_diagonal(matrix : Array[Array[Double]]) -> Array[Double] {
  let result = Array::new(capacity=matrix.length())
  for i in 0..<matrix.length() {
    if i < matrix[i].length() {
      result.push(matrix[i][i])
    }
  }
  result
}

///|
/// Computes a Frobenius norm.
pub fn advanced_matrix_frobenius(matrix : Array[Array[Double]]) -> Double {
  let mut result = 0.0
  for row in matrix {
    for value in row {
      result += value * value
    }
  }
  result.sqrt()
}

///|
/// Computes the maximum absolute matrix entry.
pub fn advanced_matrix_max_abs(matrix : Array[Array[Double]]) -> Double {
  let mut result = 0.0
  for row in matrix {
    for value in row {
      if value.abs() > result {
        result = value.abs()
      }
    }
  }
  result
}

///|
/// Computes a rank proxy from diagonal energy.
pub fn advanced_matrix_rank_proxy(
  matrix : Array[Array[Double]],
  tolerance? : Double = 1.0e-8,
) -> Int {
  let diagonal = advanced_matrix_diagonal(matrix)
  let mut result = 0
  for value in diagonal {
    if value.abs() > tolerance {
      result += 1
    }
  }
  result
}

///|
/// Computes a condition proxy from diagonal extrema.
pub fn advanced_matrix_condition_proxy(matrix : Array[Array[Double]]) -> Double {
  let diagonal = advanced_matrix_diagonal(matrix)
  let mut minimum = 1.0e300
  let mut maximum = 0.0
  for value in diagonal {
    if value.abs() > maximum {
      maximum = value.abs()
    }
    if value.abs() > 0.0 && value.abs() < minimum {
      minimum = value.abs()
    }
  }
  if minimum == 1.0e300 {
    0.0
  } else {
    maximum / minimum
  }
}

///|
/// Audits a numeric matrix.
pub fn advanced_matrix_diagnostics(
  matrix : Array[Array[Double]],
) -> MatrixDiagnostics {
  let columns = if matrix.length() == 0 { 0 } else { matrix[0].length() }
  let mut finite = true
  let mut nonzero = 0
  let mut total = 0
  for row in matrix {
    for value in row {
      if !is_finite(value) {
        finite = false
      }
      if value.abs() > 1.0e-12 {
        nonzero += 1
      }
      total += 1
    }
  }
  {
    rows: matrix.length(),
    columns,
    rank_proxy: advanced_matrix_rank_proxy(matrix),
    condition_proxy: advanced_matrix_condition_proxy(matrix),
    density: if total == 0 {
      0.0
    } else {
      nonzero.to_double() / total.to_double()
    },
    finite,
    fingerprint: matrix_checksum(matrix),
  }
}

///|
/// Returns a compact matrix diagnostic vector.
pub fn advanced_matrix_summary(
  diagnostics : MatrixDiagnostics,
) -> Array[Double] {
  [
    diagnostics.rows.to_double(),
    diagnostics.columns.to_double(),
    diagnostics.rank_proxy.to_double(),
    diagnostics.condition_proxy,
    diagnostics.density,
    if diagnostics.finite {
      1.0
    } else {
      0.0
    },
    diagnostics.fingerprint.to_double(),
  ]
}
