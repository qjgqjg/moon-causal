///|
pub fn dot(a : Array[Double], b : Array[Double]) -> Double {
  if a.length() != b.length() {
    return 0.0
  }
  let mut result = 0.0
  for i in 0..<a.length() {
    result += a[i] * b[i]
  }
  result
}

///|
pub fn add_vectors(a : Array[Double], b : Array[Double]) -> Array[Double] {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let result = Array::new(capacity=n)
  for i in 0..<n {
    result.push(a[i] + b[i])
  }
  result
}

///|
pub fn subtract_vectors(a : Array[Double], b : Array[Double]) -> Array[Double] {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let result = Array::new(capacity=n)
  for i in 0..<n {
    result.push(a[i] - b[i])
  }
  result
}

///|
pub fn scale_vector(a : Array[Double], factor : Double) -> Array[Double] {
  let result = Array::new(capacity=a.length())
  for value in a {
    result.push(value * factor)
  }
  result
}

///|
pub fn l2_norm(a : Array[Double]) -> Double {
  dot(a, a).sqrt()
}

///|
pub fn l1_norm(a : Array[Double]) -> Double {
  let mut result = 0.0
  for value in a {
    result += value.abs()
  }
  result
}

///|
pub fn matrix_zeros(rows : Int, columns : Int) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(
    capacity=if rows > 0 { rows } else { 0 },
  )
  for _ in 0..<rows {
    result.push(Array::make(columns, 0.0))
  }
  result
}

///|
pub fn identity_matrix(size : Int) -> Array[Array[Double]] {
  let result = matrix_zeros(size, size)
  for i in 0..<size {
    result[i][i] = 1.0
  }
  result
}

///|
pub fn matrix_transpose(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  if matrix.length() == 0 {
    return []
  }
  let columns = matrix[0].length()
  let result = matrix_zeros(columns, matrix.length())
  for i in 0..<matrix.length() {
    for j in 0..<columns {
      result[j][i] = matrix[i][j]
    }
  }
  result
}

///|
pub fn matrix_multiply(
  left : Array[Array[Double]],
  right : Array[Array[Double]],
) -> Array[Array[Double]] {
  if left.length() == 0 ||
    right.length() == 0 ||
    left[0].length() != right.length() {
    return []
  }
  let result = matrix_zeros(left.length(), right[0].length())
  for i in 0..<left.length() {
    for k in 0..<right.length() {
      let left_value = left[i][k]
      for j in 0..<right[0].length() {
        result[i][j] += left_value * right[k][j]
      }
    }
  }
  result
}

///|
pub fn matrix_vector_multiply(
  matrix : Array[Array[Double]],
  vector : Array[Double],
) -> Array[Double] {
  let result = Array::new(capacity=matrix.length())
  for row in matrix {
    result.push(dot(row, vector))
  }
  result
}

///|
/// Solves a dense linear system using Gaussian elimination with partial pivoting.
pub fn solve_linear_system(
  matrix : Array[Array[Double]],
  right_hand_side : Array[Double],
) -> Array[Double] {
  let n = matrix.length()
  if n == 0 || right_hand_side.length() != n {
    return []
  }
  let augmented = matrix_zeros(n, n + 1)
  for i in 0..<n {
    if matrix[i].length() != n {
      return []
    }
    for j in 0..<n {
      augmented[i][j] = matrix[i][j]
    }
    augmented[i][n] = right_hand_side[i]
  }
  for pivot in 0..<n {
    let mut best = pivot
    for row in (pivot + 1)..<n {
      if augmented[row][pivot].abs() > augmented[best][pivot].abs() {
        best = row
      }
    }
    if augmented[best][pivot].abs() < 1.0e-12 {
      return []
    }
    if best != pivot {
      let temporary = augmented[pivot]
      augmented[pivot] = augmented[best]
      augmented[best] = temporary
    }
    let divisor = augmented[pivot][pivot]
    for column in pivot..<(n + 1) {
      augmented[pivot][column] /= divisor
    }
    for row in 0..<n {
      if row != pivot {
        let factor = augmented[row][pivot]
        if factor != 0.0 {
          for column in pivot..<(n + 1) {
            augmented[row][column] -= factor * augmented[pivot][column]
          }
        }
      }
    }
  }
  let result = Array::new(capacity=n)
  for row in augmented {
    result.push(row[n])
  }
  result
}

///|
pub fn matrix_inverse(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = matrix.length()
  if n == 0 {
    return []
  }
  let result = matrix_zeros(n, n)
  for column in 0..<n {
    let rhs = Array::make(n, 0.0)
    rhs[column] = 1.0
    let solution = solve_linear_system(matrix, rhs)
    if solution.length() != n {
      return []
    }
    for row in 0..<n {
      result[row][column] = solution[row]
    }
  }
  result
}

///|
pub fn add_ridge(
  matrix : Array[Array[Double]],
  ridge : Double,
) -> Array[Array[Double]] {
  let result = matrix.copy()
  for i in 0..<result.length() {
    if i < result[i].length() {
      result[i][i] += ridge
    }
  }
  result
}

///|
pub fn covariance_matrix(
  rows : Array[Array[Double]],
  sample? : Bool = true,
) -> Array[Array[Double]] {
  if rows.length() == 0 {
    return []
  }
  let columns = rows[0].length()
  let means = Array::make(columns, 0.0)
  for row in rows {
    for column in 0..<columns {
      means[column] += row[column]
    }
  }
  for column in 0..<columns {
    means[column] /= rows.length().to_double()
  }
  let result = matrix_zeros(columns, columns)
  for row in rows {
    for i in 0..<columns {
      for j in 0..<columns {
        result[i][j] += (row[i] - means[i]) * (row[j] - means[j])
      }
    }
  }
  let denominator = if sample && rows.length() > 1 {
    rows.length() - 1
  } else {
    rows.length()
  }
  if denominator > 0 {
    for i in 0..<columns {
      for j in 0..<columns {
        result[i][j] /= denominator.to_double()
      }
    }
  }
  result
}

///|
pub fn weighted_covariance_matrix(
  rows : Array[Array[Double]],
  weights : Array[Double],
) -> Array[Array[Double]] {
  if rows.length() == 0 || rows.length() != weights.length() {
    return []
  }
  let columns = rows[0].length()
  let means = Array::make(columns, 0.0)
  let total = sum(weights)
  if total == 0.0 {
    return matrix_zeros(columns, columns)
  }
  for i in 0..<rows.length() {
    for j in 0..<columns {
      means[j] += rows[i][j] * weights[i]
    }
  }
  for j in 0..<columns {
    means[j] /= total
  }
  let result = matrix_zeros(columns, columns)
  for i in 0..<rows.length() {
    for j in 0..<columns {
      for k in 0..<columns {
        result[j][k] += weights[i] *
          (rows[i][j] - means[j]) *
          (rows[i][k] - means[k])
      }
    }
  }
  for j in 0..<columns {
    for k in 0..<columns {
      result[j][k] /= total
    }
  }
  result
}
