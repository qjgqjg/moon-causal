///|
pub struct DesignMatrix {
  values : Array[Array[Double]]
  column_names : Array[String]
}

///|
pub fn make_design_matrix(
  covariates : Array[Array[Double]],
  include_intercept : Bool,
) -> DesignMatrix {
  let width = if covariates.length() == 0 { 0 } else { covariates[0].length() }
  let intercept_columns = if include_intercept { 1 } else { 0 }
  let names : Array[String] = Array::new(capacity=width + intercept_columns)
  if include_intercept {
    names.push("intercept")
  }
  for j in 0..<width {
    names.push("x\{j + 1}")
  }
  let values : Array[Array[Double]] = Array::new(capacity=covariates.length())
  for row in covariates {
    let expanded : Array[Double] = Array::new(
      capacity=row.length() + intercept_columns,
    )
    if include_intercept {
      expanded.push(1.0)
    }
    for value in row {
      expanded.push(value)
    }
    values.push(expanded)
  }
  { values, column_names: names }
}

///|
pub fn add_treatment_interactions(
  covariates : Array[Array[Double]],
  treatment : Array[Bool],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=covariates.length())
  for i in 0..<covariates.length() {
    let row = covariates[i].copy()
    let indicator = if treatment[i] { 1.0 } else { 0.0 }
    for value in covariates[i] {
      row.push(value * indicator)
    }
    result.push(row)
  }
  result
}

///|
pub fn append_column(
  matrix : Array[Array[Double]],
  column : Array[Double],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  let n = if matrix.length() < column.length() {
    matrix.length()
  } else {
    column.length()
  }
  for i in 0..<n {
    let row = matrix[i].copy()
    row.push(column[i])
    result.push(row)
  }
  result
}

///|
pub fn rowwise_linear_predict(
  design : Array[Array[Double]],
  coefficients : Array[Double],
) -> Array[Double] {
  let result = Array::new(capacity=design.length())
  for row in design {
    result.push(dot(row, coefficients))
  }
  result
}

///|
pub fn marginal_effects(
  model : ModelFit,
  covariates : Array[Array[Double]],
) -> Array[Double] {
  let result = Array::new(capacity=model.coefficients.length())
  for coefficient in model.coefficients {
    result.push(coefficient)
  }
  let _ = covariates
  result
}

///|
pub fn interaction_average_effect(
  model : ModelFit,
  covariates : Array[Array[Double]],
  treatment_column : Int,
) -> Double {
  if covariates.length() == 0 {
    return 0.0
  }
  let treated = covariates.copy()
  let control = covariates.copy()
  for row in treated {
    if treatment_column >= 0 && treatment_column < row.length() {
      row[treatment_column] = 1.0
    }
  }
  for row in control {
    if treatment_column >= 0 && treatment_column < row.length() {
      row[treatment_column] = 0.0
    }
  }
  mean(
    subtract_vectors(
      predict_outcomes(model, treated),
      predict_outcomes(model, control),
    ),
  )
}

///|
pub fn polynomial_design(
  covariates : Array[Array[Double]],
  degree : Int,
) -> DesignMatrix {
  let values = polynomial_features(covariates, degree)
  let names : Array[String] = Array::new()
  if values.length() > 0 {
    for j in 0..<values[0].length() {
      names.push("feature-\{j + 1}")
    }
  }
  { values, column_names: names }
}
