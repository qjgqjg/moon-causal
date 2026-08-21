///|
/// Declarative validation rule for incoming analysis data.
pub struct ContractRule {
  name : String
  required : Bool
  minimum : Double
  maximum : Double
  allow_missing : Bool
}

///|
/// One data-contract violation.
pub struct ContractViolation {
  rule : String
  row : Int
  column : Int
  value : Double
  message : String
}

///|
/// Data-contract audit result.
pub struct ContractAudit {
  rows : Int
  columns : Int
  violations : Array[ContractViolation]
  checked_cells : Int
  score : Double
  passes : Bool
}

///|
/// Creates a numeric range rule.
pub fn contract_rule(
  name : String,
  minimum? : Double = -1.0e300,
  maximum? : Double = 1.0e300,
  required? : Bool = false,
  allow_missing? : Bool = true,
) -> ContractRule {
  {
    name,
    required,
    minimum: minimum.min(maximum),
    maximum: maximum.max(minimum),
    allow_missing,
  }
}

///|
/// Validates matrix cells against per-column rules.
pub fn audit_data_contract(
  matrix : Array[Array[Double]],
  rules : Array[ContractRule],
) -> ContractAudit {
  let violations : Array[ContractViolation] = Array::new()
  let mut checked = 0
  for row_index in 0..<matrix.length() {
    for column_index in 0..<rules.length() {
      let rule = rules[column_index]
      let present = column_index < matrix[row_index].length()
      let value = if present { matrix[row_index][column_index] } else { 0.0 }
      checked += 1
      if !present || (!is_finite(value) && rule.required && !rule.allow_missing) {
        violations.push({
          rule: rule.name,
          row: row_index,
          column: column_index,
          value,
          message: "required value is missing",
        })
      } else if is_finite(value) &&
        (value < rule.minimum || value > rule.maximum) {
        violations.push({
          rule: rule.name,
          row: row_index,
          column: column_index,
          value,
          message: "value outside contract range",
        })
      }
    }
  }
  let score = if checked == 0 {
    0.0
  } else {
    clamp(1.0 - violations.length().to_double() / checked.to_double(), 0.0, 1.0)
  }
  {
    rows: matrix.length(),
    columns: rules.length(),
    violations,
    checked_cells: checked,
    score,
    passes: violations.length() == 0,
  }
}

///|
/// Validates monotonic integer time within each unit.
pub fn audit_monotonic_panel(
  unit : Array[Int],
  time : Array[Int],
) -> Array[ContractViolation] {
  let result : Array[ContractViolation] = Array::new()
  let n = unit.length().min(time.length())
  for i in 1..<n {
    if unit[i] == unit[i - 1] && time[i] < time[i - 1] {
      result.push({
        rule: "monotonic-time",
        row: i,
        column: 1,
        value: time[i].to_double(),
        message: "time decreases within unit",
      })
    }
  }
  result
}

///|
/// Validates unique unit-time keys.
pub fn audit_unique_panel_keys(
  unit : Array[Int],
  time : Array[Int],
) -> Array[ContractViolation] {
  let result : Array[ContractViolation] = Array::new()
  let n = unit.length().min(time.length())
  for i in 0..<n {
    for j in 0..<i {
      if unit[i] == unit[j] && time[i] == time[j] {
        result.push({
          rule: "unique-unit-time",
          row: i,
          column: 0,
          value: time[i].to_double(),
          message: "duplicate unit-time key",
        })
        break
      }
    }
  }
  result
}

///|
/// Validates binary treatment support.
pub fn audit_treatment_support(treatment : Array[Bool]) -> ContractAudit {
  let treated = treatment.filter(fn(value) { value }).length()
  let control = treatment.length() - treated
  let violations : Array[ContractViolation] = Array::new()
  if treated == 0 {
    violations.push({
      rule: "treatment-support",
      row: -1,
      column: -1,
      value: 0.0,
      message: "no treated observations",
    })
  }
  if control == 0 {
    violations.push({
      rule: "treatment-support",
      row: -1,
      column: -1,
      value: 1.0,
      message: "no control observations",
    })
  }
  {
    rows: treatment.length(),
    columns: 1,
    violations,
    checked_cells: treatment.length(),
    score: if violations.length() == 0 {
      1.0
    } else {
      0.0
    },
    passes: violations.length() == 0,
  }
}

///|
/// Merges multiple contract violations into an audit.
pub fn merge_contract_audits(
  first : ContractAudit,
  second : ContractAudit,
) -> ContractAudit {
  let violations = first.violations.copy()
  for violation in second.violations {
    violations.push(violation)
  }
  let checked = first.checked_cells + second.checked_cells
  {
    rows: first.rows.max(second.rows),
    columns: first.columns.max(second.columns),
    violations,
    checked_cells: checked,
    score: if checked == 0 {
      0.0
    } else {
      clamp(
        1.0 - violations.length().to_double() / checked.to_double(),
        0.0,
        1.0,
      )
    },
    passes: violations.length() == 0,
  }
}

///|
/// Returns a compact data-contract summary vector.
pub fn contract_summary(audit : ContractAudit) -> Array[Double] {
  [
    audit.rows.to_double(),
    audit.columns.to_double(),
    audit.violations.length().to_double(),
    audit.checked_cells.to_double(),
    audit.score,
    if audit.passes {
      1.0
    } else {
      0.0
    },
  ]
}
