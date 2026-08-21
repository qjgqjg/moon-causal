///|
/// Columnar table schema for analysis pipelines.
pub struct TableSchema {
  names : Array[String]
  types : Array[String]
  required : Array[Bool]
  version : String
}

///|
/// Numeric table profile.
pub struct TableProfile {
  rows : Int
  columns : Int
  complete_rows : Int
  missing_cells : Int
  duplicate_rows : Int
  schema_fingerprint : UInt64
  data_fingerprint : UInt64
  passes : Bool
}

///|
/// Group-level aggregate with count, mean and dispersion.
pub struct GroupAggregate {
  key : Int
  count : Int
  total : Double
  mean : Double
  variance : Double
  minimum : Double
  maximum : Double
}

///|
/// Result of a keyed inner or left join.
pub struct TableJoinResult {
  left_indices : Array[Int]
  right_indices : Array[Int]
  matched : Int
  unmatched_left : Int
  unmatched_right : Int
  passes : Bool
}

///|
/// A deterministic rolling summary for one key and time index.
pub struct TableWindowSummary {
  key : Int
  time : Int
  count : Int
  mean : Double
  standard_deviation : Double
  minimum : Double
  maximum : Double
}

///|
fn table_unique_int(values : Array[Int]) -> Array[Int] {
  let result : Array[Int] = Array::new()
  for value in values {
    if !result.contains(value) {
      result.push(value)
    }
  }
  result
}

///|
fn table_min(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let mut result = values[0]
  for value in values[1:] {
    if value < result {
      result = value
    }
  }
  result
}

///|
fn table_max(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let mut result = values[0]
  for value in values[1:] {
    if value > result {
      result = value
    }
  }
  result
}

///|
/// Creates a schema with stable defaults.
pub fn table_schema(
  names : Array[String],
  types? : Array[String] = [],
  required? : Array[Bool] = [],
  version? : String = "1",
) -> TableSchema {
  let actual_types : Array[String] = Array::new(capacity=names.length())
  let actual_required : Array[Bool] = Array::new(capacity=names.length())
  for i in 0..<names.length() {
    actual_types.push(if i < types.length() { types[i] } else { "double" })
    actual_required.push(
      if i < required.length() {
        required[i]
      } else {
        false
      },
    )
  }
  {
    names: names.copy(),
    types: actual_types,
    required: actual_required,
    version,
  }
}

///|
/// Computes a schema fingerprint.
pub fn table_schema_fingerprint(schema : TableSchema) -> UInt64 {
  let rows : Array[Array[Double]] = Array::new()
  for i in 0..<schema.names.length() {
    rows.push([
      i.to_double(),
      schema.types[i].length().to_double(),
      if schema.required[i] {
        1.0
      } else {
        0.0
      },
    ])
  }
  matrix_checksum(rows)
}

///|
/// Selects columns from a numeric table.
pub fn table_select_columns(
  table : Array[Array[Double]],
  indexes : Array[Int],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=table.length())
  for row in table {
    let selected = Array::new(capacity=indexes.length())
    for index in indexes {
      if index >= 0 && index < row.length() {
        selected.push(row[index])
      }
    }
    result.push(selected)
  }
  result
}

///|
/// Filters rows using an explicit mask.
pub fn table_filter_rows(
  table : Array[Array[Double]],
  mask : Array[Bool],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  for i in 0..<table.length().min(mask.length()) {
    if mask[i] {
      result.push(table[i])
    }
  }
  result
}

///|
/// Filters rows using a numeric predicate expressed as a threshold.
pub fn table_filter_threshold(
  table : Array[Array[Double]],
  column : Int,
  threshold : Double,
  keep_above? : Bool = true,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  for row in table {
    if column >= 0 && column < row.length() {
      let passes = if keep_above {
        row[column] >= threshold
      } else {
        row[column] <= threshold
      }
      if passes {
        result.push(row)
      }
    }
  }
  result
}

///|
/// Adds a numeric column, truncating to the existing row count when needed.
pub fn table_add_column(
  table : Array[Array[Double]],
  column : Array[Double],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=table.length())
  for i in 0..<table.length() {
    let row = table[i].copy()
    row.push(if i < column.length() { column[i] } else { 0.0 })
    result.push(row)
  }
  result
}

///|
/// Replaces a numeric column with a transformed vector.
pub fn table_replace_column(
  table : Array[Array[Double]],
  column_index : Int,
  values : Array[Double],
) -> Array[Array[Double]] {
  let result = table.copy()
  for i in 0..<result.length() {
    if column_index >= 0 &&
      column_index < result[i].length() &&
      i < values.length() {
      result[i][column_index] = values[i]
    }
  }
  result
}

///|
/// Profiles a numeric table and its schema.
pub fn profile_table(
  table : Array[Array[Double]],
  schema : TableSchema,
) -> TableProfile {
  let missing = assess_matrix(table).missing_cells
  let complete = complete_row_count(table)
  let duplicates = duplicate_row_flags(table)
    .filter(fn(value) { value })
    .length()
  let columns = if table.length() == 0 { 0 } else { table[0].length() }
  let passes = table.length() == complete &&
    duplicates == 0 &&
    schema.names.length() == columns
  {
    rows: table.length(),
    columns,
    complete_rows: complete,
    missing_cells: missing,
    duplicate_rows: duplicates,
    schema_fingerprint: table_schema_fingerprint(schema),
    data_fingerprint: matrix_checksum(table),
    passes,
  }
}

///|
/// Counts rows without non-finite cells.
pub fn complete_row_count(table : Array[Array[Double]]) -> Int {
  let mut result = 0
  for row in table {
    let mut complete = true
    for value in row {
      if !is_finite(value) {
        complete = false
      }
    }
    if complete {
      result += 1
    }
  }
  result
}

///|
/// Aggregates a value by integer group key.
pub fn group_aggregate(
  keys : Array[Int],
  values : Array[Double],
) -> Array[GroupAggregate] {
  let n = keys.length().min(values.length())
  let groups = table_unique_int(keys[:n].to_owned())
  let result : Array[GroupAggregate] = Array::new(capacity=groups.length())
  for key in groups {
    let selected = Array::new()
    for i in 0..<n {
      if keys[i] == key {
        selected.push(values[i])
      }
    }
    result.push({
      key,
      count: selected.length(),
      total: sum(selected),
      mean: mean_or(selected, 0.0),
      variance: variance(selected),
      minimum: table_min(selected),
      maximum: table_max(selected),
    })
  }
  result
}

///|
/// Computes weighted group means.
pub fn weighted_group_means(
  keys : Array[Int],
  values : Array[Double],
  weights : Array[Double],
) -> Array[Array[Double]] {
  let n = keys.length().min(values.length()).min(weights.length())
  let groups = table_unique_int(keys[:n].to_owned())
  let result : Array[Array[Double]] = Array::new(capacity=groups.length())
  for key in groups {
    let mut numerator = 0.0
    let mut denominator = 0.0
    for i in 0..<n {
      if keys[i] == key {
        numerator += values[i] * weights[i]
        denominator += weights[i]
      }
    }
    result.push([
      key.to_double(),
      if denominator == 0.0 {
        0.0
      } else {
        numerator / denominator
      },
      denominator,
    ])
  }
  result
}

///|
/// Computes group-level differences between two treatment arms.
pub fn group_treatment_effect(
  keys : Array[Int],
  treatment : Array[Bool],
  outcomes : Array[Double],
) -> Array[Array[Double]] {
  let n = keys.length().min(treatment.length()).min(outcomes.length())
  let groups = table_unique_int(keys[:n].to_owned())
  let result : Array[Array[Double]] = Array::new(capacity=groups.length())
  for key in groups {
    let treated = Array::new()
    let control = Array::new()
    for i in 0..<n {
      if keys[i] == key {
        if treatment[i] {
          treated.push(outcomes[i])
        } else {
          control.push(outcomes[i])
        }
      }
    }
    result.push([
      key.to_double(),
      mean_or(treated, 0.0) - mean_or(control, 0.0),
      treated.length().to_double(),
      control.length().to_double(),
    ])
  }
  result
}

///|
/// Performs a one-to-one keyed join using integer identifiers.
pub fn table_inner_join(
  left_keys : Array[Int],
  right_keys : Array[Int],
) -> TableJoinResult {
  let left_indices : Array[Int] = Array::new()
  let right_indices : Array[Int] = Array::new()
  let matched_right : Array[Int] = Array::new()
  for i in 0..<left_keys.length() {
    let mut match_index = -1
    for j in 0..<right_keys.length() {
      if right_keys[j] == left_keys[i] && !matched_right.contains(j) {
        match_index = j
        break
      }
    }
    if match_index >= 0 {
      left_indices.push(i)
      right_indices.push(match_index)
      matched_right.push(match_index)
    }
  }
  {
    left_indices,
    right_indices,
    matched: left_indices.length(),
    unmatched_left: left_keys.length() - left_indices.length(),
    unmatched_right: right_keys.length() - matched_right.length(),
    passes: left_indices.length() > 0 ||
    (left_keys.length() == 0 && right_keys.length() == 0),
  }
}

///|
/// Performs a left keyed join and marks absent right rows with -1.
pub fn table_left_join(
  left_keys : Array[Int],
  right_keys : Array[Int],
) -> TableJoinResult {
  let left_indices : Array[Int] = Array::new()
  let right_indices : Array[Int] = Array::new()
  let used : Array[Int] = Array::new()
  for i in 0..<left_keys.length() {
    left_indices.push(i)
    let mut match_index = -1
    for j in 0..<right_keys.length() {
      if right_keys[j] == left_keys[i] && !used.contains(j) {
        match_index = j
        break
      }
    }
    right_indices.push(match_index)
    if match_index >= 0 {
      used.push(match_index)
    }
  }
  {
    left_indices,
    right_indices,
    matched: used.length(),
    unmatched_left: left_keys.length() - used.length(),
    unmatched_right: right_keys.length() - used.length(),
    passes: true,
  }
}

///|
/// Calculates a lag within each integer group after assuming row order is time order.
pub fn group_lag(
  keys : Array[Int],
  values : Array[Double],
  lag : Int,
  fill? : Double = 0.0,
) -> Array[Double] {
  let result = Array::make(keys.length().min(values.length()), fill)
  if lag <= 0 {
    return values[:result.length()].to_owned()
  }
  for i in lag..<result.length() {
    if keys[i] == keys[i - lag] {
      result[i] = values[i - lag]
    }
  }
  result
}

///|
/// Computes a grouped rolling mean.
pub fn group_rolling_mean(
  keys : Array[Int],
  values : Array[Double],
  window : Int,
) -> Array[Double] {
  let n = keys.length().min(values.length())
  let result = Array::make(n, 0.0)
  let width = if window > 0 { window } else { 1 }
  for i in 0..<n {
    let selected = Array::new()
    let start = if i + 1 > width { i + 1 - width } else { 0 }
    let end = i + 1
    for j in start..<end {
      if keys[j] == keys[i] {
        selected.push(values[j])
      }
    }
    result[i] = mean_or(selected, 0.0)
  }
  result
}

///|
/// Computes a grouped rolling standard deviation.
pub fn group_rolling_std(
  keys : Array[Int],
  values : Array[Double],
  window : Int,
) -> Array[Double] {
  let n = keys.length().min(values.length())
  let result = Array::make(n, 0.0)
  let width = if window > 0 { window } else { 1 }
  for i in 0..<n {
    let selected = Array::new()
    let start = if i + 1 > width { i + 1 - width } else { 0 }
    let end = i + 1
    for j in start..<end {
      if keys[j] == keys[i] {
        selected.push(values[j])
      }
    }
    result[i] = std_dev(selected)
  }
  result
}

///|
/// Returns one summary row for each key-time pair in a panel-like table.
pub fn table_window_summaries(
  keys : Array[Int],
  times : Array[Int],
  values : Array[Double],
  window : Int,
) -> Array[TableWindowSummary] {
  let n = keys.length().min(times.length()).min(values.length())
  let result : Array[TableWindowSummary] = Array::new(capacity=n)
  let width = if window > 0 { window } else { 1 }
  for i in 0..<n {
    let selected = Array::new()
    for j in 0..<n {
      if keys[j] == keys[i] &&
        times[j] <= times[i] &&
        times[j] >= times[i] - width + 1 {
        selected.push(values[j])
      }
    }
    result.push({
      key: keys[i],
      time: times[i],
      count: selected.length(),
      mean: mean_or(selected, 0.0),
      standard_deviation: std_dev(selected),
      minimum: table_min(selected),
      maximum: table_max(selected),
    })
  }
  result
}
