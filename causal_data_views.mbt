///|
/// A named view over one numeric dataset column.
pub struct CausalColumnView {
  name : String
  index : Int
  values : Array[Double]
  finite_count : Int
  missing_count : Int
  mean : Double
  standard_deviation : Double
  minimum : Double
  maximum : Double
}

///|
/// A row-level view used by scoring and review workflows.
pub struct CausalRowView {
  index : Int
  covariates : Array[Double]
  treatment : Bool
  outcome : Double
  finite : Bool
  score : Double
}

///|
/// Deterministic train/validation split indices.
pub struct CausalDataSplit {
  train : Array[Int]
  validation : Array[Int]
  train_fraction : Double
  seed : UInt64
  passes : Bool
}

///|
/// Cross-fitting fold assignment and size diagnostics.
pub struct CausalFoldAssignment {
  folds : Array[Int]
  sizes : Array[Int]
  requested_folds : Int
  passes : Bool
}

///|
/// Treatment-stratified group summary.
pub struct CausalStratumSummary {
  stratum : Int
  count : Int
  treated : Int
  control : Int
  treatment_rate : Double
  outcome_mean : Double
}

///|
/// Dataset view diagnostics used by audit and model preparation.
pub struct CausalDataViewSummary {
  rows : Int
  columns : Int
  finite_rows : Int
  treated : Int
  control : Int
  missing_cells : Int
  constant_columns : Int
  outcome_mean : Double
  outcome_standard_deviation : Double
  passes : Bool
}

///|
/// Returns a sorted copy without mutating caller-owned values.
pub fn causal_view_sorted(values : Array[Double]) -> Array[Double] {
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
/// Creates a column view from a matrix with a stable fallback for short rows.
pub fn causal_column_view(
  matrix : Array[Array[Double]],
  index : Int,
  name? : String = "",
  fill? : Double = 0.0,
) -> CausalColumnView {
  let values : Array[Double] = Array::new(capacity=matrix.length())
  for row in matrix {
    values.push(
      if index >= 0 && index < row.length() {
        row[index]
      } else {
        fill
      },
    )
  }
  let finite = causal_finite_values(values)
  {
    name: if name == "" {
      "x\{index}"
    } else {
      name
    },
    index,
    values,
    finite_count: finite.length(),
    missing_count: values.length() - finite.length(),
    mean: mean_or(finite, 0.0),
    standard_deviation: if finite.length() > 1 {
      std_dev(finite)
    } else {
      0.0
    },
    minimum: causal_minimum(finite),
    maximum: causal_maximum(finite),
  }
}

///|
/// Creates all column views in matrix order.
pub fn causal_column_views(
  matrix : Array[Array[Double]],
  names? : Array[String] = [],
) -> Array[CausalColumnView] {
  let width = if matrix.length() == 0 {
    names.length()
  } else {
    let mut maximum = 0
    for row in matrix {
      if row.length() > maximum {
        maximum = row.length()
      }
    }
    maximum.max(names.length())
  }
  let result : Array[CausalColumnView] = Array::new(capacity=width)
  for index in 0..<width {
    let name = if index < names.length() { names[index] } else { "x\{index}" }
    result.push(causal_column_view(matrix, index, name~))
  }
  result
}

///|
/// Returns one row view with a finite-value and score diagnostic.
pub fn causal_row_view(
  dataset : CausalDataset,
  index : Int,
  score? : Double = 0.0,
) -> CausalRowView {
  let valid = index >= 0 && index < dataset.n()
  let covariates = if valid { dataset.covariates[index].copy() } else { [] }
  let treatment = if valid { dataset.treatment[index] } else { false }
  let outcome = if valid { dataset.outcome[index] } else { 0.0 }
  let finite = valid && is_finite(outcome) && finite_array(covariates)
  {
    index,
    covariates,
    treatment,
    outcome,
    finite,
    score: if is_finite(score) {
      score
    } else {
      0.0
    },
  }
}

///|
/// Returns all row views in dataset order.
pub fn causal_row_views(dataset : CausalDataset) -> Array[CausalRowView] {
  let result : Array[CausalRowView] = Array::new(capacity=dataset.n())
  for index in 0..<dataset.n() {
    result.push(causal_row_view(dataset, index))
  }
  result
}

///|
/// Computes the indices that sort a score vector ascending.
pub fn causal_order_indices(values : Array[Double]) -> Array[Int] {
  let order : Array[Int] = Array::new(capacity=values.length())
  for index in 0..<values.length() {
    order.push(index)
  }
  for i in 1..<order.length() {
    let index = order[i]
    let value = if is_finite(values[index]) { values[index] } else { 1.0e300 }
    let mut cursor = i
    while cursor > 0 {
      let previous = order[cursor - 1]
      let previous_value = if is_finite(values[previous]) {
        values[previous]
      } else {
        1.0e300
      }
      if previous_value <= value {
        break
      }
      order[cursor] = previous
      cursor -= 1
    }
    order[cursor] = index
  }
  order
}

///|
/// Splits indices into train and validation partitions reproducibly.
pub fn causal_train_validation_split(
  rows : Int,
  train_fraction? : Double = 0.8,
  seed? : UInt64 = 20260819,
) -> CausalDataSplit {
  let n = rows.max(0)
  let fraction = clamp(train_fraction, 0.0, 1.0)
  let target = (fraction * n.to_double()).to_int().min(n)
  let rng = RandomState::new(seed)
  let order : Array[Int] = Array::new(capacity=n)
  for index in 0..<n {
    order.push(index)
  }
  for i in 0..<n {
    let remaining = n - i
    if remaining > 1 {
      let offset = (rng.uniform() * remaining.to_double())
        .to_int()
        .min(remaining - 1)
      let swap = i + offset
      let temporary = order[i]
      order[i] = order[swap]
      order[swap] = temporary
    }
  }
  let train : Array[Int] = Array::new(capacity=target)
  let validation : Array[Int] = Array::new(capacity=n - target)
  for i in 0..<n {
    if i < target {
      train.push(order[i])
    } else {
      validation.push(order[i])
    }
  }
  {
    train,
    validation,
    train_fraction: if n == 0 {
      0.0
    } else {
      train.length().to_double() / n.to_double()
    },
    seed,
    passes: train.length() + validation.length() == n,
  }
}

///|
/// Splits rows while preserving the treatment rate as closely as possible.
pub fn causal_stratified_split(
  treatment : Array[Bool],
  train_fraction? : Double = 0.8,
  seed? : UInt64 = 20260819,
) -> CausalDataSplit {
  let fraction = clamp(train_fraction, 0.0, 1.0)
  let treated : Array[Int] = Array::new()
  let control : Array[Int] = Array::new()
  for index in 0..<treatment.length() {
    if treatment[index] {
      treated.push(index)
    } else {
      control.push(index)
    }
  }
  let treated_split = causal_train_validation_split(
    treated.length(),
    train_fraction=fraction,
    seed~,
  )
  let control_split = causal_train_validation_split(
    control.length(),
    train_fraction=fraction,
    seed=seed + 1UL,
  )
  let train : Array[Int] = Array::new()
  let validation : Array[Int] = Array::new()
  for index in treated_split.train {
    train.push(treated[index])
  }
  for index in control_split.train {
    train.push(control[index])
  }
  for index in treated_split.validation {
    validation.push(treated[index])
  }
  for index in control_split.validation {
    validation.push(control[index])
  }
  {
    train,
    validation,
    train_fraction: if treatment.length() == 0 {
      0.0
    } else {
      train.length().to_double() / treatment.length().to_double()
    },
    seed,
    passes: train.length() + validation.length() == treatment.length(),
  }
}

///|
/// Creates deterministic fold labels by shuffled round-robin assignment.
pub fn causal_fold_assignment(
  rows : Int,
  requested_folds : Int,
  seed? : UInt64 = 20260819,
) -> CausalFoldAssignment {
  let n = rows.max(0)
  let folds_count = requested_folds.max(1).min(n.max(1))
  let split = causal_train_validation_split(n, train_fraction=1.0, seed~)
  let labels = Array::make(n, 0)
  let sizes = Array::make(folds_count, 0)
  for i in 0..<split.train.length() {
    let fold = i % folds_count
    labels[split.train[i]] = fold
    sizes[fold] += 1
  }
  {
    folds: labels,
    sizes,
    requested_folds: folds_count,
    passes: sizes.fold(init=0, fn(total, size) { total + size }) == n,
  }
}

///|
/// Returns train and validation indices for one fold.
pub fn causal_fold_indices(
  assignment : CausalFoldAssignment,
  held_out_fold : Int,
) -> CausalDataSplit {
  let fold = held_out_fold.max(0).min(assignment.requested_folds - 1)
  let train : Array[Int] = Array::new()
  let validation : Array[Int] = Array::new()
  for index in 0..<assignment.folds.length() {
    if assignment.folds[index] == fold {
      validation.push(index)
    } else {
      train.push(index)
    }
  }
  {
    train,
    validation,
    train_fraction: if assignment.folds.length() == 0 {
      0.0
    } else {
      train.length().to_double() / assignment.folds.length().to_double()
    },
    seed: 0UL,
    passes: train.length() + validation.length() == assignment.folds.length(),
  }
}

///|
/// Produces k-fold index pairs for cross-fitting.
pub fn causal_all_folds(
  rows : Int,
  requested_folds : Int,
  seed? : UInt64 = 20260819,
) -> Array[CausalDataSplit] {
  let assignment = causal_fold_assignment(rows, requested_folds, seed~)
  let result : Array[CausalDataSplit] = Array::new(
    capacity=assignment.requested_folds,
  )
  for fold in 0..<assignment.requested_folds {
    result.push(causal_fold_indices(assignment, fold))
  }
  result
}

///|
/// Creates bootstrap indices with replacement.
pub fn causal_bootstrap_indices(
  rows : Int,
  replicates : Int,
  seed? : UInt64 = 20260819,
) -> Array[Array[Int]] {
  let n = rows.max(0)
  let count = replicates.max(0)
  let rng = RandomState::new(seed)
  let result : Array[Array[Int]] = Array::new(capacity=count)
  for _ in 0..<count {
    let sample : Array[Int] = Array::new(capacity=n)
    for _ in 0..<n {
      let index = if n == 0 {
        0
      } else {
        (rng.uniform() * n.to_double()).to_int().min(n - 1)
      }
      if n > 0 {
        sample.push(index)
      }
    }
    result.push(sample)
  }
  result
}

///|
/// Resamples a causal dataset using one bootstrap index vector.
pub fn causal_bootstrap_dataset(
  dataset : CausalDataset,
  indices : Array[Int],
) -> CausalDataset {
  dataset.select(indices)
}

///|
/// Returns the finite row mask for a rectangular dataset.
pub fn causal_finite_row_mask(dataset : CausalDataset) -> Array[Bool] {
  let result : Array[Bool] = Array::new(capacity=dataset.n())
  for index in 0..<dataset.n() {
    result.push(causal_row_view(dataset, index).finite)
  }
  result
}

///|
/// Returns indices for rows that pass a finite-value mask.
pub fn causal_mask_indices(mask : Array[Bool]) -> Array[Int] {
  let result : Array[Int] = Array::new()
  for index in 0..<mask.length() {
    if mask[index] {
      result.push(index)
    }
  }
  result
}

///|
/// Filters a dataset to complete finite rows.
pub fn causal_complete_rows(dataset : CausalDataset) -> CausalDataset {
  dataset.select(causal_mask_indices(causal_finite_row_mask(dataset)))
}

///|
/// Assigns rows to integer strata using a finite score and fixed width.
pub fn causal_score_strata(
  scores : Array[Double],
  strata : Int,
  minimum? : Double = 0.0,
  maximum? : Double = 1.0,
) -> Array[Int] {
  let count = strata.max(1)
  let lower = minimum.min(maximum)
  let upper = maximum.max(minimum)
  let width = (upper - lower).max(1.0e-12)
  scores.map(fn(score) {
    if !is_finite(score) {
      0
    } else {
      ((score - lower) / width * count.to_double())
      .to_int()
      .max(0)
      .min(count - 1)
    }
  })
}

///|
/// Summarizes treatment and outcome within integer strata.
pub fn causal_stratum_summaries(
  strata : Array[Int],
  treatment : Array[Bool],
  outcome : Array[Double],
  stratum_count : Int,
) -> Array[CausalStratumSummary] {
  let count = stratum_count.max(1)
  let sizes = Array::make(count, 0)
  let treated = Array::make(count, 0)
  let outcome_sum = Array::make(count, 0.0)
  let n = strata.length().min(treatment.length()).min(outcome.length())
  for i in 0..<n {
    let stratum = strata[i].max(0).min(count - 1)
    sizes[stratum] += 1
    if treatment[i] {
      treated[stratum] += 1
    }
    if is_finite(outcome[i]) {
      outcome_sum[stratum] += outcome[i]
    }
  }
  let result : Array[CausalStratumSummary] = Array::new(capacity=count)
  for stratum in 0..<count {
    let total = sizes[stratum]
    result.push({
      stratum,
      count: total,
      treated: treated[stratum],
      control: total - treated[stratum],
      treatment_rate: if total == 0 {
        0.0
      } else {
        treated[stratum].to_double() / total.to_double()
      },
      outcome_mean: if total == 0 {
        0.0
      } else {
        outcome_sum[stratum] / total.to_double()
      },
    })
  }
  result
}

///|
/// Reorders a dataset by a score vector.
pub fn causal_reorder_dataset(
  dataset : CausalDataset,
  scores : Array[Double],
) -> CausalDataset {
  let n = dataset.n().min(scores.length())
  let order = causal_order_indices(scores[:n].to_owned())
  dataset.select(order)
}

///|
/// Selects rows by a deterministic score quantile interval.
pub fn causal_quantile_slice(
  dataset : CausalDataset,
  scores : Array[Double],
  lower_probability : Double,
  upper_probability : Double,
) -> CausalDataset {
  let finite = causal_finite_values(scores)
  if finite.length() == 0 {
    return CausalDataset::new([], [], [])
  }
  let lower = quantile(finite, clamp(lower_probability, 0.0, 1.0))
  let upper = quantile(finite, clamp(upper_probability, 0.0, 1.0))
  causal_dataset_filter_score(dataset, scores, lower, upper)
}

///|
/// Builds a complete dataset view summary.
pub fn causal_data_view_summary(
  dataset : CausalDataset,
) -> CausalDataViewSummary {
  let mask = causal_finite_row_mask(dataset)
  let finite_rows = mask.fold(init=0, fn(total, value) {
    if value {
      total + 1
    } else {
      total
    }
  })
  let missing_cells = causal_matrix_missing(dataset.covariates)
  let constant_columns = causal_constant_columns(dataset.covariates)
  {
    rows: dataset.n(),
    columns: dataset.p(),
    finite_rows,
    treated: dataset.treated_count(),
    control: dataset.control_count(),
    missing_cells,
    constant_columns,
    outcome_mean: mean_or(causal_finite_values(dataset.outcome), 0.0),
    outcome_standard_deviation: std_dev(causal_finite_values(dataset.outcome)),
    passes: dataset.is_valid() &&
    finite_rows == dataset.n() &&
    dataset.treated_count() > 0 &&
    dataset.control_count() > 0,
  }
}

///|
/// Returns a compact data-view summary vector.
pub fn causal_data_view_summary_vector(
  summary : CausalDataViewSummary,
) -> Array[Double] {
  [
    summary.rows.to_double(),
    summary.columns.to_double(),
    summary.finite_rows.to_double(),
    summary.treated.to_double(),
    summary.control.to_double(),
    summary.missing_cells.to_double(),
    summary.constant_columns.to_double(),
    summary.outcome_mean,
    summary.outcome_standard_deviation,
    if summary.passes {
      1.0
    } else {
      0.0
    },
  ]
}

///|
/// Computes a stable fingerprint for split indices.
pub fn causal_split_fingerprint(split : CausalDataSplit) -> UInt64 {
  let rows : Array[Array[Double]] = Array::new()
  rows.push([split.train_fraction, split.seed.to_double()])
  for index in split.train {
    rows.push([1.0, index.to_double()])
  }
  for index in split.validation {
    rows.push([0.0, index.to_double()])
  }
  matrix_checksum(rows)
}

///|
/// Computes a stable fingerprint for fold assignments.
pub fn causal_fold_fingerprint(assignment : CausalFoldAssignment) -> UInt64 {
  let rows : Array[Array[Double]] = Array::new()
  for index in 0..<assignment.folds.length() {
    rows.push([index.to_double(), assignment.folds[index].to_double()])
  }
  matrix_checksum(rows)
}
