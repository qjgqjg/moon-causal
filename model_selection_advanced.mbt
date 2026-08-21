///|
/// Cross-validation score for a candidate model.
pub struct CandidateScore {
  candidate : Int
  train_loss : Double
  validation_loss : Double
  validation_metric : Double
  sample_size : Int
  converged : Bool
  selected : Bool
}

///|
/// Hyperparameter grid for regularized binary models.
pub struct RegularizationGrid {
  learning_rates : Array[Double]
  l1_values : Array[Double]
  l2_values : Array[Double]
  iterations : Int
}

///|
/// Selection result with all candidates and winner index.
pub struct SelectionResult {
  scores : Array[CandidateScore]
  best_candidate : Int
  best_metric : Double
  tie_count : Int
  passes : Bool
}

///|
/// Threshold-selection result for policy deployment.
pub struct ThresholdSelection {
  threshold : Double
  metric : Double
  predicted_positive_rate : Double
  true_positive_rate : Double
  false_positive_rate : Double
  objective : String
}

///|
/// Calibration summary for a fitted score.
pub struct CalibrationSummary {
  bins : Array[Array[Double]]
  expected_calibration_error : Double
  maximum_calibration_error : Double
  brier_score : Double
  passes : Bool
}

///|
/// Creates a regularization grid with safe defaults.
pub fn regularization_grid(
  learning_rates? : Array[Double] = [0.05, 0.1],
  l1_values? : Array[Double] = [0.0],
  l2_values? : Array[Double] = [1.0e-6, 1.0e-3],
  iterations? : Int = 400,
) -> RegularizationGrid {
  {
    learning_rates,
    l1_values,
    l2_values,
    iterations: if iterations > 0 {
      iterations
    } else {
      400
    },
  }
}

///|
/// Computes binary log loss for a score vector.
pub fn candidate_log_loss(
  probabilities : Array[Double],
  treatment : Array[Bool],
) -> Double {
  let n = probabilities.length().min(treatment.length())
  if n == 0 {
    return 0.0
  }
  let mut loss = 0.0
  for i in 0..<n {
    let probability = safe_probability(probabilities[i])
    let target = if treatment[i] { 1.0 } else { 0.0 }
    loss -= target * @math.ln(probability) +
      (1.0 - target) * @math.ln(1.0 - probability)
  }
  loss / n.to_double()
}

///|
/// Computes a balanced accuracy from boolean predictions.
pub fn balanced_accuracy(
  actual : Array[Bool],
  predicted : Array[Bool],
) -> Double {
  let counts = confusion_counts(actual, predicted)
  let sensitivity = if counts[2] + counts[3] == 0 {
    0.0
  } else {
    counts[3].to_double() / (counts[2] + counts[3]).to_double()
  }
  let specificity = if counts[0] + counts[1] == 0 {
    0.0
  } else {
    counts[0].to_double() / (counts[0] + counts[1]).to_double()
  }
  (sensitivity + specificity) / 2.0
}

///|
/// Finds a classification threshold maximizing Youden's J statistic.
pub fn optimize_youden_threshold(
  probabilities : Array[Double],
  treatment : Array[Bool],
  candidates? : Array[Double] = [],
) -> ThresholdSelection {
  let grid = if candidates.length() > 0 {
    candidates
  } else {
    let result : Array[Double] = Array::new()
    for i in 0..<=100 {
      result.push(i.to_double() / 100.0)
    }
    result
  }
  let mut best_threshold = 0.5
  let mut best_metric = -1.0
  let mut best_tpr = 0.0
  let mut best_fpr = 0.0
  for threshold in grid {
    let predicted = threshold_probabilities(probabilities, threshold)
    let counts = confusion_counts(treatment, predicted)
    let tpr = if counts[2] + counts[3] == 0 {
      0.0
    } else {
      counts[3].to_double() / (counts[2] + counts[3]).to_double()
    }
    let fpr = if counts[0] + counts[1] == 0 {
      0.0
    } else {
      counts[1].to_double() / (counts[0] + counts[1]).to_double()
    }
    let metric = tpr - fpr
    if metric > best_metric {
      best_metric = metric
      best_threshold = threshold
      best_tpr = tpr
      best_fpr = fpr
    }
  }
  let predicted = threshold_probabilities(probabilities, best_threshold)
  let mut positive = 0
  for value in predicted {
    if value {
      positive += 1
    }
  }
  {
    threshold: best_threshold,
    metric: best_metric,
    predicted_positive_rate: if predicted.length() == 0 {
      0.0
    } else {
      positive.to_double() / predicted.length().to_double()
    },
    true_positive_rate: best_tpr,
    false_positive_rate: best_fpr,
    objective: "youden-j",
  }
}

///|
/// Finds a threshold maximizing F1.
pub fn optimize_f1_threshold(
  probabilities : Array[Double],
  treatment : Array[Bool],
  candidates : Array[Double],
) -> ThresholdSelection {
  let mut best_threshold = 0.5
  let mut best_metric = -1.0
  for threshold in candidates {
    let predicted = threshold_probabilities(probabilities, threshold)
    let summary = classification_summary(treatment, predicted)
    if summary[3] > best_metric {
      best_metric = summary[3]
      best_threshold = threshold
    }
  }
  let predicted = threshold_probabilities(probabilities, best_threshold)
  let counts = confusion_counts(treatment, predicted)
  let positive = counts[1] + counts[3]
  let tpr = if counts[2] + counts[3] == 0 {
    0.0
  } else {
    counts[3].to_double() / (counts[2] + counts[3]).to_double()
  }
  let fpr = if counts[0] + counts[1] == 0 {
    0.0
  } else {
    counts[1].to_double() / (counts[0] + counts[1]).to_double()
  }
  {
    threshold: best_threshold,
    metric: best_metric,
    predicted_positive_rate: if predicted.length() == 0 {
      0.0
    } else {
      positive.to_double() / predicted.length().to_double()
    },
    true_positive_rate: tpr,
    false_positive_rate: fpr,
    objective: "f1",
  }
}

///|
/// Builds equal-width calibration bins and operational error diagnostics.
pub fn calibration_summary(
  probabilities : Array[Double],
  treatment : Array[Bool],
  bins? : Int = 10,
  maximum_error? : Double = 0.1,
) -> CalibrationSummary {
  let table = calibration_bins(probabilities, treatment, bins~)
  let mut expected = 0.0
  let mut maximum = 0.0
  let total = probabilities.length().to_double().max(1.0)
  for row in table {
    let error = (row[0] - row[1]).abs()
    expected += row[2] / total * error
    if error > maximum {
      maximum = error
    }
  }
  {
    bins: table,
    expected_calibration_error: expected,
    maximum_calibration_error: maximum,
    brier_score: brier_score(probabilities, treatment),
    passes: expected <= maximum_error && maximum <= maximum_error,
  }
}

///|
/// Creates sequential folds for time-ordered validation.
pub fn time_ordered_folds(sample_size : Int, folds : Int) -> Array[Fold] {
  let n = if sample_size > 0 { sample_size } else { 0 }
  let count = if folds > 1 { folds } else { 2 }
  let result : Array[Fold] = Array::new(capacity=count)
  for fold in 0..<count {
    let training_indices : Array[Int] = Array::new()
    let validation_indices : Array[Int] = Array::new()
    let validation_start = n * fold / count
    let validation_end = n * (fold + 1) / count
    for i in 0..<validation_start {
      training_indices.push(i)
    }
    for i in validation_start..<validation_end {
      validation_indices.push(i)
    }
    result.push({ training_indices, validation_indices })
  }
  result
}

///|
/// Scores a regularized logistic candidate on a supplied train/validation split.
pub fn score_logistic_candidate(
  covariates : Array[Array[Double]],
  treatment : Array[Bool],
  train_indices : Array[Int],
  validation_indices : Array[Int],
  learning_rate : Double,
  l1 : Double,
  l2 : Double,
  candidate : Int,
) -> CandidateScore {
  let train_x : Array[Array[Double]] = Array::new(
    capacity=train_indices.length(),
  )
  let train_y : Array[Bool] = Array::new(capacity=train_indices.length())
  let validation_x : Array[Array[Double]] = Array::new(
    capacity=validation_indices.length(),
  )
  let validation_y : Array[Bool] = Array::new(
    capacity=validation_indices.length(),
  )
  for index in train_indices {
    if index >= 0 && index < covariates.length() && index < treatment.length() {
      train_x.push(covariates[index])
      train_y.push(treatment[index])
    }
  }
  for index in validation_indices {
    if index >= 0 && index < covariates.length() && index < treatment.length() {
      validation_x.push(covariates[index])
      validation_y.push(treatment[index])
    }
  }
  let model = fit_logistic_regression(
    train_x,
    train_y,
    learning_rate~,
    max_iterations=400,
    l2=l2 + l1,
  )
  let probabilities = predict_propensity(model, validation_x)
  let loss = candidate_log_loss(probabilities, validation_y)
  let predictions = threshold_probabilities(probabilities, 0.5)
  {
    candidate,
    train_loss: model.loss,
    validation_loss: loss,
    validation_metric: balanced_accuracy(validation_y, predictions),
    sample_size: validation_y.length(),
    converged: model.converged,
    selected: false,
  }
}

///|
/// Selects the best candidate by validation loss, breaking ties by balanced accuracy.
pub fn select_best_candidate(scores : Array[CandidateScore]) -> SelectionResult {
  if scores.length() == 0 {
    return {
      scores: [],
      best_candidate: -1,
      best_metric: 0.0,
      tie_count: 0,
      passes: false,
    }
  }
  let mut best = 0
  let mut ties = 1
  for i in 1..<scores.length() {
    if scores[i].validation_loss < scores[best].validation_loss - 1.0e-12 {
      best = i
      ties = 1
    } else if (scores[i].validation_loss - scores[best].validation_loss).abs() <=
      1.0e-12 {
      ties += 1
      if scores[i].validation_metric > scores[best].validation_metric {
        best = i
      }
    }
  }
  let selected = Array::new(capacity=scores.length())
  for i in 0..<scores.length() {
    selected.push({
      candidate: scores[i].candidate,
      train_loss: scores[i].train_loss,
      validation_loss: scores[i].validation_loss,
      validation_metric: scores[i].validation_metric,
      sample_size: scores[i].sample_size,
      converged: scores[i].converged,
      selected: i == best,
    })
  }
  {
    scores: selected,
    best_candidate: scores[best].candidate,
    best_metric: scores[best].validation_metric,
    tie_count: ties,
    passes: scores[best].sample_size > 0 &&
    is_finite(scores[best].validation_loss),
  }
}

///|
/// Evaluates a grid over one train-validation split.
pub fn grid_search_logistic(
  covariates : Array[Array[Double]],
  treatment : Array[Bool],
  train_indices : Array[Int],
  validation_indices : Array[Int],
  grid : RegularizationGrid,
) -> SelectionResult {
  let scores : Array[CandidateScore] = Array::new()
  let mut candidate = 0
  for learning_rate in grid.learning_rates {
    for l1 in grid.l1_values {
      for l2 in grid.l2_values {
        scores.push(
          score_logistic_candidate(
            covariates, treatment, train_indices, validation_indices, learning_rate,
            l1, l2, candidate,
          ),
        )
        candidate += 1
      }
    }
  }
  select_best_candidate(scores)
}

///|
/// Computes a cost-sensitive threshold using false-positive and false-negative costs.
pub fn cost_sensitive_threshold(
  probabilities : Array[Double],
  treatment : Array[Bool],
  false_positive_cost : Double,
  false_negative_cost : Double,
) -> ThresholdSelection {
  let threshold = clamp(
    false_positive_cost /
    (false_positive_cost + false_negative_cost).max(1.0e-12),
    0.0,
    1.0,
  )
  let predicted = threshold_probabilities(probabilities, threshold)
  let counts = confusion_counts(treatment, predicted)
  let metric = -(false_positive_cost * counts[1].to_double() +
    false_negative_cost * counts[2].to_double())
  {
    threshold,
    metric,
    predicted_positive_rate: if predicted.length() == 0 {
      0.0
    } else {
      (counts[1] + counts[3]).to_double() / predicted.length().to_double()
    },
    true_positive_rate: if counts[2] + counts[3] == 0 {
      0.0
    } else {
      counts[3].to_double() / (counts[2] + counts[3]).to_double()
    },
    false_positive_rate: if counts[0] + counts[1] == 0 {
      0.0
    } else {
      counts[1].to_double() / (counts[0] + counts[1]).to_double()
    },
    objective: "cost-sensitive",
  }
}

///|
/// Returns a compact model-selection summary vector.
pub fn selection_summary(selection : SelectionResult) -> Array[Double] {
  [
    selection.best_candidate.to_double(),
    selection.best_metric,
    selection.tie_count.to_double(),
    selection.scores.length().to_double(),
    if selection.passes {
      1.0
    } else {
      0.0
    },
  ]
}
