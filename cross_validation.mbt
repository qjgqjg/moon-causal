///|
pub struct Fold {
  training_indices : Array[Int]
  validation_indices : Array[Int]
}

///|
pub fn make_folds(sample_size : Int, folds : Int, seed : UInt64) -> Array[Fold] {
  let actual_folds = if folds < 2 {
    1
  } else if folds > sample_size {
    sample_size
  } else {
    folds
  }
  if sample_size <= 0 {
    return []
  }
  let order = Array::new(capacity=sample_size)
  for i in 0..<sample_size {
    order.push(i)
  }
  let rng = RandomState::new(seed)
  for i in 0..<sample_size {
    let j = (rng.uniform() * sample_size.to_double()).to_int()
    let temporary = order[i]
    order[i] = order[j]
    order[j] = temporary
  }
  let result : Array[Fold] = Array::new(capacity=actual_folds)
  for fold in 0..<actual_folds {
    let validation = Array::new()
    let training = Array::new()
    for position in 0..<sample_size {
      if position % actual_folds == fold {
        validation.push(order[position])
      } else {
        training.push(order[position])
      }
    }
    result.push({ training_indices: training, validation_indices: validation })
  }
  result
}

///|
pub fn accuracy_from_labels(
  predictions : Array[Bool],
  actual : Array[Bool],
) -> Double {
  let n = if predictions.length() < actual.length() {
    predictions.length()
  } else {
    actual.length()
  }
  if n == 0 {
    return 0.0
  }
  let mut correct = 0
  for i in 0..<n {
    if predictions[i] == actual[i] {
      correct += 1
    }
  }
  correct.to_double() / n.to_double()
}

///|
pub fn cross_validate_propensity(
  covariates : Array[Array[Double]],
  treatment : Array[Bool],
  folds : Int,
  seed : UInt64,
) -> Double {
  let partitions = make_folds(covariates.length(), folds, seed)
  let scores = Array::new()
  let labels = Array::new()
  for partition in partitions {
    let train_x : Array[Array[Double]] = Array::new()
    let train_t : Array[Bool] = Array::new()
    for index in partition.training_indices {
      train_x.push(covariates[index])
      train_t.push(treatment[index])
    }
    let model = fit_logistic_regression(train_x, train_t, max_iterations=200)
    for index in partition.validation_indices {
      scores.push(model_predict(model, covariates[index]))
      labels.push(treatment[index])
    }
  }
  let predicted = Array::new(capacity=scores.length())
  for score in scores {
    predicted.push(score >= 0.0)
  }
  accuracy_from_labels(predicted, labels)
}

///|
pub fn cross_validate_outcome_rmse(
  covariates : Array[Array[Double]],
  outcomes : Array[Double],
  folds : Int,
  seed : UInt64,
) -> Double {
  let partitions = make_folds(covariates.length(), folds, seed)
  let errors = Array::new()
  for partition in partitions {
    let train_x : Array[Array[Double]] = Array::new()
    let train_y : Array[Double] = Array::new()
    for index in partition.training_indices {
      train_x.push(covariates[index])
      train_y.push(outcomes[index])
    }
    let model = fit_linear_outcome_model(train_x, train_y)
    for index in partition.validation_indices {
      let prediction = model_predict(model, covariates[index])
      let difference = prediction - outcomes[index]
      errors.push(difference * difference)
    }
  }
  if errors.length() == 0 {
    0.0
  } else {
    mean(errors).sqrt()
  }
}

///|
pub fn auc_score(
  probabilities : Array[Double],
  treatment : Array[Bool],
) -> Double {
  let positive = Array::new()
  let negative = Array::new()
  let n = if probabilities.length() < treatment.length() {
    probabilities.length()
  } else {
    treatment.length()
  }
  for i in 0..<n {
    if treatment[i] {
      positive.push(probabilities[i])
    } else {
      negative.push(probabilities[i])
    }
  }
  if positive.length() == 0 || negative.length() == 0 {
    return 0.5
  }
  let mut wins = 0.0
  for p in positive {
    for n in negative {
      if p > n {
        wins += 1.0
      } else if p == n {
        wins += 0.5
      }
    }
  }
  wins / (positive.length() * negative.length()).to_double()
}
