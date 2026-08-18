///|
fn euclidean_distance(
  a : Array[Double],
  b : Array[Double],
  scales : Array[Double],
) -> Double {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let mut distance = 0.0
  for j in 0..<n {
    let scale = if j < scales.length() && scales[j] > 1.0e-12 {
      scales[j]
    } else {
      1.0
    }
    let difference = (a[j] - b[j]) / scale
    distance += difference * difference
  }
  distance.sqrt()
}

///|
fn column_scales(covariates : Array[Array[Double]]) -> Array[Double] {
  if covariates.length() == 0 {
    return []
  }
  let width = covariates[0].length()
  let result = Array::make(width, 0.0)
  for row in covariates {
    for j in 0..<width {
      result[j] += row[j] * row[j]
    }
  }
  for j in 0..<width {
    result[j] = (result[j] / covariates.length().to_double()).sqrt()
    if result[j] < 1.0e-12 {
      result[j] = 1.0
    }
  }
  result
}

///|
fn contains_int(values : Array[Int], target : Int) -> Bool {
  for value in values {
    if value == target {
      return true
    }
  }
  false
}

///|
/// Performs one-to-one nearest-neighbor matching on standardized covariates.
pub fn nearest_neighbor_match(
  covariates : Array[Array[Double]],
  treatment : Array[Bool],
  caliper? : Double = 0.0,
  replace? : Bool = false,
) -> Array[MatchPair] {
  let scales = column_scales(covariates)
  let controls = Array::new()
  for i in 0..<treatment.length() {
    if !treatment[i] {
      controls.push(i)
    }
  }
  let used = Array::new()
  let pairs : Array[MatchPair] = Array::new()
  for treated_index in 0..<treatment.length() {
    if treatment[treated_index] {
      let mut best_index = -1
      let mut best_distance = 1.0e300
      for control_index in controls {
        if replace || !contains_int(used, control_index) {
          let distance = euclidean_distance(
            covariates[treated_index],
            covariates[control_index],
            scales,
          )
          if distance < best_distance {
            best_distance = distance
            best_index = control_index
          }
        }
      }
      if best_index >= 0 && (caliper <= 0.0 || best_distance <= caliper) {
        pairs.push({
          treated_index,
          control_index: best_index,
          distance: best_distance,
          weight: 1.0,
        })
        if !replace {
          used.push(best_index)
        }
      }
    }
  }
  pairs
}

///|
/// Computes the matched-sample effect and its pair-level standard error.
pub fn estimate_matched_ate(
  outcomes : Array[Double],
  pairs : Array[MatchPair],
) -> Estimate {
  let differences = Array::new(capacity=pairs.length())
  for pair in pairs {
    if pair.treated_index >= 0 &&
      pair.control_index >= 0 &&
      pair.treated_index < outcomes.length() &&
      pair.control_index < outcomes.length() {
      differences.push(
        outcomes[pair.treated_index] - outcomes[pair.control_index],
      )
    }
  }
  let n = differences.length()
  let se = if n <= 1 {
    0.0
  } else {
    std_dev(differences, sample=true) / n.to_double().sqrt()
  }
  Estimate::from_standard_error(
    mean(differences),
    se,
    n,
    n.to_double(),
    "ATE (nearest-neighbor matching)",
  )
}

///|
pub fn matching_effects(
  outcomes : Array[Double],
  pairs : Array[MatchPair],
) -> Array[Double] {
  let result = Array::new(capacity=pairs.length())
  for pair in pairs {
    if pair.treated_index < outcomes.length() &&
      pair.control_index < outcomes.length() {
      result.push(outcomes[pair.treated_index] - outcomes[pair.control_index])
    }
  }
  result
}

///|
pub fn matching_mean_distance(pairs : Array[MatchPair]) -> Double {
  let distances = Array::new(capacity=pairs.length())
  for pair in pairs {
    distances.push(pair.distance)
  }
  mean(distances)
}
