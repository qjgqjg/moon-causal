///|
test "sigmoid function" {
  inspect(sigmoid(0.0), content="0.5")
  let s_inf = sigmoid(100.0)
  // Should be close to 1.0
  if s_inf < 0.99 {
    fail("Sigmoid of large value should be close to 1.0")
  }
}

///|
test "estimate propensity scores" {
  let covariates = [[1.0, 2.0], [1.0, -1.0], [-1.0, 0.0]]
  let weights = [0.5, -0.5]
  let scores = estimate_propensity_scores(covariates, weights)
  // score 1: z = 0.5 - 1.0 = -0.5 => sigmoid(-0.5)
  // score 2: z = 0.5 + 0.5 = 1.0 => sigmoid(1.0)
  // score 3: z = -0.5 + 0 = -0.5 => sigmoid(-0.5)
  if scores.length() != 3 {
    fail("Expected 3 scores")
  }
}
