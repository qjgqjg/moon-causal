test "stratify_by_score" {
  let scores = [0.1, 0.4, 0.6, 0.9]
  let strata = stratify_by_score(scores, 4)
  inspect(strata[0], content="0")
  inspect(strata[1], content="1")
  inspect(strata[2], content="2")
  inspect(strata[3], content="3")
}
