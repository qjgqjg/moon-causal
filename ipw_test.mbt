test "ipw ate weights" {
  let t = [true, false, true, false]
  let p = [0.8, 0.2, 0.5, 0.5]
  let weights = ipw_weights_ate(t, p)
  inspect(weights[0], content="1.25")
  inspect(weights[1], content="1.25")
  inspect(weights[2], content="2")
  inspect(weights[3], content="2")
}

test "ipw att weights" {
  let t = [true, false, true, false]
  let p = [0.8, 0.2, 0.5, 0.5]
  let weights = ipw_weights_att(t, p)
  inspect(weights[0], content="1")
  inspect(weights[1], content="0.25")
  inspect(weights[2], content="1")
  inspect(weights[3], content="1")
}
