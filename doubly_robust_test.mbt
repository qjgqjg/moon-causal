///|
test "doubly robust ate" {
  let y = [10.0, 15.0, 20.0, 25.0]
  let t = [true, false, true, false]
  let p = [0.5, 0.5, 0.5, 0.5]
  let y1 = [10.0, 10.0, 20.0, 20.0]
  let y0 = [5.0, 15.0, 15.0, 25.0]
  let ate = doubly_robust_ate(y, t, p, y1, y0)
  inspect(
    ate,
    content=(
      #|0
    ),
  )
}
