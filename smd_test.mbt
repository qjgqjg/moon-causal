///|
test "smd" {
  let cov = [1.0, 2.0, 1.0, 2.0]
  let treat = [true, true, false, false]
  let s = smd(cov, treat)
  inspect(s, content="0")
}
