test "mean and variance" {
  let arr = [1.0, 2.0, 3.0, 4.0, 5.0]
  inspect(mean(arr), content="3.0")
  inspect(variance(arr), content="2.5")
}
