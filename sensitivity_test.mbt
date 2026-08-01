test "rosenbaum bounds gamma" {
  inspect(rosenbaum_bounds_gamma(0.05, 2.0), content="0.1")
  inspect(rosenbaum_bounds_gamma(0.05, 50.0), content="1")
}
