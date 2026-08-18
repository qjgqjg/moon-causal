///|
fn main {
  let config = @causal.SyntheticConfig::new(200, 4, 1.75, 20260818UL)
  let report = @causal.run_synthetic_benchmark(config)
  println("moon-causal deterministic benchmark")
  println("sample_size=\{report.sample_size()}")
  println("true_ate=\{report.true_ate()}")
  println("estimated_ate=\{report.estimated_ate()}")
  println("absolute_error=\{report.absolute_error()}")
  println("standard_error=\{report.standard_error()}")
  println("effective_sample_size=\{report.effective_sample_size()}")
  println("propensity_auc=\{report.propensity_auc()}")
}
