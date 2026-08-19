# moon-causal

moon-causal is a pure MoonBit toolkit for causal inference on observational and experimental data. It provides composable data contracts, propensity models, treatment-effect estimators, balance diagnostics, uncertainty quantification, and deterministic simulation utilities for applications such as A/B testing, policy evaluation, and medical research.

## Core capabilities

- Data preparation: rectangular dataset validation, summaries, scaling, transformations, missing-value handling, and train/validation splits.
- Modeling: standardized logistic propensity regression, linear outcome regression, design matrices, polynomial terms, and treatment interactions.
- Estimation: IPW, stabilized IPW, ATT, ATC, overlap weighting, g-computation, doubly robust estimation, stratification, and nearest-neighbor matching.
- Diagnostics: standardized mean differences, variance ratios, overlap and positivity checks, calibration, effective sample size, residual diagnostics, and quality reports.
- Extended methods: subgroup effects, policy value, multi-arm treatments, continuous-dose response, difference-in-differences, event studies, synthetic controls, mediation, and survival curves.
- Uncertainty and evaluation: bootstrap intervals, permutation inference, jackknife influence values, cross-validation, classification metrics, and deterministic synthetic benchmarks.

## Quick start

Requirements: MoonBit stable toolchain.

```bash
moon check --deny-warn
moon test --deny-warn
moon run cmd/main
```

The library has no runtime dependency outside the MoonBit toolchain. The command-line example runs a deterministic observational benchmark and prints the estimated effect, uncertainty, effective sample size, and propensity AUC.

## CLI

Run the included example from the repository root:

```bash
moon run cmd/main
```

Expected output for the default configuration:

```text
moon-causal deterministic benchmark
sample_size=200
true_ate=1.75
estimated_ate=1.695501947510638
absolute_error=0.0544980524893619
standard_error=0.15562533888456542
effective_sample_size=177.61530018366125
propensity_auc=0.6822682268226823
```

The benchmark uses 200 observations, 4 covariates, treatment effect 1.75, and seed 20260818UL. It is intended as a reproducible regression and integration fixture, not as a claim about performance on a particular real-world dataset.

## Architecture

The root package is organized by responsibility:

- core_types.mbt, validation.mbt, and preprocessing.mbt define data contracts and transformations.
- linear_algebra.mbt, models.mbt, formula.mbt, and metrics.mbt provide numerical foundations and predictive models.
- estimators.mbt, matching.mbt, causal_methods.mbt, and advanced_causal.mbt implement treatment-effect estimators.
- balance.mbt, diagnostics.mbt, audit.mbt, and sensitivity.mbt provide diagnostics and robustness analysis.
- heterogeneity.mbt, multi_arm.mbt, continuous.mbt, panel.mbt, survival.mbt, and mediation.mbt cover extended designs.
- uncertainty.mbt, bootstrap_advanced.mbt, cross_validation.mbt, and simulation.mbt provide inference and reproducible evaluation.
- cmd/main contains a runnable benchmark application.

## Benchmark and reproducibility

The simulation API is deterministic for a fixed SyntheticConfig seed. run_synthetic_benchmark reports:

- sample size and configured ground-truth ATE;
- estimated IPW ATE and absolute error;
- standard error and effective sample size;
- propensity-model AUC.

This makes the benchmark suitable for regression testing, examples, and comparing implementation changes without relying on external data.

## Tests and quality checks

The test suite covers normal paths, empty and singleton inputs, dimension mismatches, non-finite values, extreme probabilities, repeated survival times, numerical boundaries, integration paths, and deterministic simulation.

```bash
moon fmt --check
moon info
moon check --deny-warn
moon test --deny-warn
moon build --target wasm-gc
moon build --target native
```

## CI

GitHub Actions installs the latest MoonBit stable toolchain from the official installer and runs:

- format and warning-as-error checks;
- wasm-gc and native builds;
- public-interface regeneration and clean-diff verification;
- wasm-gc tests with coverage summary;
- native type checks and tests.

## License

Apache License 2.0. See [LICENSE](LICENSE).

## Project links

- [GitHub](https://github.com/qjgqjg/moon-causal)
- [Gitlink](https://gitlink.org.cn/qjgqjg68/moon-causal)

```mbt check
///|
test "README example" {
  inspect(mean([1.0, 3.0, 5.0]), content="3")
}
```
