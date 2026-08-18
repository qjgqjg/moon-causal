///|
/// Calculates Standardized Mean Difference (SMD) for a single covariate between treatment and control groups.
pub fn smd(covariate : Array[Double], treatments : Array[Bool]) -> Double {
  let treated = Array::new()
  let control = Array::new()
  let n = if covariate.length() < treatments.length() {
    covariate.length()
  } else {
    treatments.length()
  }
  for i in 0..<n {
    if treatments[i] {
      treated.push(covariate[i])
    } else {
      control.push(covariate[i])
    }
  }
  let mean_t = mean(treated)
  let mean_c = mean(control)
  let var_t = variance(treated, sample=true)
  let var_c = variance(control, sample=true)
  let pooled_sd = ((var_t + var_c) / 2.0).sqrt()
  if pooled_sd == 0.0 {
    return 0.0
  }
  (mean_t - mean_c) / pooled_sd
}
