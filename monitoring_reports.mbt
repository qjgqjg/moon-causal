///|
/// Monitoring rule for a production causal metric.
pub struct MonitoringRule {
  name : String
  lower : Double
  upper : Double
  severity : Int
  enabled : Bool
}

///|
/// Evaluation of one monitoring rule.
pub struct MonitoringAlert {
  name : String
  value : Double
  triggered : Bool
  severity : Int
  direction : String
}

///|
/// Baseline comparison for an estimate or data-quality metric.
pub struct BaselineComparison {
  current : Double
  baseline : Double
  absolute_change : Double
  relative_change : Double
  z_score : Double
  drifted : Bool
}

///|
/// Audit entry for a reproducible analysis run.
pub struct AuditEntry {
  stage : String
  status : String
  value : Double
  fingerprint : UInt64
  message : String
}

///|
/// Monitoring report with alerts and audit trail.
pub struct MonitoringReport {
  run_id : String
  entries : Array[AuditEntry]
  alerts : Array[MonitoringAlert]
  score : Double
  passes : Bool
}

///|
/// Creates a monitoring rule.
pub fn monitoring_rule(
  name : String,
  lower : Double,
  upper : Double,
  severity? : Int = 1,
  enabled? : Bool = true,
) -> MonitoringRule {
  {
    name,
    lower: lower.min(upper),
    upper: upper.max(lower),
    severity: severity.max(1),
    enabled,
  }
}

///|
/// Evaluates one rule against a value.
pub fn evaluate_monitoring_rule(
  rule : MonitoringRule,
  value : Double,
) -> MonitoringAlert {
  let below = value < rule.lower
  let above = value > rule.upper
  let triggered = rule.enabled && (below || above)
  {
    name: rule.name,
    value,
    triggered,
    severity: rule.severity,
    direction: if below {
      "below-lower-bound"
    } else if above {
      "above-upper-bound"
    } else {
      "within-bounds"
    },
  }
}

///|
/// Evaluates multiple rules in declaration order.
pub fn evaluate_monitoring_rules(
  rules : Array[MonitoringRule],
  values : Array[Double],
) -> Array[MonitoringAlert] {
  let result : Array[MonitoringAlert] = Array::new(capacity=rules.length())
  for i in 0..<rules.length() {
    result.push(
      evaluate_monitoring_rule(
        rules[i],
        if i < values.length() {
          values[i]
        } else {
          0.0
        },
      ),
    )
  }
  result
}

///|
/// Compares a current metric with a historical baseline.
pub fn compare_baseline(
  current : Double,
  baseline : Double,
  standard_error : Double,
  threshold? : Double = 2.0,
) -> BaselineComparison {
  let absolute = current - baseline
  let relative = if baseline.abs() < 1.0e-12 {
    0.0
  } else {
    absolute / baseline.abs()
  }
  let z = if standard_error.abs() < 1.0e-12 {
    0.0
  } else {
    absolute / standard_error.abs()
  }
  {
    current,
    baseline,
    absolute_change: absolute,
    relative_change: relative,
    z_score: z,
    drifted: z.abs() > threshold,
  }
}

///|
/// Computes a weighted estimate-quality score.
pub fn estimate_monitoring_score(
  estimate : Estimate,
  overlap_fraction : Double,
  calibration_error : Double,
  quality_score : Double,
) -> Double {
  let precision = 1.0 / (1.0 + estimate.standard_error.abs())
  let overlap = clamp(overlap_fraction, 0.0, 1.0)
  let calibration = 1.0 - clamp(calibration_error, 0.0, 1.0)
  let quality = clamp(quality_score, 0.0, 1.0)
  clamp(
    0.35 * precision + 0.25 * overlap + 0.2 * calibration + 0.2 * quality,
    0.0,
    1.0,
  )
}

///|
/// Creates a deterministic audit entry.
pub fn audit_entry(
  stage : String,
  status : String,
  value : Double,
  fingerprint : UInt64,
  message? : String = "",
) -> AuditEntry {
  { stage, status, value, fingerprint, message }
}

///|
/// Builds a monitoring report from audit entries and alerts.
pub fn monitoring_report(
  run_id : String,
  entries : Array[AuditEntry],
  alerts : Array[MonitoringAlert],
  minimum_score? : Double = 0.8,
) -> MonitoringReport {
  let mut penalty = 0.0
  for alert in alerts {
    if alert.triggered {
      penalty += alert.severity.to_double() * 0.1
    }
  }
  let score = clamp(1.0 - penalty, 0.0, 1.0)
  {
    run_id,
    entries,
    alerts,
    score,
    passes: score >= minimum_score && entries.length() > 0,
  }
}

///|
/// Counts triggered alerts by severity.
pub fn alert_counts(alerts : Array[MonitoringAlert]) -> Array[Int] {
  let result = Array::make(4, 0)
  for alert in alerts {
    if alert.triggered {
      result[alert.severity.min(3)] += 1
    }
  }
  result
}

///|
/// Returns whether all critical alerts are clear.
pub fn critical_alerts_clear(
  alerts : Array[MonitoringAlert],
  critical_severity? : Int = 3,
) -> Bool {
  for alert in alerts {
    if alert.triggered && alert.severity >= critical_severity {
      return false
    }
  }
  true
}

///|
/// Computes a drift report from reference and current samples.
pub fn monitoring_drift(
  reference : Array[Double],
  current : Array[Double],
  psi_threshold? : Double = 0.2,
  ks_threshold? : Double = 0.1,
) -> Array[Double] {
  let psi = population_stability_index(reference, current)
  let ks = kolmogorov_smirnov_distance(reference, current)
  [
    psi,
    ks,
    if psi > psi_threshold {
      1.0
    } else {
      0.0
    },
    if ks > ks_threshold {
      1.0
    } else {
      0.0
    },
  ]
}

///|
/// Computes a rolling estimate trend and alert direction.
pub fn estimate_trend(
  values : Array[Double],
  window : Int,
  threshold : Double,
) -> Array[MonitoringAlert] {
  let result : Array[MonitoringAlert] = Array::new(capacity=values.length())
  let baseline = if values.length() == 0 { 0.0 } else { values[0] }
  let smooth = moving_average(values, window)
  for i in 0..<smooth.length() {
    let rule = monitoring_rule(
      "estimate-trend",
      -threshold,
      threshold,
      severity=1,
    )
    result.push(evaluate_monitoring_rule(rule, smooth[i] - baseline))
  }
  result
}

///|
/// Serializes a numeric monitoring vector as a compact CSV-like string.
pub fn monitoring_vector_text(values : Array[Double]) -> String {
  let builder = StringBuilder::new()
  for i in 0..<values.length() {
    if i > 0 {
      builder.write_string(",")
    }
    builder.write_string(values[i].to_string())
  }
  builder.to_string()
}

///|
/// Serializes alerts for log shipping.
pub fn monitoring_alert_text(alerts : Array[MonitoringAlert]) -> String {
  let builder = StringBuilder::new()
  for alert in alerts {
    builder.write_string(alert.name)
    builder.write_string("|")
    builder.write_string(alert.value.to_string())
    builder.write_string("|")
    builder.write_string(if alert.triggered { "triggered" } else { "clear" })
    builder.write_string("|")
    builder.write_string(alert.direction)
    builder.write_string("\n")
  }
  builder.to_string()
}

///|
/// Builds a compact monitoring summary vector.
pub fn monitoring_summary(report : MonitoringReport) -> Array[Double] {
  let counts = alert_counts(report.alerts)
  [
    report.entries.length().to_double(),
    report.alerts.length().to_double(),
    counts[1].to_double(),
    counts[2].to_double(),
    counts[3].to_double(),
    report.score,
    if report.passes {
      1.0
    } else {
      0.0
    },
  ]
}
