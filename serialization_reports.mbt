///|
/// Human-readable report section.
pub struct ReportSection {
  title : String
  status : String
  lines : Array[String]
  fingerprint : UInt64
}

///|
/// CSV serialization policy.
pub struct CsvPolicy {
  delimiter : String
  decimal_digits : Int
  include_header : Bool
  missing_token : String
}

///|
/// Creates a CSV policy with safe defaults.
pub fn csv_policy(
  delimiter? : String = ",",
  decimal_digits? : Int = 6,
  include_header? : Bool = true,
  missing_token? : String = "",
) -> CsvPolicy {
  {
    delimiter,
    decimal_digits: decimal_digits.max(0).min(15),
    include_header,
    missing_token,
  }
}

///|
fn sr_number(value : Double, missing_token : String) -> String {
  if is_finite(value) {
    value.to_string()
  } else {
    missing_token
  }
}

///|
/// Serializes a numeric vector as one CSV row.
pub fn serialize_vector(values : Array[Double], policy : CsvPolicy) -> String {
  let builder = StringBuilder::new()
  for i in 0..<values.length() {
    if i > 0 {
      builder.write_string(policy.delimiter)
    }
    builder.write_string(sr_number(values[i], policy.missing_token))
  }
  builder.to_string()
}

///|
/// Serializes a numeric matrix as CSV.
pub fn serialize_matrix(
  matrix : Array[Array[Double]],
  policy : CsvPolicy,
  headers? : Array[String] = [],
) -> String {
  let builder = StringBuilder::new()
  if policy.include_header && headers.length() > 0 {
    for i in 0..<headers.length() {
      if i > 0 {
        builder.write_string(policy.delimiter)
      }
      builder.write_string(headers[i])
    }
    builder.write_string("\n")
  }
  for row in matrix {
    builder.write_string(serialize_vector(row, policy))
    builder.write_string("\n")
  }
  builder.to_string()
}

///|
/// Serializes an estimate into a stable key-value line.
pub fn serialize_estimate(estimate : Estimate) -> String {
  let builder = StringBuilder::new()
  builder.write_string("estimand=")
  builder.write_string(estimate.estimand)
  builder.write_string(",estimate=")
  builder.write_string(estimate.estimate.to_string())
  builder.write_string(",standard_error=")
  builder.write_string(estimate.standard_error.to_string())
  builder.write_string(",lower=")
  builder.write_string(estimate.lower.to_string())
  builder.write_string(",upper=")
  builder.write_string(estimate.upper.to_string())
  builder.to_string()
}

///|
/// Creates a quality report section.
pub fn quality_report_section(quality : DatasetQuality) -> ReportSection {
  let lines : Array[String] = Array::new()
  lines.push("rows=\{quality.rows}")
  lines.push("columns=\{quality.columns}")
  lines.push("complete_rows=\{quality.valid_rows}")
  lines.push("missing_cells=\{quality.missing_cells}")
  lines.push("score=\{quality.score}")
  let rows : Array[Array[Double]] = [quality_summary_vector(quality)]
  {
    title: "data-quality",
    status: if quality_gate(quality) {
      "passed"
    } else {
      "warning"
    },
    lines,
    fingerprint: matrix_checksum(rows),
  }
}

///|
/// Creates a pipeline report section.
pub fn pipeline_report_section(result : PipelineResult) -> ReportSection {
  let lines : Array[String] = Array::new()
  lines.push("estimate=\{result.estimate.estimate}")
  lines.push("standard_error=\{result.estimate.standard_error}")
  lines.push("score=\{result.score}")
  lines.push("passes=\{result.passes}")
  {
    title: "causal-pipeline",
    status: if result.passes {
      "passed"
    } else {
      "warning"
    },
    lines,
    fingerprint: pipeline_fingerprint(result),
  }
}

///|
/// Serializes a report section.
pub fn serialize_report_section(section : ReportSection) -> String {
  let builder = StringBuilder::new()
  builder.write_string("[")
  builder.write_string(section.title)
  builder.write_string("] ")
  builder.write_string(section.status)
  builder.write_string("\n")
  for line in section.lines {
    builder.write_string(line)
    builder.write_string("\n")
  }
  builder.write_string("fingerprint=")
  builder.write_string(section.fingerprint.to_string())
  builder.to_string()
}

///|
/// Serializes a report composed of sections.
pub fn serialize_report(sections : Array[ReportSection]) -> String {
  let builder = StringBuilder::new()
  for section in sections {
    builder.write_string(serialize_report_section(section))
    builder.write_string("\n")
  }
  builder.to_string()
}

///|
/// Returns a compact report fingerprint.
pub fn report_fingerprint(sections : Array[ReportSection]) -> UInt64 {
  let rows : Array[Array[Double]] = Array::new()
  for section in sections {
    rows.push([
      section.lines.length().to_double(),
      section.fingerprint.to_double(),
    ])
  }
  matrix_checksum(rows)
}
