///|
/// Undirected or directed network represented by an adjacency matrix.
pub struct CausalNetwork {
  adjacency : Array[Array[Int]]
  nodes : Int
  directed : Bool
}

///|
/// Exposure mapping for one unit.
pub struct ExposureSummary {
  treated_neighbors : Int
  control_neighbors : Int
  total_neighbors : Int
  treated_fraction : Double
  exposure_level : Int
}

///|
/// Direct and spillover effect summary.
pub struct SpilloverEffect {
  direct_effect : Double
  spillover_effect : Double
  treated_exposure_count : Int
  control_exposure_count : Int
  standard_error : Double
  passes : Bool
}

///|
/// Network diagnostics.
pub struct NetworkAudit {
  nodes : Int
  edges : Int
  isolated_nodes : Int
  maximum_degree : Int
  average_degree : Double
  symmetric : Bool
  passes : Bool
}

///|
/// Creates a network from a square integer adjacency matrix.
pub fn causal_network(
  adjacency : Array[Array[Int]],
  directed? : Bool = false,
) -> CausalNetwork {
  let nodes = adjacency.length()
  let normalized : Array[Array[Int]] = Array::new(capacity=nodes)
  for row in adjacency {
    let copied = Array::make(nodes, 0)
    for i in 0..<nodes.min(row.length()) {
      copied[i] = if row[i] == 0 { 0 } else { 1 }
    }
    normalized.push(copied)
  }
  { adjacency: normalized, nodes, directed }
}

///|
/// Counts network edges.
pub fn network_edge_count(network : CausalNetwork) -> Int {
  let mut total = 0
  for row in network.adjacency {
    for value in row {
      if value != 0 {
        total += 1
      }
    }
  }
  if network.directed {
    total
  } else {
    total / 2
  }
}

///|
/// Returns node degrees.
pub fn network_degrees(network : CausalNetwork) -> Array[Int] {
  let result = Array::make(network.nodes, 0)
  for i in 0..<network.nodes {
    for value in network.adjacency[i] {
      if value != 0 {
        result[i] += 1
      }
    }
  }
  result
}

///|
/// Lists neighbors of a node.
pub fn network_neighbors(network : CausalNetwork, node : Int) -> Array[Int] {
  let result : Array[Int] = Array::new()
  if node < 0 || node >= network.nodes {
    return result
  }
  for i in 0..<network.nodes {
    if network.adjacency[node][i] != 0 {
      result.push(i)
    }
  }
  result
}

///|
/// Computes the exposure summary induced by a treatment assignment.
pub fn network_exposure(
  network : CausalNetwork,
  treatment : Array[Bool],
  node : Int,
) -> ExposureSummary {
  let neighbors = network_neighbors(network, node)
  let mut treated = 0
  for neighbor in neighbors {
    if neighbor < treatment.length() && treatment[neighbor] {
      treated += 1
    }
  }
  let total = neighbors.length()
  let fraction = if total == 0 {
    0.0
  } else {
    treated.to_double() / total.to_double()
  }
  {
    treated_neighbors: treated,
    control_neighbors: total - treated,
    total_neighbors: total,
    treated_fraction: fraction,
    exposure_level: if treated == 0 {
      0
    } else if treated == total {
      2
    } else {
      1
    },
  }
}

///|
/// Computes exposure summaries for every node.
pub fn network_exposures(
  network : CausalNetwork,
  treatment : Array[Bool],
) -> Array[ExposureSummary] {
  let result : Array[ExposureSummary] = Array::new(capacity=network.nodes)
  for node in 0..<network.nodes {
    result.push(network_exposure(network, treatment, node))
  }
  result
}

///|
/// Audits connectivity and adjacency symmetry.
pub fn audit_network(network : CausalNetwork) -> NetworkAudit {
  let degrees = network_degrees(network)
  let mut isolated = 0
  let mut maximum = 0
  let mut total = 0
  let mut symmetric = true
  for degree in degrees {
    if degree == 0 {
      isolated += 1
    }
    if degree > maximum {
      maximum = degree
    }
    total += degree
  }
  if !network.directed {
    for i in 0..<network.nodes {
      for j in 0..<network.nodes {
        if network.adjacency[i][j] != network.adjacency[j][i] {
          symmetric = false
        }
      }
    }
  }
  {
    nodes: network.nodes,
    edges: network_edge_count(network),
    isolated_nodes: isolated,
    maximum_degree: maximum,
    average_degree: if network.nodes == 0 {
      0.0
    } else {
      total.to_double() / network.nodes.to_double()
    },
    symmetric,
    passes: network.nodes > 0 && (network.directed || symmetric),
  }
}

///|
/// Computes a direct effect among units with no treated neighbors.
pub fn network_direct_effect(
  treatment : Array[Bool],
  outcome : Array[Double],
  exposures : Array[ExposureSummary],
) -> Double {
  let treated = Array::new()
  let control = Array::new()
  let n = treatment.length().min(outcome.length()).min(exposures.length())
  for i in 0..<n {
    if exposures[i].treated_neighbors == 0 {
      if treatment[i] {
        treated.push(outcome[i])
      } else {
        control.push(outcome[i])
      }
    }
  }
  mean_or(treated, 0.0) - mean_or(control, 0.0)
}

///|
/// Computes a spillover effect among untreated units by neighbor exposure.
pub fn network_spillover_effect(
  treatment : Array[Bool],
  outcome : Array[Double],
  exposures : Array[ExposureSummary],
) -> SpilloverEffect {
  let exposed = Array::new()
  let unexposed = Array::new()
  let direct_treated = Array::new()
  let direct_control = Array::new()
  let n = treatment.length().min(outcome.length()).min(exposures.length())
  for i in 0..<n {
    if treatment[i] {
      direct_treated.push(outcome[i])
    } else if exposures[i].treated_neighbors > 0 {
      exposed.push(outcome[i])
    } else {
      unexposed.push(outcome[i])
      direct_control.push(outcome[i])
    }
  }
  let direct = mean_or(direct_treated, 0.0) - mean_or(direct_control, 0.0)
  let spillover = mean_or(exposed, 0.0) - mean_or(unexposed, 0.0)
  let combined = Array::new(capacity=exposed.length() + unexposed.length())
  for value in exposed {
    combined.push(value)
  }
  for value in unexposed {
    combined.push(value)
  }
  let standard_error = if combined.length() < 2 {
    0.0
  } else {
    std_dev(combined) / combined.length().to_double().sqrt()
  }
  {
    direct_effect: direct,
    spillover_effect: spillover,
    treated_exposure_count: exposed.length(),
    control_exposure_count: unexposed.length(),
    standard_error,
    passes: exposed.length() > 1 && unexposed.length() > 1,
  }
}

///|
/// Returns integer exposure strata for network-aware adjustment.
pub fn network_exposure_strata(
  network : CausalNetwork,
  treatment : Array[Bool],
) -> Array[Int] {
  let exposures = network_exposures(network, treatment)
  let result = Array::new(capacity=exposures.length())
  for exposure in exposures {
    result.push(exposure.exposure_level)
  }
  result
}

///|
/// Assigns complete network clusters from a cluster label vector.
pub fn network_cluster_assignment(
  cluster_ids : Array[Int],
  treatment_fraction? : Double = 0.5,
  seed? : UInt64 = 20260819,
) -> Array[Bool] {
  let clusters : Array[Int] = Array::new()
  for cluster in cluster_ids {
    if !clusters.contains(cluster) {
      clusters.push(cluster)
    }
  }
  let assignment = blocked_assignment(
    clusters.length(),
    block_size=2,
    treatment_fraction~,
    seed~,
  )
  let result = Array::make(cluster_ids.length(), false)
  for i in 0..<cluster_ids.length() {
    for j in 0..<clusters.length() {
      if clusters[j] == cluster_ids[i] {
        result[i] = assignment.treatment[j]
      }
    }
  }
  result
}

///|
/// Computes a thresholded network exposure for policy simulation.
pub fn network_threshold_exposure(
  network : CausalNetwork,
  treatment : Array[Bool],
  threshold : Double,
) -> Array[Bool] {
  let exposures = network_exposures(network, treatment)
  let result = Array::new(capacity=exposures.length())
  for exposure in exposures {
    result.push(exposure.treated_fraction >= threshold)
  }
  result
}

///|
/// Returns connected components under an undirected view.
pub fn network_components(network : CausalNetwork) -> Array[Array[Int]] {
  let visited = Array::make(network.nodes, false)
  let result : Array[Array[Int]] = Array::new()
  for start in 0..<network.nodes {
    if !visited[start] {
      let component : Array[Int] = Array::new()
      let queue : Array[Int] = [start]
      visited[start] = true
      let mut cursor = 0
      while cursor < queue.length() {
        let node = queue[cursor]
        component.push(node)
        for neighbor in network_neighbors(network, node) {
          if !visited[neighbor] {
            visited[neighbor] = true
            queue.push(neighbor)
          }
        }
        cursor += 1
      }
      result.push(component)
    }
  }
  result
}

///|
/// Computes a compact network summary vector.
pub fn network_summary(network : CausalNetwork) -> Array[Double] {
  let audit = audit_network(network)
  [
    audit.nodes.to_double(),
    audit.edges.to_double(),
    audit.isolated_nodes.to_double(),
    audit.maximum_degree.to_double(),
    audit.average_degree,
    if audit.symmetric {
      1.0
    } else {
      0.0
    },
    if audit.passes {
      1.0
    } else {
      0.0
    },
  ]
}
