///|
/// Result of validating a causal graph.
pub struct GraphAudit {
  node_count : Int
  edge_count : Int
  duplicate_nodes : Int
  dangling_edges : Int
  self_edges : Int
  cyclic : Bool
  connected_components : Int
  passes : Bool
}

///|
/// A graph path represented by node names.
pub struct GraphPath {
  nodes : Array[String]
  directed : Bool
  contains_treatment : Bool
  contains_outcome : Bool
}

///|
/// Summary of d-separation or path enumeration.
pub struct SeparationResult {
  separated : Bool
  active_path_count : Int
  blocked_path_count : Int
  adjustment_candidates : Array[String]
  explanation : String
}

///|
/// A layer in a topological causal ordering.
pub struct GraphLayer {
  depth : Int
  nodes : Array[String]
}

///|
/// A compact edge table for serialization and auditing.
pub struct GraphEdgeTable {
  sources : Array[String]
  targets : Array[String]
  edge_count : Int
  fingerprint : UInt64
}

///|
fn advanced_node_index(nodes : Array[String], name : String) -> Int {
  for i in 0..<nodes.length() {
    if nodes[i] == name {
      return i
    }
  }
  -1
}

///|
fn advanced_unique_count(values : Array[String]) -> Int {
  let unique : Array[String] = Array::new()
  for value in values {
    if !unique.contains(value) {
      unique.push(value)
    }
  }
  unique.length()
}

///|

///|
/// Counts duplicate names in a graph declaration.
pub fn graph_duplicate_node_count(graph : CausalGraph) -> Int {
  graph.nodes.length() - advanced_unique_count(graph.nodes)
}

///|
/// Counts edges that refer to undeclared nodes.
pub fn graph_dangling_edge_count(graph : CausalGraph) -> Int {
  let mut result = 0
  for edge in graph.edges {
    if advanced_node_index(graph.nodes, edge.source) < 0 ||
      advanced_node_index(graph.nodes, edge.target) < 0 {
      result += 1
    }
  }
  result
}

///|
/// Counts self loops in a graph.
pub fn graph_self_edge_count(graph : CausalGraph) -> Int {
  let mut result = 0
  for edge in graph.edges {
    if edge.source == edge.target {
      result += 1
    }
  }
  result
}

///|
/// Detects directed cycles using depth-first color marking.
pub fn graph_has_cycle(graph : CausalGraph) -> Bool {
  let state = Array::make(graph.nodes.length(), 0)
  fn visit(graph : CausalGraph, state : Array[Int], index : Int) -> Bool {
    if state[index] == 1 {
      return true
    }
    if state[index] == 2 {
      return false
    }
    state[index] = 1
    let name = graph.nodes[index]
    for edge in graph.edges {
      if edge.source == name {
        let target = advanced_node_index(graph.nodes, edge.target)
        if target >= 0 && visit(graph, state, target) {
          return true
        }
      }
    }
    state[index] = 2
    false
  }
  for i in 0..<graph.nodes.length() {
    if visit(graph, state, i) {
      return true
    }
  }
  false
}

///|
/// Counts undirected connected components for graph diagnostics.
pub fn graph_component_count(graph : CausalGraph) -> Int {
  let visited = Array::make(graph.nodes.length(), false)
  let mut components = 0
  for start in 0..<graph.nodes.length() {
    if !visited[start] {
      components += 1
      let queue : Array[Int] = [start]
      visited[start] = true
      let mut cursor = 0
      while cursor < queue.length() {
        let current = graph.nodes[queue[cursor]]
        for edge in graph.edges {
          let neighbor_name = if edge.source == current {
            edge.target
          } else if edge.target == current {
            edge.source
          } else {
            ""
          }
          if neighbor_name != "" {
            let neighbor = advanced_node_index(graph.nodes, neighbor_name)
            if neighbor >= 0 && !visited[neighbor] {
              visited[neighbor] = true
              queue.push(neighbor)
            }
          }
        }
        cursor += 1
      }
    }
  }
  components
}

///|
/// Audits node declarations and directed edges.
pub fn audit_causal_graph(graph : CausalGraph) -> GraphAudit {
  let duplicate_nodes = graph_duplicate_node_count(graph)
  let dangling_edges = graph_dangling_edge_count(graph)
  let self_edges = graph_self_edge_count(graph)
  let cyclic = graph_has_cycle(graph)
  let components = graph_component_count(graph)
  {
    node_count: graph.nodes.length(),
    edge_count: graph.edges.length(),
    duplicate_nodes,
    dangling_edges,
    self_edges,
    cyclic,
    connected_components: components,
    passes: duplicate_nodes == 0 &&
    dangling_edges == 0 &&
    self_edges == 0 &&
    !cyclic,
  }
}

///|
/// Returns nodes with no incoming edges.
pub fn graph_roots(graph : CausalGraph) -> Array[String] {
  let result : Array[String] = Array::new()
  for node in graph.nodes {
    if graph.parents(node).length() == 0 {
      result.push(node)
    }
  }
  result
}

///|
/// Returns nodes with no outgoing edges.
pub fn graph_leaves(graph : CausalGraph) -> Array[String] {
  let result : Array[String] = Array::new()
  for node in graph.nodes {
    if graph.children(node).length() == 0 {
      result.push(node)
    }
  }
  result
}

///|
/// Returns descendants including the starting node when requested.
pub fn graph_descendants(
  graph : CausalGraph,
  node : String,
  include_self? : Bool = false,
) -> Array[String] {
  let result : Array[String] = Array::new()
  let queue : Array[String] = [node]
  let mut cursor = 0
  while cursor < queue.length() {
    let current = queue[cursor]
    if (include_self || current != node) && !result.contains(current) {
      result.push(current)
    }
    for child in graph.children(current) {
      if !queue.contains(child) {
        queue.push(child)
      }
    }
    cursor += 1
  }
  result
}

///|
/// Returns ancestors including the starting node when requested.
pub fn graph_ancestors(
  graph : CausalGraph,
  node : String,
  include_self? : Bool = false,
) -> Array[String] {
  let result : Array[String] = Array::new()
  let queue : Array[String] = [node]
  let mut cursor = 0
  while cursor < queue.length() {
    let current = queue[cursor]
    if (include_self || current != node) && !result.contains(current) {
      result.push(current)
    }
    for parent in graph.parents(current) {
      if !queue.contains(parent) {
        queue.push(parent)
      }
    }
    cursor += 1
  }
  result
}

///|
/// Returns whether a directed path exists between two nodes.
pub fn graph_reaches(
  graph : CausalGraph,
  source : String,
  target : String,
) -> Bool {
  if source == target {
    return true
  }
  let visited : Array[String] = Array::new()
  let queue : Array[String] = [source]
  let mut cursor = 0
  while cursor < queue.length() {
    let current = queue[cursor]
    if current == target {
      return true
    }
    if !visited.contains(current) {
      visited.push(current)
      for child in graph.children(current) {
        queue.push(child)
      }
    }
    cursor += 1
  }
  false
}

///|
/// Returns all direct mediators between treatment and outcome.
pub fn graph_mediators(
  graph : CausalGraph,
  treatment : String,
  outcome : String,
) -> Array[String] {
  let result : Array[String] = Array::new()
  for candidate in graph.children(treatment) {
    if candidate != outcome && graph_reaches(graph, candidate, outcome) {
      result.push(candidate)
    }
  }
  result
}

///|
/// Returns descendants of treatment that are not on a directed treatment-outcome path.
pub fn graph_post_treatment_confounders(
  graph : CausalGraph,
  treatment : String,
  outcome : String,
) -> Array[String] {
  let descendants = graph_descendants(graph, treatment)
  let result : Array[String] = Array::new()
  for candidate in descendants {
    if candidate != outcome && !graph_reaches(graph, candidate, outcome) {
      result.push(candidate)
    }
  }
  result
}

///|
/// Returns a conservative set of pre-treatment backdoor candidates.
pub fn graph_adjustment_candidates(
  graph : CausalGraph,
  treatment : String,
  outcome : String,
) -> Array[String] {
  let result : Array[String] = Array::new()
  let treatment_parents = graph.parents(treatment)
  for parent in treatment_parents {
    if parent != outcome && !graph_reaches(graph, parent, outcome) {
      result.push(parent)
    }
    for ancestor in graph_ancestors(graph, parent) {
      if ancestor != treatment &&
        ancestor != outcome &&
        !result.contains(ancestor) {
        result.push(ancestor)
      }
    }
  }
  for candidate in graph.backdoor_candidates(treatment, outcome) {
    if candidate != outcome && !result.contains(candidate) {
      result.push(candidate)
    }
  }
  result
}

///|
/// Finds one adjustment set by removing post-treatment variables.
pub fn graph_backdoor_adjustment_set(
  graph : CausalGraph,
  treatment : String,
  outcome : String,
) -> Array[String] {
  let candidates = graph_adjustment_candidates(graph, treatment, outcome)
  let post_treatment = graph_descendants(graph, treatment)
  let result : Array[String] = Array::new()
  for candidate in candidates {
    if !post_treatment.contains(candidate) &&
      candidate != treatment &&
      candidate != outcome {
      result.push(candidate)
    }
  }
  result
}

///|
/// Returns whether a candidate is safe under a conservative backdoor rule.
pub fn graph_is_adjustment_candidate(
  graph : CausalGraph,
  treatment : String,
  outcome : String,
  candidate : String,
) -> Bool {
  candidate != treatment &&
  candidate != outcome &&
  !graph_descendants(graph, treatment).contains(candidate) &&
  graph_adjustment_candidates(graph, treatment, outcome).contains(candidate)
}

///|
/// Computes graph node depths from roots.
pub fn graph_depths(graph : CausalGraph) -> Array[Int] {
  let depths = Array::make(graph.nodes.length(), 0)
  let order = graph.topological_order()
  for name in order {
    let index = advanced_node_index(graph.nodes, name)
    let mut depth = 0
    for parent in graph.parents(name) {
      let parent_index = advanced_node_index(graph.nodes, parent)
      if parent_index >= 0 && depths[parent_index] + 1 > depth {
        depth = depths[parent_index] + 1
      }
    }
    if index >= 0 {
      depths[index] = depth
    }
  }
  depths
}

///|
/// Groups a DAG into topological layers.
pub fn graph_layers(graph : CausalGraph) -> Array[GraphLayer] {
  let depths = graph_depths(graph)
  let result : Array[GraphLayer] = Array::new()
  let mut maximum = 0
  for depth in depths {
    if depth > maximum {
      maximum = depth
    }
  }
  for depth in 0..<=maximum {
    let names : Array[String] = Array::new()
    for i in 0..<graph.nodes.length() {
      if depths[i] == depth {
        names.push(graph.nodes[i])
      }
    }
    if names.length() > 0 {
      result.push({ depth, nodes: names })
    }
  }
  result
}

///|
/// Returns the number of paths up to a maximum depth.
pub fn graph_path_count(
  graph : CausalGraph,
  source : String,
  target : String,
  maximum_depth? : Int = 20,
) -> Int {
  fn count_paths(
    graph : CausalGraph,
    current : String,
    target : String,
    path : Array[String],
    depth : Int,
    limit : Int,
  ) -> Int {
    if current == target {
      return 1
    }
    if depth >= limit {
      return 0
    }
    let mut count = 0
    for child in graph.children(current) {
      if !path.contains(child) {
        let next_path = path.copy()
        next_path.push(child)
        count += count_paths(graph, child, target, next_path, depth + 1, limit)
      }
    }
    count
  }
  count_paths(
    graph,
    source,
    target,
    [source],
    0,
    if maximum_depth > 0 {
      maximum_depth
    } else {
      20
    },
  )
}

///|
/// Enumerates directed paths up to a bounded depth.
pub fn graph_paths(
  graph : CausalGraph,
  source : String,
  target : String,
  maximum_depth? : Int = 8,
) -> Array[GraphPath] {
  let result : Array[GraphPath] = Array::new()
  fn walk(
    graph : CausalGraph,
    current : String,
    target : String,
    path : Array[String],
    result : Array[GraphPath],
    depth : Int,
    limit : Int,
  ) -> Unit {
    if current == target {
      result.push({
        nodes: path,
        directed: true,
        contains_treatment: path.contains(current),
        contains_outcome: true,
      })
      return
    }
    if depth >= limit {
      return
    }
    for child in graph.children(current) {
      if !path.contains(child) {
        let next = path.copy()
        next.push(child)
        walk(graph, child, target, next, result, depth + 1, limit)
      }
    }
  }
  walk(
    graph,
    source,
    target,
    [source],
    result,
    0,
    if maximum_depth > 0 {
      maximum_depth
    } else {
      8
    },
  )
  result
}

///|
/// Returns a simple active-path summary using a supplied conditioning set.
pub fn graph_separation(
  graph : CausalGraph,
  source : String,
  target : String,
  conditioned : Array[String],
) -> SeparationResult {
  let paths = graph_paths(graph, source, target)
  let active = Array::new()
  let mut blocked = 0
  for path in paths {
    let has_conditioned = path.nodes.any(fn(node) { conditioned.contains(node) })
    if has_conditioned {
      blocked += 1
    } else {
      active.push(path)
    }
  }
  {
    separated: active.length() == 0,
    active_path_count: active.length(),
    blocked_path_count: blocked,
    adjustment_candidates: graph_adjustment_candidates(graph, source, target),
    explanation: if active.length() == 0 {
      "all enumerated directed paths are conditioned"
    } else {
      "at least one unconditioned directed path remains"
    },
  }
}

///|
/// Returns the total number of incoming and outgoing edges for each node.
pub fn graph_degree_table(graph : CausalGraph) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=graph.nodes.length())
  for node in graph.nodes {
    result.push([
      graph.parents(node).length().to_double(),
      graph.children(node).length().to_double(),
    ])
  }
  result
}

///|
/// Creates an edge table in declaration order.
pub fn graph_edge_table(graph : CausalGraph) -> GraphEdgeTable {
  let sources : Array[String] = Array::new()
  let targets : Array[String] = Array::new()
  let fingerprint_rows : Array[Array[Double]] = Array::new()
  for edge in graph.edges {
    sources.push(edge.source)
    targets.push(edge.target)
    fingerprint_rows.push([
      advanced_node_index(graph.nodes, edge.source).to_double(),
      advanced_node_index(graph.nodes, edge.target).to_double(),
    ])
  }
  {
    sources,
    targets,
    edge_count: sources.length(),
    fingerprint: matrix_checksum(fingerprint_rows),
  }
}

///|
/// Returns a stable adjacency matrix ordered by graph nodes.
pub fn graph_adjacency_matrix(graph : CausalGraph) -> Array[Array[Double]] {
  let matrix = Array::make(
    graph.nodes.length(),
    Array::make(graph.nodes.length(), 0.0),
  )
  for edge in graph.edges {
    let source = advanced_node_index(graph.nodes, edge.source)
    let target = advanced_node_index(graph.nodes, edge.target)
    if source >= 0 && target >= 0 {
      matrix[source][target] = 1.0
    }
  }
  matrix
}

///|
/// Reconstructs a graph from a node list and edge table.
pub fn graph_from_edge_table(
  nodes : Array[String],
  sources : Array[String],
  targets : Array[String],
) -> CausalGraph {
  let graph = CausalGraph::new(nodes)
  let n = if sources.length() < targets.length() {
    sources.length()
  } else {
    targets.length()
  }
  let mut result = graph
  for i in 0..<n {
    result = result.add_edge(sources[i], targets[i])
  }
  result
}

///|
/// Removes variables from a graph while retaining the remaining direct edges.
pub fn graph_remove_nodes(
  graph : CausalGraph,
  removed : Array[String],
) -> CausalGraph {
  let nodes = Array::new()
  for node in graph.nodes {
    if !removed.contains(node) {
      nodes.push(node)
    }
  }
  let result = CausalGraph::new(nodes)
  let mut rebuilt = result
  for edge in graph.edges {
    if !removed.contains(edge.source) && !removed.contains(edge.target) {
      rebuilt = rebuilt.add_edge(edge.source, edge.target)
    }
  }
  rebuilt
}

///|
/// Returns the induced ancestral subgraph of selected nodes.
pub fn graph_ancestral_subgraph(
  graph : CausalGraph,
  selected : Array[String],
) -> CausalGraph {
  let keep = selected.copy()
  for node in selected {
    for ancestor in graph_ancestors(graph, node) {
      if !keep.contains(ancestor) {
        keep.push(ancestor)
      }
    }
  }
  let removed = Array::new()
  for node in graph.nodes {
    if !keep.contains(node) {
      removed.push(node)
    }
  }
  graph_remove_nodes(graph, removed)
}

///|
/// Produces a compact graph diagnostic vector.
pub fn graph_summary_vector(graph : CausalGraph) -> Array[Double] {
  let audit = audit_causal_graph(graph)
  [
    audit.node_count.to_double(),
    audit.edge_count.to_double(),
    audit.duplicate_nodes.to_double(),
    audit.dangling_edges.to_double(),
    audit.self_edges.to_double(),
    if audit.cyclic {
      1.0
    } else {
      0.0
    },
    audit.connected_components.to_double(),
    if audit.passes {
      1.0
    } else {
      0.0
    },
  ]
}
