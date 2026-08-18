///|
pub struct GraphEdge {
  source : String
  target : String
}

///|
pub struct CausalGraph {
  nodes : Array[String]
  edges : Array[GraphEdge]
}

///|
pub fn CausalGraph::new(nodes : Array[String]) -> CausalGraph {
  { nodes, edges: [] }
}

///|
pub fn CausalGraph::add_edge(
  self : CausalGraph,
  source : String,
  target : String,
) -> CausalGraph {
  if self.has_edge(source, target) || source == target {
    return self
  }
  let edges = self.edges.copy()
  edges.push({ source, target })
  { nodes: self.nodes, edges }
}

///|
pub fn CausalGraph::has_edge(
  self : CausalGraph,
  source : String,
  target : String,
) -> Bool {
  for edge in self.edges {
    if edge.source == source && edge.target == target {
      return true
    }
  }
  false
}

///|
pub fn CausalGraph::edge_count(self : CausalGraph) -> Int {
  self.edges.length()
}

///|
fn node_index(nodes : Array[String], target : String) -> Int {
  for i in 0..<nodes.length() {
    if nodes[i] == target {
      return i
    }
  }
  -1
}

///|
/// Returns a topological ordering, or an empty array when a cycle is detected.
pub fn CausalGraph::topological_order(self : CausalGraph) -> Array[String] {
  let indegree = Array::make(self.nodes.length(), 0)
  for edge in self.edges {
    let target = node_index(self.nodes, edge.target)
    if target >= 0 {
      indegree[target] += 1
    }
  }
  let ready = Array::new()
  for i in 0..<indegree.length() {
    if indegree[i] == 0 {
      ready.push(i)
    }
  }
  let order : Array[String] = Array::new(capacity=self.nodes.length())
  let mut cursor = 0
  while cursor < ready.length() {
    let current = ready[cursor]
    cursor += 1
    order.push(self.nodes[current])
    for edge in self.edges {
      if edge.source == self.nodes[current] {
        let target = node_index(self.nodes, edge.target)
        if target >= 0 {
          indegree[target] -= 1
          if indegree[target] == 0 {
            ready.push(target)
          }
        }
      }
    }
  }
  if order.length() == self.nodes.length() {
    order
  } else {
    []
  }
}

///|
/// Returns direct parents of a node.
pub fn CausalGraph::parents(
  self : CausalGraph,
  target : String,
) -> Array[String] {
  let result = Array::new()
  for edge in self.edges {
    if edge.target == target {
      result.push(edge.source)
    }
  }
  result
}

///|
/// Returns direct children of a node.
pub fn CausalGraph::children(
  self : CausalGraph,
  source : String,
) -> Array[String] {
  let result = Array::new()
  for edge in self.edges {
    if edge.source == source {
      result.push(edge.target)
    }
  }
  result
}

///|
/// Returns all ancestors of a node using a breadth-first traversal.
pub fn CausalGraph::ancestors(
  self : CausalGraph,
  target : String,
) -> Array[String] {
  let result : Array[String] = Array::new()
  let queue = self.parents(target)
  let mut cursor = 0
  while cursor < queue.length() {
    let current = queue[cursor]
    cursor += 1
    if !result.contains(current) {
      result.push(current)
      for parent in self.parents(current) {
        if !result.contains(parent) {
          queue.push(parent)
        }
      }
    }
  }
  result
}

///|
/// Lists treatment parents that can be used as a simple backdoor adjustment set.
pub fn CausalGraph::backdoor_candidates(
  self : CausalGraph,
  treatment : String,
  outcome : String,
) -> Array[String] {
  let result = Array::new()
  for parent in self.parents(treatment) {
    if parent != outcome {
      result.push(parent)
    }
  }
  result
}
