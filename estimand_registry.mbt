///|
/// A documented estimand definition for an analysis registry.
pub struct EstimandDefinition {
  name : String
  description : String
  treatment : String
  outcome : String
  population : String
  identification : String
  required_diagnostics : Array[String]
  version : String
}

///|
/// Registry of named estimands used by a project.
pub struct EstimandRegistry {
  definitions : Array[EstimandDefinition]
  active : String
  version : String
}

///|
/// Validation result for an estimand registry.
pub struct EstimandAudit {
  registered : Int
  duplicate_names : Int
  missing_fields : Int
  active_found : Bool
  passes : Bool
}

///|
/// Creates a named estimand definition.
pub fn estimand_definition(
  name : String,
  description : String,
  treatment : String,
  outcome : String,
  population? : String = "all",
  identification? : String = "backdoor",
  required_diagnostics? : Array[String] = [],
  version? : String = "1",
) -> EstimandDefinition {
  {
    name,
    description,
    treatment,
    outcome,
    population,
    identification,
    required_diagnostics,
    version,
  }
}

///|
/// Creates an empty registry.
pub fn estimand_registry(version? : String = "1") -> EstimandRegistry {
  { definitions: [], active: "", version }
}

///|
/// Registers or replaces an estimand by name.
pub fn register_estimand(
  registry : EstimandRegistry,
  definition : EstimandDefinition,
  activate? : Bool = false,
) -> EstimandRegistry {
  let definitions = registry.definitions.copy()
  let mut found = -1
  for i in 0..<definitions.length() {
    if definitions[i].name == definition.name {
      found = i
      break
    }
  }
  if found >= 0 {
    definitions[found] = definition
  } else {
    definitions.push(definition)
  }
  {
    definitions,
    active: if activate || registry.active == "" {
      definition.name
    } else {
      registry.active
    },
    version: registry.version,
  }
}

///|
/// Looks up an estimand by name.
pub fn find_estimand(
  registry : EstimandRegistry,
  name : String,
) -> EstimandDefinition? {
  for definition in registry.definitions {
    if definition.name == name {
      return Some(definition)
    }
  }
  None
}

///|
/// Sets the active estimand when the name exists.
pub fn activate_estimand(
  registry : EstimandRegistry,
  name : String,
) -> EstimandRegistry {
  if find_estimand(registry, name) is Some(_) {
    {
      definitions: registry.definitions,
      active: name,
      version: registry.version,
    }
  } else {
    registry
  }
}

///|
/// Counts duplicate estimand names.
pub fn duplicate_estimand_names(registry : EstimandRegistry) -> Int {
  let mut duplicates = 0
  for i in 0..<registry.definitions.length() {
    for j in 0..<i {
      if registry.definitions[i].name == registry.definitions[j].name {
        duplicates += 1
        break
      }
    }
  }
  duplicates
}

///|
/// Audits registry completeness and active selection.
pub fn audit_estimand_registry(registry : EstimandRegistry) -> EstimandAudit {
  let mut missing = 0
  for definition in registry.definitions {
    if definition.name == "" {
      missing += 1
    }
    if definition.treatment == "" {
      missing += 1
    }
    if definition.outcome == "" {
      missing += 1
    }
    if definition.identification == "" {
      missing += 1
    }
  }
  let active_found = registry.active == "" ||
    find_estimand(registry, registry.active) is Some(_)
  let duplicates = duplicate_estimand_names(registry)
  {
    registered: registry.definitions.length(),
    duplicate_names: duplicates,
    missing_fields: missing,
    active_found,
    passes: duplicates == 0 &&
    missing == 0 &&
    active_found &&
    registry.definitions.length() > 0,
  }
}

///|
/// Returns the active estimand, if configured.
pub fn active_estimand(registry : EstimandRegistry) -> EstimandDefinition? {
  if registry.active == "" {
    None
  } else {
    find_estimand(registry, registry.active)
  }
}

///|
/// Returns all registry names in declaration order.
pub fn estimand_names(registry : EstimandRegistry) -> Array[String] {
  let result : Array[String] = Array::new(
    capacity=registry.definitions.length(),
  )
  for definition in registry.definitions {
    result.push(definition.name)
  }
  result
}

///|
/// Checks whether a required diagnostic is registered for an estimand.
pub fn requires_diagnostic(
  definition : EstimandDefinition,
  diagnostic : String,
) -> Bool {
  definition.required_diagnostics.contains(diagnostic)
}

///|
/// Adds a required diagnostic without duplicating it.
pub fn add_required_diagnostic(
  definition : EstimandDefinition,
  diagnostic : String,
) -> EstimandDefinition {
  let diagnostics = definition.required_diagnostics.copy()
  if !diagnostics.contains(diagnostic) {
    diagnostics.push(diagnostic)
  }
  {
    name: definition.name,
    description: definition.description,
    treatment: definition.treatment,
    outcome: definition.outcome,
    population: definition.population,
    identification: definition.identification,
    required_diagnostics: diagnostics,
    version: definition.version,
  }
}

///|
/// Computes a registry fingerprint from definition metadata.
pub fn estimand_registry_fingerprint(registry : EstimandRegistry) -> UInt64 {
  let rows : Array[Array[Double]] = Array::new()
  for definition in registry.definitions {
    rows.push([
      definition.name.length().to_double(),
      definition.treatment.length().to_double(),
      definition.outcome.length().to_double(),
      definition.required_diagnostics.length().to_double(),
      definition.version.length().to_double(),
    ])
  }
  matrix_checksum(rows)
}

///|
/// Creates a compact registry summary vector.
pub fn estimand_registry_summary(registry : EstimandRegistry) -> Array[Double] {
  let audit = audit_estimand_registry(registry)
  [
    audit.registered.to_double(),
    audit.duplicate_names.to_double(),
    audit.missing_fields.to_double(),
    if audit.active_found {
      1.0
    } else {
      0.0
    },
    if audit.passes {
      1.0
    } else {
      0.0
    },
    estimand_registry_fingerprint(registry).to_double(),
  ]
}
