---

# MacroIR Indexed Evaluation and `IndexedResults`

**Status:** DRAFT

**Scope:** DSL design note for indexed expansion and downstream propagation

**Derived from:** `grammar_untyped.h`, `lexer_typed.h`, `grammar_typed.h`, and design discussion

**Last updated:** 2026-03-30

---

## 1. Purpose

This document specifies a proposed indexed extension of the MacroIR DSL.

The goal is to allow some function arguments to act as batch sources while preserving:

- typed compilation
- explicit axis semantics
- downstream propagation of indexed results
- compatibility with existing scalar and vector payload types

The design is intentionally close in spirit to the indexed semantics previously explored in
`quimulun`, but adapted to the current MacroIR DSL and runtime.

---

## 2. Core Idea

An expression may evaluate either:

- as an ordinary value `T`
- or as an indexed value varying over one or more **root axes**

When a function argument expects `T` but receives a `vector<T>` expression, the vector is
interpreted as introducing a new root axis and the function is evaluated pointwise over that axis.

When a function argument expects `vector<T>`, the same expression is interpreted as ordinary data
and introduces no axis.

Therefore indexed expansion is **contextual**: it depends on the expected type at the call site.

---

## 3. Terminology

The following names are the preferred vocabulary for the indexed DSL extension.

### 3.1 Axis Terms

- `axis`
  Independent semantic dimension of variation.
  Examples: `model`, `algorithm`, `seed`.

- `root_axis`
  An axis introduced directly by an expandable argument source.
  Root axes are independent and combine by Cartesian product.

- `dependent_axis`
  An axis whose domain depends on the coordinate of one or more parent axes.
  A dependent axis does not create an independent Cartesian dimension.

- `axis_id`
  Unique internal identity of an axis.
  This is not the same thing as the human-readable label.

- `axis_label`
  Human-readable label associated with an axis.

- `axis_domain<T>`
  The runtime values associated with one axis.

- `local_axis_domain<T>`
  The values associated with a dependent axis at a particular parent coordinate.

- `axis_index`
  Integer position inside one axis domain.

### 3.2 Coordinate Terms

- `coordinate`
  One point in the product of several root axes.

- `index_space`
  The set of axes on which an expression depends.
  For rectangular cases this is just the product of root axes.
  For dependent-axis cases this may be ragged.

### 3.3 Indexed Value Terms

- `IndexedResults<T>`
  First-class DSL value representing a payload of type `T` evaluated over an `index_space`.

- `indexed payload`
  The `T` carried at each coordinate of an `IndexedResults<T>`.

### 3.4 Expression Operations

- `collect_root_axes(env)`
  Collect the root axes on which an expression depends.

- `eval_at(env, coordinate)`
  Evaluate an expression at one coordinate in an index space.

- `materialize_indexed(expr, env)`
  Evaluate an expression over its collected root axes and produce either a scalar `T` or an
  `IndexedResults<T>`.

---

## 4. Contextual Meaning of `vector<T>`

The indexed extension does not give vector literals a single fixed meaning.

Instead:

### Rule A: Expected `vector<T>`

If a function parameter expects `vector<T>`, then an actual vector expression is compiled as:

```cpp
typed_expression<vector<T>>
```

Its `index_space` is empty.

This is an ordinary vector payload.

### Rule B: Expected `T`

If a function parameter expects `T`, but the actual expression is a vector expression convertible
to `vector<T>`, then the expression is compiled as an indexed scalar expression:

```cpp
typed_expression<T>
```

with a newly introduced `root_axis`.

The vector elements become the domain of that root axis.

### Consequence

The same syntax may be interpreted either:

- as a payload vector
- or as a source of indexed expansion

depending on the formal parameter type.

This is intentional and is considered the core rule of contextual indexed evaluation.

---

## 5. Root Axes, Zip, Cartesian Product

### 5.1 Root Axes

A standalone expandable source introduces a fresh `root_axis`.

Examples:

- `models = [m1, m2]`
- `algorithms = [a1, a2, a3]`

If these are later passed where scalar `T` is expected, they become batch sources.

### 5.2 Cartesian Product

Different root axes combine by Cartesian product.

If one argument depends on root axis `model` and another depends on root axis `algorithm`, then a
scalar function over both arguments is evaluated over:

```text
model × algorithm
```

### 5.3 Zip

Zip is not a syntactic rule on adjacent vectors.

Zip means:

> two expressions are aligned because they depend on the same root axis.

So if:

- `models : model`
- `resolution : model`

then they zip naturally because both are functions on the same root axis.

### 5.4 Inherited Axes

If an expression is derived from indexed inputs, it inherits the union of their root axes.

Examples:

- `simulate(model=models)` depends on `model`
- `diagnose(simulate(model=models))` also depends on `model`

### 5.5 Dependent Axes

Not every axis-like object is a new independent Cartesian dimension.

An axis may itself depend on the coordinate of another axis.

Example:

```text
model : root axis
resolution(model) : dependent axis
```

In this case:

- `model` is independent
- `resolution` is not an independent root axis
- the domain of `resolution` is evaluated locally for each `model`

This is the dynamic analogue of variable inner vector sizes or ragged nested indexing.

For example:

```text
domain(resolution | model = m1) = [r11, r12]
domain(resolution | model = m2) = [r21, r22, r23, r24]
```

The resulting index space is not a rectangular `model × resolution` tensor.

It is a ragged or fibered structure:

- `{m1} × resolution(m1)`
- `{m2} × resolution(m2)`

and so on.

### 5.6 Auxiliary Correspondences

It should be possible to map one axis domain into values of another kind while preserving the
source axis.

Example:

```text
model -> resolution
```

This does not create a new independent root axis.

Instead it creates a value indexed by `model`.

This kind of object is an auxiliary correspondence, not a new root axis.

Important distinction:

- if `model -> scalar resolution parameter`, this is an auxiliary correspondence
- if `model -> list of resolution points`, this is a dependent axis

---

## 6. First-Class Indexed Results

Once indexed expansion happens, the result is no longer just `T`.

It becomes:

```cpp
IndexedResults<T>
```

and this indexed result must be propagated downstream.

This implies:

- `IndexedResults<T>` is a first-class DSL type
- assignment may bind identifiers to `IndexedResults<T>`
- downstream functions may:
  - consume `IndexedResults<T>` directly
  - lift scalar `f(T)` elementwise
  - explicitly reduce/project indexed results

### Important Principle

Indexed expansion is not temporary sugar used only at the call site.

It is part of the value semantics of the DSL.

### Ragged Indexed Results

`IndexedResults<T>` must not be understood as inherently rectangular.

The indexed extension should allow:

- rectangular products of independent root axes
- ragged products involving dependent axes

Therefore `IndexedResults<T>` should be semantically compatible with sparse or row-wise storage.

Nested `IndexedResults<IndexedResults<T>>` is possible in principle, but is not the preferred base
model for the DSL.

The preferred model is:

- one indexed result
- one coordinate per realized payload
- support for dependent-axis coordinates whose domains are validated relative to parent coordinates

---

## 7. Overload Policy

The indexed extension requires explicit overload-preference rules.

### Preferred Ranking

1. Exact direct match
   Examples:
   - scalar actual with formal `T`
   - vector actual with formal `vector<T>`
   - indexed actual with formal `IndexedResults<T>`

2. Exact indexed/container-aware match
   Example:
   - `write_csv(IndexedResults<T>)`

3. Implicit lifted scalar match
   Example:
   - only `f(T)` exists, but actual is indexed

### Consequence

If both `f(T)` and `f(IndexedResults<T>)` exist and the actual argument is indexed, then
`f(IndexedResults<T>)` should win.

Otherwise defining the indexed overload would not be meaningful.

---

## 8. Proposed Dynamic Evaluation Strategy

The dynamic analogue of the `quimulun` style is:

1. collect root axes from the typed argument-expression tree
2. iterate coordinates over the resulting index space
3. evaluate the expression at each coordinate
4. fill an `IndexedResults<T>`

For rectangular cases this means iterating the Cartesian product of root axes.

For dependent-axis cases this means:

- first iterate parent coordinates
- then expand dependent-axis local domains relative to those parent coordinates
- then emit realized coordinates row-wise

### Proposed Expression-Level Interface

Each `typed_expression<T>` may support:

```cpp
collect_root_axes(env) -> index_space
eval_at(env, coordinate) -> Maybe_error<T>
```

with defaults:

- `collect_root_axes(env)` returns empty
- `eval_at(env, coordinate)` delegates to ordinary `run(env)`

Specializations/overrides are expected for:

- `typed_identifier`
- `typed_identifier_ref`
- `typed_conversion`
- `typed_vector_construction`
- `typed_tuple_construction`
- `typed_function_evaluation`

### Materialization Rule

```text
if collect_root_axes(expr) is empty:
    evaluate scalar T
else:
    iterate coordinate over index_space
    fill IndexedResults<T>
```

### Axis Collection Rule

`collect_root_axes` is intentionally named conservatively.

It should collect:

- root-axis dependencies directly
- enough metadata to discover dependent axes from derived expressions

So the resulting `index_space` may later need to represent:

- a set of independent root axes
- a set of dependent-axis definitions with parent links

The document does not require these to be represented by the same runtime type.

---

## 9. Axis Registry and Immutability

The environment should own an axis registry.

### Proposed Environment Responsibilities

- allocate fresh `axis_id`
- store `axis_id -> axis metadata`
- refine display labels after assignment if desired
- keep axis identities stable across downstream evaluation
- store parent dependencies for dependent axes
- support lookup of local dependent-axis domains when materializing ragged indexed results

### Why `axis_id` and `axis_label` Must Be Separate

Human-readable labels are not guaranteed unique:

- multiple arguments may want the label `model`
- temporary expressions may have no stable assigned name

Therefore:

- uniqueness must come from `axis_id`
- readability must come from `axis_label`

### Immutability Policy

Axes should be semantically immutable.

That means the following should not change once an axis exists:

- `axis_id`
- root-axis status
- semantic alignment meaning
- domain membership

Presentation metadata may be enriched later, but semantic axis identity must remain fixed.

---

## 10. Proposed Struct and Class Names

The following names are proposed for the indexed DSL extension.

### Runtime / Semantic Names

- `AxisId`
- `AxisLabel`
- `AxisInfo`
- `AxisDomain<T>`
- `LocalAxisDomain<T>`
- `AxisIndex`
- `Coordinate`
- `IndexSpace`
- `IndexedResults<T>`

### Optional Supporting Names

- `RootAxis`
- `DependentAxis`
- `AuxiliaryCoordinate`
- `AxisRegistry`
- `IndexedValue<T>` if a lighter-weight per-coordinate representation is later needed

### Expression Method Names

- `collect_root_axes`
- `eval_at`
- `materialize_indexed`

These names are preferred over older or more overloaded alternatives such as:

- `Position`
- `Index<...>` used ambiguously for both axis identity and integer slot
- `indexed_object`

---

## 11. Naming of New Root Axes

New root axes need:

- unique internal identity
- useful display labels

### Internal Identity

Fresh root axes receive a new `AxisId` from the environment axis registry.

This identity is unique and does not depend on any display name.

### Display Label Policy

Suggested display-label precedence:

1. if introduced from an assigned identifier, use that identifier name
2. otherwise use the function argument name if available
3. otherwise use the argument position
4. later assignment may refine presentation, but not identity

Example:

- provisional label: `model`
- later qualified label: `res.model`

Only the label changes.

The underlying `AxisId` does not.

---

## 12. Minimal Structural Sketch

The document does not fix the exact in-memory layout, but the intended shape is roughly:

```cpp
struct AxisId {
    std::size_t value;
};

struct AxisIndex {
    std::size_t value;
};

struct Coordinate {
    std::vector<AxisIndex> values;
};

template<class T>
struct AxisDomain {
    std::vector<T> values;
};

struct AxisInfo {
    AxisId id;
    std::string label;
};

struct IndexSpace {
    std::vector<AxisId> root_axes;
    std::vector<AxisId> dependent_axes;
};

template<class T>
struct IndexedResults {
    IndexSpace index_space;
    std::vector<Coordinate> coordinates;
    std::vector<T> values;
};
```

This is illustrative only.

The final implementation may use:

- a richer dependency graph
- row-wise sparse storage
- parent-linked dependent-axis domain objects

and does not need to be rectangular.

---

## 13. Open Design Questions

The following points remain open:

- whether indexed evaluation is always enabled or only in selected argument positions
- whether `IndexedResults<T>` should be sparse only or also have a dense variant
- how explicit auxiliary correspondences are represented in the DSL
- how dependent axes are represented at the expression and storage layers
- whether provenance/history should later be attached to every value
- how much source text should be preserved versus reconstructed from typed expressions
- how indexed overload ranking should be integrated into the current overload-resolution machinery

---

## 14. Summary

The indexed DSL extension is based on five main commitments:

1. `vector<T>` is interpreted contextually, either as data or as an axis source.
2. Root axes combine by Cartesian product.
3. Shared root axes imply zip/alignment.
4. Dependent axes allow ragged local domains derived from parent coordinates.
5. Indexed expansion is propagated downstream as `IndexedResults<T>`.
6. Axis identity is immutable and distinct from the user-facing label.

This provides a coherent path from current typed MacroIR evaluation to a dynamic indexed model
inspired by `quimulun` while remaining compatible with the current DSL structure.
