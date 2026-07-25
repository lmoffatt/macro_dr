# Flat model declaration language (level 0)

Status: draft for discussion. Scope: `macroir` v1. Allosteric/conformational
generation (level 1) is explicitly out of scope; it goes with its own paper.

## Why

Today a kinetic scheme is declared three times in three incompatible forms
(`legacy/models_simple.h`, `legacy/models_MoffattHume_linear.h`):

1. a symbolic `Q0_formula` / `Qa_formula` / `g_formula` of strings,
2. a hand-indexed lambda `p -> Patch_Model` that is what actually runs,
3. a hand-indexed inverse lambda `Patch_Model -> p`.

Only (2) feeds the numerics, so (1) and (3) drift silently. Two drifts already
present in the tree:

- `scheme_CO` names its parameters `on` / `off` but its formula says `kon`.
- `scheme_1` sets `v_g_formula()[4] = "unitary_current"` while the running
  lambda uses `p[4] * -1.0`. The sign of the current exists only in the lambda.

Adding a scheme also means editing three hand-maintained `Models_Library(...)`
lists and recompiling (~25 min on the critical TU), because
`Model_v = decltype(get_model(...))` threads a variant over every registered
model through every instantiation.

Goal: one representation, loaded at runtime, from which Q0, Qa, g, the parameter
names, the inverse map and the printout are all derived.

## Canonical form

A model is five blocks. All sizes are runtime; nothing is a template parameter.

### 1. states

Ordered list of unique names. Length `k`. Index order defines the matrix order
and is part of the model's identity (fixtures compare index by index).

```
states: C0, C1, C2, C3, O
```

### 2. parameters

| column      | type   | notes                                   |
|-------------|--------|-----------------------------------------|
| `name`      | string | unique                                  |
| `transform` | enum   | `log10` \| `linear`                     |
| `value`     | double | in natural units, not transformed       |

Kinetic parameters are `log10` in every scheme in the tree. `linear` exists and
is used (`scheme_1` sets it for the baseline), so it must be per-parameter, not
a global switch.

### 3. rates

One directed edge per row. `i -> i` is not allowed; the diagonal is computed,
never declared.

| column   | type    | default | notes                                         |
|----------|---------|---------|-----------------------------------------------|
| `from`   | state   |         |                                               |
| `to`     | state   |         |                                               |
| `coef`   | double  | `1`     | combinatorial multiplicity (the `3` in `3*kon`) |
| `agonist`| string  | none    | name of the ligand, or none                   |
| `power`  | integer | `1`     | exponent on the agonist concentration         |

and, in a companion table, the parameter factors (one or more per edge):

| column     | type    | default | notes                        |
|------------|---------|---------|------------------------------|
| `from`     | state   |         |                              |
| `to`       | state   |         |                              |
| `param`    | string  |         | must exist in `parameters`   |
| `exponent` | integer | `1`     |                              |

so that

```
q_ij = coef * prod_p (theta_p ^ exponent_p) * [A] ^ power
```

The common case is exactly one parameter factor with exponent 1, which is every
edge in every scheme currently in the tree. Products with exponents cost nothing
extra to support and are what constrained parameterizations (detailed balance,
allosteric factors) need later.

The split `Q0` / `Qa` in the current code is exactly `power == 0` versus
`power == 1`. It is not a separate concept and does not survive into this form.

### 4. conductance

| column   | type   | default | notes                                     |
|----------|--------|---------|-------------------------------------------|
| `state`  | state  |         |                                           |
| `coef`   | double | `1`     | **carries the sign**                      |
| `param`  | string |         | must exist in `parameters`                |

States absent from the table have conductance 0.

The sign convention is explicit here and nowhere else. Inward current is
`coef = -1`, which reproduces the current `p[i] * -1.0` while keeping the
parameter itself positive and therefore `log10`-transformable.

### 5. initial distribution

One of:

- `state: <name>` — all probability in one state (what every scheme in the tree
  does today, written as `kon * (1.0/kon)` to carry the derivative type),
- `equilibrium` — `peq` of `Q(0)`, i.e. at zero agonist,
- `explicit: [p_0, ..., p_{k-1}]`.

### 6. observation model

Separate block; not kinetics.

| field                  | required | current name          |
|------------------------|----------|-----------------------|
| `n_channels`           | yes      | `N_Ch_mean`           |
| `current_noise`        | yes      | `Current_Noise`       |
| `baseline`             | yes      | `Current_Baseline`    |
| `pink_noise`           | no       | `Pink_Noise`          |
| `proportional_noise`   | no       | `Proportional_Noise`  |

## Not part of the model

These are currently written inside every model definition and must move to
solver options, where they can be set once and reported:

`min_P(1e-12)`, `Probability_error_tolerance(1e-2)`,
`Conductance_variance_error_tolerance(1e-2)`, `Binomial_magical_number(5.0)`,
`N_Ch_mean_time_segment_duration(121)`.

Also: declaring a model must not have side effects. Today declaring an
allosteric scheme opens `std::ofstream f("scheme_N_model_description.txt")` in
static initialisation and writes to the current directory. Printing is a
function you call.

## Surface syntax

The canonical form is tables, which is what R data frames, Python dicts and
JSON all express directly. For humans and for papers, a one-line-per-edge text
form parses to it. The grammar is small enough to fit here:

```
rate    := coef? factor ( '*' factor )*
coef    := number
factor  := param [ '^' integer ] | '[' agonist ']' [ '^' integer ]
```

`scheme_1` (p2x2, 5 states) in full:

```
states: C0, C1, C2, C3, O

rates:
  C0 -> C1 : 3 * kon * [ATP]
  C1 -> C0 : koff
  C1 -> C2 : 2 * kon * [ATP]
  C2 -> C1 : 2 * koff
  C2 -> C3 : kon * [ATP]
  C3 -> C2 : 3 * koff
  C3 -> O  : gating_on
  O  -> C3 : gating_off

conductance:
  O : -unitary_current

initial: state C0
```

Compare with the current source, which needs eight `v_Q0_formula()[i][j] = "..."`
assignments plus a 30-line lambda that repeats the same content with literal
indices, plus a second lambda that inverts it.

## Why this makes the derivatives easy

Every rate is a monomial in the parameters, and kinetic parameters are `log10`.
So, with `theta` the transformed parameter vector,

```
log q_ij = log(coef) + ln(10) * sum_p exponent_p * theta_p + power * log[A]
```

is **affine in theta**. Therefore `dQ/dtheta_p` is a constant sparse matrix,
computable once when the scheme is loaded, with entries `ln(10) * exponent_p * q_ij`.

Consequence: no automatic differentiation is needed to get from parameters to Q,
and no expression AST is needed. The forward sweep is seeded exactly and the AD
burden starts only at the matrix exponential. This removes a large part of what
was identified as the main technical risk of the rewrite.

Detailed balance, when it is added, is also linear in `theta` (the product of
rates around a cycle equals the product of the reverse rates), so a constrained
parameterisation is a null-space projection of an integer matrix. Still linear.

## Validation the loader must enforce

Each rule corresponds to a class of error that the current representation
cannot detect.

- every `param` referenced exists in `parameters`; every declared parameter is
  referenced (warn, do not fail: observation parameters are declared separately);
- state names unique, parameter names unique;
- no `i -> i` edge; no duplicate `(from, to)` in the transitions table;
- the transition graph is connected as an undirected graph, and every state is
  reachable, otherwise `peq` is not defined and the model is not identifiable;
- at least one state with non-zero conductance and at least one with zero,
  otherwise the recording carries no information about the kinetics;
- if `initial: equilibrium`, `Q(0)` must be irreducible;
- rows of the assembled `Q` sum to zero to machine precision (a check on the
  assembler, not on the declaration).

## Exact validation against the current code

The model layer is integers and strings, so it validates with no tolerance at
all. For each scheme currently registered (`scheme_CO`, `scheme_CCO`,
`scheme_COC`, `scheme_1`):

- same state count,
- same sparsity pattern of `Q0` and `Qa`,
- same parameter names in the same order,
- same `coef` per edge,
- same conductance vector including sign,
- and `Q` assembled at the reference parameter values equal to the binary's,
  element by element.

This is the one part of the rewrite that can be proven correct rather than
compared to a tolerance, and it should be done before any numerics are written.

## Open decisions

1. Multiple agonists. The record carries an agonist name, so the form supports
   them, but the experiment/protocol side currently assumes one. Decide whether
   v1 rejects more than one or carries it through.
2. Time-varying `n_channels`. The current code reads `N_Ch_mean` as a range
   (`p[pair(i, i)]`) with `N_Ch_mean_time_segment_duration(121)`, so the
   machinery for a per-segment channel count exists. v1 proposal: scalar only,
   with the field typed as a vector so it generalises without a format change.
3. Whether `coef` is a double or a rational. Every value in the tree is a small
   integer. A double is simpler; a rational keeps the log-linearity exact.
4. Parameter transforms beyond `log10` and `linear` (e.g. logit for a
   probability-valued parameter). Not needed by anything in the tree.
