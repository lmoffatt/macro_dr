# Error breadcrumbs — the pattern

How to attach context to `Maybe_error` failures so the user reads the right diagnostic.

This is a progressive convention. Existing code may not follow it; new code and code being touched for other reasons should adopt it.

## The principle

> **An error message is a navigation aid, not a snapshot.** It tells you *which run* and *which step* failed, so you can re-run that one case under a debugger or with extra logging. It is not the place to dump argument values.

A breadcrumb is one tiny identifying value (an index, a template choice, a coordinate label). If a breadcrumb requires more than ~30 characters to be useful, it is the wrong breadcrumb.

## What each breadcrumb should carry

**Yes** — short, indexing, identifying:

| kind | example | why |
|---|---|---|
| loop index | `k=234`, `i_step=12`, `i_state=3` | identifies *which iteration* failed |
| dispatch coordinate | `{algorithm=micro_IR, Num_ch=50}` | identifies *which lifted run* |
| template-parameter choice | `recursive=true`, `averaging=2`, `variance=true` | identifies *which compile-time variant* of the function fired |
| call-site name | `MacroR2`, `calc_micro_Qdt`, `Bayes:joint`, `to_Probability:projection` | identifies *which call inside this function* failed |
| parameter index | `param=on`, `param_idx=3` | identifies *which θ_j* the gradient/sum invariant was violated for |

**No** — bulk values, anything you'd need a viewer for:

| anti-example | why not |
|---|---|
| a `Matrix<double>` | size dominates the error string; impossible to read |
| a struct with many fields | same |
| the full `t_prior` / `t_Qdt` / `Pi` | these are exactly the things you'd need a debugger to inspect; an error string is the wrong medium |
| `repr(args)` of all arguments | always too long, never the right cutoff |

If a breadcrumb wants a bulk value, replace it with the *index* into the bulk. The index plus a re-run is enough to reach the data.

## The idiom

At every site where a `Maybe_error` is forwarded, do this inline — no helper wrapper, no macro:

```cpp
auto m = inner_call();
if (!m)
    return error_message("<my breadcrumb> | " + m.error()());
```

The breadcrumb prefix is whatever this layer adds that the *outer* layer doesn't already know. Separator is ` | ` (space-pipe-space). Read outer-to-inner, left-to-right.

### What "smart" means

The discipline is to *not* add a breadcrumb when the layer above already provides the same information.

- A function whose only context is "I exist" (no loop, no template choice, no dispatch) doesn't earn a breadcrumb. Forwarding `m.error()` raw is correct.
- A function called from one site only doesn't need to add its own name (the caller's breadcrumb already implies it).
- A loop body adds the iteration index, *only the iteration index*. Not the function name (visible from the source), not the bulk arguments.
- A dispatch over template parameters adds *which branch* fired (e.g., `MacroR2 (Binomial)` vs. `MacroR2 (non-Binomial fallback)` — both call the same primitive but with different template arguments).

If you can't decide whether to add a breadcrumb, ask: *would the chain be unambiguous without this one?* If yes, skip it.

## The recommended layering

For a likelihood-style call chain, breadcrumbs accrete in this order — outermost first:

```
{algorithm=micro_IR, Num_ch=50} | k=234 | Micror[recursive=true, averaging=2] | Bayes_Rule: zero evidence (...)
```

| layer | breadcrumb | source |
|---|---|---|
| lifted-function dispatch | `{algorithm=micro_IR, Num_ch=50}` | `Coordinate::str()` from the index space |
| per-trace (when applicable) | `trace=2` | the trace iterator |
| per-step | `k=234` | the fold loop body |
| per-template-variant call site | `Micror[recursive=true, averaging=2]` or `MacroR2 (Binomial)` | the dispatch within the step |
| per-named-call-site (only when the same step has *multiple* failure-prone calls) | `Bayes:joint`, `to_Probability:projection`, `to_Cov:M+diag` | within the function, if it has more than one Maybe-returning call |
| innermost (the actual numerical failure) | `to_Probability: \|Σ − 1\| = 0.220879` | the primitive itself |

## Templated values *are* indices

Compile-time choices like `recursive`, `averaging`, `variance`, `variance_correction` take a small finite set of values. They are *exactly* indices: each picks one branch out of a few, and "which branch" is always part of "which run failed."

When wrapping a Maybe-returning call that's templated on these, the breadcrumb should encode them — readable form preferred:

```cpp
"MacroR2[recursive=" + std::string(recursive::value ? "true" : "false") +
 ", averaging=" + std::to_string(averaging::value) + "]"
```

…or, since these are usually known by name in this codebase, lean on existing `ToString` helpers (e.g. `ToString(MacroR2<recursive, averaging, variance, variance_correction>{})`) when one exists.

## What we *don't* do (in this codebase, by choice)

1. **No `with_context(maybe, ctx)` helper.** The user prefers each function to author its own breadcrumbs inline. Helpers invite over-wrapping; inline forces the author to think about whether the breadcrumb earns its place.

2. **No automatic argument printing.** No "compiler-style" `function_name(arg1, arg2, …)` dumps. Bulk args don't fit the medium.

3. **No `std::source_location`-based file:line auto-capture (yet).** Could be added later as a free win — `error_message` could take a `source_location` defaulted to `current()` and prepend `file:line:function`. This would let breadcrumbs drop the function-name component and carry only domain values. Not done today; a future cleanup.

4. **No structured-logger split.** Errors are flat strings. A future refactor could move to a `{coord, k, fn, msg}` structure for machine-parseable diagnostics, but that's a bigger lift than the current incremental approach.

## How to migrate progressively

Don't refactor all sites at once. Touch them when:

- you're already editing the function for another reason, or
- you hit a bug whose error message was unhelpful — that's the cue that *this* function's breadcrumb is missing.

The payoff is per-call-site: each breadcrumb earns its place when an error message becomes one click more navigable.

## Worked example: figure-2 micro path

Before:

```
to_Covariance_Probability: |Σ − 1| = 0.220879
```

After current breadcrumbs (lifted dispatch + fold k):

```
{algorithm=micro_IR, Num_ch=50} | k=234 | Micror | to_Covariance_Probability: |Σ − 1| = 0.220879
```

After full pattern (template-variant + named call site):

```
{algorithm=micro_IR, Num_ch=50} | k=234 | Micror[recursive=true, averaging=2] | macro-shape projection | to_Covariance_Probability: |Σ − 1| = 0.220879
```

Each segment narrows the search by an order of magnitude. The bulk-data path (`t_prior`, `t_Qdt`, the actual matrices) is *never* in the error string — they're recovered by re-running just `{algorithm=micro_IR, Num_ch=50}` at step `k=234` under a debugger.
