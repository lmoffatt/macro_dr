# figure_2 parallelization plan — combos serial, simulations + bootstrap parallel

## Goal
Move OpenMP parallelism from the **axis-combo** level (memory ∝ concurrent combos ×
nsim × evolution → OOM) to the **per-simulation** and **bootstrap-replicate** levels,
so memory is bounded by **one cell** (the dispatcher already runs one `(Num_ch, noise)`
cell per job) and the cores are filled by the within-cell simulation/bootstrap loops.

OpenMP does not nest by default: an inner `parallel for` inside the combo `parallel for`
serializes. So the combo loop must be **serial** for the per-sim/bootstrap loops to get
the threads. We gate it behind `MACRODR_AXIS_SERIAL=1` so other scripts are unaffected.

## Site 0 — gate the combo loop serial (prerequisite)
`include/macrodr/dsl/grammar_typed.h` lines ~802, ~1175, ~1457:
`#pragma omp parallel for schedule(dynamic, 1)` → add `if(!macrodr_axis_serial())`.
Add a cached env reader:
```cpp
inline bool macrodr_axis_serial() {
    static const bool v = [] { const char* e = std::getenv("MACRODR_AXIS_SERIAL");
                               return e && std::string(e) != "0"; }();
    return v;
}
```
Dispatcher sets `MACRODR_AXIS_SERIAL=1`. Default off → combo-parallel as today.

## Site 1 — per-simulation dlikelihood loop  (src/core/likelihood.cpp ~810)
`calculate_n_simulation_mdlikelihood_predictions_impl::run_loop`:
- fork the table: `auto forks = ftbl.fork(omp_get_max_threads());`
- pre-size `results(n)` + `errs(n)`; NO push_back.
- `#pragma omp parallel for schedule(static)` over i; thread uses `forks[omp_get_thread_num()]`;
  write `results[i]` / `errs[i]`.
- after loop: first non-empty `errs[i]` → return error, else `results`.

## Site 2 — numerical Fisher loop  (include/macrodr/cmd/likelihood.h ~219)
`calculate_n_simulation_mnumerical_fisher_information`: same shape (fork table, pre-size,
parallel for, per-index error). Inner 2·n_params loop stays serial.

## Site 3 — bootstrap  (legacy/bootstrap.h ~163 bootstrap_it_two_paired; ~125 bootstrap_it_two)
RNG must stay serial for determinism:
1. serial: pre-draw all replicate index vectors from `gen` into `std::vector<std::vector<size_t>> idxs(B)`.
2. parallel: `#pragma omp parallel for schedule(dynamic)` over b → `results[b] = f(vs1, idxs[b], vs2, idxs[b], args...)`.
   Use `f`/args as lvalues (called B times — do NOT `std::forward` in a loop).
3. move results into `bootstrap<R> out`.
`samples_dir` binary write stays serial (after the parallel pass).

## Thread-safety / determinism
- per-thread `FuncMap` fork → no Qdt-cache race.
- preset diagnostic fn must be pure (reads dy + indices, returns result) — verify no static/shared state.
- bootstrap RNG drawn serially → results seed-reproducible regardless of thread count.
- no throw across omp; collect `Maybe_error` per index, combine after.
- MKL is sequential (lp64_seq) → no hidden BLAS thread level.

## Validation (before trusting speedup)
1. numerics identity: same input, `OMP_NUM_THREADS=1` vs `=N` (combos serial both) → diagnostics CSV bit-identical.
2. memory: RSS plateaus at one cell's `nsim × evo`.
3. `macrodr_tests` + `contains_dib` green.

## Order
Site 0 + Site 1 first (validate), then Site 3, then Site 2.
