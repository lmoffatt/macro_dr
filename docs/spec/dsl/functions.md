Currently implemented DSL functions in MacroIR, as registered via `push_function(...)` using `to_typed_function`.

---

# MacroIR DSL Functions

**Status:** Generated from `make_compiler()`
**Source:** `main.cpp`, `CLI_macro_dr_base.h`, `CLI_thermo_evidence_dts.h`, etc.
**Last updated:** 2025-07-17

---

## 📌 Simulation and Evidence

### `simulate`

```text
simulate(output, recording, experiment, init_seed, modelName, parameter_values, simulation_algorithm)
```

Simulates data using a model and records the result.
**Returns:** string (filename)

---

### `likelihood`

```text
likelihood(output, model, parameter_values, likelihood_algorithm, recording, experiment)
```

Computes likelihood under specified algorithm.
**Returns:** string (filename)

---

## 📌 Environment & Model Loading

### `load_Parameter`

```text
load_Parameter(filename, separator)
```

Loads parameters from CSV or TSV.
**Returns:** parameter value object

### `load_Prior`

```text
load_Prior(filename, separator)
```

Loads prior distribution from file.
**Returns:** prior object

### `get_Prior`

```text
get_Prior(prior_error, model)
```

Constructs a prior programmatically.
**Returns:** prior object

---

## 📌 Algorithm Configuration

### `simulation_algorithm`

```text
simulation_algorithm(include_N_states, number_of_substeps)
```

Defines settings for simulation.
**Returns:** simulation\_algo\_type

### `set_Likelihood_algorithm`

```text
set_Likelihood_algorithm(adaptive, recursive, averaging, correction, variance, n_sub_dt)
```

Configures the likelihood approximation strategy.
**Returns:** algorithm config

---

## 📌 Experiment Handling

### `get_Experiment`

```text
get_Experiment(filename, frequency_of_sampling, initial_ATP)
```

Parses an experiment from data.
**Returns:** experiment\_type

### `get_Observations`

```text
get_Observations(filename)
```

Loads raw experimental traces.
**Returns:** string (trace object)

### `idealize_Experiment`

```text
idealize_Experiment(experiment_filename, sep, idealized_filename)
```

Applies trace cleanup and binarization.

---

## 📌 Thermodynamic MCMC

### `set_ThermoAlgorithm_dts`

```text
set_ThermoAlgorithm_dts(num_scouts_per_ensemble, number_trials_until_give_up, max_iter_equilibrium, ...)
```

Creates a DTS-style MCMC sampler configuration.

### `thermo_evidence_dts`

```text
thermo_evidence_dts(idname, model, prior, likelihood_algorithm, data, experiment, thermo_algorithm, sampling_interval, max_samples, init_seed)
```

Runs thermodynamic MCMC for model evidence estimation.

### `thermo_evidence_dts_continuation`

```text
thermo_evidence_dts_continuation(idname, continuation_number, init_seed)
```

Continues a previous evidence run.

---

## 📌 Utility Functions

### `get_random_Id`

```text
get_random_Id(prefix)
```

Generates a random ID string.

### `get_number`

```text
get_number(n)
```

Echo function (for debugging or scripting).

### `write_text`

```text
write_text(filename, text)
```

Writes string content to a file.

### `write_script`

```text
write_script(script_name)
```

Persists the current program to disk.

---

Would you like this turned into an auto-updating reference from the compiler registry?

