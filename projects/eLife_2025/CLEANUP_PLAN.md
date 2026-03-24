# eLife 2025 Root Cleanup Plan

This document classifies the files currently sitting at the root of
`projects/eLife_2025/` and records where they should live.

## Current status

The main cleanup described here has now been applied:

- root-level figure datasets were moved into `figures/data/`
- root-level model and prior inputs were moved into `data/models/`
- debug and scratch residue were moved into `archive/debug/`

The project root is now limited to identity and coordination files.

## Guiding rule

Root-level files in this project should be limited to:

- project identity documents
- a small number of project-level coordination documents

Operational inputs, generated datasets, and debug residue should live in
owning subdirectories.

## Current project-level root files

These belong at the project root:

- `README.md`
- `PROJECT.md`
- `PROJECT.yaml`
- `CLEANUP_PLAN.md`

## Model inputs

These are stable model or prior inputs and should live under:

- `data/models/`

Files:

- `scheme_CO_par.csv`
- `scheme_CO_prior.csv`
- `scheme_CCO_par.csv`
- `scheme_CCO_prior.csv`
- `scheme_COC_par.csv`
- `scheme_COC_prior.csv`
- `scheme_1_par.csv`
- `scheme_1_prior.csv`
- `scheme_10_inact_par.csv`
- `scheme_10_inact_prior.csv`
- `scheme_5_model_description.txt`
- `scheme_6_model_description.txt`
- `scheme_7_model_description.txt`
- `scheme_8_model_description.txt`
- `scheme_9_model_description.txt`
- `scheme_10_model_description.txt`
- `scheme_11_model_description.txt`
- `scheme_12_model_description.txt`
- `scheme_13_model_description.txt`
- `scheme_14_model_description.txt`
- `scheme_15_model_description.txt`

Notes:

- `scheme_CO_*`, `scheme_CCO_*`, and `scheme_COC_*` already exist under
  `data/models/`; the root-level copies should be treated as duplicates.
- `scheme_1_*`, `scheme_10_inact_*`, and the model-description `.txt` files
  should also be moved into `data/models/` for consistency.

## Figure datasets

These are figure-facing generated datasets and should live under:

- `figures/data/`

Files:

- `figure_1_dlikelihood_predictions.csv`
- `figure_1_likelihood_diagnostic.csv`
- `figure_1_likelihood_diagnostic_IR.csv`
- `figure_1_likelihood_diagnostic_IRV.csv`
- `figure_1_likelihood_diagnostic_MNR.csv`
- `figure_1_likelihood_diagnostic_MNRV.csv`
- `figure_1_likelihood_diagnostic_MR.csv`
- `figure_1_likelihood_diagnostic_MRV.csv`
- `figure_1_likelihood_diagnostic_NR.csv`
- `figure_1_likelihood_diagnostic_R.csv`
- `figure_1_likelihood_predictions.csv`
- `figure_1_simulation.csv`
- `figure_2_dlikelihood_analysis_MRV_test.csv`
- `figure_2_dlikelihood_predictions.csv`
- `figure_2_dlikelihood_predictions_IRV.csv`
- `figure_2_dlikelihood_predictions_IRV2.csv`
- `figure_2_dlikelihood_predictions_MRV.csv`
- `figure_2_dlikelihood_predictions_MRV2.csv`
- `figure_2_dlikelihood_predictions_MRV_test.csv`
- `figure_2_simulation.csv`
- `figure_2_simulation_test_test.csv`

Notes:

- many of these files already exist under `figures/data/`; the root-level copies
  should be treated as duplicates or misplaced outputs
- test-oriented figure datasets can either remain in `figures/data/` or later be
  split into a figure-local `test/` subfolder if they accumulate

## Anomalous/generated leftovers

These should not remain at the project root.

### Scratch or debug artifacts

- `temp_script_file.txt`
- `stderr.log`
- `stderr2.log`
- `valgrind_figure1.log`

Recommended destination:

- `archive/` if kept for history
- or a future `tmp/` / `debug/` area if they are still actively useful

### Ambiguous generated artifact

- `figure_1_simulation`

This appears to be a generated output without a normal extension and should be
either:

- moved into `figures/data/` if still needed
- or archived/removed if superseded by `figure_1_simulation.csv`

## Recommended next cleanup move

1. treat `figures/data/` as the canonical location for figure datasets
2. treat `data/models/` as the canonical location for model/prior/description inputs
3. move or archive debug residue out of the root
4. leave only project identity and coordination files at the root
