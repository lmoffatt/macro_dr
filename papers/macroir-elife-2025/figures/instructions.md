# MacroIR / eLife 2025 Figure Instructions

This note extracts the figure-related guidance we need from the latest eLife author instructions checked on 2026-03-25.

Primary sources:

- https://submit.elifesciences.org/html/elife_author_instructions.html
- https://elifesciences.org/about/submit-your-research

This file is project-facing. It translates the journal rules into a concrete prep checklist for the figures produced under `projects/eLife_2025/figures/`.

## 1. Hard eLife requirements

### Initial submission PDF

- Embed figures in the manuscript PDF near where they are discussed.
- Keep each figure on a single page with its title and concise legend on the same page.
- Number figures in citation order.
- Keep composite figures compact enough that the figure and legend do not spill over multiple pages.
- Make all labels and legend text legible at normal reading size.

### Separate figure uploads

- Upload each main figure as a separate file.
- Upload each figure supplement as a separate file.
- Accepted formats: `TIFF`, `JPG`, `EPS`, `AI`, `PDF`.
- Use vector output when possible.
- If exporting `TIFF`, use 8-bit TIFF.
- Do not upload multi-page figure files.
- Minimum resolution: `300 dpi`.
- Minimum width: `10 cm`.
- If the figure is intended to span a full page, upload at `20 cm` width and `300 dpi`.
- Use `RGB` color for revised submissions.
- Remove unnecessary white margins.
- Turn off hidden layers before export.

### Legends and reporting details

- Each figure legend must state enough information to understand the panel without returning to Methods.
- Report the statistical test used where relevant.
- Report exact sample size `n`.
- Report exact `p` values where relevant.
- State replicate structure.
- State any exclusion or filtering criteria that affect the plotted data.
- If the reporting details become too long for the legend, move them into a linked source data file and mention that in the legend.

### Special content rules

- Micrographs and micrograph videos must include scale bars.
- Gels and blots must include molecular weight markers.
- Error bars must be defined explicitly in the legend.
- Source data for gels/blots must include raw unedited files and uncropped labeled versions.

### Figure supplements and source data

- Each supplement belongs to exactly one main figure.
- Name supplements as `Figure N—figure supplement M`.
- Numbering restarts from `1` for each main figure.
- Mention each supplement in the main figure legend and in the article text.
- Each source data file must correspond to exactly one figure or one table.
- Name source data as `Figure N—Source Data M`.
- Source data files should stay below `100 MB` each.

## 2. Practical standards for this paper

These are not extra journal rules. They are the standards we should follow so the figures survive review, revision, and production cleanly.

### Layout

- Prefer 1-column or 2-column compositions that can be read without zooming.
- Keep panel lettering simple and consistent: `A`, `B`, `C`, `D`.
- Align panels to a common grid.
- Avoid decorative elements that do not carry information.
- Keep axis titles, tick labels, and annotations consistent across figures.

### Typography

- Use one sans-serif family consistently across all figures.
- Do not let any label fall below a clearly readable print size.
- Keep notation identical to the manuscript text and Methods.
- Use sentence case for panel titles unless a symbol-only title is clearer.

### Color and accessibility

- Use color only when it adds information.
- Follow color-universal design: figures must remain interpretable under common color-vision deficiencies.
- Do not rely on color alone to distinguish conditions. Also vary line type, point shape, or ordering.
- Keep the same color mapping across figures for the same algorithm or condition.

### Statistical display

- Every uncertainty band or error bar must be defined in the legend.
- Distinguish clearly between variation across simulations, uncertainty in the mean, and posterior uncertainty.
- Avoid significance stars without exact values and a stated test.
- If smoothing or binning is applied, say so in the legend or source data.

### Export discipline

- Keep a vector master for each figure whenever possible.
- Rasterize only components that truly need it.
- Export a submission-ready file separately from any editable working file.
- Trim empty canvas space before final export.
- Check the exported PDF at 100% zoom before considering it finished.

## 3. Project-specific deliverables

For each main figure in this project, prepare the following bundle:

- `Figure N.pdf`
- Editable source used to generate it, for example the corresponding `Rmd`, script, or drawing source
- `Figure N legend.md` or equivalent manuscript-ready legend text
- `Figure N—Source Data 1.csv` if the figure is data-backed
- Optional supplement files as `Figure N—figure supplement M.pdf`

Suggested working mapping to the current repo:

- Figure generation code: `projects/eLife_2025/figures/`
- Storyboard/spec: `papers/macroir-elife-2025/04_figures_storyboard.md`
- Manuscript drafts: `papers/macroir-elife-2025/docs/manuscript-drafts/`

## 4. Figure-by-figure checklist

Run this checklist before calling a figure done.

- The figure answers one explicit question from the storyboard.
- Panel order matches the narrative in the manuscript.
- All symbols, abbreviations, and color encodings are defined.
- Axis limits and transformations are intentional and documented.
- Units are shown everywhere they are needed.
- Error bars, ribbons, and intervals are defined.
- The legend reports `n`, test, and exact `p` values where applicable.
- The final export is a single-page figure file.
- The exported file is at least `10 cm` wide and `300 dpi` equivalent if rasterized.
- White margins are trimmed.
- The source data file matches the plotted values.
- File naming follows eLife conventions exactly.

## 5. Current figure implications for MacroIR paper

Based on the current storyboard, the likely figure risks are:

- Figure 1: conceptual schematics can become too text-heavy; keep the punchline visual.
- Figure 2: multi-line comparisons can become unreadable; enforce consistent line styles and direct labeling where possible.
- Figure 3: diagnostics panels can drift into Methods; keep only the essential validation view in the main figure and push extras to supplements.
- Figure 4: validity maps must use a color scale that remains interpretable in grayscale and by color-blind readers.

## 6. Recommended workflow

1. Finalize the scientific question for the figure in `04_figures_storyboard.md`.
2. Generate data in `projects/eLife_2025/ops/local/` or the reproducible pipeline.
3. Produce the figure from `projects/eLife_2025/figures/`.
4. Write the legend at the same time as the plot, not afterward.
5. Export a clean single-file submission version.
6. Save or generate the matching source data file.
7. Review the figure at manuscript scale inside the draft PDF.

## 7. Open decisions for us

These are not journal blockers, but we should make them explicit as we prepare the submission:

- Choose a single project-wide font and panel-letter style.
- Choose a stable algorithm-to-color mapping and keep it fixed across all figures.
- Decide which diagnostics remain in main figures versus supplements.
- Decide whether any figure needs a full-page `20 cm` export.

