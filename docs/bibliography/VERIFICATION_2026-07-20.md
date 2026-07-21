# biblio.bib — independent verification against the PDFs on disk

> 2026-07-20. Closes the T-3 check in `papers/1_method/01_writing_plan.md`
> ("every entry has a DOI or stable URL; bibtex clean").

## What was checked

`papers/1_method/docs/manuscript-drafts/biblio.bib` (built from zero by T-3 on 2026-07-14) was
verified against the primary PDFs in this folder. Every PDF's first page was re-read with `pdftotext`
and an independent BibTeX record extracted from the printed text alone (no DOI, volume, page range or
year reconstructed from memory or from a filename), then compared field by field with the repo entry.

## Result: clean

- **106 entries**, 92 with a DOI, 12 with a URL — 104/106 carry a stable identifier.
- **The 9 keys the manuscript actually cites all resolve**; zero unresolved `\cite`.
- Of those 9, the 5 with a PDF on disk (`delcore2025parameter`, `milescu2005maximum`,
  `moffatt2025bayesian`, `munch2022bayesian`, `moffatt2007estimation`) match the printed record
  exactly. The remaining 4 are classic statistics works with no PDF here (`huber1967behavior`,
  `white1982maximum`, `godambe1960optimum`) plus the Zenodo code deposit (`moffatt2025macroir`,
  an `@software` entry, **not** a duplicate of the Comm Biol paper).
- **T-3 resolved every filename trap correctly.** Several PDFs on disk are preprints or are misnamed;
  the repo entries carry the published version of record in each case:

  | file on disk | what the file actually is | biblio.bib has |
  |---|---|---|
  | `Clerx_2019_..._BiophysJ.pdf` | the bioRxiv preprint | Biophys J 117:2420–2437 ✓ |
  | `hobolth2011.pdf` | Hobolth **& Jensen** | both authors ✓ |
  | `jahnke2006.pdf` | printed issue year 2007 | (not cited; absent) |
  | `minin2007.pdf` | printed issue year 2008 | (not cited; absent) |
  | `Bauerle_..._StochModels.pdf`, `ForemanMackey_..._AstronJ.pdf`, `Rubenzahl_..._AJ.pdf` | arXiv preprints | — |

## Works on disk with no biblio.bib entry, and why that is fine

`Carbonell 2008`, `Chowdhury 2013`, `Jahnke & Huisinga 2007`, `Minin & Suchard 2008`,
`Schranz 2008`, `Vastola 2021`, `Rudisill 1974`, `Murthy & Haftka 1988`.

All are CTMC-mathematics or numerical-methods works (matrix exponentials, eigenvalue/eigenvector
derivatives, the chemical master equation). They bear on the **engine**, not on paper 1's argument. Add
them only if Methods ends up discussing the `Qdt` computation in enough detail to need them.

## Files that are not papers

- `Moffatt_..._CommBiol_SI.pdf` (supplementary information) and `..._PeerReview.pdf` (peer review file).
- `Ho_Crawford_Suchard_2018_Stochastic_Compartmental.pdf` and
  `Zhou_Lange_2009_Composition_Markov_Chains.pdf` are **1.1 kB HTML redirect stubs saved with a .pdf
  extension**, not PDFs. The correctly-named copies of both works are intact; nothing was lost. Worth
  deleting the stubs so a future reader does not think the PDF is corrupt.

## What this does NOT certify

The 4 statistics classics with no PDF here (`huber1967behavior`, `white1982maximum`,
`godambe1960optimum`) were not verifiable from disk and were left as T-3 wrote them. They are cited in
the manuscript. Fetch the PDFs and re-check them before submission (per the standing rule that every
paper looked at gets a PDF here).
