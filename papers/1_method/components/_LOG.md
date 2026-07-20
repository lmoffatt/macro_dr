# Component scheme — does-it-work LOG

> Append-only, dated. Failure modes, observations, adjustments to the scheme (`_SPEC.md`).
> The scheme is a novel application (`_SPEC` "what is novel"); it will fail in ways we did not anticipate.
> This file is where we catch that. An entry = {date · observation · failure-mode? · rule change?}.

## 2026-07-15 — v0.1 opened
- Scheme defined (`_SPEC.md`), precursor-backed. First component built: `algorithms/MR.md` (type=reference).
- The `[LOCKED src date]` provenance tag emerged organically on MR **before** the spec existed → validates
  the nanopublication-provenance choice (we reached for it without being told; good sign it fits).
- **OBSERVATION / first failure mode.** MR (a `reference` entity) carried a **contested attribute** — the
  "MR sign" (observable variance vs parameter covariance). It did not fit the flat reference schema; it needed
  an extra "Tension" field.
  - **Failure mode:** where does a contested/unresolved attribute live — inside the entity, or as its own node?
  - **Proposed rule (v0.1):** a contested attribute becomes a `claim` node (Toulmin schema); the entity states
    only the settled part and links to the claim. NOT yet retrofitted: the owning claim (`MR-sign`) is D-4's,
    still pending Luciano approval, so MR carries the tension inline tagged `[OPEN → D-4]` until then.
  - **Do:** when D-4 is approved, split `MR-sign` into `claims/MR-sign.md` and thin MR's Tension to a link.

## 2026-07-15 — FAILURE MODE #2 (Luciano caught it): dimension vs fact not explicit
- MR (a dimension) was hoarding facts that reference OTHER dimensions ("MR worse than R" is about R + the
  window axis too). DRY implied "own node + link", but the type system did not separate **dimension** from
  **fact**, so the violation was easy to make (I made it). This is the first real design failure.
- **Precursor:** dimensional modeling / star schema (Kimball) — facts reference dimensions, dimensions hold
  no facts. Also RDF resource-vs-statement, faceted classification, nanopub entity-vs-assertion.
- **Fix (pending Luciano sign-off before applying):** two families ×2 types.
  - dimension = `entity` (algorithm, param) + `concept` (axis, term) → intrinsic coordinates + POINTERS only.
  - fact = `claim` (finding/conclusion) + `source` (figure, data, run) → references N dimensions; a fact
    referencing ≥2 dimensions MUST be its own node.
- **Retrofit MR** to coordinates {recursive, av=1} + pointers; its mechanism / over-confidence /
  non-monotonicity / sign become fact nodes.

## 2026-07-15 — REFINEMENT #3 (Luciano): functions + the bidirectional-pointer kilombo
- Missing: **functions** (code that transforms facts, not just locates it — the sandwich). Added type
  `operation` = **W3C PROV activity** (inputs → output). Precursor named; not reinvented.
- Luciano's real barrier stated: bidirectional pointers (define ↔ apply) make his "sense of complexity rise
  exponentially", so he never attempts systematization. **Dissolved by Rule 5:** edges are ONE direction
  (consumer → definition); the reverse is DERIVED by grep, never stored. Maintain N edges, not 2N/N².
- **DEMONSTRATED (works):** `diagnostics/sandwich.md` stores 0 applications; `facts/fact-MR-overconfident.md`
  carries the only edge (`computed-via: [[sandwich]]`); `grep -rl '\[\[sandwich\]\]'` derives the backlink.
  No double bookkeeping appeared. The kilombo did not materialize.
- Rules 5 (one-direction+derived backlinks), 6 (read one hop = his "recursividad limitada"), 7 (operations
  define once) added to `_SPEC`.

## 2026-07-15 — VERDICT: scheme PARKED; master_list is the spine
- Analysis of "will this work" found the killer: this scheme was a 3rd/4th copy of knowledge already owned
  by `nomenclature.md` / `decisions/` / section plans (DRY violated at birth). And the minimalist index
  Luciano wanted **already exists** = `00_master_list.md` (one-owner rule + status registry + topic index).
- **Decision:** adopt master_list as the live index (YAGNI); PARK this component scheme. Revive ONLY when a
  `_LOG` entry records master_list failing — the concrete test is the **D-0 rerun-refresh** (are the changed
  numbers a grep-and-replace or a re-derivation?).
- Actions taken: master_list got a **completeness guarantee** ("not indexed = not owned/not decided") +
  `check.sh` item 9 enforcing it mechanically (green: every working .md registered). The three demo nodes
  moved to `../archive/components_poc/`. Registered the scheme as PARKED in master_list.
- **Vocabulary kept** (Luciano's, better than the 5 types): **fact = data** (figure/CSV/code = `source`) ·
  **interpretation = reading** (single-dimension → an attribute of that dimension; multi-dimension → its own
  node). Best unit for writing = **interpretation-of-a-figure across algorithms/variables** (= a Results paragraph).
- What survived as proven value (not wasted): backlink-by-grep dissolved the double-pointer fear; MR caught a
  real regression (the sign); the `[LOCKED]` provenance tag; the fact/interpretation distinction.

### Watchlist (anticipated failure modes — confirm or refute as we go)
- A fact relevant to TWO components: does DRY "one owner + link" actually hold, or do we duplicate? (test)
- Does a `[LOCKED]` regression actually get caught by git-diff + Rule 2, or does it slip? (the real test of the terror)
- Do agents restate instead of link? (the copy-paste relapse)
- Does `_MAP.md` drift from the actual node set? (index rot)
- Type ambiguity: some things are both entity and claim (MR). Is 3 types enough, too many, wrong cut?
