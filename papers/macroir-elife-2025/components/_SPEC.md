# Component-file scheme — SPEC

> Canonical definition of how paper knowledge is stored as atomic component files. Governs everything
> under `components/`. **DRY: this is the only place the scheme is defined. Cite it, do not restate it.**
> v0.1, 2026-07-15. A novel *application* of established methods, not a new method. Whether it works: `_LOG.md`.

## Governing principle
**DRY** (Hunt & Thomas, *The Pragmatic Programmer*): *"every piece of knowledge must have a single,
authoritative, unambiguous representation within a system."* Each fact/conclusion has ONE owner file.
Everywhere else: **link, never copy.** This is the anti-regression guarantee.

## What a component is
Atomic unit of paper knowledge — one concern, one owner file (**Single Responsibility**, Martin).
Filename = ID. Path `components/<group>/<ID>.md`.

## Topic types — TWO FAMILIES (dimensional model; Kimball star schema)
**Facts reference dimensions; dimensions hold NO facts.** (see `_LOG` failure-mode #2, 2026-07-15)

**DIMENSIONS** — what you slice by. Hold intrinsic definition + typed POINTERS only. Nothing cross-cutting inside.
- **`entity`** — a thing with an identity (algorithm MR, a parameter).
- **`concept`** — an axis / term / device (the N_ch axis, boundary-state).

**FACTS** — assertions. Reference N dimensions; each referenced dimension lists them.
**Rule: a fact touching ≥2 dimensions MUST be its own node** (never inside a dimension).
- **`claim`** — a finding / conclusion / novelty. Schema = Toulmin + nanopub provenance.
- **`source`** — an artifact: figure, dataset, run, **or code location**.

**OPERATIONS** — code that transforms facts (a function, not a place). W3C PROV *activity*: **inputs → output**.
- **`operation`** — defined ONCE (the sandwich, a decomposition). Applications point to it; its applications
  list is DERIVED, never stored (see Rule 5). Precursor: PROV activity / dataflow DAG / Makefile rule.

## Schemas (Diátaxis *reference* density — no filler, fields/tables, scannable)
- **dimension** (`entity`/`concept`): **Definition** (intrinsic coordinates, `[LOCKED]`) + typed **pointer-lists** —
  `facts` · `claims` · `figures` · `data` · `code` · `neighbors` (dimensions) · `decisions` · `open`.
  A pointer to a not-yet-built node is a worklist item (Zettelkasten), not an error.
- **`claim`**: Claim · Grounds (evidence → `source` nodes) · Warrant · Qualifier (scope/strength) ·
  Rebuttal (attack surface) · Status {live|dead|conceded} · Provenance (decided where/when) · Dimensions-referenced · Links
- **`source`**: What (figure/data/run/code) · Where (path/id) · Shows (which fact) · Provenance · Dimensions-referenced

## Field provenance (tag every non-trivial field)
- `[LOCKED <src> <date>]` = a settled conclusion.
- `[OPEN → <where>]` = unresolved; pointer to the owner that will settle it.
- (= nanopublication assertion+provenance, applied.)

## Rules (the anti-regression, which is the whole point)
1. **One owner.** A contested or duplicated fact becomes its OWN node; others link to it. (DRY)
2. **No silent overwrite.** Never replace a `[LOCKED]` line without (a) a dated, sourced counter-statement
   AND (b) a `_LOG.md` entry. A regression must be a **visible downgrade**, never a quiet swap.
3. **Cite, don't copy** (DITA conref). Every number carries its source (this is LINT-SRC from `01_writing_plan.md`).
4. **Hierarchy lives in `_MAP.md`**, not in the nodes.
5. **Edges are ONE direction: dependent → dependency** (the application points to the definition/operation).
   The reverse (who applies/references X = backlinks) is **DERIVED by grep, never stored.** This is the whole
   anti-explosion move: you maintain N edges, not 2N and not N². The tool computes the "viceversa", not you.
6. **Read one hop.** Reading a node, you do not expand its pointers transitively; deep navigation is on demand.
   (This is Luciano's "recursividad limitada" — it is a rule, not a hope.)
7. **Operations define once.** An `operation` (the sandwich) is one node; consumers carry a `computed-via:
   [[op]]` edge; the operation's applications list is a derived backlink (Rule 5), never hand-maintained.

## Precursors (borrowed — look these up when the scheme fails)
DRY / SSOT (Hunt & Thomas) · Single Responsibility (Martin, SOLID) · DITA typed topics + conref (OASIS) ·
Diátaxis reference mode (Procida) · Zettelkasten atomic notes + links (Luhmann) · Frames/slots (Minsky 1974) ·
Nanopublications (Groth 2010) / Micropublications (Clark) · Toulmin argument model (1958).

## What is novel here (and therefore may fail)
The **assembly**: an agent-maintained single-source component graph for one *living* paper, with
anti-regression against LLM drift. Pieces are off-the-shelf; this composition is not a named system.
Failure modes → `_LOG.md`.
