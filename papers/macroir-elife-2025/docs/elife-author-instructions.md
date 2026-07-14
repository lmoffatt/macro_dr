# eLife Instructions for Authors (local copy)

**Captured:** 2026-06-18. **Completed and corrected: 2026-07-14.** Faithful transcription
of eLife's author guide. Figure specifications come from eLife's figure-preparation guidance.

**Retrieval note.** `https://reviewer.elifesciences.org/author-guide/full` 301-redirects to
`https://elife-rp.msubmit.net/html/elife-rp_author_instructions.html`, which is the single
full guide (all tabs on one page). **WebFetch gets HTTP 403; plain `curl` with a browser
user-agent returns 200.** Use curl to re-check this document.

**Status: complete.** All sections now captured: Editorial Process, Article Types,
**Manuscript Preparation**, Reviewed Preprints, **Revised Submissions**, **Versions of
Record**, Data Availability, **Publishing Policies**, **Publication Fees**, **Media
Policy**, **FAQs**.

> **Two contradictions inside eLife's own pages (2026-07-14).**
> 1. **Abstract length.** The live guide says **150–200 words**; the LaTeX template
>    (`elife.cls` v1.11, dated 2022-06-01) still says "no more than 150 words". The widely
>    cited "150" comes from the template. **Write to 150, treat 200 as the ceiling.**
> 2. **Fee.** The Fees tab and the public site say **$3,750**; the FAQ tab on the same page
>    still says $3,000. The FAQ is stale.

**Sources**

- Author Guide — https://reviewer.elifesciences.org/author-guide/editorial-process
- Submit your research — https://elifesciences.org/about/submit-your-research
- Full author instructions (browser only) — https://elife-rp.msubmit.net/html/elife-rp_author_instructions.html
- LaTeX template — download from the author guide on elifesciences.org

> Project-facing figure subset already lives in
> `papers/macroir-elife-2025/figures/instructions.md` (checked 2026-03-25) and is
> consistent with this document.

---

## Editorial Process

eLife combines the immediacy and openness of preprints with peer review by experts.
Key features:

- eLife **only peer reviews submissions that are available as preprints**.
- Editors prioritise submissions where reviews will be of greatest public value —
  results significant within a field, or on a topic of broad interest.
- In general eLife will **not** review papers that report incremental results, are
  obviously flawed, or for which appropriate, unbiased reviewers cannot be recruited.
- **No accept/reject decision after peer review**: every reviewed article is
  published on the eLife website as a **Reviewed Preprint** — an integrated
  publication that includes the article, an eLife assessment, public reviews, and a
  response from the authors (if available).
- Where reviewers disagree with the manuscript's claims, this is made explicit in the
  eLife assessment and public reviews.
- Editors and reviewers discuss their reviews with each other and produce a written
  assessment of the **significance of the findings** and the **strength of the
  evidence**, using a common vocabulary for consistency.
- Two further outputs: (i) **Public Reviews** on strengths/weaknesses and whether
  claims are justified by the data; (ii) **Recommendations for the authors** (revision
  requests and suggestions).
- The eLife assessment and Public Reviews are posted alongside the preprint **two
  weeks** after being sent to the authors; authors may respond or ask for factual
  corrections. Recommendations for the authors are not made public at this stage.
- Authors may submit a **revised preprint at any time**; editors decide if it warrants
  a new assessment and/or new Public Reviews. The updated Reviewed Preprint also
  includes the recommendations and the authors' response to them.
- Authors can request that eLife produces a **Version of Record (VOR)**.
- **Fee:** **$3,750** for submissions from **July 1, 2026** (corrected 2026-07-14; the old
  $3,000 figure is stale and still appears in eLife's own FAQ tab), charged at the point
  eLife commits to peer reviewing the work. Full waivers available, see Publication Fees below.

Authors may submit directly, or transfer files from **bioRxiv** or **medRxiv** if a
preprint is already posted. eLife also welcomes submissions from **Review Commons**.
Authors do **not** need to have already posted a preprint for an initial evaluation —
during submission they can indicate an existing preprint, or ask eLife to deposit the
work to bioRxiv/medRxiv. A LaTeX template is available.

### Submit your research — five stages

1. **Submission** — submit directly or transfer a preprint from bioRxiv/medRxiv.
2. **Peer review** — editors decide whether the submission is suitable for expert review.
3. **Reviewed Preprint published** — published **two weeks** after reviews are sent.
4. **Revise when ready** — revise at your own pace based on reviewer feedback.
5. **Version of Record** — you decide when to end the process; the VOR is widely
   indexed after author proofing.

Feature Articles and Scientific Correspondence use the **legacy submission system**.

## Article Types

### Research Articles
No maximum length, but **try not to exceed 5,000 words** in the main text (excluding
Materials and Methods, References, and Figure legends). **No limit on display items.**
Usual structure: Introduction; Results; Discussion; Materials and Methods (or Methods);
Acknowledgements; References; Figures (legend below each); Tables. A Methods/Model
section may appear after the Introduction where sensible.

### Short Reports
Usually **≤1,500 words** main text (excluding M&M, References, Figure legends), with
**no more than three or four** main display items. More format flexibility (e.g. a
combined Results and Discussion section).

### Tools and Resources
Must fully describe the tool/resource. Major datasets must be publicly deposited
(barring strong ethical/legal reasons); code must conform to the Open Source Definition
and be deposited in a public repository; methods comprehensively described. Need not
report major new biological insight but must clearly enable such advances. New methods
must be **benchmarked against existing methods**; minor improvements are unlikely to be
reviewed. Follow Research Article or Short Report format.

### Research Advances
For developments that **directly build on** a prior eLife Research Article, Short
Report, or Tools and Resources article. Any length, any number of display items,
minimal introductory material; flexible format; no detailed M&M when methods match the
original. Abstract should reference the original, e.g. "Previously we showed that XXXX
(author, year). Here we show that YYYYY." Linked to the original; independently indexed
and citable.

### Replication Studies
Welcomed for new insight into previously published research (eLife or elsewhere).
Encouraged to work with original authors and summarise interactions in a cover letter.
Explain the rationale and methodology; justify any deviations; discuss implications;
pre-register where possible. Title: **"Replication Study: 'Title of original
article'"**. Abstract summarises the original's central findings and states which
part(s) are the focus. Outcome does not affect whether review proceeds; if reviewed it
must be a preprint. eLife may involve an original author in review.

### Review Articles
Primarily interested in reviews that synthesize a field in an original/thought-provoking
way and/or advance it. Not interested in routine summaries, much-reviewed fields, or
specialist-only reviews. Notes:
- **Preprint** required if review proceeds; bioRxiv/medRxiv may not accept reviews, so
  authors may need another DOI-and-versioning preprint server (direct transfer may be
  unavailable).
- **Cover letter** required: why you are a suitable author, how it moves the field
  forward, how it differs from recent reviews (with references).
- Usually **no more than five authors**.
- **Length:** no restriction; short perspectives ~2–3,000 words or longer scholarly
  reviews up to **10,000 words and up to 8 figures**.
- Authors should not unduly mention their own work; ideally conclude with open
  questions / outlook.
- **Figures:** do not reuse previously published figures without explicit copyright
  permission to reproduce under a CC-BY license.

### Article types that do NOT result in Reviewed Preprints

- **Editorials** — written by eLife editors or staff.
- **Insights** — commissioned by staff, always tied to a Research Article; written by
  experts to explain significance and remaining challenges.
- **Feature Articles** — fresh insight on topics of broad interest; advised to stay
  below **2,000 words, two display items, 20 references**; active/engaging style;
  peer reviewed at editors' discretion. Propose via features@elifesciences.org.
- **Scientific Correspondence** — challenges the central findings of an eLife paper (and
  the formal response). Must first contact the original corresponding author and include
  evidence of those efforts; submit **within a year** of the original; as of Feb 2021
  no longer covers papers in other journals. Title: **"Comment on 'Title of original
  article'"**. Usually **≤1,500 words**, **≤4** main display items; measured tone
  required. Uses the **legacy submission system**.

## Manuscript Preparation

*Source: full author guide, "Versions of Record" tab, which is where eLife now keeps all
manuscript-preparation guidance. Structural note: there is **no separate initial-submission
formatting section**. The strict requirements formally bite at VOR stage; the initial
submission is essentially your preprint plus submission-form metadata.*

### Title
"Succinct (we suggest **up to 120 characters**), with a clear indication of the biological
system under investigation (if appropriate), avoiding abbreviations and unfamiliar acronyms."
**Two-part titles should be avoided for research papers.** Suggested, not enforced.

> Our title, *Information distortion in likelihood approximations for macroscopic ion-channel
> currents*, is **88 characters**. Within the guidance, and not two-part.

### Abstract
"**Normally 150–200 words.**" Subheadings discouraged (except medical submissions, which may
run to ~250 words with Background/Methods/Results/Conclusions/Funding). "Clear, measured, and
concise." **If the biological system is not in the title, it must be in the abstract.**
See the contradiction note at the top: the LaTeX template says 150. Target 150.

### Impact Statement — STILL REQUIRED
The eLife assessment did **not** replace it. It is a **submission-form field**, not part of the
article file (the LaTeX class has no `\impactstatement` macro).

- **One sentence, typically 15–30 words.**
- Summarises the single most important finding.
- Must **complement, not repeat, the title**.
- **Third person.** No "we", no "our".
- Do **not** open with "We show…", "This study…", "Our work…".
- Avoid acronyms not familiar to a broad readership.

### eLife digest
A plain-language summary, **only at Version of Record stage** and only for some papers. Opt in
by answering "Yes" to the digest query when submitting the VOR; the Features Team then sends
four questions. Not written at initial submission.

### Required sections and order
"Introduction; Results; Discussion; Materials and Methods." Fuller running order for Research
Articles: Introduction; Results; Discussion; Materials and Methods; Acknowledgements;
References; Figures (legend below each); Tables. **A Methods or Model section may appear after
the Introduction where it makes sense.** Optional: an "**Ideas and Speculation**" subsection
inside the Discussion.

Other preparation rules: personal communications inline; manufacturers named on first mention;
footnotes in parentheses in the text (no footnote list); tables editable and inside the article
file (tables-as-figures not accepted); appendices at the end or as a second Article File, with
no separate reference list; supplementary files labelled "Supplementary File 1, 2, …".

### LaTeX specifics
Upload the `.tex` plus `.bib` and any `.bst`/`.sty` as **"LaTeX Support Files"**, and a clean
compiled PDF as a Related Manuscript File.

> **Binding constraint for our figures:** "LaTeX code uploaded to our submission system
> **cannot make use of packages such as `subfig` or `subcaption`**, so **all figures with
> panels must be collected into single figure files prior to upload**." Our multi-panel
> figures must each be one file.

The template ships `vancouver-elife.bst`, so the existing bibliography setup is fine.

### References
"**Authors do not need to spend time formatting their references.**" APA is the closest match
to eLife style, but any style is accepted and eLife reformats in production. **Provide a DOI
wherever possible.** Datasets, program code, and previously published methods must all be cited
and appear in the reference list with a persistent identifier. Software citations must include a
**version**. Preprints are citable.

### Author contributions, competing interests, funding
- **Author contributions:** required, using the **CRediT** taxonomy, entered per co-author.
- **Competing interests:** required **for each author**; ICMJE window is the **36 months** before
  submission; full ICMJE forms only for medicine submissions.
- **Funding:** entered in the submission form with **grant DOIs or reference numbers** and the
  authors tied to each source. "**Do not include information about direct funding in the
  acknowledgements**" (avoids duplication). State whether funders were involved in design, data
  collection, interpretation, or the decision to submit.

### Key Resources Table
"Where appropriate, and **especially for studies including bench research**", at the **VOR**
stage, at the very beginning of Materials and Methods. A purely computational paper is not
obviously forced into a full KRT, but a **software row with version and RRID** is expected.

### Cover letter — NOT a general requirement
The current guide contains **no general cover-letter requirement** for Research Articles, Short
Reports, or Tools and Resources. It is required only for **Review Articles** and **Replication
Studies**, and it is the place to request **data-availability exceptions** or note
**preregistration**. There are also **no presubmission enquiries**.

## Reviewed Preprints

Once peer-reviewed, eLife publishes a **Reviewed Preprint**: the full preprint text, the
eLife assessment, the public peer reviews, and (if the authors wish) an author response.
It receives a **DOI and eLife citation**, and is typically posted **within two weeks** of
the conclusion of peer review (authors can request earlier posting).

Authors should alert the editorial office as soon as possible on receiving the public
reviews if there are **factual errors** or other significant concerns, so errors can be
corrected and concerns discussed before posting.

On receiving the reviews, authors can:

- **Submit a revised version** — eLife decides whether to re-review, then publishes an
  updated Reviewed Preprint with updated reviews and assessment as appropriate. Multiple
  rounds of revision/re-review are discouraged. Before resubmitting, authors must **post
  a new version of their preprint**, and include a **response to the previous letter**
  and a **tracked-changes version** of the article alongside updated files.
- **Do nothing further** — if comfortable with the current version, reviews, and
  assessment; the Reviewed Preprint is citable in its own right.
- **Proceed to a VOR without changes** — though most authors revise first.

Additional points:

- Authors are referred to *"Ten common statistical mistakes to watch out for when writing
  or reviewing a manuscript"* and should follow eLife's **MDAR Checklist** for
  interpretation/replication.
- **Author-explainer videos (optional):** short engaging overviews of the paper or
  individual figures.
  - ≤5 minutes per video; presenter introduces themself; plain language.
  - Clear audible narration (external mic, minimal background noise); good lighting,
    stable camera, appropriate background; consider subtitles/transcripts.
  - May combine presenter footage with animations/figures; **no copyrighted or
    inappropriate material**.
  - High-level videos: don't repeat the abstract — explain motivation, methodology, main
    discoveries and implications. Labelled **Video 1, Video 2, …**.
  - Figure-specific videos: the figure legend must still be complete on its own. Labelled
    **"Figure 1—video 1"**, etc., each with a concise title/legend at the end of the
    article file.
  - Formats: **AVI, WMV, MOV, MP4, or H264**; upload as **"Rich Media"** file type.
  - Reviewed Preprints do **not** support in-line videos; the VOR does.
  - Published **CC-BY** unless otherwise noted. Longer talks (e.g. conference) can be
    hosted externally and linked.
- **Striking images:** encouraged, preferably colour, **landscape**, ≥**1800 × 900 px**,
  **PNG/TIFF/JPEG** preferred; no labels/text; ≤1–2 panels; must be usable under **CC-BY**
  (so **no AI-generated images**). Upload as **"Potential Striking Image"** with a short
  caption.
- Figures may be **screened** to ensure they have not been adjusted in ways that could
  mislead.
- Reviewed Preprints are distributed under **CC-BY** (except where noted), per the BOAI
  open-access definition.
- eLife reserves the right to **suspend or terminate** review where it cannot provide
  high-quality public reviews.

## Data Availability

**Principle.** eLife requires all data associated with an article to be **freely and
widely available**, in the most useful formats and per relevant reporting standards,
unless there are compelling legal/ethical reasons to restrict access. Data provision
should comply with **FAIR** principles (Findable, Accessible, Interoperable, Reusable).

Authors must make all original data supporting (or needed to reproduce) the claims
available in the manuscript text/tables/figures/supplements, or (recommended) in a
**trusted digital repository** — including all variables, treatment conditions, and
observations. Authors must give a full account of the materials and procedures used to
collect, pre-process, clean, generate and analyse the data so it can be independently
reproduced, and must ensure **data provenance for at least 10 years**.

- **Timing.** eLife considers work published when posted as a preprint and expects
  reviewed preprints to meet these standards (acknowledging preprint infrastructure is
  still evolving).
- **Data Availability Statement.** Every eLife paper must contain one, detailing how all
  relevant data are made available, including accession codes/identifiers.
- **Laws and ethics.** Data sharing must comply with legal requirements and institutional
  standards; eLife policy never supersedes legal/ethical standards.
- **Restricted access.** If access must be limited (privacy, legal, ethical), reasons must
  be stated in the statement and a specific exemption granted by the handling editors.
  *"Data will be made available upon request" is not acceptable* — a precise mechanism
  (contacts, timeline) is required. Proprietary data needs prior owner agreements,
  clearly stated owners, and editor awareness of constraints. Publicly-available-but-
  restricted data: describe restrictions. In all restricted cases, data should be
  available to editors/reviewers unless legally/ethically impossible (disclose any such
  limits at submission).
- **External/unpublished data.** Honour all access agreements, embargoes, and citation
  norms (e.g. **Fort Lauderdale** and **Toronto** agreements for genomic data); obtain
  and state permissions in the cover letter; cite lab/website/accession numbers.
- **Mechanisms.** Prefer trusted institutional/third-party repositories with persistent
  unique identifiers and long-term preservation. Author-maintained websites are generally
  non-compliant (exceptions only when the sole option). See recommended repositories on
  the **FAIRsharing Resource**; **DataSeer** gives a free sharing report.
- **Exceptions** are requested in the cover letter at submission; rare unshareable cases
  are considered expeditiously, sometimes via confidential sharing with editors.
- **Ongoing commitment.** Readers struggling to obtain data should contact
  editorial@elifesciences.org; unresolved cases may prompt a published statement.

### Specific data repositories

- **Small-molecule crystal structures:** `.cif` file, probability-ellipsoid figure,
  structure factors with IUCR CheckCIF output → **Cambridge Structural Database**.
- **Macromolecules:** validation reports + coordinates + structure factors/intensities →
  **wwPDB** or **BMRB**; EM density maps/coordinates → **EMDB** (annotated for immediate
  release).
- **Sequences:** reference-genome reads and assembly; novel short sequences (epitopes,
  domains, markers, haplotypes) with surrounding sequence → **GenBank / DDBJ / ENA**;
  DNA/RNA sequencing → **NCBI Trace Archive** or **SRA**.
- **Genetic variants:** **dbSNP, dbVar, or EVA**.
- **Proteins / proteomics:** **UniProt** / **PRIDE**.
- **Gene expression:** **GEO** or **ArrayExpress**.
- **Human genotype/phenotype:** **dbGaP** or **EGA** (unless private/sensitive — notify
  editors in advance and explain the restriction).
- **No domain-specific archive:** **Dryad, Dataverse, or OSF**; full catalogue at the
  **FAIRsharing Resource**.

## Figure / image specifications

From eLife's figure-preparation guidance:

- **Formats:** TIFF, JPG, EPS, AI, PDF.
- **Resolution:** minimum **300 dpi**; minimum physical width **10 cm**. Full-page-width
  figures must be supplied at minimum width **20 cm** and remain at 300 dpi.
- **TIFF:** provide in **8-bit** format.
- **PDF:** must be a **vector** image format.
- **Color space:** **RGB**; minimise whitespace around figures (especially in PDFs).
- Turn off any **hidden layers** before uploading to avoid processing issues.

## Revised Submissions

- Revise **when ready**. Alert eLife if you need **longer than 12 months**.
- Before resubmitting, **update the preprint on the preprint server** to match the revision, and
  give its URL, DOI and version number in the submission form. On bioRxiv/medRxiv choose the
  option indicating the work has **not** been published in a journal.
- Upload a **response to the previous letter**; a **tracked-changes** article file is recommended.
- The revision usually goes back to the reviewers. If you are not seeking an updated eLife
  Assessment, or the changes are minor, say so at the top of the Author Response and the editors
  can assess without the reviewers.
- Version 2 of the Reviewed Preprint carries the revised preprint, updated Assessment and Public
  Reviews, and the author response. Authors get a **two-week window** to flag factual errors.
- "We are always willing to consider the **first** round of revisions; requests to consider
  subsequent rounds will need to be justified."

## Versions of Record

- The authors ask for the latest Reviewed Preprint to be declared the VOR. **Most do this at
  version 2** (after one round of revision and re-review).
- The VOR must be **requested within 12 months** of the Reviewed Preprint being published.
- It is identical to the latest Reviewed Preprint apart from minor non-scientific changes, and is
  "expected to meet **more stringent policies and standards around ethics and data availability**".
- Publication of the VOR "marks the formal end of eLife's peer review process ... and is therefore
  akin to a traditional journal paper". There is a VOR checklist.
- **Indexing (important):** Reviewed Preprints are indexed by **Google Scholar and Scopus only**.
  VORs are indexed by DOAJ, Google Scholar, **PubMed**, PubMed Central, Europe PMC, Web of Science
  and OpenAlex. **PubMed indexing requires the VOR.**
- DOIs: a version DOI (cite this), an umbrella DOI resolving to the latest version (for grants and
  CVs), and separate DOIs for the peer reviews.

## Publishing Policies

The parts that bind a computational methods paper:

- **Code and reproducibility.** Authors "must provide **program code, algorithms, scripts for
  statistical packages, and other documentation sufficient to allow an informed researcher to
  precisely reproduce all published results**."
- **Software (binding).** If new software or a new algorithm is central, authors must confirm the
  software **conforms to the Open Source Definition**, is **deposited in an appropriate public
  repository**, and is released under an **open source license**. Version control (GitHub/GitLab)
  encouraged; **eLife archives GitHub-hosted code at Software Heritage**. Binary files should be
  minimised, ideally **under 50 MB**, avoiding files over 100 MB.
- **RRIDs** encouraged in Materials and Methods for tools (software or databases), e.g.
  `RRID:SCR_007358`.
- **Use of AI tools (binding).** Generative AI "cannot be listed as or serve as co-authors of a
  paper, and **their use must be clearly described in the Materials and Methods section**."
  AI-generated images cannot be used as striking images (not CC-BY compatible).
- **MDAR.** eLife endorses the **MDAR framework**; a completed **MDAR checklist must be uploaded
  and is published as a supplementary file**.
- **Reusable and open methods.** Do not bury detailed methods in supplements or on lab websites.
  Citing a method instead of describing it is allowed only if the citation describes a very similar
  or identical method, is detailed enough to reproduce it, and is **open access**. All
  modifications must be described.
- Also in this tab: COPE ethics, editorial independence, misconduct, **appeals** (within one month,
  on grounds of factual error, competing interest, or new data; one appeal per version), authorship
  (CRediT), ORCID, funder OA compliance (CC-BY, auto-deposit to PMC at VOR), plagiarism screening
  (iThenticate), preregistration, inclusive language, licensing, nomenclature.

## Publication Fees

- **$3,750** for submissions from **1 July 2026**, charged when the preprint is sent for peer
  review. 20% VAT added if invoiced to an individual rather than a business.
  *(eLife's own FAQ tab still says $3,000. It is stale.)*
- Covers initial evaluation, staff checks, peer review, publication of the Reviewed Preprint,
  re-review, subsequent versions, and the VOR on request.
- **Applies to** Research Articles, Short Reports, **Tools and Resources**, Research Advances,
  Replication Studies, Review Articles.
- **Full waivers available**, applied for during submission with a brief justification. **No strict
  criteria** (lack of funder or institutional support, terminated grant, career stage, geographic
  location, under-funded discipline all count). **Confidential:** editors and reviewers have no
  access to waiver information, so it cannot affect the outcome.

## Media Policy

- **No embargo.** Because eLife only reviews work already public as preprints, "we do not release
  our content under embargo ... journalists can write and publish articles about an eLife paper or
  Reviewed Preprint at any time".
- Authors are free to present findings at meetings, speak to the media at any time, and share the
  preprint with journalists.

## FAQs (the useful ones)

- **Scoop protection: yes.** eLife "will not decline to peer review your submission on the grounds
  that it lacks novelty because a paper on a similar topic has been published ... or posted as a
  preprint". Applies to competing preprints posted later, and to similar manuscripts already under
  review at eLife.
- **You can take the paper elsewhere after eLife review**, including taking eLife's reviews to
  another journal. "It is their paper, not ours." A VOR, however, ends the process.
- **No presubmission enquiries. No preprint required before submission** (eLife will deposit to
  bioRxiv for you if it proceeds with review).
- eLife does **not** track or support the impact factor (DORA co-founder).
- **Authors cannot switch from the pre-2023 model** to the new one.

## Publisher

eLife Sciences Publications Ltd, 95 Regent Street, Cambridge, CB2 1AW, UK.
editorial@elifesciences.org · ISSN 2050-084X · content under CC-BY except where noted.
