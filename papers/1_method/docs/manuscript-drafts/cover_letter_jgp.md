# Carta a JGP — borrador 2026-09-04

Esqueleto: la v5 de eLife (congelada en ../elife-submitted-20260828/cover_letter.md) con el
párrafo de venue rehecho para JGP y los requisitos de su guía (conceptual advance, relacionados,
conflictos, datos accesibles para revisión, sugerencias). Dos [PENDIENTE]: editor sugerido y DOI
del preprint. El .tex/.pdf se compilan cuando se resuelvan.

---

Dear Editors,

I am submitting "Likelihood approximations distort the ion channel kinetic information in
macroscopic currents" for consideration as an Article.

We have reached the limit of what a macroscopic current shows by eye. The flip state we proposed
for P2X2 in 2007 was still on the visible side, a delay before the current rises, later confirmed.
The questions now asked of these receptors, how subunits couple and whether they do so
symmetrically, leave no such mark: they are settled by ratios of evidence between kinetic schemes,
and whether the likelihoods behind those ratios report the information a recording actually holds
had never been measured. That measurement is the conceptual advance of this manuscript. At known
parameters a correct likelihood must satisfy two classical identities, a zero-mean score and the
information equality, and channel gating can be simulated exactly even where its likelihood cannot
be evaluated exactly. We measure both identities for eight likelihoods, from least squares to a
filter conditioning each acquisition interval on both of its endpoints, over 210 design cells with
ten thousand simulated recordings each. Every likelihood that leaves the gating correlation
unmodelled misreports its information, by up to an order of magnitude, and no rescaling repairs
it; the published control of our own P2X2 comparison is among them, accumulating about twenty
times the information it holds, and the recursion is what returns that factor to one. The
boundary-conditioned member is calibrated over almost the whole plane, its own limit is measured
rather than assumed, and it departs exactly where the record begins to resolve single openings,
the frontier with the microscopic regime. For an experimentalist the map has one reading that
does not depend on where a threshold is drawn: no measured region gives both a trustworthy
least-squares error bar and information about the unitary current.

The study is entirely simulation on a two-state scheme, and we believe JGP is its natural home.
The journal publishes theoretical work grounded on established experimental evidence, and the
closest precedent in kind is recent and in these pages: the identifiability analysis of Benndorf
and Schulz (2023), likewise a statistics-of-inference study on simulated recordings. The record
this work grew from is also here: the flip state was reported in this journal. The calibrated
likelihood, its score and Fisher information, the exact simulator and the distortion diagnostic
are archived at Zenodo (code 10.5281/zenodo.22167744, data 10.5281/zenodo.22168409, the macroir
library with R and Python bindings 10.5281/zenodo.22168263) and are accessible to editors and
reviewers; a reader can run the same two identities on their own scheme.

A preliminary version of this work was deposited in bioRxiv [PENDIENTE: DOI y fecha]. No related
manuscript is under consideration elsewhere, and we declare no conflicts of interest.

Should it help, we suggest [PENDIENTE: editor, de la lista actual de JGP] as editor, and as
reviewers Feng Qin (SUNY Buffalo), Klaus Benndorf (Jena), Gary Mirams (Nottingham), Lucia
Sivilotti (UCL), Andrew Plested (HU Berlin) and Colin Kinz-Thompson (Rutgers). We request no
exclusions.

Yours sincerely,

Luciano Moffatt
INQUIMAE-CONICET, Facultad de Ciencias Exactas y Naturales, Universidad de Buenos Aires
