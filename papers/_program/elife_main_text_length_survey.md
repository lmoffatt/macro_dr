# How long a published eLife article actually is

Measured 2026-08-09. Fills the gap named in `elife-author-instructions.md`: the repo had surveyed
abstracts (n = 41), titles (n = 20), captions (7 papers) and figure supplements (n = 1,079), and had
never measured main-text length. The only number anywhere was eLife's transcribed advisory, "try not
to exceed 5,000 words in the main text (excluding Materials and Methods, References, and Figure
legends)".

Script `elife_length_survey.py`, per-article data `elife_length_survey_data.csv`, both beside this
file. Re-runnable; the API is public and unauthenticated.

## Method

413 Version-of-Record articles from the eLife API, the most recent in the two subject areas this
programme already samples: Structural Biology and Molecular Biophysics (n = 209) and Computational
and Systems Biology (n = 213), split by article type into research-article (n = 295) and
tools-resources (n = 118). 85% are 2024 or later; the distribution is stable when restricted to
2025+, so recency is not driving it.

Main text is counted eLife's own way: every top-level body section except Materials and Methods,
with figure legends, table content, equations and references excluded. Appendices are counted
separately and never inside the main text.

The count was cross-checked against a second, independent parser reading the JATS XML instead of the
JSON API, on five articles spanning 4.4k to 13.3k words: agreement within 0.1% to 2.0%, the XML
running slightly high because it picks up section titles.

## The distribution

| | median | p75 | p90 | p95 | max |
|---|---|---|---|---|---|
| main text, all 413 | **5,209** | 6,543 | 8,017 | 9,066 | 13,300 |
| research articles (295) | 5,298 | 6,554 | 8,087 | 9,102 | 13,300 |
| tools and resources (118) | 5,013 | 6,493 | 7,660 | 8,457 | 10,777 |
| Materials and Methods | 2,249 | | 4,147 | | |

55% of published articles exceed the 5,000-word advisory. 14% exceed 7,500, 2.4% exceed 10,000, and
0.5% exceed 13,000. Where a candidate target lands:

| main text | 5,000 | 6,000 | 7,000 | 8,000 | 9,000 | 10,000 | 11,515 |
|---|---|---|---|---|---|---|---|
| percentile | 45 | 66 | 82 | 90 | 95 | 98 | 99.3 |

## Section by section

Composition is as informative as the total, and it is where this manuscript is most off-norm.

| section | n | median | p75 | p90 | max |
|---|---|---|---|---|---|
| Introduction | 413 | **873** | 1,079 | 1,324 | 8,243 |
| Results | 377 | 3,066 | 3,976 | 4,941 | 10,533 |
| Discussion | 379 | **1,184** | 1,635 | 2,081 | 4,850 |
| Results and discussion (combined) | 32 | 3,977 | 5,017 | 7,278 | 8,256 |

eLife introductions are short: 873 words at the median, and only 10% reach 1,324. Discussions are
shorter still. Paper 1 currently runs Introduction 2,028 (2.3× the median, above p90) and Discussion
4,328 (3.7× the median, and within 500 words of the longest Discussion in the sample).

Methods placement: 395 of 409 articles put Materials and Methods last. Only 14 use eLife's
permission to place a Methods or Model section after the Introduction, so that route is available
but rare, and it inverts reading order for anything the Results depend on.

**The advisory is the median, not a wall.** This is the same pattern the abstract survey found, where
the transcribed 150-word rule sat 50 words below the measured median. The difference is that here the
tail is short: past about 9,000 words an article is in the top 5%, and past 13,300 there are no
examples in this sample at all.

**Article type buys nothing.** Tools and Resources is not a longer format in practice; its median
main text is 285 words *below* research articles and its maximum is 2,500 words below. Whatever D-1
is decided on, it should not be decided on length.

## Appendices

28.1% of articles carry at least one appendix (23 of those carry only figures, so 22.5% carry
appendix prose). Among the 116 that have any: median 1 appendix, p75 2, p90 5, max 11; median 1,146
words, p75 2,752, p90 6,480, max 13,791 (e94586). Appendix figures, among the 101 articles that have
any: median 6, p90 20, max 61.

Across all 413 articles, appendix prose exceeds 2,000 words in 9.7%, 5,000 words in 3.9%, 7,500 in
1.9% and 10,000 in 0.7%. Main text is no shorter in articles that use appendices (median 5,120 with,
5,277 without), so in practice an appendix is used to *add* material rather than to relocate it.
That is a fact about other authors' habits, not a constraint: eLife states no limit, and the only
stated rules are that appendices sit at the end or in a second Article File and carry no separate
reference list.

Methods is not a free dumping ground either. Its own distribution: median 2,249, p75 3,247, p90
4,147, p95 5,046, p99 6,738, max 11,580; 31% exceed 3,000 words and only 5.3% exceed 5,000. The
median Methods-to-main-text ratio is 0.44 (p90 0.81). Paper 1's Methods at 5,393 words is already
around p96, so it can absorb perhaps a thousand more words before it becomes remarkable in its own
right.

Appendix-heavy precedents, all in these two subject areas: e94586 (13,791 words / 11 appendices),
e89862 (13,258 / 2), e92497 (11,234 / 1), e79812 (8,832 / 8), e86365 (8,730 / 11), e100284
(8,426 / 1).

## The direct comparator

Münch 2022, `e62714`, the Bayesian Kalman filter for ion-channel data that this programme positions
against, measured the same way:

- main text **11,515** words (Introduction 1,842 + a combined Results and discussion 9,673)
- Materials and Methods 3,315 words
- 9 appendices carrying 5,889 words
- 16 figures, 0 figure supplements

Its main text is in the **top 0.7%** of this sample: only 3 of 413 articles reach it. So the
theory-led, appendix-carrying, very long architecture is precedented in exactly this niche and in
exactly this journal, and it is also rare. Aiming at Münch means aiming at the 99th percentile
deliberately, with the argument for why this paper is that kind of paper, and not drifting there.

## What this implies for paper 1

Counted the same way (Introduction + Theory + Diagnostics + Results prose + Discussion, captions and
Methods excluded), the manuscript is at roughly **25,100 words**, which is 4.8× the median, 2.8× the
p90 and 1.9× the largest article in the sample.

Reference points, in order of ambition:

- 5,209 — the median. Not reachable for this paper without cutting scope.
- 8,017 — p90. Reachable only if Diagnostics also leaves the body.
- 11,515 — Münch. Reachable with the measured triage: Theory to ~2,500 in the body with ~7,500 in an
  appendix, Results prose to ~3,200, Discussion to ~2,300, Diagnostics to ~1,600, Introduction to
  ~1,700, which sums to ~11,300.
- 13,300 — the largest article in the sample. Anything above this has no precedent here.

Methods is exempt from the count but its own norm is worth knowing: at 5,393 words it sits above the
p90 of 4,147, and Münch's is 3,315.
