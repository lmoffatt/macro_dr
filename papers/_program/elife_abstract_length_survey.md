# How long a published eLife abstract actually is

Measured 2026-08-11. Replaces the undated `n = 41` figure that `00_abstract.md` and the retired
`SPINE.md` both quote and that no data file ever backed. The old numbers hold up: it reported median
197 and 39% over 200, against 195 and 44.8% here on ten times the sample.

Script `elife_abstract_survey.py`, per-article data `elife_abstract_survey_data.csv`, both beside
this file. Re-runnable; the API is public and unauthenticated. Same 413 article ids as
`elife_main_text_length_survey.md`, so the two surveys describe one sample and can be joined on `id`.

## Method

The `abstract` block of each article's API record. Structured-abstract section titles ("Background",
"Methods") are dropped; eLife's separate one-sentence impact statement is a different field and is
not counted. All 413 parsed, none missing.

## The distribution

| | median | p25 | p75 | p90 | max |
|---|---|---|---|---|---|
| all 413 | 195 | 157 | 225 | 252 | 310 |
| research articles (295) | **204** | 174 | 232 | 254 | 310 |
| tools and resources (118) | 157 | 148 | 198 | 229 | 289 |

Where a candidate length lands, research articles only:

| abstract | 150 | 200 | 220 | 239 | 250 | 260 |
|---|---|---|---|---|---|---|
| percentile | 11 | 47 | 65 | 80 | 85 | 92 |

## What this settles

**The 150-word rule is a dead letter.** `elife.cls` v1.11 says "no more than 150 words" and the live
author guide says 150 to 200 (`elife-author-instructions.md`). 89% of published research articles in
these two subject areas exceed 150, and 47% exceed 200. This is the same pattern the main-text survey
found for the 5,000-word advisory, and the same conclusion applies: **the advisory is a median, not a
wall.**

**Article type matters here, unlike main-text length.** Tools and Resources abstracts run 47 words
shorter at the median and sit much closer to the stated rule. Main-text length showed no such split.

**The tail is short.** Past 260 words an abstract is in the top 8% of research articles; past 300
there are two examples in 413. Nothing in the sample is longer than 310.

## Consequence for `00_abstract.md`

The brief's constraint, 200 to 220, sits at the 47th to 65th percentile and needs no revision on the
evidence. What needs revision is that no version of the section since 2026-08-01 has been inside it,
so the number is reopened by every review round without ever binding. Either the constraint moves and
is written down with a percentile beside it, or the section is cut to it. That decision is the
abstract's, not this file's.
