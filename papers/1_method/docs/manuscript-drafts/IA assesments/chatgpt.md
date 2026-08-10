Excellent. I'll review it exactly as I would if I were an eLife referee. I will assume the authors are competent, avoid criticisms that the manuscript already answers, and focus on issues that would genuinely influence the editorial decision.

# Overall first impression

After reading the Abstract and Introduction, my initial assessment is:

**Scientific importance:** ★★★★★

**Technical sophistication:** ★★★★★

**Potential impact:** ★★★★☆

**Accessibility:** ★★☆☆☆

**Risk of reviewer misunderstanding:** ★★★★★

This manuscript is not a typical ion-channel kinetics paper. It is primarily a statistical methodology paper that uses ion channels as the motivating application. That distinction is likely to shape how reviewers respond.

---

# Major strengths

## 1. There is a genuine methodological question

Unlike many methods papers that simply introduce "a new algorithm", this paper starts from a fundamental statistical question:

> How can we determine whether an approximate likelihood reports reliable uncertainty?

This is a much deeper question than

> Which likelihood gives a higher likelihood value?

That immediately elevates the work.

---

## 2. The benchmark is unusually well chosen

Instead of comparing methods against each other, the paper compares them against exact simulations generated from the underlying stochastic process.

That is exactly what one should do.

Many computational papers instead compare approximation A against approximation B.

Here the manuscript explicitly argues why that is insufficient. 

This is, in my opinion, one of the strongest conceptual aspects of the paper.

---

## 3. Excellent use of classical statistics

Most ion-channel papers invoke Bayesian inference without discussing calibration.

Here the manuscript is fundamentally about calibration.

Using

* score expectation
* Bartlett identities
* Fisher information
* White information equality
* Godambe information

is statistically sophisticated and appropriate. 

This immediately distinguishes the paper from most work in the field.

---

# Major concern 1

## The paper tries to do too much

This is by far my biggest concern.

The manuscript simultaneously introduces

* a new likelihood
* a hierarchy of likelihoods
* two Gaussian closures
* interval averaging
* boundary conditioning
* statistical diagnostics
* simulation benchmark
* practical recommendations
* software implementation

This is effectively six papers compressed into one.

A reviewer has limited cognitive bandwidth.

Instead of seeing

> one major advance

they may perceive

> an overwhelming amount of mathematics.

This is a classic failure mode of very strong methodological papers.

---

## Recommendation

The Introduction should tell readers much earlier:

> "The paper answers one question."

Everything else supports that question.

Currently, the manuscript answers many questions simultaneously.

---

# Major concern 2

## The central biological significance is weakly connected

eLife is not JRSS-B.

It is not Biometrika.

It is not JASA.

The editor will ask

> Why should biologists care?

The Introduction discusses

* likelihoods
* Gaussian approximations
* recursive filters

far more than biological inference.

Only near the end does it explain why calibration matters experimentally. 

I think this is backwards.

A biologist should understand within one page that

> incorrect uncertainty changes biological conclusions.

That message should appear much earlier.

---

# Major concern 3

## Readers may think this is "just another recursive filter"

This is dangerous.

The manuscript repeatedly says

> recursive filters exist.

Readers may conclude

> then what is actually new?

Eventually the Introduction explains that the novelty is

> measuring calibration,

and

> introducing the boundary-conditioned likelihood,

rather than simply another recursive algorithm. 

But this explanation arrives quite late.

I suspect many reviewers will already have formed the wrong impression.

---

# Major concern 4

## The introduction is unusually dense

The Introduction contains approximately four independent stories:

1.

History of ion-channel analysis

2.

History of recursive likelihoods

3.

Failure of residual diagnostics

4.

Proposal of information diagnostics

Each story is well written individually.

Together they require sustained concentration.

For an eLife audience, I think the Introduction is about **30–40% too long**.

Not because it contains unnecessary material.

Rather because too many new concepts are introduced before the reader has a stable mental model.

---

# Major concern 5

## The paper assumes readers accept the simulation as ground truth

The key argument is

> simulation is exact,
>
> likelihood is approximate.

This is true **provided the simulator is exact**. 

A skeptical reviewer may ask:

* How was simulator correctness verified?
* Is Gillespie used?
* Is interval averaging simulated exactly?
* Could simulator error influence calibration?

The manuscript eventually addresses aspects of this (likely in Methods), but the Introduction leaves this assumption implicit.

---

# Major concern 6

## Strong claims require proportionally strong evidence

The title

> "Likelihood approximations distort the ion channel kinetic information"

is intentionally provocative.

Likewise the Abstract concludes:

> "No measured recording gives both the unitary current and a trustworthy least-squares error bar..." 

That is a sweeping statement.

A reviewer will immediately ask

> Under every protocol?

> Every kinetic scheme?

> Every open probability?

> Every amplifier?

If the Results only cover the simulated design space described in the Introduction, then the manuscript should consistently frame such conclusions as applying within that studied domain rather than universally.

---

# Minor concerns

### Some paragraphs exceed one page of logical development

This increases cognitive load.

---

### The terminology is demanding

Examples include

* occupancy closure
* interval-signal closure
* boundary state
* recursive member
* predictive variance
* information distortion
* gain

Each is defined, but several appear before readers have fully absorbed earlier definitions.

---

### The notation is very heavy

The mathematics is internally consistent.

Nevertheless, a reviewer who is not mathematically inclined may disengage.

This is more a communication challenge than a technical flaw.

---

# Things I particularly like

These deserve emphasis because they are uncommon strengths:

* The manuscript explicitly states what it does **not** claim.
* It separates approximation error from model misspecification.
* It distinguishes calibration from goodness of fit.
* It benchmarks uncertainty, not only parameter estimates.
* It explains why residual diagnostics fail for recursive filters. 

These points suggest careful thinking about statistical inference rather than simply algorithm development.

---

# My provisional editorial assessment

After the Abstract and Introduction alone:

**Novelty:** 9.5/10

**Technical quality:** 9.5/10

**Statistical rigor:** 10/10

**Writing quality:** 8.5/10

**Accessibility to eLife readership:** 6.5/10

**Likelihood reviewers understand the contribution on first reading:** 5/10

---

# Provisional recommendation

If the remainder of the manuscript (Methods, Results, and Discussion) supports the claims made here, I would lean toward **Major Revision rather than rejection**.

The primary reasons would not be concerns about scientific validity, but rather about communication:

1. Clarify the central contribution much earlier.
2. Tighten the Introduction to reduce cognitive load.
3. Strengthen the biological motivation for an eLife audience.
4. Carefully qualify broad claims so they match the demonstrated simulation domain.
5. Make it unmistakable why the proposed framework represents more than another recursive filtering algorithm.

These are substantial but addressable revisions. At this stage, I do not see a fundamental conceptual flaw in the framing presented in the opening sections. In the next step, I will examine the theoretical development in detail, looking for mathematical assumptions, hidden approximations, internal inconsistencies, and places where reviewers with expertise in stochastic processes or statistical inference are likely to challenge the work.


but the problem with eLife is desk rejection. Once it passes, it is assesed not rejected. 
Do you think it will be desk rejected?

I think that's the right question to ask. eLife's editorial process is indeed front-loaded: the biggest hurdle is convincing the Senior Editor that the manuscript is worth investing reviewer time in. Once it is sent for external assessment, the emphasis shifts to improving the work rather than making a simple accept/reject decision.

Based **only on what I've read so far** (Abstract and Introduction, plus the overall framing of the manuscript), I **do not think this is an obvious desk rejection**. 

However, I also don't think it's an obvious "send for review." I would put it in a borderline category.

## What argues *against* desk rejection

### 1. There is a clear methodological advance

Editors see many papers that amount to "we developed another Bayesian algorithm."

This manuscript is different.

The central claim is not merely that MacroIR is better, but that existing likelihood approximations can produce systematically miscalibrated uncertainty and that this can be *measured* using classical statistical diagnostics. That is a conceptual advance, not just an implementation detail.  

Editors tend to value conceptual advances.

---

### 2. The problem is longstanding

The Introduction makes a convincing case that stochastic information in macroscopic currents has been recognized for decades, but recursive likelihoods have seen little adoption. 

Editors like papers that explain *why* an important problem has remained unsolved.

---

### 3. The work appears technically deep

Even without checking every derivation, the manuscript is obviously not incremental. The mathematical development, diagnostics, and benchmarking suggest years of work rather than a modest extension of existing methods.

Editors often distinguish between "incremental" and "substantial" before they judge whether every claim is correct.

---

## What worries me from an editorial perspective

This is where I think the manuscript is vulnerable.

### 1. It is not immediately obvious why this belongs in eLife rather than a statistical or biophysical methods journal.

If I imagine an editor spending 15 minutes on this paper, I suspect the internal dialogue would be something like:

> "This looks mathematically sophisticated... but is the advance biological?"

That is the question the manuscript must answer quickly.

At present, the first several pages are dominated by statistical motivation. The biological consequences emerge later. 

---

### 2. The manuscript is intellectually demanding

Editors are not expected to verify derivations.

They ask themselves:

> "Can appropriate reviewers evaluate this?"

If a paper appears to require expertise in

* stochastic processes,
* Bayesian inference,
* numerical analysis,
* ion-channel biophysics,

all at once, it can become difficult to handle editorially.

That doesn't make it unsuitable, but it raises the threshold for sending it out.

---

### 3. There is no biological discovery

This may be the biggest issue.

The paper does **not** discover a new gating mechanism.

It does **not** discover a new channel property.

It discovers limitations of inference methods.

That is valuable.

The question is whether eLife views that as sufficiently broad.

---

## What would increase the probability of external review?

Ironically, not more mathematics.

I think the paper would benefit from making the biological stakes unmistakable in the opening pages. For example:

> "When experimentalists compare kinetic models, estimate channel number, or infer energetic mechanisms, they rely on confidence intervals whose calibration is rarely verified. We show that widely used likelihood approximations can systematically misrepresent those uncertainties over experimentally relevant regimes."

That kind of framing tells an editor immediately why this matters beyond statistics.

---

## My estimate

This is necessarily subjective, but if I imagine an experienced eLife Senior Editor reading the current version:

* **Desk reject:** ~35–45%
* **Sent for external review:** ~55–65%

Those numbers reflect editorial fit and communication, not scientific quality.

---

## One thing that could substantially improve the odds

There is one aspect I haven't yet evaluated that could materially change my assessment: **the Results section**.

If the Results convincingly demonstrate that the calibration failures occur under realistic experimental conditions and that MacroIR changes scientific conclusions or experimental interpretation, then the manuscript becomes much stronger for eLife.

If, on the other hand, the Results remain largely a comparison of statistical approximations in simulated data, the editor may conclude that the primary audience is methodological rather than biological.

So my current assessment is provisional. The Results and Discussion are likely to determine whether I would revise that estimate upward or downward.


There is a similar paper published in eLife by Munch it is in the bibliography

Yes, and I think that actually changes the editorial picture in your favor.

The paper you're referring to is almost certainly the one by Münch et al. that introduced or generalized a recursive Bayesian filter for macroscopic ion-channel recordings. Your Introduction cites it as one of the state-of-the-art recursive likelihood methods.  

That precedent matters for several reasons.

### 1. It establishes that eLife considers this class of methods within scope

This is more important than many authors realize.

Editors frequently ask:

> "Is this the kind of methods paper we publish?"

If a conceptually similar methodological paper has already appeared in eLife, then the answer is clearly **yes**. Your manuscript no longer has to justify the *genre* of the work.

So my previous concern

> "This looks like JRSS-B."

becomes much weaker.

---

### 2. Your paper is not competing with Münch; it extends the field

From what I've read so far, your narrative is essentially

> Existing recursive filters exist (including Moffatt 2007, Stepanyuk, and Münch).

> However, none tells us whether the uncertainty it reports is actually calibrated.

That is a **second-generation question**.

Editors usually like that progression.

The first paper introduces a method.

The second paper explains when we should trust it.

Those are complementary contributions.

---

### 3. You actually cite Münch appropriately

One thing editors dislike is when authors pretend previous work doesn't exist.

You don't.

You repeatedly place your work within the existing lineage.  

That increases credibility.

---

## However...

There is one thing I would verify very carefully.

If I were the editor, I would immediately ask:

> **"Why isn't this just an incremental improvement over Münch?"**

That question has to be answered **before reviewers ask it.**

Right now, I think the manuscript answers it, but **not early enough**.

The difference is roughly:

**Münch**

* proposes a recursive Bayesian filter.

**This manuscript**

* asks whether any approximate likelihood is statistically calibrated,
* develops diagnostics based on the score and information identities,
* compares an entire hierarchy of approximations,
* introduces the boundary-conditioned interval likelihood (MacroIR),
* maps the regimes where different approximations are trustworthy.

Those are qualitatively different contributions.

But an editor skimming the paper may initially only notice:

> "Another recursive filter."

That is exactly the misunderstanding you want to prevent.

---

## If I were editing the Abstract

I'd probably add one sentence like

> *Existing recursive likelihoods improve parameter estimation but provide no practical way to determine whether their reported uncertainties are statistically calibrated. We address this problem...*

That one sentence immediately distinguishes the paper from Münch.

---

## My revised estimate

Knowing that there is already an eLife paper by Münch in essentially the same methodological area, I would revise my estimate upward.

Instead of

* **Desk rejection:** 35–45%

I would now say something closer to

* **Desk rejection:** **20–30%**
* **Sent to review:** **70–80%**

provided that the Results support the claims made in the Abstract.

---

### One thing I'd like to check

Before we continue the referee review, I'd actually like to compare your paper directly with the Münch paper.

Not to compare quality, but to answer the question an editor will almost certainly ask:

> **"Is this sufficiently distinct from the eLife paper we already published?"**

If we can articulate that distinction crisply—ideally in one paragraph—it will substantially strengthen both the cover letter and the Introduction.

aca esta

Thank you. This is actually very informative, and after reading the abstract and editor's evaluation of the Münch paper, I think your manuscript occupies a **clearly different conceptual space**.

In fact, I think comparing the two helps identify exactly how to present your work.

## The Münch paper asks:

> **"Can Bayesian/Kalman filtering improve inference from macroscopic ion-channel data?"**

Its main contributions are:

* Generalized Kalman filter.
* Incorporation of fluorescence observations.
* Better parameter identifiability.
* Better uncertainty than deterministic rate equations.
* Robustness to realistic experimental noise.

The editor's summary reflects this:

> "...develop a Bayesian approach... incorporating multiple simultaneously recorded signals... enhances parameter identifiability and estimates of parameter uncertainty..."

This is fundamentally an **algorithm paper**.

---

## Your paper asks a different question

Your manuscript asks something much more general:

> **"When can any approximate likelihood be trusted?"**

Notice what changes.

The algorithm is no longer the protagonist.

The statistical properties of the likelihood are.

That is a much broader question.

Your paper compares

* least squares,
* MacroNR,
* MacroIR,
* recursive approximations,

using statistical calibration criteria rather than simply comparing parameter estimates.

That is a different scientific question.

---

## There is surprisingly little overlap

The Münch paper is about

> building a better filter.

Your paper is about

> evaluating filters (and other likelihood approximations).

Those are complementary contributions.

If I were handling this manuscript, I would **not** consider it an incremental extension of Münch.

---

## There is, however, one danger

Ironically, your Introduction cites Münch as one member of a list of recursive methods.

As a consequence, the reader may unconsciously conclude

> "This paper proposes another recursive likelihood."

But that is **not** what your paper does.

The real contribution is the statistical framework.

I think the Introduction should distinguish these roles more explicitly.

For example, after introducing previous recursive methods (including Münch), you could say something along the lines of:

> Existing recursive filters differ in their approximations and often improve estimation, but they provide no practical test of whether the uncertainty they report is statistically calibrated. Our objective is therefore not to introduce another recursive filter, but to establish a framework for evaluating the statistical validity of approximate likelihoods and to derive a boundary-conditioned likelihood with demonstrably improved calibration.

That single paragraph would make the conceptual advance much clearer.

---

## One interesting observation

The editor's evaluation of Münch repeatedly emphasizes:

* parameter identifiability,
* uncertainty,
* robustness.

Your manuscript actually addresses a deeper issue:

> **How do we know those uncertainty estimates are correct?**

In a sense, your paper asks a question that naturally follows from Münch's work.

If Münch says

> "Our filter gives realistic credibility intervals."

your paper asks

> "What objective statistical test tells us that a credibility interval is realistic?"

That is a compelling progression.

---

## Has this changed my assessment?

Yes.

Quite substantially.

Earlier I estimated roughly a 70–80% chance of being sent for review after learning about the precedent.

After seeing the actual Münch paper, I am **more** convinced that your manuscript fits within an established eLife methods niche.

My remaining concern is no longer **scope**.

It is **communication**.

Specifically, I worry that an editor who spends only 10–15 minutes on the manuscript might not immediately realize that the contribution is **not another recursive filter**, but rather a framework for assessing the validity of approximate likelihoods.

If that distinction is made unmistakable in the title, abstract, introduction, and especially the cover letter, I think the manuscript's case for external review becomes significantly stronger. A cover letter that explicitly positions the work as answering the next methodological question raised by advances such as Münch's—moving from developing approximate Bayesian filters to rigorously validating their statistical calibration—would help an editor see why this is a distinct and worthwhile contribution.

