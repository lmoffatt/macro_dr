# Reading baseline noise off a published current trace

Written 2026-07-27, after a literature sweep (13 agents, adversarially verified) established that
**typical instrumental noise per recording configuration is not published in usable form**. Every
whole-cell figure found is an amplifier floor with no cell, a model cell, or theory; for two-electrode
oocyte clamp there is no rms figure with a stated bandwidth at all; and for a true macropatch there is
none either. See `SOURCES.md` §2.

So the noise has to be read off published traces. This file is the recipe and its calibration.

## The recipe

1. **Vertical scale bar** → pA per pixel.
2. **Horizontal scale bar** + the sampling rate from Methods → **n, the number of samples that fall in
   one pixel column**. This is the step that makes the method work; see the calibration below.
3. **Ink stroke width**, in pixels, measured on **the trace's own pen**, on a near-vertical limb of the
   trace where the pen width is the horizontal extent. Use the scale bar only as a fallback.
   *Corrected 2026-07-27 after it went wrong in five papers out of eight, always in the same
   direction:* Harnau bar 3.09 px vs trace 2.44; Karasawa bar 2.0 vs trace 1.0 in one panel and 3.25
   in another panel of the same figure; Li bar 3.9–5.3 vs trace 2.80; Ivica bar 9 vs trace ≤3.76;
   Sokolov bar 12 vs trace 6; Dallas, read out of the PDF drawing commands, bar 0.648 pt vs trace
   0.300 pt. That last one would have turned a good 4.21 pA measurement into a spurious floor at
   0.92 pA. A fatter scale bar biases sigma **down**, which is the dangerous direction here.
4. Over a segment whose centreline slope is small, take the **median per-column vertical extent** of
   the trace, in pixels.
5. Then

       sigma = (extent - ink) * (pA per px) / d2(n)

   where `d2(n) = E[range of n standard normals]`, computed as `∫ (1 - F^n - (1-F)^n) dx`. Values:
   d2(2)=1.13, d2(5)=2.33, d2(7)=2.70, d2(10)=3.08, d2(20)=3.73, d2(50)=4.50, d2(85)=4.90, d2(100)=5.02.
6. Record the **filter type and bandwidth B** from Methods and convert to a power spectral density,
   `S = sigma^2 / B`, which is the only form comparable across papers. A sigma without its bandwidth
   is unusable, because patch-clamp noise is not white.

## Calibration against a record whose answer is known

Moffatt & Hume 2007 J Gen Physiol 130:183, outside-out patches from HEK293T, P2X2, −60 mV, 10 kHz
4-pole Bessel, digitized 50 kHz. The raw record of one patch is in
`projects/p2x2/data/experiments/Moffatt_Hume_2007_ATP_time.txt`, so the answer is known independently.

| source | scale bars | n per column | d2(n) | extent − ink | **sigma** |
|---|---|---|---|---|---|
| raw record, ns=1 points, first differences over 5 contiguous runs | — | — | — | — | **2.28 pA** |
| Fig. (NMDG external), panel A, late tail | 50 pA = 140 px; 50 ms = 380 px | 6.6 | 2.70 | 6.07 pA | **2.25 pA** |
| Fig. (long ATP application), panel A, plateau | 400 pA = 168 px; 500 ms = 293 px | 85.3 | 4.90 | 14.3 pA | **2.92 pA** |

Three numbers from the same study, same rig, same filter, different patches: 2.25, 2.28, 2.92 pA.
The figure-derived values bracket the data-derived one and the spread is patch-to-patch, about ±25 %.

**The 13-fold difference in n between the two panels is the real test.** A fixed "width ≈ 5 sigma"
rule of thumb would have given 1.21 pA for the first panel (47 % low) and 2.86 pA for the second
(2 % low). Counting the samples per pixel column is what makes one recipe work at both time bases.

## Check n before anything else

`n`, the samples per pixel column, decides whether the method applies at all. Measured range of
validity: **about n = 5 to a few hundred.**

- **n < 1**: the plotted trace is decimated, each column holds at most one sample, `d2(n)` is
  undefined and step 7 cannot be applied. One paper in the 2026-07-27 sweep failed this way (n = 0.49).
- **n in the thousands**: `d2` is still fine, it only grows as √(2 ln n), but one pixel column then
  spans a long stretch of record — at n = 22145 it is 1.1 s — so the drawn envelope contains mains
  hum and drift as well as in-band noise, and sigma comes out high for reasons that have nothing to
  do with the amplifier.
- The plotted bandwidth is not always the acquisition bandwidth. One paper states in Methods that an
  extra 1 kHz Gaussian filter was applied in Clampfit **for the figures only**. Read the figure
  preparation part of Methods, not just the acquisition part.

## Where it fails, and how it is biased

- **Sloping segments inflate it.** The same figure measured on its rising phase gives 3.43 pA instead
  of 2.25: within a column the signal itself moves. Only use segments with a small centreline slope,
  and report the slope.
- **A centreline-SD variant does not work.** Fitting the centreline and taking the residual SD gave
  3.86 and 7.86 pA on the same two panels, because the residual is dominated by real signal curvature,
  not noise. Use the per-column extent, which is local.
- **Floor at the ink width.** When sigma falls below roughly one pixel the trace is pure ink and the
  method returns zero. It therefore fails on clean single-channel records and works exactly where the
  literature hole is: whole cell and oocyte, where the noise is large and plainly visible.
- **Publication bias, downward.** Published traces are the best records of a set, not typical ones.
  Values obtained this way are lower bounds on the typical noise of that configuration and should be
  reported as such.
- Check the caption for averaging or additional digital filtering before taking the Methods bandwidth
  as the trace's bandwidth.

## Reproducing it

Images out of a PDF at 600 dpi with `pdfimages -png`; scale bars found as the longest contiguous
horizontal and vertical runs of dark pixels; ink width as the stroke thickness of the bar itself. The
working scripts for the calibration above are in the session scratchpad; they are short enough to
rewrite from this description, which is why they are described rather than committed.
