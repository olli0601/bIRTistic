---
bibliography: amortised_decision_making.bib
---

```{=html}
<!--
Companion to amortised_decision_making.md. Render with pandoc + citeproc so the
References section is generated from the shared .bib:
  pandoc dev/literature.md --citeproc --bibliography dev/amortised_decision_making.bib \
    --mathjax --standalone -o dev/literature.html
-->
```

# Amortising real-time clinical and humanitarian decision-making

## 1. Introduction: Related work and positioning

*Companion to `amortised_decision_making.md`. Places the amortised predictive-probability-of-success (PPS) estimator developed there within the surrounding statistical and machine-learning literature. Written as a related-work section for a statistical journal; references are collected in the shared `amortised_decision_making.bib` and rendered below by citeproc.*

### 1.1 Overview

The method of this work amortises the **interim decision quantity** of a sequential trial — the predictive probability that the study will succeed, computed from the accumulating data — so that at each interim it is returned by a single forward pass of a pre-trained network rather than by a fresh nested Monte-Carlo computation. It therefore sits at the intersection of four literatures that have so far developed largely apart: (i) the *assurance / probability-of-success* tradition in clinical-trial design, which defines the target functional but computes it by simulation; (ii) *simulation-based (likelihood-free) and amortised Bayesian inference*, which amortises the **posterior** but not a downstream decision; (iii) *amortised and sequential Bayesian optimal experimental design (BOED)*, which amortises a **design policy** under an information- or utility-based objective; and (iv) the classical theory of *posterior contraction* and the *Bernstein–von Mises* theorem, which we use to make the amortised predictive's width honest as a function of cohort size. To our knowledge no existing work amortises the assurance/PPS functional itself, with conditional predictives that are calibrated for contraction over variable cohort sizes and variable problem dimension; that is the gap this work occupies.

### 1.2 Predictive probability of success and assurance

The quantity we estimate is the Bayesian **predictive probability of success (PPS)**, known in trial design as **assurance**: the prior- or posterior-predictive probability that a future analysis will meet its success criterion, averaging the frequentist success indicator over the current uncertainty about the parameters rather than fixing a single alternative as in a power calculation. The foundational reference is @ohagan2005assurance; the posterior-data variant used for interim monitoring — "will the trial succeed given what we have seen so far?" — is developed by, among others, @ibrahim2015bayesian. These quantities are used exactly as we use them: to size a trial, to stop early for futility, or to continue.

Crucially, the standard computation is **by simulation**: the PPS is obtained by drawing parameters from the posterior of the observed data, simulating the remaining cohort, re-fitting, and averaging the success indicator — a nested Monte-Carlo loop repeated at every interim and every design. A recent representative is the simulation-based Bayesian PPoS for interim monitoring with competing-event data [@micoli2024simulation]. This nested-simulation estimator is precisely the reference ("nested-MC") baseline our method is trained to reproduce at a fraction of the deployment cost; the contribution is not a new definition of PPS but the amortisation of its computation.

### 1.3. Value of information

The decision-theoretic sibling of assurance is the **expected value of sample information (EVSI)** and the broader value-of-information calculus [@ades2004expected; @strong2014estimating; @heath2018efficient]: the expected gain in decision utility from collecting further data, again an expectation over posterior-predictive futures. EVSI shares our computational pain point — it is classically a nested Monte-Carlo integral — and much of the modern EVSI literature is about cheaper surrogates (regression / GAM approximations, moment matching). Our estimator can be read as amortising the inner object these methods approximate, specialised to the success-indicator utility.

### 1.4. Simulation-based and amortised Bayesian inference

The engine we borrow is **simulation-based inference (SBI)**, also called likelihood-free inference: fit a neural conditional density (or ratio) on simulated parameter–data pairs so that inference for a new dataset is a forward pass. The field is surveyed by @cranmer2020frontier; the neural-posterior-estimation line runs @papamakarios2016fast to @greenberg2019automatic. The defining property is **amortisation**: front-load the cost into training, then reuse the network on any new observation without further simulation — the property our deployment exploits. The `BayesFlow` framework [@radev2020bayesflow] packages this for Bayesian workflows and is the closest engineering analogue to our training loop.

Two SBI tools are load-bearing in our work. First, **simulation-based calibration (SBC)** [@talts2018validating]: the self-consistency of the joint prior-predictive draw that we use, in its probability-integral-transform form, as the conditional-calibration diagnostic (PIT–KS). Second, the observation that SBI amortises the **posterior of the parameters**, whereas we amortise a **scalar decision functional of the posterior predictive** (the success probability), a different and lower-dimensional target that lets us bring the network's capacity to bear on the decision rather than on the full parameter vector.

### 1.5. Amortised and sequential Bayesian optimal experimental design (BOED)

The literature closest to ours in *spirit* — amortising a decision that would otherwise require per-step Bayesian computation — is amortised/sequential **BOED**. Its objective, however, differs: BOED amortises a **design policy** to maximise an information- or utility-based criterion (typically the expected information gain, EIG), whereas we amortise the **evaluation of a fixed success criterion** on a fixed accrual schedule. Key references:

- @foster2019variational — the variational-EIG foundation for gradient-based BOED.
- @foster2021deep (Deep Adaptive Design, DAD) — learns a history→design policy network so an adaptive design is a millisecond forward pass; the archetype of "amortise the sequential decision."
- @ivanova2021implicit (Implicit DAD) — extends DAD to the implicit-likelihood (simulator) setting.
- @blau2022optimizing — casts adaptive BOED as reinforcement learning.
- @huang2024amortized (Amortized Decision-Aware BED) — the most decision-focused member: amortises the design under a downstream decision utility rather than pure information gain. This is the nearest neighbour to our decision-first framing, but it still amortises the *design choice*, not the success probability at a pre-specified analysis.
- @bracher2025jadai (JADAI) — amortises design **and** posterior end-to-end.
- Review: @rainforth2024modern.

The distinction from our work is consistent across this list: BOED chooses *what experiment to run* to be maximally informative; we score *whether an ongoing trial, run on a fixed schedule, is on track to succeed*. The design is not a decision variable in our setting — the accrual is given — so the amortised object is the predictive success probability and its calibrated uncertainty, not a policy.

### 1.6. Exchangeable and set-invariant architectures

Because a cohort is an unordered set of exchangeable participants, the encoder must be permutation-invariant, and its pooling must transmit the *sample size* to the head so that the posterior width can contract. Our nested-DeepSets-with-cross-attention encoder descends from **DeepSets** [@zaheer2017deep] and the **Set Transformer** [@lee2019set], and from the use of exchangeable summary networks in likelihood-free inference [@chan2018likelihood]. Our contribution at the architecture level is minor and empirical — a ragged, mask-free 1/n mean-pool that delivers the Bernstein–von Mises precision, plus item cross-attention that prices a queried component from the others — but the invariance requirement, and the finding that the pooling normalisation governs whether contraction is learned, are the relevant points of contact.

### 1.7. Posterior contraction and Bernstein–von Mises

The statistical backbone of the width model is classical asymptotics: under regularity the posterior of a smooth functional is asymptotically normal with variance of order $n^{-1}$ at the semiparametric efficiency bound (the **Bernstein–von Mises** theorem; @vandervaart1998asymptotic; @kleijn2012bernstein for the misspecified case), so the posterior standard deviation contracts as $n^{-1/2}$. An amortised predictive used for an *interim* decision must reproduce this contraction — a fixed-width head is dishonest at large $n$ — and our parametric power-law / floor heads and the between-variance self-consistency correction are exactly a device to impose the Bernstein–von Mises rate on the network's output. We are not aware of prior amortised-inference work that makes posterior contraction an explicit training and calibration target; SBI calibration is typically assessed by SBC/coverage at fixed data size, and BOED does not model it.

### 1.8. Where this work sits

The following map summarises the four neighbouring literatures against the present method along the axes that distinguish them: what is amortised, what functional is targeted, whether uncertainty is calibrated for contraction, and whether the estimator generalises across cohort sizes and problem dimension.

| Line of work                 | Representative refs                                                                                | What is amortised                                       | Target functional                                      | Contraction-calibrated                                | Varies with n, m, J                      |
| ---------------------------- | -------------------------------------------------------------------------------------------------- | ------------------------------------------------------- | ------------------------------------------------------ | ----------------------------------------------------- | ---------------------------------------- |
| Assurance / PPS (simulation) | @ohagan2005assurance; @ibrahim2015bayesian; @micoli2024simulation                                  | nothing (nested MC each time)                           | PPS / assurance                                        | —                                                     | n fixed per run                          |
| EVSI / value of information  | @ades2004expected; @strong2014estimating; @heath2018efficient                                      | inner nested-MC surrogate                               | expected decision value                                | —                                                     | design-specific                          |
| Amortised SBI / NPE          | @cranmer2020frontier; @radev2020bayesflow; @papamakarios2016fast                                   | the posterior $p(\theta\mid x)$                         | parameter posterior                                    | SBC/coverage at fixed n                               | fixed data shape                         |
| Amortised / sequential BOED  | @foster2021deep; @ivanova2021implicit; @blau2022optimizing; @huang2024amortized; @bracher2025jadai | a design policy                                         | expected information gain / utility                    | —                                                     | sequential designs                       |
| **This work**                | `amortised_decision_making.md`                                                                     | **the decision functional** $p(H_1\mid x,z^{s})\to$ PPS | **direction-aware endpoint effect** $\rho$ and its PPS | **yes — BvM power-law head + PIT/SBC + affine shift** | **yes — n, m, and item set J amortised** |

Against this background the specific claims of the present work are:

1.  **Amortising the decision, not the posterior or the design.** The network returns the conditional predictive of the endpoint effect over posterior-predictive future cohorts, from which the PPS follows directly. This is the assurance/PPS functional of §1 computed at deployment by a forward pass — the amortisation that the simulation-based PPS literature has not taken, and a different target from both SBI (posterior) and BOED (design).
2.  **Contraction as a first-class target.** The predictive width is made to obey the Bernstein–von Mises $n^{-1/2}$ law through a parametric head and a self-consistency correction, with the fit validated on the multivariate-normal benchmark where the law is exact. Honest interim uncertainty is a requirement of the decision problem that neither SBI calibration (fixed-n coverage) nor BOED addresses.
3.  **A deployment-calibration stack.** Expanding-window head fine-tuning plus a per-item affine shift correct the prior-to-operational mislocation, and PIT / marginal KS diagnostics (an SBC descendant) certify conditional and marginal calibration — the practical machinery that turns a prior-trained network into a trustworthy interim tool.
4.  **Amortisation across cohort size and problem dimension.** A single network prices arbitrary observed/future cohort sizes ($n, m$) and, via item cross-attention, an arbitrary number of components / items $J$ (demonstrated over $J = 2$–$100$), going beyond the fixed-data-shape regime of standard amortised SBI.
5.  **A substantive application class.** The estimand is an item-level effect in ordinal item-response (partial-credit) models for humanitarian and social-science trials — a setting the assurance literature has treated only with summed scores and without amortisation.

The sharpest one-line positioning: existing work amortises either the *posterior* (SBI) or the *design* (BOED), and computes assurance/PPS by *simulation*; this work amortises *assurance/PPS itself*, with contraction-calibrated conditional predictives, across cohort sizes and item sets.

## 2. Methods

*The innovations of §1.8, as a methods outline. Each advances one axis of the "This work" row: the amortised target (2.1–2.2), amortisation across $n,m,J$ (2.3–2.4), contraction calibration (2.5), the deployment stack (2.6), and generalisation across endpoint functionals (2.7). Key points only; expanded in the full companion.*

### 2.1 Target: amortise the decision, not the posterior or the design

- Estimand is the PPS / assurance $P(H_1\mid x)=\mathbb{E}_{z\mid x}\big[\mathbf 1\{P(H_1\mid x,z)>\eta_H\}\big]$ — the predictive probability a future analysis meets its success criterion.
- The network returns the conditional predictive $p(\rho\mid x,z)$ of the endpoint effect over posterior-predictive future cohorts; the PPS follows by a forward pass. Amortises the *decision functional* — not the posterior (SBI), not a design policy (BOED).
- Non-reduction to today's knowledge: $\mathbb{E}_z[\mathbf 1\{\cdot>\eta_H\}]\neq\mathbf 1\{P(H_1\mid x)>\eta_H\}$ (Jensen); the integration over future data is the point of the object.
- Reference target = the nested-MC estimator (SVI/HMC); the net is trained to reproduce it at $\sim$one forward pass.

### 2.2 Endpoint effect $\rho$ from ordinal item-response models

- $\rho$ = direction-aware **item-level** effect in a partial-credit (ordinal IRT) model — the substantive application class (humanitarian / social-science trials and vaccine immunogenicity), beyond summed scores.
- One general functional: a per-condition **reduction** (mean $\mathbb{E}[y]$ or threshold $P(y\ge c)$) combined by a **compare** (level / fold / relative change / between-group difference) — recovering seroprotection rate, GMT fold-rise, relative change and arm difference from the same fitted category probabilities.
- Several $\rho$ carried on one Monte-Carlo sample $\Rightarrow$ **composite** decisions ("at least one of SCR / SPR / GMFR met") as the exact joint posterior probability, not a product of marginals.

### 2.3 Amortised workflow

- Three steps: (i) simulate labelled (cohort, $\rho$) pairs from the model's prior-predictive; (ii) learn a neural conditional predictive; (iii) deploy by forward pass on the observed cohort — cost front-loaded into training.
- Prior-coverage requirement: the wide training prior must envelope the operational posterior slice, or the head fine-tune (2.6) cannot repair the distribution shift.

### 2.4 Set-invariant architecture over participants and items

- Cohort = exchangeable set $\Rightarrow$ permutation-invariant encoder; pooling must transmit the sample size $n$ so the width can contract.
- Ragged, mask-free $1/n$ **mean-pool** over participants — the normalisation, not the architecture, is what makes contraction learnable (delivers the BvM precision).
- **Item cross-attention** prices a queried item from the others $\Rightarrow$ amortise over the item set $J$ (2–64 items) from item *features*, not identity.
- **Wide cumulative-exceedance token** $[\mathbf 1\{k\ge c\}]_c$: its mean-pool is the empirical category CDF — *sufficient* for the mean, any threshold, and group contrasts on one encoder.
- Continuous $\rho$ predictive via a **multi-quantile head** trained with the pinball loss; amortises arbitrary $n$ (observed) and $m$ (future).

### 2.5 Contraction as a first-class target (Bernstein–von Mises)

- Interim honesty requires the width to obey the BvM $n^{-1/2}$ law (posterior SD of a smooth functional at the efficiency bound); a fixed-width head is dishonest at large $n$.
- A parametric **power-law / floor head** imposes the rate; validated on the multivariate-normal benchmark where the law is exact.
- A **between-variance self-consistency correction** rescales the deployed predictive width to the contraction law — posterior contraction made an explicit calibration target (not addressed by SBI's fixed-$n$ coverage or by BOED).

### 2.6 Deployment calibration stack

- **Expanding-window head fine-tune + per-item affine shift**: train once on the wide prior, then re-fit only the head against the application's own SVI reference — corrects the prior$\to$operational mislocation.
- **Target warp** (log$_2$ for folds, logit for rates) crosses endpoint-boundary skew; monotone $\Rightarrow$ KS-invariant; unwarped to natural units for the decision.
- **Diagnostics**: PIT–KS (conditional calibration, an SBC/PIT descendant) and marginal KS certify the deployed predictive, per interim and per item.

### 2.7 A registry keyed by the endpoint functional

- One item-general network **per functional** ($S_1$ rate, $S_2$ fold, $S_3/S_4$ relative-change, $S_5$ standardised between-group difference, $S_6$ seroconversion); the token and warp are structural to the estimand.
- **Not** keyed by item family — a single net transfers across assays/instruments (a psychometric-trained net deploys on titre families) within a $\le 64$-item, $\le 10$-level envelope; per-app calibration is the only gate.
- Between-group differences are **standardised** (Cohen's $d$) so one net is scale-free across applications.
- An application deploys a *manifest* of instances; the per-functional set is marginally sharper than a single conditioned encoder.

## 3. Results

### 3.1 Binomial and Multivariate normal simulation studies

Compare to nested-MC with SVI, importance sampling, regression-based PPS.
Point out we built one amortiser that can be flexibly deployed to arbitrary data structures up to J=64.

### 3.2 HIV vaccine immunicity

This is CAVD DataSpace / HVTN 505.
Cumulative exceedance token.
Vaccine vs placebo arm comparison of endline immunogenicity markers. 
Simple to start with; point out dilutions are ordered categorical, PCM or ordered categorical give excellent fit.

### 3.3 Flu vaccine immunicity

This is SDY312.
Pre/post.
Can inspect response to each flu strain, nicely heterogeneous PPS results, well explainable.
Point out multiple rho, can inspect joint posterior, decision roles across all of them.

### 3.4 Flu head to head vaccine comparison

This is SDY269.
Mix of pre/post and group, with one common item to compare on.

### 3.5 International red cross

Pre/post, large n, many interims.
Stratified psychometric scales, as scores often hard to interpret
Large effect across survey items; point out posterior contraction of our BvM-net, discuss Rao Blackwellization of Monte Carlo error

### 3.6 EU refugees

Pre/post.
Mainly demonstrate that amortiser also correctly reports small effects, and extends to very different reporting structure

### 3.7 Hope Group intervention

Group mean comparison.
20 items, focus on experimental platform, understand which components work well, which dont.

## References

::: {#refs}
:::