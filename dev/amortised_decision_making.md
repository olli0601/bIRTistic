# Amortized decision-focussed learning for real-time interim analyses in crisis settings

## 1. Summary

In contrast to neural posterior estimation, the focus of this work is on amortization for Bayesian decision making. The key point is that decision rules are low-dimensional, so it should be much easier to perform amortised learning by focussing on what is relevant to make real-life decision in real-time.

## 2. Target objectives

Perhaps the simplest approach to this is to consider an intervention that estimates a treatment effect \(p\). The actions are declaring success (\(A_1\)) or failure (\(A_0\)) based on the binary decision rule
\[
d(x) = 1\{P(H_1 \mid x) > \eta\},
\]
where 
\[
H_1 : p > p_0 \Leftrightarrow \rho := p/p_0 > 0 
\]

In practice, the interventions can be Cholera vaccination, our Hope group interventions in Jordan, alternative medical procedures, or similar. As generative models I am interested in various versions of item-response models such as the ordered categorical or partial credit model, but many other generative models are possible. The baseline treatment effect is without loss of generality \(p_0 = 0\), for example \(p\) can be coded as a contrast between the intervention and standard-of-care arm. The decision threshold \(\eta\) can be set to \(0.89\) in defiance to \(0.95\), but can also be derived as the Bayes optimal threshold given false positive and false negative loss functions.

There are two distinct targets.

### Target 1, current decision: “Have we already won?”

We evaluate the posterior success probability of \(H_1\), \( P(H_1 \mid x) \), and evaluate the binary decision rule

\[
d(x) = 1\{P(H_1 \mid x) > \eta\} \in \{ 0, 1 \},
\]

which declares success or failure. This answers the question: Have we already accumulated enough evidence to declare success now?


### Target 2, predictive probability of success (PPS): "How likely will we win in the future?"

We are interested in real-time decision-making, typically referred to as interim analyses to stop for efficacy or to stop for futility. This reframes the problem. 

The new questions are: If the current data \(x\) declare failure, how likely is it that future data collection could overturn this result, or should we stop for failure and safety reasons now? The other way round, if the current data \(x\) declare success, how likely is it that future data collection could overturn this result, or can we stop for efficacy early now? 

Suppose \(m\) additional observations \(z\) will be collected before the final analysis. Let the final decision rule be

\[
d(x,z) = 1\{P(H_1 \mid x,z) > \eta\}.
\]

The predictive probability of success (PPS) is

\[
PPS(x) =
\mathbb{E}_{z \sim p(z\mid x)}
[
1\{P(H_1 \mid x,z) > \eta\}
] \in [ 0, 1],
\]
where the posterior predictive data generating distribution is \[ 
    p(z|x) = \int p(z|\theta)p(\theta|x)d\theta.
\]


The PPS quantifies the probability that the intervention will declare success in the future after \(m\) further hypothetical data points \(z\), given the current data \(x\).

Thus PPS thus answers: How likely are we to win in the future if the intervention continues?

---

# 3. Examples

## 3.1 Binomial Model

Let

\[
X_i \sim \text{Bernoulli}(p)
\]

with prior

\[
p \sim \text{Beta}(a,b).
\]

The Bernoulli distribution belongs to the exponential family:
$$
\begin{aligned}
p(x_i \mid p) & 
= p^{x_i}(1-p)^{1-x_i} = \exp\{x_i \log p + (1-x_i)\log(1-p)\} 
\\
& 
= \exp\{\log\left(\frac{p}{1-p}\right) x_i + \log(1-p)\},
\end{aligned}
$$
with natural parameter \(\eta(p) = \log(p/(1-p))\), sufficient statistics \(T(x_i) = x_i\), and log-partition function \(A(p) = -\log(1-p)\). For \(n\) observations, the likelihood depends only on the sufficient statistic \(T(x_{1:n}) = \sum_{i=1}^n x_i = k^n\), the total number of successes. After observing \(x_{1:n}\) data points with $k^n$ successes the posterior is

\[
p \mid x \sim \text{Beta}(a+k^n, b+n-k^n).
\]


Suppose the future, remaining data \(z\) comprise \(m\) additional data points. Then the posterior predictive distribution of $k^m$ future successes is

$$
\begin{aligned}
k^m & \sim \int \text{Binomial}(m,p) \text{Beta}(p; a+k^n, b+n-k^n) \\
& = \text{Beta-Binomial}(m,a+k^n,b+n-k^n).
\end{aligned}
$$

If \(k^m\) additional successes occur:

\[
p \mid x,z \sim
\text{Beta}(a+k^n+k^m, b+n+m-k^n-k^m).
\]

The final posterior success probability is

\[
P(H_1 \mid x,z)
=
1 - F_{Beta}(p_0; a+k^n+k^m, b+n+m-k^n-k^m).
\]

In the future, with \(z\) additional data, success occurs when

\[
P(H_1 \mid x,z) > \eta.
\]

Since the posterior depends only on the sufficient statistic \(k^m\) (the number of successes in the future data), and \(P(H_1 \mid x, k^m)\) is monotonically increasing in \(k^m\), we can find the threshold

\[
k^m_\star = \min\{k^m \in \{0, 1, \ldots, m\} : P(H_1 \mid x, k^m) > \eta\}.
\]

Then the PPS becomes a simple tail sum:

\[
PPS(x)
=
\sum_{k^m = k^m_\star}^m
p(k^m \mid x),
\]

where \(p(k^m \mid x)\) is the Beta–Binomial predictive distribution. Explicitly,

\[
PPS(x)
=
\sum_{k^m = k^m_\star}^m
{m \choose k^m}
\frac{B(a+k^n+k^m,b+n+m-k^n-k^m)}{B(a+k^n,b+n-k^n)}.
\]

---

## 3.2 Categorical model

Let

\[
x_i \sim \text{Categorical}(p)
\]

with parameter vector \(p = (p_1, \ldots, p_K)\) where \(\sum_{k=1}^K p_k = 1\), and prior

\[
p \sim \text{Dirichlet}(\alpha),
\]

where \(\alpha = (\alpha_1, \ldots, \alpha_K)\).

The Categorical distribution belongs to the exponential family:

$$
\begin{aligned}
p(x_i \mid p) & 
= \prod_{k=1}^K p_k^{1_{x_i = k}} 
= \exp\left\{\sum_{k=1}^K 1_{x_i = k} \log p_k\right\}
\\
&
= \exp\left\{\sum_{k=1}^{K-1} \log\left(\frac{p_k}{p_K}\right) 1_{x_i = k} + \log p_K\right\},
\end{aligned}
$$
with natural parameters \(\eta_k(p) = \log(p_k/p_K)\) for \(k = 1, \ldots, K-1\), sufficient statistics \(T_k(x_i) = 1_{x_i = k}\) for \(k = 1, \ldots, K-1\), and log-partition function \(A(p) = -\log p_K\). For \(n\) observations, the likelihood depends only on the sufficient statistic \(T(x_{1:n}) = (k_1^n, \ldots, k_K^n)\), the count vector where \(k_k^n = \sum_{i=1}^n 1_{x_i = k}\) is the number of observations in category \(k\). After observing \(x_{1:n}\) data points with counts \(k_k^n\) the posterior is

\[
p \mid x \sim \text{Dirichlet}(\alpha_1 + k_1^n, \ldots, \alpha_K + k_K^n).
\]

Suppose the future, remaining data \(z\) comprise \(m\) additional data points. Then the posterior predictive distribution of the future count vector \((k_1^m, \ldots, k_K^m)\) is

$$
\begin{aligned}
(k_1^m, \ldots, k_K^m) & \sim \int \text{Multinomial}(m, p) \, \text{Dirichlet}(p; \alpha_1 + k_1^n, \ldots, \alpha_K + k_K^n) \, dp \\
& = \text{Dirichlet-Multinomial}(m, \alpha_1 + k_1^n, \ldots, \alpha_K + k_K^n).
\end{aligned}
$$

If future counts \((k_1^m, \ldots, k_K^m)\) occur:

\[
p \mid x,z \sim \text{Dirichlet}(\alpha_1 + k_1^n + k_1^m, \ldots, \alpha_K + k_K^n + k_K^m).
\]

The final posterior success probability is

\[
P(H_1 \mid x,z) = \int 1_{p \in H_1} \, \text{Dirichlet}(p; \alpha_1 + k_1^n + k_1^m, \ldots, \alpha_K + k_K^n + k_K^m) \, dp,
\]

which can be computed via numerical integration or Monte Carlo sampling from the posterior.

In the future, with \(z\) additional data, success occurs when

\[
P(H_1 \mid x,z) > \eta.
\]

Since the posterior depends only on the sufficient statistic \((k_1^m, \ldots, k_K^m)\) (the count vector of future data), we can define the success set directly in terms of this count vector. Let

\[
S = \{(k_1^m, \ldots, k_K^m) : \sum_{j=1}^K k_j^m = m, \, P(H_1 \mid x, k_1^m, \ldots, k_K^m) > \eta\},
\]

where the constraint \(\sum_{j=1}^K k_j^m = m\) ensures the counts sum to the total number of future observations. Then the PPS becomes

\[
PPS(x) = \sum_{(k_1^m, \ldots, k_K^m) \in S} p(k_1^m, \ldots, k_K^m \mid x),
\]

where \(p(k_1^m, \ldots, k_K^m \mid x)\) is the Dirichlet-Multinomial predictive distribution:

\[
p(k_1^m, \ldots, k_K^m \mid x) = \frac{m!}{k_1^m! \cdots k_K^m!} \frac{B(\alpha_1 + k_1^n + k_1^m, \ldots, \alpha_K + k_K^n + k_K^m)}{B(\alpha_1 + k_1^n, \ldots, \alpha_K + k_K^n)},
\]

where \(B(\cdot)\) is the multivariate Beta function. This is analytically computable but combinatorially more expensive than the Binomial case due to the larger state space of count vectors.

## 3.3 Multivariate Normal model with known correlation structure

A canonical high-dimensional analytically tractable PPS benchmark for continuous outcomes is the multivariate normal model with structured covariance. Let

\[
y_n \mid \mu, \sigma^2 \;\sim\; \mathrm{MVN}(\mu, \sigma^2 K), \qquad n = 1, \ldots, N,
\]

where \(K = R R^\top \in \mathbb{R}^{J \times J}\) is a known positive-definite covariance shape (for example, \(R\) a Cholesky factor of an AR(1) or arbitrary correlation matrix), and the unknowns are the mean vector \(\mu \in \mathbb{R}^J\) and the scalar variance \(\sigma^2 > 0\).

We test \(J\) component-wise alternative hypotheses,

\[
H_{1j} : \mu_j > \mu^0_j \quad \text{for all } j = 1, \ldots, J,
\]

with the baseline pinned at \(\mu^0_j := 1\) for concreteness. Using the relative-effect parameterisation,

\[
\rho_j \;:=\; 1 - \mu_j/\mu^0_j \;=\; 1 - \mu_j,
\]

the alternative is equivalently \(H_{1j} \Leftrightarrow \rho_j < 0\).

**Conjugate prior** (Normal–Inverse-Gamma with metric \(K\)):

\[
\mu \mid \sigma^2 \;\sim\; \mathrm{MVN}(\mu_0,\, \sigma^2 \Lambda_0^{-1}), \qquad \sigma^2 \;\sim\; \mathrm{Inv\text{-}Gamma}(a_0, b_0),
\]

with hyperparameters \((\mu_0, \Lambda_0, a_0, b_0)\) known.

**Exponential family and sufficient statistics.** With \(K\) known, the MVN likelihood is in the exponential family with natural parameters \(\eta(\mu, \sigma^2) = \big( K^{-1}\mu/\sigma^2,\, -1/(2\sigma^2) \big)\) and sufficient statistics for \(n\) current observations

\[
T_1(x_{1:n}) \;=\; \sum_{i=1}^n x_i, \qquad T_2(x_{1:n}) \;=\; \sum_{i=1}^n x_i^\top K^{-1} x_i,
\]

so the likelihood factors as \(p(x_{1:n} \mid \mu, \sigma^2) = h(x_{1:n})\, \exp\{\eta^\top T - n A(\mu, \sigma^2)\}\) with \(A(\mu, \sigma^2) = \tfrac{J}{2}\log\sigma^2 + \tfrac{1}{2\sigma^2}\mu^\top K^{-1} \mu\).

**Posterior** after current data \(x_{1:n}\), writing \(\bar{x}_n = T_1(x_{1:n})/n\):

$$
\begin{aligned}
\Lambda_n &= \Lambda_0 + n K^{-1},\\
\mu_n &= \Lambda_n^{-1}\!\big( \Lambda_0 \mu_0 + n K^{-1} \bar{x}_n \big),\\
a_n &= a_0 + nJ/2,\\
b_n &= b_0 + \tfrac{1}{2}\!\left( T_2(x_{1:n}) + \mu_0^\top \Lambda_0 \mu_0 - \mu_n^\top \Lambda_n \mu_n \right).
\end{aligned}
$$

Marginalising \(\sigma^2\) gives the multivariate Student-\(t\) marginal posterior of the mean,

\[
\mu \mid x \;\sim\; t_{2a_n}\!\big( \mu_n,\; (b_n/a_n)\,\Lambda_n^{-1} \big),
\]

with per-component marginals

\[
\mu_j \mid x \;\sim\; t_{2a_n}\!\big( \mu_{n,j},\; (b_n/a_n)\,[\Lambda_n^{-1}]_{jj} \big).
\]

The current posterior probability of the per-component hypothesis is therefore closed form,

\[
P(H_{1j} \mid x) \;=\; 1 - F_{t_{2a_n}}\!\left( \frac{1 - \mu_{n,j}}{\sqrt{(b_n/a_n)\,[\Lambda_n^{-1}]_{jj}}} \right),
\]

and the joint hypothesis \(P(H_1 \mid x) = P(\mu_j > 1 \;\forall j \mid x)\) is a multivariate-\(t\) orthant probability, evaluated by quasi-Monte Carlo (Genz QMC).

**Posterior predictive of future data.** Suppose \(m\) additional observations \(z_{1:m}\) will be collected. Conditional on \(x\), the future data are jointly matrix-\(t\) distributed, and their joint distribution factors through the future sufficient statistics

\[
T_1(z_{1:m}) \;=\; \sum_{i=1}^m z_i, \qquad T_2(z_{1:m}) \;=\; \sum_{i=1}^m z_i^\top K^{-1} z_i,
\]

with \(\bar{z}_m := T_1(z_{1:m})/m\). The marginal predictive of \(\bar{z}_m\) integrating out \((\mu, \sigma^2)\) is multivariate \(t\),

\[
\bar{z}_m \mid x \;\sim\; t_{2a_n}\!\left( \mu_n,\; (b_n/a_n)\big( \tfrac{1}{m} K + \Lambda_n^{-1} \big) \right),
\]

and \(T_2(z_{1:m})\) conditional on \(\bar{z}_m\) is a quadratic form whose distribution is a scaled \(F\) (the future analogue of the residual sum of squares under NIG).

**Posterior after current and future data.** With \(\bar{x}_{n+m} = (n\bar{x}_n + m \bar{z}_m)/(n+m)\):

$$
\begin{aligned}
\Lambda_{n+m} &= \Lambda_0 + (n+m) K^{-1}, & a_{n+m} &= a_0 + (n+m) J / 2,\\
\mu_{n+m} &= \Lambda_{n+m}^{-1}\!\big( \Lambda_0 \mu_0 + n K^{-1} \bar{x}_n + m K^{-1} \bar{z}_m \big),\\
b_{n+m} &= b_0 + \tfrac{1}{2}\!\left( T_2(x_{1:n}) + T_2(z_{1:m}) + \mu_0^\top \Lambda_0 \mu_0 - \mu_{n+m}^\top \Lambda_{n+m} \mu_{n+m} \right).
\end{aligned}
$$

Crucially, \(\Lambda_{n+m}\) and \(a_{n+m}\) are deterministic given the future cohort size \(m\); only \(\mu_{n+m}\) and \(b_{n+m}\) depend on the random future statistics \((\bar{z}_m, T_2(z_{1:m}))\).

The final per-component success probability is

\[
P(H_{1j} \mid x, z) \;=\; 1 - F_{t_{2a_{n+m}}}\!\left( \frac{1 - \mu_{n+m, j}}{\sqrt{(b_{n+m}/a_{n+m})\,[\Lambda_{n+m}^{-1}]_{jj}}} \right).
\]

**Predictive probability of success.** The per-component decision rule \(d_j(x, z) = 1\{P(H_{1j} \mid x, z) > \eta\}\) depends on \(z\) only through the two sufficient statistics \((\bar{z}_m, T_2(z_{1:m}))\), so the PPS reduces to an integral against the closed-form predictive of these statistics. Define the critical region in sufficient-statistic space,

\[
A_j(x) \;=\; \left\{ (\bar{z}, S) \;:\; \frac{1 - \mu_{n+m,j}(\bar{z})}{\sqrt{(b_{n+m}(\bar{z}, S)/a_{n+m})\,[\Lambda_{n+m}^{-1}]_{jj}}} \;<\; F_{t_{2a_{n+m}}}^{-1}(1 - \eta) \right\}.
\]

Then

\[
PPS_j(x) \;=\; \int 1_{A_j(x)}(\bar{z}, S)\; p(\bar{z}, S \mid x)\; d\bar{z}\, dS,
\]

a two-dimensional quadrature against the marginal predictive of \((\bar{z}_m, T_2(z_{1:m}))\) derived above. The joint-hypothesis PPS replaces the per-component CDF by the multivariate-\(t\) orthant probability \(P(\mu > \mathbf{1} \mid x, z)\) and integrates against the same predictive.

**Special case: \(\sigma^2\) known.** Fixing \(\sigma^2\), the conjugate prior collapses to MVN on \(\mu\) alone, posteriors and predictives are Gaussian, and the per-component PPS reduces to a Gaussian tail probability,

\[
PPS_j(x) \;=\; \Phi\!\left( \frac{\mu_{n,j} - 1 - z_\eta \sqrt{\sigma^2\,[\Lambda_{n+m}^{-1}]_{jj}}}{\sqrt{\sigma^2\,\big([\Lambda_n^{-1}]_{jj} - [\Lambda_{n+m}^{-1}]_{jj} + [K]_{jj}/m\big)}} \right),
\]

where \(z_\eta = \Phi^{-1}(\eta)\), i.e. a single \(\Phi\) evaluation per interim per component. The full NIG case differs only by replacing the inner Gaussian tails with Student-\(t\) tails and adding one \(F\)-distribution integration over the residual sum of squares.

**Why this benchmark.** With \(K\) arbitrary known PSD and \(J\) arbitrary (here we target \(J\) up to several hundreds or thousands), the construction supplies a high-dimensional continuous-outcome decision problem whose PPS is closed-form, with which the nested-MC / IS / SMC / regression estimators of Section 6 can be validated.

### 3.3.1 Concrete simulation setup

We fix \(\sigma^2 = 1\) (units chosen so that the per-component noise has unit variance), use the Gaussian special case throughout, and run the benchmark at three problem dimensions

\[
J \in \{50, 100, 200\},
\]

with a total cohort \(N = 500\) accrued over 10 monthly interims of 50 units each, so the future cohort size at interim \(t\) is \(m_t = N - n_t\) with \(n_t = 50 t\). The decision threshold is \(\eta = 0.89\) and the baseline is \(\mu^0 = \mathbf{1}_J\).

**Prior.** A weakly informative \(g\)-prior in the same metric \(K\) keeps the algebra closed,

\[
\mu \;\sim\; \mathrm{MVN}\!\big(\mathbf{1}_J,\; \tau_0^2\, K\big), \qquad \tau_0^2 = 100,
\]

so \(\Lambda_0^{-1} = \tau_0^2 K\) and \(\Lambda_0 = K^{-1}/\tau_0^2\). With \(\sigma^2 = 1\),

\[
\Lambda_n \;=\; (\tau_0^{-2} + n)\,K^{-1}, \qquad \Lambda_n^{-1} \;=\; \frac{K}{\tau_0^{-2} + n},
\]

so the per-component closed-form PPS specialises to a clean \(\Phi\)-tail with \(\Lambda_\bullet^{-1}\) replaced by a scaled \(K\).

**True \(\mu\).** Half the components above baseline, half below, to stress-test the per-component decision:

\[
\mu_{\text{true}, j} \;=\; \begin{cases} 1 + \Delta & j \le J/2 \\ 1 - \Delta & j > J/2 \end{cases}, \qquad \Delta = 0.3.
\]

This sets the true component-wise rejection rate to exactly \(1/2\), so the benchmark probes both efficacy and futility regimes simultaneously across the \(J\) hypotheses.

**Choice of \(R\).** We need a nontrivial known correlation structure: enough off-diagonal mass that components are correlated (so the closed-form predictive of \((\bar{z}_m, T_2(z))\) is non-trivial), but well-conditioned at all three \(J\). We propose three structures, used as a panel:

1. **AR(1) Cholesky** (primary, ordered components):

   \[
   K^{(1)}_{ij} \;=\; \rho^{|i - j|}, \qquad \rho = 0.7.
   \]

   Cholesky factor \(R^{(1)}\) with \(R^{(1)} (R^{(1)})^\top = K^{(1)}\). Banded, PSD for any \(\rho \in (-1, 1)\), condition number bounded uniformly in \(J\) (\(\kappa(K^{(1)}) = (1+\rho)/(1-\rho) \approx 5.67\) for \(\rho = 0.7\)). Reflects a natural ordering (e.g.\ time or position along a scale).

2. **Block equicorrelation** (mid difficulty, exchangeable within block):

   \[
   K^{(2)} \;=\; \begin{bmatrix} K_w & \rho_b \mathbf{1}\mathbf{1}^\top & \cdots \\ \rho_b \mathbf{1}\mathbf{1}^\top & K_w & \cdots \\ \vdots & & \ddots \end{bmatrix},
   \]

   with \(B = J/10\) blocks of size \(10\), intra-block correlation \(\rho_w = 0.8\) (so \(K_w = (1-\rho_w) I + \rho_w \mathbf{1}\mathbf{1}^\top\)), inter-block correlation \(\rho_b = 0.1\). PSD by construction (\(\rho_w > \rho_b \ge 0\)). Captures clustered (e.g.\ subscale) dependence.

3. **Low-rank-plus-diagonal factor structure** (hard, dense long-range correlation):

   \[
   K^{(3)} \;=\; \beta \beta^\top + \psi I_J, \qquad \beta \in \mathbb{R}^{J \times r},\; r = 5,\; \psi = 0.1,
   \]

   with \(\beta_{jk} \sim \mathcal{N}(0, 1/r)\) drawn once and fixed, then \(K^{(3)}\) rescaled so \(\mathrm{diag}(K^{(3)}) = \mathbf{1}\). \(R^{(3)}\) is the Cholesky factor. Five-factor latent structure, dense correlations decaying smoothly with no banding.

All three are scaled to unit diagonal so the per-component noise is comparable across \(J\) and across structures. The Cholesky factors \(R^{(\cdot)}\) are precomputed once per \((J, R\text{-type})\) cell and shared across the 10 interims.

**What we compare.**

| Cell            | Knobs                                           | Notes                                   |
| --------------- | ----------------------------------------------- | --------------------------------------- |
| Cell A (easy)   | \(R = R^{(1)}\) AR(1), \(J \in \{50,100,200\}\) | Banded, well-conditioned.               |
| Cell B (medium) | \(R = R^{(2)}\) block, same \(J\)               | Block structure; eigenvalues clustered. |
| Cell C (hard)   | \(R = R^{(3)}\) factor, same \(J\)              | Dense long-range; eigenvalue decay.     |

For each cell and each interim \(t \in \{1, \ldots, 10\}\), we compare the closed-form per-component \(PPS_j(x)\) (from the \(\Phi\)-tail above) against the four estimators of Section 6 (nested-MC, self-normalised IS, moment-matching IS, SMC resample-move, regression on \(w(z)\)). Aggregate diagnostics: per-component absolute error and bias against the closed form; per-cell timing; ESS / \(\hat{k}\) for the IS variants; tempering steps \(T\) for SMC.

---

## 3.4 Item-response model

In our applications, we are interested in Item Response Theory (IRT) models that can be fitted to real-time survey data collected about interventions in humanitarian and social science settings. 
The data consist of responses to Likert-scale survey items collected from participants at baseline and endline. The data increases in size as more individuals are surveyed. The goal is to quantify intervention effectiveness based on the current data, and to quantify the predictive probability of success (PPS).

Throughout, let $i$ index participants ($i = 1, \ldots, n$), $j$ index items/questions ($j = 1, \ldots, J$), and $k$ index response categories ($k = 1, \ldots, K$). The response of participant $i$ to item $j$ at time $t$ (baesline or endline) is $Y_{ijt}$; for simplicity we suppress $t$ in what follows. A widely-used IRT model which we use here is the partial credit model. The PCM models ordered categorical probabilities with cumulative logits,

$$
\begin{aligned}
P(Y_{ij} = k) &= \text{softmax}(\phi_{i,j,k}) = \frac{\exp(\phi_{i,j,k})}{\sum_{k'=1}^K \exp(\phi_{i,j,k'})}\\
\phi_{i,j,1} &= 0\\
\phi_{i,j,k} &= \sum_{s=1}^k \lambda_{j} \cdot \left(\theta_i + \mathbf{X}_{i,j}^T \boldsymbol{\beta} - c_{j,s}\right), \quad k=2,\dotsc,K
\end{aligned}
$$

The number of free parameters are \(N\) latent skills parameters \(\theta_i\), one for each participant; \(J(K-1)\) incremental skill thresholds \(c_{j,s}\), one for each categorical increment for each response item; \(J\) item loadings; and \(P\) fixed participants effects. The total number of parameters is thus \(N+JK + P\), comprising \(N\) local parameters that grow wich each new participant and \(JK + P\) global parameters that are shared across participants.

By construction, the PCM admits an incremental log-risk structure:

$$
\begin{aligned}
\log \frac{\Pr(Y_{i,j} = k)}{\Pr(Y_{i,j} = k-1)}
=  \lambda_{j} \cdot \left(\theta_i + \mathbf{X}_{i,j}^T \boldsymbol{\beta} - c_{j,k}\right)
\end{aligned}
$$

The model thus does not follow the proportional cumulative odds assumptions, but it has a proportional incremental risk assumption: when $\beta$ represents an intervention effect (e.g., $X_{i,j,t} = 1$ at time $t=1$ and $X_{i,j,t} = 0$ at time $t=0$), then 

$$
\begin{aligned}
&\frac{\Pr(Y_{i,j,t=1} = k \mid \eta_{i,j,t=1})}{\Pr(Y_{i,j,t=1} = k-1 \mid \eta_{i,j,t=1})} \bigg/
\frac{\Pr(Y_{i,j,t=0} = k \mid \eta_{i,j,t=0})}{\Pr(Y_{i,j,t=0} = k-1 \mid \eta_{i,j,t=0})}\\
&\quad = 
\exp\Big( \lambda_j ( \theta_i + \beta ) - \lambda_j \theta_i  \Big)
= 
\exp( \lambda_j \beta ),
\end{aligned}
$$

so when measured in incremental risks, the effect of the intervention is proportional to $\exp\beta$, and the same regardless of the category $k$.

Computationally, the PCM involves direct evaluations of category specific events ($Pr(Y=k)$ not $Pr(Y\leq k)$) and therefore tends to be much faster to evaluate and less prone to numerical issues than other IRT models. Standard priors can be attached to all free parameters.

---

# 4. No reduction of the PPS to today's knowledge

A particular confusion that often arises with the PPS is that it does not fall back to today's knowledge. 

A key concept in Bayesian forecasting is self-consistency, in that if we make posterior predictions \(z\) from today's posterior and then consider the updated posterior \(p(\theta|x,z)\), the two expansion and contraction steps cancel each other out. 

In our setting, this Bayesian prediction self-consistency property equates to

\[
\mathbb{E}_{z|x}[P(H_1|x,z)] = P(H_1|x),
\]

and it follows from replacing \( p(\theta|x,z)p(z|x) \) with \( p(\theta,z|x) \), exchanging integrals, and integrating out \(z\).

The key point is that the PPS contains the indicator decision rule

\[
1\{P(H_1 \mid x,z) > \eta\},
\]

and so

\[
\mathbb{E}_{z|x}[1\{P(H_1|x,z) > \eta\}] 
\neq 
1\{P(H_1|x) > \eta\},
\]

since \( E(f(X)) \neq f(E(X)) \). Conceptually the PPS is a probability in \( [0,1]\) whereas today's decision rule always evaluates to a binary value.

---

# 5. Links to general decision-theoretic learning

## 5.1 Bayes optimal decisions under loss

In the broader decision theoretic literature, the objective is to evaluate the expected loss (or expected utility) of all possible actions given the observed data, and then choose the action that minimizes the expected loss.

Let \(\mathcal{A}\) denote the action space. Given data \(x\), the Bayes optimal action is

\[
a^* = \arg\min_{a \in \mathcal{A}} \mathbb{E}_{\theta \mid x}[L(a, \theta)]
\]

where \(L(a, \theta)\) is the loss incurred by taking action \(a\) when the unknown parameters are \(\theta\).

In our setting, for final data \(x\), the decision theoretic problem simplifies to two actions:

* \(a_1\): declare success, adopt intervention
* \(a_0\): declare failure, reject intervention

and two partitions on the unknown true treatment effect:

* \(H_1\): success hypothesis (\(p > p_0\))
* \(H_0\): null hypothesis (\(p \leq p_0\))

We can define the loss matrix:

| Action                    | \(H_1\) true | \(H_0\) true |
| ------------------------- | ------------ | ------------ |
| Declare success (\(a_1\)) | 0            | \(L_{FP}\)   |
| Declare failure (\(a_0\)) | \(L_{FN}\)   | 0            |

where \(L_{FP}\) is the cost of a false positive decision and \(L_{FN}\) is the cost of a false negative decision.

Given observed data \(x\), the expected loss of declaring success is

\[
\mathbb{E}_{H \mid x}[L(a_1, H)]
= L_{FP} \cdot P(H_0 \mid x) + 0 \cdot P(H_1 \mid x)
= L_{FP} \cdot P(H_0 \mid x),
\]

amd similarly the expected loss of declaring failure is

\[
\mathbb{E}_{H \mid x}[L(a_0, H)]
= L_{FN} \cdot P(H_1 \mid x).
\]

The Bayes optimal action minimizes expected loss,

\[
a^*(x) = \arg\min_{a \in \{a_0, a_1\}} \mathbb{E}_{H \mid x}[L(a_1, H)].
\]

And so, to declare success we want 
$$
\begin{aligned}
& 
\mathbb{E}_{H \mid x}[L(a_1, H)] < \mathbb{E}_{H \mid x}[L(a_0, H)] \\
\Leftrightarrow &
L_{FP} \cdot (1 - P(H_1 \mid x)) < L_{FN} \cdot P(H_1 \mid x) \\
\Leftrightarrow &
P(H_1 \mid x) > L_{FP} / ( L_{FP} + L_{FN} ).
\end{aligned}
$$


This shows that the Bayes optimal decision threshold \(\eta\) can be expressed in terms of utilities or loss terms. In particular, if false positives are more costly, \(\eta > 0.5\) and if false negatives are more costly, \(\eta < 0.5\). This provides a utility-based foundation for choosing \(\eta\) in contrast to ad-hoc 0.89 or 0.95.

## 5.2 Net benefit 

A related concept is expected utility, normalised such that 1 unit corresponds to 1 true positive. Under the above action/partition set matrix, the relative cost of one false positive is \(w\) and true pos and true neg have utility/cost 0. This does not depend on the particular units of the losses, and so is more easily interpretable.

Under these utilities, the net benefit of declaring success is given by

$$
\begin{aligned}
NB(x) 
& 
= P(H_1 \mid x) - w \cdot P(H_0 \mid x) 
\\
&
= P(H_1 \mid x) - w \cdot (1 - P(H_1 \mid x) )
\end{aligned}
$$

and the net benefit decision rule is

\[
1\{ NB(x) > 0 \}.    
\]

Under Bayes optimality, we find that the relative cost of a false positive must be \( w = L_{FP} / L_{FN} \). Using the expression for \(\eta\) above, we then find \( w = \eta / (1 - \eta)\). Thus we can express the net benefit in the standard form

\[
NB(x) = P(H_1 \mid x) - \frac{\eta} {1 - \eta} \cdot (1 - P(H_1 \mid x) ).
\]

This shows that our learning problem directly connects to standard cost-benefit analyses. Our PPS is equivalent to posterior predictive net benefit over future data \(z\), and so finding a way to amortise the PPS also means a solution to amortising posterior predictive net benefit.

---

# 6. State-of-the-art approaches to estimate PPS

State-of-the-art approaches consider that we have observed a specific interim dataset \(x_{1:n}\) and have access to the posterior distribution \(p(\theta \mid x)\). The core task is to estimate

\[
   y^{(s)} := P(H_1 \mid x, z^{(s)}) = \int 1_{\theta \in H_1} \, p(\theta \mid x, z^{(s)}) \, d\theta.
\]
   
for fixed \(x\) and for each simulated future data \( z^{(s)} \sim p(\cdot \mid \theta^{(s)}) \) where \( \theta^{(s)} \sim p(\theta \mid x) \). From this, the PPS is straightforward to compute via

\[
PPS(x) = \sum_{s=1}^S 1\{ y^{(s)} > \eta\}.    
\]

### 6.1 Nested Monte Carlo

A computationally costly, but robust approach is to numerically estimate the new joint posterior \( p(\theta \mid x, z^{(s)}) \) and then simply obtain the label \( y^{(s)} \) by evaluating the above integral over posterior draws from \( p(\theta \mid x, z^{(s)}) \).

### 6.2 Importance sampling (self-normalized)

Since \( p(\theta \mid x, z^{(s)}) \propto p(\theta \mid x)\, p(z^{(s)} \mid \theta) \), we can re-use all existing posterior draws \( \theta_{k} \sim p(\theta \mid x) \) for \(k=1,\dotsc,K\) and reweight these by the future-data likelihood. 

Working in log space for stability, for each \(z^{(s)}\) separately, we compute for all \(k=1,\dotsc,K\) the importance sampling weights
$$
\begin{aligned}
\log w_k^{(s)} &= \log p(z^{(s)} \mid \theta_k) = \sum_{i=1}^m \log p\big(z^{(s)}_i \mid \theta_k\big), \\
\tilde w_k^{(s)} &= \operatorname{softmax}_k\!\big(\log w_k^{(s)}\big) = \frac{\exp\!\big(\log w_k^{(s)} - \ell^{(s)}\big)}{\sum_j \exp\!\big(\log w_j^{(s)} - \ell^{(s)}\big)}, \quad \ell^{(s)} = \log\!\textstyle\sum_j \exp \log w_j^{(s)}, \\
y^{(s)} &\approx \sum_{k=1}^K \tilde w_k^{(s)} \, 1_{\theta_k \in H_1}.
\end{aligned}
$$
The log-sum-exp / softmax map is the numerically stable form of \( w / \sum w \): subtracting the maximum log-weight before exponentiating prevents overflow while leaving the normalized weights unchanged.

Reliability is monitored by the effective sample size [@kong1992note; @liu2001monte]
\[
\mathrm{ESS} = \frac{\big(\sum_k w_k\big)^2}{\sum_k w_k^2} = \frac{1}{\sum_k (\tilde w_k)^2}, \qquad \frac{\mathrm{ESS}}{K} \in [1/K,\, 1],
\]
equivalently the second-order weight moment \( \mathbb{E}(\tilde w^2) = \tfrac{1}{K}\sum_k \tilde w_k^2 \), with \( \mathrm{ESS}/K = 1/\big(K^2\, \mathbb{E}(\tilde w^2)\big) \). 

Self-normalized IS is consistent but \(O(1/K)\) biased, and its variance is finite only when \( \mathbb{E}_{p(\theta \mid x)}\!\big[p(z \mid \theta)^2\big] < \infty \); the Pareto-smoothed importance sampling (PSIS) tail index \( \hat k \) [@vehtari2024pareto] both estimates this and stabilizes the largest weights, with \( \hat k > 0.7 \) flagging an unreliable estimate.

The main issue is that the proposal/target mismatch grows with the amount of assimilated future data: \( \mathrm{KL}\big(p(\theta \mid x, z)\,\|\,p(\theta \mid x)\big) \) increases in \(m\), so the weight mass concentrates on a single draw and \( \mathrm{ESS} \to 1 \). Empirically, at the earliest interim of our case study (\( n \approx 48 \) current vs \( m \approx 455 \) future units) the weights collapse to \( \mathrm{ESS} \approx 1 \) after even a *single* future participant, with PSIS \( \hat k = \infty \). Plain IS labels are therefore trustworthy only when \(z\) is small relative to \(x\) (late interims). The two corrections below target this regime.

### 6.3 Moment-matching importance sampling

Moment-matching IS [@paananen2021implicitly] repairs a mild proposal/target mismatch without new model fits, by transforming the draws so the transformed cloud better covers the target and reweighting with the change-of-variables Jacobian. Starting from the IS weights \( \tilde w_k \) of 6.1, compute the weighted and proposal moments
\[
\hat\mu_w = \sum_k \tilde w_k\, \theta_k, \qquad \hat\mu_q = \tfrac{1}{K}\sum_k \theta_k, \qquad (\text{optionally } \hat\Sigma_w,\ \hat\Sigma_q),
\]
and apply an invertible affine map \(T\) that matches them. The mean-match step uses
\[
T(\theta) = \theta + (\hat\mu_w - \hat\mu_q), \qquad \theta_k^* = T(\theta_k),
\]
(the covariance-match variant uses \( T(\theta) = \hat\mu_w + L_w L_q^{-1}(\theta - \hat\mu_q) \) with \( \hat\Sigma_\bullet = L_\bullet L_\bullet^\top \)). The transformed draws are reweighted against the target with the Jacobian of \(T^{-1}\),
\[
w_k^* = \frac{p(\theta_k^* \mid x, z^{(s)})}{q^*(\theta_k^*)}, \qquad q^*(\theta^*) = q\big(T^{-1}\theta^*\big)\,\big|\det \nabla T^{-1}\big|,
\]
where \(q\) is the proposal density \(p(\theta \mid x)\) (in practice a diagonal-Gaussian fit to the base draws in an unconstrained reparameterisation, with positive parameters mapped through \(\log\)). One iterates over a small family of transforms and keeps the one maximizing \(\mathrm{ESS}\) (or minimizing \( \hat k \)).

We found that when the base \( \mathrm{ESS} \approx 1 \), the weighted mean \( \hat\mu_w \) *equals* the single dominating draw, so the affine shift only relocates the whole cloud onto that point and \(\mathrm{ESS}\) does not recover. Moment matching corrects mild mismatch but cannot manufacture the support the fixed base draws lack — it never moves a particle to a region the proposal failed to sample. In our case study it leaves the early-interim \( \mathrm{ESS}/K \) unchanged at \( \approx 1/K \).

### 6.4 Sequential Monte Carlo with resample-move

To cross an arbitrarily large \( x \to (x, z) \) gap, another idea is to bridge the proposal to the target through a tempered sequence (annealed importance sampling [@neal2001annealed]; SMC samplers [@delmoral2006smc; @chopin2002sequential]),
\[
\pi_t(\theta) \;\propto\; p(\theta \mid x)\; p(z^{(s)} \mid \theta)^{\beta_t}, \qquad 0 = \beta_0 < \beta_1 < \dots < \beta_T = 1,
\]
so \( \pi_0 = p(\theta \mid x) \) and \( \pi_T = p(\theta \mid x, z^{(s)}) \) (the target). Initialise particles \( \theta_k \sim p(\theta \mid x) \) with uniform weights; at step \(t\):

1. **Reweight** by the incremental likelihood,
   \[
   \tilde w_k = \operatorname{softmax}_k\!\Big( (\beta_t - \beta_{t-1})\, \log p(z^{(s)} \mid \theta_k) \Big).
   \]
2. **Adapt** \( \Delta\beta = \beta_t - \beta_{t-1} \) by bisection so the tempering \( \mathrm{ESS}/K \) hits a target (e.g. \( \tfrac12 \)), which is an automatic schedule [@jasra2011inference; @zhou2016toward].
3. **Resample** the particles by \( \tilde w_k \) (systematic resampling) and reset weights to \(1/K\).
4. **Move** each particle with an MCMC kernel \(M_t\) leaving \(\pi_t\) invariant (resample-move [@gilks2001following]). We use Metropolis-adjusted Langevin (MALA [@roberts1996exponential]) in the unconstrained reparameterisation, adapting the step size toward the optimal acceptance \( \approx 0.574 \) [@roberts1998optimal].

Because the temperature enters only as an exponent, the move kernel's log-density \( \log p(\theta \mid x) + \beta_t \log p(z^{(s)} \mid \theta) \) only need to be compiled **once** with \( \beta_t \) a traced argument and reused across all temperatures. The final particles approximate \( \pi_T \) with uniform weights, giving the label
\[
y^{(s)} \approx \frac{1}{K} \sum_{k=1}^K 1_{\theta_k^{(T)} \in H_1}.
\]
Unlike IS and moment matching, the move step relocates particles into the target's typical set, so SMC crosses an arbitrarily large gap; the cost is \(T\) tempering steps, each a short MCMC sweep. 

We found that at the worst (earliest) interim the adaptive schedule reaches \( \beta = 1 \) in \( \approx 40 \) steps and restores \( \mathrm{ESS}/K \approx 1 \) at a wall-clock cost dominated by the per-step move rather than the one-off compile. The overall computational cost was 7-8 times larger than SVI estimation of the posterior \(p(\theta | x, z^{(s)}) \).

### 6.5 Regression-based functional inference

The estimators 6.1–6.4 all approximate the conditional posterior \(p(\theta \mid x, z^{(s)})\) once per future sample \(z^{(s)}\). The regression-based approach [@strong2014estimating] instead learns a function \(q_\psi\) with tuning parameters \(\psi\) that produces the label as a function of a low-dimensional summary of the future data directly.

Pick a summary statistic \(w:\mathcal{Z}\to\mathbb{R}^d\) of the future data and proceed as follows.

1. **Joint sampling.** For \(s=1,\dotsc,S\), draw
   $$   
   \theta^{(s)} \sim p(\theta\mid x),\qquad
   z^{(s)} \sim p(z\mid \theta^{(s)}),
   $$
   compute based on \(x\) and not \(z^{(s)}\) 
   $$
   y^{(s)} := 1_{\theta^{(s)}\in H_1},
   $$
   and also based \(z^{(s)}\) and not \(x\)
   $$
   w^{(s)} := w(z^{(s)}).
   $$
2. **Fit a regressor** \(q_\psi:\mathbb{R}^d\to[0,1]\) to the pairs \(\{(w^{(s)},y^{(s)})\}_{s=1}^S\) by minimizing the empirical binomial cross-entropy (logit link)
   \[
   \hat\psi=\arg\min_\psi -\frac{1}{S}\sum_{s=1}^S\Big[y^{(s)}\log q_\psi(w^{(s)})+(1-w^{(s)})\log(1-q_\psi(w^{(s)}))\Big],
   \]
   using a a GAM for \(d\le 6\) or a GP for higher \(d\).
3. **Predict** the label \(y^\star\) for any new future sample \(z^\star\) with the learned regressor by \(\hat y(z)=q_{\hat\psi}\big(w(z^\star)\big)\).

Instead of the binary labels \(y^{(s)}\), it is advantageous to consider the continuous \(\rho^{(s)}\) that underlie the alternative hypothesis, \(H_1 : \rho^{(s)} > 0\). This leads to the steps 

1. **Joint sampling.** For \(s=1,\dotsc,S\), draw
   $$   
   \theta^{(s)} \sim p(\theta\mid x),\qquad
   z^{(s)} \sim p(z\mid \theta^{(s)}),
   $$
   compute based on \(x\) and not \(z^{(s)}\) 
   $$
   \rho^{(s)} := \theta^{(s)}/\theta_0
   $$
   and also based \(z^{(s)}\) and not \(x\)
   $$
   w^{(s)} := w(z^{(s)}).
   $$
2. **Fit a regressor** \(q_\psi:\mathbb{R}^d\to\mathbb{R}\) to the pairs \(\{(w^{(s)},\rho^{(s)})\}_{s=1}^S\) by minimizing the least-squares loss
   \[
   \hat\psi=\arg\min_\psi \frac{1}{S}\sum_{s=1}^S\big(\rho^{(s)} - q_\psi(w^{(s)})\big)^2,
   \]
   using a a GAM for \(d\le 6\) or a GP for higher \(d\).
3. **Predict** the label \(\rho^\star\) for any new future sample \(z^\star\) with the learned regressor by \(\hat \rho(z)=q_{\hat\psi}\big(w(z^\star)\big)\), then compute \(y^\star := 1\{ \rho^\star > 0 \}\).

This scheme is guaranteed to create unbiased labels when \(w\) is a sufficient statistic for \(\theta\), as \(S \to \infty \). Indeed, in that case there are \(h\) and \(g\) such that
\[
p(z\mid \theta)=h(z)\,g\big(w(z);\theta\big),
\]
and further,
$$
\begin{aligned} 
p(\theta\mid x, z) & = 
\frac{p(\theta, z \mid x)}{p(z\mid x)} = 
\frac{p(\theta\mid x) p(z\mid\theta)}{p(z\mid x)} \\
& \propto p(\theta\mid x) p(z \mid\theta) \\
& \propto p(\theta\mid x) g(w(z) ; \theta). 
\end{aligned}
$$
Thus, we can retrieve the contracted posterior of the LHS by joint sampling as on the RHS (which is step 1 above) and then conditioning on \(w(z)\).

The label inherits the same reduction,
$$
\begin{aligned} 
y = P(H_1\mid x,z) & = \int 1_{\theta\in H_1}\,p(\theta\mid x,z)\,d\theta \\ 
& = C^{-1}\int 1_{\theta\in H_1}\,p(\theta\mid x) g(w(z) ; \theta)\,d\theta \\
& =:\;q^\star(w(z); x).
\end{aligned}
$$


 Thus the regression task in step 2 above provides a consistent estimator of \(q^\star\) under the joint sampling of step 1. Specifically, the marginal of the training pairs is
\[
p(w,y\mid x)=\int p(w\mid\theta)\,\mathrm{Bern}\big(y;1_{\theta\in H_1}\big)\,p(\theta\mid x)\,d\theta,
\]
and the population minimiser of either the MSE or the cross-entropy loss is the conditional expectation
\[
\arg\min_\tau \mathbb{E}\big[(y-q_\tau(w))^2\,\big|\,x\big]
\;=\;\mathbb{E}[y\mid w,x]
\;=\;\int 1_{\theta\in H_1}\,p(\theta\mid w,x)\,d\theta.
\]
Any consistent regression family \(q_\tau\) therefore recovers \(q^\star\) as \(S\to\infty\), and step 3 evaluates the exact label \(P(H_1\mid x,z)\) at any new \(z\) by a single regression evaluation. 

When \(w\) is only approximately sufficient, the same population minimiser applies, but it is now the projection
\[
\mathbb{E}[y\mid w,x]
\;=\;\mathbb{E}_{z\mid w,x}\!\big[P(H_1\mid x,z)\big],
\]
i.e. the smoothed label averaged over the residual \(z\)-information not captured by \(w\). The bias relative to the true \(P(H_1\mid x,z)\) is bounded by the conditional variance \(\mathrm{Var}\big(P(H_1\mid x,z)\,\big|\,w(z),x\big)\); the bias goes to zero as \(w\) approaches a sufficient statistic. 

In exponential-family models (e.g. the Bernoulli and Categorical examples of Section 3.1, 3.2) the sufficient statistic is canonical and finite-dimensional, so \(w\) can be chosen exactly. 

In partial-credit / IRT models with per-unit random effects no finite-dimensional sufficient statistic exists; a pragmatic choice for \(w\) is the per-\((\text{item},\text{time})\) response summary that the analyst would actually use post-trial. In our case, for out-of-7 outcomes, a suitable summary might just be the mean of the item responses at baseline and endline,
\[
   w_t(z^{(s)}) = m^{-1}\sum_{i=1}^m z_{it}^{(s)}.
\]
For categorical responses, a suitable summary might just be the proportion of the item responses above grade 3 at baseline and endline,
\[
   w_t(z^{(s)}) = m^{-1}\sum_{i=1}^m 1\{ z_{it}^{(s)} > 3\}.
\]


**Practical caveats.**

- *Finite \(S\).* The regression variance is \(O(S^{-1}\,\sigma^2/n_{\text{eff}}(W))\) where \(\sigma^2=\mathrm{Var}(Y\mid W)\le 1/4\) (Bernoulli) and \(n_{\text{eff}}(W)\) is the effective local sample size. @strong2014estimating recommend \(S\) in the low thousands.
- *Dimensionality.* Additive structure breaks down for \(\dim(W)\gtrsim 6\); GPs handle higher \(d\) but scale as \(O(S^3)\) without sparse approximations.
- *Uncertainty.* GP variants give a posterior credible band on \(\hat y(z)\); GAMs give standard errors from the link-scale Gaussian approximation.

---

# 7. Amortized inference workflow

Rather than learning a neural posterior of high dimensional model parameters, an initial step might be to learn the function

\[
    z \rightarrow 𝑃( H_1 \mid x,z ),
\]

because this is the key component in the decision rules and PPS above. Here \(x\) is today's fixed observed data, but it is undesirable to repeat learning for every \(x\). This motivates me to learn a neural function

\[
(x_{1:n},z_{1:m}) \rightarrow q_\phi(x_{1:n},z_{1:m}) \approx P(H_1 \mid x_{1:n},z_{1:m}) \in [0,1]
\]

for arbitrary data \(x\) with \(n\) data points today and for arbitrary future data \(z\) with remaining \(m\) data points until the end of the scheduled intervention. In my setting, the posterior always factorizes,
\[
    p(\theta \mid x,z) \propto p(\theta) \prod_{i=1}^n p(x_i \mid \theta) \prod_{i=1}^m p(z_i \mid \theta).
\]

I can re-frame this learning task as
$$
\begin{aligned}
q_\phi(x,z) 
& 
\approx \mathbb{E}_{\theta|x,z}(L(A,\theta)) 
\\
& 
=  \int L(A,\theta) p(\theta | x,z) d\theta,
\end{aligned}
$$
which makes clear that we focus on learning the loss weighted mean of the future posterior. This is decision-focussed learning, we learn what is required to make an optimal decision.

If we step back for a moment, for product likelihood models in the exponential family, we have for each data input \(\chi_i\)

\[
p(\chi_i|\theta)=h(x_i)\exp\{\eta(\theta)^\top T(\chi_i)-A(\theta)\},
\]

and so with a conjugate prior, the posterior is also in the exponential family and depends only on the sufficient statistic. In this case, if we condition on \(x\) and \(z\), we can re-express our target in terms of some function \(g\)
\[
P(H_1|x,z) = g\big(\sum_{i=1}^n T(x_i) + \sum_{i^*=1}^m T(z_{i^*})\big).
\]

We always have that the data \(x_i\) and \(z_{i^*}\) are permutation invariant, and this motivates that we want to learn end-to-end a two-stage neural model \( q_\psi, q_\tau \) that decomposes along the DeepSets theorem,

\[
(x,z) \mapsto q_\psi\bigg(\sum_{i=1}^n q_\tau(x_i) + \sum_{i^* =1}^m q_\tau(z_{i^*})\bigg) \approx P(H_1|x,z).
\]

In particular, this immediately provides generalisability to different current data sets of different sizes \(n\), and to different future data sets of different sizes \(m\).

There are three steps:

* Workflow Step 1: generating labelled training data
* Workflow Step 2: learning neural architectures
* Workflow Step 3: deployment, amortise PPS


---


# 8. Step 1: generating labelled training data

Generating training data is more challenging than in neural posterior estimation (NPE). We consider two distinct cases depending on whether the current data \(x\) is fixed or variable.

## Case A: Fixed \(x\) with access to the posterior

These are the same approaches as in Step 6.

## Case B: Variable \(x\)

In the more general setting, we want to learn a function that works for arbitrary interim data \(x\), not just a single fixed dataset. This enables true real-time decision-making as data accumulates. 

The key challenge is that for each training example, we need to approximate \(p(\theta \mid x^{(s)}, z^{(s)})\), but this posterior varies for every simulated dataset pair. We exploit that when we generate data from \(\theta^{(s)} \sim p(\theta)\), this parameter lies in a high-likelihood region for both \(x^{(s)}\) and \(z^{(s)}\) generated from it, making \(\theta^{(s)}\) a natural anchor point for approximating the joint posterior.

**Option 1: Local importance sampling**

This approach uses a local proposal distribution centered at the data-generating parameter \(\theta^{(s)}\) to approximate the joint posterior \(p(\theta \mid x^{(s)}, z^{(s)})\).

**Algorithm:**

For each \(s = 1, \ldots, S\), do:

1. Sample from prior: \(\theta^{(s)} \sim p(\theta) \). 
   
2. Generate current data: \( x^{(s)} \sim p(\cdot \mid \theta^{(s)}) \).

3. Generate future data: sample future observations from the same parameter \( z^{(s)} \sim p(\cdot \mid \theta^{(s)}) \). Both \(x^{(s)}\) and \(z^{(s)}\) are generated from \(\theta^{(s)}\), so \(\theta^{(s)}\) has high likelihood for both datasets.

4. Construct local proposal:
    \[
    q(\theta \mid \theta^{(s)}) = \mathcal{N}(\theta; \theta^{(s)}, \sigma^2 I)
    \]
    where \(\sigma^2\) controls the local exploration radius, e.g., \(\sigma = 0.1\) times the prior standard deviation.

5. Sample particles locally: Sample \(K\) particles around \(\theta^{(s)}\):
    \[
    \theta_k \sim q(\cdot \mid \theta^{(s)}), \quad k = 1, \ldots, K.
    \]
    Optionally include \(\theta^{(s)}\) itself as \(\theta_1 = \theta^{(s)}\) to ensure at least one high-likelihood particle.

6. Compute importance weights: Calculate weights for the joint posterior \(p(\theta \mid x^{(s)}, z^{(s)})\):
    $$
    \begin{aligned}
    w_k & = \frac{p(\theta_k \mid x^{(s)}, z^{(s)})}{q(\theta_k \mid \theta^{(s)})} \\
    & \propto \frac{p(\theta_k) \, p(x^{(s)} \mid \theta_k) \, p(z^{(s)} \mid \theta_k)}{q(\theta_k \mid \theta^{(s)})},
    \end{aligned}
    $$
    and normalize to
    \[
    \tilde{w}_k = w_k / \sum_{j=1}^K w_j.
    \]

7. Estimate the label by
   \[
    y^{(s)} := P(H_1 \mid x^{(s)}, z^{(s)}) \approx \sum_{k=1}^K \tilde{w}_k \, 1_{\theta_k \in H_1}.
   \]
   
**Option 2: Sequential Monte Carlo**
   
An alternative approach uses sequential importance sampling, which processes data incrementally and reuses the same \(K\) particles across both \(x^{(s)}\) and \(z^{(s)}\).
   
**Algorithm:**

For each \(s = 1, \ldots, S\), do:
   
   1. Sample generating parameter: Sample \(\theta^{(s)} \sim p(\theta)\).
   
   2. Initialize particle system: Create \(K\) particles by including \(\theta^{(s)}\) and sampling \(K-1\) additional particles from the prior:
      \[
      \{\theta_{s,1} = \theta^{(s)}, \theta_{s,2}, \ldots, \theta_{s,K}\}, \quad \theta_{s,k} \sim p(\theta) \text{ for } k \geq 2
      \]
      with initial uniform weights \(w_{s,k}^{(0)} = 1/K\) for all \(k\). Including \(\theta^{(s)}\) as a particle ensures at least one particle has high likelihood.
   
   3. Generate current data: Sample one dataset \(x^{(s)} \sim p(\cdot \mid \theta^{(s)})\).
   
   4. Update particle weights for \(x^{(s)}\): Compute importance weights to approximate \(p(\theta \mid x^{(s)})\) by evaluating how well each of the \(K\) particles explains the observed \(x^{(s)}\):
      $$
      \begin{aligned}
      w_{s,k}^{(x)} & \propto w_{s,k}^{(0)} \cdot p(x^{(s)} \mid \theta_{s,k}) \\
      \tilde{w}_{s,k}^{(x)} & = w_{s,k}^{(x)} / \sum_{j=1}^K w_{s,j}^{(x)}
      \end{aligned}
      $$
      The normalized weights \(\{\theta_{s,k}, \tilde{w}_{s,k}^{(x)}\}_{k=1}^K\) now represent \(p(\theta \mid x^{(s)})\). Resample particles if the effective sample size \(1/\sum_k (\tilde{w}_{s,k}^{(x)})^2\) is too small.
   
   5. Generate future data from posterior predictive: Sample one future dataset from the weighted particle system:
      \[
      z^{(s)} \sim \sum_{k=1}^K \tilde{w}_{s,k}^{(x)} \, p(\cdot \mid \theta_{s,k}).
      \]
      In practice, select a particle index \(\kappa \sim \text{Categorical}(\tilde{w}_{s,1}^{(x)}, \ldots, \tilde{w}_{s,K}^{(x)})\) and generate \(z^{(s)} \sim p(\cdot \mid \theta_{s,\kappa})\).
   
   6. Sequential update for \(z^{(s)}\): Update the weights again to approximate \(p(\theta \mid x^{(s)}, z^{(s)})\) by evaluating how well each particle explains the joint data:
      $$
      \begin{aligned}
      w_{s,k}^{(x,z)} & \propto \tilde{w}_{s,k}^{(x)} \cdot p(z^{(s)} \mid \theta_{s,k}) \\
      \tilde{w}_{s,k}^{(x,z)} & = w_{s,k}^{(x,z)} / \sum_{j=1}^K w_{s,j}^{(x,z)}
      \end{aligned}
      $$
   
   7. Estimate label:
      \[
      y^{(s)} = P(H_1 \mid x^{(s)}, z^{(s)}) \approx \sum_{k=1}^K \tilde{w}_{s,k}^{(x,z)} \, 1_{\theta_{s,k} \in H_1}
      \]

**Comparing Options 1 and 2**

I was motivated to consider an alternative to option 1 because in option 1, both \(x^{(s)}\) and \(z^{(s)}\) are generated from the same \(\theta^{(s)} \sim p(\theta)\). The joint distribution is
  \[
  p^{\text{prior}}(x,z) = \int p(x \mid \theta) p(z \mid \theta) p(\theta) \, d\theta.
  \]
However during deployment, for a fixed observed \(x\), we generate \(z \sim p(z \mid x)\) from the posterior predictive to evaluate PPS, and thus the corresponding joint distribution is tighter. In option 2, after generating \(x^{(s)} \sim p(\cdot \mid \theta^{(s)})\), the future data \(z^{(s)}\) is sampled from the posterior predictive \(p(z \mid x^{(s)})\) via the particle approximation in step 4. The joint distribution is
  \[
  p^{\text{post}}(z \mid x) p^{\text{prior}}(x) = \left[\int p(z \mid \theta) p(\theta \mid x) \, d\theta\right] p(x),
  \]
and this matches the deployment distribution.


---


# 9. Step 2: learning neural architectures

Given the training data \( (\theta^{(s)}, x^{(s)}, z^{(s)}, y^{(s)}) \), the current idea is to learn end-to-end two neural models \(q_\tau\) and \(q_\psi\) that respectively learn low-dimensional summaries and the decision target,

\[
(x^{(s)},z^{(s)}) \mapsto q_\psi\bigg(\sum_{i=1}^n q_\tau(x^{(s)}_i) + \sum_{i^* =1}^m q_\tau(z^{(s)}_{i^*})\bigg) \approx y^{(s)} = P(H_1|x^{(s)},z^{(s)}).
\]

Since \(y^{(s)}\) are probabilities, perhaps set the loss function to logit mean square error,

\[
\mathrm{arg\,min}_{\psi, \tau} \frac{1}{S} \sum_{s=1}^S \bigg[ q_\psi\bigg(\sum_{i=1}^n q_\tau(x^{(s)}_i) + \sum_{i^* =1}^m q_\tau(z^{(s)}_{i^*})\bigg) - \text{logit} \: y^{(s)} \bigg]^2.
\]

The specifics of our core application (Hope Groups intervention) are as follows: our data are ordered categorical responses, so each person \(i\) provides \(J\) item responses \(x_{ij} \in \{1, \ldots, K\}\) where \(K\) is the number of ordered categories. The typical scale of our data sets is \(N \approx 1000\) persons, \(J \approx 20\) items, \(K \approx 8\) categories. The data are permutation invariant over participants \(i\), but not over survey items \(j\).

**Option 1:  unrestricted embeddings**

Learn \(q_\psi\) and \(q_\tau\) without constraints, starting with a simple MLP architecture?

**Option 2:  minimum embeddings**

There are two simple, relevant models in the exponential family. First, the categorical model assuming that participants and responses are iid, and in this case the sufficient statistics are the \(K-1\) dimensional vector of counts, which can be written as the sum of \(K-1\) dimensional one-hot encodings for each \(i,j\). This suggests that the network \(q_\tau\) should be at least into \(\mathbb{R}^{K-1}\), and that the first \(K-1\) should be one-hot embeddings. 

**Option 3:  minimum restricted embeddings**

Second, the other extreme is a categorical model that assumes that participants are identical but responses are indpendent, in which case the sufficient statistics are the sum of \(JK-J\) dimensional one-hot encodings for each \(i\). This suggests that the network \(q_\tau\) should at most be into \(\mathbb{R}^{JK-J}\). 

**Additional points on the training strategy**

* Hold out a 20\% of the simulated examples \((\theta^{(s)}, x^{(s)}, z^{(s)}, y^{(s)})\) as test set
* Mini-batching: With \(S \sim 10^4\) to \(10^5\) training examples, use stochastic gradient descent with mini-batches of size 32-256.
* Since \(y^{(s)} \in [0,1]\) is a probability, set the output layer to sigmoid activation
* Vary \(n\) and \(m\) (number of persons in current and future data) within training examples to improve generalization

**Nested DeepSets architecture:**

Respect the hierarchical structure with person-level and dataset-level aggregation:
\[
q_\psi(x) = \rho_{\text{outer}}\left(\sum_{i=1}^N \rho_{\text{inner}}\left(\sum_{j=1}^J \phi_{\text{embed}}(x_{ij}, j)\right)\right)
\]
This first aggregates item responses within each person, then aggregates across persons.


---


# 10. Step 3: Deployment, amortising PPS 

After training, 

* for any observed data \(x_{1:n}\) we can compute \( \sum_{i=1}^n q_\tau(x^{(s)}_i) \). 
* give interim data \(x_{1:n}\), in practice we will estimate today's posterior \(p(\theta | x)\) with some Monte Carlo algorithm, and compute the first target objective, have we already won, \( 1\{ P(H_1 \mid x ) > \eta \). It is cheap to generate posterior predictions \(z^\star_{1:m}\). Then we can compute \( \sum_{i^\star=1}^m q_\tau(z^{(s)}_{i^\star}) \) for the corresponding \(m = N - n\).
* Next we can compute the PPS approximation

\[
q_\psi\bigg(\sum_{i=1}^n q_\tau(x^{(s)}_i) + \sum_{i^* =1}^m q_\tau(z^{(s)}_{i^*})\bigg).
\]


---


# 11. Tasks

* Generate functions to compute the exact \( PPS(x) \) in the case of simple Binomial and Categorial models
* Generate a case study that illustrates contracting posterior distributions with additional data, and evolving \( PPS(x) \) with additional data under the simple Binomial model. Illustrate stopping early for futility/safety and for efficacy. 
* Generate functions that implement Case A and B in generating labelled training data.
* Test these functions against the analytical formula under the simple Binomial and Categorial model.
* Test these functions against the analytical formula under the Categorial model.

# 12. References

<!--
Bibliography is in dev/amortised_decision_making.bib. Render with pandoc:
  pandoc dev/amortised_decision_making.md --citeproc \
    --bibliography dev/amortised_decision_making.bib \
    -o amortised_decision_making.html
The list below is auto-generated by --citeproc from the inline [@key] citations.
-->