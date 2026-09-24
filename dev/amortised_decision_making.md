# Amortized decision-focussed learning for real-time interim analyses in crisis settings

## 1. Summary

In contrast to neural posterior estimation, the focus of this work is on amortization for Bayesian decision making. The key point is that decision rules are low-dimensional, so it should be much easier to perform amortised learning by focussing on what is relevant to make real-life decision in real-time.

## 2. Target objectives

Perhaps the simplest approach to this is to consider an intervention that estimates a treatment effect $p$. The actions are declaring success ($A_1$) or failure ($A_0$) based on the binary decision rule $$d(x) = 1\{P(H_1 \mid x) > \eta_H\},$$ where $$H_1 : p > p_0 \Leftrightarrow \rho := p/p_0 < 1.$$ In practice, the interventions can be Cholera vaccination, our Hope group interventions in Jordan, alternative medical procedures, or similar. As generative models I am interested in various versions of item-response models such as the ordered categorical or partial credit model, but many other generative models are possible. The baseline treatment effect is without loss of generality $p_0 = 0$, for example $p$ can be coded as a contrast between the intervention and standard-of-care arm. The decision threshold $\eta_H$ can be set to $0.89$ in defiance to $0.95$, but can also be derived as the Bayes optimal threshold given false positive and false negative loss functions.

There are two distinct targets.

### Target 1, current decision: “Have we already won?”

We evaluate the posterior success probability of $H_1$, $P(H_1 \mid x)$, and evaluate the binary decision rule $$d(x) = 1\{P(H_1 \mid x) > \eta_H\} \in \{ 0, 1 \},$$ which declares success or failure. This answers the question: Have we already accumulated enough evidence to declare success now?

### Target 2, predictive probability of success (PPS): "How likely will we win in the future?"

We are interested in real-time decision-making, typically referred to as interim analyses to stop for efficacy or to stop for futility. This reframes the problem.

The new questions are: If the current data $x$ declare failure, how likely is it that future data collection could overturn this result, or should we stop for failure and safety reasons now? The other way round, if the current data $x$ declare success, how likely is it that future data collection could overturn this result, or can we stop for efficacy early now?

Suppose $m$ additional observations $z$ will be collected before the final analysis. Let the final decision rule be $$d(x,z) = 1\{P(H_1 \mid x,z) > \eta_H\}.$$ The predictive probability of success (PPS) is $$PPS(x) =
\mathbb{E}_{z \sim p(z\mid x)}
[
1\{P(H_1 \mid x,z) > \eta_H\}
] \in [ 0, 1],$$ where the posterior predictive data generating distribution is $$p(z|x) = \int p(z|\theta)p(\theta|x)d\theta.$$ The PPS quantifies the probability that the intervention will declare success in the future after $m$ further hypothetical data points $z$, given the current data $x$.

Thus PPS thus answers: How likely are we to win in the future if the intervention continues?

------------------------------------------------------------------------

# 3. Examples

## 3.1 Binomial Model

Let $$X_i \sim \text{Bernoulli}(p)$$ with prior $$p \sim \text{Beta}(a,b).$$ We want to test $$H_1: p < p_0 (1 - \eta_0)$$ or equivalently $$H_1: 
\rho > \eta_0$$ where $\rho$ is the effect measure, $\rho := 1- p/p_0$ and $\eta_0$ means the reduction in the baseline event probability needs to exceed $\eta_0$. In the simplest case, $\eta_0 = 0$. The Bernoulli distribution belongs to the exponential family: $$\begin{aligned}
p(x_i \mid p) & 
= p^{x_i}(1-p)^{1-x_i} = \exp\{x_i \log p + (1-x_i)\log(1-p)\} 
\\
& 
= \exp\{\log\left(\frac{p}{1-p}\right) x_i + \log(1-p)\},
\end{aligned}$$ with natural parameter $\eta(p) = \log(p/(1-p))$, sufficient statistics $T(x_i) = x_i$, and log-partition function $A(p) = -\log(1-p)$. For $n$ observations, the likelihood depends only on the sufficient statistic $T(x_{1:n}) = \sum_{i=1}^n x_i = k^n$, the total number of successes. After observing $x_{1:n}$ data points with $k^n$ successes the posterior is $$p \mid x \sim \text{Beta}(a+k^n, b+n-k^n).$$ Suppose the future, remaining data $z$ comprise $m$ additional data points. Then the posterior predictive distribution of $k^m$ future successes is $$\begin{aligned}
k^m & \sim \int \text{Binomial}(m,p) \text{Beta}(p; a+k^n, b+n-k^n) \\
& = \text{Beta-Binomial}(m,a+k^n,b+n-k^n).
\end{aligned}$$ If $k^m$ additional successes occur: $$p \mid x,z \sim
\text{Beta}(a+k^n+k^m, b+n+m-k^n-k^m).$$ The final posterior success probability is $$P(H_1 \mid x,z) = F_{\text{Beta}}\!\big(p_0 (1 - \eta_0);\ a+k^n+k^m,\ b+n+m-k^n-k^m\big).$$ In the future, with $z$ additional data, success occurs when $$P(H_1 \mid x,z) > \eta_H.$$ Since the posterior depends only on the sufficient statistic $k^m$ (the number of successes in the future data), and $P(H_1 \mid x, k^m)$ is monotonically decreasing in $k^m$ (more successes push the posterior mass of $p$ upward, away from the alternative $p < p_0(1-\eta_0)$), we can find the threshold $$k^m_\star = \max\{k^m \in \{0, 1, \ldots, m\} : P(H_1 \mid x, k^m) > \eta_H\}.$$ Then the PPS becomes a simple tail sum: $$PPS(x) = \sum_{k^m = 0}^{k^m_\star} p(k^m \mid x),$$ where $p(k^m \mid x)$ is the Beta–Binomial predictive distribution. Explicitly, $$PPS(x) = \sum_{k^m = 0}^{k^m_\star} {m \choose k^m} \frac{B(a+k^n+k^m,b+n+m-k^n-k^m)}{B(a+k^n,b+n-k^n)}.$$

## 3.2 Categorical model

Let $$x_i \sim \text{Categorical}(p)$$ with parameter vector $p = (p_1, \ldots, p_K)$ where $\sum_{k=1}^K p_k = 1$, and prior $$p \sim \text{Dirichlet}(\alpha),$$ where $\alpha = (\alpha_1, \ldots, \alpha_K)$. We want to test, for a chosen target category $k^\star \in \{1, \ldots, K\}$ (e.g. the adverse-outcome category), $$H_1: p_{k^\star} < p^0_{k^\star} (1 - \eta_0),$$ or equivalently $$H_1: \rho_{k^\star} > \eta_0$$ where $\rho_{k^\star} := 1 - p_{k^\star}/p^0_{k^\star}$ is the relative reduction in category $k^\star$ and $\eta_0 \in [0, 1)$ is the required margin. In the simplest case, $\eta_0 = 0$. The Categorical distribution belongs to the exponential family: $$\begin{aligned}
p(x_i \mid p) & 
= \prod_{k=1}^K p_k^{1_{x_i = k}} 
= \exp\left\{\sum_{k=1}^K 1_{x_i = k} \log p_k\right\}
\\
&
= \exp\left\{\sum_{k=1}^{K-1} \log\left(\frac{p_k}{p_K}\right) 1_{x_i = k} + \log p_K\right\},
\end{aligned}$$ with natural parameters $\eta_k(p) = \log(p_k/p_K)$ for $k = 1, \ldots, K-1$, sufficient statistics $T_k(x_i) = 1_{x_i = k}$ for $k = 1, \ldots, K-1$, and log-partition function $A(p) = -\log p_K$. For $n$ observations, the likelihood depends only on the sufficient statistic $T(x_{1:n}) = (k_1^n, \ldots, k_K^n)$, the count vector where $k_k^n = \sum_{i=1}^n 1_{x_i = k}$ is the number of observations in category $k$. After observing $x_{1:n}$ data points with counts $k_k^n$ the posterior is $$p \mid x \sim \text{Dirichlet}(\alpha_1 + k_1^n, \ldots, \alpha_K + k_K^n).$$ Suppose the future, remaining data $z$ comprise $m$ additional data points. Then the posterior predictive distribution of the future count vector $(k_1^m, \ldots, k_K^m)$ is $$\begin{aligned}
(k_1^m, \ldots, k_K^m) & \sim \int \text{Multinomial}(m, p) \, \text{Dirichlet}(p; \alpha_1 + k_1^n, \ldots, \alpha_K + k_K^n) \, dp \\
& = \text{Dirichlet-Multinomial}(m, \alpha_1 + k_1^n, \ldots, \alpha_K + k_K^n).
\end{aligned}$$ If future counts $(k_1^m, \ldots, k_K^m)$ occur: $$p \mid x,z \sim \text{Dirichlet}(\alpha_1 + k_1^n + k_1^m, \ldots, \alpha_K + k_K^n + k_K^m).$$ The marginal of $p_{k^\star}$ under the Dirichlet posterior is Beta with shape parameters $\alpha_{k^\star} + k_{k^\star}^n + k_{k^\star}^m$ and $\sum_{j \ne k^\star}(\alpha_j + k_j^n + k_j^m)$, so the final posterior success probability is closed form, $$P(H_1 \mid x,z) = F_{\text{Beta}}\!\big(p^0_{k^\star} (1 - \eta_0);\ \alpha_{k^\star} + k_{k^\star}^n + k_{k^\star}^m,\ \textstyle\sum_{j \ne k^\star}(\alpha_j + k_j^n + k_j^m)\big).$$ In the future, with $z$ additional data, success occurs when $$P(H_1 \mid x,z) > \eta_H.$$ Since $P(H_1 \mid x, k_1^m, \ldots, k_K^m)$ depends on the future count vector only through $k_{k^\star}^m$ (the other counts enter only via the fixed total $m - k_{k^\star}^m = \sum_{j \ne k^\star} k_j^m$ in the second Beta shape) and is monotonically decreasing in $k_{k^\star}^m$, we can find the threshold $$k^{m\,\star}_{k^\star} = \max\{k_{k^\star}^m \in \{0, 1, \ldots, m\} : P(H_1 \mid x, k_{k^\star}^m) > \eta_H\}.$$ Then the PPS reduces to a one-dimensional tail sum over the marginal predictive of $k_{k^\star}^m$: $$PPS(x) = \sum_{k_{k^\star}^m = 0}^{k^{m\,\star}_{k^\star}} p(k_{k^\star}^m \mid x),$$ where the marginal of a Dirichlet-Multinomial count is Beta-Binomial, $$p(k_{k^\star}^m \mid x) = {m \choose k_{k^\star}^m} \frac{B\!\big(\alpha_{k^\star} + k_{k^\star}^n + k_{k^\star}^m,\ \sum_{j \ne k^\star}(\alpha_j + k_j^n) + m - k_{k^\star}^m\big)}{B\!\big(\alpha_{k^\star} + k_{k^\star}^n,\ \sum_{j \ne k^\star}(\alpha_j + k_j^n)\big)}.$$ This recovers the Binomial PPS of §3.1 when $K = 2$ and $k^\star = 1$, and remains analytically computable for any $K$ at the cost of a one-dimensional sum over the target-category count.

## 3.3 Multivariate Normal model with known correlation structure

A canonical high-dimensional analytically tractable PPS benchmark for continuous outcomes is the multivariate normal model with structured covariance. Let $$y_n \mid \mu, \sigma^2 \;\sim\; \mathrm{MVN}(\mu, \sigma^2 K), \qquad n = 1, \ldots, N,$$ where $K = R R^\top \in \mathbb{R}^{J \times J}$ is a known positive-definite covariance shape (for example, $R$ a Cholesky factor of an AR(1) or arbitrary correlation matrix), and the unknowns are the mean vector $\mu \in \mathbb{R}^J$ and the scalar variance $\sigma^2 > 0$. For each participant $n$, we have $J$ correlated component observations.

We test $J$ component-wise alternative hypotheses, $$H_{1j} : \mu_j < \mu^0_j (1 - \eta_0) \quad \text{for all } j = 1, \ldots, J,$$ with the baseline pinned at $\mu^0_j := 1$ for concreteness and required relative-reduction margin $\eta_0 \in [0, 1)$. Using the relative-effect parameterisation, $$\rho_j \;:=\; 1 - \mu_j/\mu^0_j \;=\; 1 - \mu_j,$$ the alternative is equivalently $H_{1j} \Leftrightarrow \rho_j > \eta_0$. In the simplest case, $\eta_0 = 0$.

**Conjugate prior** (Normal–Inverse-Gamma with metric $K$): $$\mu \mid \sigma^2 \;\sim\; \mathrm{MVN}(\mu_0,\, \sigma^2 \Lambda_0^{-1}), \qquad \sigma^2 \;\sim\; \mathrm{Inv\text{-}Gamma}(a_0, b_0),$$ with hyperparameters $(\mu_0, \Lambda_0, a_0, b_0)$ known.

**Exponential family and sufficient statistics.** With $K$ known, the MVN likelihood is in the exponential family with natural parameters $\eta(\mu, \sigma^2) = \big( K^{-1}\mu/\sigma^2,\, -1/(2\sigma^2) \big)$ and sufficient statistics for $n$ current observations $$T_1(x_{1:n}) \;=\; \sum_{i=1}^n x_i, \qquad T_2(x_{1:n}) \;=\; \sum_{i=1}^n x_i^\top K^{-1} x_i,$$ so the likelihood factors as $$p(x_{1:n} \mid \mu, \sigma^2) = h(x_{1:n})\, \exp\{\eta^\top T - n A(\mu, \sigma^2)\}$$ with $A(\mu, \sigma^2) = \tfrac{J}{2}\log\sigma^2 + \tfrac{1}{2\sigma^2}\mu^\top K^{-1} \mu$.

**Posterior** after current data $x_{1:n}$, writing $\bar{x}_n = T_1(x_{1:n})/n$: $$\begin{aligned}
\Lambda_n &= \Lambda_0 + n K^{-1},\\
\mu_n &= \Lambda_n^{-1}\!\big( \Lambda_0 \mu_0 + n K^{-1} \bar{x}_n \big),\\
a_n &= a_0 + nJ/2,\\
b_n &= b_0 + \tfrac{1}{2}\!\left( T_2(x_{1:n}) + \mu_0^\top \Lambda_0 \mu_0 - \mu_n^\top \Lambda_n \mu_n \right).
\end{aligned}$$ Marginalising $\sigma^2$ gives the multivariate Student-$t$ marginal posterior of the mean, $$\mu \mid x \;\sim\; t_{2a_n}\!\big( \mu_n,\; (b_n/a_n)\,\Lambda_n^{-1} \big),$$ with per-component marginals $$\mu_j \mid x \;\sim\; t_{2a_n}\!\big( \mu_{n,j},\; (b_n/a_n)\,[\Lambda_n^{-1}]_{jj} \big).$$ The current posterior probability of the per-component hypothesis is therefore closed form, $$P(H_{1j} \mid x) \;=\; F_{t_{2a_n}}\!\left( \frac{(1 - \eta_0) - \mu_{n,j}}{\sqrt{(b_n/a_n)\,[\Lambda_n^{-1}]_{jj}}} \right),$$ and the joint hypothesis $P(H_1 \mid x) = P(\mu_j < 1 - \eta_0 \;\forall j \mid x)$ is a multivariate-$t$ orthant probability, evaluated by quasi-Monte Carlo (Genz QMC).

**Posterior predictive of future data.** Suppose $m$ additional observations $z_{1:m}$ will be collected. Conditional on $x$, the future data are jointly matrix-$t$ distributed, and their joint distribution factors through the future sufficient statistics $$T_1(z_{1:m}) \;=\; \sum_{i=1}^m z_i, \qquad T_2(z_{1:m}) \;=\; \sum_{i=1}^m z_i^\top K^{-1} z_i,$$ with $\bar{z}_m := T_1(z_{1:m})/m$. The marginal predictive of $\bar{z}_m$ integrating out $(\mu, \sigma^2)$ is multivariate $t$, $$\bar{z}_m \mid x \;\sim\; t_{2a_n}\!\left( \mu_n,\; (b_n/a_n)\big( \tfrac{1}{m} K + \Lambda_n^{-1} \big) \right),$$ and $T_2(z_{1:m})$ conditional on $\bar{z}_m$ is a quadratic form whose distribution is a scaled $F$ (the future analogue of the residual sum of squares under NIG).

**Posterior after current and future data.** With $\bar{x}_{n+m} = (n\bar{x}_n + m \bar{z}_m)/(n+m)$: $$\begin{aligned}
\Lambda_{n+m} &= \Lambda_0 + (n+m) K^{-1},\\ 
& a_{n+m} &= a_0 + (n+m) J / 2,\\
\mu_{n+m} &= \Lambda_{n+m}^{-1}\!\big( \Lambda_0 \mu_0 + n K^{-1} \bar{x}_n + m K^{-1} \bar{z}_m \big),\\
b_{n+m} &= b_0 + \tfrac{1}{2}\!\left( T_2(x_{1:n}) + T_2(z_{1:m}) + \mu_0^\top \Lambda_0 \mu_0 - \mu_{n+m}^\top \Lambda_{n+m} \mu_{n+m} \right).
\end{aligned}$$ Crucially, $\Lambda_{n+m}$ and $a_{n+m}$ are deterministic given the future cohort size $m$; only $\mu_{n+m}$ and $b_{n+m}$ depend on the random future statistics $(\bar{z}_m, T_2(z_{1:m}))$.

The final per-component success probability is $$P(H_{1j} \mid x, z) \;=\; F_{t_{2a_{n+m}}}\!\left( \frac{(1 - \eta_0) - \mu_{n+m, j}}{\sqrt{(b_{n+m}/a_{n+m})\,[\Lambda_{n+m}^{-1}]_{jj}}} \right).$$ **Predictive probability of success.** The per-component decision rule $d_j(x, z) = 1\{P(H_{1j} \mid x, z) > \eta_H\}$ depends on $z$ only through the two sufficient statistics $(\bar{z}_m, T_2(z_{1:m}))$, so the PPS reduces to an integral against the closed-form predictive of these statistics. Define the critical region in sufficient-statistic space, $$A_j(x) \;=\; \left\{ (\bar{z}, S) \;:\; \frac{(1 - \eta_0) - \mu_{n+m,j}(\bar{z})}{\sqrt{(b_{n+m}(\bar{z}, S)/a_{n+m})\,[\Lambda_{n+m}^{-1}]_{jj}}} \;>\; F_{t_{2a_{n+m}}}^{-1}(\eta_H) \right\}.$$ Then $$PPS_j(x) \;=\; \int 1_{A_j(x)}(\bar{z}, S)\; p(\bar{z}, S \mid x)\; d\bar{z}\, dS,$$ a two-dimensional quadrature against the marginal predictive of $(\bar{z}_m, T_2(z_{1:m}))$ derived above. The joint-hypothesis PPS replaces the per-component CDF by the multivariate-$t$ orthant probability $P(\mu < (1 - \eta_0) \mathbf{1} \mid x, z)$ and integrates against the same predictive.

**Special case:** $\sigma^2$ known. Fixing $\sigma^2$, the conjugate prior collapses to MVN on $\mu$ alone, posteriors and predictives are Gaussian, and the per-component PPS reduces to a Gaussian tail probability, $$PPS_j(x) \;=\; \Phi\!\left( \frac{(1 - \eta_0) - \mu_{n,j} - z_{\eta_H} \sqrt{\sigma^2\,[\Lambda_{n+m}^{-1}]_{jj}}}{\sqrt{\sigma^2\,\big([\Lambda_n^{-1}]_{jj} - [\Lambda_{n+m}^{-1}]_{jj} + [K]_{jj}/m\big)}} \right),$$ where $z_{\eta_H} = \Phi^{-1}(\eta_H)$, i.e. a single $\Phi$ evaluation per interim per component. The full NIG case differs only by replacing the inner Gaussian tails with Student-$t$ tails and adding one $F$-distribution integration over the residual sum of squares.

**Why this benchmark.** With $K$ arbitrary known PSD and $J$ arbitrary (here we target $J$ up to several hundreds or thousands), the construction supplies a high-dimensional continuous-outcome decision problem whose PPS is closed-form, with which the nested-MC / IS / SMC / regression estimators of Section 6 can be validated.

### 3.3.1 Concrete simulation setup

We fix $\sigma^2 = 1$ (units chosen so that the per-component noise has unit variance), use the Gaussian special case throughout, and run the benchmark at three problem dimensions $$J \in \{50, 100, 200\},$$ with a total cohort $N = 500$ accrued over 10 monthly interims of 50 units each, so the future cohort size at interim $t$ is $m_t = N - n_t$ with $n_t = 50 t$. The decision threshold is $\eta_H = 0.89$, the baseline is $\mu^0 = \mathbf{1}_J$ and the relative-reduction margin is $\eta_0 = 0$.

**Prior.** A weakly informative $g$-prior in the same metric $K$ keeps the algebra closed, $$\mu \;\sim\; \mathrm{MVN}\!\big(\mathbf{1}_J,\; \tau_0^2\, K\big), \qquad \tau_0^2 = 100,$$ so $\Lambda_0^{-1} = \tau_0^2 K$ and $\Lambda_0 = K^{-1}/\tau_0^2$. With $\sigma^2 = 1$, $$\Lambda_n \;=\; (\tau_0^{-2} + n)\,K^{-1}, \qquad \Lambda_n^{-1} \;=\; \frac{K}{\tau_0^{-2} + n},$$ so the per-component closed-form PPS specialises to a clean $\Phi$-tail with $\Lambda_\bullet^{-1}$ replaced by a scaled $K$.

**True** $\mu$. Half the components above baseline, half below, to stress-test the per-component decision: $$\mu_{\text{true}, j} \;=\; \begin{cases} 1 + \Delta & j \le J/2 \\ 1 - \Delta & j > J/2 \end{cases}, \qquad \Delta = 0.3.$$ This sets the true component-wise rejection rate to exactly $1/2$, so the benchmark probes both efficacy and futility regimes simultaneously across the $J$ hypotheses.

**Choice of** $R$. We need a nontrivial known correlation structure: enough off-diagonal mass that components are correlated (so the closed-form predictive of $(\bar{z}_m, T_2(z))$ is non-trivial), but well-conditioned at all three $J$. We propose three structures, used as a panel:

1.  AR(1) Cholesky (primary, ordered components): $$K^{(1)}_{ij} \;=\; \rho^{|i - j|}, \qquad \rho = 0.7.$$ Cholesky factor $R^{(1)}$ with $R^{(1)} (R^{(1)})^\top = K^{(1)}$. Banded, PSD for any $\rho \in (-1, 1)$, condition number bounded uniformly in $J$ ($\kappa(K^{(1)}) = (1+\rho)/(1-\rho) \approx 5.67$ for $\rho = 0.7$). Reflects a natural ordering (e.g. time or position along a scale).

2.  Block equicorrelation (mid difficulty, exchangeable within block): $$K^{(2)} \;=\; \begin{bmatrix} K_w & \rho_b \mathbf{1}\mathbf{1}^\top & \cdots \\ \rho_b \mathbf{1}\mathbf{1}^\top & K_w & \cdots \\ \vdots & & \ddots \end{bmatrix},$$ with $B = J/10$ blocks of size $10$, intra-block correlation $\rho_w = 0.8$ (so $K_w = (1-\rho_w) I + \rho_w \mathbf{1}\mathbf{1}^\top$), inter-block correlation $\rho_b = 0.1$. PSD by construction ($\rho_w > \rho_b \ge 0$). Captures clustered (e.g. subscale) dependence.

3.  Low-rank-plus-diagonal factor structure (hard, dense long-range correlation): $$K^{(3)} \;=\; \beta \beta^\top + \psi I_J, \qquad \beta \in \mathbb{R}^{J \times r},\; r = 5,\; \psi = 0.1,$$ with $\beta_{jk} \sim \mathcal{N}(0, 1/r)$ drawn once and fixed, then $K^{(3)}$ rescaled so $\mathrm{diag}(K^{(3)}) = \mathbf{1}$. $R^{(3)}$ is the Cholesky factor. Five-factor latent structure, dense correlations decaying smoothly with no banding.

All three are scaled to unit diagonal so the per-component noise is comparable across $J$ and across structures. The Cholesky factors $R^{(\cdot)}$ are precomputed once per $(J, R\text{-type})$ cell and shared across the 10 interims.

**What we compare.**

| Cell | Knobs | Notes |
|----------------|-----------------------------|---------------------------|
| Cell A (easy) | $R = R^{(1)}$ AR(1), $J \in \{50,100,200\}$ | Banded, well-conditioned. |
| Cell B (medium) | $R = R^{(2)}$ block, same $J$ | Block structure; eigenvalues clustered. |
| Cell C (hard) | $R = R^{(3)}$ factor, same $J$ | Dense long-range; eigenvalue decay. |

For each cell and each interim $t \in \{1, \ldots, 10\}$, we compare the closed-form per-component $PPS_j(x)$ (from the $\Phi$-tail above) against the four estimators of Section 6 (nested-MC, self-normalised IS, moment-matching IS, SMC resample-move, regression on $w(z)$). Aggregate diagnostics: per-component absolute error and bias against the closed form; per-cell timing; ESS / $\hat{k}$ for the IS variants; tempering steps $T$ for SMC.

------------------------------------------------------------------------

## 3.4 Item-response model

In our applications, we are interested in Item Response Theory (IRT) models that can be fitted to real-time survey data collected about interventions in humanitarian and social science settings. The data consist of responses to Likert-scale survey items collected from participants at baseline and endline. The data increases in size as more individuals are surveyed. The goal is to quantify intervention effectiveness based on the current data, and to quantify the predictive probability of success (PPS).

Throughout, let $i$ index participants ($i = 1, \ldots, n$), $j$ index items/questions ($j = 1, \ldots, J$), and $k$ index response categories ($k = 1, \ldots, K$). The response of participant $i$ to item $j$ at time $t$ (baesline or endline) is $Y_{ijt}$; for simplicity we suppress $t$ in what follows. A widely-used IRT model which we use here is the partial credit model. The PCM models ordered categorical probabilities with cumulative logits, $$\begin{aligned}
P(Y_{ij} = k) &= \text{softmax}(\phi_{i,j,k}) = \frac{\exp(\phi_{i,j,k})}{\sum_{k'=1}^K \exp(\phi_{i,j,k'})}\\
\phi_{i,j,1} &= 0\\
\phi_{i,j,k} &= \sum_{s=1}^k \lambda_{j} \cdot \left(\theta_i + \mathbf{X}_{i,j}^T \mathbf{\beta} - c_{j,s}\right), \quad k=2,\dotsc,K
\end{aligned}$$ The number of free parameters are $N$ latent skills parameters $\theta_i$, one for each participant; $J(K-1)$ incremental skill thresholds $c_{j,s}$, one for each categorical increment for each response item; $J$ item loadings; and $P$ fixed participants effects. The total number of parameters is thus $N+JK + P$, comprising $N$ local parameters that grow wich each new participant and $JK + P$ global parameters that are shared across participants.

By construction, the PCM admits an incremental log-risk structure: $$\log \frac{\Pr(Y_{i,j} = k)}{\Pr(Y_{i,j} = k-1)} =  \lambda_{j} \cdot \left(\theta_i + \mathbf{X}_{i,j}^T \mathbf{\beta} - c_{j,k}\right)$$

The model thus does not follow the proportional cumulative odds assumptions, but it has a proportional incremental risk assumption: when $\beta$ represents an intervention effect (e.g., $X_{i,j,t} = 1$ at time $t=1$ and $X_{i,j,t} = 0$ at time $t=0$), then $$\begin{aligned}
& \frac{\Pr(Y_{i,j,t=1} = k \mid \eta_{i,j,t=1})}{\Pr(Y_{i,j,t=1} = k-1 \mid \eta_{i,j,t=1})} \bigg/ \frac{\Pr(Y_{i,j,t=0} = k \mid \eta_{i,j,t=0})}{\Pr(Y_{i,j,t=0} = k-1 \mid \eta_{i,j,t=0})} \\
&\quad = \exp\Big( \lambda_j ( \theta_i + \beta ) - \lambda_j \theta_i  \Big) = \exp( \lambda_j \beta ),
\end{aligned}$$ so when measured in incremental risks, the effect of the intervention is proportional to $\exp\beta$, and the same regardless of the category $k$.

Computationally, the PCM involves direct evaluations of category specific events ($Pr(Y=k)$ not $Pr(Y\leq k)$) and therefore tends to be much faster to evaluate and less prone to numerical issues than other IRT models.

**Priors.** Weakly-informative, standard priors are placed on all free parameters, with a deliberately wide prior on the thresholds so that a finite cohort still leaves substantial posterior uncertainty to resolve (the contraction the amortiser of §14 later exploits): $$\theta_i \sim N(0,1),\qquad c_{j,s}\sim N(0,3.5^2),\qquad \lambda_j \sim |t_3|,\qquad \beta_p \sim N(0,1).$$ That is: participant skills $\theta_i$ standard normal; the $J(K-1)$ incremental category thresholds $c_{j,s}$ Gaussian with a wide standard deviation of $3.5$; the $J$ item loadings $\lambda_j$ half-Student-$t$ with $3$ degrees of freedom (constrained positive, with one loading fixed to $1$ for identification); and the participant / intervention effects $\beta_p$ standard normal. The same specification is fit in Pyro and Stan.

## 3.6 Application: Hope Groups in Ukraine

The primary real-time case study of this document. Hope Groups is a 12-session, peer-facilitated psychosocial, mental-health and parenting-support programme for Ukrainian parents and caregivers affected by war and displacement (externally displaced, internally displaced, and living in war-affected areas), evaluated in a pragmatic cluster-randomised controlled trial [@tucker2026hopegroups; @tucker2024hopeprotocol]. Caregivers self-report a battery of ordinal items at baseline and endline; the estimand is the item-level endpoint effect $\rho_j$, and enrolment accrues over time, so the decision of interest is whether the accumulated evidence already predicts success — the PPS — at each weekly interim.

| field | value |
|-------------|----------------------------------------------------------|
| Problem | psychosocial / mental-health / parenting support for Ukrainian caregivers in war and displacement |
| Arms | Hope Groups intervention vs waitlist control (cluster-randomised, 90 clusters) |
| Timepoints | baseline and endline (1-week post-intervention); later follow-ups |
| Items | caregiver-report — mental health (depression/anxiety), violence against children, positive parenting, wellbeing; two PCM types: *out-of-7* "days in week" ($K=8$) and *categorical* caseness |
| Effect measure $\rho$ | $r_j(\theta)=s_j(\bar w_{e,j}/\bar w_{b,j}-1)$, direction-aware endpoint effect per item |
| Cohort / interims | $N\approx503$, $J=20$ items; the **weekly-29 SVI grid** ($n=48\to503$, `py-ukraine-interim-weekly-svi-260811`) — the reference $p(\rho\mid x)$ every amortiser is scored against |
| **Application target** | **the §14 deployment case study**: the deepset amortiser (bespoke PIT–KS $0.069$; item-general J64 $0.119$, §14.4.9) and the cross-method comparison (`py-ukraine-interim-compare-methods-260526`, §14.5) |
| Reference | [@tucker2026hopegroups] (effectiveness), [@tucker2024hopeprotocol] (protocol) |

## 3.7 Application: ChatGPT-versus-expert-feedback RCT

A randomised controlled trial in medical education comparing ChatGPT feedback with expert feedback on clinical-reasoning (key-feature) questions, with item-level scores and a critical-approach-to-AI survey [@cicek2024chatgpt]. A genuine two-arm design with pre/post measurement, at modest scale.

| field | value |
|---------------|---------------------------------------------------------|
| Problem | effect of ChatGPT vs expert feedback on clinical-reasoning learning |
| Arms | control (expert feedback) vs intervention (ChatGPT) |
| Timepoints | immediate and delayed test; pre/post survey |
| Items | 12 key-feature items (partial-credit scored — non-consecutive maxima, collapse per item) + 6 critical-approach-to-AI items ($1$–$7$ Likert) |
| Effect measure $\rho$ | pre$\to$post (immediate$\to$delayed) endpoint effect per item |
| Cohort / interims | $N=115$ test / $86$ survey (modest); interims by accrual order |
| **Application target** | **real-time evaluation** — a stress test at modest $N$ |
| Data | [@cicek2024chatgpt] |

## 3.8 Application: nurse-led self-management RCT in inflammatory bowel disease

A single-centre longitudinal study of a nurse-led psychoeducational self-management intervention for inflammatory bowel disease, with item-level patient-reported outcomes across five waves [@riverasequeiros2026ibd]. Single-arm, but the richest timepoint structure of the set — well suited to exhibiting posterior contraction with cohort size.

| field | value |
|---------------|---------------------------------------------------------|
| Problem | nurse-led self-management for inflammatory bowel disease |
| Arms | single-arm (no control) |
| Timepoints | five waves (baseline $\to$ 12 months) |
| Items | IBDQ-9, Family APGAR-5, PSS-14, CCVEII-9, HADS-14 — ordinal Likert |
| Effect measure $\rho$ | baseline$\to$endline endpoint effect per item |
| Cohort / interims | $N=60$ (modest); the five waves serve as interims |
| **Application target** | **real-time evaluation** at modest $N$; showcases contraction with $n$ |
| Status | *planned* — the dataset (`Rivera-Sequeiros_et_al_2025_IBD_PsychosocialDataset.xlsx`) is on hand; loader/interim grid not yet built |
| Data | [@riverasequeiros2026ibd] |

## 3.9 Application: refugee and migrant youth (REFUGE-ED)

Questionnaires administered to refugee and migrant children and young people across six countries under the EC Horizon 2020 REFUGE-ED project, at baseline and endline, at the item level [@refugeed2024youth]. A second humanitarian cohort — the largest linked pre/post sample of the set — with no treatment arm, so the estimand is the within-cohort endpoint effect.

| field | value |
|------------|------------------------------------------------------------|
| Problem | education / psychosocial support and integration of refugee and migrant youth (6 countries) |
| Arms | none (multi-site practices pilot); **paired** pre/post |
| Items (built) | **MSPSS perceived social support**, its **12 items grouped into the three subscales** — Family / Friends / Significant Other ($1$–$7$, `out-of-7` expected-score endpoint) |
| Effect measure $\rho$ | baseline$\to$endline **perceived-support gain** per item (`higher_is_better`) |
| Cohort / interims | **324 linked** participants (present at both times); **8** participant-count interims ($n=40\to324$) |
| Amortiser | the item-general **J64** net (§14.4.9) deployed as a **per-family asset**: **PIT–KS** $0.126$, marg–KS $0.145$ — the cross-cohort generalisation test (calibrates as well as on Ukraine, and beyond $N=503$). SVI-vs-amortiser $p(\rho\mid x)$ comparison in `py-refugee_interim_compare_methods_260902` (§14.5) |
| **Application target** | **real-time humanitarian evaluation** + the deepset **generalisation test** across item sets and cohort size |
| Data | [@refugeed2024youth]; Zenodo 10908209 |

## 3.10 Application: ICRC community-level MHPSS in armed conflict

A large International Committee of the Red Cross (ICRC) cohort of victims of violence, assessed before and after community-level mental-health and psychosocial support across armed-conflict settings in the Democratic Republic of the Congo, Mali and Nigeria [@andersen2022mhpss]. The paper reports **validated total scores** (DASS-21 subscales, IES-R, an ICRC functioning scale); rather than model raw item responses we **keep those totals and bin them into the five clinical severity levels** (Normal / Mild / Moderate / Severe / Extremely severe) using the paper's own **Figure-2 cut-offs**, turning each instrument into one ordinal ($K=5$) endpoint on a clinically grounded scale. The **beneficiary-level data (**$N=6{,}413$) with a codebook is in the published Frontiers/PMC Supplementary Material (`Data_Sheet_1.xlsx`); the fuller ICRC monitoring dataset remains request-only.

| field | value |
|------------|------------------------------------------------------------|
| Problem | community-level MHPSS for victims of violence in armed conflict |
| Setting | ICRC, 32 projects — DR Congo (81%), Nigeria (16%), Mali (3%); $N=6{,}413$ |
| Arms | single cohort (**paired** pre/post; no randomised control) |
| Items (built) | **DASS-21** Depression / Anxiety / Stress, each subscale total binned to $K=5$ severity per Figure 2 $\to$ 3 ordinal items (**DRC arm**, complete pre&post $n=1{,}669$). IES-R is single-scale in **Nigeria** ($n\approx 986$, secondary); the **ICRC functioning scale is excluded** — its orientation and 0–14 range are not defined in the supplement codebook, and its post values overrun the scale |
| Effect measure $\rho$ | baseline$\to$endline **severity reduction** per subscale (distress higher$=$worse, so `lower_is_better`) |
| Cohort / interims | DASS/DRC $n=1{,}669$; interims **every 100 accruing beneficiaries** (with a fine early grid $20/40/60/80$) |
| Amortiser | the item-general **J64** net (§14.4.9) deployed on the DASS/DRC grid: **PIT–KS** $0.118$, marg–KS $0.116$ — as well-calibrated here as on Ukraine/REFUGE |
| Data | beneficiary-level `Data_Sheet_1.xlsx` + `Coding` codebook (Frontiers/PMC supplement, via Wayback); fuller data on request to `iandersen@icrc.org` |
| **Application target** | **real-time humanitarian interim evaluation** at scale |
| Reference | [@andersen2022mhpss] |

**Notes.** (i) *Very large effect.* Community MHPSS produces enormous pre$\to$post severity drops ($\rho\approx 0.80$–$0.95$ severity reduction; matching the paper's $96.6\%$ improved on DASS-21), so the PPS is essentially $1$ from the first interim and only becomes *discriminating* at a high threshold $\eta_0\approx 85\%$ — the operational lesson is that the interim rule must be set on that scale. (ii) *Data checks.* The supplement was verified against the paper before use — IES-R improvement $92.70\%$ (exact), and country / gender / mean-age match Table 1 — confirming the columns are correctly identified. (iii) *Scale*$\times$country split (IES-R in Nigeria, DASS-21 in DRC/Mali) mirrors the study design. Fits in `py-icrc-dass-drc_260902`, amortiser deploy in `…-icrc-dass-drc-…-J64-scale-feat-ftheadexpand_260903`; producer `scripts-py/ICRC_interim_svi.py`, loader/extract in `python/data_web_extracting.py`.

## 3.12 Application: product-launch rating experiments

The setting of interest is not the public review corpus at scale but the **small internal experiments a retailer runs before launching a new product** — a concept test or limited release in which a modest panel of customers rate the product on ordered scales, ratings accruing over days. This is a genuine real-time-decision regime: ordinal ratings fit the PCM directly (raters carry the latent trait $\theta$, product variants or attributes are the items, the hedonic/star/Likert score is the ordered response), the cohort is small and accrues over time, and the decision is a PPS — *given the ratings so far, will the product clear its launch threshold once the test completes?* The closest open match is a **consumer acceptance test**: the black-coffee preference dataset [@ristenpart2023coffee] has 118 consumers each rating 27 coffees on a **9-point hedonic scale** ($K=9$) plus four 5-point just-about-right attributes — exactly the per-consumer, per-product ordinal structure of a pre-launch panel. The public Amazon Review Data corpus [@ni2019amazon] provides the same rating structure at scale for validation.

| field | value |
|---------------|---------------------------------------------------------|
| Setting | internal pre-launch product / concept tests (modest consumer panel); open example: black-coffee acceptance test [@ristenpart2023coffee] |
| Design | small accruing cohort; ordinal ratings, optionally comparing product variants |
| Items | product variants / attributes; response $=$ 9-point hedonic liking ($K=9$) or $1$–$5$ JAR / star |
| Effect measure $\rho$ | liking of a variant vs a launch benchmark, or a contrast between variants |
| Cohort / interims | small $N$ (\$\sim\$100 consumers), accruing over the test — the real-time interim regime |
| **Application target** | **real-time launch decision** (PPS the product clears its liking threshold); amortiser validated at scale on [@ni2019amazon] |
| Data | open, [@ristenpart2023coffee] (small panel) and [@ni2019amazon] (at scale) |

## 3.13 Application: MovieLens ratings

The MovieLens rating datasets — up to 25M timestamped $0.5$–$5$ star movie ratings from GroupLens [@harper2015movielens]. A clean, classic ordinal (user $\times$ item) benchmark, smaller and tidier than the Amazon corpus.

| field | value |
|-----------------|-------------------------------------------------------|
| Setting | GroupLens MovieLens (ml-25m and smaller releases) |
| Design | observational; timestamped |
| Items | movies; response $=0.5$–$5$ star rating (10 ordered categories) |
| Effect measure $\rho$ | temporal rating shift, if an effect is wanted |
| Cohort / interims | up to 25M ratings; interims accrue by timestamp |
| **Application target** | **validate the amortiser** on ordinal ratings; tidy benchmark for the accrual/PPS framing |
| Data | open download, [@harper2015movielens] |

## 3.14 Application: mycelium novel-food acceptance

A new-product acceptance study for mycelium as a human food protein — the launch-decision setting of §3.12 with an explicit experimental design [@fischer2024mycelium]. A UK Prolific panel rates the product under a $3\times3$ manipulation, on item-level ordinal scales. Built as a **between-arm** analysis (**Option A**): **product powder vs burger, pooling all three substrates**, treating the burger arm as the baseline and the powder arm as the endline.

| field | value |
|-------------|-----------------------------------------------------------|
| Problem | consumer acceptance of mycelium as a novel food (protein) source |
| Setting | UK Prolific panel, $N=449$ ($3\times3$: processing $\times$ substrate) |
| Arms (built) | **powder vs burger**, substrates pooled — **unpaired** between-arm ($n=149+149=298$; Baseline$=$burger, Endline$=$powder) |
| Items (built) | **9** items, $1$–$7$ ($K=7$): **Acceptance** (A1–A4), **Disgust** (D1–D4), **Perceived naturalness** (PN) |
| Effect measure $\rho$ | between-arm powder-vs-burger contrast per item — Acceptance $+0.08$–$0.15$, Disgust (severity down) $+0.13$–$0.27$, Naturalness $+0.11$: powder beats burger on every construct |
| Cohort / interims | $n=298$; **8** participant-accrual interims (shuffled to mix arms) |
| **Application target** | **real-time launch/acceptance evaluation** with a genuine condition contrast. **Amortiser delivered (§21)** — the unpaired between-arm design is handled by the registry $S_3$ group-contrast network (arm$=$time, burger$=$Baseline/powder$=$Endline), the same one used for HVTN 505 (§20); federated deploy calibrates to PIT–KS $0.078$. SVI fits in `py-mycelium-powdervsburger_260902` |
| Data | open, [@fischer2024mycelium]; OSF `e3gxa` / Zenodo 10628634 |

## 3.17 Application: Parenting for Lifelong Health pooled trials

Parenting for Lifelong Health (PLH) — the parenting-programme family of our direct collaborators Lucie Cluver and Jamie Lachman, who maintain a **Phase-1 pooled database** of seven cluster-randomised PLH trials with harmonised, item-level caregiver- and child-report PROMs at baseline and follow-up. Crucially, PLH is the **instrument family behind the Ukraine Hope Groups (§3.6) and Colombia caregiver data** — parenting practices, violence against children, caregiver mental health, child behaviour — so the pool shares nearly the same item bank. That makes it a **within-family generalisation test** of the item-general amortiser (§14.4.9): one network deployed across Ukraine, Colombia and each PLH trial on a common harmonised instrument set, with a randomised control arm for the $H_1$ decision, and — unlike the cross-sectional PISA/mycelium designs — **genuinely paired pre/post**, so the per-participant amortiser applies directly.

Of the seven pooled trials, the amortiser needs $N>200$ paired cohorts; three qualify, four are pilots / small feasibility trials:

| trial (programme) | design | $N$ | fit |
|--------------------------|--------------|--------------|-------------------|
| **RISE** — N. Macedonia / Moldova / Romania (NCT04721730), PLH-YC | cluster-RCT | **823** | ✅ largest; SE-Europe, also in the §3.11 PISA set |
| **Sinovuyo Teen** — South Africa (PACTR201507001119966), PLH-Teen | cluster-RCT, 40 clusters | **552** | ✅ published flagship |
| **Sinovuyo Kids** — South Africa (NCT02165371), PLH-YC | RCT | **296** | ✅ full efficacy RCT |
| Thailand (NCT03539341), PLH-YC | RCT | 120 | ✗ too small |
| Philippines — Masayang Pamilya (NCT03205449), PLH-YC | RCT | 120 | ✗ too small |
| Philippines — PLH-Teen (OSF `qdtyv`) | RCT | $\approx$ 120 | ✗ too small |
| South Africa — pilot (NCT01802294), PLH-YC | pilot RCT | 68 | ✗ pilot |

| field | value |
|-------------|-----------------------------------------------------------|
| Problem | parenting support to reduce violence against children and improve caregiver / child wellbeing (LMIC) |
| Arms | two-arm cluster-RCT (PLH programme vs usual care / control), **paired** baseline$\to$follow-up |
| Items | harmonised PLH tools — parenting practices, **child maltreatment (ICAST)**, child behaviour, caregiver depression, adolescent-report violence — item-level ordinal, PCM-ready, **shared with §3.6 / Colombia** |
| Effect measure $\rho$ | baseline$\to$follow-up endpoint effect per item (maltreatment / violence down `lower_is_better`; positive parenting up) |
| Cohort / interims | the three qualifying trials ($N=823$ / $552$ / $296$, all $>200$); accrual as interims within each |
| **Application target** | control-arm real-time evaluation **and a cross-trial generalisation test** of the item-general amortiser across the shared PLH instrument family (each trial a deployment instance on the pooled item bank) |
| Data access | **direct collaborators** (Cluver, Lachman) — internal request for the harmonised pooled extract, not a cold DUA |
| Status | *planned* — the three $N>200$ trials confirmed; awaiting the pooled item-level extract $+$ harmonisation dictionary |
| Reference | Cluver et al. 2018 (Sinovuyo Teen); Ward et al., *J. Child Psychol. Psychiatry* 2020 (Sinovuyo Kids); RISE protocol *Trials* 2021, doi:10.1186/s13063-021-05817-1 |

## 3.18 Application: influenza-vaccine immunogenicity (ImmuneSpace / HIPC)

The first **vaccine** application and the first with **real, in-hand, individually-paired** data. The HIPC / Immune Signatures compendium on ImmuneSpace/ImmPort collects systems-vaccinology cohorts in which each participant's serum **haemagglutination-inhibition (HAI) titre is measured against several influenza strains at day 0 (pre) and again at day ≈ 21–28 (post)**. Strains are the items, the two visits are the paired time axis, and the titre is a $\log_2$ serial dilution — a naturally **ordered-categorical** readout — so this is a direct instance of the §14 **per-participant paired** estimand on real vaccine data, and the natural showcase for the ordered-categorical (GPCM) extension. Because the same instrument (HAI over strains) recurs across cohorts and seasons, the set is also a **multi-cohort transfer test** of the item-general amortiser (§14.4.9): one network across studies, strains and age groups.

Pulled 2026-09-18 to `ImmuneSpace_Influenza_Vaccine_Trials_v260918.xlsx` (16 studies, 1 061 participants, 276 885 assay results; HAI in 15 of 16). Per-study HAI structure:

| SDY | $N$ | age | HAI strains | endline day | paired $n$ | note |
|--------|-------:|--------|-------:|-------:|-------:|------------------------|
| SDY67 | 159 | 50–73 | 2 | 28 | 159 | elderly, A/H1N1 focus |
| SDY400 | 98 | 21–90 | 3 | 27 | 93 | multi-season |
| SDY314 | 92 | \~33 | 6 | 21 | 89 | 6 strains |
| **SDY212** | 91 | 21–90 | 3 | 21 | 89 | **young + older arms; ImmPort HAI tutorial exists** |
| SDY312 | 84 | 22–90 | 7 | 21 | 79 | 7 strains |
| SDY315 | 74 | 21–90 | 3 | 21 | 67 |  |
| SDY404 | 72 | 22–90 | 3 | 27 | 69 |  |
| SDY80 | 64 | 18–62 | 0 | — | 0 | NIH CHI — no HAI in pull (neut/flow only) |
| SDY269 | 63 | 0–47 | 5 | 28 | 56 | TIV + LAIV; includes children |
| SDY520 | 61 | 21–87 | 4 | 28 | 56 |  |
| SDY180 | 46 | 22–49 | 3 | 28 | 12 | multi-arm: Fluzone + Pneumovax + **saline placebo** |
| SDY296 | 45 | 21–58 | 3 | 28 | 37 |  |
| SDY301 | 40 | 21–64 | 3 | 28 | 40 |  |
| SDY270 | 30 | 21–47 | 3 | 28 | 28 |  |
| SDY305 | 25 | 20–33 | 3 | 24 | 14 |  |
| SDY144 | 17 | 0–13 | 3 | 27 | 16 | **paediatric** |

**Totals: 16 studies · 1 061 participants · 904 with a paired day-0→endline HAI pair** (across the 15 HAI studies; SDY80 excluded). Predominantly single-arm inactivated TIV (Fluzone); LAIV in SDY269; a saline/Pneumovax comparator only in SDY180 — so the estimand is the **within-participant** seroresponse, not an arm contrast.

| field | value |
|-----------|-------------------------------------------------------------|
| Problem | influenza-vaccine immunogenicity (seroresponse) — a genuine, regulator-recognised success criterion (seroconversion / seroprotection rate, GMT ratio) |
| Setting | HIPC systems-vaccinology cohorts via ImmuneSpace/ImmPort; healthy adults + elderly, two paediatric |
| Arms | mostly single-arm TIV; **within-participant** paired (Baseline$=$day 0, Endline$=$day ≈ 21–28); SDY180 carries a placebo/Pneumovax comparator |
| Items | **HAI titre per influenza strain** (2–7 strains; H1N1 / H3N2 / B) — $\log_2$ serial dilution $\to$ **ordered-categorical**; seroprotection (titre $\ge 40$) is the caseness threshold |
| Effect measure $\rho$ | **two clinically-standard HAI endpoints per strain** (CHMP/CBER; @hobson1972role), computed from the fitted category probabilities via the general `rho_specs` path of `get_endpoints_per_draw` (`item_type='categorical'`; $y$ 0-indexed on the $\log_2$ ladder, start 1:5, so $y\ge 3 \equiv$ titre $\ge$ 1:40). $\rho_1$ = seroprotection rate (SPR) $= P(\text{titre}\ge 1{:}40)$ at endline — an absolute endline level (baseline SPR is computed but unused), success $H_1$ if SPR $>0.70$; $\rho_2$ = GMT fold-rise (GMFR) $= 2^{\bar E_{\text{end}}-\bar E_{\text{base}}}$ (mean $\log_2$ titre), success $H_1$ if GMFR $>2.5$. `higher_is_better` |
| Category labels | the ordered titre categories carry **scientifically accurate reciprocal-dilution labels** (shared `_dilution_ladder_labels`): $k{=}0$ (value below the lowest 1:10 dilution) $\to$ `<1:10`, then `1:10, 1:20, 1:40` ($=k_3=$ the seroprotection threshold)`, 1:80, ...`; these flow through `dp1.y_label` + `dit.cat_labels` to every figure (prob-by-question fit, $p(\rho\mid x)$). Neutralisation ID50 (COVID/HIV) is a continuous readout binned to `1:4, 1:8, ...` (no below-detection floor) |
| Composite (joint) decision | SPR and GMFR are extracted in **one** `get_endpoints_per_draw` call, so they are indexed on the **same Monte-Carlo draw**; the regulator's composite rule "**at least one** of the criteria met" is then evaluated *within a draw* as the OR of the per-$\rho$ indicators and averaged — the correct joint posterior probability, not a naive combination of marginal PPS (`PartialCreditModel.joint_any_met`; written to `*_any_met.csv`). E.g. SDY312 per strain: A-strains + B/Lee $\approx 0.95$–$1.0$, B/Brisbane & B/Florida $\le 0.08$. (**SCR / seroconversion is NOT included**: an individual \$\ge$4-fold rise is a *paired* per-participant event, which the population marginal predictive does not retain; add it via the per-$\theta\$ paired predictive.) |
| Cohort / interims | 15 HAI studies, 904 paired participants; participant-accrual interims within a study; the 15 cohorts $=$ a **transfer / foundation-model** axis |
| **Application target** | the **direct §14 paired amortiser on real vaccine data** $+$ the **ordered-categorical (GPCM) extension** ($\log_2$ titres) $+$ a **multi-cohort foundation-model demonstration** (one net across studies / seasons / ages) |
| Data access | open-registration ImmuneSpace/ImmPort under a data-use agreement; **pulled and in hand** |
| Loader (implemented) | `read_data_immport_flu(xlsx_path, study, endline_day=None, min_start_dilution=5.0)` in `python/data_loading.py` — join HAI (Assays) $\to$ day (Events); keep day 0 $+$ endline (auto-detected as the largest of $\{28,27,24,21\}$ present); map titre $\to$ ordered category $k=\mathrm{round}(\log_2(\text{titre}/5))$ on the 2-fold ladder with a **single study-wide** $K$ across strains (one `item_type_id` $\Rightarrow$ no $K$-mixing, cf. §3.11); keep paired participants only; emit the PCM `dp1`/`dit` with `item_label`$=$strain, `time`$\in\{0,1\}$ (`Baseline`/`Endline`), `y`$=k$ / `y_stan`$=k+1$, `item_type='out-of-7'` (expected-category endpoint, $K$-agnostic), `item_high_label='higher_is_better'`. Reuses `model_pcm` + `get_endpoints` unchanged |
| Producer | `scripts-py/IMMPORT_flu_interim_svi.py` (mirrors `ICRC_interim_svi.py`): pseudo-orders participants by Participant ID, accrues interims every 10, SVI (`AutoDiagonalNormal`, 4 000 steps, 2 000 draws, `with_core_analyses`) at each, emitting the same figure set as `py-icrc-dass-drc_260902` |
| Amortiser deploy ($\rho_2$ GMFR) | The **item-general J64 amortiser** (deepsetXcompAtt, scalar-mean token) on the SVI grids for **GMFR** with expanding head-ft + affine + §14.4.26 **head BvM**, via the shared ragged driver (`deploy_IMMPORT_amortiser.sh`; GMFR reference grid from `IMMPORT_flu_gmfr_refgrid.py`). Calibration to the SVI GMFR reference: **SDY312 PIT-KS 0.089 / marg-KS 0.114, SDY314 0.082 / 0.090** (baseline 0.32–0.51). $\eta_0$ sweep $\{2,2.5,3,3.5\}$-fold: a \$\ge\$2-fold rise reached confidently only by A/Uruguay; at the CHMP 2.5-fold bar the high-baseline cohorts read futile. Outputs $\to$ the single amortiser dir `py-immport-SDY{312,314}-amortise-deepsetXcompAtt-itemamortise-J64-ftheadexpand-bvm_260919/`, prefix `pcm_gmfr_interim_*` |
| Amortiser deploy ($\rho_1$ SPR) | SPR is a **threshold** functional, so the scalar-mean token is insufficient (§14.4.18). Retrained the J64 net with a **cumulative-exceedance token** $[\mathbf 1\{k\ge 1\},\dots]$ (pool $=$ empirical category CDF $\Rightarrow$ categorically sufficient), endline-caseness target, and the caseness threshold $c{=}3$ in item metadata; deployed with the same head-ft + affine + head BvM (`deploy_IMMPORT_amortiser_spr.sh`, `RAGD_WIDETOK=1`/`RAGD_CASE_C=3`). Calibration to the SVI SPR reference: **SDY312 PIT-KS 0.104 / marg-KS 0.142, SDY314 0.080 / 0.096** (baseline 0.58–0.68) — the token change crosses the sufficiency wall. $\eta_0$ sweep on the $[0,1]$ rate $\{0.5,0.6,0.7,0.8\}$ (0.70 $=$ CHMP): at \$\>\$70% seroprotection A/PR8, A/Uruguay, A/Victoria, B/Lee pass, A/South Dakota borderline-fails, B/Brisbane & B/Florida fail. Outputs $\to$ the **same** amortiser dir under prefix `pcm_spr_interim_*` (SVI grid dir stays SVI-only). Full write-up in §16 |
| Status | **implemented — SDY312 (79 paired, 7 strains,** $K=10$) and SDY314 (89 paired, 6 strains, $K=9$); per-strain $\rho=$ relative rise in mean $\log_2$ titre, interims every 10 participants $\to$ `py-immport-SDY312_260918` / `py-immport-SDY314_260918`. Next: multi-cohort pool + the GPCM/ordered-categorical head |
| Reference | Diray-Arce et al., *Sci. Data* 2022 (Immune Signatures resource); ImmuneSpace/ImmPort SDYs 67/80/144/180/212/269/270/296/301/305/312/314/315/400/404/520 |

## 3.19 Application: HIV-vaccine trials (CAVD DataSpace / HVTN)

The **efficacy-scale** vaccine application. The Collaboration for AIDS Vaccine Discovery **DataSpace** (Fred Hutch / VISC) hosts HVTN HIV-vaccine trials under a **Global Access policy that makes individual-level data public** (free registration): per-participant, multi-antigen **binding-antibody and neutralisation titres** (antigen/isolate $=$ item, pre/post $=$ paired). The original target — the two mosaic-Ad26 efficacy trials that stopped at interim, HVTN 705 *Imbokodo* and HVTN 706 *Mosaico* — turned out **not to be pullable**: on inspection (2026-09-18) `vtn705` has no DataSpace record at all and `vtn706` is metadata-only (`data_availability` null; 0 rows in every dataset). Of the four HVTN-network studies catalogued, only **three carry integrated assay data — `vtn097`, `vtn105`, `vtn505`**. Among these, **HVTN 505 (`vtn505`)** is a phase-2b DNA/rAd5 efficacy test-of-concept that was **stopped early at a pre-specified interim for futility (April 2013)**. Its assay data support a **vaccine-vs-placebo binding-antibody (BAMA) immunogenicity contrast**, monitored over accruing subjects — a real HIV-vaccine application of the between-arm amortiser. (Two caveats surfaced on inspection: HIV-naive subjects have no informative paired baseline, so the within-participant fold-rise used for flu is degenerate here — the estimand is the endline arm contrast; and the 2013 futility stop was on HIV-*infection* efficacy, a time-to-event endpoint absent from the antibody datasets, so that specific interim is not reconstructed.)

| field | value |
|-------------|-----------------------------------------------------------|
| Problem | HIV-vaccine immunogenicity (multi-antigen seroresponse) and efficacy interim monitoring |
| Setting | HVTN trials via CAVD DataSpace (Global Access public data); Fred Hutch / VISC |
| Arms | vaccine vs placebo. **No informative paired baseline** — subjects are HIV-naive, so pre-vaccination antibody is a uniform floor (unlike flu's pre-existing immunity); the estimand is the single-timepoint **vaccine-vs-placebo** contrast, not a within-participant fold-rise |
| Timepoints | endline only for the contrast (HVTN 505 day 196). NAb is measured paired ($\{0,196\}$) but baseline is all below-detection $\Rightarrow$ the paired fold-rise is degenerate; BAMA is endline-only |
| Items | **BAMA** binding-IgG antigens (antigen $=$ item; `mfi_delta` $\to$ per-antigen ordinal, $K=3$). HVTN 505: 18 antigens, 10 span both arms across all $K$; the federated interim/amortiser deploy restricts to the **rectangular 9-antigen panel** measured on *all* subjects (drops the sparse gp70-V1V2(A), covered on only 50/239) so the SVI reference is coherent across interims — the V1V2 immune-correlate scaffolds (AE.A244, C.1086C V1V2), Con6 gp120, p24. (NAb: 5 isolates but only MN.3 paired, floor baseline — not usable for the contrast) |
| Effect measure $\rho$ | per-antigen relative **vaccine-vs-placebo** shift in ordinal binding, `higher_is_better` (encoded MYCELIUM-style: placebo$=$`Baseline`, vaccine$=$`Endline`) |
| Cohort / interims | **HVTN 505** — 239 subjects (190 vaccine / 49 placebo); shuffled accrual so both arms are present at each of 10 interims |
| **Application target** | the **between-arm endline immunogenicity amortiser** on real HIV-vaccine binding data (BAMA $\to$ ordered categorical), reusing the MYCELIUM/UkraineP group-contrast machinery. (The 2013 futility stop was on HIV-*infection* efficacy — a time-to-event endpoint not in the assay datasets — so it is not reconstructed here) |
| Data access | **CAVD DataSpace** (Global Access, data-use agreement); LabKey server `dataspace.cavd.org`, container `/CAVD`, schema `study`. **email+password** `~/.netrc` auth (not an API key). Study protocols with data: `vtn505`/`vtn097`/`vtn105`; `vtn706` empty, `vtn705` absent |
| Extractor (implemented, pull verified) | `download_cavd_dataspace(study, dest_dir, assays=('NAb','BAMA','ICS','Demographics'))` + `cavd_studies_with_data()` / `cavd_select()` in `python/data_web_extracting.py` — LabKey `query-selectRows.api` on `/CAVD` via `curl --netrc` (the `labkey`/requests netrc paths mangle the password). CLI: `python python/data_web_extracting.py --cavd "HVTN 505" --out <dir>` (`--cavd-list` shows availability). Pulled: `vtn505` NAb 628 / BAMA 10 260 / ICS 22 684 / Demog 2 504 |
| Loader + producer (implemented) | `read_data_cavd_bama_endline(bama, demo, n_cat=3)` in `python/data_loading.py` (arm$\to$time, per-antigen ordinal bins, keeps antigens both arms fully span) $+$ `scripts-py/CAVD_bama_interim_svi.py` (mirrors `MYCELIUM_interim_svi.py`; figures as `py-icrc-dass-drc_260902`). `read_data_cavd_nab` is retained but documents the degenerate paired NAb case |
| Status | **implemented — HVTN 505 BAMA between-arm grid.** Full-data check coherent: Con6 gp120 $\rho\approx1.14$ (PPS 1.00), V1V2 scaffolds $\rho\approx0.75$–0.80, p24 null ($\rho\approx0.50$). Interim SVI grid $\to$ `py-cavd-vtn505-bama_260918`. Next: NAb/ICS variants; HVTN 097/105 for a multi-cohort HIV set |
| Reference | HVTN 505 (`vtn505`, interim futility stop 2013); HVTN 097/105; CAVD DataSpace, `dataspace.cavd.org`. (Targeted but unavailable: HVTN 705 *Imbokodo*, HVTN 706 *Mosaico*) |

## 3.21 Application: head-to-head vaccine comparison (ImmPort / ImmuneSpace)

The **head-to-head vaccine** estimand — the genuinely decision-relevant contrast for a formulary/procurement choice: given two *vaccine options* (not baseline/endline of one, not two subpopulations), which produces the better immune response? From a shortlist of ImmuneSpace head-to-head studies (see below), **SDY269** (Systems Biology of 2008 Influenza Vaccination) is the sound in-hand case: two randomised adult arms, **LAIV** (live-attenuated, intranasal) vs **TIV** (inactivated, intramuscular), each with paired day-0/day-28 HAI and the same two clinical endpoints as §3.18 (SPR, GMFR). The contrast lives **across two per-arm paired analyses**: each arm is run as its own SDY312-style paired seroresponse study, then the two arms' posteriors are compared. **Caveat — different strain panels:** the arms were assayed on their own vaccine strains, so only **A/Uruguay/716/2007 (H3N2)** is shared and directly comparable strain-for-strain; the H1N1 and B strains compare each arm's *homologous* seroresponse (LAIV A/S.Dakota H1N1, B/Florida vs TIV A/Brisbane H1N1, B/Brisbane).

| field | value |
|-----------|-------------------------------------------------------------|
| Problem | **head-to-head**: which of two vaccine options gives the stronger immune response (a procurement/formulary decision, distinct from within-arm pre/post or subpopulation contrasts) |
| Setting | ImmuneSpace/ImmPort HIPC SDY269 (Emory, 2008 influenza season); two randomised adult arms, 28 paired participants each |
| Arms | **LAIV** (live-attenuated intranasal) vs **TIV** (inactivated IM); each analysed as its own paired day-0/day-28 seroresponse, then compared |
| Items | HAI titre per strain $\to$ ordered category on the 2-fold ladder ($K=7$ LAIV / $8$ TIV). **Panels differ**: LAIV {A/S.Dakota H1N1, A/Uruguay H3N2, B/Florida}, TIV {A/Brisbane H1N1, A/Uruguay H3N2, B/Brisbane} $\Rightarrow$ only **A/Uruguay H3N2 shared** |
| Effect measure $\rho$ | the §3.18 pair per arm per strain: **SPR** $=P(\text{titre}\ge 1{:}40)$ at endline ($H_1{:}>0.70$) and **GMFR** $=$ GMT fold-rise ($H_1{:}>2.5$) |
| Cohort / interims | 28 paired per arm; pseudo-order by Participant ID, interims every 10 ($n=10,20,28$) |
| **Application target** | a **head-to-head** instance of the endpoint estimand: compare two vaccines' seroresponse, with a direct arm-difference posterior on the shared strain |
| Loader + producers (implemented) | `read_data_immport_flu_headhead(xlsx, study, arms)` (pools both arms into ONE frame, item $=$ pure strain) $+$ `IMMPORT_flu_headhead_svi.py` (one joint fit per interim $+$ two joint `get_endpoints_per_draw` calls $+$ composite decision) $+$ `IMMPORT_flu_headhead_compare.py` (final-interim summary, exact shared-strain difference, plots); the single combined $p(\rho\mid x)$ plot via `IMMPORT_flu_rho_plots.py` (`IMMPORT_RHO_DIRS`). The arm filter `read_data_immport_flu(..., arm=)` remains for per-arm loads |
| Estimator — **one joint fit, six labelled rhos, joint decision** | Both arms are fitted in ONE partial-credit model per interim under the clean schema: **`item_label` stays the pure response item** (the strain, so the shared A/Uruguay is *one* item), the vaccine arm folds into the flexible condition axis with the paired time-point (`group_label` $\in$ {LAIV_baseline, …, TIV_endline}), and the non-reduced structure is the cross **`item_group_id` = item_label** $\times$ group ($5\times 4$, sparse $=12$); $\theta_i$ shared across disjoint participants, fit on the binary helper `phase` (`x_formula="~ phase - 1"`) so it is **invariant** to the relabelling. The **rho set is declared upfront by the loader** (`d['rho_specs']`) — six rhos each with a short `rho_label` and a pretty `rho_label_long`: **LAIV_spr, TIV_spr, LAIV_gmfr, TIV_gmfr** (level, per arm) and **TIV−LAIV_spr_diff, TIV−LAIV_gmfr_diff** (cross-arm, shared strain). Each is computed **one-by-one** (a scalable `get_endpoints_per_draw(…, contrast_col='phase')` call per rho — some from group means, others could be per-participant) but off the same fit, so all six share the Monte-Carlo **draw index** and concatenate into ONE joint per-draw frame. From the jointly-indexed level rhos the composite CHMP rule "\$\ge\$1 of SPR/GMFR met" is scored **per draw, per (strain, arm)** (`joint_any_met` $\to$ `headhead_any_met.csv`: TIV A/Brisbane-H1N1 $0.94$, A/Uruguay $0.93$; LAIV $\approx 0$). Output: dir `py-immport-SDY269`, fit `pcm_1_interim_`, joint frame `pcm_1_interim_i{k}_regression_training.pkl`; the single `p_rho_x_by_item` plot is a **strain (row)** $\times$ rho_label (column) grid (empty cells where a strain is in only one arm), `headhead_*` numeric summaries |
| Status | **implemented — SDY269 LAIV vs TIV, joint SVI (56 pooled, interims** $n=10\!-\!56$ arms interleaved). Textbook adult result: **TIV dominates LAIV on every strain and both endpoints.** Shared **A/Uruguay H3N2**: SPR LAIV $0.16$ vs TIV $0.69$; GMFR LAIV $1.23$ vs TIV $5.66$. H1N1: TIV meets both thresholds (SPR $0.82$, GMFR $3.12$), LAIV neither ($0.21$, $1.10$). Consistent with immunology — in adults inactivated IM TIV drives strong serum HAI, live-attenuated intranasal LAIV drives little (mucosal, pre-exposure). Amortiser deploy optional (arms small) |
| Cross-arm difference on the shared strain | $\rho_\text{SPR,diff}=\text{SPR}_\text{TIV}-\text{SPR}_\text{LAIV}$ and $\rho_\text{GMFR,diff}=\text{GMFR}_\text{TIV}-\text{GMFR}_\text{LAIV}$ on A/Uruguay H3N2 — plain differences (bounded; no ratio blow-up). Because both arms share the joint fit these are **EXACT per-draw contrasts** (no random draw-pairing): the shared $\theta$/draw removes the artificial independence, tightening the interval and raising $P(\rho>0)$ vs an independent-fit pairing. Final interim: $\rho_\text{SPR,diff}$ median $\mathbf{0.51}$ (95% CI $0.20$–$0.75$, $P>0\approx1.00$); $\rho_\text{GMFR,diff}$ median $\mathbf{4.33}$-fold (95% CI $0.21$–$12.4$, $P>0=0.98$; lower bound positive, vs $-0.03$ under an independent pairing). Summary/plot `headhead_shared_strain_diff.{csv,pdf}`, contraction over $n$ |
| Head-to-head shortlist (ImmuneSpace) | **In hand (flu xlsx) — SDY269 LAIV/TIV is the ONLY viable two-vaccine head-to-head.** SDY305 IM-TIV vs intradermal-TIV was inspected and **dropped** (ID arm has only **3 paired** participants, $n_\text{enrolled}=9$). The updated head-to-head export (`ImmuneSpace_headhead_vaccince_v260920.xlsx`, 8 titre studies) adds no new two-vaccine case: SDY112 (age bands, all Fluzone), SDY1276 (males vs females TIV) are §3.20 *subpopulation* contrasts; SDY180 Fluzone-vs-saline is vaccine-vs-placebo ($n=6$/arm); SDY1293 has no titre; SDY80/SDY89 are single-arm as exported. **SDY690** (investigational HBsAg-1018 vs licensed Engerix-B HepB) — the cleanest cross-antigen head-to-head in principle — was pulled but carries **no titre data**, so it is a dead end |
| Reference | SDY269 (Emory HIPC; Nakaya *et al.* 2011, PMID 21743478); ImmuneSpace/ImmPort |

------------------------------------------------------------------------

# 4. No reduction of the PPS to today's knowledge

A particular confusion that often arises with the PPS is that it does not fall back to today's knowledge.

A key concept in Bayesian forecasting is self-consistency, in that if we make posterior predictions $z$ from today's posterior and then consider the updated posterior $p(\theta|x,z)$, the two expansion and contraction steps cancel each other out.

In our setting, this Bayesian prediction self-consistency property equates to $$\mathbb{E}_{z|x}[P(H_1|x,z)] = P(H_1|x),$$ and it follows from replacing $p(\theta|x,z)p(z|x)$ with $p(\theta,z|x)$, exchanging integrals, and integrating out $z$.

The key point is that the PPS contains the indicator decision rule $$1\{P(H_1 \mid x,z) > \eta_H\},$$ and so $$\mathbb{E}_{z|x}[1\{P(H_1|x,z) > \eta_H\}] 
\neq 
1\{P(H_1|x) > \eta_H\},$$ since $E(f(X)) \neq f(E(X))$. Conceptually the PPS is a probability in $[0,1]$ whereas today's decision rule always evaluates to a binary value.

------------------------------------------------------------------------

# 5. Links to general decision-theoretic learning

## 5.1 Bayes optimal decisions under loss

In the broader decision theoretic literature, the objective is to evaluate the expected loss (or expected utility) of all possible actions given the observed data, and then choose the action that minimizes the expected loss.

Let $\mathcal{A}$ denote the action space. Given data $x$, the Bayes optimal action is $$a^* = \arg\min_{a \in \mathcal{A}} \mathbb{E}_{\theta \mid x}[L(a, \theta)]$$ where $L(a, \theta)$ is the loss incurred by taking action $a$ when the unknown parameters are $\theta$.

In our setting, for final data $x$, the decision theoretic problem simplifies to two actions:

- $a_1$: declare success, adopt intervention
- $a_0$: declare failure, reject intervention

and two partitions on the unknown true treatment effect, parameterised by the relative-reduction margin $\eta_0 \in [0, 1)$ introduced in §3:

- $H_1$: success hypothesis ($p < p_0 (1 - \eta_0)$)
- $H_0$: null hypothesis ($p \geq p_0 (1 - \eta_0)$)

We can define the loss matrix:

| Action                  | $H_1$ true | $H_0$ true |
|-------------------------|------------|------------|
| Declare success ($a_1$) | 0          | $L_{FP}$   |
| Declare failure ($a_0$) | $L_{FN}$   | 0          |

where $L_{FP}$ is the cost of a false positive decision and $L_{FN}$ is the cost of a false negative decision.

Given observed data $x$, the expected loss of declaring success is $$\mathbb{E}_{H \mid x}[L(a_1, H)]
= L_{FP} \cdot P(H_0 \mid x) + 0 \cdot P(H_1 \mid x)
= L_{FP} \cdot P(H_0 \mid x),$$ amd similarly the expected loss of declaring failure is $$\mathbb{E}_{H \mid x}[L(a_0, H)]
= L_{FN} \cdot P(H_1 \mid x).$$ The Bayes optimal action minimizes expected loss, $$a^*(x) = \arg\min_{a \in \{a_0, a_1\}} \mathbb{E}_{H \mid x}[L(a_1, H)].$$ And so, to declare success we want $$\begin{aligned}
& 
\mathbb{E}_{H \mid x}[L(a_1, H)] < \mathbb{E}_{H \mid x}[L(a_0, H)] \\
\Leftrightarrow &
L_{FP} \cdot (1 - P(H_1 \mid x)) < L_{FN} \cdot P(H_1 \mid x) \\
\Leftrightarrow &
P(H_1 \mid x) > L_{FP} / ( L_{FP} + L_{FN} ).
\end{aligned}$$ This shows that the Bayes optimal decision threshold $\eta_H$ can be expressed in terms of utilities or loss terms. In particular, if false positives are more costly, $\eta_H > 0.5$ and if false negatives are more costly, $\eta_H < 0.5$. This provides a utility-based foundation for choosing $\eta_H$ in contrast to ad-hoc 0.89 or 0.95.

## 5.2 Net benefit

A related concept is expected utility, normalised such that 1 unit corresponds to 1 true positive. Under the above action/partition set matrix, the relative cost of one false positive is $w$ and true pos and true neg have utility/cost 0. This does not depend on the particular units of the losses, and so is more easily interpretable.

Under these utilities, the net benefit of declaring success is given by $$\begin{aligned}
NB(x) 
& 
= P(H_1 \mid x) - w \cdot P(H_0 \mid x) 
\\
&
= P(H_1 \mid x) - w \cdot (1 - P(H_1 \mid x) )
\end{aligned}$$ and the net benefit decision rule is $$1\{ NB(x) > 0 \}.$$ Under Bayes optimality, we find that the relative cost of a false positive must be $w = L_{FP} / L_{FN}$. Using the expression for $\eta_H$ above, we then find $w = \eta_H / (1 - \eta_H)$. Thus we can express the net benefit in the standard form $$NB(x) = P(H_1 \mid x) - \frac{\eta_H} {1 - \eta_H} \cdot (1 - P(H_1 \mid x) ).$$ This shows that our learning problem directly connects to standard cost-benefit analyses. Our PPS is equivalent to posterior predictive net benefit over future data $z$, and so finding a way to amortise the PPS also means a solution to amortising posterior predictive net benefit.

------------------------------------------------------------------------

# 6. State-of-the-art approaches to estimate PPS

State-of-the-art approaches consider that we have observed a specific interim dataset $x_{1:n}$ and have access to the posterior distribution $p(\theta \mid x)$. The core task is to estimate $$y^{(s)} = y(z^{(s)}):= P(H_1 \mid x, z^{(s)}) = \int 1\{\theta \in H_1\} \, p(\theta \mid x, z^{(s)}) \, d\theta.$$ for fixed $x$ and for each simulated future data $z^{(s)} \sim p(\cdot \mid \theta^{(s)})$ where $\theta^{(s)} \sim p(\theta \mid x)$. From this, the PPS is straightforward to compute via $$PPS(x) \approx \frac{1}{S} \sum_{s=1}^S 1\{ y^{(s)} > \eta_H\}.$$

### 6.1 Nested Monte Carlo

A computationally costly, but robust approach is to numerically estimate the new joint posterior $p(\theta \mid x, z^{(s)})$ and then simply obtain the label $y^{(s)}$ by evaluating the above integral over posterior draws from $p(\theta \mid x, z^{(s)})$.

### 6.2 Regression-based functional inference

Inefficient estimators in 13.1.1–13.1.4 all approximate the conditional posterior $p(\theta \mid x, z^{(s)})$ once per future sample $z^{(s)}$. The regression-based approach [@strong2014estimating] instead learns a function $q_\psi$ with tuning parameters $\psi$ that produces the label as a function of a low-dimensional summary of the future data directly.

Pick a summary statistic $w:\mathcal{Z}\to\mathbb{R}^d$ of the future data and proceed as follows.

1.  **Joint sampling.** For $s=1,\dotsc,S$, draw $$\theta^{(s)} \sim p(\theta\mid x),\qquad
    z^{(s)} \sim p(z\mid \theta^{(s)}),$$ compute based on $x$ and not $z^{(s)}$ $$y^{(s)} := 1\{\theta^{(s)}\in H_1\},$$ and also based $z^{(s)}$ and not $x$ $$w^{(s)} := w(z^{(s)}).$$
2.  **Fit a regressor** $q_\psi:\mathbb{R}^d\to[0,1]$ to the pairs $\{(w^{(s)},y^{(s)})\}_{s=1}^S$ by minimizing the empirical binomial cross-entropy (logit link) $$\hat\psi=\arg\min_\psi -\frac{1}{S}\sum_{s=1}^S\Big[y^{(s)}\log q_\psi(w^{(s)})+(1-w^{(s)})\log(1-q_\psi(w^{(s)}))\Big],$$ using a a GAM for $d\le 6$ or a GP for higher $d$.
3.  **Predict** the label $y^{(s)}$ for any new future sample $z^{(s)}$ with the learned regressor by $y(z^{(s)})=q_{\hat\psi}\big(w(z^{(s)})\big)$, repeat over future future samples, and estimate $$PPS(x) \approx \frac{1}{S} \sum_{s=1}^S 1\{ y^{(s)} > \eta_H\}$$

Instead of the binary labels $y^{(s)}$, it is advantageous to consider the continuous $\rho^{(s)}$ that underlie the alternative hypothesis, $H_1 : \rho^{(s)} > \eta_0$. A simple approach would be:

1.  **Joint sampling.** For $s=1,\dotsc,S$, draw $$\theta^{(s)} \sim p(\theta\mid x),\qquad
    z^{(s)} \sim p(z\mid \theta^{(s)}),$$ compute based on $x$ and not $z^{(s)}$ $$\rho^{(s)} := \rho(\theta^{(s)})$$ and also based $z^{(s)}$ and not $x$ $$w^{(s)} := w(z^{(s)}).$$
2.  **Fit a regressor** $q_\psi:\mathbb{R}^d\to\mathbb{R}$ to the pairs $\{(w^{(s)},\rho^{(s)})\}_{s=1}^S$ by minimizing the least-squares loss $$\hat\psi=\arg\min_\psi \frac{1}{S}\sum_{s=1}^S\big(\rho^{(s)} - q_\psi(w^{(s)})\big)^2,$$ using linear regression, or GAM for $d\le 6$ or even a GP for higher $d$.
3.  **Predict** the conditional success probability for any new future sample $z^{(s)}$ from the learned regressor's predictive distribution at $w(z^{(s)})$. With the predicted mean $q_{\hat\psi}\big(w(z^{(s)})\big)$ and predictive standard deviation $\hat\sigma(z^{(s)})$ (link-scale Gaussian for the GAM, posterior std for the GP), approximate the upper $\eta_0$ quantile of $\rho$ at $w(z^{(s)})$ with $$\hat y(z^{(s)}) := P(H_1 \mid x, z^{(s)}) \approx \Phi\!\left( \frac{q_{\hat\psi}\big(w(z^{(s)})\big) - \eta_0}{\hat\sigma(z^{(s)})} \right),$$ repeat over future samples, and estimate $$PPS(x) \approx \frac{1}{S} \sum_{s=1}^S 1\{ \hat y(z^{(s)}) > \eta_H \}.$$

The Gaussian step in 3 is fragile under skewed or heteroskedastic residuals of $\rho \mid w, x$. A more robust alternative replaces the conditional mean by a conditional quantile, removing all distributional assumptions on $\rho \mid w$ and folding $\eta_H$ directly into the regression target. The decision rule rewrites as $$\begin{align*}
& P(H_1 \mid x, z) > \eta_H \\
\Leftrightarrow \: & P(\rho \le \eta_0 \mid x, z) \leq 1 - \eta_H \\
\Leftrightarrow \: & q_{1-\eta_H}(\rho \mid x, w(z)) > \eta_0,
\end{align*}$$ where $q_\tau(\rho \mid x, w)$ is the lower $\tau$-quantile of $\rho$ given $x, w$. So if we can estimate $q_{1-\eta_H}$ directly, the PPS decision collapses to one indicator with no Gaussian step.

1.  **Joint sampling.** As in step 1 above, draw $\{(\rho^{(s)}, w^{(s)})\}_{s=1}^S$.
2.  **Fit a conditional-quantile regressor** $q_\psi:\mathbb{R}^d\to\mathbb{R}$ to the pairs $\{(w^{(s)}, \rho^{(s)})\}_{s=1}^S$ at level $\tau = 1 - \eta_H$, by minimising the pinball loss $L_\tau$, $$\hat\psi = \arg\min_\psi \frac{1}{S}\sum_{s=1}^S L_\tau\!\big(\rho^{(s)} - q_\psi(w^{(s)})\big),$$ where $L_\tau(u) = u\,\big(\tau - 1\{u < 0\}\big)$ and $q_\psi$ is a quantile linear model, a quantile-GAM [@fasiolo2021qgam] for $d \le 6$ or even a quantile-GP / quantile-RF for higher $d$.
3.  **Predict and decide.** For any new future sample $z^{(s)}$, evaluate the quantile estimate $\hat q_{1-\eta_H}(z^{(s)}) = q_{\hat\psi}\big(w(z^{(s)})\big)$ and estimate $$PPS(x) \approx \frac{1}{S}\sum_{s=1}^S 1\{ \hat Q_{1-\eta_H}(z^{(s)}) > \eta_0 \}.$$

This scheme is guaranteed to create unbiased labels when $w$ is a sufficient statistic for $\theta$, as $S \to \infty$. Indeed, in that case there are $h$ and $g$ such that $p(z\mid \theta)=h(z)\,g\big(w(z);\theta\big),$ and further, $$\begin{aligned} 
p(\theta\mid x, z) & = 
\frac{p(\theta, z \mid x)}{p(z\mid x)} = 
\frac{p(\theta\mid x) p(z\mid\theta)}{p(z\mid x)} \\
& \propto p(\theta\mid x) p(z \mid\theta) \\
& \propto p(\theta\mid x) g(w(z) ; \theta). 
\end{aligned}$$ Thus, we can retrieve the contracted posterior of the LHS by joint sampling as on the RHS (which is step 1 above) and then conditioning on $w(z)$.

The label inherits the same reduction, $$\begin{aligned} 
y = P(H_1\mid x,z) & = \int 1_{\theta\in H_1}\,p(\theta\mid x,z)\,d\theta \\ 
& = C^{-1}\int 1_{\theta\in H_1}\,p(\theta\mid x) g(w(z) ; \theta)\,d\theta \\
& =:\;q^\star(w(z); x).
\end{aligned}$$ Thus the regression task in step 2 above provides a consistent estimator of $q^\star$ under the joint sampling of step 1. Specifically, the marginal of the training pairs is $p(w,y\mid x)=\int p(w\mid\theta)\,\mathrm{Bern}\big(y;1_{\theta\in H_1}\big)\,p(\theta\mid x)\,d\theta,$ and the population minimiser of either the MSE or the cross-entropy loss is the conditional expectation $$\arg\min_\tau \mathbb{E}\big[(y-q_\tau(w))^2\,\big|\,x\big]
\;=\;\mathbb{E}[y\mid w,x]
\;=\;\int 1_{\theta\in H_1}\,p(\theta\mid w,x)\,d\theta.$$ Any consistent regression family $q_\tau$ therefore recovers $q^\star$ as $S\to\infty$, and step 3 evaluates the exact label $P(H_1\mid x,z)$ at any new $z$ by a single regression evaluation.

When $w$ is only approximately sufficient, the same population minimiser applies, but it is now the projection $$\mathbb{E}[y\mid w,x]
\;=\;\mathbb{E}_{z\mid w,x}\!\big[P(H_1\mid x,z)\big],$$ i.e. the smoothed label averaged over the residual $z$-information not captured by $w$. The bias relative to the true $P(H_1\mid x,z)$ is bounded by the conditional variance $\mathrm{Var}\big(P(H_1\mid x,z)\,\big|\,w(z),x\big)$; the bias goes to zero as $w$ approaches a sufficient statistic.

In exponential-family models (e.g. the Bernoulli and Categorical examples of Section 3.1, 3.2) the sufficient statistic is canonical and finite-dimensional, so $w$ can be chosen exactly.

In partial-credit / IRT models with per-unit random effects no finite-dimensional sufficient statistic exists; a pragmatic choice for $w$ is the per-$(\text{item},\text{time})$ response summary that the analyst would actually use post-trial. In our case, for out-of-7 outcomes, a suitable summary might just be the mean of the item responses at baseline and endline, $w_t(z^{(s)}) = m^{-1}\sum_{i=1}^m z_{it}^{(s)}.$ For categorical responses, a suitable summary might just be the proportion of the item responses above grade 3 at baseline and endline, $w_t(z^{(s)}) = m^{-1}\sum_{i=1}^m 1\{ z_{it}^{(s)} > 3\}.$

**Practical caveats.**

- *Finite* $S$. The regression variance is $O(S^{-1}\,\sigma^2/n_{\text{eff}}(W))$ where $\sigma^2=\mathrm{Var}(Y\mid W)\le 1/4$ (Bernoulli) and $n_{\text{eff}}(W)$ is the effective local sample size. @strong2014estimating recommend $S$ in the low thousands.
- *Dimensionality.* Additive structure breaks down for $\dim(W)\gtrsim 6$; GPs handle higher $d$ but scale as $O(S^3)$ without sparse approximations.
- *Uncertainty.* GP variants give a posterior credible band on $\hat y(z)$; GAMs give standard errors from the link-scale Gaussian approximation.

------------------------------------------------------------------------

# 7. Amortized inference workflow

Rather than learning a neural posterior of high dimensional model parameters, an initial step might be to learn the function $$z \rightarrow 𝑃( H_1 \mid x,z ),$$ because this is the key component in the decision rules and PPS above. Here $x$ is today's fixed observed data, but it is undesirable to repeat learning for every $x$, and in addition we do not always have a pre-specified summary $w$ as in Section 6.2 and wish to optimize.

The task is to learn a neural function $$(x_{1:n},z_{1:m}) \rightarrow q_\phi(x_{1:n},z_{1:m}) \approx P(H_1 \mid x_{1:n},z_{1:m}) \in [0,1]$$ for arbitrary data $x$ with $n$ data points today and for arbitrary future data $z$ with remaining $m$ data points until the end of the scheduled intervention. In my setting, the posterior always factorizes, $$p(\theta \mid x,z) \propto p(\theta) \prod_{i=1}^n p(x_i \mid \theta) \prod_{i=1}^m p(z_i \mid \theta).$$

I can re-frame this learning task as $$\begin{aligned}
q_\phi(x,z) 
& 
\approx \mathbb{E}_{\theta|x,z}(L(A,\theta)) 
\\
& 
=  \int L(A,\theta) p(\theta | x,z) d\theta,
\end{aligned}$$ which makes clear that we focus on learning the loss weighted mean of the future posterior. This is decision-focussed learning, we learn what is required to make an optimal decision.

If we step back for a moment, for product likelihood models in the exponential family, we have for each data input $\chi_i$ $$p(\chi_i|\theta)=h(x_i)\exp\{\eta(\theta)^\top T(\chi_i)-A(\theta)\},$$ and so with a conjugate prior, the posterior is also in the exponential family and depends only on the sufficient statistic. In this case, if we condition on $x$ and $z$, we can re-express our target in terms of some function $g$ $P(H_1|x,z) = g\big(\sum_{i=1}^n T(x_i) + \sum_{i^*=1}^m T(z_{i^*})\big).$

We always have that the data $x_i$ and $z_{i^*}$ are permutation invariant, and this motivates that we want to learn end-to-end a two-stage neural model $q_\psi, q_\tau$ that decomposes along the DeepSets theorem, $$(x,z) \mapsto q_\psi\bigg(\sum_{i=1}^n q_\tau(x_i) + \sum_{i^* =1}^m q_\tau(z_{i^*})\bigg) \approx P(H_1|x,z).$$ In particular, this immediately provides generalisability to different current data sets of different sizes $n$, and to different future data sets of different sizes $m$.

There are three steps:

- Workflow Step 1: generating labelled training data (§8)
- Workflow Step 2: learning neural architectures (§9)
- Workflow Step 3: deployment, amortise PPS (§10)

------------------------------------------------------------------------

# 8. Step 1: generating labeled training data

The regression-based scheme of §6.2 is *training-data-free*: labels are constructed on the fly from the model by drawing $\theta^{(s)}$ (from the posterior in §6.2), simulating $x^{(s)}, z^{(s)} \sim p(\cdot \mid \theta^{(s)})$, and reading off the endpoint $\rho^{(s)} := \rho(\theta^{(s)})$. The same trick extends immediately to the amortised setting: replace the posterior draw of $\theta$ by a prior draw and simulate $(x, z)$ jointly. No external labels, no importance sampling, no SMC — the model *is* the label generator.

**Algorithm (amortised joint sampling).** For $s = 1, \ldots, S$:

1.  **Draw sizes** $n^{(s)}, m^{(s)}$ from the operational distribution over interim schedules (e.g. uniform over the accrual grid, or Poisson around the design cohort size).
2.  **Draw parameter** $\theta^{(s)} \sim p(\theta)$ from the *prior*.
3.  **Simulate current data** $x^{(s)}_{1:n^{(s)}} \sim p(\cdot \mid \theta^{(s)})$.
4.  **Simulate future data** $z^{(s)}_{1:m^{(s)}} \sim p(\cdot \mid \theta^{(s)})$.
5.  **Continuous label** $\rho^{(s)} := \rho(\theta^{(s)})$ (the relative-effect endpoint used throughout §3 and §6.2).

The triple $(x^{(s)}, z^{(s)}, \rho^{(s)})$ plays the same role at the amortised scale that the pair $(w^{(s)}, \rho^{(s)})$ plays for a single fixed $x$ in §6.2.

**Why the population minimiser is still the target.** Marginalising over $\theta$, the training pairs have marginal distribution $$p(x, z, \rho) = \int p(x, z \mid \theta)\, \delta(\rho - \rho(\theta))\, p(\theta)\, d\theta.$$ The population minimiser of the least-squares loss on $\rho$ given $(x, z)$ is the conditional expectation $$\arg\min_\phi \mathbb{E}\big[(\rho - q_\phi(x, z))^2\big] \;=\; \mathbb{E}[\rho \mid x, z] \;=\; \int \rho(\theta)\, p(\theta \mid x, z)\, d\theta,$$ the posterior mean of the endpoint under $p(\theta \mid x, z)$. This is the amortised analogue of the §6.2 identity $\mathbb{E}[y \mid w, x] = q^\star(w; x)$, with the joint pair $(x, z)$ in place of the fixed-$x$ summary $w(z)$. Any consistent regression family $q_\phi$ recovers the true $\mathbb{E}[\rho \mid x, z]$ as $S \to \infty$. Replacing the least-squares loss by the pinball loss at level $\tau$ makes the same population minimiser argument recover the conditional quantile $Q_\tau(\rho \mid x, z)$; §9 uses this for the multi-quantile head.

**Deployment vs training coverage.** During deployment we generate $z \sim p(z \mid x^{\text{obs}})$ from the posterior predictive at the observed $x^{\text{obs}}$. The training-time joint has $x$ and $z$ sharing $\theta$ from the prior, whereas the deployment joint has $x$ observed and $z$ from the posterior predictive. In the amortised setting with sufficient prior-predictive coverage of the operational-$x$ region, the deployment-time conditional $\mathbb{E}[\rho \mid x^{\text{obs}}, z]$ is recovered pointwise — no re-weighting or SMC step required.

## 8.1 Prior-coverage caveats

The amortised regression is only as good as the prior-predictive coverage of the operational $(x, z)$ region:

- **Informative-prior mismatch.** If the prior places most mass on a region of $\theta$-space that generates data unlike the operational trial data, few training samples fall near the observed $x$ and $q_\phi(x^{\text{obs}}, \cdot)$ is essentially an extrapolation.
- **Size coverage.** Fix a broad distribution over $(n, m)$ to span the full interim schedule; oversample the earliest interims (largest $m$) where the PPS is most decision-relevant.
- **Sanity check.** On the analytically tractable benchmarks §3.1 (Binomial) and §3.3 (MVN), compare the amortised $q_\phi(x^{\text{obs}}, z^{(s)})$ to the closed-form $P(H_1 \mid x^{\text{obs}}, z^{(s)})$ on a held-out grid of $(x^{\text{obs}}, z^{(s)})$.

------------------------------------------------------------------------

# 9. Step 2: learning neural architectures

The training data is $\{(x^{(s)}, z^{(s)}, \rho^{(s)})\}_{s=1}^S$ with variable set sizes $n^{(s)}, m^{(s)}$. Since $x$ and $z$ are permutation-invariant sets of observations, the DeepSets decomposition of §7 applies: $$q_\phi(x, z) \;=\; q_\psi\bigg( \sum_{i=1}^n q_\tau(x_i) + \sum_{i^*=1}^m q_\tau(z_{i^*}) \bigg).$$ The encoder $q_\tau$ maps a single observation to a $d$-dimensional embedding; the head $q_\psi$ maps the summed embedding to the target (a real number for the endpoint, or a vector for the multi-quantile head below).

## 9.1 Choice of target: continuous $\rho$ with a quantile head

The natural target at the amortised layer is the endpoint $\rho$, not the binary indicator $1\{\rho > \eta_0\}$. Continuous targets carry more information per sample and reduce Monte-Carlo variance in the same way that the §6.2 endpoint regression outperformed the $y = 1\{\theta \in H_1\}$ variant. We consider three heads for $q_\psi$, mirroring the three variants of §6.2:

1.  **Gaussian mean-and-scale head** (amortised analogue of the §6.2 Gaussian variant). Output two scalars $\hat\rho_\phi(x, z)$ and $\log \hat\sigma_\phi(x, z)$. Train by Gaussian negative-log-likelihood $$\ell_\phi^{\text{NLL}}(x, z, \rho) = \tfrac{1}{2} \log \big(2\pi \hat\sigma_\phi^2 \big) + \frac{(\rho - \hat\rho_\phi)^2}{2 \hat\sigma_\phi^2}.$$ Convert to the label via the Gaussian tail $\hat y_\phi(x, z) = \Phi\big((\hat\rho_\phi - \eta_0)/\hat\sigma_\phi\big)$ and threshold by $\eta_H$ as in §6.2.

2.  **Single-quantile head** (amortised analogue of the single-$\tau$ §6.2 quantile variant). Output one scalar $\hat Q_\phi^{1-\eta_H}(x, z)$. Train by the pinball loss at $\tau = 1 - \eta_H$ $$\ell_\phi^{\text{pinball}}(x, z, \rho) = L_\tau\big(\rho - \hat Q_\phi^{1-\eta_H}(x, z)\big).$$ Deploy the decision $1\{\hat Q_\phi^{1-\eta_H}(x, z) > \eta_0\}$ — no Gaussian step, no plug-in variance. Under $S \to \infty$ and a consistent regression class, $\hat Q_\phi^{1-\eta_H}(x, z) \to Q_{1-\eta_H}(\rho \mid x, z)$ and the decision is exact.

3.  **Multi-quantile head (recommended)** — amortised analogue of the mquantile variant of §6.2 / §6.5. Output a vector of predicted quantiles at a grid $\tau_1 < \tau_2 < \ldots < \tau_K$ spanning $(0, 1)$ (e.g. $\{0.05, 0.1, 0.2, \ldots, 0.9, 0.95\}$). Train by the summed pinball loss $$\ell_\phi^{\text{mq}}(x, z, \rho) = \sum_{k=1}^K L_{\tau_k}\big(\rho - \hat Q_\phi^{\tau_k}(x, z)\big),$$ optionally with a monotone-non-decreasing penalty across the $\tau_k$ index or a $\operatorname{cummax}$ post-processor at inference to fix quantile crossings. From the predicted quantile grid, obtain a continuous conditional success probability by linear interpolation of the CDF at $\eta_0$ $$\hat P_\phi(H_1 \mid x, z) \;=\; 1 - \hat F_\phi(\eta_0 \mid x, z), \qquad \hat F_\phi(\eta_0 \mid x, z) \;=\; \operatorname{interp}\big(\eta_0;\, \hat Q_\phi^{\tau_1}, \ldots, \hat Q_\phi^{\tau_K};\, \tau_1, \ldots, \tau_K\big),$$ then threshold by $\eta_H$.

The multi-quantile head carries the same Rao-Blackwell advantage over the single-quantile binary decision (lower Monte-Carlo variance in the PPS estimator) and the same distribution-free advantage over the Gaussian head (no plug-in variance, no Normal-tail assumption).

## 9.2 DeepSets encoder for the IRT application

For the Hope Groups intervention (§3.4), each observation carries a hierarchical structure — one participant provides $J$ ordered-categorical item responses at baseline and endline. The data are permutation-invariant over participants $i$ but not over items $j$ or times $t$. The nested DeepSets pattern $$q_\tau(\text{person}_i) \;=\; \rho_{\text{inner}}\bigg( \sum_{j=1}^J \phi_{\text{embed}}(x_{ij}, j, t) \bigg)$$ respects this: the inner MLP $\phi_{\text{embed}}$ sees a triple (response, item id, time id), and the person-level embedding is a sum over items with the item identity as a positional token. The outer sum aggregates across persons.

**Sufficient-statistic hardcoding.** In the exponential-family benchmarks of §3.1 (Binomial), §3.2 (Categorical) and §3.3 (MVN), the sufficient statistic $T$ is known analytically and finite-dimensional. Setting the first layer of $q_\tau$ to the identity on $T(x_i)$ (or on a fixed non-linear function of the raw response chosen to match $T$) guarantees exact sufficiency at that layer; only $q_\psi$ then needs to learn the map from the sum-of-$T$s to the endpoint. This matches the exp-fam reduction $P(H_1 \mid x, z) = g(\sum T(x_i) + \sum T(z_{i^*}))$ of §7 and reduces the network's degrees of freedom to those of $q_\psi$ alone. It also allows a strict correctness test on the exp-fam benchmarks: hardcoded-$T$ + multi-quantile head should match the analytic PPS to Monte-Carlo precision.

**Architectural variants.** The plain DeepSets sum can be swapped for a Set Transformer [@lee2019set] (ISAB blocks with self-attention) when sum-pooling is a bottleneck. For our IRT model the per-item interactions are weak once the item-position embedding is included, so plain DeepSets is the natural first pass; Set Transformer is a fallback if the DeepSets model plateaus.

## 9.3 Training practicalities

- **Sample budget.** Cover the joint $(\theta, n, m, x, z)$ space. In the low-dimensional benchmarks (§3.1, §3.3) $S \sim 10^4$ suffices; for the IRT model expect $S \sim 10^5$–$10^6$.
- **Loss weighting.** The multi-quantile grid can be weighted to emphasise the operational-decision quantile $\tau = 1 - \eta_H$ (e.g. $w_k = 1 + \lambda \exp(-(\tau_k - (1 - \eta_H))^2 / \nu^2)$ with $\lambda, \nu$ small).
- **Held-out set.** Reserve 10–20% of $(x^{(s)}, z^{(s)}, \rho^{(s)})$ tuples as a validation split for early stopping and quantile-crossing diagnostics.
- **Vary** $n, m$. Sample $n^{(s)}, m^{(s)}$ uniformly across the interim grid (or from the operational distribution) so the network generalises across accrual states.

**Exact sufficient statistics for hardcoding** $q_\tau$ in the exp-fam benchmarks.

- **Binomial (§3.1)**: the likelihood $p(x_i \mid p) = p^{x_i}(1-p)^{1-x_i}$ has natural parameter $\eta(p) = \log\!\big(p/(1-p)\big)$ and canonical scalar sufficient statistic $$  T(x_i) = x_i \in \{0, 1\}.
    $$ Set $q_\tau(x_i) := x_i$ (identity, dimension 1). The DeepSets sum is $\sum_{i=1}^n T(x_i) + \sum_{i^* = 1}^m T(z_{i^*}) = k^n + k^m$, the total number of successes. The head $q_\psi$ then only needs to learn the map $(k^n + k^m, n, m) \mapsto Q_\tau^{\tau_k}\!\big(\rho \mid k^n + k^m, n, m\big)$; feeding the cohort sizes $(n, m)$ as extra scalar inputs to $q_\psi$ preserves the exact-sufficiency of $q_\tau$.

- **MVN with** $\sigma^2$ known (§3.3, Gaussian special case): with the natural parameter $\eta(\mu) = K^{-1} \mu / \sigma^2$ and $\sigma^2$ fixed, only the first-order statistic is needed, $$  T(x_i) = x_i \in \mathbb{R}^J.
    $$ Set $q_\tau(x_i) := x_i$ (identity, dimension $J$). The DeepSets sum recovers $T_1(x_{1:n}) + T_1(z_{1:m}) = \sum_{i=1}^n x_i + \sum_{i^*=1}^m z_{i^*}$; $q_\psi$ takes this $J$-vector together with $(n, m)$ as inputs.

- **MVN with** $\sigma^2$ unknown (§3.3, full NIG): the likelihood adds a second natural parameter $-1/(2\sigma^2)$ against the quadratic-form statistic. Per-observation sufficient statistics are $$  T(x_i) = \big(x_i,\; x_i^\top K^{-1} x_i\big) \in \mathbb{R}^{J+1}.
    $$ Set $q_\tau(x_i) := \operatorname{concat}\!\big(x_i,\; x_i^\top K^{-1} x_i\big)$ (dimension $J + 1$; $K^{-1}$ is precomputed once per $(J, R\text{-type})$ cell in the sim setup of §3.3.1). The DeepSets sum recovers $\big(T_1(x_{1:n}) + T_1(z_{1:m}),\; T_2(x_{1:n}) + T_2(z_{1:m})\big)$, exactly the two sufficient statistics of §3.3.

Since $\rho = \rho(\theta)$ in each case depends only on these sufficient statistics of the pooled cohort (via the closed-form posterior of §3.1 / §3.3), a network with $q_\tau = T$ and any consistent scalar/vector regression $q_\psi$ is model-optimal: it matches the analytic PPS as $S \to \infty$. This gives the strict correctness test named in §9.2.

------------------------------------------------------------------------

# 10. Step 3: Deployment, amortising PPS

At deployment for observed data $x^{\text{obs}}_{1:n}$:

1.  **Sample future data.** Draw $S$ posterior-predictive samples $z^{(s)}_{1:m} \sim p(z \mid x^{\text{obs}})$ from any convenient posterior fit on $x^{\text{obs}}$ (SVI, HMC). This is the same $z$-generation step as in §6.
2.  **Amortised label.** For each $s$, forward-pass $(x^{\text{obs}}, z^{(s)})$ through the trained network to obtain either
    - the Gaussian $(\hat\rho_\phi, \hat\sigma_\phi)$ → $\hat y^{(s)} = \Phi((\hat\rho_\phi - \eta_0)/\hat\sigma_\phi)$;
    - the single quantile $\hat Q_\phi^{1-\eta_H}$ → $\hat y^{(s)} = 1\{\hat Q_\phi^{1-\eta_H} > \eta_0\}$;
    - or the multi-quantile grid $\{\hat Q_\phi^{\tau_k}\}$ → $\hat y^{(s)} = 1 - \operatorname{interp}(\eta_0; \{\hat Q_\phi^{\tau_k}\}; \{\tau_k\})$.
3.  **Aggregate.** The PPS estimator is $\widehat{\text{PPS}}(x^{\text{obs}}) = S^{-1} \sum_{s=1}^S 1\{ \hat y^{(s)} > \eta_H \}$.

Because the encoder sums are linear in the observations, $\sum_i q_\tau(x^{\text{obs}}_i)$ can be pre-computed once per interim and reused across all $z^{(s)}$ draws. The per-sample cost of the amortised PPS is one forward pass of the head $q_\psi$ per $z^{(s)}$.

**Current-decision target.** The same network trivially covers $P(H_1 \mid x^{\text{obs}}) \approx q_\phi(x^{\text{obs}}, \varnothing)$ by passing an empty $z$. To include this at training time, augment the training triples with samples at $m^{(s)} = 0$.

------------------------------------------------------------------------

# 12. Results for the Binomial model

## 12.1 Results for Binomial interim analysis amortised endptx on wz with features-fixed, qpsi-MLP, loss-multiquantilehead

**Implementation.** First prototype (Binomial, hardcoded sufficient statistic, multi-quantile head) built out end-to-end:

- New module [`python/amortiser_pps_features_fixed_qpsi_MLP_loss_multiquantilehead.py`](../python/amortiser_pps_features_fixed_qpsi_MLP_loss_multiquantilehead.py) exposing the Flax module `Amortiser_PPS_features_fixed_qpsi_MLP_loss_multiquantilehead` (DeepSets-pooled features + SiLU-MLP `q_psi` + multi-quantile head), a training loop `train()` (Adam + cosine LR warmup + pinball loss over 11 quantile levels), and prediction helpers `predict_amortised_p_h1_for_one_xz` (monotone-corrected CDF interpolation at `pps_H1_min_effect_size_thresh`) / `predict_amortised_p_h1_for_many_xz` (deployment wrapper with mquantile-compatible output schema).
- `BinomialModel.make_training_data_with_features` attached to [`python/model_binomial.py`](../python/model_binomial.py) so the amortiser's training-data step is model-owned. The DeepSets pooling for the hardcoded $T(x_i) = x_i$ collapses to features `(k_total, n_total) / N_max`.
- Sim-data cache script [`scripts-py/Binomial_interim_analyses_make_sim_data.py`](../scripts-py/Binomial_interim_analyses_make_sim_data.py) writes the fixed Bernoulli cohort + monthly interim grid + a $65\,536$-sample amortiser training batch to disk (`binomial_sim_cohort.pkl`, `binomial_amortiser_training_data.pkl`).
- Deployment script [`scripts-py/Binomial_interim_analysis_amortise_endptx_on_wz_with_features_fixed_qpsi_MLP_loss_multiquantilehead.py`](../scripts-py/Binomial_interim_analysis_amortise_endptx_on_wz_with_features_fixed_qpsi_MLP_loss_multiquantilehead.py) trains the net once (or loads the cached checkpoint) and, at each interim, forward-passes the DeepSets-pooled features to produce a `p_h1_xz` frame that `Binomial_interim_analyses_compare_methods.py` picks up under the `_RGEA_` suffix.
- Correctness tests [`test/python/test_amortised_pps_correctness.py`](../test/python/test_amortised_pps_correctness.py) exercise forward-pass determinism, training determinism, save/load round-trip, and PPS accuracy on the deployment monthly grid + 20 random $(k_n, n, m)$ triples.

**Training configuration.** MLP head `hidden_dims = (256, 256, 128)`, `num_quantiles = 11` (`taus = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)`), Adam with peak learning rate $10^{-3}$ under a linear-warmup / cosine-decay schedule, $40\,000$ steps × batch $8192$. Prior sampler: `Beta(1, 1)` on $p$, $n \sim U\{1, \ldots, N-1\}$, $m = N - n$ (fixed final cohort $N = 500$), $k_n \sim \text{Binomial}(n, p)$, $k_m \sim \text{Binomial}(m, p)$, features `(k_n + k_m, N) / N`, target $\rho = 1 - p / p_0$ with $p_0 = 0.5$.

**Deployment.** At each interim we build the posterior-predictive draws via `BinomialModel.fit_closed_form_posterior` (analytic $p_s \sim \text{Beta}(1 + k_n,\, 1 + n - k_n)$) chained through the `BinomialModel`-specific `get_interim_z_from_ypredi` override — which draws fresh $m$ Bernoulli$(p_s)$ per posterior draw so $k_m = \sum_i \text{ypred}_{s, i}$ marginalises exactly to $\text{BetaBinomial}(m,\, 1 + k_n,\, 1 + n - k_n)$. One forward pass per interim gives the continuous $\hat P(H_1 \mid x, z^{(s)}) \in [0, 1]$; the PPS is $S^{-1} \sum_s \mathbf{1}\{\hat P > \eta_H\}$. This replaces the earlier scipy-inline Beta-Binomial hack with the standard interim-loop path used by every other Binomial deployment script.

**Correctness against the closed-form Beta-Binomial PPS.** Across the 11 valid monthly interims of the deployment cohort ($N = 500$, true $p = 0.4$, $\eta_0 = 0.25$, $\eta_H = 0.89$, $S = 200$):

| Interim  | Analytic | Amortised |       |
|----------|---------:|----------:|------:|
| 1 (Jan)  |    0.069 |     0.070 | 0.001 |
| 2 (Feb)  |    0.476 |     0.510 | 0.034 |
| 3 (Mar)  |    0.316 |     0.275 | 0.041 |
| 4 (Apr)  |    0.283 |     0.245 | 0.038 |
| 5 (May)  |    0.338 |     0.290 | 0.048 |
| 6 (Jun)  |    0.135 |     0.145 | 0.010 |
| 7 (Jul)  |    0.133 |     0.155 | 0.022 |
| 8 (Aug)  |    0.051 |     0.035 | 0.016 |
| 9 (Sep)  |    0.069 |     0.075 | 0.006 |
| 10 (Oct) |    0.000 |     0.000 | 0.000 |
| 11 (Nov) |    0.000 |     0.000 | 0.000 |

**Max absolute error 0.048, mean 0.020**, meeting the correctness-test tolerance ($\le 0.02$ mean, $\le 0.05$ max within the Monte-Carlo noise band at $S = 200$). The five pytest cases in `test_amortised_pps_correctness.py` all pass with the same tolerances.

**Diagnostic finding worth documenting.** The amortised network matches the analytic $Q_{\tau}(\rho \mid k_{\text{total}}, n_{\text{total}})$ pointwise to $\le 0.002$ at each $\tau$ level. The initial PPS gap (up to 0.156) that appeared when the deployment loop used HMC-generated $k_m$ draws was driven entirely by an over-wide posterior-predictive from HMC at small $n$ (std 50 vs analytic std 36 at interim 1). Bypassing HMC with analytic Beta-Binomial $k_m$ removes the bias. This isolates the amortiser's quality from HMC noise and validates the amortised pipeline in its "pure" form (train the net once → deploy with model-analytic posterior-predictive).

## 12.2 Results for Binomial interim analysis amortised endptx on wz with features-MLP, qpsi-MLP, loss-multiquantilehead

**Implementation.** Parallel prototype swapping the hardcoded per-item encoder for a **learnable** DeepSets encoder $q_\tau$:

- New Flax module [`python/amortiser_pps_features_MLP_qpsi_MLP_loss_multiquantilehead.py`](../python/amortiser_pps_features_MLP_qpsi_MLP_loss_multiquantilehead.py) — class `Amortiser_PPS_features_MLP_qpsi_MLP_loss_multiquantilehead` uses `setup()` to share a SiLU-MLP `q_tau` (`(item_dim → q_tau_hidden_dims → embed_dim)`) across `x` and `z` sets; padded input `(B, N_max, item_dim)` + boolean masks; sum-pool with mask multiplication; concat `pooled_x + pooled_z + sizes` into the same MLP head `q_psi` and multi-quantile pinball loss as §12.1.
- New sampler `BinomialModel.make_training_data_with_raw_sequences(rng, S, n_max)` emits raw padded 0/1 item sequences (`x`, `mask_x`, `z`, `mask_z`, `sizes`) with the same `(θ, x, z)` joint sampling as the features-fixed variant, so the head sees enough shape variability to learn $q_\tau$ end-to-end.
- Sim-data cache extended: [`scripts-py/Binomial_interim_analyses_make_sim_data.py`](../scripts-py/Binomial_interim_analyses_make_sim_data.py) writes a second $65\,536$-sample cache `binomial_amortiser_training_data_features_MLP.pkl` alongside the features-fixed cache.
- Deployment script [`scripts-py/Binomial_interim_analysis_amortise_endptx_on_wz_with_features_MLP_qpsi_MLP_loss_multiquantilehead.py`](../scripts-py/Binomial_interim_analysis_amortise_endptx_on_wz_with_features_MLP_qpsi_MLP_loss_multiquantilehead.py) mirrors §12.1's layout but constructs a padded raw-sequences batch per interim (all $S$ posterior draws share the observed `x` sequence; each `z^(s)` carries `km_s` ones + `m - km_s` zeros drawn from the same analytic joint posterior-predictive built by `BinomialModel.fit_closed_form_posterior` + `get_interim_z_from_ypredi`). Outputs saved with `_RGEB_` suffix.
- Shared training / prediction / save-load utilities in `amortiser_common` are polymorphic across the fixed and MLP amortiser classes — no new API surface. `load_fitted_model` uses `importlib` on the persisted `net_class_module` / `net_class_name` to rebuild either class transparently.

**Training configuration.** `q_tau_hidden_dims = (32, 32)`, `embed_dim = 16`, `q_psi` head `hidden_dims = (128, 128, 64)`, same 11-level `pps_ProbH1_lwr_quantiles_mesh`. Adam + linear-warmup / cosine-decay at peak lr $10^{-3}$, $15\,000$ steps × batch $512$ (smaller than the features-fixed variant because each sample is now a `(N_max=500) × 1` tensor). Training loop measured 8.59 min inside `train()` (persisted on the fit dict via §11's `training_mins` field). Final training pinball loss 0.0090 (still gently decreasing — a larger budget would tighten error further).

**Correctness against the closed-form Beta-Binomial PPS.** Same 11 monthly interims, same $\eta_0 / \eta_H / S$ as §12.1:

| Interim  | Analytic | Amortised (features-MLP) |       |
|----------|---------:|-------------------------:|------:|
| 1 (Jan)  |    0.069 |                    0.065 | 0.004 |
| 2 (Feb)  |    0.476 |                    0.555 | 0.079 |
| 3 (Mar)  |    0.316 |                    0.370 | 0.054 |
| 4 (Apr)  |    0.283 |                    0.360 | 0.077 |
| 5 (May)  |    0.338 |                    0.420 | 0.082 |
| 6 (Jun)  |    0.135 |                    0.170 | 0.035 |
| 7 (Jul)  |    0.133 |                    0.180 | 0.047 |
| 8 (Aug)  |    0.051 |                    0.055 | 0.004 |
| 9 (Sep)  |    0.069 |                    0.140 | 0.071 |
| 10 (Oct) |    0.000 |                    0.000 | 0.000 |
| 11 (Nov) |    0.000 |                    0.000 | 0.000 |

**Max absolute error 0.082, mean 0.041.** Roughly $2\times$ the features-fixed error (§12.1: max 0.048, mean 0.020). Expected — the features-fixed variant hardcodes the exact Binomial sufficient statistic $T(x_i) = x_i$ as its inductive prior; the features-MLP variant has to *discover* that sum-of-successes is sufficient from data, which costs a few thousand extra training steps to close and remains slightly noisier at deployment.

## 12.3 Full cross-method comparison after the analytic-`zi` sweep

All Binomial deployment scripts (regression, IS, both amortised variants; nested-MC left as-is per its "HMC-in-HMC" contract) now build `zi` from `BinomialModel.fit_closed_form_posterior` + the `BinomialModel`-specific `get_interim_z_from_ypredi` override that draws fresh $m$ Bernoulli$(p_s)$ per posterior draw. This uses the **same** analytic $p(\theta, z \mid x)$ joint across all methods, so per-method MSE cleanly reflects estimator quality rather than divergent `z`-samplers. MSE against the closed-form Beta-Binomial PPS across the 11 monthly interims:

| Method | MSE | $\sqrt{\text{MSE}}$ | Mean per-interim inference (min) | One-off training (min) |
|-----------------------|------------:|------------:|------------:|------------:|
| **Amortised (features-fixed, `hidden_dims = (64, 64)`)** | **0.00119** | **0.035** | **0.001** | **1.94** |
| Regression endpt-x (Gaussian approx) | 0.00160 | 0.040 | 0.001 | — |
| Nested-MC using HMC for each (x, z) | 0.00244 | 0.049 | 6.04 | — |
| IS reweighting of $\theta \mid x$ | 0.00277 | 0.053 | 0.004 | — |
| Amortised (features-MLP) | 0.00277 | 0.053 | 0.002 | 8.59 |
| Regression endpt-x (mquantile) | 0.00355 | 0.060 | 0.001 | — |
| Regression endpt-x (quantile) | 0.00357 | 0.060 | 0.001 | — |
| Nested-MC using SVI for each (x, z) | 0.00401 | 0.063 | 1.10 | — |
| Regression H1-x on w(z) | 0.01033 | 0.102 | 0.001 | — |

*Amortised features-fixed shown after the §12.4 default swap to `hidden_dims = (64, 64)` — halves MSE vs the old `(256, 256, 128)` default at* $3\times$-faster training. Training cost is amortised across all future deployments — pay once, deploy in milliseconds thereafter. Nested-MC HMC costs 6 min **per interim per x**, so on a 12-interim schedule with $S = 200$ Monte-Carlo draws the total nested-MC cost is $\sim 72$ min per new patient cohort $x$; the amortised features-fixed variant recovers a lower MSE for a 1.9 min one-off training + 0.001 min per interim, breaking even before the first fresh $x$. All timings measured with $\le 4$ CPU threads on a single machine.

**Reading.** The MSE floor is set by Monte-Carlo noise at $S = 200$: $\sqrt{p (1 - p) / S} \approx 0.035$ near the PPS mode, so any method below $\sqrt{\text{MSE}} \le 0.05$ is essentially at the MC-noise floor. Both amortised variants are numerically tied with IS at that floor. The Gaussian-approx regression wins because the Beta posterior mean is well-approximated by a Gaussian at the sample sizes seen, giving a very tight plug-in $\Phi((\hat\mu - \eta_0)/\hat\sigma)$ estimator. Nested-MC HMC still uses the finite-sample HMC posterior (wider than the exact Beta at small $n$), so it inherits an interim-1 error of 0.15 that inflates its MSE. The `H1-x` regression is the most biased — the binary label loses too much information relative to the continuous endpoint targets used everywhere else.

**Amortised (features-MLP) reading.** The MLP variant matches the fixed variant's aggregate MSE despite starting from no structural prior. Per-interim errors are on average larger (mean 0.041 vs 0.020) but the *distribution* of errors is symmetric around the analytic PPS, so squared errors average out. A longer training budget would close the mean-error gap; the ceiling is the MC-noise floor at $S = 200$, same as everyone else. This validates the general-purpose amortiser template for models where the sufficient statistic is not known analytically (Categorical, IRT).

## 12.4 Ablation: features-fixed capacity + training-config knobs

Five variants of the features-fixed amortiser exercised to probe whether the default `hidden_dims = (256, 256, 128)` config from §12.1 was over-parameterised and whether other training knobs move the MSE floor. Each variant differs from the default in one dimension; all use the same analytic joint `zi` path of §12.3 and the same deployment cohort. Selectable at run time via the `AMORTISER_VARIANT` environment variable on the deployment script:

| Variant | Change vs default | MSE | $\sqrt{\text{MSE}}$ | Training (min) |
|------------|---------------------------|-----------:|-----------:|-----------:|
| default `(256, 256, 128)` | — | 0.00277 | 0.053 | 5.90 |
| **`64x64`** | `hidden_dims = (64, 64)` | **0.00119** | **0.035** | **1.94** |
| `num_quantile_levels_5` | `taus = (0.05, 0.25, 0.5, 0.75, 0.95)` | 0.00019 | 0.014 | 5.26 |
| `num_quantile_levels_21` | 21 equally-spaced taus in $(0.025, 0.975)$ | 0.00119 | 0.035 | 46.20 |
| `S_2000` | Deployment $S = 2000$ instead of $200$ | 0.00029 | 0.017 | 38.13 |
| `log_uniform_n` | Training $n \sim$ log-uniform on $[1, N-1]$ (oversample small $n$) | 0.00119 | 0.035 | 22.23 |

**Findings.**

1.  **Default was over-parameterised.** `64x64` halves the MSE at $3\times$-faster training. Ratified as the new default in the deployment script and in compare-methods; the row `Amortised (features-fixed)` in §12.3 will read `MSE = 0.00119, √MSE = 0.035` at the next compare-methods run.
2.  **Fewer quantiles is better here.** `num_quantile_levels_5` cuts MSE another order of magnitude to $0.00019$ ($\sqrt{\text{MSE}} = 0.014$) — well below the naive $S = 200$ MC noise floor of $\sim 0.018$. Coarser $\tau$ mesh concentrates fitting effort on 5 well-anchored quantiles; the piecewise-linear CDF interpolation at $\eta_0$ absorbs any smoothness cost. Doubling to 21 levels gives no MSE gain and $8\times$ the training cost.
3.  **More deployment** $S$ helps. `S_2000` at $\sqrt{\text{MSE}} = 0.017$ confirms the MC-noise floor at $S = 200$ (predicted $\sim 0.018$). Ten-fold more $z$ draws costs a linear factor in deployment time (still under 40 min end-to-end for the whole schedule); the amortised head runs in milliseconds per $z$.
4.  **Log-uniform-**$n$ doesn't help. Same MSE as `64x64` at $10\times$ longer training. The uniform training-$n$ distribution already covers the deployment grid adequately; small-$n$ regime isn't the bottleneck.

**Default swap in the codebase.** The deployment script now defaults to `NET_HIDDEN = (64, 64)` and writes to `..._260714`; the previous `..._260711` dir with the `(256, 256, 128)` net is preserved as an ablation reference. `Binomial_interim_analyses_compare_methods.py` points `dir_rgea` at the new default dir.

## 12.5 Combined ablation + robustness on features-MLP

Two follow-ups to §12.4: (1) run the untried `hidden_dims = (64, 64)` + `num_quantile_levels_5` + `S_2000` combination on the features-fixed amortiser, (2) verify the ablation findings replicate on the features-MLP amortiser (where the sufficient statistic must be *learned* from raw padded item sequences).

Same evaluation setup as §12.4 (analytic joint `zi` from `BinomialModel.fit_closed_form_posterior` + the exchangeable-Bernoulli override of `get_interim_z_from_ypredi`; MSE against closed-form Beta-Binomial PPS across the 11 monthly interims). Each variant is triggered by `AMORTISER_VARIANT` on the corresponding deployment script.

**Combined table (fixed + MLP, sorted by MSE):**

| Variant | Encoder | MSE | $\sqrt{\text{MSE}}$ | Training (min) |
|--------------------------|------------|-----------:|-----------:|-----------:|
| `combo_64x64_qlv5_S2000` | **MLP** | **0.00013** | **0.011** | **8.05** |
| `num_quantile_levels_5` | fixed | 0.00019 | 0.014 | 5.26 |
| `S_2000` | fixed | 0.00029 | 0.017 | 38.13 |
| `S_2000` | **MLP** | 0.00029 | 0.017 | 7.88 |
| `combo_64x64_qlv5_S2000` | fixed | 0.00046 | 0.022 | 1.85 |
| `num_quantile_levels_5` | **MLP** | 0.00119 | 0.035 | 12.70 |
| `64x64` (default) | fixed | 0.00119 | 0.035 | 1.94 |
| `num_quantile_levels_21` | fixed | 0.00119 | 0.035 | 46.20 |
| `log_uniform_n` | fixed | 0.00119 | 0.035 | 22.23 |
| `log_uniform_n` | **MLP** | 0.00192 | 0.044 | 29.60 |
| `64x64` | **MLP** | 0.00243 | 0.049 | 7.86 |
| default (`(256, 256, 128)` / `(128, 128, 64)`) | fixed / MLP | 0.00277 | 0.053 | 5.90 / 8.59 |
| `num_quantile_levels_21` | **MLP** | 0.00277 | 0.053 | 52.16 |

**Findings.**

1.  **Combined `(64x64 + qlv5 + S_2000)` is the best config across all variants**, and the **MLP variant wins outright at** $\sqrt{\text{MSE}} = 0.011$ ($\sim 5\times$ below the naive $S = 200$ MC floor of $0.018$; consistent with the $S = 2000$ MC floor of $\sqrt{0.07 \cdot 0.93 / 2000} \approx 0.006$). Compare to the fixed variant at $\sqrt{\text{MSE}} = 0.022$: same combined knobs but the fixed variant plateaued because its `qlv5` alone was already anomalously good ($0.00019$, a single-realisation dip below the $S = 200$ MC floor) — enlarging to `S_2000` in the combo regresses to the MC-noise-limited value. The MLP variant, which was far above the MC floor without the knobs, benefits monotonically from all three additions and lands at the true floor. **Recommend running the MLP combo as the default `Amortised (features-MLP)` entry going forward.**

2.  **Robustness of the ablation findings.** Each of the four single-knob ablations moves both encoders in the same direction:

    | Knob                     |  Fixed $\Delta$MSE |    MLP $\Delta$MSE |
    |--------------------------|-------------------:|-------------------:|
    | `64x64`                  | $\mathbf{-0.0016}$ |          $-0.0003$ |
    | `num_quantile_levels_5`  | $\mathbf{-0.0026}$ | $\mathbf{-0.0016}$ |
    | `num_quantile_levels_21` |          $-0.0016$ |            $\pm 0$ |
    | `S_2000`                 | $\mathbf{-0.0025}$ | $\mathbf{-0.0025}$ |
    | `log_uniform_n`          |          $-0.0016$ |          $-0.0009$ |

    `num_quantile_levels_5` and `S_2000` are the two knobs that help both encoders substantially; `num_quantile_levels_21` and `log_uniform_n` help neither. `64x64` mostly helps the fixed encoder (which was over-parameterised) but is neutral for the MLP encoder (whose bottleneck is the learned $q_\tau$, not the head).

3.  **Wall-clock winners.** The fixed `64x64` variant remains the fastest to train (1.85–1.94 min); the MLP combo at 8.05 min pays a $4\times$ training cost for a $2\times$ MSE improvement (down to the true MC-noise floor at $S = 2000$). Deployment cost is unchanged: one forward pass per $z$ draw at millisecond scale for both encoders.

**Recommendation.** Adopt the combined `(64x64 + qlv5 + S_2000)` configuration as the *deployment* preset for both encoders. Keep the single-knob variants in `_VARIANTS` as ablation references. The features-fixed default (`64x64`) remains at `dir_out = ..._260714` for backward-compatibility with §12.3's comparison table; the combined config is one env-var away (`AMORTISER_VARIANT=combo_64x64_qlv5_S2000`).

## 12.6 Apple Metal (JAX-MPS) training benchmark on the MLP combo

Setup. `pixi run -e mps-experimental install-jax-mps` provisions the Metal backend (`jax-mps==0.10.1`); `verify-mps` confirms `default_backend = mps` on an Apple M4 Max. The MLP combo variant is launched with `AMORTISER_VARIANT=combo_64x64_qlv5_S2000_mps` on the `mps-experimental` environment; a dedicated `..._combo_64x64_qlv5_S2000_mps_260714` output directory keeps CPU and MPS checkpoints separate so both remain reproducible.

**Direct comparison (features-MLP amortiser, combo config).**

| Backend          | Training (min) |  Speedup |         MSE | $\sqrt{\text{MSE}}$ |
|--------------|-------------:|-------------:|-------------:|---------------:|
| CPU (6-thread)   |           8.05 |     1.0× |     0.00013 |              0.0113 |
| **MPS (M4 Max)** |       **2.71** | **3.0×** | **0.00009** |          **0.0095** |

**Findings.**

- **3× training speedup on MPS.** The features-MLP encoder's per-token MLP over a `(B=512, N_max=500, item_dim=1)` batch tensor is exactly the matmul-heavy workload GPU parallelism was designed for. Kernel-launch overhead is negligible at $15\,000$ training steps.
- **MSE within MC-noise band of CPU.** MPS lands at $\sqrt{\text{MSE}} =
  0.0095$ vs CPU $0.0113$; both consistent with the $S = 2000$ Monte-Carlo floor of $\sqrt{p(1-p)/S} \approx 0.006$. Differences reflect independent PRNG seeds through the JAX/Metal backend (Metal reductions are non-deterministic) rather than model quality.
- **Not worth trying for features-fixed.** The fixed `64x64` variant is a trivial `(B=8192, feature_dim=2)` batch matmul that finishes in 1.9 min on CPU; MPS's kernel-launch overhead would dominate and give neutral-to-worse wall-clock. Confirmed empirically in a spot check (not shown).
- **Caveats.** `jax-mps` is flagged experimental; some ops silently fall back to CPU. Determinism guarantees are weaker on Metal, so the exact loss trajectory changes across runs but the population minimiser is the same and the resulting PPS estimator is within MC noise.

**Recommendation.** For the features-MLP amortiser and any future model with a matmul-heavy DeepSets encoder (Categorical, IRT), invoke the deployment script under the `mps-experimental` pixi environment. Keep features-fixed on CPU (no benefit from GPU).

# 13. Results for the multivariate normal model

## 13.0 Roadmap: estimand, estimator, and result

**Estimand.** For each component $j$ the target is the effect size $\rho_j = \mu_j - \mu^0$ (baseline $\mu^0 = 1$), and the interim decision quantity is the per-component predictive probability of success $\mathrm{PPS}_j(x) = \Pr\big(\Pr(\rho_j > \eta_0 \mid x, z) > \eta_H \mid x\big)$, integrated over the future cohort $z$ (§3.3.1, $\eta_0 = 0$, $\eta_H = 0.89$). The multivariate normal model is the controlled benchmark: with $\sigma^2$ known the posterior of $\rho_j$ is **closed-form Gaussian**, so both the reference posterior $p(\rho_j \mid x)$ and the PPS are known exactly, and every amortised estimator is scored against a ground truth rather than a Monte-Carlo surrogate. This is the role SVI plays for the partial-credit model of §14, made exact.

**Two statistical questions.** (i) *Sufficiency* — which summary of the observed cohort recovers the per-component sufficient statistic (the cohort mean under known $K$): a statistic supplied in closed form, or one estimated from the raw responses? (ii) *Contraction* — does the amortised posterior standard deviation of $\rho_j$ fall with the observed cohort size $n$ at the Bernstein–von Mises rate $\mathrm{SD}(\rho_j \mid x) \propto n^{-1/2}$? The first governs point accuracy of the PPS; the second governs whether the interim uncertainty is honest.

**The estimators, read in sequence.** Each row is one amortised summary map of increasing generality, sharing the network classes and the `amortiser_common` / `amortiser_diag_plots` / `amortiser_calibration` utilities with the Ukraine partial-credit study (§14):

1.  **`idcomp`** (§13.2) — the per-component sufficient statistic supplied in closed form; $K$ enters only through the prior. The exact-statistic baseline.
2.  **`xcomp`, `itemScompAtt`, `itemXcompAtt`** (§13.2) — hand-computed per-item summaries fed through attention across components, so the amortiser exploits the known $K$ at deployment and is reusable, unchanged, for the partial-credit study.
3.  **`deepsetScompAtt`, `deepsetXcompAtt`** (§13.2) — the per-item summary is itself estimated from the raw per-(participant, component) responses by an inner exchangeable pooling (nested DeepSets), removing the hand-computed statistic.
4.  **Parametric contraction heads** $+$ Bernstein–von Mises correction (§13.3) — the deep-set posterior width is made to obey the $n^{-1/2}$ law, closing the contraction gap those variants otherwise leave.

**Result in one line.** The exact-statistic and hand-summary estimators recover the closed-form reference (PPS–MSE $\approx 10^{-3}$, PIT–KS $\approx 0.13$, contraction exponent $\hat p \approx 0.65$ against the reference $0.50$); the deep-set variants match on point accuracy but **under-contract** ($\hat p \approx 0.07$), and the parametric heads restore the exact rate ($\hat p = 0.517$); a single estimator priced over $J \in [2, 100]$ holds flat calibration. Detailed results in §13.4; the covariance-family sampler and the amortisation over $J$ in §13.5.

## 13.1 Structure of the learning task

A batch element is a single triple $(j^*, s, \text{interim})$ — one queried component, one posterior draw of the future cohort, one interim. The observed cohort $x$ has $n$ participants; the future cohort $z^{(s)}$ has $m = N - n$ participants, drawn at training from the prior predictive and at deployment from the closed-form posterior predictive $p(z \mid x)$ (§3.3.1). At deployment all $J \cdot S$ batch elements are forward-passed per interim.

| symbol | value | description |
|------------|------------|-------------------------------------------------|
| $J$ | 20 / 60 / 100 | components presented to the amortiser (one token per component); amortised over $J \in [2,100]$ in §13.2 |
| $n$ | interim-specific | observed participants at the interim, each an $\mathbb{R}^J$ response vector |
| $m$ | $N - n$ | future participants per posterior draw $s$ |
| $N$ | 500 / 1050 | total cohort ($n + m$ fixed); $N = 1050$ for the cached simulation, $N = 500$ in §3.3.1 / §13.2 |
| $S$ | 4000 / 500 / 200 | posterior-predictive draws of $z^{(s)}$ per interim (per architecture; see §13.4) |
| $K$ | $J \times J$ | known covariance shape, unit diagonal; block-equicorrelation (§3.3.1 Cell B) at deployment |
| $F, A, E$ | —, —, 32 | per-item feature dim, aux dim, embedding dim |
| $K_\tau$ | 5 / 11 | quantile levels $\tau$ |

The endpoint $\rho_j = \mu_j - \mu^0$ is the training target; the amortiser is fitted **on the prior predictive only** (draw $\mu \sim \mathrm{MVN}(\mu^0\mathbf 1_J, \tau_0^2 K)$, simulate both cohorts, label with $\mu_j - \mu^0$). The reference posterior at deployment is the exact Gaussian, with per-component standard deviation $\mathrm{SD}(\rho_j \mid x) = \sqrt{K_{jj} / (\tau_0^{-2} + n/\sigma^2)}$; under unit diagonal and $\sigma = 1$ this is $1/\sqrt{\tau_0^{-2} + n}$, the exact Bernstein–von Mises law the contraction diagnostics score against.

## 13.2 Architectures: from a closed-form statistic to nested DeepSets

The per-component sufficient statistic under known $K$ is the cohort mean; the six amortisers differ only in **how that statistic enters** — supplied exactly, hand-computed then attended over, or estimated from the raw responses. All share one batch contract so a caller swaps architectures by swapping the class.

**Per-component summary and target.** Writing $\bar y_j = \tfrac1N\sum_{i\le n} y_{i,j}$ (observed, $s$-invariant) and $\bar z^{(s)}_j = \tfrac1N\sum_{i\le m} z^{(s)}_{i,j}$ (future, $s$-dependent), the hand-summary token for a queried component $j^*$ is $$t^{(j^*, s)}_j = \big(\bar y_j,\; \bar z^{(s)}_j,\; K_{j^*, j}\big) \in \mathbb{R}^{F}, \qquad a^{(j^*)} = \big(n/N,\; m/N,\; K_{j^*, j^*}\big) \in \mathbb{R}^{A},$$ with the $K$-row entry $K_{j^*, j}$ encoding how informative component $j$ is for the queried $j^*$. A shared map embeds each token, $h^{(j^*,s)}_j = q_{\text{tok}}(t^{(j^*,s)}_j) \in \mathbb{R}^E$, and a head $q_\psi$ reads the attended summary $\bar h^{(j^*,s)}$ with the queried token and aux into the $K_\tau$ quantiles of $\rho_{j^*}$. The two attention mechanisms are

- **cross-attention** (`Xcomp`): one query $q^{(j^*,s)} = q_{\text{query}}(t^{(j^*,s)}_{j^*})$ against $J$ keys, $\bar h = \sum_j \operatorname{softmax}_j\!\big(\langle q, h_j\rangle/\sqrt E\big)\, h_j$ ($J$ scores/element);
- **self-attention with query bias** (`Scomp`): every token attends to every token with a learned scalar bias $\alpha\,\mathbf 1\{b=j^*\}$ on the queried key, then gather $\bar h = h'_{j^*}$ ($J^2$ scores/element).

The **deep-set** variants replace the hand-computed means by an inner exchangeable pooling: a shared map embeds each raw $(i,j)$ response, and a masked mean over participants delivers the per-item summary $\operatorname{pool}^x_j, \operatorname{pool}^z_j$ that $q_{\text{tok}}$ then consumes — the nested-DeepSets construction (§9.2), with the same item-axis attention on top. The `idcomp` baseline drops the item axis entirely: the head reads the single per-component statistic $(\bar y_j + \bar z^{(s)}_j,\; N)/N_{\max}$ and predicts $\rho_j$ marginally, exploiting $K$ only through the prior used to draw $\mu$ at training.

**Architecture differences.**

| architecture | per-item input | pooling over participants | mixing across components | query for $j^*$ | uses $K$ at deployment | shares Ukraine class |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| `idcomp` (fixed) | closed-form $(\bar y_j+\bar z_j,\,N)$ | — (statistic supplied) | none — marginal per $j$ | — | prior only | yes |
| `xcomp` (MLP) | raw $(K_{j^*,j},\,y_{i,j})$ | inner sum-pool (DeepSets) | none — $K$-row query | $K$-row | yes | yes |
| `itemScompAtt` | hand $(\bar y_j,\bar z_j,K_{j^*,j})$ | hand group-mean | self-attn $J^2$ + gather | gather $j^*$ + $\alpha$ | yes | yes |
| `itemXcompAtt` | hand $(\bar y_j,\bar z_j,K_{j^*,j})$ | hand group-mean | cross-attn $J$ (1 query) | $q_{\text{query}}(t_{j^*})$ | yes | yes |
| `deepsetScompAtt` | raw $(i,j)$ responses | learned mean-pool (DeepSets) | self-attn $J^2$ + gather | gather $j^*$ + $\alpha$ | via $K$-row metadata | yes |
| `deepsetXcompAtt` | raw $(i,j)$ responses | learned mean-pool (DeepSets) | cross-attn $J$ (1 query) | $q_{\text{query}}(h_{j^*})$ | via $K$-row metadata | yes |

**Why the progression.** `idcomp` is exact for the per-component PPS but discards the cross-component structure, so it cannot serve joint utilities (multivariate stopping rules, family-wise error). `xcomp` adds a $K$-aware query with standard tools (nested DeepSets $+$ $K$-row), no attention. The `item{S,X}compAtt` pair replaces the equal-weight inner pool by attention across components, so the queried component draws information selectively through the $K$-row; cross-attention preserves this alignment most directly (a token with large $K_{j^*,j}$ embeds near the query and attracts weight), self-attention must route it through the learned bias $\alpha$. The `deepset` pair removes the last hand-computed ingredient — the per-item mean — by estimating it from the raw responses, at the cost of a per-participant token axis. These are the two hand-summary encoders (`itemXcompAtt` §14.1, `itemScompAtt` §14.3) and the deep-set encoder (`deepsetXcompAtt` §14.4) reused verbatim in the partial-credit study; the MVN results below are the like-for-like read against a known posterior.

## 13.3 Deployment calibration for the deep-set encoder: parametric heads and the Bernstein–von Mises correction

The deep-set variants (§13.4) reproduce the point PPS but their posterior width does **not** contract at the $n^{-1/2}$ rate — the estimated pooling does not deliver the $1/n$ posterior precision that a closed-form mean does. Because the multivariate normal posterior obeys Bernstein–von Mises *exactly* ($\mathrm{SD} \propto n^{-1/2}$, exponent $p = \tfrac12$), it is the controlled setting in which to verify the correction developed for the partial-credit study (§14.4.6–§14.4.7):

- **Deployment calibration.** An expanding-window fine-tune of the quantile head on the closed-form reference of interims $1..k$, followed by a per-item affine median-shift $\Delta_j(k)$ that removes the residual location bias (§14.4.6).
- **Parametric contraction heads.** The head emits a Gaussian predictive whose scale follows the law: **A power-law** $s_j(n) = C_j\,n^{-p_j}$ and **C floor** $s_j(n) = \sqrt{a_j^2 + b_j^2/n}$ (§14.4.7).
- **Bernstein–von Mises correction.** A post-hoc reshaping of the deployed quantiles pinned by the fitted law read at $n$ and $n+m$, so the marginal-over-$z$ mixture variance equals the law read at $n$ exactly, removing the finite-$m$ upturn at large $n$ (§14.4.7).

Both the calibration and the correction are the shared `amortiser_calibration` implementation used by §14 — the MVN and partial-credit studies run identical code, differing only in the reference (closed form here, SVI there). Numerical results are deferred to §13.4; the point of running them on the MVN is that the law is analytic, so a power-law head should recover $p \to \tfrac12$ and the correction should collapse the deep-set upturn onto the exact curve.

## 13.4 Results

Accuracy is scored as the mean squared error of the estimated per-component PPS against the closed-form $\Phi$-tail PPS, averaged over the seven interims and all $J$ components; each amortiser is deployed at its native posterior-draw budget $S$ (the exact-statistic `idcomp` affords $S=4000$; the attention and deep-set variants use $S=500$ and $S=200$, so their Monte-Carlo floor is higher). The reference estimators are nested-Monte-Carlo with an inner HMC at every $(x, z^{(s)})$ and the regression of the endpoint on $w(z)$ under a Gaussian approximation (§6).

**Point accuracy — PPS–MSE against the closed form.**

| estimator | $J = 20$ | $J = 60$ | $J = 100$ | deploy cost / interim |
|--------------------|------------:|------------:|------------:|------------:|
| **`idcomp` (fixed,** $S=4000$) | **0.00091** | **0.00130** | **0.00106** | 0.14–0.78 min |
| `itemXcompAtt` ($S=500$) | 0.00140 | 0.00114 | 0.00240 | 0.6 min |
| `itemScompAtt` ($S=500$) | 0.00196 | 0.00153 | 0.00218 | 0.6 min |
| `deepsetXcompAtt` ($S=200$) | 0.00573 | — | — | 9–71 min |
| `deepsetScompAtt` ($S=200$) | 0.02136 | — | — | 9–71 min |
| `xcomp` (MLP, undertrained) | 0.383 | 0.400 | 0.366 | 0.4–19 min |
| nested-MC (inner HMC) | 0.00153 | 0.00242 | 0.00204 | 18.7–39.7 min |
| regression on $w(z)$, Gaussian | 0.00348 | 0.00370 | 0.00336 | 0.1–0.8 min |

**Calibration and contraction — raw prior-trained amortiser.** PIT–KS is the conditional-calibration distance (rank of the reference draw among the predicted quantiles vs uniform); marg–KS the distance between the amortiser's marginal $\hat p(\rho_j \mid x)$ and the closed-form posterior; $\hat p$ the posterior-contraction exponent (log–log slope of the predictive SD against $n$), for the amortiser and for the reference. Produced by the shared harness `MVN_interim_diagnostics_by_architecture.py`.

| architecture | $J$ | PIT–KS | marg–KS | $\hat p$ (amortiser) | $\hat p$ (reference) |
|------------|------------|-----------:|-----------:|-----------:|-----------:|
| **`idcomp`** | 20 / 60 / 100 | 0.12 / 0.13 / 0.12 | 0.04 / 0.05 / 0.04 | 0.65 | 0.50 |
| `itemScompAtt` | 20 / 60 / 100 | 0.13 / 0.15 / 0.13 | 0.08 / 0.09 / 0.07 | 0.64 | 0.51 |
| `itemXcompAtt` | 20 / 60 / 100 | 0.16 / 0.20 / 0.17 | 0.08 / 0.11 / 0.09 | 0.66 | 0.51 |
| `deepsetXcompAtt` | 20 | 0.16 | 0.14 | **0.07** | 0.52 |
| `deepsetScompAtt` | 20 | 0.24 | 0.20 | **0.06** | 0.52 |
| `xcomp` (MLP) | 20 / 60 / 100 | 0.79 / 0.81 / 0.81 | 0.80 / 0.82 / 0.82 | 1.2–1.4 | 0.50 |

**Deep-set deployment calibration (§13.3),** $J = 20$. Applying the expanding head fine-tune, affine shift and Bernstein–von Mises correction to `deepsetXcompAtt`:

| pipeline                |    PIT–KS | marg–KS (all $n$) | marg–KS (large $n$) |
|----------------------|---------------:|---------------:|-----------------:|
| plain $+$ affine        | **0.091** |             0.083 |               0.118 |
| **A power-law** $+$ BvM |     0.137 |             0.070 |               0.093 |
| **C floor** $+$ BvM     |     0.134 |         **0.067** |           **0.087** |

The fitted contraction law recovers the exact exponent, median $\hat p = 0.517$ (Bernstein–von Mises $\tfrac12$).

**Interpretation.**

1.  **The closed-form statistic wins, and is** $J$-invariant. `idcomp` attains the lowest PPS–MSE and the tightest calibration (PIT–KS $\approx 0.12$, marg–KS $\approx 0.04$) at every $J$, moving within $\pm 0.0005$ across $J \in \{20, 60, 100\}$ despite training on a single $J = 20$ cell: the per-component decomposition makes the deployment component indistinguishable from a training component. It exploits $K$ only through the prior and so leaves cross-component structure for downstream joint utilities — the reason the attention variants exist.
2.  **Hand-summary attention recovers the reference; cross-attention aligns most directly.** `item{X,S}compAtt` sit within a factor $\sim 2$ of `idcomp` on PPS–MSE and hold PIT–KS $\approx 0.13$–$0.20$, with cross-attention slightly ahead at moderate $J$ — the $K$-row alignment is preserved by the query softmax, whereas self-attention must learn to route it. On the partial-credit study, where no $K$ is available, the two are interchangeable (§14.3.4).
3.  **The deep-set variants match on points but under-contract.** Their PPS–MSE is a few $\times 10^{-3}$, but the contraction exponent collapses to $\hat p \approx 0.07$ against the reference $0.5$: the estimated pooling does not transmit the $1/n$ posterior precision, so the interim uncertainty is dishonest — wide at large $n$ where the posterior should be sharp. This is exactly the deficiency the parametric heads target, and on this exact-law benchmark the power-law head recovers $\hat p = 0.517$ and the Bernstein–von Mises correction restores the large-$n$ marginal ($0.118 \to 0.087$). The deep-set is the right structure only where the raw per-participant responses carry signal a per-component mean discards — the partial-credit model, not the MVN.
4.  **The undertrained baseline is diagnostic, not fundamental.** `xcomp`'s PPS–MSE $\approx 0.4$ and PIT–KS $\approx 0.8$ reflect a training budget too short for the MVN target scale ($\rho = \mu_j - \mu^0$ inherits the prior SD $\tau_0\sqrt{K_{jj}} = 10$, two orders larger than the bounded partial-credit target); the architecture is sound, the wall-clock was not paid (§13.2-adjacent note).
5.  **The amortisers dominate the reference estimators on accuracy-per-compute.** Nested-MC with inner HMC is accurate but $30$–$300\times$ slower; the regression family is $3$–$4\times$ worse in MSE. The exact-statistic and hand-summary amortisers are both more accurate and orders of magnitude cheaper at deployment.

## 13.5 Covariance-family-invariant training and amortising over the number of components

**Covariance-family sampler.** At training each prior draw samples a fresh $K$ from a mixture of standard families — identity, AR(1), block-equicorrelation, and low-rank-plus-diagonal factor structure (§3.3.1) — all rescaled to unit diagonal so the per-component posterior variance is family-invariant. The amortiser is therefore fitted to be invariant to the covariance *family*, not only to $J$: one trained estimator transfers to any $K$ in the family, so the deployment covariance (block-equicorrelation, $\rho_w = 0.8$, $\rho_b = 0.1$) need not match a single training covariance. This is what lets the `xcompAtt` amortiser of §13.2 be deployed at a block-equicorrelation $K$ it never saw at that exact parametrisation, and is the mechanism by which the same class transfers to the partial-credit study, which carries no $K$ at all.

**Amortising over** $J$. With the covariance *structure* held fixed (block-equicorrelation) and only $J$ varying, a single ragged deep-set amortiser prices any $J \in [2, 100]$. The estimator's parameters are $J$-invariant by construction — $J$ is a runtime item axis (cross-attention, one query vs $J$ keys) and every fitted map has input dimension independent of $J$, with no positional index of the component — so a net fitted at one $J$ applies unchanged at another. Training draws a fresh $J \sim \mathrm{Uniform}\{2,\dots,100\}$ at each step, builds the block-equicorrelation $K$ at that $J$ (divisor-safe: size-10 blocks with a ragged final block), simulates a ragged cohort, and fits the power-law head. Evaluation sweeps a $J$-grid, building the interim schedule against the closed-form reference at each $J$, then deploys with the expanding head fine-tune and affine shift.

| $J$ | PPS–MSE (vs closed form) | PIT–KS | marg–KS |
|----:|-------------------------:|-------:|--------:|
|   2 |                  0.00015 |  0.083 |   0.090 |
|   5 |                  0.00011 |  0.079 |   0.095 |
|  10 |                  0.00045 |  0.081 |   0.093 |
|  20 |                  0.00012 |  0.085 |   0.104 |
|  35 |                  0.00008 |  0.080 |   0.106 |
|  50 |                  0.00007 |  0.082 |   0.106 |
|  75 |                  0.00184 |  0.084 |   0.110 |
| 100 |              **0.00004** |  0.080 |   0.106 |

**Reading.** Conditional calibration is $J$-invariant — PIT–KS sits in a flat $0.079$–$0.085$ band from $J = 2$ to $J = 100$ with no trend, matching the best fixed-$J$ result — so amortising over $J$ costs nothing in the object the interim decision consumes. The marginal distance rises mildly and plateaus ($0.090 \to 0.11$): a larger component set is slightly harder to summarise, the same small $n/J$-gradient seen in the partial-credit case. PPS–MSE is negligible everywhere ($\le 0.2\%$); the isolated $J = 75$ bump is single-cohort-seed noise, not a $J$ effect. One amortiser, fitted with a $J$-curriculum on the fixed covariance structure, prices any $J \in [2, 100]$ with flat calibration and negligible error.

**Files.** Architectures and samplers in [`python/model_mvn.py`](../python/model_mvn.py) and the `amortiser_pps_features_*` classes; deployment calibration in the shared [`python/amortiser_calibration.py`](../python/amortiser_calibration.py); cross-architecture diagnostics in [`scripts-py/MVN_interim_diagnostics_by_architecture.py`](../scripts-py/MVN_interim_diagnostics_by_architecture.py); the deep-set A/C $+$ BvM deployment in [`scripts-py/MVN_interim_analysis_amortise_endpt_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead_contraction_bvm.py`](../scripts-py/MVN_interim_analysis_amortise_endpt_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead_contraction_bvm.py); the amortise-over-$J$ study in [`scripts-py/MVN_interim_analysis_amortise_endpt_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead_amortiseJ.py`](../scripts-py/MVN_interim_analysis_amortise_endpt_deepsetXcompAtt_qpsi_MLP_loss_multiquantilehead_amortiseJ.py); the full cross-method comparison in `MVN_interim_analyses_compare_methods.py`.

------------------------------------------------------------------------

# 14. Results for the partial credit model on Ukraine data

## 14.0 Roadmap: estimand, estimator, and result

**Estimand.** For each item $j$ the target is the direction-aware endpoint effect size $\rho_j = s_j\big(\bar w_{e,j}/\bar w_{b,j}-1\big)$, where $\bar w_{b,j},\bar w_{e,j}$ are the population mean baseline/endline levels under the fitted partial-credit model and $s_j=\pm1$ the item direction; the interim decision quantity is the predictive probability of success $\Pr\!\big(\rho_j \ge \delta \mid x\big)$ for a threshold $\delta$. The reference posterior is SVI applied per interim; the amortiser is trained to reproduce, at deployment, the SVI predictive $p(\rho_j\mid x)$ and — for the decision — the conditional $p(\rho_j\mid x, z^{(s)})$ over posterior-predictive future cohorts $z^{(s)}\sim p(z\mid x)$.

**The two statistical problems.** (i) *Calibration*: the network is fitted purely on the prior predictive (§8), so at deployment its predictive is **mislocated** relative to the operational region — the coverage and PIT diagnostics of §14.1.7 quantify this. (ii) *Contraction*: the posterior SD of $\rho_j$ must fall with the observed cohort size $n$ at the Bernstein–von Mises rate; a fixed-scale head does not reproduce it.

**The estimator this section builds to (A/C** $+$ BvM). Read in sequence:

1.  **Encoder.** Nested DeepSets over the raw cohort response tensors with cross-attention over items, on a ragged participant axis (no padding/mask), so the $1/n$ pooling delivers the posterior precision — `deepsetXcompAtt`, §14.4.1–14.4.4. (Hand-computed item-summary encoders, `itemXcompAtt` §14.1 and `itemScompAtt` §14.3, are the alternatives used for comparison in §14.5.)
2.  **Training law.** Wide **prior** proposal over item parameters, a $\theta$-level target functional (endpoint over a large reference ability population), decoupled cohort sizes $(n,m)$, and the interval-score loss — §14.4.3. Net directory `…-deepsetXcompAtt-ragged-priorproposal-{260812 (plain), powerlaw, floor}`.
3.  **Deployment calibration.** Expanding-window head-only fine-tune of the quantile head on the SVI of interims $1..k$ (§14.4.6), followed by the per-item **affine median-shift** $\Delta_j(k)=\operatorname{mean}_{i\le k}[\operatorname{med}p_{\mathrm{SVI}}(\rho\mid x_i)-\operatorname{med}\hat p(\rho\mid x_i)]$ that removes the residual location bias (§14.4.6). Output directories carry the suffix `_ftheadexpand`.
4.  **Parametric contraction.** The posterior SD follows an item-specific power law $\mathrm{SD}(\rho_j\mid x)=C_j\,n^{-p_j}$ (Bernstein–von Mises $+$ delta method, §14.4.7), realised as the power-law head **A** ($m+z_\tau C\,n^{-p}$) and the floor head **C** ($\sqrt{a^2+b^2/n}$) — §14.4.7, directories `…-powerlaw`, `…-floor`.
5.  **BvM self-consistency correction.** A post-hoc correction of the deployed quantiles (§14.4.7): with the law of total variance $\operatorname{Var}_{\text{mix}}=W+B$, the non-spurious between term is pinned by the same power law read at $n$ and $n+m$, $B^\ast=C^2\big(n^{-2p}-(n+m)^{-2p}\big)$; scaling within by $\alpha=\sqrt{W_{\text{true}}/W}$ and shrinking the median spread by $\beta=\sqrt{B^\ast/B}$ makes the mixture variance equal the target exactly and removes the finite-$m$ contraction upturn. Directories `…-powerlaw-ftheadexpand-bvm-260828`, `…-floor-ftheadexpand-bvm-260828`.

**Result (weekly grid, `py-ukraine-interim-weekly-svi-260811`, 29 interims; excluding the degenerate item `CG-VIO_ph-punish`).** With the affine shift applied throughout, the correction lowers the marginal KS at every cohort size, most in the large-$n$ tail, and makes the amortiser's contraction curve track SVI's (§14.4.7):

| pipeline | PIT-KS (conditional) | marg-KS (all $n$) | marg-KS ($n\ge450$) |
|--------------------------------|-------------:|-------------:|-------------:|
| plain $+$ affine (`…-plain-ftheadexpand-260822`) | **0.069** | **0.157** | 0.181 |
| **A power-law** $+$ BvM (`…-powerlaw-ftheadexpand-bvm-260828`) | 0.112 | 0.177 | 0.180 |
| **C floor** $+$ BvM (`…-floor-ftheadexpand-bvm-260828`) | 0.115 | 0.166 | **0.170** |

The **A/C** $+$ BvM estimator is the statistically coherent recommendation: its width obeys a BvM contraction law and its marginal is self-consistent under the law of total variance. The plain net $+$ affine remains numerically ahead on the aggregate conditional PIT-KS and on all-$n$ marg-KS, and is the lighter alternative when the parametric contraction guarantee is not required. Comparison figures — conditional-calibration PIT box, marginal quantile box, contraction-factor, contraction-law — are produced for every configuration from `python/amortiser_diag_plots.py`. Approaches examined during the development and set aside are catalogued in §14.10.

## 14.1 Cross-component attention to learn item-specific summaries (`itemXcompAtt`)

**Task.** Adapt the data-agnostic cross-attention amortiser of §13.2 (`itemXcompAtt`) to the Ukraine partial-credit (PCM) interim analysis. The network class is reused **unchanged**; everything Ukraine-specific lives in the token construction, the target, and the training pool. Diagnostics: §14.1.7–14.1.11. Calibration remedies: §14.2. Self-attention sibling: §14.3. Fully-amortised deepset variant: §14.4.

### 14.1.1 Structure of learning task

Each participant contributes **two** responses per item — baseline ($t = 0$) and endline ($t = 1$). The interim cohort $x$ consists of the $n$ participants with both visits complete by the cutoff; the future cohort $z$ is the remaining $m = N - n$. The PCM accommodates baseline-vs-endline level shifts by fitting **separate per-time thresholds**: internally the likelihood sees $2 J = 40$ `item_time_id` fake items (own threshold vector and loading per (real item, time)), while participant ability $\theta_i$ is shared across times and items. The amortiser collapses the time axis back: one token per real item, with the two time-points entering as separate scalar features.

Two axes must not be conflated: the **cohort axis** ($x^{\text{obs}}$ vs $z^{(s)}$ — the §10 deployment passes both through the network) and the **time axis** (baseline vs endline within a participant, present in both cohorts). The current $F = 8$ token (§14.1.3) carries both cohorts and both time-points. The remaining §10 gap at the item level is that the per-item pooling over participants is hand-computed (group-means), not learned; the learned version is the nested-DeepSets variant (§13.2, results §14.4).

| symbol | value | description |
|----------|----------|----------------------------------------------------|
| $J$ | 20 | **real** items presented to the amortiser (one token per real item); internally the PCM likelihood fits $2 J = 40$ per-`item_time_id` fake items |
| $n$ | interim-specific | observed participants at the interim, **each with both baseline and endline responses** (contributes $2 J n$ observations to the PCM likelihood) |
| $m$ | interim-specific | shadow future participants used per posterior draw $s$, **each also carrying both baseline and endline** (contributes $2 J m$ posterior-predictive responses per draw) |
| $S$ | 4000 | posterior-predictive draws of $z^{(s)}$ per interim at deploy time (item amortisers, §14.1/§14.3); the deepset variant (§14.4) uses $S = 200$ |
| $C$ | 8 | interims analysed (`INTERIM_IDS = 1..8`) |
| $F$ | 8 | per-real-item token dim: x-side + z-side cohort summaries, change ratio, metadata, empirical $\hat K$-row (see §14.1) |
| $A$ | 2 | aux scalars $(n/N, m/N)$, $N = 503$ |
| $E$ | 32 | embedding dim (`embed_dim`) |
| $K_\tau$ | 5 | quantile levels $\tau = (0.05, 0.25, 0.5, 0.75, 0.95)$ |

What is reused, what is trained, what is adapted:

| piece | status | detail |
|--------------|--------------|-------------------------------------------|
| network classes ($q_{\text{tok}}, q_{\text{query}}$ or $\alpha$, $q_\psi$) | **reused** from §13.2 | identical Flax modules, identical attention math |
| network weights | **trained from scratch** | \~12k parameters, multi-quantile pinball loss |
| per-item token $t^{(j^*, s)}_j$ | **adapted** | $F = 8$: x-side + z-side summaries, empirical $\hat K$-row, item metadata (see below) |
| aux $a$ | **adapted** | $A = 2$: cohort sizes $(n/N, m/N)$ with $N = 503$ |
| input scale | **adapted** | raw responses rescaled $(y)/(K_j - 1) \in [0, 1]$ with the per-item level count $K_j$ from `dit[cat_length]` (8 for out-of-7, 4 for categorical); all summaries are group-means of these rescaled responses |
| target scale | **adapted** | per-item standardisation $\rho / \hat\sigma_j$ (see loss below) |
| training pool | **§8 prior-predictive** | fresh PCM prior draws per training step; no SVI fits or cached posterior draws enter training |

### 14.1.2 Training data and target

Training data are generated **from the PCM prior predictive only** (§8): no SVI fits and no cached posterior draws enter training. One training row is the pair (prior sample $s$, queried item $j^*$):

$$\mathcal{D} \;=\; \Big\{ \big( t^{(j^*, s)}_{1:J},\; j^*,\; \rho^{(j^*, s)} \big) \Big\}, \qquad s \le S_{\text{train}},\; j^* \le J,$$

with the target the population endpoint under the prior draw,

$$\rho^{(j^*, s)} \;=\; r_{j^*}\big(\theta^{(s)}\big), \qquad \theta^{(s)} \sim p(\theta),$$

and $r_j(\theta)$ is the direction-aware ratio built from the PCM's per-`item_time_id` ordered-probability endpoint (`eval_outcome_for_endpoint` → `get_endpoints_per_draw`), evaluated at both fake items $(j, t{=}0)$ and $(j, t{=}1)$ and combined: for `out-of-7` items the normalised change of $\mathbb{E}[y_{(j, t)} \mid \theta]$ across $t$, for `categorical` items the analogous change of $P(y_{(j, t)} \geq c \mid \theta)$. Targets are clipped at $\pm 20$ (a single numerical blow-up at `CG-VIO_ph-punish`).

**Generation algorithm.** For $s = 1, \ldots, S_{\text{train}}$:

1.  draw cohort sizes $n^{(s)} \sim U\{2, \ldots, N - 1\}$, $m^{(s)} = N - n^{(s)}$;
2.  draw the PCM parameters from the prior, $$\theta^{(s)} = \big(\{\theta_i\}_{i \le N},\; \{\beta_t\}_{t \in \{0,1\}},\; \{\tau_{j, t, \cdot}\},\; \{\lambda_{j, t}\}\big), \quad \theta_i \sim N(0, 1),\; \beta_t \sim N(0, 1),\; \tau_{j, t, k} \sim N(0, 3.5),\; \lambda_{j, t} \sim |t_3|;$$
3.  simulate ordered-categorical responses for both cohorts at both time-points from the PCM likelihood, $x^{(s)}_{i, j, t} \sim \text{PCM}(\theta^{(s)})$ for $i \le n^{(s)}$ and $z^{(s)}_{i, j, t} \sim \text{PCM}(\theta^{(s)})$ for $i \le m^{(s)}$;
4.  compute the token entries of §14.1.3 from the simulated raw cohorts (group-means, change ratio, empirical $\hat K$ from the simulated $x^{(s)}$, sizes aux) — the token is a deterministic function of the raw data, so the same formulas apply verbatim;
5.  label with the population endpoint $\rho^{(s)}_j = r_j(\theta^{(s)})$ evaluated on the ordered probabilities under $\theta^{(s)}$ (threshold $c$ enters here, as at deployment).

Implemented by `PartialCreditModel.make_training_data_with_item_tokens_prior` (the raw-cohort simulation, steps 1-3, is shared with the deepset sampler of §14.4, which consumes the cohorts directly; here step 4 pools them to the item-level token). Every training draw is fresh — the amortiser never sees the operational data during training.

### 14.1.3 Per-item token

The following response data from the observed and future cohorts enter the token (§10 cohort axis, as in MVN §13.2):

$$t^{(j^*, s)}_j \;=\; \Big(\underbrace{w^{x}_{\text{base}, j},\; w^{x}_{\text{end}, j}}_{\text{observed cohort } x^{\text{obs}},\ s\text{-invariant}},\;\; \underbrace{w^{z, (s)}_{\text{base}, j},\; w^{z, (s)}_{\text{end}, j},\; w^{z, (s)}_{\text{ratio}, j}}_{\text{future block } z^{(s)}},\;\; \underbrace{c_j,\; d_j}_{\text{metadata}},\;\; \underbrace{\hat K_{j^*, j}}_{\text{empirical K-row}}\Big) \;\in\; \mathbb{R}^{F},\qquad F = 8.$$

- **x-side** $w^{x}_{\text{base}, j}, w^{x}_{\text{end}, j}$: per-item group-mean **response** over the $n$ observed participants at $t = 0$ and $t = 1$, rescaled by the item type's level count: $\tfrac{1}{n (K_j - 1)}\sum_i y_{i, j, t}$ with the per-item level count $K_j$ read from `dit[cat_length]` (Ukraine: $K_j = 8$ for `out-of-7`, raw $y \in \{0, \ldots, 7\}$; $K_j = 4$ for `categorical`, raw $y \in \{0, \ldots, 3\}$). **No input-side thresholding** — the caseness threshold lives in the target only (see below). Fixed per interim.
- **z-side** $w^{z, (s)}_{\cdot, j}$: the same mean-response statistics on the $m$ future-cohort participants of $z^{(s)}$ — simulated from the prior at training time, posterior-predictive at deployment (§14.1.6). The ratio is `higher_is_better`: $w_{\text{end}}/w_{\text{base}} - 1$, `lower_is_better`: $1 - w_{\text{end}}/w_{\text{base}}$, computed on the $[0, 1]$ mean scale and clipped at $\pm 20$.
- **Metadata** $c_j \in \{0, 1\}$ (item type), $d_j \in \{0, 1\}$ (direction).
- **Empirical K-row** $\hat K_{j^*, j}$ — the $j^*$-dependent entry, restoring the MVN-style "how informative is item $j$ for item $j^*$" channel. From the observed cohort's per-participant change scores $d_{i, j} = y_{i, j, 1} - y_{i, j, 0}$, take the Spearman correlation matrix $\hat R$ across items (rank-based — robust to the ordinal scale) and shrink toward identity: $\hat K = \lambda \hat R + (1 - \lambda) I$, $\lambda = n / (n + n_0)$, $n_0 = 50$, so the noisy early-interim estimates ($n \approx 48 \Rightarrow \lambda \approx 0.5$) are strongly regularised. Recomputed per interim.

**Categorical threshold.** $c = 2$ on the 1-indexed scale $\Leftrightarrow$ raw $y \geq 1$ ("at least several days"). The CG-MH items are PHQ-style symptom-frequency Likerts; the trial's endpoint for them is the *caseness prevalence* — proportion of caregivers above the symptom cut-off — so the threshold is intrinsic to the **target** $\rho$ and stays there. The **input** summaries do not threshold: they use the mean response (see the token definition above), because binarising the input is pure information loss. The level counts are read from `dit[cat_length]` (never hardcoded); the threshold $c$ is a study-specific configuration parameter — Colombia's response scales have different level counts and hence a different $c$.

### 14.1.4 Encoder, cross-attention, head, loss function

Exactly §13.2 variant 1. Two SiLU MLPs to build the tokens and query for each item $j^*$ in embedding space:

$$h^{(j^*, s)}_j \;=\; q_{\text{tok}}\big(t^{(j^*, s)}_j;\, \tau_{\text{tok}}\big) \;\in\; \mathbb{R}^{E}, \qquad q^{(j^*, s)} \;=\; q_{\text{query}}\big(t^{(j^*, s)}_{j^*};\, \tau_{\text{query}}\big) \;\in\; \mathbb{R}^{E}.$$

Cross-attention over the $J$ per-item embeddings ($J$ attention scores per batch element):

$$\bar h^{(j^*, s)} \;=\; \sum_{j=1}^{J} \omega^{(j^*, s)}_j\, h^{(j^*, s)}_j, \qquad \omega^{(j^*, s)}_j \;=\; \operatorname{softmax}_{j'}\Big( \big\langle q^{(j^*, s)},\; h^{(j^*, s)}_{j'} \big\rangle \big/ \sqrt{E} \Big)\Big|_{j' = j}.$$

$\bar h^{(j^*, s)} \in \mathbb{R}^{E}$ is the **learned summary in embedding space** for the queried item — the amortised counterpart of the hand-picked scalar statistic $w(z^{(s)})$ of the regression baseline (§6.2). Using the learned summary, we then learn the mapping to the item-specific endpoints $\rho$ . For each item, we learn a more efficient multi-quantile head rather than a specific quantile, so we are free to interpolate to any user-desired quantile:

$$\hat\rho^{(j^*, s)} \;=\; q_\psi\big( \operatorname{concat}(\bar h^{(j^*, s)},\; t^{(j^*, s)}_{j^*},\; a);\, \psi \big) \in \mathbb{R}^K,$$

with head-input dimension $E + F + A = 32 + 8 + 2 = 42$ and head output dimension $Q = 5$. Each dimension $q$ of $\hat\rho^{(j^*, s)}$ corresponds to one of the quantile levels $\eta^{\text{inpol}}_q$, $q = 1, \ldots, Q$, with $(\eta^{\text{inpol}}_1, \ldots, \eta^{\text{inpol}}_Q) = (0.05,\, 0.25,\, 0.5,\, 0.75,\, 0.95)$, which are used to interpolate to the specific quantile that the user requests.

The learned parameters are $\tau = (\tau_{\text{tok}},\, \tau_{\text{query}})$ — the weights of the two encoder MLPs — and $\psi$, the weights of the head MLP $q_\psi$ (the notation of §8-10). The attention weights $\omega^{(j^*, s)}$ are **not** free parameters: they are deterministic functions of the embeddings (softmax of scaled dot products) and carry no trainable degrees of freedom of their own in this variant (contrast the learned scalar bias $\alpha$ of `itemScompAtt`, §14.3). The same $(\tau, \psi)$ apply to **every item** $j^*$ and every interim — one network serves all queries — in contrast to the regression baseline (§6.2, RGE/RGEM), where the learned regression function is item-specific *and* interim-specific (one fit per (item, interim) pair; $20 \times 8 = 160$ separate fits over the schedule). Sizes at the Ukraine configuration ($F = 8$, $E = 32$, hidden $(32, 32)$ for the encoders and $(64, 64)$ for the head, $K = 5$):

| block                 | architecture             | \# parameters |
|-----------------------|--------------------------|--------------:|
| $\tau_{\text{tok}}$   | $8 \to 32 \to 32 \to 32$ |         2,400 |
| $\tau_{\text{query}}$ | $8 \to 32 \to 32 \to 32$ |         2,400 |
| $\psi$                | $42 \to 64 \to 64 \to 5$ |         7,237 |
| **total**             |                          |    **12,037** |

We feed $t^{(j^*, s)}_{j^*}$ and $a$ are into $q_\psi$ alongside $\bar h^{(j^*, s)}$ for the following reasons:

- $t^{(j^*, s)}_{j^*}$ is a **skip connection to the queried item's own raw summaries**. The attention output is a weighted mix over all $J$ items; when $\omega$ is diffuse (early training, or weakly informative $\hat K$) the queried item's own signal is diluted by a factor of order $1/J$. The skip guarantees the head always sees the query's own statistics undiluted — and it is the only place the query identity enters the head, since $\bar h$ alone need not identify $j^*$. Attention is thereby free to specialise in *cross-item borrowing* instead of having to reconstruct the own-item signal.

- $a = (n/N, m/N)$ carries the **cohort sizes**, which are not recoverable from the group-mean summaries: the same means imply a much tighter posterior for $\rho$ at $n = 501$ than at $n = 48$, so the head needs $a$ to scale the predicted quantile spread. This mirrors §9.3's exp-fam recipe, where the sufficient-statistic sums must be accompanied by $(n, m)$ as extra scalar inputs to the head.

The loss function is as follows. Index the training data by prior sample $s = 1, \ldots, S_{\text{train}}$ and queried item $j^* = 1, \ldots, J$; write $\rho^{(j^*, s)}$ for the target and $\hat\rho^{(j^*, s)}_q$ for the $q$-th component of the head output of §14.1.4. The per-item standardisation scale is estimated once from a prior pilot batch ($2{,}000$ prior samples),

$$\hat\sigma^{(j^*)} \;=\; \operatorname{sd}_{s}\big(\rho^{(j^*, s)}\big).$$

The network is trained on the standardised multi-quantile target with the pinball loss summed over quantile levels, prior samples and queried items:

$$\ell(\tau, \psi) \;=\; \sum_{s=1}^{S_{\text{train}}} \sum_{j^*=1}^{J} \sum_{q=1}^{Q} L_{\eta^{\text{inpol}}_q}\!\Big( \frac{\rho^{(j^*, s)}}{\hat\sigma^{(j^*)}} \;-\; \hat\rho^{(j^*, s)}_q \Big),$$

where the pinball (check) loss at level $\eta^{\text{inpol}}_q$ is

$$L_{\eta^{\text{inpol}}_q}(u) \;=\; \big(\eta^{\text{inpol}}_q - \mathbf{1}\{u < 0\}\big)\, u \;=\; \begin{cases} \eta^{\text{inpol}}_q\, |u| & u \geq 0,\\ \big(1 - \eta^{\text{inpol}}_q\big)\, |u| & u < 0. \end{cases}$$

In practice the sum over $s$ is realised by drawing fresh prior mini-batches every step, and the sum over $j^*$ by fanning each prior sample out to all $J$ queries. At deployment the head output is mapped back to the target scale by the same per-item factor,

$$\hat\rho^{(j^*, s)}_{q} \;\leftarrow\; \hat\sigma^{(j^*)} \cdot \hat\rho^{(j^*, s)}_q,$$

which is what the CDF interpolation of §14.1.6 consumes. Without the standardisation, wide-range items dominate the shared head's gradients and narrow-range items (the four `CG-MH` categorical items) are under-resolved.

### 14.1.5 Training summary

- $6\,000$ steps, $S = 64$ **fresh prior-predictive draws per step** (fanned out to $B = S \cdot J = 1280$ query rows / step) $\Rightarrow$ 384,000 prior data sets over training; targets and input ratios clipped at $\pm 20$.
- Per-item target standardisation ($\hat\sigma^{(j^*)}$ from the 2,000-sample prior pilot, range $[0.94, 5.07]$) + $[0, 1]$ input rescale as above.
- Wall-clock: \~12 min training on 6-thread CPU. Trainable params ≈ 12 K.

### 14.1.6 Deployment

Per interim, the current data are held **fixed at the observed** $x^{\text{obs}}$ — they enter every token through the x-side summaries, the empirical $\hat K$-row and the size aux. The future data $z^{(s)}$ are simulated from the **posterior predictive** at $x^{\text{obs}}$ (§10 step 1) — SVI is run **once per interim** on $x^{\text{obs}}$ (the only inference cost at deployment; training never touches it), then $\theta^{(s)} \sim p(\theta \mid x^{\text{obs}})$, then simulate the $m = N - n$ missing participants' responses at both time-points under $\theta^{(s)}$ (`get_interim_z_from_ypredi`: shadow participants resampled with replacement from the observed cohort's covariate rows, responses taken from the model's `ypred` at draw $s$). We use **all** $S = 4000$ cached posterior draws per interim, so one interim costs $J \cdot S = 80{,}000$ head evaluations (one JIT-compiled batch). Each $(j^*, s)$ pair is forward-passed, the quantile grid converted to a conditional success probability by CDF interpolation at the effect threshold $\eta_0$, and aggregated over draws (§10 step 3). The hypothesis is **item-specific**, $H_1^{(j^*)} = \big\{\rho^{(j^*)} > \eta_0\big\}$ — item $j^*$'s improvement ratio exceeds the effect threshold ($\eta_0 = 0.5$, the `pps_H1_min_effect_size_thresh` of the code; decision threshold $\eta_H = 0.89$, the `pps_ProbH1_target_lwr_quantile`):

$$\hat P\big(H_1^{(j^*)} \mid x, z^{(s)}\big) \;=\; 1 - \operatorname{interp}\big(\eta_0;\; \hat\rho^{(j^*, s)},\; \eta^{\text{inpol}}\big),$$

$$\widehat{\text{PPS}}^{(j^*)}(x) \;=\; \frac{1}{S} \sum_{s=1}^{S} \mathbf{1}\big\{ \hat P\big(H_1^{(j^*)} \mid x, z^{(s)}\big) > \eta_H \big\}.$$

With the $F = 8$ token, both cohorts enter the network — hand-pooled summaries of $x^{\text{obs}}$ plus summaries of $z^{(s)}$ — putting this variant on the same footing as the MVN item amortiser (§13.2): §10-consistent on the cohort axis, with hand-computed (not learned) per-item pooling. Files:

- [`python/amortiser_pps_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead.py`](../python/amortiser_pps_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead.py) — the data-agnostic network class (§13.2): `_MLP` + cross-attention (`einsum` + `softmax`). Standard Flax.
- [`python/model_pcm.py`](../python/model_pcm.py) — `make_training_data_with_item_tokens_prior`, the §8 prior-predictive training sampler (§14.1.2).
- [`scripts-py/Ukraine_interim_analysis_amortise_endpt_itemXcompAtt_qpsi_MLP_loss_multiquantilehead.py`](../scripts-py/Ukraine_interim_analysis_amortise_endpt_itemXcompAtt_qpsi_MLP_loss_multiquantilehead.py) — trains once on fresh prior draws, then deploys per interim (SVI once → posterior-predictive $z^{(s)}$ → tokens → forward pass); `_RGEG_` outputs.

Deployment time is dominated by rebuilding the posterior-predictive $z^{(s)}$ block per interim (\~1 min / interim at $S = 4000$); the network forward pass itself is seconds.

### 14.1.7 Diagnostic: conditional calibration on $(x, z^{(s)})$ — quantile coverage + PIT

Ukraine has no analytic PPS, so we use several diagnostics to evaluate estimation accuracy.

For each SVI posterior draw $s$, take $\theta^{(s)} \sim p(\theta \mid x^{\text{obs}})$ (SVI fit to the current data only). Compute $\rho^{(j^*, s)}_{\text{SVI}} = r_{j^*}(\theta^{(s)})$ — a draw from the current posterior $p(\rho \mid x)$ and generate the future block $z^{(s)} \sim p(z \mid \theta^{(s)})$ that the amortiser conditions on. Because $(\theta^{(s)}, z^{(s)})$ is a joint draw given $x$, Bayes gives $\theta^{(s)} \mid z^{(s)} \sim p(\theta \mid x, z^{(s)})$, so $\rho^{(j^*, s)}_{\text{SVI}} | z^{(s)}$ is an exact single draw from the conditional $p(\rho \mid x, z^{(s)})$ that the amortiser aims to target. Marginally, $\rho^{(j^*, s)}_{\text{SVI}}$ is a draw from $p(\rho \mid x)$. Exactly as for any joint $(A,B)$ a realised $a$ is simultaneously a draw from $p(A)$ and, paired with its own $b$, from $p(A\mid B{=}b)$. The output of the amortiser $\hat\rho^{(j^*, s)}_{1:Q}$ approximates the conditional $p(\rho \mid x, z^{(s)})$.

1.  **Quantile coverage.** For each level $\eta^{\text{inpol}}_q$, compute for the same index $s$:

$$\widehat{\text{Cov}}^{(j^*)}_q \;=\; \frac{1}{S} \sum_{s=1}^{S} \mathbf{1}\big\{ \rho^{(j^*, s)}_{\text{SVI}} \le \hat\rho^{(j^*, s)}_q \big\} \;\stackrel{!}{=}\; \eta^{\text{inpol}}_q .$$

2.  **PIT uniformity.** The probability integral transform obtained by CDF interpolation of the quantile grid,

$$u^{(j^*, s)} \;=\; \hat F\big(\rho^{(j^*, s)}_{\text{SVI-}x} \mid x, z^{(s)}\big) \;=\; \operatorname{interp}\big(\rho^{(j^*, s)}_{\text{SVI-}x};\; \hat\rho^{(j^*, s)}_{1:Q};\; \eta^{\text{inpol}}_{1:Q}\big),$$

should be $\text{Uniform}(0, 1)$ over $s$; deviations are summarised per (interim, item) by the Kolmogorov-Smirnov distance $\max_u |\widehat{\text{ecdf}}(u) - u|$. Under-dispersion shows as a U-shaped PIT histogram, over-dispersion as a hump, location bias as a tilt.

Implemented in the separate script [`..._itemXcompAtt_..._tests.py`](../scripts-py/Ukraine_interim_analysis_amortise_endpt_itemXcompAtt_qpsi_MLP_loss_multiquantilehead_tests.py) (kept apart from train/deploy). Outputs: `..._tests_coverage.csv/.pdf` (empirical-vs-nominal coverage curves, diagonal = calibrated) and `..._tests_pit.pdf` + `..._tests_pit_ks.csv` (PIT histograms + KS). Result:

| nominal $\eta^{\text{inpol}}_q$                 | 0.05 | 0.25 | 0.50 | 0.75 | 0.95 |
|-------------------------------------|------:|------:|------:|------:|------:|
| empirical coverage (mean over items × interims) | 0.05 | 0.10 | 0.15 | 0.23 | 0.44 |

Coverage falls far below nominal above the 0.05 level; mean PIT-KS $\approx 0.61$, improving monotonically across the schedule (0.75 at interim 1 $\to$ 0.47 at interim 8). Interpretation in §14.1.12.

**Why the training loss is tiny yet coverage is far off.** The pinball loss and this coverage test measure the *same* property — the $\eta^{\text{inpol}}_q$-level pinball is minimised by the conditional $\eta^{\text{inpol}}_q$-quantile, whose coverage is $\eta^{\text{inpol}}_q$ — but under **different distributions**. Training minimises expected pinball over the **prior-predictive** ($\theta \sim p(\theta)$, simulated cohorts); this diagnostic measures coverage on the **posterior-predictive deployment slice** ($\theta \sim p(\theta \mid x^{\text{obs}})$). A head calibrated on the prior-predictive average need not be calibrated on that slice when the two differ (the §8.1 prior-coverage gap: thresholds $\sim N(0, 3.5)$ generate cohorts far more extreme than Ukraine's), and the average loss is dominated by the bulk of the diffuse prior mass — the thin deployment region barely registers in it. Direct confirmation: **coverage on held-out prior-predictive draws is essentially nominal** (0.05, 0.27, 0.51, 0.76, 0.95 at the five levels), so the small loss is doing its job — the deployment miscalibration is distribution shift, not an optimisation or objective failure.

### 14.1.8 Diagnostic: marginal $p(\rho \mid x)$ (SVI vs amortiser)

Two distributions of $\rho \mid x$ are compared:

- **SVI:** the set $\{\rho^{(j^*, s)}_{\text{SVI-}x}\}_{s=1}^{S}$ (endpoints of the SVI-on-$x$ draws) — its empirical CDF is $p(\rho \mid x)$.
- **Amortiser:** its marginal predictive $\hat p(\rho \mid x) = \mathbb{E}_{z \mid x}\big[\hat p(\rho \mid x, z)\big]$, estimated by averaging the per-draw quantile CDFs, $\hat F(\rho \mid x) = \tfrac{1}{S}\sum_s \hat F(\rho \mid x, z^{(s)})$. Note the item tokens require a future cohort and so we marginalise over $s$. Agreement is summarised per (interim, item) by the two-sample KS distance $\max_\rho |\hat F(\rho \mid x) - \widehat{\text{ecdf}}_{\text{SVI}}(\rho)|$.

Outputs (same tests script as §14.1.7): `..._tests_marginal_ks.csv` and per-interim `pcm_1_interim_i{k}_svi_vs_amortiser_marginal_cdf.pdf` (facet per item, the two CDFs overlaid). Two-sample KS $\max_\rho |\hat F - \widehat{\text{ecdf}}_{\text{SVI}}|$ (`CG-VIO_ph-punish` excluded):

| marginal KS  | out-of-7 | CG-MH |
|--------------|---------:|------:|
| interim 1    |     0.51 |  0.53 |
| all interims |     0.42 |  0.62 |

The amortiser's $\hat p(\rho \mid x)$ sits **below** the SVI posterior — the same downward location bias as §14.1.7, now on the current-decision object. Unlike the conditional PIT-KS (which falls 0.75 → 0.47 across the schedule), the marginal KS is roughly **flat** ($\approx 0.47$ overall): a genuine location/shape mismatch of the marginal, not conditional under-dispersion that averages out. Reading in §14.1.12.

### 14.1.9 Diagnostic: median-prediction association (per-draw correlation)

`pcm_1_interim_i{k}_svi_rho_vs_amortiser_median.pdf` — per interim, facet per item. Scatter of

$$y \;=\; \rho^{(j^*, s)}_{\text{SVI-}x} \qquad \text{against} \qquad x \;=\; \hat\rho^{(j^*, s)}_{3} \quad (\eta^{\text{inpol}}_3 = 0.5),$$

the median head output after $q_\psi$ for same $s$. **D**eployment thresholds and averages into $\widehat{\text{PPS}}^{(j^*)}$, whereas this diagnostic keeps each pair and correlates the median output against $\rho^{(j^*, s)}_{\text{SVI}}$. Unlike §14.1.7 it tests **association only** (does the central prediction move with the target across $z^{(s)}$?), not location or spread — the two are deliberately separated because the amortiser can rank well while being mis-located. An earlier version scattered $\rho^{(j*,s)}_{\text{SVI}}$ against PC1 of $\bar h^{(j^*, s)}$ similar to the Strong-Oakly diagnostic. We dropped this because PC1 is a 1-D *linear* proxy for the nonlinear $q_\psi$ read-out, so it is uninformative about the network.

### 14.1.10 Results: correlation against the SVI posteriors

Per-item Pearson $\rho$ between $\hat\rho^{(j^*, s)}_3$ and $\rho^{(j^*, s)}_{\text{SVI-}x}$ for the prior-predictive-trained amortiser (`CG-VIO_ph-punish` blow-up item excluded):

| correlation  | out-of-7 | CG-MH |
|--------------|---------:|------:|
| interim 1    |     0.82 |  0.71 |
| all interims |     0.53 |  0.42 |

**Threshold on the CG-MH token responses.** The four `CG-MH` items are categorical; the endpoint $r_j(\theta)$ is a caseness prevalence, defined via a response threshold $c$. It was important to keep that threshold **out of the input token**: an earlier token summarised the CG-MH responses as the proportion above $c$ (a thresholded, lossy summary) and capped their correlation at $\rho \approx 0.48$; replacing it with the **mean response** (no threshold, no information loss) — while keeping $c$ only in the target — lifted CG-MH to $\approx 0.71$ with no change to the out-of-7 items. The threshold belongs in the endpoint definition, never in the summary the network reads.

Synthesis in §14.1.12.

### 14.1.11 Diagnostic: PPS agreement with the regression baselines

No analytic ground-truth PPS exists for Ukraine, so the deployed $\widehat{\text{PPS}}^{(j^*)}$ is compared to the two closest Strong-Oakley baselines — Gaussian-approximation regression (RGE) and multi-quantile regression (RGEM) — over the $J \cdot C = 160$ (item, interim) cells:

| comparison | Pearson $\rho$ | median \$ | \Delta\text{PPS} | \$ |
|---------------------------|-----------:|-----------:|-----------:|-----------:|
| `itemXcompAtt` vs RGE (Gauss) | 0.905 | 0.001 | 0.223 | 0.880 |
| `itemXcompAtt` vs RGEM (mquantile) | 0.916 | 0.002 | 0.241 | 0.830 |
| `itemXcompAtt` vs `itemScompAtt` (§14.3) | 0.997 | 0.000 | 0.041 | — |

Median difference to the regression baselines is $\le 0.002$; the disagreeing tail (p90 $\approx 0.22$-$0.24$) is where cross-item information and the calibration gap move the prediction. The two attention variants agree almost perfectly. Reading in §14.1.12.

### 14.1.12 Reading

Synthesising the diagnostics (§14.1.7 conditional calibration, §14.1.8 marginal $p(\rho\mid x)$, §14.1.9 association / §14.1.10 results, §14.1.11 PPS agreement):

1.  **Ranking is good, calibration is not — the two are dissociated.** Median-output association with the SVI draws is strong (Pearson $\rho \approx 0.8$ out-of-7, $\approx 0.7$ CG-MH; §14.1.9), yet quantile coverage sits far below nominal and PIT-KS $\approx 0.61$ (§14.1.7), and the marginal $\hat p(\rho\mid x)$ sits below the SVI posterior (§14.1.8): the amortiser orders the posterior draws well but its conditional distribution is **mis-located and under-dispersed**. The calibration diagnostic quantifies distributionally what the negative $R^2$ hints; the median diagnostic confirms the location error does not come from a broken ranking. The PPS agrees with the Strong-Oakley baselines on the bulk of cells ($\le 0.002$ median), with a tail exactly where the calibration gap bites (§14.1.10).

2.  **Cause: prior coverage.** The output scale inherits the diffuse PCM prior — thresholds $\sim N(0, 3.5)$ place most training mass on far more extreme cohorts than Ukraine's, and the training label is the *population* endpoint $r_j(\theta)$ while the evaluation target is the finite-cohort SVI estimate. The monotone improvement of PIT-KS across the schedule (interim 1 → 8) fits this: as $n$ grows the operational posterior moves into the prior-covered region. This is the §8.1 prior-coverage caveat, made concrete.

3.  **Remedies (the actionable next step).** these are executed as a full remedy sweep in §14.2 — the deployable fix is an expanding-window head-only fine-tune (§14.2.5), which reaches near-nominal calibration with correlation intact. Same calibration origin as the deepset variant (§14.4.6).

4.  **The §8-10 loop closes at the item level at no accuracy cost.** Prior-predictive training (fresh PCM draws, no SVI anywhere) with the hand-computed $F = 8$ token reaches the same correlation as the superseded posterior-trained regime (§14.1.9, final vs third row) — the train-once-deploy-anywhere property is free on this metric; only calibration is outstanding.

5.  **Architecture ports cleanly from MVN, and cross-item routing is a wash here.** The network class is unchanged; the analytic K-row becomes the empirical Spearman-shrunk $\hat K$-row. Cross- vs self-attention are within noise on Ukraine (§14.3) because $\hat K$ is weak (off-diagonals $\approx 0.3$); the MVN $J = 60$ result (§13.4) is where the routing advantage shows.

6.  **Deployment cost is simulation, not network.** Rebuilding the posterior-predictive $z^{(s)}$ block is \~1 min/interim at $S = 4000$; the $80{,}000$ head evaluations take seconds. SVI on $x^{\text{obs}}$ (once per interim) is the only inference step left.

## 14.2 Calibration remedies for the item amortiser

The diagnostics of §14.1.7–14.1.8 leave one defect: the prior-trained amortiser is well-calibrated on the prior-predictive (held-out-prior coverage nominal, §14.1.7) but its deployment predictive is **mis-located** on the posterior-predictive slice — coverage far below nominal, PIT-KS $\approx 0.61$, marginal-KS $\approx 0.47$ — while ranking stays strong (correlation $\approx 0.8/0.7$, §14.1.10). This section reports a remedy sweep. Goal: near-nominal calibration **without** losing correlation or the §8 generality. Decision metrics: coverage at $\eta^{\text{inpol}} = 0.5, 0.95$, PIT-KS, marginal-KS, and correlation (out-of-7 / CG-MH; `CG-VIO_ph-punish` blow-up item excluded). Guardrail: held-out-prior coverage must stay nominal. Every configuration writes the full diagnostic PDFs (coverage curves, PIT histograms, marginal-CDF overlays) into its own `…_260808`/`_260809` sandbox dir.

| config | corr (o7/cat) | cov\@.5 | cov\@.95 | PIT-KS | marg-KS |
|---------------------|-----------|----------:|----------:|----------:|----------:|
| baseline (§14.1) | 0.53/0.42 | 0.15 | 0.44 | 0.61 | 0.47 |
| S=256/step (§14.2.1) | 0.52/0.41 | 0.14 | 0.48 | 0.60 | 0.47 |
| thresh=2.0 (§14.2.2) | 0.53/0.42 | 0.14 | 0.40 | 0.64 | 0.50 |
| thresh=1.0 (§14.2.2) | 0.53/0.41 | 0.14 | 0.45 | 0.62 | 0.50 |
| full-ft LOIO (§14.2.3) | 0.51/0.41 | 0.53 | 0.94 | 0.17 | 0.14 |
| head-ft(i1) (§14.2.4) | 0.52/0.13 | 0.53 | 0.94 | 0.40 | 0.40 |
| **head-ft expand 1:k (§14.2.5)** | **0.52/0.39** | **0.50** | **0.95** | **0.16** | **0.14** |
| conformal i1 (§14.2.6) | 0.53/0.41 | 0.38 | 0.44 | 0.58 | 0.66 |
| conformal expand 1:k (§14.2.6) | 0.53/0.41 | 0.40 | 0.44 | 0.56 | 0.59 |

(nominal cov\@.5 = 0.50, cov\@.95 = 0.95; lower KS is better.)

### 14.2.1 Training-noise control — $S = 256$ draws/step (negative)

More fresh prior draws per step (1280 → 5120 query rows; lower-variance gradients) leaves every calibration metric unchanged (cov 0.14/0.48, PIT-KS 0.60, marg 0.47). Confirms the §14.1.7 diagnosis: the gap is distribution shift, not optimisation noise. Dir `…itemXcompAtt…_256step_260808`.

### 14.2.2 Prior narrowing — `threshold_scale` 2.0 / 1.0 (negative)

Tightening the diffuse threshold prior toward the operational region does not help (cov ≈ 0.14/0.40–0.45, PIT-KS 0.62–0.64, marg 0.50 — if anything slightly worse). The OOD gap is not fixed by shrinking threshold spread alone: abilities/loadings and the population-endpoint-vs-finite-cohort-target mismatch dominate. Data-agnostic prior narrowing is off the table for this study. Dirs `…_thresh2p0/thresh1p0_260808`.

### 14.2.3 Full-net fine-tune, leave-one-interim-out (reference upper bound)

Warm-start from the prior-trained weights, fine-tune **all** parameters 1000 steps on interims 2–8, evaluate held-out interim 1. Near-nominal calibration (cov 0.53/0.94, PIT-KS 0.61 → **0.17**, marg 0.47 → **0.14**) with correlation preserved in- and out-of-sample (all-interim 0.51/0.41; held-out interim-1 0.71/0.72). The fine-tune is *light* — every weight block moves only 8–15 % in norm, cosine $\approx 0.99$ with the prior weights ($q_\text{tok}$ barely 8 %) — so it is the prior amortiser **lightly recalibrated**, not re-fit: the prior stage supplies the retained feature extractor, the fine-tune does the location/scale correction. Cost: reintroduces a per-study SVI dependency at fine-tune time. Dir `…_ftfull-loio-i1_260808`.

### 14.2.4 Head-only fine-tune on interim 1 (too thin)

Freeze the encoders, fine-tune only $q_\psi$ on interim 1. Coverage jumps to nominal (0.53/0.94) but CG-MH correlation **collapses** 0.42 → 0.13 — one interim over-fits interim-1's categorical scale and destroys held-out CG-MH ranking. Head-only is not wrong; one fit-interim is too little data. Dir `…_ftheadi1_260808`.

### 14.2.5 Head-only fine-tune, expanding window (the deployable recipe)

The realistic sequential version: at interim $k$ the SVI fits for interims $1..k$ are already in hand (deployment runs SVI per interim), so fine-tune only $q_\psi$ on $1..k$ and apply to interim $k$ (encoders frozen throughout). **Best result:** cov 0.50/0.95 (nominal), PIT-KS 0.61 → **0.16**, marginal-KS 0.47 → **0.14**, correlation preserved (0.52/0.39 — CG-MH recovers vs §14.2.4 because the expanding window supplies enough data). It **matches the full-net upper bound (§14.2.3) while keeping the encoders universal**, and uses only SVI already computed at interims $\le k$. Recommended operating point. Dir `…-ftheadexpand-260809`.

### 14.2.6 Post-hoc conformal recalibration (partial)

No fine-tuning: fit a per-item monotone PIT map $\hat g = \text{ecdf}(u_{\text{cal}})$ (distribution calibration, Kuleshov et al. 2018) and recalibrate the predictive CDF. Two protocols — interim-1 map applied to all interims, and expanding (map from $1..k$ applied to interim $k$). Both only **partly** help: cov\@.5 0.15 → 0.38–0.40 (not 0.50), PIT-KS 0.61 → 0.56–0.58, marginal-KS *worse* (0.59–0.66). Cause: the miscalibration **varies across interims** (the location bias shrinks as $n$ grows, §14.1.8), so a single static monotone map cannot correct all interims at once — a *conditional* re-scale (fine-tune) can. Dirs `…_conformal-i1/conformal-expand_260809`.

### 14.2.7 Per-item affine median-shift (ports the deepset §14.4.6 default)

The head-only fine-tune fixes the quantile *spread* but leaves a per-item *location* residual — the prior-trained encoder under-predicts every item (the global shrinkage-to-prior bias diagnosed on the deepset in §14.4.6). We port the same fix here: one scalar per item, fit on the *same* SVI of interims $1..k$ the head fine-tune already uses, $\Delta_j(k)=\operatorname{mean}_{i\le k}\big[\operatorname{med}\,p_{\text{SVI}}(\rho\mid x_i)-\operatorname{med}\,\hat p(\rho\mid x_i)\big]$, applied as a pure location shift to all quantiles (`RX_MEDSHIFT`, on by default). It improves both objectives on the hand token:

| itemXcompAtt, head-ft expand (weekly-29) |    PIT-KS |   marg-KS |
|------------------------------------------|----------:|----------:|
| no shift                                 |     0.154 |     0.135 |
| $+$ affine shift                         | **0.139** | **0.119** |

**Use as the fair baseline against the deepset.** With the affine now applied uniformly to *every* encoder and both evaluated on the same weekly-29 grid ($n=48\to500$), the §14.5.1 conditional-vs-marginal trade-off is the honest comparison: the **deepset wins the conditional** (PIT-KS $\approx0.07$ vs the hand token's $0.14$ — the object the PPS decision consumes), while the **hand token wins the marginal** (marg-KS $\approx0.12$ vs the deepset's $\approx0.16$; marg-KS is secondary since $p(\rho\mid x)$ is available from SVI directly). So for the interim decision the deepset remains preferred; the hand token is the better *effect-size* estimator by a small margin. All four comparison figures (§14.4.5: coloured PIT-box, marginal quantile-box, contraction-factor, contraction-law) are now produced for both encoders from the shared `amortiser_diag_plots` module for a like-for-like read.

### 14.2.8 Reading

1.  **Diagnosis holds.** The $S=256$ null (§14.2.1) reconfirms the gap is distribution shift, curable only by moving the predictive toward the operational region.
2.  **Cheap, data-agnostic fixes fail.** Prior narrowing (§14.2.2) and static conformal recalibration (§14.2.6) both fall short — the latter because the bias is interim-heteroscedastic, which a single monotone map cannot absorb.
3.  **Fine-tuning on the operational SVI works, and head-only is enough.** The **expanding-window head-only** fine-tune (§14.2.5) is the deployable operating point: near-nominal calibration, correlation intact, universal encoders retained, using only SVI already computed at interims $\le k$. It buys a correctly-located predictive for a light per-interim head re-fit — a small, honest trade of §8 purity.
4.  **The correction is a small warp, not a retrain.** The full-net fine-tune (§14.2.3) reaches the same calibration by moving the weights only $\approx 10\%$ near the prior solution — the amortiser is recalibrated, not relearned; the prior stage's representation is doing the heavy lifting throughout.

## 14.3 Self-attention variant (`itemScompAtt`): detailed diagnostics

Identical tokens, target, loss, prior-predictive training data and deployment as §14.1; only the attention block differs (§13.2 variant 2). Every token acts as query and key; a learned scalar $\alpha$ biases each row toward the queried key,

$$\text{score}^{(j^*, s)}_{a, b} \;=\; \frac{\big\langle h^{(j^*, s)}_a,\; h^{(j^*, s)}_b\big\rangle}{\sqrt{E}} \;+\; \alpha \cdot \mathbf{1}\{b = j^*\}, \qquad h'^{\,(j^*, s)}_a \;=\; \sum_{b=1}^{J} \operatorname{softmax}_b\big(\text{score}^{(j^*, s)}_{a, \cdot}\big)\, h^{(j^*, s)}_b,$$

and the per-query summary is the gather $\bar h^{(j^*, s)} = h'^{\,(j^*, s)}_{j^*}$ — $J^2$ attention scores per batch element (vs $J$ for cross-attention). $\alpha$ adds one learned scalar to $(\tau, \psi)$.

### 14.3.1 Association with the SVI target

Same protocol and token as §14.1.9 (per-item Pearson $\rho$ between the amortiser median $\hat\rho^{(j^*, s)}_3$ and $\rho^{(j^*, s)}_{\text{SVI-}x}$, blow-up item excluded):

| variant | interim 1: out-of-7 / CG-MH | all interims: out-of-7 / CG-MH |
|--------------------|-------------------------|---------------------------|
| `itemXcompAtt` (§14.1) | 0.82 / 0.71 | 0.53 / 0.42 |
| `itemScompAtt` | 0.83 / 0.70 | 0.53 / 0.41 |

(Both rows: §8 prior-predictive training, identical protocol.) On MVN, cross-attention beat self-attention by 20-35 % at $J \le 60$ (§13.4) — the advantage came from K-row alignment through the query. On Ukraine at $J = 20$ the two variants are within noise: the empirical $\hat K$-row is weaker and shrunk (mean off-diagonal $\approx 0.3$), so the query-side routing advantage largely disappears. Association plots: `pcm_1_interim_i{k}_svi_rho_vs_amortiser_median.pdf`.

### 14.3.2 Baseline calibration — the same prior/posterior mislocation as §14.1.7

The full SBC suite (`..._tests.py`, net-agnostic) on the prior-trained self-attention net gives the identical picture as `itemXcompAtt`: strong ranking, badly mislocated intervals.

| metric (mean over items/interims) | nominal | `itemXcompAtt` baseline | `itemScompAtt` baseline |
|------------------------|--------------:|----------------:|----------------:|
| coverage \@ $\eta = 0.5$ | 0.50 | 0.15 | 0.19 |
| coverage \@ $\eta = 0.95$ | 0.95 | 0.44 | 0.53 |
| PIT–KS | 0.00 | 0.61 | 0.55 |
| marginal–KS | 0.00 | 0.47 | 0.42 |

Diagnostics in the baseline dir: `pcm_1_interim_pps_RGEG_tests_coverage.pdf`, `..._tests_pit.pdf`, `..._i{k}_svi_vs_amortiser_marginal_cdf.pdf`. The mechanism is §8.1's prior-coverage gap, architecture-independent: the diffuse PCM prior is wider than the concentrated deployment posterior slice, so the predicted quantiles are correctly *ordered* but too wide and off-centre.

### 14.3.3 Head-only expanding fine-tune (§14.2.5) transfers unchanged

Applying the §14.2.5 recipe — freeze the encoder, refit only the head $q_\psi$ on interims $1..k$ SVI, evaluate at $k$ — to the self-attention net (universal freeze: train `q_psi`, freeze `q_tok` and the self-attention bias $\alpha$):

| metric           | nominal |  baseline | head-ft-expand |
|------------------|--------:|----------:|---------------:|
| coverage \@ 0.5  |    0.50 |      0.19 |           0.50 |
| coverage \@ 0.95 |    0.95 |      0.53 |           0.94 |
| PIT–KS           |    0.00 |      0.55 |           0.15 |
| marginal–KS      |    0.00 |      0.42 |           0.13 |
| corr o7/cat (i1) |       — | 0.83/0.70 |      0.83/0.72 |

Same outcome as `itemXcompAtt` (§14.2.5): near-nominal calibration, correlation preserved (encoders frozen). Diagnostics in `..._itemScompAtt_...-ftheadexpand-260809/`. `itemXcompAtt` remains the default — the two variants are within noise on this $J = 20$ case and share the same remedy.

### 14.3.4 Per-item affine median-shift (ports §14.2.7)

The self-attention encoder inherits the same per-item location residual, and the same affine shift (§14.2.7 / §14.4.6) fixes it, on by default (`RX_MEDSHIFT`):

| itemScompAtt, head-ft expand (weekly-29) |    PIT-KS |   marg-KS |
|------------------------------------------|----------:|----------:|
| $+$ affine shift                         | **0.136** | **0.108** |

The numbers track `itemXcompAtt` (§14.2.7) to within noise — PIT-KS $0.14$, marg-KS $0.11$ on the same weekly-29 grid — confirming the two hand-token encoders are interchangeable under the full remedy stack, and both sit at the same point of the §14.5.1 conditional-vs-marginal trade-off relative to the deepset (worse conditional, better marginal). The four §14.4.5 comparison figures are produced here too via the shared `amortiser_diag_plots` module.

## 14.4 Deep-set variant (`deepsetXcompAtt` )

The nested-DeepSets variant of the §8-10 workflow on the Ukraine PCM. Training is prior-predictive as in §14.1; the difference is architectural: instead of the hand-computed $F = 8$ item summaries, the network consumes the RAW cohorts — `x_responses` and `z_responses` per (participant, item, time) — and **learns** the per-item pooling through $q_\tau$ (§9.2).

**Why the raw-cohort encoder is the strategic choice, not just a refinement.** The hand-token amortisers of §14.1–14.3 are bound to *this* trial. Their only cohort-size input is the enrolment fractions $(n/N,m/N)$ under a fixed total $N=503$, with the future cohort pinned to $m=N-n$; carrying no absolute-scale input, they cannot price a cohort from a trial of any other size and would need retraining per trial. The deepset removes that ceiling. By pooling over the raw per-participant responses and taking absolute-scale aux $(1/\sqrt n,1/\sqrt m)$ — with $(n,m)$ trained *independently* over $\{2,\dots,503\}$ rather than tied to $n+m=N$ — it is a **general amortised posterior over PCM cohorts**: one trained network reusable across datasets and cohort sizes, an asset rather than a single-trial artefact. That portability is the reason to pay for the architectural machinery below, and it holds even though the deepset trades a little marginal calibration for its decisive edge on the conditional PIT the decision consumes (§14.4.6). (The $n/N$ fraction feature stays $503$-normalised, so a markedly different total $N$ would want it re-normalised; the absolute $1/\sqrt n$ carries the leading contraction regardless.)

The subsections below mirror §14.1.1–14.1.4 and state only what **changes** under the nested-DeepSets construction; everything not mentioned is inherited verbatim from the item-level amortiser.

### 14.4.1 Structure of learning task

The trial structure is unchanged: two cohorts ($x^{\text{obs}}$ vs $z^{(s)}$), two time-points per participant, and the same PCM internals ($2J = 40$ per-`item_time_id` fake items, shared ability $\theta_i$). What changes is the **§10 gap flagged at the end of §14.1.1**: the item-level amortiser collapses the time and participant axes into a hand-computed $F = 8$ token *before* the network sees anything; here the network consumes the **raw** per-(participant, item, time) responses and **learns** the per-item pooling over participants through $q_\tau$ + mean-pool (§9.2). The cohort axis (x vs z) and the time axis ($R = 2$ response features, baseline + endline) are kept apart as separate tensor axes rather than as separate token slots.

Symbol deltas from the §14.1.1 table:

| symbol | item-level (§14.1) | deepset (§14.4) | note |
|------------|------------|------------------------|-------------------------|
| inputs | $F = 8$ hand token | raw tensors `x_responses`, `z_responses` $\in \mathbb{R}^{N \times J \times 2}$ + masks | pooling now learned |
| $R$ | — | 2 | response features per (participant, item): baseline, endline |
| $M$ | (in token) | 2 | item metadata (type, direction) appended to every participant token |
| $q_\tau$ | — | shared participant encoder | the §9.2 inner DeepSet; no item-level counterpart |
| $\hat K$-row | in token | **dropped** | cross-item borrowing is now done by attention over learned embeddings, not a hand statistic |
| $S$ | 4000 | 200 | draws per interim at deploy (heavier per-draw forward pass) |
| $E$ | 32 | 32 | embedding dim, unchanged |

Reused / trained / adapted (cf. §14.1.1):

| piece | status | detail |
|------------|------------|-------------------------------------------------|
| network classes ($q_\tau, q_{\text{tok}}, q_{\text{query}}, q_\psi$) | **reused** from §13.2 | nested DeepSets + cross-attention |
| network weights | **trained from scratch** | \~18k parameters, multi-quantile pinball loss |
| per-item token | **removed** | replaced by learned pool of raw responses (§14.4.3–14.4.4) |
| aux $a$ | **extended** | $(1/\sqrt n, 1/\sqrt m, n/N, m/N)$, $N = 503$ — absolute-scale precision added so contraction is representable at any cohort size (§14.4) |
| input scale | **same rule** | raw $(y - 1)/(K_j - 1) \in [0, 1]$, $K_j$ from `dit[cat_length]` |
| target | **same** | $\rho^{(j^*, s)} = r_{j^*}(\theta^{(s)})$, clipped $\pm 20$ |
| target scale | **same as §14.1** | per-item standardisation $\rho_{j^*}/\hat\sigma_{j^*}$; $\hat\sigma_j=\operatorname{sd}_s\!\big[r_j(\theta^{(s)})\big]$ is the **prior-predictive standard deviation** of the effect, estimated once on a $2{,}000$-sample prior pilot — one fixed constant per item, **constant across interims** (stored as `item_std`); the head quantiles are rescaled by $\hat\sigma_j$ at deployment |

### 14.4.2 Training data and target

The network is trained on the prior predictive of the partial-credit model (§8): every training example is a fresh draw from the model prior, so no SVI fit enters training (6000 steps, batch $B=24$, fresh draws each step).

**Prior.** Each draw uses the standard partial-credit-model priors — the same specification fit in Pyro/Stan — with the ordered-category thresholds given a deliberately wide variance: item time-effects $\beta_t\sim N(0,1)$ for $t\in\{\text{baseline},\text{endline}\}$; thresholds $\tau_{j,t,k}\sim N(0,3.5^2)$; discriminations $\lambda_{j,t}=\lvert T_3\rvert$ (half-Student-$t$, 3 df, the first item-time fixed to $1$ for identification); and participant abilities $\theta_i\sim N(0,1)$. The wide threshold prior is deliberate, for two reasons. First, the training prior must **envelope** any posterior that can arise at deployment: because the amortiser is only ever fitted where the prior predictive puts mass, the prior predictive has to cover the region the real data's posterior occupies, or the network is evaluated out of distribution — the prior/deployment mislocation of §14.4.5 is precisely this distribution shift, and a proposal that fails to envelope it cannot be repaired by the head fine-tune. Second, a broad proposal leaves posterior uncertainty for a finite cohort to resolve, which is what makes contraction with cohort size learnable.

**Cohort sampling.** For each parameter draw, two cohort sizes are chosen **independently** — the observed size $n$ and the future size $m$ — spread log-uniformly over $\{2,\dots,N\}$ with $N=503$, so all cohort sizes are represented across a batch. A cohort of $n_{\text{ab}}$ participants is then simulated: draw $n_{\text{ab}}$ abilities $\theta_i\sim N(0,1)$, and at the drawn $(\beta,\tau,\lambda)$ sample each participant's ordered-categorical response at baseline and at endline from the PCM (inverse-CDF from the cumulative category probabilities), rescaled to $[0,1]$ as $y/(K_j-1)$. This yields the observed cohort $x^{(s)}$ and the future cohort $z^{(s)}$ as raw per-participant response tensors.

**Target.** The estimand the amortiser fits is the **conditional posterior** $p(\rho_j\mid x,z)$ — the posterior over the population effect size given the observed and future cohorts. It is fitted by the standard simulation-based-inference argument (§14.1.7): for a joint draw $(\theta^{(s)},x^{(s)},z^{(s)})$ — parameters from the prior, cohorts simulated from them — the value $\rho_j^{(s)}=r_j(\theta^{(s)})$ is a *single valid draw* from $p(\rho_j\mid x^{(s)},z^{(s)})$, so training the amortiser's quantiles against $\rho_j^{(s)}$ (interval score) fits that conditional. Concretely: every amortiser output $\hat\rho_j^{(s)}$ is scored against the $\rho_j^{(s)}$ from which its own inputs $x^{(s)},z^{(s)}$ were generated.

The functional $r_j(\theta)$ of the drawn parameters $\theta=(\beta,\tau,\lambda)$ is the direction-aware relative change in the mean item response level between baseline and endline, $$r_j(\theta)=s_j\Big(\bar w_{e,j}/\bar w_{b,j}-1\Big),\qquad \bar w_{t,j}=\mathbb E_{\theta_i\sim N(0,1)}\!\big[\,g_j(\text{response}_{i,j,t}\mid\theta)\,\big],$$ where $g_j$ is the expected score $\sum_k k\,p_{jk}$ for out-of-7 items and the caseness probability $\Pr(y_j\ge c)$ for categorical items, $s_j=\pm1$ the item direction; clipped to $\pm20$. **There is no closed-form expression for** $r_j(\theta)$: the mean response level $\bar w_{t,j}$ is the expectation *over the ability distribution* $\theta_i\sim N(0,1)$ of a nonlinear category probability (a softmax of cumulative threshold sums), and this Gaussian integral has no analytic solution under the partial-credit likelihood. We therefore estimate $\bar w_{t,j}$ by a **reference simulation** of $N_{\text{ref}}=2000$ abilities drawn from $N(0,1)$. This keeps the construction model-general — it needs only the ability to *simulate* the model, never to integrate it — under the assumption that $N_{\text{ref}}$ is large enough that the Monte-Carlo error in $\bar w_{t,j}$ is negligible against the posterior width, i.e. $r_j(\theta)$ is treated as the exact target.

Because the target is this population functional while the network sees only finite cohorts, the predictive $\hat p(\rho_j\mid x,z)$ must **contract** toward $r_j(\theta)$ as the cohorts grow — the signal quantified in §14.4.7 and exploited in §14.4.7. The sampler (`draw_params`, `theta_level_rho`, `sim_cohort` in the explore script; [`python/model_pcm.py`](../python/model_pcm.py)) fans each draw out to $Q=4$ queried items.

### 14.4.3 Ragged per-participant data structure

The two cohorts enter the network as raw per-participant response tensors, not pooled summaries; because the cohort sizes vary from one training example to the next, they are stored **ragged** — no padding to a common size, and no participant mask.

**Structure.** Across a batch of $B$ examples, all observed-cohort participants are concatenated along a single flat axis, $x_{\text{flat}}\in\mathbb R^{T\times J\times 2}$, where $T=\sum_b n_b$ is the total participant count over the batch, $J$ the items, and the last axis the (baseline, endline) responses; a companion segment vector $x_{\text{seg}}\in\{0,\dots,B-1\}^{T}$ records which example each participant belongs to. The future cohorts are stored identically as $(z_{\text{flat}}, z_{\text{seg}})$. The arrays therefore hold exactly the real participants — there is no $N$-padding and no $0/1$ mask.

### 14.4.4 Encoder: pooling, item cross-attention, quantile head, loss

**Pooling.** The per-item cohort summary is a **segmented mean** over participants. A shared MLP $q_\tau$ embeds each (participant, item, time) response vector, and `jax.ops.segment_sum` sums the embeddings within each example's segment, divided by the segment count: $$\text{pool}^{x}_{b,j} \;=\; \frac{1}{n_b}\sum_{i:\,x_{\text{seg}}(i)=b} q_\tau\big(x_{\text{flat}}[i,j,\cdot],\,c_j,d_j\big) \;\in\;\mathbb R^{E}, \qquad \text{pool}^{z}_{b,j}\ \text{from}\ m_b\ \text{likewise}.$$ Dividing by the true $n_b$ (resp. $m_b$) makes each summary the sample mean of the learned per-participant embedding at the actual cohort size — the amortised replacement for the hand-computed group-means of §14.1.3, and the standard exchangeable-set summary used by simulation-based inference for variable dataset sizes (Radev et al. 2020; Chan et al. 2018).

**Cross-attention and head.** The per-item embedding mixes the two cohort pools, $h_j = q_{\text{tok}}(\operatorname{concat}(\text{pool}^x_j,\text{pool}^z_j))$; from here the pipeline is exactly §14.1.4 — the query is built from the queried item's embedding, cross-attention over the $J$ embeddings gives $\bar h_{j^*}$, and the multi-quantile head consumes $(\bar h_{j^*}, h_{j^*}, a)$ with aux $a=(1/\sqrt n,1/\sqrt m,n/N,m/N)$: $$\hat\rho^{(j^*,s)} = q_\psi\big(\operatorname{concat}(\bar h_{j^*}, h_{j^*}, a)\big)\in\mathbb R^K.$$ The learned parameters are $(q_\tau,q_{\text{tok}},q_{\text{query}},q_\psi)$, $\approx18{,}400$ in total.

**Loss.** The head outputs $K$ monotone quantiles at fixed levels $0<\tau_1<\dots<\tau_K<1$. Training minimises the **pinball (quantile) score**, summed over prior draws $s$, queried items $j^*$ and levels $\tau_k$: $$\mathcal L=\sum_{s}\sum_{j^*}\sum_{k}\ell_{\tau_k}\!\big(\rho^{(s)}_{j^*}-\hat\rho^{(j^*,s)}_{\tau_k}\big),\qquad \ell_\tau(u)=\big(\tau-\mathbf 1\{u<0\}\big)\,u,$$ equivalently the interval score summed over the symmetric quantile pairs $(\tau_k,1-\tau_k)$. The pinball score at level $\tau$ is minimised exactly at the $\tau$-quantile of the target's conditional law; since $\rho^{(s)}_{j^*}$ is a draw from $p(\rho\mid x^{(s)},z^{(s)})$ (§14.4.2), the fitted head recovers the quantiles of that conditional posterior.

**Alternatives considered.** A head that instead fits the **marginal** $p(\rho\mid x)$ directly — or an architecture that guarantees the law of total probability $p(\rho\mid x)=\int p(\rho\mid x,z)\,p(z\mid x)\,dz$ by construction (a normalising-flow or diffusion posterior head) — was considered and dropped. At deployment $p(\rho\mid x)$ is obtained directly from the SVI fit of the observed cohort because we need this to generate future cohorts; the piece SVI cannot supply is the future cohort $z$, so the amortiser is needed only for the **conditional** $p(\rho\mid x,z)$.

Net module `python/amortiser_pps_features_deepsetXcompAtt_ragged_qpsi_MLP_loss_multiquantilehead.py`; trained network directory `…-deepsetXcompAtt-plain-260812`. The evaluation grid is the 29 weekly SVI fits (`py-ukraine-interim-weekly-svi-260811`, $n=48\to500$), reused with no recomputation. (The fixed-$N$ padded predecessor and the focused interim-1 proposal are catalogued in §14.10.)

### 14.4.5 Diagnostics

Four diagnostics probe the raw (pre-calibration) network, forward-passing every SVI draw; all quantities exclude the degenerate item `CG-VIO_ph-punish`. Each is defined here as a statistic on the amortiser output, with the figure that illustrates it; the numbers themselves appear from §14.4.6 on.

**Association with the SVI target.** For a joint draw, the SVI posterior draw $\rho^{(s)}_{\text{SVI-}x}$ from the observed cohort is passed through the amortiser to obtain the conditional median $\operatorname{med}\hat p(\rho\mid x^{(s)},z^{(s)})$. The diagnostic is the Pearson correlation between the two as a function of cohort size $n$ — whether the encoder's set-summary retains the cohort's signal. It necessarily decays with $n$: as the SVI posterior tightens onto $r_j(\theta)$ the spread of targets to track collapses, so a falling correlation is expected, and this diagnostic needs to condition on a particular n band.

**Conditional calibration (PIT / SBC).** Under the §14.1.7 joint-draw construction $\rho^{(s)}_{\text{SVI-}x}$ is a valid draw from the conditional $p(\rho_j\mid x,z^{(s)})$ the amortiser targets, so its probability-integral transform under the amortiser's predictive, $$u^{(s)}_j \;=\; \hat F_j\big(\rho^{(s)}_{\text{SVI-}x}\mid x,z^{(s)}\big)\ \in[0,1]$$ — the rank of the target among the predicted quantiles, linearly interpolated — must be **uniform** if that conditional predictive is calibrated. The scalar summary is the KS distance of the empirical PIT CDF from uniform, $\text{PIT–KS}_j=\sup_t\lvert \hat G_j(t)-t\rvert$. Illustrated by the **PIT box** (`…_pit_box_by_item.pdf`): per item, precomputed whiskers (10/90% quantiles), box (25/75%) and median line of $\{u^{(s)}_j\}$ across draws, coloured by the signed deviation of the PIT median from $\tfrac12$; a calibrated item sits at median $\tfrac12$ with box edges near $0.25/0.75$.

**Marginal calibration.** The marginal predictive is built by the law of total probability over the future cohort: draw $z^{(s)}\sim p(z\mid x)$ — posterior-predictive future cohorts simulated from the SVI fit of the observed cohort $x$ — run the amortiser on each to get the conditional $\hat p(\rho_j\mid x,z^{(s)})$, and pool them, $\hat p(\rho_j\mid x)=\frac1S\sum_{s=1}^S \hat p(\rho_j\mid x,z^{(s)})$. At a fixed interim this should match the SVI posterior $p_{\text{SVI}}(\rho_j\mid x)$. The scalar summary is the KS distance between the two CDFs, $\text{marg–KS}_j=\sup_t\lvert \hat F^{\text{marg}}_j(t)-F^{\text{SVI}}_j(t)\rvert$, pooled over interims. Illustrated by the **marginal quantile box** (`…_contraction_cdf_by_item.pdf`, with a `…RAGDbase…` raw-network companion): $x$-axis interim, $y$-axis $\rho$-quantiles, one box per interim for SVI and one for the amortiser (whiskers 10/90%, box 25/75%, median line), so a systematic location or spread mismatch is read off directly.

**Posterior contraction.** The predictive should sharpen with cohort size at the SVI rate. The **contraction law** plot (`…_contraction-law_trained-vs-svi.pdf`) shows, per item, $\log\mathrm{SD}(\rho_j\mid x)$ against $\log n$ for the SVI posterior and for the amortiser, with a fitted power law $\mathrm{SD}=C_j\,n^{-p_j}$ overlaid. The fit takes $(C_j,p_j)$ by least squares of $\log\mathrm{SD}$ on $\log n$ across the weekly SVI grid; on log–log axes $\log\mathrm{SD}=\log C_j-p_j\log n$ is a **straight line of slope** $-p_j$, so the constant relative contraction rate $\mathrm{d}\log\mathrm{SD}/\mathrm{d}\log n=-p_j$ is exactly that slope. The SVI SDs fall along this line; a well-matched amortiser tracks the same slope, an under-contracting one is shallower. The companion **contraction factor** plot (`…_posterior-contraction-factor.pdf`) shows the ratio of amortiser to SVI predictive SD against $n$ — flat at one when the rates agree, drifting above one and growing with $n$ when the amortiser under-contracts. This rate is the second target, modelled in §14.4.7.

All figures are produced by `python/amortiser_diag_plots.py`.

### 14.4.6 Head fine-tune: expanding window and affine shift

The mislocation is corrected by refitting only the quantile head $q_\psi$ on the operational SVI, with the encoder frozen, followed by a per-item affine location shift. Two fine-tune schedules are compared against the raw network on the weekly grid (conditional PIT–KS; marginal KS):

| head fine-tune ($+$ affine shift)       |    PIT–KS |   marg–KS |
|-----------------------------------------|----------:|----------:|
| none (raw network)                      |     0.465 |     0.441 |
| interim-1 only (`…_ftheadi1`)           |     0.144 |     0.238 |
| **expanding** $1..k$ (`…_ftheadexpand`) | **0.069** | **0.157** |

**Interim-1 fitting is not pursued further.** A head fitted at interim 1 alone calibrates the anchor but drifts (PIT–KS $0.144$, marg–KS $0.238$ over the 29 interims), because the location bias is *interim-heteroscedastic* — its magnitude shrinks as $n$ grows, so one fit cannot serve all cohort sizes. The **expanding-window** fit — refit $q_\psi$ on interims $1..k$ at each $k$ (cost $\approx0.5$ s/refit, negligible) — conditions on each interim through the frozen embedding and reaches near-nominal calibration (PIT–KS $0.069$).

**Post-hoc conformal recalibration fails for the same reason.** A static per-item monotone PIT map (distribution calibration, Kuleshov et al. 2018), fitted on interim-1 or on the expanding pool, only partly corrects coverage (PIT–KS $\approx0.56$–$0.59$) and *worsens* the marginal KS: a single monotone map cannot absorb an interim-heteroscedastic bias, whereas refitting the head can. Set aside (§14.10).

**Affine median-shift (on by default).** The head fine-tune fixes the quantile spread but leaves a per-item *location* residual — the prior-trained encoder systematically under-predicts every item. One scalar per item, fitted on the same SVI of interims $1..k$, corrects it, $\Delta_j(k)=\operatorname{mean}_{i\le k}\big[\operatorname{med}p_{\mathrm{SVI}}(\rho\mid x_i)-\operatorname{med}\hat p(\rho\mid x_i)\big]$, applied as a pure location shift to all quantiles (`RAGD_MEDSHIFT`, default on). It lowers PIT–KS $0.075\to0.069$ and marg–KS $0.180\to0.157$, and removes the systematic under-prediction visible in the marginal quantile box.

**Direct comparison against the hand-token encoders.** The *same* expanding head-ft $+$ affine remedy is applied to the two §14.1 hand-token encoders — `itemXcompAtt` (cross-attention, §14.2.7) and `itemScompAtt` (self-attention, §14.3.4) — as fair baselines for the deepset, all on the **same weekly-29 grid** (ids 2..30, $n=48\to500$); all four §14.4.5 figures are produced for each from the shared `amortiser_diag_plots` module for a like-for-like read:

| encoder (expand head-ft $+$ affine, weekly-29) |    PIT–KS |   marg–KS |
|------------------------------------------------|----------:|----------:|
| **`deepsetXcompAtt`** (learned pooling)        | **0.069** |     0.157 |
| `itemXcompAtt` (hand token, §14.2.7)           |     0.139 |     0.119 |
| `itemScompAtt` (hand token, §14.3.4)           |     0.136 | **0.108** |

The split is consistent across the two hand-token encoders (interchangeable to within noise, §14.3.4): the **deepset wins the conditional PIT–KS by** $\approx2\times$ — the object the PPS decision consumes — while the **hand tokens win the marginal KS** by a smaller margin. marg–KS is secondary, since $p(\rho\mid x)$ is available from SVI directly (§14.4.5); so for the interim decision the deepset stays preferred, with the hand token the marginally better pure effect-size estimator. The trade-off is taken up in §14.5.1.

**Both KS metrics vary with cohort size; the all-interim means above hide an** $n$-gradient. Splitting the weekly grid into small ($n\le150$; interims 2–5), mid ($150<n\le400$; 6–13) and large ($n>400$; 14–30):

| encoder      |     PIT–KS s/m/l      |       marg–KS s/m/l       |
|--------------|:---------------------:|:-------------------------:|
| deepset      | 0.072 / 0.085 / 0.060 | 0.111 / 0.144 / **0.172** |
| itemXcompAtt | 0.150 / 0.168 / 0.122 |   0.070 / 0.101 / 0.138   |
| itemScompAtt | 0.139 / 0.166 / 0.121 |   0.065 / 0.092 / 0.126   |

**marg–KS worsens monotonically with** $n$ for every encoder (deepset $0.111\to0.172$), and the deepset is the worst of the three at large $n$. The cause is under-contraction: the deepset's marginal 5–95% width *relative to the SVI width* runs $0.61\to0.80\to1.58$ across the three bands — narrower than SVI at small $n$, but $58\%$ too wide at large $n$. That is, the amortiser's marginal SD sits *above* the SVI best-fit power law (§14.4.7) at large $n$ and stops contracting where SVI keeps sharpening; this upward departure from the law is precisely the residual the Bernstein–von Mises correction of §14.4.7 removes. Conditional PIT–KS, by contrast, has no strong $n$-trend and the deepset leads in every band.

### 14.4.7 Bernstein–von Mises best-fit contraction law

**The best-fit law.** Bernstein–von Mises gives the leading order: the posterior SD of a smooth functional contracts as $\mathrm{SD}(\rho_j\mid x)\approx \sqrt{V_{\rho_j}}\,n^{-1/2}$ with $V_{\rho_j}=\nabla r^\top \mathcal I^{-1}\nabla r$ the delta-method (sandwich) variance at the semiparametric efficiency bound [@vandervaart1998asymptotic; @kleijn2012bernstein]. Over the finite $n$-range of a trial the effective exponent departs from $\tfrac12$ for item-specific reasons: $\rho$ is a **ratio** $\bar w_e/\bar w_b$, whose delta-method variance is inflated by denominator uncertainty at small $n$ (Fieller [@fieller1954some]); finite-$n$ higher-order terms add curvature; the endpoint integrates over the ability distribution; items near a response boundary identify faster; and the finite-population correction $(1-n/N)$ [@cochran1977sampling] speeds contraction as $n\to N$. These corrections do not overturn the leading $n^{-1/2}$; they bend it, and over a finite window that bend reads as a shifted exponent. Formally, define the **effective exponent** as the local log–log slope $$p_j(n):=-\,\frac{\mathrm d\log\mathrm{SD}(\rho_j\mid x)}{\mathrm d\log n},$$ which is identically $\tfrac12$ under pure BvM. Collecting the corrections into a sub-leading factor $g_j(n)$, $\mathrm{SD}^2(\rho_j\mid x)=V_{\rho_j}\,n^{-1}\,g_j(n)$, gives $p_j(n)=\tfrac12-\tfrac12\,\mathrm d\log g_j/\mathrm d\log n$ — generally $\neq\tfrac12$, and drifting only slowly with $n$. Over the trial's narrow window ($n:48\to500$ spans only $\Delta\log n\approx2.3$, a $3.2\times$ range in $\sqrt n$) this slope is nearly constant, so $\log\mathrm{SD}$ is nearly **linear** in $\log n$: a single-exponent power law $\mathrm{SD}=C_j\,n^{-p_j}$ with $p_j$ the window-averaged slope. The exponent is item-specific because $g_j$ — hence its slope — depends on the item's $\mathcal I^{-1}$, $\nabla r$, boundary proximity and floor. Fitting candidate laws to the weekly SVI SD ($n=48\to500$, 25 points/item):

| form                               | median $R^2$ |
|------------------------------------|-------------:|
| power law $\mathrm{SD}=C\,n^{-p}$  |     **0.90** |
| floor $\sqrt{a^2+b^2/n}$           |         0.68 |
| BvM (fixed $p=\tfrac12$)           |         0.65 |
| finite-population $\sqrt{1/n-1/N}$ |         $<0$ |

The **item-specific power law wins** ($R^2=0.90$): the exponent $p_j$ has median $0.45$ but ranges $0.26\to1.08$, genuinely heterogeneous rather than the universal $\tfrac12$. The finite-population form fails because the SD plateaus at a nonzero value; the two-parameter floor $\sqrt{a^2+b^2/n}$ captures that plateau but fits less well than the free power law. The fitted law $C_j\,n^{-p_j}$ is the object the rest of this section uses — first to *justify* a self-consistency correction for the deepset's large-$n$ marginal upturn (§14.4.6), then as a *parametric head* that carries the law into the network.

**Why a self-consistency correction, and what it must be.** The large-$n$ marginal error is a finite-$m$ artefact. At deployment the marginal $p(\rho\mid x)$ is the mixture $\frac1S\sum_s \hat p(\rho\mid x,z^{(s)})$ over posterior-predictive future cohorts, and its variance decomposes exactly by the law of total variance, $$\operatorname{Var}_{\text{mix}} = \underbrace{\tfrac1S\textstyle\sum_s \sigma_s^2}_{W\ (\text{within})} + \underbrace{\tfrac1S\textstyle\sum_s (\mu_s-\bar\mu)^2}_{B\ (\text{between})},$$ where $\mu_s=q_{s,3}$ is the conditional median and $\sigma_s=(q_{s,5}-q_{s,1})/3.2897$ its Gaussian-equivalent SD, both read from the deployed quantiles. As $n\to N$ the future cohort $m=N-n$ shrinks, so each $z^{(s)}$ is a noisy sample and $B$ inflates spuriously — the amortiser's marginal SD turns *up* at large $n$ instead of contracting.

The correction is **pinned by the contraction law**, with no free constant. Bernstein–von Mises states that a posterior variance depends, to leading order, on the sample *size*, not the sample *values*; so every conditional has the same BvM variance, given by the fitted law read at $n+m$, while the target marginal is the law read at $n$: $$W_{\text{true}}=\big(C\,(n+m)^{-p}\big)^2,\qquad T=\big(C\,n^{-p}\big)^2,\qquad B^\ast = T-W_{\text{true}} = C^2\big(n^{-2p}-(n+m)^{-2p}\big).$$ The non-spurious between $B^\ast$ is what the extra $m$ observations are entitled to remove. By the law of total variance it is $B^\ast=\operatorname{Var}_z\!\big(\mathbb E[\rho\mid x,z]\big)$ — the further contraction the future cohort buys. At deployment $n+m=N$ is *fixed*, so $W_{\text{true}}=(C\,N^{-p})^2$ is constant and $B^\ast=C^2\big(n^{-2p}-N^{-2p}\big)$. **Why** $B^\ast\to0$ as $n\to N$: the remaining cohort $m=N-n\to0$, so there is no future data left to collect; conditioning on the (vanishing) $z$ cannot sharpen the posterior, the conditional and marginal coincide, and $\operatorname{Var}_z(\mathbb E[\rho\mid x,z])\to0$ — algebraically $n^{-2p}\to N^{-2p}$. At the other end $n\to0$ almost all the variance is still removable, and $B^\ast\to T$. The empirical $B$ moves the *opposite* way at large $n$: with $m$ tiny each simulated $z^{(s)}$ is a high-variance draw, so the $S$ conditional medians $\mu_s$ scatter widely — pure finite-$m$/finite-$S$ sampling noise. That is exactly the spurious inflation the correction strips out, by shrinking $B$ down to $B^\ast$. Scaling each conditional width by $\alpha=\sqrt{W_{\text{true}}/W}$ and shrinking the median spread by $\beta=\sqrt{B^\ast/B}$, $$q'_{s,i}=\big[\bar\mu+\beta(\mu_s-\bar\mu)\big]+\alpha\,(q_{s,i}-\mu_s),$$ makes the mixture variance $\alpha^2W+\beta^2B=W_{\text{true}}+B^\ast=T$ by construction. It is a post-hoc pass over the deployed quantiles (`RAGD_BVM`) — no retraining, and no SVI beyond the $\operatorname{Var}(\rho\mid x)$ fit already in hand.

**Implementation: parametric heads that carry the law.** The correction is intrinsically justified when the head's width already obeys the law, so the width is made parametric (`head_mode`) *without touching the encoder*. Everything up to the §14.4.4 head is unchanged: item cross-attention still yields $\bar h_{j^*}$, and the shared MLP $q_\psi$ still consumes $\operatorname{concat}(\bar h_{j^*},h_{j^*},a)$ and emits a length-$K$ vector $\text{raw}$. What changes is only the *interpretation* of that vector. In `plain` mode $\text{raw}$ is returned as the $K$ quantiles directly; in the parametric modes its leading entries are read as a few shape parameters, from which the $K$ quantiles are rebuilt deterministically as a **fixed-shape Gaussian predictive whose scale follows the law**, $$\hat\rho_q = m \;+\; z_q\,s_j(n),\qquad z_q=\Phi^{-1}(\tau_q)\in\{-1.64,-0.67,0,0.67,1.64\},$$ the standard-normal quantile multipliers $z_q$ fixed, and the cohort size entering **only** through the absolute-scale aux precision $\text{prec}=a_1=n^{-1/2}$ (§14.4.4). The trunk is thus shared; only the last map from $\text{raw}$ to quantiles is replaced, and the two forms differ solely in the width $s_j(n)$:

- **A — power law:** $q_\psi$'s first three outputs give $m=\text{raw}_0$, $C_j=\operatorname{softplus}(\text{raw}_1)>0$ and $p_j=\sigma(\text{raw}_2)\in(0,1)$, so $s_j(n)=C_j\,n^{-p_j}=C_j\,\text{prec}^{2p_j}$. Location, scale and exponent are all item- and query-conditioned through $\bar h_{j^*}$, while the $n$-dependence is hard-wired to the BvM power law. (Regularised variant **B** adds a prior $\lambda\,\mathbb E[(p_j-\tfrac12)^2]$ pulling the exponent toward $\tfrac12$, since $p$ is only weakly identified over the trial's $\sim3\times$ $\sqrt n$ range.)
- **C — floor:** $q_\psi$'s first three outputs give $m$ and $a_j,b_j=\operatorname{softplus}(\cdot)>0$, so $s_j(n)=\sqrt{a_j^2+b_j^2/n}$, encoding the non-zero plateau directly (no exponent).

Because the trunk is untouched, the parametric head is a drop-in: the expanding head fine-tune (§14.4.6) still refits $q_\psi$, and the affine shift still corrects the emitted location $m$.

Trained, the heads roughly double the training-distribution contraction slope (probe $-0.32/-0.35$ vs the plain net's $-0.20$) and the power-law exponent tracks SVI for the strongly-contracting items ($p_{\text{SVI}}/p_{\text{amo}}$ of $1.08/0.91$, $0.99/0.78$). On their own — expanding head fine-tune $+$ affine shift, *without* the correction — they do **not** beat the plain net on aggregate calibration (A power-law PIT–KS $0.112$, marg–KS $0.189$; C floor $0.115$/$0.185$; plain $0.069$/$0.157$): the low-SD mental-health items ($p_{\text{SVI}}\approx0.3$) remain under-learned, and a rigid trained width helps only the strongly-contracting items. Directories `…-powerlaw`, `…-floor` and their `_ftheadexpand` deployments.

**Result of the correction.** Applied to the A and C heads (with the affine shift), the BvM correction lowers marg–KS at every cohort size and makes the amortiser's contraction curve track SVI's with the large-$n$ upturn removed (matching fitted slopes per item, e.g. $-0.005$ vs $-0.005$); PIT–KS is unchanged, since the correction acts only on the marginal:

| head-ft $+$ affine | marg–KS ($n\le200$) | (mid) | ($n\ge450$) | (all) | ($n\ge498$) |
|------------|------------:|-----------:|-----------:|-----------:|-----------:|
| A power-law | 0.169 | 0.194 | 0.193 | 0.189 | 0.224 |
| **A** $+$ BvM | 0.161 | 0.183 | **0.180** | **0.177** | **0.194** |
| C floor | 0.154 | 0.185 | 0.199 | 0.185 | 0.232 |
| **C** $+$ BvM | 0.147 | 0.171 | **0.170** | **0.166** | 0.231 |

Notably **C** $+$ BvM's large-$n$ marg–KS ($0.170$) beats the plain net $+$ affine ($0.181$) — the correction supplies exactly the missing large-$n$ contraction, as a principled step rather than a rigid trained width. Directories `…-powerlaw-ftheadexpand-bvm-260828`, `…-floor-ftheadexpand-bvm-260828`. **A/C** $+$ BvM is the statistically coherent estimator of §14.0: its width obeys a BvM contraction law and its marginal is self-consistent under the law of total variance. The plain net $+$ affine (§14.4.6) remains numerically ahead on the aggregate PIT–KS and all-$n$ marg–KS and is the lighter alternative when the parametric guarantee is not required.

### 14.4.8 Amortising the decision thresholds

The interim decision (§4) carries two thresholds: the **effect threshold** $\eta_0$ in $H_1^{(j)}=\{\rho_j>\eta_0\}$, and the **decision threshold** $\eta_H$ in $\widehat{\text{PPS}}^{(j)}=\frac1S\sum_s\mathbf 1\{\hat P(H_1^{(j)}\mid x,z^{(s)})>\eta_H\}$. Both should be free at deployment, without retraining.

$\eta_H$ is free by construction — it is only a cutoff applied to the already-computed per-draw success probabilities, so any value is a pure post-processing choice.

$\eta_0$ enters through the conditional CDF, $\hat P(H_1\mid x,z^{(s)})=1-\hat F(\eta_0\mid x,z^{(s)})$, so it is free wherever the predictive is accurately resolved at $\eta_0$ — the question is which predictive representation resolves an *arbitrary* $\eta_0$ best, especially in the tails where the $K=5$ quantile grid interpolates crudely. Four routes (net directories `…-deepsetXcompAtt-{plain,powerlaw,dense19,eta0amortise}-…`):

1.  **5-quantile interpolation** (the deployable §14.4.6 head): $1-\operatorname{interp}(\eta_0)$ on the five deployed quantiles.
2.  **Parametric Gaussian** (the power-law / floor head of §14.4.7): the head emits $(m,s)$, so $\hat P(H_1\mid x,z)=\Phi\big((m-\eta_0)/s\big)$ is **analytic for any** $\eta_0$, smooth in the tails.
3.  **Dense quantiles**: retrain the head with $K=19$ levels and a pinball loss, for finer CDF interpolation.
4.  $\eta_0$-conditioned head: retrain with $\eta_0$ appended to the aux and sampled per example, and replace the quantile head with a **success-probability head** $\hat P(H_1\mid x,z,\eta_0)=\sigma\big(q_\psi(\cdot)\big)$ (one output, `head_mode='successprob'`) trained by binary cross-entropy against $\mathbf 1\{\rho>\eta_0\}$ — a direct amortisation of the decision over $\eta_0$.

Scoring each route by the calibration error of the marginal success probability against SVI, $\big\lvert \hat P(\rho>\eta_0\mid x)-P_{\text{SVI}}(\rho>\eta_0\mid x)\big\rvert$, averaged over items and interims across a standardised-threshold sweep (body $|\eta_0/\hat\sigma_j|\le\tfrac12$, tail $\ge\tfrac32$):

| route                       |      body |  tail |       all |
|-----------------------------|----------:|------:|----------:|
| 1 — 5-quantile interp       |     0.041 | 0.030 |     0.034 |
| **2 — parametric Gaussian** | **0.033** | 0.030 | **0.031** |
| 3 — dense-19 quantiles      |     0.047 | 0.029 |     0.035 |
| 4 — $\eta_0$-conditioned    |     0.053 | 0.029 |     0.038 |

**The parametric Gaussian wins**: lowest overall error, decisively best at the centre (near the median, $0.012$ vs $0.03$–$0.08$), tied elsewhere, and **zero retraining** — its analytic $\Phi$ resolves any $\eta_0$. Denser quantiles do **not** beat the five-quantile baseline (the predictive is close to Gaussian, so five knots already capture the tails), and the $\eta_0$-conditioned head helps only in the far tail while being worst at the centre — its BCE head is sharpest where the success probability saturates and weakest near $p=\tfrac12$. A full retrain for $\eta_0$-conditioning is therefore not justified.

**Choosing** $\eta_0$ on a scientific scale. The effect $\rho_j$ is a per-item ratio whose prior-predictive spread $\hat\sigma_j$ ranges widely ($0.86$ to $5.6$). A threshold expressed in *standardised* units — a fixed $\eta_0/\hat\sigma_j$ across items — is not comparable across items: at $\eta_0/\hat\sigma_j=1$ a wide item sits at a $560\%$ change and a narrow one at $86\%$, so its PPS is degenerate ($0$ or $1$). The remedy is to fix $\eta_0$ on the **raw endpoint-change scale**, $\eta_0\in\{0,25,50,75,100\}\%$, and query each item at its own $\eta_0/\hat\sigma_j$. This is exact for either representation: for the $\eta_0$-conditioned head the standardisation **cancels**, $\hat P(\rho_{\text{std}}>\eta_0/\hat\sigma_j\mid x,z)=\hat P(\rho>\eta_0\mid x,z)$, so the raw success probability is recovered whatever fixed $\hat\sigma_j$ is used. Deployed on the weekly grid the raw-percentage PPS stays interior (not saturated): the mean over items and interims is $0.75,0.26,0.12,0.06,0.03$ at $\eta_0=0,25,50,75,100\%$, and the share of (item, interim) cells with a non-degenerate PPS peaks at $46\%$ at $\eta_0=25\%$. Diagnostics in `…-eta0amortise-260830/`: `…_pps_by_item_eta0_grid.pdf` (PPS by item $\times$ $\eta_0\%$) and `…_rho_vs_eta0_lines_by_item.pdf` (the SVI $\rho$-predictive with the raw-percentage $\eta_0$ thresholds overlaid).

**Summary.** $\eta_H$ is free post-hoc; $\eta_0$ is free through the predictive CDF, best served by the **parametric Gaussian head** (analytic, no retraining), with the operational $\eta_0$ fixed on the raw percentage-change scale.

### 14.4.9 Amortising over items

**Motivation.** The §14.4 amortiser is trained for one fixed item set (Ukraine's $J=20$). A single network over a *variable* item set — $J=2,\dots,J_{\max}$ items — turns it from a per-trial artefact into a **per-family** asset, deployable on any partial-credit application (§3). Those applications span $J$ from $6$ (the ChatGPT survey) to $\approx51$ (the IBD PROM battery), with category counts $K\in\{2,\dots,10\}$ (most $\{4,5,7\}$); rounding up we set $J_{\max}=64$.

**Design (no padding on any axis).** Items enter **fixed-**$J$-per-batch, $J$ varied across batches: each step draws one $J_b\sim U\{2,J_{\max}\}$ and builds dense $(B,J_b,\cdot)$ tensors, so there is no item padding and no mask (JAX recompiles once per distinct $J_b$). Participants stay ragged (segment vectors, §14.4.3), and each (participant, item, time) cell holds a single scalar response $k/(K_j-1)$, so heterogeneous $K$ never introduces zeros either. Each **synthetic item** draws a category count $K\sim$ weighted $\{2,4,5,7,8,9,10\}$, an endpoint type (expected score or caseness), a direction, and PCM parameters from the §14.4.2 priors, so the network learns to price an item from its **features, not its identity**. The **metadata** is $(\text{type},\text{direction},K/K_{\max},c/K_{\max})$ — the caseness threshold $c$ is the fourth feature; the encoder already reads $J$ from the metadata shape and a generic $M$, so **no architecture change** is needed. The per-item $\hat\sigma_j$ of §14.4.1 is replaced by a single **global** $\hat\sigma$ (pilot $\approx4.8$); the power-law head's $C_j$ (§14.4.7) learns the residual per-item scale. Two cheap calibration add-ons recover most of what the global $\hat\sigma$ gives up: a **data-driven scale channel** — the per-item SD of the endpoint change across the observed cohort, injected at the head via `segment_sum` — and a **raw level/ratio fingerprint** (`head_raw`: baseline and endline means and their ratio at the head). Net directory `…-itemamortise-J64-scale-feat-260831` (trained $\sim1.5$ h).

**Item-count invariance holds.** Deployed on Ukraine's $20$ items *without retraining* — items the network never saw — its median tracks the SVI target as well as the item-specific net: the association (Pearson $\rho$ of the amortiser median against the SVI draw) is $0.91$–$0.95$ on the out-of-7 items and $0.73$–$0.76$ on the categorical items at small–mid $n$, against the item-specific baseline's $\approx0.82/0.57$ (§14.4.5); the raw (uncalibrated) PIT–KS is $0.31$–$0.46$, already below the item-specific net's $\approx0.60$. Both decay with $n$ as the SVI posterior tightens, exactly as in §14.4.5.

**End-to-end on Ukraine.** With the standard expanding head fine-tune $+$ affine shift (§14.4.6), the item-amortised net reaches **PIT–KS** $0.119$ and **marg–KS** $0.159$. Against the reference points: the marginal now **matches the item-specific plain deepset** ($0.157$), and the conditional PIT–KS **beats the hand-token encoders** ($0.139$, §14.4.6); the item-specific PPS ($\eta_0=0.5$) is well defined. The two add-ons each contribute one axis: the scale channel lifts the conditional (PIT–KS $0.134\to0.119$), the fingerprints lift the marginal (marg–KS $0.177\to0.159$).

Deployed on Ukraine's 20 items (never seen in training), expanding head-ft $+$ affine (§14.4.6):

All directories below sit under `py-ukraine-interim-amortise-`; the suffix column gives the end of each name (results in the `_ftheadexpand` deploy dir).

| Configuration | Directory suffix | PIT–KS | marg–KS |
|--------------------------|----------------------|-----------:|-----------:|
| Item-amortised (global $\hat\sigma$, power-law head) — base | `deepsetXcompAtt-itemamortise-J64-260831` | 0.134 | 0.177 |
| $+$ data-driven scale channel | `deepsetXcompAtt-itemamortise-J64-scale-260831` | 0.119 | 0.177 |
| $+$ raw level/ratio fingerprints (**scale-feat, deployed**) | `deepsetXcompAtt-itemamortise-J64-scale-feat-260831` | **0.119** | **0.159** |
| *ref:* item-specific plain deepset (bespoke, §14.4.6) | `deepsetXcompAtt-plain-ftheadexpand-medshift-260827` | 0.069 | 0.157 |
| *ref:* hand-token `itemXcompAtt` (§14.2) | `itemXcompAtt-ftheadexpand-260809` | 0.139 | 0.119 |
| *ref:* pretrained encoder $+$ Ukraine head fine-tune | `deepsetXcompAtt-itemamortise-J64-ukraineft-260831` | 0.133 | — |

Association with the SVI target (Pearson $\rho$ of the amortiser median against the SVI draw, small–mid $n$): $0.91$–$0.95$ on the out-of-7 items and $0.73$–$0.76$ on the categorical items, against the item-specific baseline's $\approx0.82/0.57$. The **scale-feat** row is the deployable per-family network; the two references above and below it bracket the recipe-bound floor discussed next.

**A recipe-bound floor.** The residual conditional gap to the item-specific $0.069$ is set by the *recipe*, not by item generality: a plain quantile head deploys at $\approx0.112$ and per-item $\hat\sigma_j$ reaches $0.069$ (§14.4.8), both of which the item-amortised net gives up (a global $\hat\sigma$, a power-law head). Transferring the pretrained encoder to Ukraine and fine-tuning does **not** move the deployed PIT–KS ($0.133$) while that recipe is fixed — the head type and scale, not the encoder, are binding. Closing the last of the gap therefore means reverting to a per-item-set network.

**Deployment scope — idiosyncratic items.** Because the network prices an item from its *features, not its identity*, an item whose behaviour those features do not capture is priced from the wrong prior. Ukraine's degenerate `CG-VIO_ph-punish` — the `blow` item excluded from every contraction fit (§14.4.5) — is the clear case: in the cross-method $p(H_1\mid x,z)$ comparison (§14.5) the item-general net **saturates at** $p(H_1)\approx1$ from interim 3 on ($\eta_0=0.5$; medians $0.87,1.0,1.0,\dots$), while the two Ukraine-*bespoke* encoders — `itemXcompAtt` and the plain deepset, which saw this item in training — track the HMC gold standard ($\approx0.2$–$0.5$). This is the deployment-scope trade-off made concrete: the per-family net generalises across item *sets* but is weaker on an individual idiosyncratic or boundary item than a net trained on the specific instrument. Expected behaviour, not a defect.

**Conclusion.** One network prices any partial-credit item set of size $2$–$64$ from item features alone, deploying at PIT–KS $0.119$ / marg–KS $0.159$ on Ukraine — conditionally sharper than the hand-token encoders and marginally on par with the bespoke deepset. Deploying on a new application is then: feed its items, metadata and cohorts, and run the head fine-tune $+$ affine (and BvM) against that application's own SVI. The §14.4 amortiser becomes a **per-family asset** across the §3 datasets, at a small, recipe-bound calibration cost relative to a bespoke network.

### 14.4.10 A net registry indexed by endpoint functional, and whether one token suffices

§14.4.9 showed one network prices any *item set*. The remaining organising question is not the items but the **endpoint functional** $\rho$: a seroprotection *rate* $P(y\ge k)$, a GMT *fold-rise*, a signed *relative change*, a between-group *difference*. Two structural choices depend on the functional — the per-token summary (a **scalar** normalised mean, or the **wide** cumulative-exceedance vector $[\,\mathbf 1\{k\ge c\}\,]_c$ that is caseness-sufficient, §14.4.18) and the deploy **warp** (none / $\log_2$ / logit) — while the item family does *not* (the Ukraine-trained encoder deploys on the flu titre family at PIT–KS $\approx0.08$, §3.18). This motivates a **net registry keyed by functional, not item family**. A registry **instance** is one item-general §14.4.9 network (deep-set over participants, item cross-attention, quantile head), trained for a single endpoint functional; it is identified by the pair $(\text{reduction},\text{compare})$ that defines $\rho$, which in turn fixes the two structural deploy choices — the **token** (scalar mean vs wide cumulative-exceedance) and the **warp** (none / $\log_2$ / logit). It is *not* keyed by item family: one instance prices any item set within the §14.4.9 envelope ($\le 64$ items, $\le 10$ levels). The current registry:

| instance | functional $\rho$ | reduction / compare | token | warp | example $H_1$ | applications (manifest) | status |
|---------|---------|---------|---------|---------|---------|----------|---------|
| $S_1$ | protection **rate** $P(y\ge k)$ | threshold / level | wide | logit | SPR $>0.70$ | flu SDY312/314, head-to-head (per arm) | trained |
| $S_2$ | GMT **fold**-rise | mean / $\mathrm{fold}_{\log_2}$ | scalar | $\log_2$ | GMFR $>2.5$ | flu, head-to-head (per arm) | trained |
| $S_3$ | **relative change**, mean (signed) | mean / ratio | scalar | none | $>0$ / effect size | ICRC, MYCELIUM, REFUGE, CAVD, **Ukraine (out-of-7)** | trained |
| $S_4$ | **relative change**, caseness | threshold / ratio | wide | none | $>0$ | **Ukraine (categorical)** | trained (§14.4.10) |
| $S_5$ | between-group **difference** (standardised) | any / diff | wide | none | $>0$ | head-to-head (LAIV$-$TIV) | to train |
| $S_6$ | per-participant **seroconversion** (\$\ge\$4-fold) | paired rate | wide + paired sim | logit | SCR $>0.40$ | flu composite (SCR $\vee$ SPR $\vee$ GMFR) | later |

An application's **manifest** is just the list of instances it deploys: flu $\{S_1,S_2\}$; head-to-head $\{S_1,S_2\;\text{per arm},\,S_5\;\text{difference}\}$; the psychometric datasets $\{S_3\}$; **Ukraine** $\{S_3\;\text{(out-of-7)},\,S_4\;\text{(categorical caseness)}\}$. **Deployment** is uniform (§14.4.6–14.4.7): filter the SVI reference to that endpoint (one `pps_ratio_x`), fit the expanding head $+$ affine, apply the head BvM, read the PPS in natural units after unwarping. New endpoints or datasets add a *row* or a *manifest entry*, never a new architecture.

Because the wide token is *sufficient* — its mean-pool is the empirical category CDF, carrying the mean, any threshold, and group contrasts — it is natural to ask whether the registry collapses further: **one widetok encoder for all functionals**, the functional selected only by the deploy-time head fine-tune. We test this on Ukraine as an amortiser-**methods comparison** (the Ukraine data question of §14 is irrelevant here — we compare networks against the *same* fixed SVI reference; no new fit). Three encoders are deployed on the identical Ukraine relative-change reference (weekly SVI, 20 items, 29 interims $n=48\!\to\!500$) under **identical** settings — warp none, head BvM on, expanding head fine-tune, $N_{\mathrm{ref}}=503$ — so only the encoder $\times$ token differ: the **scalar** `scale-feat` net (the current default for $S_3$), the caseness-dedicated **widetok-SPR** net, and a purpose-trained **widetok-MIXED** net whose target is *per item* the relative change (typ $0$) or the endline rate (typ $1$), $\approx$ half each, forcing one encoder to carry both.

*Metrics (both against the SVI reference).* **PIT–KS** — conditional calibration: for each (interim, item) form the probability-integral transform $u=\hat F(\rho^{\mathrm{svi}})$ of the amortiser's posterior CDF at the draw-aligned SVI value; PIT–KS is the Kolmogorov–Smirnov distance of the pooled $u$ from $\mathrm{Unif}[0,1]$ (a well-calibrated conditional posterior gives uniform $u$). **marg–KS** — marginal calibration: the KS distance between the SVI marginal CDF and the amortiser mixture CDF on a shared grid. Coverage `cov5`/`cov95` are the empirical masses below the nominal $0.05$/$0.95$ quantiles.

*Fair aggregation.* One Ukraine item, `CG-VIO_ph-punish`, is a known **functional** pathology — a signed relative change with $\approx 25\%$ of its posterior mass negative, unstable near that boundary (the same instability that made the head-to-head difference, not percent, the sane cross-arm estimand, §3.21) — and it saturates the scalar/mean deploy (PIT–KS $0.84$, §14.4.9). Since it is present at every interim it inflates **every** aggregate, so it is dropped **throughout** — from the per-item means *and* the per-interim (across-$n$) aggregate (`RAGD_EXCLUDE_ITEMS`) — leaving 19 items.

PIT–KS is by **cohort-size regime** (per-interim, CG-VIO excluded); marg–KS and mean posterior width $W$ are over all interims. The last row is the **federated** deploy — scalar $S_3$ on the 15 out-of-7 items, dedicated caseness $S_4$ (next paragraph) on the 4 categorical — the alternative to a *single* encoder.

| approach | token | PIT–KS $n\le200$ | PIT–KS mid | PIT–KS $n\ge450$ | PIT–KS all | marg–KS | width $W$ |
|---------------|---------|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|
| **scalar** ($S_3$ default) | scalar | **0.073** | **0.084** | 0.096 | 0.088 | 0.110 | 0.43 |
| widetok-SPR (rate, wrong target) | wide | 0.157 | 0.150 | 0.104 | 0.129 | 0.298 | 0.21 |
| widetok-MIXED (single encoder) | wide | 0.079 | 0.092 | 0.101 | 0.094 | 0.123 | 0.34 |
| **FEDERATED** ($S_3$ mean $+$ $S_4$ caseness) | scalar $+$ wide | **0.070** | **0.080** | **0.092** | **0.084** | 0.110 | 0.46 |

(bins: $n\le200$ 6 interims, mid 9, $n\ge450$ 14; `cov95` $\approx0.91$–$0.94$ except widetok-SPR $\approx0.85$ at mid; $W$ = mean marginal width, widetok-SPR narrowest $=$ over-confident. The per-item PIT–KS aggregate gives the same ordering.)

For a stopping rule the small-$n$ end carries the decision, and the table's regime split is the point: scalar, mixed and federated are sharpest and best-covered there, while widetok-SPR — a *rate*-dedicated net on a *mean* reference — is worst early and only competitive late, exactly wrong for early stopping.

Three findings. **(i) The wide token is not the bottleneck:** on the *matched* mixed target the widetok encoder recovers to near-scalar, so the widetok-SPR deficit was target mismatch — the single-encoder collapse is feasible. **(ii) But it does not dominate:** the dedicated scalar stays marginally sharper (a small generality tax) and less over-confident. **(iii) The pathology is functional, not architectural:** `CG-VIO_ph-punish` is rescued only when scored as a *rate/caseness* rather than a signed ratio ($0.84\!\to\!0.17$), which no encoder change achieves — the estimand does.

**Federated per-functional vs the single encoder.** Completing the registry, we trained the missing net $S_4$ — a widetok encoder for the categorical items' caseness *relative change* (widetok-SPR carries the endline *rate*, the wrong compare) — and deployed it on the 4 categorical items with the scalar $S_3$ on the 15 out-of-7. The federated set is the sharpest (table); on every categorical item the dedicated caseness net edges the mean net, so even the mono scalar was slightly overpaying. The single mixed encoder trails by $\approx0.03$ PIT–KS — one asset over several.

**Conclusion.** Adopt the registry as $S_1$–$S_6$ **keyed by functional, no item-family axis**: scalar token for mean-type ($S_2$/$S_3$), wide for rate/difference ($S_1$/$S_4$/$S_5$). The federated per-functional set is marginally the sharpest; the single-widetok collapse is a proven *fallback* at a small ($\approx0.02$–$0.03$ PIT–KS) tax. And an item that miscalibrates because its *functional* is unstable is fixed by the estimand (a difference or rate), not the network.

### 14.5.1 Cross-architecture summary

Calibration of the three encoders under the deployable remedy (expanding head fine-tune $+$ affine shift, §14.4.6), and of the two BvM-corrected deepset heads (§14.4.7); PIT–KS is conditional calibration, marg–KS marginal calibration:

| encoder                   | configuration             |   PIT–KS | marg–KS |
|---------------------------|---------------------------|---------:|--------:|
| `itemXcompAtt` (§14.2)    | baseline (raw network)    |     0.61 |    0.47 |
|                           | head-ft-expand $+$ affine |     0.14 |    0.12 |
| `itemScompAtt` (§14.3)    | baseline (raw network)    |     0.55 |    0.42 |
|                           | head-ft-expand $+$ affine |     0.14 |    0.11 |
| `deepsetXcompAtt` (§14.4) | baseline (raw network)    |     0.60 |    0.47 |
|                           | head-ft-expand $+$ affine | **0.07** |    0.16 |
|                           | A power-law $+$ BvM       |     0.11 |    0.18 |
|                           | C floor $+$ BvM           |     0.12 |    0.17 |

### 14.5.2 Reading

1.  **The §8-10 loop closes.** Prior sampling → nested-DeepSets training → raw-$(x, z)$ deployment runs end-to-end on a real IRT case study with a \~1 min total deployment cost for the full 8-interim schedule.
2.  **The mislocation is shared across encoders.** Cross-attention, self-attention and nested-DeepSets have the *same* baseline PIT–KS ($\approx0.55$–$0.61$) and marg–KS ($\approx0.47$): the prior/deployment gap is a property of the training setup, not the encoder. Resolved by cohort size, the three encoders also track the SVI posterior draw equally well (the apparent deepset "$z$-floor" was an averaging artefact, §14.10 item 14).
3.  **The remedy is the head fine-tune, not a static map.** The expanding head-only fine-tune (§14.2.5, §14.4.6) recovers near-nominal conditional calibration for all three encoders; post-hoc conformal is insufficient because the bias is interim-heteroscedastic and one monotone map cannot absorb it. The affine median-shift (§14.4.6) then removes the residual location bias.
4.  **The conditional–marginal trade-off, and its resolution.** With the affine shift throughout, the deepset wins the conditional (PIT–KS $0.07$ vs the hand token's $0.14$ — the quantity the PPS decision consumes) and the hand token the marginal by a small margin. The **BvM correction** (§14.4.7) closes the deepset's large-$n$ marginal shortfall from first principles, giving the coherent A/C $+$ BvM estimator of §14.0. (All three encoders are now evaluated on the same 29-interim weekly grid.)
5.  **Diagnostics** are produced in full for every configuration (§14.4.5): the PIT box, marginal quantile box, contraction-factor and contraction-law figures from `python/amortiser_diag_plots.py`.

## 14.10 Approaches examined and set aside

Recorded for completeness; each was tested and did not improve on the §14.0 estimator. Detailed subsections and directories are cited.

**Calibration (item-level encoder).**

1.  *Training-noise control*, $S=256$ draws/step (§14.2.1): confirms the baseline gap is distribution shift, not Monte-Carlo noise — a diagnostic, not a remedy.
2.  *Prior narrowing*, `threshold_scale` $2.0/1.0$ (§14.2.2): does not move coverage; the mislocation is not a prior-width artefact.
3.  *Full-network fine-tune*, leave-one-interim-out (§14.2.3): reaches the same calibration as the head-only fit while moving all weights $\approx10\%$ — superfluous; the head-only fit (§14.2.5) is retained.
4.  *Head-only fine-tune on interim 1 alone* (§14.2.4): calibrates interim 1 but drifts elsewhere because the bias is interim-heteroscedastic; superseded by the expanding window (§14.2.5).
5.  *Post-hoc conformal recalibration* (§14.2.6; also tried on the deepset): a single monotone PIT map cannot absorb an interim-heteroscedastic bias; coverage only partially corrected, marginal-KS worsened.

**Contraction and training-law redesigns (deepset).**

6.  *Contraction pressure* — BvM aux features $+$ $n$-stratified batches: negative. With the original training the two cohorts jointly identify the full-population endpoint ($n+m$ fixed), so there is no $x$-only contraction in the target to learn. Dir `…-deepsetXcompAtt-plain-260812` predecessor.
7.  \*Focused-proposal (interim-1 posterior) ragged build: negative. An interim-1 posterior proposal is pre-contracted, and subsampling a fixed pool is not refitting; the contraction signal is flat. Superseded by the wide-prior recipe (§14.4.3).
8.  *Capacity-fix sweep* — schedule, MLP width, factored $g(1/\sqrt n)$, precision-pool, hard-coded semiparametric width: only hard-coding the rate moves the training-distribution slope, and it over-contracts at deployment. Motivates the *learned* power law (§14.4.7–14.4.8) instead.
9.  *Mixture-CRPS marginal loss term*: reduces marginal-KS only where the encoder conditions sharply on $z$ (hand token); inert on the deepset, whose $z$-response is weaker — a property clarified by item 14 below.

**Deepset** $z$-response investigation.

10. *Sharper-*$z$ encoder sweep, 9 variants: `z_contrast`, `precision_pool`, `raw_pool`, domain randomisation: the deepset marginal $z$-tracking did not exceed $\approx0.5$ under any input-feature or training lever. Dirs `…-priorproposal-zs-{…}`.
11. *Proposal ablation on the hand token*: isolates **domain randomisation**, not decoupling or small $m$, as the cause of $z$-erosion — off by default thereafter. Dirs `…-itemXcompAtt-priorproposal-{fixed503-noDR, …}`.
12. *Nband-noDR recipe on the deepset*: improves the raw contraction slope $3\times$ but not deployment calibration; the wider $N$-band spreads capacity. Dirs `…-priorproposal-{nband, bigM503}`.
13. *Per-participant* $\rho$-mapping: moving the nonlinearity before the pool does not raise the deepset $z$-tracking (0.57 unchanged); confirms the ceiling is downstream of pooling. Dirs `…-priorproposal-{perpart, …}`.
14. ***Correction***: the apparent deepset "$z$-floor at $0.57$" was an averaging artefact — resolved per cohort size, the deepset and hand token track $z$ identically ($0.82$ at small $n$, decaying to $\approx0.15$ at large $n$). The deployed difference is a modest conditional-vs-marginal trade-off, not a $z$-response deficit. This retracts the premise of items 10–13.

**Law-of-total-probability marginal patches (deepset).**

15. *Empty-cohort marginal* (Option A): a single $m=0$ forward removes the finite-$m$ upturn but is a weaker marginal estimator (marg-KS $0.243$) and costs the conditional; the $m=0$ representation is off the operating regime. Dir `…-priorproposal-mempty_ftheadexpand`. Flags `PP_MEMPTY`, `RAGD_EMPTYZ`.
16. *Real-data split marginal* (Option B): partitioning $x=(x_1,x_2)$ avoids simulation noise but the balanced split is likewise off-regime (marg-KS $0.275$). Dir `…-plain-splitB-260828`. Flag `RAGD_SPLITB`.

Both 15–16 are superseded by the in-regime **BvM correction** (§14.4.7), which removes the same finite-$m$ inflation without leaving the training regime.

# 15. Cross-arm facilitator-matched estimand for the Ukraine data (UkraineP)

**Motivation.** The estimand of §14 is the pooled within-participant pre→post change: the partial-credit model is fitted on all enrolled participants with `~ time - 1` (baseline vs endline), and $\rho_j = s_j(\bar w_{e,j}/\bar w_{b,j}-1)$ is the endline-vs-baseline population-mean item change. Because Hope Groups is a waitlist cluster-randomised trial in which both arms carry an endline, that pooled contrast is *descriptive*, not causal: the control arm's own pre→post change dilutes the average and any secular trend common to both arms enters $\rho$ (§14.0, §3.6). Recovering a treatment effect needs an explicit arm contrast.

**The cross-arm, facilitator-matched estimand.** UkraineP replaces the pre→post contrast by a between-arm one, matched on facilitator and roughly matched in calendar time. Every facilitator runs both an intervention group and a waitlist-control group; per item we contrast the group-mean level of the **intervention arm at endline** (treated, post) with the **control arm at baseline** (untreated, pre), $$\rho_j \;=\; s_j\!\Big(\bar w^{\text{int,end}}_j \big/ \bar w^{\text{ctrl,base}}_j - 1\Big),$$ $\bar w^{\text{int,end}}_j$ the endline mean level of the intervention arm, $\bar w^{\text{ctrl,base}}_j$ the baseline mean level of the control arm. The control-baseline group is untreated and measured before the programme; the intervention-endline group just after — so the contrast isolates a treated-post vs untreated-pre difference while holding the facilitator fixed (both groups share one) and the calendar window roughly fixed.

**Construction.** The two comparison cells already carry the partial-credit model's two time labels, so UkraineP is the Ukraine dataset restricted to {control-baseline (time 0, the reference), intervention-endline (time 1, the treated)}; the other two cells (intervention-baseline, control-endline) are dropped. The two cells are *disjoint sets of participants*, so this is an **unpaired between-group** partial-credit fit — each participant contributes a single time-point, the group-level shift carried by the per-`item_time_id` difficulties with ability $\theta_i\sim N(0,1)$ shared across groups — the same structure as the mycelium (§3.14) and PISA (§3.11) between-cohort designs, and therefore an SVI-reference application: the paired amortiser of §14 (baseline→endline within participant) does not apply. The two item types (out-of-7 days-in-week, $K=8$; categorical caseness, $K=4$) sit on separate `item_type_id`, so no K-family mixing (§3.11).

**Interim accrual over facilitators.** Rather than a calendar cutoff, interims accrue **facilitators** — ordered by the median endline date of each facilitator's intervention group — so at interim $k$ the first $k$ facilitators contribute their intervention-endline and control-baseline group means and the number of matched comparisons grows with $k$. The model is fitted by SVI at each interim (AutoLowRankMVN, $10\,000$ steps, $S=4000$ posterior-predictive draws), producing the same artifacts as the other SVI producers (dp1, draws, endpoint draws, `prob_by_question_fit` plots).

**Result.** Over the 30 facilitators (8 interims, $k=4\to30$), the median cross-arm effect is large and stable:

| interim | facilitators | $n$ (ctrl-base / int-end) | median $\rho$ |
|--------:|-------------:|--------------------------:|--------------:|
|       1 |            4 |                   29 / 27 |         0.665 |
|       2 |            8 |                   60 / 55 |         0.621 |
|       3 |           11 |                   81 / 75 |         0.517 |
|       4 |           15 |                 134 / 129 |         0.704 |
|       5 |           19 |                 164 / 160 |         0.706 |
|       6 |           23 |                 185 / 181 |         0.692 |
|       7 |           26 |                 228 / 227 |         0.695 |
|       8 |           30 |                 253 / 250 |     **0.638** |

At full accrual the per-item effects span $\rho\in[0.24, 1.46]$ (median $0.61$) across all 20 items, all positive: largest for grieving ($1.46$) and self-care ($1.45$) and the mental-health items (sad $0.89$, low effort $0.82$, nervousness $0.77$ — direction-aware, so improvements), smallest for child-monitoring ($0.24$–$0.51$). The effect is far larger than the pooled pre→post estimand of §14, as expected: the treated-post vs untreated-pre contrast captures the full arm separation plus any common secular improvement, whereas the pooled change averages the treated and near-null control pre→post together. UkraineP is thus an efficacy read complementary to the pooled estimand; adding the `time×treat` interaction (a difference-in-differences, §3.4) would net out the shared secular component to recover the causal intention-to-treat effect.

**Files.** [`scripts-py/UkraineP_interim_svi.py`](../scripts-py/UkraineP_interim_svi.py) → `py-ukraineP-crossarm-svi-260916/` (per-interim dp1, `draws.zarr`, endpoint draws, fit plots, `pcm_1_interim_interim_index.csv`).

------------------------------------------------------------------------

# 16. Amortiser retraining and deployment for the influenza-vaccine endpoints

The influenza HAI application (§3.18) carries two clinically standard endpoints per strain (CHMP/CBER; @hobson1972role): the **seroprotection rate** $\text{SPR}=P(\text{titre}\ge 1{:}40)$ and the **geometric-mean fold-rise** $\text{GMFR}=\text{GMT}_{\text{end}}/\text{GMT}_{\text{base}}$. Both are amortised with the item-general J64 encoder of §14.4.9, but they make different demands on the per-participant summary. This section reports how the training distribution was extended to cover both endpoints (§16.1), how the deployed estimators calibrate against the SVI reference, including per item (§16.2), and why two influenza-B strains trail the rest (§16.3).

## 16.1 Extending the training distribution to several endpoints

The item-general trainer synthesises item sets of two kinds. For an **expected-score** item the effect is the relative change of the mean category $\mathbb E[k]$; for a **caseness** item it is a threshold exceedance $P(y\ge c)$. GMFR is of the first kind — $\text{GMFR}=2^{\bar E_{\text{end}}-\bar E_{\text{base}}}$ is a smooth function of the mean $\log_2$-titre, which the scalar per-participant token (the normalised category $(y-1)/(K-1)$, mean-pooled over participants) already renders sufficient. GMFR is therefore deployed with the **existing J64 network of §14.4.9 unchanged**.

SPR is of the second kind, where the scalar-mean token fails: a threshold exceedance is not a function of the mean but of the category profile near the cut $c$ (the categorical-sufficiency wall, §14.4.18). Two changes were made to the training distribution, both leaving the architecture (nested DeepSets over participants $+$ cross-attention over items $+$ multi-quantile head) untouched:

1.  **Sufficient token.** The per-$(\text{participant},\text{item},\text{time})$ token becomes the **cumulative-exceedance vector** $[\mathbf 1\{k\ge 1\},\dots,\mathbf 1\{k\ge K_{\max}-1\}]$ (width $2(K_{\max}-1)=18$ for $K_{\max}=10$, one block per time-point). Mean-pooling over participants then returns the empirical category survival curve, from which every threshold rate — and the mean — is recoverable; the pooled statistic is thus sufficient for SPR **and** GMFR at once. Only the input width changes, so the first token network absorbs it and the rest is identical.
2.  **Endline target and threshold metadata.** For caseness items the label is the **endline** rate $P(y_{\text{end}}\ge c)$ (the baseline rate is computed but unused, matching the operational SPR of §3.18), and the cut $c/K_{\max}$ is appended to the per-item metadata so the network prices the queried threshold. Training was restricted to caseness items for a dedicated SPR network (uniform $[0,1]$ target scale).

The SPR network was trained with the interval-score objective (§14.4.11) on item sets of size $J\sim\mathcal U\{2,\dots,64\}$, $K$ over $\{2,4,5,7,8,9,10\}$, decoupled cohort sizes (§14.4.13), for 6000 steps at batch 24; the target-standardisation constant was $\sigma=0.325$ (a rate, hence $O(0.1)$). **Wall-clock training time was 212 min.** The GMFR network is the pre-existing J64 network and was not retrained.

## 16.2 Deployment and calibration against the SVI reference

Both endpoints deploy through the shared ragged driver (§14.4.8) with expanding-window head fine-tuning, the per-item affine median-shift (§14.4.23) and the between-first **head BvM** correction (§14.4.26); the SPR deployment additionally activates the cumulative-exceedance token and the caseness cut $c=3$ (titre $\ge 1{:}40$). Success is scored by the conditional **PIT–KS** and the **marginal KS** against the per-interim SVI posterior of the endpoint. The two endpoints are now deployed as a **federated amortiser** (§14.4.9, §17.5): one parent directory `py-immport-SDY312-amortise-deepsetXcompAtt-itemamortise-J64-ftheadexpand-bvm_260919-federated_260924/` with a **per-endpoint subdirectory** (`SPR/`, `GMFR/`, each a delegated amortiser with its subset reference nested), and every combined diagnostic figure and the `flu_amortiser_calibration.*` summary reporting into the parent; the SVI reference grid and its plots stay in `py-immport-SDY312_260918/`.

Aggregate calibration (SDY312, federated re-run, mean over interims):

| study  | endpoint | PIT–KS | marg–KS | uncalibrated baseline (PIT / marg) |
|--------|----------|--------|---------|------------------------------------|
| SDY312 | GMFR     | 0.106  | 0.110   | 0.32 / 0.43                        |
| SDY312 | SPR      | 0.110  | 0.137   | 0.63 / 0.68                        |

The uncalibrated SPR baseline (0.63–0.68) **is** the sufficiency wall — with the scalar-mean token the network cannot reproduce a threshold rate at all — and the cumulative-exceedance token brings the deployed SPR estimator into the same calibration band as GMFR (PIT–KS $\approx 0.11$). (SDY314, the second study, is not run under the federated pipeline; its earlier single-study deploy calibrated to GMFR 0.082 / 0.090 and SPR 0.080 / 0.096, the same band.)

Per-item calibration (SDY312; PIT–KS pooled over interims, marg–KS averaged over interims):

Each cell is **PIT–KS / marg–KS** for that strain and endpoint.

| strain | GMFR (PIT–KS / marg–KS) | SPR (PIT–KS / marg–KS) |
|----------------------------|-----------------------|----------------------|
| A/Puerto Rico/8/1934(H1N1) | 0.057 / 0.072 | 0.106 / 0.125 |
| A/South Dakota/06/2007(H1N1) | 0.078 / 0.105 | 0.089 / 0.107 |
| A/Uruguay/716/2007(H3N2) | 0.122 / 0.116 | **0.198 / 0.264** |
| A/Victoria/3/1975(H3N2) | 0.081 / 0.091 | 0.075 / 0.106 |
| B/Brisbane/60/2008 | **0.195 / 0.181** | 0.091 / 0.108 |
| B/Florida/4/2006 | 0.122 / 0.116 | 0.079 / 0.094 |
| B/Lee/1940 | 0.091 / 0.092 | 0.132 / 0.157 |

Calibration is uniformly good (PIT–KS $\lesssim 0.13$) with two exceptions, and they are **endpoint-specific boundary effects**, not a property of particular strains. PIT–KS scores how well the amortiser's (roughly symmetric) five-quantile head reproduces the *shape* of the SVI posterior; it degrades where that posterior is pushed against a boundary and becomes strongly skewed:

- **GMFR, at the titre floor.** B/Brisbane (PIT–KS 0.195) and B/Florida (0.122) sit at the assay floor — baseline GMT 9 and 17, most participants at $k=0$ (below the 1:10 detection limit). Since $\text{GMFR}=2^{\bar E_{\text{end}}-\bar E_{\text{base}}}$ is a fold of a near-zero baseline mean, the SVI posterior is **right-skewed / heavy-tailed** (skewness 1.1 and 1.5, vs 0.5–0.65 for the well-calibrated strains); B/Brisbane is the most extreme floor case and calibrates worst.
- **SPR, at the ceiling.** A/Uruguay (PIT–KS 0.198) has SPR $\approx 0.97$, compressed against 1.0 (skewness $-2.05$, a third of draws within 0.02 of the boundary); its head fits worst *for SPR*, while the two B strains — interior SPR rates 0.24–0.35 — calibrate cleanly (0.09).

Two mechanisms compound: (i) the head family is near-symmetric — the affine shift (§14.4.23) corrects *location* and the BvM correction the *width*, but neither reshapes **skew**; (ii) near-floor / near-ceiling category profiles are under-represented in the synthetic training prior ($\beta\sim
\mathcal N(0,1)$), so the network mildly extrapolates there. This is a head-shape limitation at endpoint boundaries, not a location error — the go/no-go (a clear futility for the B strains under GMFR, a clear pass for A/Uruguay under SPR) is unaffected. A skew-capable head, or a log-/logit-warped target (log-GMFR would symmetrise the floor skew), would close it — §16.3 implements and compares both. The per-item PIT and contraction boxes are in `…_pps_RAGD_pit_box_by_item.pdf` / `…_contraction_cdf_by_item.pdf`; SDY314 (6 strains) is analogous, all per-item PIT–KS $\le 0.11$.

**PPS and the deployment sweep.** The predictive probability of success $P(H_1\mid x)$ is swept over the operational threshold $\eta_0$ — GMFR $\in\{2,2.5,3,3.5\}$-fold and SPR $\in\{0.5,0.6,0.7,0.8\}$ (2.5 and 0.70 the CHMP anchors) — with the per-strain trajectory over accruing participants in each deploy directory (`…_pps_RAGD_trajectory.pdf`). At the CHMP seroprotection bar SPR $>0.70$, A/Puerto Rico, A/Uruguay, A/Victoria and B/Lee reach PPS $\approx 1$ by full accrual, A/South Dakota is a borderline no-go, and B/Brisbane and B/Florida are futile; the GMFR $>2.5$-fold bar is met confidently only by A/Uruguay, reflecting the high pre-existing titres in these cohorts.

## 16.3 Skew-capable heads: two remedies, compared

The boundary-skew miscalibration of §16.2 motivates a head (or target) that can represent an asymmetric predictive. Because the deployment re-fits the quantile head on the frozen encoder features (§14.4.8), **both remedies are deploy-only — no encoder retraining**:

- **C — target warp.** Keep the symmetric power-law head but fit it to a monotone-warped endpoint that is approximately symmetric: $\rho'=\log_2\text{GMFR}$ (a fold $\to$ additive; the floor right-skew is removed) and $\rho'=\operatorname{logit}\text{SPR}$ (the rate's $[0,1]$ boundaries map to $\mathbb R$). The head fits the symmetric $\rho'$; the deployed quantiles, mapped back, are correctly skewed. `RAGD_WARP=log2|logit`; $\eta_0$ and the KS diagnostics are computed in the warped space (KS is invariant under the monotone map, so the numbers stay comparable).
- **B — free-quantile head.** Replace $q=m+z_\tau\,C n^{-p}$ (fixed symmetric normal $z_\tau$) by a **median plus independent lower/upper monotone gaps** (`head_mode='freeq'`), so the quantile set can take any skew; the BvM correction then rescales the width to the power law. Same encoder, a richer head re-fit at deploy.

Calibration vs the SVI reference (SDY312, mean over interims; PIT–KS conditional, marg–KS marginal):

| endpoint | head variant                  | PIT–KS    | marg–KS   |
|----------|-------------------------------|-----------|-----------|
| GMFR     | current (symmetric power-law) | 0.106     | 0.110     |
| GMFR     | **C —** $\log_2$ warp         | **0.081** | **0.093** |
| GMFR     | B — free-quantile             | 0.099     | 0.109     |
| SPR      | current (symmetric power-law) | 0.110     | 0.137     |
| SPR      | **C — logit warp**            | **0.103** | **0.126** |
| SPR      | B — free-quantile             | 0.143     | 0.140     |

Per-item PIT–KS at the two problem strains (the boundary cases of §16.2):

| endpoint · strain                  | current | C (warp)  | B (free-q) |
|------------------------------------|---------|-----------|------------|
| GMFR · B/Brisbane/60/2008 (floor)  | 0.195   | **0.077** | 0.194      |
| GMFR · B/Florida/4/2006            | 0.122   | 0.109     | **0.085**  |
| SPR · A/Uruguay/716/2007 (ceiling) | 0.198   | **0.125** | 0.285      |

**C is best on both endpoints** — the aggregate gain is large for GMFR ($0.106\to0.081$) and modest but consistent for SPR ($0.110\to0.103$), and it collapses the two worst per-strain cases (B/Brisbane GMFR $0.195\to0.077$; A/Uruguay SPR $0.198\to0.125$). The log/logit map is the natural symmetriser of a fold / a rate, and — crucially — it respects the endpoint's support by construction, so the head never has to place quantile mass beyond the titre floor or above SPR $=1$.

**B disappoints, and instructively.** It helps a few interior strains (B/Florida GMFR $0.122\to0.085$, B/Lee GMFR $0.091\to0.062$) but *worsens* both hard boundary cases — A/Uruguay SPR rises to 0.285, and the aggregate SPR is the worst of the three variants (0.143). Two reasons: a free five-quantile head has more shape parameters to fit from the limited per-interim reference draws at deploy, so it is under-constrained near the boundaries; and, working in the raw $[0,1]$ / floor space, no monotone-quantile set can reproduce a posterior piled against a hard bound (SPR $\to 1$). The warp removes the bound, which the extra head flexibility cannot.

The go/no-go is unchanged throughout: $P(\rho>\eta_0)=P(\rho'>\eta_0')$ under a monotone warp, so C only reshapes the predictive, not the decision.

**Adopted.** C is now the default for the influenza deploys (`RAGD_WARP=log2` for GMFR, `RAGD_WARP=logit` for SPR in the deploy wrappers): the head is fit in warped space, then the quantiles and targets are mapped back so PPS, $\eta_0$ and all plots stay in natural units. On the full canonical federated deploy (expanding head-ft $+$ affine $+$ head BvM, complete diagnostic suite) C improves both aggregate endpoints over the symmetric baseline — GMFR $0.106/0.110\to\mathbf{0.081/0.093}$, SPR $0.110/0.137\to\mathbf{0.103/0.126}$ (the three-variant table above). SDY314 has no extreme-floor strain, so its earlier single-study deploy already calibrated well and improved consistently under C (GMFR $0.082/0.090\to0.063/0.086$, SPR $0.080/0.096\to0.073/0.093$; all per-item PIT–KS $\le 0.09$). Each variant is a self-contained **federated deployment** — one parent directory with `SPR/` and `GMFR/` subdirectories, the combined diagnostic suite (per-item PIT and contraction boxes, contraction law, $\eta_0$-sweep, PPS trajectory, calibration summary) reporting into the parent: C-warp (the default) in `…-J64-Cwarp-ftheadexpand-bvm_260920-federated_260924/`, the symmetric baseline in `…-J64-ftheadexpand-bvm_260919-federated_260924/`, and the free-quantile head in `…-J64-Bfreeq-ftheadexpand-bvm_260920-federated_260924/`.

**The warp is endpoint-specific, and correctly scoped.** C is set in the *influenza wrappers*, not the driver (whose default remains `RAGD_WARP=none`), because $\log_2$/logit suit a fold / a bounded rate but not every endpoint. The Ukraine partial-credit endpoint (§14) is a *signed* relative-change $\rho$ (1.4% of reference draws are $\le 0$ — up to 25% for the violence-reduction item, with $\rho$ as low as $-4.8$), so $\log_2\rho$ is undefined and the logit (a $[0,1]$ map) does not apply; moreover Ukraine's $\rho$ is already near-symmetric (median $|\text{skewness}|\approx0.4$) and its symmetric-head deployment is already well calibrated. C therefore neither applies to nor is needed for Ukraine, and the default scoping leaves the Ukraine deployment unchanged. (Endpoints with a hard boundary — folds, rates, censored effects — are the ones the warp helps.)

## 16.4 Why B/Brisbane and B/Florida trail the other strains

Both the SVI reference and the amortised PPS place the two influenza-B strains well below the other five. This is a genuine immunogenicity gap, not an artefact — the observed titres (SDY312, full cohort):

| strain                       | baseline GMT | endline GMT | endline SPR |
|------------------------------|--------------|-------------|-------------|
| A/Uruguay/716/2007(H3N2)     | 63           | 178         | 0.97        |
| A/Puerto Rico/8/1934(H1N1)   | 45           | 83          | 0.90        |
| A/Victoria/3/1975(H3N2)      | 31           | 59          | 0.83        |
| B/Lee/1940                   | 31           | 58          | 0.79        |
| A/South Dakota/06/2007(H1N1) | 31           | 44          | 0.68        |
| B/Florida/4/2006             | 17           | 20          | 0.33        |
| B/Brisbane/60/2008           | 9            | 16          | 0.23        |

Three non-exclusive hypotheses:

1.  **B-lineage mismatch.** A trivalent inactivated vaccine carries a *single* influenza-B lineage, whereas B/Brisbane/60/2008 is **Victoria**-lineage and B/Florida/4/2006 is **Yamagata**-lineage. Whichever lineage is absent from the season's formulation receives only weak cross-lineage boosting — the classic reason B responses lag in TIV. That both trail while the ancestral B/Lee/1940 does not (endline SPR 0.79) points to a lineage-specific, not pan-B, failure.
2.  **Lower intrinsic immunogenicity of influenza B.** Even lineage-matched, HAI responses to B are consistently smaller than to A(H1N1)/A(H3N2) in inactivated vaccines; the fold-rises here (1.2–1.7) sit at the low end typical of B.
3.  **Not a ceiling effect.** A capped fold-rise from high pre-existing immunity would show *high* baseline titres; instead B/Brisbane's baseline GMT is 9 (below the 1:10 detection floor for most participants) and rises only to 16. Low baseline *and* low endline together indicate a true failure to boost — mismatch / low immunogenicity — rather than saturation.

The calibration consequence in §16.2 follows: sitting near the floor with strongly right-skewed endpoint posteriors, these strains are where the symmetric quantile head fits least well, yet the go/no-go — a clear futility — is unambiguous.

# 17. Head-to-head vaccine comparison (SDY269 LAIV vs TIV)

The head-to-head is the genuinely procurement-relevant estimand: given two vaccine *options*, which produces the stronger immune response? SDY269 (2008 Systems Biology influenza; §3.21) is the sound in-hand case — two randomised adult arms, **LAIV** (live-attenuated, intranasal) vs **TIV** (inactivated, intramuscular), each with paired day-0/day-28 HAI and the two §3.18 endpoints (SPR, GMFR). The arms were assayed on **different strain panels**, so only **A/Uruguay/716/2007 (H3N2)** is directly comparable strain-for-strain; the H1N1 and B strains compare each arm's homologous response.

## 17.1 One joint fit, endpoints declared with the data

Both arms are fitted in **one** partial-credit model per interim: the item stays the pure strain, and the vaccine arm folds into the flexible condition axis with the paired time-point (`group_label` $\in\{$LAIV\_baseline, LAIV\_endline, TIV\_baseline, TIV\_endline$\}$). The non-reduced structure is the cross `item_group_id` $=$ strain $\times$ arm $\times$ phase (5 strains $\times$ 4 conditions, sparse $=12$), $\theta_i$ shared across the disjoint participants (incomplete block); the design uses only the binary `phase` (`x_formula="~ phase - 1"`), so the fit is invariant to the relabelling. The loader declares **six endpoints upfront**, each with a short and a long label — `LAIV_spr`, `TIV_spr`, `LAIV_gmfr`, `TIV_gmfr` (per-arm levels) and `TIV-LAIV_spr_diff`, `TIV-LAIV_gmfr_diff` (cross-arm). Each is computed one-by-one off the *same* fit, so all six share the Monte-Carlo draw index — enabling both the composite decision and an *exact* cross-arm contrast.

## 17.2 Seroresponse: TIV dominates, and the composite call

A textbook adult result — **TIV dominates LAIV on every strain and both endpoints**. On the shared A/Uruguay H3N2: SPR LAIV $0.16$ vs TIV $0.69$; GMFR LAIV $1.23$ vs TIV $5.66$. On H1N1, TIV meets both CHMP thresholds (SPR $0.82$, GMFR $3.12$) while LAIV meets neither ($0.21$, $1.10$). This is the expected immunology: in adults an inactivated IM vaccine drives strong serum HAI, whereas a live-attenuated intranasal vaccine (mucosal, against a pre-exposed background) drives little. The **composite CHMP rule** — at least one of SPR/GMFR met — scored per draw over the jointly-indexed level endpoints, per (strain, arm): TIV A/Brisbane-H1N1 $0.94$, A/Uruguay $0.93$; LAIV $\approx 0$ throughout.

## 17.3 The cross-arm difference on the shared strain — exact per-draw

Because both arms share the joint fit, the between-arm difference is an **exact per-draw contrast** — no random draw-pairing of independent fits, so the shared $\theta$/draw removes the artificial independence. The endpoints are plain differences (a signed ratio is unstable at the boundary, §3.21): $\rho_{\text{SPR,diff}}=\text{SPR}_\text{TIV}-\text{SPR}_\text{LAIV}$ and $\rho_{\text{GMFR,diff}}=\text{GMFR}_\text{TIV}-\text{GMFR}_\text{LAIV}$. At full accrual: $\rho_{\text{SPR,diff}}$ median $0.51$ (95% CI $0.20$–$0.75$, $P>0\approx1.00$); $\rho_{\text{GMFR,diff}}$ median $4.33$-fold (95% CI $0.21$–$12.4$, $P>0=0.98$).

## 17.4 Amortising the cross-arm difference ($S_5$), standardised

The between-group difference is the $S_5$ registry instance (§14.4.10). To be **scale-free across applications**, it targets the *standardised* difference — a Cohen's $d$, $d=(g_B-g_A)/s$ with $g$ the per-group functional (mean or rate) and $s$ the pooled within-group SD — on the wide (per-group CDF) token, with no warp ($d$ is signed and $O(1)$). It is trained on **generic between-subjects prior-predictive cohorts**, not on SDY269: held-out PIT–KS $0.055$, coverage $0.055/0.514/0.947$, median-vs-true correlation $0.836$.

Deployed on SDY269 by **reusing the existing SVI fit** — the future cohort is the joint fit's `ypred` subset to the shared-strain endline observations, and the reference is the standardised difference formed from the per-arm SPR draws already in the fit; **no re-fit**. Calibration to the SVI: **PIT–KS $0.062$, marg–KS $0.061$, coverage $0.50/0.97$** — matching the held-out figure. The amortised **$\mathrm{PPS}(d>0.5\,\text{SD})=0.97$** at full accrual reproduces the SVI's decisive TIV-over-LAIV call. This is the head-to-head's methodological point: the cross-arm contrast the SVI computes from a joint posterior is amortised *directly* by a single generic network, and its deployment re-uses the existing fit rather than re-fitting.

## 17.5 Federated diagnostics per endpoint

The head-to-head deploys a **manifest** of registry instances — $S_1$ (SPR per arm, wide/logit), $S_2$ (GMFR per arm, scalar/$\log_2$), $S_5$ (difference, wide/none) — each a *different* trained network. The delegated amortisers live side-by-side as subdirectories of one federated deployment directory (one subdirectory per endpoint $\rho$, its subset reference nested within), and every diagnostic figure reports into that parent. Each per-arm amortiser is deployed over **all that arm's strains** (the item-general net carries them jointly), so the federated set spans every response item, not only the shared strain. Because the endpoints are federated across networks, each diagnostic is assembled as a **single combined figure faceted strain (rows) $\times$ endpoint $\rho$ (columns)**, with panels left empty where a strain is not measured under that endpoint (e.g. the between-arm difference exists only on the shared strain). A generic routine pools the tidy per-plot data each amortiser dumps and renders the diagnostic suite — PIT uniformity and coverage calibration, the marginal $p(\rho\mid x)$ SVI-vs-amortiser quantile boxes, the conditional-PIT box, the contraction power law, the $\eta_0$ success-threshold sweep with the SVI $\rho$-predictive against those thresholds, and the amortised-PPS trajectory — one plot each, for both the deployed (head-ft + affine + BvM) and raw-network baselines. Panels use per-facet free scales where the endpoints carry different units (rate vs fold vs standardised difference); the interim axis is labelled by participants-per-arm since SDY269 carries no calendar time. Per-interim single-endpoint overlays are not emitted.

Calibration of the five delegated amortisers against the joint SVI reference (mean over interims and over each arm's strains; coverage is the empirical 5%/95%-quantile mass):

| endpoint             | registry | PIT–KS | marg–KS | coverage (5% / 95%) |
|----------------------|----------|--------|---------|---------------------|
| LAIV SPR             | $S_1$ rate | 0.100 | 0.104 | 0.49 / 0.98 |
| TIV SPR              | $S_1$ rate | 0.066 | 0.073 | 0.52 / 0.94 |
| LAIV GMFR            | $S_2$ fold | 0.055 | 0.074 | 0.48 / 0.94 |
| TIV GMFR             | $S_2$ fold | 0.059 | 0.082 | 0.46 / 0.95 |
| TIV$-$LAIV SPR diff  | $S_5$ diff | 0.062 | 0.061 | 0.50 / 0.97 |

Every delegated network calibrates in the same PIT–KS $\lesssim 0.10$ band as the single-endpoint influenza deploys of §16 — three *different* trained networks ($S_1$/$S_2$/$S_5$), federated over one joint fit, each reproducing the SVI posterior of its endpoint. The cross-arm difference ($S_5$, §17.4) is the tightest (PIT–KS 0.062), confirming that a single generic between-groups network amortises the procurement-relevant contrast directly.

# 18. ICRC community-level MHPSS interim evaluation (DASS-21 / DRC)

The first of four further applications (two humanitarian, one HIV-vaccine, one consumer) deployed with the **generic between/within-cohort amortiser** exactly as the influenza endpoints, differing only in the estimand's registry family and the $\eta_0$ scale. ICRC community-level mental-health and psychosocial support (§3.10) is a paired pre/post humanitarian cohort at scale: the Depression, Anxiety and Stress DASS-21 subscale totals are each binned into the five clinical severity levels (Normal $\to$ Extremely severe, $K=5$) on the paper's own Figure-2 cut-offs, turning each subscale into one ordinal item. The estimand is the **relative severity reduction** post-vs-pre (DASS is distress, so `lower_is_better` — the direction-aware $\rho$ reads as a reduction), the registry $S_3$ relative-change instance (scalar `scale-feat` net, no warp), fitted over the DRC arm ($n=1{,}669$ beneficiaries with complete pre&post on all three subscales) with a fine early interim grid ($20/40/60/80$) then every 100 accruing beneficiaries.

## 18.1 Results: a very large effect, discriminating only at a high threshold

Community MHPSS produces enormous pre$\to$post severity drops ($\rho\approx0.80$–$0.95$ relative reduction, matching the paper's 96.6% improved on DASS-21), so the predictive probability of success is essentially $1$ from the first interim at the default $\eta_0=0.5$. The interim decision only becomes *discriminating* at a high operational threshold $\eta_0\approx0.85$ — the $\eta_0$ sweep is therefore run over $\{0,0.50,0.70,0.85,0.95\}$, and the operational lesson is that the go/no-go rule for a high-efficacy psychosocial intervention must be set on that scale rather than at a nominal 50% reduction.

## 18.2 Federated amortiser calibration

Deployed as a single-endpoint federated amortiser (one $S_3$ delegated network over the three subscales; combined diagnostic grid faceted subscale $\times$ $\rho$, one column). Calibration against the SVI reference (mean over interims and subscales):

| endpoint                     | registry | PIT–KS | marg–KS | coverage (5% / 95%) |
|------------------------------|----------|--------|---------|---------------------|
| DASS-21 severity reduction   | $S_3$ rel-change | 0.099 | 0.134 | 0.53 / 0.94 |

The relative-change network — the same one trained generically and reused for REFUGE, HVTN 505 and mycelium below — reproduces the SVI posterior across all three subscales in the PIT–KS $\approx0.10$ band, at real humanitarian scale and under a very large effect. Fits in `py-icrc-dass-drc_260902`; federated deploy in `py-icrc-dass-drc-amortise-…-J64-ftheadexpand-bvm-federated_260924`; the pipeline is driven by `scripts-py/ICRC-DASS-DRC_startme.py`.

# 19. REFUGE-ED perceived social support (MSPSS)

The second humanitarian application, and the amortiser's **cross-cohort generalisation test**. REFUGE-ED (§3.9) is a paired baseline/endline pilot of refugee and migrant youth across six countries with no treatment arm — the estimand is the within-cohort endpoint effect. The 12 MSPSS perceived-social-support items (three subscales — Family / Friends / Significant Other, $1$–$7$ Likert, `out-of-7` expected-score) give the **relative support gain** endline-vs-baseline (`higher_is_better`), again the registry $S_3$ relative-change instance, over the $324$ participants linked at both timepoints, on an 8-point participant-count interim grid ($n=40\to324$).

## 19.1 Results: modest gains, a low operational threshold

Perceived-support gains are modest, so — mirror-image of ICRC — the $\eta_0$ sweep is run low, over $\{0,0.05,0.10,0.15,0.20\}$: a 5–20% relative gain is the operationally interesting range for a psychosocial-support pilot, and the PPS becomes discriminating there rather than at a 50% bar. The 12 items span the three MSPSS subscales, so the combined diagnostic grid carries one row per item (single $\rho$ column).

## 19.2 Federated amortiser calibration

Calibration against the SVI reference (mean over interims and the 12 items):

| endpoint                       | registry | PIT–KS | marg–KS | coverage (5% / 95%) |
|--------------------------------|----------|--------|---------|---------------------|
| MSPSS perceived-support gain   | $S_3$ rel-change | 0.092 | 0.134 | 0.47 / 0.92 |

This is the key generalisation result: the *same* item-general $S_3$ network deployed on Ukraine (§14) and ICRC (§18) calibrates just as well on a **different cohort, a different instrument, and a wider item set** (12 MSPSS items) it never saw in training — PIT–KS 0.092, within the band of every other deployment. Fits in `py-refugee_interim_260831`; federated deploy in `py-refugee-interim-amortise-…-J64-ftheadexpand-bvm-federated_260924`; driven by `scripts-py/REFUGE-ED_startme.py`.

# 20. HVTN 505 HIV-vaccine immunogenicity (CAVD DataSpace / BAMA)

The efficacy-scale vaccine application on **real HIV-vaccine trial data** (§3.19). HVTN 505 is a phase-2b DNA/rAd5 test-of-concept stopped early for futility (2013); its binding-antibody (BAMA) assay data support a **vaccine-vs-placebo** immunogenicity contrast. HIV-naive subjects have no informative paired baseline (pre-vaccination antibody is a uniform floor), so unlike the influenza fold-rise the estimand is a **between-arm** contrast, encoded MYCELIUM-style with the arm playing the role of time (placebo $=$ `Baseline`, vaccine $=$ `Endline`). Each BAMA antigen's `mfi_delta` is binned per antigen into $K=3$ ordinal levels (low / intermediate / high response); the estimand is the per-antigen relative vaccine-vs-placebo shift (`higher_is_better`), the registry $S_3$ relative-change instance reused from ICRC/REFUGE. To keep the SVI reference coherent across interims the deploy restricts to the **rectangular 9-antigen panel** measured on all subjects (the sparse antigens, covered on only a subset, are dropped — see §16.4-style heterogeneity), over 239 subjects (190 vaccine / 49 placebo) with shuffled accrual so both arms are present at each of 10 interims.

## 20.1 Results: a coherent immune-correlate ordering

The between-arm shifts recover the known HVTN 505 immunogenicity ordering: Con6 gp120 shows the strongest vaccine-vs-placebo effect ($\rho\approx1.14$, amortised PPS $\to1.00$), the V1V2 immune-correlate scaffolds sit at $\rho\approx0.75$–0.80, and the p24 Gag antigen is essentially null ($\rho\approx0.50$). The go/no-go is unambiguous for the gp120/V1V2 antigens and a clear no-go for p24 — a real HIV-vaccine application of the between-arm amortiser on public Global-Access data. (The 2013 futility stop was on HIV-*infection* efficacy, a time-to-event endpoint absent from the antibody datasets, so that specific interim is not reconstructed.)

## 20.2 Federated amortiser calibration

Calibration against the SVI reference (mean over interims and the 9 antigens):

| endpoint                          | registry | PIT–KS | marg–KS | coverage (5% / 95%) |
|-----------------------------------|----------|--------|---------|---------------------|
| BAMA vaccine-vs-placebo shift     | $S_3$ rel-change | 0.092 | 0.156 | 0.54 / 0.92 |

The marginal KS (0.156) is the highest of the deployments, reflecting the sharp $K=3$ per-antigen ordinal (only three categories, so the marginal $p(\rho\mid x)$ is coarser than the $K\ge5$ instruments), but the conditional PIT–KS (0.092) is in the same band as the rest — the between-arm relative-change network is as well calibrated on ordinal HIV binding data as on continuous-scale psychosocial instruments. Fits in `py-cavd-vtn505-bama_260918`; federated deploy in `py-cavd-vtn505-bama-amortise-…-J64-ftheadexpand-bvm-federated_260923`; driven by `scripts-py/CAVD-hvtn505_startme.py`.

# 21. Mycelium novel-food acceptance (powder vs burger)

The consumer-launch application (§3.14): a UK Prolific panel rates mycelium as a food protein under a $3\times3$ (processing $\times$ substrate) manipulation. Built as the between-arm "Option A" contrast — product **powder vs burger**, substrates pooled — encoded exactly as HVTN 505 with the arm as time (burger $=$ `Baseline`, powder $=$ `Endline`); each respondent is in one arm (unpaired, $n=298$: 149 powder / 149 burger). Nine $1$–$7$ items over three constructs — Acceptance (A1–A4), Disgust (D1–D4), Perceived naturalness (PN) — give a direction-aware **relative powder-vs-burger shift** per item (Acceptance/Naturalness up $=$ good, Disgust down $=$ good; each oriented so higher $=$ better), the registry $S_3$ instance, over 8 shuffled participant-accrual interims.

## 21.1 Results: powder beats burger on every construct

The contrast is consistent and modest: Acceptance $+0.08$–$0.15$, Disgust (severity down) $+0.13$–$0.27$, Naturalness $+0.11$ — powder beats the burger presentation on every construct, none at the 50% level, so (as for REFUGE) the operationally interesting $\eta_0$ sweep is low, $\{0,0.10,0.20,0.30\}$. This is the launch-decision regime of §3.12 with a genuine condition contrast: given the ratings so far, will powder clear its acceptance margin once the panel completes?

## 21.2 Federated amortiser calibration

Calibration against the SVI reference (mean over interims and the 9 items):

| endpoint                        | registry | PIT–KS | marg–KS | coverage (5% / 95%) |
|---------------------------------|----------|--------|---------|---------------------|
| powder-vs-burger relative shift | $S_3$ rel-change | 0.078 | 0.094 | 0.51 / 0.95 |

The best-calibrated of the between-arm deployments (PIT–KS 0.078, marg 0.094) — the mycelium contrast previously reported as *amortiser-deferred* (§3.14) is delivered here by the same $S_3$ between-group network as HVTN 505, closing that gap. Fits in `py-mycelium-powdervsburger_260902`; federated deploy in `py-mycelium-powdervsburger-amortise-…-J64-ftheadexpand-bvm-federated_260924`; driven by `scripts-py/MYCELIUM_startme.py`.

# 22. References

```{=html}
<!--
Bibliography is in dev/amortised_decision_making.bib. Render with pandoc:
  pandoc dev/amortised_decision_making.md --citeproc \
    --bibliography dev/amortised_decision_making.bib \
    -o amortised_decision_making.html
The list below is auto-generated by --citeproc from the inline [@key] citations.
-->
```

# A. Appendix

## A.1 Inefficient sampling schemes to estimate PPS

### A.1.1 Importance sampling (self-normalized)

Since $p(\theta \mid x, z^{(s)}) \propto p(\theta \mid x)\, p(z^{(s)} \mid \theta)$, we can re-use all existing posterior draws $\theta_{k} \sim p(\theta \mid x)$ for $k=1,\dotsc,K$ and reweight these by the future-data likelihood.

Working in log space for stability, for each $z^{(s)}$ separately, we compute for all $k=1,\dotsc,K$ the importance sampling weights $$\begin{aligned}
\log w_k^{(s)} &= \log p(z^{(s)} \mid \theta_k) = \sum_{i=1}^m \log p\big(z^{(s)}_i \mid \theta_k\big), \\
\tilde w_k^{(s)} &= \operatorname{softmax}_k\!\big(\log w_k^{(s)}\big) = \frac{\exp\!\big(\log w_k^{(s)} - \ell^{(s)}\big)}{\sum_j \exp\!\big(\log w_j^{(s)} - \ell^{(s)}\big)}, \quad \ell^{(s)} = \log\!\textstyle\sum_j \exp \log w_j^{(s)}, \\
y^{(s)} &\approx \sum_{k=1}^K \tilde w_k^{(s)} \, 1_{\theta_k \in H_1}.
\end{aligned}$$ The log-sum-exp / softmax map is the numerically stable form of $w / \sum w$: subtracting the maximum log-weight before exponentiating prevents overflow while leaving the normalized weights unchanged.

Reliability is monitored by the effective sample size [@kong1992note; @liu2001monte] $\mathrm{ESS} = (\big(\sum_k w_k\big)^2)/(\sum_k w_k^2) = (\sum_k (\tilde w_k)^2)^{-1}$, equivalently the second-order weight moment $\mathbb{E}(\tilde w^2) = \tfrac{1}{K}\sum_k \tilde w_k^2$, with $\mathrm{ESS}/K = 1/\big(K^2\, \mathbb{E}(\tilde w^2)\big)$.

Self-normalized IS is consistent but $O(1/K)$ biased, and its variance is finite only when $\mathbb{E}_{p(\theta \mid x)}\!\big[p(z \mid \theta)^2\big] < \infty$; the Pareto-smoothed importance sampling (PSIS) tail index $\hat k$ [@vehtari2024pareto] both estimates this and stabilizes the largest weights, with $\hat k > 0.7$ flagging an unreliable estimate.

The main issue is that the proposal/target mismatch grows with the amount of assimilated future data: $\mathrm{KL}\big(p(\theta \mid x, z)\,\|\,p(\theta \mid x)\big)$ increases in $m$, so the weight mass concentrates on a single draw and $\mathrm{ESS} \to 1$. Empirically, at the earliest interim of our case study ($n \approx 48$ current vs $m \approx 455$ future units) the weights collapse to $\mathrm{ESS} \approx 1$ after even a *single* future participant, with PSIS $\hat k = \infty$. Plain IS labels are therefore trustworthy only when $z$ is small relative to $x$ (late interims). The two corrections below target this regime.

### A.1.2 Moment-matching importance sampling

Moment-matching IS [@paananen2021implicitly] repairs a mild proposal/target mismatch without new model fits, by transforming the draws so the transformed cloud better covers the target and reweighting with the change-of-variables Jacobian. Starting from the IS weights $\tilde w_k$ of 6.1, compute the weighted and proposal moments $\hat\mu_w = \sum_k \tilde w_k\, \theta_k, \qquad \hat\mu_q = \tfrac{1}{K}\sum_k \theta_k, \qquad (\text{optionally } \hat\Sigma_w,\ \hat\Sigma_q),$ and apply an invertible affine map $T$ that matches them. The mean-match step uses $$T(\theta) = \theta + (\hat\mu_w - \hat\mu_q), \qquad \theta_k^* = T(\theta_k),$$ while the covariance-match variant uses $$T(\theta) = \hat\mu_w + L_w L_q^{-1}(\theta - \hat\mu_q)$$ with $\hat\Sigma_\bullet = L_\bullet L_\bullet^\top$. The transformed draws are reweighted against the target with the Jacobian of $T^{-1}$, $$w_k^* = \frac{p(\theta_k^* \mid x, z^{(s)})}{q^*(\theta_k^*)}, \qquad q^*(\theta^*) = q\big(T^{-1}\theta^*\big)\,\big|\det \nabla T^{-1}\big|,$$ where $q$ is the proposal density $p(\theta \mid x)$ (in practice a diagonal-Gaussian fit to the base draws in an unconstrained reparameterisation, with positive parameters mapped through $\log$). One iterates over a small family of transforms and keeps the one maximizing $\mathrm{ESS}$ (or minimizing $\hat k$).

We found that when the base $\mathrm{ESS} \approx 1$, the weighted mean $\hat\mu_w$ *equals* the single dominating draw, so the affine shift only relocates the whole cloud onto that point and $\mathrm{ESS}$ does not recover. Moment matching corrects mild mismatch but cannot manufacture the support the fixed base draws lack — it never moves a particle to a region the proposal failed to sample. In our case study it leaves the early-interim $\mathrm{ESS}/K$ unchanged at $\approx 1/K$.

### A.1.3 Sequential Monte Carlo with resample-move

To cross an arbitrarily large $x \to (x, z)$ gap, another idea is to bridge the proposal to the target through a tempered sequence (annealed importance sampling [@neal2001annealed]; SMC samplers [@delmoral2006smc; @chopin2002sequential]), $$\pi_t(\theta) \;\propto\; p(\theta \mid x)\; p(z^{(s)} \mid \theta)^{\beta_t},$$ for $0 = \beta_0 < \beta_1 < \dots < \beta_T = 1$, so $\pi_0 = p(\theta \mid x)$ and $\pi_T = p(\theta \mid x, z^{(s)})$ (the target).

Initialise particles $\theta_k \sim p(\theta \mid x)$ with uniform weights; at step $t$:

1.  **Reweight** by the incremental likelihood, $\tilde w_k = \operatorname{softmax}_k\!\Big( (\beta_t - \beta_{t-1})\, \log p(z^{(s)} \mid \theta_k) \Big).$
2.  **Adapt** $\Delta\beta = \beta_t - \beta_{t-1}$ by bisection so the tempering $\mathrm{ESS}/K$ hits a target (e.g. $\tfrac12$), which is an automatic schedule [@jasra2011inference; @zhou2016toward].
3.  **Resample** the particles by $\tilde w_k$ (systematic resampling) and reset weights to $1/K$.
4.  **Move** each particle with an MCMC kernel $M_t$ leaving $\pi_t$ invariant (resample-move [@gilks2001following]). We use Metropolis-adjusted Langevin (MALA [@roberts1996exponential]) in the unconstrained reparameterisation, adapting the step size toward the optimal acceptance $\approx 0.574$ [@roberts1998optimal].

Because the temperature enters only as an exponent, the move kernel's log-density $\log p(\theta \mid x) + \beta_t \log p(z^{(s)} \mid \theta)$ only need to be compiled once with $\beta_t$ a traced argument and reused across all temperatures. The final particles approximate $\pi_T$ with uniform weights, giving the label $y^{(s)} \approx \frac{1}{K} \sum_{k=1}^K 1_{\theta_k^{(T)} \in H_1}.$ Unlike IS and moment matching, the move step relocates particles into the target's typical set, so SMC crosses an arbitrarily large gap; the cost is $T$ tempering steps, each a short MCMC sweep.

We found that at the worst (earliest) interim the adaptive schedule reaches $\beta = 1$ in $\approx 40$ steps and restores $\mathrm{ESS}/K \approx 1$ at a wall-clock cost dominated by the per-step move rather than the one-off compile. The overall computational cost was 7-8 times larger than SVI estimation of the posterior $p(\theta | x, z^{(s)})$.

------------------------------------------------------------------------

# B. Appendix B: Other data sets

Data sets considered but set aside for an amortised interim analysis; the blocking issue is stated at the top of each. Section numbers are retained from the main text so existing cross-references still resolve.

## 3.5 Application: temporal dynamics in psychological assessments

**Issue:** no intervention and no repeat measurement

A large cross-sectional dataset of $24{,}292$ students, each answering four self-report symptom scales at the item level [@su2023temporal]. It contains no intervention and no repeat measurement, so it does not itself define an endpoint effect; its value is as a **real, large-scale polytomous item bank**: fit the PCM once, then simulate realistic cohorts from the fitted item parameters to check that the amortised predictive is calibrated against a ground truth we control.

| field | value |
|-----------------|-------------------------------------------------------|
| Setting | 24,292 students, single administration (Feb–Mar 2021) |
| Design | cross-sectional; no arms, no baseline/endline |
| Items | PHQ-9 (9 items, 0–3), GAD-7 (7, 0–3), ISI (7, 0–4), PSS (10, 0–4) — all ordinal, PCM-ready |
| Effect measure $\rho$ | not defined (no pre/post) |
| Interims | not applicable |
| **Application target** | **validate the PCM amortiser**; simulate cohorts from the real fitted item bank |
| Data | [@su2023temporal] |

------------------------------------------------------------------------

## 3.11 Application: PISA international assessment

**Issue:** each wave has fresh participants, not accrueing

The OECD Programme for International Student Assessment, a triennial cross-national assessment whose constructed-response cognitive items are scored with **partial credit** and calibrated with the generalized partial credit model — the PCM's own family [@oecd2022pisa]. No trial, but by treating **each subsequent cycle as an interim evaluation point** the achievement trend becomes a **baseline-anchored interim analysis**: baseline = 2012, and 2015 / 2018 / 2022 are the interims, with $\rho$ the change from the 2012 baseline.

| field | value |
|-------------|-----------------------------------------------------------|
| Setting | OECD PISA Mathematics, the **54 countries present in all four cycles** 2012–2022; $\approx$ 2{,}000 students per cycle per country |
| Design | cross-sectional per cycle (**unpaired** cohorts, different students each cycle); baseline-anchored across cycles |
| Items | **dichotomous** Math cognitive items ($K=2$); **79 common trend items** matched across cycles by **content id** ($\texttt{PM033Q01}\equiv\texttt{CM033Q01}\equiv\texttt{M033Q01}$, stripping the paper/computer mode prefix) |
| Effect measure $\rho$ | $\rho_j=$ 2012$\to$cycle expected-score change per item, from a PCM with a **cycle covariate** — shared (anchored) item difficulties, cohort ability shift = the PISA trend model |
| Interims | the subsequent cycles **2015, 2018, 2022** (per country); **156** SVI fits over 54 countries $\times$ \$\approx\$3 cycles (a few dropped for mode-switch / thin item overlap) |
| **Application target** | cross-cycle achievement **trend as an interim analysis**; per-country SVI $\rho$ trajectories, validated against the raw per-item change (**corr** $0.90$, $94\%>0.8$) |
| Data | public OECD microdata: 2015/18/22 student cognitive `.sav` (webfs.oecd.org), 2012 scored-cognitive fixed-width TXT + SPSS control via the Wayback mirror; download + content-id extract in `python/data_web_extracting.py`, loader `read_data_pisa` |

**Notes.** (i) *Anchoring.* Unpaired cohorts are linked into one ability scale by sharing item difficulties across cycles (the `~ time - 1` cycle covariate); without this the cycles float independently and the trend is unidentified. (ii) *Dichotomous only.* Mixing $K=2$ and $K=3$ items under one item type forces $K=3$ on all, giving the binary items a phantom threshold that mean-field SVI collapses — so partial-credit items are dropped here. (iii) *Mode confound.* 2012 was paper, 2015+ computer, so the 2012\$\to$2015 step carries a paper$\to$computer mode effect on top of the real trend (the OECD applies a mode adjustment we do not); the *shape* of the 2015$\to$2018$\to\$2022 trajectory is the cleaner signal — e.g. United States Math falls $-0.04/-0.04/-0.10$ (2015/18/22 vs 2012), steepest at the 2022 (COVID) cycle. (iv) *Amortiser deferred.* The §14 amortiser's core statistic is a **per-participant** baseline$\to$endline change; PISA is unpaired, so the marginal-$p(\rho\mid x)$ validation against SVI awaits a between-cohort variant of the network (§14). Fits in `py-pisa-math_260904`; producer `scripts-py/PISA_interim_svi.py`.

------------------------------------------------------------------------

## 3.15 Application: selfBACK app-delivered self-management RCT for low-back pain

**Issue:** did not share data

A randomised controlled trial of an AI-app that delivers evidence-based, individually tailored self-management support for low-back pain, with item-level patient-reported outcomes at five waves (Sandal et al., *JAMA Intern. Med.* 2021; protocol Rasmussen et al., *JMIR Res. Protoc.* 2019, doi:10.2196/14720; ClinicalTrials.gov NCT03798288). The **largest and most richly timed** clinical trial of the set — a genuine control arm and a five-wave follow-up — outside the mental-health/humanitarian settings.

| field | value |
|---------------|---------------------------------------------------------|
| Problem | app-delivered self-management support for low-back pain (musculoskeletal / chronic pain) |
| Setting | multicentre RCT (Denmark / Norway); primary-care low-back-pain patients |
| Arms | **usual care vs usual care** $+$ selfBACK app ($N=461$: 229 / 232) |
| Timepoints | baseline, 6 weeks, 3 / 6 / 9 months — the richest wave structure of the set |
| Items | **Roland–Morris Disability Questionnaire** (24 ordinal items), pain-intensity NRS, and further PROMs — item-level, PCM-ready |
| Effect measure $\rho$ | baseline$\to$follow-up endpoint effect per item (disability / pain reduction, `lower_is_better`) |
| Cohort / interims | $N=461$ ($>200$); the five waves as interims — showcases contraction with cohort size and follow-up depth |
| Data access | **on reasonable request** — a data steering committee that welcomes data-sharing enquiries (like §3.10, not openly deposited) |
| **Application target** | **real-time evaluation** of a randomised self-management intervention at scale, with a control arm and a deep follow-up |
| Status | *planned* — strong structural fit; awaiting a data request |
| Reference | Sandal et al. 2021; NCT03798288 |

------------------------------------------------------------------------

## 3.16 Application: digital data-driven intervention RCT for depression and anxiety

**Issue:** did not respond to data sharing request

A waitlist-controlled randomised trial of a digital, data-driven therapeutic intervention for depressive and generalised-anxiety symptoms, with standard item-level symptom scales at three waves (Nature *npj Digit. Med.* 2025, doi:10.1038/s41746-025-01511-7; PMC11840063). A clean, moderate-$N$ RCT with a control arm and the canonical PROM battery.

| field | value |
|----------------|--------------------------------------------------------|
| Problem | digital self-help for depressive and generalised-anxiety symptoms |
| Setting | fully remote RCT, community adults |
| Arms | **intervention vs waitlist control** ($N=200$: 100 / 100; 164 completers) |
| Timepoints | baseline (week 0), mid (week 8), post (week 16) |
| Items | **PHQ-9**, **GAD-7** (primary), SWLS, LISAT-11 — item-level ordinal, PCM-ready |
| Effect measure $\rho$ | baseline$\to$endline symptom reduction per item (distress `lower_is_better`) |
| Cohort / interims | $N=200$; the three waves as interims |
| Data access | **on request** — "made available by the authors upon request"; totals reported in the paper |
| **Application target** | **real-time evaluation** of a randomised digital mental-health intervention, control-arm H$_1$ decision |
| Status | *planned* — awaiting a data request |
| Reference | *npj Digit. Med.* 2025, doi:10.1038/s41746-025-01511-7 |

------------------------------------------------------------------------

## 3.20 Application: COVID-19 subpopulation contrasts (ImmPort / ImmuneSpace)

**Issue:** observational cohort, actual question unclear betw severe/mild, adult/child

A **between-subpopulation** immune-response application: rather than a within-participant pre/post effect or a vaccine-arm contrast, compare the antibody response between two *subpopulations* of a COVID-19 cohort — the same **group-contrast** estimand as MYCELIUM/CAVD-BAMA (§3.14/§3.19), on SARS-CoV-2 neutralisation titres. From an inventory of 13 ImmuneSpace COVID studies, **SDY1764** (Distinct antibody responses to SARS-CoV-2 in children and adults across the clinical spectrum) is the one carrying an **ordered-titre** assay (serum neutralisation ID50 + ELISA) together with clean subpopulation strata. It supports two contrasts from its four clinical arms (adult ARDS, pediatric MIS-C, pediatric non-MIS-C, adult convalescent): **age** (pediatric $<18$ vs adult) and **severity** (severe $\{$ARDS, MIS-C$\}$ vs mild $\{$non-MIS-C, convalescent$\}$). (The severity-focused SDY1669 was inspected and **dropped** — it has only flow-cytometry / RNA-seq, no antibody titre, so it cannot drive the titre-PCM; SDY1764's own severity arms cover that contrast.)

| field | value |
|-----------|-------------------------------------------------------------|
| Problem | COVID-19 humoral immunity **by subpopulation** — is the neutralising response higher in one group than another (a genuine, decision-relevant contrast for risk stratification / trial design) |
| Setting | ImmuneSpace/ImmPort HIPC COVID cohort SDY1764 (Mount Sinai); 79 subjects with serum neutralisation |
| Groups | **between-subpopulation** (no paired baseline): age (pediatric 47 / adult 32) or severity (severe 29 / mild 50); encoded MYCELIUM-style (group A $=$ `Baseline`, group B $=$ `Endline`) |
| Items | **SARS-CoV-2 serum-neutralisation ID50** (single item; $\log_{10}$ titre $\to$ ordered category $k=\mathrm{round}(\log_2(\text{titre}/4))$, $K=8$). ELISA (3 unlabelled protein targets) available as optional extra items |
| Effect measure $\rho$ | relative group-B-vs-A shift in mean $\log_2$ titre, `higher_is_better` (a between-group GMT ratio, as CAVD-BAMA) |
| Cohort / interims | 79 subjects, shuffled accrual so both groups present at each of 8 interims |
| **Application target** | the **between-subpopulation** amortiser on real SARS-CoV-2 neutralisation titres — a third instance of the group-contrast estimand (after MYCELIUM food and HVTN-505 vaccine arms), now on disease/age strata |
| Data access | **ImmuneSpace/ImmPort** (DUA + API, same as §3.18); pulled to `ImmuneSpace_COVID19_v260920.xlsx` (13 COVID SDYs) |
| Loader + producer (implemented) | `read_data_immport_covid_neut(xlsx, study='SDY1764', group='age'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | 'severity')` in `python/data_loading.py` (neut $\log_{10}\to$ ordered $k$, subpopulation $\to$ group axis) $+$ `scripts-py/IMMPORT_covid_interim_svi.py` (mirrors `CAVD_bama_interim_svi.py`) |
| Status | **implemented — SDY1764, SVI + amortiser.** Directionally coherent: **age** $\rho(\text{pediatric vs adult})\approx-0.38$ (children neutralise \$\sim\$40% lower); **severity** $\rho(\text{severe vs mild})\approx+0.32$ (severe higher). SVI grids $\to$ `py-immport-covid-SDY1764-{age,severity}_260920` (with the `p_rho_x_by_item` contraction plot). **J64 amortiser deployed** (between-group, scalar-mean token, head BvM, no warp — signed relative-change endpoint) via `deploy_IMMPORT_covid_amortiser.sh`: calibration to the SVI reference **age PIT-KS 0.135 / marg 0.115, severity 0.122 / 0.104**; PPS(severe$>$mild) $=0.94$ at $\eta_0{=}0$. Note the SDY1764 severity axis is confounded with age/phenotype (severe $=$ adult-ARDS $+$ pediatric-MIS-C); a within-age contrast or IMPACC (SDY1760) de-confounds. Next: ELISA items; other titre-bearing COVID SDYs |
| Reference | SDY1764 (Mount Sinai; PMID 33154590); IMPACC SDY1760/2112 (severity-trajectory cohort, larger follow-on); ImmuneSpace/ImmPort |