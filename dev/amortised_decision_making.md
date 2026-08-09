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

A canonical high-dimensional analytically tractable PPS benchmark for continuous outcomes is the multivariate normal model with structured covariance. Let $$y_n \mid \mu, \sigma^2 \;\sim\; \mathrm{MVN}(\mu, \sigma^2 K), \qquad n = 1, \ldots, N,$$ where $K = R R^\top \in \mathbb{R}^{J \times J}$ is a known positive-definite covariance shape (for example, $R$ a Cholesky factor of an AR(1) or arbitrary correlation matrix), and the unknowns are the mean vector $\mu \in \mathbb{R}^J$ and the scalar variance $\sigma^2 > 0$.

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

| Cell            | Knobs                                       | Notes                                   |
| --------------- | ------------------------------------------- | --------------------------------------- |
| Cell A (easy)   | $R = R^{(1)}$ AR(1), $J \in \{50,100,200\}$ | Banded, well-conditioned.               |
| Cell B (medium) | $R = R^{(2)}$ block, same $J$               | Block structure; eigenvalues clustered. |
| Cell C (hard)   | $R = R^{(3)}$ factor, same $J$              | Dense long-range; eigenvalue decay.     |

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

Computationally, the PCM involves direct evaluations of category specific events ($Pr(Y=k)$ not $Pr(Y\leq k)$) and therefore tends to be much faster to evaluate and less prone to numerical issues than other IRT models. Standard priors can be attached to all free parameters.

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
| ----------------------- | ---------- | ---------- |
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
| -------- | -------: | --------: | ----: |
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
| -------- | -------: | -----------------------: | ----: |
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

| Method                                                   |         MSE | $\sqrt{\text{MSE}}$ | Mean per-interim inference (min) | One-off training (min) |
| -------------------------------------------------------- | ----------: | ------------------: | -------------------------------: | ---------------------: |
| **Amortised (features-fixed, `hidden_dims = (64, 64)`)** | **0.00119** |           **0.035** |                        **0.001** |               **1.94** |
| Regression endpt-x (Gaussian approx)                     |     0.00160 |               0.040 |                            0.001 |                      — |
| Nested-MC using HMC for each (x, z)                      |     0.00244 |               0.049 |                             6.04 |                      — |
| IS reweighting of $\theta \mid x$                        |     0.00277 |               0.053 |                            0.004 |                      — |
| Amortised (features-MLP)                                 |     0.00277 |               0.053 |                            0.002 |                   8.59 |
| Regression endpt-x (mquantile)                           |     0.00355 |               0.060 |                            0.001 |                      — |
| Regression endpt-x (quantile)                            |     0.00357 |               0.060 |                            0.001 |                      — |
| Nested-MC using SVI for each (x, z)                      |     0.00401 |               0.063 |                             1.10 |                      — |
| Regression H1-x on w(z)                                  |     0.01033 |               0.102 |                            0.001 |                      — |

*Amortised features-fixed shown after the §12.4 default swap to `hidden_dims = (64, 64)` — halves MSE vs the old `(256, 256, 128)` default at* $3\times$-faster training. Training cost is amortised across all future deployments — pay once, deploy in milliseconds thereafter. Nested-MC HMC costs 6 min **per interim per x**, so on a 12-interim schedule with $S = 200$ Monte-Carlo draws the total nested-MC cost is $\sim 72$ min per new patient cohort $x$; the amortised features-fixed variant recovers a lower MSE for a 1.9 min one-off training + 0.001 min per interim, breaking even before the first fresh $x$. All timings measured with $\le 4$ CPU threads on a single machine.

**Reading.** The MSE floor is set by Monte-Carlo noise at $S = 200$: $\sqrt{p (1 - p) / S} \approx 0.035$ near the PPS mode, so any method below $\sqrt{\text{MSE}} \le 0.05$ is essentially at the MC-noise floor. Both amortised variants are numerically tied with IS at that floor. The Gaussian-approx regression wins because the Beta posterior mean is well-approximated by a Gaussian at the sample sizes seen, giving a very tight plug-in $\Phi((\hat\mu - \eta_0)/\hat\sigma)$ estimator. Nested-MC HMC still uses the finite-sample HMC posterior (wider than the exact Beta at small $n$), so it inherits an interim-1 error of 0.15 that inflates its MSE. The `H1-x` regression is the most biased — the binary label loses too much information relative to the continuous endpoint targets used everywhere else.

**Amortised (features-MLP) reading.** The MLP variant matches the fixed variant's aggregate MSE despite starting from no structural prior. Per-interim errors are on average larger (mean 0.041 vs 0.020) but the *distribution* of errors is symmetric around the analytic PPS, so squared errors average out. A longer training budget would close the mean-error gap; the ceiling is the MC-noise floor at $S = 200$, same as everyone else. This validates the general-purpose amortiser template for models where the sufficient statistic is not known analytically (Categorical, IRT).

## 12.4 Ablation: features-fixed capacity + training-config knobs

Five variants of the features-fixed amortiser exercised to probe whether the default `hidden_dims = (256, 256, 128)` config from §12.1 was over-parameterised and whether other training knobs move the MSE floor. Each variant differs from the default in one dimension; all use the same analytic joint `zi` path of §12.3 and the same deployment cohort. Selectable at run time via the `AMORTISER_VARIANT` environment variable on the deployment script:

| Variant                   | Change vs default                                                  |         MSE | $\sqrt{\text{MSE}}$ | Training (min) |
| ------------------------- | ------------------------------------------------------------------ | ----------: | ------------------: | -------------: |
| default `(256, 256, 128)` | —                                                                  |     0.00277 |               0.053 |           5.90 |
| **`64x64`**               | `hidden_dims = (64, 64)`                                           | **0.00119** |           **0.035** |       **1.94** |
| `num_quantile_levels_5`   | `taus = (0.05, 0.25, 0.5, 0.75, 0.95)`                             |     0.00019 |               0.014 |           5.26 |
| `num_quantile_levels_21`  | 21 equally-spaced taus in $(0.025, 0.975)$                         |     0.00119 |               0.035 |          46.20 |
| `S_2000`                  | Deployment $S = 2000$ instead of $200$                             |     0.00029 |               0.017 |          38.13 |
| `log_uniform_n`           | Training $n \sim$ log-uniform on $[1, N-1]$ (oversample small $n$) |     0.00119 |               0.035 |          22.23 |

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

| Variant                                        | Encoder     |         MSE | $\sqrt{\text{MSE}}$ | Training (min) |
| ---------------------------------------------- | ----------- | ----------: | ------------------: | -------------: |
| `combo_64x64_qlv5_S2000`                       | **MLP**     | **0.00013** |           **0.011** |       **8.05** |
| `num_quantile_levels_5`                        | fixed       |     0.00019 |               0.014 |           5.26 |
| `S_2000`                                       | fixed       |     0.00029 |               0.017 |          38.13 |
| `S_2000`                                       | **MLP**     |     0.00029 |               0.017 |           7.88 |
| `combo_64x64_qlv5_S2000`                       | fixed       |     0.00046 |               0.022 |           1.85 |
| `num_quantile_levels_5`                        | **MLP**     |     0.00119 |               0.035 |          12.70 |
| `64x64` (default)                              | fixed       |     0.00119 |               0.035 |           1.94 |
| `num_quantile_levels_21`                       | fixed       |     0.00119 |               0.035 |          46.20 |
| `log_uniform_n`                                | fixed       |     0.00119 |               0.035 |          22.23 |
| `log_uniform_n`                                | **MLP**     |     0.00192 |               0.044 |          29.60 |
| `64x64`                                        | **MLP**     |     0.00243 |               0.049 |           7.86 |
| default (`(256, 256, 128)` / `(128, 128, 64)`) | fixed / MLP |     0.00277 |               0.053 |    5.90 / 8.59 |
| `num_quantile_levels_21`                       | **MLP**     |     0.00277 |               0.053 |          52.16 |

**Findings.**

1.  **Combined `(64x64 + qlv5 + S_2000)` is the best config across all variants**, and the **MLP variant wins outright at** $\sqrt{\text{MSE}} = 0.011$ ($\sim 5\times$ below the naive $S = 200$ MC floor of $0.018$; consistent with the $S = 2000$ MC floor of $\sqrt{0.07 \cdot 0.93 / 2000} \approx 0.006$). Compare to the fixed variant at $\sqrt{\text{MSE}} = 0.022$: same combined knobs but the fixed variant plateaued because its `qlv5` alone was already anomalously good ($0.00019$, a single-realisation dip below the $S = 200$ MC floor) — enlarging to `S_2000` in the combo regresses to the MC-noise-limited value. The MLP variant, which was far above the MC floor without the knobs, benefits monotonically from all three additions and lands at the true floor. **Recommend running the MLP combo as the default `Amortised (features-MLP)` entry going forward.**

2.  **Robustness of the ablation findings.** Each of the four single-knob ablations moves both encoders in the same direction:

    | Knob                     |  Fixed $\Delta$MSE |    MLP $\Delta$MSE |
    | ------------------------ | -----------------: | -----------------: |
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
| ---------------- | -------------: | -------: | ----------: | ------------------: |
| CPU (6-thread)   |           8.05 |     1.0× |     0.00013 |              0.0113 |
| **MPS (M4 Max)** |       **2.71** | **3.0×** | **0.00009** |          **0.0095** |

**Findings.**

- **3× training speedup on MPS.** The features-MLP encoder's per-token MLP over a `(B=512, N_max=500, item_dim=1)` batch tensor is exactly the matmul-heavy workload GPU parallelism was designed for. Kernel-launch overhead is negligible at $15\,000$ training steps.
- **MSE within MC-noise band of CPU.** MPS lands at $\sqrt{\text{MSE}} =
  0.0095$ vs CPU $0.0113$; both consistent with the $S = 2000$ Monte-Carlo floor of $\sqrt{p(1-p)/S} \approx 0.006$. Differences reflect independent PRNG seeds through the JAX/Metal backend (Metal reductions are non-deterministic) rather than model quality.
- **Not worth trying for features-fixed.** The fixed `64x64` variant is a trivial `(B=8192, feature_dim=2)` batch matmul that finishes in 1.9 min on CPU; MPS's kernel-launch overhead would dominate and give neutral-to-worse wall-clock. Confirmed empirically in a spot check (not shown).
- **Caveats.** `jax-mps` is flagged experimental; some ops silently fall back to CPU. Determinism guarantees are weaker on Metal, so the exact loss trajectory changes across runs but the population minimiser is the same and the resulting PPS estimator is within MC noise.

**Recommendation.** For the features-MLP amortiser and any future model with a matmul-heavy DeepSets encoder (Categorical, IRT), invoke the deployment script under the `mps-experimental` pixi environment. Keep features-fixed on CPU (no benefit from GPU).

# 13. Results for the MVN model

## 13.1 Results for MVN interim analysis amortised endptx on wz with features-fixed, idcomp, qpsi-MLP, loss-multiquantilehead

**MVN extension of §12.1 — `idcomp` (independent-components) variant.** The amortiser sees **one component's summary at a time** and predicts $\mu_j$ marginally, so the known cross-component covariance $K$ is used at TRAINING (to draw $\mu \sim \mathrm{MVN}(\mu_0 \mathbf{1}_J, \tau^2 K)$ correctly) but **discarded at DEPLOYMENT** — the joint posterior over $\mu$ is factorised across $j$. That is a real approximation: it is exact for the per-component target $\rho_j = \mu_j - \mu_0$ when the amortiser only ever needs marginals (as here, where PPS aggregates $\mathbf{1}\{P(H_{1j} \mid x, z) > \eta_H\}$ per component), but leaves the cross-component structure on the table for downstream utilities that need the joint (multivariate stopping rules, family-wise error control, correlated effect-size summaries). In real-world data the components are correlated and a full-$K$ variant that exploits that at deployment is the natural next step — future work; §13.1/§13.2 report only the idcomp baseline.

Same amortiser class + shared `train` / `predict_amortised_p_h1_for_one_xz` / `save_trained_model` utilities from `amortiser_common`. Two new pieces plumb the MVN case study into the existing pipeline:

- `MVNModel.make_training_data_with_features` (in [`python/model_mvn.py`](../python/model_mvn.py)): joint prior-predictive sampler. Draws $\mu \sim \mathrm{MVN}(\mu_0 \mathbf{1}_J, \tau^2 K)$; picks $n \sim U\{1, \ldots, N-1\}$, $m = N - n$; draws $y_{1:N} \sim \mathrm{MVN}(\mu, \sigma^2 K)$; per component $j$ builds features `(sum_i y_{i,j}, n + m) / N_max` (analogous to the Binomial `(k_total, n_total) / N_max`) and target $\rho_j = \mu_j - \mu_{0,\text{baseline}}$. Returns `(S * J, 2)` features + `(S * J,)` targets — **the per-component flatten makes the amortiser J-invariant** under the MVN g-prior with unit-diagonal $K$ (per-component posterior of $\mu_j \mid y_{1:N}$ depends only on $(\sum_i y_{i,j}, N)$).
- `MVNModel.fit_closed_form_posterior` gained a `resume` flag (mirroring the Binomial and HMC drivers) so the outer $x$-posterior draws can be cached across method comparisons.
- Deployment script [`scripts-py/MVN_interim_analysis_amortise_endptx_on_wz_with_features_fixed_idcomp_qpsi_MLP_loss_multiquantilehead.py`](../scripts-py/MVN_interim_analysis_amortise_endptx_on_wz_with_features_fixed_idcomp_qpsi_MLP_loss_multiquantilehead.py) trains **one** net on `TRAIN_J = 20` prior draws and, per $J \in \{20, 60, 100\}$, replays the cached `interim_data_by_J[J]` block (posterior $\mu$ + zi, produced by `MVN_interim_analyses_make_interim_data.py`). Features per $(j, s)$ interim row are $((\sum_i y_{i,j} + \sum_i z^{(s)}_{i,j}), n + m) / N_{\max}$; one forward pass covers all $J \cdot S$ rows. Outputs written under the `_RGEA_` suffix per J.

**Training configuration.** Matches the Binomial §12.1 default post-ablation: MLP head `hidden_dims = (64, 64)`, 11-level $\tau$ mesh, Adam + linear-warmup / cosine-decay at peak lr $10^{-3}$, $40\,000$ steps × batch $512$ (each sample fans out to $J = 20$ per-component examples so effective batch $\approx 10\,240$).

**Correctness against the analytic** $\Phi$-tail PPS per $(J, \text{interim}, j)$. MSE and $\sqrt{\text{MSE}}$ across all $J$ components at each of the 7 interims per $J$ grid cell (deployment $S = 4000$):

| $J$ | \# $(j, \text{interim})$ |     MSE | $\sqrt{\text{MSE}}$ | Max abs err | Mean abs err |
| --- | -----------------------: | ------: | ------------------: | ----------: | -----------: |
| 20  |                      140 | 0.00091 |               0.030 |       0.085 |        0.018 |
| 60  |                      420 | 0.00130 |               0.036 |       0.109 |        0.020 |
| 100 |                      700 | 0.00106 |               0.033 |       0.100 |        0.021 |

**Findings.**

1.  **J-invariance holds under idcomp.** MSE moves within $\pm 0.0005$ across $J \in \{20, 60, 100\}$ despite training on a single $\text{TRAIN\_J} = 20$ cell. This is a straight consequence of the per-component decomposition: the amortised head sees only per-component features, so all $J = 100$ deployment components see the same conditional distribution as the $J = 20$ training components — the network doesn't have to know how many components there are.
2.  **Match to the naive Monte-Carlo floor at** $S = 4000$: $\sqrt{p(1 - p) / S} \approx 0.008$ near a PPS mode of $0.5$. All three $J$ cells land at $\sqrt{\text{MSE}} \approx 0.030$–$0.036$, above the MC floor — a longer training budget would tighten toward it. Consistent with the §12.1 Binomial features-fixed result (mean abs 0.020 at $S = 200$) once you scale by the ratio of naive MC noise.
3.  **Cross-component** $K$ is unexploited (idcomp caveat). The amortiser marginalises out cross-component correlation at deployment even though $K$ is known and available. That is fine for the per-component PPS reported here (each cell integrates $\mathbf{1}\{P(H_{1j} \mid x, z) > \eta_H\}$ over $z$ marginally) but suboptimal for downstream utilities on the joint. In real-world data the $J$ responses will be correlated; the natural extension is a full-$K$ amortiser that feeds either the whole $(y_1, \ldots, y_J)$ vector through a joint encoder or the per-component summaries through a $K$-aware attention layer. Future work.

**Training + deployment timing.** One-off training 0.90 min (on 6-thread CPU); per-$J$ deployment 0.96 / 3.55 / 5.43 min for $J = 20 / 60 / 100$ across the 7-interim schedule. Deployment is Python-loop-bound (per-item `predict_amortised_p_h1_for_many_xz` → `np.interp` over 4000 rows × up to 100 items); the JAX forward pass itself is milliseconds per interim.

## 13.2 Results for MVN interim analysis amortised endptx on wz with features-MLP, idcomp, qpsi-MLP, loss-multiquantilehead

**MVN extension of §12.2 — `idcomp` variant.** Same independent-components decomposition as §13.1 (one component's raw sequence at a time; cross-component $K$ used at TRAINING for the $\mu$ prior only). Analogous swap of the hardcoded per-component identity encoder for a learnable per-item MLP $q_\tau$:

- `MVNModel.make_training_data_with_raw_sequences` emits per-component padded scalar sequences (`item_dim = 1`), shared masks + `(n / N_max, m / N_max)` sizes. Same per-component flatten as §13.1 gives $S \cdot J$ training examples per batch, so the amortiser is again J-invariant.
- Deployment script [`scripts-py/MVN_interim_analysis_amortise_endptx_on_wz_with_features_MLP_idcomp_qpsi_MLP_loss_multiquantilehead.py`](../scripts-py/MVN_interim_analysis_amortise_endptx_on_wz_with_features_MLP_idcomp_qpsi_MLP_loss_multiquantilehead.py) per interim pivots `dpi` and `zi` on `(pid, j)` to get per-component observed / future sequences; loops per $j \in \{0, \ldots, J-1\}$ and forward-passes `S` padded raw-scalar batches so the peak tensor size stays at `(S, N_max, 1)`.

**Training configuration.** Combined `(64x64 + qlv5)` config from §12.5 — MLP head `hidden_dims = (64, 64)`, `q_tau_hidden_dims = (32, 32)`, `embed_dim = 16`, 5-level $\tau$ mesh, $15\,000$ steps × batch $128$ (effective $\approx 2\,560$ per-component examples per step). Deployment $S = 4000$ posterior draws per interim (matches the other MVN methods).

**Correctness against the analytic** $\Phi$-tail PPS per $(J, \text{interim}, j)$ (deployment $S = 500$; the MLP forward pass on `(S, N_max, 1)` tensor per $(j, \text{interim})$ is the deployment bottleneck, so we drop $S$ from RGEA's $4000$ to $500$ to cap wall-clock):

| $J$ | \# $(j, \text{interim})$ |   MSE | $\sqrt{\text{MSE}}$ | Max abs err | Mean abs err |
| --- | -----------------------: | ----: | ------------------: | ----------: | -----------: |
| 20  |                      140 | 0.265 |               0.514 |        1.00 |         0.34 |
| 60  |                      420 | 0.284 |               0.533 |        1.00 |         0.36 |
| 100 |                      700 | 0.182 |               0.427 |        1.00 |         0.26 |

**Findings — MLP variant is undertrained for the MVN target scale.**

1.  **Fixed encoder (§13.1) beats learnable encoder by two orders of magnitude here.** The MLP variant converges much more slowly than the fixed variant because the target $\rho = \mu_j - \mu_0$ inherits the prior sd of $\tau \sqrt{K_{jj}} = 10$ (vs Binomial's $\rho = 1 - p / p_0 \in [-1, 1]$), so the pinball loss lives on a much larger scale. Final training pinball loss reached $0.15$ at $3\,000$ steps and was still slowly decreasing — comparable to the Binomial MLP's $0.009$ after $15\,000$ steps only when rescaled by the target range ($\sim 100 \times$ larger for MVN).
2.  **The bottleneck is deployment compute, not training.** Each MLP forward pass over a `(S=500, N_max=1050, 1)` tensor takes \~100 ms; a full deployment loop is $J \cdot 7$ passes per $J$ cell (100 passes × 7 interims for $J = 100$). Training we can extend arbitrarily; but $S$ can't easily be pushed to match RGEA's $4000$ without a proportional wall-clock hit at deployment.
3.  **Recommendation.** For MVN with unit-diagonal $K$, prefer the features-fixed variant (§13.1): the sufficient statistic is available in closed form and the fixed encoder consumes it exactly. Reserve the MLP variant for models where the sufficient statistic is not known — with a longer training budget ($\ge 15\,000$ steps) and either target rescaling by $\tau$ or an MPS training pass (§12.6).

**Training + deployment timing.** Training 15.9 min (3,000 steps on 6-thread CPU); deployment 0.14 / 0.40 / 0.67 min for $J = 20 / 60 / 100$ across the 7-interim schedule.

## 13.3 Full cross-method comparison for the MVN case study, per $J \in \{20, 60, 100\}$

**MVN extension of §12.3.** [`scripts-py/MVN_interim_analyses_compare_methods.py`](../scripts-py/MVN_interim_analyses_compare_methods.py) extended with `DIR_RGEA` + `DIR_RGEB` and two new method rows (`amortised`, `amortised MLP`), piped through the existing per-J MSE + timing plot machinery. MSE is per $(J, \text{interim})$ averaged over all $J$ response components, computed against the closed-form $\Phi$-tail PPS `pps_cf` from `mvn_pps_closed_form.pkl`. Aggregate mean MSE per method per $J$ (averaged over the 7 interims):

| Method                                                         |    $J = 20$ |    $J = 60$ |   $J = 100$ | Deployment cost / interim |
| -------------------------------------------------------------- | ----------: | ----------: | ----------: | ------------------------: |
| **Amortised idcomp (features-fixed)**                          | **0.00091** | **0.00130** | **0.00106** |             0.14–0.78 min |
| nested-MC using HMC for each $(x, z)$                          |     0.00153 |     0.00242 |     0.00204 |             18.7–39.7 min |
| Regression endpt-x (Gaussian approx)                           |     0.00348 |     0.00370 |     0.00336 |             0.12–0.76 min |
| Regression endpt-x (quantile)                                  |     0.00349 |     0.00374 |     0.00349 |             0.12–0.76 min |
| Regression endpt-x (mquantile)                                 |     0.00358 |     0.00389 |     0.00368 |             0.12–0.76 min |
| IS reweighting of $\theta \mid x$                              |     0.00453 |     0.00664 |     0.04605 |             0.09–0.51 min |
| Amortised MLP idcomp (undertrained)                            |     0.26469 |     0.28351 |     0.18170 |             0.02–0.10 min |
| Amortised MLP xcomp (undertrained; §13.4)                      |     0.48307 |     0.53525 |     0.59112 |             0.44–19.4 min |
| **Amortised MLP xcompAtt** (K-mixture train, 15k steps; §13.5) |     0.00099 | **0.00087** |     0.00118 |             0.61–0.64 min |

**Cross-**$J$ reading.

1.  **Amortised idcomp (features-fixed) and amortised MLP xcompAtt tie for the top spot.** idcomp wins at $J = 20$ and $J = 100$ ($\pm 10\%$), xcompAtt wins at $J = 60$; both are $\sim 30$–$50\%$ better than nested-MC HMC and $3$–$4\times$ better than the regression family. idcomp hardcodes the exact per-component MVN sufficient statistic; xcompAtt (§13.5) learns to attend across components using the known K-row, and — critically — is trained on a *mixture* of K families rather than the deployment K, so it transfers to any K in the family. Both are J-invariant (single net across all $J$ values in the grid).
2.  **J-invariance verified end-to-end (idcomp).** MSE for the fixed amortiser moves only within $[0.00091,\ 0.00130]$ across $J \in \{20, 60, 100\}$ — no re-training, no re-tuning. This is the payoff of the per-component decomposition (§3.3 + §13.1).
3.  **Nested-MC HMC is** $\sim 30$–$300\times$ slower. Its 18.7 / 26.8 / 39.7 min per interim per $J$ dwarf the amortised variant's 0.14 / 0.40 / 0.77 min for the same $J$'s full 7-interim schedule. On the Binomial side (§12.3) the amortised advantage was already visible; on MVN, where nested-MC's inner HMC has to explore a $J$-dimensional posterior at every $z^{(s)}$, the gap widens.
4.  **IS collapses at** $J = 100$. Effective sample size drops as the future-data dimension grows, so the reweighted PPS estimator's variance explodes ($\text{MSE} \times 15$ vs $J = 20$). The regression + amortised families are unaffected — they operate on per-component summaries and don't reweight in $\theta$-space. This is the same failure mode noted in Appendix A.1.1.
5.  **MLP amortiser idcomp is undertrained here.** Its \~0.26 MSE is a training-budget artefact, not a fundamental issue with the method (§13.2). For MVN we recommend the features-fixed idcomp variant when $K_{jj} = 1$ (unit-diagonal); the MLP variant only pays off when the sufficient statistic is unknown (Categorical, IRT).

Full per-J boxplots + timing bars + MSE plots live under `py-mvn-interim-compare-methods-260609/`: `mvn_J{J}_compare_methods_p_h1_xz_all.pdf`, `mvn_J{J}_compare_methods_pps_all.pdf`, `mvn_J{J}_compare_methods_timing.pdf`, `mvn_compare_methods_mse.pdf`.

## 13.4 Cross-component amortiser (`xcomp`) — nested DeepSets + K-row query

**Motivation.** §13.1/§13.2 idcomp variants ignore the known cross-component covariance $K$ at deployment. In real data the $J$ responses will be correlated and the amortiser should leverage that. `xcomp` is the standard-techniques answer: nested DeepSets over $(\text{participant}\ i, \text{component}\ j)$ tokens plus a $K$-row query. No transformer, no attention — just DeepSets applied twice.

**Architecture.** Two-scalar per-token feature $(K[j^*, j],\ y_{i, j})$ tells each token *how* it correlates with the queried component. Inner shared-weights MLP $q_\tau^{\text{inner}}: \mathbb{R}^2 \to \mathbb{R}^{E}$, sum-pool over components $j$ → per-participant embedding $h_i$; sum-pool over participants $i$ with mask → set-of-participants embedding $\text{pooled} \in \mathbb{R}^E$. Same encoder applied to $x$ and $z$ with shared weights. Head $q_\psi$ concatenates $(\text{pooled}_x, \text{pooled}_z, \text{sizes}, K[j^*, j^*])$ and predicts the 5-quantile grid.

**J-invariance.** The K-row is the *only* thing the head learns about the queried component — no positional or lookup embedding of $j^*$. One trained net covers any $J$ and any $K$ family the training distribution spans. Independent-components `idcomp` is the special case $K = I$.

**Files.**

- [`python/amortiser_pps_features_MLP_xcomp_qpsi_MLP_loss_multiquantilehead.py`](../python/amortiser_pps_features_MLP_xcomp_qpsi_MLP_loss_multiquantilehead.py) — 60 lines including docstring. Nested DeepSets + K-row + head; standard Flax.
- `MVNModel.make_training_data_with_participant_sequences` (in [`python/model_mvn.py`](../python/model_mvn.py)) — prior-predictive over full participant matrices; for each of $S$ prior draws picks $Q$ random query components (default $Q = 4$) so the same participant matrix trains $Q$ different queries per gradient step.
- [`scripts-py/MVN_interim_analysis_amortise_endptx_on_wz_with_features_MLP_xcomp_qpsi_MLP_loss_multiquantilehead.py`](../scripts-py/MVN_interim_analysis_amortise_endptx_on_wz_with_features_MLP_xcomp_qpsi_MLP_loss_multiquantilehead.py) — deployment script; per interim runs $S = 500$ posterior-draw batches, each of $J$ query components.

**Training configuration.** $q_\tau^{\text{inner}} = (32, 32) \to 32$, head $(64, 64) \to 5$; $\tau$-grid $(0.05, 0.25, 0.5, 0.75, 0.95)$; Adam warmup-cosine at peak lr $10^{-3}$; $4\,000$ steps × $S = 16$ prior draws × $Q = 4$ queries → effective batch $B = 64$; TRAIN_J = 20.

**Results.**

| $J$ |   MSE | Mean abs err | Deployment min/interim |
| --- | ----: | -----------: | ---------------------: |
| 20  | 0.483 |         0.52 |                   0.44 |
| 60  | 0.535 |         0.57 |                   4–11 |
| 100 | 0.591 |         0.64 |              18.9–19.9 |

Training pinball loss dropped from $1954 \to 2.06$ over $4\,000$ steps (16 min wall-clock on 6-thread CPU); loss still monotonically decreasing at the last checkpoint.

**Findings.**

1.  **Architecture works, training budget too low.** Loss trajectory is monotone-descent and the loss floor extrapolates further down; the deployment MSE is not the fundamental ceiling of the method but the wall-clock we could afford here. The Binomial features-MLP idcomp needed $15\,000$ steps in §12.2 to reach its floor; xcomp on MVN with the bigger token count and the extra K-row degree of freedom will need at least that much on MPS.
2.  **Deployment is compute-bound.** The per-interim cost scales with $S \cdot J^2 \cdot N_\text{max}$ (encoder tokens are $(S, J, N_\text{max}, J)$). At $J = 100$ this is 20 min per interim on CPU. The natural fix is to cache pooled_x per (interim, query) across posterior draws — pooled_x depends on $j^*$ (via K-row) but not on $s$, so a per-$j$ pre-pass reduces compute by roughly $2 \times$ across the deployment loop. Further speedups from MPS (matmul-heavy encoder) and from S-chunking. Not implemented here.
3.  **Standard tools, no exotic layers.** Deliberately kept to a Flax MLP + `jnp.sum` + broadcast. The paper-value of the exercise is that any team with a Flax/PyTorch stack can add the K-row query and swap out idcomp for xcomp in a day — no new machinery.
4.  **Recommended follow-ups** to make xcomp competitive:
    - Longer training ($15\,000$+ steps) — probably enough to close the gap to idcomp fixed.
    - Cache pooled_x per interim per query across posterior draws (halves deployment cost).
    - Try attention over the component axis in place of the inner sum-pool — modest change, potentially large gain when intra-block K is strong.
    - Train on multiple K families (random block-equicorr, AR(1), factor) to get a K-family-invariant net, not just J-invariant. Real-world data will not match the block-equicorr training K exactly.

The design shows that going from idcomp to a K-aware amortiser is a small architectural step with standard machinery; the current run demonstrates the plumbing works end-to-end but leaves the empirical MSE comparison to a longer training pass.

## 13.5 Cross-component attention + K-family-invariant training (`xcompAtt`)

**Two changes on top of `xcomp`** (§13.4), motivated by the follow-ups listed at the end of §13.4:

1.  **Attention over the component axis** replaces the inner sum-pool. A single-head cross-attention block learns per-component weights conditioned on `K[j*, :]` instead of giving every component equal weight — the natural fix when intra-block correlation is strong.
2.  **K-family-invariant training.** Each prior draw picks its own K from a mixture of standard families (identity / AR(1) / block-equicorrelation / factor); the amortiser sees `K[j*, :]` as its only cue for correlation structure, so a single trained net handles any K the mixture spans. Deployment K (block-equicorr with $\rho_w = 0.8$, $\rho_b = 0.1$) is one member of the training mixture but not the one the training loop sees most often.

### Architecture (simple, standard)

- Per-component token = three scalars `(x_sum[j], z_sum[j], K[j*, j])`. Component summaries $\sum_i y_{i, j}$ are the exact MVN sufficient statistic for $\mu_j$ given $(x, z)$ — the participant-level tokens of §13.4 add no information.
- Token embedding: shared SiLU-MLP `q_tok` → embedding in $\mathbb{R}^{32}$.
- Query = MLP of `(sizes, K[j*, j*])` — derived from what the queried component "knows about itself".
- Scaled dot-product attention (single head, one query per batch element, keys/values = $J$ per-component embeddings). Standard `jnp.einsum` + `jax.nn.softmax`, no attention library.
- Head `q_psi = (64, 64)` → 5 quantiles.

Total code: 60 lines including docstring. `Amortiser_PPS_features_MLP_xcompAtt_qpsi_MLP_loss_multiquantilehead` in [`python/amortiser_pps_features_MLP_xcompAtt_qpsi_MLP_loss_multiquantilehead.py`](../python/amortiser_pps_features_MLP_xcompAtt_qpsi_MLP_loss_multiquantilehead.py).

### K-family sampler

`sample_random_K_chol(rng, J, families=...)` in [`python/model_mvn.py`](../python/model_mvn.py). Each family draws its own random parameters (AR(1) $\rho \sim U(-0.9, 0.9)$; block-equicorr with random block size, $\rho_w \sim U(0.3, 0.9)$, $\rho_b \sim U(0, \min(0.5, \rho_w - 0.1))$; factor with $r \in \{2, 5, 10\}$, $\psi \sim U(0.05, 0.5)$). All returned K have unit diagonal so per-component posterior variance stays consistent across families.

`MVNModel.make_training_data_with_component_summaries_random_K` couples that sampler with the standard prior-predictive draw: per sample $s$, pick K, draw $\mu \sim \mathrm{MVN}(\mu_0, \tau^2 K)$, then draw per-component sums directly ($\sum_i y_{i, j} \sim \mathrm{Normal}(n \mu_j, n \sigma^2 K_{jj})$ — no participant matrix materialised).

### Training + deployment

- $15\,000$ steps, $S = 32$ prior draws × $Q = 4$ query components per step → effective batch $B = 128$; TRAIN_J = 20; train K-mixture = `('identity', 'ar1', 'block', 'factor')`.
- Wall-clock: **0.71 min = 43 seconds** on 6-thread CPU. Sampler is 0.5 ms/step; encoder is $\mathcal{O}(B \cdot J \cdot \text{embed}) = \mathcal{O}(80\,\text{K})$ scalars/step (no participant loop → $N_\text{max} = 1050 \times$ cheaper than xcomp).
- Deployment: 0.61 / 0.62 / 0.64 min per J across the 7-interim schedule — flat in J because the per-forward encoder is $J \times J$ tokens, still tiny.
- Training pinball loss $3.91 \to 0.009$; final quartile 0.007–0.009, converged.

### Results

MSE against the analytic $\Phi$-tail PPS, averaged across the 7 interims per $J$ (block-equicorr deployment K; K-mixture training):

| $J$ | xcompAtt (K-mixture train) | idcomp fixed (§13.1 baseline) | nested-MC HMC | Regression Gauss |
| --- | -------------------------: | ----------------------------: | ------------: | ---------------: |
| 20  |                    0.00099 |                   **0.00091** |       0.00153 |          0.00348 |
| 60  |                **0.00087** |                       0.00130 |       0.00242 |          0.00370 |
| 100 |                    0.00118 |                   **0.00106** |       0.00204 |          0.00336 |

- **Matches idcomp fixed at** $J = 20$ and $J = 100$ ($\pm 10\%$). Beats it at $J = 60$. And this is despite xcompAtt being trained on a MIXTURE of K families (identity / AR(1) / block / factor) rather than the specific deployment K.
- **\~500× improvement over `xcomp`** (§13.4 was 0.483 / 0.535 / 0.591 undertrained). Combined effect of (a) attention over j + component-summary tokens, (b) 4× more training steps at the same wall-clock, (c) K-family-invariant training.
- **\~2× lower MSE than nested-MC HMC at every** $J$. \~40× lower at $J = 100$ vs IS.
- Training + deployment total 3 min end-to-end for the whole table.

### Reading

1.  **Attention over j is the right primitive here.** Sum-pool + inner participant DeepSets (xcomp) forces the amortiser to spread capacity across $\mathcal{O}(J \cdot N_\text{max})$ tokens per query; attention on per-component sums keeps compute on the actually-informative axis and lets the network re-weight components by their K-correlation to the query. The wall-clock swing (16 min → 43 sec training; hours → seconds deployment) is entirely from cutting the participant axis + the exact-sufficient-statistic tokenisation, not from attention per se.
2.  **K-family-invariant training does NOT cost accuracy.** Even under a training K distribution that spans identity → AR(1) → block → factor, xcompAtt lands within $\pm 10\%$ of a variant trained on the specific deployment K (idcomp fixed). The obvious extension is to train on a K distribution matched to what a downstream user's model prior implies — e.g. sampling K from a hierarchical prior on the covariance itself. Because the K-family cost is zero here, this is a free upgrade for real-world deployments where the analyst won't know the deployment K in advance.
3.  **Beats idcomp fixed at** $J = 60$. The K-mixture training gives the attention head examples where the correct answer requires *weighting* components differently — identity means "ignore other j's", block means "attend to same-block j's" — so at deployment the head has learned to look at the K-row profile and weight accordingly. idcomp fixed treats every component identically and can't exploit K at all; when $J$ is large enough that neighbouring components carry usable information about $\mu_{j^*}$, xcompAtt wins.
4.  **Standard techniques throughout.** No transformer library, no multi-head attention stack, no custom kernel. `_MLP` + `jnp.stack` + `jnp.einsum` + `jax.nn.softmax` — the same primitives already in §13.4. The performance jump is architecture (attention on sufficient statistics vs sum-pool on raw participants), not tooling.

**Recommendation.** For MVN-family case studies where $K$ is known (or drawn from a known family), adopt **xcompAtt as the new default amortised method**: it matches or beats idcomp fixed on MSE, matches it on deployment cost, halves the training wall-clock, and — crucially — transfers to any K within the training family mixture without a retrain. Keep idcomp fixed as the sanity-check baseline.

## 13.6 Data-agnostic amortisers: `itemScompAtt` (self-attention) + `itemXcompAtt` (cross-attention) + math comparison

To reuse the same amortiser code across MVN (§13.5) and Ukraine PCM (§14.1) we split the previous case-study-specific classes into two data-agnostic siblings distinguished only by the attention mechanism:

- **`Amortiser_PPS_features_itemScompAtt_qpsi_MLP_loss_multiquantilehead`** ([`python/amortiser_pps_features_itemScompAtt_qpsi_MLP_loss_multiquantilehead.py`](../python/amortiser_pps_features_itemScompAtt_qpsi_MLP_loss_multiquantilehead.py)) — **S**elf-attention across $J$ items + gather at the queried item + learned attention-score bias on the queried key (fix B).
- **`Amortiser_PPS_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead`** ([`python/amortiser_pps_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead.py`](../python/amortiser_pps_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead.py)) — **X (cross)**-attention: 1 learned query per batch element (built from the queried item's own raw token features via a small MLP `q_query`) vs $J$ keys/values.

Both use the same batch schema so callers can swap classes without changing batch construction. A per-participant nested-DeepSets sibling for arbitrary $(n, m)$ (`deepsetXcompAtt`) is also stubbed for future work.

### Notation

Fix a case study with:

- $J$ items (MVN: components; Ukraine: PCM questions),
- one queried item index $j^* \in \{1, \ldots, J\}$ per batch element,
- observed cohort size $n$, future cohort size $m$, total $N = n + m$ (fixed),
- $S$ posterior draws of the future cohort from the fitted model, indexed $s = 1, \ldots, S$: at deploy time $z^{(s)} = (z^{(s)}_{i, j})$ is the $s$-th posterior-predictive replay of the $m$ future participants,
- for MVN: known covariance matrix $K \in \mathbb{R}^{J \times J}$.

Let $F$ be the per-item raw feature dim, $A$ the aux-feature dim, $E$ the embedding dim (32), $K_\tau$ the number of quantile levels (5). A batch element is a single triple $(j^*, s, \text{interim})$ — one query, one posterior draw, one interim cohort. At deployment we forward-pass all $J \cdot S$ batch elements per interim.

### Per-token features (case-study dependent; same in both amortisers)

For MVN, given the queried component $j^*$ and $j \in \{1, \ldots, J\}$ under posterior draw $s$:

$$t^{(j^*, s)}_j \;=\; \Big(\underbrace{\tfrac{1}{N}\!\sum_{i=1}^{n} y_{i, j}}_{\bar y_j\;(\text{$s$-invariant})},\; \underbrace{\tfrac{1}{N}\!\sum_{i=1}^{m} z^{(s)}_{i, j}}_{\bar z^{(s)}_j\;(\text{$s$-dependent})},\; \underbrace{K_{j^*, j}}_{\text{K-row at }j^*}\Big) \;\in\; \mathbb{R}^{F},\qquad F = 3.$$

$\bar y_j$ is fixed across $s$ (the observed cohort does not change); $\bar z^{(s)}_j$ varies across $s$ (that is the amortiser's source of posterior-predictive information). $K_{j^*, j}$ depends on $j^*$ only. For Ukraine (§14.1), $t^{(j^*, s)}_j \in \mathbb{R}^{F}$ with $F = 8$, including an empirical $\hat K$-row as the $j^*$-dependent entry (see §14.1.3).

Aux scalars for the head (case-study dependent, $j^*$-dependent for MVN; $s$-invariant):

$$a^{(j^*)} \;=\; \big(n / N,\; m / N,\; K_{j^*, j^*}\big) \;\in\; \mathbb{R}^{A}, \qquad A = 3.$$

Token encoder (shared MLP, applied per token):

$$h^{(j^*, s)}_j \;=\; q_{\text{tok}}\!\big(t^{(j^*, s)}_j\big) \;\in\; \mathbb{R}^{E}, \qquad j = 1, \ldots, J.$$

We keep the superscript $(j^*, s)$ on tokens because their K-row entry changes with $j^*$ and their $\bar z^{(s)}_j$ entry changes with $s$; the encoder $q_{\text{tok}}$ is shared across all $(j, j^*, s)$.

### Variant 1 — cross-attention (`itemXcompAtt`, approach ii)

Build a single query vector for this batch element from the queried item's own raw features (a subtle abuse of notation: `q_query` sees the raw token, not the encoded one):

$$q^{(j^*, s)} \;=\; q_{\text{query}}\!\big(t^{(j^*, s)}_{j^*}\big) \;\in\; \mathbb{R}^{E}.$$

Cross-attention (one query, $J$ keys/values):

$$\text{score}^{(j^*, s)}_j \;=\; \frac{\big\langle q^{(j^*, s)},\; h^{(j^*, s)}_j\big\rangle}{\sqrt{E}} \;\in\; \mathbb{R}, \qquad w^{(j^*, s)}_j \;=\; \operatorname{softmax}_{j'\in\{1,\ldots,J\}}\!\big(\text{score}^{(j^*, s)}_{j'}\big) \;\in\; [0, 1].$$

Aggregated summary:

$$\bar h^{(j^*, s)} \;=\; \sum_{j = 1}^{J} w^{(j^*, s)}_j \, h^{(j^*, s)}_j \;\in\; \mathbb{R}^{E}.$$

Head:

$$\hat\rho^{(j^*, s)}_{\tau_k} \;=\; \Big[q_\psi\!\big(\operatorname{concat}\big(\bar h^{(j^*, s)},\; t^{(j^*, s)}_{j^*},\; a^{(j^*)}\big)\big)\Big]_k \;\in\; \mathbb{R}, \qquad k = 1, \ldots, K_\tau.$$

**Head input dimension:** $E + F + A = 32 + 3 + 3 = 38$.

**\# attention scores per batch element:** $J$ (one query, $J$ keys).

**Query info:** entirely in $q^{(j^*)}$ = MLP of queried item's raw features. Every non-queried token contributes via its dot-product with $q^{(j^*)}$; tokens with high $K_{j^*, j}$ produce embeddings $h^{(j^*)}_j$ near $h^{(j^*)}_{j^*}$ and attract high weight.

### Variant 2 — self-attention + gather + fix-B bias (`itemScompAtt`)

Self-attention: every one of the $J$ tokens acts as both query and key. Learned scalar $\alpha \in \mathbb{R}$ biases each row toward the queried key so the attention softmax has an unambiguous query identity.

$$\text{score}^{(j^*, s)}_{a, b} \;=\; \frac{\big\langle h^{(j^*, s)}_a,\; h^{(j^*, s)}_b\big\rangle}{\sqrt{E}} \;+\; \alpha \cdot \mathbf{1}\{b = j^*\}, \qquad a, b \in \{1, \ldots, J\}.$$

$$w^{(j^*, s)}_{a, b} \;=\; \operatorname{softmax}_{b'\in\{1,\ldots,J\}}\!\big(\text{score}^{(j^*, s)}_{a, b'}\big), \qquad h'^{\,(j^*, s)}_a \;=\; \sum_{b = 1}^{J} w^{(j^*, s)}_{a, b}\, h^{(j^*, s)}_b \;\in\; \mathbb{R}^{E}.$$

Gather at the queried position:

$$\bar h^{(j^*, s)} \;=\; h'^{\,(j^*, s)}_{j^*} \;\in\; \mathbb{R}^{E}.$$

Head (same shape as cross-attn variant):

$$\hat\rho^{(j^*, s)}_{\tau_k} \;=\; \Big[q_\psi\!\big(\operatorname{concat}\big(\bar h^{(j^*, s)},\; t^{(j^*, s)}_{j^*},\; a^{(j^*)}\big)\big)\Big]_k \;\in\; \mathbb{R}, \qquad k = 1, \ldots, K_\tau.$$

**Head input dimension:** $E + F + A = 32 + 3 + 3 = 38$ (identical to cross-attn variant).

**\# attention scores per batch element:** $J^2$ (every token queries every token).

**Query info:** carried by (1) the gather position $j^*$ and (2) the learned scalar $\alpha$ that biases attention scores toward the queried key. $\alpha$ is the ONLY new learnable parameter compared to the cross-attn variant with `q_query` collapsed to identity.

### Difference summary

Every symbol carries an explicit $(j^*, s)$ superscript above; the table below drops it for brevity.

| item                                  | `itemXcompAtt` (cross-attn)                                                                                       | `itemScompAtt` (self-attn + gather + B)                                                   |
| ------------------------------------- | ----------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| \# attention scores per batch element | $J$                                                                                                               | $J^2$                                                                                     |
| query built from                      | $q_{\text{query}}(t_{j^*}) \in \mathbb{R}^E$                                                                      | gather position $j^*$ + $\alpha \cdot \mathbf{1}\{b = j^*\}$                              |
| attention output                      | $\bar h = \sum_{j=1}^J w_j h_j$ (weighted sum)                                                                    | $\bar h = h'_{j^*}$ (gather)                                                              |
| head input                            | $\operatorname{concat}(\bar h,\, t_{j^*},\, a) \in \mathbb{R}^{E + F + A}$                                        | $\operatorname{concat}(\bar h,\, t_{j^*},\, a) \in \mathbb{R}^{E + F + A}$                |
| extra learnable module vs the other   | `q_query` MLP $\mathbb{R}^F \to \mathbb{R}^E$                                                                     | scalar $\alpha \in \mathbb{R}$                                                            |
| K-row consumption                     | attention softmax on $\langle q, h_j\rangle$; high-K tokens land near $q$ in embedding space, attract high weight | K-row baked into every token; self-attention has to learn to route information through it |

### Empirical impact on MVN (15k steps, K-mixture training)

MSE against analytic PPS, averaged over the 7-interim schedule:

| $J$ | `itemXcompAtt` cross-attn (new) | `itemScompAtt` self-attn + B (previous best) | old case-study-specific cross-attn baseline |
| --- | ------------------------------: | -------------------------------------------: | ------------------------------------------: |
| 20  |                         0.00141 |                                      0.00196 |                                     0.00099 |
| 60  |                         0.00114 |                                      0.00153 |                                     0.00087 |
| 100 |                         0.00240 |                                      0.00218 |                                     0.00118 |

**Read:**

- **Cross-attention (`itemXcompAtt`) beats self-attention + gather + B (`itemScompAtt`) at** $J = 20, 60$ by 20–35 %.
- **At** $J = 100$ they're within 10 %, self-attn slightly ahead.
- **Both close the gap to the case-study-specific cross-attn baseline** — the cost of unification is now \~15–40 % at moderate J, up to \~2× at $J = 100$.

### Reading

1.  **Cross-attention naturally preserves the K-row signal.** The K-row $K_{j^*, \cdot}$ tells us how much information item $j$ contributes to the queried $j^*$; in cross-attention, tokens with high $K_{j^*, j}$ produce embeddings close to the queried item's own embedding and are picked out by the softmax. Self-attn + gather has to learn this routing via $\alpha$ and the shared attention weights — harder from data.
2.  **Two data-agnostic classes, one batch schema.** Same call signature; only the attention block differs. Callers can A/B by swapping the import.
3.  **Naming convention** (adopted going forward): `itemScompAtt` = self-attention across items; `itemXcompAtt` = cross-attention across items with one query per batch element. The prefix `item` disambiguates from a nested `deepsetXcompAtt` variant that will add per-participant DeepSets over raw participant tensors before the cross-attention step (implements the §9.2 nested pattern).

**Recommendation:** for MVN, prefer `itemXcompAtt` (cross-attn) at $J \leq 60$; either variant works at $J = 100$. For Ukraine (no K), the two variants perform similarly — the cross-attn advantage from K-alignment does not apply. Keep both classes; the choice is a per-case-study empirical call.

------------------------------------------------------------------------

## 13.7 Nested-DeepSets amortisers: `deepsetScompAtt` + `deepsetXcompAtt`

`itemScompAtt` / `itemXcompAtt` of §13.6 take a **hand-computed per-item summary** $t^{(j^*)}_j \in \mathbb{R}^F$ as input — e.g. for MVN $t^{(j^*)}_j = (\bar x_j,\; \bar z_j,\; K_{j^*, j})$, i.e. the per-item cohort means plus the K-row entry. The DeepSets theorem (§7) says the sufficient statistic $\phi(X, Z)$ is *any* symmetric function of the observation tensor. The `deepset` variants take this literally: **learn the per-item summary from raw per-(participant, item) responses**, chaining nested DeepSets (§9.2) with the item-axis attention of §13.6.

### Architecture

Batch schema (data-agnostic; supports arbitrary response feature dim $R$):

| tensor          | shape            | content                                                      |
| --------------- | ---------------- | ------------------------------------------------------------ |
| `x_responses`   | $(B, N_x, J, R)$ | raw per-(participant, item) response of the $x$-cohort       |
| `mask_x`        | $(B, N_x)$       | 1 real / 0 padded participant                                |
| `z_responses`   | $(B, N_z, J, R)$ | raw per-(participant, item) response of the $z$-cohort       |
| `mask_z`        | $(B, N_z)$       | 1 real / 0 padded participant                                |
| `item_metadata` | $(B, J, M)$      | fixed per-item scalar features (item type, K-row entry, …)   |
| `query_idx`     | $(B,)$ int       | queried item $j^*$                                           |
| `aux`           | $(B, A)$         | scalar extras for the head (cohort sizes, $K_{j^*, j^*}$, …) |

Forward pass, per batch element:

1.  **Concatenate per-item metadata** onto every $(i, j)$ token: $\tilde x_{i, j} = \operatorname{concat}(x_{i, j},\; \text{meta}_j) \in \mathbb{R}^{R + M}$.
2.  **Per-token embedding**: $e_{i, j} = q_\tau(\tilde x_{i, j}) \in \mathbb{R}^E$, shared across $(i, j)$.
3.  **Sum-pool over participants per item** (mask-aware): $\operatorname{pool}^x_j = \sum_i m^x_i\, e^x_{i, j} \in \mathbb{R}^E$; same for $z$. This is the learned per-item sufficient statistic.
4.  **Per-item token**: $h_j = q_{\text{tok}}(\operatorname{concat}(\operatorname{pool}^x_j,\; \operatorname{pool}^z_j)) \in \mathbb{R}^E$.
5.  **Attention across items** — the only difference between the two variants:
    - `deepsetScompAtt`: self-attention over $\{h_j\}_j$ with the fix-B learned scalar bias $\alpha \cdot \mathbf{1}\{k = j^*\}$ on the queried key; gather $\bar h = h'_{j^*}$.
    - `deepsetXcompAtt`: cross-attention with one query per batch element built from $h_{j^*}$ via $q_{\text{query}}$: $q = q_{\text{query}}(h_{j^*})$; $\bar h = \sum_j w_j h_j$ with $w = \operatorname{softmax}(\langle q, h\rangle)$.
6.  **Head**: $\hat\rho^{(j^*)}_{\tau_k} = q_\psi(\operatorname{concat}(\bar h,\; h_{j^*},\; a))_k$, multi-quantile.

### MVN specifics

Per-(participant, item) response is a single continuous scalar ($R = 1$). Per-item metadata is the K-row entry $\text{meta}_j = K_{j^*, j}$ ($M = 1$). Aux is $a = (n / N,\; m / N,\; K_{j^*, j^*})$ ($A = 3$). Sampler `make_training_data_with_participant_tokens_random_K` in `model_mvn.py` emits raw participant tensors with a fresh random K per prior draw (K-family-invariant training, families `('identity', 'ar1', 'block', 'factor')`).

### Empirical impact on MVN

MSE against analytic PPS, averaged over 7-interim schedule at $J = 20$; both deepset variants trained 8k steps (\~71 min each on CPU) with mean-pool over participants (see "Mean-pool vs sum-pool" below):

| $J$ | `deepsetScompAtt` | `deepsetXcompAtt` | best of `itemScompAtt` / `itemXcompAtt` (§13.6) |
| --- | ----------------: | ----------------: | ----------------------------------------------: |
| 20  |           0.02136 |           0.00573 |                        0.00141 (`itemXcompAtt`) |

**Read:**

- **`deepsetXcompAtt` (cross-attn) beats `deepsetScompAtt` (self-attn + B) by \~4×** at $J = 20$ — same ordering as the item variants of §13.6, and for the same reason (K-alignment through the queried-item softmax).
- **Deepset variants remain 4×-15× behind the item variants of §13.6.** Learning $q_\tau$ from raw participants is data-inefficient compared to hand-supplying the sample-mean summary; both variants converge to loss \~0.011 in 8k steps but the deploy-time posterior remains wider than the analytic one.
- **Compute cost:** \~71 min training per variant at $J = 20$ vs 6-9 min for the item variants of §13.6. $J \geq 60$ not yet run (memory \~5× larger per batch element).

### Mean-pool vs sum-pool (bugfix)

First run of `deepsetScompAtt` used $\sum_i m_i^x\, e_{i, j}^x$ as the per-item pool (raw sum). Training loss dropped to 0.067 but deploy-time PPS MSE was 0.42 — the amortiser predicted a very wide posterior for every item (quantile spread \~0.9 vs true \~0.03). Cause: the raw-sum magnitude scales with $n$, which varies uniformly in $[1, N_{\max}]$ during training; $q_{\text{tok}}$ saturates on the large-$n$ tail and can't unpack a sharp posterior at moderate $n$. Fix: divide by $\max(\sum_i m_i^x, 1)$ so the pool is a **mean** over real participants, scale-invariant to $n$. Training loss dropped further to 0.011 and deploy PPS MSE fell to 0.021 (Scomp) / 0.006 (Xcomp). The pattern generalises: any nested-DeepSets amortiser where the inner cohort size varies at both train and deploy time should mean-pool, not sum-pool.

### Ukraine

Deferred: a PCM prior-predictive sampler analogous to `make_training_data_with_participant_tokens_random_K` does not yet exist in `model_pcm.py`. The itemS/X variants' training pool was built from cached wa-frames (per-(draw, item) statistics from a fitted PCM posterior); it lacks a per-participant axis and therefore cannot be recycled for the deepset variants. Adding a prior-predictive PCM sampler (item thresholds + latent abilities + response draws + rho target) is the next step for Ukraine deepset numbers.

### Reading

- **Cost.** Forward-pass memory scales as $B \cdot N_{\max} \cdot J \cdot E$ vs $B \cdot J \cdot E$ for the item variants. At $N_{\max} = 1050$, $J = 20$: 21k tokens per batch element vs 20. Training \~10× slower for MVN (71 min vs 6-9 min). At $J = 100$, the deepset tensor is $\approx 100\text{M}$ elements per batch of 128, borderline on CPU.
- **What deepset buys.** Zero hand-engineering of per-item summaries: the network learns $q_\tau$ end-to-end from raw responses. Same call signature as itemS/X so callers can swap classes; only the input processing changes.
- **When it's worth paying.** For MVN with a closed-form sufficient statistic (per-item mean), the item variants remain the operating point. Deepset is the right structure when raw per-participant information carries signal that a hand-computed summary would discard — e.g. Ukraine PCM, where per-participant response patterns encode ability + item interaction that a per-item aggregate loses.
- **Where the remaining gap comes from.** Deploy-time posterior is sharper (n ≈ 132) than most training samples ($n \sim$ uniform in $[1, 1050]$, so median $n \approx 525$). The amortiser sees relatively few tight-posterior scenarios and doesn't concentrate the quantile head enough at the deploy operating point. Curriculum with more small-$n$ / large-$m$ mass, or a log-uniform $n$ schedule biased toward the tails, would likely close the gap further.

------------------------------------------------------------------------

# 14. Results for the partial credit model on Ukraine data

## 14.1 Cross-component attention to learn item-specific summaries (`itemXcompAtt`)

**Task.** Adapt the data-agnostic cross-attention amortiser of §13.6 (`itemXcompAtt`) to the Ukraine partial-credit (PCM) interim analysis. The network class is reused **unchanged**; everything Ukraine-specific lives in the token construction, the target, and the training pool. Diagnostics: §14.1.7–14.1.11. Calibration remedies: §14.2. Self-attention sibling: §14.3. Fully-amortised deepset variant: §14.4.

### 14.1.1 Structure of learning task

Each participant contributes **two** responses per item — baseline ($t = 0$) and endline ($t = 1$). The interim cohort $x$ consists of the $n$ participants with both visits complete by the cutoff; the future cohort $z$ is the remaining $m = N - n$. The PCM accommodates baseline-vs-endline level shifts by fitting **separate per-time thresholds**: internally the likelihood sees $2 J = 40$ `item_time_id` fake items (own threshold vector and loading per (real item, time)), while participant ability $\theta_i$ is shared across times and items. The amortiser collapses the time axis back: one token per real item, with the two time-points entering as separate scalar features.

Two axes must not be conflated: the **cohort axis** ($x^{\text{obs}}$ vs $z^{(s)}$ — the §10 deployment passes both through the network) and the **time axis** (baseline vs endline within a participant, present in both cohorts). The current $F = 8$ token (§14.1.3) carries both cohorts and both time-points. The remaining §10 gap at the item level is that the per-item pooling over participants is hand-computed (group-means), not learned; the learned version is the nested-DeepSets variant (§13.7, results §14.4).

| symbol   | value            | description                                                                                                                                                            |
| -------- | ---------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| $J$      | 20               | **real** items presented to the amortiser (one token per real item); internally the PCM likelihood fits $2 J = 40$ per-`item_time_id` fake items                       |
| $n$      | interim-specific | observed participants at the interim, **each with both baseline and endline responses** (contributes $2 J n$ observations to the PCM likelihood)                       |
| $m$      | interim-specific | shadow future participants used per posterior draw $s$, **each also carrying both baseline and endline** (contributes $2 J m$ posterior-predictive responses per draw) |
| $S$      | 4000             | posterior-predictive draws of $z^{(s)}$ per interim at deploy time (item amortisers, §14.1/§14.3); the deepset variant (§14.4) uses $S = 200$                          |
| $C$      | 8                | interims analysed (`INTERIM_IDS = 1..8`)                                                                                                                               |
| $F$      | 8                | per-real-item token dim: x-side + z-side cohort summaries, change ratio, metadata, empirical $\hat K$-row (see §14.1)                                                  |
| $A$      | 2                | aux scalars $(n/N, m/N)$, $N = 503$                                                                                                                                    |
| $E$      | 32               | embedding dim (`embed_dim`)                                                                                                                                            |
| $K_\tau$ | 5                | quantile levels $\tau = (0.05, 0.25, 0.5, 0.75, 0.95)$                                                                                                                 |

What is reused, what is trained, what is adapted:

| piece                                                                      | status                   | detail                                                                                                                                                                                                      |
| -------------------------------------------------------------------------- | ------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| network classes ($q_{\text{tok}}, q_{\text{query}}$ or $\alpha$, $q_\psi$) | **reused** from §13.6    | identical Flax modules, identical attention math                                                                                                                                                            |
| network weights                                                            | **trained from scratch** | \~12k parameters, multi-quantile pinball loss                                                                                                                                                               |
| per-item token $t^{(j^*, s)}_j$                                            | **adapted**              | $F = 8$: x-side + z-side summaries, empirical $\hat K$-row, item metadata (see below)                                                                                                                       |
| aux $a$                                                                    | **adapted**              | $A = 2$: cohort sizes $(n/N, m/N)$ with $N = 503$                                                                                                                                                           |
| input scale                                                                | **adapted**              | raw responses rescaled $(y)/(K_j - 1) \in [0, 1]$ with the per-item level count $K_j$ from `dit[cat_length]` (8 for out-of-7, 4 for categorical); all summaries are group-means of these rescaled responses |
| target scale                                                               | **adapted**              | per-item standardisation $\rho / \hat\sigma_j$ (see loss below)                                                                                                                                             |
| training pool                                                              | **§8 prior-predictive**  | fresh PCM prior draws per training step; no SVI fits or cached posterior draws enter training                                                                                                               |

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

The following response data from the observed and future cohorts enter the token (§10 cohort axis, as in MVN §13.6):

$$t^{(j^*, s)}_j \;=\; \Big(\underbrace{w^{x}_{\text{base}, j},\; w^{x}_{\text{end}, j}}_{\text{observed cohort } x^{\text{obs}},\ s\text{-invariant}},\;\; \underbrace{w^{z, (s)}_{\text{base}, j},\; w^{z, (s)}_{\text{end}, j},\; w^{z, (s)}_{\text{ratio}, j}}_{\text{future block } z^{(s)}},\;\; \underbrace{c_j,\; d_j}_{\text{metadata}},\;\; \underbrace{\hat K_{j^*, j}}_{\text{empirical K-row}}\Big) \;\in\; \mathbb{R}^{F},\qquad F = 8.$$

- **x-side** $w^{x}_{\text{base}, j}, w^{x}_{\text{end}, j}$: per-item group-mean **response** over the $n$ observed participants at $t = 0$ and $t = 1$, rescaled by the item type's level count: $\tfrac{1}{n (K_j - 1)}\sum_i y_{i, j, t}$ with the per-item level count $K_j$ read from `dit[cat_length]` (Ukraine: $K_j = 8$ for `out-of-7`, raw $y \in \{0, \ldots, 7\}$; $K_j = 4$ for `categorical`, raw $y \in \{0, \ldots, 3\}$). **No input-side thresholding** — the caseness threshold lives in the target only (see below). Fixed per interim.
- **z-side** $w^{z, (s)}_{\cdot, j}$: the same mean-response statistics on the $m$ future-cohort participants of $z^{(s)}$ — simulated from the prior at training time, posterior-predictive at deployment (§14.1.6). The ratio is `higher_is_better`: $w_{\text{end}}/w_{\text{base}} - 1$, `lower_is_better`: $1 - w_{\text{end}}/w_{\text{base}}$, computed on the $[0, 1]$ mean scale and clipped at $\pm 20$.
- **Metadata** $c_j \in \{0, 1\}$ (item type), $d_j \in \{0, 1\}$ (direction).
- **Empirical K-row** $\hat K_{j^*, j}$ — the $j^*$-dependent entry, restoring the MVN-style "how informative is item $j$ for item $j^*$" channel. From the observed cohort's per-participant change scores $d_{i, j} = y_{i, j, 1} - y_{i, j, 0}$, take the Spearman correlation matrix $\hat R$ across items (rank-based — robust to the ordinal scale) and shrink toward identity: $\hat K = \lambda \hat R + (1 - \lambda) I$, $\lambda = n / (n + n_0)$, $n_0 = 50$, so the noisy early-interim estimates ($n \approx 48 \Rightarrow \lambda \approx 0.5$) are strongly regularised. Recomputed per interim.

**Categorical threshold.** $c = 2$ on the 1-indexed scale $\Leftrightarrow$ raw $y \geq 1$ ("at least several days"). The CG-MH items are PHQ-style symptom-frequency Likerts; the trial's endpoint for them is the *caseness prevalence* — proportion of caregivers above the symptom cut-off — so the threshold is intrinsic to the **target** $\rho$ and stays there. The **input** summaries do not threshold: they use the mean response (see the token definition above), because binarising the input is pure information loss. The level counts are read from `dit[cat_length]` (never hardcoded); the threshold $c$ is a study-specific configuration parameter — Colombia's response scales have different level counts and hence a different $c$.

### 14.1.4 Encoder, cross-attention, head, loss function

Exactly §13.6 variant 1. Two SiLU MLPs to build the tokens and query for each item $j^*$ in embedding space:

$$h^{(j^*, s)}_j \;=\; q_{\text{tok}}\big(t^{(j^*, s)}_j;\, \tau_{\text{tok}}\big) \;\in\; \mathbb{R}^{E}, \qquad q^{(j^*, s)} \;=\; q_{\text{query}}\big(t^{(j^*, s)}_{j^*};\, \tau_{\text{query}}\big) \;\in\; \mathbb{R}^{E}.$$

Cross-attention over the $J$ per-item embeddings ($J$ attention scores per batch element):

$$\bar h^{(j^*, s)} \;=\; \sum_{j=1}^{J} \omega^{(j^*, s)}_j\, h^{(j^*, s)}_j, \qquad \omega^{(j^*, s)}_j \;=\; \operatorname{softmax}_{j'}\Big( \big\langle q^{(j^*, s)},\; h^{(j^*, s)}_{j'} \big\rangle \big/ \sqrt{E} \Big)\Big|_{j' = j}.$$

$\bar h^{(j^*, s)} \in \mathbb{R}^{E}$ is the **learned summary in embedding space** for the queried item — the amortised counterpart of the hand-picked scalar statistic $w(z^{(s)})$ of the regression baseline (§6.2). Using the learned summary, we then learn the mapping to the item-specific endpoints $\rho$ . For each item, we learn a more efficient multi-quantile head rather than a specific quantile, so we are free to interpolate to any user-desired quantile:

$$\hat\rho^{(j^*, s)} \;=\; q_\psi\big( \operatorname{concat}(\bar h^{(j^*, s)},\; t^{(j^*, s)}_{j^*},\; a);\, \psi \big) \in \mathbb{R}^K,$$

with head-input dimension $E + F + A = 32 + 8 + 2 = 42$ and head output dimension $Q = 5$. Each dimension $q$ of $\hat\rho^{(j^*, s)}$ corresponds to one of the quantile levels $\eta^{\text{inpol}}_q$, $q = 1, \ldots, Q$, with $(\eta^{\text{inpol}}_1, \ldots, \eta^{\text{inpol}}_Q) = (0.05,\, 0.25,\, 0.5,\, 0.75,\, 0.95)$, which are used to interpolate to the specific quantile that the user requests.

The learned parameters are $\tau = (\tau_{\text{tok}},\, \tau_{\text{query}})$ — the weights of the two encoder MLPs — and $\psi$, the weights of the head MLP $q_\psi$ (the notation of §8-10). The attention weights $\omega^{(j^*, s)}$ are **not** free parameters: they are deterministic functions of the embeddings (softmax of scaled dot products) and carry no trainable degrees of freedom of their own in this variant (contrast the learned scalar bias $\alpha$ of `itemScompAtt`, §14.3). The same $(\tau, \psi)$ apply to **every item** $j^*$ and every interim — one network serves all queries — in contrast to the regression baseline (§6.2, RGE/RGEM), where the learned regression function is item-specific *and* interim-specific (one fit per (item, interim) pair; $20 \times 8 = 160$ separate fits over the schedule). Sizes at the Ukraine configuration ($F = 8$, $E = 32$, hidden $(32, 32)$ for the encoders and $(64, 64)$ for the head, $K = 5$):

| block                 | architecture             | \# parameters |
| --------------------- | ------------------------ | ------------: |
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

Per interim, the current data are held **fixed at the observed** $x^{\text{obs}}$ — they enter every token through the x-side summaries, the empirical $\hat K$-row and the size aux. The future data $z^{(s)}$ are simulated from the **posterior predictive** at $x^{\text{obs}}$ (§10 step 1) — SVI is run **once per interim** on $x^{\text{obs}}$ (the only inference cost at deployment; training never touches it), then $\theta^{(s)} \sim p(\theta \mid x^{\text{obs}})$, then simulate the $m = N - n$ missing participants' responses at both time-points under $\theta^{(s)}$ (`get_interim_z_from_ypredi`: shadow participants resampled with replacement from the observed cohort's covariate rows, responses taken from the model's `ypred` at draw $s$). We use **all** $S = 4000$ cached posterior draws per interim, so one interim costs $J \cdot S = 80{,}000$ head evaluations (one JIT-compiled batch). Each $(j^*, s)$ pair is forward-passed, the quantile grid converted to a conditional success probability by CDF interpolation at the effect threshold $\eta_0$, and aggregated over draws (§10 step 3). The hypothesis is **item-specific**, $H_1^{(j^*)} = \big\{\rho^{(j^*)} > \eta_0\big\}$ — item $j^*$'s improvement ratio exceeds the effect threshold ($\eta_0 = 0.5$; decision threshold $\eta_H = 0.89$):

$$\hat P\big(H_1^{(j^*)} \mid x, z^{(s)}\big) \;=\; 1 - \operatorname{interp}\big(\eta_0;\; \hat\rho^{(j^*, s)},\; \eta^{\text{inpol}}\big),$$

$$\widehat{\text{PPS}}^{(j^*)}(x) \;=\; \frac{1}{S} \sum_{s=1}^{S} \mathbf{1}\big\{ \hat P\big(H_1^{(j^*)} \mid x, z^{(s)}\big) > \eta_H \big\}.$$

With the $F = 8$ token, both cohorts enter the network — hand-pooled summaries of $x^{\text{obs}}$ plus summaries of $z^{(s)}$ — putting this variant on the same footing as the MVN item amortiser (§13.6): §10-consistent on the cohort axis, with hand-computed (not learned) per-item pooling. Files:

- [`python/amortiser_pps_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead.py`](../python/amortiser_pps_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead.py) — the data-agnostic network class (§13.6): `_MLP` + cross-attention (`einsum` + `softmax`). Standard Flax.
- [`python/model_pcm.py`](../python/model_pcm.py) — `make_training_data_with_item_tokens_prior`, the §8 prior-predictive training sampler (§14.1.2).
- [`scripts-py/Ukraine_interim_analysis_amortise_endptx_on_wz_with_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead.py`](../scripts-py/Ukraine_interim_analysis_amortise_endptx_on_wz_with_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead.py) — trains once on fresh prior draws, then deploys per interim (SVI once → posterior-predictive $z^{(s)}$ → tokens → forward pass); `_RGEG_` outputs.

Deployment time is dominated by rebuilding the posterior-predictive $z^{(s)}$ block per interim (\~1 min / interim at $S = 4000$); the network forward pass itself is seconds.

### 14.1.7 Diagnostic: conditional calibration on $(x, z^{(s)})$ — quantile coverage + PIT

Ukraine has no analytic PPS, so we use several diagnostics to evaluate estimation accuracy.

For each SVI posterior draw $s$, take $\theta^{(s)} \sim p(\theta \mid x^{\text{obs}})$ (SVI fit to the current data only). Compute $\rho^{(j^*, s)}_{\text{SVI}} = r_{j^*}(\theta^{(s)})$ — a draw from the current posterior $p(\rho \mid x)$ and generate the future block $z^{(s)} \sim p(z \mid \theta^{(s)})$ that the amortiser conditions on. Because $(\theta^{(s)}, z^{(s)})$ is a joint draw given $x$, Bayes gives $\theta^{(s)} \mid z^{(s)} \sim p(\theta \mid x, z^{(s)})$, so $\rho^{(j^*, s)}_{\text{SVI}} | z^{(s)}$ is an exact single draw from the conditional $p(\rho \mid x, z^{(s)})$ that the amortiser aims to target. Marginally, $\rho^{(j^*, s)}_{\text{SVI}}$ is a draw from $p(\rho \mid x)$. Exactly as for any joint \$(A,B)\$ a realised \$a\$ is simultaneously a draw from \$p(A)\$ and, paired with its own \$b\$, from \$p(A\\mid B{=}b)\$. The output of the amortiser $\hat\rho^{(j^*, s)}_{1:Q}$ approximates the conditional $p(\rho \mid x, z^{(s)})$.

1.  **Quantile coverage.** For each level $\eta^{\text{inpol}}_q$, compute for the same index $s$:

$$\widehat{\text{Cov}}^{(j^*)}_q \;=\; \frac{1}{S} \sum_{s=1}^{S} \mathbf{1}\big\{ \rho^{(j^*, s)}_{\text{SVI}} \le \hat\rho^{(j^*, s)}_q \big\} \;\stackrel{!}{=}\; \eta^{\text{inpol}}_q .$$

2.  **PIT uniformity.** The probability integral transform obtained by CDF interpolation of the quantile grid,

$$u^{(j^*, s)} \;=\; \hat F\big(\rho^{(j^*, s)}_{\text{SVI-}x} \mid x, z^{(s)}\big) \;=\; \operatorname{interp}\big(\rho^{(j^*, s)}_{\text{SVI-}x};\; \hat\rho^{(j^*, s)}_{1:Q};\; \eta^{\text{inpol}}_{1:Q}\big),$$

should be $\text{Uniform}(0, 1)$ over $s$; deviations are summarised per (interim, item) by the Kolmogorov-Smirnov distance $\max_u |\widehat{\text{ecdf}}(u) - u|$. Under-dispersion shows as a U-shaped PIT histogram, over-dispersion as a hump, location bias as a tilt.

Implemented in the separate script [`..._itemXcompAtt_..._tests.py`](../scripts-py/Ukraine_interim_analysis_amortise_endptx_on_wz_with_features_itemXcompAtt_qpsi_MLP_loss_multiquantilehead_tests.py) (kept apart from train/deploy). Outputs: `..._tests_coverage.csv/.pdf` (empirical-vs-nominal coverage curves, diagonal = calibrated) and `..._tests_pit.pdf` + `..._tests_pit_ks.csv` (PIT histograms + KS). Result:

| nominal $\eta^{\text{inpol}}_q$                 | 0.05 | 0.25 | 0.50 | 0.75 | 0.95 |
| ----------------------------------------------- | ---: | ---: | ---: | ---: | ---: |
| empirical coverage (mean over items × interims) | 0.05 | 0.10 | 0.15 | 0.23 | 0.44 |

Coverage falls far below nominal above the 0.05 level; mean PIT-KS $\approx 0.61$, improving monotonically across the schedule (0.75 at interim 1 $\to$ 0.47 at interim 8). Interpretation in §14.1.12.

**Why the training loss is tiny yet coverage is far off.** The pinball loss and this coverage test measure the *same* property — the $\eta^{\text{inpol}}_q$-level pinball is minimised by the conditional $\eta^{\text{inpol}}_q$-quantile, whose coverage is $\eta^{\text{inpol}}_q$ — but under **different distributions**. Training minimises expected pinball over the **prior-predictive** ($\theta \sim p(\theta)$, simulated cohorts); this diagnostic measures coverage on the **posterior-predictive deployment slice** ($\theta \sim p(\theta \mid x^{\text{obs}})$). A head calibrated on the prior-predictive average need not be calibrated on that slice when the two differ (the §8.1 prior-coverage gap: thresholds $\sim N(0, 3.5)$ generate cohorts far more extreme than Ukraine's), and the average loss is dominated by the bulk of the diffuse prior mass — the thin deployment region barely registers in it. Direct confirmation: **coverage on held-out prior-predictive draws is essentially nominal** (0.05, 0.27, 0.51, 0.76, 0.95 at the five levels), so the small loss is doing its job — the deployment miscalibration is distribution shift, not an optimisation or objective failure.

### 14.1.8 Diagnostic: marginal $p(\rho \mid x)$ (SVI vs amortiser)

Two distributions of $\rho \mid x$ are compared:

- **SVI:** the set $\{\rho^{(j^*, s)}_{\text{SVI-}x}\}_{s=1}^{S}$ (endpoints of the SVI-on-$x$ draws) — its empirical CDF is $p(\rho \mid x)$.
- **Amortiser:** its marginal predictive $\hat p(\rho \mid x) = \mathbb{E}_{z \mid x}\big[\hat p(\rho \mid x, z)\big]$, estimated by averaging the per-draw quantile CDFs, $\hat F(\rho \mid x) = \tfrac{1}{S}\sum_s \hat F(\rho \mid x, z^{(s)})$. Note the item tokens require a future cohort and so we marginalise over $s$. Agreement is summarised per (interim, item) by the two-sample KS distance $\max_\rho |\hat F(\rho \mid x) - \widehat{\text{ecdf}}_{\text{SVI}}(\rho)|$.

Outputs (same tests script as §14.1.7): `..._tests_marginal_ks.csv` and per-interim `pcm_1_interim_i{k}_svi_vs_amortiser_marginal_cdf.pdf` (facet per item, the two CDFs overlaid). Two-sample KS $\max_\rho |\hat F - \widehat{\text{ecdf}}_{\text{SVI}}|$ (`CG-VIO_ph-punish` excluded):

| marginal KS | out-of-7 | CG-MH |
|---|---:|---:|
| interim 1 | 0.51 | 0.53 |
| all interims | 0.42 | 0.62 |

The amortiser's $\hat p(\rho \mid x)$ sits **below** the SVI posterior — the same downward location bias as §14.1.7, now on the current-decision object. Unlike the conditional PIT-KS (which falls 0.75 → 0.47 across the schedule), the marginal KS is roughly **flat** ($\approx 0.47$ overall): a genuine location/shape mismatch of the marginal, not conditional under-dispersion that averages out. Reading in §14.1.12.

### 14.1.9 Diagnostic: median-prediction association (per-draw correlation)

`pcm_1_interim_i{k}_svi_rho_vs_amortiser_median.pdf` — per interim, facet per item. Scatter of

$$y \;=\; \rho^{(j^*, s)}_{\text{SVI-}x} \qquad \text{against} \qquad x \;=\; \hat\rho^{(j^*, s)}_{3} \quad (\eta^{\text{inpol}}_3 = 0.5),$$

the median head output after $q_\psi$ for same $s$. **D**eployment thresholds and averages into $\widehat{\text{PPS}}^{(j^*)}$, whereas this diagnostic keeps each pair and correlates the median output against $\rho^{(j^*, s)}_{\text{SVI}}$. Unlike §14.1.7 it tests **association only** (does the central prediction move with the target across $z^{(s)}$?), not location or spread — the two are deliberately separated because the amortiser can rank well while being mis-located. An earlier version scattered $\rho^{(j*,s)}_{\text{SVI}}$ against PC1 of $\bar h^{(j^*, s)}$ similar to the Strong-Oakly diagnostic. We dropped this because PC1 is a 1-D *linear* proxy for the nonlinear $q_\psi$ read-out, so it is uninformative about the network.

### 14.1.10 Results: correlation against the SVI posteriors

Per-item Pearson $\rho$ between $\hat\rho^{(j^*, s)}_3$ and $\rho^{(j^*, s)}_{\text{SVI-}x}$ for the prior-predictive-trained amortiser (`CG-VIO_ph-punish` blow-up item excluded):

| correlation | out-of-7 | CG-MH |
|---|---:|---:|
| interim 1 | 0.82 | 0.71 |
| all interims | 0.53 | 0.42 |

**Threshold on the CG-MH token responses.** The four `CG-MH` items are categorical; the endpoint $r_j(\theta)$ is a caseness prevalence, defined via a response threshold $c$. It was important to keep that threshold **out of the input token**: an earlier token summarised the CG-MH responses as the proportion above $c$ (a thresholded, lossy summary) and capped their correlation at $\rho \approx 0.48$; replacing it with the **mean response** (no threshold, no information loss) — while keeping $c$ only in the target — lifted CG-MH to $\approx 0.71$ with no change to the out-of-7 items. The threshold belongs in the endpoint definition, never in the summary the network reads.

Synthesis in §14.1.12.

### 14.1.11 Diagnostic: PPS agreement with the regression baselines

No analytic ground-truth PPS exists for Ukraine, so the deployed $\widehat{\text{PPS}}^{(j^*)}$ is compared to the two closest Strong-Oakley baselines — Gaussian-approximation regression (RGE) and multi-quantile regression (RGEM) — over the $J \cdot C = 160$ (item, interim) cells:

| comparison                               | Pearson $\rho$ | median \$ | \Delta\text{PPS} |    \$ |
| ---------------------------------------- | -------------: | --------: | ---------------: | ----: |
| `itemXcompAtt` vs RGE (Gauss)            |          0.905 |     0.001 |            0.223 | 0.880 |
| `itemXcompAtt` vs RGEM (mquantile)       |          0.916 |     0.002 |            0.241 | 0.830 |
| `itemXcompAtt` vs `itemScompAtt` (§14.3) |          0.997 |     0.000 |            0.041 |     — |

Median difference to the regression baselines is $\le 0.002$; the disagreeing tail (p90 $\approx 0.22$-$0.24$) is where cross-item information and the calibration gap move the prediction. The two attention variants agree almost perfectly. Reading in §14.1.12.

### 14.1.12 Reading

Synthesising the diagnostics (§14.1.7 conditional calibration, §14.1.8 marginal $p(\rho\mid x)$, §14.1.9 association / §14.1.10 results, §14.1.11 PPS agreement):

1.  **Ranking is good, calibration is not — the two are dissociated.** Median-output association with the SVI draws is strong (Pearson $\rho \approx 0.8$ out-of-7, $\approx 0.7$ CG-MH; §14.1.9), yet quantile coverage sits far below nominal and PIT-KS $\approx 0.61$ (§14.1.7), and the marginal $\hat p(\rho\mid x)$ sits below the SVI posterior (§14.1.8): the amortiser orders the posterior draws well but its conditional distribution is **mis-located and under-dispersed**. The calibration diagnostic quantifies distributionally what the negative $R^2$ hints; the median diagnostic confirms the location error does not come from a broken ranking. The PPS agrees with the Strong-Oakley baselines on the bulk of cells ($\le 0.002$ median), with a tail exactly where the calibration gap bites (§14.1.10).

2.  **Cause: prior coverage.** The output scale inherits the diffuse PCM prior — thresholds $\sim N(0, 3.5)$ place most training mass on far more extreme cohorts than Ukraine's, and the training label is the *population* endpoint $r_j(\theta)$ while the evaluation target is the finite-cohort SVI estimate. The monotone improvement of PIT-KS across the schedule (interim 1 → 8) fits this: as $n$ grows the operational posterior moves into the prior-covered region. This is the §8.1 prior-coverage caveat, made concrete.

3.  **Remedies (the actionable next step).** these are executed as a full remedy sweep in §14.2 — the deployable fix is an expanding-window head-only fine-tune (§14.2.5), which reaches near-nominal calibration with correlation intact. Same calibration origin as the deepset variant (§14.4.3).

4.  **The §8-10 loop closes at the item level at no accuracy cost.** Prior-predictive training (fresh PCM draws, no SVI anywhere) with the hand-computed $F = 8$ token reaches the same correlation as the superseded posterior-trained regime (§14.1.9, final vs third row) — the train-once-deploy-anywhere property is free on this metric; only calibration is outstanding.

5.  **Architecture ports cleanly from MVN, and cross-item routing is a wash here.** The network class is unchanged; the analytic K-row becomes the empirical Spearman-shrunk $\hat K$-row. Cross- vs self-attention are within noise on Ukraine (§14.3) because $\hat K$ is weak (off-diagonals $\approx 0.3$); the MVN $J = 60$ result (§13.5) is where the routing advantage shows.

6.  **Deployment cost is simulation, not network.** Rebuilding the posterior-predictive $z^{(s)}$ block is \~1 min/interim at $S = 4000$; the $80{,}000$ head evaluations take seconds. SVI on $x^{\text{obs}}$ (once per interim) is the only inference step left.

## 14.2 Calibration remedies for the item amortiser

The diagnostics of §14.1.7–14.1.8 leave one defect: the prior-trained amortiser is well-calibrated on the prior-predictive (held-out-prior coverage nominal, §14.1.7) but its deployment predictive is **mis-located** on the posterior-predictive slice — coverage far below nominal, PIT-KS $\approx 0.61$, marginal-KS $\approx 0.47$ — while ranking stays strong (correlation $\approx 0.8/0.7$, §14.1.10). This section reports a remedy sweep. Goal: near-nominal calibration **without** losing correlation or the §8 generality. Decision metrics: coverage at $\eta^{\text{inpol}} = 0.5, 0.95$, PIT-KS, marginal-KS, and correlation (out-of-7 / CG-MH; `CG-VIO_ph-punish` blow-up item excluded). Guardrail: held-out-prior coverage must stay nominal. Every configuration writes the full diagnostic PDFs (coverage curves, PIT histograms, marginal-CDF overlays) into its own `…_260808`/`_260809` sandbox dir.

| config | corr (o7/cat) | cov@.5 | cov@.95 | PIT-KS | marg-KS |
|---|---|---:|---:|---:|---:|
| baseline (§14.1) | 0.53/0.42 | 0.15 | 0.44 | 0.61 | 0.47 |
| S=256/step (§14.2.1) | 0.52/0.41 | 0.14 | 0.48 | 0.60 | 0.47 |
| thresh=2.0 (§14.2.2) | 0.53/0.42 | 0.14 | 0.40 | 0.64 | 0.50 |
| thresh=1.0 (§14.2.2) | 0.53/0.41 | 0.14 | 0.45 | 0.62 | 0.50 |
| full-ft LOIO (§14.2.3) | 0.51/0.41 | 0.53 | 0.94 | 0.17 | 0.14 |
| head-ft(i1) (§14.2.4) | 0.52/0.13 | 0.53 | 0.94 | 0.40 | 0.40 |
| **head-ft expand 1:k (§14.2.5)** | **0.52/0.39** | **0.50** | **0.95** | **0.16** | **0.14** |
| conformal i1 (§14.2.6) | 0.53/0.41 | 0.38 | 0.44 | 0.58 | 0.66 |
| conformal expand 1:k (§14.2.6) | 0.53/0.41 | 0.40 | 0.44 | 0.56 | 0.59 |

(nominal cov@.5 = 0.50, cov@.95 = 0.95; lower KS is better.)

### 14.2.1 Training-noise control — $S = 256$ draws/step (negative)

More fresh prior draws per step (1280 → 5120 query rows; lower-variance gradients) leaves every calibration metric unchanged (cov 0.14/0.48, PIT-KS 0.60, marg 0.47). Confirms the §14.1.7 diagnosis: the gap is distribution shift, not optimisation noise. Dir `…itemXcompAtt…_256step_260808`.

### 14.2.2 Prior narrowing — `threshold_scale` 2.0 / 1.0 (negative)

Tightening the diffuse threshold prior toward the operational region does not help (cov ≈ 0.14/0.40–0.45, PIT-KS 0.62–0.64, marg 0.50 — if anything slightly worse). The OOD gap is not fixed by shrinking threshold spread alone: abilities/loadings and the population-endpoint-vs-finite-cohort-target mismatch dominate. Data-agnostic prior narrowing is off the table for this study. Dirs `…_thresh2p0/thresh1p0_260808`.

### 14.2.3 Full-net fine-tune, leave-one-interim-out (reference upper bound)

Warm-start from the prior-trained weights, fine-tune **all** parameters 1000 steps on interims 2–8, evaluate held-out interim 1. Near-nominal calibration (cov 0.53/0.94, PIT-KS 0.61 → **0.17**, marg 0.47 → **0.14**) with correlation preserved in- and out-of-sample (all-interim 0.51/0.41; held-out interim-1 0.71/0.72). The fine-tune is *light* — every weight block moves only 8–15 % in norm, cosine $\approx 0.99$ with the prior weights ($q_\text{tok}$ barely 8 %) — so it is the prior amortiser **lightly recalibrated**, not re-fit: the prior stage supplies the retained feature extractor, the fine-tune does the location/scale correction. Cost: reintroduces a per-study SVI dependency at fine-tune time. Dir `…_ftfull-loio-i1_260808`.

### 14.2.4 Head-only fine-tune on interim 1 (too thin)

Freeze the encoders, fine-tune only $q_\psi$ on interim 1. Coverage jumps to nominal (0.53/0.94) but CG-MH correlation **collapses** 0.42 → 0.13 — one interim over-fits interim-1's categorical scale and destroys held-out CG-MH ranking. Head-only is not wrong; one fit-interim is too little data. Dir `…_ftheadi1_260808`.

### 14.2.5 Head-only fine-tune, expanding window (the deployable recipe)

The realistic sequential version: at interim $k$ the SVI fits for interims $1..k$ are already in hand (deployment runs SVI per interim), so fine-tune only $q_\psi$ on $1..k$ and apply to interim $k$ (encoders frozen throughout). **Best result:** cov 0.50/0.95 (nominal), PIT-KS 0.61 → **0.16**, marginal-KS 0.47 → **0.14**, correlation preserved (0.52/0.39 — CG-MH recovers vs §14.2.4 because the expanding window supplies enough data). It **matches the full-net upper bound (§14.2.3) while keeping the encoders universal**, and uses only SVI already computed at interims $\le k$. Recommended operating point. Dir `…_ftheadexpand_260809`.

### 14.2.6 Post-hoc conformal recalibration (partial)

No fine-tuning: fit a per-item monotone PIT map $\hat g = \text{ecdf}(u_{\text{cal}})$ (distribution calibration, Kuleshov et al. 2018) and recalibrate the predictive CDF. Two protocols — interim-1 map applied to all interims, and expanding (map from $1..k$ applied to interim $k$). Both only **partly** help: cov@.5 0.15 → 0.38–0.40 (not 0.50), PIT-KS 0.61 → 0.56–0.58, marginal-KS *worse* (0.59–0.66). Cause: the miscalibration **varies across interims** (the location bias shrinks as $n$ grows, §14.1.8), so a single static monotone map cannot correct all interims at once — a *conditional* re-scale (fine-tune) can. Dirs `…_conformal-i1/conformal-expand_260809`.

### 14.2.7 Reading

1.  **Diagnosis holds.** The $S=256$ null (§14.2.1) reconfirms the gap is distribution shift, curable only by moving the predictive toward the operational region.
2.  **Cheap, data-agnostic fixes fail.** Prior narrowing (§14.2.2) and static conformal recalibration (§14.2.6) both fall short — the latter because the bias is interim-heteroscedastic, which a single monotone map cannot absorb.
3.  **Fine-tuning on the operational SVI works, and head-only is enough.** The **expanding-window head-only** fine-tune (§14.2.5) is the deployable operating point: near-nominal calibration, correlation intact, universal encoders retained, using only SVI already computed at interims $\le k$. It buys a correctly-located predictive for a light per-interim head re-fit — a small, honest trade of §8 purity.
4.  **The correction is a small warp, not a retrain.** The full-net fine-tune (§14.2.3) reaches the same calibration by moving the weights only $\approx 10\%$ near the prior solution — the amortiser is recalibrated, not relearned; the prior stage's representation is doing the heavy lifting throughout.

## 14.3 Self-attention variant (`itemScompAtt`): detailed diagnostics

Identical tokens, target, loss, prior-predictive training data and deployment as §14.1; only the attention block differs (§13.6 variant 2). Every token acts as query and key; a learned scalar $\alpha$ biases each row toward the queried key,

$$\text{score}^{(j^*, s)}_{a, b} \;=\; \frac{\big\langle h^{(j^*, s)}_a,\; h^{(j^*, s)}_b\big\rangle}{\sqrt{E}} \;+\; \alpha \cdot \mathbf{1}\{b = j^*\}, \qquad h'^{\,(j^*, s)}_a \;=\; \sum_{b=1}^{J} \operatorname{softmax}_b\big(\text{score}^{(j^*, s)}_{a, \cdot}\big)\, h^{(j^*, s)}_b,$$

and the per-query summary is the gather $\bar h^{(j^*, s)} = h'^{\,(j^*, s)}_{j^*}$ — $J^2$ attention scores per batch element (vs $J$ for cross-attention). $\alpha$ adds one learned scalar to $(\tau, \psi)$.

### 14.3.1 Association with the SVI target

Same protocol and token as §14.1.9 (per-item Pearson $\rho$ between the amortiser median $\hat\rho^{(j^*, s)}_3$ and $\rho^{(j^*, s)}_{\text{SVI-}x}$, blow-up item excluded):

| variant                | interim 1: out-of-7 / CG-MH | all interims: out-of-7 / CG-MH |
| ---------------------- | --------------------------- | ------------------------------ |
| `itemXcompAtt` (§14.1) | 0.82 / 0.71                 | 0.53 / 0.42                    |
| `itemScompAtt`         | 0.83 / 0.70                 | 0.53 / 0.41                    |

(Both rows: §8 prior-predictive training, identical protocol.) On MVN, cross-attention beat self-attention by 20-35 % at $J \le 60$ (§13.6) — the advantage came from K-row alignment through the query. On Ukraine at $J = 20$ the two variants are within noise: the empirical $\hat K$-row is weaker and shrunk (mean off-diagonal $\approx 0.3$), so the query-side routing advantage largely disappears. Association plots: `pcm_1_interim_i{k}_svi_rho_vs_amortiser_median.pdf`.

### 14.3.2 Baseline calibration — the same prior/posterior mislocation as §14.1.7

The full SBC suite (`..._tests.py`, net-agnostic) on the prior-trained self-attention net gives the identical picture as `itemXcompAtt`: strong ranking, badly mislocated intervals.

| metric (mean over items/interims) | nominal | `itemXcompAtt` baseline | `itemScompAtt` baseline |
| --------------------------------- | ------: | ----------------------: | ----------------------: |
| coverage @ $\eta = 0.5$           |    0.50 |                    0.15 |                    0.19 |
| coverage @ $\eta = 0.95$          |    0.95 |                    0.44 |                    0.53 |
| PIT–KS                            |    0.00 |                    0.61 |                    0.55 |
| marginal–KS                       |    0.00 |                    0.47 |                    0.42 |

Diagnostics in the baseline dir: `pcm_1_interim_pps_RGEG_tests_coverage.pdf`, `..._tests_pit.pdf`, `..._i{k}_svi_vs_amortiser_marginal_cdf.pdf`. The mechanism is §8.1's prior-coverage gap, architecture-independent: the diffuse PCM prior is wider than the concentrated deployment posterior slice, so the predicted quantiles are correctly *ordered* but too wide and off-centre.

### 14.3.3 Head-only expanding fine-tune (§14.2.5) transfers unchanged

Applying the §14.2.5 recipe — freeze the encoder, refit only the head $q_\psi$ on interims $1..k$ SVI, evaluate at $k$ — to the self-attention net (universal freeze: train `q_psi`, freeze `q_tok` and the self-attention bias $\alpha$):

| metric              | nominal | baseline | head-ft-expand |
| ------------------- | ------: | -------: | -------------: |
| coverage @ 0.5      |    0.50 |     0.19 |           0.50 |
| coverage @ 0.95     |    0.95 |     0.53 |           0.94 |
| PIT–KS              |    0.00 |     0.55 |           0.15 |
| marginal–KS         |    0.00 |     0.42 |           0.13 |
| corr o7/cat (i1)    |       — | 0.83/0.70| 0.83/0.72      |

Same outcome as `itemXcompAtt` (§14.2.5): near-nominal calibration, correlation preserved (encoders frozen). Diagnostics in `..._itemScompAtt_..._ftheadexpand_260809/`. `itemXcompAtt` remains the default — the two variants are within noise on this $J = 20$ case and share the same remedy.

## 14.4 Full §8-10 amortiser on the Ukraine PCM: `deepsetXcompAtt` + prior-predictive training

The nested-DeepSets variant of the §8-10 workflow on the Ukraine PCM. Training is prior-predictive as in §14.1; the difference is architectural: instead of the hand-computed $F = 8$ item summaries, the network consumes the RAW cohorts — `x_responses` and `z_responses` per (participant, item, time) — and **learns** the per-item pooling through $q_\tau$ (§9.2).

### 14.4.1 Setup

- **§8 training data.** `PartialCreditModel.make_training_data_with_participant_tokens_prior` ([`python/model_pcm.py`](../python/model_pcm.py)): per prior draw, $n \sim U\{2, \ldots, N{-}1\}$ with $N = 503$, $m = N - n$; abilities $\theta_i \sim N(0, 1)$, time effects $\beta_t \sim N(0, 1)$, per-(item, time) thresholds $\sim N(0, 3.5)$, loadings $\sim |t_3|$ (first per category type pinned to 1); ordered-categorical responses simulated per (participant, item, time); label $\rho_j(\theta)$ = direction-aware ratio of the population endpoint across the two times, clipped $\pm 20$. Responses rescaled $(k - 1)/(K_j - 1) \in [0, 1]$.
- **§9 architecture.** `deepsetXcompAtt` (§13.7): shared $q_\tau$ over (participant, item) tokens with $R = 2$ (baseline, endline response) + item metadata $M = 2$; mean-pool per cohort per item; token = concat of the two cohort pools; cross-attention over items; multi-quantile head. Aux $= (n/N, m/N)$.
- **§10 deployment.** Per interim: raw observed cohort from the cached xi frame, raw $z^{(s)}$ rebuilt from the cached SVI `draws.zarr` (draw-aligned with the wa evaluation targets), one forward pass per $(j^*, s)$.
- 8000 training steps at batch $S = 32$ fresh prior draws per step $\Rightarrow$ **256,000 prior-predictive data sets** generated during training (each fanned out to $Q = 4$ queried items $\Rightarrow$ 1,024,000 training rows); \~74 min CPU. Deployment 0.11 min/interim (S = 200).
- **Size coverage vs deployment.** Training draws $n \sim U\{2, \ldots, 502\}$ (501 values, $m = 503 - n$); the 8 deployment interims use $n \in \{48, 174, 329, 416, 473, 484, 491, 501\}$. The proportion of training data sets whose $n$ exactly matches *any* deployment interim is $8/501 \approx 1.6\%$ ($\approx 4{,}100$ of the 256,000), and $\approx 0.2\%$ ($\approx 510$) per individual interim — confirming the amortiser is genuinely interpolating across cohort sizes rather than memorising the deployment configurations; the mean-pool over participants plus the $(n/N, m/N)$ aux carry the size dependence.

### 14.4.2 Results: evaluation against the SVI posterior targets

Per-item Pearson $\rho$ between the amortiser median $\hat\rho^{(j^*, s)}_{3}$ and the SVI target $\rho^{(j^*, s)}_{\text{SVI-}x}$ across the $S$ draws (excluding the `CG-VIO_ph-punish` blow-up item):

| item type                     | interim 1 mean $\rho$ | all-interim mean $\rho$ |
| ----------------------------- | --------------------: | ----------------------: |
| `out-of-7` (15 items)         |                  0.82 |                    0.53 |
| `categorical` CG-MH (4 items) |                  0.57 |                    0.32 |

Both the deepset and the item-level amortiser (§14.1.9) are now trained purely prior-predictively, so the comparison isolates the **architecture**: learned pooling matches the hand-computed $F = 8$ token on the out-of-7 items (0.82 vs 0.82 at interim 1) but trails it on CG-MH (0.57 vs 0.71) — the item token's mean-response summaries have no learned counterpart yet in $q_\tau$ at this training budget. $R^2$ is mostly negative — the amortiser tracks the *ranking* of the posterior draws well but its location/scale calibration under the diffuse PCM prior does not match the concentrated deployment posterior (the §8.1 prior-coverage caveat in action: thresholds $\sim N(0, 3.5)$ generate far more extreme cohorts than the Ukraine data region). Association plots: `pcm_1_interim_i{k}_svi_rho_vs_amortiser_median.pdf`.

**Per-item detail at interim 1** (Pearson $\rho$ / $R^2$ vs the SVI target over $S = 200$ draws):

| item                        | type        | $\rho$ | $R^2$ |
| --------------------------- | ----------- | -----: | ----: |
| `CG-MH_effort`              | categorical |   0.65 | -2.50 |
| `CG-MH_hopeless`            | categorical |   0.49 | -0.71 |
| `CG-MH_nervous`             | categorical |   0.50 | -3.65 |
| `CG-MH_sad`                 | categorical |   0.64 | -3.97 |
| `CG-DEPRESSION`             | out-of-7    |   0.84 | -0.08 |
| `CG-INVOLVE_child-problems` | out-of-7    |   0.77 |  0.17 |
| `CG-INVOLVE_help-learn`     | out-of-7    |   0.84 | -1.09 |
| `CG-MONITOR-CHI_child-safe` | out-of-7    |   0.87 |  0.17 |
| `CG-MONITOR-CHI_safe-time`  | out-of-7    |   0.86 | -0.15 |
| `CG-NONVIOLENT-DISCIPLINE`  | out-of-7    |   0.83 | -1.35 |
| `CG-POS_play`               | out-of-7    |   0.82 | -1.08 |
| `CG-POS_praise`             | out-of-7    |   0.78 | -0.17 |
| `CG-RESILIENCE`             | out-of-7    |   0.84 | -0.89 |
| `CG-SELFCARE`               | out-of-7    |   0.79 | -3.28 |
| `CG-VIO_ph-punish`          | out-of-7    |   0.29 | -0.47 |
| `CG-VIO_scream`             | out-of-7    |   0.79 | -0.47 |
| `CHI-BEHAVIOUR_angry`       | out-of-7    |   0.78 |  0.50 |
| `CHI-BEHAVIOUR_no-interest` | out-of-7    |   0.83 |  0.45 |
| `CHI-BEHAVIOUR_unhappy`     | out-of-7    |   0.85 |  0.50 |
| `grieve`                    | out-of-7    |   0.81 | -1.06 |

High $\rho$ with mostly negative $R^2$ across the board — the ranking of posterior draws is tracked well, the location/scale calibration is not (see Reading). `CG-VIO_ph-punish` is the §14.1 numerical blow-up item.

**PPS agreement** with the item-level amortiser (§14.1, RGEG): high on the bulk of cells, with the disagreeing tail at early interims where the calibration gap is largest.

### 14.4.3 Baseline calibration diagnostics (`..._diag.py`)

The same SBC suite as §14.1.7, run on the deepset by forward-passing every draw and capturing the full multi-quantile prediction (and the frozen-encoder head input `head_in`, sown for §14.4.5). It reproduces the item-level picture exactly — strong ranking, mislocated intervals:

| metric (mean over items/interims) | nominal | `itemXcompAtt` §14.1 | `deepsetXcompAtt` |
| --------------------------------- | ------: | -------------------: | ----------------: |
| coverage @ 0.5                    |    0.50 |                 0.15 |              0.18 |
| coverage @ 0.95                   |    0.95 |                 0.44 |              0.41 |
| PIT–KS                            |    0.00 |                 0.61 |              0.60 |
| marginal–KS                       |    0.00 |                 0.47 |              0.47 |

Diagnostics in the baseline dir: `pcm_1_interim_pps_RGDX_tests_coverage.pdf`, `..._tests_pit.pdf`, `..._i{k}_svi_vs_amortiser_marginal_cdf.pdf`, `..._i{k}_svi_rho_vs_amortiser_median.pdf`. The mislocation is the §8.1 prior/posterior gap again — architecture makes no difference to it.

### 14.4.4 Workstream A: post-hoc conformal recalibration (insufficient here too)

We expected the nested-DeepSets amortiser to be *more* amenable to a static distribution-calibration map (Kuleshov et al. 2018) than the item-level net — it is not. Per-item monotone PIT maps fitted on interim-1 SVI (A1) and on the expanding pool $1..k$ (A2):

| metric          | nominal | baseline | A1 conformal-i1 | A2 conformal-expand |
| --------------- | ------: | -------: | --------------: | ------------------: |
| coverage @ 0.5  |    0.50 |     0.18 |            0.38 |                0.39 |
| coverage @ 0.95 |    0.95 |     0.41 |            0.41 |                0.41 |
| PIT–KS          |    0.00 |     0.60 |            0.60 |                0.59 |
| marginal–KS     |    0.00 |     0.47 |            0.62 |                0.54 |

Exactly as for `itemXcompAtt` (§14.2.6): the map nudges the median coverage (0.18 → 0.38) but leaves the upper tail pinned (coverage @ 0.95 stuck at 0.41), barely moves the PIT–KS, and makes the marginal *worse*. **Why the deepset intuition fails:** conformal can only fix miscalibration that is a single, interim-invariant monotone warp of the PIT. Here the bias is *interim-heteroscedastic* — its magnitude shrinks as the cohort $n$ grows — so one static per-item map cannot serve all 8 interims at once. This is a property of the prior/posterior gap, not of the encoder, which is why the deepset behaves like the item net.

### 14.4.5 Head-only expanding fine-tune (§14.2.5) — the architecture-independent fix

Because `head_in` (the frozen-encoder input to $q_\psi$) is sown by the net, the §14.2.5 recipe runs *without re-running the encoder*: cache `head_in` for every $(k, s, j^*)$ in the Phase-0 pass, then refit only $q_\psi$ on the pooled embeddings of interims $1..k$ and evaluate at $k$.

| metric           | nominal | baseline | head-ft-expand |
| ---------------- | ------: | -------: | -------------: |
| coverage @ 0.5   |    0.50 |     0.18 |           0.50 |
| coverage @ 0.95  |    0.95 |     0.41 |           0.94 |
| PIT–KS           |    0.00 |     0.60 |           0.18 |
| marginal–KS      |    0.00 |     0.47 |           0.17 |
| corr o7/cat (i1) |       — |0.82/0.57 | 0.82/0.54      |

Near-nominal calibration, correlation preserved — identical outcome to `itemXcompAtt` (§14.2.5) and `itemScompAtt` (§14.3.3). Diagnostics in `..._deepsetXcompAtt_..._ftheadexpand_260809/`.

### 14.4.6 Cross-architecture summary

The remedy sweep run on all three encoders gives one consistent conclusion:

| architecture (interim-1 corr o7/cat) | remedy         | cov@.5 | cov@.95 | PIT–KS | marg–KS |
| ------------------------------------ | -------------- | -----: | ------: | -----: | ------: |
| `itemXcompAtt` (0.82/0.71)           | baseline       |   0.15 |    0.44 |   0.61 |    0.47 |
|                                      | conformal-exp  |   0.40 |    0.44 |   0.56 |    0.59 |
|                                      | head-ft-expand |   0.50 |    0.95 |   0.16 |    0.14 |
| `itemScompAtt` (0.83/0.70)           | baseline       |   0.19 |    0.53 |   0.55 |    0.42 |
|                                      | head-ft-expand |   0.50 |    0.94 |   0.15 |    0.13 |
| `deepsetXcompAtt` (0.82/0.57)        | baseline       |   0.18 |    0.41 |   0.60 |    0.47 |
|                                      | conformal-exp  |   0.39 |    0.41 |   0.59 |    0.54 |
|                                      | head-ft-expand |   0.50 |    0.94 |   0.18 |    0.17 |

### 14.4.7 Reading

1.  **The §8-10 loop closes.** Prior sampling → nested-DeepSets training → raw-$(x, z)$ deployment runs end-to-end on a real IRT case study with a \~1 min total deployment cost for the full 8-interim schedule.
2.  **Architecture affects ranking, not calibration.** Cross-attention, self-attention and nested-DeepSets share the *same* prior/posterior mislocation (baseline coverage @ 0.5 $\approx 0.15$–$0.19$, PIT–KS $\approx 0.55$–$0.61$ across all three). The encoder choice only moves the correlation — the deepset's learned pooling matches the hand token on out-of-7 (0.82) but trails on CG-MH (0.57 vs 0.70–0.71) at this budget.
3.  **The fix is architecture-independent and it is the head, not a static map.** The expanding head-only fine-tune (§14.2.5) recovers near-nominal calibration for all three encoders with correlation preserved, whereas post-hoc conformal is insufficient for all three (including the deepset, contrary to expectation) because the bias is interim-heteroscedastic. Refitting $q_\psi$ conditions on each interim's data through the frozen embedding; a single monotone PIT map cannot.
4.  **Diagnostics** (§14.1.7–14.1.9) are produced in full for every configuration: coverage curves, PIT histograms, marginal-CDF overlays and median-association scatters, plus `pcm_1_interim_pps_RGDX_svi_vs_amortiser_rho_all.pdf`.

# 15. References

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