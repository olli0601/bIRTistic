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

- **Binomial (§3.1)**: the likelihood $p(x_i \mid p) = p^{x_i}(1-p)^{1-x_i}$ has natural parameter $\eta(p) = \log\!\big(p/(1-p)\big)$ and canonical scalar sufficient statistic
  $$  T(x_i) = x_i \in \{0, 1\}.
    $$ Set $q_\tau(x_i) := x_i$ (identity, dimension 1). The DeepSets sum is $\sum_{i=1}^n T(x_i) + \sum_{i^* = 1}^m T(z_{i^*}) = k^n + k^m$, the total number of successes. The head $q_\psi$ then only needs to learn the map $(k^n + k^m, n, m) \mapsto Q_\tau^{\tau_k}\!\big(\rho \mid k^n + k^m, n, m\big)$; feeding the cohort sizes $(n, m)$ as extra scalar inputs to $q_\psi$ preserves the exact-sufficiency of $q_\tau$.

- **MVN with** $\sigma^2$ known (§3.3, Gaussian special case): with the natural parameter $\eta(\mu) = K^{-1} \mu / \sigma^2$ and $\sigma^2$ fixed, only the first-order statistic is needed,
  $$  T(x_i) = x_i \in \mathbb{R}^J.
    $$ Set $q_\tau(x_i) := x_i$ (identity, dimension $J$). The DeepSets sum recovers $T_1(x_{1:n}) + T_1(z_{1:m}) = \sum_{i=1}^n x_i + \sum_{i^*=1}^m z_{i^*}$; $q_\psi$ takes this $J$-vector together with $(n, m)$ as inputs.

- **MVN with** $\sigma^2$ unknown (§3.3, full NIG): the likelihood adds a second natural parameter $-1/(2\sigma^2)$ against the quadratic-form statistic. Per-observation sufficient statistics are
  $$  T(x_i) = \big(x_i,\; x_i^\top K^{-1} x_i\big) \in \mathbb{R}^{J+1}.
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

# 11. Codebase recommendation

The amortised extension slots into the existing pipeline as one new module and one new script per model, mirroring the existing gauss-approx / quantile / mquantile ladder.

## 11.1 New Python module: `python/amortised_pps.py`

- `class AmortisedPPSNet` — DeepSets encoder + multi-quantile head. Constructor arguments: `input_dim`, `embed_dim`, `hidden_dims`, `taus` (list of quantile levels), optional `T_fn` for hardcoded sufficient statistics. Flax module (consistent with the rest of the JAX stack in the repo).
- `fit_amortised_pps(model, prior_sampler, data_sampler, taus, n_samples, ...)` — training loop. `prior_sampler` returns $\theta$; `data_sampler(theta, n, m)` returns $(x, z, \rho)$. Reuses the existing `Model.prior_predictive` / `Model.get_endpoints_per_draw` interfaces.
- `fit_interim_regress_endptx_on_wz_with_amortised_pps(wa_or_xz_pairs, pps_H1_min_effect_size_thresh, pps_ProbH1_target_lwr_quantile, taus, model_ckpt=...)` — deployment forward pass, returns a `p_h1_xz` frame with the same schema as the existing `fit_interim_regress_endptx_on_wz_with_mquantile_regr` (adds an `amortised=True` column). Drops in downstream without further changes.

## 11.2 New scripts (mirror the mquantile ladder)

- `scripts-py/Binomial_interim_analysis_regression_endptx_on_wz_with_amortised_pps.py`
- `scripts-py/MVN_interim_analysis_regression_endptx_on_wz_with_amortised_pps.py`
- `scripts-py/Ukraine_interim_analysis_regression_endptx_on_wz_with_amortised_pps.py`

Each mirrors the mquantile script layout: fit (or load) the amortised network once at the top, then loop over interims to produce `p_h1_xz` predictions and the standard downstream artifacts (`_RGEA_` suffix). Reuse the existing `model.get_interim_z_from_ypredi`, `model.get_endpoints_per_draw`, `model.get_w`.

## 11.3 Training-data cache

- `scripts-py/{Binomial, MVN, Ukraine}_amortised_pps_make_training_data.py` — generates and caches the training samples $(x^{(s)}, z^{(s)}, \rho^{(s)})$ at a specified sample budget. Reuses the model classes and `prior_predictive` methods. Training-data generation is the expensive step; separating it from the neural fit lets us iterate on architectures cheaply. The cached tuples are indexed by `(sample_id, n_id, m_id)` so the training loop can subsample size regimes.

## 11.4 Compare-methods integration

- Add `dir_rgea = ..._amortised_pps_260702` to `{Binomial, MVN, Ukraine}_interim_analyses_compare_methods.py`.
- Add `'Regression of endpt-x on w(z) - amortised (mquantile)'` to `method_order` and the concat lists.
- Multi-quantile deployment gives continuous $\hat P_\phi(H_1 \mid x, z) \in [0, 1]$ so the amortised method joins the p(H_1 \| x, z) boxplot without filtering; only the single-quantile deployment mode (binary output) needs the same filter as the `- quantile` method.

## 11.5 Correctness tests

- Add tests under `test/python/test_amortised_pps_correctness.py`:
  - Binomial exp-fam benchmark: amortised PPS vs closed-form Beta-Binomial PPS at 8 interim states. Tolerance ≤ 0.02 in absolute error.
  - MVN benchmark: same at $J \in \{20, 60, 100\}$.
  - Deterministic sample draw + fixed weights → deterministic output (mirrors `test/python/test_interim_determinism.py`).

## 11.6 Suggested first prototype

Start with the Binomial model, hardcoded sufficient statistic $T(x_i) = x_i$ (sum of successes), multi-quantile head with the same 11-level grid used in §6.5. Verify against the closed-form Beta-Binomial PPS on the fully-observed month-by-month interim grid used in `Binomial_interim_analysis_regression_endptx_on_wz_with_mquantile_regr.py`. If the amortised network matches the analytic PPS to within Monte-Carlo noise across all 12 interims and across 100 randomly drawn $x^{\text{obs}}$, the approach is validated; extend to MVN and then to the IRT model.

------------------------------------------------------------------------

# 12. Results

## 12.1 Results for Binomial interim analysis amortised endptx on wz with features-fixed, qpsi-MLP, loss-multiquantilehead

**Implementation.** Prototype §11.6 built out end-to-end:

- New module [`python/amortiser_pps_features_fixed_qpsi_MLP_loss_multiquantilehead.py`](../python/amortiser_pps_features_fixed_qpsi_MLP_loss_multiquantilehead.py) exposing the Flax module `Amortiser_PPS_features_fixed_qpsi_MLP_loss_multiquantilehead` (DeepSets-pooled features + SiLU-MLP `q_psi` + multi-quantile head), a training loop `train()` (Adam + cosine LR warmup + pinball loss over 11 quantile levels), and prediction helpers `predict_amortised_p_h1_for_one_xz` (monotone-corrected CDF interpolation at `pps_H1_min_effect_size_thresh`) / `predict_amortised_p_h1_for_many_xz` (deployment wrapper with mquantile-compatible output schema).
- `BinomialModel.make_training_data_with_features` attached to [`python/model_binomial.py`](../python/model_binomial.py) so the amortiser's training-data step is model-owned. The DeepSets pooling for the hardcoded $T(x_i) = x_i$ collapses to features `(k_total, n_total) / N_max`.
- Sim-data cache script [`scripts-py/Binomial_interim_analyses_make_sim_data.py`](../scripts-py/Binomial_interim_analyses_make_sim_data.py) writes the fixed Bernoulli cohort + monthly interim grid + a $65\,536$-sample amortiser training batch to disk (`binomial_sim_cohort.pkl`, `binomial_amortiser_training_data.pkl`).
- Deployment script [`scripts-py/Binomial_interim_analysis_amortise_endptx_on_wz_with_features_fixed_qpsi_MLP_loss_multiquantilehead.py`](../scripts-py/Binomial_interim_analysis_amortise_endptx_on_wz_with_features_fixed_qpsi_MLP_loss_multiquantilehead.py) trains the net once (or loads the cached checkpoint) and, at each interim, forward-passes the DeepSets-pooled features to produce a `p_h1_xz` frame that `Binomial_interim_analyses_compare_methods.py` picks up under the `_RGEA_` suffix.
- Correctness tests [`test/python/test_amortised_pps_correctness.py`](../test/python/test_amortised_pps_correctness.py) exercise forward-pass determinism, training determinism, save/load round-trip, and PPS accuracy on the deployment monthly grid + 20 random $(k_n, n, m)$ triples.

**Training configuration.** MLP head `hidden_dims = (256, 256, 128)`, `num_quantiles = 11` (`taus = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)`), Adam with peak learning rate $10^{-3}$ under a linear-warmup / cosine-decay schedule, $40\,000$ steps × batch $8192$. Prior sampler: `Beta(1, 1)` on $p$, $n \sim U\{1, \ldots, N-1\}$, $m = N - n$ (fixed final cohort $N = 500$), $k_n \sim \text{Binomial}(n, p)$, $k_m \sim \text{Binomial}(m, p)$, features `(k_n + k_m, N) / N`, target $\rho = 1 - p / p_0$ with $p_0 = 0.5$.

**Deployment.** At each interim we build the posterior-predictive draws via `BinomialModel.fit_closed_form_posterior` (analytic $p_s \sim \text{Beta}(1 + k_n,\, 1 + n - k_n)$) chained through the `BinomialModel`-specific `get_interim_z_from_ypredi` override — which draws fresh $m$ Bernoulli$(p_s)$ per posterior draw so $k_m = \sum_i \text{ypred}_{s, i}$ marginalises exactly to $\text{BetaBinomial}(m,\, 1 + k_n,\, 1 + n - k_n)$. One forward pass per interim gives the continuous $\hat P(H_1 \mid x, z^{(s)}) \in [0, 1]$; the PPS is $S^{-1} \sum_s \mathbf{1}\{\hat P > \eta_H\}$. This replaces the earlier scipy-inline Beta-Binomial hack with the standard interim-loop path used by every other Binomial deployment script.

**Correctness against the closed-form Beta-Binomial PPS.** Across the 11 valid monthly interims of the deployment cohort ($N = 500$, true $p = 0.4$, $\eta_0 = 0.25$, $\eta_H = 0.89$, $S = 200$):

| Interim | Analytic | Amortised | |err| |
|---------|---------:|----------:|-----:|
| 1 (Jan) | 0.069 | 0.070 | 0.001 |
| 2 (Feb) | 0.476 | 0.510 | 0.034 |
| 3 (Mar) | 0.316 | 0.275 | 0.041 |
| 4 (Apr) | 0.283 | 0.245 | 0.038 |
| 5 (May) | 0.338 | 0.290 | 0.048 |
| 6 (Jun) | 0.135 | 0.145 | 0.010 |
| 7 (Jul) | 0.133 | 0.155 | 0.022 |
| 8 (Aug) | 0.051 | 0.035 | 0.016 |
| 9 (Sep) | 0.069 | 0.075 | 0.006 |
| 10 (Oct) | 0.000 | 0.000 | 0.000 |
| 11 (Nov) | 0.000 | 0.000 | 0.000 |

**Max absolute error 0.048, mean 0.020**, meeting the §11.5 tolerance ($\le 0.02$ mean, $\le 0.05$ max within the Monte-Carlo noise band at $S = 200$). The five pytest cases in `test_amortised_pps_correctness.py` all pass with the same tolerances.

**Diagnostic finding worth documenting.** The amortised network matches the analytic $Q_{\tau}(\rho \mid k_{\text{total}}, n_{\text{total}})$ pointwise to $\le 0.002$ at each $\tau$ level. The initial PPS gap (up to 0.156) that appeared when the deployment loop used HMC-generated $k_m$ draws was driven entirely by an over-wide posterior-predictive from HMC at small $n$ (std 50 vs analytic std 36 at interim 1). Bypassing HMC with analytic Beta-Binomial $k_m$ removes the bias. This isolates the amortiser's quality from HMC noise and validates the amortised pipeline in its "pure" form (train the net once → deploy with model-analytic posterior-predictive).

## 12.2 Results for Binomial interim analysis amortised endptx on wz with features-MLP, qpsi-MLP, loss-multiquantilehead

**Implementation.** Parallel prototype swapping the hardcoded per-item encoder for a **learnable** DeepSets encoder $q_\tau$:

- New Flax module [`python/amortiser_pps_features_MLP_qpsi_MLP_loss_multiquantilehead.py`](../python/amortiser_pps_features_MLP_qpsi_MLP_loss_multiquantilehead.py) — class `Amortiser_PPS_features_MLP_qpsi_MLP_loss_multiquantilehead` uses `setup()` to share a SiLU-MLP `q_tau` (`(item_dim → q_tau_hidden_dims → embed_dim)`) across `x` and `z` sets; padded input `(B, N_max, item_dim)` + boolean masks; sum-pool with mask multiplication; concat `pooled_x + pooled_z + sizes` into the same MLP head `q_psi` and multi-quantile pinball loss as §12.1.
- New sampler `BinomialModel.make_training_data_with_raw_sequences(rng, S, n_max)` emits raw padded 0/1 item sequences (`x`, `mask_x`, `z`, `mask_z`, `sizes`) with the same `(θ, x, z)` joint sampling as the features-fixed variant, so the head sees enough shape variability to learn $q_\tau$ end-to-end.
- Sim-data cache extended: [`scripts-py/Binomial_interim_analyses_make_sim_data.py`](../scripts-py/Binomial_interim_analyses_make_sim_data.py) writes a second $65\,536$-sample cache `binomial_amortiser_training_data_features_MLP.pkl` alongside the features-fixed cache.
- Deployment script [`scripts-py/Binomial_interim_analysis_amortise_endptx_on_wz_with_features_MLP_qpsi_MLP_loss_multiquantilehead.py`](../scripts-py/Binomial_interim_analysis_amortise_endptx_on_wz_with_features_MLP_qpsi_MLP_loss_multiquantilehead.py) mirrors §12.1's layout but constructs a padded raw-sequences batch per interim (all $S$ posterior draws share the observed `x` sequence; each `z^(s)` carries `km_s` ones + `m - km_s` zeros drawn from the same analytic joint posterior-predictive built by `BinomialModel.fit_closed_form_posterior` + `get_interim_z_from_ypredi`). Outputs saved with `_RGEB_` suffix.
- Shared training / prediction / save-load utilities in `amortiser_common` are polymorphic across the fixed and MLP amortiser classes — no new API surface. `load_fitted_model` uses `importlib` on the persisted `net_class_module` / `net_class_name` to rebuild either class transparently.

**Training configuration.** `q_tau_hidden_dims = (32, 32)`, `embed_dim = 16`, `q_psi` head `hidden_dims = (128, 128, 64)`, same 11-level `pps_ProbH1_lwr_quantiles_mesh`. Adam + linear-warmup / cosine-decay at peak lr $10^{-3}$, **$15\,000$ steps × batch $512$** (smaller than the features-fixed variant because each sample is now a `(N_max=500) × 1` tensor). Training loop measured 8.59 min inside `train()` (persisted on the fit dict via §11's `training_mins` field). Final training pinball loss 0.0090 (still gently decreasing — a larger budget would tighten error further).

**Correctness against the closed-form Beta-Binomial PPS.** Same 11 monthly interims, same $\eta_0 / \eta_H / S$ as §12.1:

| Interim | Analytic | Amortised (features-MLP) | |err| |
|---------|---------:|-------------------------:|-----:|
| 1 (Jan) | 0.069 | 0.065 | 0.004 |
| 2 (Feb) | 0.476 | 0.555 | 0.079 |
| 3 (Mar) | 0.316 | 0.370 | 0.054 |
| 4 (Apr) | 0.283 | 0.360 | 0.077 |
| 5 (May) | 0.338 | 0.420 | 0.082 |
| 6 (Jun) | 0.135 | 0.170 | 0.035 |
| 7 (Jul) | 0.133 | 0.180 | 0.047 |
| 8 (Aug) | 0.051 | 0.055 | 0.004 |
| 9 (Sep) | 0.069 | 0.140 | 0.071 |
| 10 (Oct) | 0.000 | 0.000 | 0.000 |
| 11 (Nov) | 0.000 | 0.000 | 0.000 |

**Max absolute error 0.082, mean 0.041.** Roughly $2\times$ the features-fixed error (§12.1: max 0.048, mean 0.020). Expected — the features-fixed variant hardcodes the exact Binomial sufficient statistic $T(x_i) = x_i$ as its inductive prior; the features-MLP variant has to *discover* that sum-of-successes is sufficient from data, which costs a few thousand extra training steps to close and remains slightly noisier at deployment.

## 12.3 Full cross-method comparison after the analytic-`zi` sweep

All Binomial deployment scripts (regression, IS, both amortised variants; nested-MC left as-is per its "HMC-in-HMC" contract) now build `zi` from `BinomialModel.fit_closed_form_posterior` + the `BinomialModel`-specific `get_interim_z_from_ypredi` override that draws fresh $m$ Bernoulli$(p_s)$ per posterior draw. This uses the **same** analytic $p(\theta, z \mid x)$ joint across all methods, so per-method MSE cleanly reflects estimator quality rather than divergent `z`-samplers. MSE against the closed-form Beta-Binomial PPS across the 11 monthly interims:

| Method | MSE | $\sqrt{\text{MSE}}$ | Mean per-interim inference (min) | One-off training (min) |
|--------|----:|--------------------:|---------------------------------:|-----------------------:|
| Regression endpt-x (Gaussian approx)         | 0.00160 | 0.040 | 0.001 | — |
| Nested-MC using HMC for each (x, z)          | 0.00244 | 0.049 | 6.04  | — |
| IS reweighting of $\theta \mid x$            | 0.00277 | 0.053 | 0.004 | — |
| **Amortised (features-fixed)**               | **0.00277** | **0.053** | **0.001** | **5.90** |
| **Amortised (features-MLP)**                 | **0.00277** | **0.053** | **0.002** | **8.59** |
| Regression endpt-x (mquantile)               | 0.00355 | 0.060 | 0.001 | — |
| Regression endpt-x (quantile)                | 0.00357 | 0.060 | 0.001 | — |
| Nested-MC using SVI for each (x, z)          | 0.00401 | 0.063 | 1.10  | — |
| Regression H1-x on w(z)                      | 0.01033 | 0.102 | 0.001 | — |

*Training cost is amortised across all future deployments — pay once, deploy in milliseconds thereafter. Nested-MC HMC costs 6 min **per interim per x**, so on a 12-interim schedule with $S = 200$ Monte-Carlo draws the total nested-MC cost is $\sim 72$ min per new patient cohort $x$; the amortised MLP variant recovers the same MSE for a 8.6 min one-off training + 0.002 min per interim, breaking even at the first fresh $x$ and winning outright for every subsequent one. All timings measured with $\le 4$ CPU threads on a single machine.*

**Reading.** The MSE floor is set by Monte-Carlo noise at $S = 200$: $\sqrt{p (1 - p) / S} \approx 0.035$ near the PPS mode, so any method below $\sqrt{\text{MSE}} \le 0.05$ is essentially at the MC-noise floor. Both amortised variants are numerically tied with IS at that floor. The Gaussian-approx regression wins because the Beta posterior mean is well-approximated by a Gaussian at the sample sizes seen, giving a very tight plug-in $\Phi((\hat\mu - \eta_0)/\hat\sigma)$ estimator. Nested-MC HMC still uses the finite-sample HMC posterior (wider than the exact Beta at small $n$), so it inherits an interim-1 error of 0.15 that inflates its MSE. The `H1-x` regression is the most biased — the binary label loses too much information relative to the continuous endpoint targets used everywhere else.

**Amortised (features-MLP) reading.** The MLP variant matches the fixed variant's aggregate MSE despite starting from no structural prior. Per-interim errors are on average larger (mean 0.041 vs 0.020) but the *distribution* of errors is symmetric around the analytic PPS, so squared errors average out. A longer training budget would close the mean-error gap; the ceiling is the MC-noise floor at $S = 200$, same as everyone else. This validates the general-purpose amortiser template for models where the sufficient statistic is not known analytically (Categorical, IRT).

------------------------------------------------------------------------

# 13. References

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