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
H_1 : p > p_0.
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

# 3. Analytically tractable examples

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
\begin{align*}
p(x_i \mid p) & 
= p^{x_i}(1-p)^{1-x_i} = \exp\{x_i \log p + (1-x_i)\log(1-p)\} 
\\
& 
= \exp\{\log\left(\frac{p}{1-p}\right) x_i + \log(1-p)\},
\end{align*}
with natural parameter \(\eta(p) = \log(p/(1-p))\), sufficient statistics \(T(x_i) = x_i\), and log-partition function \(A(p) = -\log(1-p)\). For \(n\) observations, the likelihood depends only on the sufficient statistic \(T(x_{1:n}) = \sum_{i=1}^n x_i = k^n\), the total number of successes. After observing \(x_{1:n}\) data points with $k^n$ successes the posterior is

\[
p \mid x \sim \text{Beta}(a+k^n, b+n-k^n).
\]


Suppose the future, remaining data \(z\) comprise \(m\) additional data points. Then the posterior predictive distribution of $k^m$ future successes is

\begin{align*}
k^m & \sim \int \text{Binomial}(m,p) \text{Beta}(p; a+k^n, b+n-k^n) \\
& = \text{Beta-Binomial}(m,a+k^n,b+n-k^n).
\end{align*}

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
\begin{align*}
p(x_i \mid p) & 
= \prod_{k=1}^K p_k^{1_{x_i = k}} 
= \exp\left\{\sum_{k=1}^K 1_{x_i = k} \log p_k\right\}
\\
&
= \exp\left\{\sum_{k=1}^{K-1} \log\left(\frac{p_k}{p_K}\right) 1_{x_i = k} + \log p_K\right\},
\end{align*}
with natural parameters \(\eta_k(p) = \log(p_k/p_K)\) for \(k = 1, \ldots, K-1\), sufficient statistics \(T_k(x_i) = 1_{x_i = k}\) for \(k = 1, \ldots, K-1\), and log-partition function \(A(p) = -\log p_K\). For \(n\) observations, the likelihood depends only on the sufficient statistic \(T(x_{1:n}) = (k_1^n, \ldots, k_K^n)\), the count vector where \(k_k^n = \sum_{i=1}^n 1_{x_i = k}\) is the number of observations in category \(k\). After observing \(x_{1:n}\) data points with counts \(k_k^n\) the posterior is

\[
p \mid x \sim \text{Dirichlet}(\alpha_1 + k_1^n, \ldots, \alpha_K + k_K^n).
\]

Suppose the future, remaining data \(z\) comprise \(m\) additional data points. Then the posterior predictive distribution of the future count vector \((k_1^m, \ldots, k_K^m)\) is

\begin{align*}
(k_1^m, \ldots, k_K^m) & \sim \int \text{Multinomial}(m, p) \, \text{Dirichlet}(p; \alpha_1 + k_1^n, \ldots, \alpha_K + k_K^n) \, dp \\
& = \text{Dirichlet-Multinomial}(m, \alpha_1 + k_1^n, \ldots, \alpha_K + k_K^n).
\end{align*}

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
\begin{align*}
& 
\mathbb{E}_{H \mid x}[L(a_1, H)] < \mathbb{E}_{H \mid x}[L(a_0, H)] \\
\Leftrightarrow &
L_{FP} \cdot (1 - P(H_1 \mid x)) < L_{FN} \cdot P(H_1 \mid x) \\
\Leftrightarrow &
P(H_1 \mid x) > L_{FP} / ( L_{FP} + L_{FN} ).
\end{align*}


This shows that the Bayes optimal decision threshold \(\eta\) can be expressed in terms of utilities or loss terms. In particular, if false positives are more costly, \(\eta > 0.5\) and if false negatives are more costly, \(\eta < 0.5\). This provides a utility-based foundation for choosing \(\eta\) in contrast to ad-hoc 0.89 or 0.95.

## 5.2 Net benefit 

A related concept is expected utility, normalised such that 1 unit corresponds to 1 true positive. Under the above action/partition set matrix, the relative cost of one false positive is \(w\) and true pos and true neg have utility/cost 0. This does not depend on the particular units of the losses, and so is more easily interpretable.

Under these utilities, the net benefit of declaring success is given by

\begin{align*}
NB(x) 
& 
= P(H_1 \mid x) - w \cdot P(H_0 \mid x) 
\\
&
= P(H_1 \mid x) - w \cdot (1 - P(H_1 \mid x) )
\end{align*}

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

# 6. Amortized inference workflow

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
\begin{align*}
q_\phi(x,z) 
& 
\approx \mathbb{E}_{\theta|x,z}(L(A,\theta)) 
\\
& 
=  \int L(A,\theta) p(\theta | x,z) d\theta,
\end{align*}
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

# 7. Workflow Step 1: generating labelled training data

Generating training data is more challenging than in neural posterior estimation (NPE). We consider two distinct cases depending on whether the current data \(x\) is fixed or variable.

### Case A: Fixed \(x\) with access to the posterior

In this simpler setting, we have observed a specific interim dataset \(x_{1:n}\) and have access to the posterior distribution \(p(\theta \mid x)\). Then it is relatively straightforward to generate labelled training data:

1. Draw posterior samples: For simulation budget \(s = 1, \dots, S\), sample \( \theta^{(s)} \sim p(\theta \mid x) \).

2. Generate future data: Sample future observations from the likelihood \( z^{(s)} \sim p(\cdot \mid \theta^{(s)}) \).

3. Compute labels: For each simulated future dataset, we need to compute the posterior success probability
   \[
   y^{(s)} := P(H_1 \mid x, z^{(s)}) = \int 1_{\theta \in H_1} \, p(\theta \mid x, z^{(s)}) \, d\theta.
   \]
   Since we have the posterior \(p(\theta \mid x)\), we can use importance reweighting:
   \begin{align*}
   w^{(s)}_k & \propto p(z^{(s)} \mid \theta_k), \quad \theta_k \sim p(\cdot \mid x) \\
   y^{(s)} & \approx \frac{\sum_k w^{(s)}_k \, 1_{\theta_k \in H_1}}{\sum_k w^{(s)}_k}.
   \end{align*}

### Case B: Variable \(x\)

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
    \begin{align*}
    w_k & = \frac{p(\theta_k \mid x^{(s)}, z^{(s)})}{q(\theta_k \mid \theta^{(s)})} \\
    & \propto \frac{p(\theta_k) \, p(x^{(s)} \mid \theta_k) \, p(z^{(s)} \mid \theta_k)}{q(\theta_k \mid \theta^{(s)})},
    \end{align*}
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
      \begin{align*}
      w_{s,k}^{(x)} & \propto w_{s,k}^{(0)} \cdot p(x^{(s)} \mid \theta_{s,k}) \\
      \tilde{w}_{s,k}^{(x)} & = w_{s,k}^{(x)} / \sum_{j=1}^K w_{s,j}^{(x)}
      \end{align*}
      The normalized weights \(\{\theta_{s,k}, \tilde{w}_{s,k}^{(x)}\}_{k=1}^K\) now represent \(p(\theta \mid x^{(s)})\). Resample particles if the effective sample size \(1/\sum_k (\tilde{w}_{s,k}^{(x)})^2\) is too small.
   
   5. Generate future data from posterior predictive: Sample one future dataset from the weighted particle system:
      \[
      z^{(s)} \sim \sum_{k=1}^K \tilde{w}_{s,k}^{(x)} \, p(\cdot \mid \theta_{s,k}).
      \]
      In practice, select a particle index \(\kappa \sim \text{Categorical}(\tilde{w}_{s,1}^{(x)}, \ldots, \tilde{w}_{s,K}^{(x)})\) and generate \(z^{(s)} \sim p(\cdot \mid \theta_{s,\kappa})\).
   
   6. Sequential update for \(z^{(s)}\): Update the weights again to approximate \(p(\theta \mid x^{(s)}, z^{(s)})\) by evaluating how well each particle explains the joint data:
      \begin{align*}
      w_{s,k}^{(x,z)} & \propto \tilde{w}_{s,k}^{(x)} \cdot p(z^{(s)} \mid \theta_{s,k}) \\
      \tilde{w}_{s,k}^{(x,z)} & = w_{s,k}^{(x,z)} / \sum_{j=1}^K w_{s,j}^{(x,z)}
      \end{align*}
   
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


# 8. Step 2: learning neural architectures

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

# 9. Step 3: Deployment, amortising PPS 

After training, 

* for any observed data \(x_{1:n}\) we can compute \( \sum_{i=1}^n q_\tau(x^{(s)}_i) \). 
* give interim data \(x_{1:n}\), in practice we will estimate today's posterior \(p(\theta | x)\) with some Monte Carlo algorithm, and compute the first target objective, have we already won, \( 1\{ P(H_1 \mid x ) > \eta \). It is cheap to generate posterior predictions \(z^\star_{1:m}\). Then we can compute \( \sum_{i^\star=1}^m q_\tau(z^{(s)}_{i^\star}) \) for the corresponding \(m = N - n\).
* Next we can compute the PPS approximation

\[
q_\psi\bigg(\sum_{i=1}^n q_\tau(x^{(s)}_i) + \sum_{i^* =1}^m q_\tau(z^{(s)}_{i^*})\bigg).
\]

# 10. Tasks

* Generate functions to compute the exact \( PPS(x) \) in the case of simple Binomial and Categorial models
* Generate a case study that illustrates contracting posterior distributions with additional data, and evolving \( PPS(x) \) with additional data under the simple Binomial model. Illustrate stopping early for futility/safety and for efficacy. 
* Generate functions that implement Case A and B in generating labelled training data.
* Test these functions against the analytical formula under the simple Binomial and Categorial model.
* Test these functions against the analytical formula under the Categorial model.