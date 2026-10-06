<!--
architecture_figure_1_text.md — EDIT THE BOX TEXT HERE, then regenerate with:
    pixi run python paper_adm/architecture_figure_1.py
One `## <box-id>` section per box. Fields (repeat `line:` for each body line):
    title:     bold heading (one line)
    subtitle:  optional second heading line (smaller)
    line:      one body line (repeatable, in order)
    hover:     tooltip shown on hover in the HTML version (one line)
Inline maths goes in $...$ (TeX): renders as mathtext in the PDF and as unicode in the HTML.
The box POSITIONS, arrows and the manifest ρ-cards live in architecture_figure_1.py (LAYOUT / CARDS).
-->

## interim_data
title: Interim outcome data
line: $J$ ordered categorical responses of $n_t$ participants
line: up to interim time $t$ out of $N$ scheduled total participants, 
line: $y_{1:n_{t}}$. Responses in groups $g$ (arms, pre/post) for 
line: effectiveness comparisons.
hover: Per participant $i$, $i=1,\dotsc,n_t$, responses $y_{i,j}$  for $j=1,\dotsc,J$ items are recorded. Responses are ordered categories of $K_j$ graded levels. Responses are structured into groups $g$ that may correspond to responses before and after an intervention for the same individual, or represent individuals in placebo or intervention arms.

## manifest
title: Pretrained federated any-J any-N decision-making amortisers
subtitle: one endpoint per $\rho$ · routed to a trained net + warp + $\eta_0$ grid · federated over ONE joint SVI fit
hover: A list of endpoints; each rho declares its trained net, warp, eta0 grid and build strategy. federated_deploy routes each to the net registry and federates them over one joint SVI fit — a single joint decision across several rho.

## item_metadata
title: Item metadata 
line: Descriptors for the $j=1,\dotsc,J$ items, $x^{\text{meta}}$ to 
line: make the generic model deployable to applications with different
line: item panels.
hover: Auxiliary information for each response item $j=1,\dotsc,J$ to facilitate any-item amortisation, including for each iten the number of categories, the type of the effect size to build appropriate tokens, and the direction of improvement in the signed effect size.


## tokeniser
title: Tokeniser
line: Prepare/scale raw and predicted future data samples into tokens $\text{tok}_{1:N}$ to facilitate learning one generic embedding into $\mathbb{R}^E$, add helpers to memorize data structure. 
hover: Concatenate participant- and item-specific raw inputs from
interim data $y_{1:n_t}$ and posterior predicted future data samples $z^s_{1:m_t}$; scale data; add index of participants, interim/future, groups; and add metadata. Here, $m_t = N - n_t$ and $s=1,\dotsc,S$ future data samples are generated from the posterior predictive distribution under the PCM model fit to the interim data.

## deepset
title: Ragged deep-set encoder
line: Standard SBI mean pool over encoded participant features in $E$-dim 
line: embedding space, separately for each item $j$, group $g$, and 
line: current/future data $b$. 
hover: Mean pool over a ragged batch of participant embeddings with no padding or masks, $\text{deepset}^s_{j,g,b} = 1/n_b \sum_{i \in g,b} q_\tau( \text{tok}^s_i )$. Mean-pooling of encoded features injects symmetry over participants. Separate mean-pools by $j,g,b$ ensure that the $E$-dimnesional encoded features can borrow information across related items, while structure over groups and interim/future is retained.

## xattn
title: Cross-component attention
line: item query $q_j(M_j)$ over components
line: deepsetXcompAtt · one net, all j
hover: Each item's query, built from its metadata, attends across the pooled components (deepsetXcompAtt). One network serves every item.

## posterior
title: Posterior effect sizes on interim data (Inference target)
line: Partial credit model fitted on interim outcome data
line: with SVI, providing numerical samples of the joint posterior
line: of $R$ effect sizes $\rho_{1:R}$ given interim
line: data, $p(\rho_{1:R} | y_{1:n_{t}})$.
hover: The partial credit model fitted with SVI on the interim data. Supplies the target draws $\rho_r^(s)$, for each of $r=1,\dotsc,R$ effects relevant for decision making. Effect sizes are deterministic functions of model parameters $\rho=f(\theta)$. The number of effects $R$ for decision making is much smaller than the number model parameters. 

## head
title: Multi-quantile head
subtitle: generating approximation to effect size posterior over interim and future data
line: 5 quantiles $\tau=\{.05,.25,.5,.75,.95\}$ of $\rho_j\,|\,x^{(t)},z^{(s)}$
line: $\times$ global $\sigma$ · warp$^{-1}$ · head-ft (expand)
hover: Monotone standardised quantiles of rho_j given x^(t) and z^(s), un-standardised by the global sigma and the inverse warp. The head is the one gradient-fine-tuned block at deploy (head-ft, expand).

## stein
title: Stein–von Mises calibration
subtitle: to match posterior contraction seen over interims
line: fit $SD_j=C_j\,n^{-p_j}$ · James–Stein shrink $p_j\to\bar p$
line: reshape var: pin at n, at $n+m=N_{ref}$
hover: Deploy-time calibration (not a net fine-tune): fit the per-item power law to the SVI SD targets, James-Stein shrink the exponent toward the pooled median, then reshape each interim's conditional so the mixture variance matches the law at n, pinned at n+m=N_ref.

## pps
title: Amortised predictive probability of success (Decision-making target)
line: Predicted probability that effect will have scientifically meaningful size
line: in the total of $N$ participants, over future response data.
line: Central quantity for interim decision-making.
hover: Scientifically meaningful effects are defined as the posterior distribution that the $r$th effect size falls into $H^(1)_r :=\{ \rho_r > \eta^{(1)}_r\}$. The scientifically meaningful effect size is denoted by $\eta^{(1)}_r$. Here, $\rho_r$ are from the posterior distribution over the interim outcome data and an instance of unseen future data, $p(\rho | y_{1:n_{t}}, z^s_{1:m_{t}})$, where $m_t$ are the number of unseen participants, $m_t = N - n_t$.
