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

## scheduler
title: Scheduler over pretrained any-J any-N decision-making amortisers
line: Each manifest specifies a trained amortiser for a specific decision target type $d$. Target types are defined by functional used for group comparison, directionality of effects, effect size warping. The scheduler directs deployment of the specified amortiser(s). Each amortiser is applicable to any-J and any-N.
card: $\rho^{(1)}$  pre-post rel improvement
cardline: net=scale-feat-paired · warp=none
cardline: bvm_shrink=Stein
card: $\rho^{(2)}$  2-arm signed difference
cardline: net=groupdiff · warp=none
card: $\rho^{(3)}$  2-arm rel improvement
cardline: net=widetok-spr · warp=log2/logit
card: $\cdots$
cardline: 
hover: The scheduler reads the yaml-manifest, and directs the deployment of pretrained amortisers. The manifest declares item types, effect size types and decision endpoints, and maps library functions to generate tokens, compute effect sizes, and fetch pretrained amortisers. Additional manifest elements are grouping instructions, directionality in item levels, effect size thresholds $\eta_0$. 


## manifest
title: Item/Effect size specifications
line: Yaml-type manifest with descriptors for the $j=1,\dotsc,J$ items, and $r=1,\dotsc,R$ effect sizes, denoted by $x^{\text{manifest}}$
hover: Auxiliary information for each response item $j=1,\dotsc,J$ to facilitate any-item amortisation, including for each iten the number of categories, the type of the effect size to build appropriate tokens, and the direction of improvement in the signed effect size. Auxiliary information for each effect size $\rho_r$, $r=1,\dotsc,R$, including library function for item-level group comparison (pre/post relative change, pre/post difference, intervention/control relative change, etc), effect size warps, scientifically meaningfully levels $\eta_0$. Effect size descriptors specify the amortiser that is deployed. Each effect-specific amortiser accepts arbitrary numbers of items and arbitrary numbers of ordered categories for each item. Item descriptors feed into the effect-specific amortiser to enable deployment to different item panels.


## itemattn
title: Item-level self-attention
line: Generate $s=1,\dotsc,S$ future data sets $z^s$ of $m_t = N-n_t$ participants, tokenize $(y_{1:n_t},z^s)$ by participant and item, embed to $e^{\text{in}}_{ijs}\in\mathbb{R}^{32}$, refine embeddings to $e^{\text{att}}_{ijs}\in\mathbb{R}^{32}$ with $J\times J$ self-attention. $(n_t + m_tS)J$ tokens, $(n_t + m_t)SJ$ embeddings.

Prepare/scale raw and predicted future data samples into tokens $\text{tok}_{1:N}$ to facilitate learning one generic embedding into $\mathbb{R}^E$, add helpers to memorize data structure. 
hover: Concatenate participant- and item-specific raw inputs from
interim data $y_{1:n_t}$ and $s=1,\dotsc,S$ posterior predicted future data samples $z^s_{1:m_t}$; scale; add index of participants, interim/future, groups, metadata. Here, $m_t = N - n_t$ and future data samples $$ are generated from the posterior predictive distribution under the PCM model fit to the interim data. Tokens are embedded to $e^{\text{in}}_{ijs}\in\mathbb{R}^{32}$ with an input MLP $q_{\text{input}}$. If items carry cross-correlation, we refine encodings to $e^{\text{att}}_{ijs}\in\mathbb{R}^32$ with $J\times J$ item-level self-attention, batched over participants $i$ and futures $s$, injecting invariance to item ordering and equivariance to participant and futures ordering.

## deepset
title: Deep-set feature encoder
line: Mean-pool participant embeddings $e^{\text{att}}_{ijs}$ to latent features $h^{\text{y}}_{j}= n^{-1}_t \sum_{i=1}^{n_t} e^{\text{att}}_{ij1}$, $h^{\text{z}}_{js}= m_t^{-1} \sum_{i=n_t+1}^{N} e^{\text{att}}_{ijs}$ from current and future data. $2SJ$ Features.
hover: Mean pool embeddings over current and future participants, $h^{\text{y}}_{j}= n^{-1}_t \sum_{i=1}^{n_t} e^{\text{att}}_{ij1}$, $h^{\text{z}}_{js}= m_t^{-1} \sum_{i=n_t+1}^{N} e^{\text{att}}_{ijs}$, injecting invariance to participant ordering. 

## featureattn
title: Feature self-attention
line: Concatenate and refine latent features $h_{js} = [h^{\text{y}}_{j}, h^{\text{z}}_{js}]$ with $JxJ$ self-attention, batched over futures $s$, to $h_{js}^{\text{att}} \in \mathbb{R}^{32}$. $SJ$ attended features.
hover: Transform $h_{js}\in \mathbb{R}^{64}$ to keys, queries and values in $\mathbb{R}^{32} with three MLPs. Refine latent features with $J\times J$ item-level self-attention, maintaining invariance to item ordering and equivariance to futures ordering. 

## posterior
title: Posterior effect sizes on interim data (Inference target)
line: Partial credit model fitted on interim outcome data with SVI, providing numerical samples of the joint posterior of $R$ effect sizes $\rho_{1:R}$ given interim data, $p(\rho_{1:R} | y_{1:n_{t}})$.
hover: The partial credit model fitted with SVI on the interim data. Supplies the target draws $\rho_r^(s)$, for each of $r=1,\dotsc,R$ effects relevant for decision making. Effect sizes are deterministic functions of model parameters $\rho=f(\theta)$. The number of effects $R$ for decision making is much smaller than the number model parameters. 

## head
title: Distributional approximation of future effect sizes
line: MLP $g_{\tau^{\text{effects}}}$: From tokens and attended features generate mesh of $q=1,\dotsc,Q$ quantiles of the future marginals $p(\rho_j | y_{1:n_t}, z^s)$ for each item $j$ and future $s$, denoted by $\hat{\rho}_{jsq}$. 
hover: To capture uncertainty in future effect sizes for each item, we train a mesh of monotone quantiles of $\rho_j | y_{1:n_t}, z^s$. 

## train
title: Train
line: Train on prior effect sizes $\rho^s$ and prior-predictive data $y^s_{1:n_t}, z^s \sim \rho^s$, any-J and any-N, ahead of outcome data. Trained with pinball loss over all quantiles $q=1,\dotsc,Q$ against the conditional view $\rho^s|y^s_{1:n_t}, z^s$.
hover: The whole amortiser including the item/feature self-attention encoders and the quantile head MLP $g_{\tau^{\text{effects}}}$ are trained end-to-end on prior-predictive simulations (parameters drawn from the PCM prior; interim and future cohorts simulated at any J and any N), so we obtain a foundational decision-making amortiser. We use a standardised version of the pinball loss function $\ell(\tau) = \sum_{s=1}^{S^\text{train}} \sum_{j=1}^J \sum_{q=1}^Q L^{\text{pinball}}_q(\rho_{js} , \hat{\rho}_{jsq})$, where $\rho_{js}$ are from the prior $p(\rho)$ used to generate the train-time interim data $y^{\text{train}}_{1:n_t}$ and the train-time future data $z^s$. $\hat{\rho}_{jsq}$ is the $q$th output of the MLP $g_{\tau^{\text{effects}}}$. Thus we train estimated $\hat{\rho}_{jsq}$ to match the $q$th quantile of the $Q$ dimensional approximation of $p(\rho_j | y_{1:n_t}, z^s)$. Training is over fresh mini-batches $s$. 


## finetune
title: Fine-tune
line: Use the interim history to refit only the MLP $g_{\tau^{\text{effects}}}$ with all upstream components frozen. Fine tuned with the same pinball loss function, using posterior effect sizes at interims $\rho^s \sim p(\rho^s|y_{1:n_t^'})$ for all $t^'=1,\dotsc,t$ in the conditional view $\rho^s|y_{1:n_t^'}, z^s$.
hover: At deployment, the MLP's parameters $\tau^{\text{effects}}$ are refit with the pinball loss function $\ell(\tau^{\text{effects}}) = \sum_{t'=1}^t \sum_{s=1}^{S} \sum_{j=1}^J \sum_{q=1}^Q L^{\text{pinball}}_q(\rho_{t^'js} , \hat{\rho}_{t^'jsq})$, where the $\rho_{t^'j,s}$ are from the $t^'$ posterior interim $p(\rho_j | y_{1:n_{t^'}})$ used to generate the future data $z^s$, and $\hat{\rho}_{t^'jsq}$ is the estimated $q$th quantile of the $Q$ dimensional approximation of $p(\rho_j | y_{1:n_t}, z^s)$.

## stein
title: von-Mises calibration
line: Estimate shifts and posterior contraction of $p(\rho_j | y_{1:n_{t^'}})$ for $t^'=1,\dotsc,t$, and use to calibrate $p(\rho_j | y_{1:n_{t^'}}, z^s)$.
hover: At deployment, a von-Mises powerlaw is fit to the marginal posteriors $p(\rho_j | y_{1:n_{t^'}})$ over the interim history $t^'=1,\dotsc,t$ for each $j$. The variance contributions in $p(\rho_j | y_{1:n_{t^'}}, z^s)$ are decomposed into between and within variance terms, and the within-variance term that is inflated due to Monte Carlo error in small future data sets is adjusted to match the von-Mises powerlaw. Affine shifts are also adjusted.

## pps
title: Amortised predictive probability of success (Decision-making target)
line: Predicted probability that effect will have scientifically meaningful size
line: in the total of $N$ participants, over future response data.
line: Central quantity for interim decision-making.
hover: Scientifically meaningful effects are defined as the posterior distribution that the $r$th effect size falls into $H^(1)_r :=\{ \rho_r > \eta^{(1)}_r\}$. The scientifically meaningful effect size is denoted by $\eta^{(1)}_r$. Here, $\rho_r$ are from the posterior distribution over the interim outcome data and an instance of unseen future data, $p(\rho | y_{1:n_{t}}, z^s_{1:m_{t}})$, where $m_t$ are the number of unseen participants, $m_t = N - n_t$.
