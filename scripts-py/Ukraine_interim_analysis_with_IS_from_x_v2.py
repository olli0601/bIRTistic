#!/usr/bin/env python3
"""
Ukraine interim analysis v2: diagnosing + mitigating importance-sampling weight
degeneracy at an EARLY interim (interim_id = 1 only).

Background
----------
``Ukraine_interim_analysis_with_IS_from_x.py`` reuses the x-posterior draws
``theta_k ~ p(theta | x)`` and reweights them by ``p(z_s | theta_k)`` (Case A,
fixed x). At early interims the future block z adds a large number of
participants (interim 1: m = N_full - n_x), so the target p(theta | x, z) sits
far from the proposal p(theta | x). The self-normalised weights then collapse
onto a single draw (ESS/N -> 1/N), as the perf plot showed.

This script experiments on interim 1 only:

  EXP A  ESS/N as a function of the number of assimilated z-participants
         (cumulative). Shows *where* the proposal/target mismatch breaks IS.
  EXP B  Pareto-smoothed importance sampling (PSIS, Vehtari et al.) on the full
         z block: report the Pareto k-hat tail diagnostic and the ESS/N before
         vs after smoothing.

Why only these two (no tempering / SMC here): in pure Case-A reweighting the
particle *locations* are frozen at the p(theta|x) draws, so no reweighting
scheme (plain, tempered, or PSIS) can manufacture effective sample size beyond
what those fixed draws support. PSIS still helps by taming the tail and giving
an honest reliability diagnostic (k-hat); the principled cure (SMC with
resample-move, or simply the MC refit path) needs to *move* particles and is
out of scope for this no-refit experiment.

Usage:
    cd /Users/or105/git/bIRTistic
    pixi run python scripts-py/Ukraine_interim_analysis_with_IS_from_x_v2.py
"""

# %%

import os
import sys
import time
from pathlib import Path

try:
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
except NameError:
    project_root = Path.cwd()
    if project_root.name == 'scripts-py':
        project_root = project_root.parent

python_path = str(project_root / 'python')
if python_path not in sys.path:
    sys.path.insert(0, python_path)

import numpy as np
import pandas as pd
import arviz as az
from scipy.special import softmax
from plotnine import (
    ggplot, aes, geom_line, geom_hline, geom_point, facet_wrap,
    scale_y_log10, scale_color_manual, theme_bw, theme, element_text, labs,
)
import jax
import jax.numpy as jnp
from jax import random
import jax.scipy.stats as jstats
import numpyro
import numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS
import warnings
warnings.filterwarnings('ignore')

from data_loading import read_data_ukraine
from fit_partial_credit_model import (
    fit_partial_credit_model_ncats_pyrosvi,
    _fit_partial_credit_make_stan_data,
    eval_loglik_partial_credit_model_ncats,
)
from fit_interim import get_interim_x, get_interim_z_from_ypredi
from utils import _futurama_palette

print("✓ Imports successful")

# %%

# =============================================================================
# Configuration
# =============================================================================

np.random.seed(42)
seed = 123

dir_data = "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data"
file_data = os.path.join(dir_data, "Ukraine_Hope_Groups_Baseline_Endline_Wide_Aug6.csv")

# Reuse the existing IS run directory so the resume=True fit re-loads the cached
# interim-1 x-posterior instead of re-fitting from scratch.
dir_out = "/Users/or105/sandbox/bIRTistic/py-ukraine-interim-with-importance-sampling-260526"
os.makedirs(dir_out, exist_ok=True)

file_prefix = "pcm_1_interim"
svi_algorithm = 'AutoLowRankMultivariateNormal'

pps_z_total = 10
INTERIM_ID = 1
x_formula = "~ time - 1"
categorical_threshold = 2

print(f"Output dir: {dir_out}")

# %%

# =============================================================================
# Load + preprocess full Ukraine data (same pipeline as the main IS script)
# =============================================================================

print("\nLoading Ukraine data...")
raw = read_data_ukraine(file_data)
dp = raw['dp'].copy()
dit = raw['dit'].copy()
dmeta = raw['dmeta'].copy()

tmp = (
    dmeta[['pid', 'time_label', 'displacement_status']]
    .drop_duplicates()
    .rename(columns={'pid': 'pid_label'})
    .dropna(subset=["displacement_status"], how='all')
)
dp = dp.merge(tmp, on=['pid_label', 'time_label'], how='inner', validate='many_to_one')

dp1 = dp[~dp['item_label'].str.contains('agg')].copy()
dp1['y_stan'] = dp1['y'] + 1
dp1 = dp1.merge(dit[['item_label', 'item_type']], on='item_label', how='left')

item_time_df = (
    dp1[['item_type', 'item_label', 'time']]
    .drop_duplicates()
    .sort_values(['item_type', 'time', 'item_label'])
    .reset_index(drop=True)
)
item_time_df['item_time_id'] = item_time_df.groupby('item_type').cumcount() + 1
dp1 = dp1.merge(item_time_df, on=['item_label', 'time', 'item_type'], how='left')
dp1 = dp1.merge(
    dit[['item_type', 'item_type_id']].drop_duplicates(),
    on='item_type', how='left',
)
dp1 = dp1.sort_values(['item_type_id', 'pid', 'time', 'item_label']).reset_index(drop=True)
dp1['oid'] = range(1, len(dp1) + 1)
dp1['oidt'] = dp1.groupby('item_type').cumcount() + 1

n_full = dp1['pid'].nunique()

_endline_dates = pd.to_datetime(
    dp1.loc[dp1['time_label'] == 'Endline', 'submission_date']
).dropna()
_start = _endline_dates.min().replace(day=1)
_end = (_endline_dates.max() + pd.offsets.MonthEnd(0)).normalize()
month_starts = pd.date_range(_start, _end, freq='MS')
di = pd.DataFrame({
    'interim_id': range(1, len(month_starts) + 1),
    'month_start': month_starts,
    'interim_date': month_starts + pd.offsets.MonthEnd(0),
})
interim_date = pd.to_datetime(di.loc[di['interim_id'] == INTERIM_ID, 'interim_date'].iloc[0])
print(f"  n_full={n_full} | interim {INTERIM_ID} date={interim_date.date()}")

# %%

# =============================================================================
# Fit (resume) interim-1 x-posterior + build the future block z
# =============================================================================

xi = get_interim_x(dp1, interim_date)
interim_m = n_full - xi['pid'].nunique()
print(f"  interim {INTERIM_ID}: n_x={xi['pid'].nunique()} | m (z participants added)={interim_m}")

interim_prefix = os.path.join(dir_out, f"{file_prefix}_{INTERIM_ID}")
fit = fit_partial_credit_model_ncats_pyrosvi(
    dit, xi,
    output_file_prefix=interim_prefix,
    algorithm=svi_algorithm,
    lr=0.01, num_steps=10000, output_samples=4000, seed=seed,
    x_formula=x_formula, resume=True,
    with_core_analyses=True, with_additional_analyses=False, verbose=False,
)
zi = get_interim_z_from_ypredi(
    xi, f"{interim_prefix}_draws.zarr", interim_m, pps_z_total=pps_z_total, seed=seed,
)

# Stack the x-posterior draws the loglik consumes (same scheme as fit_interim).
post = fit['draws'].posterior
theta_names = (
    'latent_factor_unit', 'latent_factor_beta',
    'skill_thresholds', 'skill_thresholds_1', 'skill_thresholds_incs',
    'loadings_questions_m1',
)
theta = {
    name: np.asarray(post[name].values).reshape(-1, *post[name].shape[2:])
    for name in theta_names if name in post
}
theta = {k: jnp.asarray(v) for k, v in theta.items()}
n_draw = int(theta['latent_factor_beta'].shape[0])
print(f"  posterior draws K={n_draw}")


def _z_stan_for_sample(s_idx: int):
    """Build the stan_data for future-data sample s (mirrors fit_interim)."""
    zcol = f'ypred_{s_idx}'
    z_dcati = zi.assign(pid=zi['src_pid'], y_stan=zi[zcol].astype(int))
    z_dcati['y'] = z_dcati['y_stan'] - 1
    z_dcati = z_dcati.sort_values(
        ['item_type_id', 'pid', 'time', 'item_label']
    ).reset_index(drop=True)
    z_dcati['oid'] = np.arange(1, len(z_dcati) + 1)
    z_dcati['oidt'] = z_dcati.groupby('item_type').cumcount() + 1
    return _fit_partial_credit_make_stan_data(
        dit=dit, dcati=z_dcati, x_formula=x_formula, verbose=False,
    )


def _ess_over_n(logw: np.ndarray) -> float:
    w = softmax(logw)
    return float(1.0 / (n_draw * np.sum(w ** 2)))


# %%

# =============================================================================
# EXP A: ESS/N vs number of assimilated z-participants (cumulative)
# EXP B: PSIS k-hat + ESS/N on the full z block
# =============================================================================

S_EXPERIMENT = min(3, pps_z_total)   # repeat over the first few future-data draws
curve_rows = []
psis_rows = []
for s_idx in range(S_EXPERIMENT):
    z_stan = _z_stan_for_sample(s_idx)
    # (K, N_obs_z) pointwise loglik; columns aligned with z_stan['unit_of_obs'].
    loglik = np.asarray(jax.vmap(lambda p: eval_loglik_partial_credit_model_ncats(z_stan, p))(theta))
    units = np.asarray(z_stan['unit_of_obs'])

    # Sum loglik per z-participant, then cumulate over participants (assimilation
    # order = ascending pid) to get logw after the first j participants.
    uniq = np.unique(units)
    ll_pid = np.stack([loglik[:, units == u].sum(axis=1) for u in uniq], axis=1)  # (K, M)
    cum_logw = np.cumsum(ll_pid, axis=1)                                          # (K, M)

    for j in range(cum_logw.shape[1]):
        curve_rows.append({
            's': s_idx + 1,
            'n_participants': j + 1,
            'ess_over_n': _ess_over_n(cum_logw[:, j]),
        })

    # EXP B: PSIS on the full-z log weights.
    logw_full = cum_logw[:, -1]
    ess_raw = _ess_over_n(logw_full)
    smoothed_lw, khat = az.psislw(logw_full.copy())
    w_ps = np.exp(np.asarray(smoothed_lw))
    w_ps = w_ps / w_ps.sum()
    ess_psis = float(1.0 / (n_draw * np.sum(w_ps ** 2)))
    khat = float(np.asarray(khat))
    psis_rows.append({
        's': s_idx + 1, 'M': int(interim_m),
        'ess_over_n_raw': ess_raw, 'ess_over_n_psis': ess_psis, 'pareto_khat': khat,
    })
    print(f"  s={s_idx+1}: ESS/N raw={ess_raw:.4f} | PSIS={ess_psis:.4f} | k-hat={khat:.2f}")

curve = pd.DataFrame(curve_rows)
psis = pd.DataFrame(psis_rows)
psis.to_csv(os.path.join(dir_out, f"{file_prefix}_v2_psis_i{INTERIM_ID}.csv"), index=False)
curve.to_csv(os.path.join(dir_out, f"{file_prefix}_v2_ess_curve_i{INTERIM_ID}.csv"), index=False)

# n_participants at which ESS/N first drops below 0.1 (per s), a rough viable horizon.
thresh = 0.1
horizon = (
    curve[curve['ess_over_n'] < thresh]
    .groupby('s')['n_participants'].min()
    .reset_index(name=f'first_n_below_{thresh}')
)
print("\n--- viable horizon (first #participants with ESS/N < 0.1) ---")
print(horizon.to_string(index=False))
print("\n--- PSIS summary (full z) ---")
print(psis.to_string(index=False))

# %%

# =============================================================================
# Plot: ESS/N degeneration curve vs #assimilated participants
# =============================================================================

curve['s'] = curve['s'].astype(str)
s_levels = sorted(curve['s'].unique())
s_colors = dict(zip(s_levels, _futurama_palette(len(s_levels))))

p = (
    ggplot(curve, aes(x='n_participants', y='ess_over_n', color='s'))
    + geom_line(size=0.8)
    + geom_hline(yintercept=thresh, colour='black', linetype='dashed', size=0.6)
    + geom_hline(yintercept=1.0 / n_draw, colour='grey', linetype='dotted', size=0.6)
    + scale_y_log10()
    + scale_color_manual(values=s_colors)
    + theme_bw()
    + theme(figure_size=(9, 5), axis_text_x=element_text())
    + labs(
        x='# z-participants assimilated (cumulative)',
        y='ESS / N  (log scale)',
        color='future-data\nsample s',
        title=f'IS weight degeneration at interim {INTERIM_ID} (m={interim_m}); '
              f'dashed=0.1, dotted=1/N',
    )
)
curve_pdf = os.path.join(dir_out, f"{file_prefix}_v2_ess_curve_i{INTERIM_ID}.pdf")
p.save(curve_pdf, verbose=False, limitsize=False)
print(f"\nSaved ESS/N degeneration curve to: {curve_pdf}")

# %%

# =============================================================================
# Shared pieces for the two mitigation experiments (C, D)
# =============================================================================
#
# Both need the full p(x | theta) factor, so build the x-cohort stan_data once,
# and a pointwise log-prior matching the model's priors (Normal/ZeroSumNormal/
# Normal/FoldedStudentT). loglik_x / loglik_z come from the exposed eval_loglik.

x_dcati = xi.copy()
if 'y_stan' not in x_dcati:
    x_dcati['y_stan'] = x_dcati['y'] + 1
x_dcati = x_dcati.sort_values(
    ['item_type_id', 'pid', 'time', 'item_label']
).reset_index(drop=True)
x_dcati['oid'] = np.arange(1, len(x_dcati) + 1)
x_dcati['oidt'] = x_dcati.groupby('item_type').cumcount() + 1
x_stan = _fit_partial_credit_make_stan_data(
    dit=dit, dcati=x_dcati, x_formula=x_formula, verbose=False,
)

U = int(theta['latent_factor_unit'].shape[1])
P = int(theta['latent_factor_beta'].shape[1])
L = int(theta['skill_thresholds'].shape[1])
Ld = int(theta['loadings_questions_m1'].shape[1])
s2z_sd_unit = 1.0 / np.sqrt(1.0 - 1.0 / U)
_eval = eval_loglik_partial_credit_model_ncats


def _logprior(pr: dict):
    lp = jstats.norm.logpdf(pr['latent_factor_unit'], 0.0, s2z_sd_unit).sum()
    lp += jstats.norm.logpdf(pr['latent_factor_beta'], 0.0, 1.0).sum()
    lp += jstats.norm.logpdf(pr['skill_thresholds'], 0.0, 3.5).sum()
    lp += (jnp.log(2.0) + jstats.t.logpdf(pr['loadings_questions_m1'], 3.0, 0.0, 1.0)).sum()
    return lp


# %%

# =============================================================================
# EXP C: moment-matching importance sampling (Paananen et al. 2021), mean-match
# =============================================================================
#
# Work in a "match space" u where the positive loadings are log-transformed and
# all other params are identity. Approximate the proposal p(theta|x) by a
# diagonal Gaussian fit to the x-posterior draws, q(u) = N(u; mu_p, diag sd_p^2)
# (this Gaussian approximation is intrinsic to applying MMIS without the exact
# variational density). Importance weights in match space:
#     log w(u) = logprior(theta(u)) + loglik_x(theta(u)) + loglik_z(theta(u))
#                + sum(u_loadings)            # Jacobian of the log-transform
#                - logq(u)
# Mean-match step: shift every draw by (mu_weighted - mu_proposal), then
# recompute weights. If the base weights already collapse, the weighted mean is
# essentially the single dominating draw, so the shift cannot rescue ESS.

# Stack draws into the match-space matrix Theta_u (K, D).
Theta_u = np.concatenate([
    np.asarray(theta['latent_factor_unit']),
    np.asarray(theta['latent_factor_beta']),
    np.asarray(theta['skill_thresholds']),
    np.log(np.asarray(theta['loadings_questions_m1'])),
], axis=1)
mu_p = Theta_u.mean(axis=0)
sd_p = Theta_u.std(axis=0) + 1e-8
sl_u = slice(0, U)
sl_b = slice(U, U + P)
sl_s = slice(U + P, U + P + L)
sl_l = slice(U + P + L, U + P + L + Ld)


def _u_to_params(u):
    return {
        'latent_factor_unit': u[sl_u],
        'latent_factor_beta': u[sl_b],
        'skill_thresholds': u[sl_s],
        'loadings_questions_m1': jnp.exp(u[sl_l]),
    }


def _logq(u):
    return jstats.norm.logpdf(u, jnp.asarray(mu_p), jnp.asarray(sd_p)).sum()


def _make_logw(z_stan):
    def logw(u):
        pr = _u_to_params(u)
        lt = (_logprior(pr) + _eval(x_stan, pr).sum() + _eval(z_stan, pr).sum()
              + u[sl_l].sum())          # log|d loading / d u|
        return lt - _logq(u)
    return jax.jit(jax.vmap(logw))


mm_rows = []
Theta_u_j = jnp.asarray(Theta_u)
for s_idx in range(S_EXPERIMENT):
    z_stan = _z_stan_for_sample(s_idx)
    logw_fn = _make_logw(z_stan)

    logw0 = np.asarray(logw_fn(Theta_u_j))
    w0 = softmax(logw0)
    ess0 = float(1.0 / (n_draw * np.sum(w0 ** 2)))
    # Sanity: how well the diag-Gaussian base weights track the exact IS weights
    # (exact per-draw z log-likelihood).
    llz_exact = np.asarray(jax.vmap(lambda p: _eval(z_stan, p).sum())(theta))
    corr = float(np.corrcoef(logw0, llz_exact)[0, 1])

    mu_w = (w0[:, None] * Theta_u).sum(axis=0)          # weighted mean
    shift = mu_w - mu_p
    Theta_u_star = jnp.asarray(Theta_u + shift[None, :])
    logw1 = np.asarray(logw_fn(Theta_u_star))
    w1 = softmax(logw1)
    ess1 = float(1.0 / (n_draw * np.sum(w1 ** 2)))

    mm_rows.append({
        's': s_idx + 1,
        'ess_over_n_base': ess0,
        'ess_over_n_meanmatch': ess1,
        'basew_vs_exact_corr': corr,
    })
    print(f"  [MM] s={s_idx+1}: ESS/N base={ess0:.4f} -> mean-match={ess1:.4f} "
          f"(base-weight vs exact corr={corr:.2f})")

mm = pd.DataFrame(mm_rows)
mm.to_csv(os.path.join(dir_out, f"{file_prefix}_v2_momentmatch_i{INTERIM_ID}.csv"), index=False)
print("\n--- moment-matching summary ---")
print(mm.to_string(index=False))

# %%

# =============================================================================
# EXP D: SMC sampler with resample-move (adaptive data tempering + NUTS moves)
# =============================================================================
#
# Bridge p(theta|x) -> p(theta|x, z_s) through pi_beta ∝ p(theta|x) p(z_s|theta)^beta,
# beta: 0 -> 1. Particles start at the x-posterior draws (subsampled). Each step:
#   1. adaptively pick the next beta so the tempering ESS ≈ 0.5 * K_smc,
#   2. resample particles by the incremental weights,
#   3. MOVE them with a few NUTS steps invariant to pi_beta (this is what plain
#      reweighting cannot do — it relocates particles into the target's support).
# Demonstrated for the first future-data sample only.

K_SMC = 128
ESS_FRAC_TARGET = 0.5
NUTS_WARMUP = 60
MAX_TEMPS = 15
S_SMC = 1

z_stan_smc = _z_stan_for_sample(S_SMC - 1)


def _tempered_model(beta):
    def model():
        lfu = numpyro.sample('latent_factor_unit',
                             dist.ZeroSumNormal(scale=s2z_sd_unit, event_shape=(U,)))
        beta_p = numpyro.sample('latent_factor_beta', dist.Normal(0, 1).expand([P]))
        skill = numpyro.sample('skill_thresholds', dist.Normal(0, 3.5).expand([L]))
        load = numpyro.sample(
            'loadings_questions_m1',
            dist.FoldedDistribution(dist.StudentT(3.0, 0.0, 1.0)).expand([Ld]),
        )
        pr = {'latent_factor_unit': lfu, 'latent_factor_beta': beta_p,
              'skill_thresholds': skill, 'loadings_questions_m1': load}
        numpyro.factor('llx', _eval(x_stan, pr).sum())
        numpyro.factor('llz', beta * _eval(z_stan_smc, pr).sum())
    return model


def _systematic_resample(w, key):
    n = len(w)
    positions = (random.uniform(key) + np.arange(n)) / n
    return np.searchsorted(np.cumsum(w), positions)


def _llz_at(particles):
    return np.asarray(jax.vmap(lambda p: _eval(z_stan_smc, p).sum())(particles))


# Init particles: subsample K_SMC of the x-posterior draws.
rng = np.random.default_rng(seed)
idx0 = rng.choice(n_draw, size=K_SMC, replace=False)
particles = {k: jnp.asarray(np.asarray(v)[idx0]) for k, v in theta.items()}
llz = _llz_at(particles)

key = random.PRNGKey(seed)
beta = 0.0
smc_rows = []
print(f"\n[SMC] resample-move on s={S_SMC}, K={K_SMC} particles")
for t in range(MAX_TEMPS):
    # Adaptive next beta: bisection so tempering ESS-fraction ≈ ESS_FRAC_TARGET.
    def ess_frac(db):
        wn = softmax(db * llz)
        return 1.0 / (K_SMC * np.sum(wn ** 2))
    lo, hi = 0.0, 1.0 - beta
    if ess_frac(hi) >= ESS_FRAC_TARGET:
        db = hi
    else:
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if ess_frac(mid) >= ESS_FRAC_TARGET:
                lo = mid
            else:
                hi = mid
        db = lo
    beta_new = beta + db

    wn = softmax(db * llz)
    ess_temper = float(1.0 / (K_SMC * np.sum(wn ** 2)))

    # Resample by incremental weights.
    key, ksub = random.split(key)
    anc = _systematic_resample(wn, ksub)
    particles = {k: v[anc] for k, v in particles.items()}

    # Move: NUTS invariant to pi_{beta_new}, warm-started at the particles.
    key, kmcmc = random.split(key)
    mcmc = MCMC(NUTS(_tempered_model(beta_new), max_tree_depth=6),
                num_warmup=NUTS_WARMUP, num_samples=1,
                num_chains=K_SMC, chain_method='vectorized', progress_bar=False)
    mcmc.run(kmcmc, init_params={k: np.asarray(v) for k, v in particles.items()})
    last = mcmc.get_samples(group_by_chain=True)
    particles = {k: last[k][:, -1] for k in particles}  # (K_SMC, *event)
    llz = _llz_at(particles)

    smc_rows.append({
        'temp': t + 1, 'beta': round(beta_new, 4),
        'd_beta': round(db, 4), 'ess_frac_temper': round(ess_temper, 3),
    })
    print(f"  temp {t+1}: beta {beta:.3f}->{beta_new:.3f} (dbeta={db:.3f}) "
          f"| tempering ESS/K={ess_temper:.2f}")
    beta = beta_new
    if beta >= 1.0 - 1e-9:
        break

smc = pd.DataFrame(smc_rows)
smc.to_csv(os.path.join(dir_out, f"{file_prefix}_v2_smc_i{INTERIM_ID}.csv"), index=False)
print("\n--- SMC schedule ---")
print(smc.to_string(index=False))
print(f"  reached beta={beta:.3f} in {len(smc)} tempering steps; "
      f"post-move particles are uniformly weighted (effective size = K = {K_SMC}).")

# %%

print("\n✓ IS v2 degeneracy experiment complete")
