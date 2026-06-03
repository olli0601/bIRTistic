# Plan: `Model` blueprint for model-agnostic interim algorithms

Status: design only — no code yet. Goal: split PCM-specific math out of
[python/fit_interim.py](../python/fit_interim.py) and out of
[python/fit_partial_credit_model.py](../python/fit_partial_credit_model.py)
so the same IS / MM / SMC / regression algorithms run on PCM, credit,
ordered-logit, and Binomial models with one subclass each.

## 1. Why

Two files mix generic plumbing with PCM-specific math:

- `fit_interim.py`: IS / MM / SMC / regression algorithms hard-code PCM
  loglik / prior / ordered-prob and a PCM-only param-name list.
- `fit_partial_credit_model.py`: holds the three fit drivers
  (`_pyrosvi`, `_stanadvi`, `_stanhmc`) plus loglik / prior / stan-data
  helpers. Same pattern duplicated in `fit_credit_model.py` and
  `fit_ordered_logit_model.py`.

The plan: one `Model` ABC + four subclasses (PCM, credit, ordered-logit,
Binomial). All current free-function call sites in `fit_interim.py` and
the scripts route through the model instance.

## 2. The `Model` blueprint

`python/model.py` (NEW). The base class holds **data**, exposes the three
**fitters**, and declares the **model-math callables** every algorithm
needs.

```python
class Model(ABC):
    """Generative model + its data, the unit every interim algorithm
    consumes. Holds dit, dcati, stan_data on the instance."""

    # ---- 2.1 Subclass-declared metadata ----
    param_names: tuple = ()
    """Posterior-draw parameter names this model expects. Used by
    stack_posterior_theta to flatten arviz draws into a dict of jnp arrays."""

    positive_params: tuple = ()
    """Subset of param_names that live on the positive reals -- log-
    transformed into the MM 'match space' (Paananen 2021)."""

    # ---- 2.2 Construction ----
    def __init__(self, dit: pd.DataFrame, dcati: pd.DataFrame,
                 x_formula: str = "~ time - 1", *, seed: int = 123):
        self.dit = dit
        self.dcati = dcati
        self.x_formula = x_formula
        self.seed = seed
        self.stan_data = self.make_stan_data(dcati, x_formula)

    @abstractmethod
    def make_stan_data(self, dcati: pd.DataFrame, x_formula: str) -> dict:
        """Build the data dict that eval_loglik / eval_outcome consume.
        Reused inside interim algorithms with a z-cohort dcati."""

    # ---- 2.3 Fitters (generic across models) ----
    @abstractmethod
    def fit_pyro_svi(self, output_file_prefix: str, *,
                     save_to_file: bool = True, resume: bool = False,
                     verbose: bool = True, **method_kwargs) -> dict: ...
    @abstractmethod
    def fit_stan_svi(self, output_file_prefix: str, *,
                     save_to_file: bool = True, resume: bool = False,
                     verbose: bool = True, **method_kwargs) -> dict: ...
    @abstractmethod
    def fit_stan_hmc(self, output_file_prefix: str, *,
                     save_to_file: bool = True, resume: bool = False,
                     verbose: bool = True, **method_kwargs) -> dict: ...
    def fit_closed_form(self, output_file_prefix: str, *,
                        save_to_file: bool = True, verbose: bool = True,
                        **method_kwargs) -> dict:
        raise NotImplementedError(
            f"{type(self).__name__} has no closed-form posterior."
        )

    Each fitter returns the same dict shape as today
    ({'draws': arviz_idata, 'timing': ..., 'posterior_samples': ...}).
    Method-specific tuning (lr, num_steps; iter, grad_samples; chains,
    iter_warmup) flows through ``**method_kwargs`` -- no model-specific
    leakage into the API.

    # ---- 2.4 Posterior <-> params ----
    def stack_posterior_theta(self, draws) -> dict:
        """Default implementation: flatten draws.posterior[name] for every
        ``name in self.param_names`` to (K, ...) jnp arrays. Subclass only
        needs to set param_names."""
        ...

    # ---- 2.5 Likelihood + prior ----
    @abstractmethod
    def eval_loglik(self, data: dict, params: dict) -> jnp.ndarray: ...
    @abstractmethod
    def eval_loglik_annealed(self, data: dict, params: dict) -> jnp.ndarray:
        """params['temperature'] in [0, 1] scales the loglik for SMC."""
    @abstractmethod
    def logprior(self, params: dict) -> jnp.ndarray: ...

    # ---- 2.6 Endpoint quantity per draw ----
    @abstractmethod
    def eval_outcome_for_endpoint(self, data: dict, params: dict) -> jnp.ndarray:
        """Whatever quantity the per-draw endpoint code consumes. PCM:
        ordered_prob_by_cat_qu_fit. Binomial: p(theta)."""
    @abstractmethod
    def endpoints_per_draw(self, dcati: pd.DataFrame, draws,
                           categorical_threshold: int,
                           endpoint_type: str = 'items') -> pd.DataFrame:
        """One row per (item, draw): item_label, item_type, item_high_label,
        draw, ratio."""
```

That's the whole base class. Algorithms do NOT live on it.

## 3. Subclasses

### 3.1 IRT family (PCM, credit, ordered-logit)

All three follow the same wrapper pattern. Factor the shared boilerplate
into a private `_IRTModel(Model)` mixin so each concrete IRT subclass is
~20 lines.

```python
class _IRTModel(Model):
    """Boilerplate shared by PCM, credit, ordered-logit. Subclass declares
    _module (the numpyro module functions live there) and _make_stan_data
    (the IRT-specific stan-data builder)."""

    @property
    @abstractmethod
    def _module(self): ...    # the .pyro / .py namespace

    def make_stan_data(self, dcati, x_formula):
        return self._make_stan_data(dit=self.dit, dcati=dcati,
                                    x_formula=x_formula, verbose=False)

    def fit_pyro_svi(self, output_file_prefix, **kw):
        return self._module.fit_pyrosvi(self.dit, self.dcati,
                                        output_file_prefix, **kw)
    def fit_stan_svi(self, output_file_prefix, **kw):
        return self._module.fit_stanadvi(self.dit, self.dcati,
                                         output_file_prefix, **kw)
    def fit_stan_hmc(self, output_file_prefix, **kw):
        return self._module.fit_stanhmc(self.dit, self.dcati,
                                        output_file_prefix, **kw)

    def eval_loglik(self, data, params):
        return self._module.eval_loglik(data, params)
    def eval_loglik_annealed(self, data, params):
        return self._module.eval_loglik_with_annealing(data, params)
    def logprior(self, params):
        return self._module.get_prior(params)
    def eval_outcome_for_endpoint(self, data, params):
        return self._module.get_ordered_prob(data, params)
    def endpoints_per_draw(self, dcati, draws, categorical_threshold,
                           endpoint_type='items'):
        return get_endpoints_per_draw(
            dcati=dcati, dit=self.dit, draws=draws,
            categorical_threshold=categorical_threshold,
            endpoint_type=endpoint_type,
            param_name='ordered_prob_by_cat_qu_fit', verbose=False,
        )
```

Each concrete IRT subclass: just `_module`, `_make_stan_data`,
`param_names`, `positive_params`.

```python
class PartialCreditModelNCats(_IRTModel):
    param_names = ('latent_factor_unit', 'latent_factor_beta',
                   'skill_thresholds', 'loadings_questions_m1')
    positive_params = ('loadings_questions_m1',)
    _module = pcm_module   # currently python/fit_partial_credit_model.py
    _make_stan_data = staticmethod(_fit_partial_credit_make_stan_data)

class CreditModelNCats(_IRTModel):
    param_names = (...)        # different param list
    positive_params = (...)
    _module = credit_module
    _make_stan_data = staticmethod(_fit_credit_make_stan_data)

class OrderedLogitNCats(_IRTModel):
    param_names = (...)
    positive_params = (...)
    _module = ordered_logit_module
    _make_stan_data = staticmethod(_fit_ordered_logit_make_stan_data)
```

### 3.2 Binomial

Different data shape, but still implements all three fitters. SVI / HMC
target the same Beta-Binomial generative model (numpyro + cmdstanpy), so
the user can test them against the closed-form Beta posterior.

```python
class BinomialModel(Model):
    param_names = ('p',)
    positive_params = ()                # logit-transform if MM needs it

    def make_stan_data(self, dcati, x_formula):
        return {
            'y': dcati['y'].astype(int).to_numpy(),
            'N': len(dcati),
            'item_label': dcati['item_label'].to_numpy(),
        }

    # fit_pyro_svi / fit_stan_svi / fit_stan_hmc: implement against the
    # numpyro Beta-Binomial program and the Stan equivalent.
    def fit_pyro_svi(self, output_file_prefix, **kw): ...
    def fit_stan_svi(self, output_file_prefix, **kw): ...
    def fit_stan_hmc(self, output_file_prefix, **kw): ...

    def fit_closed_form(self, output_file_prefix, **kw):
        """Beta(a + k, b + n - k) posterior, sampled directly. Returns
        the same dict shape as the SVI/HMC fitters."""
        ...

    def eval_loglik(self, data, params):       # log Bernoulli(y | p)
        ...
    def eval_loglik_annealed(self, data, params): ...
    def logprior(self, params):                # log Beta(p; a, b)
        ...
    def eval_outcome_for_endpoint(self, data, params):
        return params['p']                     # p IS the endpoint
    def endpoints_per_draw(self, dcati, draws, categorical_threshold,
                           endpoint_type='items'):
        # No items / no categorical threshold here -- single-quantity case.
        ...
```

The IRT subclasses raise `NotImplementedError` from
`fit_closed_form` (inherited default). Binomial overrides it.

## 4. Algorithms (free functions, take a `Model`)

`python/fit_interim.py` keeps the algorithms as free functions but now
each takes a `Model` instance. Functions read `model.dit`,
`model.dcati`, `model.stan_data`, and call `model.eval_loglik(...)`
etc. instead of importing the PCM module-level callables.

```python
def fit_interim_IS_reweight(model, draws, zi, *, pps_z_total, ...): ...
def fit_interim_IS_moment_matching(model, draws, zi, *, pps_z_total, ...): ...
def fit_interim_SMC_resample(model, draws, zi, *, pps_z_total, ...): ...
def fit_interim_SMC_PPS(model, draws, zi, *, pps_z_total, ...): ...
def get_interim_endpt_and_w_from_poi(model, draws, draws_file, interim_m,
                                     pps_z_total, ...): ...
# Pure, no model:
def fit_interim_regress_H1x_on_wz(wa): ...
def fit_interim_regress_endptx_on_wz(wa, pps_H1_def=0.5): ...
```

Z-cohort handling: inside the inner s_idx loop, an algorithm needs a
stan_data dict for `zi[s]`. It calls **`model.make_stan_data(dcati=zi_for_s_idx)`**
to produce one — no fresh Model instance per s_idx (would re-validate /
re-build dit-derived structures for nothing). The current PCM-specific
helpers `_interim_make_x_stan` and `_interim_make_z_stan` collapse into
that one method.

## 5. Open questions — resolved

a. **z-cohort handling**: `model.make_stan_data(dcati=zi_for_s_idx)`.
   No fresh Model per s_idx. All current `_interim_make_x_stan` /
   `_interim_make_z_stan` references go away.

b. **`stack_posterior_theta`**: each subclass declares `param_names`;
   the base class implements `stack_posterior_theta` once using it.

c. **`endpoints_per_draw`**: required on every subclass (Binomial gets
   tested with HMC/SVI against its closed-form, so it must expose
   per-draw endpoints too — single-quantity case, no items / threshold).

d. **`fit_closed_form`**: lives on the base, defaults to
   `NotImplementedError`. Binomial overrides; PCM / credit /
   ordered-logit inherit the default.

## 6. Module layout after refactor

```
python/
  model.py                       # ABC + _IRTModel mixin
  model_pcm.py                   # PartialCreditModelNCats
  model_credit.py                # CreditModelNCats
  model_ordered_logit.py         # OrderedLogitNCats
  model_binomial.py              # BinomialModel
  fit_interim.py                 # free-function algorithms taking a Model
  interim_helpers.py             # _load_ypred, get_interim_x,
                                 # get_interim_z_from_ypred{i,f}
  # GONE: fit_partial_credit_model.py, fit_credit_model.py,
  # fit_ordered_logit_model.py.
  # Their three fit drivers move onto the subclass; loglik / prior /
  # stan-data helpers stay in the underlying numpyro modules (.pyro
  # files) and are called via the subclass.
```

Scripts change from:

```python
from fit_interim import (
    get_interim_x, get_interim_endpt_and_w_from_poi,
    fit_interim_regress_H1x_on_wz,
)
from fit_partial_credit_model import fit_partial_credit_model_ncats_pyrosvi
...
fit = fit_partial_credit_model_ncats_pyrosvi(dit, xi, ...)
wa = get_interim_endpt_and_w_from_poi(xi=xi, dit=dit, draws=fit['draws'], ...)
p, perf = fit_interim_regress_H1x_on_wz(wa)
```

to:

```python
from interim_helpers import get_interim_x
from model_pcm import PartialCreditModelNCats
from fit_interim import (
    get_interim_endpt_and_w_from_poi, fit_interim_regress_H1x_on_wz,
)
...
model = PartialCreditModelNCats(dit, xi, x_formula="~ time - 1")
fit = model.fit_pyro_svi(interim_prefix, num_steps=10000, lr=0.01, ...)
wa = get_interim_endpt_and_w_from_poi(model, draws=fit['draws'], ...)
p, perf = fit_interim_regress_H1x_on_wz(wa)
```

## 7. Migration steps

Each step is one file / one concern. Two baseline tests run after every
step (must stay green):

- `test/python/test_interim_alignment.py` — 5 alignment invariants
- `test/python/test_interim_determinism.py` — 8 determinism invariants

Each step also **adds** new tests that pin down behaviour preserved by
that step's refactor. New tests stay forever; they form the regression
shield for downstream steps.

### Step 1 — move pure-pandas helpers to `interim_helpers.py`

Action: relocate `_load_ypred`, `get_interim_x`,
`get_interim_z_from_ypred{i,f}`. No logic change.

Tests added:
- `test/python/test_interim_helpers.py::test_imports_from_new_module` —
  asserts the four functions import from `interim_helpers` and resolve
  to the same callables as the old `fit_interim` re-exports.
- `test/python/test_interim_helpers.py::test_get_interim_z_from_ypredi_golden` —
  hash/snapshot the output frame for the fixture cohort (interim_m=12,
  pps_z_total=8, seed=123). Locks the exact column / row layout so any
  later refactor preserves it.

Gate: alignment + determinism + the two new tests.

### Step 2 — define `Model` ABC + `_IRTModel` mixin

Action: write `python/model.py` with the abstract base + the IRT mixin.
No subclass yet.

Tests added:
- `test/python/test_model_abc.py::test_cannot_instantiate_base` —
  `Model(dit, dcati)` raises `TypeError` for abstract methods.
- `test/python/test_model_abc.py::test_irt_mixin_abstract` — same for
  `_IRTModel` (must subclass with `_module` + `_make_stan_data`).

Gate: imports clean + the two new tests.

### Step 3 — implement `PartialCreditModelNCats` (delegating)

Action: write `python/model_pcm.py` that delegates every blueprint method
to the existing module-level functions in `fit_partial_credit_model.py`
+ the `.pyro` namespace. Zero algorithmic change.

Tests added: `test/python/test_model_pcm.py`
- `test_make_stan_data_matches_free_function` — `model.make_stan_data(dcati)`
  equals `_fit_partial_credit_make_stan_data(dit, dcati, ...)` key-by-key,
  including dtypes.
- `test_eval_loglik_matches_free_function` — at a fixed param dict
  (loaded from the fixture draws.zarr posterior mean), the two return
  bit-identical arrays.
- `test_logprior_matches_free_function`, `test_eval_outcome_matches_free_function`,
  `test_endpoints_per_draw_matches_free_function` — same pattern.
- `test_stack_posterior_theta_param_names` — flattened dict has exactly
  the keys `PartialCreditModelNCats.param_names` and array shapes
  `(K, ...)`.
- `test_fit_pyro_svi_matches_free_function` — small `num_steps=200` fit
  via `model.fit_pyro_svi(...)` returns the same `draws` array values
  as the direct free-function call with identical seed.

Gate: alignment + determinism + all `test_model_pcm.py` tests.

### Step 4 — refactor `fit_interim.py` algorithms onto a `Model`

Action: convert each algorithm function to take `model` as its first
arg, internally calling `model.eval_loglik(...)` etc. Keep the old
free-function signature as a thin shim that instantiates a default
`PartialCreditModelNCats` and forwards — one cycle of backward
compatibility, removed in Step 8.

Order: IS reweight → IS moment-match → SMC resample → SMC PPS →
`get_interim_endpt_and_w_from_poi`. One commit per algorithm.

Tests added per algorithm:
- `test/python/test_interim_algorithms.py::test_{name}_model_matches_shim` —
  on the fixture, call the new `model`-based form and the shim form;
  assert outputs bit-identical (`pd.testing.assert_frame_equal` with
  `check_exact=True`).
- `test_{name}_runs_on_alternate_model` — once Credit and OrderedLogit
  exist (after Step 6), parameterize this test over all IRT subclasses
  to prove the algorithm really is model-agnostic.

Production-script gate (per algorithm): run the corresponding script
(IS-from-x, MM-from-x, SMC-resample-from-x, regression) and diff the
output PDFs / CSVs / pkls against the cached versions in
`/Users/or105/sandbox/bIRTistic/py-ukraine-interim-with-*-260601/`. Any
non-floating-point difference fails the gate; floats must agree at
`atol=1e-10`. Wire this into `scripts-py/test_interim_outputs_diff.py`
(new, run manually before merging the algorithm step).

Gate: alignment + determinism + `test_interim_algorithms.py` for that
algorithm + script-diff gate.

### Step 5 — migrate scripts

Action: update each script to instantiate `PartialCreditModelNCats` and
call the model-aware algorithms (the new free-function form, not the
shim). One script per commit.

Tests added:
- `test/python/test_scripts_smoke.py::test_{script}_smoke` — runs each
  updated script under `pytest -m slow` with `pps_z_total=4` against
  the fixture, asserting the script exits 0 and produces the expected
  output files. Catches import / signature breakage; not a numerical
  check.

Gate: alignment + determinism + smoke test for the migrated script + the
script-diff gate from Step 4 (re-run against cached production
outputs).

### Step 6 — add `CreditModelNCats` and `OrderedLogitNCats`

Action: write `python/model_credit.py` and `python/model_ordered_logit.py`,
analogous to the PCM subclass. Existing fixtures in
`test/test_data/` cover PCM; reuse / add small fixtures for credit +
ordered-logit.

Tests added: two files mirroring `test_model_pcm.py` —
- `test/python/test_model_credit.py` — same six tests as PCM, all
  comparing the new subclass methods against the existing free-function
  call sites in `fit_credit_model.py`.
- `test/python/test_model_ordered_logit.py` — same pattern.
- Parametrize `test_interim_algorithms.py::test_{name}_runs_on_alternate_model`
  (added in Step 4) to include Credit + OrderedLogit.

Gate: all above + alignment + determinism (which now also run with the
two new subclasses as fixtures).

### Step 7 — add `BinomialModel`

Action: write `python/model_binomial.py`. New numpyro program for
SVI / HMC + Stan code if we want fit_stan_*. `fit_closed_form` samples
directly from the Beta posterior.

Tests added: `test/python/test_model_binomial.py`
- `test_closed_form_matches_analytic_pps` — synthetic data with known
  Beta(a, b) prior + binomial likelihood; assert
  `|PPS_closed_form - PPS_analytic| < 1e-8`. Analytic PPS from §3.1 of
  `dev/amortised_decision_making.md` (closed-form Beta-Binomial tail
  sum).
- `test_pyro_svi_vs_closed_form` — SVI fit at `num_steps=20_000` agrees
  with closed-form Beta within Monte-Carlo tol on posterior mean and
  PPS (`atol=0.01`).
- `test_stan_hmc_vs_closed_form` — same, HMC fit at
  `iter_sampling=2000`.
- `test_binomial_with_interim_algorithms` — run IS reweight + regression
  label estimators on synthetic Binomial data; PPS within Monte-Carlo
  tol of the closed-form answer. Proves the IRT-trained algorithms work
  unchanged on a non-IRT model.

Gate: all above. This is the first end-to-end proof that the blueprint
generalises.

### Step 8 — delete legacy modules + shims

Action: remove `fit_partial_credit_model.py`, `fit_credit_model.py`,
`fit_ordered_logit_model.py` and the free-function shims left in
`fit_interim.py` from Step 4.

Tests added:
- `test/python/test_no_legacy_imports.py::test_grep_legacy_modules` —
  greps the repo (excluding `dev/`, `htmlcov/`, `.pixi/`,
  `test/test_data/`) for any import of the three removed modules. Must
  match zero lines.
- `test/python/test_no_legacy_imports.py::test_no_shim_signatures` —
  greps `fit_interim.py` for the Step-4 shim signatures (e.g.
  `def fit_interim_IS_reweight(xi, zi, dit, draws,`). Must match zero.

Gate: full test suite green + grep tests + every production script runs
on a small smoke fixture (Step 5 smoke tests, expanded).

### Coverage summary

| Step | Files of new tests | Existing tests still gating |
|------|--------------------|------------------------------|
| 1    | `test_interim_helpers.py` | alignment, determinism |
| 2    | `test_model_abc.py` | alignment, determinism |
| 3    | `test_model_pcm.py` | alignment, determinism |
| 4    | `test_interim_algorithms.py` (grows per algorithm) | alignment, determinism, `test_model_pcm.py` |
| 5    | `test_scripts_smoke.py` | everything above |
| 6    | `test_model_credit.py`, `test_model_ordered_logit.py` | everything above |
| 7    | `test_model_binomial.py` | everything above |
| 8    | `test_no_legacy_imports.py` | everything above |
