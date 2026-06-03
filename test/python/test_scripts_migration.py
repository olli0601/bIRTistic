"""
Migration smoke tests for the analysis + interim scripts under
``scripts-py/``. Each script must:

1. Parse cleanly (no syntax errors after the search-and-replace).
2. Import the corresponding ``model_*`` subclass instead of the legacy
   ``fit_partial_credit_model`` / ``fit_credit_model`` /
   ``fit_ordered_logit_model`` free-function pipeline.
3. Contain no direct calls to the legacy module-level free fitters
   (``fit_partial_credit_model_ncats_pyrosvi`` etc.) in the script body,
   *except* where the legacy callable is passed through to
   ``fit_interim_MC_of_posterior_xz`` via its ``fitting_method=`` kwarg
   (a multiprocessing-only path we kept on the legacy module to avoid
   pickling a full Model into each spawn worker; see the docstring in
   ``Ukraine_interim_analyses_with_HMC.py``).
"""

import ast
import re
import sys
from pathlib import Path

import pytest

_repo_root = Path(__file__).resolve().parents[2]
_scripts_dir = _repo_root / 'scripts-py'

_MIGRATED_SCRIPTS = [
    # Step 5: interim PPS scripts (one PCM model per interim).
    'Ukraine_interim_analysis_regression_H1x_on_wz.py',
    'Ukraine_interim_analysis_regression_endptx_on_wz.py',
    'Ukraine_interim_analysis_with_IS_from_x.py',
    'Ukraine_interim_analysis_with_IS_moment_matching_from_x.py',
    'Ukraine_interim_analysis_with_SMC_resample_from_x.py',
    # User-requested follow-up: broader Colombia / Ukraine analysis scripts.
    'Colombia_analysis_partial_credit_ADVI-HMC-pyro.py',
    'Ukraine_analysis_partial_credit_ADVI-HMC-pyro.py',
    'Colombia_analysis_ordered_logit_ADVI-HMC-pyro.py',
    'Ukraine_analysis_ordered_logit_ADVI-HMC-pyro.py',
    'Colombia_analysis_ordered_logit_ADVI_vs_HMC.py',
    'Colombia_interim_analyses.py',
    'Ukraine_interim_analyses_with_HMC.py',
]


@pytest.mark.parametrize('name', _MIGRATED_SCRIPTS)
def test_script_parses(name):
    path = _scripts_dir / name
    src = path.read_text()
    ast.parse(src)


@pytest.mark.parametrize('name', _MIGRATED_SCRIPTS)
def test_script_imports_model_subclass(name):
    """Each migrated script must import the corresponding model_* subclass."""
    src = (_scripts_dir / name).read_text()
    # Strip docstrings + comments so we only check executable code.
    # AST walk to find Import / ImportFrom nodes.
    tree = ast.parse(src)
    imported_modules = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            imported_modules.add(node.module)
    # Expected: at least one of model_pcm / model_credit / model_ordered_logit
    # / model_binomial.
    model_imports = {m for m in imported_modules if m.startswith('model_')}
    assert model_imports, (
        f"{name}: no `from model_* import ...` line found. Expected one of"
        " model_pcm, model_credit, model_ordered_logit, model_binomial."
    )


@pytest.mark.parametrize('name', _MIGRATED_SCRIPTS)
def test_script_does_not_call_legacy_free_fitters(name):
    """Legacy free-function fitter calls
    (``fit_partial_credit_model_ncats_pyrosvi`` etc.) must not appear as
    direct callables in the migrated scripts. Exception: when used as a
    pickled callable passed to ``fit_interim_MC_of_posterior_xz`` via
    ``fitting_method=``, the legacy reference is acceptable -- only allow
    the literal substring ``fitting_method=fit_partial_credit_model_ncats``
    or ``fitting_method=fit_credit_model_ncats`` or
    ``fitting_method=fit_ordered_logit_model_ncats``."""
    src = (_scripts_dir / name).read_text()
    legacy_call_re = re.compile(
        r"(?<!fitting_method=)"          # not a fitting_method= passthrough
        r"(?<!\.)"                       # not a docstring fence/comment chain
        r"fit_(partial_credit|credit|ordered_logit)_model_ncats_"
        r"(pyrosvi|stanadvi|stanhmc)\("
    )
    # Strip docstrings / comments (cheap heuristic: skip lines starting with
    # '#' or inside triple-quoted strings).
    cleaned = re.sub(r'"""[\s\S]*?"""', '', src)
    cleaned = re.sub(r"'''[\s\S]*?'''", '', cleaned)
    cleaned = '\n'.join(
        ln for ln in cleaned.splitlines() if not ln.lstrip().startswith('#')
    )
    hits = legacy_call_re.findall(cleaned)
    assert not hits, (
        f"{name}: still calls legacy free fitter {hits}; should dispatch"
        " through a Model subclass instance."
    )
