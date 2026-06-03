"""
Step 8 of the OO-port refactor: forbid stale references to the deleted
back-compat shims and to the legacy ``fit_*_model`` free-function
fitters in user-facing code (scripts + algorithm library).

The legacy module files themselves (python/fit_partial_credit_model.py,
fit_credit_model.py, fit_ordered_logit_model.py) are still on disk:
the ``Model`` subclasses in python/model_pcm.py / model_credit.py /
model_ordered_logit.py delegate to them via the ``_module`` dict. Only
those subclass files (and python/fit_interim.py for the
``fit_interim_MC_of_posterior_xz`` multiprocessing default) are allowed
to import from the legacy modules.
"""

import re
import sys
from pathlib import Path

import pytest

_repo_root = Path(__file__).resolve().parents[2]
_python_dir = _repo_root / 'python'
_scripts_dir = _repo_root / 'scripts-py'

# Files allowed to keep an internal reference to the legacy free-function
# pipeline -- the Model subclasses (which delegate to it) and the one
# fit_interim multiprocessing default that predates the refactor.
_ALLOWED_LEGACY_IMPORTERS = {
    'model_pcm.py',
    'model_credit.py',
    'model_ordered_logit.py',
    'fit_interim.py',
}

_LEGACY_MODULES_RE = re.compile(
    r"\bfrom\s+(fit_partial_credit_model|fit_credit_model|fit_ordered_logit_model)\s+import"
)

_DELETED_SHIM_RE = re.compile(
    r"\bfit_interim_(importance_sampling_of_posterior_xz_from_x|"
    r"IS_moment_matching_of_posterior_xz_from_x|"
    r"SMC_resample_of_posterior_xz_from_x)\b"
)


def _iter_py_files(root: Path):
    for p in root.rglob('*.py'):
        if '__pycache__' in p.parts:
            continue
        yield p


def test_no_legacy_module_imports_in_scripts():
    """Scripts must not import from fit_partial_credit_model /
    fit_credit_model / fit_ordered_logit_model. They route through a Model
    subclass instead."""
    offenders = []
    for path in _iter_py_files(_scripts_dir):
        src = path.read_text()
        if _LEGACY_MODULES_RE.search(src):
            # Exception: Ukraine_interim_analyses_with_HMC.py keeps one
            # legacy import for fit_interim_MC_of_posterior_xz's
            # fitting_method= callable. See its docstring.
            if path.name == 'Ukraine_interim_analyses_with_HMC.py':
                continue
            offenders.append(str(path.relative_to(_repo_root)))
    assert not offenders, (
        f"Scripts still import legacy free fitters: {offenders}"
    )


def test_no_legacy_module_imports_in_python_lib():
    """python/*.py files must not import from the legacy modules except via
    the allow-listed Model subclasses + fit_interim's MC default."""
    offenders = []
    for path in _iter_py_files(_python_dir):
        if path.name in _ALLOWED_LEGACY_IMPORTERS:
            continue
        src = path.read_text()
        if _LEGACY_MODULES_RE.search(src):
            offenders.append(str(path.relative_to(_repo_root)))
    assert not offenders, (
        f"python/ files leak legacy imports: {offenders}"
    )


def test_no_deleted_shim_references():
    """The step-4 back-compat shim signatures
    (fit_interim_importance_sampling_of_posterior_xz_from_x etc.) were
    removed in step 8. No script or test may reference them."""
    offenders = []
    for root in (_scripts_dir, _python_dir, _repo_root / 'test' / 'python'):
        for path in _iter_py_files(root):
            if path.name == 'test_no_legacy_imports.py':
                continue   # this file documents the deleted names
            src = path.read_text()
            # Strip docstrings + comments before checking.
            stripped = re.sub(r'"""[\s\S]*?"""', '', src)
            stripped = re.sub(r"'''[\s\S]*?'''", '', stripped)
            stripped = '\n'.join(
                ln for ln in stripped.splitlines()
                if not ln.lstrip().startswith('#')
            )
            if _DELETED_SHIM_RE.search(stripped):
                offenders.append(str(path.relative_to(_repo_root)))
    assert not offenders, (
        f"Deleted shim signatures still referenced: {offenders}"
    )
