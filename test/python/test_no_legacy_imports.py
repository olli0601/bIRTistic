"""
Forbid stale references to deleted back-compat shims and to the
``fit_*_model`` legacy free-function fitters in user-facing code
(scripts + algorithm library).

After the model_*.py inlining all three legacy modules
(fit_partial_credit_model.py, fit_credit_model.py,
fit_ordered_logit_model.py) are deleted -- no code in the repo should
import them.
"""

import re
import sys
from pathlib import Path

import pytest

_repo_root = Path(__file__).resolve().parents[2]
_python_dir = _repo_root / 'python'
_scripts_dir = _repo_root / 'scripts-py'

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
    offenders = [
        str(path.relative_to(_repo_root))
        for path in _iter_py_files(_scripts_dir)
        if _LEGACY_MODULES_RE.search(path.read_text())
    ]
    assert not offenders, (
        f"Scripts still import legacy free fitters: {offenders}"
    )


def test_no_legacy_module_imports_in_python_lib():
    """python/*.py files must not import from any of the deleted legacy
    modules."""
    offenders = [
        str(path.relative_to(_repo_root))
        for path in _iter_py_files(_python_dir)
        if _LEGACY_MODULES_RE.search(path.read_text())
    ]
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
