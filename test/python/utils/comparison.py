"""Utilities for comparing Python outputs to R baselines."""
import numpy as np
import pandas as pd
from typing import Union, Tuple, Dict, Any


def compare_dataframes(
    df_python: pd.DataFrame,
    df_r: pd.DataFrame,
    rtol: float = 1e-5,
    atol: float = 1e-8,
    check_dtypes: bool = True,
    check_column_order: bool = False
) -> Tuple[bool, str]:
    """
    Compare Python and R DataFrames for equivalence.
    
    Parameters
    ----------
    df_python : pd.DataFrame
        DataFrame from Python implementation
    df_r : pd.DataFrame
        DataFrame from R baseline
    rtol : float
        Relative tolerance for numerical comparisons
    atol : float
        Absolute tolerance for numerical comparisons
    check_dtypes : bool
        Whether to check data types match
    check_column_order : bool
        Whether column order must match
        
    Returns
    -------
    matches : bool
        True if DataFrames are equivalent
    message : str
        Description of differences if not equivalent
    """
    # Shape check
    if df_python.shape != df_r.shape:
        return False, f"Shape mismatch: Python {df_python.shape} vs R {df_r.shape}"
    
    # Column check
    python_cols = set(df_python.columns)
    r_cols = set(df_r.columns)
    if python_cols != r_cols:
        missing = r_cols - python_cols
        extra = python_cols - r_cols
        return False, f"Column mismatch. Missing: {missing}, Extra: {extra}"
    
    # Align columns if order doesn't matter
    if not check_column_order:
        df_python = df_python[df_r.columns]
    
    # Value comparison
    for col in df_r.columns:
        if pd.api.types.is_numeric_dtype(df_r[col]):
            # Numerical comparison with tolerance
            try:
                if not np.allclose(
                    df_python[col].values, 
                    df_r[col].values, 
                    rtol=rtol, 
                    atol=atol, 
                    equal_nan=True
                ):
                    max_diff = np.max(np.abs(df_python[col].values - df_r[col].values))
                    return False, f"Numerical mismatch in column '{col}' (max diff: {max_diff})"
            except (TypeError, ValueError) as e:
                return False, f"Error comparing column '{col}': {str(e)}"
        else:
            # Exact comparison for non-numerical (handle NaN)
            python_filled = df_python[col].fillna('__PANDAS_NA__')
            r_filled = df_r[col].fillna('__PANDAS_NA__')
            if not (python_filled.values == r_filled.values).all():
                return False, f"Value mismatch in column '{col}'"
    
    return True, "DataFrames match"


def compare_scalars(
    value_python: float,
    value_r: float,
    rtol: float = 1e-5,
    atol: float = 1e-8,
    name: str = "value"
) -> Tuple[bool, str]:
    """
    Compare scalar values with tolerance.
    
    Parameters
    ----------
    value_python : float
        Value from Python implementation
    value_r : float
        Value from R baseline
    rtol : float
        Relative tolerance
    atol : float
        Absolute tolerance
    name : str
        Name of the value being compared (for messages)
        
    Returns
    -------
    matches : bool
        True if values match within tolerance
    message : str
        Description of comparison result
    """
    if np.isclose(value_python, value_r, rtol=rtol, atol=atol):
        return True, f"{name} matches"
    else:
        diff = abs(value_python - value_r)
        rel_diff = diff / abs(value_r) if value_r != 0 else float('inf')
        return False, (
            f"{name} mismatch: Python={value_python:.6f}, R={value_r:.6f}, "
            f"diff={diff:.6f}, rel_diff={rel_diff:.2%}"
        )


def compare_posteriors(
    posterior_python: Dict[str, Dict[str, float]],
    posterior_r: Dict[str, Dict[str, float]],
    mean_rtol: float = 0.03,  # 3% per user spec
    rhat_atol: float = 0.005  # 0.005 per user spec
) -> Tuple[bool, str]:
    """
    Compare MCMC posterior outputs.
    
    Per user specification:
    - Posterior means: within 3% relative difference
    - Rhat: within 0.005 absolute difference
    
    Parameters
    ----------
    posterior_python : dict
        Python posterior with 'means' and 'rhat' dicts
    posterior_r : dict
        R posterior with 'means' and 'rhat' dicts
    mean_rtol : float
        Relative tolerance for means (default 0.03 = 3%)
    rhat_atol : float
        Absolute tolerance for Rhat (default 0.005)
        
    Returns
    -------
    matches : bool
        True if posteriors match within tolerances
    message : str
        Description of differences if not matched
    """
    messages = []
    
    # Check means
    r_means = posterior_r.get('means', {})
    py_means = posterior_python.get('means', {})
    
    for param, mean_r in r_means.items():
        if param not in py_means:
            return False, f"Missing parameter in Python posterior: {param}"
        mean_py = py_means[param]
        rel_diff = abs(mean_py - mean_r) / abs(mean_r) if mean_r != 0 else float('inf')
        if rel_diff > mean_rtol:
            messages.append(
                f"{param} mean: rel_diff={rel_diff:.2%} > {mean_rtol:.2%} "
                f"(Python={mean_py:.6f}, R={mean_r:.6f})"
            )
    
    # Check Rhat
    r_rhat = posterior_r.get('rhat', {})
    py_rhat = posterior_python.get('rhat', {})
    
    for param, rhat_r in r_rhat.items():
        if param not in py_rhat:
            return False, f"Missing Rhat for parameter: {param}"
        rhat_py = py_rhat[param]
        abs_diff = abs(rhat_py - rhat_r)
        if abs_diff > rhat_atol:
            messages.append(
                f"{param} Rhat: abs_diff={abs_diff:.4f} > {rhat_atol:.4f} "
                f"(Python={rhat_py:.4f}, R={rhat_r:.4f})"
            )
    
    if messages:
        return False, "; ".join(messages)
    return True, "Posteriors match within tolerance"


def compare_arrays(
    arr_python: np.ndarray,
    arr_r: np.ndarray,
    rtol: float = 1e-5,
    atol: float = 1e-8,
    name: str = "array"
) -> Tuple[bool, str]:
    """
    Compare NumPy arrays with tolerance.
    
    Parameters
    ----------
    arr_python : np.ndarray
        Array from Python implementation
    arr_r : np.ndarray
        Array from R baseline
    rtol : float
        Relative tolerance
    atol : float
        Absolute tolerance
    name : str
        Name of the array (for messages)
        
    Returns
    -------
    matches : bool
        True if arrays match within tolerance
    message : str
        Description of comparison result
    """
    # Shape check
    if arr_python.shape != arr_r.shape:
        return False, f"{name} shape mismatch: Python {arr_python.shape} vs R {arr_r.shape}"
    
    # Value comparison
    if np.allclose(arr_python, arr_r, rtol=rtol, atol=atol, equal_nan=True):
        return True, f"{name} matches"
    else:
        max_diff = np.max(np.abs(arr_python - arr_r))
        max_rel_diff = np.max(np.abs(arr_python - arr_r) / np.abs(arr_r + 1e-10))
        return False, (
            f"{name} mismatch: max_abs_diff={max_diff:.6e}, "
            f"max_rel_diff={max_rel_diff:.2%}"
        )
