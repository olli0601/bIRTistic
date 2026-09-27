"""Parquet-backed I/O for amortiser interim data.

The interim-data object is a dict keyed by ``interim_id``::

    {interim_id: {'zi': DataFrame, 'mu_draws': ndarray, 'dpi': DataFrame,
                  'n_obs': int, 'interim_m': int, 'interim_date': Timestamp,
                  'interim_month_year': str, ...}}

``zi`` is ~99% of the bytes (a wide table of ``ypred_*`` predictive-draw
columns). Columnar parquet + zstd compresses it ~5x (vs pickle f8), losslessly,
because each draw column is compressed on its own. Everything else per interim is
tiny and kept in a single ``meta.pkl``.

On-disk layout (a directory, sibling of the legacy ``*.pkl``)::

    mvn_J{J}_interim_data/
        meta.pkl                 # dict {interim_id: {<all keys except 'zi'>}}
        zi_i{interim_id}.parquet # the wide zi table for that interim

``load_interim_data`` transparently accepts either the new directory or a legacy
single ``.pkl`` file, so readers can migrate with a one-line swap and old data
keeps working.
"""
from __future__ import annotations

import os
import shutil

import pandas as pd

PARQUET_COMPRESSION = 'zstd'
PARQUET_LEVEL = 7


def save_df_parquet(df: pd.DataFrame, path: str) -> None:
    """Write a DataFrame to parquet with zstd (columnar, lossless)."""
    df.to_parquet(path, engine='pyarrow', compression=PARQUET_COMPRESSION,
                  compression_level=PARQUET_LEVEL, index=False)


def load_df_parquet(path: str) -> pd.DataFrame:
    return pd.read_parquet(path, engine='pyarrow')


def _dir_for(path: str) -> str:
    """Directory path for the new format given either the dir or a legacy .pkl."""
    if path.endswith('.pkl'):
        return path[:-len('.pkl')]
    return path


def save_interim_data(d: dict, path: str) -> str:
    """Persist an interim-data dict as a parquet-backed directory.

    ``path`` may be the target directory or the legacy ``*.pkl`` path (the
    ``.pkl`` suffix is stripped to name the directory). Returns the directory.
    """
    out = _dir_for(path)
    tmp = out + '.tmp'
    if os.path.exists(tmp):
        shutil.rmtree(tmp)
    os.makedirs(tmp, exist_ok=True)
    meta = {}
    for k, e in d.items():
        e = dict(e)
        zi = e.pop('zi', None)
        if zi is not None:
            save_df_parquet(zi, os.path.join(tmp, f'zi_i{int(k)}.parquet'))
        meta[k] = e                                    # everything except zi (small)
    pd.to_pickle(meta, os.path.join(tmp, 'meta.pkl'))
    if os.path.exists(out):
        shutil.rmtree(out)
    os.replace(tmp, out)
    return out


class InterimDataWriter:
    """Stream interim blocks to a parquet-backed dir, one at a time.

    Caps peak memory at ~one interim's ``zi`` (the wide table) instead of
    holding every interim in RAM before a single write. Usage::

        w = InterimDataWriter(pkl_path)
        for interim_id, entry in ...:
            w.add(interim_id, entry)     # writes zi parquet now, keeps only meta
        w.close()                        # writes meta.pkl, swaps dir into place
    """

    def __init__(self, path: str):
        self.out = _dir_for(path)
        self.tmp = self.out + '.tmp'
        if os.path.exists(self.tmp):
            shutil.rmtree(self.tmp)
        os.makedirs(self.tmp, exist_ok=True)
        self._meta = {}

    def add(self, interim_id, entry: dict) -> None:
        entry = dict(entry)
        zi = entry.pop('zi', None)
        if zi is not None:
            save_df_parquet(zi, os.path.join(self.tmp, f'zi_i{int(interim_id)}.parquet'))
        self._meta[interim_id] = entry                 # small: no zi held

    def close(self) -> str:
        pd.to_pickle(self._meta, os.path.join(self.tmp, 'meta.pkl'))
        if os.path.exists(self.out):
            shutil.rmtree(self.out)
        os.replace(self.tmp, self.out)
        return self.out


def load_interim_data(path: str, with_zi: bool = True) -> dict:
    """Load an interim-data dict from the new parquet dir or a legacy ``*.pkl``.

    Resolution order: an existing legacy ``.pkl`` file, else the parquet
    directory (either given directly or derived by stripping ``.pkl``).

    ``with_zi=False`` returns everything EXCEPT the wide ``zi`` tables (the
    memory hog): use it with :func:`load_interim_zi` to stream one interim's
    ``zi`` at a time and keep peak memory at ~one interim.
    """
    if path.endswith('.pkl') and os.path.isfile(path):
        d = pd.read_pickle(path)                       # legacy single pickle
        if not with_zi:
            for e in d.values():
                e.pop('zi', None)                      # already materialised; free it
        return d
    out = _dir_for(path)
    if not os.path.isdir(out):
        raise FileNotFoundError(f"no interim data at {path} (checked .pkl and {out}/)")
    meta = pd.read_pickle(os.path.join(out, 'meta.pkl'))
    d = {}
    for k, e in meta.items():
        e = dict(e)
        if with_zi:
            zp = os.path.join(out, f'zi_i{int(k)}.parquet')
            if os.path.isfile(zp):
                e['zi'] = load_df_parquet(zp)
        d[k] = e
    return d


def load_interim_zi(path: str, interim_id) -> pd.DataFrame:
    """Load a single interim's ``zi`` table (lazy for the parquet format).

    For the parquet dir this reads only ``zi_i{interim_id}.parquet``. For a
    legacy ``.pkl`` there are no per-interim files, so the whole pickle is read
    and the one ``zi`` returned (kept for back-compat; parquet is the lean path).
    """
    if path.endswith('.pkl') and os.path.isfile(path):
        return pd.read_pickle(path)[interim_id]['zi']
    out = _dir_for(path)
    zp = os.path.join(out, f'zi_i{int(interim_id)}.parquet')
    return load_df_parquet(zp)


def interim_data_exists(path: str) -> bool:
    """True if either the legacy ``.pkl`` or the new parquet dir is present."""
    if path.endswith('.pkl') and os.path.isfile(path):
        return True
    return os.path.isdir(_dir_for(path))
