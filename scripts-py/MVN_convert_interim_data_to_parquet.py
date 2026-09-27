"""Task 4: convert legacy interim-data pickles to parquet-backed dirs.

Two phases (env ``CONVERT_DELETE_PKL``):
  0 (default): WRITE the parquet dir beside each ``mvn_J{J}_interim_data.pkl`` and
      verify (zi bit-exact, mu_draws close). The ``.pkl`` is KEPT, so
      ``load_interim_data`` still reads it (parquet is dormant) -- zero risk to
      running readers.
  1: after all readers are patched to ``load_interim_data`` and no reader is
      mid-run, DELETE each ``.pkl`` whose verified parquet dir exists, flipping
      everything to parquet (~5x smaller).
"""
import os
import sys
import glob

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
from amortiser_io import save_interim_data, load_interim_data, _dir_for  # noqa: E402

SB = '/Users/or105/sandbox/bIRTistic'
DIRS = [os.path.join(SB, 'py-mvn-interim-simulations-260609')]
DELETE = os.environ.get('CONVERT_DELETE_PKL', '0') == '1'


def _verify(pkl_path):
    """Load legacy pkl + the parquet dir, confirm zi bit-exact + mu close."""
    src = pd.read_pickle(pkl_path)
    dst = load_interim_data(_dir_for(pkl_path))          # force dir read
    if set(src) != set(dst):
        return False, "interim-id set mismatch"
    for k in src:
        if not src[k]['zi'].equals(dst[k]['zi']):
            # numeric-only f4 compare (skips datetime cols); tolerates f8->f4 downcast
            a = src[k]['zi'].select_dtypes('number').to_numpy('float32')
            b = dst[k]['zi'].select_dtypes('number').to_numpy('float32')
            if not np.array_equal(a, b, equal_nan=True):
                return False, f"zi mismatch at interim {k}"
        if not np.allclose(src[k]['mu_draws'], dst[k]['mu_draws'], equal_nan=True):
            return False, f"mu_draws mismatch at interim {k}"
    return True, "ok"


def main():
    for d in DIRS:
        for pkl in sorted(glob.glob(os.path.join(d, 'mvn_J*_interim_data.pkl'))):
            outdir = _dir_for(pkl)
            if not DELETE:
                b0 = os.path.getsize(pkl) / 1e9
                if not os.path.isdir(outdir):
                    save_interim_data(pd.read_pickle(pkl), pkl)
                ok, msg = _verify(pkl)
                b1 = sum(os.path.getsize(os.path.join(r, f))
                         for r, _, fs in os.walk(outdir) for f in fs) / 1e9
                print(f"WRITE {os.path.basename(pkl)}: {b0:.1f}G -> {b1:.1f}G "
                      f"({b0/b1:.2f}x) verify={msg}", flush=True)
            else:
                if os.path.isdir(outdir):
                    ok, msg = _verify(pkl)
                    if ok:
                        os.remove(pkl)
                        print(f"DELETE {os.path.basename(pkl)} (verified {msg})", flush=True)
                    else:
                        print(f"KEEP {os.path.basename(pkl)}: verify FAILED ({msg})", flush=True)
                else:
                    print(f"SKIP {os.path.basename(pkl)}: no parquet dir", flush=True)
    print("DONE convert", flush=True)


if __name__ == '__main__':
    main()
