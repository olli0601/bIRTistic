"""Web data extraction for bIRTistic applications — downloading and reshaping raw
public microdata that is not shipped with the repo. Sits next to ``data_loading.py``
(which reads local/common-format files); this module fetches from the web first.

Currently covers:
  * PISA cognitive microdata (OECD), §3.11 — download the student cognitive .sav files
    and build the common-country x common-Math-item extract for cycles 2015/2018/2022.
  * ICRC community-MHPSS supplement (Andersen et al. 2022), §3.10 — the beneficiary-level
    Data_Sheet is a published Frontiers/PMC supplement (via Wayback; www.oecd/PMC bot-block).

All large raw downloads land in a scratch/cache dir; only compact derived extracts are
meant to be persisted alongside the other application data.
"""
import os
import re
import json
import subprocess
import zipfile
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd

UA = "Mozilla/5.0"

# ---------------------------------------------------------------------------
# PISA (OECD) — §3.11
# ---------------------------------------------------------------------------
# Student cognitive (SPSS) files, same coding scheme 2015+. 2012 is paper-era
# Items are matched across cycles by CONTENT ID (strip mode prefix + variant suffix:
# PM033Q01 == CM033Q01S -> M033Q01), so the paper-era 2012 links to the computer-era
# 2015/18/22. 2012 is paper (mode confound at 2012->2015; OECD applies a mode adjustment we
# do not) and comes as fixed-width TXT + SPSS control file, both via the Wayback mirror
# (www.oecd.org is 403). rho is baseline-anchored: baseline=2012, each later cycle an interim.
PISA_COG_URLS: Dict[int, str] = {
    2015: "https://webfs.oecd.org/pisa/PUF_SPSS_COMBINED_CMB_STU_COG.zip",   # 172 MB .sav
    2018: "https://webfs.oecd.org/pisa2018/SPSS_STU_COG.zip",                # 478 MB .sav
    2022: "https://webfs.oecd.org/pisa2022/STU_COG_SPSS.zip",                # 481 MB .sav
}
PISA_2012_TXT = "https://web.archive.org/web/2id_/http://www.oecd.org/pisa/pisaproducts/INT_COG12_S_DEC03.zip"
PISA_2012_CTL = ("https://web.archive.org/web/2id_/https://www.oecd.org/content/dam/oecd/en/data/"
                 "datasets/pisa/pisa-2012-datasets/main-survey/sas-and-spss-control-files/"
                 "SPSS%20syntax%20to%20read%20in%20SCORED%20cognitive%20item%20response%20data%20file.txt")


def _curl(url: str, out: str, timeout: int = 1200) -> str:
    subprocess.run(["curl", "-sL", "-A", UA, "-m", str(timeout), url, "-o", out], check=True)
    return out


def pisa_content_id(colname: str) -> Optional[str]:
    """Cross-cycle content id: strip mode prefix (P/C) + variant suffix -> M<unit>Q<qn>."""
    m = re.match(r"[CP]?(M\d+Q\d+)", str(colname))
    return m.group(1) if m else None


def download_pisa_cognitive(cycle: int, cache_dir: str) -> str:
    """Download + unzip a computer-era cycle's student cognitive .sav; return its path."""
    os.makedirs(cache_dir, exist_ok=True)
    zp = os.path.join(cache_dir, f"cog{cycle}.zip")
    if not (os.path.exists(zp) and os.path.getsize(zp) > 10_000_000):
        _curl(PISA_COG_URLS[cycle], zp)
    dd = os.path.join(cache_dir, f"s{cycle}")
    with zipfile.ZipFile(zp) as z:
        z.extractall(dd)
    savs = [str(p) for p in Path(dd).rglob("*") if p.suffix.lower() == ".sav"]
    if not savs:
        raise FileNotFoundError(f"no .sav in {dd}")
    return max(savs, key=os.path.getsize)


def _pisa_2012_colspecs(cache_dir: str):
    """Parse the 2012 scored-cognitive SPSS control file into (names, colspecs) for CNT +
    Math (PM..Q) items. Downloads TXT (fixed-width) + control via Wayback if absent."""
    txtzip = os.path.join(cache_dir, "cog2012.zip")
    if not (os.path.exists(txtzip) and os.path.getsize(txtzip) > 1_000_000):
        _curl(PISA_2012_TXT, txtzip)
    with zipfile.ZipFile(txtzip) as z:
        z.extractall(cache_dir)
    txtpath = os.path.join(cache_dir, "INT_COG12_S_DEC03.txt")
    ctl = os.path.join(cache_dir, "ctl2012_spss.txt")
    if not os.path.exists(ctl):
        _curl(PISA_2012_CTL, ctl)
    spec = open(ctl, encoding="latin-1").read()
    found = re.findall(r"([A-Za-z]\w+)\s+(\d+)\s*-\s*(\d+)|([A-Za-z]\w+)\s+(\d+)\b(?!\s*-)", spec)
    names, colspecs = [], []
    for a, s, e, b, p in found:
        nm = a or b
        if nm == "CNT" or re.match(r"^PM\d+Q\d+", nm):
            lo = (int(s) if a else int(p)) - 1
            hi = int(e) if a else int(p)
            names.append(nm); colspecs.append((lo, hi))
    return txtpath, names, colspecs


def _read_pisa_2012(cache_dir: str) -> "pd.DataFrame":
    """Fixed-width read of the 2012 scored cognitive file -> long (CNT, item=content-id, y)."""
    txtpath, names, colspecs = _pisa_2012_colspecs(cache_dir)
    df = pd.read_fwf(txtpath, colspecs=colspecs, names=names, dtype=str)
    df["CNT"] = df["CNT"].str.strip()
    item_cols = [n for n in names if n != "CNT"]
    ren = {n: pisa_content_id(n) for n in item_cols}
    df = df.rename(columns=ren)
    df["_sid"] = np.arange(len(df))
    long = df.melt(id_vars=["CNT", "_sid"], value_vars=[ren[n] for n in item_cols],
                   var_name="item", value_name="y")
    long["y"] = pd.to_numeric(long["y"], errors="coerce")
    return long


def _scored_content_map(cols, meta):
    """content-id -> scored (credit-labelled) column in a computer-era .sav. Prefers the
    COMPUTER scored variant 'C'+cid+'S' (e.g. CM033Q01S) over paper (PM..), whose column is
    mostly-missing for the majority computer cohort."""
    from collections import defaultdict
    cand = defaultdict(list)
    for c in cols:
        cid = pisa_content_id(c)
        if cid and "credit" in " ".join(str(v).lower() for v in meta.variable_value_labels.get(c, {}).values()):
            cand[cid].append(c)
    out = {}
    for cid, cs in cand.items():
        out[cid] = min(cs, key=lambda c: (c != "C" + cid + "S",
                                          not (c.startswith("C") and c.endswith("S")),
                                          not c.startswith("C"), len(c)))
    return out


def pisa_common_math_4cyc(cache_dir: str, cycles=(2012, 2015, 2018, 2022)):
    """Content-id Math item + country intersection across the cycles (2012 paper + .sav)."""
    import pyreadstat
    items, countries = {}, {}
    for yr in cycles:
        if yr == 2012:
            _, names, _ = _pisa_2012_colspecs(cache_dir)
            items[yr] = {pisa_content_id(n) for n in names if n != "CNT"}
            df = _read_pisa_2012(cache_dir)
            countries[yr] = set(df["CNT"].dropna().unique())
        else:
            p = download_pisa_cognitive(yr, cache_dir)
            _, meta = pyreadstat.read_sav(p, metadataonly=True)
            items[yr] = set(_scored_content_map(meta.column_names, meta))
            dfc, _ = pyreadstat.read_sav(p, usecols=["CNT"])
            countries[yr] = set(dfc["CNT"].dropna().unique())
    return sorted(set.intersection(*items.values())), sorted(set.intersection(*countries.values()))


def build_pisa_math_extract(cache_dir: str, out_dir: str, cycles=(2012, 2015, 2018, 2022)):
    """4-cycle content-id-matched extract. Per cycle -> long parquet (CNT, id, item, y in
    {0,1,2}, cycle) over the common Math bank; skips/NR (5-9) -> NaN. baseline=2012."""
    import pyreadstat
    os.makedirs(out_dir, exist_ok=True)
    items, countries = pisa_common_math_4cyc(cache_dir, cycles)
    json.dump({"items": items, "countries": countries}, open(f"{out_dir}/common_math.json", "w"))
    ITEMS, COUNTRIES = set(items), set(countries)
    for yr in cycles:
        if yr == 2012:
            long = _read_pisa_2012(cache_dir)
            long = long[long["CNT"].isin(COUNTRIES) & long["item"].isin(ITEMS)].copy()
            long = long.rename(columns={"_sid": "CNTSTUID"})
        else:
            p = download_pisa_cognitive(yr, cache_dir)
            _, meta = pyreadstat.read_sav(p, metadataonly=True)
            cols = meta.column_names
            idcol = next((c for c in ("CNTSTUID", "CNTSTID") if c in cols), None)
            m = {cid: c for cid, c in _scored_content_map(cols, meta).items() if cid in ITEMS}
            df, _ = pyreadstat.read_sav(p, usecols=["CNT"] + ([idcol] if idcol else []) + list(m.values()))
            df = df[df["CNT"].isin(COUNTRIES)].rename(columns={c: cid for cid, c in m.items()})
            long = df.melt(id_vars=["CNT"] + ([idcol] if idcol else []), value_vars=sorted(m),
                           var_name="item", value_name="y")
            if idcol and idcol != "CNTSTUID":
                long = long.rename(columns={idcol: "CNTSTUID"})
        long["y"] = pd.to_numeric(long["y"], errors="coerce")
        long = long[long["y"] <= 4].dropna(subset=["y"])
        long["y"] = long["y"].astype("int8"); long["cycle"] = yr
        long.to_parquet(f"{out_dir}/pisa_math_{yr}.parquet", index=False)
        print(f"{yr}: {long['CNT'].nunique()} countries, {long['item'].nunique()} items -> {len(long)} responses")
    print(f"common: {len(items)} math items x {len(countries)} countries -> {out_dir}")


# ---------------------------------------------------------------------------
# ICRC community-MHPSS supplement (Andersen et al. 2022) — §3.10
# ---------------------------------------------------------------------------
# The beneficiary-level Data Sheet is a published Frontiers/PMC supplement. Direct PMC
# and www.oecd/Frontiers endpoints are bot-blocked; the Wayback raw mirror serves them.
ICRC_ANDERSEN_PMC = "PMC8995431"   # Front. Public Health 10:815222
ICRC_ANDERSEN_FILES = {            # PMC bin filenames (Data sheet + population table)
    "Data_Sheet_1.xlsx": f"https://pmc.ncbi.nlm.nih.gov/articles/instance/8995431/bin/Data_Sheet_1.xlsx",
    "Table_1.XLSX": f"https://pmc.ncbi.nlm.nih.gov/articles/instance/8995431/bin/Table_1.XLSX",
}


def download_icrc_andersen(dest_dir: str) -> Dict[str, str]:
    """Fetch the Andersen 2022 supplement files via the Wayback raw mirror (PMC bot-blocks
    curl). Returns {name: path}. The 'Data' sheet is beneficiary-level; 'Coding' is the
    codebook. Manual browser download is the fallback if Wayback rate-limits."""
    os.makedirs(dest_dir, exist_ok=True)
    got = {}
    for name, pmc_url in ICRC_ANDERSEN_FILES.items():
        out = os.path.join(dest_dir, f"Andersen_{name}")
        wb = f"https://web.archive.org/web/2id_/{pmc_url}"
        try:
            _curl(wb, out, timeout=120)
            if os.path.getsize(out) > 100_000 and open(out, "rb").read(2) == b"PK":
                got[name] = out
        except subprocess.CalledProcessError:
            pass
    return got


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default="/tmp/pisa_cache")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    build_pisa_math_extract(a.cache, a.out)
