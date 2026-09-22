"""Generic FEDERATED diagnostic figures for the SDY269 head-to-head (§17.5). Each endpoint (rho) is
served by a DIFFERENT trained amortiser in a DIFFERENT deploy directory; this module resolves each
declared rho to its directory via a MANIFEST and assembles the cross-endpoint diagnostics — so the
figures load pkls/csvs from several federated amortisers at once. Nothing here is SDY269-specific
except the manifest: point MANIFEST at any app's per-rho deploy dirs (as per its rho_spec).

Per deploy dir it reads the driver's standard outputs:
  pcm_1_interim_pps_RAGD_pit_by_item.csv   -> PIT-KS / marg-KS / coverage (calibration)
  pcm_1_interim_pps_RAGD_headft_coverage.csv -> per-interim (across-n) pit_ks/marg_ks
  pcm_1_interim_pps_RAGD_pps_by_item.csv   -> PPS trajectory over interims
Figures (into OUT): headhead_amortiser_calibration, headhead_amortiser_pit_vs_n, headhead_amortiser_pps."""
# ---- boilerplate ----

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
import numpy as np, pandas as pd

SB = os.environ.get('IMMPORT_SB', "/Users/or105/sandbox/bIRTistic")
OUT = f"{SB}/py-immport-SDY269"; pfx = "pcm_1_interim"
# MANIFEST: rho_label -> (deploy dir, registry instance, warp) — the federated resolution per rho_spec
MANIFEST = {
    'LAIV_spr':          (f"{SB}/py-immport-SDY269-LAIV_sprdeploy_260922",  'S1 rate',  'logit'),
    'TIV_spr':           (f"{SB}/py-immport-SDY269-TIV_sprdeploy_260922",   'S1 rate',  'logit'),
    'LAIV_gmfr':         (f"{SB}/py-immport-SDY269-LAIV_gmfrdeploy_260922", 'S2 fold',  'log2'),
    'TIV_gmfr':          (f"{SB}/py-immport-SDY269-TIV_gmfrdeploy_260922",  'S2 fold',  'log2'),
    'TIV-LAIV_spr_diff': (f"{SB}/py-immport-SDY269-S5deploy_260922",        'S5 diff',  'none'),
}
RHO_ORDER = list(MANIFEST)


def _read(dirpath, name):
    p = f"{dirpath}/{pfx}_pps_RAGD_{name}.csv"
    return pd.read_csv(p) if os.path.exists(p) else None


def load_all():
    """Load each rho's federated deploy outputs, tagged by rho_label + instance."""
    cal, traj, pps = [], [], []
    for rho, (d, inst, warp) in MANIFEST.items():
        pi = _read(d, 'pit_by_item')
        if pi is not None:
            r = pi.iloc[0]; cal.append(dict(rho=rho, instance=inst, warp=warp,
                                            pit_ks=r.pit_ks, marg_ks=r.marg_ks, cov5=r.cov5, cov95=r.cov95))
        cv = _read(d, 'headft_coverage')
        if cv is not None:
            h = cv[cv.config != 'baseline(ragged)'].copy(); h['rho'] = rho; h['instance'] = inst
            traj.append(h[['rho', 'instance', 'interim_id', 'n', 'pit_ks', 'marg_ks']])
        pp = _read(d, 'pps_by_item')
        if pp is not None:
            pp = pp.copy(); pp['rho'] = rho; pp['instance'] = inst; pps.append(pp)
    return (pd.DataFrame(cal),
            pd.concat(traj, ignore_index=True) if traj else pd.DataFrame(),
            pd.concat(pps, ignore_index=True) if pps else pd.DataFrame())


def _pal(keys):
    import ggsci, matplotlib.colors as mcolors
    pal = [mcolors.to_hex(c) for c in ggsci.pal_futurama("planetexpress")(12)]
    return {k: pal[i % len(pal)] for i, k in enumerate(keys)}


def fig_calibration(cal):
    from plotnine import (ggplot, aes, geom_col, geom_hline, facet_wrap, theme_bw, theme, labs,
                          scale_fill_manual, element_text, coord_flip)
    m = cal.melt(id_vars=['rho', 'instance'], value_vars=['pit_ks', 'marg_ks'], var_name='metric', value_name='ks')
    m['rho_f'] = pd.Categorical(m['rho'], categories=RHO_ORDER[::-1], ordered=True)
    cd = _pal(cal['instance'].unique())
    p = (ggplot(m, aes(x='rho_f', y='ks', fill='instance'))
         + geom_col(show_legend=True) + coord_flip()
         + geom_hline(yintercept=0.10, linetype='dashed', color='#555555')     # ~good-calibration guide
         + facet_wrap('metric', ncol=2)
         + scale_fill_manual(values=cd, name='amortiser')
         + theme_bw() + theme(figure_size=(10, 4), strip_text=element_text(size=9))
         + labs(x='endpoint (federated amortiser)', y='KS vs SVI (dashed = 0.10)',
                title='SDY269 head-to-head: federated amortiser calibration (each rho -> its own net)'))
    p.save(f"{OUT}/headhead_amortiser_calibration.pdf", verbose=False)
    p.save(f"{OUT}/headhead_amortiser_calibration.png", dpi=120, verbose=False)


def fig_pit_vs_n(traj):
    from plotnine import (ggplot, aes, geom_line, geom_point, geom_hline, theme_bw, theme, labs,
                          scale_color_manual, element_text)
    t = traj.copy(); t['rho_f'] = pd.Categorical(t['rho'], categories=RHO_ORDER, ordered=True)
    cd = _pal(RHO_ORDER)
    p = (ggplot(t, aes(x='n', y='pit_ks', color='rho_f'))
         + geom_hline(yintercept=0.10, linetype='dashed', color='#999999')
         + geom_line() + geom_point(size=1.3)
         + scale_color_manual(values=cd, name='endpoint')
         + theme_bw() + theme(figure_size=(9, 4.5), strip_text=element_text(size=9))
         + labs(x='interim (participants per arm)', y='PIT-KS vs SVI',
                title='SDY269 head-to-head: federated amortiser calibration across n'))
    p.save(f"{OUT}/headhead_amortiser_pit_vs_n.pdf", verbose=False)
    p.save(f"{OUT}/headhead_amortiser_pit_vs_n.png", dpi=120, verbose=False)


def fig_pps(pps):
    from plotnine import (ggplot, aes, geom_line, geom_point, theme_bw, theme, labs,
                          scale_color_manual, element_text)
    t = pps.copy(); t['rho_f'] = pd.Categorical(t['rho'], categories=RHO_ORDER, ordered=True)
    cd = _pal(RHO_ORDER)
    p = (ggplot(t, aes(x='n', y='pps', color='rho_f'))
         + geom_line() + geom_point(size=1.3)
         + scale_color_manual(values=cd, name='endpoint')
         + theme_bw() + theme(figure_size=(9, 4.5), strip_text=element_text(size=9))
         + labs(x='interim (participants per arm)', y='PPS = P(endpoint meets H1 | x)',
                title=f"SDY269 head-to-head: amortised PPS trajectory (eta0={t['eta0'].iloc[0] if 'eta0' in t else '?'})"))
    p.save(f"{OUT}/headhead_amortiser_pps.pdf", verbose=False)
    p.save(f"{OUT}/headhead_amortiser_pps.png", dpi=120, verbose=False)


if __name__ == "__main__":
    cal, traj, pps = load_all()
    cal.to_csv(f"{OUT}/headhead_amortiser_calibration.csv", index=False)
    pd.set_option('display.width', 160)
    print("=== federated amortiser calibration (per rho -> its net) ===")
    print(cal.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    if not cal.empty: fig_calibration(cal)
    if not traj.empty: fig_pit_vs_n(traj)
    if not pps.empty: fig_pps(pps)
    print(f"\nfigures -> {OUT}/headhead_amortiser_*.pdf")
