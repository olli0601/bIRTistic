#!/usr/bin/env python3
"""Composite results figure 4 for the SDY312 influenza-HAI amortised-PPS application.

A4 (portrait). Panel (a) — top-left third-width x quarter-height — is a drawn-from-source
schematic of the immunology + the three per-strain effect sizes (previous infection / vaccination
-> antibodies -> cross-specific vs monospecific -> H1N1 / H3N2 -> SPR / SCR / GMFR). The remaining
panels (b, c, ...) are added as the figure is built out.

Run:  pixi run python paper_adm/imm_sdy312_figure_4.py
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, Ellipse, Polygon, FancyBboxPatch

A4_W, A4_H = 8.27, 11.69          # A4 portrait (inches)
DPI = 200
GREY = '#9aa0a6'
AB_G, AB_B = '#4caf7d', '#4a90d9'   # antibody green / blue
VIR = '#bfe3ea'                      # virus body


# ---- schematic primitives (data coords 0..100 x 0..100) --------------------------------
def _arrow(ax, x0, y0, x1, y1, lw=1.4, scale=11):
    ax.annotate('', xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle='-|>', color='black', lw=lw, mutation_scale=scale,
                                shrinkA=0, shrinkB=0))


def _inhibit(ax, x0, y0, x1, y1, lw=1.4, bar=5):
    """T-bar (inhibition) end."""
    import numpy as np
    ax.plot([x0, x1], [y0, y1], color='black', lw=lw, solid_capstyle='round')
    dx, dy = x1 - x0, y1 - y0
    n = (dx * dx + dy * dy) ** 0.5
    px, py = -dy / n, dx / n
    ax.plot([x1 - px * bar, x1 + px * bar], [y1 - py * bar, y1 + py * bar], color='black', lw=lw, solid_capstyle='round')


def _antibody(ax, cx, cy, s=3.2, col=AB_G, lw=1.3):
    """A simple Y-shaped antibody glyph centred at (cx,cy)."""
    ax.plot([cx, cx], [cy - s, cy], color=col, lw=lw, solid_capstyle='round')          # stem
    ax.plot([cx, cx - s * 0.9], [cy, cy + s], color=col, lw=lw, solid_capstyle='round')  # left arm
    ax.plot([cx, cx + s * 0.9], [cy, cy + s], color=col, lw=lw, solid_capstyle='round')  # right arm


def _person(ax, cx, cy):
    """Simplified grey human silhouette (head + shoulders/torso)."""
    ax.add_patch(Circle((cx, cy + 9), 4.6, color=GREY, zorder=2))                       # head
    torso = FancyBboxPatch((cx - 7.5, cy - 11), 15, 17, boxstyle='round,pad=0,rounding_size=4',
                           color=GREY, zorder=2)
    ax.add_patch(torso)
    ax.add_patch(Ellipse((cx, cy + 4.5), 20, 7, color=GREY, zorder=2))                  # shoulders


def _virus(ax, cx, cy, r=5.0, label='', spikes=('#8e5fb0', '#4a90d9')):
    """Virus particle: body circle + radial spikes + label."""
    import numpy as np
    for i, th in enumerate(np.linspace(0, 2 * np.pi, 16, endpoint=False)):
        c = spikes[i % len(spikes)]
        x0, y0 = cx + r * np.cos(th), cy + r * np.sin(th)
        x1, y1 = cx + (r + 2.1) * np.cos(th), cy + (r + 2.1) * np.sin(th)
        ax.plot([x0, x1], [y0, y1], color=c, lw=1.2, solid_capstyle='round', zorder=1)
        ax.add_patch(Circle((x1, y1), 0.9, color=c, zorder=1))
    ax.add_patch(Circle((cx, cy), r, facecolor=VIR, edgecolor='#7fb8c4', lw=0.8, zorder=2))
    ax.text(cx, cy, label, ha='center', va='center', fontsize=5.0, fontweight='bold', zorder=3)


def panel_a(ax):
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis('off')
    TXT = dict(ha='center', va='center', fontsize=6.2)
    # previous infection -> person ; vaccination -> person
    ax.text(50, 97, 'Previous infection', **TXT)
    _arrow(ax, 50, 94, 50, 86)
    ax.text(12, 80, 'vaccination', ha='center', va='center', fontsize=6.2)
    _arrow(ax, 24, 78, 40, 76)
    _person(ax, 50, 74)
    # person -> antibodies
    _arrow(ax, 50, 60, 50, 52)
    ax.text(50, 49, 'Antibodies', **TXT)
    for (dx, col) in [(-9, AB_B), (-3, AB_G), (3, AB_B), (9, AB_G)]:
        _antibody(ax, 50 + dx, 43, s=2.6, col=col)
    # split: cross-specific (left) / monospecific (centre)
    _arrow(ax, 44, 40, 26, 33)
    _arrow(ax, 54, 40, 60, 33)
    ax.text(22, 31, 'Cross-specific\nantibodies', ha='center', va='center', fontsize=5.8)
    ax.text(62, 31, 'Monospecific\nantibodies', ha='center', va='center', fontsize=5.8)
    # cross-specific antibody inhibits BOTH strains; monospecific inhibits each separately
    _antibody(ax, 22, 23, s=2.8, col=AB_B)
    _inhibit(ax, 20, 20, 12, 12); _inhibit(ax, 24, 20, 33, 12)
    _antibody(ax, 56, 23, s=2.8, col=AB_G)
    _inhibit(ax, 56, 20, 54, 12)
    _antibody(ax, 70, 23, s=2.8, col=AB_B)
    _inhibit(ax, 70, 20, 72, 12)
    _virus(ax, 12, 8, label='H1N1'); _virus(ax, 33, 8, label='H3N2')
    _virus(ax, 54, 8, label='H1N1'); _virus(ax, 72, 8, label='H3N2')
    # the three effect sizes, pointed to from the monospecific strains
    _arrow(ax, 78, 16, 86, 26)
    ax.text(88, 44, '3 effect sizes per strain:', ha='left', va='center', fontsize=6.0, fontweight='bold')
    for i, t in enumerate(['Absolute level of\nprotection (SPR)',
                           'Any immunogenicity\n(SCR, seroconversion)',
                           'Magnitude of\nimmunogenicity (GMFR)']):
        ax.text(89, 37 - i * 10.5, '• ' + t, ha='left', va='top', fontsize=5.4)


def main():
    fig = plt.figure(figsize=(A4_W, A4_H))
    # panel (a): top-left, 1/3 width x 1/4 height (small margins)
    axA = fig.add_axes([0.015, 0.735, 0.333, 0.25]); panel_a(axA)
    fig.text(0.012, 0.992, 'a', fontsize=15, fontweight='bold', va='top', ha='left')
    out = os.path.join(os.environ.get('PAPER_ADM_OUT', os.path.dirname(os.path.abspath(__file__))),
                       'imm_sdy312_figure_4')
    fig.savefig(out + '.pdf', dpi=DPI)
    fig.savefig(out + '.png', dpi=DPI)
    plt.close(fig)
    print(f"Saved figure 4 -> {out}.pdf / .png")


if __name__ == '__main__':
    main()
