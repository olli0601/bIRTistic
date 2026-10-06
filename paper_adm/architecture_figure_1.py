"""architecture_figure_1 — the federated PPS-amortiser architecture sketch (d2l-style).

ONE layout spec drives two deliverables:
  * architecture_figure_1.pdf / .png  — matplotlib, WHITE background (paper figure)
  * architecture_figure_1.html        — inline SVG, light theme, hover a box for its detail

Box TEXT is edited in `architecture_figure_1_text.md` (title / subtitle / line / hover per box); the box
POSITIONS, wrappers, arrows and the manifest ρ-cards live in this file (LAYOUT / WRAPPERS / ARROWS /
CARDS). Inline maths is written $...$ (TeX) in the text file: rendered as mathtext in the PDF and as
unicode + <tspan> sub/superscripts in the HTML.

    pixi run python paper_adm/architecture_figure_1.py         # -> paper_adm/architecture_figure_1.{pdf,png,html}
"""
import os, re, math, textwrap, html as _html
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.path import Path as MplPath

HERE = os.path.dirname(os.path.abspath(__file__))
TXT = os.path.join(HERE, 'architecture_figure_1_text.md')
OUT = os.environ.get('PAPER_ADM_OUT', HERE)          # write outputs next to the script (repo paper_adm)

W, H = 1400, 900                                      # diagram canvas (svg viewBox / mpl data coords)
VB_Y0 = -80                                           # headroom for the pretrain caption above the outer box
VB_H = H - VB_Y0

# ---- palette (light / paper): role -> (fill, stroke) --------------------------------------------
INK, MUTED, LINE, BG, ACCENT = '#1a1f28', '#5b6572', '#c5ccd6', '#ffffff', '#3355d1'
SVI_S = '#8b95a4'
ROLES = {
    'input': ('#eef2f8', '#6b7a93'), 'token': ('#fdf1dc', '#cf9a3e'),
    'enc':   ('#e6f0fb', '#4d8fd6'), 'attn':  ('#efe7fb', '#8a63d0'),
    'head':  ('#dff3ee', '#2f9c82'), 'stein': ('#fce3ec', '#cf5c86'),
    'dec':   ('#e4f3e1', '#57a24f'), 'mani':  ('#e8eafb', '#5a63cf'),
    'panel': ('#ffffff', '#8b95a4'),          # item-metadata: white fill, clearly-visible border
}

# ---- geometry: box id -> (x, y, w, h, role). RIGHT column is static; the LEFT column (interim_data,
# item_metadata, posterior, pps) is sized to its text and positioned by build() below ----------------
LAYOUT = {
    'manifest':      (430, 25,  940, 200, 'mani'),
    'tokeniser':     (455, 296, 265, 104, 'token'),
    'deepset':       (745, 296, 270, 104, 'enc'),
    'xattn':         (1035,296, 275, 104, 'attn'),
    'head':          (545, 582, 400, 148, 'head'),
    'stein':         (545, 745, 400, 128, 'stein'),
}
# dashed wrappers: (x, y, w, h, stroke)
WRAPPERS = [
    (415, -10, 978, 472, INK),             # PRETRAIN phase: outer box around manifest + encoder
    (435, 275, 895, 165, ACCENT),          # encoder (any-J, any-N), nested inside the pretrain box
    (520, 560, 455, 328, INK),             # DEPLOY phase: fine-tune (specific J, N, rho)
]
# free-standing labels: (text, x, y, rotation, size, colour, anchor, family)
LABELS = [
    ('any N', 1360, 348, 90, 15, MUTED, 'middle', 'mono'),
    ('any J', 882,  466, 0,  15, MUTED, 'middle', 'mono'),
    (r'Fine-tune and deploy for specific J, N, decision-making $\rho$', 748, 548, 0, 13, INK, 'middle', 'sans'),
    (r'Add to library of decision-making amortisers: desired decision-making $\rho$ functions,',
     435, -46, 0, 13, INK, 'start', 'sans'),           # caption OUTSIDE (above) the outer pretrain box
    ('train for any-J and any-N ahead of outcome data', 435, -26, 0, 13, INK, 'start', 'sans'),
]
# arrows as directional cubic Beziers: (p0, exit-dir, p1, enter-dir, label, label_xy, k).
# exit/enter dir in {u,d,l,r} = the direction the curve leaves p0 / the side it enters p1 from, so every
# connector meets a box edge perpendicularly (smooth hooks); k = control-arm length (None = auto). Solid black.
# RIGHT/internal arrows are static; the LEFT-column arrows are built by build() (geodesic edge starts).
STATIC_ARROWS = [
    ((720, 348),  'r', (745, 348), 'l', '', None, None),   # tokeniser -> deep-set
    ((1015, 348), 'r', (1035, 348),'l', '', None, None),   # deep-set -> cross-attention
    ((1172, 400), 'd', (945, 655), 'r', '', None, 110),    # encoder -> head (down, hook left)
    ((745, 730),  'd', (745, 745), 'u', '', None, None),   # head -> stein
    ((1200, 225), 'd', (1200, 275),'u', '', None, None),   # manifest -> encoder
]
ARROWS = []                            # filled by build()
_DIRV = {'u': (0, -1), 'd': (0, 1), 'l': (-1, 0), 'r': (1, 0)}


def build(texts):
    """Size the LEFT column to its text (narrow width, ONE uniform height), evenly spaced + right-aligned,
    and route its arrows from geodesic starts on each box's right edge."""
    global ARROWS, VB_H
    R = 360                            # common RIGHT edge of all four left boxes (narrower than before)
    XF, WF = 30, R - 30                # boxes 1,3,4 left-aligned at x=30
    XM, WM = 70, R - 70                # item metadata indented (clears the x=48 arrow); right edge still R
    GAP, TOP = 44, 40
    hs = [_box_height(k, (XM if k == 'item_metadata' else XF), (WM if k == 'item_metadata' else WF), texts[k])
          for k in ('interim_data', 'item_metadata', 'posterior', 'pps')]
    H = max(hs)                        # one uniform height for all four boxes
    yid = TOP; yim = yid + H + GAP; ypo = yim + H + GAP; ypp = ypo + H + GAP
    LAYOUT['interim_data']  = (XF, yid, WF, H, 'input')
    LAYOUT['item_metadata'] = (XM, yim, WM, H, 'panel')
    LAYOUT['posterior']     = (XF, ypo, WF, H, 'mani')
    LAYOUT['pps']           = (XF, ypp, WF, H, 'dec')
    tl = 455                           # tokeniser left edge; feeds enter at these three heights
    st = LAYOUT['stein']               # stein box, for SD -> stein and the horizontal stein -> PPS
    py = (max(ypp, st[1]) + min(ypp + H, st[1] + st[3])) / 2                                 # stein<->pps overlap
    left = [
        ((48, yid + H), 'd', (48, ypo), 'u', '', None, None),                                # interim data -> SVI
        ((R, yid + H * 0.5), 'r', (tl, 326), 'l', r'$x^{(t)}$', (R + 8, yid + H * 0.5 - 9), 58),       # middle
        ((R, yim + H * 0.5), 'r', (tl, 352), 'l', r'$x^{\mathrm{meta}}$', (R + 8, yim + H * 0.5 - 9), 55),  # middle
        ((R, ypo + H / 3.0), 'r', (tl, 380), 'l', r'$z^{(s)}$', (R + 8, ypo + H / 3.0 - 9), 55),       # 1/3
        ((R, ypo + H * 2 / 3.0), 'r', (st[0], 800), 'l', r'$SD_j$', (R + 8, ypo + H * 2 / 3.0 - 9), 95),  # 2/3
        ((st[0], py), 'l', (R, py), 'r', '', None, None),                                     # stein -> PPS (horiz)
    ]
    ARROWS = left + STATIC_ARROWS
    VB_H = max(900, ypp + H + 20, st[1] + st[3] + 20) - VB_Y0


def _ctrl(p0, d0, p1, d1, k):
    """Cubic control points that leave p0 in direction d0 and arrive at p1 from side d1."""
    if k is None:
        k = min(48.0, max(14.0, 0.4 * math.hypot(p1[0] - p0[0], p1[1] - p0[1])))
    v0, v1 = _DIRV[d0], _DIRV[d1]
    return (p0[0] + v0[0] * k, p0[1] + v0[1] * k), (p1[0] + v1[0] * k, p1[1] + v1[1] * k)
# manifest ρ-cards (kept here — technical config) and net-registry chips
CARDS = [
    ('\u03c1\u2081  rel_change', ['net=scale-feat-paired \u00b7 warp=none',
                                  '\u03b7\u2080\u2208{0,10,20,30,40}% \u00b7 build=svi',
                                  'bvm_shrink=True (Stein)']),
    ('\u03c1\u2082  spr_diff',   ['net=groupdiff \u00b7 warp=none',
                                  '\u03b7\u2080 grid \u00b7 build=rho_id \u00b7 Cohen-d',
                                  '(between-arm g contrast)']),
    ('\u03c1\u2083  gmfr',       ['net=widetok-spr \u00b7 warp=log2/logit',
                                  '\u03b7\u2080 grid \u00b7 build=level',
                                  '\u2026 one joint decision']),
]
REGISTRY = ['scale-feat-paired', 'scale-feat-between', 'widetok-spr', 'groupdiff']


# ================================================================================================
# text file parsing
# ================================================================================================
def parse_text(path):
    boxes, cur = {}, None
    for raw in open(path, encoding='utf-8'):
        line = raw.rstrip('\n')
        if line.startswith('## '):
            cur = line[3:].strip(); boxes[cur] = {'title': cur, 'subtitle': None, 'lines': [], 'hover': ''}
        elif cur and ':' in line and not line.startswith('<!--') and not line.startswith('-'):
            key, _, val = line.partition(':'); key = key.strip(); val = val.strip()
            if key == 'title':    boxes[cur]['title'] = val
            elif key == 'subtitle': boxes[cur]['subtitle'] = val
            elif key == 'line':   boxes[cur]['lines'].append(val)
            elif key == 'hover':  boxes[cur]['hover'] = val
    return boxes


# ================================================================================================
# TeX inline maths -> unicode + <tspan> (for the SVG); the PDF uses matplotlib mathtext directly
# ================================================================================================
_MACRO = [(r'\bar{p}', 'p\u0304'), (r'\bar p', 'p\u0304'), (r'\dotsc', '\u2026'), (r'\ldots', '\u2026'),
          (r'\dots', '\u2026'),
          (r'\rho', '\u03c1'), (r'\theta', '\u03b8'), (r'\eta', '\u03b7'), (r'\sigma', '\u03c3'),
          (r'\tau', '\u03c4'), (r'\times', '\u00d7'), (r'\cdot', '\u00b7'), (r'\to', '\u2192'),
          (r'\in', '\u2208'), (r'\ge', '\u2265'), (r'\le', '\u2264'), (r'\mid', '|'),
          (r'\{', '{'), (r'\}', '}'), (r'\,', '\u2009'), (r'\!', ''), (r'\;', ' ')]


def _match(s, i):                     # index of the '}' matching the '{' at position i
    depth = 0
    for j in range(i, len(s)):
        depth += (s[j] == '{') - (s[j] == '}')
        if depth == 0:
            return j
    return len(s) - 1


def _mpl_math(s):
    """Normalise TeX so matplotlib mathtext accepts it (\\dotsc, \\text{} are not mathtext)."""
    return s.replace(r'\dotsc', r'\ldots').replace(r'\text{', r'\mathrm{')


def _emit_svg(s):                     # char-scan with recursive sub/superscripts (nested braces ok)
    out, i, n = [], 0, len(s)
    esc = lambda t: t.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
    while i < n:
        c = s[i]
        if c in '^_':
            shift = 'super' if c == '^' else 'sub'; i += 1
            if i < n and s[i] == '{':
                j = _match(s, i); seg = s[i + 1:j]; i = j + 1
            else:
                seg = s[i] if i < n else ''; i += 1
            out.append(f'<tspan baseline-shift="{shift}" font-size="72%">{_emit_svg(seg)}</tspan>')
        else:
            j = i
            while j < n and s[j] not in '^_':
                j += 1
            out.append(esc(s[i:j])); i = j
    return ''.join(out)


def _tex_to_svg(inner):
    s = re.sub(r'\\(?:mathrm|text)\{([^{}]*)\}', r'\1', inner)
    s = re.sub(r'\\sqrt\{([^{}]*)\}', lambda m: '\u221a' + m.group(1), s)
    for k, v in _MACRO:
        s = s.replace(k, v)
    return _emit_svg(s)


def tex_line_to_svg(line):
    """Convert a mixed text/$math$ line to SVG-safe markup (unicode + tspans)."""
    parts = re.split(r'(\$[^$]*\$)', line)
    esc = lambda t: t.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
    return ''.join(_tex_to_svg(p[1:-1]) if p.startswith('$') else esc(p) for p in parts)


# ================================================================================================
# shared text layout: yield (kind, y, size, weight) rows for a box
# ================================================================================================
TITLE, SUB, BODY, BODY_SM = 15.0, 11.5, 11.5, 9.6
SMALL = {'tokeniser', 'deepset', 'xattn'}          # narrow encoder-row boxes -> smaller body font
FLOW = {'interim_data', 'item_metadata', 'posterior', 'pps'}   # prose boxes: reflow lines as one paragraph


def _vis(s):
    """Rendered length: TeX macros/braces/$ collapse so maths counts near its glyph width."""
    t = re.sub(r'\\[a-zA-Z]+', 'x', s)
    for ch in '${}^_\\':
        t = t.replace(ch, '')
    return len(t)


def _wrap(line, maxv):
    """Word-wrap `line` to `maxv` rendered chars, keeping each $...$ span atomic."""
    spans = []
    tmp = re.sub(r'\$[^$]*\$', lambda m: (spans.append(m.group(0)), f'\x00{len(spans)-1}\x00')[1], line)
    out, cur, cl = [], [], 0
    for w in tmp.split():
        lw = sum(_vis(spans[int(i)]) for i in re.findall(r'\x00(\d+)\x00', w)) + len(re.sub(r'\x00\d+\x00', '', w))
        if cur and cl + 1 + lw > maxv:
            out.append(cur); cur, cl = [w], lw
        else:
            cl += (1 if cur else 0) + lw; cur.append(w)
    if cur:
        out.append(cur)
    rest = lambda ws: re.sub(r'\x00(\d+)\x00', lambda m: spans[int(m.group(1))], ' '.join(ws))
    return [rest(c) for c in out] or ['']


def _rows(bid, box, txt):
    x, y, w, h, role = box
    body = BODY_SM if bid in SMALL else BODY
    tmax = max(12, int((w - 20) / (0.56 * TITLE)))
    tlines = [txt['title']] if '$' in txt['title'] else (textwrap.wrap(txt['title'], tmax) or [txt['title']])
    bmax = max(16, int((w - 26) / (0.52 * body)))
    rows, yy = [], y + 22
    for tl in tlines:
        rows.append(('title', yy, TITLE, 'bold', tl)); yy += 20
    yy += 2
    if txt.get('subtitle'):
        for sl in _wrap(txt['subtitle'], max(16, int((w - 24) / (0.5 * SUB)))):
            rows.append(('sub', yy, SUB, 'italic', sl)); yy += 18
    yy += 6
    body_lines = [' '.join(txt['lines'])] if (bid in FLOW and txt['lines']) else txt['lines']
    for ln in body_lines:
        for wl in _wrap(ln, bmax):
            rows.append(('body', yy, body, 'normal', wl)); yy += 17
    return rows


def _box_height(bid, x, w, txt):
    r = _rows(bid, (x, 0, w, 0, None), txt)
    return (r[-1][1] + 16) if r else 44


# ================================================================================================
# SVG / HTML renderer
# ================================================================================================
def render_svg(boxes):
    P = []
    P.append(f'<svg viewBox="0 {VB_Y0} {W} {VB_H}" role="img" aria-label="Federated PPS amortiser architecture">')
    P.append('<defs>')
    for mid, col in (('ah', INK), ('ahr', SVI_S), ('ahs', ROLES['mani'][1])):
        P.append(f'<marker id="{mid}" markerWidth="9" markerHeight="9" refX="7" refY="3" '
                 f'orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0 L7,3 L0,6 Z" fill="{col}"/></marker>')
    P.append('</defs>')

    # wrappers (behind)
    for x, y, w, h, col in WRAPPERS:
        P.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="16" fill="none" '
                 f'stroke="{col}" stroke-width="1.6" stroke-dasharray="7 5" opacity="0.85"/>')
    # any-N / any-J double arrows
    P.append(f'<line x1="1345" y1="300" x2="1345" y2="396" stroke="{MUTED}" stroke-width="1.5" '
             f'marker-start="url(#ah)" marker-end="url(#ah)"/>')
    P.append(f'<line x1="452" y1="450" x2="1312" y2="450" stroke="{MUTED}" stroke-width="1.5" '
             f'marker-start="url(#ah)" marker-end="url(#ah)"/>')

    # arrows — solid black cubic connectors
    for p0, d0, p1, d1, lab, lxy, k in ARROWS:
        c0, c1 = _ctrl(p0, d0, p1, d1, k)
        P.append(f'<path d="M{p0[0]},{p0[1]} C{c0[0]:.0f},{c0[1]:.0f} {c1[0]:.0f},{c1[1]:.0f} '
                 f'{p1[0]},{p1[1]}" fill="none" stroke="{INK}" stroke-width="1.7" marker-end="url(#ah)"/>')
        if lab and lxy:
            P.append(f'<text x="{lxy[0]}" y="{lxy[1]}" font-size="11" fill="{MUTED}" '
                     f'font-family="IBM Plex Mono, monospace">{tex_line_to_svg(lab)}</text>')

    # labels
    for text, x, y, rot, size, col, anch, fam in LABELS:
        tr = f' transform="rotate({rot} {x} {y})"' if rot else ''
        ff = 'IBM Plex Mono, monospace' if fam == 'mono' else 'IBM Plex Sans, sans-serif'
        P.append(f'<text x="{x}" y="{y}" font-size="{size}" fill="{col}" text-anchor="{anch}" '
                 f'font-family="{ff}"{tr}>{tex_line_to_svg(text)}</text>')

    # boxes
    for bid, box in LAYOUT.items():
        x, y, w, h, role = box; fill, stroke = ROLES[role]
        txt = boxes[bid]; cx = x + w / 2
        hv = _html.escape(txt.get('hover', ''), quote=True)
        P.append(f'<g class="box" data-hover="{hv}"><title>{hv}</title>')
        P.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="11" fill="{fill}" '
                 f'stroke="{stroke}" stroke-width="1.7"/>')
        if bid == 'manifest':
            P.append(_svg_manifest(box, txt))
        else:
            for kind, yy, size, weight, s in _rows(bid, box, txt):
                cls = 'title' if kind == 'title' else ('sub' if kind == 'sub' else 'body')
                fst = ' font-style="italic"' if kind == 'sub' else ''
                fw = ' font-weight="600"' if kind == 'title' else ''
                fam = 'IBM Plex Sans, sans-serif'
                P.append(f'<text x="{cx}" y="{yy}" font-size="{size}" text-anchor="middle" '
                         f'font-family="{fam}" fill="{INK if kind=="title" else MUTED}"{fw}{fst}>'
                         f'{tex_line_to_svg(s)}</text>')
        P.append('</g>')
    P.append('</svg>')
    return '\n'.join(P)


def _svg_manifest(box, txt):
    x, y, w, h, _ = box; s, mst = ROLES['mani']; stroke = ROLES['mani'][1]
    out = [f'<text x="{x+18}" y="{y+30}" font-size="16" font-weight="600" fill="{INK}" '
           f'font-family="IBM Plex Sans">{tex_line_to_svg(txt["title"])}</text>']
    if txt.get('subtitle'):
        out.append(f'<text x="{x+18}" y="{y+52}" font-size="11" fill="{MUTED}" '
                   f'font-family="IBM Plex Sans">{tex_line_to_svg(txt["subtitle"])}</text>')
    cw, cx0 = 292, [x + 15, x + 322, x + 629]
    for cxx, (ctitle, clines) in zip(cx0, CARDS):
        out.append(f'<rect x="{cxx}" y="{y+64}" width="{cw}" height="86" rx="9" fill="#ffffff" stroke="{stroke}"/>')
        out.append(f'<text x="{cxx+14}" y="{y+86}" font-size="12" font-weight="600" fill="{INK}" '
                   f'font-family="IBM Plex Mono">{_html.escape(ctitle)}</text>')
        for i, cl in enumerate(clines):
            out.append(f'<text x="{cxx+14}" y="{y+104+i*17}" font-size="10.5" fill="{MUTED}" '
                       f'font-family="IBM Plex Mono">{_html.escape(cl)}</text>')
    rx = x + 105
    out.append(f'<text x="{x+18}" y="{y+172}" font-size="11" fill="{MUTED}" font-family="IBM Plex Sans">net registry:</text>')
    for chip in REGISTRY:
        cwd = 12 + len(chip) * 6.7
        out.append(f'<rect x="{rx}" y="{y+160}" width="{cwd:.0f}" height="20" rx="10" fill="#ffffff" stroke="{LINE}"/>')
        out.append(f'<text x="{rx+7}" y="{y+174}" font-size="10" fill="{INK}" font-family="IBM Plex Mono">{_html.escape(chip)}</text>')
        rx += cwd + 12
    return '\n'.join(out)


def write_html(boxes, path):
    svg = render_svg(boxes)
    doc = f"""<title>Federated PPS Amortiser</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;600;700&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
  body{{background:{BG};color:{INK};font-family:"IBM Plex Sans",system-ui,sans-serif;line-height:1.5}}
  .wrap{{max-width:1180px;margin:0 auto;padding-inline:16px;padding-block:32px 56px}}
  .eyebrow{{font-family:"IBM Plex Mono",monospace;font-size:12px;letter-spacing:.14em;text-transform:uppercase;color:{ACCENT};margin:0 0 8px}}
  h1{{font-size:clamp(24px,4vw,34px);margin:0 0 10px;font-weight:700;letter-spacing:-.01em}}
  .lead{{max-width:70ch;color:{MUTED};font-size:15.5px;margin:0 0 18px}}
  .frame{{background:{BG};border:1px solid {LINE};border-radius:14px;overflow-x:auto;padding:10px}}
  svg{{display:block;min-width:1020px;width:100%;height:auto}}
  #cap{{margin-top:12px;min-height:2.6em;font-size:14px;color:{INK};background:#f4f6fa;border:1px solid {LINE};
        border-left:4px solid {ACCENT};border-radius:9px;padding:11px 14px}}
  #cap b{{font-family:"IBM Plex Mono",monospace}}
  .box{{cursor:default}} .box:hover rect{{filter:brightness(0.97)}}
</style>
<div class="wrap">
  <p class="eyebrow">Architecture · figure 1</p>
  <h1>Federated PPS Amortiser</h1>
  <p class="lead">Sources at the top feed the network. A <b>pretrained federated any-J any-N</b> encoder
  (tokeniser → deep-set → cross-attention) is <b>fine-tuned</b> to specific J, N and ρ at the head, then
  <b>Stein–von Mises</b> calibrated, giving the amortised PPS over future data. Hover any block for detail.</p>
  <div class="frame">{svg}</div>
  <div id="cap">Hover a block to see what it does.</div>
</div>
<script>
  var cap=document.getElementById('cap');
  document.querySelectorAll('.box').forEach(function(g){{
    g.addEventListener('mouseenter',function(){{var t=g.getAttribute('data-hover');if(t)cap.textContent=t;}});
  }});
</script>
"""
    open(path, 'w', encoding='utf-8').write(doc)


# ================================================================================================
# matplotlib renderer (white background, paper PDF)
# ================================================================================================
def render_pdf(boxes, pdf_path, png_path):
    fig = plt.figure(figsize=(14, VB_H / 100.0), dpi=200)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, VB_Y0); ax.axis('off')
    fig.patch.set_facecolor(BG)
    sc = 0.72                                            # svg-px -> mpl-pt

    def arrow(p0, d0, p1, d1, k):
        c0, c1 = _ctrl(p0, d0, p1, d1, k)
        path = MplPath([p0, c0, c1, p1],
                       [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
        ax.add_patch(FancyArrowPatch(path=path, arrowstyle='-|>', mutation_scale=11, lw=1.5,
                     color=INK, shrinkA=0, shrinkB=1, zorder=5))

    # wrappers
    for x, y, w, h, col in WRAPPERS:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=16",
                     fill=False, edgecolor=col, lw=1.4, linestyle=(0, (7, 5)), alpha=0.85, zorder=1))
    ax.annotate('', (1345, 396), (1345, 300), arrowprops=dict(arrowstyle='<->', color=MUTED, lw=1.4))
    ax.annotate('', (1312, 450), (452, 450), arrowprops=dict(arrowstyle='<->', color=MUTED, lw=1.4))

    # arrows — solid black cubic connectors
    for p0, d0, p1, d1, lab, lxy, k in ARROWS:
        arrow(p0, d0, p1, d1, k)
        if lab and lxy:
            ax.text(lxy[0], lxy[1], _mpl_math(lab), fontsize=11 * sc, color=MUTED, family='monospace', ha='left', va='center')

    # labels
    _ha = {'middle': 'center', 'start': 'left', 'end': 'right'}
    for text, x, y, rot, size, col, anch, fam in LABELS:
        ax.text(x, y, _mpl_math(text), fontsize=size * sc, color=col,
                family='monospace' if fam == 'mono' else 'sans-serif',
                ha=_ha.get(anch, 'center'), va='center', rotation=rot)

    # boxes
    for bid, box in LAYOUT.items():
        x, y, w, h, role = box; fill, stroke = ROLES[role]; cx = x + w / 2
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=11",
                     facecolor=fill, edgecolor=stroke, lw=1.6, zorder=3))
        if bid == 'manifest':
            _mpl_manifest(ax, box, boxes[bid], sc)
            continue
        for kind, yy, size, weight, s in _rows(bid, box, boxes[bid]):
            ax.text(cx, yy, _mpl_math(s), fontsize=size * sc, ha='center', va='center', zorder=4,
                    color=INK if kind == 'title' else MUTED,
                    fontweight='bold' if kind == 'title' else 'normal',
                    fontstyle='italic' if kind == 'sub' else 'normal',
                    family='sans-serif')
    fig.savefig(pdf_path, facecolor=BG); fig.savefig(png_path, facecolor=BG); plt.close(fig)


def _mpl_manifest(ax, box, txt, sc):
    x, y, w, h, _ = box; stroke = ROLES['mani'][1]
    ax.text(x + 18, y + 30, _mpl_math(txt['title']), fontsize=16 * sc, fontweight='bold', color=INK, va='center', family='sans-serif')
    if txt.get('subtitle'):
        ax.text(x + 18, y + 52, _mpl_math(txt['subtitle']), fontsize=11 * sc, color=MUTED, va='center', family='sans-serif')
    cw, cx0 = 292, [x + 15, x + 322, x + 629]
    for cxx, (ctitle, clines) in zip(cx0, CARDS):
        ax.add_patch(FancyBboxPatch((cxx, y + 64), cw, 86, boxstyle="round,pad=0,rounding_size=9",
                     facecolor='#ffffff', edgecolor=stroke, lw=1.2, zorder=4))
        ax.text(cxx + 14, y + 84, ctitle, fontsize=12 * sc, fontweight='bold', color=INK, va='center',
                family='monospace', zorder=5)
        for i, cl in enumerate(clines):
            ax.text(cxx + 14, y + 103 + i * 17, cl, fontsize=10.5 * sc, color=MUTED, va='center',
                    family='monospace', zorder=5)
    rx = x + 105
    ax.text(x + 18, y + 172, 'net registry:', fontsize=11 * sc, color=MUTED, va='center', family='sans-serif')
    for chip in REGISTRY:
        cwd = 12 + len(chip) * 6.7
        ax.add_patch(FancyBboxPatch((rx, y + 160), cwd, 20, boxstyle="round,pad=0,rounding_size=10",
                     facecolor='#ffffff', edgecolor=LINE, lw=1.0, zorder=4))
        ax.text(rx + 7, y + 170, chip, fontsize=10 * sc, color=INK, va='center', family='monospace', zorder=5)
        rx += cwd + 12


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    boxes = parse_text(TXT)
    build(boxes)                       # size + place the left column and its arrows from the text
    write_html(boxes, os.path.join(OUT, 'architecture_figure_1.html'))
    render_pdf(boxes, os.path.join(OUT, 'architecture_figure_1.pdf'), os.path.join(OUT, 'architecture_figure_1.png'))
    print('Saved architecture_figure_1 -> ' + OUT + '/architecture_figure_1.{html,pdf,png}')
