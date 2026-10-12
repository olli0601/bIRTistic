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
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, PathPatch
from matplotlib.path import Path as MplPath

HERE = os.path.dirname(os.path.abspath(__file__))
TXT = os.path.join(HERE, 'architecture_figure_1_text.md')
OUT = os.environ.get('PAPER_ADM_OUT', HERE)          # write outputs next to the script (repo paper_adm)

W, H = 1580, 900                                      # diagram canvas (svg viewBox / mpl data coords)
VB_Y0 = -80                                           # headroom for the pretrain caption above the outer box
VB_H = H - VB_Y0

# ---- palette (light / paper): Halifax Harbor (nord), warm tan -> dark purple --------------------
#   https://emilhvitfeldt.github.io/r-color-palettes/discrete/nord/halifax_harbor/
HALIFAX = ['#e1c8af', '#c8af96', '#af967d', '#967d7d', '#644b64', '#4b324b']   # light tan -> dark purple
INK, MUTED, LINE, BG, ACCENT = '#2b2230', '#6b5f6b', '#d8cfd1', '#ffffff', '#644b64'
SVI_S = '#967d7d'


def _tint(hx, f=0.68):                         # fill = palette colour lightened toward white
    r, g, b = int(hx[1:3], 16), int(hx[3:5], 16), int(hx[5:7], 16)
    r, g, b = (round(c + (255 - c) * f) for c in (r, g, b))
    return f'#{r:02x}{g:02x}{b:02x}'


def _role(fill_i, edge_i):                     # pale tint fill + darker palette edge
    return (_tint(HALIFAX[fill_i]), HALIFAX[edge_i])


ROLES = {
    'input': _role(0, 2), 'token': _role(1, 3),
    'enc':   _role(2, 4), 'attn':  _role(3, 5),
    'head':  _role(4, 5), 'stein': _role(3, 5),
    'ft':    _role(2, 4),                       # fine-tune container (overlaps train on the head)
    'train': _role(1, 3),                       # train container (overlaps fine-tune on the head)
    'dec':   _role(4, 5), 'mani':  _role(5, 5),
    'panel': ('#ffffff', HALIFAX[3]),           # item-metadata: white fill, warm border
}

# ---- geometry: box id -> (x, y, w, h, role). RIGHT column is static; the LEFT column (interim_data,
# item_metadata, posterior, pps) is sized to its text and positioned by build() below ----------------
LAYOUT = {
    'scheduler':     (430, 25,  940, 200, 'mani'),
    'itemattn':      (455, 296, 265, 104, 'token'),
    'deepset':       (745, 296, 270, 104, 'enc'),
    'featureattn':   (1035,296, 275, 104, 'attn'),
    'head':          (1010, 600, 250, 165, 'head'),    # in the train∩finetune overlap; set in build()
    'train':         (760, 582, 500, 205, 'train'),    # overlap container (left); set in build()
    'finetune':      (1010, 582, 500, 205, 'ft'),      # overlap container (right); set in build()
    'stein':         (1010, 808, 360, 128, 'stein'),
}
OVERLAP = []           # train/finetune overlap containers (role, x,y,w,h, text-cx), drawn semi-transparent
# dashed wrappers: (x, y, w, h, stroke)
WRAPPERS = [
    (415, -10, 978, 472, INK),             # PRETRAIN phase: outer box around manifest + encoder
    (435, 275, 895, 165, ACCENT),          # encoder (any-J, any-N), nested inside the pretrain box
    (520, 560, 455, 328, INK),             # DEPLOY phase: fine-tune (specific J, N, rho)
]
# free-standing labels: (text, x, y, rotation, size, colour, anchor, family)
LABELS = [
    (r'Fine-tune and deploy for specific application', 748, 548, 0, 13, INK, 'middle', 'sans'),
    (r'Add to library of decision-making amortisers: desired decision-making $\rho$ functions,',
     435, -46, 0, 13, INK, 'start', 'sans'),           # caption OUTSIDE (above) the outer pretrain box
    ('train for any-J and any-N ahead of outcome data', 435, -26, 0, 13, INK, 'start', 'sans'),
]
# arrows as directional cubic Beziers: (p0, exit-dir, p1, enter-dir, label, label_xy, k).
# exit/enter dir in {u,d,l,r} = the direction the curve leaves p0 / the side it enters p1 from, so every
# connector meets a box edge perpendicularly (smooth hooks); k = control-arm length (None = auto). Solid black.
# RIGHT/internal arrows are static; the LEFT-column arrows are built by build() (geodesic edge starts).
STATIC_ARROWS = []     # all deploy-region arrows are built in build() (head/finetune/stein move)
ARROWS = []                            # filled by build()
DASHED = []                            # dashed scheduler-card -> itemattn arrows, filled by build()
LINES = []                             # plain beziers, NO arrowhead (e.g. the posterior fork trunk)
JUNCTIONS = []                         # (x, y) split dots where a trunk line forks into arrows
_DIRV = {'u': (0, -1), 'd': (0, 1), 'l': (-1, 0), 'r': (1, 0)}


def build(texts):
    """Size the LEFT column to its text (narrow width, ONE uniform height), evenly spaced + right-aligned,
    and route its arrows from geodesic starts on each box's right edge."""
    global ARROWS, VB_H
    R = 360                            # common RIGHT edge of all four left boxes (narrower than before)
    XF, WF = 30, R - 30                # boxes 1,3,4 left-aligned at x=30
    XM, WM = 70, R - 70                # manifest box indented; right edge still R
    GAP, TOP = 44, 40
    hs = [_box_height(k, (XM if k == 'manifest' else XF), (WM if k == 'manifest' else WF), texts[k])
          for k in ('interim_data', 'manifest', 'posterior', 'pps')]
    H = max(hs)                        # one uniform height for all four boxes
    # LEFT column: manifest, interim_data, posterior (even), then pps dropped lower-left.
    yim = TOP; yid = yim + H + GAP; ypo = yid + H + GAP
    LAYOUT['manifest']      = (XM, yim, WM, H, 'panel')
    LAYOUT['interim_data']  = (XF, yid, WF, H, 'input')
    LAYOUT['posterior']     = (XF, ypo, WF, H, 'mani')
    # pps (bottom-left) placed after stein so its arrow is horizontal, middle-of-edge to middle-of-edge.
    tl = 455                           # itemattn left edge
    cxL = XF + WF / 2.0                 # middle of the left-column box width
    # encoder row: keep each box's x / width / top (y=296) as-is, grow the (uniform) height to fit the text.
    ETOP = 296
    H_enc = max(104, *[_box_height(k, LAYOUT[k][0], LAYOUT[k][2], texts[k])
                       for k in ('itemattn', 'deepset', 'featureattn')]) + 12  # headroom for stacked ∑ limits
    for k in ('itemattn', 'deepset', 'featureattn'):
        ex, _, ew, _, erole = LAYOUT[k]
        LAYOUT[k] = (ex, ETOP, ew, H_enc, erole)
    ecy, ebot = ETOP + H_enc / 2, ETOP + H_enc
    q1, q2, q3 = ETOP + H_enc / 4, ETOP + H_enc / 2, ETOP + H_enc * 3 / 4   # itemattn entry heights
    sched_left = LAYOUT['scheduler'][0]
    # --- deploy region: a full-width BOTTOM band. Left -> right: stein (von-Mises) | fine-tune | head.
    #     TRAIN rises vertically ABOVE the head, FINE-TUNE extends in from the LEFT; their intersection IS
    #     the head box. Each container places its text in its exclusive (non-overlapping) arm. ---
    enc_bar_bot = ebot + 16                                     # bottom of the encoder top bar
    head_w, head_x, head_y = 275, 1035, 760                     # HEAD/TRAIN aligned with featureattn (x=1035,w=275)
    ft_x = 745                                                  # FINE-TUNE left edge aligned with deepset (x=745)
    st_x, st_w = 455, 265                                       # STEIN (von-Mises) aligned with itemattn (x=455,w=265)
    bh = max(150, _box_height('head', head_x, head_w, texts['head']),          # one shared height for the three
                  _box_height('stein', st_x, st_w, texts['stein']),            # bottom boxes, all on head_y (top)
                  _box_height('finetune', ft_x, head_x - ft_x, texts['finetune']))
    LAYOUT['head']     = (head_x, head_y, head_w, bh, 'head')
    LAYOUT['stein']    = (st_x, head_y, st_w, bh, 'stein')
    LAYOUT['finetune'] = (ft_x, head_y, (head_x + head_w) - ft_x, bh, 'ft')
    head_bot = head_y + bh; hcx = head_x + head_w / 2
    tr_top = head_y - 200                                       # TRAIN: vertical, exclusive text above head
    LAYOUT['train'] = (head_x, tr_top, head_w, head_bot - tr_top, 'train')
    hd, tr, ft, st = LAYOUT['head'], LAYOUT['train'], LAYOUT['finetune'], LAYOUT['stein']
    global OVERLAP   # 6th element = text sub-box (sx,sy,sw,sh): the exclusive (non-overlapping) arm
    OVERLAP = [('train',    tr[0], tr[1], tr[2], tr[3], (tr[0], tr[1], tr[2], head_y - tr[1])),       # top strip
               ('finetune', ft[0], ft[1], ft[2], ft[3], (ft[0], ft[1], head_x - ft[0], ft[3]))]      # left strip
    ftcy = head_y + bh / 2                                      # band mid-line (stein / fine-tune / head / pps centre)
    scy = ftcy
    LAYOUT['pps'] = (XF, head_y, WF, bh, 'dec')                 # pps as tall as the other bottom-band boxes
    tok_cx = tl + 265 / 2
    left = [
        ((cxL, yid + H), 'd', (cxL, ypo), 'u', '', None, None),                              # interim -> posterior
        ((R, yim + H / 3.0), 'r', (sched_left, yim + H / 3.0), 'l', '', None, None),          # manifest -> scheduler
        ((R, yim + H * 2 / 3.0), 'r', (tl, q1), 'l', '', None, 55),                           # manifest -> itemattn (1/4)
        ((R, yid + H * 0.5), 'r', (tl, q2), 'l', '', None, 58),                               # interim -> itemattn (2/4)
        ((R, ypo + H / 3.0), 'r', (tl, q3), 'l', '', None, 55),                               # posterior -> itemattn (3/4)
    ]
    enc = [
        ((720, ecy),  'r', (745, ecy), 'l', '', None, None),                                 # itemattn -> deep-set
        ((1015, ecy), 'r', (1035, ecy),'l', '', None, None),                                 # deep-set -> featureattn
    ]
    stx_r = st[0] + st[2]                                       # von-Mises right edge
    # posterior fork: ONE trunk line leaves posterior at 2/3 height, drops DOWN hugging the dashed box's
    # left edge, curves right below the box, then SPLITS (no dot) into two arrows (-> von-Mises, -> fine-tune).
    global LINES, JUNCTIONS
    jx, jy = 505, head_y - 75                                   # split point: bend earlier (higher), clear of the box edge
    py = ypo + H * 2 / 3
    LINES = [((R, py), 'r', (390, py + 42), 'u', 28),           # trunk pt1: posterior edge -> mid-gap (clear of the box edge)
             ((390, py + 42), 'd', (jx, jy), 'l', 100)]         # trunk pt2: down, bending right early to the split
    JUNCTIONS = []                                              # split read from the branching lines, no dot
    dep = [
        ((1172, ebot), 'd', (tr[0] + tr[2] * 2 / 3, tr[1]), 'u', '', None, None),             # featureattn -> head (2/3)
        ((tok_cx, ebot), 'd', (tr[0] + tr[2] / 3, tr[1]), 'u', '', None, 90),                 # SKIP tokens -> head (1/3)
        ((jx, jy), 'r', (st[0] + st[2] / 2, st[1]), 'u', '', None, 65),                        # fork -> von-Mises (bezier, enters top at 1/2 width)
        ((jx, jy), 'r', (ft[0] + (head_x - ft[0]) * 0.5, ft[1]), 'u', '', None, 90),           # fork -> fine-tune (top)
        ((ft[0], ftcy), 'l', (stx_r, scy), 'r', '', None, None),                              # fine-tune -> von-Mises
        ((st[0], scy), 'l', (R, ftcy), 'r', '', None, None),                                  # von-Mises -> pps
    ]
    ARROWS = left + enc + dep + STATIC_ARROWS
    # dashed arrows: each scheduler ρ-card dispatches its pretrained amortiser into itemattn
    global DASHED
    _dl, _cards_y, _cxs, _cw, _ch = _sched_cards(LAYOUT['scheduler'],
                                                 texts['scheduler'].get('cards') or CARDS,
                                                 texts['scheduler'].get('lines', []))
    card_bot = _cards_y + _ch
    land_x = tl + 265 / 2.0
    DASHED = [((cx + _cw / 2, card_bot), 'd', (land_x, ETOP), 'u', '', None, None) for cx in _cxs]
    sx, sy, sw, _shh, srole = LAYOUT['scheduler']
    LAYOUT['scheduler'] = (sx, sy, sw, (card_bot + 14) - sy, srole)
    # wrappers: PRETRAIN rounded rect (scheduler + encoder) + DEPLOY dashed rect (finetune + head + stein)
    global WRAPPERS
    sch_r = LAYOUT['scheduler'][0] + LAYOUT['scheduler'][2]
    pre_r = sch_r + (LAYOUT['scheduler'][0] - 415)             # symmetric gap to scheduler on left & right
    dep_x = 415; dep_r = pre_r                                 # DEPLOY frame edges aligned with PRETRAIN frame
    dep_top = head_y - 26                                       # just above the band; train pokes out the top
    dep_bot = head_bot + 24
    main_bot = enc_bar_bot + 110                                # PRETRAIN bottom lowered so encoder->head arrows sit inside
    tab_l, tab_bot = 1020, dep_top - 22                         # right column drops to tab_bot; stays ABOVE the deploy frame
    pretrain_poly = ('path', [(415, -10), (pre_r, -10), (pre_r, tab_bot),    # right edge runs full height, no jog
                              (tab_l, tab_bot), (tab_l, main_bot), (415, main_bot)], INK)
    WRAPPERS = [
        pretrain_poly,                                                             # PRETRAIN frame with TRAIN tab
        (dep_x, dep_top, dep_r - dep_x, dep_bot - dep_top, INK),                    # DEPLOY dashed rect
    ]
    LABELS[0] = (r'Fine-tune and deploy for specific application', dep_x + 20, dep_bot + 22, 0, 13, INK, 'start', 'sans')  # below the deploy frame, mirroring the top caption
    VB_H = max(900, head_bot + 20, dep_bot + 44) - VB_Y0


def _round_poly_segs(pts, r):
    """Per-vertex (p_in, p, p_out) for a rounded polygon: cut back r along each edge, round through p."""
    n = len(pts); segs = []
    for i in range(n):
        p, a, b = pts[i], pts[(i - 1) % n], pts[(i + 1) % n]
        va = (a[0] - p[0], a[1] - p[1]); la = math.hypot(*va)
        vb = (b[0] - p[0], b[1] - p[1]); lb = math.hypot(*vb)
        rr = min(r, la / 2, lb / 2)
        segs.append(((p[0] + va[0] / la * rr, p[1] + va[1] / la * rr), p,
                     (p[0] + vb[0] / lb * rr, p[1] + vb[1] / lb * rr)))
    return segs


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
            cur = line[3:].strip()
            boxes[cur] = {'title': cur, 'subtitle': None, 'lines': [], 'hover': '', 'cards': [], 'registry': []}
        elif cur and ':' in line and not line.startswith('<!--') and not line.startswith('-'):
            key, _, val = line.partition(':'); key = key.strip(); val = val.strip()
            if key == 'title':    boxes[cur]['title'] = val
            elif key == 'subtitle': boxes[cur]['subtitle'] = val
            elif key == 'line':   boxes[cur]['lines'].append(val)
            elif key == 'hover':  boxes[cur]['hover'] = val
            elif key == 'card':   boxes[cur]['cards'].append((val, []))        # manifest ρ-card title
            elif key == 'cardline' and boxes[cur]['cards']:
                boxes[cur]['cards'][-1][1].append(val)                         # a line inside the last card
            elif key == 'registry': boxes[cur]['registry'] = [c.strip() for c in val.split(',') if c.strip()]
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
TITLE, SUB, BODY, BODY_SM = 15.0, 11.5, 11.5, 11.5   # uniform body font across all boxes
SMALL = set()                                      # (was the narrow encoder boxes; now uniform body size)
FLOW = {'interim_data', 'manifest', 'posterior', 'pps'}        # prose boxes: reflow lines as one paragraph


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
            tall = ('\\sum' in wl) or ('∑' in wl)   # summation limits render ~1.6x row height
            if tall:
                yy += 7                                  # extra space ABOVE a ∑ line
            rows.append(('body', yy, body, 'normal', wl))
            yy += 28 if tall else 17
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
    P.append(f'<marker id="ahm" markerWidth="14" markerHeight="14" refX="10.5" refY="4.5" '      # 1.5x (dashed)
             f'orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0 L10.5,4.5 L0,9 Z" fill="{MUTED}"/></marker>')
    P.append('</defs>')

    # Layers are drawn back-to-front so box FILLS sit UNDER the dashed wrappers and ALL text, which read true.
    # layer 1: overlap container fills (semi-transparent) ---------------------------------------------
    for oid, ox, oy, ow, oh, sub in OVERLAP:
        fill, stroke = ROLES['train' if oid == 'train' else 'ft']
        P.append(f'<rect x="{ox}" y="{oy}" width="{ow}" height="{oh}" rx="14" fill="{fill}" '
                 f'fill-opacity="0.5" stroke="{stroke}" stroke-width="1.7"/>')
    # layer 2: box fills (+ hover group) and the scheduler manifest -----------------------------------
    for bid, box in LAYOUT.items():
        if bid in ('train', 'finetune'):                       # drawn as overlap containers above
            continue
        x, y, w, h, role = box; fill, stroke = ROLES[role]; txt = boxes[bid]
        hv = _html.escape(txt.get('hover', ''), quote=True)
        P.append(f'<g class="box" data-hover="{hv}"><title>{hv}</title>')
        P.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="11" fill="{fill}" '
                 f'stroke="{stroke}" stroke-width="1.7"/>')
        if bid == 'scheduler':
            P.append(_svg_manifest(box, txt))
        P.append('</g>')
    # layer 3: dashed wrappers, OVER the fills so frames read true black ------------------------------
    for wr in WRAPPERS:
        if len(wr) == 3 and wr[0] == 'path':
            pts, col = wr[1], wr[2]
            segs = _round_poly_segs(pts, 16)
            d = f'M{segs[-1][2][0]:.1f},{segs[-1][2][1]:.1f}'
            for pin, p, pout in segs:
                d += f' L{pin[0]:.1f},{pin[1]:.1f} Q{p[0]:.1f},{p[1]:.1f} {pout[0]:.1f},{pout[1]:.1f}'
            P.append(f'<path d="{d} Z" fill="none" stroke="{col}" stroke-width="1.6" '
                     f'stroke-dasharray="7 5" opacity="0.85"/>')
        else:
            x, y, w, h, col = wr
            P.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="16" fill="none" '
                     f'stroke="{col}" stroke-width="1.6" stroke-dasharray="7 5" opacity="0.85"/>')
    # layer 4: connectors (plain lines, arrows, dashed dispatch, split dots) --------------------------
    arrow_labels = []
    for p0, d0, p1, d1, k in LINES:
        c0, c1 = _ctrl(p0, d0, p1, d1, k)
        P.append(f'<path d="M{p0[0]},{p0[1]} C{c0[0]:.0f},{c0[1]:.0f} {c1[0]:.0f},{c1[1]:.0f} '
                 f'{p1[0]},{p1[1]}" fill="none" stroke="{INK}" stroke-width="1.7"/>')
    for p0, d0, p1, d1, lab, lxy, k in ARROWS:
        c0, c1 = _ctrl(p0, d0, p1, d1, k)
        P.append(f'<path d="M{p0[0]},{p0[1]} C{c0[0]:.0f},{c0[1]:.0f} {c1[0]:.0f},{c1[1]:.0f} '
                 f'{p1[0]},{p1[1]}" fill="none" stroke="{INK}" stroke-width="1.7" marker-end="url(#ah)"/>')
        if lab and lxy:
            arrow_labels.append((lxy, lab))
    for p0, d0, p1, d1, lab, lxy, k in DASHED:
        c0, c1 = _ctrl(p0, d0, p1, d1, k)
        P.append(f'<path d="M{p0[0]},{p0[1]} C{c0[0]:.0f},{c0[1]:.0f} {c1[0]:.0f},{c1[1]:.0f} '
                 f'{p1[0]},{p1[1]}" fill="none" stroke="{MUTED}" stroke-width="2.25" '
                 f'stroke-dasharray="9 6" marker-end="url(#ahm)"/>')
    for jxp, jyp in JUNCTIONS:
        P.append(f'<circle cx="{jxp}" cy="{jyp}" r="4.5" fill="{INK}"/>')
    # layer 5: ALL text on top (true black); pointer-events none so hover reaches the box rect ---------
    for bid, box in LAYOUT.items():
        if bid in ('train', 'finetune', 'scheduler'):          # train/ft below; scheduler text is the manifest
            continue
        x, y, w, h, role = box; cx = x + w / 2
        for kind, yy, size, weight, s in _rows(bid, box, boxes[bid]):
            fst = ' font-style="italic"' if kind == 'sub' else ''
            fw = ' font-weight="600"' if kind == 'title' else ''
            P.append(f'<text x="{cx}" y="{yy}" font-size="{size}" text-anchor="middle" pointer-events="none" '
                     f'font-family="IBM Plex Sans, sans-serif" fill="{INK}"{fw}{fst}>'
                     f'{tex_line_to_svg(s)}</text>')
    for oid, ox, oy, ow, oh, sub in OVERLAP:
        sx, sy, sw, sh = sub; otcx = sx + sw / 2
        for kind, yy, size, weight, s in _rows(oid, (sx, sy, sw, sh, None), boxes[oid]):
            fw = ' font-weight="600"' if kind == 'title' else ''
            P.append(f'<text x="{otcx:.0f}" y="{yy}" font-size="{size}" text-anchor="middle" pointer-events="none" '
                     f'font-family="IBM Plex Sans, sans-serif" fill="{INK}"{fw}>'
                     f'{tex_line_to_svg(s)}</text>')
    for lxy, lab in arrow_labels:
        P.append(f'<text x="{lxy[0]}" y="{lxy[1]}" font-size="11" fill="{MUTED}" pointer-events="none" '
                 f'font-family="IBM Plex Mono, monospace">{tex_line_to_svg(lab)}</text>')
    for text, x, y, rot, size, col, anch, fam in LABELS:
        tr = f' transform="rotate({rot} {x} {y})"' if rot else ''
        ff = 'IBM Plex Mono, monospace' if fam == 'mono' else 'IBM Plex Sans, sans-serif'
        P.append(f'<text x="{x}" y="{y}" font-size="{size}" fill="{col}" text-anchor="{anch}" pointer-events="none" '
                 f'font-family="{ff}"{tr}>{tex_line_to_svg(text)}</text>')
    P.append('</svg>')
    return '\n'.join(P)


def _sched_cards(box, cards, lines):
    """Scheduler box geometry: wrap the description `line:` text, then lay the ρ-cards in a row below it.
    Returns (wrapped desc lines, cards_y, [card x's], card width, card height)."""
    x, y, w, h, _ = box
    maxv = max(40, int((w - 40) / (0.5 * 11)))
    dl = _wrap(' '.join(lines), maxv) if lines else []
    cards_y = y + 46 + len(dl) * 16 + 12
    n = max(1, len(cards)); gap = 14
    cw = (w - 30 - (n - 1) * gap) / n
    cxs = [x + 15 + i * (cw + gap) for i in range(len(cards))]
    nl = max((len([c for c in cl if c.strip()]) for _, cl in cards), default=0)
    return dl, cards_y, cxs, cw, 24 + max(1, nl) * 15 + 8


def _svg_manifest(box, txt):
    x, y, w, h, _ = box; stroke = ROLES['mani'][1]; cards = txt.get('cards') or CARDS
    out = [f'<text x="{x+18}" y="{y+28}" font-size="15" font-weight="600" fill="{INK}" '
           f'font-family="IBM Plex Sans">{tex_line_to_svg(txt["title"])}</text>']
    dl, cards_y, cxs, cw, ch = _sched_cards(box, cards, txt.get('lines', []))
    for i, ln in enumerate(dl):
        out.append(f'<text x="{x+18}" y="{y+46+i*16}" font-size="11" fill="{INK}" '
                   f'font-family="IBM Plex Sans">{tex_line_to_svg(ln)}</text>')
    for cxx, (ctitle, clines) in zip(cxs, cards):
        out.append(f'<rect x="{cxx:.0f}" y="{cards_y:.0f}" width="{cw:.0f}" height="{ch:.0f}" rx="9" '
                   f'fill="#ffffff" stroke="{stroke}"/>')
        out.append(f'<text x="{cxx+12:.0f}" y="{cards_y+20:.0f}" font-size="10.5" font-weight="600" fill="{INK}" '
                   f'font-family="IBM Plex Mono">{tex_line_to_svg(ctitle)}</text>')
        for i, cl in enumerate(c for c in clines if c.strip()):
            out.append(f'<text x="{cxx+12:.0f}" y="{cards_y+40+i*15:.0f}" font-size="9" fill="{MUTED}" '
                       f'font-family="IBM Plex Mono">{_html.escape(cl)}</text>')
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
    fig = plt.figure(figsize=(W / 100.0, VB_H / 100.0), dpi=200)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(VB_Y0 + VB_H, VB_Y0); ax.axis('off')
    fig.patch.set_facecolor(BG)
    sc = 0.72                                            # svg-px -> mpl-pt

    def arrow(p0, d0, p1, d1, k):
        c0, c1 = _ctrl(p0, d0, p1, d1, k)
        path = MplPath([p0, c0, c1, p1],
                       [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
        ax.add_patch(FancyArrowPatch(path=path, arrowstyle='-|>', mutation_scale=11, lw=1.5,
                     color=INK, shrinkA=0, shrinkB=1, zorder=4))

    # wrappers — rect (x,y,w,h,col) or L-shape ('path',[pts],col)
    for wr in WRAPPERS:
        if len(wr) == 3 and wr[0] == 'path':
            pts, col = wr[1], wr[2]
            segs = _round_poly_segs(pts, 16)
            verts = [segs[-1][2]]; codes = [MplPath.MOVETO]
            for pin, p, pout in segs:
                verts += [pin, p, pout]; codes += [MplPath.LINETO, MplPath.CURVE3, MplPath.CURVE3]
            verts += [verts[0]]; codes += [MplPath.CLOSEPOLY]
            ax.add_patch(PathPatch(MplPath(verts, codes), fill=False, edgecolor=col, lw=1.4,
                         linestyle=(0, (7, 5)), alpha=0.85, zorder=3))
        else:
            x, y, w, h, col = wr
            ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=16",
                         fill=False, edgecolor=col, lw=1.4, linestyle=(0, (7, 5)), alpha=0.85, zorder=3))
    # overlap containers (train + fine-tune) — semi-transparent fills
    for oid, ox, oy, ow, oh, sub in OVERLAP:
        fill, stroke = ROLES['train' if oid == 'train' else 'ft']
        ax.add_patch(FancyBboxPatch((ox, oy), ow, oh, boxstyle="round,pad=0,rounding_size=14",
                     facecolor=fill, edgecolor=stroke, lw=1.6, alpha=0.5, zorder=1))
        sx, sy, sw, sh = sub; otcx = sx + sw / 2
        for kind, yy, size, weight, s in _rows(oid, (sx, sy, sw, sh, None), boxes[oid]):
            ax.text(otcx, yy, _mpl_math(s), fontsize=size * sc, ha='center', va='center', zorder=5,
                    color=INK,
                    fontweight='bold' if kind == 'title' else 'normal', family='sans-serif')
    # plain connector lines (no arrowhead) — e.g. the posterior fork trunk
    for p0, d0, p1, d1, k in LINES:
        c0, c1 = _ctrl(p0, d0, p1, d1, k)
        path = MplPath([p0, c0, c1, p1], [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
        ax.add_patch(PathPatch(path, fill=False, edgecolor=INK, lw=1.5, zorder=4))
    # arrows — solid black cubic connectors
    for p0, d0, p1, d1, lab, lxy, k in ARROWS:
        arrow(p0, d0, p1, d1, k)
        if lab and lxy:
            ax.text(lxy[0], lxy[1], _mpl_math(lab), fontsize=11 * sc, color=MUTED, family='monospace', ha='left', va='center')
    # fork split dots
    for jxp, jyp in JUNCTIONS:
        ax.add_patch(plt.Circle((jxp, jyp), 4.5, color=INK, zorder=4))
    # dashed arrows — scheduler cards dispatch pretrained amortisers into itemattn
    for p0, d0, p1, d1, lab, lxy, k in DASHED:
        c0, c1 = _ctrl(p0, d0, p1, d1, k)
        path = MplPath([p0, c0, c1, p1], [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
        ax.add_patch(FancyArrowPatch(path=path, arrowstyle='-|>', mutation_scale=15, lw=2.0,
                     color=MUTED, linestyle=(0, (7, 4)), shrinkA=0, shrinkB=1, zorder=4))

    # labels
    _ha = {'middle': 'center', 'start': 'left', 'end': 'right'}
    for text, x, y, rot, size, col, anch, fam in LABELS:
        ax.text(x, y, _mpl_math(text), fontsize=size * sc, color=col,
                family='monospace' if fam == 'mono' else 'sans-serif',
                ha=_ha.get(anch, 'center'), va='center', rotation=rot, zorder=5)

    # boxes
    for bid, box in LAYOUT.items():
        if bid in ('train', 'finetune'):                       # drawn as overlap containers above
            continue
        x, y, w, h, role = box; fill, stroke = ROLES[role]; cx = x + w / 2
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=11",
                     facecolor=fill, edgecolor=stroke, lw=1.6, zorder=2))
        if bid == 'scheduler':
            _mpl_manifest(ax, box, boxes[bid], sc)
            continue
        for kind, yy, size, weight, s in _rows(bid, box, boxes[bid]):
            ax.text(cx, yy, _mpl_math(s), fontsize=size * sc, ha='center', va='center', zorder=5,
                    color=INK,
                    fontweight='bold' if kind == 'title' else 'normal',
                    fontstyle='italic' if kind == 'sub' else 'normal',
                    family='sans-serif')
    fig.savefig(pdf_path, facecolor=BG); fig.savefig(png_path, facecolor=BG); plt.close(fig)


def _mpl_manifest(ax, box, txt, sc):
    x, y, w, h, _ = box; stroke = ROLES['mani'][1]; cards = txt.get('cards') or CARDS
    ax.text(x + 18, y + 28, _mpl_math(txt['title']), fontsize=15 * sc, fontweight='bold', color=INK, va='center', family='sans-serif', zorder=5)
    dl, cards_y, cxs, cw, ch = _sched_cards(box, cards, txt.get('lines', []))
    for i, ln in enumerate(dl):
        ax.text(x + 18, y + 46 + i * 16, _mpl_math(ln), fontsize=11 * sc, color=INK, va='center', family='sans-serif', zorder=5)
    for cxx, (ctitle, clines) in zip(cxs, cards):
        ax.add_patch(FancyBboxPatch((cxx, cards_y), cw, ch, boxstyle="round,pad=0,rounding_size=9",
                     facecolor='#ffffff', edgecolor=stroke, lw=1.2, zorder=2.2))
        ax.text(cxx + 12, cards_y + 16, _mpl_math(ctitle), fontsize=10.5 * sc, fontweight='bold', color=INK, va='center',
                family='monospace', zorder=5)
        for i, cl in enumerate(c for c in clines if c.strip()):
            ax.text(cxx + 12, cards_y + 36 + i * 15, cl, fontsize=9 * sc, color=MUTED, va='center',
                    family='monospace', zorder=5)


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    boxes = parse_text(TXT)
    build(boxes)                       # size + place the left column and its arrows from the text
    write_html(boxes, os.path.join(OUT, 'architecture_figure_1.html'))
    render_pdf(boxes, os.path.join(OUT, 'architecture_figure_1.pdf'), os.path.join(OUT, 'architecture_figure_1.png'))
    print('Saved architecture_figure_1 -> ' + OUT + '/architecture_figure_1.{html,pdf,png}')
