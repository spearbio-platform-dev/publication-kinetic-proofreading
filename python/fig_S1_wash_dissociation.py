"""
fig_S1_wash_dissociation.py
Supplementary Fig. S1: kinetic cost of wash steps in a wash-based sandwich
immunoassay (a) versus separation-free SPEAR (b), for K_D = 200 pM and 1 nM.

Model
  k_on = k_on2 = 5e5 M-1 s-1 for both antibodies, k_off = K_D * k_on.
  Antibody binding sites 400 pM in the 2 uL binding step (100 pM at 8 uL).
  Analyte 1 pM; normalised curves are independent of analyte concentration
  while the antibodies are in excess (checked at run time from 1 fM to 10 pM).
  Binding reversible and non-cooperative, epitopes independent.
  Dashed = analyte still bound to the capture antibody (or B in SPEAR);
  solid = analyte bound to both antibodies (ternary complex).

  a  wash-based: capture at binding equilibrium, then Wash 1 (0-15 min) /
     detection antibody (15-45) / Wash 2 (45-60) / enzyme (60-65, static
     incubation) / Wash 3 (65-80), read out at 80 min. During washes no free
     antibody is present, so dissociation is irreversible and released analyte
     is carried away.
  b  SPEAR: binding equilibrium in 2 uL, 4-fold dilution to 8 uL, read out at
     8 min. No wash, so dissociated species can rebind.

Normalisation
  100% = analyte bound to the capture antibody (B in SPEAR) at t = 0.
  denom = N_CA(0) for the wash assay (N_CAD(0) = 0), N_CA(0) + N_CAD(0) for SPEAR.
  The a-vs-b ratio at readout does not depend on this choice: both panels are
  divided by the same per-K_D denominator.

Usage
  python fig_S1_wash_dissociation.py      -> ../figures/Fig_S1_wash_dissociation.png
  Optional positional arguments (defaults reproduce the published figure):
    1  probe concentration at the 8 uL stage, pM          (100)
    2  'equilibrium' or 'kinetic' (1 h binding) start     (equilibrium)
    3  'percent' or 'molecules' on the y axis             (percent)
    4  SPEAR readout time, min                            (8)
    5  'lower' or 'upper' panel labels                    (lower)
    6  'noparams' or 'params' (parameter line on top)     (noparams)
  The VER table maps argument sets to the internal version tags of the
  development history; the published figure is 'wq02'.
"""


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
import numpy as np
from scipy.integrate import solve_ivp
import os, sys

# argv[1]: probe concentration at the 8 uL stage, in pM (default 100)
# argv[2]: 'equilibrium' (default) or 'kinetic' (1 h binding) for the SPEAR t=0 state
# argv[3]: 'percent' (default) or 'molecules' for the y-axis units
CB_PM = float(sys.argv[1]) if len(sys.argv) > 1 else 100.0
MODE  = sys.argv[2] if len(sys.argv) > 2 else 'equilibrium'
YUNIT = sys.argv[3] if len(sys.argv) > 3 else 'percent'
assert MODE in ('kinetic', 'equilibrium'), MODE
assert YUNIT in ('percent', 'molecules'), YUNIT
# argv[4]: SPEAR readout time in minutes (default 8)
T_READ = float(sys.argv[4]) if len(sys.argv) > 4 else 8.0
# molecules mode appends a per-K_D suffix to the prefix below
VER = {(100.0, 'kinetic', 'percent', 5.0): 'wd06', (400.0, 'kinetic', 'percent', 5.0): 'we02',
       (100.0, 'equilibrium', 'percent', 5.0): 'wf02', (400.0, 'equilibrium', 'percent', 5.0): 'wg02',
       (100.0, 'equilibrium', 'molecules', 5.0): 'wh', (400.0, 'equilibrium', 'molecules', 5.0): 'wi',
       (100.0, 'kinetic', 'percent', 10.0): 'wk02', (400.0, 'kinetic', 'percent', 10.0): 'wk04',
       (100.0, 'equilibrium', 'percent', 10.0): 'wk01', (400.0, 'equilibrium', 'percent', 10.0): 'wk03',
       (100.0, 'equilibrium', 'molecules', 10.0): 'wm',
       (100.0, 'equilibrium', 'percent', 8.0): 'wn01', (100.0, 'kinetic', 'percent', 8.0): 'wn02',
       (400.0, 'equilibrium', 'percent', 8.0): 'wn03', (400.0, 'kinetic', 'percent', 8.0): 'wn04',
       (100.0, 'equilibrium', 'molecules', 8.0): 'wp',
       }[(CB_PM, MODE, YUNIT, T_READ)]
# argv[5]: panel labels -- 'lower' (default, bold 'a'/'b') or 'upper' ('A'/'B')
# argv[6]: 'noparams' (default) or 'params' to show the parameter line above the panels
LABELS = sys.argv[5] if len(sys.argv) > 5 else 'lower'
PARAMS = sys.argv[6] if len(sys.argv) > 6 else 'noparams'
assert LABELS in ('upper', 'lower'), LABELS
assert PARAMS in ('params', 'noparams'), PARAMS
if (LABELS, PARAMS) != ('upper', 'params'):   # lowercase labels exist for wn01 only
    assert VER == 'wn01' and LABELS == 'lower', 'no version ID assigned for this combination'
    VER = {'params': 'wq01', 'noparams': 'wq02'}[PARAMS]

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, '..', 'figures')
os.makedirs(OUTDIR, exist_ok=True)
# Open Sans (SIL Open Font License, https://fonts.google.com/specimen/Open+Sans).
# Put OpenSans-Regular.ttf and OpenSans-Bold.ttf in ../fonts or point FIG_FONT_DIR
# at them; otherwise matplotlib's default sans-serif is used. Only the Regular face
# is registered, as in the published figures, so fontweight='bold' renders regular;
# the bold face is loaded explicitly where it is used.
FONT_DIR = os.environ.get('FIG_FONT_DIR', os.path.join(HERE, '..', 'fonts'))
_reg = os.path.join(FONT_DIR, 'OpenSans-Regular.ttf')
if os.path.exists(_reg):
    fm.fontManager.addfont(_reg)
    FONT = fm.FontProperties(fname=_reg).get_name()
else:
    FONT = 'DejaVu Sans'
plt.rcParams.update({'font.family': FONT, 'font.size': 12, 'axes.linewidth': 1.0,
                     'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
                     'xtick.major.size': 4, 'ytick.major.size': 4, 'figure.dpi': 300})

_bold = os.path.join(FONT_DIR, 'OpenSans-Bold.ttf')
FONT_BOLD = (fm.FontProperties(fname=_bold) if os.path.exists(_bold)
             else fm.FontProperties(family=FONT, weight='bold'))

C_HI, C_LO, C_NAVY, C_GRID = '#1795BB', '#B85042', '#0B3954', '#E8E8E8'
# annotation accent -- reserved for call-outs, never a data series
C_MARK = '#FF6600'

# ---------------------------------------------------------------- model
kon, Nav = 5e5, 6.022e23
V_init, V_dil, V_det = 2e-6, 8e-6, 2e-6
C_eff = D_eff = 4 * CB_PM * 1e-12   # binding sites at the 2 uL stage (4x the 8 uL value)
A_conc = 1e-12                   # analyte; normalised curves are insensitive to this
t_inc = 3600.0                   # binding incubation (kinetic mode), as in Fig. S5
N_A = A_conc * V_init * Nav
N_C = N_D = C_eff * V_init * Nav
STAGES = [(0, 15, 'Wash 1'), (15, 45, 'Detect Ab'), (45, 60, 'Wash 2'),
          (60, 65, 'Enzyme'), (65, 80, 'Wash 3')]
p = np.array([0, 15, 45, 60, 65, 80]) * 60.0

def site_occupancy(Kd):
    """Per-site occupancy at the 2 uL stage (independent epitopes).

    Because the epitopes are independent, the solid trace at t = 0 is exactly
    this number: [B:Ag:B2] / ([B:Ag] + [B:Ag:B2]) = p1*p2 / (p1*(1-p2) + p1*p2) = p2.
    """
    p_eq = C_eff / (C_eff + Kd)
    if MODE == 'equilibrium':
        return p_eq
    return p_eq * (1 - np.exp(-kon * (C_eff + Kd) * t_inc))

def spear(Kd):
    koff, b = Kd * kon, kon / (V_dil * Nav)
    p = site_occupancy(Kd)
    y0 = [N_A * p * (1 - p), N_A * (1 - p) * p, N_A * p * p]   # CA, AD, CAD
    def f(t, y):
        CA, AD, CAD = y
        A = max(N_A - CA - AD - CAD, 0)
        C = max(N_C - CA - CAD, 0); D = max(N_D - AD - CAD, 0)
        return [b*C*A - koff*CA - b*D*CA + koff*CAD,
                b*A*D - koff*AD - b*C*AD + koff*CAD,
                b*D*CA + b*C*AD - 2*koff*CAD]
    s = solve_ivp(f, (0, 84*60), y0, t_eval=np.linspace(0, 84*60, 1200),
                  rtol=1e-10, atol=1e-1)
    return s.t/60, s.y[0] + s.y[2], s.y[2]      # t, retained analyte, ternary

def elisa(Kd):
    koff, b = Kd * kon, kon / (V_det * Nav)
    NC, ND = C_eff * V_init * Nav, D_eff * V_det * Nav
    def wash(t, y): return [koff*y[1] - koff*y[0], -2*koff*y[1]]
    def make_det(ND):
        def det(t, y):
            CA, CAD, A, AD = y
            C = max(NC - CA - CAD, 0); D = max(ND - AD - CAD, 0)
            return [b*C*A - koff*CA - b*D*CA + koff*CAD,
                    b*D*CA + b*C*AD - 2*koff*CAD,
                    koff*CA - b*C*A - b*D*A + koff*AD,
                    b*A*D - koff*AD - b*C*AD + koff*CAD]
        return det
    T, RA, TC = [], [], []
    y2 = [N_A * C_eff / (Kd + C_eff), 0.0]
    for i, (a_, b_) in enumerate(zip(p[:-1], p[1:])):
        te = np.linspace(a_, b_, 300)
        if i in (1, 3):
            # i=1 'Detect Ab': detection antibody added fresh at D_eff.
            # i=3 'Enzyme'   : a static incubation, not a wash -- there is no flow,
            #   so anything released during the step stays in the 2 uL and can
            #   re-react. Wash 2 has already removed every species lacking capture
            #   antibody, so free Ag = Ag:DetAb = 0 at the start and the only
            #   detection antibody left in the system is the amount still bound in
            #   the sandwich, i.e. ND = CAD at the start of the step.
            ND_step = ND if i == 1 else y2[1]
            s = solve_ivp(make_det(ND_step), (a_, b_), [y2[0], y2[1], 0.0, 0.0],
                          t_eval=te, rtol=1e-10, atol=1e-1)
        else:
            s = solve_ivp(wash, (a_, b_), y2, t_eval=te, rtol=1e-8, atol=1e-6)
        T.append(s.t); RA.append(s.y[0] + s.y[1]); TC.append(s.y[1])
        y2 = [s.y[0, -1], s.y[1, -1]]
    return np.concatenate(T)/60, np.concatenate(RA), np.concatenate(TC)

pct = lambda x: np.asarray(x) / N_A * 100          # % of TOTAL analyte
KDS = [(200e-12, C_HI, '200 pM'), (1e-9, C_LO, '1 nM')]

def denom_elisa(Kd):
    """Analyte on the capture antibody at t = 0 (N_CAD(0) = 0)."""
    return N_A * C_eff / (Kd + C_eff)

def denom_spear(Kd):
    """Analyte on Ab1 at t = 0, i.e. N_CA(0) + N_CAD(0) = N_A * p."""
    return N_A * site_occupancy(Kd)

# ---------------------------------------------------------------- checks
for Kd, _, _ in KDS:
    te, rae, tce = elisa(Kd); ts, ras, tcs = spear(Kd)
    dA, dB = denom_elisa(Kd), denom_spear(Kd)
    assert abs(rae[0]/dA - 1) < 1e-9 and abs(ras[0]/dB - 1) < 1e-9   # dashed starts at 100%
    assert tce[0] == 0                                               # no detection Ab at t=0
    assert np.all(np.diff(rae) <= 1e-6)                              # capture-bound only falls
    assert abs(tcs[0]/dB - site_occupancy(Kd)) < 1e-9                # solid starts at p2
    if MODE == 'equilibrium':
        assert abs(tcs[0]/dB - C_eff/(C_eff + Kd)) < 1e-9             # ... = C/(C+KD)

# normalised curves must be insensitive to analyte concentration while the
# antibody stays in excess (1 pM analyte vs 400 pM sites = 400x). Beyond ~100 pM
# analyte the sites deplete and this no longer holds -- checked, not assumed.
_A0 = N_A
def spear5(Kd):
    ts, _, tcs = spear(Kd)
    return np.interp(T_READ, ts, tcs) / denom_spear(Kd) * 100
globals()['N_A'] = _A0 * 1e-3                       # 1 fM: the dilute limit
ref = {Kd: spear5(Kd) for Kd, _, _ in KDS}
for scale in (1.0, 10.0):                            # 1 pM and 10 pM
    globals()['N_A'] = _A0 * scale
    for Kd, _, _ in KDS:
        assert abs(spear5(Kd) - ref[Kd]) < 0.1, (scale, Kd, spear5(Kd), ref[Kd])
globals()['N_A'] = _A0
print(f'checks OK [{MODE}]: dashed starts at 100%, solid starts at p2, '
      'normalised curves independent of analyte conc from 1 fM to 10 pM')

# ---------------------------------------------------------------- figure
# ================= molecules-on-y variant (scaling check) =================
# One figure per K_D. Each is 1x2 (wash | SPEAR) with ylim = (0, N_cap), where
# N_cap is the analyte bound to capture Ab / B at t = 0 for that K_D. Both dashed
# traces must therefore start exactly at the top of the axis.
if YUNIT == 'molecules':
    SUFFIX = {200e-12: '03', 1e-9: '04'} if T_READ == 5.0 else {200e-12: '01', 1e-9: '02'}
    print(f'\n{"panel":>34} {"N_cap expected":>15} {"dashed at t=0":>15}  match')
    for Kd, col, lbl in KDS:
        te, rae, tce = elisa(Kd)
        ts, ras, tcs = spear(Kd)
        nA, nB = denom_elisa(Kd), denom_spear(Kd)
        top = max(nA, nB)
        fig, axs = plt.subplots(1, 2, figsize=(12.4, 5.0), sharey=True)
        fig.subplots_adjust(wspace=0.09)
        for cc, (ax, t_, ra, tc, ncap, ttl, rx) in enumerate([
                (axs[0], te, rae, tce, nA, 'Wash-based immunoassay', 80),
                (axs[1], ts, ras, tcs, nB, 'SPEAR (separation-free)', T_READ)]):
            if cc == 0:
                for x0, x1, sl in STAGES:
                    if sl.startswith('Wash'):
                        ax.axvspan(x0, x1, color='#8FA3AE', alpha=0.11, lw=0, zorder=0)
                    ax.text((x0+x1)/2, 1.015, sl, ha='center', va='bottom', fontsize=10,
                            color=C_NAVY, transform=ax.get_xaxis_transform())
            ax.axhline(ncap, color=col, lw=0.9, ls=(0, (2, 3)), zorder=2)
            ax.plot(t_, ra, color=col, lw=1.5, ls=(0, (5, 3)), zorder=4)
            ax.plot(t_, tc, color=col, lw=3.0, zorder=5)
            ax.axvline(rx, color=C_MARK, lw=2.2, zorder=6)
            # panel A: label inside at top (dashed trace is well below there).
            # panel B: above the axes -- at t=5 the dashed trace is still near
            # the axis top, so an in-axes label would sit on it.
            if cc == 0:
                ax.text(rx - 2, 0.965, 'Read out', ha='right', va='top', fontsize=11,
                        fontweight='bold', color=C_MARK,
                        transform=ax.get_xaxis_transform(), zorder=8)
            else:
                ax.text(rx, 1.015, 'Read out', ha='center', va='bottom', fontsize=11,
                        fontweight='bold', color=C_MARK,
                        transform=ax.get_xaxis_transform(), zorder=8)
            yv = tc[-1] if cc == 0 else np.interp(rx, t_, tc)
            ax.plot(rx, yv, 'o', ms=8, mfc='white', mec=col, mew=2.0, zorder=8)
            ax.annotate(f'{yv:,.0f}', (rx, yv), textcoords='offset points',
                        xytext=(-13 if cc == 0 else 13, 17), fontsize=12.5,
                        fontweight='bold', color=col,
                        ha='right' if cc == 0 else 'left', va='center', zorder=8)
            ax.set_xlim(0, 84); ax.set_ylim(0, top)
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{v:,.0f}'))
            ax.grid(True, color=C_GRID, lw=0.5); ax.set_axisbelow(True)
            ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
            ax.spines['left'].set_color('#444444'); ax.spines['bottom'].set_color('#444444')
            ax.set_xlabel('Time (min)', fontsize=13, color=C_NAVY)
            ax.set_title(ttl, fontsize=13, fontweight='bold', color=C_NAVY, pad=26)
            ax.text(-0.02, 1.17, 'AB'[cc], transform=ax.transAxes, fontsize=24,
                    fontweight='bold', color='black', ha='right', va='top')
            got = ra[0]
            print(f'{ttl[:24]+"  KD="+lbl:>34} {ncap:15,.0f} {got:15,.0f}  '
                  f'{"OK" if abs(got-ncap) < 1e-6*ncap else "MISMATCH"}')
            assert abs(got - ncap) < 1e-6 * ncap, (ttl, lbl, got, ncap)
        axs[0].set_ylabel('Analyte molecules', fontsize=13, color=C_NAVY)
        fig.text(0.5, 1.055,
                 f'$K_\\mathrm{{D}}$ = $K_\\mathrm{{D2}}$ = {lbl}  ·  '
                 f'analyte 1 pM in 2 $\\mu$L = {N_A:,.0f} molecules  ·  '
                 f'binding sites {C_eff*1e12:.0f} pM at 2 $\\mu$L '
                 f'({CB_PM:.0f} pM at 8 $\\mu$L)\n'
                 + ('binding at equilibrium' if MODE == 'equilibrium' else '1 h binding incubation')
                 + f'  ·  dotted line and axis top = {top:,.0f} molecules bound to '
                   'capture Ab / B at $t$ = 0',
                 ha='center', va='bottom', fontsize=10, color='#6A7A85', linespacing=1.6)
        handles = [Line2D([0], [0], color=col, lw=3.0),
                   Line2D([0], [0], color=col, lw=1.5, ls=(0, (5, 3))),
                   Line2D([0], [0], color=col, lw=0.9, ls=(0, (2, 3)))]
        leg = fig.legend(handles, ['analyte bound by both antibodies',
                                   'analyte bound by capture antibody (or antibody B)',
                                   'captured at $t$ = 0 (axis top)'],
                         loc='upper center', bbox_to_anchor=(0.5, -0.02), ncol=3, fontsize=11,
                         frameon=True, facecolor='white', edgecolor='#C4CDD3', framealpha=1.0,
                         borderpad=0.85, handlelength=3.0, columnspacing=3.0)
        leg.get_frame().set_linewidth(0.9)
        out = os.path.join(OUTDIR, f'Fig_S1_wash_dissociation_{VER}{SUFFIX[Kd]}.png')
        fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        plt.close()
        print(f'  Saved: {out.split("/")[-1]}  ({os.path.getsize(out)/1024:.0f} KB)')
    raise SystemExit

fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.4, 5.0), sharey=True)
fig.subplots_adjust(wspace=0.09)

for ax in (axA, axB):
    ax.set_xlim(0, 84); ax.set_ylim(0, 100)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{v:.0f}%'))
    ax.grid(True, color=C_GRID, lw=0.5); ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#444444'); ax.spines['bottom'].set_color('#444444')
    ax.set_xlabel('Time (min)', fontsize=13, color=C_NAVY)
axA.set_ylabel('% of analyte captured at $t$ = 0', fontsize=13, color=C_NAVY)

# wash phases shaded: no free antibody -> dissociation is irreversible
for x0, x1, lbl in STAGES:
    if lbl.startswith('Wash'):
        axA.axvspan(x0, x1, color='#8FA3AE', alpha=0.11, lw=0, zorder=0)
    axA.text((x0+x1)/2, 1.015, lbl, ha='center', va='bottom', fontsize=10,
             color=C_NAVY, transform=axA.get_xaxis_transform())

for Kd, col, lbl in KDS:
    te, rae, tce = elisa(Kd); dA = denom_elisa(Kd)
    ts, ras, tcs = spear(Kd); dB = denom_spear(Kd)
    axA.plot(te, rae/dA*100, color=col, lw=1.5, ls=(0, (5, 3)), zorder=4)
    axA.plot(te, tce/dA*100, color=col, lw=3.0, zorder=5)
    axB.plot(ts, ras/dB*100, color=col, lw=1.5, ls=(0, (5, 3)), zorder=4)
    axB.plot(ts, tcs/dB*100, color=col, lw=3.0, zorder=5)

def readout(ax, x, vals, side, above=False):
    ax.axvline(x, color=C_MARK, lw=2.2, zorder=6)
    if above:
        # panel B: at t=5 the dashed trace is still near the axis top, so an
        # in-axes label sits on it. Put it above the axes instead, on the same
        # band as panel A's stage labels.
        ax.text(x, 1.015, 'Read out', ha='center', va='bottom', fontsize=11,
                fontweight='bold', color=C_MARK,
                transform=ax.get_xaxis_transform(), zorder=8)
    else:
        ax.text(x + (2 if side > 0 else -2), 97, 'Read out',
                ha='left' if side > 0 else 'right',
                va='top', fontsize=11, fontweight='bold', color=C_MARK, zorder=8)
    for v, col, dx, dy in vals:
        ax.plot(x, v, 'o', ms=8, mfc='white', mec=col, mew=2.0, zorder=8)
        ax.annotate(f'{v:.1f}%' if v >= 1 else f'{v:.2f}%', (x, v),
                    textcoords='offset points', xytext=(side*dx, dy), fontsize=12,
                    fontweight='bold', color=col, ha='left' if side > 0 else 'right',
                    va='center', zorder=8)

vA, vB = [], []
for Kd, col, _ in KDS:
    te, _, tce = elisa(Kd); ts, _, tcs = spear(Kd)
    vA.append((tce[-1]/denom_elisa(Kd)*100, col))
    vB.append((np.interp(T_READ, ts, tcs)/denom_spear(Kd)*100, col))
readout(axA, 80, [(vA[0][0], vA[0][1], 12, 13), (vA[1][0], vA[1][1], 12, 9)], -1)
readout(axB, T_READ, [(vB[0][0], vB[0][1], 11, -30), (vB[1][0], vB[1][1], 11, 0)], +1, above=True)

axA.set_title('Wash-based immunoassay', fontsize=13,
              fontweight='bold', color=C_NAVY, pad=26)
axB.set_title('SPEAR (separation-free)', fontsize=13,
              fontweight='bold', color=C_NAVY, pad=26)
if PARAMS == 'params':
    fig.text(0.5, 1.012,
             f'binding sites {C_eff*1e12:.0f} pM at 2 $\\mu$L '
             f'({CB_PM:.0f} pM at 8 $\\mu$L)  ·  $k_\\mathrm{{on}}$ = 5$\\times$10$^{{5}}$ '
             'M$^{-1}$s$^{-1}$  ·  '
             + ('binding at equilibrium' if MODE == 'equilibrium' else '1 h binding incubation')
             + f'  ·  SPEAR read out at {T_READ:.0f} min',
             ha='center', va='bottom', fontsize=10, color='#6A7A85')
for ax, L in [(axA, 'A'), (axB, 'B')]:
    if LABELS == 'lower':
        ax.text(-0.02, 1.16, L.lower(), transform=ax.transAxes, fontsize=24,
                fontproperties=FONT_BOLD, color='black', ha='right', va='top')
    else:
        ax.text(-0.02, 1.16, L, transform=ax.transAxes, fontsize=24,
                fontweight='bold', color='black', ha='right', va='top')

handles = [Line2D([0], [0], color=C_HI, lw=3.0), Line2D([0], [0], color=C_HI, lw=1.5, ls=(0, (5, 3))),
           Line2D([0], [0], color=C_LO, lw=3.0), Line2D([0], [0], color=C_LO, lw=1.5, ls=(0, (5, 3)))]
labels = ['$K_\\mathrm{D}$ = 200 pM, analyte bound by both antibodies',
          '$K_\\mathrm{D}$ = 200 pM, analyte bound by capture antibody (or antibody B)',
          '$K_\\mathrm{D}$ = 1 nM, analyte bound by both antibodies',
          '$K_\\mathrm{D}$ = 1 nM, analyte bound by capture antibody (or antibody B)']
leg = fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.02), ncol=2,
                 fontsize=11, frameon=True, facecolor='white', edgecolor='#C4CDD3',
                 framealpha=1.0, borderpad=0.85, handlelength=3.0, columnspacing=3.0,
                 labelspacing=0.55)
leg.get_frame().set_linewidth(0.9)

OUT = os.path.join(OUTDIR, 'Fig_S1_wash_dissociation.png' if VER == 'wq02'
                   else f'Fig_S1_wash_dissociation_{VER}.png')
fig.savefig(OUT, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print(f'Saved: {OUT}  ({os.path.getsize(OUT)/1024:.0f} KB)')
print(f'\n{"KD":>8} {"washA@80":>9} {"SPEAR@5":>9} {"ratio":>8} {"SPEAR t=0":>10}')
for Kd, _, lbl in KDS:
    te, _, tce = elisa(Kd); ts, _, tcs = spear(Kd)
    a = tce[-1]/denom_elisa(Kd)*100
    b = np.interp(T_READ, ts, tcs)/denom_spear(Kd)*100
    z = tcs[0]/denom_spear(Kd)*100
    print(f'{lbl:>8} {a:8.2f}% {b:8.1f}% {b/a:7.0f}x {z:9.1f}%')
