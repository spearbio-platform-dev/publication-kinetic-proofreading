"""
fig_1d_amplicon_kinetics.py
Main-text Fig. 1d: total amplicon (signal + background) versus extension time
for a 1-step proximity amplicon-generation (PAG) assay and 2-step SPEAR, at
analyte concentrations of 0 to 1000 fM.

Model (SI equation numbers)
  1-step (PAG) signal          Eq. 2.21     r1*(1-exp(-lam t)) * occ
  1-step background            Eq. 2.10.6   C_B*(1-exp(-C_B2 k_bi t))
  2-step (SPEAR) signal        Eq. 2.22.1   r1^2*(1-(1+lam t)exp(-lam t)) * occ
  2-step background            Eq. 2.18     C_B*(1-(1+C_B2 k_bi t)exp(-C_B2 k_bi t))
  lam = k_uni + K_D k_on + K_D2 k_on2,  r1 = k_uni/lam
  occ = [C_bind/(K_D + C_bind)]^2, evaluated at the binding-step concentration
        (4 x C_B, i.e. 400 pM in 2 uL); extension and background use C_B at the
        8 uL stage. Analyte and amplicon concentrations are at the 8 uL stage.

Parameters (SI section 2.1.7)
  k_uni = 0.01 s-1, k_bi = 5e3 M-1s-1, k_on = k_on2 = 5e5 M-1s-1,
  K_D = K_D2 = 1 nM, C_B = C_B2 = 100 pM at 8 uL, t = 0-15 min,
  C_A = 0 / 0.01 / 0.1 / 1 / 10 / 100 / 1000 fM.

Style: Arial (falls back to matplotlib's default sans-serif if Arial is not
installed); 8 pt panel letter, 7 pt labels, 6 pt ticks.

Usage
  python fig_1d_amplicon_kinetics.py    -> ../figures/Fig_1d_amplicon_kinetics.{pdf,svg,png}
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import os

VER = 'hk02'
HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, '..', 'figures')
os.makedirs(OUTDIR, exist_ok=True)
_ARIAL = os.environ.get('FIG_ARIAL_DIR', '/System/Library/Fonts/Supplemental')
_have_arial = False
for f in ['Arial.ttf', 'Arial Bold.ttf']:
    if os.path.exists(os.path.join(_ARIAL, f)):
        fm.fontManager.addfont(os.path.join(_ARIAL, f))
        _have_arial = True
plt.rcParams.update({'font.family': 'Arial' if _have_arial else 'sans-serif',
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})

# ---------------------------------------------------------------- model
k_uni, k_bi, k_on = 0.01, 5e3, 5e5
KD = 1e-9
CB = 100e-12                 # 8 uL extension step
CB_BIND = 4 * CB             # 2 uL binding step
t_min = np.arange(0.1, 15.0001, 0.1)     # log axis: start at the first step, as in the original
ts = t_min * 60.0
CA = np.array([0, 0.01, 0.1, 1, 10, 100, 1000]) * 1e-15

lam = k_uni + 2 * KD * k_on
r1 = k_uni / lam
occ = (CB_BIND / (KD + CB_BIND)) ** 2
f1 = lambda x: -np.expm1(-x)
f2 = lambda x: -np.expm1(-x) - x * np.exp(-x)

bg_1  = CB * f1(CB * k_bi * ts)            # Eq. 2.10.6
bg_2  = CB * f2(CB * k_bi * ts)            # Eq. 2.18
sig_1 = r1 * f1(lam * ts) * occ            # Eq. 2.21
sig_2 = r1 ** 2 * f2(lam * ts) * occ       # Eq. 2.22.1
tot_1 = (np.outer(CA, sig_1) + bg_1) * 1e15    # fM
tot_2 = (np.outer(CA, sig_2) + bg_2) * 1e15

assert abs(r1 ** 2 * f2(k_uni * 1e9) - r1 ** 2) < 1e-12
assert np.all(np.diff(tot_1, axis=1) > 0) and np.all(np.diff(tot_2, axis=1) > 0)
i8 = np.argmin(abs(t_min - 8)); i15 = -1
print(f'occ = {occ:.4f}, r1 = {r1:.4f}, lam = {lam:.4f} s-1')
print(f'{"C_A (fM)":>9} | {"1-step 8 min":>12} {"15 min":>9} | {"2-step 8 min":>12} {"15 min":>9}   (total amplicon, fM)')
for i, ca in enumerate(CA * 1e15):
    print(f'{ca:9g} | {tot_1[i, i8]:12.4g} {tot_1[i, i15]:9.4g} | {tot_2[i, i8]:12.4g} {tot_2[i, i15]:9.4g}')

# ---------------------------------------------------------------- figure (all layout in points)
COLORS = ['#D0D4D4', '#2094BC', '#A854A0', '#40B448', '#F8A030', '#043854', '#F07C80']
LABELS = ['0 fM', '0.01 fM', '0.1 fM', '1 fM', '10 fM', '100 fM', '1000 fM']
PT = 1 / 72
AX_W, AX_H, GAP = 0.748, 1.117, 0.277          # inches
LEFT, RIGHT, TOP, BOTTOM = 0.44, 0.12, 0.30, 0.66
FIG_W = LEFT + 2 * AX_W + GAP + RIGHT
FIG_H = TOP + AX_H + BOTTOM
FS_TICK, FS_LAB, FS_PANEL = 6, 7, 8
Y_LO, Y_HI = -7, 3

fig = plt.figure(figsize=(FIG_W, FIG_H))
axes = []
for n, (data, ttl) in enumerate([(tot_1, 'PAG'), (tot_2, 'SPEAR')]):
    ax = fig.add_axes([(LEFT + n * (AX_W + GAP)) / FIG_W, BOTTOM / FIG_H, AX_W / FIG_W, AX_H / FIG_H])
    for i in range(len(CA)):
        ax.plot(t_min, data[i], color=COLORS[i], lw=0.6, zorder=3 + i, solid_capstyle='butt')
    ax.set_yscale('log')
    ax.set_xlim(0, 15); ax.set_ylim(10.0 ** Y_LO, 10.0 ** Y_HI)
    ax.set_xticks([0, 5, 10, 15])
    ax.set_yticks(10.0 ** np.arange(Y_LO, Y_HI + 1)); ax.minorticks_off()
    ax.set_yticklabels([])
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_linewidth(0.5)
    ax.tick_params(axis='both', which='major', direction='out', length=2.7, width=0.5,
                   labelsize=FS_TICK, pad=1.5, color='black')
    # decade labels drawn by hand so the exponent uses a hyphen, as in the original
    for e in range(Y_LO, Y_HI + 1):
        ax.annotate(f'{e}'.replace('-', '-'), xy=(0, 10.0 ** e), xycoords=('axes fraction', 'data'),
                    xytext=(-4.0, 1.9), textcoords='offset points', ha='right', va='center',
                    fontsize=FS_TICK * 0.62)
        # width of the exponent in pt (Arial: digit 0.556 em, hyphen 0.333 em)
        w = 0.62 * FS_TICK * (0.556 * len(str(abs(e))) + (0.333 if e < 0 else 0))
        ax.annotate('10', xy=(0, 10.0 ** e), xycoords=('axes fraction', 'data'),
                    xytext=(-4.0 - w, -0.3), textcoords='offset points', ha='right',
                    va='center', fontsize=FS_TICK)
    ax.annotate('Time (min)', xy=(0.5, 0), xycoords='axes fraction', xytext=(0, -17.1),
                textcoords='offset points', ha='center', va='baseline', fontsize=FS_LAB,
                fontweight='bold')
    ax.annotate(ttl, xy=(0.5, 1), xycoords='axes fraction', xytext=(0, 9.3),
                textcoords='offset points', ha='center', va='baseline', fontsize=FS_LAB,
                fontweight='bold')
    axes.append(ax)

axL = axes[0]
axL.annotate('[Amplicon] (fM)', xy=(0, 0.5), xycoords='axes fraction', xytext=(-22.5, 0),
             textcoords='offset points', ha='center', va='center', rotation=90,
             fontsize=FS_LAB, fontweight='bold')
axL.annotate('d', xy=(0, 1), xycoords='axes fraction', xytext=(-26.5, 10.9),
             textcoords='offset points', ha='left', va='baseline', fontsize=FS_PANEL,
             fontweight='bold')

# legend: 4 columns x 2 rows, positions (pt) measured from the original relative to the left axis
COLS = [(-9.3, 0.7), (32.1, 42.3), (66.4, 76.4), (104.3, 114.1)]     # (line start, text start)
ROWS = [-29.4, -40.4]                                                # centre, below the x axis
GRID = [[6, 4, 2, 0], [5, 3, 1, None]]                               # index into LABELS
for r, row in enumerate(GRID):
    for c, i in enumerate(row):
        if i is None:
            continue
        x0, xt = COLS[c]
        axL.plot([x0 * PT / AX_W, (x0 + 7.0) * PT / AX_W], [ROWS[r] * PT / AX_H] * 2,
                 transform=axL.transAxes, color=COLORS[i], lw=0.75, clip_on=False,
                 solid_capstyle='butt')
        axL.annotate(LABELS[i], xy=(0, 0), xycoords='axes fraction', xytext=(xt, ROWS[r]),
                     textcoords='offset points', ha='left', va='center_baseline',
                     fontsize=FS_LAB, fontweight='bold')

base = os.path.join(OUTDIR, 'Fig_1d_amplicon_kinetics')
fig.savefig(base + '.pdf', facecolor='white')
fig.savefig(base + '.svg', facecolor='white')
fig.savefig(base + '.png', dpi=900, facecolor='white')
plt.close()
for ext in ('pdf', 'svg', 'png'):
    print(f'Saved: {base}.{ext}  ({os.path.getsize(base + "." + ext)/1024:.0f} KB)')
