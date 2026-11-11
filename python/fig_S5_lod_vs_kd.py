"""
fig_S5_lod_vs_kd.py
Supplementary Fig. S5: achievable LOD versus antibody affinity (K_D = K_D2) for
allostery-based binders, proximity amplicon-generation (PAG) assays and SPEAR.

Allostery (equilibrium relations of matlab/Script1_Allostery_250108.m, solved
symbolically here):
    ABs/(A*Bs) = 1/KD     A + ABs = CA     Bs/B = K1     B + Bs + ABs = CB
    LOD condition:
      (CB*E0/(1+K1) + CB*K1/(1+K1))*(1+sigma) = (B*E0 + Bs + ABs)*(1-sigma)
    solved for CA and minimised over K1; E0 = E_in/E_act.
  At start-up the script checks that E0 = 1e-4, CB = 200 pM reproduces
  matlab/Allostery.mat (the last pass of the E0 loop in the MATLAB script).

PAG   : SI Eq. 2.22, r1*[1-exp(-lam t)]
SPEAR : SI Eq. 2.23, r1^2*[1-(1+lam t)exp(-lam t)]
  lam = k_uni + K_D k_on + K_D2 k_on2, r1 = k_uni/lam.
  Binding occupancy at the 2 uL binding step (4 x C_B); extension and background
  at the 8 uL stage. Dashed SPEAR curve: occupancy after a 1 h incubation instead
  of equilibrium, p(t) = p_eq*(1-exp(-k_on*(C+K_D)*t)) per site.

Parameters (published figure)
  t_ext = 8 min, C_B = C_B2 = 100 pM at 8 uL (400 pM at 2 uL), allostery
  C_B = 100 pM, E0 = 0.01, sigma = 0.1, k_uni = 0.01 s-1, k_bi = 5e3 M-1s-1,
  k_on = k_on2 = 5e5 M-1s-1.

Usage
  python fig_S5_lod_vs_kd.py            -> ../figures/Fig_S5_lod_vs_kd.png
  Optional positional arguments: probe pM at 8 uL (100), extension min (8),
  E_in/E_act (0.01), allostery binder pM (100).
"""


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import sympy as sp
import scipy.io as sio
import os, sys

VER = 'na01'
CB_PM      = float(sys.argv[1]) if len(sys.argv) > 1 else 100.0   # probe, pM at 8 uL
T_EXT_MIN  = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0     # extension, min
E0_ALLO    = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01    # E_in/E_act
CB_ALLO_PM = float(sys.argv[4]) if len(sys.argv) > 4 else 100.0   # allostery binder, pM

# ---------------------------------------------------------------- paths, fonts
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
C_TEAL, C_NAVY, C_GRID = '#1795BB', '#0B3954', '#E8E8E8'
C_TERRA, C_STEEL, C_TEAL_LT = '#B85042', '#465C8B', '#7FC6DC'

# ---------------------------------------------------------------- allostery, from scratch
_A, _B, _Bs, _ABs, _KD, _K1, _CA, _CB, _s, _E0 = sp.symbols(
    'A B Bs ABs KD K1 CA CB sigma E0', positive=True)
_eqs = [sp.Eq(_ABs/(_A*_Bs), 1/_KD), sp.Eq(_A+_ABs, _CA), sp.Eq(_Bs/_B, _K1),
        sp.Eq(_B+_Bs+_ABs, _CB),
        sp.Eq((_CB*_E0/(1+_K1) + _CB*_K1/(1+_K1))*(1+_s), (_B*_E0 + _Bs + _ABs)*(1-_s))]
_sol = sp.solve(_eqs, [_A, _B, _Bs, _ABs, _CA], dict=True)[0][_CA]
_lod_allo = sp.lambdify((_CB, _E0, _K1, _KD, _s), sp.simplify(_sol), 'numpy')
print('allostery LOD(CA) =', sp.simplify(_sol))

def allostery_lod(KD_arr, CB_allo, E0, sigma, n_K1=2000):
    """min over K1 of the closed-form LOD, as Script1 does."""
    K1 = np.logspace(-3, 0, n_K1)
    return np.array([np.min(_lod_allo(CB_allo, E0, K1, kd, sigma)) for kd in KD_arr])

# --- start-up check: reproduce Allostery.mat with ITS parameters
_m = sio.loadmat(os.path.join(HERE, '..', 'matlab', 'Allostery.mat'))
_kd_ref, _lod_ref = _m['vfKD_fine'].ravel(), _m['vfLOD_AllosteryFine'].ravel()
_rep = allostery_lod(_kd_ref, 0.2e-9, 1e-4, 0.1, n_K1=50)   # MATLAB: logspace(-3,0) = 50 pts
_err = np.max(np.abs(_rep/_lod_ref - 1))
assert _err < 1e-9, _err
print(f'independent re-derivation reproduces Allostery.mat (E0=1e-4, C_B=200 pM) '
      f'to {_err:.1e} relative')

# ---------------------------------------------------------------- SPEAR / proximity
k_uni, k_bi, k_on, sigma = 0.01, 5000.0, 5e5, 0.1
CB    = CB_PM * 1e-12          # at 8 uL
CB_b  = 4 * CB                 # at 2 uL
t_ext = T_EXT_MIN * 60.0
t_inc = 3600.0
NA, V = 6.022e23, 8e-6

KD  = np.logspace(-13, -7, 400)
lam = k_uni + 2 * KD * k_on
r1  = k_uni / lam
f2 = lambda x: -np.expm1(-x) - x * np.exp(-x)
f1 = lambda x: -np.expm1(-x)

p_eq  = CB_b / (CB_b + KD)
p_kin = p_eq * (1 - np.exp(-k_on * (CB_b + KD) * t_inc))
pref  = 2 * sigma / (1 - sigma) * CB

LOD_spear_eq  = pref * f2(CB*k_bi*t_ext) / (p_eq**2  * r1**2 * f2(lam*t_ext))
LOD_spear_kin = pref * f2(CB*k_bi*t_ext) / (p_kin**2 * r1**2 * f2(lam*t_ext))
LOD_pxa       = pref * f1(CB*k_bi*t_ext) / (p_eq**2  * r1    * f1(lam*t_ext))
LOD_allo      = allostery_lod(KD, CB_ALLO_PM*1e-12, E0_ALLO, sigma)

assert np.all(np.diff(LOD_spear_eq) > 0) and np.all(np.diff(LOD_allo) > 0)
assert abs((r1**2*f2(lam*t_ext))[0] - f2(k_uni*t_ext)) < 1e-4     # k_off->0 -> Eq. 2.17

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(figsize=(7.6, 5.6))
ax.axvspan(100e-12, 10e-9, color=C_TEAL, alpha=0.06, lw=0, zorder=0)
y100 = 100 / NA / V
ax.axhline(y100, color='#C9C9C9', lw=0.8, zorder=1)
ax.text(9.6e-8, y100*1.9, '100 molecules in 8 $\\mu$L', fontsize=9,
        color='#9A9A9A', ha='right', va='bottom', zorder=7)

ax.plot(KD, KD, ls=(0, (1, 3)), color='#BFBFBF', lw=1.1, zorder=2, label='$y = x$')
ax.plot(KD, LOD_allo, color=C_TERRA, lw=2.4, zorder=5, label='Allostery-based')
ax.plot(KD, LOD_pxa, color=C_STEEL, lw=2.4, zorder=5, label='Proximity Amplicon-generation')
ax.plot(KD, LOD_spear_kin, color=C_TEAL_LT, lw=2.2, ls=(0, (5, 2.5)), zorder=6,
        label='SPEAR (incubation only 1 h)')
ax.plot(KD, LOD_spear_eq, color=C_TEAL, lw=3.0, zorder=7,
        label='SPEAR (incubation achieving equilibrium)')

ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(1e-13, 1e-7); ax.set_ylim(1e-20, 1e-7)
ax.set_yticks(np.logspace(-20, -8, 13))
ax.grid(True, which='major', color=C_GRID, lw=0.5, zorder=0); ax.set_axisbelow(True)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('#444444'); ax.spines['bottom'].set_color('#444444')
ax.set_xlabel('$K_\\mathrm{D}$ = $K_\\mathrm{D2}$ (M)', fontsize=14, color=C_NAVY)
ax.set_ylabel('LOD (M, at 8 $\\mu$L)', fontsize=14, color=C_NAVY)
ax.set_title('Achievable LOD vs. antibody affinity', fontsize=15,
             fontweight='bold', color=C_NAVY, pad=44)
ax.text(0.5, 1.018,
        f'$C_\\mathrm{{B}}$ = $C_\\mathrm{{B2}}$ = {CB_PM:.0f} pM at 8 $\\mu$L '
        f'({CB_PM*4:.0f} pM at 2 $\\mu$L)  ·  allostery $C_\\mathrm{{B}}$ = {CB_ALLO_PM:.0f} pM, '
        f'$E_\\mathrm{{in}}/E_\\mathrm{{act}}$ = {E0_ALLO:g}\n'
        '$k_\\mathrm{uni}$ = 0.01 s$^{-1}$  ·  $k_\\mathrm{bi}$ = 5$\\times$10$^{3}$ M$^{-1}$s$^{-1}$  ·  '
        '$k_\\mathrm{on}$ = $k_\\mathrm{on2}$ = 5$\\times$10$^{5}$ M$^{-1}$s$^{-1}$  ·  '
        f'$\\sigma$ = 0.1  ·  $t$ = {T_EXT_MIN:.0f} min',
        transform=ax.transAxes, ha='center', va='bottom',
        fontsize=9.5, color='#6A7A85', linespacing=1.5)
ax.text(1e-9, 1.9e-20, 'typical antibody $K_\\mathrm{D}$  (100 pM – 10 nM)',
        fontsize=9.5, color='#4E8FA6', ha='center', va='bottom', style='italic', zorder=7)
leg = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.155), ncol=2, fontsize=10.5,
                frameon=True, facecolor='white', edgecolor='#C4CDD3', framealpha=1.0,
                borderpad=0.85, handlelength=2.6, columnspacing=2.4, labelspacing=0.55)
leg.get_frame().set_linewidth(0.9)

out = os.path.join(OUTDIR, 'Fig_S5_lod_vs_kd.png' if len(sys.argv) == 1
                   else f'Fig_S5_lod_vs_kd_{CB_PM:g}pM_{T_EXT_MIN:g}min.png')
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print(f'Saved: {out}  ({os.path.getsize(out)/1024:.0f} KB)\n')

print(f'{"KD":>8} {"allostery":>11} {"prox-amp":>11} {"SPEAR eq":>11} {"SPEAR 1h":>11} {"mol/8uL(1h)":>12}')
for kd in [1e-10, 2e-10, 5e-10, 1e-9, 5e-9, 1e-8]:
    i = np.argmin(abs(KD - kd))
    lbl = f'{kd*1e12:.0f} pM' if kd < 1e-9 else f'{kd*1e9:.0f} nM'
    print(f'{lbl:>8} {LOD_allo[i]*1e15:10.3g}f {LOD_pxa[i]*1e15:10.3g}f '
          f'{LOD_spear_eq[i]*1e15:10.4f}f {LOD_spear_kin[i]*1e15:10.4f}f '
          f'{LOD_spear_kin[i]*V*NA:12.0f}')
