"""
export_allostery_json.py
Writes the contents of matlab/Allostery.mat to matlab/Allostery.json, so the data
can be read without MATLAB. Values are written with full double precision
(Python repr), and the script checks that they read back bit-identical.

Usage: python export_allostery_json.py
"""
import json, os
import numpy as np
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', 'matlab', 'Allostery.mat')
DST = os.path.join(HERE, '..', 'matlab', 'Allostery.json')

m = sio.loadmat(SRC)
kd = m['vfKD_fine'].ravel()
lod = m['vfLOD_AllosteryFine'].ravel()

out = {
    'source': 'matlab/Allostery.mat, written by matlab/Script1_Allostery_250108.m',
    'description': ('Lowest achievable LOD of an allostery-based binder versus K_D, '
                    'minimised over K1 in logspace(-3, 0, 50). This is the last pass of '
                    'the E_in/E_act loop in the script that makes Fig. S10f.'),
    'parameters': {'C_B_M': 200e-12, 'E_in_over_E_act': 1e-4, 'sigma': 0.1},
    'variables': {
        'vfKD_fine': {'description': 'K_D of the binder in the active state',
                      'unit': 'M', 'values': kd.tolist()},
        'vfLOD_AllosteryFine': {'description': 'minimum LOD over K1',
                                'unit': 'M', 'values': lod.tolist()},
    },
}
with open(DST, 'w') as f:
    json.dump(out, f, indent=2)
    f.write('\n')

back = json.load(open(DST))['variables']
assert np.array_equal(np.array(back['vfKD_fine']['values']), kd)
assert np.array_equal(np.array(back['vfLOD_AllosteryFine']['values']), lod)
print(f'Saved: {DST}  ({len(kd)} points, round-trip exact)')
