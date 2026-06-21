#!/usr/bin/env python3
"""
plot_pile_section.py [rundir]
-----------------------------
Visualiza la senial BAROCLINICA (ondas internas) de la corrida de columna apilada:
secciones verticales de temperatura/densidad a lo largo del eje de la bahia, donde
el desplazamiento de las isotermas/isopicnas revela las ondas internas generadas
por la relajacion de la pila.

Genera (para rundir, por defecto ../run_expand_pile):
  fila 1: seccion de Temperatura (con isotermas) en 4 instantes;
  fila 2: anomalia de Temperatura respecto al estado inicial (resalta las ondas).

Salida: <tag>_section.png

Uso (en el cluster):
  module load herramientas/python/3.11.8
  python3 plot_pile_section.py ../run_expand_pile
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mitgcm_pile_io as io

RUNDIR = sys.argv[1] if len(sys.argv) > 1 else '../run_expand_pile'
TAG = os.path.basename(os.path.normpath(RUNDIR))
ICOL = 280                                   # eje de la bahia (cross-shelf)

nob = 'nobay' in TAG
g = io.load_grid(*( ('nobahia_01_expand_dx.bin','nobahia_01_expand_dy.bin','nobahia_01_expand_bat.bin')
                    if nob else ('bahia_01_expand_dx.bin','bahia_01_expand_dy.bin','bahia_01_expand_bat.bin') ))
yc, zc, depth = g['yc'], g['zc'], g['depth']

out = io.MITgcmOut(RUNDIR)
th = out.times / 3600.0

# region de interes en y: la pila (eta t=0>0) + campo cercano mar adentro
OFFSHORE = 60                                  # celdas offshore (~300 km)
e0 = out.var_time('Eta', 0)[:, ICOL]
pile_rows = np.where(e0 > 5e-5)[0]             # umbral 0.05 mm
j0 = max(pile_rows.min() - OFFSHORE, 0)
j1 = min(pile_rows.max() + 6, io.NY)
Y = yc[j0:j1]

# mascara de tierra/fondo para la columna (NR, j)
wet = io.wet3d(g)[:, j0:j1, ICOL]              # (NR, nj)

def section(var, it):
    s = out.transect_time(var, icol=ICOL, tidx=it)[:, j0:j1]   # (NR, nj)
    return np.ma.masked_where(~wet, s)

idxs = [0, out.nt // 6, out.nt // 3, out.nt - 1]
sec0 = section('Temp', 0)
isos = np.arange(4, 34, 2.0)                    # isotermas cada 2 C

fig, axes = plt.subplots(2, 4, figsize=(17, 8), sharex=True, sharey=True,
                         constrained_layout=True)
# escala comun para la anomalia
anoms = [section('Temp', it) - sec0 for it in idxs]
amax = max(float(np.nanmax(np.abs(a))) for a in anoms[1:]) or 1e-3

for k, it in enumerate(idxs):
    T = section('Temp', it)
    # fila 1: temperatura + isotermas
    ax = axes[0, k]
    pc = ax.pcolormesh(Y, zc, T, cmap='turbo', vmin=4, vmax=34, shading='auto')
    cs = ax.contour(Y, zc, T, levels=isos, colors='k', linewidths=0.4)
    ax.set_title('t=%.1f h' % th[it], fontsize=10)
    if k == 0: ax.set_ylabel('z / m')
    if k == 3: plt.colorbar(pc, ax=axes[0, :], label='T / C', shrink=0.8, pad=0.01)

    # fila 2: anomalia respecto a t=0
    ax = axes[1, k]
    pa = ax.pcolormesh(Y, zc, anoms[k], cmap='RdBu_r', vmin=-amax, vmax=amax, shading='auto')
    ax.contour(Y, zc, section('Temp', it), levels=isos, colors='0.3', linewidths=0.3)
    ax.set_xlabel('y / km (costa a la derecha)')
    if k == 0: ax.set_ylabel('z / m')
    if k == 3: plt.colorbar(pa, ax=axes[1, :], label="T' / C (vs t=0)", shrink=0.8, pad=0.01)

fig.suptitle('Seccion vertical en x=%.0f km (eje de la bahia) — %s\n'
             'fila 1: Temperatura+isotermas   fila 2: anomalia T (ondas internas)'
             % (g['xc'][ICOL], TAG), fontsize=12)
fn = '%s_section.png' % TAG
fig.savefig(fn, dpi=140)
print('escrito:', fn)
