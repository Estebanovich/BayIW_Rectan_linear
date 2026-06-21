#!/usr/bin/env python3
"""
plot_pile_eta.py [rundir]
-------------------------
Visualiza la senial BAROTROPICA (superficie libre Eta) de la corrida de columna
apilada: como la pila inicial se relaja y radia ondas.

Genera (para rundir, por defecto ../run_expand_pile):
  1) serie temporal de eta_max y eta en la cabeza de la bahia;
  2) mapas de Eta en 4 instantes (zoom a la region de la bahia);
  3) Hovmoller de Eta a lo largo del eje de la bahia (transecto cross-shelf).

Salida: <tag>_eta.png   (tag = nombre del run dir)

Uso (en el cluster):
  module load herramientas/python/3.11.8
  python3 plot_pile_eta.py ../run_expand_pile
  python3 plot_pile_eta.py ../run_expand_nobay_pile
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mitgcm_pile_io as io

RUNDIR = sys.argv[1] if len(sys.argv) > 1 else '../run_expand_pile'
TAG = os.path.basename(os.path.normpath(RUNDIR))
ICOL = 280                                   # columna central de la banda de la bahia

# grid acorde al caso (bahia / nobay)
nob = 'nobay' in TAG
g = io.load_grid(*( ('nobahia_01_expand_dx.bin','nobahia_01_expand_dy.bin','nobahia_01_expand_bat.bin')
                    if nob else ('bahia_01_expand_dx.bin','bahia_01_expand_dy.bin','bahia_01_expand_bat.bin') ))
xc, yc, depth = g['xc'], g['yc'], g['depth']

out = io.MITgcmOut(RUNDIR)
t, E = out.eta_all()
th = t / 3600.0
Emm = E * 1e3                                  # mm
land = io.land2d(g)
Emm_m = np.ma.masked_where(np.broadcast_to(land, Emm.shape), Emm)

# columna de la pila (donde hubo IC>0) para localizar cabeza/region
wetcol = np.where(depth[:, ICOL] < 0)[0]
jhead = wetcol.max()                           # costa/cabeza en esa columna
# region de la pila (donde eta(t=0)>0) + campo cercano mar adentro
OFFSHORE = 60                                  # celdas offshore para ver radiacion (~300 km)
pile_rows = np.where(Emm[0, :, ICOL] > 0.05)[0]   # umbral 0.05 mm
jp0, jp1 = pile_rows.min(), pile_rows.max()

# zoom: region de la pila (cabeza) hacia mar adentro
j0 = max(min(jp0, jp1) - OFFSHORE, 0)
j1 = min(max(jp0, jp1) + 6, io.NY)
i0, i1 = ICOL - 35, ICOL + 35
vmax = np.nanmax(np.abs(Emm[0]))               # escala con la pila inicial

fig = plt.figure(figsize=(16, 9))
gs = fig.add_gridspec(2, 4)

# --- (1) serie temporal ---
ax = fig.add_subplot(gs[0, 0:2])
ax.plot(th, np.nanmax(np.abs(Emm), axis=(1, 2)), '-', label='|eta|_max dominio')
ax.plot(th, Emm[:, jhead, ICOL], '-', label='eta en la cabeza (x=%.0f km)' % xc[ICOL])
ax.axhline(0, color='0.7', lw=0.6)
ax.set_xlabel('t / h'); ax.set_ylabel('eta / mm')
ax.set_title('relajacion de la pila'); ax.legend(fontsize=8); ax.grid(alpha=0.3)

# --- (2) Hovmoller a lo largo del eje de la bahia ---
ax = fig.add_subplot(gs[0, 2:4])
_, Hov = out.hovmoller_eta(ICOL)               # (nt, NY)
Hovmm = np.ma.masked_where(np.broadcast_to(land[:, ICOL], Hov.shape), Hov * 1e3)
pc = ax.pcolormesh(yc[j0:j1], th, Hovmm[:, j0:j1], cmap='RdBu_r',
                   vmin=-vmax, vmax=vmax, shading='auto')
ax.set_xlabel('y / km (cross-shelf; costa a la derecha)'); ax.set_ylabel('t / h')
ax.set_title('Hovmoller eta en x=%.0f km' % xc[ICOL])
plt.colorbar(pc, ax=ax, label='eta / mm')

# --- (3) mapas de Eta en 4 instantes (zoom) ---
idxs = [0, out.nt // 6, out.nt // 3, out.nt - 1]
for k, it in enumerate(idxs):
    ax = fig.add_subplot(gs[1, k])
    pc = ax.pcolormesh(xc[i0:i1], yc[j0:j1], Emm_m[it, j0:j1, i0:i1],
                       cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='auto')
    ax.contour(xc[i0:i1], yc[j0:j1], depth[j0:j1, i0:i1], levels=[0], colors='k', linewidths=0.8)
    ax.set_title('t=%.1f h' % th[it], fontsize=9)
    ax.set_xlabel('x/km', fontsize=8); ax.set_ylabel('y/km', fontsize=8)
    plt.colorbar(pc, ax=ax, label='eta/mm')

fig.suptitle('Senial barotropica (Eta) — %s  [%d snapshots, hasta %.1f h]'
             % (TAG, out.nt, th[-1]), fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fn = '%s_eta.png' % TAG
fig.savefig(fn, dpi=140)
print('escrito:', fn)
