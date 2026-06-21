#!/usr/bin/env python3
"""
plot_initial_pile.py
--------------------
Inspeccion visual de los campos de columna apilada (pSurfInitFile) generados por
make_initial_pile.py:  pile_init_bahia_560x352.bin y pile_init_nobahia_560x352.bin.

Para cada caso (con bahia / sin bahia) grafica:
  (a) mapa de eta (mm) con la linea de costa (bathy=0) superpuesta, zoom a la
      region de la bahia;
  (b) mapa de eta en todo el dominio (contexto);
  (c) perfil cross-shelf de eta en la columna central de la banda.

Salida: pile_init_inspeccion.png
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NX, NY = 560, 352
DT = np.dtype('>f8')
XG, YG = -1400e3, -1621e3

dx = np.fromfile('bahia_01_expand_dx.bin', DT)
dy = np.fromfile('bahia_01_expand_dy.bin', DT)
xc = (XG + np.cumsum(dx) - dx / 2) / 1e3        # km
yc = (YG + np.cumsum(dy) - dy / 2) / 1e3        # km

cases = [
    ('CON bahia',  'bahia_01_expand_bat.bin',   'pile_init_bahia_560x352.bin'),
    ('SIN bahia',  'nobahia_01_expand_bat.bin', 'pile_init_nobahia_560x352.bin'),
]

fig, axes = plt.subplots(2, 3, figsize=(16, 9))
emax_mm = 8.9

for r, (name, batf, pilef) in enumerate(cases):
    bat = np.fromfile(batf, DT).reshape(NY, NX)
    eta = np.fromfile(pilef, DT).reshape(NY, NX) * 1e3      # mm
    etam = np.ma.masked_where(eta <= 0, eta)

    # banda activa (columnas con eta>0) y su centro
    cols = np.where((eta > 0).any(axis=0))[0]
    rows = np.where((eta > 0).any(axis=1))[0]
    ic = int(round(cols.mean()))

    # --- (a) zoom a la region de la pila ---
    ax = axes[r, 0]
    i0, i1 = cols.min() - 8, cols.max() + 9
    j0, j1 = rows.min() - 8, rows.max() + 9
    pc = ax.pcolormesh(xc[i0:i1], yc[j0:j1], etam[j0:j1, i0:i1],
                       cmap='viridis', vmin=0, vmax=emax_mm, shading='auto')
    ax.contour(xc[i0:i1], yc[j0:j1], bat[j0:j1, i0:i1], levels=[0],
               colors='k', linewidths=1.2)
    ax.contour(xc[i0:i1], yc[j0:j1], bat[j0:j1, i0:i1], levels=[-300, -160],
               colors='0.6', linewidths=0.6)
    plt.colorbar(pc, ax=ax, label='eta / mm')
    ax.set_title('%s — zoom pila + costa' % name)
    ax.set_xlabel('x / km'); ax.set_ylabel('y / km (norte=costa)')
    ax.axvline(xc[ic], color='r', ls='--', lw=0.8)

    # --- (b) dominio completo ---
    ax = axes[r, 1]
    pc = ax.pcolormesh(xc, yc, etam, cmap='viridis', vmin=0, vmax=emax_mm, shading='auto')
    ax.contour(xc, yc, bat, levels=[0], colors='k', linewidths=0.5)
    plt.colorbar(pc, ax=ax, label='eta / mm')
    ax.set_title('%s — dominio completo' % name)
    ax.set_xlabel('x / km'); ax.set_ylabel('y / km')

    # --- (c) perfil cross-shelf en la columna central ---
    ax = axes[r, 2]
    col = eta[:, ic]
    m = col > 0
    ax.plot(col[m], yc[m], '-o', ms=3, color='C0')
    ax.set_title('%s — perfil eta(y) en x=%.0f km' % (name, xc[ic]))
    ax.set_xlabel('eta / mm'); ax.set_ylabel('y / km')
    ax.grid(alpha=0.3)
    jh = np.where(m)[0]
    ax.annotate('cabeza/costa\n(eta_max=%.2f mm)' % col[jh].max(),
                xy=(col[jh[-1]], yc[jh[-1]]), fontsize=8, color='r')
    ax.annotate('boca/mar adentro', xy=(col[jh[0]], yc[jh[0]]), fontsize=8, color='r')

fig.suptitle('Condicion inicial de columna apilada (pSurfInitFile) — '
             'eta_max=8.9 mm = set-up del viento', fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = 'pile_init_inspeccion.png'
fig.savefig(out, dpi=140)
print('escrito:', out)
