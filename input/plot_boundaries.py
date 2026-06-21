#!/usr/bin/env python3
"""
plot_boundaries.py
==================
Visualiza las CONDICIONES DE FRONTERA ABIERTA (OBCS) sobre el dominio.

Dibuja la batimetria en la malla REAL (con margenes estirados) y sobrepone:
  - las fronteras abiertas (Sur, Oeste, Este) leidas de data.obcs,
  - la frontera Norte cerrada (costa),
  - la capa esponja (spongeThickness celdas hacia el interior),
  - el margen estirado (NSTRETCH celdas, make_stretched_grid.py).

Lee los binarios de batimetria y de espaciado (dx/dy) que ya estan en input/.
Genera un PNG por escenario (con bahia / sin bahia).

Uso:
    python3 plot_boundaries.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# ---- Parametros del setup (SIZE.h + data.obcs + make_stretched_grid.py) ----
NX, NY          = 560, 352      # tamano de la malla
SPONGE          = 30            # spongeThickness (data.obcs, OBCS_PARM03)
NSTRETCH        = 20            # celdas estiradas por margen (make_stretched_grid.py)
RELAX_BOUND     = 60.0          # Urelax/Vrelaxobcsbound (s) en la frontera
RELAX_INNER     = 1200.0        # Urelax/Vrelaxobcsinner (s) en el borde interior
DTYPE           = np.dtype('>f8')

# Fronteras abiertas (data.obcs, OBCS_PARM01)
#   OB_Jsouth=560*1 (abierta), OB_Jnorth=560*0 (cerrada/costa),
#   OB_Iwest=352*1 (abierta), OB_Ieast=352*-1 (abierta)
OPEN = {'sur': True, 'norte': False, 'oeste': True, 'este': True}

SCENARIOS = [
    ('bahia_01_expand',   'Con bahia'),
    ('nobahia_01_expand', 'Sin bahia'),
]


def cell_edges(delta):
    """Bordes de celda (N+1) a partir del vector de espaciados (N)."""
    e = np.zeros(delta.size + 1)
    e[1:] = np.cumsum(delta)
    return e


def tau_of_sl(sl):
    """Tiempo de relajacion (s) en la capa sl de la esponja (sl=1..SPONGE),
    interpolacion lineal de obcs_sponge.F:
        tau = ((SPONGE-sl)*bound + sl*inner)/SPONGE."""
    return ((SPONGE - sl) * RELAX_BOUND + sl * RELAX_INNER) / SPONGE


def relax_timescale_field():
    """Campo 2D [NY,NX] del tiempo de relajacion de la esponja (s).
    En cada celda toma el tau MAS CORTO (relajacion mas fuerte) entre las
    fronteras abiertas que la alcanzan; NaN fuera de la esponja."""
    tau = np.full((NY, NX), np.nan)
    jj, ii = np.indices((NY, NX))
    # distancia en celdas (sl) a cada frontera abierta
    sls = []
    if OPEN['sur']:
        sls.append(jj)                 # sl crece hacia el interior desde y=0
    if OPEN['norte']:
        sls.append((NY - 1) - jj)
    if OPEN['oeste']:
        sls.append(ii)
    if OPEN['este']:
        sls.append((NX - 1) - ii)
    for sl in sls:
        m = (sl >= 1) & (sl <= SPONGE)
        cand = np.where(m, tau_of_sl(sl), np.nan)
        tau = np.fmin(tau, cand)       # fmin ignora NaN -> tau mas corto gana
    return tau


def plot_scenario(prefix, title):
    bat = np.fromfile(f'{prefix}_bat.bin', DTYPE).reshape(NY, NX)  # [j, i]
    dx  = np.fromfile(f'{prefix}_dx.bin',  DTYPE)
    dy  = np.fromfile(f'{prefix}_dy.bin',  DTYPE)

    xe = cell_edges(dx) / 1e3   # km
    ye = cell_edges(dy) / 1e3   # km
    Lx, Ly = xe[-1], ye[-1]

    fig, ax = plt.subplots(figsize=(11, 8))

    # --- batimetria (tierra = bat>=0 en gris) ---
    bat_ma = np.ma.masked_where(bat >= 0.0, bat)
    pc = ax.pcolormesh(xe, ye, bat_ma, cmap='Blues_r', shading='flat')
    cb = fig.colorbar(pc, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label('Profundidad del fondo (m)')
    ax.pcolormesh(xe, ye, np.ma.masked_where(bat < 0.0, bat),
                  cmap='Greys', shading='flat', vmin=-1, vmax=1)

    # --- capa esponja: gradiente del tiempo de relajacion (s) ---
    tau = relax_timescale_field()
    tau_ma = np.ma.masked_invalid(tau)
    sp = ax.pcolormesh(xe, ye, tau_ma, cmap='autumn', shading='flat',
                       vmin=RELAX_BOUND, vmax=RELAX_INNER, zorder=2, alpha=0.85)
    cb2 = fig.colorbar(sp, ax=ax, shrink=0.8, pad=0.12, orientation='horizontal')
    cb2.set_label('Esponja: tiempo de relajacion tau (s)\n'
                  '(corto = amortiguamiento fuerte)')

    # --- margen estirado (celdas externas) lineas punteadas ---
    st_kw = dict(color='purple', ls=':', lw=1.4, zorder=3)
    if OPEN['sur']:
        ax.axhline(ye[NSTRETCH], **st_kw)
    if OPEN['oeste']:
        ax.axvline(xe[NSTRETCH], **st_kw)
    if OPEN['este']:
        ax.axvline(xe[NX - NSTRETCH], **st_kw)

    # --- fronteras (bordes del dominio) ---
    open_kw   = dict(color='red',   lw=4, zorder=5)
    closed_kw = dict(color='black', lw=5, zorder=5)
    # Sur (y=0)
    ax.plot([0, Lx], [0, 0], **(open_kw if OPEN['sur'] else closed_kw))
    # Norte (y=Ly) - costa cerrada
    ax.plot([0, Lx], [Ly, Ly], **(open_kw if OPEN['norte'] else closed_kw))
    # Oeste (x=0)
    ax.plot([0, 0], [0, Ly], **(open_kw if OPEN['oeste'] else closed_kw))
    # Este (x=Lx)
    ax.plot([Lx, Lx], [0, Ly], **(open_kw if OPEN['este'] else closed_kw))

    # --- etiquetas de cada frontera ---
    txt = dict(ha='center', va='center', fontsize=11, fontweight='bold')
    ax.text(Lx / 2, -Ly * 0.035, 'SUR (abierta, Orlanski)', color='red', **txt)
    ax.text(Lx / 2, Ly * 1.03, 'NORTE (cerrada / costa)', color='black', **txt)
    ax.text(-Lx * 0.04, Ly / 2, 'OESTE (abierta)', color='red', rotation=90, **txt)
    ax.text(Lx * 1.04, Ly / 2, 'ESTE (abierta)', color='red', rotation=90, **txt)

    ax.set_xlim(-Lx * 0.08, Lx * 1.08)
    ax.set_ylim(-Ly * 0.08, Ly * 1.08)
    ax.set_aspect('equal')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_title(f'Condiciones de frontera OBCS — {title}\n'
                 f'dominio {Lx:.0f}x{Ly:.0f} km  ({NX}x{NY} celdas)')

    legend = [
        Line2D([0], [0], color='red',   lw=4, label='Frontera abierta (radiacion Orlanski)'),
        Line2D([0], [0], color='black', lw=5, label='Frontera cerrada (costa)'),
        Patch(facecolor='red',    alpha=0.85, label=f'Esponja borde: tau~{tau_of_sl(1):.0f} s (fuerte)'),
        Patch(facecolor='yellow', alpha=0.85, label=f'Esponja interior: tau={tau_of_sl(SPONGE):.0f} s (debil)'),
        Line2D([0], [0], color='purple', ls=':', lw=1.6,
               label=f'Inicio margen estirado ({NSTRETCH} celdas)'),
    ]
    ax.legend(handles=legend, loc='upper right', fontsize=9, framealpha=0.95)

    fig.tight_layout()
    out = f'fronteras_{prefix}.png'
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f'  escrito {out}  (dominio {Lx:.0f}x{Ly:.0f} km)')


def main():
    print('Visualizando condiciones de frontera OBCS sobre el dominio...')
    for prefix, title in SCENARIOS:
        try:
            plot_scenario(prefix, title)
        except FileNotFoundError as e:
            print(f'  [omito {prefix}] falta archivo: {e.filename}')


if __name__ == '__main__':
    main()
