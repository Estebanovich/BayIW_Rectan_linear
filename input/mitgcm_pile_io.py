#!/usr/bin/env python3
"""
mitgcm_pile_io.py
-----------------
Utilidades para leer y ensamblar la salida MNC (NetCDF en tiles) de los
experimentos de columna apilada (run_expand_pile / run_expand_nobay_pile).

La salida MITgcm con MNC se escribe en 20 tiles (5x4) bajo mnc_*/ (o OUT_pile/mnc_*/
si se archivo con el SLURM). Cada tile cubre 112x88 puntos y guarda todos los
registros de tiempo de la corrida. Este modulo:
  * localiza los tiles y los recoloca en la malla global (560x352) usando las
    coordenadas X,Y de cada archivo;
  * entrega Eta(t) completa (2D barato) o campos 3D en un instante dado;
  * provee la malla (xc,yc en km, profundidad Z de las capas, batimetria).

Uso tipico (ver plot_pile_eta.py / plot_pile_section.py):
    import mitgcm_pile_io as io
    g   = io.load_grid()
    out = io.MITgcmOut('../run_expand_pile')       # o ruta a OUT_pile
    t, eta = out.eta_all()                          # (nt,), (nt,NY,NX)
    temp   = out.var_time('Temp', tidx=-1)          # (NR,NY,NX) en el ultimo paso
"""
import os, glob, re
import numpy as np
from netCDF4 import Dataset

NX, NY, NR = 560, 352, 50
XG, YG = -1400e3, -1621e3
DT = np.dtype('>f8')
HERE = os.path.dirname(os.path.abspath(__file__))     # .../input


def load_grid(dxf='bahia_01_expand_dx.bin', dyf='bahia_01_expand_dy.bin',
              batf='bahia_01_expand_bat.bin', delR=None):
    """Devuelve dict con xc,yc (km), depth(NY,NX, <0 agua), zc (m, centros de capa)."""
    dx = np.fromfile(os.path.join(HERE, dxf), DT)
    dy = np.fromfile(os.path.join(HERE, dyf), DT)
    depth = np.fromfile(os.path.join(HERE, batf), DT).reshape(NY, NX)
    xc = (XG + np.cumsum(dx) - dx / 2) / 1e3          # km
    yc = (YG + np.cumsum(dy) - dy / 2) / 1e3          # km
    if delR is None:
        delR = np.array([1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.9,2.0,2.2,2.4,2.6,2.8,
                         3.0,3.2,3.5,3.8,4.1,4.4,4.8,5.2,5.6,6.0,6.5,7.1,7.6,8.2,
                         8.9,9.6,10.4,11.3,12.2,13.2,14.3,15.4,16.7,18.0,19.5,21.1,
                         22.8,24.6,26.6,28.8,31.1,33.7,36.4,39.3,42.5,45.8])
    zf = np.concatenate(([0.0], -np.cumsum(delR)))    # interfaces (m, <0)
    zc = 0.5 * (zf[:-1] + zf[1:])                     # centros de capa (m)
    return dict(xc=xc, yc=yc, depth=depth, zc=zc, zf=zf, dx=dx, dy=dy)


def land2d(grid):
    """Mascara (NY,NX) True en tierra (depth>=0)."""
    return grid['depth'] >= 0


def wet3d(grid):
    """Mascara (NR,NY,NX) True donde la celda es agua (centro por encima del fondo)."""
    zc = grid['zc'][:, None, None]            # (NR,1,1)
    return zc > grid['depth'][None, :, :]     # agua si centro de capa > batimetria (<0)


class MITgcmOut:
    """Lector de salida state.* en tiles, recolocada a malla global."""

    def __init__(self, rundir, prefix='state'):
        self.rundir = rundir
        # acepta mnc_* directo en rundir, o archivados en OUT_pile/mnc_*
        pats = [os.path.join(rundir, 'mnc_*', f'{prefix}.*.nc'),
                os.path.join(rundir, 'OUT_pile', 'mnc_*', f'{prefix}.*.nc')]
        self.files = []
        for p in pats:
            self.files = sorted(glob.glob(p))
            if self.files:
                break
        if not self.files:
            raise FileNotFoundError(f'sin archivos {prefix}.* bajo {rundir} (mnc_* / OUT_pile)')
        # un segmento por simplicidad (corrida desde cero). Si hay varios, toma el de iter0 menor.
        segs = sorted(set(os.path.basename(f).split('.t')[0] for f in self.files))
        self.seg = segs[0]
        self.files = [f for f in self.files if os.path.basename(f).startswith(self.seg + '.t')]
        if len(segs) > 1:
            print(f'[aviso] varios segmentos {segs}; usando {self.seg}. '
                  f'(concatenar segmentos no implementado)')
        # mapa tile -> (i0, j0) y eje temporal
        self.tiles = {}
        self._place_tiles()
        # nt comun = minimo de registros entre tiles (la corrida puede estar escribiendo EN VIVO)
        nmin = min(n for (_, _, _, _, n) in self.tiles.values())
        with Dataset(self.files[0]) as nc:
            self.times = np.array(nc.variables['T'][:nmin])  # s
            self.iters = np.array(nc.variables['iter'][:nmin])
        self._nt = nmin

    def _place_tiles(self):
        g = load_grid()
        xc, yc = g['xc'] * 1e3, g['yc'] * 1e3
        for f in self.files:
            with Dataset(f) as nc:
                X = np.array(nc.variables['X'][:]); Y = np.array(nc.variables['Y'][:])
                ntf = nc.dimensions['T'].size if 'T' in nc.dimensions else len(nc.variables['T'])
            i0 = int(np.argmin(np.abs(xc - X[0]))); j0 = int(np.argmin(np.abs(yc - Y[0])))
            self.tiles[f] = (i0, j0, len(X), len(Y), ntf)

    @property
    def nt(self):
        return self._nt

    def _abs_idx(self, tidx):
        """Convierte un indice (posiblemente negativo) a absoluto dentro de nt comun."""
        return tidx % self._nt

    def eta_all(self):
        """(times[s], Eta[nt,NY,NX])."""
        G = np.full((self.nt, NY, NX), np.nan)
        for f, (i0, j0, nxx, nyy, _) in self.tiles.items():
            with Dataset(f) as nc:
                E = np.array(nc.variables['Eta'][:self._nt])  # (T,Y,X)
            G[:, j0:j0 + nyy, i0:i0 + nxx] = E
        return self.times.copy(), G

    def var_time(self, var, tidx=-1):
        """Campo global en UN instante. var 2D->(NY,NX); var 3D->(NR,NY,NX).
        U/V se interpolan a centros de celda."""
        ti = self._abs_idx(tidx)
        out = None
        for f, (i0, j0, nxx, nyy, _) in self.tiles.items():
            with Dataset(f) as nc:
                v = nc.variables[var]
                a = np.array(v[ti])                           # (...,Y,X) sin eje T
                dims = v.dimensions
            # interpolar staggered a centros
            if 'Xp1' in dims:                                  # U
                a = 0.5 * (a[..., :-1] + a[..., 1:])
            if 'Yp1' in dims:                                  # V
                a = 0.5 * (a[..., :-1, :] + a[..., 1:, :])
            if out is None:
                shp = (NR, NY, NX) if a.ndim == 3 else (NY, NX)
                out = np.full(shp, np.nan)
            out[..., j0:j0 + nyy, i0:i0 + nxx] = a
        return out

    def transect_time(self, var, icol, tidx=-1):
        """Seccion vertical (NR, NY) de 'var' en la columna x=icol, un instante."""
        f3 = self.var_time(var, tidx=tidx)
        return f3[:, :, icol]                                  # (NR, NY)

    def hovmoller_eta(self, icol):
        """(times, NY) de Eta a lo largo de la columna x=icol (transecto cross-shelf)."""
        t, E = self.eta_all()
        return t, E[:, :, icol]
