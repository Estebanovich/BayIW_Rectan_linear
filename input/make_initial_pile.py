#!/usr/bin/env python3
"""
make_initial_pile.py
--------------------
Genera la CONDICION INICIAL de "columna de agua apilada" (pSurfInitFile) para los
experimentos BayIW_Rectan_linear, en sustitucion del forzamiento de viento.

En lugar de apilar agua contra la cabeza de la bahia con viento meridional, se
arranca la simulacion con una elevacion de superficie libre (eta) ya apilada que
crece LINEALMENTE desde la boca de la bahia (eta=0) hasta la cabeza (eta=ETA_MAX).
Al soltarse (sin forzante), esta pila desencadena el ajuste y las ondas internas.

Definicion (unificada para los casos CON y SIN bahia):
  * Banda en x = footprint de la bahia (columnas i donde 'bahia' tiene agua y
    'nobahia' es tierra), detectada automaticamente de la diferencia batimetrica.
  * Para cada columna i de la banda se localiza la LINEA DE COSTA = celda humeda
    mas al norte (mayor j con bathy<0).
  * eta(i,j) = ETA_MAX * max(0, 1 - d/L),  con d = distancia cross-shelf (m) desde
    la costa hacia mar abierto (j decreciente) y L = LARGO_BAHIA.
  * Solo en celdas humedas (bathy<0); cero fuera de la banda y en tierra.

  -> Caso bahia : la rampa llena el canal, ETA_MAX en la cabeza (norte), ~0 en la
                  boca (sur).  "crece linealmente de boca a cabeza".
  -> Caso nobay : misma banda y amplitud, anclada a la costa recta (la cabeza no
                  existe) -> control fisicamente valido (pila contra la costa).

ETA_MAX se iguala al set-up maximo del viento medido en la corrida forzada
run_expand (etapa 1).  Valor medido: |Eta|_max = 0.0089 m, localizado EN LA CABEZA
de la bahia (j=347, i=285).  Si MEASURE_ETA=True el script lo vuelve a medir de
run_expand/OUT_stage1/mnc_*/state.0000000000.t*.nc (requiere netCDF4).

Salida (big-endian float64, layout (NY,NX) en orden C == (Nx,Ny) Fortran de MITgcm,
misma convencion que bahia_01_expand_bat.bin):
  * pile_init_bahia_560x352.bin     (bathyFile = bahia_01_expand_bat.bin)
  * pile_init_nobahia_560x352.bin   (bathyFile = nobahia_01_expand_bat.bin)
"""
import numpy as np

# --------------------------- parametros ------------------------------------
NX, NY = 560, 352
DT = np.dtype('>f8')                 # big-endian float64 (nativo MITgcm)
XG, YG = -1400e3, -1621e3            # xgOrigin, ygOrigin (data &PARM04)
LARGO_BAHIA = 120e3                  # L: largo de la bahia (m) = escala de la rampa

ETA_MAX = 0.0089                     # m, set-up maximo del viento (medido en run_expand)
MEASURE_ETA = False                  # True -> re-medir de la salida forzada (netCDF4)

BAT_BAHIA  = 'bahia_01_expand_bat.bin'
BAT_NOBAY  = 'nobahia_01_expand_bat.bin'
DXF, DYF   = 'bahia_01_expand_dx.bin', 'bahia_01_expand_dy.bin'
OUT_BAHIA  = 'pile_init_bahia_560x352.bin'
OUT_NOBAY  = 'pile_init_nobahia_560x352.bin'


def measure_wind_eta():
    """Pico de |Eta| (m) al final del set-up del viento en run_expand, etapa 1."""
    import glob
    from netCDF4 import Dataset
    dx = np.fromfile(DXF, DT); dy = np.fromfile(DYF, DT)
    xc = XG + np.cumsum(dx) - dx / 2
    yc = YG + np.cumsum(dy) - dy / 2
    files = sorted(glob.glob('../run_expand/OUT_stage1/mnc_*/state.0000000000.t*.nc'))
    G = None
    for f in files:
        with Dataset(f) as nc:
            X = np.array(nc.variables['X'][:]); Y = np.array(nc.variables['Y'][:])
            E = np.array(nc.variables['Eta'][:])           # (T,Y,X)
            if G is None:
                G = np.zeros((E.shape[0], NY, NX))
            i0 = int(np.argmin(np.abs(xc - X[0])))
            j0 = int(np.argmin(np.abs(yc - Y[0])))
            G[:, j0:j0 + len(Y), i0:i0 + len(X)] = E
    return float(np.nanmax(np.abs(G)))


def build_pile(bathy, dy, eta_max, band):
    """Campo eta(NY,NX): rampa lineal anclada a la costa de cada columna de 'band'."""
    yc = YG + np.cumsum(dy) - dy / 2                        # centros de celda en y (m)
    eta = np.zeros((NY, NX))
    for i in band:
        wet = np.where(bathy[:, i] < 0)[0]                  # celdas humedas de la columna
        if wet.size == 0:
            continue
        jc = wet.max()                                      # costa = celda humeda mas al norte
        d = yc[jc] - yc[wet]                                # distancia cross-shelf (>=0 mar adentro)
        ramp = eta_max * np.clip(1.0 - d / LARGO_BAHIA, 0.0, None)
        eta[wet, i] = ramp
    return eta


def main():
    bat  = np.fromfile(BAT_BAHIA, DT).reshape(NY, NX)
    nob  = np.fromfile(BAT_NOBAY, DT).reshape(NY, NX)
    dy   = np.fromfile(DYF, DT)

    eta_max = measure_wind_eta() if MEASURE_ETA else ETA_MAX
    print('ETA_MAX usado = %.5f m  (%s)'
          % (eta_max, 'medido de run_expand' if MEASURE_ETA else 'constante set-up viento'))

    # banda en x = footprint de la bahia (agua en 'bahia', tierra en 'nobahia')
    bay_mask = (bat < 0) & (nob >= 0)
    band = np.where(bay_mask.any(axis=0))[0]
    print('banda bahia (x): i=%d..%d  (%d columnas)' % (band.min(), band.max(), band.size))

    for tag, bathy, out in (('bahia', bat, OUT_BAHIA), ('nobay', nob, OUT_NOBAY)):
        eta = build_pile(bathy, dy, eta_max, band)
        eta.astype(DT).tofile(out)
        nwet = int(np.count_nonzero(eta > 0))
        # perfil cross-shelf en la columna central de la banda
        ic = int(round(band.mean()))
        col = eta[:, ic]; jnz = np.where(col > 0)[0]
        rng = ('j=%d..%d' % (jnz.min(), jnz.max())) if jnz.size else '(seca)'
        print('  %-6s -> %s  eta_max=%.5f  celdas>0=%d  columna i=%d activa en %s'
              % (tag, out, eta.max(), nwet, ic, rng))

    print('bytes/archivo esperados = %d' % (NX * NY * 8))


if __name__ == '__main__':
    main()
