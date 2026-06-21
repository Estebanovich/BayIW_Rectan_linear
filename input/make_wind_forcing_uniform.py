#!/usr/bin/env python3
"""
make_wind_forcing_uniform.py
-----------------------------
Genera un forzamiento de viento meridional ESPACIALMENTE UNIFORME sobre TODO el
dominio (560x352), como contraparte del forzamiento LOCAL en parche
(merid_local_bay_560x352.bin) usado en los 3 escenarios BayIW_Rectan*.

Estrategia:
  * El forzamiento local es separable:  tau(t,y,x) = g(t) * S(y,x)
    (verificado: max|tau - g*S| ~ 1e-17).
  * Se REUTILIZA el pulso temporal g(t) EXACTO del archivo local (180 registros,
    periodo 1200 s, ciclo 216000 s), de modo que SOLO cambia la distribucion
    espacial: parche local  ->  uniforme en todo el dominio.
  * Amplitud renormalizada para igualar la ENERGIA inyectada.

Renormalizacion (NORM):
  - 'energy'   (defecto): conserva  E = ∬|tau|^2 dA   (proxy de trabajo del viento
                          tau·u en respuesta lineal).  A_uni = RMS_area(S_local).
  - 'momentum'          : conserva  I = ∬|tau|   dA   (impulso/momento total).
                          A_uni = media_area(S_local).
  NOTA: la energia REAL del viento es ∬tau·u dA (necesita velocidades del modelo);
  solo puede verificarse A POSTERIORI con compute_wind_energy.py + las salidas.

Salida:  merid_uniform_560x352.bin  (big-endian float64, layout (nt,ny,nx) en C,
         identico al del archivo local -> compatible con meridWindFile de MITgcm).
"""
import numpy as np

# --------------------------- parametros ------------------------------------
NX, NY, NT = 560, 352, 180
LOCAL = 'merid_local_bay_560x352.bin'
DX    = 'bahia_01_expand_dx.bin'
DY    = 'bahia_01_expand_dy.bin'
BATHY = 'bahia_01_expand_bat.bin'
OUT   = 'merid_uniform_560x352.bin'
NORM  = 'energy'     # 'energy' (∬|tau|^2 dA) o 'momentum' (∬|tau| dA)
DT    = np.dtype('>f8')

# --------------------------- leer local ------------------------------------
tau = np.fromfile(LOCAL, dtype=DT).reshape(NT, NY, NX)        # (t,y,x)
dx  = np.fromfile(DX, dtype=DT)                                # (nx,)
dy  = np.fromfile(DY, dtype=DT)                                # (ny,)
bathy = np.fromfile(BATHY, dtype=DT).reshape(NY, NX)          # (y,x)
area  = dy[:, None] * dx[None, :]                             # (y,x)
ocean = bathy < 0
Aoc   = area[ocean].sum()

# separar  tau = g(t) * S(y,x)
ipk = np.unravel_index(np.argmax(np.abs(tau)), tau.shape)
S   = tau[ipk[0]].copy()                                      # envolvente espacial (y,x)
g   = tau[:, ipk[1], ipk[2]] / S[ipk[1], ipk[2]]             # pulso temporal, g(pico)=1
assert np.abs(tau - g[:, None, None] * S[None, :, :]).max() < 1e-12, "tau no separable"

# --------------------------- amplitud uniforme -----------------------------
S_oc = np.where(ocean, S, 0.0)
if NORM == 'energy':
    A = np.sqrt(((S_oc**2) * area).sum() / Aoc)              # RMS de area
elif NORM == 'momentum':
    A = (np.abs(S_oc) * area).sum() / Aoc                    # media de area
else:
    raise ValueError(NORM)

# --------------------------- construir y escribir --------------------------
# uniforme en TODO el dominio (x,y); MITgcm ignora el viento en celdas secas
uni = (g[:, None, None] * A) * np.ones((NT, NY, NX))         # (t,y,x)
uni.astype(DT).tofile(OUT)

# --------------------------- verificacion ----------------------------------
chk = np.fromfile(OUT, dtype=DT).reshape(NT, NY, NX)
E_loc = ((S_oc**2) * area).sum(); E_uni = ((np.where(ocean, A, 0.0)**2) * area).sum()
I_loc = (np.abs(S_oc) * area).sum(); I_uni = (np.where(ocean, A, 0.0) * area).sum()
print(f"NORM = {NORM}")
print(f"  tau pico local        : {S.max():.5f} N/m^2  (cubre {100*ocean[S!=0].size/(NX*NY):.1f}% celdas, parche)")
print(f"  tau uniforme          : {A:.6f} N/m^2  ({S.max()/A:.1f}x mas debil)")
print(f"  ∬|tau|^2 dA  local/uni: {E_loc:.4e} / {E_uni:.4e}  (ratio {E_uni/E_loc:.4f})")
print(f"  ∬|tau|   dA  local/uni: {I_loc:.4e} / {I_uni:.4e}  (ratio {I_uni/I_loc:.4f})")
print(f"  registros no-cero g(t): {int((g>1e-9).sum())}  (pico en rec {ipk[0]})")
print(f"  escrito: {OUT}  ({chk.size} valores, max {chk.max():.6f})")
