# Mensaje para asesora — Fronteras abiertas BayIW_Rectan_linear

Hola profe! 👋 Le cuento los avances con las fronteras abiertas del modelo de la bahía. Ya quedaron "realmente abiertas" para que la onda barotrópica salga sin rebotar 🌊

---

**El problema** 🤔
En la malla uniforme de 5 km, la onda barotrópica viaja a c = √(gH) ≈ 76.7 m/s (con H≈600 m). Eso da una CFL = c·dt/dx ≈ 0.46, justo en el límite de la radiación de Orlanski (necesita CMAX < 0.5). O sea, la frontera reflejaba la onda de vuelta hacia la bahía 😣

**Lo que hice — 2 cosas que se complementan** ✅

*1) Malla estirada en los márgenes* (`make_stretched_grid.py`)
Estiré geométricamente las 20 celdas de cada borde abierto (sur, este, oeste), con factor 1.05 por celda. Así la celda del borde crece de 5000 m → 13266 m, y la CFL baja de **0.46 → 0.17** 🎯 Con eso Orlanski radia con holgura. El interior y la costa norte quedan intactos.

*2) Afiné los parámetros OBCS* (`data.obcs`)
- Orlanski: cvelTimeScale 1000→250, CMAX 1.0→0.45 (1.0 era inestable)
- Esponja: grosor 10→30 celdas, relajación interior 600→1200 s, borde 90→60 s

La esponja absorbe lo poco que quede después de la radiación 👍

**Por qué es seguro** 🔒
No toca la física de la bahía: T y S son uniformes, el margen sur es plano y profundo, los bordes E/O están lejos de la bahía, y el viento es nulo en los márgenes. Se conserva el tamaño del dominio (Nx=560, Ny=352).

---

Ya verifiqué los números con los binarios reales de la malla y todo cuadra ✨ Lo apliqué tanto en el caso *con bahía* (`run_expand`) como en el *sin bahía* (`run_expand_nobay`).

Cualquier duda me dice, profe! 🙌
