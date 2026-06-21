#!/bin/bash
#==============================================================================
# submit_pile_linear.slurm
# Corrida de UNA SOLA ETAPA del caso LINEAR con CONDICION INICIAL apilada
# (columna de agua apilada via pSurfInitFile), SIN viento.
#   - run_expand_pile        (con bahia,  pile_init_bahia_560x352.bin)
#   - run_expand_nobay_pile  (sin bahia,  pile_init_nobahia_560x352.bin)
#
# A diferencia de submit_2stage_linear.slurm (etapa 1 forzada por viento +
# etapa 2 libre), aqui NO hay forzante: la simulacion arranca en t=0 con la
# superficie libre ya apilada (eta crece linealmente de la boca a la cabeza de
# la bahia, eta_max = 0.0089 m = set-up maximo del viento) y se deja evolucionar.
#   startTime=0 -> endTime=432000 s (~5 d), nIter0=0, deltaT=30 -> 14400 pasos.
#
# El campo de pila se genera con input/make_initial_pile.py.
#
# Requiere el ejecutable ya compilado en build/mitgcmuv (MPI, SIZE.h 5x4 = 20
# ranks); es EL MISMO binario que el caso con viento (pSurfInitFile es del core,
# no requiere recompilar). Compilar con el optfile AVX2 (nodos Haswell):
#   source /LUSTRE/HOME/esteban/MITgcm/setup_env.sh
#   cd build && ../../../tools/genmake2 -mods=../code -mpi -of=../optfile_ometeotl_avx2
#   make depend && make -j 8
#
# Envio:  sbatch submit_pile_linear.slurm
#==============================================================================
#SBATCH --job-name=BayIW_pile
#SBATCH --partition=id
#SBATCH --nodes=1
#SBATCH --ntasks=20                # = nPx*nPy de SIZE.h (5x4)
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --output=/LUSTRE/HOME/esteban/MITgcm/verification/BayIW_Rectan_linear/logs/%x_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=estebanci@atmosfera.unam.mx

set -eo pipefail   # NB: sin '-u' (nounset): los scripts de entorno de Intel
                   # referencian variables sin asignar y abortarian el job.

EXP=/LUSTRE/HOME/esteban/MITgcm/verification/BayIW_Rectan_linear
NP=${SLURM_NTASKS:-20}
mkdir -p "$EXP/logs"

# --- Entorno Intel (ifort + Intel MPI 2021.5 + NetCDF 4.4.1.1) --------------
set +e
source /LUSTRE/HOME/esteban/MITgcm/setup_env.sh
set -e
ulimit -n 6000 || true

EXE=$EXP/build/mitgcmuv
if [[ ! -x "$EXE" ]]; then
    echo "ERROR: no existe el ejecutable $EXE"
    echo "       Compila primero el modelo con MPI (ver cabecera de este script)."
    exit 1
fi

# ---------------------------------------------------------------------------
# Corre el caso de pila en un run dir:  run_pile <run_dir>
# ---------------------------------------------------------------------------
run_pile () {
    local dir=$1
    echo "-------------------------------------------------------------------"
    echo ">>> $(basename "$dir") — pile   inicio: $(date)"
    cd "$dir"

    if grep -q "NORMAL END" output_pile.txt 2>/dev/null; then
        echo "  [skip] ya completa (NORMAL END en output_pile.txt)"
        echo "<<< $(basename "$dir") — pile   (saltada)"
        return 0
    fi

    cp -f "$EXE" .
    rm -f STDOUT.* STDERR.* output_pile.txt 2>/dev/null || true
    rm -rf mnc_* OUT_pile 2>/dev/null || true

    time mpirun -n "$NP" ./mitgcmuv > output_pile.txt 2>&1

    if ! grep -q "NORMAL END" output_pile.txt; then
        echo "ERROR: $(basename "$dir") NO termino con NORMAL END"
        echo "----- ultimas lineas de output_pile.txt -----"
        tail -40 output_pile.txt
        exit 2
    fi

    # diagnostico: eta_max en t=0 (~0.0089 m impuesto) y theta lineal
    local e0 emax tmin tmax
    e0=$(grep "dynstat_eta_max"   STDOUT.0000 2>/dev/null | head -1 | awk '{print $NF}')
    emax=$(grep "dynstat_eta_max" STDOUT.0000 2>/dev/null | tail -1 | awk '{print $NF}')
    tmax=$(grep "dynstat_theta_max" STDOUT.0000 2>/dev/null | tail -1 | awk '{print $NF}')
    tmin=$(grep "dynstat_theta_min" STDOUT.0000 2>/dev/null | tail -1 | awk '{print $NF}')
    echo "    eta_max t0/final   = ${e0:-NA} / ${emax:-NA}  (impuesto ~0.0089 m, debe decaer)"
    echo "    theta_min/max final = ${tmin:-NA} / ${tmax:-NA}  (linear: 33.6 sup, ~3.5 fondo)"

    mkdir -p OUT_pile
    for m in mnc_*; do [ -e "$m" ] && mv "$m" OUT_pile/; done
    echo "<<< $(basename "$dir") — pile   fin: $(date)"
}

for DIR in "$EXP/run_expand_pile" "$EXP/run_expand_nobay_pile"; do
    echo "==================================================================="
    echo "### CASO: $(basename "$DIR") ###"
    run_pile "$DIR"
    echo "### $(basename "$DIR"): COMPLETA ###"
done

echo "==================================================================="
echo "TODAS LAS SIMULACIONES DE PILA (con y sin bahia) FINALIZADAS: $(date)"
