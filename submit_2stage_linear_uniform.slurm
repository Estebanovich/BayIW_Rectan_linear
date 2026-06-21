#!/bin/bash
#==============================================================================
# submit_2stage_linear.slurm
# Corrida completa en 2 etapas del caso LINEAR (estratificacion lineal)
#   - run_expand        (con bahia)
#   - run_expand_nobay  (control sin bahia)
#
# Etapa 1 (forzada):  startTime=0      -> endTime=216000 s (~2.5 d), nIter0=0
#                     viento periodico ON, lee T linear (33.6 sup -> 3.5 fondo),
#                     escribe el checkpoint permanente pickup.ckptA
# Etapa 2 (libre):    startTime=216000 -> endTime=432000 s (~5 d),  nIter0=7200
#                     reinicia desde pickup.ckptA, viento OFF
#
# Requiere el ejecutable ya compilado en build/mitgcmuv (MPI, SIZE.h 5x4 = 20 ranks).
# IMPORTANTE: compilar con el optfile AVX2 (nodos de computo Haswell, sin AVX-512).
#   NO usar $MITGCM_OF (-xHost) si se compila en el login: genera binario AVX-512
#   que aborta en los nodos con "support Intel(R) AVX512F ... instructions".
# Compilacion previa:
#   source /LUSTRE/HOME/esteban/MITgcm/setup_env.sh
#   cd build && ../../../tools/genmake2 -mods=../code -mpi -of=../optfile_ometeotl_avx2
#   make depend && make -j 8
#
# Envio:  sbatch submit_2stage_linear.slurm
#==============================================================================
#SBATCH --job-name=BayIW_linear_unif
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
NP=${SLURM_NTASKS:-4}

# --- Entorno Intel (ifort + Intel MPI 2021.5 + NetCDF 4.4.1.1) --------------
# vars.sh de Intel MPI no es 'errexit/nounset clean' -> sourcearlo relajado
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
# Ejecuta una etapa en un run dir:  run_stage <run_dir> <tag> <archivo_data>
# ---------------------------------------------------------------------------
run_stage () {
    local dir=$1 tag=$2 dataf=$3
    echo "-------------------------------------------------------------------"
    echo ">>> $(basename "$dir") — $tag   inicio: $(date)"
    cd "$dir"
    cp -f "$EXE" .
    cp -f "$dataf" data
    rm -f STDOUT.* STDERR.* "output_${tag}.txt" 2>/dev/null || true

    time mpirun -n "$NP" ./mitgcmuv > "output_${tag}.txt" 2>&1

    if ! grep -q "NORMAL END" "output_${tag}.txt"; then
        echo "ERROR: $(basename "$dir") $tag NO termino con NORMAL END"
        echo "----- ultimas lineas de output_${tag}.txt -----"
        tail -40 "output_${tag}.txt"
        exit 2
    fi
    # chequeo linear: theta va de 33.6 (sup) a ~3.5 (fondo)
    local tmin tmax
    tmax=$(grep "dynstat_theta_max" STDOUT.0000 2>/dev/null | tail -1 | awk '{print $NF}')
    tmin=$(grep "dynstat_theta_min" STDOUT.0000 2>/dev/null | tail -1 | awk '{print $NF}')
    echo "    theta_min/max final = ${tmin:-NA} / ${tmax:-NA}  (linear: 33.6 sup, ~3.5 fondo)"
    echo "<<< $(basename "$dir") — $tag   fin: $(date)"
}

# iteracion final de la etapa 1 (= endTime/deltaT = 216000/30) -> pickup permanente
STAGE1_PICKUP="pickup.0000007200.data"

# ---------------------------------------------------------------------------
# Reanudable: salta una etapa si ya termino (NORMAL END) en una corrida previa.
for DIR in "$EXP/run_expand_uniform" "$EXP/run_expand_nobay_uniform"; do
    echo "==================================================================="
    echo "### CASO: $(basename "$DIR") ###"
    cd "$DIR"

    # ---- ETAPA 1: forzada, escribe pickup.0000007200 ---------------------
    if [ -f "$STAGE1_PICKUP" ] && grep -q "NORMAL END" output_stage1.txt 2>/dev/null; then
        echo "  [skip] etapa 1 ya completa ($STAGE1_PICKUP presente)"
    else
        # arranque limpio de la etapa 1
        rm -f pickup* PHref* RhoRef.* 2>/dev/null || true
        rm -rf mnc_* OUT_stage1 2>/dev/null || true
        run_stage "$DIR" stage1 "$DIR/data.stage1"
        if [ ! -f "$STAGE1_PICKUP" ]; then
            echo "ERROR: la etapa 1 en $(basename "$DIR") no produjo $STAGE1_PICKUP"
            exit 3
        fi
    fi

    # archivar salida NetCDF de la etapa 1 (idempotente: tambien si se salto arriba)
    if ls mnc_* >/dev/null 2>&1; then
        mkdir -p OUT_stage1
        for m in mnc_*; do mv "$m" OUT_stage1/; done
    fi

    # ---- ETAPA 2: libre, restart desde pickup.0000007200, viento OFF -----
    if grep -q "NORMAL END" output_stage2.txt 2>/dev/null; then
        echo "  [skip] etapa 2 ya completa"
    else
        run_stage "$DIR" stage2 "$DIR/data.stage2"
        mkdir -p OUT_stage2
        for m in mnc_*; do [ -e "$m" ] && mv "$m" OUT_stage2/; done
    fi

    echo "### $(basename "$DIR"): 2 etapas COMPLETAS ###"
done

echo "==================================================================="
echo "TODAS LAS SIMULACIONES (con y sin bahia, 2 etapas) FINALIZADAS: $(date)"
