# BayIW_Rectan_linear — Internal Waves in a Rectangular Bay (Linear Stratification)

**Master's Degree Thesis Project**

This project uses the [MITgcm](https://mitgcm.org/) ocean general circulation model to simulate the generation and propagation of **internal waves** in a rectangular bay using a **linear stratification** profile. It is one of three stratification cases studied in the thesis:

| Case | Repository | Stratification |
|---|---|---|
| Realistic | `BayIW_Rectan/` | February climatological T/S profiles |
| **Linear** | **`BayIW_Rectan_linear/`** | Linearly varying density profile |
| Two-layer | `BayIW_Rectan_2layer/` | Idealized two-layer temperature profile |

The simulations are run with and without the bay geometry to isolate the bay's influence on internal wave dynamics.

---

## Project Structure

```
BayIW_Rectan_linear/
├── build/                        # MITgcm compilation directory
├── code/                         # Model configuration and size files
│   ├── SIZE.h                    # Grid decomposition parameters
│   ├── packages.conf             # Enabled MITgcm packages
│   ├── OBCS_OPTIONS.h            # Open boundary condition options
│   └── DIAGNOSTICS_SIZE.h        # Diagnostics buffer sizes
├── input/                        # Forcing, bathymetry, and preprocessing scripts
│   ├── *.bin                     # Binary input files (bathymetry, grid spacing, T/S, wind)
│   ├── profiles_2layer/          # Two-layer profile binaries (for reference)
│   ├── profiles_gradRho/         # Gradient density profile files
│   ├── STATE/                    # Reference state generation scripts
│   ├── *.ipynb                   # Jupyter notebooks for input generation and analysis
│   └── *.py                      # Python utility scripts
├── run_expand/                   # Run directory — with bay geometry
│   ├── data                      # Main model parameter file
│   ├── data.diagnostics          # Diagnostics configuration
│   ├── data.obcs                 # Open boundary condition parameters
│   ├── data.pkg                  # Package switches
│   ├── data.mnc                  # NetCDF output settings
│   └── mnc_*/                    # NetCDF output directories (one per MPI rank)
├── run_expand_nobay/             # Run directory — without bay (control case)
├── compile_and_run_expand.sh     # Script to compile and run both simulations
├── run_expands.sh                # Run-only script (no recompilation)
├── Merge_MPI_STDOUT.sh           # Merge STDOUT files from all MPI ranks
└── compress_nc.sh                # Compress NetCDF output files
```

---

## Model Configuration

| Parameter | Value |
|---|---|
| Model | MITgcm |
| Grid type | Cartesian |
| Horizontal resolution | Variable (expanded grid: 560 × 352 points) |
| Vertical levels | 50 (stretched: 1 m near surface to ~46 m at depth) |
| Time step | 30 s |
| Equation of state | Linear |
| Free surface | Implicit |
| Lateral viscosity | 100 m² s⁻¹ |
| Vertical viscosity | 1 × 10⁻⁵ m² s⁻¹ |
| Coriolis parameter (f₀) | 6.97 × 10⁻⁵ s⁻¹ |

### Packages used

| Package | Purpose |
|---|---|
| OBCS | Open boundary conditions (Orlanski radiation + sponge layers) |
| Diagnostics | Output of density anomaly fields (`RHOAnoma`) |
| MNC | NetCDF output |

### Boundary conditions

Open boundaries are applied on the south, west, and east edges using **Orlanski radiation** conditions. A sponge layer of 10 grid cells damps spurious reflections near the boundaries. Barotropic velocity balance is applied to maintain mass conservation.

### Forcing

- **Wind**: Periodic meridional wind forcing applied in stage 1 only (period = 1200 s, cycle = 216 000 s)(`make_wind_forcing_local.ipynb`)
- **Initial temperature**: Linear stratification profile (`linear_temp_50zlev_560x352.bin`)
- **Initial salinity**: Constant (`linear_salt_50zlev_560x352.bin`)

---

## Bay Geometry

The domain features a **rectangular bay** carved into a coastal shelf. The bay geometry is defined in `input/bahia_rectan_impar_func.ipynb` with the following dimensions:

| Parameter | Value |
|---|---|
| Bay width | ~60 km (from −L/4 to L/4, where L = 119 km) |
| Bay length | ~119 km |
| Bay depth | 164 m (constant, flat bottom) |
| Open ocean depth | ~1000 m |

The domain uses a **variable horizontal grid spacing**: uniform 200 m resolution in the central region containing the bay, expanding outward to reduce the domain size while minimizing boundary reflections. A previous attempt used a parabolic (semi-circular) bay profile, which produced numerical instabilities; the rectangular geometry was adopted for all stratification scenarios.

---

## Linear Stratification

The stratification is defined by a **linearly varying density profile** from surface to bottom. Both temperature and salinity may vary linearly with depth, producing a constant buoyancy frequency N² throughout the water column. This is the classic idealization used in analytical internal wave theory and allows direct comparison of simulated wave characteristics (phase speed, dispersion relation) against known theoretical solutions.

---

## Simulation Cases

| Run directory | Bay geometry | Description |
|---|---|---|
| `run_expand/` | With bay | Primary case with rectangular bay bathymetry |
| `run_expand_nobay/` | Without bay | Control case — open coastal domain |

Both cases share the same model executable and physical parameters, and follow the same two-stage run strategy described below.

### Two-stage run strategy

Each simulation is run in two consecutive stages:

| Stage | `startTime` | `endTime` | `nIter0` | Wind forcing | Purpose |
|---|---|---|---|---|---|
| 1 — Forced | 0 s | 216 000 s (~2.5 days) | 0 | ON (periodic, 1200 s period) | Wind-driven internal wave generation |
| 2 — Free | 216 000 s | 432 000 s (~5 days total) | 7200 | OFF | Free propagation and decay after forcing stops |

Stage 1 starts from rest with the linear T/S initial conditions and applies periodic wind forcing. It writes a checkpoint (`ckptA`) at the end. Stage 2 restarts from that checkpoint (`pickupSuff='ckptA'`) with the wind forcing commented out, allowing the internal wave field to propagate and decay freely.

The `data` files committed to the repository correspond to **stage 2**. To run stage 1, update `&PARM03` to:
```
startTime=0., endTime=216000., nIter0=0,
periodicExternalForcing=.TRUE., externForcingPeriod=1200., externForcingCycle=216000.,
```
and uncomment the `meridWindFile` and `hydrogThetaFile` lines in `&PARM05`.

---

## Differences from the Realistic Stratification Case

| Parameter | Realistic (`BayIW_Rectan`) | Linear (`BayIW_Rectan_linear`) |
|---|---|---|
| Temperature profile | February climatology | Linearly varying with depth |
| Salinity profile | February climatology (33.6 → 34.7 PSU) | Constant or linearly varying |
| Buoyancy frequency N² | Variable with depth | Constant (uniform N²) |
| T/S input files | `feb_temp/salt_50zlev_560x352.bin` | `linear_temp/salt_50zlev_560x352.bin` |
| Run strategy | Single stage (fresh run, wind-forced) | Two stages: forced (0–216 000 s) then free (216 000–432 000 s) |

---

## Input File Generation

The `input/` directory contains Jupyter notebooks used to prepare all binary input files:

| Notebook | Purpose |
|---|---|
| `bahia_rectan_impar_func.ipynb` | Rectangular bay bathymetry generation |
| `make_binary_*.ipynb` | Convert bathymetry to MITgcm binary format |
| `make_T_S_binary_files_linear.ipynb` | Generate linear stratification T/S initial conditions |
| `make_wind_forcing*.ipynb` | Generate periodic wind forcing fields |
| `check_output*.ipynb` | Post-processing and visualization of model output |
| `check_output_scenarios_func.ipynb` | Cross-scenario comparison with other stratification cases |

---

## How to Compile and Run

### Requirements

- MITgcm source tree (expected at `../../../` relative to this directory)
- Fortran compiler (`gfortran`) — build options configured for `darwin_arm64_gfortran`
- MPI (optional, for parallel runs)
- Python 3 with `numpy`, `scipy`, `matplotlib`, `xarray`, `netCDF4` (for preprocessing and analysis)
- Jupyter Notebook

### Compilation and execution

```bash
bash compile_and_run_expand.sh
```

The script will interactively ask whether to:
1. Clean the `build/` directory before compiling
2. Enable MPI compilation
3. Clean run directories before execution
4. Run with MPI (and how many cores)

Both `run_expand/` (with bay) and `run_expand_nobay/` (without bay) are run sequentially.

### Run only (no recompilation)

```bash
bash run_expands.sh
```

### Run a single case

```bash
cd run_expand
bash run_expand.sh
```

Or manually:

```bash
cp ../build/mitgcmuv .
mpirun -np 4 ./mitgcmuv > output.txt
```

---

## Output

Model output is written as **NetCDF** files into `mnc_*/` subdirectories (one per MPI rank). The primary diagnostic field is:

- `diag_rho` — density anomaly (`RHOAnoma`), output every 900 s (15 min)

### Post-processing utilities

```bash
bash Merge_MPI_STDOUT.sh    # Merge STDOUT logs from all MPI ranks into one file
bash compress_nc.sh         # Compress NetCDF output to reduce disk usage
```

---

## Scientific Context

Internal waves are gravity waves that propagate along density interfaces in the ocean interior. The linear stratification case produces a **constant buoyancy frequency N²**, which simplifies the theoretical treatment of internal waves: phase speeds, vertical mode structures, and dispersion relations all have clean analytical forms under this assumption. This makes the linear case ideal for:

1. Validating the numerical model against theoretical predictions
2. Isolating the effect of bay geometry from stratification complexity
3. Providing a reference against which the two-layer and realistic stratification results can be compared

---

## Author

Esteban Cruz Isidro
