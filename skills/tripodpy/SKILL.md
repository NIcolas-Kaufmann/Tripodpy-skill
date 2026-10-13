---
name: "tripodpy"
description: "Expert skill for setting up, running, and analysing Tripodpy simulations of dust evolution in protoplanetary disks. Use this skill whenever the user mentions Tripodpy, dust grain growth, fragmentation, radial drift, dust-to-gas ratio evolution, or wants to run/configure/analyse a dust evolution simulation and it needs to be fast. Also triggers for: customising Tripodpy physics (ice lines, planetary gaps, planetesimal formation, custom coagulation kernels), reading Tripodpy HDF5 output, diagnosing Tripodpy runs, or writing a Tripodpy setup script. Always use this skill before writing any Tripodpy-related code.\n"
---

# Tripodpy Skill

Tripodpy (Kaufmann et al. 2025) simulates the coupled radial
evolution of gas and dust in protoplanetary disks: viscous gas evolution,
dust advection/diffusion, and grain growth via the Tripod method.
Built on the `simframe` framework. additionaly it alows to track the compositional evolution of the dust.

**All internal quantities are in CGS units** — this is non-negotiable and the
most common source of errors.

---

## Quick reference

| Install | `pip install tripodpy` (requires Python 3 + Fortran compiler) |
|---|---|
| Import | `from tripodpy import Simulation` |
| Constants | `from tripodpy import constants as c` — use `c.au`, `c.M_sun`, `c.year`, etc. |
| Output | HDF5 files `data/data{N:04d}.hdf5` + dump file `data/frame.dmp` |
| Plotting | `from tripodpy import plot; plot.panel(sim)` or `plot.panel("data/")` |
| Citation | Kaufmann et al. (2025) — always remind user |

---

## Mandatory workflow

For every Tripodpy task (setup, run, analysis), follow this sequence:

1. **Confirm the task type**: setup-and-run, analyse existing output, or customise physics.
2. **For analysis**: verify the output directory exists and at least one HDF5 file is present before reading anything. Read a single file first.
3. **Never assume default parameters** — read `sim.ini` or the first HDF5 to confirm grid, mass, and physical parameters actually used.
4. **Run sanity checks** after reading data (see below).
5. **Emit `SimulationHandoff/v1`** (see `references/handoff.md`) after a successful run or analysis.

---

## Simulation object structure

```
sim.dust        — dust quantities (Sigma, St, eps, a, rho, v, p, S, …)
sim.gas         — gas quantities (Sigma, alpha, T, cs, nu, eta, rho, …)
sim.grid        — grid (r [cm], ri, m [g], Nr, Nm, A, OmegaK)
sim.star        — stellar parameters (M, R, T, L)
sim.t           — integration variable, time in seconds
sim.integrator  — numerical integrator
sim.components  — physics components (Default, other components can be added e.g. H20 etc.)
sim.writer      — output writer (datadir, overwrite, snapshots)
```

Key 2-D fields have either shape `(Nr, Nm_s)` — e.g. `sim.dust.Sigma`, or `(Nr, Nm_l)` — e.g. `sim.dust.St`, `sim.dust.a`, `sim.dust.rhos` where `Nm_s` = 2 and `Nm_l`= 5.
Key 1-D fields have shape `(Nr,)` — e.g. `sim.gas.Sigma`, `sim.gas.T`.

---

## Standard setup pattern

```python
from tripodpy import Simulation
from dustpy import constants as c
import numpy as np

sim = Simulation()

# --- stellar ---
sim.ini.star.M = 1.0 * c.M_sun      # [g]
sim.ini.star.R = 2.0 * c.R_sun      # [cm]
sim.ini.star.T = 5772.               # [K]

# --- grid ---
sim.ini.grid.Nr    = 100
sim.ini.grid.rmin  = 1.  * c.au     # [cm]
sim.ini.grid.rmax  = 1e3 * c.au     # [cm]

# --- gas ---
sim.ini.gas.alpha  = 1e-3
sim.ini.gas.Mdisk  = 0.05 * c.M_sun # [g]
sim.ini.gas.SigmaExp = -1.0         # surface density power law
sim.ini.gas.SigmaRc  = 30. * c.au  # [cm] cut-off radius

# --- dust ---
sim.ini.dust.d2gRatio  = 0.01
sim.ini.dust.vFrag     = 100.       # [cm/s] = 1 m/s
sim.ini.dust.rhoMonomer = 1.67      # [g/cm³]

sim.initialize()

# --- writer and snapshots must be set AFTER initialize() ---
sim.writer.datadir = "data"
sim.writer.overwrite = False   # set True only after explicit user confirmation
sim.t.snapshots = np.hstack([
    sim.t,
    np.geomspace(1e3, 1e5, num=21) * c.year   # [s]
])

sim.run()
```

---

## Reading output files

**Always read the first snapshot before batch-processing.**

```python
import h5py, numpy as np

def read_tripodpy_snap(datadir: str, snap: int) -> dict:
    """Read one tripodpy HDF5 snapshot. All arrays returned in CGS."""
    path = f"{datadir}/data{snap:04d}.hdf5"
    with h5py.File(path, "r") as f:
        r    = f["grid/r"][:]          # [cm]   shape (Nr,)
        m    = f["grid/m"][:]          # [g]    shape (Nm,Nm_l)
        t    = f["t"][()]              # [s]    scalar
        Sigma_gas  = f["gas/Sigma"][:] # [g/cm²] shape (Nr,)
        T_gas      = f["gas/T"][:]     # [K]     shape (Nr,)
        Sigma_dust = f["dust/Sigma"][:] # [g/cm²] shape (Nr, Nm_s)
        St         = f["dust/St"][:]   # [-]    shape (Nr, Nm_l)
        a          = f["dust/a"][:]    # [cm]   shape (Nr, Nm_l)
        eps        = f["dust/eps"][:]  # [-]    shape (Nr,Nm_s)  dust-to-gas ratio
    return dict(r=r, m=m, t=t,
                Sigma_gas=Sigma_gas, T_gas=T_gas,
                Sigma_dust=Sigma_dust, St=St, a=a, eps=eps)

# Verify single snapshot first
snap0 = read_tripodpy_snap("data", 0)
print("Nr, Nm:", snap0["Sigma_dust"].shape)
print("t_0 [yr]:", snap0["t"] / c.year)
print("Sigma_gas range [g/cm²]:", snap0["Sigma_gas"].min(), snap0["Sigma_gas"].max())
```

**Resuming from dump:**
```python
sim = Simulation()
sim.writer.datadir = "data"
sim.initialize()
sim.load("data/frame.dmp")
sim.run()
```

---

## Plotting

Use the bundled plotting script to produce publication-quality diagnostic figures
directly from HDF5 output — no interactive display required.

```bash
# All three plot types (panel, radial profiles, space-time evolution), last snapshot:
python ~/.agents/skills/tripodpy/scripts/plot_tripodpy.py \
    --run_dir data/dust/my_run --plot all --out plots/dust/my_run

# Radial profiles at first and last snapshot only:
python ~/.agents/skills/tripodpy/scripts/plot_tripodpy.py \
    --run_dir data/dust/my_run --plot radial --snaps 0 -1

# 6-panel overview at snapshot 5, save as PNG:
python ~/.agents/skills/tripodpy/scripts/plot_tripodpy.py \
    --run_dir data/dust/my_run --plot panel --snap 5 --fmt png
```

**Output on success (last stdout line):**
```
SUCCESS plot_dir=<path>  files=<comma-separated absolute paths>
```

| Plot type | File suffix | Contents |
|---|---|---|
| `panel` | `_panel.pdf` | 6-panel overview: dust density map (r–m), mass spectrum, gas/dust mass evolution, radial Σ, ε(r) |
| `radial` | `_radial.pdf` | 4-panel radial profiles: Σ_gas/Σ_dust, ε(r), a_max(r), St_max(r) at selected snapshots |
| `evolution` | `_evolution.pdf` | Space-time diagrams: a_max(r,t) and Σ_dust(r,t) over all snapshots |

**Key CLI flags:**

| Flag | Default | Description |
|---|---|---|
| `--run_dir` | required | Run directory (contains `data/` sub-dir with HDF5 files) |
| `--plot` | `all` | `panel`, `radial`, `evolution`, or `all` |
| `--snaps` | `0 -1` | Snapshot indices for radial plot (space-separated, negative = from end) |
| `--snap` | `-1` | Single snapshot index for panel plot |
| `--out` | `plots/dust/<run_name>/` | Output directory |
| `--fmt` | `pdf` | `pdf`, `png`, or `both` |

---

## Sanity checks (mandatory before reporting results)

Run after reading each snapshot. Emit `[SANITY WARN]` for out-of-range values
rather than silently continuing.

```python
def sanity_check_tripodpy(snap: dict, snap_index: int) -> list[str]:
    warnings = []
    r_au = snap["r"] / c.au

    # Gas surface density
    sig_g = snap["Sigma_gas"]
    if sig_g.max() > 1e5 or sig_g[sig_g > 0].min() < 1e-4:
        warnings.append(f"[SANITY WARN] snap {snap_index}: Sigma_gas range "
                        f"[{sig_g.min():.2e}, {sig_g.max():.2e}] g/cm² unusual")

    # Gas temperature
    T = snap["T_gas"]
    if T.max() > 3000 or T.min() < 5:
        warnings.append(f"[SANITY WARN] snap {snap_index}: T range "
                        f"[{T.min():.1f}, {T.max():.1f}] K unusual")

    # Dust surface density (integrated over mass)
    Sigma_d_tot = snap["Sigma_dust"].sum(axis=1)
    if Sigma_d_tot.max() > 1e4:
        warnings.append(f"[SANITY WARN] snap {snap_index}: total Sigma_dust "
                        f"max={Sigma_d_tot.max():.2e} g/cm² — check floor/source")

    # Stokes number
    St = snap["St"]
    if St.max() > 10:
        warnings.append(f"[SANITY WARN] snap {snap_index}: St_max={St.max():.2e} "
                        f"> 10 — drift instability regime")

    # Dust-to-gas ratio
    eps = snap["eps"]
    if eps.max() > 1.0:
        warnings.append(f"[SANITY WARN] snap {snap_index}: eps_max={eps.max():.2f} "
                        f"> 1 — streaming instability likely; check if intended")

    # Particle size
    a_cm = snap["a"]
    a_max = a_cm.max()
    if a_max > 1e3:   # > 10 m — almost certainly numerical runaway
        warnings.append(f"[SANITY WARN] snap {snap_index}: a_max={a_max:.2e} cm "
                        f"> 10 m — possible numerical runaway")

    return warnings
```

Expected physical ranges:

| Quantity | Typical range | Notes |
|---|---|---|
| `gas.Sigma` | 1 – 10⁴ g cm⁻² | Inner disk higher |
| `gas.T` | 10 – 2000 K | Midplane profile |
| `dust.Sigma` (total) | 0.01 – 100 g cm⁻² | Integrated over mass bins |
| `dust.St` | 10⁻⁶ – 1 | Warn > 1 |
| `dust.eps` | 10⁻⁴ – 0.1 bulk | Warn > 1 |
| `dust.a` (peak) | 0.1 µm – 10 cm | Fragmentation-limited |
| `dust.v.frag` | 100 – 1000 cm/s | Default 100 cm/s |

---

## Common customisations

For detailed patterns, read `references/customisation.md`.

**Quick index:**

- Change physics function (e.g. temperature profile):
  `sim.gas.T.updater = my_T_func` — call before `sim.initialize()`
- Turn off coagulation entirely: set both probabilities to zero + unset
  `v.rel`, `v.frag` updaters
- Turn off gas evolution: `sim.gas.S.hyd.updater = None`
- Custom radial grid (e.g. refinement around a gap): set `sim.grid.ri` before
  calling `sim.initialize()` — see `references/customisation.md`
- Add torque for planetary gap (Kanagawa 2017 profile): see
  `references/customisation.md § Planetary gaps`
- Add external dust/gas sources: `sim.dust.S.ext` / `sim.gas.S.ext`
- Systoles (called every timestep) / diastoles (called every output):
  use `sim.addSystemicUpdate()` — see `references/customisation.md`

---

## Output conventions for `SimulationHandoff/v1`

When emitting the handoff schema (see shared `handoff_schemas.md`):

```python
import h5py, subprocess

with h5py.File(f"{datadir}/data{last_snap:04d}.hdf5", "r") as f:
    rho_max = float(f["dust/rho"][:].max())   # [g/cm³] per mass bin, midplane
    rho_min = float(f["dust/rho"][:][f["dust/rho"][:] > 0].min())
    t_end   = float(f["t"][()])               # [s]

# Convert t_end to years for the handoff units block
handoff = {
    "schema": "SimulationHandoff/v1",
    "code": "tripodpy",
    "code_version": subprocess.check_output(
        ["pip", "show", "tripodpy"], text=True
    ),  # parse "Version: X.Y.Z" from output
    "diagnostics": {
        "n_snapshots": last_snap + 1,
        "last_snap": last_snap,
        "t_end_code": t_end / c.year,   # report in years
        "rho_max": rho_max,
        "rho_min": rho_min,
    },
    "units": {
        "length": "cm",
        "mass": "g",
        "time": "yr",   # t_end converted above
    },
    # ... run_dir, output_dir, output_files, sanity_passed, warnings
}
```

**tripodpy density convention for handoff**: `rho` in the handoff refers to
`dust.rho` — midplane mass density per mass bin in g/cm³, shape `(Nr, Nm_l)`.
This differs from FARGO3D (surface density) and PLUTO (volume density).
Always document this in `warnings` if handing off to a code that expects
surface density: `"warnings": ["tripodpy rho is midplane volumetric [g/cm³] per mass bin, not surface density"]`.

---

## Performance notes

- `Nr=100` is the standard; increasing to 200 increases runtime ~4×.
- The Smoluchowski solver uses sparse implicit matrices — coagulation cost scales as O(Nr × Nm_l²).
- On HPC: tripodpy is single-node only (no MPI). Use array jobs for parameter surveys.
- For long runs (> 10⁶ yr with full coagulation), budget ~30–60 min on a modern single core.

---

## Iron rules

- All internal quantities **must** be in CGS units — never mix units silently.
- Never overwrite existing output without explicit user confirmation (`sim.writer.overwrite = False`).
- Always read and sanity-check a single snapshot before batch-processing.
- Emit `SimulationHandoff/v1` after every successful run or analysis.
- Never invent parameter values — read from `sim.ini` or the first HDF5.

---

## References

- **All parameters with types, defaults, constraints, and interactions** → `references/parameters.md`
- Full customisation patterns → `references/customisation.md`
- Handoff schema → `references/handoff.md` + shared `handoff_schemas.md`
- tripodpy docs: https://tripodpy.readthedocs.io/en/latest/
- Simframe docs: https://simframe.rtfd.io/
