---
name: "tripodpy"
description: "Expert skill for setting up, running, and analysing TriPoDPy simulations of dust evolution in protoplanetary disks (TriPoD: parametric tri-population dust size distribution, Pfeil et al. 2024; a fast alternative to full-coagulation DustPy). Use this skill whenever the user mentions TriPoDPy / tripodpy / TriPoD, or wants to run/configure/analyse a fast 1D dust evolution simulation (grain growth, fragmentation, radial drift, dust-to-gas ratio, maximum particle size s_max, size-distribution exponent). Also triggers for: customising TriPoDPy physics (ice lines, planetary gaps, external sources/sinks, boundaries), compositional tracking with `sim.addcomponent` (gas/dust tracers, active dust species, sublimation/condensation of volatiles such as CO or H2O), reading TriPoDPy HDF5 output, diagnosing TriPoDPy runs, or writing a TriPoDPy setup script. Always use this skill before writing any TriPoDPy-related code.\n"
---

# TriPoDPy Skill

TriPoDPy (Kaufmann et al. 2025) simulates the coupled radial evolution of gas
and dust in protoplanetary disks: viscous gas evolution, dust advection/diffusion,
and grain growth/fragmentation with the **TriPoD** parametric size-distribution
model (Pfeil et al. 2024). Instead of a mass grid it evolves, per radial cell,
two dust surface densities and the maximum particle size `s_max` of a truncated
power law. It can additionally track the composition of gas and dust
(`sim.addcomponent`).

`tripodpy.Simulation` is a child of `dustpy.Simulation` (which is a
`simframe.Frame`), so much of the DustPy structure carries over — **but not all
of it**. Several DustPy recipes silently do nothing in TriPoDPy (see
"DustPy habits that break" below).

**All internal quantities are in CGS units** — the most common source of errors.

---

## Quick reference

| | |
|---|---|
| Install | `pip install tripodpy` (Python ≥ 3.11, needs a Fortran compiler for source builds; depends on `dustpy`) |
| Import | `from tripodpy import Simulation` |
| Constants | `import tripodpy.constants as c` — `c.au`, `c.year`, `c.M_sun`, `c.R_sun`, `c.M_earth`, `c.M_jup`, `c.m_p`, `c.k_B`, `c.G`, `c.pi`, `c.sigma_sb`, `c.sigma_H2` (same values as `dustpy.constants`) |
| Output | HDF5 files `<datadir>/data{N:04d}.hdf5` + dump `<datadir>/frame.dmp` (default `datadir="data"`) |
| Reading | `from tripodpy.utils import read_data` · `from tripodpy import hdf5writer, readdump` |
| Plotting | `from tripodpy import plot; plot.panel(sim)` or `plot.panel("data", it=10, ir=20)` |
| Docs | https://tripodpy.readthedocs.io/ (notebooks `docs/source/*.ipynb` in the repo) |
| Citation | Kaufmann et al. (2025) for the code **and** Pfeil et al. (2024, A&A 691, A45) for the TriPoD method — always remind the user |

---

## Mandatory workflow

For every TriPoDPy task (setup, run, analysis):

1. **Confirm the task type**: setup-and-run, analyse existing output, customise physics, or compositional tracking.
2. **For analysis**: verify the output directory exists and contains at least one `data*.hdf5`. Read a single file first.
3. **Never assume default parameters** — read `sim.ini` or the first HDF5 to confirm grid, disk mass and physical parameters actually used.
4. **Run the sanity checks** below after reading data.
5. **Emit `SimulationHandoff/v1`** (see `references/handoff.md`) after a successful run or analysis.

---

## Model in one paragraph (what the arrays mean)

At each radius the dust follows a power law n(a) ∝ a^q between `s.min` and
`s.max`, split at `s_int = sqrt(s.min * s.max)`:

- `dust.Sigma[:, 0]` = Σ₀, mass in small grains (`s.min` … `s_int`)
- `dust.Sigma[:, 1]` = Σ₁, mass in large grains (`s_int` … `s.max`)
- `dust.qrec` = q reconstructed from Σ₁/Σ₀ via `q = log(Σ₁/Σ₀) / log(s.max/s_int) − 4` (negative, typically −3 … −4; n(a) ∝ a^q)
- `dust.a` holds 5 characteristic sizes `[a0, f.dv*a1, a1, f.dv*s.max, s.max]` (`f.dv = 0.4` by default), so `dust.a[:, -1] == dust.s.max`
- Integrated state: Σ₀, Σ₁ and `s.max · Σ₁` are evolved together by one implicit solver; gas is evolved per component (the `Default` component holds all gas by default).

So dust fields have shape **`(Nr, 2)`** (surface-density-like: `Sigma`, `SigmaFloor`, `S.*`, `Fi.*` with `Nr+1`) or **`(Nr, 5)`** (size-like: `a`, `m`, `St`, `H`, `rho`, `rhos`, `fill`, `D`, `v.rad`), or **`(Nr, 5, 5)`** for `v.rel.*`. `eps`, `s.max`, `s.min`, `qrec`, `q.eff`, `q.frag`, `p.*`, `v.frag` are **`(Nr,)`**. There is **no mass grid** (`grid.m`, `grid.Nm`, `ini.grid.Nmbpd`, `ini.grid.mmax` do not exist).

---

## Simulation object structure

```
sim.dust        — Sigma (Nr,2), s.{min,max,lim}, qrec, q.{eff,frag,sweep,turb1,turb2,drfrag},
                  a/m/St/H/rho/rhos/fill/D (Nr,5), eps (Nr,), v.{frag,rad,rad_flux,driftmax,rel.*},
                  p.{frag,stick,fragtrans,driftfrag}, f.{crit,drift,dv}, delta.{rad,turb,vert},
                  S.{ext,hyd,coag,compo,smax_hyd,tot}, Fi.{adv,diff,tot}, boundary, s.boundary
sim.gas         — Sigma, T, cs, Hp, nu, alpha, mu, eta, rho, P, n, mfp, v.{rad,visc}, Fi,
                  S.{ext,hyd,tot}, torque.{Lambda,v}, boundary, SigmaFloor   (all (Nr,))
sim.grid        — r, ri, A, Nr, OmegaK  (+ hidden _Nm_short=2, _Nm_long=5)
sim.star        — M, R, T, L
sim.components  — compositional components; always contains `Default` (all gas)
sim.t           — time [s] (integration variable; exists only after initialize())
sim.integrator  — implicit dust solver + one implicit solver per component
sim.writer      — hdf5writer (exists only after initialize())
```

---

## Standard setup pattern

```python
import numpy as np
from tripodpy import Simulation
import tripodpy.constants as c

sim = Simulation()

# --- stellar ---
sim.ini.star.M = 1.0 * c.M_sun        # [g]
sim.ini.star.R = 2.0 * c.R_sun        # [cm]
sim.ini.star.T = 5772.                # [K]

# --- grid ---
sim.ini.grid.Nr   = 100
sim.ini.grid.rmin = 1.   * c.au       # [cm]
sim.ini.grid.rmax = 1000 * c.au       # [cm]
sim.ini.grid.mmin = 1e-12             # [g] -> sets s.min = (3 mmin / (4 pi rhos))^(1/3)

# --- gas ---
sim.ini.gas.alpha    = 1e-3
sim.ini.gas.Mdisk    = 0.05 * c.M_sun # [g]
sim.ini.gas.SigmaExp = -1.0
sim.ini.gas.SigmaRc  = 60. * c.au     # [cm] (default)

# --- dust ---
sim.ini.dust.d2gRatio   = 0.01
sim.ini.dust.vFrag      = 100.        # [cm/s] = 1 m/s
sim.ini.dust.rhoMonomer = 1.67        # [g/cm³]
sim.ini.dust.aIniMax    = 1e-4        # [cm] initial s.max

sim.initialize()

# --- writer and snapshots must be set AFTER initialize() ---
# (before initialize(), sim.writer and sim.t are None)
sim.writer.datadir   = "data"
sim.writer.overwrite = False          # True only after explicit user confirmation
sim.t.snapshots = np.hstack([
    sim.t,                                       # t = 0 -> data0000.hdf5 is the initial state
    np.geomspace(1e3, 1e5, num=21) * c.year,     # [s]
])

sim.run()
```

Default snapshots (if not set): 21 log-spaced times 10³–10⁵ yr (no t = 0 file).
Measured on one core: the default model (Nr = 100, to 10⁵ yr) takes ≈ 25–30 s.

---

## Reading output

**Always read the first snapshot before batch-processing.**

```python
import h5py, numpy as np
import tripodpy.constants as c

def read_tripodpy_snap(datadir: str, snap: int) -> dict:
    """Read one TriPoDPy HDF5 snapshot. All arrays in CGS."""
    with h5py.File(f"{datadir}/data{snap:04d}.hdf5", "r") as f:
        return dict(
            t          = f["t"][()],              # [s]      scalar
            r          = f["grid/r"][:],          # [cm]     (Nr,)
            Sigma_gas  = f["gas/Sigma"][:],       # [g/cm²]  (Nr,)
            T_gas      = f["gas/T"][:],           # [K]      (Nr,)
            Sigma_dust = f["dust/Sigma"][:],      # [g/cm²]  (Nr, 2)  [Sigma_0, Sigma_1]
            SigmaFloor = f["dust/SigmaFloor"][:], # [g/cm²]  (Nr, 2)
            St         = f["dust/St"][:],         # [-]      (Nr, 5)
            a          = f["dust/a"][:],          # [cm]     (Nr, 5)
            smax       = f["dust/s/max"][:],      # [cm]     (Nr,)
            smin       = f["dust/s/min"][:],      # [cm]     (Nr,)
            qrec       = f["dust/qrec"][:],       # [-]      (Nr,)
            eps        = f["dust/eps"][:],        # [-]      (Nr,)  total dust-to-gas ratio
        )

snap0 = read_tripodpy_snap("data", 0)
print("Nr:", snap0["r"].size, " t_0 [yr]:", snap0["t"] / c.year)
```

HDF5 layout mirrors the object: `dust/s/max`, `dust/q/eff`, `gas/torque/Lambda`,
`components/<name>/gas/Sigma`, `components/<name>/dust/Sigma`, `components/<name>/gas/value`, …

Convenience readers:

```python
from tripodpy import hdf5writer
w = hdf5writer(datadir="data")
smax = w.read.sequence("dust.s.max")      # (Nt, Nr) — reads only one field
data = w.read.all()                        # everything, extra leading time axis
d3   = w.read.output(3)                    # one snapshot as namespace

from tripodpy.utils import read_data      # adds reconstructed size distributions
d = read_data("data", Na=50)               # SimpleNamespace: d.t, d.dust, d.gas, d.grid
d.dust.Sigma_recon   # (Nt, Nr, Na) g/cm² per size bin on d.grid.a_recon (Nt, Na)
d.dust.M, d.gas.M    # total masses vs time [g]
d.dust.StFr, d.dust.StDr   # fragmentation / drift-limit Stokes numbers
d.dust.q             # (Nt, Nr) exponent with POSITIVE sign (= -qrec)
```

`read_data` has no `Nt` attribute; use `d.t.size`.

**Resuming from a dump** (dumps are pickles — only load your own):
```python
from tripodpy import readdump
sim = readdump("data/frame.dmp")           # full frame incl. writer & custom updaters
# APPEND to the existing snapshots — files are numbered by index in t.snapshots
sim.t.snapshots = np.concatenate([sim.t.snapshots, [2e5 * c.year]])
sim.run()
```
There is no `sim.load()`.

---

## Bundled scripts

This skill's directory is `${CLAUDE_SKILL_DIR}`; the scripts live in
`${CLAUDE_SKILL_DIR}/scripts/`. Both print one machine-readable last line
(`SUCCESS ...` or `ERROR message=...`).

### Running (`run_tripodpy.py`)

```bash
# Default model into runs/my_run/data (+ run_manifest.json):
python ${CLAUDE_SKILL_DIR}/scripts/run_tripodpy.py --run_dir runs/my_run

# Parameter overrides (unit suffixes are converted; everything else CGS):
python ${CLAUDE_SKILL_DIR}/scripts/run_tripodpy.py --run_dir runs/a1e-4 \
    --json '{"alpha": 1e-4, "vfrag_cms": 1000, "t_end_yr": 1e6}'

# Continue an existing run from data/frame.dmp to a later time:
python ${CLAUDE_SKILL_DIR}/scripts/run_tripodpy.py --run_dir runs/my_run --resume --t_end_yr 1e6

# Validate only / write a SLURM jobscript without submitting:
python ${CLAUDE_SKILL_DIR}/scripts/run_tripodpy.py --run_dir runs/my_run --dry_run
python ${CLAUDE_SKILL_DIR}/scripts/run_tripodpy.py --run_dir runs/my_run --hpc --dry_run
```
Success line: `SUCCESS run_dir=<path> n_snapshots=<N> t_end_yr=<float> wall_clock_s=<float>`.
Existing output is never overwritten unless `--overwrite` is passed (ask the user first).
See `python ${CLAUDE_SKILL_DIR}/scripts/run_tripodpy.py --help` for all flags.

### Plotting (`plot_tripodpy.py`)

```bash
# All three plot types for a run made with run_tripodpy.py (<run_dir>/data/):
python ${CLAUDE_SKILL_DIR}/scripts/plot_tripodpy.py --run_dir runs/my_run --plot all --out plots/my_run

# Radial profiles at first and last snapshot:
python ${CLAUDE_SKILL_DIR}/scripts/plot_tripodpy.py --run_dir runs/my_run --plot radial --snaps 0 -1

# Overview panel at snapshot 5, size distribution at radial index 30, PNG:
python ${CLAUDE_SKILL_DIR}/scripts/plot_tripodpy.py --run_dir runs/my_run --plot panel --snap 5 --ir 30 --fmt png
```

`--run_dir` may also be a plain data directory containing `data*.hdf5`.
Last stdout line: `SUCCESS plot_dir=<path>  files=<comma-separated absolute paths>` or `ERROR message=...`.

| Plot type | File suffix | Contents |
|---|---|---|
| `panel` | `_panel.pdf` | `tripodpy.plot.panel`: reconstructed size distribution vs r (with drift/frag limits), distribution at `ir`, gas/dust mass vs time, midplane density + q, Σ_gas/Σ_dust + ε |
| `radial` | `_radial.pdf` | Σ_gas/Σ_dust, ε, s_max, St(s_max), q_rec, Σ₀/Σ₁ at selected snapshots (empty cells masked) |
| `evolution` | `_evolution.pdf` | Space-time diagrams s_max(r,t) and Σ_dust(r,t) |

| Flag | Default | Description |
|---|---|---|
| `--run_dir` | required | Run dir (with `data/`) or data dir |
| `--plot` | `all` | `panel`, `radial`, `evolution`, `all` |
| `--snaps` | `0 -1` | Snapshot indices for radial plot (negative = from end) |
| `--snap` | `-1` | Snapshot for panel plot (negative = from end) |
| `--ir` | `20` | Radial index for the size-distribution panel |
| `--out` | `plots/<run_name>/` | Output directory |
| `--fmt` | `pdf` | `pdf`, `png`, `both` |

In your own code: `plot.panel(...)` calls `plt.show()` and returns nothing; to save,
use the Agg backend and `plt.gcf().savefig(...)`. `plot.panel` clips `it < 0` to 0 —
pass the explicit last index, not `-1`.

---

## Sanity checks (mandatory before reporting results)

Empty cells (dust at `0.1 * SigmaFloor`, gas at the `1e-100` floor — typically the
outer disk beyond the gas cut-off) have meaningless `eps`/`St` of order 10⁹⁰.
**Always mask them** before checking. Emit `[SANITY WARN]` instead of silently continuing.

```python
def sanity_check_tripodpy(snap: dict, snap_index: int) -> list[str]:
    warnings = []
    Sd_tot = snap["Sigma_dust"].sum(axis=1)
    m = Sd_tot > snap["SigmaFloor"].sum(axis=1)          # cells that contain dust
    if not m.any():
        return [f"[SANITY WARN] snap {snap_index}: no cell has dust above the floor value"]

    sig_g = snap["Sigma_gas"][m]
    if sig_g.max() > 1e5:
        warnings.append(f"[SANITY WARN] snap {snap_index}: Sigma_gas max={sig_g.max():.2e} g/cm² unusual")
    T = snap["T_gas"]
    if T.max() > 3000 or T.min() < 5:
        warnings.append(f"[SANITY WARN] snap {snap_index}: T range [{T.min():.1f}, {T.max():.1f}] K unusual")
    if Sd_tot[m].max() > 1e4:
        warnings.append(f"[SANITY WARN] snap {snap_index}: Sigma_dust max={Sd_tot[m].max():.2e} g/cm²")
    St_max = snap["St"][m, -1]                            # Stokes number of s_max
    if St_max.max() > 10:
        warnings.append(f"[SANITY WARN] snap {snap_index}: St(s_max)={St_max.max():.2e} > 10")
    eps = snap["eps"][m]
    if eps.max() > 1.0:
        warnings.append(f"[SANITY WARN] snap {snap_index}: eps_max={eps.max():.2f} > 1 — check if intended")
    smax = snap["smax"][m]
    if smax.max() > 1e3:
        warnings.append(f"[SANITY WARN] snap {snap_index}: s_max={smax.max():.2e} cm > 10 m — possible runaway")
    if np.any(smax < 1.5 * snap["smin"][m] * (1 - 1e-6)):
        warnings.append(f"[SANITY WARN] snap {snap_index}: s_max < 1.5 s_min in dusty cells")
    if not np.all(np.isfinite(snap["qrec"][m])):
        warnings.append(f"[SANITY WARN] snap {snap_index}: non-finite qrec in dusty cells")
    return warnings
```

Values seen in the default model (Nr = 100, 0–4×10⁵ yr), dusty cells only:

| Quantity | Default-model range | Notes |
|---|---|---|
| `gas.Sigma` | 0.016 – 2.2×10³ g cm⁻² | outer cells drop to floor 10⁻¹⁰⁰ |
| `gas.T` | ~7 – 220 K | passive irradiation, φ = 0.05 |
| `dust.eps` | ≤ 0.012 (≪ at the dust outer edge) | starts at `d2gRatio` |
| `dust.s.max` | 10⁻⁴ cm (initial) … 0.29 cm | fragmentation/drift limited; ≥ 1.5 `s.min` enforced |
| `dust.St[:, -1]` | 10⁻⁷ … 0.026 | warn > 1 |
| `dust.qrec` | −3.2 … −3.6 in the bulk | down to ≈ −4.7 at the dust outer edge |
| `dust.v.frag` | `ini.dust.vFrag` (100 cm/s) | spatially varying only if you set an updater |

---

## DustPy habits that break in TriPoDPy (all verified)

| DustPy habit | In TriPoDPy |
|---|---|
| `sim.ini.grid.Nmbpd`, `ini.grid.mmax`, `sim.grid.m`, `grid.Nm`, `dust.kernel`, `dust.coagulation` | deleted — raise `AttributeError` |
| `ini.dust.distExp` | exists but **has no effect** (initial q comes from `dust.q.eff`) |
| `sim.gas.S.ext.updater = f` (gas sources) | **silently ignored**: `gas.S.ext` is recomputed as the sum of component sources. Use `sim.components.Default.gas._Sigma_dot` (see customisation.md §4) |
| Freeze gas with `sim.gas.Sigma.updater = None` | **inconsistent**: total frozen, but the `Default` component keeps evolving. Use `sim.gas.nu[:] = 0; sim.gas.nu.updater = None` |
| `sim.gas.boundary.inner/outer.setcondition(...)` | **no effect** on the solver — set it on `sim.components.Default.gas.boundary` (and other active gas components) |
| Freeze dust transport via `dust.Fi.*.updater = None` | no effect (fluxes are diagnostic) — zero `dust.v.rad_flux` and `dust.D` |
| Turn off coagulation via `dust.S.coag` | not enough — zero `dust.v.rel.tot` |
| `sim.dust.S.ext.updater = f` **before** `initialize()` | `AttributeError` (field is `None`) — set it after `initialize()` |
| `sim.dust.S.ext` with **active dust components** | ignored — `dust.Sigma` is overwritten by the sum of components; use the component's `dust._Sigma_dot` |
| `sim.load(...)`, `sim.addSystole`, `sim.addDiastole`, `sim.addSystemicUpdate` | do not exist — use `tripodpy.readdump` and `updater.systole/diastole` |
| `sim.addcomponent(...)` before `initialize()` | `AttributeError` — only exists after `initialize()` |
| Dust size/mass grid analysis on `dust.Sigma` | `dust.Sigma` has only 2 populations; use `tripodpy.utils.read_data` / `sim_size_distribution` for a reconstructed size distribution |

---

## Common customisations

Detailed, tested recipes in `references/customisation.md`; compositional tracking in
`references/components.md`.

- Replace a physics function: `sim.dust.v.frag.updater = f` (after `initialize()`), then `sim.update()`
- Ice lines: `v.frag` updater depending on `sim.gas.T`
- Freeze gas: `sim.gas.nu[:] = 0.; sim.gas.nu.updater = None`
- Planetary gap: α-profile `sim.gas.alpha[:] = alpha0 / f_gap` (e.g. `dustpylib` `duffell2020`) or a torque via `sim.gas.torque.Lambda.updater`
- Dust sources/sinks: `sim.dust.S.ext.updater` (after `initialize()`); gas sources: `sim.components.Default.gas._Sigma_dot`
- Custom radial grid: set `sim.grid.ri` before `initialize()`
- Per-timestep hooks: `sim.updater.systole = f` / `sim.updater.diastole = f` (called on every `sim.update()`, i.e. every step)
- TriPoD-specific knobs: `sim.dust.f.drift` (0.8), `sim.dust.f.dv` (0.4), `sim.dust.f.crit` (0.425), `sim.dust.s.lim` (1e-4 cm), `sim.dust.q.{turb1,turb2,drfrag,sweep}` (−3.5, −3.75, −3.75, −3)
- Components (tracers, multiple dust species, sublimation/condensation): `sim.addcomponent(...)` after `initialize()` — see `references/components.md`

---

## Output conventions for `SimulationHandoff/v1`

See `references/handoff.md`. Key point: `dust/rho` is the **midplane volume density
per characteristic size** in g/cm³, shape `(Nr, 5)`, computed from
`Sigma[:, [0,1,1,1,1]]` — columns 1–4 all share Σ₁, so **summing `dust/rho` over
axis 1 double counts**. For column density use `dust/Sigma` (Nr, 2).

---

## Performance notes (measured, single core)

- Default model (Nr = 100, 10³–10⁵ yr, 21 snapshots): ≈ 27 s; extending to 4×10⁵ yr: ≈ 45 s more.
- Nr = 40–50, 10⁴ yr: ≈ 2–4 s; adding a few components adds a few seconds.
- Cost does not depend on a mass grid (there is none): only Nr, the number of components, and the timestep matter.
- No MPI; run parameter surveys as independent jobs (`run_tripodpy.py --hpc`).

---

## Iron rules

- All quantities in CGS — never mix units silently.
- Never overwrite existing output without explicit user confirmation (`sim.writer.overwrite = False`; the writer raises `RuntimeError: File ... already exists`).
- Set `writer.datadir` and `t.snapshots` after `initialize()`.
- Always read and sanity-check a single snapshot before batch-processing; mask empty cells.
- Emit `SimulationHandoff/v1` after every successful run or analysis.
- Never invent parameter values — read them from `sim.ini` or the first HDF5.
- When in doubt whether a DustPy feature works in TriPoDPy, check the table above or test it on a small grid (Nr ≈ 40, 10³ yr) before relying on it.

---

## References

- All parameters, defaults and shapes → `references/parameters.md`
- Customisation recipes → `references/customisation.md`
- Compositional tracking (`addcomponent`) → `references/components.md`
- Handoff schema → `references/handoff.md` (TriPoDPy-specific fields; the generic `SimulationHandoff/v1` schema `handoff_schemas.md` is **not bundled** with this skill — it belongs to the surrounding multi-agent setup; if it is not available, follow `references/handoff.md` alone)
- TriPoDPy docs: https://tripodpy.readthedocs.io/ · DustPy docs (shared framework): https://stammler.github.io/dustpy/ · Simframe: https://simframe.rtfd.io/
