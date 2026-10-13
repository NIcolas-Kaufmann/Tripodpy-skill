# TriPoDPy Handoff Reference

This file documents how to emit a valid `SimulationHandoff/v1` JSON
after a TriPoDPy run or analysis session.
The generic schema is defined in `handoff_schemas.md`, which belongs to the
surrounding multi-agent setup and is **not bundled with this skill**. If it is
not available, the fields below and the emitter function are sufficient.

---

## TriPoDPy-specific field conventions

| Schema field | TriPoDPy source | Notes |
|---|---|---|
| `code` | `"tripodpy"` | Always this string |
| `code_version` | `tripodpy.__version__` | e.g. `"1.0.0"` |
| `output_dir` | `sim.writer.datadir` | e.g. `"runs/my_run/data"` |
| `output_files` | `data{N:04d}.hdf5` + `frame.dmp` | List all present files |
| `diagnostics.n_snapshots` | count of `data*.hdf5` files | |
| `diagnostics.last_snap` | highest N in `data{N:04d}.hdf5` | |
| `diagnostics.t_end_code` | `f["t"][()]` from last HDF5, **converted to years** | Units block says `"time": "yr"` |
| `diagnostics.rho_max` | `f["dust/rho"]` max over cells with dust above the floor | Midplane volume density, g/cm³, per characteristic size |
| `diagnostics.rho_min` | smallest positive `dust/rho` in those cells | |
| `diagnostics.smax_max` | `f["dust/s/max"]` max over dusty cells [cm] | TriPoD-specific |
| `diagnostics.wall_clock_s` | from `run_tripodpy.py` SUCCESS line | |
| `diagnostics.plot_paths` | from `plot_tripodpy.py` SUCCESS line | comma-split; verify each path exists |
| `units.length` / `mass` / `time` | `"cm"` / `"g"` / `"yr"` | raw HDF5 time is in seconds |
| `sanity_passed` | result of the sanity checks in SKILL.md | `false` → stop, do not hand off |
| `warnings` | list of `[SANITY WARN]` strings + density note | |

---

## Density convention warning (mandatory)

TriPoDPy's `dust/rho` is the **midplane volume mass density** in g/cm³ with shape
`(Nr, 5)` — one column per characteristic size `[a0, f.dv·a1, a1, f.dv·s_max, s_max]`,
computed as `Sigma[:, [0, 1, 1, 1, 1]] / (sqrt(2π) H)`. Columns 1–4 all use Σ₁, so
**`dust/rho` must not be summed over axis 1** (that counts Σ₁ four times). There is no
mass grid. For column densities use `dust/Sigma` `(Nr, 2)` = [Σ₀, Σ₁]; for a size
distribution use `tripodpy.utils.read_data(...).dust.Sigma_recon`.

This differs from DustPy (`rho` per mass bin, `(Nr, Nm)`), FARGO3D (gas surface
density in code units) and PLUTO (volume density in code units).

Always include the note in `warnings` when handing off to `@spectral-agent` or
`@mcmc-agent`:

```json
"warnings": [
  "TriPoDPy dust/rho is midplane volume density [g/cm³] per characteristic size (Nr x 5; columns 1-4 share Sigma_1 — do not sum over sizes). Use dust/Sigma (Nr x 2) for surface density and read_data().dust.Sigma_recon for a size distribution."
]
```

If `sanity_passed` is `false`, add the failing check:
```json
"warnings": [
  "SANITY FAIL: eps_max=1.43 > 1.0 at snap 18 — results may be unreliable"
]
```

---

## Full emission pattern

```python
import hashlib, json, os, glob, subprocess, h5py
import numpy as np
import tripodpy
import tripodpy.constants as c

# Path of the run script inside this skill. ${CLAUDE_SKILL_DIR} is only substituted
# in SKILL.md, so copy the resolved path from there.
SKILL_DIR    = "<this skill's directory>"   # insert the ${CLAUDE_SKILL_DIR} path shown in SKILL.md
SKILL_SCRIPT = os.path.join(SKILL_DIR, "scripts", "run_tripodpy.py")

DENSITY_NOTE = (
    "TriPoDPy dust/rho is midplane volume density [g/cm³] per characteristic size "
    "(Nr x 5; columns 1-4 share Sigma_1 — do not sum over sizes). Use dust/Sigma "
    "(Nr x 2) for surface density and read_data().dust.Sigma_recon for a size distribution."
)

def emit_tripodpy_handoff(
    task_id: str,
    run_dir: str,
    datadir: str,
    param_file: str,
    sanity_warnings: list[str],
    wall_clock_s: float,
    plot_paths: list[str] | None = None,
) -> dict:
    """
    Build, write (run_dir/handoff.json) and return a SimulationHandoff/v1 dict.

    task_id         : snake_case label (e.g. 'tripodpy_gap_test')
    run_dir         : run directory (contains run_manifest.json from run_tripodpy.py)
    datadir         : HDF5 output directory
    param_file      : setup script / manifest used to configure the run
    sanity_warnings : list of [SANITY WARN] strings from sanity_check_tripodpy()
    wall_clock_s    : elapsed wall time from the run_tripodpy.py SUCCESS line
    plot_paths      : paths from the plot_tripodpy.py SUCCESS line (optional)
    """
    run_dir    = os.path.abspath(run_dir)
    datadir    = os.path.abspath(datadir)
    param_file = os.path.abspath(param_file)

    hdf5_files   = sorted(glob.glob(os.path.join(datadir, "data*.hdf5")))
    dump_file    = os.path.join(datadir, "frame.dmp")
    manifest     = os.path.join(run_dir, "run_manifest.json")
    output_files = hdf5_files + ([dump_file] if os.path.exists(dump_file) else [])

    if not hdf5_files:
        raise RuntimeError(f"[DATA MISSING] No HDF5 files in {datadir}")
    if not os.path.exists(manifest):
        raise RuntimeError(f"[DATA MISSING] run_manifest.json not found in {run_dir}")
    inaccessible = [p for p in output_files if not os.access(p, os.R_OK)]
    if inaccessible:
        raise RuntimeError(f"[DATA MISSING] Inaccessible files: {inaccessible}")

    last_file = hdf5_files[-1]
    last_snap = int(os.path.basename(last_file)[4:8])

    with h5py.File(last_file, "r") as f:
        t_end_s = float(f["t"][()])
        Sd      = f["dust/Sigma"][:]
        dusty   = Sd.sum(-1) > f["dust/SigmaFloor"][:].sum(-1)   # mask empty cells
        rho     = f["dust/rho"][:][dusty]
        smax    = f["dust/s/max"][:][dusty]
    rho_pos = rho[rho > 0]

    try:
        skill_version = subprocess.check_output(
            ["git", "-C", SKILL_DIR, "log", "-1", "--format=%H", "--", SKILL_SCRIPT],
            text=True, stderr=subprocess.DEVNULL,
        ).strip() or "unknown"
    except Exception:
        skill_version = "unknown"

    try:
        param_md5 = hashlib.md5(open(param_file, "rb").read()).hexdigest()
    except FileNotFoundError:
        param_md5 = "[DATA MISSING]"

    sanity_passed = not any(
        "SANITY WARN" in w and any(kw in w for kw in ["eps_max", "St(s_max)", "runaway", "s_max <", "qrec"])
        for w in sanity_warnings
    )

    warnings = list(sanity_warnings)
    if DENSITY_NOTE not in warnings:
        warnings.append(DENSITY_NOTE)

    handoff = {
        "schema":               "SimulationHandoff/v1",
        "task_id":              task_id,
        "run_dir":              run_dir,
        "run_manifest":         manifest,
        "code":                 "tripodpy",
        "code_version":         tripodpy.__version__,
        "skill_script":         SKILL_SCRIPT,
        "skill_script_version": skill_version,
        "param_file":           param_file,
        "param_file_md5":       param_md5,
        "output_dir":           datadir,
        "output_files":         [os.path.abspath(p) for p in output_files],
        "diagnostics": {
            "n_snapshots":  len(hdf5_files),
            "last_snap":    last_snap,
            "t_end_code":   t_end_s / c.year,
            "rho_field":    "dust/rho",
            "rho_units":    "g/cm³ per characteristic size (midplane volumetric, Nr x 5)",
            "rho_max":      float(rho.max()) if rho.size else 0.0,
            "rho_min":      float(rho_pos.min()) if rho_pos.size else 0.0,
            "smax_max_cm":  float(smax.max()) if smax.size else 0.0,
            "wall_clock_s": wall_clock_s,
        },
        "units": {
            "length":  "cm",
            "mass":    "g",
            "time":    "yr",
            "density": "g/cm³ per characteristic size (midplane volumetric)",
        },
        "plot_paths":    [os.path.abspath(p) for p in (plot_paths or [])],
        "sanity_passed": sanity_passed,
        "warnings":      warnings,
    }

    with open(os.path.join(run_dir, "handoff.json"), "w") as fp:
        json.dump(handoff, fp, indent=2)
    return handoff
```

---

## Validation checklist before emitting

- [ ] `sanity_passed` is `True` — if not, stop and report to user
- [ ] All paths in `output_files` exist and are readable
- [ ] `run_manifest` path exists on disk
- [ ] `t_end_code` is in **years** (not seconds)
- [ ] `diagnostics.rho_field` = `"dust/rho"` and `rho_units` explicitly set; rho values from dusty cells only
- [ ] `skill_script`, `skill_script_version`, `param_file`, `param_file_md5` populated
- [ ] `task_id` matches the `prompts/` log filename
- [ ] Density convention note is in `warnings`
- [ ] `plot_paths` populated and all paths exist on disk (run `plot_tripodpy.py` first)
- [ ] `handoff.json` written to `run_dir`
