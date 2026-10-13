# TriPoDPy Handoff Reference

This file documents how to emit a valid `SimulationHandoff/v1` JSON
after a TriPoDPy run or analysis session.
The full schema is defined in the shared `handoff_schemas.md`.

---

## DustPy-specific field conventions

| Schema field | DustPy source | Notes |
|---|---|---|
| `code` | `"tripodpy"` | Always this string |
| `code_version` | `pip show dustpy` | Parse "Version: X.Y.Z" |
| `output_dir` | `sim.writer.datadir` | e.g. `"runs/my_run/data"` |
| `output_files` | `data{N:04d}.hdf5` + `frame.dmp` | List all present files |
| `diagnostics.n_snapshots` | count of `data*.hdf5` files | |
| `diagnostics.last_snap` | highest N in `data{N:04d}.hdf5` | |
| `diagnostics.t_end_code` | `f["t"][()]` from last HDF5, **converted to years** | Units block says `"time": "yr"` |
| `diagnostics.rho_max` | `f["dust/rho"][:].max()` | Midplane volumetric density, g/cm³ per mass bin |
| `diagnostics.rho_min` | `f["dust/rho"][:][f["dust/rho"][:]>0].min()` | Smallest non-zero value |
| `diagnostics.wall_clock_s` | recorded from `run_dustpy.py` output | |
| `diagnostics.plot_paths` | output of `plot_dustpy.py` (SUCCESS line) | comma-split; verify each path exists |
| `units.length` | `"cm"` | Always CGS |
| `units.mass` | `"g"` | Always CGS |
| `units.time` | `"yr"` | t_end_code converted; raw HDF5 stores seconds |
| `sanity_passed` | result of sanity checks in SKILL.md | `false` → stop, do not hand off |
| `warnings` | list of `[SANITY WARN]` strings | Empty list if all clear |

---

## Density convention warning (mandatory)

TriPoDPy's `dust.rho` is the **midplane volumetric mass density per mass bin**
in g/cm³, shape `(Nr, Nm_l)`. This differs from:
- FARGO3D: `rho` is gas surface density in code units → g/cm² after conversion
- PLUTO: `rho` is volume density but in code units (not per mass bin)

Always include this in the `warnings` list when handing off to
`@spectral-agent` or `@mcmc-agent`:

```json
"warnings": [
  "DustPy rho field is midplane volumetric density [g/cm³] per mass bin (Nr x Nm), not surface density. Use dust/Sigma [g/cm²] for surface density."
]
```

If `sanity_passed` is `false`, add the specific failing check:
```json
"warnings": [
  "SANITY FAIL: dust.eps_max=1.43 > 1.0 at snap 18 — streaming instability regime, results may be unreliable"
]
```

---

## Full emission pattern

```python
import hashlib, json, os, glob, subprocess, h5py
import dustpy.constants as c
import numpy as np

SKILL_SCRIPT = os.path.expanduser(
    "~/.agents/skills/dustpy/scripts/run_dustpy.py"
)

def emit_dustpy_handoff(
    task_id: str,
    run_dir: str,
    datadir: str,
    param_file: str,
    sanity_warnings: list[str],
    wall_clock_s: float,
    plot_paths: list[str] | None = None,
) -> dict:
    """
    Build and return a SimulationHandoff/v1 dict for a DustPy run.

    Parameters
    ----------
    task_id         : snake_case label matching prompts/ log (e.g. 'dustpy_skill_test')
    run_dir         : absolute path to the run directory (contains run_manifest.json)
    datadir         : absolute path to the HDF5 output directory
    param_file      : absolute path to the setup script / .par file used
    sanity_warnings : list of [SANITY WARN] strings from sanity_check_dustpy()
    wall_clock_s    : elapsed wall time in seconds from run_dustpy.py output
    plot_paths      : list of absolute paths from plot_dustpy.py SUCCESS line (optional)
    """
    run_dir    = os.path.abspath(run_dir)
    datadir    = os.path.abspath(datadir)
    param_file = os.path.abspath(param_file)

    # Collect output files
    hdf5_files   = sorted(glob.glob(os.path.join(datadir, "data*.hdf5")))
    dump_file    = os.path.join(datadir, "frame.dmp")
    manifest     = os.path.join(run_dir, "run_manifest.json")
    output_files = hdf5_files + ([dump_file] if os.path.exists(dump_file) else [])

    if not hdf5_files:
        raise RuntimeError(f"[DATA MISSING] No HDF5 files in {datadir}")
    if not os.path.exists(manifest):
        raise RuntimeError(f"[DATA MISSING] run_manifest.json not found in {run_dir}")

    # Verify all output files accessible
    inaccessible = [p for p in output_files
                    if not os.path.exists(p) or not os.access(p, os.R_OK)]
    if inaccessible:
        raise RuntimeError(f"[DATA MISSING] Inaccessible files: {inaccessible}")

    last_snap = len(hdf5_files) - 1
    last_file = hdf5_files[-1]

    with h5py.File(last_file, "r") as f:
        t_end_s = float(f["t"][()])
        rho     = f["dust/rho"][:]
        rho_pos = rho[rho > 0]
        rho_max = float(rho.max())
        rho_min = float(rho_pos.min()) if rho_pos.size > 0 else 0.0

    # Code version
    try:
        pip_out = subprocess.check_output(
            ["pip", "show", "dustpy"], text=True, stderr=subprocess.DEVNULL
        )
        version_line = next(l for l in pip_out.splitlines()
                            if l.startswith("Version:"))
        code_version = version_line.split(":")[1].strip()
    except Exception:
        code_version = "unknown"

    # Skill script version (git hash of the script file)
    try:
        skill_version = subprocess.check_output(
            ["git", "log", "-1", "--format=%H", "--", SKILL_SCRIPT],
            text=True, stderr=subprocess.DEVNULL,
        ).strip() or "unknown"
    except Exception:
        skill_version = "unknown"

    # param_file MD5
    try:
        param_md5 = hashlib.md5(
            open(param_file, "rb").read()
        ).hexdigest()
    except FileNotFoundError:
        param_md5 = "[DATA MISSING]"

    sanity_passed = not any(
        "SANITY WARN" in w and any(kw in w for kw in ["eps_max", "St_max", "runaway"])
        for w in sanity_warnings
    )

    # Always add the density convention advisory
    warnings = list(sanity_warnings)
    density_note = (
        "DustPy rho field is midplane volumetric density [g/cm³] "
        "per mass bin (Nr x Nm), not surface density. "
        "Use dust/Sigma [g/cm²] for column density."
    )
    if density_note not in warnings:
        warnings.append(density_note)

    handoff = {
        "schema":               "SimulationHandoff/v1",
        "task_id":              task_id,
        "run_dir":              run_dir,
        "run_manifest":         manifest,
        "code":                 "DustPy",
        "code_version":         code_version,
        "skill_script":         SKILL_SCRIPT,
        "skill_script_version": skill_version,
        "param_file":           param_file,
        "param_file_md5":       param_md5,
        "output_dir":           datadir,
        "output_files":         [os.path.abspath(p) for p in output_files],
        "diagnostics": {
            "n_snapshots":  last_snap + 1,
            "last_snap":    last_snap,
            "t_end_code":   t_end_s / c.year,
            "rho_field":    "dust/rho",
            "rho_units":    "g/cm³ per mass bin (midplane volumetric)",
            "rho_max":      rho_max,
            "rho_min":      rho_min,
            "wall_clock_s": wall_clock_s,
        },
        "units": {
            "length":  "cm",
            "mass":    "g",
            "time":    "yr",
            "density": "g/cm³ per mass bin (midplane volumetric)",
        },
        "plot_paths":    [os.path.abspath(p) for p in (plot_paths or [])],
        "sanity_passed": sanity_passed,
        "warnings":      warnings,
    }

    # Write to disk alongside the run
    out_path = os.path.join(run_dir, "handoff.json")
    with open(out_path, "w") as fp:
        json.dump(handoff, fp, indent=2)

    return handoff
```

---

## Validation checklist before emitting

- [ ] `sanity_passed` is `True` — if not, stop and report to user
- [ ] All paths in `output_files` exist and are readable
- [ ] `run_manifest` path exists on disk
- [ ] `t_end_code` is in **years** (not seconds)
- [ ] `diagnostics.rho_field` = `"dust/rho"` and `rho_units` explicitly set
- [ ] `rho_max` / `rho_min` are from `dust/rho` (volumetric, not surface density)
- [ ] `skill_script`, `skill_script_version`, `param_file`, `param_file_md5` populated
- [ ] `task_id` matches the `prompts/` log filename
- [ ] Density convention warning is in `warnings`
- [ ] `plot_paths` populated and all paths exist on disk (run `plot_dustpy.py` first)
- [ ] `handoff.json` written to `run_dir`
