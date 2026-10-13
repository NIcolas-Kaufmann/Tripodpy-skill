#!/usr/bin/env python3
"""
run_tripodpy.py — TriPoDPy launch script for the simulation agent.

Usage (local):
    python run_tripodpy.py --run_dir runs/my_run [options]
    python run_tripodpy.py --run_dir runs/my_run --json '{"alpha": 1e-3, ...}'

Usage (HPC / SLURM dry-run):
    python run_tripodpy.py --run_dir runs/my_run --hpc --dry_run

Output (stdout, last line):
    SUCCESS run_dir=<path> n_snapshots=<N> t_end_yr=<float> wall_clock_s=<float>
    ERROR   message=<description>

All physical parameters must be in CGS unless the argument name ends in
a unit suffix (_msun, _rsun, _au, _yr, _cms, _cm, _g) — those are converted internally.

TriPoDPy has no mass grid: there are always 2 surface-density populations
and 5 characteristic sizes. ``mmin_g`` only sets the minimum particle size
s_min = (3 mmin / (4 pi rho_s))^(1/3).
"""

import argparse
import json
import os
import sys
import time
import traceback

import numpy as np

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------


def parse_args():
    p = argparse.ArgumentParser(description="Launch a TriPoDPy dust-evolution simulation.")

    # --- run control ---
    p.add_argument(
        "--run_dir",
        required=True,
        help="Output directory for this run (created if absent).",
    )
    p.add_argument("--json", default=None, help="JSON string of parameter overrides (see below).")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing data directory.")
    p.add_argument(
        "--resume",
        action="store_true",
        help="Resume from existing frame.dmp dump file.",
    )
    p.add_argument(
        "--dry_run",
        action="store_true",
        help="Validate parameters and write setup script only; " "do not run the simulation.",
    )

    # --- HPC / SLURM ---
    p.add_argument(
        "--hpc",
        action="store_true",
        help="Write a SLURM jobscript instead of running locally.",
    )
    p.add_argument("--partition", default="serial", help="SLURM partition name.")
    p.add_argument("--walltime", default="04:00:00", help="SLURM wall time (HH:MM:SS).")
    p.add_argument("--memory_gb", type=float, default=8.0, help="Memory per node in GB for SLURM.")

    # --- stellar ---
    p.add_argument("--star_mass_msun", type=float, default=1.0)
    p.add_argument("--star_radius_rsun", type=float, default=2.0)
    p.add_argument("--star_teff", type=float, default=5772.0)

    # --- grid ---
    p.add_argument("--nr", type=int, default=100)
    p.add_argument("--rmin_au", type=float, default=1.0)
    p.add_argument("--rmax_au", type=float, default=1000.0)
    p.add_argument("--mmin_g", type=float, default=1e-12, help="Minimum particle mass [g]; sets s_min.")

    # --- gas ---
    p.add_argument("--alpha", type=float, default=1e-3)
    p.add_argument("--mdisk_msun", type=float, default=0.05)
    p.add_argument("--sigma_exp", type=float, default=-1.0)
    p.add_argument("--sigma_rc_au", type=float, default=60.0)

    # --- dust ---
    p.add_argument("--d2g", type=float, default=0.01)
    p.add_argument(
        "--vfrag_cms",
        type=float,
        default=100.0,
        help="Fragmentation velocity in cm/s (default 100 = 1 m/s).",
    )
    p.add_argument("--rho_monomer", type=float, default=1.67, help="Monomer bulk density [g/cm³].")
    p.add_argument(
        "--a_ini_max_cm",
        type=float,
        default=1e-4,
        help="Initial maximum particle size s_max [cm] (default 1 µm).",
    )
    p.add_argument(
        "--allow_drifting",
        action="store_true",
        help="Set ini.dust.allowDriftingParticles=True (default False).",
    )

    # --- time ---
    p.add_argument("--t_start_yr", type=float, default=1e3)
    p.add_argument("--t_end_yr", type=float, default=1e5)
    p.add_argument("--snaps_per_decade", type=int, default=10)

    return p.parse_args()


# ---------------------------------------------------------------------------
# Parameter override via --json
# ---------------------------------------------------------------------------

KNOWN_JSON_KEYS = {
    "star_mass_msun",
    "star_radius_rsun",
    "star_teff",
    "nr",
    "rmin_au",
    "rmax_au",
    "mmin_g",
    "alpha",
    "mdisk_msun",
    "sigma_exp",
    "sigma_rc_au",
    "d2g",
    "vfrag_cms",
    "rho_monomer",
    "a_ini_max_cm",
    "allow_drifting",
    "t_start_yr",
    "t_end_yr",
    "snaps_per_decade",
}


def apply_json_overrides(args, json_str: str):
    try:
        overrides = json.loads(json_str)
    except json.JSONDecodeError as e:
        print(f"ERROR message=Invalid --json: {e}", flush=True)
        sys.exit(1)
    unknown = set(overrides) - KNOWN_JSON_KEYS
    if unknown:
        print(
            f"ERROR message=Unknown JSON keys: {sorted(unknown)}. "
            f"Known keys: {sorted(KNOWN_JSON_KEYS)}",
            flush=True,
        )
        sys.exit(1)
    for key, val in overrides.items():
        setattr(args, key, val)
    return args


# ---------------------------------------------------------------------------
# SLURM jobscript writer
# ---------------------------------------------------------------------------

SLURM_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name=tripodpy_{run_name}
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem={memory_mb}M
#SBATCH --time={walltime}
#SBATCH --output={run_dir}/slurm_%j.out
#SBATCH --error={run_dir}/slurm_%j.err

# Activate environment if needed:
# source /path/to/venv/bin/activate

python {script_path} \\
    --run_dir {run_dir} \\
    --json '{params_json}'
"""


def write_slurm_script(args, params_json: str, run_dir: str) -> str:
    script_path = os.path.abspath(__file__)
    run_name = os.path.basename(run_dir.rstrip("/"))
    slurm_path = os.path.join(run_dir, "submit.sh")
    os.makedirs(run_dir, exist_ok=True)
    content = SLURM_TEMPLATE.format(
        run_name=run_name,
        partition=args.partition,
        memory_mb=int(args.memory_gb * 1024),
        walltime=args.walltime,
        run_dir=os.path.abspath(run_dir),
        script_path=script_path,
        params_json=params_json.replace("'", '"'),
    )
    with open(slurm_path, "w") as f:
        f.write(content)
    return slurm_path


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_params(args):
    errors = []
    if args.rmin_au <= 0 or args.rmax_au <= args.rmin_au:
        errors.append(f"Invalid grid: rmin={args.rmin_au} AU, " f"rmax={args.rmax_au} AU")
    if args.t_end_yr <= args.t_start_yr:
        errors.append(f"t_end_yr={args.t_end_yr} must be > " f"t_start_yr={args.t_start_yr}")
    if args.mmin_g <= 0:
        errors.append(f"mmin_g={args.mmin_g} must be positive")
    if args.t_start_yr <= 0:
        errors.append(f"t_start_yr={args.t_start_yr} must be > 0 (log-spaced snapshots)")
    if args.alpha <= 0 or args.alpha > 0.1:
        errors.append(f"alpha={args.alpha} out of physical range (0, 0.1]")
    if args.vfrag_cms <= 0:
        errors.append(f"vfrag_cms={args.vfrag_cms} must be positive")
    if errors:
        for e in errors:
            print(f"ERROR message={e}", flush=True)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Simulation setup and run
# ---------------------------------------------------------------------------


def build_and_run(args):
    from tripodpy import Simulation, readdump
    from tripodpy import constants as c

    run_dir = os.path.abspath(args.run_dir)
    datadir = os.path.join(run_dir, "data")
    os.makedirs(run_dir, exist_ok=True)

    sim = Simulation()

    # stellar
    sim.ini.star.M = args.star_mass_msun * c.M_sun
    sim.ini.star.R = args.star_radius_rsun * c.R_sun
    sim.ini.star.T = args.star_teff

    # grid
    sim.ini.grid.Nr = args.nr
    sim.ini.grid.rmin = args.rmin_au * c.au
    sim.ini.grid.rmax = args.rmax_au * c.au
    sim.ini.grid.mmin = args.mmin_g

    # gas
    sim.ini.gas.alpha = args.alpha
    sim.ini.gas.Mdisk = args.mdisk_msun * c.M_sun
    sim.ini.gas.SigmaExp = args.sigma_exp
    sim.ini.gas.SigmaRc = args.sigma_rc_au * c.au

    # dust
    sim.ini.dust.d2gRatio = args.d2g
    sim.ini.dust.vFrag = args.vfrag_cms
    sim.ini.dust.rhoMonomer = args.rho_monomer
    sim.ini.dust.aIniMax = args.a_ini_max_cm
    sim.ini.dust.allowDriftingParticles = bool(args.allow_drifting)

    # write a manifest so analysis workflow never depends on AGENTS.md alone
    manifest = {
        "code": "tripodpy",
        "run_dir": run_dir,
        "data_dir": datadir,
        "parameters": {
            "star_mass_msun": args.star_mass_msun,
            "alpha": args.alpha,
            "mdisk_msun": args.mdisk_msun,
            "nr": args.nr,
            "rmin_au": args.rmin_au,
            "rmax_au": args.rmax_au,
            "mmin_g": args.mmin_g,
            "a_ini_max_cm": args.a_ini_max_cm,
            "d2g": args.d2g,
            "vfrag_cms": args.vfrag_cms,
            "t_start_yr": args.t_start_yr,
            "t_end_yr": args.t_end_yr,
        },
    }
    import json as _json

    with open(os.path.join(run_dir, "run_manifest.json"), "w") as mf:
        _json.dump(manifest, mf, indent=2)

    if args.dry_run:
        n_snaps_dry = int(args.snaps_per_decade * np.log10(args.t_end_yr / args.t_start_yr)) + 1
        print(
            f"SUCCESS run_dir={run_dir} dry_run=True "
            f"n_snapshots={n_snaps_dry} t_end_yr={args.t_end_yr}",
            flush=True,
        )
        return

    snapshots = (
        np.geomspace(
            args.t_start_yr,
            args.t_end_yr,
            num=int(args.snaps_per_decade * np.log10(args.t_end_yr / args.t_start_yr)) + 1,
        )
        * c.year
    )

    if args.resume:
        # Dump files are pickles: only load dumps you created yourself.
        dump_path = os.path.join(datadir, "frame.dmp")
        if not os.path.exists(dump_path):
            print(
                f"ERROR message=Resume requested but dump file not found: " f"{dump_path}",
                flush=True,
            )
            sys.exit(1)
        sim = readdump(dump_path)  # restores the full frame incl. writer
        sim.writer.datadir = datadir
        # Append new snapshot times to the existing ones: simframe numbers the
        # output files by their index in sim.t.snapshots, so replacing the
        # array would restart at data0000.hdf5 and collide with existing files.
        old = np.asarray(sim.t.snapshots)
        t_now = max(float(sim.t), float(old.max()))
        t_end = args.t_end_yr * c.year
        if t_end > t_now:
            num = max(int(args.snaps_per_decade * np.log10(t_end / t_now)), 1) + 1
            sim.t.snapshots = np.concatenate([old, np.geomspace(t_now, t_end, num=num)[1:]])
        if not np.any(sim.t.snapshots > sim.t):
            print(f"ERROR message=Dump is already at t={float(sim.t)/c.year:.4g} yr >= t_end_yr", flush=True)
            sys.exit(1)
    else:
        sim.initialize()
        # writer and snapshots are set AFTER initialize(): sim.writer is None
        # and sim.t does not exist before it
        sim.writer.datadir = datadir
        sim.writer.overwrite = args.overwrite
        # include t=0 so data0000.hdf5 holds the initial condition
        sim.t.snapshots = np.hstack([sim.t, snapshots])

    t0 = time.perf_counter()
    sim.run()
    wall_clock_s = time.perf_counter() - t0

    # Count output files
    import glob as _glob

    hdf5_files = sorted(_glob.glob(os.path.join(datadir, "data*.hdf5")))
    n_snaps = len(hdf5_files)

    if n_snaps == 0:
        print(
            f"ERROR message=Run completed but no HDF5 files found in " f"{datadir}",
            flush=True,
        )
        sys.exit(1)

    # Verify last snapshot is readable
    import h5py as _h5

    last_file = hdf5_files[-1]
    try:
        with _h5.File(last_file, "r") as f:
            t_end_s = float(f["t"][()])
    except Exception as e:
        print(f"ERROR message=Last snapshot unreadable: {e}", flush=True)
        sys.exit(1)

    t_end_yr = t_end_s / c.year

    print(
        f"SUCCESS run_dir={run_dir} n_snapshots={n_snaps} "
        f"t_end_yr={t_end_yr:.4g} wall_clock_s={wall_clock_s:.1f}",
        flush=True,
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main():
    args = parse_args()

    if args.json:
        args = apply_json_overrides(args, args.json)

    validate_params(args)

    # Collect params for SLURM jobscript / manifest
    param_keys = [
        "star_mass_msun",
        "star_radius_rsun",
        "star_teff",
        "nr",
        "rmin_au",
        "rmax_au",
        "mmin_g",
        "alpha",
        "mdisk_msun",
        "sigma_exp",
        "sigma_rc_au",
        "d2g",
        "vfrag_cms",
        "rho_monomer",
        "a_ini_max_cm",
        "allow_drifting",
        "t_start_yr",
        "t_end_yr",
        "snaps_per_decade",
    ]
    params_dict = {k: getattr(args, k) for k in param_keys}
    params_json = json.dumps(params_dict)

    if args.hpc:
        slurm_path = write_slurm_script(args, params_json, args.run_dir)
        if args.dry_run:
            print(
                f"SUCCESS run_dir={os.path.abspath(args.run_dir)} "
                f"slurm_script={slurm_path} dry_run=True",
                flush=True,
            )
        else:
            import subprocess

            result = subprocess.run(["sbatch", slurm_path], capture_output=True, text=True)
            if result.returncode != 0:
                print(f"ERROR message=sbatch failed: {result.stderr.strip()}", flush=True)
                sys.exit(1)
            job_id = result.stdout.strip().split()[-1]
            print(
                f"SUCCESS run_dir={os.path.abspath(args.run_dir)} "
                f"slurm_job_id={job_id} slurm_script={slurm_path}",
                flush=True,
            )
        return

    try:
        build_and_run(args)
    except Exception:
        tb = traceback.format_exc()
        print(f"ERROR message={tb[-1800:]}", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
