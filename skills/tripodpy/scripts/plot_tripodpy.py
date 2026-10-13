#!/usr/bin/env python3
"""
plot_tripodpy.py — Diagnostic plots for TriPoDPy HDF5 output.

Produces three plot types (selectable via --plot):
  panel     : overview figure produced by ``tripodpy.plot.panel`` (reconstructed
              size-distribution map, size distribution at one radius, gas/dust mass
              evolution, midplane density + exponent, surface densities + eps)
  radial    : radial profiles of Sigma_gas / Sigma_dust, eps, s_max, St(s_max),
              q_rec and the two populations Sigma_0 / Sigma_1 at selected snapshots
  evolution : space-time diagrams of s_max(r, t) and Sigma_dust(r, t)

Usage
-----
  python plot_tripodpy.py --run_dir runs/my_run [options]

  # All three plot types:
  python plot_tripodpy.py --run_dir runs/my_run --plot all --out plots/my_run

  # Radial profiles at first and last snapshot:
  python plot_tripodpy.py --run_dir runs/my_run --plot radial --snaps 0 -1

Data directory: ``<run_dir>/data`` if it exists, otherwise ``<run_dir>`` itself
(must contain ``data*.hdf5``).

Output (stdout, last line):
  SUCCESS plot_dir=<path>  files=<comma-separated list>
  ERROR   message=<description>
"""

import argparse
import glob
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # must be set before pyplot is imported

import h5py  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import tripodpy.constants as c  # noqa: E402

# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def _hdf5_files(datadir: str) -> list[str]:
    files = sorted(glob.glob(os.path.join(datadir, "data*.hdf5")))
    if not files:
        raise FileNotFoundError(f"No data*.hdf5 files in {datadir}")
    return files


def _read_snap(path: str) -> dict:
    """Read one TriPoDPy HDF5 snapshot. Returns dict of CGS arrays.

    Shapes: Nr radial cells, 2 surface-density populations, 5 particle sizes.
    """
    with h5py.File(path, "r") as f:
        r_cm = f["grid/r"][:]  # (Nr,) cm
        t_s = float(f["t"][()])  # s
        Sigma_gas = f["gas/Sigma"][:]  # (Nr,) g/cm²
        Sigma_dust = f["dust/Sigma"][:]  # (Nr, 2) g/cm²  [Sigma_0, Sigma_1]
        SigmaFloor = f["dust/SigmaFloor"][:]  # (Nr, 2) g/cm²
        St = f["dust/St"][:]  # (Nr, 5)
        smax = f["dust/s/max"][:]  # (Nr,) cm
        eps = f["dust/eps"][:]  # (Nr,)
        qrec = f["dust/qrec"][:]  # (Nr,)
    Sigma_dust_tot = Sigma_dust.sum(axis=1)
    # Cells where the dust sits at/below its floor value carry no physical
    # information (eps/St there can be ~1e90) — mask them for plotting.
    has_dust = Sigma_dust_tot > SigmaFloor.sum(axis=1)
    return dict(
        r_au=r_cm / c.au,
        t_yr=t_s / c.year,
        Sigma_gas=Sigma_gas,
        Sigma_dust=Sigma_dust,
        Sigma_dust_tot=Sigma_dust_tot,
        has_dust=has_dust,
        St_max=St[:, -1],  # Stokes number of s_max
        smax_cm=smax,
        eps=eps,
        qrec=qrec,
    )


def _masked(x, mask):
    return np.where(mask, x, np.nan)


# ---------------------------------------------------------------------------
# Plot: overview panel (tripodpy.plot.panel, saved to file)
# ---------------------------------------------------------------------------


def plot_panel(datadir: str, snap_index: int, ir: int, out_path: str) -> str:
    """Save the tripodpy.plot.panel figure to a file without displaying it."""
    from tripodpy import plot

    n = len(_hdf5_files(datadir))
    # panel() clips negative indices to 0, so resolve them here
    it = snap_index % n
    plt.close("all")
    plot.panel(datadir, it=it, ir=ir, show_limits=True)  # calls plt.show(): no-op under Agg
    fig = plt.gcf()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# Plot: radial profiles
# ---------------------------------------------------------------------------


def plot_radial(datadir: str, snap_indices: list[int], out_path: str) -> str:
    """Radial profiles at selected snapshots."""
    files = _hdf5_files(datadir)
    n = len(files)
    snaps = [i % n for i in snap_indices]

    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    ax_sig, ax_eps, ax_smax, ax_st, ax_q, ax_pop = axes.flat

    cmap = plt.get_cmap("viridis")
    colors = [cmap(i / max(len(snaps) - 1, 1)) for i in range(len(snaps))]

    for color, idx in zip(colors, snaps):
        s = _read_snap(files[idx])
        lbl = f"t = {s['t_yr']:.2e} yr"
        r = s["r_au"]
        m = s["has_dust"]

        ax_sig.loglog(r, s["Sigma_gas"], color=color, ls="-", label=lbl)
        ax_sig.loglog(r, _masked(s["Sigma_dust_tot"], m), color=color, ls="--")

        ax_eps.loglog(r, _masked(s["eps"], m), color=color, label=lbl)
        ax_smax.loglog(r, _masked(s["smax_cm"], m), color=color, label=lbl)
        ax_st.loglog(r, _masked(s["St_max"], m), color=color, label=lbl)
        ax_q.semilogx(r, _masked(s["qrec"], m), color=color, label=lbl)

        ax_pop.loglog(r, _masked(s["Sigma_dust"][:, 0], m), color=color, ls=":")
        ax_pop.loglog(r, _masked(s["Sigma_dust"][:, 1], m), color=color, ls="-", label=lbl)

    # the gas floor (1e-100) would otherwise stretch the axis over ~100 decades
    ax_sig.set_ylim(bottom=1e-12 * ax_sig.get_ylim()[1])
    ax_sig.set_ylabel(r"$\Sigma$ [g cm$^{-2}$]")
    ax_sig.set_title("Surface density (solid=gas, dashed=dust)")
    ax_sig.legend(fontsize=7)

    ax_eps.set_ylabel(r"$\varepsilon = \Sigma_\mathrm{d}/\Sigma_\mathrm{g}$")
    ax_eps.set_title("Dust-to-gas ratio")

    ax_smax.set_ylabel(r"$s_\mathrm{max}$ [cm]")
    ax_smax.set_title("Maximum particle size")

    ax_st.set_ylabel(r"St($s_\mathrm{max}$)")
    ax_st.axhline(1.0, color="r", ls=":", lw=1)
    ax_st.set_title("Stokes number of largest grains")

    ax_q.set_ylabel(r"$q_\mathrm{rec}$")
    ax_q.set_title("Size-distribution exponent")

    ax_pop.set_ylabel(r"$\Sigma$ [g cm$^{-2}$]")
    ax_pop.set_title(r"$\Sigma_0$ (dotted) and $\Sigma_1$ (solid)")

    for ax in axes.flat:
        ax.set_xlabel("Distance [au]")

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# Plot: space-time evolution
# ---------------------------------------------------------------------------


def plot_evolution(datadir: str, out_path: str) -> str:
    """Space-time diagrams of s_max(r,t) and Sigma_dust(r,t)."""
    files = _hdf5_files(datadir)
    snaps = [_read_snap(f) for f in files]
    t_yr = np.array([s["t_yr"] for s in snaps])
    if t_yr[0] <= 0.0 and len(t_yr) > 1:  # drop t=0 so the log time axis works
        snaps, t_yr = snaps[1:], t_yr[1:]
    r_au = snaps[-1]["r_au"]

    smax_grid = np.array([_masked(s["smax_cm"], s["has_dust"]) for s in snaps])
    Sigma_d_grid = np.array([_masked(s["Sigma_dust_tot"], s["has_dust"]) for s in snaps])

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12, 5))

    def _pcolormesh_log(ax, z, cmap, label):
        z_safe = np.where(z > 0, z, np.nan)
        pm = ax.pcolormesh(r_au, t_yr, np.log10(z_safe), cmap=cmap, shading="auto")
        cb = fig.colorbar(pm, ax=ax)
        cb.set_label(label)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Distance [au]")
        ax.set_ylabel("Time [yr]")

    _pcolormesh_log(ax0, smax_grid, "magma", r"$\log_{10}\,s_\mathrm{max}$ [cm]")
    ax0.set_title("Maximum particle size")
    _pcolormesh_log(ax1, Sigma_d_grid, "viridis", r"$\log_{10}\,\Sigma_\mathrm{dust}$ [g cm$^{-2}$]")
    ax1.set_title("Dust surface density")

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--run_dir",
        required=True,
        help="Run directory (contains data/ with HDF5 files) or the data directory itself",
    )
    p.add_argument("--plot", default="all", choices=["panel", "radial", "evolution", "all"])
    p.add_argument(
        "--snaps", type=int, nargs="+", default=[0, -1],
        help="Snapshot indices for the radial plot (negative = from end)",
    )
    p.add_argument("--snap", type=int, default=-1, help="Snapshot index for the panel plot (negative = from end)")
    p.add_argument("--ir", type=int, default=20, help="Radial index for the size-distribution panel")
    p.add_argument("--out", default=None, help="Output directory (default: plots/<run_name>)")
    p.add_argument("--fmt", default="pdf", choices=["pdf", "png", "both"])
    return p.parse_args()


def main() -> None:
    args = parse_args()

    run_dir = os.path.abspath(args.run_dir)
    datadir = os.path.join(run_dir, "data")
    if not os.path.isdir(datadir):
        datadir = run_dir
    if not glob.glob(os.path.join(datadir, "data*.hdf5")):
        print(f"ERROR message=no data*.hdf5 files in {run_dir} or {run_dir}/data", flush=True)
        sys.exit(1)

    run_name = Path(run_dir).name
    out_dir = os.path.abspath(args.out or os.path.join("plots", run_name))
    os.makedirs(out_dir, exist_ok=True)

    fmts = ["pdf", "png"] if args.fmt == "both" else [args.fmt]
    made = []

    try:
        plots_requested = ["panel", "radial", "evolution"] if args.plot == "all" else [args.plot]
        for fmt in fmts:
            if "panel" in plots_requested:
                out = os.path.join(out_dir, f"{run_name}_panel.{fmt}")
                made.append(plot_panel(datadir, args.snap, args.ir, out))
            if "radial" in plots_requested:
                out = os.path.join(out_dir, f"{run_name}_radial.{fmt}")
                made.append(plot_radial(datadir, args.snaps, out))
            if "evolution" in plots_requested:
                out = os.path.join(out_dir, f"{run_name}_evolution.{fmt}")
                made.append(plot_evolution(datadir, out))
        print(f"SUCCESS plot_dir={out_dir}  files={','.join(made)}", flush=True)
    except Exception:
        import traceback

        print(f"ERROR message={traceback.format_exc()}", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
