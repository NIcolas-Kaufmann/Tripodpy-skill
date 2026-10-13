# TriPoDPy skill for Claude Code

A Claude Code plugin with one Agent Skill, `tripodpy`. The skill helps you set up, run,
customise and analyse [TriPoDPy](https://github.com/tripod-code/tripodpy) simulations of
dust evolution in protoplanetary disks.

The skill covers:

- standard setups, snapshots, output and restarts
- reading TriPoDPy HDF5 output, plotting, and sanity checks
- customisation: ice lines, planetary gaps, sources and sinks, boundaries, and switching
  off transport, coagulation or gas evolution
- compositional tracking with `sim.addcomponent`: tracers, several dust species, and
  sublimation/condensation
- known pitfalls when carrying DustPy recipes over to TriPoDPy

Everything in the skill was checked against `tripodpy 1.0.0`.

## Installation

Install from GitHub inside Claude Code:

```
/plugin marketplace add NIcolas-Kaufmann/Tripodpy-skill
/plugin install tripodpy@tripodpy-skill
```

Alternatively, copy `skills/tripodpy/` into `~/.claude/skills/` to use it as a personal skill.

## Requirements

The bundled scripts need a Python environment with:

- `tripodpy` (Python ≥ 3.11; it depends on `dustpy`), installed with `pip install tripodpy`
- `numpy`, `h5py`, `matplotlib`
- optional: `dustpylib`, used by the planetary-gap recipe

## Layout

```
.claude-plugin/
  plugin.json          plugin manifest
  marketplace.json     single-plugin marketplace (this repo)
skills/tripodpy/
  SKILL.md             entry point loaded by Claude
  references/          parameters, customisation, components, handoff
  scripts/
    run_tripodpy.py    run, resume, or write a SLURM job for a simulation
    plot_tripodpy.py   overview, radial and space-time plots from HDF5 output
```

## Citation

If you use TriPoDPy, cite:

- Kaufmann et al. (2025), the TriPoDPy code
- Pfeil et al. (2024), A&A 691, A45, the TriPoD method

## License

GPL-3.0, the same license as TriPoDPy. See `LICENSE`.
