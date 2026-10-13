# TriPoDPy Parameters Reference

All parameters in `sim.ini` must be set **before** `sim.initialize()`.
All quantities are in **CGS units**. Defaults below were read from
`tripodpy 1.0.0` (`Simulation().ini`) and from an initialised default model.

Table of contents:
1. [Stellar parameters — `sim.ini.star`](#1-stellar-parameters)
2. [Grid parameters — `sim.ini.grid`](#2-grid-parameters)
3. [Gas parameters — `sim.ini.gas`](#3-gas-parameters)
4. [Dust parameters — `sim.ini.dust`](#4-dust-parameters)
5. [TriPoD model parameters (set after `initialize()`)](#5-tripod-model-parameters)
6. [Time / snapshot parameters — `sim.t`](#6-time--snapshot-parameters)
7. [Writer parameters — `sim.writer`](#7-writer-parameters)
8. [Cross-parameter constraints and interactions](#8-cross-parameter-constraints-and-interactions)

---

## 1. Stellar parameters

Namespace: `sim.ini.star`

| Parameter | Default (CGS) | Physical | Notes |
|---|---|---|---|
| `M` | `1.9884e+33` g | 1 M☉ | Keplerian frequency, dynamical time scales |
| `R` | `1.3914e+11` cm | 2 R☉ | Enters the default irradiation temperature profile (via L) |
| `T` | `5772.0` K | — | Enters the default temperature profile (via L) |

Derived: `sim.star.L` = 4πR²σ_SB T⁴ [erg/s]. Default `gas.T` = (0.05 L / (4π r² σ_SB))^¼.

---

## 2. Grid parameters

Namespace: `sim.ini.grid` — contains **only** `Nr`, `rmin`, `rmax`, `mmin`.
(`Nmbpd` and `mmax` from DustPy are deleted; there is no mass grid.)

| Parameter | Default (CGS) | Physical | Notes |
|---|---|---|---|
| `Nr` | `100` | — | Number of radial cells (log-spaced between `rmin` and `rmax`) |
| `rmin` | `1.4960e+13` cm | 1 au | Inner boundary (interface) |
| `rmax` | `1.4960e+16` cm | 1000 au | Outer boundary (interface) |
| `mmin` | `1e-12` g | — | Minimum particle mass. Only used to set `dust.s.min = (3 mmin / (4π ρ_s fill))^(1/3)` → 5.23×10⁻⁵ cm for ρ_s = 1.67 |

Derived (after `makegrids()` / `initialize()`):

| Field | Shape | Unit | Description |
|---|---|---|---|
| `sim.grid.r` | `(Nr,)` | cm | Cell centres = arithmetic mean of adjacent interfaces |
| `sim.grid.ri` | `(Nr+1,)` | cm | Interfaces; set manually before `initialize()` for custom grids |
| `sim.grid.A` | `(Nr,)` | cm² | Annulus area π(ri₊² − ri₋²) |
| `sim.grid.OmegaK` | `(Nr,)` | 1/s | Keplerian frequency |
| `sim.grid._Nm_short` | scalar | — | 2: number of surface-density populations |
| `sim.grid._Nm_long` | scalar | — | 5: number of characteristic sizes |

---

## 3. Gas parameters

Namespace: `sim.ini.gas`

| Parameter | Default (CGS) | Physical | Notes |
|---|---|---|---|
| `alpha` | `1e-3` | — | Turbulent α. Initial value of `gas.alpha` **and** of `dust.delta.rad/turb/vert` |
| `Mdisk` | `9.942e+31` g | 0.05 M☉ | Initial gas disk mass (Lynden-Bell & Pringle 1974 profile) |
| `mu` | `3.847e-24` g | 2.3 m_p | Mean molecular weight; becomes `components.Default.gas.pars.mu` |
| `SigmaExp` | `-1.0` | — | Power-law index of initial Σ_gas |
| `SigmaRc` | `8.976e+14` cm | 60 au | Exponential cut-off radius |

> Note: `SigmaRc` default is 8.976e14 cm = **60 au** (not 30 au as in some older text).

Gas fields after `initialize()` (all shape `(Nr,)`):

| Field | Unit | Description |
|---|---|---|
| `sim.gas.Sigma` | g/cm² | **Sum of all components'** `gas.Sigma` (updater `tripodpy.std.gas.Sigma_tot`) |
| `sim.gas.SigmaFloor` | g/cm² | 1e-100 |
| `sim.gas.mu` | g | Σ / Σᵢ(Σᵢ/μᵢ) over components (updater `std.gas.mu`) |
| `sim.gas.T`, `cs`, `Hp`, `nu`, `eta`, `rho`, `P`, `n`, `mfp` | CGS | as in DustPy |
| `sim.gas.alpha` | — | No updater by default; can be overwritten with an array (e.g. gap profile) |
| `sim.gas.v.rad`, `v.visc`, `Fi` | cm/s, g/s | Computed after the implicit step (diagnostic) |
| `sim.gas.S.ext` | g/cm²/s | **Sum of components' `gas.Sigma_dot`** — assigning an updater here has no effect |
| `sim.gas.torque.Lambda`, `torque.v` | cm²/s², cm/s | Torque; `torque.v` enters the gas Jacobian (works in TriPoDPy) |

Default component (`sim.components.Default`, created by `initialize()`):

| Field | Shape | Description |
|---|---|---|
| `.gas.Sigma` | `(Nr,)` | Gas of this component (initially all gas) — this is what is integrated |
| `.gas._Sigma_dot` / `.gas.Sigma_dot` | `(Nr,)` | External gas source of this component [g/cm²/s]; zero by default |
| `.gas.pars.mu` | scalar | μ of this component (= `ini.gas.mu`) |
| `.gas.Fi`, `.gas.S.{ext,hyd,tot}` | `(Nr+1,)`, `(Nr,)` | Fluxes/sources (diagnostic) |

---

## 4. Dust parameters

Namespace: `sim.ini.dust`

| Parameter | Default | Notes |
|---|---|---|
| `d2gRatio` | `0.01` | Initial dust-to-gas ratio |
| `aIniMax` | `1e-4` cm | Initial `dust.s.max` (unless reduced by the drift limit, see below) |
| `allowDriftingParticles` | `False` | If `False`, initial `s.max` is capped at a drift-limited size (with safety margin) so no initial drift wave forms |
| `vFrag` | `100.0` cm/s | Initial `dust.v.frag` (constant, no updater) |
| `rhoMonomer` | `1.67` g/cm³ | Initial `dust.rhos` (all 5 sizes) |
| `distExp` | `-3.5` | **Not used by TriPoDPy** (inherited from DustPy; changing it does not change the initial state) |

### 4a. Dust fields after `initialize()`

| Field | Shape | Unit | Description |
|---|---|---|---|
| `dust.Sigma` | `(Nr, 2)` | g/cm² | [Σ₀ small (s.min–s_int), Σ₁ large (s_int–s.max)], s_int = √(s.min·s.max) |
| `dust.SigmaFloor` | `(Nr, 2)` | g/cm² | `1e-10` everywhere (empty cells sit at `0.1*SigmaFloor`) |
| `dust.s.min` | `(Nr,)` | cm | Minimum size, constant in time |
| `dust.s.max` | `(Nr,)` | cm | Maximum size — evolved; kept ≥ 1.5 `s.min` |
| `dust.s.lim` | scalar | cm | `1e-4`: limiting size for the shrinking term (Pfeil et al. 2024) |
| `dust.qrec` | `(Nr,)` | — | Exponent reconstructed from Σ₁/Σ₀: log(Σ₁/Σ₀)/log(s.max/s_int) − 4 |
| `dust.q.eff` | `(Nr,)` | — | Equilibrium exponent = q.frag·p.frag + q.sweep·(1−p.frag) |
| `dust.q.frag` | `(Nr,)` | — | Fragmentation exponent (between q.turb1, q.turb2, q.drfrag) |
| `dust.q.sweep`, `q.turb1`, `q.turb2`, `q.drfrag` | scalars | — | −3.0, −3.5, −3.75, −3.75 |
| `dust.a` | `(Nr, 5)` | cm | `[a0, f.dv·a1, a1, f.dv·s.max, s.max]` — a0, a1 mass-averaged sizes of the two populations |
| `dust.m` | `(Nr, 5)` | g | (4/3)π ρ_s fill a³ |
| `dust.St` | `(Nr, 5)` | — | Stokes numbers (Epstein + Stokes I) of the 5 sizes |
| `dust.H` | `(Nr, 5)` | cm | Dust scale heights (Dubrulle et al. 1995) |
| `dust.rho` | `(Nr, 5)` | g/cm³ | Midplane density, `Sigma[:, [0,1,1,1,1]] / (√(2π) H)` |
| `dust.rhos`, `dust.fill` | `(Nr, 5)` | g/cm³, — | Bulk density and filling factor. **Keep both size-independent** — TriPoD assumes ρ(a) = const |
| `dust.D` | `(Nr, 5)` | cm²/s | Diffusivity (uses `St·f.drift`; zero in the first/last two cells) |
| `dust.eps` | `(Nr,)` | — | Σ_d,tot / Σ_gas |
| `dust.v.frag` | `(Nr,)` | cm/s | Fragmentation velocity |
| `dust.v.rad` | `(Nr, 5)` | cm/s | Radial velocity; `v.rad_flux` includes `f.drift` and is used for the fluxes |
| `dust.v.rel.{azi,brown,rad,turb,vert,tot}` | `(Nr, 5, 5)` | cm/s | Relative velocities |
| `dust.v.driftmax` | `(Nr,)` | cm/s | Max drift velocity |
| `dust.p.frag`, `p.stick`, `p.fragtrans`, `p.driftfrag` | `(Nr,)` | — | Probabilities/transition functions (p.frag uses `v.rel.tot[:, -2, -1]`, i.e. 0.4 s.max vs s.max) |
| `dust.S.{ext,hyd,coag,compo,tot}` | `(Nr, 2)` | g/cm²/s | Source terms; `S.ext` zero by default |
| `dust.S.smax_hyd` | `(Nr,)` | cm/s | Hydrodynamic source of s.max |
| `dust.s.sdot_coag` | `(Nr,)` | cm/s | Coagulation growth rate of s.max |
| `dust.Fi.{adv,diff,tot}` | `(Nr+1, 2)` | g/cm/s | Interface fluxes (computed after the implicit step; diagnostic only) |
| `dust.backreaction.A`, `.B` | `(Nr,)` | — | 1 and 0 (backreaction off) |

### 4b. Mixing parameters

Set directly on the frame (no updater by default), e.g. after `initialize()`:

| Field | Default | Used for |
|---|---|---|
| `sim.dust.delta.rad` | `ini.gas.alpha` | radial diffusion |
| `sim.dust.delta.turb` | `ini.gas.alpha` | turbulent relative velocities |
| `sim.dust.delta.vert` | `ini.gas.alpha` | settling / scale height |

```python
sim.initialize()
sim.dust.delta.rad[:]  = 1e-4
sim.dust.delta.turb[:] = 1e-4
sim.update()
```

---

## 5. TriPoD model parameters

Exist only after `initialize()`; change them afterwards (then `sim.update()`).
These calibrations come from Pfeil et al. (2024); change only with good reason.

| Field | Default | Meaning |
|---|---|---|
| `sim.dust.f.drift` | `0.8` | Drift velocity calibration factor (multiplies St in `v.rad_flux` and `D`) |
| `sim.dust.f.dv` | `0.4` | Effective collision-speed parameter; also the size ratio in `dust.a` columns 1 and 3 |
| `sim.dust.f.crit` | `0.425` | Critical mass fraction in the large population (shrinking criterion) |
| `sim.dust.s.lim` | `1e-4` cm | Limiting size for shrinking |
| `sim.dust.q.turb1/turb2/drfrag/sweep` | −3.5 / −3.75 / −3.75 / −3.0 | Regime exponents |
| `sim.dust.s.boundary.inner/outer` | const. gradient / value | Boundary conditions for `s.max · Σ₁` |

---

## 6. Time / snapshot parameters

`sim.t` is created in `initialize()` (it is `None` before). Set snapshots afterwards.

| Attribute | Default | Notes |
|---|---|---|
| `sim.t.snapshots` | 21 log-spaced times 10³–10⁵ yr (in s) | Output times; include `sim.t` (= 0) to also write the initial state |
| `sim.t.cfl` | `0.1` | Safety factor of the timestep |

```python
import numpy as np
import tripodpy.constants as c
sim.initialize()
sim.t.snapshots = np.hstack([sim.t, np.geomspace(1e3, 1e6, num=31) * c.year])  # 10 per decade
```

Output files are numbered by the index in `sim.t.snapshots`. When extending a
restarted run, **append** to the existing array (replacing it restarts numbering at
`data0000.hdf5` and fails with "File ... already exists").

---

## 7. Writer parameters

`sim.writer` is `None` before `initialize()`; afterwards it is a `simframe` `hdf5writer`.

| Attribute | Default | Notes |
|---|---|---|
| `sim.writer.datadir` | `"data"` | Created automatically |
| `sim.writer.overwrite` | `False` | If `False`, writing an existing file raises `RuntimeError`. Set `True` only after explicit user confirmation |
| `sim.writer = None` | — | Run without writing any output (useful for tests) |
| `<field>.save = False` | `True` | Exclude a field/group from the files, e.g. `sim.dust.v.save = False` |

Files: `data{N:04d}.hdf5` (one per snapshot) and `frame.dmp` (latest dump, for `tripodpy.readdump`).

---

## 8. Cross-parameter constraints and interactions

**`alpha` sets four things by default**: gas viscosity and the three `dust.delta.*`.
Changing `sim.gas.alpha` after `initialize()` does **not** change `dust.delta.*`.

**`vFrag`, `alpha`, `T` set the fragmentation limit** St_frag ≈ v_frag² / (3 δ_turb c_s²)
(`read_data` returns it as `dust.StFr`). Larger `vFrag` → larger `s.max`.

**`aIniMax` vs `allowDriftingParticles`**: with the default `False`, `s.max` starts at
min(`aIniMax`, drift-limited size) — check `sim.dust.s.max` after `initialize()`.

**`mmin` and `rhoMonomer` set `s.min`**; `s.max` is never allowed below 1.5 `s.min`.

**`SigmaRc` vs `rmax`**: beyond a few `SigmaRc` the gas falls to the floor (10⁻¹⁰⁰) and the
dust to `0.1·SigmaFloor`. Those cells have meaningless `eps`/`St` — mask them in analysis.

**Bulk density / filling factor must stay size-independent** (TriPoD assumes ρ(a) = const).
Active dust components change `dust.rhos` radially (same value for all 5 sizes) — that is fine.
