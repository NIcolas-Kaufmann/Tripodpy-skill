# TriPoDPy Customisation Reference

Every recipe below was run against `tripodpy 1.0.0` on a small grid
(Nr = 40, ≤ 5×10⁴ yr) and its effect checked. Where a DustPy recipe does **not**
work in TriPoDPy, that is stated explicitly.

Table of contents:
1. [Replacing a physics updater](#1-replacing-a-physics-updater)
2. [Custom radial grid / grid refinement](#2-custom-radial-grid--grid-refinement)
3. [Turning off physics](#3-turning-off-physics)
4. [External sources and sinks](#4-external-sources-and-sinks)
5. [Planetary gaps](#5-planetary-gaps)
6. [Ice lines / spatially varying fragmentation velocity](#6-ice-lines--spatially-varying-fragmentation-velocity)
7. [Planetesimal formation sink + mass bookkeeping](#7-planetesimal-formation-sink)
8. [Per-timestep hooks (systoles / diastoles) and custom fields](#8-per-timestep-hooks-and-custom-fields)
9. [Boundary conditions](#9-boundary-conditions)
10. [Snapshots, output and restarts](#10-snapshots-output-and-restarts)

**Order of operations.** Most fields (everything in `sim.dust`, `sim.gas`, `sim.t`,
`sim.writer`, `sim.components`) only exist after `sim.initialize()`. So the usual
pattern is:

```python
sim = Simulation()
# 1. sim.ini.* parameters and (optionally) sim.grid.ri
sim.initialize()
# 2. assign updaters / overwrite fields / add components
sim.update()          # propagate changes to dependent fields
# 3. writer + snapshots, then sim.run()
```

Assigning `sim.dust.S.ext.updater = ...` *before* `initialize()` raises
`AttributeError: 'NoneType' object has no attribute 'updater'`.

---

## 1. Replacing a physics updater

Every `Field` has an `.updater`. Assigning a function `f(sim) -> ndarray` (correct
shape!) replaces the physics; assigning `None` freezes the field at its current value.

```python
import numpy as np
from tripodpy import Simulation
import tripodpy.constants as c

def my_temperature(sim):
    """Power-law midplane temperature, shape (Nr,)."""
    return 150. * (sim.grid.r / c.au) ** -0.5      # [K]

sim = Simulation()
sim.initialize()
sim.gas.T.updater = my_temperature
sim.update()
```

Rules:
- The function is called every timestep — keep it cheap.
- To set a field once, assign values and then `.updater = None`.
- Assigning a new updater replaces the heartbeat incl. its systole/diastole; to keep
  them, assign to `field.updater.updater` instead (as in the docs notebooks).
- Shapes: dust surface-density-like fields `(Nr, 2)`, size-like fields `(Nr, 5)`,
  relative velocities `(Nr, 5, 5)`, everything else `(Nr,)`.
- Keep `dust.rhos` and `dust.fill` **independent of particle size** (the TriPoD
  reconstruction assumes ρ(a) = const).

---

## 2. Custom radial grid / grid refinement

Set the interfaces `sim.grid.ri` (plain ndarray, `Nr+1` values) before
`sim.initialize()` (or `sim.makegrids()`); `Nr`, `r` (arithmetic mid-points) and `A`
are derived. `ini.grid.Nr/rmin/rmax` are then ignored.

```python
def refinegrid(ri, r0, num=3):
    """Refine the radial grid around r0 (from the TriPoDPy docs, notebook 3)."""
    if num == 0:
        return ri
    for _ in range(num):
        ind = np.argmin(r0 > ri) - 1
        ind_left, ind_right = ind - num, ind + num + 1
        n = (2 * num + 1) * 2
        ri_mid = np.empty(n)
        for i in range(0, n, 2):
            j = ind - num + int(i / 2)
            ri_mid[i] = ri[j]
            ri_mid[i + 1] = 0.5 * (ri[j] + ri[j + 1])
        ri = np.concatenate((ri[:ind_left], ri_mid, ri[ind_right:]))
    return ri

sim = Simulation()
sim.grid.ri = refinegrid(np.logspace(0., 3., num=100) * c.au, 4.5 * c.au, num=3)
sim.initialize()
print(int(sim.grid.Nr))     # > 99: refined cells were added
```

There is no mass grid to customise: the dust always has 2 populations and 5
characteristic sizes. To change the smallest grains use `ini.grid.mmin`
(sets `dust.s.min`).

---

## 3. Turning off physics

### Freeze the gas (verified: Σ_gas changes < 0.2 % in 3×10³ yr)
```python
sim.initialize()
sim.gas.nu[:] = 0.
sim.gas.nu.updater = None
```
Do **not** use the DustPy recipe `sim.gas.Sigma.updater = None`: it freezes the
*total* `gas.Sigma` while the `Default` component keeps evolving underneath, so the
two become inconsistent.

### Turn off dust transport (drift + diffusion), keep growth
The implicit dust solver reads `dust.v.rad_flux` and `dust.D` directly, so
freezing `dust.Fi.*` (DustPy recipe) does **nothing** — the fluxes are only diagnostics.
```python
sim.initialize()
sim.dust.v.rad_flux[:] = 0.; sim.dust.v.rad_flux.updater = None
sim.dust.D[:] = 0.;          sim.dust.D.updater = None
sim.update()
```

### Turn off coagulation (growth and fragmentation), keep transport
Setting only `dust.S.coag` to zero is **not** enough (the Jacobian builds its own
coagulation terms). Zero the relative velocities instead — this stops `s.max` growth
and the Σ₀ ↔ Σ₁ exchange (verified: changes < 0.1 %):
```python
sim.initialize()
sim.dust.v.rel.tot[:] = 0.
sim.dust.v.rel.tot.updater = None
sim.update()
```

---

## 4. External sources and sinks

Units g cm⁻² s⁻¹; negative = sink.

### Dust
`sim.dust.S.ext` has shape `(Nr, 2)` and is added to both populations. Assign the
updater **after** `initialize()`:

```python
def dust_sink(sim):
    """Remove dust inside 5 au on a 100-yr timescale."""
    S = np.zeros_like(sim.dust.Sigma)
    m = sim.grid.r < 5. * c.au
    S[m, :] = -sim.dust.Sigma[m, :] / (1e2 * c.year)
    return S

sim.initialize()
sim.dust.S.ext.updater = dust_sink
sim.update()
```

If the simulation has **active dust components** (`dust_active=True`, see
`components.md`), `sim.dust.Sigma` is overwritten by the sum of the components each
step and `sim.dust.S.ext` has no effect. Put the source on a component instead:

```python
bg = sim.components.bg.dust                        # an active dust component
bg._Sigma_dot.updater = dust_sink                  # (Nr, 2) [g/cm²/s]
bg.updater = ["_Sigma_dot"] + list(bg.updateorder)
sim.update()
```

### Gas
`sim.gas.S.ext.updater = ...` (DustPy way) is **silently ignored**: TriPoDPy sets
`gas.S.ext` to the sum of the components' `gas.Sigma_dot`, and only those enter the
solver. Put gas sources on a component — normally `Default`:

```python
def disk_wind(sim):
    """Mass loss of the Default gas component, (Nr,) [g/cm²/s]."""
    return -1e-5 * sim.components.Default.gas.Sigma * sim.grid.OmegaK

sim.initialize()
g = sim.components.Default.gas
g._Sigma_dot.updater = disk_wind
g.updater = ["_Sigma_dot"] + list(g.updateorder)   # make sure it is updated every step
sim.update()
```
(Verified: an e-folding sink of 10⁴ yr gave Σ/Σ_ref = 0.741 = e^−0.3 after 3×10³ yr.)
`sim.gas.S.ext` then shows the total for diagnostics.

---

## 5. Planetary gaps

### A. α-profile gap (recommended; verified)
Lowering the viscosity locally, α(r) = α₀ / f(r), produces a gap with
Σ ∝ f(r) in steady state. `dustpylib` provides f(r) (e.g. Duffell 2020; Kanagawa
et al. 2017: `kanagawa2017` from the same module).

```python
from dustpylib.substructures.gaps import duffell2020

sim = Simulation()
sim.initialize()

a_p   = 10. * c.au
q_p   = 1e-3                                   # planet/star mass ratio
h_p   = np.interp(a_p, sim.grid.r, sim.gas.Hp / sim.grid.r)
alpha0 = sim.ini.gas.alpha

f_gap = duffell2020(sim.grid.r, a_p, q_p, h_p, alpha0)   # Sigma perturbation f <= 1
sim.gas.alpha[:] = alpha0 / f_gap
sim.update()
```
Verified: Σ at a_p dropped to ~0.3 % of the unperturbed profile after 5×10⁴ yr
(Nr = 60). `dust.delta.*` stay at α₀ unless you change them too.
**Note**: `f_gap` is a *density* perturbation — use `alpha0 / f_gap`, not
`alpha0 * f_gap`.

To start with the gap already open, also multiply `sim.gas.Sigma`,
`sim.components.Default.gas.Sigma` and `sim.dust.Sigma` by `f_gap` (dust:
`[:, None]`) before running.

### B. Torque
`sim.gas.torque.Lambda` (specific angular-momentum injection rate, cm² s⁻²) is
converted to `sim.gas.torque.v` by its default updater, and `torque.v` enters the gas
Jacobian. Any Λ(r) can be supplied:

```python
sim.initialize()
def Lambda(sim):
    r = sim.grid.r
    return my_torque_profile(r)            # [cm²/s²], shape (Nr,) — user-supplied physics
sim.gas.torque.Lambda.updater = Lambda
sim.update()
```
(Verified with a test profile: Σ near a_p fell to ~8 % after 2×10⁴ yr.) The physical
torque prescription is up to the user; do not invent one — ask or cite a source.

---

## 6. Ice lines / spatially varying fragmentation velocity

```python
def vfrag_icelines(sim):
    """Pinilla et al. (2017)-style ice lines, (Nr,) [cm/s]."""
    T = sim.gas.T
    v = np.where(T < 150., 1000., 100.)    # water
    v = np.where(T <  80.,  700., v)       # NH3
    v = np.where(T <  44.,  100., v)       # CO2
    return v

sim.initialize()
sim.dust.v.frag.updater = vfrag_icelines   # evolves with T
sim.update()
```
For a static profile: `sim.dust.v.frag[:] = vfrag_icelines(sim)` (default updater is
`None`). Verified: `s.max` inside/outside the snow line responds as expected.

---

## 7. Planetesimal formation sink

```python
eps_crit = 1.0        # dust-to-gas threshold
zeta     = 0.01       # conversion efficiency per orbit

sim.initialize()
sim.grid.addfield("Mplts", 0., description="Cumulative planetesimal mass [g]")

def planetesimal_sink(sim):
    S = np.zeros_like(sim.dust.Sigma)
    over = sim.dust.eps > eps_crit                       # eps has shape (Nr,)
    S[over, :] = -zeta * sim.dust.Sigma[over, :] * sim.grid.OmegaK[over, None]
    return S
sim.dust.S.ext.updater = planetesimal_sink

def accumulate(sim):
    dt = sim.t.prevstepsize
    if dt:
        sim.grid.Mplts += -float((sim.dust.S.ext.sum(-1) * sim.grid.A).sum()) * dt
sim.updater.diastole = accumulate    # runs after every sim.update(), i.e. every step

sim.update()
```
`grid.Mplts` is a normal field, so it is written to every HDF5 file
(`grid/Mplts`). Verified in a test with a low threshold: booked mass agreed with the
dust mass removed to ~2 % (remainder = flux through the inner boundary).
With active dust components, put the sink on the components (§4).

---

## 8. Per-timestep hooks and custom fields

There is **no** `sim.addSystole` / `addDiastole` / `addSystemicUpdate`. Every updater
(Heartbeat) has a `systole` (before) and `diastole` (after):

```python
sim.updater.systole  = f    # before the full-frame update (each step)
sim.updater.diastole = g    # after the full-frame update (each step)
sim.gas.T.updater.diastole = h   # right after gas.T is updated
```
Verified: `sim.updater.systole/diastole` fire once per integration step.

Custom fields (written to the HDF5 output automatically):
```python
sim.grid.addfield("rsnow", 0., description="Snowline location [cm]")
sim.grid.rsnow.updater = lambda sim: sim.grid.ri[np.argmax(sim.gas.T < 150.)]
sim.gas.T.updater.diastole = lambda sim: sim.grid.rsnow.update()   # needs T first
```
Assigning a list to `group.updater` sets the update order (and clears its
systole/diastole), e.g. `sim.grid.updater = ["OmegaK", "rsnow"]`.

---

## 9. Boundary conditions

Defaults: gas and dust inner = constant gradient, outer = value (floor);
`dust.s.boundary` (acts on `s.max · Σ₁`) inner = constant gradient, outer = value.

```python
sim.initialize()
# dust: sim.dust.boundary acts on the solver (verified)
sim.dust.boundary.inner.setcondition("grad", -1)          # needs a value
# gas: set it on the COMPONENT(s) — sim.gas.boundary is NOT used by the solver
sim.components.Default.gas.boundary.inner.setcondition("val", 1e3)
```
Conditions (`dustpy.utils.Boundary.setcondition`): `"const_grad"`, `"const_pow"`,
`"const_val"`, `"val"`, `"grad"`, `"pow"` (the last three need a value), `None`.

**Gas boundaries:** the gas is integrated per component, and only
`sim.components.<name>.gas.boundary` enters the solver. Verified: setting
`sim.gas.boundary.inner` to `"val", 1e4` had no effect, setting
`sim.components.Default.gas.boundary.inner` to the same pinned Σ_gas[0] = 10⁴.
With several active gas components, set the boundary on each of them.

For `s.max`, replace the boundary object as in the repo tests:

```python
from dustpy.utils import Boundary
sim.dust.s.boundary.outer = Boundary(sim.grid.r[::-1], sim.grid.ri[::-1],
                                     sim.dust.Sigma[::-1, 1] * sim.dust.s.max[::-1],
                                     condition="const_pow")
```

---

## 10. Snapshots, output and restarts

```python
sim.initialize()
sim.t.snapshots = np.hstack([sim.t, np.geomspace(1e3, 1e6, num=16) * c.year])  # 5/decade
sim.writer.datadir  = "runs/my_run/data"
sim.writer.overwrite = False
sim.dust.v.save = False      # don't write the (large) velocity group
sim.run()
```

Run without any output (tests / quick experiments): `sim.writer = None`.

Restart from the dump file (pickle — trusted files only):
```python
from tripodpy import readdump
sim = readdump("runs/my_run/data/frame.dmp")
sim.t.snapshots = np.concatenate([sim.t.snapshots, [2e6 * c.year]])   # APPEND, don't replace
sim.run()
```
Files are numbered by the index in `t.snapshots`; replacing the array makes the writer
try to write `data0000.hdf5` again (`RuntimeError: ... already exists`).
