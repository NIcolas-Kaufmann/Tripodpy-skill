# tripodpy Customisation Reference

Table of contents:
1. [Replacing a physics updater](#1-replacing-a-physics-updater)
2. [Custom radial grid / grid refinement](#2-custom-radial-grid--grid-refinement)
3. [Turning off physics modules](#3-turning-off-physics-modules)
4. [External sources and sinks](#4-external-sources-and-sinks)
5. [Planetary gaps (Kanagawa 2017 torque)](#5-planetary-gaps-kanagawa-2017-torque)
6. [Ice lines / spatially varying fragmentation velocity](#6-ice-lines--spatially-varying-fragmentation-velocity)
7. [Planetesimal formation sink](#7-planetesimal-formation-sink)
8. [Systoles and diastoles (per-timestep / per-output hooks)](#8-systoles-and-diastoles)
9. [Boundary conditions](#9-boundary-conditions)
10. [Snapshot customisation](#10-snapshot-customisation)

All customisation must happen **after** `sim = Simulation()` and **before**
`sim.initialize()` unless marked otherwise.

---

## 1. Replacing a physics updater

Every `Field` has an `.updater` attribute. Setting it to a callable replaces
the default physics function for that field. Setting it to `None` freezes the
field at its current value.

```python
def my_temperature(sim):
    """Power-law midplane temperature profile."""
    T0  = 150.   # [K] at r0
    r0  = 1. * c.au
    q   = -0.5   # power-law index
    return T0 * (sim.grid.r / r0) ** q

sim.gas.T.updater = my_temperature
sim.initialize()
```

**Rules:**
- The function signature must be `f(sim) -> np.ndarray` with the correct shape.
- The function is called **every timestep** — keep it fast.
- To set a field once (not every step), set the value after `sim.initialize()`,
  then set `.updater = None`.
- After replacing updaters, call `sim.update()` to propagate changes to
  dependent fields before running.

---

## 2. Custom radial grid / grid refinement

Set `sim.grid.ri` (cell *interfaces*, shape `Nr+1`) before calling
`sim.initialize()`. The cell centres `sim.grid.r` are derived automatically.

```python
import numpy as np
from tripodpy import Simulation
import dustpy.constants as c

sim = Simulation()

# Start with a regular log grid, then refine around a gap location
ri_base = np.logspace(0, 3, num=100, base=10.) * c.au

def refinegrid(ri, r0, num=3):
    """Double resolution in `num` cells around r0, recursively."""
    if num == 0:
        return ri
    ind = np.argmin(r0 > ri) - 1
    indl, indr = ind - num, ind + num + 1
    ril, rir = ri[:indl], ri[indr:]
    N = (2 * num + 1) * 2
    rim = np.empty(N)
    for i in range(0, N, 2):
        j = ind - num + int(i / 2)
        rim[i]     = ri[j]
        rim[i + 1] = 0.5 * (ri[j] + ri[j + 1])
    return refinegrid(np.concatenate([ril, rim, rir]), r0, num=num-1)

ri = refinegrid(ri_base, 10. * c.au, num=3)   # refine around 10 AU

sim.grid.ri = ri   # assign before initialize
sim.initialize()
```

**Constraints:**
- the mass dimension of the dust fileds is either 2 (for ´sim.dust.Sigma´ or 5 for ´sim.dust.a´, ´sim.dust.St´ )

---

## 3. Turning off physics modules



### Turn off dust hydrodynamics (keep coagulation)
```python
# Set advective and diffusive fluxes to zero and freeze them
sim.initialize()
sim.dust.Fi.adv[:]  = 0.; sim.dust.Fi.adv.updater  = None
sim.dust.Fi.diff[:] = 0.; sim.dust.Fi.diff.updater  = None
sim.dust.Fi.tot[:]  = 0.; sim.dust.Fi.tot.updater   = None
sim.dust.S.hyd[:]   = 0.; sim.dust.S.hyd.updater    = None
```

### Turn off gas evolution (freeze gas)
```python
# Method A: set viscosity to zero (no backreaction assumed)
sim.initialize()
sim.gas.nu[:] = 0.
sim.gas.nu.updater = None

# Method B: freeze gas surface density entirely
sim.initialize()
sim.gas.Sigma.updater = None
sim.gas.S.hyd.updater = None
```

---

## 4. External sources and sinks

External source terms are added via `sim.dust.S.ext` and `sim.gas.S.ext`.
Units: g cm⁻² s⁻¹. Negative values = sinks.

```python
def dust_sink(sim):
    """Remove dust at a fixed fractional rate inside 5 AU."""
    rate = np.zeros_like(sim.dust.Sigma)
    mask = sim.grid.r < 5. * c.au
    rate[mask, :] = -sim.dust.Sigma[mask, :] / (1e3 * c.year)
    return rate

sim.dust.S.ext.updater = dust_sink
sim.initialize()
```

For gas external sources (e.g. disk wind mass loss):
```python
def disk_wind(sim):
    """Simple disk wind: lose a fraction of gas per dynamical time."""
    Omega = sim.grid.OmegaK
    return -1e-5 * sim.gas.Sigma * Omega   # [g/cm²/s]

sim.components.Default.Sigma_dot.updater = disk_wind
sim.initialize()
```

---

## 5. Planetary gaps (Kanagawa 2017 torque)

Impose a gap via a torque profile derived from the Kanagawa et al. (2017)
gap shape. This modifies gas evolution without freezing the gas.

```python
import numpy as np
from tripodpy import Simulation
import dustpy.constants as c

def kanagawa_gap(r, a_p, q, h, alpha):
    """Gap perturbation f = Sigma_g(r) / Sigma_0 (Kanagawa et al. 2017)."""
    f   = np.ones_like(r)
    K   = q**2 / (h**5 * alpha)
    Kp  = q**2 / (h**3 * alpha)
    Kp4 = Kp**0.25
    Smin = 1. / (1. + 0.04 * K)
    dr1  = (0.25 * Smin + 0.08) * Kp4
    dr2  = 0.33 * Kp4
    dist = np.abs(r - a_p) / a_p
    f    = np.where(dist < dr2, 4. / Kp4 * dist - 0.32, f)
    f    = np.where(dist < dr1, Smin, f)
    return np.clip(f, Smin, 1.)

def make_torque_updater(a_p, q, h, alpha):
    """Returns an updater function for sim.gas.torque.Lambda."""
    def torque(sim):
        r  = sim.grid.r
        nu = sim.gas.nu
        Ok = sim.grid.OmegaK
        f  = kanagawa_gap(r, a_p, q, h, alpha)
        # Lau (2024) torque prescription
        dlnf_dlnr = np.gradient(np.log(f), np.log(r))
        Lambda = -1.5 * nu * Ok * (0.5 / f - dlnf_dlnr - 0.5)
        return Lambda
    return torque

sim = Simulation()
a_p   = 10.  * c.au
q     = 30.  * c.M_earth / c.M_sun   # planet/star mass ratio
h     = 0.05                          # aspect ratio at planet
alpha = 1e-3

sim.gas.torque.Lambda.updater = make_torque_updater(a_p, q, h, alpha)
sim.initialize()
```

**Note:** The torque enters via `sim.gas.torque.Lambda` (specific angular
momentum injection rate, cm² s⁻²) and `sim.gas.torque.v` (effective velocity,
cm s⁻¹). Both are already wired into the gas evolution Jacobian in the
standard model — you only need to set the updater.

---

## 6. Ice lines / spatially varying fragmentation velocity

Change `sim.dust.v.frag` to depend on the local gas temperature.
If the temperature evolves, assign an updater; if it is static, set the
value once after `initialize()` and freeze it.

```python
def vfrag_icelines(sim):
    """Pinilla et al. (2017) ice line prescription."""
    T = sim.gas.T
    v = np.where(T < 150., 1000., 100.)   # water ice line at 150 K
    v = np.where(T <  80.,  700., v)      # NH3 at 80 K
    v = np.where(T <  44.,  100., v)      # CO2 at 44 K
    return v   # [cm/s]

sim = Simulation()
sim.initialize()

# Static temperature profile → set once, freeze updater
sim.dust.v.frag[:] = vfrag_icelines(sim)
sim.dust.v.frag.updater = None
sim.update()
```

If you have a custom evolving temperature, assign as an updater instead:
```python
sim.dust.v.frag.updater = vfrag_icelines
```

---

## 7. Planetesimal formation sink

Remove dust mass that crosses a streaming instability threshold by adding
a negative external source term and tracking the converted mass.

```python
# Threshold: eps > eps_crit triggers planetesimal formation
eps_crit = 1.0

def planetesimal_sink(sim):
    """Convert dust to planetesimals where eps > eps_crit."""
    sink = np.zeros_like(sim.dust.Sigma)
    eps  = sim.dust.eps   # shape (Nr,)
    # Apply sink in cells over threshold, distributed over mass bins
    over = eps > eps_crit
    if over.any():
        tau_form = 1. / sim.grid.OmegaK   # one orbital period timescale [s]
        # Remove dust proportionally across all mass bins
        sink[over, :] = -sim.dust.Sigma[over, :] / tau_form[over, np.newaxis]
    return sink

sim.dust.S.ext.updater = planetesimal_sink
sim.initialize()
```

Track total planetesimal mass by adding a diastole (see §8).

---

## 8. Systoles and diastoles

- **Systole**: called at every integration timestep, before the integration.
- **Diastole**: called at every snapshot output.

```python
# Track total dust mass at every output (diastole)
M_dust_history = []

def record_dust_mass(sim):
    r   = sim.grid.r       # [cm]
    ri  = sim.grid.ri
    dr  = ri[1:] - ri[:-1]
    A   = 2. * np.pi * r * dr   # annulus area [cm²]
    Sigma_d_tot = sim.dust.Sigma.sum(axis=1)   # [g/cm²]
    M_dust_history.append((float(sim.t), float((Sigma_d_tot * A).sum())))

sim.addDiastole(record_dust_mass)
sim.initialize()
```

```python
# Systole example: update a time-dependent external parameter every step
def update_planet_mass(sim):
    # Planet grows at a fixed rate
    sim._planet_mass = 1e-5 * c.M_sun * (sim.t / (1e4 * c.year))

sim.addSystole(update_planet_mass)
```

---

## 9. Boundary conditions

Default: inner = constant gradient, outer = floor value.

```python
# Inner boundary: constant value (zero-flux / closed inner wall)
from simframe.integration import Scheme
sim.initialize()
sim.gas.boundary.inner = None   # sets to zero-gradient (no-flux)

# Outer boundary: set to a fixed value
sim.gas.boundary.outer.value = sim.gas.SigmaFloor
```

For dust boundaries, the same pattern applies via
`sim.dust.boundary.inner` and `sim.dust.boundary.outer`.

---

## 10. Snapshot customisation

```python
import numpy as np
import dustpy.constants as c

sim = Simulation()

# Run for 1 Myr with 5 snapshots per decade starting at 1 kyr
sim.t.snapshots = np.hstack([
    sim.t,
    np.geomspace(1e3, 1e6, num=16) * c.year   # [s]
])

# Change output directory
sim.writer.datadir  = "runs/my_run"
sim.writer.overwrite = False   # raise if directory exists — safe default

sim.initialize()
```

To write only specific fields to the HDF5 output (reduces file size):
```python
# After initialize(), prune the writer's field list
# (advanced — consult simframe docs for Writer.toc customisation)
```
