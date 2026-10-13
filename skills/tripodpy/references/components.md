# TriPoDPy Compositional Tracking (`sim.addcomponent`)

TriPoDPy can split the gas and the dust into **components** (chemical species,
tracers of origin, …). Source: `tripodpy/std/addcomponent.py`, `component.py`,
`compo.py`; docs notebooks `5_compositional_tracking.ipynb` and
`Example_disk_composition.ipynb`. Everything below was run against `tripodpy 1.0.0`.

---

## 1. Basics

- `sim.components` always contains **`Default`** after `initialize()`: an active gas
  component holding all the gas (`pars.mu = ini.gas.mu`).
- `sim.gas.Sigma` = Σ of all components' `gas.Sigma`; `sim.gas.mu` = Σ / Σᵢ(Σᵢ/μᵢ).
- If any component is `dust_active`, `sim.dust.Sigma` is **overwritten each step by the sum of
  the active dust components** (floored at `dust.SigmaFloor`) and `sim.dust.rhos` becomes the
  mass-weighted harmonic mean of the components' `dust.pars.rhos`.
- Each component is integrated by its own implicit instruction in `sim.integrator`.
- Output: `components/<name>/gas/Sigma`, `.../gas/value`, `.../dust/Sigma`,
  `.../dust/value`, `.../gas/pars/mu`, … in every HDF5 file.

### Signature

```python
sim.addcomponent(name, gas_value, mu, dust_value=None,
                 dust_active=False, dust_tracer=False,
                 gas_active=False, gas_tracer=False,
                 description="", rhos=1.0)
```

- Only exists **after** `sim.initialize()` (before: `AttributeError`).
- `name` must be new (`"Default"` again → `RuntimeError: Component with name ... already exists`).
- `gas_value` `(Nr,)`: gas surface density [g/cm²] if `gas_active`, tracer value if `gas_tracer`.
  For pure dust components pass `sim.gas.SigmaFloor` (as the docs do).
- `mu` [g]: mean molecular weight (CGS! e.g. `28 * c.m_p`), used if `gas_active`.
- `dust_value` `(Nr, 2)`: dust surface density [g/cm²] if `dust_active`, tracer value if `dust_tracer`.
- `rhos` [g/cm³]: bulk density of the solid, used if `dust_active`.

### Supported mode combinations

| `gas_active` | `gas_tracer` | `dust_active` | `dust_tracer` | Meaning | Status |
|---|---|---|---|---|---|
| ✓ | | | | extra gas species (like `Default`) | works |
| | ✓ | | | passive gas tracer (`gas.value`) | works |
| | | | ✓ | passive dust tracer (`dust.value`, (Nr,2)) | works |
| | | ✓ | | dust species (contributes to `dust.Sigma`, `rhos`) | works |
| ✓ | | ✓ | | volatile in gas **and** solid, sublimation/condensation | works (see §4 pitfall) |
| ✓ | | | ✓ | gas species + its solid phase traced as `value·dust.Sigma` (exchanges via sublimation/condensation, does not change dust mass) | runs; set `pars.nu/Tsub/mu` |
| | ✓ | | ✓ | gas and dust tracer in one component | **fails**: `RuntimeError: Component type not recognized` — use two components |
| active+tracer in the same phase | | | | | `AssertionError` |

Check a component's mode with `comp.gas._active`, `comp.gas._tracer`,
`comp.dust._active`, `comp.dust._tracer`. Properties of inactive phases return zeros and
refuse assignment (`RuntimeError: Do not set Sigma for inactive ... species`).

---

## 2. Tracers

```python
import numpy as np
from tripodpy import Simulation
import tripodpy.constants as c

sim = Simulation()
sim.initialize()

tr_gas  = np.where(sim.grid.r < 10 * c.au, 1.0, 1e-5)           # (Nr,)
tr_dust = tr_gas[:, None] * np.ones_like(sim.dust.Sigma)         # (Nr, 2)
sim.addcomponent("tr_gas",  tr_gas,  1., gas_tracer=True)
sim.addcomponent("tr_dust", tr_dust, 1., dust_value=tr_dust, dust_tracer=True)

sim.writer.datadir = "data_tracers"
sim.run()

sim.components.tr_gas.gas.value     # (Nr,)   advected with the gas
sim.components.tr_dust.dust.value   # (Nr, 2) advected with each dust population
```
(For the dust tracer, `gas_value` is not used — the docs pass the same array.)

Tracer values are integrated as `value × Σ` and divided back, so they are
mass-weighted fractions in [min, max] of the initial values (verified: a uniform tracer stays
exactly uniform). Boundary conditions default to the initial edge values.

Reading tracers from disk:
```python
from tripodpy import hdf5writer
data = hdf5writer(datadir="data_tracers").read.all()
data.components.tr_gas.gas.value      # (Nt, Nr)
data.components.tr_dust.dust.value    # (Nt, Nr, 2)
```

---

## 3. Active dust components (several dust species)

Once one component is `dust_active`, **all** dust must live in active components,
because `sim.dust.Sigma` is replaced by their sum. Always add a background component:

```python
sim = Simulation()
sim.initialize()

d_out = np.array(sim.dust.Sigma); d_out[sim.grid.r <  10 * c.au] = 1e-11
d_in  = np.array(sim.dust.Sigma); d_in [sim.grid.r >= 10 * c.au] = 1e-11
rho0  = float(sim.dust.rhos[0, 0])
sim.addcomponent("dust_out", sim.gas.SigmaFloor, 1., dust_value=d_out, dust_active=True, rhos=rho0)
sim.addcomponent("dust_in",  sim.gas.SigmaFloor, 1., dust_value=d_in,  dust_active=True, rhos=rho0)
```
Verified: Σ(components) equals `sim.dust.Sigma` in every cell above the floor.

- Different `rhos` per component → `sim.dust.rhos` varies with r (same for all 5 sizes, as TriPoD requires).
  After changing `comp.dust.pars.rhos`, call `sim.dust.rhos.update()`.
- External dust sources/sinks must go on the components (`comp.dust._Sigma_dot`, shape (Nr,2));
  `sim.dust.S.ext` is ignored in this mode (see customisation.md §4).
- Per-component size distribution: `tripodpy.utils.sim_size_distribution(sim, "dust_out")`
  (without a name: total dust). Returns `(a, a_i, sig_da)` with `sig_da` shape (Nr, Na).

---

## 4. Volatiles: sublimation and condensation (`gas_active=True, dust_active=True`)

Rates follow Penzlin et al. (2024) (`compo.L_sublimation`, `compo.L_condensation`).
Required per species (defaults are placeholders: `mu = 2.3 m_p`, `nu = 1.4`, `Tsub = 10 K`):

| Parameter | Meaning |
|---|---|
| `comp.gas.pars.mu` | molecular mass [g] |
| `comp.gas.pars.nu` | attempt frequency [1/s] |
| `comp.gas.pars.Tsub` | binding energy in Kelvin [K] |
| `comp.dust.pars.rhos` | bulk density of the ice [g/cm³] |

```python
sim = Simulation()
sim.ini.grid.Nr = 100
sim.ini.gas.SigmaRc = 100 * c.au
sim.initialize()

# background (refractory) dust carrying all current dust
sim.addcomponent("default_dust", sim.gas.SigmaFloor, 1.,
                 dust_value=sim.dust.Sigma, dust_active=True, rhos=float(sim.dust.rhos[0, 0]))

# CO: start mostly in the ice, with a small but non-floor gas seed (see pitfall!)
co_ice = np.minimum(5e-3 * sim.gas.Sigma, sim.dust.Sigma.sum(-1))[:, None] * 0.5 * np.ones((1, 2))
co_gas = 1e-10 * np.asarray(sim.gas.Sigma)
sim.addcomponent("CO", co_gas, 28 * c.m_p, dust_value=co_ice,
                 dust_active=True, gas_active=True, rhos=1.)
co = sim.components.CO
co.gas.pars.mu   = 28 * c.m_p
co.gas.pars.nu   = 7.0e11
co.gas.pars.Tsub = 1180.
co.dust.pars.rhos = 1.0
sim.dust.rhos.update()
```

**Pitfall (verified, important):** in `compo.jacobian_compo` sublimation *and*
condensation are switched off in every cell where `comp.gas.Sigma < 10 * sim.gas.SigmaFloor`
(= 10⁻⁹⁹) or where `comp.dust.Sigma == sim.dust.SigmaFloor`. If the gas phase is
initialised **at** `sim.gas.SigmaFloor` (as in docs notebook 5), the species never
sublimates — CO stayed entirely solid for 10⁴ yr even at 200 K. Seeding the gas with
`1e-10 * sim.gas.Sigma` made CO sublimate inside ~120 au within 100 yr.
Starting with everything in the gas phase (as in `Example_disk_composition.ipynb`) also works.

Other pitfalls (from the docs + code):
- Truncate volatile initial profiles where the dust is at its floor
  (`np.where(sim.dust.Sigma[:,0] > sim.dust.SigmaFloor[:,0], 1., 0.)`).
- When adding gas species, subtract their mass from `sim.components.Default.gas.Sigma`
  (otherwise total gas mass increases), and set `Default.gas.pars.mu` to the remaining mixture,
  then `sim.update()`.
- Everything in CGS — especially `mu` (grams, not amu).
- `sim.gas.mu` in cells where the gas is at the floor (outer edge) can take odd values — mask them.

---

## 5. Analysing components

```python
for name, comp in sim.components.__dict__.items():
    if name.startswith("_"):
        continue                          # hidden flags live in the same dict
    if comp.gas._active:
        print(name, "gas mass [M_earth]:", float((comp.gas.Sigma * sim.grid.A).sum()) / c.M_earth)
    if comp.dust._active:
        print(name, "dust mass [M_earth]:", float((comp.dust.Sigma.sum(-1) * sim.grid.A).sum()) / c.M_earth)
# (Fields are numpy subclasses; wrap reductions in float() for clean printing)
```
From files: `hdf5writer(datadir=...).read.sequence("components.CO.gas.Sigma")` → (Nt, Nr).
