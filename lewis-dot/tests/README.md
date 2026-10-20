# Lewis Engine Tests

These tests validate deterministic Lewis generation behavior for:

- `H2`, `O2`, `N2`, `CH4`, `NH3` (basic baseline molecules)
- `CO2`
- `CSe2` (center selection and linear geometry)
- `NO2-` (two-form resonance)
- `ClO3-`, `ClO4-`, `SO4^2-`, `PO4^3-` (octet-first: single bonds, charged center, one form)
- `SO2` and `SO3` (octet-first resonance: 2 and 3 forms)
- `NO3-`
- `CO3^2-`
- `NH4+`, `H3O+`, `BH4-`, `BF4-` (four-coordinate period-2 centers)
- `BF3` (incomplete-octet acceptance)
- `SF6` (expanded-valence central atom)
- `PCl5` (expanded-valence central atom)
- `ICl5` (expanded-valence halogen central atom)
- `XeF2` and `XeF4` (hypervalent noble-gas centers)
- `IF7` (7-domain pentagonal-bipyramidal VSEPR mapping)
- resonance away from a single symmetric center: `HNO3`, `HCO3-`, `N3-`, `OCN-`, `SCN-`, `N2O`, `NO2Cl`, allyl cation `C3H5+`, oxalate `C2O4^2-`, `N2O4`
- multi-center skeletons with no charge separation: `HCOOH`, `CH3COOH`, `C2H3Cl`, `CH3CN`, `C2N2`
- skeleton placement: H on O in `HClO3`/`H2SO3`, `S2O3^2-` center, O-centered `Cl2O`
- atom entry order does not change the result (`HNO3` entered scrambled)
- VSEPR lookup mapping for `CO2`, `NO3-`, `NH4+`, `H2O`, `PCl5`, and `SF6`
- VSEPR center for chain molecules (`HNO2` reports N, not O) and no bond angle for diatomics
- VSEPR non-null handling for valid structures (including `H2`) and invalid-input guards
- no-atoms rejection
- negative-electron rejection (invalid charge)
- skeleton-build rejection (`He2`)
- odd-electron rejection (`NO`)

The suite also checks shared invariants for successful generations:

- formal charge sum equals molecular charge in every resonance form
- resonance forms are deduplicated (no duplicate bond/lone-pair states)
- expected resonance counts and central-atom bond orders for representative ions

Run from PowerShell:

```powershell
./tests/run_tests.ps1
```

The script builds a small host executable using `clang`, `gcc`, or `zig cc`.
