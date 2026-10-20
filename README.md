# Lewis Dot Structure Generator (TI-84 Plus CE)

This repository is dedicated to the `lewis-dot` CEdev program for TI-84 Plus CE calculators.

`lewis-dot` is an interactive Lewis dot structure generator with charge handling and resonance navigation.

Only `lewis-dot` project files are tracked in Git for this repo. This README documents that project only.

## Screenshot

![Periodic Table Menu](lewis-dot/demoImages/periodicTableMenu.png)

## Requirements

- [CEdev](https://ce-programming.github.io/toolchain/) installed
- `CEDEV` environment variable set to your CEdev install path
- `make` available in your shell

## Build

Build from the `lewis-dot` directory:

```powershell
cd lewis-dot
make
```

Output files are generated in `lewis-dot/bin/`.

To clean generated files:

```powershell
make clean
```

## Host Tests

Run the engine tests from repository root:

```powershell
./lewis-dot/tests/run_tests.ps1
```

The script compiles and runs a host executable with `clang`, `gcc`, or `zig cc`.

## Controls

- Arrow keys: move periodic-table cursor
- `Enter`: add highlighted element
- `Del`: remove last atom
- `Alpha`: cycle charge (`0`, `+1`, `+2`, `-1`, `-2`)
- `2nd`: generate Lewis structure
- Left/Right: cycle resonance structures
- `Clear`: return to periodic table view
- `Mode`: quit

## File Purpose Map

`lewis-dot/Makefile`
- CEdev build configuration for the `LEWIS` target.

`lewis-dot/src/main.c`
- App entry point and runtime loop.
- Coordinates screen mode switching, key handling, and warning overlays.

`lewis-dot/src/lewis_model.h`
- Shared constants and core data structures (`Element`, `Molecule`, `LewisStructure`, `InvalidReason`).

`lewis-dot/src/lewis_model.c`
- Element table definitions and periodic table grid initialization.
- Model reset helper (`molecule_reset`).

`lewis-dot/src/lewis_engine.h`
- Public API for structure generation and invalid-reason messaging.

`lewis-dot/src/lewis_engine.c`
- Lewis generation logic:
- skeleton search (backbone chain, terminal-atom and hydrogen placement)
- exhaustive pi-bond placement for each skeleton, ranked octet-first, then by formal charge
- resonance forms: every best-ranked form of the chosen skeleton
- central-atom choice for display and VSEPR
- invalid-reason classification

## Structure Conventions

- Octet first: atoms get exactly an octet when possible. Period 3+ atoms exceed it only when unavoidable (`SF6`, `PCl5`, `XeF4`, `I3-`), never just to cancel formal charges. `SO4^2-` is drawn with four single bonds (S +2), and `SO2` has two resonance forms.
- Be/B/Al may keep incomplete octets (`BF3`); carbon may be a sextet only when no octet structure exists (carbocations such as the allyl cation).
- A formula with several isomers produces one of them. For example, `C2H6O` gives ethanol and `CH3NO2` gives carbamic acid (all formal charges zero), not nitromethane.
- Rings are not supported.

`lewis-dot/src/layout.h`
- Public API for atom coordinate layout helpers.

`lewis-dot/src/layout.c`
- Connectivity-aware atom coordinate placement:
- linear-chain layout for path-like graphs
- tree-from-central layout fallback

`lewis-dot/src/ui_text.h`
- Shared UI text helper declarations.

`lewis-dot/src/ui_text.c`
- Safe clipped text drawing, integer/string append helpers, and text contrast helper.

`lewis-dot/src/ui_periodic.h`
- Public API for periodic-table cursor movement and rendering.

`lewis-dot/src/ui_periodic.c`
- Periodic-table screen rendering (header, selected atoms bar, element card, help text).
- Sparse-grid cursor movement behavior.

`lewis-dot/tests/lewis_engine_tests.c`
- Host-side deterministic engine tests for representative molecules, resonance behavior, and invalid-input paths.

`lewis-dot/tests/run_tests.ps1`
- PowerShell script to compile and run host tests (`clang`, `gcc`, or `zig cc`).

`lewis-dot/tests/README.md`
- Test scope and usage instructions.
