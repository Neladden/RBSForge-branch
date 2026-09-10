# Operon Calculator

A from-scratch recreation of the Salis Lab Operon Calculator. RBSForge
(`rbsforge/`) is one module of this project — the translation-initiation-
rate engine — not the whole of it. Everything else lives under `operon/`,
one subpackage per independent module of the design spec.

**Each subproject has its own `CLAUDE.md`.** Read the one in the
directory you're actually editing before making changes — it has that
module's spec section, status, dependencies, and gotchas. This file is
only the cross-cutting index.

## Read before doing anything non-trivial

- `docs/OPERON_CALCULATOR.md` — architecture diagram, module status table,
  implementation order.
- `docs/OPERON_CALCULATOR_SPEC.md` — the full design spec every module
  implements against. This is the source of truth for algorithms,
  citations, constants, and known failure modes. When a subproject's
  `CLAUDE.md` cites a section number, it's a section in this file.

## Repo-wide conventions

- **No required third-party dependencies**, in `rbsforge` or in any
  `operon.*` subpackage. Optional dependencies (ViennaRNA/`RNAfold`,
  eventually) are opt-in fallback paths, never a hard import.
- **Don't invent unpublished coefficients.** If a paper says a term is
  "not implemented" or "unpublished," leave it at an explicit `0.0` or
  raise `NotImplementedError` — don't guess a value. RBSForge's
  `delta_G_stacking = 0.0` and the promoter calculator's dead spacer-
  categorical term are both deliberate; don't "complete" them.
- **A published constant fit on one au/energy scale doesn't transfer to
  another without re-fitting.** Tian 2015's `C` and `k_P` are the
  standing example (spec section 4.5, 23). Ratios (like
  `k_reinitiation(d)`) do transfer; absolute rates and occupancy
  thresholds usually don't. When adding a new module that consumes a
  published constant, check spec section 23 before wiring it in.
- **Don't put one subpackage's organism-specific fields into another's
  shared type.** `operon.core.OperonHost` wraps `rbsforge.HostPack`; it
  does not fork it, and `HostPack` does not grow codon tables or IS
  motifs. See `operon/core/CLAUDE.md`.
- Tests are plain `unittest`, one file per module under `tests/`, run
  with `python -m unittest discover -s tests` from the repo root (or
  `pytest`). Every module that ships real code ships tests with it, in
  the same commit.
- A "prepared subsection" (see `docs/OPERON_CALCULATOR.md`'s status
  table) is a package with a docstring and stub function signatures that
  raise `NotImplementedError` — real scaffolding, not empty. Filling one
  in means replacing the stub body and adding tests, not restructuring
  the package.

## Layout

```
rbsforge/              Module: TIR engine (Predict mode) — implemented
operon/
  core/                shared types (OperonHost, CDS, Operon, Assembled)
  assembly/            Module A — operon assembly
  coupling/            Module B — translational coupling
  elongation/          Module C — TER + codon recoding
  stability/           Module D — mRNA stability
  promoter_calculator/ Module E — sigma70 promoter scan — implemented
  htisc/               Module F — highly translated internal start codons
  pauses/              Module G — ribosomal pause sites
  terminators/         Module H — intrinsic + rho-dependent terminators
  repeats/             Module I — repeats, IS/att sites
  synthesis/           Module J — synthesis complexity, RE sites
  design/              design_rbs + NSGA-II operon design
docs/                  architecture overview + full spec
examples/
tests/
```
