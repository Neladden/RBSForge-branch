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
rbsforge/              Module: TIR engine (Predict mode), predict_one + footprint — implemented
operon/
  core/                shared types (OperonHost, CDS, Operon, Assembled) — implemented
  assembly/            Module A — operon assembly — implemented
  coupling/            Module B — translational coupling — implemented
  elongation/          Module C — TER + codon recoding — implemented (placeholder codon table)
  stability/           Module D — mRNA stability — features implemented, rate gated on real coefficients
  promoter_calculator/ Module E — sigma70 promoter scan — implemented
  htisc/               Module F — highly translated internal start codons — implemented
  pauses/              Module G — ribosomal pause sites — implemented (3 of 4 signals)
  terminators/         Module H — intrinsic + rho-dependent terminators — implemented (presence layer)
  repeats/             Module I — repeats, IS/att sites — implemented (seed-length, not extended)
  synthesis/           Module J — synthesis complexity, RE sites — implemented (numeric hard rules)
  design/              design_rbs + NSGA-II operon design — implemented
docs/                  architecture overview + full spec
examples/
tests/
```

Every module above is implemented with real tests, but "implemented"
does not mean "spec-complete" — most have a documented scope boundary
(a placeholder data table, a deferred optional signal, a presence-only
layer standing in for a trained classifier). Read the subpackage's own
`CLAUDE.md` before assuming a number from it is production-grade; that's
exactly what each one's "What's not implemented" / "Load-bearing
constraints" section is for.
