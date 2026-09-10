# operon.core

Shared types every other `operon.*` subpackage reads and populates:
`OperonHost`, `CDS`, `Operon`, `Assembled`, `Junction`, `StandbyParams`.
**Status: implemented** (as plain dataclasses — no behavior). Spec
section 2.1.

## Rules for this package specifically

- **These types carry no behavior.** No methods beyond what a
  `@dataclass` gives for free. Logic belongs in the module that owns it
  (`operon.assembly.assemble`, `operon.coupling.coupled_tirs`, ...), not
  here. If you're tempted to add a method, it probably belongs in the
  calling module instead.
- **`OperonHost.pack` wraps an `rbsforge.HostPack`; it does not extend or
  subclass it.** Fields that belong to translation physics itself go in
  `HostPack` (`rbsforge/hostpack.py`); fields that are Operon-level and
  organism-specific but not translation-physics (codon tables, IS/att
  motifs, forbidden RE sites, `r_elong_nominal_nt_s`, coupling's
  `c_unfold`/`k_p`) go in `OperonHost`. See spec Appendix F for the exact
  split — when adding a field, check that table before deciding which
  type it goes on.
- **Coordinate convention is load-bearing across every scanner:**
  assembled DNA is 5' to 3' on the coding strand; mRNA is the same
  sequence with T to U from the promoter's TSS to the terminator's 3'
  end; every scanner reports `[start, end)`, half-open, on that molecule,
  plus a `strand` field. `operon.promoter_calculator` already follows
  this (see its `PromoterHit`); any new scanner should match its span
  convention, not invent its own.
- **`OperonHost.c_unfold` and `k_p` default to Tian 2015's published
  values (0.81, 10) as a documented prior, not a validated constant.**
  Anything that reads them (`operon.coupling`) must re-fit before
  trusting an absolute number — see `operon/coupling/CLAUDE.md`.
- Don't add a `predict`/`assemble`/`scan_*` free function to this
  package. If it doesn't belong to any specific module yet, it's not
  ready to be "shared" — leave it where it's actually used until a
  second caller needs it too.
