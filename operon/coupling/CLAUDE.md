# operon.coupling — Module B

Translational coupling: the only translation physics not already in
RBSForge. **Status: implemented.** Spec section 7, source Tian & Salis,
NAR 2015 (https://doi.org/10.1093/nar/gkv635).

## Load-bearing constraints

- **`require_vienna=True` is the default and raises `RuntimeError` if
  `calc.folder.name != "vienna"`.** The builtin folder cannot represent
  bulge-containing intergenic hairpins; running coupling on it under-
  stabilizes inhibitory helices and flattens the sigmoid (known failure
  mode 11). `require_vienna=False` only warns and proceeds — that path
  exists for tests/development in environments without `RNAfold` on
  `PATH` (this repo's own test suite uses it that way, since the sandbox
  this was built in has no Vienna installed); never use it to produce
  numbers anyone will actually trust.
- `to_physical(rate_au, host)` **is the identity function, deliberately.**
  The spec is explicit that the au->s^-1 conversion is absorbed into a
  properly re-fit `host.c_unfold`, not a separately-published formula —
  inventing one here would silently double-convert. Don't add a
  conversion formula to this function without a cited source.
- `host.c_unfold` (`C`) and `host.k_p` still default to Tian's published
  values (0.81, 10) as an unfit prior — re-fit both before trusting an
  absolute coupled TIR (see `operon/core/CLAUDE.md`). `k_reinitiation(d)`
  is the one piece that transfers across au scales as-is; don't touch it
  without a reason from the spec.
- `extra_unpaired_global` passed to `predict_one` for the "unfolded"
  score is **the upstream CDS's entire span**
  (`range(assembled.starts[i-1], assembled.cds_end[i-1])`), not just a
  window near its stop — `predict_one` clips it to whatever falls inside
  the downstream start's own folded window itself (spec Appendix C: "every
  index < start_i that is part of CDS_{i-1}... clipped to the window").
  Don't pre-clip it here; that duplicates logic `predict_one` already owns.
- Leaderless downstream start -> de novo term is exactly `0`, keeping
  only `r_reinit` (`LeaderlessCouplingError` is only raised for a
  leaderless *first* CDS, where there is nothing to couple from at all —
  a downstream leaderless start is a normal, real architecture, not an
  error).
- There is no closed-form inverse of this module — `operon.design` must
  sample-and-score across junctions, never invert `coupled_tirs` (spec
  section 7.5).

## Tests

`tests/test_coupling.py`: `k_reinitiation`'s published fixed points and
interpolation shape, the Vienna-required refuse/warn behavior, first-CDS
mono-cistronic equivalence, `d=-4` vs. insulating `d=-25` producing a
measurably larger re-initiation leak, both leaderless cases (first CDS
raises, downstream CDS keeps only re-init), and a trivial single-CDS
operon.
