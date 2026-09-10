# operon.pauses — Module G

Ribosomal pause site scanning inside CDSs. **Status: implemented** for
three of the four documented signals (internal SD/anti-SD, slow-codon
runs, polyproline/stall motifs). Spec section 12. Independent motif +
structure scanner; no standalone Salis "Pause Calculator" was ever
published, so this is an operational definition assembled from adjacent
tools, not a single model to port.

## What's not implemented

Stable mRNA hairpins in the ribosomal E/P site window are not scanned —
the spec itself flags this signal as "optional, weaker evidence for
bacteria." Don't assume `scan_pauses` reports hairpin-based pauses.

## A bug already caught here, worth knowing about

`best_hybridization` can return a duplex that satisfies the minimum
alignment length but is thermodynamically **unfavorable** (positive
`delta_g_hybrid` — a weak, coincidental sequence match, not something
that would actually form). An early version of `_scan_internal_sd`
reported these as pause hits; a clean all-preferred-codon test CDS
(`GCGGCGGCGGCG`, no real SD-like structure) still triggered false
positives this way. `_INTERNAL_SD_MAX_DELTA_G = -1.0` gates on top of
the length check now — if you touch `_scan_internal_sd`, keep both
conditions (`duplex.length >= _INTERNAL_SD_MIN_DUPLEX_LENGTH and
duplex.delta_g_hybrid <= _INTERNAL_SD_MAX_DELTA_G`), not just one.

## Other things to know

## What "pause site" means here (union of four signals)

1. **Internal SD / anti-SD hybridization during elongation** — hexamers
   complementary to the 16S tail, in-frame or out-of-frame, inside a CDS.
   Reuse `rbsforge.thermo.duplex.best_hybridization` on a sliding window
   — it's already a sliding anti-SD search, just not wrapped for this use
   (spec section 4.6 lists this as an "already works" building block).
   **Use `host.pack.asd`, never a hardcoded *E. coli* tail** — this is
   the one place a naive port would silently assume *E. coli* and be
   wrong for every other HostPack.
2. Slow codon runs (consecutive low-`w` codons; a run of >=2-3 counts).
3. Polyproline/stall motifs (PPP, PPG at minimum; SecM/TnaC-class arrest
   peptides are an optional expansion).
4. Stable mRNA hairpins in the ribosomal E/P site window — optional,
   weaker evidence for bacteria than the internal-SD signal; don't weight
   it as heavily as (1).

## Gotchas

- Score internal-SD hits weighted by predicted duplex ΔG, not as a
  binary presence/absence count — a marginal hexamer match and a
  near-perfect anti-SD duplex are not the same pause strength.
- The 2016 SEED abstract explicitly allows keeping a pause on purpose
  (e.g. to allow co-translational folding) — implement any per-CDS
  whitelist/exclusion as a first-class option, not an afterthought filed
  under "TODO."
- Design-mode kill (synonym swap that breaks the SD-like hexamer or
  replaces the slow codon) belongs in `operon.elongation`'s recoding
  mutators, not here — this module only detects and scores.

## Interface to build toward

```python
def scan_pauses(cds_nt: str, host, calc) -> list: ...
```
