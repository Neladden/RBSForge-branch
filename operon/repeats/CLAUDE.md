# operon.repeats — Module I

Repeat finding (>=12 bp) plus transposon-insertion / phage att-site
motif scanning. **Status: implemented.** Spec section 14. Source:
Hossain et al. 2020 (https://doi.org/10.1038/s41587-020-0584-2,
Nonrepetitive Parts Calculator).

## What "implemented" means here — a deliberate scope boundary

`find_repeats`'s direct/inverted hits are reported at **seed (k-mer)
length, not maximally extended**. Two 12-mer-identical loci that actually
share 20 identical bp are reported as a 12 bp hit, not a 20 bp one.
Extending each seed to its true maximal match (the "extend, merge
overlapping" part of spec section 14.1) is a real follow-up, not done —
don't assume `RepeatHit.length` is the true shared-sequence length for
direct/inverted hits; it's always exactly `k`. Tandem/terminal hits
*are* reported at their true span (they come from an independent
run-length scan, not the k-mer index), and periods are capped at 6
(homopolymer through hexanucleotide) — longer-period tandem structure
falls through to the direct-repeat scan instead once its span reaches
`k`.

## `find_repeats`

Four repeat types to report separately (spec section 14.1): direct,
reverse-complement (inverted), tandem, and terminal (tandem of length
>=5 at either sequence end — synthesis + recombination risk). `L_max =
12` is the published cutoff — a conservative design rule below the ~20-25
bp where homologous recombination becomes efficient in *E. coli*, chosen
partly because it also helps `operon.synthesis`.

**Implement this as O(L), not O(L^2).** Seed with 12-mers (canonical
strand + its reverse complement as the same key), extend matches, merge
overlapping spans, keep unique pairs — a dict keyed by up-to-4^12 12-mers
is the whole trick; don't reach for an all-pairs comparison.

## `scan_motifs` (IS/att sites — spec section 14.3)

A host-specific motif table (IS element inverted-repeat cores for the
target organism; phage attB/attP cores for relevant phages; optionally
Chi sites for RecBCD-aware hosts) scanned both strands, exact or
1-mismatch. This table is genuinely host-specific — *E. coli* K-12's IS
list is not B. subtilis's or P. putida's. Take it from
`operon.core.OperonHost.is_motifs`/`att_motifs`, don't hardcode one
organism's table as a default inside this module.

## Design-mode use (this module drives a hard constraint, not just a score)

- Maintain a `used_kmers` set across the whole operon (plus an optional
  background: genome, plasmid backbone, already-accepted cistrons in a
  multi-operon design). A codon swap or RBS mutation proposed by
  `operon.design`'s GA is **illegal**, not just penalized, if it
  introduces a 12-mer already in `used_kmers` — unless it only duplicates
  within the same local tandem already being broken.
- IS/att motifs are rare enough that a synonym swap almost never trades
  off against TIR — treat a hit as cheap to fix, not as an objective to
  balance against others.

## Interface to build toward

```python
def find_repeats(dna: str, k: int = 12): ...
def scan_motifs(dna: str, motifs: list) -> list: ...
```
