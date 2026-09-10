"""Example: scan a DNA sequence for its strongest sigma70 promoter and
list any other (cryptic) promoters found nearby."""
from operon.promoter_calculator import predict, scan_promoters

# A consensus sigma70 promoter: -35 (TTGACA), a 17 bp spacer,
# -10 (TATAAT), a 6 nt discriminator, and a 20 nt ITR.
SEQUENCE = "A" * 30 + "TTGACA" + "A" * 17 + "TATAAT" + "C" * 6 + "G" * 20 + "A" * 10

if __name__ == "__main__":
    result = predict(SEQUENCE, organism="ecoli")
    best = result.best
    print(f"strongest TSS: {best.tss} (strand {best.strand:+d})")
    print(f"  Tx_rate = {best.tx_rate:.1f}")
    print(f"  dG_total = {best.dg_total:.3f}")
    print(f"  -35 = {best.hex35}  spacer = {len(best.spacer)} nt  -10 = {best.hex10}")
    print(f"  {best.breakdown.terms()}")
    print()

    cryptic = scan_promoters(SEQUENCE, intended_tss=best.tss, intended_strand=best.strand, tau_tx=1.0)
    print(f"{len(cryptic)} cryptic promoter(s) above tau_tx=1.0:")
    for hit in cryptic[:5]:
        print(f"  TSS {hit.tss} (strand {hit.strand:+d}): Tx_rate={hit.tx_rate:.2f}")
