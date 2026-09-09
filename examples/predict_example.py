"""Example: predict translation initiation rate across a few HostPacks."""
from rbsforge import predict

SEQUENCE = "AGGAGGACAACTAAATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT"

if __name__ == "__main__":
    for species, temperature_c in [
        ("ecoli", 37.0),
        ("b_subtilis", 37.0),
        ("geobacillus", 60.0),
        ("thermus", 70.0),
    ]:
        result = predict(SEQUENCE, species=species, temperature_c=temperature_c)
        best = result.best
        print(f"=== {species} @ {temperature_c} C ===")
        if best is None:
            print("  no scorable start codon found (leaderless)")
            continue
        print(f"  best start: {best.codon} at index {best.index}")
        print(f"  delta_G_total = {best.delta_g_total:.2f} kcal/mol")
        print(f"  TIR (proportional) = {best.translation_initiation_rate:.1f}")
        print(best.breakdown)
        print()
