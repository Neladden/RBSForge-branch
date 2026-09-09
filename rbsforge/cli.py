"""Command-line interface: mRNA sequence + species + temperature -> TIR table."""
from __future__ import annotations

import argparse
import json
import sys

from .hostpack import BUILTIN_HOSTPACKS
from .rbs.calculator import RBSCalculator


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rbsforge",
        description=(
            "Thermodynamics-based ribosome binding site (RBS) strength "
            "prediction, reimplementing the Salis Lab RBS Calculator's "
            "free-energy model for arbitrary organisms."
        ),
    )
    parser.add_argument("sequence", help="mRNA (or DNA) sequence, 5' UTR through some of the CDS")
    parser.add_argument(
        "--species",
        "-s",
        default="ecoli",
        help=f"HostPack name. Built in: {', '.join(sorted(BUILTIN_HOSTPACKS))}",
    )
    parser.add_argument(
        "--temperature",
        "-t",
        type=float,
        default=None,
        help="Expression temperature in Celsius (default: the HostPack's own t_growth_c)",
    )
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON instead of text")
    return parser


def _result_to_dict(r) -> dict:
    return {
        "index": r.index,
        "codon": r.codon,
        "leaderless": r.leaderless,
        "delta_g_total": r.delta_g_total,
        "translation_initiation_rate": r.translation_initiation_rate,
        "v1_style_rate": r.v1_style_rate,
        "aligned_spacing": r.aligned_spacing,
        "raw_spacing": r.raw_spacing,
        "terms": dict(r.breakdown.terms) if r.breakdown else None,
    }


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    calculator = RBSCalculator.for_species(args.species, temperature_c=args.temperature)
    result = calculator.predict(args.sequence)

    if args.json:
        payload = {
            "species": calculator.hostpack.name,
            "temperature_c": calculator.temperature_c,
            "results": [_result_to_dict(r) for r in result.results],
        }
        print(json.dumps(payload, indent=2))
        return 0

    print(f"Species:     {calculator.hostpack.name}")
    print(f"Temperature: {calculator.temperature_c:.1f} C")
    print(f"Sequence:    {len(result.sequence)} nt")
    print()

    ranked = result.ranked()
    if not ranked:
        print("No scorable start codons found (all candidates were leaderless).")
    for r in ranked:
        marker = "*" if r is result.best else " "
        print(
            f"{marker} pos {r.index:>5}  {r.codon}  "
            f"delta_G_total = {r.delta_g_total:+7.2f} kcal/mol  "
            f"TIR = {r.translation_initiation_rate:10.1f} (proportional)"
        )
        for name, value in r.breakdown.terms.items():
            note = r.breakdown.notes.get(name, "")
            print(f"      {name:<10} {value:+7.2f}   {note}")
        print()

    skipped = [r for r in result.results if r.leaderless]
    if skipped:
        print(f"{len(skipped)} start codon(s) skipped as leaderless (no SD-like sequence found):")
        for r in skipped:
            print(f"  pos {r.index:>5}  {r.codon}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
