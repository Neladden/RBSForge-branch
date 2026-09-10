"""
Operon Calculator -- Module A: operon assembly (spec section 3).

The only place sequence is concatenated. Independent of every biophysical
model. Layout::

    [promoter][5' UTR / RBS_1][CDS_1][intergenic_1][RBS_2][CDS_2] ... [terminator]

``Operon.intergenic_policy`` governs every junction after the first
(uniformly, for the whole operon -- a documented simplification; a
per-junction policy would need a richer ``Operon`` type):

- ``"overlap-N"`` (e.g. ``"overlap-4"`` for the canonical AUGA junction,
  ``"overlap-25"`` for insulation): CDS_i is spliced onto the last N nt
  already placed -- those N nt are *replaced* by CDS_i's own leading N nt
  (CDS_i's provided sequence is authoritative for the shared span; making
  that span consistent with the upstream reading frame, e.g. a joint
  N-terminal/C-terminal recode, is a Design-mode step -- spec section
  18.3 -- not this module's job). Produces ``d_nt == -N``.
- ``"abut"``: CDS_i is placed immediately after CDS_{i-1}, no RBS content
  in between. ``d_nt == 0``. ``rbs_list[i]`` is ignored for this junction.
- ``"spacer-N"``: ``rbs_list[i]`` is inserted as literal new intergenic
  sequence and must be exactly N nt. ``d_nt == N``.
- ``"free"`` (default): ``rbs_list[i]`` (any length, including empty) is
  inserted as literal new sequence before CDS_i with no overlap trimming.
  ``d_nt`` falls out of the resulting positions.

``rbs_list[0]`` is always literal new content (the 5' UTR / RBS_1),
regardless of ``intergenic_policy``.

Coordinates: ``Assembled.features`` spans are on the full assembled DNA
(0-based from the start of ``promoter_dna``) -- promoter and terminator
only exist there. ``Assembled.starts``/``cds_end`` are on
``Assembled.mrna`` instead (0-based from the TSS), matching
``CDS.annotated_start``'s documented convention and what
``operon.coupling`` (spec section 7.3) reads.
"""

from typing import List, Tuple

from operon.core import Assembled, Junction, Operon


def _parse_policy(policy: str) -> Tuple[str, int]:
    if policy == "free":
        return "free", 0
    if policy == "abut":
        return "abut", 0
    if policy.startswith("overlap-"):
        return "overlap", int(policy[len("overlap-") :])
    if policy.startswith("spacer-"):
        return "spacer", int(policy[len("spacer-") :])
    raise ValueError(
        f"Unknown intergenic_policy {policy!r}; expected 'free', 'abut', "
        "'overlap-N', or 'spacer-N'"
    )


def assemble(operon: Operon) -> Assembled:
    """Concatenate promoter + [RBS_i + CDS_i] + terminator into one
    assembled DNA/mRNA molecule, computing each junction's intergenic
    distance ``d`` (section 3.1) along the way.

    ``d = start_{i+1} - stop_i_end``; negative means overlap.
    """
    if not operon.cds_list:
        raise ValueError("Operon.cds_list must contain at least one CDS")

    kind, n = _parse_policy(operon.intergenic_policy)

    dna_parts: List[str] = [operon.promoter_dna]
    features: List[dict] = [
        {"id": "promoter", "type": "promoter", "start": 0, "end": len(operon.promoter_dna), "strand": 1}
    ]
    junctions: List[Junction] = []
    starts_dna: List[int] = []
    cds_ends_dna: List[int] = []

    pos = len(operon.promoter_dna)

    for i, cds in enumerate(operon.cds_list):
        rbs = operon.rbs_list[i] if i < len(operon.rbs_list) else ""
        cds_nt = cds.aa_or_nt

        if i == 0:
            dna_parts.append(rbs)
            rbs_start = pos
            pos += len(rbs)
            if rbs:
                features.append(
                    {"id": f"RBS_{cds.id}", "type": "rbs", "start": rbs_start, "end": pos, "strand": 1}
                )
            start_i = pos
            d = None  # no junction before the first CDS
        else:
            prev_stop_end = cds_ends_dna[-1]
            if kind == "free":
                dna_parts.append(rbs)
                rbs_start = pos
                pos += len(rbs)
                if rbs:
                    features.append(
                        {"id": f"RBS_{cds.id}", "type": "rbs", "start": rbs_start, "end": pos, "strand": 1}
                    )
                start_i = pos
            elif kind == "abut":
                start_i = pos
            elif kind == "spacer":
                if len(rbs) != n:
                    raise ValueError(
                        f"intergenic_policy 'spacer-{n}' requires rbs_list[{i}] to be exactly "
                        f"{n} nt, got {len(rbs)}"
                    )
                dna_parts.append(rbs)
                rbs_start = pos
                pos += len(rbs)
                features.append(
                    {"id": f"intergenic_{i}", "type": "intergenic", "start": rbs_start, "end": pos, "strand": 1}
                )
                start_i = pos
            elif kind == "overlap":
                if n <= 0:
                    raise ValueError("intergenic_policy 'overlap-N' needs a positive N")
                if n > len(cds_nt):
                    raise ValueError(f"overlap-{n} is longer than CDS {cds.id!r}'s own sequence")
                combined = "".join(dna_parts)
                if n > len(combined) - len(operon.promoter_dna):
                    raise ValueError(f"overlap-{n} exceeds the sequence available upstream of CDS {cds.id!r}")
                trimmed = combined[:-n]
                dna_parts = [trimmed]
                pos = len(trimmed)
                start_i = pos
            else:  # pragma: no cover -- _parse_policy already restricts `kind`
                raise AssertionError(kind)

            d = start_i - prev_stop_end
            junctions.append(
                Junction(
                    between=(operon.cds_list[i - 1].id, cds.id),
                    stop_i=prev_stop_end,
                    start_i_plus_1=start_i,
                    d_nt=d,
                )
            )

        dna_parts.append(cds_nt)
        pos += len(cds_nt)
        stop_i_end = pos
        features.append({"id": cds.id, "type": "CDS", "start": start_i, "end": stop_i_end, "strand": 1})
        starts_dna.append(start_i)
        cds_ends_dna.append(stop_i_end)

    terminator_start = pos
    dna_parts.append(operon.terminator_dna)
    pos += len(operon.terminator_dna)
    if operon.terminator_dna:
        features.append(
            {"id": "terminator", "type": "terminator", "start": terminator_start, "end": pos, "strand": 1}
        )

    dna = "".join(dna_parts)
    tss = operon.tss if operon.tss is not None else len(operon.promoter_dna)
    mrna = dna[tss:].upper().replace("T", "U")

    return Assembled(
        dna=dna,
        mrna=mrna,
        features=features,
        junctions=junctions,
        starts=[s - tss for s in starts_dna],
        cds_end=[e - tss for e in cds_ends_dna],
    )
