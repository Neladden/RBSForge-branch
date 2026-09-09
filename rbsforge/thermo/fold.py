"""Single-strand mRNA secondary structure (minimum free energy) folding.

Provides delta_G_mRNA (the MFE of a window of mRNA sequence) and, via the
`forced_unpaired` option, the constrained folds used by the v1.0-style
standby-site penalty (`rbs/standby.py`): the MFE of the same window with a
given set of positions forbidden from pairing.

Two backends:

* `BuiltinFolder` -- a from-scratch, dependency-free simplified
  Zuker/Nussinov-style dynamic program, temperature-scaled through
  `thermo.nn_params`. To keep the implementation compact and low-risk, it
  models only perfectly nested helices (contiguous Watson-Crick/G-U
  stacks) and hairpin loops; it does not model internal loops, bulges, or
  explicit multiloop initiation penalties. This makes it a conservative
  (somewhat stability-underestimating) approximation, not a substitute
  for a full nearest-neighbor folding engine -- but it requires no
  external dependencies and is adequate for comparative scoring of short
  (<150 nt) mRNA windows.
* `ViennaFolder` -- a thin wrapper around the `RNAfold` binary from the
  ViennaRNA package (the model's reference engine from v1.1 onward), used
  automatically when available on PATH for higher-accuracy, full-featured
  structure prediction at the requested temperature. Never required.

`best_available_folder()` picks ViennaRNA when present, else the builtin.
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from typing import FrozenSet, List, Tuple

from ..constants import REFERENCE_TEMPERATURE_K
from ..seqtools import is_watson_crick, is_wobble, validate_rna
from .nn_params import hairpin_loop_energy, watson_crick_stack

_INF = float("inf")


@dataclass
class FoldResult:
    delta_g: float
    pairs: List[Tuple[int, int]] = field(default_factory=list)
    backend: str = "builtin"

    def to_dot_bracket(self, length: int) -> str:
        chars = ["."] * length
        for i, j in self.pairs:
            chars[i] = "("
            chars[j] = ")"
        return "".join(chars)


def _pairable(a: str, b: str) -> bool:
    return is_watson_crick(a, b) or is_wobble(a, b)


class BuiltinFolder:
    """Zuker/Nussinov-lite MFE folder (stacks + hairpin loops only)."""

    name = "builtin"

    def __init__(self, temperature_k: float = REFERENCE_TEMPERATURE_K):
        self.temperature_k = temperature_k

    def fold(self, sequence: str, forced_unpaired: FrozenSet[int] = frozenset()) -> FoldResult:
        """MFE fold of `sequence`, with positions in `forced_unpaired`
        (0-based) forbidden from participating in any base pair -- used to
        compute the v1.0 standby-site penalty (the energetic cost of
        forcing the standby-site nucleotides to stay single-stranded)."""
        seq = validate_rna(sequence)
        n = len(seq)
        if n < 4:
            return FoldResult(delta_g=0.0, pairs=[], backend=self.name)

        v = [[_INF] * n for _ in range(n)]  # V[i][j]: best energy with i,j paired
        v_choice = [[None] * n for _ in range(n)]  # "hairpin" or "stack"

        for span in range(4, n):  # span = j - i; minimum loop needs j-i>=4
            for i in range(0, n - span):
                j = i + span
                if i in forced_unpaired or j in forced_unpaired:
                    continue
                if not _pairable(seq[i], seq[j]):
                    continue
                best_dg, best_kind = _INF, None
                loop_size = j - i - 1
                if loop_size >= 3:
                    dg = hairpin_loop_energy(loop_size, self.temperature_k)
                    if dg < best_dg:
                        best_dg, best_kind = dg, "hairpin"
                if j - 1 > i + 1 and v[i + 1][j - 1] < _INF:
                    key = f"{seq[i]}{seq[i + 1]}/{seq[j - 1]}{seq[j]}"
                    dg = watson_crick_stack(key, self.temperature_k) + v[i + 1][j - 1]
                    if dg < best_dg:
                        best_dg, best_kind = dg, "stack"
                if best_kind is not None:
                    v[i][j] = best_dg
                    v_choice[i][j] = best_kind

        w = [0.0] * (n + 1)  # w[j+1] == W(j); w[0] == W(-1) == 0
        traceback: List[Tuple[str, int, int]] = [None] * (n + 1)  # type: ignore
        for j in range(0, n):
            best = w[j]  # j unpaired: W(j-1)
            best_choice = ("unpaired", j, j)
            for i in range(0, j - 3):
                if v[i][j] < _INF:
                    candidate = w[i] + v[i][j]
                    if candidate < best:
                        best = candidate
                        best_choice = ("paired", i, j)
            w[j + 1] = best
            traceback[j + 1] = best_choice

        pairs = self._traceback(traceback, v_choice, n)
        return FoldResult(delta_g=w[n], pairs=pairs, backend=self.name)

    def _traceback(self, traceback, v_choice, n) -> List[Tuple[int, int]]:
        pairs: List[Tuple[int, int]] = []
        j_idx = n
        while j_idx > 0:
            kind, i, j = traceback[j_idx]
            if kind == "unpaired":
                j_idx = j  # W(j) reduced to W(j-1); j here equals old j (== j_idx-1)
            else:
                self._traceback_v(i, j, v_choice, pairs)
                j_idx = i
        return sorted(pairs)

    def _traceback_v(self, i, j, v_choice, pairs) -> None:
        while True:
            pairs.append((i, j))
            if v_choice[i][j] == "stack":
                i, j = i + 1, j - 1
            else:  # hairpin leaf (or, defensively, an unset choice)
                break


class ViennaFolder:
    """Wraps the external `RNAfold` binary (ViennaRNA), if installed."""

    name = "vienna"

    def __init__(self, temperature_k: float = REFERENCE_TEMPERATURE_K):
        self.temperature_k = temperature_k

    @staticmethod
    def available() -> bool:
        return shutil.which("RNAfold") is not None

    def fold(self, sequence: str, forced_unpaired: FrozenSet[int] = frozenset()) -> FoldResult:
        seq = validate_rna(sequence)
        temp_c = self.temperature_k - 273.15
        args = ["RNAfold", "--noPS", "-T", f"{temp_c:.2f}"]
        stdin = seq + "\n"
        if forced_unpaired:
            constraint = "".join("x" if i in forced_unpaired else "." for i in range(len(seq)))
            args.append("-C")
            stdin = seq + "\n" + constraint + "\n"
        proc = subprocess.run(
            args, input=stdin, capture_output=True, text=True, timeout=10, check=True
        )
        lines = proc.stdout.strip().splitlines()
        structure_line = lines[1]
        dot_bracket, energy_str = structure_line.rsplit("(", 1)
        dot_bracket = dot_bracket.strip()
        delta_g = float(energy_str.rstrip(")"))
        pairs = _pairs_from_dot_bracket(dot_bracket)
        return FoldResult(delta_g=delta_g, pairs=pairs, backend=self.name)


def _pairs_from_dot_bracket(dot_bracket: str) -> List[Tuple[int, int]]:
    stack: List[int] = []
    pairs: List[Tuple[int, int]] = []
    for idx, ch in enumerate(dot_bracket):
        if ch == "(":
            stack.append(idx)
        elif ch == ")":
            i = stack.pop()
            pairs.append((i, idx))
    return sorted(pairs)


def best_available_folder(temperature_k: float = REFERENCE_TEMPERATURE_K):
    """Return a ViennaFolder if RNAfold is on PATH, else BuiltinFolder."""
    if ViennaFolder.available():
        return ViennaFolder(temperature_k)
    return BuiltinFolder(temperature_k)
