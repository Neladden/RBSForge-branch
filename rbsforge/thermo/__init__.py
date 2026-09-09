"""Nearest-neighbor RNA thermodynamics used by the RBS calculator.

A from-scratch, dependency-free reimplementation of the nearest-neighbor
(NN) free-energy model used to score RNA:RNA duplex hybridization
(mRNA:anti-SD binding, `duplex.py`) and single-strand mRNA secondary
structure (`fold.py`). If ViennaRNA's `RNAfold` binary is present on
PATH, `fold.best_available_folder()` uses it as a higher-accuracy
drop-in backend; everything also works with the built-in model alone.
"""
