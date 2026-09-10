"""
Operon Calculator -- a recreation of the Salis Lab Operon Calculator.

RBSForge (the ``rbsforge`` package, a sibling of this one) is the
translation-initiation-rate engine for the first module of this stack and
predates this package; it is a subset of the Operon Calculator, not the
whole of it. See ``docs/OPERON_CALCULATOR.md`` for the full architecture
and ``docs/OPERON_CALCULATOR_SPEC.md`` for the complete design spec these
subpackages implement against.

Subpackages, one per independent module of the spec:

===================================  =======  =========================================  ==========
Subpackage                           Module   What                                       Status
===================================  =======  =========================================  ==========
``rbsforge`` (sibling)                --      TIR engine (Predict mode)                  implemented
``operon.core``                       --      shared types (HostPack wrapper, Operon)    implemented
``operon.assembly``                   A       operon assembly, intergenic distance d     planned
``operon.coupling``                   B       translational coupling                     planned
``operon.elongation``                 C       TER + synonymous codon recoding            planned
``operon.stability``                  D       mRNA stability                             planned
``operon.promoter_calculator``        E       sigma70 promoter / cryptic-promoter scan   implemented
``operon.htisc``                      F       highly translated internal start codons    planned
``operon.pauses``                     G       ribosomal pause sites                      planned
``operon.terminators``                H       intrinsic + rho-dependent terminators      planned
``operon.repeats``                    I       repeats, IS/att sites                      planned
``operon.synthesis``                  J       synthesis complexity, restriction sites    planned
``operon.design``                     --      design_rbs + NSGA-II design search         planned
===================================  =======  =========================================  ==========

Module K (system-level RNAP/ribosome load) is optional bookkeeping over
the other modules' outputs and has no subpackage of its own yet.
"""
