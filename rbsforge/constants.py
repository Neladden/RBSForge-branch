"""Physical constants and reference values for the RBS Calculator model.

Values here reproduce the public parameter set of the Salis Lab RBS
Calculator (v1.0/v2.0-class), as reconstructed from the primary papers and
the v1.0 GPL source (Salis, Mirsky & Voigt, Nat Biotechnol 2009; see
docs/MODEL.md for full citations). All free energies are kcal/mol,
temperatures in Kelvin unless noted otherwise.
"""

# Universal gas constant, kcal/(mol*K)
GAS_CONSTANT_KCAL = 0.0019872041

REFERENCE_TEMPERATURE_K = 310.15  # 37 degrees C, the model's default/fit temperature

# --- v1.0 Boltzmann-map constants (Salis 2009) ---------------------------
#
#   r = K * exp(-delta_G_total / RT_eff)   ==   r = exp(-beta * delta_G_total)
#
# with RT_eff = 1/beta. beta is an *empirical* conversion from computed
# free energy to in-cell rate, fit almost entirely on E. coli near 37 C --
# it is NOT 1/RT and should not be replaced by 1/RT(T) without a new fit
# (see docs/MODEL.md section "Temperature and beta").
V1_PREFACTOR_K = 2500.0
V1_LOG_K = 7.824
BETA_DEFAULT = 0.45  # mol/kcal, fit +/- 0.05 (Salis 2009)
RT_EFF_DEFAULT = 1.0 / BETA_DEFAULT  # ~2.222 kcal/mol


def celsius_to_kelvin(temperature_c: float) -> float:
    return temperature_c + 273.15


def rt(temperature_k: float = REFERENCE_TEMPERATURE_K) -> float:
    """Ideal-solution RT in kcal/mol at the given temperature. Distinct
    from RT_eff / beta, which are empirical, not ideal-solution, values --
    see module docstring."""
    return GAS_CONSTANT_KCAL * temperature_k
