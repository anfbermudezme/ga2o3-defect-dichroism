from __future__ import annotations

# Physical constants / conversions used throughout the plotting scripts.
# Keep these in one place so tests and plots share the same values.

HARTREE_TO_EV: float = 27.211386245988  # Hartree energy in eV
HC_EVN_M: float = 1239.841984          # (h*c) in eV*nm, useful for lambda = hc/E
NM_TO_CM: float = 1e-7                 # 1 nm = 1e-7 cm
