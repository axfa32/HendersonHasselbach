# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 14:55:24 2025

@author: Frank
"""

# Henderson–Hasselbalch titration curve for Histidine (triprotic)
# y-axis: pH, x-axis: equivalents of OH− added (average deprotonation n̄)

import numpy as np
import matplotlib.pyplot as plt

# ---- Inputs ----
# Histidine pKa values and isoelectric point
pKa1 = 2.1
pKa2 = 3.9
pKa3 = 9.8
pI = 3.0

# ---- Calculations ----
# Sort pKa ascending to ensure Ka1 >= Ka2 >= Ka3 in the polyprotic formulas
pKas = np.sort(np.array([pKa1, pKa2, pKa3], dtype=float))
Ka1, Ka2, Ka3 = 10**(-pKas[0]), 10**(-pKas[1]), 10**(-pKas[2])

# Sweep pH and compute distribution fractions
pH = np.linspace(0.0, 14.0, 2000)
H = 10**(-pH)

# Triprotic acid distribution denominator:
# D = [H+]^3 + Ka1[H+]^2 + Ka1Ka2[H+] + Ka1Ka2Ka3
D = (H**3) + (Ka1 * H**2) + (Ka1 * Ka2 * H) + (Ka1 * Ka2 * Ka3)

alpha0 = (H**3) / D                  # H3A
alpha1 = (Ka1 * H**2) / D            # H2A−
alpha2 = (Ka1 * Ka2 * H) / D         # HA2−
alpha3 = (Ka1 * Ka2 * Ka3) / D       # A3−

# Average number of protons lost per molecule (≈ equivalents of OH− added)
nbar = alpha1 + 2*alpha2 + 3*alpha3  # 0 → 3

# ---- Plot ----
plt.figure(figsize=(7,5))
plt.plot(nbar, pH, linewidth=2, label="Histidine titration curve")

# Vertical guides at mid-equivalence points (where pH ≈ each pKa)
mid_eqs = [0.5, 1.5, 2.5]
for x, pka in zip(mid_eqs, pKas):
    plt.axvline(x=x, linestyle="--", linewidth=1)
    # place a small label near the corresponding pH = pKa height
    plt.text(x, pka, f" pKa ≈ {pka:.2f}", va="bottom", ha="left")

# Isoelectric point (horizontal guide)
plt.axhline(y=pI, linestyle=":", linewidth=1, label=f"Isoelectric point (pI ≈ {pI:.2f})")

plt.xlabel("Equivalents of OH⁻ added per mole of histidine (n̄)")
plt.ylabel("pH")
plt.title("Henderson–Hasselbalch Titration Curve: Histidine (triprotic)")
plt.xlim(-0.05, 3.05)
plt.ylim(0, 14)
plt.grid(True, which="both", linestyle=":")
plt.legend(loc="lower right")

plt.tight_layout()
plt.show()
