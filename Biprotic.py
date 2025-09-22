# Diprotic Henderson–Hasselbalch titration curve
# y-axis: pH, x-axis: equivalents of OH− added (average deprotonation n̄)

import numpy as np
import matplotlib.pyplot as plt

# ---- Inputs (edit these) ----
pKa1 = 2.46      # example
pKa2 = 9.41      # example
pI   = 5.88      # set to a float to draw a pI line, e.g., 5.97; or leave as None

# ---- Calculations ----
# Ensure pKa1 <= pKa2 for correct ordering
pKa1, pKa2 = sorted([float(pKa1), float(pKa2)])
Ka1, Ka2 = 10**(-pKa1), 10**(-pKa2)

# pH sweep
pH = np.linspace(0.0, 14.0, 2000)
H  = 10**(-pH)

# Diprotic distributions:
# D = [H+]^2 + Ka1[H+] + Ka1Ka2
D     = (H**2) + (Ka1*H) + (Ka1*Ka2)
alpha0 = (H**2)      / D            # H2A
alpha1 = (Ka1*H)     / D            # HA−
alpha2 = (Ka1*Ka2)   / D            # A2−

# Average deprotonation (≈ equivalents of OH− added)
nbar = alpha1 + 2*alpha2            # ranges 0 → 2

# ---- Plot ----
plt.figure(figsize=(7,5))
plt.plot(nbar, pH, linewidth=2, label="Diprotic titration curve")

# Mid-equivalence guides (where pH ≈ pKa1 and pKa2)
for x, pka in zip([0.5, 1.5], [pKa1, pKa2]):
    plt.axvline(x=x, linestyle="--", linewidth=1)
    plt.text(x, pka, f" pKa ≈ {pka:.2f}", va="bottom", ha="left")

# Optional isoelectric point
if pI is not None:
    plt.axhline(y=float(pI), linestyle=":", linewidth=1, label=f"pI ≈ {float(pI):.2f}")

plt.xlabel("Equivalents of OH⁻ added per mole ")
plt.ylabel("pH")
plt.title("Henderson–Hasselbalch Titration Curve (Diprotic): Tryptophan")
plt.xlim(-0.05, 2.05)
plt.ylim(0, 14)
plt.grid(True, which="both", linestyle=":")
plt.legend(loc="lower right")
plt.tight_layout()
plt.show()
