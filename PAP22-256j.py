import measure as ms
from measure import npfarray as npf
from measure import sqrt
import numpy as np
import scipy.constants as cs

ms.plt.rc('text', usetex=True)
ms.plt.rc('font', family='serif')

titles = [
  r'Abhängigkeit der Energie $E_\alpha$ der $K_\alpha$-Strahlung der Elemente in Abhängigkeit der Kernladungszahl $Z$.',
  r'Abhängigkeit der Energie $E_\beta$ der $K_\beta$-Strahlung der Elemente in Abhängigkeit der Kernladungszahl $Z$.'
]

# Determination of the Rydberg-energy and the screening constant for the K_α-Radiation
Z1 = npf([26, 29, 47, 22, 40, 30, 28, 42])
E_alpha = npf([6.42, 8.05, 21.93, 4.64, 15.80, 8.65, 7.48, 17.47]) * cs.kilo * cs.e
d_E_alpha = npf([0.15, 0.17, 0.20, 0.16, 0.17, 0.17, 0.17, 0.18]) * cs.kilo * cs.e
sr_E_alpha = sqrt(E_alpha)
d_sr_E_alpha = d_E_alpha / (2 * sr_E_alpha)

ms.pltext.initplot(num=1, title=titles[0], xlabel=r'$Z$', ylabel=r'$\sqrt{E_\alpha} / \sqrt{\mathrm{eV}}$')
s1, d_s1, b1, d_b1 = ms.linreg(Z1, sr_E_alpha / sqrt(cs.e), d_sr_E_alpha / sqrt(cs.e), plot=True)
s1 *= sqrt(cs.e)
d_s1 *= sqrt(cs.e)
b1 *= sqrt(cs.e)
d_b1 *= sqrt(cs.e)

E1_R = 4 / 3 * s1**2
d_E1_R = E1_R * (2 * d_s1 / s1)
sigma1_12 = -sqrt(4 / 3) * b1 / sqrt(E1_R)
d_sigma1_12 = sigma1_12 * sqrt((d_b1 / b1)**2 + (d_E1_R / (2 * E1_R))**2)

print()
print("K_α Radiation:")
print(ms.val("s", s1, d_s1))
print(ms.val("b", b1, d_b1))
print()
print(ms.val("E_R", E1_R / cs.e, d_E1_R / cs.e, unit='eV'))
print(ms.sig("E_R,l", E1_R, d_E1_R, cs.Rydberg * cs.h * cs.c))
print(ms.val("σ_12", sigma1_12, d_sigma1_12))
print()

# Determination of the Rydberg-energy and the screening constant for the K_β-Radiation
Z2 = npf([Z1[i] for i in range(len(Z1)) if i != 3])
E_beta = npf([7.06, 8.91, 24.59, 17.69, 9.59, 8.26, 19.58]) * cs.kilo * cs.e
d_E_beta = npf([0.17, 0.17, 0.15, 0.18, 0.17, 0.18, 0.24]) * cs.kilo * cs.e
sr_E_beta = sqrt(E_beta)
d_sr_E_beta = d_E_beta / (2 * sr_E_beta)

ms.pltext.initplot(num=2, title=titles[1], xlabel=r'$Z$', ylabel=r'$\sqrt{E_\beta} / \sqrt{\mathrm{eV}}$')
s2, d_s2, b2, d_b2 = ms.linreg(Z2, sr_E_beta / sqrt(cs.e), d_sr_E_beta / sqrt(cs.e), plot=True)
s2 *= sqrt(cs.e)
d_s2 *= sqrt(cs.e)
b2 *= sqrt(cs.e)
d_b2 *= sqrt(cs.e)

E2_R = 9 / 8 * s2**2
d_E2_R = E2_R * (2 * d_s2 / s2)
sigma2_12 = -sqrt(9 / 8) * b2 / sqrt(E2_R)
d_sigma2_12 = sigma2_12 * sqrt((d_b2 / b2)**2 + (d_E2_R / (2 * E2_R))**2)

print("K_β Radiation:")
print(ms.val("s", s2, d_s2))
print(ms.val("b", b2, d_b2))
print()
print(ms.val("E_R", E2_R / cs.e, d_E2_R / cs.e, unit='eV'))
print(ms.sig("E_R,l", E2_R, d_E2_R, cs.Rydberg * cs.h * cs.c))
print(ms.val("σ_12", sigma2_12, d_sigma2_12))
print()

print("Vergleich")
print(ms.sig("E_R", E1_R, d_E1_R, E2_R, d_E2_R))
print(ms.sig("σ_12", sigma1_12, d_sigma1_12, sigma2_12, d_sigma2_12))
print()

ms.pltext.savefigs('figures/256')
ms.plt.show()
