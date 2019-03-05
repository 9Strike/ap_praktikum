import measure as ms
import numpy as np
import scipy.constants as cs

from measure import npfarray as fa
from measure import pi
from measure import sqrt, sin, cos, tan, exp, ln

ms.plt.rc('text', usetex=True)
ms.plt.rc('font', family='serif')

# (1) Determination of the response time of a RC-element
C1 = fa([470, 4.7, 47]) * cs.nano
d_C1 = 0.10 * C1
R1 = fa([1, 10, 1]) * cs.kilo
d_R1 = 0.05 * R1
T1_12 = fa([0.32, 0.04, 0.04]) * cs.milli
d_T1_12 = fa([0.03, 0.01, 0.01]) * cs.milli
f1 = fa([110, 600, 600])
U1_pp = fa([0.95, 0.95, 0.95])
d_U1_pp = fa([0.02, 0.02, 0.02])

tau1_O = T1_12 / ln(2)
d_tau1_O = d_T1_12 / ln(2)
tau1_T = R1 * C1
d_tau1_T = tau1_T * sqrt((d_R1 / R1)**2 + (d_C1 / C1)**2)

print()
print('1. Determination of the response time of a RC-element:')
print(ms.tbl([ms.lst(C1, d_C1, name='C', unit='F'),
              ms.lst(R1, d_R1, name='R', unit='Ω'),
              ms.lst(f1, name='f', unit='Hz'),
              ms.lst(tau1_O, d_tau1_O, name='τ', unit='s'),
              ms.lst(tau1_T, d_tau1_T, name='τ', unit='s'), ms.dev(tau1_O, d_tau1_O, tau1_T, d_tau1_T, name='τ')]))
print()

# (2) Frequency and phase of a RC-element
f2 = np.arange(1, 11, 1) * cs.kilo
delta_t2 = fa([0.20, 0.08, 0.042, 0.027, 0.019, 0.013, 0.010, 0.007, 0.007, 0.005]) * cs.milli
d_delta_t2 = fa([0.03, 0.02, 0.015, 0.015, 0.010, 0.010, 0.005, 0.005, 0.005, 0.004]) * cs.milli
phi2 = 2 * pi * f2 * delta_t2
d_phi2 = 2 * pi * f2 * d_delta_t2

ms.pltext.initplot(num=1, xlabel=r'$f$ / Hz', ylabel=r'$\varphi$ / $^\circ$')
s2, d_s2, i2, d_i2 = ms.linreg(f2, phi2 / cs.degree, d_phi2 / cs.degree, fit_range=range(3), plot=True)
s2 *= cs.degree
d_s2 *= cs.degree
i2 *= cs.degree
d_i2 *= cs.degree

f2_G = (pi - 4 * i2) / (4 * s2)
d_f2_G = f2_G * sqrt((d_i2 / (pi - 4 * i2))**2 + (d_s2 / s2)**2)

print('2. Frequency and phase of a RC-element:')
print(ms.val('s', s2, d_s2))
print(ms.val('i', i2, d_i2))
print(ms.val('f_G', f2_G, d_f2_G, unit='Hz'))
print(ms.dev(f2_G, d_f2_G, 3.0e3, 0.3e3, name='f_G'))
print()


# Save and show plots
ms.pltext.savefigs('figures/241')
ms.plt.show()