import measure as ms
import numpy as np
import scipy.constants as cs

from measure import npfarray as npfarray
from measure import pi
from measure import sqrt, sin, cos, tan, exp, ln

ms.plt.rc('text', usetex=True)
ms.plt.rc('font', family='serif')

titles = [
  r'Phasenverschiebung $\varphi$ eines Hochpassfilters ($C = 47$ nF, $R = 1$ k$\Omega$) in Abhängigkeit der Frequenz $f$.'
]

# (1) Determination of the response time of a RC-element
C1 = npfarray([470, 4.7, 47]) * cs.nano
d_C1 = 0.10 * C1
R1 = npfarray([1, 10, 1]) * cs.kilo
d_R1 = 0.05 * R1
T1_12 = npfarray([0.32, 0.04, 0.04]) * cs.milli
d_T1_12 = npfarray([0.03, 0.01, 0.01]) * cs.milli
f1 = npfarray([110, 600, 600])
U1_pp = npfarray([0.95, 0.95, 0.95])
d_U1_pp = npfarray([0.02, 0.02, 0.02])

tau1_O = T1_12 / ln(2)
d_tau1_O = d_T1_12 / ln(2)
tau1_T = R1 * C1
d_tau1_T = tau1_T * sqrt((d_R1 / R1)**2 + (d_C1 / C1)**2)

print()
print('1. Determination of the response time of a RC-element:')
print(ms.tbl([ms.lst(C1, d_C1, name='C', unit='F'),
              ms.lst(R1, d_R1, name='R', unit='Ω'),
              ms.lst(tau1_O, d_tau1_O, name='τ', unit='s'),
              ms.lst(tau1_T, d_tau1_T, name='τ', unit='s'), ms.dev(tau1_O, d_tau1_O, tau1_T, d_tau1_T, name='τ')]))

# (3) Frequency and phase of a RC-element
R3 = cs.kilo
d_R3 = 0.05 * R3
C3 = 47 * cs.nano
d_C3 = 0.10 * C3
f3_G_low = 3.0 * cs.kilo
d_f3_G_low = 0.3 * cs.kilo
f3_G_high = 3.1 * cs.kilo
d_f3_G_high = 0.3 * cs.kilo

f3 = np.arange(1, 11, 1) * cs.kilo
delta_t3 = npfarray([0.20, 0.08, 0.042, 0.027, 0.019, 0.013, 0.010, 0.007, 0.007, 0.005]) * cs.milli
d_delta_t3 = npfarray([0.03, 0.02, 0.015, 0.015, 0.010, 0.010, 0.005, 0.005, 0.005, 0.004]) * cs.milli
phi3 = 2 * pi * f3 * delta_t3
d_phi3 = 2 * pi * f3 * d_delta_t3

ms.pltext.initplot(num=7, title=titles[0], xlabel=r'$f$ / kHz', ylabel=r'$\varphi$ / $^\circ$', fignum=True)
s3, d_s3, i3, d_i3 = ms.linreg(f3 / cs.kilo, phi3 / cs.degree, d_phi3 / cs.degree, fit_range=range(3), plot=True)
s3 *= cs.degree / cs.kilo
d_s3 *= cs.degree / cs.kilo
i3 *= cs.degree
d_i3 *= cs.degree

f3_G_low_p1 = 2.9 * cs.kilo
d_f3_G_low_p1 = 0.3 * cs.kilo
f3_G_high_p1 = 3.1 * cs.kilo
d_f3_G_high_p1 = 0.3 * cs.kilo

f3_G_high_p2 = (pi - 4 * i3) / (4 * s3)
d_f3_G_high_p2 = f3_G_high_p2 * sqrt((4 * d_i3 / (pi - 4 * i3))**2 + (d_s3 / s3)**2)

f3_G_T = 1 / (2 * pi * R3 * C3)
d_f3_G_T = f3_G_T * sqrt((d_R3 / R3)**2 + (d_C3 / C3)**2)

print('3. Frequency and phase of a RC-element:')
print(ms.val('s', s3, d_s3))
print(ms.val('i', i3, d_i3))
print(ms.val('low: f_G', f3_G_low, d_f3_G_low, unit='Hz'))
print(ms.val('high: f_G', f3_G_high, d_f3_G_high, unit='Hz'))
print(ms.val('low: f_G', f3_G_low_p1, d_f3_G_low_p1, unit='Hz'))
print(ms.val('high: f_G', f3_G_high_p1, d_f3_G_high_p1, unit='Hz'))
print(ms.val('high: f_G', f3_G_high_p2, d_f3_G_high_p2, unit='Hz'))
print(ms.val('f_G', f3_G_T, d_f3_G_T, unit='Hz'))
print(ms.dev(f3_G_low, d_f3_G_low, f3_G_T, d_f3_G_T, name='low: f_G'))
print(ms.dev(f3_G_high, d_f3_G_high, f3_G_T, d_f3_G_T, name='high: f_G'))
print(ms.dev(f3_G_low_p1, d_f3_G_low_p1, f3_G_T, d_f3_G_T, name='low: f_G'))
print(ms.dev(f3_G_high_p1, d_f3_G_high_p1, f3_G_T, d_f3_G_T, name='high: f_G'))
print(ms.dev(f3_G_high_p2, d_f3_G_high_p2, f3_G_T, d_f3_G_T, name='high: f_G'))
print()

# (4) Frequency of an oscillating circuit
C4 = 47 * cs.nano
d_C4 = 0.10 * C4
R4 = npfarray([1000, 220, 47])
d_R4 = 0.05 * R4

U4_E = 1.00
d_U4_E = 0.01
U4_A = npfarray([0.94, 0.73, 0.31])
d_U4_A = npfarray([0.01, 0.01, 0.01])
f4_R = npfarray([3.85, 3.70, 3.70]) * cs.kilo
d_f4_R = npfarray([0.05, 0.05, 0.05]) * cs.kilo
f4_1 = npfarray([2.15, 3.18, 3.37]) * cs.kilo
d_f4_1 = npfarray([0.05, 0.05, 0.05]) * cs.kilo
f4_2 = npfarray([7.08, 4.41, 4.00]) * cs.kilo
d_f4_2 = npfarray([0.05, 0.05, 0.05]) * cs.kilo
delta_omega4 = 2 * pi * (f4_2 - f4_1)
d_delta_omega4 = 2 * pi * sqrt(d_f4_2**2 + d_f4_1**2)

L4 = 1 / ((2 * pi * f4_R)**2 * C4)
d_L4 = L4 * sqrt((2 * d_f4_R / f4_R)**2 + (d_C4 / C4)**2)
d_L4 = 1 / len(L4) * sqrt(np.sum(d_L4**2))
L4 = ms.mv(L4)

R4_G = delta_omega4 * L4
d_R4_G = R4_G * sqrt((d_delta_omega4 / delta_omega4)**2 + (d_L4 / L4)**2)

R4_V1 = (U4_E / U4_A - 1) * R4
d_R4_V1 = sqrt((R4 / U4_A)**2 * (d_U4_E**2 + (U4_E * d_U4_A / U4_A)**2) + ((U4_E / U4_A - 1) * d_R4)**2)

R4_V2 = R4_G - R4
d_R4_V2 = sqrt(d_R4_G**2 + d_R4**2)

print('4. Frequency of an oscillating circuit')
print(ms.val('L', L4, d_L4, unit='H'))
print()

print(ms.tbl([ms.lst(R4, d_R4, name='R', unit='Ω'),
              ms.lst(delta_omega4, d_delta_omega4, name='Δω', unit='Hz'),
              ms.lst(R4_G, d_R4_G, name='R_G', unit='Ω'),
              ms.lst(R4_V1, d_R4_V1, name='R_V', unit='Ω'),
              ms.lst(R4_V2, d_R4_V2, name='R_V', unit='Ω'),
              ms.dev(R4_V1, d_R4_V1, R4_V2, d_R4_V2, name='R_V')]))

# (5) Determination of the dampting constant of a free, damped oscillating circuit
L5 = L4
d_L5 = d_L4
A5 = npfarray([1.83, 1.17, 0.80, 0.56, 0.38])
d_A5 = npfarray([0.10, 0.10, 0.10, 0.10, 0.10])
T5 = npfarray([0.26, 0.26, 0.26, 0.26]) * cs.milli
d_T5 = npfarray([0.03, 0.03, 0.03, 0.03]) * cs.milli

Ld5 = ln(A5[:-1] / A5[1:])
d_Ld5 = np.sum((d_A5[:-1] / A5[:-1])**2 + (d_A5[1:] / A5[1:])**2) / (len(A5) - 1)
Ld5 = ms.mv(Ld5)

d_T5 = np.sum(d_T5**2) / len(T5)
T5 = ms.mv(T5)

delta5 = Ld5 / T5
d_delta5 = delta5 * sqrt((d_Ld5 / Ld5)**2 + (d_T5 / T5)**2)

R5_G = 2 * L5 * delta5
d_R5_G = R5_G * sqrt((d_L5 / L5)**2 + (d_delta5 / delta5)**2)

print('5. Determination of the dampting constant of a free, damped oscillating circuit')
print(ms.val('Λ', Ld5, d_Ld5))
print(ms.val('T', T5, d_T5, unit='s'))
print(ms.val('δ', delta5, d_delta5, unit='1/s', prefix=False))
print(ms.val('R_G', R5_G, d_R5_G, unit='Ω'))
print(ms.dev(R4_G[2], d_R4_G[2], R5_G, d_R5_G, name='R_G'))
print()

# (6) Resonance magnification
C6 = 47 * cs.nano
d_C6 = 0.10 * C6
R6 = 220
d_R6 = 0.05 * R6
L6 = L4
d_L6 = d_L4
omega6_R = 2 * pi * 3.84 * cs.kilo
d_omega6_R = 2 * pi * 0.05 * cs.kilo
omega6_C = 2 * pi * 3.75 * cs.kilo
d_omega6_C = 2 * pi * 0.05 * cs.kilo
omega6_L = 2 * pi * 3.93 * cs.kilo
d_omega6_L = 2 * pi * 0.05 * cs.kilo

delta6 = R6 / (2 * L6)
d_delta6 = delta6 * sqrt((d_R6 / R6)**2 + (d_L6 / L6)**2)
omega6_R_T = 1 / sqrt(L6 * C6)
d_omega6_R_T = omega6_R_T / 2 * sqrt((d_L6 / L6)**2 + (d_C6 / C6)**2)
omega6_C_T = sqrt(omega6_R_T**2 - 2 * delta6**2)
d_omega6_C_T = sqrt((2 * omega6_R_T * d_omega6_R_T)**2 + (4 * delta6 * d_delta6)**2) / (2 * omega6_C_T)
omega6_L_T = sqrt(omega6_R_T**2 + 2 * delta6**2)
d_omega6_L_T = sqrt((2 * omega6_R_T * d_omega6_R_T)**2 + (4 * delta6 * d_delta6)**2) / (2 * omega6_L_T)

print('6. Resonance magnification')
print(ms.val('ω_R', omega6_R, d_omega6_R, unit='Hz'))
print(ms.val('ω_R', omega6_R_T, d_omega6_R_T, unit='Hz'))
print(ms.dev(omega6_R, d_omega6_R, omega6_R_T, d_omega6_R_T, name='ω_R'))
print(ms.val('ω_C', omega6_C, d_omega6_C, unit='Hz'))
print(ms.val('ω_C', omega6_C_T, d_omega6_C_T, unit='Hz'))
print(ms.dev(omega6_C, d_omega6_C, omega6_C_T, d_omega6_C_T, name='ω_C'))
print(ms.val('ω_L', omega6_L, d_omega6_L, unit='Hz'))
print(ms.val('ω_L', omega6_L_T, d_omega6_L_T, unit='Hz'))
print(ms.dev(omega6_L, d_omega6_L, omega6_L_T, d_omega6_L_T, name='ω_L'))
print()

# (7) Band-stop filter
R7 = cs.kilo
d_R7 = 0.05 * R7
C7 = 47 * cs.nano
d_C7 = 0.10 * C7
L7 = L4
d_L7 = d_L4
omega7 = 2 * pi * 3.89 * cs.kilo
d_omega7 = 2 * pi * 0.15 * cs.kilo

omega7_T = 1 / sqrt(L7 * C7)
d_omega7_T = omega7_T / 2 * sqrt((d_L7 / L7)**2 + (d_C7 / C7)**2)

print('7. band-stop filter')
print(ms.val('ω', omega7, d_omega7, unit='Hz'))
print(ms.val('ω', omega7_T, d_omega7_T, unit='Hz'))
print(ms.dev(omega7, d_omega7, omega7_T, d_omega7_T, name='ω'))
print()

# Save and show plots
ms.pltext.savefigs('figures/241')
#ms.plt.show()
