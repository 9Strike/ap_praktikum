import measure as ms
from measure import pi, sqrt, sin, cos, arcsin, arccos, arctan
from measure import npfarray as npf
import numpy as np
import scipy.constants as sc
from scipy.optimize import curve_fit

ms.plt.rc('figure', figsize=(11.69, 8.27))

# Data of the helmholtz-coil
d_h = 0.295
r_h = 0.147
n_h = 124

# Data of the induction coil
n_i = 4000
A_i = 41.7 * sc.centi**2

# Measured values
# Resistance and maximum allowed supply voltage
I_max = 5.0
I_test = 1.83
d_I_test = 0.01
U_test = 2.2
d_U_test = 0.1
R_sp = U_test / I_test
d_R_sp = R_sp * ms.sqrt((d_U_test / U_test)**2 + (d_I_test / I_test)**2)
U_max = R_sp * I_max

# Frequency and induced voltage at constant current
I1 = 8.0 / 2
d_I1 = 0.1 / 2
f1 = npf([3.0, 5.6, 9.3, 11.93, 14.90])
d_f1 = npf([0.1, 0.1, 0.1, 0.40, 0.40])
U1_i = npf([0.81, 2.40, 5.08, 6.96, 9.04]) / 2
d_U1_i = npf([0.02, 0.03, 0.04, 0.01, 0.01]) / 2

# Current, voltage and induced voltage at constant frequency
f2 = 9.96
d_f2 = 0.05
I2 = npf([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
d_I2 = npf([0.1, 0.1, 0.1, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1])
U2 = npf([0.6, 1.2, 1.8, 2.4, 3.0, 3.6, 4.2, 4.8, 5.4])
d_U2 = npf([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
U2_i = npf([0.855, 1.540, 2.240, 2.940, 3.560, 4.240, 4.920, 5.620, 6.260]) / 2
d_U2_i = npf([0.010, 0.010, 0.020, 0.020, 0.040, 0.020, 0.040, 0.020, 0.020]) / 2

# Rotation angle and induced voltage at constant a.c. voltage
f3 = 98.03
d_f3 = 1.0
α3 = npf([0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330])
d_α3 = npf([3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3])
U3_i = npf([2.52, 2.24, 1.36, 0.046, 1.25, 2.22, 2.51, 2.14, 1.24, 0.028, 1.19, 2.14]) / 2
d_U3_i = npf([0.01, 0.01, 0.01, 0.002, 0.01, 0.01, 0.01, 0.01, 0.01, 0.002, 0.01, 0.01]) / 2

# Current, voltage, a.c. frequency and induced voltage at a.c. voltage
I4 = npf([0.68, 0.44, 0.31, 0.24, 0.185, 0.15, 0.13, 0.11, 0.095, 0.08, 0.0509, 0.035175, 0.02667, 0.021515, 0.01792, 0.01533, 0.013475, 0.012125, 0.01070])
d_I4 = npf([0.02, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.02, 0.01, 0.0001, 0.00002, 0.00002, 0.00002, 0.00002, 0.00002, 0.00002, 0.00002, 0.00002])
U4 = npf([2.8, 4.20, 4.16, 4.08, 4.08, 4.04, 4.00, 4.00, 4.00, 3.96, 3.92, 3.80, 3.80, 3.80, 3.84, 3.80, 3.76, 3.80, 3.80])
d_U4 = npf([0.2, 0.04, 0.02, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04])
f4 = npf([20.60, 40.20, 60.00, 79.87, 101.2, 123.02, 140.82, 161.13, 179.8, 200.0, 405.4, 603.26, 803.8, 1000, 1205, 1412, 1609, 1805, 2033])
d_f4 = npf([0.03, 0.04, 0.10, 0.20, 0.10, 0.3, 0.3, 0.3, 0.3, 0.4, 0.1, 0.10, 0.4, 2, 2, 3, 3, 3, 3])
U4_i = npf([3.30, 4.40, 4.76, 4.92, 5.00, 5.00, 5.04, 5.04, 5.08, 5.08, 5.08, 5.04, 5.08, 5.08, 5.16, 5.16, 5.16, 5.20, 5.32]) / 2
d_U4_i = npf([0.01, 0.01, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) / 2

# Frequency and induced voltage at no current
f5 = 15.06
d_f5 = 0.10
U5_i = 0.158 / 2
d_U5_i = 0.002 / 2

# Frequency, compensating current and induced voltage
f6 = 15.1
d_f6 = 0.2
I6_c = 0.018205
d_I6_c = 0.00001
U6_i = 0.072 / 2
d_U6_i = 0.001 / 2

# Determination of the magnetic field at the center of the helmholtz coil
ms.pltext.initplot(num=2, title='Induzierte Spannung in Abhängigkeit der Rotationsfrequenz bei konstantem Strom und Gleichspannung', xlabel='f / Hz', ylabel='U / V', fignum=True)
[s1, d_s1, tmp, tmp] = ms.linreg(f1, U1_i, d_U1_i, d_f1, plot=True)
ms.pltext.initplot(num=3, title='Induzierte Spannung in Abhängigkeit des Spulenstroms bei konstanter Rotationsfrequenz und Gleichspannung', xlabel='I / A', ylabel='U / V', fignum=True)
ms.linreg(I2, U2_i, d_U2_i, d_I2, plot=True)

B1 = s1 / (2 * pi * n_i * A_i)
d_B1 = d_s1 / (2 * pi * n_i * A_i)

B1_ = (8 / sqrt(125)) * sc.mu_0 * I1 * n_h / r_h
d_B1_ = (8 / sqrt(125)) * sc.mu_0 * d_I1 * n_h / r_h

print(ms.val("s1", s1, d_s1))
print(ms.val("B1", B1, d_B1, unit='T'))
print(ms.val("B1'", B1_, d_B1_, unit='T'))
print(ms.sig("B1,B1'", B1, d_B1, B1_, d_B1_))
print()

# Determination of the inductance of the helmholtz coil
ms.pltext.initplot(num=4, title='Amplitude der Induktionsspannung in Abhängigkeit des Winkels bei Wechselspannung', xlabel='α / °', ylabel='U / V', fignum=True)
ms.pltext.plotdata(α3, U3_i, d_U3_i, d_α3, connect=True)

Ui_U_ratio = U4_i / U4
d_Ui_U_ratio = Ui_U_ratio * sqrt((d_U4_i / U4_i)**2 + (d_U4 / U4)**2)
ms.pltext.initplot(num=5, title='Verhältnis der Induktionsspannungsamplitude zur Spulenspannungsamplitude\nin Abhängigkeit der Spannungsfrequenz bei Wechselspannung', xlabel='f / Hz', ylabel='Ui / U', fignum=True)
ms.pltext.plotdata(f4, Ui_U_ratio, d_Ui_U_ratio, d_f4, connect=True)

U_I_ratio = U4 / I4
d_U_I_ratio = U_I_ratio * sqrt((d_U4 / U4)**2 + (d_I4 / I4)**2)
ms.pltext.initplot(num=6, title='Widerstand der Helmholtzspule in Abhängigkeit der Spannungsfrequenz bei Welchselspannung', xlabel='f / Hz', ylabel='R / Ω', fignum=True)
[s4, d_s4, i_4, d_i4] = ms.linreg(f4, U_I_ratio, d_U_I_ratio, d_f4, plot=True, fit_range=range(10,len(f4)))

L = s4 / (2 * pi)
d_L = d_s4 / (2 * pi)

print(ms.val("s4", s4, d_s4))
print(ms.val("L", L, d_L, unit='H'))
print()

# Determination of the absolute value of the earths magnetic field
B5 = U5_i / (2 * pi * n_i * A_i * f5)
d_B5 = B5 * sqrt((d_U5_i / U5_i)**2 + (d_f5 / f5)**2)

print(ms.val("B5", B5, d_B5, unit='T'))
print(ms.sig("B5,B5_l", B5, d_B5, 4.87e-5))
print()

# Determination of the angle of inclination of the earths magnetic field
B6_v = (8 / sqrt(125)) * sc.mu_0 * I6_c * n_h / r_h
d_B6_v = (8 / sqrt(125)) * sc.mu_0 * d_I6_c * n_h / r_h
B6_h = U6_i / (2 * pi * n_i * A_i * f6)
d_B6_h = B6_h * sqrt((d_U6_i / U6_i)**2 + (d_f6 / f6)**2)
α = arctan(B6_v / B6_h)
d_α = B6_v * B6_h / (B6_v**2 + B6_h**2) * sqrt((d_B6_v / B6_v)**2 + (d_B6_h / B6_h)**2)
α2 = arcsin(B6_v / B5)
d_α2 = B6_v / sqrt(B5**2 - B6_v**2) * sqrt((d_B6_v / B6_v)**2 + (d_B5 / B5)**2)
α3 = arccos(B6_h / B5)
d_α3 = B6_h / sqrt(B5**2 - B6_h**2) * sqrt((d_B6_h / B6_h)**2 + (d_B5 / B5)**2)

B5_ = sqrt(B6_v**2 + B6_h**2)
d_B5_ = 1 / B5_ * sqrt(B6_v**2 * d_B6_v**2 + B6_h**2 * d_B6_h**2)

α /= sc.degree
d_α /= sc.degree
α2 /= sc.degree
d_α2 /= sc.degree
α3 /= sc.degree
d_α3 /= sc.degree

print(ms.val("B6_v", B6_v, d_B6_v, unit='T'))
print(ms.val("B6_h", B6_h, d_B6_h, unit='T'))
print(ms.val("B6_v'", B5 * sin(66 * sc.degree), unit='T'))
print(ms.val("B6_h'", B5 * cos(66 * sc.degree), unit='T'))
print(ms.val("B5'", B5_, d_B5_, unit='T'))
print(ms.val("α", α, d_α, unit='°', prefix=False))
print(ms.val("α2", α2, d_α2, unit='°', prefix=False))
print(ms.val("α3", α3, d_α3, unit='°', prefix=False))
print(ms.sig("α", α, d_α, 66))
print(ms.sig("α2", α2, d_α2, 66))
print(ms.sig("α3", α3, d_α3, 66))
print()

#ms.plt.show()

#for i in range(2, 7):
#  ms.plt.figure(i).savefig('figures/fig' + str(i) +'.pdf', papertype='a4', orientation='landscape', bbox_inches='tight', format='pdf')
