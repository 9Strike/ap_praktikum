import datplot as dp
import datstr as ds
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as cs
from numpy import pi, sqrt

### General

print()

fa = lambda a: np.array(a, dtype=float)



### Measured data

## General
rho = 7.5 * cs.gram / cs.centi**3
l = fa([393.0, 293.0, 193.0]) * cs.milli
d_l = fa([3.0, 3.0, 3.0]) * cs.milli


## Without coupling
tl = fa([1.03, 17.23])
d_tl = fa([0.1, 0.1])

tr = fa([1.41, 17.56])
d_tr = fa([0.1, 0.1])


## Symmetric oscillation
tl_s = fa([[1.78, 17.92], [2.91, 19.08], [1.84, 17.92]])
tr_s = fa([[1.79, 17.93], [1.34, 17.43], [1.85, 17.88]])
d_tl_s = fa([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
d_tr_s = fa([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
f_s = fa([0.620, 0.620, 0.620])
d_f_s = fa([0.008, 0.005, 0.008])


## Antisymmetric oscillation
tl_a = fa([[1.23, 14.81], [1.26, 15.81], [1.45, 16.82]])
tr_a = fa([[1.92, 15.49], [1.96, 16.53], [2.20, 17.58]])
d_tl_a = fa([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
d_tr_a = fa([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
f_a = fa([0.736, 0.686, 0.649])
d_f_a = fa([0.009, 0.01, 0.009])


## Mixed excitation
tl_b = fa([[0.65, 14.06], [2.75, 17.29], [2.09, 17.86]]) # 5 Periods and 10 with beats
tr_b = fa([[1.10, 15.01], [2.34, 16.90], [1.71, 18.07]])
tl_b_b = fa([[8.78, 52.24], [8.81, 84.62], [17.11, 187.32]])
tr_b_b = fa([[4.68, 47.50], [16.08, 91.78], [33.98, 204.40]])
d_tl_b = fa([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
d_tr_b = fa([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
d_tl_b_b = fa([[0.5, 0.5], [0.5, 0.5], [0.5, 0.5]])
d_tr_b_b = fa([[0.5, 0.5], [0.5, 0.5], [0.5, 0.5]])
f_b = fa([0.735, 0.686, 0.649])
d_f_b = fa([0.009, 0.006, 0.003])
f_b_b = fa([0.620, 0.620, 0.619])
d_f_b_b = fa([0.009, 0.006, 0.003])



### Data preparation

## Without coupling
Tl = (tl[1] - tl[0]) / 10.0
d_Tl = sqrt(d_tl[0]**2 + d_tl[1]**2) / 10.0

Tr = (tr[1] - tr[0]) / 10.0
d_Tr = sqrt(d_tr[0]**2 + d_tr[1]**2) / 10.0

omega_l = 2.0 * pi / Tl
d_omega_l = omega_l * d_Tl / Tl
omega_r = 2.0 * pi / Tr
d_omega_r = omega_r * d_Tr / Tr

omega = 0.5 * (omega_l + omega_r)
d_omega = 0.5 * sqrt(d_omega_l**2 + d_omega_r**2)


## Symmetric oscillation
Tl_s = np.zeros(3)
d_Tl_s = np.zeros(3)
Tr_s = np.zeros(3)
d_Tr_s = np.zeros(3)
for i in range(3):
  Tl_s[i] = (tl_s[i][1] - tl_s[i][0]) / 10.0
  d_Tl_s[i] = sqrt(d_tl_s[i][1]**2 + d_tl_s[i][0]**2) / 10.0
  Tr_s[i] = (tr_s[i][1] - tr_s[i][0]) / 10.0
  d_Tr_s[i] = sqrt(d_tr_s[i][1]**2 + d_tr_s[i][0]**2) / 10.0

omega_l_s = 2.0 * pi / Tl_s
d_omega_l_s = omega_l_s * d_Tl_s / Tl_s
omega_r_s = 2.0 * pi / Tr_s
d_omega_r_s = omega_r_s * d_Tr_s / Tr_s

omega_s = 0.5 * (omega_l_s + omega_r_s)
d_omega_s = 0.5 * sqrt(d_omega_l_s**2 + d_omega_r_s**2)

omega2_s = 2.0 * pi * f_s
d_omega2_s = 2.0 * pi * d_f_s


## Antisymmetric oscillation
Tl_a = np.zeros(3)
d_Tl_a = np.zeros(3)
Tr_a = np.zeros(3)
d_Tr_a = np.zeros(3)
for i in range(3):
  Tl_a[i] = (tl_a[i][1] - tl_a[i][0]) / 10.0
  d_Tl_a[i] = sqrt(d_tl_a[i][1]**2 + d_tl_a[i][0]**2) / 10.0
  Tr_a[i] = (tr_a[i][1] - tr_a[i][0]) / 10.0
  d_Tr_a[i] = sqrt(d_tr_a[i][1]**2 + d_tr_a[i][0]**2) / 10.0

omega_l_a = 2.0 * pi / Tl_a
d_omega_l_a = omega_l_a * d_Tl_a / Tl_a
omega_r_a = 2.0 * pi / Tr_a
d_omega_r_a = omega_r_a * d_Tr_a / Tr_a

omega_a = 0.5 * (omega_l_a + omega_r_a)
d_omega_a = 0.5 * sqrt(d_omega_l_a**2 + d_omega_r_a**2)

omega2_a = 2.0 * pi * f_a
d_omega2_a = 2.0 * pi * d_f_a


## Mixed excitation
Tl_b, d_Tl_b = np.zeros(3), np.zeros(3)
Tr_b, d_Tr_b = np.zeros(3), np.zeros(3)
Tl_b_b, d_Tl_b_b = np.zeros(3), np.zeros(3)
Tr_b_b, d_Tr_b_b = np.zeros(3), np.zeros(3)
for i in range(3):
  Tl_b[i] = (tl_b[i][1] - tl_b[i][0]) / 10.0
  d_Tl_b[i] = sqrt(d_tl_b[i][0]**2 + d_tl_b[i][1]**2) / 10.0
  Tr_b[i] = (tr_b[i][1] - tr_b[i][0]) / 10.0
  d_Tr_b[i] = sqrt(d_tr_b[i][0]**2 + d_tr_b[i][1]**2) / 10.0
  Tl_b_b[i] = 2.0 * (tl_b_b[i][1] - tl_b_b[i][0]) / 5.0
  d_Tl_b_b[i] = 2.0 * sqrt(d_tl_b_b[i][0]**2 + d_tl_b_b[i][1]**2) / 5.0
  Tr_b_b[i] = 2.0 * (tr_b_b[i][1] - tr_b_b[i][0]) / 5.0
  d_Tr_b_b[i] = 2.0 * sqrt(d_tr_b_b[i][0]**2 + d_tr_b_b[i][1]**2) / 5.0

omega_l_b = 2.0 * pi / Tl_b
d_omega_l_b = omega_l_b * d_Tl_b / Tl_b
omega_r_b = 2.0 * pi / Tr_b
d_omega_r_b = omega_r_b * d_Tr_b / Tr_b
omega_l_b_b = 2.0 * pi / Tl_b_b
d_omega_l_b_b = omega_l_b_b * d_Tl_b_b / Tl_b_b
omega_r_b_b = 2.0 * pi / Tr_b_b
d_omega_r_b_b = omega_r_b_b * d_Tr_b_b / Tr_b_b

omega_b = 0.5 * (omega_l_b + omega_r_b)
d_omega_b = 0.5 * sqrt(d_omega_l_b**2 + d_omega_r_b**2)
omega_b_b = 0.5 * (omega_l_b_b + omega_r_b_b)
d_omega_b_b = 0.5 * sqrt(d_omega_l_b_b**2 + d_omega_r_b_b**2)

omega2_2_b = 2.0 * pi * f_b
d_omega2_2_b = 2.0 * pi * d_f_b
omega2_1_b = 2.0 * pi * f_b_b
d_omega2_1_b = 2.0 * pi * d_f_b_b



### Evaluation

## Without coupling
print(ds.val('ω', omega, d_omega))
print()


## Symmetric Oscillation

omega_array = np.full_like(omega_s, omega)
d_omega_array = np.full_like(omega_s, d_omega)
print(ds.tbl([
  ds.lst(omega_s, d_omega_s, name='ω1', prefix=False, unit='1/s'),
  ds.lst(omega2_s, d_omega2_s, name='ω2', prefix=False, unit='1/s'),
  ds.dev(omega_s, d_omega_s, omega2_s, d_omega2_s, name='ω1, ω2'),
  ds.dev(omega_s, d_omega_s, omega_array, d_omega_array, name='ω1, ω'),
  ds.dev(omega2_s, d_omega2_s, omega_array, d_omega_array, name='ω2, ω')
], name='Symmetric oscillation frequencys'))
print()


## Antisymmetric oscillation

print(ds.tbl([
  ds.lst(omega_a, d_omega_a, name='ω1', prefix=False, unit='1/s'),
  ds.lst(omega2_a, d_omega2_a, name='ω2', prefix=False, unit='1/s'),
  ds.dev(omega_a, d_omega_a, omega2_a, d_omega2_a, name='ω1, ω2')
], name='Antisymmetric oscillation frequencys'))
print()


## Mixed excitation

omega2_b = 0.5 * (omega2_1_b + omega2_2_b)
d_omega2_b = 0.5 * sqrt(d_omega2_1_b**2 + d_omega2_2_b**2)
omega2_b_b = 0.5 * (omega2_2_b - omega2_1_b)
d_omega2_b_b = 0.5 * sqrt(d_omega2_1_b**2 + d_omega2_2_b**2)

omega_theo_b = 0.5 * (omega_s + omega_a)
d_omega_theo_b = 0.5 * sqrt(d_omega_s**2 + d_omega_a**2)
omega_theo_b_b = 0.5 * (omega_a - omega_s)
d_omega_theo_b_b = 0.5 * sqrt(d_omega_s**2 + d_omega_a**2)

print(ds.tbl([
  ds.lst(omega_b, d_omega_b, name='ω1', prefix=False, unit='1/s'),
  ds.lst(omega2_b, d_omega2_b, name='ω2', prefix=False, unit='1/s'),
  ds.lst(omega_theo_b, d_omega_theo_b, name='ω_E', prefix=False, unit='1/s'),
  ds.dev(omega_b, d_omega_b, omega2_b, d_omega2_b, name='ω1, ω2'),
  ds.dev(omega_b, d_omega_b, omega_theo_b, d_omega_theo_b, name='ω1, ω_E'),
  ds.dev(omega2_b, d_omega2_b, omega_theo_b, d_omega_theo_b, name='ω2, ω_E'),
], name='Mixed excitation frequencys'))
print(ds.tbl([
  ds.lst(omega_b_b, d_omega_b_b, name='ω1', prefix=False, unit='1/s'),
  ds.lst(omega2_b_b, d_omega2_b_b, name='ω2', prefix=False, unit='1/s'),
  ds.lst(omega_theo_b_b, d_omega_theo_b_b, name='ω_E', prefix=False, unit='1/s'),
  ds.dev(omega_b_b, d_omega_b_b, omega2_b_b, d_omega2_b_b, name='ω1, ω2'),
    ds.dev(omega_b_b, d_omega_b_b, omega_theo_b_b, d_omega_theo_b_b, name='ω1, ω_E'),
  ds.dev(omega2_b_b, d_omega2_b_b, omega_theo_b_b, d_omega_theo_b_b, name='ω2, ω_E'),
], name='Mixed excitation beats frequencys'))
print()

kappa = (omega_a**2 - omega_s**2) / (omega_s**2 + omega_a**2)
d_kappa = kappa * sqrt(4 * (omega_a**2 * d_omega_a**2 + omega_s**2 * d_omega_s**2)
                        * (1 / (omega_a**2 - omega_s**2)**2 + 1 / (omega_s**2 + omega_a**2)**2))

print(ds.tbl([
  ds.lst(kappa, d_kappa, name='κ', prefix=False)
], name='coupling factors'))

kappa_ratio = np.zeros(len(kappa) - 1)
d_kappa_ratio = np.zeros(len(d_kappa) - 1)
l2_ratio = np.zeros(len(l) - 1)
d_l2_ratio = np.zeros(len(d_l) - 1)
for i in range(len(kappa_ratio)):
  kappa_ratio[i] = kappa[i + 1] / kappa[i]
  d_kappa_ratio[i] = kappa_ratio[i] * sqrt((d_kappa[i + 1] / kappa[i + 1])**2 + (d_kappa[i] / kappa[i])**2)
  l2_ratio[i] = l[i + 1]**2 / l[i]**2
  d_l2_ratio[i] = l2_ratio[i] * sqrt((2 * d_l[i + 1] / l[i + 1])**2 + (2 * d_l[i] / l[i])**2)

print(ds.tbl([
  ds.lst(kappa_ratio, d_kappa_ratio, name='κ_ratio', prefix=False),
  ds.lst(l2_ratio, d_l2_ratio, name='l²_ratio', prefix=False),
  ds.dev(kappa_ratio, d_kappa_ratio, l2_ratio, d_l2_ratio, name='κ_ratio, l²_ratio')
], name='coupling factor ratios'))
