import datplot as dp
import datstr as ds
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as cs
from numpy import sqrt

## General

fa = lambda a: np.array(a, dtype=float)

## Measured data

# General
rho = 7.5 * cs.gram / cs.centi**3

# 3.
tl = [1.03, 17.23]
d_tl = [0.1, 0.1]
Tl = (tl[1] - tl[0]) / 10.0
d_Tl = sqrt(d_tl[0]**2 + d_tl[1]**2) / 10.0

tr = [1.41, 17.56]
d_tr = [0.1, 0.1]
Tr = (tr[1] - tr[0]) / 10.0
d_Tr = sqrt(d_tr[0]**2 + d_tr[1]**2) / 10.0

print(ds.val('Tl', Tl, d_Tl, unit='s'))
print(ds.val('Tr', Tr, d_Tr, unit='s'))
print()

# 4.
tl_s = fa([[1.78, 17.92], [2.91, 19.08], [1.84, 17.92]])
tr_s = fa([[1.79, 17.93], [1.34, 17.43], [1.85, 17.88]])
d_tl_s = fa([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
d_tr_s = fa([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
Tl_s = np.zeros(3)
d_Tl_s = np.zeros(3)
Tr_s = np.zeros(3)
d_Tr_s = np.zeros(3)
for i in range(3):
  Tl_s[i] = (tl_s[i][1] - tl_s[i][0]) / 10.0
  d_Tl_s[i] = sqrt(d_tl_s[i][1]**2 + d_tl_s[i][0]**2) / 10.0
  Tr_s[i] = (tr_s[i][1] - tr_s[i][0]) / 10.0
  d_Tr_s[i] = sqrt(d_tr_s[i][1]**2 + d_tr_s[i][0]**2) / 10.0

print(ds.tbl([
  ds.lst(Tl_s, d_Tl_s, name='Tl_s', unit='s'),
  ds.lst(Tr_s, d_Tr_s, name='Tr_s', unit='s')
]))
print()

tl_a = fa([[1.23, 14.81], [1.26, 15.81], [1.45, 16.82]])
tr_a = fa([[1.92, 15.49], [1.96, 16.53], [2.20, 17.58]])
d_tl_a = fa([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
d_tr_a = fa([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
Tl_a = np.zeros(3)
d_Tl_a = np.zeros(3)
Tr_a = np.zeros(3)
d_Tr_a = np.zeros(3)
for i in range(3):
  Tl_a[i] = (tl_a[i][1] - tl_a[i][0]) / 10.0
  d_Tl_a[i] = sqrt(d_tl_a[i][1]**2 + d_tl_a[i][0]**2) / 10.0
  Tr_a[i] = (tr_a[i][1] - tr_a[i][0]) / 10.0
  d_Tr_a[i] = sqrt(d_tr_a[i][1]**2 + d_tr_a[i][0]**2) / 10.0

print(ds.tbl([
  ds.lst(Tl_a, d_Tl_a, name='Tl_a', unit='s'),
  ds.lst(Tr_a, d_Tr_a, name='Tr_a', unit='s')
]))
print()

# 4.
tl_b = fa([[0.65, 14.06], [2.75, 17.29], [2.09, 17.86]]) # 5 Perioden davor 10
tr_b = fa([[1.10, 15.01], [2.34, 16.90], [1.71, 18.07]])
tl_b_b = fa([[8.78, 52.24], [8.81, 84.62], [17.11, 187.32]])
tr_b_b = fa([[4.68, 47.50], [16.08, 91.78], [33.98, 204.40]]) # Perioden
d_tl_b = fa([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
d_tr_b = fa([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
d_tl_b_b = fa([[0.5, 0.5], [0.5, 0.5], [0.5, 0.5]])
d_tr_b_b = fa([[0.5, 0.5], [0.5, 0.5], [0.5, 0.5]])
Tl_b = np.zeros(3)
d_Tl_b = np.zeros(3)
Tr_b = np.zeros(3)
d_Tr_b = np.zeros(3)
Tl_b_b = np.zeros(3)
d_Tl_b_b = np.zeros(3)
Tr_b_b = np.zeros(3)
d_Tr_b_b = np.zeros(3)
for i in range(3):
  Tl_b[i] = (tl_b[i][1] - tl_b[i][0]) / 10.0
  d_Tl_b[i] = sqrt(d_tl_b[i][0]**2 + d_tl_b[i][1]**2) / 10.0
  Tr_b[i] = (tr_b[i][1] - tr_b[i][0]) / 10.0
  d_Tr_b[i] = sqrt(d_tr_b[i][0]**2 + d_tr_b[i][1]**2) / 10.0
  Tl_b_b[i] = (tl_b_b[i][1] - tl_b_b[i][0]) / 5.0
  d_Tl_b_b[i] = sqrt(d_tl_b_b[i][0]**2 + d_tl_b_b[i][1]**2) / 5.0
  Tr_b_b[i] = (tr_b_b[i][1] - tr_b_b[i][0]) / 5.0
  d_Tr_b_b[i] = sqrt(d_tr_b_b[i][0]**2 + d_tr_b_b[i][1]**2) / 5.0

print(ds.tbl([
  ds.lst(Tl_b, d_Tl_b, name='Tl_b', unit='s'),
  ds.lst(Tr_b, d_Tr_b, name='Tr_b', unit='s'),
  ds.lst(Tl_b_b, d_Tl_b_b, name='Tl_b_b', unit='s'),
  ds.lst(Tr_b_b, d_Tr_b_b, name='Tr_b_b', unit='s')
]))
