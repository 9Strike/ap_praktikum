import datplot as dp
import datstr as ds
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sqrt, sin, cos, sinc
from scipy.integrate import quad

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fa = lambda a: np.array(a, dtype=float)

# The cut is performed after cut as index. lshift = 0 returns original array.
def dat_overlap(arr, d_arr, cut, lshift):
  c = cut
  s = lshift
  arr_ = np.concatenate((arr[:c-s], np.mean([arr[c-s:c], arr[c:c+s]], axis=0)))
  arr_ = np.concatenate((arr_, arr[c+s:]))
  d_arr_ = np.concatenate((d_arr[:c-s], sqrt(d_arr[c-s:c]**2 + d_arr[c:c+s]**2) / 2))
  d_arr_ = np.concatenate((d_arr_, d_arr[c+s:]))
  return arr_, d_arr_

## Measured data
# Maxima, minima of the single slit fourier image
n_sf = fa([0, 1, 1, 2, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6])
I_ug_sf = 72.0
d_I_ug_sf = 14.0
I_sf = fa([3700, 90, 290, 84, 3340, 130, 1330, 110, 640, 100, 410, 96, 287, 99])
d_I_sf = fa([50, 10, 20, 10, 50, 20, 40, 20, 30, 20, 30, 10, 30, 10])
x_sf = fa([1073, 992, 955, 913, 957, 912, 877, 834, 794, 755, 717, 676, 639, 598])
d_x_sf = fa([5] * len(x_sf))

I_sf -= I_ug_sf
d_I_sf = sqrt(d_I_sf**2 + d_I_ug_sf**2)
I_sf_ = np.concatenate((
  I_sf[:4] / I_sf[0],
  I_sf[4:] / I_sf[4] * I_sf[2] / I_sf[0]
))
d_I_sf_ = np.concatenate((
  I_sf_[:4] * sqrt((d_I_sf[:4] / I_sf[:4])**2 + (d_I_sf[0] / I_sf[0])**2),
  I_sf_[4:] * sqrt((d_I_sf[4:] / I_sf[4:])**2 + (d_I_sf[4] / I_sf[4])**2 + (d_I_sf[2] / I_sf[2])**2 + (d_I_sf[0] / I_sf[0])**2)
))

n_sf = np.concatenate((n_sf[:4], n_sf[6:]))
I_sf, d_I_sf = dat_overlap(I_sf_, d_I_sf_, 4, 2)
x_sf, d_x_sf = dat_overlap(x_sf, d_x_sf, 4, 2)

n_min_sf = n_sf[1::2]
I_min_sf = I_sf[1::2]
d_I_min_sf = d_I_sf[1::2]
x_min_sf = x_sf[1::2]
d_x_min_sf = d_x_sf[1::2]

# n_max_sf = n_sf[::2]
I_max_sf = I_sf[::2]
d_I_max_sf = d_I_sf[::2]
x_max_sf = x_sf[::2]
d_x_max_sf = d_x_sf[::2]

# Maxima, minima of the double slit fourier image
n_df = fa([0, 1, 1, 2, 1, 2, 2, 3, 3, 4])
I_ug_df = 72.0
d_I_ug_df = 14.0
I_df = fa([3530, 120, 2230, 90, 3590, 100, 390, 90, 435, 90])
d_I_df = fa([50, 30, 50, 20, 50, 20, 40, 20, 30, 20])
x_df = fa([1219, 1202, 1188, 1169, 1184, 1164, 1155, 1140, 1114, 1098])
d_x_df = fa([5, 5, 5, 5, 5, 5, 5, 30, 5, 5])

I_df -= I_ug_df
d_I_df = sqrt(d_I_df**2 + d_I_ug_df**2)
I_df_ = np.concatenate((
  I_df[:4] / I_df[0],
  I_df[4:] / I_df[4] * I_df[2] / I_df[0]
))
d_I_df_ = np.concatenate((
  I_df_[:4] * sqrt((d_I_df[:4] / I_df[:4])**2 + (d_I_df[0] / I_df[0])**2),
  I_df_[4:] * sqrt((d_I_df[4:] / I_df[4:])**2 + (d_I_df[4] / I_df[4])**2 + (d_I_df[2] / I_df[2])**2 + (d_I_df[0] / I_df[0])**2)
))

n_df = np.concatenate((n_df[:4], n_df[6:]))
I_df, d_I_df = dat_overlap(I_df_, d_I_df_, 4, 2)
x_df, d_x_df = dat_overlap(x_df, d_x_df, 4, 2)

n_min_df = n_df[1::2]
I_min_df = I_df[1::2]
d_I_min_df = d_I_df[1::2]
x_min_df = x_df[1::2]
d_x_min_df = d_x_df[1::2]

# n_max_df = n_df[::2]
I_max_df = I_df[::2]
d_I_max_df = d_I_df[::2]
x_max_df = x_df[::2]
d_x_max_df = d_x_df[::2]

# Abscissa calibration with the single slit
n_sc = fa([1, 2, 3, 4])
d_sc = fa([0.235, 0.440, 0.650, 0.860])
d_d_sc = fa([0.005, 0.005, 0.005, 0.005])

# Single slit object image
I1_so = fa([1580, 90])
I2_so = fa([1080, 1580, 86])
I3_so = fa([1440, 1120, 1530, 89])
d_I1_so = fa([40, 3])
d_I2_so = fa([20, 20, 5])
d_I3_so = fa([30, 30, 50, 5])

I1_i_max = np.argmax(I1_so)[0]
I2_i_max = np.argmax(I2_so)[0]
I3_i_max = np.argmax(I3_so)[0]
d_I1_so = I1_so / I1_so[I1_i_max] * sqrt((d_I1_so / I1_so)**2 + (d_I1_so[I1_i_max] / I1_so[I1_i_max])**2)
d_I2_so = I2_so / I2_so[I2_i_max] * sqrt((d_I2_so / I2_so)**2 + (d_I2_so[I2_i_max] / I2_so[I2_i_max])**2)
d_I3_so = I3_so / I3_so[I3_i_max] * sqrt((d_I3_so / I3_so)**2 + (d_I3_so[I3_i_max] / I3_so[I3_i_max])**2)
I1_so = I1_so / I1_so[I1_i_max]
I2_so = I2_so / I2_so[I2_i_max]
I3_so = I3_so / I3_so[I3_i_max]

x1_so = 939
d_x1_so = 3
x2_so = 1070
d_x2_so = 3
b_so = 682
d_b_so = 5
f_so = 80

## Evaluation
# Single slit fourier image
dp.initplot(num=1, title=r'Positionen $x$ der Minima und Maxima eines Einzelspaltes in Abhängigkeit der Ordnung $n$.', xlabel=r'$n$', ylabel=r'$x$ / px')
s, d_s, b, d_b = dp.linreg(n_min_sf, x_min_sf, d_x_min_sf, plot=True)
n_max_sf = (x_max_sf - b) / s
d_n_max_sf = abs(n_max_sf) * sqrt((d_x_max_sf**2 + d_b**2) / (x_max_sf - b)**2 + (d_s / s)**2)
dp.plotdata(n_max_sf, x_max_sf, d_x_max_sf, d_n_max_sf)

n_max_sf_theo = np.arange(0.5, 0.5 + len(n_max_sf))
n_max_sf_theo[0] = 0.0
I_max_sf_theo = sinc(n_max_sf)**2
d_I_max_sf_theo = 2 * abs((sinc(n_max_sf) - cos(pi * n_max_sf)) * sinc(n_max_sf) * d_n_max_sf / n_max_sf)

print()
print(ds.tbl([
  ds.lst(n_max_sf, d_n_max_sf, name='n_o', expToFix=0),
  ds.lst(n_max_sf_theo, name='n_t', expToFix=0),
  ds.sig('n_o, n_t', n_max_sf, d_n_max_sf, n_max_sf_theo, perc=True)
]))
print(ds.tbl([
  ds.lst(I_max_sf, d_I_max_sf, name='I_o', unit='I0', prefix=False, expToFix=0),
  ds.lst(I_max_sf_theo, d_I_max_sf_theo, name='I_t', unit='I0', prefix=False, expToFix=0),
  ds.sig('I_o, I_t', I_max_sf, d_I_max_sf, I_max_sf_theo, d_I_max_sf_theo, perc=True)
]))

# Double slit fourier image


# Single slit object image
def slit(n, y):
  return 2 * sinc(n) * cos(2 * pi * n * y)

n_so = 3
y_so = np.linspace(-1, 1, 1000)

I_so = np.array([quad(lambda n: slit(n, y), 0, n_so) for y in y_so])
I_so = np.array([x[0]**2 for x in I_so])
I_so = I_so / np.max(I_so)

dp.initplot(num=2, title=r'', xlabel=r'$y$ / $d$', ylabel=r'$I$ / b.E.')
dp.plot(y_so, I_so)

# Double slit object image
def dslit(n, y, g):
  return 4 * sinc(n) * cos(pi * n * g) * cos(2 * pi * n * y)

g_do = 2
n_do = 3
y_do = np.linspace(-2, 2, 1000)

I_do = np.array([quad(lambda n: dslit(n, y, g_do), 0, n_do) for y in y_do])
I_do = np.array([x[0]**2 for x in I_do])
I_do = I_do / np.max(I_do)

dp.initplot(num=3, title=r'', xlabel=r'$y$ / $d$', ylabel=r'$I$ / b.E.')
dp.plot(y_do, I_do)

# Show plots
plt.show()
