import measure as ms
from measure import pi, sqrt, sin, cos, exp
from measure import npfarray as npf
import numpy as np
import scipy.constants as cs

ms.plt.rc('text', usetex=True)
ms.plt.rc('font', family='serif')
conv = [lambda s: float(s.replace(',', '.'))]

titles = [
  r'Zählrate $n$ in Abhängigkeit der Spannung der Röntgenröhre $U$ bei einem konstanten Winkel von 7.5$^\circ$'
]

# Constants
d_LiF = 201.4 * cs.pico

# Analysis of the spectrum of the LiF-crystal
# Determination of planck's constant
U1 = 35.0 * cs.kilo
t1 = 5.0
beta1, n1 = np.loadtxt('data/255/data1.txt', unpack=True)
d_n1 = sqrt(n1 * t1) / t1
n0 = ms.mv(n1[0:7])
d_n0 = ms.dsto_mv(n1[0:7])

ms.pltext.initplot(num=1, xlabel=r'$\beta$ / °', ylabel=r'$n$ / (1/s)')
s1, d_s1, i1, d_i1 = ms.linreg(beta1[:20], n1[:20], d_n1[:20], fit_range=range(10, 13), plot=True)
beta1_G = (n0 - i1) / s1
d_beta1_G = beta1_G * sqrt((d_n0**2 + d_i1**2) / (n0 - i1)**2 + (d_s1 / s1)**2)
beta1_G *= cs.degree
d_beta1_G *= cs.degree
ld1_G = 2 * d_LiF * sin(beta1_G)
d_ld1_G = 2 * d_LiF * cos(beta1_G) * d_beta1_G
h = (cs.e * U1 / cs.c) * ld1_G
d_h = (cs.e * U1 / cs.c) * d_ld1_G

print(ms.val("h", h, d_h, unit='Js', prefix=False))
print(ms.sig("h", h, d_h, cs.h))

# Analysis of the K_α, K_β peaks in first and second order
def gauss(x, mu, sigma, A, B):
  return A / sqrt(2 * pi * sigma**2) * exp(-(x - mu)**2 / (2 * sigma**2)) + B

tp = 20.0
betap1, np1 = np.loadtxt('data/255/data2.txt', unpack=True)
betap2, np2 = np.loadtxt('data/255/data3.txt', unpack=True)
betap3, np3 = np.loadtxt('data/255/data4.txt', unpack=True)
betap4, np4 = np.loadtxt('data/255/data5.txt', unpack=True)
d_np1 = sqrt(np1 * tp) / tp
d_np2 = sqrt(np2 * tp) / tp
d_np3 = sqrt(np3 * tp) / tp
d_np4 = sqrt(np4 * tp) / tp

ms.pltext.initplot(num=2, xlabel=r'$\beta / °$', ylabel=r'$n$ / (1/s)')
ms.fit(betap1, np1, d_np1, gauss, p0=[9.0, 0.15, 2000, 375], plot=True)

ms.pltext.initplot(num=3, xlabel=r'$\beta / °$', ylabel=r'$n$ / (1/s)')
ms.fit(betap2, np2, d_np2, gauss, p0=[9.0, 0.15, 2000, 375], plot=True)

ms.pltext.initplot(num=4, xlabel=r'$\beta / °$', ylabel=r'$n$ / (1/s)')
ms.fit(betap3, np3, d_np3, gauss, p0=[9.0, 0.15, 2000, 375], plot=True)

ms.pltext.initplot(num=5, xlabel=r'$\beta / °$', ylabel=r'$n$ / (1/s)')
ms.fit(betap4, np4, d_np4, gauss, p0=[9.0, 0.15, 2000, 375], plot=True)

# Counting rate - Voltage dependency Measurement
t = 20
φ = 7.5 * cs.degree
U = np.arange(20.0, 36.0, 1.0) * cs.kilo
n = npf([1.35, 1.35, 2.75, 5.55, 32.95, 78.35, 122.8, 163.3, 200.6, 237.0, 270.2, 307.6, 337.1, 374.7, 403.7, 433.3])
d_n = sqrt(n * t) / t

#ms.pltext.initplot(num=2, title=titles[0], xlabel=r'$U$ / V', ylabel=r'$n$ / (1/s)')
#s, d_s, i, d_i = ms.linreg(U, n, d_n, fit_range=range(3, len(U)), plot=True)

ms.plt.show()
