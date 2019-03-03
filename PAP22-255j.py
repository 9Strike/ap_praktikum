import measure as ms
from measure import pi, sqrt, sin, tan, cos, arcsin, arccos, exp, ln
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
rho_NaCl = 2.164 * cs.gram / cs.centi**3
M_NaCl = 58.44 * cs.gram

# (1) Analysis of the spectrum of the LiF-crystal
# Determination of planck's constant
U1 = 35.0 * cs.kilo
t1 = 5.0
beta1, n1 = np.loadtxt('data/255/data1.txt', unpack=True)
d_n1 = sqrt(n1 * t1) / t1
n1_0 = ms.mv(n1[0:7])
d_n1_0 = ms.dsto_mv(n1[0:7])

ms.pltext.initplot(num=1, xlabel=r'$\beta$ / °', ylabel=r'$n$ / (1/s)')
s1, d_s1, b1, d_b1 = ms.linreg(beta1[:20], n1[:20], d_n1[:20], fit_range=range(10, 13), plot=True)
beta1_G = (n1_0 - b1) / s1
d_beta1_G = beta1_G * sqrt((d_n1_0**2 + d_b1**2) / (n1_0 - b1)**2 + (d_s1 / s1)**2)
beta1_G *= cs.degree
d_beta1_G *= cs.degree
ld1_G = 2 * d_LiF * sin(beta1_G)
d_ld1_G = 2 * d_LiF * cos(beta1_G) * d_beta1_G
h1 = (cs.e * U1 / cs.c) * ld1_G
d_h1 = (cs.e * U1 / cs.c) * d_ld1_G
beta1_G2 = arcsin(ld1_G / d_LiF)
d_beta1_G2 = d_ld1_G / sqrt(d_LiF**2 - ld1_G**2)

print()
print(ms.val("n0", n1_0, d_n1_0, unit='1/s', prefix=False))
print(ms.val("s1", s1, d_s1))
print(ms.val("b1", b1, d_b1))
print(ms.val("β_G", beta1_G / cs.degree, d_beta1_G / cs.degree, unit='°', prefix=False))
print(ms.val("λ_G", ld1_G, d_ld1_G, unit='m'))
print(ms.val("h", h1, d_h1, unit='Js', prefix=False))
print(ms.sig("h,l", h1, d_h1, cs.h))
print(ms.val("β_G2", beta1_G2 / cs.degree, d_beta1_G2 / cs.degree, unit='°', prefix=False))
print()

# (p1) Analysis of the K_α, K_β peaks in first and second order
def gauss(x, mu, sigma, A):
  return A / sqrt(2 * pi * sigma**2) * exp(-(x - mu)**2 / (2 * sigma**2))

t_p1 = 20.0
beta1_p1, n1_p1 = np.loadtxt('data/255/data2.txt', unpack=True)
beta2_p1, n2_p1 = np.loadtxt('data/255/data3.txt', unpack=True)
beta3_p1, n3_p1 = np.loadtxt('data/255/data4.txt', unpack=True)
beta4_p1, n4_p1 = np.loadtxt('data/255/data5.txt', unpack=True)
d_n1_p1 = sqrt(n1_p1 * t_p1) / t_p1
d_n2_p1 = sqrt(n2_p1 * t_p1) / t_p1
d_n3_p1 = sqrt(n3_p1 * t_p1) / t_p1
d_n4_p1 = sqrt(n4_p1 * t_p1) / t_p1

ms.pltext.initplot(num=2, nrows=2, ncols=2, xlabel=r'$\beta$ / $^\circ$', ylabel=r'$n$ / (1/s)')
ms.pltext.set_axis(0)
[mu1_p1, sigma1_p1, A1_p1], [d_mu1_p1, d_sigma1_p1, d_A1_p1] = ms.fit(beta1_p1, n1_p1, d_n1_p1, gauss, p0=[9.0, 0.2, 450], plot=True, fit_range=range(3, 7))
ms.pltext.set_axis(1)
[mu2_p1, sigma2_p1, A2_p1], [d_mu2_p1, d_sigma2_p1, d_A2_p1] = ms.fit(beta2_p1, n2_p1, d_n2_p1, gauss, p0=[10.15, 0.15, 500], plot=True, fit_range=range(3, 7))
ms.pltext.set_axis(2)
[mu3_p1, sigma3_p1, A3_p1], [d_mu3_p1, d_sigma3_p1, d_A3_p1] = ms.fit(beta3_p1, n3_p1, d_n3_p1, gauss, p0=[18.3, 0.2, 50], plot=True, fit_range=range(5, 9))
ms.pltext.set_axis(3)
[mu4_p1, sigma4_p1, A4_p1], [d_mu4_p1, d_sigma4_p1, d_A4_p1] = ms.fit(beta4_p1, n4_p1, d_n4_p1, gauss, p0=[20.7, 0.15, 100], plot=True, fit_range=range(4, 9))
ld1_p1 = 2 * d_LiF * sin(mu1_p1 * cs.degree)
d_ld1_p1 = 2 * d_LiF * cos(mu1_p1 * cs.degree) * d_mu1_p1 * cs.degree
ld2_p1 = 2 * d_LiF * sin(mu2_p1 * cs.degree)
d_ld2_p1 = 2 * d_LiF * cos(mu2_p1 * cs.degree) * d_mu2_p1 * cs.degree
ld3_p1 = 2 * d_LiF * sin(mu3_p1 * cs.degree)
d_ld3_p1 = 2 * d_LiF * cos(mu3_p1 * cs.degree) * d_mu3_p1 * cs.degree
ld4_p1 = 2 * d_LiF * sin(mu4_p1 * cs.degree)
d_ld4_p1 = 2 * d_LiF * cos(mu4_p1 * cs.degree) * d_mu4_p1 * cs.degree

print(ms.val("β1", mu1_p1, d_mu1_p1, unit='°', prefix=False))
print(ms.val("β2", mu2_p1, d_mu2_p1, unit='°', prefix=False))
print(ms.val("β3", mu3_p1, d_mu3_p1, unit='°', prefix=False))
print(ms.val("β4", mu4_p1, d_mu4_p1, unit='°', prefix=False))
print(ms.val("Δβ_2", 2 * sqrt(2 * ln(2)) * sigma2_p1, 2 * sqrt(2 * ln(2)) * d_sigma2_p1, unit='°', prefix=False))
print() 
print(ms.val("λ1", ld1_p1, d_ld1_p1, unit='m'))
print(ms.val("λ2", ld2_p1, d_ld2_p1, unit='m'))
print(ms.val("λ3", ld3_p1, d_ld3_p1, unit='m'))
print(ms.val("λ4", ld4_p1, d_ld4_p1, unit='m'))
print()
print(ms.sig("λ_α", ld2_p1, d_ld2_p1, 71.1e-12, perc=True))
print(ms.sig("λ_β", ld1_p1, d_ld1_p1, 63.1e-12, perc=True))
print()

# Counting rate - Voltage dependency Measurement
t = 20
beta = 7.5 * cs.degree
d_beta = 0.1 * cs.degree
U = np.arange(20.0, 36.0, 1.0) * cs.kilo
n = npf([1.35, 1.35, 2.75, 5.55, 32.95, 78.35, 122.8, 163.3, 200.6, 237.0, 270.2, 307.6, 337.1, 374.7, 403.7, 433.3])
d_n = sqrt(n * t) / t
n_0 = ms.mv(n[:3])
d_n_0 = ms.dsto_mv(n[:3])

ms.pltext.initplot(num=3, title=titles[0], xlabel=r'$U$ / V', ylabel=r'$n$ / (1/s)')
s, d_s, b, d_b = ms.linreg(U, n, d_n, fit_range=range(3, len(U)), plot=True)

U_G = (n_0 - b) / s
d_U_G = U_G * sqrt((d_n_0**2 + d_b**2) / (n_0 - b)**2 + (d_s / s)**2)
h = (2 * cs.e * d_LiF / cs.c) * sin(beta) * U_G
d_h = h * sqrt((d_U_G / U_G)**2 + (d_beta / tan(beta))**2)

print(ms.val("n0", n_0, d_n_0, unit='1/s', prefix=False))
print(ms.val("s", s, d_s))
print(ms.val("b", b, d_b))
print(ms.val("U_G", U_G, d_U_G, unit='V'))
print(ms.val("h", h, d_h, unit='Js', prefix=False))
print(ms.sig("h,l", h, d_h, cs.h))
print(ms.sig("h", h1, d_h1, h, d_h))
print()

# (2) Analysis of the spectrum of the NaCl-crystal
# Analysis of the K_α, K_β peaks in first and second order
# Determination of the lattice constant of NaCl and Avogadro's constant
U2 = 35 * cs.kilo
t2 = 5.0
beta2, n2 = np.loadtxt('data/255/data6.txt', unpack=True)
d_n2 = sqrt(n2 * t2) / t2

ms.pltext.initplot(num=4, ncols=2, nrows=2, xlabel=r'$\beta$ / $^\circ$', ylabel=r'$n$ / (1/s)')
ms.pltext.set_axis(0)
[mu1_2, sigma1_2, A1_2], [d_mu1_2, d_sigma1_2, d_A1_2] = ms.fit(beta2, n2, d_n2, gauss, p0=[6.5, 0.2, 650], fit_range=range(15, 19), plot=True)
ms.pltext.set_axis(1)
[mu2_2, sigma2_2, A2_2], [d_mu2_2, d_sigma2_2, d_A2_2] = ms.fit(beta2, n2, d_n2, gauss, p0=[7.25, 0.2, 1000], fit_range=range(20, 24), plot=True)
ms.pltext.set_axis(2)
[mu3_2, sigma3_2, A3_2], [d_mu3_2, d_sigma3_2, d_A3_2] = ms.fit(beta2, n2, d_n2, gauss, p0=[13, 0.15, 100], fit_range=range(49, 53), plot=True)
ms.pltext.set_axis(3)
[mu4_2, sigma4_2, A4_2], [d_mu4_2, d_sigma4_2, d_A4_2] =ms.fit(beta2, n2, d_n2, gauss, p0=[14.5, 0.15, 200], fit_range=range(57, 61), plot=True)

print(ms.val("μ1", mu1_2, d_mu1_2, unit='°', prefix=False))
print(ms.val("μ2", mu2_2, d_mu2_2, unit='°', prefix=False))
print(ms.val("μ3", mu3_2, d_mu3_2, unit='°', prefix=False))
print(ms.val("μ4", mu4_2, d_mu4_2, unit='°', prefix=False))
print()

beta1_2 = mu1_2 * cs.degree
d_beta1_2 = d_mu1_2 * cs.degree
beta2_2 = mu2_2 * cs.degree
d_beta2_2 = d_mu2_2 * cs.degree

l_alpha_2 = ld2_p1 / sin(beta2_2)
d_l_alpha_2 = l_alpha_2 * sqrt((d_ld2_p1 / ld2_p1)**2 + (d_beta2_2 / tan(beta2_2))**2)
l_beta_2 = ld1_p1 / sin(beta1_2)
d_l_beta_2 = l_beta_2 * sqrt((d_ld1_p1 / ld1_p1)**2 + (d_beta1_2 / tan(beta1_2))**2)

N_A_alpha_2 = 4 * M_NaCl / (rho_NaCl * l_alpha_2**3)
d_N_A_alpha_2 = N_A_alpha_2 * (3 * d_l_alpha_2 / l_alpha_2)
N_A_beta_2 = 4 * M_NaCl / (rho_NaCl * l_beta_2**3)
d_N_A_beta_2 = N_A_beta_2 * (3 * d_l_beta_2 / l_beta_2)

print(ms.val("a", l_alpha_2, d_l_alpha_2, unit='m'))
print(ms.val("a", l_beta_2, d_l_beta_2, unit='m'))
print(ms.sig("a", l_alpha_2, d_l_alpha_2, l_beta_2, d_l_beta_2))
print()
print(ms.val("N_A", N_A_alpha_2, d_N_A_alpha_2, unit='1/mol', prefix=False))
print(ms.sig("N_A,l", N_A_alpha_2, d_N_A_alpha_2, cs.N_A))
print(ms.val("N_A", N_A_beta_2, d_N_A_beta_2, unit='1/mol', prefix=False))
print(ms.sig("N_A,l", N_A_beta_2, d_N_A_beta_2, cs.N_A))
print(ms.sig("N_A", N_A_alpha_2, d_N_A_alpha_2, N_A_beta_2, d_N_A_beta_2))
print()

# Show plots
#ms.plt.show()
