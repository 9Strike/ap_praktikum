import datetime as dt
import matplotlib.pyplot as plt
import measure as ms
from measure import npfarray as npf
from measure import pi
from measure import sqrt, exp, log10, ln
import numpy as np
import scipy.constants as cs
from scipy.optimize import curve_fit

ms.plt.rc('text', usetex=True)
ms.plt.rc('font', family='serif')

titles = [
  r'Zählrate $n$ der $\beta$-Strahlung eines $^{90}$Sr Präparats in Abhängigkeit der Absorberdicke $d$ von Aluminium',
  r'Zählrate $n$ der $\gamma$-Strahlung eines $^{60}$Co Präparats in Abhängigkeit der Absorberdicke $d$ von Blei',
  r'Zählrate $n$ der $\alpha$-Strahlung eines $^{241}$Am Präparats in Abhängigkeit des Luftdrucks $p$'
]

# Constants
rho_Pb = 11.34 * cs.gram / cs.centi**3
rho_Al = 2.699 * cs.gram / cs.centi**3

# Measurements of the counter tube
r_c = 14 * cs.milli / 2
l_c = 4 * cs.centi
U0_c = 520
d_U0_c = 10

# Measurement of the background
t0 = 5 * cs.minute
n0 = 122
d_n0 = sqrt(n0) / t0
n0 = n0 / t0

print(ms.val("n0", n0, d_n0, unit='1/s', prefix=False))

# Measurement of β-Radiation absorption, Sr 90, GS 527
A_Sr = 74 * cs.kilo
a_Sr = 60 * cs.milli
d_a_Sr = 2 * cs.milli
t1_Sr = 30
t2_Sr = 2 * cs.minute
t3_Sr = 5 * cs.minute
n_Sr = npf([1136, 665, 412, 273, 180, 110, 306, 208, 131, 76, 60, 42])
d_n_Sr = sqrt(n_Sr)
n_Sr = np.append(n_Sr[:6] / t1_Sr, n_Sr[6:] / t2_Sr)
d_n_Sr = np.append(d_n_Sr[:6] / t1_Sr, d_n_Sr[6:] / t2_Sr)
d_Sr = np.arange(0.0, 0.3 * cs.milli * len(n_Sr), 0.3 * cs.milli)
n0_Sr = 51
d_n0_Sr = sqrt(n0_Sr) / t3_Sr
n0_Sr = n0_Sr / t3_Sr

n_Sr -= n0_Sr
d_n_Sr = sqrt(d_n_Sr**2 + d_n0_Sr**2)

ms.pltext.initplot(num=1, title=titles[0], xlabel=r'$d$ / mm', ylabel=r'$\lg(n)$ / (1/s)', scale='linlog', fignum=True)
ms.expreg(d_Sr / cs.milli, n_Sr, d_n_Sr, fit_range=range(len(d_Sr) - 3), plot=True)
ms.pltext.initplot(num=2, title=titles[0], xlabel=r'$d$ / mm', ylabel=r'$n$ / (1/s)', scale='linlin', fignum=True)
sl_Sr, d_sl_Sr, i_Sr, d_i_Sr = ms.linreg(d_Sr / cs.milli, n_Sr, d_n_Sr, fit_range=range(len(d_Sr) - 3, len(d_Sr)), plot=True)
sl_Sr /= cs.milli
d_sl_Sr /= cs.milli

R_Sr = -i_Sr / sl_Sr
d_R_Sr = R_Sr * sqrt((d_i_Sr / i_Sr)**2 + (d_sl_Sr / sl_Sr)**2)

sigma_ES = 0.130 * cs.gram / cs.centi**2
sigma_Sr = R_Sr * rho_Al + sigma_ES      # => E_Sr = (2.3 ± 0.8) MeV
d_sigma_Sr = d_R_Sr * rho_Al

print("Absorption of β-Radiation")
print(ms.val("sl", sl_Sr * cs.milli, d_sl_Sr * cs.milli, unit='1/(mm s)', prefix=False))
print(ms.val("i", i_Sr, d_i_Sr, unit='1/s', prefix=False))
print(ms.val("R", R_Sr, d_R_Sr, unit='m'))
print(ms.val("σ", sigma_Sr * cs.centi**2 / cs.gram, d_sigma_Sr * cs.centi**2 / cs.gram, unit='g/cm²'))
print(ms.val("E", 2.3e6, 0.8e6, unit='eV'))
print()

# Measurement of γ-Radiation absorption, Co 60, UB 595
a_Co = 150 * cs.milli
d_a_Co = 2 * cs.milli
t_Co = cs.minute
n_Co = npf([2447, 1760, 1279, 908, 714, 541, 412, 312, 223, 195, 164])
d_n_Co = sqrt(n_Co) / t_Co
n_Co = n_Co / t_Co
d_Co = np.arange(0.0, 5 * cs.milli * len(n_Co), 5 * cs.milli)

n_Co = n_Co - n0
d_n_Co = sqrt(d_n_Co**2 + d_n0**2)

ms.pltext.initplot(num=3, title=titles[1], xlabel=r'$d$ / mm', ylabel=r'$\lg(n)$ / (1/s)', scale='linlog', fignum=True)
sl_Co, d_sl_Co, _, _ = ms.expreg(d_Co / cs.milli, n_Co, d_n_Co, plot=True)
mu_Co = -sl_Co / cs.milli
d_mu_Co = d_sl_Co / cs.milli
mu_rho_Co = mu_Co / rho_Pb      # => E_Co = (1.45 ± 0.05) MeV
d_mu_rho_Co = d_mu_Co / rho_Pb

print("Absorption of γ-Radiation")
print(ms.val("μ", mu_Co * cs.milli, d_mu_Co * cs.milli, unit='1/mm', prefix=False))
print(ms.val("μ/ρ", mu_rho_Co * cs.gram / cs.centi**2, d_mu_rho_Co * cs.gram / cs.centi**2, unit='cm²/g', prefix=False))
print(ms.val("E", 1.45e6, 0.05e6, unit='eV'))
print()

# Measurement of γ-Radiation activity, Co 60, UB 595
A_N_CoA = 3.7 * cs.mega
T_CoA = (dt.datetime(2019, 2, 21) - dt.datetime(2012, 2, 2)).total_seconds()
T_H_CoA = 5.27 * cs.year
eps_CoA = 0.04
rho_abs_CoA = 7.9 * cs.gram / cs.centi**3
d_abs_CoA = 1.4 * cs.milli
t_CoA = cs.minute
a_CoA = npf([50, 105, 190]) * cs.milli
d_a_CoA = npf([2, 2, 2]) * cs.milli
n_CoA = npf([33865, 8266, 2171])
d_n_CoA = sqrt(n_CoA) / t_CoA
n_CoA = n_CoA / t_CoA

mu_abs_CoA = mu_rho_Co * rho_abs_CoA
d_mu_abs_CoA = d_mu_rho_Co * rho_abs_CoA
A1_CoA = 4 * n_CoA * a_CoA**2 / (eps_CoA * r_c**2)
d_A1_CoA = A1_CoA * sqrt((d_n_CoA / (n_CoA * a_CoA))**2 + (2 * d_a_CoA / a_CoA)**2)
A2_CoA = 4 * n_CoA * (a_CoA + l_c / 2)**2 / (eps_CoA * r_c**2)
d_A2_CoA = A2_CoA * sqrt((d_n_CoA / (n_CoA * a_CoA))**2 + (2 * d_a_CoA / (a_CoA + l_c / 2))**2)
A3_CoA = A2_CoA * exp(-mu_abs_CoA * d_abs_CoA)
d_A3_CoA = A3_CoA * sqrt((d_A2_CoA / A2_CoA)**2 + (d_abs_CoA * d_mu_abs_CoA)**2)
A_l_CoA = A_N_CoA * exp(-ln(2) * T_CoA / T_H_CoA)

print(ms.tbl([ms.lst(A1_CoA, d_A1_CoA, name='A', unit='Bq'),
              ms.lst(A2_CoA, d_A2_CoA, name='A', unit='Bq'),
              ms.lst(A3_CoA, d_A3_CoA, name='A', unit='Bq')]))
print(ms.val("A", A_l_CoA))
print(ms.tbl([ms.dev(A1_CoA, d_A1_CoA, [A_l_CoA] * 3, name='A'),
              ms.dev(A2_CoA, d_A2_CoA, [A_l_CoA] * 3, name='A'),
              ms.dev(A3_CoA, d_A3_CoA, [A_l_CoA] * 3, name='A')]))

# Measurement of Am-Radiation absorption and energy, Am 241, AP 15.2
s_c = 4.2 * cs.centi
sigma_c = 2.25 * cs.milli / cs.centi**2
A_Am = 90 * cs.kilo
a_Am = 3.95 * cs.centi
d_a_Am = 0.05 * cs.centi
t_Am = cs.minute
p1_Am = npf([18, 98, 120, 225, 324, 383, 416, 450, 475, 515, 617, 721, 813, 911, 1013]) * cs.milli * cs.bar
p2_Am = npf([20, 98, 120, 222, 324, 392, 413, 450, 473, 512, 614, 717, 809, 908, 1009]) * cs.milli * cs.bar
d_p_Am = cs.milli * cs.bar
n_Am = npf([13144, 13142, 13131, 12933, 12615, 9883, 7101, 3491, 1680, 451, 239, 205, 229, 225, 212])
d_n_Am = sqrt(n_Am) / t_Am
n_Am = n_Am / t_Am

p_Am = npf([0.5 * (p1_Am[i] + p2_Am[i]) for i in range(len(p1_Am))])
d_p_Am = npf([np.max([0.5 * np.abs(p2_Am[i] - p1_Am[i]), 1]) for i in range(len(p1_Am))])

ms.pltext.initplot(num=4, title=titles[2], xlabel=r'$p$ / Pa', ylabel=r'$n$ / (1/s)', fignum=True)
sl_Am, d_sl_Am, i_Am, d_i_Am = ms.linreg(p_Am, n_Am, d_n_Am, d_p_Am, fit_range=range(5,9), plot=True)
p_H = (n_Am[0] / 2 - i_Am) / sl_Am
d_p_H = p_H * sqrt((d_n_Am[0]**2 + 4 * d_i_Am**2) / (n_Am[0] - 2 * i_Am)**2 + (d_sl_Am / sl_Am)**2)
s_Am = p_H / ms.p0 * a_Am + sigma_c / (1.43 * cs.milli / cs.centi**2) * cs.centi + 0.68 * cs.centi      # => E_Am = (6.0 ± 0.5) MeV
d_s_Am = sqrt((a_Am * d_p_H)**2 + (p_H * d_a_Am)**2) / ms.p0

print("Absorption of α-Radiation:")
print(ms.val("sl", sl_Am, d_sl_Am, unit='1/(s Pa)', prefix=False))
print(ms.val("i", i_Am, d_i_Am, unit='1/s', prefix=False))
print(ms.val("p", p_H, d_p_H, unit='Pa'))
print(ms.val("s", s_Am, d_s_Am, unit='m'))

ms.plt.show()
