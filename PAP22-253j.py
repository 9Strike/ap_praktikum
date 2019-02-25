import matplotlib.pyplot as plt
from measure import npfarray as npf
from measure import pi as π
from measure import sqrt, log10
import measure as ms
import numpy as np
import scipy.constants as cs

ms.plt.rc('text', usetex=True)
ms.plt.rc('font', family='serif')

titles = [
  r'Zählrate $n$ der $\beta$-Strahlung eines $^{90}$Sr Präparats in Abhängigkeit der Absorberdicke $d$ von Aluminium',
  r'Zählrate $n$ der $\gamma$-Strahlung eines $^{60}$Co Präparats in Abhängigkeit der Absorberdicke $d$ von Blei',
  r'Zählrate $n$ der $\alpha$-Strahlung eines $^{241}$Am Präparats in Abhängigkeit des Luftdrucks $p$'
]

# Constants
ρ_Pb = 11.34 * cs.gram / cs.centi**3

# Measurements of the counter tube
r_c = 14 * cs.milli / 2
U0_c = 520
d_U0_c = 10

# Measurement of the background
t0 = 5 * cs.minute
n0 = 122
d_n0 = sqrt(n0) / t0
n0 = n0 / t0

print(ms.val("n0", n0 * 2 * cs.minute))

# Measurement of β-Radiation absorption, Sr 90, GS 527
A_β = 74 * cs.kilo
a_β = 60 * cs.milli
d_a_β = 2 * cs.milli
t1_β = 30
t2_β = 2 * cs.minute
t3_β = 5 * cs.minute
n_β = npf([1136, 665, 412, 273, 180, 110, 306, 208, 131, 76, 60, 42])
d_n_β = sqrt(n_β)
n_β = np.append(n_β[:6] / t1_β, n_β[6:] / t2_β)
d_n_β = np.append(d_n_β[:6] / t1_β, d_n_β[6:] / t2_β)
d_β = np.arange(0.0, 0.3 * cs.milli * len(n_β), 0.3 * cs.milli)
n0_β = 51
d_n0_β = sqrt(n0_β) / t3_β
n0_β = n0_β / t3_β

n_β = n_β - n0_β
d_n_β = sqrt(d_n_β**2 + d_n0_β**2)

ms.pltext.initplot(num=1, title=titles[0], xlabel=r'$d$ / mm', ylabel=r'$\lg(n)$ / (1/s)', scale='linlog')
ms.pltext.plotdata(d_β / cs.milli, n_β, d_n_β)

# Measurement of γ-Radiation absorption, Co 60, UB 595
A_γ = 3.7 * cs.mega
a_γ = 150 * cs.milli
d_a_γ = 2 * cs.milli
t_γ = cs.minute
n_γ = npf([2447, 1760, 1279, 908, 714, 541, 412, 312, 223, 195, 164])
d_n_γ = sqrt(n_γ) / t_γ
n_γ = n_γ / t_γ
d_γ = np.arange(0.0, 5 * cs.milli * len(n_γ), 5 * cs.milli)

n_γ = n_γ - n0
d_n_γ = sqrt(d_n_γ**2 + d_n0**2)

ms.pltext.initplot(num=2, title=titles[1], xlabel=r'$d$ / mm', ylabel=r'$\lg(n)$ / (1/s)', scale='linlog')
sl_γ, d_sl_γ, _, _ = ms.expreg(d_γ / cs.milli, n_γ, d_n_γ, plot=True)
μ_γ = -sl_γ / cs.milli
d_μ_γ = d_sl_γ / cs.milli
μ_ρ_γ = μ_γ / ρ_Pb      # => E_γ = (1.45 ± 0.05) MeV
d_μ_ρ_γ = d_μ_γ / ρ_Pb

print(ms.val("μ/ρ", μ_ρ_γ * cs.gram / cs.centi**2, d_μ_ρ_γ * cs.gram / cs.centi**2, unit='cm²/g', prefix=False))
print()

# Measurement of γ-Radiation activity, Co 60, UB 595
t_Aγ = cs.minute
a_Aγ = npf([50, 105, 190]) * cs.milli
d_a_Aγ = npf([2, 2, 2]) * cs.milli
n_Aγ = npf([33865, 8266, 2171])
d_n_Aγ = sqrt(n_Aγ) / t_Aγ
n_Aγ = n_Aγ / t_Aγ

ε_Aγ = 0.04
A_Aγ = 4 * n_Aγ * a_Aγ**2 / (ε_Aγ * r_c**2)
d_A_Aγ = A_Aγ * sqrt((d_n_Aγ / (n_Aγ * a_Aγ))**2 + (2 * d_a_Aγ / a_Aγ)**2)

print(ms.tbl([ms.lst(A_Aγ, d_A_Aγ, name='A')]))
print()

# Measurement of α-Radiation absorption and energy, Am 241, AP 15.2
s_c = 4.2 * cs.centi
σ_c = 2.25 * cs.milli / cs.centi**2
A_α = 90 * cs.kilo
a_α = 3.95 * cs.centi
d_a_α = 0.05 * cs.centi
t_α = cs.minute
p1_α = npf([18, 98, 120, 225, 324, 383, 416, 450, 475, 515, 617, 721, 813, 911, 1013]) * cs.milli * cs.bar
p2_α = npf([20, 98, 120, 222, 324, 392, 413, 450, 473, 512, 614, 717, 809, 908, 1009]) * cs.milli * cs.bar
d_p_α = cs.milli * cs.bar
n_α = npf([13144, 13142, 13131, 12933, 12615, 9883, 7101, 3491, 1680, 451, 239, 205, 229, 225, 212])
d_n_α = sqrt(n_α) / t_α
n_α = n_α / t_α

p_α = npf([0.5 * (p1_α[i] + p2_α[i]) for i in range(len(p1_α))])
d_p_α = npf([np.max([0.5 * np.abs(p2_α[i] - p1_α[i]), 1]) for i in range(len(p1_α))])

ms.pltext.initplot(num=3, title=titles[2], xlabel=r'$p$ / Pa', ylabel=r'$n$ / (1/s)')
sl_α, d_sl_α, i_α, d_i_α = ms.linreg(p_α, n_α, d_n_α, d_p_α, fit_range=range(5,9), plot=True)
p_H = (n_α[0] / 2 - i_α) / sl_α
d_p_H = p_H * sqrt((d_n_α[0]**2 + 4 * d_i_α**2) / (n_α[0] - 2 * i_α)**2 + (d_sl_α / sl_α)**2)
s_α = p_H / ms.p0 * a_α + σ_c / (1.43 * cs.milli / cs.centi**2) * cs.centi + 0.68 * cs.centi      # => E_α = ( ± ) eV
d_s_α = sqrt((a_α * d_p_H)**2 + (p_H * d_a_α)**2) / ms.p0

print("Absorption of α-Radiation:")
print(ms.val("sl", sl_α, d_sl_α))
print(ms.val("i", i_α, d_i_α))
print(ms.val("p", p_H, d_p_H))
print(ms.val("s", s_α, d_s_α))

ms.plt.show()
