import measure as ms
from measure import npfarray as npf
import numpy as np
from numpy import pi as π
from numpy import sqrt, exp
import scipy.constants as sc
from scipy.optimize import curve_fit
from scipy.stats import chi2
from scipy.special import gamma

# Measured values
# Measurement of the plateau-sector
V_E = 440
t0 = 30.0
U = npf([500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750])
n = npf([2462, 2606, 2636, 2738, 2753, 2671, 2618, 2685, 2742, 2715, 2792])
d_n = sqrt(n) / t0
n = n / t0

# Determination of the plateau-sector
ms.pltext.initplot(num=1, title='Plateaubereich, Abhängigkeit der Zählrate n von der Zählrohrspannung U', xlabel='U/V', ylabel='n', fignum=True)
[s, d_s, i, d_i] = ms.linreg(U, n, d_n, fit_range=range(1, len(U)), plot=True)

U_0 = 600

# Measurement of the plateau-slope
t1 = sc.minute
U1 = npf([U_0, U_0 + 100])
n1 = npf([12022, 12172])
d_n1 = sqrt(n1) / t1
n1 = n1 / t1

t2 = 3 * sc.minute
U2 = npf([U_0, U_0 + 100])
n2 = npf([35870, 36579])
d_n2 = sqrt(n2) / t2
n2 = n2 / t2

# Determination of the plateau-slope
s1 = n1[1] - n1[0]
d_s1 = sqrt(d_n1[1]**2 + d_n1[0]**2)
s1_p = s1 / n1[0]
d_s1_p = s1_p * sqrt((d_s1 / s1)**2 + (d_n1[0] / n1[0])**2)

T1 = 1/0.01 * (n1[1] + n1[0]) / (n1[1] - n1[0])

s2 = n2[1] - n2[0]
d_s2 = sqrt(d_n2[1]**2 + d_n2[0]**2)
s2_p = s2 / n2[0]
d_s2_p = s2_p * sqrt((d_s2 / s2)**2 + (d_n2[0] / n2[0])**2)

T2 = 1/0.01 * (n2[1] + n2[0]) / (n2[1] - n2[0])

print()
print("Plateau-slope:")
print(ms.val('s1', s1, d_s1))
print(ms.val('s1_p', s1_p * 100, d_s1_p * 100, unit='%', prefix=False))
print()
print(ms.sig('s1,0', s1, d_s1, 0.0))
print(ms.val('T1', T1 / sc.hour, unit='hours', prefix=False))
print(ms.val('s1_p + Δs1_p', (s1_p + d_s1_p) * 100, unit='%', prefix=False))
print(ms.val('s1_p + 2 Δs1_p', (s1_p + 2 * d_s1_p) * 100, unit='%', prefix=False))
print()
print(ms.val('s2', s2, d_s2))
print(ms.val('s2_p', s2_p * 100, d_s2_p * 100, unit='%', prefix=False))
print()
print(ms.sig('s2,0', s2, d_s2, 0.0))
print(ms.val('T2', T2 / sc.hour, unit='hours', prefix=False))
print(ms.val('s2_p + Δs2_p', (s2_p + d_s2_p) * 100, unit='%', prefix=False))
print(ms.val('s2_p + 2 Δs2_p', (s2_p + 2 * d_s2_p) * 100, unit='%', prefix=False))
print()


# Propability distributions
def g_pdf(x, μ, σ, A):
  return A / sqrt(2 * π * σ**2) * exp(-(x - μ)**2 / (2 * σ**2))
def p_pdf(k, μ, A):
  return A * (μ**k / gamma(k + 1)) * exp(-μ)

# 1.Determination of the statistics
n3, f3 = np.loadtxt('data/PAP22/251_1j.dat', unpack=True)
d_f3 = sqrt(f3)

[μ3_g, σ3_g, A3_g], pcov_g = curve_fit(g_pdf, n3[1:-2], f3[1:-2], sigma=d_f3[1:-2], p0=[25, 5, 2000])
[μ3_p, A3_p], pcov_p = curve_fit(p_pdf, n3[1:-2], f3[1:-2], sigma=d_f3[1:-2], p0=[25, 2000])
[d_μ3_g, d_σ3_g, d_A3_g] = sqrt(np.diag(pcov_g))
[d_μ3_p, d_A3_p] = sqrt(np.diag(pcov_p))
n3_gp = np.linspace(n3[0], n3[-1], 1000)
f3_g = g_pdf(n3_gp, μ3_g, σ3_g, A3_g)
f3_p = p_pdf(n3_gp, μ3_p, A3_p)

dof_g3 = len(n3[1:-2]) - 3
dof_p3 = len(n3[1:-2]) - 2
χ2_g3 = ms.chi2(f3[1:-2], d_f3[1:-2], g_pdf(n3[1:-2], μ3_g, σ3_g, A3_g))
χ2_p3 = ms.chi2(f3[1:-2], d_f3[1:-2], p_pdf(n3[1:-2], μ3_p, A3_p))
χ2_red_g3 = χ2_g3 / dof_g3
χ2_red_p3 = χ2_p3 / dof_p3
fitprob_g3 = 1 - chi2.cdf(χ2_g3, dof_g3)
fitprob_p3 = 1 - chi2.cdf(χ2_p3, dof_p3)

print("Statistics 1:")
print("Gaussian fit:")
print(ms.val("μ", μ3_g, d_μ3_g))
print(ms.val("σ", σ3_g, d_σ3_g))
print(ms.val("A", A3_g, d_A3_g))
print(ms.val("χ²", χ2_g3))
print(ms.val("χ_r²", χ2_red_g3))
print(ms.val("p_fit", round(fitprob_g3, 2) * 100, unit='%', prefix=False))
print("Poisson fit:")
print(ms.val("μ", μ3_p, d_μ3_p))
print(ms.val("A", A3_p, d_A3_p))
print(ms.val("χ²", χ2_p3))
print(ms.val("χ_r²", χ2_red_p3))
print(ms.val("p_fit", round(fitprob_p3, 2) * 100, unit='%', prefix=False))
print()

ms.pltext.initplot(num=2, title='Häufigkeit f der Zählrate n bei 2200 Messungen der Zählraten von Co60 mit einer Messzeit von 0,5s', xlabel='n/(1/s)', ylabel='f', fignum=True)
ms.pltext.plotdata(n3, f3, d_f3, color='gray')
ms.plt.plot(n3_gp, f3_g)
ms.plt.plot(n3_gp, f3_p)
ms.plt.legend(["Gauss fit", "Poisson fit"])

# 2.Determination of the statistics
n4, f4 = np.loadtxt('data/PAP22/251_2j.dat', unpack=True)
d_f4 = sqrt(f4)

[μ4_g, σ4_g, A4_g], _ = curve_fit(g_pdf, n4[:-2], f4[:-2], sigma=d_f4[:-2], p0=[4, 1.5, 4000])
[μ4_p, A4_p], _ = curve_fit(p_pdf, n4[:-2], f4[:-2], sigma=d_f4[:-2], p0=[4, 4000])
[d_μ4_g, d_σ4_g, d_A4_g] = sqrt(np.diag(pcov_g))
[d_μ4_p, d_A4_p] = sqrt(np.diag(pcov_p))
n4_gp = np.linspace(n4[0], n4[-1], 1000)
f4_g = g_pdf(n4_gp, μ4_g, σ4_g, A4_g)
f4_p = p_pdf(n4_gp, μ4_p, A4_p)

dof_g4 = len(n4[:-2]) - 3
dof_p4 = len(n4[:-2]) - 2
χ2_g4 = ms.chi2(f4[:-2], d_f4[:-2], g_pdf(n4[:-2], μ4_g, σ4_g, A4_g))
χ2_p4 = ms.chi2(f4[:-2], d_f4[:-2], p_pdf(n4[:-2], μ4_p, A4_p))
χ2_red_g4 = χ2_g4 / dof_g4
χ2_red_p4 = χ2_p4 / dof_p4
fitprob_g4 = 1 - chi2.cdf(χ2_g4, dof_g4)
fitprob_p4 = 1 - chi2.cdf(χ2_p4, dof_p4)

print("Statistics 2:")
print("Gaussian fit:")
print(ms.val("μ", μ4_g, d_μ4_g))
print(ms.val("σ", σ4_g, d_σ4_g))
print(ms.val("A", A4_g, d_A4_g))
print(ms.val("χ²", χ2_g4))
print(ms.val("χ_r²", χ2_red_g4))
print(ms.val("p_fit", round(fitprob_g4, 2) * 100, unit='%', prefix=False))
print("Poisson fit:")
print(ms.val("μ", μ4_p, d_μ4_p))
print(ms.val("A", A4_p, d_A4_p))
print(ms.val("χ²", χ2_p4))
print(ms.val("χ_r²", χ2_red_p4))
print(ms.val("p_fit", round(fitprob_p4, 2) * 100, unit='%', prefix=False))

ms.pltext.initplot(num=3, title='Häufigkeit f der Zählrate n bei 5100 Messungen der Zählraten von Co60 mit einer Messzeit von 0,1s', xlabel='n/(1/s)', ylabel='lg(f)', scale='linlog', fignum=True)
ms.pltext.plotdata(n4, f4, d_f4, color='gray')
ms.plt.plot(n4_gp, f4_g)
ms.plt.plot(n4_gp, f4_p)
ms.plt.legend(["Gauss fit", "Poisson fit"])

ms.pltext.savefigs('figures')

ms.plt.show()
