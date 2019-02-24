import measure as ms
from measure import sqrt, exp, ln
import numpy as np
import scipy.constants as cs
from scipy.optimize import curve_fit
from scipy.stats import chi2

ms.plt.rc('figure', figsize=(11.69, 8.27))

# Measured data
U0 = 520
d_U0 = 10
ΔT0 = 10.0
ΔT_Ag = 10.0
ΔT_In = 120.0
T0 = 8 * cs.minute
T_Ag = 400
T_In = 50 * cs.minute

n0 = np.loadtxt('data/PAP22/252_1j.dat', usecols=[1], unpack=True)
n_Ag_1 = np.loadtxt('data/PAP22/252_2j.dat', usecols=[1], unpack=True)
n_Ag_2 = np.loadtxt('data/PAP22/252_3j.dat', usecols=[1], unpack=True)
n_Ag_3 = np.loadtxt('data/PAP22/252_4j.dat', usecols=[1], unpack=True)
n_Ag_4 = np.loadtxt('data/PAP22/252_5j.dat', usecols=[1], unpack=True)
n_In = np.loadtxt('data/PAP22/252_6j.dat', usecols=[1], unpack=True)

# Background radiation
n0_m = ms.mv(n0 / ΔT0)
d_n0_m = ms.dsto_mv(n0 / ΔT0)

# Fit function
def f_Ag(x, A1, λ1, A2, λ2):
  return A1 * exp(-λ1 * x) + A2 * exp(-λ2 * x)
def f_In(x, A, λ):
  return A * exp(-λ * x)

# Ag decay
t_Ag = np.arange(ΔT_Ag / 2, T_Ag + ΔT_Ag / 2, ΔT_Ag)
N_Ag = (n_Ag_1 + n_Ag_2 + n_Ag_3 + n_Ag_4)
d_N_Ag = sqrt(N_Ag) / (4 * ΔT_Ag)
N_Ag = N_Ag / (4 * ΔT_Ag)
N_Ag = N_Ag - n0_m
d_N_Ag = sqrt(d_N_Ag**2 + d_n0_m**2)

cut = 0
for i in range(len(N_Ag)):
  if N_Ag[i] / d_N_Ag[i] < 3.0:
    cut = i
    break
[A1_Ag, λ1_Ag, A2_Ag, λ2_Ag], pcov = curve_fit(f_Ag, t_Ag[:cut], N_Ag[:cut], sigma=d_N_Ag[:cut], p0=[30, ln(2) / 24.6, 5, ln(2) / (2.41 * cs.minute)])
[d_A1_Ag, d_λ1_Ag, d_A2_Ag, d_λ2_Ag] = sqrt(np.diag(pcov))

t_Ag_fit = np.linspace(t_Ag[0], t_Ag[-1], 1000)
N_Ag_fit = f_Ag(t_Ag_fit, A1_Ag, λ1_Ag, A2_Ag, λ2_Ag)
ms.pltext.initplot(num=1, title='Mittelwert der Zerfallsrate A aus 4 Messungen von Ag 108 und Ag 110 zu der Zeit t', xlabel='t / s', ylabel='A / Bq', fignum=True)
ms.pltext.plotdata(t_Ag, N_Ag, d_N_Ag, color='gray')
ms.plt.plot(t_Ag_fit, N_Ag_fit)

dof_Ag_fit = len(N_Ag[:cut]) - 4
χ2_Ag = ms.chi2(N_Ag[:cut], d_N_Ag[:cut], f_Ag(t_Ag[:cut], A1_Ag, λ1_Ag, A2_Ag, λ2_Ag))
χ2_Ag_red = χ2_Ag / dof_Ag_fit

p_Ag_fit = 1 - chi2.cdf(χ2_Ag, dof_Ag_fit)

τ1_Ag = ln(2) / λ1_Ag
d_τ1_Ag = ln(2) * d_λ1_Ag / λ1_Ag**2
τ2_Ag = ln(2) / λ2_Ag
d_τ2_Ag = ln(2) * d_λ2_Ag / λ2_Ag**2

print("Ag108, Ag110 decay:")
print(ms.val("A1", A1_Ag, d_A1_Ag, unit='Bq'))
print(ms.val("λ1", λ1_Ag, d_λ1_Ag, unit='1/s', prefix=False))
print(ms.val("A2", A2_Ag, d_A2_Ag, unit='Bq'))
print(ms.val("λ2", λ2_Ag, d_λ2_Ag, unit='1/s', prefix=False))
print(ms.val("χ²", χ2_Ag))
print(ms.val("χr²", χ2_Ag_red))
print(ms.val("pfit", 100 * p_Ag_fit, unit='%'))
print()
print("Half life of Ag110:")
print(ms.val("τ", τ1_Ag, d_τ1_Ag, unit='s'))
print(ms.sig("τ", τ1_Ag, d_τ1_Ag, 24.6))
print("Half life of Ag108:")
print(ms.val("τ", τ2_Ag / cs.minute, d_τ2_Ag / cs.minute, unit='min'))
print(ms.sig("τ", τ2_Ag, d_τ2_Ag, 2.37 * cs.minute))
print()

# In decay
t_In = np.arange(ΔT_In / 2, T_In + ΔT_In / 2, ΔT_In)
d_n_In = sqrt(n_In) / ΔT_In
n_In = n_In / ΔT_In
n_In = n_In - n0_m
d_n_In = sqrt(d_n_In**2 + d_n0_m**2)

[A_In, λ_In], pcov = curve_fit(f_In, t_In[1:], n_In[1:], sigma=d_n_In[1:], p0=[4, 0.5e-3])
[d_A_In, d_λ_In] = sqrt(np.diag(pcov))

t_In_fit = np.linspace(t_In[0], t_In[-1], 1000)
n_In_fit = f_In(t_In_fit, A_In, λ_In)
ms.pltext.initplot(num=2, title='Zerfallsrate A von In 116m (und In 116) in Abhängigkeit zu der Zeit t.', xlabel='t / s', ylabel='A / Bq', fignum=True)
ms.pltext.plotdata(t_In, n_In, d_n_In, color='gray')
ms.plt.plot(t_In_fit, n_In_fit)

dof_In_fit = len(n_In[1:]) - 2
χ2_In = ms.chi2(n_In[1:], d_n_In[1:], f_In(t_In[1:], A_In, λ_In))
χ2_In_red = χ2_In / (len(n_In[1:]) - 2)

p_In_fit = 1 - chi2.cdf(χ2_In, dof_In_fit)

τ_In = ln(2) / λ_In
d_τ_In = ln(2) * d_λ_In / λ_In**2

print("In decay:")
print(ms.val("A", A_In, d_A_In, unit='Bq'))
print(ms.val("λ", λ_In, d_λ_In, unit='1/s', prefix=False))
print(ms.val("χ²", χ2_In))
print(ms.val("χr²", χ2_In_red))
print(ms.val("pfit", 100 * p_In_fit, unit='%'))
print()
print("Half life of In116:")
print(ms.val("τ", τ_In / cs.minute, d_τ_In / cs.minute, unit='min'))
print(ms.sig("τ", τ_In, d_τ_In, 54.29 * cs.minute))

ms.pltext.savefigs('figures')

ms.plt.show()
