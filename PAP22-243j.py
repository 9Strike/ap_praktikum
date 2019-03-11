import measure as ms
import numpy as np
import scipy.constants as cs

from measure import npfarray as fa
from measure import sqrt
from scipy.integrate import quad

ms.plt.rc('text', usetex=True)
ms.plt.rc('font', family='serif')

titles = [
  r'Differenz der Quadrate der effektiven Rauschspannung $U_N$ und der zusätzlichen Rauschspannung $U_V$ des Verstärkers' '\n'
    r' in Abhängigkeit des Widerstandes $R$',
  r'Frequenzgang des Messsystems.'
]

# Measurement of the noise voltage in dependency of the resistance
rd_V = 0.003
rd_R = 0.005 # maximal

# (1) Measurement of the noise voltage in dependency of the resistance
R1 = np.arange(5, 35, 5) * cs.kilo
U1 = fa([])
d_U1 = fa([]) / sqrt(100)
d_U1 = sqrt(d_U1**2 + (rd_V * U1)**2)
T1 =  + cs.zero_Celsius
d_T1 =  + cs.zero_Celsius

# (2) Measurement of the inherent noise voltage of the multiplier
U2 = 
d_U2 =  / sqrt(100)
d_U2 = sqrt(d_U1**2 + (rd_V * U1)**2)

U_diff = U1**2 - U2**2
d_U_diff = sqrt((2 * U1 * d_U1)**2 + (2 * U2 * d_U2)**2)

ms.initplot(num=1, title=titles[0], xlabel=r'$R$ / k$\Omega$', ylabel=r'$(U_N**2 - U_V**2)$ / mV$**2$', scale='loglog', fignum=True)
[s, d_s, b, d_b] = ms.linreg(R1, U_diff, d_U_diff, plot=True)

# Measurement of the measurement system's frequency slope
D = 0.001
d_D = 0.002 * D

def g_fit(f, V, Omg1, Omg2, n1, n2):
  return V / sqrt((1 + (Omg1 / f)**(2 * n1)) * (1 + (f / Omg2)**(2 * n2)))
def g2_fit(f, V, Omg1, Omg2, n1, n2):
  return g(f, V, Omg1, Omg2, n1, n2)**2

# (3) Measurement of the measurement system's frequency slope and the equivalent noise bandwidth
U3_E = 0.2
f3, U3_A = np.loadtxt('data/243/data2.dat', skiprows=1, usecols=(0, 1), unpack=True)
# d_f3, d_U3_A
g3 = U3_A / (D * U3_E)
d_g3 = g3 * d_D / D

ms.pltext.initplot(num=2, title=titles[1], xlabel=r'$f$ / Hz', ylabel=r'$g(f)$', scale='loglog', fignum=True)
param3, d_param = ms.fit(f3, g3, d_g3, g_fit, p0=[1000, 1000, 50000, 5, 5], fit_range=range(len(f3)), plot=True)
[V3, Omg3_1, Omg3_2, n3_1, n3_2] = param3
[d_V3, d_Omg3_1, d_Omg3_2, d_n3_1, d_n3_2] = d_param

B3 = quad(g2_fit, f3[0], f3[-1], args=tuple(param3))
d_B3 = 0.02 * B3

# Determination of boltzmann's constant
k = s / (4 * T1 * B3)
d_k_stat = k * sqrt((d_s / s)**2 + (d_T1 / T1)**2)
d_k_sys = k * d_B3 / B3
d_k = sqrt(d_k_stat**2 + d_k_sys**2)

k_theo = cs.physical_constants["Boltzmann constant"]
d_k_theo = k_theo[2]
k_theo = k_theo[0]

print(ms.val('k', k, d_k_stat, d_k_sys, unit='J/K', prefix=False))
print(ms.val('k', k_theo, k_theo, unit='J/K', prefix=False))
print(ms.dev(k, d_k, k_theo, d_k_theo, name='k'))
