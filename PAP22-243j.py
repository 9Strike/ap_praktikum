import measure as ms
import numpy as np
import scipy.constants as cs

from measure import npfarray as fa
from measure import sqrt, ln
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
n1 = fa([103, 196, 104, 104, 105, 106])
U1 = fa([2.4233, 3.1416, 3.7295, 4.238, 4.6961, 5.1225]) * cs.milli
d_U1 = fa([0.01, 0.0113, 0.0149, 0.0165, 0.0173, 0.0171]) * cs.milli / sqrt(n1)
d_U1 = sqrt(d_U1**2 + (rd_V * U1)**2)
T1 = 23.4 + cs.zero_Celsius
d_T1 = 0.1

# (2) Measurement of the inherent noise voltage of the multiplier
n2 = 104
U2 = 1.3698 * cs.milli
d_U2 = 0.00544 * cs.milli / sqrt(n2)
d_U2 = sqrt(d_U1**2 + (rd_V * U1)**2)

U_diff = U1**2 - U2**2
d_U_diff = sqrt((2 * U1 * d_U1)**2 + (2 * U2 * d_U2)**2)

ms.pltext.initplot(num=1, title=titles[0], xlabel=r'$R$ / k$\Omega$', ylabel=r'$(U_N**2 - U_V**2)$ / mV$**2$', scale='loglog', fignum=True)
[s, d_s, b, d_b] = ms.linreg(R1, U_diff, d_U_diff, plot=True)

# Measurement of the measurement system's frequency slope
D = 0.001
d_D = 0.002 * D

def g_fit(f, V, Omg1, Omg2, n1, n2):
  return V / sqrt((1 + (Omg1 / f)**(2 * n1)) * (1 + (f / Omg2)**(2 * n2)))
def g2_fit(f, V, Omg1, Omg2, n1, n2):
  return g_fit(f, V, Omg1, Omg2, n1, n2)**2

def g2_fit_partial1(f, V, Omg1, Omg2, n1, n2):
  return (2*V)/(((Omg1/f)**(2*n1)+1)*(f/Omg2)**(2*n2)+(Omg1/f)**(2*n1)+1)
def g2_fit_partial2(f, V, Omg1, Omg2, n1, n2):
  return -(2*V**2*(Omg1/f)**(2*n1)*n1)/((Omg1*(Omg1/f)**(4*n1)+2*Omg1*(Omg1/f)**(2*n1)+Omg1)*(f/Omg2)**(2*n2)+Omg1*(Omg1/f)**(4*n1)+2*Omg1*(Omg1/f)**(2*n1)+Omg1)
def g2_fit_partial3(f, V, Omg1, Omg2, n1, n2):
  return (2*V**2*(f/Omg2)**(2*n2)*n2)/((Omg2*(Omg1/f)**(2*n1)+Omg2)*(f/Omg2)**(4*n2)+(2*Omg2*(Omg1/f)**(2*n1)+2*Omg2)*(f/Omg2)**(2*n2)+Omg2*(Omg1/f)**(2*n1)+Omg2)
def g2_fit_partial4(f, V, Omg1, Omg2, n1, n2):
  return -(2*V**2*ln(Omg1/f)*(Omg1/f)**(2*n1))/(((Omg1/f)**(4*n1)+2*(Omg1/f)**(2*n1)+1)*(f/Omg2)**(2*n2)+(Omg1/f)**(4*n1)+2*(Omg1/f)**(2*n1)+1)
def g2_fit_partial5(f, V, Omg1, Omg2, n1, n2):
  return -(2*V**2*(f/Omg2)**(2*n2)*ln(f/Omg2))/(((Omg1/f)**(2*n1)+1)*(f/Omg2)**(4*n2)+(2*(Omg1/f)**(2*n1)+2)*(f/Omg2)**(2*n2)+(Omg1/f)**(2*n1)+1)

# (3) Measurement of the measurement system's frequency slope and the equivalent noise bandwidth
U3_E = 0.2
f3, U3_A = np.loadtxt('data/243/data1.txt', skiprows=1, usecols=(0, 1), unpack=True)
g3 = U3_A / (D * U3_E)
d_g3 = g3 * d_D / D

ms.pltext.initplot(num=2, title=titles[1], xlabel=r'$f$ / Hz', ylabel=r'$g(f)$', scale='loglog', fignum=True)
param3, d_param3 = ms.fit(f3, g3, d_g3, g_fit, p0=[1000, 1000, 50000, 5, 5], fit_range=range(17, 138), plot=True)
[V3, Omg3_1, Omg3_2, n3_1, n3_2] = param3
[d_V3, d_Omg3_1, d_Omg3_2, d_n3_1, d_n3_2] = d_param3

B3, d_B3 = quad(g2_fit, f3[0], f3[-1], args=tuple(param3))
#d_B3 = 0.02 * B3
d_B3 = d_B3**2
d_B3 += quad(g2_fit_partial1, f3[0], f3[-1], args=tuple(param3))[0]**2 * d_param3[0]**2
d_B3 += quad(g2_fit_partial2, f3[0], f3[-1], args=tuple(param3))[0]**2 * d_param3[1]**2
d_B3 += quad(g2_fit_partial3, f3[0], f3[-1], args=tuple(param3))[0]**2 * d_param3[2]**2
d_B3 += quad(g2_fit_partial4, f3[0], f3[-1], args=tuple(param3))[0]**2 * d_param3[3]**2
d_B3 += quad(g2_fit_partial5, f3[0], f3[-1], args=tuple(param3))[0]**2 * d_param3[4]**2
d_B3 = sqrt(d_B3)

# (4)
R4_1 = fa([4792.6, 5559.6, 6309.4, 7062.9, 7704])
R4_2 = fa([4794.6, 5562.9, 6312.9, 7065.6, 7725.4])
T4_1 = fa([51.084, 101.27, 151.1, 201.95, 245.86])
T4_2 = fa([51.214, 101.49, 151.33, 202.13, 247.34])
n4 = fa([108, 106, 104, 113, 105])
U4 = fa([3.1927, 2.9036, 3.0386, 3.6315, 3.7054]) * cs.milli
d_U4 = fa([0.469, 0.0162, 0.0146, 0.0283, 0.0226]) * cs.milli / sqrt(n4)

# Determination of boltzmann's constant
k = s / (4 * T1 * B3)
d_k_stat = k * sqrt((d_s / s)**2 + (d_T1 / T1)**2)
d_k_sys = k * d_B3 / B3
d_k = sqrt(d_k_stat**2 + d_k_sys**2)

k_theo = cs.physical_constants["Boltzmann constant"]
d_k_theo = k_theo[2]
k_theo = k_theo[0]

print(ms.val('k', k, d_k_stat, d_k_sys, unit='J/K', prefix=False))
print(ms.val('k', k_theo, d_k_theo, unit='J/K', prefix=False))
print(ms.dev(k, d_k, k_theo, d_k_theo, name='k'))

# Show plots
ms.plt.show()
