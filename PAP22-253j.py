import matplotlib.pyplot as plt
from measure import npfarray as npf
from measure import sqrt
import measure as ms
import numpy as np
import scipy.constants as cs

# Measured data
r_c = 14 * cs.milli / 2
U0_c = 520

# Measurement of the zero-effect
t0 = 5 * cs.minute
n0 = 122 / t0

print(ms.val("n0", n0 * 2 * cs.minute))

# Measurement of β-Radiation absorption, 90Sr, GS 527
A_Sr = 74 * cs.kilo
s_Sr = 60 * cs.milli
d_s_Sr = 2 * cs.milli
t1_Sr = 30
t2_Sr = 2 * cs.minute
t3_Sr = 5 * cs.minute
n_Sr = npf([1136, 665, 412, 273, 180, 110, 306, 208, 131, 76, 60, 42])
d_n_Sr = sqrt(n_Sr)
n_Sr = np.append(n_Sr[:6] / t1_Sr, n_Sr[6:] / t2_Sr)
d_n_Sr = np.append(d_n_Sr[:6] / t1_Sr, d_n_Sr[6:] / t2_Sr)
d_Sr = np.arange(0.0, 0.3 * cs.milli * len(n_Sr), 0.3 * cs.milli)
n0_β = 51 / t3_Sr

ms.pltext.initplot(num=1, xlabel='d / mm', ylabel='n / (1/s)')
ms.pltext.plotdata(d_Sr / cs.milli, n_Sr, d_n_Sr)

# Measurement of γ-Radiation absorption, 60Co, UB 595
A_Co = 3.7 * cs.mega
s_Co = 150 * cs.milli
d_s_Co = 2 * cs.milli
t_Co = cs.minute
n_Co = npf([2447, 1760, 1279, 908, 714, 541, 412, 312, 223, 195, 164])
d_n_Co = sqrt(n_Co) / t_Co
n_Co = n_Co / t_Co
d_Co = np.arange(0.0, 5 * cs.milli * len(n_Co), 5 * cs.milli)

ms.pltext.initplot(num=2, xlabel='d / mm', ylabel='n / (1/s)')
ms.pltext.plotdata(d_Co / cs.milli, n_Co, d_n_Co)

# Measurement of γ-Radiation activity, 60Co, UB 595
t_A = cs.minute
s_A = npf([50, 105, 190]) * cs.milli
d_s_A = npf([2, 2, 2]) * cs.milli
n_A = npf([33865, 8266, 2171])
d_n_A = sqrt(n_A) / t_A
n_A = n_A / t_A

# Measurement of α-Radiation absorption and energy, 241Am, AP 15.2
s_c = 4.2 * cs.centi
σ_c = 2.25 * cs.milli / cs.centi**2
A_Am = 90 * cs.kilo
t_Am = cs.minute
p1_Am = npf([18, 98, 120, 225, 324, 383, 416, 450, 475, 515, 617, 721, 813, 911, 1013]) * cs.milli * cs.bar
p2_Am = npf([20, 98, 120, 222, 324, 392, 413, 450, 473, 512, 614, 717, 809, 908, 1009]) * cs.milli * cs.bar
d_p_Am = cs.milli * cs.bar
n_Am = npf([13144, 13142, 13131, 12933, 12615, 9883, 7101, 3491, 1680, 451, 239, 205, 229, 225, 212])
d_n_Am = sqrt(n_Am) / t_Am
n_Am = n_Am / t_Am

ms.pltext.initplot(num=3, xlabel='p / Pa', ylabel='n / (1/s)')
ms.pltext.plotdata(p1_Am, n_Am, d_n_Am)

ms.plt.show()
