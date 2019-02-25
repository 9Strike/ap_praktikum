import measure as ms
from measure import sqrt
from measure import npfarray as npf
import numpy as np
import scipy.constants as cs

ms.plt.rc('text', usetex=True)
ms.plt.rc('font', family='serif')

titles = [
  r'Zählrate $n$ in Abhängigkeit der Spannung der Röntgenröhre $U$ bei einem konstanten Winkel von 7.5$^\circ$'
]

# Counting rate - Voltage dependency Measurement
t = 20
φ = 7.5 * cs.degree
U = np.arange(20.0, 36.0, 1.0) * cs.kilo
n = npf([1.35, 1.35, 2.75, 5.55, 32.95, 78.35, 122.8, 163.3, 200.6, 237.0, 270.2, 307.6, 337.1, 374.7, 403.7, 433.3])
d_n = sqrt(n * t) / t

ms.pltext.initplot(num=1, title=titles[0], xlabel=r'$U$ / V', ylabel=r'$n$ / (1/s)')
s, d_s, i, d_i = ms.linreg(U, n, d_n, fit_range=range(3, len(U)), plot=True)

ms.plt.show()
