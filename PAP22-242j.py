import measure as ms
import numpy as np
import scipy.constants as cs

from measure import npfarray as fa
from measure import sqrt

ms.plt.rc('text', usetex=True)
ms.plt.rc('font', family='serif')

# (2) Frequency slope
U2_G = 300 * cs.milli
f2 = fa([0.1, 0.3, 0.6, 1.0, 3.0, 6.0, 10.0, 30.0, 60.0, 100.0, 300.0]) * cs.kilo
U2_E = fa([176, 250, 248, 250, 250, 250, 250, 252, 250, 252, 250]) * cs.milli / 10.0
d_U2_E = fa([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]) * cs.milli / 10.0
U2_A = fa([4.12, 5.76, 5.72, 5.52, 4.12, 2.60, 1.68, 0.572, 0.296, 0.184, 0.060])
d_U2_A = fa([0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.004, 0.004, 0.004, 0.004])

f3 = fa([0.1, 0.3, 0.6, 1.0, 3.0, 6.0, 10.0, 30.0, 60.0, 100.0, 300.0]) * cs.kilo
U3_E = fa([250, 248, 250, 250, 250, 250, 248, 252, 254, 250, 250]) * cs.milli / 10.0
d_U3_E = fa([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]) * cs.milli / 10.0
U3_A = fa([2.30, 2.32, 2.28, 2.30, 2.12, 1.78, 1.38, 0.560, 0.296, 0.174, 0.062])
d_U3_A = fa([0.04, 0.04, 0.04, 0.02, 0.02, 0.02, 0.02, 0.004, 0.004, 0.004, 0.004])

f4 = fa([0.1, 0.3, 0.6, 1.0, 3.0, 6.0, 10.0, 30.0, 60.0, 100.0, 300.0]) * cs.kilo
U4_E = fa([832, 832, 832, 832, 840, 848, 840, 848, 848, 840, 848]) * cs.milli / 10.0
d_U4_E = fa([2, 2, 2, 2, 2, 2, 8, 8, 8, 8, 8]) * cs.milli / 10.0
U4_A = fa([7.60, 7.60, 7.52, 7.52, 7.12, 6.0, 4.56, 1.84, 0.98, 0.56, 0.19])
d_U4_A = fa([0.02, 0.02, 0.04, 0.04, 0.04, 0.8, 0.04, 0.02, 0.02, 0.02, 0.02])

V2 = U2_A / U2_E
d_V2 = V2 * sqrt((d_U2_A / U2_A)**2 + (d_U2_E / U2_E)**2)
V3 = U3_A / U3_E
d_V3 = V3 * sqrt((d_U3_A / U3_A)**2 + (d_U3_E / U3_E)**2)
V4 = U4_A / U4_E
d_V4 = V4 * sqrt((d_U4_A / U4_A)**2 + (d_U4_E / U4_E)**2)

ms.pltext.initplot(num=1, xlabel=r'$f$ / Hz', ylabel=r'$V$', scale='loglog')
ms.pltext.plotdata(f2, V2, d_V2, connect=True)
ms.pltext.plotdata(f3, V3, d_V3, connect=True)
ms.pltext.plotdata(f4, V4, d_V4, connect=True)

# Show plots
ms.plt.show()