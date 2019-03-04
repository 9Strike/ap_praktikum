import measure as ms
from measure import npfarray as npf
from measure import pi
import numpy as np
import scipy.constants as cs

f = np.arange(1, 11, 1) * cs.kilo
t = npf([0.20, 0.08, 0.042, 0.027, 0.019, 0.013, 0.010, 0.007, 0.007, 0.005]) * cs.milli
d_t = npf([0.03, 0.02, 0.015, 0.015, 0.010, 0.010, 0.005, 0.005, 0.005, 0.004]) * cs.milli
phi = 2 * pi * f * t
d_phi = 2 * pi * f * d_t

print(ms.tbl([ms.lst(phi / cs.degree, d_phi / cs.degree, name='φ', unit='°', prefix=False)]))
