import measure as ms
import numpy as np
import scipy.constants as cs

from measure import npfarray as npfarray
from measure import sqrt

ms.plt.rc('text', usetex=True)
ms.plt.rc('font', family='serif')

titles = [
  r'Ausgangsspannung $U_A$ in Abhängigkeit der Eingangsgleichspannung $U_E$ für verschiedene Gegenkopplungswiderstände',
  r'Ausgangsspannung $U_A$ in Abhängigkeit der Eingangswechselspannung $U_E$ für verschiedene Gegenkopplungswiderstände',
  r'Verstärkung $V$ in Abhängigkeit der Wechselspannungsfrequenz $f$ für verschiedene Schaltungskonfigurationen'
]

# (1) Dependency of input and output voltage
# DC
R1_G = 48.7 * cs.kilo
R1_E = 3 * cs.kilo
U1_E = npfarray([-250, -190, -120, -60.3, -1.3, 60.7, 120, 250]) * cs.milli / 2
d_U1_E = npfarray([1, 1, 1, 0.2, 0.2, 0.1, 1, 1]) * cs.milli / 2
U1_A = npfarray([4.24, 3.19, 2.11, 1.10, 0.16, -0.85, -1.80, -3.96]) / 2
d_U1_A = npfarray([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]) / 2

R2_G = 274 * cs.kilo
R2_E = 3 * cs.kilo
U2_E = npfarray([-250, -190, -120, -60.1, -1.5, 80.0, 160, 250]) * cs.milli / 2
d_U2_E = npfarray([1, 1, 1, 0.1, 0.1, 0.1, 1, 1]) * cs.milli / 2
U2_A = npfarray([14.4, 14.4, 10.9, 5.28, 0.82, -6.80, -13.0, -13.0]) / 2
d_U2_A = npfarray([0.1, 0.1, 0.1, 0.01, 0.01, 0.01, 0.1, 0.1]) / 2

# AC
R3_G = 274 * cs.kilo
R3_E = 3 * cs.kilo
U3_G = npfarray([150, 300, 450, 600, 750, 900, 950, 1000]) * cs.milli / 2
U3_E = npfarray([128, 250, 372, 500, 620, 744, 792, 832]) * cs.milli / 20.0
d_U3_E = npfarray([1, 1, 1, 1, 1, 1, 1, 1]) * cs.milli / 20.0
U3_A = npfarray([1.16, 2.30, 3.40, 4.52, 5.68, 6.80, 7.12, 7.52]) / 2
d_U3_A = npfarray([0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) / 2

R4_G = 680 * cs.kilo
R4_E = 3 * cs.kilo
U4_G = npfarray([150, 300, 450, 600, 750, 900, 950, 1000]) * cs.milli / 2
U4_E = npfarray([136, 250, 372, 500, 620, 744, 792, 832]) * cs.milli / 20.0
d_U4_E = npfarray([1, 1, 1, 1, 1, 1, 1, 1]) * cs.milli / 20.0
U4_A = npfarray([3.00, 5.48, 8.24, 10.90, 13.60, 16.60, 17.40, 18.20]) / 2
d_U4_A = npfarray([0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) / 2

ms.pltext.initplot(num=1, title=titles[0], xlabel=r'$U_E$ / V', ylabel=r'$U_A$ / V', fignum=True)
s1, d_s1, b1, d_b1 = ms.linreg(U1_E, U1_A, d_U1_A, d_U1_E, plot=True, graphname=r'$R_G$ = 48.7 k$\Omega$')
s2, d_s2, b2, d_b2 = ms.linreg(U2_E, U2_A, d_U2_A, d_U2_E, fit_range=range(2, 6), plot=True, graphname=r'$R_G$ = 274 k$\Omega$')

ms.pltext.initplot(num=2, title=titles[1], xlabel=r'$U_E$ / V', ylabel=r'$U_A$ / V', fignum=True)
s3, d_s3, b3, d_b3 = ms.linreg(U3_E, U3_A, d_U3_A, d_U3_E, plot=True, graphname=r'$R_G$ = 274 k$\Omega$')
s4, d_s4, b4, d_b4 = ms.linreg(U4_E, U4_A, d_U4_A, d_U4_E, plot=True, graphname=r'$R_G$ = 680 k$\Omega$')
1
V1 = -s1
d_V1 = d_s1
V2 = -s2
d_V2 = d_s2
V3 = s3
d_V3 = d_s3
V4 = s4
d_V4 = d_s4

V1_T = R1_G / R1_E
V2_T = R2_G / R2_E
V3_T = R3_G / R3_E
V4_T = R4_G / R4_E

print()
print(ms.tbl([ms.lst(npfarray([R1_E, R2_E]), name='R_E', unit='Ω'),
              ms.lst(npfarray([R1_G, R2_G]), name='R_G', unit='Ω'),
              ms.lst(npfarray([V1, V2]), npfarray([d_V1, d_V2]), name='V'),
              ms.lst(npfarray([V1_T, V2_T]), name='V'),
              ms.dev(npfarray([V1, V2]), npfarray([d_V1, d_V2]) ,npfarray([V1_T, V2_T]), name='V', perc=True)]))
print()
print(ms.tbl([ms.lst(npfarray([R3_E, R4_E]), name='R_E', unit='Ω'),
              ms.lst(npfarray([R3_G, R4_G]), name='R_G', unit='Ω'),
              ms.lst(npfarray([V3, V4]), npfarray([d_V3, d_V4]), name='V'),
              ms.lst(npfarray([V3_T, V4_T]), name='V'),
              ms.dev(npfarray([V3, V4]), npfarray([d_V3, d_V4]) ,npfarray([V3_T, V4_T]), name='V', perc=True)]))
print()

# (2) Frequency slope of the amplification
R5_G = 680 * cs.kilo
U5_G = 0.3 / 2
f5 = npfarray([0.1, 0.3, 0.6, 1.0, 3.0, 6.0, 10.0, 30.0, 60.0, 100.0, 300.0]) * cs.kilo
U5_E = npfarray([176, 250, 248, 250, 250, 250, 250, 252, 250, 252, 250]) * cs.milli / 20.0
d_U5_E = npfarray([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]) * cs.milli / 20.0
U5_A = npfarray([4.12, 5.76, 5.72, 5.52, 4.12, 2.60, 1.68, 0.572, 0.296, 0.184, 0.060]) / 2
d_U5_A = npfarray([0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.004, 0.004, 0.004, 0.004]) / 2

R6_G = 274 * cs.kilo
U6_G = 0.3 / 2
f6 = npfarray([0.1, 0.3, 0.6, 1.0, 3.0, 6.0, 10.0, 30.0, 60.0, 100.0, 300.0]) * cs.kilo
U6_E = npfarray([250, 248, 250, 250, 250, 250, 248, 252, 254, 250, 250]) * cs.milli / 20.0
d_U6_E = npfarray([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]) * cs.milli / 20.0
U6_A = npfarray([2.30, 2.32, 2.28, 2.30, 2.12, 1.78, 1.38, 0.560, 0.296, 0.174, 0.062]) / 2
d_U6_A = npfarray([0.04, 0.04, 0.04, 0.02, 0.02, 0.02, 0.02, 0.004, 0.004, 0.004, 0.004]) / 2

R7_G = 48.7 * cs.kilo
U7_G = 0.3 / 2
f7 = npfarray([0.1, 0.3, 0.6, 1.0, 3.0, 6.0, 10.0, 30.0, 60.0, 100.0, 300.0]) * cs.kilo
U7_E = npfarray([832, 832, 832, 832, 840, 848, 840, 848, 848, 840, 848]) * cs.milli / 20.0
d_U7_E = npfarray([2, 2, 2, 2, 2, 2, 8, 8, 8, 8, 8]) * cs.milli / 20.0
U7_A = npfarray([7.60, 7.60, 7.52, 7.52, 7.12, 6.0, 4.56, 1.84, 0.98, 0.56, 0.19]) / 2
d_U7_A = npfarray([0.02, 0.02, 0.04, 0.04, 0.04, 0.8, 0.04, 0.02, 0.02, 0.02, 0.02]) / 2

R8_G = 48.7 * cs.kilo
C8_G = 560 * cs.pico
U8_G = 1.0 / 2
f8 = npfarray([0.1, 0.3, 0.6, 1.0, 3.0, 6.0, 10.0, 30.0, 60.0, 100.0, 300.0]) * cs.kilo
U8_E = npfarray([832, 832, 832, 840, 848, 840, 840, 848, 848, 840, 840]) * cs.milli / 20.0
d_U8_E = npfarray([2, 2, 2, 2, 8, 8, 8, 8, 8, 8, 8]) * cs.milli / 20.0
U8_A = npfarray([1.36, 1.36, 1.36, 1.36, 1.26, 0.98, 0.68, 0.27, 0.140, 0.83, 0.33]) / 2
d_U8_A = npfarray([0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.004, 0.02, 0.02]) / 2

R9_G = 48.7 * cs.kilo
C9_E = 47 * cs.nano
U9_G = 1.0 / 2
f9 = npfarray([0.3, 0.6, 1.0, 3.0, 6.0, 10.0, 20.0]) * cs.kilo
U9_E = npfarray([832, 832, 832, 848, 848, 840, 848]) * cs.milli / 20.0
d_U9_E = npfarray([8, 8, 8, 8, 8, 8, 8]) * cs.milli / 20.0
U9_A = npfarray([0.388, 0.680, 0.936, 1.34, 1.38, 1.36, 1.26]) / 2
d_U9_A = npfarray([0.004, 0.004, 0.004, 0.02, 0.02, 0.02, 0.02]) / 2

V5 = U5_A / U5_E
d_V5 = V5 * sqrt((d_U5_A / U5_A)**2 + (d_U5_E / U5_E)**2)
V6 = U6_A / U6_E
d_V6 = V6 * sqrt((d_U6_A / U6_A)**2 + (d_U6_E / U6_E)**2)
V7 = U7_A / U7_E
d_V7 = V7 * sqrt((d_U7_A / U7_A)**2 + (d_U7_E / U7_E)**2)
V8 = U8_A / U8_E
d_V8 = V8 * sqrt((d_U8_A / U8_A)**2 + (d_U8_E / U8_E)**2)
V9 = U9_A / U9_E
d_V9 = V9 * sqrt((d_U9_A / U9_A)**2 + (d_U9_E / U9_E)**2)

ms.pltext.initplot(num=3, title=titles[2], xlabel=r'$f$ / Hz', ylabel=r'$V$', scale='loglog', fignum=True)
ms.pltext.plotdata(f5, V5, d_V5, label= r'$R_G$ = 680 k$\Omega$', connect=True)
ms.pltext.plotdata(f6, V6, d_V6, label= r'$R_G$ = 274 k$\Omega$', connect=True)
ms.pltext.plotdata(f7, V7, d_V7, label= r'$R_G$ = 48.7 k$\Omega$', connect=True)
ms.pltext.plotdata(f8, V8, d_V8, label= r'$R_G$ = 48.7 k$\Omega$, $C_G$ = 560 pF', connect=True)
ms.pltext.plotdata(f9, V9, d_V9, label= r'$R_G$ = 48.7 k$\Omega$, $C_E$ = 47 nF', connect=True)

# Show plots
#ms.pltext.savefigs('figures/242')
#ms.plt.show()
