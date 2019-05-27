import datplot as dp
import datstat as dst
import datstr as ds
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as cs

from numpy import sqrt



### General

print()

titles = [
  r'',
  r'',
  r'',
  r'',
  r'',
  r'',
  r''
]

dataConverters = {
  0 : lambda s: s.decode('utf-8').replace(',', '.'),
  1 : lambda s: s.decode('utf-8').replace(',', '.')
}

ER = cs.Rydberg * cs.h * cs.c

lda_fl_t = np.array([
  759.4, 686.7, 656.3, 589.6, 589.0, 587.6, 527.0, 518.4,
  486.1, 430.8, 396.8, 393.4
]) * cs.nano
lda_fl_t_ = np.array([
  759.4, 686.7, 656.3, 589.3, 527.0, 518.4, 486.1, 430.8,
  396.8, 393.4
]) * cs.nano



### Measured data

## Sun light
lda_g, I_g = np.loadtxt('data/234/sky_glass.txt', skiprows=17,
                        converters=dataConverters, comments='>', unpack=True)
lda_sk, I_sk = np.loadtxt('data/234/sky.txt', skiprows=17,
                        converters=dataConverters, comments='>', unpack=True)
lda_s, I_s = np.loadtxt('data/234/sun.txt', skiprows=17,
                        converters=dataConverters, comments='>', unpack=True)

lda1_fl_s = np.array([
  759.74, 686.95, 656.03, 589.44, 527.27, 517.44, 486.93, 431.61,
  398.27, 394.87
]) * cs.nano
lda2_fl_s = np.array([
  760.91, 687.70, 656.96, 590.45, 528.13, 518.25, 487.62, 432.60,
  398.99, 395.83
]) * cs.nano

## Natrium

lda_nw, I_nw = np.loadtxt('data/234/natrium_weak.txt', skiprows=17,
                          converters=dataConverters, comments='>', unpack=True)
lda_ns, I_ns = np.loadtxt('data/234/natrium_strong.txt', skiprows=17,
                          converters=dataConverters, comments='>', unpack=True)

lda1_ol_na = np.array([
  332.43, 498.67, 568.46, 588.86, 615.67, 695.97, 706.12, 726.66,
  737.61, 749.57, 762.57, 765.56, 769.04, 771.46, 793.81, 799.81,
  810.30, 817.11, 825.34, 839.73, 841.18, 306.88, 332.18, 396.27,
  405.89, 417.00, 420.83, 427.60, 430.99, 434.43, 451.35, 456.27,
  467.49, 475.88, 493.96, 498.47, 515.54, 641.47, 649.41, 667.42,
  674.92, 686.70, 695.82, 702.57, 705.96, 714.08, 720.05, 726.49,
  737.46, 749.38, 762.39, 765.54, 768.93, 771.20, 793.66, 799.78,
  809.76, 817.36, 825.22, 839.47, 840.79
]) * cs.nano
lda2_ol_na = np.array([
  333.89, 499.97, 570.34, 591.08, 617.00, 697.12, 707.23, 727.62,
  738.72, 751.74, 763.82, 766.75, 770.07, 772.65, 794.90, 801.62,
  811.65, 819.94, 826.45, 840.75, 842.45, 308.32, 334.55, 397.47,
  406.89, 418.70, 422.61, 429.31, 432.16, 436.04, 453.84, 457.50,
  468.85, 477.14, 495.22, 500.47, 517.20, 642.50, 650.45, 668.48,
  676.08, 687.99, 697.65, 703.68, 707.81, 715.47, 721.41, 728.14,
  739.40, 752.75, 764.73, 766.91, 770.64, 773.46, 795.52, 802.09,
  812.45, 820.60, 826.91, 841.26, 843.28
]) * cs.nano



### Data Preparation

## Sun light
lda_g *= cs.nano
lda_sk *= cs.nano
lda_s *= cs.nano

lda_fl_s = 0.5 * (lda1_fl_s + lda2_fl_s)
d_lda_fl_s = 0.5 * (lda2_fl_s - lda1_fl_s)


## Natrium

lda_ol_na = 0.5 * (lda1_ol_na + lda2_ol_na)
d_lda_ol_na = 0.5 * (lda2_ol_na - lda1_ol_na)

i_lda_ol_sorted_na = np.argsort(lda_ol_na)
lda_ol_na = lda_ol_na[i_lda_ol_sorted_na]
d_lda_ol_na = d_lda_ol_na[i_lda_ol_sorted_na]



### Evaluation

## Sun light
dp.initplot(num=1, title=titles[0], xlabel=r'$\lambda$ / nm', ylabel=r'$I$ / b.E.')
plt.xlim(250, 900)
dp.plot(lda_g / cs.nano, I_g, label='Spektrum des Himmels durch Glas')
dp.plot(lda_sk / cs.nano, I_sk, label='Spektrum des Himmels')

A_g = 1 - I_g / I_sk
dp.initplot(num=2, title=titles[1], xlabel=r'$\lambda$ / nm', ylabel=r'$A$ / b.E.')
plt.xlim(320, 800)
plt.ylim(-0.1, 1.0)
dp.plot(lda_g / cs.nano, A_g)

dp.initplot(num=3, title=titles[2], xlabel=r'$\lambda$ / nm', ylabel=r'$I$ / b.E.')
plt.xlim(350, 800)
dp.plot(lda_s / cs.nano, I_s)

# Print
print(ds.tbl([
  ds.lst(lda_fl_s, d_lda_fl_s, name='λ', unit='m'),
  ds.lst(lda_fl_t_, name='λ_lit', unit='m'),
  ds.dev(lda_fl_s, d_lda_fl_s, lda_fl_t_, name='λ, λ_lit', perc=True)
]))
print()


## Natrium
dp.initplot(num=4, title=titles[3], xlabel=r'$\lambda$ / nm', ylabel=r'$I$ / b.E.', scale='linlog')
plt.xlim(300, 850)
dp.plot(lda_ns, I_ns)

dp.initplot(num=5, title=titles[4], xlabel=r'$\lambda$ / nm', ylabel=r'$I$ / b.E.', scale='linlog')
plt.xlim(300, 540)
dp.plot(lda_nw, I_nw)

dp.initplot(num=5, title=titles[4], xlabel=r'$\lambda$ / nm', ylabel=r'$I$ / b.E.', scale='linlog')
plt.xlim(600, 850)
dp.plot(lda_nw, I_nw)

# 1. Side series
m1_ref_na = 3
lda1_ref_na = 819 * cs.nano
m1_na = np.arange(3, 13)

i1_ref_na = np.argmin(np.abs(lda_ol_na - lda1_ref_na))
E_3p_na = -ER / m1_ref_na**2 - cs.h * cs.c / lda_ol_na[i1_ref_na]
lda1_tl_na = cs.c * cs.h / (-ER / m1_na**2 - E_3p_na)

# 2. Side series
lda2_ref_na = 589 * cs.nano
m2_na = np.arange(4, 10)

i2_ref_na = np.argmin(np.abs(lda_ol_na - lda2_ref_na))
E_3s_na = E_3p_na - cs.h * cs.c / lda_ol_na[i2_ref_na]
delta_s_na = 3.0 - sqrt(-ER / E_3s_na)
lda2_tl_na = cs.c * cs.h / (-ER / (m2_na - delta_s_na)**2 - E_3p_na)

# Main series
m3_na = np.arange(4, 6)

delta_p_na = 3.0 - sqrt(-ER / E_3p_na)
lda3_tl_na = cs.c * cs.h / (-ER / (m3_na - delta_p_na)**2 - E_3s_na)

# Arrange theoretical wavelengths
m_na = np.concatenate((m1_na, m2_na, m3_na))
lda_tl_na = np.concatenate((lda1_tl_na, lda2_tl_na, lda3_tl_na))

i_lda_tl_sorted_na = np.argsort(lda_tl_na)
m_na = m_na[i_lda_tl_sorted_na]
lda_tl_na = lda_tl_na[i_lda_tl_sorted_na]

# Print
print(ds.tbl([
  ds.lst(lda_ol_na, d_lda_ol_na, name='λ_o', unit='m')
]))
print(ds.tbl([
  ds.lst(lda_tl_na, name='λ_t', unit='m')
]))

i_assign_na = np.array([np.argmin(np.abs(lda_ol_na - lda)) for lda in lda_tl_na])
d_lda_ol_na = d_lda_ol_na[i_assign_na]
lda_ol_na = lda_ol_na[i_assign_na]

print(ds.tbl([
  ds.lst(lda_ol_na, d_lda_ol_na, name='λ_o', unit='m'),
  ds.lst(lda_tl_na, name='λ_t', unit='m'),
  ds.dev(lda_ol_na, d_lda_ol_na, lda_tl_na, name='λ_o, λ_t', perc=True)
]))

isInAccScope_na = dst.dev(lda_ol_na, d_lda_ol_na, lda_tl_na) <= 3.0
m_na = m_na[isInAccScope_na]
lda_ol_na = lda_ol_na[isInAccScope_na]
d_lda_ol_na = d_lda_ol_na[isInAccScope_na]
lda_tl_na = lda_tl_na[isInAccScope_na]

print(ds.tbl([
  ds.lst(lda_ol_na, d_lda_ol_na, name='λ_o', unit='m'),
  ds.lst(lda_tl_na, name='λ_t', unit='m'),
  ds.dev(lda_ol_na, d_lda_ol_na, lda_tl_na, name='λ_o, λ_t', perc=True)
]))

# Plot of the first side series
def lda1_func(m, E_Ry, E_3p, delta_d):
  return cs.c * cs.h / (E_Ry / (m - delta_d)**2 - E_3p)

isInLda1_na = np.in1d(lda_tl_na, lda1_tl_na)
m1_na = m_na[isInLda1_na]
lda1_ol_na = lda_ol_na[isInLda1_na]
d_lda1_ol_na = d_lda_ol_na[isInLda1_na]

dp.initplot(num=6, title=titles[5], xlabel=r'$m$', ylabel=r'$\lambda$ / nm')
popt1_na, d_popt1_na = dp.fit(m1_na, lda1_ol_na, d_lda1_ol_na, lda1_func, p0=[-13.6 * cs.e, -3 * cs.e, -0.02], plot=True)

chi2_1_na = dst.chi2(lda1_ol_na, d_lda1_ol_na, lda1_func(m1_na, *popt1_na))
chi2_1_red_na = chi2_1_na / (len(lda1_ol_na) - 3)

# Plot of the second side series
def lda2_func(m, E_Ry, E_3p, delta_s):
  return cs.c * cs.h / (E_Ry / (m - delta_s)**2 - E_3p)

isInLda2_na = np.in1d(lda_tl_na, lda2_tl_na)
m2_na = m_na[isInLda2_na]
lda2_ol_na = lda_ol_na[isInLda2_na]
d_lda2_ol_na = d_lda_ol_na[isInLda2_na]

dp.initplot(num=7, title=titles[6], xlabel=r'$m$', ylabel=r'$\lambda$ / nm')
popt2_na, d_popt2_na = dp.fit(m2_na, lda2_ol_na, d_lda2_ol_na, lda2_func, p0=[-13.6 * cs.e, -3 * cs.e, -0.02], plot=True)

chi2_2_na = dst.chi2(lda2_ol_na, d_lda2_ol_na, lda2_func(m2_na, *popt2_na))
chi2_2_red_na = chi2_2_na / (len(lda2_ol_na) - 3)

# Print
print(ds.val('E_Ry', popt1_na[0], d_popt1_na[0]))
print(ds.val('E_3p', popt2_na[1], d_popt2_na[1]))
print(ds.val('E_delta_d'))



### Show plots
plt.show()
