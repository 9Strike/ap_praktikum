### measure libraby version 1.8.1
import numpy as np
import matplotlib.pyplot as plt

# Settings
linreg_change = 0.00001 # min relative change per step to end linear regression
minfloat = 1e-80 # replaces zeros in linreg
default_dpi = 144.0 # default dpi for plots

# Variables for export
sqrt = np.sqrt
exp = np.exp
ln = np.log
log10 = np.log10
sin = np.sin
cos = np.cos
tan = np.tan
arctan = np.arctan
pi = np.pi
euler_e = np.e
c = 2.99792458e8 # Speed of light
h = 6.62607015e-34 # Planck's constant
e = 1.602176634e-19 # Elementary charge
kB = 1.380649e-23 # Boltzmann constant
NA = 6.02214076e23 # Avogadro constant
T0 = 273.15 # Zero Celsius in Kelvin
g = 9.80984 # Gravitantional acceleration in Heidelberg 
dg = 2e-5 # Uncertanty of the gravitational acceleration

def npfarray(x):
  return np.array(x, dtype='float')

def mv(x):
  s = 0.0
  for i in range(len(x)):
    s += x[i]
  return s / len(x)

def dsto(x):
  s = 0.0
  for i in range(len(x)):
    s += (x[i] - mv(x))**2
  return sqrt(s / (len(x) - 1))

def dsto_mv(x):
  return dsto(x) / sqrt(len(x))

def dsys_mv(x):
  return sqrt(np.sum(x**2)) / len(x)

def dtot(dsys, dsto):
  return sqrt(dsys**2 + dsto**2)

def signval(val, err=0.0):
  if (err == 0.0):
    return ['{:g}'.format(val), '']
  if ('{:.1e}'.format(err)[0] == '9' and '{:.0e}'.format(err)[0] == '1'):
    err = float('{:.0e}'.format(err))
  firstdigit = int('{:.1e}'.format(err)[0])
  if (firstdigit <= 2):
    round2 = 1
    errstr = '{:.1e}'.format(err)
  else:
    round2 = 0
    errstr = '{:.0e}'.format(err)
  expdiff = int(np.floor(np.log10(abs(val))) - np.floor(np.log10(err)))
  if (expdiff < 0):
    sdigits = 0 #bug: try val('',0.999,1.1) and val('',1.0,1.1)
    if (round2 != 1 or expdiff != -1):
      val = 0.0
  else:
    sdigits = expdiff + round2
  valstr = '{:.{digits}e}'.format(val, digits=sdigits)
  return [valstr, errstr]

def val(name, val, err=0.0):
  out = ''
  if (name != ''):
    out += name + ' = '
  tmp = signval(val, err)
  out += tmp[0]
  if (tmp[1] != ''):
    out += ' ± ' + tmp[1]
  return out

def lst(val, err=[], name=''):
  """
  Parameters

  val: Array of floats with length N, which represent measured values
  err: Array of floats with length N, which represent the errors of the corresponding values
  name: String, which is added before the list
  ----------
  Returns

  Array of strings, where each string is in the form val[i] ± err[i]
  """
  if (err == []):
    err = [0.0 for i in range(len(val))]
  N = len(val)
  valmaxlen = 0
  errmaxlen = 0
  for i in range(N):
    tmp = signval(val[i], err[i])
    if (len(tmp[0]) > valmaxlen):
      valmaxlen = len(tmp[0])
    if (len(tmp[1]) > errmaxlen):
      errmaxlen = len(tmp[1])
  out = []
  if (name != ''):
    out.append(name)
  for i in range(len(val)):
    tmp = signval(val[i], err[i])
    tmp2 =  tmp[0].ljust(valmaxlen)
    if (tmp[1] != ''):
      tmp2 += ' ± ' + tmp[1].ljust(errmaxlen)
    elif (errmaxlen != 0):
      tmp2 += ''.ljust(errmaxlen + 3)
    out.append(tmp2)
  return out

def tbl(lists, title=''):
  """
  Parameters

  lists: Array of rowarrays with length N, which should be arrays with length M of the column strings.
  title: String, which is added before the table
  ----------
  Returns

  String of the MxN array.
  """
  M = len(lists[0])
  N = len(lists)
  out = ''
  if (title != ''):
    out += title + ':\n'
  lens = [int(npfarray([len(lists[i][j]) for j in range(M)]).max()) for i in range(N)]
  for j in range(M):
    for i in range(N):
      suffix = ' | '
      if (i == N-1):
        suffix = '\n' 
      out += lists[i][j].ljust(lens[i]) + suffix
  return out

def sig(name, val1, dVal1, val2, dVal2=0.0):
  nominator = abs(val1 - val2)
  denominator = np.sqrt(dVal1**2 + dVal2**2)
  if (nominator == 0.0):
    sigstr = '0'
  elif (denominator == 0.0):
    sigstr = '∞'
  else:
    sigma = nominator / denominator
    if (sigma < 0.95):
      digits = int(abs(np.floor(np.log10(sigma))))
    elif (sigma < 3.95):
      digits = 1
    else:
      digits = 0
    sigstr = '{:.{digits}f}'.format(sigma, digits = digits)
  prefix = ''
  if (name != ''):
    prefix = name + ': '
  return prefix + sigstr + 'σ'

def chi2(yo, dyo, ye, dye=[]):
  if (dye == []):
    dye = [0.0 for i in range(len(ye))]
  chi2 = 0.0
  for i in range(len(yo)):
    chi2 += (yo[i] - ye[i])**2 / (dyo[i]**2 + dye[i]**2)
  return chi2

def chi2_red(yo, dyo, ye, dye=[], dof=0):
  if (dof == 0):
    dof = len(ye)
  return chi2(yo, dyo, ye, dye) / dof

class pltext:
  @staticmethod
  def initplot(num=0, title='', xlabel='', ylabel='', scale='linlin', dpi=default_dpi):
    fig = plt.figure(num)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    fig.dpi = dpi
    fig.set_size_inches(11.69,8.27)
    plt.tight_layout()
    if (scale == 'linlin'):
      plt.ticklabel_format(style='sci', axis='both', scilimits=(-2,3))
    elif (scale == 'linlog'):
      plt.yscale('log')
      plt.ticklabel_format(style='sci', axis='x', scilimits=(-2,3))
    elif (scale == 'loglin'):
      plt.xscale('log')
      plt.ticklabel_format(style='sci', axis='y', scilimits=(-2,3))
    elif (scale == 'loglog'):
      plt.yscale('log')
      plt.xscale('log')

  @staticmethod
  def plotdata(x, y, dy=[], dx=[], label='', color=None, connect=False):
    if (dx == []):
      dx = [0.0 for i in range(len(x))]
    if (dy == []):
      dy = [0.0 for i in range(len(y))]
    plot = plt.errorbar(x=x, y=y, yerr=dy, xerr=dx, label=label, color=color, fmt='o', markersize=3)
    if (connect == True):
      if (color == None):
        color = plot[0].get_color()
      plt.plot(x, y, color=color)
    if (label != ''):
      plt.legend()

def linreg(x, y, dy, dx=[], plot=False, graphname=''):
  def linreg_iter(x, y, dy):
    [s0, s1, s2, s3, s4] = [0.0, 0.0, 0.0, 0.0, 0.0]
    for i in range(len(x)):
      if (dy[i] == 0.0):
        dy[i] = minfloat
      s0 += dy[i]**-2
      s1 += x[i] * dy[i]**-2
      s2 += y[i] * dy[i]**-2
      s3 += x[i]**2 * dy[i]**-2
      s4 += x[i] * y[i] * dy[i]**-2
    eta = s0 * s3 - s1**2
    g = (s0 * s4 - s1 * s2) / eta
    dg = np.sqrt(s0 / eta)
    b = (s3 * s2 - s1 * s4) / eta
    db = np.sqrt(s3 / eta)
    return [g, dg, b, db]

  iter0 = linreg_iter(x, y, dy)
  result = []
  if (dx == []):
    dx = [0.0 for i in range(len(x))]
    result = iter0
  else:
    g = iter0[0]
    g_old = g * (1 - 2 * linreg_change)
    while (abs(1 - g_old / g) >= linreg_change):
      g_old = g
      dy = [np.sqrt((g * dx[i])**2 + dy[i]**2) for i in range(len(dy))]
      g = linreg_iter(x, y, dy)[0]
    result = linreg_iter(x, y, dy)
  if (plot == True):
    [g, dg, b, db] = result
    min_x = np.argmin(x)
    max_x = np.argmax(x)
    xint = [x[min_x] - dx[min_x], x[max_x] + dx[max_x]]
    yfit = [g * xint[i] + b for i in range(2)]
    yerr = [(g + dg) * xint[i] + (b - db) for i in range(2)]
    fitfunc = plt.plot(xint, yfit, marker='')
    color = fitfunc[0].get_color()
    plt.plot(xint, yerr, marker='', linestyle='dashed', color=color)
    pltext.plotdata(x=x, y=y, dy=dy, dx=dx, label=graphname, color=color)
    if (graphname != ''):
      plt.legend()
  return result

def lin_yerr(x, dx, y, dy):
  g = linreg(x, y, dx, dy)
  new_dy = [np.sqrt(dy[i]**2 + (g * dx[i])**2) for i in range(len(dy))]
  return new_dy
