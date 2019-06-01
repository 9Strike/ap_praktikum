import numpy as np
from sigval import sigval

# Unit prefixes for val and lst
unitPrefixes = "kMGTPEZYyzafpnμm"

# Table chars
singleFrameChars = ['│', '─', '┼', '┌', '┐', '└', '┘', '├', '┬', '┤' ,'┴']
doubleFrameChars = ['║', '═', '╬', '╔', '╗', '╚', '╝', '╠', '╦', '╣', '╩']
tableChars = singleFrameChars

def val(name, val, err=0.0, syserr=0.0, unit='', prefix=True):
  """
  Parameters

  val: float
  err: float, uncertainty of val
  name: string, name of val
  ----------
  Returns

  string, format: "name = val ± err" with two significant digits
  """

  if syserr != 0.0 and err == 0.0:
    raise ValueError('If syserr is specified one must also specify err')

  if err < 0.0:
    raise ValueError('The Uncertainty must be greater than zero')

  if abs(val) < err:
    print('Warning: The Uncertainty is greater than the value itself.')

  out = ''
  if name != '':
    out += name + ' = '
  
  syserrstr = None
  if syserr != 0.0:
    if syserr > err:
      valstr, syserrstr, expstr = sigval(val, syserr, unit != '' and prefix)
      exp = int(expstr)
      _, errstr, _ = sigval(val, err, True, exp)
    else:
      valstr, errstr, expstr = sigval(val, err, unit != '' and prefix)
      exp = int(expstr)
      _, syserrstr, _ = sigval(val, syserr, True, exp)
  else:
    valstr, errstr, expstr = sigval(val, err, unit != '' and prefix)

  if err != 0.0 and (expstr[0] != '0' or unit != ''):
    out += '('
  out += valstr
  if err != 0.0:
    out += ' ± ' + errstr
  if syserr != 0.0:
    out += ' stat. ± ' + syserrstr + ' syst.'
  if err != 0.0 and (expstr[0] != '0' or unit != ''):
    out += ')'
  if expstr[0] != '0':
    exp = int(expstr)
    if unit != '' and prefix and abs(exp) <= 3 * len(unitPrefixes) / 2:
      p = exp // 3
      if p > 0:
        p -= 1
      out += ' ' + unitPrefixes[p] + unit
    else:
      out += 'e' + expstr
      if unit != '':
        out += ' ' + unit
  else:
    out += ' ' + unit

  return out

def lst(val, err=[], name='', unit='', prefix=True, expToFix=None):
  """
  Parameters

  val: array of floats with length N
  err: array of floats with length N, uncertainties of val
  name: string, name of the list
  ----------
  Returns

  array of strings, format "val[i] ± err[i]" with significant digits
  """
  # Use zeros in case of empty err
  if (err == []):
    err = [0.0 for i in range(len(val))]

  # Use most frequent exponent (multiple of 3)
  N = len(val)
  lstExp = expToFix
  if expToFix == None or prefix:
    exps = np.zeros(N)
    for i in range(N):
      _, _, exps[i] = sigval(val[i], err[i], True)
    exps, counts = np.unique(exps, return_counts=True)
    lstExp = int(exps[np.argmax(counts)])

  # Determine maximal val and err lengths
  valmaxlen = 0
  errmaxlen = 0
  for i in range(N):
    tmp = sigval(val[i], err[i], True, lstExp)
    if (len(tmp[0]) > valmaxlen):
      valmaxlen = len(tmp[0])
    if (len(tmp[1]) > errmaxlen):
      errmaxlen = len(tmp[1])
  colWidth = valmaxlen + errmaxlen + 3 if errmaxlen > 0 else valmaxlen
  
  # Create title, center title and write to out
  out = []
  title = ''
  if name != '':
    title += name
  if unit != '' and lstExp != 0:
    title += ' / '
    if prefix:
      p = lstExp // 3
      uPrefix = ''
      if p > 0:
        uPrefix = unitPrefixes[p - 1]
      elif p < 0:
        uPrefix = unitPrefixes[p]
      title += uPrefix + unit
    else:
      title += '(' + 'e' + str(lstExp) + ' ' + unit + ')'
  elif unit != '':
    title += ' / ' + unit
  elif lstExp != 0:
    title += ' / ' + 'e' + str(lstExp)
  colWidth = max(colWidth, len(title))
  adjust = (colWidth + len(title)) // 2
  out.append(title.rjust(adjust))

  # Write and adjust value error strings to out
  for i in range(len(val)):
    tmp = sigval(val[i], err[i], True, lstExp)
    entry = tmp[0].rjust(valmaxlen)
    if (tmp[1] != ''):
      entry += ' ± ' + tmp[1].ljust(errmaxlen)
    elif (errmaxlen != 0):
      entry += ''.ljust(errmaxlen + 3)
    adjust = (colWidth + len(entry)) // 2
    out.append(entry.rjust(adjust))
  
  return out

def tbl(lists, name='', endl=True):
  """
  Parameters

  lists: array of rowarrays with length N, which should be arrays with length M of the column strings
  name: string, which is added before the table
  ----------
  Returns

  string of the MxN array
  """
  out = ''
  colWidths = [max([len(lists[i][j]) for j in range(len(lists[i]))]) for i in range(len(lists))]
  titles = [lists[i][0] for i in range(len(lists))]
  cols = [lists[i][1:] for i in range(len(lists))]
  nRows = len(cols[0])
  nCols = len(cols)

  # Print column titles
  for i in range(len(titles) - 1):
    out += titles[i].ljust(colWidths[i]) + ' ' + tableChars[0] + ' '
  out += titles[-1].ljust(colWidths[-1]) + '\n'

  # Print crossbar
  for i in range(len(titles) - 1):
    out += tableChars[1] * colWidths[i] + tableChars[1] + tableChars[2] + tableChars[1]
  out += tableChars[1] * colWidths[-1] + tableChars[1] + '\n'

  # Print tabular rows, by column entries
  for j in range(nRows - 1):
    for i in range(nCols - 1):
      out += cols[i][j].ljust(colWidths[i]) + ' ' + tableChars[0] + ' '
    out += cols[-1][j].ljust(colWidths[-1]) + '\n'
  for i in range(nCols - 1):
    out += cols[i][-1].ljust(colWidths[i]) + ' ' + tableChars[0] + ' '
  out += cols[-1][-1].ljust(colWidths[-1])

  # Connect extra column, which might be generated by dev
  rows = out.split('\n')
  for i in range(1, len(rows)):
    inds = [s for s in range(len(rows[i])) if rows[i][s] == tableChars[0]]
    for s in inds:
      upperRow = list(rows[i - 1])
      if upperRow[s] == tableChars[1]:
        upperRow[s] = tableChars[8]
        rows[i - 1] = ''.join(upperRow)
  out = ''
  for i in range(len(rows) - 1):
    out += rows[i] + '\n'
  out += rows[-1]

  # Print subtitle
  if (name != ''):
    out += '\n' + name
  return out + ('\n' if endl else '')

def sig(name, val1, err1, val2, err2=0.0, perc=False):
  ### deprecated, use dev instead
  return dev(val1,err1,val2,err2,name=name,perc=perc)

def dev(val1, err1, val2, err2=0.0, name='', perc=False):
  # Returns deviation string
  def get_dev(nom, denom):
    if (nom == 0.0):
      sigstr = '0'
    elif (denom == 0.0):
      sigstr = '∞ '
    else:
      sigma = nom / denom
      if (sigma < 0.95):
        digits = int(abs(np.floor(np.log10(sigma))))
      elif (sigma < 3.95):
        digits = 1
      else:
        digits = 0
      sigstr = '{:.{digits}f}'.format(sigma,digits=digits)
    sigstr += 'σ'
    return sigstr

  # Returns percental deviation string
  def get_perc(val1,val2,pformat='{:.2f}'):
    percval = abs(val1 - val2) / abs(val2) * 100
    percstr = pformat.format(percval) + '%'
    return percstr

  # Gets deviation of the difference from zero
  out = None
  nom = abs(val1 - val2)
  denom = np.sqrt(err1**2 + err2**2)
  if type(nom) is np.ndarray or type(denom) is np.ndarray:
    # Reconcile argument types
    out = []
    N = len(val1)
    if type(val1) is not np.ndarray:
      val1 = np.array([val1] * N)
    if type(val2) is not np.ndarray:
      val2 = np.array([val2] * N)
    if type(err1) is not np.ndarray:
      err1 = np.array([err1] * N)
    if type(err2) is not np.ndarray:
      err2 = np.array([err2] * N)
    if type(nom) is not np.ndarray:
      nom = np.array([nom] * N)
    if type(denom) is not np.ndarray:
      denom = np.array([denom] * N)
    
    if perc:
      # Get deviation and percental deviation strings and determine max lengths
      devs = []
      percs = []
      devmaxlen = 0
      percmaxlen = 0
      for i in range(N):
        # Get deviation string
        devs.append(get_dev(nom[i], denom[i]))
        siglen = len(devs[i])
        if (siglen > devmaxlen):
          devmaxlen = siglen
        # Get percental deviation string
        percs.append(get_perc(val1[i], val2[i]))
        perclen = len(percs[i])
        if (perclen > percmaxlen):
          percmaxlen = perclen
      colWidth = devmaxlen + 3 + percmaxlen if percmaxlen > 0 else devmaxlen

      if (name != ''):
        # Center name and write to out
        colWidth = max(colWidth, len(name))
        adjust = (colWidth + len(name)) // 2
        out.append(name.rjust(adjust))
      for i in range(N):
        # Center entry and write to out
        entry = devs[i].rjust(devmaxlen) + ' ' + tableChars[0] + ' ' + percs[i].rjust(percmaxlen)
        adjust = (colWidth + len(entry)) // 2
        out.append(entry.rjust(adjust))
    else:
      devs = []
      devmaxlen = 0
      for i in range(N):
        # Get deviation string
        devs.append(get_dev(nom[i], denom[i]))
        siglen = len(devs[i])
        if (siglen > devmaxlen):
          devmaxlen = siglen
      colWidth = devmaxlen

      if (name != ''):
        # Center name and write to out
        colWidth = max(colWidth, len(name))
        adjust = (colWidth + len(name)) // 2
        out.append(name.rjust(adjust))
      for i in range(N):
        # Center entry and write to out
        out.append(get_dev(nom[i], denom[i]).rjust(colWidth))
  else:
    out = ''
    prefix = ''
    if (name != ''):
      prefix = name + ': '
    out += prefix + get_dev(nom, denom)
    if perc:
      out += ' ≙ ' + get_perc(val1, val2, pformat='{:.2g}')
  return out
