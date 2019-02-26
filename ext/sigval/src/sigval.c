#include "sigval.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

const char* digits = "0123456789";

double dround(double val, int d);
char* dtostr(double val, unsigned int d, char* valstr);
char* shiftPoint(char* val, int shift);
char* rmZerosbp(char* val);
int isDigit(char c);
int getDigit(char c);

void _sigval(double val, double err, char* valstr, char* errstr, char* expstr) {
  // Check for invalid uncertainty
  if (err < 0)
    return;
  // Check for sign
  int sign = val < 0;
  if (sign)
    val = -val;
  // Check for zero uncertainty
  if (err == 0.0) {
    const int d = 3;

    // get Exponent
    sprintf(valstr, "%.0e", val);
    int valExp = atoi(valstr + 2);

    // Round and get cutted string
    val = dround(val, -valExp + d);
    dtostr(val, d, valstr);

    // Shift
    int fixShift = 0;
    if (valExp > -3 && valExp < 3) {
      fixShift = valExp;
      shiftPoint(valstr, -fixShift);
      valExp = 0;
    }
    else {
      shiftPoint(valstr, 0);
    }

    // Cut the exponent
    valstr[2 + d] = '\0';
    errstr[0] = '\0';
    sprintf(expstr, "%i", valExp);

    // Add sign if val < 0 at the beginning
    if (sign && val != 0) {
      char val__[0x40];
      strcpy(val__, (char*)"-");
      strcat(val__, valstr);
      strcpy(valstr, val__);
    }
    return;
  }
  // Get exponents
  sprintf(valstr, "%.0e", val);
  sprintf(errstr, "%.0e", err);
  int valExp = atoi(valstr + 2);
  int errExp = atoi(errstr + 2);

  // Round to second nonzero place of err
  int d_exp = valExp - errExp;
  int shift = -valExp + d_exp + 1;
  val = dround(val, shift);
  err = dround(err, shift);

  // Get cutted strings
  int d = MAX(0, d_exp + 1);
  dtostr(val, d, valstr);
  dtostr(err, 1, errstr);

  // Round to first nonzero place of err
  if (getDigit(errstr[0]) > 2) {
    --shift; --d;
    val = dround(val, shift);
    err = dround(err, shift);
    dtostr(val, d, valstr);
    dtostr(err, 0, errstr);
  }

  // Shift errstr, so the exponend matches val
  int multDigitsbp_val;
  int multDigitsbp_err;
  int fixShift = 0;
  if (valExp > -3 && valExp < 3) {
    fixShift = valExp;
    shiftPoint(valstr, -fixShift);
    shiftPoint(errstr, -fixShift + d_exp);
    multDigitsbp_val = valExp > 0;
    multDigitsbp_err = fixShift - d_exp > 0;
    valExp = 0;
  }
  else {
    shiftPoint(valstr, 0);
    shiftPoint(errstr, d_exp);
    multDigitsbp_val = 0;
    multDigitsbp_err = -d_exp > 0;
  }

  // Cut the exponent and write it seperatly in expstr
  int c = 0;
  for (int i = 0; i < strlen(errstr); ++i) {
    if (getDigit(errstr[i]) > 0)
      c = i;
    else if (errstr[i] == '.')
      c = i - 1;
    else if (errstr[i] == 'e')
      break;
  }
  valstr[c + multDigitsbp_val * fixShift - multDigitsbp_err * (fixShift - d_exp) + 1] = '\0';
  errstr[c + 1] = '\0';
  sprintf(expstr, "%i", valExp);

  // Add sign if val < 0 at the beginning
  if (sign && val != 0) {
    char val__[0x40];
    strcpy(val__, (char*)"-");
    strcat(val__, valstr);
    strcpy(valstr, val__);
  }
}
void _sigval_fix(double val, double err, double fixExp, char* valstr, char* errstr, char* expstr) {
    // Check for invalid uncertainty
  if (err < 0)
    return;
  // Check for sign
  int sign = val < 0;
  if (sign)
    val = -val;
  // Check for zero uncertainty
  if (err == 0.0) {
    const int d = 3;

    // get Exponent
    sprintf(valstr, "%.0e", val);
    int valExp = atoi(valstr + 2);

    // Round and get cutted string
    val = dround(val, -valExp + d);
    dtostr(val, d, valstr);

    // Shift
    int fixShift = 0;
    if (valExp > -3 && valExp < 3) {
      fixShift = valExp;
      shiftPoint(valstr, -fixShift);
      valExp = 0;
    }
    else {
      shiftPoint(valstr, 0);
    }

    // Cut the exponent
    valstr[2 + d] = '\0';
    errstr[0] = '\0';
    sprintf(expstr, "%i", valExp);

    // Add sign if val < 0 at the beginning
    if (sign && val != 0) {
      char val__[0x40];
      strcpy(val__, (char*)"-");
      strcat(val__, valstr);
      strcpy(valstr, val__);
    }
    return;
  }
  // Get exponents
  sprintf(valstr, "%.0e", val);
  sprintf(errstr, "%.0e", err);
  int valExp = atoi(valstr + 2);
  int errExp = atoi(errstr + 2);

  // Round to second nonzero place of err
  int d_exp = valExp - errExp;
  int shift = -valExp + d_exp + 1;
  val = dround(val, shift);
  err = dround(err, shift);

  // Get cutted strings
  int d = MAX(0, d_exp + 1);
  dtostr(val, d, valstr);
  dtostr(err, 1, errstr);

  // Round to first nonzero place of err
  if (getDigit(errstr[0]) > 2) {
    --shift; --d;
    val = dround(val, shift);
    err = dround(err, shift);
    dtostr(val, d, valstr);
    dtostr(err, 0, errstr);
  }

  // Shift errstr, so the exponend matches val
  int multDigitsbp_val;
  int multDigitsbp_err;
  int fixShift = 0;
  fixShift = valExp - fixExp;
  shiftPoint(valstr, -fixShift);
  shiftPoint(errstr, -fixShift + d_exp);
  valExp -= fixShift;
  multDigitsbp_val = fixShift > 0;
  multDigitsbp_err = fixShift - d_exp > 0;

  // Cut the exponent and write it seperatly in expstr
  int c = 0;
  for (int i = 0; i < strlen(errstr); ++i) {
    if (getDigit(errstr[i]) > 0)
      c = i;
    else if (errstr[i] == '.')
      c = i - 1;
    else if (errstr[i] == 'e')
      break;
  }
  valstr[c + multDigitsbp_val * fixShift - multDigitsbp_err * (fixShift - d_exp) + 1] = '\0';
  errstr[c + 1] = '\0';
  sprintf(expstr, "%i", valExp);

  // Add sign if val < 0 at the beginning
  if (sign && val != 0) {
    char val__[0x40];
    strcpy(val__, (char*)"-");
    strcat(val__, valstr);
    strcpy(valstr, val__);
  }
}
void _sigval_fix_mul3(double val, double err, char* valstr, char* errstr, char* expstr) {
  // Check for invalid uncertainty
  if (err < 0)
    return;
  // Check for sign
  int sign = val < 0;
  if (sign)
    val = -val;
  // Check for zero uncertainty
  if (err == 0.0) {
    const int d = 3;

    // get Exponent
    sprintf(valstr, "%.0e", val);
    int valExp = atoi(valstr + 2);

    // Round and get cutted string
    val = dround(val, -valExp + d);
    dtostr(val, d, valstr);

    // Shift
    int fixShift = 0;
    if (valExp > -3 && valExp < 3) {
      fixShift = valExp;
      shiftPoint(valstr, -fixShift);
      valExp = 0;
    }
    else {
      shiftPoint(valstr, 0);
    }

    // Cut the exponent
    valstr[2 + d] = '\0';
    errstr[0] = '\0';
    sprintf(expstr, "%i", valExp);

    // Add sign if val < 0 at the beginning
    if (sign && val != 0) {
      char val__[0x40];
      strcpy(val__, (char*)"-");
      strcat(val__, valstr);
      strcpy(valstr, val__);
    }
    return;
  }
  // Get exponents
  sprintf(valstr, "%.0e", val);
  sprintf(errstr, "%.0e", err);
  int valExp = atoi(valstr + 2);
  int errExp = atoi(errstr + 2);

  // Round to second nonzero place of err
  int d_exp = valExp - errExp;
  int shift = -valExp + d_exp + 1;
  val = dround(val, shift);
  err = dround(err, shift);

  // Get cutted strings
  int d = MAX(0, d_exp + 1);
  dtostr(val, d, valstr);
  dtostr(err, 1, errstr);

  // Round to first nonzero place of err
  if (getDigit(errstr[0]) > 2) {
    --shift; --d;
    val = dround(val, shift);
    err = dround(err, shift);
    dtostr(val, d, valstr);
    dtostr(err, 0, errstr);
  }

  // Shift errstr, so the exponend matches val
  int multDigitsbp_val;
  int multDigitsbp_err;
  int fixShift = 0;
  fixShift = (valExp % 3 + (valExp < 0 ? 3 : 0)) % 3;
  shiftPoint(valstr, -fixShift);
  shiftPoint(errstr, -fixShift + d_exp);
  valExp -= fixShift;
  multDigitsbp_val = fixShift > 0;
  multDigitsbp_err = fixShift - d_exp > 0;

  // Cut the exponent and write it seperatly in expstr
  int c = 0;
  for (int i = 0; i < strlen(errstr); ++i) {
    if (getDigit(errstr[i]) > 0)
      c = i;
    else if (errstr[i] == '.')
      c = i - 1;
    else if (errstr[i] == 'e')
      break;
  }
  valstr[c + multDigitsbp_val * fixShift - multDigitsbp_err * (fixShift - d_exp) + 1] = '\0';
  errstr[c + 1] = '\0';
  sprintf(expstr, "%i", valExp);

  // Add sign if val < 0 at the beginning
  if (sign && val != 0) {
    char val__[0x40];
    strcpy(val__, (char*)"-");
    strcat(val__, valstr);
    strcpy(valstr, val__);
  }
}

double dround(double val, int d) {
  long shiftFactor;
  if (d < 0) {
    shiftFactor = pow(10, -d);
    val = round(val / shiftFactor) * shiftFactor;
  }
  else {
    shiftFactor = pow(10, d);
    val = round(val * shiftFactor) / shiftFactor;
  }
  return val;
}
char* dtostr(double val, unsigned int d, char* valstr) {
  char prec[0x40];
  char format[0x40];

  sprintf(prec, "%i", d);
  strcpy(format, "%.");
  strcat(format, prec);
  strcat(format, "e");
  sprintf(valstr, format, val);
  return valstr;
}
char* shiftPoint(char* val, int shift) {
  // Split val and exp
  char* valExp[2];
  valExp[0] = strtok(val, "e");
  valExp[1] = strtok(NULL, "e");

  // Split val at point
  char* valSplit[2];
  valSplit[0] = strtok(val, ".");
  valSplit[1] = strtok(NULL, ".");
  if (valSplit[1] == NULL)
    valSplit[1] = (char*)"";

  // Create new val string without point
  char newVal[0x40];
  int lenbp = strlen(valSplit[0]);
  int lenap = strlen(valSplit[1]);
  int iMin = MIN(lenbp - shift - 1, 0);
  int iMaxPlus1 = MAX(-shift - lenap + 1, 0) + lenbp + lenap;
  // Add zeros at the beginning
  for (int i = 0; i < -iMin; ++i) {
    newVal[i] = digits[0];
  }
  // Add the digits
  newVal[-iMin] = '\0';
  strcat(newVal, valSplit[0]);
  strcat(newVal, valSplit[1]);
  // Add zeros at the end
  for (int i = lenbp + lenap - iMin; i < iMaxPlus1 - iMin; ++i) {
    newVal[i] = digits[0];
  }
  newVal[iMaxPlus1 - iMin] = '\0';

  // Add point to string
  int p = -iMin + lenbp - shift;
  newVal[strlen(newVal) + 1] = '\0';
  for (int i = strlen(newVal); i >= 0; --i) {
    if (i == p) {
      newVal[i] = '.';
      break;
    }
    newVal[i] = newVal[i - 1];
  }

  // Add exponent
  int ePos = strlen(newVal);
  int newExp = (valExp[1] == NULL ? 0 : atoi(valExp[1])) + shift;
  if (newExp != 0) {
    newVal[ePos] = 'e';
    sprintf(newVal + ePos + 1, "%i", newExp);
  }
  strcpy(val, newVal);
  return val;
}
char* rmZerosbp(char* val) {
  // Check for sign
  int sign = val[0] == '-';
  int i = sign ? 1 : 0;
  
  // Find first nonzero digit
  while (val[i] == digits[0]) { ++i; };

  // Remove zeros before point
  char valstr[0x40];
  if (val[i] == '.') {
    // Remove all zeros but one before point
    strcpy(valstr, (char*)(sign ? "-0" : "0"));
    strcat(valstr, val + i);
  }
  else if (i == strlen(val)) {
    // All digits were zeros. Leave the last one
    strcpy(valstr, &val[i - 1]);
  }
  else {
    strcpy(valstr, (char*)(sign ? "-" : ""));
    strcat(valstr, val + i);
  }
  strcpy(val, valstr);
  return val;
}
int isDigit(char c) {
  for (unsigned int i = 0; i < strlen(digits); ++i) {
    if (digits[i] == 'c') {
      return 1;
    }
  }
  return 0;
}
int getDigit(char c) {
  for (unsigned int i = 0; i < strlen(digits); ++i) {
    if (c == digits[i]) {
      return i;
    }
  }
  return -1;
}
