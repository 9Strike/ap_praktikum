#ifndef SIGVAL_H
#define SIGVAL_H

#include <python3.7/Python.h>

void _sigval(double val, double err, char* valstr, char* errstr, char* expstr);
void _sigval_fix(double val, double err, double fixExp, char* valstr, char* errstr, char* expstr);
void _sigval_fix_mul3(double val, double err, char* valstr, char* errstr, char* expstr);

static PyObject* sigval(PyObject* self, PyObject* args) {
  double val;
  double err;
  int fix = 0;
  int fixExp = INT_MIN;
  char valstr[0x40];
  char errstr[0x40];
  char expstr[0x40];

  if (!PyArg_ParseTuple(args, "dd|ii", &val, &err, &fix, &fixExp)) {
    return NULL;
  }
  if (fix && fixExp != INT_MIN)
    _sigval_fix(val, err, fixExp, valstr, errstr, expstr);
  else if (fix)
    _sigval_fix_mul3(val, err, valstr, errstr, expstr);
  else _sigval(val, err, valstr, errstr, expstr);
  return Py_BuildValue("sss", valstr, errstr, expstr);
}
static PyMethodDef svMethods[] = {
  {"sigval", sigval, METH_VARARGS, "Displays the significant parts of a defective value"},
  {NULL, NULL, 0, NULL}
};
static struct PyModuleDef sv = {
  PyModuleDef_HEAD_INIT,
  "sv",
  "Parses a value and its error to a string with only the significant part",
  -1,
  svMethods
};
PyMODINIT_FUNC PyInit_sigval() {
  return PyModule_Create(&sv);
}

#endif 