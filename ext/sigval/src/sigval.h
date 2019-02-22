#ifndef SIGVAL_H
#define SIGVAL_H

#include <python3.7/Python.h>

void _sigval(double val, double err, char* valstr, char* errstr, char* expstr);
static PyObject* sigval(PyObject* self, PyObject* args) {
  double val;
  double err;
  char valstr[0x40];
  char errstr[0x40];
  char expstr[0x40];

  if (!PyArg_ParseTuple(args, "dd", &val, &err)) {
    return NULL;
  }
  _sigval(val, err, valstr, errstr, expstr);
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