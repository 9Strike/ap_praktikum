#include <iostream>

#include "sigval.h"

int main(int argc, char** argv) {
  char valstr[0x40];
  char errstr[0x40];
  char expstr[0x40];
  _sigval_fix(992.0, 5.0, 0.0, valstr, errstr, expstr);

  std::cout << valstr << " +- " << errstr << " e" << expstr << std::endl;
  return 0;
}
