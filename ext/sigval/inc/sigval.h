#ifndef SIGVAL_H
#define SIGVAL_H

void _sigval(double val, double err, char* valstr, char* errstr, char* expstr);
void _sigval_fix(double val, double err, double fixExp, char* valstr, char* errstr, char* expstr);
void _sigval_fix_mul3(double val, double err, char* valstr, char* errstr, char* expstr);

#endif
