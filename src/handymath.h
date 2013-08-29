// -*- mode: C++; -*-

#pragma once

#include <stdio.h>
#include "Faddeeva.hh"

inline double erfi(double x) {
  return Faddeeva::erfi(x);
}

inline double sqr(double x) { return x*x; }

inline double uipow(double x, unsigned int n) {
  double sofar = 1.0;
  while (n > 0) {
    if (n & 1) sofar *= x;
    x *= x;
    n = n >> 1;
  }
  return sofar;
}

inline double min(double a, double b) { return (a<b) ? a : b; }
inline double max(double a, double b) { return (a>b) ? a : b; }
inline double heaviside(double x) { return (x>0) ? 1.0 : 0.0; }

inline double ipow(double x, int n) {
  if (n >= 0) return uipow(x,n);
  else return 1.0/uipow(x,-n);
}

inline void print_double(const char *prefix, double x, int width=26, int digits = 14, int min_digits = 11) {
  double max_f, min_f;
  if (min_digits == 0 || min_digits > digits) min_digits = digits;
  width -= strlen(prefix);
  if (width < digits + 8) width = digits + 8;

  if (x > 0) {
    max_f = min(uipow(10, width - digits - 3), 1e6);
    min_f = max(ipow(10, min_digits - digits - 4), 1e-5);
  }
  else {
    min_f = max(-uipow(10, width - digits - 5), -1e6);
    max_f = min(-ipow(10, min_digits - digits - 4), -1e-5);
  }
  if ((x < max_f && x > min_f) || x == 0.0)
    printf("%s %*.*f", prefix, width, digits + 3, x);
  else printf("%s %*.*e", prefix, width, digits, x);
}
