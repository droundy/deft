// Deft is a density functional package developed by the research
// group of Professor David Roundy
//
// Copyright 2010 The Deft Authors
//
// Deft is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with deft; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
// Please see the file AUTHORS for a list of authors.

#include "Expression.h"
#include <sstream>
#include <iostream>
#include <cassert>
#include <stdio.h>
#include <stdlib.h>

Expression::Expression() {
  name = "UNKNOWN";
  kind = "NOKIND";
  type = "Grid";
  arg1 = arg2 = arg3 = 0;
}

Expression::Expression(const Expression &e) {
  arg1 = 0;
  arg2 = 0;
  arg3 = 0;
  *this = e;
}

void Expression::operator=(const Expression &e) {
  name = e.name;
  value = e.value;
  kind = e.kind;
  type = e.type;
  delete arg1;
  delete arg2;
  if (e.arg1) arg1 = new Expression(*e.arg1);
  else arg1 = 0;
  if (e.arg2) arg2 = new Expression(*e.arg2);
  else arg2 = 0;
  if (e.arg3) arg3 = new Expression(*e.arg3);
  else arg3 = 0;
}

Expression::Expression(std::string n) {
  arg1 = 0;
  arg2 = 0;
  arg3 = 0;
  name = n;
  type = "Grid";
  kind = "variable";
}

Expression::Expression(double c) {
  arg1 = arg2 = arg3 = 0;
  std::ostringstream oss;
  oss.precision(16);
  oss << c;
  name = oss.str();
  kind = "constant";
  type = "double";
  value = c;
}

Expression Expression::method(const char *n) const {
  Expression out;
  out.kind = "method";
  out.name = "." + std::string(n) + "(";
  out.arg1 = new Expression(*this);
  out.type = type; // default to methods not changing types
  return out;
}

Expression Expression::method(const char *n, const Expression &a) const {
  Expression out = method(n);
  out.arg2 = new Expression(a);
  return out;
}

Expression Expression::method(const char *n, const Expression &a, const Expression &b) const {
  Expression out = method(n, a);
  out.arg3 = new Expression(b);
  return out;
}

Expression Expression::operator()(const Expression &e) const {
  Expression out;
  out.name = "(";
  out.kind = "method";
  out.arg1 = new Expression(*this);
  out.arg2 = new Expression(e);
  return out;
}

Expression Expression::operator+(const Expression &e) const {
  if (kind == "constant" && value == 0) {
    return e;
  } else if (e.kind == "constant" && e.value == 0) {
    return *this;
  } else if (e.kind == "unary" && e.name == "-") {
    return *this - *e.arg1;
  } else if (kind == "unary" && name == "-") {
    return *e.arg1 - *this;
  } else if (e.kind == "constant" && kind == "constant") {
    return Expression(value+e.value);
  }
  Expression out;
  out.name = "+";
  out.kind = "+-";
  out.arg1 = new Expression(*this);
  out.arg2 = new Expression(e);
  return out;
}

Expression Expression::operator-(const Expression &e) const {
  if (kind == "constant" && value == 0) {
    return e;
  } else if (e.kind == "constant" && e.value == 0) {
    return *this;
  } else if (e.kind == "constant" && kind == "constant") {
    return Expression(value+e.value);
  }
  Expression out;
  out.name = "-";
  out.kind = "+-";
  out.arg1 = new Expression(*this);
  out.arg2 = new Expression(e);
  return out;
}

Expression Expression::operator-() const {
  if (kind == "constant") {
    return -value;
  } else if (kind == "*/" && arg1->kind == "constant") {
    return (-arg1->value) * *arg2;
  } else if (kind == "*/" && arg2->kind == "constant") {
    return *arg1 * (-arg2->value);
  }
  Expression out;
  out.name = "-";
  out.kind = "unary";
  out.arg1 = new Expression(*this);
  return out;
}

Expression Expression::operator*(const Expression &e) const {
  Expression larg = *this;
  Expression rarg = e;
  if (rarg.kind == "constant") {
    // prefer to have scalar on left.
    //rarg = e;
    //larg = *this;
  }
  // First, let's make a few optimizations...
  if (larg.kind == "constant" && larg.value == 1) {
    return rarg;
  } else if (rarg.kind == "constant" && rarg.value == 1) {
    return larg;
  } else if (larg.kind == "constant" && larg.value == -1) {
    return -rarg;
  } else if (rarg.kind == "constant" && rarg.value == -1) {
    return -larg;
  } else if (larg.kind == "constant" && larg.value == 0) {
    return larg;
  } else if (rarg.kind == "constant" && rarg.value == 0) {
    return rarg;
  } else if (rarg.kind == "constant" && larg.kind == "constant") {
    return Expression(larg.value*rarg.value);
  }
  Expression out;
  if (larg.type == "ReciprocalGrid") {
    if (rarg.type == "ReciprocalGrid") larg = larg.method("cwise");
    assert(rarg.type != "Grid");
    out.type = "ReciprocalGrid";
  } else if (larg.type == "Grid") {
    if (rarg.type == "Grid") larg = larg.method("cwise");
    assert(rarg.type != "ReciprocalGrid");
  }
  out.name = "*";
  out.kind = "*/";
  out.arg1 = new Expression(larg);
  out.arg2 = new Expression(e);
  return out;
}

Expression Expression::operator/(const Expression &e) const {
  // First, let's make a few optimizations...
  if (kind == "constant" && value == 0) {
    return *this;
  } else if (e.kind == "constant" && e.value == 1) {
    return *this;
  } else if (e.kind == "constant" && e.value == -1) {
    return - *this;
  } else if (e.kind == "constant" && kind == "constant") {
    return Expression(value/e.value);
  }
  Expression larg = *this;
  Expression rarg = e;
  Expression out;
  if (larg.type == "ReciprocalGrid") {
    if (rarg.type == "ReciprocalGrid") larg = larg.method("cwise");
    assert(rarg.type != "Grid");
    out.type = "ReciprocalGrid";
  } else if (larg.type == "Grid") {
    if (rarg.type == "Grid") larg = larg.method("cwise");
    assert(rarg.type != "ReciprocalGrid");
  }
  out.name = "/";
  out.kind = "*/";
  out.arg1 = new Expression(*this);
  out.arg2 = new Expression(e);
  return out;
}

Expression funexpr(const char *n) {
  Expression out;
  out.name = std::string(n) + "("; // slightly weird...
  out.kind = "function";
  return out;
}

Expression funexpr(const char *n, const Expression &arg) {
  Expression out = funexpr(n);
  out.arg1 = new Expression(arg);
  return out;
}

Expression funexpr(const char *n, const Expression &arg, const Expression &a2) {
  Expression out = funexpr(n, arg);
  out.arg2 = new Expression(a2);
  return out;
}

Expression funexpr(const char *n, const Expression &arg, const Expression &a2, const Expression &a3) {
  Expression out = funexpr(n, arg, a2);
  out.arg3 = new Expression(a3);
  return out;
}

std::string Expression::printme() const {
  if (kind == "+-") {
    // Addition
    std::string a1 = arg1->printme();
    std::string a2 = arg2->printme();
    if (arg2->kind == "+-") a2 = "(" + a2 + ")";
    return a1 + " " + name + " " + a2;
  } else if (kind == "*/") {
    // Multiplication
    std::string a1 = arg1->printme();
    if (arg1->kind == "+-") a1 = "(" + a1 + ")";
    std::string a2 = arg2->printme() ;
    if (arg2->kind == "+-" || arg2->kind == "*/") a2 = "(" + a2 + ")";
    return a1 + name + a2;
  } else if (kind == "unary" && name == "-") {
    std::string arg = arg1->printme();
    if (arg1->kind == "+-") arg = "(" + arg + ")";
    return "-" + arg;
  } else if (kind == "function") {
    // Function calls
    std::string out = name;
    if (arg1) {
      out += arg1->printme();
      if (arg2) {
        out += ", " + arg2->printme();
        if (arg3) out += ", " + arg3->printme();
      }
    }
    out += ")";
    return out;
  } else if (kind == "method") {
    std::string out = arg1->printme();
    if (arg1->kind == "+-" || arg1->kind == "*/" || arg1->kind == "unary") out = "(" + out + ")";
    out += name;
    if (arg2) {
      out += arg2->printme();
      if (arg3) out += ", " + arg3->printme();
    }
    out += ")";
    return out;
  }
  return name;
}

Expression fft(const Expression &g) {
  if (g.type != "Grid") {
    printf("fft: Expression %s should have type Grid.\n", g.printme().c_str());
    exit(1);
  }
  Expression out = funexpr("fft", Expression("gd"), g);
  out.type = "ReciprocalGrid";
  return out;
}

Expression ifft(const Expression &g) {
  if (g.type != "ReciprocalGrid") {
    printf("ifft: Expression '%s' should have type ReciprocalGrid but instead has type '%s'.\n",
           g.printme().c_str(), g.type.c_str());
    exit(1);
  }
  Expression out = funexpr("ifft", Expression("gd"), g);
  out.type = "Grid";
  return out;
}
