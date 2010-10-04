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
  alias = "";
  kind = "NOKIND";
  type = "Grid";
  arg1 = arg2 = arg3 = 0;
  unlazy = false;
}

Expression::Expression(const Expression &e) {
  arg1 = 0;
  arg2 = 0;
  arg3 = 0;
  *this = e;
}

void Expression::operator=(const Expression &e) {
  name = e.name;
  alias = e.alias;
  value = e.value;
  kind = e.kind;
  type = e.type;
  unlazy = e.unlazy;
  delete arg1;
  delete arg2;
  delete arg3;
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
  alias = "";
  type = "Grid";
  kind = "variable";
  unlazy = false;
}

Expression::Expression(double c) {
  arg1 = arg2 = arg3 = 0;
  std::ostringstream oss;
  oss.precision(16);
  oss << c;
  name = oss.str();
  alias = "";
  kind = "constant";
  type = "double";
  value = c;
  unlazy = false;
}

Expression Expression::method(const char *n) const {
  Expression out;
  out.kind = "method";
  out.name = "." + std::string(n) + "(";
  out.arg1 = new Expression(*this);
  out.type = type; // default to methods not changing types
  // check for special cases...
  if (out.name == ".square(" && type == "double") {
    return *this * *this;
  }
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

Expression Expression::operator()(const Expression &e, const Expression &f) const {
  Expression out;
  out.name = "(";
  out.kind = "method";
  out.arg1 = new Expression(*this);
  out.arg2 = new Expression(e);
  out.arg3 = new Expression(f);
  return out;
}

Expression Expression::set_alias(std::string a) {
  if (kind == "constant") {
    kind = "variable";
    name = a;
  }
  alias = a;
  return *this;
}

Expression Expression::cwise() const {
  if (iscwise()) return *this;
  if (type == "double") return *this;
  if (name == ".cwise(") return *this;
  return method("cwise");
}

bool Expression::iscwise() const {
  return name == ".cwise(" || type == "double";
}

Expression grid_ones = Expression("VectorXd::Ones(gd.NxNyNz)");
Expression recip_ones = Expression("VectorXcd::Ones(gd.NxNyNzOver2)").set_type("ReciprocalGrid");

Expression Expression::operator+(const Expression &e) const {
  if (kind == "constant" && value == 0) {
    return e;
  } else if (e.kind == "constant" && e.value == 0) {
    return *this;
  } else if (e.kind == "constant" && e.value < 0) {
    return *this - (-e);
  } else if (e.kind == "unary" && e.name == "-") {
    return *this - *e.arg1;
  } else if (kind == "unary" && name == "-") {
    return e - *arg1;
  } else if (e.kind == "constant" && kind == "constant") {
    return Expression(value+e.value);
  }
  Expression out;
  if (type == "ReciprocalGrid") {
    assert(e.type != "Grid");
    if (e.type == "double") return *this + e*recip_ones;
    out.type = "ReciprocalGrid";
  } else if (type == "Grid") {
    assert(e.type != "ReciprocalGrid");
    if (e.type == "double") return *this + e*grid_ones;
  } else if (type == "double") {
    if (e.type == "ReciprocalGrid") return (*this)*recip_ones + e;
    else if (e.type == "Grid") return (*this)*grid_ones + e;
    assert(e.type == "double");
    out.type = "double";
  }
  out.name = "+";
  out.kind = "+-";
  out.arg1 = new Expression(*this);
  out.arg2 = new Expression(e);
  return out;
}

Expression Expression::operator-(const Expression &e) const {
  if (kind == "constant" && value == 0) {
    return -e; // 0 - b = -b
  } else if (e.kind == "constant" && e.value == 0) {
    return *this; // a - 0 = a
  } else if (e.kind == "constant" && e.value < 0) {
    return *this + (-e.value); // a - -3.0 = a + 3.0
  } else if (e.kind == "unary" && e.name == "-") {
    return *this + *e.arg1; // a - -b = a + b
  } else if (e.kind == "constant" && kind == "constant") {
    return Expression(value-e.value);
  }
  Expression out;

  if (type == "ReciprocalGrid") {
    assert(e.type != "Grid");
    if (e.type == "double") return *this - e*recip_ones;
    out.type = "ReciprocalGrid";
  } else if (type == "Grid") {
    assert(e.type != "ReciprocalGrid");
    if (e.type == "double") return *this - e*grid_ones;
  } else if (type == "double") {
    if (e.type == "ReciprocalGrid") {
      return (*this)*recip_ones - e;
    } else if (e.type == "Grid") {
      return (*this)*grid_ones - e;
    }
    assert(e.type == "double");
    out.type = "double";
  }
  out.name = "-";
  out.kind = "+-";
  out.arg1 = new Expression(*this);
  out.arg2 = new Expression(e);
  out.checkWellFormed();
  return out;
}

Expression Expression::operator-() const {
  checkWellFormed();
  if (kind == "constant") {
    return -value;
  } else if (kind == "*/" && name == "*" && arg1->kind == "constant") {
    return (-arg1->value) * *arg2;
  } else if (kind == "*/" && name == "/" && arg1->kind == "constant") {
    return (-arg1->value) / *arg2;
  } else if (kind == "*/" && name == "*" && arg2->kind == "constant") {
    return *arg1 * (-arg2->value);
  } else if (kind == "*/" && name == "/" && arg2->kind == "constant") {
    return *arg1 / (-arg2->value);
  } else if (kind == "+-" && name == "+") {
    return (-*arg1) - *arg2;
  } else if (kind == "+/" && name == "-") {
    return (-*arg1) - (-*arg2);
  } else if (kind == "unary" && name == "-") {
    return *arg1;
  }
  Expression out;
  out.name = "-";
  out.kind = "unary";
  out.type = type;
  out.arg1 = new Expression(*this);
  out.checkWellFormed();
  return out;
}

Expression Expression::operator*(const Expression &e) const {
  checkWellFormed();
  e.checkWellFormed();
  if (kind == "constant" && e.kind == "constant") {
    return Expression(value*e.value);
  } else if (e.kind == "constant") {
    // prefer to have scalar on left.
    return e*(*this);
  }
  // First, let's make a few optimizations...
  if (kind == "constant" && value == 1) {
    return e;
  } else if (kind == "constant" && value == -1) {
    return -e;
  } else if (kind == "constant" && value == 0) {
    return *this;
  } else if (kind == "constant" && e.kind == "constant") {
    return Expression(value*e.value);
  }
  Expression out;
  if (e.type == "ReciprocalGrid") {
    assert(type != "Grid");
    out.type = "ReciprocalGrid";
  } else if (type == "ReciprocalGrid") {
    out.type = "ReciprocalGrid";
  } else if (e.type == "Grid") {
    assert(e.type != "ReciprocalGrid");
  } else if (e.type == "double" && type == "double") {
    out.type = "double";
  }
  out.name = "*";
  out.kind = "*/";
  out.arg1 = new Expression(*this);
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
  Expression out;
  if (type == "ReciprocalGrid") {
    assert(e.type != "Grid");
    out.type = "ReciprocalGrid";
  } else if (type == "Grid") {
    assert(e.type != "ReciprocalGrid");
  } else if (type == "double" && kind == "constant") {
    if (e.type == "Grid") return (*this)*grid_ones/e;
    if (e.type == "ReciprocalGrid") return (*this)*recip_ones/e;
    if (e.type == "double") out.type = "double";
  } else if (type == "ReciprocalGrid") {
    out.type = "ReciprocalGrid";
  } else if (e.type == "double" && type == "double") {
    out.type = "double";
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


Expression linearfunexpr(const char *n, const Expression &arg) {
  if (arg.kind == "*/") {
    if (arg.arg1->type == "double" && arg.name == "*") {
      return *arg.arg1 * linearfunexpr(n, *arg.arg2);
    } else if (arg.arg2->type == "double" && arg.name == "*") {
      return *arg.arg2 * linearfunexpr(n, *arg.arg1);
    } else if (arg.arg2->type == "double" && arg.name == "/") {
      return linearfunexpr(n, *arg.arg1) / *arg.arg2;
    }
  }
  Expression out = funexpr(n, arg);
  out.kind = "linear function";
  return out;
}

Expression linearfunexprgd(const char *n, const char *type, const Expression &arg) {
  if (arg.kind == "*/") {
    if (arg.arg1->type == "double" && arg.name == "*") {
      return *arg.arg1 * linearfunexprgd(n, type, *arg.arg2);
    } else if (arg.arg2->type == "double" && arg.name == "*") {
      return *arg.arg2 * linearfunexprgd(n, type, *arg.arg1);
    } else if (arg.arg2->type == "double" && arg.name == "/") {
      return linearfunexprgd(n, type, *arg.arg1) / *arg.arg2;
    }
  } else if (arg.kind == "unary" && arg.name == "-") {
    return - linearfunexprgd(n, type, *arg.arg1);
  }
  Expression out = funexpr(n, arg);
  out.kind = "linear function";
  if (arg.arg1 && arg.arg2)
    out.name = std::string(n) + "(gd, /* " + arg.arg1->type + "  " + arg.name + "  " + arg.arg2->type + " */";
  else
    out.name = std::string(n) + "(gd, /* " + arg.name + " */";
  out.name = std::string(n) + "(gd, ";
  out.type = type;
  return out;
}

Expression fft(const Expression &g) {
  if (g.type != "Grid") {
    printf("fft: Expression %s should have type Grid.\n", g.printme().c_str());
    exit(1);
  }
  Expression out = funexpr("fft", Expression("gd"), g);
  out.type = "ReciprocalGrid";
  out.unlazy = true;
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
  out.unlazy = true;
  return out;
}

std::string Expression::printme() const {
  if (kind == "+-") {
    // Addition and subtraction
    std::string a1 = arg1->printme();
    std::string a2 = arg2->printme();
    if (arg2->kind == "+-") a2 = "(" + a2 + ")";
    return a1 + " " + name + " " + a2;
  } else if (kind == "*/") {
    // Multiplication and division
    Expression myarg1 = *arg1;
    if (arg2->type != "double" && !arg1->iscwise()) myarg1 = myarg1.cwise();
    std::string a1 = myarg1.printme();
    if (myarg1.kind == "+-") a1 = "(" + a1 + ")";
    std::string a2 = arg2->printme() ;
    if (arg2->kind == "+-" || arg2->kind == "*/") a2 = "(" + a2 + ")";
    return a1 + name + a2;
  } else if (kind == "unary" && name == "-") {
    // Unary negation
    std::string arg = arg1->printme();
    if (arg1->kind == "+-" || arg1->kind == "*/") arg = "(" + arg + ")";
    return "-" + arg;
  } else if (kind == "function" || kind == "linear function") {
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
  } else if (kind == "constant" && alias != "") {
    return alias;
  }
  return name;
}

int Expression::checkWellFormed() const {
  int retval = 0;
  if (kind == "+-") {
    if (arg1->type != type || arg2->type != type) {
      printf("Types don't match on addition/subtraction!\n");
      retval++;
    }
  } else if (kind == "*/") {
    if (arg1->type != type && arg1->type != "double") {
      printf("Types don't match on multiplication!\n");
      retval++;
    }
    if (arg2->type != type && arg2->type != "double") {
      printf("Types don't match on multiplication!\n");
      retval++;
    }
  } else if (kind == "unary" && name == "-") {
    if (arg1->type != type) {
      printf("Types don't match on unary minus: %s\n", printme().c_str());
      printf("Type of %s is %s\n", arg1->printme().c_str(), arg1->type.c_str());
      printf("\n\n");
      retval++;
    }
  }
  if (arg1) retval += arg1->checkWellFormed();
  if (arg2) retval += arg2->checkWellFormed();
  if (arg3) retval += arg3->checkWellFormed();
  return retval;
}

Expression Expression::simplify() const {
  if (name == "*" && arg2->name == "*") {
    if (arg2->arg1->kind == "constant")
      return (*arg2->arg1 * *arg1 * *arg2->arg2).simplify();
    if (arg2->arg2->kind == "constant")
      return (*arg2->arg2 * *arg1 * *arg2->arg1).simplify();
    return (*arg1 * *arg2->arg1 * *arg2->arg2).simplify();
  }
  if (kind == "unary" && name == "-") {
    Expression minusthis = arg1->simplify();
    if (minusthis.kind == "*/" && minusthis.arg1->kind == "constant") {
      minusthis.arg1 = new Expression(-minusthis.arg1->value);
      return minusthis;
    }
    if (minusthis.kind == "*/" && minusthis.arg1->kind == "*/" && minusthis.arg1->arg1->kind == "constant") {
      minusthis.arg1->arg1 = new Expression(-minusthis.arg1->arg1->value);
      return minusthis;
    }
    return -minusthis;
  }
  Expression out;
  out.name = name;
  out.type = type;
  out.kind = kind;
  out.alias = alias;
  if (arg1) out.arg1 = new Expression(arg1->simplify());
  if (arg2) out.arg2 = new Expression(arg2->simplify());
  if (arg3) out.arg3 = new Expression(arg3->simplify());
  return *this;
}

void Expression::generate_code(FILE *o, const char *fmt, std::set<std::string> allvars) const {
  // Set of variables that need to be deleted.
  std::set<std::string> vars;

  Expression e = *this;
  Expression s = e.FindCommonSubexpression();
  while (s.unlazy) {
    std::string a = s.alias + "_v";
    if (a == "_v") a = "var";
    if (allvars.count(a)) {
      // need to make this thing unique...
      for (int varnum=0; true; varnum++) {
        std::ostringstream oss;
        oss << varnum;
        std::string newa = a + "_" + oss.str();
        if (!allvars.count(newa)) {
          a = newa;
          break;
        }
      }
    }
    fprintf(o, "    %s %s(%s);\n", s.ctype(), a.c_str(), s.printme().c_str());
    //fprintf(o, "    %s *%s_ptr = new %s(%s);\n", s.ctype(), a.c_str(), s.ctype(), s.printme().c_str());
    //fprintf(o, "    %s %s = *%s_ptr;\n", s.ctype(), a.c_str(), a.c_str());
    vars.insert(a);
    allvars.insert(a);
    while (e.EliminateThisSubexpression(s, a));
    s = e.FindCommonSubexpression();
    for (std::set<std::string>::iterator i = vars.begin(); i != vars.end(); ++i) {
      if (!e.FindVariable(*i)) {
        //fprintf(o, "    // Couldn't find %s in:  %s\n", i->c_str(), e.printme().c_str());
        //fprintf(o, "    delete %s_ptr;\n", i->c_str());
        fprintf(o, "    %s.resize(0); // We're done with this...\n", i->c_str());
        vars.erase(i);
      }
    }
  }
  fprintf(o, fmt, e.printme().c_str());
}
