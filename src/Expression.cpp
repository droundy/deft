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
  depth = 1;
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
  depth = e.depth;
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
  depth = 1;
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
  if (c == -0) value = 0; // I don't want to bother with minus zero...
  unlazy = false;
  depth = 1;
}

Expression Expression::method(const std::string n) const {
  Expression out;
  out.kind = "method";
  out.name = "." + n + "(";
  out.arg1 = new Expression(*this);
  out.type = type; // default to methods not changing types
  out.SetDepth();
  return out;
}

Expression Expression::method(const char *n, const Expression &a) const {
  Expression out = method(n);
  out.arg2 = new Expression(a);
  out.SetDepth();
  return out;
}

Expression Expression::method(const char *n, const Expression &a, const Expression &b) const {
  Expression out = method(n, a);
  out.arg3 = new Expression(b);
  out.SetDepth();
  return out;
}

Expression abs(const Expression &x) {
  if (x.typeIs("double")) return funexpr("fabs", x).set_type("double");
  return x.cwisemethod("abs");
}

Expression log(const Expression &x) {
  if (x.typeIs("double")) return funexpr("log", x).set_type("double");
  return x.cwisemethod("log");
}

Expression exp(const Expression &x) {
  if (x.typeIs("double")) return funexpr("exp", x).set_type("double");
  return x.cwisemethod("exp");
}

Expression sqrt(const Expression &x) {
  if (x.typeIs("double")) return funexpr("sqrt", x).set_type("double");
  return x.cwisemethod("sqrt");
}

Expression sqr(const Expression &x) {
  if (x.typeIs("double")) return x*x;
  return x.cwisemethod("square");
}

static inline int max(int a, int b) {
  return (a>b) ? a : b;
}

void Expression::SetDepth() {
  depth = 1;
  if (arg1) depth = max(depth, 1 + arg1->depth);
  if (arg2) depth = max(depth, 1 + arg2->depth);
  if (arg3) depth = max(depth, 1 + arg3->depth);
}

Expression Expression::operator()(const Expression &e) const {
  Expression out;
  out.name = "(";
  out.kind = "method";
  out.arg1 = new Expression(*this);
  out.arg2 = new Expression(e);
  out.SetDepth();
  return out;
}

Expression Expression::operator()(const Expression &e, const Expression &f) const {
  Expression out;
  out.name = "(";
  out.kind = "method";
  out.arg1 = new Expression(*this);
  out.arg2 = new Expression(e);
  out.arg3 = new Expression(f);
  out.SetDepth();
  return out;
}

Expression Expression::set_alias(std::string a) {
  if (kindIs("constant") && a != "literal") {
    kind = "variable";
    name = a;
  }
  alias = a;
  return *this;
}

Expression Expression::append_to_alias(const std::string a) {
  if (alias != "" && !kindIs("variable")) {
    alias = alias + a;
  }
  if (arg1) arg1->append_to_alias(a);
  if (arg2) arg2->append_to_alias(a);
  if (arg3) arg3->append_to_alias(a);
  return *this;
}

Expression Expression::cwise() const {
  if (iscwise()) return *this;
  if (typeIs("double")) return *this;
  if (name == ".cwise(") return *this;
  return method("cwise");
}

Expression Expression::cwisemethod(const char *m) const {
  return method("cwise()." + std::string(m));
}

bool Expression::iscwise() const {
  return name == ".cwise(" || typeIs("double");
}

Expression grid_ones = Expression("VectorXd::Ones(gd.NxNyNz)");
Expression recip_ones = Expression("VectorXcd::Ones(gd.NxNyNzOver2)").set_type("ReciprocalGrid");

Expression Expression::operator+(const Expression &e) const {
  if (kindIs("constant") && value == 0) {
    return e;
  } else if (e.kindIs("constant") && e.value == 0) {
    return *this;
  } else if (e.kindIs("constant") && kindIs("constant")) {
    return Expression(value+e.value);
  }
  Expression thispost(*this), epost(e);
  Expression thispre = thispost.ScalarFactor(), epre = epost.ScalarFactor();
  if (thispost.kindIs("linear function") && epost.kindIs("linear function") &&
      thispost.name == epost.name) {
    if (thispre == epre) {
      Expression out = linearfunexprgd("oops", thispost.type, *thispost.arg1 + *epost.arg1);
      out.name = thispost.name;
      return epre*out;
    }
    if (thispre == -epre && thispre != Expression(1)) {
      Expression out = linearfunexprgd("oops", thispost.type, *thispost.arg1 - *epost.arg1);
      out.name = thispost.name;
      return thispre*out;
    }
    Expression out = linearfunexprgd("oops", thispost.type,
                                     thispre * *thispost.arg1 + epre * *epost.arg1);
    out.name = thispost.name;
    return out;
  }
  Expression out;
  if (typeIs("ReciprocalGrid")) {
    assert(!e.typeIs("Grid"));
    out.type = "ReciprocalGrid";
  } else if (typeIs("Grid")) {
    assert(!e.typeIs("ReciprocalGrid"));
  } else if (e.typeIs("double")) {
    out.type = "double";
  } else if (e.typeIs("ReciprocalGrid")) {
    out.type = "ReciprocalGrid";
  }
  out.name = "+";
  out.kind = "+-";
  out.arg1 = new Expression(*this);
  out.arg2 = new Expression(e);
  out.SetDepth();
  return out;
}

Expression Expression::operator-(const Expression &e) const {
  return *this + (-e);
}

Expression Expression::operator-() const {
  // checkWellFormed();
  if (kindIs("constant")) {
    return -value;
  } else if (kindIs("*/") && name == "*" && arg1->kindIs("constant")) {
    return (-arg1->value) * *arg2;
  } else if (kindIs("*/") && name == "/" && arg1->kindIs("constant")) {
    return (-arg1->value) / *arg2;
  } else if (kindIs("*/") && name == "*" && arg2->kindIs("constant")) {
    return *arg1 * (-arg2->value);
  } else if (kindIs("*/") && name == "/" && arg2->kindIs("constant")) {
    return *arg1 / (-arg2->value);
  //} else if (kindIs("+-") && name == "+") {
  //  return (-*arg1) + (-*arg2);
  } else if (kindIs("unary") && name == "-") {
    return *arg1;
  }
  Expression out;
  out.name = "-";
  out.kind = "unary";
  out.type = type;
  out.arg1 = new Expression(*this);
  out.SetDepth();
  // out.checkWellFormed();
  return out;
}

Expression Expression::operator*(const Expression &e) const {
  // checkWellFormed();
  // e.checkWellFormed();
  if (kindIs("constant") && e.kindIs("constant")) {
    return Expression(value*e.value);
  } else if (e.kindIs("constant")) {
    // prefer to have scalar on left.
    return e*(*this);
  }
  // First, let's make a few optimizations...
  if (kindIs("constant") && value == 1) {
    return e;
  } else if (kindIs("constant") && value == -1) {
    return -e;
  } else if (kindIs("constant") && value == 0) {
    return *this;
  } else if (e.kindIs("unary") && e.name == "-") {
    return - (*this * *e.arg1);
  }
  Expression out;
  if (e.typeIs("ReciprocalGrid")) {
    assert(!typeIs("Grid"));
    out.type = "ReciprocalGrid";
  } else if (typeIs("ReciprocalGrid")) {
    out.type = "ReciprocalGrid";
  } else if (e.typeIs("Grid")) {
    assert(!e.typeIs("ReciprocalGrid"));
  } else if (e.typeIs("double") && typeIs("double")) {
    out.type = "double";
  }
  out.name = "*";
  out.kind = "*/";
  out.arg1 = new Expression(*this);
  out.arg2 = new Expression(e);
  out.SetDepth();
  return out;
}

Expression Expression::operator/(const Expression &e) const {
  // First, let's make a few optimizations...
  if (kindIs("constant") && value == 0) {
    return *this;
  } else if (e.kindIs("constant") && e.value == 1) {
    return *this;
  } else if (e.kindIs("constant") && e.value == -1) {
    return - *this;
  } else if (e.kindIs("constant") && kindIs("constant")) {
    return Expression(value/e.value);
  } else if (e.kindIs("*/") && e.name == "/") {
    return *this * *e.arg2 / *e.arg1; // inverse of inverse
  }
  //Expression postfactor = e;
  //Expression prefactor = postfactor.ScalarFactor();
  //if (postfactor.kindIs("*/") && postfactor.name == "/") {
  //  return (*this / prefactor) * *postfactor.arg2 / *postfactor.arg1;
  //}
  Expression out;
  if (typeIs("ReciprocalGrid")) {
    assert(!e.typeIs("Grid"));
    out.type = "ReciprocalGrid";
  } else if (typeIs("Grid")) {
    assert(!e.typeIs("ReciprocalGrid"));
  } else if (typeIs("ReciprocalGrid")) {
    out.type = "ReciprocalGrid";
  } else if (e.typeIs("double") && typeIs("double")) {
    out.type = "double";
  } else if (e.typeIs("ReciprocalGrid")) {
    out.type = "ReciprocalGrid";
  }
  out.name = "/";
  out.kind = "*/";
  out.arg1 = new Expression(*this);
  out.arg2 = new Expression(e);
  out.SetDepth();
  return out;
}

Expression Expression::ScalarFactor() {
  if (kindIs("constant")) {
    Expression out(value);
    value = 1;
    return out;
  }
  if (typeIs("double")) {
    // If we have type "double" then we *are* a scalar!
    Expression out(*this);
    *this = Expression(1);
    return out;
  }
  if (kindIs("*/") && name == "*") {
    Expression sf1 = arg1->ScalarFactor();
    Expression sf2 = arg2->ScalarFactor();
    *this = *arg1 * *arg2; // to handle cases where arg1 or arg2 becomes 1.
    return sf1*sf2;
  }
  if (kindIs("*/") && name == "/") {
    Expression sf1 = arg1->ScalarFactor();
    Expression sf2 = arg2->ScalarFactor();
    *this = *arg1 / *arg2; // to handle cases where arg1 or arg2 becomes 1.
    return sf1/sf2;
  }
  if (kindIs("linear function")) {
    return arg1->ScalarFactor();
  }
  if (kindIs("unary") && name == "-") {
    Expression postfactor = *arg1;
    Expression prefactor = postfactor.ScalarFactor();
    *this = postfactor;
    return -prefactor;
  }
  return Expression(1);
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
  out.SetDepth();
  return out;
}

Expression funexpr(const char *n, const Expression &arg, const Expression &a2) {
  Expression out = funexpr(n, arg);
  out.arg2 = new Expression(a2);
  out.SetDepth();
  return out;
}

Expression funexpr(const char *n, const Expression &arg, const Expression &a2, const Expression &a3) {
  Expression out = funexpr(n, arg, a2);
  out.arg3 = new Expression(a3);
  out.SetDepth();
  return out;
}

Expression linearfunexprgd(const char *n, const char *type, const Expression &arg) {
  if (arg.kindIs("constant") && arg.value == 0)
    return Expression(0); // Nice optimization!  :)
  Expression newarg(arg);
  Expression prefactor = newarg.ScalarFactor();

  Expression out = funexpr(n, newarg);
  out.kind = "linear function";
  out.name = std::string(n) + "(gd, ";
  out.unlazy = true;
  out.type = type;
  out.SetDepth();
  return prefactor*out;
}

Expression fft(const Expression &g) {
  if (!g.typeIs("Grid")) {
    if (g.typeIs("double") && g.kindIs("constant") && g.value == 0) {
      return g;
    } else if (g.typeIs("double")) {
      // The fft of a constant is a delta function in G space!
      printf("fft: Expression should have type Grid, but accepting double anyhow.\n");
      printf("fft: Expression is:  %s\n", g.printme().c_str());
      // FIXME: We shouldn't need to do an actual FFT here!
      return fft(g * grid_ones);
    }
    printf("fft: Expression %s should have type Grid.\n", g.printme().c_str());
    exit(1);
  }
  return linearfunexprgd("fft", "ReciprocalGrid", g);
}

Expression ifft(const Expression &g) {
  if (!g.typeIs("ReciprocalGrid")) {
    if (g.typeIs("double") && g.kindIs("constant") && g.value == 0) {
      return g;
    } else {
      printf("ifft: Expression '%s' should have type ReciprocalGrid but instead has type '%s'.\n",
             g.printme().c_str(), g.type);
      exit(1);
    }
  }
  return linearfunexprgd("ifft", "Grid", g);
}

std::string Expression::printme() const {
  if (kindIs("+-")) {
    // Addition and subtraction
    std::string a1 = arg1->printme();
    std::string a2 = arg2->printme();
    if (arg1->typeIs("double") && arg2->typeIs("Grid"))
      a1 = (*arg1 * grid_ones).printme();
    if (arg1->typeIs("double") && arg2->typeIs("ReciprocalGrid"))
      a1 = (*arg1 * recip_ones).printme();
    if (arg2->typeIs("double") && arg1->typeIs("Grid"))
      a2 = (*arg2 * grid_ones).printme();
    if (arg2->typeIs("double") && arg1->typeIs("ReciprocalGrid"))
      a2 = (*arg2 * recip_ones).printme();
    if (arg2->kindIs("+-")) a2 = "(" + a2 + ")";
    return a1 + " " + name + " " + a2;
  } else if (kindIs("*/")) {
    // Multiplication and division
    Expression myarg1 = *arg1;
    if (myarg1.typeIs("double") && arg2->typeIs("Grid") && name == "/")
      myarg1 = myarg1*grid_ones;
    if (myarg1.typeIs("double") && arg2->typeIs("ReciprocalGrid") && name == "/")
      myarg1 = myarg1*recip_ones;
    if (!arg2->typeIs("double") && !myarg1.iscwise()) myarg1 = myarg1.cwise();
    std::string a1 = myarg1.printme();
    if (myarg1.kindIs("+-")) a1 = "(" + a1 + ")";
    std::string a2 = arg2->printme() ;
    if (arg2->kindIs("+-") || arg2->kindIs("*/")) a2 = "(" + a2 + ")";
    return a1 + name + a2;
  } else if (kindIs("unary") && name == "-") {
    // Unary negation
    std::string arg = arg1->printme();
    if (arg1->kindIs("+-") || arg1->kindIs("*/")) arg = "(" + arg + ")";
    return "-" + arg;
  } else if (kindIs("function") || kindIs("linear function")) {
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
  } else if (kindIs("method")) {
    std::string out = arg1->printme();
    if (arg1->kindIs("+-") || arg1->kindIs("*/") || arg1->kindIs("unary")) out = "(" + out + ")";
    out += name;
    if (arg2) {
      out += arg2->printme();
      if (arg3) out += ", " + arg3->printme();
    }
    out += ")";
    return out;
  } else if (kindIs("constant") && alias != "" && alias != "literal") {
    return alias + "xxx";
  }
  return name;
}

int Expression::checkWellFormed() const {
  int retval = 0;
  if (kindIs("+-")) {
    if ((!arg1->typeIs(type) && !arg1->typeIs("double")) ||
        (!arg2->typeIs(type) && !arg2->typeIs("double"))) {
      printf("Types don't match on addition!\n");
      retval++;
    }
  } else if (kindIs("*/")) {
    if (!arg1->typeIs(type) && !arg1->typeIs("double")) {
      printf("Types don't match on multiplication!\n");
      retval++;
    }
    if (!arg2->typeIs(type) && !arg2->typeIs("double")) {
      printf("Types don't match on multiplication!\n");
      retval++;
    }
  } else if (kindIs("unary") && name == "-") {
    if (!arg1->typeIs(type)) {
      printf("Types don't match on unary minus: %s\n", printme().c_str());
      printf("Type of %s is %s\n", arg1->printme().c_str(), arg1->type);
      printf("\n\n");
      retval++;
    }
  }
  if (arg1) retval += arg1->checkWellFormed();
  if (arg2) retval += arg2->checkWellFormed();
  if (arg3) retval += arg3->checkWellFormed();
  return retval;
}

void Expression::multiplyOut() {
  //if (!IsUnlazy()) return;
  if (kindIs("+-")) {
    arg1->multiplyOut();
    arg2->multiplyOut();
  } else {
    multiplyOutHelper();
  }
}

void Expression::multiplyOutHelper() {
  if (kindIs("*/") && name == "*") {
    arg1->multiplyOutHelper();
    arg2->multiplyOutHelper();
    if (arg1->name == "+") {
      *this = (*arg1->arg1 * *arg2) + (*arg1->arg2 * *arg2);
      multiplyOut();
    } else if (arg2->name == "+") {
      *this = (*arg2->arg1 * *arg1) + (*arg2->arg2 * *arg1);
      multiplyOut();
    }
    return;
  } else if (name == "-" && kindIs("unary")) {
    arg1->multiplyOutHelper();
    if (arg1->name == "+") {
      *this = (- *arg1->arg1) + (- *arg1->arg2);
      multiplyOut();
    }
  } else if (kindIs("linearfunexprgd")) {
    arg1->multiplyOutHelper();
    if (arg1->name == "+") {
      *this = linearfunexprgd(name.c_str(), type, *arg1->arg1)
        + linearfunexprgd(name.c_str(), type, *arg1->arg2);
      multiplyOut();
    }
  }
}

Expression Expression::simplify() const {
  Expression postfactor = *this;
  Expression prefactor = postfactor.ScalarFactor();
  if (prefactor != Expression(1)) {
    return prefactor * postfactor.simplify();
  }
  if (name == "*" && arg2->name == "*") {
    if (arg2->arg1->kindIs("constant"))
      return (*arg2->arg1 * *arg1 * *arg2->arg2).simplify();
    if (arg2->arg2->kindIs("constant"))
      return (*arg2->arg2 * *arg1 * *arg2->arg1).simplify();
    return (*arg1 * *arg2->arg1 * *arg2->arg2).simplify();
  }
  if (kindIs("unary") && name == "-") {
    Expression minusthis = arg1->simplify();
    if (minusthis.kindIs("*/") && minusthis.arg1->kindIs("constant")) {
      minusthis.arg1 = new Expression(-minusthis.arg1->value);
      return minusthis;
    }
    if (minusthis.kindIs("*/") && minusthis.arg1->kindIs("*/") && minusthis.arg1->arg1->kindIs("constant")) {
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

std::set<std::string> Expression::top_level_vars(std::set<std::string> *allvars) {
  std::set<std::string> out;
  //printf("Looking for toplevel vars in %s\n", alias.c_str());
  if (kindIs("+-") && name == "+") {
    std::set<std::string> out1 = arg1->top_level_vars(allvars);
    for (std::set<std::string>::iterator i = out1.begin(); i != out1.end(); ++i) {
      //printf("First side has %s\n", i->c_str());
      out.insert(*i);
    }
    std::set<std::string> out2 = arg2->top_level_vars(allvars);
    for (std::set<std::string>::iterator i = out2.begin(); i != out2.end(); ++i) {
      //printf("Second side has %s\n", i->c_str());
      out.insert(*i);
    }
    return out;
  }
  if (alias == "") alias = "toplevel";
  if (allvars->count(alias)) {
    // need to make this thing unique...
    for (int varnum=0; true; varnum++) {
      std::ostringstream oss;
      oss << varnum;
      std::string newa = alias + "_" + oss.str();
      if (!allvars->count(newa)) {
        alias = newa;
        break;
      }
    }
  }
  //printf("Adding %s as %s\n", printme().c_str(), alias.c_str());
  out.insert(alias);
  allvars->insert(alias);
  return out;
}

Expression Expression::FindNamedSubexpression(const std::string n) const {
  Expression e;
  if (alias == n) return *this;
  if (arg1) {
    e = arg1->FindNamedSubexpression(n);
    if (e.alias == n) return e;
  }
  if (arg2) {
    e = arg2->FindNamedSubexpression(n);
    if (e.alias == n) return e;
  }
  if (arg3) {
    e = arg3->FindNamedSubexpression(n);
    if (e.alias == n) return e;
  }
  return *this;
}

std::string Expression::FindASubexpressionName(const std::string not_this) const {
  if (arg1 && alias != "" && alias != not_this) return alias;
  if (arg1) {
    std::string n = arg1->FindASubexpressionName(not_this);
    if (n != "") return n;
  }
  if (arg2) {
    std::string n = arg2->FindASubexpressionName(not_this);
    if (n != "") return n;
  }
  if (arg3) {
    std::string n = arg3->FindASubexpressionName(not_this);
    if (n != "") return n;
  }
  return "";
}

void Expression::generate_free_code(FILE *o,  std::set<std::string> *myvars) const {
  // Free any newly unused variables...
  for (std::set<std::string>::iterator i = myvars->begin(); i != myvars->end(); ++i) {
    if (!FindVariable(*i)) {
      //fprintf(o, "    // Couldn't find %s in:  %s\n", i->c_str(), e.printme().c_str());
      //fprintf(o, "    delete %s_ptr;\n", i->c_str());
      fprintf(o, "    %s.resize(0); // We're done with this...\n", i->c_str());
      fflush(o);
      myvars->erase(i);
    }
  }
}

void Expression::generate_increment_code(FILE *o, const char *fmt, const std::string thisvar,
                                         std::set<std::string> important,
                                         std::set<std::string> *allvars,
                                         std::set<std::string> *myvars) {

  Expression e = *this;
  if (thisvar != "") e = FindNamedSubexpression(thisvar);

  if (e.kindIs("+-")) {
    arg1->generate_increment_code(o, fmt, "", important, allvars, myvars);
    arg2->generate_increment_code(o, fmt, "", important, allvars, myvars);
  } else {
    fprintf(o, "\t{\n");
    generate_code(o, fmt, thisvar, important, allvars, myvars);
    fprintf(o, "\t}\n\n");
  }
}


void Expression::generate_code(FILE *o, const char *fmt, const std::string thisvar,
                               std::set<std::string> important,
                               std::set<std::string> *allvars, std::set<std::string> *myvars) {
  std::set<std::string> otherallvars, othermyvars;
  if (!allvars) allvars = &otherallvars;
  if (!myvars) myvars = &othermyvars;
  // Set of variables that need to be deleted.
  std::set<std::string> gridvars;
  std::set<std::string> recipvars;

  Expression e = *this;
  if (thisvar != "") e = FindNamedSubexpression(thisvar);

  //printf("Starting %s...\n", thisvar.c_str());
  {
    // First, let's work out any "explicitly-named" subexpressions, common or not...
    std::string n = e.FindASubexpressionName(e.alias);
    Expression s = e.FindNamedSubexpression(n);
    //fprintf(o, "\t\t// Next to consider is %s (for %s)\n", n.c_str(), thisvar.c_str());
    while (s.alias == n && n != "" && n != e.alias) {
      {
        Expression mysube = FindNamedSubexpression(n);
        if (mysube != s) {
          printf("Nasty situation with %s.  The two options are:\n", n.c_str());
          printf("a = %s\n", s.printme().c_str());
          printf("b = %s\n", mysube.printme().c_str());
          printf("alias of b is %s, and its name is %s compared with %s\n",
                 mysube.alias.c_str(), mysube.name.c_str(), n.c_str());
          exit(1);
        }
      }
      std::string a = "my_" + n;
      if (allvars->count(a)) {
        // need to make this thing unique...
        for (int varnum=0; true; varnum++) {
          std::ostringstream oss;
          oss << varnum;
          a = n + "_" + oss.str();
          if (!allvars->count(a)) break;
        }
      }
      std::string myfmt = "\t\t" +
        (s.ctype() + (" " + a + " = %s; // explicitly-named\n"));
      //printf("I am having fun with %s", myfmt.c_str());
      //printf("I am %s\n", e.printme().c_str());
      //printf("%s is %s\n", a.c_str(), s.printme().c_str());

      generate_code(o, myfmt.c_str(), n, important, allvars, myvars);
      fflush(o);

      if (s.typeIs("Grid")) gridvars.insert(a);
      if (s.typeIs("ReciprocalGrid")) recipvars.insert(a);
      if (!s.typeIs("double")) myvars->insert(a);
      allvars->insert(a);
      s = FindNamedSubexpression(n); // Find out if the expression has since changed...
      //fprintf(o, "// expression was %s\n", printme().c_str());
      //fprintf(o, "// I eliminated %s\n", s.printme().c_str());
      EliminateThisSubexpression(s, a);
      //fprintf(o, "// expression is now %s\n", printme().c_str());
      
      // We need the following in order to update e with any changes
      // that generate_code might have made...
      if (thisvar != "") e = FindNamedSubexpression(thisvar);
      else e = *this;
      n = e.FindASubexpressionName(e.alias);
      s = e.FindNamedSubexpression(n);

      // Free unused variables...
      generate_free_code(o, myvars);
    }
  }
  //printf("All done with this one %s!\n", thisvar.c_str());

  // I don't perform the following simplification, since it makes the
  // generated code harder to read, and also takes a fair amount of
  // time to perform...
  if (false) {
    // Second, let's see what double expressions we can simplify...
    std::set<std::string> doublevars;
    Expression s = e.FindDoubleSubexpression();
    while (!s.kindIs("constant")) {
      std::string a = s.alias + "_v";
      if (a == "_v") a = "d";
      if (allvars->count(a)) {
        // need to make this thing unique...
        for (int varnum=0; true; varnum++) {
          std::ostringstream oss;
          oss << varnum;
          std::string newa = a + "_" + oss.str();
          if (!allvars->count(newa)) {
            a = newa;
            break;
          }
        }
      }
      doublevars.insert(a);
      allvars->insert(a);
      e.EliminateThisSubexpression(s, a);
      EliminateThisSubexpression(s, a);
      fprintf(o, "    const double %s = %s;\n", a.c_str(), s.printme().c_str());
      s = e.FindDoubleSubexpression();
    }
  } // Done simplifying double-valued expressions.

  Expression s = e.FindCommonSubexpression();
  //Expression ssmaller = s;
  //ssmaller.ScalarFactor(); // This is s without any prefactor.
  // Use the smaller expression if it is more common.
  //if (e.CountThisSubexpression(ssmaller) > e.CountThisSubexpression(s)) s = ssmaller;
  // Let's lump in any "cheap" calculations we can get done at the same time!
  //if (e.CountThisSubexpression(e.EasyParentOfThisSubexpression(s)) >= e.CountThisSubexpression(s))
  //  s = e.EasyParentOfThisSubexpression(s);
  int counter = 0;
  while (s != e) {
    std::string a = s.alias + "_v";
    if (a == "_v") {
      if (s.typeIs("ReciprocalGrid")) a = "recip";
      else a = "var";
    }
    if (allvars->count(a)) {
      // need to make this thing unique...
      for (int varnum=0; true; varnum++) {
        std::ostringstream oss;
        oss << varnum;
        std::string newa = a + "_" + oss.str();
        if (!allvars->count(newa)) {
          a = newa;
          break;
        }
      }
    }
    //fprintf(o, "    %s *%s_ptr = new %s(%s);\n", s.ctype(), a.c_str(), s.ctype(), s.printme().c_str());
    //fprintf(o, "    %s %s = *%s_ptr;\n", s.ctype(), a.c_str(), a.c_str());
    if (s.typeIs("Grid")) gridvars.insert(a);
    if (s.typeIs("ReciprocalGrid")) recipvars.insert(a);
    if (!s.typeIs("double")) myvars->insert(a);
    allvars->insert(a);
    e.EliminateThisSubexpression(s, a);
    EliminateThisSubexpression(s, a);
    fprintf(o, "    %s %s = %s;\n", s.ctype(), a.c_str(), s.printme().c_str());
    fprintf(o, "    //printf(\"Memory use %x is %%g with peak %%g\\n\", current_memory()/1024.0/1024, peak_memory()/1024.0/1024);\n", counter++);
    //fprintf(o, "    // expr = %s\n", e.printme().c_str());

    // Free unused variables...
    generate_free_code(o, myvars);
    for (std::set<std::string>::iterator i = gridvars.begin(); i != gridvars.end(); ++i)
      if (!FindVariable(*i)) gridvars.erase(i);
    for (std::set<std::string>::iterator i = recipvars.begin(); i != recipvars.end(); ++i)
      if (!FindVariable(*i)) recipvars.erase(i);
    
    bool have_changed = true;
    while (have_changed) {
      have_changed = false;
      // Now we look for a variable that could be simplified a bit...
      for (std::set<std::string>::iterator i = gridvars.begin(); i != gridvars.end(); ++i) {
        Expression easy = EasyParentOfThisSubexpression(Expression(*i), important);
        if (easy != Expression(*i)) {
          if (!easy.typeIs("Grid")) easy = EasyParentOfThisSubexpression(easy, important);
          if (easy.typeIs("Grid")) {
            //printf("I am reusing Grid variable %s!!!\n", i->c_str());
            e.EliminateThisSubexpression(easy, *i);
            EliminateThisSubexpression(easy, *i);
            fprintf(o, "    %s = %s; // We can reuse this variable\n", i->c_str(), easy.printme().c_str());
            fprintf(o, "    //printf(\"Memory use %x is %%g with peak %%g\\n\", current_memory()/1024.0/1024, peak_memory()/1024.0/1024);\n", counter++);
            //fprintf(o, "    // expr = %s\n", e.printme().c_str());
            fflush(o);
            have_changed = true;
          }
        }
      }
      for (std::set<std::string>::iterator i = recipvars.begin(); i != recipvars.end(); ++i) {
        Expression expi = Expression(*i).set_type("ReciprocalGrid");
        Expression easy = EasyParentOfThisSubexpression(expi, important);
        if (easy != expi) {
          if (!easy.typeIs("ReciprocalGrid")) easy = EasyParentOfThisSubexpression(easy, important);
          if (easy.typeIs("ReciprocalGrid")) {
            //printf("I am reusing ReciprocalGrid variable %s!!!\n", i->c_str());
            e.EliminateThisSubexpression(easy, *i);
            EliminateThisSubexpression(easy, *i);
            fprintf(o, "    %s = %s; // We can reuse this variable\n", i->c_str(), easy.printme().c_str());
            fprintf(o, "    //printf(\"Memory use %x is %%g with peak %%g\\n\", current_memory()/1024.0/1024, peak_memory()/1024.0/1024);\n", counter++);
            //fprintf(o, "    // expr = %s\n", e.printme().c_str());
            fflush(o);
            have_changed = true;
          }
        }
      }

      // Free any newly unused variables...
      generate_free_code(o, myvars);
      for (std::set<std::string>::iterator i = gridvars.begin(); i != gridvars.end(); ++i)
        if (!FindVariable(*i)) gridvars.erase(i);
      for (std::set<std::string>::iterator i = recipvars.begin(); i != recipvars.end(); ++i)
        if (!FindVariable(*i)) recipvars.erase(i);
    }

    // Now pick our next subexpression!
    s = e.FindCommonSubexpression();
    //Expression ssmaller = s;
    //ssmaller.ScalarFactor(); // This is s without any prefactor.
    // Use the smaller expression if it is more common.
    //if (e.CountThisSubexpression(ssmaller) > e.CountThisSubexpression(s)) s = ssmaller;
    // Let's lump in any "cheap" calculations we can get done at the same time!
    //if (e.CountThisSubexpression(e.EasyParentOfThisSubexpression(s)) >= e.CountThisSubexpression(s))
    //  s = e.EasyParentOfThisSubexpression(s);
  }
  int numvars = 0;
  for (int i=0; fmt[i]; i++) {
    if (fmt[i] == '%') numvars++;
  }
  if (numvars == 1) fprintf(o, fmt, e.printme().c_str());
  else if (numvars == 2) fprintf(o, fmt, e.printme().c_str(), e.printme().c_str());
  else {
    fprintf(o, "Looks like there is trouble here!\n");
    fprintf(o, fmt, e.printme().c_str());
  }
  fflush(o);
}
