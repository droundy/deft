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
  Expression thispost(*this), epost(e);
  Expression thispre = thispost.ScalarFactor(), epre = epost.ScalarFactor();
  if (thispost.kind == "linear function" && epost.kind == "linear function" &&
      thispost.name == epost.name) {
    if (thispre == epre) {
      Expression out = linearfunexprgd("oops", thispost.type.c_str(), *thispost.arg1 + *epost.arg1);
      out.name = thispost.name;
      return epre*out;
    }
    if (thispre == -epre) {
      Expression out = linearfunexprgd("oops", thispost.type.c_str(), *thispost.arg1 - *epost.arg1);
      out.name = thispost.name;
      return thispre*out;
    }
    Expression out = linearfunexprgd("oops", thispost.type.c_str(),
                                     thispre * *thispost.arg1 + epre * *epost.arg1);
    out.name = thispost.name;
    return out;
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

  Expression thispost(*this), epost(e);
  Expression thispre = thispost.ScalarFactor(), epre = epost.ScalarFactor();
  // for some reason the following makes memory use worse!!!
  if (thispost.kind == "linear function" && epost.kind == "linear function" &&
      thispost.name == epost.name) {
    if (thispre == epre) {
      Expression out = linearfunexprgd("oops", thispost.type.c_str(), *thispost.arg1 - *epost.arg1);
      out.name = thispost.name;
      return epre*out;
    }
    if (thispre == -epre) {
      Expression out = linearfunexprgd("oops", thispost.type.c_str(), *thispost.arg1 + *epost.arg1);
      out.name = thispost.name;
      return thispre*out;
    }
    Expression out = linearfunexprgd("oops", thispost.type.c_str(),
                                     thispre * *thispost.arg1 - epre * *epost.arg1);
    out.name = thispost.name;
    return out;
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
  } else if (e.kind == "unary" && e.name == "-") {
    return - (*this * *e.arg1);
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

Expression Expression::ScalarFactor() {
  if (kind == "constant") {
    Expression out(value);
    value = 1;
    return out;
  }
  if (type == "double") {
    // If we have type "double" then we *are* a scalar!
    Expression out(*this);
    *this = Expression(1);
    return out;
  }
  if (kind == "*/" && name == "*") {
    Expression sf1 = arg1->ScalarFactor();
    Expression sf2 = arg2->ScalarFactor();
    *this = *arg1 * *arg2; // to handle cases where arg1 or arg2 becomes 1.
    return sf1*sf2;
  }
  if (kind == "*/" && name == "/") {
    Expression sf1 = arg1->ScalarFactor();
    Expression sf2 = arg2->ScalarFactor();
    *this = *arg1 / *arg2; // to handle cases where arg1 or arg2 becomes 1.
    return sf1/sf2;
  }
  if (kind == "linear function") {
    return arg1->ScalarFactor();
  }
  if (kind == "unary" && name == "-") {
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

Expression linearfunexprgd(const char *n, const char *type, const Expression &arg) {
  if (arg.kind == "constant" && arg.value == 0)
    return Expression(0); // Nice optimization!  :)
  Expression newarg(arg);
  Expression prefactor = newarg.ScalarFactor();

  Expression out = funexpr(n, newarg);
  out.kind = "linear function";
  out.name = std::string(n) + "(gd, ";
  out.unlazy = true;
  out.type = type;
  return prefactor*out;
}

Expression fft(const Expression &g) {
  if (g.type != "Grid") {
    if (g.type == "double") {
      // The fft of a constant is a delta function in G space!
      printf("fft: Expression should have type Grid, but accepting double anyhow.\n");
      // FIXME: W shouldn't need to do an actual FFT here!
      return fft(g * grid_ones);
    }
    printf("fft: Expression %s should have type Grid.\n", g.printme().c_str());
    exit(1);
  }
  return linearfunexprgd("fft", "ReciprocalGrid", g);
}

Expression ifft(const Expression &g) {
  if (g.type != "ReciprocalGrid") {
    printf("ifft: Expression '%s' should have type ReciprocalGrid but instead has type '%s'.\n",
           g.printme().c_str(), g.type.c_str());
    exit(1);
  }
  return linearfunexprgd("ifft", "Grid", g);
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
  Expression postfactor = *this;
  Expression prefactor = postfactor.ScalarFactor();
  if (prefactor != Expression(1)) {
    return prefactor * postfactor.simplify();
  }
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
  std::set<std::string> gridvars;
  std::set<std::string> recipvars;

  Expression e = *this;
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
      if (s.type == "ReciprocalGrid") a = "recip";
      else a = "var";
    }
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
    //fprintf(o, "    %s *%s_ptr = new %s(%s);\n", s.ctype(), a.c_str(), s.ctype(), s.printme().c_str());
    //fprintf(o, "    %s %s = *%s_ptr;\n", s.ctype(), a.c_str(), a.c_str());
    if (s.type == "Grid") gridvars.insert(a);
    if (s.type == "ReciprocalGrid") recipvars.insert(a);
    allvars.insert(a);
    while (e.EliminateThisSubexpression(s, a));
    fprintf(o, "    %s %s = %s;\n", s.ctype(), a.c_str(), s.printme().c_str());
    fprintf(o, "    printf(\"Memory use %x is %%g with peak %%g\\n\", current_memory()/1024.0/1024, peak_memory()/1024.0/1024);\n", counter++);
    //fprintf(o, "    // expr = %s\n", e.printme().c_str());
    fflush(o);

    // Free unused variables...
    for (std::set<std::string>::iterator i = gridvars.begin(); i != gridvars.end(); ++i) {
      if (!e.FindVariable(*i)) {
        //fprintf(o, "    // Couldn't find %s in:  %s\n", i->c_str(), e.printme().c_str());
        //fprintf(o, "    delete %s_ptr;\n", i->c_str());
        fprintf(o, "    %s.resize(0); // We're done with this...\n", i->c_str());
        fflush(o);
        gridvars.erase(i);
      }
    }
    for (std::set<std::string>::iterator i = recipvars.begin(); i != recipvars.end(); ++i) {
      if (!e.FindVariable(*i)) {
        //fprintf(o, "    // Couldn't find %s in:  %s\n", i->c_str(), e.printme().c_str());
        //fprintf(o, "    delete %s_ptr;\n", i->c_str());
        fprintf(o, "    %s.resize(0); // We're done with this...\n", i->c_str());
        fflush(o);
        recipvars.erase(i);
      }
    }
    
    bool have_changed = true;
    while (have_changed) {
      have_changed = false;
      // Now we look for a variable that could be simplified a bit...
      for (std::set<std::string>::iterator i = gridvars.begin(); i != gridvars.end(); ++i) {
        Expression easy = e.EasyParentOfThisSubexpression(Expression(*i));
        if (easy != Expression(*i)) {
          if (easy.type != "Grid") easy = e.EasyParentOfThisSubexpression(easy);
          if (easy.type == "Grid") {
            //printf("I am reusing Grid variable %s!!!\n", i->c_str());
            while (e.EliminateThisSubexpression(easy, "THIS IS A UNIQUE NAME"));
            while (e.EliminateThisSubexpression(Expression("THIS IS A UNIQUE NAME"), *i));
            fprintf(o, "    %s = %s; // We can reuse this variable\n", i->c_str(), easy.printme().c_str());
            fprintf(o, "    printf(\"Memory use %x is %%g with peak %%g\\n\", current_memory()/1024.0/1024, peak_memory()/1024.0/1024);\n", counter++);
            //fprintf(o, "    // expr = %s\n", e.printme().c_str());
            fflush(o);
            have_changed = true;
          }
        }
      }
      for (std::set<std::string>::iterator i = recipvars.begin(); i != recipvars.end(); ++i) {
        Expression expi = Expression(*i).set_type("ReciprocalGrid");
        Expression easy = e.EasyParentOfThisSubexpression(expi);
        if (easy != expi) {
          if (easy.type != "ReciprocalGrid") easy = e.EasyParentOfThisSubexpression(easy);
          if (easy.type == "ReciprocalGrid") {
            //printf("I am reusing ReciprocalGrid variable %s!!!\n", i->c_str());
            while (e.EliminateThisSubexpression(easy, "THIS IS A UNIQUE NAME"));
            while (e.EliminateThisSubexpression(Expression("THIS IS A UNIQUE NAME").set_type("ReciprocalGrid"), *i));
            fprintf(o, "    %s = %s; // We can reuse this variable\n", i->c_str(), easy.printme().c_str());
            fprintf(o, "    printf(\"Memory use %x is %%g with peak %%g\\n\", current_memory()/1024.0/1024, peak_memory()/1024.0/1024);\n", counter++);
            //fprintf(o, "    // expr = %s\n", e.printme().c_str());
            fflush(o);
            have_changed = true;
          }
        }
      }

      // Free any newly unused variables...
      for (std::set<std::string>::iterator i = gridvars.begin(); i != gridvars.end(); ++i) {
        if (!e.FindVariable(*i)) {
          //fprintf(o, "    // Couldn't find %s in:  %s\n", i->c_str(), e.printme().c_str());
          //fprintf(o, "    delete %s_ptr;\n", i->c_str());
          fprintf(o, "    %s.resize(0); // We're done with this...\n", i->c_str());
          fflush(o);
          gridvars.erase(i);
        }
      }
      for (std::set<std::string>::iterator i = recipvars.begin(); i != recipvars.end(); ++i) {
        if (!e.FindVariable(*i)) {
          //fprintf(o, "    // Couldn't find %s in:  %s\n", i->c_str(), e.printme().c_str());
          //fprintf(o, "    delete %s_ptr;\n", i->c_str());
          fprintf(o, "    %s.resize(0); // We're done with this...\n", i->c_str());
          fflush(o);
          recipvars.erase(i);
        }
      }
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
  fprintf(o, fmt, e.printme().c_str());
  fflush(o);
}
