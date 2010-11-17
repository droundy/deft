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

#include "Functionals.h"
#include "handymath.h"
#include "Grid.h"

Expression FunctionalInterface::cwiseprintme(const Expression &x) const {
  return printme(x);
}

bool FunctionalInterface::I_have_analytic_grad() const {
  return true;
}

void FunctionalInterface::pgrad(const GridDescription &gd, const VectorXd &x, const VectorXd &ingrad,
                                VectorXd *outpgrad) const {
  Grid trash(gd);
  grad(gd, x, ingrad, &trash, outpgrad);
}

void FunctionalInterface::print_summary(const char *prefix, double e, const char *name) const {
  if (name) printf("%s%25s =", prefix, name);
  else printf("%s%25s =", prefix, "UNKNOWN");
  print_double("", e);
  printf("\n");
}

double FunctionalInterface::integral(const GridDescription &gd, const VectorXd &x) const {
  return transform(gd, x).sum()*gd.dvolume;
}

// The following is a "fake" functional, used for dumping code to
// generate the gradient.
class PretendIngradType : public FunctionalInterface {
public:
  PretendIngradType() {}

  VectorXd transform(const GridDescription &, const VectorXd &x) const {
    return x;
  }
  double transform(double) const {
    return 1;
  }
  double derive(double) const {
    return 1; // hokey!
  }
  Functional grad(const Functional &, const Functional &, bool) const {
    return 0;
  }
  void grad(const GridDescription &, const VectorXd &, const VectorXd &, VectorXd *, VectorXd *) const {
  }
  Expression printme(const Expression &) const {
    return Expression("ingrad");
  }
};

/*
double average(int count, ...)
{
    va_list ap;
    int j;
    double tot = 0;
    va_start(ap, count); //Requires the last fixed parameter (to get the address)
    for(j=0; j<count; j++)
        tot+=va_arg(ap, double); //Requires the type to cast to. Increments ap to the next argument.
    va_end(ap);
    return tot/count;
}
*/
void Functional::create_source(const std::string filename, const std::string classname,
                               const char *a1, const char *a2, const char *a3,
                               const char *a4, const char *a5, const char *a6,
                               const char *a7, bool isheader) const {
  FILE *o = fopen(filename.c_str(), "w");
  if (isheader) fprintf(o, "// -*- mode: C++; -*-\n\n#pragma once\n\n");
  //printf("Generating %s\n", classname.c_str());
  fprintf(o, "#include \"Functionals.h\"\n");
  fprintf(o, "#include \"utilities.h\"\n");
  fprintf(o, "#include \"handymath.h\"\n\n");

  fprintf(o, "class %s_type : public FunctionalInterface {\n", classname.c_str());
  fprintf(o, "public:\n");
  fprintf(o, "  %s_type(", classname.c_str());
  if (a1) fprintf(o, "double %s_arg", a1);
  if (a2) fprintf(o, ", double %s_arg", a2);
  if (a3) fprintf(o, ", double %s_arg", a3);
  if (a4) fprintf(o, ", double %s_arg", a4);
  if (a5) fprintf(o, ", double %s_arg", a5);
  if (a6) fprintf(o, ", double %s_arg", a6);
  if (a7) fprintf(o, ", double %s_arg", a7);
  fprintf(o, ") ");
  if (a1) fprintf(o, ": %s(%s_arg)", a1, a1);
  if (a2) fprintf(o, ", %s(%s_arg)", a2, a2);
  if (a3) fprintf(o, ", %s(%s_arg)", a3, a3);
  if (a4) fprintf(o, ", %s(%s_arg)", a4, a4);
  if (a5) fprintf(o, ", %s(%s_arg)", a5, a5);
  if (a6) fprintf(o, ", %s(%s_arg)", a6, a6);
  if (a7) fprintf(o, ", %s(%s_arg)", a7, a7);
  fprintf(o, " {\n");
  fprintf(o, "    have_integral = true;\n");
  fprintf(o, "  }\n");
  fprintf(o, "  bool I_have_analytic_grad() const {\n");
  fprintf(o, "    return false;\n");
  fprintf(o, "  }\n");
  fprintf(o, "  double integral(const GridDescription &gd, const VectorXd &x) const {\n");
  fprintf(o, "    assert(&gd); // to avoid an unused parameter error\n");
  fprintf(o, "    assert(&x); // to avoid an unused parameter error\n");
  std::set<std::string> allvars;
  if (a1) allvars.insert(a1);
  if (a2) allvars.insert(a2);
  if (a3) allvars.insert(a3);
  if (a4) allvars.insert(a4);
  if (a5) allvars.insert(a5);
  if (a6) allvars.insert(a6);
  if (a7) allvars.insert(a7);
  Expression myself = printme(Expression("x"));
  std::set<std::string> toplevel = myself.top_level_vars(&allvars);
  {
    std::set<std::string> myvars;
    for (std::set<std::string>::iterator i = toplevel.begin(); i != toplevel.end(); ++i) {
      fprintf(o, "    // examining %s\n", i->c_str());
      Expression thisguy = myself.FindNamedSubexpression(*i);
      //fprintf(o, "    // It actually has expression %s\n", thisguy.printme().c_str());
      if (thisguy.alias != *i) {
        printf("I am erasing variable %s\n", i->c_str());
        printf("It actually has alias %s\n", thisguy.alias.c_str());
        toplevel.erase(i);
        continue;
      }
      char *buf = new char[300];
      if (thisguy.typeIs("double"))
        snprintf(buf, 300, "    %s = %%s*gd.Lat.volume();\n", i->c_str());
      else
        snprintf(buf, 300, "    %s = (%%s).sum()*gd.dvolume;\n", i->c_str());
      myself.generate_code(o, buf, *i, toplevel, &allvars, &myvars);
      //fprintf(o, "    // Myself starts as %s\n", myself.printme().c_str());
      //fprintf(o, "    // Substituting %s\n", thisguy.printme().c_str());
      myself.EliminateThisDouble(myself.FindNamedSubexpression(*i), *i);
      //fprintf(o, "    // Myself is now %s\n", myself.printme().c_str());
      delete[] buf;
      myself.generate_free_code(o, &myvars);
    }
  }
  //myself.method("sum").generate_code(o, "    return %s*gd.dvolume;\n", "", toplevel, &allvars, &myvars);
  {
    bool amfirst = true;
    fprintf(o, "    return ");
    for (std::set<std::string>::iterator i = toplevel.begin(); i != toplevel.end(); ++i) {
      if (amfirst) {
        fprintf(o, "%s", i->c_str());
        amfirst = false;
      } else {
        fprintf(o, " + %s", i->c_str());
      }
    }
    fprintf(o, ";\n");
  }
  fprintf(o, "  }\n");
  fprintf(o, "  double transform(double) const {\n");
  fprintf(o, "    assert(false);\n");
  fprintf(o, "    return 0;\n");
  fprintf(o, "  }\n");
  fprintf(o, "  double derive(double) const {\n");
  fprintf(o, "    assert(false);\n");
  fprintf(o, "    return 0;\n");
  fprintf(o, "  }\n\n");
  fprintf(o, "  VectorXd transform(const GridDescription &gd, const VectorXd &x) const {\n");
  fprintf(o, "    assert(&gd); // to avoid an unused parameter error\n");
  fprintf(o, "    assert(&x); // to avoid an unused parameter error\n");
  printme(Expression("x")).generate_code(o, "    return %s;\n");
  fprintf(o, "  }\n");
  fprintf(o, "  Functional grad(const Functional &, const Functional &, bool) const {\n");
  fprintf(o, "    assert(false);\n");
  fprintf(o, "    return 0;\n");
  fprintf(o, "  }\n\n");

  fprintf(o, "  void pgrad(const GridDescription &gd, const VectorXd &x, const VectorXd &ingrad, ");
  fprintf(o,              "VectorXd *outpgrad) const {\n");
  // The following requires that "R" always be defined!
  Functional r(1.0);
  r.set_name("R");
  Functional smoothed = StepConvolve(1.0)/Functional(4*M_PI*r*r*r);
  Expression curvature =
    grad(dV, Identity(), false).grad(dV, Identity(), false).cwiseprintme(smoothed.printme(Expression("x")));
  //Expression epg = grad(Functional(new PretendIngradType()), Identity(), true).printme(Expression("x"));
  Expression eg = grad(Functional(new PretendIngradType()), Identity(), false).printme(Expression("x"));
  Expression epg = eg * Expression("invcurvature");
  if (true || curvature.typeIs("double")) fprintf(o, "    grad(gd, x, ingrad, outpgrad, 0);\n");
  else {
    fprintf(o, "    assert(&gd); // to avoid an unused parameter error\n");
    fprintf(o, "    assert(&x); // to avoid an unused parameter error\n");
    fprintf(o, "    // curvature is %s\n", curvature.printme().c_str());
    fflush(o);
    std::set<std::string> allvars, myvars;
    (Expression(1)/curvature).generate_code(o, "    VectorXd invcurvature = %s;\n", "",
                                            std::set<std::string>(), &allvars, &myvars);
    fprintf(o, "    invcurvature = (invcurvature.cwise() > 0).select(invcurvature, 1);\n");
    fprintf(o, "    invcurvature = (invcurvature.cwise() < 1e30).select(invcurvature, 1);\n");
    //fprintf(o, "    printf(\"Our invcurvature is %%g\\n\", invcurvature.sum());\n");
    epg.generate_code(o, "    (*outpgrad) += %s;\n", "", std::set<std::string>(), &allvars, &myvars);
  }
  fprintf(o, "  }\n");

  fprintf(o, "  void grad(const GridDescription &gd, const VectorXd &x, const VectorXd &ingrad, ");
  fprintf(o,                                         "VectorXd *outgrad, VectorXd *outpgrad) const {\n");
  fprintf(o, "    assert(&gd); // to avoid an unused parameter error\n");
  fprintf(o, "    assert(&x); // to avoid an unused parameter error\n");
  fprintf(o, "    if (outpgrad) {\n");
  if (true || curvature.typeIs("double")) {
    eg.generate_code(o, "    (*outgrad) += %s;\n    (*outpgrad) += %s;\n");
  } else {
    fprintf(o, "      assert(&gd); // to avoid an unused parameter error\n");
    fprintf(o, "      assert(&x); // to avoid an unused parameter error\n");
    fprintf(o, "      // curvature is %s\n", curvature.printme().c_str());
    std::set<std::string> allvars, myvars;
    (Expression(1)/curvature).generate_code(o, "    VectorXd invcurvature = %s;\n", "",
                                            std::set<std::string>(), &allvars, &myvars);
    fprintf(o, "      invcurvature = (invcurvature.cwise() > 0).select(invcurvature, 1);\n");
    fprintf(o, "      invcurvature = (invcurvature.cwise() < 1e30).select(invcurvature, 1);\n");
    //fprintf(o, "      printf(\"Our invcurvature is %%g\\n\", invcurvature.sum());\n");
    eg.generate_code(o, "    (*outgrad) += %s;\n      (*outpgrad) += (%s).cwise() * invcurvature;\n",
                     "", std::set<std::string>(), &allvars, &myvars);
  }
  fprintf(o, "    } else {\n");
  eg = grad(Functional(new PretendIngradType()), Identity(), false).printme(Expression("x"));
  eg.generate_code(o, "      (*outgrad) += %s;\n");
  fprintf(o, "    }\n");
  fprintf(o, "  }\n");

  fprintf(o, "  Expression printme(const Expression &) const {\n");
  fprintf(o, "    return Expression(\"Can't print optimized Functionals\");\n");
  fprintf(o, "  }\n");
  fprintf(o, "  void print_summary(const char *prefix, double energy, const char *name=0) const {\n");
  if (toplevel.size() > 1) {
    for (std::set<std::string>::iterator i = toplevel.begin(); i != toplevel.end(); ++i) {
      fprintf(o, "    printf(\"%%s%25s =\", prefix);\n", i->c_str());
      fprintf(o, "    print_double(\"\", %s);\n", i->c_str());
      fprintf(o, "    printf(\"\\n\");\n");
    }
  }
  fprintf(o, "    FunctionalInterface::print_summary(prefix, energy, name);\n");
  fprintf(o, "  }\n");
  fprintf(o, "private:\n");
  if (a1) fprintf(o, "  double %s;\n", a1);
  if (a2) fprintf(o, "  double %s;\n", a2);
  if (a3) fprintf(o, "  double %s;\n", a3);
  if (a4) fprintf(o, "  double %s;\n", a4);
  if (a5) fprintf(o, "  double %s;\n", a5);
  if (a6) fprintf(o, "  double %s;\n", a6);
  if (a7) fprintf(o, "  double %s;\n", a7);
  for (std::set<std::string>::iterator i = toplevel.begin(); i != toplevel.end(); ++i) {
    fprintf(o, "  mutable double %s;\n", i->c_str());
  }
  if ((!a1 || std::string(a1) != "R") &&
      (!a2 || std::string(a2) != "R") &&
      (!a3 || std::string(a3) != "R") &&
      (!a4 || std::string(a4) != "R") &&
      (!a5 || std::string(a5) != "R") &&
      (!a6 || std::string(a6) != "R") &&
      (!a7 || std::string(a7) != "R"))
    fprintf(o, "  double R;\n");
  fprintf(o, "};\n\n");
  if (isheader) fprintf(o, "inline ");
  fprintf(o, "Functional %s(", classname.c_str());
  if (a1) fprintf(o, "double %s", a1);
  if (a2) fprintf(o, ", double %s", a2);
  if (a3) fprintf(o, ", double %s", a3);
  if (a4) fprintf(o, ", double %s", a4);
  if (a5) fprintf(o, ", double %s", a5);
  if (a6) fprintf(o, ", double %s", a6);
  if (a7) fprintf(o, ", double %s", a7);
  fprintf(o, ") {\n");
  fprintf(o, "  return Functional(new %s_type(", classname.c_str());
  if (a1) fprintf(o, "%s", a1);
  if (a2) fprintf(o, ", %s", a2);
  if (a3) fprintf(o, ", %s", a3);
  if (a4) fprintf(o, ", %s", a4);
  if (a5) fprintf(o, ", %s", a5);
  if (a6) fprintf(o, ", %s", a6);
  if (a7) fprintf(o, ", %s", a7);
  fprintf(o, "), \"%s\");\n", classname.c_str());
  fprintf(o, "}\n");
  fclose(o);
}

Expression Functional::printme(const Expression &x) const {
  Expression myself = itsCounter->ptr->printme(x);
  if (get_name() && myself.get_alias() != "literal") myself.set_alias(get_name());
  if (next()) {
    // Get associativity right...
    if (next()->next()) {
      Expression nextguy = next()->itsCounter->ptr->printme(x);
      if (next()->get_name() && nextguy.get_alias() != "literal")
        nextguy.set_alias(next()->get_name());
      return myself + nextguy + next()->next()->printme(x);
    } else {
      return myself + next()->printme(x);
    }
  } else {
    return myself;
  }
}

Expression Functional::cwiseprintme(const Expression &x) const {
  Expression myself = itsCounter->ptr->cwiseprintme(x);
  if (get_name() && myself.get_alias() != "literal") myself.set_alias(get_name());
  if (next()) {
    // Get associativity right...
    if (next()->next()) {
      Expression nextguy = next()->itsCounter->ptr->cwiseprintme(x);
      if (next()->get_name() && nextguy.get_alias() != "literal")
        nextguy.set_alias(next()->get_name());
      return myself + nextguy + next()->next()->cwiseprintme(x);
    } else {
      return myself + next()->cwiseprintme(x);
    }
  } else {
    return myself;
  }
}

int Functional::run_finite_difference_test(const char *testname, const Grid &x,
                                                 const VectorXd *direction) const {
  printf("\nRunning finite difference test on %s:\n", testname);
  const GridDescription gd(x.description());

  int retval = 0;

  VectorXd my_grad(x);
  double Eold = integral(x);
  my_grad.setZero();
  integralgrad(x, &my_grad);
  VectorXd my_direction(my_grad);
  if (direction) my_direction = *direction;
  my_direction /= my_direction.norm();

  const double lderiv = my_direction.dot(my_grad);
  if (lderiv == 0.0) {
    printf("FAIL: Gradient is zero...\n");
    retval++;
  }
  {
    // compare the gradient computed by the two functions grad() and
    // grad_and_pgrad()
    VectorXd my_pgrad(x);
    my_grad.setZero();
    my_pgrad.setZero();
    integralgrad(x, &my_grad, &my_pgrad);
    const double lderiv_new = my_direction.dot(my_grad);
    if (fabs(lderiv_new/lderiv - 1) > 1e-12) {
      printf("\n*** WARNING!!! INCONSISTENT GRADIENTS! ***\n");
      printf("Different gradient in energy_and_grad() and grad_and_pgrad()\n");
      printf("FAIL: Fractional error is %g\n\n", lderiv_new/lderiv - 1);
      retval++;
    }

    if (I_have_analytic_grad()) {
      // We have an analytic grad, so let's make sure it matches the
      // other one...
      VectorXd othergrad(grad(dV, Identity(), false)(x));
      double maxgraderr = (othergrad-my_grad).cwise().abs().maxCoeff();
      double maxgrad = my_grad.cwise().abs().maxCoeff();
      if (maxgraderr/maxgrad > 1e-12) {
        printf("func itself is %s\n", printme(Expression("x")).printme().c_str());
        printf("othergrad[0] = %g\n", othergrad[0]);
        printf("othergrad itself is %s\n",
               grad(dV, Identity(), false).printme(Expression("x")).printme().c_str());
        printf("my_grad[0] = %g\n", my_grad[0]);
        printf("maxgraderr is %g while maxgrad is %g\n", maxgraderr, maxgrad);
        printf("FAIL: Discrepancy in the gradient is too big: %g\n\n", maxgraderr/maxgrad);      
        retval++;
      }
    }
  }

  // We try to choose as small values for epsilon as are consistent with
  // getting meaningful results:
  const double sigfigs_best = 1e-15;
  const double epsilon_max = 1e-15*fabs(Eold)/fabs(lderiv)/sigfigs_best;
  const int min_p = (int) -floor(log10(epsilon_max) + 0.5);
  const int max_p = min_p + 7;
  printf("choosing delta of 1e%d based on directional grad %g and energy %g\n",
         -min_p, lderiv, Eold);
  double min = 1e300, best_ratio_error = 1.0;
  VectorXd grads(max_p + 1 - min_p);
  VectorXd diff_grads(max_p - min_p);
  // Take steps in the grad direction.
  for(int p=min_p; p <= max_p; p++) {
    const double eps_ratio = 10.0;
    const double epsilon = pow(eps_ratio, -p);
    // The following is a little wasteful of memory...
    const double Eplus= integral(gd, x + epsilon*my_direction);
    const double Eminus=integral(gd, x - epsilon*my_direction);

    grads[p-min_p] = (Eplus-Eminus)/(2*epsilon);
    printf("    eps^2 = %25.16f deltaE %.12g\n",
           epsilon*epsilon*pow(eps_ratio, 2.0*min_p), Eplus-Eminus);
    printf("FD   Ratio: %25.16f (grad ~ %g)\n", grads[p-min_p]/lderiv,
           grads[p-min_p]);
    printf("FD sigfigs: %25.16f\n", 1e-15*fabs(Eold/(Eplus-Eminus)));
    fflush(stdout);

    if (fabs(grads[p-min_p]/lderiv - 1.0) < best_ratio_error)
      best_ratio_error = fabs(grads[p-min_p]/lderiv - 1.0);
    if (p > min_p) {
      diff_grads[p-min_p-1]=(grads[p-min_p]-grads[p-min_p-1]);
      double diff = (grads[p-min_p-1]-lderiv)*
        (-1+1./(eps_ratio*eps_ratio))/diff_grads[p-min_p-1];
      if ( min > fabs(diff-1) ) min = fabs(diff-1);
      //dft_log("FD diff: %25.16lf\n\n", diff);
    }
  }

  if (min < 1e-3 && best_ratio_error < 1e-7) {
    printf("Passed on basis of reasonable scaling (%g) and accuracy (%g).\n",
           min, best_ratio_error);
  } else if (best_ratio_error < 1e-10) {
    printf("Passed on basis of a gradient accuracy of (%g) (with scaling %g).\n",
           best_ratio_error, min);
  } else if (min < 1e-5 && best_ratio_error < 0.01) {
    printf("Passed on basis of seriously nice scaling %g (and low accuracy %g).\n",
           min, best_ratio_error);
  } else {
    printf("FAIL: Failed with scaling ratio of %g and gradient accuracy of %g\n",
           min, best_ratio_error);
    retval++;
  }
  return retval;
}

double Functional::print_iteration(const char *prefix, int iter) const {
  printf("%s==============\n", prefix);
  printf("%sIteration %4d\n", prefix, iter);
  printf("%s==============\n", prefix);

  print_summary(prefix, itsCounter->last_energy, get_name());
  const Functional *nxt = this;
  double etot = 0;
  while (nxt) {
    etot += nxt->itsCounter->last_energy;
    //printf("In print_iteration summarizing %p\n", nxt);
    //nxt->print_summary(prefix, nxt->itsCounter->last_energy);
    nxt = nxt->next();
  }
  printf("%s%25s =", prefix, "total energy");
  print_double("", etot);
  printf("\n");
  return etot;
}

void Functional::print_summary(const char *prefix, double energy, const char *name) const {
  if (get_name()) name = get_name();
  if (!next()) {
    itsCounter->ptr->print_summary(prefix, energy, name);
    return;
  }
  const Functional *nxt = this;
  double etot = 0;
  while (nxt) {
    if (nxt->get_name()) name = nxt->get_name();
    etot += nxt->itsCounter->last_energy;
    nxt->itsCounter->ptr->print_summary(prefix, nxt->itsCounter->last_energy, name);
    nxt = nxt->next();
  }
  //assert(fabs(etot - energy) < 1e-6);
}

Functional Identity() { return Pow(1); }

class dVType : public FunctionalInterface {
public:
  dVType() {}

  VectorXd transform(const GridDescription &gd, const VectorXd &) const {
    return gd.dvolume*VectorXd::Ones(gd.NxNyNz);
  }
  double transform(double) const {
    return 1;
  }
  double derive(double) const {
    return 1; // hokey!
  }
  Functional grad(const Functional &, const Functional &, bool) const {
    return 0;
  }
  void grad(const GridDescription &, const VectorXd &, const VectorXd &, VectorXd *, VectorXd *) const {
  }
  Expression printme(const Expression &) const {
    return Expression("gd.dvolume").set_type("double");
  }
};

Functional dV = Functional(new dVType());

class Constant : public FunctionalInterface {
public:
  Constant(double x, const char *n) : c(x), name(n) {}

  VectorXd transform(const GridDescription &, const VectorXd &data) const {
    return c*VectorXd::Ones(data.rows());
  }
  double transform(double) const {
    return c;
  }
  double derive(double) const {
    return 0;
  }
  Functional grad(const Functional &, const Functional &, bool) const {
    return 0;
  }
  void grad(const GridDescription &, const VectorXd &, const VectorXd &, VectorXd *, VectorXd *) const {
  }
  Expression printme(const Expression &) const {
    if (name) return Expression(name).set_type("double");
    return c;
  }
private:
  double c;
  const char *name;
};

Functional::Functional(double x, const char *n) : itsCounter(0) {
  init(new Constant(x, n), 0);
}

class ConstantField : public FunctionalInterface {
public:
  ConstantField(const VectorXd &x) : c(x) {}

  VectorXd transform(const GridDescription &, const VectorXd &) const {
    return c;
  }
  double transform(double) const {
    return c.sum()/c.rows();
  }
  double derive(double) const {
    return 0;
  }
  Functional grad(const Functional &, const Functional &, bool) const {
    return 0;
  }
  void grad(const GridDescription &, const VectorXd &, const VectorXd &, VectorXd *, VectorXd *) const {
  }
  Expression printme(const Expression &) const {
    return Expression("an unknown field...");
  }
private:
  VectorXd c;
};

Functional::Functional(const VectorXd &x) : itsCounter(0) {
  init(new ConstantField(x), 0);
}

class ChainRuleType : public FunctionalInterface {
public:
  ChainRuleType(const Functional &fa, const Functional &fb) : f1(fa), f2(fb) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    if (!f1.next()) {
      // This is the simple, efficient case!
      return f1(gd, f2(gd, data));
    }
    // This does some extra work to save the energies of each term in
    // the sum, just in case we want to print it!
    VectorXd f2data(f2(gd,data));
    VectorXd f1f2data(f1.justMe(gd, f2data));
    double e = gd.dvolume*f1f2data.sum();
    f1.set_last_energy(e);
    Functional *nxt = f1.next();
    while (nxt) {
      f1f2data += nxt->justMe(gd, f2data);
      double etot = gd.dvolume*f1f2data.sum();
      nxt->set_last_energy(etot - e);
      e = etot;
      nxt = nxt->next();
    }
    return f1f2data;
  }
  double transform(double n) const {
    return f1(f2(n));
  }
  double derive(double n) const {
    return f1.derive(f2(n))*f2.derive(n);
  }
  Functional grad(const Functional &ingrad, const Functional &x, bool ispgrad) const {
    return f2.grad(f1.grad(ingrad, f2(x), ispgrad), x, ispgrad);
    return f1.grad(f2.grad(ingrad, x, ispgrad), f2, ispgrad);
  }
  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    Grid outgrad1(gd);
    outgrad1.setZero();
    f1.grad(gd, f2(gd, data), ingrad, &outgrad1, 0);
    f2.grad(gd, data, outgrad1, outgrad, outpgrad);
  }
  void print_summary(const char *prefix, double e, const char *name) const {
    f1.print_summary(prefix, e, name);
  }
  Expression printme(const Expression &x) const {
    return f1.printme(f2.printme(x));
  }
  Expression cwiseprintme(const Expression &x) const {
    return f1.cwiseprintme(f2.cwiseprintme(x));
  }
  bool I_have_analytic_grad() const {
    return f1.I_have_analytic_grad() && f2.I_have_analytic_grad();
  }
private:
  Functional f1, f2;
};

Functional Functional::operator()(const Functional &f) const {
  return Functional(new ChainRuleType(*this, f));
}

class QuotientRuleType : public FunctionalInterface {
public:
  QuotientRuleType(const Functional &fa, const Functional &fb) : f1(fa), f2(fb) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    return f1(gd, data).cwise()/f2(gd, data);
  }
  double transform(double n) const {
    return f1(n)/f2(n);
  }
  double derive(double n) const {
    return f1.derive(n)/f2(n) - f1(n)*f2.derive(n)/f2(n)/f2(n);
  }
  Functional grad(const Functional &ingrad, const Functional &x, bool ispgrad) const {
    return f1.grad(ingrad/f2(x), x, ispgrad) - f2.grad(f1(x)*ingrad/sqr(f2(x)), x, ispgrad);
  }
  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    VectorXd out2 = f2(gd, data);
    f1.grad(gd, data, ingrad.cwise()/out2, outgrad, outpgrad);
    f2.grad(gd, data, (ingrad.cwise()*f1(gd, data)).cwise()/((-out2).cwise()*out2), outgrad, outpgrad);
  }
  Expression printme(const Expression &x) const {
    return f1.printme(x) / f2.printme(x);
  }
  Expression cwiseprintme(const Expression &x) const {
    return f1.cwiseprintme(x) / f2.cwiseprintme(x);
  }
  bool I_have_analytic_grad() const {
    return f1.I_have_analytic_grad() && f2.I_have_analytic_grad();
  }
private:
  Functional f1, f2;
};

Functional Functional::operator/(const Functional &f) const {
  return Functional(new QuotientRuleType(*this, f));
}

class ProductRuleType : public FunctionalInterface {
public:
  ProductRuleType(const Functional &fa, const Functional &fb) : f1(fa), f2(fb) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    return f1(gd, data).cwise()*f2(gd, data);
  }
  double transform(double n) const {
    return f1(n)*f2(n);
  }
  double derive(double n) const {
    return f1(n)*f2.derive(n) + f1.derive(n)*f2(n);
  }
  Functional grad(const Functional &ingrad, const Functional &x, bool ispgrad) const {
    return f2.grad(f1(x)*ingrad, x, ispgrad) + f1.grad(f2(x)*ingrad, x, ispgrad);
  }
  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    f1.grad(gd, data, ingrad.cwise()*f2(gd, data), outgrad, outpgrad);
    f2.grad(gd, data, ingrad.cwise()*f1(gd, data), outgrad, outpgrad);
  }
  Expression printme(const Expression &x) const {
    return f1.printme(x)*f2.printme(x);
  }
  Expression cwiseprintme(const Expression &x) const {
    return f1.cwiseprintme(x)*f2.cwiseprintme(x);
  }
  bool I_have_analytic_grad() const {
    return f1.I_have_analytic_grad() && f2.I_have_analytic_grad();
  }
private:
  Functional f1, f2;
};

Functional Functional::operator*(const Functional &f) const {
  return Functional(new ProductRuleType(*this, f));
}


class LogType : public FunctionalInterface {
public:
  LogType() {}

  VectorXd transform(const GridDescription &, const VectorXd &data) const {
    return data.cwise().log();
  }
  double transform(double n) const {
    return log(n);
  }
  double derive(double n) const {
    return 1/n;
  }
  Functional grad(const Functional &ingrad, const Functional &x, bool) const {
    return ingrad/x;
  }
  void grad(const GridDescription &, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    *outgrad += ingrad.cwise()/data;
    if (outpgrad) *outpgrad += ingrad.cwise()/data;
  }
  Expression printme(const Expression &x) const {
    return log(x);
  }
};

Functional log(const Functional &f) {
  return Functional(new LogType())(f);
}


class ExpType : public FunctionalInterface {
public:
  ExpType() {}

  VectorXd transform(const GridDescription &, const VectorXd &data) const {
    return data.cwise().exp();
  }
  double transform(double n) const {
    return exp(n);
  }
  double derive(double n) const {
    return exp(n);
  }
  Functional grad(const Functional &ingrad, const Functional &x, bool) const {
    return ingrad*Functional(new ExpType())(x);
  }
  void grad(const GridDescription &, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    *outgrad += ingrad.cwise() * data.cwise().exp();
    if (outpgrad) *outpgrad += ingrad.cwise() * data.cwise().exp();
  }
  Expression printme(const Expression &x) const {
    return exp(x);
  }
};

Functional exp(const Functional &f) {
  return Functional(new ExpType())(f);
}


Functional Functional::operator-(const Functional &f) const {
  return *this + -f;
}

Functional Functional::operator-() const {
  return -1* *this;
}

Functional sqr(const Functional &f) {
  return Pow(2)(f);
}


class Constraint : public FunctionalInterface {
public:
  Constraint(const Grid &g, const Functional &y) : constraint(g), f(y) {};

  double integral(const GridDescription &gd, const VectorXd &x) const {
    return f.integral(gd, x);
  }
  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    return f(gd, data);
  }
  double transform(double n) const {
    return f(n);
  }
  double derive(double n) const {
    return f.derive(n);
  }
  Functional grad(const Functional &ingrad, const Functional &x, bool ispgrad) const {
    return Functional(constraint)*f.grad(ingrad, x, ispgrad);
  }
  void pgrad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
             VectorXd *outpgrad) const {
    VectorXd mypgrad(data);
    mypgrad.setZero();
    f.pgrad(gd, data, ingrad, &mypgrad);
    *outpgrad += constraint.cwise() * mypgrad;
  }
  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    if (outpgrad) {
      VectorXd mygrad(data), mypgrad(data);
      mygrad.setZero();
      mypgrad.setZero();
      f.grad(gd, data, ingrad, &mygrad, &mypgrad);
      *outgrad += constraint.cwise() * mygrad;
      *outpgrad += constraint.cwise() * mypgrad;
    } else {
      VectorXd mygrad(data);
      mygrad.setZero();
      f.grad(gd, data, ingrad, &mygrad, 0);
      *outgrad += constraint.cwise() * mygrad;
    }
  }
  Expression printme(const Expression &x) const {
    return f.printme(x);
  }
  bool I_have_analytic_grad() const {
    return f.I_have_analytic_grad();
  }
private:
  const Grid constraint;
  const Functional f;
};

Functional constrain(const Grid &g, Functional f) {
  return Functional(new Constraint(g, f));
}
