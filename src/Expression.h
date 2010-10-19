// -*- mode: C++; -*-

#pragma once

#include <string>
#include <set>

class Expression {
public:
  Expression();
  Expression(double); // Constant literal
  Expression(std::string); // Variable
  Expression(const Expression &);
  void operator=(const Expression &);
  ~Expression() {
    delete arg1;
    delete arg2;
  }

  std::string name, alias, kind, type;
  Expression *arg1, *arg2, *arg3;
  double value;
  bool unlazy;
  Expression set_type(const char *t) {
    type = t;
    return *this;
  }
  const char *ctype() const {
    if (type == "double") return "double";
    if (type == "ReciprocalGrid") return "VectorXcd";
    return "VectorXd";
  }

  Expression operator+(const Expression &) const;
  Expression operator-() const;
  Expression operator-(const Expression &) const;
  Expression operator*(const Expression &) const;
  Expression operator/(const Expression &) const;
  Expression operator()(const Expression &) const;
  Expression operator()(const Expression &, const Expression &) const;
  Expression cwise() const; // creates a coefficient-wise version...
  bool iscwise() const; // true if we're already coefficient-wise
  Expression method(const char *) const;
  Expression method(const char *, const Expression &) const;
  Expression method(const char *, const Expression &, const Expression &) const;

  Expression set_alias(std::string a);
  std::string get_alias() const {
    return alias;
  }

  // Modifies this expression in-place, and returns a declaration of
  // the common subexpression.
  std::string EliminateCommonSubexpression();
  bool EliminateThisSubexpression(const Expression &, const std::string alias);
  bool EliminateThisDouble(const Expression &, const std::string alias);
  Expression EasyParentOfThisSubexpression(const Expression &, std::set<std::string> important) const;
  bool FindVariable(const std::string n) const;
  Expression FindNamedSubexpression(const std::string n) const;

  Expression FindCommonSubexpression() const;
  int CountThisSubexpression(const Expression &) const;
  int Depth() const;
  bool operator==(const Expression &) const;
  bool operator!=(const Expression &e) const {
    return !(*this == e);
  }
  bool IsUnlazy() const;
  bool IsUnlazyApartFrom(const Expression &, std::set<std::string> important) const;

  Expression ScalarFactor(); // Removes any scalar factors from *this

  std::string printme() const;
  int checkWellFormed() const;

  Expression simplify() const;
  void generate_code(FILE *outfile, const char *fmt, const std::string thisvar = "",
                     std::set<std::string> important = std::set<std::string>(),
                     std::set<std::string> *scopevars = 0, std::set<std::string> *myvars = 0);
  void generate_free_code(FILE *o,  std::set<std::string> *myvars) const;
  std::set<std::string> top_level_vars(std::set<std::string> *allvars);
};

Expression funexpr(const char *name);
Expression funexpr(const char *name, const Expression &);
Expression funexpr(const char *name, const Expression &, const Expression &);
Expression funexpr(const char *name, const Expression &, const Expression &, const Expression &);

Expression linearfunexprgd(const char *name, const char *type, const Expression &);

inline Expression operator+(double a, const Expression &b) {
  return Expression(a) + b;
}

inline Expression operator-(double a, const Expression &b) {
  return Expression(a) - b;
}

inline Expression operator*(double a, const Expression &b) {
  return Expression(a) * b;
}

inline Expression operator/(double a, const Expression &b) {
  return Expression(a) / b;
}

Expression fft(const Expression &x);
Expression ifft(const Expression &x);

Expression abs(const Expression &);
Expression log(const Expression &);
Expression exp(const Expression &);
Expression sqr(const Expression &);
