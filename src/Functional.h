// -*- mode: C++; -*-

#pragma once

#include "ReciprocalGrid.h"
#include "Expression.h"

class Functional;

class FunctionalInterface {
public:
  FunctionalInterface() : have_integral(false) {}
  virtual ~FunctionalInterface() {}

  // A functional mapping one field onto another...
  virtual double integral(const GridDescription &gd, double kT, const VectorXd &data) const;
  virtual VectorXd transform(const GridDescription &gd, double kT, const VectorXd &data) const = 0;
  virtual double transform(double kT, double x) const = 0;

  virtual void pgrad(const GridDescription &gd, double kT, const VectorXd &x,
                     const VectorXd &ingrad, VectorXd *outpgrad) const;

  // This computes the gradient of the functional, given a gradient of
  // its output field (i.e. it applies the chain rule).
  virtual void grad(const GridDescription &gd, double kT, const VectorXd &data,
                    const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const = 0;
  virtual double derive(double kT, double data) const = 0;
  virtual Expression derive_homogeneous(const Expression &x) const = 0;
  virtual double d_by_dT(double kT, double data) const = 0;
  virtual Functional grad(const Functional &ingrad, const Functional &x, bool ispgrad) const = 0;
  virtual Functional grad_T(const Functional &ingradT) const = 0;
  virtual Expression printme(const Expression &) const = 0;

  virtual void print_summary(const char *prefix, double energy, std::string name) const;
  virtual bool I_have_analytic_grad() const;
  virtual bool I_am_constant_wrt_x() const;
  virtual bool I_preserve_homogeneous() const;
  virtual bool I_am_homogeneous() const;
  virtual bool I_am_zero() const;
  virtual bool I_am_one() const;
  virtual bool I_give_zero_for_zero() const;
  virtual bool I_am_local() const; // WARNING:  this defaults to true!

  bool have_integral;
};

template<typename Derived, typename extra> class ConvolveWith;

class Functional {
public:
  // Handle reference counting so we can pass these things around freely...
  explicit Functional(double, const char *name=0); // This handles constants!
  explicit Functional(const VectorXd &); // This handles constant fields!
  template<typename Derived, typename extra>
  explicit Functional(Derived (*f)(const GridDescription &, extra), extra e,
                      const Expression &R, const Expression &gzero, bool iseven)
    : itsCounter(0) {
    Lattice lat(Cartesian(1,0,0), Cartesian(0,1,0), Cartesian(0,0,1));
    GridDescription gd(lat, 2, 2, 2);
    // This handles constant ephemeral fields!
    std::string n = f(gd, e).name();
    if (R != Expression("R")) n = "";
    init(new ConvolveWith<Derived,extra>(f,e,R,gzero,iseven), n);
  }
  explicit Functional(FunctionalInterface* p = 0, const std::string name = "") // allocate a new counter
    : itsCounter(0) {
    init(p, name);
  }
  ~Functional() {
    release();
    delete mynext;
  }
  Functional(const Functional& r) {
    acquire(r.itsCounter);
    if (r.mynext) mynext = new Functional(*r.mynext);
    else mynext = 0;
  }
  Functional& operator=(const Functional& r) {
    if (this != &r) {
      release();
      acquire(r.itsCounter);
    }
    if (r.mynext) mynext = new Functional(*r.mynext);
    else {
      if (mynext) mynext->release();
      mynext = 0;
    }
    return *this;
  }

  Functional operator()(const Functional &) const;
  Functional operator+(const Functional &x) const {
    Functional out(*this);
    return out += x;
  }
  Functional operator+=(const Functional &x) {
    if (I_am_zero()) {
      *this = x;
      return *this;
    }
    if (x.I_am_zero()) {
      return *this;
    }
    if (mynext) *mynext += x;
    else mynext = new Functional(x);
    return *this;
  }
  Functional operator-() const;
  Functional operator-(const Functional &) const;
  Functional operator/(const Functional &) const;
  Functional operator*(const Functional &) const;
  VectorXd justMe(const GridDescription &gd, double kT, const VectorXd &data) const {
    return itsCounter->ptr->transform(gd, kT, data);
  }
  double justMeIntegral(const GridDescription &gd, double kT, const VectorXd &data) const {
    return itsCounter->ptr->integral(gd, kT, data);
  }
  VectorXd operator()(double kT, const GridDescription &gd, const VectorXd &data) const {
    return (*this)(gd, kT, data);
  }
  VectorXd operator()(const GridDescription &gd, double kT, const VectorXd &data) const {
    VectorXd out = itsCounter->ptr->transform(gd, kT, data);
    if (mynext) out += (*mynext)(gd, kT, data);
    return out;
  }
  VectorXd operator()(double kT, const Grid &g) const {
    return (*this)(kT, g.description(), g);
  }
  double integral(double kT, const Grid &g) const {
    return integral(g.description(), kT, g);
  }
  double integral(double kT, const GridDescription &gd, const VectorXd &data) const {
    return integral(gd, kT, data);
  }
  double integral(const GridDescription &gd, double kT, const VectorXd &data) const {
    // This takes care to save the energies of each term in the sum.
    double e = itsCounter->ptr->integral(gd, kT, data);
    set_last_energy(e);
    Functional *nxt = next();
    while (nxt) {
      double enext = nxt->itsCounter->ptr->integral(gd, kT, data);
      nxt->set_last_energy(enext);
      e += enext;
      nxt = nxt->next();
    }
    return e;
  }
  void integralgrad(double kT, const Grid &g, VectorXd *gr, VectorXd *pg=0) const {
    grad(kT, g.description(), g, g.description().dvolume*VectorXd::Ones(g.description().NxNyNz), gr, pg);
  }
  void integralgrad(double kT, const GridDescription &gd, const VectorXd &x, VectorXd *g, VectorXd *pg=0) const {
    grad(kT, gd, x, gd.dvolume*VectorXd::Ones(gd.NxNyNz), g, pg);
  }
  void integralpgrad(double kT, const Grid &g, VectorXd *gr) const {
    pgrad(kT, g.description(), g, g.description().dvolume*VectorXd::Ones(g.description().NxNyNz), gr);
  }
  void integralpgrad(double kT, const GridDescription &gd, const VectorXd &x, VectorXd *g) const {
    pgrad(kT, gd, x, gd.dvolume*VectorXd::Ones(gd.NxNyNz), g);
  }
  double operator()(double kT, double data) const {
    assert(itsCounter);
    assert(itsCounter->ptr);
    double out = itsCounter->ptr->transform(kT, data);
    if (mynext) out += (*mynext)(kT, data);
    return out;
  }
  Functional grad_T(const Functional &ingrad) const {
    if (mynext) {
      return itsCounter->ptr->grad_T(ingrad) + mynext->grad_T(ingrad);
    } else {
      return itsCounter->ptr->grad_T(ingrad);
    }
  }
  Functional grad(const Functional &ingrad, const Functional &x, bool ispgrad) const {
    if (mynext) {
      return itsCounter->ptr->grad(ingrad, x, ispgrad) + mynext->grad(ingrad, x, ispgrad);
    } else {
      return itsCounter->ptr->grad(ingrad, x, ispgrad);
    }
  }
  Expression derive_homogeneous(const Expression &x) const {
    if (mynext)
      return itsCounter->ptr->derive_homogeneous(x) + mynext->derive_homogeneous(x);
    else return itsCounter->ptr->derive_homogeneous(x);
  }
  Functional pgrad(const Functional &ingrad, const Functional &x) const {
    return grad(ingrad, x, true);
  }
  void grad(double kT, const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    grad(gd, kT, data, ingrad, outgrad, outpgrad);
  }
  void grad(const GridDescription &gd, double kT, const VectorXd &data,
            const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    itsCounter->ptr->grad(gd, kT, data, ingrad, outgrad, outpgrad);
    if (mynext) mynext->grad(gd, kT, data, ingrad, outgrad, outpgrad);
  }
  void pgrad(double kT, const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
             VectorXd *outpgrad) const {
    pgrad(gd, kT, data, ingrad, outpgrad);
  }
  void pgrad(const GridDescription &gd, double kT, const VectorXd &data,
             const VectorXd &ingrad, VectorXd *outpgrad) const {
    itsCounter->ptr->pgrad(gd, kT, data, ingrad, outpgrad);
    if (mynext) mynext->pgrad(gd, kT, data, ingrad, outpgrad);
  }
  double derive(double kT, double data) const {
    double out = itsCounter->ptr->derive(kT, data);
    if (mynext) out += mynext->derive(kT, data);
    return out;
  }
  double d_by_dT(double kT, double data) const {
    double out = itsCounter->ptr->d_by_dT(kT, data);
    if (mynext) out += mynext->d_by_dT(kT, data);
    return out;
  }
  const std::string get_name() const { return itsCounter->name; }
  Functional set_name_stdstring(const std::string &n) { itsCounter->name = n; return *this; }
  Functional set_name(const char *n) {
    try {
      if (n) itsCounter->name = n;
      else itsCounter->name = "";
    }
    catch (...) {
      assert(0);
    }
    return *this;
  }
  Functional *next() const {
    return mynext;
  }
  Functional set_last_energy(double e) const { itsCounter->last_energy = e; return *this; }

  void print_summary(const char *prefix, double energy, std::string name="") const;
  // The following utility methods do not need to be overridden.  Its
  // return value is the total energy that is printed.
  double print_iteration(const char *prefix, int iter) const;
  // run_finite_difference_test returns non-zero when the test fails.
  int run_finite_difference_test(const char *testname,
                                 double kT, const Grid &data,
                                 const VectorXd *direction = 0) const;
  int run_homogeneous_finite_difference_test(const char *testname,
                                             double kT, double data) const;
  Expression printme(const Expression &) const;
  void create_source(const std::string filename, const std::string classname,
                     const char *args[], bool isheader=false) const;
  void create_header(const std::string filename, const std::string classname,
                     const char *args[]) const {
    create_source(filename, classname, args, true);
  }
  bool I_have_analytic_grad() const {
    const Functional *nxt = this;
    while (nxt) {
      if (!nxt->itsCounter->ptr->I_have_analytic_grad()) return false;
      nxt = nxt->next();
    }
    return true;
  }
  bool I_am_homogeneous() const {
    const Functional *nxt = this;
    while (nxt) {
      if (!nxt->itsCounter->ptr->I_am_homogeneous()) return false;
      nxt = nxt->next();
    }
    return true;
  }
  bool I_preserve_homogeneous() const {
    const Functional *nxt = this;
    while (nxt) {
      if (!nxt->itsCounter->ptr->I_am_homogeneous() &&
          !nxt->itsCounter->ptr->I_preserve_homogeneous()) return false;
      nxt = nxt->next();
    }
    return true;
  }
  bool I_am_constant_wrt_x() const {
    const Functional *nxt = this;
    while (nxt) {
      if (!nxt->itsCounter->ptr->I_am_constant_wrt_x()) return false;
      nxt = nxt->next();
    }
    return true;
  }
  bool I_give_zero_for_zero() const {
    const Functional *nxt = this;
    while (nxt) {
      if (!nxt->itsCounter->ptr->I_give_zero_for_zero()) return false;
      nxt = nxt->next();
    }
    return true;
  }
  bool I_am_zero() const {
    const Functional *nxt = this;
    while (nxt) {
      assert(nxt->itsCounter);
      if (!nxt->itsCounter->ptr->I_am_zero()) return false;
      nxt = nxt->next();
    }
    return true;
  }
  bool I_am_one() const {
    const Functional *nxt = this;
    if (!itsCounter->ptr->I_am_one() || nxt) return false;
    return true;
  }
  bool I_am_local() const {
    const Functional *nxt = this;
    while (nxt) {
      if (!nxt->itsCounter->ptr->I_am_local()) return false;
      nxt = nxt->next();
    }
    return true;
  }
private:
  void init(FunctionalInterface *p, const std::string &name) {
    if (p) {
      itsCounter = new counter(p);
      itsCounter->name = name;
    }
    mynext = 0;
  }
  Functional *mynext;
  struct counter {
    counter(FunctionalInterface* p = 0, unsigned c = 1) : ptr(p), count(c) {}
    FunctionalInterface* ptr;
    unsigned count;
    std::string name;
    Functional *next;
    double last_energy;
  };
  counter *itsCounter;
  void acquire(counter* c) {
    // increment the count
    itsCounter = c;
    if (c) (c->count)++;
  }
  void release() {
    // decrement the count, delete if it is 0
    if (itsCounter) {
      if (--itsCounter->count == 0) {
        delete itsCounter->ptr;
        delete itsCounter;
      }
      itsCounter = 0;
    }
  }
};

inline Functional operator*(double x, const Functional &f) {
//  return f*x;  modified 3/14
  return f*Functional(x);
}
inline Functional operator-(double x, const Functional &f) {
  return Functional(x) - f;
}

Functional log(const Functional &);
Functional exp(const Functional &);
Functional abs(const Functional &);
Functional sqr(const Functional &);
Functional sqrt(const Functional &);
Functional constrain(const Grid &, Functional);

Functional dV();



template<typename Derived, typename extra>
class ConvolveWith : public FunctionalInterface {
public:
  ConvolveWith(Derived (*ff)(const GridDescription &, extra),
               extra e, const Expression &R, const Expression &gzerovv, bool isev)
    : f(ff), radexpr(R), gzerov(gzerovv), data(e), iseven(isev) {}
  ConvolveWith(const ConvolveWith &cw)
    : f(cw.f), radexpr(cw.radexpr), gzerov(cw.gzerov), data(cw.data), iseven(cw.iseven) {}
  bool I_am_local() const {
    return false;
  }
  bool I_preserve_homogeneous() const {
    return true;
  }
  bool I_give_zero_for_zero() const {
    return true;
  }

  EIGEN_STRONG_INLINE VectorXd transform(const GridDescription &gd, double, const VectorXd &x) const {
    Grid out(gd, x);
    ReciprocalGrid recip = out.fft();
    recip.cwise() *= Eigen::CwiseNullaryOp<Derived, VectorXcd>(gd.NxNyNzOver2, 1, f(gd, data));
    return recip.ifft();
  }
  double gzero() const {
    Lattice lat(Cartesian(1,0,0), Cartesian(0,1,0), Cartesian(0,0,1));
    GridDescription gd(lat, 2, 2, 2);
    return f(gd, data).func(Reciprocal(0,0,0)).real();
  }
  double transform(double, double n) const {
    return n*gzero();
  }
  double derive(double, double) const {
    return gzero();
  }
  double d_by_dT(double, double) const {
    return 0;
  }
  Expression derive_homogeneous(const Expression &) const {
    return gzerov;
  }
  Functional grad(const Functional &ingrad, const Functional &, bool) const {
    if (iseven)
      return Functional(new ConvolveWith(f, data, radexpr, gzerov, iseven))(ingrad);
    else
      return Functional(new ConvolveWith(f, data, radexpr, gzerov, iseven))(-1*ingrad);
  }
  Functional grad_T(const Functional &) const {
    // FIXME: I assume here that the convolution kernel itself doesn't
    // depend on temperature, which may not be the case.
    return Functional(0.0);
  }
  EIGEN_STRONG_INLINE void grad(const GridDescription &gd, double, const VectorXd &,
                                const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    Grid out(gd, ingrad);
    ReciprocalGrid recip = out.fft();
    recip.cwise() *= Eigen::CwiseNullaryOp<Derived, VectorXcd>(gd.NxNyNzOver2, 1, f(gd, data));
    if (iseven) out = recip.ifft();
    else out = -recip.ifft();
    *outgrad += out;
    // FIXME: we will want to propogate preexisting preconditioning
    if (outpgrad) *outpgrad += out;
  }
  Expression printme(const Expression &x) const {
    if (x.typeIs("double")) {
      return gzerov*x;
    } else {
      Lattice lat(Cartesian(1,0,0), Cartesian(0,1,0), Cartesian(0,0,1));
      Derived c(GridDescription(lat, 2, 2, 2), data);
      return ifft(funexpr(c.name(), Expression("gd"), radexpr).set_type("ReciprocalGrid") * fft(x));
    }
  }
private:
  Derived (*f)(const GridDescription &, extra);
  Expression radexpr, gzerov;
  extra data;
  bool iseven;
};
