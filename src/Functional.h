// -*- mode: C++; -*-

#pragma once

#include "ReciprocalGrid.h"
#include "Expression.h"

class Functional;

class FunctionalInterface {
public:
  FunctionalInterface() : have_analytic_grad(true), have_integral(false) {}
  virtual ~FunctionalInterface() {}

  // A functional mapping one field onto another...
  virtual double integral(const GridDescription &gd, const VectorXd &data) const;
  virtual VectorXd transform(const GridDescription &gd, const VectorXd &data) const = 0;
  virtual double transform(double) const = 0;

  // This computes the gradient of the functional, given a gradient of
  // its output field (i.e. it applies the chain rule).
  virtual void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
                    VectorXd *outgrad, VectorXd *outpgrad) const = 0;
  virtual double derive(double data) const = 0;
  virtual Functional grad(const Functional &ingrad, const Functional &x, bool ispgrad) const = 0;
  virtual Expression printme(const Expression &) const = 0;

  virtual void print_summary(const char *prefix, double energy, const char *name) const;

  bool have_analytic_grad, have_integral;
};

template<typename Derived, typename extra> class ConvolveWith;

class Functional {
public:
  // Handle reference counting so we can pass these things around freely...
  Functional(double, const char *name=0); // This handles constants!
  explicit Functional(const VectorXd &); // This handles constant fields!
  template<typename Derived, typename extra>
  explicit Functional(Derived (*f)(const GridDescription &, extra), extra e, bool iseven)
    : itsCounter(0) {
    Lattice lat(Cartesian(1,0,0), Cartesian(0,1,0), Cartesian(0,0,1));
    GridDescription gd(lat, 2, 2, 2);
    // This handles constant ephemeral fields!
    init(new ConvolveWith<Derived,extra>(f,e,iseven), f(gd, e).name());
  }
  explicit Functional(FunctionalInterface* p = 0, const char *name = 0) // allocate a new counter
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
    return *this;
  }

  Functional operator()(const Functional &) const;
  Functional operator+(const Functional &x) const {
    Functional out(*this);
    return out += x;
  }
  Functional operator+=(const Functional &x) {
    if (mynext) *mynext += x;
    else mynext = new Functional(x);
    return *this;
  }
  Functional operator-() const;
  Functional operator-(const Functional &) const;
  Functional operator/(const Functional &) const;
  Functional operator*(const Functional &) const;
  VectorXd justMe(const GridDescription &gd, const VectorXd &data) const {
    return itsCounter->ptr->transform(gd, data);
  }
  VectorXd operator()(const GridDescription &gd, const VectorXd &data) const {
    VectorXd out = itsCounter->ptr->transform(gd, data);
    if (mynext) out += (*mynext)(gd, data);
    return out;
  }
  VectorXd operator()(const Grid &g) const {
    VectorXd out = itsCounter->ptr->transform(g.description(), g);
    if (mynext) out += (*mynext)(g.description(), g);
    return out;
  }
  double integral(const Grid &g) const {
    return integral(g.description(), g);
  }
  double integral(const GridDescription &gd, const VectorXd &data) const {
    if (itsCounter->ptr->have_integral && !next()) {
      return itsCounter->ptr->integral(gd, data);
    }
    // This does some extra work to save the energies of each term in
    // the sum.
    VectorXd fdata(justMe(gd, data));
    double e = gd.dvolume*fdata.sum();
    set_last_energy(e);
    Functional *nxt = next();
    while (nxt) {
      fdata += nxt->justMe(gd, data);
      double etot = gd.dvolume*fdata.sum();
      nxt->set_last_energy(etot - e);
      e = etot;
      nxt = nxt->next();
    }
    return e;
  }
  void integralgrad(const Grid &g, VectorXd *gr, VectorXd *pg=0) const {
    grad(g.description(), g, g.description().dvolume*VectorXd::Ones(g.description().NxNyNz), gr, pg);
  }
  void integralgrad(const GridDescription &gd, const VectorXd &x, VectorXd *g, VectorXd *pg=0) const {
    grad(gd, x, gd.dvolume*VectorXd::Ones(gd.NxNyNz), g, pg);
  }
  double operator()(double data) const {
    double out = itsCounter->ptr->transform(data);
    if (mynext) out += (*mynext)(data);
    return out;
  }
  Functional grad(const Functional &ingrad, const Functional &x, bool ispgrad) const {
    if (mynext) {
      return itsCounter->ptr->grad(ingrad, x, ispgrad) + mynext->grad(ingrad, x, ispgrad);
    } else {
      return itsCounter->ptr->grad(ingrad, x, ispgrad);
    }
  }
  Functional pgrad(const Functional &ingrad, const Functional &x) const {
    return grad(ingrad, x, true);
  }
  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    itsCounter->ptr->grad(gd, data, ingrad, outgrad, outpgrad);
    if (mynext) mynext->grad(gd, data, ingrad, outgrad, outpgrad);
  }
  double derive(double data) const {
    double out = itsCounter->ptr->derive(data);
    if (mynext) out += mynext->derive(data);
    return out;
  }
  const char *get_name() const { return itsCounter->name; }
  Functional set_name(const char *n) { itsCounter->name = n; return *this; }
  Functional *next() const {
    return mynext;
  }
  Functional set_last_energy(double e) const { itsCounter->last_energy = e; return *this; }

  void print_summary(const char *prefix, double energy, const char *name=0) const;
  // The following utility methods do not need to be overridden.  Its
  // return value is the total energy that is printed.
  double print_iteration(const char *prefix, int iter) const;
  // run_finite_difference_test returns non-zero when the test fails.
  int run_finite_difference_test(const char *testname,
                                  const Grid &data,
                                  const VectorXd *direction = 0) const;
  Expression printme(const Expression &) const;
  void create_source(const std::string filename, const std::string classname,
                     const char *a1 = 0, const char *a2 = 0, bool isheader=false) const;
  void create_header(const std::string filename, const std::string classname,
                     const char *a1 = 0, const char *a2 = 0) const {
    create_source(filename, classname, a1, a2, true);
  }
private:
  void init(FunctionalInterface *p, const char *name) {
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
    const char *name;
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
  return f*x;
}
inline Functional operator-(double x, const Functional &f) {
  return Functional(x) - f;
}

Functional log(const Functional &);
Functional exp(const Functional &);
Functional sqr(const Functional &);
Functional constrain(const Grid &, Functional);

// choose combines two local functionals, with which one is applied
// depending on the local value of the field.
Functional choose(double, const Functional &lower, const Functional &higher);

extern Functional dV;



template<typename Derived, typename extra>
class ConvolveWith : public FunctionalInterface {
public:
  ConvolveWith(Derived (*ff)(const GridDescription &, extra), extra e, bool isev)
    : f(ff), data(e), iseven(isev) {}
  ConvolveWith(const ConvolveWith &cw) : f(cw.f), data(cw.data), iseven(cw.iseven) {}

  EIGEN_STRONG_INLINE VectorXd transform(const GridDescription &gd, const VectorXd &x) const {
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
  double transform(double n) const {
    return n*gzero();
  }
  double derive(double) const {
    return gzero();
  }
  Functional grad(const Functional &ingrad, const Functional &, bool) const {
    if (iseven)
      return Functional(new ConvolveWith(f, data, iseven))(ingrad);
    else
      return Functional(new ConvolveWith(f, data, iseven))(-1*ingrad);
  }
  EIGEN_STRONG_INLINE void grad(const GridDescription &gd, const VectorXd &, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
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
    Lattice lat(Cartesian(1,0,0), Cartesian(0,1,0), Cartesian(0,0,1));
    Derived c(GridDescription(lat, 2, 2, 2), data);
    return ifft(funexpr(c.name(), Expression("gd"), Expression("R")).set_type("ReciprocalGrid") * fft(x));
  }
private:
  Derived (*f)(const GridDescription &, extra);
  extra data;
  bool iseven;
};
