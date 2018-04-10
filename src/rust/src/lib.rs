//! This is the deft crate!

extern crate internment;
extern crate tinyset;

use tinyset::{Map64, Fits64};
use internment::Intern;

/// This is an expression.
///
/// It implements `Copy` and so is efficient to use.
#[derive(Debug, Clone, Copy, Eq, Hash)]
pub struct Expr {
    inner: Intern<InnerExpr>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct InnerExpr {
    et: ExprType,
    op: ExprOp,
}

/// This is the type an `Expr` will have in generated code.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ExprType {
    Scalar,
    Vector,
    RealSpace,
    RealSpaceVector,
    KSpace,
    KSpaceVector,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum ExprOp {
    Var(&'static str),
    Exp(Expr),
    Log(Expr),
    Sum(CommutativeMap),
    Mul(CommutativeMap),
}

#[derive(Debug, Clone)]
struct CommutativeMap {
    inner: Map64<Expr, f64>,
}

impl Expr {
    /// Express a variable.
    pub fn var(sym: &'static str, et: ExprType) -> Expr {
        Expr::from_inner(InnerExpr {
            et: et,
            op: ExprOp::Var(sym),
        })
    }

    /// Express an exponential.
    pub fn exp(self) -> Expr {
        Expr::from_inner(InnerExpr {
            et: self.get_type(),
            op: ExprOp::Exp(self),
        })
    }

    /// Express a logarithm.
    pub fn log(self) -> Expr {
        Expr::from_inner(InnerExpr {
            et: self.get_type(),
            op: ExprOp::Log(self),
        })
    }

    /// Raise an expression to a constant power.
    pub fn pow<T: Into<f64>>(self, power: T) -> Expr {
        let power = power.into();
        if power == 0.0 {
            Expr::one(self.get_type())
        } else if power == 1.0 {
            self
        } else {
            Expr::mul_from_map(CommutativeMap::from_pair(self, power), self.get_type())
        }
    }

    /// Differentiate symbolically with respect to `wrt`.
    pub fn deriv(&self, wrt: Expr) -> Expr {
        if *self == wrt {
            return Expr::one(self.get_type());
        }
        match &self.inner.op {
            &ExprOp::Var(..) => Expr::zero(self.get_type()),
            &ExprOp::Exp(a) => *self * a.deriv(wrt),
            &ExprOp::Log(a) => a.deriv(wrt) / a,
            &ExprOp::Sum(ref a) =>
                a.inner.iter()
                       .fold(Expr::zero(self.get_type()),
                             |a, (b, &c)| a + Expr::pow(b, c - 1.0) * c * b.deriv(wrt)),
            &ExprOp::Mul(ref a) =>
                a.inner.iter()
                       .fold(Expr::zero(self.get_type()),
                             |a, (b, &c)| a + b.deriv(wrt) * c * *self / b),
        }
    }

    fn zero(et: ExprType) -> Expr {
        Expr::sum_from_map(CommutativeMap::new(), et)
    }

    fn one(et: ExprType) -> Expr {
        Expr::mul_from_map(CommutativeMap::new(), et)
    }

    fn from_inner(i: InnerExpr) -> Expr {
        Expr { inner: Intern::new(i), }
    }

    fn sum_from_map(map: CommutativeMap, et: ExprType) -> Expr {
        Expr::from_inner(InnerExpr {
            et: et,
            op: ExprOp::Sum(map),
        })
    }

    fn mul_from_map(map: CommutativeMap, et: ExprType) -> Expr {
        Expr::from_inner(InnerExpr {
            et: et,
            op: ExprOp::Mul(map),
        })
    }

    fn get_type(&self) -> ExprType {
        self.inner.et
    }

    /// Creates a fragment of C++ code which evalutes the expression.
    ///
    /// If `assignto` is not `""`, the C++ code will also instantiate a variable
    /// named the value of `assignto` and assign the result of the expression to
    /// it.
    pub fn cpp(&self, assignto: &'static str) -> String {
        if assignto != "" {
            match self.get_type() {
                ExprType::Scalar =>
                    String::from("double ") + assignto + &" = " + &self.cpp(""),
                ExprType::Vector =>
                    String::from("double ") + assignto
                        + &"[3]; for (int i = 0; i < 3; ++i) { "
                        + assignto + &"[i] = " + &self.cpp("") + &"; }",
                _ =>
                    String::from("\n#error unimplemented"),
            }
        } else {
            match &self.inner.op {
                &ExprOp::Var(s) =>
                    match self.get_type() {
                        ExprType::Vector => String::from(s) + &"[i]",
                        _ => String::from(s),
                    },
                &ExprOp::Exp(a) => String::from("exp(") + &a.cpp("") + &")",
                &ExprOp::Log(a) => String::from("log(") + &a.cpp("") + &")",
                &ExprOp::Sum(ref m) => m.sum_string(),
                &ExprOp::Mul(ref m) => m.mul_string(),
            }
        }
    }
}

impl<RHS: Into<Expr>> std::ops::Add<RHS> for Expr {
    type Output = Self;

    fn add(self, other: RHS) -> Self {
        let other = other.into();

        assert_eq!(self.get_type(), other.get_type());

        let mut sum: CommutativeMap;
        match (&self.inner.op, &other.inner.op) {
            (&ExprOp::Sum(ref lhs), &ExprOp::Sum(ref rhs)) => {
                sum = lhs.clone();
                sum.union(&rhs);
            },

            (&ExprOp::Sum(ref lhs), _) => {
                sum = lhs.clone();
                if other != Expr::zero(other.get_type()) {
                    sum.insert(other, 1.0);
                }
            },

            (_, &ExprOp::Sum(ref rhs)) => {
                sum = rhs.clone();
                if self != Expr::zero(self.get_type()) {
                    sum.insert(self, 1.0);
                }
            },

            (_, _) => {
                if self != Expr::zero(self.get_type()) {
                    sum = CommutativeMap::from_pair(self, 1.0);
                } else {
                    sum = CommutativeMap::new()
                }
                if other != Expr::zero(other.get_type()) {
                    sum.insert(other, 1.0);
                }
            },
        }

        if sum.inner.len() == 1 {
            let (x, &c) = sum.inner.iter().next().unwrap();
            if c == 1.0 {
                return x
            }
        }
        Expr::sum_from_map(sum, self.get_type())
    }
}

impl<RHS: Into<Expr>> std::ops::Mul<RHS> for Expr {
    type Output = Self;

    fn mul(self, other: RHS) -> Self {
        let other = other.into();

        assert_eq!(self.get_type(), other.get_type());

        let mut mul: CommutativeMap;
        match (&self.inner.op, &other.inner.op) {
            (&ExprOp::Mul(ref lhs), &ExprOp::Mul(ref rhs)) => {
                mul = lhs.clone();
                mul.union(&rhs);
            },

            (&ExprOp::Mul(ref lhs), _) => {
                mul = lhs.clone();
                if other != Expr::one(other.get_type()) {
                    if other == Expr::zero(other.get_type()) {
                        return Expr::zero(other.get_type());
                    }
                    mul.insert(other, 1.0);
                }
            },

            (_, &ExprOp::Mul(ref rhs)) => {
                mul = rhs.clone();
                if self != Expr::one(self.get_type()) {
                    if self == Expr::zero(self.get_type()) {
                        return Expr::zero(self.get_type());
                    }
                    mul.insert(self, 1.0);
                }
            },

            (_, _) => {
                if self == Expr::zero(self.get_type()) || other == Expr::zero(other.get_type()) {
                    return Expr::zero(self.get_type());
                }
                if self != Expr::one(self.get_type()) {
                    mul = CommutativeMap::from_pair(self, 1.0);
                } else {
                    mul = CommutativeMap::new()
                }
                if other != Expr::one(other.get_type()) {
                    mul.insert(other, 1.0);
                }
            },
        }

        if mul.inner.len() == 1 {
            let (x, &p) = mul.inner.iter().next().unwrap();
            if p == 1.0 {
                return x
            }
        }
        Expr::mul_from_map(mul, self.get_type())
    }
}

impl<RHS: Into<Expr>> std::ops::Sub<RHS> for Expr {
    type Output = Self;

    fn sub(self, other: RHS) -> Self {
        let other = other.into();
        if other == Expr::zero(other.get_type()) {
            return self;
        }
        match &other.inner.op {
            &ExprOp::Sum(ref map) =>
                self + Expr::sum_from_map(map.negate(), other.get_type()),
            _ =>
                self + Expr::sum_from_map(CommutativeMap::from_pair(other, -1.0), other.get_type()),
        }
    }
}

impl<RHS: Into<Expr>> std::ops::Div<RHS> for Expr {
    type Output = Self;

    fn div(self, other: RHS) -> Self {
        let other = other.into();
        assert!(other != Expr::zero(other.get_type()));
        if other == Expr::one(other.get_type()) {
            return self;
        }
        match &other.inner.op {
            &ExprOp::Mul(ref map) =>
                self * Expr::mul_from_map(map.negate(), other.get_type()),
            _ =>
                self * Expr::mul_from_map(CommutativeMap::from_pair(other, -1.0), other.get_type()),
        }
    }
}

impl std::ops::Neg for Expr {
    type Output = Self;

    fn neg(self) -> Self {
        match &self.inner.op {
            &ExprOp::Sum(ref sum) =>
                Expr::sum_from_map(sum.negate(), self.get_type()),
            _ =>
                Expr::sum_from_map(CommutativeMap::from_pair(self, -1.0), self.get_type()),
        }
    }
}

impl std::ops::Sub<Expr> for f64 {
    type Output = Expr;

    fn sub(self, other: Expr) -> Expr {
        let x: Expr = self.into();
        x - other
    }
}

impl std::ops::Sub<Expr> for i32 {
    type Output = Expr;

    fn sub(self, other: Expr) -> Expr {
        let x: Expr = self.into();
        x - other
    }
}

impl std::ops::Div<Expr> for f64 {
    type Output = Expr;

    fn div(self, other: Expr) -> Expr {
        let x: Expr = self.into();
        x / other
    }
}

impl std::ops::Div<Expr> for i32 {
    type Output = Expr;

    fn div(self, other: Expr) -> Expr {
        let x: Expr = self.into();
        x / other
    }
}

impl CommutativeMap {
    fn new() -> Self {
        CommutativeMap { inner: Map64::new(), }
    }

    fn as_vec(&self) -> Vec<(u64, u64)> {
        self.inner.iter().map(|(k, &v)| (k.inner.to_u64(), v as u64)).collect()
    }

    // This lets you create a `CommutativeMap` without wasting time checking if
    // the first term was already there.
    fn from_pair(k: Expr, v: f64) -> CommutativeMap {
        let mut map = CommutativeMap::new();
        map.inner.insert(k, v);
        map
    }

    fn insert(&mut self, k: Expr, v: f64) {
        let v = self.inner.get(&k).unwrap_or(&0.0) + v;
        if v == 0.0 {
            self.inner.remove(&k);
        } else {
            self.inner.insert(k, v);
        }
    }

    fn union(&mut self, other: &CommutativeMap) {
        for (k, &v) in other.inner.iter() {
            self.insert(k, v);
        }
    }

    // This is used for negation, subtraction, reciprocals, division.
    fn negate(&self) -> Self {
        let mut neg = CommutativeMap::new();
        for (k, &v) in self.inner.iter() {
            neg.inner.insert(k, -v);
        }
        neg
    }

    // Returns a result useful for generating code.
    // Splits the map by sign, sorts each by c++ string, and returns the two vectors;
    // .0 are the positives, .1 are the negatives
    fn pre_string(&self) -> (Vec<(String, f64)>, Vec<(String, f64)>) {
        let mut plusses: Vec<(String, f64)>;
        let mut minuses: Vec<(String, f64)>;
        plusses = self.inner.iter().filter_map(|(k, &v)| {
            if v > 0.0 {
                Some((k.cpp(""), v))
            } else {
                None
            }
        }).collect();
        minuses = self.inner.iter().filter_map(|(k, &v)| {
            if v < 0.0 {
                Some((k.cpp(""), v))
            } else {
                None
            }
        }).collect();
        plusses.sort_by(|&(ref p, _), &(ref q, _)| p.cmp(&q));
        minuses.sort_by(|&(ref p, _), &(ref q, _)| p.cmp(&q));
        // Would use .sort() here but can't sort f64
        (plusses, minuses)
    }

    fn sum_string(&self) -> String {
        let (p, n) = self.pre_string();
        let positive_coefficient = |&(ref x, ref c): &(String, f64)| -> String {
            if x == "1" {
                c.to_string()
            } else if *c == 1.0 {
                x.clone()
            } else {
                c.to_string() + &" * " + &x
            }
        };
        let negative_coefficient = |&(ref x, ref c): &(String, f64)| -> String {
            if x == "1" {
                c.abs().to_string()
            } else if *c == -1.0 {
                x.clone()
            } else {
                c.abs().to_string() + &" * " + &x
            }
        };
        match (p.len(), n.len()) {
            (0, 0) => String::from("0"),
            (_, 0) => p.iter()
                       .map(positive_coefficient)
                       .collect::<Vec<String>>()
                       .join(" + "),
            (0, _) => String::from("-")
                      + &n.iter()
                          .map(negative_coefficient)
                          .collect::<Vec<String>>()
                          .join(" - "),
            (_, _) => p.iter()
                       .map(positive_coefficient)
                       .collect::<Vec<String>>()
                       .join(" + ")
                      + &" - "
                      + &n.iter()
                          .map(negative_coefficient)
                          .collect::<Vec<String>>()
                          .join(" - "),
        }
    }

    fn mul_string(&self) -> String {
        let (n, d) = self.pre_string();
        let power = |&(ref x, ref p): &(String, f64)| -> String {
            if x == "1" || p.abs() == 1.0 {
                x.clone()
            } else if p.abs() == 2.0 {
                x.clone() + &" * " + &x
            } else {
                String::from("pow(") + &x + &", " + &p.abs().to_string() + &")"
            }
        };
        match (n.len(), d.len()) {
            (0, 0) => String::from("1"),
            (_, 0) => n.iter()
                       .map(power)
                       .collect::<Vec<String>>()
                       .join(" * "),
            (0, _) => String::from("1 / (")
                      + &d.iter()
                          .map(power)
                          .collect::<Vec<String>>()
                          .join(" * ")
                      + &")",
            (_, _) => n.iter()
                       .map(&power)
                       .collect::<Vec<String>>()
                       .join(" * ")
                      + &" / "
                      + &d.iter()
                          .map(&power)
                          .collect::<Vec<String>>()
                          .join(" * ")
                      + &")",
        }
    }
}

/// This lets numbers easily become `Expr`s of type `ExprType::Scalar`.
impl<N: Into<f64>> From<N> for Expr {
    fn from(n: N) -> Self {
        let n = n.into();
        if n == 0.0 {
            Expr::sum_from_map(CommutativeMap::new(), ExprType::Scalar)
        } else if n == 1.0 {
            Expr::mul_from_map(CommutativeMap::new(), ExprType::Scalar)
        } else {
            Expr::sum_from_map(CommutativeMap::from_pair(1.into(), n), ExprType::Scalar)
        }
    }
}

/// This lets you compare `Expr`s for equality.
impl PartialEq for Expr {
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

// This lets `Expr`s be keys in `CommutativeMap`s.
impl Fits64 for Expr {
    fn to_u64(self) -> u64 {
        self.inner.to_u64()
    }

    unsafe fn from_u64(x: u64) -> Self {
        Expr { inner: Intern::from_u64(x), }
    }
}

// Everything in `InnerExpr` must be `Hash` so that `InnerExpr` can be
// interned, but only `CommutativeMap` needs a manual implementation.
impl std::hash::Hash for CommutativeMap {
    fn hash<H: std::hash::Hasher>(&self, hasher: &mut H) {
        let mut terms = self.as_vec();
        terms.sort();
        for (k, v) in terms {
            k.hash(hasher);
            v.hash(hasher);
        }
    }
}

// Everything in `InnerExpr` must be `Eq` so that `InnerExpr` can be
// interned, but only `CommutativeMap` needs a manual implementation.
impl PartialEq for CommutativeMap {
    fn eq(&self, other: &Self) -> bool {
        let (mut lhs, mut rhs) = (self.as_vec(), other.as_vec());
        lhs.sort();
        rhs.sort();
        lhs == rhs
    }
}

impl Eq for CommutativeMap {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basics() {
        let a = Expr::var("a", ExprType::Scalar);
        let b = Expr::var("b", ExprType::Scalar);

        assert_eq!(Expr::from(1) + 0, Expr::from(1));
        assert_eq!(a + 0, a);
        assert_eq!(a, a);
        assert!(a != b);
        assert_eq!(a + a, a + a);
        assert_eq!(a + b, a + b);
        assert_eq!(a + b, b + a);
        assert_eq!(a * a, a * a);
        assert_eq!(a * b, a * b);
        assert_eq!(a * b, b * a);
        assert!(a + a != b + b);
        assert!(a * a != b * b);
        assert_eq!(a - a, 0.into());
        assert_eq!(a / a, 1.into());
        assert!(a - b != 0.into());
        assert!(a / b != 0.into());
        assert_eq!(a - a - a, Expr::from(0) - a);
        assert_eq!(a - a - a, 0 - a);
        assert_eq!(a - a - a, -a);
        assert_eq!(a / a, 1.into());
        assert_eq!(a / a / a, Expr::from(1) / a);
        assert_eq!(a / a / a, 1 / a);
    }

    #[test]
    fn strings() {
        let a = Expr::var("a", ExprType::Scalar);
        let u = Expr::var("u", ExprType::Vector);

        assert_eq!(a.cpp(""), "a");
        assert_eq!(u.cpp("t"), "double t[3]; for (int i = 0; i < 3; ++i) { t[i] = u[i]; }");
    }

    #[test]
    fn derivatives() {
        let a = Expr::var("a", ExprType::Scalar);
        let b = Expr::var("b", ExprType::Scalar);

        assert_eq!(a.deriv(a), Expr::from(1));
        assert_eq!((a * a).deriv(a), a * 2);
        assert_eq!((a * a * a).deriv(a), a*a*3);
        assert_eq!(a.deriv(b), 0.into());
        assert_eq!(a.deriv(b), a * 0);
        assert_eq!((a + b).deriv(a), Expr::from(1));

        assert_eq!(Expr::log(a).deriv(a), 1.0 / a);
        assert_eq!(Expr::log(a * a).deriv(a), 2 / a);
        assert_eq!(Expr::exp(a).deriv(a), Expr::exp(a));
        assert_eq!(Expr::exp(a * a).deriv(a), a * Expr::exp(Expr::pow(a, 2)) * 2);
        assert_eq!(Expr::pow(Expr::log(a) + a, 3).deriv(a),
                   (a.log() + a).pow(2) * 3 * (Expr::from(1) / a + 1));
        assert_eq!((a.log() + Expr::log(a) * Expr::log(a)).deriv(Expr::log(a)),
                   Expr::from(1) + Expr::log(a) * 2);
    }
}
