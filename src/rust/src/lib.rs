//! This is the deft crate!

extern crate internment;
extern crate tinyset;

use tinyset::{Set64, Map64, Fits64};
use internment::Intern;

#[derive(Clone, Debug)]
struct CommutativeMap {
    map: Map64<Expr, f64>,
}

impl std::hash::Hash for CommutativeMap {
    fn hash<H>(&self, hasher: &mut H) where H: std::hash::Hasher {
        let mut terms: Vec<(u64, u64)>
            = self.map.iter()
                      .map(|(k, &v)| (k.inner.to_u64(), v as u64))
                      .collect();
        terms.sort();
        for (k, v) in terms {
            k.hash(hasher);
            v.hash(hasher);
        }
    }
}

impl std::cmp::PartialEq for CommutativeMap {
    fn eq(&self, other: &CommutativeMap) -> bool {
        for (k, v) in self.map.iter() {
            if other.get(k) != Some(v) {
                return false
            }
        }
        true
    }
}

impl std::cmp::Eq for CommutativeMap {}

impl CommutativeMap {
    fn new() -> CommutativeMap {
        CommutativeMap { map: Map64::new() }
    }

    fn len(&self) -> usize {
        self.map.len()
    }

    fn get(&self, k: Expr) -> Option<&f64> {
        self.map.get(&k)
    }

    fn insert(&mut self, k: Expr, v: f64) {
        let v = self.map.get(&k).unwrap_or(&0.0) + v;
        if v == 0.0 {
            self.map.remove(&k);
        } else {
            self.map.insert(k, v);
        }
    }

    fn union(&mut self, other: &CommutativeMap) {
        for (k, &v) in other.map.iter() {
            self.insert(k, v);
        }
    }

    fn negate(&self) -> CommutativeMap {
        let mut neg = CommutativeMap::new();
        for (k, &v) in self.map.iter() {
            neg.insert(k, -v);
        }
        neg
    }

    fn to_string<F>(&self, operator: &'static str, coeff: F) -> String where
        F: FnMut(&(String, f64)) -> String {
        let mut terms: Vec<(String, f64)>
            = self.map.iter()
                      .map(|(k, &v)| (k.cpp(), v))
                      .collect();
        terms.sort_by(|&(ref p, _), &(ref q, _)| p.cmp(&q));
        terms.iter()
             .map(coeff)
             .collect::<Vec<String>>()
             .join(operator)
    }

    fn split_by_sign(&self) -> (CommutativeMap, CommutativeMap) {
        let mut plusses = CommutativeMap::new();
        let mut minuses = CommutativeMap::new();

        for (k, &v) in self.map.iter() {
            if v > 0.0 {
                plusses.insert(k, v)
            } else if v < 0.0 {
                minuses.insert(k, v)
            }
        }

        (plusses, minuses)
    }
}

#[derive(Clone, PartialEq, Eq, Debug)]
struct Subexprs {
    set: Set64<Expr>,
}

impl std::hash::Hash for Subexprs {
    fn hash<H>(&self, hasher: &mut H) where H: std::hash::Hasher {
        let mut membs: Vec<u64>
            = self.set.iter()
                      .map(|i| i.inner.to_u64())
                      .collect();
        membs.sort();
        for i in membs {
            i.hash(hasher);
        }
    }
}

impl Subexprs {
    fn new_interned() -> Intern<Self> {
        Intern::new(Self { set: Set64::new() })
    }
}

impl From<Set64<Expr>> for Subexprs {
    fn from(set: Set64<Expr>) -> Self {
        Self { set: set }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
enum InnerExpr {
    Var(Intern<&'static str>),
    Exp(Expr),
    Log(Expr),
    Sum(Intern<CommutativeMap>),
    Mul(Intern<CommutativeMap>),
}

/// This is an expression. You can use it to do arithmetic.
///
/// # Example
///
/// ```
/// use deft::Expr;
/// assert_eq!(Expr::var("a"), Expr::var("a"));
/// ```
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct Expr {
    inner: Intern<InnerExpr>,
    subx: Intern<Subexprs>,
}

impl Fits64 for Expr {
    fn to_u64(self) -> u64 {
        self.inner.to_u64()
    }

    unsafe fn from_u64(x: u64) -> Self {
        Expr {
            inner: Intern::<InnerExpr>::from_u64(x),
            subx: Subexprs::new_interned(),
        }
    }
}

impl Expr {
    fn from_inner(i: InnerExpr) -> Expr {
        Expr {
            inner: Intern::new(i),
            subx: match i {
                InnerExpr::Var(_) => Subexprs::new_interned(),
                InnerExpr::Exp(a) | InnerExpr::Log(a) => a.subx.clone(),
                InnerExpr::Sum(m) | InnerExpr::Mul(m) => {
                    let subx = m.map.iter()
                                    .fold(Set64::new(), |a, (e, _)| &a | &e.subx.set);
                    Intern::new(Subexprs::from(subx))
                },
            }
        }
    }

    fn from_sum_map(m: CommutativeMap) -> Expr {
        if m.len() == 1 {
            let (k, &v) = m.map.iter().next().unwrap();
            if v == 1.0 {
                return k;
            }
        }
        Expr::from_inner(InnerExpr::Sum(Intern::new(m)))
    }

    fn from_mul_map(m: CommutativeMap) -> Expr {
        if m.get(Expr::zero()).is_some() {
            return Expr::zero();
        }
        if m.len() == 1 {
            let (k, &v) = m.map.iter().next().unwrap();
            if v == 1.0 {
                return k;
            }
        }
        Expr::from_inner(InnerExpr::Mul(Intern::new(m)))
    }

    /// Express a variable.
    pub fn var(sym: &'static str) -> Expr {
        Expr::from_inner(InnerExpr::Var(Intern::new(sym)))
    }

    /// Express an exponential.
    pub fn exp(arg: Expr) -> Expr {
        Expr::from_inner(InnerExpr::Exp(arg))
    }

    /// Express a logarithm.
    pub fn log(arg: Expr) -> Expr {
        Expr::from_inner(InnerExpr::Log(arg))
    }

    /// Express unity/one.
    pub fn one() -> Expr {
        Expr::from_inner(InnerExpr::Mul(Intern::new(CommutativeMap::new())))
    }

    /// Express zero.
    pub fn zero() -> Expr {
        Expr::from_inner(InnerExpr::Sum(Intern::new(CommutativeMap::new())))
    }

    /// Raise an expression to a constant power.
    ///
    /// # Example
    ///
    /// ```
    /// use deft::Expr;
    /// let a = Expr::var("a");
    /// assert_eq!(Expr::pow(a, 2) / a, a);
    /// ```
    pub fn pow<P: Into<f64>>(base: Expr, power: P) -> Expr {
        let mut mul = CommutativeMap::new();
        mul.insert(base, power.into());
        Expr::from_mul_map(mul)
    }

    /// Differentiate symbolically.
    ///
    /// # Example
    ///
    /// ```
    /// use deft::Expr;
    /// let x = Expr::var("x");
    /// assert_eq!((x * x).deriv("x"), x * 2);
    /// ```
    pub fn deriv(&self, wrt: &'static str) -> Expr {
        match *self.inner {
            InnerExpr::Var(s) => if *s == wrt { Expr::one() } else { Expr::zero() },
            InnerExpr::Exp(a) => *self * a.deriv(wrt),
            InnerExpr::Log(a) => a.deriv(wrt) / a,
            InnerExpr::Sum(a)
                => a.map.iter()
                        .fold(Expr::zero(),
                              |a, (b, &c)| a + Expr::pow(b, c - 1.0) * c * b.deriv(wrt)),
            InnerExpr::Mul(a)
                => a.map.iter()
                        .fold(Expr::zero(),
                              |a, (b, &c)| a + b.deriv(wrt) * c * *self / b),
            // InnerExpr::Mul(a)
            //     => a.map.iter()
            //             .map(|(expr, &power)| expr.deriv(wrt) * power * *self / expr).sum(),
        }
    }

    /// Creates a fragment of C++ code which evaluates the expression.
    pub fn cpp(&self) -> String {
        match *self.inner {
            InnerExpr::Var(s) => String::from(*s),
            InnerExpr::Exp(a) => String::from("exp(") + &a.cpp() + &")",
            InnerExpr::Log(a) => String::from("log(") + &a.cpp() + &")",
            InnerExpr::Sum(m) => {
                if (*m).len() == 0 {
                    return String::from("0");
                }

                let coeff = |&(ref s, c): &(String, f64)| -> String {
                    if s == "1" {
                        c.abs().to_string()
                    } else if c.abs() == 1.0 {
                        s.clone()
                    } else {
                        c.abs().to_string() + &" * " + &s
                    }
                };

                let (p, n) = (*m).split_by_sign();
                let p = p.to_string(" + ", &coeff);
                if n.len() != 0 {
                    let n = n.to_string(" - ", &coeff);
                    p + &" - " + &n
                } else {
                    p
                }
            },
            InnerExpr::Mul(m) => {
                if (*m).len() == 0 {
                    return String::from("1");
                }

                let (n, d) = (*m).split_by_sign();
                let power = |&(ref s, p): &(String, f64)| -> String {
                    if p.abs() == 1.0 {
                        s.clone()
                    } else {
                        String::from("pow(") + &s + &", " + &p.abs().to_string() + &")"
                    }
                };

                let n = if n.len() != 0 {
                    n.to_string(" * ", &power)
                } else {
                    String::from("1")
                };

                match d.len() {
                    0 => n,
                    1 => n + &" / " + &d.to_string(" * ", &power),
                    _ => n + &" / (" + &d.to_string(" * ", &power) + &")",
                }
            },
        }
    }
}

impl<N: Into<f64>> From<N> for Expr {
    fn from(n: N) -> Self {
        let mut m = CommutativeMap::new();
        let n = n.into();
        assert!(!n.is_nan());
        m.insert(Expr::one(), n);
        Expr::from_sum_map(m)
    }
}

impl<RHS: Into<Expr>> std::ops::Add<RHS> for Expr {
    type Output = Self;

    fn add(self, other: RHS) -> Self {
        let other = other.into();
        match (*self.inner, *other.inner) {
            (InnerExpr::Sum(a), InnerExpr::Sum(b)) => {
                let mut sum = (*a).clone();
                sum.union(&*b);
                Expr::from_sum_map(sum)
            },
            (InnerExpr::Sum(m), _) => {
                let mut sum = (*m).clone();
                sum.insert(other, 1.0);
                Expr::from_sum_map(sum)
            },
            (_, InnerExpr::Sum(m)) => {
                let mut sum = (*m).clone();
                sum.insert(self, 1.0);
                Expr::from_sum_map(sum)
            },
            _ => {
                let mut sum = CommutativeMap::new();
                sum.insert(self, 1.0);
                sum.insert(other, 1.0);
                Expr::from_sum_map(sum)
            },
        }
    }
}

impl<RHS: Into<Expr>> std::ops::Mul<RHS> for Expr {
    type Output = Self;

    fn mul(self, other: RHS) -> Self {
        let mut mul = CommutativeMap::new();
        let other = other.into();

        match *self.inner {
            InnerExpr::Mul(m) => mul.union(&*m),
            _ => mul.insert(self, 1.0),
        };

        match *other.inner {
            InnerExpr::Mul(m) => mul.union(&*m),
            _ => mul.insert(other, 1.0),
        };

        Expr::from_mul_map(mul)
    }
}

impl<RHS: Into<Expr>> std::ops::Sub<RHS> for Expr {
    type Output = Self;

    fn sub(self, other: RHS) -> Self {
        let mut subtrahend = CommutativeMap::new();
        let other = other.into();

        match *other.inner {
            InnerExpr::Sum(m) => subtrahend.union(&(*m).negate()),
            _ => subtrahend.insert(other, -1.0),
        };

        self + Expr::from_sum_map(subtrahend)
    }
}

impl<RHS: Into<Expr>> std::ops::Div<RHS> for Expr {
    type Output = Self;

    fn div(self, other: RHS) -> Self {
        let mut divisor = CommutativeMap::new();
        let other = other.into();

        match *other.inner {
            InnerExpr::Mul(m) => divisor.union(&(*m).negate()),
            _ => divisor.insert(other, -1.0),
        };

        self * Expr::from_mul_map(divisor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basics() {
        let a = Expr::var("a");

        assert_eq!(a.cpp(), "a");
        assert_eq!(Expr::var("\\aleph").cpp(), "\\aleph");
        assert_eq!(Expr::one().cpp(), "1");
        assert_eq!(Expr::zero().cpp(), "0");
    }

    #[test]
    fn commutatives() {
        let a = Expr::var("a");
        let b = Expr::var("b");

        assert_eq!((a + a).cpp(), "2 * a");
        assert_eq!((a + b).cpp(), "a + b");
        assert_eq!((a - a).cpp(), "0");
        assert_eq!((a * a).cpp(), "pow(a, 2)");
        assert_eq!((a * b).cpp(), "a * b");
        assert_eq!((a / a).cpp(), "1");
        assert_eq!((a / b).cpp(), "a / b");
    }

    #[test]
    fn constants() {
        assert_eq!((Expr::zero() + 0.0).cpp(), "0");
        assert_eq!((Expr::one() + 0.0).cpp(), "1");
        assert_eq!((Expr::zero() + 1.0).cpp(), "1");
        assert_eq!((Expr::one() + 1.0).cpp(), "2");
        assert_eq!((Expr::zero() * 0.0).cpp(), "0");
        assert_eq!((Expr::one() * 0.0).cpp(), "0");
        assert_eq!((Expr::zero() * 1.0).cpp(), "0");
        assert_eq!((Expr::one() * 1.0).cpp(), "1");
    }

    #[test]
    fn derivatives() {
        let a = Expr::var("a");
        let b = Expr::var("b");
        assert_eq!(a.deriv("a"), Expr::one());
        assert_eq!((a * a).deriv("a"), a*2);
        assert_eq!((a * a * a).deriv("a").cpp(), "3 * pow(a, 2)");
        assert_eq!(a.deriv("b").cpp(), "0");
        assert_eq!((a + b).deriv("a").cpp(), "1");
        assert_eq!(Expr::log(a).deriv("a").cpp(), "1 / a");
        assert_eq!(Expr::log(a * a).deriv("a").cpp(), "2 / a");
        assert_eq!(Expr::exp(a).deriv("a").cpp(), "exp(a)");
        assert_eq!(Expr::exp(a * a).deriv("a").cpp(), "2 * a * exp(pow(a, 2))");
        assert_eq!(Expr::pow(Expr::log(a) + a, 3).deriv("a"),
                   Expr::pow(Expr::log(a) + a, 2) * 3 * (Expr::one() / a + 1));
    }
}
