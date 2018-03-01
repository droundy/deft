//! This is the deft crate!

extern crate internment;
extern crate tinyset;

use tinyset::{Map64, Fits64};
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
}

impl Fits64 for Expr {
    fn to_u64(self) -> u64 {
        self.inner.to_u64()
    }

    unsafe fn from_u64(x: u64) -> Self {
        Expr { inner: Intern::<InnerExpr>::from_u64(x) }
    }
}

impl Expr {
    fn from_inner(i: InnerExpr) -> Expr {
        Expr { inner: Intern::new(i) }
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
                    if c.abs() == 1.0 {
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
                let n = n.to_string(" * ", &power);
                match d.len() {
                    0 => n,
                    1 => n + &" / " + &d.to_string(" * ", &power),
                    _ => n + &" / (" + &d.to_string(" * ", &power) + &")",
                }
            },
        }
    }
}

impl std::ops::Add for Expr {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let mut sum = CommutativeMap::new();

        match *self.inner {
            InnerExpr::Sum(m) => sum.union(&*m),
            _ => sum.insert(self, 1.0),
        };

        match *other.inner {
            InnerExpr::Sum(m) => sum.union(&*m),
            _ => sum.insert(other, 1.0),
        };

        Expr::from_sum_map(sum)
    }
}

impl std::ops::Mul for Expr {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let mut mul = CommutativeMap::new();

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

impl std::ops::Sub for Expr {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let mut subtrahend = CommutativeMap::new();

        match *other.inner {
            InnerExpr::Sum(m) => subtrahend.union(&(*m).negate()),
            _ => subtrahend.insert(other, -1.0),
        };

        self + Expr::from_sum_map(subtrahend)
    }
}

impl std::ops::Div for Expr {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        let mut divisor = CommutativeMap::new();

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
}
