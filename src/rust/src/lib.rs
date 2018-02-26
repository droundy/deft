extern crate internment;
extern crate tinyset;

use tinyset::TinyMap;
use internment::Intern;
use std::cmp::{PartialEq, Eq};
use std::hash::{Hash, Hasher};

#[derive(Clone)]
struct TinyMapWrapper {
    map: TinyMap<Expr, f64>,
}

impl Hash for TinyMapWrapper {
    fn hash<H: Hasher>(&self, state: &mut H) {
    let mut terms: Vec<(*const InnerExpr, u64)>
        = self.map.iter()
                  .map(|(&t, &c)| (t.inner.as_ref() as *const InnerExpr, c as u64))
                  .collect();
    terms.sort();
        for term in terms {
            term.0.hash(state);
            term.1.hash(state);
        }
    }
}

impl PartialEq for TinyMapWrapper {
    fn eq(&self, other: &TinyMapWrapper) -> bool {
        for term in self.map.iter() {
            if other.map.get(term.0) != Some(term.1) {
                return false;
            }
        }
        true
    }
}

impl Eq for TinyMapWrapper {}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum InnerExpr {
    Var(Intern<&'static str>),
    Exp(Expr),
    Log(Expr),
    Sum(Intern<TinyMapWrapper>),
    Mul(Intern<TinyMapWrapper>),
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Expr {
    inner: Intern<InnerExpr>,
}

impl Expr {
    fn from_inner(i: InnerExpr) -> Expr {
        Expr{ inner: Intern::new(i) }
    }

    pub fn var(sym: &'static str) -> Expr {
        Expr::from_inner(InnerExpr::Var(Intern::new(sym)))
    }

    fn exp(arg: Expr) -> Expr {
        Expr::from_inner(InnerExpr::Exp(arg))
    }

    fn log(arg: Expr) -> Expr {
        Expr::from_inner(InnerExpr::Log(arg))
    }

    fn zero() -> Expr {
        Expr::from_inner(InnerExpr::Sum(Intern::new(TinyMapWrapper { map: TinyMap::new() })))
    }

    fn one() -> Expr {
        Expr::from_inner(InnerExpr::Mul(Intern::new(TinyMapWrapper { map: TinyMap::new() })))
    }

    fn cnst(x: f64) -> Expr {
        let mut term = TinyMap::new();
        term.insert(Expr::one(), x);
        Expr::from_inner(InnerExpr::Sum(Intern::new(TinyMapWrapper { map: term })))
    }

    pub fn cpp(&self) -> String {
        use InnerExpr::*;

        match self.inner.as_ref() {
            &Var(sym)
                => String::from(*sym),
            &Exp(arg)
                => String::from("exp(") + &arg.cpp() + &")",
            &Log(arg)
                => String::from("log(") + &arg.cpp() + &")",
            &Sum(wrapper) => {
                let mut terms: Vec<(String, f64)>
                    = wrapper.map.iter()
                    .map(|(s, c)| (s.cpp(), *c)).collect();

                // cannot just .sort() because f64 is not Ord
                terms.sort_by(|&(ref p, _), &(ref q, _)| p.cmp(q));

                let coeff = |t: (String, f64)| -> String {
                    if t.1.abs() == 1.0 {
                        t.0.clone()
                    } else if t.0 == "1" {
                        t.1.to_string()
                    } else {
                        t.1.abs().to_string() + &" * " + &t.0
                    }
                };

                let mut plusses: Vec<String> = Vec::new();
                let mut minuses: Vec<String> = Vec::new();

                for term in terms {
                    if term.1 > 0.0 {
                        plusses.push(coeff(term));
                    } else if term.1 < 0.0 {
                        minuses.push(coeff(term));
                    }
                    // else if term.1 == 0, discard it
                };

                if plusses.len() == 0 && minuses.len() == 0 {
                    String::from("0")
                } else {
                    plusses.join(" + ") + &minuses.iter()
                                                  .fold(String::new(), |a, s| a + &" - " + &s)
                }
            },
            &Mul(wrapper) => {
                let parens = |e: &Expr| -> String {
                    match e.inner.as_ref() {
                        &Sum(..) => String::from("(") + &e.cpp() + &")",
                        _ => e.cpp(),
                    }
                };

                let mut terms: Vec<(String, f64)>
                    = wrapper.map.iter()
                                 .map(|(s, c)| (parens(s), *c))
                                 .collect();
                terms.sort_by(|&(ref p, _), &(ref q, _)| p.cmp(q));

                let power = |t: (String, f64)| {
                    if t.1.abs() == 1.0 {
                        t.0.clone()
                    } else if t.1.abs() == 2.0 {
                        t.0.clone() + &" * " + &t.0
                    } else {
                        String::from("pow(") + &t.0 + &", " + &t.1.abs().to_string() + &")"
                    }
                };

                let mut numer: Vec<String> = Vec::new();
                let mut denom: Vec<String> = Vec::new();

                for term in terms {
                    if term.1 > 0.0 {
                        numer.push(power(term));
                    } else if term.1 < 0.0 {
                        denom.push(power(term));
                    }
                    // else if term.1 == 0, it equals one; discard
                }

                match (numer.len(), denom.len()) {
                    (0, 0) => String::from("1"),
                    (0, _) => String::from("1 / (") + &denom.join(" * ") + &")",
                    (_, 0) => numer.join(" * "),
                    (_, _) => numer.join(" * ") + &" / (" + &denom.join(" * ") + &")",
                }
            },
        }
    }
}

impl std::ops::Add for Expr {
    type Output = Expr;

    fn add(self, addend: Expr) -> Expr {
        use InnerExpr::*;

        let mut sum: TinyMap<Expr, f64> = TinyMap::new();

        match *self.inner {
            Sum(terms) => for term in terms.map.iter() {
                sum.insert(*term.0, *term.1);
            },
            _ => {
                sum.insert(self, 1.0);
            },
        };

        match *addend.inner {
            Sum(terms) => for term in terms.map.iter() {
                if sum.contains_key(term.0) {
                    let coeff = *sum.get(term.0).unwrap() + *term.1;
                    sum.insert(*term.0, coeff);
                } else {
                    sum.insert(*term.0, *term.1);
                }
            },
            _ => {
                if sum.contains_key(&addend) {
                    let coeff = *sum.get(&addend).unwrap() + 1.0;
                    sum.insert(addend, coeff);
                } else {
                    sum.insert(addend, 1.0);
                }
            },
        };

        Expr::from_inner(InnerExpr::Sum(Intern::new(TinyMapWrapper { map: sum })))
    }
}

impl std::ops::Mul for Expr {
    type Output = Expr;

    fn mul(self, multiplier: Expr) -> Expr {
        use InnerExpr::*;

        let mut mul: TinyMap<Expr, f64> = TinyMap::new();

        match *self.inner {
            Mul(terms) => for term in terms.map.iter() {
                mul.insert(*term.0, *term.1);
            },
            _ => {
                mul.insert(self, 1.0);
            },
        };

        match *multiplier.inner {
            Mul(terms) => for term in terms.map.iter() {
                if mul.contains_key(term.0) {
                    let power = *mul.get(term.0).unwrap() + *term.1;
                    mul.insert(*term.0, power);
                } else {
                    mul.insert(*term.0, *term.1);
                }
            },
            _ => {
                if mul.contains_key(&multiplier) {
                    let power = *mul.get(&multiplier).unwrap() + 1.0;
                    mul.insert(multiplier, power);
                } else {
                    mul.insert(multiplier, 1.0);
                }
            },
        };

        Expr::from_inner(InnerExpr::Mul(Intern::new(TinyMapWrapper { map: mul })))
    }
}

impl std::ops::Sub for Expr {
    type Output = Expr;

    fn sub(self, subtrahend: Expr) -> Expr {
        let mut term: TinyMap<Expr, f64> = TinyMap::new();
        term.insert(subtrahend, -1.0);
        self + Expr::from_inner(InnerExpr::Sum(Intern::new(TinyMapWrapper { map: term })))
    }
}

impl std::ops::Div for Expr {
    type Output = Expr;

    fn div(self, divisor: Expr) -> Expr {
        let mut term: TinyMap<Expr, f64> = TinyMap::new();
        term.insert(divisor, -1.0);
        self * Expr::from_inner(InnerExpr::Mul(Intern::new(TinyMapWrapper { map: term })))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn var() {
        assert_eq!(Expr::var("x").cpp(), "x");
        assert_eq!(Expr::var("\\alpha").cpp(), "\\alpha");
    }

    #[test]
    fn unaries() {
        assert_eq!(Expr::exp(Expr::var("x")).cpp(), "exp(x)");
        assert_eq!(Expr::exp(Expr::exp(Expr::var("x"))).cpp(), "exp(exp(x))");
        assert_eq!(Expr::log(Expr::var("x")).cpp(), "log(x)");
        assert_eq!(Expr::log(Expr::log(Expr::var("x"))).cpp(), "log(log(x))");
    }

    #[test]
    fn add() {
        let a = Expr::var("a");
        let b = Expr::var("b");
        let z = Expr::var("z");
        assert_eq!((a+b).cpp(), "a + b");
        assert_eq!((a+z+b).cpp(), "a + b + z");
        assert_eq!((a+a).cpp(), "2 * a");
    }

    #[test]
    fn mul() {
        let k = Expr::var("k");
        let p = Expr::var("p");
        let q = Expr::var("q");
        assert_eq!((p * q).cpp(), "p * q");
        assert_eq!((k * p * q).cpp(), "k * p * q");
    }

    #[test]
    fn mul_add() {
        let a = Expr::var("a");
        let b = Expr::var("b");
        let c = Expr::var("c");
        let d = Expr::var("d");
        let n = Expr::var("n");
        let z = Expr::var("z");
        assert_eq!((a + n * z).cpp(), "a + n * z");
        assert_eq!((a * (n + z)).cpp(), "(n + z) * a");
        assert_eq!((a * (b + c) - b * (d - z + n)).cpp(), "(b + c) * a - (d + n - z) * b");
    }

    #[test]
    fn zero() {
        assert_eq!(Expr::zero().cpp(), "0");
    }

    #[test]
    fn one() {
        assert_eq!(Expr::one().cpp(), "1");
    }

    #[test]
    fn cnst() {
        assert_eq!(Expr::cnst(6.283185).cpp(), "6.283185")
    }

    #[test]
    fn minus() {
        let a = Expr::var("a");
        let b = Expr::var("b");
        let c = Expr::var("c");
        let d = Expr::var("d");
        assert_eq!((a - a).cpp(), "0");
        assert_eq!((a - b).cpp(), "a - b");
        assert_eq!((a - d + b - c).cpp(), "a + b - c - d");
    }

    #[test]
    fn fraction() {
        let a = Expr::var("a");
        let b = Expr::var("b");
        let c = Expr::var("c");
        let d = Expr::var("d");
        assert_eq!((a / b * d / c).cpp(), "a * d / (b * c)");
        assert_eq!((a / (b + c) * d / c).cpp(), "a * d / ((b + c) * c)");
    }
}
