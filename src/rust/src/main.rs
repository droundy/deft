extern crate internment;
extern crate tinyset;

use tinyset::TinyMap;
use internment::Intern;
use std::cmp::{PartialEq, Eq};
use std::hash::{Hash, Hasher};
use std::mem::discriminant;

#[derive(Clone)]
struct TinyMapWrapper {
    map: TinyMap<Intern<Expr>, f64>,
}

impl Hash for TinyMapWrapper {
    fn hash<H: Hasher>(&self, state: &mut H) {
        for term in self.map.iter() {
            (*term.0).hash(state);
            (*term.1 as u64).hash(state);
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

#[derive(Clone, Copy)]
enum Expr {
    Var(Intern<&'static str>),
    Exp(Intern<Expr>),
    Log(Intern<Expr>),
    Sum(Intern<TinyMapWrapper>),
    Mul(Intern<TinyMapWrapper>),
}

impl PartialEq for Expr {
    fn eq(&self, other: &Expr) -> bool {
        use Expr::*;

        match self {
            &Var(sym_l) => match other { &Var(sym_r) => sym_l == sym_r, _ => false, },
            &Exp(arg_l) => match other { &Exp(arg_r) => arg_l == arg_r, _ => false, },
            &Log(arg_l) => match other { &Log(arg_r) => arg_l == arg_r, _ => false, },
            &Sum(map_l) => match other { &Sum(map_r) => map_l == map_r, _ => false, },
            &Mul(map_l) => match other { &Mul(map_r) => map_l == map_r, _ => false, },
        }
    }
}

impl Eq for Expr {}

impl Hash for Expr {
    fn hash<H: Hasher>(&self, state: &mut H) {
        use Expr::*;

        discriminant(self).hash(state);

        match self {
            &Var(ref sym)
                => sym.hash(state),
            &Exp(ref arg) | &Log(ref arg)
                => arg.hash(state),
            &Sum(ref map) | &Mul(ref map)
                => map.hash(state),
        };
    }
}

impl Expr {
    fn var(sym: &'static str) -> Intern<Expr> {
        Intern::new(Expr::Var(Intern::new(sym)))
    }

    fn exp(arg: Intern<Expr>) -> Intern<Expr> {
        Intern::new(Expr::Exp(arg))
    }

    fn log(arg: Intern<Expr>) -> Intern<Expr> {
        Intern::new(Expr::Log(arg))
    }

    fn sum(augend: Intern<Expr>, addend: Intern<Expr>) -> Intern<Expr> {
        use Expr::*;

        let mut sum: TinyMap<Intern<Expr>, f64> = TinyMap::new();

        match *augend {
            Sum(terms) => for term in terms.map.iter() {
                sum.insert(*term.0, *term.1);
            },
            _ => {
                sum.insert(augend, 1.0);
            },
        };

        match *addend {
            Sum(terms) => for term in terms.map.iter() {
                if sum.contains_key(term.0) {
                    let coeff = *sum.get(term.0).unwrap() + 1.0;
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

        Intern::new(Expr::Sum(Intern::new(TinyMapWrapper { map: sum })))
    }

    fn mul(multiplicand: Intern<Expr>, multiplier: Intern<Expr>) -> Intern<Expr> {
        use Expr::*;

        let mut mul: TinyMap<Intern<Expr>, f64> = TinyMap::new();

        match *multiplicand {
            Mul(terms) => for term in terms.map.iter() {
                mul.insert(*term.0, *term.1);
            },
            _ => {
                mul.insert(multiplicand, 1.0);
            },
        };

        match *multiplier {
            Mul(terms) => for term in terms.map.iter() {
                if mul.contains_key(term.0) {
                    let power = *mul.get(term.0).unwrap() + 1.0;
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

        Intern::new(Expr::Mul(Intern::new(TinyMapWrapper { map: mul })))
    }

    fn cpp(&self) -> String {
        use Expr::*;

        match self {
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
                terms.sort_by(|&(ref p, _), &(ref q, _)| p.cmp(q));

                let coeff = |t: &(String, f64)| {
                    if t.1 == 1.0 {
                        t.0.clone()
                    } else {
                        t.1.to_string() + &" * " + &t.0
                    }
                };

                match terms.first() {
                    None => String::from("0"),
                    Some(first) => String::from("(") + &terms.iter().skip(1)
                        .fold(coeff(first), |acc, term| acc + &" + " + &coeff(term)) + &")"
                }
            },
            &Mul(wrapper) => {
                let mut terms: Vec<(String, f64)>
                    = wrapper.map.iter()
                    .map(|(s, c)| (s.cpp(), *c)).collect();
                terms.sort_by(|&(ref p, _), &(ref q, _)| p.cmp(q));

                let power = |t: &(String, f64)| {
                    if t.1 == 1.0 {
                        t.0.clone()
                    } else if t.1 == 2.0 {
                        t.0.clone() + &" * " + &t.0
                    } else {
                        String::from("pow(") + &t.0 + &", " + &t.1.to_string() + &")"
                    }
                };

                match terms.first() {
                    None => String::from("1"),
                    Some(first) => String::from("(") + &terms.iter().skip(1)
                        .fold(power(first), |acc, term| acc + &" * " + &power(term)) + &")"
                }
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn var() {
        assert_eq!(Expr::var("x").cpp(), String::from("x"));
        assert_eq!(Expr::var("\\alpha").cpp(), String::from("\\alpha"));
    }

    #[test]
    fn unaries() {
        assert_eq!(Expr::exp(Expr::var("x")).cpp(), String::from("exp(x)"));
        assert_eq!(Expr::exp(Expr::exp(Expr::var("x"))).cpp(), String::from("exp(exp(x))"));
        assert_eq!(Expr::log(Expr::var("x")).cpp(), String::from("log(x)"));
        assert_eq!(Expr::log(Expr::log(Expr::var("x"))).cpp(), String::from("log(log(x))"));
    }

    #[test]
    fn sum() {
        assert_eq!(Expr::sum(Expr::var("a"), Expr::var("b")).cpp(), String::from("(a + b)"));
        assert_eq!(Expr::sum(Expr::sum(Expr::var("a"), Expr::var("z")), Expr::var("b")).cpp(),
                   String::from("(a + b + z)"));
        assert_eq!(Expr::sum(Expr::var("a"), Expr::var("a")).cpp(), String::from("(2 * a)"));
    }

    #[test]
    fn mul() {
        assert_eq!(Expr::mul(Expr::var("p"), Expr::var("q")).cpp(), String::from("(p * q)"));
        assert_eq!(Expr::mul(Expr::mul(Expr::var("p"), Expr::var("k")), Expr::var("q")).cpp(),
                   String::from("(k * p * q)"));
    }
}
