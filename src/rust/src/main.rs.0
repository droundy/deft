extern crate internment;
extern crate tinyset;

use tinyset::{TinyMap};
use internment::Intern;
use std::cmp::{PartialEq, Eq};
use std::hash::{Hash, Hasher};

#[derive(Clone)]
struct TemporaryTinyMapWrapper {
    tm: TinyMap<Intern<Expr>, f64>,
}

impl Hash for TemporaryTinyMapWrapper {
    fn hash<H: Hasher>(&self, state: &mut H) {
        for term in self.tm.iter() {
            (*term.0).hash(state);
            (*term.1 as u64).hash(state);
        }
    }
}

impl PartialEq for TemporaryTinyMapWrapper {
    fn eq(&self, other: &TemporaryTinyMapWrapper) -> bool {
        for term in self.tm.iter() {
            if other.tm.get(term.0) != Some(term.1) {
                return false;
            }
        }
        true
    }
}

impl Eq for TemporaryTinyMapWrapper {}

#[derive(Clone, Copy)]
enum Expr {
    Var(Intern<&'static str>),
    Exp(Intern<Expr>),
    Log(Intern<Expr>),
    Sum(Intern<TemporaryTinyMapWrapper>),
    Prod(Intern<TemporaryTinyMapWrapper>),
}

impl PartialEq for Expr {
    fn eq(&self, other: &Expr) -> bool {
        use Expr::*;

        match self {
            &Var(sym1) => match other {
                &Var(sym2) => sym1 == sym2,
                _ => false,
            },
            &Exp(arg1) => match other {
                &Exp(arg2) => arg1 == arg2,
                _ => false,
            },
            &Log(arg1) => match other {
                &Log(arg2) => arg1 == arg2,
                _ => false,
            },
            &Sum(terms1) => match other {
                &Sum(terms2) => terms1 == terms2,
                _ => false,
            },
            &Prod(terms1) => match other {
                &Prod(terms2) => terms1 == terms2,
                _ => false,
            },
        }
    }
}

impl Eq for Expr {}

impl Hash for Expr {
    fn hash<H: Hasher>(&self, state: &mut H) {
        use Expr::*;

        match self {
            &Var(..)  => 1.hash(state),
            &Exp(..)  => 2.hash(state),
            &Log(..)  => 3.hash(state),
            &Sum(..)  => 4.hash(state),
            &Prod(..) => 5.hash(state),
        }

        match self {
            &Var(ref sym) => sym.hash(state),
            &Exp(ref arg) | &Log(ref arg) => arg.hash(state),
            &Sum(ref map) | &Prod(ref map) => map.hash(state),
        }
    }
}

impl Expr {
    fn cpp(&self) -> String {
        use Expr::*;

        match self {
            &Var(sym) => String::from(*sym),
            &Exp(arg) => String::from("exp(") + &arg.cpp() + &")",
            &Log(arg) => String::from("log(") + &arg.cpp() + &")",
            &Sum(term_map) => {
                let mut term_vec: Vec<(String, f64)>
                    = term_map.tm.iter()
                              .map(|(s, c)| (s.cpp(), *c))
                              .collect();
                term_vec.sort_by(|&(ref s1, _), &(ref s2, _)| s1.cmp(s2));

                let coeff = |t: &(String, f64)| {
                    if t.1 == 1.0 {
                        t.0.clone()
                    } else {
                        t.1.to_string() + &" * " + &t.0
                    }
                };

                match term_vec.first() {
                    None => String::from("0"),
                    Some(term1) => term_vec.iter().skip(1).fold(coeff(term1), |acc, term| acc + &" + " + &coeff(term))
                }
            },
            &Prod(term_map) => {
                let mut term_vec: Vec<(String, f64)>
                    = term_map.tm.iter()
                              .map(|(s, c)| (s.cpp(), *c))
                              .collect();
                term_vec.sort_by(|&(ref s1, _), &(ref s2, _)| s1.cmp(s2));
                term_vec.iter()
                         .fold(String::new(), |acc, &(ref s, ref c)|
                                   acc + &" * pow(" + &s + &", " + &c.to_string() + &")")
            },
        }
    }

    fn var(sym: &'static str) -> Intern<Expr> {
        Intern::new(Expr::Var(Intern::new(sym)))
    }

    fn sum(augend: Intern<Expr>, addend: Intern<Expr>) -> Intern<Expr> {
        let mut sum: TinyMap<Intern<Expr>, f64> = TinyMap::new();

        match *augend {
            Expr::Var(..) | Expr::Exp(..) | Expr::Log(..) | Expr::Prod(..) => {
                sum.insert(augend, 1.0);
            },
            Expr::Sum(terms) => for term in terms.tm.iter() {
                sum.insert(*term.0, *term.1);
            },
        };

        match *addend {
            Expr::Var(..) | Expr::Exp(..) | Expr::Log(..) | Expr::Prod(..) => {
                if sum.contains_key(&addend) {
                    let coeff = *sum.get(&addend).unwrap() + 1.0;
                    sum.insert(addend, coeff);
                } else {
                    sum.insert(addend, 1.0);
                }
            },
            Expr::Sum(terms) => for term in terms.tm.iter() {
                if sum.contains_key(term.0) {
                    let coeff = *sum.get(term.0).unwrap() + 1.0;
                    sum.insert(*term.0, coeff);
                } else {
                    sum.insert(*term.0, *term.1);
                }
            },
        };

        Intern::new(Expr::Sum(Intern::new(TemporaryTinyMapWrapper { tm: sum })))
    }
}

#[cfg(test)]
mod tests {
}

fn main() {
    println!("{}", Expr::sum(Expr::sum(Intern::new(Expr::Log(Expr::var("z"))), Expr::var("k")), Expr::var("m")).cpp());
}
