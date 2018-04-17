extern crate internment;
extern crate tinyset;

use tinyset::{Map64, Fits64};
use internment::Intern;

pub trait InnerExpr: 'static + Send + Clone + Eq + std::fmt::Debug + std::hash::Hash {
    fn cpp(&self) -> String;
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Expr<T: InnerExpr> {
    inner: Intern<T>,
}

impl<T: InnerExpr> Expr<T> {
    fn cpp(&self) -> String {
        self.inner.cpp()
    }
}

impl<T: InnerExpr> Copy for Expr<T> {}

impl<T: InnerExpr> Fits64 for Expr<T> {
    fn to_u64(self) -> u64 {
        self.inner.to_u64()
    }

    unsafe fn from_u64(x: u64) -> Self {
        Expr { inner: Intern::from_u64(x), }
    }
}

impl<T: InnerExpr> From<T> for Expr<T> {
    fn from(inner: T) -> Self {
        Expr { inner: Intern::new(inner) }
    }
}

trait ExprAdd: InnerExpr {
    fn add(&self, other: &Self) -> Self;
    fn neg(&self) -> Self;
    fn zero() -> Self;
}

trait ExprMul {
    fn mul(&self, other: &Self) -> Self;
    fn one() -> Self;
}

trait ExprField: ExprAdd + ExprMul {
    fn reciprocal(&self) -> Self;
}

impl<T: InnerExpr + ExprAdd> std::ops::Add for Expr<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        (*self.inner).add(&*other.inner).into()
    }
}

impl<T: InnerExpr + ExprAdd> std::ops::Sub for Expr<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        (*self.inner).add(&other.inner.neg()).into()
    }
}

impl<T: InnerExpr + ExprMul> std::ops::Mul for Expr<T> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        (*self.inner).mul(&*other.inner).into()
    }
}

impl<T: InnerExpr + ExprField> std::ops::Div for Expr<T> {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        (*self.inner).mul(&other.inner.reciprocal()).into()
    }
}

#[derive(PartialEq, Eq, Hash, Clone, Debug)]
enum Scalar {
    Var(&'static str),
    Exp(Expr<Scalar>),
    Log(Expr<Scalar>),
    Add(AbelianMap<Scalar>),
    Mul(AbelianMap<Scalar>),
}

impl InnerExpr for Scalar {
    fn cpp(&self) -> String { unimplemented!() }
}

impl ExprAdd for Scalar {
    fn add(&self, other: &Self) -> Self {
        let mut sum: AbelianMap<Scalar>;
        match (self, other) {
            (&Scalar::Add(ref lhs), &Scalar::Add(ref rhs)) => {
                sum = lhs.clone();
                sum.union(&rhs);
            },
            (&Scalar::Add(ref lhs), _) => {
                sum = lhs.clone();
                if *other != Scalar::zero() {
                    sum.insert(other.clone().into(), 1.0);
                }
            },
            (_, &Scalar::Add(ref rhs)) => {
                sum = rhs.clone();
                if *self != Scalar::zero() {
                    sum.insert(self.clone().into(), 1.0);
                }
            },
            (_, _) => {
                sum = (self.clone().into(), 1.0).into();
                sum.insert(other.clone().into(), 1.0);
            }
        }
        if sum.inner.len() == 1 {
            let (k, &v) = sum.inner.iter().next().unwrap();
            if v == 1.0 {
                return (*k.inner).clone();
            }
        }
        Scalar::Add(sum.into())
    }

    fn neg(&self) -> Self {
        match self {
            &Scalar::Add(ref map) =>
                Scalar::Add(map.inner.iter().map(|(k, &v)| (k, -v)).collect()),
            _ => {
                if *self == Scalar::zero() {
                    Scalar::zero()
                } else {
                    Scalar::Add((self.clone().into(), -1.0).into())
                }
            },
        }
    }

    fn zero() -> Self {
        Scalar::Mul(AbelianMap::new())
    }
}

#[derive(Debug, Clone)]
struct AbelianMap<T: InnerExpr> {
    inner: Map64<Expr<T>, f64>,
}

impl<T: InnerExpr> AbelianMap<T> {
    fn new() -> Self {
        Self { inner: Map64::new(), }
    }

    fn insert(&mut self, k: Expr<T>, v: f64) {
        let v = self.inner.get(&k).unwrap_or(&0.0) + v;
        if v == 0.0 {
            self.inner.remove(&k);
        } else {
            self.inner.insert(k, v);
        }
    }

    fn union(&mut self, other: &Self) {
        for (k, &v) in other.inner.iter() {
            self.insert(k, v);
        }
    }
}

impl<T: InnerExpr> std::iter::FromIterator<(Expr<T>, f64)> for AbelianMap<T> {
    fn from_iter<I: IntoIterator<Item = (Expr<T>, f64)>>(iter: I) -> Self {
        let mut map = AbelianMap::new();
        for (k, v) in iter {
            map.insert(k, v);
        }
        map
    }
}

impl<T: InnerExpr> From<(Expr<T>, f64)> for AbelianMap<T> {
    fn from((k, v): (Expr<T>, f64)) -> Self {
        let mut inner = Map64::new();
        inner.insert(k, v);
        Self { inner, }
    }
}

impl<T: InnerExpr> From<Map64<Expr<T>, f64>> for AbelianMap<T> {
    fn from(inner: Map64<Expr<T>, f64>) -> Self {
        Self { inner, }
    }
}

impl<T: InnerExpr> std::hash::Hash for AbelianMap<T> {
    fn hash<H: std::hash::Hasher>(&self, hasher: &mut H) {
        let mut terms: Vec<_>
            = self.inner.iter()
                        .map(|(k, &v)| (k.inner.to_u64(), v as u64))
                        .collect();
        terms.sort();
        for (k, v) in terms {
            k.hash(hasher);
            v.hash(hasher);
        }
    }
}

impl<T: InnerExpr> PartialEq for AbelianMap<T> {
    fn eq(&self, other: &Self) -> bool {
        let mut lhs: Vec<_>
            = self.inner.iter()
                        .map(|(k, &v)| (k.inner.to_u64(), v as u64))
                        .collect();
        let mut rhs: Vec<_>
            = other.inner.iter()
                         .map(|(k, &v)| (k.inner.to_u64(), v as u64))
                         .collect();
        lhs.sort();
        rhs.sort();
        lhs == rhs
    }
}

impl<T: InnerExpr> Eq for AbelianMap<T> {}
