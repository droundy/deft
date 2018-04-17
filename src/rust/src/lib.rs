extern crate internment;
extern crate tinyset;

use tinyset::{Map64, Fits64};
use internment::Intern;

pub trait Kind: 'static + Send + Clone + Eq + std::fmt::Debug + std::hash::Hash {
    fn cpp(&self) -> String;
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Expr<T: Kind> {
    inner: Intern<T>,
}

impl<T: Kind> Expr<T> {
    fn cpp(&self) -> String {
        self.inner.cpp()
    }
}

impl<T: Kind> Copy for Expr<T> {}

impl<T: Kind> Fits64 for Expr<T> {
    fn to_u64(self) -> u64 {
        self.inner.to_u64()
    }

    unsafe fn from_u64(x: u64) -> Self {
        Expr { inner: Intern::from_u64(x), }
    }
}

impl<T: Kind> From<T> for Expr<T> {
    fn from(inner: T) -> Self {
        Expr { inner: Intern::new(inner) }
    }
}

trait ClosedAdd: Kind {
    fn add(&self, other: &Self) -> Self {
        let mut sum: AbelianMap<Self>;
        match (self.borrow_sum_map(), other.borrow_sum_map()) {
            (Some(lhs), Some(rhs)) => {
                sum = lhs.clone();
                sum.union(&rhs);
            },
            (Some(lhs), _) => {
                sum = lhs.clone();
                if *other != Self::zero() {
                    sum.insert(other.clone().into(), 1.0);
                }
            },
            (_, Some(rhs)) => {
                sum = rhs.clone();
                if *self != Self::zero() {
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
        Self::sum_from_map(sum.into())
    }

    fn neg(&self) -> Self {
        match self.borrow_sum_map() {
            Some(ref map) =>
                Self::sum_from_map(map.inner.iter().map(|(k, &v)| (k, -v)).collect()),
            _ => {
                if *self == Self::zero() {
                    Self::zero()
                } else {
                    Self::sum_from_map((self.clone().into(), -1.0).into())
                }
            },
        }
    }

    fn zero() -> Self {
        Self::sum_from_map(AbelianMap::new())
    }

    fn sum_from_map(x: AbelianMap<Self>) -> Self;
    fn borrow_sum_map(&self) -> Option<&AbelianMap<Self>>;
}

trait ClosedMul {
    fn mul(&self, other: &Self) -> Self;
    fn one() -> Self;
}

trait ClosedArithmetic: ClosedAdd + ClosedMul {
    fn reciprocal(&self) -> Self;
}

impl<T: Kind + ClosedAdd> std::ops::Add for Expr<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        (*self.inner).add(&*other.inner).into()
    }
}

impl<T: Kind + ClosedAdd> std::ops::Sub for Expr<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        (*self.inner).add(&other.inner.neg()).into()
    }
}

impl<T: Kind + ClosedMul> std::ops::Mul for Expr<T> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        (*self.inner).mul(&*other.inner).into()
    }
}

impl<T: Kind + ClosedArithmetic> std::ops::Div for Expr<T> {
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

impl Kind for Scalar {
    fn cpp(&self) -> String { unimplemented!() }
}

impl ClosedAdd for Scalar {
    fn sum_from_map(m: AbelianMap<Self>) -> Self {
        Scalar::Add(m)
    }
    fn borrow_sum_map(&self) -> Option<&AbelianMap<Self>> {
        if let &Scalar::Add(ref m) = self {
            Some(m)
        } else {
            None
        }
    }
}

#[derive(PartialEq, Eq, Hash, Clone, Debug)]
enum RealSpaceScalar {
    Var(&'static str),
    Exp(Expr<RealSpaceScalar>),
    Log(Expr<RealSpaceScalar>),
    Add(AbelianMap<RealSpaceScalar>),
    Mul(AbelianMap<RealSpaceScalar>),
}

impl Kind for RealSpaceScalar {
    fn cpp(&self) -> String { unimplemented!() }
}

impl ClosedAdd for RealSpaceScalar {
    fn sum_from_map(m: AbelianMap<Self>) -> Self {
        RealSpaceScalar::Add(m)
    }
    fn borrow_sum_map(&self) -> Option<&AbelianMap<Self>> {
        if let &RealSpaceScalar::Add(ref m) = self {
            Some(m)
        } else {
            None
        }
    }
}

#[derive(Debug, Clone)]
struct AbelianMap<T: Kind> {
    inner: Map64<Expr<T>, f64>,
}

impl<T: Kind> AbelianMap<T> {
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

impl<T: Kind> std::iter::FromIterator<(Expr<T>, f64)> for AbelianMap<T> {
    fn from_iter<I: IntoIterator<Item = (Expr<T>, f64)>>(iter: I) -> Self {
        let mut map = AbelianMap::new();
        for (k, v) in iter {
            map.insert(k, v);
        }
        map
    }
}

impl<T: Kind> From<(Expr<T>, f64)> for AbelianMap<T> {
    fn from((k, v): (Expr<T>, f64)) -> Self {
        let mut inner = Map64::new();
        inner.insert(k, v);
        Self { inner, }
    }
}

impl<T: Kind> From<Map64<Expr<T>, f64>> for AbelianMap<T> {
    fn from(inner: Map64<Expr<T>, f64>) -> Self {
        Self { inner, }
    }
}

impl<T: Kind> std::hash::Hash for AbelianMap<T> {
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

impl<T: Kind> PartialEq for AbelianMap<T> {
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

impl<T: Kind> Eq for AbelianMap<T> {}
