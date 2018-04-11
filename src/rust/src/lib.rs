//! This is the deft crate!

extern crate internment;
extern crate tinyset;

use tinyset::{Map64, Fits64};
use internment::Intern;

/// This is an expression.
///
/// It implements `Copy`, making it efficient to use.
#[derive(Debug, Hash, Clone, PartialEq, Eq)]
pub struct Expr<T: ExprKind + Cpp<T>> {
    inner: Intern<ExprOp<T>>,
}

impl<T: ExprKind + Cpp<T>> Copy for Expr<T> {}

impl<T: ExprKind + Cpp<T>> Fits64 for Expr<T> {
    fn to_u64(self) -> u64 {
        self.inner.to_u64()
    }
    unsafe fn from_u64(x: u64) -> Self {
        Expr { inner: Intern::from_u64(x), }
    }
}

pub trait Cpp<T: ExprKind + Cpp<T>> {
    fn cpp(op: Expr<T>, assignto: &'static str) -> String;
}

pub trait ExprKind: PartialEq + Eq + Clone + std::fmt::Debug {
    type SpecificOps: std::hash::Hash + PartialEq + Eq + Clone + std::fmt::Debug;
}

#[derive(Debug, Hash)]
enum ExprOp<T: ExprKind + Cpp<T>> {
    Add(CommutativeMap<T>),
    ScalarMul(Scalar, Expr<T>),
    Specific(T::SpecificOps),
}

#[derive(PartialEq, Eq, Clone, Debug, Hash)]
struct Scalar;
impl ExprKind for Scalar { type SpecificOps = ScalarOps; }
#[derive(Hash, PartialEq, Eq, Clone, Debug)]
enum ScalarOps {
    Mul(CommutativeMap<Scalar>),
    Exp(Expr<Scalar>),
    Log(Expr<Scalar>),
}
impl Cpp<Scalar> for Scalar {
    fn cpp(op: Expr<Scalar>, assignto: &'static str) -> String {
        match &*op.inner {
            &ExprOp::Add(ref map) => { unimplemented!() },
            &ExprOp::Specific(ref spec_op) =>
                match spec_op {
                    &ScalarOps::Mul(ref map) => { unimplemented!() },
                    otherwise => { unimplemented!() },
                },
            otherwise => { unimplemented!() }
        }
    }
}

#[derive(Debug, Clone)]
struct CommutativeMap<T: ExprKind + Cpp<T>> {
    inner: Map64<Expr<T>, f64>,
}

impl<T: ExprKind + Cpp<T>> CommutativeMap<T> {
    fn new() -> Self {
        CommutativeMap { inner: Map64::new(), }
    }

    fn as_vec(&self) -> Vec<(u64, u64)> {
        self.inner.iter().map(|(k, &v)| (k.inner.to_u64(), v as u64)).collect()
    }

    // This lets you create a `CommutativeMap` without wasting time checking if
    // the first term was already there.
    fn from_pair(k: Expr<T>, v: f64) -> Self {
        let mut map = CommutativeMap::new();
        map.inner.insert(k, v);
        map
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

    // This is used for negation, subtraction, reciprocals, division.
    fn negate(&self) -> Self {
        let mut neg = CommutativeMap::new();
        for (k, &v) in self.inner.iter() {
            neg.inner.insert(k, -v);
        }
        neg
    }
}

impl<T: ExprKind + Cpp<T>> std::hash::Hash for CommutativeMap<T> {
    fn hash<H: std::hash::Hasher>(&self, hasher: &mut H) {
        let mut terms = self.as_vec();
        terms.sort();
        for (k, v) in terms {
            k.hash(hasher);
            v.hash(hasher);
        }
    }
}

impl<T: ExprKind + Cpp<T>> PartialEq for CommutativeMap<T> {
    fn eq(&self, other: &Self) -> bool {
        let (mut lhs, mut rhs) = (self.as_vec(), other.as_vec());
        lhs.sort();
        rhs.sort();
        lhs == rhs
    }
}

impl<T: ExprKind + Cpp<T>> Eq for CommutativeMap<T> {}
