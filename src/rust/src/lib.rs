
extern crate internment;
extern crate tinyset;

use tinyset::{Map64, Fits64};
use internment::Intern;

pub trait InnerExpr: 'static + Send + Eq + std::hash::Hash {
    fn cpp(&self) -> String;
}

pub struct Expr<T: InnerExpr> {
    inner: Intern<T>,
}

impl<T: InnerExpr> Expr<T> {
    fn cpp(&self) -> String {
        self.inner.cpp()
    }
}

impl<T: InnerExpr> From<T> for Expr<T> {
    fn from(inner: T) -> Self {
        Expr { inner: Intern::new(inner) }
    }
}

#[derive(PartialEq, Eq, Hash)]
enum Scalar {
    Var(&'static str),
    Add(
        // unimplemented!()
    ),
}

impl InnerExpr for Scalar {
    fn cpp(&self) -> String { unimplemented!() }
}

impl std::ops::Add for Scalar {
    type Output = Self;

    fn add(self, other: Self) -> Self { unimplemented!() }
}

trait ClosedAdd {
    fn add(&self, other: &Self) -> Self;
}

impl<T: InnerExpr + ClosedAdd> std::ops::Add for Expr<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        ((*self.inner).add(&*other.inner)).into()
    }
}
