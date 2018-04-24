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
    fn sum_from_map(x: AbelianMap<Self>) -> Self;
    fn borrow_sum_map(&self) -> Option<&AbelianMap<Self>>;

    fn add(&self, other: &Self) -> Self {
        let mut sum: AbelianMap<Self>;
        match (self.borrow_sum_map(), other.borrow_sum_map()) {
            (Some(lhs), Some(rhs)) => {
                sum = lhs.clone();
                sum.union(&rhs);
            },
            (Some(lhs), _) => {
                sum = lhs.clone();
                sum.insert(other.clone().into(), 1.0);
            },
            (_, Some(rhs)) => {
                sum = rhs.clone();
                sum.insert(self.clone().into(), 1.0);
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

}

trait ClosedMul: Kind {
    fn product_from_map(x: AbelianMap<Self>) -> Self;
    fn borrow_product_map(&self) -> Option<&AbelianMap<Self>>;

    fn mul(&self, other: &Self) -> Self {
        let mut product: AbelianMap<Self>;
        match (self.borrow_product_map(), other.borrow_product_map()) {
            (Some(lhs), Some(rhs)) => {
                product = lhs.clone();
                product.union(&rhs);
            },
            (Some(lhs), _) => {
                product = lhs.clone();
                product.insert(other.clone().into(), 1.0);
            },
            (_, Some(rhs)) => {
                product = rhs.clone();
                product.insert(self.clone().into(), 1.0);
            },
            (_, _) => {
                product = (self.clone().into(), 1.0).into();
                product.insert(other.clone().into(), 1.0);
            },
        }
        if product.inner.len() == 1 {
            let (k, &v) = product.inner.iter().next().unwrap();
            if v == 1.0 {
                return (*k.inner).clone();
            }
        }
        Self::product_from_map(product)
    }

    fn one() -> Self {
        Self::product_from_map(AbelianMap::new())
    }
}

trait ClosedArithmetic: ClosedAdd + ClosedMul {
    fn reciprocal(&self) -> Self {
        match self.borrow_product_map() {
            Some(ref map) =>
                Self::product_from_map(map.inner.iter().map(|(k, &v)| (k, -v)).collect()),
            _ => {
                if *self == Self::zero() {
                    Self::zero()
                } else {
                    Self::product_from_map((self.clone().into(), -1.0).into())
                }
            },
        }
    }
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
    fn cpp(&self) -> String {
        match self {
            &Scalar::Var(sym) => String::from(sym),
            &Scalar::Exp(arg) => String::from("exp(") + &arg.cpp() + &")",
            &Scalar::Log(arg) => String::from("log(") + &arg.cpp() + &")",
            &Scalar::Add(ref map) => {
                let pos_coeff = |&(ref x, ref c): &(String, f64)| -> String {
                    if x == "1" {
                        c.to_string()
                    } else if *c == 1.0 {
                        x.clone()
                    } else {
                        c.to_string() + &" * " + &x
                    }
                };
                let neg_coeff = |&(ref x, ref c): &(String, f64)| -> String {
                    if x == "1" {
                        c.abs().to_string()
                    } else if *c == -1.0 {
                        x.clone()
                    } else {
                        c.abs().to_string() + &" * " + &x
                    }
                };
                let (p, n) = map.split_cpp_sort();
                match (p.len(), n.len()) {
                    (0, 0) => String::from("0"),
                    (_, 0) => p.iter()
                               .map(pos_coeff)
                               .collect::<Vec<String>>()
                               .join(" + "),
                    (0, _) => String::from("-")
                              + &n.iter()
                                  .map(neg_coeff)
                                  .collect::<Vec<String>>()
                                  .join(" - "),
                    (_, _) => p.iter()
                               .map(pos_coeff)
                               .collect::<Vec<String>>()
                               .join(" + ")
                              + &" - "
                              + &n.iter()
                                  .map(neg_coeff)
                                  .collect::<Vec<String>>()
                                  .join(" - "),
                }
            },
            &Scalar::Mul(ref map) => {
                let ref power = |&(ref x, ref p): &(String, f64)| -> String {
                    if x == "1" || p.abs() == 1.0 {
                        x.clone()
                    } else if p.abs() == 2.0 {
                        x.clone() + &" * " + &x
                    } else {
                        String::from("pow(") + &x + &", " + &p.abs().to_string() + &")"
                    }
                };
                let (n, d) = map.split_cpp_sort();
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
                               .map(power)
                               .collect::<Vec<String>>()
                               .join(" * ")
                              + &" / ("
                              + &d.iter()
                                  .map(power)
                                  .collect::<Vec<String>>()
                                  .join(" * ")
                              + &")",
                }
            },
        }
    }
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

impl ClosedMul for Scalar {
    fn product_from_map(m: AbelianMap<Self>) -> Self {
        Scalar::Mul(m)
    }

    fn borrow_product_map(&self) -> Option<&AbelianMap<Self>> {
        if let &Scalar::Mul(ref m) = self {
            Some(m)
        } else {
            None
        }
    }
}

impl ClosedArithmetic for Scalar {}

#[derive(PartialEq, Eq, Hash, Clone, Debug)]
enum RealSpaceScalar {
    Var(&'static str),
    ScalarVar(&'static str),
    Exp(Expr<RealSpaceScalar>),
    Log(Expr<RealSpaceScalar>),
    Add(AbelianMap<RealSpaceScalar>),
    Mul(AbelianMap<RealSpaceScalar>),
    FFT(Expr<RealSpaceScalar>),
}

impl From<Scalar> for RealSpaceScalar {
    fn from(s: Scalar) -> Self {
        match s {
            Scalar::Var(sym) => RealSpaceScalar::ScalarVar(sym),
            _ => panic!(),
        }
    }
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

    fn split_cpp_sort(&self) -> (Vec<(String, f64)>, Vec<(String, f64)>) {
        let (mut p, mut n): (Vec<(String, f64)>, Vec<(String, f64)>)
            = self.inner
                  .iter()
                  .map(|(k, &v)| (k.cpp(), v))
                  .partition(|&(_, v)| v > 0.0);
        let ref by_key = |&(ref p, _): &(String, f64), &(ref q, _): &(String, f64)| {
            p.cmp(&q)
        };
        p.sort_by(by_key);
        n.sort_by(by_key);
        (p, n)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scalar() {
        let a: Expr<Scalar> = Scalar::Var("a").into();
        let b: Expr<Scalar> = Scalar::Var("b").into();

        let zero: Expr<Scalar> = Scalar::zero().into();
        let one: Expr<Scalar> = Scalar::zero().into();

        assert_eq!(a + b, a + b);
        assert_eq!(a + b, b + a);
        assert!(a + a != b + b);
        assert_eq!(a + zero, a);

        if let Scalar::Mul(ref map) = *(a * one).inner {
            for (k, &v) in map.inner.iter() {
                println!("{} {}", k.cpp(), v);
            }
        }

        assert_eq!(a * one, a);

        assert_eq!(a - a, zero);
        assert_eq!(a + a - a, a);
        assert_eq!(a + a - a - a, zero);
        assert_eq!(a + b - a - b, zero);
        assert_eq!(a + b - b - a, zero);
        assert!(b - a != a - b);

        assert_eq!(a.cpp(), "a");
        assert!(b.cpp() != "a");
        assert_eq!((a + b).cpp(), "a + b");
        assert_eq!((b + a).cpp(), "a + b");
        assert_eq!((a + a).cpp(), "2 * a");
        assert_eq!((a + a - a).cpp(), a.cpp());
    }
}
