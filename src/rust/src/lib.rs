extern crate internment;
extern crate tinyset;

use tinyset::{Map64, Fits64};
use internment::Intern;

pub enum SingleStatement<T: Kind> {
    Init(&'static str, Expr<T>),
    Assign(&'static str, Expr<T>),
    Free(&'static str),
}

pub enum Statement {
    Block(Vec<Statement>),
    S(SingleStatement<Scalar>),
    RSS(SingleStatement<RealSpaceScalar>),
    InitScalar(&'static str, Expr<Scalar>),
    InitRealSpaceScalar(&'static str, Expr<RealSpaceScalar>),
    AssignRealSpaceScalar(&'static str, Expr<RealSpaceScalar>),
    FreeRealSpaceScalar(&'static str),
}

pub trait Kind: 'static + Send + Clone + Eq + std::fmt::Debug + std::hash::Hash {
    fn cpp(&self) -> String;
}

#[derive(Debug, Clone, Eq, Hash)]
pub struct Expr<T: Kind> {
    inner: Intern<T>,
}

impl<T: Kind> Expr<T> {
    fn cpp(&self) -> String {
        self.inner.cpp()
    }
}

impl<T: Kind, RHS: Clone + Into<Expr<T>>> PartialEq<RHS> for Expr<T> {
    fn eq(&self, other: &RHS) -> bool {
        self.inner == other.clone().into().inner
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
        Self::sum_from_map(sum)
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

    fn mul<T: ClosedMul + Into<Self>>(&self, other: &T) -> Self {
        let mut product: AbelianMap<Self>;
        let other: &Self = &other.clone().into();
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

    fn from_f64<F: Into<f64>>(f: F) -> Self {
        let f = f.into();
        if f == 0.0 {
            Self::zero()
        } else if f == 1.0 {
            Self::one()
        } else {
            Self::sum_from_map((Expr::from(Self::one()), f).into())
        }
    }
}

impl<F: Into<f64>> From<F> for Expr<Scalar> {
    fn from(f: F) -> Self {
        Scalar::from_f64(f).into()
    }
}

// impl<T: Kind, U: Kind + From<T>> From<Expr<T>> for Expr<U> {
//     fn from(f: Expr<T>) -> Self {
//         f.into().into()
//     }
// }

impl<T: Kind + ClosedAdd, RHS: Into<Expr<T>>> std::ops::Add<RHS> for Expr<T> {
    type Output = Self;

    fn add(self, other: RHS) -> Self {
        let other = other.into();
        (*self.inner).add(&*other.inner).into()
    }
}

impl<T: Kind + ClosedAdd, RHS: Into<Expr<T>>> std::ops::Sub<RHS> for Expr<T> {
    type Output = Self;

    fn sub(self, other: RHS) -> Self {
        let other = other.into();
        (*self.inner).add(&other.inner.neg()).into()
    }
}

impl<T: Kind + ClosedMul, RHS: Into<Expr<T>>> std::ops::Mul<RHS> for Expr<T> {
    type Output = Self;

    fn mul(self, other: RHS) -> Self {
        let other = other.into();
        (*self.inner).mul(&*other.inner).into()
    }
}

impl<T: Kind + ClosedArithmetic, RHS: Into<Expr<T>>> std::ops::Div<RHS> for Expr<T> {
    type Output = Self;

    fn div(self, other: RHS) -> Self {
        let other = other.into();
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
        if m.inner.len() == 1 {
            let (k, &v) = m.inner.iter().next().unwrap();
            if v == 1.0 {
                return (*k.inner).clone();
            }
        }
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
            Scalar::Exp(arg) => RealSpaceScalar::Exp(RealSpaceScalar::from((*arg.inner).clone()).into()),
            Scalar::Log(arg) => RealSpaceScalar::Log(RealSpaceScalar::from((*arg.inner).clone()).into()),
            Scalar::Add(map) =>
                RealSpaceScalar::sum_from_map(map.inner.iter()
                                              .map(|(k,&v)| (k.into(), v)).collect()),
            _ => panic!(),
        }
    }
}

// impl Into<Scalar> for RealSpaceScalar {
//     fn into(s: Self) -> Scalar {
// 
//     }
// }

impl Kind for RealSpaceScalar {
    fn cpp(&self) -> String {
        unimplemented!()
    }
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

impl ClosedMul for RealSpaceScalar {
    fn product_from_map(m: AbelianMap<Self>) -> Self {
        RealSpaceScalar::Mul(m)
    }

    fn borrow_product_map(&self) -> Option<&AbelianMap<Self>> {
        if let &RealSpaceScalar::Mul(ref m) = self {
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

// impl<F: Kind, T: Kind + From<F>> From<AbelianMap<F>> for AbelianMap<T> {
//     fn from(f: AbelianMap<F>) -> Self {
//          <AbelianMap<T>>::from(f.inner.iter().map(|(k, &v)| (k.into(), v)).collect());
//     }
// }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scalar() {
        let a: Expr<Scalar> = Scalar::Var("a").into();
        let b: Expr<Scalar> = Scalar::Var("b").into();

        assert_eq!(a + b, a + b);
        assert_eq!(a + b, b + a);
        assert!(a + a != b + b);
        assert_eq!(a + 0, a);

        assert_eq!(a * 1, a);
        assert_eq!(a * 1, a);

        assert_eq!(a - a, 0);
        assert_eq!(a + a - a, a);
        assert_eq!(a + a - a - a, 0);
        assert_eq!(a + b - a - b, 0);
        assert_eq!(a + b - b - a, 0);
        assert!(b - a != a - b);

        assert_eq!(a.cpp(), "a");
        assert!(b.cpp() != "a");
        assert_eq!((a + b).cpp(), "a + b");
        assert_eq!((b + a).cpp(), "a + b");
        assert_eq!((a + a).cpp(), "2 * a");
        assert_eq!((a + a - a).cpp(), a.cpp());
    }

    #[test]
    fn realspace() {
        let a: Expr<RealSpaceScalar> = RealSpaceScalar::Var("a").into();
        let s: Expr<Scalar> = <Expr<Scalar>>::from(Scalar::Var("s"));

        assert_eq!(a * s, a * <Expr<RealSpaceScalar>>::from(RealSpaceScalar::ScalarVar("s")));
    }
}

impl From<Expr<Scalar>> for Expr<RealSpaceScalar> {
    fn from(f: Expr<Scalar>) -> Self {
        let f: RealSpaceScalar = (*f.inner).clone().into();
        f.into()
    }
}

// impl<F: Kind, T: Kind + Into<F>> From<Expr<F>> for Expr<T> {
//     fn from(f: Expr<F>) -> Self {
//         let f: T = (*f.inner).clone().into();
//         f.into()
//     }
// }
