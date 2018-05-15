extern crate internment;
extern crate tinyset;

use internment::Intern;
use tinyset::{Map64, Fits64};
use tinyset::u64set::Map64Iter;
use std::fmt::Debug;
use std::hash::Hash;
use std::iter::FromIterator;
use std::ops::Deref;

pub trait ExprType: 'static + Send + Clone + Eq + Debug + Hash {
    fn cpp(&self) -> String;
}

#[derive(Clone, Debug, Eq, Hash)]
pub struct Expr<T: ExprType> {
    inner: Intern<T>,
}

impl<T: ExprType> Expr<T> {
    fn new(inner: &T) -> Expr<T> {
        Expr { inner: Intern::new(inner.clone()), }
    }

    fn cpp(&self) -> String {
        self.inner.cpp()
    }

    fn cast<U: ExprType + From<T>>(&self) -> Expr<U> {
        Expr::new(&self.inner.deref().clone().into())
    }
}

impl<T: ExprType, U: Clone + Into<Expr<T>>> PartialEq<U> for Expr<T> {
    fn eq(&self, other: &U) -> bool {
        self.inner == other.clone().into().inner
    }
}

impl<T: ExprType> Fits64 for Expr<T> {
    fn to_u64(self) -> u64 {
        self.inner.to_u64()
    }

    unsafe fn from_u64(x: u64) -> Self {
        Expr { inner: Intern::from_u64(x), }
    }
}

impl<T: ExprType> Copy for Expr<T> {}

impl<T: ExprAdd + ExprMul, N: Into<f64>> From<N> for Expr<T> {
    fn from(n: N) -> Self {
        let n = n.into();
        let inner = if n == 0.0 {
            T::zero()
        } else if n == 1.0 {
            T::one()
        } else {
            let mut map = AbelianMap::new();
            map.insert(Expr::new(&T::one()), n);
            T::sum_from_map(map)
        };
        Expr::new(&inner)
    }
}

pub trait ExprAdd: ExprType {
    fn sum_from_map(x: AbelianMap<Self>) -> Self;
    fn map_from_sum(&self) -> Option<&AbelianMap<Self>>;

    fn add(&self, other: &Self) -> Self {
        let mut sum: AbelianMap<Self>;
        match (self.map_from_sum(), other.map_from_sum()) {
            (Some(l), Some(r)) => {
                sum = l.clone();
                sum.union(r);
            },
            (Some(s), _) => {
                sum = s.clone();
                sum.insert(Expr::new(other), 1.0);
            },
            (_, Some(s)) => {
                sum = s.clone();
                sum.insert(Expr::new(self), 1.0);
            },
            (_, _) => {
                sum = AbelianMap::new();
                sum.insert(Expr::new(self), 1.0);
                sum.insert(Expr::new(other), 1.0);
            },
        }
        Self::sum_from_map(sum)
    }

    fn neg(&self) -> Self {
        match self.map_from_sum() {
            Some(ref map) =>
                Self::sum_from_map(map.iter().map(|(k, &v)| (k, -v)).collect()),
            _ => if *self == Self::zero() {
                    Self::zero()
                } else {
                    let mut map = AbelianMap::new();
                    map.insert(Expr::new(self), -1.0);
                    Self::sum_from_map(map)
                },
        }
    }

    fn zero() -> Self {
        Self::sum_from_map(AbelianMap::new())
    }
}

impl<T: ExprAdd, U: Into<Expr<T>>> std::ops::Add<U> for Expr<T> {
    type Output = Self;

    fn add(self, other: U) -> Self {
        let other = other.into();
        Expr::new(&self.inner.deref().add(other.inner.deref()))
    }
}

impl<T: ExprAdd, U: Into<Expr<T>>> std::ops::Sub<U> for Expr<T> {
    type Output = Self;

    fn sub(self, other: U) -> Self {
        let other = other.into();
        Expr::new(&self.inner.deref().add(&other.inner.deref().neg()))
    }
}

impl<T: ExprAdd> std::ops::Neg for Expr<T> {
    type Output = Self;

    fn neg(self) -> Self {
        Expr::new(&self.inner.deref().neg())
    }
}

pub trait ExprMul: ExprType {
    fn mul_from_map(x: AbelianMap<Self>) -> Self;
    fn map_from_mul(&self) -> Option<&AbelianMap<Self>>;

    fn mul(&self, other: &Self) -> Self {
        let mut prod: AbelianMap<Self>;
        match (self.map_from_mul(), other.map_from_mul()) {
            (Some(l), Some(r)) => {
                prod = l.clone();
                prod.union(r);
            },
            (Some(s), _) => {
                prod = s.clone();
                prod.insert(Expr::new(other), 1.0);
            },
            (_, Some(s)) => {
                prod = s.clone();
                prod.insert(Expr::new(self), 1.0);
            },
            (_, _) => {
                prod = AbelianMap::new();
                prod.insert(Expr::new(self), 1.0);
                prod.insert(Expr::new(other), 1.0);
            },
        }
        Self::mul_from_map(prod)
    }

    fn recip(&self) -> Self {
        match self.map_from_mul() {
            Some(ref map) =>
                Self::mul_from_map(map.iter().map(|(k, &v)| (k, -v)).collect()),
            _ => if *self == Self::one() {
                    Self::one()
                } else {
                    let mut map = AbelianMap::new();
                    map.insert(Expr::new(self), -1.0);
                    Self::mul_from_map(map)
                },
        }
    }

    fn one() -> Self {
        Self::mul_from_map(AbelianMap::new())
    }
}

impl<T: ExprMul, U: Into<Expr<T>>> std::ops::Mul<U> for Expr<T> {
    type Output = Self;

    fn mul(self, other: U) -> Self {
        let other = other.into();
        Expr::new(&self.inner.deref().mul(other.inner.deref()))
    }
}

impl<T: ExprMul, U: Into<Expr<T>>> std::ops::Div<U> for Expr<T> {
    type Output = Self;

    fn div(self, other: U) -> Self {
        let other = other.into();
        Expr::new(&self.inner.deref().mul(&other.inner.deref().recip()))
    }
}

macro_rules! impl_expr_add {
    ( $name:ident ) => {
        impl ExprAdd for $name {
            fn sum_from_map(map: AbelianMap<Self>) -> Self {
                if map.len() == 1 {
                    let (k, &v) = map.iter().next().unwrap();
                    if v == 1.0 {
                        return k.inner.deref().clone();
                    }
                }
                $name::Add(map)
            }

            fn map_from_sum(&self) -> Option<&AbelianMap<Self>> {
                if let &$name::Add(ref map) = self {
                    Some(&map)
                } else {
                    None
                }
            }
        }
    }
}

macro_rules! impl_expr_mul {
    ( $name:ident ) => {
        impl ExprMul for $name {
            fn mul_from_map(map: AbelianMap<Self>) -> Self {
                if map.len() == 1 {
                    let (k, &v) = map.iter().next().unwrap();
                    if v == 1.0 {
                        return k.inner.deref().clone();
                    }
                }
                $name::Mul(map)
            }

            fn map_from_mul(&self) -> Option<&AbelianMap<Self>> {
                if let &$name::Mul(ref map) = self {
                    Some(&map)
                } else {
                    None
                }
            }
        }
    }
}

#[derive(Clone, PartialEq, Eq, Debug, Hash)]
pub enum Scalar {
    Var(Intern<String>),
    Exp(Expr<Scalar>),
    Log(Expr<Scalar>),
    Add(AbelianMap<Scalar>),
    Mul(AbelianMap<Scalar>),
}

impl_expr_add!(Scalar);
impl_expr_mul!(Scalar);
impl_expr_add!(RealSpaceScalar);
impl_expr_mul!(RealSpaceScalar);
impl_expr_add!(KSpaceScalar);
impl_expr_mul!(KSpaceScalar);

impl Scalar {
    fn exp(&self) -> Self { Scalar::Exp(Expr::new(self)) }
    fn log(&self) -> Self { Scalar::Exp(Expr::new(self)) }
    fn var(name: &str) -> Self { Scalar::Var(Intern::new(String::from(name))) }
}

impl ExprType for Scalar {
    fn cpp(&self) -> String {
        match self {
            &Scalar::Var(s) =>
                (*s).clone(),
            &Scalar::Exp(a) =>
                String::from("exp(") + &a.cpp() + &")",
            &Scalar::Log(a) =>
                String::from("log(") + &a.cpp() + &")",
            &Scalar::Add(ref m) => {
                let pcoeff = |&(ref x, ref c): &(String, f64)| -> String {
                    if x == "1" {
                        c.to_string()
                    } else if *c == 1.0 {
                        x.clone()
                    } else {
                        c.to_string() + &" * " + &x
                    }
                };
                let ncoeff = |&(ref x, ref c): &(String, f64)| -> String {
                    if x == "1" {
                        c.abs().to_string()
                    } else if *c == -1.0 {
                        x.clone()
                    } else {
                        c.abs().to_string() + &" * " + &x
                    }
                };
                let (p, n) = m.split_cpp_sort();
                match (p.len(), n.len()) {
                    (0, 0) =>
                        String::from("0"),
                    (_, 0) =>
                        p.iter().map(pcoeff).collect::<Vec<String>>().join(" + "),
                    (0, _) =>
                        String::from("-")
                            + &n.iter().map(ncoeff).collect::<Vec<String>>().join(" + "),
                    (_, _) =>
                        p.iter().map(pcoeff).collect::<Vec<String>>().join(" + ")
                            + &" - "
                            + &n.iter().map(ncoeff).collect::<Vec<String>>().join(" + "),
                }
            },
            &Scalar::Mul(ref m) => {
                let ref power = |&(ref x, ref p): &(String, f64)| -> String {
                    if x == "1" || p.abs() == 1.0 {
                        x.clone()
                    } else if p.abs() == 2.0 {
                        x.clone() + &" * " + &x
                    } else {
                        String::from("pow(") + &x + &", " + &p.abs().to_string() + &")"
                    }
                };
                let (n, d) = m.split_cpp_sort();
                match (n.len(), d.len()) {
                    (0, 0) =>
                        String::from("1"),
                    (_, 0) =>
                        n.iter().map(power).collect::<Vec<String>>().join(" * "),
                    (0, _) =>
                        String::from("1 / (")
                            + &d.iter().map(power).collect::<Vec<String>>().join(" * ")
                            + &")",
                    (_, _) =>
                        n.iter().map(power).collect::<Vec<String>>().join(" * ")
                            + &" / ("
                            + &d.iter().map(power).collect::<Vec<String>>().join(" * ")
                            + &")",
                }
            },
        }
    }
}

#[derive(Clone, PartialEq, Eq, Debug, Hash)]
pub enum RealSpaceScalar {
    Var(Intern<String>),
    ScalarVar(Intern<String>),
    Exp(Expr<RealSpaceScalar>),
    Log(Expr<RealSpaceScalar>),
    Add(AbelianMap<RealSpaceScalar>),
    Mul(AbelianMap<RealSpaceScalar>),
    IFFT(Expr<KSpaceScalar>),
}

impl RealSpaceScalar {
    fn exp(&self) -> Self { RealSpaceScalar::Exp(Expr::new(self)) }
    fn log(&self) -> Self { RealSpaceScalar::Exp(Expr::new(self)) }
    fn var(name: &str) -> Self { RealSpaceScalar::Var(Intern::new(String::from(name))) }
    fn scalar_var(name: &str) -> Self { RealSpaceScalar::ScalarVar(Intern::new(String::from(name))) }
}

impl ExprType for RealSpaceScalar { fn cpp(&self) -> String { unimplemented!() } }

#[derive(Clone, PartialEq, Eq, Debug, Hash)]
pub enum KSpaceScalar {
    Var(Intern<String>),
    ScalarVar(Intern<String>),
    Exp(Expr<KSpaceScalar>),
    Log(Expr<KSpaceScalar>),
    Add(AbelianMap<KSpaceScalar>),
    Mul(AbelianMap<KSpaceScalar>),
    FFT(Expr<RealSpaceScalar>),
}

impl KSpaceScalar {
    fn exp(&self) -> Self { KSpaceScalar::Exp(Expr::new(self)) }
    fn log(&self) -> Self { KSpaceScalar::Exp(Expr::new(self)) }
    fn var(name: &str) -> Self { KSpaceScalar::Var(Intern::new(String::from(name))) }
    fn scalar_var(name: &str) -> Self { KSpaceScalar::ScalarVar(Intern::new(String::from(name))) }
}

impl ExprType for KSpaceScalar { fn cpp(&self) -> String { unimplemented!() } }

impl From<Scalar> for RealSpaceScalar {
    fn from(s: Scalar) -> Self {
        match s {
            Scalar::Var(name) =>
                RealSpaceScalar::ScalarVar(name),
            Scalar::Exp(a) =>
                RealSpaceScalar::Exp(Expr::new(&a.inner.deref().clone().into())),
            Scalar::Log(a) =>
                RealSpaceScalar::Log(Expr::new(&a.inner.deref().clone().into())),
            Scalar::Add(m) =>
                RealSpaceScalar::Add(
                    m.iter()
                     .map(|(k, &v)| (Expr::new(&k.inner.deref().clone().into()), v))
                     .collect()),
            Scalar::Mul(m) =>
                RealSpaceScalar::Mul(
                    m.iter()
                     .map(|(k, &v)| (Expr::new(&k.inner.deref().clone().into()), v))
                     .collect()),
        }
    }
}

// impl RealSpaceScalar {
//     fn fake_scalar(&self) -> Scalar {
//         match s {
//             RealSpaceScalar::Var(name) =>
//                 Scalar::ScalarVar(
//                     // Leak `s` because `Var`s are static because
//                     // they are in `Intern`ed `Expr::inner`s.
//                     // See <https://stackoverflow.com/a/30527289>.
//                     unsafe {
//                         let s = String::from(s) + &"[i]";
//                         let ss = std::mem::transmute(&s as &str);
//                         std::mem::forget(s);
//                         ss
//                     }),
//             RealSpaceScalar::Exp(a) =>
//                 Scalar::Exp(Expr::new(&a.inner.deref().clone().into())),
//             RealSpaceScalar::Log(a) =>
//                 Scalar::Log(Expr::new(&a.inner.deref().clone().into())),
//             RealSpaceScalar::Add(m) =>
//                 Scalar::Add(
//                     m.iter()
//                      .map(|(k, &v)| (Expr::new(&k.inner.deref().clone().into()), v))
//                      .collect()),
//             RealSpaceScalar::Mul(m) =>
//                 Scalar::Mul(
//                     m.iter()
//                      .map(|(k, &v)| (Expr::new(&k.inner.deref().clone().into()), v))
//                      .collect()),
//         }
//     }
// }

#[derive(Clone)]
pub struct AbelianMap<T: ExprType> {
    inner: Map64<Expr<T>, f64>,
}

impl<T: ExprType> AbelianMap<T> {
    fn new() -> Self {
        Self { inner: Map64::new(), }
    }

    fn iter(&self) -> <&Self as IntoIterator>::IntoIter {
        self.into_iter()
    }

    fn len(&self) -> usize {
        self.inner.len()
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
        for (k, &v) in other {
            self.insert(k, v);
        }
    }

    fn split_cpp_sort(&self) -> (Vec<(String, f64)>, Vec<(String, f64)>) {
        let (mut p, mut n): (Vec<(String, f64)>, Vec<(String, f64)>)
            = self.iter()
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

impl<T: ExprType> Debug for AbelianMap<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut res = write!(f, "AbelianMap [ ");
        for term in self {
            res = res.and(write!(f, "{:?}, ", term));
        }
        res.and(write!(f, "] "))
    }
}

impl<'a, T: ExprType> IntoIterator for &'a AbelianMap<T> {
    type Item = (Expr<T>, &'a f64);
    type IntoIter = Map64Iter<'a, Expr<T>, f64>;

    fn into_iter(self) -> Self::IntoIter {
        self.inner.iter()
    }
}

// May not be needed.
//
// impl<'a, T: ExprType> FromIterator<(Expr<T>, &'a f64)> for AbelianMap<T> {
//     fn from_iter<I: IntoIterator<Item = (Expr<T>, &'a f64)>>(iter: I) -> Self {
//         let mut map = AbelianMap::new();
//         for (k, &v) in iter {
//             map.insert(k, v);
//         }
//         map
//     }
// }

impl<'a, T: ExprType> FromIterator<(Expr<T>, f64)> for AbelianMap<T> {
    fn from_iter<I: IntoIterator<Item = (Expr<T>, f64)>>(iter: I) -> Self {
        let mut map = AbelianMap::new();
        for (k, v) in iter {
            map.insert(k, v);
        }
        map
    }
}

// May not be needed.
//
// impl<'a, T: ExprType> FromIterator<Expr<T>> for AbelianMap<T> {
//     fn from_iter<I: IntoIterator<Item = Expr<T>>>(iter: I) -> Self {
//         let mut map = AbelianMap::new();
//         for k in iter {
//             map.insert(k, 1.0);
//         }
//         map
//     }
// }

impl<'a, T: ExprType> FromIterator<T> for AbelianMap<T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let mut map = AbelianMap::new();
        for k in iter {
            map.insert(Expr::new(&k), 1.0);
        }
        map
    }
}

impl<T: ExprType> Hash for AbelianMap<T> {
    fn hash<H: std::hash::Hasher>(&self, hasher: &mut H) {
        let mut terms: Vec<_>
            = self.iter()
                  .map(|(k, &v)| (k.inner.to_u64(), v as u64))
                  .collect();
        terms.sort();
        for (k, v) in terms {
            k.hash(hasher);
            v.hash(hasher);
        }
    }
}

impl<T: ExprType> PartialEq for AbelianMap<T> {
    fn eq(&self, other: &Self) -> bool {
        let mut lhs: Vec<_>
            = self.iter()
                  .map(|(k, &v)| (k.inner.to_u64(), v as u64))
                  .collect();
        let mut rhs: Vec<_>
            = other.iter()
                   .map(|(k, &v)| (k.inner.to_u64(), v as u64))
                   .collect();
        lhs.sort();
        rhs.sort();
        lhs == rhs
    }
}

impl<T: ExprType> Eq for AbelianMap<T> {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn abelian() {
        let t = vec![Scalar::var("x"), Scalar::var("x"), Scalar::var("y"), Scalar::var("z")];
        let s: AbelianMap<Scalar> = t.into_iter().collect();
        assert_eq!(s, s);
        assert_eq!(s, s.iter().map(|(k, &v)| (k, v)).collect());
        assert!(s != vec![Scalar::var("x"), Scalar::var("y"), Scalar::var("z")].into_iter().collect());
    }

    #[test]
    fn arith() {
        let a = Expr::new(&Scalar::var("a"));
        let b = Expr::new(&Scalar::var("b"));

        assert_eq!(a, a);
        assert!(a != b);

        assert_eq!(a + 0, a);
        assert_eq!(a + a, a + a);
        assert_eq!(a + b, a + b);
        assert_eq!(a + b, b + a);
        assert_eq!(a + a - a, a);
        assert_eq!(a - a, 0);
        assert_eq!(a + a - a - a, 0);

        assert_eq!(a * 1, a);
        assert_eq!(a * a, a * a);
        assert_eq!(a * b, a * b);
        assert_eq!(a * b, b * a);
        assert_eq!(a * a / a, a);
        assert_eq!(a / a, 1);
        assert_eq!(a * a / a / a, 1);
    }

    #[test]
    fn realspace() {
        let a: Expr<RealSpaceScalar> = <Expr<RealSpaceScalar>>::new(&RealSpaceScalar::var("a"));
        let s: Expr<Scalar> = <Expr<Scalar>>::new(&Scalar::var("s"));
        let rs_s: Expr<RealSpaceScalar> =
            <Expr<RealSpaceScalar>>::new(&RealSpaceScalar::scalar_var("s"));

        assert_eq!(a * rs_s, a * rs_s);
        // FIXME the following is broken:
        // assert_eq!(a * s, a * rs_s);
    }

    #[test]
    fn cpp() {
        let a = Expr::new(&Scalar::var("a"));
        let b = Expr::new(&Scalar::var("b"));
        let c = Expr::new(&Scalar::var("c"));


        assert_eq!((c + b + a).cpp(), "a + b + c");
        assert_eq!((c * b * a).cpp(), "a * b * c");
        assert_eq!((a / c / b).cpp(), "a / (b * c)");
    }
}
