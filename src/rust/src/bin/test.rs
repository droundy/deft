//! This is the deft crate!

extern crate deft;

use deft::Expr;

fn main() {
    let a = Expr::var("a");
    let b = Expr::var("b");
    let z = Expr::var("z");
    println!("a: {:?}", a);
    println!("b: {:?}", b);
    let sum = a+b;
    println!("a+b: {:?}", sum);
    println!("a+b: {}", sum.cpp());
    assert_eq!((a+b).cpp(), "a + b");
    assert_eq!((a+z+b).cpp(), "a + b + z");
    assert_eq!((a+a).cpp(), "2 * a");
}
