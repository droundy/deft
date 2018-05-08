extern crate proc_macro;
extern crate syn;
#[macro_use]
extern crate quote;

use proc_macro::TokenStream;

#[proc_macro_derive(ExprAdd)]
pub fn derive_expr_add(input: TokenStream) -> TokenStream {
    let s = input.to_string();
    let ast = syn::parse_derive_input(&s).unwrap();
    let gen = impl_expr_add(&ast);
    gen.parse().unwrap()
}

fn impl_expr_add(ast: &syn::DeriveInput) -> quote::Tokens {
    let name = &ast.ident;
    quote! {
        impl ExprAdd for #name {
            fn sum_from_map(map: AbelianMap<Self>) -> Self {
                if map.len() == 1 {
                    let (k, &v) = map.iter().next().unwrap();
                    if v == 1.0 {
                        return k.inner.deref().clone();
                    }
                }
                #name::Add(map)
            }

            fn map_from_sum(&self) -> Option<&AbelianMap<Self>> {
                if let &#name::Add(ref map) = self {
                    Some(&map)
                } else {
                    None
                }
            }
        }
    }
}

#[proc_macro_derive(ExprMul)]
pub fn derive_expr_mul(input: TokenStream) -> TokenStream {
    let s = input.to_string();
    let ast = syn::parse_derive_input(&s).unwrap();
    let gen = impl_expr_mul(&ast);
    gen.parse().unwrap()
}

fn impl_expr_mul(ast: &syn::DeriveInput) -> quote::Tokens {
    let name = &ast.ident;
    quote! {
        impl ExprMul for #name {
            fn mul_from_map(map: AbelianMap<Self>) -> Self {
                if map.len() == 1 {
                    let (k, &v) = map.iter().next().unwrap();
                    if v == 1.0 {
                        return k.inner.deref().clone();
                    }
                }
                #name::Mul(map)
            }

            fn map_from_mul(&self) -> Option<&AbelianMap<Self>> {
                if let &#name::Mul(ref map) = self {
                    Some(&map)
                } else {
                    None
                }
            }
        }
    }
}
