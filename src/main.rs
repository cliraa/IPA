use lambdaworks_math::unsigned_integer::element::U384;
use lambdaworks_math::field::{
    element::FieldElement,
    fields::montgomery_backed_prime_fields::{IsModulus, U384PrimeField}};

#[derive(Clone, Debug)]
pub struct U384ModulusP;

impl IsModulus<U384> for U384ModulusP {
    const MODULUS: U384 = U384::from_hex_unchecked("47");
}

type U384FP = U384PrimeField<U384ModulusP>;
type U384FPElement = FieldElement<U384FP>;

fn calculate_p(
    a: [U384FPElement; 2], 
    b: [U384FPElement; 2], 
    g: [U384FPElement; 2], 
    h: [U384FPElement; 2]) -> U384FPElement {
    let p: U384FPElement = (a[0].clone() * g[0].clone()) + (b[0].clone() * h[0].clone()) + (a[1].clone() * g[1].clone()) + (b[1].clone() * h[1].clone());
    p
}

fn calculate_l(
    a: [U384FPElement; 2], 
    b: [U384FPElement; 2], 
    g: [U384FPElement; 2], 
    h: [U384FPElement; 2], u: U384FPElement) -> U384FPElement {
    let l: U384FPElement = (a[0].clone() * g[1].clone()) + (b[1].clone() * h[0].clone()) + (a[0].clone() * b[1].clone() * u.clone());
    l
}

fn calculate_r(
    a: [U384FPElement; 2],
    b: [U384FPElement; 2],
    g: [U384FPElement; 2],
    h: [U384FPElement; 2],
    u: U384FPElement) -> U384FPElement {
        let r: U384FPElement = (a[1].clone() * g[0].clone()) + (b[0].clone() * h[1].clone()) + (a[1].clone() * b[0].clone() * u.clone());
        r
}

fn calculate_a_prime(
    a: [U384FPElement; 2],
    x: U384FPElement) -> U384FPElement {
        let a_prime: U384FPElement = (a[0].clone() * x.clone()) + (a[1].clone() * mod_inverse(x.clone()));
        a_prime
}

fn calculate_b_prime(
    b: [U384FPElement; 2],
    x: U384FPElement) -> U384FPElement {
        let b_prime: U384FPElement = (b[0].clone() * mod_inverse(x.clone())) + (b[1].clone() * x.clone());
        b_prime
}

fn calculate_c(
    a: [U384FPElement; 2], 
    b: [U384FPElement; 2]) -> U384FPElement {
        let c: U384FPElement = (a[0].clone() * b[0].clone()) + (a[1].clone() * b[1].clone());
        c
}

fn verify_n_2(
    a: [U384FPElement; 2],
    b: [U384FPElement; 2],
    g: [U384FPElement; 2],
    h: [U384FPElement; 2],
    u: U384FPElement, x: U384FPElement) -> bool {

        let p: U384FPElement = calculate_p(a.clone(), b.clone(), g.clone(), h.clone());
        let l: U384FPElement = calculate_l(a.clone(), b.clone(), g.clone(), h.clone(), u.clone());
        let r: U384FPElement = calculate_r(a.clone(), b.clone(), g.clone(), h.clone(), u.clone());
        let a_prime: U384FPElement = calculate_a_prime(a.clone(), x.clone());
        let b_prime: U384FPElement = calculate_b_prime(b.clone(), x.clone());
        
        // Inner Product Between A and B
        let c: U384FPElement = calculate_c(a, b);

        // Verification
        let left: U384FPElement = (x.clone().pow(2_u64)) * l.clone() + p.clone() + (c.clone() * u.clone()) + (mod_inverse(x.clone()).pow(2_u64)) * r.clone();
        let right: U384FPElement = (mod_inverse(x.clone()) * a_prime.clone() * g[0].clone()) + (x.clone() * a_prime.clone() * g[1].clone()) + (x.clone() * b_prime.clone() * h[0].clone()) + (mod_inverse(x) * b_prime.clone() * h[1].clone()) + (a_prime.clone() * b_prime.clone() * u.clone());

        left == right
}

fn mod_inverse(a: U384FPElement) -> U384FPElement {
    let mod_minustwo: U384FPElement = U384FPElement::from_hex_unchecked("45");
    let result: U384FPElement = a.pow(mod_minustwo.representative());
    result
}

#[cfg(test)]
mod tests {

    use crate::verify_n_2;
    use crate::U384FPElement;

    #[test]
    fn verify_test_1() {
        let a: [U384FPElement; 2] = [U384FPElement::from_hex_unchecked("2"), U384FPElement::from_hex_unchecked("3")];
        let b: [U384FPElement; 2] = [U384FPElement::from_hex_unchecked("4"), U384FPElement::from_hex_unchecked("6")];
        let g: [U384FPElement; 2] = [U384FPElement::from_hex_unchecked("5"), U384FPElement::from_hex_unchecked("7")];
        let h: [U384FPElement; 2] = [U384FPElement::from_hex_unchecked("b"), U384FPElement::from_hex_unchecked("a")];
        let u: U384FPElement = U384FPElement::from_hex_unchecked("5");
        let x: U384FPElement = U384FPElement::from_hex_unchecked("3");
    
        let result: bool = verify_n_2(a, b, g, h, u, x);
        assert_eq!(result, true);
    }

    #[test]
    fn verify_test_2() {
        let a: [U384FPElement; 2] = [U384FPElement::from_hex_unchecked("21"), U384FPElement::from_hex_unchecked("13")];
        let b: [U384FPElement; 2] = [U384FPElement::from_hex_unchecked("3"), U384FPElement::from_hex_unchecked("24")];
        let g: [U384FPElement; 2] = [U384FPElement::from_hex_unchecked("10"), U384FPElement::from_hex_unchecked("17")];
        let h: [U384FPElement; 2] = [U384FPElement::from_hex_unchecked("5"), U384FPElement::from_hex_unchecked("3")];
        let u: U384FPElement = U384FPElement::from_hex_unchecked("f");
        let x: U384FPElement = U384FPElement::from_hex_unchecked("d");
    
        let result: bool = verify_n_2(a, b, g, h, u, x);
        assert_eq!(result, true);
    }
}
