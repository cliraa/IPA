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
    a: Vec<U384FPElement>,
    b: Vec<U384FPElement>,
    g: Vec<U384FPElement>,
    h: Vec<U384FPElement>) -> U384FPElement {
    let p: U384FPElement = (a[0].clone() * g[0].clone()) + (b[0].clone() * h[0].clone()) + (a[1].clone() * g[1].clone()) + (b[1].clone() * h[1].clone());
    p
}

fn calculate_l(
    a: Vec<U384FPElement>,
    b: Vec<U384FPElement>,
    g: Vec<U384FPElement>,
    h: Vec<U384FPElement>, u: U384FPElement) -> U384FPElement {
    let l: U384FPElement = (a[0].clone() * g[1].clone()) + (b[1].clone() * h[0].clone()) + (a[0].clone() * b[1].clone() * u.clone());
    l
}

fn calculate_r(
    a: Vec<U384FPElement>,
    b: Vec<U384FPElement>,
    g: Vec<U384FPElement>,
    h: Vec<U384FPElement>,
    u: U384FPElement) -> U384FPElement {
        let r: U384FPElement = (a[1].clone() * g[0].clone()) + (b[0].clone() * h[1].clone()) + (a[1].clone() * b[0].clone() * u.clone());
        r
}

fn calculate_a_prime(
    a: Vec<U384FPElement>,
    x: U384FPElement) -> U384FPElement {
        let a_prime: U384FPElement = (a[0].clone() * x.clone()) + (a[1].clone() * mod_inverse(x.clone()));
        a_prime
}

fn calculate_b_prime(
    b: Vec<U384FPElement>,
    x: U384FPElement) -> U384FPElement {
        let b_prime: U384FPElement = (b[0].clone() * mod_inverse(x.clone())) + (b[1].clone() * x.clone());
        b_prime
}

fn calculate_c(
    a: &mut Vec<U384FPElement>,
    b: &mut Vec<U384FPElement>) -> U384FPElement {
        let c: U384FPElement = (a[0].clone() * b[0].clone()) + (a[1].clone() * b[1].clone());
        c
}

fn verify_n_2(
    a: &mut Vec<U384FPElement>,
    b: &mut Vec<U384FPElement>,
    g: &mut Vec<U384FPElement>,
    h: &mut Vec<U384FPElement>,
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

fn inner_product(v1: &[U384FPElement], v2: &[U384FPElement]) -> U384FPElement {
    if v1.len() != v2.len() {
        panic!("Vectors must have the same size.");
    }
    let result = v1.iter().zip(v2.iter()).map(|(x, y)| x * y).sum::<U384FPElement>();
    result
}

fn vec_low(vector: &[U384FPElement]) -> Vec<U384FPElement> {
    let mid_idx = vector.len() / 2;
    vector[..mid_idx].to_vec()
}

fn vec_hi(vector: &[U384FPElement]) -> Vec<U384FPElement> {
    let mid_idx = vector.len() / 2;
    vector[mid_idx..].to_vec()
}

fn computek(vector: &[U384FPElement]) -> u32 {
    let mut k = 0;
    let mut length = vector.len();
    while length != 1 {
        length /= 2;
        k += 1;
    }
    k
}

fn multiply_vector(vector: &[U384FPElement], s: U384FPElement) -> Vec<U384FPElement> {
    vector.iter().map(|x| x * s.clone()).collect()
}

fn vector_addition(a: &[U384FPElement], b: &[U384FPElement]) -> Vec<U384FPElement> {
    if a.len() != b.len() {
        panic!("Error: Vectors must have the same length");
    }
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}

fn verify(a: &mut Vec<U384FPElement>, b: &mut Vec<U384FPElement>, g: &mut Vec<U384FPElement>, h: &mut Vec<U384FPElement>, u: U384FPElement, s: U384FPElement, s_prime: U384FPElement, x: U384FPElement) -> bool {
    if a.len() == 2 {
        verify_n_2(a, b, g, h, u, x)
    } else {
    let mut k = computek(a);
    while k > 1 {
        let a_lo = vec_low(a);
        let a_hi = vec_hi(a);
        let b_lo = vec_low(b);
        let b_hi = vec_hi(b);
        let g_lo = vec_low(g);
        let g_hi = vec_hi(g);
        let h_lo = vec_low(h);
        let h_hi = vec_hi(h);
        
        let L = inner_product(&a_lo, &g_hi) + inner_product(&b_hi, &h_lo) + inner_product(&a_lo, &b_hi) * u.clone();
        let R = inner_product(&a_hi, &g_lo) + inner_product(&b_lo, &h_hi) + inner_product(&a_hi, &b_lo) * u.clone();

        let a_prime = vector_addition(&multiply_vector(&a_lo, x.clone()), &multiply_vector(&a_hi, mod_inverse(x.clone())));
        let b_prime = vector_addition(&multiply_vector(&b_lo, mod_inverse(x.clone())),&multiply_vector(&b_hi, x.clone()));
        let g_prime = vector_addition(&multiply_vector(&g_lo, mod_inverse(x.clone())),&multiply_vector(&g_hi, x.clone()));
        let h_prime = vector_addition(&multiply_vector(&h_lo, x.clone()),&multiply_vector(&h_hi, mod_inverse(x.clone())));
        
        *a = a_prime;
        *b = b_prime;
        *g = g_prime;
        *h = h_prime;
        k -= 1;
    }
    verify_n_2(a, b, g, h, u, x)
    }
}

#[cfg(test)]
mod tests {

    use crate::verify_n_2;
    use crate::verify;
    use crate::U384FPElement;

    #[test]
    fn verify_test_1() {
        let mut a = vec![U384FPElement::from_hex_unchecked("2"), U384FPElement::from_hex_unchecked("3")];
        let mut b = vec![U384FPElement::from_hex_unchecked("4"), U384FPElement::from_hex_unchecked("6")];
        let mut g = vec![U384FPElement::from_hex_unchecked("5"), U384FPElement::from_hex_unchecked("7")];
        let mut h = vec![U384FPElement::from_hex_unchecked("b"), U384FPElement::from_hex_unchecked("a")];
        let u: U384FPElement = U384FPElement::from_hex_unchecked("5");
        let x: U384FPElement = U384FPElement::from_hex_unchecked("3");
    
        let result: bool = verify_n_2(&mut a, &mut b, &mut g, &mut h, u, x);
        assert_eq!(result, true);
    }    

    #[test]
    fn verify_test_2() {
        let mut a = vec![U384FPElement::from_hex_unchecked("21"), U384FPElement::from_hex_unchecked("13")];
        let mut b = vec![U384FPElement::from_hex_unchecked("3"), U384FPElement::from_hex_unchecked("24")];
        let mut g = vec![U384FPElement::from_hex_unchecked("10"), U384FPElement::from_hex_unchecked("17")];
        let mut h = vec![U384FPElement::from_hex_unchecked("5"), U384FPElement::from_hex_unchecked("3")];
        let u: U384FPElement = U384FPElement::from_hex_unchecked("f");
        let x: U384FPElement = U384FPElement::from_hex_unchecked("d");
    
        let result: bool = verify_n_2(&mut a, &mut b, &mut g, &mut h, u, x);
        assert_eq!(result, true);
    }

    #[test]
    fn verify_test_3() {
        let mut a = 
            vec![U384FPElement::from_hex_unchecked("21"), U384FPElement::from_hex_unchecked("13"), U384FPElement::from_hex_unchecked("2"),U384FPElement::from_hex_unchecked("3")];
        let mut b = 
            vec![U384FPElement::from_hex_unchecked("10"), U384FPElement::from_hex_unchecked("3"), U384FPElement::from_hex_unchecked("a"),U384FPElement::from_hex_unchecked("b")];
        let mut g = 
            vec![U384FPElement::from_hex_unchecked("21"), U384FPElement::from_hex_unchecked("12"), U384FPElement::from_hex_unchecked("3"),U384FPElement::from_hex_unchecked("6")];
        let mut h = 
            vec![U384FPElement::from_hex_unchecked("20"), U384FPElement::from_hex_unchecked("11"), U384FPElement::from_hex_unchecked("4"),U384FPElement::from_hex_unchecked("5")];
        let u: U384FPElement = U384FPElement::from_hex_unchecked("f");
        let s: U384FPElement = U384FPElement::from_hex_unchecked("e");
        let s_prime: U384FPElement = U384FPElement::from_hex_unchecked("c");
        let x: U384FPElement = U384FPElement::from_hex_unchecked("d");
    
        let result = verify(&mut a, &mut b, &mut g, &mut h, u, s, s_prime, x);
        assert_eq!(result, true);
    }
    
}
