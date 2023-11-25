use lambdaworks_math::field::element::FieldElement;

fn calculate_p<F: lambdaworks_math::field::traits::IsField>(
    a: Vec<FieldElement<F>>,
    b: Vec<FieldElement<F>>,
    g: Vec<FieldElement<F>>,
    h: Vec<FieldElement<F>>) -> FieldElement<F> {
    let p: FieldElement<F> = (a[0].clone() * g[0].clone()) + (b[0].clone() * h[0].clone()) + (a[1].clone() * g[1].clone()) + (b[1].clone() * h[1].clone());
    p
}

fn calculate_l<F: lambdaworks_math::field::traits::IsField>(
    a: Vec<FieldElement<F>>,
    b: Vec<FieldElement<F>>,
    g: Vec<FieldElement<F>>,
    h: Vec<FieldElement<F>>,
    u: FieldElement<F>) -> FieldElement<F> {
    let l: FieldElement<F> = (a[0].clone() * g[1].clone()) + (b[1].clone() * h[0].clone()) + (a[0].clone() * b[1].clone() * u.clone());
    l
}

fn calculate_r<F: lambdaworks_math::field::traits::IsField>(
    a: Vec<FieldElement<F>>,
    b: Vec<FieldElement<F>>,
    g: Vec<FieldElement<F>>,
    h: Vec<FieldElement<F>>,
    u: FieldElement<F>) -> FieldElement<F> {
        let r: FieldElement<F> = (a[1].clone() * g[0].clone()) + (b[0].clone() * h[1].clone()) + (a[1].clone() * b[0].clone() * u.clone());
        r
}

fn calculate_a_prime<F: lambdaworks_math::field::traits::IsField>(
    a: Vec<FieldElement<F>>,
    x: FieldElement<F>) -> FieldElement<F> {
        let a_prime: FieldElement<F> = (a[0].clone() * x.clone()) + (a[1].clone() * mod_inverse(x.clone()));
        a_prime
}

fn calculate_b_prime<F: lambdaworks_math::field::traits::IsField>(
    b: Vec<FieldElement<F>>,
    x: FieldElement<F>) -> FieldElement<F> {
        let b_prime: FieldElement<F> = (b[0].clone() * mod_inverse(x.clone())) + (b[1].clone() * x.clone());
        b_prime
}

fn calculate_c<F: lambdaworks_math::field::traits::IsField>(
    a: &mut Vec<FieldElement<F>>,
    b: &mut Vec<FieldElement<F>>) -> FieldElement<F> {
        let c: FieldElement<F> = (a[0].clone() * b[0].clone()) + (a[1].clone() * b[1].clone());
        c
}

fn verify_n_2<F: lambdaworks_math::field::traits::IsField>(
    a: &mut Vec<FieldElement<F>>,
    b: &mut Vec<FieldElement<F>>,
    g: &mut Vec<FieldElement<F>>,
    h: &mut Vec<FieldElement<F>>,
    u: FieldElement<F>, x: FieldElement<F>) -> bool {

        let p: FieldElement<F> = calculate_p(a.clone(), b.clone(), g.clone(), h.clone());
        let l: FieldElement<F> = calculate_l(a.clone(), b.clone(), g.clone(), h.clone(), u.clone());
        let r: FieldElement<F> = calculate_r(a.clone(), b.clone(), g.clone(), h.clone(), u.clone());
        let a_prime: FieldElement<F> = calculate_a_prime(a.clone(), x.clone());
        let b_prime: FieldElement<F> = calculate_b_prime(b.clone(), x.clone());
        
        let c: FieldElement<F> = calculate_c(a, b);

        let left: FieldElement<F> = (x.clone().pow(2_u64)) * l.clone() + p.clone() + (c.clone() * u.clone()) + (mod_inverse(x.clone()).pow(2_u64)) * r.clone();
        let right: FieldElement<F> = (mod_inverse(x.clone()) * a_prime.clone() * g[0].clone()) + (x.clone() * a_prime.clone() * g[1].clone()) + (x.clone() * b_prime.clone() * h[0].clone()) + (mod_inverse(x) * b_prime.clone() * h[1].clone()) + (a_prime.clone() * b_prime.clone() * u.clone());

        left == right
}

fn mod_inverse<F: lambdaworks_math::field::traits::IsField>(a: FieldElement<F>) -> FieldElement<F> {
    let result: FieldElement<F> = a.inv().unwrap();
    result
}

fn inner_product<F: lambdaworks_math::field::traits::IsField>(v1: &[FieldElement<F>], v2: &[FieldElement<F>]) -> FieldElement<F> {
    if v1.len() != v2.len() {
        panic!("Vectors must have the same size.");
    }
    let result = v1.iter().zip(v2.iter()).map(|(x, y)| x * y).sum::<FieldElement<F>>();
    result
}

fn vec_low<F: lambdaworks_math::field::traits::IsField>(vector: &[FieldElement<F>]) -> Vec<FieldElement<F>> {
    let mid_idx = vector.len() / 2;
    vector[..mid_idx].to_vec()
}

fn vec_hi<F: lambdaworks_math::field::traits::IsField>(vector: &[FieldElement<F>]) -> Vec<FieldElement<F>> {
    let mid_idx = vector.len() / 2;
    vector[mid_idx..].to_vec()
}

fn computek<F: lambdaworks_math::field::traits::IsField>(vector: &[FieldElement<F>]) -> u128 {
    let mut k = 0;
    let mut length = vector.len();
    while length != 1 {
        length /= 2;
        k += 1;
    }
    k
}

fn vector_multiply_scalar<F: lambdaworks_math::field::traits::IsField>(vector: &[FieldElement<F>], s: FieldElement<F>) -> Vec<FieldElement<F>> {
    vector.iter().map(|x| x * s.clone()).collect()
}

fn vector_addition<F: lambdaworks_math::field::traits::IsField>(a: &[FieldElement<F>], b: &[FieldElement<F>]) -> Vec<FieldElement<F>> {
    if a.len() != b.len() {
        panic!("Error: Vectors must have the same length");
    }
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}

fn verify<F: lambdaworks_math::field::traits::IsField>(
    a: &mut Vec<FieldElement<F>>,
    b: &mut Vec<FieldElement<F>>,
    g: &mut Vec<FieldElement<F>>,
    h: &mut Vec<FieldElement<F>>,
    u: FieldElement<F>,
    s: FieldElement<F>,
    s_prime: FieldElement<F>,
    x: FieldElement<F>) -> bool {
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

                let a_prime = vector_addition(&vector_multiply_scalar(&a_lo, x.clone()), &vector_multiply_scalar(&a_hi, mod_inverse(x.clone())));
                let b_prime = vector_addition(&vector_multiply_scalar(&b_lo, mod_inverse(x.clone())),&vector_multiply_scalar(&b_hi, x.clone()));
                let g_prime = vector_addition(&vector_multiply_scalar(&g_lo, mod_inverse(x.clone())),&vector_multiply_scalar(&g_hi, x.clone()));
                let h_prime = vector_addition(&vector_multiply_scalar(&h_lo, x.clone()),&vector_multiply_scalar(&h_hi, mod_inverse(x.clone())));
                
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
    
    use lambdaworks_math::unsigned_integer::element::U384;
    use lambdaworks_math::field::{
        element::FieldElement,
        fields::montgomery_backed_prime_fields::{IsModulus, U384PrimeField}};

    #[derive(Clone, Debug)]
    pub struct U384ModulusP;

    impl IsModulus<U384> for U384ModulusP {
        const MODULUS: U384 = U384::from_hex_unchecked("6b");
    }

    type U384FP = U384PrimeField<U384ModulusP>;
    type U384FPElement = FieldElement<U384FP>;

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
