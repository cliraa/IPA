use lambdaworks_math::field::element::FieldElement;

// Given the coefficients of a polynomial p, a vector g, and values r and h, let's produce the commitment [p]:

fn commitment_p<F: lambdaworks_math::field::traits::IsField>(
    a: &mut Vec<FieldElement<F>>,
    g: &mut Vec<FieldElement<F>>,
    blinding_factor: FieldElement<F>,
    h: FieldElement<F>) -> FieldElement<F> {
        let result = inner_product(&a, &g) + (blinding_factor*h);
        result
    }

// Open:

// Case n = 2:

fn calculate_l<F: lambdaworks_math::field::traits::IsField>(
    a: Vec<FieldElement<F>>,
    b: Vec<FieldElement<F>>,
    g: Vec<FieldElement<F>>,
    h: FieldElement<F>,
    u: FieldElement<F>,
    s: FieldElement<F>) -> FieldElement<F> {
        let l: FieldElement<F> = (a[0].clone() * g[1].clone()) + (s.clone() * h.clone()) + (a[0].clone() * b[1].clone() * u.clone());
        l
}

fn calculate_r<F: lambdaworks_math::field::traits::IsField>(
    a: Vec<FieldElement<F>>,
    b: Vec<FieldElement<F>>,
    g: Vec<FieldElement<F>>,
    h: FieldElement<F>,
    u: FieldElement<F>,
    s_prime: FieldElement<F>) -> FieldElement<F> {
        let r: FieldElement<F> = (a[1].clone() * g[0].clone()) + (s_prime.clone() * h.clone()) + (a[1].clone() * b[0].clone() * u.clone());
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
    h: FieldElement<F>,
    blinding_factor: FieldElement<F>,
    u: FieldElement<F>,
    x: FieldElement<F>,
    s: FieldElement<F>,
    s_prime: FieldElement<F>) -> bool {

        let p: FieldElement<F> = commitment_p(a, g, blinding_factor.clone(), h.clone());
        let l: FieldElement<F> = calculate_l(a.clone(), b.clone(), g.clone(), h.clone(), u.clone(), s.clone());
        let r: FieldElement<F> = calculate_r(a.clone(), b.clone(), g.clone(), h.clone(), u.clone(), s_prime.clone());
        let a_prime: FieldElement<F> = calculate_a_prime(a.clone(), x.clone());
        let b_prime: FieldElement<F> = calculate_b_prime(b.clone(), x.clone());
        
        let c: FieldElement<F> = calculate_c(a, b);
        let r_prime: FieldElement<F> = (s.clone() * x.pow(2_u64)) + (blinding_factor.clone()) + ((mod_inverse(x.clone()).pow(2_u64)) * s_prime.clone());

        let left: FieldElement<F> = (x.clone().pow(2_u64)) * l.clone() + p.clone() + (c.clone() * u.clone()) + (mod_inverse(x.clone()).pow(2_u64)) * r.clone();
        let right: FieldElement<F> = (mod_inverse(x.clone()) * a_prime.clone() * g[0].clone()) + (x.clone() * a_prime.clone() * g[1].clone()) + (r_prime.clone() * h) + (a_prime.clone() * b_prime.clone() * u.clone());

        left == right
}

// General case n = 2^k, where k > 1

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

fn verify<F: lambdaworks_math::field::traits::IsField>( // This function must be fixed
    a: &mut Vec<FieldElement<F>>,
    b: &mut Vec<FieldElement<F>>,
    g: &mut Vec<FieldElement<F>>,
    h: FieldElement<F>,
    blinding_factor: FieldElement<F>,
    u: FieldElement<F>,
    x: FieldElement<F>,
    s: FieldElement<F>,
    s_prime: FieldElement<F>) -> bool {
        if a.len() == 2 {
            verify_n_2(a, b, g, h, blinding_factor, u, x, s, s_prime)
        } else {
            let mut k = computek(a);
            while k > 1 {
                let a_lo: Vec<FieldElement<F>> = vec_low(a);
                let a_hi: Vec<FieldElement<F>> = vec_hi(a);
                let b_lo: Vec<FieldElement<F>> = vec_low(b);
                let b_hi: Vec<FieldElement<F>> = vec_hi(b);
                let g_lo: Vec<FieldElement<F>> = vec_low(g);
                let g_hi: Vec<FieldElement<F>> = vec_hi(g);

                let commitment: FieldElement<F> = commitment_p(a, g, blinding_factor.clone(), h.clone());
                let c: FieldElement<F> = inner_product(a,b);
                let r_prime: FieldElement<F> = (s.clone() * x.pow(2_u64)) + (blinding_factor.clone()) + ((mod_inverse(x.clone()).pow(2_u64)) * s_prime.clone());
                
                let L: FieldElement<F> = inner_product(&a_lo, &g_hi) + (s.clone() * h.clone()) + inner_product(&a_lo, &b_hi) * u.clone();
                let R: FieldElement<F> = inner_product(&a_hi, &g_lo) + (s_prime.clone() * h.clone()) + inner_product(&a_hi, &b_lo) * u.clone();

                // Verification:

                let a_prime: FieldElement<F> = calculate_a_prime(a.clone(), x.clone());
                let b_prime: FieldElement<F> = calculate_b_prime(b.clone(), x.clone());

                let left: FieldElement<F> = (x.clone().pow(2_u64)) * L.clone() + commitment.clone() + (c.clone() * u.clone()) + (mod_inverse(x.clone()).pow(2_u64)) * R.clone();
                let right: FieldElement<F> = (mod_inverse(x.clone()) * a_prime.clone() * g[0].clone()) + (x.clone() * a_prime.clone() * g[1].clone()) + (r_prime.clone() * h.clone()) + (a_prime.clone() * b_prime.clone() * u.clone());
                
                let result: bool = left == right;
                assert_eq!(result, true);

                // Next step:
                
                let a_prime_vector: Vec<FieldElement<F>> = vector_addition(&vector_multiply_scalar(&a_lo, x.clone()), &vector_multiply_scalar(&a_hi, mod_inverse(x.clone())));
                let b_prime_vector: Vec<FieldElement<F>> = vector_addition(&vector_multiply_scalar(&b_lo, mod_inverse(x.clone())),&vector_multiply_scalar(&b_hi, x.clone()));
                
                let g_prime_vector: Vec<FieldElement<F>> = vector_addition(&vector_multiply_scalar(&g_lo, mod_inverse(x.clone())),&vector_multiply_scalar(&g_hi, x.clone()));
                                
                *a = a_prime_vector;
                *b = b_prime_vector;
                *g = g_prime_vector;
                k -= 1;
            };
            true
        }
    }
    
#[cfg(test)]
mod tests {

    use crate::commitment_p;
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
    fn commitment_p_test() {
        let mut p = 
            vec![U384FPElement::from_hex_unchecked("2"), U384FPElement::from_hex_unchecked("4"), U384FPElement::from_hex_unchecked("6"), U384FPElement::from_hex_unchecked("8")];
        let mut g = 
            vec![U384FPElement::from_hex_unchecked("a"), U384FPElement::from_hex_unchecked("b"), U384FPElement::from_hex_unchecked("c"), U384FPElement::from_hex_unchecked("d")];
        
        let r: U384FPElement = U384FPElement::from_hex_unchecked("f");
        let h: U384FPElement = U384FPElement::from_hex_unchecked("e");

        let expected: U384FPElement = U384FPElement::from_hex_unchecked("16");
    
        let result = commitment_p(&mut p, &mut g, r, h);
        assert_eq!(result, expected);
    }

    #[test]
    fn verify_n_2_test_1() {
        let mut a = vec![U384FPElement::from_hex_unchecked("2"), U384FPElement::from_hex_unchecked("3")];
        let mut b = vec![U384FPElement::from_hex_unchecked("1"), U384FPElement::from_hex_unchecked("6")];
        let mut g = vec![U384FPElement::from_hex_unchecked("5"), U384FPElement::from_hex_unchecked("7")];
        let h: U384FPElement = U384FPElement::from_hex_unchecked("9");
        let bf: U384FPElement = U384FPElement::from_hex_unchecked("11");
        let u: U384FPElement = U384FPElement::from_hex_unchecked("5");

        let x: U384FPElement = U384FPElement::from_hex_unchecked("6");
        let s: U384FPElement = U384FPElement::from_hex_unchecked("7");
        let s_prime: U384FPElement = U384FPElement::from_hex_unchecked("8");
    
        let result: bool = verify_n_2(&mut a, &mut b, &mut g, h, bf, u, x, s, s_prime);
        assert_eq!(result, true);
    }

    #[test]
    fn verify_n_2_test_2() {
        let mut a = vec![U384FPElement::from_hex_unchecked("21"), U384FPElement::from_hex_unchecked("13")];
        let mut b = vec![U384FPElement::from_hex_unchecked("1"), U384FPElement::from_hex_unchecked("24")];
        let mut g = vec![U384FPElement::from_hex_unchecked("10"), U384FPElement::from_hex_unchecked("17")];
        let h: U384FPElement = U384FPElement::from_hex_unchecked("15");
        let bf: U384FPElement = U384FPElement::from_hex_unchecked("7");
        let u: U384FPElement = U384FPElement::from_hex_unchecked("9");

        let x: U384FPElement = U384FPElement::from_hex_unchecked("3");
        let s: U384FPElement = U384FPElement::from_hex_unchecked("6");
        let s_prime: U384FPElement = U384FPElement::from_hex_unchecked("2");
    
        let result: bool = verify_n_2(&mut a, &mut b, &mut g, h, bf, u, x, s, s_prime);
        assert_eq!(result, true);
    }

    #[test]
    fn verify_test_3() {
        let mut a = 
            vec![U384FPElement::from_hex_unchecked("21"), U384FPElement::from_hex_unchecked("13")];
        let mut b = 
            vec![U384FPElement::from_hex_unchecked("1"), U384FPElement::from_hex_unchecked("3")];
        let mut g = 
            vec![U384FPElement::from_hex_unchecked("21"), U384FPElement::from_hex_unchecked("12")];
        let h: U384FPElement = U384FPElement::from_hex_unchecked("d");
        let bf: U384FPElement = U384FPElement::from_hex_unchecked("17");
        
        let u: U384FPElement = U384FPElement::from_hex_unchecked("f");
        let s: U384FPElement = U384FPElement::from_hex_unchecked("e");
        let s_prime: U384FPElement = U384FPElement::from_hex_unchecked("c");
        let x: U384FPElement = U384FPElement::from_hex_unchecked("d");
    
        let result = verify(&mut a, &mut b, &mut g, h, bf, u, x, s, s_prime);
        assert_eq!(result, true);
    }

    #[test]
    fn verify_test_4() { // This test is failing:
        let mut a = 
            vec![U384FPElement::from_hex_unchecked("21"), U384FPElement::from_hex_unchecked("13"), U384FPElement::from_hex_unchecked("2"),U384FPElement::from_hex_unchecked("3")];
        let mut b = 
            vec![U384FPElement::from_hex_unchecked("1"), U384FPElement::from_hex_unchecked("3"), U384FPElement::from_hex_unchecked("9"),U384FPElement::from_hex_unchecked("1b")];
        let mut g = 
            vec![U384FPElement::from_hex_unchecked("21"), U384FPElement::from_hex_unchecked("12"), U384FPElement::from_hex_unchecked("3"),U384FPElement::from_hex_unchecked("6")];
        let h: U384FPElement = U384FPElement::from_hex_unchecked("d");
        let bf: U384FPElement = U384FPElement::from_hex_unchecked("17");
        
        let u: U384FPElement = U384FPElement::from_hex_unchecked("f");
        let s: U384FPElement = U384FPElement::from_hex_unchecked("e");
        let s_prime: U384FPElement = U384FPElement::from_hex_unchecked("c");
        let x: U384FPElement = U384FPElement::from_hex_unchecked("d");
    
        let result = verify(&mut a, &mut b, &mut g, h, bf, u, x, s, s_prime);
        assert_eq!(result, true);
    }
}
