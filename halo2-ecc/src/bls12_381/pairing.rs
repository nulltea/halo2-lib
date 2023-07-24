#![allow(non_snake_case)]
use super::{Fp12Chip, Fp2Chip, FpChip, FpPoint, Fq, FqPoint};
use crate::ecc;
use crate::fields::vector::{FieldVector, FieldVectorChip};
use crate::{
    ecc::{EcPoint, EccChip},
    fields::fp12::mul_no_carry_w6,
    fields::{FieldChip, PrimeField},
};
use halo2_base::utils::ScalarField;
use halo2_base::Context;
use halo2curves::bls12_381::{Fq2, BLS_X_IS_NEGATIVE};
use halo2curves::{
    bls12_381::{BitIterator, Fq12, G1Affine, G2Affine, BLS_X, FROBENIUS_COEFF_FQ12_C1},
    bn256::SIX_U_PLUS_2_NAF,
};
use itertools::Itertools;
use rayon::vec;

const XI_0: i64 = 1;

// Inputs:
//  Q0 = (x_1, y_1) and Q1 = (x_2, y_2) are points in E(Fp2)
//  P is point (X, Y) in E(Fp)
// Assuming Q0 != Q1
// Output:
//  line_{Psi(Q0), Psi(Q1)}(P) where Psi(x,y) = (w^2 x, w^3 y)
//  - equals w^3 (y_1 - y_2) X + w^2 (x_2 - x_1) Y + w^5 (x_1 y_2 - x_2 y_1) =: out3 * w^3 + out2 * w^2 + out5 * w^5 where out2, out3, out5 are Fp2 points
// Output is [None, None, out2, out3, None, out5] as vector of `Option<FqPoint>`s
pub fn sparse_line_function_unequal<F: PrimeField>(
    fp2_chip: &Fp2Chip<F>,
    ctx: &mut Context<F>,
    Q: (&EcPoint<F, FqPoint<F>>, &EcPoint<F, FqPoint<F>>),
    P: &EcPoint<F, FpPoint<F>>,
) -> Vec<Option<FqPoint<F>>> {
    let (x_1, y_1) = (&Q.0.x, &Q.0.y);
    let (x_2, y_2) = (&Q.1.x, &Q.1.y);
    let (X, Y) = (&P.x, &P.y);
    assert_eq!(x_1.0.len(), 2);
    assert_eq!(y_1.0.len(), 2);
    assert_eq!(x_2.0.len(), 2);
    assert_eq!(y_2.0.len(), 2);

    let y1_minus_y2 = fp2_chip.sub_no_carry(ctx, y_1, y_2);
    let x2_minus_x1 = fp2_chip.sub_no_carry(ctx, x_2, x_1);
    let x1y2 = fp2_chip.mul_no_carry(ctx, x_1, y_2);
    let x2y1 = fp2_chip.mul_no_carry(ctx, x_2, y_1);

    let out3 = fp2_chip.0.fp_mul_no_carry(ctx, y1_minus_y2, X);
    let out4 = fp2_chip.0.fp_mul_no_carry(ctx, x2_minus_x1, Y);
    let out2 = fp2_chip.sub_no_carry(ctx, &x1y2, &x2y1);

    // so far we have not "carried mod p" for any of the outputs
    // we do this below
    [None, Some(out2), None, Some(out3), Some(out4), None]
        .into_iter()
        .map(|option_nc| option_nc.map(|nocarry| fp2_chip.carry_mod(ctx, nocarry)))
        .collect()
}

// Assuming curve is of form Y^2 = X^3 + b (a = 0) to save operations
// Inputs:
//  Q = (x, y) is a point in E(Fp)
//  P = (P.x, P.y) in E(Fp2)
// Output:
//  line_{Psi(Q), Psi(Q)}(P) where Psi(x,y) = (w^2 x, w^3 y)
//  - equals (3x^3 - 2y^2)(XI_0 + u) + w^4 (-3 x^2 * Q.x) + w^3 (2 y * Q.y) =: out0 + out4 * w^4 + out3 * w^3 where out0, out3, out4 are Fp2 points
// Output is [out0, None, out2, out3, None, None] as vector of `Option<FqPoint>`s
pub fn sparse_line_function_equal<F: PrimeField>(
    fp2_chip: &Fp2Chip<F>,
    ctx: &mut Context<F>,
    Q: &EcPoint<F, FpPoint<F>>,
    P: &EcPoint<F, FqPoint<F>>,
) -> Vec<Option<FqPoint<F>>> {
    let (x, y) = (&P.x, &P.y);
    assert_eq!(x.0.len(), 2);
    assert_eq!(y.0.len(), 2);

    let x_sq = fp2_chip.mul(ctx, x, x);

    let x_cube = fp2_chip.mul_no_carry(ctx, &x_sq, x);
    let three_x_cu = fp2_chip.scalar_mul_no_carry(ctx, &x_cube, 3);
    let y_sq = fp2_chip.mul_no_carry(ctx, y, y);
    let two_y_sq = fp2_chip.scalar_mul_no_carry(ctx, &y_sq, 2);
    let out0 = fp2_chip.sub_no_carry(ctx, &three_x_cu, &two_y_sq);

    let x_sq_Px = fp2_chip.0.fp_mul_no_carry(ctx, x_sq, &Q.x);
    let out2 = fp2_chip.scalar_mul_no_carry(ctx, x_sq_Px, -3);

    let y_Py = fp2_chip.0.fp_mul_no_carry(ctx, y.clone(), &Q.y);
    let out3 = fp2_chip.scalar_mul_no_carry(ctx, &y_Py, 2);

    // so far we have not "carried mod p" for any of the outputs
    // we do this below
    [Some(out0), None, Some(out2), Some(out3), None, None]
        .into_iter()
        .map(|option_nc| option_nc.map(|nocarry| fp2_chip.carry_mod(ctx, nocarry)))
        .collect()
}

// multiply Fp12 point `a` with Fp12 point `b` where `b` is len 6 vector of Fp2 points, where some are `None` to represent zero.
// Assumes `b` is not vector of all `None`s
pub fn sparse_fp12_multiply<F: PrimeField>(
    fp2_chip: &Fp2Chip<F>,
    ctx: &mut Context<F>,
    a: &FqPoint<F>,
    b_fp2_coeffs: &[Option<FqPoint<F>>],
) -> FqPoint<F> {
    assert_eq!(a.0.len(), 12);
    assert_eq!(b_fp2_coeffs.len(), 6);
    let mut a_fp2_coeffs = Vec::with_capacity(6);
    for i in 0..6 {
        a_fp2_coeffs.push(FieldVector(vec![a[i].clone(), a[i + 6].clone()]));
    }
    // a * b as element of Fp2[w] without evaluating w^6 = (XI_0 + u)
    let mut prod_2d = vec![None; 11];
    for i in 0..6 {
        for j in 0..6 {
            prod_2d[i + j] =
                match (prod_2d[i + j].clone(), &a_fp2_coeffs[i], b_fp2_coeffs[j].as_ref()) {
                    (a, _, None) => a,
                    (None, a, Some(b)) => {
                        let ab = fp2_chip.mul_no_carry(ctx, a, b);
                        Some(ab)
                    }
                    (Some(a), b, Some(c)) => {
                        let bc = fp2_chip.mul_no_carry(ctx, b, c);
                        let out = fp2_chip.add_no_carry(ctx, &a, &bc);
                        Some(out)
                    }
                };
        }
    }

    let mut out_fp2 = Vec::with_capacity(6);
    for i in 0..6 {
        // prod_2d[i] + prod_2d[i+6] * w^6
        let prod_nocarry = if i != 5 {
            let eval_w6 = prod_2d[i + 6]
                .as_ref()
                .map(|a| mul_no_carry_w6::<_, _, XI_0>(fp2_chip.fp_chip(), ctx, a.clone()));
            match (prod_2d[i].as_ref(), eval_w6) {
                (None, b) => b.unwrap(), // Our current use cases of 235 and 034 sparse multiplication always result in non-None value
                (Some(a), None) => a.clone(),
                (Some(a), Some(b)) => fp2_chip.add_no_carry(ctx, a, &b),
            }
        } else {
            prod_2d[i].clone().unwrap()
        };
        let prod = fp2_chip.carry_mod(ctx, prod_nocarry);
        out_fp2.push(prod);
    }

    let mut out_coeffs = Vec::with_capacity(12);
    for fp2_coeff in &out_fp2 {
        out_coeffs.push(fp2_coeff[0].clone());
    }
    for fp2_coeff in &out_fp2 {
        out_coeffs.push(fp2_coeff[1].clone());
    }
    FieldVector(out_coeffs)
}

// Input:
// - g is Fp12 point
// - P = (P0, P1) with Q0, Q1 points in E(Fp2)
// - Q is point in E(Fp)
// Output:
// - out = g * l_{Psi(Q0), Psi(Q1)}(P) as Fp12 point
pub fn fp12_multiply_with_line_unequal<F: PrimeField>(
    fp2_chip: &Fp2Chip<F>,
    ctx: &mut Context<F>,
    g: &FqPoint<F>,
    P: (&EcPoint<F, FqPoint<F>>, &EcPoint<F, FqPoint<F>>),
    Q: &EcPoint<F, FpPoint<F>>,
    debug: bool,
) -> FqPoint<F> {
    let line = sparse_line_function_unequal::<F>(fp2_chip, ctx, P, Q);
    // if debug {
    //     println!(
    //         "line unequal: {:#?}",
    //         line.iter()
    //             .map(|lf| lf.clone().map(|lfv| fp2_chip.get_assigned_value(&lfv.into())))
    //             .collect_vec()
    //     );
    // }
    sparse_fp12_multiply::<F>(fp2_chip, ctx, g, &line)
}

// Input:
// - g is Fp12 point
// - P is point in E(Fp2)
// - Q is point in E(Fp)
// Output:
// - out = g * l_{Psi(Q), Psi(Q)}(P) as Fp12 point
pub fn fp12_multiply_with_line_equal<F: PrimeField>(
    fp2_chip: &Fp2Chip<F>,
    ctx: &mut Context<F>,
    g: &FqPoint<F>,
    P: &EcPoint<F, FqPoint<F>>,
    Q: &EcPoint<F, FpPoint<F>>,
    debug: bool,
) -> FqPoint<F> {
    let line = sparse_line_function_equal::<F>(fp2_chip, ctx, Q, P);
    if debug {
        println!(
            "line: {:#?}",
            line.iter()
                .map(|lf| lf.clone().map(|lfv| fp2_chip.get_assigned_value(&lfv.into())))
                .collect_vec()
        );
    }
    sparse_fp12_multiply::<F>(fp2_chip, ctx, g, &line)
}

// Assuming curve is of form `y^2 = x^3 + b` for now (a = 0) for less operations
// Value of `b` is never used
// Inputs:
// - Q = (x, y) is a point in E(Fp2)
// - P is a point in E(Fp)
// - `pseudo_binary_encoding` is fixed vector consisting of {-1, 0, 1} entries such that `loop_count = sum pseudo_binary_encoding[i] * 2^i`
// Output:
//  - f_{loop_count}(Q,P) * l_{[loop_count] Q', Frob_p(Q')}(P) * l_{[loop_count] Q' + Frob_p(Q'), -Frob_p^2(Q')}(P)
//  - where we start with `f_1(Q,P) = 1` and use Miller's algorithm f_{i+j} = f_i * f_j * l_{i,j}(Q,P)
//  - Q' = Psi(Q) in E(Fp12)
//  - Frob_p(x,y) = (x^p, y^p)
//  - Above formula is specific to BN curves
// Assume:
//  - Q != O and the order of Q in E(Fp2) is r
//  - r is prime, so [i]Q != [j]Q for i != j in Z/r
//  - `0 <= loop_count < r` and `loop_count < p` (to avoid [loop_count]Q' = Frob_p(Q'))
//  - x^3 + b = 0 has no solution in Fp2, i.e., the y-coordinate of Q cannot be 0.
pub fn miller_loop<F: PrimeField>(
    ecc_chip: &EccChip<F, Fp2Chip<F>>,
    ctx: &mut Context<F>,
    Q: &EcPoint<F, FqPoint<F>>,
    P: &EcPoint<F, FpPoint<F>>,
    pseudo_binary_encoding: &[i8],
) -> FqPoint<F> {
    unimplemented!()
}

// let pairs = [(a_i, b_i)], a_i in G_1, b_i in G_2
// output is Prod_i e'(a_i, b_i), where e'(a_i, b_i) is the output of `miller_loop_BN(b_i, a_i)`
pub fn multi_miller_loop<F: PrimeField>(
    ecc_chip: &EccChip<F, Fp2Chip<F>>,
    ctx: &mut Context<F>,
    pairs: Vec<(&EcPoint<F, FpPoint<F>>, &EcPoint<F, FqPoint<F>>)>,
) -> FqPoint<F> {
    println!("---------- multi_miller_loop_BN ----------");
    println!(
        "[circuit] P1: ({:?} {:?})",
        ecc_chip.field_chip().get_assigned_value(&pairs[0].1.x.clone().into()),
        ecc_chip.field_chip().get_assigned_value(&pairs[0].1.y.clone().into())
    );
    let fp_chip = ecc_chip.field_chip.fp_chip();
    let fp12_chip = Fp12Chip::<F>::new(fp_chip);

    // let mut f = fp12_chip.load_private(ctx, Fq12::one());
    // initialize the first line function into Fq12 point
    let mut f = {
        let sparse_f =
            sparse_line_function_equal::<F>(ecc_chip.field_chip(), ctx, pairs[0].0, pairs[0].1);
        assert_eq!(sparse_f.len(), 6);

        let zero_fp = fp_chip.load_constant(ctx, Fq::zero());
        let mut f_coeffs = Vec::with_capacity(12);
        for coeff in &sparse_f {
            if let Some(fp2_point) = coeff {
                f_coeffs.push(fp2_point[0].clone());
            } else {
                f_coeffs.push(zero_fp.clone());
            }
        }
        for coeff in &sparse_f {
            if let Some(fp2_point) = coeff {
                f_coeffs.push(fp2_point[1].clone());
            } else {
                f_coeffs.push(zero_fp.clone());
            }
        }
        FieldVector(f_coeffs) // this is just to skip first mul leveraging the fact that f = 1
    };
    println!(
        "[circuit] f * line equal 0 = {:#?}",
        fp12_chip.get_assigned_value(&f.clone().into())
    );
    for &(q, p) in pairs.iter().skip(1) {
        f = fp12_multiply_with_line_equal::<F>(ecc_chip.field_chip(), ctx, &f, p, q, true);
    }

    println!(
        "[circuit] f * line equal 1 = {:#?}",
        fp12_chip.get_assigned_value(&f.clone().into())
    );

    let fp12_chip = Fp12Chip::<F>::new(fp_chip);

    const DBG_N: usize = 2;

    let mut r = pairs.iter().map(|pair| pair.1.clone()).collect::<Vec<_>>();
    let mut found_one = true;
    let mut j = 0;
    let mut prevBit = true;

    for r in r.iter_mut() {
        *r = ecc_chip.double(ctx, r.clone());
    }
    for (i, bit) in (0..62).rev().map(|i| (i as usize, ((BLS_X >> i) & 1) == 1)) {
        if prevBit {
            for (r, &(q, p)) in r.iter_mut().zip(pairs.iter()) {
                f = fp12_multiply_with_line_unequal::<F>(
                    ecc_chip.field_chip(),
                    ctx,
                    &f,
                    (r, p),
                    q,
                    j < DBG_N,
                );
                *r = ecc_chip.add_unequal(ctx, r.clone(), p.clone(), false);
                if j < DBG_N {
                    println!(
                        "[circuit] f * line unequal = {:#?}",
                        fp12_chip.get_assigned_value(&f.clone().into())
                    );
                    println!(
                        "[circuit] add(r): ({:?} {:?})",
                        ecc_chip.field_chip().get_assigned_value(&r.x.clone().into()),
                        ecc_chip.field_chip().get_assigned_value(&r.y.clone().into())
                    );
                    println!("---------pair");
                }
            }

            if j < DBG_N {
                println!("---------------add");
            }
        }

        prevBit = bit;

        println!("---------------bit {} = {}", i, bit);

        if !found_one {
            println!("---------------skipping");
            found_one = bit;
            continue;
        }
        

        f = fp12_chip.mul(ctx, &f, &f);

        for (r, &(q, p)) in r.iter_mut().zip(pairs.iter()) {
            f = fp12_multiply_with_line_equal::<F>(ecc_chip.field_chip(), ctx, &f, r, q, j < DBG_N);
            *r = ecc_chip.double(ctx, r.clone());
            if j < DBG_N {
                println!(
                    "[circuit] double(r): ({:?} {:?})",
                    ecc_chip.field_chip().get_assigned_value(&r.x.clone().into()),
                    ecc_chip.field_chip().get_assigned_value(&r.y.clone().into())
                );
                println!(
                    "[circuit] f * line equal = {:#?}",
                    fp12_chip.get_assigned_value(&f.clone().into())
                );

                println!("---------pair");
            }
        }
        if j < DBG_N {
            println!("---------------double");
        }

        

        j += 1
    }

    f
}

// let pairs = [(a_i, b_i)], a_i in G_1, b_i in G_2
// output is Prod_i e'(a_i, b_i), where e'(a_i, b_i) is the output of `miller_loop_BN(b_i, a_i)`
pub fn multi_miller_loop_BLS<F: PrimeField>(
    ecc_chip: &EccChip<F, Fp2Chip<F>>,
    ctx: &mut Context<F>,
    pairs: Vec<(&EcPoint<F, FpPoint<F>>, &EcPoint<F, FqPoint<F>>)>,
) -> FqPoint<F> {
    let fp_chip = ecc_chip.field_chip.fp_chip();

    let fp12_chip = Fp12Chip::<F>::new(fp_chip);
    let mut f = fp12_chip.load_private(ctx, Fq12::one());

    let z_one = ecc_chip.field_chip().load_constant(ctx, Fq2::one());
    let mut coeffs = vec![];

    let mut r = pairs
        .iter()
        .map(|pair| EcProjective::from_affine(pair.1.clone(), ecc_chip.field_chip(), ctx))
        .collect::<Vec<_>>();
    let mut found_one = false;
    let mut j = 0;
    for bit in (0..64).rev().map(|b| (((BLS_X >> 1) >> b) & 1) == 1) {
        if !found_one {
            found_one = bit;
            continue;
        }

        for (r, &(a, _)) in r.iter_mut().zip(pairs.iter()) {
            coeffs = doubling_step::<F>(ecc_chip.field_chip(), ctx, r, true);
            f = ell::<F>(f, &coeffs, a, &fp12_chip, ctx);
        }

        if bit {
            for (r, &(a, b)) in r.iter_mut().zip(pairs.iter()) {
                coeffs = addition_step::<F>(ecc_chip.field_chip(), ctx, r, b, true);
                f = ell::<F>(f, &coeffs, a, &fp12_chip, ctx);
            }
        }

        f = fp12_chip.mul(ctx, &f, &f);

        j += 1;
    }

    println!("[circuit] f -> loop = {:?}", fp12_chip.get_assigned_value(&f.clone().into()));

    for (r, &(a, _)) in r.iter_mut().zip(pairs.iter()) {
        coeffs = doubling_step::<F>(ecc_chip.field_chip(), ctx, r, true);
        f = ell::<F>(f, &coeffs, a, &fp12_chip, ctx);
    }

    if BLS_X_IS_NEGATIVE {
        f = fp12_chip.conjugate(ctx, f);
    }

    println!("[circuit] f -> end = {:?}", fp12_chip.get_assigned_value(&f.clone().into()));

    f
}

fn ell<'chip, F: PrimeField>(
    f: FqPoint<F>,
    coeffs: &Vec<FqPoint<F>>,
    p: &EcPoint<F, FpPoint<F>>,
    fp12_chip: &Fp12Chip<'chip, F>,
    ctx: &mut Context<F>,
) -> FqPoint<F> {
    let c00 = &coeffs[0].0[0];
    let c01 = &coeffs[0].0[1];
    let c10 = &coeffs[1].0[0];
    let c11 = &coeffs[1].0[1];

    let c00 = fp12_chip.fp_chip().mul(ctx, c00, p.y.clone());
    let c01 = fp12_chip.fp_chip().mul(ctx, c01, p.y.clone());
    let c10 = fp12_chip.fp_chip().mul(ctx, c10, p.x.clone());
    let c11 = fp12_chip.fp_chip().mul(ctx, c11, p.x.clone());

    let fp2_chip = Fp2Chip::<F>::new(fp12_chip.fp_chip());

    fq12_mul_by_014(
        f,
        coeffs[2].clone(),
        FieldVector(vec![c10, c11]),
        FieldVector(vec![c00, c01]),
        &fp2_chip,
        ctx,
    )
}

type Fp6Chip<'chip, F> = FieldVectorChip<'chip, F, FpChip<'chip, F>>;

fn fq12_mul_by_014<'chip, F: PrimeField>(
    f: FqPoint<F>,
    c0: FqPoint<F>,
    c1: FqPoint<F>,
    c4: FqPoint<F>,
    fp2_chip: &Fp2Chip<'chip, F>,
    ctx: &mut Context<F>,
) -> FqPoint<F> {
    let fp6_chip = Fp6Chip::<F>::new(fp2_chip.fp_chip());

    let fc0 = FieldVector(permute_vector(&f.0, &[0, 6, 2, 8, 4, 10]));
    let fc1 = FieldVector(permute_vector(&f.0, &[1, 7, 3, 9, 5, 11]));

    let t0 = fq6_mul_by_01(fc0.clone(), c0.clone(), c1.clone(), fp2_chip, ctx);
    let t1 = fq6_mul_by_1(fc1.clone(), c4.clone(), fp2_chip, ctx);
    let o = fp2_chip.add(ctx, c1, c4);

    let x0 = fq6_mul_by_nonresidue(&t1, fp2_chip, ctx);
    let x0 = fp6_chip.add_no_carry(ctx, x0, t0.clone());
    let x0 = fp6_chip.carry_mod(ctx, x0);

    let x1 = fp6_chip.add_no_carry(ctx, fc0.clone(), fc1.clone());
    let x1 = fp6_chip.carry_mod(ctx, x1);
    let x1 = fq6_mul_by_01(x1, c0, o, fp2_chip, ctx);
    let x1 = fp6_chip.sub_no_carry(ctx, x1, t0);
    let x1 = fp6_chip.sub_no_carry(ctx, x1, t1);
    let x1 = fp6_chip.carry_mod(ctx, x1);

    let c = x0.0.clone().into_iter().chain(x1.0.clone()).collect_vec();

    FieldVector(
        permute_vector(&c, &[0, 6, 2, 8, 4, 10])
            .into_iter()
            .chain(permute_vector(&c, &[1, 7, 3, 9, 5, 11]))
            .collect_vec(),
    )
}

pub fn permute_vector<T: Clone>(v1: &Vec<T>, indexes: &[usize]) -> Vec<T> {
    indexes.iter().map(|i| v1[*i].clone()).collect()
}

pub fn fq6_mul<'chip, F: ScalarField>(
    fp2_chip: &Fp2Chip<'chip, F>,
    ctx: &mut Context<F>,
    a: &FqPoint<F>, // fq6
    b: &FqPoint<F>, // fq6
) -> FqPoint<F> {
    let a = (
        FieldVector(a.0[..2].to_vec()),
        FieldVector(a.0[2..4].to_vec()),
        FieldVector(a.0[4..6].to_vec()),
    );
    let b = (
        FieldVector(b.0[..2].to_vec()),
        FieldVector(b.0[2..4].to_vec()),
        FieldVector(b.0[4..6].to_vec()),
    );
    let ab00 = fp2_chip.mul(ctx, &a.0, &b.0);
    let ab11 = fp2_chip.mul(ctx, &a.1, &b.1);
    let ab22 = fp2_chip.mul(ctx, &a.2, &b.2);

    let c0 = {
        let b12 = fp2_chip.add(ctx, b.1.clone(), b.2.clone());
        let a12 = fp2_chip.add(ctx, a.1.clone(), a.2.clone());
        let t = fp2_chip.mul(ctx, a12, b12);
        let t = fp2_chip.sub_no_carry(ctx, t, ab11.clone());
        let t = fp2_chip.sub_no_carry(ctx, t, ab22.clone());
        let t = fp2_chip.carry_mod(ctx, t);
        let t = fq2_mul_by_nonresidue(&t, fp2_chip.fp_chip(), ctx);
        fp2_chip.add(ctx, t, ab00.clone())
    };

    let c1 = {
        let b01 = fp2_chip.add(ctx, b.0.clone(), b.1.clone());
        let a01 = fp2_chip.add(ctx, a.0.clone(), a.1.clone());
        let t = fp2_chip.mul(ctx, a01, b01);
        let t = fp2_chip.sub_no_carry(ctx, t, ab00.clone());
        let t = fp2_chip.sub_no_carry(ctx, t, ab11.clone());
        let ab22 = fq2_mul_by_nonresidue(&ab22, fp2_chip.fp_chip(), ctx);
        fp2_chip.add(ctx, t, &ab22)
    };

    let c2 = {
        let b02 = fp2_chip.add(ctx, b.0.clone(), b.2.clone());
        let a02 = fp2_chip.add(ctx, a.0.clone(), a.2.clone());
        let t = fp2_chip.mul(ctx, a02, b02);
        let t = fp2_chip.sub_no_carry(ctx, t, ab00);
        let t = fp2_chip.add_no_carry(ctx, t, ab11);
        fp2_chip.sub(ctx, t, ab22)
    };

    FieldVector(c0.0.into_iter().chain(c1.0).chain(c2.0).collect())
}

fn fq6_mul_by_01<'chip, F: PrimeField>(
    a: FqPoint<F>,  // fq6
    b0: FqPoint<F>, // fq2
    b1: FqPoint<F>, // fq2
    fp2_chip: &Fp2Chip<'chip, F>,
    ctx: &mut Context<F>,
) -> FqPoint<F> {
    let ab00 = fp2_chip.mul(ctx, FieldVector(a.0[..2].to_vec()), b0.clone());
    let ab11 = fp2_chip.mul(ctx, FieldVector(a.0[2..4].to_vec()), b1.clone());

    let c0 = {
        let b12 = b1.clone();
        let a12 =
            fp2_chip.add(ctx, FieldVector(a.0[2..4].to_vec()), FieldVector(a.0[4..6].to_vec()));
        let t = fp2_chip.mul(ctx, a12, b12);
        let t = fp2_chip.sub_no_carry(ctx, t, ab11.clone());
        let t = fp2_chip.carry_mod(ctx, t);
        let t = fq2_mul_by_nonresidue(&t, fp2_chip.fp_chip(), ctx);
        fp2_chip.add(ctx, t, ab00.clone())
    };

    let c1 = {
        let b01 = fp2_chip.add(ctx, b0.clone(), b1);
        let a01 =
            fp2_chip.add(ctx, FieldVector(a.0[..2].to_vec()), FieldVector(a.0[2..4].to_vec()));
        let t = fp2_chip.mul(ctx, a01, b01);
        let t = fp2_chip.sub_no_carry(ctx, t, ab00.clone());
        fp2_chip.sub(ctx, t, ab11.clone())
    };

    let c2 = {
        let b02 = b0.clone();
        let a02 =
            fp2_chip.add(ctx, FieldVector(a.0[..2].to_vec()), FieldVector(a.0[4..6].to_vec()));
        let t = fp2_chip.mul(ctx, a02, b02);
        let t = fp2_chip.sub_no_carry(ctx, t, ab00);
        fp2_chip.add(ctx, t, ab11)
    };

    FieldVector(c0.0.into_iter().chain(c1.0).chain(c2.0).collect())
}

fn fq6_mul_by_1<'chip, F: PrimeField>(
    a: FqPoint<F>,  // fq6
    b1: FqPoint<F>, // fq2
    fp2_chip: &Fp2Chip<'chip, F>,
    ctx: &mut Context<F>,
) -> FqPoint<F> {
    let ab11 = fp2_chip.mul(ctx, FieldVector(a.0[2..4].to_vec()), b1.clone());

    let c0 = {
        let b12 = b1.clone();
        let a12 =
            fp2_chip.add(ctx, FieldVector(a.0[2..4].to_vec()), FieldVector(a.0[4..6].to_vec()));
        let t = fp2_chip.mul(ctx, a12, b12);
        let t = fp2_chip.sub(ctx, t, ab11.clone());
        fq2_mul_by_nonresidue(&t, fp2_chip.fp_chip(), ctx)
    };

    let c1 = {
        let b01 = b1;
        let a01 =
            fp2_chip.add(ctx, FieldVector(a.0[..2].to_vec()), FieldVector(a.0[2..4].to_vec()));
        let t = fp2_chip.mul(ctx, a01, b01);
        fp2_chip.sub(ctx, t, ab11.clone())
    };

    let c2 = ab11;

    FieldVector(c0.0.into_iter().chain(c1.0).chain(c2.0).collect())
}

pub fn fq2_mul_by_nonresidue<'chip, F: PrimeField>(
    a: &FqPoint<F>,
    fp_chip: &FpChip<'chip, F>,
    ctx: &mut Context<F>,
) -> FqPoint<F> {
    let c0 = fp_chip.sub(ctx, &a.0[0], &a.0[1]);
    let c1 = fp_chip.add(ctx, &a.0[0], &a.0[1]);

    FieldVector(vec![c0, c1])
}

pub fn fq6_mul_by_nonresidue<'chip, F: PrimeField>(
    a: &FqPoint<F>,
    fp2_chip: &Fp2Chip<'chip, F>,
    ctx: &mut Context<F>,
) -> FqPoint<F> {
    let c0 = fq2_mul_by_nonresidue(&FieldVector(a.0[4..6].to_vec()), fp2_chip.fp_chip(), ctx);
    FieldVector(c0.0.into_iter().chain(a.0[..2].to_vec()).chain(a.0[2..4].to_vec()).collect())
}

struct EcProjective<F: PrimeField> {
    x: FqPoint<F>,
    y: FqPoint<F>,
    z: FqPoint<F>,
}

impl<F: PrimeField> EcProjective<F> {
    fn from_affine<'chip>(
        affine: EcPoint<F, FqPoint<F>>,
        fp2_chip: &Fp2Chip<'chip, F>,
        ctx: &mut Context<F>,
    ) -> Self {
        let z = fp2_chip.load_constant(ctx, Fq2::one());
        Self { x: affine.x, y: affine.y, z }
    }

    fn to_affine<'chip>(
        &self,
        fp2_chip: &Fp2Chip<'chip, F>,
        ctx: &mut Context<F>,
    ) -> EcPoint<F, FqPoint<F>> {
        let x = fp2_chip.divide_unsafe(ctx, &self.x, &self.z);
        let y = fp2_chip.divide_unsafe(ctx, &self.y, &self.z);
        EcPoint::new(x, y)
    }
}

fn addition_step<'chip, F: PrimeField>(
    fp2_chip: &Fp2Chip<'chip, F>,
    ctx: &mut Context<F>,
    r: &mut EcProjective<F>,
    p: &EcPoint<F, FqPoint<F>>,
    with_coeffs: bool,
) -> Vec<FqPoint<F>> {
    let zsq = fp2_chip.mul(ctx, &r.z, &r.z);
    let ysq = fp2_chip.mul(ctx, p.y(), p.y());
    let tv0 = fp2_chip.mul(ctx, &zsq, p.x());
    let tv1 = {
        let tv = fp2_chip.add_no_carry(ctx, p.y(), &r.z);
        let tv = fp2_chip.carry_mod(ctx, tv);
        let tv = fp2_chip.mul(ctx, &tv, &tv);
        let tv = fp2_chip.sub_no_carry(ctx, tv, &ysq);
        let tv = fp2_chip.sub_no_carry(ctx, tv, &zsq);
        let tv = fp2_chip.carry_mod(ctx, tv);
        fp2_chip.mul(ctx, &tv, &zsq)
    };
    let tv2 = fp2_chip.sub_no_carry(ctx, &tv0, &r.x);
    let tv3 = fp2_chip.mul(ctx, &tv2, &tv2);
    let tv4 = fp2_chip.add_no_carry(ctx, &tv3, &tv3);
    let tv4 = fp2_chip.add_no_carry(ctx, &tv4, &tv4);
    let tv4 = fp2_chip.carry_mod(ctx, tv4);
    let tv5 = fp2_chip.mul(ctx, &tv4, &tv2);
    let tv6 = {
        let tv = fp2_chip.sub_no_carry(ctx, &tv1, &r.y);
        let tv = fp2_chip.sub_no_carry(ctx, tv, &r.y);
        fp2_chip.carry_mod(ctx, tv)
    };
    let tv9 = fp2_chip.mul(ctx, &tv6, p.x());
    let tv7 = fp2_chip.mul(ctx, &tv4, &r.x);
    let xout = {
        let tv = fp2_chip.mul(ctx, &tv6, &tv6);
        let tv = fp2_chip.sub_no_carry(ctx, tv, &tv5);
        let tv = fp2_chip.sub_no_carry(ctx, tv, &tv7);
        let tv = fp2_chip.sub_no_carry(ctx, tv, &tv7);
        fp2_chip.carry_mod(ctx, tv)
    };
    let zout = {
        let tv = fp2_chip.add_no_carry(ctx, &r.z, &tv2);
        let tv = fp2_chip.carry_mod(ctx, tv);
        let tv = fp2_chip.mul(ctx, &tv, &tv);
        let tv = fp2_chip.sub_no_carry(ctx, tv, &zsq);
        let tv = fp2_chip.sub_no_carry(ctx, tv, &tv3);
        fp2_chip.carry_mod(ctx, tv)
    };
    let tv10 = fp2_chip.add_no_carry(ctx, p.y(), &zout);
    let tv8 = {
        let tv = fp2_chip.sub_no_carry(ctx, &tv7, &xout);
        let tv = fp2_chip.carry_mod(ctx, tv);
        fp2_chip.mul(ctx, &tv, &tv6)
    };
    let tv0 = fp2_chip.mul(ctx, &r.y, &tv5);
    let tv0 = fp2_chip.add_no_carry(ctx, &tv0, &tv0);
    let yout = fp2_chip.sub(ctx, &tv8, &tv0);

    *r = EcProjective { x: xout.clone(), y: yout, z: zout.clone() };

    if !with_coeffs {
        return vec![];
    }

    let tv10 = {
        let tv = fp2_chip.mul(ctx, &tv10, &tv10);
        fp2_chip.sub(ctx, tv, &ysq)
    };
    let zsq = fp2_chip.mul(ctx, &r.z, &r.z);
    let tv10 = fp2_chip.sub_no_carry(ctx, &tv10, &zsq);
    let tv9 = fp2_chip.add_no_carry(ctx, &tv9, &tv9);
    let tv9 = fp2_chip.sub(ctx, &tv9, &tv10);
    let tv10 = fp2_chip.add(ctx, &r.z, &r.z);
    let tv6 = fp2_chip.negate(ctx, tv6);
    let tv1 = fp2_chip.add(ctx, &tv6, &tv6);

    vec![tv10, tv1, tv9]
}

fn doubling_step<'chip, F: PrimeField>(
    fp2_chip: &Fp2Chip<'chip, F>,
    ctx: &mut Context<F>,
    P: &mut EcProjective<F>,
    with_coeffs: bool,
) -> Vec<FqPoint<F>> {
    let x = P.x.clone();
    let y = P.y.clone();
    let z = P.z.clone();

    let tv0 = fp2_chip.mul(ctx, &x, &x);
    let tv1 = fp2_chip.mul(ctx, &y, &y);
    let tv2 = fp2_chip.mul(ctx, &tv1, &tv1);
    let tv3 = {
        let tv = fp2_chip.add_no_carry(ctx, &tv1, &x);
        let tv = fp2_chip.mul(ctx, &tv, &tv);
        let tv = fp2_chip.sub_no_carry(ctx, &tv, &tv0);
        let tv = fp2_chip.sub_no_carry(ctx, &tv, &tv2);
        fp2_chip.carry_mod(ctx, tv)
    };
    let tv3 = fp2_chip.add_no_carry(ctx, &tv3, &tv3);
    let tv4 = {
        let tv = fp2_chip.add_no_carry(ctx, &tv0, &tv0);
        let tv = fp2_chip.add_no_carry(ctx, tv, &tv0);
        fp2_chip.carry_mod(ctx, tv)
    };
    let tv6 = fp2_chip.add(ctx, &x, &tv4);
    let tv5 = fp2_chip.mul(ctx, &tv4, &tv4);
    let zsq = fp2_chip.mul(ctx, &z, &z);
    let xout = {
        let tv = fp2_chip.sub_no_carry(ctx, &tv5, &tv3);
        let tv = fp2_chip.sub_no_carry(ctx, &tv, &tv3);
        fp2_chip.carry_mod(ctx, tv)
    };
    let zout = {
        let tv = fp2_chip.add_no_carry(ctx, &z, &y);
        let tv = fp2_chip.carry_mod(ctx, tv);
        let tv = fp2_chip.mul(ctx, &tv, &tv);
        let tv = fp2_chip.sub_no_carry(ctx, &tv, &tv1);
        let tv = fp2_chip.sub_no_carry(ctx, &tv, &zsq);
        fp2_chip.carry_mod(ctx, tv)
    };
    let yout = {
        let tv = fp2_chip.sub_no_carry(ctx, &tv3, &xout);
        let tv = fp2_chip.carry_mod(ctx, tv);
        fp2_chip.mul(ctx, tv, &tv4)
    };

    let tv2 = fp2_chip.add_no_carry(ctx, &tv2, &tv2);
    let tv2 = fp2_chip.add_no_carry(ctx, &tv2, &tv2);
    let tv2 = fp2_chip.add_no_carry(ctx, &tv2, &tv2);
    let yout = fp2_chip.sub_no_carry(ctx, &yout, &tv2);
    let yout = fp2_chip.carry_mod(ctx, yout);

    *P = EcProjective { x: xout, y: yout, z: zout };

    if !with_coeffs {
        return vec![];
    }

    let tv3 = fp2_chip.mul(ctx, &tv4, &zsq);
    let tv3 = fp2_chip.add(ctx, &tv3, &tv3);
    let tv3 = fp2_chip.negate(ctx, tv3);
    let tv6 = {
        let tv = fp2_chip.mul(ctx, &tv6, &tv6);
        let tv = fp2_chip.sub_no_carry(ctx, tv, tv0);
        fp2_chip.sub(ctx, tv, &tv5)
    };
    let tv1 = fp2_chip.add_no_carry(ctx, &tv1, &tv1);
    let tv1 = fp2_chip.add_no_carry(ctx, &tv1, &tv1);
    let tv6 = fp2_chip.sub(ctx, tv6, &tv1);
    let tv0 = fp2_chip.mul(ctx, &P.z, zsq);
    let tv0 = fp2_chip.add(ctx, &tv0, &tv0);

    vec![tv0, tv3, tv6]
}

// Frobenius coefficient coeff[1][j] = ((9+u)^{(p-1)/6})^j
// Frob_p( twist(Q) ) = ( (w^2 x)^p, (w^3 y)^p ) = twist( coeff[1][2] * x^p, coeff[1][3] * y^p )
// Input:
// - Q = (x, y) point in E(Fp2)
// - coeff[1][2], coeff[1][3] as assigned cells: this is an optimization to avoid loading new constants
// Output:
// - (coeff[1][2] * x^p, coeff[1][3] * y^p) point in E(Fp2)
pub fn twisted_frobenius<F: PrimeField>(
    ecc_chip: &EccChip<F, Fp2Chip<F>>,
    ctx: &mut Context<F>,
    Q: impl Into<EcPoint<F, FqPoint<F>>>,
    c2: impl Into<FqPoint<F>>,
    c3: impl Into<FqPoint<F>>,
) -> EcPoint<F, FqPoint<F>> {
    let Q = Q.into();
    let c2 = c2.into();
    let c3 = c3.into();
    assert_eq!(c2.0.len(), 2);
    assert_eq!(c3.0.len(), 2);

    let frob_x = ecc_chip.field_chip.conjugate(ctx, Q.x);
    let frob_y = ecc_chip.field_chip.conjugate(ctx, Q.y);
    let out_x = ecc_chip.field_chip.mul(ctx, c2, frob_x);
    let out_y = ecc_chip.field_chip.mul(ctx, c3, frob_y);
    EcPoint::new(out_x, out_y)
}

// Frobenius coefficient coeff[1][j] = ((9+u)^{(p-1)/6})^j
// -Frob_p( twist(Q) ) = ( (w^2 x)^p, -(w^3 y)^p ) = twist( coeff[1][2] * x^p, coeff[1][3] * -y^p )
// Input:
// - Q = (x, y) point in E(Fp2)
// Output:
// - (coeff[1][2] * x^p, coeff[1][3] * -y^p) point in E(Fp2)
pub fn neg_twisted_frobenius<F: PrimeField>(
    ecc_chip: &EccChip<F, Fp2Chip<F>>,
    ctx: &mut Context<F>,
    Q: impl Into<EcPoint<F, FqPoint<F>>>,
    c2: impl Into<FqPoint<F>>,
    c3: impl Into<FqPoint<F>>,
) -> EcPoint<F, FqPoint<F>> {
    let Q = Q.into();
    let c2 = c2.into();
    let c3 = c3.into();
    assert_eq!(c2.0.len(), 2);
    assert_eq!(c3.0.len(), 2);

    let frob_x = ecc_chip.field_chip.conjugate(ctx, Q.x);
    let neg_frob_y = ecc_chip.field_chip.neg_conjugate(ctx, Q.y);
    let out_x = ecc_chip.field_chip.mul(ctx, c2, frob_x);
    let out_y = ecc_chip.field_chip.mul(ctx, c3, neg_frob_y);
    EcPoint::new(out_x, out_y)
}

// To avoid issues with mutably borrowing twice (not allowed in Rust), we only store fp_chip and construct g2_chip and fp12_chip in scope when needed for temporary mutable borrows
pub struct PairingChip<'chip, F: PrimeField> {
    pub fp_chip: &'chip FpChip<'chip, F>,
}

impl<'chip, F: PrimeField> PairingChip<'chip, F> {
    pub fn new(fp_chip: &'chip FpChip<F>) -> Self {
        Self { fp_chip }
    }

    pub fn load_private_g1_unchecked(
        &self,
        ctx: &mut Context<F>,
        point: G1Affine,
    ) -> EcPoint<F, FpPoint<F>> {
        let g1_chip = EccChip::new(self.fp_chip);
        g1_chip.load_private_unchecked(ctx, (point.x, point.y))
    }

    pub fn load_private_g2_unchecked(
        &self,
        ctx: &mut Context<F>,
        point: G2Affine,
    ) -> EcPoint<F, FqPoint<F>> {
        let fp2_chip = Fp2Chip::new(self.fp_chip);
        let g2_chip = EccChip::new(&fp2_chip);
        g2_chip.load_private_unchecked(ctx, (point.x, point.y))
    }

    pub fn miller_loop(
        &self,
        ctx: &mut Context<F>,
        Q: &EcPoint<F, FqPoint<F>>,
        P: &EcPoint<F, FpPoint<F>>,
    ) -> FqPoint<F> {
        let fp2_chip = Fp2Chip::<F>::new(self.fp_chip);
        let g2_chip = EccChip::new(&fp2_chip);
        // miller_loop_BN::<F>(
        //     &g2_chip,
        //     ctx,
        //     Q,
        //     P,
        //     &SIX_U_PLUS_2_NAF, // pseudo binary encoding for BN254
        // )
        unimplemented!()
    }

    pub fn multi_miller_loop(
        &self,
        ctx: &mut Context<F>,
        pairs: Vec<(&EcPoint<F, FpPoint<F>>, &EcPoint<F, FqPoint<F>>)>,
    ) -> FqPoint<F> {
        let fp2_chip = Fp2Chip::<F>::new(self.fp_chip);
        let g2_chip = EccChip::new(&fp2_chip);
        let f = multi_miller_loop::<F>(&g2_chip, ctx, pairs);
        let fp12_chip = Fp12Chip::<F>::new(self.fp_chip);

        f
    }

    pub fn final_exp(&self, ctx: &mut Context<F>, f: FqPoint<F>) -> FqPoint<F> {
        let fp12_chip = Fp12Chip::<F>::new(self.fp_chip);
        fp12_chip.final_exp(ctx, f)
    }

    // optimal Ate pairing
    pub fn pairing(
        &self,
        ctx: &mut Context<F>,
        Q: &EcPoint<F, FqPoint<F>>,
        P: &EcPoint<F, FpPoint<F>>,
    ) -> FqPoint<F> {
        let f0 = self.miller_loop(ctx, Q, P);
        let fp12_chip = Fp12Chip::<F>::new(self.fp_chip);
        // final_exp implemented in final_exp module
        fp12_chip.final_exp(ctx, f0)
    }

    /*
     * Conducts an efficient pairing check e(P, Q) = e(S, T) using only one
     * final exponentiation. In particular, this constraints
     * (e'(-P, Q)e'(S, T))^x = 1, where e' is the optimal ate pairing without
     * the final exponentiation. Reduces number of necessary advice cells by
     * ~30%.
     */
    pub fn pairing_check(
        &self,
        ctx: &mut Context<F>,
        Q: &EcPoint<F, FqPoint<F>>,
        P: &EcPoint<F, FpPoint<F>>,
        T: &EcPoint<F, FqPoint<F>>,
        S: &EcPoint<F, FpPoint<F>>,
    ) {
        let ecc_chip_fp = EccChip::new(self.fp_chip);
        let negated_P = ecc_chip_fp.negate(ctx, P);
        let mml = self.multi_miller_loop(ctx, vec![(&negated_P, Q), (S, T)]);
        let fp12_chip = Fp12Chip::<F>::new(self.fp_chip);
        let fe = fp12_chip.final_exp(ctx, mml);
        let fp12_one = fp12_chip.load_constant(ctx, Fq12::one());
        fp12_chip.assert_equal(ctx, fe, fp12_one);
    }
}
