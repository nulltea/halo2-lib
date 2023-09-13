use std::marker::PhantomData;

use halo2_base::{utils::modulus, AssignedValue, Context};
use num_bigint::BigUint;

use super::{
    vector::{FieldVector, FieldVectorChip},
    FieldChip, FieldExtConstructor, PrimeField, PrimeFieldChip,
};
use crate::impl_field_ext_chip_common;
use ff::PrimeField as _;

/// Represent Fp12 point as FqPoint with degree = 12
/// `Fp12 = Fp2[w] / (w^6 - u - xi)`
/// This implementation assumes p = 3 (mod 4) in order for the polynomial u^2 + 1 to
/// be irreducible over Fp; i.e., in order for -1 to not be a square (quadratic residue) in Fp
/// This means we store an Fp12 point as `\sum_{i = 0}^6 (a_{i0} + a_{i1} * u) * w^i`
/// This is encoded in an FqPoint of degree 12 as `(a_{00}, ..., a_{50}, a_{01}, ..., a_{51})`
#[derive(Clone, Copy, Debug)]
pub struct Fp12Chip<'a, F: PrimeField, FpChip: FieldChip<F>, Fp12, const XI_0: i64>(
    pub FieldVectorChip<'a, F, FpChip>,
    PhantomData<Fp12>,
);

impl<'a, F, FpChip, Fp12, const XI_0: i64> Fp12Chip<'a, F, FpChip, Fp12, XI_0>
where
    F: PrimeField,
    FpChip: PrimeFieldChip<F>,
    FpChip::FieldType: PrimeField,
    Fp12: ff::Field,
{
    /// User must construct an `FpChip` first using a config. This is intended so everything shares a single `FlexGateChip`, which is needed for the column allocation to work.
    pub fn new(fp_chip: &'a FpChip) -> Self {
        assert_eq!(
            modulus::<FpChip::FieldType>() % 4usize,
            BigUint::from(3u64),
            "p must be 3 (mod 4) for the polynomial u^2 + 1 to be irreducible"
        );
        Self(FieldVectorChip::new(fp_chip), PhantomData)
    }

    pub fn fp_chip(&self) -> &FpChip {
        self.0.fp_chip
    }

    pub fn fp2_mul_no_carry(
        &self,
        ctx: &mut Context<F>,
        fp12_pt: FieldVector<FpChip::UnsafeFieldPoint>,
        fp2_pt: FieldVector<FpChip::UnsafeFieldPoint>,
    ) -> FieldVector<FpChip::UnsafeFieldPoint> {
        let fp12_pt = fp12_pt.0;
        let fp2_pt = fp2_pt.0;
        assert_eq!(fp12_pt.len(), 12);
        assert_eq!(fp2_pt.len(), 2);

        let fp_chip = self.fp_chip();
        let mut out_coeffs = Vec::with_capacity(12);
        for i in 0..6 {
            let coeff1 = fp_chip.mul_no_carry(ctx, fp12_pt[i].clone(), fp2_pt[0].clone());
            let coeff2 = fp_chip.mul_no_carry(ctx, fp12_pt[i + 6].clone(), fp2_pt[1].clone());
            let coeff = fp_chip.sub_no_carry(ctx, coeff1, coeff2);
            out_coeffs.push(coeff);
        }
        for i in 0..6 {
            let coeff1 = fp_chip.mul_no_carry(ctx, fp12_pt[i + 6].clone(), fp2_pt[0].clone());
            let coeff2 = fp_chip.mul_no_carry(ctx, fp12_pt[i].clone(), fp2_pt[1].clone());
            let coeff = fp_chip.add_no_carry(ctx, coeff1, coeff2);
            out_coeffs.push(coeff);
        }
        FieldVector(out_coeffs)
    }

    // for \sum_i (a_i + b_i u) w^i, returns \sum_i (-1)^i (a_i + b_i u) w^i
    pub fn conjugate(
        &self,
        ctx: &mut Context<F>,
        a: FieldVector<FpChip::FieldPoint>,
    ) -> FieldVector<FpChip::FieldPoint> {
        let a = a.0;
        assert_eq!(a.len(), 12);

        let coeffs = a
            .into_iter()
            .enumerate()
            .map(|(i, c)| if i % 2 == 0 { c } else { self.fp_chip().negate(ctx, c) })
            .collect();
        FieldVector(coeffs)
    }
}

/// multiply Fp2 elts: (a0 + a1 * u) * (XI0 + u) without carry
///
/// # Assumptions
/// * `a` is `Fp2` point represented as `FieldVector` with degree = 2
pub fn mul_no_carry_w6<F: PrimeField, FC: FieldChip<F>, const XI_0: i64>(
    fp_chip: &FC,
    ctx: &mut Context<F>,
    a: FieldVector<FC::UnsafeFieldPoint>,
) -> FieldVector<FC::UnsafeFieldPoint> {
    let [a0, a1]: [_; 2] = a.0.try_into().unwrap();
    // (a0 + a1 u) * (XI_0 + u) = (a0 * XI_0 - a1) + (a1 * XI_0 + a0) u     with u^2 = -1
    // This should fit in the overflow representation if limb_bits is large enough
    let a0_xi0 = fp_chip.scalar_mul_no_carry(ctx, a0.clone(), XI_0);
    let out0_0_nocarry = fp_chip.sub_no_carry(ctx, a0_xi0, a1.clone());
    let out0_1_nocarry = fp_chip.scalar_mul_and_add_no_carry(ctx, a1, a0, XI_0);
    FieldVector(vec![out0_0_nocarry, out0_1_nocarry])
}

// a lot of this is common to any field extension (lots of for loops), but due to the way rust traits work, it is hard to create a common generic trait that does this. The main problem is that if you had a `FieldExtCommon` trait and wanted to implement `FieldChip` for anything with `FieldExtCommon`, rust will stop you because someone could implement `FieldExtCommon` and `FieldChip` for the same type, causing a conflict.
// partially solved using macro

impl<'a, F, FpChip, Fp12, const XI_0: i64> FieldChip<F> for Fp12Chip<'a, F, FpChip, Fp12, XI_0>
where
    F: PrimeField,
    FpChip: PrimeFieldChip<F>,
    FpChip::FieldType: PrimeField,
    Fp12: ff::Field + FieldExtConstructor<FpChip::FieldType, 12>,
    FieldVector<FpChip::UnsafeFieldPoint>: From<FieldVector<FpChip::FieldPoint>>,
    FieldVector<FpChip::FieldPoint>: From<FieldVector<FpChip::ReducedFieldPoint>>,
{
    const PRIME_FIELD_NUM_BITS: u32 = FpChip::FieldType::NUM_BITS;
    type UnsafeFieldPoint = FieldVector<FpChip::UnsafeFieldPoint>;
    type FieldPoint = FieldVector<FpChip::FieldPoint>;
    type ReducedFieldPoint = FieldVector<FpChip::ReducedFieldPoint>;
    type FieldType = Fp12;
    type RangeChip = FpChip::RangeChip;

    fn get_assigned_value(&self, x: &Self::UnsafeFieldPoint) -> Fp12 {
        assert_eq!(x.0.len(), 12);
        let values = x.0.iter().map(|v| self.fp_chip().get_assigned_value(v)).collect::<Vec<_>>();
        Fp12::new(values.try_into().unwrap())
    }

    // w^6 = u + xi for xi = 9
    fn mul_no_carry(
        &self,
        ctx: &mut Context<F>,
        a: impl Into<Self::UnsafeFieldPoint>,
        b: impl Into<Self::UnsafeFieldPoint>,
    ) -> Self::UnsafeFieldPoint {
        let a = a.into().0;
        let b = b.into().0;
        assert_eq!(a.len(), 12);
        assert_eq!(b.len(), 12);

        let fp_chip = self.fp_chip();
        // a = \sum_{i = 0}^5 (a_i * w^i + a_{i + 6} * w^i * u)
        // b = \sum_{i = 0}^5 (b_i * w^i + b_{i + 6} * w^i * u)
        let mut a0b0_coeffs: Vec<FpChip::UnsafeFieldPoint> = Vec::with_capacity(11);
        let mut a0b1_coeffs: Vec<FpChip::UnsafeFieldPoint> = Vec::with_capacity(11);
        let mut a1b0_coeffs: Vec<FpChip::UnsafeFieldPoint> = Vec::with_capacity(11);
        let mut a1b1_coeffs: Vec<FpChip::UnsafeFieldPoint> = Vec::with_capacity(11);
        for i in 0..6 {
            for j in 0..6 {
                let coeff00 = fp_chip.mul_no_carry(ctx, &a[i], &b[j]);
                let coeff01 = fp_chip.mul_no_carry(ctx, &a[i], &b[j + 6]);
                let coeff10 = fp_chip.mul_no_carry(ctx, &a[i + 6], &b[j]);
                let coeff11 = fp_chip.mul_no_carry(ctx, &a[i + 6], &b[j + 6]);
                if i + j < a0b0_coeffs.len() {
                    a0b0_coeffs[i + j] = fp_chip.add_no_carry(ctx, &a0b0_coeffs[i + j], coeff00);
                    a0b1_coeffs[i + j] = fp_chip.add_no_carry(ctx, &a0b1_coeffs[i + j], coeff01);
                    a1b0_coeffs[i + j] = fp_chip.add_no_carry(ctx, &a1b0_coeffs[i + j], coeff10);
                    a1b1_coeffs[i + j] = fp_chip.add_no_carry(ctx, &a1b1_coeffs[i + j], coeff11);
                } else {
                    a0b0_coeffs.push(coeff00);
                    a0b1_coeffs.push(coeff01);
                    a1b0_coeffs.push(coeff10);
                    a1b1_coeffs.push(coeff11);
                }
            }
        }

        let mut a0b0_minus_a1b1 = Vec::with_capacity(11);
        let mut a0b1_plus_a1b0 = Vec::with_capacity(11);
        for i in 0..11 {
            let a0b0_minus_a1b1_entry = fp_chip.sub_no_carry(ctx, &a0b0_coeffs[i], &a1b1_coeffs[i]);
            let a0b1_plus_a1b0_entry = fp_chip.add_no_carry(ctx, &a0b1_coeffs[i], &a1b0_coeffs[i]);

            a0b0_minus_a1b1.push(a0b0_minus_a1b1_entry);
            a0b1_plus_a1b0.push(a0b1_plus_a1b0_entry);
        }

        // out_i       = a0b0_minus_a1b1_i + XI_0 * a0b0_minus_a1b1_{i + 6} - a0b1_plus_a1b0_{i + 6}
        // out_{i + 6} = a0b1_plus_a1b0_{i} + a0b0_minus_a1b1_{i + 6} + XI_0 * a0b1_plus_a1b0_{i + 6}
        let mut out_coeffs = Vec::with_capacity(12);
        for i in 0..6 {
            if i < 5 {
                let mut coeff = fp_chip.scalar_mul_and_add_no_carry(
                    ctx,
                    &a0b0_minus_a1b1[i + 6],
                    &a0b0_minus_a1b1[i],
                    XI_0,
                );
                coeff = fp_chip.sub_no_carry(ctx, coeff, &a0b1_plus_a1b0[i + 6]);
                out_coeffs.push(coeff);
            } else {
                out_coeffs.push(a0b0_minus_a1b1[i].clone());
            }
        }
        for i in 0..6 {
            if i < 5 {
                let mut coeff =
                    fp_chip.add_no_carry(ctx, &a0b1_plus_a1b0[i], &a0b0_minus_a1b1[i + 6]);
                coeff =
                    fp_chip.scalar_mul_and_add_no_carry(ctx, &a0b1_plus_a1b0[i + 6], coeff, XI_0);
                out_coeffs.push(coeff);
            } else {
                out_coeffs.push(a0b1_plus_a1b0[i].clone());
            }
        }
        FieldVector(out_coeffs)
    }

    impl_field_ext_chip_common!();
}

mod bn254 {
    use crate::fields::FieldExtConstructor;
    use crate::halo2_proofs::halo2curves::bn256::{Fq, Fq12, Fq2, Fq6};
    // This means we store an Fp12 point as `\sum_{i = 0}^6 (a_{i0} + a_{i1} * u) * w^i`
    // This is encoded in an FqPoint of degree 12 as `(a_{00}, ..., a_{50}, a_{01}, ..., a_{51})`
    impl FieldExtConstructor<Fq, 12> for Fq12 {
        fn new(c: [Fq; 12]) -> Self {
            Fq12 {
                c0: Fq6 {
                    c0: Fq2 { c0: c[0], c1: c[6] },
                    c1: Fq2 { c0: c[2], c1: c[8] },
                    c2: Fq2 { c0: c[4], c1: c[10] },
                },
                c1: Fq6 {
                    c0: Fq2 { c0: c[1], c1: c[7] },
                    c1: Fq2 { c0: c[3], c1: c[9] },
                    c2: Fq2 { c0: c[5], c1: c[11] },
                },
            }
        }

        fn coeffs(&self) -> Vec<Fq> {
            let x = self;
            vec![
                x.c0.c0.c0, x.c1.c0.c0, x.c0.c1.c0, x.c1.c1.c0, x.c0.c2.c0, x.c1.c2.c0, x.c0.c0.c1,
                x.c1.c0.c1, x.c0.c1.c1, x.c1.c1.c1, x.c0.c2.c1, x.c1.c2.c1,
            ]
        }
    }
}

mod bls12_381 {
    use crate::{
        bls12_381::{
            pairing::{fq2_mul_by_nonresidue, fq6_mul, fq6_mul_by_nonresidue},
            Fp12Chip, Fp2Chip, FpChip, FqPoint,
        },
        fields::{
            vector::{FieldVector, FieldVectorChip},
            FieldChip, FieldExtConstructor,
        },
    };
    use halo2_base::{
        halo2_proofs::halo2curves::bn256::Fr, safe_types::RangeChip, utils::ScalarField, Context,
    };
    use halo2curves::bls12_381::{Fq, Fq12, Fq2, Fq6};
    use itertools::Itertools;
    // This means we store an Fp12 point as `\sum_{i = 0}^6 (a_{i0} + a_{i1} * u) * w^i`
    // This is encoded in an FqPoint of degree 12 as `(a_{00}, ..., a_{50}, a_{01}, ..., a_{51})`
    impl FieldExtConstructor<Fq, 12> for Fq12 {
        fn new(c: [Fq; 12]) -> Self {
            Fq12 {
                c0: Fq6 {
                    c0: Fq2 { c0: c[0], c1: c[6] },
                    c1: Fq2 { c0: c[2], c1: c[8] },
                    c2: Fq2 { c0: c[4], c1: c[10] },
                },
                c1: Fq6 {
                    c0: Fq2 { c0: c[1], c1: c[7] },
                    c1: Fq2 { c0: c[3], c1: c[9] },
                    c2: Fq2 { c0: c[5], c1: c[11] },
                },
            }
        }

        fn coeffs(&self) -> Vec<Fq> {
            let x = self;
            vec![
                x.c0.c0.c0, x.c1.c0.c0, x.c0.c1.c0, x.c1.c1.c0, x.c0.c2.c0, x.c1.c2.c0, x.c0.c0.c1,
                x.c1.c0.c1, x.c0.c1.c1, x.c1.c1.c1, x.c0.c2.c1, x.c1.c2.c1,
            ]
        }
    }

    #[test]
    fn fq12_mul_test() {
        let a = Fq12::new([
            "095f4afdc084d57eee56b4202c532a0df4c3ddd95462597b984a365b728e9c81fb2bb317f5dfb9711e30875e2165aff3",
            "1435572b6fb6067632a4182a3bdd595dd4621ccc11cf4000111ed2c1b6b66422e911c96e6d069e542d258785d021ac00",
            "0668c42fd3ba46102620bab4712ab825fa3922cc2ff5d8ef1f6f8063caf12a97f0bd9926063ba5b4816b7234aa2febcd",
            "0de0baf382f95f2b82fdae92e937e69b2fb1edbfbfd1a6e4de4ace750418307725f43ca854c898b763aee2656da078d5",
            "0be98d82c63cf9984bb08fdbfb10171a12f4fc97d9daf2a592cce173802c73ea934fd3d3f3c1697b8a54978b0312ef90",
            "173f49f514603120824d50acb1bc8d379086bcbeeb9d3b044cfa52b244bb0bf7516c01277ec34fb804957a4402b7feb7",
            "036a5c0c6bd528e41bfbf01fbb1ae0b1f07c68e891b1177e1cffa7873fb0744d01b8948e58158be9496fcbcecc9c5c41",
            "0fa00a6dd146ee3f380fc1da81f0328f8d4db5be7edc08c71031896621a9f547acc8885f60bed73f80db2c0df30d670d",
            "088cbdf6007360acd11df5a86c4af64b1b04e644efd10f6224354558edf664d209c042f299e4c4ac495f947f68045dbd",
            "01fc612df463dca68ee17dfb03a3d2930a03ff1d5ebe7e89ab018805a418fc8fe82f289d09f37ec468da09bd8a69fd7d",
            "0348b0220c9154b68c5ee57e1c8da961c7666ed36b043cd0192677778fb676e92c96a20be40b7118e1bbf6edc4a64b12",
            "180419099d423d6a7a0c6f7562e1f2c0926f4ae08a5bb1b3deb5e23e19884e39f06870186a9e4bc7354f45abfcbd4e37",
        ].map(|e| Fq::from_bytes(&hex::decode(e).unwrap().into_iter().rev().collect_vec().try_into().unwrap()).unwrap()));

        let b = Fq12::new([
            "0e067ba6baee481ab2cf183e940f4a0d0a0a227b689a0a1a22fe842bbd6968578f39e68e10c9da589e187e0c0dac4724",
            "00454f0072513adce8c7b3a8265d549a0f3e23679dbea1ca2e57b9fd547eda1c3891933d2f8af06240ebad3d13fa896e",
            "164676476a575b26801c395b81d24e911c05997dd94871c36e9b801b6a588319af8ad536a610ae30a7d503bccb328783",
            "1327595233db66e264d694b398e073403aab1359eb24dfd3940c46ed5373bef283e2454750c718eda16c0a28e3e238c1",
            "18870dbe1314a1c093a5db4762daf457c74be3740eeff23a8cba16eb221c95157e858443f07558d2cc17ade048e775e2",
            "106a97367211091a444381f22b3294c0a1c02a097c2cf32e7dbbf66ee96a5cabd96f9cec2a02e3bfa598218726205ebf",
            "05a4442e5c49306100db33bcdde02d82023ad218f0c60ab1a50c3e2549ac79faf8d5dbc45698249a40232c40d5c15b49",
            "0622b4a66469f188a595dc048f6a292db7bbe1644dbdd8f2f3cd358ccdab63e201bb6046f09236a8f38a8fd461964a24",
            "0462f7612d700d43404d00bfe6d5047bf4f0584c49d2f7aa39d0fabad95803629675050261f5965afc29c7995c74f3dd",
            "04ae33b025312bcfa99ae2e0006fc98bf5d2c6e693bb2ef6e36250f6a2c2b4ed70e5561d307a0e56ee4b5beb279bf3eb",
            "0a79a216dcb1fdf6e68b59396a38f4e6ce2074c8ef9ccec7a0b9c7e0d5ec1e87d78a7f84e751ab83cb8badb85687d713",
            "160412c671989f8dcdc99bda8ff56e0e103207de43049e53dd5382d6d406bd0b63e75ad49a41f2e5289dcbb970b35ee6",
        ].map(|e| Fq::from_bytes(&hex::decode(e).unwrap().into_iter().rev().collect_vec().try_into().unwrap()).unwrap()));

        let mut ctx = Context::<Fr>::new(true, 0);
        let range_chip = RangeChip::<Fr>::new(halo2_base::gates::range::RangeStrategy::Vertical, 8);
        let fp_chip = FpChip::<Fr>::new(&range_chip, 112, 4);
        let fp12_chip = super::Fp12Chip::<Fr, FpChip<Fr>, Fq12, 1>::new(&fp_chip);

        let a_assigned = fp12_chip.load_private(&mut ctx, a);
        let b_assigned = fp12_chip.load_private(&mut ctx, b);

        // let c = fq12_mul(&fp12_chip, &mut ctx, &a_assigned, &b_assigned);
        let c = fp12_chip.mul(&mut ctx, &a_assigned, &b_assigned);
        let c_value = fp12_chip.get_assigned_value(&c.into());

        let c_check = a * b;

        assert_eq!(c_value, c_check);
    }
}
