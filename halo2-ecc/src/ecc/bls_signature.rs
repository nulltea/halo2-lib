use halo2_base::{Context, utils::ScalarField};

use crate::{fields::{fp2::Fp2Chip, fp::FpChip, PrimeField, vector::FieldVector}, bigint::ProperCrtUint};

use super::EcPoint;

pub trait BlsSignatureChip<'chip, F: ScalarField> {
    // fn new(fp_chip: &'chip FpChip<F, Fp>) -> Self;

    // Verifies that e(g1, signature) = e(pubkey, H(m)) by checking e(g1, signature)*e(pubkey, -H(m)) === 1
    // where e(,) is optimal Ate pairing
    // G1: {g1, pubkey}, G2: {signature, message}
   fn verify_signature(
        &self,
        signature: EcPoint<F, FieldVector<ProperCrtUint<F>>>,
        msghash: EcPoint<F, FieldVector<ProperCrtUint<F>>>,
        pubkey: EcPoint<F, ProperCrtUint<F>>,
        g1_gen: EcPoint<F, ProperCrtUint<F>>,
        ctx: &mut Context<F>,
        is_strict: bool,
    ) -> FieldVector<ProperCrtUint<F>>;
}
