use halo2_base::{utils::ScalarField, Context};

use crate::{
    bigint::ProperCrtUint,
    fields::{fp::FpChip, fp2::Fp2Chip, vector::FieldVector, PrimeField},
};

use super::EcPoint;

pub trait BlsSignatureChip<'chip, F: ScalarField> {
    // Verifies that e(g1, signature) = e(pubkey, H(m)) by checking e(g1, signature)*e(pubkey, -H(m)) === 1
    // where e(,) is optimal Ate pairing
    // G1: {g1, pubkey}, G2: {signature, message}
    fn verify_signature(
        &self,
        ctx: &mut Context<F>,
        signature: EcPoint<F, FieldVector<ProperCrtUint<F>>>,
        msghash: EcPoint<F, FieldVector<ProperCrtUint<F>>>,
        pubkey: EcPoint<F, ProperCrtUint<F>>,
        g1_gen: EcPoint<F, ProperCrtUint<F>>,
        is_strict: bool,
    ) -> FieldVector<ProperCrtUint<F>>;

    fn verify_pairing(
        &self,
        ctx: &mut Context<F>,
        signature: EcPoint<F, FieldVector<ProperCrtUint<F>>>,
        msghash: EcPoint<F, FieldVector<ProperCrtUint<F>>>,
        pubkey: EcPoint<F, ProperCrtUint<F>>,
        g1_neg: EcPoint<F, ProperCrtUint<F>>,
    ) -> FieldVector<ProperCrtUint<F>>;
}
