use std::cell::RefCell;

use crate::{
    sha256::vanilla::{
        columns::Sha256CircuitConfig, param::*, util::get_sha2_capacity,
        witness::AssignedSha256Block,
    },
    util::{eth_types::Field, word::Word},
};
use getset::{CopyGetters, Getters};
use halo2_base::{
    gates::{
        circuit::{builder::BaseCircuitBuilder, BaseCircuitParams, BaseConfig},
        flex_gate::MultiPhaseThreadBreakPoints,
    },
    halo2_proofs::{
        circuit::{Layouter, SimpleFloorPlanner},
        plonk::{Circuit, ConstraintSystem, Error},
    },
    AssignedValue,
};
use snark_verifier_sdk::CircuitExt;

/// Sha Coprocessor Leaf Circuit
#[derive(Getters)]
pub struct ShaCoprocessorLeafCircuit<F: Field> {
    inputs: Vec<Vec<u8>>,
    num_rows: Option<usize>,

    /// Parameters of this circuit. The same parameters always construct the same circuit.
    #[getset(get = "pub")]
    params: ShaCoprocessorLeafCircuitParams,

    base_circuit_builder: RefCell<BaseCircuitBuilder<F>>,
    // hasher: RefCell<PoseidonHasher<F, POSEIDON_T, POSEIDON_RATE>>,
}

/// Parameters of KeccakCoprocessorLeafCircuit.
#[derive(Default, Clone, CopyGetters)]
pub struct ShaCoprocessorLeafCircuitParams {
    /// This circuit has 2^k rows.
    #[getset(get_copy = "pub")]
    k: usize,
    // Number of unusable rows withhold by Halo2.
    #[getset(get_copy = "pub")]
    num_unusable_row: usize,
    /// The bits of lookup table for RangeChip.
    #[getset(get_copy = "pub")]
    lookup_bits: usize,
    /// Max keccak_f this circuits can aceept. The circuit can at most process <capacity> of inputs
    /// with < NUM_BYTES_TO_ABSORB bytes or an input with <capacity> * NUM_BYTES_TO_ABSORB - 1 bytes.
    #[getset(get_copy = "pub")]
    capacity: usize,
    // If true, publish raw outputs. Otherwise, publish Poseidon commitment of raw outputs.
    #[getset(get_copy = "pub")]
    publish_raw_outputs: bool,

    // Derived parameters of sub-circuits.
    // pub sha_circuit_params: KeccakConfigParams,
    pub base_circuit_params: BaseCircuitParams,
}

impl ShaCoprocessorLeafCircuitParams {
    /// Create a new KeccakCoprocessorLeafCircuitParams.
    pub fn new(
        k: usize,
        num_unusable_row: usize,
        lookup_bits: usize,
        capacity: usize,
        publish_raw_outputs: bool,
    ) -> Self {
        assert!(1 << k > num_unusable_row, "Number of unusable rows must be less than 2^k");
        let base_circuit_params = BaseCircuitParams {
            k,
            lookup_bits: Some(lookup_bits),
            num_instance_columns: 0,
            // if publish_raw_outputs {
            //     OUTPUT_NUM_COL_RAW
            // } else {
            //     OUTPUT_NUM_COL_COMMIT
            // },
            ..Default::default()
        };
        Self {
            k,
            num_unusable_row,
            lookup_bits,
            capacity,
            publish_raw_outputs,
            base_circuit_params,
        }
    }
}

/// Circuit::Config for Keccak Coprocessor Leaf Circuit.
#[derive(Clone)]
pub struct ShaCoprocessorLeafConfig<F: Field> {
    pub base_circuit_config: BaseConfig<F>,
    pub sha_circuit_config: Sha256CircuitConfig<F>,
}

impl<F: Field> Circuit<F> for ShaCoprocessorLeafCircuit<F> {
    type Config = ShaCoprocessorLeafConfig<F>;
    type FloorPlanner = SimpleFloorPlanner;
    type Params = ShaCoprocessorLeafCircuitParams;

    fn params(&self) -> Self::Params {
        self.params.clone()
    }

    /// Creates a new instance of the [RangeCircuitBuilder] without witnesses by setting the witness_gen_only flag to false
    fn without_witnesses(&self) -> Self {
        unimplemented!()
    }

    /// Configures a new circuit using [`BaseConfigParams`]
    fn configure_with_params(meta: &mut ConstraintSystem<F>, params: Self::Params) -> Self::Config {
        let sha_circuit_config = Sha256CircuitConfig::new(meta);
        let base_circuit_params = params.base_circuit_params;
        // BaseCircuitBuilder::configure_with_params must be called in the end in order to get the correct
        // unusable_rows.
        let base_circuit_config =
            BaseCircuitBuilder::configure_with_params(meta, base_circuit_params.clone());
        Self::Config { base_circuit_config, sha_circuit_config }
    }

    fn configure(_: &mut ConstraintSystem<F>) -> Self::Config {
        unreachable!("You must use configure_with_params");
    }

    fn synthesize(
        &self,
        config: Self::Config,
        mut layouter: impl Layouter<F>,
    ) -> Result<(), Error> {
        let k = self.params.k;
        let sha_assigned_blocks = layouter.assign_region(
            || "sha256 bit circuit",
            |mut region| {
                Ok(config.sha_circuit_config.multi_sha256(
                    &mut region,
                    self.inputs.clone(),
                    self.num_rows.map(get_sha2_capacity),
                ))
            },
        )?;

        // Base circuit witness generation.
        // let loaded_keccak_fs = self.load_sha_assigned_blocks(sha_assigned_blocks);
        // self.generate_base_circuit_witnesses(&loaded_keccak_fs);

        // self.base_circuit_builder.borrow().synthesize(config.base_circuit_config, layouter)?;

        // // Reset the circuit to the initial state so synthesize could be called multiple times.
        // self.base_circuit_builder.borrow_mut().clear();
        // self.hasher.borrow_mut().clear();
        Ok(())
    }
}

impl<F: Field> CircuitExt<F> for ShaCoprocessorLeafCircuit<F> {
    fn num_instance(&self) -> Vec<usize> {
        vec![]
    }

    fn instances(&self) -> Vec<Vec<F>> {
        vec![]
    }
}

/// Witnesses of a sha256 which are necessary to be loaded into halo2-lib.
#[derive(Clone, Copy, Debug, CopyGetters, Getters)]
pub struct LoadedSha256<F: Field> {
    /// bytes_left of the first row of the first round of this keccak_f. This could be used to determine the length of the input.
    // #[getset(get_copy = "pub")]
    // pub(crate) bytes_left: AssignedValue<F>,

    /// The output of this sha256. is_final/hash_lo/hash_hi come from the first row of the last round(NUM_ROUNDS).
    #[getset(get_copy = "pub")]
    pub is_final: AssignedValue<F>,

    // Hash word consisting of two limbs - lower 16 bits and the high 16 bits, in big-endian.
    #[getset(get = "pub")]
    pub hash: Word<AssignedValue<F>>,

    /// Input words (u64) of this keccak_f.
    #[getset(get = "pub")]
    pub word_values: [AssignedValue<F>; NUM_WORDS_TO_ABSORB],
    // pub(crate) length: AssignedValue<F>,
}

impl<F: Field> ShaCoprocessorLeafCircuit<F> {
    /// Create a new KeccakCoprocessorLeafCircuit.
    pub fn new(
        inputs: Vec<Vec<u8>>,
        num_rows: Option<usize>,
        params: ShaCoprocessorLeafCircuitParams,
        witness_gen_only: bool,
    ) -> Self {
        // let input_size = inputs.iter().map(|input| get_num_keccak_f(input.len())).sum::<usize>();
        // assert!(input_size < params.capacity, "Input size exceeds capacity");
        let mut base_circuit_builder = BaseCircuitBuilder::new(witness_gen_only);
        base_circuit_builder.set_params(params.base_circuit_params.clone());
        Self {
            inputs,
            num_rows,
            params,
            base_circuit_builder: RefCell::new(base_circuit_builder),
            // hasher: RefCell::new(create_hasher()),
        }
    }

    /// Get break points of BaseCircuitBuilder.
    pub fn base_circuit_break_points(&self) -> MultiPhaseThreadBreakPoints {
        self.base_circuit_builder.borrow().break_points()
    }

    /// Set break points of BaseCircuitBuilder.
    pub fn set_base_circuit_break_points(&self, break_points: MultiPhaseThreadBreakPoints) {
        self.base_circuit_builder.borrow_mut().set_break_points(break_points);
    }

    pub fn update_base_circuit_params(&mut self, params: &BaseCircuitParams) {
        self.params.base_circuit_params = params.clone();
        self.base_circuit_builder.borrow_mut().set_params(params.clone());
    }

    // /// Simulate witness generation of the base circuit to determine BaseCircuitParams because the number of columns
    // /// of the base circuit can only be known after witness generation.
    // pub fn calculate_base_circuit_params(
    //     params: &ShaCoprocessorLeafCircuitParams,
    // ) -> BaseCircuitParams {
    //     // Create a simulation circuit to calculate base circuit parameters.
    //     let simulation_circuit = Self::new(vec![], params.clone(), false);
    //     let loaded_keccak_fs = simulation_circuit.mock_load_keccak_assigned_rows();
    //     simulation_circuit.generate_base_circuit_witnesses(&loaded_keccak_fs);

    //     let base_circuit_params = simulation_circuit
    //         .base_circuit_builder
    //         .borrow_mut()
    //         .calculate_params(Some(params.num_unusable_row));
    //     // prevent drop warnings
    //     simulation_circuit.base_circuit_builder.borrow_mut().clear();

    //     base_circuit_params
    // }

    // /// Mock loading Keccak assigned rows from Keccak circuit. This function doesn't create any witnesses/constraints.
    // fn mock_load_keccak_assigned_rows(&self) -> Vec<LoadedShaF<F>> {
    //     let base_circuit_builder = self.base_circuit_builder.borrow();
    //     let mut copy_manager = base_circuit_builder.core().copy_manager.lock().unwrap();
    //     (0..self.params.capacity)
    //         .map(|_| LoadedShaF {
    //             bytes_left: copy_manager.mock_external_assigned(F::ZERO),
    //             word_values: core::array::from_fn(|_| copy_manager.mock_external_assigned(F::ZERO)),
    //             is_final: SafeTypeChip::unsafe_to_bool(
    //                 copy_manager.mock_external_assigned(F::ZERO),
    //             ),
    //             hash_lo: copy_manager.mock_external_assigned(F::ZERO),
    //             hash_hi: copy_manager.mock_external_assigned(F::ZERO),
    //         })
    //         .collect_vec()
    // }

    // /// Load needed witnesses into halo2-lib from sha assigned rows. This function doesn't create any witnesses/constraints.
    fn load_sha_assigned_blocks(
        &self,
        assigned_blocks: Vec<AssignedSha256Block<'_, F>>,
    ) -> Vec<LoadedSha256<F>> {
        // let rows_per_round = self.params.sha_circuit_params.rows_per_round;
        let base_circuit_builder = self.base_circuit_builder.borrow();
        let mut copy_manager = base_circuit_builder.core().copy_manager.lock().unwrap();
        assigned_blocks
            .into_iter()
            .map(|block| {
                // let mut rounds = block.collect_vec();
                // assert_eq!(rounds.len(), NUM_ROUNDS + 1);
                // let bytes_left = copy_manager.load_external_assigned(rounds[0].bytes_left.clone());
                // let output_row = rounds.pop().unwrap();
                let word_values = core::array::from_fn(|i| {
                    let cell = &block.word_values[i];
                    copy_manager.load_external_assigned(cell.clone())
                });
                let is_final = copy_manager.load_external_assigned(block.is_final);
                let hash_lo = copy_manager.load_external_assigned(block.output.lo());
                let hash_hi = copy_manager.load_external_assigned(block.output.hi());
                LoadedSha256 { word_values, is_final, hash: Word::new([hash_lo, hash_hi]) }
            })
            .collect()
    }

    // /// Generate witnesses of the base circuit.
    // fn generate_base_circuit_witnesses(&self, loaded_keccak_fs: &[LoadedShaF<F>]) {
    //     let range = self.base_circuit_builder.borrow().range_chip();
    //     let gate = range.gate();
    //     let circuit_final_outputs = {
    //         let mut base_circuit_builder_mut = self.base_circuit_builder.borrow_mut();
    //         let ctx = base_circuit_builder_mut.main(0);
    //         let mut hasher = self.hasher.borrow_mut();
    //         hasher.initialize_consts(ctx, gate);

    //         let lookup_key_per_keccak_f =
    //             encode_inputs_from_keccak_fs(ctx, gate, &hasher, loaded_keccak_fs);
    //         Self::generate_circuit_final_outputs(
    //             ctx,
    //             gate,
    //             &lookup_key_per_keccak_f,
    //             loaded_keccak_fs,
    //         )
    //     };
    //     self.publish_outputs(&circuit_final_outputs);
    // }

    // /// Combine lookup keys and Keccak results to generate final outputs of the circuit.
    // pub fn generate_circuit_final_outputs(
    //     ctx: &mut Context<F>,
    //     gate: &impl GateInstructions<F>,
    //     lookup_key_per_keccak_f: &[PoseidonCompactOutput<F>],
    //     loaded_keccak_fs: &[LoadedShaF<F>],
    // ) -> Vec<KeccakCircuitOutput<AssignedValue<F>>> {
    //     let KeccakCircuitOutput {
    //         key: dummy_key_val,
    //         hash_lo: dummy_keccak_val_lo,
    //         hash_hi: dummy_keccak_val_hi,
    //     } = dummy_circuit_output::<F>();

    //     // Dummy row for keccak_fs with is_final = false. The corresponding logical input is empty.
    //     let dummy_key_witness = ctx.load_constant(dummy_key_val);
    //     let dummy_keccak_lo_witness = ctx.load_constant(dummy_keccak_val_lo);
    //     let dummy_keccak_hi_witness = ctx.load_constant(dummy_keccak_val_hi);

    //     let mut circuit_final_outputs = Vec::with_capacity(loaded_keccak_fs.len());
    //     for (compact_output, loaded_keccak_f) in
    //         lookup_key_per_keccak_f.iter().zip_eq(loaded_keccak_fs)
    //     {
    //         let is_final = AssignedValue::from(loaded_keccak_f.is_final);
    //         let key = gate.select(ctx, *compact_output.hash(), dummy_key_witness, is_final);
    //         let hash_lo =
    //             gate.select(ctx, loaded_keccak_f.hash_lo, dummy_keccak_lo_witness, is_final);
    //         let hash_hi =
    //             gate.select(ctx, loaded_keccak_f.hash_hi, dummy_keccak_hi_witness, is_final);
    //         circuit_final_outputs.push(KeccakCircuitOutput { key, hash_lo, hash_hi });
    //     }
    //     circuit_final_outputs
    // }

    // /// Publish outputs of the circuit as public instances.
    // fn publish_outputs(&self, outputs: &[KeccakCircuitOutput<AssignedValue<F>>]) {
    //     // The length of outputs should always equal to params.capacity.
    //     assert_eq!(outputs.len(), self.params.capacity);
    //     if !self.params.publish_raw_outputs {
    //         let range_chip = self.base_circuit_builder.borrow().range_chip();
    //         let gate = range_chip.gate();
    //         let mut base_circuit_builder_mut = self.base_circuit_builder.borrow_mut();
    //         let ctx = base_circuit_builder_mut.main(0);

    //         // TODO: wrap this into a function which should be shared wiht App circuits.
    //         let output_commitment = self.hasher.borrow().hash_fix_len_array(
    //             ctx,
    //             gate,
    //             &outputs
    //                 .iter()
    //                 .flat_map(|output| [output.key, output.hash_lo, output.hash_hi])
    //                 .collect_vec(),
    //         );

    //         let assigned_instances = &mut base_circuit_builder_mut.assigned_instances;
    //         // The commitment should be in the first row.
    //         assert!(assigned_instances[OUTPUT_COL_IDX_COMMIT].is_empty());
    //         assigned_instances[OUTPUT_COL_IDX_COMMIT].push(output_commitment);
    //     } else {
    //         let assigned_instances = &mut self.base_circuit_builder.borrow_mut().assigned_instances;

    //         // Outputs should be in the top of instance columns.
    //         assert!(assigned_instances[OUTPUT_COL_IDX_KEY].is_empty());
    //         assert!(assigned_instances[OUTPUT_COL_IDX_HASH_LO].is_empty());
    //         assert!(assigned_instances[OUTPUT_COL_IDX_HASH_HI].is_empty());
    //         for output in outputs {
    //             assigned_instances[OUTPUT_COL_IDX_KEY].push(output.key);
    //             assigned_instances[OUTPUT_COL_IDX_HASH_LO].push(output.hash_lo);
    //             assigned_instances[OUTPUT_COL_IDX_HASH_HI].push(output.hash_hi);
    //         }
    //     }
    // }
}

// /// Encode raw inputs from Keccak circuit witnesses into lookup keys.
// ///
// /// Each element in the return value corrresponds to a Keccak chunk. If is_final = true, this element is the lookup key of the corresponding logical input.
// pub fn encode_inputs_from_keccak_fs<F: Field>(
//     ctx: &mut Context<F>,
//     gate: &impl GateInstructions<F>,
//     initialized_hasher: &PoseidonHasher<F, POSEIDON_T, POSEIDON_RATE>,
//     loaded_keccak_fs: &[LoadedShaF<F>],
// ) -> Vec<PoseidonCompactOutput<F>> {
//     // Circuit parameters
//     let num_poseidon_absorb_per_keccak_f = num_poseidon_absorb_per_keccak_f::<F>();
//     let num_word_per_witness = num_word_per_witness::<F>();
//     let num_witness_per_keccak_f = POSEIDON_RATE * num_poseidon_absorb_per_keccak_f;

//     // Constant witnesses
//     let one_const = ctx.load_constant(F::ONE);
//     let zero_const = ctx.load_zero();
//     let multipliers_val = get_words_to_witness_multipliers::<F>()
//         .into_iter()
//         .map(|multiplier| Constant(multiplier))
//         .collect_vec();

//     let mut compact_chunk_inputs = Vec::with_capacity(loaded_keccak_fs.len());
//     let mut last_is_final = one_const;
//     for loaded_keccak_f in loaded_keccak_fs {
//         // If this keccak_f is the last of a logical input.
//         let is_final = loaded_keccak_f.is_final;
//         let mut poseidon_absorb_data = Vec::with_capacity(num_witness_per_keccak_f);

//         // First witness of a keccak_f: [<length_placeholder>, word_values[0], word_values[1], ...]
//         // <length_placeholder> is the length of the input if this is the first keccak_f of a logical input. Otherwise 0.
//         let mut words = Vec::with_capacity(num_word_per_witness);
//         let input_bytes_len = gate.mul(ctx, loaded_keccak_f.bytes_left, last_is_final);
//         words.push(input_bytes_len);
//         words.extend_from_slice(&loaded_keccak_f.word_values);

//         // Turn every num_word_per_witness words later into a witness.
//         for words in words.chunks(num_word_per_witness) {
//             let mut words = words.to_vec();
//             words.resize(num_word_per_witness, zero_const);
//             let witness = gate.inner_product(ctx, words, multipliers_val.clone());
//             poseidon_absorb_data.push(witness);
//         }
//         // Pad 0s to make sure poseidon_absorb_data.len() % RATE == 0.
//         poseidon_absorb_data.resize(num_witness_per_keccak_f, zero_const);
//         let compact_inputs: Vec<_> = poseidon_absorb_data
//             .chunks_exact(POSEIDON_RATE)
//             .map(|chunk| chunk.to_vec().try_into().unwrap())
//             .collect_vec();
//         debug_assert_eq!(compact_inputs.len(), num_poseidon_absorb_per_keccak_f);
//         compact_chunk_inputs.push(PoseidonCompactChunkInput::new(compact_inputs, is_final));
//         last_is_final = is_final.into();
//     }

//     initialized_hasher.hash_compact_chunk_inputs(ctx, gate, &compact_chunk_inputs)
// }

#[cfg(test)]
mod tests {
    use std::path::Path;

    use halo2_base::{
        gates::circuit::CircuitBuilderStage,
        halo2_proofs::plonk::{keygen_pk, keygen_vk},
        utils::{
            fs::gen_srs,
            testing::{check_proof_with_instances, gen_proof_with_instances},
        },
    };
    use itertools::Itertools;

    use super::*;
    use crate::halo2_proofs::{dev::MockProver, halo2curves::bn256::Fr};

    #[test]
    fn test_mock_leaf_circuit_commit() {
        let k: usize = 18;
        let num_unusable_row: usize = 109;
        let lookup_bits: usize = 4;
        let capacity: usize = 1024;
        let publish_raw_outputs: bool = false;

        let inputs = vec![(0u8..64).collect::<Vec<_>>(); 1024];

        let params = ShaCoprocessorLeafCircuitParams::new(
            k,
            num_unusable_row,
            lookup_bits,
            capacity,
            publish_raw_outputs,
        );
        // let base_circuit_params =
        //     ShaCoprocessorLeafCircuit::<Fr>::calculate_base_circuit_params(&params);
        // params.base_circuit_params = base_circuit_params;
        let circuit =
            ShaCoprocessorLeafCircuit::<Fr>::new(inputs.clone(), None, params.clone(), false);
        // let circuit_outputs = multi_inputs_to_circuit_outputs::<Fr>(&inputs, params.capacity());

        // let instances = vec![vec![calculate_circuit_outputs_commit(&circuit_outputs)]];
        let instances = vec![];

        let prover = MockProver::<Fr>::run(k as u32, &circuit, instances).unwrap();
        prover.assert_satisfied();
    }

    #[test]
    fn test_prove_leaf_circuit_commit() {
        let _ = env_logger::builder().is_test(true).try_init();

        let k: usize = 18;
        let num_unusable_row: usize = 109;
        let lookup_bits: usize = 4;
        let capacity: usize = 1024;
        let publish_raw_outputs: bool = false;

        let inputs = vec![];
        let circuit_params = ShaCoprocessorLeafCircuitParams::new(
            k,
            num_unusable_row,
            lookup_bits,
            capacity,
            publish_raw_outputs,
        );
        // let base_circuit_params =
        //     ShaCoprocessorLeafCircuit::<Fr>::calculate_base_circuit_params(&circuit_params);
        // circuit_params.base_circuit_params = base_circuit_params;
        let circuit =
            ShaCoprocessorLeafCircuit::<Fr>::new(inputs, None, circuit_params.clone(), false);

        let params = gen_srs(k as u32);

        let vk = keygen_vk(&params, &circuit).unwrap();
        let pk = keygen_pk(&params, vk, &circuit).unwrap();

        let inputs = vec![(0u8..64).collect::<Vec<_>>(); 1024];

        // let break_points = circuit.base_circuit_break_points();
        let circuit = ShaCoprocessorLeafCircuit::<Fr>::new(
            inputs.clone(),
            None,
            circuit_params.clone(),
            true,
        );
        // circuit.set_base_circuit_break_points(break_points);

        // let circuit_outputs =
        //     multi_inputs_to_circuit_outputs::<Fr>(&inputs, circuit_params.capacity());
        // let instances = vec![vec![calculate_circuit_outputs_commit(&circuit_outputs)]];
        let instances: Vec<Vec<Fr>> = vec![];

        let proof = gen_proof_with_instances(
            &params,
            &pk,
            circuit,
            instances.iter().map(|f| f.as_slice()).collect_vec().as_slice(),
        );
        check_proof_with_instances(
            &params,
            pk.get_vk(),
            &proof,
            instances.iter().map(|f| f.as_slice()).collect_vec().as_slice(),
            true,
        );
    }

    use snark_verifier_sdk::{
        evm::{evm_verify, gen_evm_proof_shplonk, gen_evm_verifier_shplonk},
        halo2::{
            aggregation::{AggregationCircuit, AggregationConfigParams},
            gen_snark_shplonk,
        },
        SHPLONK,
    };

    #[test]
    fn test_aggregate_leaf_circuit_commit() {
        let _ = env_logger::builder().is_test(true).try_init();

        let k: usize = 18;
        let num_unusable_row: usize = 109;
        let lookup_bits: usize = 4;
        let capacity: usize = 1024;
        let publish_raw_outputs: bool = false;

        let inputs = vec![];
        let circuit_params = ShaCoprocessorLeafCircuitParams::new(
            k,
            num_unusable_row,
            lookup_bits,
            capacity,
            publish_raw_outputs,
        );
        // let base_circuit_params =
        //     ShaCoprocessorLeafCircuit::<Fr>::calculate_base_circuit_params(&circuit_params);
        // circuit_params.base_circuit_params = base_circuit_params;
        let circuit =
            ShaCoprocessorLeafCircuit::<Fr>::new(inputs, None, circuit_params.clone(), false);

        let params = gen_srs(k as u32);

        let vk = keygen_vk(&params, &circuit).unwrap();
        let pk = keygen_pk(&params, vk, &circuit).unwrap();

        let inputs = vec![(0u8..64).collect::<Vec<_>>(); 1024];

        // let break_points = circuit.base_circuit_break_points();
        let circuit = ShaCoprocessorLeafCircuit::<Fr>::new(
            inputs.clone(),
            None,
            circuit_params.clone(),
            true,
        );

        // let circuit_outputs =
        //     multi_inputs_to_circuit_outputs::<Fr>(&inputs, circuit_params.capacity());
        // let instances = vec![vec![calculate_circuit_outputs_commit(&circuit_outputs)]];

        println!("gen_snark_shplonk");
        let snark = gen_snark_shplonk(&params, &pk, circuit, None::<String>);

        const K_AGG: u32 = 22;
        let agg_params = gen_srs(K_AGG);

        let circuit_params =
            AggregationConfigParams { degree: K_AGG, lookup_bits: 8, ..Default::default() };

        let mut agg_circuit = AggregationCircuit::new::<SHPLONK>(
            CircuitBuilderStage::Keygen,
            circuit_params,
            &agg_params,
            vec![snark.clone()],
            Default::default(),
        );
        agg_circuit.expose_previous_instances(false);
        let circuit_params = agg_circuit.calculate_params(None);

        let agg_vk = keygen_vk(&agg_params, &agg_circuit).unwrap();
        let agg_pk = keygen_pk(&agg_params, agg_vk.clone(), &agg_circuit).unwrap();

        let break_points = agg_circuit.break_points();

        let mut agg_circuit = AggregationCircuit::new::<SHPLONK>(
            CircuitBuilderStage::Prover,
            circuit_params,
            &params,
            vec![snark.clone()],
            Default::default(),
        );
        agg_circuit.expose_previous_instances(false);
        agg_circuit.set_params(circuit_params);
        agg_circuit.set_break_points(break_points);

        let instances = agg_circuit.instances();
        let num_instances = agg_circuit.num_instance();

        println!("num_instances: {:?}", num_instances);
        println!("instances: {:?}", instances);

        let proof = gen_evm_proof_shplonk(&agg_params, &agg_pk, agg_circuit, instances.clone());
        println!("proof size: {}", proof.len());
        let deployment_code = gen_evm_verifier_shplonk::<AggregationCircuit>(
            &agg_params,
            &agg_vk,
            num_instances,
            Some(Path::new("contractyul")),
        );
        println!("deployment_code size: {}", deployment_code.len());
        evm_verify(deployment_code, instances, proof);
    }
}
