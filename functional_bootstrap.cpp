#include "lib.cpp"

enum fbt_operation {
  fbt_op_f_bootstrap, 
  fbt_op_tlwe_keyswitch,
  fbt_op_lwe_keyswitch,
  fbt_op_mv_f_bootstrap,
  fbt_op_mv_f_bootstrap_init,
  fbt_op_trivial,
  fbt_op_load_literal,
  fbt_op_load_hardcoded_LUT,
  fbt_op_generate_mv_b_LUT,
  fbt_op_generate_b_LUT,
  fbt_op_generate_LUT_from_f,
  fbt_op_next_tree_level} fbt_operation;

typedef struct fbt_context {
  const TFheGateBootstrappingParameterSet * params;
  const LweBootstrappingKeyFFT * bootstrap_key;
  const TLweKeySwitchKey * tlweKS;
  LweSample * input, * output;
  TLweSample * tlwe_temp;
  LweSample * lwe_temp_in;
  int log_torus_base, torus_base;
} fbt_context;

typedef void (*fbt_op_trivial_function)(int, int, fbt_context *);
int lut_test1[192] = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,
                      0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
int lut_test2[64] = {1,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2};
int lut_int_comp_base4[4] = {0,1,1,1};
int lut_int_comp_base8[8] = {0,1,1,1,1,1,1};
int lut_add_base4[4] = {0,1,2,3};
int lut_relu_base4[4] = {0,1,0,0};
int * FBT_HARDCODED_LUT_LIST[6] = {lut_test1, lut_test2, lut_int_comp_base4, lut_int_comp_base8, lut_add_base4, lut_relu_base4};
fbt_op_trivial_function * fbt_op_trivial_FUNCTIONS;

fbt_context * fbt_init(const TFheGateBootstrappingParameterSet * params,
                       const LweBootstrappingKeyFFT * bootstrap_key,
                       const TLweKeySwitchKey * tlweKS, int log_torus_base,
                       const int size_of_temps){
  fbt_context * context = (fbt_context *) malloc(sizeof(fbt_context));
  context->params = params;
  context->bootstrap_key = bootstrap_key;
  context->tlweKS = tlweKS;

  context->log_torus_base = log_torus_base;
  context->torus_base = 1 << log_torus_base;

  const TLweParams * tlwe_params = context->params->tgsw_params->tlwe_params;
  const LweParams * extracted_params = &tlwe_params->extracted_lweparams;
  context->input = new_LweSample_array(size_of_temps, extracted_params);
  context->output = new_LweSample_array(size_of_temps, extracted_params);
  context->tlwe_temp = new_TLweSample_array(size_of_temps, tlwe_params);
  context->lwe_temp_in = new_LweSample_array(size_of_temps, params->in_out_params);
  return context;
}

void fbt_main_loop(int ** function, const int size, LweSample * input, fbt_context * context){
  const TLweParams * tlwe_params = context->params->tgsw_params->tlwe_params;
  const LweParams * lwe_params = context->params->in_out_params;
  const LweParams * extracted_params = &tlwe_params->extracted_lweparams;
  int current_input = 0;
  const int log_carpov_base = 1;
  
  TLweSample * TV_0 = new_TLweSample(tlwe_params);
  LweSample * temp;
  IntPolynomial * lut_clear;
  int * lut_input, number_of_luts;

  for (size_t i = 0; i < size; i++)
  {
    switch (function[i][0])
    {
    case fbt_op_f_bootstrap:
      for (size_t j = 0; j < function[i][3]; j++){
        tfhe_functional_bootstrap_woKS_FFT(&context->output[function[i][1] + j], context->bootstrap_key, &context->tlwe_temp[function[i][2] + j], &input[current_input], context->torus_base);
      }
      break;

    case fbt_op_mv_f_bootstrap:
      carpov_second_phase_with_base_composition(&context->output[function[i][1]], lut_clear, TV_0, function[i][3], tlwe_params, context->log_torus_base, log_carpov_base);
      delete_IntPolynomial_array(number_of_luts, lut_clear);
      break;
    
    case fbt_op_mv_f_bootstrap_init:
      tfhe_carpov_multivalue_bootstrap_init(TV_0, context->bootstrap_key, &input[current_input], context->torus_base);
      break;

    case fbt_op_generate_mv_b_LUT:
      IntPolynomial * luts;
      number_of_luts = generate_polynomial_LUT(&luts, lut_input, function[i][3], context->torus_base, tlwe_params->N);
      lut_clear = new_IntPolynomial_array((context->log_torus_base/log_carpov_base)*(function[i][3]/context->torus_base), tlwe_params->N);
      for (size_t j = 0; j < function[i][3]/context->torus_base; j++){
        carpov_factorization(&lut_clear[(context->log_torus_base/log_carpov_base)*j], &luts[j], context->log_torus_base, log_carpov_base);
      }
      delete_IntPolynomial_array(number_of_luts, luts);
      break;

    case fbt_op_generate_b_LUT:
      generate_noiseless_trivial_LUT(&context->tlwe_temp[function[i][1]], lut_input, function[i][3], context->torus_base, tlwe_params);
      break;
    
    case fbt_op_tlwe_keyswitch:
      for (size_t j = 0; j < function[i][3]; j++){
        public_functional_tlweKeySwitch(&context->tlwe_temp[function[i][1] + j], context->tlweKS, &context->input[function[i][2] + j*context->torus_base], context->torus_base);
      }
      break;

    case fbt_op_load_literal:
      for (size_t j = 0; j < function[i][3]; j++){
        lweNoiselessTrivial(&context->output[function[i][1] + j], modSwitchToTorus32(function[i][2], 2*context->torus_base), extracted_params);
      }
      break;

    case fbt_op_load_hardcoded_LUT:
      lut_input = FBT_HARDCODED_LUT_LIST[function[i][2]];
      break;

    case fbt_op_trivial:
      fbt_op_trivial_FUNCTIONS[function[i][2]](function[i][1], function[i][3], context);
      break;

    case fbt_op_lwe_keyswitch:
      // TODO
      break;
      
    case fbt_op_generate_LUT_from_f:
      // TODO
      for (size_t j = 0; j < function[i][3]; j++)
      {
        
      }
      break;
      

    case fbt_op_next_tree_level:
      current_input+=1;
      temp = context->input;
      context->input = context->output;
      context->output = temp;
      break;
      
    default:
      break;
    }
  }
}

fbt_integer * fbt_new_encrypted_integer(int value, int digits, LweKey * key, fbt_context * context){
  fbt_integer * result = (fbt_integer*)malloc(sizeof(fbt_integer));
  result->lwe_samples = encryptInteger(value, digits, context->log_torus_base, key);
  result->digits = digits;
  result->lwe_params = key->params;
  result->log_torus_base = context->log_torus_base;
  return result;
}

fbt_integer * fbt_integer_new_array(int digits, int n, fbt_context * context){
  fbt_integer * result = (fbt_integer*)malloc(sizeof(fbt_integer) * n);
  for (size_t i = 0; i < n; i++){
    result[i].lwe_samples = new_LweSample_array(digits, &context->params->tgsw_params->tlwe_params->extracted_lweparams);
    assert(result[i].lwe_samples != 0);
    result[i].digits = digits;
    result[i].lwe_params =  &context->params->tgsw_params->tlwe_params->extracted_lweparams;
    result[i].log_torus_base = context->log_torus_base;
  }
  return result;
}

void fbt_integer_delete(fbt_integer * integer){
  delete_LweSample_array(integer->digits, integer->lwe_samples);
  free(integer);
}

void fbt_integer_encrypt(fbt_integer * output, int value, LweKey * key, fbt_context * context){
  const int torus_base = 1 << context->log_torus_base;
  const int mask = torus_base - 1;
  for (size_t i = 0; i < output->digits; i++){
    lweSymEncrypt(&output->lwe_samples[i], modSwitchToTorus32(value & mask, torus_base*2), key->params->alpha_min, key);
    value >>= context->log_torus_base;
  }
  output->lwe_params = key->params;
  output->log_torus_base = context->log_torus_base;
}

void fbt_integer_trivial(fbt_integer * output, int value){
  const int torus_base = 1 << output->log_torus_base;
  const int mask = torus_base - 1;
  for (size_t i = 0; i < output->digits; i++){
    lweNoiselessTrivial(&output->lwe_samples[i], modSwitchToTorus32(value & mask, torus_base*2), output->lwe_params);
    value >>= output->log_torus_base;
  }
  output->log_torus_base = output->log_torus_base;
}

void fbt_integer_trivial_with_scale(fbt_integer * output, int value, int scale){
  const int torus_base = 1 << output->log_torus_base;
  const int mask = torus_base - 1;
  for (size_t i = 0; i < output->digits; i++){
    lweNoiselessTrivial(&output->lwe_samples[i], modSwitchToTorus32(value & mask, scale*torus_base*2), output->lwe_params);
    value >>= output->log_torus_base;
  }
  output->log_torus_base = output->log_torus_base;
}

void fbt_integer_encrypt_with_scale(fbt_integer * output, int value, int scale, LweKey * key, fbt_context * context){
  const int torus_base = 1 << context->log_torus_base;
  const int mask = torus_base - 1;
  for (size_t i = 0; i < output->digits; i++){
    lweSymEncrypt(&output->lwe_samples[i], modSwitchToTorus32(value & mask, scale*torus_base*2), key->params->alpha_min, key);
    value >>= context->log_torus_base;
  }
  output->lwe_params = key->params;
  output->log_torus_base = context->log_torus_base;
}

int fbt_integer_decrypt_with_scale(fbt_integer * input, int scale, LweKey * key, fbt_context * context){
  int result = 0;
  for (int i = input->digits - 1; i >= 0; i--){
    result <<= context->log_torus_base;
    const int tmp = modSwitchFromTorus32(lweSymDecrypt(&input->lwe_samples[i], key, scale*context->torus_base*2), scale*context->torus_base*2);
    // printf("DEBUG: %d * %d^%d\n", tmp, torus_base, i);
    result += tmp;
  }
  // return result;
  return (result << (32 - context->log_torus_base*input->digits)) >> (32 - context->log_torus_base*input->digits);
}

void fbt_integer_bootstrap(fbt_integer * integer, fbt_context * context){
  const int torus_base = 1 << integer->log_torus_base;
  LweSample * input = new_LweSample(context->params->in_out_params);

  int lut[torus_base];
  for (int i = 0; i < torus_base; i++){
    lut[i] = i;
  }

  generate_noiseless_trivial_LUT(context->tlwe_temp, (const int *) lut, context->torus_base, context->torus_base, context->params->tgsw_params->tlwe_params);
  for (size_t j = 0; j < integer->digits; j++){
    lweKeySwitch(input, context->bootstrap_key->ks, &integer->lwe_samples[j]);
    tfhe_functional_bootstrap_woKS_FFT(&integer->lwe_samples[j], context->bootstrap_key, context->tlwe_temp, input, context->torus_base);
  }

  delete_LweSample(input);
}

void fbt_integer_bootstrap_with_scale(fbt_integer * integer, int scale, fbt_context * context){
  const int torus_base = 1 << integer->log_torus_base;
  LweSample * input = new_LweSample(context->params->in_out_params);

  int lut[torus_base];
  for (int i = 0; i < torus_base; i++){
    lut[i] = i;
  }

  generate_noiseless_trivial_LUT_with_scaling(context->tlwe_temp, (const int *) lut, context->torus_base, context->torus_base, scale, context->params->tgsw_params->tlwe_params);
  for (size_t j = 0; j < integer->digits; j++){
    lweKeySwitch(input, context->bootstrap_key->ks, &integer->lwe_samples[j]);
    tfhe_functional_bootstrap_woKS_FFT(&integer->lwe_samples[j], context->bootstrap_key, context->tlwe_temp, input, context->torus_base);
  }

  delete_LweSample(input);
}

int fbt_integer_decrypt(fbt_integer * input, LweKey * key, fbt_context * context){
  return decryptInteger(input->lwe_samples, input->digits, context->log_torus_base, key);
}

void fbt_integer_add(fbt_integer * output, fbt_integer * input1, fbt_integer * input2, fbt_context * context){
  const TFheGateBootstrappingParameterSet * params = context->params;
  const LweParams * extracted_params = &params->tgsw_params->tlwe_params->extracted_lweparams;
  const int torus_base = 1 << input1->log_torus_base;
  int lut[torus_base];
  for (int i = 0; i < torus_base; i++){
    lut[i] = -1;
  }
  generate_noiseless_trivial_LUT_with_scaling(context->tlwe_temp, (const int *) lut, torus_base, torus_base, 2, params->tgsw_params->tlwe_params);
  LweSample * sum = new_LweSample(params->in_out_params);


  lweNoiselessTrivial(&output->lwe_samples[0], 0, extracted_params);
  for (int i = 0; i < input1->digits; i++){
    lweAddTo(&output->lwe_samples[i], &input1->lwe_samples[i], &params->tgsw_params->tlwe_params->extracted_lweparams);
    lweAddTo(&output->lwe_samples[i], &input2->lwe_samples[i], &params->tgsw_params->tlwe_params->extracted_lweparams);
    lweKeySwitch(sum, context->bootstrap_key->ks, &output->lwe_samples[i]);
    tfhe_functional_bootstrap_with_multiple_extract(context->output, context->torus_base, context->bootstrap_key, context->tlwe_temp, sum, context->torus_base);
    for (int j = 0; j < context->torus_base; j++){
      lweSubTo(&output->lwe_samples[i], &context->output[i], extracted_params);
    }
    output->lwe_samples[i].b -= modSwitchToTorus32(2, context->torus_base*2);
    if(i != output->digits - 1){ // carry
      lweNoiselessTrivial(&output->lwe_samples[i + 1], modSwitchToTorus32(1, 2*context->torus_base*2), extracted_params);
      lweAddTo(&output->lwe_samples[i + 1], context->output, extracted_params);
    }
  }
  
  
  delete_LweSample(sum);
}

void fbt_integer_propagate_carry(fbt_integer * input, fbt_context * context){
  const TFheGateBootstrappingParameterSet * params = context->params;
  const LweParams * extracted_params = &params->tgsw_params->tlwe_params->extracted_lweparams;
  const int torus_base = 1 << input->log_torus_base;
  int lut[torus_base];
  for (int i = 0; i < torus_base; i++){
    lut[i] = -1;
  }
  generate_noiseless_trivial_LUT_with_scaling(context->tlwe_temp, (const int *) lut, torus_base, torus_base, 2, params->tgsw_params->tlwe_params);
  LweSample * sum = new_LweSample(params->in_out_params);


  for (int i = 0; i < input->digits; i++){
    lweKeySwitch(sum, context->bootstrap_key->ks, &input->lwe_samples[i]);
    tfhe_functional_bootstrap_with_multiple_extract(context->output, context->torus_base, context->bootstrap_key, context->tlwe_temp, sum, context->torus_base);
    for (int j = 0; j < context->torus_base; j++){
      lweSubTo(&input->lwe_samples[i], &context->output[i], extracted_params);
    }
    input->lwe_samples[i].b -= modSwitchToTorus32(2, context->torus_base*2);
    if(i != input->digits - 1){ // carry
      input->lwe_samples[i + 1].b +=  modSwitchToTorus32(1, 2*context->torus_base*2);
      lweAddTo(&input->lwe_samples[i + 1], context->output, extracted_params);
    }
  }
  delete_LweSample(sum);
}

void fbt_integer_add_to_wo_carry(fbt_integer * output, fbt_integer * input, fbt_context * context){
  for (size_t i = 0; i < input->digits; i++){
    lweAddTo(&output->lwe_samples[i], &input->lwe_samples[i], output->lwe_params);
  }
}

void fbt_integer_addMul_to_wo_carry(fbt_integer * output, fbt_integer * input, int value, fbt_context * context){
  const int mask = context->torus_base - 1;
  int dec_value = 0;
  for (size_t j = 0; j < input->digits; j++){
    dec_value = value&mask;
    for (size_t i = 0; i <= j; i++){
      lweAddMulTo(&output->lwe_samples[i], dec_value, &input->lwe_samples[context->torus_base - 1 + i - j], output->lwe_params);
    }
    value >>= context->log_torus_base;
  }
}

void fbt_integer_negate_wo_carry(fbt_integer * number, fbt_context * context){
  for (size_t i = 0; i < number->digits; i++){
    lweCopy(context->input, &number->lwe_samples[i], number->lwe_params);
    lweNoiselessTrivial(&number->lwe_samples[i], modSwitchToTorus32(context->torus_base - 1, context->torus_base*2), number->lwe_params);
    lweSubTo(&number->lwe_samples[i], context->input, number->lwe_params);
  }
  number->lwe_samples[0].b += 1; // 2-complement conversion
}

void fbt_create_LUT_from_literals(fbt_lut * lut, const int * literals, int num_of_numbers, fbt_context * context){
  lut->literal_promisse = literals;
  lut->size = num_of_numbers/context->torus_base;
}

void fbt_create_LUT_from_lwe_samples(fbt_lut * lut, LweSample * samples, int num_of_samples, fbt_context * context){
  for (size_t j = 0; j < num_of_samples/context->torus_base; j++){
    public_functional_tlweKeySwitch(&lut->tlwe_samples[j], context->tlweKS, &samples[j*context->torus_base], context->torus_base);
  }
  lut->size = num_of_samples/context->torus_base;
}

void fbt_f_bootstrap_LUT(LweSample * output, LweSample * input, fbt_lut * lut, fbt_context * context){
  if(lut->tlwe_samples){
    for (size_t j = 0; j < lut->size; j++){
      tfhe_functional_bootstrap_woKS_FFT(&output[j], context->bootstrap_key, &lut->tlwe_samples[j], input, context->torus_base);
    }
  }else if(lut->literal_promisse){
    generate_noiseless_trivial_LUT(context->tlwe_temp, lut->literal_promisse, lut->size*context->torus_base, context->torus_base, context->params->tgsw_params->tlwe_params);
    for (size_t j = 0; j < lut->size; j++){
      tfhe_functional_bootstrap_woKS_FFT(&output[j], context->bootstrap_key, &context->tlwe_temp[j], input, context->torus_base);
    }
  }else{
    printf("Error: undefined LUT");
    abort();
  }
}


void fbt_mux(fbt_integer * output, LweSample * select, fbt_integer * input, fbt_context * context){

  for (size_t i = 0; i < input->digits; i++){
    for (size_t j = 0; j < context->torus_base; j++){
      lweCopy(&context->input[j], &input[j].lwe_samples[i], &context->params->tgsw_params->tlwe_params->extracted_lweparams);
    }
    public_functional_tlweKeySwitch(context->tlwe_temp, context->tlweKS, context->input, context->torus_base);
    tfhe_functional_bootstrap_woKS_FFT(&output->lwe_samples[i], context->bootstrap_key, context->tlwe_temp, select, context->torus_base);
  }
}


void fbt_run_bootstrap_level(){

}
