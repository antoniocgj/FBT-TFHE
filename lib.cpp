#include "functional_bootstrap.h"

uint64_t __g_clock_begin, __g_clock_end;
#define MEASURE_TIME_AND_PRINT(X) \
  __g_clock_begin = get_time(); \
  X; \
  __g_clock_end = get_time(); \
  printf("Time: %lu,%03lu,%03lu μs\n", (__g_clock_end - __g_clock_begin)/1000000, ((__g_clock_end - __g_clock_begin)/1000)%1000, (__g_clock_end - __g_clock_begin)%1000);



void generate_random_bytes(uint64_t amount, uint8_t * pointer){
  unsigned long long rnd;
  int i = 0;
  while (i < amount) {
    if( i%8 == 0 && 0 == _rdrand64_step (&rnd)){
      printf("Random Generation Failed\n");
      continue;
    }
    pointer[i] = ((uint8_t *) &rnd)[i%8];
    i++;
  }
}

TFheGateBootstrappingParameterSet *new_TFHE_parameters(const int32_t bk_l, const int32_t bk_Bgbit, const int32_t ks_basebit, const int32_t ks_t) {
  static const int32_t N = 1024;
  static const int32_t k = 1;
  static const int32_t n = 630;
  // static const int32_t bk_l = 4;
  // static const int32_t bk_Bgbit = 5; 
  // static const int32_t ks_basebit = 5;
  // static const int32_t ks_length = 2; // t
  static const double ks_stdev = pow(2, -15); //standard deviation
  static const double bk_stdev = pow(2, -25); //standard deviation
  static const double max_stdev = 0.012467; //max standard deviation for a 1/4 msg space

  double a = n*(k+1)*N*bk_l*pow((1 << bk_Bgbit)/2, 2)*pow(bk_stdev, 2); // _n(¯k+1)¯l¯N¯β^2¯alpha^2
  double b = n*(k*N+1)*pow(1/(2*pow((1 << bk_Bgbit), bk_l)), 2); // _n(1 + ¯k¯N)¯epsilon^2
  double c = N*pow(1 << ks_basebit, -2*(ks_t+1)); // ¯n*base^(−2(t+1))
  double d = N*ks_t*pow(ks_stdev, 2); // t¯n*ks_stdev^2
  // ¯n = 1024, β = Bg/2, Bg = 1 << bk_Bgbit
  // printf("Error <= %12.5e + %12.5e + %12.5e + %12.5e = %12.5e\n", a, b, c, d, a+b+c+d);
  
  LweParams *params_in = new_LweParams(n, ks_stdev, max_stdev);
  TLweParams *params_accum = new_TLweParams(N, k, bk_stdev, max_stdev);
  TGswParams *params_bk = new_TGswParams(bk_l, bk_Bgbit, params_accum);

  TfheGarbageCollector::register_param(params_in);
  TfheGarbageCollector::register_param(params_accum);
  TfheGarbageCollector::register_param(params_bk);

  return new TFheGateBootstrappingParameterSet(ks_t, ks_basebit, params_in, params_bk);
}
 

TFheGateBootstrappingParameterSet *parameters_setup() {
  static const int32_t N = 1024;
  static const int32_t k = 1;
  static const int32_t n = 630;
  static const int32_t bk_l = 4;
  static const int32_t bk_Bgbit = 5; 
  static const int32_t ks_basebit = 6;
  static const int32_t ks_length = 2; // t
  static const double ks_stdev = pow(2, -15); //standard deviation
  static const double bk_stdev = pow(2, -25); //standard deviation
  static const double max_stdev = 0.012467; //max standard deviation for a 1/4 msg space

  double a = n*(k+1)*N*bk_l*pow((1 << bk_Bgbit)/2, 2)*pow(bk_stdev, 2); // _n(¯k+1)¯l¯N¯β^2¯alpha^2
  double b = n*(k*N+1)*pow(1/(2*pow((1 << bk_Bgbit), bk_l)), 2); // _n(1 + ¯k¯N)¯epsilon^2
  double c = N*pow(1 << ks_basebit, -2*(ks_length+1)); // ¯n*base^(−2(t+1))
  double d = N*ks_length*pow(ks_stdev, 2); // t¯n*ks_stdev^2
  // ¯n = 1024, β = Bg/2, Bg = 1 << bk_Bgbit
  // printf("Error <= %12.5e + %12.5e + %12.5e + %12.5e = %12.5e\n", a, b, c, d, a+b+c+d);
  
  LweParams *params_in = new_LweParams(n, ks_stdev, max_stdev);
  TLweParams *params_accum = new_TLweParams(N, k, bk_stdev, max_stdev);
  TGswParams *params_bk = new_TGswParams(bk_l, bk_Bgbit, params_accum);

  TfheGarbageCollector::register_param(params_in);
  TfheGarbageCollector::register_param(params_accum);
  TfheGarbageCollector::register_param(params_bk);

  return new TFheGateBootstrappingParameterSet(ks_length, ks_basebit, params_in, params_bk);
}

LweSample * encryptInteger(int number, int digits, const int log_torus_base, LweKey * key){
  LweSample * result = new_LweSample_array(digits, key->params);
  const int torus_base = 1 << log_torus_base;
  const int mask = torus_base - 1;
  for (size_t i = 0; i < digits; i++){
    lweSymEncrypt(&result[i], modSwitchToTorus32(number & mask, torus_base*2), key->params->alpha_min, key);
    number >>= log_torus_base;
  }
  return result;
}


int decryptInteger(LweSample * sample, int digits, const int log_torus_base, LweKey * key){
  int result = 0;
  const int torus_base = 1 << log_torus_base;
  for (int i = digits - 1; i >= 0; i--){
    result <<= log_torus_base;
    const int tmp = modSwitchFromTorus32(lweSymDecrypt(&sample[i], key, torus_base*2), torus_base*2);
    // printf("DEBUG: %d * %d^%d\n", tmp, torus_base, i);
    result += tmp;
  }
  // return result;
  return (result << (32 - log_torus_base*digits)) >> (32 - log_torus_base*digits);
}

typedef struct TFheFunctionalBootstrappingSecretKeySet{
  TFheGateBootstrappingSecretKeySet* tfhe_keys;
  const TFheGateBootstrappingCloudKeySet* cloud_key;
  LweKey *extracted_key;
  TLweKeySwitchKey * tlweKS;
  int log_torus_base;
} TFheFunctionalBootstrappingSecretKeySet;

TFheFunctionalBootstrappingSecretKeySet* new_random_functional_bootstrapping_secret_keyset(TFheGateBootstrappingParameterSet* params, const int log_torus_base){
  TFheFunctionalBootstrappingSecretKeySet * result = (TFheFunctionalBootstrappingSecretKeySet *) malloc(sizeof(TFheFunctionalBootstrappingSecretKeySet));
  uint32_t seed[8];
  generate_random_bytes(32, (uint8_t * ) seed);
  result->log_torus_base = log_torus_base;
  const int torus_base = 1 << log_torus_base;
  tfhe_random_generator_setSeed(seed, 8);
  result->tfhe_keys = new_random_gate_bootstrapping_secret_keyset(params);
  result->cloud_key = &result->tfhe_keys->cloud;
  const TLweKey *tlwe_key = &result->tfhe_keys->tgsw_key->tlwe_key;
  result->extracted_key = new_LweKey(&result->tfhe_keys->tgsw_key->tlwe_params->extracted_lweparams);
  tLweExtractKey(result->extracted_key, tlwe_key);  
  result->tlweKS = new_TLweKeySwitchKey(result->tfhe_keys->cloud.bkFFT->accum_params->N, result->tfhe_keys->cloud.bkFFT->accum_params->k, params->ks_t, params->ks_basebit, result->tfhe_keys->cloud.bkFFT->bk_params->tlwe_params, torus_base);
  tLweCreateKeySwitchKey(result->tlweKS, result->extracted_key, tlwe_key, torus_base);
  return result;
}

void export_tfhe_functional_bootstrapping_public_keyset_to_file(FILE* F, const TFheFunctionalBootstrappingSecretKeySet* key){
  export_tfheGateBootstrappingCloudKeySet_toFile(F, key->cloud_key);
  fwrite(&key->log_torus_base, 4, 1, F);
  const int ks_base = 1<<key->tlweKS->basebit;
  const int torus_base = 1 << key->log_torus_base;
  for (size_t i = 0; i < (ks_base - 1)*key->tlweKS->N*key->tlweKS->t*torus_base; i++){
    const TLweSample * sample = &key->tlweKS->ks0_raw[i];
    fwrite(sample->a->coefsT, sizeof(Torus32), key->tlweKS->N, F);
    fwrite(sample->b->coefsT, sizeof(Torus32), key->tlweKS->N, F);
  }
}

void export_tfhe_functional_bootstrapping_secret_keyset_to_file(FILE* F, const TFheFunctionalBootstrappingSecretKeySet* key){
  export_tfheGateBootstrappingSecretKeySet_toFile(F, key->tfhe_keys);
  fwrite(&key->log_torus_base, 4, 1, F);
  const int ks_base = 1<<key->tlweKS->basebit;
  const int torus_base = 1 << key->log_torus_base;
  for (size_t i = 0; i < (ks_base - 1)*key->tlweKS->N*key->tlweKS->t*torus_base; i++){
    const TLweSample * sample = &key->tlweKS->ks0_raw[i];
    fwrite(sample->a->coefsT, sizeof(Torus32), key->tlweKS->N, F);
    fwrite(sample->b->coefsT, sizeof(Torus32), key->tlweKS->N, F);
  }
}

TFheFunctionalBootstrappingSecretKeySet* new_functional_bootstrapping_public_keyset_from_file(FILE * F){
  TFheFunctionalBootstrappingSecretKeySet * result = (TFheFunctionalBootstrappingSecretKeySet *) malloc(sizeof(TFheFunctionalBootstrappingSecretKeySet));
  result->cloud_key = new_tfheGateBootstrappingCloudKeySet_fromFile(F);
  int read = fread((void *)&result->log_torus_base, 4, 1, F);
  if(read != 4) printf("Error while reading log_torus_base\n");
  const TLweParams * tlwe_params = result->cloud_key->params->tgsw_params->tlwe_params;
  const int torus_base = 1 << result->log_torus_base;
  const int ks_base = 1<< result->cloud_key->params->ks_basebit;
  result->tlweKS = new_TLweKeySwitchKey(tlwe_params->N, tlwe_params->k, result->cloud_key->params->ks_t, result->cloud_key->params->ks_basebit, tlwe_params, torus_base);
  for (size_t i = 0; i < (ks_base - 1)*result->tlweKS->N*result->tlweKS->t*torus_base; i++){
    const TLweSample * sample = &result->tlweKS->ks0_raw[i];
    read = fread((void *)sample->a->coefsT, sizeof(Torus32), result->tlweKS->N, F);
    read += fread((void *)sample->b->coefsT, sizeof(Torus32), result->tlweKS->N, F);
    if(read != result->tlweKS->N*sizeof(Torus32)*2) printf("Error while reading TLwe KS key\n");
  }
  return result;
}

TFheFunctionalBootstrappingSecretKeySet* new_functional_bootstrapping_secret_keyset_from_file(FILE * F){
  uint32_t seed[8];
  generate_random_bytes(32, (uint8_t * ) seed);
  tfhe_random_generator_setSeed(seed, 8);
  TFheFunctionalBootstrappingSecretKeySet * result = (TFheFunctionalBootstrappingSecretKeySet *) malloc(sizeof(TFheFunctionalBootstrappingSecretKeySet));
  result->tfhe_keys = new_tfheGateBootstrappingSecretKeySet_fromFile(F);
  result->cloud_key = &result->tfhe_keys->cloud;
  int read = fread((void *)&result->log_torus_base, 4, 1, F);
  if(read != 1) printf("Error while reading log_torus_base: %d\n", result->log_torus_base);
  const TLweParams * tlwe_params = result->cloud_key->bkFFT->bk_params->tlwe_params;
  const int torus_base = 1 << result->log_torus_base;
  const int ks_base = 1<< result->cloud_key->params->ks_basebit;
  result->tlweKS = new_TLweKeySwitchKey(tlwe_params->N, tlwe_params->k, result->cloud_key->params->ks_t, result->cloud_key->params->ks_basebit, tlwe_params, torus_base);
  for (size_t i = 0; i < (ks_base - 1)*result->tlweKS->N*result->tlweKS->t*torus_base; i++){
    const TLweSample * sample = &result->tlweKS->ks0_raw[i];
    read = fread((void *)sample->a->coefsT, sizeof(Torus32), result->tlweKS->N, F);
    read += fread((void *)sample->b->coefsT, sizeof(Torus32), result->tlweKS->N, F);
    if(read != result->tlweKS->N*2) printf("Error while reading TLwe KS key\n");
  }
  result->extracted_key = new_LweKey(&result->tfhe_keys->tgsw_key->tlwe_params->extracted_lweparams);
  tLweExtractKey(result->extracted_key, &result->tfhe_keys->tgsw_key->tlwe_key);  
  return result;
}


void decryptAndPrint(LweSample * ct, const LweKey * key, int Msize){
  Torus32 ai = lweSymDecrypt(ct, key, Msize);
  Torus32 phase = lwePhase(ct, key);
  printf("Phase: %lf, Round(Phase): %lf, Torus32: 0x%08x, Round(Torus32): 0x%08x, ModSwitch: %d\n", t32tod(phase), t32tod(ai), phase, ai, modSwitchFromTorus32(phase, Msize));
}

void decryptAndPrintPolynomial(TLweSample * ct, const TLweKey * key, int Msize){
  TorusPolynomial * result = new_TorusPolynomial(ct->b->N);
  tLweSymDecrypt(result, ct, key, Msize);
  // tLwePhase(result, ct, key);
  for (size_t i = 0; i < ct->b->N; i++)
  {
    printf("%d (%lf), ", modSwitchFromTorus32(result->coefsT[i], Msize), t32tod(result->coefsT[i]));
  }
  printf("\n");
  delete_TorusPolynomial(result);
}

void decryptAndPrintPolynomialNoRepeat(TLweSample * ct, const TLweKey * key, int Msize){
  TorusPolynomial * result = new_TorusPolynomial(ct->b->N);
  tLweSymDecrypt(result, ct, key, Msize);
  // tLwePhase(result, ct, key);
  int last, at, repeatedCount = 1;
  at = modSwitchFromTorus32(result->coefsT[0], Msize);
  printf("%d (%lf) ", at, t32tod(result->coefsT[0]));
  last = at;
  for (size_t i = 1; i < ct->b->N; i++)
  {
    at = modSwitchFromTorus32(result->coefsT[i], Msize);
    if(at == last){
      repeatedCount++;
      continue;
    }
    printf("x%d, %d (%lf) ", repeatedCount, at, t32tod(result->coefsT[i]));
    repeatedCount = 1;
    last = at;
  }
  printf("x%d \n", repeatedCount);
  delete_TorusPolynomial(result);
}


/**
 * result = LWE(v_p) where p=barb-sum(bara_i.s_i) mod 2N
 * @param result the output LWE sample
 * @param v a 2N-elt anticyclic function (represented by a TorusPolynomial)
 * @param bk An array of n TGSW FFT samples where bk_i encodes s_i
 * @param barb A coefficients between 0 and 2N-1
 * @param bara An array of n coefficients between 0 and 2N-1
 * @param bk_params The parameters of bk
 */
void tfhe_functional_blindRotateAndExtract_FFT(LweSample *result,
                                           const TLweSample *v,
                                           const TGswSampleFFT *bk,
                                           const int32_t barb,
                                           const int32_t *bara,
                                           const int32_t n,
                                           const TGswParams *bk_params) {

    const TLweParams *accum_params = bk_params->tlwe_params;
    const LweParams *extract_params = &accum_params->extracted_lweparams;
    const int32_t N = accum_params->N;
    const int32_t _2N = 2 * N;

    // Test polynomial 
    TorusPolynomial *testvectbis = new_TorusPolynomial(N);
    // Accumulator
    TLweSample *acc = new_TLweSample(accum_params);

    // testvector = X^{2N-barb}*v
    if (barb != 0){
      for (size_t i = 0; i <= accum_params->k; i++){
        torusPolynomialMulByXai(&acc->a[i], _2N - barb, &v->a[i]);
      }
    }else{
      for (size_t i = 0; i <= accum_params->k; i++){
        torusPolynomialCopy(&acc->a[i], &v->a[i]);
      }
    } 
    // Blind rotation
    tfhe_blindRotate_FFT(acc, bk, bara, n, bk_params);
    // Extraction
    tLweExtractLweSample(result, acc, extract_params, accum_params);

    delete_TLweSample(acc);
    delete_TorusPolynomial(testvectbis);
}

/**
 * result = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
 * @param result The resulting LweSample
 * @param bk The bootstrapping + keyswitch key
 * @param mu The output message (if phase(x)>0)
 * @param x The input sample
 */
void tfhe_functional_bootstrap_woKS_FFT(LweSample *result,
                                    const LweBootstrappingKeyFFT *bk,
                                    TLweSample * testvect,
                                    const LweSample *x, const int torus_base) {

  const TGswParams *bk_params = bk->bk_params;
  const TLweParams *accum_params = bk->accum_params;
  const LweParams *in_params = bk->in_out_params;
  const int32_t N = accum_params->N;
  const int32_t Nx2 = 2 * N;
  const int32_t n = in_params->n;

  // precision offset
  const Torus32 prec_offset = modSwitchToTorus32(N/(2*torus_base), Nx2);

  int32_t *bara = new int32_t[N];


  // Modulus switching
  int32_t barb = modSwitchFromTorus32(x->b + prec_offset, Nx2);
  for (int32_t i = 0; i < n; i++) {
      bara[i] = modSwitchFromTorus32(x->a[i], Nx2);
  }

  // the initial testvec = [mu,mu,mu,...,mu]
  // for (int32_t i = 0; i < N; i++) testvect->coefsT[i] = mu[i];

  // Bootstrapping rotation and extraction
  tfhe_functional_blindRotateAndExtract_FFT(result, testvect, bk->bkFFT, barb, bara, n, bk_params);


  delete[] bara;

}

void tfhe_functional_bootstrap_with_multiple_extract(LweSample *result, int num_of_extracts, const LweBootstrappingKeyFFT *bk, TLweSample * testvect, const LweSample *x, const int torus_base) {
  const TGswParams *bk_params = bk->bk_params;
  const LweParams *in_params = bk->in_out_params;
  const int32_t n = in_params->n;
  const TLweParams *accum_params = bk_params->tlwe_params;
  const LweParams *extract_params = &accum_params->extracted_lweparams;
  const int32_t N = accum_params->N;
  const int32_t _2N = 2 * N;

  const Torus32 prec_offset = modSwitchToTorus32(N/(2*torus_base), _2N);

  int32_t *bara = new int32_t[N];


  // Modulus switching
  int32_t barb = modSwitchFromTorus32(x->b + prec_offset, _2N);
  for (int32_t i = 0; i < n; i++) {
      bara[i] = modSwitchFromTorus32(x->a[i], _2N);
  }

  // Test polynomial 
  TorusPolynomial *testvectbis = new_TorusPolynomial(N);
  // Accumulator
  TLweSample *acc = new_TLweSample(accum_params);

  // testvector = X^{2N-barb}*v
  if (barb != 0){
    for (size_t i = 0; i <= accum_params->k; i++){
      torusPolynomialMulByXai(&acc->a[i], _2N - barb, &testvect->a[i]);
    }
  }else{
    for (size_t i = 0; i <= accum_params->k; i++){
      torusPolynomialCopy(&acc->a[i], &testvect->a[i]);
    }
  } 
  // Blind rotation
  tfhe_blindRotate_FFT(acc, bk->bkFFT, bara, n, bk_params);
  // Extraction
  LweSample * res_temp = new_LweSample(extract_params);
  for (size_t i = 0; i < num_of_extracts/2; i++){
    tLweExtractLweSampleIndex(&result[i], acc, i, extract_params, accum_params);
  }
  for (size_t i = num_of_extracts/2; i < num_of_extracts; i++){
    tLweExtractLweSampleIndex(res_temp, acc, N - 1 - (i - num_of_extracts/2), extract_params, accum_params);
    lweNegate(&result[i], res_temp, extract_params);
  }
  
  // tLweExtractLweSample(result, acc, extract_params, accum_params);

  delete_LweSample(res_temp);
  delete_TLweSample(acc);
  delete_TorusPolynomial(testvectbis);
  delete[] bara;
}

void tfhe_functional_bootstrap_wo_extract(TLweSample *result, const LweBootstrappingKeyFFT *bk, TLweSample * testvect, const LweSample *x, const int torus_base) {
  const TGswParams *bk_params = bk->bk_params;
  const LweParams *in_params = bk->in_out_params;
  const int32_t n = in_params->n;
  const TLweParams *accum_params = bk_params->tlwe_params;
  const LweParams *extract_params = &accum_params->extracted_lweparams;
  const int32_t N = accum_params->N;
  const int32_t _2N = 2 * N;

  const Torus32 prec_offset = modSwitchToTorus32(N/(2*torus_base), _2N);

  int32_t *bara = new int32_t[N];


  // Modulus switching
  int32_t barb = modSwitchFromTorus32(x->b + prec_offset, _2N);
  for (int32_t i = 0; i < n; i++) {
      bara[i] = modSwitchFromTorus32(x->a[i], _2N);
  }
  
  // The accumulator is the result


  // testvector = X^{2N-barb}*v
  if (barb != 0){
    for (size_t i = 0; i <= accum_params->k; i++){
      torusPolynomialMulByXai(&result->a[i], _2N - barb, &testvect->a[i]);
    }
  }else{
    for (size_t i = 0; i <= accum_params->k; i++){
      torusPolynomialCopy(&result->a[i], &testvect->a[i]);
    }
  } 
  // Blind rotation
  tfhe_blindRotate_FFT(result, bk->bkFFT, bara, n, bk_params);
  delete[] bara;
}


/**
 * result = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
 * @param result The resulting LweSample
 * @param bk The bootstrapping + keyswitch key
 * @param mu The output message (if phase(x)>0)
 * @param x The input sample
 */
void tfhe_functional_bootstrap_FFT(LweSample *result,
                               const LweBootstrappingKeyFFT *bk,
                               Torus32 * mu,
                               const LweSample *x, const int torus_base) {

  LweSample *u = new_LweSample(&bk->accum_params->extracted_lweparams);

  const int32_t N = bk->accum_params->N;
  TLweSample * testvect = new_TLweSample(bk->bk_params->tlwe_params);

  for (int32_t i = 0; i < N; i++){
    testvect->a[0].coefsT[i] = 0;
    testvect->b->coefsT[i] = mu[i];
  } 

  tfhe_functional_bootstrap_woKS_FFT(u, bk, testvect, x, torus_base);
  
  
  // Key switching
  lweKeySwitch(result, bk->ks, u);

  delete_LweSample(u);
  delete_TLweSample(testvect);
}

void tfhe_carpov_multivalue_bootstrap_init(TLweSample *result, const LweBootstrappingKeyFFT *bk, const LweSample *x, const int torus_base){
  const TGswParams *bk_params = bk->bk_params;
  const TLweParams *accum_params = bk->accum_params;
  const LweParams *in_params = bk->in_out_params;
  const int32_t N = accum_params->N;
  const int32_t Nx2 = 2 * N;
  const int32_t n = in_params->n;

  // precision offset
  const Torus32 prec_offset = modSwitchToTorus32(N/(2*torus_base), Nx2);

  TorusPolynomial *testvect = new_TorusPolynomial(N);
  for (int32_t i = 0; i < N; i++) testvect->coefsT[i] = modSwitchToTorus32(N/torus_base, 2*Nx2);
  
  int32_t *bara = new int32_t[N];
  // Modulus switching
  int32_t barb = modSwitchFromTorus32(x->b + prec_offset, Nx2);
  for (int32_t i = 0; i < n; i++) {
    bara[i] = modSwitchFromTorus32(x->a[i], Nx2);
  }

  accum_params = bk_params->tlwe_params; ///< Params of each row
  const LweParams *extract_params = &accum_params->extracted_lweparams; 

  TorusPolynomial *testvectbis = new_TorusPolynomial(N);

  if (barb != 0) torusPolynomialMulByXai(testvectbis, Nx2 - barb, testvect);
  else torusPolynomialCopy(testvectbis, testvect);
  tLweNoiselessTrivial(result, testvectbis, accum_params);
  tfhe_blindRotate_FFT(result, bk->bkFFT, bara, n, bk_params);

  delete_TorusPolynomial(testvectbis);
  delete[] bara;
  delete_TorusPolynomial(testvect);
}

void carpov_factorization(IntPolynomial * TV_1, IntPolynomial *TV, const int log_torus_base, const int log_carpov_base){
  const int bitMask = (1<<log_carpov_base) - 1;
  
  for (size_t j = 0; j < log_torus_base; j+=log_carpov_base){
    TV_1[j/log_carpov_base].coefs[0] = ((TV->coefs[0] >> j) & bitMask) + ((TV->coefs[TV->N - 1] >> j) & bitMask);
  }
  
  for (size_t i = 1; i < TV->N; i++){
    for (size_t j = 0; j < log_torus_base; j+=log_carpov_base){
      TV_1[j/log_carpov_base].coefs[i] = ((TV->coefs[i] >> j) & bitMask) - ((TV->coefs[i-1] >> j) & bitMask);
    }
  }
}

void extract_and_keyswitch(LweSample *result, TLweSample * in, const LweBootstrappingKeyFFT *bk){
  LweSample *u = new_LweSample(&bk->accum_params->extracted_lweparams);

  tLweExtractLweSample(u, in, &bk->bk_params->tlwe_params->extracted_lweparams, bk->bk_params->tlwe_params);
  // Key switching
  lweKeySwitch(result, bk->ks, u);
  delete_LweSample(u);
}

void carpov_second_phase(LweSample * results, IntPolynomial * tv_f, TLweSample * rotated_tv_0, int q, const TLweParams * params){
  TLweSample * result_tlwe = new_TLweSample(params);
  LweSample * u = new_LweSample(&params->extracted_lweparams);

  for (size_t i = 0; i < q; i++){
    tLweNoiselessTrivialT(result_tlwe, (const Torus32) 0, params);
    tLweAddMulRTo(result_tlwe, &tv_f[i], rotated_tv_0, params);
    tLweExtractLweSample(&results[i], result_tlwe, &params->extracted_lweparams, params);
  }

  delete_TLweSample(result_tlwe);
  delete_LweSample(u);
}

void carpov_second_phase_with_base_composition(LweSample * results, IntPolynomial * tv_f, TLweSample * rotated_tv_0, int q, const TLweParams * params, int log_torus_base, const int log_carpov_base){
  TLweSample * result_tlwe = new_TLweSample(params);
  LweSample * res_temp = new_LweSample(&params->extracted_lweparams);

  const int number_of_elements = log_torus_base/log_carpov_base;

  for (size_t i = 0; i < q; i++){
    lweNoiselessTrivial(&results[i], 0, &params->extracted_lweparams);
    for (size_t j = 0; j < number_of_elements; j++){
      tLweNoiselessTrivialT(result_tlwe, (const Torus32) 0, params);
      tLweAddMulRTo(result_tlwe, &tv_f[i*number_of_elements + j], rotated_tv_0, params);
      for (size_t s = 0; s < (1<<(j*log_carpov_base)); s++){
        tLweExtractLweSampleIndex(res_temp, result_tlwe, s, &params->extracted_lweparams, params);
        lweAddTo(&results[i], res_temp, &params->extracted_lweparams);
      }
    }    
  }
  delete_TLweSample(result_tlwe);
  delete_LweSample(res_temp);
}


void carpov_base_compose(LweSample * result, LweSample * in, const LweParams * params, int number_of_elements, const int log_carpov_base){
  lweNoiselessTrivial(result, 0, params);
  for (size_t i = 0; i < number_of_elements; i++){
    lweAddMulTo(result, 1<<(i*log_carpov_base), &in[i], params);
  }
}

// TRLWE KeySwitch considering bases

void tLweKeySwitchTranslate_fromArray(TLweSample* result, 
	TLweSample**** ks, const TLweParams* params, 
	const TorusPolynomial* ai, 
	const int32_t n, const int32_t t, const int32_t basebit, const int torus_base){
  const int32_t base=1<<basebit;       // base=2 in [CGGI16]
  const int32_t prec_offset=1<<(32-(1+basebit*t)); //precision
  const int32_t mask=base-1;

  for (int32_t i=0;i<n;i++){
    for (size_t e = 0; e < torus_base; e++){
      const uint32_t aibar=ai[i].coefsT[e]+prec_offset;
      for (int32_t j=0;j<t;j++){
        const uint32_t aij=(aibar>>(32-(j+1)*basebit)) & mask;
        if(aij != 0) {
          for (int q = 0; q <= 1; q++) {
            for (int p = 0; p < n; p++)
              result->a[q].coefsT[p] -= ks[i][e][j][aij - 1].a[q].coefsT[p];                
          }
        }
      }
    }     
  }
}

void public_functional_tlweKeySwitch(TLweSample* result, const TLweKeySwitchKey* ks, LweSample* samples, const int torus_base){
  const TLweParams* params=ks->out_params;
  const int32_t n=ks->N;
  const int32_t basebit=ks->basebit;
  const int32_t t=ks->t;
  TorusPolynomial * f_a = new_TorusPolynomial_array(n, torus_base);

  for (size_t i = 0; i < n; i++) {
    result->a[0].coefsT[i] = 0;
    result->b->coefsT[i] = samples[i/(n/torus_base)].b;
  }

  for (size_t i = 0; i < n; i++){ 
    for (size_t j = 0; j < torus_base; j++){
      f_a[i].coefsT[j] = samples[j].a[i];
    }
  }

  tLweKeySwitchTranslate_fromArray(result, ks->ks, params, f_a, n, t, basebit, torus_base);

  delete_TorusPolynomial_array(n, f_a);
}


TLweKeySwitchKey * new_TLweKeySwitchKey(int32_t N, int32_t k, int32_t t, int32_t basebit, const TLweParams * out_params, const int torus_base){
  TLweKeySwitchKey * tlweKS = (TLweKeySwitchKey *) malloc(sizeof(TLweKeySwitchKey));
  tlweKS->N = N;
  tlweKS->k = k;
  tlweKS->t = t;
  tlweKS->basebit = basebit;
  tlweKS->out_params = out_params;
  int32_t base = 1<<basebit;

  double c = N*pow(base, -2*(t+1)); // ¯n*base^(−2(t+1))
  double d = torus_base*N*t*pow(out_params->alpha_min, 2); // t¯n*ks_stdev^2
  // ¯n = 1024, β = Bg/2, Bg = 1 << bk_Bgbit
  // printf("KS Error <= %12.5e + %12.5e = %12.5e\n", c, d, c+d);

  tlweKS->ks0_raw = new_TLweSample_array((base - 1)*N*t*torus_base, out_params);
  tlweKS->ks1_raw = new TLweSample*[N*t*torus_base];
  tlweKS->ks2_raw = new TLweSample**[N*torus_base];
  tlweKS->ks = new TLweSample***[N];
   
  for (int32_t p = 0; p < N*t*torus_base; ++p)
	  tlweKS->ks1_raw[p] = tlweKS->ks0_raw + (base - 1)*p;
	for (int32_t p = 0; p < N*torus_base; ++p)
	  tlweKS->ks2_raw[p] = tlweKS->ks1_raw + t*p;
	for (int32_t p = 0; p < N; ++p)
	  tlweKS->ks[p] = tlweKS->ks2_raw + torus_base*p;


  return tlweKS;
}

void tLweCreateKeySwitchKey(TLweKeySwitchKey* result, const LweKey* in_key, const TLweKey* out_key, const int torus_base){
  const int32_t n = result->N;
  const int32_t t = result->t;
  const int32_t basebit = result->basebit;
  const int32_t base = 1<<basebit;
  const double alpha = out_key->params->alpha_min;
  const int32_t sizeks = n*t*(base-1);
  //const int32_t n_out = out_key->params->n;

  // double err = 0;

  // // chose a random vector of gaussian noises
  // double* noise = new double[sizeks];
  // for (int32_t i = 0; i < sizeks; ++i){
  //     normal_distribution<double> distribution(0.,alpha); 
  //     noise[i] = distribution(generator);
  //     err += noise[i];
  // }
  // // recenter the noises
  // err = err/sizeks;
  // for (int32_t i = 0; i < sizeks; ++i) noise[i] -= err;

  // generate the ks
  int32_t index = 0; 
  for (int32_t i = 0; i < n; ++i) {
    for (size_t e = 0; e < torus_base; e++){
      for (int32_t j = 0; j < t; ++j){
        for (int32_t h = 1; h < base; ++h) { 
          Torus32 mess = (in_key->key[i]*h)*(1<<(32-(j+1)*basebit));
          tLweSymEncryptZero(&result->ks[i][e][j][h - 1], alpha, out_key);
          for (size_t q = e*(n/torus_base); q < (e+1)*(n/torus_base); q++) result->ks[i][e][j][h - 1].b->coefsT[q] += mess;
        }
      }
    }
  }
  // delete[] noise; 
}

// Util functions

int new_noiseless_trivial_LUT(TLweSample ** result, int * LUT, int size, int torus_base, const TLweParams * params){
  *result = new_TLweSample_array(size/torus_base, params);
  const int N = params->N;
  const int slice_size = N/torus_base;
  for (size_t i = 0; i < size/torus_base; i++){
    for (size_t j = 0; j < N; j++){
      (*result)[i].a->coefsT[j] = 0;
      (*result)[i].b->coefsT[j] = modSwitchToTorus32(LUT[i*torus_base + (j/slice_size)]*slice_size, 2*N);
    }
  }
  return (size/torus_base);
}

int new_encrypted_LUT_from_literals(TLweSample ** result, int * LUT, int size, int torus_base, const TLweKey * key, const TLweParams * params){
  *result = new_TLweSample_array(size/torus_base, params);
  const int N = params->N;
  const int slice_size = N/torus_base;
  for (size_t i = 0; i < size/torus_base; i++){
    tLweSymEncryptZero(*result, key->params->alpha_min, key);
    for (size_t j = 0; j < N; j++){
      (*result)[i].b->coefsT[j] += modSwitchToTorus32(LUT[i*torus_base + (j/slice_size)]*slice_size, 2*N);
    }
  }
  return (size/torus_base);
}

int generate_noiseless_trivial_LUT(TLweSample * result, const int * LUT, int size, int torus_base, const TLweParams * params){
  const int N = params->N;
  const int slice_size = N/torus_base;
  for (size_t i = 0; i < size/torus_base; i++){
    for (size_t j = 0; j < N; j++){
      result[i].a->coefsT[j] = 0;
      result[i].b->coefsT[j] = modSwitchToTorus32(LUT[i*torus_base + (j/slice_size)]*slice_size, 2*N);
    }
  }
  return (size/torus_base);
}

int generate_noiseless_trivial_LUT_with_scaling(TLweSample * result, const int * LUT, int size, int torus_base, int scaling, const TLweParams * params){
  const int N = params->N;
  const int slice_size = N/torus_base;
  for (size_t i = 0; i < size/torus_base; i++){
    for (size_t j = 0; j < N; j++){
      result[i].a->coefsT[j] = 0;
      result[i].b->coefsT[j] = modSwitchToTorus32(LUT[i*torus_base + (j/slice_size)]*slice_size, 2*N*scaling);
    }
  }
  return (size/torus_base);
}

// int generate_noiseless_trivial_LUT_with_scaling(TLweSample * result, int * LUT, int size, int torus_base, int scaling, const TLweParams * params){
//   const int N = params->N;
//   const int slice_size = N/torus_base;
//   for (size_t i = 0; i < size/torus_base; i++){
//     for (size_t j = 0; j < N; j++){
//       result[i].a->coefsT[j] = 0;
//       result[i].b->coefsT[j] = modSwitchToTorus32(LUT[i*torus_base + (j/slice_size)]*slice_size, 2*N*scaling);
//     }
//   }
//   return (size/torus_base);
// }


int generate_polynomial_LUT(IntPolynomial ** result, int * LUT, int size, int torus_base, int N){
  *result = new_IntPolynomial_array(size/torus_base, N);
  const int slice_size = N/torus_base;
  for (size_t i = 0; i < size/torus_base; i++){
    (*result)[i].coefs[0] = LUT[i*torus_base];
    for (size_t j = 1; j < N; j++){
      (*result)[i].coefs[j] = LUT[i*torus_base + (j/slice_size)];
    }
  }
  return (size/torus_base);
}

