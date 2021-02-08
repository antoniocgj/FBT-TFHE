#include <stdio.h>
#include <math.h>
#include <x86intrin.h>
#include <sys/time.h>
#include <assert.h>
#include "tfhe.h"
#include "tfhe_garbage_collector.h"

typedef struct TLweKeySwitchKey
{
  TLweSample *  ks0_raw;
  TLweSample ** ks1_raw;
  TLweSample *** ks2_raw;
  TLweSample **** ks;
  int32_t N, k, basebit, t;
  const TLweParams* out_params;
} TLweKeySwitchKey;

typedef struct fbt_integer{
  LweSample * lwe_samples;
  const LweParams * lwe_params;
  int log_torus_base;
  int digits;
} fbt_integer;

typedef struct fbt_lut{
  TLweSample * tlwe_samples = 0;
  IntPolynomial * int_polynomials = 0;
  const int * literal_promisse = 0;
  int size;
} fbt_lut;

uint64_t get_time() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (tv.tv_usec) + (tv.tv_sec * 1000000);
}

TLweKeySwitchKey * new_TLweKeySwitchKey(int32_t N, int32_t k, int32_t t, int32_t basebit, const TLweParams * out_params, const int torus_base);
void tLweCreateKeySwitchKey(TLweKeySwitchKey* result, const LweKey* in_key, const TLweKey* out_key, const int torus_base);
EXPORT void tLweNoiselessTrivialT(TLweSample *result, const Torus32 mu, const TLweParams *params);