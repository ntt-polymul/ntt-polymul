#include "api.h"
#include "hal.h"
#include "randombytes.h"
#include "sendfn.h"

#include <stdint.h>
#include <string.h>

#ifdef SABER_TYPE
#include "SABER_indcpa.h"
#include "SABER_params.h"
#include "poly.h"
void MatrixVectorMul(polyvec *a, uint16_t skpv[SABER_K][SABER_N], uint16_t res[SABER_K][SABER_N], uint16_t mod, int16_t transpose);
void InnerProd(uint16_t pkcl[SABER_K][SABER_N],uint16_t skpv[SABER_K][SABER_N],uint16_t mod,uint16_t res[SABER_N]);
#elif defined(NTRU_N)
#include "owcpa.h"
#include "poly.h"
#else
#include "bin-lwe.h"
#endif



#define printcycles(S, U) send_unsignedll((S), (U))

int main(void)
{
  unsigned char key_a[CRYPTO_BYTES], key_b[CRYPTO_BYTES];
  unsigned char sk[CRYPTO_SECRETKEYBYTES];
  unsigned char pk[CRYPTO_PUBLICKEYBYTES];
  unsigned char ct[CRYPTO_CIPHERTEXTBYTES];
  unsigned char cpa_r[32];

  unsigned long long t0, t1;
  int i;

  hal_setup(CLOCK_BENCHMARK);

  for(i=0;i<10;i++)
    hal_send_str("==========================");



  for(i=0;i<CRYPTO_ITERATIONS;i++) {
    memset(key_a, 0, sizeof key_a);
    memset(key_b, 0, sizeof key_a);
    memset(sk, 0, sizeof sk);
    memset(pk, 0, sizeof pk);
    memset(ct, 0, sizeof ct);

    // Key-pair generation
    t0 = hal_get_time();
    crypto_kem_keypair(pk, sk);
    t1 = hal_get_time();
    printcycles("cca keypair cycles:", t1-t0);

    // Encapsulation
    t0 = hal_get_time();
    crypto_kem_enc(ct, key_a, pk);
    t1 = hal_get_time();
    printcycles("encaps cycles:", t1-t0);

    // Decapsulation
    t0 = hal_get_time();
    crypto_kem_dec(key_b, ct, sk);
    t1 = hal_get_time();
    printcycles("decaps cycles:", t1-t0);

    if(memcmp(key_a, key_b, CRYPTO_BYTES)) {
      hal_send_str("ERROR KEYS\n");
    }
    else {
      hal_send_str("OK KEYS\n");
    }


    #ifdef SABER_TYPE

    t0 = hal_get_time();
    indcpa_kem_keypair(pk, sk);
    t1 = hal_get_time();
    printcycles("cpa keypair cycles:", t1-t0);

    randombytes(key_a, sizeof key_a);
    randombytes(cpa_r, sizeof cpa_r);

    t0 = hal_get_time();
    indcpa_kem_enc(key_a, cpa_r, pk, ct);
    t1 = hal_get_time();
    printcycles("cpa enc cycles:", t1-t0);

    t0 = hal_get_time();
    indcpa_kem_dec(sk, ct, key_b);
    t1 = hal_get_time();
    printcycles("cpa dec cycles:", t1-t0);

    if(memcmp(key_a, key_b, CRYPTO_BYTES)) {
      hal_send_str("ERROR KEYS\n");
    }
    else {
      hal_send_str("OK KEYS\n");
    }



    polyvec a;
    uint16_t skpv[SABER_K][SABER_N];
    uint16_t res[SABER_K][SABER_N];

    t0 = hal_get_time();
    MatrixVectorMul(&a, skpv, res, SABER_Q-1, 0);
    t1 = hal_get_time();
    printcycles("matrix vector mul cycles:", t1-t0);


    uint16_t pkcl[SABER_K][SABER_N];
    uint16_t skpv2[SABER_K][SABER_N];
    uint16_t res2[SABER_N];

    t0 = hal_get_time();
    InnerProd(pkcl, skpv2, SABER_Q-1, res2);
    t1 = hal_get_time();
    printcycles("inner prod cycles:", t1-t0);

    #elif defined(NTRU_N)

    unsigned char seed[NTRU_SEEDBYTES];
    randombytes(seed, sizeof seed);

    t0 = hal_get_time();
    owcpa_keypair(pk, sk, seed);
    t1 = hal_get_time();
    printcycles("cpa keypair cycles:", t1-t0);

    uint8_t rm1[NTRU_OWCPA_MSGBYTES];
    uint8_t rm_seed[NTRU_SAMPLE_RM_BYTES];

    randombytes(rm_seed, NTRU_SAMPLE_RM_BYTES);
    owcpa_samplemsg(rm1, rm_seed);

    t0 = hal_get_time();
    owcpa_enc(ct, rm1, pk);
    t1 = hal_get_time();
    printcycles("cpa enc cycles:", t1-t0);

    uint8_t rm2[NTRU_OWCPA_MSGBYTES];
    t0 = hal_get_time();
    owcpa_dec(rm2, ct, sk);
    t1 = hal_get_time();
    printcycles("cpa dec cycles:", t1-t0);

    if(memcmp(rm1, rm2, sizeof rm1)) {
      hal_send_str("ERROR KEYS\n");
    }
    else {
      hal_send_str("OK KEYS\n");
    }

    #ifndef TOOM
    poly r, a, b;
    randombytes((unsigned char *)&a, sizeof(poly));
    randombytes((unsigned char *)&b, sizeof(poly));
    t0 = hal_get_time();
    poly_SignedZ3_Rq_mul(&r, &a, &b);
    t1 = hal_get_time();
    printcycles("polymul cycles:", t1-t0);
    #endif


    #else // LAC

    unsigned long long clen, mlen;
    uint8_t buf[MESSAGE_LEN];


    t0 = hal_get_time();
    crypto_encrypt_keypair(pk, sk);
    t1 = hal_get_time();
    printcycles("cpa keypair cycles:", t1-t0);

    randombytes(key_a, MESSAGE_LEN);

    t0 = hal_get_time();
    crypto_encrypt(ct, &clen, key_a, MESSAGE_LEN, pk);
    t1 = hal_get_time();
    printcycles("cpa enc cycles:", t1-t0);

    t0 = hal_get_time();
    crypto_encrypt_open(key_b, &mlen, ct, clen, sk);
    t1 = hal_get_time();
    printcycles("cpa dec cycles:", t1-t0);

    if(memcmp(key_a, key_b, MESSAGE_LEN)) {
      hal_send_str("ERROR KEYS\n");
    }
    else {
      hal_send_str("OK KEYS\n");
    }

    uint8_t a[DIM_N], s[DIM_N], b[DIM_N];
    uint8_t seed[2*SAMPLE_LEN];
    randombytes(b, sizeof b);
    randombytes(seed, sizeof seed);
    gen_r(s, seed);

    t0 = hal_get_time();
    poly_mul(a, s, b, DIM_N);
    t1 = hal_get_time();
    printcycles("polymul cycles:", t1-t0);



    #endif



    hal_send_str("#");
  }
  while(1);
  return 0;
}
