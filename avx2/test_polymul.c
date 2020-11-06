#include <stdio.h>
#include "params.h"
#include "randombytes.h"
#include "poly.h"
#include "consts.h"

static void poly_naivemul(poly *c, const poly *a, const poly *b) {
  unsigned int i,j;
  int16_t t[2*KEM_N] = {0};

  for(i = 0; i < KEM_N; ++i)
    for(j = 0; j < KEM_N; ++j)
      t[i+j] = (t[i+j] + (int32_t)a->coeffs[i]*b->coeffs[j]) % KEM_Q;

  for(i = KEM_N; i < 2*KEM_N; i++)
    c->coeffs[i - KEM_N] = (t[i - KEM_N] + t[i]) % KEM_Q;
}

int main(void) {
  int i;
  uint16_t nonce = 0;
  uint8_t seed[POLYMUL_SYMBYTES];
  poly a, b, c, d;
  nttpoly ahat, bhat, chat;

  randombytes(seed,POLYMUL_SYMBYTES);
  poly_uniform(&a,seed,nonce++);
  poly_noise(&b,seed,nonce++);

  poly_naivemul(&c,&a,&b);
  ntru_poly_mul(&d,&a,&b);
  for(i=0;i<KEM_N;i++)
    if((c.coeffs[i] - d.coeffs[i]) % KEM_Q)
      printf("ERROR1: %d\n", i);

  poly_ntt(&ahat,&a,PDATA0);
  poly_ntt(&bhat,&b,PDATA0);
  poly_basemul_montgomery(&ahat,&ahat,&bhat,PDATA0);
  poly_invntt_tomont(&ahat,&ahat,PDATA0);
  poly_ntt(&bhat,&a,PDATA1);
  poly_ntt(&chat,&b,PDATA1);
  poly_basemul_montgomery(&bhat,&bhat,&chat,PDATA1);
  poly_invntt_tomont(&bhat,&bhat,PDATA1);
  poly_crt(&d,&ahat,&bhat);
  for(i=0;i<KEM_N;i++)
    if((c.coeffs[i] - d.coeffs[i]) % KEM_Q)
      printf("ERROR2: %d, %d, %d\n", i, c.coeffs[i], d.coeffs[i]);

  return 0;
}
