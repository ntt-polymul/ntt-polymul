#include <stdio.h>
#include "params.h"
#include "randombytes.h"
#include "poly.h"
#include "consts.h"

#define P P1
#define PDATA PDATA1

static void poly_naivemul(poly *c, const poly *a, const poly *b) {
  unsigned int i,j;
  int16_t buf[2*KEM_N] = {0};

  for(i = 0; i < KEM_N; ++i)
    for(j = 0; j < KEM_N; ++j)
      buf[i+j] = (buf[i+j] + (int32_t)a->coeffs[i]*b->coeffs[j]) % P;

  for(i=0;i<KEM_N;i++)
    c->coeffs[i] = (buf[i] - buf[KEM_N+i]) % P;
}

int main(void) {
  int j, err;
  uint16_t nonce = 0;
  uint8_t seed[POLYMUL_SYMBYTES];
  poly a, b, c;
  nttpoly ahat, bhat;

  randombytes(seed,POLYMUL_SYMBYTES);
  poly_uniform(&a,seed,nonce++);
  poly_uniform(&b,seed,nonce++);

  poly_naivemul(&c,&a,&b);

  poly_ntt(&ahat,&a,PDATA);
  poly_ntt(&bhat,&b,PDATA);
  poly_basemul_montgomery(&ahat,&ahat,&bhat,PDATA);
  poly_invntt_tomont(&ahat,&ahat,PDATA);

  err = 0;
  for(j=0;j<KEM_N;j++)
    if((c.coeffs[j] - ahat.coeffs[j]) % P) {
      printf("ERROR: %d, %d, %d\n", j, c.coeffs[j], ahat.coeffs[j]);
      err = 1;
    }
  if(!err) printf("ALL GOOD.\n");
  return 0;
}
