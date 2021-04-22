#include <stdio.h>
#include "params.h"
#include "randombytes.h"
#include "polyvec.h"
#include "poly.h"
#include "consts.h"

static void poly_naivemul(poly *c, const poly *a, const poly *b) {
  unsigned int i,j;
  int16_t t[2*KEM_N] = {0};

  for(i = 0; i < KEM_N; ++i)
    for(j = 0; j < KEM_N; ++j)
      t[i+j] = (t[i+j] + (int32_t)a->coeffs[i]*b->coeffs[j]) % KEM_Q;

  for(i = KEM_N; i < 2*KEM_N; i++)
    c->coeffs[i - KEM_N] = (t[i - KEM_N] - t[i]) % KEM_Q;
}

int main(void) {
  int i,j,err;
  uint16_t nonce = 0;
  uint8_t seed[POLYMUL_SYMBYTES];
  polyvec a[KEM_K];
  polyvec s, t, u;
  poly tmp;
  err = 0;

  randombytes(seed,POLYMUL_SYMBYTES);
  for(i=0;i<KEM_K;i++)
    polyvec_uniform(&a[i],seed,nonce++);
  polyvec_noise(&s,seed,nonce++);

  for(i=0;i<KEM_K;i++) {
    poly_naivemul(&t.vec[i],&a[i].vec[0],&s.vec[0]);
    for(j=1;j<KEM_K;j++) {
      poly_naivemul(&tmp,&a[i].vec[j],&s.vec[j]);
      poly_add(&t.vec[i],&t.vec[i],&tmp);
    }
  }

  saber_matrix_vector_mul(&u,a,&s);
  for(i=0;i<KEM_K;i++)
    for(j=0;j<KEM_N;j++)
      if((u.vec[i].coeffs[j] - t.vec[i].coeffs[j]) % KEM_Q) {
        printf("ERROR1: %d, %d, %d, %d\n", i, j, t.vec[i].coeffs[j], u.vec[i].coeffs[j]);
        err = 1;
      }

  polyvec_matrix_vector_mul(&u,a,&s,0);
  for(i=0;i<KEM_K;i++)
    for(j=0;j<KEM_N;j++)
      if((u.vec[i].coeffs[j] - t.vec[i].coeffs[j]) % KEM_Q) {
        printf("ERROR2: %d, %d, %d, %d\n", i, j, t.vec[i].coeffs[j], u.vec[i].coeffs[j]);
        err = 1;
      }
  if(!err) printf("ALL GOOD.\n");

  return 0;
}
