#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include "randombytes.h"
#include "params.h"
#include "poly.h"
#include "cpucycles.h"
#include "consts.h"

#define P P1
#define PDATA PDATA1

static int16_t pow_simple(int16_t a, unsigned int e) {
  int16_t r;
  if(e == 0) return 1;
  else if(e == 1) return a;

  r = pow_simple(a,e/2);
  r = (int32_t)r*r % P;
  if(e&1) r = (int32_t)r*a % P;

  if(r > (P-1)/2) r -= P;
  else if(r < -(P-1)/2) r += P;

  return r;
}

#if NTT_N == 256
static int idx(int i) {
  int r;
  r  = i/128*128;
  i %= 128;
  r += i/8;
  i %= 8;
  r += 16*i;
  return r;
}
#elif NTT_N == 512
static int idx(int i) {
  int r;
  r  = i/32*32;
  i %= 32;
  r += i/4;
  i %= 4;
  r += i/2*16;
  i %= 2;
  r += i*8;
  return r;
}
#elif NTT_N == 1024
static int idx(int i) {
  int r;
  r  = i/128*128;
  i %= 128;
  r += i/8;
  i %= 8;
  r += i/4*64;
  i %= 4;
  r += i*16;
  return r;
}
#endif

static int16_t zeta[NTT_N/NTT_F];
static int16_t zetapow[NTT_N/NTT_F][POLY_N/NTT_F];

int main(void) {
  int i,j, err;
  //uint64_t t[20], overhead;
  int16_t out[NTT_N/NTT_F];
  uint8_t seed[POLYMUL_SYMBYTES];
  poly a;
  nttpoly b;
  err = 0;

  //overhead = cpucycles_overhead();
  randombytes(seed,POLYMUL_SYMBYTES);

  for(i=0;i<POLY_N;i++)
    a.coeffs[i] = 0;
  a.coeffs[NTT_F] = 1;
  poly_ntt(&b,&a,PDATA);
  for(i=0;i<NTT_N/NTT_F;i++) {
    zeta[i] = b.coeffs[idx(NTT_F*i)] % P;
    if((pow_simple(zeta[i],NTT_N/NTT_F) + 1) % P) {
      fprintf(stderr, "ERROR1: %d, %d\n", NTT_F*i, zeta[i]);
      err = 1;
    }
    for(j=0;j<i;j++)
      if((zeta[j] - zeta[i]) % P == 0) {
        fprintf(stderr, "ERROR2: %d, %d, %d, %d\n", NTT_F*i, NTT_F*j, zeta[i], zeta[j]);
        err = 1;
      }
  }

  for(i=0;i<NTT_N/NTT_F;i++)
    for(j=0;j<POLY_N/NTT_F;j++)
      zetapow[i][j] = pow_simple(zeta[i],j);

  for(i=0;i<POLY_N/NTT_F;i++) {
    for(j=0;j<POLY_N;j++)
      a.coeffs[j] = 0;
    a.coeffs[NTT_F*i] = 1;
    poly_ntt(&b,&a,PDATA);
    for(j=0;j<NTT_N/NTT_F;j++)
      if((b.coeffs[idx(NTT_F*j)] - zetapow[j][i]) % P) {
        fprintf(stderr,"ERROR3: %d, %d, %d %d\n", i, j, b.coeffs[idx(NTT_F*j)], zetapow[j][i]);
        err = 1;
      }
  }

  poly_uniform(&a,seed,0);
  for(i=0;i<NTT_N/NTT_F;i++) {
    out[i] = 0;
    for(j=0;j<POLY_N/NTT_F;j++)
      out[i] = (out[i] + (int32_t)a.coeffs[NTT_F*j]*zetapow[i][j]) % P;
  }
  poly_ntt(&b,&a,PDATA);
  for(i=0;i<NTT_N/NTT_F;i++)
    if((b.coeffs[idx(NTT_F*i)] - out[i]) % P) {
      fprintf(stderr,"ERROR4: %d, %d %d\n", i, b.coeffs[idx(NTT_F*i)], out[i]);
      err = 1;
    }

  poly_uniform(&a,seed,1);
  poly_ntt(&b,&a,PDATA);
  const int16_t div = pow_simple(NTT_N/NTT_F,P-2);
  for(i=0;i<NTT_N;i++)
    b.coeffs[i] = (int32_t)b.coeffs[i]*div % P;
  for(i=0;i<NTT_N;i++)
    if(b.coeffs[i] > (P-1)/2) b.coeffs[i] -= P;
  for(i=0;i<NTT_N;i++)
    if(b.coeffs[i] < -(P-1)/2) b.coeffs[i] += P;
  poly_invntt_tomont(&b,&b,PDATA);
  for(i=0;i<KEM_N;++i)
    if((a.coeffs[i] - b.coeffs[i]) % P) {
      fprintf(stderr, "ERROR5: %d, %d, %d\n", i, a.coeffs[i], b.coeffs[i]);
      err = 1;
    }
  for(i=KEM_N;i<NTT_N;++i)
    if(b.coeffs[i] % P) {
      fprintf(stderr, "ERROR6: %d, %d \n", i, b.coeffs[i]);
      err = 1;
    }

  if(!err) printf("ALL GOOD.\n");
  return 0;

/*
  for(i=0;i<20;i++) {
    t[i] = cpucycles();
    poly_ntt(&b,&a,PDATA);
  }
  for(i=0;i<19;i++)
    printf("ntt: %2d: %lu\n", i+1, t[i+1] - t[i] - overhead);
  for(i=0;i<20;i++) {
    t[i] = cpucycles();
    poly_invntt_tomont(&b,&b,PDATA);
  }
  for(i=0;i<19;i++)
    printf("invntt: %2d: %lu\n", i, t[i+1] - t[i] - overhead);
*/
}
