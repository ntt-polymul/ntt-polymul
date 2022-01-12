#include <stdint.h>
#include <immintrin.h>
#include <assert.h>
#include "params.h"
#include "consts1024.h"
#include "nttrange.h"

void ntt_range(int16_t *coeffs, uint32_t *bounds) {
  int32_t j,k,l;
  int16_t zeta,t;
  char s[] = "Level 0";

  for(l=NTT_N;l>NTT_F;l/=2) {
    for(j=0;j<NTT_N;j+=l) {
      zeta = pdata[_ZETAS+j/l];
      for(k=j;k<j+l/2;k++) {
        if( (l == NTT_N/32 && (j < 64 || j >= 512))
            || (l == NTT_N/64 && j >= 64 && j < 512) ) // extra reduction levels 5+6
        {
          coeffs[k] = mulmod(coeffs[k],pdata[_16XMONT]);
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT]);
        }

        if(j > 0 || l <= NTT_N/32) {
          coeffs[l/2+k] = mulmod(coeffs[l/2+k],zeta);
          bounds[l/2+k] = maxmulmod(bounds[l/2+k],zeta);
        }
        t = coeffs[l/2+k];
        coeffs[l/2+k] = coeffs[k] - t;
        coeffs[k] = coeffs[k] + t;
        bounds[k] = bounds[l/2+k] = bounds[k] + bounds[l/2+k];
      }
    }
    print_bounds(bounds,s);
    s[6] += 1;
  }
}

void invntt_range(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,k,l;
  int16_t zeta,t;
  char s[] = "Inv Level 0";

  for(l=NTT_F;l<NTT_N;l*=2) {
    for(j=0;j<NTT_N;j+=2*l) {
      if(j > 0) {
        i = j/(2*l);
        i = (3 << (31 - _lzcnt_u32(i))) - i - 1;
        zeta = -pdata[_ZETAS+i];
      }
      else
        zeta = pdata[_16XMONT];
      for(k=j;k<j+l;k++) {
        t = coeffs[l+k];
        coeffs[l+k] = coeffs[k] - t;
        coeffs[k] = coeffs[k] + t;
        bounds[k] = bounds[l+k] = bounds[k] + bounds[l+k];
        if( (j > 0 || l <= 8*NTT_F)
            || (l == 64*NTT_F && k >= 32 && k < 96)
            || (l == 128*NTT_F && k >= 32) )
        {
          coeffs[l+k] = mulmod(coeffs[l+k],zeta);
          bounds[l+k] = maxmulmod(bounds[l+k],zeta);
        }

        if( (l == 2*NTT_F || l == 8*NTT_F)
            || (l == 32*NTT_F && k < j+32)
            || (l == 64*NTT_F && k >= j+32 && (j == 0 || k < j+64))
            || (l == 128*NTT_F && j == 0 && k >= 32)
            || (l == 128*NTT_F && j > 0 && (k < j+32 || (k >= j+64 && k < j+128))) )
        {
          coeffs[k] = mulmod(coeffs[k],pdata[_16XMONT]);
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT]);
        }
      }
    }
    print_bounds(bounds,s);
    s[10] += 1;
  }
}

void basemul_range(int16_t *coeffs0, const int16_t *coeffs1, uint32_t *bounds0, const uint32_t *bounds1) {
  int32_t i;
  int16_t zeta;
  int16_t t,u;
  uint32_t e,f;

  for(i=0;i<NTT_N/2;i++) {
    if(i%2 == 0) zeta = pdata[_ZETAS+i/2];
    else zeta = -pdata[_ZETAS+i/2];

    coeffs0[2*i+0] = mulmod(coeffs0[2*i+0],pdata[_16XMONT]);
    coeffs0[2*i+1] = mulmod(coeffs0[2*i+1],pdata[_16XMONT]);
    t  = mulmod(coeffs0[2*i+1],coeffs1[2*i+1]);
    t  = mulmod(t,zeta);
    t += mulmod(coeffs0[2*i+0],coeffs1[2*i+0]);
    u  = mulmod(coeffs0[2*i+0],coeffs1[2*i+1]);
    u += mulmod(coeffs0[2*i+1],coeffs1[2*i+0]);
    t  = mulmod(t,pdata[_16XF]);
    u  = mulmod(u,pdata[_16XF]);
    coeffs0[2*i+0] = t;
    coeffs0[2*i+1] = u;

    bounds0[2*i+0] = maxmulmod(bounds0[2*i+0],pdata[_16XMONT]);
    bounds0[2*i+1] = maxmulmod(bounds0[2*i+1],pdata[_16XMONT]);
    e  = maxmulmod(bounds0[2*i+1],bounds1[2*i+1]);
    e  = maxmulmod(e,zeta);
    e += maxmulmod(bounds0[2*i+0],bounds1[2*i+0]);
    f  = maxmulmod(bounds0[2*i+0],bounds1[2*i+1]);
    f += maxmulmod(bounds0[2*i+1],bounds1[2*i+0]);
    e  = maxmulmod(e,pdata[_16XF]);
    f  = maxmulmod(f,pdata[_16XF]);
    bounds0[2*i+0] = e;
    bounds0[2*i+1] = f;
  }
  print_bounds(bounds0,"Basemul");
}

void crt_range(int16_t *coeffs0, const int16_t *coeffs1, uint32_t *bounds0, const uint32_t *bounds1) {
  int32_t i;
  int16_t t;
  uint32_t u;

  pdata = PDATA1;
  for(i=0;i<2*KEM_N;i++) {
    t = mulmod(coeffs1[i] - coeffs0[i],CRT_U);
    t *= P0;
    coeffs0[i] = coeffs0[i] + t;

    u = maxmulmod(bounds1[i] + (P0-1)/2,CRT_U);  // extra reduction
    u *= P0;
    bounds0[i] = bounds0[i] + u;
    assert(bounds0[i] < P0*P1 - KEM_Q/2*KEM_N);
  }
}
