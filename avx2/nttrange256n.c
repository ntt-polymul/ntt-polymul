#include <stdint.h>
#include <assert.h>
#include "params.h"
#include "consts256.h"
#include "nttrange.h"

static void ntt0t2(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,k,l;
  int16_t t,zeta;
  char s[] = "Level 0";

  i = 0;
  for(l=NTT_N;l>NTT_N/8;l/=2) {
    for(j=0;j<NTT_N;j+=l) {
      const int32_t idx[7] = {0,32,96,64,72,128,136};
      zeta = pdata[_ZETAS+16+idx[i++]];
      for(k=j;k<j+l/2;k++) {
        t = mulmod(coeffs[l/2+k],zeta);
        coeffs[l/2+k] = coeffs[k] - t;
        coeffs[k] = coeffs[k] + t;

        bounds[l/2+k] = maxmulmod(bounds[l/2+k],zeta);
        bounds[k] = bounds[l/2+k] = bounds[k] + bounds[l/2+k];
      }
    }
    print_bounds(bounds,s);
    s[6] += 1;
  }
}

static void twist32(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,idx;
  int16_t zeta;

  for(i=0;i<NTT_N;i++) {
    j  = i;
    idx  = j/128*256;
    j %= 128;
    idx += j/32*4;
    j %= 32;
    idx += j/4*32;
    j %= 4;
    idx += j;
    zeta = pdata[_TWIST32+16+idx];
    coeffs[i] = mulmod(coeffs[i],zeta);
    bounds[i] = maxmulmod(bounds[i],zeta);
  }
  print_bounds(bounds,"Twist 32");
}

static void ntt3t5(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,k,l;
  int16_t t,zeta;
  char s[] = "Level 3";

  for(l=NTT_N/8;l>NTT_N/64;l/=2) {
    for(j=0;j<NTT_N;j+=l) {
      i = j%(NTT_N/8) / l;
      const int32_t idx[4] = {-1,0,32,96};
      if(i > 0) zeta = pdata[_ZETAS+16+idx[i]];
      else zeta = pdata[_16XMONT];
      for(k=j;k<j+l/2;k++) {
        if(i > 0 || l <= NTT_N/16) t = mulmod(coeffs[l/2+k],zeta);  // extra reduction
        else t = coeffs[l/2+k];
        coeffs[l/2+k] = coeffs[k] - t;
        coeffs[k] = coeffs[k] + t;

        if(i > 0 || l <= NTT_N/16) bounds[l/2+k] = maxmulmod(bounds[l/2+k],zeta);
        bounds[k] = bounds[l/2+k] = bounds[k] + bounds[l/2+k];
      }
    }
    print_bounds(bounds,s);
    s[6] += 1;
  }
}

static void twist4(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,idx;
  int16_t zeta;

  for(i=0;i<NTT_N;i++) {
    j  = i % 32;
    idx = j/4*8;
    j %= 4;
    idx += j;
    if(idx < 4) zeta = pdata[_16XMONT];
    else zeta = pdata[_TWIST4+4+idx];
    coeffs[i] = mulmod(coeffs[i],zeta);
    bounds[i] = maxmulmod(bounds[i],zeta);
  }
  print_bounds(bounds,"Twist 4");
}

static void ntt6t7(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,k,l;
  int16_t t,zeta;
  char s[] = "Level 6";

  zeta = pdata[_ZETAS+16];
  for(l=NTT_N/64;l>NTT_N/256;l/=2) {
    for(j=0;j<NTT_N;j+=l) {
      i = j%(NTT_N/64) / l;
      for(k=j;k<j+l/2;k++) {
        if(i > 0) t = mulmod(coeffs[l/2+k],zeta);
        else t = coeffs[l/2+k];
        coeffs[l/2+k] = coeffs[k] - t;
        coeffs[k] = coeffs[k] + t;

        if(i > 0) bounds[l/2+k] = maxmulmod(bounds[l/2+k],zeta);
        bounds[k] = bounds[l/2+k] = bounds[k] + bounds[l/2+k];
      }
    }
    print_bounds(bounds,s);
    s[6] += 1;
  }
}

void ntt_range(int16_t *coeffs, uint32_t *bounds) {
  ntt0t2(coeffs,bounds);
  twist32(coeffs,bounds);
  ntt3t5(coeffs,bounds);
  twist4(coeffs,bounds);
  ntt6t7(coeffs,bounds);
}

static void invntt0t1(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,k,l;
  int16_t zeta,t;
  char s[] = "Inv Level 0";

  zeta = -pdata[_ZETAS+16];
  for(l=NTT_N/128;l<=NTT_N/64;l*=2) {
    for(j=0;j<NTT_N;j+=l) {
      i = j%(NTT_N/64) / l;
      for(k=j;k<j+l/2;k++) {
        t = coeffs[l/2+k];
        coeffs[l/2+k] = coeffs[k] - t;
        coeffs[k] = coeffs[k] + t;
        if(i > 0) coeffs[l/2+k] = mulmod(coeffs[l/2+k],zeta);

        bounds[k] = bounds[l/2+k] = bounds[k] + bounds[l/2+k];
        if(i > 0) bounds[l/2+k] = maxmulmod(bounds[l/2+k],zeta);
      }
    }
    print_bounds(bounds,s);
    s[10] += 1;
  }
}

static void invtwist4(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,idx;
  int16_t zeta;

  for(i=0;i<NTT_N;i++) {
    const int32_t lut[8] = {-1, 0, 3, 2, 7, 6, 5, 4};
    j  = i % 32;
    idx = lut[j/4]*8;
    j %= 4;
    idx += j;
    if(idx < 0) zeta = pdata[_16XMONT];
    else zeta = pdata[_TWIST4+4+idx];
    coeffs[i] = mulmod(coeffs[i],zeta);
    bounds[i] = maxmulmod(bounds[i],zeta);
  }
  print_bounds(bounds,"Inv Twist 4");
}

static void invntt2t4(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,k,l;
  int16_t zeta,t;
  char s[] = "Inv Level 2";

  for(l=NTT_N/32;l<=NTT_N/8;l*=2) {
    for(j=0;j<NTT_N;j+=l) {
      i = j%(NTT_N/8) / l;
      const int32_t lut[4] = {-1,0,96,32};
      zeta = -pdata[_ZETAS+16+lut[i]];
      for(k=j;k<j+l/2;k++) {
        if(l == NTT_N/8 && (k%(NTT_N/8)) < NTT_N/64) {  // extra reduction
          coeffs[k] = mulmod(coeffs[k],pdata[_16XMONT]);
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT]);
        }

        t = coeffs[l/2+k];
        coeffs[l/2+k] = coeffs[k] - t;
        coeffs[k] = coeffs[k] + t;
        if(i > 0) coeffs[l/2+k] = mulmod(coeffs[l/2+k],zeta);

        bounds[k] = bounds[l/2+k] = bounds[k] + bounds[l/2+k];
        if(i > 0) bounds[l/2+k] = maxmulmod(bounds[l/2+k],zeta);
      }
    }
    print_bounds(bounds,s);
    s[10] += 1;
  }
}

static void invtwist32(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,idx;
  int16_t zeta;

  for(i=0;i<NTT_N;i++) {
    j  = i;
    idx  = (1-j/128)*256;
    j %= 128;
    idx += (3-j/32)*4;
    j %= 32;
    idx += j/4*32;
    j %= 4;
    idx += j;
    zeta = pdata[_TWIST32+16+idx];
    coeffs[i] = mulmod(coeffs[i],zeta);
    bounds[i] = maxmulmod(bounds[i],zeta);
  }
  print_bounds(bounds,"Inv Twist 32");
}

static void invntt5t7(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,k,l;
  int16_t zeta,t;
  char s[] = "Inv Level 5";

  i = 6;
  for(l=NTT_N/4;l<=NTT_N;l*=2) {
    for(j=0;j<NTT_N;j+=l) {
      const int32_t idx[7] = {0,32,96,64,72,128,136};
      zeta = -pdata[_ZETAS+16+idx[i--]];
      for(k=j;k<j+l/2;k++) {
        t = coeffs[l/2+k];
        coeffs[l/2+k] = coeffs[k] - t;
        coeffs[k] = coeffs[k] + t;
        coeffs[l/2+k] = mulmod(coeffs[l/2+k],zeta);

        bounds[k] = bounds[l/2+k] = bounds[k] + bounds[l/2+k];
        bounds[l/2+k] = maxmulmod(bounds[l/2+k],zeta);
      }
      if(l == NTT_N/2) {  // extra reduction
        for(k=j;k<j+l/4;k++) {
          coeffs[k] = mulmod(coeffs[k],pdata[_16XMONT]);
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT]);
        }
      }
    }
    print_bounds(bounds,s);
    s[10] += 1;
  }
}

void invntt_range(int16_t *coeffs, uint32_t *bounds) {
  invntt0t1(coeffs,bounds);
  invtwist4(coeffs,bounds);
  invntt2t4(coeffs,bounds);
  invtwist32(coeffs,bounds);
  invntt5t7(coeffs,bounds);
}

void basemul_range(int16_t *coeffs0, const int16_t *coeffs1, uint32_t *bounds0, const uint32_t *bounds1) {
  int32_t i;

  for(i=0;i<NTT_N;i++) {
    coeffs0[i] = mulmod(coeffs0[i],coeffs1[i]);
    coeffs0[i] = mulmod(coeffs0[i],pdata[_16XF]);

    bounds0[i] = 2*maxmulmod2(bounds0[i],bounds1[i]);  // accumulation
    bounds0[i] = maxmulmod(bounds0[i],pdata[_16XF]);
  }

  print_bounds(bounds0,"Basemul");
}

void crt_range(int16_t *coeffs0, const int16_t *coeffs1, uint32_t *bounds0, const uint32_t *bounds1) {
  int32_t i;
  int16_t t;
  uint32_t u;

  for(i=0;i<NTT_N;i++) {
    pdata = PDATA0;
    coeffs0[i] = mulmod(coeffs0[i],PDATA0[_16XMONT]);
    bounds0[i] = maxmulmod(bounds0[i],PDATA0[_16XMONT]);
    pdata = PDATA1;
    t = mulmod(coeffs1[i] - coeffs0[i],CRT_U);
    u = maxmulmod(bounds1[i] + bounds0[i],CRT_U);
    t *= P0;
    u *= P0;
    coeffs0[i] = coeffs0[i] + t;
    bounds0[i] = bounds0[i] + u;
    assert(bounds0[i] < P0*P1-5*4096*256*4);
  }
}
