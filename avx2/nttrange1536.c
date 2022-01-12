#include <stdint.h>
#include <assert.h>
#include "params.h"
#include "consts1536.h"
#include "nttrange.h"

static void ntt0t3(int16_t *coeffs, uint32_t *bounds) {
  int32_t j,k,l;
  int16_t zeta,t;
  char s[] = "Level 0";

  for(l=NTT_N;l>=NTT_N/8;l/=2) {
    for(j=0;j<NTT_N;j+=l) {
      zeta = pdata[_ZETAS+2*(j/l-1)];
      for(k=j;k<j+l/2;k++) {
        if(j > 0) {
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

static void twist96(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,idx;
  int16_t zeta;

  for(i=0;i<NTT_N;i++) {
    j  = i;
    idx  = j/16*32;
    j %= 16;
    idx += j;
    if(idx < 192) zeta = pdata[_16XMONT];
    else zeta = pdata[_TWIST96+16+idx];
    coeffs[i] = mulmod(coeffs[i],zeta);
    bounds[i] = maxmulmod(bounds[i],zeta);
  }
  print_bounds(bounds,"Twist 96");
}

static void ntt4t6(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,k,l;
  int16_t zeta,t;
  char s[] = "Level 4";

  for(l=NTT_N/16;l>=NTT_N/64;l/=2) {
    for(j=0;j<NTT_N;j+=l) {
      i = j%(NTT_N/16) / l;
      if(i > 0) zeta = pdata[_ZETAS+2*(i-1)];
      else zeta = pdata[_16XMONT];
      for(k=j;k<j+l/2;k++) {
        if(l == NTT_N/64 && i <= 1 && k >= j+8) {
          coeffs[k] = mulmod(coeffs[k],pdata[_16XMONT]);
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT]);
        }

        if(i > 0 || (l == NTT_N/32 && k < j+8) || l == NTT_N/64) {
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

static void twist12(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,idx;
  int16_t zeta;

  for(i=0;i<NTT_N;i++) {
    j  = i % 96;
    idx  = j/48*96;
    j %= 48;
    idx += j/24*4;
    j %= 24;
    idx += j/4*16;
    j %= 4;
    idx += j;
    if(idx < 48 && idx % 8 < 4) zeta = pdata[_16XMONT];
    else zeta = pdata[_TWIST12+8+idx];
    coeffs[i] = mulmod(coeffs[i],zeta);
    bounds[i] = maxmulmod(bounds[i],zeta);
  }
  print_bounds(bounds,"Twist 12");
}

static void ntt7t8(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,k,l;
  int16_t zeta,t;
  char s[] = "Level 7";

  for(l=NTT_N/128;l>NTT_F;l/=2) {
    for(j=0;j<NTT_N;j+=l) {
      i = j%(NTT_N/128) / l;
      zeta = pdata[_ZETAS+2*(i-1)];
      for(k=j;k<j+l/2;k++) {
        if(i > 0) {
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

void ntt_range(int16_t *coeffs, uint32_t *bounds) {
  ntt0t3(coeffs,bounds);
  twist96(coeffs,bounds);
  ntt4t6(coeffs,bounds);
  twist12(coeffs,bounds);
  ntt7t8(coeffs,bounds);
}

static void invntt0t1(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,k,l;
  int16_t zeta,t;
  char s[] = "Inv Level 0";

  for(l=NTT_F;l<=2*NTT_F;l*=2) {
    for(j=0;j<NTT_N;j+=2*l) {
      const int32_t lut[4] = {-1, 0, 4, 2};
      i = j%(4*NTT_F) / (2*l);
      zeta = -pdata[_ZETAS+lut[i]];
      for(k=j;k<j+l;k++) {
        t = coeffs[l+k];
        coeffs[l+k] = coeffs[k] - t;
        coeffs[k] = coeffs[k] + t;
        bounds[k] = bounds[l+k] = bounds[k] + bounds[l+k];
        if(i > 0) {
          coeffs[l+k] = mulmod(coeffs[l+k],zeta);
          bounds[l+k] = maxmulmod(bounds[l+k],zeta);
        }
      }
    }
    print_bounds(bounds,s);
    s[10] += 1;
  }
}

static void invtwist12(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,idx;
  int16_t zeta;

  for(i=0;i<NTT_N;i++) {
    j  = i % 96;
    const int32_t lut[] = {-48,0,52,4,148,100,144,96};
    idx  = lut[j/12];
    j %= 12;
    idx += j/4*16;
    j %= 4;
    idx += j;
    if(idx < 0) zeta = pdata[_16XMONT];
    else zeta = pdata[_TWIST12+8+idx];
    coeffs[i] = mulmod(coeffs[i],zeta);
    bounds[i] = maxmulmod(bounds[i],zeta);
  }
  print_bounds(bounds,"Inv Twist 12");
}

static void invntt2t4(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,k,l;
  int16_t zeta,t;
  char s[] = "Inv Level 2";

  for(l=4*NTT_F;l<=16*NTT_F;l*=2) {
    for(j=0;j<NTT_N;j+=2*l) {
      const int32_t lut[4] = {-1, 0, 4, 2};
      i = j%(32*NTT_F) / (2*l);
      zeta = -pdata[_ZETAS+lut[i]];
      for(k=j;k<j+l;k++) {
        t = coeffs[l+k];
        coeffs[l+k] = coeffs[k] - t;
        coeffs[k] = coeffs[k] + t;
        bounds[k] = bounds[l+k] = bounds[k] + bounds[l+k];
        if(i > 0) {
          coeffs[l+k] = mulmod(coeffs[l+k],zeta);
          bounds[l+k] = maxmulmod(bounds[l+k],zeta);
        }

        if( (l == 4*NTT_F && i == 1)
            || (l == 8*NTT_F && i == 1 && k < j+16) )
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

static void invtwist96(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,idx;
  int16_t zeta;

  for(i=0;i<NTT_N;i++) {
    const int32_t lut[16] = {-1,0,3,2,7,6,5,4,15,14,13,12,11,10,9,8};
    j  = i;
    idx  = lut[j/96]*192;
    j %= 96;
    idx += j/16*32;
    j %= 16;
    idx += j;
    if(idx < 0) zeta = pdata[_16XMONT];
    else zeta = pdata[_TWIST96+16+idx];
    coeffs[i] = mulmod(coeffs[i],zeta);
    bounds[i] = maxmulmod(bounds[i],zeta);
  }
  print_bounds(bounds,"Inv Twist 96");
}

static void invntt5t8(int16_t *coeffs, uint32_t *bounds) {
  int32_t j,k,l;
  int16_t zeta,t;
  char s[] = "Inv Level 5";

  for(l=32*NTT_F;l<=256*NTT_F;l*=2) {
    for(j=0;j<NTT_N;j+=2*l) {
      const int32_t lut[8] = {-1, 0, 4, 2, 12, 10, 8, 6};
      if(j > 0) zeta = -pdata[_ZETAS+lut[j/(2*l)]];
      else zeta = pdata[_16XMONT];
      for(k=j;k<j+l;k++) {
        t = coeffs[l+k];
        coeffs[l+k] = coeffs[k] - t;
        coeffs[k] = coeffs[k] + t;
        bounds[k] = bounds[l+k] = bounds[k] + bounds[l+k];
        if(j > 0 || l == 64*NTT_F) {
          coeffs[l+k] = mulmod(coeffs[l+k],zeta);
          bounds[l+k] = maxmulmod(bounds[l+k],zeta);
        }

        if( (l == 64*NTT_F && (j == 0 || k < j+96))
            || (l == 128*NTT_F && k >= 768+96 && k < 768+192) )
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

void invntt_range(int16_t *coeffs, uint32_t *bounds) {
  invntt0t1(coeffs,bounds);
  invtwist12(coeffs,bounds);
  invntt2t4(coeffs,bounds);
  invtwist96(coeffs,bounds);
  invntt5t8(coeffs,bounds);
}

void basemul_range(int16_t *coeffs0, const int16_t *coeffs1, uint32_t *bounds0, const uint32_t *bounds1) {
  int32_t i;
  int16_t zeta;
  int16_t t[5];
  int32_t u[5];

  for(i=0;i<NTT_N/3;i++) {
    zeta = pdata[_ZETAS+2*(i/2-1)];
    zeta = (i%2==0) ? zeta : -zeta;
    coeffs0[3*i+0] = maxmulmod(coeffs0[3*i+0],pdata[_16XMONT]);  // FIXME
    coeffs0[3*i+1] = maxmulmod(coeffs0[3*i+1],pdata[_16XMONT]);
    coeffs0[3*i+2] = maxmulmod(coeffs0[3*i+2],pdata[_16XMONT]);
    t[0]  = mulmod(coeffs0[3*i+0],coeffs1[3*i+0]);
    t[1]  = mulmod(coeffs0[3*i+0],coeffs1[3*i+1]);
    t[2]  = mulmod(coeffs0[3*i+0],coeffs1[3*i+2]);
    t[1] += mulmod(coeffs0[3*i+1],coeffs1[3*i+0]);
    t[2] += mulmod(coeffs0[3*i+1],coeffs1[3*i+1]);
    t[3]  = mulmod(coeffs0[3*i+1],coeffs1[3*i+2]);
    t[2] += mulmod(coeffs0[3*i+2],coeffs1[3*i+0]);
    t[3] += mulmod(coeffs0[3*i+2],coeffs1[3*i+1]);
    t[4]  = mulmod(coeffs0[3*i+2],coeffs1[3*i+2]);
    if(i == 0) {
      t[0] += t[3];
      t[1] += t[4];
    }
    else if(i == 1) {
      t[0] -= t[3];
      t[1] -= t[4];
    }
    else {
      t[0] += mulmod(t[3],zeta);
      t[1] += mulmod(t[4],zeta);
    }
    t[0] = mulmod(t[0],pdata[_16XF]);
    t[1] = mulmod(t[1],pdata[_16XF]);
    t[2] = mulmod(t[2],pdata[_16XF]);
    coeffs0[3*i+0] = t[0];
    coeffs0[3*i+1] = t[1];
    coeffs0[3*i+2] = t[2];

    bounds0[3*i+0] = maxmulmod(bounds0[3*i+0],pdata[_16XMONT]);
    bounds0[3*i+1] = maxmulmod(bounds0[3*i+1],pdata[_16XMONT]);
    bounds0[3*i+2] = maxmulmod(bounds0[3*i+2],pdata[_16XMONT]);
    u[0]  = maxmulmod2(bounds0[3*i+0],bounds1[3*i+0]);
    u[1]  = maxmulmod2(bounds0[3*i+0],bounds1[3*i+1]);
    u[2]  = maxmulmod2(bounds0[3*i+0],bounds1[3*i+2]);
    u[1] += maxmulmod2(bounds0[3*i+1],bounds1[3*i+0]);
    u[2] += maxmulmod2(bounds0[3*i+1],bounds1[3*i+1]);
    u[3]  = maxmulmod2(bounds0[3*i+1],bounds1[3*i+2]);
    u[2] += maxmulmod2(bounds0[3*i+2],bounds1[3*i+0]);
    u[3] += maxmulmod2(bounds0[3*i+2],bounds1[3*i+1]);
    u[4]  = maxmulmod2(bounds0[3*i+2],bounds1[3*i+2]);
    if(i == 0) {
      u[0] += u[3];
      u[1] += u[4];
    }
    else if(i == 1) {
      u[0] += u[3];
      u[1] += u[4];
    }
    else {
      u[0] += maxmulmod(u[3],zeta);
      u[1] += maxmulmod(u[4],zeta);
    }
    u[0] = maxmulmod(u[0],pdata[_16XF]);
    u[1] = maxmulmod(u[1],pdata[_16XF]);
    u[2] = maxmulmod(u[2],pdata[_16XF]);
    bounds0[3*i+0] = u[0];
    bounds0[3*i+1] = u[1];
    bounds0[3*i+2] = u[2];
  }
  print_bounds(bounds0,"Basemul");
}

void crt_range(int16_t *coeffs0, const int16_t *coeffs1, uint32_t *bounds0, const uint32_t *bounds1) {
  int32_t i;
  int16_t t;
  uint32_t u;

  for(i=0;i<2*KEM_N;i++) {
    pdata = PDATA0;
    coeffs0[i] = mulmod(coeffs0[i],pdata[_16XMONT]);
    bounds0[i] = maxmulmod(bounds0[i],pdata[_16XMONT]);
    pdata = PDATA1;
    t = mulmod(coeffs1[i] - coeffs0[i],CRT_U);
    u = maxmulmod(bounds1[i] + bounds0[i],CRT_U);
    t *= P0;
    u *= P0;
    coeffs0[i] = coeffs0[i] + t;
    bounds0[i] = bounds0[i] + u;
    assert(bounds0[i] < P0*P1 - KEM_Q/2*KEM_N);
  }
}

