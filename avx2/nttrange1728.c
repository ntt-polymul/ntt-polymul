#include <stdint.h>
#include <assert.h>
#include "params.h"
#include "consts1728.h"
#include "nttrange.h"

static void ntt0t1(int16_t *coeffs, uint32_t *bounds) {
  int32_t j,k,l;
  int16_t zeta,zetasq;
  const int16_t omega = pdata[_ZETAS3];
  int32_t t,u,v,w;
  char s[] = "Level 0";

  for(l=NTT_N;l>=NTT_N/3;l/=3) {
    for(j=0;j<NTT_N;j+=l) {
      zeta = pdata[_ZETAS3+2*j/l];
      zetasq = pdata[_ZETAS3+2*(j/l+1)];
      for(k=j;k<j+l/3;k++) {
        if(j > 0) {
          coeffs[l/3+k] = mulmod(coeffs[l/3+k],zeta);
          bounds[l/3+k] = maxmulmod(bounds[l/3+k],zeta);
          coeffs[2*l/3+k] = mulmod(coeffs[2*l/3+k],zetasq);
          bounds[2*l/3+k] = maxmulmod(bounds[2*l/3+k],zetasq);
        }

        t = mulmod(coeffs[l/3+k]-coeffs[2*l/3+k],omega);
        u = coeffs[k] + coeffs[l/3+k] + coeffs[2*l/3+k];
        v = coeffs[k] + t - coeffs[2*l/3+k];
        w = coeffs[k] - coeffs[l/3+k] - t;
        coeffs[k] = u;
        coeffs[l/3+k] = v;
        coeffs[2*l/3+k] = w;

        t = maxmulmod(bounds[l/3+k]+bounds[2*l/3+k],omega);
        u = bounds[k] + bounds[l/3+k] + bounds[2*l/3+k];
        v = bounds[k] + t + bounds[2*l/3+k];
        w = bounds[k] + bounds[l/3+k] + t;
        bounds[k] = u;
        bounds[l/3+k] = v;
        bounds[2*l/3+k] = w;
      }
    }
    print_bounds(bounds,s);
    s[6] += 1;
  }
}

static void twist192(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,idx;
  int16_t zeta;

  for(i=0;i<NTT_N;i++) {
    j  = i;
    idx  = j/192*32;
    j %= 192;
    idx += j/16*288;
    j %= 16;
    idx += j;
    zeta = pdata[_TWIST192+16+idx];
    coeffs[i] = mulmod(coeffs[i],zeta);
    bounds[i] = maxmulmod(bounds[i],zeta);
  }
  print_bounds(bounds,"Twist 192");
}

static void ntt2t4(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,k,l;
  int16_t zeta,t;
  char s[] = "Level 2";

  for(l=NTT_N/9;l>=NTT_N/36;l/=2) {
    for(j=0;j<NTT_N;j+=l) {
      i = j%(NTT_N/9) / l;
      if(i > 0) zeta = pdata[_ZETAS+2*(i-1)];
      else zeta = pdata[_16XMONT];
      for(k=j;k<j+l/2;k++) {
        if(i > 0 || l == NTT_N/36) {
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

static void twist24(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,idx;
  int16_t zeta;

  for(i=0;i<NTT_N;i++) {
    j  = i % 192;
    idx  = j/24*2;
    j %= 24;
    idx += j/2*32;
    j %= 2;
    idx += j;
    zeta = pdata[_TWIST24+16+idx];
    coeffs[i] = mulmod(coeffs[i],zeta);
    bounds[i] = maxmulmod(bounds[i],zeta);
  }
  print_bounds(bounds,"Twist 24");
}

static void ntt5t7(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,k,l;
  int16_t zeta,t;
  char s[] = "Level 5";

  for(l=NTT_N/72;l>=NTT_N/288;l/=2) {
    for(j=0;j<NTT_N;j+=l) {
      i = j%(NTT_N/72) / l;
      if(i > 0) zeta = pdata[_ZETAS+2*(i-1)];
      else zeta = pdata[_16XMONT];
      for(k=j;k<j+l/2;k++) {
        if(i > 0 || l == NTT_N/288) {
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
  ntt0t1(coeffs,bounds);
  twist192(coeffs,bounds);
  ntt2t4(coeffs,bounds);
  twist24(coeffs,bounds);
  ntt5t7(coeffs,bounds);
}

static void invntt0t2(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,k,l;
  int16_t zeta,t;
  char s[] = "Inv Level 0";

  for(l=NTT_F;l<=4*NTT_F;l*=2) {
    for(j=0;j<NTT_N;j+=2*l) {
      const int32_t lut[4] = {_16XMONT-_ZETAS, 0, 4, 2};
      i = j%(8*NTT_F) / (2*l);
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
        if( (l == NTT_F && i == 2 && k >= j+2*l/3)
            || (l == 2*NTT_F && i == 1 && k < j+2*l/3) )
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

static void invtwist24(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,idx;
  int16_t zeta;

  for(i=0;i<NTT_N;i++) {
    j  = i % 192;
    idx  = j/24*2;
    j %= 24;
    idx += j/2*32;
    j %= 2;
    idx += j;
    zeta = pdata[_TWIST24INV+16+idx];
    coeffs[i] = mulmod(coeffs[i],zeta);
    bounds[i] = maxmulmod(bounds[i],zeta);
  }
  print_bounds(bounds,"Inv Twist 24");
}

static void invntt3t5(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,k,l;
  int16_t zeta,t;
  char s[] = "Inv Level 3";

  for(l=8*NTT_F;l<=32*NTT_F;l*=2) {
    for(j=0;j<NTT_N;j+=2*l) {
      const int32_t lut[4] = {_16XMONT-_ZETAS, 0, 4, 2};
      i = j%(64*NTT_F) / (2*l);
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
        if( (l == 8*NTT_F && i == 2 && k >= j+2*l/3)
            || (l == 16*NTT_F && i == 1 && k < j+2*l/3) )
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

static void invtwist192(int16_t *coeffs, uint32_t *bounds) {
  int32_t i,j,idx;
  int16_t zeta;

  for(i=0;i<NTT_N;i++) {
    const int32_t lut[9] = {0,2,1,8,7,6,5,4,3};
    j  = i;
    idx  = lut[j/192]*32;
    j %= 192;
    idx += j/16*288;
    j %= 16;
    idx += j;
    zeta = pdata[_TWIST192+16+idx];
    coeffs[i] = mulmod(coeffs[i],zeta);
    bounds[i] = maxmulmod(bounds[i],zeta);
  }
  print_bounds(bounds,"Inv Twist 192");
}

static void invntt6t7(int16_t *coeffs, uint32_t *bounds) {
  int32_t j,k,l;
  int16_t zeta,zetasq;
  const int16_t omega = pdata[_ZETAS3];
  int32_t t,u,v,w;
  char s[] = "Inv Level 6";

  for(l=64*NTT_F;l<NTT_N;l*=3) {
    for(j=0;j<NTT_N;j+=3*l) {
      zeta = pdata[_ZETAS3_INV+2*j/(3*l)];
      zetasq = pdata[_ZETAS3_INV+2*(j/(3*l)+1)];
      for(k=j;k<j+l;k++) {
        if(l == 192*NTT_F)
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT]);
        if(l == 192*NTT_F && k < l/3) {
          bounds[l+k] = maxmulmod(bounds[l+k],pdata[_16XMONT]);
          bounds[2*l+k] = maxmulmod(bounds[2*l+k],pdata[_16XMONT]);
        }

        t = mulmod(coeffs[l+k]-coeffs[2*l+k],omega);
        u = coeffs[k] + coeffs[l+k] + coeffs[2*l+k];
        v = coeffs[k] + t - coeffs[2*l+k];
        w = coeffs[k] - coeffs[l+k] - t;
        coeffs[k] = u;
        coeffs[l+k] = w;
        coeffs[2*l+k] = v;

        t = maxmulmod(bounds[l+k]+bounds[2*l+k],omega);
        u = bounds[k] + bounds[l+k] + bounds[2*l+k];
        v = bounds[k] + t + bounds[2*l+k];
        w = bounds[k] + bounds[l+k] + t;
        bounds[k] = u;
        bounds[l+k] = w;
        bounds[2*l+k] = v;

        if(j > 0) {
          coeffs[l+k] = mulmod(coeffs[l+k],zeta);
          coeffs[2*l+k] = mulmod(coeffs[2*l+k],zetasq);
          bounds[l+k] = maxmulmod(bounds[l+k],zeta);
          bounds[2*l+k] = maxmulmod(bounds[2*l+k],zetasq);
        }
      }
    }
    print_bounds(bounds,s);
    s[10] += 1;
  }
}

void invntt_range(int16_t *coeffs, uint32_t *bounds) {
  invntt0t2(coeffs,bounds);
  invtwist24(coeffs,bounds);
  invntt3t5(coeffs,bounds);
  invtwist192(coeffs,bounds);
  invntt6t7(coeffs,bounds);
}

void basemul_range(int16_t *coeffs0, const int16_t *coeffs1, uint32_t *bounds0, const uint32_t *bounds1) {
  int32_t i;
  int16_t zeta;
  int16_t t[5];
  uint32_t u[5];

  for(i=0;i<NTT_N/3;i++) {
    zeta = pdata[_ZETAS+2*(i/2-1)];
    zeta = (i%2==0) ? zeta : -zeta;
    coeffs0[3*i+0] = maxmulmod(coeffs0[3*i+0],pdata[_16XMONT]);
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
    //pdata = PDATA0;
    //coeffs0[i] = mulmod(coeffs0[i],pdata[_16XMONT]);
    //bounds0[i] = maxmulmod(bounds0[i],pdata[_16XMONT]);
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
