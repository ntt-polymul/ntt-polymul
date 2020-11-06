#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <db.h>
#include "params.h"
#include "consts1728.h"

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

static DB *btree = NULL;

static int16_t maxmulmod(int32_t bound, int16_t r, const int16_t *pdata) {
  int32_t i, t;
  int16_t a, b, m;
  uint64_t ikey;
  DBT key, value;
  const int16_t p = pdata[_16XP];
  const int16_t pinv = pdata[_16XPINV];

  assert(bound >= 0);
  if(r != 1) assert(bound <= 32768);

  memset(&key, 0, sizeof(key));
  memset(&value, 0, sizeof(value));
  ikey = ((uint64_t)r << 48) | ((uint64_t)bound << 16) | ((uint64_t)p << 0);
  key.data = &ikey;
  key.size = 8;
  if(!btree->get(btree, NULL, &key, &value, 0))
    return *(int16_t *)value.data;

  a = b = 0;
  for(i=-bound;i<bound;i++) {
    t = (int32_t)i*r;
    m = (int16_t)t*pinv;
    t = (t - (int32_t)m*p) >> 16;
    a = min(a,t);
    b = max(b,t);
  }

  a = max(-a,b);
  value.data = &a;
  value.size = 2;
  btree->put(btree, NULL, &key, &value, 0);
  return a;
}

static void ntt(int32_t *bounds, int32_t n, const int16_t *pdata) {
  int32_t j,k,l;
  int32_t t,u,x,y,z;

  for(l=n/3;l>=n/9;l/=3) {
    for(j=0;j<n/(3*l);j++) {
      for(k=3*l*j;k<3*l*j+l;k++) {
        if(j > 0) {
          bounds[l+k] = maxmulmod(bounds[l+k],pdata[_ZETAS3+2*j],pdata);
          bounds[2*l+k] = maxmulmod(bounds[2*l+k],pdata[_ZETAS3+2*(j+1)],pdata);
        }
        t = maxmulmod(bounds[l+k],pdata[_ZETAS3],pdata);
        u = maxmulmod(bounds[2*l+k],pdata[_ZETAS3],pdata);
        x = bounds[k] + bounds[l+k] + bounds[2*l+k];
        y = bounds[k] + t + (bounds[2*l+k] + u);
        z = bounds[k] + (bounds[l+k] + t) + u;
        bounds[k] = x;
        bounds[l+k] = y;
        bounds[2*l+k] = z;
      }
    }
  }
}

static void ntt2(int32_t *bounds, int32_t n, const int16_t *pdata) {
  int32_t j,k,l;

  for(l=n/2;l>=n/8;l/=2) {
    for(j=0;j<n/(2*l);j++) {
      for(k=2*l*j;k<2*l*j+l;k++) {
        if(l == n/8 && j == 0)
          bounds[l+k] = maxmulmod(bounds[l+k],pdata[_16XMONT],pdata);
        if(j > 0) bounds[l+k] = maxmulmod(bounds[l+k],pdata[_ZETAS+2*(j-1)],pdata);
        bounds[k] = bounds[l+k] = bounds[k] + bounds[l+k];
      }
    }
  }
}

static void invntt(int32_t *bounds, int32_t n, const int16_t *pdata) {
  int32_t j,k,l;

  for(l=n/8;l<=n/2;l*=2) {
    for(j=0;j<n/(2*l);j++) {
      for(k=2*l*j;k<2*l*j+l;k++) {
        if(l == n/2 && k < n/8)
          bounds[l+k] = maxmulmod(bounds[l+k],pdata[_16XMONT],pdata);
        bounds[k] = bounds[l+k] = bounds[k] + bounds[l+k];
        if(j > 0) {
          const int32_t lut[4] = {0, 0, 4, 2};
          bounds[l+k] = maxmulmod(bounds[l+k],pdata[_ZETAS+lut[j]],pdata);
        }
      }
    }
  }
}

static void invntt2(int32_t *bounds, int32_t n, const int16_t *pdata) {
  int32_t j,k,l;
  int32_t t,u,x,y,z;

  for(l=n/9;l<=n/3;l*=3) {
    for(j=0;j<n/(3*l);j++) {
      for(k=3*l*j;k<3*l*j+l;k++) {
        if(l == n/3)
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT],pdata);
        if(l == n/3 && k < n/9) {
          bounds[l+k] = maxmulmod(bounds[l+k],pdata[_16XMONT],pdata);
          bounds[2*l+k] = maxmulmod(bounds[2*l+k],pdata[_16XMONT],pdata);
        }
        t = maxmulmod(bounds[l+k],pdata[_ZETAS3],pdata);
        u = maxmulmod(bounds[2*l+k],pdata[_ZETAS3],pdata);
        x = bounds[k] + bounds[l+k] + bounds[2*l+k];
        y = bounds[k] + t + (bounds[2*l+k] + u);
        z = bounds[k] + (bounds[l+k] + t) + u;
        bounds[k] = x;
        bounds[l+k] = y;
        bounds[2*l+k] = z;
        if(j > 0) {
          bounds[l+k] = maxmulmod(bounds[l+k],pdata[_ZETAS3_INV+2*j],pdata);
          bounds[2*l+k] = maxmulmod(bounds[2*l+k],pdata[_ZETAS3_INV+2*(j+1)],pdata);
        }
      }
    }
  }
}

static void basemul(int32_t *bounds, const int16_t *pdata) {
  int32_t i, t[5], x, y, z;

  for(i=0;i<8;i++) {
    x = maxmulmod(bounds[3*i],pdata[_16XMONT],pdata);
    y = maxmulmod(bounds[3*i+1],pdata[_16XMONT],pdata);
    z = maxmulmod(bounds[3*i+2],pdata[_16XMONT],pdata);
    t[0]  = maxmulmod(bounds[3*i]*x,1,pdata);
    t[1]  = 2*maxmulmod(bounds[3*i]*y,1,pdata);
    t[2]  = 2*maxmulmod(bounds[3*i]*z,1,pdata);
    t[2] += maxmulmod(bounds[3*i+1]*y,1,pdata);
    t[3]  = 2*maxmulmod(bounds[3*i+1]*z,1,pdata);
    t[4]  = maxmulmod(bounds[3*i+2]*z,1,pdata);
    if(i >= 2) {
      t[3] = maxmulmod(t[3],pdata[_ZETAS+2*(i/2-1)],pdata);
      t[4] = maxmulmod(t[3],pdata[_ZETAS+2*(i/2-1)],pdata);
    }
    t[0] += t[3];
    t[1] += t[4];
    t[0] = maxmulmod(t[0],pdata[_16XF],pdata);
    t[1] = maxmulmod(t[1],pdata[_16XF],pdata);
    t[2] = maxmulmod(t[2],pdata[_16XF],pdata);
    bounds[3*i+0] = t[0];
    bounds[3*i+1] = t[1];
    bounds[3*i+2] = t[2];
  }
}

static void poly_mul_modp(int32_t *bounds, const int16_t *pdata) {
  int32_t i,j,idx,t;

  /* forward ntt */
  ntt(bounds,NTT_N,pdata);
  for(i=0;i<NTT_N;i++) {
    j  = i;
    idx  = j/192*32;
    j %= 192;
    idx += j/16*288;
    j %= 16;
    idx += j;
    bounds[i] = maxmulmod(bounds[i],pdata[_TWIST192+16+idx],pdata);
  }
  for(i=0;i<NTT_N;i+=NTT_N/9)
    ntt2(bounds+i,NTT_N/9,pdata);
  for(i=0;i<NTT_N;i++) {
    j  = i % 192;
    idx  = j/24*2;
    j %= 24;
    idx += j/2*32;
    j %= 2;
    idx += j;
    bounds[i] = maxmulmod(bounds[i],pdata[_TWIST24+16+idx],pdata);
  }
  for(i=0;i<NTT_N;i+=NTT_N/72)
    ntt2(bounds+i,NTT_N/72,pdata);

  for(i=t=0;i<NTT_N;i++)
    t = max(bounds[i],t);
  printf("Maximum output coefficient forward NTT: %d\n",t);

  /* basemul */
  for(i=0;i<NTT_N;i+=24)
    basemul(bounds+i,pdata);

  //for(i=0;i<NTT_N;i++) {
  //  int16_t p = pdata[_16XP];
  //  bounds[i] = maxmulmod((1<<15)*(p-1)/2,1,pdata);
  //}

  for(i=t=0;i<NTT_N;i++)
    t = max(bounds[i],t);
  printf("Maximum output coefficient basemul: %d\n",t);

  /* inverse ntt */
  for(i=0;i<NTT_N;i+=NTT_N/72)
    invntt(bounds+i,NTT_N/72,pdata);
  for(i=0;i<NTT_N;i++) {
    j  = i % 192;
    idx  = j/24*2;
    j %= 24;
    idx += j/2*32;
    j %= 2;
    idx += j;
    bounds[i] = maxmulmod(bounds[i],pdata[_TWIST24INV+16+idx],pdata);
  }
  for(i=0;i<NTT_N;i+=NTT_N/9)
    invntt(bounds+i,NTT_N/9,pdata);
  for(i=0;i<NTT_N;i++) {
    const int32_t lut[9] = {0,2,1,8,7,6,5,4,3};
    j  = i;
    idx  = lut[j/192]*32;
    j %= 192;
    idx += j/16*288;
    j %= 16;
    idx += j;
    bounds[i] = maxmulmod(bounds[i],pdata[_TWIST192+16+idx],pdata);
  }
  for(i=0;i<NTT_N;i+=NTT_N)
    invntt2(bounds+i,NTT_N,pdata);

  for(i=t=0;i<NTT_N;i++)
    t = max(bounds[i],t);
  printf("Maximum output coefficient inverse NTT: %d\n",t);
}

static void crt(int32_t *bounds2, const int32_t *bounds0, const int32_t *bounds1) {
  int32_t i,t;

  for(i=0;i<2*KEM_N;i++) {
    t = maxmulmod(bounds1[i] + bounds0[i],CRT_U,PDATA1);
    t *= P0;
    bounds2[i] = bounds0[i] + t;
  }

  for(i=t=0;i<2*KEM_N;i++)
    t = max(bounds2[i],t);
  printf("Maximum output coefficient CRT: %d\n",t);
}

int main(void) {
  int32_t i;
  int32_t bounds0[NTT_N], bounds1[NTT_N];

  db_create(&btree, NULL, 0);
  btree->open(btree, NULL, "access.db", NULL, DB_BTREE, DB_CREATE, 0664);

  for(i=0;i<KEM_N;i++)
    bounds0[i] = bounds1[i] = 4096;
  for(i=KEM_N;i<NTT_N;i++)
    bounds0[i] = bounds1[i] = 0;

  poly_mul_modp(bounds0,PDATA0);
  poly_mul_modp(bounds1,PDATA1);
  crt(bounds0,bounds0,bounds1);

  btree->close(btree,0);
  return 0;
}
