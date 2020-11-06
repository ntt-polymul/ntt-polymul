#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <db.h>
#include "params.h"
#include "consts256.h"

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
  if(!btree->get(btree, NULL, &key, &value, 0)) {
    a  = ((uint8_t *)value.data)[0];
    a |= (int16_t)((uint8_t *)value.data)[1] << 8;
    return a;
  }

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

static void ntt_negacyclic(int32_t *bounds, int32_t n, const int16_t *pdata) {
  int32_t i=0,j,k,l;
  int16_t zeta;

  for(l=n/2;l>=n/8;l/=2) {
    for(j=0;j<n/(2*l);j++) {
      const int32_t idx[7] = {0,32,96,64,72,128,136};
      zeta = pdata[_ZETAS+16+idx[i++]];
      for(k=2*l*j;k<2*l*j+l;k++) {
        bounds[l+k] = maxmulmod(bounds[l+k],zeta,pdata);
        bounds[k] = bounds[l+k] = bounds[k] + bounds[l+k];
      }
    }
  }
}

static void ntt(int32_t *bounds, int32_t n, const int16_t *pdata) {
  int32_t j,k,l;

  for(l=n/2;l>=n/8;l/=2) {
    for(j=0;j<n/(2*l);j++) {
      for(k=2*l*j;k<2*l*j+l;k++) {
        const int32_t idx[4] = {-1,0,32,96};
        if(j > 0) bounds[l+k] = maxmulmod(bounds[l+k],pdata[_ZETAS+16+idx[j]],pdata);
        else if(l <= n/4) bounds[l+k] = maxmulmod(bounds[l+k],pdata[_16XMONT],pdata);
        bounds[k] = bounds[l+k] = bounds[k] + bounds[l+k];
      }
    }
  }
}

static void ntt2(int32_t *bounds, int32_t n, const int16_t *pdata) {
  int32_t j,k,l;

  for(l=n/2;l>=n/4;l/=2) {
    for(j=0;j<n/(2*l);j++) {
      for(k=2*l*j;k<2*l*j+l;k++) {
        if(j > 0) bounds[l+k] = maxmulmod(bounds[l+k],pdata[_ZETAS+16],pdata);
        bounds[k] = bounds[l+k] = bounds[k] + bounds[l+k];
      }
    }
  }
}

static void invntt(int32_t *bounds, int32_t n, const int16_t *pdata) {
  int32_t j,k,l;

  for(l=n/4;l<=n/2;l*=2) {
    for(j=0;j<n/(2*l);j++) {
      for(k=2*l*j;k<2*l*j+l;k++) {
        bounds[k] = bounds[l+k] = bounds[k] + bounds[l+k];
        if(j > 0) bounds[l+k] = maxmulmod(bounds[l+k],pdata[_ZETAS+16],pdata);
      }
    }
  }
}

static void invntt2(int32_t *bounds, int32_t n, const int16_t *pdata) {
  int32_t j,k,l;

  for(l=n/8;l<=n/2;l*=2) {
    for(j=0;j<n/(2*l);j++) {
      for(k=2*l*j;k<2*l*j+l;k++) {
        if(l == n/2 && k < n/8)
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT],pdata);
        bounds[k] = bounds[l+k] = bounds[k] + bounds[l+k];
        if(j > 0) {
          const int32_t lut[8] = {-1,0,96,32,136,128,72,64};
          bounds[l+k] = maxmulmod(bounds[l+k],pdata[_ZETAS+16+lut[j]],pdata);
        }
      }
    }
  }
}

static void invntt_negacyclic(int32_t *bounds, int32_t n, const int16_t *pdata) {
  int32_t i,j,k,l;
  int16_t zeta;

  i = 6;
  for(l=n/8;l<=n/2;l*=2) {
    for(j=0;j<n/(2*l);j++) {
      const int32_t idx[7] = {0,32,96,64,72,128,136};
      zeta = pdata[_ZETAS+16+idx[i--]];
      for(k=2*l*j;k<2*l*j+l;k++) {
        if(l == n/2 && k < n/8) {
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT],pdata);
          bounds[l+k] = maxmulmod(bounds[l+k],pdata[_16XMONT],pdata);
        }
        bounds[k] = bounds[l+k] = bounds[k] + bounds[l+k];
        bounds[l+k] = maxmulmod(bounds[l+k],zeta,pdata);
      }
    }
  }
}

static void basemul_acc(int32_t *bounds, int32_t n, const int16_t *pdata) {
  int32_t i, x;

  for(i=0;i<n;i++) {
    x = maxmulmod(bounds[i],pdata[_16XMONT],pdata);
    bounds[i] = 4*maxmulmod(bounds[i]*x,1,pdata);
    bounds[i] = maxmulmod(bounds[i],pdata[_16XF],pdata);
  }
}

static void poly_mul_modp(int32_t *bounds, const int16_t *pdata) {
  int32_t i,j,idx,t;

  /* forward ntt */
  ntt_negacyclic(bounds,NTT_N,pdata);
  for(i=0;i<NTT_N;i++) {
    j  = i;
    idx  = j/128*256;
    j %= 128;
    idx += j/32*4;
    j %= 32;
    idx += j/4*32;
    j %= 4;
    idx += j;
    bounds[i] = maxmulmod(bounds[i],pdata[_TWIST32+16+idx],pdata);
  }
  for(i=0;i<NTT_N;i+=NTT_N/8)
    ntt(bounds+i,NTT_N/8,pdata);
  for(i=0;i<NTT_N;i++) {
    j  = i % 32;
    idx = j/4*8;
    j %= 4;
    idx += j;
    if(idx < 4) bounds[i] = maxmulmod(bounds[i],pdata[_16XMONT],pdata);
    else bounds[i] = maxmulmod(bounds[i],pdata[_TWIST4+4+idx],pdata);
  }
  for(i=0;i<NTT_N;i+=NTT_N/64)
    ntt2(bounds+i,NTT_N/64,pdata);

  for(i=t=0;i<NTT_N;i++)
    t = max(bounds[i],t);
  printf("Maximum output coefficient forward NTT: %d\n",t);

  /* basemul */
  basemul_acc(bounds,NTT_N,pdata);

  //for(i=0;i<NTT_N;i++)
  //  bounds[i] = maxmulmod((1<<15)*(p-1)/2,1,pdata);

  for(i=t=0;i<NTT_N;i++)
    t = max(bounds[i],t);
  printf("Maximum output coefficient basemul: %d\n",t);

  /* inverse ntt */
  for(i=0;i<NTT_N;i+=NTT_N/64)
    invntt(bounds+i,NTT_N/64,pdata);
  for(i=0;i<NTT_N;i++) {
    const int32_t lut[8] = {-1, 0, 3, 2, 7, 6, 5, 4};
    j  = i % 32;
    idx  = lut[j/4]*8;
    j %= 4;
    idx += j;
    if(idx < 0) bounds[i] = maxmulmod(bounds[i],pdata[_16XMONT],pdata);
    else bounds[i] = maxmulmod(bounds[i],pdata[_TWIST4+4+idx],pdata);
  }
  for(i=0;i<NTT_N;i+=NTT_N/8)
    invntt2(bounds+i,NTT_N/8,pdata);
  for(i=0;i<NTT_N;i++) {
    j  = i;
    idx  = j/128*256;
    j %= 128;
    idx += j/32*4;
    j %= 32;
    idx += (7-j/4)*32;
    j %= 4;
    idx += j;
    bounds[i] = maxmulmod(bounds[i],pdata[_TWIST32+16+idx],pdata);
  }
  invntt_negacyclic(bounds,NTT_N,pdata);

  for(i=t=0;i<NTT_N;i++)
    t = max(bounds[i],t);
  printf("Maximum output coefficient inverse NTT: %d\n",t);
}

static void crt(int32_t *bounds2, const int32_t *bounds0, const int32_t *bounds1) {
  int32_t i,t,u;

  for(i=0;i<KEM_N;i++) {
    t = maxmulmod(bounds0[i],PDATA0[_16XMONT],PDATA0);
    u = maxmulmod(bounds1[i] + t,CRT_U,PDATA1);
    u *= P0;
    bounds2[i] = t + u;
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
