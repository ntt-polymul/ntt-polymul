#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <immintrin.h>
#include <db.h>
#include "params.h"
#include "consts1024.h"

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

static void ntt(int32_t *bounds, int32_t n, const int16_t *pdata) {
  int32_t j,k,l;

  for(l=n/2;l>=n/512;l/=2) {
    for(j=0;j<n/(2*l);j++) {
      for(k=2*l*j;k<2*l*j+l;k++) {
        if(l == n/64 && (j < 2 || j >= 16))
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT],pdata);
        if(l == n/128 && j >= 4 && j < 32)
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT],pdata);
        if(j > 0 || l <= n/128 || l == n/64)
          bounds[l+k] = maxmulmod(bounds[l+k],pdata[_ZETAS+j],pdata);
        bounds[k] = bounds[l+k] = bounds[k] + bounds[l+k];
      }
    }
  }
}

static void invntt(int32_t *bounds, int32_t n, const int16_t *pdata) {
  int32_t j,k,l,i;

  for(l=n/512;l<=n/2;l*=2) {
    for(j=0;j<n/(2*l);j++) {
      for(k=2*l*j;k<2*l*j+l;k++) {
        bounds[k] = bounds[l+k] = bounds[k] + bounds[l+k];
        if(j > 0 || l <= n/128 || l == n/64) {
          i = (j==0) ? -1431655765 : 1 << (31 - _lzcnt_u32(j));
          i = 3*i - j - 1;
          bounds[l+k] = maxmulmod(bounds[l+k],pdata[_ZETAS+i],pdata);
        }
        i = k - 2*l*j;
        if(l == n/256 || l == n/64)
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT],pdata);
        if(l == n/16 && i < l/2)
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT],pdata);
        if(l == n/8 && j == 0 && i >= l/4 && i < 3*l/4) {
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT],pdata);
          bounds[l+k] = maxmulmod(bounds[l+k],pdata[_16XMONT],pdata);
        }
        if(l == n/8 && j == 0 && i >= 3*l/4)
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT],pdata);
        if(l == n/8 && j > 0 && i >= l/4 && i < 2*l/4)
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT],pdata);
        if(l == n/4 && j == 0 && i >= l/8) {
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT],pdata);
          bounds[l+k] = maxmulmod(bounds[l+k],pdata[_16XMONT],pdata);
        }
        if(l == n/4 && j > 0 && (i < l/8 || (i >= 2*l/8 && i < 4*l/8)))
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT],pdata);
      }
    }
  }
}

static void basemul(int32_t *bounds, const int16_t *pdata) {
  int32_t i, t[3], x, y;

  for(i=0;i<NTT_N/2;i++) {
    x = maxmulmod(bounds[2*i],pdata[_16XMONT],pdata);
    y = maxmulmod(bounds[2*i+1],pdata[_16XMONT],pdata);
    t[0] = maxmulmod(bounds[2*i]*x,1,pdata);
    t[1] = 2*maxmulmod(bounds[2*i]*y,1,pdata);
    t[2] = maxmulmod(bounds[2*i+1]*y,1,pdata);
    t[2] = maxmulmod(t[2],pdata[_ZETAS+i/2],pdata);
    t[0] += t[2];
    t[0] = maxmulmod(t[0],pdata[_16XF],pdata);
    t[1] = maxmulmod(t[1],pdata[_16XF],pdata);
    bounds[2*i+0] = t[0];
    bounds[2*i+1] = t[1];
  }
}

static void poly_mul_modp(int32_t *bounds, const int16_t *pdata) {
  int32_t i,t;

  /* forward ntt */
  ntt(bounds,NTT_N,pdata);
  for(i=t=0;i<NTT_N;i++)
    t = max(bounds[i],t);
  printf("Maximum output coefficient forward NTT: %d\n",t);

  /* basemul */
  basemul(bounds,pdata);


  for(i=t=0;i<NTT_N;i++)
    t = max(bounds[i],t);
  printf("Maximum output coefficient basemul: %d\n",t);

  /* inverse ntt */
  invntt(bounds,NTT_N,pdata);
  for(i=t=0;i<NTT_N;i++)
    t = max(bounds[i],t);
  printf("Maximum output coefficient inverse NTT: %d\n",t);
}

static void crt(int32_t *bounds2, const int32_t *bounds0, const int32_t *bounds1) {
  // first input is reduced, i.e., bounds0 does not matter.
  (void) bounds0;
  int32_t i,t,u;

  for(i=0;i<2*KEM_N;i++) {
    t = (P0-1)/2; //  reduction of bounds0[i]
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
    bounds0[i] = bounds1[i] = 2048;
  for(i=KEM_N;i<NTT_N;i++)
    bounds0[i] = bounds1[i] = 0;

  poly_mul_modp(bounds0,PDATA0);
  poly_mul_modp(bounds1,PDATA1);
  crt(bounds0,bounds0,bounds1);

  btree->close(btree,0);
  return 0;
}
