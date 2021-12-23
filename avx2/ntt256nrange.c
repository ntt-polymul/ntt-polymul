#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <db.h>
#include "params.h"
#include "consts256.h"
#include "poly.h"

void range_mul(uint32_t *bounds0, uint32_t *bounds1);

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

static DB *btree = NULL;
static const int16_t *pdata;

static int16_t mulmod(int16_t a, int16_t b) {
  int32_t t;
  const int16_t p = pdata[_16XP];
  const int16_t pinv = pdata[_16XPINV];

  t = (int32_t)a*b;
  a = (int16_t)t*pinv;
  a = (t - (int32_t)a*p) >> 16;
  return a;
}

static int16_t maxmulmod(uint32_t bound, int16_t r) {
  int32_t i;
  int16_t a,b,t,buf[3];
  DBT key, value;
  const int16_t p = pdata[_16XP];

  assert(bound < 32768);

  memset(&key, 0, sizeof(key));
  memset(&value, 0, sizeof(value));
  buf[0] = p;
  buf[1] = bound;
  buf[2] = r;
  key.data = buf;
  key.size = 6;
  if(!btree->get(btree, NULL, &key, &value, 0)) {
    a  = ((uint8_t *)value.data)[0];
    a |= (int16_t)((uint8_t *)value.data)[1] << 8;
    return a;
  }

  a = b = 0;
  for(i=-bound;i<(int32_t)bound;i++) {
    t = mulmod(i,r);
    a = min(a,t);
    b = max(b,t);
  }

  a = max(-a,b);
  value.data = &a;
  value.size = 2;
  btree->put(btree, NULL, &key, &value, 0);
  return a;
}

static int16_t maxmulmod2(uint32_t bound0, uint32_t bound1) {
  int32_t i,j;
  int16_t a,b,t,buf[3];
  DBT key, value;
  const int16_t p = pdata[_16XP];

  assert(bound0 < 32768);
  assert(bound1 < 32768);

  memset(&key, 0, sizeof(key));
  memset(&value, 0, sizeof(value));
  buf[0] = -p;
  buf[1] = bound0;
  buf[2] = bound1;
  key.data = buf;
  key.size = 6;
  if(!btree->get(btree, NULL, &key, &value, 0)) {
    a  = ((uint8_t *)value.data)[0];
    a |= (int16_t)((uint8_t *)value.data)[1] << 8;
    return a;
  }

  a = b = 0;
  for(i=-bound0;i<=(int32_t)bound0;i++) {
    for(j=-bound1;j<=(int32_t)bound1;j++) {
      t = mulmod(i,j);
      a = min(a,t);
      b = max(b,t);
    }
  }

  a = max(-a,b);
  value.data = &a;
  value.size = 2;
  btree->put(btree, NULL, &key, &value, 0);
  return a;
}

static void print_bounds(uint32_t *bounds, char *s) {
  int32_t i;

  for(i=0;i<NTT_N;i++)
    assert(bounds[i] < 32768);

  printf("%s:\n  %d", s, bounds[0]);
  for(i=1;i<NTT_N;i++)
    printf(", %d", bounds[i]);
  printf("\n");
}

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

static void ntt(int16_t *coeffs, uint32_t *bounds) {
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
        if(l == NTT_N && k < NTT_N/8) {  // extra reduction
          coeffs[k] = mulmod(coeffs[k],pdata[_16XMONT]);
          coeffs[l/2+k] = mulmod(coeffs[l/2+k],pdata[_16XMONT]);
          bounds[k] = maxmulmod(bounds[k],pdata[_16XMONT]);
          bounds[l/2+k] = maxmulmod(bounds[l/2+k],pdata[_16XMONT]);
        }

        t = coeffs[l/2+k];
        coeffs[l/2+k] = coeffs[k] - t;
        coeffs[k] = coeffs[k] + t;
        coeffs[l/2+k] = mulmod(coeffs[l/2+k],zeta);

        bounds[k] = bounds[l/2+k] = bounds[k] + bounds[l/2+k];
        bounds[l/2+k] = maxmulmod(bounds[l/2+k],zeta);
      }
    }
    print_bounds(bounds,s);
    s[10] += 1;
  }
}

static void invntt(int16_t *coeffs, uint32_t *bounds) {
  invntt0t1(coeffs,bounds);
  invtwist4(coeffs,bounds);
  invntt2t4(coeffs,bounds);
  invtwist32(coeffs,bounds);
  invntt5t7(coeffs,bounds);
}

static void basemul(int16_t *coeffs0, int16_t *coeffs1, uint32_t *bounds0, uint32_t *bounds1) {
  int32_t i;

  for(i=0;i<NTT_N;i++) {
    coeffs0[i] = mulmod(coeffs0[i],coeffs1[i]);
    coeffs0[i] = mulmod(coeffs0[i],pdata[_16XF]);

    bounds0[i] = 2*maxmulmod2(bounds0[i],bounds1[i]);  // accumulation
    bounds0[i] = maxmulmod(bounds0[i],pdata[_16XF]);
  }
  print_bounds(bounds0,"Basemul");
}

static void crt(uint32_t *bounds2, const uint32_t *bounds0, const uint32_t *bounds1) {
  int32_t i,t,u;

  for(i=0;i<NTT_N;i++) {
    pdata = PDATA0;
    t = maxmulmod(bounds0[i],PDATA0[_16XMONT]);
    pdata = PDATA1;
    u = maxmulmod(bounds1[i] + t,CRT_U);
    u *= P0;
    bounds2[i] = t + u;
    assert(bounds2[i] < P0*P1-5*4096*256*4);
  }
}

void poly_ntt(nttpoly *r, const poly *a, const int16_t *mypdata) {
  int32_t i;
  uint32_t bounds[NTT_N];

  for(i=0;i<NTT_N;i++) {
    r->coeffs[i] = a->coeffs[i];
    bounds[i] = abs(a->coeffs[i]);
  }

  db_create(&btree, NULL, 0);
  btree->open(btree, NULL, "access.db", NULL, DB_BTREE, DB_CREATE, 0664);
  pdata = mypdata;
  ntt(r->coeffs,bounds);
  btree->close(btree,0);
}

void poly_invntt_tomont(nttpoly *r, const nttpoly *a, const int16_t *mypdata) {
  int32_t i;
  uint32_t bounds[NTT_N];

  for(i=0;i<NTT_N;i++) {
    r->coeffs[i] = a->coeffs[i];
    bounds[i] = abs(a->coeffs[i]);
  }

  db_create(&btree, NULL, 0);
  btree->open(btree, NULL, "access.db", NULL, DB_BTREE, DB_CREATE, 0664);
  pdata = mypdata;
  invntt(r->coeffs,bounds);
  btree->close(btree,0);
}

void range_mul(uint32_t *bounds0, uint32_t *bounds1) {
  int32_t i;
  int16_t coeffs[NTT_N];
  uint32_t bounds2[NTT_N], bounds3[NTT_N];

  for(i=0;i<NTT_N;i++) {
    bounds2[i] = bounds0[i];
    bounds3[i] = bounds1[i];
  }

  db_create(&btree, NULL, 0);
  btree->open(btree, NULL, "access.db", NULL, DB_BTREE, DB_CREATE, 0664);

  pdata = PDATA0;
  ntt(coeffs,bounds0);
  ntt(coeffs,bounds1);
  basemul(coeffs,coeffs,bounds0,bounds1);
  invntt(coeffs,bounds0);

  pdata = PDATA1;
  ntt(coeffs,bounds2);
  ntt(coeffs,bounds3);
  basemul(coeffs,coeffs,bounds2,bounds3);
  invntt(coeffs,bounds2);

  crt(bounds0,bounds0,bounds2);

  btree->close(btree,0);
}
