#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <db.h>
#include "consts.h"
#include "nttrange.h"
#include "poly.h"

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

static DB *btree = NULL;
const int16_t *pdata;

int16_t mulmod(int16_t a, int16_t b) {
  int32_t t;
  const int16_t p = pdata[_16XP];
  const int16_t pinv = pdata[_16XPINV];

  t = (int32_t)a*b;
  a = (int16_t)t*pinv;
  a = (t - (int32_t)a*p) >> 16;
  return a;
}

int16_t maxmulmod(uint32_t bound, int16_t r) {
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

int16_t maxmulmod2(uint32_t bound0, uint32_t bound1) {
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

void print_bounds(uint32_t *bounds, char *s) {
  int32_t i,j;

  printf("%s:\n", s);
  for(i=0;i<NTT_N/16;i++) {
    //printf(" ");
    for(j=0;j<16;j++)
      printf(" %5d,",bounds[16*i+j]);
    printf("\n");
  }

  for(i=0;i<NTT_N;i++)
    assert(bounds[i] < 32768);
}

void poly_ntt(nttpoly *r, const poly *a, const int16_t *mypdata) {
  int32_t i;
  uint32_t bounds[NTT_N];

  for(i=0;i<POLY_N;i++) {
    r->coeffs[i] = a->coeffs[i];
    bounds[i] = abs(a->coeffs[i]);
  }
  for(i=POLY_N;i<NTT_N;i++) {
    r->coeffs[i] = 0;
    bounds[i] = 0;
  }

  db_create(&btree, NULL, 0);
  btree->open(btree, NULL, "range.db", NULL, DB_BTREE, DB_CREATE, 0664);
  pdata = mypdata;
  ntt_range(r->coeffs,bounds);
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
  btree->open(btree, NULL, "range.db", NULL, DB_BTREE, DB_CREATE, 0664);
  pdata = mypdata;
  invntt_range(r->coeffs,bounds);
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
  btree->open(btree, NULL, "range.db", NULL, DB_BTREE, DB_CREATE, 0664);

  pdata = PDATA0;
  ntt_range(coeffs,bounds0);
  ntt_range(coeffs,bounds1);
  basemul_range(coeffs,coeffs,bounds0,bounds1);
  invntt_range(coeffs,bounds0);

  pdata = PDATA1;
  ntt_range(coeffs,bounds2);
  ntt_range(coeffs,bounds3);
  basemul_range(coeffs,coeffs,bounds2,bounds3);
  invntt_range(coeffs,bounds2);
  crt_range(coeffs,coeffs,bounds0,bounds2);

  btree->close(btree,0);
}
