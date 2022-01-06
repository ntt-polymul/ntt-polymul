#ifndef RANGE_H
#define RANGE_H

#include <stdint.h>

extern const int16_t *pdata;

int16_t mulmod(int16_t a, int16_t b);
int16_t maxmulmod(uint32_t bound, int16_t r);
int16_t maxmulmod2(uint32_t bound0, uint32_t bound1);
void print_bounds(uint32_t *bounds, char *s);
void range_mul(uint32_t *bounds0, uint32_t *bounds1);

void ntt_range(int16_t *coeffs, uint32_t *bounds);
void invntt_range(int16_t *coeffs, uint32_t *bounds);
void basemul_range(int16_t *coeffs0, const int16_t *coeffs1, uint32_t *bounds0, const uint32_t *bounds1);
void crt_range(int16_t *coeffs0, const int16_t *coeffs1, uint32_t *bounds0, const uint32_t *bounds1);

#endif
