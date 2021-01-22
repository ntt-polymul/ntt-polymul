#include <stdint.h>

//generate the public parameter a from seed
int gen_a(uint8_t *a,  const uint8_t *seed);
//generate the small random vector for secret and error
int gen_e(uint8_t *e, uint8_t *seed);
int gen_r(uint8_t *r, uint8_t *seed);

// poly_mul  b=[as]
int poly_mul(const uint8_t *a, const uint8_t *s, uint8_t *b, unsigned int  vec_num);
// poly_aff  b=as+e 
int poly_aff(const uint8_t *a, const  uint8_t *s, uint8_t *e, uint8_t *b, unsigned int vec_num);
// compress: cut the low 4bit
int poly_compress(const uint8_t *in,  uint8_t *out, const unsigned int vec_num);
// de-compress: set the low 4bit to be zero
int poly_decompress(const uint8_t *in, uint8_t *out, const unsigned int vec_num);
