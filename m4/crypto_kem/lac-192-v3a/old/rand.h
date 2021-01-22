#include "randombytes.h"
#include "sha2.h"

//random bytes
#define random_bytes(R, LEN) randombytes(R, LEN)
//hash
#define hash(IN, INBYTES, OUT) sha256(OUT, IN, INBYTES)
//generate seed
#define gen_seed(IN, INBYTES, OUT) sha256(OUT, IN, INBYTES)

//pseudo-random bytes
int pseudo_random_bytes(uint8_t *r, unsigned int len, const uint8_t *seed);
//hash
int hash_to_k(const unsigned char *in, unsigned int len_in, unsigned char * out);
