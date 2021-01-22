#include "rand.h"
#include "lac_param.h"
#include "aes.h"
#include "sha2.h"

#include <string.h>

//pseudo-random bytes
int pseudo_random_bytes(uint8_t *r, unsigned int len, const uint8_t *seed)
{
	uint8_t data[12] = {0};
	aes256ctx ctx;
	aes256_ctr_keyexp(&ctx, seed);
	aes256_ctr(r, len, data, &ctx);
	return 0;
}

//hash
int hash_to_k(const unsigned char *in, unsigned int len_in, unsigned char * out)
{
	unsigned char tmp_out[32];
	sha256(tmp_out, in, len_in);
	memcpy(out,tmp_out,MESSAGE_LEN);
	return 0;
}
