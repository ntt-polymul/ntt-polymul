#include "lac_param.h"

#if defined(LAC_LIGHT)
//bch(255,128,7)
#define DATA_LEN 16
#define ECC_LEN 1
#define CODE_LEN 17
#endif

#if defined(LAC128)
//bch(255,128,15)
#define DATA_LEN 16
#define ECC_LEN 8
#define CODE_LEN 24
#endif

#if defined(LAC192)
//bch(511,256,15)
#define DATA_LEN 32
#define ECC_LEN 9
#define CODE_LEN 41
#endif

#if defined(LAC256)
//D2+bch(511,264,33)
#define DATA_LEN 32
#define ECC_LEN 23
#define CODE_LEN 55 
#endif

//error correction encode
int ecc_enc(const uint8_t *d, uint8_t *c);

//error corrction decode
int ecc_dec(uint8_t *d, const uint8_t *c);

