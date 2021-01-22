//select os
#define LINUX
//#define WIN
//security level
#define LAC128
//modulus
#define Q 251
#define BIG_Q 257024//1024*Q 

#if defined(LAC_LIGHT)
//parameter for LAC_LIGHT
#define STRENGTH "LAC_LIGHT"
#define DIM_N 512
#define SEED_LEN 32
#define PK_LEN 544 //N+SEED_LEN
#define SK_LEN 256 //NUM_ONE*4
#define MESSAGE_LEN 16
#define CIPHER_LEN 648//N+C2_VEC_NUM/2
#define C2_VEC_NUM 272//CODE_LEN*8*2
#define NUM_ONE 64 //number of 1 or -1 in noise vector
#define HASH_TYPE "SHA256"
#define SAMPLE_LEN 245
#endif

#if defined(LAC128)
//parameter for LAC_LIGHT
#define STRENGTH "LAC128"
#define DIM_N 512
#define SEED_LEN 32
#define PK_LEN 544 //N+SEED_LEN
#define SK_LEN 512 //NUM_ONE*4
#define MESSAGE_LEN 16
#define CIPHER_LEN 704//N+C2_VEC_NUM/2
#define C2_VEC_NUM 384//(DATA_LEN+ECC_LEN)*8*2
#define NUM_ONE 128 //number of 1 or -1 in noise vector
#define HASH_TYPE "SHA256"
#define SAMPLE_LEN 590
#endif

#if defined(LAC192)
//parameter for LAC_STANDARD
#define STRENGTH "LAC192"
#define DIM_N 1024
#define SEED_LEN 32
#define PK_LEN 1056//N+SEED_LEN
#define SK_LEN 512 //NUM_ONE*4
#define MESSAGE_LEN 32
#define CIPHER_LEN 1352//N+C2_VEC_NUM/2=1183
#define C2_VEC_NUM 656//(DATA_LEN+ECC_LEN)*8*2
#define NUM_ONE 128 //number of 1 or -1 in noise vector
#define HASH_TYPE "SHA256"
#define SAMPLE_LEN 495
#endif

#if defined(LAC256)
//parameter for LAC_HIGH
#define STRENGTH "LAC256"
#define DIM_N 1024
#define SEED_LEN 32
#define PK_LEN 1056//(N+SEED_LEN)
#define SK_LEN 768 //NUM_ONE*4
#define MESSAGE_LEN 32
#define CIPHER_LEN 1464//(N+C2_VEC_NUM/2)
#define C2_VEC_NUM 880//2*(DATA_LEN+ECC_LEN)*8
#define NUM_ONE 192 //number of 1 or -1 in noise vector
#define HASH_TYPE "SHA256"
#define SAMPLE_LEN 815
#endif
