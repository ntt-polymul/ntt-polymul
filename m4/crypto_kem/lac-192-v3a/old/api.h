#ifndef api_h
#define api_h

#include <stdint.h>

#include  "lac_param.h"
//  Set these three values apropriately for your algorithm
#define CRYPTO_SECRETKEYBYTES SK_LEN+PK_LEN
#define CRYPTO_PUBLICKEYBYTES PK_LEN
#define CRYPTO_BYTES MESSAGE_LEN
#define CRYPTO_CIPHERTEXTBYTES CIPHER_LEN

// Change the algorithm name
#define CRYPTO_ALGNAME STRENGTH
//functions for pke
int crypto_encrypt_keypair( uint8_t *pk, uint8_t *sk);
int crypto_encrypt( uint8_t *c, unsigned long long *clen, const uint8_t *m, unsigned long long mlen, const uint8_t *pk);
int crypto_encrypt_open(uint8_t *m, unsigned long long *mlen,const uint8_t *c, unsigned long long clen,const uint8_t *sk);
//key generation
int kg(uint8_t *pk, uint8_t *sk);
//key generation with seed
int kg_seed(uint8_t *pk, uint8_t *sk, uint8_t *seed);
// encryption
int pke_enc(const uint8_t *pk, const uint8_t *m, unsigned long long mlen, uint8_t *c, unsigned long long *clen);
// encryption with seed
int pke_enc_seed(const uint8_t *pk, const uint8_t *m, unsigned long long mlen, uint8_t *c, unsigned long long *clen, uint8_t *seed);

// decrypt
int pke_dec(const uint8_t *sk, const uint8_t *c, unsigned long long clen, uint8_t *m, unsigned long long *mlen);

//functions for kem
int crypto_kem_keypair( uint8_t *pk, uint8_t *sk);
int crypto_kem_enc( uint8_t *ct, uint8_t *ss, const uint8_t *pk);
int crypto_kem_dec( uint8_t *ss, const uint8_t *ct, const uint8_t *sk);

int kem_enc_fo(const uint8_t *pk, uint8_t *k, uint8_t *c);
// fo encryption for cca security with seed
int kem_enc_fo_seed(const uint8_t *pk, uint8_t *k, uint8_t *c, uint8_t *seed);
// decrypt of fo mode
int kem_dec_fo(const uint8_t *pk, const uint8_t *sk, const  uint8_t *c, uint8_t *k);

//functions for ke
//Alice send: generate pk and sk, and send pk to Bob
// int crypto_ke_alice_send(uint8_t *pk,uint8_t *sk);
// Bob receive: receive  pk, randomly choose m, and encryrpt m with pk to generate c, k=HASH(pk,m).
// int crypto_ke_bob_receive(uint8_t *pk, uint8_t *c, uint8_t *k);
//Alice receive: receive c, and decrypt to get m and comute k=HASH(pk,m)
// int crypto_ke_alice_receive(uint8_t *pk, uint8_t *sk, uint8_t *c, uint8_t *k);

//functions for ake
//Alice send: generate pk and sk, and send pk to Bob
// int crypto_ake_alice_send(uint8_t *pk,uint8_t *sk, uint8_t *pk_b, uint8_t *sk_a, uint8_t *c, uint8_t *k1);
// Bob receive: receive  pk, randomly choose m, and encryrpt m with pk to generate c1 c2, k=HASH(pk,m).
// int crypto_ake_bob_receive(uint8_t *pk_b, uint8_t *sk_b, uint8_t *pk_a, uint8_t *pk, uint8_t *c_in, uint8_t *c_out, uint8_t *k);
//Alice receive: receive c1,c2, and decrypt to get m and comute k=HASH(pk,m)
// int crypto_ake_alice_receive(uint8_t *pk_a, uint8_t *sk_a,uint8_t *pk_b, uint8_t *pk, uint8_t *sk, uint8_t *c1, uint8_t *c_in, uint8_t *k1, uint8_t *k);

#endif /* api_h */
