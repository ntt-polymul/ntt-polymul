#include <string.h>
#include "api.h"
#include "rand.h"
#include "bch.h"
#include "ecc.h"

//generate keypair
int crypto_kem_keypair(uint8_t *pk, uint8_t *sk)
{
	//call the key generation algorithm of pke
	crypto_encrypt_keypair(pk, sk);
	return 0;
}
int crypto_kem_enc( uint8_t *ct, uint8_t *ss, const uint8_t *pk)
{
	kem_enc_fo(pk,ss,ct);
	return 0;
}
int crypto_kem_dec( uint8_t *ss, const uint8_t *ct, const uint8_t *sk)
{
	const uint8_t *pk;
	pk=sk+SK_LEN;
	kem_dec_fo(pk,sk,ct,ss);
	return 0;
}
// fo encryption for cca security 
int kem_enc_fo(const uint8_t *pk, uint8_t *k, uint8_t *c)
{
	uint8_t buf[MESSAGE_LEN],seed[SEED_LEN],seed_buf[MESSAGE_LEN+SEED_LEN];
	unsigned long long clen;
	
	//generate random message m, stored in buf
	random_bytes(buf,MESSAGE_LEN);
	//compute seed=gen_seed(m|pk), add pk for multi key attack protection
	memcpy(seed_buf,buf,MESSAGE_LEN);
	memcpy(seed_buf+MESSAGE_LEN,pk,SEED_LEN);
	gen_seed(seed_buf,MESSAGE_LEN+SEED_LEN,seed);
	//encrypt m with seed
	pke_enc_seed(pk,buf,MESSAGE_LEN,c,&clen,seed);
	
	//compute k=hash(m|c)
	hash_to_k(buf,MESSAGE_LEN,k);
	
	return 0;
}

// fo encryption for cca security with seed
int kem_enc_fo_seed(const uint8_t *pk, uint8_t *k, uint8_t *c, uint8_t *seed)
{
	uint8_t buf[MESSAGE_LEN],local_seed[SEED_LEN],seed_buf[MESSAGE_LEN+SEED_LEN];
	unsigned long long clen;
	
	//generate random message m from seed, stored in buf
	pseudo_random_bytes(buf,MESSAGE_LEN,seed);
	//compute loacal_seed=gen_seed(m|pk), add pk for multi key attack protection
	memcpy(seed_buf,buf,MESSAGE_LEN);
	memcpy(seed_buf+MESSAGE_LEN,pk,SEED_LEN);
	gen_seed(seed_buf,MESSAGE_LEN+SEED_LEN,local_seed);
	//encrypt m with local_seed
	pke_enc_seed(pk,buf,MESSAGE_LEN,c,&clen,local_seed);
	
	//compute k=hash(m|c)
	hash_to_k(buf,MESSAGE_LEN,k);
	
	return 0;
}

// decrypt of fo mode
int kem_dec_fo(const uint8_t *pk, const uint8_t *sk, const uint8_t *c, uint8_t *k)
{
	uint8_t buf[MESSAGE_LEN+CIPHER_LEN],seed[SEED_LEN],seed_buf[MESSAGE_LEN+SEED_LEN];
	unsigned long long mlen,clen;
	uint8_t c_v[CIPHER_LEN];
	
	//compute m from c
	pke_dec(sk,c,CIPHER_LEN, buf,&mlen);
	//compte k=hash(m|c)
	hash_to_k(buf,MESSAGE_LEN,k);
	//re-encryption with seed=gen_seed(m|pk), add pk for multi key attack protection
	memcpy(seed_buf,buf,MESSAGE_LEN);
	memcpy(seed_buf+MESSAGE_LEN,pk,SEED_LEN);
	gen_seed(seed_buf,MESSAGE_LEN+SEED_LEN,seed);
	pke_enc_seed(pk,buf,MESSAGE_LEN,c_v,&clen,seed);
	
	//verify
	if(memcmp(c,c_v,CIPHER_LEN)!=0)
	{
		//k=hash(hash(sk)|c)
		hash((uint8_t*)sk,SK_LEN,buf);
		hash(buf,MESSAGE_LEN+CIPHER_LEN,k);
	}
	
	return 0;
}

