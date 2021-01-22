#include <string.h>
#include "api.h"
#include "rand.h"
#include "bch.h"
#include "ecc.h"
#include "bin-lwe.h"

#define RATIO 125

static int d2_decode(uint8_t *v1, uint8_t *v2, uint8_t *data, int len)
{
	int i;
	int temp1,temp2;
	int d2_bound=Q/2;
	int center_point=Q/2;
	int vec_bound=len/2;
	
	for(i=0;i<vec_bound;i++)
	{
		//D2 decoding:compute m*q/2+e1 + m*q/2+e2 in [0,2*Q]
		temp1=(v1[i]-v2[i]+Q)%Q;
		temp2=(v1[i+vec_bound]-v2[i+vec_bound]+Q)%Q;
		
		//shift
		if(temp1<center_point)
		{
			temp1=center_point-temp1+center_point;//mirror around Q/2
		}
		if(temp2<center_point)
		{
			temp2=center_point-temp2+center_point;//mirror around Q/2
		}
		//merge erors
		temp1+=temp2-Q;

		//recover m from m*q/2+e1 + m*q/2+e2, RATIO=q/2
		if(temp1<d2_bound)
		{
			data[i/8]=data[i/8]^(1<<(i%8));
		}
	}
	
	return 0;
}

//key generation
int crypto_encrypt_keypair( uint8_t *pk, uint8_t *sk)
{
	kg(pk,sk);
	
	return 0;
}

//encryption
int crypto_encrypt( uint8_t *c, unsigned long long *clen, const uint8_t *m, unsigned long long mlen, const uint8_t *pk)
{
	if(mlen>MESSAGE_LEN)
	{
		return -1;
	}
	
	//call pke encryption function
	pke_enc(pk,m, mlen,c,clen);

	return 0;
}
//decryption
int crypto_encrypt_open(uint8_t *m, unsigned long long *mlen,const uint8_t *c, unsigned long long clen,const uint8_t *sk)
{
	//call pke decryption function
	pke_dec(sk,c,clen,m,mlen);

	return 0;
}

//key generation with seed
int kg_seed(uint8_t *pk, uint8_t *sk, uint8_t *seed)
{
	uint8_t seeds[3*SEED_LEN];
	uint8_t a[DIM_N];
	uint8_t e[DIM_N];
	//char outs[128];
	
	//generate three seeds for a,sk,e
	pseudo_random_bytes(seeds,3*SEED_LEN,seed);	
	
	//copy the seed of a to pk
	memcpy(pk,seeds,SEED_LEN);
	
	//generate a
	gen_a(a,pk);//print_bytes(sk,CRYPTO_SECRETKEYBYTES);
	
	//generate  sk,e
	gen_r(sk,seeds+SEED_LEN);
	gen_e(e,seeds+2*SEED_LEN);
	//compute pk=a*sk+e
	//unsigned long long int t1 = hal_get_time();
	poly_aff(a,sk,e,pk+SEED_LEN,DIM_N);
	//unsigned long long int t2 = hal_get_time();
	//sprintf(outs, "poly_aff(a,sk,e,pk+SEED_LEN,DIM_N):%llu\n",t2-t1 );
	//hal_send_str(outs);
	//copy pk=as+e to the second part of sk, now sk=s|pk
	memcpy(sk+SK_LEN,pk,PK_LEN);
	
	return 0;
}

//key generation
int kg(uint8_t *pk, uint8_t *sk)
{
	uint8_t seed[SEED_LEN];
	
	//generate seed
	random_bytes(seed,SEED_LEN);		
	//key generation with seed 
	kg_seed(pk,sk,seed);	
	
	return 0;
}
// encryption
int pke_enc(const uint8_t *pk, const uint8_t *m, unsigned long long mlen, uint8_t *c, unsigned long long *clen)
{
	uint8_t seed[SEED_LEN];
	
	//generate seed
	random_bytes(seed,SEED_LEN);
	//encrypt with seed 
	pke_enc_seed(pk,m,mlen,c,clen,seed);	

	return 0;
}


// decrypt
int pke_dec(const uint8_t *sk, const uint8_t *c,unsigned long long clen, uint8_t *m, unsigned long long *mlen)
{
	//char outs[128];
	uint8_t out[DIM_N];
	uint8_t c2[C2_VEC_NUM];
	int c2_len=(clen-DIM_N)*2;
	uint8_t code[CODE_LEN],*p_code;
	uint8_t m_buf[MESSAGE_LEN];

	//c2 decompress
	poly_decompress(c+DIM_N,c2,c2_len);
	//c1*sk
	//poly_mul(c,sk,out,c2_len);	
	//unsigned long long int t1 = hal_get_time();
	poly_mul(c,sk,out,c2_len);	
	//unsigned long long int t2 = hal_get_time();
	//sprintf(outs, "poly_mul(c,sk,out,c2_len):%llu\n",t2-t1 );
	//hal_send_str(outs);

	//no d2
	
	//compute mlen
	*mlen=c2_len/16-ECC_LEN;
	//shif the pointer of ecc data
	p_code=code+(DATA_LEN-(*mlen));


	//init code
	memset(code,0,CODE_LEN);


	//D2 decode
	d2_decode(c2,out,p_code,c2_len);

	//bch decode to recover m
	ecc_dec(m_buf,code);

	//get plaintext
	memcpy(m,m_buf+(DATA_LEN-(*mlen)),*mlen);
	
	return 0;
}

static int encode_to_e2(uint8_t *e2, const uint8_t *m, unsigned long long mlen, int *c2_len)
{
	int i;
	int8_t message;
	int vec_bound;
	uint8_t m_buf[MESSAGE_LEN];
	uint8_t code[CODE_LEN],*p_code;
	//BCH encoding
	//package m_buf
	memset(m_buf,0,MESSAGE_LEN);
	//set data
	memcpy(m_buf+(DATA_LEN-mlen),m,mlen);
	//encode m with ecc code
	ecc_enc(m_buf,code);
	//set p_code
	p_code=code+(DATA_LEN-mlen);
	
	//compute the length of c2
	*c2_len=(mlen+ECC_LEN)*8*2;
	vec_bound=*c2_len/2;
	//compute  code*q/2+e2, 
	for(i=0;i<vec_bound;i++)
	{
		//RATIO=q/2. add code*q/2 to e2
		message=RATIO*((p_code[i/8]>>(i%8))&1);
		e2[i]=e2[i]+message;
		//D2 encode, repeat at i+vec_bound
		e2[i+vec_bound]=e2[i+vec_bound]+message;
	}
	
	return 0;
}

// encryption with seed
int pke_enc_seed(const uint8_t *pk, const uint8_t *m, unsigned long long mlen, uint8_t *c, unsigned long long *clen, uint8_t *seed)
{
	uint8_t seeds[3*SEED_LEN];
	uint8_t r[DIM_N];
	uint8_t e1[DIM_N],e2[DIM_N];
	uint8_t c2[C2_VEC_NUM];
	uint8_t a[DIM_N];
	//char outs[128];

	int c2_len;

	if(mlen>MESSAGE_LEN)
	{
		return -1;
	}
	
	//generate  a from seed in the first part of pk
	gen_a(a,pk);
	//generate three seeds for r,e1,e2
	pseudo_random_bytes(seeds,3*SEED_LEN,seed);
	//generate random vector r
	gen_r(r,seeds);
	//generate error vector e1
	gen_e(e1,seeds+SEED_LEN);
	//generate c1: c1=a*r+e1
	//unsigned long long int t1 = hal_get_time();
	poly_aff(a,r,e1,c,DIM_N);
	//unsigned long long int t2 = hal_get_time();
	//sprintf(outs, "poly_aff(a,r,e1,c,DIM_N):%llu\n",t2-t1 );
	//hal_send_str(outs);
	//generate error vector e2
	gen_e(e2,seeds+2*SEED_LEN);
	//encode message to e2
	encode_to_e2(e2,m,mlen, &c2_len);
	//c2=b*r+e2+m*[q/2]
	//t1 = hal_get_time();
	poly_aff(pk+SEED_LEN,r,e2,c2,c2_len);
	//t2 = hal_get_time();
	//sprintf(outs, "poly_aff(pk+SEED_LEN,r,e2,c2,c2_len):%llu\n",t2-t1 );
	//hal_send_str(outs);
	//compress c2
	poly_compress(c2,c+DIM_N,c2_len);
	*clen=DIM_N+c2_len/2;

	return 0;

}

