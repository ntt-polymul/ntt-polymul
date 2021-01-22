#include "bch.h"
#include "ecc.h"
#include <string.h>
#include <stdlib.h>

//error corretion encode
int ecc_enc(const uint8_t *d, uint8_t *c)
{
	uint8_t ecc[ECC_LEN];
	
	//encoode
	encode_bch(d,DATA_LEN,ecc);
	//copy data to the first part of code
	memcpy(c,d,DATA_LEN);
	// compy ecc to the second part of code
	memcpy(c+DATA_LEN,ecc,ECC_LEN);
	
	return 0;
}

//error corrction decode
int ecc_dec(uint8_t *d, const uint8_t *c)
{
	int error_num;
	uint8_t ecc[ECC_LEN];
	
	//copy data from c to_d
	memcpy(d,c,DATA_LEN);
	memcpy(ecc,c+DATA_LEN,ECC_LEN);
	//decode
	error_num=decode_bch(d,DATA_LEN,ecc);
	
	return error_num;
}