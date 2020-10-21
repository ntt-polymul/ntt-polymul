#include "SABER_params.h"

#define N_SB (SABER_N >> 2)
#define N_SB_RES (2*N_SB-1)


#define N_SB_16 (N_SB >> 2)
#define N_SB_16_RES (2*N_SB_16-1)

#define TC_TC 1
#define TC_KARA 2

//#define MUL_TYPE TC_TC
#define MUL_TYPE TC_KARA

#if MUL_TYPE == TC_TC
	#define NUM_POLY_MID 7
#elif MUL_TYPE == TC_KARA
	#define NUM_POLY_MID 9
#endif


void pol_mul(uint16_t* a, uint16_t* b, uint16_t* res, uint16_t p, uint32_t n);

void pol_mul_sb(int16_t* a, int16_t* b, int16_t* res, uint16_t p, uint32_t n,uint32_t start);

void toom_cook_4way(const uint16_t* a1, const uint16_t* b1, uint16_t* result);

void pol_mul(uint16_t* a, uint16_t* b, uint16_t* res, uint16_t p, uint32_t n);

#if MUL_TYPE == TC_TC
	void evaluation_single(const uint16_t *b, uint32_t bw_ar[7][NUM_POLY_MID][N_SB_16]);

	void TC_evaluation_64_unrolled(const uint16_t* a1, uint32_t bw_ar[7][NUM_POLY_MID][N_SB_16], uint32_t w_ar[7][NUM_POLY_MID][N_SB_16_RES]);

	void TC_interpol_64_unrolled(uint32_t w_ar[7][NUM_POLY_MID][N_SB_16_RES], uint16_t *result);

#elif MUL_TYPE == TC_KARA

	void evaluation_single_kara(const uint16_t *b, uint16_t bw_ar[7][9][N_SB_16]);
	void TC_evaluation_unrolled_kara(const uint16_t* a1, uint16_t bw_ar[7][9][N_SB_16], uint16_t w_ar[7][NUM_POLY_MID][N_SB_16_RES]);

	void TC_interpol1_kara(uint16_t w_ar[7][NUM_POLY_MID][N_SB_16_RES], uint16_t *result);

#endif
