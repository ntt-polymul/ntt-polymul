#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "poly_mul.h"

void TC_evaluation_16(const uint32_t* a1, uint32_t bw_ar[NUM_POLY_MID][N_SB_16], uint32_t *w1,  uint32_t *w2,  uint32_t *w3,  uint32_t *w4,  uint32_t *w5,  uint32_t *w6,  uint32_t *w7);

void TC_interpol_16(uint32_t *w1, uint32_t *w2, uint32_t *w3, uint32_t *w4, uint32_t *w5, uint32_t *w6, uint32_t *w7, uint32_t *result);

static inline int16_t reduce(int16_t a, int64_t p);

#if MUL_TYPE == TC_TC
//----------------------------------routines for lazy interpolation TC-TC--------------------------------

/*
void pol_mul_noreduce_32(const uint32_t* a, const uint32_t* b, uint32_t* res, uint32_t n) //simple school book without polynomial reduction. For 32 bits 
{ 
	uint32_t i, j;
	
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			res[i+j]=res[i+j] + (a[i] * b[j]);
		}
	}
}
*/

void pol_mul_noreduce_32(const uint32_t* a, const uint32_t* b, uint32_t* res, uint32_t n) //simple school book without polynomial reduction. For 32 bits 
{ 
	uint32_t i, j;
	uint32_t a0, a1;
	
	for (i = 0; i < n/2; i++) {
		a0=a[2*i];
		a1=a[2*i+1];
		for (j = 0; j < n; j++) {
			res[2*i+j]=res[2*i+j] + (a0 * b[j]);
			res[2*i+j+1]=res[2*i+j+1] + ( a1 * b[j]);
		}
	}
}


void evaluation_single(const uint16_t *b, uint32_t bw_ar[7][NUM_POLY_MID][N_SB_16]){ // for precomputing B

	//printf("\nEvaluation Single TC-TC\n");
	uint32_t r0, r1, r2, r3, r4, r5, r6, r7;

	uint32_t bw1[N_SB], bw2[N_SB], bw3[N_SB], bw4[N_SB], bw5[N_SB], bw6[N_SB], bw7[N_SB]; //can be reduced

	uint16_t *B0, *B1, *B2, *B3;

	int j;

	B0 = (uint16_t*)b;
	B1 = (uint16_t*)&b[N_SB];
	B2 = (uint16_t*)&b[2*N_SB];
	B3 = (uint16_t*)&b[3*N_SB];

	for (j = 0; j < N_SB; ++j) {
		r0 = B0[j];
		r1 = B1[j];
		r2 = B2[j];
		r3 = B3[j];
		r4 = r0 + r2;
		r5 = r1 + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		bw3[j] = r6;
		bw4[j] = r7;
		r4 = ((r0 << 2)+r2) << 1;
		r5 = (r1 << 2) + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		bw5[j] = r6;
		bw6[j] = r7;
		r4 = (r3 << 3) + (r2 << 2) + (r1 << 1) + r0;
		bw2[j] = r4; 
		bw7[j] = r0;
		bw1[j] = r3;
	}

	for (j = 0; j < N_SB_16; ++j) {
		r0 = bw1[j];
		r1 = bw1[j+N_SB_16];
		r2 = bw1[j+2*N_SB_16];
		r3 = bw1[j+3*N_SB_16];
		r4 = r0 + r2;
		r5 = r1 + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		bw_ar[0][2][j] = r6;
		bw_ar[0][3][j] = r7;
		r4 = ((r0 << 2)+r2) << 1;
		r5 = (r1 << 2) + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		bw_ar[0][4][j] = r6;
		bw_ar[0][5][j] = r7;
		r4 = (r3 << 3) + (r2 << 2) + (r1 << 1) + r0;
		bw_ar[0][1][j] = r4; 
		bw_ar[0][6][j] = r0;
		bw_ar[0][0][j] = r3;
	}


	for (j = 0; j < N_SB_16; ++j) {
		r0 = bw2[j];
		r1 = bw2[j+N_SB_16];
		r2 = bw2[j+2*N_SB_16];
		r3 = bw2[j+3*N_SB_16];
		r4 = r0 + r2;
		r5 = r1 + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		bw_ar[1][2][j] = r6;
		bw_ar[1][3][j] = r7;
		r4 = ((r0 << 2)+r2) << 1;
		r5 = (r1 << 2) + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		bw_ar[1][4][j] = r6;
		bw_ar[1][5][j] = r7;
		r4 = (r3 << 3) + (r2 << 2) + (r1 << 1) + r0;
		bw_ar[1][1][j] = r4; 
		bw_ar[1][6][j] = r0;
		bw_ar[1][0][j] = r3;
	}

	for (j = 0; j < N_SB_16; ++j) {
		r0 = bw3[j];
		r1 = bw3[j+N_SB_16];
		r2 = bw3[j+2*N_SB_16];
		r3 = bw3[j+3*N_SB_16];
		r4 = r0 + r2;
		r5 = r1 + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		bw_ar[2][2][j] = r6;
		bw_ar[2][3][j] = r7;
		r4 = ((r0 << 2)+r2) << 1;
		r5 = (r1 << 2) + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		bw_ar[2][4][j] = r6;
		bw_ar[2][5][j] = r7;
		r4 = (r3 << 3) + (r2 << 2) + (r1 << 1) + r0;
		bw_ar[2][1][j] = r4; 
		bw_ar[2][6][j] = r0;
		bw_ar[2][0][j] = r3;
	}


	for (j = 0; j < N_SB_16; ++j) {
		r0 = bw4[j];
		r1 = bw4[j+N_SB_16];
		r2 = bw4[j+2*N_SB_16];
		r3 = bw4[j+3*N_SB_16];
		r4 = r0 + r2;
		r5 = r1 + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		bw_ar[3][2][j] = r6;
		bw_ar[3][3][j] = r7;
		r4 = ((r0 << 2)+r2) << 1;
		r5 = (r1 << 2) + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		bw_ar[3][4][j] = r6;
		bw_ar[3][5][j] = r7;
		r4 = (r3 << 3) + (r2 << 2) + (r1 << 1) + r0;
		bw_ar[3][1][j] = r4; 
		bw_ar[3][6][j] = r0;
		bw_ar[3][0][j] = r3;
	}

	for (j = 0; j < N_SB_16; ++j) {
		r0 = bw5[j];
		r1 = bw5[j+N_SB_16];
		r2 = bw5[j+2*N_SB_16];
		r3 = bw5[j+3*N_SB_16];
		r4 = r0 + r2;
		r5 = r1 + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		bw_ar[4][2][j] = r6;
		bw_ar[4][3][j] = r7;
		r4 = ((r0 << 2)+r2) << 1;
		r5 = (r1 << 2) + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		bw_ar[4][4][j] = r6;
		bw_ar[4][5][j] = r7;
		r4 = (r3 << 3) + (r2 << 2) + (r1 << 1) + r0;
		bw_ar[4][1][j] = r4; 
		bw_ar[4][6][j] = r0;
		bw_ar[4][0][j] = r3;
	}


	for (j = 0; j < N_SB_16; ++j) {
		r0 = bw6[j];
		r1 = bw6[j+N_SB_16];
		r2 = bw6[j+2*N_SB_16];
		r3 = bw6[j+3*N_SB_16];
		r4 = r0 + r2;
		r5 = r1 + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		bw_ar[5][2][j] = r6;
		bw_ar[5][3][j] = r7;
		r4 = ((r0 << 2)+r2) << 1;
		r5 = (r1 << 2) + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		bw_ar[5][4][j] = r6;
		bw_ar[5][5][j] = r7;
		r4 = (r3 << 3) + (r2 << 2) + (r1 << 1) + r0;
		bw_ar[5][1][j] = r4; 
		bw_ar[5][6][j] = r0;
		bw_ar[5][0][j] = r3;
	}

	for (j = 0; j < N_SB_16; ++j) {
		r0 = bw7[j];
		r1 = bw7[j+N_SB_16];
		r2 = bw7[j+2*N_SB_16];
		r3 = bw7[j+3*N_SB_16];
		r4 = r0 + r2;
		r5 = r1 + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		bw_ar[6][2][j] = r6;
		bw_ar[6][3][j] = r7;
		r4 = ((r0 << 2)+r2) << 1;
		r5 = (r1 << 2) + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		bw_ar[6][4][j] = r6;
		bw_ar[6][5][j] = r7;
		r4 = (r3 << 3) + (r2 << 2) + (r1 << 1) + r0;
		bw_ar[6][1][j] = r4; 
		bw_ar[6][6][j] = r0;
		bw_ar[6][0][j] = r3;
	}


}

void TC_evaluation_64_unrolled(const uint16_t* a1, uint32_t bw_ar[7][NUM_POLY_MID][N_SB_16], uint32_t w_ar[7][NUM_POLY_MID][N_SB_16_RES])//TC+TC unrolled
{

	//printf("\nTC 1st level evaluation TC-TC\n");
	uint32_t aw1[N_SB], aw2[N_SB], aw3[N_SB], aw4[N_SB], aw5[N_SB], aw6[N_SB], aw7[N_SB];

	uint32_t r0, r1, r2, r3, r4, r5, r6, r7;
	uint16_t *A0, *A1, *A2, *A3;

	A0 = (uint16_t*)a1;
	A1 = (uint16_t*)&a1[N_SB];
	A2 = (uint16_t*)&a1[2*N_SB];
	A3 = (uint16_t*)&a1[3*N_SB];

	int j;

// EVALUATION
	for (j = 0; j < N_SB; ++j) {
		r0 = A0[j];
		r1 = A1[j];
		r2 = A2[j];
		r3 = A3[j];
		r4 = r0 + r2;
		r5 = r1 + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		aw3[j] = r6;
		aw4[j] = r7;
		r4 = ((r0 << 2)+r2) << 1;
		r5 = (r1 << 2) + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		aw5[j] = r6;
		aw6[j] = r7;
		r4 = (r3 << 3) + (r2 << 2) + (r1 << 1) + r0;
		aw2[j] = r4; aw7[j] = r0;
		aw1[j] = r3;
	}

// MULTIPLICATION

	TC_evaluation_16(aw1, bw_ar[0], w_ar[0][0], w_ar[0][1], w_ar[0][2], w_ar[0][3], w_ar[0][4], w_ar[0][5], w_ar[0][6]);
	TC_evaluation_16(aw2, bw_ar[1], w_ar[1][0], w_ar[1][1], w_ar[1][2], w_ar[1][3], w_ar[1][4], w_ar[1][5], w_ar[1][6]);
	TC_evaluation_16(aw3, bw_ar[2], w_ar[2][0], w_ar[2][1], w_ar[2][2], w_ar[2][3], w_ar[2][4], w_ar[2][5], w_ar[2][6]);
	TC_evaluation_16(aw4, bw_ar[3], w_ar[3][0], w_ar[3][1], w_ar[3][2], w_ar[3][3], w_ar[3][4], w_ar[3][5], w_ar[3][6]);	
	TC_evaluation_16(aw5, bw_ar[4], w_ar[4][0], w_ar[4][1], w_ar[4][2], w_ar[4][3], w_ar[4][4], w_ar[4][5], w_ar[4][6]);
	TC_evaluation_16(aw6, bw_ar[5], w_ar[5][0], w_ar[5][1], w_ar[5][2], w_ar[5][3], w_ar[5][4], w_ar[5][5], w_ar[5][6]);
	TC_evaluation_16(aw7, bw_ar[6], w_ar[6][0], w_ar[6][1], w_ar[6][2], w_ar[6][3], w_ar[6][4], w_ar[6][5], w_ar[6][6]);

}

void TC_interpol_64_unrolled(uint32_t w_ar[7][NUM_POLY_MID][N_SB_16_RES], uint16_t *result){ //unrolled

	//printf("\nInterpolation called\n");

	//printf("\nTC 1st level interpol TC-TC\n");

	uint32_t r0, r1, r2, r3, r4, r5, r6;
	uint16_t inv3 = 43691, inv9 = 36409, inv15 = 61167;

	uint32_t w1[N_SB_RES] = {0}, w2[N_SB_RES] = {0}, w3[N_SB_RES] = {0}, w4[N_SB_RES] = {0}, w5[N_SB_RES] = {0}, w6[N_SB_RES] = {0}, w7[N_SB_RES] = {0};

	int i;
	
	uint16_t * C;
	C = result;

	TC_interpol_16(w_ar[0][0], w_ar[0][1], w_ar[0][2], w_ar[0][3], w_ar[0][4], w_ar[0][5], w_ar[0][6], w1);
	TC_interpol_16(w_ar[1][0], w_ar[1][1], w_ar[1][2], w_ar[1][3], w_ar[1][4], w_ar[1][5], w_ar[1][6], w2);
	TC_interpol_16(w_ar[2][0], w_ar[2][1], w_ar[2][2], w_ar[2][3], w_ar[2][4], w_ar[2][5], w_ar[2][6], w3);
	TC_interpol_16(w_ar[3][0], w_ar[3][1], w_ar[3][2], w_ar[3][3], w_ar[3][4], w_ar[3][5], w_ar[3][6], w4);
	TC_interpol_16(w_ar[4][0], w_ar[4][1], w_ar[4][2], w_ar[4][3], w_ar[4][4], w_ar[4][5], w_ar[4][6], w5);
	TC_interpol_16(w_ar[5][0], w_ar[5][1], w_ar[5][2], w_ar[5][3], w_ar[5][4], w_ar[5][5], w_ar[5][6], w6);
	TC_interpol_16(w_ar[6][0], w_ar[6][1], w_ar[6][2], w_ar[6][3], w_ar[6][4], w_ar[6][5], w_ar[6][6], w7);	

	for (i = 0; i < N_SB_RES; ++i) {
		r0 = w1[i];
		r1 = w2[i];
		r2 = w3[i];
		r3 = w4[i];
		r4 = w5[i];
		r5 = w6[i];
		r6 = w7[i];

		r1 = r1 + r4;
		r5 = r5 - r4;
		r3 = ((r3-r2) >> 1);
		r4 = r4 - r0;
		r4 = r4 - (r6 << 6);
		r4 = (r4 << 1) + r5;
		r2 = r2 + r3;
		r1 = r1 - (r2 << 6) - r2;
		r2 = r2 - r6;
		r2 = r2 - r0;
		r1 = r1 + 45*r2;
		r4 = (((r4 - (r2 << 3))*inv3) >> 3);
		r5 = r5 + r1;
		r1 = (((r1 + (r3 << 4))*inv9) >> 1);
		r3 = -(r3 + r1);
		r5 = (((30*r1 - r5)*inv15) >> 2);
		r2 = r2 - r4;
		r1 = r1 - r5;

		C[i]     += r6;
		C[i+64]  += r5;
		C[i+128] += r4;
		C[i+192] += r3;
		C[i+256] += r2;
		C[i+320] += r1;
		C[i+384] += r0;
	}



}


void TC_evaluation_16(const uint32_t* a1, uint32_t bw_ar[NUM_POLY_MID][N_SB_16], uint32_t *w1,  uint32_t *w2,  uint32_t *w3,  uint32_t *w4,  uint32_t *w5,  uint32_t *w6,  uint32_t *w7){

	//printf("\nTC 2nd level evaluation TC-TC\n");

	uint32_t aw1[N_SB_16], aw2[N_SB_16], aw3[N_SB_16], aw4[N_SB_16], aw5[N_SB_16], aw6[N_SB_16], aw7[N_SB_16];

	uint32_t r0, r1, r2, r3, r4, r5, r6, r7;
	uint32_t *A0, *A1, *A2, *A3;
	A0 = (uint32_t*)a1;
	A1 = (uint32_t*)&a1[N_SB_16];
	A2 = (uint32_t*)&a1[2*N_SB_16];
	A3 = (uint32_t*)&a1[3*N_SB_16];
	/*
	B0 = (uint32_t*)b1;
	B1 = (uint32_t*)&b1[N_SB_16];
	B2 = (uint32_t*)&b1[2*N_SB_16];
	B3 = (uint32_t*)&b1[3*N_SB_16];
	*/
	int j;

// EVALUATION
	for (j = 0; j < N_SB_16; ++j) {
		r0 = A0[j];
		r1 = A1[j];
		r2 = A2[j];
		r3 = A3[j];
		r4 = r0 + r2;
		r5 = r1 + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		aw3[j] = r6;
		aw4[j] = r7;
		r4 = ((r0 << 2)+r2) << 1;
		r5 = (r1 << 2) + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		aw5[j] = r6;
		aw6[j] = r7;
		r4 = (r3 << 3) + (r2 << 2) + (r1 << 1) + r0;
		aw2[j] = r4; aw7[j] = r0;
		aw1[j] = r3;
	}
// MULTIPLICATION

	pol_mul_noreduce_32(aw1,bw_ar[0],w1, 16);
	pol_mul_noreduce_32(aw2,bw_ar[1],w2, 16);
	pol_mul_noreduce_32(aw3,bw_ar[2],w3, 16);
	pol_mul_noreduce_32(aw4,bw_ar[3],w4, 16);
	pol_mul_noreduce_32(aw5,bw_ar[4],w5, 16);
	pol_mul_noreduce_32(aw6,bw_ar[5],w6, 16);
	pol_mul_noreduce_32(aw7,bw_ar[6],w7, 16);
}

void TC_interpol_16(uint32_t *w1, uint32_t *w2, uint32_t *w3, uint32_t *w4, uint32_t *w5, uint32_t *w6, uint32_t *w7, uint32_t *result){

	//printf("\nTC 2nd level interpol TC-TC\n");
	uint32_t r0, r1, r2, r3, r4, r5, r6;
	uint32_t inv3 = 174763, inv9 = 233017, inv15 = 454383;
	uint32_t * C;
	C = result;
	int i;

		for (i = 0; i < N_SB_16_RES; ++i) {
		r0 = w1[i];
		r1 = w2[i];
		r2 = w3[i];
		r3 = w4[i];
		r4 = w5[i];
		r5 = w6[i];
		r6 = w7[i];

		r1 = r1 + r4;
		r5 = r5 - r4;
		r3 = ((r3-r2) >> 1);
		r4 = r4 - r0;
		r4 = r4 - (r6 << 6);
		r4 = (r4 << 1) + r5;
		r2 = r2 + r3;
		r1 = r1 - (r2 << 6) - r2;
		r2 = r2 - r6;
		r2 = r2 - r0;
		r1 = r1 + 45*r2;
		r4 = (((r4 - (r2 << 3))*inv3) >> 3);
		r5 = r5 + r1;
		r1 = (((r1 + (r3 << 4))*inv9) >> 1);
		r3 = -(r3 + r1);
		r5 = (((30*r1 - r5)*inv15) >> 2);
		r2 = r2 - r4;
		r1 = r1 - r5;
		
		C[i]    +=r6;
		C[i+16] +=r5;
		C[i+32] +=r4;
		C[i+48] +=r3;
		C[i+64] +=r2;
		C[i+80] +=r1;
		C[i+96] +=r0;
		
	}



}



//----------------------------------routines for lazy interpolation TC-TC ends---------------------------
#elif MUL_TYPE == TC_KARA

// void pol_mul_noreduce(const int16_t* a, const int16_t* b, int16_t* res, uint16_t p, uint32_t n) //simple school book //without polynomial reduction
// { 
// 	uint32_t i, j;
// 	uint16_t a0, a1;	

// 	for (i = 0; i < n/2; i++) {
// 		a0=a[2*i];
// 		a1=a[2*i+1];
// 		for (j = 0; j < n; j++) {
// 			res[2*i+j]=res[2*i+j] + (a0 * b[j]);
// 			res[2*i+j+1]=res[2*i+j+1] + ( a1 * b[j]);

// 		}
// 	}
// }


//----------------------------------routines for lazy interpolation TC-KARA -----------------------------
// void evaluation_single_kara(const uint16_t *b, uint16_t bw_ar[7][9][N_SB_16]){ // for precomputing B

// 	uint16_t r0, r1, r2, r3, r4, r5, r6, r7;

// 	//uint16_t bw1[N_SB], bw2[N_SB], bw3[N_SB], bw4[N_SB], bw5[N_SB], bw6[N_SB], bw7[N_SB]; //can be reduced
// 	uint16_t bw[N_SB*5];
// 	uint16_t *bw1, *bw2, *bw3, *bw4, *bw5, *bw6, *bw7;
// 	bw1 = (uint16_t*)&b[3*N_SB];
// 	bw2 = (uint16_t*)bw;
// 	bw3 = (uint16_t*)&bw[N_SB];
// 	bw4 = (uint16_t*)&bw[2*N_SB];
// 	bw5 = (uint16_t*)&bw[3*N_SB];
// 	bw6 = (uint16_t*)&bw[4*N_SB];
// 	bw7 = (uint16_t*)b;

// 	//evaluation_single_kara_asm(b, bw_ar);

// 	// kara_eval_precomp(&b[3*N_SB],  bw_ar[0]);
// 	// kara_eval_precomp(bw,          bw_ar[1]);
// 	// kara_eval_precomp(&bw[N_SB],   bw_ar[2]);
// 	// kara_eval_precomp(&bw[2*N_SB], bw_ar[3]);
// 	// kara_eval_precomp(&bw[3*N_SB], bw_ar[4]);
// 	// kara_eval_precomp(&bw[4*N_SB], bw_ar[5]);
// 	// kara_eval_precomp(bw7,         bw_ar[6]);

// 	uint16_t *B0, *B1, *B2, *B3;

// 	int i;

// 	B0 = (uint16_t*)b;
// 	B1 = (uint16_t*)&b[N_SB];
// 	B2 = (uint16_t*)&b[2*N_SB];
// 	B3 = (uint16_t*)&b[3*N_SB];

// 	for (i = 0; i < N_SB; ++i) {
// 		r0 = B0[i];
// 		r1 = B1[i];
// 		r2 = B2[i];
// 		r3 = B3[i];
// 		r4 = r0 + r2;
// 		r5 = r1 + r3;
// 		r6 = r4 + r5; r7 = r4 - r5;
// 		bw3[i] = r6;
// 		bw4[i] = r7;
// 		r4 = ((r0 << 2)+r2) << 1;
// 		r5 = (r1 << 2) + r3;
// 		r6 = r4 + r5; r7 = r4 - r5;
// 		bw5[i] = r6;
// 		bw6[i] = r7;
// 		r4 = (r3 << 3) + (r2 << 2) + (r1 << 1) + r0;
// 		bw2[i] = r4; 
// 		bw7[i] = r0;
// 		bw1[i] = r3;
// 	}

// 	for(i=0;i<16;i++){
// 		bw_ar[0][0][i]=bw1[i]+bw1[i+16];
// 		bw_ar[0][1][i]=bw1[i+32]+bw1[i+16+32];
// 		bw_ar[0][2][i]=bw1[i]+bw1[i+32];
// 		bw_ar[0][3][i]=bw1[i+16]+bw1[i+32+16];
// 		bw_ar[0][4][i]=bw1[i]+bw1[i+16]+bw1[i+32]+bw1[i+32+16];
// 		bw_ar[0][5][i]=bw1[i];
// 		bw_ar[0][6][i]=bw1[i+16];
// 		bw_ar[0][7][i]=bw1[i+32];
// 		bw_ar[0][8][i]=bw1[i+32+16];
// 	}

// 	for(i=0;i<16;i++){
// 		bw_ar[1][0][i]=bw2[i]+bw2[i+16];
// 		bw_ar[1][1][i]=bw2[i+32]+bw2[i+16+32];
// 		bw_ar[1][2][i]=bw2[i]+bw2[i+32];
// 		bw_ar[1][3][i]=bw2[i+16]+bw2[i+32+16];
// 		bw_ar[1][4][i]=bw2[i]+bw2[i+16]+bw2[i+32]+bw2[i+32+16];
// 		bw_ar[1][5][i]=bw2[i];
// 		bw_ar[1][6][i]=bw2[i+16];
// 		bw_ar[1][7][i]=bw2[i+32];
// 		bw_ar[1][8][i]=bw2[i+32+16];
// 	}

// 	for(i=0;i<16;i++){
// 		bw_ar[2][0][i]=bw3[i]+bw3[i+16];
// 		bw_ar[2][1][i]=bw3[i+32]+bw3[i+16+32];
// 		bw_ar[2][2][i]=bw3[i]+bw3[i+32];
// 		bw_ar[2][3][i]=bw3[i+16]+bw3[i+32+16];
// 		bw_ar[2][4][i]=bw3[i]+bw3[i+16]+bw3[i+32]+bw3[i+32+16];
// 		bw_ar[2][5][i]=bw3[i];
// 		bw_ar[2][6][i]=bw3[i+16];
// 		bw_ar[2][7][i]=bw3[i+32];
// 		bw_ar[2][8][i]=bw3[i+32+16];
// 	}

// 	for(i=0;i<16;i++){
// 		bw_ar[3][0][i]=bw4[i]+bw4[i+16];
// 		bw_ar[3][1][i]=bw4[i+32]+bw4[i+16+32];
// 		bw_ar[3][2][i]=bw4[i]+bw4[i+32];
// 		bw_ar[3][3][i]=bw4[i+16]+bw4[i+32+16];
// 		bw_ar[3][4][i]=bw4[i]+bw4[i+16]+bw4[i+32]+bw4[i+32+16];
// 		bw_ar[3][5][i]=bw4[i];
// 		bw_ar[3][6][i]=bw4[i+16];
// 		bw_ar[3][7][i]=bw4[i+32];
// 		bw_ar[3][8][i]=bw4[i+32+16];
// 	}

// 	for(i=0;i<16;i++){
// 		bw_ar[4][0][i]=bw5[i]+bw5[i+16];
// 		bw_ar[4][1][i]=bw5[i+32]+bw5[i+16+32];
// 		bw_ar[4][2][i]=bw5[i]+bw5[i+32];
// 		bw_ar[4][3][i]=bw5[i+16]+bw5[i+32+16];
// 		bw_ar[4][4][i]=bw5[i]+bw5[i+16]+bw5[i+32]+bw5[i+32+16];
// 		bw_ar[4][5][i]=bw5[i];
// 		bw_ar[4][6][i]=bw5[i+16];
// 		bw_ar[4][7][i]=bw5[i+32];
// 		bw_ar[4][8][i]=bw5[i+32+16];
// 	}

// 	for(i=0;i<16;i++){
// 		bw_ar[5][0][i]=bw6[i]+bw6[i+16];
// 		bw_ar[5][1][i]=bw6[i+32]+bw6[i+16+32];
// 		bw_ar[5][2][i]=bw6[i]+bw6[i+32];
// 		bw_ar[5][3][i]=bw6[i+16]+bw6[i+32+16];
// 		bw_ar[5][4][i]=bw6[i]+bw6[i+16]+bw6[i+32]+bw6[i+32+16];
// 		bw_ar[5][5][i]=bw6[i];
// 		bw_ar[5][6][i]=bw6[i+16];
// 		bw_ar[5][7][i]=bw6[i+32];
// 		bw_ar[5][8][i]=bw6[i+32+16];
// 	}

// 	for(i=0;i<16;i++){
// 		bw_ar[6][0][i]=bw7[i]+bw7[i+16];
// 		bw_ar[6][1][i]=bw7[i+32]+bw7[i+16+32];
// 		bw_ar[6][2][i]=bw7[i]+bw7[i+32];
// 		bw_ar[6][3][i]=bw7[i+16]+bw7[i+32+16];
// 		bw_ar[6][4][i]=bw7[i]+bw7[i+16]+bw7[i+32]+bw7[i+32+16];
// 		bw_ar[6][5][i]=bw7[i];
// 		bw_ar[6][6][i]=bw7[i+16];
// 		bw_ar[6][7][i]=bw7[i+32];
// 		bw_ar[6][8][i]=bw7[i+32+16];
// 	}

// }

// void TC_evaluation_unrolled_kara(const uint16_t* a1, uint16_t bw_ar[7][9][N_SB_16], uint16_t w_ar[7][NUM_POLY_MID][N_SB_16_RES])
// {

// //	TC_evaluation_unrolled_kara_asm(a1, bw_ar, w_ar);

// 	//printf("\nEvaluation called\n");
// 	//uint16_t aw1[N_SB], aw2[N_SB], aw3[N_SB], aw4[N_SB], aw5[N_SB], aw6[N_SB], aw7[N_SB];
// 	// uint16_t bw1[N_SB], bw2[N_SB], bw3[N_SB], bw4[N_SB], bw5[N_SB], bw6[N_SB], bw7[N_SB];

// 	uint16_t r0, r1, r2, r3, r4, r5, r6, r7;
// 	uint16_t *A0, *A1, *A2, *A3, *B0, *B1, *B2, *B3;
// 	A0 = (uint16_t*)a1;
// 	A1 = (uint16_t*)&a1[N_SB];
// 	A2 = (uint16_t*)&a1[2*N_SB];
// 	A3 = (uint16_t*)&a1[3*N_SB];
// 	// B0 = (uint16_t*)b1;
// 	// B1 = (uint16_t*)&b1[N_SB];
// 	// B2 = (uint16_t*)&b1[2*N_SB];
// 	// B3 = (uint16_t*)&b1[3*N_SB];

// 	//uint16_t * C;
// 	//C = result;

// 	int j;

// // EVALUATION
// 	for (j = 0; j < N_SB; ++j) {
// 		r0 = A0[j];
// 		r1 = A1[j];
// 		r2 = A2[j];
// 		r3 = A3[j];
// 		r4 = r0 + r2;
// 		r5 = r1 + r3;
// 		r6 = r4 + r5; r7 = r4 - r5;
// 		aw3[j] = r6;
// 		aw4[j] = r7;
// 		r4 = ((r0 << 2)+r2) << 1;
// 		r5 = (r1 << 2) + r3;
// 		r6 = r4 + r5; r7 = r4 - r5;
// 		aw5[j] = r6;
// 		aw6[j] = r7;
// 		r4 = (r3 << 3) + (r2 << 2) + (r1 << 1) + r0;
// 		aw2[j] = r4; //aw7[j] = r0;
// 		//aw1[j] = r3;
// 	}

	
// 	// for (j = 0; j < N_SB; ++j) { //this is needed as we are passing the bw too. Can be put in evaluation single.(TODO)
// 	// 	r0 = B0[j];
// 	// 	r1 = B1[j];
// 	// 	r2 = B2[j];
// 	// 	r3 = B3[j];
// 	// 	r4 = r0 + r2;
// 	// 	r5 = r1 + r3;
// 	// 	r6 = r4 + r5; r7 = r4 - r5;
// 	// 	bw3[j] = r6;
// 	// 	bw4[j] = r7;
// 	// 	r4 = ((r0 << 2)+r2) << 1;
// 	// 	r5 = (r1 << 2) + r3;
// 	// 	r6 = r4 + r5; r7 = r4 - r5;
// 	// 	bw5[j] = r6;
// 	// 	bw6[j] = r7;
// 	// 	r4 = (r3 << 3) + (r2 << 2) + (r1 << 1) + r0;
// 	// 	bw2[j] = r4; bw7[j] = r0;
// 	// 	bw1[j] = r3;
// 	// }
	


// // MULTIPLICATION

// 	eval_kara(&a1[3*N_SB], bw_ar[0], w_ar[0]);
// 	eval_kara(aw2, bw_ar[1], w_ar[1]);
// 	eval_kara(aw3, bw_ar[2], w_ar[2]);
// 	eval_kara(aw4, bw_ar[3], w_ar[3]);
// 	eval_kara(aw5, bw_ar[4], w_ar[4]);
// 	eval_kara(aw6, bw_ar[5], w_ar[5]);
// 	eval_kara(a1, bw_ar[6], w_ar[6]);
// }



// void TC_interpol1_kara(uint16_t w_ar[7][NUM_POLY_MID][N_SB_16_RES], uint16_t *result){

// 	//printf("\nInterpolation called\n");

// 	// uint16_t w1[N_SB_RES], w2[N_SB_RES], w3[N_SB_RES], w4[N_SB_RES], w5[N_SB_RES], w6[N_SB_RES], w7[N_SB_RES]; //can be reduced. use the result directly
// 	// uint16_t w[7*(N_SB_RES+1)];
// 	// uint16_t *w1, *w2, *w3, *w4, *w5, *w6, *w7;
// 	// w1 = (uint16_t *)w;
// 	// w2 = (uint16_t *)&w[N_SB_RES+1];
// 	// w3 = (uint16_t *)&w[2*(N_SB_RES+1)];
// 	// w4 = (uint16_t *)&w[3*(N_SB_RES+1)];
// 	// w5 = (uint16_t *)&w[4*(N_SB_RES+1)];
// 	// w6 = (uint16_t *)&w[5*(N_SB_RES+1)];
// 	// w7 = (uint16_t *)&w[6*(N_SB_RES+1)];

// 	// w1[31]=w1[63]=w1[95]=0; //three places are enough. others are over written in interpol_kara
// 	// w2[31]=w2[63]=w2[95]=0; 
// 	// w3[31]=w3[63]=w3[95]=0; 
// 	// w4[31]=w4[63]=w4[95]=0; 
// 	// w5[31]=w5[63]=w5[95]=0; 
// 	// w6[31]=w6[63]=w6[95]=0; 
// 	// w7[31]=w7[63]=w7[95]=0;
// 	// w1[127]=w2[127]=w3[127]=w4[127]=w5[127]=w6[127]=w7[127]=0;//to simplify asm

// 	// uint16_t r2, r3, r4, r5, r6, r7, r8, r9;
// 	// uint16_t inv3 = 43691, inv9 = 36409, inv15 = 61167;

// 	// int i;

// 	// uint16_t * C;
// 	// C = result;

// 	// interpol_kara(w_ar[0][0], w_ar[0][1], w_ar[0][2], w_ar[0][3], w_ar[0][4], w_ar[0][5], w_ar[0][6], w_ar[0][7], w_ar[0][8], w1);
// 	// interpol_kara(w_ar[1][0], w_ar[1][1], w_ar[1][2], w_ar[1][3], w_ar[1][4], w_ar[1][5], w_ar[1][6], w_ar[1][7], w_ar[1][8], w2);
// 	// interpol_kara(w_ar[2][0], w_ar[2][1], w_ar[2][2], w_ar[2][3], w_ar[2][4], w_ar[2][5], w_ar[2][6], w_ar[2][7], w_ar[2][8], w3);
// 	// interpol_kara(w_ar[3][0], w_ar[3][1], w_ar[3][2], w_ar[3][3], w_ar[3][4], w_ar[3][5], w_ar[3][6], w_ar[3][7], w_ar[3][8], w4);
// 	// interpol_kara(w_ar[4][0], w_ar[4][1], w_ar[4][2], w_ar[4][3], w_ar[4][4], w_ar[4][5], w_ar[4][6], w_ar[4][7], w_ar[4][8], w5);
// 	// interpol_kara(w_ar[5][0], w_ar[5][1], w_ar[5][2], w_ar[5][3], w_ar[5][4], w_ar[5][5], w_ar[5][6], w_ar[5][7], w_ar[5][8], w6);
// 	// interpol_kara(w_ar[6][0], w_ar[6][1], w_ar[6][2], w_ar[6][3], w_ar[6][4], w_ar[6][5], w_ar[6][6], w_ar[6][7], w_ar[6][8], w7);

// 	// interpol_kara(w_ar[0][0], w1);
// 	// interpol_kara(w_ar[1][0], w2);
// 	// interpol_kara(w_ar[2][0], w3);
// 	// interpol_kara(w_ar[3][0], w4);
// 	// interpol_kara(w_ar[4][0], w5);
// 	// interpol_kara(w_ar[5][0], w6);
// 	// interpol_kara(w_ar[6][0], w7);

// 	TC_interpol1_kara_asm(w_ar[0][0], result);

// 	// for (i = 0; i < N_SB_RES+1; ++i) {
// 	// 	TC_interpol1_kara_asm_round(&w1[i], &result[i], i);
// 	// 	r2 = w1[i];
// 	// 	r3 = w2[i];
// 	// 	r4 = w3[i];
// 	// 	r5 = w4[i];
// 	// 	r6 = w5[i];
// 	// 	r7 = w6[i];
// 	// 	r8 = w7[i];

// 	// 	r3 = r3 + r6;
// 	// 	r7 = r7 - r6;
// 	// 	r5 = ((r5-r4) >> 1);
// 	// 	r6 = r6 - r2;
// 	// 	r6 = r6 - (r8 << 6);
// 	// 	r6 = (r6 << 1) + r7;
// 	// 	r4 = r4 + r5;
// 	// 	r3 = r3 - (r4 << 6) - r4;
// 	// 	r4 = r4 - r8;
// 	// 	r4 = r4 - r2;
// 	// 	r3 = r3 + 45*r4;
// 	// 	r6 = (((r6 - (r4 << 3))*inv3) >> 3);
// 	// 	r7 = r7 + r3;
// 	// 	r3 = (((r3 + (r5 << 4))*inv9) >> 1);
// 	// 	r5 = -(r5 + r3);
// 	// 	r7 = (((30*r3 - r7)*inv15) >> 2);
// 	// 	r4 = r4 - r6;
// 	// 	r3 = r3 - r7;

// 	// 	// C[i]     += r8;
// 	// 	// C[i+64]  += r7;
// 	// 	// C[i+128] += r6;
// 	// 	// C[i+192] += r5;
// 	// 	// C[i+256] += r4;
// 	// 	// C[i+320] += r3;
// 	// 	// C[i+384] += r2;
// 	// 	//C[i]     = (C[i]    + r8 - r4) & (SABER_Q-1);
// 	// 	//C[i+64]  = (C[i+64] + r7 - r3) & (SABER_Q-1);
// 	// 	//C[i+128] = (C[i+128]+ r6 - r2) & (SABER_Q-1);

// 	// 	call_nop();

// 	// 	if (i < 64) {
// 	// 		//C[i+192] = (C[i+192] + r5) & (SABER_Q-1);
// 	// 	} else {
// 	// 		//C[i-64]  = (C[i-64]  - r5) & (SABER_Q-1);
// 	// 	}
// 	// 	call_nop();
// 	// }
// // */

// }

// void eval_kara(const uint16_t *a, const uint16_t *b, uint16_t bw_ar[9][N_SB_16], uint16_t res[NUM_POLY_MID][N_SB_16_RES]){


// 	//eval_kara_asm(a, b, bw_ar, res);

// 	int16_t i;

// 	//uint16_t a0[16], a1[16], a2[16], a3[16], a4[16];
// 	//uint16_t b0[16], b1[16], b2[16], b3[16], b4[16];

// 	uint16_t aw[80];
// 	uint16_t *a0, *a1, *a2, *a3, *a4;
// 	a0 = (uint16_t*)&aw[0];
// 	a1 = (uint16_t*)&aw[16];
// 	a2 = (uint16_t*)&aw[32];
// 	a3 = (uint16_t*)&aw[48];
// 	a4 = (uint16_t*)&aw[64];

// 	for(i=0;i<16;i++) {
// 		a0[i]=a[i]+a[i+16];
// 		a1[i]=a[i+32]+a[i+16+32];
// 		a2[i]=a[i]+a[i+32];
// 		a3[i]=a[i+16]+a[i+32+16];
// 		a4[i]=a[i]+a[i+16]+a[i+32]+a[i+32+16];
// 	}

// //----------------------------------------------------------------------

//     school_book_mul2_16_fast(a,b,res[3]);
// 	school_book_mul2_16_fast(a+16,b+16,res[4]);
// 	school_book_mul2_16_fast(a0,bw_ar[0],res[0]);
// 	//pol_mul_noreduce(a,b,res[3], 8192,16);
// 	//pol_mul_noreduce(a+16,b+16,res[4], 8192,16);
// 	//pol_mul_noreduce(a0,bw_ar[0],res[0], 8192,16);

// 	school_book_mul2_16_fast(a+32,b+32,res[5]);
// 	school_book_mul2_16_fast(a+48,b+48,res[6]);
// 	school_book_mul2_16_fast(a1,bw_ar[1],res[1]);	
// 	school_book_mul2_16_fast(a2,bw_ar[2],res[7]);	
// 	school_book_mul2_16_fast(a3,bw_ar[3],res[8]);
// 	school_book_mul2_16_fast(a4,bw_ar[4],res[2]);
// 	//pol_mul_noreduce(a+32,b+32,res[5], 8192,16);
// 	//pol_mul_noreduce(a+48,b+48,res[6], 8192,16);
// 	//pol_mul_noreduce(a1,bw_ar[1],res[1], 8192,16);	
// 	//pol_mul_noreduce(a2,bw_ar[2],res[7], 8192,16);	
// 	//pol_mul_noreduce(a3,bw_ar[3],res[8], 8192,16);
// 	//pol_mul_noreduce(a4,bw_ar[4],res[2], 8192,16);
// 	//printf("Poly mul called\n");


// }

// //void interpol_kara(uint16_t *res0, uint16_t *res1, uint16_t *res2, uint16_t *res3, uint16_t *res4, uint16_t *res5, uint16_t *res6, uint16_t *res7, uint16_t *res8, uint16_t *result){
// void interpol_kara(uint16_t *res_ptr, uint16_t *result){

// 	// uint16_t *res0, *res1, *res2, *res3, *res4, *res5, *res6, *res7, *res8;
// 	// res0 = res_ptr;
// 	// res1 = &res_ptr[31];
// 	// res2 = &res_ptr[2*31];
// 	// res3 = &res_ptr[3*31];
// 	// res4 = &res_ptr[4*31];
// 	// res5 = &res_ptr[5*31];
// 	// res6 = &res_ptr[6*31];
// 	// res7 = &res_ptr[7*31];
// 	// res8 = &res_ptr[8*31];


// 	// uint32_t i;
// 	uint16_t mid_result[63];

// 	mid_result[31]=0; //only one place is enough since other places are overwritten

// 	interpol_kara_asm(res_ptr, result, mid_result);

// 	//loop2
// 	// for(i=0;i<31;i++){ //first subtraction
// 	// 	res0[i]=res0[i]-res3[i]-res4[i];
// 	// 	res1[i]=res1[i]-res5[i]-res6[i];
// 	// 	res2[i]=res2[i]-res7[i]-res8[i];
// 	// // }

// 	// // for(i=0;i<31;i++){// copying to result directory to save sapce
// 	// 	result[i]=res3[i]; 
// 	// 	result[i+32]=res4[i];
// 	// 	result[i+64]=res5[i];
// 	// 	result[i+96]=res6[i];
// 	// 	mid_result[i]=res7[i];
// 	// 	mid_result[i+32]=res8[i];
// 	// }
// 	//loop3
// 	// for(i=0;i<31;i++){
// 	// 	result[i+16]=result[i+16]+res0[i];
// 	// 	result[i+16+64]=result[i+16+64]+res1[i];
// 	// 	mid_result[i+16]=mid_result[i+16]+res2[i];
// 	// }

// 	//------------------------------------------------------------------------

// 	// for(i=0;i<63;i++){ // first subtraction in the second loop
// 	// 	mid_result[i]=(mid_result[i]-result[i]-result[i+64]);
// 	// }

// 	// for(i=0;i<63;i++){
// 	// 	result[i+32]=(result[i+32]+mid_result[i]);
// 	// }
// }

//----------------------------------routines for lazy interpolation TC-Kara ends---------------------------

#endif
