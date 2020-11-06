#include "consts1024.h"
.include "shuffle.inc"
.include "fq.inc"

.macro update rln,rl0,rl1,rl2,rl3,rh0,rh1,rh2,rh3
vpaddw		%ymm\rh0,%ymm\rl0,%ymm\rln
vpsubw		%ymm\rh0,%ymm\rl0,%ymm\rh0

vpaddw		%ymm\rh1,%ymm\rl1,%ymm\rl0
vpsubw		%ymm\rh1,%ymm\rl1,%ymm\rh1

vpaddw		%ymm\rh2,%ymm\rl2,%ymm\rl1
vpsubw		%ymm\rh2,%ymm\rl2,%ymm\rh2

vpaddw		%ymm\rh3,%ymm\rl3,%ymm\rl2
vpsubw		%ymm\rh3,%ymm\rl3,%ymm\rh3
.endm

.macro levels0t2 off
/* level 0 */
vpbroadcastw	(_ZETAS_PINV+1)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+1)*2(%rdx),%ymm2
vmovdqa		(512+16*\off)*2(%rsi),%ymm8
vmovdqa		(640+16*\off)*2(%rsi),%ymm9
vmovdqa		(768+16*\off)*2(%rsi),%ymm10
vmovdqa		(896+16*\off)*2(%rsi),%ymm11

fqmulprecomp	1,2,8,x=3
fqmulprecomp	1,2,9,x=3
fqmulprecomp	1,2,10,x=3
fqmulprecomp	1,2,11,x=3

vmovdqa		(  0+16*\off)*2(%rsi),%ymm4
vmovdqa		(128+16*\off)*2(%rsi),%ymm5
vmovdqa		(256+16*\off)*2(%rsi),%ymm6
vmovdqa		(384+16*\off)*2(%rsi),%ymm7

update		3,4,5,6,7,8,9,10,11

/* level 1 */
vpbroadcastw	(_ZETAS_PINV+2)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+2)*2(%rdx),%ymm2
fqmulprecomp	1,2,5,x=7
fqmulprecomp	1,2,6,x=7
vpbroadcastw	(_ZETAS_PINV+3)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+3)*2(%rdx),%ymm2
fqmulprecomp	1,2,10,x=7
fqmulprecomp	1,2,11,x=7

update		7,3,4,8,9,5,6,10,11

/* level 2 */
vpbroadcastw	(_ZETAS_PINV+4)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+4)*2(%rdx),%ymm2
fqmulprecomp	1,2,3,x=9
vpbroadcastw	(_ZETAS_PINV+5)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+5)*2(%rdx),%ymm2
fqmulprecomp	1,2,6,x=9
vpbroadcastw	(_ZETAS_PINV+6)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+6)*2(%rdx),%ymm2
fqmulprecomp	1,2,8,x=9
vpbroadcastw	(_ZETAS_PINV+7)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+7)*2(%rdx),%ymm2
fqmulprecomp	1,2,11,x=9

update		9,7,5,4,10,3,6,8,11

vmovdqa		%ymm9,(  0+16*\off)*2(%rdi)
vmovdqa		%ymm3,(128+16*\off)*2(%rdi)
vmovdqa		%ymm7,(256+16*\off)*2(%rdi)
vmovdqa		%ymm6,(384+16*\off)*2(%rdi)
vmovdqa		%ymm5,(512+16*\off)*2(%rdi)
vmovdqa		%ymm8,(640+16*\off)*2(%rdi)
vmovdqa		%ymm4,(768+16*\off)*2(%rdi)
vmovdqa		%ymm11,(896+16*\off)*2(%rdi)
.endm

.macro levels3t7 off
vmovdqa		(128*\off+  0)*2(%rdi),%ymm4
vmovdqa		(128*\off+ 16)*2(%rdi),%ymm5
vmovdqa		(128*\off+ 32)*2(%rdi),%ymm6
vmovdqa		(128*\off+ 48)*2(%rdi),%ymm7
vmovdqa		(128*\off+ 64)*2(%rdi),%ymm8
vmovdqa		(128*\off+ 80)*2(%rdi),%ymm9
vmovdqa		(128*\off+ 96)*2(%rdi),%ymm10
vmovdqa		(128*\off+112)*2(%rdi),%ymm11

/* level 3 */
vpbroadcastw	(_ZETAS_PINV+8+\off)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+8+\off)*2(%rdx),%ymm2
fqmulprecomp	1,2,8,x=3
fqmulprecomp	1,2,9,x=3
fqmulprecomp	1,2,10,x=3
fqmulprecomp	1,2,11,x=3

update		3,4,5,6,7,8,9,10,11

/* level 4 */
vpbroadcastw	(_ZETAS_PINV+16+2*\off)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+16+2*\off)*2(%rdx),%ymm2
fqmulprecomp	1,2,5,x=7
fqmulprecomp	1,2,6,x=7
vpbroadcastw	(_ZETAS_PINV+16+2*\off+1)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+16+2*\off+1)*2(%rdx),%ymm2
fqmulprecomp	1,2,10,x=7
fqmulprecomp	1,2,11,x=7

update		7,3,4,8,9,5,6,10,11

/* level 5 */
vpbroadcastw	(_ZETAS_PINV+32+4*\off)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+32+4*\off)*2(%rdx),%ymm2
fqmulprecomp	1,2,3,x=9
vpbroadcastw	(_ZETAS_PINV+32+4*\off+1)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+32+4*\off+1)*2(%rdx),%ymm2
fqmulprecomp	1,2,6,x=9
vpbroadcastw	(_ZETAS_PINV+32+4*\off+2)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+32+4*\off+2)*2(%rdx),%ymm2
fqmulprecomp	1,2,8,x=9
vpbroadcastw	(_ZETAS_PINV+32+4*\off+3)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+32+4*\off+3)*2(%rdx),%ymm2
fqmulprecomp	1,2,11,x=9
vmovdqa		(_16XMONT_PINV)*2(%rdx),%ymm1
vmovdqa		(_16XMONT)*2(%rdx),%ymm2
fqmulprecomp	1,2,7,x=9 // extra reduction
fqmulprecomp	1,2,5,x=9 // extra reduction
fqmulprecomp	1,2,4,x=9 // extra reduction
fqmulprecomp	1,2,10,x=9 // extra reduction

update		9,7,5,4,10,3,6,8,11

/* level 6 */
vmovdqa		cdecl(idxdata)+_PACKWIDX(%rip),%ymm1
vpshufb		%ymm1,%ymm9,%ymm9
vpermq		$0xD8,%ymm9,%ymm9
vpshufb		%ymm1,%ymm3,%ymm3
vpermq		$0xD8,%ymm3,%ymm3
vpshufb		%ymm1,%ymm7,%ymm7
vpermq		$0xD8,%ymm7,%ymm7
vpshufb		%ymm1,%ymm6,%ymm6
vpermq		$0xD8,%ymm6,%ymm6
vpshufb		%ymm1,%ymm5,%ymm5
vpermq		$0xD8,%ymm5,%ymm5
vpshufb		%ymm1,%ymm8,%ymm8
vpermq		$0xD8,%ymm8,%ymm8
vpshufb		%ymm1,%ymm4,%ymm4
vpermq		$0xD8,%ymm4,%ymm4
vpshufb		%ymm1,%ymm11,%ymm11
vpermq		$0xD8,%ymm11,%ymm11

unpack		9,3,10,3
unpack		7,6,9,6
unpack		5,8,7,8
unpack		4,11,5,11

vpbroadcastd	(_ZETAS_PINV+64+8*\off)*2(%rdx),%ymm1
vpbroadcastd	(_ZETAS+64+8*\off)*2(%rdx),%ymm2
fqmulprecomp	1,2,3,x=4
vpbroadcastd	(_ZETAS_PINV+64+8*\off+2)*2(%rdx),%ymm1
vpbroadcastd	(_ZETAS+64+8*\off+2)*2(%rdx),%ymm2
fqmulprecomp	1,2,6,x=4
vpbroadcastd	(_ZETAS_PINV+64+8*\off+4)*2(%rdx),%ymm1
vpbroadcastd	(_ZETAS+64+8*\off+4)*2(%rdx),%ymm2
fqmulprecomp	1,2,8,x=4
vpbroadcastd	(_ZETAS_PINV+64+8*\off+6)*2(%rdx),%ymm1
vpbroadcastd	(_ZETAS+64+8*\off+6)*2(%rdx),%ymm2
fqmulprecomp	1,2,11,x=4

update		4,10,9,7,5,3,6,8,11

/* level 7 */
unpack		4,3,5,3
unpack		10,6,4,6
unpack		9,8,10,8
unpack		7,11,9,11

vpbroadcastq	(_ZETAS_PINV+128+16*\off)*2(%rdx),%ymm1
vpbroadcastq	(_ZETAS+128+16*\off)*2(%rdx),%ymm2
fqmulprecomp	1,2,3,x=7
vpbroadcastq	(_ZETAS_PINV+128+16*\off+4)*2(%rdx),%ymm1
vpbroadcastq	(_ZETAS+128+16*\off+4)*2(%rdx),%ymm2
fqmulprecomp	1,2,6,x=7
vpbroadcastq	(_ZETAS_PINV+128+16*\off+8)*2(%rdx),%ymm1
vpbroadcastq	(_ZETAS+128+16*\off+8)*2(%rdx),%ymm2
fqmulprecomp	1,2,8,x=7
vpbroadcastq	(_ZETAS_PINV+128+16*\off+12)*2(%rdx),%ymm1
vpbroadcastq	(_ZETAS+128+16*\off+12)*2(%rdx),%ymm2
fqmulprecomp	1,2,11,x=7

update		7,5,4,10,9,3,6,8,11

shuffle4	7,5,9,5
shuffle4	3,6,7,6
shuffle4	4,10,3,10
shuffle4	8,11,4,11

shuffle8	9,3,8,3
shuffle8	7,4,9,4
shuffle8	5,10,7,10
shuffle8	6,11,5,11

vmovdqa		%ymm8,(128*\off+  0)*2(%rdi)
vmovdqa		%ymm3,(128*\off+ 16)*2(%rdi)
vmovdqa		%ymm7,(128*\off+ 32)*2(%rdi)
vmovdqa		%ymm10,(128*\off+ 48)*2(%rdi)
vmovdqa		%ymm9,(128*\off+ 64)*2(%rdi)
vmovdqa		%ymm4,(128*\off+ 80)*2(%rdi)
vmovdqa		%ymm5,(128*\off+ 96)*2(%rdi)
vmovdqa		%ymm11,(128*\off+112)*2(%rdi)
.endm

.global cdecl(poly_ntt)
cdecl(poly_ntt):
vmovdqa		_16XP*2(%rdx),%ymm0

levels0t2	0
levels0t2	1
levels0t2	2
levels0t2	3
levels0t2	4
levels0t2	5
levels0t2	6
levels0t2	7

levels3t7	0
levels3t7	1
levels3t7	2
levels3t7	3
levels3t7	4
levels3t7	5
levels3t7	6
levels3t7	7

ret
