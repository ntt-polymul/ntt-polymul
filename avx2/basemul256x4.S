#include "consts1024.h"
.include "fq.inc"

.macro schoolbook off
vmovdqa		(64*\off+ 0)*2(%rsi),%ymm2
vmovdqa		(64*\off+16)*2(%rsi),%ymm3
vmovdqa		(64*\off+32)*2(%rsi),%ymm4
vmovdqa		(64*\off+48)*2(%rsi),%ymm5

vmovdqa		(_16XMONT_PINV)*2(%rcx),%ymm6
vmovdqa		(_16XMONT)*2(%rcx),%ymm7
fqmulprecomp	6,7,2,x=15
fqmulprecomp	6,7,3,x=15
fqmulprecomp	6,7,4,x=15
fqmulprecomp	6,7,5,x=15

vmovdqa		(64*\off+ 0)*2(%rdx),%ymm6
vpmullw		%ymm1,%ymm6,%ymm7
fqmulprecomp2	7,6,2,8,x=15
fqmulprecomp2	7,6,3,9,x=15
fqmulprecomp2	7,6,4,10,x=15
fqmulprecomp2	7,6,5,11,x=15

vmovdqa		(64*\off+16)*2(%rdx),%ymm6
vpmullw		%ymm1,%ymm6,%ymm7
fqmulprecomp2	7,6,2,14,x=15
vpaddw		%ymm14,%ymm9,%ymm9
fqmulprecomp2	7,6,3,14,x=15
vpaddw		%ymm14,%ymm10,%ymm10
fqmulprecomp2	7,6,4,14,x=15
vpaddw		%ymm14,%ymm11,%ymm11
fqmulprecomp2	7,6,5,12,x=15

vmovdqa		(64*\off+32)*2(%rdx),%ymm6
vpmullw		%ymm1,%ymm6,%ymm7
fqmulprecomp2	7,6,2,14,x=15
vpaddw		%ymm14,%ymm10,%ymm10
fqmulprecomp2	7,6,3,14,x=15
vpaddw		%ymm14,%ymm11,%ymm11
fqmulprecomp2	7,6,4,14,x=15
vpaddw		%ymm14,%ymm12,%ymm12
fqmulprecomp2	7,6,5,13,x=15

vmovdqa		(64*\off+48)*2(%rdx),%ymm6
vpmullw		%ymm1,%ymm6,%ymm7
fqmulprecomp2	7,6,2,14,x=15
vpaddw		%ymm14,%ymm11,%ymm11
fqmulprecomp2	7,6,3,14,x=15
vpaddw		%ymm14,%ymm12,%ymm12
fqmulprecomp2	7,6,4,14,x=15
vpaddw		%ymm14,%ymm13,%ymm13
fqmulprecomp2	7,6,5,14,x=15

vmovdqa		(_ZETAS_PINV+128+16*(\off/2))*2(%rcx),%ymm6
vmovdqa		(_ZETAS+128+16*(\off/2))*2(%rcx),%ymm7
fqmulprecomp	6,7,12,x=15
fqmulprecomp	6,7,13,x=15
fqmulprecomp	6,7,14,x=15
.if \off % 2 == 0
vpaddw		%ymm12,%ymm8,%ymm8
vpaddw		%ymm13,%ymm9,%ymm9
vpaddw		%ymm14,%ymm10,%ymm10
.else
vpsubw		%ymm12,%ymm8,%ymm8
vpsubw		%ymm13,%ymm9,%ymm9
vpsubw		%ymm14,%ymm10,%ymm10
.endif

vmovdqa		(_16XF_PINV)*2(%rcx),%ymm6
vmovdqa		(_16XF)*2(%rcx),%ymm7
fqmulprecomp	6,7,8,x=15
fqmulprecomp	6,7,9,x=15
fqmulprecomp	6,7,10,x=15
fqmulprecomp	6,7,11,x=15

vmovdqa		%ymm8,(64*\off+ 0)*2(%rdi)
vmovdqa		%ymm9,(64*\off+16)*2(%rdi)
vmovdqa		%ymm10,(64*\off+32)*2(%rdi)
vmovdqa		%ymm11,(64*\off+48)*2(%rdi)
.endm

.text
.global cdecl(poly_basemul_montgomery)
cdecl(poly_basemul_montgomery):
vmovdqa		(_16XP)*2(%rcx),%ymm0
vmovdqa		(_16XPINV)*2(%rcx),%ymm1

schoolbook	0
schoolbook	1
schoolbook	2
schoolbook	3
schoolbook	4
schoolbook	5
schoolbook	6
schoolbook	7
schoolbook	8
schoolbook	9
schoolbook	10
schoolbook	11
schoolbook	12
schoolbook	13
schoolbook	14
schoolbook	15

ret
