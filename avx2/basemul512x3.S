#include "consts1536.h"
.include "fq.inc"

.macro schoolbook off
vmovdqa		(48*\off+ 0)*2(%rsi),%ymm2		# a0
vmovdqa		(48*\off+16)*2(%rsi),%ymm3		# a1
vmovdqa		(48*\off+32)*2(%rsi),%ymm4		# a2

vpmullw		%ymm1,%ymm2,%ymm5			# a0.pinv
vpmullw		%ymm1,%ymm3,%ymm6			# a1.pinv
vpmullw		%ymm1,%ymm4,%ymm7			# a2.pinv

vmovdqa		(48*\off+ 0)*2(%rdx),%ymm8		# b0
vmovdqa		(48*\off+16)*2(%rdx),%ymm9		# b1
vmovdqa		(48*\off+32)*2(%rdx),%ymm10		# b2

fqmulprecomp2	5,2,8,11,x=15				# a0b0
fqmulprecomp2	5,2,9,12,x=15				# a0b1
fqmulprecomp2	5,2,10,13,x=15				# a0b2
fqmulprecomp2	6,3,8,14,x=15				# a1b0
fqmulprecomp2	6,3,9,2,x=15				# a1b1
fqmulprecomp2	6,3,10,5,x=15				# a1b2
fqmulprecomp2	7,4,8,3,x=15				# a2b0
fqmulprecomp2	7,4,9,6,x=15				# a2b1
fqmulprecomp2	7,4,10,8,x=15				# a2b2

vpaddw		%ymm14,%ymm12,%ymm12			# a0b1 + a1b0
vpaddw		%ymm3,%ymm13,%ymm13			# a0b2 + a2b0
vpaddw		%ymm2,%ymm13,%ymm13			# a0b2 + a2b0 + a1b1
vpaddw		%ymm6,%ymm5,%ymm5			# a1b2 + a2b1

.if (\off % 4) == 0
vpaddw		%ymm5,%ymm11,%ymm11			# c0
vpaddw		%ymm8,%ymm12,%ymm12			# c1
.elseif (\off % 4) == 1
vpsubw		%ymm5,%ymm11,%ymm11			# c0
vpsubw		%ymm8,%ymm12,%ymm12			# c1
.elseif (\off % 4) == 2
vpbroadcastd	(_ZETAS_PINV)*2(%rcx),%ymm9
vpbroadcastd	(_ZETAS)*2(%rcx),%ymm10
fqmulprecomp	9,10,5,x=15
fqmulprecomp	9,10,8,x=15
vpaddw		%ymm5,%ymm11,%ymm11			# c0
vpaddw		%ymm8,%ymm12,%ymm12			# c1
.elseif (\off % 4) == 3
vpbroadcastd	(_ZETAS_PINV)*2(%rcx),%ymm9
vpbroadcastd	(_ZETAS)*2(%rcx),%ymm10
fqmulprecomp	9,10,5,x=15
fqmulprecomp	9,10,8,x=15
vpsubw		%ymm5,%ymm11,%ymm11			# c0
vpsubw		%ymm8,%ymm12,%ymm12			# c1
.endif

vmovdqa		(_16XF_PINV)*2(%rcx),%ymm2
vmovdqa		(_16XF)*2(%rcx),%ymm3
fqmulprecomp	2,3,11,x=15
fqmulprecomp	2,3,12,x=15
fqmulprecomp	2,3,13,x=15

vmovdqa		%ymm11,(48*\off+ 0)*2(%rdi)
vmovdqa		%ymm12,(48*\off+16)*2(%rdi)
vmovdqa		%ymm13,(48*\off+32)*2(%rdi)
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
schoolbook	16
schoolbook	17
schoolbook	18
schoolbook	19
schoolbook	20
schoolbook	21
schoolbook	22
schoolbook	23
schoolbook	24
schoolbook	25
schoolbook	26
schoolbook	27
schoolbook	28
schoolbook	29
schoolbook	30
schoolbook	31

ret
