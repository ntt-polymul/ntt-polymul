	.thumb
	.syntax unified
.text

.global school_book_mul2_16_fast
.func school_book_mul2_16_fast, school_book_mul2_16_fast
.type school_book_mul2_16_fast, %function
.align 2
school_book_mul2_16_fast:
	stmdb	sp!, {lr}
	ldr.w r6, [r1, #0]
	ldr.w ip, [r1, #4]
	ldr.w r3, [r1, #8]
	ldr.w sl, [r1, #12]
	ldrh.w r9, [r2]
	ldrh.w fp, [r2, #2]
	ldr.w r7, [r0, #0]
	ldr.w r8, [r0, #4]
	ldr.w r4, [r0, #8]
	ldr.w lr, [r0, #12]
	smlabb r9, r7, r6, r9
	smladx fp, r7, r6, fp
	pkhbt r9, r9, fp, lsl #16
	ldrh.w fp, [r2, #6]
	ldrh.w r5, [r2, #4]
	str.w r9, [r2]
	smladx fp, r7, ip, fp
	smlabb r5, r7, ip, r5
	pkhbt r9, r8, r7
	smladx fp, r8, r6, fp
	smlad r5, r9, r6, r5
	pkhbt fp, r5, fp, lsl #16
	str.w fp, [r2, #4]
	ldrh.w r5, [r2, #8]
	ldrh.w fp, [r2, #10]
	smlabb r5, r3, r7, r5
	smladx fp, r3, r7, fp
	smlad r5, r9, ip, r5
	pkhbt r9, r4, r8
	smlad r5, r9, r6, r5
	smladx fp, r8, ip, fp
	smladx fp, r4, r6, fp
	pkhbt fp, r5, fp, lsl #16
	str.w fp, [r2, #8]
	ldrh.w r5, [r2, #12]
	ldrh.w fp, [r2, #14]
	smlabb r5, sl, r7, r5
	smladx fp, sl, r7, fp
	smlad r5, ip, r9, r5
	pkhbt r9, r8, r7
	smlad r5, r3, r9, r5
	pkhbt r9, lr, r4
	smlad r5, r6, r9, r5
	smladx fp, r3, r8, fp
	smladx fp, ip, r4, fp
	smladx fp, r6, lr, fp
	pkhbt fp, r5, fp, lsl #16
	str.w fp, [r2, #12]
	ldrh.w r5, [r2, #16]
	ldrh.w fp, [r2, #18]
	smlatt r5, r6, lr, r5
	smladx fp, sl, r8, fp
	smlad r5, ip, r9, r5
	pkhbt r9, r8, r7
	smlad r5, sl, r9, r5
	pkhbt r9, r4, r8
	smlad r5, r3, r9, r5
	smladx fp, r3, r4, fp
	smladx fp, ip, lr, fp
	pkhbt fp, r5, fp, lsl #16
	str.w fp, [r2, #16]
	ldrh.w r5, [r2, #20]
	ldrh.w fp, [r2, #22]
	smlatt r5, ip, lr, r5
	smladx fp, sl, r4, fp
	smlad r5, sl, r9, r5
	pkhbt r9, lr, r4
	smlad r5, r3, r9, r5
	smladx fp, r3, lr, fp
	pkhbt fp, r5, fp, lsl #16
	str.w fp, [r2, #20]
	ldrh.w r5, [r2, #24]
	ldrh.w fp, [r2, #26]
	smlad r5, sl, r9, r5
	smladx fp, sl, lr, fp
	smlatt r5, r3, lr, r5
	pkhbt fp, r5, fp, lsl #16
	ldr.w r9, [r2, #28]
	str.w fp, [r2, #24]
	smlatt fp, sl, lr, r9
	pkhbt fp, fp, r9
	str.w fp, [r2, #28]
	ldr.w r7, [r0, #16]
	ldr.w r8, [r0, #20]
	ldr.w r4, [r0, #24]
	ldr.w lr, [r0, #28]
	ldr.w r9, [r2, #16]
	ldr.w r5, [r2, #20]
	mov.w fp, r9, lsr #16
	smlabb r9, r7, r6, r9
	smladx fp, r7, r6, fp
	pkhbt r9, r9, fp, lsl #16
	str.w r9, [r2, #16]
	mov.w fp, r5, lsr #16
	smladx fp, r7, ip, fp
	smlabb r5, r7, ip, r5
	pkhbt r9, r8, r7
	smladx fp, r8, r6, fp
	smlad r5, r9, r6, r5
	pkhbt fp, r5, fp, lsl #16
	ldr.w r5, [r2, #24]
	str.w fp, [r2, #20]
	mov.w fp, r5, lsr #16
	smlabb r5, r3, r7, r5
	smladx fp, r3, r7, fp
	smlad r5, r9, ip, r5
	pkhbt r9, r4, r8
	smlad r5, r9, r6, r5
	smladx fp, r8, ip, fp
	smladx fp, r4, r6, fp
	pkhbt fp, r5, fp, lsl #16
	ldr.w r5, [r2, #28]
	str.w fp, [r2, #24]
	mov.w fp, r5, lsr #16
	smlabb r5, sl, r7, r5
	smladx fp, sl, r7, fp
	smlad r5, ip, r9, r5
	pkhbt r9, r8, r7
	smlad r5, r3, r9, r5
	pkhbt r9, lr, r4
	smlad r5, r6, r9, r5
	smladx fp, r3, r8, fp
	smladx fp, ip, r4, fp
	smladx fp, r6, lr, fp
	pkhbt fp, r5, fp, lsl #16
	str.w fp, [r2, #28]
	ldrh.w r5, [r2, #32]
	ldrh.w fp, [r2, #34]
	smlatt r5, r6, lr, r5
	smladx fp, sl, r8, fp
	smlad r5, ip, r9, r5
	pkhbt r9, r8, r7
	smlad r5, sl, r9, r5
	pkhbt r9, r4, r8
	smlad r5, r3, r9, r5
	smladx fp, r3, r4, fp
	smladx fp, ip, lr, fp
	pkhbt fp, r5, fp, lsl #16
	str.w fp, [r2, #32]
	ldrh.w r5, [r2, #36]
	ldrh.w fp, [r2, #38]
	smlatt r5, ip, lr, r5
	smladx fp, sl, r4, fp
	smlad r5, sl, r9, r5
	pkhbt r9, lr, r4
	smlad r5, r3, r9, r5
	smladx fp, r3, lr, fp
	pkhbt fp, r5, fp, lsl #16
	str.w fp, [r2, #36]
	ldrh.w r5, [r2, #40]
	ldrh.w fp, [r2, #42]
	smlad r5, sl, r9, r5
	smladx fp, sl, lr, fp
	smlatt r5, r3, lr, r5
	pkhbt fp, r5, fp, lsl #16
	ldr.w r9, [r2, #44]
	str.w fp, [r2, #40]
	smlatt fp, sl, lr, r9
	pkhbt fp, fp, r9
	str.w fp, [r2, #44]
	ldr.w r6, [r1, #16]
	ldr.w ip, [r1, #20]
	ldr.w r3, [r1, #24]
	ldr.w sl, [r1, #28]
	ldr.w r9, [r2, #32]
	ldr.w r5, [r2, #36]
	mov.w fp, r9, lsr #16
	smlabb r9, r7, r6, r9
	smladx fp, r7, r6, fp
	pkhbt r9, r9, fp, lsl #16
	str.w r9, [r2, #32]
	mov.w fp, r5, lsr #16
	smladx fp, r7, ip, fp
	smlabb r5, r7, ip, r5
	pkhbt r9, r8, r7
	smladx fp, r8, r6, fp
	smlad r5, r9, r6, r5
	pkhbt fp, r5, fp, lsl #16
	ldr.w r5, [r2, #40]
	str.w fp, [r2, #36]
	mov.w fp, r5, lsr #16
	smlabb r5, r3, r7, r5
	smladx fp, r3, r7, fp
	smlad r5, r9, ip, r5
	pkhbt r9, r4, r8
	smlad r5, r9, r6, r5
	smladx fp, r8, ip, fp
	smladx fp, r4, r6, fp
	pkhbt fp, r5, fp, lsl #16
	ldr.w r5, [r2, #44]
	str.w fp, [r2, #40]
	mov.w fp, r5, lsr #16
	smlabb r5, sl, r7, r5
	smladx fp, sl, r7, fp
	smlad r5, ip, r9, r5
	pkhbt r9, r8, r7
	smlad r5, r3, r9, r5
	pkhbt r9, lr, r4
	smlad r5, r6, r9, r5
	smladx fp, r3, r8, fp
	smladx fp, ip, r4, fp
	smladx fp, r6, lr, fp
	pkhbt fp, r5, fp, lsl #16
	str.w fp, [r2, #44]
	ldrh.w r5, [r2, #48]
	ldrh.w fp, [r2, #50]
	smlatt r5, r6, lr, r5
	smladx fp, sl, r8, fp
	smlad r5, ip, r9, r5
	pkhbt r9, r8, r7
	smlad r5, sl, r9, r5
	pkhbt r9, r4, r8
	smlad r5, r3, r9, r5
	smladx fp, r3, r4, fp
	smladx fp, ip, lr, fp
	pkhbt fp, r5, fp, lsl #16
	str.w fp, [r2, #48]
	ldrh.w r5, [r2, #52]
	ldrh.w fp, [r2, #54]
	smlatt r5, ip, lr, r5
	smladx fp, sl, r4, fp
	smlad r5, sl, r9, r5
	pkhbt r9, lr, r4
	smlad r5, r3, r9, r5
	smladx fp, r3, lr, fp
	pkhbt fp, r5, fp, lsl #16
	str.w fp, [r2, #52]
	ldrh.w r5, [r2, #56]
	ldrh.w fp, [r2, #58]
	smlad r5, sl, r9, r5
	smladx fp, sl, lr, fp
	smlatt r5, r3, lr, r5
	ldrh.w r9, [r2, #60]
	pkhbt fp, r5, fp, lsl #16
	str.w fp, [r2, #56]
	smlatt fp, sl, lr, r9
	strh.w fp, [r2, #60]
	ldr.w r7, [r0, #0]
	ldr.w r8, [r0, #4]
	ldr.w r4, [r0, #8]
	ldr.w lr, [r0, #12]
	ldr.w r9, [r2, #16]
	ldr.w r5, [r2, #20]
	mov.w fp, r9, lsr #16
	smlabb r9, r7, r6, r9
	smladx fp, r7, r6, fp
	pkhbt r9, r9, fp, lsl #16
	str.w r9, [r2, #16]
	mov.w fp, r5, lsr #16
	smladx fp, r7, ip, fp
	smlabb r5, r7, ip, r5
	pkhbt r9, r8, r7
	smladx fp, r8, r6, fp
	smlad r5, r9, r6, r5
	pkhbt fp, r5, fp, lsl #16
	ldr.w r5, [r2, #24]
	str.w fp, [r2, #20]
	mov.w fp, r5, lsr #16
	smlabb r5, r3, r7, r5
	smladx fp, r3, r7, fp
	smlad r5, r9, ip, r5
	pkhbt r9, r4, r8
	smlad r5, r9, r6, r5
	smladx fp, r8, ip, fp
	smladx fp, r4, r6, fp
	pkhbt fp, r5, fp, lsl #16
	ldr.w r5, [r2, #28]
	str.w fp, [r2, #24]
	mov.w fp, r5, lsr #16
	smlabb r5, sl, r7, r5
	smladx fp, sl, r7, fp
	smlad r5, ip, r9, r5
	pkhbt r9, r8, r7
	smlad r5, r3, r9, r5
	pkhbt r9, lr, r4
	smlad r5, r6, r9, r5
	smladx fp, r3, r8, fp
	smladx fp, ip, r4, fp
	smladx fp, r6, lr, fp
	pkhbt fp, r5, fp, lsl #16
	ldr.w r5, [r2, #32]
	str.w fp, [r2, #28]
	mov.w fp, r5, lsr #16
	smlatt r5, r6, lr, r5
	smladx fp, sl, r8, fp
	smlad r5, ip, r9, r5
	pkhbt r9, r8, r7
	smlad r5, sl, r9, r5
	pkhbt r9, r4, r8
	smlad r5, r3, r9, r5
	smladx fp, r3, r4, fp
	smladx fp, ip, lr, fp
	pkhbt fp, r5, fp, lsl #16
	ldr.w r5, [r2, #36]
	str.w fp, [r2, #32]
	mov.w fp, r5, lsr #16
	smlatt r5, ip, lr, r5
	smladx fp, sl, r4, fp
	smlad r5, sl, r9, r5
	pkhbt r9, lr, r4
	smlad r5, r3, r9, r5
	smladx fp, r3, lr, fp
	pkhbt fp, r5, fp, lsl #16
	ldr.w r5, [r2, #40]
	str.w fp, [r2, #36]
	mov.w fp, r5, lsr #16
	smlad r5, sl, r9, r5
	smladx fp, sl, lr, fp
	smlatt r5, r3, lr, r5
	pkhbt fp, r5, fp, lsl #16
	ldrh.w r5, [r2, #44]
	str.w fp, [r2, #40]
	smlatt fp, sl, lr, r5
	strh.w fp, [r2, #44]
	ldmia.w	sp!, {pc}
.endfunc

.global eval_kara
.func eval_kara, eval_kara
.type eval_kara, %function
.align 2
eval_kara:
	stmdb sp!, {r2, r3, lr}
	sub.w sp, #164
//iter 1/4
	ldr r4, [r0]
	ldr r5, [r0, #32]
	ldr r6, [r0, #64]
	ldr r7, [r0, #96]
	ldr r8, [r0, #4]
	ldr r9, [r0, #36]
	ldr r10, [r0, #68]
	ldr r11, [r0, #100]
	uadd16 r2, r4, r5
	uadd16 r3, r8, r9
	uadd16 r12, r6, r7
	uadd16 r14, r10, r11
	uadd16 r4, r4, r6
	uadd16 r8, r8, r10
	uadd16 r5, r5, r7
	uadd16 r9, r9, r11
	uadd16 r6, r2, r12
	uadd16 r7, r3, r14
	str r2, [sp]
	str r3, [sp, #4]
	str r12, [sp, #32]
	str r14, [sp, #36]
	str r4, [sp, #64]
	str r5, [sp, #96]
	str r8, [sp, #68]
	str r9, [sp, #100]
	str r6, [sp, #128]
	str r7, [sp, #132]
//iter 2/4
	ldr r4, [r0, #8]
	ldr r5, [r0, #40]
	ldr r6, [r0, #72]
	ldr r7, [r0, #104]
	ldr r8, [r0, #12]
	ldr r9, [r0, #44]
	ldr r10, [r0, #76]
	ldr r11, [r0, #108]
	uadd16 r2, r4, r5
	uadd16 r3, r8, r9
	uadd16 r12, r6, r7
	uadd16 r14, r10, r11
	uadd16 r4, r4, r6
	uadd16 r8, r8, r10
	uadd16 r5, r5, r7
	uadd16 r9, r9, r11
	uadd16 r6, r2, r12
	uadd16 r7, r3, r14
	str r2, [sp, #8]
	str r3, [sp, #12]
	str r12, [sp, #40]
	str r14, [sp, #44]
	str r4, [sp, #72]
	str r5, [sp, #104]
	str r8, [sp, #76]
	str r9, [sp, #108]
	str r6, [sp, #136]
	str r7, [sp, #140]
//iter 3/4
	ldr r4, [r0, #16]
	ldr r5, [r0, #48]
	ldr r6, [r0, #80]
	ldr r7, [r0, #112]
	ldr r8, [r0, #20]
	ldr r9, [r0, #52]
	ldr r10, [r0, #84]
	ldr r11, [r0, #116]
	uadd16 r2, r4, r5
	uadd16 r3, r8, r9
	uadd16 r12, r6, r7
	uadd16 r14, r10, r11
	uadd16 r4, r4, r6
	uadd16 r8, r8, r10
	uadd16 r5, r5, r7
	uadd16 r9, r9, r11
	uadd16 r6, r2, r12
	uadd16 r7, r3, r14
	str r2, [sp, #16]
	str r3, [sp, #20]
	str r12, [sp, #48]
	str r14, [sp, #52]
	str r4, [sp, #80]
	str r5, [sp, #112]
	str r8, [sp, #84]
	str r9, [sp, #116]
	str r6, [sp, #144]
	str r7, [sp, #148]
//iter 4/4
	ldr r4, [r0, #24]
	ldr r5, [r0, #56]
	ldr r6, [r0, #88]
	ldr r7, [r0, #120]
	ldr r8, [r0, #28]
	ldr r9, [r0, #60]
	ldr r10, [r0, #92]
	ldr r11, [r0, #124]
	uadd16 r2, r4, r5
	uadd16 r3, r8, r9
	uadd16 r12, r6, r7
	uadd16 r14, r10, r11
	uadd16 r4, r4, r6
	uadd16 r8, r8, r10
	uadd16 r5, r5, r7
	uadd16 r9, r9, r11
	uadd16 r6, r2, r12
	uadd16 r7, r3, r14
	str r2, [sp, #24]
	str r3, [sp, #28]
	str r12, [sp, #56]
	str r14, [sp, #60]
	str r4, [sp, #88]
	str r5, [sp, #120]
	ldr.w r2, [sp, #164]
	str r8, [sp, #92]
	str r9, [sp, #124]
	add.w r2, #186
	add.w r1, #160
	str r6, [sp, #152]
	str r7, [sp, #156]
//end iter
// r0 = a   r1 = bw_ar   r2 = res   sp=aw
	bl school_book_mul2_16_fast
	add r2, #62
	add r1, #32
	add r0, #32
	bl school_book_mul2_16_fast
	add r2, #62
	add r1, #32
	add r0, #32
	bl school_book_mul2_16_fast
	add r2, #62
	add r1, #32
	add r0, #32
	bl school_book_mul2_16_fast
	sub r2, #372
	sub r1, #256
	mov r0, sp
	bl school_book_mul2_16_fast
	add r2, #62
	add r1, #32
	add r0, #32
	bl school_book_mul2_16_fast
	add r2, #372
	add r1, #32
	add r0, #32
	bl school_book_mul2_16_fast
	add r2, #62
	add r1, #32
	add r0, #32
	bl school_book_mul2_16_fast
	sub r2, #372
	add r1, #32
	add r0, #32
	bl school_book_mul2_16_fast
	add.w sp, #164
	ldmia.w sp!, {r2, r3, pc}
.endfunc

.global TC_evaluation_unrolled_kara
.func TC_evaluation_unrolled_kara, TC_evaluation_unrolled_kara
.type TC_evaluation_unrolled_kara, %function
.align 2
TC_evaluation_unrolled_kara:
	stmdb sp!, {r4, r5, r6, r7, r8, r9, r10, r11, r12, lr}
	sub.w sp, #644
//sp = aw   r0 = a  r1 = bw_wr   r2 = w_ar
//loop1
	add 	r12, r0, #128
TC_evaluation_unrolled_kara_loop1:
	ldr 	r6, [r0, #384]
	ldr 	r4, [r0, #128]
	ldr 	r5, [r0, #256]
	ldr 	r3, [r0], #4
	cmp 	r0, r12
	uadd16 	r8, r4, r6
	uadd16 	r7, r3, r5
	uadd16 	r9, r7, r8
	usub16 	r10, r7, r8
	str 	r9, [sp, #128]
	str 	r10, [sp, #256]
	uadd16 	r7, r3, r3
	uadd16	r7, r7, r7
	uadd16	r7, r7, r5
	uadd16 	r7, r7, r7
	uadd16 	r8, r4, r4
	uadd16 	r8, r8, r8
	uadd16 	r8, r8, r6
	uadd16 	r9, r7, r8
	usub16 	r10, r7, r8
	str 	r9, [sp, #384]
	str 	r10, [sp, #512]
	uadd16 	r7, r6, r6
	uadd16 	r7, r7, r5
	uadd16 	r7, r7, r7
	uadd16 	r7, r7, r4
	uadd16 	r7, r7, r7
	uadd16 	r7, r7 ,r3
	str 	r7, [sp], #4
	bne 	TC_evaluation_unrolled_kara_loop1
//sp = aw[64]   r0 = a[64]  r1 = bw_wr   r2 = w_ar
	add r0, #256
	sub sp, #128
	sub r3, r0, #384
//sp = aw   r0 = a[192]  r1 = bw_wr   r2 = w_ar
//bw_wr -> +288(16*9*2) -> 288-128=+160  w_ar -> +558(31*9*2)
	bl eval_kara
//sp = aw   r0=corrupted   r1=passed bw_wr + 128   r2=passed w_ar   r3=a
	add r2, #558
	add r1, #160
	mov r0, sp
	bl eval_kara
	add r2, #558
	add r1, #160
	add r0, sp, #128
	bl eval_kara
	add r2, #558
	add r1, #160
	add r0, sp, #256
	bl eval_kara
	add r2, #558
	add r1, #160
	add r0, sp, #384
	bl eval_kara
	add r2, #558
	add r1, #160
	add r0, sp, #512
	bl eval_kara
	add r2, #558
	add r1, #160
	mov r0, r3
	bl eval_kara
	add.w sp, #644
	ldmia.w sp!, {r4, r5, r6, r7, r8, r9, r10, r11, r12, pc}
.endfunc

.global evaluation_single_kara
.func evaluation_single_kara, evaluation_single_kara
.type evaluation_single_kara, %function
.align 2
evaluation_single_kara:
	stmdb sp!, {r4, r5, r6, r7, r8, r9, r10, r11, r12, lr}
	sub.w sp, #644
//r0 = b  r1 = bw_ar   sp = bw
//loop1
	add 	r12, r0, #128
evaluation_single_kara_asm_loop1:
	ldr 	r6, [r0, #384]
	ldr 	r4, [r0, #128]
	ldr 	r5, [r0, #256]
	ldr 	r3, [r0], #4
	cmp 	r0, r12
	uadd16 	r8, r4, r6
	uadd16 	r7, r3, r5
	uadd16 	r9, r7, r8
	usub16 	r10, r7, r8
	str 	r9, [sp, #128]
	str 	r10, [sp, #256]
	uadd16 	r7, r3, r3
	uadd16	r7, r7, r7
	uadd16	r7, r7, r5
	uadd16 	r7, r7, r7
	uadd16 	r8, r4, r4
	uadd16 	r8, r8, r8
	uadd16 	r8, r8, r6
	uadd16 	r9, r7, r8
	usub16 	r10, r7, r8
	str 	r9, [sp, #384]
	str 	r10, [sp, #512]
	uadd16 	r7, r6, r6
	uadd16 	r7, r7, r5
	uadd16 	r7, r7, r7
	uadd16 	r7, r7, r4
	uadd16 	r7, r7, r7
	uadd16 	r7, r7 ,r3
	str 	r7, [sp], #4
	bne 	evaluation_single_kara_asm_loop1
//r0 = b[64]  r1 = bw_ar   r2 = bw[64]
	add r0, #256
	sub sp, #128
	bl kara_eval_precomp
	add r1, #1728
	sub r0, #384
	bl kara_eval_precomp
	sub r1, #1440
	mov r0, sp
	bl kara_eval_precomp
	add r1, #288
	add r0, #128
	bl kara_eval_precomp
	add r1, #288
	add r0, #128
	bl kara_eval_precomp
	add r1, #288
	add r0, #128
	bl kara_eval_precomp
	add r1, #288
	add r0, #128
	bl kara_eval_precomp
// kara_eval_precomp(&b[3*N_SB],  bw_ar[0]);
// kara_eval_precomp(&b[0],       bw_ar[6]);
// kara_eval_precomp(bw,          bw_ar[1]);
// kara_eval_precomp(&bw[N_SB],   bw_ar[2]);
// kara_eval_precomp(&bw[2*N_SB], bw_ar[3]);
// kara_eval_precomp(&bw[3*N_SB], bw_ar[4]);
// kara_eval_precomp(&bw[4*N_SB], bw_ar[5]);
	add.w sp, #644
	ldmia.w sp!, {r4, r5, r6, r7, r8, r9, r10, r11, r12, pc}
.endfunc

.global kara_eval_precomp
.func kara_eval_precomp, kara_eval_precomp
.type kara_eval_precomp, %function
.align 2
kara_eval_precomp:
	stmdb sp!, {lr}
//iter 1/4
	ldr r4, [r0]
	ldr r5, [r0, #32]
	ldr r6, [r0, #64]
	ldr r7, [r0, #96]
	ldr r8, [r0, #4]
	ldr r9, [r0, #36]
	ldr r10, [r0, #68]
	ldr r11, [r0, #100]
	uadd16 r2, r4, r5
	uadd16 r3, r8, r9
	uadd16 r12, r6, r7
	uadd16 r14, r10, r11
	str r4, [r1, #160]
	str r5, [r1, #192]
	str r6, [r1, #224]
	str r7, [r1, #256]
	str r8, [r1, #164]
	str r9, [r1, #196]
	str r10, [r1, #228]
	str r11, [r1, #260]
	uadd16 r4, r4, r6
	uadd16 r8, r8, r10
	uadd16 r5, r5, r7
	uadd16 r9, r9, r11
	uadd16 r6, r2, r12
	uadd16 r7, r3, r14
	str r2, [r1]
	str r3, [r1, #4]
	str r12, [r1, #32]
	str r14, [r1, #36]
	str r4, [r1, #64]
	str r5, [r1, #96]
	str r8, [r1, #68]
	str r9, [r1, #100]
	str r6, [r1, #128]
	str r7, [r1, #132]
//iter 2/4
	ldr r4, [r0, #8]
	ldr r5, [r0, #40]
	ldr r6, [r0, #72]
	ldr r7, [r0, #104]
	ldr r8, [r0, #12]
	ldr r9, [r0, #44]
	ldr r10, [r0, #76]
	ldr r11, [r0, #108]
	uadd16 r2, r4, r5
	uadd16 r3, r8, r9
	uadd16 r12, r6, r7
	uadd16 r14, r10, r11
	str r4, [r1, #168]
	str r5, [r1, #200]
	str r6, [r1, #232]
	str r7, [r1, #264]
	str r8, [r1, #172]
	str r9, [r1, #204]
	str r10, [r1, #236]
	str r11, [r1, #268]
	uadd16 r4, r4, r6
	uadd16 r8, r8, r10
	uadd16 r5, r5, r7
	uadd16 r9, r9, r11
	uadd16 r6, r2, r12
	uadd16 r7, r3, r14
	str r2, [r1, #8]
	str r3, [r1, #12]
	str r12, [r1, #40]
	str r14, [r1, #44]
	str r4, [r1, #72]
	str r5, [r1, #104]
	str r8, [r1, #76]
	str r9, [r1, #108]
	str r6, [r1, #136]
	str r7, [r1, #140]
//iter 3/4
	ldr r4, [r0, #16]
	ldr r5, [r0, #48]
	ldr r6, [r0, #80]
	ldr r7, [r0, #112]
	ldr r8, [r0, #20]
	ldr r9, [r0, #52]
	ldr r10, [r0, #84]
	ldr r11, [r0, #116]
	uadd16 r2, r4, r5
	uadd16 r3, r8, r9
	uadd16 r12, r6, r7
	uadd16 r14, r10, r11
	str r4, [r1, #176]
	str r5, [r1, #208]
	str r6, [r1, #240]
	str r7, [r1, #272]
	str r8, [r1, #180]
	str r9, [r1, #212]
	str r10, [r1, #244]
	str r11, [r1, #276]
	uadd16 r4, r4, r6
	uadd16 r8, r8, r10
	uadd16 r5, r5, r7
	uadd16 r9, r9, r11
	uadd16 r6, r2, r12
	uadd16 r7, r3, r14
	str r2, [r1, #16]
	str r3, [r1, #20]
	str r12, [r1, #48]
	str r14, [r1, #52]
	str r4, [r1, #80]
	str r5, [r1, #112]
	str r8, [r1, #84]
	str r9, [r1, #116]
	str r6, [r1, #144]
	str r7, [r1, #148]
//iter 4/4
	ldr r4, [r0, #24]
	ldr r5, [r0, #56]
	ldr r6, [r0, #88]
	ldr r7, [r0, #120]
	ldr r8, [r0, #28]
	ldr r9, [r0, #60]
	ldr r10, [r0, #92]
	ldr r11, [r0, #124]
	uadd16 r2, r4, r5
	uadd16 r3, r8, r9
	uadd16 r12, r6, r7
	uadd16 r14, r10, r11
	str r4, [r1, #184]
	str r5, [r1, #216]
	str r6, [r1, #248]
	str r7, [r1, #280]
	str r8, [r1, #188]
	str r9, [r1, #220]
	str r10, [r1, #252]
	str r11, [r1, #284]
	uadd16 r4, r4, r6
	uadd16 r8, r8, r10
	uadd16 r5, r5, r7
	uadd16 r9, r9, r11
	uadd16 r6, r2, r12
	uadd16 r7, r3, r14
	str r2, [r1, #24]
	str r3, [r1, #28]
	str r12, [r1, #56]
	str r14, [r1, #60]
	str r4, [r1, #88]
	str r5, [r1, #120]
	str r8, [r1, #92]
	str r9, [r1, #124]
	str r6, [r1, #152]
	str r7, [r1, #156]
//end iter
	ldmia.w sp!, {pc}
.endfunc

.global TC_interpol1_kara
.func TC_interpol1_kara, TC_interpol1_kara
.type TC_interpol1_kara, %function
.align 2
TC_interpol1_kara:
	stmdb sp!, {r4, r5, r6, r7, r8, r9, r10, r11, r12, lr}
	subw	sp, sp, #1796
	eor r2, r2, r2
	strh.w	r2, [sp, #190]
	strh.w	r2, [sp, #126]
	strh.w	r2, [sp, #62]
	strh.w	r2, [sp, #446]
	strh.w	r2, [sp, #382]
	strh.w	r2, [sp, #318]
	strh.w	r2, [sp, #702]
	strh.w	r2, [sp, #638]
	strh.w	r2, [sp, #574]
	strh.w	r2, [sp, #958]
	strh.w	r2, [sp, #894]
	strh.w	r2, [sp, #830]
	strh.w	r2, [sp, #1214]
	strh.w	r2, [sp, #1150]
	strh.w	r2, [sp, #1086]
	strh.w	r2, [sp, #1470]
	strh.w	r2, [sp, #1406]
	strh.w	r2, [sp, #1342]
	strh.w	r2, [sp, #1726]
	strh.w	r2, [sp, #1662]
	strh.w	r2, [sp, #1598]
	strh.w	r2, [sp, #1790]
	strh.w	r2, [sp, #1534]
	strh.w	r2, [sp, #1278]
	strh.w	r2, [sp, #1022]
	strh.w	r2, [sp, #766]
	strh.w	r2, [sp, #510]
	strh.w	r2, [sp, #254]

	mov r4, r0
	mov r5, r1

	mov r1, sp
	bl interpol_kara
	add	r1, sp, #256
	addw	r0, r4, #558
	bl interpol_kara
	add	r1, sp, #512
	addw	r0, r4, #1116
	bl interpol_kara
	add	r1, sp, #768
	addw	r0, r4, #1674
	bl interpol_kara
	add.w	r1, sp, #1024
	addw	r0, r4, #2232
	bl interpol_kara
	add.w	r1, sp, #1280
	addw	r0, r4, #2790
	bl interpol_kara
	add.w	r1, sp, #1536
	addw	r0, r4, #3348
	bl interpol_kara

	mov r1, r5
	mov r0, sp

	eor r14, r14, r14
TC_interpol1_kara_asm_loop:
//r0=w_m r1=result r14=i
	ldr 	r8, [r0, #1536]
	ldr 	r7, [r0, #1280]
	ldr 	r6, [r0, #1024]
	ldr 	r5, [r0, #768]
	ldr 	r4, [r0, #512]
	ldr 	r3, [r0, #256]
	ldr 	r2, [r0], #4
	uadd16 r3, r3, r6
	usub16 r7, r7, r6
	usub16 r5, r5, r4
	lsr r5, r5, #1
	and r5, r5, #0x7FFF7FFF
	usub16 r6, r6, r2
	lsl r9, r8, #6
	and r9, r9, #0xFFC0FFC0
	usub16 r6, r6, r9
	lsl r6, r6, #1
	and r6, r6, #0xFFFEFFFE
	uadd16 r6, r6, r7
	uadd16 r4, r4, r5
	lsl r9, r4, #6
	and r9, r9, #0xFFC0FFC0
	usub16 r3, r3, r9
	usub16 r3, r3, r4
	usub16 r4, r4, r8
	usub16 r4, r4, r2
	mov r9, #45
	smulbb r10, r4, r9
	smultb r11, r4, r9
	pkhbt r9, r10, r11, lsl #16
	uadd16 r3, r3, r9
	lsl r9, r4, #3
	and r9, r9, #0xFFF8FFF8
	usub16 r6, r6, r9
	movw r9, #43691
	smulbb r10, r6, r9
	smultb r11, r6, r9
	pkhbt r9, r10, r11, lsl #16
	lsr r6, r9, #3
	and r6, r6, #0x1FFF1FFF
	uadd16 r7, r7, r3
	lsl r9, r5, #4
	and r9, r9, #0xFFF0FFF0
	uadd16 r3, r3, r9
	movw r9, #36409
	smulbb r10, r3, r9
	smultb r11, r3, r9
	pkhbt r9, r10, r11, lsl #16
	lsr r3, r9, #1
	and r3, r3, #0x7FFF7FFF
	eor r9, r9, r9
	usub16 r9, r9, r5
	usub16 r5, r9, r3
	mov r9, #30
	smulbb r10, r3, r9
	smultb r11, r3, r9
	pkhbt r9, r10, r11, lsl #16
	usub16 r7, r9, r7
	movw r9, #61167
	smulbb r10, r7, r9
	smultb r11, r7, r9
	pkhbt r9, r10, r11, lsl #16
	lsr r7, r9, #2
	and r7, r7, #0x3FFF3FFF
	usub16 r4, r4, r6
	usub16 r3, r3, r7
//stores
	ldr r9, [r1, #256]
	ldr r10, [r1, #128]
	ldr r11, [r1]
	uadd16 r9, r9, r6
	usub16 r9, r9, r2
	and r9, r9, #0x1FFF1FFF
	uadd16 r10, r10, r7
	usub16 r10, r10, r3
	and r10, r10, #0x1FFF1FFF
	uadd16 r11, r11, r8
	usub16 r11, r11, r4
	and r11, r11, #0x1FFF1FFF
	str r9, [r1, #256]
	str r10, [r1, #128]
	str r11, [r1], #4
	ldr r10, [r1, #-132]
	ldr r11, [r1, #380]
	usub16 r10, r10, r5
	uadd16 r11, r11, r5
	and r10, r10, #0x1FFF1FFF
	and r11, r11, #0x1FFF1FFF
	cmp r14, #64
	bge greater_equal_64
	str r11, [r1, #380]
	b cond_done
greater_equal_64:
	str r10, [r1, #-132]
cond_done:
	add r14, #2
	cmp r14, #128
	bne TC_interpol1_kara_asm_loop

	addw sp, sp, #1796

	ldmia.w sp!, {r4, r5, r6, r7, r8, r9, r10, r11, r12, pc}
.endfunc

.global interpol_kara
.func interpol_kara, interpol_kara
.type interpol_kara, %function
.align 2
interpol_kara:
	stmdb sp!, {r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, lr}
	sub sp, #132
//// OJO USAR R12 SOLO TEMPORALMENTE
	eor r12, r12, r12
	strh.w	r12, [sp, #62]
	mov r2, sp
////(res_ptr, result, mid_result)
interpol_kara_loop2:
	add r12, r12, #2
	ldr r3, [r0]
	ldr r4, [r0, #62]
	ldr r5, [r0, #124]
	ldr r6, [r0, #186]
	ldr r7, [r0, #248]
	ldr r8, [r0, #310]
	ldr r9, [r0, #372]
	ldr r10, [r0, #434]
	ldr r11, [r0, #496]
	cmp r12, #32
	usub16 r3, r3, r6
	usub16 r3, r3, r7
	usub16 r4, r4, r8
	usub16 r4, r4, r9
	usub16 r5, r5, r10
	usub16 r5, r5, r11
	beq exit_loop2
	str r5, [r0, #124]
	str r4, [r0, #62]
	str r3, [r0], #4
	str r9, [r1, #192]
	str r8, [r1, #128]
	str r7, [r1, #64]
	str r6, [r1], #4
	str r11, [r2, #64]
	str r10, [r2], #4
	b interpol_kara_loop2
exit_loop2:
	strh r5, [r0, #124]
	strh r4, [r0, #62]
	strh r3, [r0], #-60
	strh r9, [r1, #192]
	strh r8, [r1, #128]
	strh r7, [r1, #64]
	strh r6, [r1], #-28
	strh r11, [r2, #64]
	strh r10, [r2], #-28////(res_ptr, result, mid_result)
////(res_ptr, result[16], mid_result[16])
//// OJO USAR R12 SOLO TEMPORALMENTE
	eor r12, r12, r12
interpol_kara_loop3:
	add r12, r12, #2
	ldr r5, [r0, #124]
	ldr r4, [r0, #62]
	ldr r3, [r0], #4
	ldr r6, [r1]
	ldr r7, [r1, #128]
	ldr r8, [r2]
	cmp r12, #32
	uadd16 r6, r6, r3
	uadd16 r7, r7, r4
	uadd16 r8, r8, r5
	beq exit_loop3
	str r7, [r1, #128]
	str r6, [r1], #4
	str r8, [r2], #4
	b interpol_kara_loop3
exit_loop3:
	strh r7, [r1, #128]
	strh r6, [r1], #-92
	strh r8, [r2], #-92
////(res_ptr, result, mid_result)
//// OJO USAR R12 SOLO TEMPORALMENTE
	eor r12, r12, r12
interpol_kara_sub_loop:
	add r12, r12, #2
	ldr r8, [r1, #132]
	ldr r7, [r1, #4]
	ldr r6, [r2, #4]
	ldr r5, [r1, #128]
	ldr r4, [r1], #8
	ldr r3, [r2]
	cmp r12, #32
	usub16 r6, r6, r7
	usub16 r6, r6, r8
	usub16 r3, r3, r4
	usub16 r3, r3, r5
	beq exit_sub_loop
	str r6, [r2, #4]
	str r3, [r2], #8
	b interpol_kara_sub_loop
exit_sub_loop:
	strh r6, [r2, #4]
	str r3, [r2], #-120
////(res_ptr, result, mid_result[32])
	sub r1, r1, #64
//// OJO USAR R12 SOLO TEMPORALMENTE
	eor r12, r12, r12
interpol_kara_add_loop:
// 	result[i+32]=(result[i+32]+mid_result[i]);
	add r12, r12, #2
	ldr r3, [r1]
	ldr r4, [r1, #4]
	ldr r7, [r2], #8
	ldr r8, [r2, #-4]
	cmp r12, #32
	uadd16 r3, r3, r7
	uadd16 r4, r4, r8
	beq exit_add_loop
	str r4, [r1, #4]
	str r3, [r1], #8
	b interpol_kara_add_loop
exit_add_loop:
	strh r4, [r1, #4]
	str r3, [r1], #8
	add sp, #132
	ldmia.w sp!, {r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, pc}
.endfunc

//***************************
//* FUNCTIONS FOR DEBUGGING *
//***************************

.global TC_interpol1_kara_asm_round
.func TC_interpol1_kara_asm_round, TC_interpol1_kara_asm_round
.type TC_interpol1_kara_asm_round, %function
.align 2
TC_interpol1_kara_asm_round:
	stmdb sp!, {r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, lr}

	mov r14, r2

//r0=w_m r1=result r14=i
	ldr 	r8, [r0, #1536]
	ldr 	r7, [r0, #1280]
	ldr 	r6, [r0, #1024]
	ldr 	r5, [r0, #768]
	ldr 	r4, [r0, #512]
	ldr 	r3, [r0, #256]
	ldr 	r2, [r0], #2

	uadd16 r3, r3, r6
	usub16 r7, r7, r6
	usub16 r5, r5, r4
	lsr r5, r5, #1
	and r5, r5, #0x7FFF7FFF
	usub16 r6, r6, r2
	lsl r9, r8, #6
	and r9, r9, #0xFFC0FFC0
	usub16 r6, r6, r9
	lsl r6, r6, #1
	and r6, r6, #0xFFFEFFFE
	uadd16 r6, r6, r7
	uadd16 r4, r4, r5
	lsl r9, r4, #6
	and r9, r9, #0xFFC0FFC0
	usub16 r3, r3, r9
	usub16 r3, r3, r4
	usub16 r4, r4, r8
	usub16 r4, r4, r2
	mov r9, #45
	smulbb r10, r4, r9
	smultb r11, r4, r9
	pkhbt r9, r10, r11, lsl #16
	uadd16 r3, r3, r9
	lsl r9, r4, #3
	and r9, r9, #0xFFF8FFF8
	usub16 r6, r6, r9
	movw r9, #43691
	smulbb r10, r6, r9
	smultb r11, r6, r9
	pkhbt r9, r10, r11, lsl #16
	lsr r6, r9, #3
	and r6, r6, #0x1FFF1FFF
	uadd16 r7, r7, r3
	lsl r9, r5, #4
	and r9, r9, #0xFFF0FFF0
	uadd16 r3, r3, r9
	movw r9, #36409
	smulbb r10, r3, r9
	smultb r11, r3, r9
	pkhbt r9, r10, r11, lsl #16
	lsr r3, r9, #1
	and r3, r3, #0x7FFF7FFF
	eor r9, r9, r9
	usub16 r9, r9, r5
	usub16 r5, r9, r3
	mov r9, #30
	smulbb r10, r3, r9
	smultb r11, r3, r9
	pkhbt r9, r10, r11, lsl #16
	usub16 r7, r9, r7
	movw r9, #61167
	smulbb r10, r7, r9
	smultb r11, r7, r9
	pkhbt r9, r10, r11, lsl #16
	lsr r7, r9, #2
	and r7, r7, #0x3FFF3FFF
	usub16 r4, r4, r6
	usub16 r3, r3, r7

	ldr r9, [r1, #256]
	ldr r10, [r1, #128]
	ldr r11, [r1]
	uadd16 r9, r9, r6
	usub16 r9, r9, r2
	and r9, r9, #0x1FFF1FFF
	uadd16 r10, r10, r7
	usub16 r10, r10, r3
	and r10, r10, #0x1FFF1FFF
	uadd16 r11, r11, r8
	usub16 r11, r11, r4
	and r11, r11, #0x1FFF1FFF
	strh r9, [r1, #256]
	strh r10, [r1, #128]
	strh r11, [r1], #2

	ldr r10, [r1, #-130]
	ldr r11, [r1, #382]
	uadd16 r11, r11, r5
	usub16 r10, r10, r5
	and r10, r10, #0x1FFF1FFF
	and r11, r11, #0x1FFF1FFF
	cmp r14, #64
	bge greater_equal_64_round
	strh r11, [r1, #382]
	b cond_done_round
greater_equal_64_round:
	strh r10, [r1, #-130]
cond_done_round:

	ldmia.w sp!, {r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, pc}
.endfunc

.global call_nop
.func call_nop, call_nop
.type call_nop, %function
.align 2
call_nop:
	stmdb sp!, {r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, lr}
	ldmia.w sp!, {r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, pc}
.endfunc

