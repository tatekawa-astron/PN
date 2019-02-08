	.file	"mainbody.bc"


	.section	.rodata.cst8,"aM",@progbits,8
	.align	8
.LCPI1_0:					
	.quad	4613937818241073152	# double value: 3.000000e+00
.LCPI1_1:					
	.quad	4616189618054758400	# double value: 4.000000e+00
.LCPI1_2:					
	.quad	4598175219545276416	# double value: 2.500000e-01
.LCPI1_3:					
	.quad	4608308318706860032	# double value: 1.250000e+00
.LCPI1_4:					
	.quad	4615063718147915776	# double value: 3.500000e+00
	.text
	.align	16
	.globl	main_Limit_2D_sticky9_1
	.type	main_Limit_2D_sticky9_1,@function
main_Limit_2D_sticky9_1:
.Leh_func_begin1:
.Llabel1:
	subq	$184, %rsp
	movsd	%xmm7, 96(%rsp)
	movsd	%xmm6, 24(%rsp)
	movsd	%xmm5, 72(%rsp)
	movsd	%xmm4, 64(%rsp)
	movsd	%xmm3, 8(%rsp)
	movsd	%xmm2, 176(%rsp)
	movsd	%xmm1, 168(%rsp)
	movsd	%xmm0, 128(%rsp)
	movsd	x_0_1__Limit_2D_sticky9_1, %xmm0
	movsd	%xmm0, 136(%rsp)
	movapd	%xmm1, %xmm0
	subsd	136(%rsp), %xmm0
	movsd	%xmm0, 144(%rsp)
	mulsd	%xmm0, %xmm0
	movsd	x_0_0__Limit_2D_sticky9_1, %xmm1
	movsd	%xmm1, 120(%rsp)
	movsd	128(%rsp), %xmm1
	subsd	120(%rsp), %xmm1
	movsd	%xmm1, 160(%rsp)
	mulsd	%xmm1, %xmm1
	addsd	%xmm0, %xmm1
	movsd	x_0_2__Limit_2D_sticky9_1, %xmm0
	movsd	%xmm0, 112(%rsp)
	movapd	%xmm2, %xmm0
	subsd	112(%rsp), %xmm0
	movsd	%xmm0, 88(%rsp)
	mulsd	%xmm0, %xmm0
	addsd	%xmm1, %xmm0
	movsd	%xmm0, 152(%rsp)
	addsd	eps2_Limit_2D_sticky9_1, %xmm0
	call	rsqrt
	movsd	%xmm0, 48(%rsp)
	movsd	192(%rsp), %xmm0
	subsd	136(%rsp), %xmm0
	mulsd	%xmm0, %xmm0
	movsd	96(%rsp), %xmm1
	subsd	120(%rsp), %xmm1
	mulsd	%xmm1, %xmm1
	addsd	%xmm0, %xmm1
	movsd	200(%rsp), %xmm0
	subsd	112(%rsp), %xmm0
	mulsd	%xmm0, %xmm0
	addsd	%xmm1, %xmm0
	movsd	%xmm0, 136(%rsp)
	addsd	eps2_Limit_2D_sticky9_1, %xmm0
	call	rsqrt
	movsd	%xmm0, 16(%rsp)
	movsd	192(%rsp), %xmm0
	subsd	168(%rsp), %xmm0
	movsd	%xmm0, 80(%rsp)
	mulsd	%xmm0, %xmm0
	movsd	96(%rsp), %xmm1
	subsd	128(%rsp), %xmm1
	movsd	%xmm1, 96(%rsp)
	mulsd	%xmm1, %xmm1
	addsd	%xmm0, %xmm1
	movsd	200(%rsp), %xmm0
	subsd	176(%rsp), %xmm0
	movsd	%xmm0, 40(%rsp)
	mulsd	%xmm0, %xmm0
	addsd	%xmm1, %xmm0
	addsd	eps2_Limit_2D_sticky9_1, %xmm0
	movsd	%xmm0, 112(%rsp)
	call	rsqrt
	movsd	48(%rsp), %xmm1
	mulsd	8(%rsp), %xmm1
	mulsd	48(%rsp), %xmm1
	mulsd	48(%rsp), %xmm1
	movsd	%xmm1, 176(%rsp)
	movsd	208(%rsp), %xmm1
	movsd	176(%rsp), %xmm2
	mulsd	%xmm1, %xmm2
	movsd	%xmm2, 176(%rsp)
	mulsd	160(%rsp), %xmm2
	movsd	176(%rsp), %xmm3
	mulsd	144(%rsp), %xmm3
	movsd	%xmm3, 104(%rsp)
	movsd	176(%rsp), %xmm3
	divsd	m_0__Limit_2D_sticky9_1, %xmm3
	movsd	%xmm3, 56(%rsp)
	movsd	144(%rsp), %xmm3
	mulsd	72(%rsp), %xmm3
	movsd	160(%rsp), %xmm4
	mulsd	64(%rsp), %xmm4
	addsd	%xmm3, %xmm4
	movsd	88(%rsp), %xmm3
	mulsd	24(%rsp), %xmm3
	addsd	%xmm4, %xmm3
	mulsd	.LCPI1_0(%rip), %xmm3
	movsd	%xmm3, 128(%rsp)
	movsd	216(%rsp), %xmm3
	movsd	128(%rsp), %xmm4
	mulsd	%xmm3, %xmm4
	movsd	72(%rsp), %xmm5
	mulsd	224(%rsp), %xmm5
	movsd	64(%rsp), %xmm6
	mulsd	%xmm3, %xmm6
	addsd	%xmm5, %xmm6
	movsd	24(%rsp), %xmm5
	mulsd	232(%rsp), %xmm5
	addsd	%xmm6, %xmm5
	movsd	%xmm5, 168(%rsp)
	movsd	.LCPI1_1(%rip), %xmm5
	movsd	168(%rsp), %xmm6
	mulsd	%xmm5, %xmm6
	movsd	%xmm6, 168(%rsp)
	mulsd	160(%rsp), %xmm6
	subsd	%xmm4, %xmm6
	movsd	144(%rsp), %xmm4
	mulsd	224(%rsp), %xmm4
	movsd	160(%rsp), %xmm7
	mulsd	%xmm3, %xmm7
	addsd	%xmm4, %xmm7
	movsd	%xmm7, 160(%rsp)
	movsd	88(%rsp), %xmm3
	mulsd	232(%rsp), %xmm3
	addsd	%xmm7, %xmm3
	mulsd	%xmm5, %xmm3
	movsd	%xmm3, 120(%rsp)
	movsd	64(%rsp), %xmm3
	mulsd	120(%rsp), %xmm3
	movsd	%xmm3, 64(%rsp)
	subsd	%xmm3, %xmm6
	mulsd	56(%rsp), %xmm6
	movsd	.LCPI1_2(%rip), %xmm3
	movsd	136(%rsp), %xmm4
	mulsd	%xmm3, %xmm4
	divsd	112(%rsp), %xmm4
	mulsd	%xmm0, %xmm4
	movsd	%xmm4, 136(%rsp)
	movapd	%xmm0, %xmm4
	mulsd	%xmm4, %xmm4
	mulsd	%xmm0, %xmm4
	movsd	%xmm4, 160(%rsp)
	movsd	16(%rsp), %xmm4
	mulsd	%xmm5, %xmm4
	movsd	%xmm4, 16(%rsp)
	movapd	%xmm0, %xmm4
	mulsd	.LCPI1_3(%rip), %xmm4
	addsd	16(%rsp), %xmm4
	subsd	136(%rsp), %xmm4
	movsd	152(%rsp), %xmm5
	mulsd	%xmm3, %xmm5
	divsd	112(%rsp), %xmm5
	mulsd	%xmm0, %xmm5
	addsd	%xmm4, %xmm5
	movsd	%xmm5, 152(%rsp)
	mulsd	%xmm5, %xmm2
	movsd	160(%rsp), %xmm0
	mulsd	.LCPI1_4(%rip), %xmm0
	mulsd	48(%rsp), %xmm0
	mulsd	8(%rsp), %xmm0
	mulsd	%xmm1, %xmm0
	movsd	%xmm0, 160(%rsp)
	movsd	96(%rsp), %xmm0
	mulsd	160(%rsp), %xmm0
	movsd	%xmm0, 96(%rsp)
	subsd	%xmm0, %xmm2
	subsd	%xmm6, %xmm2
	movapd	%xmm2, %xmm0
	pxor	%xmm1, %xmm1
	addsd	%xmm1, %xmm0
	call	accum
	movsd	%xmm0, 136(%rsp)
	movsd	224(%rsp), %xmm0
	mulsd	128(%rsp), %xmm0
	movsd	%xmm0, 32(%rsp)
	movsd	144(%rsp), %xmm0
	mulsd	168(%rsp), %xmm0
	subsd	32(%rsp), %xmm0
	movsd	%xmm0, 144(%rsp)
	movsd	72(%rsp), %xmm0
	mulsd	120(%rsp), %xmm0
	movsd	%xmm0, 72(%rsp)
	movsd	144(%rsp), %xmm0
	subsd	72(%rsp), %xmm0
	mulsd	56(%rsp), %xmm0
	movsd	%xmm0, 144(%rsp)
	movsd	104(%rsp), %xmm0
	mulsd	152(%rsp), %xmm0
	movsd	%xmm0, 104(%rsp)
	movsd	80(%rsp), %xmm0
	mulsd	160(%rsp), %xmm0
	movsd	%xmm0, 80(%rsp)
	movsd	104(%rsp), %xmm0
	subsd	80(%rsp), %xmm0
	subsd	144(%rsp), %xmm0
	movsd	%xmm0, 104(%rsp)
	pxor	%xmm1, %xmm1
	addsd	%xmm1, %xmm0
	call	accum
	movsd	%xmm0, 144(%rsp)
	movsd	128(%rsp), %xmm0
	mulsd	232(%rsp), %xmm0
	movsd	%xmm0, 128(%rsp)
	movsd	168(%rsp), %xmm0
	mulsd	88(%rsp), %xmm0
	subsd	128(%rsp), %xmm0
	movsd	%xmm0, 168(%rsp)
	movsd	120(%rsp), %xmm0
	mulsd	24(%rsp), %xmm0
	movsd	%xmm0, 120(%rsp)
	movsd	168(%rsp), %xmm0
	subsd	120(%rsp), %xmm0
	mulsd	56(%rsp), %xmm0
	movsd	%xmm0, 168(%rsp)
	movsd	176(%rsp), %xmm0
	mulsd	88(%rsp), %xmm0
	mulsd	152(%rsp), %xmm0
	movsd	%xmm0, 176(%rsp)
	movsd	160(%rsp), %xmm0
	mulsd	40(%rsp), %xmm0
	movsd	%xmm0, 160(%rsp)
	movsd	176(%rsp), %xmm0
	subsd	160(%rsp), %xmm0
	subsd	168(%rsp), %xmm0
	movsd	%xmm0, 176(%rsp)
	pxor	%xmm1, %xmm1
	addsd	%xmm1, %xmm0
	call	accum
	movsd	136(%rsp), %xmm1
	movsd	%xmm1, a1c_0__Limit_2D_sticky9_1
	movsd	144(%rsp), %xmm1
	movsd	%xmm1, a1c_1__Limit_2D_sticky9_1
	movsd	%xmm0, a1c_2__Limit_2D_sticky9_1
	pxor	%xmm0, %xmm0
	addq	$184, %rsp
	ret
	.size	main_Limit_2D_sticky9_1, .-main_Limit_2D_sticky9_1
.Leh_func_end1:
	.type	a1c_0__Limit_2D_sticky9_1,@object
	.bss
	.globl a1c_0__Limit_2D_sticky9_1
	.align	8
a1c_0__Limit_2D_sticky9_1:				# a1c_0__Limit-sticky9_1
	.size	a1c_0__Limit_2D_sticky9_1, 8
	.zero	8
	.type	a1c_1__Limit_2D_sticky9_1,@object
	.globl a1c_1__Limit_2D_sticky9_1
	.align	8
a1c_1__Limit_2D_sticky9_1:				# a1c_1__Limit-sticky9_1
	.size	a1c_1__Limit_2D_sticky9_1, 8
	.zero	8
	.type	a1c_2__Limit_2D_sticky9_1,@object
	.globl a1c_2__Limit_2D_sticky9_1
	.align	8
a1c_2__Limit_2D_sticky9_1:				# a1c_2__Limit-sticky9_1
	.size	a1c_2__Limit_2D_sticky9_1, 8
	.zero	8
	.section	.eh_frame,"aw",@progbits
.LEH_frame0:
.Lsection_eh_frame:
.Leh_frame_common:
	.long	.Leh_frame_common_end-.Leh_frame_common_begin
.Leh_frame_common_begin:
	.long	0x0
	.byte	0x1
	.asciz	"zR"
	.uleb128	1
	.sleb128	-8
	.byte	0x10
	.uleb128	1
	.byte	0x1B
	.byte	0xC
	.uleb128	7
	.uleb128	8
	.byte	0x90
	.uleb128	1
	.align	8
.Leh_frame_common_end:

.Lmain_Limit_2D_sticky9_1.eh:
	.long	.Leh_frame_end1-.Leh_frame_begin1
.Leh_frame_begin1:
	.long	.Leh_frame_begin1-.Leh_frame_common
	.long	.Leh_func_begin1-.
	.long	.Leh_func_end1-.Leh_func_begin1
	.uleb128	0
	.byte	0xE
	.uleb128	192
	.byte	0x4
	.long	.Llabel1-.Leh_func_begin1
	.byte	0xD
	.uleb128	7
	.align	8
.Leh_frame_end1:

	.section	.note.GNU-stack,"",@progbits
