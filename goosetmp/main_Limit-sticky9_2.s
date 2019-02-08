	.file	"mainbody.bc"


	.section	.rodata.cst8,"aM",@progbits,8
	.align	8
.LCPI1_0:					
	.quad	4598175219545276416	# double value: 2.500000e-01
.LCPI1_1:					
	.quad	4615063718147915776	# double value: 3.500000e+00
.LCPI1_2:					
	.quad	4608308318706860032	# double value: 1.250000e+00
.LCPI1_3:					
	.quad	4616189618054758400	# double value: 4.000000e+00
.LCPI1_4:					
	.quad	4613937818241073152	# double value: 3.000000e+00
.LCPI1_5:					
	.quad	4609434218613702656	# double value: 1.500000e+00
.LCPI1_6:					
	.quad	13835058055282163712	# double value: -2.000000e+00
	.text
	.align	16
	.globl	main_Limit_2D_sticky9_2
	.type	main_Limit_2D_sticky9_2,@function
main_Limit_2D_sticky9_2:
.Leh_func_begin1:
.Llabel1:
	subq	$184, %rsp
	movsd	%xmm7, 88(%rsp)
	movsd	%xmm6, 80(%rsp)
	movsd	%xmm5, 16(%rsp)
	movsd	%xmm4, 168(%rsp)
	movsd	%xmm3, 72(%rsp)
	movsd	%xmm2, 160(%rsp)
	movsd	%xmm1, 136(%rsp)
	movsd	%xmm0, 104(%rsp)
	movsd	x_0_1__Limit_2D_sticky9_2, %xmm0
	movsd	%xmm0, 64(%rsp)
	subsd	%xmm2, %xmm0
	movsd	%xmm0, 144(%rsp)
	mulsd	%xmm0, %xmm0
	movsd	x_0_0__Limit_2D_sticky9_2, %xmm1
	movsd	%xmm1, 112(%rsp)
	subsd	104(%rsp), %xmm1
	movsd	%xmm1, 152(%rsp)
	mulsd	%xmm1, %xmm1
	addsd	%xmm0, %xmm1
	movsd	x_0_2__Limit_2D_sticky9_2, %xmm0
	movsd	%xmm0, 8(%rsp)
	subsd	%xmm4, %xmm0
	movsd	%xmm0, 56(%rsp)
	mulsd	%xmm0, %xmm0
	addsd	%xmm1, %xmm0
	movsd	%xmm0, 176(%rsp)
	addsd	eps2_Limit_2D_sticky9_2, %xmm0
	call	rsqrt
	movsd	%xmm0, 40(%rsp)
	movsd	64(%rsp), %xmm0
	subsd	192(%rsp), %xmm0
	movsd	%xmm0, 64(%rsp)
	mulsd	%xmm0, %xmm0
	movsd	112(%rsp), %xmm1
	subsd	80(%rsp), %xmm1
	movsd	%xmm1, 112(%rsp)
	mulsd	%xmm1, %xmm1
	addsd	%xmm0, %xmm1
	movsd	8(%rsp), %xmm0
	subsd	208(%rsp), %xmm0
	movsd	%xmm0, 8(%rsp)
	mulsd	%xmm0, %xmm0
	addsd	%xmm1, %xmm0
	addsd	eps2_Limit_2D_sticky9_2, %xmm0
	call	rsqrt
	movsd	%xmm0, 128(%rsp)
	movsd	192(%rsp), %xmm0
	subsd	160(%rsp), %xmm0
	movsd	%xmm0, 120(%rsp)
	mulsd	%xmm0, %xmm0
	movsd	80(%rsp), %xmm1
	subsd	104(%rsp), %xmm1
	movsd	%xmm1, 80(%rsp)
	mulsd	%xmm1, %xmm1
	addsd	%xmm0, %xmm1
	movsd	208(%rsp), %xmm0
	subsd	168(%rsp), %xmm0
	movsd	%xmm0, 48(%rsp)
	mulsd	%xmm0, %xmm0
	addsd	%xmm1, %xmm0
	movsd	%xmm0, 168(%rsp)
	addsd	eps2_Limit_2D_sticky9_2, %xmm0
	movsd	%xmm0, 24(%rsp)
	call	rsqrt
	movsd	m_0__Limit_2D_sticky9_2, %xmm1
	movsd	%xmm1, 160(%rsp)
	movsd	224(%rsp), %xmm1
	movsd	40(%rsp), %xmm2
	mulsd	%xmm1, %xmm2
	mulsd	40(%rsp), %xmm2
	mulsd	40(%rsp), %xmm2
	movsd	%xmm2, 104(%rsp)
	movsd	176(%rsp), %xmm2
	subsd	168(%rsp), %xmm2
	movsd	168(%rsp), %xmm3
	subsd	176(%rsp), %xmm3
	movsd	%xmm3, 168(%rsp)
	movsd	.LCPI1_0(%rip), %xmm3
	movsd	168(%rsp), %xmm4
	mulsd	%xmm3, %xmm4
	mulsd	128(%rsp), %xmm4
	mulsd	128(%rsp), %xmm4
	mulsd	128(%rsp), %xmm4
	movsd	%xmm4, 168(%rsp)
	mulsd	%xmm3, %xmm2
	mulsd	128(%rsp), %xmm2
	mulsd	128(%rsp), %xmm2
	mulsd	128(%rsp), %xmm2
	movapd	%xmm0, %xmm3
	mulsd	%xmm1, %xmm3
	mulsd	.LCPI1_1(%rip), %xmm1
	mulsd	160(%rsp), %xmm1
	mulsd	128(%rsp), %xmm1
	movsd	128(%rsp), %xmm4
	mulsd	.LCPI1_2(%rip), %xmm4
	movsd	40(%rsp), %xmm5
	mulsd	.LCPI1_3(%rip), %xmm5
	addsd	%xmm4, %xmm5
	addsd	168(%rsp), %xmm5
	mulsd	%xmm0, %xmm3
	mulsd	%xmm0, %xmm3
	movapd	%xmm3, %xmm6
	mulsd	160(%rsp), %xmm6
	mulsd	%xmm5, %xmm6
	movsd	%xmm6, 96(%rsp)
	movapd	%xmm6, %xmm5
	mulsd	80(%rsp), %xmm5
	movapd	%xmm0, %xmm6
	mulsd	.LCPI1_3(%rip), %xmm6
	addsd	%xmm4, %xmm6
	addsd	%xmm2, %xmm6
	movsd	160(%rsp), %xmm2
	mulsd	104(%rsp), %xmm2
	mulsd	%xmm6, %xmm2
	movsd	%xmm2, 160(%rsp)
	mulsd	152(%rsp), %xmm2
	addsd	%xmm5, %xmm2
	mulsd	128(%rsp), %xmm1
	mulsd	128(%rsp), %xmm1
	subsd	40(%rsp), %xmm0
	mulsd	%xmm1, %xmm0
	movsd	%xmm0, 40(%rsp)
	movsd	112(%rsp), %xmm0
	mulsd	40(%rsp), %xmm0
	movsd	%xmm0, 112(%rsp)
	subsd	%xmm0, %xmm2
	movsd	144(%rsp), %xmm0
	mulsd	72(%rsp), %xmm0
	movsd	72(%rsp), %xmm1
	mulsd	200(%rsp), %xmm1
	movsd	152(%rsp), %xmm4
	mulsd	136(%rsp), %xmm4
	movsd	%xmm4, 128(%rsp)
	movsd	136(%rsp), %xmm4
	mulsd	88(%rsp), %xmm4
	addsd	%xmm1, %xmm4
	movsd	56(%rsp), %xmm1
	mulsd	16(%rsp), %xmm1
	movsd	16(%rsp), %xmm5
	mulsd	216(%rsp), %xmm5
	addsd	%xmm4, %xmm5
	mulsd	.LCPI1_3(%rip), %xmm5
	movsd	%xmm5, 168(%rsp)
	movapd	%xmm5, %xmm4
	mulsd	152(%rsp), %xmm4
	movsd	144(%rsp), %xmm5
	mulsd	200(%rsp), %xmm5
	movsd	152(%rsp), %xmm6
	mulsd	88(%rsp), %xmm6
	addsd	%xmm5, %xmm6
	movsd	%xmm6, 152(%rsp)
	movsd	56(%rsp), %xmm5
	mulsd	216(%rsp), %xmm5
	addsd	%xmm6, %xmm5
	movsd	%xmm5, 112(%rsp)
	movsd	.LCPI1_4(%rip), %xmm5
	movsd	112(%rsp), %xmm6
	mulsd	%xmm5, %xmm6
	movsd	%xmm6, 112(%rsp)
	mulsd	136(%rsp), %xmm6
	subsd	%xmm6, %xmm4
	movsd	128(%rsp), %xmm6
	addsd	%xmm0, %xmm6
	addsd	%xmm1, %xmm6
	mulsd	.LCPI1_3(%rip), %xmm6
	movsd	%xmm6, 128(%rsp)
	movapd	%xmm6, %xmm0
	mulsd	88(%rsp), %xmm0
	subsd	%xmm0, %xmm4
	mulsd	104(%rsp), %xmm4
	subsd	%xmm4, %xmm2
	movsd	120(%rsp), %xmm0
	mulsd	200(%rsp), %xmm0
	movsd	80(%rsp), %xmm1
	mulsd	88(%rsp), %xmm1
	addsd	%xmm0, %xmm1
	movsd	48(%rsp), %xmm0
	mulsd	216(%rsp), %xmm0
	addsd	%xmm1, %xmm0
	movapd	%xmm0, %xmm1
	mulsd	.LCPI1_5(%rip), %xmm1
	mulsd	%xmm0, %xmm1
	divsd	24(%rsp), %xmm1
	movsd	200(%rsp), %xmm0
	subsd	72(%rsp), %xmm0
	movsd	%xmm0, 24(%rsp)
	movapd	%xmm0, %xmm4
	mulsd	%xmm4, %xmm4
	movsd	88(%rsp), %xmm0
	subsd	136(%rsp), %xmm0
	movapd	%xmm0, %xmm6
	mulsd	%xmm6, %xmm6
	addsd	%xmm4, %xmm6
	movsd	216(%rsp), %xmm4
	subsd	16(%rsp), %xmm4
	movsd	%xmm4, (%rsp)
	mulsd	%xmm4, %xmm4
	addsd	%xmm6, %xmm4
	mulsd	.LCPI1_6(%rip), %xmm4
	movsd	72(%rsp), %xmm6
	mulsd	%xmm6, %xmm6
	movsd	136(%rsp), %xmm7
	mulsd	%xmm7, %xmm7
	addsd	%xmm6, %xmm7
	movsd	16(%rsp), %xmm6
	mulsd	%xmm6, %xmm6
	addsd	%xmm7, %xmm6
	addsd	%xmm4, %xmm6
	addsd	%xmm1, %xmm6
	mulsd	%xmm3, %xmm6
	movsd	%xmm6, 152(%rsp)
	movapd	%xmm6, %xmm1
	mulsd	80(%rsp), %xmm1
	addsd	%xmm2, %xmm1
	movsd	200(%rsp), %xmm2
	mulsd	%xmm5, %xmm2
	movsd	72(%rsp), %xmm4
	mulsd	.LCPI1_3(%rip), %xmm4
	subsd	%xmm2, %xmm4
	mulsd	120(%rsp), %xmm4
	movsd	88(%rsp), %xmm2
	mulsd	%xmm5, %xmm2
	movsd	%xmm2, 88(%rsp)
	movsd	136(%rsp), %xmm2
	mulsd	.LCPI1_3(%rip), %xmm2
	subsd	88(%rsp), %xmm2
	mulsd	80(%rsp), %xmm2
	addsd	%xmm4, %xmm2
	movsd	%xmm2, 136(%rsp)
	mulsd	216(%rsp), %xmm5
	movsd	.LCPI1_3(%rip), %xmm2
	mulsd	16(%rsp), %xmm2
	subsd	%xmm5, %xmm2
	mulsd	48(%rsp), %xmm2
	addsd	136(%rsp), %xmm2
	mulsd	%xmm3, %xmm2
	movsd	%xmm2, 176(%rsp)
	mulsd	%xmm2, %xmm0
	addsd	%xmm1, %xmm0
	pxor	%xmm1, %xmm1
	addsd	%xmm1, %xmm0
	call	accum
	movsd	%xmm0, 136(%rsp)
	movsd	96(%rsp), %xmm0
	mulsd	120(%rsp), %xmm0
	movsd	160(%rsp), %xmm1
	mulsd	144(%rsp), %xmm1
	addsd	%xmm0, %xmm1
	movsd	64(%rsp), %xmm0
	mulsd	40(%rsp), %xmm0
	movsd	%xmm0, 64(%rsp)
	subsd	%xmm0, %xmm1
	movsd	144(%rsp), %xmm0
	mulsd	168(%rsp), %xmm0
	movsd	%xmm0, 144(%rsp)
	movsd	72(%rsp), %xmm0
	mulsd	112(%rsp), %xmm0
	movsd	%xmm0, 72(%rsp)
	movsd	144(%rsp), %xmm0
	subsd	72(%rsp), %xmm0
	movsd	%xmm0, 144(%rsp)
	movsd	200(%rsp), %xmm0
	mulsd	128(%rsp), %xmm0
	movsd	%xmm0, 32(%rsp)
	movsd	144(%rsp), %xmm0
	subsd	32(%rsp), %xmm0
	mulsd	104(%rsp), %xmm0
	movsd	%xmm0, 144(%rsp)
	subsd	%xmm0, %xmm1
	movsd	120(%rsp), %xmm0
	mulsd	152(%rsp), %xmm0
	addsd	%xmm1, %xmm0
	movsd	%xmm0, 120(%rsp)
	movsd	24(%rsp), %xmm0
	mulsd	176(%rsp), %xmm0
	addsd	120(%rsp), %xmm0
	movsd	%xmm0, 24(%rsp)
	pxor	%xmm1, %xmm1
	addsd	%xmm1, %xmm0
	call	accum
	movsd	%xmm0, 144(%rsp)
	movsd	96(%rsp), %xmm0
	mulsd	48(%rsp), %xmm0
	movsd	%xmm0, 96(%rsp)
	movsd	160(%rsp), %xmm0
	mulsd	56(%rsp), %xmm0
	addsd	96(%rsp), %xmm0
	movsd	%xmm0, 160(%rsp)
	movsd	40(%rsp), %xmm0
	mulsd	8(%rsp), %xmm0
	movsd	%xmm0, 40(%rsp)
	movsd	160(%rsp), %xmm0
	subsd	40(%rsp), %xmm0
	movsd	%xmm0, 160(%rsp)
	movsd	168(%rsp), %xmm0
	mulsd	56(%rsp), %xmm0
	movsd	%xmm0, 168(%rsp)
	movsd	112(%rsp), %xmm0
	mulsd	16(%rsp), %xmm0
	movsd	%xmm0, 112(%rsp)
	movsd	168(%rsp), %xmm0
	subsd	112(%rsp), %xmm0
	movsd	%xmm0, 168(%rsp)
	movsd	128(%rsp), %xmm0
	mulsd	216(%rsp), %xmm0
	movsd	%xmm0, 128(%rsp)
	movsd	168(%rsp), %xmm0
	subsd	128(%rsp), %xmm0
	mulsd	104(%rsp), %xmm0
	movsd	%xmm0, 168(%rsp)
	movsd	160(%rsp), %xmm0
	subsd	168(%rsp), %xmm0
	movsd	%xmm0, 160(%rsp)
	movsd	152(%rsp), %xmm0
	mulsd	48(%rsp), %xmm0
	addsd	160(%rsp), %xmm0
	movsd	%xmm0, 152(%rsp)
	movsd	176(%rsp), %xmm0
	mulsd	(%rsp), %xmm0
	addsd	152(%rsp), %xmm0
	movsd	%xmm0, 176(%rsp)
	pxor	%xmm1, %xmm1
	addsd	%xmm1, %xmm0
	call	accum
	movsd	136(%rsp), %xmm1
	movsd	%xmm1, ac_i_0__Limit_2D_sticky9_2
	movsd	144(%rsp), %xmm1
	movsd	%xmm1, ac_i_1__Limit_2D_sticky9_2
	movsd	%xmm0, ac_i_2__Limit_2D_sticky9_2
	pxor	%xmm0, %xmm0
	addq	$184, %rsp
	ret
	.size	main_Limit_2D_sticky9_2, .-main_Limit_2D_sticky9_2
.Leh_func_end1:
	.type	ac_i_0__Limit_2D_sticky9_2,@object
	.bss
	.globl ac_i_0__Limit_2D_sticky9_2
	.align	8
ac_i_0__Limit_2D_sticky9_2:				# ac_i_0__Limit-sticky9_2
	.size	ac_i_0__Limit_2D_sticky9_2, 8
	.zero	8
	.type	ac_i_1__Limit_2D_sticky9_2,@object
	.globl ac_i_1__Limit_2D_sticky9_2
	.align	8
ac_i_1__Limit_2D_sticky9_2:				# ac_i_1__Limit-sticky9_2
	.size	ac_i_1__Limit_2D_sticky9_2, 8
	.zero	8
	.type	ac_i_2__Limit_2D_sticky9_2,@object
	.globl ac_i_2__Limit_2D_sticky9_2
	.align	8
ac_i_2__Limit_2D_sticky9_2:				# ac_i_2__Limit-sticky9_2
	.size	ac_i_2__Limit_2D_sticky9_2, 8
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

.Lmain_Limit_2D_sticky9_2.eh:
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
