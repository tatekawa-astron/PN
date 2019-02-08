	.file	"mainbody.bc"


	.text
	.align	16
	.globl	main_Limit_2D_sticky9_0
	.type	main_Limit_2D_sticky9_0,@function
main_Limit_2D_sticky9_0:
.Leh_func_begin1:
.Llabel1:
	subq	$40, %rsp
	movsd	%xmm6, 32(%rsp)
	subsd	%xmm1, %xmm4
	movsd	%xmm4, 24(%rsp)
	movapd	%xmm4, %xmm1
	mulsd	%xmm1, %xmm1
	subsd	%xmm0, %xmm3
	movsd	%xmm3, 16(%rsp)
	mulsd	%xmm3, %xmm3
	addsd	%xmm1, %xmm3
	movapd	%xmm5, %xmm0
	subsd	%xmm2, %xmm0
	movsd	%xmm0, 8(%rsp)
	mulsd	%xmm0, %xmm0
	addsd	%xmm3, %xmm0
	addsd	eps2_Limit_2D_sticky9_0, %xmm0
	call	rsqrt
	movsd	32(%rsp), %xmm1
	mulsd	%xmm0, %xmm1
	mulsd	%xmm0, %xmm1
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, 32(%rsp)
	movsd	16(%rsp), %xmm0
	mulsd	%xmm1, %xmm0
	movsd	%xmm0, 16(%rsp)
	pxor	%xmm1, %xmm1
	addsd	%xmm1, %xmm0
	call	accum
	movsd	%xmm0, 16(%rsp)
	movsd	24(%rsp), %xmm0
	mulsd	32(%rsp), %xmm0
	movsd	%xmm0, 24(%rsp)
	pxor	%xmm1, %xmm1
	addsd	%xmm1, %xmm0
	call	accum
	movsd	%xmm0, 24(%rsp)
	movsd	32(%rsp), %xmm0
	mulsd	8(%rsp), %xmm0
	movsd	%xmm0, 32(%rsp)
	pxor	%xmm1, %xmm1
	addsd	%xmm1, %xmm0
	call	accum
	movsd	16(%rsp), %xmm1
	movsd	%xmm1, a_i_0__Limit_2D_sticky9_0
	movsd	24(%rsp), %xmm1
	movsd	%xmm1, a_i_1__Limit_2D_sticky9_0
	movsd	%xmm0, a_i_2__Limit_2D_sticky9_0
	pxor	%xmm0, %xmm0
	addq	$40, %rsp
	ret
	.size	main_Limit_2D_sticky9_0, .-main_Limit_2D_sticky9_0
.Leh_func_end1:
	.type	a_i_0__Limit_2D_sticky9_0,@object
	.bss
	.globl a_i_0__Limit_2D_sticky9_0
	.align	8
a_i_0__Limit_2D_sticky9_0:				# a_i_0__Limit-sticky9_0
	.size	a_i_0__Limit_2D_sticky9_0, 8
	.zero	8
	.type	a_i_1__Limit_2D_sticky9_0,@object
	.globl a_i_1__Limit_2D_sticky9_0
	.align	8
a_i_1__Limit_2D_sticky9_0:				# a_i_1__Limit-sticky9_0
	.size	a_i_1__Limit_2D_sticky9_0, 8
	.zero	8
	.type	a_i_2__Limit_2D_sticky9_0,@object
	.globl a_i_2__Limit_2D_sticky9_0
	.align	8
a_i_2__Limit_2D_sticky9_0:				# a_i_2__Limit-sticky9_0
	.size	a_i_2__Limit_2D_sticky9_0, 8
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

.Lmain_Limit_2D_sticky9_0.eh:
	.long	.Leh_frame_end1-.Leh_frame_begin1
.Leh_frame_begin1:
	.long	.Leh_frame_begin1-.Leh_frame_common
	.long	.Leh_func_begin1-.
	.long	.Leh_func_end1-.Leh_func_begin1
	.uleb128	0
	.byte	0xE
	.uleb128	48
	.byte	0x4
	.long	.Llabel1-.Leh_func_begin1
	.byte	0xD
	.uleb128	7
	.align	8
.Leh_frame_end1:

	.section	.note.GNU-stack,"",@progbits
