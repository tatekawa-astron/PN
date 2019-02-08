il_ps_2_0
dcl_input_interp(linear) v0.xy
;
;################################################################################
; generated on 2017-04-20 20:51:49
; generated from Limit-sticky9_1.q
; VARI  x_i_0_ x_i_1_ x_i_2_ m_i_ v_i_0_ v_i_1_ v_i_2_
; VARJ  x_j_0_ x_j_1_ x_j_2_ m_j_ v_j_0_ v_j_1_ v_j_2_
; VARF  a1c_0_ a1c_1_ a1c_2_
; CONST x_0_0_ x_0_1_ x_0_2_ eps2 m_0_
;################################################################################
;x_i_0_ -> res0
;x_i_1_ -> res1
;x_i_2_ -> res2
;m_i_ -> res3
;v_i_0_ -> res4
;v_i_1_ -> res5
;v_i_2_ -> res6
;x_j_0_ -> res7
;x_j_1_ -> res8
;x_j_2_ -> res9
;m_j_ -> res10
;v_j_0_ -> res11
;v_j_1_ -> res12
;v_j_2_ -> res13
;a1c_0_ -> o0
;a1c_1_ -> o1
;a1c_2_ -> o2
dcl_resource_id(0)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(1)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(2)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(3)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(4)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(5)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(6)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(7)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(8)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(9)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(10)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(11)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(12)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(13)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_output_generic o0
dcl_output_generic o1
dcl_output_generic o2
; const literal 0.000000e+00
dcl_literal l20, 0x0, 0x0, 0x0, 0x0
; const literal 2.500000e-01
dcl_literal l21, 0x0, 0x3fd00000, 0x0, 0x3fd00000
; const literal 1.250000e+00
dcl_literal l22, 0x0, 0x3ff40000, 0x0, 0x3ff40000
; const literal 3.000000e+00
dcl_literal l23, 0x0, 0x40080000, 0x0, 0x40080000
; const literal 3.500000e+00
dcl_literal l24, 0x0, 0x400c0000, 0x0, 0x400c0000
; const literal 4.000000e+00
dcl_literal l25, 0x0, 0x40100000, 0x0, 0x40100000
dcl_cb cb0[5]
; const x_0_0_ -> r506 cb0[0]
; const x_0_1_ -> r507 cb0[1]
; const x_0_2_ -> r508 cb0[2]
; const eps2 -> r509 cb0[3]
; const m_0_ -> r510 cb0[4]

; nj
dcl_cb cb10[1]

; float zero
dcl_literal l10, 0.0, 0.0, 0.0, 0.0
; for loop counter
dcl_literal l11, 0.0, 0.0, 1.0, 0x1
; for rsq
dcl_literal l4, 0x0, 0x3fe00000, 0x0, 0x3ff80000
; for negation
dcl_literal l8, 0x0, 0x80000000, 0x0, 0x80000000
;;
dcl_literal l50, 0x0, 0x0, 0x0, 0x0
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;load VARI
sample_resource(0)_sampler(0) r100, v0.xy
sample_resource(1)_sampler(1) r101, v0.xy
sample_resource(2)_sampler(2) r102, v0.xy
sample_resource(3)_sampler(3) r103, v0.xy
sample_resource(4)_sampler(4) r104, v0.xy
sample_resource(5)_sampler(5) r105, v0.xy
sample_resource(6)_sampler(6) r106, v0.xy
;load CB
mov r506, cb0[0]
mov r507, cb0[1]
mov r508, cb0[2]
mov r509, cb0[3]
mov r510, cb0[4]
mov r500, l20
mov r501, l21
mov r502, l22
mov r503, l23
mov r504, l24
mov r505, l25
; zero accum
mov r400, l10
mov r401, l10
mov r402, l10



mov r200.xy, l11.xy
mov r201.w, cb10[0].w
mov r201.z, cb10[0].z
ixor r201.x, r201.x, r201.x
itof r203.x, r201.z

whileloop
  sample_resource(7)_sampler(7) r300, r200.xy
  sample_resource(8)_sampler(8) r301, r200.xy
  sample_resource(9)_sampler(9) r302, r200.xy
  sample_resource(10)_sampler(10) r303, r200.xy
  sample_resource(11)_sampler(11) r304, r200.xy
  sample_resource(12)_sampler(12) r305, r200.xy
  sample_resource(13)_sampler(13) r306, r200.xy

  call 0
  call 1
  call 0


  add r200.x, r200.x, l11.z
  iadd r201.x, r201.x, l11.w
  ige r201.y, r201.x, r201.w
  break_logicalnz r201.y

;  umod r202.x, r201.x, r201.z
;  if_logicalz r202.x
;    add r200.0y, r200.0y, l11.z
;  endif

  eq r202.x, r200.x, r203.x
  if_logicalnz r202.x
    add r200.0y, r200.0y, l11.z
  endif
endloop

mov o0, r400
mov o1, r401
mov o2, r402


ret_dyn
endmain
; end of the main function

;; force
func 0
  dadd r799.xy, r100.xy, r506_neg(yw).xy
  dadd r799.zw, r100.zw, r506_neg(yw).zw
  dadd r811.xy, r300.xy, r506_neg(yw).xy
  dadd r811.zw, r300.zw, r506_neg(yw).zw
  dadd r812.xy, r300.xy, r100_neg(yw).xy
  dadd r812.zw, r300.zw, r100_neg(yw).zw
  dadd r814.xy, r101.xy, r507_neg(yw).xy
  dadd r814.zw, r101.zw, r507_neg(yw).zw
  dadd r816.xy, r301.xy, r507_neg(yw).xy
  dadd r816.zw, r301.zw, r507_neg(yw).zw
  dadd r818.xy, r301.xy, r101_neg(yw).xy
  dadd r818.zw, r301.zw, r101_neg(yw).zw
  dadd r800.xy, r102.xy, r508_neg(yw).xy
  dadd r800.zw, r102.zw, r508_neg(yw).zw
  dadd r802.xy, r302.xy, r508_neg(yw).xy
  dadd r802.zw, r302.zw, r508_neg(yw).zw
  dadd r804.xy, r302.xy, r102_neg(yw).xy
  dadd r804.zw, r302.zw, r102_neg(yw).zw
  dmul r740.xy, r799.xy, r799.xy
  dmul r740.zw, r799.zw, r799.zw
  dmul r747.xy, r814.xy, r814.xy
  dmul r747.zw, r814.zw, r814.zw
  dadd r703.xy, r740.xy, r747.xy
  dadd r703.zw, r740.zw, r747.zw
  dmul r749.xy, r800.xy, r800.xy
  dmul r749.zw, r800.zw, r800.zw
  dadd r705.xy, r703.xy, r749.xy
  dadd r705.zw, r703.zw, r749.zw
  dadd r706.xy, r705.xy, r509.xy
  dadd r706.zw, r705.zw, r509.zw
  mov r10.xy, r706.xy
  call 15
  mov r723.xy, r11.xy
  mov r10.xy, r706.zw
  call 15
  mov r723.__zw, r11.xyxy
  dmul r756.xy, r811.xy, r811.xy
  dmul r756.zw, r811.zw, r811.zw
  dmul r757.xy, r816.xy, r816.xy
  dmul r757.zw, r816.zw, r816.zw
  dadd r708.xy, r756.xy, r757.xy
  dadd r708.zw, r756.zw, r757.zw
  dmul r758.xy, r802.xy, r802.xy
  dmul r758.zw, r802.zw, r802.zw
  dadd r709.xy, r708.xy, r758.xy
  dadd r709.zw, r708.zw, r758.zw
  dadd r710.xy, r709.xy, r509.xy
  dadd r710.zw, r709.zw, r509.zw
  mov r10.xy, r710.xy
  call 15
  mov r724.xy, r11.xy
  mov r10.xy, r710.zw
  call 15
  mov r724.__zw, r11.xyxy
  dmul r759.xy, r812.xy, r812.xy
  dmul r759.zw, r812.zw, r812.zw
  dmul r760.xy, r818.xy, r818.xy
  dmul r760.zw, r818.zw, r818.zw
  dadd r711.xy, r759.xy, r760.xy
  dadd r711.zw, r759.zw, r760.zw
  dmul r761.xy, r804.xy, r804.xy
  dmul r761.zw, r804.zw, r804.zw
  dadd r712.xy, r711.xy, r761.xy
  dadd r712.zw, r711.zw, r761.zw
  dadd r713.xy, r712.xy, r509.xy
  dadd r713.zw, r712.zw, r509.zw
  mov r10.xy, r713.xy
  call 15
  mov r725.xy, r11.xy
  mov r10.xy, r713.zw
  call 15
  mov r725.__zw, r11.xyxy
  dmul r762.xy, r723.xy, r103.xy
  dmul r762.zw, r723.zw, r103.zw
  dmul r763.xy, r762.xy, r723.xy
  dmul r763.zw, r762.zw, r723.zw
  dmul r764.xy, r763.xy, r723.xy
  dmul r764.zw, r763.zw, r723.zw
  dmul r765.xy, r764.xy, r303.xy
  dmul r765.zw, r764.zw, r303.zw
  dmul r766.xy, r765.xy, r799.xy
  dmul r766.zw, r765.zw, r799.zw
  dmul r767.xy, r724.xy, r505.xy
  dmul r767.zw, r724.zw, r505.zw
  dmul r768.xy, r725.xy, r502.xy
  dmul r768.zw, r725.zw, r502.zw
  dadd r714.xy, r767.xy, r768.xy
  dadd r714.zw, r767.zw, r768.zw
  dmul r769.xy, r709.xy, r501.xy
  dmul r769.zw, r709.zw, r501.zw
  ddiv r737.xy, r769.xy, r713.xy
  ddiv r737.zw, r769.zw, r713.zw
  dmul r770.xy, r737.xy, r725.xy
  dmul r770.zw, r737.zw, r725.zw
  dadd r813.xy, r714.xy, r770_neg(yw).xy
  dadd r813.zw, r714.zw, r770_neg(yw).zw
  dmul r771.xy, r705.xy, r501.xy
  dmul r771.zw, r705.zw, r501.zw
  ddiv r738.xy, r771.xy, r713.xy
  ddiv r738.zw, r771.zw, r713.zw
  dmul r772.xy, r738.xy, r725.xy
  dmul r772.zw, r738.zw, r725.zw
  dadd r715.xy, r813.xy, r772.xy
  dadd r715.zw, r813.zw, r772.zw
  dmul r773.xy, r766.xy, r715.xy
  dmul r773.zw, r766.zw, r715.zw
  dmul r774.xy, r725.xy, r725.xy
  dmul r774.zw, r725.zw, r725.zw
  dmul r775.xy, r774.xy, r725.xy
  dmul r775.zw, r774.zw, r725.zw
  dmul r776.xy, r775.xy, r504.xy
  dmul r776.zw, r775.zw, r504.zw
  dmul r777.xy, r776.xy, r723.xy
  dmul r777.zw, r776.zw, r723.zw
  dmul r778.xy, r777.xy, r103.xy
  dmul r778.zw, r777.zw, r103.zw
  dmul r779.xy, r778.xy, r303.xy
  dmul r779.zw, r778.zw, r303.zw
  dmul r780.xy, r779.xy, r812.xy
  dmul r780.zw, r779.zw, r812.zw
  dadd r815.xy, r773.xy, r780_neg(yw).xy
  dadd r815.zw, r773.zw, r780_neg(yw).zw
  ddiv r739.xy, r765.xy, r510.xy
  ddiv r739.zw, r765.zw, r510.zw
  dmul r781.xy, r104.xy, r304.xy
  dmul r781.zw, r104.zw, r304.zw
  dmul r782.xy, r105.xy, r305.xy
  dmul r782.zw, r105.zw, r305.zw
  dadd r716.xy, r781.xy, r782.xy
  dadd r716.zw, r781.zw, r782.zw
  dmul r783.xy, r106.xy, r306.xy
  dmul r783.zw, r106.zw, r306.zw
  dadd r717.xy, r716.xy, r783.xy
  dadd r717.zw, r716.zw, r783.zw
  dmul r784.xy, r717.xy, r505.xy
  dmul r784.zw, r717.zw, r505.zw
  dmul r785.xy, r784.xy, r799.xy
  dmul r785.zw, r784.zw, r799.zw
  dmul r786.xy, r799.xy, r104.xy
  dmul r786.zw, r799.zw, r104.zw
  dmul r787.xy, r814.xy, r105.xy
  dmul r787.zw, r814.zw, r105.zw
  dadd r718.xy, r786.xy, r787.xy
  dadd r718.zw, r786.zw, r787.zw
  dmul r788.xy, r800.xy, r106.xy
  dmul r788.zw, r800.zw, r106.zw
  dadd r719.xy, r718.xy, r788.xy
  dadd r719.zw, r718.zw, r788.zw
  dmul r789.xy, r719.xy, r503.xy
  dmul r789.zw, r719.zw, r503.zw
  dmul r790.xy, r789.xy, r304.xy
  dmul r790.zw, r789.zw, r304.zw
  dadd r817.xy, r785.xy, r790_neg(yw).xy
  dadd r817.zw, r785.zw, r790_neg(yw).zw
  dmul r791.xy, r799.xy, r304.xy
  dmul r791.zw, r799.zw, r304.zw
  dmul r792.xy, r814.xy, r305.xy
  dmul r792.zw, r814.zw, r305.zw
  dadd r720.xy, r791.xy, r792.xy
  dadd r720.zw, r791.zw, r792.zw
  dmul r793.xy, r800.xy, r306.xy
  dmul r793.zw, r800.zw, r306.zw
  dadd r721.xy, r720.xy, r793.xy
  dadd r721.zw, r720.zw, r793.zw
  dmul r794.xy, r721.xy, r505.xy
  dmul r794.zw, r721.zw, r505.zw
  dmul r795.xy, r794.xy, r104.xy
  dmul r795.zw, r794.zw, r104.zw
  dadd r819.xy, r817.xy, r795_neg(yw).xy
  dadd r819.zw, r817.zw, r795_neg(yw).zw
  dmul r796.xy, r739.xy, r819.xy
  dmul r796.zw, r739.zw, r819.zw
  dadd r820.xy, r815.xy, r796_neg(yw).xy
  dadd r820.zw, r815.zw, r796_neg(yw).zw
  dadd r722.xy, r820.xy, r500.xy
  dadd r722.zw, r820.zw, r500.zw
  dmul r797.xy, r765.xy, r814.xy
  dmul r797.zw, r765.zw, r814.zw
  dmul r741.xy, r797.xy, r715.xy
  dmul r741.zw, r797.zw, r715.zw
  dmul r742.xy, r779.xy, r818.xy
  dmul r742.zw, r779.zw, r818.zw
  dadd r801.xy, r741.xy, r742_neg(yw).xy
  dadd r801.zw, r741.zw, r742_neg(yw).zw
  dmul r743.xy, r784.xy, r814.xy
  dmul r743.zw, r784.zw, r814.zw
  dmul r744.xy, r789.xy, r305.xy
  dmul r744.zw, r789.zw, r305.zw
  dadd r803.xy, r743.xy, r744_neg(yw).xy
  dadd r803.zw, r743.zw, r744_neg(yw).zw
  dmul r745.xy, r794.xy, r105.xy
  dmul r745.zw, r794.zw, r105.zw
  dadd r805.xy, r803.xy, r745_neg(yw).xy
  dadd r805.zw, r803.zw, r745_neg(yw).zw
  dmul r746.xy, r739.xy, r805.xy
  dmul r746.zw, r739.zw, r805.zw
  dadd r806.xy, r801.xy, r746_neg(yw).xy
  dadd r806.zw, r801.zw, r746_neg(yw).zw
  dadd r704.xy, r806.xy, r500.xy
  dadd r704.zw, r806.zw, r500.zw
  dmul r748.xy, r765.xy, r800.xy
  dmul r748.zw, r765.zw, r800.zw
  dmul r750.xy, r748.xy, r715.xy
  dmul r750.zw, r748.zw, r715.zw
  dmul r751.xy, r779.xy, r804.xy
  dmul r751.zw, r779.zw, r804.zw
  dadd r807.xy, r750.xy, r751_neg(yw).xy
  dadd r807.zw, r750.zw, r751_neg(yw).zw
  dmul r752.xy, r784.xy, r800.xy
  dmul r752.zw, r784.zw, r800.zw
  dmul r753.xy, r789.xy, r306.xy
  dmul r753.zw, r789.zw, r306.zw
  dadd r808.xy, r752.xy, r753_neg(yw).xy
  dadd r808.zw, r752.zw, r753_neg(yw).zw
  dmul r754.xy, r794.xy, r106.xy
  dmul r754.zw, r794.zw, r106.zw
  dadd r809.xy, r808.xy, r754_neg(yw).xy
  dadd r809.zw, r808.zw, r754_neg(yw).zw
  dmul r755.xy, r739.xy, r809.xy
  dmul r755.zw, r739.zw, r809.zw
  dadd r810.xy, r807.xy, r755_neg(yw).xy
  dadd r810.zw, r807.zw, r755_neg(yw).zw
  dadd r707.xy, r810.xy, r500.xy
  dadd r707.zw, r810.zw, r500.zw
  dadd r400.xy, r400.xy, r722.xy
  dadd r400.zw, r400.zw, r722.zw
  dadd r401.xy, r401.xy, r704.xy
  dadd r401.zw, r401.zw, r704.zw
  dadd r402.xy, r402.xy, r707.xy
  dadd r402.zw, r402.zw, r707.zw

  ret
endfunc

;; permutation
func 1
  mov r300, r300.zwxy
  mov r301, r301.zwxy
  mov r302, r302.zwxy
  mov r303, r303.zwxy
  mov r304, r304.zwxy
  mov r305, r305.zwxy
  mov r306, r306.zwxy

  ret
endfunc

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; rsq double
; in r10.xy
; out r11.xy
func 15
  ; x0 = rsq(x)
  d2f r11, r10.xy
  rsq r11, r11
  f2d r12, r11

  ; x*0.5
  dmul r13.xy, r10.xy, l4.xy

  ; b = 1.5*x0
  dmul r14.xy, r12.xy, l4.zw

  ; a = x0*x0*x0*x*0.5
  dmul r16.xy, r12.xy, r12.xy
  dmul r16.xy, r16.xy, r12.xy
  dmul r16.xy, r16.xy, r13.xy

  ; x1 = b - a
  ixor r999, r16, l8
  dadd r11.xy, r14.xy, r999.xy
  ret
endfunc

end

