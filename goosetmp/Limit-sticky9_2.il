il_ps_2_0
dcl_input_interp(linear) v0.xy
;
;################################################################################
; generated on 2017-04-20 20:51:49
; generated from Limit-sticky9_2.q
; VARI  x_i_0_ v_i_0_ x_i_1_ v_i_1_ x_i_2_ v_i_2_
; VARJ  x_j_0_ v_j_0_ x_j_1_ v_j_1_ x_j_2_ v_j_2_ m_j_
; VARF  ac_i_0_ ac_i_1_ ac_i_2_
; CONST x_0_0_ x_0_1_ x_0_2_ eps2 m_0_
;################################################################################
;x_i_0_ -> res0
;v_i_0_ -> res1
;x_i_1_ -> res2
;v_i_1_ -> res3
;x_i_2_ -> res4
;v_i_2_ -> res5
;x_j_0_ -> res6
;v_j_0_ -> res7
;x_j_1_ -> res8
;v_j_1_ -> res9
;x_j_2_ -> res10
;v_j_2_ -> res11
;m_j_ -> res12
;ac_i_0_ -> o0
;ac_i_1_ -> o1
;ac_i_2_ -> o2
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
dcl_output_generic o0
dcl_output_generic o1
dcl_output_generic o2
; const literal 0.000000e+00
dcl_literal l20, 0x0, 0x0, 0x0, 0x0
; const literal 2.500000e-01
dcl_literal l21, 0x0, 0x3fd00000, 0x0, 0x3fd00000
; const literal 1.250000e+00
dcl_literal l22, 0x0, 0x3ff40000, 0x0, 0x3ff40000
; const literal 1.500000e+00
dcl_literal l23, 0x0, 0x3ff80000, 0x0, 0x3ff80000
; const literal 2.000000e+00
dcl_literal l24, 0x0, 0x40000000, 0x0, 0x40000000
; const literal 3.000000e+00
dcl_literal l25, 0x0, 0x40080000, 0x0, 0x40080000
; const literal 3.500000e+00
dcl_literal l26, 0x0, 0x400c0000, 0x0, 0x400c0000
; const literal 4.000000e+00
dcl_literal l27, 0x0, 0x40100000, 0x0, 0x40100000
dcl_cb cb0[5]
; const x_0_0_ -> r508 cb0[0]
; const x_0_1_ -> r509 cb0[1]
; const x_0_2_ -> r510 cb0[2]
; const eps2 -> r511 cb0[3]
; const m_0_ -> r512 cb0[4]

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
;load CB
mov r508, cb0[0]
mov r509, cb0[1]
mov r510, cb0[2]
mov r511, cb0[3]
mov r512, cb0[4]
mov r500, l20
mov r501, l21
mov r502, l22
mov r503, l23
mov r504, l24
mov r505, l25
mov r506, l26
mov r507, l27
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
  sample_resource(6)_sampler(6) r300, r200.xy
  sample_resource(7)_sampler(7) r301, r200.xy
  sample_resource(8)_sampler(8) r302, r200.xy
  sample_resource(9)_sampler(9) r303, r200.xy
  sample_resource(10)_sampler(10) r304, r200.xy
  sample_resource(11)_sampler(11) r305, r200.xy
  sample_resource(12)_sampler(12) r306, r200.xy

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
  dadd r867.xy, r508.xy, r100_neg(yw).xy
  dadd r867.zw, r508.zw, r100_neg(yw).zw
  dadd r880.xy, r508.xy, r300_neg(yw).xy
  dadd r880.zw, r508.zw, r300_neg(yw).zw
  dadd r884.xy, r300.xy, r100_neg(yw).xy
  dadd r884.zw, r300.zw, r100_neg(yw).zw
  dadd r889.xy, r301.xy, r101_neg(yw).xy
  dadd r889.zw, r301.zw, r101_neg(yw).zw
  dadd r891.xy, r509.xy, r102_neg(yw).xy
  dadd r891.zw, r509.zw, r102_neg(yw).zw
  dadd r893.xy, r509.xy, r302_neg(yw).xy
  dadd r893.zw, r509.zw, r302_neg(yw).zw
  dadd r895.xy, r302.xy, r102_neg(yw).xy
  dadd r895.zw, r302.zw, r102_neg(yw).zw
  dadd r868.xy, r303.xy, r103_neg(yw).xy
  dadd r868.zw, r303.zw, r103_neg(yw).zw
  dadd r871.xy, r510.xy, r104_neg(yw).xy
  dadd r871.zw, r510.zw, r104_neg(yw).zw
  dadd r874.xy, r510.xy, r304_neg(yw).xy
  dadd r874.zw, r510.zw, r304_neg(yw).zw
  dadd r877.xy, r304.xy, r104_neg(yw).xy
  dadd r877.zw, r304.zw, r104_neg(yw).zw
  dadd r878.xy, r305.xy, r105_neg(yw).xy
  dadd r878.zw, r305.zw, r105_neg(yw).zw
  dmul r768.xy, r867.xy, r867.xy
  dmul r768.zw, r867.zw, r867.zw
  dmul r799.xy, r891.xy, r891.xy
  dmul r799.zw, r891.zw, r891.zw
  dadd r703.xy, r768.xy, r799.xy
  dadd r703.zw, r768.zw, r799.zw
  dmul r801.xy, r871.xy, r871.xy
  dmul r801.zw, r871.zw, r871.zw
  dadd r719.xy, r703.xy, r801.xy
  dadd r719.zw, r703.zw, r801.zw
  dadd r720.xy, r719.xy, r511.xy
  dadd r720.zw, r719.zw, r511.zw
  mov r10.xy, r720.xy
  call 15
  mov r744.xy, r11.xy
  mov r10.xy, r720.zw
  call 15
  mov r744.__zw, r11.xyxy
  dmul r806.xy, r880.xy, r880.xy
  dmul r806.zw, r880.zw, r880.zw
  dmul r809.xy, r893.xy, r893.xy
  dmul r809.zw, r893.zw, r893.zw
  dadd r721.xy, r806.xy, r809.xy
  dadd r721.zw, r806.zw, r809.zw
  dmul r810.xy, r874.xy, r874.xy
  dmul r810.zw, r874.zw, r874.zw
  dadd r723.xy, r721.xy, r810.xy
  dadd r723.zw, r721.zw, r810.zw
  dadd r726.xy, r723.xy, r511.xy
  dadd r726.zw, r723.zw, r511.zw
  mov r10.xy, r726.xy
  call 15
  mov r745.xy, r11.xy
  mov r10.xy, r726.zw
  call 15
  mov r745.__zw, r11.xyxy
  dmul r814.xy, r884.xy, r884.xy
  dmul r814.zw, r884.zw, r884.zw
  dmul r816.xy, r895.xy, r895.xy
  dmul r816.zw, r895.zw, r895.zw
  dadd r728.xy, r814.xy, r816.xy
  dadd r728.zw, r814.zw, r816.zw
  dmul r820.xy, r877.xy, r877.xy
  dmul r820.zw, r877.zw, r877.zw
  dadd r729.xy, r728.xy, r820.xy
  dadd r729.zw, r728.zw, r820.zw
  dadd r730.xy, r729.xy, r511.xy
  dadd r730.zw, r729.zw, r511.zw
  mov r10.xy, r730.xy
  call 15
  mov r746.xy, r11.xy
  mov r10.xy, r730.zw
  call 15
  mov r746.__zw, r11.xyxy
  dmul r824.xy, r746.xy, r306.xy
  dmul r824.zw, r746.zw, r306.zw
  dmul r826.xy, r824.xy, r746.xy
  dmul r826.zw, r824.zw, r746.zw
  dmul r827.xy, r826.xy, r746.xy
  dmul r827.zw, r826.zw, r746.zw
  dmul r828.xy, r744.xy, r306.xy
  dmul r828.zw, r744.zw, r306.zw
  dmul r829.xy, r828.xy, r744.xy
  dmul r829.zw, r828.zw, r744.zw
  dmul r830.xy, r829.xy, r744.xy
  dmul r830.zw, r829.zw, r744.zw
  dmul r831.xy, r830.xy, r512.xy
  dmul r831.zw, r830.zw, r512.zw
  dmul r832.xy, r746.xy, r507.xy
  dmul r832.zw, r746.zw, r507.zw
  dmul r833.xy, r745.xy, r502.xy
  dmul r833.zw, r745.zw, r502.zw
  dadd r734.xy, r832.xy, r833.xy
  dadd r734.zw, r832.zw, r833.zw
  dadd r890.xy, r719.xy, r729_neg(yw).xy
  dadd r890.zw, r719.zw, r729_neg(yw).zw
  dmul r834.xy, r890.xy, r501.xy
  dmul r834.zw, r890.zw, r501.zw
  dmul r835.xy, r834.xy, r745.xy
  dmul r835.zw, r834.zw, r745.zw
  dmul r836.xy, r835.xy, r745.xy
  dmul r836.zw, r835.zw, r745.zw
  dmul r837.xy, r836.xy, r745.xy
  dmul r837.zw, r836.zw, r745.zw
  dadd r735.xy, r734.xy, r837.xy
  dadd r735.zw, r734.zw, r837.zw
  dmul r838.xy, r831.xy, r735.xy
  dmul r838.zw, r831.zw, r735.zw
  dmul r839.xy, r838.xy, r867.xy
  dmul r839.zw, r838.zw, r867.zw
  dmul r840.xy, r827.xy, r512.xy
  dmul r840.zw, r827.zw, r512.zw
  dmul r841.xy, r744.xy, r507.xy
  dmul r841.zw, r744.zw, r507.zw
  dadd r736.xy, r841.xy, r833.xy
  dadd r736.zw, r841.zw, r833.zw
  dadd r737.xy, r729.xy, r719_neg(yw).xy
  dadd r737.zw, r729.zw, r719_neg(yw).zw
  dmul r842.xy, r737.xy, r501.xy
  dmul r842.zw, r737.zw, r501.zw
  dmul r843.xy, r842.xy, r745.xy
  dmul r843.zw, r842.zw, r745.zw
  dmul r844.xy, r843.xy, r745.xy
  dmul r844.zw, r843.zw, r745.zw
  dmul r845.xy, r844.xy, r745.xy
  dmul r845.zw, r844.zw, r745.zw
  dadd r738.xy, r736.xy, r845.xy
  dadd r738.zw, r736.zw, r845.zw
  dmul r846.xy, r840.xy, r738.xy
  dmul r846.zw, r840.zw, r738.zw
  dmul r847.xy, r846.xy, r884.xy
  dmul r847.zw, r846.zw, r884.zw
  dadd r739.xy, r839.xy, r847.xy
  dadd r739.zw, r839.zw, r847.zw
  dmul r848.xy, r306.xy, r506.xy
  dmul r848.zw, r306.zw, r506.zw
  dmul r849.xy, r848.xy, r512.xy
  dmul r849.zw, r848.zw, r512.zw
  dmul r850.xy, r849.xy, r745.xy
  dmul r850.zw, r849.zw, r745.zw
  dmul r851.xy, r850.xy, r745.xy
  dmul r851.zw, r850.zw, r745.zw
  dmul r852.xy, r851.xy, r745.xy
  dmul r852.zw, r851.zw, r745.zw
  dadd r892.xy, r746.xy, r744_neg(yw).xy
  dadd r892.zw, r746.zw, r744_neg(yw).zw
  dmul r853.xy, r852.xy, r892.xy
  dmul r853.zw, r852.zw, r892.zw
  dmul r854.xy, r853.xy, r880.xy
  dmul r854.zw, r853.zw, r880.zw
  dadd r894.xy, r739.xy, r854_neg(yw).xy
  dadd r894.zw, r739.zw, r854_neg(yw).zw
  dmul r855.xy, r101.xy, r301.xy
  dmul r855.zw, r101.zw, r301.zw
  dmul r856.xy, r103.xy, r303.xy
  dmul r856.zw, r103.zw, r303.zw
  dadd r740.xy, r855.xy, r856.xy
  dadd r740.zw, r855.zw, r856.zw
  dmul r857.xy, r105.xy, r305.xy
  dmul r857.zw, r105.zw, r305.zw
  dadd r741.xy, r740.xy, r857.xy
  dadd r741.zw, r740.zw, r857.zw
  dmul r858.xy, r741.xy, r507.xy
  dmul r858.zw, r741.zw, r507.zw
  dmul r859.xy, r858.xy, r867.xy
  dmul r859.zw, r858.zw, r867.zw
  dmul r860.xy, r867.xy, r301.xy
  dmul r860.zw, r867.zw, r301.zw
  dmul r861.xy, r891.xy, r303.xy
  dmul r861.zw, r891.zw, r303.zw
  dadd r742.xy, r860.xy, r861.xy
  dadd r742.zw, r860.zw, r861.zw
  dmul r862.xy, r871.xy, r305.xy
  dmul r862.zw, r871.zw, r305.zw
  dadd r743.xy, r742.xy, r862.xy
  dadd r743.zw, r742.zw, r862.zw
  dmul r863.xy, r743.xy, r505.xy
  dmul r863.zw, r743.zw, r505.zw
  dmul r864.xy, r863.xy, r101.xy
  dmul r864.zw, r863.zw, r101.zw
  dadd r896.xy, r859.xy, r864_neg(yw).xy
  dadd r896.zw, r859.zw, r864_neg(yw).zw
  dmul r865.xy, r867.xy, r101.xy
  dmul r865.zw, r867.zw, r101.zw
  dmul r769.xy, r891.xy, r103.xy
  dmul r769.zw, r891.zw, r103.zw
  dadd r704.xy, r865.xy, r769.xy
  dadd r704.zw, r865.zw, r769.zw
  dmul r770.xy, r871.xy, r105.xy
  dmul r770.zw, r871.zw, r105.zw
  dadd r705.xy, r704.xy, r770.xy
  dadd r705.zw, r704.zw, r770.zw
  dmul r771.xy, r705.xy, r507.xy
  dmul r771.zw, r705.zw, r507.zw
  dmul r772.xy, r771.xy, r301.xy
  dmul r772.zw, r771.zw, r301.zw
  dadd r869.xy, r896.xy, r772_neg(yw).xy
  dadd r869.zw, r896.zw, r772_neg(yw).zw
  dmul r773.xy, r830.xy, r869.xy
  dmul r773.zw, r830.zw, r869.zw
  dadd r870.xy, r894.xy, r773_neg(yw).xy
  dadd r870.zw, r894.zw, r773_neg(yw).zw
  dmul r774.xy, r101.xy, r101.xy
  dmul r774.zw, r101.zw, r101.zw
  dmul r775.xy, r103.xy, r103.xy
  dmul r775.zw, r103.zw, r103.zw
  dadd r706.xy, r774.xy, r775.xy
  dadd r706.zw, r774.zw, r775.zw
  dmul r776.xy, r105.xy, r105.xy
  dmul r776.zw, r105.zw, r105.zw
  dadd r707.xy, r706.xy, r776.xy
  dadd r707.zw, r706.zw, r776.zw
  dmul r777.xy, r889.xy, r889.xy
  dmul r777.zw, r889.zw, r889.zw
  dmul r778.xy, r868.xy, r868.xy
  dmul r778.zw, r868.zw, r868.zw
  dadd r708.xy, r777.xy, r778.xy
  dadd r708.zw, r777.zw, r778.zw
  dmul r779.xy, r878.xy, r878.xy
  dmul r779.zw, r878.zw, r878.zw
  dadd r709.xy, r708.xy, r779.xy
  dadd r709.zw, r708.zw, r779.zw
  dmul r780.xy, r709.xy, r504.xy
  dmul r780.zw, r709.zw, r504.zw
  dadd r872.xy, r707.xy, r780_neg(yw).xy
  dadd r872.zw, r707.zw, r780_neg(yw).zw
  dmul r781.xy, r884.xy, r301.xy
  dmul r781.zw, r884.zw, r301.zw
  dmul r782.xy, r895.xy, r303.xy
  dmul r782.zw, r895.zw, r303.zw
  dadd r710.xy, r781.xy, r782.xy
  dadd r710.zw, r781.zw, r782.zw
  dmul r783.xy, r877.xy, r305.xy
  dmul r783.zw, r877.zw, r305.zw
  dadd r711.xy, r710.xy, r783.xy
  dadd r711.zw, r710.zw, r783.zw
  dmul r784.xy, r711.xy, r503.xy
  dmul r784.zw, r711.zw, r503.zw
  dmul r785.xy, r784.xy, r711.xy
  dmul r785.zw, r784.zw, r711.zw
  ddiv r767.xy, r785.xy, r730.xy
  ddiv r767.zw, r785.zw, r730.zw
  dadd r712.xy, r872.xy, r767.xy
  dadd r712.zw, r872.zw, r767.zw
  dmul r786.xy, r827.xy, r712.xy
  dmul r786.zw, r827.zw, r712.zw
  dmul r787.xy, r786.xy, r884.xy
  dmul r787.zw, r786.zw, r884.zw
  dadd r713.xy, r870.xy, r787.xy
  dadd r713.zw, r870.zw, r787.zw
  dmul r788.xy, r101.xy, r507.xy
  dmul r788.zw, r101.zw, r507.zw
  dmul r789.xy, r301.xy, r505.xy
  dmul r789.zw, r301.zw, r505.zw
  dadd r873.xy, r788.xy, r789_neg(yw).xy
  dadd r873.zw, r788.zw, r789_neg(yw).zw
  dmul r790.xy, r884.xy, r873.xy
  dmul r790.zw, r884.zw, r873.zw
  dmul r791.xy, r103.xy, r507.xy
  dmul r791.zw, r103.zw, r507.zw
  dmul r792.xy, r303.xy, r505.xy
  dmul r792.zw, r303.zw, r505.zw
  dadd r875.xy, r791.xy, r792_neg(yw).xy
  dadd r875.zw, r791.zw, r792_neg(yw).zw
  dmul r793.xy, r895.xy, r875.xy
  dmul r793.zw, r895.zw, r875.zw
  dadd r714.xy, r790.xy, r793.xy
  dadd r714.zw, r790.zw, r793.zw
  dmul r794.xy, r105.xy, r507.xy
  dmul r794.zw, r105.zw, r507.zw
  dmul r795.xy, r305.xy, r505.xy
  dmul r795.zw, r305.zw, r505.zw
  dadd r876.xy, r794.xy, r795_neg(yw).xy
  dadd r876.zw, r794.zw, r795_neg(yw).zw
  dmul r796.xy, r877.xy, r876.xy
  dmul r796.zw, r877.zw, r876.zw
  dadd r715.xy, r714.xy, r796.xy
  dadd r715.zw, r714.zw, r796.zw
  dmul r797.xy, r827.xy, r715.xy
  dmul r797.zw, r827.zw, r715.zw
  dmul r798.xy, r797.xy, r889.xy
  dmul r798.zw, r797.zw, r889.zw
  dadd r716.xy, r713.xy, r798.xy
  dadd r716.zw, r713.zw, r798.zw
  dadd r717.xy, r716.xy, r500.xy
  dadd r717.zw, r716.zw, r500.zw
  dmul r800.xy, r838.xy, r891.xy
  dmul r800.zw, r838.zw, r891.zw
  dmul r802.xy, r846.xy, r895.xy
  dmul r802.zw, r846.zw, r895.zw
  dadd r718.xy, r800.xy, r802.xy
  dadd r718.zw, r800.zw, r802.zw
  dmul r803.xy, r853.xy, r893.xy
  dmul r803.zw, r853.zw, r893.zw
  dadd r879.xy, r718.xy, r803_neg(yw).xy
  dadd r879.zw, r718.zw, r803_neg(yw).zw
  dmul r804.xy, r858.xy, r891.xy
  dmul r804.zw, r858.zw, r891.zw
  dmul r805.xy, r863.xy, r103.xy
  dmul r805.zw, r863.zw, r103.zw
  dadd r881.xy, r804.xy, r805_neg(yw).xy
  dadd r881.zw, r804.zw, r805_neg(yw).zw
  dmul r807.xy, r771.xy, r303.xy
  dmul r807.zw, r771.zw, r303.zw
  dadd r882.xy, r881.xy, r807_neg(yw).xy
  dadd r882.zw, r881.zw, r807_neg(yw).zw
  dmul r808.xy, r830.xy, r882.xy
  dmul r808.zw, r830.zw, r882.zw
  dadd r883.xy, r879.xy, r808_neg(yw).xy
  dadd r883.zw, r879.zw, r808_neg(yw).zw
  dmul r811.xy, r786.xy, r895.xy
  dmul r811.zw, r786.zw, r895.zw
  dadd r722.xy, r883.xy, r811.xy
  dadd r722.zw, r883.zw, r811.zw
  dmul r812.xy, r797.xy, r868.xy
  dmul r812.zw, r797.zw, r868.zw
  dadd r724.xy, r722.xy, r812.xy
  dadd r724.zw, r722.zw, r812.zw
  dadd r725.xy, r724.xy, r500.xy
  dadd r725.zw, r724.zw, r500.zw
  dmul r813.xy, r838.xy, r871.xy
  dmul r813.zw, r838.zw, r871.zw
  dmul r815.xy, r846.xy, r877.xy
  dmul r815.zw, r846.zw, r877.zw
  dadd r727.xy, r813.xy, r815.xy
  dadd r727.zw, r813.zw, r815.zw
  dmul r817.xy, r853.xy, r874.xy
  dmul r817.zw, r853.zw, r874.zw
  dadd r885.xy, r727.xy, r817_neg(yw).xy
  dadd r885.zw, r727.zw, r817_neg(yw).zw
  dmul r818.xy, r858.xy, r871.xy
  dmul r818.zw, r858.zw, r871.zw
  dmul r819.xy, r863.xy, r105.xy
  dmul r819.zw, r863.zw, r105.zw
  dadd r886.xy, r818.xy, r819_neg(yw).xy
  dadd r886.zw, r818.zw, r819_neg(yw).zw
  dmul r821.xy, r771.xy, r305.xy
  dmul r821.zw, r771.zw, r305.zw
  dadd r887.xy, r886.xy, r821_neg(yw).xy
  dadd r887.zw, r886.zw, r821_neg(yw).zw
  dmul r822.xy, r830.xy, r887.xy
  dmul r822.zw, r830.zw, r887.zw
  dadd r888.xy, r885.xy, r822_neg(yw).xy
  dadd r888.zw, r885.zw, r822_neg(yw).zw
  dmul r823.xy, r786.xy, r877.xy
  dmul r823.zw, r786.zw, r877.zw
  dadd r731.xy, r888.xy, r823.xy
  dadd r731.zw, r888.zw, r823.zw
  dmul r825.xy, r797.xy, r878.xy
  dmul r825.zw, r797.zw, r878.zw
  dadd r732.xy, r731.xy, r825.xy
  dadd r732.zw, r731.zw, r825.zw
  dadd r733.xy, r732.xy, r500.xy
  dadd r733.zw, r732.zw, r500.zw
  dadd r400.xy, r400.xy, r717.xy
  dadd r400.zw, r400.zw, r717.zw
  dadd r401.xy, r401.xy, r725.xy
  dadd r401.zw, r401.zw, r725.zw
  dadd r402.xy, r402.xy, r733.xy
  dadd r402.zw, r402.zw, r733.zw

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

