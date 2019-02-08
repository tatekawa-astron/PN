CONST x_0_0_, x_0_1_, x_0_2_, eps2, m_0_;
VARI x_i_0_, x_i_1_, x_i_2_, m_i_, v_i_0_, v_i_1_, v_i_2_;
VARJ x_j_0_, x_j_1_, x_j_2_, m_j_, v_j_0_, v_j_1_, v_j_2_;
VARF a1c_0_, a1c_1_, a1c_2_;
a1c_2__ = 0.0;
a1c_1__ = 0.0;
a1c_0__ = 0.0;

dxb_0_ = 0.0;  // the value 0.0 is meaningless.
dxc_0_ = 0.0;  // the value 0.0 is meaningless.
dxbc_0_ = 0.0;  // the value 0.0 is meaningless.
dxb_1_ = 0.0;  // the value 0.0 is meaningless.
dxc_1_ = 0.0;  // the value 0.0 is meaningless.
dxbc_1_ = 0.0;  // the value 0.0 is meaningless.
dxb_2_ = 0.0;  // the value 0.0 is meaningless.
dxc_2_ = 0.0;  // the value 0.0 is meaningless.
dxbc_2_ = 0.0;  // the value 0.0 is meaningless.
r1b2 = 0.0;  // the value 0.0 is meaningless.
r1b2e = 0.0;  // the value 0.0 is meaningless.
r1be = 0.0;  // the value 0.0 is meaningless.
r1c2 = 0.0;  // the value 0.0 is meaningless.
r1c2e = 0.0;  // the value 0.0 is meaningless.
r1ce = 0.0;  // the value 0.0 is meaningless.
rbc2e = 0.0;  // the value 0.0 is meaningless.
rbce = 0.0;  // the value 0.0 is meaningless.
mr1b3e = 0.0;  // the value 0.0 is meaningless.
dxb_0_ = x_i_0_ - x_0_0_;
dxc_0_ = x_j_0_ - x_0_0_;
dxbc_0_ = x_j_0_ - x_i_0_;
dxb_1_ = x_i_1_ - x_0_1_;
dxc_1_ = x_j_1_ - x_0_1_;
dxbc_1_ = x_j_1_ - x_i_1_;
dxb_2_ = x_i_2_ - x_0_2_;
dxc_2_ = x_j_2_ - x_0_2_;
dxbc_2_ = x_j_2_ - x_i_2_;
r1b2 = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_;
r1b2e = r1b2 + eps2;
r1be = rsqrt(r1b2e);
r1c2 = dxc_0_ * dxc_0_ + dxc_1_ * dxc_1_ + dxc_2_ * dxc_2_;
r1c2e = r1c2 + eps2;
r1ce = rsqrt(r1c2e);
rbc2e = dxbc_0_ * dxbc_0_ + dxbc_1_ * dxbc_1_ + dxbc_2_ * dxbc_2_ + eps2;
rbce = rsqrt(rbc2e);
mr1b3e = m_i_ * r1be * r1be * r1be;
a1c_0__ = a1c_0__ + (mr1b3e * m_j_ * dxb_0_ * (4.0 * r1ce + 1.25 * rbce - 0.25 * r1c2 / rbc2e * rbce + 0.25 * r1b2 / rbc2e * rbce) - 3.5 * (rbce * rbce * rbce) * r1be * m_i_ * m_j_ * dxbc_0_ - mr1b3e * m_j_ / m_0_ * (4.0 * (v_i_0_ * v_j_0_ + v_i_1_ * v_j_1_ + v_i_2_ * v_j_2_) * dxb_0_ - 3.0 * (v_i_0_ * dxb_0_ + v_i_1_ * dxb_1_ + v_i_2_ * dxb_2_) * v_j_0_ - 4.0 * (v_j_0_ * dxb_0_ + v_j_1_ * dxb_1_ + v_j_2_ * dxb_2_) * v_i_0_));
a1c_1__ = a1c_1__ + (mr1b3e * m_j_ * dxb_1_ * (4.0 * r1ce + 1.25 * rbce - 0.25 * r1c2 / rbc2e * rbce + 0.25 * r1b2 / rbc2e * rbce) - 3.5 * (rbce * rbce * rbce) * r1be * m_i_ * m_j_ * dxbc_1_ - mr1b3e * m_j_ / m_0_ * (4.0 * (v_i_0_ * v_j_0_ + v_i_1_ * v_j_1_ + v_i_2_ * v_j_2_) * dxb_1_ - 3.0 * (v_i_0_ * dxb_0_ + v_i_1_ * dxb_1_ + v_i_2_ * dxb_2_) * v_j_1_ - 4.0 * (v_j_0_ * dxb_0_ + v_j_1_ * dxb_1_ + v_j_2_ * dxb_2_) * v_i_1_));
a1c_2__ = a1c_2__ + (mr1b3e * m_j_ * dxb_2_ * (4.0 * r1ce + 1.25 * rbce - 0.25 * r1c2 / rbc2e * rbce + 0.25 * r1b2 / rbc2e * rbce) - 3.5 * (rbce * rbce * rbce) * r1be * m_i_ * m_j_ * dxbc_2_ - mr1b3e * m_j_ / m_0_ * (4.0 * (v_i_0_ * v_j_0_ + v_i_1_ * v_j_1_ + v_i_2_ * v_j_2_) * dxb_2_ - 3.0 * (v_i_0_ * dxb_0_ + v_i_1_ * dxb_1_ + v_i_2_ * dxb_2_) * v_j_2_ - 4.0 * (v_j_0_ * dxb_0_ + v_j_1_ * dxb_1_ + v_j_2_ * dxb_2_) * v_i_2_));
a1c_0_ += a1c_0__;
a1c_1_ += a1c_1__;
a1c_2_ += a1c_2__;
