CONST x_0_0_, x_0_1_, x_0_2_, eps2, m_0_;
VARI x_i_0_, v_i_0_, x_i_1_, v_i_1_, x_i_2_, v_i_2_;
VARJ x_j_0_, v_j_0_, x_j_1_, v_j_1_, x_j_2_, v_j_2_, m_j_;
VARF ac_i_0_, ac_i_1_, ac_i_2_;
ac_i_2__ = 0.0;
ac_i_1__ = 0.0;
ac_i_0__ = 0.0;

dxb_0_ = 0.0;  // the value 0.0 is meaningless.
dxc_0_ = 0.0;  // the value 0.0 is meaningless.
dxbc_0_ = 0.0;  // the value 0.0 is meaningless.
dvbc_0_ = 0.0;  // the value 0.0 is meaningless.
dxb_1_ = 0.0;  // the value 0.0 is meaningless.
dxc_1_ = 0.0;  // the value 0.0 is meaningless.
dxbc_1_ = 0.0;  // the value 0.0 is meaningless.
dvbc_1_ = 0.0;  // the value 0.0 is meaningless.
dxb_2_ = 0.0;  // the value 0.0 is meaningless.
dxc_2_ = 0.0;  // the value 0.0 is meaningless.
dxbc_2_ = 0.0;  // the value 0.0 is meaningless.
dvbc_2_ = 0.0;  // the value 0.0 is meaningless.
r1b2 = 0.0;  // the value 0.0 is meaningless.
r1b2e = 0.0;  // the value 0.0 is meaningless.
r1be = 0.0;  // the value 0.0 is meaningless.
r1c2 = 0.0;  // the value 0.0 is meaningless.
r1c2e = 0.0;  // the value 0.0 is meaningless.
r1ce = 0.0;  // the value 0.0 is meaningless.
rbc2 = 0.0;  // the value 0.0 is meaningless.
rbc2e = 0.0;  // the value 0.0 is meaningless.
rbce = 0.0;  // the value 0.0 is meaningless.
mr3inv = 0.0;  // the value 0.0 is meaningless.
dxb_0_ = x_0_0_ - x_i_0_;
dxc_0_ = x_0_0_ - x_j_0_;
dxbc_0_ = x_j_0_ - x_i_0_;
dvbc_0_ = v_j_0_ - v_i_0_;
dxb_1_ = x_0_1_ - x_i_1_;
dxc_1_ = x_0_1_ - x_j_1_;
dxbc_1_ = x_j_1_ - x_i_1_;
dvbc_1_ = v_j_1_ - v_i_1_;
dxb_2_ = x_0_2_ - x_i_2_;
dxc_2_ = x_0_2_ - x_j_2_;
dxbc_2_ = x_j_2_ - x_i_2_;
dvbc_2_ = v_j_2_ - v_i_2_;
r1b2 = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_;
r1b2e = r1b2 + eps2;
r1be = rsqrt(r1b2e);
r1c2 = dxc_0_ * dxc_0_ + dxc_1_ * dxc_1_ + dxc_2_ * dxc_2_;
r1c2e = r1c2 + eps2;
r1ce = rsqrt(r1c2e);
rbc2 = dxbc_0_ * dxbc_0_ + dxbc_1_ * dxbc_1_ + dxbc_2_ * dxbc_2_;
rbc2e = rbc2 + eps2;
rbce = rsqrt(rbc2e);
mr3inv = m_j_ * rbce * rbce * rbce;
ac_i_0__ = ac_i_0__ + (m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_0_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_0_ - 3.5 * m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_0_ - m_j_ * r1be * r1be * r1be * (4.0 * (v_i_0_ * v_j_0_ + v_i_1_ * v_j_1_ + v_i_2_ * v_j_2_) * dxb_0_ - 3.0 * (v_j_0_ * dxb_0_ + v_j_1_ * dxb_1_ + v_j_2_ * dxb_2_) * v_i_0_ - 4.0 * (v_i_0_ * dxb_0_ + v_i_1_ * dxb_1_ + v_i_2_ * dxb_2_) * v_j_0_) + mr3inv * (v_i_0_ * v_i_0_ + v_i_1_ * v_i_1_ + v_i_2_ * v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (v_j_0_ * dxbc_0_ + v_j_1_ * dxbc_1_ + v_j_2_ * dxbc_2_) * (v_j_0_ * dxbc_0_ + v_j_1_ * dxbc_1_ + v_j_2_ * dxbc_2_) / rbc2e) * dxbc_0_ + mr3inv * (dxbc_0_ * (4.0 * v_i_0_ - 3.0 * v_j_0_) + dxbc_1_ * (4.0 * v_i_1_ - 3.0 * v_j_1_) + dxbc_2_ * (4.0 * v_i_2_ - 3.0 * v_j_2_)) * dvbc_0_);
ac_i_1__ = ac_i_1__ + (m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_1_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_1_ - 3.5 * m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_1_ - m_j_ * r1be * r1be * r1be * (4.0 * (v_i_0_ * v_j_0_ + v_i_1_ * v_j_1_ + v_i_2_ * v_j_2_) * dxb_1_ - 3.0 * (v_j_0_ * dxb_0_ + v_j_1_ * dxb_1_ + v_j_2_ * dxb_2_) * v_i_1_ - 4.0 * (v_i_0_ * dxb_0_ + v_i_1_ * dxb_1_ + v_i_2_ * dxb_2_) * v_j_1_) + mr3inv * (v_i_0_ * v_i_0_ + v_i_1_ * v_i_1_ + v_i_2_ * v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (v_j_0_ * dxbc_0_ + v_j_1_ * dxbc_1_ + v_j_2_ * dxbc_2_) * (v_j_0_ * dxbc_0_ + v_j_1_ * dxbc_1_ + v_j_2_ * dxbc_2_) / rbc2e) * dxbc_1_ + mr3inv * (dxbc_0_ * (4.0 * v_i_0_ - 3.0 * v_j_0_) + dxbc_1_ * (4.0 * v_i_1_ - 3.0 * v_j_1_) + dxbc_2_ * (4.0 * v_i_2_ - 3.0 * v_j_2_)) * dvbc_1_);
ac_i_2__ = ac_i_2__ + (m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_2_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_2_ - 3.5 * m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_2_ - m_j_ * r1be * r1be * r1be * (4.0 * (v_i_0_ * v_j_0_ + v_i_1_ * v_j_1_ + v_i_2_ * v_j_2_) * dxb_2_ - 3.0 * (v_j_0_ * dxb_0_ + v_j_1_ * dxb_1_ + v_j_2_ * dxb_2_) * v_i_2_ - 4.0 * (v_i_0_ * dxb_0_ + v_i_1_ * dxb_1_ + v_i_2_ * dxb_2_) * v_j_2_) + mr3inv * (v_i_0_ * v_i_0_ + v_i_1_ * v_i_1_ + v_i_2_ * v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (v_j_0_ * dxbc_0_ + v_j_1_ * dxbc_1_ + v_j_2_ * dxbc_2_) * (v_j_0_ * dxbc_0_ + v_j_1_ * dxbc_1_ + v_j_2_ * dxbc_2_) / rbc2e) * dxbc_2_ + mr3inv * (dxbc_0_ * (4.0 * v_i_0_ - 3.0 * v_j_0_) + dxbc_1_ * (4.0 * v_i_1_ - 3.0 * v_j_1_) + dxbc_2_ * (4.0 * v_i_2_ - 3.0 * v_j_2_)) * dvbc_2_);
ac_i_0_ += ac_i_0__;
ac_i_1_ += ac_i_1__;
ac_i_2_ += ac_i_2__;
