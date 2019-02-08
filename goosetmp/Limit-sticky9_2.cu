               __global__ void
               f02_calculator(int ioff_, int ni_, int nj_, double x_0_0_, double x_0_1_, double x_0_2_, double eps2, double m_0_, f02_jp_t *f02_jp_, f02_ip_t *f02_ip_, f02_result_t *f02_result_)
               {
                   extern __shared__ char smembuf_[];
                   int kbdim_ = blockDim.x;
                   
                   f02_jp_t * f02_jp_smem_ = (f02_jp_t *)smembuf_;
                   
                   f02_result_t * f02_result_smem_ = (f02_result_t *)smembuf_;
                   double dxb_0_, dxc_0_, dxbc_0_, dvbc_0_, dxb_1_, dxc_1_, dxbc_1_, dvbc_1_, dxb_2_, dxc_2_, dxbc_2_, dvbc_2_, r1b2, r1b2e, r1be, r1c2, r1c2e, r1ce, rbc2, rbc2e, rbce, mr3inv;
                   double ac_i_0_wcache_ = 0.0f;
double ac_i_1_wcache_ = 0.0f;
double ac_i_2_wcache_ = 0.0f;

                   int njdiv_ = gridDim.x;
                   int jbid_ = blockIdx.x;
                   int ibid_ = blockIdx.y;
                   int tid_ = threadIdx.x;
                   int npipe_ = (ni_ - 1) / 1 + 1;
                   int nvalidthread_ = npipe_ - kbdim_ * ibid_;
                   if (nvalidthread_ > kbdim_) {
                       nvalidthread_ = kbdim_;
                   }
                   int njsub_ = (nj_ - 1) / njdiv_ + 1;
                   int joff0_ = njsub_ * jbid_;
                   int joff1_ = joff0_ + njsub_;
                   if (joff1_ > nj_) {
                       joff1_ = nj_;
                   }
                   int jstride_ = 1;
                   if (nvalidthread_ <= kbdim_ / 2) jstride_ = 2;
if (nvalidthread_ <= kbdim_ / 4) jstride_ = 4;

                   int njhsub_ = kbdim_ / jstride_;
                   int jstart_ = tid_ / njhsub_;
                   int isrc_ = kbdim_ * ibid_ + tid_ % njhsub_;
                   int idst_ = njdiv_ * isrc_ + jbid_;
                   int icnt_ = ioff_ + isrc_;
                   for (int joff_ = joff0_; joff_ < joff1_; joff_ += kbdim_) {
                       int jsrc_ = tid_+joff_;
                       __syncthreads();
                       #if 1
                     {
                         float4 *srcbuf_ = (float4 *)(f02_jp_ + joff_);
                         float4 *dstbuf_ = (float4 *)smembuf_;
                         for (int icpy = 0; icpy < sizeof(f02_jp_t)/sizeof(float4); icpy++) {
                             dstbuf_[tid_] = srcbuf_[tid_];
                             dstbuf_ += kbdim_;
                             srcbuf_ += kbdim_;
                         }
                     }
#else
                     f02_jp_smem_[tid_] = f02_jp_[jsrc_];
#endif

                       __syncthreads();
                       int jsup_ = kbdim_;
                       if (joff_ + jsup_ > joff1_) {
                           jsup_ = joff1_ - joff_;
                       }
                       if (jsup_ < kbdim_) {
                           for (int j_ = jstart_; j_ < jsup_; j_+= jstride_) {
                               dxb_0_ = x_0_0_ - f02_ip_[isrc_].x_i_0_;
dxc_0_ = x_0_0_ - f02_jp_smem_[j_].x_j_0_;
dxbc_0_ = f02_jp_smem_[j_].x_j_0_ - f02_ip_[isrc_].x_i_0_;
dvbc_0_ = f02_jp_smem_[j_].v_j_0_ - f02_ip_[isrc_].v_i_0_;
dxb_1_ = x_0_1_ - f02_ip_[isrc_].x_i_1_;
dxc_1_ = x_0_1_ - f02_jp_smem_[j_].x_j_1_;
dxbc_1_ = f02_jp_smem_[j_].x_j_1_ - f02_ip_[isrc_].x_i_1_;
dvbc_1_ = f02_jp_smem_[j_].v_j_1_ - f02_ip_[isrc_].v_i_1_;
dxb_2_ = x_0_2_ - f02_ip_[isrc_].x_i_2_;
dxc_2_ = x_0_2_ - f02_jp_smem_[j_].x_j_2_;
dxbc_2_ = f02_jp_smem_[j_].x_j_2_ - f02_ip_[isrc_].x_i_2_;
dvbc_2_ = f02_jp_smem_[j_].v_j_2_ - f02_ip_[isrc_].v_i_2_;
r1b2 = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_;
r1b2e = r1b2 + eps2;
r1be = rsqrt(r1b2e);
r1c2 = dxc_0_ * dxc_0_ + dxc_1_ * dxc_1_ + dxc_2_ * dxc_2_;
r1c2e = r1c2 + eps2;
r1ce = rsqrt(r1c2e);
rbc2 = dxbc_0_ * dxbc_0_ + dxbc_1_ * dxbc_1_ + dxbc_2_ * dxbc_2_;
rbc2e = rbc2 + eps2;
rbce = rsqrt(rbc2e);
mr3inv = f02_jp_smem_[j_].m_j_ * rbce * rbce * rbce;
ac_i_0_wcache_ += f02_jp_smem_[j_].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_0_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_0_ - 3.5 * f02_jp_smem_[j_].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_0_ - f02_jp_smem_[j_].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_].v_j_2_) * dxb_0_ - 3.0 * (f02_jp_smem_[j_].v_j_0_ * dxb_0_ + f02_jp_smem_[j_].v_j_1_ * dxb_1_ + f02_jp_smem_[j_].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_0_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_].v_j_0_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_0_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_].v_j_2_)) * dvbc_0_;
ac_i_1_wcache_ += f02_jp_smem_[j_].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_1_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_1_ - 3.5 * f02_jp_smem_[j_].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_1_ - f02_jp_smem_[j_].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_].v_j_2_) * dxb_1_ - 3.0 * (f02_jp_smem_[j_].v_j_0_ * dxb_0_ + f02_jp_smem_[j_].v_j_1_ * dxb_1_ + f02_jp_smem_[j_].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_1_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_].v_j_1_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_1_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_].v_j_2_)) * dvbc_1_;
ac_i_2_wcache_ += f02_jp_smem_[j_].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_2_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_2_ - 3.5 * f02_jp_smem_[j_].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_2_ - f02_jp_smem_[j_].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_].v_j_2_) * dxb_2_ - 3.0 * (f02_jp_smem_[j_].v_j_0_ * dxb_0_ + f02_jp_smem_[j_].v_j_1_ * dxb_1_ + f02_jp_smem_[j_].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_2_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_].v_j_2_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_2_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_].v_j_2_)) * dvbc_2_;

                           }
                       }
                       else {    
                           for (int j_ = jstart_; j_ < kbdim_; j_+= jstride_ * 8) {
                                       // loop 0
dxb_0_ = x_0_0_ - f02_ip_[isrc_].x_i_0_;
dxc_0_ = x_0_0_ - f02_jp_smem_[j_ + jstride_ * 0].x_j_0_;
dxbc_0_ = f02_jp_smem_[j_ + jstride_ * 0].x_j_0_ - f02_ip_[isrc_].x_i_0_;
dvbc_0_ = f02_jp_smem_[j_ + jstride_ * 0].v_j_0_ - f02_ip_[isrc_].v_i_0_;
dxb_1_ = x_0_1_ - f02_ip_[isrc_].x_i_1_;
dxc_1_ = x_0_1_ - f02_jp_smem_[j_ + jstride_ * 0].x_j_1_;
dxbc_1_ = f02_jp_smem_[j_ + jstride_ * 0].x_j_1_ - f02_ip_[isrc_].x_i_1_;
dvbc_1_ = f02_jp_smem_[j_ + jstride_ * 0].v_j_1_ - f02_ip_[isrc_].v_i_1_;
dxb_2_ = x_0_2_ - f02_ip_[isrc_].x_i_2_;
dxc_2_ = x_0_2_ - f02_jp_smem_[j_ + jstride_ * 0].x_j_2_;
dxbc_2_ = f02_jp_smem_[j_ + jstride_ * 0].x_j_2_ - f02_ip_[isrc_].x_i_2_;
dvbc_2_ = f02_jp_smem_[j_ + jstride_ * 0].v_j_2_ - f02_ip_[isrc_].v_i_2_;
r1b2 = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_;
r1b2e = r1b2 + eps2;
r1be = rsqrt(r1b2e);
r1c2 = dxc_0_ * dxc_0_ + dxc_1_ * dxc_1_ + dxc_2_ * dxc_2_;
r1c2e = r1c2 + eps2;
r1ce = rsqrt(r1c2e);
rbc2 = dxbc_0_ * dxbc_0_ + dxbc_1_ * dxbc_1_ + dxbc_2_ * dxbc_2_;
rbc2e = rbc2 + eps2;
rbce = rsqrt(rbc2e);
mr3inv = f02_jp_smem_[j_ + jstride_ * 0].m_j_ * rbce * rbce * rbce;
ac_i_0_wcache_ += f02_jp_smem_[j_ + jstride_ * 0].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_0_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_0_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 0].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_0_ - f02_jp_smem_[j_ + jstride_ * 0].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 0].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 0].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 0].v_j_2_) * dxb_0_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 0].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_0_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 0].v_j_0_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 0].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 0].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_0_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 0].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 0].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 0].v_j_2_)) * dvbc_0_;
ac_i_1_wcache_ += f02_jp_smem_[j_ + jstride_ * 0].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_1_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_1_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 0].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_1_ - f02_jp_smem_[j_ + jstride_ * 0].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 0].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 0].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 0].v_j_2_) * dxb_1_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 0].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_1_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 0].v_j_1_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 0].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 0].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_1_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 0].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 0].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 0].v_j_2_)) * dvbc_1_;
ac_i_2_wcache_ += f02_jp_smem_[j_ + jstride_ * 0].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_2_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_2_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 0].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_2_ - f02_jp_smem_[j_ + jstride_ * 0].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 0].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 0].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 0].v_j_2_) * dxb_2_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 0].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_2_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 0].v_j_2_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 0].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 0].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 0].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_2_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 0].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 0].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 0].v_j_2_)) * dvbc_2_;
        // loop 1
dxb_0_ = x_0_0_ - f02_ip_[isrc_].x_i_0_;
dxc_0_ = x_0_0_ - f02_jp_smem_[j_ + jstride_ * 1].x_j_0_;
dxbc_0_ = f02_jp_smem_[j_ + jstride_ * 1].x_j_0_ - f02_ip_[isrc_].x_i_0_;
dvbc_0_ = f02_jp_smem_[j_ + jstride_ * 1].v_j_0_ - f02_ip_[isrc_].v_i_0_;
dxb_1_ = x_0_1_ - f02_ip_[isrc_].x_i_1_;
dxc_1_ = x_0_1_ - f02_jp_smem_[j_ + jstride_ * 1].x_j_1_;
dxbc_1_ = f02_jp_smem_[j_ + jstride_ * 1].x_j_1_ - f02_ip_[isrc_].x_i_1_;
dvbc_1_ = f02_jp_smem_[j_ + jstride_ * 1].v_j_1_ - f02_ip_[isrc_].v_i_1_;
dxb_2_ = x_0_2_ - f02_ip_[isrc_].x_i_2_;
dxc_2_ = x_0_2_ - f02_jp_smem_[j_ + jstride_ * 1].x_j_2_;
dxbc_2_ = f02_jp_smem_[j_ + jstride_ * 1].x_j_2_ - f02_ip_[isrc_].x_i_2_;
dvbc_2_ = f02_jp_smem_[j_ + jstride_ * 1].v_j_2_ - f02_ip_[isrc_].v_i_2_;
r1b2 = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_;
r1b2e = r1b2 + eps2;
r1be = rsqrt(r1b2e);
r1c2 = dxc_0_ * dxc_0_ + dxc_1_ * dxc_1_ + dxc_2_ * dxc_2_;
r1c2e = r1c2 + eps2;
r1ce = rsqrt(r1c2e);
rbc2 = dxbc_0_ * dxbc_0_ + dxbc_1_ * dxbc_1_ + dxbc_2_ * dxbc_2_;
rbc2e = rbc2 + eps2;
rbce = rsqrt(rbc2e);
mr3inv = f02_jp_smem_[j_ + jstride_ * 1].m_j_ * rbce * rbce * rbce;
ac_i_0_wcache_ += f02_jp_smem_[j_ + jstride_ * 1].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_0_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_0_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 1].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_0_ - f02_jp_smem_[j_ + jstride_ * 1].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 1].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 1].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 1].v_j_2_) * dxb_0_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 1].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_0_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 1].v_j_0_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 1].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 1].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_0_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 1].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 1].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 1].v_j_2_)) * dvbc_0_;
ac_i_1_wcache_ += f02_jp_smem_[j_ + jstride_ * 1].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_1_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_1_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 1].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_1_ - f02_jp_smem_[j_ + jstride_ * 1].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 1].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 1].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 1].v_j_2_) * dxb_1_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 1].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_1_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 1].v_j_1_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 1].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 1].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_1_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 1].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 1].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 1].v_j_2_)) * dvbc_1_;
ac_i_2_wcache_ += f02_jp_smem_[j_ + jstride_ * 1].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_2_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_2_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 1].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_2_ - f02_jp_smem_[j_ + jstride_ * 1].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 1].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 1].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 1].v_j_2_) * dxb_2_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 1].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_2_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 1].v_j_2_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 1].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 1].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 1].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_2_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 1].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 1].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 1].v_j_2_)) * dvbc_2_;
        // loop 2
dxb_0_ = x_0_0_ - f02_ip_[isrc_].x_i_0_;
dxc_0_ = x_0_0_ - f02_jp_smem_[j_ + jstride_ * 2].x_j_0_;
dxbc_0_ = f02_jp_smem_[j_ + jstride_ * 2].x_j_0_ - f02_ip_[isrc_].x_i_0_;
dvbc_0_ = f02_jp_smem_[j_ + jstride_ * 2].v_j_0_ - f02_ip_[isrc_].v_i_0_;
dxb_1_ = x_0_1_ - f02_ip_[isrc_].x_i_1_;
dxc_1_ = x_0_1_ - f02_jp_smem_[j_ + jstride_ * 2].x_j_1_;
dxbc_1_ = f02_jp_smem_[j_ + jstride_ * 2].x_j_1_ - f02_ip_[isrc_].x_i_1_;
dvbc_1_ = f02_jp_smem_[j_ + jstride_ * 2].v_j_1_ - f02_ip_[isrc_].v_i_1_;
dxb_2_ = x_0_2_ - f02_ip_[isrc_].x_i_2_;
dxc_2_ = x_0_2_ - f02_jp_smem_[j_ + jstride_ * 2].x_j_2_;
dxbc_2_ = f02_jp_smem_[j_ + jstride_ * 2].x_j_2_ - f02_ip_[isrc_].x_i_2_;
dvbc_2_ = f02_jp_smem_[j_ + jstride_ * 2].v_j_2_ - f02_ip_[isrc_].v_i_2_;
r1b2 = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_;
r1b2e = r1b2 + eps2;
r1be = rsqrt(r1b2e);
r1c2 = dxc_0_ * dxc_0_ + dxc_1_ * dxc_1_ + dxc_2_ * dxc_2_;
r1c2e = r1c2 + eps2;
r1ce = rsqrt(r1c2e);
rbc2 = dxbc_0_ * dxbc_0_ + dxbc_1_ * dxbc_1_ + dxbc_2_ * dxbc_2_;
rbc2e = rbc2 + eps2;
rbce = rsqrt(rbc2e);
mr3inv = f02_jp_smem_[j_ + jstride_ * 2].m_j_ * rbce * rbce * rbce;
ac_i_0_wcache_ += f02_jp_smem_[j_ + jstride_ * 2].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_0_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_0_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 2].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_0_ - f02_jp_smem_[j_ + jstride_ * 2].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 2].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 2].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 2].v_j_2_) * dxb_0_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 2].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_0_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 2].v_j_0_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 2].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 2].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_0_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 2].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 2].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 2].v_j_2_)) * dvbc_0_;
ac_i_1_wcache_ += f02_jp_smem_[j_ + jstride_ * 2].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_1_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_1_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 2].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_1_ - f02_jp_smem_[j_ + jstride_ * 2].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 2].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 2].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 2].v_j_2_) * dxb_1_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 2].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_1_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 2].v_j_1_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 2].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 2].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_1_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 2].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 2].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 2].v_j_2_)) * dvbc_1_;
ac_i_2_wcache_ += f02_jp_smem_[j_ + jstride_ * 2].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_2_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_2_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 2].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_2_ - f02_jp_smem_[j_ + jstride_ * 2].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 2].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 2].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 2].v_j_2_) * dxb_2_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 2].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_2_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 2].v_j_2_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 2].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 2].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 2].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_2_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 2].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 2].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 2].v_j_2_)) * dvbc_2_;
        // loop 3
dxb_0_ = x_0_0_ - f02_ip_[isrc_].x_i_0_;
dxc_0_ = x_0_0_ - f02_jp_smem_[j_ + jstride_ * 3].x_j_0_;
dxbc_0_ = f02_jp_smem_[j_ + jstride_ * 3].x_j_0_ - f02_ip_[isrc_].x_i_0_;
dvbc_0_ = f02_jp_smem_[j_ + jstride_ * 3].v_j_0_ - f02_ip_[isrc_].v_i_0_;
dxb_1_ = x_0_1_ - f02_ip_[isrc_].x_i_1_;
dxc_1_ = x_0_1_ - f02_jp_smem_[j_ + jstride_ * 3].x_j_1_;
dxbc_1_ = f02_jp_smem_[j_ + jstride_ * 3].x_j_1_ - f02_ip_[isrc_].x_i_1_;
dvbc_1_ = f02_jp_smem_[j_ + jstride_ * 3].v_j_1_ - f02_ip_[isrc_].v_i_1_;
dxb_2_ = x_0_2_ - f02_ip_[isrc_].x_i_2_;
dxc_2_ = x_0_2_ - f02_jp_smem_[j_ + jstride_ * 3].x_j_2_;
dxbc_2_ = f02_jp_smem_[j_ + jstride_ * 3].x_j_2_ - f02_ip_[isrc_].x_i_2_;
dvbc_2_ = f02_jp_smem_[j_ + jstride_ * 3].v_j_2_ - f02_ip_[isrc_].v_i_2_;
r1b2 = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_;
r1b2e = r1b2 + eps2;
r1be = rsqrt(r1b2e);
r1c2 = dxc_0_ * dxc_0_ + dxc_1_ * dxc_1_ + dxc_2_ * dxc_2_;
r1c2e = r1c2 + eps2;
r1ce = rsqrt(r1c2e);
rbc2 = dxbc_0_ * dxbc_0_ + dxbc_1_ * dxbc_1_ + dxbc_2_ * dxbc_2_;
rbc2e = rbc2 + eps2;
rbce = rsqrt(rbc2e);
mr3inv = f02_jp_smem_[j_ + jstride_ * 3].m_j_ * rbce * rbce * rbce;
ac_i_0_wcache_ += f02_jp_smem_[j_ + jstride_ * 3].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_0_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_0_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 3].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_0_ - f02_jp_smem_[j_ + jstride_ * 3].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 3].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 3].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 3].v_j_2_) * dxb_0_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 3].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_0_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 3].v_j_0_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 3].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 3].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_0_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 3].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 3].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 3].v_j_2_)) * dvbc_0_;
ac_i_1_wcache_ += f02_jp_smem_[j_ + jstride_ * 3].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_1_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_1_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 3].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_1_ - f02_jp_smem_[j_ + jstride_ * 3].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 3].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 3].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 3].v_j_2_) * dxb_1_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 3].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_1_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 3].v_j_1_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 3].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 3].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_1_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 3].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 3].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 3].v_j_2_)) * dvbc_1_;
ac_i_2_wcache_ += f02_jp_smem_[j_ + jstride_ * 3].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_2_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_2_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 3].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_2_ - f02_jp_smem_[j_ + jstride_ * 3].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 3].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 3].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 3].v_j_2_) * dxb_2_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 3].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_2_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 3].v_j_2_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 3].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 3].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 3].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_2_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 3].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 3].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 3].v_j_2_)) * dvbc_2_;
        // loop 4
dxb_0_ = x_0_0_ - f02_ip_[isrc_].x_i_0_;
dxc_0_ = x_0_0_ - f02_jp_smem_[j_ + jstride_ * 4].x_j_0_;
dxbc_0_ = f02_jp_smem_[j_ + jstride_ * 4].x_j_0_ - f02_ip_[isrc_].x_i_0_;
dvbc_0_ = f02_jp_smem_[j_ + jstride_ * 4].v_j_0_ - f02_ip_[isrc_].v_i_0_;
dxb_1_ = x_0_1_ - f02_ip_[isrc_].x_i_1_;
dxc_1_ = x_0_1_ - f02_jp_smem_[j_ + jstride_ * 4].x_j_1_;
dxbc_1_ = f02_jp_smem_[j_ + jstride_ * 4].x_j_1_ - f02_ip_[isrc_].x_i_1_;
dvbc_1_ = f02_jp_smem_[j_ + jstride_ * 4].v_j_1_ - f02_ip_[isrc_].v_i_1_;
dxb_2_ = x_0_2_ - f02_ip_[isrc_].x_i_2_;
dxc_2_ = x_0_2_ - f02_jp_smem_[j_ + jstride_ * 4].x_j_2_;
dxbc_2_ = f02_jp_smem_[j_ + jstride_ * 4].x_j_2_ - f02_ip_[isrc_].x_i_2_;
dvbc_2_ = f02_jp_smem_[j_ + jstride_ * 4].v_j_2_ - f02_ip_[isrc_].v_i_2_;
r1b2 = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_;
r1b2e = r1b2 + eps2;
r1be = rsqrt(r1b2e);
r1c2 = dxc_0_ * dxc_0_ + dxc_1_ * dxc_1_ + dxc_2_ * dxc_2_;
r1c2e = r1c2 + eps2;
r1ce = rsqrt(r1c2e);
rbc2 = dxbc_0_ * dxbc_0_ + dxbc_1_ * dxbc_1_ + dxbc_2_ * dxbc_2_;
rbc2e = rbc2 + eps2;
rbce = rsqrt(rbc2e);
mr3inv = f02_jp_smem_[j_ + jstride_ * 4].m_j_ * rbce * rbce * rbce;
ac_i_0_wcache_ += f02_jp_smem_[j_ + jstride_ * 4].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_0_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_0_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 4].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_0_ - f02_jp_smem_[j_ + jstride_ * 4].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 4].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 4].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 4].v_j_2_) * dxb_0_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 4].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_0_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 4].v_j_0_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 4].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 4].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_0_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 4].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 4].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 4].v_j_2_)) * dvbc_0_;
ac_i_1_wcache_ += f02_jp_smem_[j_ + jstride_ * 4].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_1_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_1_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 4].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_1_ - f02_jp_smem_[j_ + jstride_ * 4].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 4].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 4].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 4].v_j_2_) * dxb_1_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 4].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_1_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 4].v_j_1_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 4].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 4].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_1_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 4].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 4].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 4].v_j_2_)) * dvbc_1_;
ac_i_2_wcache_ += f02_jp_smem_[j_ + jstride_ * 4].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_2_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_2_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 4].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_2_ - f02_jp_smem_[j_ + jstride_ * 4].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 4].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 4].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 4].v_j_2_) * dxb_2_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 4].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_2_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 4].v_j_2_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 4].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 4].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 4].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_2_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 4].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 4].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 4].v_j_2_)) * dvbc_2_;
        // loop 5
dxb_0_ = x_0_0_ - f02_ip_[isrc_].x_i_0_;
dxc_0_ = x_0_0_ - f02_jp_smem_[j_ + jstride_ * 5].x_j_0_;
dxbc_0_ = f02_jp_smem_[j_ + jstride_ * 5].x_j_0_ - f02_ip_[isrc_].x_i_0_;
dvbc_0_ = f02_jp_smem_[j_ + jstride_ * 5].v_j_0_ - f02_ip_[isrc_].v_i_0_;
dxb_1_ = x_0_1_ - f02_ip_[isrc_].x_i_1_;
dxc_1_ = x_0_1_ - f02_jp_smem_[j_ + jstride_ * 5].x_j_1_;
dxbc_1_ = f02_jp_smem_[j_ + jstride_ * 5].x_j_1_ - f02_ip_[isrc_].x_i_1_;
dvbc_1_ = f02_jp_smem_[j_ + jstride_ * 5].v_j_1_ - f02_ip_[isrc_].v_i_1_;
dxb_2_ = x_0_2_ - f02_ip_[isrc_].x_i_2_;
dxc_2_ = x_0_2_ - f02_jp_smem_[j_ + jstride_ * 5].x_j_2_;
dxbc_2_ = f02_jp_smem_[j_ + jstride_ * 5].x_j_2_ - f02_ip_[isrc_].x_i_2_;
dvbc_2_ = f02_jp_smem_[j_ + jstride_ * 5].v_j_2_ - f02_ip_[isrc_].v_i_2_;
r1b2 = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_;
r1b2e = r1b2 + eps2;
r1be = rsqrt(r1b2e);
r1c2 = dxc_0_ * dxc_0_ + dxc_1_ * dxc_1_ + dxc_2_ * dxc_2_;
r1c2e = r1c2 + eps2;
r1ce = rsqrt(r1c2e);
rbc2 = dxbc_0_ * dxbc_0_ + dxbc_1_ * dxbc_1_ + dxbc_2_ * dxbc_2_;
rbc2e = rbc2 + eps2;
rbce = rsqrt(rbc2e);
mr3inv = f02_jp_smem_[j_ + jstride_ * 5].m_j_ * rbce * rbce * rbce;
ac_i_0_wcache_ += f02_jp_smem_[j_ + jstride_ * 5].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_0_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_0_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 5].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_0_ - f02_jp_smem_[j_ + jstride_ * 5].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 5].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 5].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 5].v_j_2_) * dxb_0_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 5].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_0_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 5].v_j_0_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 5].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 5].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_0_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 5].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 5].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 5].v_j_2_)) * dvbc_0_;
ac_i_1_wcache_ += f02_jp_smem_[j_ + jstride_ * 5].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_1_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_1_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 5].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_1_ - f02_jp_smem_[j_ + jstride_ * 5].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 5].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 5].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 5].v_j_2_) * dxb_1_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 5].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_1_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 5].v_j_1_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 5].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 5].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_1_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 5].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 5].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 5].v_j_2_)) * dvbc_1_;
ac_i_2_wcache_ += f02_jp_smem_[j_ + jstride_ * 5].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_2_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_2_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 5].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_2_ - f02_jp_smem_[j_ + jstride_ * 5].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 5].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 5].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 5].v_j_2_) * dxb_2_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 5].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_2_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 5].v_j_2_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 5].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 5].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 5].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_2_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 5].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 5].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 5].v_j_2_)) * dvbc_2_;
        // loop 6
dxb_0_ = x_0_0_ - f02_ip_[isrc_].x_i_0_;
dxc_0_ = x_0_0_ - f02_jp_smem_[j_ + jstride_ * 6].x_j_0_;
dxbc_0_ = f02_jp_smem_[j_ + jstride_ * 6].x_j_0_ - f02_ip_[isrc_].x_i_0_;
dvbc_0_ = f02_jp_smem_[j_ + jstride_ * 6].v_j_0_ - f02_ip_[isrc_].v_i_0_;
dxb_1_ = x_0_1_ - f02_ip_[isrc_].x_i_1_;
dxc_1_ = x_0_1_ - f02_jp_smem_[j_ + jstride_ * 6].x_j_1_;
dxbc_1_ = f02_jp_smem_[j_ + jstride_ * 6].x_j_1_ - f02_ip_[isrc_].x_i_1_;
dvbc_1_ = f02_jp_smem_[j_ + jstride_ * 6].v_j_1_ - f02_ip_[isrc_].v_i_1_;
dxb_2_ = x_0_2_ - f02_ip_[isrc_].x_i_2_;
dxc_2_ = x_0_2_ - f02_jp_smem_[j_ + jstride_ * 6].x_j_2_;
dxbc_2_ = f02_jp_smem_[j_ + jstride_ * 6].x_j_2_ - f02_ip_[isrc_].x_i_2_;
dvbc_2_ = f02_jp_smem_[j_ + jstride_ * 6].v_j_2_ - f02_ip_[isrc_].v_i_2_;
r1b2 = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_;
r1b2e = r1b2 + eps2;
r1be = rsqrt(r1b2e);
r1c2 = dxc_0_ * dxc_0_ + dxc_1_ * dxc_1_ + dxc_2_ * dxc_2_;
r1c2e = r1c2 + eps2;
r1ce = rsqrt(r1c2e);
rbc2 = dxbc_0_ * dxbc_0_ + dxbc_1_ * dxbc_1_ + dxbc_2_ * dxbc_2_;
rbc2e = rbc2 + eps2;
rbce = rsqrt(rbc2e);
mr3inv = f02_jp_smem_[j_ + jstride_ * 6].m_j_ * rbce * rbce * rbce;
ac_i_0_wcache_ += f02_jp_smem_[j_ + jstride_ * 6].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_0_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_0_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 6].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_0_ - f02_jp_smem_[j_ + jstride_ * 6].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 6].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 6].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 6].v_j_2_) * dxb_0_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 6].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_0_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 6].v_j_0_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 6].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 6].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_0_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 6].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 6].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 6].v_j_2_)) * dvbc_0_;
ac_i_1_wcache_ += f02_jp_smem_[j_ + jstride_ * 6].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_1_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_1_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 6].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_1_ - f02_jp_smem_[j_ + jstride_ * 6].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 6].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 6].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 6].v_j_2_) * dxb_1_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 6].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_1_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 6].v_j_1_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 6].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 6].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_1_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 6].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 6].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 6].v_j_2_)) * dvbc_1_;
ac_i_2_wcache_ += f02_jp_smem_[j_ + jstride_ * 6].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_2_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_2_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 6].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_2_ - f02_jp_smem_[j_ + jstride_ * 6].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 6].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 6].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 6].v_j_2_) * dxb_2_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 6].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_2_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 6].v_j_2_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 6].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 6].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 6].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_2_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 6].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 6].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 6].v_j_2_)) * dvbc_2_;
        // loop 7
dxb_0_ = x_0_0_ - f02_ip_[isrc_].x_i_0_;
dxc_0_ = x_0_0_ - f02_jp_smem_[j_ + jstride_ * 7].x_j_0_;
dxbc_0_ = f02_jp_smem_[j_ + jstride_ * 7].x_j_0_ - f02_ip_[isrc_].x_i_0_;
dvbc_0_ = f02_jp_smem_[j_ + jstride_ * 7].v_j_0_ - f02_ip_[isrc_].v_i_0_;
dxb_1_ = x_0_1_ - f02_ip_[isrc_].x_i_1_;
dxc_1_ = x_0_1_ - f02_jp_smem_[j_ + jstride_ * 7].x_j_1_;
dxbc_1_ = f02_jp_smem_[j_ + jstride_ * 7].x_j_1_ - f02_ip_[isrc_].x_i_1_;
dvbc_1_ = f02_jp_smem_[j_ + jstride_ * 7].v_j_1_ - f02_ip_[isrc_].v_i_1_;
dxb_2_ = x_0_2_ - f02_ip_[isrc_].x_i_2_;
dxc_2_ = x_0_2_ - f02_jp_smem_[j_ + jstride_ * 7].x_j_2_;
dxbc_2_ = f02_jp_smem_[j_ + jstride_ * 7].x_j_2_ - f02_ip_[isrc_].x_i_2_;
dvbc_2_ = f02_jp_smem_[j_ + jstride_ * 7].v_j_2_ - f02_ip_[isrc_].v_i_2_;
r1b2 = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_;
r1b2e = r1b2 + eps2;
r1be = rsqrt(r1b2e);
r1c2 = dxc_0_ * dxc_0_ + dxc_1_ * dxc_1_ + dxc_2_ * dxc_2_;
r1c2e = r1c2 + eps2;
r1ce = rsqrt(r1c2e);
rbc2 = dxbc_0_ * dxbc_0_ + dxbc_1_ * dxbc_1_ + dxbc_2_ * dxbc_2_;
rbc2e = rbc2 + eps2;
rbce = rsqrt(rbc2e);
mr3inv = f02_jp_smem_[j_ + jstride_ * 7].m_j_ * rbce * rbce * rbce;
ac_i_0_wcache_ += f02_jp_smem_[j_ + jstride_ * 7].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_0_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_0_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 7].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_0_ - f02_jp_smem_[j_ + jstride_ * 7].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 7].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 7].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 7].v_j_2_) * dxb_0_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 7].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_0_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 7].v_j_0_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 7].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 7].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_0_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 7].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 7].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 7].v_j_2_)) * dvbc_0_;
ac_i_1_wcache_ += f02_jp_smem_[j_ + jstride_ * 7].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_1_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_1_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 7].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_1_ - f02_jp_smem_[j_ + jstride_ * 7].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 7].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 7].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 7].v_j_2_) * dxb_1_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 7].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_1_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 7].v_j_1_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 7].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 7].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_1_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 7].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 7].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 7].v_j_2_)) * dvbc_1_;
ac_i_2_wcache_ += f02_jp_smem_[j_ + jstride_ * 7].m_j_ * r1be * r1be * r1be * m_0_ * (4.0 * rbce + 1.25 * r1ce + 0.25 * (r1b2 - rbc2) * r1ce * r1ce * r1ce) * dxb_2_ + mr3inv * m_0_ * (4.0 * r1be + 1.25 * r1ce + 0.25 * (-r1b2 + rbc2) * r1ce * r1ce * r1ce) * dxbc_2_ - 3.5 * f02_jp_smem_[j_ + jstride_ * 7].m_j_ * m_0_ * r1ce * r1ce * r1ce * (rbce - r1be) * dxc_2_ - f02_jp_smem_[j_ + jstride_ * 7].m_j_ * r1be * r1be * r1be * (4.0 * (f02_ip_[isrc_].v_i_0_ * f02_jp_smem_[j_ + jstride_ * 7].v_j_0_ + f02_ip_[isrc_].v_i_1_ * f02_jp_smem_[j_ + jstride_ * 7].v_j_1_ + f02_ip_[isrc_].v_i_2_ * f02_jp_smem_[j_ + jstride_ * 7].v_j_2_) * dxb_2_ - 3.0 * (f02_jp_smem_[j_ + jstride_ * 7].v_j_0_ * dxb_0_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_1_ * dxb_1_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_2_ * dxb_2_) * f02_ip_[isrc_].v_i_2_ - 4.0 * (f02_ip_[isrc_].v_i_0_ * dxb_0_ + f02_ip_[isrc_].v_i_1_ * dxb_1_ + f02_ip_[isrc_].v_i_2_ * dxb_2_) * f02_jp_smem_[j_ + jstride_ * 7].v_j_2_) + mr3inv * (f02_ip_[isrc_].v_i_0_ * f02_ip_[isrc_].v_i_0_ + f02_ip_[isrc_].v_i_1_ * f02_ip_[isrc_].v_i_1_ + f02_ip_[isrc_].v_i_2_ * f02_ip_[isrc_].v_i_2_ - 2.0 * (dvbc_0_ * dvbc_0_ + dvbc_1_ * dvbc_1_ + dvbc_2_ * dvbc_2_) + 1.5 * (f02_jp_smem_[j_ + jstride_ * 7].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_2_ * dxbc_2_) * (f02_jp_smem_[j_ + jstride_ * 7].v_j_0_ * dxbc_0_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_1_ * dxbc_1_ + f02_jp_smem_[j_ + jstride_ * 7].v_j_2_ * dxbc_2_) / rbc2e) * dxbc_2_ + mr3inv * (dxbc_0_ * (4.0 * f02_ip_[isrc_].v_i_0_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 7].v_j_0_) + dxbc_1_ * (4.0 * f02_ip_[isrc_].v_i_1_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 7].v_j_1_) + dxbc_2_ * (4.0 * f02_ip_[isrc_].v_i_2_ - 3.0 * f02_jp_smem_[j_ + jstride_ * 7].v_j_2_)) * dvbc_2_;

                           }
                       }
                   }
                   __syncthreads();
                   f02_result_smem_[tid_].ac_i_0_ = ac_i_0_wcache_;f02_result_smem_[tid_].ac_i_1_ = ac_i_1_wcache_;f02_result_smem_[tid_].ac_i_2_ = ac_i_2_wcache_;
                   __syncthreads();
                   if (jstride_ > 1) {if (tid_ < kbdim_ / 2) {f02_result_smem_[tid_].ac_i_0_ += f02_result_smem_[tid_ + kbdim_ / 2].ac_i_0_;f02_result_smem_[tid_].ac_i_1_ += f02_result_smem_[tid_ + kbdim_ / 2].ac_i_1_;f02_result_smem_[tid_].ac_i_2_ += f02_result_smem_[tid_ + kbdim_ / 2].ac_i_2_;}// __syncthreads(); // this is not necessary since kbdim_ / 2 <= warp size.
}if (jstride_ > 2) {if (tid_ < kbdim_ / 4) {f02_result_smem_[tid_].ac_i_0_ += f02_result_smem_[tid_ + kbdim_ / 4].ac_i_0_;f02_result_smem_[tid_].ac_i_1_ += f02_result_smem_[tid_ + kbdim_ / 4].ac_i_1_;f02_result_smem_[tid_].ac_i_2_ += f02_result_smem_[tid_ + kbdim_ / 4].ac_i_2_;}// __syncthreads(); // this is not necessary since kbdim_ / 4 <= warp size.
}
                   __syncthreads();

#if 1
                   if (tid_ < nvalidthread_) {
                       int idstoff_ = njdiv_ * kbdim_ * ibid_ + jbid_ + njdiv_ * (tid_ % njhsub_);
                       float4 *srcbuf_ = (float4 *) (f02_result_smem_ + tid_);
                       float4 *dstbuf_ = (float4 *) (f02_result_ + idstoff_);
                       for (int icpy_ = 0; icpy_ < sizeof(f02_result_t) / sizeof(float4); icpy_++) {
                           dstbuf_[icpy_] = srcbuf_[icpy_];
                       }
                   }
#else
                   if (tid_ < nvalidthread_) {
                       f02_result_[idst_] = f02_result_smem_[tid_];
                   }
#endif
               }

/*
 * njdiv_    : # of result fragments per result packet.
 * njdiv_ru_ : njdiv_ rounded up to a power of two.
 * rbdim_    : # of result fragments to be reduced to (rbdim_ / njdiv_) result packets.
 */
               __global__ void
               f02_reducer(int njdiv_, int njdiv_ru_, f02_result_t *f02_result_, f02_result_t *f02_result_sub_)
               {
                   extern __shared__ char smembuf_[];
                   int rbdim_ = blockDim.x;
                   f02_result_t * f02_result_smem_ = (f02_result_t *)smembuf_;
                   f02_result_t * f02_result_smem_packed_ = (f02_result_t *)(smembuf_ + rbdim_ * sizeof(f02_result_t));
                   int tid_ = threadIdx.x;
                   int bid_ = blockIdx.x;
                   int isrc_ = rbdim_ * bid_ + tid_;
                   int ndst_ = rbdim_ / njdiv_;
                   int idst_ = ndst_ * bid_ + tid_;
                   __syncthreads();
                   #if 1
                     {
                         float4 *srcbuf_ = (float4 *)(f02_result_sub_ + rbdim_ * bid_);
                         float4 *dstbuf_ = (float4 *)smembuf_;
                         for (int icpy = 0; icpy < sizeof(f02_result_t)/sizeof(float4); icpy++) {
                             dstbuf_[tid_] = srcbuf_[tid_];
                             dstbuf_ += rbdim_;
                             srcbuf_ += rbdim_;
                         }
                     }
#else
                     f02_result_smem_[tid_] = f02_result_sub_[isrc_];
#endif

                   __syncthreads();

                   int n_ = njdiv_ru_;
                   while (n_ > 1) {
                       n_ /= 2;
                       int ipartner_ = tid_ + n_;
                       if (tid_ % njdiv_ < n_ && ipartner_ % njdiv_ru_ < njdiv_) {
                           f02_result_smem_[tid_].ac_i_0_ += f02_result_smem_[ipartner_].ac_i_0_;f02_result_smem_[tid_].ac_i_1_ += f02_result_smem_[ipartner_].ac_i_1_;f02_result_smem_[tid_].ac_i_2_ += f02_result_smem_[ipartner_].ac_i_2_;
                       }
                       __syncthreads(); // this is not necessary if rbdim_ <= warp size.
                   }
                   __syncthreads();
                   if (tid_ % njdiv_ == 0) {
                       int ipack_ = tid_ / njdiv_;
                       f02_result_smem_packed_[ipack_] = f02_result_smem_[tid_];
                   }
                   __syncthreads();
#if 1
                   {
                       float4 *srcbuf_ = (float4 *) f02_result_smem_packed_;
                       float4 *dstbuf_ = (float4 *) (f02_result_ + ndst_ * bid_);
                       if (tid_ < ndst_) {
                           for (int icpy = 0; icpy < sizeof(f02_result_t) / sizeof(float4); icpy++) {
                               dstbuf_[tid_] = srcbuf_[tid_];
                               dstbuf_ += ndst_;
                               srcbuf_ += ndst_;
                           }
                       }
                   }
#else
                   if (tid_ < ndst_) {
                       f02_result_[idst_] = f02_result_smem_packed_[tid_];
                   }
#endif
               }
