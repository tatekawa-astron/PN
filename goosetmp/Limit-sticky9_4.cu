               __global__ void
               f04_calculator(int ioff_, int ni_, int nj_, double x_0_0_, double x_0_1_, double x_0_2_, double eps2, double m_0_, f04_jp_t *f04_jp_, f04_ip_t *f04_ip_, f04_result_t *f04_result_)
               {
                   extern __shared__ char smembuf_[];
                   int kbdim_ = blockDim.x;
                   
                   f04_jp_t * f04_jp_smem_ = (f04_jp_t *)smembuf_;
                   
                   f04_result_t * f04_result_smem_ = (f04_result_t *)smembuf_;
                   double dxa_0_, dxb_0_, dxab_0_, dxa_1_, dxb_1_, dxab_1_, dxa_2_, dxb_2_, dxab_2_, r1a2e, r1ainv, r1b2e, rab2e, r1binv, rabinv, vi2;
                   double pot_pn2_i_wcache_ = 0.0f;

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
                         float4 *srcbuf_ = (float4 *)(f04_jp_ + joff_);
                         float4 *dstbuf_ = (float4 *)smembuf_;
                         for (int icpy = 0; icpy < sizeof(f04_jp_t)/sizeof(float4); icpy++) {
                             dstbuf_[tid_] = srcbuf_[tid_];
                             dstbuf_ += kbdim_;
                             srcbuf_ += kbdim_;
                         }
                     }
#else
                     f04_jp_smem_[tid_] = f04_jp_[jsrc_];
#endif

                       __syncthreads();
                       int jsup_ = kbdim_;
                       if (joff_ + jsup_ > joff1_) {
                           jsup_ = joff1_ - joff_;
                       }
                       if (jsup_ < kbdim_) {
                           for (int j_ = jstart_; j_ < jsup_; j_+= jstride_) {
                               dxa_0_ = f04_ip_[isrc_].x_i_0_ - x_0_0_;
dxb_0_ = f04_jp_smem_[j_].x_j_0_ - x_0_0_;
dxab_0_ = f04_jp_smem_[j_].x_j_0_ - f04_ip_[isrc_].x_i_0_;
dxa_1_ = f04_ip_[isrc_].x_i_1_ - x_0_1_;
dxb_1_ = f04_jp_smem_[j_].x_j_1_ - x_0_1_;
dxab_1_ = f04_jp_smem_[j_].x_j_1_ - f04_ip_[isrc_].x_i_1_;
dxa_2_ = f04_ip_[isrc_].x_i_2_ - x_0_2_;
dxb_2_ = f04_jp_smem_[j_].x_j_2_ - x_0_2_;
dxab_2_ = f04_jp_smem_[j_].x_j_2_ - f04_ip_[isrc_].x_i_2_;
r1a2e = dxa_0_ * dxa_0_ + dxa_1_ * dxa_1_ + dxa_2_ * dxa_2_ + eps2;
r1ainv = rsqrt(r1a2e);
r1b2e = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_ + eps2;
rab2e = dxab_0_ * dxab_0_ + dxab_1_ * dxab_1_ + dxab_2_ * dxab_2_ + eps2;
r1binv = rsqrt(r1b2e);
rabinv = rsqrt(rab2e);
vi2 = f04_ip_[isrc_].v_i_0_ * f04_ip_[isrc_].v_i_0_ + f04_ip_[isrc_].v_i_1_ * f04_ip_[isrc_].v_i_1_ + f04_ip_[isrc_].v_i_2_ * f04_ip_[isrc_].v_i_2_;
pot_pn2_i_wcache_ += 0.25 * rabinv * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_].m_j_ * (6.0 * vi2 - 7.0 * (f04_ip_[isrc_].v_i_0_ * f04_jp_smem_[j_].v_j_0_ + f04_ip_[isrc_].v_i_1_ * f04_jp_smem_[j_].v_j_1_ + f04_ip_[isrc_].v_i_2_ * f04_jp_smem_[j_].v_j_2_) - (dxab_0_ * f04_ip_[isrc_].v_i_0_ + dxab_1_ * f04_ip_[isrc_].v_i_1_ + dxab_2_ * f04_ip_[isrc_].v_i_2_) * (dxab_0_ * f04_jp_smem_[j_].v_j_0_ + dxab_1_ * f04_jp_smem_[j_].v_j_1_ + dxab_2_ * f04_jp_smem_[j_].v_j_2_) / rab2e) + m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_].m_j_ * r1ainv * rabinv + 0.5 * m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_].m_j_ * r1ainv * r1binv;

                           }
                       }
                       else {    
                           for (int j_ = jstart_; j_ < kbdim_; j_+= jstride_ * 8) {
                                       // loop 0
dxa_0_ = f04_ip_[isrc_].x_i_0_ - x_0_0_;
dxb_0_ = f04_jp_smem_[j_ + jstride_ * 0].x_j_0_ - x_0_0_;
dxab_0_ = f04_jp_smem_[j_ + jstride_ * 0].x_j_0_ - f04_ip_[isrc_].x_i_0_;
dxa_1_ = f04_ip_[isrc_].x_i_1_ - x_0_1_;
dxb_1_ = f04_jp_smem_[j_ + jstride_ * 0].x_j_1_ - x_0_1_;
dxab_1_ = f04_jp_smem_[j_ + jstride_ * 0].x_j_1_ - f04_ip_[isrc_].x_i_1_;
dxa_2_ = f04_ip_[isrc_].x_i_2_ - x_0_2_;
dxb_2_ = f04_jp_smem_[j_ + jstride_ * 0].x_j_2_ - x_0_2_;
dxab_2_ = f04_jp_smem_[j_ + jstride_ * 0].x_j_2_ - f04_ip_[isrc_].x_i_2_;
r1a2e = dxa_0_ * dxa_0_ + dxa_1_ * dxa_1_ + dxa_2_ * dxa_2_ + eps2;
r1ainv = rsqrt(r1a2e);
r1b2e = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_ + eps2;
rab2e = dxab_0_ * dxab_0_ + dxab_1_ * dxab_1_ + dxab_2_ * dxab_2_ + eps2;
r1binv = rsqrt(r1b2e);
rabinv = rsqrt(rab2e);
vi2 = f04_ip_[isrc_].v_i_0_ * f04_ip_[isrc_].v_i_0_ + f04_ip_[isrc_].v_i_1_ * f04_ip_[isrc_].v_i_1_ + f04_ip_[isrc_].v_i_2_ * f04_ip_[isrc_].v_i_2_;
pot_pn2_i_wcache_ += 0.25 * rabinv * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 0].m_j_ * (6.0 * vi2 - 7.0 * (f04_ip_[isrc_].v_i_0_ * f04_jp_smem_[j_ + jstride_ * 0].v_j_0_ + f04_ip_[isrc_].v_i_1_ * f04_jp_smem_[j_ + jstride_ * 0].v_j_1_ + f04_ip_[isrc_].v_i_2_ * f04_jp_smem_[j_ + jstride_ * 0].v_j_2_) - (dxab_0_ * f04_ip_[isrc_].v_i_0_ + dxab_1_ * f04_ip_[isrc_].v_i_1_ + dxab_2_ * f04_ip_[isrc_].v_i_2_) * (dxab_0_ * f04_jp_smem_[j_ + jstride_ * 0].v_j_0_ + dxab_1_ * f04_jp_smem_[j_ + jstride_ * 0].v_j_1_ + dxab_2_ * f04_jp_smem_[j_ + jstride_ * 0].v_j_2_) / rab2e) + m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 0].m_j_ * r1ainv * rabinv + 0.5 * m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 0].m_j_ * r1ainv * r1binv;
        // loop 1
dxa_0_ = f04_ip_[isrc_].x_i_0_ - x_0_0_;
dxb_0_ = f04_jp_smem_[j_ + jstride_ * 1].x_j_0_ - x_0_0_;
dxab_0_ = f04_jp_smem_[j_ + jstride_ * 1].x_j_0_ - f04_ip_[isrc_].x_i_0_;
dxa_1_ = f04_ip_[isrc_].x_i_1_ - x_0_1_;
dxb_1_ = f04_jp_smem_[j_ + jstride_ * 1].x_j_1_ - x_0_1_;
dxab_1_ = f04_jp_smem_[j_ + jstride_ * 1].x_j_1_ - f04_ip_[isrc_].x_i_1_;
dxa_2_ = f04_ip_[isrc_].x_i_2_ - x_0_2_;
dxb_2_ = f04_jp_smem_[j_ + jstride_ * 1].x_j_2_ - x_0_2_;
dxab_2_ = f04_jp_smem_[j_ + jstride_ * 1].x_j_2_ - f04_ip_[isrc_].x_i_2_;
r1a2e = dxa_0_ * dxa_0_ + dxa_1_ * dxa_1_ + dxa_2_ * dxa_2_ + eps2;
r1ainv = rsqrt(r1a2e);
r1b2e = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_ + eps2;
rab2e = dxab_0_ * dxab_0_ + dxab_1_ * dxab_1_ + dxab_2_ * dxab_2_ + eps2;
r1binv = rsqrt(r1b2e);
rabinv = rsqrt(rab2e);
vi2 = f04_ip_[isrc_].v_i_0_ * f04_ip_[isrc_].v_i_0_ + f04_ip_[isrc_].v_i_1_ * f04_ip_[isrc_].v_i_1_ + f04_ip_[isrc_].v_i_2_ * f04_ip_[isrc_].v_i_2_;
pot_pn2_i_wcache_ += 0.25 * rabinv * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 1].m_j_ * (6.0 * vi2 - 7.0 * (f04_ip_[isrc_].v_i_0_ * f04_jp_smem_[j_ + jstride_ * 1].v_j_0_ + f04_ip_[isrc_].v_i_1_ * f04_jp_smem_[j_ + jstride_ * 1].v_j_1_ + f04_ip_[isrc_].v_i_2_ * f04_jp_smem_[j_ + jstride_ * 1].v_j_2_) - (dxab_0_ * f04_ip_[isrc_].v_i_0_ + dxab_1_ * f04_ip_[isrc_].v_i_1_ + dxab_2_ * f04_ip_[isrc_].v_i_2_) * (dxab_0_ * f04_jp_smem_[j_ + jstride_ * 1].v_j_0_ + dxab_1_ * f04_jp_smem_[j_ + jstride_ * 1].v_j_1_ + dxab_2_ * f04_jp_smem_[j_ + jstride_ * 1].v_j_2_) / rab2e) + m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 1].m_j_ * r1ainv * rabinv + 0.5 * m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 1].m_j_ * r1ainv * r1binv;
        // loop 2
dxa_0_ = f04_ip_[isrc_].x_i_0_ - x_0_0_;
dxb_0_ = f04_jp_smem_[j_ + jstride_ * 2].x_j_0_ - x_0_0_;
dxab_0_ = f04_jp_smem_[j_ + jstride_ * 2].x_j_0_ - f04_ip_[isrc_].x_i_0_;
dxa_1_ = f04_ip_[isrc_].x_i_1_ - x_0_1_;
dxb_1_ = f04_jp_smem_[j_ + jstride_ * 2].x_j_1_ - x_0_1_;
dxab_1_ = f04_jp_smem_[j_ + jstride_ * 2].x_j_1_ - f04_ip_[isrc_].x_i_1_;
dxa_2_ = f04_ip_[isrc_].x_i_2_ - x_0_2_;
dxb_2_ = f04_jp_smem_[j_ + jstride_ * 2].x_j_2_ - x_0_2_;
dxab_2_ = f04_jp_smem_[j_ + jstride_ * 2].x_j_2_ - f04_ip_[isrc_].x_i_2_;
r1a2e = dxa_0_ * dxa_0_ + dxa_1_ * dxa_1_ + dxa_2_ * dxa_2_ + eps2;
r1ainv = rsqrt(r1a2e);
r1b2e = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_ + eps2;
rab2e = dxab_0_ * dxab_0_ + dxab_1_ * dxab_1_ + dxab_2_ * dxab_2_ + eps2;
r1binv = rsqrt(r1b2e);
rabinv = rsqrt(rab2e);
vi2 = f04_ip_[isrc_].v_i_0_ * f04_ip_[isrc_].v_i_0_ + f04_ip_[isrc_].v_i_1_ * f04_ip_[isrc_].v_i_1_ + f04_ip_[isrc_].v_i_2_ * f04_ip_[isrc_].v_i_2_;
pot_pn2_i_wcache_ += 0.25 * rabinv * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 2].m_j_ * (6.0 * vi2 - 7.0 * (f04_ip_[isrc_].v_i_0_ * f04_jp_smem_[j_ + jstride_ * 2].v_j_0_ + f04_ip_[isrc_].v_i_1_ * f04_jp_smem_[j_ + jstride_ * 2].v_j_1_ + f04_ip_[isrc_].v_i_2_ * f04_jp_smem_[j_ + jstride_ * 2].v_j_2_) - (dxab_0_ * f04_ip_[isrc_].v_i_0_ + dxab_1_ * f04_ip_[isrc_].v_i_1_ + dxab_2_ * f04_ip_[isrc_].v_i_2_) * (dxab_0_ * f04_jp_smem_[j_ + jstride_ * 2].v_j_0_ + dxab_1_ * f04_jp_smem_[j_ + jstride_ * 2].v_j_1_ + dxab_2_ * f04_jp_smem_[j_ + jstride_ * 2].v_j_2_) / rab2e) + m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 2].m_j_ * r1ainv * rabinv + 0.5 * m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 2].m_j_ * r1ainv * r1binv;
        // loop 3
dxa_0_ = f04_ip_[isrc_].x_i_0_ - x_0_0_;
dxb_0_ = f04_jp_smem_[j_ + jstride_ * 3].x_j_0_ - x_0_0_;
dxab_0_ = f04_jp_smem_[j_ + jstride_ * 3].x_j_0_ - f04_ip_[isrc_].x_i_0_;
dxa_1_ = f04_ip_[isrc_].x_i_1_ - x_0_1_;
dxb_1_ = f04_jp_smem_[j_ + jstride_ * 3].x_j_1_ - x_0_1_;
dxab_1_ = f04_jp_smem_[j_ + jstride_ * 3].x_j_1_ - f04_ip_[isrc_].x_i_1_;
dxa_2_ = f04_ip_[isrc_].x_i_2_ - x_0_2_;
dxb_2_ = f04_jp_smem_[j_ + jstride_ * 3].x_j_2_ - x_0_2_;
dxab_2_ = f04_jp_smem_[j_ + jstride_ * 3].x_j_2_ - f04_ip_[isrc_].x_i_2_;
r1a2e = dxa_0_ * dxa_0_ + dxa_1_ * dxa_1_ + dxa_2_ * dxa_2_ + eps2;
r1ainv = rsqrt(r1a2e);
r1b2e = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_ + eps2;
rab2e = dxab_0_ * dxab_0_ + dxab_1_ * dxab_1_ + dxab_2_ * dxab_2_ + eps2;
r1binv = rsqrt(r1b2e);
rabinv = rsqrt(rab2e);
vi2 = f04_ip_[isrc_].v_i_0_ * f04_ip_[isrc_].v_i_0_ + f04_ip_[isrc_].v_i_1_ * f04_ip_[isrc_].v_i_1_ + f04_ip_[isrc_].v_i_2_ * f04_ip_[isrc_].v_i_2_;
pot_pn2_i_wcache_ += 0.25 * rabinv * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 3].m_j_ * (6.0 * vi2 - 7.0 * (f04_ip_[isrc_].v_i_0_ * f04_jp_smem_[j_ + jstride_ * 3].v_j_0_ + f04_ip_[isrc_].v_i_1_ * f04_jp_smem_[j_ + jstride_ * 3].v_j_1_ + f04_ip_[isrc_].v_i_2_ * f04_jp_smem_[j_ + jstride_ * 3].v_j_2_) - (dxab_0_ * f04_ip_[isrc_].v_i_0_ + dxab_1_ * f04_ip_[isrc_].v_i_1_ + dxab_2_ * f04_ip_[isrc_].v_i_2_) * (dxab_0_ * f04_jp_smem_[j_ + jstride_ * 3].v_j_0_ + dxab_1_ * f04_jp_smem_[j_ + jstride_ * 3].v_j_1_ + dxab_2_ * f04_jp_smem_[j_ + jstride_ * 3].v_j_2_) / rab2e) + m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 3].m_j_ * r1ainv * rabinv + 0.5 * m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 3].m_j_ * r1ainv * r1binv;
        // loop 4
dxa_0_ = f04_ip_[isrc_].x_i_0_ - x_0_0_;
dxb_0_ = f04_jp_smem_[j_ + jstride_ * 4].x_j_0_ - x_0_0_;
dxab_0_ = f04_jp_smem_[j_ + jstride_ * 4].x_j_0_ - f04_ip_[isrc_].x_i_0_;
dxa_1_ = f04_ip_[isrc_].x_i_1_ - x_0_1_;
dxb_1_ = f04_jp_smem_[j_ + jstride_ * 4].x_j_1_ - x_0_1_;
dxab_1_ = f04_jp_smem_[j_ + jstride_ * 4].x_j_1_ - f04_ip_[isrc_].x_i_1_;
dxa_2_ = f04_ip_[isrc_].x_i_2_ - x_0_2_;
dxb_2_ = f04_jp_smem_[j_ + jstride_ * 4].x_j_2_ - x_0_2_;
dxab_2_ = f04_jp_smem_[j_ + jstride_ * 4].x_j_2_ - f04_ip_[isrc_].x_i_2_;
r1a2e = dxa_0_ * dxa_0_ + dxa_1_ * dxa_1_ + dxa_2_ * dxa_2_ + eps2;
r1ainv = rsqrt(r1a2e);
r1b2e = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_ + eps2;
rab2e = dxab_0_ * dxab_0_ + dxab_1_ * dxab_1_ + dxab_2_ * dxab_2_ + eps2;
r1binv = rsqrt(r1b2e);
rabinv = rsqrt(rab2e);
vi2 = f04_ip_[isrc_].v_i_0_ * f04_ip_[isrc_].v_i_0_ + f04_ip_[isrc_].v_i_1_ * f04_ip_[isrc_].v_i_1_ + f04_ip_[isrc_].v_i_2_ * f04_ip_[isrc_].v_i_2_;
pot_pn2_i_wcache_ += 0.25 * rabinv * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 4].m_j_ * (6.0 * vi2 - 7.0 * (f04_ip_[isrc_].v_i_0_ * f04_jp_smem_[j_ + jstride_ * 4].v_j_0_ + f04_ip_[isrc_].v_i_1_ * f04_jp_smem_[j_ + jstride_ * 4].v_j_1_ + f04_ip_[isrc_].v_i_2_ * f04_jp_smem_[j_ + jstride_ * 4].v_j_2_) - (dxab_0_ * f04_ip_[isrc_].v_i_0_ + dxab_1_ * f04_ip_[isrc_].v_i_1_ + dxab_2_ * f04_ip_[isrc_].v_i_2_) * (dxab_0_ * f04_jp_smem_[j_ + jstride_ * 4].v_j_0_ + dxab_1_ * f04_jp_smem_[j_ + jstride_ * 4].v_j_1_ + dxab_2_ * f04_jp_smem_[j_ + jstride_ * 4].v_j_2_) / rab2e) + m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 4].m_j_ * r1ainv * rabinv + 0.5 * m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 4].m_j_ * r1ainv * r1binv;
        // loop 5
dxa_0_ = f04_ip_[isrc_].x_i_0_ - x_0_0_;
dxb_0_ = f04_jp_smem_[j_ + jstride_ * 5].x_j_0_ - x_0_0_;
dxab_0_ = f04_jp_smem_[j_ + jstride_ * 5].x_j_0_ - f04_ip_[isrc_].x_i_0_;
dxa_1_ = f04_ip_[isrc_].x_i_1_ - x_0_1_;
dxb_1_ = f04_jp_smem_[j_ + jstride_ * 5].x_j_1_ - x_0_1_;
dxab_1_ = f04_jp_smem_[j_ + jstride_ * 5].x_j_1_ - f04_ip_[isrc_].x_i_1_;
dxa_2_ = f04_ip_[isrc_].x_i_2_ - x_0_2_;
dxb_2_ = f04_jp_smem_[j_ + jstride_ * 5].x_j_2_ - x_0_2_;
dxab_2_ = f04_jp_smem_[j_ + jstride_ * 5].x_j_2_ - f04_ip_[isrc_].x_i_2_;
r1a2e = dxa_0_ * dxa_0_ + dxa_1_ * dxa_1_ + dxa_2_ * dxa_2_ + eps2;
r1ainv = rsqrt(r1a2e);
r1b2e = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_ + eps2;
rab2e = dxab_0_ * dxab_0_ + dxab_1_ * dxab_1_ + dxab_2_ * dxab_2_ + eps2;
r1binv = rsqrt(r1b2e);
rabinv = rsqrt(rab2e);
vi2 = f04_ip_[isrc_].v_i_0_ * f04_ip_[isrc_].v_i_0_ + f04_ip_[isrc_].v_i_1_ * f04_ip_[isrc_].v_i_1_ + f04_ip_[isrc_].v_i_2_ * f04_ip_[isrc_].v_i_2_;
pot_pn2_i_wcache_ += 0.25 * rabinv * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 5].m_j_ * (6.0 * vi2 - 7.0 * (f04_ip_[isrc_].v_i_0_ * f04_jp_smem_[j_ + jstride_ * 5].v_j_0_ + f04_ip_[isrc_].v_i_1_ * f04_jp_smem_[j_ + jstride_ * 5].v_j_1_ + f04_ip_[isrc_].v_i_2_ * f04_jp_smem_[j_ + jstride_ * 5].v_j_2_) - (dxab_0_ * f04_ip_[isrc_].v_i_0_ + dxab_1_ * f04_ip_[isrc_].v_i_1_ + dxab_2_ * f04_ip_[isrc_].v_i_2_) * (dxab_0_ * f04_jp_smem_[j_ + jstride_ * 5].v_j_0_ + dxab_1_ * f04_jp_smem_[j_ + jstride_ * 5].v_j_1_ + dxab_2_ * f04_jp_smem_[j_ + jstride_ * 5].v_j_2_) / rab2e) + m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 5].m_j_ * r1ainv * rabinv + 0.5 * m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 5].m_j_ * r1ainv * r1binv;
        // loop 6
dxa_0_ = f04_ip_[isrc_].x_i_0_ - x_0_0_;
dxb_0_ = f04_jp_smem_[j_ + jstride_ * 6].x_j_0_ - x_0_0_;
dxab_0_ = f04_jp_smem_[j_ + jstride_ * 6].x_j_0_ - f04_ip_[isrc_].x_i_0_;
dxa_1_ = f04_ip_[isrc_].x_i_1_ - x_0_1_;
dxb_1_ = f04_jp_smem_[j_ + jstride_ * 6].x_j_1_ - x_0_1_;
dxab_1_ = f04_jp_smem_[j_ + jstride_ * 6].x_j_1_ - f04_ip_[isrc_].x_i_1_;
dxa_2_ = f04_ip_[isrc_].x_i_2_ - x_0_2_;
dxb_2_ = f04_jp_smem_[j_ + jstride_ * 6].x_j_2_ - x_0_2_;
dxab_2_ = f04_jp_smem_[j_ + jstride_ * 6].x_j_2_ - f04_ip_[isrc_].x_i_2_;
r1a2e = dxa_0_ * dxa_0_ + dxa_1_ * dxa_1_ + dxa_2_ * dxa_2_ + eps2;
r1ainv = rsqrt(r1a2e);
r1b2e = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_ + eps2;
rab2e = dxab_0_ * dxab_0_ + dxab_1_ * dxab_1_ + dxab_2_ * dxab_2_ + eps2;
r1binv = rsqrt(r1b2e);
rabinv = rsqrt(rab2e);
vi2 = f04_ip_[isrc_].v_i_0_ * f04_ip_[isrc_].v_i_0_ + f04_ip_[isrc_].v_i_1_ * f04_ip_[isrc_].v_i_1_ + f04_ip_[isrc_].v_i_2_ * f04_ip_[isrc_].v_i_2_;
pot_pn2_i_wcache_ += 0.25 * rabinv * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 6].m_j_ * (6.0 * vi2 - 7.0 * (f04_ip_[isrc_].v_i_0_ * f04_jp_smem_[j_ + jstride_ * 6].v_j_0_ + f04_ip_[isrc_].v_i_1_ * f04_jp_smem_[j_ + jstride_ * 6].v_j_1_ + f04_ip_[isrc_].v_i_2_ * f04_jp_smem_[j_ + jstride_ * 6].v_j_2_) - (dxab_0_ * f04_ip_[isrc_].v_i_0_ + dxab_1_ * f04_ip_[isrc_].v_i_1_ + dxab_2_ * f04_ip_[isrc_].v_i_2_) * (dxab_0_ * f04_jp_smem_[j_ + jstride_ * 6].v_j_0_ + dxab_1_ * f04_jp_smem_[j_ + jstride_ * 6].v_j_1_ + dxab_2_ * f04_jp_smem_[j_ + jstride_ * 6].v_j_2_) / rab2e) + m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 6].m_j_ * r1ainv * rabinv + 0.5 * m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 6].m_j_ * r1ainv * r1binv;
        // loop 7
dxa_0_ = f04_ip_[isrc_].x_i_0_ - x_0_0_;
dxb_0_ = f04_jp_smem_[j_ + jstride_ * 7].x_j_0_ - x_0_0_;
dxab_0_ = f04_jp_smem_[j_ + jstride_ * 7].x_j_0_ - f04_ip_[isrc_].x_i_0_;
dxa_1_ = f04_ip_[isrc_].x_i_1_ - x_0_1_;
dxb_1_ = f04_jp_smem_[j_ + jstride_ * 7].x_j_1_ - x_0_1_;
dxab_1_ = f04_jp_smem_[j_ + jstride_ * 7].x_j_1_ - f04_ip_[isrc_].x_i_1_;
dxa_2_ = f04_ip_[isrc_].x_i_2_ - x_0_2_;
dxb_2_ = f04_jp_smem_[j_ + jstride_ * 7].x_j_2_ - x_0_2_;
dxab_2_ = f04_jp_smem_[j_ + jstride_ * 7].x_j_2_ - f04_ip_[isrc_].x_i_2_;
r1a2e = dxa_0_ * dxa_0_ + dxa_1_ * dxa_1_ + dxa_2_ * dxa_2_ + eps2;
r1ainv = rsqrt(r1a2e);
r1b2e = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_ + eps2;
rab2e = dxab_0_ * dxab_0_ + dxab_1_ * dxab_1_ + dxab_2_ * dxab_2_ + eps2;
r1binv = rsqrt(r1b2e);
rabinv = rsqrt(rab2e);
vi2 = f04_ip_[isrc_].v_i_0_ * f04_ip_[isrc_].v_i_0_ + f04_ip_[isrc_].v_i_1_ * f04_ip_[isrc_].v_i_1_ + f04_ip_[isrc_].v_i_2_ * f04_ip_[isrc_].v_i_2_;
pot_pn2_i_wcache_ += 0.25 * rabinv * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 7].m_j_ * (6.0 * vi2 - 7.0 * (f04_ip_[isrc_].v_i_0_ * f04_jp_smem_[j_ + jstride_ * 7].v_j_0_ + f04_ip_[isrc_].v_i_1_ * f04_jp_smem_[j_ + jstride_ * 7].v_j_1_ + f04_ip_[isrc_].v_i_2_ * f04_jp_smem_[j_ + jstride_ * 7].v_j_2_) - (dxab_0_ * f04_ip_[isrc_].v_i_0_ + dxab_1_ * f04_ip_[isrc_].v_i_1_ + dxab_2_ * f04_ip_[isrc_].v_i_2_) * (dxab_0_ * f04_jp_smem_[j_ + jstride_ * 7].v_j_0_ + dxab_1_ * f04_jp_smem_[j_ + jstride_ * 7].v_j_1_ + dxab_2_ * f04_jp_smem_[j_ + jstride_ * 7].v_j_2_) / rab2e) + m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 7].m_j_ * r1ainv * rabinv + 0.5 * m_0_ * f04_ip_[isrc_].m_i_ * f04_jp_smem_[j_ + jstride_ * 7].m_j_ * r1ainv * r1binv;

                           }
                       }
                   }
                   __syncthreads();
                   f04_result_smem_[tid_].pot_pn2_i_ = pot_pn2_i_wcache_;
                   __syncthreads();
                   if (jstride_ > 1) {if (tid_ < kbdim_ / 2) {f04_result_smem_[tid_].pot_pn2_i_ += f04_result_smem_[tid_ + kbdim_ / 2].pot_pn2_i_;}// __syncthreads(); // this is not necessary since kbdim_ / 2 <= warp size.
}if (jstride_ > 2) {if (tid_ < kbdim_ / 4) {f04_result_smem_[tid_].pot_pn2_i_ += f04_result_smem_[tid_ + kbdim_ / 4].pot_pn2_i_;}// __syncthreads(); // this is not necessary since kbdim_ / 4 <= warp size.
}
                   __syncthreads();

#if 1
                   if (tid_ < nvalidthread_) {
                       int idstoff_ = njdiv_ * kbdim_ * ibid_ + jbid_ + njdiv_ * (tid_ % njhsub_);
                       float4 *srcbuf_ = (float4 *) (f04_result_smem_ + tid_);
                       float4 *dstbuf_ = (float4 *) (f04_result_ + idstoff_);
                       for (int icpy_ = 0; icpy_ < sizeof(f04_result_t) / sizeof(float4); icpy_++) {
                           dstbuf_[icpy_] = srcbuf_[icpy_];
                       }
                   }
#else
                   if (tid_ < nvalidthread_) {
                       f04_result_[idst_] = f04_result_smem_[tid_];
                   }
#endif
               }

/*
 * njdiv_    : # of result fragments per result packet.
 * njdiv_ru_ : njdiv_ rounded up to a power of two.
 * rbdim_    : # of result fragments to be reduced to (rbdim_ / njdiv_) result packets.
 */
               __global__ void
               f04_reducer(int njdiv_, int njdiv_ru_, f04_result_t *f04_result_, f04_result_t *f04_result_sub_)
               {
                   extern __shared__ char smembuf_[];
                   int rbdim_ = blockDim.x;
                   f04_result_t * f04_result_smem_ = (f04_result_t *)smembuf_;
                   f04_result_t * f04_result_smem_packed_ = (f04_result_t *)(smembuf_ + rbdim_ * sizeof(f04_result_t));
                   int tid_ = threadIdx.x;
                   int bid_ = blockIdx.x;
                   int isrc_ = rbdim_ * bid_ + tid_;
                   int ndst_ = rbdim_ / njdiv_;
                   int idst_ = ndst_ * bid_ + tid_;
                   __syncthreads();
                   #if 1
                     {
                         float4 *srcbuf_ = (float4 *)(f04_result_sub_ + rbdim_ * bid_);
                         float4 *dstbuf_ = (float4 *)smembuf_;
                         for (int icpy = 0; icpy < sizeof(f04_result_t)/sizeof(float4); icpy++) {
                             dstbuf_[tid_] = srcbuf_[tid_];
                             dstbuf_ += rbdim_;
                             srcbuf_ += rbdim_;
                         }
                     }
#else
                     f04_result_smem_[tid_] = f04_result_sub_[isrc_];
#endif

                   __syncthreads();

                   int n_ = njdiv_ru_;
                   while (n_ > 1) {
                       n_ /= 2;
                       int ipartner_ = tid_ + n_;
                       if (tid_ % njdiv_ < n_ && ipartner_ % njdiv_ru_ < njdiv_) {
                           f04_result_smem_[tid_].pot_pn2_i_ += f04_result_smem_[ipartner_].pot_pn2_i_;
                       }
                       __syncthreads(); // this is not necessary if rbdim_ <= warp size.
                   }
                   __syncthreads();
                   if (tid_ % njdiv_ == 0) {
                       int ipack_ = tid_ / njdiv_;
                       f04_result_smem_packed_[ipack_] = f04_result_smem_[tid_];
                   }
                   __syncthreads();
#if 1
                   {
                       float4 *srcbuf_ = (float4 *) f04_result_smem_packed_;
                       float4 *dstbuf_ = (float4 *) (f04_result_ + ndst_ * bid_);
                       if (tid_ < ndst_) {
                           for (int icpy = 0; icpy < sizeof(f04_result_t) / sizeof(float4); icpy++) {
                               dstbuf_[tid_] = srcbuf_[tid_];
                               dstbuf_ += ndst_;
                               srcbuf_ += ndst_;
                           }
                       }
                   }
#else
                   if (tid_ < ndst_) {
                       f04_result_[idst_] = f04_result_smem_packed_[tid_];
                   }
#endif
               }
