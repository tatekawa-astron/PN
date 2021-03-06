               __global__ void
               f00_calculator(int ioff_, int ni_, int nj_, double eps2, f00_jp_t *f00_jp_, f00_ip_t *f00_ip_, f00_result_t *f00_result_)
               {
                   extern __shared__ char smembuf_[];
                   int kbdim_ = blockDim.x;
                   
                   f00_jp_t * f00_jp_smem_ = (f00_jp_t *)smembuf_;
                   
                   f00_result_t * f00_result_smem_ = (f00_result_t *)smembuf_;
                   double dx_0_, dx_1_, dx_2_, r2, rinv, mrinv, mr3inv;
                   double a_i_0_wcache_ = 0.0f;
double a_i_1_wcache_ = 0.0f;
double a_i_2_wcache_ = 0.0f;

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
                         float4 *srcbuf_ = (float4 *)(f00_jp_ + joff_);
                         float4 *dstbuf_ = (float4 *)smembuf_;
                         for (int icpy = 0; icpy < sizeof(f00_jp_t)/sizeof(float4); icpy++) {
                             dstbuf_[tid_] = srcbuf_[tid_];
                             dstbuf_ += kbdim_;
                             srcbuf_ += kbdim_;
                         }
                     }
#else
                     f00_jp_smem_[tid_] = f00_jp_[jsrc_];
#endif

                       __syncthreads();
                       int jsup_ = kbdim_;
                       if (joff_ + jsup_ > joff1_) {
                           jsup_ = joff1_ - joff_;
                       }
                       if (jsup_ < kbdim_) {
                           for (int j_ = jstart_; j_ < jsup_; j_+= jstride_) {
                               dx_0_ = f00_jp_smem_[j_].x_j_0_ - f00_ip_[isrc_].x_i_0_;
dx_1_ = f00_jp_smem_[j_].x_j_1_ - f00_ip_[isrc_].x_i_1_;
dx_2_ = f00_jp_smem_[j_].x_j_2_ - f00_ip_[isrc_].x_i_2_;
r2 = dx_0_ * dx_0_ + dx_1_ * dx_1_ + dx_2_ * dx_2_ + eps2;
rinv = rsqrt(r2);
mrinv = f00_jp_smem_[j_].m_j_ * rinv;
mr3inv = mrinv * rinv * rinv;
a_i_0_wcache_ += mr3inv * dx_0_;
a_i_1_wcache_ += mr3inv * dx_1_;
a_i_2_wcache_ += mr3inv * dx_2_;

                           }
                       }
                       else {    
                           for (int j_ = jstart_; j_ < kbdim_; j_+= jstride_ * 8) {
                                       // loop 0
dx_0_ = f00_jp_smem_[j_ + jstride_ * 0].x_j_0_ - f00_ip_[isrc_].x_i_0_;
dx_1_ = f00_jp_smem_[j_ + jstride_ * 0].x_j_1_ - f00_ip_[isrc_].x_i_1_;
dx_2_ = f00_jp_smem_[j_ + jstride_ * 0].x_j_2_ - f00_ip_[isrc_].x_i_2_;
r2 = dx_0_ * dx_0_ + dx_1_ * dx_1_ + dx_2_ * dx_2_ + eps2;
rinv = rsqrt(r2);
mrinv = f00_jp_smem_[j_ + jstride_ * 0].m_j_ * rinv;
mr3inv = mrinv * rinv * rinv;
a_i_0_wcache_ += mr3inv * dx_0_;
a_i_1_wcache_ += mr3inv * dx_1_;
a_i_2_wcache_ += mr3inv * dx_2_;
        // loop 1
dx_0_ = f00_jp_smem_[j_ + jstride_ * 1].x_j_0_ - f00_ip_[isrc_].x_i_0_;
dx_1_ = f00_jp_smem_[j_ + jstride_ * 1].x_j_1_ - f00_ip_[isrc_].x_i_1_;
dx_2_ = f00_jp_smem_[j_ + jstride_ * 1].x_j_2_ - f00_ip_[isrc_].x_i_2_;
r2 = dx_0_ * dx_0_ + dx_1_ * dx_1_ + dx_2_ * dx_2_ + eps2;
rinv = rsqrt(r2);
mrinv = f00_jp_smem_[j_ + jstride_ * 1].m_j_ * rinv;
mr3inv = mrinv * rinv * rinv;
a_i_0_wcache_ += mr3inv * dx_0_;
a_i_1_wcache_ += mr3inv * dx_1_;
a_i_2_wcache_ += mr3inv * dx_2_;
        // loop 2
dx_0_ = f00_jp_smem_[j_ + jstride_ * 2].x_j_0_ - f00_ip_[isrc_].x_i_0_;
dx_1_ = f00_jp_smem_[j_ + jstride_ * 2].x_j_1_ - f00_ip_[isrc_].x_i_1_;
dx_2_ = f00_jp_smem_[j_ + jstride_ * 2].x_j_2_ - f00_ip_[isrc_].x_i_2_;
r2 = dx_0_ * dx_0_ + dx_1_ * dx_1_ + dx_2_ * dx_2_ + eps2;
rinv = rsqrt(r2);
mrinv = f00_jp_smem_[j_ + jstride_ * 2].m_j_ * rinv;
mr3inv = mrinv * rinv * rinv;
a_i_0_wcache_ += mr3inv * dx_0_;
a_i_1_wcache_ += mr3inv * dx_1_;
a_i_2_wcache_ += mr3inv * dx_2_;
        // loop 3
dx_0_ = f00_jp_smem_[j_ + jstride_ * 3].x_j_0_ - f00_ip_[isrc_].x_i_0_;
dx_1_ = f00_jp_smem_[j_ + jstride_ * 3].x_j_1_ - f00_ip_[isrc_].x_i_1_;
dx_2_ = f00_jp_smem_[j_ + jstride_ * 3].x_j_2_ - f00_ip_[isrc_].x_i_2_;
r2 = dx_0_ * dx_0_ + dx_1_ * dx_1_ + dx_2_ * dx_2_ + eps2;
rinv = rsqrt(r2);
mrinv = f00_jp_smem_[j_ + jstride_ * 3].m_j_ * rinv;
mr3inv = mrinv * rinv * rinv;
a_i_0_wcache_ += mr3inv * dx_0_;
a_i_1_wcache_ += mr3inv * dx_1_;
a_i_2_wcache_ += mr3inv * dx_2_;
        // loop 4
dx_0_ = f00_jp_smem_[j_ + jstride_ * 4].x_j_0_ - f00_ip_[isrc_].x_i_0_;
dx_1_ = f00_jp_smem_[j_ + jstride_ * 4].x_j_1_ - f00_ip_[isrc_].x_i_1_;
dx_2_ = f00_jp_smem_[j_ + jstride_ * 4].x_j_2_ - f00_ip_[isrc_].x_i_2_;
r2 = dx_0_ * dx_0_ + dx_1_ * dx_1_ + dx_2_ * dx_2_ + eps2;
rinv = rsqrt(r2);
mrinv = f00_jp_smem_[j_ + jstride_ * 4].m_j_ * rinv;
mr3inv = mrinv * rinv * rinv;
a_i_0_wcache_ += mr3inv * dx_0_;
a_i_1_wcache_ += mr3inv * dx_1_;
a_i_2_wcache_ += mr3inv * dx_2_;
        // loop 5
dx_0_ = f00_jp_smem_[j_ + jstride_ * 5].x_j_0_ - f00_ip_[isrc_].x_i_0_;
dx_1_ = f00_jp_smem_[j_ + jstride_ * 5].x_j_1_ - f00_ip_[isrc_].x_i_1_;
dx_2_ = f00_jp_smem_[j_ + jstride_ * 5].x_j_2_ - f00_ip_[isrc_].x_i_2_;
r2 = dx_0_ * dx_0_ + dx_1_ * dx_1_ + dx_2_ * dx_2_ + eps2;
rinv = rsqrt(r2);
mrinv = f00_jp_smem_[j_ + jstride_ * 5].m_j_ * rinv;
mr3inv = mrinv * rinv * rinv;
a_i_0_wcache_ += mr3inv * dx_0_;
a_i_1_wcache_ += mr3inv * dx_1_;
a_i_2_wcache_ += mr3inv * dx_2_;
        // loop 6
dx_0_ = f00_jp_smem_[j_ + jstride_ * 6].x_j_0_ - f00_ip_[isrc_].x_i_0_;
dx_1_ = f00_jp_smem_[j_ + jstride_ * 6].x_j_1_ - f00_ip_[isrc_].x_i_1_;
dx_2_ = f00_jp_smem_[j_ + jstride_ * 6].x_j_2_ - f00_ip_[isrc_].x_i_2_;
r2 = dx_0_ * dx_0_ + dx_1_ * dx_1_ + dx_2_ * dx_2_ + eps2;
rinv = rsqrt(r2);
mrinv = f00_jp_smem_[j_ + jstride_ * 6].m_j_ * rinv;
mr3inv = mrinv * rinv * rinv;
a_i_0_wcache_ += mr3inv * dx_0_;
a_i_1_wcache_ += mr3inv * dx_1_;
a_i_2_wcache_ += mr3inv * dx_2_;
        // loop 7
dx_0_ = f00_jp_smem_[j_ + jstride_ * 7].x_j_0_ - f00_ip_[isrc_].x_i_0_;
dx_1_ = f00_jp_smem_[j_ + jstride_ * 7].x_j_1_ - f00_ip_[isrc_].x_i_1_;
dx_2_ = f00_jp_smem_[j_ + jstride_ * 7].x_j_2_ - f00_ip_[isrc_].x_i_2_;
r2 = dx_0_ * dx_0_ + dx_1_ * dx_1_ + dx_2_ * dx_2_ + eps2;
rinv = rsqrt(r2);
mrinv = f00_jp_smem_[j_ + jstride_ * 7].m_j_ * rinv;
mr3inv = mrinv * rinv * rinv;
a_i_0_wcache_ += mr3inv * dx_0_;
a_i_1_wcache_ += mr3inv * dx_1_;
a_i_2_wcache_ += mr3inv * dx_2_;

                           }
                       }
                   }
                   __syncthreads();
                   f00_result_smem_[tid_].a_i_0_ = a_i_0_wcache_;f00_result_smem_[tid_].a_i_1_ = a_i_1_wcache_;f00_result_smem_[tid_].a_i_2_ = a_i_2_wcache_;
                   __syncthreads();
                   if (jstride_ > 1) {if (tid_ < kbdim_ / 2) {f00_result_smem_[tid_].a_i_0_ += f00_result_smem_[tid_ + kbdim_ / 2].a_i_0_;f00_result_smem_[tid_].a_i_1_ += f00_result_smem_[tid_ + kbdim_ / 2].a_i_1_;f00_result_smem_[tid_].a_i_2_ += f00_result_smem_[tid_ + kbdim_ / 2].a_i_2_;}// __syncthreads(); // this is not necessary since kbdim_ / 2 <= warp size.
}if (jstride_ > 2) {if (tid_ < kbdim_ / 4) {f00_result_smem_[tid_].a_i_0_ += f00_result_smem_[tid_ + kbdim_ / 4].a_i_0_;f00_result_smem_[tid_].a_i_1_ += f00_result_smem_[tid_ + kbdim_ / 4].a_i_1_;f00_result_smem_[tid_].a_i_2_ += f00_result_smem_[tid_ + kbdim_ / 4].a_i_2_;}// __syncthreads(); // this is not necessary since kbdim_ / 4 <= warp size.
}
                   __syncthreads();

#if 1
                   if (tid_ < nvalidthread_) {
                       int idstoff_ = njdiv_ * kbdim_ * ibid_ + jbid_ + njdiv_ * (tid_ % njhsub_);
                       float4 *srcbuf_ = (float4 *) (f00_result_smem_ + tid_);
                       float4 *dstbuf_ = (float4 *) (f00_result_ + idstoff_);
                       for (int icpy_ = 0; icpy_ < sizeof(f00_result_t) / sizeof(float4); icpy_++) {
                           dstbuf_[icpy_] = srcbuf_[icpy_];
                       }
                   }
#else
                   if (tid_ < nvalidthread_) {
                       f00_result_[idst_] = f00_result_smem_[tid_];
                   }
#endif
               }

/*
 * njdiv_    : # of result fragments per result packet.
 * njdiv_ru_ : njdiv_ rounded up to a power of two.
 * rbdim_    : # of result fragments to be reduced to (rbdim_ / njdiv_) result packets.
 */
               __global__ void
               f00_reducer(int njdiv_, int njdiv_ru_, f00_result_t *f00_result_, f00_result_t *f00_result_sub_)
               {
                   extern __shared__ char smembuf_[];
                   int rbdim_ = blockDim.x;
                   f00_result_t * f00_result_smem_ = (f00_result_t *)smembuf_;
                   f00_result_t * f00_result_smem_packed_ = (f00_result_t *)(smembuf_ + rbdim_ * sizeof(f00_result_t));
                   int tid_ = threadIdx.x;
                   int bid_ = blockIdx.x;
                   int isrc_ = rbdim_ * bid_ + tid_;
                   int ndst_ = rbdim_ / njdiv_;
                   int idst_ = ndst_ * bid_ + tid_;
                   __syncthreads();
                   #if 1
                     {
                         float4 *srcbuf_ = (float4 *)(f00_result_sub_ + rbdim_ * bid_);
                         float4 *dstbuf_ = (float4 *)smembuf_;
                         for (int icpy = 0; icpy < sizeof(f00_result_t)/sizeof(float4); icpy++) {
                             dstbuf_[tid_] = srcbuf_[tid_];
                             dstbuf_ += rbdim_;
                             srcbuf_ += rbdim_;
                         }
                     }
#else
                     f00_result_smem_[tid_] = f00_result_sub_[isrc_];
#endif

                   __syncthreads();

                   int n_ = njdiv_ru_;
                   while (n_ > 1) {
                       n_ /= 2;
                       int ipartner_ = tid_ + n_;
                       if (tid_ % njdiv_ < n_ && ipartner_ % njdiv_ru_ < njdiv_) {
                           f00_result_smem_[tid_].a_i_0_ += f00_result_smem_[ipartner_].a_i_0_;f00_result_smem_[tid_].a_i_1_ += f00_result_smem_[ipartner_].a_i_1_;f00_result_smem_[tid_].a_i_2_ += f00_result_smem_[ipartner_].a_i_2_;
                       }
                       __syncthreads(); // this is not necessary if rbdim_ <= warp size.
                   }
                   __syncthreads();
                   if (tid_ % njdiv_ == 0) {
                       int ipack_ = tid_ / njdiv_;
                       f00_result_smem_packed_[ipack_] = f00_result_smem_[tid_];
                   }
                   __syncthreads();
#if 1
                   {
                       float4 *srcbuf_ = (float4 *) f00_result_smem_packed_;
                       float4 *dstbuf_ = (float4 *) (f00_result_ + ndst_ * bid_);
                       if (tid_ < ndst_) {
                           for (int icpy = 0; icpy < sizeof(f00_result_t) / sizeof(float4); icpy++) {
                               dstbuf_[tid_] = srcbuf_[tid_];
                               dstbuf_ += ndst_;
                               srcbuf_ += ndst_;
                           }
                       }
                   }
#else
                   if (tid_ < ndst_) {
                       f00_result_[idst_] = f00_result_smem_packed_[tid_];
                   }
#endif
               }
