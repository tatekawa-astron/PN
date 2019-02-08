#include <cutil.h>
#include <cutil_inline.h>
#include <gcutil.h>
#include "Limit-sticky9_0.h"
#include "Limit-sticky9_1.h"
#include "Limit-sticky9_2.h"
#include "Limit-sticky9_3.h"
#include "Limit-sticky9_4.h"
#include "Limit-sticky9_5.h"
#include "Limit-sticky9_6.h"




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#define DIM 3


#define REAL double


#define NMAX 270000





void
force(double (*x)[3], double (*v)[3], double *m, double eps,
     double (*a)[3], double *pot, int n)
{
  double r, r2e, r2, reinv, rinv, mrinv, mr3inv, dx[3], a1[3], a1c[3];
  double ac[270000][3];
  double dxb[3], dxc[3], dxbc[3], dvbc[3];
  double r1b2, r1b2e, r1c2, r1c2e, rbc2, rbc2e, r1be, mr1b3e, r1ce, rbce;
  double eps2;
  int i, j, k;

  eps2=eps*eps;



for (i=0;i<n;i++) {
    
for (k=0;k<3;k++) { a[i][k] = 0.0;}

    }
                      /*
                       * dispatcher for a kernel defined in 'Limit-sticky9_0.cu'.
                       *
                       * kbdim_  : # of threads per block. fixed to a multiple of warp size.
                       * nip_    : # of IPs handled in one block.
                       * njdiv_  : # of JP fragments. determined in a heuristic manner as a function of
                       *           ni, nj and device specification such as max grid size and warp size,
                       * nblock_ : # of blocks (ni * njdiv / nip) rounded up to a multiple of nip.
                       * npipe_  : len of each IP array. ni/nip_pack rounded up to a multiple of (nth * nip_pack).
                       * nj_ru_  : nj rounded up to a multiple of (nth * njdiv).
                       */
                      {
#if SHARED_HOSTBUF
                          double * GlobalMem<double>::hostbuf = NULL;
                          int GlobalMem<double>::nbytemax = 0;
#endif
                          // Below may have a room for optimization yet.
                          const int kbdim_ = 64; 
                          const int rbdim_ = 64; 
                          const int nimax_ = kbdim_ * 128; // necessary amount of main mem is larger for larger value.
                          const size_t ksmemsize_ = kbdim_ * sizeof(f00_jp_t);
                          const size_t rsmemsize_ = rbdim_ * sizeof(f00_result_t) * 2;
                          const int nip_ = kbdim_ * 1;
                          const GoosePrecision gprec_ = kGoosePrecisionDouble;
                          static int njdiv_, njdiv_ru_;
                          static int njold = 0, niold = 0;
                          static int nblockmax_, nsp_;
                          static size_t smemsizemax_;
                          static int firstcall_ = 1;
                          int bufidx_, ioff_, nisub_;
                          int ni_ru_ = ((n - 0 - 1) / nimax_ + 1) * nimax_;
                          cudaError cuerr_;
                          if (firstcall_) {
                              int device_;
                              cudaDeviceProp prop;
                              firstcall_ = 0;
                              cuerr_ = cudaGetDevice(&device_);
                              cutilSafeCall(cuerr_);
                              cuerr_ = cudaGetDeviceProperties(&prop, device_);
                              cutilSafeCall(cuerr_);
                              if (gprec_ >= kGoosePrecisionDouble) {
                                  if (prop.major < 2 && prop.minor < 3) {
                                      fprintf(stderr, "GPU architecture sm_%d%d does not support double-precision arithmetic.\n",
                                              prop.major, prop.minor);
                                      exit(1);
                                  }
                              }
                              smemsizemax_ = prop.sharedMemPerBlock;
                              if (ksmemsize_ >= smemsizemax_) {
                                  fprintf(stderr,
                                          "Shared memory consumption in f00_calculator (=%dbyte) "
                                          "reached the limit (=%dbyte). Too many jvars or fvars used.\n",
                                          ksmemsize_, smemsizemax_);
                                  exit(1);
                              }
                              if (rsmemsize_ >= smemsizemax_) {
                                  fprintf(stderr,
                                          "Shared memory consumption in f00_reducer (=%dbyte) "
                                          "reached the limit (=%dbyte). Too many fvars used.\n",
                                          rsmemsize_, smemsizemax_);
                                  exit(1);
                              }
                              nblockmax_ = prop.maxGridSize[0];
                              nsp_ = prop.multiProcessorCount * 8; // # of spream processors.
                              // fprintf(stderr, "nblockmax:%d  nsp:%d\n", nblockmax_, nsp_);
                          }
                          nisub_ = nimax_;
                          for (ioff_ = 0; ioff_ < n; ioff_ += nisub_) {
                              if (ioff_ + nisub_ > n) {
                                  nisub_ = n - ioff_;
                              }

                              int nisub_ru_ = ((nisub_ - 1) / nip_ + 1) * nip_;
                              if (njold != n || niold != nisub_) {
                                  // Adjust # of JP fragments so that large enough # of threads (to fill all SPs, 
                                  // and to hide latency of the global memory) are dispatched, at long as each
                                  // JP fragments has several times warp size.
                                  // You may want to hand tune the following part to obtain the optimal 'njdiv_'
                                  // value for a given hardware configuration.
                                  njdiv_ = njdiv_ru_ = 1;
                                  while (nisub_ru_ / nip_ * njdiv_ * 2 <= nblockmax_ &&
                                         nisub_ * njdiv_ < nsp_ * 200 &&
                                         n / njdiv_ > kbdim_ * 6 &&
                                         njdiv_ * 2 <= rbdim_) {
                                      njdiv_ *= 2;
                                  }
                                  // njdiv_ru_ is set to njdiv_ rounded up to a power of two.
                                  while (njdiv_ru_ < njdiv_) {
                                      njdiv_ru_ *= 2;
                                  }
                                  njold = n;
                                  niold = nisub_;
#if 0
                                  fprintf(stderr, "\n");
                                  fprintf(stderr, "nj:%d  njdiv:%d  nj/njdiv:%d  ni:%d  ni * njdiv:%d\n",
                                          n, njdiv_, n/njdiv_, nisub_, nisub_ * njdiv_);
                                  fprintf(stderr, "\n");
#endif
                              }
                              dim3 kthreads_(kbdim_, 1, 1);
                              dim3 kgrids_(njdiv_, nisub_ru_ / nip_, 1);
                              dim3 rthreads_(rbdim_, 1, 1);
                              int npipe_ = nisub_ru_ / 1;
                              int nrblock_ = npipe_ * njdiv_ / rbdim_;
                              dim3 rgrids_(nrblock_,1, 1);
                              int nj_ru_ = ((n - 0 - 1) / (kbdim_ * njdiv_) + 1) * (kbdim_ * njdiv_);
                              static int jbufsize_ = 0, ibufsize_ = 0, rsubbufsize_ = 0, rbufsize_ = 0;
                              static GlobalMem<f00_jp_t> f00_jp_;
                              static GlobalMem<f00_ip_t> f00_ip_;
                              static GlobalMem<f00_result_t> f00_result_;
                              static GlobalMem<f00_result_t> f00_result_sub_;
                              if (nj_ru_ > jbufsize_) {
                                  jbufsize_ = nj_ru_;
                                  f00_jp_.realloc(jbufsize_);
                              }

                              // Here we need to alloc n IPs. Note that npipe_ IPs would not be enough,
                              // since IPs on the device memory should not be overwritten during the ioff_ loop.
                              if (ni_ru_ > ibufsize_) {
                                  ibufsize_ = ni_ru_;
                                  f00_ip_.realloc(ibufsize_);
                              }
                              if (ni_ru_ > rbufsize_) {
                                  rbufsize_ = ni_ru_;
                                  f00_result_.realloc(rbufsize_);
                              }
                              if (npipe_ * njdiv_ > rsubbufsize_) {
                                  rsubbufsize_ = npipe_ * njdiv_ ;
                                  f00_result_sub_.realloc(rsubbufsize_);
                              }

                              for (j = 0, bufidx_ = 0 ; j <n; j++, bufidx_++) {f00_jp_[bufidx_].x_j_0_ = x[j][0];f00_jp_[bufidx_].x_j_1_ = x[j][1];f00_jp_[bufidx_].x_j_2_ = x[j][2];f00_jp_[bufidx_].m_j_ = m[j];}f00_jp_.htod(0, bufidx_);

                              if (ioff_ < n) { // otherwise nothing to send.
                                  int nipsend_ = npipe_;
                                  if (ioff_ + nipsend_ >= n) {
                                      nipsend_ = n - ioff_;
                                  }
                                  for (i = 0 ; i < nipsend_; i++) {f00_ip_[ioff_ + i].x_i_0_ = x[ioff_ + i][0];f00_ip_[ioff_ + i].x_i_1_ = x[ioff_ + i][1];f00_ip_[ioff_ + i].x_i_2_ = x[ioff_ + i][2];}f00_ip_.htod(ioff_, nipsend_);

                              }
                              f00_calculator<<<kgrids_, kthreads_, ksmemsize_>>>(ioff_, nisub_, n - 0, (double)eps2, f00_jp_, (f00_ip_ + ioff_), f00_result_sub_);
                              f00_reducer<<<rgrids_, rthreads_, rsmemsize_>>>(njdiv_, njdiv_ru_, (f00_result_ + ioff_), f00_result_sub_);

                              if (ioff_ < n) { // otherwise nothing to receive.
                                  int niprecv_ = npipe_;
                                  if (ioff_ + niprecv_ >= n) {
                                      niprecv_ = n - ioff_;
                                  }
                                  f00_result_.dtoh(ioff_, niprecv_);
for (i = 0 ; i < niprecv_; i++) {a[ioff_ + i][0] = f00_result_[ioff_ + i].a_i_0_;a[ioff_ + i][1] = f00_result_[ioff_ + i].a_i_1_;a[ioff_ + i][2] = f00_result_[ioff_ + i].a_i_2_;}
                              }
                          } // end of ioff loop.
                      } // end of api calls.



  for(k=0;k<3;k++) a1[k] = 0.0;



  for(i=1;i<n;i++) {
    for(k=0;k<3;k++) {
      dx[k] = x[i][k] - x[0][k];
    }
    r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    rinv = rsqrt(r2);
    r2e = r2 + eps2;
    reinv = rsqrt(r2e);
    mrinv = m[i]*reinv;
    mr3inv = mrinv*reinv*reinv;

    a1[0] += mr3inv * dx[0]
      * (5.0*m[0]*reinv
  -2.0*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2])
  +1.5*(v[i][0]*dx[0]+v[i][1]*dx[1]+v[i][2]*dx[2])
  *(v[i][0]*dx[0]+v[i][1]*dx[1]+v[i][2]*dx[2])
  /r2)
      +3.0*mr3inv*v[i][0]
      * (v[i][0]*dx[0]+v[i][1]*dx[1]+v[i][2]*dx[2]);

    a1[1] += mr3inv * dx[1]
      * (5.0*m[0]*reinv
  -2.0*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2])
  +1.5*(v[i][0]*dx[0]+v[i][1]*dx[1]+v[i][2]*dx[2])
  *(v[i][0]*dx[0]+v[i][1]*dx[1]+v[i][2]*dx[2])
  /r2)
      +3.0*mr3inv*v[i][1]
      * (v[i][0]*dx[0]+v[i][1]*dx[1]+v[i][2]*dx[2]);

    a1[2] += mr3inv * dx[2]
      * (5.0*m[0]*reinv
  -2.0*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2])
  +1.5*(v[i][0]*dx[0]+v[i][1]*dx[1]+v[i][2]*dx[2])
  *(v[i][0]*dx[0]+v[i][1]*dx[1]+v[i][2]*dx[2])
  /r2)
      +3.0*mr3inv*v[i][2]
      * (v[i][0]*dx[0]+v[i][1]*dx[1]+v[i][2]*dx[2]);
  }

  for(k=0;k<3;k++) a[0][k]+=a1[k];



  for(k=0;k<3;k++) a1c[k]=0.0;
  for(i=1;i<n;i++) {
    for(k=0;k<3;k++) {
   dxb[k]=x[i][k]-x[0][k];
    }
      r1b2=dxb[0]*dxb[0]+dxb[1]*dxb[1]+dxb[2]*dxb[2];
      r1b2e=r1b2+eps2;
      r1be=rsqrt(r1b2e);
      mr1b3e=m[i]*r1be*r1be*r1be;

    a1c[0]+=4.0*mr1b3e*r1be*m[i]*dxb[0];
    a1c[1]+=4.0*mr1b3e*r1be*m[i]*dxb[1];
    a1c[2]+=4.0*mr1b3e*r1be*m[i]*dxb[2];

  }
                      /*
                       * dispatcher for a kernel defined in 'Limit-sticky9_1.cu'.
                       *
                       * kbdim_  : # of threads per block. fixed to a multiple of warp size.
                       * nip_    : # of IPs handled in one block.
                       * njdiv_  : # of JP fragments. determined in a heuristic manner as a function of
                       *           ni, nj and device specification such as max grid size and warp size,
                       * nblock_ : # of blocks (ni * njdiv / nip) rounded up to a multiple of nip.
                       * npipe_  : len of each IP array. ni/nip_pack rounded up to a multiple of (nth * nip_pack).
                       * nj_ru_  : nj rounded up to a multiple of (nth * njdiv).
                       */
                      {
#if SHARED_HOSTBUF
                          double * GlobalMem<double>::hostbuf = NULL;
                          int GlobalMem<double>::nbytemax = 0;
#endif
                          // Below may have a room for optimization yet.
                          const int kbdim_ = 64; 
                          const int rbdim_ = 64; 
                          const int nimax_ = kbdim_ * 128; // necessary amount of main mem is larger for larger value.
                          const size_t ksmemsize_ = kbdim_ * sizeof(f01_jp_t);
                          const size_t rsmemsize_ = rbdim_ * sizeof(f01_result_t) * 2;
                          const int nip_ = kbdim_ * 1;
                          const GoosePrecision gprec_ = kGoosePrecisionDouble;
                          static int njdiv_, njdiv_ru_;
                          static int njold = 0, niold = 0;
                          static int nblockmax_, nsp_;
                          static size_t smemsizemax_;
                          static int firstcall_ = 1;
                          int bufidx_, ioff_, nisub_;
                          int ni_ru_ = ((n - 1 - 1) / nimax_ + 1) * nimax_;
                          cudaError cuerr_;
                          if (firstcall_) {
                              int device_;
                              cudaDeviceProp prop;
                              firstcall_ = 0;
                              cuerr_ = cudaGetDevice(&device_);
                              cutilSafeCall(cuerr_);
                              cuerr_ = cudaGetDeviceProperties(&prop, device_);
                              cutilSafeCall(cuerr_);
                              if (gprec_ >= kGoosePrecisionDouble) {
                                  if (prop.major < 2 && prop.minor < 3) {
                                      fprintf(stderr, "GPU architecture sm_%d%d does not support double-precision arithmetic.\n",
                                              prop.major, prop.minor);
                                      exit(1);
                                  }
                              }
                              smemsizemax_ = prop.sharedMemPerBlock;
                              if (ksmemsize_ >= smemsizemax_) {
                                  fprintf(stderr,
                                          "Shared memory consumption in f01_calculator (=%dbyte) "
                                          "reached the limit (=%dbyte). Too many jvars or fvars used.\n",
                                          ksmemsize_, smemsizemax_);
                                  exit(1);
                              }
                              if (rsmemsize_ >= smemsizemax_) {
                                  fprintf(stderr,
                                          "Shared memory consumption in f01_reducer (=%dbyte) "
                                          "reached the limit (=%dbyte). Too many fvars used.\n",
                                          rsmemsize_, smemsizemax_);
                                  exit(1);
                              }
                              nblockmax_ = prop.maxGridSize[0];
                              nsp_ = prop.multiProcessorCount * 8; // # of spream processors.
                              // fprintf(stderr, "nblockmax:%d  nsp:%d\n", nblockmax_, nsp_);
                          }
                          nisub_ = nimax_;
                          for (ioff_ = 1; ioff_ < n; ioff_ += nisub_) {
                              if (ioff_ + nisub_ > n) {
                                  nisub_ = n - ioff_;
                              }

                              int nisub_ru_ = ((nisub_ - 1) / nip_ + 1) * nip_;
                              if (njold != n || niold != nisub_) {
                                  // Adjust # of JP fragments so that large enough # of threads (to fill all SPs, 
                                  // and to hide latency of the global memory) are dispatched, at long as each
                                  // JP fragments has several times warp size.
                                  // You may want to hand tune the following part to obtain the optimal 'njdiv_'
                                  // value for a given hardware configuration.
                                  njdiv_ = njdiv_ru_ = 1;
                                  while (nisub_ru_ / nip_ * njdiv_ * 2 <= nblockmax_ &&
                                         nisub_ * njdiv_ < nsp_ * 200 &&
                                         n / njdiv_ > kbdim_ * 6 &&
                                         njdiv_ * 2 <= rbdim_) {
                                      njdiv_ *= 2;
                                  }
                                  // njdiv_ru_ is set to njdiv_ rounded up to a power of two.
                                  while (njdiv_ru_ < njdiv_) {
                                      njdiv_ru_ *= 2;
                                  }
                                  njold = n;
                                  niold = nisub_;
#if 0
                                  fprintf(stderr, "\n");
                                  fprintf(stderr, "nj:%d  njdiv:%d  nj/njdiv:%d  ni:%d  ni * njdiv:%d\n",
                                          n, njdiv_, n/njdiv_, nisub_, nisub_ * njdiv_);
                                  fprintf(stderr, "\n");
#endif
                              }
                              dim3 kthreads_(kbdim_, 1, 1);
                              dim3 kgrids_(njdiv_, nisub_ru_ / nip_, 1);
                              dim3 rthreads_(rbdim_, 1, 1);
                              int npipe_ = nisub_ru_ / 1;
                              int nrblock_ = npipe_ * njdiv_ / rbdim_;
                              dim3 rgrids_(nrblock_,1, 1);
                              int nj_ru_ = ((n - 1 - 1) / (kbdim_ * njdiv_) + 1) * (kbdim_ * njdiv_);
                              static int jbufsize_ = 0, ibufsize_ = 0, rsubbufsize_ = 0, rbufsize_ = 0;
                              static GlobalMem<f01_jp_t> f01_jp_;
                              static GlobalMem<f01_ip_t> f01_ip_;
                              static GlobalMem<f01_result_t> f01_result_;
                              static GlobalMem<f01_result_t> f01_result_sub_;
                              if (nj_ru_ > jbufsize_) {
                                  jbufsize_ = nj_ru_;
                                  f01_jp_.realloc(jbufsize_);
                              }

                              // Here we need to alloc n IPs. Note that npipe_ IPs would not be enough,
                              // since IPs on the device memory should not be overwritten during the ioff_ loop.
                              if (ni_ru_ > ibufsize_) {
                                  ibufsize_ = ni_ru_;
                                  f01_ip_.realloc(ibufsize_);
                              }
                              if (ni_ru_ > rbufsize_) {
                                  rbufsize_ = ni_ru_;
                                  f01_result_.realloc(rbufsize_);
                              }
                              if (npipe_ * njdiv_ > rsubbufsize_) {
                                  rsubbufsize_ = npipe_ * njdiv_ ;
                                  f01_result_sub_.realloc(rsubbufsize_);
                              }

                              for (j = 1, bufidx_ = 0 ; j <n; j++, bufidx_++) {f01_jp_[bufidx_].x_j_0_ = x[j][0];f01_jp_[bufidx_].x_j_1_ = x[j][1];f01_jp_[bufidx_].x_j_2_ = x[j][2];f01_jp_[bufidx_].m_j_ = m[j];f01_jp_[bufidx_].v_j_0_ = v[j][0];f01_jp_[bufidx_].v_j_1_ = v[j][1];f01_jp_[bufidx_].v_j_2_ = v[j][2];}f01_jp_.htod(0, bufidx_);

                              if (ioff_ < n) { // otherwise nothing to send.
                                  int nipsend_ = npipe_;
                                  if (ioff_ + nipsend_ >= n) {
                                      nipsend_ = n - ioff_;
                                  }
                                  for (i = 0 ; i < nipsend_; i++) {f01_ip_[ioff_ + i].x_i_0_ = x[ioff_ + i][0];f01_ip_[ioff_ + i].x_i_1_ = x[ioff_ + i][1];f01_ip_[ioff_ + i].x_i_2_ = x[ioff_ + i][2];f01_ip_[ioff_ + i].m_i_ = m[ioff_ + i];f01_ip_[ioff_ + i].v_i_0_ = v[ioff_ + i][0];f01_ip_[ioff_ + i].v_i_1_ = v[ioff_ + i][1];f01_ip_[ioff_ + i].v_i_2_ = v[ioff_ + i][2];}f01_ip_.htod(ioff_, nipsend_);

                              }
                              f01_calculator<<<kgrids_, kthreads_, ksmemsize_>>>(ioff_, nisub_, n - 1, (double)x[0][0], (double)x[0][1], (double)x[0][2], (double)eps2, (double)m[0], f01_jp_, (f01_ip_ + ioff_), f01_result_sub_);
                              f01_reducer<<<rgrids_, rthreads_, rsmemsize_>>>(njdiv_, njdiv_ru_, (f01_result_ + ioff_), f01_result_sub_);

                              if (ioff_ < n) { // otherwise nothing to receive.
                                  int niprecv_ = npipe_;
                                  if (ioff_ + niprecv_ >= n) {
                                      niprecv_ = n - ioff_;
                                  }
                                  f01_result_.dtoh(ioff_, niprecv_);
for (i = 0 ; i < niprecv_; i++) {a1c[0] += f01_result_[ioff_ + i].a1c_0_;a1c[1] += f01_result_[ioff_ + i].a1c_1_;a1c[2] += f01_result_[ioff_ + i].a1c_2_;}
                              }
                          } // end of ioff loop.
                      } // end of api calls.





  for(i=1;i<n;i++) {
    for(k=0;k<3;k++) {
      dxb[k]=x[0][k]-x[i][k];
    }

    r1b2=dxb[0]*dxb[0]+dxb[1]*dxb[1]+dxb[2]*dxb[2];
    r1b2e=r1b2+eps2;
    r1be=rsqrt(r1b2e);

    mr1b3e=m[i]*r1be*r1be*r1be;

    a1c[0]-=
      mr1b3e*m[i]*dxb[0]* (4.0*r1be+1.25*eps)
      -mr1b3e*m[i]/m[0]
      *(4.0*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]) *dxb[0]
 -7.0*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])*v[i][0]);

    a1c[1]-=
      mr1b3e*m[i]*dxb[1]* (4.0*r1be+1.25*eps)
      -mr1b3e*m[i]/m[0]
      *(4.0*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]) *dxb[1]
     -7.0*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])*v[i][1]);

    a1c[2]-=
      mr1b3e*m[i]*dxb[2]* (4.0*r1be+1.25*eps)
      -mr1b3e*m[i]/m[0]
      *(4.0*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]) *dxb[2]
 -7.0*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])*v[i][2]);
  }

  for(k=0;k<3;k++) {
    a[0][k]+=a1c[k];
  }



  for(i=1;i<n;i++) {
    for(k=0;k<3;k++) {
      dxb[k]=x[0][k]-x[i][k];
    }

    r1b2=dxb[0]*dxb[0]+dxb[1]*dxb[1]+dxb[2]*dxb[2];
    r1b2e=r1b2+eps2;
    r1be=rsqrt(r1b2e);

    mrinv=m[0]*r1be;
    mr3inv=mrinv*r1be*r1be;

    a[i][0]+=mr3inv*dxb[0]
      * (4.0*mrinv-(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]))
      + 4.0*mr3inv*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])*v[i][0];

    a[i][1]+=mr3inv*dxb[1]
      * (4.0*mrinv-(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]))
      + 4.0*mr3inv*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])*v[i][1];

    a[i][2]+=mr3inv*dxb[2]
      * (4.0*mrinv-(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]))
      + 4.0*mr3inv*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])*v[i][2];
  }



  for(i=1;i<n;i++){
    for(k=0;k<3;k++){
      ac[i][k]=0.0;
    }
  }


  for(i=1;i<n;i++){
    for(k=0;k<3;k++) {
      dxb[k]=x[i][k]-x[0][k];
    }

    r1b2=dxb[0]*dxb[0]+dxb[1]*dxb[1]+dxb[2]*dxb[2];
    r1b2e=r1b2+eps2;
    r1be=rsqrt(r1b2e);

    mr3inv=m[0]*r1be*r1be*r1be;

    a[i][0]+=5.0*mr3inv*m[i]*r1be*dxb[0]
      -m[i]*r1be*r1be*r1be
      *(4.0*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2])*dxb[0]
       -7.0*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])*v[i][0]);

    a[i][1]+=5.0*mr3inv*m[i]*r1be*dxb[1]
      -m[i]*r1be*r1be*r1be
      *(4.0*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2])*dxb[1]
       -7.0*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])*v[i][1]);

    a[i][2]+=5.0*mr3inv*m[i]*r1be*dxb[2]
      -m[i]*r1be*r1be*r1be
      *(4.0*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2])*dxb[2]
       -7.0*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])*v[i][2]);
  }

                      /*
                       * dispatcher for a kernel defined in 'Limit-sticky9_2.cu'.
                       *
                       * kbdim_  : # of threads per block. fixed to a multiple of warp size.
                       * nip_    : # of IPs handled in one block.
                       * njdiv_  : # of JP fragments. determined in a heuristic manner as a function of
                       *           ni, nj and device specification such as max grid size and warp size,
                       * nblock_ : # of blocks (ni * njdiv / nip) rounded up to a multiple of nip.
                       * npipe_  : len of each IP array. ni/nip_pack rounded up to a multiple of (nth * nip_pack).
                       * nj_ru_  : nj rounded up to a multiple of (nth * njdiv).
                       */
                      {
#if SHARED_HOSTBUF
                          double * GlobalMem<double>::hostbuf = NULL;
                          int GlobalMem<double>::nbytemax = 0;
#endif
                          // Below may have a room for optimization yet.
                          const int kbdim_ = 64; 
                          const int rbdim_ = 64; 
                          const int nimax_ = kbdim_ * 128; // necessary amount of main mem is larger for larger value.
                          const size_t ksmemsize_ = kbdim_ * sizeof(f02_jp_t);
                          const size_t rsmemsize_ = rbdim_ * sizeof(f02_result_t) * 2;
                          const int nip_ = kbdim_ * 1;
                          const GoosePrecision gprec_ = kGoosePrecisionDouble;
                          static int njdiv_, njdiv_ru_;
                          static int njold = 0, niold = 0;
                          static int nblockmax_, nsp_;
                          static size_t smemsizemax_;
                          static int firstcall_ = 1;
                          int bufidx_, ioff_, nisub_;
                          int ni_ru_ = ((n - 1 - 1) / nimax_ + 1) * nimax_;
                          cudaError cuerr_;
                          if (firstcall_) {
                              int device_;
                              cudaDeviceProp prop;
                              firstcall_ = 0;
                              cuerr_ = cudaGetDevice(&device_);
                              cutilSafeCall(cuerr_);
                              cuerr_ = cudaGetDeviceProperties(&prop, device_);
                              cutilSafeCall(cuerr_);
                              if (gprec_ >= kGoosePrecisionDouble) {
                                  if (prop.major < 2 && prop.minor < 3) {
                                      fprintf(stderr, "GPU architecture sm_%d%d does not support double-precision arithmetic.\n",
                                              prop.major, prop.minor);
                                      exit(1);
                                  }
                              }
                              smemsizemax_ = prop.sharedMemPerBlock;
                              if (ksmemsize_ >= smemsizemax_) {
                                  fprintf(stderr,
                                          "Shared memory consumption in f02_calculator (=%dbyte) "
                                          "reached the limit (=%dbyte). Too many jvars or fvars used.\n",
                                          ksmemsize_, smemsizemax_);
                                  exit(1);
                              }
                              if (rsmemsize_ >= smemsizemax_) {
                                  fprintf(stderr,
                                          "Shared memory consumption in f02_reducer (=%dbyte) "
                                          "reached the limit (=%dbyte). Too many fvars used.\n",
                                          rsmemsize_, smemsizemax_);
                                  exit(1);
                              }
                              nblockmax_ = prop.maxGridSize[0];
                              nsp_ = prop.multiProcessorCount * 8; // # of spream processors.
                              // fprintf(stderr, "nblockmax:%d  nsp:%d\n", nblockmax_, nsp_);
                          }
                          nisub_ = nimax_;
                          for (ioff_ = 1; ioff_ < n; ioff_ += nisub_) {
                              if (ioff_ + nisub_ > n) {
                                  nisub_ = n - ioff_;
                              }

                              int nisub_ru_ = ((nisub_ - 1) / nip_ + 1) * nip_;
                              if (njold != n || niold != nisub_) {
                                  // Adjust # of JP fragments so that large enough # of threads (to fill all SPs, 
                                  // and to hide latency of the global memory) are dispatched, at long as each
                                  // JP fragments has several times warp size.
                                  // You may want to hand tune the following part to obtain the optimal 'njdiv_'
                                  // value for a given hardware configuration.
                                  njdiv_ = njdiv_ru_ = 1;
                                  while (nisub_ru_ / nip_ * njdiv_ * 2 <= nblockmax_ &&
                                         nisub_ * njdiv_ < nsp_ * 200 &&
                                         n / njdiv_ > kbdim_ * 6 &&
                                         njdiv_ * 2 <= rbdim_) {
                                      njdiv_ *= 2;
                                  }
                                  // njdiv_ru_ is set to njdiv_ rounded up to a power of two.
                                  while (njdiv_ru_ < njdiv_) {
                                      njdiv_ru_ *= 2;
                                  }
                                  njold = n;
                                  niold = nisub_;
#if 0
                                  fprintf(stderr, "\n");
                                  fprintf(stderr, "nj:%d  njdiv:%d  nj/njdiv:%d  ni:%d  ni * njdiv:%d\n",
                                          n, njdiv_, n/njdiv_, nisub_, nisub_ * njdiv_);
                                  fprintf(stderr, "\n");
#endif
                              }
                              dim3 kthreads_(kbdim_, 1, 1);
                              dim3 kgrids_(njdiv_, nisub_ru_ / nip_, 1);
                              dim3 rthreads_(rbdim_, 1, 1);
                              int npipe_ = nisub_ru_ / 1;
                              int nrblock_ = npipe_ * njdiv_ / rbdim_;
                              dim3 rgrids_(nrblock_,1, 1);
                              int nj_ru_ = ((n - 1 - 1) / (kbdim_ * njdiv_) + 1) * (kbdim_ * njdiv_);
                              static int jbufsize_ = 0, ibufsize_ = 0, rsubbufsize_ = 0, rbufsize_ = 0;
                              static GlobalMem<f02_jp_t> f02_jp_;
                              static GlobalMem<f02_ip_t> f02_ip_;
                              static GlobalMem<f02_result_t> f02_result_;
                              static GlobalMem<f02_result_t> f02_result_sub_;
                              if (nj_ru_ > jbufsize_) {
                                  jbufsize_ = nj_ru_;
                                  f02_jp_.realloc(jbufsize_);
                              }

                              // Here we need to alloc n IPs. Note that npipe_ IPs would not be enough,
                              // since IPs on the device memory should not be overwritten during the ioff_ loop.
                              if (ni_ru_ > ibufsize_) {
                                  ibufsize_ = ni_ru_;
                                  f02_ip_.realloc(ibufsize_);
                              }
                              if (ni_ru_ > rbufsize_) {
                                  rbufsize_ = ni_ru_;
                                  f02_result_.realloc(rbufsize_);
                              }
                              if (npipe_ * njdiv_ > rsubbufsize_) {
                                  rsubbufsize_ = npipe_ * njdiv_ ;
                                  f02_result_sub_.realloc(rsubbufsize_);
                              }

                              for (j = 1, bufidx_ = 0 ; j <n; j++, bufidx_++) {f02_jp_[bufidx_].x_j_0_ = x[j][0];f02_jp_[bufidx_].v_j_0_ = v[j][0];f02_jp_[bufidx_].x_j_1_ = x[j][1];f02_jp_[bufidx_].v_j_1_ = v[j][1];f02_jp_[bufidx_].x_j_2_ = x[j][2];f02_jp_[bufidx_].v_j_2_ = v[j][2];f02_jp_[bufidx_].m_j_ = m[j];}f02_jp_.htod(0, bufidx_);

                              if (ioff_ < n) { // otherwise nothing to send.
                                  int nipsend_ = npipe_;
                                  if (ioff_ + nipsend_ >= n) {
                                      nipsend_ = n - ioff_;
                                  }
                                  for (i = 0 ; i < nipsend_; i++) {f02_ip_[ioff_ + i].x_i_0_ = x[ioff_ + i][0];f02_ip_[ioff_ + i].v_i_0_ = v[ioff_ + i][0];f02_ip_[ioff_ + i].x_i_1_ = x[ioff_ + i][1];f02_ip_[ioff_ + i].v_i_1_ = v[ioff_ + i][1];f02_ip_[ioff_ + i].x_i_2_ = x[ioff_ + i][2];f02_ip_[ioff_ + i].v_i_2_ = v[ioff_ + i][2];}f02_ip_.htod(ioff_, nipsend_);

                              }
                              f02_calculator<<<kgrids_, kthreads_, ksmemsize_>>>(ioff_, nisub_, n - 1, (double)x[0][0], (double)x[0][1], (double)x[0][2], (double)eps2, (double)m[0], f02_jp_, (f02_ip_ + ioff_), f02_result_sub_);
                              f02_reducer<<<rgrids_, rthreads_, rsmemsize_>>>(njdiv_, njdiv_ru_, (f02_result_ + ioff_), f02_result_sub_);

                              if (ioff_ < n) { // otherwise nothing to receive.
                                  int niprecv_ = npipe_;
                                  if (ioff_ + niprecv_ >= n) {
                                      niprecv_ = n - ioff_;
                                  }
                                  f02_result_.dtoh(ioff_, niprecv_);
for (i = 0 ; i < niprecv_; i++) {ac[ioff_ + i][0] = f02_result_[ioff_ + i].ac_i_0_;ac[ioff_ + i][1] = f02_result_[ioff_ + i].ac_i_1_;ac[ioff_ + i][2] = f02_result_[ioff_ + i].ac_i_2_;}
                              }
                          } // end of ioff loop.
                      } // end of api calls.





  for(i=1;i<n;i++){
    for(k=0;k<3;k++){
      a[i][k]+=ac[i][k];
    }
  }

  for(i=1;i<n;i++){
    for(k=0;k<3;k++) {
      dxb[k]=x[0][k]-x[i][k];
    }

    r1b2=dxb[0]*dxb[0]+dxb[1]*dxb[1]+dxb[2]*dxb[2];
    r1b2e=r1b2+eps2;
    r1be=rsqrt(r1b2e);
    mr3inv=m[i]*r1be*r1be*r1be;

    a[i][0]-=mr3inv*m[0]*dxb[0]
      *(4.0/eps+1.25*r1be+0.25*r1b2*r1be*r1be*r1be)
      -3.5*mr3inv*m[0]*(1.0/eps-r1be)*dxb[0]
      -mr3inv
      *(4.0*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2])*dxb[0]
       -7.0*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])*v[i][0]);

    a[i][1]-=mr3inv*m[0]*dxb[1]
      *(4.0/eps+1.25*r1be+0.25*r1b2*r1be*r1be*r1be)
      -3.5*mr3inv*m[0]*(1.0/eps-r1be)*dxb[1]
      -mr3inv
      *(4.0*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2])*dxb[1]
       -7.0*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])*v[i][1]);

    a[i][2]-=mr3inv*m[0]*dxb[2]
      *(4.0/eps+1.25*r1be+0.25*r1b2*r1be*r1be*r1be)
      -3.5*mr3inv*m[0]*(1.0/eps-r1be)*dxb[2]
      -mr3inv
      *(4.0*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2])*dxb[2]
       -7.0*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])*v[i][2]);
  }

}


void energy(double t,
            double x[270000][3],
            double v[270000][3],
            double m[270000],
            int n,
            double init_ene,
            double eps)
{
  double pot[270000], DE, ene, eps2;
  double dx[3], dxa[3], dxb[3], dxab[3];
  double r2, rinv, mrinv, vi2;
  double r1a2e, r1ainv, mr1ainv, r1b2e, r1binv, mr1binv;
  double rab2e, rabinv;
  double kin_n, pot_n;
  double pot_pn, pot_pn2[270000];
  int i,j,k;
  double cm[3];

  eps2=eps*eps;
  kin_n=0.0;
  pot_n=0.0;
  pot_pn=0.0;



  for (i=0;i<n;i++) {
    kin_n+=m[i]*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]);
  }

  kin_n*=0.5;

  for (i=0;i<n;i++) {
    pot[i]=0.0;
  }

                      /*
                       * dispatcher for a kernel defined in 'Limit-sticky9_3.cu'.
                       *
                       * kbdim_  : # of threads per block. fixed to a multiple of warp size.
                       * nip_    : # of IPs handled in one block.
                       * njdiv_  : # of JP fragments. determined in a heuristic manner as a function of
                       *           ni, nj and device specification such as max grid size and warp size,
                       * nblock_ : # of blocks (ni * njdiv / nip) rounded up to a multiple of nip.
                       * npipe_  : len of each IP array. ni/nip_pack rounded up to a multiple of (nth * nip_pack).
                       * nj_ru_  : nj rounded up to a multiple of (nth * njdiv).
                       */
                      {
#if SHARED_HOSTBUF
                          double * GlobalMem<double>::hostbuf = NULL;
                          int GlobalMem<double>::nbytemax = 0;
#endif
                          // Below may have a room for optimization yet.
                          const int kbdim_ = 64; 
                          const int rbdim_ = 64; 
                          const int nimax_ = kbdim_ * 128; // necessary amount of main mem is larger for larger value.
                          const size_t ksmemsize_ = kbdim_ * sizeof(f03_jp_t);
                          const size_t rsmemsize_ = rbdim_ * sizeof(f03_result_t) * 2;
                          const int nip_ = kbdim_ * 1;
                          const GoosePrecision gprec_ = kGoosePrecisionDouble;
                          static int njdiv_, njdiv_ru_;
                          static int njold = 0, niold = 0;
                          static int nblockmax_, nsp_;
                          static size_t smemsizemax_;
                          static int firstcall_ = 1;
                          int bufidx_, ioff_, nisub_;
                          int ni_ru_ = ((n - 0 - 1) / nimax_ + 1) * nimax_;
                          cudaError cuerr_;
                          if (firstcall_) {
                              int device_;
                              cudaDeviceProp prop;
                              firstcall_ = 0;
                              cuerr_ = cudaGetDevice(&device_);
                              cutilSafeCall(cuerr_);
                              cuerr_ = cudaGetDeviceProperties(&prop, device_);
                              cutilSafeCall(cuerr_);
                              if (gprec_ >= kGoosePrecisionDouble) {
                                  if (prop.major < 2 && prop.minor < 3) {
                                      fprintf(stderr, "GPU architecture sm_%d%d does not support double-precision arithmetic.\n",
                                              prop.major, prop.minor);
                                      exit(1);
                                  }
                              }
                              smemsizemax_ = prop.sharedMemPerBlock;
                              if (ksmemsize_ >= smemsizemax_) {
                                  fprintf(stderr,
                                          "Shared memory consumption in f03_calculator (=%dbyte) "
                                          "reached the limit (=%dbyte). Too many jvars or fvars used.\n",
                                          ksmemsize_, smemsizemax_);
                                  exit(1);
                              }
                              if (rsmemsize_ >= smemsizemax_) {
                                  fprintf(stderr,
                                          "Shared memory consumption in f03_reducer (=%dbyte) "
                                          "reached the limit (=%dbyte). Too many fvars used.\n",
                                          rsmemsize_, smemsizemax_);
                                  exit(1);
                              }
                              nblockmax_ = prop.maxGridSize[0];
                              nsp_ = prop.multiProcessorCount * 8; // # of spream processors.
                              // fprintf(stderr, "nblockmax:%d  nsp:%d\n", nblockmax_, nsp_);
                          }
                          nisub_ = nimax_;
                          for (ioff_ = 0; ioff_ < n; ioff_ += nisub_) {
                              if (ioff_ + nisub_ > n) {
                                  nisub_ = n - ioff_;
                              }

                              int nisub_ru_ = ((nisub_ - 1) / nip_ + 1) * nip_;
                              if (njold != n || niold != nisub_) {
                                  // Adjust # of JP fragments so that large enough # of threads (to fill all SPs, 
                                  // and to hide latency of the global memory) are dispatched, at long as each
                                  // JP fragments has several times warp size.
                                  // You may want to hand tune the following part to obtain the optimal 'njdiv_'
                                  // value for a given hardware configuration.
                                  njdiv_ = njdiv_ru_ = 1;
                                  while (nisub_ru_ / nip_ * njdiv_ * 2 <= nblockmax_ &&
                                         nisub_ * njdiv_ < nsp_ * 200 &&
                                         n / njdiv_ > kbdim_ * 6 &&
                                         njdiv_ * 2 <= rbdim_) {
                                      njdiv_ *= 2;
                                  }
                                  // njdiv_ru_ is set to njdiv_ rounded up to a power of two.
                                  while (njdiv_ru_ < njdiv_) {
                                      njdiv_ru_ *= 2;
                                  }
                                  njold = n;
                                  niold = nisub_;
#if 0
                                  fprintf(stderr, "\n");
                                  fprintf(stderr, "nj:%d  njdiv:%d  nj/njdiv:%d  ni:%d  ni * njdiv:%d\n",
                                          n, njdiv_, n/njdiv_, nisub_, nisub_ * njdiv_);
                                  fprintf(stderr, "\n");
#endif
                              }
                              dim3 kthreads_(kbdim_, 1, 1);
                              dim3 kgrids_(njdiv_, nisub_ru_ / nip_, 1);
                              dim3 rthreads_(rbdim_, 1, 1);
                              int npipe_ = nisub_ru_ / 1;
                              int nrblock_ = npipe_ * njdiv_ / rbdim_;
                              dim3 rgrids_(nrblock_,1, 1);
                              int nj_ru_ = ((n - 0 - 1) / (kbdim_ * njdiv_) + 1) * (kbdim_ * njdiv_);
                              static int jbufsize_ = 0, ibufsize_ = 0, rsubbufsize_ = 0, rbufsize_ = 0;
                              static GlobalMem<f03_jp_t> f03_jp_;
                              static GlobalMem<f03_ip_t> f03_ip_;
                              static GlobalMem<f03_result_t> f03_result_;
                              static GlobalMem<f03_result_t> f03_result_sub_;
                              if (nj_ru_ > jbufsize_) {
                                  jbufsize_ = nj_ru_;
                                  f03_jp_.realloc(jbufsize_);
                              }

                              // Here we need to alloc n IPs. Note that npipe_ IPs would not be enough,
                              // since IPs on the device memory should not be overwritten during the ioff_ loop.
                              if (ni_ru_ > ibufsize_) {
                                  ibufsize_ = ni_ru_;
                                  f03_ip_.realloc(ibufsize_);
                              }
                              if (ni_ru_ > rbufsize_) {
                                  rbufsize_ = ni_ru_;
                                  f03_result_.realloc(rbufsize_);
                              }
                              if (npipe_ * njdiv_ > rsubbufsize_) {
                                  rsubbufsize_ = npipe_ * njdiv_ ;
                                  f03_result_sub_.realloc(rsubbufsize_);
                              }

                              for (j = 0, bufidx_ = 0 ; j <n; j++, bufidx_++) {f03_jp_[bufidx_].x_j_0_ = x[j][0];f03_jp_[bufidx_].x_j_1_ = x[j][1];f03_jp_[bufidx_].x_j_2_ = x[j][2];f03_jp_[bufidx_].m_j_ = m[j];}f03_jp_.htod(0, bufidx_);

                              if (ioff_ < n) { // otherwise nothing to send.
                                  int nipsend_ = npipe_;
                                  if (ioff_ + nipsend_ >= n) {
                                      nipsend_ = n - ioff_;
                                  }
                                  for (i = 0 ; i < nipsend_; i++) {f03_ip_[ioff_ + i].x_i_0_ = x[ioff_ + i][0];f03_ip_[ioff_ + i].x_i_1_ = x[ioff_ + i][1];f03_ip_[ioff_ + i].x_i_2_ = x[ioff_ + i][2];}f03_ip_.htod(ioff_, nipsend_);

                              }
                              f03_calculator<<<kgrids_, kthreads_, ksmemsize_>>>(ioff_, nisub_, n - 0, (double)eps2, f03_jp_, (f03_ip_ + ioff_), f03_result_sub_);
                              f03_reducer<<<rgrids_, rthreads_, rsmemsize_>>>(njdiv_, njdiv_ru_, (f03_result_ + ioff_), f03_result_sub_);

                              if (ioff_ < n) { // otherwise nothing to receive.
                                  int niprecv_ = npipe_;
                                  if (ioff_ + niprecv_ >= n) {
                                      niprecv_ = n - ioff_;
                                  }
                                  f03_result_.dtoh(ioff_, niprecv_);
for (i = 0 ; i < niprecv_; i++) {pot[ioff_ + i] = f03_result_[ioff_ + i].pot_i_;}
                              }
                          } // end of ioff loop.
                      } // end of api calls.
for (i=0;i<n;i++) {
    pot[i]+=m[i]/sqrt(eps2);
    pot[i]*=m[i];
  }


  for (i=0;i<n;i++) {
    pot_n+=pot[i];
  }

  pot_n*=0.5;



  pot_pn=0.0;

  for (i=1;i<n;i++) {
    for (k=0;k<3;k++) dx[k]=x[i][k]-x[0][k];

    r1a2e=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+eps2;
    r1ainv=rsqrt(r1a2e);
    mr1ainv=m[0]*r1ainv;

    vi2=v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2];

    pot_pn+=0.375*m[i]*vi2*vi2 + 1.5*mr1ainv*m[i]*vi2
      +0.5*mr1ainv*mr1ainv*m[i]
      +0.5*mr1ainv*m[i]
      *(3.0*vi2-7.0*(v[0][0]*v[i][0]+v[0][1]*v[i][1]+v[0][2]*v[i][2])
 -(dx[0]*v[0][0]+dx[1]*v[0][1]+dx[2]*v[0][2])
 *(dx[0]*v[i][0]+dx[1]*v[i][1]+dx[2]*v[i][2])
 *r1ainv*r1ainv);
  }

  for (i=1;i<n;i++) {
    pot_pn2[i] = 0.0;
  }

                      /*
                       * dispatcher for a kernel defined in 'Limit-sticky9_4.cu'.
                       *
                       * kbdim_  : # of threads per block. fixed to a multiple of warp size.
                       * nip_    : # of IPs handled in one block.
                       * njdiv_  : # of JP fragments. determined in a heuristic manner as a function of
                       *           ni, nj and device specification such as max grid size and warp size,
                       * nblock_ : # of blocks (ni * njdiv / nip) rounded up to a multiple of nip.
                       * npipe_  : len of each IP array. ni/nip_pack rounded up to a multiple of (nth * nip_pack).
                       * nj_ru_  : nj rounded up to a multiple of (nth * njdiv).
                       */
                      {
#if SHARED_HOSTBUF
                          double * GlobalMem<double>::hostbuf = NULL;
                          int GlobalMem<double>::nbytemax = 0;
#endif
                          // Below may have a room for optimization yet.
                          const int kbdim_ = 64; 
                          const int rbdim_ = 64; 
                          const int nimax_ = kbdim_ * 128; // necessary amount of main mem is larger for larger value.
                          const size_t ksmemsize_ = kbdim_ * sizeof(f04_jp_t);
                          const size_t rsmemsize_ = rbdim_ * sizeof(f04_result_t) * 2;
                          const int nip_ = kbdim_ * 1;
                          const GoosePrecision gprec_ = kGoosePrecisionDouble;
                          static int njdiv_, njdiv_ru_;
                          static int njold = 0, niold = 0;
                          static int nblockmax_, nsp_;
                          static size_t smemsizemax_;
                          static int firstcall_ = 1;
                          int bufidx_, ioff_, nisub_;
                          int ni_ru_ = ((n - 1 - 1) / nimax_ + 1) * nimax_;
                          cudaError cuerr_;
                          if (firstcall_) {
                              int device_;
                              cudaDeviceProp prop;
                              firstcall_ = 0;
                              cuerr_ = cudaGetDevice(&device_);
                              cutilSafeCall(cuerr_);
                              cuerr_ = cudaGetDeviceProperties(&prop, device_);
                              cutilSafeCall(cuerr_);
                              if (gprec_ >= kGoosePrecisionDouble) {
                                  if (prop.major < 2 && prop.minor < 3) {
                                      fprintf(stderr, "GPU architecture sm_%d%d does not support double-precision arithmetic.\n",
                                              prop.major, prop.minor);
                                      exit(1);
                                  }
                              }
                              smemsizemax_ = prop.sharedMemPerBlock;
                              if (ksmemsize_ >= smemsizemax_) {
                                  fprintf(stderr,
                                          "Shared memory consumption in f04_calculator (=%dbyte) "
                                          "reached the limit (=%dbyte). Too many jvars or fvars used.\n",
                                          ksmemsize_, smemsizemax_);
                                  exit(1);
                              }
                              if (rsmemsize_ >= smemsizemax_) {
                                  fprintf(stderr,
                                          "Shared memory consumption in f04_reducer (=%dbyte) "
                                          "reached the limit (=%dbyte). Too many fvars used.\n",
                                          rsmemsize_, smemsizemax_);
                                  exit(1);
                              }
                              nblockmax_ = prop.maxGridSize[0];
                              nsp_ = prop.multiProcessorCount * 8; // # of spream processors.
                              // fprintf(stderr, "nblockmax:%d  nsp:%d\n", nblockmax_, nsp_);
                          }
                          nisub_ = nimax_;
                          for (ioff_ = 1; ioff_ < n; ioff_ += nisub_) {
                              if (ioff_ + nisub_ > n) {
                                  nisub_ = n - ioff_;
                              }

                              int nisub_ru_ = ((nisub_ - 1) / nip_ + 1) * nip_;
                              if (njold != n || niold != nisub_) {
                                  // Adjust # of JP fragments so that large enough # of threads (to fill all SPs, 
                                  // and to hide latency of the global memory) are dispatched, at long as each
                                  // JP fragments has several times warp size.
                                  // You may want to hand tune the following part to obtain the optimal 'njdiv_'
                                  // value for a given hardware configuration.
                                  njdiv_ = njdiv_ru_ = 1;
                                  while (nisub_ru_ / nip_ * njdiv_ * 2 <= nblockmax_ &&
                                         nisub_ * njdiv_ < nsp_ * 200 &&
                                         n / njdiv_ > kbdim_ * 6 &&
                                         njdiv_ * 2 <= rbdim_) {
                                      njdiv_ *= 2;
                                  }
                                  // njdiv_ru_ is set to njdiv_ rounded up to a power of two.
                                  while (njdiv_ru_ < njdiv_) {
                                      njdiv_ru_ *= 2;
                                  }
                                  njold = n;
                                  niold = nisub_;
#if 0
                                  fprintf(stderr, "\n");
                                  fprintf(stderr, "nj:%d  njdiv:%d  nj/njdiv:%d  ni:%d  ni * njdiv:%d\n",
                                          n, njdiv_, n/njdiv_, nisub_, nisub_ * njdiv_);
                                  fprintf(stderr, "\n");
#endif
                              }
                              dim3 kthreads_(kbdim_, 1, 1);
                              dim3 kgrids_(njdiv_, nisub_ru_ / nip_, 1);
                              dim3 rthreads_(rbdim_, 1, 1);
                              int npipe_ = nisub_ru_ / 1;
                              int nrblock_ = npipe_ * njdiv_ / rbdim_;
                              dim3 rgrids_(nrblock_,1, 1);
                              int nj_ru_ = ((n - 1 - 1) / (kbdim_ * njdiv_) + 1) * (kbdim_ * njdiv_);
                              static int jbufsize_ = 0, ibufsize_ = 0, rsubbufsize_ = 0, rbufsize_ = 0;
                              static GlobalMem<f04_jp_t> f04_jp_;
                              static GlobalMem<f04_ip_t> f04_ip_;
                              static GlobalMem<f04_result_t> f04_result_;
                              static GlobalMem<f04_result_t> f04_result_sub_;
                              if (nj_ru_ > jbufsize_) {
                                  jbufsize_ = nj_ru_;
                                  f04_jp_.realloc(jbufsize_);
                              }

                              // Here we need to alloc n IPs. Note that npipe_ IPs would not be enough,
                              // since IPs on the device memory should not be overwritten during the ioff_ loop.
                              if (ni_ru_ > ibufsize_) {
                                  ibufsize_ = ni_ru_;
                                  f04_ip_.realloc(ibufsize_);
                              }
                              if (ni_ru_ > rbufsize_) {
                                  rbufsize_ = ni_ru_;
                                  f04_result_.realloc(rbufsize_);
                              }
                              if (npipe_ * njdiv_ > rsubbufsize_) {
                                  rsubbufsize_ = npipe_ * njdiv_ ;
                                  f04_result_sub_.realloc(rsubbufsize_);
                              }

                              for (j = 1, bufidx_ = 0 ; j <n; j++, bufidx_++) {f04_jp_[bufidx_].x_j_0_ = x[j][0];f04_jp_[bufidx_].x_j_1_ = x[j][1];f04_jp_[bufidx_].x_j_2_ = x[j][2];f04_jp_[bufidx_].m_j_ = m[j];f04_jp_[bufidx_].v_j_0_ = v[j][0];f04_jp_[bufidx_].v_j_1_ = v[j][1];f04_jp_[bufidx_].v_j_2_ = v[j][2];}f04_jp_.htod(0, bufidx_);

                              if (ioff_ < n) { // otherwise nothing to send.
                                  int nipsend_ = npipe_;
                                  if (ioff_ + nipsend_ >= n) {
                                      nipsend_ = n - ioff_;
                                  }
                                  for (i = 0 ; i < nipsend_; i++) {f04_ip_[ioff_ + i].x_i_0_ = x[ioff_ + i][0];f04_ip_[ioff_ + i].x_i_1_ = x[ioff_ + i][1];f04_ip_[ioff_ + i].x_i_2_ = x[ioff_ + i][2];f04_ip_[ioff_ + i].v_i_0_ = v[ioff_ + i][0];f04_ip_[ioff_ + i].v_i_1_ = v[ioff_ + i][1];f04_ip_[ioff_ + i].v_i_2_ = v[ioff_ + i][2];f04_ip_[ioff_ + i].m_i_ = m[ioff_ + i];}f04_ip_.htod(ioff_, nipsend_);

                              }
                              f04_calculator<<<kgrids_, kthreads_, ksmemsize_>>>(ioff_, nisub_, n - 1, (double)x[0][0], (double)x[0][1], (double)x[0][2], (double)eps2, (double)m[0], f04_jp_, (f04_ip_ + ioff_), f04_result_sub_);
                              f04_reducer<<<rgrids_, rthreads_, rsmemsize_>>>(njdiv_, njdiv_ru_, (f04_result_ + ioff_), f04_result_sub_);

                              if (ioff_ < n) { // otherwise nothing to receive.
                                  int niprecv_ = npipe_;
                                  if (ioff_ + niprecv_ >= n) {
                                      niprecv_ = n - ioff_;
                                  }
                                  f04_result_.dtoh(ioff_, niprecv_);
for (i = 0 ; i < niprecv_; i++) {pot_pn2[ioff_ + i] = f04_result_[ioff_ + i].pot_pn2_i_;}
                              }
                          } // end of ioff loop.
                      } // end of api calls.


  for (i=1;i<n;i++) {
     for (k=0;k<3;k++) {
       dxa[k]=x[i][k]-x[0][k];
     }

      r1a2e=dxa[0]*dxa[0]+dxa[1]*dxa[1]+dxa[2]*dxa[2]+eps2;
      r1ainv=rsqrt(r1a2e);
     vi2=v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2];

     pot_pn2[i]-=0.25*m[i]*m[i]/eps*(-vi2)
                 +m[0]*m[i]*m[i]*(r1ainv/eps+0.5*r1ainv*r1ainv);
   }

  ene=kin_n+pot_n+pot_pn;

  for (i=1;i<n;i++) {
    ene+=pot_pn2[i];
  }

  DE=(init_ene-ene)/ene;

  printf("time = %g\n",t);
  printf("pot = %22.15e kin = %22.15e \n pot_pn = %22.15e \n  total= %22.15e ratio = %e\n",
  pot_n, kin_n, pot_pn, ene, kin_n/pot_n);
  printf(" DE = %e %g\n",DE,t);
}

void initial_energy(double x[270000][3],
                    double v[270000][3],
                    double m[270000],
                    int n,
                    double *init_ene,
                    double eps)
{
  double pot[270000], ene, eps2;
  double dx[3], dxa[3], dxb[3], dxab[3];
  double r2, rinv, mrinv, vi2;
  double r1a2e, r1ainv, mr1ainv, r1b2e, r1binv, mr1binv;
  double rab2e, rabinv;
  double kin_n, pot_n;
  double pot_pn, pot_pn2[270000];
  int i,j,k;
  double cm[3];

  eps2=eps*eps;
  kin_n=0.0;
  pot_n=0.0;
  pot_pn=0.0;



  for (i=0;i<n;i++) {
    kin_n+=m[i]*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]);
  }

  kin_n*=0.5;

  for (i=0;i<n;i++) {
    pot[i]=0.0;
  }

                      /*
                       * dispatcher for a kernel defined in 'Limit-sticky9_5.cu'.
                       *
                       * kbdim_  : # of threads per block. fixed to a multiple of warp size.
                       * nip_    : # of IPs handled in one block.
                       * njdiv_  : # of JP fragments. determined in a heuristic manner as a function of
                       *           ni, nj and device specification such as max grid size and warp size,
                       * nblock_ : # of blocks (ni * njdiv / nip) rounded up to a multiple of nip.
                       * npipe_  : len of each IP array. ni/nip_pack rounded up to a multiple of (nth * nip_pack).
                       * nj_ru_  : nj rounded up to a multiple of (nth * njdiv).
                       */
                      {
#if SHARED_HOSTBUF
                          double * GlobalMem<double>::hostbuf = NULL;
                          int GlobalMem<double>::nbytemax = 0;
#endif
                          // Below may have a room for optimization yet.
                          const int kbdim_ = 64; 
                          const int rbdim_ = 64; 
                          const int nimax_ = kbdim_ * 128; // necessary amount of main mem is larger for larger value.
                          const size_t ksmemsize_ = kbdim_ * sizeof(f05_jp_t);
                          const size_t rsmemsize_ = rbdim_ * sizeof(f05_result_t) * 2;
                          const int nip_ = kbdim_ * 1;
                          const GoosePrecision gprec_ = kGoosePrecisionDouble;
                          static int njdiv_, njdiv_ru_;
                          static int njold = 0, niold = 0;
                          static int nblockmax_, nsp_;
                          static size_t smemsizemax_;
                          static int firstcall_ = 1;
                          int bufidx_, ioff_, nisub_;
                          int ni_ru_ = ((n - 0 - 1) / nimax_ + 1) * nimax_;
                          cudaError cuerr_;
                          if (firstcall_) {
                              int device_;
                              cudaDeviceProp prop;
                              firstcall_ = 0;
                              cuerr_ = cudaGetDevice(&device_);
                              cutilSafeCall(cuerr_);
                              cuerr_ = cudaGetDeviceProperties(&prop, device_);
                              cutilSafeCall(cuerr_);
                              if (gprec_ >= kGoosePrecisionDouble) {
                                  if (prop.major < 2 && prop.minor < 3) {
                                      fprintf(stderr, "GPU architecture sm_%d%d does not support double-precision arithmetic.\n",
                                              prop.major, prop.minor);
                                      exit(1);
                                  }
                              }
                              smemsizemax_ = prop.sharedMemPerBlock;
                              if (ksmemsize_ >= smemsizemax_) {
                                  fprintf(stderr,
                                          "Shared memory consumption in f05_calculator (=%dbyte) "
                                          "reached the limit (=%dbyte). Too many jvars or fvars used.\n",
                                          ksmemsize_, smemsizemax_);
                                  exit(1);
                              }
                              if (rsmemsize_ >= smemsizemax_) {
                                  fprintf(stderr,
                                          "Shared memory consumption in f05_reducer (=%dbyte) "
                                          "reached the limit (=%dbyte). Too many fvars used.\n",
                                          rsmemsize_, smemsizemax_);
                                  exit(1);
                              }
                              nblockmax_ = prop.maxGridSize[0];
                              nsp_ = prop.multiProcessorCount * 8; // # of spream processors.
                              // fprintf(stderr, "nblockmax:%d  nsp:%d\n", nblockmax_, nsp_);
                          }
                          nisub_ = nimax_;
                          for (ioff_ = 0; ioff_ < n; ioff_ += nisub_) {
                              if (ioff_ + nisub_ > n) {
                                  nisub_ = n - ioff_;
                              }

                              int nisub_ru_ = ((nisub_ - 1) / nip_ + 1) * nip_;
                              if (njold != n || niold != nisub_) {
                                  // Adjust # of JP fragments so that large enough # of threads (to fill all SPs, 
                                  // and to hide latency of the global memory) are dispatched, at long as each
                                  // JP fragments has several times warp size.
                                  // You may want to hand tune the following part to obtain the optimal 'njdiv_'
                                  // value for a given hardware configuration.
                                  njdiv_ = njdiv_ru_ = 1;
                                  while (nisub_ru_ / nip_ * njdiv_ * 2 <= nblockmax_ &&
                                         nisub_ * njdiv_ < nsp_ * 200 &&
                                         n / njdiv_ > kbdim_ * 6 &&
                                         njdiv_ * 2 <= rbdim_) {
                                      njdiv_ *= 2;
                                  }
                                  // njdiv_ru_ is set to njdiv_ rounded up to a power of two.
                                  while (njdiv_ru_ < njdiv_) {
                                      njdiv_ru_ *= 2;
                                  }
                                  njold = n;
                                  niold = nisub_;
#if 0
                                  fprintf(stderr, "\n");
                                  fprintf(stderr, "nj:%d  njdiv:%d  nj/njdiv:%d  ni:%d  ni * njdiv:%d\n",
                                          n, njdiv_, n/njdiv_, nisub_, nisub_ * njdiv_);
                                  fprintf(stderr, "\n");
#endif
                              }
                              dim3 kthreads_(kbdim_, 1, 1);
                              dim3 kgrids_(njdiv_, nisub_ru_ / nip_, 1);
                              dim3 rthreads_(rbdim_, 1, 1);
                              int npipe_ = nisub_ru_ / 1;
                              int nrblock_ = npipe_ * njdiv_ / rbdim_;
                              dim3 rgrids_(nrblock_,1, 1);
                              int nj_ru_ = ((n - 0 - 1) / (kbdim_ * njdiv_) + 1) * (kbdim_ * njdiv_);
                              static int jbufsize_ = 0, ibufsize_ = 0, rsubbufsize_ = 0, rbufsize_ = 0;
                              static GlobalMem<f05_jp_t> f05_jp_;
                              static GlobalMem<f05_ip_t> f05_ip_;
                              static GlobalMem<f05_result_t> f05_result_;
                              static GlobalMem<f05_result_t> f05_result_sub_;
                              if (nj_ru_ > jbufsize_) {
                                  jbufsize_ = nj_ru_;
                                  f05_jp_.realloc(jbufsize_);
                              }

                              // Here we need to alloc n IPs. Note that npipe_ IPs would not be enough,
                              // since IPs on the device memory should not be overwritten during the ioff_ loop.
                              if (ni_ru_ > ibufsize_) {
                                  ibufsize_ = ni_ru_;
                                  f05_ip_.realloc(ibufsize_);
                              }
                              if (ni_ru_ > rbufsize_) {
                                  rbufsize_ = ni_ru_;
                                  f05_result_.realloc(rbufsize_);
                              }
                              if (npipe_ * njdiv_ > rsubbufsize_) {
                                  rsubbufsize_ = npipe_ * njdiv_ ;
                                  f05_result_sub_.realloc(rsubbufsize_);
                              }

                              for (j = 0, bufidx_ = 0 ; j <n; j++, bufidx_++) {f05_jp_[bufidx_].x_j_0_ = x[j][0];f05_jp_[bufidx_].x_j_1_ = x[j][1];f05_jp_[bufidx_].x_j_2_ = x[j][2];f05_jp_[bufidx_].m_j_ = m[j];}f05_jp_.htod(0, bufidx_);

                              if (ioff_ < n) { // otherwise nothing to send.
                                  int nipsend_ = npipe_;
                                  if (ioff_ + nipsend_ >= n) {
                                      nipsend_ = n - ioff_;
                                  }
                                  for (i = 0 ; i < nipsend_; i++) {f05_ip_[ioff_ + i].x_i_0_ = x[ioff_ + i][0];f05_ip_[ioff_ + i].x_i_1_ = x[ioff_ + i][1];f05_ip_[ioff_ + i].x_i_2_ = x[ioff_ + i][2];}f05_ip_.htod(ioff_, nipsend_);

                              }
                              f05_calculator<<<kgrids_, kthreads_, ksmemsize_>>>(ioff_, nisub_, n - 0, (double)eps2, f05_jp_, (f05_ip_ + ioff_), f05_result_sub_);
                              f05_reducer<<<rgrids_, rthreads_, rsmemsize_>>>(njdiv_, njdiv_ru_, (f05_result_ + ioff_), f05_result_sub_);

                              if (ioff_ < n) { // otherwise nothing to receive.
                                  int niprecv_ = npipe_;
                                  if (ioff_ + niprecv_ >= n) {
                                      niprecv_ = n - ioff_;
                                  }
                                  f05_result_.dtoh(ioff_, niprecv_);
for (i = 0 ; i < niprecv_; i++) {pot[ioff_ + i] = f05_result_[ioff_ + i].pot_i_;}
                              }
                          } // end of ioff loop.
                      } // end of api calls.
for (i=0;i<n;i++) {
    pot[i]+=m[i]/sqrt(eps2);
    pot[i]*=m[i];
  }


  for (i=0;i<n;i++) {
    pot_n+=pot[i];
  }

  pot_n*=0.5;



  pot_pn=0.0;

  for (i=1;i<n;i++) {
    for (k=0;k<3;k++) dx[k]=x[i][k]-x[0][k];

    r1a2e=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+eps2;
    r1ainv=rsqrt(r1a2e);
    mr1ainv=m[0]*r1ainv;

    vi2=v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2];

    pot_pn+=0.375*m[i]*vi2*vi2 + 1.5*mr1ainv*m[i]*vi2
      +0.5*mr1ainv*mr1ainv*m[i]
      +0.5*mr1ainv*m[i]
      *(3.0*vi2-7.0*(v[0][0]*v[i][0]+v[0][1]*v[i][1]+v[0][2]*v[i][2])
 -(dx[0]*v[0][0]+dx[1]*v[0][1]+dx[2]*v[0][2])
 *(dx[0]*v[i][0]+dx[1]*v[i][1]+dx[2]*v[i][2])
 *r1ainv*r1ainv);
  }





  for (i=1;i<n;i++) {
    pot_pn2[i] = 0.0;
  }

                      /*
                       * dispatcher for a kernel defined in 'Limit-sticky9_6.cu'.
                       *
                       * kbdim_  : # of threads per block. fixed to a multiple of warp size.
                       * nip_    : # of IPs handled in one block.
                       * njdiv_  : # of JP fragments. determined in a heuristic manner as a function of
                       *           ni, nj and device specification such as max grid size and warp size,
                       * nblock_ : # of blocks (ni * njdiv / nip) rounded up to a multiple of nip.
                       * npipe_  : len of each IP array. ni/nip_pack rounded up to a multiple of (nth * nip_pack).
                       * nj_ru_  : nj rounded up to a multiple of (nth * njdiv).
                       */
                      {
#if SHARED_HOSTBUF
                          double * GlobalMem<double>::hostbuf = NULL;
                          int GlobalMem<double>::nbytemax = 0;
#endif
                          // Below may have a room for optimization yet.
                          const int kbdim_ = 64; 
                          const int rbdim_ = 64; 
                          const int nimax_ = kbdim_ * 128; // necessary amount of main mem is larger for larger value.
                          const size_t ksmemsize_ = kbdim_ * sizeof(f06_jp_t);
                          const size_t rsmemsize_ = rbdim_ * sizeof(f06_result_t) * 2;
                          const int nip_ = kbdim_ * 1;
                          const GoosePrecision gprec_ = kGoosePrecisionDouble;
                          static int njdiv_, njdiv_ru_;
                          static int njold = 0, niold = 0;
                          static int nblockmax_, nsp_;
                          static size_t smemsizemax_;
                          static int firstcall_ = 1;
                          int bufidx_, ioff_, nisub_;
                          int ni_ru_ = ((n - 1 - 1) / nimax_ + 1) * nimax_;
                          cudaError cuerr_;
                          if (firstcall_) {
                              int device_;
                              cudaDeviceProp prop;
                              firstcall_ = 0;
                              cuerr_ = cudaGetDevice(&device_);
                              cutilSafeCall(cuerr_);
                              cuerr_ = cudaGetDeviceProperties(&prop, device_);
                              cutilSafeCall(cuerr_);
                              if (gprec_ >= kGoosePrecisionDouble) {
                                  if (prop.major < 2 && prop.minor < 3) {
                                      fprintf(stderr, "GPU architecture sm_%d%d does not support double-precision arithmetic.\n",
                                              prop.major, prop.minor);
                                      exit(1);
                                  }
                              }
                              smemsizemax_ = prop.sharedMemPerBlock;
                              if (ksmemsize_ >= smemsizemax_) {
                                  fprintf(stderr,
                                          "Shared memory consumption in f06_calculator (=%dbyte) "
                                          "reached the limit (=%dbyte). Too many jvars or fvars used.\n",
                                          ksmemsize_, smemsizemax_);
                                  exit(1);
                              }
                              if (rsmemsize_ >= smemsizemax_) {
                                  fprintf(stderr,
                                          "Shared memory consumption in f06_reducer (=%dbyte) "
                                          "reached the limit (=%dbyte). Too many fvars used.\n",
                                          rsmemsize_, smemsizemax_);
                                  exit(1);
                              }
                              nblockmax_ = prop.maxGridSize[0];
                              nsp_ = prop.multiProcessorCount * 8; // # of spream processors.
                              // fprintf(stderr, "nblockmax:%d  nsp:%d\n", nblockmax_, nsp_);
                          }
                          nisub_ = nimax_;
                          for (ioff_ = 1; ioff_ < n; ioff_ += nisub_) {
                              if (ioff_ + nisub_ > n) {
                                  nisub_ = n - ioff_;
                              }

                              int nisub_ru_ = ((nisub_ - 1) / nip_ + 1) * nip_;
                              if (njold != n || niold != nisub_) {
                                  // Adjust # of JP fragments so that large enough # of threads (to fill all SPs, 
                                  // and to hide latency of the global memory) are dispatched, at long as each
                                  // JP fragments has several times warp size.
                                  // You may want to hand tune the following part to obtain the optimal 'njdiv_'
                                  // value for a given hardware configuration.
                                  njdiv_ = njdiv_ru_ = 1;
                                  while (nisub_ru_ / nip_ * njdiv_ * 2 <= nblockmax_ &&
                                         nisub_ * njdiv_ < nsp_ * 200 &&
                                         n / njdiv_ > kbdim_ * 6 &&
                                         njdiv_ * 2 <= rbdim_) {
                                      njdiv_ *= 2;
                                  }
                                  // njdiv_ru_ is set to njdiv_ rounded up to a power of two.
                                  while (njdiv_ru_ < njdiv_) {
                                      njdiv_ru_ *= 2;
                                  }
                                  njold = n;
                                  niold = nisub_;
#if 0
                                  fprintf(stderr, "\n");
                                  fprintf(stderr, "nj:%d  njdiv:%d  nj/njdiv:%d  ni:%d  ni * njdiv:%d\n",
                                          n, njdiv_, n/njdiv_, nisub_, nisub_ * njdiv_);
                                  fprintf(stderr, "\n");
#endif
                              }
                              dim3 kthreads_(kbdim_, 1, 1);
                              dim3 kgrids_(njdiv_, nisub_ru_ / nip_, 1);
                              dim3 rthreads_(rbdim_, 1, 1);
                              int npipe_ = nisub_ru_ / 1;
                              int nrblock_ = npipe_ * njdiv_ / rbdim_;
                              dim3 rgrids_(nrblock_,1, 1);
                              int nj_ru_ = ((n - 1 - 1) / (kbdim_ * njdiv_) + 1) * (kbdim_ * njdiv_);
                              static int jbufsize_ = 0, ibufsize_ = 0, rsubbufsize_ = 0, rbufsize_ = 0;
                              static GlobalMem<f06_jp_t> f06_jp_;
                              static GlobalMem<f06_ip_t> f06_ip_;
                              static GlobalMem<f06_result_t> f06_result_;
                              static GlobalMem<f06_result_t> f06_result_sub_;
                              if (nj_ru_ > jbufsize_) {
                                  jbufsize_ = nj_ru_;
                                  f06_jp_.realloc(jbufsize_);
                              }

                              // Here we need to alloc n IPs. Note that npipe_ IPs would not be enough,
                              // since IPs on the device memory should not be overwritten during the ioff_ loop.
                              if (ni_ru_ > ibufsize_) {
                                  ibufsize_ = ni_ru_;
                                  f06_ip_.realloc(ibufsize_);
                              }
                              if (ni_ru_ > rbufsize_) {
                                  rbufsize_ = ni_ru_;
                                  f06_result_.realloc(rbufsize_);
                              }
                              if (npipe_ * njdiv_ > rsubbufsize_) {
                                  rsubbufsize_ = npipe_ * njdiv_ ;
                                  f06_result_sub_.realloc(rsubbufsize_);
                              }

                              for (j = 1, bufidx_ = 0 ; j <n; j++, bufidx_++) {f06_jp_[bufidx_].x_j_0_ = x[j][0];f06_jp_[bufidx_].x_j_1_ = x[j][1];f06_jp_[bufidx_].x_j_2_ = x[j][2];f06_jp_[bufidx_].m_j_ = m[j];f06_jp_[bufidx_].v_j_0_ = v[j][0];f06_jp_[bufidx_].v_j_1_ = v[j][1];f06_jp_[bufidx_].v_j_2_ = v[j][2];}f06_jp_.htod(0, bufidx_);

                              if (ioff_ < n) { // otherwise nothing to send.
                                  int nipsend_ = npipe_;
                                  if (ioff_ + nipsend_ >= n) {
                                      nipsend_ = n - ioff_;
                                  }
                                  for (i = 0 ; i < nipsend_; i++) {f06_ip_[ioff_ + i].x_i_0_ = x[ioff_ + i][0];f06_ip_[ioff_ + i].x_i_1_ = x[ioff_ + i][1];f06_ip_[ioff_ + i].x_i_2_ = x[ioff_ + i][2];f06_ip_[ioff_ + i].v_i_0_ = v[ioff_ + i][0];f06_ip_[ioff_ + i].v_i_1_ = v[ioff_ + i][1];f06_ip_[ioff_ + i].v_i_2_ = v[ioff_ + i][2];f06_ip_[ioff_ + i].m_i_ = m[ioff_ + i];}f06_ip_.htod(ioff_, nipsend_);

                              }
                              f06_calculator<<<kgrids_, kthreads_, ksmemsize_>>>(ioff_, nisub_, n - 1, (double)x[0][0], (double)x[0][1], (double)x[0][2], (double)eps2, (double)m[0], f06_jp_, (f06_ip_ + ioff_), f06_result_sub_);
                              f06_reducer<<<rgrids_, rthreads_, rsmemsize_>>>(njdiv_, njdiv_ru_, (f06_result_ + ioff_), f06_result_sub_);

                              if (ioff_ < n) { // otherwise nothing to receive.
                                  int niprecv_ = npipe_;
                                  if (ioff_ + niprecv_ >= n) {
                                      niprecv_ = n - ioff_;
                                  }
                                  f06_result_.dtoh(ioff_, niprecv_);
for (i = 0 ; i < niprecv_; i++) {pot_pn2[ioff_ + i] = f06_result_[ioff_ + i].pot_pn2_i_;}
                              }
                          } // end of ioff loop.
                      } // end of api calls.


  for (i=1;i<n;i++) {
     for (k=0;k<3;k++) {
       dxa[k]=x[i][k]-x[0][k];
     }

      r1a2e=dxa[0]*dxa[0]+dxa[1]*dxa[1]+dxa[2]*dxa[2]+eps2;
      r1ainv=rsqrt(r1a2e);
     vi2=v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2];

     pot_pn2[i]-=0.25*m[i]*m[i]/eps*(-vi2)
                 +m[0]*m[i]*m[i]*(r1ainv/eps+0.5*r1ainv*r1ainv);
   }

  ene=kin_n+pot_n+pot_pn;

  for (i=1;i<n;i++) {
    ene+=pot_pn2[i];
  }

  printf("time = %g\n",0.0);
  printf("pot = %22.15e kin = %22.15e \n pot_pn = %22.15e \n  total= %22.15e ratio = %e\n",
  pot_n, kin_n, pot_pn, ene, kin_n/pot_n);

  *init_ene = ene;
}
void push_velocity(double v[270000][3],
                   double a[270000][3],
                   double dt,
                   int n)
{
  int j,k;
  for(j=0;j<n;j++){
    for(k=0;k<3;k++) v[j][k] += dt*a[j][k];
  }
}
void push_position(double x[270000][3],
                   double v[270000][3],
                   double dt,
                   int n)
{
  int j,k;
  for(j=0;j<n;j++){
    for(k=0;k<3;k++) x[j][k] += dt*v[j][k];
  }
}
main()
{
  static double x[270000][3];
  static double v[270000][3];
  static double m[270000];
  static double a[270000][3];
  static double ah[270000][3];
  static double pot[270000];
  static double hx[4][270000][3];
  static double hv[4][270000][3];
  static double x1[270000][3];
  static double v1[270000][3];
  double dt,eps,init_ene,time;
  double deouttime,eouttime,idtinv,endtime,epsinv;
  double icm[3],ratio;
  FILE *fp2;
  int n,i,k,dim,symp;
  double lt=0.0, st=0.0;
  double hlt=0.0, hst=0.0;
  double holdtime;
  double xsize=2.0,rotint=2.0,sustained=0.0;
  int simid;
  double gintrps;
  double peak;
  FILE *fpinput, *fpout;
  double rr;
  float xtmp,ytmp,ztmp,vxtmp,vytmp,vztmp;
  double cub2;
  time=0.0;
  eps=0.001;
  dt=1.0/131072.0;
  n=10001;
  endtime=10.0;
  deouttime=0.005;
  fpout = fopen("Sphere-r10-1PN.dat","w");
  fpinput = fopen("Init10k-sphere.dat","r");
  if (!fpinput) {
    perror("data_input");
    exit(1);
  }
  x[0][0]=0.0;
  x[0][1]=0.0;
  x[0][2]=0.0;
  v[0][0]=0.0;
  v[0][1]=0.0;
  v[0][2]=0.0;
  printf ("Init \n");
  for(i=1; i<n; i++) {
    fscanf (fpinput, "%f %f %f \n", &xtmp,&ytmp,&ztmp);
    x[i][0]=xtmp;
    x[i][1]=ytmp;
    x[i][2]=ztmp;
  }
  for(i=1; i<n; i++) {
    fscanf (fpinput, "%f %f %f \n", &vxtmp,&vytmp,&vztmp);
    v[i][0]=vxtmp;
    v[i][1]=vytmp;
    v[i][2]=vztmp;
  }
  fclose(fpinput);
  m[0] =100.0/(double)n;
  for (i=1; i<n; i++) {
    m[i] = 1.0/(double)n;
  }
  printf ("start \n");
  printf("initialdata end\n");
  for(i=0;i<n;i++) {
    fprintf(fpout,"%lf %lf %lf \n",x[i][0],x[i][1],x[i][2]);
  }
  for(i=0;i<n;i++) {
    fprintf(fpout,"%lf %lf %lf \n",v[i][0],v[i][1],v[i][2]);
  }
  eouttime= time+deouttime;
  initial_energy(x,v,m,n,&init_ene,eps);
  while(time < endtime){
    static int step = 0;
    for(i=0;i<n;i++) {
      for(k=0;k<3;k++) {
 hx[0][i][k]=x[i][k];
 hv[0][i][k]=v[i][k];
      }
    }
    force(x,v,m,eps,a,pot,n);
    for(i=0;i<n;i++) {
      for(k=0;k<3;k++) {
 x1[i][k]=x[i][k]+0.5*v[i][k]*dt;
 v1[i][k]=v[i][k]+0.5*a[i][k]*dt;
 hx[1][i][k]=x1[i][k];
 hv[1][i][k]=v1[i][k];
      }
    }
    force(x1,v1,m,eps,a,pot,n);
    for(i=0;i<n;i++) {
      for(k=0;k<3;k++) {
        x1[i][k]=x[i][k]+0.5*v1[i][k]*dt;
 v1[i][k]=v[i][k]+0.5*a[i][k]*dt;
 hx[2][i][k]=x1[i][k];
 hv[2][i][k]=v1[i][k];
      }
    }
    force(x1,v1,m,eps,a,pot,n);
    for(i=0;i<n;i++) {
      for(k=0;k<3;k++) {
        x1[i][k]=x[i][k]+v1[i][k]*dt;
 v1[i][k]=v[i][k]+a[i][k]*dt;
 hx[3][i][k]=x1[i][k];
 hv[3][i][k]=v1[i][k];
      }
    }
    for(i=0;i<n;i++) {
      for(k=0;k<3;k++) {
 x[i][k]=(hx[0][i][k]+hx[3][i][k]
   +2.0*(hx[1][i][k]+hx[2][i][k]))/6.0;
 v[i][k]=(hv[0][i][k]+hv[3][i][k]
                 +2.0*(hv[1][i][k]+hv[2][i][k]))/6.0;
      }
    }
    time += dt;
    if( time >= eouttime) {
      energy(time,x,v,m,n,init_ene,eps);
      eouttime += deouttime;
      for(i=0;i<n;i++) {
 fprintf(fpout,"%lf %lf %lf \n",x[i][0],x[i][1],x[i][2]);
      }
      for(i=0;i<n;i++) {
 fprintf(fpout,"%lf %lf %lf \n",v[i][0],v[i][1],v[i][2]);
      }
    }
  }
  for(i=0;i<n;i++) {
    fprintf(fpout,"%lf %lf %lf \n",x[i][0],x[i][1],x[i][2]);
  }
  for(i=0;i<n;i++) {
    fprintf(fpout,"%lf %lf %lf \n",v[i][0],v[i][1],v[i][2]);
  }
  fclose(fpout);
  printf("%lf \n", time);
}
#include "Limit-sticky9_0.cu"
#include "Limit-sticky9_1.cu"
#include "Limit-sticky9_2.cu"
#include "Limit-sticky9_3.cu"
#include "Limit-sticky9_4.cu"
#include "Limit-sticky9_5.cu"
#include "Limit-sticky9_6.cu"
