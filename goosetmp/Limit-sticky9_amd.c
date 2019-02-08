#include <gcalutil.h>
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

#pragma goose func
double
rsqrt(double r2)
{
    return 1.0 / sqrt(r2);
}

void
force(double (*x)[3], double (*v)[3], double *m, double eps, double (*a)[3], double *pot, int n)
{
    double r, r2e, r2, reinv, rinv, mrinv, mr3inv, dx[3], a1[3], a1c[3];
    double ac[270000][3];
    double dxb[3], dxc[3], dxbc[3], dvbc[3];
    double r1b2, r1b2e, r1c2, r1c2e, rbc2, rbc2e, r1be, mr1b3e, r1ce, rbce;
    double eps2;
    int i, j, k;

    eps2 = eps * eps;

    for (i = 0; i < n; i++) {

        for (k = 0; k < 3; k++) {
            a[i][k] = 0.0;
        }

    }
    {
        // the following parameters can be adjusted to achive better performance.
        const int idomainw_ = GCAL_DATA_GRANULARITY * 4;
        const int jdomainw_ = 256;
        const int nimax_ = 65536;

        const int mindomainh_ = 2;      // lower limit for the domain height.
        const int looplen_ = 2;
        CALdomain jdomain_, idomain_;
        int bufidx_, ioff_, nisub_, npipe_;
        int hbufsize_req_;
        int moduleid_ = 0;      // and uses one kernel-module per device so far.

        int njcnt_ru_ = looplen_;
        int nj_write_ru_ = looplen_;

        njcnt_ru_ = (n - 0 - 1) / 1 + 1;        // j-loop count.
        njcnt_ru_ = ((njcnt_ru_ - 1) / looplen_ + 1) * looplen_;        // round it up to a multiple of looplen_.

        nj_write_ru_ = (n - 0 - 1) / 1 + 1;     // do the same for nj_write.
        nj_write_ru_ = ((nj_write_ru_ - 1) / looplen_ + 1) * looplen_;
        nj_write_ru_ = 0 + nj_write_ru_ * 1;    // convert j-loop count to j-index.

        static int hbufsize_ = 0;
        static int jbufsize_ = 0;
        static int ibufsize_ = 0;
        static double *doublebuf_ = NULL;

        // initialize the device.
        GCAL_open(moduleid_, f0_0_kernelsrc_, 11, 7);

        double constbuf_[] = { eps2, eps2, };
        static int constbuf_alloced_ = 0;
        if (!constbuf_alloced_) {
            GCAL_mallocj(moduleid_, 10, gcalCtypeDouble2, gcalMemConst, 1, 1);  // height must be 1.
            constbuf_alloced_ = 1;
        }
        GCAL_memcpyj(10, sizeof(double) * 2, constbuf_, gcalMemConst);

        // each domain handles looplen_ JPs, and thus we need (njcnt_ru_/looplen_) domains.
        jdomain_.width = jdomainw_;
        jdomain_.height = (njcnt_ru_ - 1) / (jdomain_.width * looplen_) + 1;
        if (jdomain_.height < mindomainh_) {
            jdomain_.height = mindomainh_;
        }

        // (re)allocate device-memory for jp.
        if (njcnt_ru_ > jbufsize_) {
            if (jbufsize_ > 0) {
                GCAL_free(3);
                GCAL_free(4);
                GCAL_free(5);
                GCAL_free(6);

            }
            jbufsize_ = njcnt_ru_;
            GCAL_mallocj(moduleid_, 3, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][0]
            GCAL_mallocj(moduleid_, 4, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][1]
            GCAL_mallocj(moduleid_, 5, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][2]
            GCAL_mallocj(moduleid_, 6, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // m[j]

        }
        // (re)alloc a buffer on the host memory.
        hbufsize_req_ = jbufsize_;
        if (hbufsize_req_ > hbufsize_) {
            hbufsize_ = hbufsize_req_;
            doublebuf_ = (double *) realloc(doublebuf_, sizeof(double) * hbufsize_);
            if (!doublebuf_) {
                perror("doublebuf_.");
                exit(1);
            }
        }

        // copy JPs from the host memory to the device memory.
        for (j = 0, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][0];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(3, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 0, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][1];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(4, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 0, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][2];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(5, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 0, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) m[j];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(6, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        nisub_ = nimax_ * GCAL_ndevice();
        for (ioff_ = 0; ioff_ < n; ioff_ += nisub_) {
            if (ioff_ + nisub_ > n) {
                nisub_ = n - ioff_;
            }
            npipe_ = (nisub_ - 1) / (1 * 1) + 1;        // # of threads to be dispatched.

            // each domain handles looplen_ IPs, and thus we need (npipe_/looplen_) domains.
            idomain_.width = idomainw_;
            idomain_.height = (npipe_ - 1) / (idomain_.width * looplen_) + 1;
            if (idomain_.height < mindomainh_) {
                idomain_.height = mindomainh_;
            }

            // (re)allocate device-memory for ip.
            if (npipe_ > ibufsize_) {
                if (ibufsize_ > 0) {
                    GCAL_free(0);
                    GCAL_free(1);
                    GCAL_free(2);
                    GCAL_free(7);
                    GCAL_free(8);
                    GCAL_free(9);

                }
                ibufsize_ = npipe_;
                GCAL_malloci(moduleid_, 0, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][0]
                GCAL_malloci(moduleid_, 1, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][1]
                GCAL_malloci(moduleid_, 2, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][2]
                GCAL_malloci(moduleid_, 7, gcalCtypeDouble2, gcalMemDeviceToHost, idomain_.width, idomain_.height);     // a[i][0]
                GCAL_malloci(moduleid_, 8, gcalCtypeDouble2, gcalMemDeviceToHost, idomain_.width, idomain_.height);     // a[i][1]
                GCAL_malloci(moduleid_, 9, gcalCtypeDouble2, gcalMemDeviceToHost, idomain_.width, idomain_.height);     // a[i][2]

            }
            // (re)alloc a buffer on the host memory.
            hbufsize_req_ = ibufsize_;
            if (hbufsize_req_ > hbufsize_) {
                hbufsize_ = hbufsize_req_;
                doublebuf_ = (double *) realloc(doublebuf_, sizeof(double) * hbufsize_);
                if (!doublebuf_) {
                    perror("doublebuf_.");
                    exit(1);
                }
            }

            // copy IPs from the host memory to the device memory.
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][0];
            }
            GCAL_memcpyi(0, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][1];
            }
            GCAL_memcpyi(1, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][2];
            }
            GCAL_memcpyi(2, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            // launch the kernel.
            GCAL_launch(moduleid_, 11, njcnt_ru_, looplen_, &jdomain_, &idomain_);

            // copy results from the device memory to the host memory.
            GCAL_memcpyi(7, sizeof(double) * npipe_, doublebuf_, gcalMemDeviceToHost);
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                a[ioff_ + i][0] = doublebuf_[bufidx_];
            }

            GCAL_memcpyi(8, sizeof(double) * npipe_, doublebuf_, gcalMemDeviceToHost);
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                a[ioff_ + i][1] = doublebuf_[bufidx_];
            }

            GCAL_memcpyi(9, sizeof(double) * npipe_, doublebuf_, gcalMemDeviceToHost);
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                a[ioff_ + i][2] = doublebuf_[bufidx_];
            }

        }                       // end of 'ioff_' loop.
        // finalize the device.
        GCAL_close();
        hbufsize_ = 0;
        jbufsize_ = 0;
        ibufsize_ = 0;
    }

    for (k = 0; k < 3; k++)
        a1[k] = 0.0;

    for (i = 1; i < n; i++) {
        for (k = 0; k < 3; k++) {
            dx[k] = x[i][k] - x[0][k];
        }
        r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
        rinv = rsqrt(r2);
        r2e = r2 + eps2;
        reinv = rsqrt(r2e);
        mrinv = m[i] * reinv;
        mr3inv = mrinv * reinv * reinv;

        a1[0] += mr3inv * dx[0]
            * (5.0 * m[0] * reinv - 2.0 * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2])
               + 1.5 * (v[i][0] * dx[0] + v[i][1] * dx[1] + v[i][2] * dx[2])
               * (v[i][0] * dx[0] + v[i][1] * dx[1] + v[i][2] * dx[2])
               / r2)
            + 3.0 * mr3inv * v[i][0]
            * (v[i][0] * dx[0] + v[i][1] * dx[1] + v[i][2] * dx[2]);

        a1[1] += mr3inv * dx[1]
            * (5.0 * m[0] * reinv - 2.0 * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2])
               + 1.5 * (v[i][0] * dx[0] + v[i][1] * dx[1] + v[i][2] * dx[2])
               * (v[i][0] * dx[0] + v[i][1] * dx[1] + v[i][2] * dx[2])
               / r2)
            + 3.0 * mr3inv * v[i][1]
            * (v[i][0] * dx[0] + v[i][1] * dx[1] + v[i][2] * dx[2]);

        a1[2] += mr3inv * dx[2]
            * (5.0 * m[0] * reinv - 2.0 * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2])
               + 1.5 * (v[i][0] * dx[0] + v[i][1] * dx[1] + v[i][2] * dx[2])
               * (v[i][0] * dx[0] + v[i][1] * dx[1] + v[i][2] * dx[2])
               / r2)
            + 3.0 * mr3inv * v[i][2]
            * (v[i][0] * dx[0] + v[i][1] * dx[1] + v[i][2] * dx[2]);
    }

    for (k = 0; k < 3; k++)
        a[0][k] += a1[k];

    for (k = 0; k < 3; k++)
        a1c[k] = 0.0;
    for (i = 1; i < n; i++) {
        for (k = 0; k < 3; k++) {
            dxb[k] = x[i][k] - x[0][k];
        }
        r1b2 = dxb[0] * dxb[0] + dxb[1] * dxb[1] + dxb[2] * dxb[2];
        r1b2e = r1b2 + eps2;
        r1be = rsqrt(r1b2e);
        mr1b3e = m[i] * r1be * r1be * r1be;

        a1c[0] += 4.0 * mr1b3e * r1be * m[i] * dxb[0];
        a1c[1] += 4.0 * mr1b3e * r1be * m[i] * dxb[1];
        a1c[2] += 4.0 * mr1b3e * r1be * m[i] * dxb[2];

    }
    {
        // the following parameters can be adjusted to achive better performance.
        const int idomainw_ = GCAL_DATA_GRANULARITY * 4;
        const int jdomainw_ = 256;
        const int nimax_ = 65536;

        const int mindomainh_ = 2;      // lower limit for the domain height.
        const int looplen_ = 2;
        CALdomain jdomain_, idomain_;
        int bufidx_, ioff_, nisub_, npipe_;
        int hbufsize_req_;
        int moduleid_ = 0;      // and uses one kernel-module per device so far.

        int njcnt_ru_ = looplen_;
        int nj_write_ru_ = looplen_;

        njcnt_ru_ = (n - 1 - 1) / 1 + 1;        // j-loop count.
        njcnt_ru_ = ((njcnt_ru_ - 1) / looplen_ + 1) * looplen_;        // round it up to a multiple of looplen_.

        nj_write_ru_ = (n - 1 - 1) / 1 + 1;     // do the same for nj_write.
        nj_write_ru_ = ((nj_write_ru_ - 1) / looplen_ + 1) * looplen_;
        nj_write_ru_ = 1 + nj_write_ru_ * 1;    // convert j-loop count to j-index.

        static int hbufsize_ = 0;
        static int jbufsize_ = 0;
        static int ibufsize_ = 0;
        static double *doublebuf_ = NULL;

        // initialize the device.
        GCAL_open(moduleid_, f0_1_kernelsrc_, 22, 14);

        double constbuf_[] = { x[0][0], x[0][0], x[0][1], x[0][1], x[0][2], x[0][2], eps2, eps2, m[0], m[0], };
        static int constbuf_alloced_ = 0;
        if (!constbuf_alloced_) {
            GCAL_mallocj(moduleid_, 17, gcalCtypeDouble2, gcalMemConst, 5, 1);  // height must be 1.
            constbuf_alloced_ = 1;
        }
        GCAL_memcpyj(17, sizeof(double) * 10, constbuf_, gcalMemConst);

        // each domain handles looplen_ JPs, and thus we need (njcnt_ru_/looplen_) domains.
        jdomain_.width = jdomainw_;
        jdomain_.height = (njcnt_ru_ - 1) / (jdomain_.width * looplen_) + 1;
        if (jdomain_.height < mindomainh_) {
            jdomain_.height = mindomainh_;
        }

        // (re)allocate device-memory for jp.
        if (njcnt_ru_ > jbufsize_) {
            if (jbufsize_ > 0) {
                GCAL_free(7);
                GCAL_free(8);
                GCAL_free(9);
                GCAL_free(10);
                GCAL_free(11);
                GCAL_free(12);
                GCAL_free(13);

            }
            jbufsize_ = njcnt_ru_;
            GCAL_mallocj(moduleid_, 7, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][0]
            GCAL_mallocj(moduleid_, 8, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][1]
            GCAL_mallocj(moduleid_, 9, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][2]
            GCAL_mallocj(moduleid_, 10, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height);        // m[j]
            GCAL_mallocj(moduleid_, 11, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height);        // v[j][0]
            GCAL_mallocj(moduleid_, 12, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height);        // v[j][1]
            GCAL_mallocj(moduleid_, 13, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height);        // v[j][2]

        }
        // (re)alloc a buffer on the host memory.
        hbufsize_req_ = jbufsize_;
        if (hbufsize_req_ > hbufsize_) {
            hbufsize_ = hbufsize_req_;
            doublebuf_ = (double *) realloc(doublebuf_, sizeof(double) * hbufsize_);
            if (!doublebuf_) {
                perror("doublebuf_.");
                exit(1);
            }
        }

        // copy JPs from the host memory to the device memory.
        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][0];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(7, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][1];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(8, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][2];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(9, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) m[j];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(10, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) v[j][0];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(11, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) v[j][1];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(12, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) v[j][2];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(13, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        nisub_ = nimax_ * GCAL_ndevice();
        for (ioff_ = 1; ioff_ < n; ioff_ += nisub_) {
            if (ioff_ + nisub_ > n) {
                nisub_ = n - ioff_;
            }
            npipe_ = (nisub_ - 1) / (1 * 1) + 1;        // # of threads to be dispatched.

            // each domain handles looplen_ IPs, and thus we need (npipe_/looplen_) domains.
            idomain_.width = idomainw_;
            idomain_.height = (npipe_ - 1) / (idomain_.width * looplen_) + 1;
            if (idomain_.height < mindomainh_) {
                idomain_.height = mindomainh_;
            }

            // (re)allocate device-memory for ip.
            if (npipe_ > ibufsize_) {
                if (ibufsize_ > 0) {
                    GCAL_free(0);
                    GCAL_free(1);
                    GCAL_free(2);
                    GCAL_free(3);
                    GCAL_free(4);
                    GCAL_free(5);
                    GCAL_free(6);
                    GCAL_free(14);
                    GCAL_free(15);
                    GCAL_free(16);

                }
                ibufsize_ = npipe_;
                GCAL_malloci(moduleid_, 0, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][0]
                GCAL_malloci(moduleid_, 1, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][1]
                GCAL_malloci(moduleid_, 2, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][2]
                GCAL_malloci(moduleid_, 3, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // m[i]
                GCAL_malloci(moduleid_, 4, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // v[i][0]
                GCAL_malloci(moduleid_, 5, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // v[i][1]
                GCAL_malloci(moduleid_, 6, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // v[i][2]
                GCAL_malloci(moduleid_, 14, gcalCtypeDouble2, gcalMemDeviceToHost, idomain_.width, idomain_.height);    // a1c[0]
                GCAL_malloci(moduleid_, 15, gcalCtypeDouble2, gcalMemDeviceToHost, idomain_.width, idomain_.height);    // a1c[1]
                GCAL_malloci(moduleid_, 16, gcalCtypeDouble2, gcalMemDeviceToHost, idomain_.width, idomain_.height);    // a1c[2]

            }
            // (re)alloc a buffer on the host memory.
            hbufsize_req_ = ibufsize_;
            if (hbufsize_req_ > hbufsize_) {
                hbufsize_ = hbufsize_req_;
                doublebuf_ = (double *) realloc(doublebuf_, sizeof(double) * hbufsize_);
                if (!doublebuf_) {
                    perror("doublebuf_.");
                    exit(1);
                }
            }

            // copy IPs from the host memory to the device memory.
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][0];
            }
            GCAL_memcpyi(0, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][1];
            }
            GCAL_memcpyi(1, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][2];
            }
            GCAL_memcpyi(2, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) m[ioff_ + i];
            }
            GCAL_memcpyi(3, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) v[ioff_ + i][0];
            }
            GCAL_memcpyi(4, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) v[ioff_ + i][1];
            }
            GCAL_memcpyi(5, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) v[ioff_ + i][2];
            }
            GCAL_memcpyi(6, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            // launch the kernel.
            GCAL_launch(moduleid_, 22, njcnt_ru_, looplen_, &jdomain_, &idomain_);

            // copy results from the device memory to the host memory.
            GCAL_memcpyi(14, sizeof(double) * npipe_, doublebuf_, gcalMemDeviceToHost);
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                a1c[0] += doublebuf_[bufidx_];
            }

            GCAL_memcpyi(15, sizeof(double) * npipe_, doublebuf_, gcalMemDeviceToHost);
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                a1c[1] += doublebuf_[bufidx_];
            }

            GCAL_memcpyi(16, sizeof(double) * npipe_, doublebuf_, gcalMemDeviceToHost);
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                a1c[2] += doublebuf_[bufidx_];
            }

        }                       // end of 'ioff_' loop.
        // finalize the device.
        GCAL_close();
        hbufsize_ = 0;
        jbufsize_ = 0;
        ibufsize_ = 0;
    }

    for (i = 1; i < n; i++) {
        for (k = 0; k < 3; k++) {
            dxb[k] = x[0][k] - x[i][k];
        }

        r1b2 = dxb[0] * dxb[0] + dxb[1] * dxb[1] + dxb[2] * dxb[2];
        r1b2e = r1b2 + eps2;
        r1be = rsqrt(r1b2e);

        mr1b3e = m[i] * r1be * r1be * r1be;

        a1c[0] -= mr1b3e * m[i] * dxb[0] * (4.0 * r1be + 1.25 * eps)
            - mr1b3e * m[i] / m[0]
            * (4.0 * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * dxb[0]
               - 7.0 * (v[i][0] * dxb[0] + v[i][1] * dxb[1] + v[i][2] * dxb[2]) * v[i][0]);

        a1c[1] -= mr1b3e * m[i] * dxb[1] * (4.0 * r1be + 1.25 * eps)
            - mr1b3e * m[i] / m[0]
            * (4.0 * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * dxb[1]
               - 7.0 * (v[i][0] * dxb[0] + v[i][1] * dxb[1] + v[i][2] * dxb[2]) * v[i][1]);

        a1c[2] -= mr1b3e * m[i] * dxb[2] * (4.0 * r1be + 1.25 * eps)
            - mr1b3e * m[i] / m[0]
            * (4.0 * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * dxb[2]
               - 7.0 * (v[i][0] * dxb[0] + v[i][1] * dxb[1] + v[i][2] * dxb[2]) * v[i][2]);
    }

    for (k = 0; k < 3; k++) {
        a[0][k] += a1c[k];
    }

    for (i = 1; i < n; i++) {
        for (k = 0; k < 3; k++) {
            dxb[k] = x[0][k] - x[i][k];
        }

        r1b2 = dxb[0] * dxb[0] + dxb[1] * dxb[1] + dxb[2] * dxb[2];
        r1b2e = r1b2 + eps2;
        r1be = rsqrt(r1b2e);

        mrinv = m[0] * r1be;
        mr3inv = mrinv * r1be * r1be;

        a[i][0] += mr3inv * dxb[0]
            * (4.0 * mrinv - (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]))
            + 4.0 * mr3inv * (v[i][0] * dxb[0] + v[i][1] * dxb[1] + v[i][2] * dxb[2]) * v[i][0];

        a[i][1] += mr3inv * dxb[1]
            * (4.0 * mrinv - (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]))
            + 4.0 * mr3inv * (v[i][0] * dxb[0] + v[i][1] * dxb[1] + v[i][2] * dxb[2]) * v[i][1];

        a[i][2] += mr3inv * dxb[2]
            * (4.0 * mrinv - (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]))
            + 4.0 * mr3inv * (v[i][0] * dxb[0] + v[i][1] * dxb[1] + v[i][2] * dxb[2]) * v[i][2];
    }

    for (i = 1; i < n; i++) {
        for (k = 0; k < 3; k++) {
            ac[i][k] = 0.0;
        }
    }

    for (i = 1; i < n; i++) {
        for (k = 0; k < 3; k++) {
            dxb[k] = x[i][k] - x[0][k];
        }

        r1b2 = dxb[0] * dxb[0] + dxb[1] * dxb[1] + dxb[2] * dxb[2];
        r1b2e = r1b2 + eps2;
        r1be = rsqrt(r1b2e);

        mr3inv = m[0] * r1be * r1be * r1be;

        a[i][0] += 5.0 * mr3inv * m[i] * r1be * dxb[0]
            - m[i] * r1be * r1be * r1be * (4.0 * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * dxb[0]
                                           - 7.0 * (v[i][0] * dxb[0] + v[i][1] * dxb[1] + v[i][2] * dxb[2]) * v[i][0]);

        a[i][1] += 5.0 * mr3inv * m[i] * r1be * dxb[1]
            - m[i] * r1be * r1be * r1be * (4.0 * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * dxb[1]
                                           - 7.0 * (v[i][0] * dxb[0] + v[i][1] * dxb[1] + v[i][2] * dxb[2]) * v[i][1]);

        a[i][2] += 5.0 * mr3inv * m[i] * r1be * dxb[2]
            - m[i] * r1be * r1be * r1be * (4.0 * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * dxb[2]
                                           - 7.0 * (v[i][0] * dxb[0] + v[i][1] * dxb[1] + v[i][2] * dxb[2]) * v[i][2]);
    }

    {
        // the following parameters can be adjusted to achive better performance.
        const int idomainw_ = GCAL_DATA_GRANULARITY * 4;
        const int jdomainw_ = 256;
        const int nimax_ = 65536;

        const int mindomainh_ = 2;      // lower limit for the domain height.
        const int looplen_ = 2;
        CALdomain jdomain_, idomain_;
        int bufidx_, ioff_, nisub_, npipe_;
        int hbufsize_req_;
        int moduleid_ = 0;      // and uses one kernel-module per device so far.

        int njcnt_ru_ = looplen_;
        int nj_write_ru_ = looplen_;

        njcnt_ru_ = (n - 1 - 1) / 1 + 1;        // j-loop count.
        njcnt_ru_ = ((njcnt_ru_ - 1) / looplen_ + 1) * looplen_;        // round it up to a multiple of looplen_.

        nj_write_ru_ = (n - 1 - 1) / 1 + 1;     // do the same for nj_write.
        nj_write_ru_ = ((nj_write_ru_ - 1) / looplen_ + 1) * looplen_;
        nj_write_ru_ = 1 + nj_write_ru_ * 1;    // convert j-loop count to j-index.

        static int hbufsize_ = 0;
        static int jbufsize_ = 0;
        static int ibufsize_ = 0;
        static double *doublebuf_ = NULL;

        // initialize the device.
        GCAL_open(moduleid_, f0_2_kernelsrc_, 21, 13);

        double constbuf_[] = { x[0][0], x[0][0], x[0][1], x[0][1], x[0][2], x[0][2], eps2, eps2, m[0], m[0], };
        static int constbuf_alloced_ = 0;
        if (!constbuf_alloced_) {
            GCAL_mallocj(moduleid_, 16, gcalCtypeDouble2, gcalMemConst, 5, 1);  // height must be 1.
            constbuf_alloced_ = 1;
        }
        GCAL_memcpyj(16, sizeof(double) * 10, constbuf_, gcalMemConst);

        // each domain handles looplen_ JPs, and thus we need (njcnt_ru_/looplen_) domains.
        jdomain_.width = jdomainw_;
        jdomain_.height = (njcnt_ru_ - 1) / (jdomain_.width * looplen_) + 1;
        if (jdomain_.height < mindomainh_) {
            jdomain_.height = mindomainh_;
        }

        // (re)allocate device-memory for jp.
        if (njcnt_ru_ > jbufsize_) {
            if (jbufsize_ > 0) {
                GCAL_free(6);
                GCAL_free(7);
                GCAL_free(8);
                GCAL_free(9);
                GCAL_free(10);
                GCAL_free(11);
                GCAL_free(12);

            }
            jbufsize_ = njcnt_ru_;
            GCAL_mallocj(moduleid_, 6, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][0]
            GCAL_mallocj(moduleid_, 7, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // v[j][0]
            GCAL_mallocj(moduleid_, 8, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][1]
            GCAL_mallocj(moduleid_, 9, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // v[j][1]
            GCAL_mallocj(moduleid_, 10, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height);        // x[j][2]
            GCAL_mallocj(moduleid_, 11, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height);        // v[j][2]
            GCAL_mallocj(moduleid_, 12, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height);        // m[j]

        }
        // (re)alloc a buffer on the host memory.
        hbufsize_req_ = jbufsize_;
        if (hbufsize_req_ > hbufsize_) {
            hbufsize_ = hbufsize_req_;
            doublebuf_ = (double *) realloc(doublebuf_, sizeof(double) * hbufsize_);
            if (!doublebuf_) {
                perror("doublebuf_.");
                exit(1);
            }
        }

        // copy JPs from the host memory to the device memory.
        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][0];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(6, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) v[j][0];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(7, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][1];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(8, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) v[j][1];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(9, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][2];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(10, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) v[j][2];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(11, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) m[j];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(12, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        nisub_ = nimax_ * GCAL_ndevice();
        for (ioff_ = 1; ioff_ < n; ioff_ += nisub_) {
            if (ioff_ + nisub_ > n) {
                nisub_ = n - ioff_;
            }
            npipe_ = (nisub_ - 1) / (1 * 1) + 1;        // # of threads to be dispatched.

            // each domain handles looplen_ IPs, and thus we need (npipe_/looplen_) domains.
            idomain_.width = idomainw_;
            idomain_.height = (npipe_ - 1) / (idomain_.width * looplen_) + 1;
            if (idomain_.height < mindomainh_) {
                idomain_.height = mindomainh_;
            }

            // (re)allocate device-memory for ip.
            if (npipe_ > ibufsize_) {
                if (ibufsize_ > 0) {
                    GCAL_free(0);
                    GCAL_free(1);
                    GCAL_free(2);
                    GCAL_free(3);
                    GCAL_free(4);
                    GCAL_free(5);
                    GCAL_free(13);
                    GCAL_free(14);
                    GCAL_free(15);

                }
                ibufsize_ = npipe_;
                GCAL_malloci(moduleid_, 0, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][0]
                GCAL_malloci(moduleid_, 1, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // v[i][0]
                GCAL_malloci(moduleid_, 2, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][1]
                GCAL_malloci(moduleid_, 3, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // v[i][1]
                GCAL_malloci(moduleid_, 4, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][2]
                GCAL_malloci(moduleid_, 5, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // v[i][2]
                GCAL_malloci(moduleid_, 13, gcalCtypeDouble2, gcalMemDeviceToHost, idomain_.width, idomain_.height);    // ac[i][0]
                GCAL_malloci(moduleid_, 14, gcalCtypeDouble2, gcalMemDeviceToHost, idomain_.width, idomain_.height);    // ac[i][1]
                GCAL_malloci(moduleid_, 15, gcalCtypeDouble2, gcalMemDeviceToHost, idomain_.width, idomain_.height);    // ac[i][2]

            }
            // (re)alloc a buffer on the host memory.
            hbufsize_req_ = ibufsize_;
            if (hbufsize_req_ > hbufsize_) {
                hbufsize_ = hbufsize_req_;
                doublebuf_ = (double *) realloc(doublebuf_, sizeof(double) * hbufsize_);
                if (!doublebuf_) {
                    perror("doublebuf_.");
                    exit(1);
                }
            }

            // copy IPs from the host memory to the device memory.
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][0];
            }
            GCAL_memcpyi(0, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) v[ioff_ + i][0];
            }
            GCAL_memcpyi(1, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][1];
            }
            GCAL_memcpyi(2, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) v[ioff_ + i][1];
            }
            GCAL_memcpyi(3, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][2];
            }
            GCAL_memcpyi(4, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) v[ioff_ + i][2];
            }
            GCAL_memcpyi(5, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            // launch the kernel.
            GCAL_launch(moduleid_, 21, njcnt_ru_, looplen_, &jdomain_, &idomain_);

            // copy results from the device memory to the host memory.
            GCAL_memcpyi(13, sizeof(double) * npipe_, doublebuf_, gcalMemDeviceToHost);
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                ac[ioff_ + i][0] = doublebuf_[bufidx_];
            }

            GCAL_memcpyi(14, sizeof(double) * npipe_, doublebuf_, gcalMemDeviceToHost);
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                ac[ioff_ + i][1] = doublebuf_[bufidx_];
            }

            GCAL_memcpyi(15, sizeof(double) * npipe_, doublebuf_, gcalMemDeviceToHost);
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                ac[ioff_ + i][2] = doublebuf_[bufidx_];
            }

        }                       // end of 'ioff_' loop.
        // finalize the device.
        GCAL_close();
        hbufsize_ = 0;
        jbufsize_ = 0;
        ibufsize_ = 0;
    }

    for (i = 1; i < n; i++) {
        for (k = 0; k < 3; k++) {
            a[i][k] += ac[i][k];
        }
    }

    for (i = 1; i < n; i++) {
        for (k = 0; k < 3; k++) {
            dxb[k] = x[0][k] - x[i][k];
        }

        r1b2 = dxb[0] * dxb[0] + dxb[1] * dxb[1] + dxb[2] * dxb[2];
        r1b2e = r1b2 + eps2;
        r1be = rsqrt(r1b2e);
        mr3inv = m[i] * r1be * r1be * r1be;

        a[i][0] -= mr3inv * m[0] * dxb[0]
            * (4.0 / eps + 1.25 * r1be + 0.25 * r1b2 * r1be * r1be * r1be)
            - 3.5 * mr3inv * m[0] * (1.0 / eps - r1be) * dxb[0]
            - mr3inv * (4.0 * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * dxb[0]
                        - 7.0 * (v[i][0] * dxb[0] + v[i][1] * dxb[1] + v[i][2] * dxb[2]) * v[i][0]);

        a[i][1] -= mr3inv * m[0] * dxb[1]
            * (4.0 / eps + 1.25 * r1be + 0.25 * r1b2 * r1be * r1be * r1be)
            - 3.5 * mr3inv * m[0] * (1.0 / eps - r1be) * dxb[1]
            - mr3inv * (4.0 * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * dxb[1]
                        - 7.0 * (v[i][0] * dxb[0] + v[i][1] * dxb[1] + v[i][2] * dxb[2]) * v[i][1]);

        a[i][2] -= mr3inv * m[0] * dxb[2]
            * (4.0 / eps + 1.25 * r1be + 0.25 * r1b2 * r1be * r1be * r1be)
            - 3.5 * mr3inv * m[0] * (1.0 / eps - r1be) * dxb[2]
            - mr3inv * (4.0 * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * dxb[2]
                        - 7.0 * (v[i][0] * dxb[0] + v[i][1] * dxb[1] + v[i][2] * dxb[2]) * v[i][2]);
    }

}

void
energy(double t, double x[270000][3], double v[270000][3], double m[270000], int n, double init_ene, double eps)
{
    double pot[270000], DE, ene, eps2;
    double dx[3], dxa[3], dxb[3], dxab[3];
    double r2, rinv, mrinv, vi2;
    double r1a2e, r1ainv, mr1ainv, r1b2e, r1binv, mr1binv;
    double rab2e, rabinv;
    double kin_n, pot_n;
    double pot_pn, pot_pn2[270000];
    int i, j, k;
    double cm[3];

    eps2 = eps * eps;
    kin_n = 0.0;
    pot_n = 0.0;
    pot_pn = 0.0;

    for (i = 0; i < n; i++) {
        kin_n += m[i] * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
    }

    kin_n *= 0.5;

    for (i = 0; i < n; i++) {
        pot[i] = 0.0;
    }

    {
        // the following parameters can be adjusted to achive better performance.
        const int idomainw_ = GCAL_DATA_GRANULARITY * 4;
        const int jdomainw_ = 256;
        const int nimax_ = 65536;

        const int mindomainh_ = 2;      // lower limit for the domain height.
        const int looplen_ = 2;
        CALdomain jdomain_, idomain_;
        int bufidx_, ioff_, nisub_, npipe_;
        int hbufsize_req_;
        int moduleid_ = 0;      // and uses one kernel-module per device so far.

        int njcnt_ru_ = looplen_;
        int nj_write_ru_ = looplen_;

        njcnt_ru_ = (n - 0 - 1) / 1 + 1;        // j-loop count.
        njcnt_ru_ = ((njcnt_ru_ - 1) / looplen_ + 1) * looplen_;        // round it up to a multiple of looplen_.

        nj_write_ru_ = (n - 0 - 1) / 1 + 1;     // do the same for nj_write.
        nj_write_ru_ = ((nj_write_ru_ - 1) / looplen_ + 1) * looplen_;
        nj_write_ru_ = 0 + nj_write_ru_ * 1;    // convert j-loop count to j-index.

        static int hbufsize_ = 0;
        static int jbufsize_ = 0;
        static int ibufsize_ = 0;
        static double *doublebuf_ = NULL;

        // initialize the device.
        GCAL_open(moduleid_, f0_3_kernelsrc_, 9, 7);

        double constbuf_[] = { eps2, eps2, };
        static int constbuf_alloced_ = 0;
        if (!constbuf_alloced_) {
            GCAL_mallocj(moduleid_, 8, gcalCtypeDouble2, gcalMemConst, 1, 1);   // height must be 1.
            constbuf_alloced_ = 1;
        }
        GCAL_memcpyj(8, sizeof(double) * 2, constbuf_, gcalMemConst);

        // each domain handles looplen_ JPs, and thus we need (njcnt_ru_/looplen_) domains.
        jdomain_.width = jdomainw_;
        jdomain_.height = (njcnt_ru_ - 1) / (jdomain_.width * looplen_) + 1;
        if (jdomain_.height < mindomainh_) {
            jdomain_.height = mindomainh_;
        }

        // (re)allocate device-memory for jp.
        if (njcnt_ru_ > jbufsize_) {
            if (jbufsize_ > 0) {
                GCAL_free(3);
                GCAL_free(4);
                GCAL_free(5);
                GCAL_free(6);

            }
            jbufsize_ = njcnt_ru_;
            GCAL_mallocj(moduleid_, 3, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][0]
            GCAL_mallocj(moduleid_, 4, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][1]
            GCAL_mallocj(moduleid_, 5, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][2]
            GCAL_mallocj(moduleid_, 6, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // m[j]

        }
        // (re)alloc a buffer on the host memory.
        hbufsize_req_ = jbufsize_;
        if (hbufsize_req_ > hbufsize_) {
            hbufsize_ = hbufsize_req_;
            doublebuf_ = (double *) realloc(doublebuf_, sizeof(double) * hbufsize_);
            if (!doublebuf_) {
                perror("doublebuf_.");
                exit(1);
            }
        }

        // copy JPs from the host memory to the device memory.
        for (j = 0, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][0];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(3, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 0, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][1];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(4, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 0, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][2];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(5, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 0, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) m[j];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(6, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        nisub_ = nimax_ * GCAL_ndevice();
        for (ioff_ = 0; ioff_ < n; ioff_ += nisub_) {
            if (ioff_ + nisub_ > n) {
                nisub_ = n - ioff_;
            }
            npipe_ = (nisub_ - 1) / (1 * 1) + 1;        // # of threads to be dispatched.

            // each domain handles looplen_ IPs, and thus we need (npipe_/looplen_) domains.
            idomain_.width = idomainw_;
            idomain_.height = (npipe_ - 1) / (idomain_.width * looplen_) + 1;
            if (idomain_.height < mindomainh_) {
                idomain_.height = mindomainh_;
            }

            // (re)allocate device-memory for ip.
            if (npipe_ > ibufsize_) {
                if (ibufsize_ > 0) {
                    GCAL_free(0);
                    GCAL_free(1);
                    GCAL_free(2);
                    GCAL_free(7);

                }
                ibufsize_ = npipe_;
                GCAL_malloci(moduleid_, 0, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][0]
                GCAL_malloci(moduleid_, 1, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][1]
                GCAL_malloci(moduleid_, 2, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][2]
                GCAL_malloci(moduleid_, 7, gcalCtypeDouble2, gcalMemDeviceToHost, idomain_.width, idomain_.height);     // pot[i]

            }
            // (re)alloc a buffer on the host memory.
            hbufsize_req_ = ibufsize_;
            if (hbufsize_req_ > hbufsize_) {
                hbufsize_ = hbufsize_req_;
                doublebuf_ = (double *) realloc(doublebuf_, sizeof(double) * hbufsize_);
                if (!doublebuf_) {
                    perror("doublebuf_.");
                    exit(1);
                }
            }

            // copy IPs from the host memory to the device memory.
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][0];
            }
            GCAL_memcpyi(0, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][1];
            }
            GCAL_memcpyi(1, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][2];
            }
            GCAL_memcpyi(2, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            // launch the kernel.
            GCAL_launch(moduleid_, 9, njcnt_ru_, looplen_, &jdomain_, &idomain_);

            // copy results from the device memory to the host memory.
            GCAL_memcpyi(7, sizeof(double) * npipe_, doublebuf_, gcalMemDeviceToHost);
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                pot[ioff_ + i] = doublebuf_[bufidx_];
            }

        }                       // end of 'ioff_' loop.
        // finalize the device.
        GCAL_close();
        hbufsize_ = 0;
        jbufsize_ = 0;
        ibufsize_ = 0;
    }
    for (i = 0; i < n; i++) {
        pot[i] += m[i] / sqrt(eps2);
        pot[i] *= m[i];
    }

    for (i = 0; i < n; i++) {
        pot_n += pot[i];
    }

    pot_n *= 0.5;

    pot_pn = 0.0;

    for (i = 1; i < n; i++) {
        for (k = 0; k < 3; k++)
            dx[k] = x[i][k] - x[0][k];

        r1a2e = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2] + eps2;
        r1ainv = rsqrt(r1a2e);
        mr1ainv = m[0] * r1ainv;

        vi2 = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];

        pot_pn += 0.375 * m[i] * vi2 * vi2 + 1.5 * mr1ainv * m[i] * vi2 + 0.5 * mr1ainv * mr1ainv * m[i]
            + 0.5 * mr1ainv * m[i]
            * (3.0 * vi2 - 7.0 * (v[0][0] * v[i][0] + v[0][1] * v[i][1] + v[0][2] * v[i][2])
               - (dx[0] * v[0][0] + dx[1] * v[0][1] + dx[2] * v[0][2])
               * (dx[0] * v[i][0] + dx[1] * v[i][1] + dx[2] * v[i][2])
               * r1ainv * r1ainv);
    }

    for (i = 1; i < n; i++) {
        pot_pn2[i] = 0.0;
    }

    {
        // the following parameters can be adjusted to achive better performance.
        const int idomainw_ = GCAL_DATA_GRANULARITY * 4;
        const int jdomainw_ = 256;
        const int nimax_ = 65536;

        const int mindomainh_ = 2;      // lower limit for the domain height.
        const int looplen_ = 2;
        CALdomain jdomain_, idomain_;
        int bufidx_, ioff_, nisub_, npipe_;
        int hbufsize_req_;
        int moduleid_ = 0;      // and uses one kernel-module per device so far.

        int njcnt_ru_ = looplen_;
        int nj_write_ru_ = looplen_;

        njcnt_ru_ = (n - 1 - 1) / 1 + 1;        // j-loop count.
        njcnt_ru_ = ((njcnt_ru_ - 1) / looplen_ + 1) * looplen_;        // round it up to a multiple of looplen_.

        nj_write_ru_ = (n - 1 - 1) / 1 + 1;     // do the same for nj_write.
        nj_write_ru_ = ((nj_write_ru_ - 1) / looplen_ + 1) * looplen_;
        nj_write_ru_ = 1 + nj_write_ru_ * 1;    // convert j-loop count to j-index.

        static int hbufsize_ = 0;
        static int jbufsize_ = 0;
        static int ibufsize_ = 0;
        static double *doublebuf_ = NULL;

        // initialize the device.
        GCAL_open(moduleid_, f0_4_kernelsrc_, 20, 14);

        double constbuf_[] = { x[0][0], x[0][0], x[0][1], x[0][1], x[0][2], x[0][2], eps2, eps2, m[0], m[0], };
        static int constbuf_alloced_ = 0;
        if (!constbuf_alloced_) {
            GCAL_mallocj(moduleid_, 15, gcalCtypeDouble2, gcalMemConst, 5, 1);  // height must be 1.
            constbuf_alloced_ = 1;
        }
        GCAL_memcpyj(15, sizeof(double) * 10, constbuf_, gcalMemConst);

        // each domain handles looplen_ JPs, and thus we need (njcnt_ru_/looplen_) domains.
        jdomain_.width = jdomainw_;
        jdomain_.height = (njcnt_ru_ - 1) / (jdomain_.width * looplen_) + 1;
        if (jdomain_.height < mindomainh_) {
            jdomain_.height = mindomainh_;
        }

        // (re)allocate device-memory for jp.
        if (njcnt_ru_ > jbufsize_) {
            if (jbufsize_ > 0) {
                GCAL_free(7);
                GCAL_free(8);
                GCAL_free(9);
                GCAL_free(10);
                GCAL_free(11);
                GCAL_free(12);
                GCAL_free(13);

            }
            jbufsize_ = njcnt_ru_;
            GCAL_mallocj(moduleid_, 7, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][0]
            GCAL_mallocj(moduleid_, 8, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][1]
            GCAL_mallocj(moduleid_, 9, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][2]
            GCAL_mallocj(moduleid_, 10, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height);        // m[j]
            GCAL_mallocj(moduleid_, 11, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height);        // v[j][0]
            GCAL_mallocj(moduleid_, 12, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height);        // v[j][1]
            GCAL_mallocj(moduleid_, 13, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height);        // v[j][2]

        }
        // (re)alloc a buffer on the host memory.
        hbufsize_req_ = jbufsize_;
        if (hbufsize_req_ > hbufsize_) {
            hbufsize_ = hbufsize_req_;
            doublebuf_ = (double *) realloc(doublebuf_, sizeof(double) * hbufsize_);
            if (!doublebuf_) {
                perror("doublebuf_.");
                exit(1);
            }
        }

        // copy JPs from the host memory to the device memory.
        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][0];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(7, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][1];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(8, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][2];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(9, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) m[j];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(10, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) v[j][0];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(11, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) v[j][1];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(12, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) v[j][2];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(13, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        nisub_ = nimax_ * GCAL_ndevice();
        for (ioff_ = 1; ioff_ < n; ioff_ += nisub_) {
            if (ioff_ + nisub_ > n) {
                nisub_ = n - ioff_;
            }
            npipe_ = (nisub_ - 1) / (1 * 1) + 1;        // # of threads to be dispatched.

            // each domain handles looplen_ IPs, and thus we need (npipe_/looplen_) domains.
            idomain_.width = idomainw_;
            idomain_.height = (npipe_ - 1) / (idomain_.width * looplen_) + 1;
            if (idomain_.height < mindomainh_) {
                idomain_.height = mindomainh_;
            }

            // (re)allocate device-memory for ip.
            if (npipe_ > ibufsize_) {
                if (ibufsize_ > 0) {
                    GCAL_free(0);
                    GCAL_free(1);
                    GCAL_free(2);
                    GCAL_free(3);
                    GCAL_free(4);
                    GCAL_free(5);
                    GCAL_free(6);
                    GCAL_free(14);

                }
                ibufsize_ = npipe_;
                GCAL_malloci(moduleid_, 0, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][0]
                GCAL_malloci(moduleid_, 1, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][1]
                GCAL_malloci(moduleid_, 2, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][2]
                GCAL_malloci(moduleid_, 3, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // v[i][0]
                GCAL_malloci(moduleid_, 4, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // v[i][1]
                GCAL_malloci(moduleid_, 5, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // v[i][2]
                GCAL_malloci(moduleid_, 6, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // m[i]
                GCAL_malloci(moduleid_, 14, gcalCtypeDouble2, gcalMemDeviceToHost, idomain_.width, idomain_.height);    // pot_pn2[i]

            }
            // (re)alloc a buffer on the host memory.
            hbufsize_req_ = ibufsize_;
            if (hbufsize_req_ > hbufsize_) {
                hbufsize_ = hbufsize_req_;
                doublebuf_ = (double *) realloc(doublebuf_, sizeof(double) * hbufsize_);
                if (!doublebuf_) {
                    perror("doublebuf_.");
                    exit(1);
                }
            }

            // copy IPs from the host memory to the device memory.
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][0];
            }
            GCAL_memcpyi(0, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][1];
            }
            GCAL_memcpyi(1, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][2];
            }
            GCAL_memcpyi(2, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) v[ioff_ + i][0];
            }
            GCAL_memcpyi(3, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) v[ioff_ + i][1];
            }
            GCAL_memcpyi(4, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) v[ioff_ + i][2];
            }
            GCAL_memcpyi(5, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) m[ioff_ + i];
            }
            GCAL_memcpyi(6, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            // launch the kernel.
            GCAL_launch(moduleid_, 20, njcnt_ru_, looplen_, &jdomain_, &idomain_);

            // copy results from the device memory to the host memory.
            GCAL_memcpyi(14, sizeof(double) * npipe_, doublebuf_, gcalMemDeviceToHost);
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                pot_pn2[ioff_ + i] = doublebuf_[bufidx_];
            }

        }                       // end of 'ioff_' loop.
        // finalize the device.
        GCAL_close();
        hbufsize_ = 0;
        jbufsize_ = 0;
        ibufsize_ = 0;
    }

    for (i = 1; i < n; i++) {
        for (k = 0; k < 3; k++) {
            dxa[k] = x[i][k] - x[0][k];
        }

        r1a2e = dxa[0] * dxa[0] + dxa[1] * dxa[1] + dxa[2] * dxa[2] + eps2;
        r1ainv = rsqrt(r1a2e);
        vi2 = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];

        pot_pn2[i] -= 0.25 * m[i] * m[i] / eps * (-vi2)
            + m[0] * m[i] * m[i] * (r1ainv / eps + 0.5 * r1ainv * r1ainv);
    }

    ene = kin_n + pot_n + pot_pn;

    for (i = 1; i < n; i++) {
        ene += pot_pn2[i];
    }

    DE = (init_ene - ene) / ene;

    printf("time = %g\n", t);
    printf("pot = %22.15e kin = %22.15e \n pot_pn = %22.15e \n  total= %22.15e ratio = %e\n", pot_n, kin_n, pot_pn, ene, kin_n / pot_n);
    printf(" DE = %e %g\n", DE, t);
}

void
initial_energy(double x[270000][3], double v[270000][3], double m[270000], int n, double *init_ene, double eps)
{
    double pot[270000], ene, eps2;
    double dx[3], dxa[3], dxb[3], dxab[3];
    double r2, rinv, mrinv, vi2;
    double r1a2e, r1ainv, mr1ainv, r1b2e, r1binv, mr1binv;
    double rab2e, rabinv;
    double kin_n, pot_n;
    double pot_pn, pot_pn2[270000];
    int i, j, k;
    double cm[3];

    eps2 = eps * eps;
    kin_n = 0.0;
    pot_n = 0.0;
    pot_pn = 0.0;

    for (i = 0; i < n; i++) {
        kin_n += m[i] * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
    }

    kin_n *= 0.5;

    for (i = 0; i < n; i++) {
        pot[i] = 0.0;
    }

    {
        // the following parameters can be adjusted to achive better performance.
        const int idomainw_ = GCAL_DATA_GRANULARITY * 4;
        const int jdomainw_ = 256;
        const int nimax_ = 65536;

        const int mindomainh_ = 2;      // lower limit for the domain height.
        const int looplen_ = 2;
        CALdomain jdomain_, idomain_;
        int bufidx_, ioff_, nisub_, npipe_;
        int hbufsize_req_;
        int moduleid_ = 0;      // and uses one kernel-module per device so far.

        int njcnt_ru_ = looplen_;
        int nj_write_ru_ = looplen_;

        njcnt_ru_ = (n - 0 - 1) / 1 + 1;        // j-loop count.
        njcnt_ru_ = ((njcnt_ru_ - 1) / looplen_ + 1) * looplen_;        // round it up to a multiple of looplen_.

        nj_write_ru_ = (n - 0 - 1) / 1 + 1;     // do the same for nj_write.
        nj_write_ru_ = ((nj_write_ru_ - 1) / looplen_ + 1) * looplen_;
        nj_write_ru_ = 0 + nj_write_ru_ * 1;    // convert j-loop count to j-index.

        static int hbufsize_ = 0;
        static int jbufsize_ = 0;
        static int ibufsize_ = 0;
        static double *doublebuf_ = NULL;

        // initialize the device.
        GCAL_open(moduleid_, f0_5_kernelsrc_, 9, 7);

        double constbuf_[] = { eps2, eps2, };
        static int constbuf_alloced_ = 0;
        if (!constbuf_alloced_) {
            GCAL_mallocj(moduleid_, 8, gcalCtypeDouble2, gcalMemConst, 1, 1);   // height must be 1.
            constbuf_alloced_ = 1;
        }
        GCAL_memcpyj(8, sizeof(double) * 2, constbuf_, gcalMemConst);

        // each domain handles looplen_ JPs, and thus we need (njcnt_ru_/looplen_) domains.
        jdomain_.width = jdomainw_;
        jdomain_.height = (njcnt_ru_ - 1) / (jdomain_.width * looplen_) + 1;
        if (jdomain_.height < mindomainh_) {
            jdomain_.height = mindomainh_;
        }

        // (re)allocate device-memory for jp.
        if (njcnt_ru_ > jbufsize_) {
            if (jbufsize_ > 0) {
                GCAL_free(3);
                GCAL_free(4);
                GCAL_free(5);
                GCAL_free(6);

            }
            jbufsize_ = njcnt_ru_;
            GCAL_mallocj(moduleid_, 3, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][0]
            GCAL_mallocj(moduleid_, 4, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][1]
            GCAL_mallocj(moduleid_, 5, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][2]
            GCAL_mallocj(moduleid_, 6, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // m[j]

        }
        // (re)alloc a buffer on the host memory.
        hbufsize_req_ = jbufsize_;
        if (hbufsize_req_ > hbufsize_) {
            hbufsize_ = hbufsize_req_;
            doublebuf_ = (double *) realloc(doublebuf_, sizeof(double) * hbufsize_);
            if (!doublebuf_) {
                perror("doublebuf_.");
                exit(1);
            }
        }

        // copy JPs from the host memory to the device memory.
        for (j = 0, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][0];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(3, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 0, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][1];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(4, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 0, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][2];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(5, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 0, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) m[j];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(6, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        nisub_ = nimax_ * GCAL_ndevice();
        for (ioff_ = 0; ioff_ < n; ioff_ += nisub_) {
            if (ioff_ + nisub_ > n) {
                nisub_ = n - ioff_;
            }
            npipe_ = (nisub_ - 1) / (1 * 1) + 1;        // # of threads to be dispatched.

            // each domain handles looplen_ IPs, and thus we need (npipe_/looplen_) domains.
            idomain_.width = idomainw_;
            idomain_.height = (npipe_ - 1) / (idomain_.width * looplen_) + 1;
            if (idomain_.height < mindomainh_) {
                idomain_.height = mindomainh_;
            }

            // (re)allocate device-memory for ip.
            if (npipe_ > ibufsize_) {
                if (ibufsize_ > 0) {
                    GCAL_free(0);
                    GCAL_free(1);
                    GCAL_free(2);
                    GCAL_free(7);

                }
                ibufsize_ = npipe_;
                GCAL_malloci(moduleid_, 0, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][0]
                GCAL_malloci(moduleid_, 1, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][1]
                GCAL_malloci(moduleid_, 2, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][2]
                GCAL_malloci(moduleid_, 7, gcalCtypeDouble2, gcalMemDeviceToHost, idomain_.width, idomain_.height);     // pot[i]

            }
            // (re)alloc a buffer on the host memory.
            hbufsize_req_ = ibufsize_;
            if (hbufsize_req_ > hbufsize_) {
                hbufsize_ = hbufsize_req_;
                doublebuf_ = (double *) realloc(doublebuf_, sizeof(double) * hbufsize_);
                if (!doublebuf_) {
                    perror("doublebuf_.");
                    exit(1);
                }
            }

            // copy IPs from the host memory to the device memory.
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][0];
            }
            GCAL_memcpyi(0, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][1];
            }
            GCAL_memcpyi(1, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][2];
            }
            GCAL_memcpyi(2, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            // launch the kernel.
            GCAL_launch(moduleid_, 9, njcnt_ru_, looplen_, &jdomain_, &idomain_);

            // copy results from the device memory to the host memory.
            GCAL_memcpyi(7, sizeof(double) * npipe_, doublebuf_, gcalMemDeviceToHost);
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                pot[ioff_ + i] = doublebuf_[bufidx_];
            }

        }                       // end of 'ioff_' loop.
        // finalize the device.
        GCAL_close();
        hbufsize_ = 0;
        jbufsize_ = 0;
        ibufsize_ = 0;
    }
    for (i = 0; i < n; i++) {
        pot[i] += m[i] / sqrt(eps2);
        pot[i] *= m[i];
    }

    for (i = 0; i < n; i++) {
        pot_n += pot[i];
    }

    pot_n *= 0.5;

    pot_pn = 0.0;

    for (i = 1; i < n; i++) {
        for (k = 0; k < 3; k++)
            dx[k] = x[i][k] - x[0][k];

        r1a2e = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2] + eps2;
        r1ainv = rsqrt(r1a2e);
        mr1ainv = m[0] * r1ainv;

        vi2 = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];

        pot_pn += 0.375 * m[i] * vi2 * vi2 + 1.5 * mr1ainv * m[i] * vi2 + 0.5 * mr1ainv * mr1ainv * m[i]
            + 0.5 * mr1ainv * m[i]
            * (3.0 * vi2 - 7.0 * (v[0][0] * v[i][0] + v[0][1] * v[i][1] + v[0][2] * v[i][2])
               - (dx[0] * v[0][0] + dx[1] * v[0][1] + dx[2] * v[0][2])
               * (dx[0] * v[i][0] + dx[1] * v[i][1] + dx[2] * v[i][2])
               * r1ainv * r1ainv);
    }

    for (i = 1; i < n; i++) {
        pot_pn2[i] = 0.0;
    }

    {
        // the following parameters can be adjusted to achive better performance.
        const int idomainw_ = GCAL_DATA_GRANULARITY * 4;
        const int jdomainw_ = 256;
        const int nimax_ = 65536;

        const int mindomainh_ = 2;      // lower limit for the domain height.
        const int looplen_ = 2;
        CALdomain jdomain_, idomain_;
        int bufidx_, ioff_, nisub_, npipe_;
        int hbufsize_req_;
        int moduleid_ = 0;      // and uses one kernel-module per device so far.

        int njcnt_ru_ = looplen_;
        int nj_write_ru_ = looplen_;

        njcnt_ru_ = (n - 1 - 1) / 1 + 1;        // j-loop count.
        njcnt_ru_ = ((njcnt_ru_ - 1) / looplen_ + 1) * looplen_;        // round it up to a multiple of looplen_.

        nj_write_ru_ = (n - 1 - 1) / 1 + 1;     // do the same for nj_write.
        nj_write_ru_ = ((nj_write_ru_ - 1) / looplen_ + 1) * looplen_;
        nj_write_ru_ = 1 + nj_write_ru_ * 1;    // convert j-loop count to j-index.

        static int hbufsize_ = 0;
        static int jbufsize_ = 0;
        static int ibufsize_ = 0;
        static double *doublebuf_ = NULL;

        // initialize the device.
        GCAL_open(moduleid_, f0_6_kernelsrc_, 20, 14);

        double constbuf_[] = { x[0][0], x[0][0], x[0][1], x[0][1], x[0][2], x[0][2], eps2, eps2, m[0], m[0], };
        static int constbuf_alloced_ = 0;
        if (!constbuf_alloced_) {
            GCAL_mallocj(moduleid_, 15, gcalCtypeDouble2, gcalMemConst, 5, 1);  // height must be 1.
            constbuf_alloced_ = 1;
        }
        GCAL_memcpyj(15, sizeof(double) * 10, constbuf_, gcalMemConst);

        // each domain handles looplen_ JPs, and thus we need (njcnt_ru_/looplen_) domains.
        jdomain_.width = jdomainw_;
        jdomain_.height = (njcnt_ru_ - 1) / (jdomain_.width * looplen_) + 1;
        if (jdomain_.height < mindomainh_) {
            jdomain_.height = mindomainh_;
        }

        // (re)allocate device-memory for jp.
        if (njcnt_ru_ > jbufsize_) {
            if (jbufsize_ > 0) {
                GCAL_free(7);
                GCAL_free(8);
                GCAL_free(9);
                GCAL_free(10);
                GCAL_free(11);
                GCAL_free(12);
                GCAL_free(13);

            }
            jbufsize_ = njcnt_ru_;
            GCAL_mallocj(moduleid_, 7, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][0]
            GCAL_mallocj(moduleid_, 8, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][1]
            GCAL_mallocj(moduleid_, 9, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height); // x[j][2]
            GCAL_mallocj(moduleid_, 10, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height);        // m[j]
            GCAL_mallocj(moduleid_, 11, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height);        // v[j][0]
            GCAL_mallocj(moduleid_, 12, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height);        // v[j][1]
            GCAL_mallocj(moduleid_, 13, gcalCtypeDouble2, gcalMemHostToDevice, jdomain_.width, jdomain_.height);        // v[j][2]

        }
        // (re)alloc a buffer on the host memory.
        hbufsize_req_ = jbufsize_;
        if (hbufsize_req_ > hbufsize_) {
            hbufsize_ = hbufsize_req_;
            doublebuf_ = (double *) realloc(doublebuf_, sizeof(double) * hbufsize_);
            if (!doublebuf_) {
                perror("doublebuf_.");
                exit(1);
            }
        }

        // copy JPs from the host memory to the device memory.
        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][0];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(7, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][1];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(8, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) x[j][2];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(9, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) m[j];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(10, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) v[j][0];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(11, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) v[j][1];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(12, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        for (j = 1, bufidx_ = 0; j < n; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) v[j][2];
        }
        for (j = n; n == n && j < nj_write_ru_; j += 1, bufidx_++) {
            doublebuf_[bufidx_] = (double) 0.0;
        }
        GCAL_memcpyj(13, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

        nisub_ = nimax_ * GCAL_ndevice();
        for (ioff_ = 1; ioff_ < n; ioff_ += nisub_) {
            if (ioff_ + nisub_ > n) {
                nisub_ = n - ioff_;
            }
            npipe_ = (nisub_ - 1) / (1 * 1) + 1;        // # of threads to be dispatched.

            // each domain handles looplen_ IPs, and thus we need (npipe_/looplen_) domains.
            idomain_.width = idomainw_;
            idomain_.height = (npipe_ - 1) / (idomain_.width * looplen_) + 1;
            if (idomain_.height < mindomainh_) {
                idomain_.height = mindomainh_;
            }

            // (re)allocate device-memory for ip.
            if (npipe_ > ibufsize_) {
                if (ibufsize_ > 0) {
                    GCAL_free(0);
                    GCAL_free(1);
                    GCAL_free(2);
                    GCAL_free(3);
                    GCAL_free(4);
                    GCAL_free(5);
                    GCAL_free(6);
                    GCAL_free(14);

                }
                ibufsize_ = npipe_;
                GCAL_malloci(moduleid_, 0, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][0]
                GCAL_malloci(moduleid_, 1, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][1]
                GCAL_malloci(moduleid_, 2, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // x[i][2]
                GCAL_malloci(moduleid_, 3, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // v[i][0]
                GCAL_malloci(moduleid_, 4, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // v[i][1]
                GCAL_malloci(moduleid_, 5, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // v[i][2]
                GCAL_malloci(moduleid_, 6, gcalCtypeDouble2, gcalMemHostToDevice, idomain_.width, idomain_.height);     // m[i]
                GCAL_malloci(moduleid_, 14, gcalCtypeDouble2, gcalMemDeviceToHost, idomain_.width, idomain_.height);    // pot_pn2[i]

            }
            // (re)alloc a buffer on the host memory.
            hbufsize_req_ = ibufsize_;
            if (hbufsize_req_ > hbufsize_) {
                hbufsize_ = hbufsize_req_;
                doublebuf_ = (double *) realloc(doublebuf_, sizeof(double) * hbufsize_);
                if (!doublebuf_) {
                    perror("doublebuf_.");
                    exit(1);
                }
            }

            // copy IPs from the host memory to the device memory.
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][0];
            }
            GCAL_memcpyi(0, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][1];
            }
            GCAL_memcpyi(1, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) x[ioff_ + i][2];
            }
            GCAL_memcpyi(2, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) v[ioff_ + i][0];
            }
            GCAL_memcpyi(3, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) v[ioff_ + i][1];
            }
            GCAL_memcpyi(4, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) v[ioff_ + i][2];
            }
            GCAL_memcpyi(5, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                doublebuf_[bufidx_] = (double) m[ioff_ + i];
            }
            GCAL_memcpyi(6, sizeof(double) * bufidx_, doublebuf_, gcalMemHostToDevice);

            // launch the kernel.
            GCAL_launch(moduleid_, 20, njcnt_ru_, looplen_, &jdomain_, &idomain_);

            // copy results from the device memory to the host memory.
            GCAL_memcpyi(14, sizeof(double) * npipe_, doublebuf_, gcalMemDeviceToHost);
            for (i = 0, bufidx_ = 0; i < nisub_; i += 1, bufidx_++) {
                pot_pn2[ioff_ + i] = doublebuf_[bufidx_];
            }

        }                       // end of 'ioff_' loop.
        // finalize the device.
        GCAL_close();
        hbufsize_ = 0;
        jbufsize_ = 0;
        ibufsize_ = 0;
    }

    for (i = 1; i < n; i++) {
        for (k = 0; k < 3; k++) {
            dxa[k] = x[i][k] - x[0][k];
        }

        r1a2e = dxa[0] * dxa[0] + dxa[1] * dxa[1] + dxa[2] * dxa[2] + eps2;
        r1ainv = rsqrt(r1a2e);
        vi2 = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];

        pot_pn2[i] -= 0.25 * m[i] * m[i] / eps * (-vi2)
            + m[0] * m[i] * m[i] * (r1ainv / eps + 0.5 * r1ainv * r1ainv);
    }

    ene = kin_n + pot_n + pot_pn;

    for (i = 1; i < n; i++) {
        ene += pot_pn2[i];
    }

    printf("time = %g\n", 0.0);
    printf("pot = %22.15e kin = %22.15e \n pot_pn = %22.15e \n  total= %22.15e ratio = %e\n", pot_n, kin_n, pot_pn, ene, kin_n / pot_n);

    *init_ene = ene;
}

void
push_velocity(double v[270000][3], double a[270000][3], double dt, int n)
{
    int j, k;
    for (j = 0; j < n; j++) {
        for (k = 0; k < 3; k++)
            v[j][k] += dt * a[j][k];
    }
}

void
push_position(double x[270000][3], double v[270000][3], double dt, int n)
{
    int j, k;
    for (j = 0; j < n; j++) {
        for (k = 0; k < 3; k++)
            x[j][k] += dt * v[j][k];
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
    double dt, eps, init_ene, time;
    double deouttime, eouttime, idtinv, endtime, epsinv;
    double icm[3], ratio;
    FILE *fp2;
    int n, i, k, dim, symp;
    double lt = 0.0, st = 0.0;
    double hlt = 0.0, hst = 0.0;
    double holdtime;
    double xsize = 2.0, rotint = 2.0, sustained = 0.0;
    int simid;
    double gintrps;
    double peak;
    FILE *fpinput, *fpout;
    double rr;
    float xtmp, ytmp, ztmp, vxtmp, vytmp, vztmp;
    double cub2;
    time = 0.0;
    eps = 0.001;
    dt = 1.0 / 131072.0;
    n = 10001;
    endtime = 1.0;
    deouttime = 0.005;
    fpout = fopen("Disk-r10-100-1PN.dat", "w");
    fpinput = fopen("Init10k-disk-r10-100.dat", "r");
    if (!fpinput) {
        perror("data_input");
        exit(1);
    }
    x[0][0] = 0.0;
    x[0][1] = 0.0;
    x[0][2] = 0.0;
    v[0][0] = 0.0;
    v[0][1] = 0.0;
    v[0][2] = 0.0;
    printf("Init \n");
    for (i = 1; i < n; i++) {
        fscanf(fpinput, "%f %f %f \n", &xtmp, &ytmp, &ztmp);
        x[i][0] = xtmp;
        x[i][1] = ytmp;
        x[i][2] = ztmp;
    }
    for (i = 1; i < n; i++) {
        fscanf(fpinput, "%f %f %f \n", &vxtmp, &vytmp, &vztmp);
        v[i][0] = vxtmp;
        v[i][1] = vytmp;
        v[i][2] = vztmp;
    }
    fclose(fpinput);
    m[0] = 100.0 / (double) n;
    for (i = 1; i < n; i++) {
        m[i] = 1.0 / (double) n;
    }
    printf("start \n");
    printf("initialdata end\n");
    for (i = 0; i < n; i++) {
        fprintf(fpout, "%lf %lf %lf \n", x[i][0], x[i][1], x[i][2]);
    }
    for (i = 0; i < n; i++) {
        fprintf(fpout, "%lf %lf %lf \n", v[i][0], v[i][1], v[i][2]);
    }
    eouttime = time + deouttime;
    initial_energy(x, v, m, n, &init_ene, eps);
    while (time < endtime) {
        static int step = 0;
        for (i = 0; i < n; i++) {
            for (k = 0; k < 3; k++) {
                hx[0][i][k] = x[i][k];
                hv[0][i][k] = v[i][k];
            }
        }
        force(x, v, m, eps, a, pot, n);
        for (i = 0; i < n; i++) {
            for (k = 0; k < 3; k++) {
                x1[i][k] = x[i][k] + 0.5 * v[i][k] * dt;
                v1[i][k] = v[i][k] + 0.5 * a[i][k] * dt;
                hx[1][i][k] = x1[i][k];
                hv[1][i][k] = v1[i][k];
            }
        }
        force(x1, v1, m, eps, a, pot, n);
        for (i = 0; i < n; i++) {
            for (k = 0; k < 3; k++) {
                x1[i][k] = x[i][k] + 0.5 * v1[i][k] * dt;
                v1[i][k] = v[i][k] + 0.5 * a[i][k] * dt;
                hx[2][i][k] = x1[i][k];
                hv[2][i][k] = v1[i][k];
            }
        }
        force(x1, v1, m, eps, a, pot, n);
        for (i = 0; i < n; i++) {
            for (k = 0; k < 3; k++) {
                x1[i][k] = x[i][k] + v1[i][k] * dt;
                v1[i][k] = v[i][k] + a[i][k] * dt;
                hx[3][i][k] = x1[i][k];
                hv[3][i][k] = v1[i][k];
            }
        }
        for (i = 0; i < n; i++) {
            for (k = 0; k < 3; k++) {
                x[i][k] = (hx[0][i][k] + hx[3][i][k]
                           + 2.0 * (hx[1][i][k] + hx[2][i][k])) / 6.0;
                v[i][k] = (hv[0][i][k] + hv[3][i][k]
                           + 2.0 * (hv[1][i][k] + hv[2][i][k])) / 6.0;
            }
        }
        time += dt;
        if (time >= eouttime) {
            energy(time, x, v, m, n, init_ene, eps);
            eouttime += deouttime;
            for (i = 0; i < n; i++) {
                fprintf(fpout, "%lf %lf %lf \n", x[i][0], x[i][1], x[i][2]);
            }
            for (i = 0; i < n; i++) {
                fprintf(fpout, "%lf %lf %lf \n", v[i][0], v[i][1], v[i][2]);
            }
        }
    }
    for (i = 0; i < n; i++) {
        fprintf(fpout, "%lf %lf %lf \n", x[i][0], x[i][1], x[i][2]);
    }
    for (i = 0; i < n; i++) {
        fprintf(fpout, "%lf %lf %lf \n", v[i][0], v[i][1], v[i][2]);
    }
    fclose(fpout);
    printf("%lf \n", time);
}
