/* code for nbody : Limit-sticky9.c 
   (with GRAPE, Runge-Kutta scheme)
   2017.2.16
*/
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
double rsqrt(double r2)
{
  return 1.0/sqrt(r2);
}

void
force(double (*x)[3], double (*v)[3], double *m, double eps,
	    double (*a)[3], double *pot, int n)
{
  double r, r2e, r2, reinv, rinv, mrinv, mr3inv, dx[3], a1[3], a1c[3];
  double ac[NMAX][3];
  double dxb[3], dxc[3], dxbc[3], dvbc[3];
  double r1b2, r1b2e, r1c2, r1c2e, rbc2, rbc2e, r1be, mr1b3e, r1ce, rbce;
  double eps2;
  int i, j, k;

  eps2=eps*eps;

// #pragma omp parallel for private(j,dx,r2,rinv,mrinv,mr3inv,k)

#pragma goose parallel for precision ("double") loopcounter(i, j) \
nip_pack(1) // pack 4 IPs for better performance on GRAPE-DR, 1 on GTX280 and HD4850.

  for(i=0;i<n;i++) {
    for(k=0;k<3;k++) a[i][k] = 0.0;
    for (j=0;j<n;j++) {
      for(k=0;k<3;k++) dx[k] = x[j][k] - x[i][k];
      r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] + eps2;
      rinv = rsqrt(r2);
      mrinv = m[j]*rinv;
      mr3inv = mrinv*rinv*rinv;
      a[i][0] += mr3inv * dx[0];
      a[i][1] += mr3inv * dx[1];
      a[i][2] += mr3inv * dx[2];      
    }
  }


  for(k=0;k<3;k++) a1[k] = 0.0;

//#pragma goose parallel for precision ("double") result (a1[0..2]) loopcounter(i) nip_pack(1) // pack 4 IPs for better performance on GRAPE-DR, 1 on GTX280 and HD4850.

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

// Post-Newtonian a_{1-cross}

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
#pragma goose parallel for precision ("double") result (a1c[0..2]) loopcounter(i, j) nip_pack(1) // pack 4 IPs for better performance on GRAPE-DR, 1 on GTX280 and HD4850.

  for(i=1;i<n;i++) {
    for(j=1;j<n;j++) {
      for(k=0;k<3;k++) {
	dxb[k]=x[i][k]-x[0][k];
	dxc[k]=x[j][k]-x[0][k];
	dxbc[k]=x[j][k]-x[i][k];
      }

      r1b2=dxb[0]*dxb[0]+dxb[1]*dxb[1]+dxb[2]*dxb[2];
      r1b2e=r1b2+eps2;
      r1be=rsqrt(r1b2e);

      r1c2=dxc[0]*dxc[0]+dxc[1]*dxc[1]+dxc[2]*dxc[2];
      r1c2e=r1c2+eps2;
      r1ce=rsqrt(r1c2e);

      rbc2e=dxbc[0]*dxbc[0]+dxbc[1]*dxbc[1]+dxbc[2]*dxbc[2]+eps2;
      rbce=rsqrt(rbc2e);
      mr1b3e=m[i]*r1be*r1be*r1be;

      a1c[0]+=mr1b3e*m[j]*dxb[0]
	* (4.0*r1ce + 1.25*rbce - 0.25*r1c2/rbc2e*rbce
	   + 0.25*r1b2/rbc2e*rbce)
	-3.5*(rbce*rbce*rbce)*r1be*m[i]*m[j]*dxbc[0]
	-mr1b3e*m[j]/m[0]
	*(4.0*(v[i][0]*v[j][0]+v[i][1]*v[j][1]+v[i][2]*v[j][2])
	  *dxb[0]
	  -3.0*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])
	  *v[j][0]
	  -4.0*(v[j][0]*dxb[0]+v[j][1]*dxb[1]+v[j][2]*dxb[2])
	  *v[i][0]);

      a1c[1]+=mr1b3e*m[j]*dxb[1]
	* (4.0*r1ce + 1.25*rbce - 0.25*r1c2/rbc2e*rbce
	   + 0.25*r1b2/rbc2e*rbce)
	-3.5*(rbce*rbce*rbce)*r1be*m[i]*m[j]*dxbc[1]
	-mr1b3e*m[j]/m[0]
	*(4.0*(v[i][0]*v[j][0]+v[i][1]*v[j][1]+v[i][2]*v[j][2])
	  *dxb[1]
	  -3.0*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])
	  *v[j][1]
	  -4.0*(v[j][0]*dxb[0]+v[j][1]*dxb[1]+v[j][2]*dxb[2])
	  *v[i][1]);

      a1c[2]+=mr1b3e*m[j]*dxb[2]
	* (4.0*r1ce + 1.25*rbce - 0.25*r1c2/rbc2e*rbce
	   + 0.25*r1b2/rbc2e*rbce)
	-3.5*(rbce*rbce*rbce)*r1be*m[i]*m[j]*dxbc[2]
	-mr1b3e*m[j]/m[0]
	*(4.0*(v[i][0]*v[j][0]+v[i][1]*v[j][1]+v[i][2]*v[j][2])
	  *dxb[2]
	  -3.0*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])
	  *v[j][2]
	  -4.0*(v[j][0]*dxb[0]+v[j][1]*dxb[1]+v[j][2]*dxb[2])
	  *v[i][2]);
    }
  }

// Remove double-counting
// 2017.2.27

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

//a_{a-BH}

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

  //aa_{cross}

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

#pragma goose parallel for precision ("double") result (ac) loopcounter(i, j) nip_pack(1) // pack 4 IPs for better performance on GRAPE-DR, 1 on GTX280 and HD4850.

  for(i=1;i<n;i++){
    for (j=1;j<n;j++) {
      for(k=0;k<3;k++) {
	dxb[k]=x[0][k]-x[i][k];
	dxc[k]=x[0][k]-x[j][k];
        dxbc[k]=x[j][k]-x[i][k];
        dvbc[k]=v[j][k]-v[i][k];
      }

      r1b2=dxb[0]*dxb[0]+dxb[1]*dxb[1]+dxb[2]*dxb[2];
      r1b2e=r1b2+eps2;
      r1be=rsqrt(r1b2e);

      r1c2=dxc[0]*dxc[0]+dxc[1]*dxc[1]+dxc[2]*dxc[2];
      r1c2e=r1c2+eps2;
      r1ce=rsqrt(r1c2e);

      rbc2=dxbc[0]*dxbc[0]+dxbc[1]*dxbc[1]+dxbc[2]*dxbc[2];
      rbc2e=rbc2+eps2;
      rbce=rsqrt(rbc2e);

      mr3inv=m[j]*rbce*rbce*rbce;

      ac[i][0]+=m[j]*r1be*r1be*r1be*m[0]
	*(4.0*rbce+1.25*r1ce+0.25*(r1b2-rbc2)*r1ce*r1ce*r1ce)
	*dxb[0]
	+mr3inv*m[0]
	*(4.0*r1be+1.25*r1ce+0.25*(-r1b2+rbc2)*r1ce*r1ce*r1ce)
	*dxbc[0]
	-3.5*m[j]*m[0]*r1ce*r1ce*r1ce*(rbce-r1be)
	*dxc[0]
	-m[j]*r1be*r1be*r1be
	*(4.0*(v[i][0]*v[j][0]+v[i][1]*v[j][1]+v[i][2]*v[j][2])
	 *dxb[0]
	 -3.0*(v[j][0]*dxb[0]+v[j][1]*dxb[1]+v[j][2]*dxb[2])
	 *v[i][0]
	 -4.0*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])
	 *v[j][0])
	+mr3inv
	*((v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2])
	 -2.0*(dvbc[0]*dvbc[0]+dvbc[1]*dvbc[1]+dvbc[2]*dvbc[2])
	 +1.5*(v[j][0]*dxbc[0]+v[j][1]*dxbc[1]+v[j][2]*dxbc[2])
	 *(v[j][0]*dxbc[0]+v[j][1]*dxbc[1]+v[j][2]*dxbc[2])
	 /rbc2e)*dxbc[0]
	+mr3inv
	*(dxbc[0]*(4.0*v[i][0]-3.0*v[j][0])
	 +dxbc[1]*(4.0*v[i][1]-3.0*v[j][1])
	 +dxbc[2]*(4.0*v[i][2]-3.0*v[j][2]))
	*dvbc[0];

      ac[i][1]+=m[j]*r1be*r1be*r1be*m[0]
        *(4.0*rbce+1.25*r1ce+0.25*(r1b2-rbc2)*r1ce*r1ce*r1ce)
        *dxb[1]
        +mr3inv*m[0]
        *(4.0*r1be+1.25*r1ce+0.25*(-r1b2+rbc2)*r1ce*r1ce*r1ce)
        *dxbc[1]
        -3.5*m[j]*m[0]*r1ce*r1ce*r1ce*(rbce-r1be)
        *dxc[1]
        -m[j]*r1be*r1be*r1be
        *(4.0*(v[i][0]*v[j][0]+v[i][1]*v[j][1]+v[i][2]*v[j][2])
         *dxb[1]
         -3.0*(v[j][0]*dxb[0]+v[j][1]*dxb[1]+v[j][2]*dxb[2])
         *v[i][1]
         -4.0*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])
         *v[j][1])
        +mr3inv
	    *((v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2])
         -2.0*(dvbc[0]*dvbc[0]+dvbc[1]*dvbc[1]+dvbc[2]*dvbc[2])
         +1.5*(v[j][0]*dxbc[0]+v[j][1]*dxbc[1]+v[j][2]*dxbc[2])
         *(v[j][0]*dxbc[0]+v[j][1]*dxbc[1]+v[j][2]*dxbc[2])
     	 /rbc2e)*dxbc[1]
        +mr3inv
        *(dxbc[0]*(4.0*v[i][0]-3.0*v[j][0])
         +dxbc[1]*(4.0*v[i][1]-3.0*v[j][1])
         +dxbc[2]*(4.0*v[i][2]-3.0*v[j][2]))
        *dvbc[1];

      ac[i][2]+=m[j]*r1be*r1be*r1be*m[0]
        *(4.0*rbce+1.25*r1ce+0.25*(r1b2-rbc2)*r1ce*r1ce*r1ce)
        *dxb[2]
        +mr3inv*m[0]
        *(4.0*r1be+1.25*r1ce+0.25*(-r1b2+rbc2)*r1ce*r1ce*r1ce)
        *dxbc[2]
        -3.5*m[j]*m[0]*r1ce*r1ce*r1ce*(rbce-r1be)
        *dxc[2]
        -m[j]*r1be*r1be*r1be
        *(4.0*(v[i][0]*v[j][0]+v[i][1]*v[j][1]+v[i][2]*v[j][2])
         *dxb[2]
         -3.0*(v[j][0]*dxb[0]+v[j][1]*dxb[1]+v[j][2]*dxb[2])
         *v[i][2]
         -4.0*(v[i][0]*dxb[0]+v[i][1]*dxb[1]+v[i][2]*dxb[2])
         *v[j][2])
        +mr3inv
	    *((v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2])
         -2.0*(dvbc[0]*dvbc[0]+dvbc[1]*dvbc[1]+dvbc[2]*dvbc[2])
         +1.5*(v[j][0]*dxbc[0]+v[j][1]*dxbc[1]+v[j][2]*dxbc[2])
         *(v[j][0]*dxbc[0]+v[j][1]*dxbc[1]+v[j][2]*dxbc[2])
	 /rbc2e)*dxbc[2]
        +mr3inv
        *(dxbc[0]*(4.0*v[i][0]-3.0*v[j][0])
         +dxbc[1]*(4.0*v[i][1]-3.0*v[j][1])
         +dxbc[2]*(4.0*v[i][2]-3.0*v[j][2]))
        *dvbc[2];
    }
  }

// Remove double-counting

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


void energy(REAL t,
            REAL x[NMAX][DIM],
            REAL v[NMAX][DIM],
            REAL m[NMAX],
            int n,
            REAL init_ene,
            double eps)
{
  REAL pot[NMAX], DE, ene, eps2;
  REAL dx[3], dxa[3], dxb[3], dxab[3];
  REAL r2, rinv, mrinv, vi2;
  REAL r1a2e, r1ainv, mr1ainv, r1b2e, r1binv, mr1binv;
  REAL rab2e, rabinv;
  REAL kin_n, pot_n;
  REAL pot_pn, pot_pn2[NMAX];
  int i,j,k;	
  REAL cm[DIM];

  eps2=eps*eps;
  kin_n=0.0;
  pot_n=0.0;
  pot_pn=0.0;

  //   Newtonian

  for (i=0;i<n;i++) {
    kin_n+=m[i]*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]);
  }

  kin_n*=0.5;
  
  for (i=0;i<n;i++) {
    pot[i]=0.0;
  }

#pragma goose parallel for precision ("double") result (pot_n) loopcounter(i,j) nip_pack(1) // pack 4 IPs for better performance on GRAPE-DR, 1 on GTX280 and HD4850
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++){
      for (k=0;k<3;k++) dx[k]=x[j][k]-x[i][k];
      r2=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+eps2;

      rinv=rsqrt(r2);
      mrinv=rinv*m[j];
      pot[i]-=mrinv;
    }
    pot[i]+=m[i]/sqrt(eps2);
    pot[i]*=m[i];
  }

  for (i=0;i<n;i++) {
    pot_n+=pot[i];
  }

  pot_n*=0.5;

  // PN

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

#pragma goose parallel for precision ("double") result (pot_pn2) loopcounter(i,j) nip_pack(1) // pack 4 IPs for better performance on GRAPE-DR, 1 on GTX280 and  HD4850.
  for (i=1;i<n;i++) {
    for (j=1;j<n;j++) {
      for (k=0;k<3;k++) {
   	dxa[k]=x[i][k]-x[0][k];
        dxb[k]=x[j][k]-x[0][k];
        dxab[k]=x[j][k]-x[i][k];
      }

      r1a2e=dxa[0]*dxa[0]+dxa[1]*dxa[1]+dxa[2]*dxa[2]+eps2;
      r1ainv=rsqrt(r1a2e);

      r1b2e=dxb[0]*dxb[0]+dxb[1]*dxb[1]+dxb[2]*dxb[2]+eps2;
      rab2e=dxab[0]*dxab[0]+dxab[1]*dxab[1]+dxab[2]*dxab[2]+eps2;

      r1binv=rsqrt(r1b2e);
      rabinv=rsqrt(rab2e);

      vi2=v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2];

      pot_pn2[i]+=0.25*rabinv*m[i]*m[j]
	*(6.0*vi2
	  -7.0*(v[i][0]*v[j][0]+v[i][1]*v[j][1]+v[i][2]*v[j][2])
	  -(dxab[0]*v[i][0]+dxab[1]*v[i][1]+dxab[2]*v[i][2])
	  *(dxab[0]*v[j][0]+dxab[1]*v[j][1]+dxab[2]*v[j][2]) 
       /rab2e)
	+m[0]*m[i]*m[j]*r1ainv*rabinv
	+0.5*m[0]*m[i]*m[j]*r1ainv*r1binv;
    }
  }

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

void initial_energy(REAL x[NMAX][DIM],
                    REAL v[NMAX][DIM],
                    REAL m[NMAX],
                    int n,
                    REAL *init_ene,
                    double eps)
{
  REAL pot[NMAX], ene, eps2;
  REAL dx[3], dxa[3], dxb[3], dxab[3];
  REAL r2, rinv, mrinv, vi2;
  REAL r1a2e, r1ainv, mr1ainv, r1b2e, r1binv, mr1binv;
  REAL rab2e, rabinv;
  REAL kin_n, pot_n;
  REAL pot_pn, pot_pn2[NMAX];
  int i,j,k;	
  REAL cm[DIM];

  eps2=eps*eps;
  kin_n=0.0;
  pot_n=0.0;
  pot_pn=0.0;

  //   Newtonian

  for (i=0;i<n;i++) {
    kin_n+=m[i]*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]);
  }

  kin_n*=0.5;
  
  for (i=0;i<n;i++) {
    pot[i]=0.0;
  }

#pragma goose parallel for precision ("double") result (pot_n) loopcounter(i,j) nip_pack(1) // pack 4 IPs for better performance on GRAPE-DR, 1 on GTX280 and HD4850
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++){
      for (k=0;k<3;k++) dx[k]=x[j][k]-x[i][k];
      r2=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+eps2;

      rinv=rsqrt(r2);
      mrinv=rinv*m[j];
      pot[i]-=mrinv;
    }
    pot[i]+=m[i]/sqrt(eps2);
    pot[i]*=m[i];
  }

  for (i=0;i<n;i++) {
    pot_n+=pot[i];
  }

  pot_n*=0.5;

  // PN

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

  // printf("pot_pn1 = %22.15e \n", pot_pn);

  //  pot_pn2=0.0;

  for (i=1;i<n;i++) {
    pot_pn2[i] = 0.0;
  }

#pragma goose parallel for precision ("double") private(dxa,r1a2e,ra1inv,dxa,dxb,dxab,r1b2e,rab2e,r1binv,rabinv) result (pot_pn2) loopcounter(i,j) nip_pack(1) // pack 4 IPs for better performance on GRAPE-DR, 1 on GTX280 and  HD4850.
  for (i=1;i<n;i++) {
    for (j=1;j<n;j++) {
      for (k=0;k<3;k++) {
	dxa[k]=x[i][k]-x[0][k];
        dxb[k]=x[j][k]-x[0][k];
        dxab[k]=x[j][k]-x[i][k];
      }

      r1a2e=dxa[0]*dxa[0]+dxa[1]*dxa[1]+dxa[2]*dxa[2]+eps2;
      r1ainv=rsqrt(r1a2e);

      r1b2e=dxb[0]*dxb[0]+dxb[1]*dxb[1]+dxb[2]*dxb[2]+eps2;
      rab2e=dxab[0]*dxab[0]+dxab[1]*dxab[1]+dxab[2]*dxab[2]+eps2;

      r1binv=rsqrt(r1b2e);
      rabinv=rsqrt(rab2e);

      vi2=v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2];

      pot_pn2[i]+=0.25*rabinv*m[i]*m[j]
	*(6.0*vi2
	  -7.0*(v[i][0]*v[j][0]+v[i][1]*v[j][1]+v[i][2]*v[j][2])
	  -(dxab[0]*v[i][0]+dxab[1]*v[i][1]+dxab[2]*v[i][2])
	  *(dxab[0]*v[j][0]+dxab[1]*v[j][1]+dxab[2]*v[j][2]) 
       /rab2e)
	+m[0]*m[i]*m[j]*r1ainv*rabinv
	+0.5*m[0]*m[i]*m[j]*r1ainv*r1binv;

    }
  }

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

/*
  for(k=0;k<DIM;k++) cm[k] = 0.0;
  for(i=0;i<n;i++){
    for(k=0;k<DIM;k++) cm[k] += x[i][k] * m[i];
  }
*/

}

void push_velocity(REAL v[NMAX][DIM],
                   REAL a[NMAX][DIM],
                   REAL dt,
                   int n)
{
  int j,k;

  for(j=0;j<n;j++){
    for(k=0;k<DIM;k++) v[j][k] += dt*a[j][k];
  }
}

void push_position(REAL x[NMAX][DIM],
                   REAL v[NMAX][DIM],
                   REAL dt,
                   int n)
{
  int j,k;	

  for(j=0;j<n;j++){
    for(k=0;k<DIM;k++) x[j][k] += dt*v[j][k];
  }
}

main()
{
  static REAL x[NMAX][DIM];
  static REAL v[NMAX][DIM];
  static REAL m[NMAX];
  static REAL a[NMAX][DIM];
  static REAL ah[NMAX][DIM];
  static REAL pot[NMAX];
  static REAL hx[4][NMAX][DIM];
  static REAL hv[4][NMAX][DIM];
  static REAL x1[NMAX][DIM];
  static REAL v1[NMAX][DIM];
  REAL dt,eps,init_ene,time;
  REAL deouttime,eouttime,idtinv,endtime,epsinv;
  REAL icm[DIM],ratio;
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
    m[i] = 1.0/(double)(n-1);
  }

  printf ("start \n");

  for(i=0;i<n;i++) {
    fprintf(fpout,"%lf %lf %lf \n",x[i][0],x[i][1],x[i][2]);
  }

  for(i=0;i<n;i++) {
    fprintf(fpout,"%lf %lf %lf \n",v[i][0],v[i][1],v[i][2]);
  }

  eouttime= time+deouttime;
  //eouttime=endtime;

//    get_cputime(&hlt,&hst);
//    holdtime = hlt;

//  force(x,m,eps,a,pot,n);
//  a1_PN_force(x,v,m,eps,a,pot,n);

  initial_energy(x,v,m,n,&init_ene,eps);
//    fflush(stdout);

//    get_cputime(&lt,&st);
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
