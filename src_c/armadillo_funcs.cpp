#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

extern "C" void kgain_(float *xENS, float *yEns, float *dy, int *nx, int *ny, int *nEns, float *s,
		       float *dx)
{
  int i, j;
  mat xENSa(*nx,*nEns),yENSa(*ny,*nEns);
  vec dya(*ny);
  //    kGain=dot(covXY,linalg.inv(covYY+R*eye(ny)))
  //    xRet=xEns.mean()+dot(kGain,dy)
	
  //vec dy=vec(*n),sol;
  for(i=0;i<*nEns;i++)
    for(j=0;j<*nx;j++)
      {
	xENSa(j,i)=xENS[j*(*nEns)+i];
      }
  for(j=0;j<*nx;j++)
    {
      float xmean=0;
      for(i=0;i<*nEns;i++)
	xmean+=xENSa(j,i);
      xENSa(j,span(0,*nEns-1))-=(xmean/(*nEns));
    }
  for(i=0;i<*nEns;i++)
    for(j=0;j<*ny;j++)
      {
	yENSa(j,i)=yEns[j*(*nEns)+i];
      }

  for(j=0;j<*ny;j++)
    {
      float ymean=0;
      for(i=0;i<*nEns;i++)
	ymean+=yENSa(j,i);
      yENSa(j,span(0,*nEns-1))-=(ymean/(*nEns));
    }
  mat covYY=(yENSa*yENSa.t())/(*nEns-1);
  mat covXY=(xENSa*yENSa.t())/(*nEns-1);
  //cout<<covYY;
  for(j=0;j<*ny;j++)
    {
      covYY(j,j)+=*s;
      dya(j)=dy[j];
    }
  vec sol,dxa;
  sol=solve(covYY,dya);
  dxa=covXY*sol;
  for(i=0;i<*nx;i++)
    dx[i]=dxa(i);

  /*
  for(i=0;i<*n;i++)
    printf("%6.2f ",zobs[i]);
  printf("\n");
  for(i=0;i<*n;i++)
    printf("%6.2f ",z[i]);
  printf("\n");
  for(i=0;i<*n;i++)
    printf("%6.2f ",sol[i]);
  printf("\n");
  printf("****\n");
  */
}

extern "C" void gauss_newton_(float *dzdn, float *z, float *zobs, int *n, float *s,
                            float *dn)
{
  int i, j;
  mat gradZ(*n,*n),temp;
  vec dy=vec(*n),sol;
  for(i=0;i<*n;i++)
    for(j=0;j<*n;j++)
    {
      gradZ(i,j)=dzdn[j*(*n)+i];
    }

  temp=gradZ.t()*gradZ;
  for(i=0;i<*n;i++)
    {
      temp(i,i)+=*s;
      if(zobs[i]>10 && z[i]>5)
	dy(i)=zobs[i]-z[i];
      else
	dy(i)=0;
    }
  vec graddy=gradZ.t()*dy;
  sol=solve(temp,graddy);
  for(i=0;i<*n;i++)
    dn[i]=sol(i);

  /*
  for(i=0;i<*n;i++)
    printf("%6.2f ",zobs[i]);
  printf("\n");
  for(i=0;i<*n;i++)
    printf("%6.2f ",z[i]);
  printf("\n");
  for(i=0;i<*n;i++)
    printf("%6.2f ",sol[i]);
  printf("\n");
  printf("****\n");
  */
}


extern "C" void interp_arm_(float *x, float *y, int *n,
                            float *xi, float *yi, int *ni)
{
  int i;
  vec xa=vec(*n);
  vec ya=vec(*n);
  vec xia=vec(*ni);
  vec yia=vec(*ni);
  for(i=0;i<*n;i++)
    {
      xa(i)=x[i];
      ya(i)=y[i];
    }
  for(i=0;i<*ni;i++)
    xia(i)=xi[i];
  
  
  interp1(xa, ya, xia, yia);
  for(i=0;i<*ni;i++)
    {
      yi[i]=yia(i);
      printf("%i %g \n",i,yia(i));
    }
}

extern "C" void kgainc_(float *dtb, int *n, float *s, float *kgain)
{
  mat A(*n,*n);
  int i,j;
  for(i=0;i<*n;i++)
    {
      for(j=0;j<*n;j++)
        A(i,j)=dtb[i]*dtb[j];
      A(i,i)=A(i,i)+(*s);
    }
  //A.print("A=:");
  mat B=pinv(A,0.0001);
  //B.print("B=:");
  
  for(i=0;i<*n;i++)
    {
      kgain[i]=0;
      for(j=0;j<*n;j++)
        kgain[i]+=dtb[j]*B(j,i);
    }
}
//pinv=linalg.pinv(dot(dtb.T,dtb)+eye(6)*4)
//  kgain=dot(dtb,pinv)

extern "C" void interp_armi_(int *x, float *y, int *n,
                             int *xi, float *yi, int *ni)
{
  int i;
  vec xa=vec(*n);
  vec ya=vec(*n);
  vec xia=vec(*ni);
  vec yia;
  for(i=0;i<*n;i++)
    {
      xa(i)=x[i];
      ya(i)=y[i];
    }
  for(i=0;i<*ni;i++)
    {
      xia(i)=xi[i];
      //printf("%g \n",xi[i]);
    }
  //  xia.print("xi:");
  
  interp1(xa, ya, xia, yia);
  //printf(" %i %i \n",*n,*ni);
  for(i=0;i<*ni;i++)
    {
      yi[i]=yia(i);
      //  printf("%i %lg %lg \n",i,xia(i),yia(i));
    }
}
