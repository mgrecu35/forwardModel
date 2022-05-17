//  SFM 04/06/2013  Code module added in merge from M.Grecu's code
//
#include <stdlib.h>
#include <stdio.h>

unsigned char i, dataN[1441][1441];
unsigned char i, dataS[1441][1441];
short int dataElevN[1441][1441];
short int dataElevS[1441][1441];

void getwfractnc_(float *r,float *s,float *fract)
{
  int i=*r;
  int j=*s;
  float f1, f2;
  f1=*r-i;
  f2=*s-j;
  int i0,j0;
  i0=i;
  j0=j;
  *fract=(1-f1)*(1-f2)*dataN[j0][i0]+
    (f1)*(1-f2)*dataN[j0][i0+1]+
    (f1)*(f2)*dataN[j0+1][i0+1]+
    (1-f1)*(f2)*dataN[j0+1][i0];
  if(*fract>100.) *fract=100.;
}

void getelevnc_(float *r,float *s,float *elev)
{
  int i=*r;
  int j=*s;
  float f1, f2;
  f1=*r-i;
  f2=*s-j;
  int i0,j0;
  i0=i;
  j0=j;
  *elev=(1-f1)*(1-f2)*dataElevN[j0][i0]+
    (f1)*(1-f2)*dataElevN[j0][i0+1]+
    (f1)*(f2)*dataElevN[j0+1][i0+1]+
    (1-f1)*(f2)*dataElevN[j0+1][i0];
  //if(*fract>100.) *fract=100.;
}

void getelevsc_(float *r,float *s,float *elev)
{
  int i=*r;
  int j=*s;
  float f1, f2;
  f1=*r-i;
  f2=*s-j;
  int i0,j0;
  i0=i;
  j0=j;
  *elev=(1-f1)*(1-f2)*dataElevS[j0][i0]+
    (f1)*(1-f2)*dataElevS[j0][i0+1]+
    (f1)*(f2)*dataElevS[j0+1][i0+1]+
    (1-f1)*(f2)*dataElevS[j0+1][i0];
  //if(*fract>100.) *fract=100.;
}
void readnc_(void)
{
  FILE *f;
 
  f=fopen("AncData/Nh.gl_g.sds01.v4.17.bin","rb");
  fread(&dataN[0][0],1,1441*1441,f);
  fclose(f);
  f=fopen("AncData/Nhelev.pc","rb");
  fread(&dataElevN[0][0],2,1441*1441,f);
  fclose(f);
}


void getwfractsc_(float *r,float *s,float *fract)
{
  int i=*r;
  int j=*s;
  float f1, f2;
  f1=*r-i;
  f2=*s-j;
  int i0,j0;
  i0=i;
  j0=j;
  *fract=(1-f1)*(1-f2)*dataS[j0][i0]+
    (f1)*(1-f2)*dataS[j0][i0+1]+
    (f1)*(f2)*dataS[j0+1][i0+1]+
    (1-f1)*(f2)*dataS[j0+1][i0];
  if(*fract>100.) *fract=100.;
}
void readsc_(void)
{
  FILE *f;
 
  f=fopen("AncData/Sh.gl_g.sds01.v4.17.bin","rb");
  fread(&dataS[0][0],1,1441*1441,f);
  fclose(f);
  f=fopen("AncData/Shelev.pc","rb");
  fread(&dataElevS[0][0],2,1441*1441,f);
  fclose(f);
}
