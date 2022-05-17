#include <stdio.h>
#include <stdlib.h>
FILE *fenkf, *fenkf2;

void clearsc_()
{
 printf("\033[2J");
}
int openenkffile_(int *orbNumber)
{
  char fname[100],s[6];
  sprintf(s,"%6.6i",*orbNumber);
  strcpy(fname,"outdirITE109/enk2File");
  strcat(fname,s);
  fenkf=fopen(fname,"wb");
}

int openenkffilensl_(int *orbNumber,int *ichunk)
{
  char fname[100],s[9];
  sprintf(s,"%6.6i.%2.2i",*orbNumber,*ichunk);
  strcpy(fname,"outdirITE109/enk2File");
  strcat(fname,s);
  fenkf=fopen(fname,"wb");
}

int closeenkffile_()
{
  fclose(fenkf);
}

int closeenkffile2_()
{
  fclose(fenkf);
}

void enkfwi5_(int *n5)
{
  fwrite(n5,sizeof(int),5,fenkf);
}

void enkfwi1_(int *n1)
{
  fwrite(n1,sizeof(int),1,fenkf);
}

void enkfwf_(float *f)
{
  fwrite(f,sizeof(float),1,fenkf);
}
