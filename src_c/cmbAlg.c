int mainj (int i, char *pfname, int *ifs, char* jobname);
#include "stdlib.h"
#include <stdio.h>
#include <string.h>
void do_chunk_(int *i,int *one, int *idir);
void do_chunkx_(int *i,int *one, int *idir);
void dealloc_chunk_(int *i);
void radarretsub2_(int *nmu2,  int *nmfreq2,   int *icL, float *tbRgrid,  
		  float *dprrain, int *ichunk, int *orbNumb, int *ialg, int *idir);
void radarretsub3_(int *nmu2,  int *nmfreq2,   int *icL, float *tbRgrid,  
		  float *dprrain, int *ichunk, int *orbNumb, int *ialg, int *idir);
void radarretsub4_(int *nmu2,  int *nmfreq2,   int *icL, float *tbRgrid,  
		  float *dprrain, int *ichunk, int *orbNumb, int *ialg, int *idir);
void radarretsub4_fs_(int *nmu2,  int *nmfreq2,   int *icL, float *tbRgrid,  
		  float *dprrain, int *ichunk, int *orbNumb, int *ialg, int *idir);
void dealloc_struct_(int *i);
void close_files_(int *i);
int main(int argc, char *argv[])
{
  char  fname[100];
  char jobname[255];
  int ifs;
  if(argc != 3)
    {fprintf(stderr, 
	     "\nCommand Line ERROR-should be 2 arguments (jobname, parameterFile)\n");
      exit(1);}
  
  strcpy(jobname, argv[1]);
  
  strcpy(&fname[0],argv[2]);

  int ndpr=mainj(1,fname,&ifs,&jobname[0]);
  int ny=49;
  int nx=300;
  int nz=88;
  printf("%i %i\n",ndpr, ifs);
  //exit(0);
  int i,one=1;
  int nmu=5, nmfreq=8, orbNumb=0;
  float *dprrain, *tbRgrid;
  tbRgrid=(float*) malloc(sizeof(float)*9300*49*14);
  dprrain=(float*) malloc(sizeof(float)*49*300);
  int ialg=1;
  int idir;
  int icL;
  int nchunk=ndpr/300;
  nchunk=6;
  for(i=0;i<=nchunk;i++)
    {
      if(ifs==1)
	do_chunkx_(&i,&one,&idir);
      else
	do_chunk_(&i,&one,&idir);
      icL=i*300;
      if(i>=6)
	{
	  radarretsub2_(&nmu,  &nmfreq,   &icL, tbRgrid,  
			dprrain, &i, &orbNumb, &ialg, &idir);
	  radarretsub3_(&nmu,  &nmfreq,   &icL, tbRgrid,  
			dprrain, &i, &orbNumb, &ialg, &idir);
	  if(ifs==1)
	    radarretsub4_fs_(&nmu,  &nmfreq,   &icL, tbRgrid,  
			     dprrain, &i, &orbNumb, &ialg, &idir);
	  else
	    radarretsub4_(&nmu,  &nmfreq,   &icL, tbRgrid,  
			     dprrain, &i, &orbNumb, &ialg, &idir);
	  
	  dealloc_struct_(&i);
	}
      
      dealloc_chunk_(&i);
    }
closefiles_(&one);
}

