//  SFM  04/06/2013  Added string.h
//  SFM  05/06/2013  Modifications by LW to allow passing in jobname
//
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//  SFM  04/06/2013  Changed file name lengths to 1000
int inputParser(char *p2arg, char f1ctmi[3][1000], 
		char f2aKu[1000], char f2aDPR[1000], char f2AKuENV[1000], 
		char fSNOW[1000], char fSEAICE[1000], char f2CMB[1000], 
		int *rseed1,int *rseed2, int *ialg, int *ifs);
void mainfort_(char jobname[255], char f1ctmi1[1000],char f1ctmi2[1000],
               char f1ctmi3[1000],
	       char f2AKu[1000], char f2aDPR[1000], char f2AKuENV[1000], 
	       char fSNOW[1000], char fSEAICE[1000], char f2CMB[1000], 
	       int *rseed1, int *rseed2,
	       int *igmi1, int *igmi2, int *igm3, int *i2AKu, int *i2AKuENV, 
	       int *i2aDPR, int *iSNOW, int *iSEAICE, int *i2CMB, int *ialg, int *ndpr, int *ifs);

int main2 (int argc, char *argv[])
{
  char f1ctmi[3][1000],  f2AKu[1000], f2ADPR[1000], f2AKuENV[1000] ;
  char f2CMB[1000], fSNOW[1000], fSEAICE[1000] ;
  char **fname;
  char jobname[255], args[500], *p2arg;
  int igmi1, igmi2, igmi3, i2AKu, i2ADPR, i2AKuENV, iSNOW, iSEAICE, i2CMB;
  int rseed1, rseed2, ialg, ifs;

  if(argc != 3)
    {fprintf(stderr, 
     "\nCommand Line ERROR-should be 2 arguments (jobname, parameterFile)\n");
     exit(1);}
  strcpy(jobname, argv[1]);
  inputParser(argv[2],f1ctmi,f2AKu,f2ADPR,f2AKuENV,fSNOW,fSEAICE, 
              f2CMB,&rseed1,&rseed2,&ialg,&ifs);

  igmi1=strlen(&f1ctmi[0][0]);
  igmi2=strlen(&f1ctmi[1][0]);
  igmi3=strlen(&f1ctmi[2][0]);
  printf(" File name lengths, igmi1,igmi2,igmi3 : %i  %i  %i  \n",
          igmi1,igmi2,igmi3) ;
  
  i2AKu=strlen(f2AKu);
  i2ADPR=strlen(f2ADPR);
  i2AKuENV=strlen(f2AKuENV);
  iSNOW=strlen(fSNOW);
  iSEAICE=strlen(fSEAICE);
  i2CMB=strlen(f2CMB);
//begin  WSO 12/28/16 suppress output file
//  FILE *fout;
//  fout=fopen("runStatus.txt","w");
//  fprintf(fout,"%i\n",-1);
//  fclose(fout);
//end    WSO 12/28/16
  printf(" File name lengths, i2AKu,i2ADPR,i2AKuENV : %i %i %i \n",
          i2AKu,i2ADPR,i2AKuENV,i2CMB) ;
  printf(" File name lengths, iSNOW,iSEAICE,i2CMB : %i %i %i \n",
          iSNOW,iSEAICE,i2CMB) ;

  /*  mainfort_(jobname, &f1ctmi[0][0],&f1ctmi[1][0],&f1ctmi[2][0], f2AKu, 
            f2ADPR, f2AKuENV, fSNOW, fSEAICE, f2CMB, &rseed1, &rseed2, 
	    &igmi1, &igmi2, &igmi3, &i2AKu, &i2ADPR, &i2AKuENV, &iSNOW, 
	    &iSEAICE, &i2CMB, &ialg);
  */
  printf(" COMPLETION \n") ;

//begin  WSO 12/28/16 suppress output file
//  fout=fopen("runStatus.txt","w");
//  fprintf(fout,"%i\n",0);
//  fclose(fout);
//end    WSO 12/28/16

  exit(0);
}

int mainj (int i, char *pfname, int *ifs, char *jobname)
{
  char f1ctmi[3][1000],  f2AKu[1000], f2ADPR[1000], f2AKuENV[1000] ;
  char f2CMB[1000], fSNOW[1000], fSEAICE[1000] ;
  char **fname, pfname2[1000];
  char args[500], *p2arg;
  int igmi1, igmi2, igmi3, i2AKu, i2ADPR, i2AKuENV, iSNOW, iSEAICE, i2CMB;
  int rseed1, rseed2, ialg;

  printf("%s \n",pfname);
  strcpy(&pfname2[0],pfname);
  //  exit(0);
  //strcpy(jobname, "junk");
  // strcpy(jobname, argv[1]);
  inputParser(&pfname2[0],f1ctmi,f2AKu,f2ADPR,f2AKuENV,fSNOW,fSEAICE, 
              f2CMB,&rseed1,&rseed2,&ialg,ifs);
 
  igmi1=strlen(&f1ctmi[0][0]);
  igmi2=strlen(&f1ctmi[1][0]);
  igmi3=strlen(&f1ctmi[2][0]);
  printf(" File name lengths, igmi1,igmi2,igmi3 : %i  %i  %i  \n",
          igmi1,igmi2,igmi3) ;
  
  i2AKu=strlen(f2AKu);
  i2ADPR=strlen(f2ADPR);
  i2AKuENV=strlen(f2AKuENV);
  iSNOW=strlen(fSNOW);
  iSEAICE=strlen(fSEAICE);
  i2CMB=strlen(f2CMB);
//begin  WSO 12/28/16 suppress output file
//  FILE *fout;
//  fout=fopen("runStatus.txt","w");
//  fprintf(fout,"%i\n",-1);
//  fclose(fout);
//end    WSO 12/28/16
  printf(" File name lengths, i2AKu,i2ADPR,i2AKuENV : %i %i %i \n",
          i2AKu,i2ADPR,i2AKuENV,i2CMB) ;
  printf(" File name lengths, iSNOW,iSEAICE,i2CMB : %i %i %i \n",
          iSNOW,iSEAICE,i2CMB) ;
  int ndpr;
  printf("%s \n",f2ADPR);
  printf("%s \n",f2AKu);

  mainfort_(jobname, &f1ctmi[0][0],&f1ctmi[1][0],&f1ctmi[2][0], f2AKu, 
            f2ADPR, f2AKuENV, fSNOW, fSEAICE, f2CMB, &rseed1, &rseed2, 
	    &igmi1, &igmi2, &igmi3, &i2AKu, &i2ADPR, &i2AKuENV, &iSNOW, 
	    &iSEAICE, &i2CMB, &ialg, &ndpr, ifs);
  
  printf(" COMPLETION %i\n",*ifs) ;
  
//begin  WSO 12/28/16 suppress output file
//  fout=fopen("runStatus.txt","w");
//  fprintf(fout,"%i\n",0);
//  fclose(fout);
//end    WSO 12/28/16

  return ndpr;  
}

void printjobname_(char* jobname)
{

  printf("%s \n",jobname);
}
