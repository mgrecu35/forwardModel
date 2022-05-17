#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define iiad local_rd_var_mp_iiad_

int inputParser(char *p2arg, char f1CGMI[3][1000], char f2AKu[1000],
		char f2ADPR[1000], char f2KuENV[1000], char fSNOW[1000], 
		char fSEAICE[1000], char fCMB[1000], int *rseed1, int *rseed2,
		int *ialg, int *ifs)
{
  FILE *fn;
  char str[1000], key[20] = "", vstr[1000], *rs;
  int i, cc = 0, cr = 0, ck = 0, cd = 0, ce = 0, csnow = 0, csice = 0, 
      cm = 0, rv;
  extern int iiad;
  *rseed1 = -1;
  *rseed2 = -1;
  strcpy(fSNOW, "");
  strcpy(fSEAICE, "");
  fn = fopen(p2arg, "r");
  iiad=0;
  if(fn == NULL) {fprintf(stderr, "\nERROR opening parameter file.\n"); exit(1);}
  while(1)
  {
    rs = fgets(str, 1000, fn);
    if(rs == NULL || strlen(str) < 5) break;
    for(i=0;i<strlen(str);++i)
      if(str[i] == '=') {str[i] = ' '; break;}
    rv = sscanf(str, "%s %s", key, vstr);
    if(strcmp(key, "f1CGMI") == 0)
    {
      if(cc == 3)
      {
        fprintf(stderr, "\nERROR in Parameter file - more than 3 1CGMI keywords.\n");
        fclose(fn);
        exit(1);
      }
      strcpy(f1CGMI[cc], vstr);
      ++cc;
    }
    else if(strcmp(key, "f2AKu") == 0)
      {strcpy(f2AKu, vstr); ++ck;}
    else if(strcmp(key, "f2ADPR") == 0)
      {strcpy(f2ADPR, vstr); ++cd;}
    else if(strcmp(key, "f2KuENV") == 0)
      {strcpy(f2KuENV, vstr); ++ce;}
    else if(strcmp(key, "fSNOW") == 0)
      continue;
    else if(strcmp(key, "fSEAICE") == 0)
      continue;
    else if(strcmp(key, "fCMB") == 0)
      {strcpy(fCMB, vstr); ++cm;}
    else if(strcmp(key, "rseed") == 0)
    {
      if(cr == 0) *rseed1 = atoi(vstr);
      if(cr == 1) *rseed2 = atoi(vstr);
      if(cr == 2)
      {
        fprintf(stderr, "\nERROR in parameter file - more than 2 rseeds.\n");
        fclose(fn);
        exit(1);
      }
      ++cr;
    }
    else if(strcmp(key,"ialg")==0)
      { *ialg=atoi(vstr);}
    else if(strcmp(key,"ifs")==0)
      { *ifs=atoi(vstr);}
    else if(strcmp(key,"iiad")==0)
      { iiad=atoi(vstr);}
    else
    {
      fprintf(stderr, "\nError in parameter file - invalid keyword %s.\n", key);
      fprintf(stderr,"  %s\n",key);
      fclose(fn);
      exit(1);
    }
  }
  fclose(fn);
  if(cc < 3) {fprintf(stderr, "\nERROR in parameter file - less than 3 f1CGMI keywords.\n"); exit(1);}
  if(cr > 2) {fprintf(stderr, "\nERROR in parameter file - more than 2 rseeds.\n"); exit(1);}
  if(ck != 1) {fprintf(stderr, "\nERROR in parameter file - wrong number of f2AKu keywords.\n"); exit(1);}
  if(cd != 1) {fprintf(stderr, "\nERROR in parameter file - wrong number of f2ADPR keywords.\n"); exit(1);}
  if(ce != 1) {fprintf(stderr, "\nERROR in parameter file - wrong number of f2KuENV keywords.\n"); exit(1);}
  //if(csnow != 1) {fprintf(stderr, "\nERROR in parameter file - wrong number of fSNOW keywords.\n"); exit(1);}
  //if(csice != 1) {fprintf(stderr, "\nERROR in parameter file - wrong number of fSEAICE keywords.\n"); exit(1);}
  if(cm != 1) {fprintf(stderr, "\nERROR in parameter file - wrong number of fCMB keywords.\n"); exit(1);}
  return 0;
}
