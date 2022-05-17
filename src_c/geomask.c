#include <stdlib.h>
#include <stdio.h>


//[2160][1080]
//[144][91][12]
void read_geodat_(char *lsflag, short int *sstdata)
{
  FILE *fp;
  int    i, j, k, ic;
  int b;

  fp=fopen("AncData/geomask_petty.asc","r");
  ic=0;
  for (j=0; j<1080; j++)
	for(i=0; i<2160; i++)
	  {
		fscanf(fp,"%2i ",&b);
		lsflag[ic]=b;
		ic++;
	  }
  
  fclose(fp);
  fp=fopen("AncData/sst.dat","r");
  ic=0;
  
    
  for(k=0;k<12;k++)
    for (j=0; j<144; j++)
      for(i=0; i<91; i++)
	{
	  fscanf(fp,"%i ",&b);
	  sstdata[ic]=b;
	  ic++;
	  //	  printf("%i ",sstdata[j][i][k]);
	}
  
  fclose(fp);

}

//[2160][1080]
//[144][91][12]
int igetlandsea_(float *rlat1, float *rlon1,char *lsflag) 
{
  int i,j, iflag, igetlandsea;

  float rlat;
  float rlon;

  rlat=*rlat1+0.;
  rlon=*rlon1+0.;

  if (rlon < 0.0) rlon = rlon + 360.0;
  if (rlon > 360.0) rlon = rlon - 360.0;
  i = rlon*6.0  + 1;
  j = (rlat + 90.0)*6.0 + 1;
  //  printf("%g %g %g %g %i %i \n",rlon,rlat,*rlon1,*rlat1, i, j);

  iflag =  lsflag[i+2160*j];
  
  if (iflag == 5 || iflag == 4)
	igetlandsea = 0;
  else
	{
	  if (iflag == 0) 
		igetlandsea = 1;
	  else
		{
		  if(iflag == 0 || iflag == 8)
			igetlandsea = 2;
		  else
			igetlandsea = 3;
		  
		}
	}
  return igetlandsea;
  
}



/* ! ocean
   igetlandsea = 0
   else if (iflag .eq. 0) then                ! land
   igetlandsea = 1
   else if (iflag .eq. 6 .or. iflag .eq. 8) then   ! coast
   igetlandsea = 2
*/
