//  SFM 04/06/2013  Code module added in merge from M.Grecu's code
//  SFM  06/27/2013  Parameter name changes from W.Olson; reduce unused code
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf.h>
#include <mfhdf.h>
#include "TKheaders.h"
#include "TK_1CGMI.h"
//begin WSO 04/07/2013
#include <math.h>
//end WSO 04/07/2013
#ifdef GFOR 
extern int __nbinmod_MOD_nbin;
#define nbins __nbinmod_MOD_nbin
#endif

#ifdef IFORT 
extern int nbinmod_mp_nbin_;
#define nbins nbinmod_mp_nbin_
#endif

//  begin  SFM  09/12/2013
extern TKINFO ctkfile;
//  end  SFM  09/12/2013
//  begin  SFM  09/04/2013
L2AKuENV_NS        L2AKuENVData;
L2AKuENVX_FS        L2AKuENVDatax;
//L2AKuENV_NS_aa1        L2AKuENVData;
//  end  SFM  09/04/2013

TKINFO       granuleHandle2AKuENV;

//  SFM  begin 12/09/2013; to pass out TK status message
//following line  WSO 9/1/13 add skin temperature to subroutine call

//LW 05/04/18 - add alg
int readenv_(char *jobname, int *ialg, char *f2akuENV, int *n1c21, int *ic, 
             int *nBSize, float *envQv, float *envTemp,
	     float *envPress, float *envSfcWind, float *envSfcWindU, float *envSfcWindV, float * envSknTemp, 
	     float *envSfcTemp, float *envSfcPress, float *envCloud)
{
  int             status_alpha ;
//  SFM  end 12/09/2013
  int status,i;
  char emptygranuletext[100];

//  SFM  begin  12/13/2013; to pass out TK status message
  char envfname[1000];
  strcpy(envfname,f2akuENV);
//  SFM  begin  12/13/2013

  status_alpha = 0 ;
  if(*ic==0)
    {
//  SFM  begin  12/09/2013; to pass out TK status message
//  begin  SFM  09/04/2013
//      status = TKopen(&f2akuENV[0], "2AKuENV_aa1", TKREAD, "HDF5", 
//                      jobname, &granuleHandle2AKuENV, 1);

      //LW 05/04/18 - add alg, 2APRENV
      if(*ialg == 2)
        status_alpha = TKopen(&f2akuENV[0], "2APRENV", TKREAD, "HDF5",  jobname, &granuleHandle2AKuENV, 1);
      else
        status_alpha = TKopen(&f2akuENV[0], "2AKuENV", TKREAD, "HDF5",  jobname, &granuleHandle2AKuENV, 1);

    if (status_alpha != 0)
      {
       printf("WARNING: Unable to access 2AKuENV %i \n",status_alpha);
       return status_alpha ;
      }
//  end  SFM  09/04/2013
//  SFM  end  12/09/2013
//  begin  SFM  09/12/2013
    TKtransferMetaData(&granuleHandle2AKuENV, &ctkfile);
    printf("  STATUS, meta transfer 2AKuENV %i  \n",status_alpha);
//  end  SFM  09/12/2013
    }

//  SFM  begin   12/13/2013
  if(strncmp(envfname,"nil",3) == 0) status_alpha = -1 ;
//  SFM  end     12/13/2013

  TKgetMetaString(&granuleHandle2AKuENV, "FileHeader", "EmptyGranule", 
                  emptygranuletext);	  
  if (strncmp(emptygranuletext,"EMPTY",5) == 0)
     {
      status_alpha = -3;
      printf("  STATUS, 2AKuENV has empty granule \n",status_alpha);
      return status_alpha ;
      }
      
  float ucomp, vcomp, speed; //WSO 04/07/2013

  int i3d=0;
  int i2d=0;
  int j,k;
//begin  WSO 9/19/13
  float ray_angle;
  float grav = 9.8;
  float rd = 287.;
  float pi = 3.14159265;
//end    WSO 9/19/13
  i=0;

  while((TKendOfFile (&granuleHandle2AKuENV) != TK_EOF) && i<*nBSize)
    {

      status=TKreadScan(&granuleHandle2AKuENV,&L2AKuENVData);

      for(j=0;j<49;j++)
	{
//begin  WSO 9/19/13 assign approximate ray angle
     ray_angle = abs(24 - j) * 17.0 / 24.;
//end    WSO 9/19/13
	  for(k=0;k<nbins;k++)
	    {
	      //printf("%i %i %i \n",j,k,i3d);
//begin  WSO 9/19/13 interpolate all upper-level parameters
//begin WSO 04/07/2013
//...resample air temperature profile
//
//TEST********
//  L2AKuENVData.VERENV.airTemperature[j][2*k] = -9999.9;
//  L2AKuENVData.VERENV.airPressure[j][2*k] = -9999.9;
//  L2AKuENVData.VERENV.waterVapor[j][2*k][0] = -9999.9;
//  L2AKuENVData.VERENV.waterVapor[j][2*k+1][0] = -9999.9;
//  L2AKuENVData.VERENV.cloudLiquidWater[j][2*k][0] = -9999.9;
//  L2AKuENVData.VERENV.cloudLiquidWater[j][2*k+1][0] = -9999.9;
//END TEST****

         if(L2AKuENVData.VERENV.airTemperature[j][2*k] > 0. && L2AKuENVData.VERENV.airTemperature[j][2*k+1] > 0.)
	        envTemp[i3d]= 0.5 *
    	      	(L2AKuENVData.VERENV.airTemperature[j][2*k] + L2AKuENVData.VERENV.airTemperature[j][2*k+1]);
         else if(L2AKuENVData.VERENV.airTemperature[j][2*k] > 0.)
           envTemp[i3d] = L2AKuENVData.VERENV.airTemperature[j][2*k];
         else if(L2AKuENVData.VERENV.airTemperature[j][2*k+1] > 0.)
            envTemp[i3d] = L2AKuENVData.VERENV.airTemperature[j][2*k+1];
         else
            envTemp[i3d] = L2AKuENVData.VERENV.airTemperature[j][2*k];
//assume that variable contains flag value as final option

         if(L2AKuENVData.VERENV.airPressure[j][2*k] > 0. && L2AKuENVData.VERENV.airPressure[j][2*k+1] > 0.)
           envPress[i3d] = L2AKuENVData.VERENV.airPressure[j][2*k] * 
            pow((L2AKuENVData.VERENV.airPressure[j][2*k+1]/L2AKuENVData.VERENV.airPressure[j][2*k]), 0.5);
         else if(L2AKuENVData.VERENV.airPressure[j][2*k] > 0. && envTemp[i3d] > 0.)
           envPress[i3d] = L2AKuENVData.VERENV.airPressure[j][2*k] *
            exp(grav * (125. * cos(ray_angle * pi / 180.) / 2.) / (rd * envTemp[i3d]));
         else if(L2AKuENVData.VERENV.airPressure[j][2*k+1] > 0. && envTemp[i3d] > 0.)
           envPress[i3d] = L2AKuENVData.VERENV.airPressure[j][2*k+1] * 
            exp(-grav * (125. * cos(ray_angle * pi / 180.) / 2.) / (rd * envTemp[i3d]));
         else
           envPress[i3d] = L2AKuENVData.VERENV.airPressure[j][2*k];
//assume that variable contains flag value as final option

//begin WSO 04/07/2013
//...for proper units of [g/m3]
//   note here, Mircea opts to use the water vapor adjusted by the algorithm
//   in rain areas [0]
         if(L2AKuENVData.VERENV.waterVapor[j][2*k][0] >= 0. && L2AKuENVData.VERENV.waterVapor[j][2*k+1][0] >= 0.)
	        envQv[i3d] = 0.5 * (L2AKuENVData.VERENV.waterVapor[j][2*k][0] + L2AKuENVData.VERENV.waterVapor[j][2*k+1][0]) * 1.e+3;
         else if (L2AKuENVData.VERENV.waterVapor[j][2*k][0] >= 0.)
           envQv[i3d] = L2AKuENVData.VERENV.waterVapor[j][2*k][0] * 1.e+3;
         else if (L2AKuENVData.VERENV.waterVapor[j][2*k+1][0] >= 0.)
           envQv[i3d] = L2AKuENVData.VERENV.waterVapor[j][2*k+1][0] * 1.e+3;
         else
           envQv[i3d] = L2AKuENVData.VERENV.waterVapor[j][2*k][0];
//assume that variable contains flag value as final option

         if(L2AKuENVData.VERENV.cloudLiquidWater[j][2*k][0] >= 0. && L2AKuENVData.VERENV.cloudLiquidWater[j][2*k+1][0] >= 0.)
  	        envCloud[i3d] = 0.5 * (L2AKuENVData.VERENV.cloudLiquidWater[j][2*k][0] + 
                 L2AKuENVData.VERENV.cloudLiquidWater[j][2*k+1][0]) * 1.e+3;
         else if(L2AKuENVData.VERENV.cloudLiquidWater[j][2*k][0] >= 0.)
           envCloud[i3d] = L2AKuENVData.VERENV.cloudLiquidWater[j][2*k][0] * 1.e+3;
         else if(L2AKuENVData.VERENV.cloudLiquidWater[j][2*k+1][0] >= 0.)
           envCloud[i3d] = L2AKuENVData.VERENV.cloudLiquidWater[j][2*k+1][0] * 1.e+3;
         else
           envCloud[i3d] = L2AKuENVData.VERENV.cloudLiquidWater[j][2*k][0];
//assume that variable contains flag value as final option

//diagnostic
//        if(j == 24)
//          printf("j: %5i, angle:  %8.4f,  k: %5i,  envTemp: %10.2f  %10.2f %10.2f,  envPress: %12.2f %12.2f %12.2f,  envQv: %12.4f %12.4f %12.4f,  envCloud: %12.4f %12.4f %12.4f\n",
//           j, ray_angle, k, L2AKuENVData.VERENV.airTemperature[j][2*k], L2AKuENVData.VERENV.airTemperature[j][2*k+1], envTemp[i3d], 
//           L2AKuENVData.VERENV.airPressure[j][2*k], L2AKuENVData.VERENV.airPressure[j][2*k+1], envPress[i3d], 
//           L2AKuENVData.VERENV.waterVapor[j][2*k][0] * 1.e+3, L2AKuENVData.VERENV.waterVapor[j][2*k+1][0] *1.e+3, envQv[i3d], 
//           L2AKuENVData.VERENV.cloudLiquidWater[j][2*k][0] * 1.e+3, L2AKuENVData.VERENV.cloudLiquidWater[j][2*k+1][0] *1.e+3, envCloud[i3d]);
//end diagnostic

	      i3d++;
//end WSO 04/07/2013
//end    WSO 9/19/13
	      }
//  07/16/2013  SFM  New parameter name
//begin WSO 9/1/13 add skin temperature
     envSknTemp[i2d]=L2AKuENVData.VERENV.skinTemperature[j];
//end   WSO 9/1/13
//aa	  envSfcTemp[i2d]=L2AKuENVData.VERENV.groundTemperature[j];
	  envSfcTemp[i2d]=L2AKuENVData.VERENV.surfaceTemperature[j];

//begin  WSO 9/9/13 surface pressure in hundredths of hPa
     envSfcPress[i2d]=L2AKuENVData.VERENV.surfacePressure[j];
//end    WSO 9/9/13

//begin WSO 04/07/2013
//...sum u and v components of surface wind
	  envSfcWind[i2d]=L2AKuENVData.VERENV.surfaceWind[j][0];
	  //envSfcWind[i2d]=0.;
          ucomp = L2AKuENVData.VERENV.surfaceWind[j][0];
          envSfcWindU[i2d]= L2AKuENVData.VERENV.surfaceWind[j][0];
          vcomp = L2AKuENVData.VERENV.surfaceWind[j][1];
          envSfcWindV[i2d]= L2AKuENVData.VERENV.surfaceWind[j][1];
          speed = ucomp * ucomp + vcomp * vcomp;
          if(speed > 0.) envSfcWind[i2d] = sqrt(speed);
//end WSO 04/07/2013
	  i2d++;
	}

//   exit(0);
      i++;
    }

//  SFM  begin 12/09/2013; to pass out TK status message
  return status_alpha ;
//  SFM  end   12/09/2013
}


int readenvx_(char *jobname, int *ialg, char *f2akuENV, int *n1c21, int *ic, 
             int *nBSize, float *envQv, float *envTemp,
	     float *envPress, float *envSfcWind, float *envSfcWindU, float *envSfcWindV, float * envSknTemp, 
	     float *envSfcTemp, float *envSfcPress, float *envCloud)
{
  int             status_alpha ;
//  SFM  end 12/09/2013
  int status,i;
  char emptygranuletext[100];

//  SFM  begin  12/13/2013; to pass out TK status message
  char envfname[1000];
  strcpy(envfname,f2akuENV);
//  SFM  begin  12/13/2013

  status_alpha = 0 ;
  if(*ic==0)
    {
//  SFM  begin  12/09/2013; to pass out TK status message
//  begin  SFM  09/04/2013
//      status = TKopen(&f2akuENV[0], "2AKuENV_aa1", TKREAD, "HDF5", 
//                      jobname, &granuleHandle2AKuENV, 1);

      //LW 05/04/18 - add alg, 2APRENV
      if(*ialg == 2)
        status_alpha = TKopen(&f2akuENV[0], "2APRENV", TKREAD, "HDF5",  jobname, &granuleHandle2AKuENV, 1);
      else
        status_alpha = TKopen(&f2akuENV[0], "2AKuENVX", TKREAD, "HDF5",  jobname, &granuleHandle2AKuENV, 1);

    if (status_alpha != 0)
      {
       printf("WARNING: Unable to access 2AKuENV %i \n",status_alpha);
       return status_alpha ;
      }
//  end  SFM  09/04/2013
//  SFM  end  12/09/2013
//  begin  SFM  09/12/2013
    TKtransferMetaData(&granuleHandle2AKuENV, &ctkfile);
    printf("  STATUS, meta transfer 2AKuENV %i  \n",status_alpha);
//  end  SFM  09/12/2013
    }

//  SFM  begin   12/13/2013
  if(strncmp(envfname,"nil",3) == 0) status_alpha = -1 ;
//  SFM  end     12/13/2013

  TKgetMetaString(&granuleHandle2AKuENV, "FileHeader", "EmptyGranule", 
                  emptygranuletext);	  
  if (strncmp(emptygranuletext,"EMPTY",5) == 0)
     {
      status_alpha = -3;
      printf("  STATUS, 2AKuENV has empty granule \n",status_alpha);
      return status_alpha ;
      }
      
  float ucomp, vcomp, speed; //WSO 04/07/2013

  int i3d=0;
  int i2d=0;
  int j,k;
//begin  WSO 9/19/13
  float ray_angle;
  float grav = 9.8;
  float rd = 287.;
  float pi = 3.14159265;
//end    WSO 9/19/13
  i=0;

  while((TKendOfFile (&granuleHandle2AKuENV) != TK_EOF) && i<*nBSize)
    {

      status=TKreadScan(&granuleHandle2AKuENV,&L2AKuENVDatax);

      for(j=0;j<49;j++)
	{
//begin  WSO 9/19/13 assign approximate ray angle
     ray_angle = abs(24 - j) * 17.0 / 24.;
//end    WSO 9/19/13
	  for(k=0;k<nbins;k++)
	    {
	      //printf("%i %i %i \n",j,k,i3d);
//begin  WSO 9/19/13 interpolate all upper-level parameters
//begin WSO 04/07/2013
//...resample air temperature profile
//
//TEST********
//  L2AKuENVData.VERENV.airTemperature[j][2*k] = -9999.9;
//  L2AKuENVData.VERENV.airPressure[j][2*k] = -9999.9;
//  L2AKuENVData.VERENV.waterVapor[j][2*k][0] = -9999.9;
//  L2AKuENVData.VERENV.waterVapor[j][2*k+1][0] = -9999.9;
//  L2AKuENVData.VERENV.cloudLiquidWater[j][2*k][0] = -9999.9;
//  L2AKuENVData.VERENV.cloudLiquidWater[j][2*k+1][0] = -9999.9;
//END TEST****

/*
  if(L2AKuENVDatax.VERENV.airTemperature[j][2*k] > 0. && L2AKuENVDatax.VERENV.airTemperature[j][2*k+1] > 0.)
  envTemp[i3d]= 0.5 *
  (L2AKuENVDatax.VERENV.airTemperature[j][2*k] + L2AKuENVDatax.VERENV.airTemperature[j][2*k+1]);
  else if(L2AKuENVDatax.VERENV.airTemperature[j][2*k] > 0.)
  envTemp[i3d] = L2AKuENVDatax.VERENV.airTemperature[j][2*k];
  else if(L2AKuENVDatax.VERENV.airTemperature[j][2*k+1] > 0.)
  envTemp[i3d] = L2AKuENVDatax.VERENV.airTemperature[j][2*k+1];
  else
  envTemp[i3d] = L2AKuENVDatax.VERENV.airTemperature[j][2*k];*/
//assume that variable contains flag value as final option
	 envSfcTemp[i2d]=L2AKuENVDatax.VERENV.surfaceTemperature[j];
         if(L2AKuENVDatax.VERENV.airPressure[j][2*k] > 0. && L2AKuENVDatax.VERENV.airPressure[j][2*k+1] > 0.)
           envPress[i3d] = L2AKuENVDatax.VERENV.airPressure[j][2*k] * 
            pow((L2AKuENVDatax.VERENV.airPressure[j][2*k+1]/L2AKuENVDatax.VERENV.airPressure[j][2*k]), 0.5);
         else if(L2AKuENVDatax.VERENV.airPressure[j][2*k] > 0. && envTemp[i3d] > 0.)
           envPress[i3d] = L2AKuENVDatax.VERENV.airPressure[j][2*k] *
            exp(grav * (125. * cos(ray_angle * pi / 180.) / 2.) / (rd * envTemp[i3d]));
         else if(L2AKuENVDatax.VERENV.airPressure[j][2*k+1] > 0. && envTemp[i3d] > 0.)
           envPress[i3d] = L2AKuENVDatax.VERENV.airPressure[j][2*k+1] * 
            exp(-grav * (125. * cos(ray_angle * pi / 180.) / 2.) / (rd * envTemp[i3d]));
         else
           envPress[i3d] = L2AKuENVDatax.VERENV.airPressure[j][2*k];
//assume that variable contains flag value as final option

//begin WSO 04/07/2013
//...for proper units of [g/m3]
//   note here, Mircea opts to use the water vapor adjusted by the algorithm
//   in rain areas [0]
         if(L2AKuENVDatax.VERENV.waterVapor[j][2*k][0] >= 0. && L2AKuENVDatax.VERENV.waterVapor[j][2*k+1][0] >= 0.)
	        envQv[i3d] = 0.5 * (L2AKuENVDatax.VERENV.waterVapor[j][2*k][0] + L2AKuENVDatax.VERENV.waterVapor[j][2*k+1][0]) * 1.e+3;
         else if (L2AKuENVDatax.VERENV.waterVapor[j][2*k][0] >= 0.)
           envQv[i3d] = L2AKuENVDatax.VERENV.waterVapor[j][2*k][0] * 1.e+3;
         else if (L2AKuENVDatax.VERENV.waterVapor[j][2*k+1][0] >= 0.)
           envQv[i3d] = L2AKuENVDatax.VERENV.waterVapor[j][2*k+1][0] * 1.e+3;
         else
           envQv[i3d] = L2AKuENVDatax.VERENV.waterVapor[j][2*k][0];
//assume that variable contains flag value as final option

         if(L2AKuENVDatax.VERENV.cloudLiquidWater[j][2*k][0] >= 0. && L2AKuENVDatax.VERENV.cloudLiquidWater[j][2*k+1][0] >= 0.)
  	        envCloud[i3d] = 0.5 * (L2AKuENVDatax.VERENV.cloudLiquidWater[j][2*k][0] + 
                 L2AKuENVDatax.VERENV.cloudLiquidWater[j][2*k+1][0]) * 1.e+3;
         else if(L2AKuENVDatax.VERENV.cloudLiquidWater[j][2*k][0] >= 0.)
           envCloud[i3d] = L2AKuENVDatax.VERENV.cloudLiquidWater[j][2*k][0] * 1.e+3;
         else if(L2AKuENVDatax.VERENV.cloudLiquidWater[j][2*k+1][0] >= 0.)
           envCloud[i3d] = L2AKuENVDatax.VERENV.cloudLiquidWater[j][2*k+1][0] * 1.e+3;
         else
           envCloud[i3d] = L2AKuENVDatax.VERENV.cloudLiquidWater[j][2*k][0];
//assume that variable contains flag value as final option

//diagnostic
//        if(j == 24)
//          printf("j: %5i, angle:  %8.4f,  k: %5i,  envTemp: %10.2f  %10.2f %10.2f,  envPress: %12.2f %12.2f %12.2f,  envQv: %12.4f %12.4f %12.4f,  envCloud: %12.4f %12.4f %12.4f\n",
//           j, ray_angle, k, L2AKuENVDatax.VERENV.airTemperature[j][2*k], L2AKuENVDatax.VERENV.airTemperature[j][2*k+1], envTemp[i3d], 
//           L2AKuENVDatax.VERENV.airPressure[j][2*k], L2AKuENVDatax.VERENV.airPressure[j][2*k+1], envPress[i3d], 
//           L2AKuENVDatax.VERENV.waterVapor[j][2*k][0] * 1.e+3, L2AKuENVDatax.VERENV.waterVapor[j][2*k+1][0] *1.e+3, envQv[i3d], 
//           L2AKuENVDatax.VERENV.cloudLiquidWater[j][2*k][0] * 1.e+3, L2AKuENVDatax.VERENV.cloudLiquidWater[j][2*k+1][0] *1.e+3, envCloud[i3d]);
//end diagnostic

	      i3d++;
//end WSO 04/07/2013
//end    WSO 9/19/13
	      }
//  07/16/2013  SFM  New parameter name
//begin WSO 9/1/13 add skin temperature
     envSknTemp[i2d]=L2AKuENVDatax.VERENV.skinTemperature[j];
//end   WSO 9/1/13
//aa	  envSfcTemp[i2d]=L2AKuENVDatax.VERENV.groundTemperature[j];
     
     //printf("%f ", envSfcTemp[i2d]);
//begin  WSO 9/9/13 surface pressure in hundredths of hPa
     envSfcPress[i2d]=L2AKuENVDatax.VERENV.surfacePressure[j];
//end    WSO 9/9/13

//begin WSO 04/07/2013
//...sum u and v components of surface wind
	  envSfcWind[i2d]=L2AKuENVDatax.VERENV.surfaceWind[j][0];
	  //envSfcWind[i2d]=0.;
          ucomp = L2AKuENVDatax.VERENV.surfaceWind[j][0];
          envSfcWindU[i2d]= L2AKuENVDatax.VERENV.surfaceWind[j][0];
          vcomp = L2AKuENVDatax.VERENV.surfaceWind[j][1];
          envSfcWindV[i2d]= L2AKuENVDatax.VERENV.surfaceWind[j][1];
          speed = ucomp * ucomp + vcomp * vcomp;
          if(speed > 0.) envSfcWind[i2d] = sqrt(speed);
//end WSO 04/07/2013
	  i2d++;
	}

//   exit(0);
      i++;
    }

//  SFM  begin 12/09/2013; to pass out TK status message
  return status_alpha ;
//  SFM  end   12/09/2013
}
