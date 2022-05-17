//  SFM 04/06/2013  Code module added in merge from M.Grecu's code
//  SFM 05/06/2013  Modifications from LW to facilitate using job names
//  SFM 06/27/2013  Parameter name changes from W.Olson; reduce unused code
//  SFM 07/19/2013  Large volume of code added for M.Grecu
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf.h>
#include <mfhdf.h>
#include "TKheaders.h"
#include "TK_2BCMBT_hdf5.h"
#ifdef GFOR 
extern int __nbinmod_MOD_imemb;
#define nbins __nbinmod_MOD_imemb
//begin  WSO 9/15/13 
extern float __missingmod_MOD_missing_r4;
#define missing_r4c __missingmod_MOD_missing_r4
extern short __missingmod_MOD_missing_i2;
#define missing_i2c __missingmod_MOD_missing_i2
extern long __missingmod_MOD_missing_i4;
#define missing_i4c __missingmod_MOD_missing_i4
extern int __nbinmod_MOD_ntransition;
#define ntransitions __nbinmod_MOD_ntransition
//end    WSO 9/15/13
#endif

#ifdef IFORT 
extern int nbinmod_mp_nbin_;
//begin  WSO 8/8/13
extern int nbinmod_mp_ntransition_;
//end    WSO 8/8/13
#define nbins nbinmod_mp_nbin_
//begin  WSO 8/8/13
#define ntransitions nbinmod_mp_ntransition_
//end    WSO 8/8/13
//begin  WSO 9/15/13
extern float missingmod_mp_missing_r4_;
#define missing_r4c missingmod_mp_missing_r4_
extern short missingmod_mp_missing_i2_;
#define missing_i2c missingmod_mp_missing_i2_
extern long  missingmod_mp_missing_i4_;
#define missing_i4c missingmod_mp_missing_i4_
//end    WSO 9/15/13
#endif

//begin WSO 04/07/2013
//Note that there were many structure/variable name changes in this
//version to be compatible with TKIO 3.50.8
//All S1 and S2 were changed to NS and MS, respectively
//The variable ending "Out" was removed because a separate Input structure
//was created
//end WSO 04/07/2013

extern TKINFO dprtkfile;
TKINFO ctkfile;
TKINFO ctkfileIn;

L2BCMBT_SWATHS swath1;
L2ADPR_SWATHS dprswathT;
L2AKu_NS      L2AKuData;

void openoutputfilet_(char *jobname, char *fname)
{
  char ranstring[80], ranMsg[80];
  int status;
  printf(" Output file = %s \n",fname);

  int ret;
//  SFM  04/06/2013  Changed file type to 2BCMB   WSO 04/07/2013
//  SFM  09/04/2013  Changed jobid to jobname
//  SFM  09/11/2013  Moved metadta settings to centralized location

  ret = TKopen(fname, "2BCMBT", TKWRITE, "HDF5", jobname, &ctkfile,1); //WSO 04/07/2013
}

void openinputfilet_(char *jobname, char *fname)
{
  int ret;

//  SFM  04/06/2013  Changed file type to 2BCMB   WSO 04/07/2013
  ret = TKopen(fname, "2BCMB", TKREAD, "HDF5", jobname, &ctkfileIn, 1);
  printf("%s %i \n",fname,ret);
}

void readscant_(void)
{
  int ret, i;
  ret= TKreadScan(&ctkfileIn,&swath1);

}

void closeoutputfilet_(void)
{
  int ret;
  ret=TKclose(&ctkfile);
}

void closedproutputfilet_(void)
{
  int ret;
  ret=TKclose(&dprtkfile);
}

void setlatlons1t_(float *lat, float *lon, float *sfcPrecip, 
                  float *sfcPrecipStd, float *piaOut)
{
  int i;
  extern L2BCMBT_SWATHS swath1;

  for(i=0;i<49;i++)
    {
      swath1.NS.Latitude[i]=lat[i];
      swath1.NS.Longitude[i]=lon[i];
      if(swath1.NS.Longitude[i]>180)
	swath1.NS.Longitude[i]-=360;
      swath1.NS.surfPrecipTotRate[i]=sfcPrecip[i];
      swath1.NS.surfPrecipTotRateSigma[i]=sfcPrecipStd[i];
      swath1.NS.pia[i]=piaOut[i];
    }
}

void setlatlons2t_(float *lat, float *lon, float *sfcPrecip, 
                  float *sfcPrecipStd, float *piaOutKu, float *piaOutKa)
{
  int i;
  //  extern L2BCMBT_SWATHS swath;

  for(i=12;i<37;i++)
    {
      //  swath.MS.Latitude[i-12]=lat[i];
      // swath.MS.Longitude[i-12]=lon[i];
      // if(swath.MS.Longitude[i-12]>180)
      //	swath.MS.Longitude[i-12]-=360;
      // swath.MS.surfPrecipTotRate[i-12]=sfcPrecip[i];
      // swath.MS.surfPrecipTotRateSigma[i-12]=sfcPrecipStd[i];
      //swath.MS.pia[i-12][0]=piaOutKu[i];
      //swath.MS.pia[i-12][1]=piaOutKa[i];
    }
}

void copyrrates1t_(float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;

  for(k=0;k<nbins;k++)
    {
      swath1.NS.precipTotRate[*i][k]=rrate[k];
      swath1.NS.precipTotRateSigma[*i][k]=rratestd[k];
    }
}

void copysflfractt_(float *lfract, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;

  swath1.NS.surfLiqRateFrac[*i]=*lfract;
  
}

void copyzkat_(float *zka, int *i)
{
  int k;
  extern L2ADPR_SWATHS dprswathT;

  for(k=0;k<nbins;k++)
    {
      dprswathT.MS.PRE.zFactorMeasured[*i][2*k]=(int)(zka[k]*100);
      dprswathT.MS.PRE.zFactorMeasured[*i][2*k+1]=(int)(zka[k]*100);
      //      printf("%g ",zka[k]);
    }
}

void copypiakat_(float *piaKa, int *i)
{
  int k;
  extern L2ADPR_SWATHS dprswathT;
 
  for(k=0;k<nbins;k++)
    {
      dprswathT.MS.SRT.pathAtten[*i]=*piaKa;
    }
}

void copytruerratet_(float *rrate, int *i)
{
  int k;
  extern L2ADPR_SWATHS dprswathT;

  for(k=0;k<nbins;k++)
    {
      dprswathT.NS.SLV.precipRate[*i][2*k]=(int)(rrate[k]*100);
      dprswathT.NS.SLV.precipRate[*i][2*k+1]=(int)(rrate[k]*100);
    }
}

//begin  WSO 8/30/13
void copyenvsfqvs1t_(float *envQv, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  swath1.NS.surfaceVaporDensity[*i]=envQv[nbins-1];
}

void copyenvsfqvs2t_(float *envQv, int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;
  //swath1.MS.surfaceVaporDensity[*i]=envQv[nbins-1];
}

void copyenvqvs1t_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  for(k=0;k<10;k++)
      swath1.NS.vaporDensity[*i][k]=envQv[envnodes[k]-1];
}

void copyenvqvs2t_(float *envQv, short *envnodes, int *i)
{
  int k;
  //extern L2BCMBT_SWATHS swath;
  //for(k=0;k<10;k++)
  //   swath1.MS.vaporDensity[*i][k]=envQv[envnodes[k]-1];
}

void copyenvpresss1t_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  for(k=0;k<10;k++)
    swath1.NS.airPressure[*i][k]=envQv[envnodes[k]-1];
}

void copyenvpresss2t_(float *envQv, short *envnodes, int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;
  //for(k=0;k<10;k++)
  //  swath1.MS.airPressure[*i][k]=envQv[envnodes[k]-1];
}

void copyenvtemps1t_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  for(k=0;k<10;k++)
     {
      swath1.NS.envParamNode[*i][k]=envnodes[k]-1;
      swath1.NS.airTemperature[*i][k]=envQv[envnodes[k]-1];
     }
}

void copyenvtemps2t_(float *envQv, short *envnodes, int *i)
{
  int k;
  // extern L2BCMBT_SWATHS swath;
  for(k=0;k<10;k++)
    {
      //   swath1.MS.envParamNode[*i][k]=envnodes[k]-1;
      //swath1.MS.airTemperature[*i][k]=envQv[envnodes[k]-1];
    }
}

void copyenvsftemps1t_(float *envQv, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  swath1.NS.surfaceAirTemperature[*i]=envQv[nbins-1];
}

void copyenvsftemps2t_(float *envQv, int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;
  // swath1.MS.surfaceAirTemperature[*i]=envQv[nbins-1];
}

//end    WSO 8/30/13

void copypwcs1t_(float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;

  for(k=0;k<nbins;k++)
    {
      swath1.NS.precipTotWaterCont[*i][k]=rrate[k];
      swath1.NS.precipTotWaterContSigma[*i][k]=rratestd[k];
    }
}

//begin  WSO 8/7/13
void copylwcfracs1t_(float *mlwc_frac, float *mrate_frac, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;

  for(k=0;k<ntransitions;k++)
    {
      swath1.NS.liqMassFracTrans[*i][k]=mlwc_frac[k];
      swath1.NS.liqRateFracTrans[*i][k]=mrate_frac[k];
    }
}

void copysfcrainliqfracs1t_(float *sfcrainliq_frac, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  swath1.NS.surfLiqRateFrac[*i]=*sfcrainliq_frac;

}
//end    WSO 8/7/13

void copyd0s1t_(float *dm, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;

  for(k=0;k<nbins;k++)
    {
      swath1.NS.precipTotPSDparamHigh[*i][k]=dm[k];
    }
}

void copyzckus1t_(float *zc, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;

   for(k=0;k<nbins;k++)
    {
      if(zc[k] > -90.)
        swath1.NS.correctedReflectFactor[*i][k] = zc[k];
      else
//begin  WSO 9/17/13 standardized missing flags
        swath1.NS.correctedReflectFactor[*i][k] = missing_r4c;
//end    WSO 9/17/13
    }
}

void copynodess1t_(int *node, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  for(k=0;k<5;k++)
    {
      swath1.NS.phaseBinNodes[*i][k]=node[k];
    }
}

void copyrrates2t_(float *rrate, float *rratestd, int *i)
{
  int k;
  //extern L2BCMBT_SWATHS swath;

  for(k=0;k<nbins;k++)
    {
      //swath.MS.precipTotRate[*i][k]=rrate[k];
      //  swath.MS.precipTotRateSigma[*i][k]=rratestd[k];
    }
}

void copypwcs2t_(float *rrate, float *rratestd, int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;

  for(k=0;k<nbins;k++)
    {
//begin WSO 4/18/2013
//changed NS to MS
      //  swath.MS.precipTotWaterCont[*i][k]=rrate[k];
      //swath.MS.precipTotWaterContSigma[*i][k]=rratestd[k];
//end  WSO 4/18/2013
    }
}

//begin  WSO 8/7/13
void copylwcfracs2t_(float *mlwc_frac, float *mrate_frac, int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;

  for(k=0;k<ntransitions;k++)
    {
      //  swath.MS.liqMassFracTrans[*i][k]=mlwc_frac[k];
      //swath.MS.liqRateFracTrans[*i][k]=mrate_frac[k];
    }
}

void copysfcrainliqfracs2t_(float *sfcrainliq_frac, int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;
  //swath.MS.surfLiqRateFrac[*i]=*sfcrainliq_frac;

}
//end   WSO 8/7/13


void copyd0s2t_(float *dm, int *i)
{
  int k;
  //extern L2BCMBT_SWATHS swath;

  for(k=0;k<nbins;k++)
    {
// SFM 05/06/2013 Changed NS to MS to match M.Grecu code from 04/19/2013
      //  swath.MS.precipTotPSDparamHigh[*i][k]=dm[k];
    }
}

void copyzckus2t_(float *zku, float *zka, int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;

   for(k=0;k<nbins;k++)
    {
      //  if(zku[k] > -90.)
      //  swath.MS.correctedReflectFactor[*i][k][0] = zku[k];
      //else
//begin  WSO 9/17/13 standardized missing flags
      // swath.MS.correctedReflectFactor[*i][k][0] = missing_r4c;
//end    WSO 9/17/13
      //if(zka[k] > -90.)
      //  swath.MS.correctedReflectFactor[*i][k][1] = zka[k];
      //else
//begin  WSO 9/17/13 standardized missing flags
      //swath.MS.correctedReflectFactor[*i][k][1] = missing_r4c;
//end    WSO 9/17/13
    }
}
void copynodess2t_(int *node, int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;
  for(k=0;k<5;k++)
    {
      //  swath.MS.phaseBinNodes[*i][k]=node[k];
    }
}

void rewindt_(int *ic)
{
  extern TKINFO       granuleHandle2AKu;
  int status = TKseek(&granuleHandle2AKu, *ic, TK_ABS_SCAN_OFF); 
}

//begin WSO 9/8/13 rewind DPR file
void rewind_dprt_(int *ic)
{
  extern TKINFO       dprtkfile;
    int status_dpr = TKseek(&dprtkfile, *ic, TK_ABS_SCAN_OFF);
}
//end WSO 9/8/13

//  SFM  begin  12/13/2013; add flag to call sequence
void frominputt_(long *st_2adpr)
{
//  SFM  begin  12/13/2013
  extern TKINFO       granuleHandle2AKu;
  extern L2AKu_NS        L2AKuData;
  extern L2BCMBT_SWATHS swath1;
//begin  WSO 9/1/13
  extern L2ADPR_SWATHS dprswathT;
//end    WSO 9/1/13
  int j;
  int status, status_dpr ;
//for diagnostic
  float dummyPIA[49];
  int k, printPIA[49];
//end for diagnostic
//

//  SFM  begin  12/13/2013; add conditional to dpr read
  status=TKreadScan(&granuleHandle2AKu,&L2AKuData);
  if (*st_2adpr == 0) status_dpr=TKreadScan(&dprtkfile,&dprswathT);
//  SFM  begin  12/13/2013

  for( j=0; j<49; j++)
    {
      //swath1.NS.Input.piaEffective[j]=L2AKuData.SRT.pathAtten[j];
      swath1.NS.Input.piaEffective[j]=L2AKuData.SRT.PIAhybrid[j]; //MG 7/31/18, use hybrid PIA
      //begin  WSO 9/5/13 remove flag assignment
      //       swath1.NS.Input.piaEffectiveSigma[j]=-99;
      //end    WSO 9/5/13
      //swath1.NS.Input.piaEffectiveReliabFlag[j]=
	//L2AKuData.SRT.reliabFlag[j];
      swath1.NS.Input.piaEffectiveReliabFlag[j]=
	L2AKuData.SRT.reliabFlagHY[j];                      //WSO  8/2/18, use hybrid flag
      swath1.NS.Input.precipitationType[j]=
	L2AKuData.CSF.typePrecip[j];
      swath1.NS.Input.precipTypeQualityFlag[j]=
	L2AKuData.CSF.qualityTypePrecip[j];
      swath1.NS.Input.surfaceElevation[j]=L2AKuData.PRE.elevation[j];
      swath1.NS.Input.localZenithAngle[j]=L2AKuData.PRE.localZenithAngle[j];
      swath1.NS.Input.surfaceType[j]=L2AKuData.PRE.landSurfaceType[j];
      //begin  WSO 9/28/13 use alternate rain flag that includes missing for bad scans
      //      swath1.NS.Input.precipitationFlag[j]=L2AKuData.PRE.flagPrecip[j];
      //end    WSO 9/28/13
      swath1.NS.Input.surfaceRangeBin[j]=(L2AKuData.PRE.binRealSurface[j]-1)/2; // MG 04/11/2014
      swath1.NS.Input.stormTopBin[j]=(L2AKuData.PRE.binStormTop[j]-1)/2; // MG 04/11/2014
      if(swath1.NS.Input.stormTopBin[j]<0)
	swath1.NS.Input.stormTopBin[j]=missing_i2c;
      swath1.NS.Input.stormTopAltitude[j]=L2AKuData.PRE.heightStormTop[j];
      //begin  WSO 09/30/15 add one bin to the binClutterFreeBottom to temporarily compensate for
      //the subtraction of one bin by the radar team
      //      swath1.NS.Input.lowestClutterFreeBin[j]=
      //	(L2AKuData.PRE.binClutterFreeBottom[j]-1)/2; // MG 04/11/2014
      //begin  WSO 10/19/15 subtract one 125 m bin from binClutterFreeBottom, and 
      //restore V3 definition of lowestClutterFreeBin as a test
      //      swath1.NS.Input.lowestClutterFreeBin[j]=
      //	(L2AKuData.PRE.binClutterFreeBottom[j])/2;
      swath1.NS.Input.lowestClutterFreeBin[j]=
	(L2AKuData.PRE.binClutterFreeBottom[j] - 2)/2;
      //end    WSO 10/15/15
      //end    WSO 09/30/15
      //begin  WSO 9/17/13 correction for two bin average location in combined
      swath1.NS.Input.ellipsoidBinOffset[j]=
	L2AKuData.PRE.ellipsoidBinOffset[j] + 0.125/2.;
      //end    WSO 9/17/13
      //begin  WSO 8/19/13
      swath1.NS.Input.zeroDegAltitude[j] = L2AKuData.VER.heightZeroDeg[j];
      swath1.NS.Input.zeroDegBin[j] = (L2AKuData.VER.binZeroDeg[j]-1)/2; // MG 04/11/2014
      //end    WSO 8/19/13
      //         WSO 8/19/13
    }


  swath1.NS.scanStatus.FractionalGranuleNumber =    
     L2AKuData.scanStatus.FractionalGranuleNumber;
    swath1.NS.scanStatus.SCorientation =
     L2AKuData.scanStatus.SCorientation;
    swath1.NS.scanStatus.acsModeMidScan =
     L2AKuData.scanStatus.acsModeMidScan;
    swath1.NS.scanStatus.dataQuality =
     L2AKuData.scanStatus.dataQuality;
    swath1.NS.scanStatus.dataWarning =
     L2AKuData.scanStatus.dataWarning;
    swath1.NS.scanStatus.geoError =
     L2AKuData.scanStatus.geoError;
    swath1.NS.scanStatus.geoWarning =
     L2AKuData.scanStatus.geoWarning;
    swath1.NS.scanStatus.limitErrorFlag =
     L2AKuData.scanStatus.limitErrorFlag;
    swath1.NS.scanStatus.missing =
     L2AKuData.scanStatus.missing;
    swath1.NS.scanStatus.modeStatus =
     L2AKuData.scanStatus.modeStatus;
    swath1.NS.scanStatus.operationalMode =
     L2AKuData.scanStatus.operationalMode;
    swath1.NS.scanStatus.pointingStatus =
     L2AKuData.scanStatus.pointingStatus;
    swath1.NS.scanStatus.targetSelectionMidScan =
     L2AKuData.scanStatus.targetSelectionMidScan;
//from 2ADPR
   

}

void copyscantimet_(int *i)
{
  extern L2BCMBT_SWATHS swath1;
  extern int DayOfMonth[300], DayOfYear[300], Hour[300], MilliSecond[300],
    Minute[300], Month[300], Second[300], Year[300], SecondOfDay[300];
  extern NAVIGATION navigation[300];

 swath1.NS.ScanTime.DayOfMonth=DayOfMonth[*i];
 swath1.NS.ScanTime.DayOfYear=DayOfYear[*i];
 swath1.NS.ScanTime.Hour=Hour[*i];
 swath1.NS.ScanTime.MilliSecond=MilliSecond[*i];
 swath1.NS.ScanTime.Minute=Minute[*i];
 swath1.NS.ScanTime.Month=Month[*i];
 swath1.NS.ScanTime.Second=Second[*i];
 swath1.NS.ScanTime.Year=Year[*i];
 swath1.NS.ScanTime.SecondOfDay=SecondOfDay[*i];
 memcpy(&swath1.NS.navigation, &navigation[*i], sizeof(NAVIGATION));

//begin WSO 04/07/2013
//added MS swath scantimes
//end WSO 04/07/2013
}

void copypreciptypet_(int *ptype, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  //swath.S1.precipitationType[*i]=*ptype;
}

void copyw10t_(float *w10, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  swath1.NS.tenMeterWindSpeed[*i]=*w10;
}

void copyw10sigmat_(float *w10s, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  swath1.NS.tenMeterWindSigma[*i]=*w10s;
}

void copyw10smallt_(float *w10, int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;
  // swath.MS.tenMeterWindSpeed[*i]=*w10;
}

void copyw10smallsigmat_(float *w10s, int *i)
{
  int k;
  //extern L2BCMBT_SWATHS swath;
  //swath.MS.tenMeterWindSigma[*i]=*w10s;
}

void writedprscant_(void)
{
  int ret;
  ret= TKwriteScan(&dprtkfile,&dprswathT);
}

//  begin  SFM  12/26/2013
void write_emptyt_(void)

//    brief utility to put empty keyword into output file header
//    when needed
{
  char emptygranuletext[100];

  strcpy(emptygranuletext,"EMPTY") ;
  TKsetMetaString(&ctkfile, "FileHeader", "EmptyGranule", 
                  emptygranuletext);	  

  //LW 05/03/18
  TKsetMetaString(&ctkfile, "FileHeader", "SatelliteName", "TRMM");

}
//  end    SFM  12/26/2013

//  begin  SFM  11/27/2013
void writescant_(void)
{
  int ret;
  char emptygranuletext[100];

  TKgetMetaString(&ctkfile, "FileHeader", "EmptyGranule", 
                  emptygranuletext);	  
  if (strncmp(emptygranuletext,"NOT_EMPTY",9) == 0)
      ret= TKwriteScan(&ctkfile,&swath1);
}
//  end    SFM  11/27/2013

void copysfcairtemps1t_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  swath1.NS.surfaceAirTemperature[*i]=*sfcVar;
}

void copysfcairtemps2t_(float *sfcVar, int *i)
{
  int k;
  //extern L2BCMBT_SWATHS swath;
  // swath.MS.surfaceAirTemperature[*i]=*sfcVar;
}

void copysfcairpresss1t_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  swath1.NS.surfaceAirPressure[*i]=*sfcVar;
}

void copysfcairpresss2t_(float *sfcVar, int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;
  // swath.MS.surfaceAirPressure[*i]=*sfcVar;
}

void copyskintemps1t_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  swath1.NS.skinTemperature[*i]=*sfcVar;
}

void copyskintemps2t_(float *sfcVar, int *i)
{  
  int k;
  //  extern L2BCMBT_SWATHS swath;
  // swath.MS.skinTemperature[*i]=*sfcVar;
}

//write skin temperature estimate uncertainty
void copyskintempsigmas1t_(float *skinsigma, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;

  swath1.NS.skinTempSigma[*i] = *skinsigma;
}

void copyskintempsigmas2t_(float *skinsigma, int *i)
{  
  int k;
  //  extern L2BCMBT_SWATHS swath;

  //swath.MS.skinTempSigma[*i] = *skinsigma;
}

//write column vapor estimate uncertainty
void copycolumnvaporsigmas1t_(float *colvaporsigma, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;

  swath1.NS.columnVaporSigma[*i] = *colvaporsigma;
}

void copycolumnvaporsigmas2t_(float *colvaporsigma, int *i)
{  
  int k;
  // extern L2BCMBT_SWATHS swath;

  //swath.MS.columnVaporSigma[*i] = *colvaporsigma;
}

//write column cloud liquid estimate uncerainty
void copycolumncloudliqsigmas1t_(float *colcldsigma, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  
  swath1.NS.columnCloudLiqSigma[*i] = *colcldsigma;
}

void copycolumncloudliqsigmas2t_(float *colcldsigma, int *i)
{  
  int k;
  // extern L2BCMBT_SWATHS swath;

  // swath.MS.columnCloudLiqSigma[*i] = *colcldsigma;
}

//write algorithm type flag
void copyalgotypes1t_(int *algotype, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;

  swath1.NS.FLG.algoType[*i] = *algotype;
}

void copyalgotypes2t_(int *algotype, int *i)
{
  int k;
  //extern L2BCMBT_SWATHS swath;

  //  swath.MS.FLG.algoType[*i] = *algotype;
}


//write error of non-raining data fit
void copyerrorofdatafits1t_(float *erroroffit, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  
  swath1.NS.errorOfDataFit[*i] = *erroroffit;
}

void copyerrorofdatafits2t_(float *erroroffit, int *i)
{  
  int k;
  //extern L2BCMBT_SWATHS swath;

  //swath.MS.errorOfDataFit[*i] = *erroroffit;
}

void copysfcemissouts1t_(float *tbout, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<9;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swath1.NS.surfEmissivity[*i][k]=tbout[k];
    else
      swath1.NS.surfEmissivity[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//  for(k=13;k<15;k++)
//    swath1.NS.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copysfcemissouts1sigmat_(float *tbout, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<9;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swath1.NS.surfEmissSigma[*i][k]=tbout[k];
    else
      swath1.NS.surfEmissSigma[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//  for(k=13;k<15;k++)
//    swath1.NS.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copysfcemissouts2t_(float *tbout, int *i)
{
  int k;
  //extern L2BCMBT_SWATHS swath;
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  // for(k=0;k<13;k++)
//begin  WSO 9/16/13
  // if(tbout[k] > 0.)
  //    swath.MS.surfEmissivity[*i][k]=tbout[k];
  //  else
  //    swath.MS.surfEmissivity[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//    for(k=13;k<15;k++)
//      swath.MS.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copysfcemissouts2sigmat_(float *tbout, int *i)
{
  int k;
  // extern L2BCMBT_SWATHS swath;
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  // for(k=0;k<13;k++)
//begin  WSO 9/16/13
  // if(tbout[k] > 0.)
  //    swath.MS.surfEmissSigma[*i][k]=tbout[k];
  //  else
  //    swath.MS.surfEmissSigma[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//    for(k=13;k<15;k++)
//      swath.MS.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copytbouts1t_(float *tbout, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
//begin  WSO 9/16/13
  for(k=0;k<9;k++)
    if(tbout[k] > -90.)
      swath1.NS.simulatedBrightTemp[*i][k]=tbout[k];
    else
      swath1.NS.simulatedBrightTemp[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels  
//  for(k=13;k<15;k++)
//    swath1.NS.simulatedBrightTemp[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
}

void copytbouts2t_(float *tbout, int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;
//begin  WSO 9/16/13
  // for(k=0;k<13;k++)
  //   if(tbout[k] > -90.)
  //    swath.MS.simulatedBrightTemp[*i][k]=tbout[k];
  //   else
  //    swath.MS.simulatedBrightTemp[*i][k]=missing_r4c;
  //for(k=0;k<2;k++)
  //  swath.MS.simulatedBrightTemp[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//  for(k=13;k<15;k++)
//    swath.MS.simulatedBrightTemp[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end  WSO 9/16/13
}

void copyrainflags1t_(int *sfcVar, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  swath1.NS.Input.precipitationFlag[*i]=*sfcVar;
}

void copyrainflags2t_(int *sfcVar, int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;
  //swath.MS.Input.precipitationFlag[*i][0]=*sfcVar;
  //swath.MS.Input.precipitationFlag[*i][1]=*sfcVar;
}

//begin  WSO 8/20/14 write new ioquality flags
void copyioqualitys1t_(int *sfcVar, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  swath1.NS.FLG.ioQuality[*i]=*sfcVar;
}

void copyioqualitys2t_(int *sfcVar, int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;
  // swath.MS.FLG.ioQuality[*i]=*sfcVar;
}
//end    WSO 8/20/14
//
//begin  WSO 3/17/17 write snow ice cover flags
void copysnowices1t_(int *sfcVar, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  swath1.NS.Input.snowIceCover[*i]=*sfcVar;
}

void copysnowices2t_(int *sfcVar, int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;
  // swath.MS.Input.snowIceCover[*i]=*sfcVar;
}
//end    WSO 3/17/17

void copysfcliqfracts1t_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;
  swath1.NS.surfLiqRateFrac[*i]=*sfcVar;
}

void copysfcliqfracts2t_(float *sfcVar, int *i)
{
  int k;
  // extern L2BCMBT_SWATHS swath;
  // swath.MS.surfLiqRateFrac[*i]=*sfcVar;
}

void copycldwaters1t_(float *var1d, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;

  for(k=0;k<nbins;k++)
    {
        swath1.NS.cloudLiqWaterCont[*i][k]=var1d[k];
    }
}

void copycldwaters2t_(float *var1d, int *i)
{
  int k;
  // extern L2BCMBT_SWATHS swath; 

  for(k=0;k<nbins;k++)
    {
      //    swath.MS.cloudLiqWaterCont[*i][k]=var1d[k];
    }
}

void copycldices1t_(float *var1d, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;

  for(k=0;k<nbins;k++)
    {
      swath1.NS.cloudIceWaterCont[*i][k]=var1d[k];
    }
}

void copycldices2t_(float *var1d, int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;

  for(k=0;k<nbins;k++)
    {
      //  swath.MS.cloudIceWaterCont[*i][k]=var1d[k];
    }
}

//begin  WSO 9/5/13 new copy routine for SRT and DSRT pia effective sigma's
void copysigmapias1t_(float *sigmapia, int *i)
{
   extern L2BCMBT_SWATHS swath1;
//diagnostic
//  printf("in writeCMBOut i: %5i,  sigmapia: %10.4f\n", *sigmapia, *i);
//end diagnostic
  swath1.NS.Input.piaEffectiveSigma[*i] = *sigmapia;
}
void copysigmapias2t_(float *sigmapiaku, float *sigmapiaka, int *i)
{
  //  extern L2BCMBT_SWATHS swath;
//    diagnostic
//  printf("in writeCMBOut i: %5i,  sigmapiaku: %10.4f,  sigmapiaka: %10.4f\n", 
//     *sigmapiaku, *sigmapiaka, *i);
//end diagnostic
  //  swath.MS.Input.piaEffectiveSigma[*i][0] = *sigmapiaku;
  //    swath.MS.Input.piaEffectiveSigma[*i][1] = *sigmapiaka;
}
//end    WSO 9/5/13

//write principal components
void copyprincomps1t_(float *princomp, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;

  for(k=0;k<5;k++)
    {
      swath1.NS.aPriori.prinComp[*i][k] = princomp[k];
    }
}
//
void copyprincomps2t_(float *princomp, int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;

  for(k=0;k<5;k++)
    {
      // swath.MS.aPriori.prinComp[*i][k] = princomp[k];
    }
}

//write profile class
void copyprofclasss1t_(int *profclass, int *i)
{
  extern L2BCMBT_SWATHS swath1;

  swath1.NS.aPriori.profClass[*i] = *profclass;
}

void copyprofclasss2t_(int *profclass, int *i)
{
  //  extern L2BCMBT_SWATHS swath;

  //swath.MS.aPriori.profClass[*i] = *profclass;
}

//write surface precip bias ratio
void copysurfprecipbiasratios1t_(float *biasratio, int *i)
{ 
  extern L2BCMBT_SWATHS swath1;

  swath1.NS.aPriori.surfPrecipBiasRatio[*i] = *biasratio;
}

void copysurfprecipbiasratios2t_(float *biasratio, int *i)
{
  //  extern L2BCMBT_SWATHS swath;

  // swath.MS.aPriori.surfPrecipBiasRatio[*i] = *biasratio;
} 


//write initial log10 of the PSD intercept
void copyinitnws1t_(float *initlogNw, int *n9, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;

  for(k=0;k<9;k++)
    {
     if(n9[k]>1 && n9[k]<=88)
       {
         swath1.NS.aPriori.initNw[*i][k] = initlogNw[n9[k]-1];
       }
     else
       {
         swath1.NS.aPriori.initNw[*i][k] = missing_r4c;
       }
    }
}

void copyinitnws2t_(float *initlogNw, int *n9, int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;

  for(k=0;k<9;k++)
    {
     if(n9[k]>1 && n9[k]<=88)
       {
	 //     swath.MS.aPriori.initNw[*i][k] = initlogNw[n9[k]-1];
       }
     else
       { 
         //swath.MS.aPriori.initNw[*i][k] = missing_r4c;
       }
    }
}

//write sub-footprint variability parameter
void copysubfootvariabilitys1t_(float *subfoot, int *i)
{ 
  extern L2BCMBT_SWATHS swath1;

  swath1.NS.nubfPIAfactor[*i] = *subfoot;
} 
  
void copysubfootvariabilitys2t_(float *subfoot, int *i)
{
  //  extern L2BCMBT_SWATHS swath;
    
  //swath.MS.nubfPIAfactor[*i] = *subfoot;
}

//write multiple scattering flag
void copymultiscatcalcs1t_(int *multiscat, int *i)
{ 
  extern L2BCMBT_SWATHS swath1;
  
  swath1.NS.FLG.multiScatCalc[*i] = *multiscat;
}

void copymultiscatcalcs2t_(int *multiscat, int *i)
{ 
  // extern L2BCMBT_SWATHS swath;
  
  //swath.MS.FLG.multiScatCalc[*i] = *multiscat;
}

//write multiple scattering surface parameter
void copymultiscatsurfaces1t_(float *multisfc, int *i)
{
  extern L2BCMBT_SWATHS swath1;

  swath1.NS.multiScatMaxContrib[*i] = *multisfc;
}

void copymultiscatsurfaces2t_(float *multisfc, int *i)
{
  //  extern L2BCMBT_SWATHS swath;

  //swath.MS.multiScatMaxContrib[*i] = *multisfc;
}

//
//begin  WSO 2/8/17 copy routine for measured sigma-zeros
void copysigmazeros1t_(float *sigmazeroku, int *i)
{
  extern L2BCMBT_SWATHS swath1;
  swath1.NS.Input.sigmaZeroMeasured[*i] = *sigmazeroku;
}
void copysigmazeros2t_(float *sigmazeroku, float *sigmazeroka, int *i)
{
  //    extern L2BCMBT_SWATHS swath;
//        swath.MS.Input.sigmaZeroMeasured[*i][0] = *sigmazeroku;
//            swath.MS.Input.sigmaZeroMeasured[*i][1] = *sigmazeroka;
//          swath.MS.Input.sigmaZeroMeasured[*i] = *sigmazeroka;
}
//end    WSO 2/8/17

//begin  WSO 8/19/13 modified copy routines to include nodes
void copylognws1t_(float *logNw, int *n9, int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;

  for(k=0;k<9;k++)
    {
//begin WSO 9/7/14 added upper limit for n9 in next line
      if(n9[k]>1 && n9[k]<=88)
        {
         swath1.NS.PSDparamLowNode[*i][k] = n9[k]-1;
	 swath1.NS.precipTotPSDparamLow[*i][k][0]=logNw[n9[k]-1];
        }
      else
        {
         swath1.NS.PSDparamLowNode[*i][k] = missing_i2c;
         swath1.NS.precipTotPSDparamLow[*i][k][0]= missing_r4c;
        }
    }
}

void copylognws2t_(float *logNw, int *n9,int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;
  
  for(k=0;k<9;k++)
    {
      if(n9[k]>0)
        {
	  //         swath.MS.PSDparamLowNode[*i][k] = n9[k]-1;
	 //	 swath.MS.precipTotPSDparamLow[*i][k][0]=logNw[n9[k]-1];
        }
      else
        {
	  //         swath.MS.PSDparamLowNode[*i][k] = missing_i2c;
	  //         swath.MS.precipTotPSDparamLow[*i][k][0]= missing_r4c;
        }
    }
}
//end    WSO 8/19/13

//begin  WSO 8/19/13 add mu as second low-resolution parameter
void copymus1t_(float *mu_prof, int *n9,int *i)
{
  int k;
  extern L2BCMBT_SWATHS swath1;

  for(k=0;k<9;k++)
    {
      if(n9[k]>1)
        {
         swath1.NS.precipTotPSDparamLow[*i][k][1]=mu_prof[n9[k]-1];
        }
      else
        {
         swath1.NS.precipTotPSDparamLow[*i][k][1]=missing_r4c;
        }
    }
}

void copymus2t_(float *mu_prof, int *n9,int *i)
{
  int k;
  //  extern L2BCMBT_SWATHS swath;

  for(k=0;k<9;k++)
    {
      if(n9[k]>0)
        {
	  //   swath.MS.precipTotPSDparamLow[*i][k][1]=mu_prof[n9[k]-1];
        }
      else
        {
	  // swath.MS.precipTotPSDparamLow[*i][k][1]=missing_r4c;
        }
    }
}
//end    WSO 8/19/13


void idiot_checkt_(int *number, char *ident)
{
printf(" sfm idiot check %i %s \n",number,ident);
}
