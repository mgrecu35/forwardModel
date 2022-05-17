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
#include "TK_2BCMB_hdf5.h"
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

L2BCMB_SWATHS swath;
L2ADPR_SWATHS dprswath;
L2ADPRX_SWATHS dprxswath;
L2BCMB_SWATHS swath1;
L2AKu_NS      L2AKuData;
L2AKuX_FS      L2AKuDataX;

void openoutputfile_(char *jobname, char *fname)
{
  char ranstring[80], ranMsg[80];
  int status;
  printf(" Output file = %s \n",fname);

  int ret;
//  SFM  04/06/2013  Changed file type to 2BCMB   WSO 04/07/2013
//  SFM  09/04/2013  Changed jobid to jobname
//  SFM  09/11/2013  Moved metadta settings to centralized location

  ret = TKopen(fname, "2BCMB", TKWRITE, "HDF5", jobname, &ctkfile,1); //WSO 04/07/2013
}

void openinputfile_(char *jobname, char *fname)
{
  int ret;

//  SFM  04/06/2013  Changed file type to 2BCMB   WSO 04/07/2013
  ret = TKopen(fname, "2BCMB", TKREAD, "HDF5", jobname, &ctkfileIn, 1);
  printf("%s %i \n",fname,ret);
}

void readscan_(void)
{
  int ret, i;
  ret= TKreadScan(&ctkfileIn,&swath1);

}

void closeoutputfile_(void)
{
  int ret;
  ret=TKclose(&ctkfile);
}

void closedproutputfile_(void)
{
  int ret;
  ret=TKclose(&dprtkfile);
}

void setlatlons1_(float *lat, float *lon, float *sfcPrecip, 
                  float *sfcPrecipStd, float *piaOut)
{
  int i;
  extern L2BCMB_SWATHS swath;

  for(i=0;i<49;i++)
    {
      swath.NS.Latitude[i]=lat[i];
      swath.NS.Longitude[i]=lon[i];
      if(swath.NS.Longitude[i]>180)
	swath.NS.Longitude[i]-=360;
      swath.NS.surfPrecipTotRate[i]=sfcPrecip[i];
      swath.NS.surfPrecipTotRateSigma[i]=sfcPrecipStd[i];
      swath.NS.pia[i]=piaOut[i];
    }
}

void setlatlons2_(float *lat, float *lon, float *sfcPrecip, 
                  float *sfcPrecipStd, float *piaOutKu, float *piaOutKa)
{
  int i;
  extern L2BCMB_SWATHS swath;

  for(i=12;i<37;i++)
    {
      swath.MS.Latitude[i-12]=lat[i];
      swath.MS.Longitude[i-12]=lon[i];
      if(swath.MS.Longitude[i-12]>180)
	swath.MS.Longitude[i-12]-=360;
      swath.MS.surfPrecipTotRate[i-12]=sfcPrecip[i];
      swath.MS.surfPrecipTotRateSigma[i-12]=sfcPrecipStd[i];
      swath.MS.pia[i-12][0]=piaOutKu[i];
      swath.MS.pia[i-12][1]=piaOutKa[i];
    }
}

void copyrrates1_(float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<nbins;k++)
    {
      swath.NS.precipTotRate[*i][k]=rrate[k];
      swath.NS.precipTotRateSigma[*i][k]=rratestd[k];
    }
}

void copysflfract_(float *lfract, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  swath.NS.surfLiqRateFrac[*i]=*lfract;
  if(*i>=12 && *i<=37)
    swath.MS.surfLiqRateFrac[*i-12]=*lfract;
  
}

void copyzka_(float *zka, int *i)
{
  int k;
  extern L2ADPR_SWATHS dprswath;

  for(k=0;k<nbins;k++)
    {
      dprswath.MS.PRE.zFactorMeasured[*i][2*k]=(int)(zka[k]*100);
      dprswath.MS.PRE.zFactorMeasured[*i][2*k+1]=(int)(zka[k]*100);
      //      printf("%g ",zka[k]);
    }
}

void copypiaka_(float *piaKa, int *i)
{
  int k;
  extern L2ADPR_SWATHS dprswath;
 
  for(k=0;k<nbins;k++)
    {
      dprswath.MS.SRT.pathAtten[*i]=*piaKa;
    }
}

void copytruerrate_(float *rrate, int *i)
{
  int k;
  extern L2ADPR_SWATHS dprswath;

  for(k=0;k<nbins;k++)
    {
      dprswath.NS.SLV.precipRate[*i][2*k]=(int)(rrate[k]*100);
      dprswath.NS.SLV.precipRate[*i][2*k+1]=(int)(rrate[k]*100);
    }
}

//begin  WSO 8/30/13
void copyenvsfqvs1_(float *envQv, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.NS.surfaceVaporDensity[*i]=envQv[nbins-1];
}

void copyenvsfqvs2_(float *envQv, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.MS.surfaceVaporDensity[*i]=envQv[nbins-1];
}

void copyenvqvs1_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  for(k=0;k<10;k++)
      swath.NS.vaporDensity[*i][k]=envQv[envnodes[k]-1];
}

void copyenvqvs2_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  for(k=0;k<10;k++)
      swath.MS.vaporDensity[*i][k]=envQv[envnodes[k]-1];
}

void copyenvpresss1_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  for(k=0;k<10;k++)
    swath.NS.airPressure[*i][k]=envQv[envnodes[k]-1];
}

void copyenvpresss2_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  for(k=0;k<10;k++)
    swath.MS.airPressure[*i][k]=envQv[envnodes[k]-1];
}

void copyenvtemps1_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  for(k=0;k<10;k++)
     {
      swath.NS.envParamNode[*i][k]=envnodes[k]-1;
      swath.NS.airTemperature[*i][k]=envQv[envnodes[k]-1];
     }
}

void copyenvtemps2_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  for(k=0;k<10;k++)
    {
     swath.MS.envParamNode[*i][k]=envnodes[k]-1;
     swath.MS.airTemperature[*i][k]=envQv[envnodes[k]-1];
    }
}

void copyenvsftemps1_(float *envQv, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.NS.surfaceAirTemperature[*i]=envQv[nbins-1];
}

void copyenvsftemps2_(float *envQv, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.MS.surfaceAirTemperature[*i]=envQv[nbins-1];
}

//end    WSO 8/30/13

void copypwcs1_(float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<nbins;k++)
    {
      swath.NS.precipTotWaterCont[*i][k]=rrate[k];
      swath.NS.precipTotWaterContSigma[*i][k]=rratestd[k];
    }
}

//begin  WSO 8/7/13
void copylwcfracs1_(float *mlwc_frac, float *mrate_frac, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<ntransitions;k++)
    {
      swath.NS.liqMassFracTrans[*i][k]=mlwc_frac[k];
      swath.NS.liqRateFracTrans[*i][k]=mrate_frac[k];
    }
}

void copysfcrainliqfracs1_(float *sfcrainliq_frac, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.NS.surfLiqRateFrac[*i]=*sfcrainliq_frac;

}
//end    WSO 8/7/13

void copyd0s1_(float *dm, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<nbins;k++)
    {
      swath.NS.precipTotPSDparamHigh[*i][k]=dm[k];
    }
}

void copyzckus1_(float *zc, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

   for(k=0;k<nbins;k++)
    {
      if(zc[k] > -90.)
        swath.NS.correctedReflectFactor[*i][k] = zc[k];
      else
//begin  WSO 9/17/13 standardized missing flags
        swath.NS.correctedReflectFactor[*i][k] = missing_r4c;
//end    WSO 9/17/13
    }
}

void copynodess1_(int *node, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  for(k=0;k<5;k++)
    {
      swath.NS.phaseBinNodes[*i][k]=node[k];
    }
}

void copyrrates2_(float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<nbins;k++)
    {
      swath.MS.precipTotRate[*i][k]=rrate[k];
      swath.MS.precipTotRateSigma[*i][k]=rratestd[k];
    }
}

void copypwcs2_(float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<nbins;k++)
    {
//begin WSO 4/18/2013
//changed NS to MS
      swath.MS.precipTotWaterCont[*i][k]=rrate[k];
      swath.MS.precipTotWaterContSigma[*i][k]=rratestd[k];
//end  WSO 4/18/2013
    }
}

//begin  WSO 8/7/13
void copylwcfracs2_(float *mlwc_frac, float *mrate_frac, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<ntransitions;k++)
    {
      swath.MS.liqMassFracTrans[*i][k]=mlwc_frac[k];
      swath.MS.liqRateFracTrans[*i][k]=mrate_frac[k];
    }
}

void copysfcrainliqfracs2_(float *sfcrainliq_frac, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.MS.surfLiqRateFrac[*i]=*sfcrainliq_frac;

}
//end   WSO 8/7/13


void copyd0s2_(float *dm, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<nbins;k++)
    {
// SFM 05/06/2013 Changed NS to MS to match M.Grecu code from 04/19/2013
      swath.MS.precipTotPSDparamHigh[*i][k]=dm[k];
    }
}

void copyzckus2_(float *zku, float *zka, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

   for(k=0;k<nbins;k++)
    {
      if(zku[k] > -90.)
        swath.MS.correctedReflectFactor[*i][k][0] = zku[k];
      else
//begin  WSO 9/17/13 standardized missing flags
        swath.MS.correctedReflectFactor[*i][k][0] = missing_r4c;
//end    WSO 9/17/13
      if(zka[k] > -90.)
        swath.MS.correctedReflectFactor[*i][k][1] = zka[k];
      else
//begin  WSO 9/17/13 standardized missing flags
        swath.MS.correctedReflectFactor[*i][k][1] = missing_r4c;
//end    WSO 9/17/13
    }
}
void copynodess2_(int *node, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  for(k=0;k<5;k++)
    {
      swath.MS.phaseBinNodes[*i][k]=node[k];
    }
}

void rewind_(int *ic)
{
  extern TKINFO       granuleHandle2AKu;
  int status = TKseek(&granuleHandle2AKu, *ic, TK_ABS_SCAN_OFF); 
}

//begin WSO 9/8/13 rewind DPR file
void rewind_dpr_(int *ic)
{
  extern TKINFO       dprtkfile;
    int status_dpr = TKseek(&dprtkfile, *ic, TK_ABS_SCAN_OFF);
}
//end WSO 9/8/13

//  SFM  begin  12/13/2013; add flag to call sequence
void frominput_(long *st_2adpr)
{
//  SFM  begin  12/13/2013
  extern TKINFO       granuleHandle2AKu;
  extern L2AKu_NS        L2AKuData;
  extern L2BCMB_SWATHS swath;
//begin  WSO 9/1/13
  extern L2ADPR_SWATHS dprswath;
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
  if (*st_2adpr == 0) status_dpr=TKreadScan(&dprtkfile,&dprswath);
//  SFM  begin  12/13/2013

  for( j=0; j<49; j++)
    {
      //swath.NS.Input.piaEffective[j]=L2AKuData.SRT.pathAtten[j];
      swath.NS.Input.piaEffective[j]=L2AKuData.SRT.PIAhybrid[j];  //MG  7/31/18, use hybrid PIA
//begin  WSO 9/5/13 remove flag assignment
//       swath.NS.Input.piaEffectiveSigma[j]=-99;
//end    WSO 9/5/13
   //   swath.NS.Input.piaEffectiveReliabFlag[j]=
	// L2AKuData.SRT.reliabFlag[j];
      swath.NS.Input.piaEffectiveReliabFlag[j]=
	L2AKuData.SRT.reliabFlagHY[j];                              //WSO  8/2/18 use hybrid flag
      swath.NS.Input.precipitationType[j]=
	L2AKuData.CSF.typePrecip[j];
      swath.NS.Input.precipTypeQualityFlag[j]=
	L2AKuData.CSF.qualityTypePrecip[j];
      swath.NS.Input.surfaceElevation[j]=L2AKuData.PRE.elevation[j];
      swath.NS.Input.localZenithAngle[j]=L2AKuData.PRE.localZenithAngle[j];
      swath.NS.Input.surfaceType[j]=L2AKuData.PRE.landSurfaceType[j];
//begin  WSO 9/28/13 use alternate rain flag that includes missing for bad scans
//      swath.NS.Input.precipitationFlag[j]=L2AKuData.PRE.flagPrecip[j];
//end    WSO 9/28/13
      swath.NS.Input.surfaceRangeBin[j]=(L2AKuData.PRE.binRealSurface[j]-1)/2; // MG 04/11/2014
      swath.NS.Input.stormTopBin[j]=(L2AKuData.PRE.binStormTop[j]-1)/2; // MG 04/11/2014
      if(swath.NS.Input.stormTopBin[j]<0)
	swath.NS.Input.stormTopBin[j]=missing_i2c;
      swath.NS.Input.stormTopAltitude[j]=L2AKuData.PRE.heightStormTop[j];
//begin  WSO 09/30/15 add one bin to the binClutterFreeBottom to temporarily compensate for
//the subtraction of one bin by the radar team
//      swath.NS.Input.lowestClutterFreeBin[j]=
//	(L2AKuData.PRE.binClutterFreeBottom[j]-1)/2; // MG 04/11/2014
//begin  WSO 10/19/15 subtract one 125 m bin from binClutterFreeBottom, and 
//restore V3 definition of lowestClutterFreeBin as a test
//      swath.NS.Input.lowestClutterFreeBin[j]=
//	(L2AKuData.PRE.binClutterFreeBottom[j])/2;
      swath.NS.Input.lowestClutterFreeBin[j]=
	(L2AKuData.PRE.binClutterFreeBottom[j] - 2)/2;
//end    WSO 10/15/15
//end    WSO 09/30/15
//begin  WSO 9/17/13 correction for two bin average location in combined
      swath.NS.Input.ellipsoidBinOffset[j]=
	    L2AKuData.PRE.ellipsoidBinOffset[j] + 0.125/2.;
//end    WSO 9/17/13
//begin  WSO 8/19/13
      swath.NS.Input.zeroDegAltitude[j] = L2AKuData.VER.heightZeroDeg[j];
      swath.NS.Input.zeroDegBin[j] = (L2AKuData.VER.binZeroDeg[j]-1)/2; // MG 04/11/2014
//end    WSO 8/19/13
      if(j>=12 && j<37)
	{
	  swath.MS.Input.surfaceElevation[j-12]=
	    L2AKuData.PRE.elevation[j];
	  swath.MS.Input.localZenithAngle[j-12]=
	    L2AKuData.PRE.localZenithAngle[j];
	  swath.MS.Input.surfaceType[j-12]=
	    L2AKuData.PRE.landSurfaceType[j];
//begin  WSO 9/28/13 use alternate rain flag that includes missing for bad scans
//	  swath.MS.Input.precipitationFlag[j-12][0]=
//	    L2AKuData.PRE.flagPrecip[j];
//	  swath.MS.Input.precipitationFlag[j-12][1]=
//	    L2AKuData.PRE.flagPrecip[j];
//end    WSO 9/28/13
	  swath.MS.Input.surfaceRangeBin[j-12][0]=
	    (L2AKuData.PRE.binRealSurface[j]-1)/2; // MG 04/11/2014
	  swath.MS.Input.surfaceRangeBin[j-12][1]=
	    (L2AKuData.PRE.binRealSurface[j]-1)/2; // MG 04/11/2014
	  swath.MS.Input.stormTopBin[j-12][0]=
	    (L2AKuData.PRE.binStormTop[j]-1)/2; // MG 04/11/2014
	  swath.MS.Input.stormTopBin[j-12][1]=  // MG 04/11/2014
	    (L2AKuData.PRE.binStormTop[j]-1)/2; // MG 04/11/2014
	  if(swath.MS.Input.stormTopBin[j-12][0]<0)
	    swath.MS.Input.stormTopBin[j-12][0]=missing_i2c;
	  if(swath.MS.Input.stormTopBin[j-12][1]<0)
	    swath.MS.Input.stormTopBin[j-12][1]=missing_i2c;
	  swath.MS.Input.stormTopAltitude[j-12][0]=
	    L2AKuData.PRE.heightStormTop[j];
	  swath.MS.Input.stormTopAltitude[j-12][1]=
	    L2AKuData.PRE.heightStormTop[j];
//begin  WSO 09/30/15 add one bin to the binClutterFreeBottom to temporarily compensate for
//the subtraction of one bin by the radar team
//	  swath.MS.Input.lowestClutterFreeBin[j-12][0]=
//    (L2AKuData.PRE.binClutterFreeBottom[j]-1)/2; // MG 04/11/2014
//	  swath.MS.Input.lowestClutterFreeBin[j-12][1]=
//	    (L2AKuData.PRE.binClutterFreeBottom[j]-1)/2; // MG 04/11/2014
//begin  WSO 10/19/15 subtract one 125 m bin from binClutterFreeBottom, and 
// restore V3 definition of lowestClutterFreeBin in test
//	  swath.MS.Input.lowestClutterFreeBin[j-12][0]=
//	    (L2AKuData.PRE.binClutterFreeBottom[j])/2;
//	  swath.MS.Input.lowestClutterFreeBin[j-12][1]=
//	    (L2AKuData.PRE.binClutterFreeBottom[j])/2;
	  swath.MS.Input.lowestClutterFreeBin[j-12][0]=
	    (L2AKuData.PRE.binClutterFreeBottom[j] - 2)/2;
	  swath.MS.Input.lowestClutterFreeBin[j-12][1]=
	    (L2AKuData.PRE.binClutterFreeBottom[j] - 2)/2;
//end    WSO 10/19/15
//end    WSO 09/30/15
//begin  WSO 9/17/13 correction for two bin average location in combined
	  swath.MS.Input.ellipsoidBinOffset[j-12][0]=
	    L2AKuData.PRE.ellipsoidBinOffset[j] + 0.125/2.;
	  swath.MS.Input.ellipsoidBinOffset[j-12][1]=
	    L2AKuData.PRE.ellipsoidBinOffset[j] + 0.125/2.;
//end    WSO 9/17/13
//begin  WSO 9/5/13 reset pia's using DPR output
	  swath.MS.Input.piaEffective[j-12][0]=  
	    dprswath.NS.SRT.pathAtten[j];
	  swath.MS.Input.piaEffective[j-12][1]=
	    dprswath.MS.SRT.pathAtten[j-12];
//begin  WSO 9/5/13 remove flag assignment
//	  swath.MS.Input.piaEffectiveSigma[j-12][0]=-99;
//end    WSO 9/5/13
	  swath.MS.Input.piaEffectiveReliabFlag[j-12][0]=
	    dprswath.NS.SRT.reliabFlag[j];
	  swath.MS.Input.piaEffectiveReliabFlag[j-12][1]=
	    dprswath.MS.SRT.reliabFlag[j-12];
//end    WSO 9/5/13
	  swath.MS.Input.precipitationType[j-12]=
	    L2AKuData.CSF.typePrecip[j];
	  swath.MS.Input.precipTypeQualityFlag[j-12]=
	    L2AKuData.CSF.qualityTypePrecip[j];
//begin  WSO 8/19/13 need to update toolkit
          swath.MS.Input.zeroDegAltitude[j-12] =
            L2AKuData.VER.heightZeroDeg[j];
          swath.MS.Input.zeroDegBin[j-12][0] =
            (L2AKuData.VER.binZeroDeg[j]-1)/2;   // MG 04/11/2014
//end    WSO 8/19/13
	}
//diagnostic assignment
          dummyPIA[j] = dprswath.NS.SRT.pathAtten[j];
//end diagnostic 
    }

//diagnostic
//       if(L2AKuData.Latitude[24] > 30. &&  L2AKuData.Latitude[24] < 40. && L2AKuData.Longitude[24] > -165. && L2AKuData.Longitude[24] <-155.)
//         {
//           for(k=0;k<49;k++)
//             if(dummyPIA[k] < -99.)
//               {
//                 printPIA[k] = 99;
//               }
//             else
//               printPIA[k] = dummyPIA[k]*10.;
//           printf("lon: %10.2f,  ", L2AKuData.Longitude[24]);
//           for(k=0;k<49;k++)
//             printf("%2i", printPIA[k]);
//           printf("\n");
//         }
//end diagnostic

//begin  WSO 9/1/13 scanStatus variables copied from 2AKu
    swath.NS.scanStatus.FractionalGranuleNumber =    
     L2AKuData.scanStatus.FractionalGranuleNumber;
    swath.NS.scanStatus.SCorientation =
     L2AKuData.scanStatus.SCorientation;
    swath.NS.scanStatus.acsModeMidScan =
     L2AKuData.scanStatus.acsModeMidScan;
    swath.NS.scanStatus.dataQuality =
     L2AKuData.scanStatus.dataQuality;
    swath.NS.scanStatus.dataWarning =
     L2AKuData.scanStatus.dataWarning;
    swath.NS.scanStatus.geoError =
     L2AKuData.scanStatus.geoError;
    swath.NS.scanStatus.geoWarning =
     L2AKuData.scanStatus.geoWarning;
    swath.NS.scanStatus.limitErrorFlag =
     L2AKuData.scanStatus.limitErrorFlag;
    swath.NS.scanStatus.missing =
     L2AKuData.scanStatus.missing;
    swath.NS.scanStatus.modeStatus =
     L2AKuData.scanStatus.modeStatus;
    swath.NS.scanStatus.operationalMode =
     L2AKuData.scanStatus.operationalMode;
    swath.NS.scanStatus.pointingStatus =
     L2AKuData.scanStatus.pointingStatus;
    swath.NS.scanStatus.targetSelectionMidScan =
     L2AKuData.scanStatus.targetSelectionMidScan;
//from 2ADPR
    swath.MS.scanStatus.FractionalGranuleNumber =
     dprswath.MS.scanStatus.FractionalGranuleNumber;
    swath.MS.scanStatus.SCorientation =
     dprswath.MS.scanStatus.SCorientation;
    swath.MS.scanStatus.acsModeMidScan =
     dprswath.MS.scanStatus.acsModeMidScan;
    swath.MS.scanStatus.dataQuality =
     dprswath.MS.scanStatus.dataQuality;
    swath.MS.scanStatus.dataWarning =
     dprswath.MS.scanStatus.dataWarning;
    swath.MS.scanStatus.geoError =
     dprswath.MS.scanStatus.geoError;
    swath.MS.scanStatus.geoWarning =
     dprswath.MS.scanStatus.geoWarning;
    swath.MS.scanStatus.limitErrorFlag =
     dprswath.MS.scanStatus.limitErrorFlag;
    swath.MS.scanStatus.missing =
     dprswath.MS.scanStatus.missing;
    swath.MS.scanStatus.modeStatus =
     dprswath.MS.scanStatus.modeStatus;
    swath.MS.scanStatus.operationalMode =
     dprswath.MS.scanStatus.operationalMode;
    swath.MS.scanStatus.pointingStatus =
     dprswath.MS.scanStatus.pointingStatus;
    swath.MS.scanStatus.targetSelectionMidScan =
     dprswath.MS.scanStatus.targetSelectionMidScan;
//end    WSO 9/1/13

}

void copyscantime_(int *i)
{
  extern L2BCMB_SWATHS swath;
  extern int DayOfMonth[300], DayOfYear[300], Hour[300], MilliSecond[300],
    Minute[300], Month[300], Second[300], Year[300], SecondOfDay[300];
  extern NAVIGATION navigation[300];

 swath.NS.ScanTime.DayOfMonth=DayOfMonth[*i];
 swath.NS.ScanTime.DayOfYear=DayOfYear[*i];
 swath.NS.ScanTime.Hour=Hour[*i];
 swath.NS.ScanTime.MilliSecond=MilliSecond[*i];
 swath.NS.ScanTime.Minute=Minute[*i];
 swath.NS.ScanTime.Month=Month[*i];
 swath.NS.ScanTime.Second=Second[*i];
 swath.NS.ScanTime.Year=Year[*i];
 swath.NS.ScanTime.SecondOfDay=SecondOfDay[*i];
 memcpy(&swath.NS.navigation, &navigation[*i], sizeof(NAVIGATION));

//begin WSO 04/07/2013
//added MS swath scantimes
 swath.MS.ScanTime.DayOfMonth=DayOfMonth[*i];
 swath.MS.ScanTime.DayOfYear=DayOfYear[*i];
 swath.MS.ScanTime.Hour=Hour[*i];
 swath.MS.ScanTime.MilliSecond=MilliSecond[*i];
 swath.MS.ScanTime.Minute=Minute[*i];
 swath.MS.ScanTime.Month=Month[*i];
 swath.MS.ScanTime.Second=Second[*i];
 swath.MS.ScanTime.Year=Year[*i];
 swath.MS.ScanTime.SecondOfDay=SecondOfDay[*i];
 memcpy(&swath.MS.navigation, &navigation[*i], sizeof(NAVIGATION));
//end WSO 04/07/2013
}

void copypreciptype_(int *ptype, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  //swath.S1.precipitationType[*i]=*ptype;
}

void copyw10_(float *w10, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.NS.tenMeterWindSpeed[*i]=*w10;
}

void copyw10sigma_(float *w10s, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.NS.tenMeterWindSigma[*i]=*w10s;
}

void copyw10small_(float *w10, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.MS.tenMeterWindSpeed[*i]=*w10;
}

void copyw10smallsigma_(float *w10s, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.MS.tenMeterWindSigma[*i]=*w10s;
}

void writedprscan_(void)
{
  int ret;
  ret= TKwriteScan(&dprtkfile,&dprswath);
}

//  begin  SFM  12/26/2013
void write_empty_(void)

//    brief utility to put empty keyword into output file header
//    when needed
{
  char emptygranuletext[100];

  strcpy(emptygranuletext,"EMPTY") ;
  TKsetMetaString(&ctkfile, "FileHeader", "EmptyGranule", 
                  emptygranuletext);	  
}
//  end    SFM  12/26/2013

//  begin  SFM  11/27/2013
void writescan_(void)
{
  int ret;
  char emptygranuletext[100];

  // TKgetMetaString(&ctkfile, "FileHeader", "EmptyGranule", 
  //                emptygranuletext);	  
  // if (strncmp(emptygranuletext,"NOT_EMPTY",9) == 0)
  ret= TKwriteScan(&ctkfile,&swath);
}
//  end    SFM  11/27/2013

void copysfcairtemps1_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.NS.surfaceAirTemperature[*i]=*sfcVar;
}

void copysfcairtemps2_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.MS.surfaceAirTemperature[*i]=*sfcVar;
}

void copysfcairpresss1_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.NS.surfaceAirPressure[*i]=*sfcVar;
}

void copysfcairpresss2_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.MS.surfaceAirPressure[*i]=*sfcVar;
}

void copyskintemps1_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.NS.skinTemperature[*i]=*sfcVar;
}

void copyskintemps2_(float *sfcVar, int *i)
{  
  int k;
  extern L2BCMB_SWATHS swath;
  swath.MS.skinTemperature[*i]=*sfcVar;
}

//write skin temperature estimate uncertainty
void copyskintempsigmas1_(float *skinsigma, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  swath.NS.skinTempSigma[*i] = *skinsigma;
}

void copyskintempsigmas2_(float *skinsigma, int *i)
{  
  int k;
  extern L2BCMB_SWATHS swath;

  swath.MS.skinTempSigma[*i] = *skinsigma;
}

//write column vapor estimate uncertainty
void copycolumnvaporsigmas1_(float *colvaporsigma, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  swath.NS.columnVaporSigma[*i] = *colvaporsigma;
}

void copycolumnvaporsigmas2_(float *colvaporsigma, int *i)
{  
  int k;
  extern L2BCMB_SWATHS swath;

  swath.MS.columnVaporSigma[*i] = *colvaporsigma;
}

//write column cloud liquid estimate uncerainty
void copycolumncloudliqsigmas1_(float *colcldsigma, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  
  swath.NS.columnCloudLiqSigma[*i] = *colcldsigma;
}

void copycolumncloudliqsigmas2_(float *colcldsigma, int *i)
{  
  int k;
  extern L2BCMB_SWATHS swath;

  swath.MS.columnCloudLiqSigma[*i] = *colcldsigma;
}

//write algorithm type flag
void copyalgotypes1_(int *algotype, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  swath.NS.FLG.algoType[*i] = *algotype;
}

void copyalgotypes2_(int *algotype, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  swath.MS.FLG.algoType[*i] = *algotype;
}


//write error of non-raining data fit
void copyerrorofdatafits1_(float *erroroffit, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  
  swath.NS.errorOfDataFit[*i] = *erroroffit;
}

void copyerrorofdatafits2_(float *erroroffit, int *i)
{  
  int k;
  extern L2BCMB_SWATHS swath;

  swath.MS.errorOfDataFit[*i] = *erroroffit;
}

void copysfcemissouts1_(float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swath.NS.surfEmissivity[*i][k]=tbout[k];
    else
      swath.NS.surfEmissivity[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//  for(k=13;k<15;k++)
//    swath.NS.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copysfcemissouts1sigma_(float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swath.NS.surfEmissSigma[*i][k]=tbout[k];
    else
      swath.NS.surfEmissSigma[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//  for(k=13;k<15;k++)
//    swath.NS.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copysfcemissouts2_(float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swath.MS.surfEmissivity[*i][k]=tbout[k];
    else
      swath.MS.surfEmissivity[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//    for(k=13;k<15;k++)
//      swath.MS.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copysfcemissouts2sigma_(float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swath.MS.surfEmissSigma[*i][k]=tbout[k];
    else
      swath.MS.surfEmissSigma[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//    for(k=13;k<15;k++)
//      swath.MS.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copytbouts1_(float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
//begin  WSO 9/16/13
  for(k=0;k<13;k++)
    if(tbout[k] > -90.)
      swath.NS.simulatedBrightTemp[*i][k]=tbout[k];
    else
      swath.NS.simulatedBrightTemp[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels  
//  for(k=13;k<15;k++)
//    swath.NS.simulatedBrightTemp[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
}

void copytbouts2_(float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
//begin  WSO 9/16/13
  for(k=0;k<13;k++)
    if(tbout[k] > -90.)
      swath.MS.simulatedBrightTemp[*i][k]=tbout[k];
    else
      swath.MS.simulatedBrightTemp[*i][k]=missing_r4c;
  //for(k=0;k<2;k++)
  //  swath.MS.simulatedBrightTemp[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//  for(k=13;k<15;k++)
//    swath.MS.simulatedBrightTemp[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end  WSO 9/16/13
}

void copyrainflags1_(int *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.NS.Input.precipitationFlag[*i]=*sfcVar;
}

void copyrainflags2_(int *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.MS.Input.precipitationFlag[*i][0]=*sfcVar;
  swath.MS.Input.precipitationFlag[*i][1]=*sfcVar;
}

//begin  WSO 8/20/14 write new ioquality flags
void copyioqualitys1_(int *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.NS.FLG.ioQuality[*i]=*sfcVar;
}

void copyioqualitys2_(int *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.MS.FLG.ioQuality[*i]=*sfcVar;
}
//end    WSO 8/20/14
//
//begin  WSO 3/17/17 write snow ice cover flags
void copysnowices1_(int *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.NS.Input.snowIceCover[*i]=*sfcVar;
}

void copysnowices2_(int *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.MS.Input.snowIceCover[*i]=*sfcVar;
}
//end    WSO 3/17/17

void copysfcliqfracts1_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.NS.surfLiqRateFrac[*i]=*sfcVar;
}

void copysfcliqfracts2_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  swath.MS.surfLiqRateFrac[*i]=*sfcVar;
}

void copycldwaters1_(float *var1d, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<nbins;k++)
    {
        swath.NS.cloudLiqWaterCont[*i][k]=var1d[k];
    }
}

void copycldwaters2_(float *var1d, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath; 

  for(k=0;k<nbins;k++)
    {
        swath.MS.cloudLiqWaterCont[*i][k]=var1d[k];
    }
}

void copycldices1_(float *var1d, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<nbins;k++)
    {
      swath.NS.cloudIceWaterCont[*i][k]=var1d[k];
    }
}

void copycldices2_(float *var1d, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<nbins;k++)
    {
      swath.MS.cloudIceWaterCont[*i][k]=var1d[k];
    }
}

//begin  WSO 9/5/13 new copy routine for SRT and DSRT pia effective sigma's
void copysigmapias1_(float *sigmapia, int *i)
{
  extern L2BCMB_SWATHS swath;
//diagnostic
//  printf("in writeCMBOut i: %5i,  sigmapia: %10.4f\n", *sigmapia, *i);
//end diagnostic
  swath.NS.Input.piaEffectiveSigma[*i] = *sigmapia;
}
void copysigmapias2_(float *sigmapiaku, float *sigmapiaka, int *i)
{
  extern L2BCMB_SWATHS swath;
//    diagnostic
//  printf("in writeCMBOut i: %5i,  sigmapiaku: %10.4f,  sigmapiaka: %10.4f\n", 
//     *sigmapiaku, *sigmapiaka, *i);
//end diagnostic
    swath.MS.Input.piaEffectiveSigma[*i][0] = *sigmapiaku;
    swath.MS.Input.piaEffectiveSigma[*i][1] = *sigmapiaka;
}
//end    WSO 9/5/13

//write principal components
void copyprincomps1_(float *princomp, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<5;k++)
    {
      swath.NS.aPriori.prinComp[*i][k] = princomp[k];
    }
}
//
void copyprincomps2_(float *princomp, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<5;k++)
    {
     swath.MS.aPriori.prinComp[*i][k] = princomp[k];
    }
}

//write profile class
void copyprofclasss1_(int *profclass, int *i)
{
  extern L2BCMB_SWATHS swath;

  swath.NS.aPriori.profClass[*i] = *profclass;
}

void copyprofclasss2_(int *profclass, int *i)
{
  extern L2BCMB_SWATHS swath;

  swath.MS.aPriori.profClass[*i] = *profclass;
}

//write surface precip bias ratio
void copysurfprecipbiasratios1_(float *biasratio, int *i)
{ 
  extern L2BCMB_SWATHS swath;

  swath.NS.aPriori.surfPrecipBiasRatio[*i] = *biasratio;
}

void copysurfprecipbiasratios2_(float *biasratio, int *i)
{
  extern L2BCMB_SWATHS swath;

  swath.MS.aPriori.surfPrecipBiasRatio[*i] = *biasratio;
} 


//write initial log10 of the PSD intercept
void copyinitnws1_(float *initlogNw, int *n9, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<9;k++)
    {
     if(n9[k]>1 && n9[k]<=88)
       {
         swath.NS.aPriori.initNw[*i][k] = initlogNw[n9[k]-1];
       }
     else
       {
         swath.NS.aPriori.initNw[*i][k] = missing_r4c;
       }
    }
}

void copyinitnws2_(float *initlogNw, int *n9, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<9;k++)
    {
     if(n9[k]>1 && n9[k]<=88)
       {
         swath.MS.aPriori.initNw[*i][k] = initlogNw[n9[k]-1];
       }
     else
       { 
         swath.MS.aPriori.initNw[*i][k] = missing_r4c;
       }
    }
}

//write sub-footprint variability parameter
void copysubfootvariabilitys1_(float *subfoot, int *i)
{ 
  extern L2BCMB_SWATHS swath;

  swath.NS.nubfPIAfactor[*i] = *subfoot;
} 
  
void copysubfootvariabilitys2_(float *subfoot, int *i)
{
  extern L2BCMB_SWATHS swath;
    
  swath.MS.nubfPIAfactor[*i] = *subfoot;
}

//write multiple scattering flag
void copymultiscatcalcs1_(int *multiscat, int *i)
{ 
  extern L2BCMB_SWATHS swath;
  
  swath.NS.FLG.multiScatCalc[*i] = *multiscat;
}

void copymultiscatcalcs2_(int *multiscat, int *i)
{ 
  extern L2BCMB_SWATHS swath;
  
  swath.MS.FLG.multiScatCalc[*i] = *multiscat;
}

//write multiple scattering surface parameter
void copymultiscatsurfaces1_(float *multisfc, int *i)
{
  extern L2BCMB_SWATHS swath;

  swath.NS.multiScatMaxContrib[*i] = *multisfc;
}

void copymultiscatsurfaces2_(float *multisfc, int *i)
{
  extern L2BCMB_SWATHS swath;

  swath.MS.multiScatMaxContrib[*i] = *multisfc;
}

//
//begin  WSO 2/8/17 copy routine for measured sigma-zeros
void copysigmazeros1_(float *sigmazeroku, int *i)
{
  extern L2BCMB_SWATHS swath;
  swath.NS.Input.sigmaZeroMeasured[*i] = *sigmazeroku;
}
void copysigmazeros2_(float *sigmazeroku, float *sigmazeroka, int *i)
{
    extern L2BCMB_SWATHS swath;
//        swath.MS.Input.sigmaZeroMeasured[*i][0] = *sigmazeroku;
//            swath.MS.Input.sigmaZeroMeasured[*i][1] = *sigmazeroka;
          swath.MS.Input.sigmaZeroMeasured[*i] = *sigmazeroka;
}
//end    WSO 2/8/17

//begin  WSO 8/19/13 modified copy routines to include nodes
void copylognws1_(float *logNw, int *n9, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<9;k++)
    {
//begin WSO 9/7/14 added upper limit for n9 in next line
      if(n9[k]>1 && n9[k]<=88)
        {
         swath.NS.PSDparamLowNode[*i][k] = n9[k]-1;
	 swath.NS.precipTotPSDparamLow[*i][k][0]=logNw[n9[k]-1];
        }
      else
        {
         swath.NS.PSDparamLowNode[*i][k] = missing_i2c;
         swath.NS.precipTotPSDparamLow[*i][k][0]= missing_r4c;
        }
    }
}

void copylognws2_(float *logNw, int *n9,int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  
  for(k=0;k<9;k++)
    {
      if(n9[k]>0)
        {
         swath.MS.PSDparamLowNode[*i][k] = n9[k]-1;
	 swath.MS.precipTotPSDparamLow[*i][k][0]=logNw[n9[k]-1];
        }
      else
        {
         swath.MS.PSDparamLowNode[*i][k] = missing_i2c;
	 swath.MS.precipTotPSDparamLow[*i][k][0]= missing_r4c;
        }
    }
}
//end    WSO 8/19/13

//begin  WSO 8/19/13 add mu as second low-resolution parameter
void copymus1_(float *mu_prof, int *n9,int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<9;k++)
    {
      if(n9[k]>1)
        {
         swath.NS.precipTotPSDparamLow[*i][k][1]=mu_prof[n9[k]-1];
        }
      else
        {
         swath.NS.precipTotPSDparamLow[*i][k][1]=missing_r4c;
        }
    }
}

void copymus2_(float *mu_prof, int *n9,int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

  for(k=0;k<9;k++)
    {
      if(n9[k]>0)
        {
         swath.MS.precipTotPSDparamLow[*i][k][1]=mu_prof[n9[k]-1];
        }
      else
        {
         swath.MS.precipTotPSDparamLow[*i][k][1]=missing_r4c;
        }
    }
}
//end    WSO 8/19/13


void idiot_check_(int *number, char *ident)
{
printf(" sfm idiot check %i %s \n",number,ident);
}
