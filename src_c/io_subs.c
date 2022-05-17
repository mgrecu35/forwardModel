//  SFM  04/06/2013  Code changes from M.Grecu
//  SFM  06/27/2013  Parameter name changes from W.Olson; reduce unused code
//  SFM  09/04/2013  Eliminated folllowing routines : openoutputfile32_ ;
//                    openinputfile32_ ; closedproutputfile32_ ; read1c21dot2 ;
//                    writedprscan32_ ; readscan32_ ;  writescan32_ ;
//                    dsphere_ ; copypreciptype32_ ; copytbout32_
//                    They are not used. Do not have jobname required for PPS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <hdf.h>
#include <mfhdf.h>
#include "TKheaders.h"
#include "TK_2BCMB.h"
#ifdef GFOR 
extern int __nbinmod_MOD_imemb;
#define nbins __nbinmod_MOD_imemb
#endif

#ifdef IFORT 
extern int nbinmod_mp_nbin_;
#define nbins nbinmod_mp_nbin_
#endif

//begin WSO 04/07/2013
//Note that many changes of the structures were made to conform
//with TKIO 3.31.0
//All S1 and S2 were changed to NS and MS, respectively
//Also, an Input structure was created under TKIO 3.31.0; see below
//end WSO 04/07/2013 

TKINFO ctkfile32;
extern TKINFO dprtkfile;
TKINFO ctkfileIn32;
L2BCMB_SWATHS swath;
L2BCMBX_SWATHS swathx;
L2ADPR_SWATHS dprswath;
L2BCMB_SWATHS swath1;

void setlatlon2adpr32_(float *lat, float *lon, float *sfcPrecip, 
                       float *sfcPrecipStd, float *piaOut)
{
  int i;
  extern L2BCMB_SWATHS swath;


  for(i=0;i<49;i++)
    {
      dprswath.NS.Latitude[i]=lat[i];
      dprswath.NS.Longitude[i]=lon[i];
      if(dprswath.NS.Longitude[i]>180)
	dprswath.NS.Longitude[i]-=360;
      dprswath.NS.SLV.precipRateNearSurface[i]=sfcPrecip[i];
      dprswath.NS.SRT.pathAtten[i]=piaOut[i];
     
    }
}

void copyrrates132_(float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

 
  for(k=0;k<nbins;k++)
    {
      swath.NS.precipTotRate[*i][k]=rrate[k];
      swath.NS.precipTotRateSigma[*i][k]=rratestd[k];
    }
}

void copyzka32_(float *zka, int *i)
{
  int k;
  extern L2ADPR_SWATHS dprswath;


  for(k=0;k<nbins;k++)
    {
      dprswath.MS.PRE.zFactorMeasured[*i][2*k]=(int)(zka[k]*100);
      dprswath.MS.PRE.zFactorMeasured[*i][2*k+1]=(int)(zka[k]*100);
    }
}

void copypiaka32_(float *piaKa, int *i)
{
  int k;
  extern L2ADPR_SWATHS dprswath;

 
  for(k=0;k<nbins;k++)
    {
      dprswath.MS.SRT.pathAtten[*i]=*piaKa;
    }
}

void copytruerrate32_(float *rrate, int *i)
{
  int k;
  extern L2ADPR_SWATHS dprswath;

 
  for(k=0;k<nbins;k++)
    {
      dprswath.NS.SLV.precipRate[*i][2*k]=(int)(rrate[k]*100);
      dprswath.NS.SLV.precipRate[*i][2*k+1]=(int)(rrate[k]*100);
    }
}

void copypwcs132_(float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

 

  for(k=0;k<nbins;k++)
    {
      swath.NS.precipTotWaterCont[*i][k]=rrate[k];
      swath.NS.precipTotWaterContSigma[*i][k]=rratestd[k];
    }
}

void copyd0s132_(float *dm, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

 

  for(k=0;k<nbins;k++)
    {
      swath.NS.precipTotPSDparamHigh[*i][k]=dm[k];
    }
}

void copyzckus132_(float *zc, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;

   for(k=0;k<nbins;k++)
    {
      swath.NS.correctedReflectFactor[*i][k]=zc[k];
    }
}

void copynodess132_(int *node, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath;
  for(k=0;k<5;k++)
    {
      swath.NS.phaseBinNodes[*i][k]=node[k];
    }
}

void copyscantime32_(int *i)
{
  extern L2BCMB_SWATHS swath;
  extern int DayOfMonth[300], DayOfYear[300], Hour[300], MilliSecond[300],
    Minute[300], Month[300], Second[300], Year[300], SecondOfDay[300];

 swath.NS.ScanTime.DayOfMonth=DayOfMonth[*i];
 swath.NS.ScanTime.DayOfYear=DayOfYear[*i];
 swath.NS.ScanTime.Hour=Hour[*i];
 swath.NS.ScanTime.MilliSecond=MilliSecond[*i];
 swath.NS.ScanTime.Minute=Minute[*i];
 swath.NS.ScanTime.Month=Month[*i];
 swath.NS.ScanTime.Second=Second[*i];
 swath.NS.ScanTime.Year=Year[*i];
}

void get_scAngle(float *scAngle)
{
  int i;
  float sc_angle[49]={17.0340,16.3270,15.6170,14.9050,14.1990,
		      13.4870,12.7780,12.0670,11.3600,10.6480,9.93700,9.22900,
		      8.52000,7.80400,7.09800,6.38900,5.67700,4.96900,4.25900,
		      3.55100,2.84100, 2.13100,1.42200,0.711000,0.000000,
		      -0.711000,-1.42000,-2.12900,-2.84100,-3.54800,
		      -4.25800,-4.97000,-5.67700,
		      -6.38900,-7.09700,-7.81000,-8.51700,
		      -9.23100,-9.93700,
		      -10.6480,-11.3590,-12.0680,
		      -12.7770,-13.4870,-14.1950,
		      -14.9040,-15.6130,-16.3280,-17.0320};
  for(i=0;i<49;i++)
    scAngle[i]=sc_angle[i];
}
