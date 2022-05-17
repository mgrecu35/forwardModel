//  SFM 04/06/2013  Code module added in merge from M.Grecu's code
//  SFM 04/16/2013  Modified file name lengths from 80, 100, and 256 to 
//                    1000 for all occurrences; including "malloc" useage
//  SFM 05/06/2013  Code from LW to pass jobname into modules
//  SFM 06/19/2013  Undocumented code changes from M.Grecu
//  SFM  6/27/2013  Parameter name changes from W.Olson; reduce unused code
//  SFM  7/30/2013  zFactorMeasured value reduced by factor 100 to accommodate
//                    changes in 2AKu data format
//  LAW  5/02/2014  Received update readpr.GPM.model_05_01_14.c.gz from Bill Olson

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf.h>
#include <mfhdf.h>
//begin  WSO 9/5/13 include math library for power calculation
#include <math.h>
//end    WSO 9/5/13
#include "TKheaders.h"
//begin  aaa LAW 9/4/13 1GMI to 1CGMI_r1
#include "TK_1CGMI.h"
//end    aaa LAW 9/4/13
#ifdef GFOR 
extern int __nbinmod_MOD_nbin;
#define nbins __nbinmod_MOD_nbin
//begin  WSO 9/15/13 
extern float __missingmod_MOD_missing_r4;
#define missing_r4c __missingmod_MOD_missing_r4
extern short __missingmod_MOD_missing_i2;
#define missing_i2c __missingmod_MOD_missing_i2
extern long __missingmod_MOD_missing_i4;
#define missing_i4c __missingmod_MOD_missing_i4
//end    WSO 9/15/13
#endif

int DayOfMonth[300], DayOfYear[300], Hour[300], MilliSecond[300],
  Minute[300], Month[300], Second[300], Year[300], SecondOfDay[300];
//  begin  LAW  05/08/2014; add navigation data to output file
NAVIGATION navigation[300];
//  end    LAW  05/08/2014

#ifdef IFORT 
extern int nbinmod_mp_nbin_;
#define nbins nbinmod_mp_nbin_
//begin  WSO 9/15/13
extern float missingmod_mp_missing_r4_;
#define missing_r4c missingmod_mp_missing_r4_
extern short missingmod_mp_missing_i2_;
#define missing_i2c missingmod_mp_missing_i2_
extern long missingmod_mp_missing_i4_;
#define missing_i4c missingmod_mp_missing_i4_
//end    WSO 9/15/13
#endif
TKINFO granuleHandle2AKu;
TKINFO tkfileinfo;

TKINFO dprtkfile;
extern L2ADPR_SWATHS dprswath;

//  begin  SFM  09/12/2013
extern TKINFO ctkfile;
int iYearSet=0, iDaySet=0;  // 4/14/MG
//  end  SFM  09/12/2013

//  begin WSO 8/14/2014
unsigned long getbits(unsigned long x, int p, int n);
unsigned long temp_byte;
//  end   WSO 8/14/2014

//  SFM  begin 12/06/2013; to pass out TK status message
int readgmi_(char *jobname, char *f1cgmi,  int *n1b11, float *tmilow,
             float *tmihigh, float *xlon,float *xlat, float *gmilon, 
	     float *gmilat, int *mm, int *yyyy, int *jday, 
	     int *dd, int *ifile, float *secondOfDay,
	     float *SCLon, float *SCLat, 
	     float *xlonS2, float *xlatS2,
	     float *eiaS1, float *eiaS2)  //4/14/14 MG
{
  int             status_alpha ;
//  SFM  end   12/06/2013
  char            *granuleID;
  int             dataType;
  char            filemode;
  int             status, orbNumber,i, ind, ipia, indh,j, k, nscans, nscanmax;

//  TKINFO          tkfileinfo;
  L1CGMI_SWATHS   data;
  char emptygranuletext[100];
 
  i=0;

  granuleID=(char*)malloc(sizeof(char)*1000);
  strcpy(&granuleID[0],f1cgmi);
  printf(" readgmi: f1cgmi %s \n",f1cgmi);
  
//  SFM  begin 12/06/2013; to pass out TK status message
//begin  LAW  9/4/13 1CGMI to 1CGMI_r1
//aa   status = TKopen(granuleID,"1CGMI", TKREAD, "HDF5", jobname, &tkfileinfo, 1);
    status_alpha = TKopen(granuleID,"1CGMI", TKREAD, "HDF5", jobname, &tkfileinfo, 1);
    if (status_alpha != 0)
      {
       printf("WARNING: Unable to access 1CGMI %i \n",status_alpha);
       return status_alpha ;
      }
//  SFM  end   12/06/2013   
//end    LAW  9/4/13

//  begin  SFM  09/12/2013
    TKtransferMetaData(&tkfileinfo, &ctkfile);
    printf("  STATUS, meta transfer 1CGMI %i  \n",status_alpha);
//  end  SFM  09/12/2013
 
  TKgetMetaString(&tkfileinfo, "FileHeader", "EmptyGranule", 
                  emptygranuletext);	  
  if (strncmp(emptygranuletext,"EMPTY",5) == 0)
     {
      status_alpha = -3;
      printf(" STATUS, 1CGMI has empty granule \n",status_alpha);
      return status_alpha ;
      }
      
  status = TKgetScanCount ( &tkfileinfo, "S1" );
    
  *n1b11=0;
  
  status = TKseek(&tkfileinfo, 0,
                  TK_ABS_SCAN_OFF);
  nscans=0;
  while (TKendOfFile (&tkfileinfo) != TK_EOF)
    {
      status=TKreadScan(&tkfileinfo,&data);
      nscans++;
    }

  printf("   number scans = %i \n",nscans);
  if(*ifile==1)
    {
      status = TKseek(&tkfileinfo, nscans-200,
		      TK_ABS_SCAN_OFF);
    }
  else
    {
      status = TKseek(&tkfileinfo, 0,
		      TK_ABS_SCAN_OFF);
    }
  if(*ifile==3)
    nscanmax=200;
  else
    nscanmax=3300;
      
  ind=0;
  ipia=0;
  indh=0;

  int dday;  //4/14/14 MG
  while (TKendOfFile (&tkfileinfo) != TK_EOF && *n1b11<nscanmax)
    {
      status=TKreadScan(&tkfileinfo,&data);

      *mm=data.S1.ScanTime.Month;
      *yyyy=data.S1.ScanTime.Year;
      *jday=data.S1.ScanTime.DayOfYear;
      *dd=data.S1.ScanTime.DayOfMonth;
      
//  SFM  begin  04/16/2014; for M.Grecu, revision of nodes processing
      if(iYearSet==0)                                     //4/14/14 MG Begin
	{
	  iYearSet=*yyyy;
	  iDaySet=*jday;
	  dday=0;
	}
      else
	{
	  if(iYearSet/4*4==iYearSet)
	    dday=(iYearSet-*yyyy)*366*24*3600;
	  else
	    dday=(iYearSet-*yyyy)*365*24*3600.;
	  dday+=(iDaySet-*jday)*24*3600.;
	}
      secondOfDay[*n1b11]=dday+data.S1.ScanTime.SecondOfDay;     //4/14/14 MG End
      SCLat[*n1b11]=data.S1.SCstatus.SClatitude;
      SCLon[*n1b11]=data.S1.SCstatus.SClongitude;

//  SFM  end  04/16/2014

      (*n1b11)+=1;

      for(j=0; j<221; j++)
	{
	  for(k=0;k<4;k++)
	    {
	      tmihigh[indh]=data.S2.Tc[j][k]/1.;
	      indh++;
	    }
	  for(k=0;k<9;k++)
	    {
	      tmilow[ind]=data.S1.Tc[j][k]/1.;
//  SFM  begin  12/24/2013; diagnostic to insert temperature defaults
//              if (j<110) tmilow[ind] = -999.0 ;
//  SFM  end    12/24/2013
	      ind++;
	    }
	}
      for(j=0;j<221;j++)
	{
	  eiaS1[ipia]=*data.S1.incidenceAngle[j];
	  eiaS2[ipia]=*data.S2.incidenceAngle[j];
	  xlon[ipia]=data.S1.Longitude[j];
	  xlat[ipia]=data.S1.Latitude[j];
	  xlonS2[ipia]=data.S2.Longitude[j];
	  xlatS2[ipia]=data.S2.Latitude[j];
	  ipia++;
	}

    }

//  SFM  begin 12/06/2013; to pass out TK status message
  status_alpha = TKclose(&tkfileinfo);
  free(granuleID);    
  return status_alpha ;
//  SFM  end   12/06/2013; to pass out TK status message
}

void get_scAngle(float *scAngle);
TKINFO       granuleHandle2AKu;
char fileName[1000];

//  SFM  begin 12/12/2013; to pass out TK status message
//  SFM  begin  04/16/2014; for M.Grecu, revision of nodes processing
//begin  WSO 08/18/14; add ioqualityflags pass to main routine
void read2aku_(char *jobname, char *f2aku, int *n1c21, 
	       float *zKuObs, float *zKaObs, float *snrRatioku, float *snrRatioka, 
               float *srtPIAku, float *dsrtPIAku, float *dsrtPIAka,
               float *srtsigmaPIAku, float *dsrtsigmaPIAku, float *dsrtsigmaPIAka,
               float *sigmaZeroKu, float *sigmaZeroKa, float *sclon, float *sclat,//SJM 2/3/15 //SJM 2/3/15
	       float *xlon, float *xlat, int *badRayFlag, int *rainFlagBad,
	       int *nodes, int *raintype, float *scAngle, int *ic,
	       int *nBSize, float *freezH, float *surfaceZku, 
	       int *iLandOcean, float *srtrelPIAku, float *dsrtrelPIA, 
	       float *piaHB, int *ioqualityflagku, int *ioqualityflagdpr,
	       char *f2ADPR, float *dprrain, int *BBbin, int *binRealSurface,
	       float *localZenithAngle, float *elevation, int *status_alpha,
	       float *secondOfDay, int *NSRelibFlag, int *MSRelibFlag,
	       int *snowIceCover, float *seaIceConcentration, float *cBEst)
//end    WSO 08/18/14
//  SFM  end    04/16/2014
//  SFM  end   12/12/2013
{

  extern int DayOfMonth[300], DayOfYear[300], Hour[300], MilliSecond[300],
    Minute[300], Month[300], Second[300], Year[300], SecondOfDay[300];
  int             dataType;
  char            filemode;
  int             status, orbNumber;
  char            foutname[1000], npcafname[1000];
  float           maxr, bb1x1[360][90];
  int i, j, nscan, k, kk, indClutf, k1, izc, inodes, ipia;
  float           q80[nbins], aw, bw, xs;
  int16           firstbin[49];
  float           z80[nbins],x,y;

//begin  SFM  9/04/2013
  L2AKu_NS        L2AKuData;
//  L2AKu_NS_aa1     L2AKuData;
//end   SFM  9/04/2013

//begin  SFM  3/13/2014
  int missing_flag ;
  int dataquality_flag ;
  int modestatus_flag ;
//end    SFM  3/13/2014
    
  int inode=0, ii, nscans, rtype;
  char fname[1000];
  int32 sd_id, sds_id, n_datasets, n_file_attrs,
           rain_index, geo_index, index, status32;
  float pia;
  int32 dim_sizes[10];
  int32 start[3];
  int32 rank, num_type, attributes, istat;
  int itop, ibb;
  int iSurf,iBinEllip,iClutf, ibb1, ibb2; 
//begin    WSO 9/5/13 pia effective sigma variables 
  float				piaAltern, rfactAltern, sigmaAltern, suminvvar;
  float           dummyPIAku[49], dummyPIAka[49];
  int             printPIA[49];
//end      WSO 9/5/13
 
  FILE *fout1,*fout2,*fout3;
  int isfzku;
    char startTime[100], stopTime[100];
//  begin  SFM  11/25/2013
    char emptygranuletext[100];
//  end    SFM  11/25/2013

  get_scAngle(scAngle);
 
  char *ifDPR;
  char dummychar[] = "N"; //WSO 04/07/2013
  
// begin WSO 04/07/2013
//...suppress writing of DPR file
  ifDPR = &dummychar[0];
//  ifDPR=getenv("IF2ADPROUT");
// end WSO 04/07/2013

  strcpy(&fname[0],&f2aku[0]);

  get_scAngle(scAngle);
  int iwrite=0;

  char dprfname[1000];
  strcpy(dprfname,f2ADPR);

//...First pass through, check need to open files
  if(*ic==0)
  {
//aa  begin  SFM  09/04/2013
    status = TKopen(&fname[0], "2AKu", TKREAD, "HDF5", jobname,
		    &granuleHandle2AKu, 1);
//  status = TKopen(&fname[0], "2AKu_aa1", TKREAD, "HDF5", jobname,
//		  &granuleHandle2AKu, 1);
//aa  end  SFM  09/04/2013
//  begin  SFM  09/12/2013
    TKtransferMetaData(&granuleHandle2AKu, &ctkfile);
    printf("  status, meta transfer 2AKu,a  %i  \n",status);
    TKgetMetaString(&granuleHandle2AKu, "FileHeader", "StartGranuleDateTime", 
                    startTime);
    TKgetMetaString(&granuleHandle2AKu, "FileHeader", "StopGranuleDateTime", 
                    stopTime);
    TKsetMetaString(&ctkfile, "FileHeader", "StartGranuleDateTime", startTime);
    TKsetMetaString(&ctkfile, "FileHeader", "StopGranuleDateTime", stopTime);
//  begin  SFM  11/25/2013
    TKgetMetaString(&granuleHandle2AKu, "FileHeader", "EmptyGranule", 
                    emptygranuletext);
    TKsetMetaString(&ctkfile, "FileHeader", "EmptyGranule", emptygranuletext);

    
  if (strncmp(emptygranuletext,"EMPTY",5) == 0)
     {
//      status_alpha = -3;
      printf(" STATUS, 2AKu has empty granule \n",status);
      return  ;
      }

//  SFM  begin 12/12/2013; to pass out TK status message
//  end    SFM  11/25/2013
//  end  SFM  09/12/2013
    if(ifDPR[0]=='Y' && *status_alpha == 0)
      {
//aa  begin  SFM  09/04/2013
//  SFM  start  06/22/2014; for M.Grecu (no longer needed)
//    status = TKopen(&fname[0], "2AKu", TKREAD, "HDF5", jobname,
//		    &granuleHandle2AKu, 1);
//  SFM  end    06/22/2014; for M.Grecu
//  SFM  end   12/12/2013
//  status = TKopen(&fname[0], "2AKu_aa1", TKREAD, "HDF5", jobname,
//		  &granuleHandle2AKu, 1);
//aa  end  SFM  09/04/2013
//  begin  SFM  09/12/2013
    TKtransferMetaData(&granuleHandle2AKu, &ctkfile);
    printf("  status, meta transfer 2AKu,b  %i  \n",status);
//  end  SFM  09/12/2013

      }
    else
      {
//  SFM  begin 12/12/2013; to pass out TK status message
	if(strncmp(dprfname,"nil",3)!=0 && *status_alpha == 0)
	  {
//aa  begin  SFM  09/04/2013
	    *status_alpha = TKopen(dprfname,"2ADPR",TKREAD,"HDF5",jobname,&dprtkfile,1);
//	    status = TKopen(dprfname,"2ADPR_aa1",TKREAD,"HDF5",jobname,&dprtkfile,1);
//aa  end  SFM  09/04/2013
//  begin  SFM  09/12/2013

            if (*status_alpha == 0) 
	      {
               TKtransferMetaData(&dprtkfile, &ctkfile);
               printf("  STATUS, meta transfer 2ADPR %i  \n",*status_alpha);
//  SFM  end   12/12/2013
               int nscans;
               TKgetMetaInt(&dprtkfile,"NS_SwathHeader","NumberScansGranule",&nscans);
               printf(" routine read2aku: nscans, from 2ADPR   %i \n",nscans);

               TKgetMetaString(&dprtkfile, "FileHeader", "EmptyGranule", 
                               emptygranuletext);	  
               if (strncmp(emptygranuletext,"EMPTY",5) == 0)
                  {
                   *status_alpha = -3;
                   printf(" STATUS, 2ADPR has empty granule \n",*status_alpha);
                   }
	      }
	  }
      }
  }

//  SFM  begin   12/13/2013
  if(strncmp(dprfname,"nil",3) == 0) *status_alpha = -2 ;
//  SFM  end     12/13/2013

////  TKgetMetaString(&dprtkfile, "FileHeader", "EmptyGranule", 
////                  emptygranuletext);	  
//  if (
//  if (strncmp(emptygranuletext,"EMPTY",5) == 0)
//     {
//      *status_alpha = -3;
//      printf(" STATUS, 2ADPR has empty granule \n",*status_alpha);
//      }
//  else
  
  i=0;
  izc=0;
  ipia=0;
  inodes=0;
  isfzku=0;
  int iland=0;

  if(*ic==0)
    {
     status = TKseek(&granuleHandle2AKu, *ic, TK_ABS_SCAN_OFF); 
//  SFM  begin 12/13/2013; to pass out TK status message
     if(ifDPR[0]!='Y' && strncmp(dprfname,"nil",3)!=0 && *status_alpha == 0)
        {
         if (*status_alpha == 0) status = TKseek(&dprtkfile, *ic, TK_ABS_SCAN_OFF); 
	 }
//  SFM  end   12/13/2013
    }
  int ifact=1;
  //  float dZku=-0.7, dZka=-0.5;
   float dZku=0.0, dZka=0.0;
  int yyyy, jday, dday; //4/14/14 MG
  while((TKendOfFile (&granuleHandle2AKu) != TK_EOF) && i<*nBSize)
    {
      status=TKreadScan(&granuleHandle2AKu,&L2AKuData);
//  SFM  begin 12/12/2013; to pass out TK status message
      if(ifDPR[0]!='Y' && strncmp(dprfname,"nil",3)!=0 && *status_alpha == 0)
//  SFM  end   12/12/2013
        {
         if (*status_alpha == 0) status = TKreadScan(&dprtkfile,&dprswath);
	 //printf("%s \n","reading DPR");
	 }
//  SFM  end   12/13/2013
      
      DayOfMonth[i]=L2AKuData.ScanTime.DayOfMonth;
      DayOfYear[i]=L2AKuData.ScanTime.DayOfYear;
      Hour[i]=L2AKuData.ScanTime.Hour;
      MilliSecond[i]=L2AKuData.ScanTime.MilliSecond;
      Minute[i]=L2AKuData.ScanTime.Minute;
      Month[i]=L2AKuData.ScanTime.Month;
      Second[i]=L2AKuData.ScanTime.Second;
      SecondOfDay[i]=L2AKuData.ScanTime.SecondOfDay;
      Year[i]=L2AKuData.ScanTime.Year;
//  begin  LAW  05/08/2014; add navigation data to output file
      memcpy(&navigation[i], &L2AKuData.navigation, sizeof(NAVIGATION));
//  end    LAW  05/08/2014
      sclon[i] = L2AKuData.navigation.scLon;
      sclat[i] = L2AKuData.navigation.scLat;

      for(j=0;j<49;j++)
	{

//begin  WSO 8/18/14 initialize ioquality flags
          ioqualityflagku[ipia] = 0;
          ioqualityflagdpr[ipia] = 0;
//end    WSO 8/18/14
//begin  WSO 8/14/14 extract first byte from qualityData, which contain Zku scan quality info
          badRayFlag[ipia] = getbits(L2AKuData.FLG.qualityData[j], 7, 8);
	  BBbin[ipia]=L2AKuData.CSF.binBBPeak[j];
//end    WSO 8/14/14

//begin  WSO 8/20/14 set ioquality flags for freezing level detection (based on Ku only)
          if(badRayFlag[ipia] != 0)  //bad scan
            {
              ioqualityflagku[ipia] = ioqualityflagku[ipia] + 9000;   //add 9000 for bad Ku input data
              if(j>=12 && j<37)
                {
                  ioqualityflagdpr[ipia] = ioqualityflagdpr[ipia] + 9000; //add 9000 for bad Ku input data
                }
            }
          else if(L2AKuData.CSF.binBBPeak[j]-1 > 0)  //scan OK, but and freezing level determined from bright band
            {
              ioqualityflagku[ipia] = ioqualityflagku[ipia] + 0;   //add 0 for freezing level from bright band
              if(j>=12 && j<37)
               {
                  ioqualityflagdpr[ipia] = ioqualityflagdpr[ipia] + 0; //add 0 for freezing level from bright band
               }
            }
          else  //scan OK, but no bright band, so freezing level from analysis
            {
              ioqualityflagku[ipia] = ioqualityflagku[ipia] + 1000;  // add 1000 for freezing level from analysis
              if(j>=12 && j<37)
               {
                  ioqualityflagdpr[ipia] = ioqualityflagdpr[ipia] + 1000; //add 1000 for freezing level from analysis
               }
            }
//end    WSO 8/20/14


//begin  WSO 8/20/14 set ioquality flags for precipitation type classification (based on Ku only)
          if(badRayFlag[ipia] != 0)  //bad scan
            {
              ioqualityflagku[ipia] = ioqualityflagku[ipia] + 90000;   //add 90000 for bad Ku input data
              if(j>=12 && j<37)
                {
                  ioqualityflagdpr[ipia] = ioqualityflagdpr[ipia] + 90000; //add 90000 for bad Ku input data
                }
            }
          else if(L2AKuData.CSF.typePrecip[j] <= 0)  //scan OK, but no feature for classification (no rain)
            {
              ioqualityflagku[ipia] = ioqualityflagku[ipia] + 20000;   //add 20000 for no feature (no rain)
              if(j>=12 && j<37)
               {
                  ioqualityflagdpr[ipia] = ioqualityflagdpr[ipia] + 20000; //add 20000 for no feature (no rain)
               }
            }
          else if(L2AKuData.CSF.typePrecip[j]/10000000 >= 3)  //scan OK, but indeterminate classification
            {
              ioqualityflagku[ipia] = ioqualityflagku[ipia] + 10000;  // add 10000 for indeterminate classification
              if(j>=12 && j<37)
               {
                  ioqualityflagdpr[ipia] = ioqualityflagdpr[ipia] + 10000; //add 10000 for indeterminate classification
               }
            }
//end    WSO 8/20/14



//  SFM  begin 04/16/2014; for M.Grecu, revision procssing of nodes
//begin  WSO 05/01/2014 restore index offset for 1-176 input
	  nodes[inodes]= (L2AKuData.PRE.binStormTop[j]-1)/2; //WSO 05/01/2014
	  inodes++;

	  if(L2AKuData.CSF.binBBPeak[j]-1>0)                 //WSO 05/01/2014 
	    {
	      nodes[inodes]=(L2AKuData.CSF.binBBPeak[j]-1)/2-3; //WSO 05/01/14
	      inodes++;
	      nodes[inodes]=(L2AKuData.CSF.binBBPeak[j]-1)/2;   //WSO 05/01/2014
	      inodes++;
	      nodes[inodes]=(L2AKuData.CSF.binBBPeak[j]-1)/2+2; //WSO 05/01/2014
	      inodes++;
	    }
	  else
	    {
	      nodes[inodes]=(L2AKuData.VER.binZeroDeg[j]-1)/2-3; //WSO 05/01/2014
	      inodes++;
	      nodes[inodes]=(L2AKuData.VER.binZeroDeg[j]-1)/2;   //WSO 05/01/2014
	      inodes++;
	      nodes[inodes]=(L2AKuData.VER.binZeroDeg[j]-1)/2+2; //WSO 05/01/2014
	      inodes++;
	    }
//end   WSO 05/01/2014 restore index offset for 1-176 input
//  SFM  end   04/16/2014; for M.Grecu, revision procssing of nodes
//     if(L2AKuData.PRE.binClutterFreeBottom[j]-1 < 0)  //WSO 4/12/14
//       {
//        printf("lat: %10.2f,  lon: %10.2f  binCluttFree: %5i,  nodes: %5i %5i %5i %5i\n", 
//         L2AKuData.Latitude[j], L2AKuData.Longitude[j], L2AKuData.PRE.binClutterFreeBottom[j]-1, 
//         nodes[inodes-4], nodes[inodes-3], nodes[inodes-2], nodes[inodes-1]);  //WSO 4/12/14
//       }

//begin  WSO 09/30/15 add one 125 m bin to the lowest node to temporarily
//compensate for the subtraction of one bin from L2AKuData.PRE.binClutterFreeBottom by the radar team
//	  nodes[inodes]=(L2AKuData.PRE.binClutterFreeBottom[j]-1)/2-1;      //Feb 2015 MG
//
//begin WSO 10/19/15 subtract one 125 m bin from binClutterFreeBottom,and then use
// V3 node definition in test
//	  nodes[inodes]=(L2AKuData.PRE.binClutterFreeBottom[j])/2-1;      
	  nodes[inodes]=(L2AKuData.PRE.binClutterFreeBottom[j] - 2)/2;      
//end   WSO 10/19/15
//end    WSO 9/30/15
	  inodes++;

//begin  WSO 9/28/13 account for "bad scans"
//begin  WSO 8/18/14 set ioqualityflags; assume quality hinges on presence of Ku input and rain
//detected by Ku, only
          if(badRayFlag[ipia] != 0)  //bad scan
            {
              ioqualityflagku[ipia] = ioqualityflagku[ipia] + 90;   //add 90 for bad Ku input data
              ioqualityflagku[ipia] = ioqualityflagku[ipia] + 9;    //add 9 for no valid CMB output
              if(j>=12 && j<37)
                {
                  ioqualityflagdpr[ipia] = ioqualityflagdpr[ipia] + 90; //add 90 for bad Ku input data
                  ioqualityflagdpr[ipia] = ioqualityflagdpr[ipia] + 9;  //add 9 for no valid CMB output
                }
            }
          else if(L2AKuData.PRE.flagPrecip[j] == 0)  //scan OK, but no rain at Ku
            {
              ioqualityflagku[ipia] = ioqualityflagku[ipia] + 10;   //add 10 for no rain at Ku
              ioqualityflagku[ipia] = ioqualityflagku[ipia] + 0;    //add 0 for valid CMB output
              if(j>=12 && j<37)
               {
                  ioqualityflagdpr[ipia] = ioqualityflagdpr[ipia] + 10; //add 10 for no rain at Ku
                  ioqualityflagdpr[ipia] = ioqualityflagdpr[ipia] + 0;  //add 0 for valid CMB output
               }
            }
//end    WSO 8/18/14

          rainFlagBad[ipia] = L2AKuData.PRE.flagPrecip[j];

//begin  WSO 6/29/14 change criteria from >0 to != 0 in next line
     if(badRayFlag[ipia] != 0)     //indicates bad scan
//end    WSO 6/29/14 
       {
         rainFlagBad[ipia] = missing_i4c;
         for(k=inodes-5;k<inodes;k++) nodes[k] = missing_i4c;
//        printf("lat: %10.2f,  lon: %10.2f  badRayFlag:  %5i,  binCluttFree: %5i,  missing: %5i,  nodes: %5i %5i %5i %5i %5i\n",
//         L2AKuData.Latitude[j], L2AKuData.Longitude[j], badRayFlag[ipia], L2AKuData.PRE.binClutterFreeBottom[j]-1, missing_i4c, 
//         nodes[inodes-5], nodes[inodes-4], nodes[inodes-3], nodes[inodes-2], nodes[inodes-1]); //WSO 4/12/14
        }
     else if(nodes[inodes-5]<0)  //no storm top, but VER gives freezing level
        {
        nodes[inodes-5]=nodes[inodes-4]-1;
//        printf("lat: %10.2f,  lon: %10.2f  binCluttFree: %5i,  missing: %5i,  nodes: %5i %5i %5i %5i %5i\n",
//         L2AKuData.Latitude[j], L2AKuData.Longitude[j], L2AKuData.PRE.binClutterFreeBottom[j]-1, missing_i4c,
//         nodes[inodes-5], nodes[inodes-4], nodes[inodes-3], nodes[inodes-2], nodes[inodes-1]);  //WSO 4/12/14
        }
//end    WSO 9/28/13

//begin  WSO 8/20/14 test SRT to see if it's available and useable
     if(badRayFlag[ipia] != 0)  //bad scan
       {
         ioqualityflagku[ipia] = ioqualityflagku[ipia] + 900;
       }
     else if(L2AKuData.PRE.snRatioAtRealSurface[j] <= 2.)  //Ku channel attenuated
       {
         ioqualityflagku[ipia] = ioqualityflagku[ipia] + 200;
       }
     else if(L2AKuData.SRT.reliabFactorHY[j] <= 2.9) //Ku sigma-zero in noise; WSO 8/2/18 hybrid
       {
         ioqualityflagku[ipia] = ioqualityflagku[ipia] + 100;
       }
     if(j>=12 && j<37)

       {
         if(badRayFlag[ipia] != 0)  //bad scan
           {
             ioqualityflagdpr[ipia] = ioqualityflagdpr[ipia] + 900;
           }
         else if(dprswath.MS.PRE.snRatioAtRealSurface[j-12] <= 2.)  //Ka channel attenuated
           {
             ioqualityflagdpr[ipia] = ioqualityflagdpr[ipia] + 200;
           }
         else if(dprswath.MS.SRT.reliabFactor[j-12] <= 2.9) //Ka sigma-zero in noise of background
           {
             ioqualityflagdpr[ipia] = ioqualityflagdpr[ipia] + 100;
           }
       }
//end    WSO 8/20/14

//begin  WSO 9/5/13 read individual PIA estimates and 
////       reliability factors, and create composite
////       uncertainty estimate
////begin  WSO 9/10/13 read surface SNR estimate at Ku
     snowIceCover[ipia]=L2AKuData.PRE.snowIceCover[j];
     seaIceConcentration[ipia]=L2AKuData.Experimental.seaIceConcentration[j];
     //printf("%5i %8.2f", snowIceCover[ipia], seaIceConcentration[ipia]);
     snrRatioku[ipia] = L2AKuData.PRE.snRatioAtRealSurface[j];
    // srtPIAku[ipia]=L2AKuData.SRT.pathAtten[j];
     srtPIAku[ipia]=L2AKuData.SRT.PIAhybrid[j];  //MG 7/31/18 use hybrid PIA
     srtsigmaPIAku[ipia] = missing_r4c;
     // begin SJM 7/9/14 read Ku sigma_zero
     sigmaZeroKu[ipia] = L2AKuData.PRE.sigmaZeroMeasured[j];
     //printf("%5i %5i %8.2f\n",j,ipia,sigmaZeroKu[ipia]);
     // end SJM 7/9/14
//begin  WSO 11/17/15 exclude cases where reliabFactor is zero
//begin WSO 8/2/18 calculate hybrid PIA uncertainty
     if(L2AKuData.SRT.PIAhybrid[j] > -9990. && (L2AKuData.SRT.reliabFactorHY[j] > 1.e-5 || L2AKuData.SRT.reliabFactorHY[j] < -1.e-5))
       srtsigmaPIAku[ipia] = L2AKuData.SRT.PIAhybrid[j] / L2AKuData.SRT.reliabFactorHY[j];
     else
       srtsigmaPIAku[ipia] = missing_r4c;
//end    WSO 8/2/18
//end    WSO 11/17/15
//end    WSO 9/10/13
//end    WSO 9/5/13
     
	  piaHB[ipia]=L2AKuData.SLV.piaFinal[j];
	 // srtrelPIAku[ipia]=L2AKuData.SRT.reliabFactor[j];
	  srtrelPIAku[ipia]=L2AKuData.SRT.reliabFactorHY[j];  //WSO 8/2/18 use hybrid factor
	 // NSRelibFlag[ipia]=L2AKuData.SRT.reliabFlag[j];
	  NSRelibFlag[ipia]=L2AKuData.SRT.reliabFlagHY[j];    //WSO 8/2/18 use hybrid flag
	  raintype[ipia]=L2AKuData.CSF.typePrecip[j];
	  binRealSurface[ipia]=(L2AKuData.PRE.binRealSurface[j]-1)/2;  //WSO 05/01/2014
	  localZenithAngle[ipia]=L2AKuData.PRE.localZenithAngle[j]; 
	  elevation[ipia]=L2AKuData.PRE.elevation[j];
	  
	  if(L2AKuData.PRE.binStormTop[j]-1<0)                       //WSO 05/01/2014
	    raintype[ipia]=0;
//begin  WSO 09/30/15 temporarily compensate for the subtraction of one
//bin from L2AKuData.PRE.binClutterFreeBottom by the radar team
//	  if(L2AKuData.PRE.binClutterFreeBottom[j]-1==0)             //WSO 05/01/2014
//	     raintype[ipia]=0;
	  if(L2AKuData.PRE.binClutterFreeBottom[j]==0)
	     raintype[ipia]=0;
//end    WSO 09/30/15
	  
	  ifact=1;
//begin  WSO 9/13/13 changed > to >=
	  while(raintype[ipia]/ifact>=10)
	    {
	      ifact=ifact*10;
	    }
	  raintype[ipia]=raintype[ipia]/ifact*100;
	  if(L2AKuData.CSF.binBBPeak[j]>10)
	    raintype[ipia]=100;
	  freezH[ipia]=L2AKuData.VER.heightZeroDeg[j];
	  dprrain[ipia]=dprswath.NS.SLV.precipRateNearSurface[j];
//end    WSO 9/13/13
//begin  WSO 9/5/13 calculate uncertainty estimates at Ku and Ka
//begin  WSO 9/10/13 add surface SNR estimate at Ka
	  if(j>=12 && j<37)
           {
           //printf("%i5 %i5 %8.2f\n",j,ipia,dprswath.MS.PRE.sigmaZeroMeasured[j-12]);
           sigmaZeroKa[ipia] = dprswath.MS.PRE.sigmaZeroMeasured[j-12]; //SJM 7/9/14
           snrRatioka[ipia] = dprswath.MS.PRE.snRatioAtRealSurface[j-12];
             
	   dsrtPIAku[ipia]=dprswath.NS.SRT.pathAtten[j];
	   dsrtPIAka[ipia]=dprswath.MS.SRT.pathAtten[j-12];
	   dsrtrelPIA[ipia]=dprswath.MS.SRT.reliabFactor[j-12];
	   MSRelibFlag[ipia]=dprswath.MS.SRT.reliabFlag[j-12];
//begin  WSO 11/17/15  exclude cases where reliabFactor is zero
      dsrtsigmaPIAku[ipia] = missing_r4c;
      if(dprswath.NS.SRT.pathAtten[j] > -9990. && (dprswath.NS.SRT.reliabFactor[j] > 1.e-5 || dprswath.NS.SRT.reliabFactor[j] < -1.e-5))
         dsrtsigmaPIAku[ipia] = dprswath.NS.SRT.pathAtten[j] / dprswath.NS.SRT.reliabFactor[j];
      else
         dsrtsigmaPIAku[ipia] = missing_r4c;

      dsrtsigmaPIAka[ipia] = missing_r4c;
      if(dprswath.MS.SRT.pathAtten[j-12] > -9990. && (dprswath.MS.SRT.reliabFactor[j-12] > 1.e-5 || dprswath.MS.SRT.reliabFactor[j-12] < -1.e-5))
        dsrtsigmaPIAka[ipia] = dprswath.MS.SRT.pathAtten[j-12] / dprswath.MS.SRT.reliabFactor[j-12];
      else
        dsrtsigmaPIAka[ipia] = missing_r4c;
//end    WSO 11/17/15

           dummyPIAku[j] = dprswath.NS.SRT.pathAtten[j];
           dummyPIAka[j] = dprswath.MS.SRT.pathAtten[j-12];
           }
	  else
           {
          snrRatioka[ipia] = missing_r4c;
	       dsrtPIAku[ipia] = missing_r4c;
	       dsrtPIAka[ipia] = missing_r4c;
          dummyPIAku[j] = missing_r4c;
          dummyPIAka[j] = missing_r4c;
          dsrtsigmaPIAku[ipia] = missing_r4c;
          dsrtsigmaPIAka[ipia] = missing_r4c;
          sigmaZeroKa[ipia] = missing_r4c; //SJM 7/9/14
	  dsrtrelPIA[ipia]=missing_r4c;
	  MSRelibFlag[ipia]=4;
           }

//  SFM  begin 12/12/2013; to pass out TK status message
	  if(strncmp(dprfname,"nil",3)==0 && *status_alpha == 0)
//  SFM  end   12/12/2013
           {
          snrRatioka[ipia] = missing_r4c;
	       dsrtPIAku[ipia] = missing_r4c;
	       dsrtPIAku[ipia] = missing_r4c;
	       dsrtPIAka[ipia] = missing_r4c;
          dummyPIAku[j] = missing_r4c;
          dummyPIAka[j] = missing_r4c;
          dsrtsigmaPIAku[ipia] = missing_r4c;
          dsrtsigmaPIAka[ipia] = missing_r4c;
           }
//end    WSO 9/10/13
//end    WSO 9/5/13
	  xlon[ipia]=L2AKuData.Longitude[j];
	  xlat[ipia]=L2AKuData.Latitude[j];
	  ipia++;
	  for(k=0;k<nbins;k++)
	    {
	      if(L2AKuData.PRE.zFactorMeasured[j][2*k]>0 && 
		 L2AKuData.PRE.zFactorMeasured[j][2*k+1]>0)
		zKuObs[izc]= 0.5*(L2AKuData.PRE.zFactorMeasured[j][2*k]
				  +L2AKuData.PRE.zFactorMeasured
				  [j][2*k+1])+dZku ;
	      else
		if(L2AKuData.PRE.zFactorMeasured[j][2*k]>0)
		  zKuObs[izc]= L2AKuData.PRE.zFactorMeasured[j][2*k]+dZku ;
		else
		  if(L2AKuData.PRE.zFactorMeasured[j][2*k+1]>0)
		    zKuObs[izc]= L2AKuData.PRE.zFactorMeasured[j][2*k+1]+dZku ;
		  else
		    zKuObs[izc]=0;
	      
	      //  SFM  begin 12/12/2013; to pass out TK status message
	      //printf("%i %i z=%g,",j,k,(float) dprswath.MS.PRE.zFactorMeasured[j-12][2*k]);
	      if(strncmp(dprfname,"nil",3)!=0 && *status_alpha == 0)
		{
		  if(j>=12 && j<37)
		    {
		      if(dprswath.MS.PRE.zFactorMeasured[j-12][2*k]>0 && 
			 dprswath.MS.PRE.zFactorMeasured[j-12][2*k+1]>0)
			zKaObs[izc]= 0.5*(dprswath.MS.PRE.zFactorMeasured[j-12][2*k]
					  +dprswath.MS.PRE.zFactorMeasured
					  [j-12][2*k+1])+dZka ;
		      else
			if(dprswath.MS.PRE.zFactorMeasured[j-12][2*k]>0)
			  zKaObs[izc]= dprswath.MS.PRE.zFactorMeasured[j-12][2*k]+dZka ;
			else
			  if(dprswath.MS.PRE.zFactorMeasured[j-12][2*k+1]>0)
			    zKaObs[izc]= dprswath.MS.PRE.zFactorMeasured[j-12][2*k+1]+dZka ;
			  else
			    zKaObs[izc]=0;
		    }
		  else
		    zKaObs[izc]=-99;
		}
	      else
		zKaObs[izc]=-99;

	      izc++;
	    }
	}
//diagnostic
//       if(L2AKuData.Latitude[24] > 30. &&  L2AKuData.Latitude[24] < 40. && L2AKuData.Longitude[24] > -165. && L2AKuData.Longitude[24] <-155.)
//         {
//           for(k=0;k<49;k++)
//             if(dummyPIAku[k] < -99.)
//               {
//                 printPIA[k] = 99;
//               }
//             else
//               printPIA[k] = dummyPIAku[k]*10.;
//           printf("i: %5i, lon: %10.2f,  ", i, L2AKuData.Longitude[24]);
//           for(k=0;k<49;k++)
//             printf("%2i", printPIA[k]);
//           printf("\n");
//         }
//end diagnostic

//  SFM  begin 3/13/2014; diagnostics
      missing_flag     = L2AKuData.scanStatus.missing;
      dataquality_flag = L2AKuData.scanStatus.dataQuality;
      modestatus_flag  = L2AKuData.scanStatus.modeStatus;
//      printf("flag values %i %i %i  %i %i %i \n",i,j,*ic,
//              missing_flag,dataquality_flag,modestatus_flag);
//  SFM  end   3/13/2014

//  SFM  begin  04/16/2014; for M.Grecu, revision of nodes processing
      yyyy=L2AKuData.ScanTime.Year;                 //4/14/14 MG Begin
      jday=L2AKuData.ScanTime.DayOfYear;
      if(iYearSet==0)                  
	{
	  iYearSet=yyyy;
	  iDaySet=jday;
	  dday=0;
	}
      else
	{
	  if(iYearSet/4*4==iYearSet)
	    dday=(iYearSet-yyyy)*366*24*3600;
	  else
	    dday=(iYearSet-yyyy)*365*24*3600.;
	  dday+=(iDaySet-jday)*24*3600.;
	}
      secondOfDay[i]=L2AKuData.ScanTime.SecondOfDay;  //4/14/14 MG End
//  SFM  end  04/16/2014
      i++;
    }

    *n1c21=i;

    return;
}

//  SFM  begin 12/09/2013; to pass out TK status message
int readdprpflag_(char *jobname, char *f2aku, int *n1c21, 
		  float *xlon, float *xlat, int *iLandOcean, 
		  int *flagPrecip, float *sfcRain)
{
  int             status_alpha ;
//  SFM  begin 12/09/2013; to pass out TK status message

  TKINFO       granuleHandle2AKu;
  int             dataType;
  char            filemode;
  int             status, orbNumber;
  char            foutname[1000], npcafname[1000];
  float           maxr, bb1x1[360][90];
  int i, j, nscan, k, kk, indClutf, k1, izc, inodes, ipia;
  float           q80[nbins], aw, bw, xs;
  int16           firstbin[49];
  float           z80[nbins],x,y;
//  begin  LW  05/06/2014; update metadata in File Header
  char startGranuleTime[100], stopGranuleTime[100];
  int granuleNumber;
//  end    LW  05/06/2014

//  begin  SFM  09/04/2013
  L2AKu_NS     L2AKuData;
//  L2AKu_NS_aa1     L2AKuData;
//  end  SFM  09/04/2013
 
  int inode=0, ii, nscans, rtype;
  char fname[1000];
  char emptygranuletext[100];
  int32 sd_id, sds_id, n_datasets, n_file_attrs,
           rain_index, geo_index, index, status32;
  float pia;
  int32 dim_sizes[10];
  int32 start[3];
  int32 rank, num_type, attributes, istat;
  int itop, ibb;
  int iSurf,iBinEllip,iClutf, ibb1, ibb2; 
 
  FILE *fout1,*fout2,*fout3;
  int isfzku;

  strcpy(&fname[0],&f2aku[0]);

//  SFM  begin  12/09/2013; to pass out TK status message
//  begin  SFM  09/04/2013
//  status = TKopen(&fname[0], "2AKu_aa1", TKREAD, "HDF5", jobname,
//		  &granuleHandle2AKu, 1);
  printf("%s\n",&fname[0]);
 
  status_alpha = TKopen(&fname[0], "2AKu", TKREAD, "HDF5", jobname,
		  &granuleHandle2AKu, 1);
  if (status_alpha != 0)
      {
	printf("WARNING: Unable to access 2AKu %i \n",status_alpha);
	return status_alpha ;
      }
  //  end  SFM  09/04/2013
  
  //  SFM  begin  12/17/2013; to avoid empty granule
  TKgetMetaString(&granuleHandle2AKu, "FileHeader", "EmptyGranule", 
		  emptygranuletext);
  TKsetMetaString(&ctkfile, "FileHeader", "EmptyGranule", emptygranuletext);
  //  SFM  end  12/17/2013
  
  //  begin  LW  05/06/2014; update metadata in File Header
  //  begin  WSO 10/01/2015; comment out following line to prevent doubling of output
  //         in InputRecord
  //    TKtransferMetaData(&granuleHandle2AKu, &ctkfile);
  //  end    WSO 10/01/2015
  TKgetMetaString(&granuleHandle2AKu, "FileHeader", "StartGranuleDateTime",
		  startGranuleTime);
  TKsetMetaString(&ctkfile, "FileHeader", "StartGranuleDateTime",
		  startGranuleTime);
  TKgetMetaString(&granuleHandle2AKu, "FileHeader", "StopGranuleDateTime",
		  stopGranuleTime);
  TKsetMetaString(&ctkfile, "FileHeader", "StopGranuleDateTime",
		  stopGranuleTime);
  TKgetMetaInt(&granuleHandle2AKu, "FileHeader", "GranuleNumber",
	       &granuleNumber);
  TKsetMetaInt(&ctkfile, "FileHeader", "GranuleNumber",
	       granuleNumber);
  //  end    LW  05/06/2014
  
  if (strncmp(emptygranuletext,"EMPTY",5) == 0)
    {
      status_alpha = -3;
      printf(" WARNING: 2AKu has empty granule \n",status);
      return status_alpha ;
    }
  

  //  SFM  end  12/17/2013
    
  i=0;
  izc=0;
  ipia=0;
  inodes=0;
  isfzku=0;
  int iland=0;
  
  while((TKendOfFile (&granuleHandle2AKu) != TK_EOF))
    {
      status=TKreadScan(&granuleHandle2AKu,&L2AKuData);
      
      for(j=0;j<49;j++)
	{  
	  iLandOcean[ipia]=L2AKuData.PRE.landSurfaceType[j];
	  flagPrecip[ipia]=L2AKuData.PRE.flagPrecip[j];
	  sfcRain[ipia]=L2AKuData.SLV.precipRateNearSurface[j];
	  xlon[ipia]=L2AKuData.Longitude[j];
	  xlat[ipia]=L2AKuData.Latitude[j];
	  ipia++;
	}
      i++;
    }

  *n1c21=i;

//  SFM  begin 12/09/2013; to pass out TK status message
  status_alpha = TKclose(&granuleHandle2AKu);
  return status_alpha ;
//  SFM  end   12/09/2013

}
