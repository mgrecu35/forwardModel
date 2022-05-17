#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf.h>
#include <mfhdf.h>
#include "TKheaders.h"
#include "TK_1CGMI.h"

extern TKINFO ctkfile;

// New routine  SFM  09/09/2013; transfers miscellaneous metadata to 
//              output file
// New routine  SFM  04/14/2014; added meta_mini

//  SFM  start  12/11/2013  add file status flags
//  SFM  start  09/25/2013

//  SFM  start  08/22/2014 ; modify cumulative status flag
int status_trans_ (int status_flag) ;
//  SFM  end    08/22/2014  

void meta_for_outputfile_(char* fSNOW, char* fSEAICE,
                          int *orbitNumber, int *rseed1, int *rseed2,
                          char* algorithmVersion, long *st_1, long *st_2, 
			  long *st_3, long *st_2akuenv, long *st_2aku, 
			  long *st_2adpr, long *st_snow, long *st_seaice)
//  SFM  end    09/25/2013
{

    long status_cumulative ;
//  SFM  end    12/11/2013 
    int  status,i;
    char ranstring[80], AlgRunInfo[1000], flag_string[8] ;

    printf(" meta transfer %s  \n",fSNOW) ;
    printf(" meta transfer %s  \n",fSEAICE) ;
    printf(" meta transfer %i  %i  %i  \n",*orbitNumber,*rseed1,*rseed2) ;
//  SFM  start  12/11/2013  add file status flags
    printf(" Status flags: Ku/DPR/GMIa/GMIb/GMIc/ENV/Snow/Sea: %i %i %i %i %i %i %i %i \n",
           *st_2aku, *st_2adpr, *st_1, *st_2, *st_3, *st_2akuenv, *st_snow,
	   *st_seaice ) ;

//  SFM  start  08/22/2014 ; modify cumulative status flag
//  reassign values from call sequence
    int int_st_2aku, int_st_2adpr, int_st_1, int_st_2, int_st_3 ;
    int int_st_2akuenv, int_st_snow, int_st_seaice ;
    int_st_2aku = *st_2aku ;
    int_st_2adpr = *st_2adpr ;
    int_st_1 = *st_1 ;
    int_st_2 = *st_2 ;
    int_st_3 = *st_3 ;
    int_st_2akuenv = *st_2akuenv ;
    int_st_snow = *st_snow ;
    int_st_seaice = *st_seaice ;

//  translate to new format flags
    int  new_flag ;
    status_cumulative = 0 ;

    new_flag = status_trans_ (int_st_2aku) ;
    status_cumulative = status_cumulative + new_flag * 10000000 ;
    new_flag = status_trans_ (int_st_2adpr) ;
    status_cumulative = status_cumulative + new_flag * 1000000 ;
    new_flag = status_trans_ (int_st_1) ;
    status_cumulative = status_cumulative + new_flag * 100000 ;
    new_flag = status_trans_ (int_st_2) ;
    status_cumulative = status_cumulative + new_flag * 10000 ;
    new_flag = status_trans_ (int_st_3) ;
    status_cumulative = status_cumulative + new_flag * 1000 ;
    new_flag = status_trans_ (int_st_2akuenv) ;
    status_cumulative = status_cumulative + new_flag * 100 ;
    //new_flag = status_trans_ (int_st_snow) ;
    new_flag = 7;
    status_cumulative = status_cumulative + new_flag * 10 ;
    //new_flag = status_trans_ (int_st_seaice) ;
    new_flag = 7;
    status_cumulative = status_cumulative + new_flag * 1 ;
    printf(" Cumulative status flag  %08i \n",status_cumulative) ;
//  SFM  end    08/22/2014  

//  SFM  end    12/11/2013 
    status = TKsetMetaInt (&ctkfile,"FileHeader", 
			  "GranuleNumber", *orbitNumber );

    strcpy(AlgRunInfo, "Random seeds: ");
    sprintf(ranstring, "%d %d", *rseed1, *rseed2);
    strcat(AlgRunInfo, ranstring);

//  SFM  begin  08/21/2014 ;  add protection so snow & seaice files aren't reported if not used
//  SFM  start  10/18/2013  for LW
    if (int_st_snow == 0) 
    {
        strcat(AlgRunInfo, "  Snow data: ");
        strcat(AlgRunInfo, fSNOW);
    }
//  SFM  end    10/18/2013  for LW
    if (int_st_seaice == 0) 
    {
        strcat(AlgRunInfo, "  Sea ice data: ");
        strcat(AlgRunInfo, fSEAICE);
    }
//  SFM  end    08/21/2014

//  SFM  start  12/11/2013  add file status flags
    strcat(AlgRunInfo, "  StatusCum: ");
    sprintf(ranstring, "%08i", status_cumulative);
    strcat(AlgRunInfo, ranstring);
//  SFM  end    12/11/2013 

    TKsetMetaString(&ctkfile, "AlgorithmRuntimeInfo", "", AlgRunInfo);

//  SFM  start  09/25/2013
    TKsetMetaString(&ctkfile, "FileHeader", "AlgorithmVersion",
                        algorithmVersion) ;
//  SFM  end    09/25/2013

}

//  SFM  begin    04/14/2014

void meta_mini_(char* algorithmVersion)
{
    TKsetMetaString(&ctkfile, "FileHeader", "AlgorithmVersion",
                        algorithmVersion) ;
//  SFM  end    04/14/2014
}

//  subroutine translates the output from file-read programs into
//  status flags to be written to cmb file header data

int status_trans_ (int status_flag)
    {
    int new_flag ;
    switch (status_flag)
        {
        //  file ok
        case 0 :
            new_flag = 0 ;
	    break ;
        //  file read error
        case 1 :
            new_flag = 7 ;
	    break ;
        //  intentional "nil" file found
        case 2 :
            new_flag = 8 ;
	    break ;
        //  empty granule
        case 3 :
            new_flag = 9 ;
	    break ;
        }
    return new_flag ;
    }
