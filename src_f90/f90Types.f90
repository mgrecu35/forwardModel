
module f90Types
   
  type stormStructType
     integer,dimension(:),pointer :: nodes     
     ! 5 nodes that define the storm structure 
     ! for stratiform, with BB, nodes[0] is the storm top, 
     ! nodes[1] is the BB top,
     ! nodes[2] is the BB peak, nodes[3] is the BB bottom, 
     ! and nodes[4] is the lowest
     ! clutter free gate (0<=nodes[0]<nodes[1]<nodes[2]
     ! <nodes[3]<nodes[4]<ngates)
     ! for convective, nodes[3] is not used
     ! the C convention is used the gates are 
     ! numbered from 0 to ngates-1
     integer   ::  iSurf      ! surface gate number
     real      ::  freezH     ! freezing height - not used
     integer   ::  rainType;  ! rainType -1 stratiform
  end type stormStructType

  type radarDataType
     integer                    :: ngates     ! number of gates
     real,dimension(:),pointer  :: z13obs     ! observed Ku-reflectivities
     real,dimension(:),pointer  :: z35obs     ! observed Ka-reflectivities
     real                       :: xlong      ! longitude -not used
     real                       :: xlat       ! latitude -not used
     real                       :: pia13srt   ! Ku-band SRT PIA 
     real                       :: relPia13srt   ! Ku-band SRT PIA 
     real                       :: pia35srt   ! Ka-band SRT PIA
     real                       :: relPia35srt   ! Ka-band SRT PIA
     real                       :: dr         ! gate size
     real,dimension(:),pointer  :: hh;        ! hh[i] is the height of gate [i]
     real                       :: hfreez
     real			:: sigmaZeroKu
     real			:: sigmaZeroKa
  end type radarDataType

  type  radarRetType
    
     integer   ::  ngates     ! number of gates
     integer   ::  nMemb      ! number of members
     integer   ::  nmfreq     ! # of simulated passive microwave frequencies
     real,dimension(:),pointer :: z13c        ! effective reflectivity 
                                              ! (attenuation corrected)
                                              ! at Ku-band
     real,dimension(:),pointer :: z35mod0     ! simulated observations at Ka-band
!  SFM  begin  06/16/2014; for M. Grecu, multiple scattering
     real,dimension(:),pointer :: dz35ms      !  multiple scattering effect
!  SFM  end    06/16/2014
     real,dimension(:),pointer :: z35         ! simulated observations at Ka-band
     real,dimension(:),pointer :: pwc         ! precipitation water content 
                                              ! (g/m3) 
     real,dimension(:),pointer :: rrate       ! rain rate (mm/h)
     real,dimension(:),pointer :: d0          ! median diameter (mm)
     real,dimension(:),pointer :: log10dNw    !
     real,dimension(:),pointer :: tb          ! simulated brightness temperatures
     real,dimension(:),pointer :: emTb          ! simulated brightness temperatures
     real,dimension(:),pointer :: emis          ! emissivity
     integer,dimension(:),pointer  ::  imu    ! mu index of the look up table
     integer,dimension(:),pointer  ::  iwc    
     integer,dimension(:),pointer  ::  icc        
                                              ! index of RH profile (from 1 to nc)
                                              ! nc is the number of possible 
                                              ! RH profiles see cloud.f90
     integer,dimension(:),pointer   ::  jcc       
                                              ! index of cloud profile (from 1 to nc)
     real,dimension(:),pointer      ::  sfc_wind, sfc_windU, sfc_windV   
                                              ! surface wind speed
     real, dimension(:),pointer     ::  pia13
     real, dimension(:),pointer     ::  pia35
     real, dimension(:),pointer     ::  simSigmaZeroKu
     real, dimension(:),pointer     ::  simSigmaZeroKa
     real, dimension(:),pointer     ::  z35mMean
     real, dimension(:),pointer     ::  z35mSig
     real                           ::  pia13mMean, pia35mMean, pia13mSig, pia35mSig
     integer                        ::  ntpw
     real, dimension(:),pointer     ::  tpw
     real, dimension(:),pointer     ::  tpwCldMod
     real, dimension(:),pointer     ::  logdNw
                                             ! the last four variables are not retrieved
                                             ! they are randomly set and 
                                             ! specify the conditions 
     ! in which the retrievals and 
     ! associated brightness temperatures are derived
  end type radarRetType

  type  retParamType
     
  !
  !  F=wz*SUM(zsim,ka-zobs,ka)**2+w13*(pia13-pia13srt)**2+
  !  w35*(pia35,sim-pia35srt)**2
  !  wz the weight of the reflectivity term
  !  w13 the weight of the pia squared difference at ku-band
  !  w35 the weight of the pia squared difference at ka-band
  !  z13thresh - threshold (dBZ) to determine 
  !  the Ku-band observations used in the algorithm
  !  z35thresh - threshold (dBZ) to determine the Ka-band observations 
  !  used in the algorithm
  
  
     real :: wz, w13,  w35,  z13thresh,  z35thresh
  end type retParamType


end module f90Types
