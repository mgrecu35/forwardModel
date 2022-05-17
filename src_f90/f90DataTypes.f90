!  SFM  04/06/2013  Added Lv2AKuENV_DataType

module f90DataTypes
  !...Level 2AKuENV Ku Environment File Data
  !---Programmer:  S. McLaughlin   9/13/2012

  TYPE Lv2AKuENV_DataType
                     ! ***  Allocation sizes, internal storage *** 
    INTEGER*4 nscan_2akuenv  ! Number scans in granule
    INTEGER*4 nray2akuenv    ! Number angle bins in each NS scan
    INTEGER*4 nbin2akuenv    ! Number range bins in each NS and MS ray
    INTEGER*4 nwater         ! TBD
    INTEGER*4 nwind          ! Number of wind components: u,v

    REAL*8,    POINTER :: secondOfDay(:)    ! Seconds of day, time associated w 
                                            !   scan, 0 to 86400
    INTEGER*4, POINTER :: dayofyear(:)      ! Day of year, DDD format, 1 to 365
    INTEGER*2, POINTER :: year(:)           ! Year, YYYY format, 1950 to 2100
    INTEGER*2, POINTER :: millisecond(:)    ! Thousandths of second, MMM, 0 to 999
    INTEGER*1, POINTER :: month(:)          ! Month of year, MM format, 1 to 12
    INTEGER*1, POINTER :: dayofmonth(:)     ! Day of month, DD format, 1 to 31
                                            !   UTC time
    INTEGER*1, POINTER :: hour(:)           ! Hour of Day, HH format, 0 to 23
    INTEGER*1, POINTER :: minute(:)         ! Minute of hour, MM format, 0 to 59
    INTEGER*1, POINTER :: second(:)         ! Second of minute, SS format, 0 to 60
    REAL*4,    POINTER :: latitude(:,:)     ! Eath lat @ center FOV @ ellipsoid
                                            !   -90 to 90
    REAL*4,    POINTER :: longitude(:,:)    ! Eath lon @ center FOV @ ellipsoid
                                            !   -180 to 180
    REAL*4,  POINTER :: airTemp(:,:,:)      ! air temperature, [K]
    REAL*4,  POINTER :: airPressure(:,:,:)  ! air pressure, [hPa]
    REAL*4,  POINTER :: waterVapor(:,:,:,:) ! water vapor, [kg/m**3]
    REAL*4,  POINTER :: cldLiqWat(:,:,:,:)  ! cloud liquid water, [kg/m**3]
!begin  WSO 9/1/13 add skin temperature
    REAL*4,  POINTER :: skinTemperature(:,:)! skin temperature, [K]
!end    WSO 9/1/13
    REAL*4,  POINTER :: surfPressure(:,:)   ! surface pressure, [hPa]
    REAL*4,  POINTER :: groundTemp(:,:)     ! ground temperature, [K]
    REAL*4,  POINTER :: surfWind(:,:,:)     ! surface wind, [m/s]	    
  END TYPE Lv2AKuENV_DataType

  type geoDataType
     
     byte:: lsflag(1080,2160)
     integer*2 :: sstdata(91,144,12)
     integer   :: month
  end type geoDataType
  type dPRDataType
     integer :: n1c21                                   ! number of scans
     integer :: nmax                                    ! max number of scans
     integer :: nrays                                   ! max number of scans
     integer :: ngates                                  ! max number of scans
     real  :: scAngle(49)                               ! DPR view angle
     real, dimension(:,:,:), pointer    :: zku1c21      ! ku-band reflectivity
     real, dimension(:,:,:), pointer    :: zka1c21      ! ka-band reflectivity
     real, dimension(:,:), pointer      :: xlon,xlat    ! longitude and latitude
     real, dimension(:), pointer 	:: sclon, sclat ! sub-satellite lat/lon SJM 8/4/16
!begin  WSO 9/28/13 add bad ray flag
     integer, dimension(:,:), pointer   :: badRayFlag   ! indicates unacceptable
                                                        ! data quality
     integer, dimension(:,:), pointer   :: rainFlagBad  ! precip flag with bad ray info
!end    WSO 9/28/13
!begin  WSO 9/10/13 add snr ratios
     real, dimension(:,:), pointer      :: snrRatioku   ! surface SNR at Ku-band
     real, dimension(:,:), pointer      :: snrRatioka   ! surface SNR at Ka-band
!end    WSO 9/10/13
!begin SJM 7/9/14 add sigma_zero
     real, dimension(:,:), pointer      :: sigmaZeroKu   ! observed sigma_zero at Ku-band
     real, dimension(:,:), pointer      :: sigmaZeroKa   ! observed sigma_zero at Ka-band
     !real, dimension(:,:), pointer      :: dummy
!end SJM 7/9/14
!begin  WSO 9/5/13 add dsrt PIA's, relability factors
     real, dimension(:,:), pointer      :: srtPIAku     ! SRT PIA at Ku-band
     real, dimension(:,:), pointer      :: dsrtPIAku    ! DSRT PIA at Ku-band
     real, dimension(:,:), pointer      :: dsrtPIAka    ! DSRT PIA at Ka-band
     real, dimension(:,:), pointer      :: srtrelPIAku  ! SRT PIA at Ku-band
     real, dimension(:,:), pointer      :: dsrtrelPIA 
     integer, dimension(:,:), pointer      :: NSrelibFlag
     integer, dimension(:,:), pointer      :: MSrelibFlag

!end    WSO 9/5/13
!begin  WSO 9/5/13 add composite sigma
     real, dimension(:,:), pointer      :: srtsigmaPIAku   ! uncertainty of SRT PIA at Ku-band
     real, dimension(:,:), pointer      :: dsrtsigmaPIAku   ! uncertainty of DSRT PIA at Ku-band
     real, dimension(:,:), pointer      :: dsrtsigmaPIAka   ! uncertainty of DSRT PIA at Ka-band
!end    WSO 9/5/13 
!begin  WSO 8/15/14 add ioquality flag variables
     integer, dimension(:,:), pointer   :: ioqualityflagku  ! ioquality flag for NS (Ku+GMI) product
     integer, dimension(:,:), pointer   :: ioqualityflagdpr ! ioquality flag for MS (Ku+Ka+GMI) product
!end    WSO 8/15/14
     real, dimension(:,:), pointer      :: piaHB        ! SRT PIA at Ka-band
     integer, dimension(:,:,:), pointer :: node         ! storm struct
     integer, dimension(:,:), pointer   :: rainType     ! rain type
     integer, dimension(:,:), pointer   :: BBbin        ! rain type
     real, dimension(:,:), pointer      :: freezH          ! freezingH
     integer, dimension(:,:), pointer   :: binRealSurface
     real, dimension(:,:), pointer      :: localZenithAngle
     real, dimension(:,:), pointer      :: elevation
     integer, dimension(:,:), pointer   :: ig, jg
     real, dimension(:,:),  pointer   :: surfaceZku
     real, dimension(:,:),  pointer   :: surfaceZka
     integer, dimension(:,:), pointer   :: iLandOcean
     real, dimension(:,:,:), pointer :: envQv
     real, dimension(:,:,:), pointer :: envPress
     real, dimension(:,:,:), pointer :: envTemp
     real, dimension(:,:,:), pointer :: envCloud
     real, dimension(:,:), pointer   :: envSfcWind
     real, dimension(:,:), pointer   :: envSfcWindU !SJM 8/4/16
     real, dimension(:,:), pointer   :: envSfcWindV !SJM 8/4/16
     integer, dimension(:,:), pointer   :: snowIceCover !SJM 11/16/16
     real, dimension(:,:), pointer   :: seaIceConcentration !SJM 11/16/16
!begin  WSO 9/1/13 add skin temperature
     real, dimension(:,:), pointer   :: envSknTemp
!end    WSO 9/1/13
     real, dimension(:,:), pointer   :: envSfcTemp
     real, dimension(:,:), pointer   :: envSfcPress
     real, dimension(:,:), pointer   :: cBest
     real, dimension(:), pointer :: SecondOfDay ! 4/14/14 MG
  end type dPRDataType

  type gMIDataType
     integer :: n1b11, mm
     integer :: n1b11H
     integer :: nmax
     integer :: nlowf, nhf, nrays
 
     real,pointer :: gmilow(:,:,:), gmihigh(:,:,:), &
          txlon(:,:),  txlat(:,:)
     real,pointer :: sst(:,:), tpw(:,:), lwp(:,:), sfc_wind(:,:)
     real,pointer :: gmilon(:), gmilat(:)
     real,pointer :: sScan(:), eScan(:)
     byte, pointer :: landSea(:,:), rainFlag(:,:)
     byte, pointer :: rainFlagDPR(:,:)
     real, pointer :: tb10hC(:,:), tb10vC(:,:), tb19hC(:,:), tb19vC(:,:), &
          tb21vC(:,:),tb37hC(:,:), tb37vC(:,:)
     real, pointer :: tb10hC3(:,:), tb10vC3(:,:), tb19hC3(:,:), tb19vC3(:,:), &
          tb21vC3(:,:),tb37hC3(:,:), tb37vC3(:,:)
     real,pointer  ::  gmilow3(:,:,:), gmilow3H(:,:,:), &
          txlon3(:,:), txlat3(:,:),gmihigh3(:,:,:), &
          sst3(:,:), tpw3(:,:), lwp3(:,:), sfc_wind3(:,:) 
     real,pointer  ::  simTb(:,:,:,:,:), simTbMax(:,:,:,:), &
          simTbMin(:,:,:,:), simTbMean(:,:,:,:), simTbStd(:,:,:,:)

     real,pointer :: dPRsfcRain(:,:)
     real,pointer :: dPRsfcRainEns(:,:,:)
     real,pointer :: tpwEns(:,:,:)

  end type gMIDataType

 type cGMIDataType
     integer :: n1b11, mm
     integer :: n1b11H
     integer :: nmax
     integer :: nS1f, nS2f, nrays
     real,pointer :: gmiS1(:,:,:), gmiS2(:,:,:), &
          S1lon(:,:),  S1lat(:,:), S2lon(:,:),  S2lat(:,:), &
          emissS1(:,:,:), emissS2(:,:,:)
     real, pointer :: sst(:,:)
     integer*8, pointer :: landSea(:,:), rainFlag(:,:)
     real, pointer :: gmiLon(:), gmiLat(:)
     real, pointer :: tb10hC(:,:), tb10vC(:,:), tb19hC(:,:), tb19vC(:,:), &
          tb21vC(:,:),tb37hC(:,:), tb37vC(:,:),  tpw(:,:), lwp(:,:), &
          sfc_wind(:,:)

     real,pointer :: gmiS13(:,:,:), gmiS23(:,:,:)
     real, pointer :: tb10hC3(:,:), tb10vC3(:,:), tb19hC3(:,:), tb19vC3(:,:), &
          tb21vC3(:,:),tb37hC3(:,:), tb37vC3(:,:),  S1lon3(:,:),  S1lat3(:,:), sst3(:,:), &
          sfc_wind3(:,:), tpw3(:,:), gmiS13H(:,:,:), SCLon3(:), SCLat3(:), S1eia3(:,:), S2eia3(:,:) !SJM 8/4/16
     real, pointer :: emissS13(:,:,:)
     logical, pointer :: rainFlagDPR(:,:)
     real, pointer :: SecondOfDay(:) !4/14/14 MG
     real, pointer :: SCLon(:) !4/14/14 MG
     real, pointer :: SCLat(:) !4/14/14 MG
     real, pointer :: S1eia(:,:) !SJM 8/4/16
     real, pointer :: S2eia(:,:) !SJM 8/4/16
  end type cGMIDataType

  type gmi2GridType
     real         :: xmin, ymin, dx
     integer      :: nx, ny
     real,allocatable :: tpw(:,:)
     integer,allocatable :: actOb(:,:)
     integer,allocatable :: ig(:,:), jg(:,:)
  end type gmi2GridType

  type gridDataType
     integer:: nxg, nyg, ig, jg
     integer:: nMemb
     real:: dx,xmin,xmax,ymin,ymax 
     integer, pointer :: indigPR(:,:), &
       indjgPR(:,:), indiPR(:,:), indjPR(:,:)
  end type gridDataType

  type gridEnvDataType
     integer:: nxg, nyg
     integer:: nMemb, nmfreq
     real:: dx,xmin,xmax,ymin,ymax 
     real, dimension(:,:,:,:,:), pointer :: tb
     real, dimension(:,:,:), pointer :: sfc_wind
     integer, dimension(:,:,:), pointer :: icc
     integer, dimension(:,:,:), pointer :: jcc
     integer, dimension(:,:), pointer ::   yesNo
  end type gridEnvDataType

  type MSRetType
     real,dimension(:,:,:,:),pointer :: pwc       ! precipitation water content 
     ! (g/m3) 
     real,dimension(:,:,:,:),pointer :: rrate     ! rain rate (mm/h)
     real,dimension(:,:), pointer    :: sfcRain     !
     real,dimension(:,:,:), pointer    :: pia13mod  !
     real,dimension(:,:,:), pointer    :: pia35mod  !
     real,dimension(:,:,:), pointer    :: sfcRainEns ! 
     real,dimension(:,:,:), pointer    :: sfcd0Ens   !
     real,dimension(:,:,:), pointer    :: sfcNwEns   ! 
     real,dimension(:,:,:), pointer    :: sfcWindEns !SJM 12/4/2014 
     real,dimension(:,:,:,:),pointer :: d0           ! median diameter (mm)
     real,dimension(:,:,:,:),pointer :: log10dNw     ! 
     real,dimension(:,:,:,:,:),pointer :: tb         ! simulated brightness temperatures
     real,dimension(:,:,:,:,:),pointer :: emTb         ! simulated brightness temperatures
     real,dimension(:,:,:,:,:),pointer :: emis         ! simulated emissivity !SJM 8/4/16
     real,dimension(:,:,:,:,:),pointer :: convtb         ! simulated brightness temperatures
     real,dimension(:,:,:,:),pointer :: zkuEns         
     real,dimension(:,:,:,:),pointer :: zkaEns
     real, dimension(:,:,:), pointer :: simSigmaZeroKu !Simulate Sigma-Zero at Ku band
     real, dimension(:,:,:), pointer :: simSigmaZeroKa !Simulated Sigma-Zero at Ka band
  end type MSRetType

  type dPRRetType
     type (MSRetType) :: MS
     integer :: n1c21                             ! number of scans
     integer :: nMemb
     integer :: ngates
     integer :: nmfreq
     real,dimension(:,:,:,:),pointer :: z13c      ! effective reflectivity 
     ! (attenuation corrected)
     ! at Ku-band
     real,dimension(:,:,:,:),pointer :: z35mod0   ! simulated observations at Ka-band
!  SFM  start  06/16/2014;  for M.Grecu, multiple scattering
     real,dimension(:,:,:,:),pointer :: dz35ms   ! simulated observations at Ka-band
!  SFM  start  06/16/2014;
     real,dimension(:,:,:,:),pointer :: z35   ! simulated observations at Ka-band
     real,dimension(:,:,:,:),pointer :: pwc       ! precipitation water content 
                                                  ! (g/m3) 
     real,dimension(:,:,:,:),pointer :: rrate     ! rain rate (mm/h)
     real,dimension(:,:), pointer    :: sfcRain     !
     real,dimension(:,:), pointer    :: sfcRainRef  !
     real,dimension(:,:,:), pointer    :: pia13mod  !
     real,dimension(:,:,:), pointer    :: pia35mod  !
     real,dimension(:,:,:), pointer    :: sfcRainEns ! 
     real,dimension(:,:,:), pointer    :: sfcd0Ens   !
     real,dimension(:,:,:), pointer    :: sfcNwEns   ! 
     real,dimension(:,:,:), pointer    :: sfcWindEns !SJM 12/4/2014 
     real,dimension(:,:,:,:),pointer :: d0           ! median diameter (mm)
     real,dimension(:,:,:,:),pointer :: log10dNw     ! 
     integer, dimension(:,:,:),pointer :: n9         ! 
     real,dimension(:,:,:,:,:),pointer :: tb         ! simulated brightness temperatures
     real,dimension(:,:,:,:,:),pointer :: emTb         ! simulated brightness temperatures
     real,dimension(:,:,:,:,:),pointer :: emis         ! simulated emissivity SJM 8/4/16
     real,dimension(:,:,:,:,:),pointer :: convtb         ! simulated brightness temperatures
     !begin SJM 7/9/2014 add sigma_zero simulated variables
     real, dimension(:,:,:), pointer :: simSigmaZeroKu !Simulate Sigma-Zero at Ku band
     real, dimension(:,:,:), pointer :: simSigmaZeroKa !Simulated Sigma-Zero at Ka band
     !end SJM 7/9/2014
     integer,dimension(:,:,:),pointer  :: imu        ! mu index of the look up table
     integer,dimension(:,:,:),pointer  :: icc        
     real, dimension(:,:,:,:),pointer  :: cldwcoeff
     ! index of RH profile (from 1 to nc)
     ! nc is the number of possible 
     ! RH profiles see cloud.f90
     integer,dimension(:,:,:),pointer   ::  jcc       
     ! index of cloud profile (from 1 to nc)
     real,dimension(:),pointer      ::  sfc_wind
     real,dimension(:),pointer      ::  tpw
     real,dimension(:,:,:), pointer      ::  tpwEns
!     real,dimension(:,:,:,:), pointer :: tbMin
!     real,dimension(:,:,:,:), pointer :: tbMax
!     real,dimension(:,:,:,:), pointer :: tbStd
!     real,dimension(:,:,:,:), pointer :: tbMean
!     real,dimension(:,:,:,:), pointer :: tbCMean
!     real,dimension(:,:,:,:), pointer :: tbSim
!     real,dimension(:,:,:), pointer ::   obsSig
!     real,dimension(:,:,:), pointer ::   simSig
   


     ! surface wind speed
  end type dPRRetType
end module f90DataTypes

