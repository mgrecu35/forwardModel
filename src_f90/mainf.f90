!! SFM 05/06/2013  Modifications from LW to facilitate using job names
!! SFM 07/19/2013  Modifications from ngniM.Grecu for convergence problem
!! SFM 08/09/2013  Eliminate IF2ADPROUT dependences; few extra diagnostics
!!
!microwave table freq 19.0 37. 10. 85
!TMIchannel freq  10.0, 19.0, 21, 37, 85! V,H
!npol 0 H, npol 1 V


SUBROUTINE init_random_seed(rseed1, rseed2)
  INTEGER :: rseed1, rseed2, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed

  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))

  if(rseed1.le.-1) then
    ! Don't set seed - what delivered code does.
  elseif(rseed1.eq.0) then
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)
  else
    seed = rseed1 + 37 * (/ (i - 1, i = 1, n) /)
    if(rseed2.gt.0) then
      seed(2) = rseed2
    endif
    CALL RANDOM_SEED(PUT = seed)
  endif

  CALL RANDOM_SEED(GET = seed)
  rseed1 = seed(1)
  rseed2 = seed(2)
  DEALLOCATE(seed)
END SUBROUTINE init_random_seed

SUBROUTINE init_random_seed2()
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  
  CALL SYSTEM_CLOCK(COUNT=clock)
  
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)
  
  DEALLOCATE(seed)
END SUBROUTINE init_random_seed2

subroutine mainfort(jobname, f1ctmi1,f1ctmi2,f1ctmi3,                    &
                    f2AKu, f2aDPR, f2AkuENV, fSNOW, fSEAICE, f2CMB,      &
                    rseed1, rseed2, igmi1, igmi2, igmi3, i2AKu,          &
                    i2ADPR, i2AkuENV, iSNOW, iSEAICE, i2CMB, ialg, ndpr1, ifs)

  use iso_c_binding    !iso c binding statement
  use globalData
  !use gfsmodel        ! SFM  04/16/2014  deleted
  use emissmod
  use nbinmod
!begin  WSO 9/14/13
  use missingMod
!end    WSO 9/14/13
!begin WSO 8/7/13
  use Tables_frac
!  use writeENKF
!end WSO 8/7/13

  implicit none
  integer :: rseed1, rseed2, ifs
  real :: tbRgrid(14,49,9300)  ! resampled brightness temperatures on same grid 
                              !   as 2adpr data; (channel/ray/scan)


!  SFM  04/06/2013  Changed file name lengths to 1000
  character(c_char) :: jobname(255), f1ctmi1(1000), f1ctmi2(1000), &
       f1ctmi3(1000), f2AKu(1000), f2aDPR(1000), f2AKuENV(1000), f2CMB(1000)
  character(c_char) :: fSNOW(1000), fSEAICE(1000)
  character(c_char) :: outcmb(1000) !iso c binding

  integer :: igmi1, igmi2, igmi3, i2AKu, i2ADPR, i2AKuENV, i2CMB, ialg
  integer :: iSNOW, iSEAICE
  integer :: i1ctmi, i1c21, i2a23
  real    :: mu
  integer :: nGMIS1, nGMIS2, nS1f, nS2f, nchunks
  character*3 :: ifdpr
  INTEGER*4 :: orbitNumber    ! retrieved orbit number
  integer :: ngmi1,ngmi2,ngmi3, j
  integer :: isnow1, iseaice1, ndpr1
  real    :: lastGood          !  SFM  04/16/2014   added for M.Grecu
  
  character*90 :: lut_file !SJM 7/9/2014

!  SFM  start  09/25/2013
  character(16), parameter :: algorithmVersion = "2BCMB_20180919"
  !  SFM  end  09/25/2013
!  SFM  start  09/27/2013
!...Variables required for implementation of autosnow option
!  SFM  start  10/25/2013  for LW
    integer    :: date(3)=(/0,0,0/), istat,    kk1, kk2
!  SFM  end    10/25/2013  for LW
    integer(1) :: autosnow(9000,4500)                 !4km global grid
!  SFM  end    09/27/2013

!...Allocate storage space and array size parameters for new data files

    TYPE (Lv2AKuENV_DataType)  :: Lv2AKuENV_scan
    INTEGER :: nscan_2akuenv, nray2akuenv, nbin2akuenv, nwater, nwind

!  SFM  start  01/02/2014
   
    INTEGER*4     date_number      ! file date in format yyyymmdd
    INTEGER*4     month_get        ! computed file month
    INTEGER*4     st_1, st_2, st_3 ! status codes for 1CGMI files
 
    INTEGER*4     readgmi, readtmi         ! declare function
    INTEGER*4     readenv, readenvx          ! declare function
    INTEGER*4     readdprpflag, readdprpflagx     ! declare function
    !LW 05/04/18
    INTEGER*4     readdprtpflag     ! declare function

    integer*4     i, ifract, idir
    real :: x1L, x2L
!  SFM  end    01/02/2014
  
!  SFM  04/06/2013  Added code and supporting datasets for query_2akuenv, 
!                     allocate_2AKuENV_space, and read_2kauenv 
!  SFM  06/19/2013  Code from M.Grecu integrated

  rec=0
  
  mu=cos(52.8/180*3.14159)
  do i=1,7
     !call mie2(i-3.)
  enddo
  
  nmfreq=8
  nmu=5
  call init_random_seed(rseed1,rseed2)
  print *, 'Random Seeds : ',rseed1,rseed2
  print*, ialg
  print*, f2AKu(1:i2AKu) 
  print*, f2aDPR(1:i2adpr)
  !stop
  call readtablesLiang2(nmu,nmfreq)         !This option reads in Liang's mu=2 table
  !call readtablesDSDWG(nmu,nmfreq) !Replaces Liang's DSD table with DSD working group relationships (rain only)
  
  call makeHashTables()
  call cloud_init(nmfreq) 
  call initWFlag(nmfreq)
  call init_nbin
!begin  WSO 9/14/13 initial missing flags in missingMod
  call init_missing_flags
!end    WSO 9/14/13

!begin WSO 8/7/13
  call read_melt_percentages
!end WSO 8/7/13

  !begin SJM 7/9/2014 read emissivity and sigma_zero LUTs
  lut_file = 'Emiss/LUT.GMI.emis.MW_ITE101.bin'
  call read_LUTwatemis(lut_file)
  lut_file = 'Emiss/LUT.DPR.sigma0.MWopt.GANAL_V5.bin'
  call read_lutwatsigma0(lut_file)
  !end SJM 7/9/2014
  !read in emissivity-sigma0 EOFs
  lut_file = 'GANAL_ITE030.V04'
  call read_LUTlandclass(lut_file)
  !end SJM 7/9/2014
  nmF=4
  if(.not.(allocated(ip))) &      
       allocate(ip(nMF), iGMIChan(nMF), iSimF(nMF), iGMIf(nMF))
  
  nPrScanM = 300; ngates=nbin; ndPRrays=49;
  nGMIMax  = 3700; nGMIS1=9; nGMIS2=4; nGMIrays=221;
  gridData%dx=0.05
!  call readdmnw()

!  SFM  start  11/18/2013
  nMemb=50
!  SFM  end    11/18/2013
  ifdpr(1:1)='N'
  
!...Check file availability for processing
    !print*, 
    CALL file_check(igmi1, igmi2, igmi3, i2aku, i2adpr, i2akuenv, isnow, iseaice,  &
                    i2cmb, f1ctmi1(1:igmi1), f1ctmi2(1:igmi2), f1ctmi3(1:igmi3),   &
		    f2adpr(1:i2adpr), f2akuenv(1:i2akuenv), f2aku(1:i2aku),        &
		    fsnow(1:isnow), fseaice(1:iseaice), f2cmb(1:i2cmb),            &
		    st_1, st_2, st_3, st_2akuenv, st_2adpr, st_2aku, st_snow,      &
		    st_seaice, st_cmb, date_number, orbitNumber)

!  SFM  start  01/02/2014  reworked checks for nil and dead files

!  Open output file
  IF (st_cmb .NE. 0) THEN
      IF (st_cmb .EQ. -1) PRINT *,'WARNING: 2CMB directory NOT available'
      IF (st_cmb .EQ. -2) PRINT *,'WARNING: 2CMB file designated nil'
      RETURN
  ELSE
      outcmb(1:i2cmb+1)=f2CMB(1:i2cmb)//char(0)
      if(ialg==2) then
         CALL openoutputfilet(jobname,outcmb)
      else
         if(ifs==1) then
            CALL openoutputfile_fs(jobname,outcmb)
         else
            CALL openoutputfile(jobname,outcmb)
         endif
      endif
  ENDIF
 
!  SFM  start  04/10/2014  guarantee metadata on early exit
!.Look for Cases of Bad 2Aku
  IF (st_2aku .LT. 0) THEN
      IF (st_2aku .EQ. -1) PRINT *,'WARNING: 2AKu file is NOT available'
      IF (st_2aku .EQ. -2) PRINT *,'WARNING: 2AKu file is "nil" '
      IF (st_2aku .EQ. -4) PRINT *,'WARNING: 2AKu file is NOT HDF format'

      !LW 05/03/18
      !CALL write_empty()
      if(ialg==2) then
         CALL write_emptyt()
      else
         CALL write_empty()
      endif

      CALL meta_mini(algorithmVersion)
      CALL closeoutputfilet()
      RETURN
  ENDIF
  
  !.Initialize the 2Aku Data
  PRINT *,'File read status 2AKu  : ',st_2aku

  if(ialg==2) then
    st_2aku = readdprtpflag(jobname, f2AKu(1:i2AKu),ndpr,rlon,rlat, &
             iLandOcean, rPrecipFlag, rSfcPrecip)
  else
     if(ifs==1) then
        st_2aku = readdprpflagx(jobname, f2AKu(1:i2AKu),ndpr,rlon,rlat, &
             iLandOcean, rPrecipFlag, rSfcPrecip)
     else
        st_2aku = readdprpflag(jobname, f2AKu(1:i2AKu),ndpr,rlon,rlat, &
             iLandOcean, rPrecipFlag, rSfcPrecip)
     end if

  endif

!.Process Case of Empty 2Aku Granule
  IF (st_2aku .EQ. -3) THEN
      PRINT *,'WARNING: 2AKu file is empty granule'

      !LW 05/03/18
      !CALL write_empty()
      if(ialg==2) then
         CALL write_emptyt()
      else
         CALL write_empty()
      endif

      CALL meta_mini(algorithmVersion)
      CALL closeoutputfilet()
      RETURN
  ENDIF
!  SFM  end    01/02/2014
!  SFM  end    04/10/2014

!  iwENKF=0
  call read_geodat(geoData%lsflag, geoData%sstdata)

!  SFM  start  12/06/2013  reworked checks for nil and dead files
  nS1f=9
  nS2f=4
  ngmi_total = 0
  
  call allocatecGMISpace(gMIData, nS1f, nS2f, nGMIrays, nGMIMax, nmfreq, nmemb)

  ngmi1=0
  if(ialg==2) then
     IF (st_1 .GE. 0) THEN
        st_1=readtmi(jobname,f1ctmi1(1:igmi1),ngmi1,gMIData%gmiS1(:,:,ngmi_total+1:ngmimax), &
             gMIData%gmiS2(:,:,ngmi_total+1:ngmimax),gMIData%S1lon(:,ngmi_total+1:ngmimax),  &
             gMIData%S1lat(:,ngmi_total+1:ngmimax),                                          &
             gMIData%gmilon(ngmi_total+1:ngmimax), gMIData%gmilat(ngmi_total+1:ngmimax),     &
             gmiData%mm,year,jday,dd,1,gmiData%secondOfDay(ngmi_total+1:ngmimax),&
             gmiData%SCLon(ngmi_total+1:ngmimax),&
             gmiData%SCLat(ngmi_total+1:ngmimax),&
             gMIData%S2lon(:,ngmi_total+1:ngmimax),  &
             gMIData%S2lat(:,ngmi_total+1:ngmimax), &
             GMIdata%S1eia(:,ngmi_total+1:ngmimax), &
             GMIdata%S2eia(:,ngmi_total+1:ngmimax)) !SJM 3/15/16
     ENDIF
  !stop
     ngmi_total=ngmi1+ngmi_total

     ngmi2=0
     IF (st_2 .GE. 0) THEN
        st_2=readtmi(jobname,f1ctmi2(1:igmi2),ngmi2,gMIData%gmiS1(:,:,ngmi_total+1:ngmimax), &
             gMIData%gmiS2(:,:,ngmi_total+1:ngmimax),gMIData%S1lon(:,ngmi_total+1:ngmimax),  &
             gMIData%S1lat(:,ngmi_total+1:ngmimax),                                          &
             gMIData%gmilon(ngmi_total+1:ngmimax), gMIData%gmilat(ngmi_total+1:ngmimax),     &
             gmiData%mm,year,jday,dd,2,gmiData%secondOfDay(ngmi_total+1:ngmimax),&
             gmiData%SCLon(ngmi_total+1:ngmimax),&
             gmiData%SCLat(ngmi_total+1:ngmimax),&
             gMIData%S2lon(:,ngmi_total+1:ngmimax),  &
             gMIData%S2lat(:,ngmi_total+1:ngmimax), &
             GMIdata%S1eia(:,ngmi_total+1:ngmimax), &
             GMIdata%S2eia(:,ngmi_total+1:ngmimax)) !SJM 3/15/16
     ENDIF
     ngmi_total=ngmi2+ngmi_total
     
     ngmi3=0
     IF (st_3 .GE. 0) THEN
        st_3=readtmi(jobname,f1ctmi3(1:igmi2),ngmi3,gMIData%gmiS1(:,:,ngmi_total+1:ngmimax), &
             gMIData%gmiS2(:,:,ngmi_total+1:ngmimax),gMIData%S1lon(:,ngmi_total+1:ngmimax),  &
             gMIData%S1lat(:,ngmi_total+1:ngmimax),                                          &
             gMIData%gmilon(ngmi_total+1:ngmimax), gMIData%gmilat(ngmi_total+1:ngmimax),     &
             gmiData%mm,year,jday,dd,3,gmiData%secondOfDay(ngmi_total+1:ngmimax),&
             gmiData%SCLon(ngmi_total+1:ngmimax),&
             gmiData%SCLat(ngmi_total+1:ngmimax),&
             gMIData%S2lon(:,ngmi_total+1:ngmimax),  &
             gMIData%S2lat(:,ngmi_total+1:ngmimax), &
             GMIdata%S1eia(:,ngmi_total+1:ngmimax), &
             GMIdata%S2eia(:,ngmi_total+1:ngmimax)) ! SJM
     ENDIF
  else
     IF (st_1 .GE. 0) THEN
        st_1=readgmi(jobname,f1ctmi1(1:igmi1),ngmi1,gMIData%gmiS1(:,:,ngmi_total+1:ngmimax), &
             gMIData%gmiS2(:,:,ngmi_total+1:ngmimax),gMIData%S1lon(:,ngmi_total+1:ngmimax),  &
             gMIData%S1lat(:,ngmi_total+1:ngmimax),                                          &
             gMIData%gmilon(ngmi_total+1:ngmimax), gMIData%gmilat(ngmi_total+1:ngmimax),     &
             gmiData%mm,year,jday,dd,1,gmiData%secondOfDay(ngmi_total+1:ngmimax),&
             gmiData%SCLon(ngmi_total+1:ngmimax),&
             gmiData%SCLat(ngmi_total+1:ngmimax),&
             gMIData%S2lon(:,ngmi_total+1:ngmimax),  &
             gMIData%S2lat(:,ngmi_total+1:ngmimax), &
             GMIdata%S1eia(:,ngmi_total+1:ngmimax), &
             GMIdata%S2eia(:,ngmi_total+1:ngmimax)) !SJM 3/15/16
     ENDIF
  !stop
     ngmi_total=ngmi1+ngmi_total

     ngmi2=0
     IF (st_2 .GE. 0) THEN
        st_2=readgmi(jobname,f1ctmi2(1:igmi2),ngmi2,gMIData%gmiS1(:,:,ngmi_total+1:ngmimax), &
             gMIData%gmiS2(:,:,ngmi_total+1:ngmimax),gMIData%S1lon(:,ngmi_total+1:ngmimax),  &
             gMIData%S1lat(:,ngmi_total+1:ngmimax),                                          &
             gMIData%gmilon(ngmi_total+1:ngmimax), gMIData%gmilat(ngmi_total+1:ngmimax),     &
             gmiData%mm,year,jday,dd,2,gmiData%secondOfDay(ngmi_total+1:ngmimax),&
             gmiData%SCLon(ngmi_total+1:ngmimax),&
             gmiData%SCLat(ngmi_total+1:ngmimax),&
             gMIData%S2lon(:,ngmi_total+1:ngmimax),  &
             gMIData%S2lat(:,ngmi_total+1:ngmimax), &
             GMIdata%S1eia(:,ngmi_total+1:ngmimax), &
             GMIdata%S2eia(:,ngmi_total+1:ngmimax)) !SJM 3/15/16
     ENDIF
     ngmi_total=ngmi2+ngmi_total
     
     ngmi3=0
     IF (st_3 .GE. 0) THEN
        st_3=readgmi(jobname,f1ctmi3(1:igmi2),ngmi3,gMIData%gmiS1(:,:,ngmi_total+1:ngmimax), &
             gMIData%gmiS2(:,:,ngmi_total+1:ngmimax),gMIData%S1lon(:,ngmi_total+1:ngmimax),  &
             gMIData%S1lat(:,ngmi_total+1:ngmimax),                                          &
             gMIData%gmilon(ngmi_total+1:ngmimax), gMIData%gmilat(ngmi_total+1:ngmimax),     &
             gmiData%mm,year,jday,dd,3,gmiData%secondOfDay(ngmi_total+1:ngmimax),&
             gmiData%SCLon(ngmi_total+1:ngmimax),&
             gmiData%SCLat(ngmi_total+1:ngmimax),&
             gMIData%S2lon(:,ngmi_total+1:ngmimax),  &
             gMIData%S2lat(:,ngmi_total+1:ngmimax), &
             GMIdata%S1eia(:,ngmi_total+1:ngmimax), &
             GMIdata%S2eia(:,ngmi_total+1:ngmimax)) ! SJM
     ENDIF
  endif
  ngmi_total=ngmi3+ngmi_total

  gMIData%n1b11=ngmi_total
  gMIData%n1b11=gMIData%n1b11-4
  ngmi_total=ngmi_total-4
  !stop
  PRINT *,'Number 1CGMI reads     : ',ngmi1, ngmi2, ngmi3, '   ', ngmi_total
  !print*, gmiData%secondOfDay(ngmi1+1:(ngmi_total-ngmi3))
  PRINT *,'File read status 1CGMI : ',st_1, st_2, st_3

!  SFM  end    12/06/2013  reworked checks for nil and dead files

  ifdpr(1:1)='N'

!  SFM  start  12/06/2013
!...Compute orbit & date data
  month_get = (date_number/100) - (date_number/10000) * 100
  year = date_number/10000
  dd   = date_number - year*10000 - month_get * 100
  print*, month_get
  call read_emis_s0_map(month_get)
  write(lut_file,'(A22,I2.2,A4)') 'gpm_surfmaps/surfmaps.',month_get,'.dat'
  call read_surfmap_month(lut_file)
  call readwfract()
  
!-----------------------
  nBSize=300
  call allocateDPRSpace(dPRData, ngates, ndPRrays, nPRScanM)

  call param_set_BMCV(5) 
 iStart=1

  iEnd=1
  dx=0

  tbRgrid=-99             ! set to default

  nchunks=int(ndpr/300)
  if(nchunks*300==ndpr) nchunks=nchunks-1
  print*, 'ialg=', ialg, ndpr

  call sst2(gmiData,geoData,gMIData%mm,gMIData%n1b11)
  !stop
  do i=0,3!nchunks
    
  enddo
  f2AKuc(:)=f2AKu(:)
  f2ADPRc(:)=f2ADPR(:)
  f2AkuENVc(:)=f2AkuENV(:)
  i2ADPRc=i2ADPR
  i2akuc=i2aku
  i2AKuENVc=i2AKuENV
  jobnamec(:)=jobname
  ndpr1=ndpr
  !return
  !call deallocatecGMISpace(gmiData)

!  SFM  start  10/22/2013

!...Isolate seaice file name w/o prefixes
    iSEAICE1 = 1
    DO i=iSEAICE,1,-1
      IF (fSEAICE(i) .EQ. '/') THEN
        iSEAICE1 = i+1
        EXIT
      ENDIF
    ENDDO

!...Isolate snow file names w/o prefixes
    iSNOW1 = 1
!=====================================
!  SFM  This segment of code gets commented out if using autosnow
!       split files. If using the ims data or autosnow global data,
!       leave it in.
    DO i=iSNOW,1,-1
      IF (fSNOW(i) .EQ. '/') THEN
        iSNOW1 = i+1
        EXIT
      ENDIF
    ENDDO
!=====================================

!  SFM  end  10/22/2013
!  SFM  start  01/02/2014

!..Transfer remaining meta data to output file
   WRITE(UNIT=*,FMT=500) st_2aku, st_2adpr, st_1, st_2, st_3, st_2akuenv,      &
                         st_snow, st_seaice
500 FORMAT('File read status ALLof : ',8I5)
   call meta_for_outputfile(fSNOW(iSNOW1:iSNOW), fSEAICE(iSEAICE1:iSEAICE),    &
                            orbitNumber, rseed1, rseed2, algorithmVersion,     &
			    ABS(st_1), ABS(st_2), ABS(st_3), ABS(st_2akuenv),  &
			    ABS(st_2aku), ABS(st_2adpr), ABS(st_snow),         &
			    ABS(st_seaice))
!  SFM  end  01/02/2014

   
   
!stop
 
end subroutine mainfort


subroutine mainfortpy()

  use iso_c_binding    !iso c binding statement
  use globalData
  use emissmod
  use nbinmod
  use missingMod
  use Tables_frac


  implicit none
  integer :: rseed1, rseed2, ifs
  real :: tbRgrid(14,49,9300)  ! resampled brightness temperatures on same grid 
                              !   as 2adpr data; (channel/ray/scan)


  character(c_char) :: jobname(255), f1ctmi1(1000), f1ctmi2(1000), &
       f1ctmi3(1000), f2AKu(1000), f2aDPR(1000), f2AKuENV(1000), f2CMB(1000)
  character(c_char) :: fSNOW(1000), fSEAICE(1000)
  character(c_char) :: outcmb(1000) !iso c binding

  integer :: igmi1, igmi2, igmi3, i2AKu, i2ADPR, i2AKuENV, i2CMB, ialg
  integer :: iSNOW, iSEAICE
  integer :: i1ctmi, i1c21, i2a23
  real    :: mu
  integer :: nGMIS1, nGMIS2, nS1f, nS2f, nchunks
  character*3 :: ifdpr
  INTEGER*4 :: orbitNumber    ! retrieved orbit number
  integer :: ngmi1,ngmi2,ngmi3, j
  integer :: isnow1, iseaice1, ndpr1
  real    :: lastGood          !  SFM  04/16/2014   added for M.Grecu
  
  character*90 :: lut_file !SJM 7/9/2014

  character(16), parameter :: algorithmVersion = "2BCMB_20180919"
  integer    :: date(3)=(/0,0,0/), istat,    kk1, kk2
  integer(1) :: autosnow(9000,4500)                 !4km global grid
  
  INTEGER :: nscan_2akuenv, nray2akuenv, nbin2akuenv, nwater, nwind

  INTEGER*4     date_number      ! file date in format yyyymmdd
  INTEGER*4     month_get        ! computed file month
  INTEGER*4     st_1, st_2, st_3 ! status codes for 1CGMI files
  
  INTEGER*4     readgmi, readtmi         ! declare function
  INTEGER*4     readenv, readenvx          ! declare function
  INTEGER*4     readdprpflag, readdprpflagx     ! declare function
  INTEGER*4     readdprtpflag     ! declare function
  
  integer*4     i, ifract, idir
  real :: x1L, x2L

  rec=0
  
  mu=cos(52.8/180*3.14159)
  do i=1,7
     !call mie2(i-3.)
  enddo
  
  nmfreq=8
  nmu=5
  call init_random_seed(rseed1,rseed2)
  print *, 'Random Seeds : ',rseed1,rseed2
  print*, ialg

  call readtablesLiang2(nmu,nmfreq)         !This option reads in Liang's mu=2 table
  
  call makeHashTables()
  call cloud_init(nmfreq) 
  call initWFlag(nmfreq)
  call init_nbin
  call init_missing_flags


  nmF=4
  if(.not.(allocated(ip))) &      
       allocate(ip(nMF), iGMIChan(nMF), iSimF(nMF), iGMIf(nMF))
  
  nPrScanM = 300; ngates=nbin; ndPRrays=49;
  nGMIMax  = 3700; nGMIS1=9; nGMIS2=4; nGMIrays=221;
  gridData%dx=0.05

  nMemb=50
  ifdpr(1:1)='N'
  
  call read_geodat(geoData%lsflag, geoData%sstdata)

  nS1f=9
  nS2f=4
  ngmi_total = 0
  
  call allocatecGMISpace(gMIData, nS1f, nS2f, nGMIrays, nGMIMax, nmfreq, nmemb)

  call readwfract()
  
!-----------------------
  nBSize=300
  call allocateDPRSpace(dPRData, ngates, ndPRrays, nPRScanM)

  call param_set_BMCV(5) 
  iStart=1

  iEnd=1
  dx=0
  
  tbRgrid=-99             ! set to default
  
  nchunks=int(ndpr/300)
  if(nchunks*300==ndpr) nchunks=nchunks-1
end subroutine mainfortpy

subroutine do_chunk(i,ialg, idir)
  use globalData
  implicit none

  integer :: i, ialg
  integer :: ifract, idir, j
  real :: x1L, x2L
  INTEGER*4     readenv, readenvx 
  real :: dprrain(49,300)
  print*,'ichunk=', i, nMemb, istart
  ic=i*nBSize
  !jobnamec(1)='j'
  !jobnamec(2)='u'
  !jobnamec(3)='n'
  !jobnamec(4)='k'
  
  call printjobname(jobnamec)
  print*, f2AKuc(1:i2akuc)
  print*, f2ADPRc(1:i2adprc)
  print*, ic
  print*, 'do_chunk'
  
  !return
  if(ialg.eq.2) then
     call read2akut(jobnamec, f2AKuc(1:i2akuc),                             &
          dPRData%n1c21,dPRData%zku1c21, dPRData%zka1c21,                &
          dPRData%snrRatioku, dPRData%snrRatioka,                        &
          dPRData%srtPIAku,dPRData%dsrtPIAku,dPRData%dsrtPIAka,          &
          dPRData%srtsigmaPIAku, dPRData%dsrtsigmaPIAku,                 &
          dPRData%dsrtsigmaPIAka,                                        &
          dPRData%sigmaZeroKu, dPRData%sigmaZeroKa,                      & !SJM 12/3/14
          dPRdata%sclon, DPRdata%sclat,                                  & !SJM 3/31/16
          dPRData%xlon,dPRData%xlat,dPRData%badRayFlag,                  &
          dPRData%rainFlagBad,dPRData%node,dPRData%rainType,             &
          dPRData%scAngle, ic,                                           &
          nBSize,dPRData%freezH, dPRData%surfaceZKu,                     &
          dPRData%iLandOcean,dPRData%srtrelPIAku,dPRData%dsrtrelPIA,     &
          dPRData%piaHB,                                                 &
          dPRData%ioqualityflagku, dPRData%ioqualityflagdpr,             &
          f2ADPRc(1:i2adprc),dprrain,dPRData%BBbin,dPRData%binRealSurface, &
          dPRData%localZenithAngle, dPRData%elevation, st_2adpr,         &
          dPRData%secondOfDay,dPRData%NSRelibFlag,dPRData%MSRelibFlag,   &
          dPRdata%snowIceCover, dPRdata%seaIceConcentration, dPRdata%cBEst)
  else
     call read2aku(jobnamec, f2AKuc(1:i2akuc),                             &
          dPRData%n1c21,dPRData%zku1c21, dPRData%zka1c21,                &
          dPRData%snrRatioku, dPRData%snrRatioka,                        &
          dPRData%srtPIAku,dPRData%dsrtPIAku,dPRData%dsrtPIAka,          &
          dPRData%srtsigmaPIAku, dPRData%dsrtsigmaPIAku,                 &
          dPRData%dsrtsigmaPIAka,                                        &
          dPRData%sigmaZeroKu, dPRData%sigmaZeroKa,                      & !SJM 12/3/14
          dPRdata%sclon, DPRdata%sclat,                                  & !SJM 3/31/16
          dPRData%xlon,dPRData%xlat,dPRData%badRayFlag,                  &
          dPRData%rainFlagBad,dPRData%node,dPRData%rainType,             &
          dPRData%scAngle, ic,                                           &
          nBSize,dPRData%freezH, dPRData%surfaceZKu,                     &
          dPRData%iLandOcean,dPRData%srtrelPIAku,dPRData%dsrtrelPIA,     &
          dPRData%piaHB,                                                 &
          dPRData%ioqualityflagku, dPRData%ioqualityflagdpr,             &
          f2ADPRc(1:i2adprc),dprrain,dPRData%BBbin,dPRData%binRealSurface, &
          dPRData%localZenithAngle, dPRData%elevation, st_2adpr,         &
          dPRData%secondOfDay,dPRData%NSRelibFlag,dPRData%MSRelibFlag,   &
          dPRdata%snowIceCover, dPRdata%seaIceConcentration, dPRdata%cBEst)
  endif
  print*, maxval(dPRData%xlat)
  IF (i .EQ. 1) PRINT *,'File read status 2ADPR : ',st_2adpr
  print*, dPRData%n1c21, st_2akuenv
  
  IF (st_2akuenv .GE. 0) THEN
     st_2akuenv = readenv(jobnamec, ialg, f2AkuENVc(1:i2AkuENVc), &
          dPRData%n1c21,    &
          ic, nBSize, dPRData%envQv, dPRData%envTemp,                  &
          dPRData%envPress, dPRData%envSfcWind,dPRData%envSfcWindU, &
          dPRData%envSfcWindV, dPRData%envSknTemp,    &
          dPRData%envSfcTemp, dPRData%envSfcPress,dPRData%envCloud)
  ENDIF
  
  IF (i .EQ. 1) PRINT *,'File read status 2ADPR : ',st_2adpr
  IF (i .EQ. 1) PRINT *,'File read status ENV   : ',st_2akuenv
  
  !begin SJM 12/12/2014 fix missing GMI data when orbit spans 2 days (unwrap times)
  if(i .eq. 0) then
     dprstart_sec = dPRData%secondOfDay(1)
     gmistart_sec = gmiData%secondOfDay(1)
  endif
  !print*, dprstart_sec, gmistart_sec
  do j=1,ngmi_total
     if(gmiData%secondOfDay(j) .lt. gmistart_sec) gmiData%secondOfDay(j)=gmiData%secondOfDay(j)+86400
  end do
  
  do j=1,dPRData%n1c21
     if(dPRData%secondOfDay(j) .lt. dprstart_sec) dPRData%secondOfDay(j)=dPRData%secondOfDay(j)+86400
  end do
  !end SJM 12/12/2014
  
  !  SFM  start  04/16/2014;  for M.Grecu, node processing change issues
  iStart=1
  do while(gmiData%secondOfDay(istart)< &
       dPRData%secondOfDay(1).and. &
       gmiData%secondOfDay(istart+1)<&
       dPRData%secondOfDay(1).and. &
       istart<ngmi_total) 
     istart=istart+1
  end do
  iEnd=iStart
  !begin  WSO 4/24/14 replaced istart with iEnd as limiter of do while
  do while(gmiData%secondOfDay(iEnd)< &
       dPRData%secondOfDay(dPRData%n1c21).and. &
       gmiData%secondOfDay(iEnd+1)<&
       dPRData%secondOfDay(dPRData%n1c21).and. &
       iEnd<ngmi_total) 
     iEnd=iEnd+1
  end do
  if(iEnd<ngmi_total) iEnd=iEnd+1
  iEnd=min(ngmi_total,iEnd+70)
  iStart=max(1,iStart-70)
  print*, 'istart=',iStart,iEnd, i
  call reSampleGMI(gmiData,iStart,iEnd,gmi2Grid,i)
  
  sysdN=0.0
  sysdN=-.25 !new calibration !-0.25 ITE
  print*, 'allocate dPRRet'
  call allocateDPRRetSpace(dPRRet, nMemb, nmfreq, ngates, &
       ndPRrays, dPRData%n1c21)
  ifract=0
  
  x1L=gmiData%sclon(istart)
  x2L=gmiData%s1lon(110,istart)
  if(x1L>180) x1L=x1L-180*2.
  if(x2L>180) x2L=x2L-180*2.
  if(x1L>x2L .and. x1L-x2L<15) then
     idir=1
  else
     if(x2L-x1L<15) then
        idir=-1
     else
        idir=1
     endif
  endif
  print*, 'lon_istart=', x1L,x2L, 'idir=',idir
  
  !call radarRetSubEmpty(geoData,dPRData, &
  !     dPRRet,gridEnvData, gmi2Grid,        &
  !     gmiData, nmu, nMemb, nmfreq, ngates, sysdn, ic, tbRgrid, dprrain, &
  !     i, st_2adpr,orbitNumber,iStart,iEnd)
  
  
  
  ! iStart=iEnd       SFM  04/16/2014   for M.Grecu; node processing
  print*, 'ichunks=',i, ran1()
end subroutine do_chunk

subroutine closefiles(ialg)
  implicit none
  integer :: ialg
  if(ialg==2) then
      call closeoutputfilet()
   else
      call closeoutputfile()
   end if
 end subroutine closefiles

subroutine getzs(zku,zka,nx,ny,nz)
  use globalData
  real :: zku(nz,ny,nx),zka(nz,ny,nx)
  integer :: nx,ny,nz
  print*, nx, ny, nz
  print*, shape(dPRData%zku1c21)
  zku=dPRData%zku1c21
  zka=dPRData%zka1c21
end subroutine getzs


subroutine do_chunkx(i,ialg, idir)
  use globalData
  implicit none

  integer :: i, ialg
  integer :: ifract, idir, j
  real :: x1L, x2L
  INTEGER*4     readenv, readenvx 
  real :: dprrain(49,300)
  print*,'ichunk=', i, nMemb, istart
  ic=i*nBSize
  !jobnamec(1)='j'
  !jobnamec(2)='u'
  !jobnamec(3)='n'
  !jobnamec(4)='k'
  !print*, jobnamec(:10)
  !print*, f2AKuc(1:i2akuc)
  !print*, f2ADPRc(1:i2adprc)
  !print*, ic
  !stop
  !return
  if(ialg.eq.2) then
     call read2akut(jobnamec, f2AKuc(1:i2aDPRc),                             &
          dPRData%n1c21,dPRData%zku1c21, dPRData%zka1c21,                &
          dPRData%snrRatioku, dPRData%snrRatioka,                        &
          dPRData%srtPIAku,dPRData%dsrtPIAku,dPRData%dsrtPIAka,          &
          dPRData%srtsigmaPIAku, dPRData%dsrtsigmaPIAku,                 &
          dPRData%dsrtsigmaPIAka,                                        &
          dPRData%sigmaZeroKu, dPRData%sigmaZeroKa,                      & !SJM 12/3/14
          dPRdata%sclon, DPRdata%sclat,                                  & !SJM 3/31/16
          dPRData%xlon,dPRData%xlat,dPRData%badRayFlag,                  &
          dPRData%rainFlagBad,dPRData%node,dPRData%rainType,             &
          dPRData%scAngle, ic,                                           &
          nBSize,dPRData%freezH, dPRData%surfaceZKu,                     &
          dPRData%iLandOcean,dPRData%srtrelPIAku,dPRData%dsrtrelPIA,     &
          dPRData%piaHB,                                                 &
          dPRData%ioqualityflagku, dPRData%ioqualityflagdpr,             &
          f2ADPRc(1:i2adprc),dprrain,dPRData%BBbin,dPRData%binRealSurface, &
          dPRData%localZenithAngle, dPRData%elevation, st_2adpr,         &
          dPRData%secondOfDay,dPRData%NSRelibFlag,dPRData%MSRelibFlag,   &
          dPRdata%snowIceCover, dPRdata%seaIceConcentration, dPRdata%cBEst)
  else
     call read2akux(jobnamec, f2AKuc(1:i2aDPRc),                             &
          dPRData%n1c21,dPRData%zku1c21, dPRData%zka1c21,                &
          dPRData%snrRatioku, dPRData%snrRatioka,                        &
          dPRData%srtPIAku,dPRData%dsrtPIAku,dPRData%dsrtPIAka,          &
          dPRData%srtsigmaPIAku, dPRData%dsrtsigmaPIAku,                 &
          dPRData%dsrtsigmaPIAka,                                        &
          dPRData%sigmaZeroKu, dPRData%sigmaZeroKa,                      & !SJM 12/3/14
          dPRdata%sclon, DPRdata%sclat,                                  & !SJM 3/31/16
          dPRData%xlon,dPRData%xlat,dPRData%badRayFlag,                  &
          dPRData%rainFlagBad,dPRData%node,dPRData%rainType,             &
          dPRData%scAngle, ic,                                           &
          nBSize,dPRData%freezH, dPRData%surfaceZKu,                     &
          dPRData%iLandOcean,dPRData%srtrelPIAku,dPRData%dsrtrelPIA,     &
          dPRData%piaHB,                                                 &
          dPRData%ioqualityflagku, dPRData%ioqualityflagdpr,             &
          f2ADPRc(1:i2adprc),dprrain,dPRData%BBbin,dPRData%binRealSurface, &
          dPRData%localZenithAngle, dPRData%elevation, st_2adpr,         &
          dPRData%secondOfDay,dPRData%NSRelibFlag,dPRData%MSRelibFlag,   &
          dPRdata%snowIceCover, dPRdata%seaIceConcentration, dPRdata%cBEst, &
          dPRData%envTemp)
  endif
  print*, maxval(dPRData%xlat)
  IF (i .EQ. 1) PRINT *,'File read status 2ADPR : ',st_2adpr
  print*, dPRData%n1c21, st_2akuenv
  
  IF (st_2akuenv .GE. 0) THEN
     st_2akuenv = readenvx(jobnamec, ialg, f2AkuENVc(1:i2AkuENVc), &
          dPRData%n1c21,    &
          ic, nBSize, dPRData%envQv, dPRData%envTemp,                  &
          dPRData%envPress, dPRData%envSfcWind,dPRData%envSfcWindU, &
          dPRData%envSfcWindV, dPRData%envSknTemp,    &
          dPRData%envSfcTemp, dPRData%envSfcPress,dPRData%envCloud)
  ENDIF
  
  IF (i .EQ. 1) PRINT *,'File read status 2ADPR : ',st_2adpr
  IF (i .EQ. 1) PRINT *,'File read status ENV   : ',st_2akuenv
  
  !begin SJM 12/12/2014 fix missing GMI data when orbit spans 2 days (unwrap times)
  if(i .eq. 0) then
     dprstart_sec = dPRData%secondOfDay(1)
     gmistart_sec = gmiData%secondOfDay(1)
  endif
  !print*, dprstart_sec, gmistart_sec
  do j=1,ngmi_total
     if(gmiData%secondOfDay(j) .lt. gmistart_sec) gmiData%secondOfDay(j)=gmiData%secondOfDay(j)+86400
  end do
  
  do j=1,dPRData%n1c21
     if(dPRData%secondOfDay(j) .lt. dprstart_sec) dPRData%secondOfDay(j)=dPRData%secondOfDay(j)+86400
  end do
  !end SJM 12/12/2014
  
  !  SFM  start  04/16/2014;  for M.Grecu, node processing change issues
  iStart=1
  do while(gmiData%secondOfDay(istart)< &
       dPRData%secondOfDay(1).and. &
       gmiData%secondOfDay(istart+1)<&
       dPRData%secondOfDay(1).and. &
       istart<ngmi_total) 
     istart=istart+1
  end do
  iEnd=iStart
  !begin  WSO 4/24/14 replaced istart with iEnd as limiter of do while
  do while(gmiData%secondOfDay(iEnd)< &
       dPRData%secondOfDay(dPRData%n1c21).and. &
       gmiData%secondOfDay(iEnd+1)<&
       dPRData%secondOfDay(dPRData%n1c21).and. &
       iEnd<ngmi_total) 
     iEnd=iEnd+1
  end do
  if(iEnd<ngmi_total) iEnd=iEnd+1
  iEnd=min(ngmi_total,iEnd+70)
  iStart=max(1,iStart-70)
  print*,  'istart=',iStart,iEnd, i
  call reSampleGMI(gmiData,iStart,iEnd,gmi2Grid,i)
  
  sysdN=0.0
  sysdN=-.25 !new calibration !-0.25 ITE
  print*, 'allocate dPRRet'
  call allocateDPRRetSpace(dPRRet, nMemb, nmfreq, ngates, &
       ndPRrays, dPRData%n1c21)
  ifract=0
  
  x1L=gmiData%sclon(istart)
  x2L=gmiData%s1lon(110,istart)
  if(x1L>180) x1L=x1L-180*2.
  if(x2L>180) x2L=x2L-180*2.
  if(x1L>x2L .and. x1L-x2L<15) then
     idir=1
  else
     if(x2L-x1L<15) then
        idir=-1
     else
        idir=1
     endif
  endif
  print*, 'lon_istart=', x1L,x2L, 'idir=',idir
  
  !call radarRetSubEmpty(geoData,dPRData, &
  !     dPRRet,gridEnvData, gmi2Grid,        &
  !     gmiData, nmu, nMemb, nmfreq, ngates, sysdn, ic, tbRgrid, dprrain, &
  !     i, st_2adpr,orbitNumber,iStart,iEnd)
  
  
  
  ! iStart=iEnd       SFM  04/16/2014   for M.Grecu; node processing
  print*, 'ichunks=',i, ran1()
end subroutine do_chunkx
