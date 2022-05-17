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
  
  !call read_geodat(geoData%lsflag, geoData%sstdata)

  nS1f=9
  nS2f=4
  ngmi_total = 0
  
  call allocatecGMISpace(gMIData, nS1f, nS2f, nGMIrays, nGMIMax, nmfreq, nmemb)

  !call readwfract()
  
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

subroutine do_chunkpy(i,ialg, idir)
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
  

end subroutine do_chunkpy



subroutine getzs(zku,zka,nx,ny,nz)
  use globalData
  real :: zku(nz,ny,nx),zka(nz,ny,nx)
  integer :: nx,ny,nz
  print*, nx, ny, nz
  print*, shape(dPRData%zku1c21)
  zku=dPRData%zku1c21
  zka=dPRData%zka1c21
end subroutine getzs


