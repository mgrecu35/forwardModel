module globalData
  use iso_c_binding    !iso c binding statement
  use f90DataTypes
  use ran_mod
  use gif_util
  use weight
  use nbinMod
  implicit none

  integer :: nmfreq,nPrScanM, nGMIMax, ngates, nMemb, nmu, ndPRrays
  integer :: nGMIlowf, nGMIhf, nGMIrays

  type (dPRDataType)    :: dPRData
  type (dPRRetType)     :: dPRRet 
  type (cgMIDataType)    :: gMIData
  type (gridDataType)   :: gridData
  type (gridEnvDataType):: gridEnvData
  type (gmi2GridType)   :: gmi2Grid
  integer, allocatable  :: ip(:), iGMIChan(:), iSimF(:), iGMIf(:)
  integer               :: nmF, k
  real                  :: fobj
  type (geoDataType) :: geoData
  integer:: i0,j0!, ilandSea, igetlandsea
  real:: m1,m2, x, y, sstv, w1,w2, wt
  integer :: row, col, ic, nBSize
  integer*4 :: iStart, iEnd
  real ::  mval, x1(49), x2(49), dx(49)
  character*50 :: fname 
  real :: sysdn
  integer :: ncid, rec, iret, create_netcdf_file
  integer :: ndpr, year, jday, dd
!!  SFM  04/06/2013  Changed file name lengths to 1000
  real    :: rlon(49,10000), rlat(49,10000)
  integer :: iLandOcean(49,10000), rPrecipFlag(49,10000), imsFlag(49,10000), &
             reynoldsICE(49,10000)
  real    :: rSfcPrecip(49,10000), reynoldsSST(49,10000)
  real    :: emissOut(221,3300)
  character(c_char) :: jobnamec(255),&
      f2AKuc(1000), f2aDPRc(1000), f2AKuENVc(1000), f2CMBc(1000)
  INTEGER*4     st_2akuenv, st_2aku, st_2adpr  ! additional status codes
  INTEGER*4     st_snow, st_seaice, st_cmb     ! additional status codes
  INTEGER*4     ngmi_total       ! Accumulated 1CGMI read total
  real :: gmistart_sec, dprstart_sec !SJM 12/12/2014
  integer :: i2ADPRc, i2AKuc,i2AKuENVc
end module globalData
