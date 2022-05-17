!  SFM  04/06/2013  Code changes from M.Grecui
!
subroutine setemtbm(emtbm,emtb,nmfreqm,nmemb)
real :: emtb(2,nmfreqm,nmemb)
real :: emtbm(9)
integer :: ifreqG(9), ipolG(9)

ifreqG(1:9)=(/1,1,2,2,3,4,4,5,5/)
ipolG(1:9)=(/1,2,1,2,1,1,2,1,2/)
do i=1,9
   emtbm(i)=sum(emtb(ipolG(i),ifreqG(i),1:nmemb))/nmemb
enddo
end subroutine setemtbm
real function corrcoef(x,y,n)
  integer :: n
  real :: x(n), y(n)
  sumx=0
  sumx2=0
  sumxy=0
  sumy=0
  sumy2=0
  do i=1,n     
     sumx  = sumx + x(i) 
     sumx2 = sumx2 + x(i) * x(i)
     sumxy = sumxy + x(i) * y(i)
     sumy  = sumy + y(i)                                
     sumy2 = sumy2 + y(i) * y(i)                        
  end do
  r = (sumxy - sumx * sumy / n) /    &           
       sqrt((sumx2 - sumx**2/n) * (sumy2 - sumy**2/n))
  corrcoef=r
end function corrcoef

real function covar(x,y,n)
  integer :: n
  real :: x(n), y(n)
  sumx=0
  sumx2=0
  sumxy=0
  sumy=0
  sumy2=0
  do i=1,n     
     sumx  = sumx + x(i) 
     sumy  = sumy + y(i)                                
  end do
  sumx=sumx/n
  sumy=sumy/n
  do i=1,n     
     sumxy  = sumxy+ (x(i)-sumx)*(y(i)-sumy)                                
  end do
  covar=sumxy/(n-1)
  if(covar.gt.300) then
     !do i=1,n
     !   write(*,*) x(i),y(i)
     !enddo
     !stop
  endif
end function covar


real function stddev(x,n)
  integer :: n
  real :: x(n)
  xm=sum(x(1:n))/n
  xs=sqrt(sum((x(1:n)-xm)**2)/(n-1))
  stddev=xs
end function stddev

real function srtpiaf(x,n)
  integer :: n
  real :: x(n)
  xm=0
  xm2=0
  do i=1,n
     xm=xm+10.**(-0.1*x(i))
     xm2=xm2+x(i)
  enddo
  xm=-10.*log10(xm/n)
  xm2=xm2/n
  !print*, xm2, xm
  if(xm2>0.) then
     srtpiaf=xm/xm2
  else
     srtpiaf=1
  endif
end function srtpiaf

subroutine DPRLonUnfold(x1,dPRData)
  use f90DataTypes
  implicit none
  real                   :: x1(49), dx(49)
  type (dPRDataType)     :: dPRData
  integer                :: i, j, k
  real                   :: dxl
 
  do i=1,49
     dx(i)=x1(i)-dPRData%xlon(i,1)
     dPRData%xlon(i,1)=dPRData%xlon(i,1)+dx(i) 
     do j=2,dPRData%n1c21
        dPRData%xlon(i,j)=dPRData%xlon(i,j)+dx(i)
        if(dPRData%xlon(i,j)<dPRData%xlon(i,j-1)) then
           dPRData%xlon(i,j)=dPRData%xlon(i,j)+360         
           !longitude unfolding
           dx(i)=dx(i)+360
        endif
     enddo
  enddo
end subroutine DPRLonUnfold


!  SFM  begin  12/13/2013; add file status flag
subroutine radarRetSub(nmu2,  nmfreq2,   icL, tbRgrid,               &
      dprrain,ichunk,orbNumb,ialg,idir)
!  SFM  end    12/13/2013
  use globalData
  use f90DataTypes
  use f90Types
  use cldclass
  use ran_mod
  use geophysEns
  use nbinMod
  !use tables2
  use weight
  Use BMCVparameters
  use emissMod
!begin  MG 10/29/15 add gEnv module
  use gEnv
!end    MG 10/29/15
!begin  WSO 9/14/13 incorporate missing flags
  use missingMod
!end    WSO 9/14/13
!begin  WSO 6/5/18 add limits to output variables
  use outputminmax
!end    WSO 6/5/18
  use LUT_def !SJM 7/9/2015
  implicit none
!  type (geoDataType) :: geoData
!  type (gridEnvDataType) :: gridEnvData
!  type (dPRDataType)     :: dPRData
!  type (dPRRetType)      :: dPRRet
!  type (gmi2GridType)    :: gmi2Grid
!  type (cgMIDataType)     :: gmiData
  integer :: nmu2, nmfreq2
  integer*4 :: ichunk
!  integer :: st_2adpr              ! file open status for 2adpr file
  real :: tbRgrid(14,49,9300), dprrain(49,300)
 ! integer :: tbRgridIn(9,49,9300)

  real                   :: meansfcRain,stddevSfcRain, tbout(14)
  type (radarRetType)    :: radarRet
  type (radarDataType)   :: radarData
  type (stormStructType) :: stormStruct
  type (retParamType)    :: retParam
  real :: xin363(363)
  real :: sfcRain(49,300),sfcRainStd(49,300)
  real :: rRate3D(nbin,49,300),  rRate3Dstd(nbin,49,300)
  real :: pwc3D(nbin,49,300),  pwc3Dstd(nbin,49,300)
  real :: zcKu3D(nbin,49,300), d03D(nbin,49,300), piaOut(49,300)

  real :: sfcRainMS(49,300),sfcRainStdMS(49,300),pia35m(49,300)
  real :: rRate3DMS(nbin,49,300),  rRate3DstdMS(nbin,49,300)
  real :: pwc3DMS(nbin,49,300),  pwc3DstdMS(nbin,49,300)
  real :: zcKu3DMS(nbin,49,300), zcKa3DMS(nbin,49,300), d03DMS(nbin,49,300), &
          piaOutKuMS(49,300), piaOutKaMS(49,300)

  integer :: ii, jj, iGMI, jGMI
  integer :: di(8), dj(8)  
  integer :: i, j, ig, jg, ntpw, nmemb1, itop, irand
  real    :: pia13m, rms1, rms2,  unSortedRR(200), corrcoef, sfcRain2, tpw_ij
  integer :: iy(200), kflag, it
  real    :: sysdNl, pia13mean
  integer :: iLandSea(5,5), i1, j1, igetlandsea, ic2, nobs, iit
  real :: a0(2,8), emiss(2), tb19v, tb19h, tb37v, tb37h, tb22v, tb85v, tb85h
  real :: meanTb85
  real :: stdTb85,kgain,kgain37, ymean(3)
!...Ensemble parameters
  real, allocatable ::  Yens(:,:), Xens(:,:), Yobs(:), Xup(:)
  integer :: ibatch
  real    :: stddev, srtpiaf
  real    :: FWHMx, FWHMy, tbconv(2), tbconvEns(2,100)
  integer :: dnx,dny, ik
  integer :: ipol(15), ifreq(15), iobs(15), ifreq1
  real, allocatable :: ndn(:), ndnp(:,:), xscalev(:), logdNwf(:), randemiss(:), dwind(:)
  real, allocatable  :: rhPCij(:,:), cldwPCij(:,:)
  real :: cldw(nlayer), rh(nlayer), pia13s
  integer :: nx,ny, icount, imin
  real ::  xm1,xs,rmsmin, prob, probtot, rmstot
  real :: piaR(100), fPIA, z13m
  integer :: ntbpix, ntbpix2
  real :: emtbm(9)
  real :: zminsc
  real :: realOut(49)
  real :: w10(49,300), w10_out_NS(49,300), w10_out_MS(49,300), w10_min, w10_max, emis, relAz
  real :: w10_rms_NS(49,300), emis_rms_NS(49,300,13), w10_rms_MS(49,300), emis_rms_MS(49,300,13)
  real :: dZms(49,300) !! MS addition Feb 10, 2017
  integer :: msFlag(49, 300) !!WSO addition Feb 11, 2017
!begin  WSO 2/8/17 new variables
  integer :: multiscatcalc_NS(49, 300), multiscatcalc_MS(49, 300)
  integer :: algotype_NS(49, 300), algotype_MS(49, 300)
  integer :: profclass_NS(49, 300), profclass_MS(49, 300)
  real :: subfootvariability_NS(49, 300), subfootvariability_MS(49, 300)
  real :: multiscatsurface_NS(49, 300), multiscatsurface_MS(49, 300)
  real :: skintempsigma_NS(49, 300), skintempsigma_MS(49, 300)
  real :: columnvaporsigma_NS(49, 300), columnvaporsigma_MS(49, 300)
  real :: columncloudliqsigma_NS(49, 300), columncloudliqsigma_MS(49, 300)
  real :: errorofdatafit_NS(49, 300), errorofdatafit_MS(49, 300)
  real :: initnw_NS(nbin, 49, 300), initnw_MS(nbin, 49, 300) 
  real :: princomp_NS(5, 49, 300), princomp_MS(5, 49, 300)
  real :: surfprecipbiasratio_NS(49, 300), surfprecipbiasratio_MS(49, 300)
!end    WSO 2/8/17 
  integer :: l, ipias(2)
  character*3 :: ifdpr, iftest
  character*90 :: outfile
  integer :: ink
  integer :: ifreqG(15), ipolG(15) 
  DOUBLE PRECISION input(6)
  DOUBLE PRECISION output(2)
  real :: wfract(5,5), wfractm, wfractsd
  real                    emissv(n_chan)
  real                    emissh(n_chan)
  real                    emissv_std(n_chan)
  real                    emissh_std(n_chan)
  real :: emis_eofs(nmemb,12) !SJM 7/9/2015
  integer :: stype!SJM 7/9/2015
  real  :: vLand(18,18), vOcean(10,10)
  real  :: pMLand(18), pMOcean(10)
  real  :: mTbLand(9), mTbOcean(9)
  real  :: stTbLand(9), stTbOcean(9)
  double precision  :: xin(18), xpred(18), yout(9)
  real, allocatable :: emissoutL(:,:,:), emis_out_NS(:,:,:), emis_out_MS(:,:,:) !sjm 8/10/15
!begin  WSO 8/19/13 change Nw variable name (not dN) and add mu
  real :: cldwprof(88), cldiprof(88), log10NwMean(88), mu_mean_prof(88)
  integer *2 :: env_nodes(10, 49)
  real :: env_levs(10), ray_angle, pi
!end    WSO 8/19/13
  real :: lFract(49,300), sprobs, probs(100), rmsS(100)
  real :: covar, xf
  integer  :: orbNumb
  !begin SJM 7/25/14
  real :: s0Ku, s0Ka, s0stdKu, s0stdKa, s0corr, ds0Ku, ds0Ka
  real :: sigmaZeroVarKu(49,300), sigmaZeroVarKa(49,300), sigmaZeroCov(49,300)
  !end SJM 7/25/2014
!begin WSO 8/8/13
  real :: gatelength
  real :: depthBB, depthML, depth
  real :: mu_mean(49, 300)
  real :: mu_meanMS(49, 300)
  real :: scLatPR(49,300),scLonPR(49,300),wfmap(49,300), fpmap(49,300,15), fpmapN(49,300,15)
  real :: S1eiaPR(49,300), S2eiaPR(49,300)
  real :: mlwc_frac(10, 49, 300)
  real :: mrate_frac(10, 49, 300)
  real :: mlwc_fracMS(10, 49, 300)
  real :: mrate_fracMS(10, 49, 300)
  real :: sfcRainLiqFrac(49, 300)
  real :: sfcRainLiqFracMS(49, 300)
  real :: tbMax1(15), tbMin1(15)

!  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
!  SFM  begin  06/22/2014
  real :: wfractPix, windPert(100), windPertU(100), windPertV(100), qvPert(100), dnqv
!  SFM  end    06/22/2014
!  SFM  end    07/29/2014
!end   WSO 8/8/13
  real :: covTb(49,300,15,15), tbMax(49,300,15), tbMin(49,300,15), &
       tbMean(49,300,15)
  real :: invCovTb(49,300,15,15)
  real :: tbout2D(49,300,15), tb(49,300,15), tbNoOcean(49,300,15), &
       tbout2DNoOcean(49,300,15), tbObs(49,300,15)
  real :: dfdtb(49,300,15), rerr(15), tb0(49,300,15), fem(15) , &
       tb0MS(49,300,15), tbNoOceanMS(49,300,15), tbout2DNoOceanMS(49,300,15),&
       tbout2DMS(49,300,15)

  integer :: actOb(49,300), iactOb
  integer :: jk, nf
  integer :: dig               ! SFM  04/16/2014  for M.Grecu
  real   :: cl(9,25), xin25(25),dtb(9)
  real   :: ebar, minl
  real, allocatable :: geoloc(:), hFreqTbs(:,:), PRgeoloc(:), hFreqPRg(:,:,:)
  !integer*4 :: istart, iend
  integer :: iconv, ialg, icL
  real :: nubfc, stdpia35

  integer,parameter :: nscans=300, npixs=25, nlev=88, nchans=13
  integer :: nfreq, idir
  integer :: pType(nscans,npixs)
  real :: sfcTemp(nscans,npixs), cldw3d(nscans,npixs,nlev)
  integer :: clutFree(nscans,npixs)
  real :: pRate(nscans,npixs,nlev), swc3d(nscans,npixs,nlev), tbobsT(nscans,npixs,nchans)
  real :: z13(nscans,npixs,nlev),emiss2d(nscans,npixs,nchans)
  real :: nw3d(nscans,npixs,nlev), press3d(nscans,npixs,nlev), &
       airTemp3d(nscans,npixs,nlev),qv3d(nscans,npixs,nlev)
  integer :: binNodes(nscans,npixs,5)
  integer :: envNode(nscans,npixs,10)
  real    :: pRateOut(nscans,npixs,nlev), swcOut(nscans,npixs,nlev), nwOut(nscans,npixs,nlev)
  integer :: sfcBin(nscans,npixs)
  real    :: tbsim(nscans,npixs,nchans)
!begin  WSO 8/30/13 prescribe levels for environmental parameters
  data env_levs/18., 14., 10., 8., 6., 4., 2., 1., 0.5, 0./
  data pi/3.14159265/
  tbMax1(1:9)=(/ 285.96795654,284.71334839,286.23388672,&
       284.92977905,286.37451172,&
       287.90933228,286.43515015,284.8597412,284.51959229/)
  tbMin1(1:9)=(/175.59274292,95.61299896,&
       215.97111511,161.71897888,205.35340881,&
       143.96331787,143.96331787,87.26193237,87.26193237/)
  tbMin1(1:9)=(/167.43344116,   86.27742767, &
       200.01838684,  134.10232544, &
       233.02462769,&
       221.44563293,  159.45098877,  267.19430542,  242.64634705/)
  !print*, orbNumb,ichunk
  !call openascii(orbNumb,ichunk)
  !stop
  call readpcoeff(vLand,vOcean,pMLand,pMOcean,&
       mTbLand,mTbOcean,stTbLand,stTbOcean)
  dPRRet%cldwcoeff=0
  ifreqG(1:13)=(/1,1,2,2,3,4,4,5,5,6,6,7,8/)
  ipolG(1:13)=(/1,2,1,2,1,1,2,1,2,1,2,1,1/)
! 1 is V
! 2 is H
  !call readclust()
!begin  MG Sept 15 2015
  !call openenkffilensl(orbNumb,ichunk)
  allocate(hFreqPRg(49,dPRData%n1c21,4))
  hFreqPRg=missing_r4
!end  MG
  dZms=missing_r4 !! MS addition Feb 10, 2017
  msFlag = missing_i2 !! WSO addition Feb 11, 2017
  print*, iEnd, iStart, ichunk
  if(iEnd-iStart>3) then
     allocate(geoloc(2*(iEnd+1-iStart)*81))
     allocate(hFreqTbs((iEnd+1-iStart)*81,4))
     ic2=1

     do i=iStart,iEnd
        do j=70,150
           geoloc(ic2)=gmiData%S2lat(j,i)
           ic2=ic2+1
           geoloc(ic2)=gmiData%S2lon(j,i)
           ic2=ic2+1
        enddo
     enddo
     print*, minval(gmiData%S2lat(70:150,iStart:iEnd))
     print*, maxval(gmiData%S2lat(70:150,iStart:iEnd))
     print*, minval(gmiData%S2lon(70:150,iStart:iEnd))
     print*, maxval(gmiData%S2lon(70:150,iStart:iEnd))
     minl=minval(gmiData%S2lon(70:150,iStart:iEnd))
     do k=1,4
        ic2=1
        do i=iStart,iEnd
           do j=70,150
              hFreqTbs(ic2,k)=gmiData%gmiS2(k,j,i)
              ic2=ic2+1
           enddo
        enddo
     enddo
     
!begin  MG Sept 15 2015
!     allocate(hFreqPRg(49,dPRData%n1c21,4))
!     hFreqPRg=-99
!end    MG Sept 15 2015
     allocate(prgeoloc(dPRData%n1c21*49*2))
     ic2=1
     do i=1,dPRData%n1c21
        do j=1,49
           prgeoloc(ic2)=dPRData%xlat(j,i)
           ic2=ic2+1
           if(dPRData%xlon(j,i)+180>minl) then
              prgeoloc(ic2)=dPRData%xlon(j,i)
           else
              prgeoloc(ic2)=dPRData%xlon(j,i)!+360
           endif
           ic2=ic2+1
        enddo
     enddo
     print*, minval(dPRData%xlat)
     print*, maxval(dPRData%xlat)
     print*, minval(dPRData%xlon)
     print*, maxval(dPRData%xlon)

     !allocate(hFreqPRg(dPRData%n1c21,49,4))
     !print*, 
     !print*, geoloc(1:49*2)
     !print*, geoloc((iEnd+1-iStart-5)*81*2:(iEnd+1-iStart)*81*2)
     !print*, prgeoloc(1:49*2)
     do k=1,4
        call flannint(geoloc,  prgeoloc, hFreqTbs(:,k), hFreqPRg(:,:,k), &
             (iEnd+1-iStart)*81, 2, &
             dPRData%n1c21*49)
     enddo
     !print*, hFreqPRg(:,1,1)
     !print*, hFreqPRg(:,1,2)
     !print*, hFreqPRg(:,1,3)
     !print*, hFreqPRg(:,50,1)
     !print*, hFreqPRg(:,50,2)
     !print*, hFreqPRg(:,50,3)
     !stop
     deallocate(geoloc)
     deallocate(prgeoloc)
     deallocate(hFreqTbs)
  endif

  di=(/0,  0, -1, 1, -2, 2, 0, 0/)
  dj=(/-1, 1,  0, 0, 0, 0, 2, -2/)

  a0(1,1:8)=(/0.5098,4.4664E-3,-6.0427E-6,&
       -2.5285E-3,-2.3725E-3,9.8163E-4,-2.2269E-3,-1.3193E-3/)
  a0(2,1:8)=(/0.3186,-1.5225E-3,1.7213E-3,&
       -3.7164E-4,6.5607E-3,8.1213E-4,-1.7678E-3,-1.7250E-3/)
  
  nmemb1=nmemb
  nobs=2
  allocate(ndn(nmemb1),xscalev(nmemb1),&
       logdNwf(9*nmemb1), randemiss(nmfreq2*nmemb1*2)) 
  allocate(ndnp(10,nmemb1)) 
  allocate(emissoutL(49,dPRData%n1c21,13))

  ntbpix=0
  ntbpix2=0
!  begin  SFM  07/29/2014; for M.Grecu,  eliminate NANs
  do k=0,nmemb1-1
     windPert(k+1) = normal2(0.,.15) !SJM 2/4/15
     windPertU(k+1) = normal2(0.,.15) !SJM 2/4/15
     windPertV(k+1) = normal2(0.,.15)
     qvPert(k+1)=ran1()
  end do
!  end    SFM  07/29/2014
  iactOb=0
  do k=0,nmemb1-1
     ndn(k+1)=normal2(0.,1.0)
     xscalev(k+1)=ran1()
     do i=1,9
        ndnp(i,k+1)=normal2(0.,1.0)
     enddo
  enddo
  do k=1,2*nmfreq2*nmemb1
     randemiss(k)=.5*ran1()
  enddo
  !set random emissivity/sigm0 PCs
  do k=0,nmemb1-1
     do i=1,12
       emis_eofs(k+1,i) = normal2(0.,1.)
     end do
  end do
 

  rRate3D=0.
  rRate3Dstd=0.
  pwc3D=0.
  pwc3Dstd=0.
!begin  WSO 9/15/13 set to flag instead of 0
  d03D = missing_r4
  zcKu3D = missing_r4
  piaOut = missing_r4
!end    WSO 9/15/13 
  sfcRain=0.
  sfcRainStd=0.
  dPRRet%z35mod0=0.
!begin WSO 04/07/13
  rRate3DMS=0.
  rRate3DstdMS=0.
  pwc3DMS=0.
  pwc3DstdMS=0.
!begin  WSO 9/15/13 set to flag instead of 0
  d03DMS = missing_r4
  zcKu3DMS = missing_r4
  zcKa3DMS = missing_r4
  piaOutKuMS = missing_r4
  piaOutKaMS = missing_r4
!end    WSO 9/15/13 
  sfcRainMS=0.
  sfcRainStdMS=0.
!begin  WSO 9/15/13 set to flag instead of 0
  w10 = missing_r4
  w10_rms_NS = missing_r4
  emis_rms_NS = missing_r4
  w10_rms_MS = missing_r4
  emis_rms_MS = missing_r4
!end    WSO 9/15/13
!end WSO 04/07/13
!begin WSO 2/8/17 initialize new output variables
  multiscatcalc_NS = missing_i4
  algotype_NS = missing_i4
  subfootvariability_NS = missing_r4
  multiscatsurface_NS = missing_r4
  skintempsigma_NS = missing_r4
  columnvaporsigma_NS = missing_r4
  columncloudliqsigma_NS = missing_r4
  errorofdatafit_NS = missing_r4
  multiscatcalc_MS = missing_i4
  algotype_MS = missing_i4
  subfootvariability_MS = missing_r4
  multiscatsurface_MS = missing_r4
  skintempsigma_MS = missing_r4
  columnvaporsigma_MS = missing_r4
  columncloudliqsigma_MS = missing_r4
  errorofdatafit_MS = missing_r4
  profclass_NS = missing_i4
  initnw_NS = missing_r4
  princomp_NS = missing_r4
  surfprecipbiasratio_NS = missing_r4
  profclass_MS = missing_i4
  initnw_MS = missing_r4
  princomp_MS = missing_r4
  surfprecipbiasratio_MS = missing_r4
!end WSO 2/8/17
!begin WSO 9/15/13 set to flag instead of 0
!begin  WSO 8/8/13
  mlwc_frac = missing_r4
  mrate_frac = missing_r4
  mlwc_fracMS = missing_r4
  mrate_fracMS = missing_r4
  sfcRainLiqFrac = missing_r4
  sfcRainLiqFracMS = missing_r4
!end    WSO 8/8/13
!begin  WSO 8/19/13
  mu_mean = missing_r4
  mu_meanMS = missing_r4
!end    WSO 8/19/13
!end    WSO 9/15/13
  print*, nmfreq2
  call allocateDPRProfRet(radarRet,nmfreq2,nmemb1,ngates, 9)   ! allocates memory for the 1D 
  
  !...retrieval structures
  radarRet%rrate = 0.0  

  radarRet%tb=-99

  call allocateDPRProfData(radarData, ngates)                 ! allocates memory for the 
                                                              ! 1-D DPR observations
  call allocateStormStructData(stormStruct)                   ! allocates memory for the 5-node 
                                                              ! storm structure
  call setRetParam(retParam)
  call setrandclass(radarRet, nmu2)
  
  dPRRet%sfc_wind(1:nmemb)=radarRet%sfc_wind(1:nmemb)
  
  dPRRet%sfcRainEns=0
  stormStruct%iSurf=ngates
  radarData%ngates=ngates
  radarData%dr=0.25
  dPRRet%tb=-99
  dPRRet%emtb=-99
  dPRRet%emis=-99
  dPRRet%n9=0
  call allocGeophys(6,61,9,nmemb1,nmfreq2*nmemb1*2)
  call setdNwIcJcL(sysdN,nmemb1)
  nx=30+nbin*7
!begin  MG 9/18/13 changed 110 to 130
  ny=130
!end    MG 9/18/
! SFM  begin  03/27/2014; execution protections
  IF (ALLOCATED(xens)) deallocate(xens)
  IF (ALLOCATED(yens)) deallocate(yens)
  IF (ALLOCATED(rhPCij)) deallocate(rhPCij)
  IF (ALLOCATED(cldwPCij)) deallocate(cldwPCij)
! SFM  end    03/27/2014

  allocate(Xens(nx,nmemb),Yens(ny,nmemb),Yobs(ny), Xup(nx))
  allocate(rhPCij(nmemb1,nRhEofs), cldwPCij(nmemb1,nCldwEofs))
!00000000000000000000000000000000000000000000000000000000000000000000000000000
  lFract=-99.

  tb19v=0
  tb19h=0
  tb22v=0
  tb37v=0
  tb37h=0
  tb85v=0
  tb85h=0
  ic2=0
!begin  WSO 9/16/13
  emissoutL=missing_r4
!end    WSO 9/16/13
!begin WSO 04/07/13
   ifdpr(1:1) = 'N'
   iftest(1:1) = 'N'


  sfcRain=0
  !begin SJM 12/9/2014
  w10=dPRData%envSfcWind
  w10_out_NS=dPRData%envSfcWind
  w10_out_MS=dPRData%envSfcWind
  
  sigmaZeroVarKu = 0.
  sigmaZeroVarKa = 0.
  sigmaZeroCov = 0.
  S1eiaPR = 52.7
  S2eiaPR = 49.1
!   write(outfile,'(A64,I6.6,A4)') '/PANFS/user/home/smunchak/data/combAlg-tests/ensdata/InitialEns.',orbNumb,'.bin'
!   if(ichunk .eq. 0) then
!     open(31,file=outfile,form='unformatted') 
!   else
!     open(31,file=outfile,access='append',form='unformatted') 
!   endif
!   write(outfile,'(A64,I6.6,A4)') '/PANFS/user/home/smunchak/data/combAlg-tests/ensdata/NSFiltered.',orbNumb,'.bin'
!   if(ichunk .eq. 0) then
!     open(32,file=outfile,form='unformatted') 
!   else
!     open(32,file=outfile,access='append',form='unformatted') 
!   endif
!   write(outfile,'(A64,I6.6,A4)') '/PANFS/user/home/smunchak/data/combAlg-tests/ensdata/MSFiltered.',orbNumb,'.bin'
!   if(ichunk .eq. 0) then
!     open(33,file=outfile,form='unformatted') 
!   else
!     open(33,file=outfile,access='append',form='unformatted') 
!   endif
!   if(ichunk .eq. 0) then
!     open(34,file='deconvTb.bin',form='unformatted')
!   else
!     open(34,file='deconvTb.bin',access='append',form='unformatted') 
!   endif
  !end SJM 12/9/2014
  w10=dPRData%envSfcWind
  actOb=0

  dPRRet%convtb=-99
  scLonPR=-99
  scLatPR=-99
  tbMean=-99
  tbMax=-99
  tbMin=-99
  !call startprofs()
  rrate3D=-99
  tbNoOcean=-99
  pia35m=0.

  do j=1,dPRData%n1c21
     do i=1,49
        eLon=dPRData%xlon(i,j)
        eLat=dPRData%xlat(i,j)
        call getwfraction(eLat,&
             eLon,wfmap(i,j))
        wfmap(i,j)=wfmap(i,j)/100.
        
!begin  WSO 8/21/14 initialize ioquality flag to bad
       dPRData%ioqualityflagku(i, j) = dPRData%ioqualityflagku(i,j) + 900000
       if(i > 12 .and. i < 38) then
         dPRData%ioqualityflagdpr(i, j) = dPRData%ioqualityflagdpr(i, j) + 900000
       endif
!end    WSO 8/21/14
!  SFM  begin  04/16/2014; for M.Grecu, revision of nodes processing
        if(dPRData%xlon(i,j)>-998) then  !4/15/14 MG begin
        iLandSea=-2
        iGMI=-99  !4/22/14 MG
        jGMI=-99  !4/22/14 MG
!begin  WSO 3/16/17 initialize ig, jg
        ig = -99
        jg = -99
!end    WSO 3/16/17
        if(gmi2Grid%xmin>-998 .and. gmi2Grid%ymin>-998) then
           ig=(dPRData%xlon(i,j)-gmi2Grid%xmin)/gmi2Grid%dx+1
           jg=(dPRData%xlat(i,j)-gmi2Grid%ymin)/gmi2Grid%dx+1
           if(ig>0 .and. ig<gmi2Grid%nx .and. jg>0 .and. jg<gmi2Grid%ny) then
              iGMI=gmi2Grid%ig(ig,jg)
              jGMI=gmi2Grid%jg(ig,jg)
              if(gmi2Grid%actOb(ig,jg)==1) then
                 actOb(i,j)=1
                 iactob=iactob+1
              endif
              
           else
              if(jg>0 .and. jg<gmi2Grid%ny) then
                 dig=int(360/gmi2Grid%dx)
                 if(ig+dig>0 .and.  ig+dig<gmi2Grid%nx) then
                    ig=ig+dig
                    iGMI=gmi2Grid%ig(ig,jg)
                    jGMI=gmi2Grid%jg(ig,jg)
                 else
                    if(ig+2*dig>0 .and.  ig+2*dig<gmi2Grid%nx) then
                       ig=ig+2*dig
                       iGMI=gmi2Grid%ig(ig,jg)
                       jGMI=gmi2Grid%jg(ig,jg)
                    else
                       iGMI=-99
                       jGMI=-99
                    endif
                 endif
              endif
           endif
        else
           iGMI=-99
           jGMI=-99
        endif   ! 4/15/14 MG End
        !print*, iGMI, jGMI, i, j
!  SFM  end    04/16/2104
        k=1
        if(iGMI<0) then
           do while(iGMI<0 .and. k<8)
              if( ig+di(k)>0 .and. ig+di(k)<=gmi2Grid%nx .and.                 &
                   jg+dj(k)>0 .and. jg+dj(k)<=gmi2Grid%ny) then
                 iGMI=gmi2Grid%ig(ig+di(k),jg+dj(k))
                 jGMI=gmi2Grid%jg(ig+di(k),jg+dj(k))
                 
              endif
              k=k+1
           end do
        endif
        dPRData%ig(i,j)=iGMI
        dPRData%jg(i,j)=jGMI       
        if(jGMI>-99) then
           scLonPR(i,j)=gmiData%SCLon3(jGMI)
           scLatPR(i,j)=gmiData%SCLat3(jGMI)
           S1eiaPR(i,j)=gmidata%S1eia3(iGMI,jGMI)
           S2eiaPR(i,j)=gmidata%S2eia3(iGMI,jGMI)
           !print*,i,j,S1eiaPR(i,j), S2eiaPR(i,j)
        endif
        
        if(iGMI>0 .and. jGMI>0 .and. gmi2Grid%xmin>-99 ) then !4/22/14 MG
           if(gmiData%tpw3(iGMI,jGMI)>0) then
              if(gmiData%sfc_wind3(iGMI,jGMI)>0) then
                 radarRet%sfc_wind=gmiData%sfc_wind3(iGMI,jGMI)
                 tbRgrid(14,i,j+icL)=gmiData%tpw3(iGMI,jGMI) 
              else
                 radarRet%sfc_wind=dPRData%envSfcWind(i,j)
                 tbRgrid(14,i,j+icL)=-99
              endif

              if(ifdpr(1:1).ne.'Y'.and.iftest(1:1).ne.'Y') then
                 tbRgrid(1:9,i,j+icL)=gmiData%gmiS13(1:9,iGMI,jGMI)

!begin  WSO 8/21/2014 for at least some good data, re-set quality to nominal

                 if(maxval(tbRgrid(1:9, i, j + icL)) > 0) then
                   dPRData%ioqualityflagku(i, j) = dPRData%ioqualityflagku(i, j) - 900000
                   if(i > 12 .and. i < 38) then
                     dPRData%ioqualityflagdpr(i, j) = dPRData%ioqualityflagdpr(i, j) - 900000
                   endif
                 endif

              endif
           else
              if(ifdpr(1:1).ne.'Y'.and.iftest(1:1).ne.'Y') then
                 tbRgrid(1:9,i,j+icL)=gmiData%gmiS13(1:9,iGMI,jGMI)
                 
                 if(maxval(tbRgrid(1:9, i, j + icL)) > 0) then
                    dPRData%ioqualityflagku(i, j) = &
                         dPRData%ioqualityflagku(i, j) - 900000
                    if(i > 12 .and. i < 38) then
                       dPRData%ioqualityflagdpr(i, j) = &
                            dPRData%ioqualityflagdpr(i, j) - 900000
                    endif
                 endif


              endif
              radarRet%sfc_wind=dPRData%envSfcWind(i,j)
           endif
        endif

!begin  WSO 9/5/13 rename SRT PIA
        if(dPRData%rainType(i,j)>=100 ) then
!end    WSO 9/5/13 
           !print*, dPRData%node(:,i,j)
           dPRData%node(5,i,j)=dPRData%node(5,i,j)

!end    WSO 9/5/13 
           if (dPRData%rainType(i,j)/100==2) then
              dPRData%node(2,i,j)=dPRData%node(2,i,j)-1
           endif
           if (dPRData%rainType(i,j)/100==1) then
              dPRData%node(2,i,j)=dPRData%node(2,i,j)-1
              dPRData%node(4,i,j)=dPRData%node(4,i,j)+1
           endif
           stormStruct%nodes  = dPRData%node(:,i,j)
           radarData%z13obs   = dPRData%zku1c21(:,i,j)
           radarData%z35obs   = dPRData%zka1c21(:,i,j)
           !radarData%z35obs   = -99
!begin  WSO 9/5/13 rename SRT and DSRT PIA's and reliability factor
           radarData%pia13srt = dPRData%srtPIAku(i,j)
           radarData%relpia13srt = dPRData%srtrelPIAku(i,j)
           radarData%pia35srt = -99 
           radarData%pia35srt = dPRData%dsrtPIAka(i,j)
!end    WSO 9/5/13
           !begin SJM 7/10/14 add sigma-zero
           radarData%sigmaZeroKu = dPRData%sigmaZeroKu(i,j)
           radarData%sigmaZeroKa = dPRData%sigmaZeroKa(i,j)
           !print*, radarData%sigmaZeroKu, radarData%sigmaZeroKa
           !end SJM 7/10/14
           radarData%hfreez   = dPRData%freezH(i,j) /1000. 
           stormStruct%iSurf = dPRData%binRealSurface(i,j)+1

           do k=1,nbin
              if(radarData%z13obs(k)<-99) radarData%z13obs(k)=-99
           enddo

           !if(iGMI>0 .and. jGMI>0 ) then
           !   if(gmiData%tpw3(iGMI,jGMI)<0) then
           !      call landEmiss(emissout(i,j,1:9),dPRData,i,j,&
           !           dPRData%envSfcWind(i,j),dPRData%envSfcTemp(i,j))
           !   endif
           !endif

           stormStruct%rainType=dPRData%rainType(i,j)
           stormStruct%rainType=stormStruct%rainType/100
           
           itop=1
           radarRet%sfc_wind(1)=dPRData%envSfcWind(i,j)
           radarRet%sfc_windU(1)=dPRData%envSfcWindU(i,j)
           radarRet%sfc_windV(1)=dPRData%envSfcWindV(i,j)
           w10(i,j)=radarRet%sfc_wind(1)
           !print*, i,j,radarRet%sfc_wind(1),radarRet%sfc_windU(1),radarRet%sfc_windV(1)
           !w10(i,j)=0.5*(radarRet%sfc_wind(1)+dPRData%envSfcWind(i,j))
          
           if(dPRData%rainType(i,j)/100>=1) then
              do ibatch=1,1
                 eLon=dPRData%xlon(i,j)
                 eLat=dPRData%xlat(i,j)
                 call getwfraction(eLat,&
                      eLon,wfractPix)
                 call interpoldNw(i,j, logdNwf)

                 if(stormStruct%rainType == 1) then
                    do k=0,nmemb1-1 
                       radarRet%logdNw(k*9+1:(k+1)*9)=sysdn-0.1+ & !Sept 17, 2015 MG
                            logdNwf(k*9+1:(k+1)*9)
                    enddo
                 else
                    if(dPRData%rainType(i,j)/100==2) then
                       do k=0,nmemb1-1
!  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
!  SFM  begin  06/22/2014; for M.Grecu  (unknown justification)
                          radarRet%logdNw(k*9+1:(k+1)*9)=sysdn+0.0+            &
                               logdNwf(k*9+1:(k+1)*9)
!  SFM  end    06/22/2014
!  SFM  end    07/29/2014
                       enddo
                    else
                       do k=0,nmemb1-1
                          radarRet%logdNw(k*9+1:(k+1)*9)=sysdn+0.0+            &
                               logdNwf(k*9+1:(k+1)*9)
                       enddo
                    endif
                 endif
                 if(stormStruct%rainType .eq. 1 .and. dPRData%BBbin(i,j)<=0)   &
		 then
                    stormStruct%rainType = 3
                 endif
                 iRad=i
                 jRad=j
                 if(wfractPix<90 .and.  stormStruct%rainType ==1) then
                    do k=0,nmemb1-1 
                       radarRet%logdNw(k*9+1:(k+1)*9)=&
                            radarRet%logdNw(k*9+1:(k+1)*9)-0.1
                    enddo
                 endif
                 do k=0,nmemb1-1
                    xscalev(k+1)=1
                 enddo

		 
                 !radarRet%sfc_wind(1:nmemb1)=&
                 !    2*windPert(1:nmemb1)
                 

            
                      
                 !begin SJM 10/16/15
                 !print*, dPRData%snowIceCover(i,j)
                 jj=2880-floor((elat+90.)/180.*2880.)
                 if(jj .lt. 1) jj=1
                 if(jj .gt. 2880) jj = 2880
                 ii=floor((elon+180.)/360.*5760.)+1
                 if(ii .lt. 1) ii = 1
                 if(ii .gt. 5760) ii = 5760
                 !print*, elat,elon,ii,jj
                 if(dPRData%snowIceCover(i,j) .eq. 0) then
                   stype = LUT%land_class_map_bare(ii,jj) !SJM 9/9/15
                 else
                   stype = LUT%land_class_map_snow(ii,jj) !SJM 9/9/15
                 endif
                 !print*, eLon, eLat, stype
                 !begin SJM 12/9/2014
                 !Determine emissivity and surface backscatter here. Surface classification will be as follows:
                 !If stype=1 (GPROF water class), use water
                 !If stype=2 (GPROF sea ice), use sea ice
                 !For all other stype values, if wfract is greater than 50, use a weighted mix of the stype class and water
                 !If wfract is < 50, just use the stype class as is.
                 !Perturb wind in all pixels
                 w10_min = 9999.9
                 w10_max = -9999.9
                 do ii=max(1,i-25),min(49,i+25)
                   do jj=max(1,j-25),min(dPRdata%n1c21,j+25)
                     if((ii-i)**2+(jj-j)**2 .gt. 25**2) cycle
                     if(dPRData%envsfcWind(ii,jj) .ge. 0. .and. dPRData%envsfcWind(ii,jj) .lt. w10_min) w10_min = dPRData%envsfcWind(ii,jj)
                     if(dPRData%envsfcWind(ii,jj) .ge. 0. .and. dPRData%envsfcWind(ii,jj) .gt. w10_max) w10_max = dPRData%envsfcWind(ii,jj)
                   end do
                 end do
                 !radarRet%sfc_wind(1:nmemb1)=dPRData%envsfcWind(i,j)+max(12.5,w10_max-w10_min,dPRData%envsfcWind(i,j))*windPertU(1:nmemb1)
                 radarRet%sfc_windU(1:nmemb)=dprData%envSfcWindU(i,j)+max(12.5,w10_max-w10_min,dPRData%envsfcWind(i,j))*windPertU(1:nmemb1)
                 !radarRet%sfc_windU(1:nmemb)=7.5+25.*windPertU(1:nmemb1)
                 radarRet%sfc_windV(1:nmemb)=dprData%envSfcWindV(i,j)+max(12.5,w10_max-w10_min,dPRData%envsfcWind(i,j))*windPertV(1:nmemb1)
                 !radarRet%sfc_windV(1:nmemb)=0.+25.*windPertV(1:nmemb1)
                 !where(radarRet%sfc_wind .lt. 0.) radarRet%sfc_wind = 0.
                 !calculate TPW to choose land class EOFS (8 or 10-channel)
                 tpw_ij=0.
                 do k = 1, dprData%binRealSurface(i,j)
                   !there are many approximations in this calc but for purposes of choosing 8- or 10-channel EOFs it is ok
                   if(dprData%envQv(k,i,j) .gt. 0.) tpw_ij=tpw_ij+250.*100.*dprData%EnvPress(k,i,j)/(dprData%EnvTemp(k,i,j)*(461.5+286.9/(0.001*dprData%envQv(k,i,j))))
                 end do
                 
                 if(stype .eq. 1) then
                   do k=0,nmemb1-1
                     radarRet%sfc_wind(k+1) = sqrt(radarRet%sfc_windU(k+1)**2+radarRet%sfc_windV(k+1)**2)
                     call calc_relAz(DPRData%sclon(j), DPRData%sclat(j), elon, elat, radarRet%sfc_windU(k+1), radarRet%sfc_windV(k+1), relAz)
                     !relAz=0.
                     !print*, radarRet%sfc_windU(k+1), radarRet%sfc_windV(k+1), relAz
                     call intplte_water_sigma0(i,radarRet%sfc_wind(k+1),relAz,s0Ku, s0Ka, s0stdKu, s0stdKa, s0corr)
                     ds0Ku = normal2(0.,1.)
                     ds0Ka = (s0corr**2)*ds0Ku+(1.-s0corr**2)*normal2(0.,1.)
                     ds0Ku = ds0Ku*s0stdKu
                     if(s0Ka .ne. -99.) ds0Ka = ds0Ka*s0stdKa
                     !print*, i, k, radarRet%sfc_wind(k),s0Ku, s0Ka
                     !add (correlated) noise to obs
                     radarRet%simSigmaZeroKu(k+1) = s0Ku+ds0Ku!-radarRet%pia13(k)!+ds0Ku
                     if(s0Ka .ne. -99.) then
                       radarRet%simSigmaZeroKa(k+1) = s0Ka+ds0Ka!-radarRet%pia35(k)!+ds0Ka !need to add multiple scattering contribution as well
                     else
                       radarRet%simSigmaZeroKa(k+1) = -99.
                     endif
                     sigmaZeroVarKu(i,j) = s0stdKu**2
                     if(s0Ka .ne. -99.) sigmaZeroVarKa(i,j) = s0stdKa**2
                     if(s0Ka .ne. -99.) sigmaZeroCov(i,j) = s0corr*s0stdKu*s0stdKa
                     call calc_relAz(scLonPR(i,j), scLatPR(i,j), elon, elat, radarRet%sfc_windU(k+1), radarRet%sfc_windV(k+1), relAz)
                     !relAz=0.
                     !print*, i,j,radarRet%sfc_wind(k+1), '1'
                     
                     do ii=1,5
                       call intplte_emis(ii,0,dPRData%envSfcTemp(i,j),radarRet%sfc_wind(k+1),relAz,S1eiaPR(i,j),emis,ebar)
                       emissv(ii)=emis
                       call intplte_emis(ii,1,dPRData%envSfcTemp(i,j),radarRet%sfc_wind(k+1),relAz,S1eiaPR(i,j),emis,ebar)
                       emissh(ii)=emis
                     end do
                     
                     do ii=6,6
                       call intplte_emis(ii,0,dPRData%envSfcTemp(i,j),radarRet%sfc_wind(k+1),relAz,S2eiaPR(i,j),emis,ebar)
                       emissv(ii)=emis
                       call intplte_emis(ii,1,dPRData%envSfcTemp(i,j),radarRet%sfc_wind(k+1),relAz,S2eiaPR(i,j),emis,ebar)
                       emissh(ii)=emis
                     end do
                     !print '(15F8.3)', relAz,S1eiaPR(i,j),S2eiaPR(i,j), emissv, emissh
                     radarRet%emis(k*2*nmfreq2+1) = emissv(1)
                     radarRet%emis(k*2*nmfreq2+2) = emissh(1)
                     radarRet%emis(k*2*nmfreq2+3) = emissv(2)
                     radarRet%emis(k*2*nmfreq2+4) = emissh(2)
                     radarRet%emis(k*2*nmfreq2+5) = emissv(3)
                     radarRet%emis(k*2*nmfreq2+6) = emissv(4)
                     radarRet%emis(k*2*nmfreq2+7) = emissh(4)
                     radarRet%emis(k*2*nmfreq2+8) = emissv(5)
                     radarRet%emis(k*2*nmfreq2+9) = emissh(5)
                     radarRet%emis(k*2*nmfreq2+10) = emissv(6)
                     radarRet%emis(k*2*nmfreq2+11) = emissh(6)
                     radarRet%emis(k*2*nmfreq2+12:k*2*nmfreq2+16) = radarRet%emis(k*2*nmfreq2+10)
                     !print '(I3,14F8.3,2F8.2)', k, radarRet%sfc_wind(k+1), radarRet%emis(k*2*nmfreq+1:k*2*nmfreq+13), radarRet%simSigmaZeroKu(k+1), radarRet%simSigmaZeroKa(k+1)
                   end do
                   
                   
                   !if(max(w10_max-w10_min,dPRData%envsfcWind(i,j)) .ge. 10.) radarRet%sfc_wind(1:nmemb1)=dPRData%envsfcWind(i,j)+max(w10_max-w10_min,dPRData%envsfcWind(i,j))*windPert(1:nmemb1)
                 else 

                   call getemiss(elat,elon,dPRData%snowIceCover(i,j),emissv,emissh,emissv_std,emissh_std)
                   !print '(2F8.2,I5, 12F8.3)', elat, elon, stype, emissv, emissh
                   !stop
                   !call get_s0(elat,elon,i,stype,s0Ku,s0Ka, s0stdKu, s0stdKa)
                   !alternative: get reference sigma_zero from observed sigma_zero +srt PIA
                   s0Ku = DPRData%sigmaZeroKu(i,j)+DPRData%srtPIAKu(i,j)
                   s0stdKu = DPRData%srtsigmaPIAKu(i,j)
                   s0Ka = DPRData%sigmaZeroKa(i,j)+DPRData%dsrtPIAKa(i,j)
                   s0stdKa = DPRData%dsrtsigmaPIAKa(i,j)
                   !print '(2F8.2,I5,11F8.2)', elat, elon, stype, s0Ku, s0stdKu,LUT%land_class_sigma0Ku_std(stype-1,i), DPRData%sigmaZeroKu(i,j),DPRData%srtsigmaPIAKu(i,j),DPRData%dsrtsigmaPIAKu(i,j), s0Ka, s0stdKa, DPRData%sigmaZeroKa(i,j),DPRData%dsrtsigmaPIAKa(i,j),LUT%land_class_sigma0Ka_std(stype-1,i)
                   !reset PIA for forward model
                   !radarData%pia13srt = s0Ku-DPRData%sigmaZeroKu(i,j)
                   !radarData%pia35srt = s0Ka-DPRData%sigmaZeroKa(i,j)
                   sigmaZeroVarKu(i,j) = s0stdKu**2
                   sigmaZeroVarKa(i,j) = s0stdKa**2
                   sigmaZeroCov(i,j) = 0.5*s0stdKu*s0stdKa !need to replace w/ actual correlation; this is a conservative value
                   !use same sfc reference as DPR
                   !if(dPRData%snrRatioku(i, j) > 2.) s0Ku = dPRData%sigmaZeroKu(i,j)+dPRData%srtPIAKu(i,j)
                   !if(dPRData%snrRatioka(i, j) > 2.) s0Ka = dPRData%sigmaZeroKa(i,j)+dPRData%dsrtPIAKa(i,j)
                   !Set ensemble emissivities
                   do k=0,nmemb1-1
                     radarRet%emis(k*2*nmfreq2+1) = emissv(1)
                     radarRet%emis(k*2*nmfreq2+2) = emissh(1)
                     radarRet%emis(k*2*nmfreq2+3) = emissv(2)
                     radarRet%emis(k*2*nmfreq2+4) = emissh(2)
                     radarRet%emis(k*2*nmfreq2+6) = emissv(4)
                     radarRet%emis(k*2*nmfreq2+7) = emissh(4)
                     radarRet%emis(k*2*nmfreq2+8) = emissv(5)
                     radarRet%emis(k*2*nmfreq2+9) = emissh(5)
                     radarRet%emis(k*2*nmfreq2+10) = emissv(6)
                     radarRet%emis(k*2*nmfreq2+11) = emissh(6)
                     radarRet%simSigmaZeroKu(k+1) = s0Ku
                     radarRet%simSigmaZeroKa(k+1) = s0Ka
                     if(tpw_ij .gt. 10.) then
                       do ii = 1,10
                         radarRet%emis(k*2*nmfreq2+1) = radarRet%emis(k*2*nmfreq2+1) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_LF(stype-1,ii,1)
                         radarRet%emis(k*2*nmfreq2+2) = radarRet%emis(k*2*nmfreq2+2) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_LF(stype-1,ii,2)
                         radarRet%emis(k*2*nmfreq2+3) = radarRet%emis(k*2*nmfreq2+3) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_LF(stype-1,ii,3)
                         radarRet%emis(k*2*nmfreq2+4) = radarRet%emis(k*2*nmfreq2+4) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_LF(stype-1,ii,4)
                         radarRet%emis(k*2*nmfreq2+6) = radarRet%emis(k*2*nmfreq2+6) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_LF(stype-1,ii,5)
                         radarRet%emis(k*2*nmfreq2+7) = radarRet%emis(k*2*nmfreq2+7) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_LF(stype-1,ii,6)
                         radarRet%emis(k*2*nmfreq2+8) = radarRet%emis(k*2*nmfreq2+8) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_LF(stype-1,ii,7)
                         radarRet%emis(k*2*nmfreq2+9) = radarRet%emis(k*2*nmfreq2+9) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_LF(stype-1,ii,8)
                         !radarRet%emis(k*2*nmfreq2+10) = radarRet%emis(k*2*nmfreq2+10) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_LF(stype-1,ii,9)
                         !radarRet%emis(k*2*nmfreq2+11) = radarRet%emis(k*2*nmfreq2+11) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_LF(stype-1,ii,10)
                         radarRet%simSigmaZeroKu(k+1) = radarRet%simSigmaZeroKu(k+1)+emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_LF(stype-1,ii,10+i)*s0stdKu/LUT%land_class_sigma0Ku_std(stype-1,i)
                         radarRet%simSigmaZeroKa(k+1) = radarRet%simSigmaZeroKa(k+1)+emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_LF(stype-1,ii,10+49+i)*s0stdKa/LUT%land_class_sigma0Ka_std(stype-1,i)
                         !print*, emis_eofs(k+1,ii), LUT%land_class_emis_eofs_NS(stype,ii,1)
                       end do
                     else
                       do ii = 1,12
                         radarRet%emis(k*2*nmfreq2+1) = radarRet%emis(k*2*nmfreq2+1) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_AF(stype-1,ii,1)
                         radarRet%emis(k*2*nmfreq2+2) = radarRet%emis(k*2*nmfreq2+2) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_AF(stype-1,ii,2)
                         radarRet%emis(k*2*nmfreq2+3) = radarRet%emis(k*2*nmfreq2+3) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_AF(stype-1,ii,3)
                         radarRet%emis(k*2*nmfreq2+4) = radarRet%emis(k*2*nmfreq2+4) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_AF(stype-1,ii,4)
                         radarRet%emis(k*2*nmfreq2+6) = radarRet%emis(k*2*nmfreq2+6) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_AF(stype-1,ii,5)
                         radarRet%emis(k*2*nmfreq2+7) = radarRet%emis(k*2*nmfreq2+7) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_AF(stype-1,ii,6)
                         radarRet%emis(k*2*nmfreq2+8) = radarRet%emis(k*2*nmfreq2+8) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_AF(stype-1,ii,7)
                         radarRet%emis(k*2*nmfreq2+9) = radarRet%emis(k*2*nmfreq2+9) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_AF(stype-1,ii,8)
                         radarRet%emis(k*2*nmfreq2+10) = radarRet%emis(k*2*nmfreq2+10) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_AF(stype-1,ii,9)
                         radarRet%emis(k*2*nmfreq2+11) = radarRet%emis(k*2*nmfreq2+11) + emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_AF(stype-1,ii,10)
                         radarRet%simSigmaZeroKu(k+1) = radarRet%simSigmaZeroKu(k+1)+emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_AF(stype-1,ii,10+i)*s0stdKu/LUT%land_class_sigma0Ku_std(stype-1,i)
                         radarRet%simSigmaZeroKa(k+1) = radarRet%simSigmaZeroKa(k+1)+emis_eofs(k+1,ii)*LUT%land_class_emis_eofs_AF(stype-1,ii,10+49+i)*s0stdKa/LUT%land_class_sigma0Ka_std(stype-1,i)
                         !print*, emis_eofs(k+1,ii), LUT%land_class_emis_eofs_NS(stype,ii,1)
                       end do
                     endif
                     radarRet%emis(k*2*nmfreq2+5) = 0.64*radarRet%emis(k*2*nmfreq2+4)+0.36*radarRet%emis(k*2*nmfreq2+6)
                     radarRet%emis(k*2*nmfreq2+12:k*2*nmfreq2+16) = radarRet%emis(k*2*nmfreq2+10)
                     !print '(I3,13F8.3,2F8.2)', k, radarRet%emis(k*2*nmfreq+1:k*2*nmfreq+13), radarRet%simSigmaZeroKu(k+1), radarRet%simSigmaZeroKa(k+1)
                   end do
                   
                   !print*, s0Ku, dPRData%sigmaZeroKu(i,j)+dPRData%srtPIAKu(i,j), dPRData%sigmaZeroKa(i,j)+dPRData%dsrtPIAKa(i,j)
                   if(stype .gt. 2 .and. wfractPix .gt. 1) then
                     do k=0,nmemb1-1
                       radarRet%sfc_wind(k+1) = sqrt(radarRet%sfc_windU(k+1)**2+radarRet%sfc_windV(k+1)**2)
                       !relAz=0.
                       call calc_relAz(DPRData%sclon(j), DPRData%sclat(j), elon, elat, radarRet%sfc_windU(k+1), radarRet%sfc_windV(k+1), relAz)
                       call intplte_water_sigma0(i,radarRet%sfc_wind(k+1),relAz,s0Ku, s0Ka, s0stdKu, s0stdKa, s0corr)
                       ds0Ku = normal2(0.,1.)
                       ds0Ka = (s0corr**2)*ds0Ku+(1.-s0corr**2)*normal2(0.,1.)
                       ds0Ku = ds0Ku*s0stdKu
                       if(s0Ka .ne. -99.) ds0Ka = ds0Ka*s0stdKa
                       radarRet%simSigmaZeroKu(k+1) = (wfractPix/100.*s0Ku+ds0Ku)+(1.-wfractPix/100.)*radarRet%simSigmaZeroKu(k+1)!-radarRet%pia13(k)!+ds0Ku
                       radarRet%simSigmaZeroKa(k+1) = (wfractPix/100.*s0Ka+ds0Ka)+(1.-wfractPix/100.)*radarRet%simSigmaZeroKa(k+1)!-radarRet%pia35(k)!+ds0Ka !need to add multiple scattering contribution as well
                       sigmaZeroVarKu(i,j) = (wfractPix/100.*s0stdKu)+(1.-wfractPix/100.)*sigmaZeroVarKu(i,j)
                       sigmaZeroVarKa(i,j) = (wfractPix/100.*s0stdKa)+(1.-wfractPix/100.)*sigmaZeroVarKa(i,j)
                       call calc_relAz(scLonPR(i,j), scLatPR(i,j), elon, elat, radarRet%sfc_windU(k+1), radarRet%sfc_windV(k+1), relAz)
                       !relAz=0.
                       !print*, i,j,radarRet%sfc_wind(k+1),'2'
                       do ii=1,5
                         call intplte_emis(ii,0,dPRData%envSfcTemp(i,j),radarRet%sfc_wind(k+1),relAz,S1eiaPR(i,j),emis,ebar)
                         emissv(ii)=emis
                         call intplte_emis(ii,1,dPRData%envSfcTemp(i,j),radarRet%sfc_wind(k+1),relAz,S1eiaPR(i,j),emis,ebar)
                         emissh(ii)=emis
                       end do
                       do ii=6,6
                         call intplte_emis(ii,0,dPRData%envSfcTemp(i,j),radarRet%sfc_wind(k+1),relAz,S2eiaPR(i,j),emis,ebar)
                         emissv(ii)=emis
                         call intplte_emis(ii,1,dPRData%envSfcTemp(i,j),radarRet%sfc_wind(k+1),relAz,S2eiaPR(i,j),emis,ebar)
                         emissh(ii)=emis
                       end do
                       
                       radarRet%emis(k*2*nmfreq2+1) = emissv(1)*wfractPix/100.+(1.-wfractPix/100.)*radarRet%emis(k*2*nmfreq2+1)
                       radarRet%emis(k*2*nmfreq2+2) = emissh(1)*wfractPix/100.+(1.-wfractPix/100.)*radarRet%emis(k*2*nmfreq2+2)
                       radarRet%emis(k*2*nmfreq2+3) = emissv(2)*wfractPix/100.+(1.-wfractPix/100.)*radarRet%emis(k*2*nmfreq2+3)
                       radarRet%emis(k*2*nmfreq2+4) = emissh(2)*wfractPix/100.+(1.-wfractPix/100.)*radarRet%emis(k*2*nmfreq2+4)
                       radarRet%emis(k*2*nmfreq2+5) = emissv(3)*wfractPix/100.+(1.-wfractPix/100.)*radarRet%emis(k*2*nmfreq2+5)
                       radarRet%emis(k*2*nmfreq2+6) = emissv(4)*wfractPix/100.+(1.-wfractPix/100.)*radarRet%emis(k*2*nmfreq2+6)
                       radarRet%emis(k*2*nmfreq2+7) = emissh(4)*wfractPix/100.+(1.-wfractPix/100.)*radarRet%emis(k*2*nmfreq2+7)
                       radarRet%emis(k*2*nmfreq2+8) = emissv(5)*wfractPix/100.+(1.-wfractPix/100.)*radarRet%emis(k*2*nmfreq2+8)
                       radarRet%emis(k*2*nmfreq2+9) = emissh(5)*wfractPix/100.+(1.-wfractPix/100.)*radarRet%emis(k*2*nmfreq2+9)
                       radarRet%emis(k*2*nmfreq2+10) = emissv(6)*wfractPix/100.+(1.-wfractPix/100.)*radarRet%emis(k*2*nmfreq2+10)
                       radarRet%emis(k*2*nmfreq2+11) = emissh(6)*wfractPix/100.+(1.-wfractPix/100.)*radarRet%emis(k*2*nmfreq2+11)
                     end do
                   else
                     radarRet%sfc_wind(1:nmemb1)=dPRData%envsfcWind(i,j) !remove wind perturbations over land
                   endif
                   
                 endif
                 !end SJM 10/16/15
                 call setEnv(dPRData%envQv(:,i,j),dPRData%envTemp(:,i,j),&
                      dPRData%envPress(:,i,j),dPRData%envSfcTemp(i,j),&
                      dPRData%envSknTemp(i,j))
!begin  MG 10/29/15 assign reliability flag
                 reliabFlag=dPRData%NSRelibFlag(i,j)
!end    MG 10/20/15
                 call  ensRadRetStCvKu(radarData,stormStruct,                  &
                      retParam, nmu2,radarRet, itop, rms1, rms2, sysdN, iit, &
                      xscalev, randemiss, dPRData%localZenithAngle(i,j), &
                      wfractPix, ichunk, i, j, dZms(i,j), msFlag(i, j)) !! MS&WSO addition Feb 11, 2017

!begin WSO 2/11/17 assign dZms to output variable
                 multiscatsurface_MS(i, j) = dZms(i, j)
                 multiscatcalc_MS(i, j) = msFlag(i, j)
!end   WSO 2/11/17
                 
                 !begin SJM 10/16/15
                      !generate simulated sigma_zero over water
                      !print*, wfractPix
                 
                 if(stype .ne. 1) then
                    call getemiss(elat,elon,dPRData%snowIceCover(i,j),emissv,emissh,emissv_std,emissh_std)
                    emissoutL(i,j,1:13)=(/emissv(1),emissh(1),&
                         emissv(2),emissh(2),&
                         emissv(3),emissv(4),emissh(4),&
                         emissv(5),emissh(5),emissv(6),emissh(6),emissv(6),emissv(6)/)
                    if(maxval(emissoutL(i,j,1:13)) .gt. 2.) print'(A6,2I5,13F8.3)', 'EnsRet', i,j,emissoutL(i,j,1:13)
                 endif
!                  if(j>1 .and. j<dPRData%n1c21 .and. pia35m(i,j)>0.1) then
!                    nubfc= 1./srtpiaf(pia35m(i-1:i+1,j-1:j+1),9)
!                  else
!                    nubfc=1.
!                  endif
!                  print*, i,j, nubfc
                 do k=0,nmemb1-1
                   radarRet%simSigmaZeroKu(k+1) = radarRet%simSigmaZeroKu(k+1)-radarRet%pia13(k+1)!+ds0Ku
                   !if(radarRet%simSigmaZeroKa(k+1) .ne. -99.) radarRet%simSigmaZeroKa(k+1) = radarRet%simSigmaZeroKa(k+1)-radarRet%pia35(k+1)/nubfc!+ds0Ka 
                   !print '(I3,13F8.3,2F8.2)', k, radarRet%tb(k*2*nmfreq+1:k*2*nmfreq+13), radarRet%simSigmaZeroKu(k+1), radarRet%simSigmaZeroKa(k+1)
                 end do
                 !end sjm 10/16/15
                 dPRRet%n9(:,i,j)=n9+1

                 dPRRet%cldwcoeff(i,j,:,:)=cldwcoeff(1:10,1:nmemb)


                 if(radarRet%tb(2)>0) ntbpix=ntbpix+1
                 
                 do k=0,nmemb1-1
                    iy(k+1)=k+1
                 enddo
                 do k=0,nmemb1-1
                    dPRRet%tb(i,j,2,:,k+1+(ibatch-1)*nmemb1)=                  &
                         radarRet%tb((iy(k+1)-1)*2*radarRet%nmfreq+1:          &
                         2*(iy(k+1)-1)*radarRet%nmfreq+radarRet%nmfreq)
                    dPRRet%tb(i,j,1,:,k+1+(ibatch-1)*nmemb1)=                  &
                         radarRet%tb(((iy(k+1)-1)*2+1)*radarRet%nmfreq+1:      &
                         2*((iy(k+1)-1)+1)*radarRet%nmfreq)
                    dPRRet%emtb(i,j,2,:,k+1+(ibatch-1)*nmemb1)=           &
                         radarRet%emtb((iy(k+1)-1)*2*radarRet%nmfreq+1:    &
                         2*(iy(k+1)-1)*radarRet%nmfreq+radarRet%nmfreq)
                    dPRRet%emtb(i,j,1,:,k+1+(ibatch-1)*nmemb1)=             &
                         radarRet%emtb(((iy(k+1)-1)*2+1)*radarRet%nmfreq+1:   &
                         2*((iy(k+1)-1)+1)*radarRet%nmfreq)
                    !begin SJM 10/16/15
                    dPRRet%emis(i,j,1,:,k+1+(ibatch-1)*nmemb1)=             &
                         radarRet%emis(((iy(k+1)-1)*2+1)*radarRet%nmfreq+1:   &
                         2*((iy(k+1)-1)+1)*radarRet%nmfreq)
                    dPRRet%emis(i,j,2,:,k+1+(ibatch-1)*nmemb1)=           &
                         radarRet%emis((iy(k+1)-1)*2*radarRet%nmfreq+1:    &
                         2*(iy(k+1)-1)*radarRet%nmfreq+radarRet%nmfreq)
                    !end SJM 10/16/15
                    dPRRet%log10dNw (k+1+(ibatch-1)*nmemb1,:,i,j)=-99
                    dPRRet%d0 (k+1+(ibatch-1)*nmemb1,:,i,j)=-99
                 enddo
!                 print*, dPRRet%tb(i,j,1,:,1)
!  ifreqG(1:9)=(/1,1,2,2,3,4,4,5,5/)
!  ipolG(1:9)=(/1,2,1,2,1,1,2,1,2/)
                 do ik=1,9
                    tbMean(i,j,ik)=&
                         sum(dPRRet%tb(i,j,ipolG(ik),ifreqG(ik),1:nmemb1))/&
                         nmemb1
                    tbMin(i,j,ik)=&
                         minval(dPRRet%tb(i,j,ipolG(ik),ifreqG(ik),1:nmemb1))
                    tbMax(i,j,ik)=&
                         maxval(dPRRet%tb(i,j,ipolG(ik),ifreqG(ik),1:nmemb1))
                    tbNoOcean(i,j,ik)=tbMean(i,j,ik)
                 enddo
                 
                 do ik=1,9
                    do jk=1,9
                       covTb(i,j,ik,jk)=&
                            covar(dPRRet%tb(i,j,ipolG(ik),ifreqG(ik),1:nmemb1),&
                            dPRRet%tb(i,j,ipolG(jk),ifreqG(jk),1:nmemb1),nmemb1)
                    enddo
                   
                 enddo
121 format(81(F8.2,1x))
                 if(wfractPix>90) then
                    !write(*,121) wfmap(i,j), tbMean(i,j,:)
                 endif
                 if(ifdpr(1:1)=='Y') then
                    do k=1,9
                       !tbRgrid(k,i,j+icL)= dPRRet%tb(i,j,ipolG(k),ifreqG(k),nmemb1)
                    enddo
                 endif

                 do k=0,nmemb1-1
                    dPRRet%log10dNw (k+1+(ibatch-1)*nmemb1,1:ngates,i,j)=      &
                         radarRet%log10dNw((iy(k+1)-1)*ngates+1:(iy(k+1))*ngates)
                    dPRRet%d0 (k+1+(ibatch-1)*nmemb1,1:ngates,i,j)=            &
                         radarRet%d0((iy(k+1)-1)*ngates+1:(iy(k+1))*ngates)
                    dPRRet%rrate (k+1+(ibatch-1)*nmemb1,1:ngates,i,j)=         &
                         radarRet%rrate((iy(k+1)-1)*ngates+1:(iy(k+1))*ngates)
                    dPRRet%pwc (k+1+(ibatch-1)*nmemb1,1:ngates,i,j)=           &
                         (radarRet%pwc((iy(k+1)-1)*ngates+1:(iy(k+1))*ngates))
                    dPRRet%z13c(k+1+(ibatch-1)*nmemb1,1:ngates,i,j)=           &
                         radarRet%z13c((iy(k+1)-1)*ngates+1:(iy(k+1))*ngates)
                    dPRRet%z35mod0(k+1+(ibatch-1)*nmemb1,1:ngates,i,j)=        &
                         radarRet%z35mod0((iy(k+1)-1)*ngates+1:(iy(k+1))*ngates)
                    dPRRet%z35(k+1+(ibatch-1)*nmemb1,1:ngates,i,j)=            &
                         radarRet%z35((iy(k+1)-1)*ngates+1:(iy(k+1))*ngates)
                 enddo
                 
                 meansfcRain=0
                 do k=0,nmemb1-1
                    dPRRet%sfcRainEns(i,j,k+1+(ibatch-1)*nmemb1)=              &
                         radarRet%rrate((iy(k+1)-1)*ngates+dPRData%node(5,i,j))
                    dPRRet%sfcd0Ens(i,j,k+1+(ibatch-1)*nmemb1)=                &
                         radarRet%d0((iy(k+1)-1)*ngates+1+dPRData%node(5,i,j))
                    meansfcRain=meansfcRain+ dPRRet%sfcRainEns(i,j,k+1)
                    dPRRet%sfcNwEns(i,j,k+1+(ibatch-1)*nmemb1)=                &
                         radarRet%log10dNw((iy(k+1)-1)*ngates+1                &
                         +dPRData%node(5,i,j))
                    dPRRet%sfcWindEns(i,j,k+1+(ibatch-1)*nmemb1)=radarRet%sfc_wind(k+1) !SJM 12/4/2014
                    dPRRet%pia13mod(i,j,k+1+(ibatch-1)*nmemb1)=                &
                         radarRet%pia13(iy(k+1))
                    dPRRet%pia35mod(i,j,k+1+(ibatch-1)*nmemb1)=&
                         radarRet%pia35(iy(k+1))
                    dPRRet%simSigmaZeroKu(i,j,k+1+(ibatch-1)*nmemb1)=                & !SJM 2/4/2015
                         radarRet%simSigmaZeroKu(iy(k+1))
                    dPRRet%simSigmaZeroKa(i,j,k+1+(ibatch-1)*nmemb1)=                & !SJM 2/4/2015
                         radarRet%simSigmaZeroKa(iy(k+1))
                 enddo
                 pia13m=sum(dPRRet%pia13mod(i,j,1:nmemb1))/nmemb1
                 pia35m(i,j)=sum(dPRRet%pia35mod(i,j,1:nmemb1))/nmemb1
                 pia13s=(sum((dPRRet%pia13mod(i,j,1:nmemb1)-pia13m)**2)/&
                      nmemb1)**.5
                 sfcRain(i,j)=meansfcRain/nmemb1

                 sfcRainStd(i,j)=sqrt(sum((dPRRet%sfcRainEns(i,j,1:nmemb1)- &
                      sfcRain(i,j))**2)/(nmemb1-1))
                 piaOut(i,j)=pia13m
              enddo
           endif
        else
           do i1=1,5
              do j1=1,5
                 eLon=dPRData%xlon(i,j)
                 eLat=dPRData%xlat(i,j)
                 call getwfraction(eLat+(i1-2)*0.15,&
                      eLon+(j1-2)*0.15,wfract(i1,j1))
              enddo
           enddo
           
           wfractm=sum(wfract(:,:))/25.
           wfractsd=sqrt(sum((wfract(:,:)-wfractm)**2)/25)
           !begin SJM 10/16/15
           !stype = LUT%land_class_map(mod(floor((dPRData%xlon(i,j)+180.)/360.*5760.),5760)+1, 2880-floor((dPRData%xlat(i,j)+90.)/180.*2880.)) !SJM 9/9/15
           jj=2880-floor((DPRData%xlat(i,j)+90.)/180.*2880.)
           if(jj .lt. 1) jj=1
           if(jj .gt. 2880) jj = 2880
           ii=floor((DPRData%xlon(i,j)+180.)/360.*5760.)+1
           if(ii .lt. 1) ii = 1
           if(ii .gt. 5760) ii = 5760
           if(dPRData%snowIceCover(i,j) .eq. 0) then
             stype = LUT%land_class_map_bare(ii,jj) !SJM 9/9/15
           else
             stype = LUT%land_class_map_snow(ii,jj) !SJM 9/9/15
           endif
           if(stype .eq. 1) then
             !print*, i,j, radarRet%sfc_wind(k+1),'3'
             call calc_relAz(scLonPR(i,j), scLatPR(i,j), dPRData%xlon(i,j), dPRData%xlat(i,j), dprData%envSfcWindU(i,j), dprData%envSfcWindV(i,j), relAz)
             do ii=1,5
               call intplte_emis(ii,0,dPRData%envSfcTemp(i,j),dprData%envSfcWind(i,j),relAz,S1eiaPR(i,j),emis,ebar)
               emissv(ii)=emis
               call intplte_emis(ii,1,dPRData%envSfcTemp(i,j),dprData%envSfcWind(i,j),relAz,S1eiaPR(i,j),emis,ebar)
               emissh(ii)=emis
             end do
             do ii=6,6
               call intplte_emis(ii,0,dPRData%envSfcTemp(i,j),dprData%envSfcWind(i,j),relAz,S2eiaPR(i,j),emis,ebar)
               emissv(ii)=emis
               call intplte_emis(ii,1,dPRData%envSfcTemp(i,j),dprData%envSfcWind(i,j),relAz,S2eiaPR(i,j),emis,ebar)
               emissh(ii)=emis
             end do
!begin  SJM 2/14/17 emissivity standard deviation code
             emissv_std = missing_r4
             emissh_std = missing_r4
           else if(stype .gt. 2 .and. wfract(2,2) .gt. 1) then !SJM 2/15/17 index changes
             call getemiss(dPRData%xlat(i,j),dPRData%xlon(i,j),dPRData%snowIceCover(i,j),emissv,emissh,emissv_std,emissh_std)
             call calc_relAz(scLonPR(i,j), scLatPR(i,j), dPRData%xlon(i,j), dPRData%xlat(i,j), dprData%envSfcWindU(i,j), dprData%envSfcWindV(i,j), relAz)
             do ii=1,5
               call intplte_emis(ii,0,dPRData%envSfcTemp(i,j),dprData%envSfcWind(i,j),relAz,S1eiaPR(i,j),emis,ebar)
               emissv(ii)=0.01*wfract(2,2)*emis+(1.-0.01*wfract(2,2))*emissv(ii) !SJM 2/15/17 index changes
               call intplte_emis(ii,1,dPRData%envSfcTemp(i,j),dprData%envSfcWind(i,j),relAz,S1eiaPR(i,j),emis,ebar)
               emissh(ii)=0.01*wfract(2,2)*emis+(1.-0.01*wfract(2,2))*emissh(ii) !SJM 2/15/17 index changes
             end do
             do ii=6,6
               call intplte_emis(ii,0,dPRData%envSfcTemp(i,j),dprData%envSfcWind(i,j),relAz,S2eiaPR(i,j),emis,ebar)
               emissv(ii)=0.01*wfract(2,2)*emis+(1.-0.01*wfract(2,2))*emissv(ii) !SJM 2/15/17 index changes
               call intplte_emis(ii,1,dPRData%envSfcTemp(i,j),dprData%envSfcWind(i,j),relAz,S2eiaPR(i,j),emis,ebar)
               emissh(ii)=0.01*wfract(2,2)*emis+(1.-0.01*wfract(2,2))*emissh(ii) !SJM 2/15/17 index changes
             end do
!end  SJM 2/14/17
           else
              call getemiss(dPRData%xlat(i,j),dPRData%xlon(i,j),dPRData%snowIceCover(i,j),emissv,emissh,emissv_std,emissh_std)
           endif
           emissoutL(i,j,1:13)=(/emissv(1),emissh(1),&
                                emissv(2),emissh(2),&
                                emissv(3), &
                                emissv(4),emissh(4),&
                                emissv(5),emissh(5),&
                                emissv(6),emissh(6),&
                                emissv(6),emissv(6)/)
!begin  SJM 2/14/17 emissivity rms vectors
           emis_rms_NS(i,j,1:13)=(/emissv_std(1),emissh_std(1),&
                                emissv_std(2),emissh_std(2),&
                                emissv_std(3), &
                                emissv_std(4),emissh_std(4),&
                                emissv_std(5),emissh_std(5),&
                                emissv_std(6),emissh_std(6),&
                                emissv_std(6),emissv_std(6)/)
           emis_rms_MS(i,j,1:13) = emis_rms_NS(i,j,1:13)
!end   SJM 2/14/17
           if(maxval(emissoutL(i,j,1:13)) .gt. 2.) print'(A6,2I5,13F8.3)','nopcp ', i,j,emissoutL(i,j,1:13)
           !end SJM 10/16/15
           do ik=1,9
              tbNoOcean(i,j,ik)=tbRgrid(ik,i,j+icL)
           enddo
           do ik=1,9
              if(wfract(3,3)<5) then
                 tbNoOcean(i,j,ik)=tbRgrid(ik,i,j+icL)
              endif
           enddo
        endif
     endif
     !if(maxval(emissout(i,j,1:13)) .gt. 2.) print*, i,j, emissout(i,j,1:13)
     !if(maxval(dPRRet%emis(i,j,:,:,:)) .gt. 2.) print*, dPRRet%emis(i,j,:,:,:)
  enddo
enddo

!print*, scLonPR
!print*, scLatPR

print*, maxval(dZms)  !!MG addition Feb 10, 2017
!print*,'iactOb=',iactob
!stop
!print*, maxval(wfmap)

! footprintmap(ifreq,wfmap,fpmap,n,ns,lat,lon,scLon,scLat)
do ifreq1=1,9
   call footprintmap2(ifreqG(ifreq1),wfmap(1:49,1:dPRData%n1c21),&
        fpmap(1:49,1:dPRData%n1c21,ifreq1),&
        49,dPRData%n1c21,dPRData%xlat(1:49,1:dPRData%n1c21), &
        dPRData%xlon(1:49,1:dPRData%n1c21),&
        scLonPR(1:49,1:dPRData%n1c21),scLatPR(1:49,1:dPRData%n1c21))
enddo
fpmapN=1
tb=tbMean
!tbRgrid(1:9,i,j+icL)
nf=9
!call endprofs()
!print*, maxval(tbRgrid(1:9,1:49,ic+1:ic+dPRData%n1c21))
call setoptvars(covTb(1:49,1:dPRData%n1c21,1:9,1:9),&
     invCovTb(1:49,1:dPRData%n1c21,1:9,1:9),tbMean(1:49,1:dPRData%n1c21,1:9),&
     tbMin(1:49,1:dPRData%n1c21,1:9),&
     tbMax(1:49,1:dPRData%n1c21,1:9),&
     tb(1:49,1:dPRData%n1c21,1:9),&
     tbRgrid(1:9,1:49,ic+1:ic+dPRData%n1c21),&
     tbObs(1:49,1:dPRData%n1c21,1:9),wfmap(1:49,1:dPRData%n1c21),&
     49,dPRData%n1c21,nf)

tb0=tb
tbout2d=-99
tbout2dNoOcean=-99
iconv=1
print*, 'before conv'
if(iconv==1) then
   call convallfreq(actOb,tb0(:,:,1:9),tbMean(:,:,1:9),&
        invCovTb(:,:,1:9,1:9),&
        tbObs(:,:,1:9),tbout2D(:,:,1:9),dfdtb(:,:,1:9),49,dPRData%n1c21,&
        dPRData%xlat(:,1:dPRData%n1c21), dPRData%xlon(1:49,1:dPRData%n1c21),&
        scLonPR(1:49,1:dPRData%n1c21),scLatPR(1:49,1:dPRData%n1c21),&
        wfmap(1:49,1:dPRData%n1c21),&
        fpmap(1:49,1:dPRData%n1c21,1:9), &
        nf,fobj,ifreqG(1:9),sfcRain(1:49,1:dPRData%n1c21),ialg)
   
   call convallfreq(actOb,tbNoOcean(:,:,1:9),tbNoOcean(:,:,1:9),&
        invCovTb(:,:,1:9,1:9),&
        tbObs(:,:,1:9),tbout2DNoOcean(:,:,1:9),dfdtb(:,:,1:9),&
        49,dPRData%n1c21,&
        dPRData%xlat(:,1:dPRData%n1c21), dPRData%xlon(1:49,1:dPRData%n1c21),&
        scLonPR(1:49,1:dPRData%n1c21),scLatPR(1:49,1:dPRData%n1c21),&
        wfmap(1:49,1:dPRData%n1c21),&
        fpmapN(1:49,1:dPRData%n1c21,1:9), &
        nf,fobj,ifreqG(1:9),sfcRain(1:49,1:dPRData%n1c21),ialg)
   print*, 'after conv'
endif
!1 sqrt(117/ 17)+  1 sqrt(127/473)= 3.1415927329


131 format(18(F6.2,1x))

tb0=tb

print*, maxval(dprret%tb)

do j=3,dPRData%n1c21-3
   do i=3,49-2
      if(dPRData%raintype(i,j)>0 .and. &
           minval(tbout2dNoOcean(i,j,1:9))>0 .and. &
           minval(tbRgrid(1:9,i,ic+j))>0. .and. wfmap(i,j)<.1) then
         !write(*,131) tbRgrid(1:9,i,j+icL),tbout2dNoOcean(i,j,1:9)
         !!write
      endif
   enddo
enddo

call clearsc()
call asciiplot(sfcRain(:,1:150),49,150,2,1,1e-1,100.)
print*,''
print*,''
call asciiplot(sfcRain(:,151:300),49,150,2,1,1e-1,100.)

if(iconv==1) then
   call updateTbs(dPRData%n1c21,tbObs(:,:,1:9),&
        tbout2D(:,:,1:9),tb(:,:,1:9),fpmap(:,:,1:9),ialg)
   call updateTbsL(dPRData%n1c21,tbObs(:,:,1:9),tbout2DNoOcean(:,:,1:9),&
        tbNoOcean(:,:,1:9),fpmapN(:,:,1:9),ialg)
  
!begin  MG 9/17/18 added following code to adjust deconvolution
   if(ialg==2) then
      do i=1,2
         call updateTbs(dPRData%n1c21,tbObs(:,:,1:9),&
              tbout2D(:,:,1:9),tb(:,:,1:9),fpmap(:,:,1:9),ialg)
         call updateTbsL(dPRData%n1c21,tbObs(:,:,1:9),tbout2DNoOcean(:,:,1:9),&
              tbNoOcean(:,:,1:9),fpmapN(:,:,1:9),ialg)
      enddo
   endif
!end    MG 9/17/18
endif

print*,'ialg=',ialg
print*, maxval(dprret%tb)
print*, maxval(tb)
!stop
tb0=tb
print*, 'after update'

do j=3,-dPRData%n1c21-3
   do i=3,49-2
      if(dPRData%raintype(i,j)>0 .and. minval(tbout2d(i,j,1:9))>0 .and. &
         minval(tbRgrid(1:9,i,ic+j))>0.) then
         call setemtbm(emtbm(1:9),dPRRet%emtb(i,j,:,:,1:nmemb1),nmfreq2,nmemb1)
         call dboux2(tb(i,j,1:9),tb0(i,j,1:9),&
              emtbm,ifreqG,ipolG,fem)

         rms1=1e7
         do jk=1,7
            fem(jk)=fem(jk)/10.
            dPRRet%tb(i,j,ipolG(jk),ifreqG(jk),:)= &
                 fem(jk)*dPRRet%emtb(i,j,ipolG(jk),ifreqG(jk),:)+&
                 (1-fem(jk))*dPRRet%tb(i,j,ipolG(jk),ifreqG(jk),:)
         enddo

         !do ik=1,nmemb1
         !   rms2=0
            
         !   do jk=1,5
         !      rms2=rms2+(tb(i,j,jk)- &
         !           dPRRet%tb(i,j,ipolG(jk),ifreqG(jk),ik))**2
         !   enddo
         !
         !   rmsS(ik)=rms2
         !   if(rms2<rms1) then
         !      rms1=rms2
         !      do jk=1,9
         !         tb0(i,j,jk)= &
         !              dPRRet%tb(i,j,ipolG(jk),ifreqG(jk),ik)
         !      enddo
         !   end if
         !enddo
         !tb0(i,j,1:9)=tb(i,j,1:9)
      else
         if (minval(tbRgrid(1:9,i,ic+j))>0.) then
            tb0(i,j,1:9)=tb(i,j,1:9)
         endif
      endif
   enddo
enddo





!do j=3,dPRData%n1c21-3
!   do i=3,49-2
!      if(dPRData%raintype(i,j)>0 .and. minval(tbout2d(i,j,:))>0 .and. &
!           minval(tbRgrid(:,i,ic+j))>0.) then
!         write(*,131) tbRgrid(1:9,i,j+icL),tbout2d(i,j,1:9)
!      endif
!   enddo
!enddo
!call endtbs()
!call closeascii()
132 format(10(F6.2,1x),I3)
105 format(51(F6.2,1x))

close(10)
close(20)
103 format(13(F8.3,1x))
104 format(50(F7.2,1x))
itcount=itcount+ntbpix
!include 'convolveTbs.f90'
!print*, maxval(dPRRet%convtb)
!  SFM  begin  07/29/2014; for M. Grecu, elminate NANs


do j=1,dPRData%n1c21
   do i=1,49
      !write(34) j,i,dPRdata%xlat(i,j), DPRData%xlon(i,j), tb(i,j,:)
      if(dPRData%rainType(i,j)>=100) then
         do k=max(1,1+dPRData%node(1,i,j)-10),1+dPRData%node(5,i,j)
            if(maxval(dPRRet%rrate(1:nmemb1,k,i,j))>0.01 .and. &
                 minval(dPRRet%rrate(1:nmemb1,k,i,j))<0) then
               do imemb=1,nmemb1
                  if(dPRRet%rrate(imemb,k,i,j)<0) then
                     dPRRet%rrate(imemb,k,i,j)=0.
                     dPRRet%pwc(imemb,k,i,j)=0.
                     dPRRet%d0(imemb,k,i,j)=0.
                     dPRRet%z35mod0(imemb,k,i,j)=0.
                  endif
               enddo
            endif
         enddo
         !print*, minval(dPRRet%rrate(:,dPRData%node(5,i,j):dPRData%node(5,i,j),i,j)),  maxval(dPRRet%rrate(:,dPRData%node(5,i,j):dPRData%node(5,i,j),i,j))
      endif
   enddo
enddo
!close(34)

!  SFM  end    07/29/2014
dPRRet%MS%pwc=dPRRet%pwc
dPRRet%MS%rrate=dPRRET%rrate
dPRRet%MS%d0=dPRRet%d0
dPRRet%MS%log10dNw=dPRRet%log10dNw
dPRRet%MS%sfcRainEns=dPRRet%sfcRainEns
dPRRet%MS%zkuEns=dPRRet%z13c
dPRRet%MS%zkaEns=dPRRet%z35
dPRRet%MS%convtb=dPRRet%convtb
dPRRet%MS%tb=dPRRet%tb
dPRRet%MS%pia13mod=dPRRet%pia13mod
dPRRet%MS%pia35mod=dPRRet%pia35mod
!Start SJM 12/9/2014
dPRRet%MS%sfcWindEns=dPRRet%sfcWindEns
dPRRet%MS%simSigmaZeroKu=dPRRet%simSigmaZeroKu
dPRRet%MS%simSigmaZeroKa=dPRRet%simSigmaZeroKa
dPRRet%MS%emis=dPRRet%emis
!end SJM 12/9/2014

!  SFM  begin  06/22/2014; for M.Grecu  (unknown justification)
print*, maxval(dprret%tb)
print*, 'before Kalman'
tb0MS=tb0
tbNoOceanMS=tbNoOcean
!
!tb=-999
!hfreqPRg=-999
!
print*, ic
print*, maxval(hFreqPRg(:,:,1))
!print*, maxval(dprRet%emis), maxval(dprRet%MS%emis)
do j=1,dPRData%n1c21
   do i=1,49
!  SFM  begin  07/29/2014; for M.Grecu, elminate NANs
      do i1=1,5
         do j1=1,5
            eLon=dPRData%xlon(i,j)
            eLat=dPRData%xlat(i,j)
            call getwfraction(eLat+(i1-2)*0.15,&
                 eLon+(j1-2)*0.15,wfract(i1,j1))
            if(i1.eq.2.and.j1.eq.2) then
               wfmap(i1,j1)=wfract(i1,j1)/100.
            endif
         enddo
      enddo
      wfractm=sum(wfract(:,:))/25.
!      if(i .eq. 25) write(*, '("scan: ", i5, "  lat: ", f10.4, "  lon: ", f10.4, &
!        "  wfractm: ", f10.4)'), j, dPRData%xlat(i,j), dPRData%xlon(i,j), wfractm
      tbRgrid(10:13,i,j+icL)= hFreqPRg(i,j,1:4)
      !stype = LUT%land_class_map(mod(floor((dPRData%xlon(i,j)+180.)/360.*5760.),5760)+1, 2880-floor((dPRData%xlat(i,j)+90.)/180.*2880.)) !SJM 9/9/15
      jj=2880-floor((DPRData%xlat(i,j)+90.)/180.*2880.)
      if(jj .lt. 1) jj=1
      if(jj .gt. 2880) jj = 2880
      ii=floor((DPRData%xlon(i,j)+180.)/360.*5760.)+1
      if(ii .lt. 1) ii = 1
      if(ii .gt. 5760) ii = 5760
      if(dPRData%snowIceCover(i,j) .eq. 0) then
        stype = LUT%land_class_map_bare(ii,jj) !SJM 9/9/15
      else
        stype = LUT%land_class_map_snow(ii,jj) !SJM 9/9/15
      endif
      if( dPRData%rainType(i,j)>0) then ! &
 !          .and. j>0 .and. j<dPRData%n1c21-2 )           &
!      then
         ntbpix2=ntbpix2+1
         !print*, i, j
         if(ifdpr(1:1)=='Y') goto 30
      !tbout(10)=sum(dPRRet%tb(i,j,1,6,1:1*nmemb1))/(nmemb1)
      !tbout(11)=sum(dPRRet%tb(i,j,2,6,1:1*nmemb1))/(nmemb1)
      !tbout(12)=sum(dPRRet%tb(i,j,1,7,1:1*nmemb1))/(nmemb1)
      !tbout(13)=sum(dPRRet%tb(i,j,1,8,1:1*nmemb1))/(nmemb1)
      !do k=10,13
      !   tbout(k)=hFreqPRg(i,j,k-9)
      !enddo
!          write(31) i,j,ichunk,dPRData%sigmaZeroKu(i,j), dPRData%sigmaZeroKa(i,j), sigmaZeroVarKu(i,j),sigmaZeroVarKa(i,j),sigmaZeroCov(i,j), tb(i,j,1:9),hFreqPRg(i,j,:), &!tbRgrid(1:9,i,j+icL), &!scan/ray number, observed sigma_zero
!          dPRRet%sfcRainEns(i,j,1:nmemb1), & !Sfc rain for each ens member
!          dPRRet%sfcd0Ens(i,j,1:nmemb1), & !Sfc D0 for each ens member
!          dPRRet%sfcNwEns(i,j,1:nmemb1), &!Sfc Nw for each ens member
!          dPRRet%sfcWindEns(i,j,1:nmemb1), &!wind for each ens member
!          dPRRet%pia13mod(i,j,1:nmemb1), &!PIAKu for each ens member
!          dPRRet%pia35mod(i,j,1:nmemb1), &!PIAKa for each end member
!          dPRRet%simSigmaZeroKu(i,j,1:nmemb1), &!sigma_zero_Ku for each ens member
!          dPRRet%simSigmaZeroKa(i,j,1:nmemb1), &!sigma_zero_Ka for each ens member
!          dPRRet%tb(i,j,:,:,1:nmemb1), &!Pixel Tb for each ens member (10-19, at least)
!          dPRRet%emis(i,j,:,:,1:nmemb1)
         if(wfractm>99 .and. stype .eq. 1) then
            call filterUpNS(dPRData,dPRRet, Xens,Yens,Yobs,Xup,  &
                 tb,dprRain,sfcRain,nmemb1,ic,i,j,nx,ny,wfractm,&
                 sigmaZeroVarKu, sigmaZeroVarKa, sigmaZeroCov,&
                 hFreqPRg(i,j,:))
            do k=1,9
               tb0(i,j,k)=sum(dPRRet%tb(i,j,ipolG(k),ifreqG(k),1:nmemb1))/nmemb1 
            enddo
            !print*, tb0(i,j,1:9)

         else
            call filterUpNSLand(dPRData,dPRRet, Xens,Yens,Yobs,Xup,  &
                 tbNoOcean,dprRain,sfcRain,nmemb1,ic,i,j,nx,ny,wfractm,&
                 sigmaZeroVarKu, sigmaZeroVarKa, sigmaZeroCov,&
                 hFreqPRg(i,j,:))
            do k=1,9
               tbNoOcean(i,j,k)=&
                    sum(dPRRet%tb(i,j,ipolG(k),ifreqG(k),1:nmemb1))/nmemb1
            enddo
            !print*, tbNoOcean(i,j,1:9)
            
         endif
!                   write(32) i,j,ichunk,dPRData%sigmaZeroKu(i,j), dPRData%sigmaZeroKa(i,j), dPRData%srtPIAKu(i,j), dPRData%dsrtPIAKu(i,j), dPRData%dsrtPIAKa(i,j), tb(i,j,1:9),hFreqPRg(i,j,:), &!tbRgrid(1:9,i,j+icL), &!scan/ray number, observed sigma_zero
!          dPRRet%sfcRainEns(i,j,1:nmemb1), & !Sfc rain for each ens member
!          dPRRet%sfcd0Ens(i,j,1:nmemb1), & !Sfc D0 for each ens member
!          dPRRet%sfcNwEns(i,j,1:nmemb1), &!Sfc Nw for each ens member
!          dPRRet%sfcWindEns(i,j,1:nmemb1), &!wind for each ens member
!          dPRRet%pia13mod(i,j,1:nmemb1), &!PIAKu for each ens member
!          dPRRet%pia35mod(i,j,1:nmemb1), &!PIAKa for each end member
!          dPRRet%simSigmaZeroKu(i,j,1:nmemb1), &!sigma_zero_Ku for each ens member
!          dPRRet%simSigmaZeroKa(i,j,1:nmemb1), &!sigma_zero_Ka for each ens member
!          dPRRet%tb(i,j,:,:,1:nmemb1), &!Pixel Tb for each ens member (10-19, at least)
!          dPRRet%emis(i,j,:,:,1:nmemb1)
         
         if(i>=13 .and. i<=37) then
            stdpia35=0.
            if(j>1 .and. j<dPRData%n1c21 .and. pia35m(i,j)>0.1) then
               nubfc= 1./srtpiaf(pia35m(i-1:i+1,j-1:j+1),9)
            else
               nubfc=1.
            endif
!begin  WSO 2/10/17 assigne 1/nubfc to output variable
            subfootvariability_MS(i, j) = 1./nubfc 
!end    WSO 2/10/17
            !nubfc=1
            !print*, i,j,nubfc
            do k=0,nmemb1-1
              if(dPRRet%simSigmaZeroKa(i,j,k+1) .ne. -99.) dPRRet%simSigmaZeroKa(i,j,k+1) = dPRRet%simSigmaZeroKa(i,j,k+1)-dPRRet%pia35mod(i,j,k+1)/nubfc!+ds0Ka 
              !print '(I3,13F8.3,2F8.2)', k, radarRet%tb(k*2*nmfreq2+1:k*2*nmfreq2+13), radarRet%simSigmaZeroKu(k+1), radarRet%simSigmaZeroKa(k+1)
            end do
            !nubfc=0.1-0.45*stdpia35+0.42*stdpia35*stdpia35
            !print*, dPRData%dsrtPIAka(i,j),nubfc
            if (nubfc*dPRData%dsrtPIAka(i,j)<50 .and. nubfc>=1.) then
               !dPRData%dsrtPIAka(i,j)=dPRData%dsrtPIAka(i,j)*1.1*nubfc
            endif
            call filterUpMS(dPRData,dPRRet, Xens,Yens,Yobs,Xup,  &
                 tb,dprRain,sfcRain,nmemb1,ic,i,j,nx,ny,wfractm,&
                 sigmaZeroVarKu, sigmaZeroVarKa, sigmaZeroCov,&
                 hFreqPRg(i,j,:),nubfc)
!                         write(33) i,j,ichunk,dPRData%sigmaZeroKu(i,j), dPRData%sigmaZeroKa(i,j), dPRData%srtPIAKu(i,j), dPRData%dsrtPIAKu(i,j), dPRData%dsrtPIAKa(i,j), tb(i,j,1:9), hFreqPRg(i,j,:),&!scan/ray number, observed sigma_zero
!             dPRRet%MS%sfcRainEns(i,j,1:nmemb1), & !Sfc rain for each ens member
!             dPRRet%MS%sfcd0Ens(i,j,1:nmemb1), & !Sfc D0 for each ens member
!             dPRRet%MS%sfcNwEns(i,j,1:nmemb1), &!Sfc Nw for each ens member
!             dPRRet%MS%sfcWindEns(i,j,1:nmemb1), &!wind for each ens member
!             dPRRet%MS%pia13mod(i,j,1:nmemb1), &!PIAKu for each ens member
!             dPRRet%MS%pia35mod(i,j,1:nmemb1), &!PIAKa for each end member
!             dPRRet%MS%simSigmaZeroKu(i,j,1:nmemb1), &!sigma_zero_Ku for each ens member
!             dPRRet%MS%simSigmaZeroKa(i,j,1:nmemb1), &!sigma_zero_Ka for each ens member
!             dPRRet%MS%tb(i,j,:,:,1:nmemb1), &!Pixel Tb for each ens member (10-19, at least)
!             dPRRet%MS%emis(i,j,:,:,1:nmemb1)
            do k=1,9
               tb0MS(i,j,k)=sum(dPRRet%MS%tb(i,j,ipolG(k),&
                    ifreqG(k),1:nmemb1))/nmemb1
               tbNoOceanMS(i,j,k)=sum(dPRRet%MS%tb(i,j,ipolG(k),&
                    ifreqG(k),1:nmemb1))/nmemb1
            enddo
         endif
         !  SFM  end    07/29/2014
30       continue
         !  SFM  begin  07/29/2014; for M.Gecu, elminate NANs
      else
         !if(wfractm>80) &
         !     call filterUpNSClearSky(dPRData,dPRRet, Xens,Yens,Yobs,Xup,  &
         !     tbRgrid,dprRain,sfcRain,nmemb1,ic,i,j,nx,ny)
!  SFM  end    07/29/2014
!simulate clear-sky sigma_zero
        
!         if(stype .eq. 1) then
!           call calc_relAz(DPRData%sclon(j), DPRData%sclat(j), DPRData%xlon(i,j), DPRData%xlat(i,j), DPRData%envSfcWindU(i,j), DPRData%envSfcWIndV(i,j), relAz)
!           call intplte_water_sigma0(i,DPRData%envSfcWind(i,j),relAz,s0Ku, s0Ka, s0stdKu, s0stdKa, s0corr)
!         else
!           call get_s0(DPRData%xlat(i,j), DPRData%xlon(i,j),i,stype,s0Ku,s0Ka, s0stdKu, s0stdKa)
!         endif
!         !print*, i,j,s0Ku,s0Ka
!         dPRRet%simSigmaZeroKu(i,j,1:nmemb1) = s0Ku
!         dPRRet%simSigmaZeroKa(i,j,1:nmemb1) = s0Ka
!         dPRRet%sfcWindEns(i,j,1:nmemb1) = DPRData%envSfcWind(i,j)
!         write(31) i,j,ichunk,dPRData%sigmaZeroKu(i,j), dPRData%sigmaZeroKa(i,j), s0stdKu**2,s0stdKa**2,s0corr*s0stdKu*s0stdKa, tb(i,j,1:9),hFreqPRg(i,j,:), &!tbRgrid(1:9,i,j+icL), &!scan/ray number, observed sigma_zero
!          dPRRet%sfcRainEns(i,j,1:nmemb1), & !Sfc rain for each ens member
!          dPRRet%sfcd0Ens(i,j,1:nmemb1), & !Sfc D0 for each ens member
!          dPRRet%sfcNwEns(i,j,1:nmemb1), &!Sfc Nw for each ens member
!          dPRRet%sfcWindEns(i,j,1:nmemb1), &!wind for each ens member
!          dPRRet%pia13mod(i,j,1:nmemb1), &!PIAKu for each ens member
!          dPRRet%pia35mod(i,j,1:nmemb1), &!PIAKa for each end member
!          dPRRet%simSigmaZeroKu(i,j,1:nmemb1), &!sigma_zero_Ku for each ens member
!          dPRRet%simSigmaZeroKa(i,j,1:nmemb1), &!sigma_zero_Ka for each ens member
!          dPRRet%tb(i,j,:,:,1:nmemb1), &!Pixel Tb for each ens member (10-19, at least)
!          dPRRet%emis(i,j,:,:,1:nmemb1)
      end if
   enddo
enddo

!call closeenkffile2()
! close(31)
! close(32)
! close(33)
print*, 'after filtering'
print*, maxval(dprret%tb)
!print*, maxval(dprRet%emis), maxval(dprRet%MS%emis)
!stop
!print*, maxval(dPRRet%convtb)
31 format(53(F7.2,1x))
!  SFM  end    06/22/2014
do j=1,dPRData%n1c21
   do i=1,49
      if(dPRData%rainType(i,j)>=100) then
!diagnostic
!       write(*, '("early scan: ", i5, "  scene: ", i5, "  nodes: ", 5i5)') j, i, dPRData%node(1:5, i, j)
!end diagnostic
         do k=1+dPRData%node(1,i,j),1+dPRData%node(5,i,j)
!            if(k < 1 .or. k > 88) write(*, '("bad k at: ", 3i10, "  min k: ", i10, "  max k: ", i10)') &
!             i, j, k, 1+dPRData%node(1,i,j), 1+dPRData%node(5,i,j)
            rrate3D(k,i,j)=sum(dPRRet%rrate(1:nmemb1,k,i,j))/nmemb1
            rrate3Dstd(k,i,j)=sqrt(sum((dPRRet%rrate(1:nmemb1,k,i,j)-          &
                 rrate3D(k,i,j))**2)/(nmemb1-1))
            pwc3D(k,i,j)=sum(dPRRet%pwc(1:nmemb1,k,i,j))/nmemb1
            d03D(k,i,j)=sum(dPRRet%d0(1:nmemb1,k,i,j))/nmemb1
            pwc3Dstd(k,i,j)=sqrt(sum((dPRRet%pwc(1:nmemb1,k,i,j)-              &
                 pwc3D(k,i,j))**2)/(nmemb1-1))
            zcKu3D(k,i,j)=sum(dPRRet%z13c(1:nmemb1,k,i,j))/nmemb1
!begin  WSO 8/28/14 if rrate3D is -99 then
            if(abs(rrate3D(k, i, j) - (-99)) < 1.) then
              rrate3D(k, i, j) = 0.
              rrate3Dstd(k, i, j) =  0.
              pwc3D(k, i, j) = 0.
              pwc3Dstd(k, i, j) = 0.
            endif
!end    WSO 8/28/14
         enddo
!begin  WSO 9/15/13 set to missing instead of -99
!diagnostic
!       write(*, '("early scan: ", i5, "  scene: ", i5, "  nbin: ", i5)') nbin
!end diagnostic
!         if(2+dPRData%node(5,i,j) < 1 .or. 2+dPRData%node(5,i,j) > 88) &
!          write(*, '("bad node 5+2at: ", 2i10, "  node5+2: ", i10, "  nbin: ", i10, "  nodes: ", 5i10)') &
!             i, j, 2+dPRData%node(5,i,j), nbin, (dPRData%node(l, i, j), l=1,5)
         rrate3D(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
         rrate3Dstd(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
         pwc3D(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
         pwc3Dstd(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
         d03D(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
!end    WSO 9/15/13
         pia13m=sum(dPRRet%pia13mod(i,j,1:nmemb1))/nmemb1
         pia13s=(sum((dPRRet%pia13mod(i,j,1:nmemb1)-pia13m)**2)/nmemb1)**.5
         sfcRain(i,j)=sum(dPRRet%sfcRainEns(i,j,1:nmemb1))/nmemb1
         sfcRainStd(i,j)=sqrt(sum((dPRRet%sfcRainEns(i,j,1:nmemb1)-            &
              sfcRain(i,j))**2)/(nmemb1-1))
!begin SJM 12/4/2014
         w10_out_NS(i,j) = sum(dPRRet%sfcWindEns(i,j,1:nmemb1))/nmemb1
         w10_rms_NS(i,j) = sqrt(sum((dPRRet%sfcWindEns(i,j,1:nmemb1)-w10_out_NS(i,j))**2)/nmemb1)
         !print*, i,j, w10_out_NS(i,j), w10_rms_NS(i,j)
         !w10_rms_NS(i,j) = sqrt(sum((dPRRet%sfcWindEns(i,j,1:nmemb1)-)/nmemb1
         !print '(2I5,2F8.2)', i,j,dPRData%envSfcWind(i,j), w10_out_NS(i,j)
!end SJM 12/4/2014

!begin  WSO 8/28/14 if sfcRain is -99, set sfcRain and its uncertainty to zero
         if(abs(sfcRain(i, j) - (-99)) < 1.) then
           sfcRain(i, j) = 0.
           sfcRainStd(i, j) = 0.
         endif
!end    WSO 8/28/24 
         piaOut(i,j)=pia13m

!begin  WSO 8/7/13 add liquid fractions in transition layer
         gatelength = 250.
!begin  WSO 12/30/13 set mu_mean = 2 regardless of values
!since all tables represent only mu = 2
!         mu_mean(i, j) = sum(mu_tab(radarRet%imu(1:nmemb1)))/nmemb1
         mu_mean(i, j) = 2.
!end    WSO 12/30/13

         if(dPRData%node(5, i, j) > dPRData%node(4,i,j)) then   !rain at surface

           sfcRainLiqFrac(i, j) = 1.0

         else if(dPRData%node(5, i, j) < dPRData%node(2,i,j)) then !snow at surface

           sfcRainLiqFrac(i, j) = 0.

         endif
         do  k=1+dPRData%node(2,i,j),1+dPRData%node(4,i,j)
!begin  WSO 9/5/13 add logic to account for lowest bin in mixed-phase layer
           if(k .le. 1 + dPRData%node(5, i, j)) then  !valid mixed-phase bin
!end    WSO 9/5/13
             depthBB = (dPRData%node(3,i,j) - dPRData%node(2,i,j)) * gatelength
             depthML = (dPRData%node(4,i,j) - dPRData%node(2,i,j)) * gatelength
             depth = (k - 1 - dPRData%node(2,i,j)) * gatelength

             if(dPRData%rainType(i,j)<200) then  !stratiform
               call interp_melt_percentages(depthBB, depthML, &
                mu_mean(i, j), d03D(k, i, j), depth, mlwc_frac(k - dPRData%node(2,i,j), i, j), &
                mrate_frac(k - dPRData%node(2,i,j), i, j))
             else  !convective or undefined
                mlwc_frac(k - dPRData%node(2,i,j), i, j) = depth / depthML
                mrate_frac(k - dPRData%node(2,i,j), i, j) = depth / depthML
             endif

             if(k == 1 + dPRData%node(5,i,j)) then
               sfcRainLiqFrac(i, j) = mrate_frac(k - dPRData%node(2,i,j), i, j)

               go to 200

             endif

           endif

         enddo
  200    continue

!end  WSO 8/7/13
      else

!begin  WSO 9/22/13 set to missing in clutter zone
!diagnostic
!       write(*, '("later scan: ", i5, "  scene: ", i5, "  nodes: ", 5i5, "  nbin: ", i5)') &
!        j, i, dPRData%node(1:5, i, j), nbin
!end diagnostic

         !Change from MG, 5/2/18
         !if(dPRData%badRayFlag(i,j) == 0)
         if(dPRData%badRayFlag(i,j) == 0 .and. &
             2+dPRData%node(5,i,j)>1) & 
         then   !valid non-raining point

!           if(2+dPRData%node(5,i,j) < 1 .or. 2+dPRData%node(5,i,j) > 88) &
!             write(*, '("bad node 5+2at: ", 2i10, "  node5+2: ", i10, "  nbin: ", i10, " nodes: ", 5i10, "  lat/lon: ", 2f10.4)') &
!             i, j, 2+dPRData%node(5,i,j), nbin, (dPRData%node(l, i, j), l=1,5), dPRData%xlat(i,j), dPRData%xlon(i,j)
           rrate3D(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
           rrate3Dstd(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
           pwc3D(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
           pwc3Dstd(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
           d03D(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
         else                               !invalid nodes
           sfcRain(i,j)=missing_r4
           sfcRainStd(i,j)=missing_r4
           rrate3D(1:nbin,i,j) = missing_r4
           rrate3Dstd(1:nbin,i,j) = missing_r4
           pwc3D(1:nbin,i,j) = missing_r4
           pwc3Dstd(1:nbin,i,j) = missing_r4
           d03D(1:nbin,i,j) = missing_r4
         endif

!end    WSO 9/22/13

      endif
   end do
end do
print*, 'before mlw'
do j=1,dPRData%n1c21
   do i=13,37
      if(dPRData%rainType(i,j)>=100) then
         do k=1+dPRData%node(1,i,j),1+dPRData%node(5,i,j)
            rrate3DMS(k,i,j)=sum(dPRRet%MS%rrate(1:nmemb1,k,i,j))/nmemb1
            rrate3DstdMS(k,i,j)=sqrt(sum((dPRRet%MS%rrate(1:nmemb1,k,i,j)-     &
                 rrate3DMS(k,i,j))**2)/(nmemb1-1))
            pwc3DMS(k,i,j)=sum(dPRRet%MS%pwc(1:nmemb1,k,i,j))/nmemb1
            d03DMS(k,i,j)=sum(dPRRet%MS%d0(1:nmemb1,k,i,j))/nmemb1
            pwc3DstdMS(k,i,j)=sqrt(sum((dPRRet%pwc(1:nmemb1,k,i,j)-            &
                 pwc3D(k,i,j))**2)/(nmemb1-1))
            zcKu3DMS(k,i,j)=sum(dPRRet%MS%zkuEns(1:nmemb1,k,i,j))/nmemb1
            zcKa3DMS(k,i,j)=sum(dPRRet%MS%zkaEns(1:nmemb1,k,i,j))/nmemb1
            !zcKa3DMS(k,i,j)=sum(dPRRet%z35mod0(1:nmemb1,k,i,j))/nmemb1
!begin  WSO 8/28/14 if rrate3DMS is -99 then
            if(abs(rrate3DMS(k, i, j) - (-99)) < 1.) then
              rrate3DMS(k, i, j) = 0.
              rrate3DstdMS(k, i, j) =  0.
              pwc3DMS(k, i, j) = 0.
              pwc3DstdMS(k, i, j) = 0.
            endif
!end    WSO 8/28/14
         enddo
!begin  WSO 9/15/13 set to missing instead of -99
         rrate3DMS(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
         rrate3DstdMS(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
         pwc3DMS(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
         pwc3DstdMS(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
         d03DMS(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
!end    WSO 9/15/13
         sfcRainMS(i,j)=sum(dPRRet%MS%sfcRainEns(i,j,1:nmemb1))/nmemb1
         sfcRainStdMS(i,j)=sqrt(sum((dPRRet%MS%sfcRainEns(i,j,1:nmemb1)-       &
           sfcRainMS(i,j))**2)/(nmemb1-1))
!begin  WSO 8/28/14 if sfcRainMS is -99, set to sfcRainMS and its uncertainty to zero
         if(abs(sfcRainMS(i, j) - (-99)) < 1.) then
           sfcRainMS(i, j) = 0.
           sfcRainStdMS(i, j) = 0.
         endif
!end    WSO 8/28/14
         piaOutKuMS(i,j)=sum(dPRRet%MS%pia13mod(i,j,1:nmemb1))/nmemb1
         piaOutKaMS(i,j)=sum(dPRRet%MS%pia35mod(i,j,1:nmemb1))/nmemb1
         w10_out_MS(i,j) = sum(dPRRet%MS%sfcWindEns(i,j,1:nmemb1))/nmemb1 !SJM 2/5/15
         w10_rms_MS(i,j) = sqrt(sum((dPRRet%MS%sfcWindEns(i,j,1:nmemb1)-w10_out_MS(i,j))**2)/nmemb1)
         !print*, i,j, w10_out_MS(i,j), w10_rms_MS(i,j)
!begin  WSO 8/7/13 add liquid fractions in transition layer
         gatelength = 250.
!begin  WSO 12/30/13 set mu_meanMS = 2, since
!all tables only represent mu = 2
!         mu_meanMS(i, j) = sum(mu_tab(radarRet%imu(1:nmemb1)))/nmemb1
         mu_meanMS(i, j) = 2.
!end    WSO 12/30/13

         if(dPRData%node(5, i, j) > dPRData%node(4,i,j)) then   !rain at surface

           sfcRainLiqFracMS(i, j) = 1.0


         else if(dPRData%node(5, i, j) < dPRData%node(2,i,j)) then !snow at surface

           sfcRainLiqFracMS(i, j) = 0.

         endif

!diagnostic
!           if(dPRData%node(5,i,j) .le. dPRData%node(4, i,j) .and. &
!            dPRData%node(5,i,j) .ge. dPRData%node(2, i,j)) then
!             write(6, '(/,"  Level  i: ", i5, "  j: ", i5, "  lat: ", f10.4, " lon: ", f10.4, &
!             "  dPRnodes: ", 5i5)') &
!                 i, j, dPRData%xlat(i, j), dPRData%xlon(i, j), dPRData%node(1:5, i, j) + 1
!           endif
!end diagnostic

         do  k=1+dPRData%node(2,i,j),1+dPRData%node(4,i,j)
!begin  WSO 9/5/13 add logic to account for lowest bin in mixed-phase layer
           if(k .le. 1 + dPRData%node(5, i, j)) then  !valid mixed-phase bin
!end    WSO 9/5/13

             depthBB = (dPRData%node(3,i,j) - dPRData%node(2,i,j)) * gatelength
             depthML = (dPRData%node(4,i,j) - dPRData%node(2,i,j)) * gatelength
             depth = (k - 1 - dPRData%node(2,i,j)) * gatelength

             if(dPRData%rainType(i,j)<200) then  !stratiform
               call interp_melt_percentages(depthBB, depthML, &
                mu_meanMS(i, j), d03DMS(k, i, j), depth, mlwc_fracMS(k - dPRData%node(2,i,j), i, j), &
                mrate_fracMS(k - dPRData%node(2,i,j), i, j))
             else  !convective or undefined
                mlwc_fracMS(k - dPRData%node(2,i,j), i, j) = depth / depthML
                mrate_fracMS(k - dPRData%node(2,i,j), i, j) = depth / depthML
             endif


             if(k - 1 == dPRData%node(5,i,j)) then
               sfcRainLiqFracMS(i, j) = mrate_fracMS(k - dPRData%node(2,i,j), i, j)

               go to 300

             endif

           endif

         enddo
  300    continue

!end  WSO 8/7/13
      else

!begin  WSO 9/22/13 set to missing instead of -99

        !Change from MG, 5/2/18
        !if(dPRData%badRayFlag(i,j) == 0) 
         if(dPRData%badRayFlag(i,j) == 0 .and. &
             2+dPRData%node(5,i,j)>1) &
        then    !valid non-raining footprint


          rrate3DMS(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
          rrate3DstdMS(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
          pwc3DMS(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
          pwc3DstdMS(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
          d03DMS(2+dPRData%node(5,i,j):nbin,i,j) = missing_r4
        else
          sfcRainMS(i,j)=missing_r4
          sfcRainStdMS(i,j)=missing_r4
          rrate3DMS(1:nbin,i,j) = missing_r4
          rrate3DstdMS(1:nbin,i,j) = missing_r4
          pwc3DMS(1:nbin,i,j) = missing_r4
          pwc3DstdMS(1:nbin,i,j) = missing_r4
          d03DMS(1:nbin,i,j) = missing_r4
        endif
!end    WSO 9/22/13

      endif
   enddo
enddo
print*, 'before out'
print*, maxval(tbout2D), maxval(tbout2dnoocean)
40 format(14(F7.2,1x))


if(iconv==1) then
   call convallfreq(actOb,tb0(:,:,1:9),tbMean(:,:,1:9),&
     invCovTb(:,:,1:9,1:9),&
     tbObs(:,:,1:9),tbout2D(:,:,1:9),dfdtb(:,:,1:9),49,dPRData%n1c21,&
     dPRData%xlat(:,1:dPRData%n1c21), dPRData%xlon(1:49,1:dPRData%n1c21),&
     scLonPR(1:49,1:dPRData%n1c21),scLatPR(1:49,1:dPRData%n1c21),&
     wfmap(1:49,1:dPRData%n1c21),&
     fpmap(1:49,1:dPRData%n1c21,1:9), &
     nf,fobj,ifreqG,sfcRain(1:49,1:dPRData%n1c21),ialg)

   call convallfreq(actOb,tbNoOcean(:,:,1:9),tbNoOcean(:,:,1:9),&
        invCovTb(:,:,1:9,1:9),&
        tbObs(:,:,1:9),tbout2DNoOcean(:,:,1:9),dfdtb(:,:,1:9),&
        49,dPRData%n1c21,&
        dPRData%xlat(:,1:dPRData%n1c21), dPRData%xlon(1:49,1:dPRData%n1c21),&
        scLonPR(1:49,1:dPRData%n1c21),scLatPR(1:49,1:dPRData%n1c21),&
        wfmap(1:49,1:dPRData%n1c21),&
        fpmapN(1:49,1:dPRData%n1c21,1:9), &
        nf,fobj,ifreqG(1:9),sfcRain(1:49,1:dPRData%n1c21),ialg)
   
   call convallfreq(actOb,tb0MS(:,:,1:9),tbMean(:,:,1:9),&
        invCovTb(:,:,1:9,1:9),&
        tbObs(:,:,1:9),tbout2DMS(:,:,1:9),dfdtb(:,:,1:9),49,dPRData%n1c21,&
        dPRData%xlat(:,1:dPRData%n1c21), dPRData%xlon(1:49,1:dPRData%n1c21),&
        scLonPR(1:49,1:dPRData%n1c21),scLatPR(1:49,1:dPRData%n1c21),&
        wfmap(1:49,1:dPRData%n1c21),&
        fpmap(1:49,1:dPRData%n1c21,1:9), &
        nf,fobj,ifreqG,sfcRain(1:49,1:dPRData%n1c21),ialg)
   
   call convallfreq(actOb,tbNoOceanMS(:,:,1:9),tbNoOceanMS(:,:,1:9),&
        invCovTb(:,:,1:9,1:9),&
        tbObs(:,:,1:9),tbout2DNoOceanMS(:,:,1:9),dfdtb(:,:,1:9),&
        49,dPRData%n1c21,&
        dPRData%xlat(:,1:dPRData%n1c21), dPRData%xlon(1:49,1:dPRData%n1c21),&
        scLonPR(1:49,1:dPRData%n1c21),scLatPR(1:49,1:dPRData%n1c21),&
        wfmap(1:49,1:dPRData%n1c21),&
        fpmapN(1:49,1:dPRData%n1c21,1:9), &
        nf,fobj,ifreqG(1:9),sfcRain(1:49,1:dPRData%n1c21),ialg)
endif

print*, 'ialg=',ialg
!  SFM  begin  12/13/2013; conditional to rewind dpr
call rewind(ic)
IF (st_2adpr .EQ. 0) call rewind_dpr(ic)

 
do i=1,25
   do j=1,dPRData%n1c21
      emiss2d(j,i,1) = sum(dprRet%MS%emis(i+12,j,1,1,1:nmemb1))/nmemb1
      emiss2d(j,i,2) = sum(dprRet%MS%emis(i+12,j,2,1,1:nmemb1))/nmemb1
      emiss2d(j,i,3) = sum(dprRet%MS%emis(i+12,j,1,2,1:nmemb1))/nmemb1
      emiss2d(j,i,4) = sum(dprRet%MS%emis(i+12,j,2,2,1:nmemb1))/nmemb1
      emiss2d(j,i,5) = sum(dprRet%MS%emis(i+12,j,1,3,1:nmemb1))/nmemb1
      emiss2d(j,i,6) = sum(dprRet%MS%emis(i+12,j,1,4,1:nmemb1))/nmemb1
      emiss2d(j,i,7) = sum(dprRet%MS%emis(i+12,j,2,4,1:nmemb1))/nmemb1
      emiss2d(j,i,8) = sum(dprRet%MS%emis(i+12,j,1,5,1:nmemb1))/nmemb1
      emiss2d(j,i,9) = sum(dprRet%MS%emis(i+12,j,2,5,1:nmemb1))/nmemb1
      emiss2d(j,i,10) = sum(dprRet%MS%emis(i+12,j,1,6,1:nmemb1))/nmemb1
      emiss2d(j,i,11) = sum(dprRet%MS%emis(i+12,j,2,6,1:nmemb1))/nmemb1
      pRate(j,i,:)=rrate3DMS(:,i+12,j)
      swc3D(j,i,:)=pwc3DMS(:,i+12,j)
      airTemp3D(j,i,:)=dPRData%envTemp(:,i+12,j)
      sfcTemp(j,i)=dPRData%envSfcTemp(i+12,j)
      press3D(j,i,:)=dPRData%envPress(:,i+12,j)
      qv3D(j,i,:)=dPRData%envQv(:,i+12,j)
      z13(j,i,:)=zcKu3DMS(:,i+12,j)
      do k=1,88
         if(dPRData%rainType(i+12,j)>=100) then
            nw3d(j,i,k)=sum(dPRRet%MS%log10dNw(1:nmemb1,k,i+12,j))/nmemb1
         else
            nw3d(j,i,k)=-0.5
         endif
      enddo
      sfcBin(j,i)= dprData%binRealSurface(i+12,j)
      cldw3d(j,i,:)=0.
      if(dPRData%rainType(i+12,j)>=100) then
         call getcldwfromcoeff(dPRRet%cldwcoeff(i+12,j,:,:),                      &
              dPRData%freezH(i+12,j)/1000.,cldwprof, dPRData%node(:,i+12,j),         &
              radarData%dr,dPRData%localZenithAngle(i+12,j),nmemb1)
         do k = 1, 88
            if(cldwprof(k) < 0.) then
               cldwprof(k) = 0.
            endif
         enddo
         cldw3d(j,i,:)=cldwprof
      endif
      tbobsT(j,i,1:9)=tbobs(i+12,j,1:9)
      tbobsT(j,i,10:13)=hFreqPRg(i+12,j,1:4)
      pType(j,i)=dPRData%rainType(i+12,j)
      binNodes(j,i,:)=dPRData%node(:,i+12,j)
      clutFree(j,i)=binNodes(j,i,5)
      do k = 1, 10
         envNode(j,i,k)= 88 - nint((env_levs(k)/&
              cos(dPRData%localZenithAngle(i+12, j) * pi / 180.)) * 4.)-1
      enddo
   enddo
enddo
nfreq=7
!idir=1
call frteprep(binNodes,pRate,swc3d,pRateOut,swcOut,nwOut,z13,emiss2d,&
     envNode,pType,&
     qv3D,press3D,airTemp3D,nw3d,tbsim,nscans,npixs,nlev,nchans,idir,&
     sfcBin,sfcTemp,cldw3d,tbobsT,nfreq,clutFree,dPRData%n1c21)

!pRateOut=pRate
do i=1,-25
   do j=1,dPRData%n1c21
      if(tbSim(j,i,10)>10) then
         print*, tbSim(j,i,10),tbobsT(j,i,10),&
              sum(dPRRet%tb(i+12,j,1,6,1:1*nmemb1))/(nmemb1)
      endif
   enddo
enddo
do i=1,25
   do j=1,dPRData%n1c21
      do k=1,min(binNodes(j,i,5),binNodes(j,i,4))+1
         rrate3DMS(k,i+12,j)=pRateOut(j,i,k)
         pwc3DMS(k,i+12,j)=swcOut(j,i,k)
      enddo
      do k=min(binNodes(j,i,5),binNodes(j,i,2))+1,&
           min(binNodes(j,i,5),binNodes(j,i,3))+1
         !xf=1-(k-binNodes(j,i,2)-1.)/(binNodes(j,i,3)-binNodes(j,i,2)+1e-5)
         !rrate3DMS(k,i+12,j)=xf*pRateOut(j,i,k)+(1-xf)*pRate(j,i,k)
         !pwc3DMS(k,i+12,j)=xf*swcOut(j,i,k)+(1-xf)*swc3D(j,i,k)
      enddo
      do k=1,88
         if(dPRData%rainType(i+12,j)>=100) then
            dPRRet%MS%log10dNw(1:nmemb1,k,i+12,j)=&
                 dPRRet%MS%log10dNw(1:nmemb1,k,i+12,j)-nw3d(j,i,k)+nwOut(j,i,k)
         endif
      enddo
   enddo
enddo
!stop
do j=1,dPRData%n1c21
   if(ialg==1) then
      do i=1,25
         call copyw10small(w10_out_MS(i+12,j),i-1) !SJM 12/9/2014
         call copyw10smallsigma(w10_rms_MS(i+12,j),i-1) !SJM 12/9/2014
      enddo
   endif
   if(ialg==2) then
      call setlatlons1t( dPRData%xlat(:,j), dPRData%xlon(:,j),                &
           sfcRain(:,j),sfcRainStd(:,j),piaOut(:,j))
   else
      call setlatlons1( dPRData%xlat(:,j), dPRData%xlon(:,j),                &
           sfcRain(:,j),sfcRainStd(:,j),piaOut(:,j))
      call setlatlons2( dPRData%xlat(:,j), dPRData%xlon(:,j),                 &
           sfcRainMS(:,j),sfcRainStdMS(:,j),piaOutKuMS(:,j),piaOutKaMS(:,j))
   endif


   do i=1,49
      if(ialg==1) then
         !begin  WSO 9/28/13 added rain flag that includes bad scan as missing
         call copyrainflags1(dPRData%rainFlagBad(i,j), i-1)
         !end    WSO 9/28/13
         !begin  WSO 8/15/2014 added ioquality flag 
         call copyioqualitys1(dPRData%ioqualityflagku(i,j), i-1)
         !end    WSO 8/15/2014
         !begin  WSO 3/17/17 add snow ice cover flag to output
         call copysnowices1(dPRData%snowIceCover(i,j), i-1)
         !end    WSO 3/17/17
         
         !begin  WSO 2/8/17 new output variables
         call copyinitnws1(initnw_NS(:, i, j), dPRRet%n9(:, i, j), i-1)
         call copyprincomps1(princomp_NS(:, i, j), i-1)
         call copyprofclasss1(profclass_NS(i, j), i-1)
         call copysurfprecipbiasratios1(surfprecipbiasratio_NS(i, j), i-1)
         call copysubfootvariabilitys1(subfootvariability_NS(i, j), i-1)
         call copymultiscatcalcs1(multiscatcalc_NS(i, j), i-1)
         call copymultiscatsurfaces1(multiscatsurface_NS(i, j), i-1)
         !end    WSO 2/8/17
      else
         call copyrainflags1t(dPRData%rainFlagBad(i,j), i-1)
         call copyioqualitys1t(dPRData%ioqualityflagku(i,j), i-1)
         call copysnowices1t(dPRData%snowIceCover(i,j), i-1)
         call copyinitnws1t(initnw_NS(:, i, j), dPRRet%n9(:, i, j), i-1)
         call copyprincomps1t(princomp_NS(:, i, j), i-1)
         call copyprofclasss1t(profclass_NS(i, j), i-1)
         call copysurfprecipbiasratios1t(surfprecipbiasratio_NS(i, j), i-1)
         call copysubfootvariabilitys1t(subfootvariability_NS(i, j), i-1)
         call copymultiscatcalcs1t(multiscatcalc_NS(i, j), i-1)
         call copymultiscatsurfaces1t(multiscatsurface_NS(i, j), i-1)
      endif
      !begin  WSO 9/19/13 more accurate environmental bin calculation
      do k = 1, 10
         env_nodes(k, i) = 88 - nint((env_levs(k)/cos(dPRData%localZenithAngle(i, j) * pi / 180.)) * 4.)
         !       write(*, '("i: ", i5, "  k: ", i5, "  zenangl: ", f8.4, "  env_node: ", &
         !        i5)') i, k, dPRData%localZenithAngle(i, j), env_nodes(k, i)
      enddo
      !end    WSO 9/19/13
      
      if(ialg==2) then
         call copysigmapias1t(dPRData%srtsigmaPIAku(i, j), i-1)
         call copysigmazeros1t(dPRData%sigmaZeroKu(i, j), i-1)
      else
         call copysigmapias1(dPRData%srtsigmaPIAku(i, j), i-1)
         call copysigmazeros1(dPRData%sigmaZeroKu(i, j), i-1)
      endif

      if(dPRData%rainType(i,j)>=100) then
         !print*,minval(dPRRet%cldwcoeff(i,j,:,:))  !MG 09/26/13
         call getcldwfromcoeff(dPRRet%cldwcoeff(i,j,:,:),                      &
              dPRData%freezH(i,j)/1000.,cldwprof, dPRData%node(:,i,j),         &
              radarData%dr,dPRData%localZenithAngle(i,j),nmemb1)
!begin  WSO 9/17/13 remove negative cloud water contents and standardize missing
!                  values
        do k = 1, 88
          if(cldwprof(k) < 0. .and. cldwprof(k) > -90.) then
            cldwprof(k) = 0.
          else if(cldwprof(k) < -90.) then
            cldwprof(k) = missing_r4
          endif
        end do
!end    WSO 9/17/13
      else
         cldwprof=0
      endif
      if(ialg==2) then
         call copycldwaters1t(cldwprof,i-1)
         !begin  WSO 9/15/13 add missing for cloud ice profiles
         cldiprof = missing_r4
         call copycldices1t(cldiprof,i-1)
         !end  WSO 9/15/13 
         call copyrrates1t(rrate3D(:,i,j),rrate3Dstd(:,i,j),i-1)
         call copypwcs1t(pwc3D(:,i,j),pwc3Dstd(:,i,j),i-1)
         !begin  WSO 8/7/13
         call copylwcfracs1t(mlwc_frac(:,i,j),mrate_frac(:,i,j),i-1)
         call copysfcrainliqfracs1t(sfcRainLiqFrac(i, j), i-1)
         !end    WSO 8/7/13
         call copyd0s1t(d03D(:,i,j),i-1)
         call copyzckus1t(zcKu3D(:,i,j),i-1)
         call copynodess1t(dPRData%node(:,i,j),i-1)
      else
         call copycldwaters1(cldwprof,i-1)
         cldiprof = missing_r4
         call copycldices1(cldiprof,i-1)
         call copyrrates1(rrate3D(:,i,j),rrate3Dstd(:,i,j),i-1)
         call copypwcs1(pwc3D(:,i,j),pwc3Dstd(:,i,j),i-1)
         call copylwcfracs1(mlwc_frac(:,i,j),mrate_frac(:,i,j),i-1)
         call copysfcrainliqfracs1(sfcRainLiqFrac(i, j), i-1)
         call copyd0s1(d03D(:,i,j),i-1)
         call copyzckus1(zcKu3D(:,i,j),i-1)
         call copynodess1(dPRData%node(:,i,j),i-1)
      endif
         
!begin  WSO 8/19/13 change dNw to Nw and add mu
      if(dPRData%rainType(i,j)>=100) then
         do k=1,88
            log10NwMean(k)=sum(dPRRet%log10dNw(1:nmemb1,k,i,j))/nmemb1 + &
             log10(8.e+6)
            mu_mean_prof(k) = mu_mean(i, j)
         enddo
      else
         log10NwMean = missing_r4
         mu_mean_prof = missing_r4
      endif
!end    WSO 8/19/13
      !print*, dPRRet%n9(:,i,j)
      if(ialg==2) then 
         call copylognws1t( log10NwMean,dPRRet%n9(:,i,j),i-1)
         call copymus1t( mu_mean_prof, dPRRet%n9(:,i,j),i-1)
         call copypreciptypet(dPRData%raintype(i,j),i-1)
         call copyw10t(w10_out_NS(i,j),i-1) !modified SJM 12/4/2014
         call copyw10sigmat(w10_rms_NS(i,j),i-1) !modified SJM 12/4/2014
         !begin  WSO 8/30/13 update to specify environmental nodes
         call copyenvtemps1t(dPRData%envTemp(:,i,j),env_nodes(:,i),i-1)
         call copysfcairtemps1t(dPRData%envSfcTemp(i,j),i-1)
         call copyenvpresss1t(dPRData%envPress(:,i,j),env_nodes(:,i),i-1)
         call copysfcairpresss1t(dPRData%envSfcPress(i,j),i-1)
         call copyenvqvs1t(dPRData%envQv(:,i,j),env_nodes(:,i),i-1)
         call copyenvsfqvs1t(dPRData%envQv(:,i,j),i-1)
         call copyskintemps1t(dPRData%envSknTemp(i,j),i-1)
         call copyskintempsigmas1t(skintempsigma_NS(i, j), i-1)
         call copycolumnvaporsigmas1t(columnvaporsigma_NS(i, j), i-1)
         call copycolumncloudliqsigmas1t(columncloudliqsigma_NS(i, j), i-1)
         call copyalgotypes1t(algotype_NS(i, j), i-1)
         call copyerrorofdatafits1t(errorofdatafit_NS(i, j), i-1)
      else
         call copylognws1( log10NwMean,dPRRet%n9(:,i,j),i-1)
         call copymus1( mu_mean_prof, dPRRet%n9(:,i,j),i-1)
         call copypreciptype(dPRData%raintype(i,j),i-1)
         call copyw10(w10_out_NS(i,j),i-1) !modified SJM 12/4/2014
         call copyw10sigma(w10_rms_NS(i,j),i-1) !modified SJM 12/4/2014
         call copyenvtemps1(dPRData%envTemp(:,i,j),env_nodes(:,i),i-1)
         call copysfcairtemps1(dPRData%envSfcTemp(i,j),i-1)
         call copyenvpresss1(dPRData%envPress(:,i,j),env_nodes(:,i),i-1)
         call copysfcairpresss1(dPRData%envSfcPress(i,j),i-1)
         call copyenvqvs1(dPRData%envQv(:,i,j),env_nodes(:,i),i-1)
         call copyenvsfqvs1(dPRData%envQv(:,i,j),i-1)
         call copyskintemps1(dPRData%envSknTemp(i,j),i-1)
         call copyskintempsigmas1(skintempsigma_NS(i, j), i-1)
         call copycolumnvaporsigmas1(columnvaporsigma_NS(i, j), i-1)
         call copycolumncloudliqsigmas1(columncloudliqsigma_NS(i, j), i-1)
         call copyalgotypes1(algotype_NS(i, j), i-1)
         call copyerrorofdatafits1(errorofdatafit_NS(i, j), i-1)
      endif
!end    WSO 8/30/13
      !stype = LUT%land_class_map(mod(floor((DPRData%xlon(i,j)+180.)/360.*5760.),5760)+1, 2880-floor((DPRData%xlat(i,j)+90.)/180.*2880.)) !SJM 9/9/15'
      jj=2880-floor((DPRData%xlat(i,j)+90.)/180.*2880.)
      if(jj .lt. 1) jj=1
      if(jj .gt. 2880) jj = 2880
      ii=floor((DPRData%xlon(i,j)+180.)/360.*5760.)+1
      if(ii .lt. 1) ii = 1
      if(ii .gt. 5760) ii = 5760
      if(dPRData%snowIceCover(i,j) .eq. 0) then
        stype = LUT%land_class_map_bare(ii,jj) !SJM 9/9/15
      else
        stype = LUT%land_class_map_snow(ii,jj) !SJM 9/9/15
      endif
      !print*, DPRData%xlon(i,j), DPRData%xlat(i,j), stype
      if(w10_out_NS(i,j)>0 .and. stype .eq. 1) then!wfmap(i,j)>0.9) then
        call calc_relAz(sclonPR(i,j), sclatPR(i,j), DPRData%xlon(i,j), DPRData%xlat(i,j), dprData%envSfcWindU(i,j), dprData%envSfcWindV(i,j), relAz)
        !print*, i,j,w10_out_NS(i,j),'4'
         do jk=1,9
!begin  WSO 10/14/15 flip polarization indices
!            !ipolG(jk),ifreqG(jk)
!            call intplte_emis(ifreqG(jk),1-(ipolG(jk)-1),dPRData%envSknTemp(i,j),w10_out_NS(i,j),emissout(i,j,jk),ebar)
            !ipolG(jk),ifreqG(jk)
            call intplte_emis(ifreqG(jk),ipolG(jk)-1,dPRData%envSknTemp(i,j),w10_out_NS(i,j),relAz,S1eiaPR(i,j),emissoutL(i,j,jk),ebar)
!end    WSO 10/14/14
         enddo
         call intplte_emis(ifreqG(10),0,dPRData%envSknTemp(i,j),w10_out_NS(i,j),relAz,S2eiaPR(i,j),emissoutL(i,j,10),ebar)
         call intplte_emis(ifreqG(11),1,dPRData%envSknTemp(i,j),w10_out_NS(i,j),relAz,S2eiaPR(i,j),emissoutL(i,j,11),ebar)
!          if(dPRData%rainType(i,j)<1 .and.  abs(dPRData%xlat(i, j))>55) then
!             call setEnv(dPRData%envQv(:,i,j),dPRData%envTemp(:,i,j),&
!                  dPRData%envPress(:,i,j),dPRData%envSfcTemp(i,j),&
!                  dPRData%envSknTemp(i,j))
!             call getemissout2(tbRgrid(1:9,i,j+icL),emissout(i,j,1:9))
!             emissout(i,j,10)=emissout(i,j,8)
!             emissout(i,j,11)=emissout(i,j,9)
!             emissout(i,j,12)=emissout(i,j,8)
!             emissout(i,j,13)=emissout(i,j,8)
!          endif
!begin  WSO 10/14/15 replicate emissivities to HF
           emissoutL(i,j,12:13)=emissoutL(i,j,10)
!end    WSO 10/14/15
      else if(dprRet%emis(i,j,1,1,1) .gt. 0.) then
        !print*, dprRet%emis(i,j,1,1,1:nmemb1)
        emissoutL(i,j,1) = sum(dprRet%emis(i,j,1,1,1:nmemb1))/nmemb1
        emissoutL(i,j,2) = sum(dprRet%emis(i,j,2,1,1:nmemb1))/nmemb1
        emissoutL(i,j,3) = sum(dprRet%emis(i,j,1,2,1:nmemb1))/nmemb1
        emissoutL(i,j,4) = sum(dprRet%emis(i,j,2,2,1:nmemb1))/nmemb1
        emissoutL(i,j,5) = sum(dprRet%emis(i,j,1,3,1:nmemb1))/nmemb1
        emissoutL(i,j,6) = sum(dprRet%emis(i,j,1,4,1:nmemb1))/nmemb1
        emissoutL(i,j,7) = sum(dprRet%emis(i,j,2,4,1:nmemb1))/nmemb1
        emissoutL(i,j,8) = sum(dprRet%emis(i,j,1,5,1:nmemb1))/nmemb1
        emissoutL(i,j,9) = sum(dprRet%emis(i,j,2,5,1:nmemb1))/nmemb1
        emissoutL(i,j,10) = sum(dprRet%emis(i,j,1,6,1:nmemb1))/nmemb1
        emissoutL(i,j,11) = sum(dprRet%emis(i,j,2,6,1:nmemb1))/nmemb1
        emissoutL(i,j,12:13) = emissoutL(i,j,10)
        
        emis_rms_NS(i,j,1) = sqrt(sum((dprRet%emis(i,j,1,1,1:nmemb1)-emissoutL(i,j,1))**2)/nmemb1)
        emis_rms_NS(i,j,2) = sqrt(sum((dprRet%emis(i,j,2,1,1:nmemb1)-emissoutL(i,j,2))**2)/nmemb1)
        emis_rms_NS(i,j,3) = sqrt(sum((dprRet%emis(i,j,1,2,1:nmemb1)-emissoutL(i,j,3))**2)/nmemb1)
        emis_rms_NS(i,j,4) = sqrt(sum((dprRet%emis(i,j,2,2,1:nmemb1)-emissoutL(i,j,4))**2)/nmemb1)
        emis_rms_NS(i,j,5) = sqrt(sum((dprRet%emis(i,j,1,3,1:nmemb1)-emissoutL(i,j,5))**2)/nmemb1)
        emis_rms_NS(i,j,6) = sqrt(sum((dprRet%emis(i,j,1,4,1:nmemb1)-emissoutL(i,j,6))**2)/nmemb1)
        emis_rms_NS(i,j,7) = sqrt(sum((dprRet%emis(i,j,2,4,1:nmemb1)-emissoutL(i,j,7))**2)/nmemb1)
        emis_rms_NS(i,j,8) = sqrt(sum((dprRet%emis(i,j,1,5,1:nmemb1)-emissoutL(i,j,8))**2)/nmemb1)
        emis_rms_NS(i,j,9) = sqrt(sum((dprRet%emis(i,j,2,5,1:nmemb1)-emissoutL(i,j,9))**2)/nmemb1)
        emis_rms_NS(i,j,10) = sqrt(sum((dprRet%emis(i,j,1,6,1:nmemb1)-emissoutL(i,j,10))**2)/nmemb1)
        emis_rms_NS(i,j,11) = sqrt(sum((dprRet%emis(i,j,2,6,1:nmemb1)-emissoutL(i,j,11))**2)/nmemb1)
        emis_rms_NS(i,j,12:13) = emis_rms_NS(i,j,10)
        !print '(2i5, 13F8.3)', i, j, emis_rms_NS(i,j,1:13)
      endif
      !if(maxval(emissout(i,j,1:13)) .gt. 2.) then 
      !  print'(A6,2I5,13F8.3)', 'NS Out', i,j,emissout(i,j,1:13)
      !  print*, dprRet%emis(i,j,1,1,1:nmemb1)
      !endif
      !print '(2I5,13F8.3)',i,j,emissout(i,j,1:13)
!begin  WSO 6/5/18 limit emissivities
      do k = 1, 13
        if(emissoutL(i,j,k) .gt. missing_r4) then
          emissoutL(i,j,k) = min(max(emissoutL(i,j,k), surfEmissivity_min), surfEmissivity_max)
        endif
      enddo
!end    WSO 6/5/18 
      if(ialg==2) then
         call copysfcemissouts1t(emissoutL(i,j,:),i-1)
         call copysfcemissouts1sigmat(emis_rms_NS(i,j,:),i-1)
      else
         call copysfcemissouts1(emissoutL(i,j,:),i-1)
         call copysfcemissouts1sigma(emis_rms_NS(i,j,:),i-1)
      endif
      if(i>12 .and. i<38 .and. ialg==1) then
!begin  WSO 9/28/13 add rain flag including missing for bad scans
        call copyrainflags2(dPRData%rainFlagBad(i,j), i-13)
!end    WSO 9/28/13
!begin  WSO 8/15/2014 add ioquality flag
        call copyioqualitys2(dPRData%ioqualityflagdpr(i,j), i-13)
!end    WSO 8/15/2014
!begin  WSO 3/17/17 add snow ice cover flag to output
        call copysnowices2(dPRData%snowIceCover(i,j), i-13)
!end    WSO 3/17/17

!begin  WSO 2/8/17 new output variables
        call copyinitnws2(initnw_MS(:, i, j), dPRRet%n9(:, i, j), i-13)
        call copyprincomps2(princomp_MS(:, i, j), i-13)
        call copyprofclasss2(profclass_MS(i, j), i-13)
        call copysurfprecipbiasratios2(surfprecipbiasratio_MS(i, j), i-13)
        call copysubfootvariabilitys2(subfootvariability_MS(i, j), i-13)
        call copymultiscatcalcs2(multiscatcalc_MS(i, j), i-13)
        call copymultiscatsurfaces2(multiscatsurface_MS(i, j), i-13)
!end    WSO 2/8/17
!begin  WSO 9/15/13 add missing for cloud ice profiles
        call copycldwaters2(cldwprof, i-13)
        cldiprof = missing_r4
        call copycldices2(cldiprof,i-13)
!end  WSO 9/15/13 
         call copyzckus2(zcKu3DMS(:,i,j), zcKa3DMS(:,i,j), i-13)
!begin  WSO 9/5/13 add DSRT sigma PIA information
         call copysigmapias2(dPRData%dsrtsigmaPIAku(i, j), dPRData%dsrtsigmaPIAka(i, j), i-13)
!end    WSO 9/5/13
!begin  WSO 2/8/17 output sigmaZero at Ka band
         call copysigmazeros2(dPRData%sigmaZeroKu(i, j), dPRData%sigmaZeroKa(i, j), i-13)
!end    WSO 2/8/17
!begin  WSO 8/19/13 change dNw to Nw and add mu
         if(dPRData%rainType(i,j)>=100) then
            do k=1,88
               log10NwMean(k)=sum(dPRRet%MS%log10dNw(1:nmemb1,k,i,j))/nmemb1 + &
                log10(8.e+6)
               mu_mean_prof(k) = mu_meanMS(i, j)
            enddo
         else
            log10NwMean = missing_r4
            mu_mean_prof = missing_r4
         endif
!end     WSO 8/19/13
         call copylognws2( log10NwMean,dPRRet%n9(:,i,j),i-13)
         call copymus2( mu_mean_prof,dPRRet%n9(:,i,j),i-13)
         call copyrrates2(rrate3DMS(:,i,j),rrate3DstdMS(:,i,j),i-13)
         call copypwcs2(pwc3DMS(:,i,j),pwc3DstdMS(:,i,j),i-13)
!begin  WSO 8/7/13
         call copylwcfracs2(mlwc_fracMS(:,i,j),mrate_fracMS(:,i,j),i-13)
         call copysfcrainliqfracs2(sfcRainLiqFracMS(i, j), i-13)
!end    WSO 8/7/13
         call copyd0s2(d03DMS(:,i,j),i-13)
         call copynodess2(dPRData%node(:,i,j),i-13)
!begin  WSO 8/30/13 update to specify environmental nodes
         call copyenvtemps2(dPRData%envTemp(:,i,j), env_nodes(:,i), i-13)
!         call copyenvsftemps2(dPRData%envTemp(i,j),i-13)
         call copysfcairtemps2(dPRData%envSfcTemp(i,j),i-13)
         call copyenvpresss2(dPRData%envPress(:,i,j), env_nodes(:,i), i-13)
         call copysfcairpresss2(dPRData%envSfcPress(i,j),i-13)
         call copyenvqvs2(dPRData%envQv(:,i,j), env_nodes(:,i), i-13)
         call copyenvsfqvs2(dPRData%envQv(:,i,j),i-13)
         call copyskintemps2(dPRData%envSknTemp(i,j),i-13)
         call copyskintempsigmas2(skintempsigma_MS(i, j), i-13)
         call copycolumnvaporsigmas2(columnvaporsigma_MS(i, j), i-13)
         call copycolumncloudliqsigmas2(columncloudliqsigma_MS(i, j), i-13)
         call copyalgotypes2(algotype_MS(i, j), i-13)
         call copyerrorofdatafits2(errorofdatafit_MS(i, j), i-13)
!end    WSO 8/30/13
         if(w10_out_MS(i,j)>0 .and. wfmap(i,j)>0.9) then
            !print*, i,j,w10_out_MS(i,j),'5'
            call calc_relAz(sclonPR(i,j),sclatPR(i,j), DPRData%xlon(i,j), DPRData%xlat(i,j), dprData%envSfcWindU(i,j), dprData%envSfcWindV(i,j), relAz)
            do jk=1,9
               !ipolG(jk),ifreqG(jk)
!begin  WSO 10/14/15 flip polarization index
!               call intplte_emis(ifreqG(jk),&
!                    1-(ipolG(jk)-1),dPRData%envSknTemp(i,j),&
!                    w10_out_MS(i,j),emissout(i,j,jk),ebar)
               call intplte_emis(ifreqG(jk),&
                    ipolG(jk)-1,dPRData%envSknTemp(i,j),&
                    w10_out_MS(i,j),relAz,S1eiaPR(i,j),emissoutL(i,j,jk),ebar)
!end    WSO 10/14/15
            enddo
            call intplte_emis(ifreqG(10),0,dPRData%envSknTemp(i,j),w10_out_NS(i,j),relAz,S2eiaPR(i,j),emissoutL(i,j,10),ebar)
            call intplte_emis(ifreqG(11),1,dPRData%envSknTemp(i,j),w10_out_NS(i,j),relAz,S2eiaPR(i,j),emissoutL(i,j,11),ebar)
!begin  WSO 10/14/15 replicate emissivities at HF
            emissoutL(i,j,12:13)=emissoutL(i,j,10)
!end    WSO 10/14/15
         else if(dprRet%MS%emis(i,j,1,1,1) .gt. 0.) then
           emissoutL(i,j,1) = sum(dprRet%MS%emis(i,j,1,1,1:nmemb1))/nmemb1
           emissoutL(i,j,2) = sum(dprRet%MS%emis(i,j,2,1,1:nmemb1))/nmemb1
           emissoutL(i,j,3) = sum(dprRet%MS%emis(i,j,1,2,1:nmemb1))/nmemb1
           emissoutL(i,j,4) = sum(dprRet%MS%emis(i,j,2,2,1:nmemb1))/nmemb1
           emissoutL(i,j,5) = sum(dprRet%MS%emis(i,j,1,3,1:nmemb1))/nmemb1
           emissoutL(i,j,6) = sum(dprRet%MS%emis(i,j,1,4,1:nmemb1))/nmemb1
           emissoutL(i,j,7) = sum(dprRet%MS%emis(i,j,2,4,1:nmemb1))/nmemb1
           emissoutL(i,j,8) = sum(dprRet%MS%emis(i,j,1,5,1:nmemb1))/nmemb1
           emissoutL(i,j,9) = sum(dprRet%MS%emis(i,j,2,5,1:nmemb1))/nmemb1
           emissoutL(i,j,10) = sum(dprRet%MS%emis(i,j,1,6,1:nmemb1))/nmemb1
           emissoutL(i,j,11) = sum(dprRet%MS%emis(i,j,2,6,1:nmemb1))/nmemb1
           emissoutL(i,j,12:13) = emissoutL(i,j,10)
           
           emis_rms_MS(i,j,1) = sqrt(sum((dprRet%MS%emis(i,j,1,1,1:nmemb1)-emissoutL(i,j,1))**2)/nmemb1)
           emis_rms_MS(i,j,2) = sqrt(sum((dprRet%MS%emis(i,j,2,1,1:nmemb1)-emissoutL(i,j,2))**2)/nmemb1)
           emis_rms_MS(i,j,3) = sqrt(sum((dprRet%MS%emis(i,j,1,2,1:nmemb1)-emissoutL(i,j,3))**2)/nmemb1)
           emis_rms_MS(i,j,4) = sqrt(sum((dprRet%MS%emis(i,j,2,2,1:nmemb1)-emissoutL(i,j,4))**2)/nmemb1)
           emis_rms_MS(i,j,5) = sqrt(sum((dprRet%MS%emis(i,j,1,3,1:nmemb1)-emissoutL(i,j,5))**2)/nmemb1)
           emis_rms_MS(i,j,6) = sqrt(sum((dprRet%MS%emis(i,j,1,4,1:nmemb1)-emissoutL(i,j,6))**2)/nmemb1)
           emis_rms_MS(i,j,7) = sqrt(sum((dprRet%MS%emis(i,j,2,4,1:nmemb1)-emissoutL(i,j,7))**2)/nmemb1)
           emis_rms_MS(i,j,8) = sqrt(sum((dprRet%MS%emis(i,j,1,5,1:nmemb1)-emissoutL(i,j,8))**2)/nmemb1)
           emis_rms_MS(i,j,9) = sqrt(sum((dprRet%MS%emis(i,j,2,5,1:nmemb1)-emissoutL(i,j,9))**2)/nmemb1)
           emis_rms_MS(i,j,10) = sqrt(sum((dprRet%MS%emis(i,j,1,6,1:nmemb1)-emissoutL(i,j,10))**2)/nmemb1)
           emis_rms_MS(i,j,11) = sqrt(sum((dprRet%MS%emis(i,j,2,6,1:nmemb1)-emissoutL(i,j,11))**2)/nmemb1)
           emis_rms_MS(i,j,12:13) = emis_rms_MS(i,j,10)
           !print '(2i5, 13F8.3)', i, j, emis_rms_MS(i,j,1:13)
         endif
         !print '(2I5,13F8.3)',i,j,emissout(i,j,1:13)
         !if(maxval(emissout(i,j,1:13)) .gt. 2.) print'(A6,2I5,13F8.3)', 'MS Out ',i,j,emissout(i,j,1:13)
!begin  WSO 6/5/18 limit emissivities
         do k = 1, 13
           if(emissoutL(i,j,k) .gt. missing_r4) then
             emissoutL(i,j,k) = min(max(emissoutL(i,j,k), surfEmissivity_min), surfEmissivity_max)
           endif
         enddo
!end    WSO 6/5/18 
         call copysfcemissouts2(emissoutL(i,j,:),i-13)
         call copysfcemissouts2sigma(emis_rms_MS(i,j,:),i-13)
      endif

      do k=1,9
         if(tbout2d(i,j,k)>0) then
            tbout(k)=tbout2d(i,j,k)!-tbRgrid(k,i,j+icL)
         else
            tbout(k)=tbout2dNoOcean(i,j,k)!-tbRgrid(k,i,j+icL)
         endif
      enddo
      tbout(10)=sum(dPRRet%tb(i,j,1,6,1:1*nmemb1))/(nmemb1)
      tbout(11)=sum(dPRRet%tb(i,j,2,6,1:1*nmemb1))/(nmemb1)
      tbout(12)=sum(dPRRet%tb(i,j,1,7,1:1*nmemb1))/(nmemb1)
      tbout(13)=sum(dPRRet%tb(i,j,1,8,1:1*nmemb1))/(nmemb1)
      do k=10,13
         !tbout(k)=hFreqPRg(i,j,k-9)
      enddo
      !print*, tbout
      if(dPRData%rainType(i,j)>=100) then
         !print*, tbout
      endif
      if(tbout(7)>0 .and. tbRgrid(7,i,j+icL)>0) then
         !write(*,10) tbout(6:9), tbRgrid(6:9,i,j+icL)
      endif
      if(ialg==2) then
         call copytbouts1t(tbout,i-1)
      else
         call copytbouts1(tbout,i-1)
      endif
!!begin MG 09172013
      if(i>12 .and. i<38) then
         do k=1,9
            if(tbout2dMS(i,j,k)>0) then
               tbout(k)=tbout2dMS(i,j,k)!-tbRgrid(k,i,j+icL)
            else
               tbout(k)=tbout2dNoOceanMS(i,j,k)!-tbRgrid(k,i,j+icL)
            endif
         enddo
         tbout(10)=sum(dPRRet%MS%tb(i,j,1,6,1:1*nmemb1))/(nmemb1)
         tbout(11)=sum(dPRRet%MS%tb(i,j,2,6,1:1*nmemb1))/(nmemb1)
         tbout(12)=sum(dPRRet%MS%tb(i,j,1,7,1:1*nmemb1))/(nmemb1)
         tbout(13)=sum(dPRRet%MS%tb(i,j,1,8,1:1*nmemb1))/(nmemb1)
         call copytbouts2(tbout,i-13)
      endif
!!end MG 09172013
   enddo

   if(ialg==2) then
      call frominputt(st_2adpr)
      call copyscantimet(j-1)
      call writescant()
   else
      call frominput(st_2adpr)
      call copyscantime(j-1)
      call writescan()
   endif
  
   
 
!begin  WSO 9/5/13 suppress writing to 2ADPR
!   if(ifdpr(1:1)=='Y') then
!      call writedprscan()
!      !print*, ifdpr
!   endif
!end    WSO 9/5/13 
enddo


!do j=3,dPRData%n1c21-3
!   do i=3,49-2
!      if(dPRData%raintype(i,j)>0 .and. minval(tbout2d(i,j,:))>0 .and. &
!           minval(tbRgrid(:,i,ic+j))>0.) then
!         write(*,131) tbRgrid(1:9,i,j+ic),tbout2d(i,j,1:9)
!      endif
!   enddo
!enddo
! SFM  begin  03/27/2014; execution protections
  IF (ALLOCATED(emissoutL)) deallocate(emissoutL)
  if(allocated(hFreqPRg)) deallocate(hFreqPRg)
! SFM  end    03/27/2014
!print*, ntbpix, ntbpix2
!print*, ifdpr(1:1)=='Y'
101 format(10(F8.3,1x),I3,2(F8.3,1x))

!000000000000000000000000000000000000000000000000000000000000000000000000000000
! SFM  begin  03/27/2014; execution protections
  IF (ALLOCATED(Yobs)) deallocate(Yobs)
  IF (ALLOCATED(Xup)) deallocate(Xup)
  IF (ALLOCATED(randemiss)) deallocate(randemiss)
  IF (ALLOCATED(Xens)) deallocate(Xens)
  IF (ALLOCATED(Yens)) deallocate(Yens)
! SFM  end    03/27/2014

  call deallocGeophys()
  call deallocateStormStructData(stormStruct)
  call deallocateDPRProfRet(radarRet)
  call deallocateDPRProfData(radarData)
  
! SFM  begin  03/27/2014; execution protections
  IF (ALLOCATED(ndn)) deallocate(ndn) 
  IF (ALLOCATED(xscalev)) deallocate(xscalev) 
  IF (ALLOCATED(logdnwf)) deallocate(logdnwf) 
  IF (ALLOCATED(ndnp)) deallocate(ndnp) 
  IF (ALLOCATED(rhPCij)) deallocate(rhPCij)
  IF (ALLOCATED(cldwPCij)) deallocate(cldwPCij)
! SFM  end    03/27/2014

end subroutine radarRetSub

subroutine setRetParam(retParam)
  use f90Types
  implicit none
  type (retParamType)    :: retParam
  retParam%wz=1.
  retParam%w13=1.
  retParam%w35=1.
  retParam%z13thresh=15
  retParam%z35thresh=15
end subroutine setRetParam

subroutine readpcoeff(vLand,vOcean,pMLand,pMOcean,mTbLand,mTbOcean,            &
                      stTbLand,stTbOcean)
  real  :: vLand(18,18), vOcean(10,10)
  real  :: pMLand(18), pMOcean(10)
  real  :: mTbLand(9), mTbOcean(9)
  real  :: stTbLand(9), stTbOcean(9)

open(10,file='AncData/predEofs.Land.txt')
do i=1,18
   read(10,*) (vLand(i,j),j=1,18)
enddo
close(10)
open(10,file='AncData/predMeans.Land.txt')
read(10,*) (pMLand(i),i=1,18)
close(10)
open(10,file='AncData/stdCSTb.Land.txt')
read(10,*) (stTbLand(i),i=1,9)
close(10)
open(10,file='AncData/meanCSTb.Land.txt')
read(10,*) (mTbLand(i),i=1,9)
close(10)

open(10,file='AncData/envpredEofs.Ocean.txt')
do i=1,10
   read(10,*) (vOcean(i,j),j=1,10)
enddo
close(10)
open(10,file='AncData/envpredMeans.Ocean.txt')
read(10,*) (pMOcean(i),i=1,10)
close(10)
open(10,file='AncData/envstdCSTb.Ocean.txt')
read(10,*) (stTbOcean(i),i=1,9)
close(10)
open(10,file='AncData/envmeanCSTb.Ocean.txt')
read(10,*) (mTbOcean(i),i=1,9)
close(10)
end subroutine readpcoeff
