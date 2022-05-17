!  SFM  04/06/2013  Code changes from M.Grecui
!
subroutine radarRetSubEmpty(geoData, dPRData, dPRRet, gridEnvData, gmi2Grid,        &
      gmiData, nmu2, nMemb, nmfreq2, ngates, sysdN, ic, tbRgrid,               &
      dprrain,ichunk,st_2adpr,orbNumb,istart,iend)
!  SFM  end    12/13/2013
  use f90DataTypes
  use f90Types
  use cldclass
  use ran_mod
  use geophysEns
  use nbinMod
  use tables2
  use weight
  Use BMCVparameters
  use emissMod
!begin  MG 10/29/15 add gEnv module
  use gEnv
!end    MG 10/29/15
!begin  WSO 9/14/13 incorporate missing flags
  use missingMod
!end    WSO 9/14/13

  implicit none
  type (geoDataType) :: geoData
  type (gridEnvDataType) :: gridEnvData
  type (dPRDataType)     :: dPRData
  type (dPRRetType)      :: dPRRet
  type (gmi2GridType)    :: gmi2Grid
  type (cgMIDataType)     :: gmiData
  integer :: nmu2, nMemb, nmfreq2, ngates
  integer*4 :: ichunk
  integer :: st_2adpr              ! file open status for 2adpr file
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

  integer :: ic, ii, jj, iGMI, jGMI
  integer :: di(8), dj(8)  
  integer :: i, j, k, ig, jg, ntpw, nmemb1, itop, irand
  real    :: pia13m, rms1, rms2, sysdN, unSortedRR(200), corrcoef, sfcRain2
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
  real :: w10(49,300), w10_out_NS(49,300), w10_out_MS(49,300), w10_min, w10_max
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
  real  :: vLand(18,18), vOcean(10,10)
  real  :: pMLand(18), pMOcean(10)
  real  :: mTbLand(9), mTbOcean(9)
  real  :: stTbLand(9), stTbOcean(9)
  double precision  :: xin(18), xpred(18), yout(9)
  real, allocatable :: emissout(:,:,:)
!begin  WSO 8/19/13 change Nw variable name (not dN) and add mu
  real :: cldwprof(88), cldiprof(88), log10NwMean(88), mu_mean_prof(88)
  integer *2 :: env_nodes(10, 49)
  real :: env_levs(10), ray_angle, pi
!end    WSO 8/19/13
  real :: lFract(49,300), sprobs, probs(100), rmsS(100)
  real :: covar
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
  real :: mlwc_frac(10, 49, 300)
  real :: mrate_frac(10, 49, 300)
  real :: mlwc_fracMS(10, 49, 300)
  real :: mrate_fracMS(10, 49, 300)
  real :: sfcRainLiqFrac(49, 300)
  real :: sfcRainLiqFracMS(49, 300)
  real :: tbMax1(15), tbMin1(15)

!  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
!  SFM  begin  06/22/2014
  real :: wfractPix, windPert(100), qvPert(100), dnqv
!  SFM  end    06/22/2014
!  SFM  end    07/29/2014
!end   WSO 8/8/13
  real :: covTb(49,300,15,15), tbMax(49,300,15), tbMin(49,300,15), &
       tbMean(49,300,15)
  real :: invCovTb(49,300,15,15)
  real :: tbout2D(49,300,15), tb(49,300,15), tbNoOcean(49,300,15), &
       tbout2DNoOcean(49,300,15), tbObs(49,300,15), fobj
  real :: dfdtb(49,300,15), rerr(15), tb0(49,300,15), fem(15) , &
       tb0MS(49,300,15), tbNoOceanMS(49,300,15), tbout2DNoOceanMS(49,300,15),&
       tbout2DMS(49,300,15)

  integer :: actOb(49,300), iactOb
  integer :: jk, nf
  integer :: dig               ! SFM  04/16/2014  for M.Grecu
  real   :: cl(9,25), xin25(25),dtb(9)
  real   :: ebar, minl
  real, allocatable :: geoloc(:), hFreqTbs(:,:), PRgeoloc(:), hFreqPRg(:,:,:)
  integer*4 :: istart, iend
  integer :: iconv
  real :: nubfc, stdpia35
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
  dPRRet%cldwcoeff=0
  ifreqG(1:13)=(/1,1,2,2,3,4,4,5,5,6,6,7,8/)
  ipolG(1:13)=(/1,2,1,2,1,1,2,1,2,1,2,1,1/)
! 1 is V
! 2 is H
  !call readclust()
!begin  MG Sept 15 2015
  allocate(hFreqPRg(49,dPRData%n1c21,4))
  hFreqPRg=missing_r4
!end  MG
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
       logdNwf(9*nmemb1), randemiss(nmfreq*nmemb1*2)) 
  allocate(ndnp(10,nmemb1)) 
  allocate(emissout(49,dPRData%n1c21,15))

  ntbpix=0
  ntbpix2=0
!  begin  SFM  07/29/2014; for M.Grecu,  eliminate NANs
  do k=0,nmemb1-1
     !windPert(k+1)=20*ran1()
     windPert(k+1) = normal2(0.,.15) !SJM 2/4/15
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
  do k=1,2*nmfreq*nmemb1
     randemiss(k)=.5*ran1()
  enddo

 

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
!end    WSO 9/15/13
!end WSO 04/07/13
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
  print*, nmfreq
  call allocateDPRProfRet(radarRet,nmfreq,nmemb1,ngates, 9)   ! allocates memory for the 1D 

  !...retrieval structures
  radarRet%rrate = 0.0  

  radarRet%tb=-99

  call allocateDPRProfData(radarData, ngates)                 ! allocates memory for the 
                                                              ! 1-D DPR observations
  call allocateStormStructData(stormStruct)                   ! allocates memory for the 5-node 
                                                              ! storm structure
  
  dPRRet%sfc_wind(1:nmemb)=radarRet%sfc_wind(1:nmemb)
  
  dPRRet%sfcRainEns=0
  stormStruct%iSurf=ngates
  radarData%ngates=ngates
  radarData%dr=0.25
  dPRRet%tb=-99
  dPRRet%emtb=-99
  dPRRet%n9=0
  call allocGeophys(6,61,9,nmemb1,nmfreq*nmemb1*2)
  call setdNwIcJcL(sysdN,nmemb1)
  nx=1+nbin*7
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
  emissout=missing_r4
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
!end SJM 12/9/2014

!  SFM  end    06/22/2014
do j=1,dPRData%n1c21
   do i=1,49
      sfcRain(i,j)=missing_r4
      sfcRainStd(i,j)=missing_r4
      rrate3D(1:nbin,i,j) = missing_r4
      rrate3Dstd(1:nbin,i,j) = missing_r4
      pwc3D(1:nbin,i,j) = missing_r4
      pwc3Dstd(1:nbin,i,j) = missing_r4
      d03D(1:nbin,i,j) = missing_r4
   end do
end do

do j=1,dPRData%n1c21
   do i=13,37
      sfcRainMS(i,j)=missing_r4
      sfcRainStdMS(i,j)=missing_r4
      rrate3DMS(1:nbin,i,j) = missing_r4
      rrate3DstdMS(1:nbin,i,j) = missing_r4
      pwc3DMS(1:nbin,i,j) = missing_r4
      pwc3DstdMS(1:nbin,i,j) = missing_r4
      d03DMS(1:nbin,i,j) = missing_r4
   enddo
enddo


call rewind(ic)
IF (st_2adpr .EQ. 0) call rewind_dpr(ic)
!print*, minval( dPRRet%cldwcoeff), maxval( dPRRet%cldwcoeff)  !MG 09/26/13

do j=1,dPRData%n1c21
   do i=1,25
      call copyw10small(w10_out_MS(i+12,j),i-1) !SJM 12/9/2014
   enddo

   call setlatlons1( dPRData%xlat(:,j), dPRData%xlon(:,j),                     &
        sfcRain(:,j),sfcRainStd(:,j),piaOut(:,j))
  
   call setlatlons2( dPRData%xlat(:,j), dPRData%xlon(:,j),                     &
        sfcRainMS(:,j),sfcRainStdMS(:,j),piaOutKuMS(:,j),piaOutKaMS(:,j))


   do i=1,49

     call copyrainflags1(dPRData%rainFlagBad(i,j), i-1)

     call copyioqualitys1(dPRData%ioqualityflagku(i,j), i-1)

     call copysigmapias1(dPRData%srtsigmaPIAku(i, j), i-1)
     call copyrrates1(rrate3D(:,i,j),rrate3Dstd(:,i,j),i-1)
     call copypwcs1(pwc3D(:,i,j),pwc3Dstd(:,i,j),i-1)
     call copylwcfracs1(mlwc_frac(:,i,j),mrate_frac(:,i,j),i-1)
     call copysfcrainliqfracs1(sfcRainLiqFrac(i, j), i-1)
     call copyd0s1(d03D(:,i,j),i-1)
     call copyzckus1(zcKu3D(:,i,j),i-1)
     call copynodess1(dPRData%node(:,i,j),i-1)
     call copylognws1( log10NwMean,dPRRet%n9(:,i,j),i-1)
     call copymus1( mu_mean_prof, dPRRet%n9(:,i,j),i-1)
     call copypreciptype(dPRData%raintype(i,j),i-1)
     call copyw10(w10_out_NS(i,j),i-1) !modified SJM 12/4/2014
     call copyenvtemps1(dPRData%envTemp(:,i,j),env_nodes(:,i),i-1)
     call copysfcairtemps1(dPRData%envSfcTemp(i,j),i-1)
     call copyenvpresss1(dPRData%envPress(:,i,j),env_nodes(:,i),i-1)
     call copysfcairpresss1(dPRData%envSfcPress(i,j),i-1)
     call copyenvqvs1(dPRData%envQv(:,i,j),env_nodes(:,i),i-1)
     call copyenvsfqvs1(dPRData%envQv(:,i,j),i-1)
     call copyskintemps1(dPRData%envSknTemp(i,j),i-1)
     
   enddo
!  SFM  begin  12/12/2013, add file flag to calling sequence
   call frominput(st_2adpr)
!  SFM  end    12/12/2013
   call copyscantime(j-1)
   call writescan()

enddo


  IF (ALLOCATED(emissout)) deallocate(emissout)
  if(allocated(hFreqPRg)) deallocate(hFreqPRg)


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

end subroutine radarRetSubEmpty

