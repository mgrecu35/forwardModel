module local_RD_var
  use f90Types
  real :: covTb(49,300,15,15), tbMax(49,300,15), tbMin(49,300,15), &
       tbMean(49,300,15)
  real :: invCovTb(49,300,15,15)
  real :: tbout2D(49,300,15), tb(49,300,15), tbNoOcean(49,300,15), &
       tbout2DNoOcean(49,300,15), tbObs(49,300,15)
  real :: dfdtb(49,300,15), rerr(15), tb0(49,300,15), fem(15) , &
       tb0MS(49,300,15), tbNoOceanMS(49,300,15), tbout2DNoOceanMS(49,300,15),&
       tbout2DMS(49,300,15)
  integer :: ipol(15), ifreq(15), iobs(15), ifreq1
  integer :: ifreqG(15), ipolG(15) 
  integer,parameter :: nbinL=88
  real :: sfcRain(49,300),sfcRainStd(49,300)
  real :: rRate3D(nbinL,49,300),  rRate3Dstd(nbinL,49,300)
  real :: pwc3D(nbinL,49,300),  pwc3Dstd(nbinL,49,300)
  real :: zcKu3D(nbinL,49,300), d03D(nbinL,49,300), piaOut(49,300)

  real :: sfcRainMS(49,300),sfcRainStdMS(49,300),pia35m(49,300)
  real :: rRate3DMS(nbinL,49,300),  rRate3DstdMS(nbinL,49,300)
  real :: pwc3DMS(nbinL,49,300),  pwc3DstdMS(nbinL,49,300)
  real :: zcKu3DMS(nbinL,49,300), zcKa3DMS(nbinL,49,300), d03DMS(nbinL,49,300), &
          piaOutKuMS(49,300), piaOutKaMS(49,300)
  real, allocatable :: emissoutL(:,:,:), emis_out_NS(:,:,:), emis_out_MS(:,:,:) !sjm 8/10/15
  type (stormStructType) :: stormStruct
  real, allocatable ::  Yens(:,:), Xens(:,:), Yobs(:), Xup(:)
  real, allocatable :: geoloc(:), hFreqTbs(:,:), PRgeoloc(:), hFreqPRg(:,:,:)
  real, allocatable :: ndn(:), ndnp(:,:), xscalev(:), logdNwf(:), randemiss(:), dwind(:)
  real, allocatable  :: rhPCij(:,:), cldwPCij(:,:)
  type (radarRetType)    :: radarRet
  type (radarDataType)   :: radarData
  real :: emis_eofs(100,12) !SJM 7/9/2015
  real :: scLatPR(49,300),scLonPR(49,300),wfmap(49,300), fpmap(49,300,15), fpmapN(49,300,15)
  real :: w10(49,300), w10_out_NS(49,300), w10_out_MS(49,300), w10_min, w10_max, emis, relAz
  real :: w10_rms_NS(49,300), emis_rms_NS(49,300,13), w10_rms_MS(49,300), emis_rms_MS(49,300,13)
  real :: S1eiaPR(49,300), S2eiaPR(49,300)
  real :: sigmaZeroVarKu(49,300), sigmaZeroVarKa(49,300), sigmaZeroCov(49,300)
  integer :: iiad
end module local_RD_var

subroutine get_rain_type(raintype, nscans)
  use globalData
  integer :: raintype(49,300)
  integer :: nscans
  raintype=dPRData%rainType
  nscans=dPRData%n1c21
  print*, nscans
end subroutine get_rain_type

subroutine testr(it)
  integer :: it
  it=136
end subroutine testr

subroutine radarRetSub2(nmu2,  nmfreq2,   icL, tbRgrid,               &
      dprrain,ichunk,orbNumb,ialg,idir)
!  SFM  end    12/13/2013
  use local_RD_var
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
 
  !type (stormStructType) :: stormStruct
  type (retParamType)    :: retParam
  real :: xin363(363)
 
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

  integer :: ibatch
  real    :: stddev, srtpiaf
  real    :: FWHMx, FWHMy, tbconv(2), tbconvEns(2,100)
  integer :: dnx,dny, ik
  
  real :: cldw(nlayer), rh(nlayer), pia13s
  integer :: nx,ny, icount, imin
  real ::  xm1,xs,rmsmin, prob, probtot, rmstot
  real :: piaR(100), fPIA, z13m
  integer :: ntbpix, ntbpix2
  real :: emtbm(9)
  real :: zminsc
  real :: realOut(49)
 
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
  !integer :: ifreqG(15), ipolG(15) 
  DOUBLE PRECISION input(6)
  DOUBLE PRECISION output(2)
  real :: wfract(5,5), wfractm, wfractsd
  real                    emissv(n_chan)
  real                    emissh(n_chan)
  real                    emissv_std(n_chan)
  real                    emissh_std(n_chan)
  integer :: stype!SJM 7/9/2015
  real  :: vLand(18,18), vOcean(10,10)
  real  :: pMLand(18), pMOcean(10)
  real  :: mTbLand(9), mTbOcean(9)
  real  :: stTbLand(9), stTbOcean(9)
  double precision  :: xin(18), xpred(18), yout(9)

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
 
  !end SJM 7/25/2014
!begin WSO 8/8/13
  real :: gatelength
  real :: depthBB, depthML, depth
  real :: mu_mean(49, 300)
  real :: mu_meanMS(49, 300)

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
 

  integer :: actOb(49,300), iactOb
  integer :: jk, nf
  integer :: dig               ! SFM  04/16/2014  for M.Grecu
  real   :: cl(9,25), xin25(25),dtb(9)
  real   :: ebar, minl
  !real, allocatable :: geoloc(:), hFreqTbs(:,:), PRgeoloc(:), hFreqPRg(:,:,:)
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
  print*, orbNumb,ichunk
  !call openascii(orbNumb,ichunk)
  
  print*, dPRData%n1c21
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
     allocate(prgeoloc(dPRData%n1c21*49*2))
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


     do k=1,4
        call flannint(geoloc,  prgeoloc, hFreqTbs(:,k), hFreqPRg(:,:,k), &
             (iEnd+1-iStart)*81, 2, &
             dPRData%n1c21*49)
     enddo

     deallocate(geoloc)
     deallocate(prgeoloc)
     deallocate(hFreqTbs)
  end if

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

! SFM  end    03/27/2014

  print*, nx, nmemb1
  allocate(Xens(nx,nmemb1),Yens(ny,nmemb1),Yobs(ny), Xup(nx))
  
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
  print*,'max_gmiData%tpw3=', maxval(gmiData%tpw3)
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
  rrate3D=0!-99
  tbNoOcean=-99
  pia35m=0.
  print*, gmi2Grid%xmin, gmi2Grid%ymin
  !return 
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
           !print*, iGMI,jGMI, i, j, gmiData%tpw3(iGMI,jGMI), gmi2Grid%xmin, icL  
           if(jGMI>-99) then
              scLonPR(i,j)=gmiData%SCLon3(jGMI)
              scLatPR(i,j)=gmiData%SCLat3(jGMI)
              S1eiaPR(i,j)=gmidata%S1eia3(iGMI,jGMI)
              S2eiaPR(i,j)=gmidata%S2eia3(iGMI,jGMI)
              !print*,i,j,S1eiaPR(i,j), S2eiaPR(i,j)
           endif
           
           if(iGMI>0 .and. jGMI>0 .and. gmi2Grid%xmin>-998 ) then !4/22/14 MG
              if(gmiData%tpw3(iGMI,jGMI)>0) then
                 if(gmiData%sfc_wind3(iGMI,jGMI)>0) then
                    radarRet%sfc_wind=gmiData%sfc_wind3(iGMI,jGMI)
                    tbRgrid(14,i,j+icL)=gmiData%tpw3(iGMI,jGMI) 
                 else
                    radarRet%sfc_wind=dPRData%envSfcWind(i,j)
                    tbRgrid(14,i,j+icL)=-99
                 endif
                 !print*, ifdpr(1_, iftest(1), icL
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
                    if(dprData%binRealSurface(i,j)>88) then
                       print*, dprData%binRealSurface(i,j), i, j
                       stop
                    end if
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
                    
                    !print*,dPRData%node(:,i,j)
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
121                 format(81(F8.2,1x))
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
  print*, maxval(dPRRet%log10dNw)
  print*, minval(dPRRet%log10dNw)
  !print*,tbRgrid(7,:,1:300)
  !stop
end subroutine radarRetSub2

subroutine radarRetSub3(nmu2,  nmfreq2,   icL, tbRgrid,               &
      dprrain,ichunk,orbNumb,ialg,idir)
!  SFM  end    12/13/2013
  use local_RD_var
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
  !type (radarRetType)    :: radarRet
  !type (radarDataType)   :: radarData
  type (retParamType)    :: retParam
  real :: xin363(363)
 
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
  !real, allocatable ::  Yens(:,:), Xens(:,:), Yobs(:), Xup(:)
  integer :: ibatch
  real    :: stddev, srtpiaf
  real    :: FWHMx, FWHMy, tbconv(2), tbconvEns(2,100)
  integer :: dnx,dny, ik
  !real, allocatable :: ndn(:), ndnp(:,:), xscalev(:), logdNwf(:), randemiss(:), dwind(:)
  !real, allocatable  :: rhPCij(:,:), cldwPCij(:,:)
  real :: cldw(nlayer), rh(nlayer), pia13s
  integer :: nx,ny, icount, imin
  real ::  xm1,xs,rmsmin, prob, probtot, rmstot
  real :: piaR(100), fPIA, z13m
  integer :: ntbpix, ntbpix2
  real :: emtbm(9)
  real :: zminsc
  real :: realOut(49)
  !real :: w10(49,300), w10_out_NS(49,300), w10_out_MS(49,300), w10_min, w10_max, emis, relAz
  !real :: w10_rms_NS(49,300), emis_rms_NS(49,300,13), w10_rms_MS(49,300), emis_rms_MS(49,300,13)
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
  DOUBLE PRECISION input(6)
  DOUBLE PRECISION output(2)
  real :: wfract(5,5), wfractm, wfractsd
  real                    emissv(n_chan)
  real                    emissh(n_chan)
  real                    emissv_std(n_chan)
  real                    emissh_std(n_chan)
  integer :: stype!SJM 7/9/2015
  real  :: vLand(18,18), vOcean(10,10)
  real  :: pMLand(18), pMOcean(10)
  real  :: mTbLand(9), mTbOcean(9)
  real  :: stTbLand(9), stTbOcean(9)
  double precision  :: xin(18), xpred(18), yout(9)
  !real, allocatable :: emissoutL(:,:,:), emis_out_NS(:,:,:), emis_out_MS(:,:,:) !sjm 8/10/15
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
  !real :: sigmaZeroVarKu(49,300), sigmaZeroVarKa(49,300), sigmaZeroCov(49,300)
  !end SJM 7/25/2014
!begin WSO 8/8/13
  real :: gatelength
  real :: depthBB, depthML, depth
  real :: mu_mean(49, 300)
  real :: mu_meanMS(49, 300)
!  real :: scLatPR(49,300),scLonPR(49,300),wfmap(49,300), fpmap(49,300,15), fpmapN(49,300,15)
!  real :: S1eiaPR(49,300), S2eiaPR(49,300)
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
  

  integer :: actOb(49,300), iactOb
  integer :: jk, nf
  integer :: dig               ! SFM  04/16/2014  for M.Grecu
  real   :: cl(9,25), xin25(25),dtb(9)
  real   :: ebar, minl

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
  !print*, tbout2D(24,:,1)
  !stop
  call clearsc()
  call asciiplot(sfcRain(:,1:150),49,150,2,1,1e-1,100.)
  print*,''
  print*,''
  call asciiplot(sfcRain(:,151:300),49,150,2,1,1e-1,100.)
  print*, maxval(dPRRet%log10dNw)
  print*, minval(dPRRet%log10dNw)
  !stop
end subroutine radarRetSub3



subroutine dealloc_struct(i)
  use local_RD_var
  use geophysEns
  print*,i
  IF (ALLOCATED(emissoutL)) deallocate(emissoutL)
  if(allocated(hFreqPRg)) deallocate(hFreqPRg)

  IF (ALLOCATED(Yobs)) deallocate(Yobs)
  IF (ALLOCATED(Xup)) deallocate(Xup)
  IF (ALLOCATED(randemiss)) deallocate(randemiss)
  IF (ALLOCATED(Xens)) deallocate(Xens)
  IF (ALLOCATED(Yens)) deallocate(Yens)

  call deallocGeophys()
  call deallocateStormStructData(stormStruct)
  call deallocateDPRProfRet(radarRet)
  call deallocateDPRProfData(radarData)
  
  IF (ALLOCATED(ndn)) deallocate(ndn) 
  IF (ALLOCATED(xscalev)) deallocate(xscalev) 
  IF (ALLOCATED(logdnwf)) deallocate(logdnwf) 
  IF (ALLOCATED(ndnp)) deallocate(ndnp) 
  IF (ALLOCATED(rhPCij)) deallocate(rhPCij)
  IF (ALLOCATED(cldwPCij)) deallocate(cldwPCij)

end subroutine dealloc_struct

subroutine dealloc_chunk(i)
  use globalData
  integer :: i
  print*, i
  call deallocateDPRRetSpace(dPRRet)
  call deallocateHRescGMI(gmiData,gmi2Grid)
end subroutine dealloc_chunk
