!  SFM  04/06/2013  Added allocate_2AKuENV_space

!*******************************************************************************
!2345678 1 2345678 2 2345678 3 2345678 4 2345678 5 2345678 6 2345678 7 2345678 8
!*******************************************************************************
!
!---NAME:     allocate_2AKuENV_space
!
!---PURPOSE:  This subroutine allocates space needed to read into program
!             memory the values of interest to the combined algorithm 
!             processing that are to be extracted from a 2AKuENV file
!             Ku environment variables.  See f90DataTypes for explanation
!             of parameter contents.  The size parameters are also stored
!             into the structure holding the allocated array space.
!
!---Programmer:  S. McLaughlin   9/13/2012
!
!---COMMENTS: None.
!
!-------------------------------------------------------------------------------
SUBROUTINE allocate_2AKuENV_space (Lv2AKuENV_scan, nscan_2akuenv, nray2akuenv,  &
                                   nbin2akuenv, nwater, nwind)

    USE f90DataTypes

    IMPLICIT NONE

!...Calling Sequence Parameters
    TYPE (Lv2AKuENV_DataType),intent(out) :: Lv2AKuENV_scan  ! allocated structure

    INTEGER*4,intent(in) :: nscan_2akuenv  ! Number scans in granule
    INTEGER*4,intent(in) :: nray2akuenv    ! Number angle bins in each NS scan
    INTEGER*4,intent(in) :: nbin2akuenv    ! Number range bins in each NS and MS ray
    INTEGER*4,intent(in) :: nwater         ! TBD
    INTEGER*4,intent(in) :: nwind          ! Number of wind components: u,v

!---Begin Procedure

!!    PRINT *,'allocate_2AKuENV_space allocation params ',nscan_2akuenv,         &
!!             nray2akuenv, nbin2akuenv, nwater, nwind

    Lv2AKuENV_scan%nscan_2akuenv = nscan_2akuenv
    Lv2AKuENV_scan%nray2akuenv   = nray2akuenv
    Lv2AKuENV_scan%nbin2akuenv   = nbin2akuenv
    Lv2AKuENV_scan%nwater        = nwater
    Lv2AKuENV_scan%nwind         = nwind

    ALLOCATE (Lv2AKuENV_scan%secondofday(nscan_2akuenv))
    ALLOCATE (Lv2AKuENV_scan%dayofyear(nscan_2akuenv))
    ALLOCATE (Lv2AKuENV_scan%year(nscan_2akuenv))
    ALLOCATE (Lv2AKuENV_scan%millisecond(nscan_2akuenv))
    ALLOCATE (Lv2AKuENV_scan%month(nscan_2akuenv))
    ALLOCATE (Lv2AKuENV_scan%dayofmonth(nscan_2akuenv))
    ALLOCATE (Lv2AKuENV_scan%hour(nscan_2akuenv))
    ALLOCATE (Lv2AKuENV_scan%minute(nscan_2akuenv))
    ALLOCATE (Lv2AKuENV_scan%second(nscan_2akuenv))

    ALLOCATE (Lv2AKuENV_scan%latitude(nray2akuenv,nscan_2akuenv))
    ALLOCATE (Lv2AKuENV_scan%longitude(nray2akuenv,nscan_2akuenv))
    ALLOCATE (Lv2AKuENV_scan%airTemp(nbin2akuenv,nray2akuenv,nscan_2akuenv))
    ALLOCATE (Lv2AKuENV_scan%airPressure(nbin2akuenv,nray2akuenv,nscan_2akuenv))
    ALLOCATE (Lv2AKuENV_scan%waterVapor(nwater,nbin2akuenv,nray2akuenv,nscan_2akuenv))
    ALLOCATE (Lv2AKuENV_scan%cldLiqWat(nwater,nbin2akuenv,nray2akuenv,nscan_2akuenv))
    ALLOCATE (Lv2AKuENV_scan%surfPressure(nray2akuenv,nscan_2akuenv))
    ALLOCATE (Lv2AKuENV_scan%groundTemp(nray2akuenv,nscan_2akuenv))
!begin WSO 9/1/13 add skin temperature
    ALLOCATE (Lv2AKuENV_scan%skinTemperature(nray2akuenv,nscan_2akuenv))
!end   WSO 9/1/13
    ALLOCATE (Lv2AKuENV_scan%surfWind(nwind,nray2akuenv,nscan_2akuenv))

!---End of Subroutine

    RETURN

END SUBROUTINE allocate_2AKuENV_space

subroutine allocateDPRSpace(this, ngates, nrays, n)
  use f90DataTypes
  integer            :: ngates, nrays
  type (dPRDataType) :: this 
  integer            :: n
  this%nmax=n
  this%ngates=ngates
!begin  WSO 9/5/13 add composite sigmapia; make variable name changes
!begin  WSO 9/10/13 add surface SNR variables
!begin  WSO 9/28/13 added bad ray flag variable and rain flag variable
  allocate(this%zku1c21(ngates,nrays,n), this%snrRatioku(nrays,n), &
       this%srtPIAku(nrays,n), this%srtsigmaPIAku(nrays,n), &
       this%xlon(nrays,n),  this%xlat(nrays,n), this%badRayFlag(nrays, n), &
       this%rainFlagBad(nrays, n), this%freezH(nrays,n), this%piaHB(nrays,n))
  allocate(this%sclon(n), this%sclat(n)) !Sjm 8/5/16
  allocate(this%zka1c21(ngates,nrays,n), this%snrRatioka(nrays, n), &
       this%dsrtPIAku(nrays,n), this%dsrtPIAka(nrays,n), &
       this%dsrtsigmaPIAku(nrays,n), this%dsrtsigmaPIAka(nrays,n))
  allocate(this%srtrelPIAku(nrays,n), this%dsrtrelPIA(nrays,n))
  allocate(this%MSRelibFlag(nrays,n), this%NSRelibFlag(nrays,n))
!end    WSO 9/28/13
!end    WSO 9/10/13
!end    WSO 9/5/13
!begin  WSO 8/15/14 added ioquality variable that summarizes status of input/output parameters
!using digital indices
  allocate(this%ioqualityflagku(nrays, n))  !for NS (Ku+GMI) mode
  allocate(this%ioqualityflagdpr(nrays, n)) !for MS (Ku+Ka+GMI) model
!end    WSO 8/15/14
!begin SJM 7/9/14 add sigma_zero
  allocate(this%sigmaZeroKu(nrays,n))
  allocate(this%sigmaZeroKa(nrays,n))
!end SJM 7/9/14
  allocate(this%node(5,nrays,n), this%rainType(nrays,n), this%BBbin(nrays,n))
  allocate(this%surfaceZku(nrays,n))
  allocate(this%surfaceZka(nrays,n))
  allocate(this%ilandOcean(nrays,n))
  this%nrays=nrays
  allocate(this%ig(nrays,n))
  allocate(this%jg(nrays,n))
  allocate(this%envQv(ngates,nrays,n))
  allocate(this%envPress(ngates,nrays,n))
  allocate(this%envTemp(ngates,nrays,n))
  allocate(this%envCloud(ngates,nrays,n))
  allocate(this%envSfcWind(nrays,n))
  allocate(this%envSfcWindU(nrays,n))
  allocate(this%envSfcWindV(nrays,n))
  allocate(this%snowIceCover(nrays,n))
  allocate(this%seaIceConcentration(nrays,n))
!begin  WSO 9/1/13 add skin temperature
  allocate(this%envSknTemp(nrays,n))
!end    WSO 9/1/13
  allocate(this%envSfcTemp(nrays,n))
  allocate(this%envSfcPress(nrays,n))
  allocate(this%binRealSurface(nrays,n))
  allocate(this%localZenithAngle(nrays,n))
  allocate(this%elevation(nrays,n))
  allocate(this%cBEst(nrays,n))
!  SFM  begin  04/16/2014; for nodes revision - in lieu G.Merca
  allocate(this%secondOfDay(n))
!  SFM  end    04/16/2014
end subroutine allocateDPRSpace

subroutine allocatecGMISpace(this, nS1f, nS2f, nrays, n, nmfreq, nmemb)
  use f90DataTypes
  type (cgMIDataType) :: this 
  integer            :: n, nS1f, nS2f, nrays, nmfreq, nmemb
  this%nmax=n
  this%n1b11=n
  this%nS1f=nS1f
  this%ns2f=nS2f
  this%nrays=nrays
  !print*, nrays, n, nS1f, nS2f
  allocate(this%gmiS1(nS1f,nrays,n))
  allocate(this%emissS1(nS1f,nrays,n))
  allocate(this%emissS2(nS1f,nrays,n))
  allocate(this%gmiS2(nS2f,nrays,n))
  allocate(this%gmilon(n))
  allocate(this%gmilat(n))
  allocate(this%S1lon(nrays,n))
  allocate(this%S1lat(nrays,n))
  allocate(this%S2lon(nrays,n))
  allocate(this%S2lat(nrays,n))
  allocate(this%sst(nrays,n))
  allocate(this%tpw(nrays,n))
  allocate(this%lwp(nrays,n))
  allocate(this%sfc_wind(nrays,n))
  allocate(this%landSea(nrays,n))
  allocate(this%rainFlag(nrays,n))
  allocate(this%tb10hC(nrays,n))
  allocate(this%tb10vC(nrays,n))
  allocate(this%tb19hC(nrays,n))
  allocate(this%tb19vC(nrays,n))
  allocate(this%tb21vC(nrays,n))
  allocate(this%tb37hC(nrays,n))
  allocate(this%tb37vC(nrays,n))
  allocate(this%rainflagdpr(nrays,n))
!  SFM  begin  04/16/2014; for nodes revision - in lieu G.Merca
  allocate(this%secondOfDay(n))
  allocate(this%SCLon(n))
  allocate(this%SCLat(n))
  allocate(this%S1eia(nrays,n))
  allocate(this%S2eia(nrays,n))
!  SFM  end    04/16/2014
  !write(*,*) 'simTb size=', size(this%simTb)
end subroutine allocatecGMISpace

subroutine allocateGMISimTb(this,nmfreq,nmemb)
  use f90DataTypes
  type (gMIDataType) :: this 
  allocate(this%simTb(2*this%nrays,3*this%n1b11-2,2,nmfreq,nmemb)) 
  allocate(this%simTbMin(2*this%nrays,3*this%n1b11-2,2,nmfreq)) 
  allocate(this%simTbMax(2*this%nrays,3*this%n1b11-2,2,nmfreq)) 
  allocate(this%simTbMean(2*this%nrays,3*this%n1b11-2,2,nmfreq)) 
  allocate(this%simTbStd(2*this%nrays,3*this%n1b11-2,2,nmfreq)) 
  allocate(this%dPRSfcRain(2*this%nrays,3*this%n1b11-2)) 
  allocate(this%dPRSfcRainEns(2*this%nrays,3*this%n1b11-2,nmemb))
  allocate(this%tpw3(2*this%nrays,3*this%n1b11-2)) 
  allocate(this%tpwEns(2*this%nrays,3*this%n1b11-2,nmemb)) 
end subroutine allocateGMISimTb

subroutine allocateDPRRetSpace(this, nmemb, nmfreq, ngates, nrays, n)
  use f90DataTypes
  use cldclass
  implicit none
  type (dPRRetType) :: this 
  integer :: nmemb, nmfreq, ngates,n,nrays

  this%nmemb=nmemb
  this%ngates=ngates
  this%n1c21=n
  this%nmfreq=nmfreq
  allocate(this%z13c(nmemb,ngates,nrays,n))
  allocate(this%z35mod0(nmemb,ngates,nrays,n))
!  SFM  begin  06/16/2014; for M. Grecu, multiple scattering
  allocate(this%dz35ms(nmemb,ngates,nrays,n))
!  SFM  end    06/16/2014
  allocate(this%z35(nmemb,ngates,nrays,n))
  allocate(this%pwc(nmemb,ngates,nrays,n))
  allocate(this%rrate(nmemb,ngates,nrays,n))
  allocate(this%sfcRain(nrays,n))
  allocate(this%sfcRainEns(nrays,n,nmemb))
  allocate(this%sfcd0Ens(nrays,n,nmemb))
  allocate(this%sfcNwEns(nrays,n,nmemb))
  allocate(this%sfcWindEns(nrays,n,nmemb)) !SJM 12/4/2014
  allocate(this%d0(nmemb,ngates,nrays,n))
  allocate(this%log10dNw(nmemb,ngates,nrays,n))
  allocate(this%n9(9,nrays,n))
  allocate(this%tb(nrays,n,2,nmfreq,nmemb))
  allocate(this%emTb(nrays,n,2,nmfreq,nmemb))
  allocate(this%emis(nrays,n,2,nmfreq,nmemb)) !SJM 5/28/15
  allocate(this%cldwcoeff(nrays,n,10,nmemb))
  allocate(this%MS%pwc(nmemb,ngates,nrays,n))
  allocate(this%MS%rrate(nmemb,ngates,nrays,n))
  allocate(this%MS%sfcRain(nrays,n))
  allocate(this%MS%sfcRainEns(nrays,n,nmemb))
  allocate(this%MS%sfcd0Ens(nrays,n,nmemb))
  allocate(this%MS%sfcNwEns(nrays,n,nmemb))
  allocate(this%MS%sfcWindEns(nrays,n,nmemb)) !SJM 12/4/2014
  allocate(this%MS%d0(nmemb,ngates,nrays,n))
  allocate(this%MS%zkuEns(nmemb,ngates,nrays,n))
  allocate(this%MS%zkaEns(nmemb,ngates,nrays,n))
  allocate(this%MS%log10dNw(nmemb,ngates,nrays,n))
  allocate(this%MS%tb(nrays,n,2,nmfreq,nmemb))
  allocate(this%MS%emTb(nrays,n,2,nmfreq,nmemb))
  allocate(this%MS%emis(nrays,n,2,nmfreq,nmemb)) !SJM 5/28/15
  allocate(this%MS%convtb(nrays,n,2,nmfreq,nmemb))
  allocate(this%MS%pia13mod(nrays,n,nmemb))
  allocate(this%MS%pia35mod(nrays,n,nmemb))
  allocate(this%MS%simSigmaZeroKu(nrays,n,nmemb)) !SJM 12/9/2014
  allocate(this%MS%simSigmaZeroKa(nrays,n,nmemb)) !SJM 12/9/2014
  allocate(this%simSigmaZeroKu(nrays,n,nmemb)) !SJM 12/9/2014
  allocate(this%simSigmaZeroKa(nrays,n,nmemb)) !SJM 12/9/2014
!  allocate(this%tbMin(nrays,n,2,nmfreq))
!  allocate(this%tbMax(nrays,n,2,nmfreq))
!  allocate(this%tbMean(nrays,n,2,nmfreq))
!  allocate(this%tbStd(nrays,n,2,nmfreq))
  allocate(this%imu(nmemb,nrays,n))
  allocate(this%icc(nmemb,nrays,n))
  allocate(this%jcc(nmemb,nrays,n))
  allocate(this%sfc_wind(nmemb))
  allocate(this%tpwEns(nrays,n,nmemb))
!  allocate(this%tbCMean(nrays,n,2,nmfreq))
!  allocate(this%tbSim(nrays,n,2,nmfreq))
!  allocate(this%obsSig(nrays,n,9))
!  allocate(this%simSig(nrays,n,9))
  allocate(this%convtb(nrays,n,2,nmfreq,nmemb))
  allocate(this%sfcRainRef(nrays,n))
  allocate(this%pia13mod(nrays,n,nmemb))
  allocate(this%pia35mod(nrays,n,nmemb))

end subroutine allocateDPRRetSpace

subroutine allocateDPRProfRet(this,nmfreq,nmemb,ngates, nNodes)
  use f90Types
  use cldclass
  implicit none
  integer :: nmemb, ngates, nmfreq, nNodes
  type (radarRetType) :: this 
  !print*, nmemb, ngates
  this%nmemb=nmemb
  this%ngates=ngates
  this%nmfreq=nmfreq
  allocate(this%z13c(nmemb*ngates))
  allocate(this%z35mod0(nmemb*ngates))
!  SFM  begin  06/16/2014; for M. Grecu, multiple scattering
  allocate(this%dz35ms(nmemb*ngates))
!  SFM  end    06/16/2014
  allocate(this%z35(nmemb*ngates))
  allocate(this%pwc(nmemb*ngates))
  allocate(this%rrate(nmemb*ngates))
  allocate(this%d0(nmemb*ngates))
  allocate(this%log10dNw(nmemb*ngates))
  allocate(this%logdNw(nmemb*nNodes))
  allocate(this%tb(nmemb*2*nmfreq))
  allocate(this%emtb(nmemb*2*nmfreq))
  allocate(this%emis(nmemb*2*nmfreq))
  allocate(this%imu(nmemb))
  allocate(this%icc(nmemb))
  allocate(this%jcc(nmemb))
  allocate(this%sfc_wind(nmemb))
  allocate(this%sfc_windU(nmemb))
  allocate(this%sfc_windV(nmemb))
  allocate(this%pia13(nmemb))
  allocate(this%pia35(nmemb))
  allocate(this%simSigmaZeroKu(nmemb)) !SJM 12/3/2014
  allocate(this%simSigmaZeroKa(nmemb)) !SJM 12/3/2014
  allocate(this%z35mMean(ngates))
  allocate(this%z35mSig(ngates))
  allocate(this%tpwCldMod(nc))
end subroutine allocateDPRProfRet

subroutine allocateDPRProfData(this, ngates)
  use f90Types
  implicit none
  integer :: ngates
  type (radarDataType) :: this
  allocate(this%z13obs(ngates))
  allocate(this%z35obs(ngates))
  allocate(this%hh(ngates))
  this%ngates=ngates
end subroutine allocateDPRProfData

subroutine deallocateDPRProfData(this)
  use f90Types
  implicit none
  type (radarDataType) :: this
  deallocate(this%z13obs)
  deallocate(this%z35obs)
  deallocate(this%hh)
end subroutine deallocateDPRProfData

subroutine allocateStormStructData(this)
  use f90Types
  type (stormStructType) :: this
  allocate(this%nodes(5))
end subroutine allocateStormStructData
 
subroutine deallocateDPRRetSpace(this)
  use f90DataTypes
  use cldclass
  implicit none
  type (dPRRetType) :: this 
  integer :: nmemb, nmfreq, ngates,n,nrays

  deallocate(this%z13c)
  deallocate(this%z35mod0)
!  SFM  begin  06/16/2014; for M. Grecu, multiple scattering
  deallocate(this%dz35ms)
!  SFM  end    06/16/2014
  deallocate(this%z35)
  deallocate(this%pwc)
  deallocate(this%rrate)
  deallocate(this%sfcRain)
  deallocate(this%sfcRainEns)
  deallocate(this%sfcd0Ens)
  deallocate(this%sfcNwEns)
  deallocate(this%sfcWindEns)!SJM 12/4/2014
  deallocate(this%d0)
  deallocate(this%log10dNw)
  deallocate(this%n9)
  deallocate(this%tb)
  deallocate(this%emTb)
  deallocate(this%emis)
  deallocate(this%cldwcoeff)
 ! deallocate(this%tbMin)
 ! deallocate(this%tbMax)
 ! deallocate(this%tbMean)
 ! deallocate(this%tbStd)
  deallocate(this%imu)
  deallocate(this%icc)
  deallocate(this%jcc)
  deallocate(this%sfc_wind)
  deallocate(this%tpwEns)
! deallocate(this%tbCMean)
! deallocate(this%tbSim)
  deallocate(this%convtb)
  deallocate(this%sfcRainRef)
  deallocate(this%pia13mod)
  deallocate(this%pia35mod)
  deallocate(this%simSigmaZeroKu)
  deallocate(this%simSigmaZeroKa)
  deallocate(this%MS%pia13mod)
  deallocate(this%MS%pia35mod)
  deallocate(this%MS%pwc)
  deallocate(this%MS%rrate)
  deallocate(this%MS%sfcRain)
  deallocate(this%MS%sfcRainEns)
  deallocate(this%MS%sfcd0Ens)
  deallocate(this%MS%sfcNwEns)
  deallocate(this%MS%sfcWindEns) !SJM 12/4/2014
  deallocate(this%MS%simSigmaZeroKu)!SJM 12/9/2014
  deallocate(this%MS%simSigmaZeroKa)!SJM 12/9/2014
  deallocate(this%MS%d0)
  deallocate(this%MS%log10dNw)
  deallocate(this%MS%tb)
  deallocate(this%MS%emtb)
  deallocate(this%MS%emis)
  deallocate(this%MS%convtb)
  deallocate(this%MS%zkuEns)
  deallocate(this%MS%zkaEns)
end subroutine deallocateDPRRetSpace

subroutine deallocateDPRProfRet(this)
  use f90Types
  use cldclass
  implicit none
  integer :: nmemb, ngates, nmfreq, ntpw
  type (radarRetType) :: this 
  !print*, nmemb, ngates


  deallocate(this%z13c)
  deallocate(this%z35mod0)
!  SFM  begin  06/16/2014; for M. Grecu, multiple scattering
  deallocate(this%dz35ms)
!  SFM  end    06/16/2014
  deallocate(this%z35)
  deallocate(this%pwc)
  deallocate(this%rrate)
  deallocate(this%d0)
  deallocate(this%log10dNw)
  deallocate(this%logdNw)
  deallocate(this%tb)
  deallocate(this%emtb)
  deallocate(this%emis)
  deallocate(this%imu)
  deallocate(this%icc)
  deallocate(this%jcc)
  deallocate(this%sfc_wind,this%sfc_windU, this%sfc_windV)
  deallocate(this%pia13)
  deallocate(this%pia35)
  deallocate(this%simSigmaZeroKu) !SJM 12/3/2014
  deallocate(this%simSigmaZeroKa) !SJM 12/3/2014
  deallocate(this%z35mMean)
  deallocate(this%z35mSig)
  deallocate(this%tpwCldMod)
end subroutine deallocateDPRProfRet


subroutine deallocateStormStructData(this)
  use f90Types
  type (stormStructType) :: this
  deallocate(this%nodes)
end subroutine deallocateStormStructData


subroutine deallocateDPRSpace()
  use f90DataTypes
  integer            :: ngates, nrays
  type (dPRDataType) :: this 
  integer            :: n
 
!begin  WSO 9/5/13 deallocate with new variable names
!begin  WSO 9/10/13 deallocate surface SNR variables
!begin  WSO 9/28/13 deallocate bad ray flag and rain flag variable
  deallocate(this%zku1c21, this%snrRatioku, & 
       this%srtPIAku, this%srtsigmaPIAku, &
       this%xlon,  this%xlat, this%badRayFlag, this%rainFlagBad, this%freezH, this%piaHB)
  deallocate(this%sclon, this%sclat)
  deallocate(this%zka1c21, this%snrRatioka, & 
       this%dsrtPIAku, this%dsrtPIAka, this%dsrtsigmaPIAku, this%dsrtsigmaPIAka)
  deallocate(this%srtrelPIAku, this%dsrtrelPIA)
  deallocate(this%MSRelibFlag, this%NSRelibFlag)
  deallocate(this%sigmaZeroKu, this%sigmaZeroKa)
!end    WSO 9/28/13
!end    WSO 9/10/13
!end    WSO 9/5/13
!begin  WSO 8/15/14 added ioquality variable that summarizes status of
!input/output parameters
!using digital indices
  deallocate(this%ioqualityflagku, this%ioqualityflagdpr)
!end    WSO 8/15/14
  deallocate(this%node, this%rainType, this%BBbin)
  deallocate(this%surfaceZku)
  deallocate(this%surfaceZka)
  deallocate(this%ilandOcean)
  
  deallocate(this%ig)
  deallocate(this%jg)
  deallocate(this%envQv)
  deallocate(this%envPress)
  deallocate(this%envTemp)
  deallocate(this%envCloud)
  deallocate(this%envSfcWind)
  deallocate(this%envSfcWindU)
  deallocate(this%envSfcWindV)
  deallocate(this%snowIceCover)
  deallocate(this%seaIceConcentration)
  deallocate(this%envSfcTemp)
  deallocate(this%envSfcPress)
  deallocate(this%binRealSurface)
  deallocate(this%localZenithAngle)
  deallocate(this%elevation)
  deallocate(this%cBEst)
!  SFM  begin  04/16/2014; for nodes revision - in lieu G.Merca
  deallocate(this%secondOfDay)
!  SFM  end    04/16/2014
end subroutine deallocateDPRSpace

subroutine deallocateGMISpace()
  use f90DataTypes
  type (gMIDataType) :: this 
  integer            :: n, nlowf, nhf, nrays, nmfreq, nmemb
 
  if(associated(this%gmilow)) deallocate(this%gmilow)
  if(associated(this%gmihigh)) deallocate(this%gmihigh)
  if(associated(this%gmilon)) deallocate(this%gmilon)
  if(associated(this%gmilat)) deallocate(this%gmilat)
  if(associated(this%sScan)) deallocate(this%sScan)
  if(associated(this%eScan)) deallocate(this%eScan)
  if(associated(this%txlon)) deallocate(this%txlon)
  if(associated(this%txlat)) deallocate(this%txlat)
  if(associated(this%sfc_wind))  deallocate(this%sfc_wind)
  if(associated(this%sst)) deallocate(this%sst)
  if(associated(this%tpw)) deallocate(this%tpw)
  if(associated(this%tb10hC)) deallocate(this%tb10hC)
  if(associated(this%tb10vC))  deallocate(this%tb10vC)
  if(associated(this%tb19hC)) deallocate(this%tb19hC)
  if(associated(this%tb19vC)) deallocate(this%tb19vC)
  if(associated(this%tb21vC)) deallocate(this%tb21vC)
  if(associated(this%tb37hC))   deallocate(this%tb37hC)
  if(associated(this%tb37vC)) deallocate(this%tb37vC)
  if(associated(this%landSea)) deallocate(this%landSea)
  if(associated(this%rainFlag)) deallocate(this%rainFlag)
  if(associated(this%rainFlagDPR)) deallocate(this%rainFlagDPR)
  !write(*,*) 'simTb size=', size(this%simTb)
  
end subroutine deallocateGMISpace

subroutine deallocatecGMISpace(this)
  use f90DataTypes
  implicit none
  type (cgMIDataType) :: this 



  deallocate(this%gmiS1)
  deallocate(this%emissS1)
  deallocate(this%emissS2)
  deallocate(this%gmiS2)
  deallocate(this%gmilon)
  deallocate(this%gmilat)
  deallocate(this%S1lon)
  deallocate(this%S1lat)
  deallocate(this%S2lon)
  deallocate(this%S2lat)
  deallocate(this%sst)
  deallocate(this%tpw)
  deallocate(this%lwp)
  deallocate(this%sfc_wind)
  deallocate(this%landSea)
  deallocate(this%rainFlag)
  deallocate(this%tb10hC)
  deallocate(this%tb10vC)
  deallocate(this%tb19hC)
  deallocate(this%tb19vC)
  deallocate(this%tb21vC)
  deallocate(this%tb37hC)
  deallocate(this%tb37vC)
  deallocate(this%rainflagdpr)
!  SFM  begin  04/16/2014; for nodes revision - in lieu G.Merca
  deallocate(this%secondOfDay)
  deallocate(this%SCLon)
  deallocate(this%SCLat)
  deallocate(this%S1eia)
  deallocate(this%S2eia)
!  SFM  end    04/16/2014
  !write(*,*) 'simTb size=', size(this%simTb)
end subroutine deallocatecGMISpace
