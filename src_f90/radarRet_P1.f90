subroutine sfc(i,j,scLatPr1,scLonPR1,wfractPix,elon,elat,nmfreq2,nmemb1,n_chan)
  use local_RD_var
  use globalData
  use f90DataTypes
  use f90Types
  use ran_mod
  use LUT_def !SJM 7/9/2015
  implicit none

  integer :: stype
  integer :: i, j, ii, jj, nmemb1, nmfreq2, n_chan
  real :: elon, elat,  tpw_ij
  real :: s0Ku, s0Ka, s0stdKu, s0stdKa, s0corr, ds0Ku, ds0Ka
  !real :: sigmaZeroVarKu(49,300), sigmaZeroVarKa(49,300), sigmaZeroCov(49,300)
  real :: wfract(5,5), wfractm, wfractsd
  real                    emissv(n_chan)
  real                    emissh(n_chan)
  real                    emissv_std(n_chan)
  real                    emissh_std(n_chan)
  real :: windPert(100), windPertU(100), windPertV(100), qvPert(100), wfractpix
  real :: scLatPR1(49,300),scLonPR1(49,300), emis1, ebar1
!  real :: S1eiaPR(49,300), S2eiaPR(49,300)
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
        call calc_relAz(scLonPR1(i,j), scLatPR1(i,j), elon, elat, radarRet%sfc_windU(k+1), radarRet%sfc_windV(k+1), relAz)
        !relAz=0.
        !print*, i,j,radarRet%sfc_wind(k+1), '1'
        
        do ii=1,5
           call intplte_emis(ii,0,dPRData%envSfcTemp(i,j),radarRet%sfc_wind(k+1),relAz,S1eiaPR(i,j),emis,ebar1)
           emissv(ii)=emis1
           call intplte_emis(ii,1,dPRData%envSfcTemp(i,j),radarRet%sfc_wind(k+1),relAz,S1eiaPR(i,j),emis,ebar1)
           emissh(ii)=emis1
        end do
        
        do ii=6,6
           call intplte_emis(ii,0,dPRData%envSfcTemp(i,j),radarRet%sfc_wind(k+1),relAz,S2eiaPR(i,j),emis1,ebar1)
           emissv(ii)=emis1
           call intplte_emis(ii,1,dPRData%envSfcTemp(i,j),radarRet%sfc_wind(k+1),relAz,S2eiaPR(i,j),emis1,ebar1)
           emissh(ii)=emis1
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
           call calc_relAz(scLonPR1(i,j), scLatPR1(i,j), elon, elat, radarRet%sfc_windU(k+1), radarRet%sfc_windV(k+1), relAz)
           !relAz=0.
           !print*, i,j,radarRet%sfc_wind(k+1),'2'
           do ii=1,5
              call intplte_emis(ii,0,dPRData%envSfcTemp(i,j),radarRet%sfc_wind(k+1),relAz,S1eiaPR(i,j),emis1,ebar1)
              emissv(ii)=emis1
              call intplte_emis(ii,1,dPRData%envSfcTemp(i,j),radarRet%sfc_wind(k+1),relAz,S1eiaPR(i,j),emis1,ebar1)
              emissh(ii)=emis1
           end do
           do ii=6,6
              call intplte_emis(ii,0,dPRData%envSfcTemp(i,j),radarRet%sfc_wind(k+1),relAz,S2eiaPR(i,j),emis1,ebar1)
              emissv(ii)=emis1
              call intplte_emis(ii,1,dPRData%envSfcTemp(i,j),radarRet%sfc_wind(k+1),relAz,S2eiaPR(i,j),emis1,ebar1)
              emissh(ii)=emis1
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
end subroutine sfc
subroutine radarRetSub1(nmu2,  nmfreq2,   icL, tbRgrid,               &
      dprrain,ichunk,orbNumb,ialg,idir, i, j)
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
 ! real :: w10(49,300), w10_out_NS(49,300), w10_out_MS(49,300), w10_min, w10_max, emis, relAz
 ! real :: w10_rms_NS(49,300), emis_rms_NS(49,300,13), w10_rms_MS(49,300), emis_rms_MS(49,300,13)
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
!  real :: sigmaZeroVarKu(49,300), sigmaZeroVarKa(49,300), sigmaZeroCov(49,300)
  !end SJM 7/25/2014
!begin WSO 8/8/13
  real :: gatelength
  real :: depthBB, depthML, depth
  real :: mu_mean(49, 300)
  real :: mu_meanMS(49, 300)
  !real :: scLatPR(49,300),scLonPR(49,300),wfmap(49,300), fpmap(49,300,15), fpmapN(49,300,15)
  !real :: S1eiaPR(49,300), S2eiaPR(49,300)
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
  !print*, orbNumb,ichunk
  !call openascii(orbNumb,ichunk)
  
  !print*, dPRData%n1c21

  ifreqG(1:13)=(/1,1,2,2,3,4,4,5,5,6,6,7,8/)
  ipolG(1:13)=(/1,2,1,2,1,1,2,1,2,1,2,1,1/)

  
  
  eLon=dPRData%xlon(i,j)
  eLat=dPRData%xlat(i,j)
  call getwfraction(eLat,&
       eLon,wfmap(i,j))
  wfmap(i,j)=wfmap(i,j)/100.
  call setEnv(dPRData%envQv(:,i,j),dPRData%envTemp(:,i,j),&
       dPRData%envPress(:,i,j),dPRData%envSfcTemp(i,j),&
       dPRData%envSknTemp(i,j))
  dPRData%ioqualityflagku(i, j) = dPRData%ioqualityflagku(i,j) + 900000
  if(i > 12 .and. i < 38) then
     dPRData%ioqualityflagdpr(i, j) = dPRData%ioqualityflagdpr(i, j) + 900000
  endif
        
  nmemb1=nmemb

  if(dPRData%xlon(i,j)>-998) then  !4/15/14 MG begin
     iLandSea=-2
     iGMI=-99  !4/22/14 MG
     jGMI=-99  !4/22/14 MG
     ig = -99
     jg = -99
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
     endif
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
        endif
     endif
  end if
  
  if(dPRData%rainType(i,j)>=100 ) then
     dPRData%node(5,i,j)=dPRData%node(5,i,j)
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
     
     stormStruct%rainType=dPRData%rainType(i,j)
     stormStruct%rainType=stormStruct%rainType/100
     itop=1
     radarRet%sfc_wind(1)=dPRData%envSfcWind(i,j)
     radarRet%sfc_windU(1)=dPRData%envSfcWindU(i,j)
     radarRet%sfc_windV(1)=dPRData%envSfcWindV(i,j)
     w10(i,j)=radarRet%sfc_wind(1)
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
           call sfc(i,j,scLatPr,scLonPR,S1eiaPR,S2eiaPR,wfractPix,elon,elat,nmfreq2,nmemb1,n_chan)
           reliabFlag=dPRData%NSRelibFlag(i,j)
           !print*,dPRData%node(:,i,j)
           call  ensRadRetStCvKu(radarData,stormStruct,                &
                retParam, nmu2,radarRet, itop, rms1, rms2, sysdN, iit, &
                xscalev, randemiss, dPRData%localZenithAngle(i,j), &
                wfractPix, ichunk, i, j, dZms(i,j), msFlag(i, j)) 
           
           !begin WSO 2/11/17 assign dZms to output variable
           multiscatsurface_MS(i, j) = dZms(i, j)
           multiscatcalc_MS(i, j) = msFlag(i, j)
           
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
              dPRRet%log10dNw (k+1+(ibatch-1)*nmemb1,:,i,j)=-99
              dPRRet%d0 (k+1+(ibatch-1)*nmemb1,:,i,j)=-99
           end do
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
     
  endif
  
  

end subroutine radarRetSub1
