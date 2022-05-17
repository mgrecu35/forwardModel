subroutine kgainf(dtb,n,s,kgain)
  integer :: n
  real :: dtb(n), s
  real, intent(out) :: kgain(n)
  call kgainc(dtb, n, s, kgain)
end subroutine kgainf

module oetbs
  real oe_tbs(300,49,13,4)
  real scattProp(13,300,49,88,3)
end module oetbs

subroutine get_oetbs(tbs,scattPropOut)
  use oetbs
  real :: tbs(300,49,13,4)
  real :: scattPropOut(13,300,49,88,3)
  tbs=oe_tbs
  scattPropOut=scattProp
end subroutine get_oetbs
subroutine getscattH(binNodes,pRate,z13,emiss,envNode,pType,&
     qv3D,press3D,airTemp3D,nw3d,nscans,npixs,nlev,nchans,idir,&
     i,j,nfreq,df,kextH,salbH,asymH,dsrate,clutFree,swc3d)
  implicit none
  integer :: nscans,npixs,nlev, nchans
  integer :: nfreq, idir
  real :: pRate(nscans,npixs,nlev)
  real :: z13(nscans,npixs,nlev),emiss(nscans,npixs,nchans)
  integer :: clutFree(nscans,npixs)
  real :: nw3d(nscans,npixs,nlev), press3d(nscans,npixs,nlev), &
       airTemp3d(nscans,npixs,nlev),qv3d(nscans,npixs,nlev),swc3d(nscans,npixs,nlev)
  integer :: binNodes(nscans,npixs,5)
  integer :: envNode(nscans,npixs,10)
  integer :: pType(nscans,npixs)
  integer :: i, j, k, ni, i1
  real :: kext(nfreq), salb(nfreq), asym(nfreq)
  real :: kext2(nfreq), salb2(nfreq), asym2(nfreq)
  real :: srate, f, swc, df, n0w1
  real :: kextH(88,nfreq),salbH(88,nfreq),asymH(88,nfreq)
  real :: kexttot(88), salbtot(88), asymtot(88)
  real :: sfcTemp(nscans,npixs), cldw(nscans,npixs,nlev)
  real :: freqs(13)
  integer :: iFreq(13), ifreqit, ireturn
  real :: absair,abswv,dsrate
  integer :: maxlyr, nlyr, ik, ifreqsit
  real :: tbout, umu, fisot, z_clw
  integer ::  sfcBin(nscans,npixs)
  real dx, dx0, swc2

  dx=0.25
  dx0=5.
  kextH=0
  salbH=0
  asymH=0
  kext2=0
  salb2=0
  asym2=0
  do k=11, min(binNodes(i,j,3)-1,clutFree(i,j))
     i1=i+int((88-k)*dx/dx0+0.5)*idir
     if(i1<1) i1=1
     if(i1>nscans) i1=nscans
     if (pType(i1,j)>0 .and. z13(i1,j,k+1)>10 ) then
        !if(i1>nscans) i1=nscans
        !if(i1<1) i1=1
        f= (k-binNodes(i1,j,1))/ &
             (1e-3+binNodes(i1,j,3)-binNodes(i1,j,1))
        if(f>1) f=1
        if(f<0) f=0.
        !f=0.0
        n0w1=nw3d(i1,j,k+1)
        if(n0w1>3.5) n0w1=3.5
        if(n0w1<-3.5) n0w1=-3.5
10      continue
        call getsnowp_z2(n0w1,z13(i1,j,k+1),kext,salb,asym,srate,swc,nfreq)
        if(dsrate>1e-5) then
           call getsnowp(n0w1,dsrate,swc2,kext2,salb2,asym2,nfreq)
           !print*, dsrate, swc2, kext2
        else
           swc2=0
           kext2=0
           salb2=0
           asym2=0
        end if
        if( isnan(kext(1)) .or. maxval(kext)>1.9e5) then
           print*, z13(i1,j,k+1), nw3d(i1,j,k+1), n0w1,binNodes(i1,j,:), f, df
           print*, i, j, i1, k
           print*, kext
           print*, salb
           print*, asym
           !stop
           if(n0w1>1) then 
              n0w1=1
           else
              n0w1=n0w1/2
           endif
           goto 10
           
        endif
     else
        kext=0
        salb=0
        asym=0
        kext2=0
        salb2=0
        asym2=0
     endif
   
     kextH(k+1,:)=kext+kext2
     salbH(k+1,:)=(salb*kext+salb2*kext2)/(kext+kext2+1e-10)
     asymH(k+1,:)=(asym*salb*kext+asym2*salb2*kext2)*&
          (salb*kext+salb2*kext2)/((salb*kext+salb2*kext2)**2+1e-10)

     if (isnan(salbH(k+1,1))) then
        print*, kext2, kext
        print*, salb2
        print*, n0w1
        print*, k+1
        !STOP
     endif
  enddo
  do k=binNodes(i,j,3), binNodes(i,j,5)
     i1=i+int((88-k)*dx/dx0+0.5)*idir
     if(i1>nscans) i1=nscans
     if(i1<1) i1=1
     if (pType(i1,j)>0) then
        if(nw3d(i1,j,k+1)>3.5) nw3d(i1,j,k+1)=3.5
        if(nw3d(i1,j,k+1)<-3.5) nw3d(i1,j,k+1)=-3.5
        call getrainp(nw3d(i1,j,k+1),pRate(i1,j,k+1),kext,salb,asym,nfreq)
        if(k<binNodes(i,j,4)) then
           f=(k-binNodes(i,j,3))/(binNodes(i,j,4)-binNodes(i,j,3)+1e-3)
           n0w1=nw3d(i1,j,k+1)
           if(n0w1>3.5) n0w1=3.5
           if(n0w1<-3.5) n0w1=-3.5
           if(f*dsrate>1e-5) then
              call getsnowp(n0w1,(1-f)*dsrate,swc2,kext2,salb2,asym2,nfreq)
           else
              kext2=0
              salb2=0
              asym2=0
           endif
        else
           swc2=0
           kext2=0
           salb2=0
           asym2=0
        end if
     else
        kext=0
        salb=0
        asym=0
     endif
    
     kextH(k+1,:)=kext+kext2
     salbH(k+1,:)=(salb*kext+salb2*kext2)/(kext+kext2+1e-10)
     asymH(k+1,:)=(asym*salb*kext+asym2*salb2*kext2)*&
          (salb*kext+salb2*kext2)/((salb*kext+salb2*kext2)**2+1e-10)
     if(maxval(kextH(k+1,:))>197089.1-1) then
        print*, kext
        print*, kext2
        print*, nw3d(i1,j,k+1),pRate(i1,j,k+1)
        !stop
     end if
101 format(4(I4)," ",4(F9.5))
     if(j.eq.24) then
        !print*, idir
        !print*,i,i1,j,k, swc3d(i1,j,k+1)+swc2
        if(idir==1) then
           !write(*,101) i, i1,j,k+1, kextH(k+1,7),salbH(k+1,7),asymH(k+1,7),swc3d(i1,j,k+1)+swc2
        endif
     endif
     if (isnan(salbH(k+1,1))) then
        print*, 'rain'
        print*, k+1
        print*, kext2, kext
        print*, salb2
        print*, n0w1
        !stop
     endif
  enddo
  do k=binNodes(i,j,5)+1,88
     kextH(k,:)=kext
     salbH(k,:)=salb
     asymH(k,:)=asym
  enddo
end subroutine getscattH

subroutine updatepfields(binNodes,pRate,swc3d,z13,emiss,envNode,pType,&
     qv3D,press3D,airTemp3D,nw3d,nscans,npixs,nlev,nchans,idir,&
     i,j,nfreq,df,dsrate,clutFree)
  implicit none
  integer :: nscans,npixs,nlev, nchans
  integer :: nfreq, idir
  integer :: clutFree(nscans,npixs)
  real :: pRate(nscans,npixs,nlev),swc3d(nscans,npixs,nlev)
  real :: z13(nscans,npixs,nlev),emiss(nscans,npixs,nchans)
  real :: nw3d(nscans,npixs,nlev), press3d(nscans,npixs,nlev), &
       airTemp3d(nscans,npixs,nlev),qv3d(nscans,npixs,nlev)
  integer :: binNodes(nscans,npixs,5)
  integer :: envNode(nscans,npixs,10)
  integer :: pType(nscans,npixs)
  integer :: i, j, k, ni, i1
  real :: kext(nfreq), salb(nfreq), asym(nfreq)
  real :: kext2(nfreq), salb2(nfreq), asym2(nfreq)
  real :: srate, f, swc, df, n0w1
  real :: kextH(88,nfreq),salbH(88,nfreq),asymH(88,nfreq)
  real :: kexttot(88), salbtot(88), asymtot(88)
  real :: sfcTemp(nscans,npixs), cldw(nscans,npixs,nlev)
  real :: freqs(13)
  integer :: iFreq(13), ifreqit, ireturn
  real :: absair,abswv
  integer :: maxlyr, nlyr, ik, ifreqsit
  real :: tbout, umu, fisot, z_clw
  integer ::  sfcBin(nscans,npixs)
  real dx, dx0
  real :: dsrate, swc2
  dx=0.25
  dx0=5.
  kextH=0
  salbH=0
  asymH=0
  do k=11, min(binNodes(i,j,3)-1,clutFree(i,j))
     i1=i+int((88-k)*dx/dx0+0.5)*idir
     if(i1>nscans) i1=nscans
     if(i1<1) i1=1
     if (pType(i1,j)>0 .and. z13(i1,j,k+1)>0 ) then
        f= (k-binNodes(i1,j,1))/ &
             (1e-3+binNodes(i1,j,3)-binNodes(i1,j,1))
        if(f>1) f=1
        if(f<0) f=0
        !f=0.0
        n0w1=nw3d(i1,j,k+1)
        if(n0w1>3.5) n0w1=3.5
        if(n0w1<-3.5) n0w1=-3.5
        call getsnowp_z2(n0w1,z13(i1,j,k+1),kext,salb,asym,srate,swc,nfreq)
        call getsnowp(n0w1,dsrate,swc2,kext2,salb2,asym2,nfreq)
        if(swc2>1e2) then
           print*, 'above'
           print*, dsrate, swc2, n0w1, z13(i1,j,k+1)
           !stop
        endif
        swc3d(i1,j,k+1)=swc+swc2
        pRate(i1,j,k+1)=srate+dsrate
        nw3d(i1,j,k+1)=n0w1
     endif
  enddo
  do k=binNodes(i,j,3), binNodes(i,j,5)
     i1=i+int((88-k)*dx/dx0+0.5)*idir
     if(i1>nscans) i1=nscans
     if(i1<1) i1=1
     if (pType(i1,j)>0) then
        call getrainp(nw3d(i1,j,k+1),pRate(i1,j,k+1),kext,salb,asym,nfreq)
        if(k<binNodes(i,j,4)) then
           n0w1=nw3d(i1,j,k+1)
           f=(k-binNodes(i,j,3))/(binNodes(i,j,4)-binNodes(i,j,3)+1e-3)
           call getsnowp(n0w1,(1-f)*dsrate,swc2,kext2,salb2,asym2,nfreq)
           if(swc2>1e2) then
              print*, 'below'
              print*, dsrate, swc2, n0w1, z13(i1,j,k+1)
              !stop

           endif
        else
           swc2=0
           f=1.
        end if
        swc3d(i1,j,k+1)=swc3D(i1,j,k+1)+swc2
        pRate(i1,j,k+1)=pRate(i1,j,k+1)+(1-f)*dsrate
     end if
  end do
end subroutine updatepfields
subroutine slant1d(binNodes,pRate,z13,emiss,envNode,pType,&
     qv3D,press3D,airTemp3D,nw3d,tbsim,nscans,npixs,nlev,nchans,idir,&
     sfcBin,sfcTemp,cldw,df,i,j,nfreq, tbout_all, f2, dsrate, clutFree, swc3d)
  use oetbs
  implicit none
  integer :: nscans,npixs,nlev, nchans
  integer :: nfreq, idir
  real :: pRate(nscans,npixs,nlev)
  integer :: clutFree(nscans,npixs)
  real :: z13(nscans,npixs,nlev),emiss(nscans,npixs,nchans)
  real :: nw3d(nscans,npixs,nlev), press3d(nscans,npixs,nlev), &
       airTemp3d(nscans,npixs,nlev),qv3d(nscans,npixs,nlev),swc3d(nscans,npixs,nlev)
  integer :: binNodes(nscans,npixs,5)
  integer :: envNode(nscans,npixs,10)
  real, intent(out) :: tbsim(nscans,npixs,nchans)
  integer :: pType(nscans,npixs)
  integer :: i, j, k, ni, i1
  integer :: rangei(nlev)
  real :: kext(nfreq), salb(nfreq), asym(nfreq)
  real :: srate, f, swc, df, n0w1
  real :: kextH_u(88,nfreq),salbH_u(88,nfreq),asymH_u(88,nfreq)
  real :: kextH_d(88,nfreq),salbH_d(88,nfreq),asymH_d(88,nfreq)
  real :: kexttot_d(88), salbtot_d(88), asymtot_d(88)
  real :: kexttot_u(88), salbtot_u(88), asymtot_u(88)
  real :: sfcTemp(nscans,npixs), cldw(nscans,npixs,nlev)
  real :: freqs(13)
  integer :: iFreq(13), ifreqit, ireturn
  real :: absair,abswv
  real :: asym1d_u(88), kext1d_u(88), salb1d_u(88), lyr_hgt(89), lyr_temp(89)
  real :: asym1d_d(88), kext1d_d(88), salb1d_d(88)
  integer :: maxlyr, nlyr, ik, ifreqsit
  real :: tbout, umu, fisot, z_clw
  integer ::  sfcBin(nscans,npixs)
  real :: tbout_all(nchans), f2, tbout2
  real :: dsrate
  freqs=(/10.6,10.6,18.7,18.7,23.,37.,37.,89.,89.,166.,166.,186.35,190.35/)
  iFreq=(/0,0,1,1,2,3,3,4,4,5,5,6,7/)
  !ifreqG=(/1,1,2,2,3,4,4,5,5,6,6,7,8/)
  call getscattH(binNodes,pRate,z13,emiss,envNode,pType,&
       qv3D,press3D,airTemp3D,nw3d,nscans,npixs,nlev,nchans,idir,&
       i,j,nfreq,df,kextH_u,salbH_u,asymH_u, dsrate, clutFree, swc3d)
  call getscattH(binNodes,pRate,z13,emiss,envNode,pType,&
       qv3D,press3D,airTemp3D,nw3d,nscans,npixs,nlev,nchans,-idir,&
       i,j,nfreq,df,kextH_d,salbH_d,asymH_d, dsrate, clutFree, swc3d)

  !print*, binNodes(i,j,:)
  !print*, z13(i,j,:)
  !print*, pRate(i,j,:)

  !print*, envNode(i,j,:)
  do ifreqit=1,13
     kexttot_d=0
     salbtot_d=0
     asymtot_d=0
     kexttot_u=0
     salbtot_u=0
     asymtot_u=0
    do k=envNode(i,j,1)+1,envNode(i,j,10)+1
       call GasabsR98(freqs(ifreqit),airTemp3d(i,j,k),&
            (1+0.1*f2)*qv3d(i,j,k)*1e-3,&
            press3d(i,j,k)*1e2,&
            absair,abswv,ireturn)
       !absair=0
       !call monortm_lut(ifreqit, press3d(i,j,k), airTemp3d(i,j,k),&
       !     (1+0.1*f2)*qv3d(i,j,k), abswv)

       if(absair+abswv>1e2) then
          print*, f2, qv3d(i,j,k)
          print*, absair+abswv
          !stop
       end if
       if(cldw(i,j,k)>0) then
          call gcloud(freqs(ifreqit),airTemp3d(i,j,k),cldw(i,j,k),z_clw)
       else
          z_clw=0
       end if
       kexttot_u(k)=absair+abswv+z_clw+kextH_u(k,iFreq(ifreqit)+1)
       if(kexttot_u(k)>197089.1-1) then
          print*, kextH_u(k,iFreq(ifreqit)+1), iFreq(ifreqit)+1, ifreqit
          !stop
       end if
       salbtot_u(k)=salbH_u(k,iFreq(ifreqit)+1)*&
            kextH_u(k,iFreq(ifreqit)+1)/kexttot_u(k)
       asymtot_u(k)=asymH_u(k,iFreq(ifreqit)+1)
       kexttot_d(k)=absair+abswv+z_clw+kextH_d(k,iFreq(ifreqit)+1)
       salbtot_d(k)=salbH_d(k,iFreq(ifreqit)+1)*&
            kextH_d(k,iFreq(ifreqit)+1)/kexttot_d(k)
       asymtot_d(k)=asymH_d(k,iFreq(ifreqit)+1)
       scattProp(iFreq(ifreqit)+1,i,j,k,1)=kexttot_u(k)
       scattProp(iFreq(ifreqit)+1,i,j,k,2)=salbtot_u(k)
       scattProp(iFreq(ifreqit)+1,i,j,k,3)=asymtot_u(k)
       if(isnan(salbtot_u(k))) then
          print*, absair, abswv, z_clw
          print*, kextH_u(k,iFreq(ifreqit)+1)
          print*, k
          print*, binNodes(i,j,:)
          !stop
       endif
    enddo
    ik=1
    lyr_hgt(1)=0
    lyr_temp(1)=sfcTemp(i,j)
    
    !print*, 'sfcBin=',sfcBin(i,j), 'nodes=',envNode(i,j,1) ,envNode(i,j,10)
    !print*,'stop' 
    !stop
    do k=sfcBin(i,j)+1,envNode(i,j,1)+1,-1
       kext1d_u(ik)=kexttot_u(k)
       salb1d_u(ik)=salbtot_u(k)
       asym1d_u(ik)=asymtot_u(k)
       kext1d_d(ik)=kexttot_d(k)
       salb1d_d(ik)=salbtot_d(k)
       asym1d_d(ik)=asymtot_d(k)       
       lyr_hgt(ik+1)=ik*0.25*cos(abs(j-13)/12.*9./180.*3.1415)
       lyr_temp(ik+1)=airTemp3d(i,j,k)
       ik=ik+1
    enddo
    nlyr=ik-1
    umu=cos(53/180.*3.1415)
    fisot=2.7
    !print*,'nlyr=',nlyr, binNodes(i,j,:), envNode(i,j,:)
    call radtran(umu, nlyr, tbout, sfcTemp(i,j), lyr_temp(1:nlyr+1), lyr_hgt(1:nlyr+1), &
         kext1d_u(1:nlyr), &
         salb1d_u(1:nlyr), asym1d_u(1:nlyr), &
         fisot, emiss(i,j,ifreqit), emiss(i,j,ifreqit), &
         nlyr)

    if (isnan(tbout)) then
       print*, nlyr, f2
       print*, kext1d_u(1:nlyr)
       print*, emiss(i,j,ifreqit), ifreqit
       print*, salb1d_u(1:nlyr)
       print*, asym1d_u(1:nlyr)
       print*, sfcTemp(i,j)
       print*, lyr_hgt(1:nlyr+1)
       print*, lyr_temp(1:nlyr+1)
       !stop '"tbout" is a NaN'

    endif
    !print*, tbout
    !stop
    tbsim(i,j,ifreqit)=tbout
    tbout_all(ifreqit)=tbout
  enddo
end subroutine slant1d

subroutine  frteprep(binNodes,pRate,swc3d,pRateOut,swcOut,nwOut,z13,emiss,envNode,pType,&
     qv3D,press3D,airTemp3D,nw3d,tbsim,nscans,npixs,nlev,nchans,idir,&
     sfcBin,sfcTemp,cldw,tbobs,nfreq,clutFree,nscans_act)
  use oetbs
  implicit none
  integer :: nscans,nscans_act,npixs,nlev, nchans
  integer :: nfreq, idir
  integer :: clutFree(nscans,npixs)
  real :: pRate(nscans,npixs,nlev), swc3d(nscans,npixs,nlev), tbobs(nscans,npixs,nchans)
  real :: z13(nscans,npixs,nlev),emiss(nscans,npixs,nchans)
  real :: nw3d(nscans,npixs,nlev), press3d(nscans,npixs,nlev), &
       airTemp3d(nscans,npixs,nlev),qv3d(nscans,npixs,nlev)
  integer :: binNodes(nscans,npixs,5)
  integer :: envNode(nscans,npixs,10)
  real, intent(out) :: tbsim(nscans,npixs,nchans)
  integer :: pType(nscans,npixs)
  integer :: i, j, k, ni, i1
  integer :: rangei(nlev)
  real :: kext(nfreq), salb(nfreq), asym(nfreq)
  real :: srate, f, swc, df, n0w1
  real :: kextH(88,nfreq),salbH(88,nfreq),asymH(88,nfreq)
  real :: sfcTemp(nscans,npixs), cldw(nscans,npixs,nlev)
  integer :: sfcBin(nscans,npixs), n6
  real :: tbout1(nchans), tbout2(nchans), dtb(6), kgain(6), s, e(6)
  real,intent(out) :: pRateOut(nscans,npixs,nlev), swcOut(nscans,npixs,nlev), nwOut(nscans,npixs,nlev)
  real :: f2, dsrate
  print*, nscans, npixs, nlev
  scattProp=-99.9
  !idir=0
  !stop
  !slant1d(binNodes,pRate,z13,emiss,envNode,pType,&
  !     qv3D,press3D,airTemp3D,nw3d,tbsim,nscans,npixs,nlev,nchans,idir,&
  !     sfcBin,sfcTemp,cldw,df,i,j,nfreq, tbout_all, f2)
  

  do i=1,88
     rangei(i)=i-1
  end do
  s=16.
  df=1.
  n6=6
  pRateOut=pRate
  swcOut=swc3d
  nwOut=nw3d
101 format(2(I6,3x),6(F7.2,1x),6(F7.2,1x))
  open(11,file='diagnost.txt',form='formatted',&
       status='unknown',action='write')
  oe_tbs=-99.
  do i=1,nscans_act
     do j=2,48!npixs
        kextH=0
        salbH=0
        asymH=0
        if (pType(i,j)>0) then
           df=0.
           f2=0.
           dsrate=0.
           !print*, nchans
           call slant1d(binNodes,pRate,z13,emiss,envNode,pType,&
                qv3D,press3D,airTemp3D,nw3d,tbsim,nscans,npixs,nlev,nchans,idir,&
                sfcBin,sfcTemp,cldw,df,i,j,nfreq, tbout1, f2, dsrate, clutFree, swc3d)
           df=1.
           dsrate=0.1
           oe_tbs(i,j,:,3)=tbout1(:)
           call slant1d(binNodes,pRate,z13,emiss,envNode,pType,&
                qv3D,press3D,airTemp3D,nw3d,tbsim,nscans,npixs,nlev,nchans,idir,&
                sfcBin,sfcTemp,cldw,df,i,j,nfreq, tbout2, f2, dsrate, clutFree, swc3d)
           
           do k=1,6
              dtb(k)=(tbout2(7+k)-tbout1(7+k))/dsrate
              
              if (tbobs(i,j,7+k)>90) then
                 e(k)=tbobs(i,j,7+k)-tbout1(7+k)
              else
                 e(k)=0
              endif
           enddo
           !print*, tbout1(7:13)
           !print*, tbout2(7:13)
           !print*, dtb
           call kgainc(dtb, n6, s, kgain)
           df=0
           do k=1,6
              dsrate=dsrate+kgain(k)*e(k)
           end do
           f2=1.
           if(dsrate<0.0001) dsrate=0.0001
           call slant1d(binNodes,pRate,z13,emiss,envNode,pType,&
                qv3D,press3D,airTemp3D,nw3d,tbsim,nscans,npixs,nlev,nchans,idir,&
                sfcBin,sfcTemp,cldw,df*0.,i,j,nfreq, tbout2, f2, dsrate, clutFree, swc3d)
           do k=1,6
              dtb(k)=(tbout2(7+k)-tbout1(7+k))/0.1
              e(k)=tbobs(i,j,7+k)-tbout1(7+k)
           enddo
           if(tbobs(i,j,10)<150) then
              !print*, tbout1(9:13)
           endif
           call kgainc(dtb, n6, s, kgain)
           f2=0.
           do k=1,6
              f2=f2+kgain(k)*e(k)
           end do
           if(f2<-2.5) f2=-2.5
           if(f2>2.5) f2=2.5
           if(df>2) df=2
           if(df<-2) df=-2
           !df=0.
           !print*, f2, df
           f2=0
           !df=1.
           df=0.0
           dsrate=dsrate*0.5
           call slant1d(binNodes,pRate,z13,emiss,envNode,pType,&
                qv3D,press3D,airTemp3D,nw3d,tbsim,nscans,npixs,nlev,nchans,idir,&
                sfcBin,sfcTemp,cldw,df,i,j,nfreq, tbout2, f2, dsrate, clutFree, swc3d)
           oe_tbs(i,j,:,1)=tbobs(i,j,:)
           oe_tbs(i,j,:,2)=tbout2(:)
202        format(I7,13(F7.2,1x))
           if(j.eq.24) then
              if(pType(i,j)>0) then
                 write(*,202) i,tbout2
              end if
              !stop
           end if
           !write(11,101) i,j,(tbout2(7+k),k=1,6),(tbobs(i,j,7+k),k=1,6)
           !dtb=array([array(tbSim)[-6:]-array(tbSim_noad)[-6:]])
           
           !pinv=linalg.pinv(dot(dtb.T,dtb)+eye(6)*4)
           !kgain=dot(dtb,pinv)
           !e=tb[i,j,7:13]-array(tbSim_noad)[-6:]
           !e[tb[i,j,7:13]<0]=0
           !df=dot(kgain,e)[0]
           call updatepfields(binNodes,pRateOut,swcOut,z13,emiss,envNode,pType,&
                qv3D,press3D,airTemp3D,nwOut,nscans,npixs,nlev,nchans,idir,&
                i,j,nfreq,df, dsrate, clutFree)
        end if
     enddo
  enddo
  close(11)
  !   pRateOut=pRate
  !   swcOut=swc3d
  !   nwOut=nw3d
end subroutine frteprep

subroutine  frteprep0(binNodes,pRate,swc3d,z13,emiss,envNode,pType,&
      qv3D,press3D,airTemp3D,nw3d,tbsim,nscans,npixs,nlev,nchans,idir,&
      sfcBin,sfcTemp,cldw,tbobs,nfreq,clutFree)
   implicit none
   integer :: nscans,npixs,nlev, nchans
   integer :: nfreq, idir
   integer :: clutFree(nscans,npixs)
   real :: pRate(nscans,npixs,nlev), swc3d(nscans,npixs,nlev), tbobs(nscans,npixs,nchans)
   real :: z13(nscans,npixs,nlev),emiss(nscans,npixs,nchans)
   real :: nw3d(nscans,npixs,nlev), press3d(nscans,npixs,nlev), &
        airTemp3d(nscans,npixs,nlev),qv3d(nscans,npixs,nlev)
   integer :: binNodes(nscans,npixs,5)
   integer :: envNode(nscans,npixs,10)
   real, intent(out) :: tbsim(nscans,npixs,nchans)
   integer :: pType(nscans,npixs)
   integer :: i, j, k, ni, i1
   integer :: rangei(nlev)
   real :: kext(nfreq), salb(nfreq), asym(nfreq)
   real :: srate, f, swc, df, n0w1
   real :: kextH(88,nfreq),salbH(88,nfreq),asymH(88,nfreq)
   real :: sfcTemp(nscans,npixs), cldw(nscans,npixs,nlev)
   integer :: sfcBin(nscans,npixs), n6
   real :: tbout1(nchans), tbout2(nchans), dtb(6), kgain(6), s, e(6)
   real :: f2
   !slant1d(binNodes,pRate,z13,emiss,envNode,pType,&
   !     qv3D,press3D,airTemp3D,nw3d,tbsim,nscans,npixs,nlev,nchans,idir,&
   !     sfcBin,sfcTemp,cldw,df,i,j,nfreq, tbout_all, f2)


   do i=1,88
      rangei(i)=i-1
   end do
   s=4.
   df=1.
   n6=6
   do i=1,nscans
      do j=1,npixs
         kextH=0
         salbH=0
         asymH=0
         if (pType(i,j)>0) then
            df=0.
            f2=0.
            call slant1d(binNodes,pRate,z13,emiss,envNode,pType,&
                 qv3D,press3D,airTemp3D,nw3d,tbsim,nscans,npixs,nlev,nchans,idir,&
                 sfcBin,sfcTemp,cldw,df,i,j,nfreq, tbout1, f2, clutFree)

            !dtb=array([array(tbSim)[-6:]-array(tbSim_noad)[-6:]])
            !df=0.
            !pinv=linalg.pinv(dot(dtb.T,dtb)+eye(6)*4)
            !kgain=dot(dtb,pinv)
            !e=tb[i,j,7:13]-array(tbSim_noad)[-6:]
            !e[tb[i,j,7:13]<0]=0
            !df=dot(kgain,e)[0]
         end if
      enddo
   enddo
   !   pRateOut=pRate
   !   swcOut=swc3d
   !   nwOut=nw3d
 end subroutine frteprep0
 
