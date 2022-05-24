subroutine reflectivity(rwc1,swc1,wv1,dn1,temp,press,dr,nz,zka_m,&
     zka_out,attka_out,dzka_m, piaka, dpiaka, kext, salb, asym, pRate)
  use tablep2
  use tables2
  implicit none
  real :: rwc1(nz),swc1(nz),wv1(nz),temp(nz),press(nz),dr
  real :: zkar_1,attkar_1,zkas_1,attkas_1
  real :: n1
  real :: dn1(nz)
  integer :: nz, k, ibin
  real :: graupCoeff(2)
  real,intent(out) :: zka_m(nz)
  real :: attka
  real,intent(out) :: piaka, dzka_m(nz,nz,2), dpiaka(nz,2)
  real :: zkar_11, attkar_11, zkas_11, attkas_11
  integer :: k1
  real :: drwc
  real, intent(out) :: zka_out(nz),attka_out(nz)
  integer :: ireturn
  real :: absair, abswv, piaKaAtm
  real,intent(out) :: kext(nz), salb(nz), asym(nz), pRate(nz)
  real :: salbr1,asymr1,salbs1,asyms1,kextr1,kexts1
  real :: kextt, salbt, asymt
  graupCoeff=(/13.63604457, 28.58466471/)

  piaka=0
  !print*,rwc1
  !print*,swc1
  !print*,temp
  !print*,press
  !print*,nz,dr
  !return
  dzka_m=0
  dpiaka=0
  piaKaAtm=0
  do k=nz,1,-1
     call GasabsR98(35.5,temp(k),wv1(k),press(k),absair,abswv,ireturn)
     !print*, absair+abswv
     piaKaAtm=piaKaAtm+(absair+abswv)*4.343*0.125*2
     kextt=(absair+abswv)
     salbt=0
     asymt=0
     pRate(k)=0
     if(rwc1(k)>1e-2) then
        call bisection2(rwc(1:289),289,rwc1(k),ibin)
        pRate(k)=pRate(k)+rainrate(ibin)
        zkar_1=zkar(ibin)
        attkar_1=attkar(ibin)
        salbr1 = salbTable(ibin,4,1)*n1
        asymr1 = asymTable(ibin,4,1)
        kextr1 = attkar_1/4.343
        kextt = kextt + kextr1
        salbt = salbt + kextr1*salbr1
        asymt = asymt + kextr1*salbr1*asymr1
        call bisection2(rwc(1:289),289,1.1*rwc1(k),ibin)
        if (ibin<289) then
           zkar_11=zkar(ibin+1)
           attkar_11=attkar(ibin+1)
           drwc=rwc(ibin+1)-rwc(ibin)
        else
           zkar_11=zkar(ibin-1)
           attkar_11=attkar(ibin-1)
           drwc=rwc(ibin-1)-rwc(ibin)
        endif
     else
        zkar_1=0
        attkar_1=0
        zkar_11=0
        attkar_11=0
        drwc=0.1
     endif
     
     if (swc1(k)>1e-2) then
        if (temp(k).lt.(273.15)) then
           n1=exp(-0.122*(temp(k)-273.15))
        else
           n1=1
        endif
        if(n1>50) n1=50
        call bisection2(gwc(1:253),253,swc1(k)/n1,ibin)
        pRate(k)=pRate(k)+snowrate(ibin)*n1
        zkas_1=zkag(ibin)+10*log10(n1)
        salbS1 = salbTableG(ibin,4,1)
        asymS1 = asymTableG(ibin,4,1)
        kextS1 = attkag(ibin)*n1/4.343
        kextt = kextt + kexts1
        salbt = salbt + kexts1*salbs1
        asymt = asymt + kexts1*salbs1*asyms1
        
        attkas_1=attkag(ibin)*n1
        call bisection2(gwc(1:253),253,(1.1*swc1(k))/n1,ibin)
        zkas_11=zkag(ibin)+10*log10(n1)
        attkas_11=attkag(ibin)*n1
        
        if(swc1(k)/n1 < gwc(1)) then
           zkas_1=log10(swc1(k)/n1)*graupCoeff(1)+ &
                graupCoeff(2)+10*log10(n1)
           zkas_11=log10(1.1*swc1(k)/n1)*graupCoeff(1)+ &
                graupCoeff(2)+10*log10(n1)
        endif
     else
        zkas_1=0
        attkas_1=0
        zkas_11=0
        attkas_11=0
     endif
     !goto 10
     piaka=piaka+(attkas_1+attkar_1)*dr
     zka_m(k)=10*log10(10**(0.1*zkas_1)+10**(0.1*zkar_1))-piaka
     zka_out(k)=10*log10(10**(0.1*zkas_1)+10**(0.1*zkar_1))
     if(swc1(k)>1e-2) then
        dzka_m(k,k,2)=dzka_m(k,k,2)+&
             (10*log10(10**(0.1*zkas_11)+10**(0.1*zkar_1))-&
             10*log10(10**(0.1*zkas_1)+10**(0.1*zkar_1)))/(0.1*swc1(k))-&
             (attkas_11-attkas_1)/(0.1*swc1(k))*dr
        dpiaka(k,2)=dpiaka(k,2)+(attkas_11-attkas_1)/(0.1*swc1(k))*2*dr
        !do k1=k-1,1,-1
           !dzka_m(k1,2)=dzka_m(k1,2)-(attkas_11-attkas_1)/(0.1*swc1(k))
        !enddo
     endif
     if(rwc1(k)>1e-2) then
        dzka_m(k,k,1)=dzka_m(k,k,1)+&
             (10*log10(10**(0.1*zkar_11)+10**(0.1*zkas_1))-&
             10*log10(10**(0.1*zkas_1)+10**(0.1*zkar_1)))/(0.1*rwc1(k))-&
             (attkar_11-attkar_1)/(drwc)*dr
        !print*, zkar_11,zkar_1,  dzka_m(k,k,1), rwc1(k), drwc
        piaka=piaka+(attkas_1+attkar_1)*dr
        dpiaka(k,1)=dpiaka(k,1)+(attkar_11-attkar_1)/(drwc)*2*dr
        !do k1=k-1,1,-1
           !dzka_m(k1,1)=dzka_m(k1,1)-(attkar_11-attkar_1)/(0.1*rwc1(k))
        !enddo
     endif
     attka_out(k)=attkar_1+attkas_1
     kext(k)=kextt
     salb(k)=salbt/kextt
     if (salbt>1e-8) then
        asym(k)=asymt/salbt
     else
        asym(k)=0
     end if
     !10 continue
  end do
  !  print*, piaKaAtm
  !stop
end subroutine reflectivity


subroutine reflectivity_ka(rwc1,swc1,wv1,dn1,temp,press,dr,sl,nz,zka_m,&
     zka_out,attka_out, piaka, kext, salb, asym, kext1, salb1, asym1, pRate)
  use tablep2
  use tables2
  implicit none
  real :: rwc1(nz),swc1(nz),wv1(nz),temp(nz),press(nz),dr
  real :: zkar_1,attkar_1,zkas_1,attkas_1
  real :: n1, sl
  real :: dn1(nz)
  integer :: nz, k, ibin
  real :: graupCoeff(2)
  real,intent(out) :: zka_m(nz)
  real :: attka
  real,intent(out) :: piaka
  real :: zkar_11, attkar_11, zkas_11, attkas_11
  integer :: k1
  real :: drwc
  real, intent(out) :: zka_out(nz),attka_out(nz)
  integer :: ireturn
  real :: absair, abswv, piaKaAtm
  real,intent(out) :: kext(nz), salb(nz), asym(nz), pRate(nz)
  real,intent(out) :: kext1(nz), salb1(nz), asym1(nz)
  real :: salbr1,asymr1,salbs1,asyms1,kextr1,kexts1
  real :: salbr11,asymr11,kextr11
  real :: kextt, salbt, asymt
  real :: kextt1, salbt1, asymt1
  integer :: ibin2
  graupCoeff=(/13.63604457, 28.58466471/)

  piaka=0
  !print*,rwc1
  !print*,swc1
  !print*,temp
  !print*,press
  !print*,nz,dr
  !return
  
  piaKaAtm=0
  do k=nz,1,-1
     call GasabsR98(35.5,temp(k),wv1(k),press(k),absair,abswv,ireturn)
     !print*, absair+abswv
     piaKaAtm=piaKaAtm+(absair+abswv)*4.343*0.125*2
     kextt=(absair+abswv)
     kextt1=(absair+abswv)
     salbt=0
     asymt=0
     salbt1=0
     asymt1=0
     pRate(k)=0
     if(rwc1(k)>1e-2) then
        call bisection2(rwc(1:289),289,rwc1(k)/10**dn1(k),ibin)
        pRate(k)=pRate(k)+rainrate(ibin)*10**dn1(k)
        zkar_1=zkar(ibin)+10*dn1(k)
        attkar_1=attkar(ibin)*10**dn1(k)
        salbr1 = salbTable(ibin,4,1)
        asymr1 = asymTable(ibin,4,1)
        kextr1 = attkar_1/4.343
        kextt = kextt + kextr1
        salbt = salbt + kextr1*salbr1
        asymt = asymt + kextr1*salbr1*asymr1
     else
        zkar_1=0
        attkar_1=0
     endif
     
     if (swc1(k)>1e-2) then
        if (temp(k).lt.(283.15)) then
           n1=exp(-0.122*(temp(k)-273.15)*sl)
        else
           n1=exp(-0.122*(10)*sl)
        endif
        if(n1>20) n1=20
        call bisection2(gwc(1:253),253,swc1(k)/n1/10**dn1(k),ibin)
        pRate(k)=pRate(k)+snowrate(ibin)*n1*10**dn1(k)
        zkas_1=zkag(ibin)+10*log10(n1)+10*dn1(k)
        salbS1 = salbTableG(ibin,4,1)
        asymS1 = asymTableG(ibin,4,1)
        kextS1 = attkag(ibin)*n1/4.343*10**dn1(k)
        kextt = kextt + kexts1
        salbt = salbt + kexts1*salbs1
        asymt = asymt + kexts1*salbs1*asyms1
        kextt1 = kextt1 + kexts1
        salbt1 = salbt1 + kexts1*salbs1
        asymt1 = asymt1 + kexts1*salbs1*asyms1

        attkas_1=attkag(ibin)*n1*10**dn1(k)
        
        if(swc1(k)/n1/10**dn1(k) < gwc(1)) then
           zkas_1=log10(swc1(k)/n1)*graupCoeff(1)+ &
                graupCoeff(2)+10*log10(n1)+10*dn1(k)
        endif
     else
        zkas_1=0
        attkas_1=0
     endif

     piaka=piaka+(attkas_1+attkar_1)*dr
     zka_m(k)=10*log10(10**(0.1*zkas_1)+10**(0.1*zkar_1))-piaka
     piaka=piaka+(attkas_1+attkar_1)*dr
     zka_out(k)=10*log10(10**(0.1*zkas_1)+10**(0.1*zkar_1))

     attka_out(k)=attkar_1+attkas_1
     kext(k)=kextt
     salb(k)=salbt/kextt
     if (salbt>1e-8) then
        asym(k)=asymt/salbt
     else
        asym(k)=0
     end if
     
     kext1(k)=kextt1
     salb1(k)=salbt1/kextt1
     if (salbt1>1e-8) then
        asym1(k)=asymt1/salbt1
     else
        asym1(k)=0
     end if

     !10 continue
  end do
  !  print*, piaKaAtm
  !stop
end subroutine reflectivity_ka


subroutine reflectivity_ku(rwc1,swc1,wv1,dn1,temp,press,dr,sl,nz,zku_m,&
     zku_out,attku_out, piaku, kext, salb, asym, kext1, salb1, asym1, pRate,&
     dm_out)
  use tablep2
  use tables2
  implicit none
  real :: rwc1(nz),swc1(nz),wv1(nz),temp(nz),press(nz),dr
  real :: zkur_1,attkur_1,zkus_1,attkus_1
  real :: n1,sl
  real :: dn1(nz)
  integer :: nz, k, ibin
  real :: graupCoeff(2)
  real,intent(out) :: zku_m(nz)
  real :: attku
  real,intent(out) :: piaku, dm_out(nz)
  real :: zkur_11, attkur_11, zkus_11, attkus_11
  integer :: k1
  real :: drwc
  real, intent(out) :: zku_out(nz),attku_out(nz)
  integer :: ireturn
  real :: absair, abswv, piaKuAtm
  real,intent(out) :: kext(nz), salb(nz), asym(nz), pRate(nz)
  real,intent(out) :: kext1(nz), salb1(nz), asym1(nz)
  real :: salbr1,asymr1,salbs1,asyms1,kextr1,kexts1
  real :: salbr11,asymr11,kextr11
  real :: kextt, salbt, asymt
  real :: kextt1, salbt1, asymt1
  integer :: ibin2
  real :: dm_r, dm_s, dm_t
  graupCoeff=(/13.63604457, 28.58466471/)

  piaku=0
  !print*,rwc1
  !print*,swc1
  !print*,temp
  !print*,press
  !print*,nz,dr
  !return
  
  piaKuAtm=0
  do k=nz,1,-1
     call GasabsR98(35.5,temp(k),wv1(k),press(k),absair,abswv,ireturn)
     !print*, absair+abswv
     piaKuAtm=piaKuAtm+(absair+abswv)*4.343*0.125*2
     kextt=(absair+abswv)
     kextt1=(absair+abswv)
     salbt=0
     asymt=0
     salbt1=0
     asymt1=0
     pRate(k)=0
     if(rwc1(k)>1e-2) then
        call bisection2(rwc(1:289),289,rwc1(k)/10**dn1(k),ibin)
        pRate(k)=pRate(k)+rainrate(ibin)*10**dn1(k)
        zkur_1=zkur(ibin)+10*dn1(k)
        attkur_1=attkur(ibin)*10**dn1(k)
        salbr1 = 0.6*salbTable(ibin,1,1)+0.4*salbTable(ibin,2,1)
        asymr1 = 0.6*asymTable(ibin,1,1)+0.4*asymTable(ibin,2,1)
        kextr1 = attkur_1/4.343
        kextt = kextt + kextr1
        salbt = salbt + kextr1*salbr1
        asymt = asymt + kextr1*salbr1*asymr1
        dm_t=dmR(ibin)*rwc1(k)
     else
        zkur_1=0
        attkur_1=0
        dm_t=0.
     endif
     
     if (swc1(k)>1e-2) then
        if (temp(k).lt.(283.15)) then
           n1=exp(-0.122*(temp(k)-273.15)*sl)
        else
           n1=exp(-0.122*(10)*sl)
        endif
        if(n1>20) n1=20
        call bisection2(gwc(1:253),253,swc1(k)/n1/10**dn1(k),ibin)
        pRate(k)=pRate(k)+snowrate(ibin)*n1*10**dn1(k)
        zkus_1=zkug(ibin)+10*log10(n1)+10*dn1(k)
        salbS1 = salbTableG(ibin,4,1)
        asymS1 = asymTableG(ibin,4,1)
        salbS1 = 0.6*salbTableG(ibin,1,1)+0.4*salbTable(ibin,2,1)
        asymS1 = 0.6*asymTableG(ibin,1,1)+0.4*asymTable(ibin,2,1)
        kextS1 = attkug(ibin)*n1/4.343*10**dn1(k)
        kextt = kextt + kexts1
        salbt = salbt + kexts1*salbs1
        asymt = asymt + kexts1*salbs1*asyms1
        kextt1 = kextt1 + kexts1
        salbt1 = salbt1 + kexts1*salbs1
        asymt1 = asymt1 + kexts1*salbs1*asyms1

        attkus_1=attkug(ibin)*n1*10**dn1(k)
        !dm_s=dmG(imbin)
        dm_t=dm_t+dmR(ibin)*swc1(k)
        if(swc1(k)/n1 < gwc(1)) then
           zkus_1=log10(swc1(k)/n1)*graupCoeff(1)+ &
                graupCoeff(2)+10*log10(n1)+10*dn1(k)
        endif
     else
        zkus_1=0
        attkus_1=0
     endif
     dm_out(k)=dm_t/(rwc1(k)+swc1(k)+1e-5)
     piaku=piaku+(attkus_1+attkur_1)*dr
     zku_m(k)=10*log10(10**(0.1*zkus_1)+10**(0.1*zkur_1))-piaku
     piaku=piaku+(attkus_1+attkur_1)*dr
     zku_out(k)=10*log10(10**(0.1*zkus_1)+10**(0.1*zkur_1))

     attku_out(k)=attkur_1+attkus_1
     kext(k)=kextt
     salb(k)=salbt/kextt
     if (salbt>1e-8) then
        asym(k)=asymt/salbt
     else
        asym(k)=0
     end if
     
     kext1(k)=kextt1
     salb1(k)=salbt1/kextt1
     if (salbt1>1e-8) then
        asym1(k)=asymt1/salbt1
     else
        asym1(k)=0
     end if

     !10 continue
  end do
  !  print*, piaKuAtm
  !stop
end subroutine reflectivity_ku



