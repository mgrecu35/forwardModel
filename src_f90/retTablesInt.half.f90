


subroutine integratecvHB(z13,z35,i1,i2,pia13,pia35,z35mod,pwc,n0w,&
     ngates,nfreq,node,dr,imu,dpia13,dpia35,kext,salb,asym,rrate,d0)
  use tables2
  implicit none
  integer :: nfreq
  integer :: ngates,node(5),imu(ngates),i1,i2
  real :: z35mod(ngates)
  real :: z13(ngates), z35(ngates),n0w(ngates), pwc(ngates)
  integer :: i, j, i0, j0
  real :: pia13, pia35, dpia13, dpia35, zb

  real :: kext(ngates,nfreq),salb(ngates,nfreq),asym(ngates,nfreq)
  real :: rrate(ngates), d0(ngates)
  real :: dr,f 
  integer :: igoto, iNoAd
  if(n0w(i1)<-3.0) then
     n0w(i1)=-3.0
  endif
  if(n0w(i1)>3.0) then
     n0w(i1)=3.0
  endif
 
  igoto=1
  iNoAd=0
10 continue  
  if(isnan(n0w(i1))) then
     n0w(i1)=-1
     iNoAd=1
  endif
!  write(*,*)  i1,i2,node
  if(n0w(i1)<-5) then
     dpia13=0
     dpia35=0
     kext(i1,1:nfreq)=0
     salb(i1,1:nfreq)=0
     asym(i1,1:nfreq)=0
     print*, 'here'
     return
  endif
  do i=max(i1,node(1)),min(i2,node(2)-1)
     dpia13=0.
     dpia35=0.
     i0=(z13(i)-10.*n0w(i)-zmin)/dzbin+1
     if(i0<1) i0=1
     if(i0>nbinS2) i0=nbinS2
     dpia13=att13TableS2(i0,imu(i))*10.**n0w(i)*dr
     dpia35=att35TableS2(i0,imu(i))*10.**n0w(i)*dr
     pwc(i)=pwc13TableS2(i0,imu(i))+n0w(i)
     rrate(i)=pr13TableS2(i0,imu(i))*10.**n0w(i)
     rrate(i)=pr13TableS2(i0,imu(i))*10.**n0w(i)
     d0(i)=d013TableS2(i0,imu(i))
     !print*, pwc(i)
     if(10.**pwc(i)>0.2 .and. z13(i)-10.*n0w(i)<zmin .and. igoto<3) then
        n0w(i)=n0w(i)-0.3
        igoto=igoto+1
        goto 10
     endif
     kext(i,1:nfreq)=kextTableS2(i0,1:nfreq,imu(i))*10.**n0w(i)
     salb(i,1:nfreq)=salbTableS2(i0,1:nfreq,imu(i))
     asym(i,1:nfreq)=asymTableS2(i0,1:nfreq,imu(i))
!  SFM  begin  06/16/2014; for M.Grecu, multiple scattering
     z35mod(i)=z35TableS2(i0,imu(i))+10*n0w(i)-(pia35+dpia35)
     z35(i)=z35TableS2(i0,imu(i))+10*n0w(i)
!  SFM  end    06/16/2014
     !pia13=pia13+dpia13
     !pia35=pia35+dpia35
  enddo
 
  zb=0
  do i=max(i1,node(2)),min(i2,node(4)-1)
     dpia13=0.
     dpia35=0.
     f=((i-node(2)+0.)/(node(4)-node(2)+0.))**2
     !f=1.
     i0=(z13(i)-10.*n0w(i)-zmin)/dzbin+1
     if(i0<1) i0=1
     if(i0>nbinS2) i0=nbinS2
     dpia13=((1-f)*att13TableS2(i0,imu(i))+&
          f*att13Table(i0,imu(i)))*10.**n0w(i)*dr
     dpia35=((1-f)*att35TableS2(i0,imu(i))+&
          f*att35Table(i0,imu(i)))*10.**n0w(i)*dr
!  S2FM  begin  06/16/2014; for M.Grecu, multiple scattering
     z35mod(i)=((1-f)*z35TableS2(i0,imu(i))+&
          f*(z35Table(i0,imu(i))))+10*n0w(i)-(pia35+dpia35)
!  S2FM  end  06/16/2014
     z35(i)=((1-f)*z35TableS2(i0,imu(i))+f*(z35Table(i0,imu(i))))+10*n0w(i)
     pwc(i)=((1-f)*pwc13TableS2(i0,imu(i))+f*pwc13Table(i0,imu(i)))+n0w(i)
     rrate(i)=pr13TableS2(i0,imu(i))*10.**n0w(i)*(1-f)+&
          pr13Table(i0,imu(i))*10.**n0w(i)*(f)
     d0(i)=d013TableS2(i0,imu(i))*(1-f)+f*d013Table(i0,imu(i))
     kext(i,1:nfreq)=(1-f)*kextTableS2(i0,1:nfreq,imu(i))*10.**n0w(i)+ &
          f*kextTable(i0,1:nfreq,imu(i))*10.**n0w(i)
     salb(i,1:nfreq)=(1-f)*kextTableS2(i0,1:nfreq,imu(i))*&
          10.**n0w(i)*salbTableS2(i0,1:nfreq,imu(i))+&
          f*kextTable(i0,1:nfreq,imu(i))*10.**n0w(i)*&
          salbTable(i0,1:nfreq,imu(i))
     asym(i,1:nfreq)=(1-f)*kextTableS2(i0,1:nfreq,imu(i))*&
          10.**n0w(i)*salbTableS2(i0,1:nfreq,imu(i))*&
          asymTableS2(i0,1:nfreq,imu(i))+ &
          f*kextTable(i0,1:nfreq,imu(i))*10.**n0w(i)*&
          salbTable(i0,1:nfreq,imu(i))*&
          asymTable(i0,1:nfreq,imu(i))
     
     asym(i,1:nfreq)=asym(i,1:nfreq)/salb(i,1:nfreq)
     salb(i,1:nfreq)=salb(i,1:nfreq)/kext(i,1:nfreq)
     salb(i,1:nfreq)=(1-f)*salbTableS2(i0,1:nfreq,imu(i))+&
          f*salbTable(i0,1:nfreq,imu(i))
     asym(i,1:nfreq)=(1-f)*&
          asymTableS2(i0,1:nfreq,imu(i))+ &
          f*asymTable(i0,1:nfreq,imu(i))
     !pia13=pia13+dpia13
     !pia35=pia35+dpia35
  enddo
  
  do i=max(i1,node(4)),min(i2,node(5))
     dpia13=0.
     dpia35=0.
     i0=(z13(i)-10.*n0w(i)-zmin)/dzbin+1
     if(i0<1) i0=1
     if(i0>nbins) i0=nbins
     dpia13=att13Table(i0,imu(i))*10.**n0w(i)*dr
     !pia13=pia13+dpia13
     dpia35=att35Table(i0,imu(i))*10.**n0w(i)*dr
     pwc(i)=pwc13Table(i0,imu(i))+n0w(i)
     d0(i)=d013Table(i0,imu(i))
     rrate(i)=pr13Table(i0,imu(i))*10.**n0w(i)
     kext(i,1:nfreq)=kextTable(i0,1:nfreq,imu(i))*10.**n0w(i)
     if(kext(i,nfreq)>60 .and. igoto<3 .and. iNoAd==0) then
        n0w(i)=n0w(i)+log(60/kext(i,nfreq))
        igoto=igoto+1
        goto 10
        !print*, rrate(i), kext(i,nfreq), 10.**n0w(i), pia13, i, node(5)
     endif
     salb(i,1:nfreq)=salbTable(i0,1:nfreq,imu(i))
     asym(i,1:nfreq)=asymTable(i0,1:nfreq,imu(i))
!  S2FM  begin  06/16/2014; for M.Grecu, multiple scattering
     z35mod(i)=z35Table(i0,imu(i))+10*n0w(i)-(pia35+dpia35)
!  S2FM  end    06/16/2014
     z35(i)=z35Table(i0,imu(i))+10*n0w(i)
     !print*, z35mod(i), z35(i), i1,i2,i
     !pia35=pia35+dpia35
  enddo
  if(rrate(i1)>=201 .and. igoto<3 .and. iNoAd==0) then
     !write(*,*) rrate(i1), n0w(i1), z13(i1), i0, pr13Table(i0,imu(i1))
     n0w(i1)=n0w(i1)+log10(200/rrate(i1))
     igoto=igoto+1
     goto 10
  endif
 
  dpia13=dpia13/dr
  dpia35=dpia35/dr
  pwc(i1)=10.**pwc(i1)
end subroutine integratecvHB


subroutine integratestHB(z13,z13obs,z35,i1,i2,pia13,pia35,z35mod,pwc,n0w,&
     ngates,nfreq,node,dr,imu, dpia13, dpia35,kext,salb,asym,rrate,d0,dz1,dz2) !Sept 17, 2015 MG begin 
  use tables2
  implicit none 
  integer :: nfreq
  real :: dr, f
  integer :: ngates,node(5),imu(ngates),i1,i2
  real :: z35mod(ngates), z13obs(ngates)
  real :: z13(ngates), z35(ngates),n0w(ngates), pwc(ngates)
  integer :: i, j, i0, j0
  real :: pia13, pia35, dpia13, dpia35, zb
  real :: kext(ngates,nfreq),salb(ngates,nfreq),asym(ngates,nfreq)
  real :: rrate(ngates),d0(ngates)
  real :: n0bb, dz1, dz2
  integer :: igoto, iNoAd
  igoto=1
  iNoAd=0

10 continue 
  
  if(n0w(i1)<-3.0) then
     n0w(i1)=-3.0
  endif
  if(n0w(i1)>3.0) then
     n0w(i1)=3.0
  endif

  if(isnan(n0w(i1))) then
     n0w(i1)=-1
     iNoAd=1
  endif

  if(n0w(i1)<-5) then
     dpia13=0
     dpia35=0
     kext(i1,1:nfreq)=0
     salb(i1,1:nfreq)=0
     asym(i1,1:nfreq)=0
     return
  endif
  do i=max(i1,node(1)),min(i2,node(2)-1)
        dpia13=0.
        dpia35=0.
        i0=int((z13(i)-10.*n0w(i)-zmin)/dzbin)+1
        if(i0<1) i0=1
        if(i0>nbinS2) i0=nbinS2
        dpia13=att13TableS2(i0,imu(i))*10.**n0w(i)*dr
        dpia35=att35TableS2(i0,imu(i))*10.**n0w(i)*dr
        pwc(i)=pwc13TableS2(i0,imu(i))+n0w(i)
        rrate(i)=pr13TableS2(i0,imu(i))*10.**n0w(i)
        d0(i)=d013TableS2(i0,imu(i))
     
        kext(i,1:nfreq)=kextTableS2(i0,1:nfreq,imu(i))*10.**n0w(i)
        salb(i,1:nfreq)=salbTableS2(i0,1:nfreq,imu(i))
        asym(i,1:nfreq)=asymTableS2(i0,1:nfreq,imu(i))
!  S2FM  begin  06/16/2014; for M.Grecu, multiple scattering
        z35mod(i)=z35TableS2(i0,imu(i))+10*n0w(i)-(pia35+dpia35)
        z35(i)=z35TableS2(i0,imu(i))+10*n0w(i)
!  S2FM  end    06/16/2014
  enddo
  zb=0
  n0bb=-0.
  do i=max(i1,node(2)),min(i2,node(3)-1)
     dpia13=0.
     dpia35=0.
     !f=1.
     f=(i-node(2)+0.)/(node(3)-node(2)+0.)
     f=1-(maxval(z13obs(node(2):node(4)))-z13obs(i))/(dz1+0.1)
     if(f<0) f=0
     if(f>1) f=1
     i0=(z13(i)-10.*n0w(i)-10.*n0bb-zmin)/dzbin+1
     if(i0<1) i0=1
     if(i0>nbinS2) i0=nbinS2
     dpia13=((1-f)*att13TableS2(i0,imu(i))+f*att13TableBB(i0,imu(i)))*10.**(n0w(i)+n0bb)*dr
     dpia35=((1-f)*att35TableS2(i0,imu(i))+f*att35TableBB(i0,imu(i)))*10.**(n0w(i)+n0bb)*dr
!  S2FM  begin  06/16/2014; for M.Grecu, multiple scattering
     z35mod(i)=((1-f)*z35TableS2(i0,imu(i))+f*(z35TableBB(i0,imu(i))))+10*(n0w(i)+n0bb)-(pia35+dpia35)
     z35(i)=((1-f)*z35TableS2(i0,imu(i))+f*(z35TableBB(i0,imu(i))))+10*(n0w(i)+n0bb)
!  S2FM  end    06/16/2014
     pwc(i)=((1-f)*pwc13TableS2(i0,imu(i))+f*pwc13TableBB(i0,imu(i)))+(n0w(i)+n0bb)
     rrate(i)=pr13TableS2(i0,imu(i))*10.**(n0w(i)+n0bb)*(1-f)+pr13TableBB(i0,imu(i))*10.**(n0w(i)+n0bb)*(f)
     d0(i)=d013TableS2(i0,imu(i))*(1-f)+f*d013TableBB(i0,imu(i))
     kext(i,1:nfreq)=(1-f)*kextTableS2(i0,1:nfreq,imu(i))*10.**(n0w(i)+n0bb)+ &
          f*kextTableBB(i0,1:nfreq,imu(i))*10.**(n0w(i)+n0bb)
     salb(i,1:nfreq)=(1-f)*kextTableS2(i0,1:nfreq,imu(i))*10.**(n0w(i)+n0bb)* &
          salbTableS2(i0,1:nfreq,imu(i))+&
             f*kextTableBB(i0,1:nfreq,imu(i))*10.**(n0w(i)+n0bb)*salbTableBB(i0,1:nfreq,imu(i))
     asym(i,1:nfreq)=(1-f)*kextTableS2(i0,1:nfreq,imu(i))*&
          10.**(n0w(i)+n0bb)*salbTableS2(i0,1:nfreq,imu(i))*&
          asymTableS2(i0,1:nfreq,imu(i))+ &
          f*kextTableBB(i0,1:nfreq,imu(i))*10.**(n0w(i)+n0bb)*salbTableBB(i0,1:nfreq,imu(i))*&
          asymTableBB(i0,1:nfreq,imu(i))
        
     asym(i,1:nfreq)=asym(i,1:nfreq)/salb(i,1:nfreq)
     salb(i,1:nfreq)=salb(i,1:nfreq)/kext(i,1:nfreq)

  
  enddo

  do i=max(i1,node(3)),min(i2,node(4)-1)
        dpia13=0.
        dpia35=0.
        f=(i-node(3)+0.)/(node(4)-node(3)+0.)
        f=(maxval(z13obs(node(2):node(4)))-z13obs(i))/(dz2+0.1)
        if(f<0) f=0
        if(f>1) f=1
        !f=1.
        i0=(z13(i)-10.*(n0w(i)+n0bb)-zmin)/dzbin+1
        if(i0<1) i0=1
        if(i0>i0max(imu(i))) i0=i0max(imu(i))
        dpia13=((1-f)*att13TableBB(i0,imu(i))+f*att13Table(i0,imu(i)))*10.**(n0w(i)+n0bb)*dr

        dpia35=((1-f)*att35TableBB(i0,imu(i))+f*att35Table(i0,imu(i)))*10.**(n0w(i)+n0bb)*dr
        kext(i,1:nfreq)=(1-f)*kextTableBB(i0,1:nfreq,imu(i))*10.**(n0w(i)+n0bb)+ &
             f*kextTable(i0,1:nfreq,imu(i))*10.**(n0w(i)+n0bb)
        salb(i,1:nfreq)=(1-f)*kextTableBB(i0,1:nfreq,imu(i))*&
             10.**(n0w(i)+n0bb)*salbTableBB(i0,1:nfreq,imu(i))+&
             f*kextTable(i0,1:nfreq,imu(i))*10.**(n0w(i)+n0bb)*salbTable(i0,1:nfreq,imu(i))
        asym(i,1:nfreq)=(1-f)*kextTableBB(i0,1:nfreq,imu(i))*&
             10.**(n0w(i)+n0bb)*salbTableBB(i0,1:nfreq,imu(i))*&
             asymTableBB(i0,1:nfreq,imu(i))+ &
             f*kextTable(i0,1:nfreq,imu(i))*10.**(n0w(i)+n0bb)*salbTable(i0,1:nfreq,imu(i))*&
             asymTable(i0,1:nfreq,imu(i))
        asym(i,1:nfreq)=asym(i,1:nfreq)/salb(i,1:nfreq)
        salb(i,1:nfreq)=salb(i,1:nfreq)/kext(i,1:nfreq)
!        salb(i,1:nfreq)=(1-f)*salbTableBB(i0,1:nfreq,imu(i))+&
!             f*salbTable(i0,1:nfreq,imu(i))
!        asym(i,1:nfreq)=(1-f)*&
!             asymTableBB(i0,1:nfreq,imu(i))+ &
!             f*asymTable(i0,1:nfreq,imu(i))
        pwc(i)=((1-f)*pwc13TableBB(i0,imu(i))+f*pwc13Table(i0,imu(i)))+(n0w(i)+n0bb)
        rrate(i)=pr13TableBB(i0,imu(i))*10.**(n0w(i)+n0bb)*(1-f)+pr13Table(i0,imu(i))*10.**(n0w(i)+n0bb)*(f)
        d0(i)=d013TableBB(i0,imu(i))*(1-f)+f*d013Table(i0,imu(i))
!  S2FM  begin  06/16/2014; for M.Grecu, multiple scattering
        z35mod(i)=((1-f)*z35TableBB(i0,imu(i))+f*z35Table(i0,imu(i)))+10*(n0w(i)+n0bb)-(pia35+dpia35)
        z35(i)=((1-f)*z35TableBB(i0,imu(i))+f*z35Table(i0,imu(i)))+10*(n0w(i)+n0bb)
!  S2FM  end    06/16/2014
  enddo

  do i=max(i1,node(4)),min(i2,node(5))
     dpia13=0.
     dpia35=0.
     i0=(z13(i)-10.*n0w(i)-zmin)/dzbin+1
     if(i0<1) i0=1
     if(i0>nbins) i0=nbins
     dpia13=att13Table(i0,imu(i))*10.**n0w(i)*dr
     
     dpia35=att35Table(i0,imu(i))*10.**n0w(i)*dr
     pwc(i)=pwc13Table(i0,imu(i))+n0w(i)
     rrate(i)=pr13Table(i0,imu(i))*10.**n0w(i)
     d0(i)=d013Table(i0,imu(i))
     kext(i,1:nfreq)=kextTable(i0,1:nfreq,imu(i))*10.**n0w(i)
     salb(i,1:nfreq)=salbTable(i0,1:nfreq,imu(i))
     asym(i,1:nfreq)=asymTable(i0,1:nfreq,imu(i))
!  SFM  begin  06/16/2014; for M.Grecu, multiple scattering
     z35mod(i)=z35Table(i0,imu(i))+10*n0w(i)-(pia35+dpia35)
     z35(i)=z35Table(i0,imu(i))+10*n0w(i)
!  SFM  end    06/16/2014
     if(kext(i,nfreq)>60 .and. igoto<3 .and. iNoAd==0) then
        n0w(i)=n0w(i)+log(60/kext(i,nfreq))
        igoto=igoto+1
        goto 10
        !print*, rrate(i), kext(i,nfreq), 10.**n0w(i), pia13, i, node(5)
     endif
  enddo

  dpia35=dpia35/dr
  dpia13=dpia13/dr
  pwc(i1)=10.**pwc(i1)
!  if(salb(i1,6)>kext(i1,6)) then
!     print*, 'here we go', z13(i1), i1, node
!     stop
!  endif
end subroutine integratestHB !Sept 17, 2015 MG end

subroutine integrateanvHB(z13,z35,i1,i2,pia13,pia35,z35mod,pwc,n0w,&
     ngates,nfreq,node,dr,imu, dpia13, dpia35,kext,salb,asym,rrate,d0)
  use tables2
  implicit none 
  integer :: nfreq
  real :: dr, f
  integer :: ngates,node(5),imu(ngates),i1,i2
  real :: z35mod(ngates)
  real :: z13(ngates), z35(ngates),n0w(ngates), pwc(ngates)
  integer :: i, j, i0, j0
  real :: pia13, pia35, dpia13, dpia35, zb
  real :: kext(ngates,nfreq),salb(ngates,nfreq),asym(ngates,nfreq)
  real :: rrate(ngates),d0(ngates)
  integer :: igoto
  igoto=1
10 continue 
  
  if(n0w(i1)<-3.0) then
     n0w(i1)=-3.0
  endif
  if(n0w(i1)>3.0) then
     n0w(i1)=3.0
  endif
  if(n0w(i1)<-5) then
     dpia13=0
     dpia35=0
     kext(i1,1:nfreq)=0
     salb(i1,1:nfreq)=0
     asym(i1,1:nfreq)=0
     return
  endif
  do i=i1,i2
     dpia13=0.
     dpia35=0.
     i0=int((z13(i)-10.*n0w(i)-zmin)/dzbin)+1
     if(i0<1) i0=1
     if(i0>i0max(imu(i))-1) i0=i0max(imu(i))-1
     dpia13=att13TableS(i0,imu(i))*10.**n0w(i)*dr
     dpia35=att35TableS(i0,imu(i))*10.**n0w(i)*dr
     pwc(i)=pwc13TableS(i0,imu(i))+n0w(i)
     rrate(i)=pr13TableS(i0,imu(i))*10.**n0w(i)
     d0(i)=d013TableS(i0,imu(i))
     
     kext(i,1:nfreq)=kextTableS(i0,1:nfreq,imu(i))*10.**n0w(i)
     salb(i,1:nfreq)=salbTableS(i0,1:nfreq,imu(i))
     asym(i,1:nfreq)=asymTableS(i0,1:nfreq,imu(i))
     !  SFM  begin  06/16/2014; for M.Grecu, multiple scattering
     z35mod(i)=z35TableS(i0,imu(i))+10*n0w(i)-(pia35+dpia35)
     z35(i)=z35TableS(i0,imu(i))+10*n0w(i)
     !  SFM  end    06/16/2014
  enddo
end subroutine integrateanvHB
