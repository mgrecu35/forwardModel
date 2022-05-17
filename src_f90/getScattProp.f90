subroutine getsnowp(n0w,srate,swc,kext,salb,asym,nfreq)
  use tables2
  real :: n0w,srate
  real, intent(out) :: kext(nfreq),salb(nfreq), asym(nfreq),swc
  integer :: imu
  imu=3
  call bisection2(pr13TableS2(:,imu),nbinS2,srate/10**n0w, i0)
  kext(1:nfreq)=kextTableS2(i0,1:nfreq,imu)*10.**n0w
  salb(1:nfreq)=salbTableS2(i0,1:nfreq,imu)
  asym(1:nfreq)=asymTableS2(i0,1:nfreq,imu)
  swc=10.**pwc13TableS2(i0,imu)*10.**n0w
end subroutine getsnowp

subroutine getrainp(n0w,rrate,kext,salb,asym,nfreq)
  use tables2
  real :: n0w,srate
  real, intent(out) :: kext(nfreq),salb(nfreq), asym(nfreq)
  integer :: imu
  imu=3
  call bisection2(pr13Table(:,imu),nbins,rrate/10**n0w, i0)
  kext(1:nfreq)=kextTable(i0,1:nfreq,imu)*10.**n0w
  salb(1:nfreq)=salbTable(i0,1:nfreq,imu)
  asym(1:nfreq)=asymTable(i0,1:nfreq,imu)
end subroutine getrainp


subroutine getsnowp2(n0w,iwc,kext,salb,asym,dpia13,dpia35,z13,z35,z94,dm,nfreq)
  use tables2
  real :: n0w,iwc
  real, intent(out) :: kext(nfreq),salb(nfreq), asym(nfreq)
  real, intent(out) :: dpia13,dpia35,z35,z13,z94,dm
  integer :: imu
  imu=3
  call bisection2(pwc13TableS2(:,imu),nbinS2,log10(iwc)-n0w, i0)
  kext(1:nfreq)=kextTableS2(i0,1:nfreq,imu)*10.**n0w
  salb(1:nfreq)=salbTableS2(i0,1:nfreq,imu)
  asym(1:nfreq)=asymTableS2(i0,1:nfreq,imu)
  dpia13=att13TableS2(i0,imu)*10.**n0w
  dpia35=att35TableS2(i0,imu)*10.**n0w
  z35=z35TableS2(i0,imu)+10*n0w
  !z94=z94TableS2(i0,imu)+10*n0w
  z13=zmin+i0*dzbin+10.*n0w
  dm=d013TableS2(i0,imu)
end subroutine getsnowp2

subroutine getsnowp_2(n0w,iwc,kext,salb,asym,nfreq)
  use tables2
  real :: n0w,iwc
  real, intent(out) :: kext(nfreq),salb(nfreq), asym(nfreq)
  integer :: imu
  imu=3
  call bisection2(pwc13TableS2(:,imu),nbinS2,log10(iwc)-n0w, i0)
  kext(1:nfreq)=kextTableS2(i0,1:nfreq,imu)*10.**n0w
  salb(1:nfreq)=salbTableS2(i0,1:nfreq,imu)
  asym(1:nfreq)=asymTableS2(i0,1:nfreq,imu)
end subroutine getsnowp_2


subroutine getrainp2(n0w,pwc,kext,salb,asym,dpia13,dpia35,z13,z35,z94,dm,nfreq)
  use tables2
  real :: n0w,pwc
  real, intent(out) :: kext(nfreq),salb(nfreq), asym(nfreq)
  real, intent(out) :: dpia13,dpia35,z13,z35,z94,dm
  integer :: imu
  imu=3
  call bisection2(pwc13Table(:,imu),nbins,log10(pwc)-n0w, i0)
  kext(1:nfreq)=kextTable(i0,1:nfreq,imu)*10.**n0w
  salb(1:nfreq)=salbTable(i0,1:nfreq,imu)
  asym(1:nfreq)=asymTable(i0,1:nfreq,imu)
  dpia13=att13Table(i0,imu)*10.**n0w
  dpia35=att35Table(i0,imu)*10.**n0w
  z35=z35Table(i0,imu)+10*n0w
  !z94=z94Table(i0,imu)+10*n0w
  z13=zmin+i0*dzbin+10.*n0w
  dm=d013Table(i0,imu)
end subroutine getrainp2


subroutine getrainp_2(n0w,pwc,kext,salb,asym,nfreq)
  use tables2
  real :: n0w,pwc
  real, intent(out) :: kext(nfreq),salb(nfreq), asym(nfreq)
  integer :: imu
  imu=3
  call bisection2(pwc13Table(:,imu),nbins,log10(pwc)-n0w, i0)
  kext(1:nfreq)=kextTable(i0,1:nfreq,imu)*10.**n0w
  salb(1:nfreq)=salbTable(i0,1:nfreq,imu)
  asym(1:nfreq)=asymTable(i0,1:nfreq,imu)
end subroutine getrainp_2


subroutine getsnowp_z(n0w,z13,kext,salb,asym,srate,nfreq)
  use tables2
  real :: n0w, z13
  real, intent(out) :: kext(nfreq),salb(nfreq), asym(nfreq),srate
  integer :: imu
  imu=3
  i0=(z13-10.*n0w-zmin)/dzbin+1
  if(i0<1) i0=1
  if(i0>nbinS2) i0=nbinS2
  kext(1:nfreq)=kextTableS2(i0,1:nfreq,imu)*10.**n0w
  salb(1:nfreq)=salbTableS2(i0,1:nfreq,imu)
  asym(1:nfreq)=asymTableS2(i0,1:nfreq,imu)
  srate=pr13TableS2(i0,imu)*10.**n0w
end subroutine getsnowp_z

subroutine getsnowp_z2(n0w,z13,kext,salb,asym,srate,swc,nfreq)
  use tables2
  real :: n0w, z13
  real, intent(out) :: kext(nfreq),salb(nfreq), asym(nfreq),srate,swc
  integer :: imu
  imu=3
  i0=(z13-10.*n0w-zmin)/dzbin+1
  if(i0<1) i0=1
  if(i0>nbinS2) i0=nbinS2
  kext(1:nfreq)=kextTableS2(i0,1:nfreq,imu)*10.**n0w
  salb(1:nfreq)=salbTableS2(i0,1:nfreq,imu)
  asym(1:nfreq)=asymTableS2(i0,1:nfreq,imu)
  srate=pr13TableS2(i0,imu)*10.**n0w
  swc=10.**pwc13TableS2(i0,imu)*10.**n0w
end subroutine getsnowp_z2

subroutine getrainp_z(n0w,z13,kext,salb,asym,nfreq)
  use tables2
  real :: n0w,srate,z13
  real, intent(out) :: kext(nfreq),salb(nfreq), asym(nfreq)
  integer :: imu
  imu=3
  i0=(z13-10.*n0w-zmin)/dzbin+1
  if(i0<1) i0=1
  if(i0>nbins) i0=nbins
  kext(1:nfreq)=kextTable(i0,1:nfreq,imu)*10.**n0w
  salb(1:nfreq)=salbTable(i0,1:nfreq,imu)
  asym(1:nfreq)=asymTable(i0,1:nfreq,imu)
end subroutine getrainp_z


subroutine getsnow(z13,n0w,swc)
  use tables2
  real :: z13, n0w
  !real, intent(out) :: d0,rrate
  real :: swc
  integer :: imu, i0
  imu=3
  i0=(z13-10.*n0w-zmin)/dzbin+1
  if(i0<1) i0=1
  if(i0>nbinS2) i0=nbinS2
  !rrate=pr13TableS2(i0,imu)*10.**n0w
  swc=10.**(pwc13TableS2(i0,imu)+n0w)
  !print*,z13, n0w, i0, dzbin, swc
  !i0=(z13(i)-10.*n0w(i)-10.*n0bb-zmin)/dzbin+1
  dpia13=att13TableS2(i0,imu)*10.**n0w
  dpia35=att35TableS2(i0,imu)*10.**n0w
end subroutine getsnow !Sept 17, 2015 MG end

subroutine getmsnow(z13,n0w,mswc,dpia13,dpia35)
  use tables2
  real :: z13, n0w
  !real, intent(out) :: d0,rrate
  real :: mswc, dpia13, dpia35
  integer :: imu, i0
  imu=3

  i0=(z13-10.*n0w-zmin)/dzbin+1
  if(i0<1) i0=1
  if(i0>i0max(imu)) i0=i0max(imu)
  !if(i0>nbins) i0=nbins
  ! print*, z13, n0w, i0
  !rrate=pr13TableS2(i0,imu)*10.**n0w
  mswc=10.**(pwc13TableBB(i0,imu)+n0w)
  !print*,z13, n0w, i0, dzbin, swc
  !i0=(z13(i)-10.*n0w(i)-10.*n0bb-zmin)/dzbin+1
  dpia13=att13TableBB(i0,imu)*10.**n0w
  dpia35=att35TableBB(i0,imu)*10.**n0w
end subroutine getmsnow
