module cldclass
  integer:: nlayer, nc
  integer:: nRhEofs, nCldwEofs
  real,allocatable :: temp(:), & ! temperature profile used in the radiative transfer module
       rlhm(:), &     ! mean relative humidity profile, currently, not used
       height(:), &   ! height
       press(:), &    ! pressure
       rlhmrc(:,:), & ! relative humidity profiles in precipitation regions, 
       rlhmcs(:,:), & ! relative humidity profiles in clear sky regions
       cc(:,:), &      ! cloud profiles
       rheofs(:,:), cldweofs(:,:), cldwm(:), stdrhPC(:), stdcldwPC(:)
  
 
  real, allocatable :: tpw(:), tpwR(:)
  real, allocatable :: atm_extKa(:,:), atm_extKag(:,:), cld_extKa(:,:), cld_extKag(:,:),   atm_extm(:), &
       cld_extm(:), atm_exts(:), cld_exts(:) 
                      ! atm_extKa(:,:), cld_extKa(:,:) water vapor and cloud extinction at Ka-band
                      ! atm_extm(:), cld_extm(:), atm_exts(:), cld_exts(:) 
                      ! were intended to hold statistics, i.e. mean, 
                      ! standard deviation, but ended up not used
  real,allocatable :: atm_extMw(:,:,:),  atm_extMwRg(:,:,:) 
! precomputed extinction at various heights and frequencies
  real,allocatable :: atm_extMwr(:,:,:) 
! precomputed extinction at various heights and frequencies
  real :: drrte       ! this is the layer thickness in the RT calculations, technically it is not needed
                      ! but used to expedite the calculation of the inclusion of radar derived 
                      ! electromagnetic properties into calculations

end module cldclass

subroutine cloud_init(nmfreq)
! this subroutines sets the height, temperature, and nlayer variables
! for the radiative transfer subroutine
! it also reads 50 Relative Humidity  profiles and 50 cloud profiles
! the extinction associated with these profiles is precomputed in used in radiative calculations
! the attenuation at Ka-band is also pre-computed 
use cldclass
use microwFreq
use ran_mod
implicit none
integer :: n, i, j, ifreq, nmfreq
real :: tsfc,dr,tavg
real :: psfc, freq35
real :: junk(40)
real :: pav, qvsat, es,  pwvsat, H(100), qv0(100), qv0r(100)
freq35=35.
n=40
drrte=0.5
psfc=1000.
nlayer=n
nc=50
nRhEofs=20
nCldwEofs=20

if(allocated(temp)) return
allocate(temp(0:n),height(0:n),press(0:n),rlhm(n))
allocate(rlhmrc(n,nc), rlhmcs(n,nc), cc(n,nc))
allocate(atm_extKa(n,nc), cld_extKa(n,nc), atm_extKag(n,50), cld_extKag(n,50))
allocate(atm_extm(n), cld_extm(n))
allocate(atm_exts(n), cld_exts(n))
allocate(atm_extMw(nmfreq,n,nc))
allocate(atm_extMwR(nmfreq,n,nc),atm_extMwRg(nmfreq,nc,50))
allocate(tpw(nc))
allocate(tpwR(nc))
allocate(rheofs(nRhEofs,n), cldweofs(nCldwEofs,n), cldwm(n), stdrhPC(nRhEofs), stdcldwPC(nCldwEofs))

open(10,file='AncData/rhCloud.dat')
read(10,*) (rlhm(i),i=1,n)

do i=1,nRhEofs
   read(10,*) (rheofs(i,j),j=1,n)
enddo
read(10,*) (stdrhPC(j),j=1,nRhEofs)

read(10,*) (cldwm(i),i=1,n)
do i=1,ncldwEofs
   read(10,*) (cldweofs(i,j),j=1,n)
enddo
read(10,*) (stdcldwPC(j),j=1,ncldwEofs)
!print*, rlhm
!stop
!print*, stdcldwPC(1:10)
close(10)
mfreq(3)=21.3
open(10,file='AncData/clusters80.dat')
do i=1,nc
   read(10,*) (rlhmrc(j,i),j=1,nlayer), (junk(j),j=n+1,40)
   read(10,*) (cc(j,i),j=1,nlayer), (junk(j),j=n+1,40)
   read(10,*) (rlhmcs(j,i),j=1,nlayer), (junk(j),j=n+1,40)

enddo
cc=cc*0.5
101 format(50(F7.3,1x))


tpw=0
read(10,*) (temp(j),j=0,nlayer)
do i=0,n
   height(i)=i*drrte
enddo

rlhmrc=rlhmrc*0.975
press(0)=psfc*( 1. - height(0)*6.5/288. )**5.253

do i=1,n
   Tavg= (temp(i)+temp(i-1))*.5
   press(i)  = psfc*( 1. - height(i)*6.5/288. )**5.253
   pav=0.5*(press(i)+press(i-1))
   do j=1,nc
      call absorb(freq35, tavg, press(i), &
           rlhmrc(i,j), atm_extKa(i,j))
      call gcloud(freq35,tavg,cc(i,j),cld_extKa(i,j))
   enddo
   do j=1,50
      call absorb(freq35, tavg, press(i), &
           j*2., atm_extKag(i,j))
   enddo
   do ifreq=1,nmfreq
      do j=1,nc
         call absorb(mfreq(ifreq),tavg,press(i),rlhmcs(i,j),&
              atm_extMw(ifreq,i,j))
         call absorb(mfreq(ifreq),tavg,press(i),rlhmrc(i,j),&
              atm_extMwR(ifreq,i,j))
      enddo
      do j=1,50
         call absorb(mfreq(ifreq),tavg,press(i),j*2.,&
              atm_extMwRg(ifreq,i,j))
      enddo
   enddo
   atm_extm(i)=sum(atm_extKa(i,:))/nc
   cld_extm(i)=sum(cld_extKa(i,:))/nc
   atm_exts(i)=sqrt(sum((atm_extKa(i,:)-atm_extm(i))**2)/nc)
   cld_exts(i)=sqrt(sum((cld_extKa(i,:)-cld_extm(i))**2)/nc)
enddo

atm_extKa=atm_extKa*4.43
atm_extKag=atm_extKag*4.43
cld_extKa=cld_extKa*4.43

end subroutine cloud_init




 Real Function pwvsat ( T )
   Implicit NONE
!
   Real :: x,T
!
! Goff Gratch (1946) equation
!
   x = 373.16/T
   pwvsat = -7.90298 * (x - 1.) + 5.02808 * log10(x) -     &
             1.3816e-7 * (10.**(11.344*(1.-1./x)) - 1.) +  &
             8.1328e-3 * (10.**(-3.49149*(x-1.)) - 1.) + log10(1013.246)
   pwvsat = 10.**pwvsat
!
 End Function pwvsat
