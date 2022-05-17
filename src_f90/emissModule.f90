module emissMod
 integer *4                n_lat
 parameter(n_lat = 720)
 integer *4                n_lon
 parameter(n_lon = 1440)
 integer *4                n_chan
 parameter(n_chan = 6) !SJM 9/9/2015
 real  :: gmiobslow(7)
 real  :: gmiobshigh(2)
 real *4                   emissv_map(n_lon, n_lat, n_chan)
 real *4                   emissh_map(n_lon, n_lat, n_chan)
 real *4                   emissv_std_map(n_lon, n_lat, n_chan)
 real *4                   emissh_std_map(n_lon, n_lat, n_chan)
 real eLat, eLon
 real ::  kexttotEns(7,40,100), salbtotEns(7,40,100), asymtotEns(7,40,100)
end module emissMod

subroutine reademiss(imonth)
  use emissMod
  implicit none

  integer imonth
  select case(imonth)
  case(1)
     open(20,file='Emiss/emiss1.bin',form='unformatted')
  case(2)
     open(20,file='Emiss/emiss2.bin',form='unformatted')
  case(3)
     open(20,file='Emiss/emiss3.bin',form='unformatted')
  case(4)
     open(20,file='Emiss/emiss4.bin',form='unformatted')
  case(5)
     open(20,file='Emiss/emiss5.bin',form='unformatted')
  case(6)
     open(20,file='Emiss/emiss6.bin',form='unformatted')
  case(7)
     open(20,file='Emiss/emiss7.bin',form='unformatted')
  case(8)
     open(20,file='Emiss/emiss8.bin',form='unformatted')
  case(9)
     open(20,file='Emiss/emiss9.bin',form='unformatted')
  case(10)
     open(20,file='Emiss/emiss10.bin',form='unformatted')
  case(11)
     open(20,file='Emiss/emiss11.bin',form='unformatted')
  case(12)
     open(20,file='Emiss/emiss12.bin',form='unformatted')
  end select
  
  read(20) emissv_map(:,:,1:5)
  read(20) emissh_map(:,:,1:5)
  read(20) emissv_std_map(:,:,1:5)
  read(20) emissh_std_map(:,:,1:5)
  close(20)
  emissv_std_map(:,:,6) = emissv_std_map(:,:,5)
  emissh_std_map(:,:,6) = emissh_std_map(:,:,5)
end subroutine reademiss

subroutine getemiss(lat,lon,snowicecover,emissv,emissh,emissv_std,emissh_std)
  use emissMod
  use LUT_def
  implicit none
  real                    lat, lon
  integer                 snowicecover
  real                    emissv(n_chan)
  real                    emissh(n_chan)
  real                    emissv_std(n_chan)
  real                    emissh_std(n_chan)
  integer :: ii, jj, stype, i

!get surface type from map
  !print*, lat, lon
  ii=2880-floor((lat+90.)/180.*2880.)
  if(ii .lt. 1) ii=1
  if(ii .gt. 2880) ii = 2880
  jj=floor((lon+180.)/360.*5760.)+1
  if(jj .lt. 1) jj = 1
  if(jj .gt. 5760) jj = 5760
  !print*, ii,jj
  !stype = LUT%land_class_map_bare(jj,ii)
  if(snowicecover .eq. 2 .or. snowicecover .eq. 3) then
    stype = LUT%land_class_map_snow(jj,ii) !SJM 9/9/15
  else 
    stype = LUT%land_class_map_snow(jj,ii) !SJM 9/9/15
  endif
  !print*, lat, lon, stype
!begin  WSO 11/7/13 added tests for indices out of bounds
!modified by SJM for use with GMI-derived maps
  ii=min(int((lat+70.)/0.25) + 1, 560)
  jj=min(int((lon+180.)/.25)+1, 1440)
  
!end  WSO  11/7/13
   if(stype .eq. 1) then
     emissv=-99.
     emissv_std = 0.
     emissh = -99.
     emissh_std=0.
     return
   endif
   !if(stype .eq. 2 .or. (stype .ge. 8 .and. stype .le. 11) .or. stype .eq. 14) then !sea ice or snow cover
   if(snowicecover .eq. 2 .or. snowicecover .eq. 3) then
    emissv(1) = LUT%emis_map_snow(jj,ii,1)
    emissv(2) = LUT%emis_map_snow(jj,ii,3)
    emissv(4) = LUT%emis_map_snow(jj,ii,5)
    emissv(5) = LUT%emis_map_snow(jj,ii,7)
    emissv(6) = LUT%emis_map_snow(jj,ii,9)
    emissh(1) = LUT%emis_map_snow(jj,ii,2)
    emissh(2) = LUT%emis_map_snow(jj,ii,4)
    emissh(4) = LUT%emis_map_snow(jj,ii,6)
    emissh(5) = LUT%emis_map_snow(jj,ii,8)
    emissh(6) = LUT%emis_map_snow(jj,ii,10)
!begin  SJM 2/14/17 read emissivity std deviations
    emissv_std(1) = LUT%estd_map_snow(jj,ii,1)
    emissv_std(2) = LUT%estd_map_snow(jj,ii,3)
    emissv_std(4) = LUT%estd_map_snow(jj,ii,5)
    emissv_std(5) = LUT%estd_map_snow(jj,ii,7)
    emissv_std(6) = LUT%estd_map_snow(jj,ii,9)
    emissh_std(1) = LUT%estd_map_snow(jj,ii,2)
    emissh_std(2) = LUT%estd_map_snow(jj,ii,4)
    emissh_std(4) = LUT%estd_map_snow(jj,ii,6)
    emissh_std(5) = LUT%estd_map_snow(jj,ii,8)
    emissh_std(6) = LUT%estd_map_snow(jj,ii,10)
!end  SJM 2/14/17
  else if(snowicecover .eq. 0 .or. snowicecover .eq. 1) then
    emissv(1) = LUT%emis_map_bare(jj,ii,1)
    emissv(2) = LUT%emis_map_bare(jj,ii,3)
    emissv(4) = LUT%emis_map_bare(jj,ii,5)
    emissv(5) = LUT%emis_map_bare(jj,ii,7)
    emissv(6) = LUT%emis_map_bare(jj,ii,9)
    emissh(1) = LUT%emis_map_bare(jj,ii,2)
    emissh(2) = LUT%emis_map_bare(jj,ii,4)
    emissh(4) = LUT%emis_map_bare(jj,ii,6)
    emissh(5) = LUT%emis_map_bare(jj,ii,8)
    emissh(6) = LUT%emis_map_bare(jj,ii,10)
!begin  SJM 2/14/17 read emissivity std deviations
    emissv_std(1) = LUT%estd_map_bare(jj,ii,1)
    emissv_std(2) = LUT%estd_map_bare(jj,ii,3)
    emissv_std(4) = LUT%estd_map_bare(jj,ii,5)
    emissv_std(5) = LUT%estd_map_bare(jj,ii,7)
    emissv_std(6) = LUT%estd_map_bare(jj,ii,9)
    emissh_std(1) = LUT%estd_map_bare(jj,ii,2)
    emissh_std(2) = LUT%estd_map_bare(jj,ii,4)
    emissh_std(4) = LUT%estd_map_bare(jj,ii,6)
    emissh_std(5) = LUT%estd_map_bare(jj,ii,8)
    emissh_std(6) = LUT%estd_map_bare(jj,ii,10)
!end  SJM 2/14/17
  else
    emissv(1) = LUT%emis_map_all(jj,ii,1)
    emissv(2) = LUT%emis_map_all(jj,ii,3)
    emissv(4) = LUT%emis_map_all(jj,ii,5)
    emissv(5) = LUT%emis_map_all(jj,ii,7)
    emissv(6) = LUT%emis_map_all(jj,ii,9)
    emissh(1) = LUT%emis_map_all(jj,ii,2)
    emissh(2) = LUT%emis_map_all(jj,ii,4)
    emissh(4) = LUT%emis_map_all(jj,ii,6)
    emissh(5) = LUT%emis_map_all(jj,ii,8)
    emissh(6) = LUT%emis_map_all(jj,ii,10)
!begin  SJM 2/14/17 read emissivity std deviations
    emissv_std(1) = LUT%estd_map_all(jj,ii,1)
    emissv_std(2) = LUT%estd_map_all(jj,ii,3)
    emissv_std(4) = LUT%estd_map_all(jj,ii,5)
    emissv_std(5) = LUT%estd_map_all(jj,ii,7)
    emissv_std(6) = LUT%estd_map_all(jj,ii,9)
    emissh_std(1) = LUT%estd_map_all(jj,ii,2)
    emissh_std(2) = LUT%estd_map_all(jj,ii,4)
    emissh_std(4) = LUT%estd_map_all(jj,ii,6)
    emissh_std(5) = LUT%estd_map_all(jj,ii,8)
    emissh_std(6) = LUT%estd_map_all(jj,ii,10)
!end  SJM 2/14/17
  endif

  if(emissv(1) .le. 0.) emissv(1) = LUT%land_class_emis_mean(stype-1,1)
  if(emissv(2) .le. 0.) emissv(2) = LUT%land_class_emis_mean(stype-1,3)
  if(emissv(4) .le. 0.) emissv(4) = LUT%land_class_emis_mean(stype-1,5)
  if(emissv(5) .le. 0.) emissv(5) = LUT%land_class_emis_mean(stype-1,7)
  if(emissv(6) .le. 0.) emissv(6) = LUT%land_class_emis_mean(stype-1,9)
  emissv(3) = 0.64*emissv(2)+0.36*emissv(4)
  if(emissh(1) .le. 0.) emissh(1) = LUT%land_class_emis_mean(stype-1,2)
  if(emissh(2) .le. 0.) emissh(2) = LUT%land_class_emis_mean(stype-1,4)
  if(emissh(4) .le. 0.) emissh(4) = LUT%land_class_emis_mean(stype-1,6)
  if(emissh(5) .le. 0.) emissh(5) = LUT%land_class_emis_mean(stype-1,8)
  if(emissh(6) .le. 0.) emissh(6) = LUT%land_class_emis_mean(stype-1,10)
  emissh(3) = 0.64*emissh(2)+0.36*emissh(4)
  
  
  !emissv_std=emissv_std_map(jj,ii,:)
  !emissh_std=emissh_std_map(jj,ii,:)
  !if(minval(emissv)<0) then
   ! print*, lat, lon
   ! print*, ii, jj
   ! print*, emissv
   ! stop
  !endif
end subroutine getemiss
! subroutine getemiss(lat,lon,emissv,emissh,emissv_std,emissh_std)
!   use emissMod
!   implicit none
!   real                    lat, lon
!   real                    emissv(n_chan)
!   real                    emissh(n_chan)
!   real                    emissv_std(n_chan)
!   real                    emissh_std(n_chan)
!   integer :: ii, jj
! 
! !begin  WSO 11/7/13 added tests for indices out of bounds
!   ii=min(int((lat+90.)/0.25) + 1, 720)
!   if(lon<0) then
!      jj=min(int((lon+360)/.25)+1, 1440)
!   else
!      if(lon>=360) then
!         jj=min(int((lon-360)/.25) + 1, 1440)
!      else
!         jj=min(int(lon/.25) + 1, 1440)
!      endif
!   endif
! !end  WSO  11/7/13
! 
!   emissv=emissv_map(jj,ii,:)
!   emissh=emissh_map(jj,ii,:)
!   emissv_std=emissv_std_map(jj,ii,:)
!   emissh_std=emissh_std_map(jj,ii,:)
!   if(minval(emissv)<0) then
!    ! print*, lat, lon
!    ! print*, ii, jj
!    ! print*, emissv
!    ! stop
!   endif
! end subroutine getemiss
