!  SFM  08/09/2013  Add use of routines interpol3sfmvc and interpol3sst to
!                    fix problems with default values and to add "best
!                    availale" functionality

subroutine reSampleGMI(gmiData,i1,i2,gmi2Grid,nchunk)
  use f90DataTypes
  implicit none
  type (gmi2GridType) :: gmi2Grid
  type (cgMIDataType) :: gmiData
  real ::  minvalx(5)
  integer :: Ni,Nh1,Nh2,Nh3, i, j, Nw, ii, jj, ic
  real :: x(15), w(600), temp, sumtpw
  integer :: i1, i2
  integer :: nchunk, k, iii
  real    :: xmin, xmax, ymin, ymax, lastGood
  integer :: nx, ny, iside1, iside2, iside
!!print*,'  reSampleGMI',i1,i2,gMIData%nS2f,gMIData%n1b11H
  minvalx=(/80.,170.,140.,200.,230./)
  gMIData%n1b11H=3*(i2-i1+1)-2
!!print*,'  reSampleGMI',i1,i2,gMIData%nS2f,gMIData%n1b11H
  !print*,  gMIData%n1b11H, gMIData%nS1f,gMIData%nrays
  !stop
  allocate(gMIData%gmiS13(gMIData%nS1f,gMIData%nrays,gMIData%n1b11H)) 
  allocate(gMIData%emissS13(gMIData%nS1f,gMIData%nrays,gMIData%n1b11H)) 
  allocate(gMIData%gmiS13H(gMIData%nS1f,gMIData%nrays,gMIData%n1b11H))
  allocate(gMIData%gmiS23(gMIData%nS2f,gMIData%nrays,gMIData%n1b11H))
  allocate(gMIData%S1lon3(gMIData%nrays,gMIData%n1b11H)) 
  allocate(gMIData%SClat3(gMIData%n1b11H)) 
  allocate(gMIData%SClon3(gMIData%n1b11H)) 
  allocate(gMIData%S1lat3(gMIData%nrays,gMIData%n1b11H)) 
  allocate(gMIData%S1eia3(gMIData%nrays,gMIData%n1b11H)) 
  allocate(gMIData%S2eia3(gMIData%nrays,gMIData%n1b11H)) 

  allocate(gMIData%sfc_wind3(gMIData%nrays,gMIData%n1b11H)) 
  allocate(gMIData%sst3(gMIData%nrays,gMIData%n1b11H)) 
  allocate(gMIData%tpw3(gMIData%nrays,gMIData%n1b11H)) 
  !print*, 'reSampleGMI',gMIData%nrays,gMIData%n1b11H, i1, i2
  allocate(gMIData%tb10hC3(gMIData%nrays,gMIData%n1b11H)) 
  allocate(gMIData%tb10vC3(gMIData%nrays,gMIData%n1b11H)) 
  allocate(gMIData%tb19hC3(gMIData%nrays,gMIData%n1b11H)) 
  allocate(gMIData%tb19vC3(gMIData%nrays,gMIData%n1b11H)) 
  allocate(gMIData%tb21vC3(gMIData%nrays,gMIData%n1b11H)) 
  allocate(gMIData%tb37hC3(gMIData%nrays,gMIData%n1b11H)) 
  allocate(gMIData%tb37vC3(gMIData%nrays,gMIData%n1b11H)) 

!begin  MG 4/25/14 add code to avoid error in empty GMI file
if(i2<i1) then
    gmi2Grid%ny=1
    gmi2Grid%ymin=-999.
    gmi2Grid%nx=1
    gmi2Grid%xmin=-999.
    allocate(gmi2Grid%ig(gmi2Grid%nx,gmi2Grid%ny))
    allocate(gmi2Grid%jg(gmi2Grid%nx,gmi2Grid%ny))
    allocate(gmi2Grid%actOb(gmi2Grid%nx,gmi2Grid%ny))
    gmi2Grid%ig=-99
    gmi2Grid%jg=-99
    return 
 endif
!end    MG 4/25/14

  CALL interpol3sfmvc(gMIData%gmiS1(:,:,i1:i2),gMIData%gmiS13,                    &
                   gMIData%nS1f,i2-i1+1)

  gMIData%sfc_wind3=-99
  
  CALL interpol3sc(gMIData%sfc_wind(:,i1:i2), gMIData%sfc_wind3,               &
                   gMIData%landSea(:,i1:i2),i2-i1+1) 

  CALL interpol3sst(gMIData%sst(:,i1:i2), gMIData%sst3,                        &
                   gMIData%landSea(:,i1:i2),i2-i1+1) 

  CALL interpol3sc(gMIData%tpw(:,i1:i2), gMIData%tpw3,                         &
                   gMIData%landSea(:,i1:i2),i2-i1+1) 

  CALL interpol3sc(gMIData%tb10hC(:,i1:i2), gMIData%tb10hC3,                   &
                   gMIData%landSea(:,i1:i2),i2-i1+1) 

  CALL interpol3sc(gMIData%tb10vC(:,i1:i2),gMIData%tb10vC3,                    &
                   gMIData%landSea(:,i1:i2),i2-i1+1) 

  CALL interpol3sc(gMIData%tb19hC(:,i1:i2),gMIData%tb19hC3,                    &
                   gMIData%landSea(:,i1:i2),i2-i1+1) 

  CALL interpol3sc(gMIData%tb19vC(:,i1:i2),gMIData%tb19vC3,                    &
                   gMIData%landSea(:,i1:i2),i2-i1+1) 

  CALL interpol3sc(gMIData%tb21vC(:,i1:i2),gMIData%tb21vC3,                    &
                   gMIData%landSea(:,i1:i2),i2-i1+1) 

!  print*, "tb21vC3", maxval(gMIData%tb21vC3)
!  print*, "tb21vC", maxval(gMIData%tb21vC(:,i1:i2))

  call interpol3sc(gMIData%tb37hC(:,i1:i2),gMIData%tb37hC3,&
       gMIData%landSea(:,i1:i2),i2-i1+1) 

!  print*, "tb37hC3", maxval(gMIData%tb37hC3)
!  print*, "tb37hC", maxval(gMIData%tb37hC(:,i1:i2))

  call interpol3sc(gMIData%tb37vC(:,i1:i2),gMIData%tb37vC3,&
       gMIData%landSea(:,i1:i2),i2-i1+1) 

!  print*, "tb37vC3", maxval(gMIData%tb37vC3)
!  print*, "tb37vC", maxval(gMIData%tb37vC(:,i1:i2))


  call interpol3xc(gMIData%S1lon(:,i1:i2),gMIData%S1lon3,i2-i1+1) 
  call interpol3xc(gMIData%S1lat(:,i1:i2),gMIData%S1lat3,i2-i1+1) 
  
  call interpol3xc(gMIData%S1eia(:,i1:i2),gMIData%S1eia3,i2-i1+1) 
  call interpol3xc(gMIData%S2eia(:,i1:i2),gMIData%S2eia3,i2-i1+1) 
  
  call interpol3sc2(gMIData%SClon(i1:i2),gMIData%SClon3,i2-i1+1) 
  call interpol3sc2(gMIData%SClat(i1:i2),gMIData%SClat3,i2-i1+1) 

!  print*, gMIData%SClon3
!  print*, gMIData%SClat3
!  stop
!  print*, gmiData%S1lon3(70,:)
  call interpol3x(gMIData%gmiS2(1,:,i1:i2),gMIData%gmiS23(1,:,:),i2-i1+1)  
  call interpol3x(gMIData%gmiS2(2,:,i1:i2),gMIData%gmiS23(2,:,:),i2-i1+1)  

 !4/15/14 MG Begin
  gmi2Grid%dx=0.05
  gmi2Grid%xmin=+999.9
!  print*, 3*(i2-i1+1)-2
!  nx= 1.1*(3*(i2-i1+1)-2)
  ymax=-999.9
  ymin=999.9
  xmax=-999.9
  xmin=999.9

  do i=1,3*(i2-i1+1)-2
     do j=70,150
        if(gMIData%S1lat3(j,i)>-998) then
           if(gMIData%S1lat3(j,i)>ymax) then
              ymax=gMIData%S1lat3(j,i)
           endif
           if(gMIData%S1lat3(j,i)<ymin) then
              ymin=gMIData%S1lat3(j,i)
           endif
        endif
        if(gMIData%S1lon3(j,i)>-998) then
           if(gMIData%S1lon3(j,i)>xmax) then
              xmax=gMIData%S1lon3(j,i)
              !print*, i, j, gMIData%S1lon3(j,i)
           endif
           if(gMIData%S1lon3(j,i)<xmin) then
              xmin=gMIData%S1lon3(j,i)
           endif
        endif
     enddo
  enddo
 
        
  !print*, ymin, ymax, xmin, xmax
  if(ymax.gt.-998) then
     gmi2Grid%ny=int((ymax-ymin)/gmi2Grid%dx)+2
     gmi2Grid%ymin=ymin
     !print*, ymin,ymax,gmi2Grid%ny
  else
     gmi2Grid%ny=1
     gmi2Grid%ymin=-999.
  endif
 if(xmax.gt.-998) then
     gmi2Grid%nx=int((xmax-xmin)/gmi2Grid%dx)+2
     gmi2Grid%xmin=xmin
     !print*, xmin, xmax,gmi2Grid%nx
  else
     gmi2Grid%nx=1
     gmi2Grid%xmin=-999.
  endif
!4/15/14 MG End
  allocate(gmi2Grid%ig(gmi2Grid%nx,gmi2Grid%ny))
  allocate(gmi2Grid%jg(gmi2Grid%nx,gmi2Grid%ny))
  allocate(gmi2Grid%actOb(gmi2Grid%nx,gmi2Grid%ny))
!  return
  gmi2Grid%ig=-99
  gmi2Grid%jg=-99
  gmi2Grid%actOb=0
!  SFM  begin  12/06/2013; accomodate conditions of all-empty 1CGMI files
  IF (i1 .gt. 1 .and. i2 .gt. 1) THEN
     do i=70,150
        do j=1,gMIData%n1b11H
           if(gMIData%S1lon3(i,j)>-999 .and.  gMIData%S1lat3(i,j)>-999) then
              ii=(gMIData%S1lon3(i,j)-gmi2Grid%xmin)/gmi2Grid%dx+1
              jj=(gMIData%S1lat3(i,j)-gmi2Grid%ymin)/gmi2Grid%dx+1
              gmi2Grid%ig(ii,jj)=i
              gmi2Grid%jg(ii,jj)=j
           endif
        enddo
     enddo
     do i=70,150
        do j=i1,i2
           if(gMIData%S1lon(i,j)>-999 .and.  gMIData%S1lat(i,j)>-999) then
              ii=(gMIData%S1lon(i,j)-gmi2Grid%xmin)/gmi2Grid%dx+1
              jj=(gMIData%S1lat(i,j)-gmi2Grid%ymin)/gmi2Grid%dx+1
              if(ii.gt.0.and.jj.gt.0.and.ii.lt.gmi2Grid%nx.and.&
                   jj.lt.gmi2Grid%ny)then
                 gmi2grid%actOb(ii,jj)=1
              endif
           endif
        enddo
     enddo
  ENDIF
!  SFM  end    12/06/2013

end subroutine reSampleGMI


subroutine deallocateHResGMI(gmiData,gmi2Grid)
  use f90DataTypes
  implicit none
  type (gmi2GridType) :: gmi2Grid
  type (gMIDataType) :: gmiData
  real ::  minvalx(5)
  integer :: Ni,Nh1,Nh2,Nh3, i, j, Nw, ii, jj, ic
  real :: x(15), w(600), temp, sumtpw
  real :: i1, i2
 
 
  deallocate(gMIData%gmilow3) 
  deallocate(gMIData%gmilow3H)
  deallocate(gMIData%txlon3) 
  deallocate(gMIData%txlat3) 
  deallocate(gmi2Grid%ig)
  deallocate(gmi2Grid%jg)
  deallocate(gmi2Grid%actOb)
  deallocate(gMIData%sfc_wind3) 
  deallocate(gMIData%sst3) 
  deallocate(gMIData%tpw3) 
  deallocate(gMIData%tb10hC3) 
  deallocate(gMIData%tb10vC3) 
  deallocate(gMIData%tb19hC3) 
  deallocate(gMIData%tb19vC3) 
  deallocate(gMIData%tb21vC3) 
  deallocate(gMIData%tb37hC3) 
  deallocate(gMIData%tb37vC3) 
  deallocate(gMIData%gmihigh3)
  return  
 
end subroutine deallocateHResGMI


subroutine deallocateHRescGMI(gmiData,gmi2Grid)
  use f90DataTypes
  implicit none
  type (gmi2GridType) :: gmi2Grid
  type (cgMIDataType) :: gmiData
  real ::  minvalx(5)
  integer :: Ni,Nh1,Nh2,Nh3, i, j, Nw, ii, jj, ic
  real :: x(15), w(600), temp, sumtpw
  real :: i1, i2
 
 
  deallocate(gMIData%gmiS13) 
  deallocate(gMIData%emissS13)
  deallocate(gMIData%gmiS13H)
  deallocate(gMIData%S1lon3) 
  deallocate(gMIData%S1lat3) 
  deallocate(gMIData%S1eia3)
  deallocate(gMIData%S2eia3)
  deallocate(gMIData%SClon3) 
  deallocate(gMIData%SClat3)
  if(allocated(gmi2grid%actOb)) deallocate(gmi2grid%actOb) 
  if(allocated(gmi2grid%ig)) deallocate(gmi2Grid%ig)
  if(allocated(gmi2grid%jg)) deallocate(gmi2Grid%jg)
  deallocate(gMIData%sfc_wind3) 
  deallocate(gMIData%sst3) 
  deallocate(gMIData%tpw3) 
  deallocate(gMIData%tb10hC3) 
  deallocate(gMIData%tb10vC3) 
  deallocate(gMIData%tb19hC3) 
  deallocate(gMIData%tb19vC3) 
  deallocate(gMIData%tb21vC3) 
  deallocate(gMIData%tb37hC3) 
  deallocate(gMIData%tb37vC3) 
  deallocate(gMIData%gmiS23)
  return  

  
 
end subroutine deallocateHRescGMI
