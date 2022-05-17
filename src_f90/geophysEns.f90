module geophysEns
  real, allocatable :: nwCoarse(:,:,:,:)
  real, allocatable :: emissCoarse(:,:,:)
  real, allocatable :: rhPCCoarse(:,:,:,:)
  real, allocatable :: cldwPCCoarse(:,:,:,:)
  real :: cldwcoeff(20,100)
  integer, allocatable :: icCoarse(:,:,:)
  integer, allocatable :: jcCoarse(:,:,:)
  integer, allocatable :: iindex(:)
  integer, allocatable :: jindex(:)
  integer  ::  nc, nt, nh, nm, nemiss
  integer  :: iRad, jRad
  contains 
    subroutine allocGeophys(nc1,nt1,nh1,nm1,nemiss1)
      integer :: nc1, nt1, nh1,nm1
      nc=nc1
      nt=nt1
      nh=nh1 
      nm=nm1
      nemiss=nemiss1
      allocate( nwCoarse(nc,nt,nh,nm))
      allocate( emissCoarse(nc,nt,nemiss))
      allocate( rhPCCoarse(nc,nt,nm,20))
      allocate( cldwPCCoarse(nc,nt,nm,20))
      allocate( iindex(nc))
      allocate( jindex(nt))
      allocate( icCoarse(nc,nt,nm))
      allocate( jcCoarse(nc,nt,nm))
    end subroutine allocGeophys
    subroutine deallocGeophys()
      integer :: ny, nx, nh
! SFM  begin  03/27/2014; execution protections
      IF (ALLOCATED(nwCoarse))     deallocate(nwCoarse)
      IF (ALLOCATED(emissCoarse))  deallocate(emissCoarse)
      IF (ALLOCATED(icCoarse))     deallocate(icCoarse)
      IF (ALLOCATED(jcCoarse))     deallocate(jcCoarse)
      IF (ALLOCATED(iindex))       deallocate(iindex)
      IF (ALLOCATED(jindex))       deallocate(jindex)
      IF (ALLOCATED(rhPCCoarse))   deallocate(rhPCCoarse)
      IF (ALLOCATED(cldwPCCoarse)) deallocate(cldwPCCoarse)
! SFM  end    03/27/2014
    end subroutine deallocGeophys
    subroutine setdNwIcJc(dn,nmemb)
      use ran_mod
      implicit none
      integer :: i, j, k, l, nmemb
      real :: dn
      iindex=(/1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 49/)
      do i=1,60
         jindex(i+1)=(i)*5
      enddo
      jindex(1)=1
      do i=1,11
         do j=1,61
            do l=1,nmemb
               icCoarse(i,j,l)= ran1()*(49-0.1)+1
               if(icCoarse(i,j,l)>49) icCoarse(i,j,l)=49
               jcCoarse(i,j,l)= ran1()*(49-0.1)+1
               if(jcCoarse(i,j,l)>49) jcCoarse(i,j,l)=49
            enddo
            do l=1,nmemb
               do k=1,9
                  if(k==1) then
                     nwCoarse(i,j,k,l)=(ran1()-0.5)*0.25+dn
                  else
                     nwCoarse(i,j,k,l)=((ran1()-0.5)*0.25+dn)*.25+0.75*nwCoarse(i,j,k-1,l)
                  endif
                  if(nwCoarse(i,j,k,l)>1.5) nwCoarse(i,j,k,l)=1.5
                  if(nwCoarse(i,j,k,l)<-1.5) nwCoarse(i,j,k,l)=-1.5
               enddo
            enddo
         enddo
      enddo
    end subroutine setdNwIcJc
    subroutine setdNwIcJcL(dn,nmemb)
      use ran_mod
      implicit none
      integer :: i, j, k, l, nmemb
      integer :: i0, j0, i1, j1
      real :: dn
      iindex=(/1,  10,  20,  30,  40,  49/)
      do i=1,50
         jindex(i+1)=(i)*10
      enddo
      jindex(1)=1
      do i=1,6
         do j=1,61     !MG 9/25/13
            do l=1,20
               do k=1,nmemb
                  rhPCCoarse(i,j,k,l)=normal2(0.,1.0)
                  cldwPCCoarse(i,j,k,l)=normal2(0.,1.0)
               enddo
            enddo
            do l=1,nemiss
               emissCoarse(i,j,l)=normal2(0.,1.)
            enddo
            do l=1,nmemb
               icCoarse(i,j,l)= ran1()*(49-0.1)+1
               if(icCoarse(i,j,l)>49) icCoarse(i,j,l)=49
               jcCoarse(i,j,l)= ran1()*(49-0.1)+1
               if(jcCoarse(i,j,l)>49) jcCoarse(i,j,l)=49
            enddo
            do l=1,nmemb
               do k=1,9
                  if(k==1) then
                     nwCoarse(i,j,k,l)=normal2(0.,1.)*0.5
                  else
!  SFM  begin  06222014; for M.Grecu (unknown justification)
                     nwCoarse(i,j,k,l)=normal2(0.,1.)*0.5*0.5+0.5*nwCoarse(i,j,k-1,l)
!  SFM  end    06222014
                  endif
                  if(nwCoarse(i,j,k,l)>2.5) nwCoarse(i,j,k,l)=2.5
                  if(nwCoarse(i,j,k,l)<-2.5) nwCoarse(i,j,k,l)=-2.5
               enddo
!  SFM  begin  06222014; for M.Grecu (unknown justification)
!               nwCoarse(i,j,1:5,l)=0.5*normal2(0.,1.)+0.5
!               nwCoarse(i,j,6:9,l)=0.5*normal2(0.,1.)
!  SFM  end    06222014
            enddo
         enddo
      enddo
    end subroutine setdNwIcJcL

    subroutine interpoldNw(i0,j0, logNw)
      use ran_mod
      implicit none
      integer :: i, j, k
      integer :: i1,j1,i0,j0
      real :: dn
      real :: dx, dy, f1, f2
      real, dimension(:) :: logNw
      i1=1
      j1=1

      do while (i0>iindex(i1))
         i1=i1+1
      end do
      do while (j0>jindex(j1))
         j1=j1+1
      end do
      if (i1 == 1) i1=2
      if (j1 == 1) j1=2
      i1=i1-1
      j1=j1-1
      
      dx=iindex(i1+1)-iindex(i1)
      dy=jindex(j1+1)-jindex(j1)
      f1=(i0-iindex(i1))/dx
      f2=(j0-jindex(j1))/dy
      if(f1<0 .or. f1>1) then
         write(*,*) iindex(i1),i0, iindex(i1+1), f1
         stop
      endif
      if(f2<0 .or. f2>1) then
         write(*,*) jindex(j1),j0, jindex(j1+1), f2
         stop
      endif
      do i=1,nm
         logNw((i-1)*nh+1:i*nh)=(1-f1)*(1-f2)*nwCoarse(i1,j1,:,i)+ &
              (1-f1)*f2*nwCoarse(i1,j1+1,:,i)+ &
              f1*(1-f2)*nwCoarse(i1+1,j1,:,i)+ &
              f1*f2*nwCoarse(i1+1,j1+1,:,i)
      enddo
      !write(*,*) logNw, f1, f2, i1, i2, j1, j2
    end subroutine interpoldNw
    
    subroutine interpolEmiss(i0,j0, emiss)
      use ran_mod
      implicit none
      integer :: i, j, k
      real :: dn
      real :: dx, dy, f1, f2
      integer :: i0,j0, i1, j1
      real, dimension(:) :: emiss
      i1=1
      j1=1

      do while (i0>iindex(i1))
         i1=i1+1
      end do
      do while (j0>jindex(j1))
         j1=j1+1
      end do
      if (i1 == 1) i1=2
      if (j1 == 1) j1=2
      i1=i1-1
      j1=j1-1
      
      dx=iindex(i1+1)-iindex(i1)
      dy=jindex(j1+1)-jindex(j1)
      f1=(i0-iindex(i1))/dx
      f2=(j0-jindex(j1))/dy
     
      if(f1<0 .or. f1>1) then
         write(*,*) iindex(i1),i0, iindex(i1+1), f1
         stop
      endif
      if(f2<0 .or. f2>1) then
         write(*,*) jindex(j1),j0, jindex(j1+1), f2
         stop
      endif
      emiss(:)=(1-f1)*(1-f2)*emissCoarse(i1,j1,:)+ &
              (1-f1)*f2*emissCoarse(i1,j1+1,:)+ &
              f1*(1-f2)*emissCoarse(i1+1,j1,:)+ &
              f1*f2*emissCoarse(i1+1,j1+1,:)
    end subroutine interpolEmiss


    subroutine interpolPC(imemb, rhPCij, cldwPCij,cldw, rh)
      use ran_mod
      use cldclass
      implicit none
      integer :: i, j, k, imemb
      real :: dn
      real :: dx, dy, f1, f2
      real :: rhPCij(nRhEofs), cldwPCij(nCldwEofs)
      real :: cldw(nlayer), rh(nlayer)
      integer:: i1, j1, i0, j0
      i1=1
      j1=1
      i0=iRad
      j0=jRad
      do while (i0>iindex(i1))
         i1=i1+1
      end do
      do while (j0>jindex(j1))
         j1=j1+1
      end do
      if (i1 == 1) i1=2
      if (j1 == 1) j1=2
      i1=i1-1
      j1=j1-1
      
      dx=iindex(i1+1)-iindex(i1)
      dy=jindex(j1+1)-jindex(j1)
      f1=(i0-iindex(i1))/dx
      f2=(j0-jindex(j1))/dy
     
      if(f1<0 .or. f1>1) then
         write(*,*) iindex(i1),i0, iindex(i1+1), f1
         stop
      endif
      if(f2<0 .or. f2>1) then
         write(*,*) jindex(j1),j0, jindex(j1+1), f2
         stop
      endif
      rhPCij(1:nrhEOFs)=(1-f1)*(1-f2)*rhPCCoarse(i1,j1,imemb,1:nrhEOFs)+ &
           (1-f1)*f2*rhPCCoarse(i1,j1+1,imemb,1:nrhEOFs)+ &
           f1*(1-f2)*rhPCCoarse(i1+1,j1,imemb,1:nrhEOFs)+ &
           f1*f2*rhPCCoarse(i1+1,j1+1,imemb,1:nrhEOFs)
      
      cldwPCij(1:ncldwEOFs)=(1-f1)*(1-f2)*cldwPCCoarse(i1,j1,imemb,1:ncldwEOFs)+ &
           (1-f1)*f2*cldwPCCoarse(i1,j1+1,imemb,1:ncldwEOFs)+ &
           f1*(1-f2)*cldwPCCoarse(i1+1,j1,imemb,1:ncldwEOFs)+ &
           f1*f2*cldwPCCoarse(i1+1,j1+1,imemb,1:ncldwEOFs)
      cldw=cldwm
      rh=rlhm
      !if(imemb==1) then
      !   print*, rlhm
      !   stop
      !endif
      do i=1,nRhEofs
         rh(:)=rh(:)+rhPCij(i)*stdrhPC(i)*rheofs(i,:)
      enddo
      !write(*,*) rlhm
      !write(*,*) rh
      !write(*,*) rhPCij
      !write(*,*) stdrhPC
      !stop
      do i=1,nCldwEofs
         cldw(:)=cldw(:)+cldwPCij(i)*stdcldwPC(i)*cldweofs(i,:)
      enddo
      
    end subroutine interpolPC

end module geophysEns

subroutine getcldwfromcoeff(cldwPCij,hfreez,cldwprof,node,dr,localZAngle,nmemb)
  use cldclass
  use nbinMod
  use geophysens    !MG 09/25/13
  implicit none
  integer :: i, j, k, node(5), nmemb
  real :: dn, hfreez, localZAngle
  real :: dx, dy, f1, f2, dr
  real :: cldwPCij(10,nmemb), cldwprof(nbin) !MG 09/25/2013
  real :: cldw(nlayer), rh(nlayer)
  integer:: i1, j1, i0, j0
  real :: hh(nbin)
  
  cldw=cldwm
  cldwprof=0
  do i=1,nbin
     hh(i)=(nbin-i-0.5)*dr*cos(localZAngle*3.1415/180.)
  enddo
  do i=1,10
     cldw(:)=cldw(:)+sum(cldwPCij(i,1:nmemb))/nmemb*stdcldwPC(i)*cldweofs(i,:)
  enddo
  !begin  MG 9/26/13
  if(maxval(cldw)>2) then
     cldw=cldw*2/maxval(cldw)
  endif
  if(maxval(cldw)>50) then
     print*, cldw
     print*,nmemb
     print*,cldweofs(1:10,3)
     print*,stdcldwPC(1:10)
     do i=1,10
        print*,sum(cldwPCij(i,1:nmemb))/nmemb
     enddo
     !print*,'here'
     !print*, cldwPCCoarse(:,:,1:10,1:10)
     print*, 'in cldw'
     stop
  endif
!end    MG 9/26/13
  do i=node(1)+1,node(5)+1
     i0=((hh(i)+4.5-hfreez)/drrte)+1
     if(i0.gt.0 .and. i0.lt.40) then
        cldwprof(i)=cldw(i0)
     else
!begin  MG 9/26/13
        if(i0<1) cldwprof(i)=cldw(1)
        if(i0>=40) cldwprof(i)=cldw(40)
!end    MG 9/26/13
        !print*, node
        !print*, hh(i)
        !print*, drrte
        !print*, dr
        !print*, hfreez
        !stop
     endif
  enddo
  cldwprof(node(5)+2:nbin)=-99

end subroutine getcldwfromcoeff
