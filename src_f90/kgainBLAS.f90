!  SFM  04/06/2013  Code changes from M.Grecu
!  SFM  05/06/2013  Code changes from M.Grecu, via L.Woltz
!
subroutine covYYp(covyy,yy,ny,nm)
integer ny, nm
real ::   yy(ny,nm), covyy(ny,ny)
real ::   my(ny)

do i=1, ny
   my(i)=sum(yy(i,1:nm))/nm
enddo

do i=1, ny
   do j=1, ny
      covyy(i,j)=sum((yy(i,1:nm)-my(i))*(yy(j,1:nm)-my(j)))/nm
   enddo
enddo

end subroutine covYYp

subroutine covXYp(covxy,yy,xx,ny,nx,nm)
integer ny, nm, nx
real ::   yy(ny,nm), xx(nx,nm), covxy(nx,ny)
real ::   mx(nx), my(ny)

do i=1, nx
   mx(i)=sum(xx(i,1:nm))/nm
enddo
do i=1, ny
   my(i)=sum(yy(i,1:nm))/nm
enddo

do i=1, nx
   do j=1, ny
      covxy(i,j)=sum((xx(i,1:nm)-mx(i))*(yy(j,1:nm)-my(j)))/nm
   enddo
enddo

end subroutine covXYp
!enkF1d(Xens(:,1:nmemb1-1),Yens(:,1:nmemb1-1),Yobs,nx,ny,nMemb-1,xup)
subroutine enkF1d(xens,yens,yobs,nx,ny,nm,xup, ipias,s0KuVar,s0KaVar,S0Cov,&
  piaKuVar,piaKaVar)
  use ran_mod
  implicit none
  integer :: ny,nx,nm, ipias(2)
  real :: xens(nx,nm), yens(ny,nm), yobs(ny), ym(ny), ydiff(ny,nm+1)
  real :: covyy(ny,ny), covxy(nx,ny)
  character :: uplo
  integer   :: info, lda, ldb, n, nrhs, i, j, irec, ii
  integer   :: nmemb1, nmemb2,  k
  real      :: xm(nx), xup(nx), maxx, minx, stdx, s1, s0KuVar,S0KaVar,s0Cov,&
       piaKuVar,piaKaVar

! WSO  begin 08/20/2015 increased array dimension to 800 to avoid overflow for
! "deep"
! profiles
  real      :: dxup(800)
! WSO  end 08/20/2015

  !print*, nx, ny, nm
  !print*, Yobs

  nmemb1=nm
!   print '(11F10.4)', yens(1:ny,1:nm)
  call covYYp(covyy(1:ny,1:ny),yens(1:ny,1:nm),ny,nmemb1)
!   print '(11F10.4)', covyy(1:ny,1:ny)
!   print*, ny
!   stop

  do i=1,ipias(1)-1 !Ka dBZ obs
     if(i.ne.ipias(1)) then
        covyy(i,i)=covyy(i,i)+4.
     else
        covyy(i,i)=covyy(i,i)+4.
     endif
  enddo
  
  !   print '(11F10.4)', covyy(1:ny,1:ny)
!   stop
  call covXYp(covxy(1:nx,1:ny),yens(1:ny,1:nm),xens(1:nx,1:nm),ny,nx,nm)
  
  uplo='U'

  do i=1,ny
     ym(i)=sum(yens(i,1:nm))/nm
     ydiff(i,1)=yobs(i)-ym(i)
     do j=1,nm
        ydiff(i,1+j)=yobs(i)-yens(i,j)
     enddo
  enddo

  
  covyy(ipias(1),ipias(1))=covyy(ipias(1),ipias(1))+s0KuVar
  covyy(ipias(2),ipias(2))=covyy(ipias(2),ipias(2))+s0KaVar
  covyy(ipias(1)+2,ipias(1)+2)=covyy(ipias(1)+2,ipias(1)+2)+piaKuVar
  covyy(ipias(1)+3,ipias(1)+3)=covyy(ipias(1)+3,ipias(1)+3)+piaKaVar

  do i=ipias(2)+3,ny
     covyy(i,i)=covyy(i,i)+100/2./2.
  enddo
  do i=ipias(2)+10,ny
     covyy(i,i)=covyy(i,i)+150/2./2.
  enddo

  call sposv(uplo, ny, 1+nm, covyy(1:ny,1:ny), ny, ydiff(1:ny,1:1+nm), ny, info)

  s1=stdx(xens(1,1:nm),nm)
  !print*, info
  if(info .ne. 0) then
     !print*, 'info1',info
     do i=1,-ny
        write(*,10) covyy(i,1:ny)
     enddo
     !stop
     10 format(70(F8.2,1x))
     !print*, s0kuVar,s0KaVar
     !print*, yens(1,1:nm)
     !print*, yens(2,1:nm)
  else
!     print*, 'info0',info 
  endif
  do j=1,nx
     xm(j)=sum(xens(j,1:nm))/nm
     maxx=maxval(xens(j,1:nm))
     minx=minval(xens(j,1:nm))
     xup(j)=xm(j)
  enddo
  if (maxval(yobs(1:ny)).gt.1e3) then
     return
  endif
  if (minval(yobs(1:ny)).lt.-1e3) then
     return
  endif

  if(info==0) then
!! begin MG 09132013  ! GPM instruments left for Japan
     do j=1,nx
        dxup(j)=0
        do i=1,ny
           dxup(j)=dxup(j)+covxy(j,i)*ydiff(i,1)
        enddo
        xm(j)=sum(xens(j,1:nm))/nm

     enddo
!! end MG 09182013     
     !print*,  xm(1), dxup
     !print*, nx
!! begin MG 09182013
     !print*, xm(1)+dxup(1)
     !print*, isnan(dxup(1))
     if(xm(1)+dxup(1) .lt. 300 .and. (isnan(dxup(1)) .eqv. .false.)) then
 !       print*, 'update'
        do j=1,nx
           xm(j)=sum(xens(j,1:nm))/nm
           maxx=maxval(xens(j,1:nm))
           minx=minval(xens(j,1:nm))
           xup(j)=xm(j)
           
           do i=1,ny
              xup(j)=xup(j)+covxy(j,i)*ydiff(i,1)
              !print*, j, i, covxy(j,i), ydiff(i,1), xup(j)
              !        if(xup(j)>maxx) xup(j)=maxx
              do k=1,nm
                 xens(j,k)=xens(j,k)+covxy(j,i)*ydiff(i,1+k)
              enddo
           enddo
           if(xup(j)<0.5*minx) xup(j)=0.5*minx
           do k=1,nm
              if(minx>0) then
                 if(xens(j,k)<0.2*minx) xens(j,k)=0.2*minx
              endif
              if(maxx>0) then
                 if(xens(j,k)>1.2*maxx) xens(j,k)=1.2*maxx
              endif
           enddo
        enddo
     endif
     
!     xup(1:nx)=xm(1:nx)+dxup(1:nx)
!! end MG 09182013                
     if(xup(1) .gt. 1110 .or. isnan(xup(1))) then
        print*, xens(1,1:nm)
        print*, yens(1:ny,1)
        print*, yobs(1:ny)
        print*, nm,ny,nx
        print*, 'Unc=',s1,stdx(xens(1,1:nm),nm)
	print *, 'COVXY',(covxy(1,ii),ii=1,ny),'   ',xup(1),isnan(xup(1))
        print*, ydiff(1:ny,1)
        !print*, dxup(1), xm(1)
        !print*, dxup
        print*, 'kgain'
        stop
     end if
  endif
end subroutine enkF1d

function  stdx(xn,n)
integer :: n
real :: stdx
real :: xn(n)

xm=sum(xn)/n
stdx=sqrt(sum((xn-xm)**2/(n-1)))
end function stdx
