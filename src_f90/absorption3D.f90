! GasabsR98(F,Tk,Rhowv,Pa,absair,abswv,ireturn)
subroutine absorption3D(nz,ny,nx,T,qv,qc,rho,press,abs_air,abs_wv,abs_clw,f)
  implicit none
  integer :: nx,ny,nz
  real :: T(nz,ny,nx),qv(nz,ny,nx),press(nz,ny,nx),rho(nz,ny,nx),qc(nz,ny,nx)
  real, intent(out) :: abs_air(nz,ny,nx), abs_wv(nz,ny,nx), abs_clw(nz,ny,nx)
  integer ireturn,i,j,k
  real  :: f

  do k=1,nz
     do j=1,ny
        do i=1,nx
           ireturn=0
           call gasabsr98(f,T(k,j,i),qv(k,j,i)*rho(k,j,i),press(k,j,i),abs_air(k,j,i),abs_wv(k,j,i),ireturn)
           call gcloud(f,T(k,j,i),qc(k,j,i)*1e3*rho(k,j,i),abs_clw(k,j,i))
        enddo
     enddo
  enddo

end subroutine absorption3D


subroutine radtran3D(nz,ny,nx,T,tsk,kext,salb,asym,emiss_2d,height,tb2d,umu)
  implicit none
  integer :: nx,ny,nz
  real :: T(nz,ny,nx), tsk(ny,nx), height(nz+1,ny,nx),umu, emiss_2d(2,ny,nx)
  !real :: tsk(ny,nx)
  real :: kext(nz,ny,nx),salb(nz,ny,nx),asym(nz,ny,nx)
  real :: tlayer(0:nz)
  real, intent(out) :: tb2d (2,ny,nx)
  integer i,j,k
  real  :: tb1, fisot
  logical lambert
  lambert=.true.
  fisot=2.7
  do j=1,ny
     do i=1,nx
        tlayer(0)=tsk(j,i)
        do k=1,nz-1
           tlayer(k)=0.5*(t(k,j,i)+t(k+1,j,i))
        enddo
        tlayer(nz)=t(nz,j,i)
        call radtran(umu,nz,tb1,tsk(j,i),tLayer,height(:,j,i), &
             kext(:,j,i),salb(:,j,i),asym(:,j,i), &
             fisot,emiss_2d(1,j,i),emiss_2d(1,j,i),lambert)
        tb2d(1,j,i)=tb1
        call radtran(umu,nz,tb1,tsk(j,i),tLayer,height(:,j,i), &
             kext(:,j,i),salb(:,j,i),asym(:,j,i), &
             fisot,emiss_2d(2,j,i),emiss_2d(2,j,i),lambert)
        tb2d(2,j,i)=tb1
     enddo
  enddo
end subroutine radtran3D
