subroutine background(tb,tbmean,invcovTb,dfdtb,n,ns,nf,fobj,wfmap)
implicit none
  integer :: n,nfreq,ns,nf
  real    :: tb(n,ns,nf), tbmean(n,ns,nf), invcovTb(n,ns,nf,nf), dfdtb(n,ns,nf)
  real    :: wfmap(n,ns)
  real    :: fobj, fobjold
  integer :: i, j, k, l,i1,j1
  fobjold=fobj
  !invCovTb=0.
  do i=1,n
     do j=1,ns
        do k=1,nf
           do l=1,nf
              if(tbmean(i,j,k)>0 .and. tbmean(i,j,l)>0 .and. &
                   wfmap(i,j)>0.9) then
                 fobj=fobj+0.5*invcovTb(i,j,k,l)*(tb(i,j,k)-tbmean(i,j,k))*&
                      (tb(i,j,l)-tbmean(i,j,l))
                 dfdtb(i,j,k)=dfdtb(i,j,k)+0.5*invcovTb(i,j,k,l)*&
                      (tb(i,j,l)-tbmean(i,j,l))
                 dfdtb(i,j,l)=dfdtb(i,j,l)+0.5*invcovTb(i,j,k,l)*&
                      (tb(i,j,k)-tbmean(i,j,k))
              endif
           enddo
        enddo
     enddo
  enddo

110 format(9(F9.4,1x))
end subroutine background

subroutine convallfreq(actOb,tb,tbmean,invCovTb,&
     tbObs,tbout,dfdtb,n,ns,lat,lon,scLon,scLat,&
     wfmap,fpmap,nf,fobj,ifreqs,sfcRain,ialg)
  implicit none

  integer :: n,nfreq,ns,nf,ifreqs(nf)
  real    :: tb(n,ns,nf),tbmean(n,ns,nf), sfcRain(n,ns)
  real    :: invCovTb(n,ns,nf,nf)

  real,intent(out)    :: tbout(n,ns,nf), dfdtb(n,ns,nf)
 
  real, intent(in)::  lat(n,ns), lon(n,ns), &
       scLon(n,ns), scLat(n,ns), tbObs(n,ns,nf), wfmap(n,ns), fpmap(n,ns,nf)
  real :: fobj, lambda, tbm, ic, w(9)
  integer :: i, ialg
  real   :: fobj1
  integer :: actOb(49,300)
  w=(/1.,1.,1.,1.,0.7,0.3,0.3,0.1,0.1/)
  fobj=0
  do i=1,9
     call convifreq(actOb,w(i),ifreqs(i),tb(:,:,i),tbObs(:,:,i),tbout(:,:,i),&
          dfdtb(:,:,i),n,ns,lat,lon,scLon,scLat,&
          wfmap,fpmap(:,:,i),fobj1,sfcRain,ialg)
     !print*, fobj1, i
     fobj=fobj+fobj1
     !print*, maxval(tb(:,:,i)),minval(tb(:,:,i))
  enddo
  !print*, 'Fobj1=',fobj
  !call background(tb,tbmean,invcovTb,dfdtb,n,ns,nf,fobj,wfmap)
  !print*, 'Fobj2=',fobj
end subroutine convallfreq
