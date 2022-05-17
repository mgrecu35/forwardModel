 subroutine boundTbs(tb,tbout,wfmap,tableTbL,pia,n,ns,npia,ich)
  implicit none
  integer :: n,ns, ich,npia
  real, intent(in) :: tb(n,ns), wfmap(n,ns)
  real, intent(out) :: tbout(n,ns) 
  real, intent(in) :: pia(n,ns)
  real, intent(in) :: tableTbL(npia,9,2)
  integer :: i, j, ipia, ipia1, ipia2
  real :: m1Tb,m2Tb
  tbout=tb
  do i=1,n
     do j=1,ns
        if(pia(i,j)>=0) then
           ipia=pia(i,j)+1
           ipia2=ipia+1
           ipia1=ipia-1
           if(ipia2>npia) ipia2=npia
           if(ipia1>npia) ipia1=npia
           if(ipia1<1) ipia1=1
           if(ipia2<1) ipia2=1
           if(wfmap(i,j)>0.5) then
              m1Tb=max(tableTbl(ipia2,ich,1),tableTbl(ipia1,ich,1))
              if(tb(i,j)>m1Tb+5) tbout(i,j)=m1Tb+5
              m2Tb=min(tableTbl(ipia2,ich,2),tableTbl(ipia1,ich,2))
              if(tb(i,j)<m2Tb-5) tbout(i,j)=m2Tb-5
           else
              !if(tb(i,j)>290) tbout(i,j)=290
              !if(tb(i,j)<260) tbout(i,j)=260
           endif
        endif
     end do
  end do

end subroutine boundTbs

subroutine dsphere(lats, latf, lons, lonf, d)
  implicit none
  
  real lats, latf, lons, lonf
  real lonsr,lonfr, dlat, dlon, d, pi, latsr, latfr
  real d2rad, y
  pi=3.1415927
  d2rad=0.017453294 !1./180*pi
  latsr=lats*d2rad
  latfr=latf*d2rad
  lonfr=lonf*d2rad
  lonsr=lons*d2rad
  dlat=(latsr-latfr)/2
  dlon=(lonsr-lonfr)/2
  y=(sin(dlat)**2+ &
       cos(latsr)*cos(latfr)*sin(dlon)**2)**0.5
  d=2*asin(y)
  d=d*6356.
  
  
end subroutine dsphere
subroutine convifreq(actOb,&
     w,ifreq,tb,tbObs,tbout,dfdtb,n,ns,lat,lon,scLon,scLat,&
     wfmap,fpmap,fobj,sfcRain,ialg)
  Use BMCVparameters
  
  implicit none
  integer :: n,nfreq,ns
  real    :: tb(n,ns)
  real,intent(out)    :: tbout(n,ns), dfdtb(n,ns)
  integer :: i, j,dnx, dny, ifreq, ialg
  integer :: ll, kk, k
  real    :: FWHMx, FWHMy, tbconv
  real, intent(in)::  lat(n,ns), lon(n,ns), &
       scLon(n,ns), scLat(n,ns), tbObs(n,ns), wfmap(n,ns), fpmap(n,ns), &
       sfcRain(n,ns),w
  real :: fobj, lambda, tbm, ic
  integer :: actOb(49,300)
  dfdtb=0
  lambda=0.01
  FWHMx=fov_ct_microw(ifreq)
  FWHMY=fov_dt_microw(ifreq)
  if(ialg==2) then
     FWHMx=fov_ct_microw(ifreq)*2
     FWHMY=fov_dt_microw(ifreq)*2
  endif
  dny = fov_dt_microw(ifreq)/5.+1
  dnx = fov_ct_microw(ifreq)/5.+1
!  print*, FWHMx, FWHMy
!  print*, dnx, dny, n, ns
  fobj=0

!!$OMP PARALLEL DO  private(tbconv)
  do j=1,n
     if(j>dny .and. j<=n-dny) then
        do i=1,ns
           if(i>dnx .and. i<=ns-dnx .and. fpmap(j,i)>0.95 ) then
              if(tbObs(j,i)>50 .and. tbObs(j,i)<300 .and. tb(j,i)>0 ) then
                 call convPixel(w,i, j, dnx, dny, &
                      FWHMx, FWHMy, tbconv, &
                      tb,dfdtb,tbObs,n,ns,lat,lon,scLon(j,i),scLat(j,i))
                 tbout(j,i)=tbconv
                 if(tbconv>0) then
                    fobj=fobj+1/2.*w*(tbObs(j,i)-tbconv)**2
                 endif
              endif
           else
              tbout(j,i)=-99
           endif
        enddo
     else
        tbout(j,:)=-99
     endif
  enddo
!!$OMP END PARALLEL DO  
  
  
end subroutine convifreq




subroutine convPixel(w,i, j, dnx, dny, &
     FWHMx, FWHMy, tbconv, &
     tb,dfdtb,tbObs,n,ns,lat,lon,scLon,scLat)
  Use BMCVparameters

  implicit none
  integer :: n,ns
  real    :: tb(n,ns)
  integer :: i, j,dnx, dny, ifreq
  integer :: ll, kk
  real    :: FWHMx, FWHMy, tbconv
  real    :: sumw, weightc, y2, sumt, ddx, ddy
  integer :: iflag, imemb, nmemb
  real, intent(in)::  lat(n,ns), lon(n,ns), tbObs(n,ns), scLon, scLat
  real, intent(out):: dfdtb(n,ns)
  real :: xsc,ysc,zsc, xpr, ypr, zpr, xpix, ypix, zpix
  real :: A1,B1,C1, A2, B2, C2, angle, ddr
  real :: torad, w
  integer :: ipix, indi(500), indj(500), wij(500),k
  tbconv=0
  sumw=0.
  sumt=0.
  torad=1./180.*3.141592653589793

  call fromlatlon(lat(j,i)*torad,lon(j,i)*torad,xpr,ypr,zpr)
  call fromlatlon(scLon*torad,scLat*torad,xsc,ysc,zsc)

  call getplane(xsc,ysc,zsc,xpr,ypr,zpr,A1,B1,C1)
  ipix=0
  
  do ll=-dny*1,3*dny 
     do kk=-dnx*1,3*dnx 
        if(i-dnx+kk>=1 .and. i-dnx+kk<=ns .and. j-dny+ll>=1 .and. &
             j-dny+ll<=n) then

           call fromlatlon(lat(j+ll-dny,i+kk-dnx)*torad,&
                lon(j+ll-dny,i+kk-dnx)*torad,&
                xpix,ypix,zpix)

           call dsphere(lat(j,i), lat(j+ll-dny,i+kk-dnx), lon(j,i), &
                lon(j+ll-dny,i+kk-dnx), ddr)
           if(ddr>1e-3) then
              call getplane(xpr,ypr,zpr,xpix,ypix,zpix,A2,B2,C2)
              call getAnglePlanes(A1,B1,C1,A2,B2,C2,angle)
              ddx=ddr*sin(angle)
              ddy=ddr*cos(angle)
           else
              ddx=0
              ddy=0.
           endif
           y2=(((ddx/FWHMx)**2+(ddy/FWHMy)**2)* &
                &                    4*log(2.))
           

           weightc=exp(-y2)


           if(tb(j-dny+ll,i-dnx+kk)>0) then
              tbconv=tbconv+&
                   tb(j-dny+ll,i-dnx+kk)*weightc
              sumw=sumw+weightc
              if(weightc>1e-4) then
                 ipix=ipix+1
                 indi(ipix)=i-dnx+kk
                 indj(ipix)=j-dny+ll
                 wij(ipix)=weightc
              endif
           endif
           sumt=sumt+weightc
        endif
     enddo
  enddo


  if(sumw/sumt>0.95 .and. sumt> .1) then
     
     if(tbObs(j,i)>50 .and. tbObs(j,i)<300) then
        tbconv=tbconv/sumw
        do k=1,ipix
           dfdtb(indj(k),indi(k))=dfdtb(indj(k),indi(k))+&
                wij(k)/sumw*(tbconv-tbObs(j,i))*w
        enddo
       
     else
        tbconv=-99
     endif
  else
     !!print*, sumw, sumt, j, i, dnx, dny, FWHMx, FWHMy 
     tbconv=-99
  endif
  !print*, ipix
end subroutine convPixel

subroutine footprintC(i, j, dnx, dny, &
     FWHMx, FWHMy, wfmap, fpmap,n,ns,lat,lon,scLon,scLat)
  Use BMCVparameters

  implicit none
  integer :: n,ns
  integer :: i, j,dnx, dny, ifreq
  integer :: ll, kk
  real    :: FWHMx, FWHMy, tbconv
  real    :: sumw, weightc, y2, sumt, ddx, ddy
  integer :: iflag, imemb, nmemb
  real, intent(in)::  lat(n,ns), lon(n,ns), wfmap(n,ns), scLon, &
     scLat
  real, intent(out):: fpmap(n,ns)
  real :: xsc,ysc,zsc, xpr, ypr, zpr, xpix, ypix, zpix
  real :: A1,B1,C1, A2, B2, C2, angle, ddr
  real :: torad
  integer :: ipix, indi(500), indj(500), wij(500),k
  tbconv=0
  sumw=0.
  sumt=0.
  torad=1./180.*3.141592653589793
 ! print*, lat(j,i), lon(j,i)
 ! print*, scLat, scLon

  call fromlatlon(lat(j,i)*torad,lon(j,i)*torad,xpr,ypr,zpr)
  call fromlatlon(scLon*torad,scLat*torad,xsc,ysc,zsc)

  call getplane(xsc,ysc,zsc,xpr,ypr,zpr,A1,B1,C1)
  ipix=0
  fpmap(j,i)=1
  do ll=0,2*dny 
     do kk=0,2*dnx 
        if(i-dnx+kk>=1 .and. i-dnx+kk<=ns .and. j-dny+ll>=1 .and. &
             j-dny+ll<=n) then

           call fromlatlon(lat(j+ll-dny,i+kk-dnx)*torad,&
                lon(j+ll-dny,i+kk-dnx)*torad,&
                xpix,ypix,zpix)

           call dsphere(lat(j,i), lat(j+ll-dny,i+kk-dnx), lon(j,i), &
                lon(j+ll-dny,i+kk-dnx), ddr)
           !print*, ddr
           if(ddr>1e-3) then
              call getplane(xpr,ypr,zpr,xpix,ypix,zpix,A2,B2,C2)
              call getAnglePlanes(A1,B1,C1,A2,B2,C2,angle)
              ddx=ddr*sin(angle)
              ddy=ddr*cos(angle)
           else
              ddx=0
              ddy=0.
           endif
           y2=(((ddx/FWHMx)**2+(ddy/FWHMy)**2)* &
                &                    4*log(2.))
           

           weightc=exp(-y2)

           
           if(weightc>1e-3) then
              if(wfmap(j-dny+ll,i-dnx+kk)<0.95) then
                 fpmap(j,i)=0
              endif
           endif
           sumt=sumt+weightc
        endif
     enddo
  enddo


end subroutine footprintC

subroutine footprintmap(ifreq,wfmap,fpmap,n,ns,lat,lon,scLon,scLat)
  Use BMCVparameters
  
  implicit none
  integer :: n,nfreq,ns
  real    :: fpmap(n,ns)
  integer :: i, j,dnx, dny, ifreq
  integer :: ll, kk, k
  real    :: FWHMx, FWHMy, tbconv
  real    ::  lat(n,ns), lon(n,ns), &
       wfmap(n,ns), scLon(n,ns), scLat(n,ns)
  real :: fobj, lambda, tbm, ic

!  FWHMx=fov_ct_microw(ifreq)
!  FWHMY=fov_dt_microw(ifreq)
  
!  dny = fov_dt_microw(ifreq)/5.+1
!  dnx = fov_ct_microw(ifreq)/5.+1

  dnx=5
  dny=5
  fpmap=0
  !print*, dnx, dny
  !print*, n, ns
  !print*, maxval(wfmap)

  do j=1,n
     if(j>dny .and. j<=n-dny) then
        do i=1,ns
           if(i>dnx .and. i<=ns-dnx) then
              !print*, j,i
              !print*, lat(j,i), lon(j,i)
              !print*, scLon(j,i), scLat(j,i), wfmap(j,i)
              !call footprintC(i, j, dnx, dny, &
              !     FWHMx, FWHMy, wfmap, fpmap,n,ns,lat,lon,scLon(j,i),&
              !     scLat(j,i))
              !stop
           endif
        enddo
     endif
  enddo

  
end subroutine footprintmap

subroutine footprintmap2(ifreq,wfmap,fpmap,n,ns,lat,lon,scLon,scLat)
  Use BMCVparameters
  implicit none
  integer :: n,nfreq,ns
  real    :: fpmap(n,ns)
  integer :: i, j,dnx, dny, ifreq
  integer :: ll, kk, k
  real    :: FWHMx, FWHMy, tbconv
  real    ::  lat(n,ns), lon(n,ns), &
       wfmap(n,ns), scLon(n,ns), scLat(n,ns)
  real :: fobj, lambda, tbm, ic


  dny = fov_dt_microw(ifreq)/5.+1
  dnx = fov_ct_microw(ifreq)/5.+1
  !print*, n, ns
  !print*, lat(1,1), lon(1,1)
  !print*, scLon(1,1),scLat(1,1)
  FWHMx=fov_ct_microw(ifreq)
  FWHMY=fov_dt_microw(ifreq)
  fpmap=0
  do j=1,n
     if(j>dny .and. j<=n-dny) then
        do i=1,ns
           if(i>dnx .and. i<=ns-dnx) then
              !print*, j,i
              !print*, lat(j,i), lon(j,i)
              !print*, scLon(j,i), scLat(j,i), wfmap(j,i)
              call footprintC(i, j, dnx, dny, &
                   FWHMx, FWHMy, wfmap, fpmap,n,ns,lat,lon,scLon(j,i),&
                   scLat(j,i))

              !stop
           endif
        enddo
     endif
  enddo

end subroutine footprintmap2
