subroutine getplane(x1,y1,z1,x2,y2,z2,A,B,C)
  implicit none
  real:: x1,y1,z1,x2,y2,z2, A,B,C
  A=y1*z2-z1*y2
  B=z1*x2-x1*z2
  C=x1*y2-x2*y1 
  !print*, A,B,C 
end subroutine getplane

subroutine getAnglePlanes(A1,B1,C1,A2,B2,C2,angle)
  implicit none
  real:: A1,B1,C1,A2,B2,C2,angle
  real:: a11,a22
  a11=abs(A1*A2+B1*B2+C1*C2)
  a22=sqrt(A1**2+B1**2+C1**2)*sqrt(A2**2+B2**2+C2**2)
  if(a22>a11) then
     angle=acos(a11/a22)
  else
     angle=0.
  endif
end subroutine getAnglePlanes
!def hdist(lats,latf,lons,lonf):
!    dphi=abs(lats-latf)
!    dlam=abs(lons-lonf)
!    d=2*asin(sqrt(sin(dphi/2)**2+cos(lats)*cos(latf)*sin(dlam/2)**2))
!    return d
!def rotation(x,y,z,incangle):
!    x1=x
!    y1=y*cos(incangle)-z*sin(incangle)
!    z1=y*sin(incangle)+z*cos(incangle)
!    return x1,y1,z1

subroutine fromlatlon(lat,lon,x,y,z)
  implicit none
  real:: x,y,z,lat,lon
  x=cos(lon)*sin(lat)
  y=sin(lon)*sin(lat)
  z=cos(lat)
end subroutine fromlatlon

!def tolatlon(x,y,z):
!    lat=atan2(sqrt(x**2+y**2),z)
!    lon=atan2(y,x)
!    return lat,lon
!def dist(xn,yn,zn,xs,ys,zs):
!    return sqrt((xn-xs)**2+(yn-ys)**2+(zn-zs)**2)
