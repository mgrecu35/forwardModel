subroutine asciiplot(imag,nx,ny,dny,ilog,imin,imax)
real:: imag(nx,ny), imin,imax
character*16 ascii_chars
character*16 ascii_chars2


character*100 line
integer ::  ilev, dny, ilog
ascii_chars="MNHQ$OC?7>!:-;. "
ascii_chars2=".-+*$#@&"
line=''//ascii_chars(1:1)//ascii_chars(2:2)
line=line(1:2)//ascii_chars(3:3)
!print*, nx, ny, dny, ilog, imin, imax
!return
if (ilog==1) then
   do i=1,nx
      line=''
      ij=0
      do j=1,ny,dny
         ij=ij+1
         ilev=1
         if(imag(i,j)>0.) then
            ilev=(8.*(log(imag(i,j))-log(imin))/(log(imax)-log(imin)))+1
         end if
         if(ilev>8) ilev=8
         if(ilev<1) ilev=1
         i1=ilev
         if(ij==1) then
            line=ascii_chars2(i1:i1)
         else
            line(1:ij)=line(1:ij-1)//ascii_chars2(i1:i1)
         endif
      enddo
      write(*,101) line
   enddo
101 format(100A)
endif


end subroutine asciiplot
