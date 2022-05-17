!  SFM  08/09/2013 Added routoines interpol3sst and interpol3sfmvc to
!                   correct problems with default values and to allow
!                   search for "best available"; removed some unused
!                   modules


subroutine project2tmi(tbmax,tbmin,tbmean,tbstd, &  !Oct06
     tbgridmax,tbgridmin,tbgridmean,tbgridstd, xmin, ymin, dx, & !Oct06
     n1b11,nxg,nyg) !Oct06
  implicit none !Oct06
  integer :: i, j, k, n1b11, nxg, nyg, ig, jg            !Oct06
  real :: tbgridmax(7,nyg,nxg),tbgridmin(7,nyg,nxg), &   !Oct06
       tbgridmean(7,nyg,nxg), tbgridstd(7,nyg,nxg)       !Oct06
  real :: tbmax(7,208,3*n1b11-2) !Oct06
  real :: tbmin(7,208,3*n1b11-2) !Oct06
  real :: tbmean(7,208,3*n1b11-2) !Oct06
  real :: tbstd(7,208,3*n1b11-2)  !Oct06
  real :: xmin, ymin, dx          !Oct06
  
end subroutine project2tmi                 !Oct06

subroutine interpol3vc(tmilow,tmilow3,nf,n1b11) !Oct06
  implicit none                                !Oct06
  integer :: n1b11, nf                         !Oct06
  real :: tmilow3(nf,221,3*n1b11-2), tmilow(nf,221,n1b11) !Oct06
  integer :: i, j, k                                      !Oct06
  do i=1,nf                                               !Oct06
     do k=70,150               !MG 09132013  
        do j=1,n1b11-1                                    !Oct06
           tmilow3(i,k,3*j-2)=tmilow(i,k,j)
           tmilow3(i,k,3*j-1)=2/3.*tmilow(i,k,j)+1./3.*tmilow(i,k,j+1)
           tmilow3(i,k,3*j)=1/3.*tmilow(i,k,j)+2./3.*tmilow(i,k,j+1) 
        enddo 
        tmilow3(i,k,3*n1b11-2)=tmilow(i,k,n1b11) !Oct06
     enddo                                           !Oct06
  enddo
end subroutine interpol3vc

subroutine interpol3sc(tmilow,tmilow3,landSea,n1b11)        
  implicit none                        
  integer :: n1b11                     
  real :: tmilow3(221,3*n1b11-2), tmilow(221,n1b11) , landSea(221,n1b11) 
  integer :: i, j, k                                   
  
  do k=70,150               !MG 09132013                                    
     do j=1,n1b11-1                          
        tmilow3(k,3*j-2)=tmilow(k,j)  
        if(tmilow(k,j)<-98 .and. tmilow(k,j+1)<-98 ) then
           tmilow3(k,3*j-1)=-99
           tmilow3(k,3*j)  =-99
        else
           if(tmilow(k,j)< -98 ) then
              tmilow3(k,3*j-1)=tmilow(k,j+1) 
              tmilow3(k,3*j)=tmilow(k,j+1)  
           else
              if(tmilow(k,j+1)<-98) then
                 tmilow3(k,3*j-1)=tmilow(k,j) 
                 tmilow3(k,3*j)=tmilow(k,j)
              else
                 tmilow3(k,3*j-1)=2/3.*tmilow(k,j)+1./3.*tmilow(k,j+1) 
                 tmilow3(k,3*j)=1/3.*tmilow(k,j)+2./3.*tmilow(k,j+1) 
                 if(tmilow3(k,3*j-1)<0) then
                    print*, tmilow3(k,3*j-1),tmilow3(k,3*j), tmilow(k,j), &
		            tmilow(k,j+1), k, j
                    !print*, n1b11
                    print*, landSea(k,j),landSea(k,j+1)
                    stop 'STOP a FROM interpol3sc'
                 endif
                 if(tmilow3(k,3*j)<0) then
                    print*, tmilow3(k,3*j-1),tmilow3(k,3*j) 
                    print*, tmilow3(k,3*j-1),tmilow3(k,3*j) 
                    print*, landSea(k,j),landSea(k,j+1)
                    stop 'STOP b FROM interpol3sc'
                 endif
              endif
           endif
        endif
     enddo
  enddo
end subroutine interpol3sc


subroutine interpol3sfmvc(tmilow,tmilow3,nf,n1b11)        

  implicit none                        
  integer :: n1b11,nf                
  real :: tmilow3(nf,221,3*n1b11-2), tmilow(nf,221,n1b11)
  integer :: i, j, k                                   

  INTEGER*4 nn,mm,low_n,high_n,low_m,high_m, nsumm
  REAL*4    summ, substitute
  
  DO i=1,nf 
  do k=70,150               !MG 09132013                                    
     do j=1,n1b11-1                          
        tmilow3(i,k,3*j-2)=tmilow(i,k,j)  
        if(tmilow(i,k,j)<-98 .and. tmilow(i,k,j+1)<-98 ) then
           tmilow3(i,k,3*j-1)=-999
           tmilow3(i,k,3*j)  =-999
        else
           if(tmilow(i,k,j)< -98 ) then
              tmilow3(i,k,3*j-1)=tmilow(i,k,j+1) 
              tmilow3(i,k,3*j)=tmilow(i,k,j+1)  
           else
              if(tmilow(i,k,j+1)<-98) then
                 tmilow3(i,k,3*j-1)=tmilow(i,k,j) 
                 tmilow3(i,k,3*j)=tmilow(i,k,j)
              else
                 tmilow3(i,k,3*j-1)=2/3.*tmilow(i,k,j)+1./3.*tmilow(i,k,j+1) 
                 tmilow3(i,k,3*j)=1/3.*tmilow(i,k,j)+2./3.*tmilow(i,k,j+1) 
                 IF(tmilow3(i,k,3*j-1)<0) THEN  ! seek nearby mn substitute in 5x5
                    low_m = MAX(j-2,1)
                    high_m= MIN(j+2,n1b11)
                    low_n = MAX(k-2,1)
                    high_n= MIN(k+2,221)
                    summ = 0.0
                    nsumm = 0
                    DO nn=low_n,high_n
                       DO mm=low_m,high_m
                          IF (tmilow(i,nn,mm) .GT. 0.10 .AND.                  &
			      tmilow(i,nn,mm) .LT. 500.0)                      &
			  THEN
                             summ = summ + tmilow(i,nn,mm)
                             nsumm = nsumm + 1
                          ENDIF
                       ENDDO
                    ENDDO
                    IF (nsumm .GE. 1) THEN
                       substitute = summ / nsumm
                       tmilow3(i,k,3*j-1) = substitute
                       tmilow3(i,k,3*j)   = substitute
                    ELSE
                       tmilow3(i,k,3*j-1) = 300.0
                       tmilow3(i,k,3*j)   = 300.0
                    ENDIF
                 ENDIF
                 IF (tmilow3(i,k,3*j)<0) THEN
                    low_m = MAX(j-2,1)
                    high_m= MIN(j+2,n1b11)
                    low_n = MAX(k-2,1)
                    high_n= MIN(k+2,221)
                    summ = 0.0
                    nsumm = 0
                    DO nn=low_n,high_n
                       DO mm=low_m,high_m
                          IF (tmilow(i,nn,mm) .GT. 0.10 .AND.                  &
			      tmilow(i,nn,mm) .LT. 500.0)                      &
			  THEN
                             summ = summ + tmilow(i,nn,mm)
                             nsumm = nsumm + 1
                          ENDIF
                       ENDDO
                    ENDDO
                    IF (nsumm .GE. 1) THEN  ! seek nearby mn substitute in 5x5
                       substitute = summ / nsumm
                       tmilow3(i,k,3*j-1) = substitute
                       tmilow3(i,k,3*j)   = substitute
                    ELSE
                       tmilow3(i,k,3*j-1) = 300.0
                       tmilow3(i,k,3*j)   = 300.0
                    ENDIF
                 endif
              endif
           endif
        endif
     enddo
     tmilow3(i,k,3*n1b11-2)=tmilow(i,k,n1b11)
  enddo
  ENDDO

end subroutine interpol3sfmvc


subroutine interpol3sst(tmilow,tmilow3,landSea,n1b11)        
!
!  SFM  08/08/2013  adapted from interpol3sc, specifically for sst 
!                   (sea surface temperature); adds provisions for seeking
!                   substitute value fromm mean of nearby 5x5 grid

  implicit none                        
  integer :: n1b11                     
  real :: tmilow3(221,3*n1b11-2), tmilow(221,n1b11) , landSea(221,n1b11) 
  integer :: i, j, k                                   

  INTEGER*4 nn,mm,low_n,high_n,low_m,high_m, nsumm
  REAL*4    summ, substitute
  
  do k=70,150               !MG 09132013                                    
     do j=1,n1b11-1                          
        tmilow3(k,3*j-2)=tmilow(k,j)  
        if(tmilow(k,j)<-98 .and. tmilow(k,j+1)<-98 ) then
           tmilow3(k,3*j-1)=-999
           tmilow3(k,3*j)  =-999
        else
           if(tmilow(k,j)< -98 ) then
              tmilow3(k,3*j-1)=tmilow(k,j+1) 
              tmilow3(k,3*j)=tmilow(k,j+1)  
           else
              if(tmilow(k,j+1)<-98) then
                 tmilow3(k,3*j-1)=tmilow(k,j) 
                 tmilow3(k,3*j)=tmilow(k,j)
              else
                 tmilow3(k,3*j-1)=2/3.*tmilow(k,j)+1./3.*tmilow(k,j+1) 
                 tmilow3(k,3*j)=1/3.*tmilow(k,j)+2./3.*tmilow(k,j+1) 
                 IF(tmilow3(k,3*j-1)<0) THEN  ! seek nearby mn substitute in 5x5
                    low_m = MAX(j-2,1)
                    high_m= MIN(j+2,n1b11)
                    low_n = MAX(k-2,1)
                    high_n= MIN(k+2,221)
                    summ = 0.0
                    nsumm = 0
                    DO nn=low_n,high_n
                       DO mm=low_m,high_m
                          IF (tmilow(nn,mm) .GT. 0.10 .AND.                    &
			      tmilow(nn,mm) .LT. 500.0)                        &
			  THEN
                             summ = summ + tmilow(nn,mm)
                             nsumm = nsumm + 1
                          ENDIF
                       ENDDO
                    ENDDO
                    IF (nsumm .GE. 1) THEN
                       substitute = summ / nsumm
                       tmilow3(k,3*j-1) = substitute
                       tmilow3(k,3*j)   = substitute
                    ELSE
                       tmilow3(k,3*j-1) = 20.0
                       tmilow3(k,3*j)   = 20.0
                    ENDIF
                 ENDIF
                 IF (tmilow3(k,3*j)<0) THEN
                    low_m = MAX(j-2,1)
                    high_m= MIN(j+2,n1b11)
                    low_n = MAX(k-2,1)
                    high_n= MIN(k+2,221)
                    summ = 0.0
                    nsumm = 0
                    DO nn=low_n,high_n
                       DO mm=low_m,high_m
                          IF (tmilow(nn,mm) .GT. 0.10 .AND.                    &
			      tmilow(nn,mm) .LT. 500.0)                        &
			  THEN
                             summ = summ + tmilow(nn,mm)
                             nsumm = nsumm + 1
                          ENDIF
                       ENDDO
                    ENDDO
                    IF (nsumm .GE. 1) THEN  ! seek nearby mn substitute in 5x5
                       substitute = summ / nsumm
                       tmilow3(k,3*j-1) = substitute
                       tmilow3(k,3*j)   = substitute
                    ELSE
                       tmilow3(k,3*j-1) = 20.0
                       tmilow3(k,3*j)   = 20.0
                    ENDIF
                 endif
              endif
           endif
        endif
     enddo
  enddo

end subroutine interpol3sst

subroutine interpol3x(txlon,txlon3,n1b11)            !Oct06
  implicit none                                      !Oct06
  integer :: n1b11                                   !Oct06
  real :: txlon3(208,3*n1b11-2), txlon(208,n1b11)    !Oct06
  integer :: i, j, k                                 !Oct06
  
  do k=1,208                                         !Oct06
     do j=1,n1b11-1                                  !Oct06
        txlon3(k,3*j-2)=txlon(k,j)                   !Oct06
        txlon3(k,3*j-1)=2./3.*txlon(k,j)+1./3.*txlon(k,j+1) !Oct06
        txlon3(k,3*j)=1./3.*txlon(k,j)+2./3.*txlon(k,j+1)   !Oct06
     enddo    
     txlon3(k,3*n1b11-2)=txlon(k,n1b11)                   !Oct06
  enddo                                                   !Oct06
end subroutine interpol3x                                 !Oct06

subroutine interpol3sc2(txlon,txlon3,n1b11)            !Oct06
  implicit none                                      !Oct06
  integer :: n1b11                                   !Oct06
  real :: txlon3(3*n1b11-2), txlon(n1b11)    !Oct06
  integer :: i, j, k                                 !Oct06
  

  do j=1,n1b11-1                                  !Oct06
     txlon3(3*j-2)=txlon(j)                   !Oct06
     if(txlon(j)>-998 .and. txlon(j+1)>-998) then !4/15/14 MG Begin
        txlon3(3*j-1)=2./3.*txlon(j)+1./3.*txlon(j+1) !Oct06
        txlon3(3*j)=1./3.*txlon(j)+2./3.*(txlon(j+1))
     else
        txlon3(3*j)=-999.9
        txlon3(3*j-1)=-999.9
     endif !4/15/14 MG End
  enddo
  txlon3(3*n1b11-2)=txlon(n1b11)                   !Oct06
  
  
end subroutine interpol3sc2

subroutine interpol3xc(txlon,txlon3,n1b11)            !Oct06
  implicit none                                      !Oct06
  integer :: n1b11                                   !Oct06
  real :: txlon3(221,3*n1b11-2), txlon(221,n1b11)    !Oct06
  integer :: i, j, k                                 !Oct06
  
  do k=70,150               !MG 09132013
     do j=1,n1b11-1                                  !Oct06
        txlon3(k,3*j-2)=txlon(k,j)                   !Oct06
        if(txlon(k,j)>-998 .and. txlon(k,j+1)>-998) then !4/15/14 MG Begin
           txlon3(k,3*j-1)=2./3.*txlon(k,j)+1./3.*txlon(k,j+1) !Oct06
           txlon3(k,3*j)=1./3.*txlon(k,j)+2./3.*(txlon(k,j+1))
        else
           txlon3(k,3*j)=-999.9
           txlon3(k,3*j-1)=-999.9
        endif !4/15/14 MG End
     enddo
     txlon3(k,3*n1b11-2)=txlon(k,n1b11)                   !Oct06
  enddo                                                   !Oct06
end subroutine interpol3xc                                 !Oct06

subroutine getgridind(xlon2a25,ylat2a25,indig,indjg,xmin,ymin,dx,n1c21,nxg,nyg)
implicit none
integer :: nxg,nyg,n1c21
integer :: indig(nxg,nyg), indjg(nxg,nyg)
real :: xlon2a25(49,n1c21), ylat2a25(49,n1c21)
real :: xmin,ymin,dx
integer :: ig, jg, i, j
indig=-99
indjg=-99
do i=1,49
   do j=1,n1c21
      ig=(xlon2a25(i,j)-xmin)/dx+1
      jg=(ylat2a25(i,j)-ymin)/dx+1
      indig(ig,jg)=i
      indjg(ig,jg)=j
   enddo
enddo
end subroutine getgridind

subroutine getPRind(xlon,ylat,indig,indjg,xmin,ymin,dx,indi,indj,n1b11,nxg,nyg)
implicit none
integer :: nxg,nyg,n1b11
integer :: indig(nxg,nyg), indjg(nxg,nyg)
real :: xlon(208,3*n1b11-2), ylat(208,3*n1b11-2)
integer :: indi(208,3*n1b11-2), indj(208,3*n1b11-2)
real :: xmin,ymin,dx
integer :: ig, jg, i, j
indi=0
indj=0
open(10,file='igjg.dat')
do i=1,208
   do j=1,3*n1b11-2
      ig=(xlon(i,j)-xmin)/dx+1
      jg=(ylat(i,j)-ymin)/dx+1
      if(ig>0 .and. ig<=nxg  .and. jg>0 .and.jg<=nyg) then
         indi(i,j)=indig(ig,jg)
         indj(i,j)=indjg(ig,jg)
      endif
      write(10,*) xlon(i,j), ylat(i,j), indi(i,j), indj(i,j)
   enddo
   write(10,*)
enddo
close(10)
!stop
end subroutine getPRind
