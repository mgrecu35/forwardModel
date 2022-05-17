module gEnv
  
  real:: qvenvG(88)
  real:: tempG(88)
  real:: pressEnvG(88), sfcTempEnvG
  real:: wvExt(88,8)
  real:: tb4emiss(8,4)
!begin  MG 10/29/15 add reliabFlag
  integer :: reliabFlag
!end    MG 10/29/15
end module gEnv

real function  getextka(i)
  use gEnv
  getextka=wvExt(i+1,3)
end function getextka
subroutine setEnv(qvEnv,tEnv,pressEnv,sfcTempEnv,sknTempEnv)
  use gEnv
  implicit none
  real :: qvEnv(88), tEnv(88), pressEnv(88), sfcTempEnv,sknTempEnv
  real :: freqs(8), absair, abswv, rhowv, rho
  integer :: i, j, ireturn, isfc
  real :: dboux,umu
  sfcTempEnvG=sfcTempEnv
  freqs=(/10.000000,19.000000,22.000000,37.000000,85.000000,166.,186.3,190.3/)
  !if(minval(tenv)<100) then
  !   print*, tenv
  !   stop 'tenv'
  !endif
  isfc=88
!begin  MG  10/28/15 new code to properly filter environmental data
  if(maxval(qvEnv(40:88))<0 .or. maxval(pressEnv(40:88))<0 .or. &
       maxval(tEnv(40:88))<0) then
     tb4emiss=-99
     return   !!skip the rest if qvEnv,..., are negative below 10 km
  endif
!end    MG 10/28/15


  do while(qvEnv(isfc)<0 .or. pressEnv(isfc)>1200)
     isfc=isfc-1
  end do
  if(isfc<88) then
     qvEnv(isfc+1:88)= qvEnv(isfc)
     pressEnv(isfc+1:88)=pressEnv(isfc)
     tEnv(isfc+1:88)=tEnv(isfc)
     sfcTempEnvG=tEnv(isfc)
  end if
  if (sfcTempEnvG<-99) then
     sfcTempEnvG=tEnv(isfc)
     print*, sfcTempEnvG, tEnv(isfc)
     print*, tEnv
     print*, qvEnv
     print*, pressEnv
     print*, 'setEnv',isfc
  end if

  if(minval(pressEnv)<0 .or. minval(qvEnv)<0 .or. minval(tEnv)<0) then
     tb4emiss=-99
     return   !!skip the rest 
  endif

  do i=88,1,-1
     qvenvG(i)=qvEnv(89-i)
     tempG(i)=tEnv(89-i)
     pressEnvG(i)=pressEnv(89-i)
     
     rho=pressEnv(89-i)*100./(tEnv(89-i)*287.)
     
     rhowv=qvEnv(89-i)*rho*1e-3
    
     do j=1,8
        if(tenv(89-i)>100) then
           call GasabsR98(freqs(j),tenv(89-i),Rhowv,pressEnv(89-i)*100.,&
                absair,abswv,ireturn)
           if (ireturn==1) then
              print*, qvEnv
              print*, tEnv
              print*, pressEnv
              print*, qvenvG(i), tEnv(i), pressEnvG(i), sfcTempEnvG, i
              stop
           endif
           wvExt(i,j)=absair+abswv
        else
           if(i<88) then
              wvExt(i,j)=wvExt(i+1,j)
           endif
        endif
     enddo
     !print*, rho, pressEnv(89-i), rhowv, qvEnv(89-i), wvExt(i,4)
  enddo
  umu=cos(53./180.*3.1415)
  do j=1,8
     do i=1,4
        tb4emiss(j,i)=dboux( wvExt(88:1:-1,1),  &
             tenv(88:1:-1), 0.25, i*0.25, umu, &
             sfcTempEnv,88)
     enddo
  enddo
  !stop
end subroutine setEnv

subroutine getemissout2(tbObs,emissout)
  use gEnv
  real :: tbObs(9), emissout(9)
  real :: emissRet(9)
  integer :: i,j
  integer :: ifreq(9)
  ifreq=(/1,1,2,2,3,4,4,5,5/)
!begin  MG 10/28/15 new code to skip code
  if(maxval(tb4emiss)<10) return !!skip the rest if no tb4emiss values are set in setEnv
!end    MG 10/28/15
  do i=1,9
     if(tbObs(i)>0) then
        j=1
        do while(.not.(tb4emiss(ifreq(i),j)<=tbObs(i) .and. &
             tb4emiss(ifreq(i),j+1)>=tbObs(i)) .and. j<3)
           j=j+1
        enddo
        f=(tbObs(i)-tb4emiss(ifreq(i),j))/&
             (tb4emiss(ifreq(i),j+1)+0.1-tb4emiss(ifreq(i),j))
        !if(i.eq.1) &
        !     print*, tbObs(i),tb4emiss(ifreq(i),j),tb4emiss(ifreq(i),j+1)
        emissOut(i)=j*0.25+f*0.25
     endif
  enddo
end subroutine getemissout2
