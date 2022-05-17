

subroutine calctb(sfc_wind,umu,kext,salb,asym,node,ic,jc,&
     ngates,nmfreqm,hh,tb,emtb,emis_out,tpwGMI,hfreez,randemiss,imemb,nmemb,&
     iprof,jprof,w1rand,jrand)
  use cldclass 
  use microwFreq
  use ran_mod
  use geophysEns
  use emissMod
  use gENV
  use LUT_def
  implicit none
  
  integer:: ngates, nmfreqm, npart
  real :: umu, sfc_temp,  dr, w1rand(nmemb), jrand(nmemb*50)
  integer :: i, ifreq, npol,  node(5)
  real :: hh(ngates)
  real :: kext(ngates,nmfreqm), salb(ngates,nmfreqm),&
       asym(ngates,nmfreqm)
  real :: atm_ext, tavg, kextcw
  real :: tb3, fisot, emis, ebar 
  real :: tb(2*nmfreqm),  tb2(2,nmfreqm),  emtb(2*nmfreqm), emis_out(2*nmfreqm)
  logical :: prnt(3), lambert
  real ::  kexttot(40), salbtot(40), asymtot(40)
  real ::  kexttotp(40), salbtotp(40), asymtotp(40), tbp
  real ::  kextrn(40), salbrn(40), asymrn(40), temp2(0:40)
  integer :: countrn(40)
  integer ::  i0, j, ic, jc, ic0, icmin(1), iprof, jprof, stype
  real :: sfc_wind
  real:: kext_clw, tpwGMI, hfreez
  integer :: itop, ibot, nlayerA, imemb, nmemb
  real    :: randemiss(nmemb,nmfreqm,2)
  real :: rhPCij(nRhEofs), cldwPCij(nCldwEofs), cldw(nlayer), rh(nlayer)
  real ::  wv_extMemb(nmfreqm,nlayer), cld_extMemb(nmfreqm,nlayer), w1
  logical :: isnan
  integer :: ifirst
  real  ::                  emissv(n_chan)
  real  ::                  emissh(n_chan)
  real  ::                  emissv_std(n_chan)
  real  ::                  emissh_std(n_chan)
  real  ::                  wfract, sfc_wind0, dboux, emisold
  lambert = .false.
  
  call getwfraction(elat,elon,wfract)

  sfc_temp=temp(0)-(4.5-hfreez)*6
!MG Feb 2015
  temp2=temp
  sfc_temp=sfcTempEnvG
  do i=1,40
     temp2(i)=0.5*(tempG(2*i)+tempG(2*i+1))
  enddo
  temp2(0)=tempG(1)-0.125*6.

!MG Feb 2015
  w1=0.1
  sfc_wind0=sfc_wind
  if(sfc_wind0<2) sfc_wind0=2.
  !SJM 10/19/15 Emissivities for all surface types are now set in the RadarRet subroutine.
  emissv(1) = emis_out(1)
  emissv(2) = emis_out(3)
  emissv(3) = emis_out(5)
  emissv(4) = emis_out(6)
  emissv(5) = emis_out(8)
  emissv(6) = emis_out(10)
  emissh(1) = emis_out(2)
  emissh(2) = emis_out(4)
  emissh(3) = emis_out(4)
  emissh(4) = emis_out(7)
  emissh(5) = emis_out(9)
  emissh(6) = emis_out(11)
  !end SJM 10/19/15
  
  call interpolPC(imemb+1, rhPCij, cldwPCij, cldw, rh)

!begin  MG 9/20/13 add statement to force passing of variable values to
!       cldwcoeff
  cldwcoeff(1:10, imemb + 1) = cldwPCij(1:10)
!end    MG 9/20/13

  if(imemb==0) then
     do i= node(1), -node(5)
        !print*, kext(i,1:5)
     enddo
     !print*, 'In RTE'
     !print*, iRad,jRad,rhPCij(1:2), rh
  endif
  w1=0.25*w1rand(imemb+1)
  !w1=0.125
  do i=1,nlayer
    do ifreq=1,nmfreqm
       j=1.0*rh(i)/2.
       tavg=(temp2(i)+temp2(i-1))*0.5
       !if(tavg<267 .and. imemb/2*2==imemb) j=j/2
       if(j<=0) j=1
       if(j>50) j=50
       !j=1
       wv_extMemb(ifreq,i)=atm_extMWRg(ifreq,i,j)
       wv_extMemb(ifreq,i)=w1*atm_extMWRg(ifreq,i,j)+&
            (1-w1)*0.5*(wvExt(2*i,ifreq)+wvExt(2*i-1,ifreq))
       
       if(imemb/2*2==imemb) then
          if(i>0.5*(88-node(3))) then
             j=j*jrand(imemb*50+i)+1
!begin WSO 5/12/15
             if(j<=0) j=1
             if(j>50) j=50
!end   WSO 5/12/15
             wv_extMemb(ifreq,i)=atm_extMWRg(ifreq,i,j)
          endif
       endif
       if(cldw(i)<0) cldw(i)=0.
       call gcloud(mfreq(ifreq),tavg,0.5*cldw(i),cld_extMemb(ifreq,i))
    enddo
 enddo

  if(node(5)>100) node(5)=100


  fisot   = 2.7
  tb=0
  !if(imemb==1 .and. iprof==38 .and. jprof==190) then

  do i=node(1),node(5)
     salb(i,8)=salb(i,7)
     kext(i,8)=kext(i,7)
     asym(i,8)=asym(i,7)
  enddo

  !endif
121 format(12(F7.2,1x))
  do ifreq=1,nmfreqm
     kextrn=0
     salbrn=0
     asymrn=0
     countrn=0
     do i=1,node(5)
        i0=((hh(i)+4.5-hfreez)/drrte)+1
        if(imemb==-1 .and. iprof==14 .and. jprof==232) then
           !print*, i, ifreq, salb(i,ifreq),asym(i,ifreq),kext(i,ifreq)
        endif
        if(i0>0 .and. i0<=nlayer .and. kext(i,ifreq)>-9) then
           asymrn(i0)=asymrn(i0)+salb(i,ifreq)*asym(i,ifreq)*kext(i,ifreq)
           kextrn(i0)=kextrn(i0)+kext(i,ifreq)
           salbrn(i0)=salbrn(i0)+salb(i,ifreq)*kext(i,ifreq)
           if(isnan(asymrn(i0))) then
              print*, asymrn(i0), kextrn(i0), salbrn(i0)
              print*,'rte'
              stop
           endif
           countrn(i0)=countrn(i0)+1
           if(imemb==111) &
                print*, 'skt=', salb(i,ifreq), &
                kext(i,ifreq), salbrn(i0)/kextrn(i0), i, i0, &
                countrn(i0)
        endif
     enddo
     do i=1,nlayer
        if(countrn(i)>0) then
           asymrn(i)=asymrn(i)/countrn(i)
           kextrn(i)=kextrn(i)/countrn(i)
           salbrn(i)=salbrn(i)/countrn(i)
           if(imemb==0) then
              if(salbrn(i) .gt. kextrn(i)) then
                 print*, 'sk=', salbrn(i), kextrn(i), i
                 print*, 'rte'
                 stop
              endif
           endif
           if(isnan(salbrn(i))) then
           endif
           if(isnan(asymrn(i))) then
              print*, asymrn(i)
           endif
        endif
     enddo
     ifirst=1
     if(countrn(1)==0 ) then
        do while(countrn(ifirst)==0 .and. ifirst<10)
           ifirst=ifirst+1
        enddo
        asymrn(1:ifirst-1)=asymrn(ifirst)
        kextrn(1:ifirst-1)=kextrn(ifirst)
        salbrn(1:ifirst-1)=salbrn(ifirst)
     endif

     do i=1,nlayer
        tavg=(temp2(i)+temp2(i-1))*0.5
        call gcloud(mfreq(ifreq),tavg,0.5*cc(i,jc),kext_clw)
        atm_ext=wv_extMemb(ifreq,i)
        !kext_clw=cld_extMemb(ifreq,i)
        kexttot(i)=atm_ext+kextrn(i)+kext_clw
     
        if(salbrn(i)> 1e-8) then
           asymtot(i)=(asymrn(i))/(salbrn(i))
        else
           asymtot(i)=0.
        endif
        salbtot(i)=(salbrn(i))/kexttot(i)
        if(salbtot(i)>1) then
           if(imemb==0) then
              print*, kexttot(i),  kextrn(i), salbrn(i), atm_ext, kext_clw
              print*, 'rte'
              stop
           endif
           salbtot(i)=1.
        endif
        if(asymtot(i)>1) then
           
           asymtot(i)=1.
        endif
        if(asymtot(i)<-1) asymtot(i)=-1.
     enddo
     ibot=0
     
     do while(height(ibot)<4.5-hfreez .and. ibot<10)
        ibot=ibot+1
     enddo

     if(ibot<1) ibot=1
     nlayerA=nlayer+1-ibot
     do npol=0,1

        if(npol==1) then
           if(ifreq<=6) then
              emis=emissv(ifreq)!+ randemiss(imemb+1,ifreq,npol+1)*0.5*0.15/2.
           else
              emis=emissv(6)!+ randemiss(imemb+1,ifreq,npol+1)*0.5*0.15/2.
           endif
        else
           if(ifreq<=6) then
              emis=emissh(ifreq)!+ randemiss(imemb+1,ifreq,npol+1)*0.5*0.15/2.
           else
              emis=emissh(6)!+ randemiss(imemb+1,ifreq,npol+1)*0.5*0.15/2.
           endif
        endif
        if(emis<0.1) emis=0.1
        !if(emis>0.999) emis=0.999
        !if(sfc_wind>=0 .and. wfract .gt. 50) then
        !   call emit(mfreq(ifreq), npol, sfc_temp, sfc_wind, &
        !        umu, emis, ebar)
        !   emisold=emis
        !   if(ifreq<=5) then 
        !      call intplte_emis(ifreq,1-NPOL,sfc_temp,sfc_wind,0.,52.8,emis,ebar)
        !      if(imemb==-1) then
        !         print*, emisold, emis, npol, sfc_wind, mfreq(ifreq)
        !      endif
        !   endif
        !endif
        if(emissv(1)<0) then
           tb=-99
           return
        endif
        !if(emis>0.999) emis=0.999
        if(emis<0.1) emis=0.1
        ebar=emis
        emis_out(npol*nmfreqm+ifreq) = emis
        !endif
        if(imemb==-1 .and. iprof==38 .and. jprof==190) then
           if(ifreq>5 .and. npol==0) then
              print*, ifreq
              print*, hfreez
              do i=ibot,nlayer
              write(*,201), temp(i), height(i),  kexttot(i), salbtot(i), &
                   asymtot(i),wv_extMemb(ifreq,i), qvEnvG(i*2)
           enddo
           endif
        endif
        call radtran(umu,nlayerA, tb(npol*nmfreqm+ifreq), sfc_temp, &
             temp(ibot-1:nlayer), height(ibot-1:nlayer)-height(ibot-1), &
             kexttot(ibot:nlayer), salbtot(ibot:nlayer), &
             asymtot(ibot:nlayer), fisot, emis, &
             ebar, lambert, prnt)
        if(ifreq>5 .and. imemb==-1 .and. iprof==38 .and. jprof==190) &
             print*, tb(npol*nmfreqm+ifreq), ifreq, npol
        emtb(npol*nmfreqm+ifreq)=&
             dboux( kexttot(ibot:nlayer),  temp(ibot-1:nlayer), &
             0.5, emis,  umu, sfc_temp, nlayerA) 
        if(isnan(tb(npol*nmfreqm+ifreq)).and.imemb==0) then   ! sanity check
           print*, 'Exception in rterain.f90'
           print*, mfreq(ifreq)
           print*, sfc_wind
           print*, emis, ifreq, emissv(ifreq)
           print*, sfc_temp, hfreez
           do i=ibot,nlayer
              print*, temp(i), height(i),  kexttot(i), salbtot(i), &
                   asymtot(i), fisot
           enddo
           print*, umu
           print*, 'rte'
           stop
        endif
!        if((tb(npol*nmfreqm+ifreq))<240 .and. mfreq(ifreq)>800) then  
!           write(31,*) emis, sfc_temp, nlayer+1-ibot
!           do i=ibot,nlayer
!              write(31,100)  i, kexttot(i), salbtot(i), asymtot(i), &
!                   height(i-1)-height(ibot-1),temp(i-1)
!           enddo
!           write(31,*) tb(npol*nmfreqm+ifreq),mfreq(ifreq)
!           do i=ibot,nlayer
!              kexttotp=kexttot
!              kexttotp(i)=1.1*kexttot(i)
!              call radtran(umu,nlayerA, tbp, sfc_temp, &
!                   temp(ibot-1:nlayer), height(ibot-1:nlayer)-height(ibot-1), &
!                   kexttotp(ibot:nlayer), salbtot(ibot:nlayer), &
!                   asymtot(ibot:nlayer), fisot, emis, &
!                   ebar, lambert, prnt)
!              write(31,101) i-ibot, tbp- tb(npol*nmfreqm+ifreq)
!           enddo
!        endif
     enddo
  enddo
201 format(7(F10.4,1x))
if(imemb==0) then
   !print*, (tb(ifreq),ifreq=5,nmfreqm)
endif
!print*, tb
!stop
101 format(i3,1x,2(F7.3,1x))
100 format(i3,8(F9.4,1x))


end subroutine calctb


subroutine avgScProp(sfc_wind,umu,kext,salb,asym,node,ic,jc,&
     ngates,nmfreqm,hh,tb,tpwGMI,hfreez,randemiss,imemb,nmemb)
  use cldclass 
  use microwFreq
  use ran_mod
  use geophysEns
  use emissMod

  implicit none
  
  integer:: ngates, nmfreqm, npart
  real :: umu, sfc_temp,  dr
  integer :: i, ifreq, npol,  node(5)
  real :: hh(ngates)
  real :: kext(ngates,nmfreqm), salb(ngates,nmfreqm),&
       asym(ngates,nmfreqm)
  real :: atm_ext, tavg, kextcw
  real :: tb3, fisot, emis, ebar 
  real :: tb(2*nmfreqm),  tb2(2,nmfreqm)
  logical :: prnt(3), lambert
  real ::  kexttot(40), salbtot(40), asymtot(40)
  real ::  kexttotp(40), salbtotp(40), asymtotp(40), tbp
  real ::  kextrn(40), salbrn(40), asymrn(40)
  integer :: countrn(40)
  integer ::  i0, j, ic, jc, ic0, icmin(1)
  real :: sfc_wind
  real:: kext_clw, tpwGMI, hfreez
  integer :: itop, ibot, nlayerA, imemb, nmemb
  real    :: randemiss(nmemb,nmfreqm,2)
  real :: rhPCij(nRhEofs), cldwPCij(nCldwEofs), cldw(nlayer), rh(nlayer)
  real ::  wv_extMemb(nmfreqm,nlayer), cld_extMemb(nmfreqm,nlayer)
  logical :: isnan
  integer :: ifirst
  real  ::                  emissv(n_chan)
  real  ::                  emissh(n_chan)
  real  ::                  emissv_std(n_chan)
  real  ::                  emissh_std(n_chan)
  lambert = .false.

 

  call interpolPC(imemb+1, rhPCij, cldwPCij, cldw, rh)

  do i=1,nlayer
     do ifreq=1,nmfreqm
       j=rh(i)/2.
       if(j<=0) j=1
       if(j>50) j=1
       wv_extMemb(ifreq,i)=atm_extMWRg(ifreq,i,j)
       tavg=(temp(i)+temp(i-1))*0.5
       if(cldw(i)<0) cldw(i)=0.
       call gcloud(mfreq(ifreq),tavg,cldw(i),cld_extMemb(ifreq,i))
    enddo
 enddo

  if(node(5)>100) node(5)=100

 
  fisot   = 2.7

  tb=0
  
  do ifreq=1,nmfreqm
     kextrn=0
     salbrn=0
     asymrn=0
     countrn=0
     do i=node(1),node(5)
        i0=((hh(i)+4.5-hfreez)/drrte)+1
        if(i0>0 .and. i0<=nlayer .and. kext(i,ifreq)>-9) then
           asymrn(i0)=asymrn(i0)+salb(i,ifreq)*asym(i,ifreq)*kext(i,ifreq)
           kextrn(i0)=kextrn(i0)+kext(i,ifreq)
           salbrn(i0)=salbrn(i0)+salb(i,ifreq)*kext(i,ifreq)
           countrn(i0)=countrn(i0)+1
        endif
     enddo
     do i=1,nlayer
        if(countrn(i)>0) then
           asymrn(i)=asymrn(i)/countrn(i)
           kextrn(i)=kextrn(i)/countrn(i)
           salbrn(i)=salbrn(i)/countrn(i)
           if(isnan(salbrn(i))) then
           endif
        endif
     enddo
     ifirst=1
     if(countrn(1)==0 ) then
        do while(countrn(ifirst)==0 .and. ifirst<10)
           ifirst=ifirst+1
        enddo
        asymrn(1:ifirst-1)=asymrn(ifirst)
        kextrn(1:ifirst-1)=kextrn(ifirst)
        salbrn(1:ifirst-1)=salbrn(ifirst)
     endif
     do i=1,nlayer
        tavg=(temp(i)+temp(i-1))*0.5
        atm_ext=wv_extMemb(ifreq,i)
        kext_clw=cld_extMemb(ifreq,i)
        kexttot(i)=atm_ext+kextrn(i)+kext_clw
        salbtot(i)=(salbrn(i))/kexttot(i)
        if(salbtot(i)*kexttot(i) > 1e-8) then
           asymtot(i)=(asymrn(i))/(salbtot(i))
        else
           asymtot(i)=0.
        endif
       
        if(salbtot(i)>1) salbtot(i)=1.
        if(asymtot(i)>1) asymtot(i)=1.
        if(asymtot(i)<-1) asymtot(i)=-1.
     enddo
     do i=1,nlayer
        kexttotEns(ifreq,i,imemb+1)=kexttot(i)
        salbtotEns(ifreq,i,imemb+1)=salbtot(i)
        asymtotEns(ifreq,i,imemb+1)=asymtot(i)
     enddo
     ibot=0
 
   
  enddo



end subroutine avgScProp
