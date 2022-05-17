subroutine rainprofst(n1,zku_obs,zka_obs,dpiaSRT,piakus,piakas,&
     reldpia,nc,dr,wzku,wzka,wpia,rrate_in,dn_in,rrate_out,dn_out)
  use tablep2
  use tables2
  implicit none
  integer :: n1, reldpia, nc
  real :: w1, w2, wzku, wzka, wpia
  real :: zku_obs(n1),zka_obs(n1),dpiaSRT, dr
  real :: piakus,piakas
  real :: piaku(n1), piaka(n1)
  real :: lam1(n1), lam2(n1)
  real :: df_dlam1(n1), df_dlam2(n1)
  real :: df_dlogr(n1), df_dn(n1)
  real :: rrate_in(n1),dn_in(n1), logrrate(n1)
  real,intent(out) :: rrate_out(n1),dn_out(n1)
  integer :: i, n1j, it, i1
  real :: piaKut,piaKat, dzku, dzka, dzkudr, dzkadr, &
       dattkudr(n1), dattkadr(n1), &
       dattkudn(n1), dattkadn(n1), dpia, attkuL, attKaL, wdn, alpha
  real :: fobj
  df_dlam1=0
  df_dlam2=0
 
  wdn=0.1
  alpha=0.05
  dn_out=dn_in
  do i=1,n1
     if(rrate_in(i).gt.0) then
        logrrate(i)=log10(rrate_in(i))
     else
        logrrate(i)=-999
     end if
  end do
  print*,zku_obs
  print*, logrrate
  do it=1,10
     fobj=0.0
     df_dlogr=0
     df_dn=0
     piaKut=piaKuS
     piaKat=piaKaS
     attKuL=0
     attKaL=0
     do i=1,n1
        if (zku_obs(i).gt.10) then
           call bisection2(logRJ,nbins,logrrate(i)-dn_out(i),n1j)
           piaKut=piaKut+attKuJ(n1j)*10**dn_out(i)*dr
           piaKat=piaKat+attKaJ(n1j)*10**dn_out(i)*dr
           attKuL=attKuJ(n1j)*10**dn_out(i)
           attKaL=attKaJ(n1j)*10**dn_out(i)
           if (n1j.gt.1.and.n1j.lt.nbins) then
              dzkudr=(zkusj(n1j)-zkusj(n1j-1))/(logRJ(n1j)-logRJ(n1j-1))
              dattkudr(i)=(attKuJ(n1j)-attKuJ(n1j-1))*10**dn_out(i)/&
                   (logRJ(n1j)-logRJ(n1j-1))
              dattkudn(i)=attKuJ(n1j)*10**dn_out(i)*log(10.0)
              dzku=zku_obs(i)-zkusj(n1j)-10*dn_out(i)+piaKut
              !print*,'zku','dzku',zkusj(n1j)+10*dn_out(i)-piaKut, dzku
              fobj=fobj+dzku**2
              df_dlogr(i)=df_dlogr(i)+wzku*dzku*(-dzkudr+dattkudr(i)*dr)
              df_dn(i)=df_dn(i)+wzku*dzku*(-10+dattkudn(i)*dr)
              do i1=1,i-1
                 df_dlogr(i1)=df_dlogr(i1)+wzku*dzku*(2*dattkudr(i1)*dr)
                 df_dn(i1)=df_dn(i1)+wzku*dzku*(2*dattkudn(i1)*dr)
              enddo
              if (zka_obs(i).gt.10) then
                 dzkadr=(zkasj(n1j)-zkasj(n1j-1))/(logRJ(n1j)-logRJ(n1j-1))
                 dzka=zka_obs(i)-zkasj(n1j)-10*dn_out(i)+piaKat
                 !print*,'zka','dzka',zkasj(n1j)+10*dn_out(i)-piaKat, dzka
                 dattkadr(i)=(attKaJ(n1j)-attKaJ(n1j-1))*10**dn_out(i)&
                      /(logRJ(n1j)-logRJ(n1j-1))
                 dattkadn(i)=attKaJ(n1j)*10**dn_out(i)*log(10.0)
                 df_dlogr(i)=df_dlogr(i)+wzka*dzka*(-dzkadr+dattkadr(i)*dr)
                 df_dn(i)=df_dn(i)+wzka*dzka*(-10+dattkadn(i)*dr)
                 fobj=fobj+dzka**2
                 do i1=1,i-1
                    df_dlogr(i1)=df_dlogr(i1)+wzka*dzka*(2*dattkadr(i1)*dr)
                    df_dn(i1)=df_dn(i1)+wzka*dzka*(2*dattkadn(i1)*dr)
              enddo
              endif
           end if
           piaKut=piaKut+attKuJ(n1j)*10**dn_out(i)*dr
           piaKat=piaKat+attKaJ(n1j)*10**dn_out(i)*dr
        end if
     end do
     piaKut=piaKut+attKuL*2*nc*dr
     piaKat=piaKat+attKaL*2*nc*dr
     print*, it, 'Fobj=',fobj
     if(reldpia.eq.1) then
        dpia=dpiaSRT-(piaKat-piaKut)
        do i=1,n1
           if (zku_obs(i).gt.10.and.zka_obs(i).gt.10) then
              df_dlogr(i)=df_dlogr(i)+wpia*dpia*(-dattkadr(i)+dattkudr(i))*dr*2
              df_dn(i)=df_dn(i)+wpia*dpia*(-dattkadn(i)+dattkudn(i))*dr*2
           end if
        end do
     end if
     !print*, df_dn
     !print*, df_dlogr
     do i=1,n1
        if (zku_obs(i).gt.10.) then
           !df_dn(i)=df_dn(i)+wdn*dn_out(i)
           dn_out(i)=dn_out(i)-alpha*df_dn(i)
           logrrate(i)=logrrate(i)-alpha*df_dlogr(i)
        end if
     end do
     print*, piaKat-piaKut,dpiaSRT
  end do
  rrate_out=10**logrrate
end subroutine rainprofst


subroutine rainprofstg(n1,zku_obs,zka_obs,dpiaSRT,piakus,piakas,&
     reldpia,nc,dr,wzku,wzka,wpia,rrate_in,dn_in,nens,rrate_out,dn_out,zkusim,zkasim,&
     zku_out,zka_out,rrens,yEns,xEns,dy,pia_out,dm_out,attOut)
  use tablep2
  use tables2
  use ran_mod
  implicit none
  integer :: n1, reldpia, nc
  real :: w1, w2, wzku, wzka, wpia
  real :: zku_obs(n1),zka_obs(n1),dpiaSRT, dr
  real :: piakus,piakas
  real :: piaku(n1), piaka(n1)
  real :: lam1(n1), lam2(n1)
  real :: df_dlam1(n1), df_dlam2(n1)
  real :: df_dlogr(n1), df_dn(n1)
  real :: rrate_in(n1),dn_in(n1), logrrate(n1)
  real,intent(out) :: rrate_out(n1),dn_out(n1),zkusim(nens,n1),zkasim(nens,n1),&
       rrens(nens,n1),pia_out(2),yEns(nens,2*n1+1),xEns(nens,2*n1+1), dy(2*n1+1)
  real, intent(out) :: attOut(n1,2)
  real,intent(out) :: zku_out(n1), zka_out(n1), dm_out(n1)
  integer :: i, n1j, it, i1, nens
  real :: piaKut,piaKat, dzku, dzka, dzkudr, dzkadr, &
       dattkudr(n1), dattkadr(n1), &
       dattkudn(n1), dattkadn(n1), dpia, attkuL, attKaL, wdn, alpha
  real :: fobj
  real :: rn,dnp,rrn
  real :: dx(2*n1+1), sigma
  df_dlam1=0
  df_dlam2=0
 
  wdn=0.1
  alpha=0.01
  dn_out=dn_in
  do i=1,n1
     if(rrate_in(i).gt.0) then
        logrrate(i)=log10(rrate_in(i))
     else
        logrrate(i)=-999
     end if
  end do
  do it=1,nens
     piaKut=piaKuS
     piaKat=piaKaS
     attKuL=0
     attKaL=0
     rn=normal2(0.0,1.0)
     do i=1,n1
        if (zku_obs(i).gt.10) then
           rn=0.5*rn+0.5*normal2(0.0,1.0)
           dnp=dn_out(i)-0.1+rn*0.5
           rrn=normal2(0.0,0.2)
           call bisection2(logRJ,nbins,rrn+logrrate(i)-dnp,n1j)
           rrens(it,i)=10**(logrrate(i)+rrn)
           attKuL=attKuJ(n1j)*10**dnp
           attKaL=attKaJ(n1j)*10**dnp
           piaKut=piaKut+attKuJ(n1j)*10**dnp*dr
           piaKat=piaKat+attKaJ(n1j)*10**dnp*dr
           zkusim(it,i)=zkusj(n1j)+10*dnp-piaKut
           zkasim(it,i)=zkasj(n1j)+10*dnp-piaKat
           yEns(it,i)=zkusim(it,i)
           yEns(it,n1+i)=zkasim(it,i)
           xEns(it,i)=rrens(it,i)
           xEns(it,i+n1)=dnp
           piaKut=piaKut+attKuJ(n1j)*10**dnp*dr
           piaKat=piaKat+attKaJ(n1j)*10**dnp*dr
        end if
     end do
     piaKut=piaKut+attKuL*2*nc*dr
     piaKat=piaKat+attKaL*2*nc*dr
     pia_out(1)=piaKut
     pia_out(2)=piaKat
     yEns(it,2*n1+1)=piaKat-piaKut
     xEns(it,2*n1+1)=piaKat-piaKut
  enddo
  do i=1,n1
     if (zku_obs(i).gt.10) then
        dy(i)=zku_obs(i)-sum(zkusim(:,i))/nens
        if (zka_obs(i).gt.10) then
           dy(i+n1)=zka_obs(i)-sum(zkasim(:,i))/nens
        end if
     end if
  end do
  if(reldpia.eq.1) then
     dy(2*n1+1)=dpiaSRT-sum(yEns(:,2*n1+1))/nens
  end if
  sigma=4
  call kgain(xENS, yEns, dy, 2*n1+1, 2*n1+1, nEns, sigma, dx)
  !print*, dx
  !stop
  !La Solution
  piaKut=piaKuS
  piaKat=piaKaS
  do i=1,n1
     if (zku_obs(i).gt.10) then
        rrate_out(i)=sum(xEns(:,i))/nEns+0.95*dx(i)
        if(rrate_out(i).lt.minval(xEns(:,i))) then
           rrate_out(i)=minval(xEns(:,i))
        end if
        dnp=sum(xEns(:,i+n1))/nEns+0.95*dx(i+n1)
        call bisection2(logRJ,nbins,log10(rrate_out(i))-dnp,n1j)
        dn_out(i)=dnp
        attKuL=attKuJ(n1j)*10**dnp
        attKaL=attKaJ(n1j)*10**dnp
        piaKut=piaKut+attKuJ(n1j)*10**dnp*dr
        piaKat=piaKat+attKaJ(n1j)*10**dnp*dr
        zku_out(i)=zkusj(n1j)+10*dnp-piaKut
        zka_out(i)=zkasj(n1j)+10*dnp-piaKat
        dm_out(i)=dmJ(n1j)
        attOut(i,1)=attKuJ(n1j)*10**dnp*dr
        attOut(i,2)=attKaJ(n1j)*10**dnp*dr
        piaKut=piaKut+attKuJ(n1j)*10**dnp*dr
        piaKat=piaKat+attKaJ(n1j)*10**dnp*dr
     end if
  end do
  piaKut=piaKut+attKuL*2*nc*dr
  piaKat=piaKat+attKaL*2*nc*dr
  pia_out(1)=piaKut
  pia_out(2)=piaKat
end subroutine rainprofstg
