subroutine iter_profcv2(btop,bzd,bcf,bsfc,zKuL,zKaL,dr,n1d,eps,imu,&
     dn1d,dm1d,rrate1d,zKuC,zKaSim,epst,piaKu,piaKa,itype,dnCoeff_new,dncv,dnp,dzdn,piaSRTKu,relPIASRTKu,&
     dt1,dt2, rrate1d_sub, dn_sub, dm_sub, zkuc_sub, piahb_sub, piaKa_sub,zetaS)
  use tables2
  use tableP2
  use ran_mod
  implicit none
  integer :: btop, bzd, bcf, bsfc, n1d, imu
  real :: zKuL(n1d), dr, zKaL(n1d), dnp(n1d)
  !real :: dmCoeff(2)=(/0.027773772993636318,-0.6624076086959071/)
  real,intent(out) :: dm1d(n1d)
  real,intent(out) :: epst,piaKu,piaKa
  real,intent(out) :: dn1d(n1d), rrate1d(n1d), zKuC(n1d), zKaSim(n1d), dzdn(n1d,n1d)
  real :: dnCoeff_new(2)
  real :: dzKa(n1d), dns(n1d)
  real :: eps, dm, pia, dncv
  integer :: it, ik, k, n1, n1H, itype
  real :: ztrueS,ztrueH,ztrue,f, attKuH, attKaH, attKuS, attKaS
  real :: ftran, probs(31), sumprob,  zeta1d(n1d,31), q
  real,intent(out) :: rrate1d_sub(n1d,31), dn_sub(n1d,31), dm_sub(n1d,31), zkuc_sub(n1d,31), piahb_sub(31)
  real,intent(out) :: piaKa_sub(31),zetaS(31)
  real :: piaKuS, piaKaS, beta, piamax, attKu, attKa, dn, dni
  real :: dnCoeff(2)=(/-0.00570364,  0.13319214/)
  real :: dmCoeff(3)= (/-2.56885571e-04,  7.18909743e-02, -1.60889523e+00/)
  !(/-4.57599678e-04,  7.852247958e-02, -1.78316499e+00/)
  !(/-4.57599678e-04,  8.52247958e-02, -1.78316499e+00/)
  real :: rn, dm_old!, dm_sub(n1d,31)
  real :: zka1,zka2,attka2
  real,intent(out) :: dt1, dt2
  real :: start, finish1, finish2
  integer :: kk, imc, isub
  real :: dmm
  real :: piaSRTKu, dzdn_sub2(n1d,n1d,31), dzdn_sub1(n1d,n1d,31), zkasim_sub(n1d,31)
  integer :: relPIASRTKu
  real :: ddn
  piaKu=0
  !dmCoeff=array([ 0.02893781, -0.69481455])
  dzdn=0.0
  !if(itype.eq.1) then
  dnCoeff=dnCoeff_new
  !end if
  zKuC=zKuL
  dns=0
  dZKa=0
  rrate1d=0
  call cpu_time(start)
  do it=1,2
     piaKu=0.0
     piaKa=0.0
     do k=btop,min(bcf,bzd-1)
        if (zKuC(k+1).gt.10) then
           if (k.lt.bzd) then
              ztrue=zKuC(k+1)
           end if
           if(ztrue.lt.40) then
              ztrueS=ztrue
              ztrueH=0.0
           else
              ztrueS=40
              ztrueH=10*log10(10**(0.1*ztrue)-10**(0.1*39.999))
           endif
           !n1,n2=bisection(zKuSJ,ztrueS-10*dns[k+1])
           call bisection2(zKuSJ,nbinS2,ztrueS-10*dns(k+1), n1)
           attKuS=att13TableS2(n1,imu)*10**dns(k+1)
           attKaS=att35TableS2(n1,imu)*10**dns(k+1)
           !n1H,n2H=bisection(zKuHJ,ztrueH)
           call bisection2(zKuSJ,nbinS2,ztrueH, n1H)
           attKuH=att13TableH(n1H,imu)
           attKaH=att35TableH(n1H,imu)
           piaKu=piaKu+(attKuS+attKuH)*dr
           piaKa=piaKa+(attKaS+attKaH)*dr
           zKuC(k+1)=zKuL(k+1)+piaKu
           zKaSim(k+1)=10*log10(10**(0.1*z35TableS2(n1,imu)+dns(k+1))+&
                10**(0.1*z35TableH(n1H,imu)))-piaKa
           if(zKaL(k+1).gt.10) then
              dZKa(k+1)=zKaSim(k+1)-zKaL(k+1)
           end if
           piaKu=piaKu+(attKuS+attKuH)*dr
           piaKa=piaKa+(attKaS+attKaH)*dr
           dm1d(k+1)=0.5*(d013TableS2(n1,imu)+d013TableH(n1H,imu))
           dn1d(k+1)=dns(k+1)
           rrate1d(k+1)=pr13TableS2(n1,imu)*10**dns(k+1)+pr13TableH(n1H,imu)
        end if
     end do
     if(it.eq.1) then
        do k=btop,min(bcf,bzd-1)
           dns(k+1)=dns(k+1)-0.6*dZKa(k+1)
           if(dns(k+1).gt.0.975) dns(k+1)=0.975
           if(dns(k+1).lt.-1.0) dns(k+1)=-1.0
        end do
     end if
  end do
  
  piaKuS=piaKu+0.0
  piaKaS=piaKa+0.0

  zeta1d=0
  beta=0.71
  
  piamax=54-zKuL(bcf+1)
  if(piamax<0) piamax=56-zKuL(bcf+1)
  if(piamax<0) piamax=4
  q=0.2*log(10.0)
  zKuC(bzd+1:bcf+1)=zKuL(bzd+1:bcf+1)+piaKuS
  epst=eps
  attKu=0.0
  attKa=0.0
  f=1.0
  zeta1d=0
  sumprob=0
  piaKu=0.0
  rrate1d(bzd+1:bcf+1)=0.0
  dm1d(bzd+1:bcf+1)=0.0
  dn1d(bzd+1:bcf+1)=0.0
  rrate1d_sub=0
  dm_sub=0
  dn_sub=0
  zkasim_sub=0
  dzdn_sub2=0
  dzdn_sub1=0
  zkasim(bzd+1:bcf+1)=0
  dzdn(bzd+1:bcf+1,bzd+1:bcf+1)=0
  call cpu_time(finish1)
  ddn=0.025
!!$OMP PARALLEL DO  PRIVATE(k,dn,n1,attKu,attKa,zKa1,zKa2)
  do isub=1,31
     zetaS(isub)=0
     do k=bzd,bcf
        if (zKuL(k+1).gt.10) then
           call bisection2(zKudN,nbins,zKuL(k+1)+piaKuS-10*(isub-16)*ddn-10*dncv, n1)
           dn_sub(k+1,isub)=(zKudN(n1)-zKuSJ(n1))/10.0+(isub-16)*ddn+dncv
           dn=dn_sub(k+1,isub)
           attKu=att13Table(n1,imu)*10**dn
           attKa=att35Table(n1,imu)*10**dn
           zetaS(isub)=zetaS(isub)+attKu*dr
           zeta1d(k+1,isub)=zetaS(isub)
        else
           attKu=0.0
           attKa=0.0
        endif
     end do
     if (q*beta*zetaS(isub).lt.0.985) then
        piaKa_sub(isub)=piaKas
        piaHB_sub(isub)=-10/beta*log10(1-q*beta*zetaS(isub))+attKu*(bsfc-bcf)*2*dr
        do k=bzd,bcf
           if(zKuL(k+1).gt.10) then
              zKuC_sub(k+1,isub)=zKuL(k+1)-10/beta*log10(1-q*beta*zeta1d(k+1,isub))
              n1=int((zKuC_sub(k+1,isub)-zmin-10*dn_sub(k+1,isub))/dzbin)
              dn=dn_sub(k+1,isub)
              attKa=att35Table(n1,imu)*10**dn
              piaKa_sub(isub)=piaKa_sub(isub)+attKa*dr
              if(n1.lt.1) n1=1
              if(n1.gt.nbins) n1=nbins
              dm_sub(k+1,isub)=d013Table(n1,imu)
              rrate1d_sub(k+1,isub)=pr13Table(n1,imu)*10**dn_sub(k+1,isub)
              zKa1=z35Table(n1,imu)+10*dn
              zkasim_sub(k+1,isub)= zKa1-piaKa_sub(isub)
              piaKa_sub(isub)=piaKa_sub(isub)+attKa*dr
              if(n1.lt.nbins-4) then
                 zKa2=z35Table(n1+4,imu)+10*(dn-0.1)
                 attKa2=att35Table(n1+4,imu)*10**(dn-0.1)
                 dzdn_sub2(k+1,k+1,isub)=dzdn_sub2(k+1,k+1,isub)+(zka2-attKa2*dr)-((zka1-attKa*dr))
                 dzdn_sub1(k+1,k+1,isub)=dzdn_sub1(k+1,k+1,isub)+zka1-attKa*dr
                 do kk=k+1,bcf
                    dzdn_sub2(kk+1,k+1,isub)=dzdn_sub2(kk+1,k+1,isub)-(attKa2*dr-attKa*dr)
                    dzdn_sub1(kk+1,k+1,isub)=dzdn_sub1(kk+1,k+1,isub)-attKa*dr
                 end do
              end if
           else
              attKa=0
           end if
        end do
        piaKa_sub(isub)=piaKa_sub(isub)+attKa*(bsfc-bcf)*2*dr
        probs(isub)=exp(-((isub-16)*ddn)**2/0.5**2)
        if(rrate1d_sub(bcf+1,isub).gt.150) probs(isub)=probs(isub)*0.9
     else
        probs(isub)=0.0
     end if
  end do
  !print*,zetaS,bzd,bcf
  !!$OMP END PARALLEL DO
  piaKu=0.0
  piaKa=0.0
  sumprob=0.0
  do isub=1,31
     do k=bzd,bcf
        if(zkuL(k+1).gt.10) then
           rrate1d(k+1)= rrate1d(k+1)+probs(isub)*rrate1d_sub(k+1,isub)
           dm1d(k+1)=dm1d(k+1)+probs(isub)*dm_sub(k+1,isub)
           dn1d(k+1)=dn1d(k+1)+probs(isub)*dn_sub(k+1,isub)
           zkasim(k+1)=zkasim(k+1)+10**(0.1*zkasim_sub(k+1,isub))*probs(isub)
           dzdn(k+1,k+1)=dzdn(k+1,k+1)+10**(0.1*zkasim_sub(k+1,isub)+0.1*dzdn_sub2(k+1,k+1,isub))*probs(isub)
           do kk=k+1,bcf
              if(zkuL(kk+1).gt.10) then
                 dzdn(kk+1,k+1)=dzdn(kk+1,k+1)+10**(0.1*zkasim_sub(kk+1,isub)+0.1*dzdn_sub2(kk+1,k+1,isub))*probs(isub)
              end if
           end do
        end if
     end do
     piaKu=piaKu+probs(isub)*piaHB_sub(isub)
     piaKa=piaKa+probs(isub)*piaKa_sub(isub)
     sumprob=sumprob+probs(isub)
  end do
 
  piaKu=piaKu/(sumprob+1e-10)+piaKuS
  piaKa=piaKa/(sumprob+1e-10)
  do k=bzd,bcf
     if(zkuL(k+1).gt.10) then
        rrate1d(k+1)= rrate1d(k+1)/(sumprob+1e-10)
        dm1d(k+1)=dm1d(k+1)/(sumprob+1e-10)
        dn1d(k+1)=dn1d(k+1)/(sumprob+1e-10)
        zkasim(k+1)=log10(zkasim(k+1)/(sumprob+1e-10))*10.0
        dzdn(k+1,k+1)=(log10(dzdn(k+1,k+1)/(sumprob+1e-10))*10.0-zkasim(k+1))/(-0.1)
        do kk=k+1,bcf
           if(zkuL(kk+1).gt.10) then
              dzdn(kk+1,k+1)=(log10(dzdn(kk+1,k+1)/(sumprob+1e-10))*10.0-log10(zkasim(kk+1)/(sumprob+1e-9))*10)/(-0.1)
           end if
        end do
     end if
  end do
  call cpu_time(finish2)
  dt1=finish1-start
  dt2=finish2-start
  if(rrate1d(bcf+1).gt.150) then
     print*, rrate1d_sub(bcf+1,1:31:3)
     print*, 'probs=',probs(1:31:3)
     print*, 'zetaS=',zetaS(1:31:3)
  end if
end subroutine iter_profcv2
