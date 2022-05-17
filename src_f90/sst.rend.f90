!  SFM  04/06/2013  Code changes from M.Grecu
!  SFM  08/09/2013  Modified deafault values

subroutine sst2(gmiData,geoData,mm,n1b11)
  use f90DataTypes
  implicit none
  type (cgMIDataType)    :: gMIData
  type (geoDataType) :: geoData
  integer :: k, i, row, col, mm, n1b11
  real :: tb(10)
  integer :: nmfreqm, iLandSea(5,5)
  real :: sstclust(1200,10)
  integer :: left(n1b11), right(n1b11), irflag(n1b11), ilsea(n1b11)
  integer :: rc, lc, j, i1, j1
  integer :: isol, ic, irainflag
  real :: tb10v, tb10h, tb19v, tb19h, tb21v, tb37v, tb37h, tb85h, tb85v, hfreez
  real :: sstv, rms, rmsmin
  real :: w1, w2, wfract
  integer :: igetlandsea


  open(10,file='AncData/sst.t.clust')
  do i=1,1200
     read(10,*) (sstclust(i,j),j=1,10)
  enddo

  gMIData%sfc_wind=-999
  gMIData%sst= -999

  nmfreqm=5
  do k=70,150
     left=-1
     right=-1
     ilsea=1
     lc=-1
     do i=1,gmiData%n1b11
        tb10v=gmiData%gmiS1(1,k,i)
        tb10h=gmiData%gmiS1(2,k,i)
        tb19v=gmiData%gmiS1(3,k,i)
        tb19h=gmiData%gmiS1(4,k,i)
        tb21v=gmiData%gmiS1(5,k,i)
        tb37v=gmiData%gmiS1(6,k,i)
        tb37h=gmiData%gmiS1(7,k,i)
        tb85h=gmiData%gmiS1(9,k,i)
        do i1=1,5 !4/14/14 MG Begin
           do j1=1,5
              if(gmiData%S1lat(k,i)>-999) then
                 call getwfraction(gmiData%S1lat(k,i)+(i1-2)*0.15,&
                      gmiData%S1lon(k,i)+(j1-2)*0.15,wfract)
                 iLandSea(i1,j1)=igetlandsea(gmiData%S1lat(k,i)+(i1-2)*0.15,&
                      gmiData%S1lon(k,i)+(j1-2)*0.15, &
                      geoData%lsflag)
                 if(wfract>75) then
                    iLandSea(i1,j1)=0
                 else
                    iLandSea(i1,j1)=1
                 endif
              else
                 iLandSea(i1,j1)=-99
              endif
           enddo
        enddo !4/14/14 MG End
        row = ((gmiData%S1lon(k,i)+ 180.0)/2.5 )
        col = ((gmiData%S1lat(k,i)+ 90.0)/2.0)

!begin WSO 04/07/2013 
!...fix for lon = 180.0 or lat = 90.0
        if(row > 143) row = 143
        if(col > 89) col = 89
!end WSO 04/07/2013
        if(gmiData%S1lat(k,i)>-999) then
           sstv=0.1*geoData%sstdata(col+1,row+1,mm)-273.15
           hfreez=(0.1*geoData%sstdata(col+1,row+1,mm)-273.15)/6.5;
        else
           sstv=-999
           hfreez=-999
        endif
        if(maxval(iLandSea)==0 .and. gmiData%S1lat(k,i)>-999) then
           ilsea(i)=0
           rmsmin=10e5
           isol=-1
           do ic=1,1200
              rms=0
              do j=1,7
                 rms=rms+(gMIData%gmiS1(j,k,i)-sstclust(ic,j))**2
              enddo
              if(rmsmin>rms)  then
                 rmsmin=rms
                 isol=ic
              endif
           enddo
           if(isol>-1) then
              gMIData%sfc_wind(k,i)=sstclust(isol,8)
              gMIData%sst(k,i)=sstclust(isol,9)
              gMIData%tpw(k,i)=sstclust(isol,10)
              gMIData%tb10hC(k,i)=tb10h
              gMIData%tb10vC(k,i)=tb10v
              gMIData%tb19hC(k,i)=tb19h
              gMIData%tb19vC(k,i)=tb19v
              gMIData%tb21vC(k,i)=tb21v
              gMIData%tb37hC(k,i)=tb37h
              gMIData%tb37vC(k,i)=tb37v
 !             write(31,131) gMIData%sfc_wind(k,i), gMIData%sst(k,i), &
 !                  gMIData%tpw(k,i), tb10h,tb10v,tb19h,tb19v, tb37h, tb37v
           else
              !gMIData%sfc_wind(k,i)=0.
              gMIData%sst(k,i)=sstv
              !gMIData%tpw(k,i)=sstclust(1,10)
           endif
           if(isol>0 .and. tb10v>0) then
              call pix_screen_ocean( tb10h, tb10v,tb19h, tb19v, tb21v, &
                   tb37h, tb37v, tb85h, &
                   gmiData%S1lat(k,i),gmiData%S1lon(k,i),&
                   geoData%sstdata, &
                   mm, irainflag,  gMIData%sst(k,i), &
                   gMIData%sfc_wind(k,i),gMIData%tpw(k,i) ) 
           else
              irainflag=1
           endif
200        format(9(F6.2,1x))
           if(irainflag==0) then
              left(i)=0
              irflag(i)=0
              lc=0
           else
              if(lc >-1)  then
                 lc=lc+1
              endif
              left(i)=lc
              irflag(i)=1
           endif
           if(irainflag==1) then
              gMIData%sst(k,i)=sstv
           endif
        else
           if(iLandSea(2,2)==0) then
              irflag(i)=1
              ilsea(i)=0
           endif
           left(i)=-1
           right(i)=-1
           lc=-1
        endif
     enddo
!131 format(15(F7.2,1x))
     rc=-1
     do i=gmiData%n1b11,1,-1
        gMIData%landSea(k,i)=ilsea(i)
        gMIData%rainFlag(k,i)=irflag(i)
        if(ilsea(i)==0) then
           if(irflag(i)==0) then
              right(i)=0
              rc=0
           else
              if(rc>=0) then
                 rc=rc+1
              endif
              gMIData%sfc_wind(k,i)=0
              gMIData%sst(k,i)=0
              gMIData%tpw(k,i)=0
              gMIData%tb10hC(k,i)=0
              gMIData%tb10vC(k,i)=0
              gMIData%tb19hC(k,i)=0
              gMIData%tb19vC(k,i)=0
              gMIData%tb21vC(k,i)=0
              gMIData%tb37hC(k,i)=0
              gMIData%tb37vC(k,i)=0
              right(i)=rc
              if(right(i)>0) then
                 w1=1./right(i)
                 gMIData%sfc_wind(k,i)=w1*gMIData%sfc_wind(k,i+right(i))
                 gMIData%sst(k,i)=w1*gMIData%sst(k,i+right(i))
                 gMIData%tpw(k,i)=w1*gMIData%tpw(k,i+right(i))
                 gMIData%tb10hC(k,i)=w1*gMIData%tb10hC(k,i+right(i))
                 gMIData%tb10vC(k,i)=w1*gMIData%tb10vC(k,i+right(i))
                 gMIData%tb19hC(k,i)=w1*gMIData%tb19hC(k,i+right(i))
                 gMIData%tb19vC(k,i)=w1*gMIData%tb19vC(k,i+right(i))
                 gMIData%tb21vC(k,i)=w1*gMIData%tb21vC(k,i+right(i))
                 gMIData%tb37hC(k,i)=w1*gMIData%tb37hC(k,i+right(i))
                 gMIData%tb37vC(k,i)=w1*gMIData%tb37vC(k,i+right(i))
              else
                 w1=0.
                 gMIData%sfc_wind(k,i)=0.
                 gMIData%sst(k,i)=0.
                 gMIData%tpw(k,i)=0.
              endif
              if(left(i)>0) then
                 w2=1./left(i)
                 gMIData%sfc_wind(k,i)=gMIData%sfc_wind(k,i)+ &
                      w2*gMIData%sfc_wind(k,i-left(i))
                 gMIData%sst(k,i)=gMIData%sst(k,i)+ &
                      w2*gMIData%sst(k,i-left(i))
                 gMIData%tpw(k,i)=gMIData%tpw(k,i)+ &
                      w2*gMIData%tpw(k,i-left(i))
                 gMIData%tb10hC(k,i)=gMIData%tb10hC(k,i)+ &
                      w2*gMIData%tb10hC(k,i-left(i))
                 gMIData%tb10vC(k,i)=gMIData%tb10vC(k,i)+ &
                      w2*gMIData%tb10vC(k,i-left(i))
                 gMIData%tb19hC(k,i)=gMIData%tb19hC(k,i)+ &
                      w2*gMIData%tb19hC(k,i-left(i))
                 gMIData%tb19vC(k,i)=gMIData%tb19vC(k,i)+ &
                      w2*gMIData%tb19vC(k,i-left(i))
                 gMIData%tb21vC(k,i)=gMIData%tb21vC(k,i)+ &
                      w2*gMIData%tb21vC(k,i-left(i))
                 gMIData%tb37hC(k,i)=gMIData%tb37hC(k,i)+ &
                      w2*gMIData%tb37hC(k,i-left(i))
                 gMIData%tb37vC(k,i)=gMIData%tb37vC(k,i)+ &
                      w2*gMIData%tb37vC(k,i-left(i))
              else
                 w2=0.
              endif
              if(w1+w2>0) then
                 gMIData%sfc_wind(k,i)=gMIData%sfc_wind(k,i)/(w1+w2)
                 gMIData%sst(k,i)=gMIData%sst(k,i)/(w1+w2)
                 gMIData%tpw(k,i)=gMIData%tpw(k,i)/(w1+w2)
                 gMIData%tb10hC(k,i)=gMIData%tb10hC(k,i)/(w1+w2)
                 gMIData%tb10vC(k,i)=gMIData%tb10vC(k,i)/(w1+w2)
                 gMIData%tb19hC(k,i)=gMIData%tb19hC(k,i)/(w1+w2)
                 gMIData%tb19vC(k,i)=gMIData%tb19vC(k,i)/(w1+w2)
                 gMIData%tb21vC(k,i)=gMIData%tb21vC(k,i)/(w1+w2)
                 gMIData%tb37hC(k,i)=gMIData%tb37hC(k,i)/(w1+w2)
                 gMIData%tb37vC(k,i)=gMIData%tb37vC(k,i)/(w1+w2)
                 
              else
                 gMIData%sfc_wind(k,i)=-99
                 gMIData%sst(k,i)=-999
                 gMIData%tpw(k,i)=-99
                 gMIData%tb10hC(k,i)=-99
                 gMIData%tb10vC(k,i)=-99
                 gMIData%tb19hC(k,i)=-99
                 gMIData%tb19vC(k,i)=-99
                 gMIData%tb21vC(k,i)=-99
                 gMIData%tb37hC(k,i)=-99
                 gMIData%tb37vC(k,i)=-99
              endif
           endif
        else
           right(i)=-1
           rc=-1
           gMIData%sfc_wind(k,i)=-999
           gMIData%sst(k,i)=-999
           gMIData%tpw(k,i)=-999
           gMIData%tb10hC(k,i)=-999
           gMIData%tb10vC(k,i)=-999
           gMIData%tb19hC(k,i)=-999
           gMIData%tb19vC(k,i)=-999
           gMIData%tb21vC(k,i)=-999
           gMIData%tb37hC(k,i)=-999
           gMIData%tb37vC(k,i)=-999
        endif

     enddo
  enddo
  do k=70,150  !4/14/14 MG
     do i=1,gmiData%n1b11
        if(gmiData%S1lat(k,i)<-99) then
           gMIData%sfc_wind(k,i)=-999
           gMIData%sst(k,i)=-999
           gMIData%tpw(k,i)=-999
           gMIData%tb10hC(k,i)=-999
           gMIData%tb10vC(k,i)=-999
           gMIData%tb19hC(k,i)=-999
           gMIData%tb19vC(k,i)=-999
           gMIData%tb21vC(k,i)=-999
           gMIData%tb37hC(k,i)=-999
           gMIData%tb37vC(k,i)=-999
        endif
     enddo
  enddo !4/14/14 MG
!begin WSO 04/07/2013
!  write(*, '("in sst rend: ")')
!  do i=1,500
!    write(*, '(81i2)') (int(gMIData%sfc_wind(k, i)), k=70,150)
!  enddo
!  print*, gmiData%n1b11, maxval(gMIData%sfc_wind(70:150,1:500))
!  stop
!end WSO 04/07/2013
  close(10)
10 format(4(F6.2,1x),2(I5))
20 format(17(F6.2,1x))
end subroutine sst2


subroutine ocean_screening
end subroutine ocean_screening
