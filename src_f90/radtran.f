
      SUBROUTINE RADTRAN(UMU, NLYR, TB, BTEMP, LYRTEMP, LYRHGT, KEXT,
     $     SALB, ASYM, FISOT, EMIS, EBAR, lambert, prnt)
cf2py  real, intent(out) :: tb
C
C     CHRIS KUMMEROW
C     INCLUDES ASYMPTOTIC EXPRESSIONS FOR TERM3, TERM4, AND TERM5 IN
C     THE LIMIT OF SMALL EFFECTIVE OPTICAL DEPTHS; BILL OLSON FEB, 1995.
C
      PARAMETER    ( MAXLYR = 90 )
      CHARACTER*1  POLN
      LOGICAL LAMBERT, PRNT(3)
      integer ier, imode, icode
      REAL  LYRTEMP(0:MAXLYR), LYRHGT(0:MAXLYR), KEXT(MAXLYR),
     $      SALB(MAXLYR), ASYM(MAXLYR), L(MAXLYR), H(MAXLYR),
     $      B0(MAXLYR), B1(MAXLYR), W(2*MAXLYR,2*MAXLYR),
     $      BB(2*MAXLYR), DP(MAXLYR), DM(MAXLYR), Z(0:MAXLYR)
      REAL  IOUT(0:MAXLYR), I_IN(MAXLYR+1,100), MU, NU
      real x(2*maxlyr)
C
c      write(6,'(1x,"umu:  ",f7.2,"  nlyr: ",i5,"  tb: ",f7.2,
c     + "  btemp: ",f7.2,"  fisot: ",f7.2,
c     + /1x,"emis: ",f10.3,"  ebar: ",f10.3,"  lamb: ",l5, 
c     + /1x,"  lyrtemp(0): ",f7.2,"  lyrhgt(0): ",f7.2)')
c     + umu,nlyr,tb,btemp,fisot,emis,ebar,lambert,lyrtemp(0),
c     + lyrhgt(0)
c      do j=nlyr,1,-1
c      write(6,'(1x,"lyrtemp: ",f7.2,"  lyrhgt: ",f7.2,
c     + "  kext: ",f10.4,"  salb: ",f10.4,"  asym: ",f10.4)')
c     + lyrtemp(j),lyrhgt(j),kext(j),salb(j),asym(j)
c      end do
C                                CALCULATE SOME CONSTANTS
      Z(0) = LYRHGT(0)
      DO 20  J = 1,NLYR
        Z(J)  = LYRHGT(J)
        B0(J) = LYRTEMP(J-1)
        B1(J) = (LYRTEMP(J) - LYRTEMP(J-1))/(LYRHGT(J) - LYRHGT(J-1))
        L(J) = SQRT(3.*KEXT(J)*KEXT(J)*(1. - SALB(J))*
     $         (1. - SALB(J)*ASYM(J)))
        H(J) = 1.5*KEXT(J)*(1. - SALB(J)*ASYM(J))
  20  CONTINUE


C     FILL IN THE MATRIX ELEMENTS WHICH FORM THE BOUNDARY CONDITIONS
C     AT THE TOP, BOTTOM AND LAYER INTERFACES OF THE CLOUD.  THERE ARE
C     TWO QUANTITIES, D+ "DP", AND D- "DM" AT EACH BOUNDARY, SO THE 
C     MATRIX HAS DIMENSIONS  2*NLYR BY 2*NLYR.
C     ORDER OF D'S:  D+(1),D-(1),D+(2),D-(2), .... , D+(NLYR),D-(NLYR)

C     SET ALL MATRIZ ELEMENTS TO ZERO	
      DO 45  I = 1,2*NLYR
       DO 45  J = 1,2*NLYR
        W(I,J) = 0.0
  45  CONTINUE	

C     FILL IN THE NON-ZERO MATRIX ELEMENTS
      W(1,1)   = ((EBAR - 2.)*L(1)/H(1)) + EBAR
      W(1,2)   = ((2. - EBAR)*L(1)/H(1)) + EBAR
      DO 50  I = 2,2*(NLYR-1),2
       W(I,I-1)   =  (1. - L(I/2)/H(I/2))*EXP(+L(I/2)*(Z(I/2)-Z(I/2-1)))
       W(I,I  )   =  (1. + L(I/2)/H(I/2))*EXP(-L(I/2)*(Z(I/2)-Z(I/2-1)))
       W(I,I+1)   = -(1. - L(I/2+1)/H(I/2+1))
       W(I,I+2)   = -(1. + L(I/2+1)/H(I/2+1))

       W(I+1,I-1) =  (1. + L(I/2)/H(I/2))*EXP(+L(I/2)*(Z(I/2)-Z(I/2-1)))
       W(I+1,I)   =  (1. - L(I/2)/H(I/2))*EXP(-L(I/2)*(Z(I/2)-Z(I/2-1)))
       W(I+1,I+1) = -(1. + L(I/2+1)/H(I/2+1))
       W(I+1,I+2) = -(1. - L(I/2+1)/H(I/2+1))
c        write(*,*) (w(i,i+j),j=-2,2)
  50  CONTINUE
      W(2*NLYR,2*NLYR-1) =  (1. + L(NLYR)/H(NLYR))*EXP(+L(NLYR)*
     $                       (Z(NLYR)-Z(NLYR-1)))
      W(2*NLYR,2*NLYR)   =  (1. - L(NLYR)/H(NLYR))*EXP(-L(NLYR)*
     $                       (Z(NLYR)-Z(NLYR-1)))

C     FILL IN THE ROW OF CONSTANTS IN THE LINEAR EQUATIONS
      BB(1)    = EBAR*BTEMP - EBAR*B0(1) - (EBAR - 2.)*B1(1)/H(1)
      DO 55  I = 2,2*(NLYR-1),2
       BB(I)   =  + B1(I/2)/H(I/2) - B1(I/2+1)/H(I/2+1)
       BB(I+1) =  - B1(I/2)/H(I/2) + B1(I/2+1)/H(I/2+1)
  55  CONTINUE
      BB(2*NLYR)  =  FISOT - B0(NLYR) - B1(NLYR)*(Z(NLYR) 
     $               - Z(NLYR-1) + 1/H(NLYR))
C
C
C     MATRIX INVERSION IN DONE IN SUBROUTINE LINPAK
c      write(*,*) bb
c      write(*,*) b1
c      stop
      imode=0
      call band(imode,2*nlyr,w,bb,icode)
!      call linsysa(w,bb,x,2*nlyr,2*maxlyr,ier)
      
!     CALL LINPAK(NLYR, W, BB, RCOND)
!      do i=1,2*nlyr
!         write(*,*) x(i), bb(i)
!      enddo
!      stop
!

      do i=1,2*nlyr
!         bb(i)=x(i)
      enddo
c     write(*,*) bb
c     stop
C     
c     WRITE(*,*) ' '
c     IF ( PRNT(2) ) WRITE(4,*) ' '
      DO 60  I = 1,NLYR
       DP(I) = BB(2*I-1)
       DM(I) = BB(2*I)
c       WRITE(6,653) I, DP(I), I, DM(I)      
c       IF ( PRNT(2) )  WRITE(4,653) I, DP(I), I, DM(I)
  60  CONTINUE
 653  FORMAT(10X,'D+ (',I2,') = ',E11.4,6X,'D-(',I2,') = ',E11.4)

C     AFTER D'S ARE KNOWN, CALCULATE SURFACE RADIANCE

      MU = UMU
      NU = -UMU
      lambert = .false.
C     FOR THE FOLLOWING CALCULATIONS, REFER TO APPENDIX B OF THESIS
C     *********************************************************************
      
      IF ( LAMBERT ) THEN 
C      CALCULATE THE DOWNWELLING FLUX THROUGH THE ATMOSPHERE AT 81 ANGLES
       NANG = 81
       DO 997  NN = 1,NANG
        XNU = -(2.*NN - 1.)/(NANG*2.)
        I_IN(NLYR+1,NN) = FISOT
C       LOOP THROUGH THE REMAINING LAYERS
        DO 100  J = NLYR,1,-1
C        CALCULATE RADIANCE FROM TOP OF LAYER "J"
         XA = B0(J) - 1.5*SALB(J)*ASYM(J)*XNU*B1(J)/H(J)
         XB = B1(J)
         XC = SALB(J)*DP(J)*(1. - 1.5*ASYM(J)*XNU*L(J)/H(J))               
         XD = SALB(J)*DM(J)*(1. + 1.5*ASYM(J)*XNU*L(J)/H(J))
         YA = KEXT(J)/XNU
         YB = YA + L(J)
         YC = YA - L(J)
         DZ = Z(J) - Z(J-1)

         TERM1 = I_IN(J+1,NN)*EXP(YA*DZ)
         TERM2 = XA*(1. - EXP(YA*DZ))
c     new tests
         if(abs(ya*dz) .lt. 1.e-5) then
           term3=-xb*ya*dz*dz
         else
           TERM3 = XB/YA*(EXP(YA*DZ)*(1. - YA*DZ) - 1.)
         end if
         if(abs(yb*dz) .lt. 1.e-5) then
           term4=-xc*ya*dz
         else
           TERM4 = XC*YA/YB*(1. - EXP(YB*DZ))
         end if
         if(abs(yc*dz) .lt. 1.e-5) then
           term5=-xd*ya*dz
         else
           TERM5 = XD*YA/YC*(1. - EXP(YC*DZ))
         end if
         I_IN(J,NN) = TERM1 + TERM2 + TERM3 + TERM4 + TERM5
c        write(6,'(1x,"lev: ",i5,"  ang: ",i5,"  term1: ",f10.4,
c     +            "  term2: ",f10.4,"  term3: ",f10.4,"  term4: ",f10.4,
c     +            "  term5: ",f10.4,"  tbdn: ",f7.2)')
c     +   j,nn,term1,term2,term3,term4,term5,i_in(j,nn)
  100   CONTINUE
  997  CONTINUE
C
C      CALCULATE THE TOTAL DOWNWELLING FLUX REACHING THE SURFACE
       XIUP = 0.
       DO 47  NN = 1,NANG
        XIUP = XIUP + I_IN(1,NN)*(1./NANG)*(2.*NN-1.)/(2.*NANG)
  47   CONTINUE
       XIUP = 2.*XIUP
c       write(6,'(1x,"xiup: ",f7.2)') xiup
C
      ELSE

C      CALCULATE THE DOWNWELLING FLUX THROUGH THE ATMOSPHERE AT ANGLE MU
       NN = 22                      ! THIS IS A DUMMY INDEX
       I_IN(NLYR+1,NN) = FISOT
C      LOOP THROUGH THE REMAINING LAYERS
       DO 110  J = NLYR,1,-1
C       CALCULATE RADIANCE FROM TOP OF LAYER "J"
        XA = B0(J) - 1.5*SALB(J)*ASYM(J)*NU*B1(J)/H(J)
        XB = B1(J)
        XC = SALB(J)*DP(J)*(1. - 1.5*ASYM(J)*NU*L(J)/H(J))               
        XD = SALB(J)*DM(J)*(1. + 1.5*ASYM(J)*NU*L(J)/H(J))
        YA = KEXT(J)/NU
        YB = YA + L(J)
        YC = YA - L(J)
        DZ = Z(J) - Z(J-1)

        TERM1 = I_IN(J+1,NN)*EXP(YA*DZ)
        TERM2 = XA*(1. - EXP(YA*DZ))
c     new tests
        if(abs(ya*dz) .lt. 1.e-5) then
          term3=-xb*ya*dz*dz
        else
          TERM3 = XB/YA*(EXP(YA*DZ)*(1. - YA*DZ) - 1.)
        end if
        if(abs(yb*dz) .lt. 1.e-5) then
          term4=-xc*ya*dz
        else
          TERM4 = XC*YA/YB*(1. - EXP(YB*DZ))
        end if
        if(abs(yc*dz) .lt. 1.e-5) then
          term5=-xd*ya*dz
        else
          TERM5 = XD*YA/YC*(1. - EXP(YC*DZ))
        end if
        I_IN(J,NN) = TERM1 + TERM2 + TERM3 + TERM4 + TERM5
  110  CONTINUE
       XIUP = I_IN(1,22)
C
      ENDIF

C    
C
      IOUT(0) = EMIS*BTEMP + (1. - EMIS)*XIUP
      DO 101  J = 1,NLYR
C      CALCULATE THE UPWELLING RADIANCES AT THE TOP OF EACH LAYER J
       XA = B0(J) - 1.5*SALB(J)*ASYM(J)*MU*B1(J)/H(J)
       XB = B1(J)
       XC = SALB(J)*DP(J)*(1. - 1.5*ASYM(J)*MU*L(J)/H(J))               
       XD = SALB(J)*DM(J)*(1. + 1.5*ASYM(J)*MU*L(J)/H(J))
       YA = KEXT(J)/MU
       YB = YA + L(J)
       YC = YA - L(J)
       DZ = Z(J) - Z(J-1)
 
       TERM1 = IOUT(J-1)*EXP(-YA*DZ)
       TERM2 = XA*(1. - EXP(-YA*DZ))
c     new tests
       if(abs(ya*dz) .lt. 1.e-5) then
         term3=0.
       else
         TERM3 = XB/YA*(YA*DZ - 1. + EXP(-YA*DZ))
       end if
       if(abs(yb*dz) .lt. 1.e-5) then
         term4=xc*ya*dz*exp(-ya*dz)
       else
         TERM4 = XC*YA/YB*(EXP( (YB-YA)*DZ ) - EXP(-YA*DZ) )
       end if
       if(abs(yc*dz) .lt. 1.e-5) then
         term5=xd*ya*dz*exp(-ya*dz)
       else
         TERM5 = XD*YA/YC*EXP(-YA*DZ)*(EXP(YC*DZ) - 1.)
       end if
       IOUT(J) = TERM1 + TERM2 + TERM3 + TERM4 + TERM5
  101 CONTINUE
C
C
c      WRITE(*,*) ' '
c      IF ( PRNT(2) )  WRITE(4,*) ' '
      DO 44  J = 0,NLYR
c       WRITE(*,77) NLYR+1-J, I_IN(NLYR+1-J,22), NLYR-J, IOUT(NLYR-J)
c       IF ( PRNT(2) ) WRITE(4,77) NLYR+1-J, I_IN(NLYR+1-J,22), 
c     $                NLYR-J, IOUT(NLYR-J)
  44  CONTINUE
  77  FORMAT(10X,'Iin(',I2,') = ',F7.1,10X,'Iout(',I2,') = ',F7.1)

C
      TB = IOUT(NLYR)
!      write(*,*) tb
      RETURN
      END
