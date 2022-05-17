      SUBROUTINE  EMIT( F, NPOL, TS, W, UMU, EMIS, EBAR )
cf2py real, intent(out) :: emis, ebar
      real f, ts, w, umu, emis, ebar
      PARAMETER ( NANG = 21 )
      REAL *4 ANG(NANG), MU(NANG), ESUM(NANG)

      PI = 2.*ASIN(1.0)
      S = 45.                                    ! SALINITY IN PPM
      ANGLE = ACOS(UMU)*180./PI
C      
C     CALCULATE EMIS AT GIVEN ANGLE
      CALL EMISS (F,ANGLE,S,TS,W,EV,EH,EMISH,EMISV)
      if(abs(f-37)<3) then
         eh=eh-0.002
         ev=ev-0.002
      endif
      IF ( NPOL .EQ. 0 ) EMIS = EH
      IF ( NPOL .EQ. 1 ) EMIS = EV
C
C     CALCULATE EMIS AT VARIOUS ANGLES
      DO 58  I = 1,NANG
       ANG(I) = 4.*( I - 1 )
       MU(I) = COS(ANG(I)*PI/180.)
       ANGLES = ANG(I)
       CALL EMISS (F,ANGLES,S,TS,W,EV,EH,EMISH,EMISV)
       if(abs(f-37)<3) then
          eh=eh-0.002
          ev=ev-0.002
       endif
       ESUM(I) =  EV + EH 
  58  CONTINUE
C
C     CALCULATE EBAR
      SUM = 0.0
      DO 59  I = 1,NANG-1
       EAVG = 0.5*( ESUM(I) + ESUM(I+1) )
       DMU = MU(I) - MU(I+1)
       AVMU = 0.5*( MU(I) + MU(I+1) )
       SUM = SUM + EAVG*AVMU*DMU
  59  CONTINUE
      EBAR = SUM
      RETURN
      END
C
C
      SUBROUTINE  DIECON(S,T,FREQ,E1,E2)

      save sold, told, st, s2, t2, sst, stt, sstt, es, tau, sigma
      DATA TWOPI /6.283185307/ , EINF /4.9/ , SOLD /0.0/ , TOLD /-99./
      IF (S .EQ. SOLD .AND. T .EQ. TOLD) GO TO 10
      ST = S*T
      S2 = S*S
      T2 = T*T
      SST = S2*T
      STT = T2*S
      SSTT = S2*T2
      ES = 88.-4.339E-01*S+1.71E-03*S2-4.035E-01*T+8.065E-04*T2+6.170
     $  E-03 * ST-8.910E-05*SST-6.934E-05*STT+1.439E-06*SSTT
C
      TAU = (18.70-7.924E-02*S+6.35E-04*S2-5.489E-01*T+5.758E-03*T2+
     $1.889E-03*ST-7.209E-06*SST-5.299E-07*STT-2.101E-07*SSTT)*1.0E-12
C
      SIGMA = (7.788E-03*S-1.672E-06*S2-8.570E-15*T+2.996E-16*T2+4.059E
     $     -04 * ST-3.215E-06*SST-1.423E-06*STT+3.229E-08*SSTT)*1.0E11
C
   10 ZNU = FREQ*1.E09
      OMEGA = TWOPI*ZNU
      DEN = 1. + (OMEGA * TAU) ** 2
      E1 = (ES-EINF)/DEN+EINF
      E2 = (ES-EINF)*OMEGA*TAU/DEN+2.*SIGMA/ZNU
      SOLD = S
      TOLD = T
C
      RETURN
      END
C
C
      SUBROUTINE  EMISS(F,ANGLE,S,TS,W,EV,EH,EMISH,EMISV)
      DATA DTR / 0.01745329252 /
      T = TS-273.16
      THETA = ANGLE*DTR
C
      CALL DIECON(S,T,F,E1,E2)
C
      CMHU = COS (THETA)
      CMHU2 = CMHU*CMHU
      FACT1 = E1+CMHU2-1.
      B = (FACT1*FACT1+E2*E2)**0.25
      B2 = B*B
      ARG = E2/FACT1
      PSIZ = 0.5*ATAN(ARG)
      COSPZ = COS(PSIZ)*CMHU*B
      F11 = CMHU2+B2+2.*COSPZ
      PSIY = ATAN(E2/E1)
      G = CMHU*SQRT(E1*E1+E2*E2)
      COSPY = COS(PSIY-PSIZ)*G*B
      F22 = G*G+B2+2.*COSPY
C
C     FOR SPECULAR SURFACES THE HORIZONTAL AND VERTICAL EMISSIVITY
C     CALCULATION
C
      EMISH = 4.*COSPZ/F11
      EMISV = 4.*COSPY/F22
C
C     FROM HERE THE EFFECT OF SURFACE ROUGHNESS AND FOAM ARE INCLUDED
C     BASED ON HOLLINGER MODEL AND FOAM MODEL OF STOGRYN
C
      RSH = 1.-EMISH
      RSV = 1.-EMISV
C
C     P0 = 1.707476E-2+8.560329E-4*F+1.120024E-5*F*F
C     P1 = -1.500792E-2+1.820672E-3*F-4.633806E-5*F*F
C     P2 = 2.442217E-4-2.282022E-6*F+4.194212E-7*F*F
C
C     FOAM = (P0+P1*W+P2*W*W)
C
      FOAM = 10.751E-06 * W ** 3.231
C
      GH = 1.-1.748E-3*ANGLE-7.336E-5*ANGLE**2+1.044E-7*ANGLE**3
      GV = 1.-9.946E-4*ANGLE+3.218E-5*ANGLE**2-1.187E-6*ANGLE**3
     $   +7.0E-20*ANGLE**10
C
      A1 = (208.0+1.29*F)/TS
C
C     RFV = 1.-A1*GV-0.005*F
C     RFH = 1.-A1*GH-0.005*F
C
      RFV = 1. - A1 * GV
      RFH = 1. - A1 * GH
C
      Y = 7.32E-02*ANGLE
C
C     TS SURFACE TEMP IS IN DEGREE KELVIN
C
      SQRTF = SQRT(F)
C
C     CORRV = (W*(1.17E-01-2.09E-03*EXP(Y))*SQRTF/TS)-0.00065*F
C     CORRH = (W*(1.15E-01+3.80E-05*ANGLE**2)*SQRTF/TS)-0.00065*F
C
      CORRV = (W*(1.17E-01-2.09E-03*EXP(Y))*SQRTF/TS)
      CORRH = (W*(1.15E-01+3.80E-05*ANGLE**2)*SQRTF/TS)
C
      RRV = RSV-CORRV
      RRH = RSH-CORRH
C
      RV = RRV*(1.-FOAM)+RFV*FOAM
      RH = RRH*(1.-FOAM)+RFH*FOAM
C
      EH = 1.-RH
      EV = 1.-RV
C
      RETURN
      END


      SUBROUTINE  EMIT2( F, NPOL, TS, W, UMU, EMIS, EBAR )

      PARAMETER ( NANG = 21 )
      REAL  ANG(NANG), MU(NANG), ESUM(NANG)

      PI = 2.*ASIN(1.0)
      S = 35.                                    ! SALINITY IN PPM
      ANGLE = ACOS(UMU)*180./PI
C      
C     CALCULATE EMIS AT GIVEN ANGLE
      CALL EMISS2 (F,ANGLE,S,TS,W,EV,EH,EMISH,EMISV)
      IF ( NPOL .EQ. 0 ) EMIS = EH
      IF ( NPOL .EQ. 1 ) EMIS = EV
C
C     CALCULATE EMIS AT VARIOUS ANGLES
      DO 58  I = 1,NANG
       ANG(I) = 4.*( I - 1 )
       MU(I) = COS(ANG(I)*PI/180.)
       ANGLES = ANG(I)
       CALL EMISS2 (F,ANGLES,S,TS,W,EV,EH,EMISH,EMISV)
       ESUM(I) =  EV + EH 
  58  CONTINUE
C
C     CALCULATE EBAR
      SUM = 0.0
      DO 59  I = 1,NANG-1
       EAVG = 0.5*( ESUM(I) + ESUM(I+1) )
       DMU = MU(I) - MU(I+1)
       AVMU = 0.5*( MU(I) + MU(I+1) )
       SUM = SUM + EAVG*AVMU*DMU
  59  CONTINUE
      EBAR = SUM
      RETURN
      END
C
C
      SUBROUTINE  DIECON2(S,T,FREQ,E1,E2)

      save sold, told, st, s2, t2, sst, stt, sstt, es, tau, sigma

      DATA TWOPI /6.283185307/ , EINF /4.9/ , SOLD /0.0/ , TOLD /-99./
      IF (S .EQ. SOLD .AND. T .EQ. TOLD) GO TO 10
      ST = S*T
      S2 = S*S
      T2 = T*T
      SST = S2*T
      STT = T2*S
      SSTT = S2*T2
      ES = 88.-4.339E-01*S+1.71E-03*S2-4.035E-01*T+8.065E-04*T2+6.170
     $  E-03 * ST-8.910E-05*SST-6.934E-05*STT+1.439E-06*SSTT
C
      TAU = (18.70-7.924E-02*S+6.35E-04*S2-5.489E-01*T+5.758E-03*T2+
     $1.889E-03*ST-7.209E-06*SST-5.299E-07*STT-2.101E-07*SSTT)*1.0E-12
C
      SIGMA = (7.788E-03*S-1.672E-06*S2-8.570E-15*T+2.996E-16*T2+4.059E
     $     -04 * ST-3.215E-06*SST-1.423E-06*STT+3.229E-08*SSTT)*1.0E11
C
   10 ZNU = FREQ*1.E09
      OMEGA = TWOPI*ZNU
      DEN = 1. + (OMEGA * TAU) ** 2
      E1 = (ES-EINF)/DEN+EINF
      E2 = (ES-EINF)*OMEGA*TAU/DEN+2.*SIGMA/ZNU
      SOLD = S
      TOLD = T
C
      RETURN
      END
C
C
      SUBROUTINE  EMISS2(F,ANGLE,S,TS,W,EV,EH,EMISH,EMISV)
      DATA DTR / 0.01745329252 /
      T = TS-273.16
      THETA = ANGLE*DTR
C
      CALL DIECON2(S,T,F,E1,E2)
C
      CMHU = COS (THETA)
      CMHU2 = CMHU*CMHU
      FACT1 = E1+CMHU2-1.
      B = (FACT1*FACT1+E2*E2)**0.25
      B2 = B*B
      ARG = E2/FACT1
      PSIZ = 0.5*ATAN(ARG)
      COSPZ = COS(PSIZ)*CMHU*B
      F11 = CMHU2+B2+2.*COSPZ
      PSIY = ATAN(E2/E1)
      G = CMHU*SQRT(E1*E1+E2*E2)
      COSPY = COS(PSIY-PSIZ)*G*B
      F22 = G*G+B2+2.*COSPY
C
C     FOR SPECULAR SURFACES THE HORIZONTAL AND VERTICAL EMISSIVITY
C     CALCULATION
C
      EMISH = 4.*COSPZ/F11
      EMISV = 4.*COSPY/F22
C
C     FROM HERE THE EFFECT OF SURFACE ROUGHNESS AND FOAM ARE INCLUDED
C     BASED ON HOLLINGER MODEL AND FOAM MODEL OF STOGRYN
C
      RSH = 1.-EMISH
      RSV = 1.-EMISV
C
C     P0 = 1.707476E-2+8.560329E-4*F+1.120024E-5*F*F
C     P1 = -1.500792E-2+1.820672E-3*F-4.633806E-5*F*F
C     P2 = 2.442217E-4-2.282022E-6*F+4.194212E-7*F*F
C
C     FOAM = (P0+P1*W+P2*W*W)
C
      FOAM = 7.751E-06 * W ** 3.231
C
      GH = 1.-1.748E-3*ANGLE-7.336E-5*ANGLE**2+1.044E-7*ANGLE**3
      GV = 1.-9.946E-4*ANGLE+3.218E-5*ANGLE**2-1.187E-6*ANGLE**3
     $   +7.0E-20*ANGLE**10
C
      A1 = (208.0+1.29*F)/TS
C
C     RFV = 1.-A1*GV-0.005*F
C     RFH = 1.-A1*GH-0.005*F
C
      RFV = 1. - A1 * GV
      RFH = 1. - A1 * GH
C
      Y = 7.32E-02*ANGLE
C
C     TS SURFACE TEMP IS IN DEGREE KELVIN
C
      SQRTF = SQRT(F)
C
C     CORRV = (W*(1.17E-01-2.09E-03*EXP(Y))*SQRTF/TS)-0.00065*F
C     CORRH = (W*(1.15E-01+3.80E-05*ANGLE**2)*SQRTF/TS)-0.00065*F
C
      CORRV = (W*(1.17E-01-2.09E-03*EXP(Y))*SQRTF/TS)
      CORRH = (W*(1.15E-01+3.80E-05*ANGLE**2)*SQRTF/TS)
C
      RRV = RSV-CORRV
      RRH = RSH-CORRH
C
      RV = RRV*(1.-FOAM)+RFV*FOAM
      RH = RRH*(1.-FOAM)+RFH*FOAM
C
      EH = 1.-RH
      EV = 1.-RV
C
      RETURN
      END
