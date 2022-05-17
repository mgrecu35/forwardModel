      SUBROUTINE GCLOUD(FREQY,T,CLW,Z_CLW)
cf2py real, intent(out) :: z_clw      
      ES = 191.4 - 0.378*T
      TAU = 0.00199*EXP(2140/T)/T
      X = 6.283*FREQY*TAU
      EZ = 4.9
      RPE = EZ + (ES-EZ)/(1+X*X)
      XIPE = -(EZ -ES)*X/(1+X*X)
      Z_CLW = 0.188*CLW*FREQY*XIPE/((RPE+2)**2 + (XIPE)**2)
      RETURN
      END
      

