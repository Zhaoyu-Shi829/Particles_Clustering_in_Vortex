      SUBROUTINE PP_STOKES_DRAG(UPART,UPART_F,DIAM,GMOL,RHO,FACT_A_W,
     $   DRAG_ACCEL)
C     comptes the stokes drag acting on a particle
      
      IMPLICIT NONE
      REAL UPART(3),UPART_F(3),DIAM,GMOL,DRAG_ACCEL(3),RHO,FACT_A_W
      



      REAL VEL_DIFF(3), NORM_VEL_DIFF, RE_P, C_D_MAL_NORM_VEL_DIFF

      VEL_DIFF = UPART - UPART_F


      NORM_VEL_DIFF = SQRT(VEL_DIFF(1)**2 + VEL_DIFF(2)**2
     $   + VEL_DIFF(3)**2)
      RE_P = NORM_VEL_DIFF * DIAM / GMOL * RHO




C     C_D AUS CLIFT ET. AL

      IF (RE_P .LT. 0.01) THEN
C        C_D = 3.0/16.0 + 24/RE_P
         C_D_MAL_NORM_VEL_DIFF = 3.0/16.0*NORM_VEL_DIFF 
     $      + 24./DIAM*GMOL/RHO
      ELSEIF (RE_P .LT. 3.0E5) THEN
         C_D_MAL_NORM_VEL_DIFF = (24 / RE_P * (1+0.15*RE_P**0.687)
     $      + 0.42/(1+4.25E4/RE_P**1.16))*NORM_VEL_DIFF
      ELSEIF (RE_P .LT. 4.0E5) THEN
         C_D_MAL_NORM_VEL_DIFF = (29.78-5.3*LOG10(RE_P))*NORM_VEL_DIFF
      ELSEIF (RE_P .LT. 1.0E6) THEN
         C_D_MAL_NORM_VEL_DIFF = (0.1*LOG10(RE_P)-0.49)*NORM_VEL_DIFF
      ELSE
         C_D_MAL_NORM_VEL_DIFF= 0.19-8.0E4/RE_P*NORM_VEL_DIFF
      ENDIF



      DRAG_ACCEL = FACT_A_W
     $    * C_D_MAL_NORM_VEL_DIFF * VEL_DIFF
      
      
      RETURN
      END
