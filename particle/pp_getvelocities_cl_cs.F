      SUBROUTINE PP_GETVELOCITIES_CL_CS (KK,JJ,II,X,Y,Z,DX,DY,DZ,
     $      DDX,  DDY,  DDZ,
     $     RDX, RDY, RDZ,
     $     RDDX, RDDY, RDDZ,
     $                          U,V,W,ALPHA,BETA,GAMMA,DELTA,
     $                          NPART,XP,YP,ZP,
     $                          IINDEX,JINDEX,KINDEX,UPART)
C------------------------------------------------------------
C
C     GETS NEIGHBOURING VELOCITIES AROUND PARTICLES - COMPLETE LINEAR AND CONSERVATIVE IN DIFFERENTIAL SENSE
C
C     030306 (CG)  ORIGINAL
C
C------------------------------------------------------------
      IMPLICIT NONE
      INTEGER KK,JJ,II,NPART

      REAL U(KK,JJ,II), V(KK,JJ,II), W(KK,JJ,II),
     $     ALPHA(KK,JJ,II),BETA(KK,JJ,II),GAMMA(KK,JJ,II),
     $     DELTA(KK,JJ,II)


      REAL    X(II),    Y(JJ),    Z(KK),
     $       DX(II),   DY(JJ),   DZ(KK),
     $      DDX(II),  DDY(JJ),  DDZ(KK),
     $      RDX(II),  RDY(JJ),  RDZ(KK),
     $     RDDX(II), RDDY(JJ), RDDZ(KK)


      REAL XP(NPART), YP(NPART), ZP(NPART)
      REAL UPART(3,NPART)

      INTEGER IINDEX(NPART), JINDEX(NPART), KINDEX(NPART)


      INTEGER I,J,K, IPART

      IF (NPART .LE. 0) THEN
         RETURN
      ENDIF


C------------------------------- INTERPOLATION OF U

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1
               ALPHA(K,J,I) = (U(K,J,I)-U(K,J,I-1))*RDDX(I)
               BETA(K,J,I) = 0.25*((U(K,J+1,I)+U(K,J+1,I-1)
     $            -U(K,J,I)-U(K,J,I-1))*RDY(J)
     $            +(U(K,J,I)+U(K,J,I-1)
     $            -U(K,J-1,I)-U(K,J-1,I-1))*RDY(J-1))
               GAMMA(K,J,I) = 0.25*((U(K+1,J,I)+U(K+1,J,I-1)
     $            -U(K,J,I)-U(K,J,I-1))*RDZ(K)
     $            +(U(K,J,I)+U(K,J,I-1)
     $            -U(K-1,J,I)-U(K-1,J,I-1))*RDZ(K-1))
               DELTA(K,J,I) = 0.5*(U(K,J,I)+U(K,J,I-1)
     $            -ALPHA(K,J,I)*(DDX(I)-DX(I-1))
     $            -BETA(K,J,I)*(DDY(J) - DY(J-1))
     $            -GAMMA(K,J,I)*(DDZ(K) - DZ(K-1)))
            ENDDO
         ENDDO
      ENDDO


C$OMP PARALLEL DO PRIVATE(IPART,I,J,K)
      DO IPART=1,NPART
         I=IINDEX(IPART)
         J=JINDEX(IPART)
         K=KINDEX(IPART)

!#ifdef _DEBUG_
         IF ((I .LT. 3) .OR. (I .GT. II-2) .OR.
     $      (J .LT. 3) .OR. (J .GT. JJ-2) .OR.
     $      (K .LT. 3) .OR. (K .GT. KK-2)) THEN
            WRITE(6,*) IPART, XP(IPART), YP(IPART), ZP(IPART),
     $         I, II, J, JJ, K, KK
            CALL ERRR(501, 'GETVEL_CL_CS')
         ENDIF
!#endif

         UPART(1,IPART)=ALPHA(K,J,I)*(XP(IPART)-X(I))
     $      + BETA(K,J,I)*(YP(IPART)-Y(J))
     $      +GAMMA(K,J,I)*(ZP(IPART)-Z(K))+DELTA(K,J,I)
      ENDDO
C$OMP END PARALLEL DO




C------------------------------- INTERPOLATION OF V

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1
               ALPHA(K,J,I) = (V(K,J,I)-V(K,J-1,I))*RDDY(J)
               BETA(K,J,I) = 0.25*((V(K,J,I+1)+V(K,J-1,I+1)
     $            -V(K,J,I)-V(K,J-1,I))*RDX(I)
     $            +(V(K,J,I)+V(K,J-1,I)
     $            -V(K,J,I-1)-V(K,J-1,I-1))*RDX(I-1))
               GAMMA(K,J,I) = 0.25*((V(K+1,J,I)+V(K+1,J-1,I)
     $            -V(K,J,I)-V(K,J-1,I))*RDZ(K)
     $            +(V(K,J,I)+V(K,J-1,I)
     $            -V(K-1,J,I)-V(K-1,J-1,I))*RDZ(K-1))
               DELTA(K,J,I) = 0.5*(V(K,J,I)+V(K,J-1,I)
     $            -ALPHA(K,J,I)*(DDY(J)-DY(J-1))
     $            -BETA(K,J,I)*(DDX(I) - DX(I-1))
     $            -GAMMA(K,J,I)*(DDZ(K) - DZ(K-1)))
            ENDDO
         ENDDO
      ENDDO



C$OMP PARALLEL DO PRIVATE(IPART,I,J,K)
      DO IPART=1,NPART
         I=IINDEX(IPART)
         J=JINDEX(IPART)
         K=KINDEX(IPART)

         UPART(2,IPART)=ALPHA(K,J,I)*(YP(IPART)-Y(J))
     $      + BETA(K,J,I)*(XP(IPART)-X(I))
     $      +GAMMA(K,J,I)*(ZP(IPART)-Z(K))+DELTA(K,J,I)
      ENDDO
C$OMP END PARALLEL DO





C------------------------------- INTERPOLATION OF W

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1
               ALPHA(K,J,I) = (V(K,J,I)-V(K,J-1,I))*RDDY(J)
               BETA(K,J,I) = 0.25*((V(K,J,I+1)+V(K,J-1,I+1)
     $            -V(K,J,I)-V(K,J-1,I))*RDX(I)
     $            +(V(K,J,I)+V(K,J-1,I)
     $            -V(K,J,I-1)-V(K,J-1,I-1))*RDX(I-1))
               GAMMA(K,J,I) = 0.25*((V(K+1,J,I)+V(K+1,J-1,I)
     $            -V(K,J,I)-V(K,J-1,I))*RDZ(K)
     $            +(V(K,J,I)+V(K,J-1,I)
     $            -V(K-1,J,I)-V(K-1,J-1,I))*RDZ(K-1))
               DELTA(K,J,I) = 0.5*(V(K,J,I)+V(K,J-1,I)
     $            -ALPHA(K,J,I)*(DDX(I)-DX(I-1))
     $            -BETA(K,J,I)*(DDY(J) - DY(J-1))
     $            -GAMMA(K,J,I)*(DDZ(K) - DZ(K-1)))

               ALPHA(K,J,I) = (W(K,J,I)-W(K-1,J,I))*RDDZ(K)
               BETA(K,J,I) = 0.25*((W(K,J,I+1)+W(K-1,J,I+1)
     $            -W(K,J,I)-W(K-1,J,I))*RDX(I)
     $            +(W(K,J,I)+W(K-1,J,I)
     $            -W(K,J,I-1)-W(K-1,J,I-1))*RDX(I-1))
               GAMMA(K,J,I) = 0.25*((W(K,J+1,I)+W(K-1,J+1,I)
     $            -W(K,J,I)-W(K-1,J,I))*RDY(J)
     $            +(W(K,J,I)+W(K-1,J,I)
     $            -W(K,J-1,I)-W(K-1,J-1,I))*RDY(J-1))
               DELTA(K,J,I) = 0.5*(W(K,J,I)+W(K-1,J,I)
     $            -ALPHA(K,J,I)*(DDZ(K)-DZ(K-1))
     $            -BETA(K,J,I)*(DDX(I) - DX(I-1))
     $            -GAMMA(K,J,I)*(DDY(J) - DY(J-1)))
            ENDDO
         ENDDO
      ENDDO



C$OMP PARALLEL DO PRIVATE(IPART,I,J,K)
      DO IPART=1,NPART
         I=IINDEX(IPART)
         J=JINDEX(IPART)
         K=KINDEX(IPART)

         UPART(3,IPART)=ALPHA(K,J,I)*(ZP(IPART)-Z(K))
     $      + BETA(K,J,I)*(XP(IPART)-X(I))
     $      +GAMMA(K,J,I)*(YP(IPART)-Y(J))+DELTA(K,J,I)
      ENDDO
C$OMP END PARALLEL DO



      RETURN

      END









