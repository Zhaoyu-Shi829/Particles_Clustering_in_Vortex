      SUBROUTINE PP_NEW_PARTS(
     $   IDIM3D,IDIM2D, IDIM1D, 
     $   II,JJ,KK,
     $   X,Y,Z,
     $   DX,DY,DZ,
     $   DDX,DDY,DDZ,
     $   RDX, RDY, RDZ,
     $   RDDX,RDDY,RDDZ,
     $   NXGRAE,NYGRAE,NZGRAE,
     $   XMIN, XMAX, YMIN, YMAX,
     $   ZMIN, ZMAX,
     $   U,V,W,
     $   HILF1,HILF2,HILF3,HILF4,
     $   DT,
     $   GMOL, RHO, 
     $   MAXBOCONDS,MYNBOCONDS,MYITYPBOCONDS,
     $   MYLBOGRIDS,
     $   FLOWTYP,
     $   XBANF,XBEND,
     $   YBANF,YBEND,
     $   ZBANF,ZBEND,
     $   AUB,AVB,AWB,
     $   NPART_GRID, 
     $   NPART_MAX, 
     $   PP_ARRSIZE2D,
     $   XPART,
     $   IINDEX,JINDEX,KINDEX,
     $   IINDEX_N,JINDEX_N,KINDEX_N,
     $   START_VEC,
     $   PART_IDX, 
     $   DIAMCLASS,
     $   UPART,
     $   UPART_F,
     $   OLD_UPART_F,
     $   PP_FLOWCASE,
     $   DS
#if defined _IB_
     $   ,BP
     $   ,PP_BODY_NCELLS_MAX, PP_BODY_NVEC
     $   ,PP_BODY_NVECNO
#endif
     $   ,NDIAM_CLASSES,IGRID)

      
      IMPLICIT NONE


      INTEGER IDIM3D,IDIM2D,IDIM1D
      INTEGER II,JJ,KK,IGRID,NPART_MAX
      INTEGER NXGRAE,NYGRAE,NZGRAE
      REAL X(II),Y(JJ),Z(KK)
      REAL DX(II),DY(JJ),DZ(KK)
      REAL DDX(II),DDY(JJ),DDZ(KK)
      REAL RDX(II), RDY(JJ), RDZ(KK)
      REAL RDDX(II),RDDY(JJ),RDDZ(KK)
      REAL XMIN, XMAX, YMIN, YMAX
      REAL ZMIN, ZMAX
      REAL U(KK,JJ,II),V(KK,JJ,II),W(KK,JJ,II)
      REAL HILF1(KK,JJ,II),HILF2(KK,JJ,II),
     $     HILF3(KK,JJ,II),HILF4(KK,JJ,II)

      REAL DT
      REAL GMOL, RHO
      INTEGER MAXBOCONDS,MYNBOCONDS(9), MYITYPBOCONDS(MAXBOCONDS,9)
      INTEGER MYLBOGRIDS(MAXBOCONDS,9)
      CHARACTER (LEN=16) FLOWTYP(MAXBOCONDS,9)
      REAL XBANF(MAXBOCONDS,9),XBEND(MAXBOCONDS,9),
     $    YBANF(MAXBOCONDS,9),YBEND(MAXBOCONDS,9),
     $    ZBANF(MAXBOCONDS,9),ZBEND(MAXBOCONDS,9),
     $    AUB(MAXBOCONDS,9),AVB(MAXBOCONDS,9),AWB(MAXBOCONDS,9)

      INTEGER NPART_GRID
      INTEGER PP_ARRSIZE2D,NDIAM_CLASSES
      REAL XPART(3,NPART_MAX)
      INTEGER IINDEX(NPART_GRID),JINDEX(NPART_GRID),KINDEX(NPART_GRID)
      INTEGER IINDEX_N(NPART_GRID),JINDEX_N(NPART_GRID),
     $   KINDEX_N(NPART_GRID)

      REAL START_VEC(3,NPART_MAX)
      INTEGER PART_IDX(NPART_MAX)
      INTEGER DIAMCLASS(NPART_MAX)
      REAL UPART(3,NPART_MAX)
      REAL UPART_F(3,NPART_MAX)
      REAL OLD_UPART_F(3,NPART_MAX)
      INTEGER PP_FLOWCASE
      REAL DS(NPART_MAX)

#if defined _IB_
      REAL BP(KK,JJ,II)
      INTEGER PP_BODY_NCELLS_MAX
      REAL PP_BODY_NVEC(4,PP_BODY_NCELLS_MAX)
      INTEGER PP_BODY_NVECNO(IDIM3D)
#endif

C     LOCAL VARS

      INTEGER IPART,I,J,K,NPART_NEW_FRO,
     $ NEW_PARTS, 
C     $ NEW_PARTS_TOP,NEW_PARTS_BOT,
C     $ NEW_PARTS_LFT,NEW_PARTS_RGT,
C     $ NEW_PARTS_FRO,NEW_PARTS_BAC,
     $ NPART_NEW,NPART_NEW_BAC,
     $ NPART_NEW_LFT,NPART_NEW_RGT,
     $ NPART_NEW_TOP,NPART_NEW_BOT

      INTEGER NEW_DIAMCLASS,NP,PP_NEWINDEX
      INTEGER KSTART(6),KSTOP(6),JSTART(6),JSTOP(6),
     $        ISTART(6),ISTOP(6)
      INTEGER NEW_IINDEX,NEW_JINDEX,NEW_KINDEX

      REAL X_NEW_PART,U_NEW_PART
      REAL VFLUX,AFLUX,ranfmglet
      
      REAL :: rrand(2)
      INTEGER :: idx1, idx2

      ALLOCATABLE
     $   NEW_PARTS(:),
C     $   NEW_PARTS_BOT(:),
C     $   NEW_PARTS_RGT(:),
C     $   NEW_PARTS_LFT(:),
C     $   NEW_PARTS_FRO(:),
C     $   NEW_PARTS_BAC(:),
     $   NEW_DIAMCLASS(:),
     $   NEW_IINDEX(:),
     $   NEW_JINDEX(:),
     $   NEW_KINDEX(:),
     $   X_NEW_PART(:,:),
     $   U_NEW_PART(:,:)

      INTEGER NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,L
      INTEGER IBOCOND_PP(6),COUNTER
      INTEGER COUNTERD1,COUNTERD2,COUNTERD3

      CALL PP_BASB(MAXBOCONDS, MYNBOCONDS,
     $   MYITYPBOCONDS,NFRO,NBAC,NRGT,NLFT,NBOT,NTOP)
      CALL PP_IBASB(MAXBOCONDS, MYNBOCONDS,
     $   MYITYPBOCONDS,IBOCOND_PP)

      NPART_NEW_FRO = 0
      NPART_NEW_BAC = 0
      NPART_NEW_LFT = 0
      NPART_NEW_RGT = 0
      NPART_NEW_BOT = 0
      NPART_NEW_TOP = 0
      NPART_NEW = 0
    
!      ALLOCATE(NEW_PARTS_TOP(PP_ARRSIZE2D))
!      ALLOCATE(NEW_PARTS_BOT(PP_ARRSIZE2D))
!      ALLOCATE(NEW_PARTS_RGT(PP_ARRSIZE2D))
      ALLOCATE(NEW_PARTS(PP_ARRSIZE2D))
!      ALLOCATE(NEW_PARTS(PP_ARRSIZE2D))
!      ALLOCATE(NEW_PARTS_BAC(PP_ARRSIZE2D))
      ALLOCATE(NEW_DIAMCLASS(PP_ARRSIZE2D))

      DO L=1,6
#if 0
        DO K=1,KK
         IF(Z(K).LE.ZBANF(IBOCOND_PP(L),L).AND.
     $      Z(K+1).GT.ZBANF(IBOCOND_PP(L),L)) KSTART(L) = K
         IF(Z(K).LE.ZBEND(IBOCOND_PP(L),L).AND.
     $      Z(K+1).GT.ZBEND(IBOCOND_PP(L),L)) KSTOP(L) = K
        ENDDO
        DO J=1,JJ
         IF(Y(J).LE.YBANF(IBOCOND_PP(L),L).AND.
     $      Y(J+1).GT.YBANF(IBOCOND_PP(L),L)) JSTART(L) = J
         IF(Y(J).LE.YBEND(IBOCOND_PP(L),L).AND.
     $      Y(J+1).GT.YBEND(IBOCOND_PP(L),L)) JSTOP(L) = J
        ENDDO
        DO I=1,II
         IF(X(I).LE.XBANF(IBOCOND_PP(L),L).AND.
     $      X(I+1).GT.XBANF(IBOCOND_PP(L),L)) ISTART(L) = I
         IF(X(I).LE.XBEND(IBOCOND_PP(L),L).AND.
     $      X(I+1).GT.XBEND(IBOCOND_PP(L),L)) ISTOP(L) = I
        ENDDO
#endif
        KSTART(L) = 3
        KSTOP(L) = KK-2
        
        JSTART(L) = 3
        JSTOP(L) = JJ-2
        
        ISTART(L) = 3
        ISTOP(L) = II-2
      ENDDO
C      WRITE(6,*)ISTART(1),ISTOP(1)
      IF (NFRO .EQ. 32) THEN
        VFLUX = 0.0
        AFLUX = 0.0

        DO K=KSTART(1),KSTOP(1)
         DO J=JSTART(1),JSTOP(1)
C          DO I=2,2
          DO I=ISTART(1),ISTOP(1)
          AFLUX = AFLUX + DDY(J)*DDZ(K)
          VFLUX = VFLUX + U(K,J,I)*DDY(J)*DDZ(K)*DT
          ENDDO
         ENDDO
        ENDDO
        NPART_NEW_FRO = INT(AUB(IBOCOND_PP(1),1)*VFLUX)
        WRITE(6,*) 'NPART_NEW_FRO: ',AFLUX,VFLUX,NPART_NEW_FRO,IGRID
      ENDIF
      IF (NBAC .EQ. 32) THEN
        VFLUX = 0.0
        AFLUX = 0.0

        DO K=KSTART(2),KSTOP(2)
         DO J=JSTART(2),JSTOP(2)
C          DO I=(II-3),(II-3)
          DO I=ISTART(2),ISTOP(2)
          AFLUX = AFLUX + DDY(J)*DDZ(K)
          VFLUX = VFLUX + U(K,J,I)*DDY(J)*DDZ(K)*DT
          ENDDO
         ENDDO
        ENDDO
        NPART_NEW_BAC = INT(AUB(IBOCOND_PP(2),2)*VFLUX)
C        WRITE(6,*) 'NPART_NEW_FRO: ',AFLUX,VFLUX,NPART_NEW_FRO,IGRID
      ENDIF
      IF (NRGT .EQ. 32) THEN
        VFLUX = 0.0
        AFLUX = 0.0

        DO K=KSTART(3),KSTOP(3)
C         DO J=2,2
         DO J=JSTART(3),JSTOP(3)
          DO I=ISTART(3),ISTOP(3)
          AFLUX = AFLUX + DDX(I)*DDZ(K)
          VFLUX = VFLUX + V(K,J,I)*DDX(I)*DDZ(K)*DT
          ENDDO
         ENDDO
        ENDDO
        NPART_NEW_RGT = INT(AVB(IBOCOND_PP(3),3)*VFLUX)
      ENDIF
      IF (NLFT .EQ. 32) THEN
        VFLUX = 0.0
        AFLUX = 0.0

        DO K=KSTART(4),KSTOP(4)
C         DO J=JJ-3,JJ-3
         DO J=JSTART(4),JSTOP(4)
          DO I=ISTART(4),ISTOP(4)
          AFLUX = AFLUX + DDX(I)*DDZ(K)
          VFLUX = VFLUX + V(K,J,I)*DDX(I)*DDZ(K)*DT
          ENDDO
         ENDDO
        ENDDO
        NPART_NEW_LFT = INT(AVB(IBOCOND_PP(4),4)*VFLUX)
      ENDIF
      IF (NBOT .EQ. 32) THEN
        VFLUX = 0.0
        AFLUX = 0.0

C        DO K=2,2
        DO K=KSTART(5),KSTOP(5)
         DO J=JSTART(5),JSTOP(5)
          DO I=ISTART(5),ISTOP(5)
          AFLUX = AFLUX + DDX(I)*DDY(J)
          VFLUX = VFLUX + W(K,J,I)*DDX(I)*DDY(J)*DT
          ENDDO
         ENDDO
        ENDDO
        NPART_NEW_BOT = INT(AWB(IBOCOND_PP(5),5)*VFLUX)
      ENDIF
      IF (NTOP .EQ. 32) THEN
        VFLUX = 0.0
        AFLUX = 0.0

C        DO K=KK-3,KK-3
        DO K=KSTART(6),KSTOP(6)
         DO J=JSTART(6),JSTOP(6)
          DO I=ISTART(6),ISTOP(6)
          AFLUX = AFLUX + DDX(I)*DDY(J)
          VFLUX = VFLUX + W(K,J,I)*DDX(I)*DDY(J)*DT
          ENDDO
         ENDDO
        ENDDO
        NPART_NEW_BOT = INT(AWB(IBOCOND_PP(6),6)*VFLUX)
      ENDIF
      NPART_NEW = NPART_NEW_FRO + NPART_NEW_BAC + 
     $            NPART_NEW_RGT + NPART_NEW_LFT +
     $            NPART_NEW_BOT + NPART_NEW_TOP

      ALLOCATE(NEW_IINDEX(NPART_NEW))
      ALLOCATE(NEW_JINDEX(NPART_NEW))
      ALLOCATE(NEW_KINDEX(NPART_NEW))
      ALLOCATE(X_NEW_PART(3,NPART_NEW))
      ALLOCATE(U_NEW_PART(3,NPART_NEW))
      COUNTER = 0
      DO NP=1,NPART_NEW_FRO
       COUNTER = COUNTER + 1
       NEW_PARTS(COUNTER) = PP_NEWINDEX()
       NEW_DIAMCLASS(COUNTER) = FLOOR(ranfmglet()*NDIAM_CLASSES)+1 
C       WRITE(6,*)'NEW INDEX: ',NP,NEW_PARTS_FRO(NP),NEW_DIAMCLASS(NP)
C       X_NEW_PART(1,COUNTER) = X(2) + DX(2)
       X_NEW_PART(1,COUNTER) = X(ISTART(1))
       
       DO
           CALL RANDOM_NUMBER(rrand)
           idx1 = INT(rrand(1)*(JSTOP(1) - JSTART(1)) + JSTART(1))
           idx2 = INT(rrand(2)*(KSTOP(1) - KSTART(1)) + KSTART(1))
           IF (bp(idx2, idx1, istart(1)) > 0.5) THEN
               EXIT
           END IF
       END DO
       CALL RANDOM_NUMBER(rrand)
       X_NEW_PART(2,COUNTER) = Y(idx1) + (rrand(1) - 0.5)*DY(idx1)
       X_NEW_PART(3,COUNTER) = Z(idx2) + (rrand(2) - 0.5)*DZ(idx2)
       
!       X_NEW_PART(2,COUNTER) = (Y(JSTART(1))+ 0.5*DY(JSTART(1)-1))+
!     $   ranfmglet()*(Y(JSTOP(1)) + 0.5*DY(JSTOP(1))
!     $          -Y(JSTART(1))+ 0.5*DY(JSTART(1)-1))
!       X_NEW_PART(3,COUNTER) = (Z(KSTART(1))+ 0.5*DZ(KSTART(1)-1))+
!     $   ranfmglet()*(Z(KSTOP(1)) + 0.5*DZ(KSTOP(1))
!     $          -Z(KSTART(1))+ 0.5*DZ(KSTART(1)-1))
C       WRITE(6,*)'NEW POS1: ',NP,X_NEW_PART(1,NP),X_NEW_PART(2,NP),
C     $                          X_NEW_PART(3,NP)
      ENDDO
      DO NP=1,NPART_NEW_BAC
       COUNTER = COUNTER + 1
       NEW_PARTS(COUNTER) = PP_NEWINDEX()
       NEW_DIAMCLASS(COUNTER) = FLOOR(ranfmglet()*NDIAM_CLASSES)+1
C       WRITE(6,*)'NEW INDEX: ',NP,NEW_PARTS_FRO(NP),NEW_DIAMCLASS(NP)
C       X_NEW_PART(1,COUNTER) = X(II-3) 
       X_NEW_PART(1,COUNTER) = X(ISTOP(2)) 
       X_NEW_PART(2,COUNTER) = (Y(JSTART(2))+ 0.5*DY(JSTART(2)-1))+
     $   ranfmglet()*(Y(JSTOP(2)) + 0.5*DY(JSTOP(2))
     $          -Y(JSTART(2))+ 0.5*DY(JSTART(2)-1))
       X_NEW_PART(3,COUNTER) = (Z(KSTART(2))+ 0.5*DZ(KSTART(2)-1))+
     $   ranfmglet()*(Z(KSTOP(2)) + 0.5*DZ(KSTOP(2))
     $          -Z(KSTART(2))+ 0.5*DZ(KSTART(2)-1))
C       WRITE(6,*)'NEW POS1: ',NP,X_NEW_PART(1,NP),X_NEW_PART(2,NP),
C     $                          X_NEW_PART(3,NP)
      ENDDO
      DO NP=1,NPART_NEW_RGT
       COUNTER = COUNTER + 1
       NEW_PARTS(COUNTER) = PP_NEWINDEX()
       NEW_DIAMCLASS(COUNTER) = FLOOR(ranfmglet()*NDIAM_CLASSES)+1
       X_NEW_PART(1,COUNTER) = (X(ISTART(3))+ 0.5*DX(ISTART(3)-1))+
     $   ranfmglet()*(X(ISTOP(3)) + 0.5*DX(ISTOP(3))
     $          -X(ISTART(3))+ 0.5*DX(ISTART(3)-1))
       X_NEW_PART(2,COUNTER) = Y(JSTART(3))
C       X_NEW_PART(2,COUNTER) = Y(3)
       X_NEW_PART(3,COUNTER) = (Z(KSTART(3))+ 0.5*DZ(KSTART(3)-1))+
     $   ranfmglet()*(Z(KSTOP(3)) + 0.5*DZ(KSTOP(3))
     $          -Z(KSTART(3))+ 0.5*DZ(KSTART(3)-1))
      ENDDO
      DO NP=1,NPART_NEW_LFT
       COUNTER = COUNTER + 1
       NEW_PARTS(COUNTER) = PP_NEWINDEX()
       NEW_DIAMCLASS(COUNTER) = FLOOR(ranfmglet()*NDIAM_CLASSES)+1
       X_NEW_PART(1,COUNTER) = (X(ISTART(4))+ 0.5*DX(ISTART(4)-1))+
     $   ranfmglet()*(X(ISTOP(4)) + 0.5*DX(ISTOP(4))
     $          -X(ISTART(4))+ 0.5*DX(ISTART(4)-1))
C       X_NEW_PART(2,COUNTER) = Y(JJ-3)
       X_NEW_PART(2,COUNTER) = Y(JSTART(4))
       X_NEW_PART(3,COUNTER) = (Z(KSTART(4))+ 0.5*DZ(KSTART(4)-1))+
     $   ranfmglet()*(Z(KSTOP(4)) + 0.5*DZ(KSTOP(4))
     $          -Z(KSTART(4))+ 0.5*DZ(KSTART(4)-1))
      ENDDO
      DO NP=1,NPART_NEW_BOT
       COUNTER = COUNTER + 1
       NEW_PARTS(COUNTER) = PP_NEWINDEX()
       NEW_DIAMCLASS(COUNTER) = FLOOR(ranfmglet()*NDIAM_CLASSES)+1
       X_NEW_PART(1,COUNTER) = (X(ISTART(5))+ 0.5*DX(ISTART(5)-1))+
     $   ranfmglet()*(X(ISTOP(5)) + 0.5*DX(ISTOP(5))
     $          -X(ISTART(5))+ 0.5*DX(ISTART(5)-1))
       X_NEW_PART(2,COUNTER) = (Y(JSTART(5))+ 0.5*DY(JSTART(5)-1))+
     $   ranfmglet()*(Y(JSTOP(5)) + 0.5*DY(JSTOP(5))
     $          -Y(JSTART(5))+ 0.5*DY(JSTART(5)-1))
       X_NEW_PART(3,COUNTER) = Z(KSTART(5))
C       X_NEW_PART(3,COUNTER) = Z(3)
      ENDDO
      DO NP=1,NPART_NEW_TOP
       COUNTER = COUNTER + 1
       NEW_PARTS(COUNTER) = PP_NEWINDEX()
       NEW_DIAMCLASS(COUNTER) = FLOOR(ranfmglet()*NDIAM_CLASSES)+1
       X_NEW_PART(1,COUNTER) = (X(ISTART(6))+ 0.5*DX(ISTART(6)-1))+
     $   ranfmglet()*(X(ISTOP(6)) + 0.5*DX(ISTOP(6))
     $          -X(ISTART(6))+ 0.5*DX(ISTART(6)-1))
       X_NEW_PART(2,COUNTER) = (Y(JSTART(6))+ 0.5*DY(JSTART(6)-1))+
     $   ranfmglet()*(Y(JSTOP(6)) + 0.5*DY(JSTOP(6))
     $          -Y(JSTART(6))+ 0.5*DY(JSTART(6)-1))
C       X_NEW_PART(3,COUNTER) = Z(KK-3)
       X_NEW_PART(3,COUNTER) = Z(KSTART(6))
      ENDDO
      COUNTERD1 = 0
      COUNTERD2 = 0
      COUNTERD3 = 0
      DO NP=1,NPART_NEW
       IF(NEW_DIAMCLASS(NP).EQ.1) THEN
        COUNTERD1 = COUNTERD1 + 1
       ELSEIF(NEW_DIAMCLASS(NP).EQ.2) THEN
        COUNTERD2 = COUNTERD2 + 1
       ELSEIF(NEW_DIAMCLASS(NP).EQ.3) THEN
        COUNTERD3 = COUNTERD3 + 1
       ENDIF
      ENDDO
      IF(NPART_NEW.GE.1) THEN
      WRITE(6,*)'NEW PARTS D1:',COUNTERD1
      WRITE(6,*)'NEW PARTS D2:',COUNTERD2
      WRITE(6,*)'NEW PARTS D3:',COUNTERD3
      ENDIF

      CALL SEARCHINDEX(KK,Z,DZ,
     $           NZGRAE,NPART_NEW,3,
     $           X_NEW_PART,
     $           NEW_KINDEX,HILF1,0)

      CALL SEARCHINDEX(JJ,Y,DY,
     $           NYGRAE,NPART_NEW,2,
     $           X_NEW_PART,
     $           NEW_JINDEX,HILF1,0)

      CALL SEARCHINDEX(II,X,DX,
     $           NXGRAE,NPART_NEW,1,
     $           X_NEW_PART,
     $           NEW_IINDEX,HILF1,0)

!      DO NP=1,NPART_NEW_FRO
!       WRITE(6,*)'NEW POS2: ',NP,NEW_IINDEX(NP),NEW_JINDEX(NP),
!     $                          NEW_KINDEX(NP)
!      ENDDO
C      DO NP=1,NPART_NEW_FRO
C       WRITE(6,*)'NEW POS2: ',NP,X_NEW_PART(1,NP),X_NEW_PART(2,NP),
C     $                          X_NEW_PART(3,NP)
C      ENDDO

      CALL PP_GETVELOCITIES_CL_CS (KK,JJ,II,
     $   X,Y,Z,
     $   DX,  DY,  DZ,
     $   DDX,  DDY,  DDZ,
     $   RDX, RDY, RDZ,
     $   RDDX, RDDY, RDDZ,
     $   U,V,W,HILF1, HILF2, HILF3, HILF4,
     $   NPART_NEW,
     $   X_NEW_PART(1,:),
     $   X_NEW_PART(2,:),
     $   X_NEW_PART(3,:),
     $   NEW_IINDEX,
     $   NEW_JINDEX,
     $   NEW_KINDEX,
     $   U_NEW_PART)

C      DO NP=1,NPART_NEW_FRO
C       WRITE(6,*)'NEW INDEX: ',NP,NEW_IINDEX(NP),NEW_JINDEX(NP),
C     $                          NEW_KINDEX(NP),IGRID
C       WRITE(6,*)'NEW VEL: ',NP,U_NEW_PART(1,NP),U_NEW_PART(2,NP),
C     $                          U_NEW_PART(3,NP),IGRID
C      ENDDO

      DO NP=1,NPART_NEW
       NPART_GRID = NPART_GRID + 1
       XPART(1,NPART_GRID) = X_NEW_PART(1,NP)
       XPART(2,NPART_GRID) = X_NEW_PART(2,NP)
       XPART(3,NPART_GRID) = X_NEW_PART(3,NP)
       UPART(1,NPART_GRID) = U_NEW_PART(1,NP)
       UPART(2,NPART_GRID) = U_NEW_PART(2,NP)
       UPART(3,NPART_GRID) = U_NEW_PART(3,NP)
       PART_IDX(NPART_GRID) = NEW_PARTS(NP)
       DIAMCLASS(NPART_GRID) = NEW_DIAMCLASS(NP)
       START_VEC(1,NPART_GRID) = X_NEW_PART(1,NP)
       START_VEC(2,NPART_GRID) = X_NEW_PART(2,NP)
       START_VEC(3,NPART_GRID) = X_NEW_PART(3,NP)
       UPART_F(1,NPART_GRID) = U_NEW_PART(1,NP)
       UPART_F(2,NPART_GRID) = U_NEW_PART(2,NP)
       UPART_F(3,NPART_GRID) = U_NEW_PART(3,NP)
       OLD_UPART_F(1,NPART_GRID) = U_NEW_PART(1,NP)
       OLD_UPART_F(2,NPART_GRID) = U_NEW_PART(2,NP)
       OLD_UPART_F(3,NPART_GRID) = U_NEW_PART(3,NP)
       DS(NPART_GRID) = 0.0
      ENDDO

!      DEALLOCATE(NEW_PARTS_TOP)
!      DEALLOCATE(NEW_PARTS_BOT)
!      DEALLOCATE(NEW_PARTS_RGT)
!      DEALLOCATE(NEW_PARTS_LFT)
      DEALLOCATE(NEW_PARTS)
!      DEALLOCATE(NEW_PARTS_BAC)
      DEALLOCATE(NEW_DIAMCLASS)
      DEALLOCATE(NEW_IINDEX)
      DEALLOCATE(NEW_JINDEX)
      DEALLOCATE(NEW_KINDEX)
      DEALLOCATE(X_NEW_PART)
      DEALLOCATE(U_NEW_PART)

      RETURN
      END
