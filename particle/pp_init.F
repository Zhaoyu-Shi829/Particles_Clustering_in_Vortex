       SUBROUTINE PP_INIT (IDIM3D,IDIM2D,IDIM1D,idim1dl,
     $      U,V,W,
     $      X, Y, Z, DX,DY,DZ,
     $      DDX,DDY,DDZ,
     $      RDX, RDY, RDZ,
     $      RDDX,RDDY,RDDZ,
     $      NXGRAE,NYGRAE,NZGRAE,
     $      XMIN,YMIN,ZMIN
#ifdef _IB_
     $    ,BP
#endif
     $     ,RHO,GMOL
     $     )
C--------------------------------------------------
C
C     INITIALISES PARTICLE POSITIONS:
C
C         PP_FLOWCASE == 1 > T-MISCHER
C         PP_FLOWCASE == 2 > CHANNEL
C         PP_FLOWCASE == 3 > T-Mischer + FRED BODY
C         PP_FLOWCASE == 4 > Benchmark T-Mischer
C         PP_FLOWCASE == 5 > Gleichverteilt �ber gesamte Geometrie
C         PP_FLOWCASE == 6 > CHANNEL Emscher, TEILKREIS (D.H. + FRED BODY)
C         PP_FLOWCASE == 7 > CHANNEL RECHTECKIG, f�r armenio, PARTIKEL WERDEN BEI Z=CONST INITIALISIERT
C         PP_FLOWCASE == 8 > CHANNEL RECHTECKIG, f�r rouse, PARTIKEL WERDEN MIT EXP KONZENTRATIONSVERTEILUNG INITIALISIERT
C         PP_FLOWCASE == 9 > Kolk (Geometrie Oscar), Partikel werden nur im Bereich
C                             3 <= x <= 5.29    X    0.68 <= y <= 1.32    X  0 < z < 0.3 initialisiert
C         PP_FLOWCASE == 10 > particles are initialized in a box specified by pp_init_xmin, pp_init_xmax, pp_init_ymin,
C                             pp_init_ymax, pp_init_zmin, pp_init_zmax, +fred_body moeglich
C
C     26.05.03 (FS): ORIGINAL
C     14.04.05 (CG): Benchmark T-Mischer
C
C     21.04.05 (CG)    GRID_OFFSET(I):    OFFSET ZU GITTER I IN XPART / PART_IDX - GRID_OFFSET(1)=0
C                                NPART_GRID(I):          ANZAHL DER (BELEGTEN PL�TZE F�R) PARTIKEL IM GITTER I
C                                PART_IDX:              ARRAY DERSELBEN STRUKTUR WIE XPART, WELCHES DEN (GLOBALEN) PARTIKELINDEX ENTH�LT
C
C--------------------------------------------------
      USE pp_fields_mod
      USE pp_settings_mod, ONLY: pp_set_precision
      USE pointer_mod, ONLY: get_ip1, get_ip1l
      USE pp_statistics_mod, ONLY: pp_statistics_init

      IMPLICIT NONE

#include "../colevel.h"
#include "../mgpar.h"
#include "pp_particles.h"

C
C               MPI-INCLUDE
C
#include "../compi.h"



      REAL X(IDIM1D),Y(IDIM1D),Z(IDIM1D),
     $   U(IDIM3D),V(IDIM3D),W(IDIM3D)
      REAL DX(IDIM1D), DY(IDIM1D), DZ(IDIM1D),
     $   DDX(IDIM1D),DDY(IDIM1D),DDZ(IDIM1D)
      REAL RDX(IDIM1DL), RDY(IDIM1DL), RDZ(IDIM1DL)
      REAL     ranfmglet
      REAL RDDX(IDIM1DL), RDDY(IDIM1DL), RDDZ(IDIM1DL)
      INTEGER IDIM3D, IDIM2D, IDIM1D, idim1dl
      INTEGER NXGRAE(MAXGRIDS),NYGRAE(MAXGRIDS),NZGRAE(MAXGRIDS)
      REAL XMIN(MAXGRIDS), YMIN(MAXGRIDS), ZMIN(MAXGRIDS)
#ifdef _IB_
      REAL BP(IDIM3D)
#endif
      REAL RHO,GMOL







      REAL XMAX_HERE, YMAX_HERE, ZMAX_HERE, SINKVEL, PP_PART_DIAMETER
      REAL XMIN_HERE, YMIN_HERE, ZMIN_HERE
      REAL XTOT_HERE, YTOT_HERE, ZTOT_HERE
      ALLOCATABLE XTOT_HERE(:), YTOT_HERE(:), ZTOT_HERE(:)

      INTEGER NPART_GRID_SAVE
      ALLOCATABLE NPART_GRID_SAVE(:)

#if defined _PP_SHARP_CUTOFF_PARTVEL_
      REAL PP_FILTERWIDTH
#endif


      INTEGER STARTINDEX, STOPINDEX, ILEVEL, IGRID_UNSL, IGRID_FL, IPART
      INTEGER II,JJ,KK
      INTEGER IP3,IP2,IP1,IBB,IBU,
     $        NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
     $        IP3T,IP2T,IP1T,IBBT
      INTEGER :: ip1l
      INTEGER IINDEX, JINDEX, KINDEX

      ALLOCATABLE IINDEX(:),JINDEX(:),KINDEX(:)

      INTEGER PP_NEWINDEX

      CALL pp_set_precision()

      CALL ppFieldsAlloc(idim3d, idim2d)
      IF (pp_stat) THEN
          CALL pp_statistics_init(idim3d)
      END IF

#ifdef _COMPUTE_MATERIAL_DERIV_
      CALL ppMatderFieldsAlloc(idim3d)
#endif

      IF (PP_LREAD) THEN
c#if defined _PP_WRITE_BED_ONLY_ && defined _PP_DIAM_FILE_
C        SET DIAMETERS OF ALL PARTICLES SUCH THAT DIAMETER OF PARTICLES
C        WHICH WILL BE INITILAIZED LATER IS SET
c         CALL PP_DIAM_FROM_FILE(PP_NXNPART,
c     $         DIAMCLASS, NDIAM_CLASSES)
c#endif
         RETURN
      ENDIF



      DS = 0

      MEAN_SUSPART_DS=0
      NSUSPART_DS=0


#ifdef _PP_SCHMEECKLE_BOUND_
      NCONTACTS=0
#endif
      NPART_REST=0
      NBR_IN_REST=-1








      IF ((PP_FLOWCASE .EQ. 5) .OR. (PP_FLOWCASE .EQ. 6) .OR.
     $   (PP_FLOWCASE .EQ. 7) .OR. (PP_FLOWCASE .EQ. 8)
     $    .OR. (PP_FLOWCASE .EQ. 9) .OR. (PP_FLOWCASE .EQ. 10)) THEN

       NPART=PP_NXNPART

C       ILEVEL = MAXLEVEL


      IF ((PP_FLOWCASE .EQ. 9) .OR. (PP_FLOWCASE .EQ. 10)) THEN
C        COMPUTE NPART_GRID
         ALLOCATE(XTOT_HERE(MAXGRIDS))
         ALLOCATE(YTOT_HERE(MAXGRIDS))
         ALLOCATE(ZTOT_HERE(MAXGRIDS))
         IF (PP_FLOWCASE .EQ. 9) THEN
C                             3 <= x <= 5.29    X    0.68 <= y <= 1.32    X   0 < z < 0.3
            PP_INIT_XMIN=3.0
            PP_INIT_XMAX=5.29
            PP_INIT_YMIN=0.68
            PP_INIT_YMAX= 1.32
            PP_INIT_ZMIN=0.0
            PP_INIT_ZMAX=0.3
         ENDIF
         DO ILEVEL=MINLEVEL,MAXLEVEL
         DO IGRID_UNSL = 1,NOFTST(ILEVEL)
            IGRID_FL = IGRDOFTST(IGRID_UNSL,ILEVEL)

            CALL get_ip1(ip1, igrid_fl)
            CALL mgdims(kk, jj, ii, igrid_fl)

            XMAX_HERE=MIN(PP_INIT_XMAX, x(ip1+ii-3) + dx(ip1+ii-3)/2.0)
            XMIN_HERE=MAX(PP_INIT_XMIN, x(ip1+1) + dx(ip1+1)/2.0)
            XTOT_HERE(IGRID_FL)=MAX(0.0, XMAX_HERE - XMIN_HERE)

            YMAX_HERE=MIN(PP_INIT_YMAX, y(ip1+jj-3) + dy(ip1+jj-3)/2.0)
            YMIN_HERE=MAX(PP_INIT_YMIN, y(ip1+1) + dy(ip1+1)/2.0)
            YTOT_HERE(IGRID_FL)=MAX(0.0, YMAX_HERE - YMIN_HERE)

            ZMAX_HERE=MIN(PP_INIT_ZMAX, z(ip1+kk-3) + dz(ip1+kk-3)/2.0)
            ZMIN_HERE=MAX(PP_INIT_ZMIN, z(ip1+1) + dz(ip1+1)/2.0)
            ZTOT_HERE(IGRID_FL)=MAX(0.0,ZMAX_HERE - ZMIN_HERE)

C            ALLOW SOME SPACE FOR ROUNDOFF ERRORS
C            ZTOT_HERE(IGRID_FL)=MAX(ZTOT_HERE(IGRID_FL),0.05)
         ENDDO
         ENDDO

         CALL PP_DIM_ARRAYS(
     $      XTOT_HERE,YTOT_HERE,ZTOT_HERE
#if defined _IB_
     $      ,IDIM3D
#endif
#if defined _PP_INIT_FROM_BED_
     $      ,IDIM2D
#endif
     $   )

C         WRITE(6,*)'NACH ERSTEM AUFRUF VON PP_DIM_ARRAYS'

#ifdef _IB_
C         WRITE(6,*)'VOR DEALLOCATE IN PP_INIT',MYID
         DEALLOCATE(PP_BODY_NVEC)
         DEALLOCATE(PP_BODY_NVECNO)
#endif

#ifdef _PP_INIT_FROM_BED_
         DEALLOCATE(PP_BEDSUM)
#endif

         ALLOCATE(NPART_GRID_SAVE(MAXGRIDS))
         NPART_GRID_SAVE=NPART_GRID


         DEALLOCATE(XTOT_HERE)
         DEALLOCATE(YTOT_HERE)
         DEALLOCATE(ZTOT_HERE)
      ENDIF



      ALLOCATE(XTOT_HERE(MAXGRIDS))
      ALLOCATE(YTOT_HERE(MAXGRIDS))
      ALLOCATE(ZTOT_HERE(MAXGRIDS))
      DO ILEVEL=MINLEVEL,MAXLEVEL
          DO IGRID_UNSL = 1,NOFTST(ILEVEL)
              IGRID_FL = IGRDOFTST(IGRID_UNSL,ILEVEL)

              CALL get_ip1(ip1, igrid_fl)
              CALL mgdims(kk, jj, ii, igrid_fl)

              XMAX_HERE = x(ip1+ii-3) + dx(ip1+ii-3)/2.0
              XMIN_HERE = x(ip1+1) + dx(ip1+1)/2.0
              XTOT_HERE(IGRID_FL) = XMAX_HERE - XMIN_HERE

              YMAX_HERE = y(ip1+jj-3) + dy(ip1+jj-3)/2.0
              YMIN_HERE = y(ip1+1) + dy(ip1+1)/2.0
              YTOT_HERE(IGRID_FL) = YMAX_HERE - YMIN_HERE

              ZMAX_HERE = z(ip1+kk-3) + dz(ip1+kk-3)/2.0
              ZMIN_HERE = z(ip1+1) + dz(ip1+1)/2.0
              ZTOT_HERE(IGRID_FL) = ZMAX_HERE - ZMIN_HERE
           ENDDO
       ENDDO

      CALL PP_DIM_ARRAYS(
     $   XTOT_HERE, YTOT_HERE, ZTOT_HERE
#if defined _IB_
     $      ,IDIM3D
#endif
#if defined _PP_INIT_FROM_BED_
     $      ,IDIM2D
#endif
     $   )
      DEALLOCATE(XTOT_HERE)
      DEALLOCATE(YTOT_HERE)
      DEALLOCATE(ZTOT_HERE)
C         WRITE(6,*)'NACH ZEITEM AUFRUF VON PP_DIM_ARRAYS'

      IF ((PP_FLOWCASE .EQ. 9) .OR. (PP_FLOWCASE .EQ. 10)) THEN
         NPART_GRID=NPART_GRID_SAVE
         DEALLOCATE(NPART_GRID_SAVE)
      ENDIF



C      WRITE(6,*)'ILEVEL',ILEVEL,MINLEVEL,MAXLEVEL
C      ILEVEL = MAXLEVEL
      ALLOCATE(IINDEX( NPART_PER_GRID_MAX ))
      ALLOCATE(JINDEX( NPART_PER_GRID_MAX ))
      ALLOCATE(KINDEX( NPART_PER_GRID_MAX ))
      DO ILEVEL=MINLEVEL,MAXLEVEL
      DO IGRID_UNSL = 1,NOFTST(ILEVEL)

      IGRID_FL = IGRDOFTST(IGRID_UNSL,ILEVEL)

      IF(MYID .EQ. IDPROCOFGRD(IGRID_FL)) THEN


      CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $        NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID_FL,
     $        IP3T,IP2T,IP1T,IBBT)
      CALL get_ip1l(ip1l, igrid_fl)

      STARTINDEX = GRID_OFFSET(IGRID_FL)+1
      STOPINDEX = GRID_OFFSET(IGRID_FL)+NPART_GRID(IGRID_FL)

      IF (STOPINDEX .GE. STARTINDEX) THEN
C     CHECK IF THERE IS A PARTICLE IN THIS GRID


      DO IPART=STARTINDEX,STOPINDEX
            PART_IDX(IPART)=PP_NEWINDEX()
#ifdef _PP_DIAM_FILE_
            CALL PP_DIAM_FROM_FILE(NPART_GRID(IGRID_FL),
     $         DIAMCLASS(STARTINDEX), NDIAM_CLASSES)
#else
            DIAMCLASS(IPART)=FLOOR(ranfmglet()*NDIAM_CLASSES)+1
#endif
      ENDDO


C      BELEGUNG PARTIKELPOSITIONEN
       XMIN_HERE = x(ip1+1) + dx(ip1+1)/2.0
       YMIN_HERE = y(ip1+1) + dy(ip1+1)/2.0
       ZMIN_HERE = z(ip1+1) + dz(ip1+1)/2.0
       XMAX_HERE = x(ip1+ii-3) + dx(ip1+ii-3)/2.0
       YMAX_HERE = y(ip1+jj-3) + dy(ip1+jj-3)/2.0
       ZMAX_HERE = z(ip1+kk-3) + dz(ip1+kk-3)/2.0
      IF (PP_FLOWCASE .EQ. 10) THEN
            XMIN_HERE=MAX(PP_INIT_XMIN,XMIN_HERE)
            XMIN_HERE=MIN(XMAX_HERE,XMIN_HERE)
            XMAX_HERE=MIN(PP_INIT_XMAX,XMAX_HERE)
            XMAX_HERE=MAX(XMIN_HERE,XMAX_HERE)

            YMIN_HERE=MAX(PP_INIT_YMIN,YMIN_HERE)
            YMIN_HERE=MIN(YMAX_HERE,YMIN_HERE)
            YMAX_HERE=MIN(PP_INIT_YMAX,YMAX_HERE)
            YMAX_HERE=MAX(YMIN_HERE,YMAX_HERE)

            ZMIN_HERE=MAX(PP_INIT_ZMIN,ZMIN_HERE)
            ZMIN_HERE=MIN(ZMAX_HERE,ZMIN_HERE)
            ZMAX_HERE=MIN(PP_INIT_ZMAX,ZMAX_HERE)
            ZMAX_HERE=MAX(ZMIN_HERE,ZMAX_HERE)
      ENDIF

      IF (PP_FLOWCASE .EQ. 7) THEN
         ZMIN_HERE = 0.025
         ZMAX_HERE = 0.025
         WRITE(6,*) 'PP_FLOWCASE 7: IGRID, XMIN_HERE, XMAX_HERE:',
     $      IGRID_FL,XMIN_HERE,XMAX_HERE
         WRITE(6,*) 'PP_FLOWCASE 7: IGRID, YMIN_HERE, YMAX_HERE:',
     $      IGRID_FL,YMIN_HERE,YMAX_HERE
         WRITE(6,*) 'PP_FLOWCASE 7: IGRID, ZMIN_HERE, ZMAX_HERE:',
     $      IGRID_FL,ZMIN_HERE,ZMAX_HERE
      ENDIF

      IF (PP_FLOWCASE .EQ. 8) THEN
         ZMIN_HERE = 0.0
         ZMAX_HERE = 1.0
      ENDIF

       CALL PP_RAND_PART_POS(NPART_GRID(IGRID_FL),
     $      XPART(:,STARTINDEX:STOPINDEX),
     $      XMIN_HERE,XMAX_HERE,YMIN_HERE,YMAX_HERE,ZMIN_HERE,ZMAX_HERE,
     $      IINDEX,
     $      JINDEX,
     $      KINDEX,KK,JJ,II,
     $      X(IP1),Y(IP1),Z(IP1),
     $      DX(IP1),DY(IP1),DZ(IP1), DDX(IP1),DDY(IP1),DDZ(IP1),
     $      NXGRAE(IGRID_FL),NYGRAE(IGRID_FL),NZGRAE(IGRID_FL), HILF1
#if defined _PP_SCHMEECKLE_BOUND_ || defined _PP_COLLIDE_
     $      ,DIAMCLASS(STARTINDEX:STOPINDEX)
#endif
#ifdef _IB_
     $      ,BP(IP3)
#endif
     $      )


!     WARNING:
!     This does obviously not work for general cases
!     because z(2) is hardcoded
      IF (PP_FLOWCASE .EQ. 8) THEN
C        EXPONENTIALVERTEILUNG �BER Z AUS GLEICHVERTEILUNG ERZEUGEN
         DO IPART=1,NPART
            SINKVEL = 1./18.*(RHOPART-RHO)*
     $         (PP_PART_DIAMETER(DIAMCLASS(IPART)))**2/GMOL*GRAVITATION
            XPART(3,IPART)=-DISPERS_FOR_INIT/SINKVEL
     $         *LOG(1-XPART(3,IPART))+
     $         Z(2)+0.5*DZ(2)
            XPART(3,IPART)=MAX(Z(2)+0.5*DZ(2),XPART(3,IPART))
            XPART(3,IPART)=MIN(Z(KK-2)+0.5*DZ(KK-2),XPART(3,IPART))
         ENDDO
      ENDIF

      START_VEC=XPART



      RAND_WALK_MAX = MIN(MIN(DX(1),DY(1)),DZ(1))/2.
#ifdef _PP_SHARP_CUTOFF_PARTVEL_
      RAND_WALK_MAX = PP_FILTERWIDTH(KK,JJ,II,3,3,3,DX(1),DY(1),DZ(1))
#endif




#ifdef _PP_INTERPOL_SL_CS_
      CALL PP_GETVELOCITIES_SL_CS(KK,JJ,II,
     $   X(IP1),Y(IP1),Z(IP1),
     $   DX(IP1),  DY(IP1),  DZ(IP1),
     $   RDDX(IP1l),RDDY(IP1l),RDDZ(IP1l),
     $   U(IP3),V(IP3),W(IP3),
     $   NPART_GRID(IGRID_FL),
     $   XPART(1,STARTINDEX:STOPINDEX),
     $   XPART(2,STARTINDEX:STOPINDEX),
     $   XPART(3,STARTINDEX:STOPINDEX),
     $   IINDEX,
     $   JINDEX,
     $   KINDEX,
     $   UPART(:,STARTINDEX:STOPINDEX)
     $   )
#elif defined _PP_INTERPOL_CL_CS_ || defined _PP_INTERPOL_4THORDER_
      CALL PP_GETVELOCITIES_CL_CS (KK,JJ,II,
     $   X(IP1),Y(IP1),Z(IP1),
     $   DX(IP1),  DY(IP1),  DZ(IP1),
     $   DDX(IP1),  DDY(IP1),  DDZ(IP1),
     $   RDX(IP1l), RDY(IP1l), RDZ(IP1l),
     $   RDDX(IP1l), RDDY(IP1l), RDDZ(IP1l),
     $   U(IP3),V(IP3),W(IP3),HILF1, HILF2, HILF3, HILF4,
     $   NPART_GRID(IGRID_FL),
     $   XPART(1,STARTINDEX:STOPINDEX),
     $   XPART(2,STARTINDEX:STOPINDEX),
     $   XPART(3,STARTINDEX:STOPINDEX),
     $   IINDEX,
     $   JINDEX,
     $   KINDEX,
     $   UPART(:,STARTINDEX:STOPINDEX)
     $   )
#else
            WRITE(6,*) 'GETVELOCITIES NOT IMPLEMENTED IN ',
     $         'PP_INIT'
            CALL ERRR(504, 'PP_INIT')
#endif

c      UPART(1,STARTINDEX:STOPINDEX)=UPART(1,STARTINDEX:STOPINDEX)
c     $   +(RHOPART-RHO)*DIAMPART_MIN**2/18./GMOL
c     $   *GRAVITATION*SIN(SLOPE/180.0*3.1415927)
c      UPART(3,STARTINDEX:STOPINDEX)=UPART(3,STARTINDEX:STOPINDEX)
c     $   -(RHOPART-RHO)*DIAMPART_MIN**2/18./GMOL
c     $   *GRAVITATION*COS(SLOPE/180.0*3.1415927)



c      WRITE(0,*) 'BEEP: SETTING FIXED PARTICLE VEL'
c      UPART(2,:)=1.0

      OLD_UPART_F(:,STARTINDEX:STOPINDEX)=UPART(:,STARTINDEX:STOPINDEX)
#if defined _PP_SHOT_SGS_ || defined _PP_SIM_SGS_
      UPART_F_EFF(:,STARTINDEX:STOPINDEX)=UPART(:,STARTINDEX:STOPINDEX)
#endif
#ifdef _PP_VEL_FILTER_POS_UNFILTER_
      UPART_UNFIL(:,STARTINDEX:STOPINDEX)=UPART(:,STARTINDEX:STOPINDEX)
      OLD_UPART_F_UNFIL(:,STARTINDEX:STOPINDEX)=
     $   UPART(:,STARTINDEX:STOPINDEX)
#endif

      ENDIF
c      ENDIF (STOPINDEX .GE. STARTINDEX) THEN


         ENDIF

      ENDDO
      ENDDO




      DEALLOCATE(IINDEX)
      DEALLOCATE(JINDEX)
      DEALLOCATE(KINDEX)









      ELSE
         WRITE (6,*) 'PP_FLOWCASE',PP_FLOWCASE
         CALL ERRR (501,' PP_INIT ')
      ENDIF




c#ifdef _PP_WRITE_BED_ONLY_
c      NPART_GRID=1
c#endif




      RETURN
 1000 FORMAT (A9,I3,3(2X,E12.4E2),I3,1X,A5,(E12.4E2))
      END
