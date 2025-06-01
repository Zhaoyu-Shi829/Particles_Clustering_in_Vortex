      SUBROUTINE PP_READ(CHANNEL, IDIM2D, IDIM3D)
C--------------------------------------------------------
C
C     READS PARTIClE-FILE FOR RESTART
C
C     28.05.03 (FS): ORIGINAL
C     27.04.05 (CG): ERWEITERT UM PART_IDX(I),START_GRID(I), SCALAR_VAL(I)
C
C--------------------------------------------------------

      USE pp_settings_mod, ONLY: pp_mpi_integer
      USE pointer_mod, ONLY: get_ip1
      USE fields_mod, ONLY: x, y, z, dx, dy, dz

      IMPLICIT NONE
#include "../colevel.h"
#include "../mgpar.h"
#include "pp_particles.h"

#include "mpif.h"
#include "../compi.h"



      INTEGER CHANNEL
      INTEGER IDIM3D     
      INTEGER IDIM2D



C     LOCAL VARS

      INTEGER I,IGRID_FL,IGRID_UNSL,NGRID, IDTAG, VERSION,ILEVEL
            
      LOGICAL LSORT


      INTEGER J,NGRID_OLD,
     $      GRID_OFFSET_OLD(MAXGRIDS)
      INTEGER  NPART_REST_N(MAXGRIDS)

      INTEGER PP_NUM_VARS

      INTEGER NSECTIONS, ISECTION
      CHARACTER (LEN=32) SECTIONTYPE


      INTEGER PROC_READ
      INTEGER IERR

      INTEGER NUMVARS_SOLL
      INTEGER NTOTGRDS,IGRID_NTOTGRDS,IGRID
      ALLOCATABLE IGRID_NTOTGRDS(:)

      INTEGER :: II, JJ, KK, IP1
      REAL :: XMAX_HERE, YMAX_HERE, ZMAX_HERE
      REAL, ALLOCATABLE :: XMIN_HERE(:), YMIN_HERE(:), ZMIN_HERE(:)
      REAL, ALLOCATABLE :: XTOT_HERE(:), YTOT_HERE(:), ZTOT_HERE(:)

      ALLOCATE(IGRID_NTOTGRDS(MAXGRIDS))
      NTOTGRDS = 0
      IGRID_NTOTGRDS = 0
      DO ILEVEL=MINLEVEL,MAXLEVEL
       NGRID = NOFTST(ILEVEL)
       DO IGRID=1,NGRID
       NTOTGRDS = NTOTGRDS + 1
       IGRID_NTOTGRDS(NTOTGRDS) = IGRDOFTST(IGRID,ILEVEL)
       ENDDO
      ENDDO



C      ILEVEL=MAXLEVEL

      IF (.NOT. (PP_LREAD)) RETURN




      

      NUMVARS_SOLL=18

#if defined _PP_SHOT_SGS_ || defined _PP_SIM_SGS_
      NUMVARS_SOLL=NUMVARS_SOLL+3
#endif

#ifdef  _PP_WRITE_UPART_EFF_
      NUMVARS_SOLL=NUMVARS_SOLL+3
#endif

#ifdef _PP_FILTER_APRIORI_
      NUMVARS_SOLL=NUMVARS_SOLL+3
#endif

#ifdef _PP_VEL_FILTER_POS_UNFILTER_
      NUMVARS_SOLL=NUMVARS_SOLL+6
#endif
      
      IF (PP_NUM_VARS() .NE. NUMVARS_SOLL) 
     $   CALL ERRR(508,' PP_READ ')



C      NGRID=NOFTST(MAXLEVEL)



C     nur der prozessor, der das letzte gitter hat, liest
C      PROC_READ = IDPROCOFGRD(IGRDOFTST(NGRID,MAXLEVEL))
      PROC_READ = IDPROCOFGRD(IGRID_NTOTGRDS(NTOTGRDS))
!      PROC_READ = 0
      IF (MYID .EQ. PROC_READ) THEN



C     READ HEADER
      READ(CHANNEL,ERR=2000) IDTAG
      
      IF (IDTAG .NE. 987519784) THEN
         VERSION = 0
         NPART=IDTAG
      ELSE
         READ(CHANNEL,ERR=2000) VERSION
         READ(CHANNEL,ERR=2000) NPART
      ENDIF
      
      READ(CHANNEL,ERR=2000) RAND_WALK_MAX

      READ(CHANNEL,ERR=2000) NGRID_OLD


      ENDIF

      CALL MPI_BCAST(
     $   VERSION,
     $   1,
     $   pp_mpi_integer,
     $   PROC_READ,
     $   MPI_WHICH_COMM,
     $   IERR)

      CALL MPI_BCAST(
     $   NGRID_OLD,
     $   1,
     $   pp_mpi_integer,
     $   PROC_READ,
     $   MPI_WHICH_COMM,
     $   IERR)


      LSORT=.FALSE.
      IF ((NTOTGRDS .NE. NGRID_OLD) .AND. (VERSION .LE. 0)) THEN
            WRITE(6,*) 'IN VERSION 0 NGRID MUST EQUAL NGRID_OLD'
            CALL ERRR (502,' PP_READ ')
      ENDIF
      IF ((NTOTGRDS .LT. NGRID_OLD) .AND. (VERSION .GT. 0)) THEN
            WRITE(6,*) 'NGRID MUST GREATER OR EQUAL NGRID_OLD'
            CALL ERRR (502,' PP_READ ')
      ENDIF
      IF ((NTOTGRDS .GT. NGRID_OLD) .AND. (VERSION .GT. 0)) THEN
C        OTHER GRIDS - MAYBE DIFFERENT SLICING? SORT PARTICLES IN GRIDS!
         LSORT=.TRUE.
      ENDIF


      CALL MPI_BCAST(
     $   NPART,
     $   1,
     $   pp_mpi_integer,
     $   PROC_READ,
     $   MPI_WHICH_COMM,
     $   IERR)


      CALL MPI_BCAST(
     $   RAND_WALK_MAX,
     $   1,
     $   MPI_DOUBLE_OR_SINGLE,
     $   PROC_READ,
     $   MPI_WHICH_COMM,
     $   IERR)



      ALLOCATE(XMIN_HERE(MAXGRIDS))
      ALLOCATE(YMIN_HERE(MAXGRIDS))
      ALLOCATE(ZMIN_HERE(MAXGRIDS))
      ALLOCATE(XTOT_HERE(MAXGRIDS))
      ALLOCATE(YTOT_HERE(MAXGRIDS))
      ALLOCATE(ZTOT_HERE(MAXGRIDS))
      DO ILEVEL=MINLEVEL,MAXLEVEL
          DO I = 1, NOFTST(ILEVEL)
              IGRID = IGRDOFTST(I,ILEVEL)
              CALL get_ip1(ip1, igrid)
              CALL mgdims(kk, jj, ii, igrid)
              
              XMAX_HERE = x(ip1+ii-3) + dx(ip1+ii-3)/2.0
              XMIN_HERE(igrid) = x(ip1+1) + dx(ip1+1)/2.0
              XTOT_HERE(IGRID) = XMAX_HERE - XMIN_HERE(igrid)

              YMAX_HERE = y(ip1+jj-3) + dy(ip1+jj-3)/2.0
              YMIN_HERE(igrid) = y(ip1+1) + dy(ip1+1)/2.0
              YTOT_HERE(IGRID) = YMAX_HERE - YMIN_HERE(igrid)

              ZMAX_HERE = z(ip1+kk-3) + dz(ip1+kk-3)/2.0
              ZMIN_HERE(igrid) = z(ip1+1) + dz(ip1+1)/2.0
              ZTOT_HERE(IGRID) = ZMAX_HERE - ZMIN_HERE(igrid)
          END DO
      END DO
      CALL PP_DIM_ARRAYS(
     $   XTOT_HERE,YTOT_HERE,ZTOT_HERE
#if defined _IB_
     $   ,IDIM3D
#endif      
#if defined _PP_INIT_FROM_BED_
     $   ,IDIM2D
#endif      
     $   )

c        folgendes statement evtl. noetig falls NGRID .GT. NGRID_OLD
      NPART_GRID=0

         IF (MYID .EQ. PROC_READ) THEN
         READ(CHANNEL,ERR=2000) NSECTIONS
         ENDIF

         CALL MPI_BCAST(
     $      NSECTIONS,
     $      1,
     $      pp_mpi_integer,
     $      PROC_READ,
     $      MPI_WHICH_COMM,
     $      IERR)

         
         DO ISECTION=1,NSECTIONS

            IF (MYID .EQ. PROC_READ) THEN

            READ(CHANNEL,ERR=2000) SECTIONTYPE

            ENDIF
   
            CALL MPI_BCAST(
     $         SECTIONTYPE,
     $         32,
     $         MPI_CHARACTER,
     $         PROC_READ,
     $         MPI_WHICH_COMM,
     $         IERR)
         
            IF (SECTIONTYPE .EQ. 'PARTICLE_ORDERED                ') 
     $         THEN
               CALL PP_READ_SECTION(CHANNEL, PROC_READ, NGRID_OLD,
     $            ISECTION)
            ELSE
               WRITE(6,*) 'UNKNOWN SECTIONTYPE', SECTIONTYPE
               CALL ERRR(512, ' PP_READ ')
            ENDIF
         ENDDO



      IF (LSORT) THEN


         CALL PP_SORT_IN_GRIDS(NPART_MAX,NPART_2D_MAX,NGRID,
     $   MAXBOCONDS,MAXGRIDS,GRID_OFFSET,NPART_GRID,
     $   XPART,
     $   PP_SHIFT_X,PP_SHIFT_Y,PP_SHIFT_Z,
     $   IGRDOFTST(1:NGRID,MAXLEVEL),
     $   MYID,IDPROCOFGRD,
     $   XMIN_HERE,XTOT_HERE,YMIN_HERE,
     &   YTOT_HERE,ZMIN_HERE,ZTOT_HERE
     $   )
      ENDIF




C     FIND RESTING PARTICLES
      NBR_IN_REST = -1
      DO IGRID_UNSL=1,NTOTGRDS
         IGRID_FL = IGRID_NTOTGRDS(IGRID_UNSL)
C         IGRID_FL = IGRDOFTST(IGRID_UNSL,ILEVEL)
         IF (MYID .EQ. IDPROCOFGRD(IGRID_FL)) THEN

         CALL PP_FIND_RESTING_PARTICLES(NPART_GRID(IGRID_FL), 
     $      NPART_REST(IGRID_FL), PART_IDX(GRID_OFFSET(IGRID_FL)+1), 
     $      IDX_REST(GRID_OFFSET2D(IGRID_FL)+1), 
     $      NBR_REST(GRID_OFFSET2D(IGRID_FL)+1), 
     $      NBR_IN_REST(GRID_OFFSET(IGRID_FL)+1))
      
         ENDIF
      ENDDO

      DEALLOCATE(XMIN_HERE)
      DEALLOCATE(YMIN_HERE)
      DEALLOCATE(ZMIN_HERE)
      DEALLOCATE(XTOT_HERE)
      DEALLOCATE(YTOT_HERE)
      DEALLOCATE(ZTOT_HERE)

      DEALLOCATE(IGRID_NTOTGRDS)

      RETURN

 2000 CALL ERRR (501,'PP_READ')

      END

