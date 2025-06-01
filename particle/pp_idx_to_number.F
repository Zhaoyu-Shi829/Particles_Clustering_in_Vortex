      SUBROUTINE PP_IDX_TO_NUMBER(NPART, PART_IDX, NPART_IDX, IDX, NBR)
      IMPLICIT NONE
      
      INTEGER NPART, NPART_IDX, IDX(NPART_IDX), NBR(NPART_IDX), 
     $   PART_IDX(NPART)
      
C      LOCAL VARS
      
      INTEGER IPART1, IPART2, IPART3
      
      
c      DEBUG
c      NBR=0
      
      DO IPART1=1,NPART
         DO IPART2=1,NPART_IDX
            IF (PART_IDX(IPART1) .EQ. IDX(IPART2)) THEN
               NBR(IPART2)=IPART1
            ENDIF
         ENDDO
      ENDDO
      
      
      
C     DEBUG
c      DO IPART2=1,NPART_IDX
c      IF (NBR(IPART2) .EQ. 0) THEN
c         WRITE (6,*) 'ERROR IN PP_IDX_TO_NUMBER'
c         CALL ERRR(501,'PP_IDX_TO_NUMBER')
c      ENDIF
c      ENDDO
      
      RETURN 
      END
      
