MODULE pp_statistics_mod

    USE precision_mod, ONLY: real64

    IMPLICIT NONE
    PRIVATE

    REAL, ALLOCATABLE, TARGET ::     &
        npart_g(:), a_npart_g(:), &
        upart_g(:), a_upart_g(:), &
        vpart_g(:), a_vpart_g(:), &
        wpart_g(:), a_wpart_g(:)

    REAL(kind=real64) :: avg_time_tot = 0.0
    INTEGER :: nclasses

    PUBLIC :: pp_statistics, pp_statistics_init, pp_statistics_finish

CONTAINS

    SUBROUTINE pp_statistics(ii, jj, kk, igrid, nxgrae, nygrae, nzgrae, &
        npart_max, npart_grid, iindex, jindex, kindex, &
        part_idx, xpart, diamclass, upart, x, y, z, dx, dy, dz, hilf, dt)

        USE precision_mod, ONLY: realk
        USE pointer_mod, ONLY: get_ip3n
        USE timer_mod, ONLY: startTimer, stopTimer
        IMPLICIT NONE

        ! Subroutine arguments
        INTEGER, INTENT(in) :: ii, jj, kk
        INTEGER, INTENT(in) :: igrid
        INTEGER, INTENT(in) :: nxgrae, nygrae, nzgrae
        INTEGER, INTENT(in) :: npart_max
        INTEGER, INTENT(in) :: npart_grid
        INTEGER, INTENT(inout) :: iindex(npart_grid), jindex(npart_grid), &
            kindex(npart_grid)
        INTEGER, INTENT(in) :: part_idx(npart_max)
        REAL, INTENT(in) :: xpart(3,npart_max)
        INTEGER, INTENT(in) :: diamclass(npart_max)
        REAL, INTENT(in) :: upart(3,npart_max)
        REAL, INTENT(in) :: x(ii), y(jj), z(kk)
        REAL, INTENT(in) :: dx(ii),dy(jj),dz(kk)
        REAL, INTENT(inout) :: hilf(kk, jj, ii)
        REAL, INTENT(in) :: dt

        ! Local variables
        INTEGER :: i, j, k, l, m
        INTEGER :: ip3n
        REAL(kind=real64) :: alpha, beta

        REAL, POINTER :: npart_gp(:, :, :, :)
        REAL, POINTER :: upart_gp(:, :, :, :)
        REAL, POINTER :: vpart_gp(:, :, :, :)
        REAL, POINTER :: wpart_gp(:, :, :, :)

        REAL, POINTER :: a_npart_gp(:, :, :, :)
        REAL, POINTER :: a_upart_gp(:, :, :, :)
        REAL, POINTER :: a_vpart_gp(:, :, :, :)
        REAL, POINTER :: a_wpart_gp(:, :, :, :)

        CALL searchindex(kk, z, dz, nzgrae, npart_grid, 3, xpart, kindex, hilf, 0)
        CALL searchindex(jj, y, dy, nygrae, npart_grid, 2, xpart, jindex, hilf, 0)
        CALL searchindex(ii, x, dx, nxgrae, npart_grid, 1, xpart, iindex, hilf, 0)

        CALL get_ip3n(ip3n, nclasses, igrid)
        npart_gp(1:kk, 1:jj, 1:ii, 1:nclasses) => npart_g(ip3n:ip3n+kk*jj*ii*nclasses-1)
        upart_gp(1:kk, 1:jj, 1:ii, 1:nclasses) => upart_g(ip3n:ip3n+kk*jj*ii*nclasses-1)
        vpart_gp(1:kk, 1:jj, 1:ii, 1:nclasses) => vpart_g(ip3n:ip3n+kk*jj*ii*nclasses-1)
        wpart_gp(1:kk, 1:jj, 1:ii, 1:nclasses) => wpart_g(ip3n:ip3n+kk*jj*ii*nclasses-1)
        a_npart_gp(1:kk, 1:jj, 1:ii, 1:nclasses) => a_npart_g(ip3n:ip3n+kk*jj*ii*nclasses-1)
        a_upart_gp(1:kk, 1:jj, 1:ii, 1:nclasses) => a_upart_g(ip3n:ip3n+kk*jj*ii*nclasses-1)
        a_vpart_gp(1:kk, 1:jj, 1:ii, 1:nclasses) => a_vpart_g(ip3n:ip3n+kk*jj*ii*nclasses-1)
        a_wpart_gp(1:kk, 1:jj, 1:ii, 1:nclasses) => a_wpart_g(ip3n:ip3n+kk*jj*ii*nclasses-1)

        ! Set arrays to zero
        CALL dphi0n(kk*jj*ii*nclasses, npart_gp)
        CALL dphi0n(kk*jj*ii*nclasses, upart_gp)
        CALL dphi0n(kk*jj*ii*nclasses, vpart_gp)
        CALL dphi0n(kk*jj*ii*nclasses, wpart_gp)

        ! Number of particles per diameter class
        DO m = 1, npart_grid
            i = iindex(m)
            j = jindex(m)
            k = kindex(m)
            l = diamclass(m)
            npart_gp(k, j, i, l) = npart_gp(k, j, i, l) + 1.0
        END DO

        ! Average velocity per diameter class
        DO m = 1, npart_grid
            i = iindex(m)
            j = jindex(m)
            k = kindex(m)
            l = diamclass(m)

            upart_gp(k, j, i, l) = upart_gp(k, j, i, l) + upart(1, m)
            vpart_gp(k, j, i, l) = vpart_gp(k, j, i, l) + upart(2, m)
            wpart_gp(k, j, i, l) = wpart_gp(k, j, i, l) + upart(3, m)
        END DO

        ! Normalize
        DO l = 1, nclasses
            DO i = 1, ii
                DO j = 1, jj
                    DO k = 1, kk
                        upart_gp(k, j, i, l) = upart_gp(k, j, i, l)/MAX(1.0, npart_gp(k, j, i, l))
                        vpart_gp(k, j, i, l) = vpart_gp(k, j, i, l)/MAX(1.0, npart_gp(k, j, i, l))
                        wpart_gp(k, j, i, l) = wpart_gp(k, j, i, l)/MAX(1.0, npart_gp(k, j, i, l))
                    ENDDO
                ENDDO
            ENDDO
        END DO

        ! Time-weigted averages
        avg_time_tot = avg_time_tot + dt
        beta = dt/avg_time_tot
        alpha = 1.0 - beta
        CALL phiadd2(nclasses*ii*jj*kk, a_npart_gp, npart_gp, REAL(alpha, realk), REAL(beta, realk))
        CALL phiadd2(nclasses*ii*jj*kk, a_upart_gp, upart_gp, REAL(alpha, realk), REAL(beta, realk))
        CALL phiadd2(nclasses*ii*jj*kk, a_vpart_gp, vpart_gp, REAL(alpha, realk), REAL(beta, realk))
        CALL phiadd2(nclasses*ii*jj*kk, a_wpart_gp, wpart_gp, REAL(alpha, realk), REAL(beta, realk))
    END SUBROUTINE pp_statistics


    SUBROUTINE pp_statistics_init(idim3d)
        USE allocator_mod, ONLY: realFieldAlloc
        USE h5fieldio_mod, ONLY: h5fieldio_read_3dfield
        USE precision_mod, ONLY: realk

        ! Subroutine arguments
        INTEGER, INTENT(in) :: idim3d

#       include "mgpar.h"
#       include "pp_particles.h"

        ! Local variables
        CHARACTER(LEN=16) :: dsetname
        INTEGER :: l
        REAL :: rattr

        nclasses = ndiam_classes
        avg_time_tot = 0.0

        ! Allocate storage
        CALL realFieldAlloc(  npart_g, idim3d*nclasses)
        CALL realFieldAlloc(a_npart_g, idim3d*nclasses)

        CALL realFieldAlloc(  upart_g, idim3d*nclasses)
        CALL realFieldAlloc(  vpart_g, idim3d*nclasses)
        CALL realFieldAlloc(  wpart_g, idim3d*nclasses)

        CALL realFieldAlloc(a_upart_g, idim3d*nclasses)
        CALL realFieldAlloc(a_vpart_g, idim3d*nclasses)
        CALL realFieldAlloc(a_wpart_g, idim3d*nclasses)

        ! Read in stat fields
        IF (pp_lread .AND. pp_stat) THEN
            DO l = 1, nclasses
                WRITE(dsetname, "(A, I0, A)") "ANPART_", l, "_G"
                CALL h5fieldio_read_3dfield(a_npart_g, dsetname, &
                    rattr=rattr, comp=l, ncomp=nclasses)
            END DO

            DO l = 1, nclasses
                WRITE(dsetname, "(A, I0, A)") "AUPART_", l, "_G"
                CALL h5fieldio_read_3dfield(a_upart_g, dsetname, &
                    rattr=rattr, comp=l, ncomp=nclasses)
            END DO

            DO l = 1, nclasses
                WRITE(dsetname, "(A, I0, A)") "AVPART_", l, "_G"
                CALL h5fieldio_read_3dfield(a_vpart_g, dsetname, &
                    rattr=rattr, comp=l, ncomp=nclasses)
            END DO

            DO l = 1, nclasses
                WRITE(dsetname, "(A, I0, A)") "AWPART_", l, "_G"
                CALL h5fieldio_read_3dfield(a_wpart_g, dsetname, &
                    rattr=rattr, comp=l, ncomp=nclasses)
            END DO

            ! Assume value stored to each field is the same
            avg_time_tot = rattr
        END IF

    END SUBROUTINE pp_statistics_init


    SUBROUTINE pp_statistics_finish
        USE h5fieldio_mod, ONLY: h5fieldio_write_3dfield
        USE precision_mod, ONLY: realk

#       include "mgpar.h"
#       include "pp_particles.h"

        CHARACTER(LEN=16) :: dsetname, fieldname
        INTEGER :: stagi, stagj, stagk
        INTEGER :: l

        IF (pp_lwrite) THEN
            DO l = 1, nclasses
                WRITE(dsetname, "(A, I0, A)") "NPART_", l, "_G"
                WRITE(fieldname, "(A, I0, A)") "NPART", l, " G"
                stagi = 0; stagj = 0; stagk = 0

                CALL h5fieldio_write_3dfield(npart_g, dsetname, &
                    stagi, stagj, stagk, fieldname, comp=l, ncomp=nclasses)
            END DO

            DO l = 1, nclasses
                WRITE(dsetname, "(A, I0, A)") "UPART_", l, "_G"
                WRITE(fieldname, "(A, I0, A)") "UPART", l, " G"
                stagi = 0; stagj = 0; stagk = 0

                CALL h5fieldio_write_3dfield(upart_g, dsetname, &
                    stagi, stagj, stagk, fieldname, comp=l, ncomp=nclasses)
            END DO

            DO l = 1, nclasses
                WRITE(dsetname, "(A, I0, A)") "VPART_", l, "_G"
                WRITE(fieldname, "(A, I0, A)") "VPART", l, " G"
                stagi = 0; stagj = 0; stagk = 0

                CALL h5fieldio_write_3dfield(vpart_g, dsetname, &
                    stagi, stagj, stagk, fieldname, comp=l, ncomp=nclasses)
            END DO

            DO l = 1, nclasses
                WRITE(dsetname, "(A, I0, A)") "WPART_", l, "_G"
                WRITE(fieldname, "(A, I0, A)") "WPART", l, " G"
                stagi = 0; stagj = 0; stagk = 0

                CALL h5fieldio_write_3dfield(wpart_g, dsetname, &
                    stagi, stagj, stagk, fieldname, comp=l, ncomp=nclasses)
            END DO
        END IF

        IF (pp_lwrite .AND. pp_stat) THEN
            DO l = 1, nclasses
                WRITE(dsetname, "(A, I0, A)") "ANPART_", l, "_G"
                WRITE(fieldname, "(A, I0, A)") "<NPART", l, "> G"
                stagi = 0; stagj = 0; stagk = 0

                CALL h5fieldio_write_3dfield(a_npart_g, dsetname, &
                    stagi, stagj, stagk, fieldname, &
                    rattr=REAL(avg_time_tot, realk), comp=l, &
                    ncomp=nclasses)
            END DO

            DO l = 1, nclasses
                WRITE(dsetname, "(A, I0, A)") "AUPART_", l, "_G"
                WRITE(fieldname, "(A, I0, A)") "<UPART", l, "> G"
                stagi = 0; stagj = 0; stagk = 0

                CALL h5fieldio_write_3dfield(a_upart_g, dsetname, &
                    stagi, stagj, stagk, fieldname, &
                    rattr=REAL(avg_time_tot, realk), comp=l, &
                    ncomp=nclasses)
            END DO

            DO l = 1, nclasses
                WRITE(dsetname, "(A, I0, A)") "AVPART_", l, "_G"
                WRITE(fieldname, "(A, I0, A)") "<VPART", l, "> G"
                stagi = 0; stagj = 0; stagk = 0

                CALL h5fieldio_write_3dfield(a_vpart_g, dsetname, &
                    stagi, stagj, stagk, fieldname, &
                    rattr=REAL(avg_time_tot, realk), comp=l, &
                    ncomp=nclasses)
            END DO

            DO l = 1, nclasses
                WRITE(dsetname, "(A, I0, A)") "AWPART_", l, "_G"
                WRITE(fieldname, "(A, I0, A)") "<WPART", l, "> G"
                stagi = 0; stagj = 0; stagk = 0

                CALL h5fieldio_write_3dfield(a_wpart_g, dsetname, &
                    stagi, stagj, stagk, fieldname, &
                    rattr=REAL(avg_time_tot, realk), comp=l, &
                    ncomp=nclasses)
            END DO
        END IF
    END SUBROUTINE pp_statistics_finish

END MODULE pp_statistics_mod


     SUBROUTINE pp_errcheck(kk, jj, ii, igrid, nxgrae, nygrae, nzgrae, &
        npart_max, npart_grid, iindex, jindex, kindex, &
        part_idx, xpart, x, y, z, dx, dy, dz, hilf)

        USE precision_mod, ONLY: realk
        USE pointer_mod, ONLY: get_ip3n
        IMPLICIT NONE

        ! Subroutine arguments
        INTEGER, INTENT(in) :: kk, jj, ii, igrid
        INTEGER, INTENT(in) :: nxgrae, nygrae, nzgrae
        INTEGER, INTENT(in) :: npart_max
        INTEGER, INTENT(in) :: npart_grid
        INTEGER, INTENT(inout) :: iindex(npart_grid), jindex(npart_grid), &
            kindex(npart_grid)
        INTEGER, INTENT(in) :: part_idx(npart_max)
        REAL, INTENT(in) :: xpart(3,npart_max)
        REAL, INTENT(in) :: x(ii), y(jj), z(kk)
        REAL, INTENT(in) :: dx(ii),dy(jj),dz(kk)
        REAL, INTENT(inout) :: hilf(kk, jj, ii)

        ! Local variables
        INTEGER :: i, j, k, m

        CALL searchindex(kk, z, dz, nzgrae, npart_grid, xpart(3,:), kindex, hilf, 0)
        CALL searchindex(jj, y, dy, nygrae, npart_grid, xpart(2,:), jindex, hilf, 0)
        CALL searchindex(ii, x, dx, nxgrae, npart_grid, xpart(1,:), iindex, hilf, 0)

        ! Number of particles per diameter class
        DO m = 1, npart_grid
            i = iindex(m)
            j = jindex(m)
            k = kindex(m)
!            IF (igrid == 1) THEN
!                WRITE(*,*) m, igrid, xpart(1, m), xpart(2, m), xpart(3, m)
!            END IF
            IF (i == 0 .OR. j == 0 .OR. k == 0) THEN
                WRITE(*,*) "PP_CHECK ERROR:"
                WRITE(*,*) m, igrid, part_idx(m), i, j, k
                WRITE(*,*) m, igrid, xpart(1, m), xpart(2, m), xpart(3, m)
                CALL errr(501, 'pp_errcheck')
            END IF
        END DO
    END SUBROUTINE pp_errcheck


    SUBROUTINE pp_write_vtk(ii, jj, kk, igrid, nxgrae, nygrae, nzgrae, &
        npart_max, npart_grid, iindex, jindex, kindex, &
        part_idx, xpart, diamclass, upart, x, y, z, dx, dy, dz, hilf, ittot)

        USE precision_mod, ONLY: realk
        USE pointer_mod, ONLY: get_ip3n
        IMPLICIT NONE

        ! Subroutine arguments
        INTEGER, INTENT(in) :: ii, jj, kk
        INTEGER, INTENT(in) :: igrid
        INTEGER, INTENT(in) :: nxgrae, nygrae, nzgrae
        INTEGER, INTENT(in) :: npart_max
        INTEGER, INTENT(in) :: npart_grid
        INTEGER, INTENT(inout) :: iindex(npart_grid), jindex(npart_grid), &
            kindex(npart_grid)
        INTEGER, INTENT(in) :: part_idx(npart_max)
        REAL, INTENT(in) :: xpart(3,npart_max)
        INTEGER, INTENT(in) :: diamclass(npart_max)
        REAL, INTENT(in) :: upart(3,npart_max)
        REAL, INTENT(in) :: x(ii), y(jj), z(kk)
        REAL, INTENT(in) :: dx(ii),dy(jj),dz(kk)
        REAL, INTENT(inout) :: hilf(kk, jj, ii)
        INTEGER, INTENT(in) :: ittot

        ! Local variables
        INTEGER :: m
        INTEGER :: iou = 123
        CHARACTER(len=32) :: filename


        !IF (igrid /= 1) THEN
        !    RETURN
        !END IF
        IF (npart_grid == 0) THEN
            RETURN
        END IF

        WRITE(filename, '(A, I4.4, A, I4.4, A)') "PARTVTK-", igrid, "/part-", ittot, ".vtk"
        WRITE(*,*) ittot, filename
        OPEN(unit=iou, file=filename)

        WRITE(iou, '(A)') "# vtk DataFile Version 3.0"
        WRITE(iou, '(A)') "MGLET Particles"
        WRITE(iou, '(A)') "ASCII"
        WRITE(iou, '(A)') ""
        WRITE(iou, '(A)') "DATASET POLYDATA"

        WRITE(iou, '(A, I0, A)') "POINTS ", npart_grid, " float"
        DO m = 1, npart_grid
            WRITE(iou, *) xpart(1,m), xpart(2,m), xpart(3,m)
        END DO

        WRITE(iou, '(A, I0)') "VERTICES 1 ", npart_grid+1
        WRITE(iou, *)  npart_grid
        DO m = 1, npart_grid
            WRITE(iou, *) m-1
        END DO

!        WRITE(iou, '(A, I0, A, I0)') "CELLS ", npart_grid, " ", npart_grid*2
!        DO m = 1, npart_grid
!            WRITE(iou, *) 1, m
!        END DO

        WRITE(iou, '(A, I0)') "POINT_DATA ", npart_grid
        WRITE(iou, '(A)') "SCALARS idx float 1"
        WRITE(iou, '(A)') "LOOKUP_TABLE default"
        DO m = 1, npart_grid
            WRITE(iou, *) part_idx(m)
        END DO

        WRITE(iou, '(A)') "VECTORS velocity float"
        DO m = 1, npart_grid
            WRITE(iou, *) upart(1,m), upart(2,m), upart(3,m)
        END DO

        CLOSE(iou)
    END SUBROUTINE pp_write_vtk
