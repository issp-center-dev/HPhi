MODULE corplot_val
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & nline, &
  & itarget, &
  & nwfc, & ! Number of state
  & nk_row,  & ! number row of total k
  & nk         ! Number of k to be computed
  !
  REAL(8),SAVE :: &
  & koff(3), &
  & bragg(3,26), &
  & braggnorm(26), &
  & bz_line(3,2,26*26), &
  & recipr(3,3)    ! Reciprocal lattice vector
  !
  LOGICAL,SAVE :: &
  & errbar, &
  & rpart
  !
  REAL(8),ALLOCATABLE,SAVE :: &
  & kvec(:,:) ! (2,nk) k-vector in the 1st BZ
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & cor_err(:,:,:), & ! (nk,6,nwfc) Correlation function in the k-space (See below)
  & cor_k(:,:,:)    ! (nk,6,nwfc) Correlation function in the k-space (See below)
  !
  ! Kind of Correlation for cor_k(:,1:6)
  !
  ! (1) C_{i up  }A_{j up  }
  ! (2) C_{i down}A_{j down}
  ! (3) N_{i}N_{J} = N_{i up}N_{j up} + N_{i up}N_{j down} + N_{i down}N_{j up} + N_{i down}N_{j down}
  ! (4) Sz_{i}Sz_{J} = 0.25*(N_{i up}N_{j up} - N_{i up}N_{j down} - N_{i down}N_{j up} + N_{i down}N_{j down})
  ! (5) S_{i +   }S_{j -   }
  ! (6) S.S
  !
END MODULE corplot_val
!
!
!
MODULE corplot_routine
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Output Fourier component of Correlation function
!
SUBROUTINE read_cor()
  !
#if defined(FUJITSU)
  USE service_routines, ONLY : IARGC
#endif
  USE corplot_val, ONLY : nwfc, nk_row, recipr, koff, &
  &                       cor_k, cor_err, kvec, nk
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ik, iwfc, nk0, idim
  REAL(8) :: rtmp(3)
  CHARACTER(256) :: filename, ctmp1, ctmp2
#if defined(SR)
  INTEGER,INTRINSIC :: IARGC
#endif
  !
  WRITE(*,*) 
  WRITE(*,*) "#####  Read Files  #####" 
  WRITE(*,*) 
  !
  nwfc = IARGC()
  WRITE(*,*) "  Number of Files : ", nwfc
  !
  DO iwfc = 1, nwfc
     !
     CALL GETARG(iwfc, filename)
     !
     OPEN(fi, file = TRIM(filename))
     !
     READ(fi,*) ctmp1, nk0, nk_row
     DO idim = 1, 3
        READ(fi,*) ctmp2, recipr(1:3,idim)
     END DO
     READ(fi,*) ctmp2
     READ(fi,*) ctmp2
     READ(fi,*) ctmp2, koff(1:3)
     !
     IF(iwfc == 1) THEN
        !
        nk = nk0
        ALLOCATE(cor_k(1:nk,1:6,1:nwfc), cor_err(1:nk,1:6,1:nwfc), kvec(3,nk))
        cor_k(1:nk,1:6,1:nwfc) = CMPLX(0d0, 0d0, KIND(0d0))
        cor_err(1:nk,1:6,1:nwfc) = 0d0
        !
        WRITE(*,*) "    Number of k : ", nk
        WRITE(*,*) "    k-point offset :"
        WRITE(*,'(4x3f15.10)') koff(1:3)
        !
     END IF
     !
     IF(TRIM(ctmp1) == "#HPhi") THEN
        WRITE(*,'(a,a,a)') "  Read ", TRIM(filename), " as HPhi Correlation File"
        DO ik = 1, nk
           READ(fi,'(15e15.5)') kvec(1:3, ik), cor_k(ik,1:6,iwfc)
        END DO
     ELSE ! mVMC
        WRITE(*,'(a,a,a)') "  Read ", TRIM(filename), " as mVMC Correlation File"
        DO ik = 1, nk
           READ(fi,'(27e15.5)') kvec(1:3, ik), cor_k(ik,1:6,iwfc), cor_err(ik,1:6,iwfc)
        END DO
     END IF
     !
     CLOSE(fi)
     !
  END DO ! iwfc = 1, nwfc
  !
END SUBROUTINE read_cor
!
! Set vectors to difine Bragg's plane
!
SUBROUTINE set_bragg_vector()
  !
  USE corplot_val, ONLY : recipr, bragg, braggnorm
  IMPLICIT NONE
  !
  INTEGER :: i0, i1, i2, ibr
  !
  ibr = 0
  !
  DO i0 = -1, 1
     DO i1 = -1, 1
        DO i2 = -1, 1
           !
           IF(ALL((/i0, i1, i2/) == 0)) CYCLE
           !
           ibr = ibr + 1
           bragg(1:3,ibr) = MATMUL(recipr(1:3,1:3), DBLE((/i0, i1, i2/))) * 0.5d0
           !
           braggnorm(ibr) = DOT_PRODUCT(bragg(1:3,ibr), bragg(1:3,ibr))
           !
        END DO
     END DO
  END DO
  !
END SUBROUTINE set_bragg_vector
!
! Solve linear system
!
FUNCTION solve3(a, b) RESULT(det)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: a(3,3)
  REAL(8),INTENT(INOUT) :: b(3)
  !
  REAL(8) :: det, c(3)
  !
  det = a(1, 1) * (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2)) &
  &   + a(1, 2) * (a(2, 3) * a(3, 1) - a(2, 1) * a(3, 3)) &
  &   + a(1, 3) * (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1))
  !
  c(1) = b(1) * (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2)) &
  &    + b(2) * (a(1, 3) * a(3, 2) - a(1, 2) * a(3, 3)) &
  &    + b(3) * (a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2))
  !
  c(2) = b(1) * (a(2, 3) * a(3, 1) - a(2, 1) * a(3, 3)) &
  &    + b(2) * (a(1, 1) * a(3, 3) - a(1, 3) * a(3, 1)) &
  &    + b(3) * (a(1, 3) * a(2, 1) - a(1, 1) * a(2, 3))
  !
  c(3) = b(1) * (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1)) &
  &    + b(2) * (a(1, 2) * a(3, 1) - a(1, 1) * a(3, 2)) &
  &    + b(3) * (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1))
  !
  b(1:3) = c(1:3) / det
  !
  RETURN
  !
END FUNCTION solve3
!
! Judge wheser this line is the edge of 1st BZ
!
FUNCTION bz_corners(ibr, jbr, nbr, corner, corner2) RESULT(lbr)
  !
  USE corplot_val, ONLY : bragg, braggnorm
  IMPLICIT NONE
  !
  INTEGER :: ibr, jbr, nbr
  REAL(8) :: corner(3)
  REAL(8) :: corner2(3)
  !
  INTEGER :: kbr, i, lbr, nbr0
  REAL(8) :: bmat(3,3), rhs(3), prod, thr = 1d-4, det
  !
  nbr0 = nbr
  !
  DO kbr = nbr0, 26
     !
     bmat(1,1:3) = bragg(1:3,ibr)
     bmat(2,1:3) = bragg(1:3,jbr)
     bmat(3,1:3) = bragg(1:3,kbr)
     !
     rhs(1) = braggnorm(ibr)
     rhs(2) = braggnorm(jbr)
     rhs(3) = braggnorm(kbr)
     !
     ! if Bragg planes do not cross, roop next kbr
     !
     det = solve3(bmat, rhs)
     IF (ABS(det) < thr) CYCLE
     !
     ! if vert0 = vert1, roop next kbr
     !
     prod = DOT_PRODUCT(corner2(1:3) - rhs(1:3), corner2(1:3) - rhs(1:3))
     IF (prod < thr) CYCLE
     !
     ! is this corner really in 1st BZ ?
     !
     i = 0
     DO lbr = 1, 26
        !
        prod = DOT_PRODUCT(bragg(1:3,lbr), rhs(1:3))
        !
        IF (prod > braggnorm(lbr) + thr) THEN
           i = 1
           EXIT
        END if
        !
     END DO
     !
     IF(i /= 1) THEN
        corner(1:3) = rhs(1:3)
        lbr = kbr + 1
        RETURN
     END IF
     !
  END DO
  !
  !  this line is not a BZ boundary
  !
  lbr = 1
  RETURN
  !
END FUNCTION bz_corners
!
! Compute Brillouin zone boundariy lines
!
SUBROUTINE set_bz_line()
  !
  USE corplot_val, ONLY : bz_line, nline
  IMPLICIT NONE
  !
  INTEGER :: ibr, jbr, nbr, lvert
  REAL(8) :: corner(3,2)
  !
  CALL set_bragg_vector()
  !
  nline = 0
  !
  DO ibr = 1, 26
     DO jbr = 1, 26
        !
        corner(1:3,1:2) = 0d0
        nbr = 1
        lvert = bz_corners(ibr, jbr, nbr, corner(1:3,1), corner(1:3,2))
        IF(lvert == 1) CYCLE
        nbr = lvert
        !
        lvert = bz_corners(ibr, jbr, nbr, corner(1:3,2), corner(1:3,1))
        IF(lvert == 1) CYCLE
        !
        nline = nline + 1
        !
        ! Sort
        !
        IF(corner(1,1) - corner(1,2) > 0.000001) THEN
           bz_line(1:3,1:2,nline) = corner(1:3,1:2)
        ELSE IF (corner(1,2) - corner(1,1) > 0.000001) THEN
           bz_line(1:3,1,nline) = corner(1:3,2)
           bz_line(1:3,2,nline) = corner(1:3,1)
        ELSE
           IF(corner(2,1) - corner(2,2) > 0.000001) THEN
              bz_line(1:3,1:2,nline) = corner(1:3,1:2)
           ELSE IF (corner(2,2) - corner(2,1) > 0.000001) THEN
              bz_line(1:3,1,nline) = corner(1:3,2)
              bz_line(1:3,2,nline) = corner(1:3,1)
           ELSE
              IF(corner(3,1) - corner(3,2) > 0.000001) THEN
                 bz_line(1:3,1:2,nline) = corner(1:3,1:2)
              ELSE 
                 bz_line(1:3,1,nline) = corner(1:3,2)
                 bz_line(1:3,2,nline) = corner(1:3,1)
              END IF
           END IF
        END IF
        !
     END DO
  END DO
  !
END SUBROUTINE set_bz_line
!
! Compute Unique BZ line
!
SUBROUTINE uniq_bz_line(minz,maxz,nline2,bz_line2)
  !
  USE corplot_val, ONLY : bz_line, nline
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: minz, maxz
  INTEGER,INTENT(OUT) :: nline2
  REAL(8),INTENT(OUT) :: bz_line2(3,2,nline*4)
  !
  INTEGER :: iline, jline
  REAL(8) :: bz_line0(3,2)
  !
  nline2 = 0
  DO iline = 1, nline
     !
     bz_line0(1:3,1) = (/bz_line(1,1,iline), bz_line(2,1,iline), maxz/)
     bz_line0(1:3,2) = (/bz_line(1,1,iline), bz_line(2,1,iline), minz/)
     !
     DO jline = 1, nline2
        IF(ALL(ABS(bz_line0(1:3,1:2) - bz_line2(1:3,1:2,jline)) < 0.000001)) GOTO 10
     END DO
     !
     nline2 = nline2 + 1
     bz_line2(1:3,1:2,nline2) = bz_line0(1:3,1:2)
     !
10   CONTINUE
     !
     bz_line0(1:3,1) = (/bz_line(1,2,iline), bz_line(2,2,iline), maxz/)
     bz_line0(1:3,2) = (/bz_line(1,2,iline), bz_line(2,2,iline), minz/)
     !
     DO jline = 1, nline2
        IF(ALL(ABS(bz_line0(1:3,1:2) - bz_line2(1:3,1:2,jline)) < 0.000001)) GOTO 20
     END DO
     !
     nline2 = nline2 + 1
     bz_line2(1:3,1:2,nline2) = bz_line0(1:3,1:2)
     !
20   CONTINUE
     !
     IF(ALL(ABS(bz_line(1:2,1,iline) - bz_line(1:2,2,iline)) < 0.000001)) CYCLE
     !
     bz_line0(1:3,1) = (/bz_line(1,1,iline), bz_line(2,1,iline), minz/)
     bz_line0(1:3,2) = (/bz_line(1,2,iline), bz_line(2,2,iline), minz/)
     !
     DO jline = 1, nline2
        IF(ALL(ABS(bz_line0(1:3,1:2) - bz_line2(1:3,1:2,jline)) < 0.000001)) GOTO 30
     END DO
     !
     nline2 = nline2 + 1
     bz_line2(1:3,1:2,nline2) = bz_line0(1:3,1:2)
     !
30   CONTINUE
     !
     bz_line0(1:3,1) = (/bz_line(1,1,iline), bz_line(2,1,iline), maxz/)
     bz_line0(1:3,2) = (/bz_line(1,2,iline), bz_line(2,2,iline), maxz/)
     !
     DO jline = 1, nline2
        IF(ALL(ABS(bz_line0(1:3,1:2) - bz_line2(1:3,1:2,jline)) < 0.000001)) GOTO 40
     END DO
     !
     nline2 = nline2 + 1
     bz_line2(1:3,1:2,nline2) = bz_line0(1:3,1:2)
     !
40   CONTINUE
     !
  END DO
  !
END SUBROUTINE uniq_bz_line
!
! Write gnuplot script
!
SUBROUTINE write_gnuplot()
  !
  USE corplot_val, ONLY : itarget, rpart, errbar, nwfc, &
  &                       cor_k, nk, nwfc, nline
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, iwfc, iline, nline2
  REAL(8) :: maxz, minz, bz_line2(3,2,nline*4)
  !
  IF(rpart) THEN
     maxz = MAXVAL(DBLE( cor_k(1:nk, itarget, 1:nwfc)))
     minz = MINVAL(DBLE( cor_k(1:nk, itarget, 1:nwfc)))
  ELSE
     maxz = MAXVAL(AIMAG(cor_k(1:nk, itarget, 1:nwfc)))
     minz = MINVAL(AIMAG(cor_k(1:nk, itarget, 1:nwfc)))
  END IF
  CALL uniq_bz_line(minz,maxz,nline2,bz_line2)
  !
  OPEN(fo, file = "correlation.gp")
  !
  WRITE(fo,'(a)') "#set terminal pdf color enhanced \" !"
  WRITE(fo,'(a)') "#dashed dl 1.0 size 20.0cm, 20.0cm"
  WRITE(fo,'(a)') "#set output 'correlation.pdf'"
  WRITE(fo,'(a)') "#set view 60.0, 30.0"
  !
  WRITE(fo,*) 
  WRITE(fo,'(a)') "#set view equal xy"
  WRITE(fo,'(a)') "set ticslevel 0"
  WRITE(fo,'(a)') "set hidden3d"
  WRITE(fo,'(a)') "set xlabel 'kx'"
  WRITE(fo,'(a)') "set ylabel 'ky'"
  WRITE(fo,'(a,e15.5,a,e15.5,a)') "set zrange [", minz, ":", maxz, "]"
  WRITE(fo,*) 
  WRITE(fo,'(a)') "set pm3d"
  WRITE(fo,'(a)') "set pm3d interpolate 5, 5"
  WRITE(fo,'(a)') "#set contour"
  WRITE(fo,'(a)') "set view 0.0, 0.0"
  WRITE(fo,*) 
  WRITE(fo,*) "#####  Set Brillouin-Zone Boundary  #####"
  WRITE(fo,*) 
  DO iline = 1, nline2
     WRITE(fo,'(a,e15.5,a,e15.5,a,e15.5,a,e15.5,a,e15.5,a,e15.5,a)') &
     &       "set arrow from ", bz_line2(1,1,iline), ",", bz_line2(2,1,iline), ",", bz_line2(3,1,iline), &
     &                  " to ", bz_line2(1,2,iline), ",", bz_line2(2,2,iline), ",", bz_line2(3,2,iline), " nohead front"
  END DO ! iline = 1, nline
  !
  WRITE(fo,*) 
  WRITE(fo,*) "#####  End Set Brillouin-Zone Boundary  #####"
  WRITE(fo,*) 
  WRITE(fo,'(a)') "splot \" !"
  IF(errbar) THEN
     DO iwfc = 1, nwfc
        !
        WRITE(fo, '(a,i0,a,i0,a,i0,a)',advance='no') &
        & "'correlation.dat' u 1:2:($", iwfc+2, "-$", iwfc+2+nwfc, ") w l tit '", iwfc, "-', "
        WRITE(fo, '(a,i0,a,i0,a,i0,a)',advance='no') &
        & "'correlation.dat' u 1:2:($", iwfc+2, "+$", iwfc+2+nwfc, ") w l tit '", iwfc, "+'"
        !
        IF(iwfc /= nwfc) WRITE(fo,'(a)') ", \" !"
        !
     END DO
  ELSE
     DO iwfc = 1, nwfc
        !
        WRITE(fo, '(a,i0,a,i0,a)',advance='no') &
        & "'correlation.dat' u 1:2:", iwfc+2, " w l tit '", iwfc, "'"
        !
        IF(iwfc /= nwfc) WRITE(fo,'(a)') ", \" !"
        !
     END DO
  END IF
  WRITE(fo,*)
  WRITE(fo,*) "pause -1"
  !
  CLOSE(fo)
  !
END SUBROUTINE write_gnuplot
!
! Output data to be plotted
!
SUBROUTINE write_data()
  !
  USE corplot_val, ONLY : itarget, rpart, nwfc, nk_row, &
  &                       cor_k, cor_err, kvec, koff, nk
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, ik
  REAL(8) :: koff2(2)
  CHARACTER(100) :: form
  !
  OPEN(fo, file = "correlation.dat")
  WRITE(form,'(a,i0,a)') "(", 2 + nwfc*2, "e15.5)"
  !
  IF(itarget <= 2) THEN
     koff2(1:2) = koff(1:2)
  ELSE
     koff2(1:2) = 0d0
  END IF
  !
  DO ik = 1, nk
     !
     IF(rpart) THEN
        WRITE(fo,TRIM(form)) kvec(1:2,ik) + koff2(1:2), &
        &  DBLE(cor_k(  ik, itarget, 1:nwfc)), &
        &  DBLE(cor_err(ik, itarget, 1:nwfc))
     ELSE
        WRITE(fo,TRIM(form)) kvec(1:2,ik) + koff2(1:2), &
        &  AIMAG(cor_k(  ik, itarget, 1:nwfc)), &
        &  AIMAG(cor_err(ik, itarget, 1:nwfc))
     END IF
     !
     IF(MOD(ik, nk_row) == 0) WRITE(fo,*)
     !
  END DO
  !
  CLOSE(fo)
  !
END SUBROUTINE write_data
!
END MODULE corplot_routine
!
! Main routine
!
PROGRAM corplot
  !
  USE corplot_routine, ONLY : set_bz_line, read_cor, write_gnuplot, write_data
  USE corplot_val, ONLY : itarget, errbar, rpart
  IMPLICIT NONE
  !
  CALL read_cor()
  CALL set_bz_line()
  !
  WRITE(*,*) 
  WRITE(*,*) "#####  Plot Start  #####"
  WRITE(*,*) 
  WRITE(*,*) "  Please specify target number from below (0 or Ctrl-C to exit): "
  WRITE(*,*) 
  WRITE(*,*) "  Real Part Without ErrorBar"
  WRITE(*,*) "    [ 1] Up-Up [ 2] Down-Down [ 3] Density-Density [ 4] SzSz [ 5] S+S- [ 6] S.S"
  WRITE(*,*) "  Imaginary Part Without ErrorBar"
  WRITE(*,*) "    [11] Up-Up [12] Down-Down [13] Density-Density [14] SzSz [15] S+S- [16] S.S"
  WRITE(*,*) "  Real Part With ErrorBar"
  WRITE(*,*) "    [21] Up-Up [22] Down-Down [23] Density-Density [24] SzSz [25] S+S- [26] S.S"
  WRITE(*,*) "  Imaginary Part With ErrorBar"
  WRITE(*,*) "    [31] Up-Up [32] Down-Down [33] Density-Density [34] SzSz [35] S+S- [36] S.S"
  WRITE(*,*) 
  !
  DO
     WRITE(*,'(a)',ADVANCE='NO') "  Target : "
     READ(*,*) itarget
     errbar = itarget / 20 >= 1
     rpart = MOD(itarget / 10, 2) == 0
     itarget = MOD(itarget, 10)
     IF(itarget < 1 .OR. 6 < itarget) EXIT
     CALL write_data()
     CALL write_gnuplot()
     CALL SYSTEM("gnuplot correlation.gp")
  END DO
  !
  WRITE(*,*) 
  WRITE(*,*) "#####  Done  #####"
  WRITE(*,*) 
  !
END PROGRAM corplot
