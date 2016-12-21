MODULE corplot_val
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & itarget, &
  & nwfc, & ! Number of state
  & nktot,   & ! Total number of k
  & nk_row,  & ! number row of total k
  & nk         ! Number of k to be computed
  !
  REAL(8),SAVE :: &
  & recipr(2,2)    ! Reciprocal lattice vector
  !
  LOGICAL,SAVE :: &
  & errbar, &
  & rpart
  !
  INTEGER,ALLOCATABLE,SAVE :: &
  & equiv(:)      ! (nktot) Equivalent k point in the 1st BZ
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
  ! (6) S_{i -   }S_{j +   }
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
  USE corplot_val, ONLY : nwfc, nktot, nk_row, &
  &                       cor_k, cor_err, equiv, kvec
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ik, iwfc, nk0, nk
  REAL(8) :: rtmp(2)
  CHARACTER(256) :: filename, ctmp1, ctmp2
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
     READ(fi,*) ctmp1, nk0
     READ(fi,*) ctmp2
     READ(fi,*) ctmp2
     IF(iwfc == 1) THEN
        nk = nk0
        ALLOCATE(cor_k(1:nk,1:6,1:nwfc), cor_err(1:nk,1:6,1:nwfc))
        cor_k(1:nk,1:6,1:nwfc) = CMPLX(0d0, 0d0, KIND(0d0))
        cor_err(1:nk,1:6,1:nwfc) = 0d0
     END IF
     !
     IF(TRIM(ctmp1) == "#HPhi") THEN
        WRITE(*,'(a,a,a)') "  Read ", TRIM(filename), " as HPhi Correlation File"
        DO ik = 1, nk
           READ(fi,'(14e15.5)') rtmp(1:2), cor_k(ik,1:6,iwfc)
        END DO
     ELSE ! mVMC
        WRITE(*,'(a,a,a)') "  Read ", TRIM(filename), " as mVMC Correlation File"
        DO ik = 1, nk
           READ(fi,'(26e15.5)') rtmp(1:2), cor_k(ik,1:6,iwfc), cor_err(ik,1:6,iwfc)
        END DO
     END IF
     !
     CLOSE(fi)
     !
  END DO ! iwfc = 1, nwfc
  !
  ! k-points in the larger area
  !
  OPEN(fi, file = "kpoints.dat")
  !
  READ(fi,*) nktot, nk_row
  ALLOCATE(kvec(2,nktot), equiv(nktot))
  WRITE(*,*) "  Total Number of k : ", nktot
  WRITE(*,*) "  Row for k : ", nk_row
  !
  DO ik = 1, nktot
     READ(fi,*) kvec(1:2,ik), equiv(ik)
  END DO
  !
  CLOSE(fi)
  !
END SUBROUTINE read_cor
!
! Write gnuplot script
!
SUBROUTINE write_gnuplot()
  !
  USE corplot_val, ONLY : itarget, rpart, errbar, nwfc
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, ik, iwfc
  !
  OPEN(fo, file = "correlation.gp")
  !
  WRITE(fo,*) "set ticslevel 0"
  WRITE(fo,*) "set hidden3d"
  WRITE(fo,*) "set xlabel 'kx'"
  WRITE(fo,*) "set ylabel 'ky'"
  WRITE(fo,*) "#set pm3d"
  WRITE(fo,*) "#set contour"
  WRITE(fo,'(a)') "splot \" !"
  IF(errbar) THEN
     DO iwfc = 1, nwfc
        WRITE(fo, '(a,i0,a,i0,a,i0a)') &
        & "'correlation.dat' u 1:2:($", iwfc+2, "-$", iwfc+2+nwfc, ") w l tit '", iwfc, "-', \" !"
        WRITE(fo, '(a,i0,a,i0,a,i0a)') &
        & "'correlation.dat' u 1:2:($", iwfc+2, "+$", iwfc+2+nwfc, ") w l tit '", iwfc, "+', \" !"
     END DO
  ELSE
     DO iwfc = 1, nwfc
        WRITE(fo, '(a,i0,a,i0a)') &
        & "'correlation.dat' u 1:2:", iwfc+2, " w l tit '", iwfc, "', \" !"
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
  USE corplot_val, ONLY : itarget, rpart, nwfc, nktot, nk_row, &
  &                       cor_k, cor_err, equiv, kvec
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, ik
  CHARACTER(100) :: form
  !
  OPEN(fo, file = "correlation.dat")
  WRITE(form,'(a,i0,a)') "(", 2 + nwfc*2, "e15.5)"
  !
  DO ik = 1, nktot
     !
     IF(rpart) THEN
        WRITE(fo,TRIM(form)) kvec(1:2,ik),  DBLE(cor_k(equiv(ik), itarget, 1:nwfc)), &
        &                                 DBLE(cor_err(equiv(ik), itarget, 1:nwfc))
     ELSE
        WRITE(fo,TRIM(form)) kvec(1:2,ik), AIMAG(cor_k(equiv(ik), itarget, 1:nwfc)), &
        &                                AIMAG(cor_err(equiv(ik), itarget, 1:nwfc))
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
  USE corplot_routine, ONLY : read_cor, write_gnuplot, write_data
  USE corplot_val, ONLY : itarget, errbar, rpart
  IMPLICIT NONE
  !
  CALL read_cor()
  !
  WRITE(*,*) 
  WRITE(*,*) "#####  Plot Start  #####"
  WRITE(*,*) 
  WRITE(*,*) "  Please specify target number from below (0 or Ctrl-C to exit): "
  WRITE(*,*) 
  WRITE(*,*) "  Real Part Without ErrorBar"
  WRITE(*,*) "    [ 1] Up-Up [ 2] Down-Down [ 3] Density-Density [ 4] SzSz [ 5] S+S- [ 6] S-S+"
  WRITE(*,*) "  Imaginary Part Without ErrorBar"
  WRITE(*,*) "    [11] Up-Up [12] Down-Down [13] Density-Density [14] SzSz [15] S+S- [16] S-S+"
  WRITE(*,*) "  Real Part With ErrorBar"
  WRITE(*,*) "    [21] Up-Up [22] Down-Down [23] Density-Density [24] SzSz [25] S+S- [26] S-S+"
  WRITE(*,*) "  Imaginary Part With ErrorBar"
  WRITE(*,*) "    [31] Up-Up [32] Down-Down [33] Density-Density [34] SzSz [35] S+S- [36] S-S+"
  WRITE(*,*) 
  !
  DO
     WRITE(*,'(a)',ADVANCE='NO') "  Target : "
     READ(*,*) itarget
     errbar = itarget / 20 >= 1
     rpart = MOD(itarget / 10, 2) == 0
     itarget = MOD(itarget, 10)
     IF(itarget < 1 .OR. 6 < itarget) EXIT
     CALL write_gnuplot()
     CALL write_data()
     CALL SYSTEM("gnuplot correlation.gp")
  END DO
  !
  WRITE(*,*) 
  WRITE(*,*) "#####  Done  #####"
  WRITE(*,*) 
  !
END PROGRAM corplot
