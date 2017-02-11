MODULE fourier_val
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & calctype, & ! (0) Lanczos (1) TPQ (2) FullDiag (3) LOBCG (4) mVMC
  & nwfc, & ! Number of state
  & box(3,3), & ! Supercell index
  & rbox(3,3), & ! Reciplocal Superlattice Vector times nk
  & nsite,   & ! Number of sites
  & ncor1,   & ! Nomber of One-body Correlation function
  & ncor2,   & ! Number of Two-body Correlation function
  & ncor(8), & ! Number of Correlation function for each index(See below)
  & nktot,   & ! Total number of k
  & nk_row,  & ! number row of total k
  & nk         ! Number of k to be computed
  !
  REAL(8),SAVE :: &
  & direct(3,3), & ! Direct lattice vector
  & recipr(3,3), & ! Reciprocal lattice vector
  & koff(3)        ! k-point offset for the boundary phase
  !
  CHARACTER(256),SAVE :: &
  & filehead, & ! Filename header for correlation functions
  & file_one, & ! Filename for One-body Correlation
  & file_two    ! Filename for Two-body Correlation
  !
  INTEGER,ALLOCATABLE,SAVE :: &
  & indx(:,:,:), & ! (nsite:nsite,8) Mapping index for each Correlation function
  & equiv(:)      ! (nktot) Equivalent k point in the 1st BZ
  !
  REAL(8),ALLOCATABLE,SAVE :: &
  & kvec(:,:), & ! (3,nk) k-vector in the 1st BZ
  & kvec_tot(:,:), & ! (3,nktot) k-vector in the lerger area
  & site(:,:)    ! (3,nsite) Site geometry in the fractional coordinate
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & cor(:,:,:,:), & ! (nsite,nsite,6,nwfc) Correlation function in real space (See below)
  & cor_k(:,:,:)    ! (nk,6,nwfc) Correlation function in the k-space (See below)
  !
  CHARACTER(256),ALLOCATABLE :: &
  & filetail(:) ! (nwfc) Index (run, step, etc.) in a file name
  !
  ! Kind of Correlation for ncor(1:8) and index(:,:,1:8)
  !
  ! (1) C_{i up  }A_{j up  }
  ! (2) C_{i down}A_{j down}
  ! (3) N_{i up  }N_{j up  } = C_{i up  }A_{i up  }C_{j up  }A_{j up  }
  ! (4) N_{i up  }N_{j down} = C_{i up  }A_{i up  }C_{j down}A_{j down}
  ! (5) N_{i down}N_{j up  } = C_{i down}A_{i down}C_{j up  }A_{j up  }
  ! (6) N_{i down}N_{j down} = C_{i down}A_{i down}C_{j down}A_{j down}
  ! (7) S_{i +   }S_{j -   } = C_{i up  }A_{i down}C_{j down}A_{j up  }
  ! (8) S_{i -   }S_{j +   } = C_{i down}A_{i up  }C_{j up  }A_{j down}
  !
  ! Kind of Correlation for cor(:,:,1:6) and cor_k(:,1:6)
  !
  ! (1) C_{i up  }A_{j up  }
  ! (2) C_{i down}A_{j down}
  ! (3) N_{i}N_{J} = N_{i up}N_{j up} + N_{i up}N_{j down} + N_{i down}N_{j up} + N_{i down}N_{j down}
  ! (4) Sz_{i}Sz_{J} = 0.25*(N_{i up}N_{j up} - N_{i up}N_{j down} - N_{i down}N_{j up} + N_{i down}N_{j down})
  ! (5) S_{i +   }S_{j -   }
  ! (6) S_{i -   }S_{j +   }
  !
END MODULE fourier_val
!
!
!
MODULE fourier_routine
  !
  IMPLICIT NONE
  !
  INTERFACE
     SUBROUTINE key2lower(key) BIND(c)
       USE,INTRINSIC :: iso_c_binding
       CHARACTER(KIND=C_CHAR) :: key(*)
     END SUBROUTINE key2lower
  END INTERFACE
  !
CONTAINS
!
! Read from HPhi/mVMC input files
!
SUBROUTINE read_filename()
  !
  USE fourier_val, ONLY : file_one, file_two, filehead, nsite, nwfc, filetail, calctype
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, lanczos_max, numave, interval, irun, istep, iwfc, idx_start
  CHARACTER(256) :: modpara, calcmod, keyname, namelist
  !
  WRITE(*,*) 
  WRITE(*,*) "#####  Read HPhi/mVMC Input Files  #####" 
  WRITE(*,*) 
  !
  namelist = ""
  CALL GETARG(1, namelist)
  !
  ! Read from NameList file
  !
  OPEN(fi,file = TRIM(namelist))
  !
  calcmod = ""
  DO
     READ(fi,*,END=10) keyname
     BACKSPACE(fi)
     CALL key2lower(keyname)
     !
     IF(TRIM(ADJUSTL(keyname)) == "onebodyg") THEN
        READ(fi,*) keyname, file_one
     ELSE IF(TRIM(ADJUSTL(keyname)) == "twobodyg") THEN
        READ(fi,*) keyname, file_two
     ELSE IF(TRIM(ADJUSTL(keyname)) == "modpara") THEN
        READ(fi,*) keyname, modpara
     ELSE IF(TRIM(ADJUSTL(keyname)) == "calcmod") THEN
        READ(fi,*) keyname, calcmod
     ELSE
        READ(fi,*) keyname
     END IF
  END DO
  !
10 CONTINUE
  WRITE(*,*) "  Read from ", TRIM(namelist)
  CLOSE(FI)
  !
  WRITE(*,*) "    OneBodyG file : ", TRIM(ADJUSTL(file_one))
  WRITE(*,*) "    TwoBodyG file : ", TRIM(ADJUSTL(file_two))
  WRITE(*,*) "    ModPara file : ", TRIM(ADJUSTL(modpara)) 
  WRITE(*,*) "    CalcMod file : ", TRIM(ADJUSTL(calcmod))
  !
  ! Read from Modpara file
  !
  OPEN(fi,file = TRIM(modpara))
  !
  DO
     READ(fi,*,END=20) keyname
     BACKSPACE(fi)
     CALL key2lower(keyname)
     !
     IF(TRIM(ADJUSTL(keyname)) == "nsite") THEN
        READ(fi,*) keyname, nsite
     ELSE IF(TRIM(ADJUSTL(keyname)) == "cdatafilehead") THEN
        READ(fi,*) keyname, filehead
     ELSE IF(TRIM(ADJUSTL(keyname)) == "numave") THEN
        READ(fi,*) keyname, numave
     ELSE IF(TRIM(ADJUSTL(keyname)) == "lanczos_max") THEN
        READ(fi,*) keyname, lanczos_max
     ELSE IF(TRIM(ADJUSTL(keyname)) == "expecinterval") THEN
        READ(fi,*) keyname, interval
     ELSE IF(TRIM(ADJUSTL(keyname)) == "exct") THEN
        READ(fi,*) keyname, nwfc
     ELSE IF(TRIM(ADJUSTL(keyname)) == "ndataidxstart") THEN
        READ(fi,*) keyname, idx_start
     ELSE IF(TRIM(ADJUSTL(keyname)) == "ndataqtysmp") THEN
        READ(fi,*) keyname, numave
     ELSE
        READ(fi,*) keyname
     END IF
  END DO
  !
20 CONTINUE
  WRITE(*,*) "  Read from ", TRIM(modpara)
  WRITE(*,*) "    FileHead : ", TRIM(ADJUSTL(filehead))
  WRITE(*,*) "    Number of site : ", nsite
  CLOSE(FI)
  !
  filehead = "output/" // TRIM(ADJUSTL(filehead))
  !
  ! Read from CalcMod file
  !
  IF(calcmod == "") THEN
     !
     ! When mVMC
     !
     calctype = 4
     !
  ELSE
     !
     OPEN(fi,file = TRIM(calcmod))
     !
     DO
        READ(fi,*,END=30) keyname
        BACKSPACE(fi)
        CALL key2lower(keyname)
        !
        IF(TRIM(ADJUSTL(keyname)) == "calctype") THEN
           READ(fi,*) keyname, calctype
        ELSE
           READ(fi,*) keyname
        END IF
     END DO
     !
30   CONTINUE
     WRITE(*,*) "  Read from ", TRIM(calcmod)
     WRITE(*,*) "    CalcType : ", calctype
     CLOSE(FI)
     !
  END IF ! IF(calcmod == "") THEN
  !
  ! Spefify the tail of file name
  !
  IF(calctype == 0) THEN
     !
     ! Lanczos
     !
     nwfc = 1
     ALLOCATE(filetail(nwfc))
     filetail = ".dat"
     !
     WRITE(*,*) "    Method : Lanczos"
     !
  ELSE IF (calctype == 1) THEN
     !
     ! TPQ
     !
     nwfc = numave * (1 + (lanczos_max - 1) / interval)
     ALLOCATE(filetail(nwfc))
     !
     iwfc = 0
     DO irun = 0, numave - 1
        DO istep = 0, lanczos_max - 1
           IF(MOD(istep, interval) == 0) THEN
              iwfc = iwfc + 1
              WRITE(filetail(iwfc),'(a,i0,a,i0,a)') &
              & "_set", irun, "step", istep, ".dat"
           END IF
        END DO
     END DO
     !
     WRITE(*,*) "    Method : TPQ"
     WRITE(*,*) "    Number of Run : ", numave
     WRITE(*,*) "    Maximum Iteration : ", lanczos_max
     WRITE(*,*) "    Expectation Interval : ", interval
     !
  ELSE IF (calctype == 2) THEN
     !
     ! FullDiag
     !
     OPEN(fi, file = "output/CHECK_Memory.dat")
     READ(fi, '("  MAX DIMENSION idim_max=1", i16)') nwfc
     CLOSE(fi)
     ALLOCATE(filetail(nwfc))
     !
     DO iwfc = 1, nwfc
        WRITE(filetail(iwfc),'(a,i0,a)') "_eigen", iwfc - 1, ".dat"
     END DO
     !
     WRITE(*,*) "    Method : Full Diagonalization"
     !
  ELSE IF (calctype == 3) THEN
     !
     ! LOBCG
     !
     ALLOCATE(filetail(nwfc))
     DO iwfc = 1, nwfc
        WRITE(filetail(iwfc),'(a,i0,a)') "_eigen", iwfc - 1, ".dat"
     END DO
     !
     WRITE(*,*) "    Method : LOBCG"
     !
  ELSE
     !
     ! mVMC
     !
     nwfc = numave
     ALLOCATE(filetail(nwfc))
     DO iwfc = idx_start, idx_start + numave - 1
        WRITE(filetail(iwfc),'(a,i3.3,a)') "_", iwfc, ".dat"
     END DO
     !
     WRITE(*,*) "    Method : mVMC"
     WRITE(*,*) "    Start Index : ", idx_start
     !
  END IF ! (clactype == ??)
  !
  WRITE(*,*) "    Number of States : ", nwfc
  !
END SUBROUTINE read_filename
!
! Read geometry from file
!
SUBROUTINE read_geometry()
  !
  USE fourier_val, ONLY : direct, recipr, box, rbox, nsite, site, nk, koff
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, isite, ii, jj
  REAL(8) :: det, phase(3), pi180 = ACOS(-1d0)/180d0
  CHARACTER(256) :: filename
  !
  ALLOCATE(site(3,nsite))
  !
  WRITE(*,*) 
  WRITE(*,*) "#####  Read Geometry Input File  #####" 
  WRITE(*,*) 
  !
  ! "nd argument should be geometry file name
  !
  CALL GETARG(2, filename)
  !
  WRITE(*,*) "  Read from ", TRIM(filename)
  OPEN(fi, file = TRIM(filename))
  !
  ! Direct lattice vector in arbitraly unit
  !
  DO ii = 1, 3
     READ(fi,*) direct(1:3,ii)
  END DO
  WRITE(*,*) "    Direct LATTICE VECTOR :"
  WRITE(*,'(4x3f15.10)') direct(1:3, 1:3)
  !
  ! Bondary phase
  !
  READ(fi,*) phase(1:3)  
  WRITE(*,*) "    Boundary phase[degree] : "
  WRITE(*,'(4x3f15.10)') phase(1:3)
  !
  ! Supercell index (a0w, a0l, a1w, a1l)
  !
  DO ii = 1, 3
     READ(fi,*) box(1:3,ii)
  END DO
  WRITE(*,*) "    Supercell Index :"
  WRITE(*,'(3i8)') box(1:3, 1:3)
  !
  DO isite = 1, nsite
     READ(fi,*) site(1:3,isite)
  END DO
  !
  CLOSE(fi)
  !
  ! The number of k in the 1st BZ
  !
  nk = 0
  DO ii = 1, 3
     nk = nk &
     &  + box(ii,1) * box(MOD(ii,     3) + 1, 2) * box(MOD(ii + 1, 3) + 1, 3) &
     &  - box(ii,1) * box(MOD(ii + 1, 3) + 1, 2) * box(MOD(ii,     3) + 1, 3)
  END DO
  !
  ! Compute reciprocal SuperLattice Vector
  !
  DO ii = 1, 3
     DO jj = 1, 3
        rbox(jj, ii) = &
        &   box(MOD(jj,     3) + 1, MOD(ii, 3) + 1) * box(MOD(jj + 1, 3) + 1, MOD(ii + 1, 3) + 1) &
        & - box(MOD(jj + 1, 3) + 1, MOD(ii, 3) + 1) * box(MOD(jj,     3) + 1, MOD(ii + 1, 3) + 1)
     END DO
  END DO
  !
  IF(nk < 0) THEN
     nk = -nk
     rbox(1:3,1:3) = - rbox(1:3,1:3)
  END IF
  WRITE(*,*) "    Number of k point : ", nk
  WRITE(*,*) "    Reciprocal superlattice vector (times nk) :"
  WRITE(*,'(3i8)') rbox(1:3, 1:3)
  !
  ! Compute Reciprocal Lattice Vector
  !
  det = 0d0
  DO ii = 1, 3
     det = det &
     &  + direct(ii,1) * direct(MOD(ii,     3) + 1, 2) * direct(MOD(ii + 1, 3) + 1, 3) &
     &  - direct(ii,1) * direct(MOD(ii + 1, 3) + 1, 2) * direct(MOD(ii,     3) + 1, 3)
  END DO
  !
  ! Compute reciprocal SuperLattice Vector
  !
  DO ii = 1, 3
     DO jj = 1, 3
        recipr(jj, ii) = &
        &   direct(MOD(jj,     3) + 1, MOD(ii, 3) + 1) * direct(MOD(jj + 1, 3) + 1, MOD(ii + 1, 3) + 1) &
        & - direct(MOD(jj + 1, 3) + 1, MOD(ii, 3) + 1) * direct(MOD(jj,     3) + 1, MOD(ii + 1, 3) + 1)
     END DO
  END DO
  recipr(1:3, 1:3) = recipr(1:3, 1:3) / det
  WRITE(*,*) "    Reciplocal lattice vector :"
  WRITE(*,'(4x3f15.10)') recipr(1:3, 1:3)
  !
  ! k-point offset for the boundary phase
  !
  phase(1:3) = phase(1:3) * pi180
  koff(1:3) = MATMUL(DBLE(rbox(1:3,1:3)), phase(1:3))
  koff(1:3) = koff(1:3) / DBLE(nk)
  koff(1:3) = MATMUL(recipr(1:3,1:3), koff(1:3))
  WRITE(*,*) "    k-point offset :"
  WRITE(*,'(4x3f15.10)') koff(1:3)
  !  
END SUBROUTINE read_geometry
!
! Set k points
!
SUBROUTINE set_kpoints()
  !
  USE fourier_val, ONLY : box, rbox, nk, nktot, kvec, kvec_tot, equiv, nk_row
  !
  IMPLICIT NONE
  !
  INTEGER :: imax(3), imin(3), i1, i2, i3, ii, edge(3,8), nk0, ik, jk, idim, ikvec0(3)
  INTEGER,ALLOCATABLE :: ikvec(:,:), ikvec_tot(:,:)
  REAL(8) :: kvec0(3)
  !
  ! Define range of k-grid index spanning [-1:1] in fractional BZ
  !
  ii = 0
  DO i3 = -1, 1, 2
     DO i2 = -1, 1, 2
        DO i1 = -1, 1, 2
           ii = ii + 1
           edge(1:3,ii) = MATMUL((/i1, i2, i3/), box(1:3,1:3))
        END DO
     END DO
  END DO
  !
  DO idim = 1, 3
     imin(idim) = MINVAL(edge(idim,1:8))
     imax(idim) = MAXVAL(edge(idim,1:8))
  END DO
  !
  nktot = PRODUCT(imax(1:3) - imin(1:3) + 1)
  nk_row = imax(1) - imin(1) + 1
  ALLOCATE(kvec(3,nk), ikvec(3,nk))
  ALLOCATE(kvec_tot(3,nktot), ikvec_tot(3,nktot), equiv(nktot))
  !
  nk0 = 0
  nktot = 0
  DO i3 = imin(3), imax(3)
     DO i2 = imin(2), imax(2)
        DO i1 = imin(1), imax(1)
           !
           ikvec0(1:3) = MATMUL(rbox(1:3,1:3), (/i1, i2, i3/))
           kvec0(1:3) = DBLE(ikvec0(1:3)) / DBLE(nk)
           !
           ! Only k-vectors in the 1st BZ is used in the fourier trans.
           !
           IF(ALL(0 <= ikvec0(1:3)) .AND. ALL(ikvec0(1:3) < nk)) THEN
              nk0 = nk0 + 1
              ikvec(1:3,nk0) = ikvec0(1:3)
              kvec( 1:3,nk0) = kvec0(1:3)
           END IF
           !
           IF(i3 == 0) THEN
              nktot = nktot + 1
              ikvec_tot(1:3,nktot) = ikvec0(1:3)
              kvec_tot( 1:3,nktot) = kvec0(1:3)
           END IF
           !
        END DO
     END DO
  END DO
  WRITE(*,*) "  Number of k : ", nk0 ! Must be the same as nk
  !
  ! Search equivalent k-point
  !
  DO ik = 1, nktot
     !
     DO jk = 1, nk
        !
        IF(ALL(MODULO(ikvec_tot(1:3,ik), nk) == ikvec(1:3,jk))) THEN
           equiv(ik) = jk
           EXIT
        END IF
        !
     END DO ! jk = 1, nk
     !
  END DO ! ik = 1, nktot
  !
END SUBROUTINE set_kpoints
!
! Read Correlation Function
!
SUBROUTINE read_corrindx()
  !
  USE fourier_val, ONLY : file_one, file_two, ncor1, ncor2, ncor, indx, nsite
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, itmp(8), icor, itmp0(2)
  CHARACTER(100) :: ctmp
  !
  WRITE(*,*) 
  WRITE(*,*) "#####  Read Correlation Index File  #####" 
  WRITE(*,*) 
  !
  ALLOCATE(indx(nsite,nsite,8))
  indx(1:nsite,1:nsite,1:8) = 0
  !
  ! Read index for the One-Body Correlation
  !
  OPEN(fi, file = TRIM(file_one))
  WRITE(*,*) "  Read from ", TRIM(file_one)
  READ(fi,*) ctmp
  READ(fi,*) ctmp, ncor1
  READ(fi,*) ctmp
  READ(fi,*) ctmp
  READ(fi,*) ctmp
  WRITE(*,*) "    Number of Correlation Function (One body) : ", ncor1
  !
  ! Search mapping index for up-up and down-down correlation
  !
  ncor(1:2) = 0
  DO icor = 1, ncor1
     READ(fi,*) itmp(1:4)
     IF(itmp(2) == 0 .AND. itmp(4) == 0) THEN
        !
        ! Up-Up correlation
        !
        indx(itmp(1) + 1, itmp(3) + 1, 1) = icor 
     ELSE IF (itmp(2) == 1 .AND. itmp(4) == 1) THEN
        !
        ! Down-Down correlation
        !
        indx(itmp(1) + 1, itmp(3) + 1, 2) = icor 
     END IF
  END DO
  !
  WRITE(*,*) "    Number of Up-Up Index : ", COUNT(indx(1:nsite, 1:nsite, 1) /= 0)
  WRITE(*,*) "    Number of Down-Down Index : ", COUNT(indx(1:nsite, 1:nsite, 2) /= 0)
  !
  CLOSE(fi)
  !
  ! Read index for the Two-Body Correlation
  !
  OPEN(fi, file = TRIM(file_two))
  WRITE(*,*) "  Read from ", TRIM(file_two)
  READ(fi,*) ctmp
  READ(fi,*) ctmp, ncor2
  READ(fi,*) ctmp
  READ(fi,*) ctmp
  READ(fi,*) ctmp
  WRITE(*,*) "    Number of Correlation Function (Two body) : ", ncor2
  !
  ! Search index for up-up, up-up-down-down, etc.
  !
  ncor(3:8) = 0
  DO icor = 1, ncor2
     !
     READ(fi,*) itmp(1:8)
     !
     IF(itmp(1) == itmp(3) .AND. itmp(5) == itmp(7)) THEN
        !
        IF(itmp(2) == 0 .AND. itmp(4) == 0) THEN
           !
           IF(itmp(6) == 0 .AND. itmp(8) == 0) THEN
              !
              ! UpUpUpUp
              !
              indx(itmp(1) + 1, itmp(5) + 1, 3) = icor 
           ELSE IF(itmp(6) == 1 .AND. itmp(8) == 1) THEN
              !
              ! UpUpDownDown
              !
              indx(itmp(1) + 1, itmp(5) + 1, 4) = icor 
           END IF
        ELSE IF(itmp(2) == 1 .AND. itmp(4) == 1) THEN
           !
           IF(itmp(6) == 0 .AND. itmp(8) == 0) THEN
              !
              ! DownDownUpUp
              !
              indx(itmp(1) + 1, itmp(5) + 1, 5) = icor 
           ELSE IF(itmp(6) == 1 .AND. itmp(8) == 1) THEN
              !
              ! DownDownDownDown
              !
              indx(itmp(1) + 1, itmp(5) + 1, 6) = icor 
           END IF
        END IF
        !
     END IF
     !
     IF(itmp(1) == itmp(7) .AND. itmp(3) == itmp(5)) THEN
        !
        IF((itmp(2) == 0 .AND. itmp(4) == 0) .AND. (itmp(6) == 1 .AND. itmp(8) == 1)) THEN
           !
           ! Up-Down-Down-Up = S+S-
           !
           indx(itmp(1) + 1, itmp(5) + 1, 7) = icor 
        ELSE IF((itmp(2) == 1 .AND. itmp(4) == 1) .AND. (itmp(6) == 0 .AND. itmp(8) == 0)) THEN
           !
           ! Down-Up-Up-Down = S-S+
           !
           indx(itmp(1) + 1, itmp(5) + 1, 8) = icor 
        END IF
        !
     END IF
     !
  END DO
  !
  WRITE(*,*) "    Number of UpUpUpUp         Index : ", COUNT(indx(1:nsite, 1:nsite, 3) /= 0)
  WRITE(*,*) "    Number of UpUpDownDown     Index : ", COUNT(indx(1:nsite, 1:nsite, 4) /= 0)
  WRITE(*,*) "    Number of DownDownUpUp     Index : ", COUNT(indx(1:nsite, 1:nsite, 5) /= 0)
  WRITE(*,*) "    Number of DownDownDownDown Index : ", COUNT(indx(1:nsite, 1:nsite, 6) /= 0)
  WRITE(*,*) "    Number of Plus-Minus       Index : ", COUNT(indx(1:nsite, 1:nsite, 7) /= 0)
  WRITE(*,*) "    Number of Minus-Plus       Index : ", COUNT(indx(1:nsite, 1:nsite, 8) /= 0)
  !
  CLOSE(fi)
  !
END SUBROUTINE read_corrindx
!
! Read Correlation Function
!
SUBROUTINE read_corrfile()
  !
  USE fourier_val, ONLY : filehead, filetail, nwfc, calctype, &  
  &                       ncor1, ncor2, ncor, indx, cor, nsite
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, icor, itmp(8), ii, iwfc, isite, jsite
  COMPLEX(8),ALLOCATABLE :: cor0(:)
  REAL(8),ALLOCATABLE :: cor0_r(:,:)
  CHARACTER(256) :: filename
  !
  ALLOCATE(cor(nsite,nsite,6,nwfc))
  ALLOCATE(cor0(0:MAX(ncor1,ncor2)), cor0_r(2,MAX(ncor1,ncor2)))
  cor(1:nsite,1:nsite,1:6,1:nwfc) = CMPLX(0d0, 0d0, KIND(1d0))
  cor0(0) = CMPLX(0d0, 0d0, KIND(1d0))
  !
  DO iwfc = 1, nwfc
     !
     ! One Body Correlation Function
     !
     filename = TRIM(filehead) // "_cisajs" // TRIM(filetail(iwfc))
     OPEN(fi, file = TRIM(filename))
     !
     IF(calctype == 4) THEN
        !
        ! mVMC
        !
        READ(fi,*) cor0_r(1:2, 1:ncor1)
     ELSE
        !
        ! HPhi
        !
        DO icor = 1, ncor1
           READ(fi,*) itmp(1:4), cor0_r(1:2, icor)
        END DO
     END IF ! (calctype == 4)
     !
     CLOSE(fi)
     !
     cor0(1:ncor1) = CMPLX(cor0_r(1,1:ncor1), cor0_r(2,1:ncor1), KIND(1d0))
     !
     ! Map it into Up-Up(1) and Down-Down(2) Correlation
     !
     DO ii = 1, 2
        DO jsite = 1, nsite
           DO isite = 1, nsite
              cor(isite, jsite, ii, iwfc) = cor0(indx(isite,jsite,ii))
           END DO
        END DO
     END DO
     !
     ! Two Body Correlation function
     !
     filename = TRIM(filehead) // "_cisajscktalt" // TRIM(filetail(iwfc))
     OPEN(fi, file = TRIM(filename))
     !
     IF(calctype == 4) THEN
        !
        ! mVMC
        !
        READ(fi,*) cor0_r(1:2, 1:ncor2)
     ELSE
        !
        ! HPhi
        !
        DO icor = 1, ncor2
           READ(fi,*) itmp(1:8), cor0_r(1:2, icor)
        END DO
     END IF
     !
     CLOSE(fi)
     !
     cor0(1:ncor2) = CMPLX(cor0_r(1,1:ncor2), cor0_r(2,1:ncor2), KIND(1d0))
     !
     ! Map it into Density-Density(3), Sz-Sz(4), S+S-(5), S-S+(6) Correlation
     !
     ! Up-Up-Up-Up and Down-Down-Down-Down into Density-Density & Sz-Sz
     !
     DO jsite = 1, nsite
        DO isite = 1, nsite
           !
           cor(isite, jsite, 3, iwfc) = cor0(indx(isite,jsite,3)) &
           &                          + cor0(indx(isite,jsite,4)) &
           &                          + cor0(indx(isite,jsite,5)) &
           &                          + cor0(indx(isite,jsite,6))
           !
           cor(isite, jsite, 4, iwfc) = cor0(indx(isite,jsite,3)) &
           &                          - cor0(indx(isite,jsite,4)) &
           &                          - cor0(indx(isite,jsite,5)) &
           &                          + cor0(indx(isite,jsite,6))
           !
        END DO
     END DO
     !
     cor(1:nsite,1:nsite, 4, iwfc) = cor(1:nsite,1:nsite, 4, iwfc) * 0.25d0
     !
     ! Up-Down-Down-Up(S+S-) and Down-Up-Up-Down(S-S+)
     !
     DO ii = 5, 6
        DO jsite = 1, nsite
           DO isite = 1, nsite
              cor(isite, jsite, ii, iwfc) = cor0(indx(isite,jsite,ii + 2))
           END DO
           !
           ! Ciu+ Cid Cjd+ Cju = delta_{ij} Ciu+ Ciu - Ciu+ Cju Cjd+ Cid
           ! Cid+ Ciu Cju+ Cjd = delta_{ij} Cid+ Cid - Cid+ Cjd Cju+ Ciu
           !
           cor(jsite, jsite, ii, iwfc) = cor(isite, isite, ii, iwfc) &
           & + cor(isite, isite, ii - 4, iwfc) 
           !
        END DO
     END DO
     !
  END DO ! iwfc = 1, nwfc
  !
  DEALLOCATE(cor0, cor0_r, indx)
  !
END SUBROUTINE read_corrfile
!
! Fourier transformation
!
SUBROUTINE fourier_cor()
  !
  USE fourier_val, ONLY : nsite, cor, cor_k, site, kvec, nwfc, nk
  IMPLICIT NONE
  !
  INTEGER :: isite, jsite, ik
  REAL(8) :: tpi = 2.0 * ACOS(-1d0), theta
  COMPLEX(8),ALLOCATABLE :: fmat(:,:,:)
  !
  ALLOCATE(fmat(nk,nsite,nsite), cor_k(nk,6,nwfc))
  !
  ! Matirx for Fourier trans. exp(-i k R)
  !
  DO jsite = 1, nsite
     DO isite = 1, nsite
        DO ik = 1, nk
           theta = - tpi * DOT_PRODUCT(kvec(1:3,ik), (site(1:3,isite) - site(1:3,jsite)))
           fmat(ik,isite,jsite) = CMPLX(COS(theta), SIN(theta), KIND(1d0))
        END DO
     END DO
  END DO
  !
  CALL zgemm('N', 'N', nk, 6*nwfc, nsite*nsite, CMPLX(1d0, 0d0, KIND(1d0)), fmat, nk, &
  &          cor, nsite*nsite, CMPLX(0d0,0d0,KIND(1d0)), cor_k, nk)
  !
  cor_k(1:nk,1:6,1:nwfc) = cor_k(1:nk,1:6,1:nwfc) / dble(nsite**2)
  !
  DEALLOCATE(fmat, cor, site)
  !
END SUBROUTINE fourier_cor
!
! Output Fourier component of Correlation function
!
SUBROUTINE output_cor()
  !
  USE fourier_val, ONLY : cor_k, nk, nktot, nk_row, kvec, kvec_tot, koff, &
  &                       equiv, nwfc, recipr, filehead, filetail, calctype
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, ik, iwfc, idim
  REAL(8) :: tpi = 2.0 * ACOS(-1d0)
  CHARACTER(256) :: filename
  COMPLEX(8),ALLOCATABLE :: cor_ave(:,:), cor_err(:,:)
  !
  ! Output Correlation function in the 1st BZ
  !
  WRITE(*,*) 
  WRITE(*,*) "#####  Output Files  #####" 
  WRITE(*,*) 
  !
  WRITE(*,*) "  Correlation in k-space : ", TRIM(filehead) // "_corr", "*.dat"
  !
  IF(calctype == 4) THEN
     !
     ! mVMC
     !
     ALLOCATE(cor_ave(nk,6), cor_err(nk,6))
     !
     ! Average
     !
     cor_ave(1:nk,1:6) = SUM(cor_k(1:nk,1:6,1:nwfc), 3) / DBLE(nwfc)
     !
     ! Variance
     !
     cor_err(1:nk,1:6) = 0d0
     DO iwfc = 1, nwfc
        cor_err(1:nk,1:6) = CMPLX( DBLE(cor_k(1:nk,1:6,iwfc) - cor_ave(1:nk,1:6))**2, &
        &                         AIMAG(cor_k(1:nk,1:6,iwfc) - cor_ave(1:nk,1:6))**2, &
        &                         KIND(0d0))
     END DO
     !
     ! Standard Error
     !
     IF(nwfc == 1) THEN
        cor_err(1:nk,1:6) = CMPLX(0d0, 0d0, KIND(0d0))
     ELSE
        cor_err(1:nk,1:6) = CMPLX(SQRT( DBLE(cor_err(1:nk,1:6))), &
        &                         SQRT(AIMAG(cor_err(1:nk,1:6))), KIND(0d0)) &
        &                 / (DBLE(nwfc - 1) * SQRT(DBLE(nwfc)))
     END IF
     !
     filename = TRIM(filehead) // "_corr.dat"
     OPEN(fo, file = TRIM(filename))
     !
     WRITE(fo,*) "#mVMC", nk
     WRITE(fo,*) "# kx[1] ky[2] kz[3](Cart.) UpUp[4,5,16,17] (Re. Im. Err.) DownDown[6,7,18,19]"
     WRITE(fo,*) "# Density[8,9,20,21] SzSz[10,11,22,23] S+S-[12,13,24,25] S-S+[14,15,26.27]"
     WRITE(fo,'(a,3f15.7)') " #k-offset", koff(1:3)
     !
     DO ik = 1, nk
        WRITE(fo,'(27e15.5)') tpi * MATMUL(recipr(1:3,1:3), kvec(1:3,ik)), &
        &                     cor_ave(ik,1:6), cor_err(ik,1:6)
     END DO
     !
     CLOSE(fo)
     !
     DEALLOCATE(cor_ave, cor_err)
     !
  ELSE
     !
     ! HPhi
     !
     DO iwfc = 1, nwfc
        !
        filename = TRIM(filehead) // "_corr" // TRIM(filetail(iwfc))
        OPEN(fo, file = TRIM(filename))
        !
        WRITE(fo,*) "#HPhi", nk
        WRITE(fo,*) "# kx[1] ky[2] kz[3](Cart.) UpUp[4,5] (Re. Im.) DownDown[6,7]"
        WRITE(fo,*) "# Density[8,9] SzSz[10,11] S+S-[12,13] S-S+[14,15]"
        WRITE(fo,'(a,3f15.7)') " #k-offset", koff(1:3)
        !
        DO ik = 1, nk
           WRITE(fo,'(15e15.5)') tpi * MATMUL(recipr(1:3,1:3), kvec(1:3,ik)), cor_k(ik,1:6,iwfc)
        END DO
        !
        CLOSE(fo)
        !
     END DO
     !
  END IF ! IF(calctype == 4)
  !
  ! k-points in the larger area
  !
  OPEN(fo, file = "kpoints.dat")
  !
  WRITE(fo,*) nktot, nk_row
  DO idim = 1, 3
     WRITE(fo,'(3e15.5)') tpi * recipr(1:3,idim)
  END DO
  !
  DO ik = 1, nktot
     WRITE(fo,'(3e15.5,i7)') tpi * MATMUL(recipr(1:3,1:3), kvec_tot(1:3,ik)), equiv(ik)
  END DO
  !
  CLOSE(fo)
  !
  DEALLOCATE(cor_k,kvec,kvec_tot,equiv)
  !
END SUBROUTINE output_cor
!
END MODULE fourier_routine
!
! Main routine
!
PROGRAM fourier
  !
  USE fourier_routine, ONLY : read_filename, read_geometry, set_kpoints, &
  &                           read_corrindx, read_corrfile, fourier_cor, output_cor
  IMPLICIT NONE
  !
  CALL read_filename()
  CALL read_geometry()
  CALL set_kpoints()
  CALL read_corrindx()
  CALL read_corrfile()
  CALL fourier_cor()
  CALL output_cor()
  !
  WRITE(*,*) 
  WRITE(*,*) "#####  Done  #####" 
  WRITE(*,*) 
  !
END PROGRAM fourier
