MODULE fourier_val
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & calctype, & ! (0) Lanczos (1) TPQ (2) FullDiag (3) LOBCG (4) mVMC
  & nwfc, & ! Number of state
  & avec(2,2), & ! Supercell index
  & bvec(2,2), & ! Reciplocal Superlattice Vector times nk
  & nsite,   & ! Number of sites
  & ncor1,   & ! Nomber of One-body Correlation function
  & ncor2,   & ! Number of Two-body Correlation function
  & ncor(8), & ! Number of Correlation function for each index(See below)
  & nktot,   & ! Total number of k
  & nk_row,  & ! number row of total k
  & nk         ! Number of k to be computed
  !
  REAL(8),SAVE :: &
  & direct(2,2), & ! Direct lattice vector
  & recipr(2,2)    ! Reciprocal lattice vector
  !
  CHARACTER(256),SAVE :: &
  & filehead, & ! Filename header for correlation functions
  & file_one, & ! Filename for One-body Correlation
  & file_two    ! Filename for Two-body Correlation
  !
  INTEGER,ALLOCATABLE,SAVE :: &
  & indx(:,:,:), & ! (3,nsite*nsite,8) Mapping index for each Correlation function
  & equiv(:)      ! (nktot) Equivalent k point in the 1st BZ
  !
  REAL(8),ALLOCATABLE,SAVE :: &
  & kvec(:,:), & ! (2,nk) k-vector in the 1st BZ
  & kvec_tot(:,:), & ! (2,nktot) k-vector in the lerger area
  & site(:,:)    ! (2,nsite) Site geometry in the fractional coordinate
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
     filehead = "output/" // TRIM(ADJUSTL(filehead))
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
  IF(calctype == 0) THEN ! Lanczos
     !
     nwfc = 1
     ALLOCATE(filetail(nwfc))
     filetail = ".dat"
     !
     WRITE(*,*) "    Method : Lanczos"
     !
  ELSE IF (calctype == 1) THEN ! TPQ
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
  ELSE IF (calctype == 2) THEN ! FullDiag
     !
     OPEN(fi, file = "output/CHECK_Memory.dat")
     READ(fi, '("  MAX DIMENSION idim_max=1", i16)') nwfc
     CLOSE(fi)
     ALLOCATE(filetail(nwfc))
     !
     DO iwfc = 1, nwfc
        WRITE(filetail(iwfc),'(a,i0,a)') "_eigen", iwfc, ".dat"
     END DO
     !
     WRITE(*,*) "    Method : Full Diagonalization"
     !
  ELSE IF (calctype == 3) THEN ! LOBCG
     !
     ALLOCATE(filetail(nwfc))
     DO iwfc = 1, nwfc
        WRITE(filetail(iwfc),'(a,i0,a)') "_eigen", iwfc, ".dat"
     END DO
     !
     WRITE(*,*) "    Method : LOBCG"
     !
  ELSE ! mVMC
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
  USE fourier_val, ONLY : direct, recipr, avec, bvec, nsite, site, nk
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, isite
  REAL(8) :: det
  CHARACTER(256) :: filename
  !
  ALLOCATE(site(2,nsite))
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
  READ(fi,*) direct(1:2,1)
  READ(fi,*) direct(1:2,2)
  WRITE(*,*) "    Direct LATTICE VECTOR :"
  WRITE(*,'(4x2f15.10)') direct(1:2, 1:2)
  !
  ! Supercell index (a0w, a0l, a1w, a1l)
  !
  READ(fi,*) avec(1:2,1)
  READ(fi,*) avec(1:2,2)
  WRITE(*,*) "    Supercell Index :"
  WRITE(*,'(2i8)') avec(1:2, 1:2)
  !
  DO isite = 1, nsite
     READ(fi,*) site(1:2,isite)
  END DO
  !
  CLOSE(fi)
  !
  ! Compute reciprocal SuperLattice Vector
  !
  nk = avec(1,1) * avec(2,2) - avec(1,2) * avec(2,1)
  bvec(1:2,1) = (/  avec(2,2), - avec(1,2)/)
  bvec(1:2,2) = (/- avec(2,1),   avec(1,1)/)
  IF(nk < 0) THEN
     nk = -nk
     bvec(1:2,1:2) = - bvec(1:2,1:2)
  END IF
  WRITE(*,*) "    Number of k point : ", nk
  WRITE(*,*) "    Reciprocal superlattice vector (times nk) :"
  WRITE(*,'(2i8)') bvec(1:2, 1:2)
  !
  ! Compute Reciprocal Lattice Vector
  !
  det = direct(1,1) * direct(2,2) - direct(1,2) * direct(2,1)
  recipr(1:2,1) = (/  direct(2,2), - direct(1,2)/) / det
  recipr(1:2,2) = (/- direct(2,1),   direct(1,1)/) / det
  WRITE(*,*) "    Reciplocal lattice vector :"
  WRITE(*,'(4x2f15.10)') recipr(1:2, 1:2)
  !
  RETURN
  !
10 WRITE(*,*) "Stop in reading namelist FOURIER !"
  STOP
  !
END SUBROUTINE read_geometry
!
! Set k points
!
SUBROUTINE set_kpoints()
  !
  USE fourier_val, ONLY : avec, bvec, nk, nktot, kvec, kvec_tot, equiv, nk_row
  !
  IMPLICIT NONE
  !
  INTEGER :: imax(2), imin(2), i1, i2, ii, edge(2,4), nk0, ik, jk, idim
  INTEGER,ALLOCATABLE :: ikvec(:,:), ikvec_tot(:,:)
  !
  ! Define range of k-grid index spanning [-1:1] in fractional BZ
  !
  ii = 0
  DO i2 = -1, 1, 2
     DO i1 = -1, 1, 2
        ii = ii + 1
        edge(1:2,ii) = MATMUL((/i1, i2/), avec(1:2,1:2))
     END DO
  END DO
  !
  DO idim = 1, 2
     imin(idim) = MINVAL(edge(idim,1:4))
     imax(idim) = MAXVAL(edge(idim,1:4))
  END DO
  !
  nktot = PRODUCT(imax(1:2) - imin(1:2) + 1)
  nk_row = imax(1) - imin(1) + 1
  ALLOCATE(kvec(2,nk), ikvec(2,nk))
  ALLOCATE(kvec_tot(2,nktot), ikvec_tot(2,nktot), equiv(nktot))
  !
  nk0 = 0
  nktot = 0
  DO i2 = imin(2), imax(2)
     DO i1 = imin(1), imax(1)
        !
        nktot = nktot + 1
        ikvec_tot(1:2,nktot) = MATMUL(bvec(1:2,1:2), (/i1, i2/))
        kvec_tot( 1:2,nktot) = DBLE(ikvec_tot(1:2,nktot)) / DBLE(nk)
        !
        ! Only k-vectors in the 1st BZ is used in the fourier trans.
        !
        IF(ALL(0 <= ikvec_tot(1:2,nktot)) .AND. ALL(ikvec_tot(1:2,nktot) < nk)) THEN
           nk0 = nk0 + 1
           ikvec(1:2,nk0) = ikvec_tot(1:2,nktot)
           kvec( 1:2,nk0) = kvec_tot( 1:2,nktot)
        END IF
        !
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
        IF(ALL(MODULO(ikvec_tot(1:2,ik), nk) == ikvec(1:2,jk))) THEN
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
  INTEGER :: fi = 10, itmp(8), icor
  CHARACTER(100) :: ctmp
  !
  WRITE(*,*) 
  WRITE(*,*) "#####  Read Correlation Index File  #####" 
  WRITE(*,*) 
  !
  ALLOCATE(indx(3,nsite*nsite,8))
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
     READ(fi,*) itmp(1), itmp(2), itmp(3), itmp(4)
     IF(itmp(2) == 0 .AND. itmp(4) == 0) THEN ! Up-Up correlation
        ncor(1) = ncor(1) + 1
        indx(1:3,ncor(1),1) = (/itmp(1) + 1, itmp(3) + 1, icor/)
     ELSE IF (itmp(2) == 1 .AND. itmp(4) == 1) THEN ! Down-Down correlation
        ncor(2) = ncor(2) + 1
        indx(1:3,ncor(2),2) = (/itmp(1) + 1, itmp(3) + 1, icor/)
     END IF
  END DO
  !
  WRITE(*,*) "    Number of Up-Up Index : ", ncor(1)
  WRITE(*,*) "    Number of Down-Down Index : ", ncor(2)
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
     IF(itmp(1) /= itmp(3) .OR. itmp(5) /= itmp(7)) CYCLE
     !
     IF(itmp(2) == 0 .AND. itmp(4) == 0) THEN
        IF(itmp(6) == 0 .AND. itmp(8) == 0) THEN ! UpUpUpUp
           ncor(3) = ncor(3) + 1
           indx(1:3,ncor(3),3) = (/itmp(1) + 1, itmp(5) + 1, icor/)
        ELSE IF(itmp(6) == 1 .AND. itmp(8) == 1) THEN ! UpUpDownDown
           ncor(4) = ncor(4) + 1
           indx(1:3,ncor(4),4) = (/itmp(1) + 1, itmp(5) + 1, icor/)
        END IF
     ELSE IF(itmp(2) == 1 .AND. itmp(4) == 1) THEN
        IF(itmp(6) == 0 .AND. itmp(8) == 0) THEN ! DownDownUpUp
           ncor(5) = ncor(5) + 1
           indx(1:3,ncor(5),5) = (/itmp(1) + 1, itmp(5) + 1, icor/)
        ELSE IF(itmp(6) == 1 .AND. itmp(8) == 1) THEN ! DownDownDownDown
           ncor(6) = ncor(6) + 1
           indx(1:3,ncor(6),6) = (/itmp(1) + 1, itmp(5) + 1, icor/)
        END IF
     ELSE IF((itmp(2) == 0 .AND. itmp(4) == 1) .AND. (itmp(6) == 1 .AND. itmp(8) == 0)) THEN
        !
        ! Up-Down-Down-Up = S+S-
        !
        ncor(7) = ncor(7) + 1
        indx(1:3,ncor(7),7) = (/itmp(1) + 1, itmp(5) + 1, icor/)
     ELSE IF((itmp(2) == 1 .AND. itmp(4) == 0) .AND. (itmp(6) == 0 .AND. itmp(8) == 1)) THEN
        !
        ! Down-Up-Up-Down = S-S+
        !
        ncor(8) = ncor(8) + 1
        indx(1:3,ncor(8),8) = (/itmp(1) + 1, itmp(5) + 1, icor/)
     END IF
  END DO
  !
  CLOSE(fi)
  !
  WRITE(*,*) "    Number of UpUpUpUp         Index : ", ncor(3)
  WRITE(*,*) "    Number of UpUpDownDown     Index : ", ncor(4)
  WRITE(*,*) "    Number of DownDownUpUp     Index : ", ncor(5)
  WRITE(*,*) "    Number of DownDownDownDown Index : ", ncor(6)
  WRITE(*,*) "    Number of Plus-Minus       Index : ", ncor(7)
  WRITE(*,*) "    Number of Minus-Plus       Index : ", ncor(8)
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
  INTEGER :: fi = 10, icor, itmp(8), ii, iwfc
  COMPLEX(8),ALLOCATABLE :: cor0(:)
  REAL(8),ALLOCATABLE :: cor0_r(:,:)
  CHARACTER(256) :: filename, set, step
  !
  ALLOCATE(cor(nsite,nsite,6,nwfc))
  ALLOCATE(cor0(MAX(ncor1,ncor2)), cor0_r(2,MAX(ncor1,ncor2)))
  cor(1:nsite,1:nsite,1:6,1:nwfc) = CMPLX(0d0, 0d0, KIND(1d0))
  !
  ! One Body Correlation Function
  !
  DO iwfc = 1, nwfc
     !
     filename = TRIM(filehead) // "_cisajs" // TRIM(filetail(iwfc))
     OPEN(fi, file = TRIM(filename))
     IF(calctype == 4) THEN ! mVMC
        READ(fi,*) cor0_r(1:2, 1:ncor1)
     ELSE ! HPhi
        DO icor = 1, ncor1
           READ(fi,*) itmp(1:4), cor0_r(1:2, icor)
        END DO
     END IF ! IF(calctype == 4)
     CLOSE(fi)
     !
     DO icor = 1, ncor1
        cor0(icor) = CMPLX(cor0_r(1,icor), cor0_r(2,icor), KIND(1d0))
     END DO
     !
     ! Map it into Up-Up(1) and Down-Down(2) Correlation
     !
     DO ii = 1, 2
        DO icor = 1, ncor(ii)
           cor(indx(1,icor,ii), indx(2,icor,ii), ii, iwfc) = cor0(indx(3,icor,ii))
        END DO
     END DO
     !
     ! Two Body Correlation function
     !
     filename = TRIM(filehead) // "_cisajscktalt" // TRIM(filetail(iwfc))
     OPEN(fi, file = TRIM(filename))
     IF(calctype == 4) THEN ! mVMC
        READ(fi,*) cor0_r(1:2, 1:ncor2)
     ELSE ! HPhi
        DO icor = 1, ncor2
           READ(fi,*) itmp(1:8), cor0_r(1:2, icor)
        END DO
     END IF
     CLOSE(fi)
     !
     DO icor = 1, ncor2
        cor0(icor) = CMPLX(cor0_r(1,icor), cor0_r(2,icor), KIND(1d0))
     END DO
     !
     ! Map it into Density-Density(3), Sz-Sz(4), S+S-(5), S-S+(6) Correlation
     !
     DO ii = 3, 6, 3
        DO icor = 1, ncor(ii)
           cor(indx(1,icor,ii), indx(2,icor,ii), 3, iwfc) &
           &  = cor( indx(1,icor,ii), indx(2,icor,ii), 3, iwfc) &
           &  + cor0(indx(3,icor,ii))
           cor(indx(1,icor,ii), indx(2,icor,ii), 4, iwfc) &
           &  = cor( indx(1,icor,ii), indx(2,icor,ii), 4, iwfc) &
           &  + cor0(indx(3,icor,ii))
        END DO
     END DO
     DO ii = 4, 5
        DO icor = 1, ncor(ii)
           cor(indx(1,icor,ii), indx(2,icor,ii), 3, iwfc) &
           &  = cor( indx(1,icor,ii), indx(2,icor,ii), 3, iwfc) &
           &  + cor0(indx(3,icor,ii))
           cor(indx(1,icor,ii), indx(2,icor,ii), 4, iwfc) &
           &  = cor( indx(1,icor,ii), indx(2,icor,ii), 4, iwfc) &
           &  - cor0(indx(3,icor,ii))
        END DO
     END DO
     cor(1:nsite,1:nsite, 4, iwfc) = cor(1:nsite,1:nsite, 4, iwfc) * 0.25d0
     DO ii = 7, 8
        DO icor = 1, ncor(ii)
           cor(indx(1,icor,ii), indx(2,icor,ii), ii-2, iwfc) = cor0(indx(3,icor,ii))
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
           theta = - tpi * DOT_PRODUCT(kvec(1:2,ik), (site(1:2,isite) - site(1:2,jsite)))
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
  USE fourier_val, ONLY : cor_k, nk, nktot, nk_row, kvec, kvec_tot, &
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
  WRITE(*,*) "  Correlation in k-space : ", "output/" // TRIM(filehead) // "_corr", "*.dat"
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
     WRITE(fo,*) "# k1[1] k2[2](Cart.) UpUp[3,4,15,16] (Re. Im. Err.) DownDown[5,6,17,18]"
     WRITE(fo,*) "# Density[7:8,19,20] SzSz[9,10,21,22] S+S-[11,12,23,24] S-S+[13,14,25.26]"
     !
     DO ik = 1, nk
        WRITE(fo,'(26e15.5)') tpi * MATMUL(recipr(1:2,1:2), kvec(1:2,ik)), &
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
        WRITE(fo,*) "# k1[1] k2[2](Cart.) UpUp[3,4] (Re. Im.) DownDown[5,6]"
        WRITE(fo,*) "# Density[7:8] SzSz[9,10] S+S-[11,12] S-S+[13,14]"
        !
        DO ik = 1, nk
           WRITE(fo,'(14e15.5)') tpi * MATMUL(recipr(1:2,1:2), kvec(1:2,ik)), cor_k(ik,1:6,iwfc)
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
  DO idim = 1, 2
     WRITE(fo,'(3e15.5)') tpi * recipr(1:2,idim)
  END DO
  !
  DO ik = 1, nktot
     WRITE(fo,'(2e15.5,i7)') tpi * MATMUL(recipr(1:2,1:2), kvec_tot(1:2,ik)), equiv(ik)
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
