MODULE fourier_val
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & interval, &
  & numave, &
  & nkg(3), & ! k-grid for momentum ditribution
  & nk_line, & ! Numberof along each k line
  & nnode, & ! Number of node of k-path
  & calctype, & ! (0) Lanczos (1) TPQ (2) FullDiag (3) LOBCG (4) mVMC
  & nwfc, & ! Number of state
  & nr,   & ! Number of R-vector
  & norb, & ! Number of orbitals per unit cell
  & box(3,3), & ! Supercell index
  & nsite,   &  ! Number of sites
  & ncor1,   &  ! Nomber of One-body Correlation function
  & ncor2,   &  ! Number of Two-body Correlation function
  & ncor(8), &  ! Number of Correlation function for each index(See below)
  & nk          ! Number of k to be computed
  !
  REAL(8),SAVE :: &
  & recipr(3,3)    ! Reciprocal lattice vector
  !
  CHARACTER(256),SAVE :: &
  & filehead, & ! Filename header for correlation functions
  & file_one, & ! Filename for One-body Correlation
  & file_two    ! Filename for Two-body Correlation
  !
  INTEGER,ALLOCATABLE,SAVE :: &
  & nreq(:),    & ! (nr) Number of equivalent R-vector for each R
  & irv(:,:,:), & ! (3,125,nr) R-vector
  & rindx(:),   & ! (nsite) Index of R
  & orb(:),     & ! (nsite) Index of orbital
  & indx(:,:,:,:) ! (nr,8,norb,norb) Mapping index for each Correlation function
  !
  REAL(8),ALLOCATABLE,SAVE :: &
  & knode(:,:), & ! (3,nnode) Nodes of k path
  & phase(:,:), & ! (125,nr) Boundary phase 
  & kvec(:,:)     ! (3,nk) k-vector in the 1st BZ
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & cor(:,:,:,:,:), & ! (nr,6,norb,norb,nwfc) Correlation function in real space (See below)
  & cor_k(:,:,:,:,:)  ! (nk,6,norb,norb,nwfc) Correlation function in the k-space (See below)
  !
  CHARACTER(256),ALLOCATABLE :: &
  & kname(:), & ! (nnode) Label of k-point node
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
CONTAINS
!
SUBROUTINE key2lower(key)
  CHARACTER(*) :: key
  !
  INTEGER :: ii, acode
  !
  DO ii = 1, LEN(TRIM(key))
     acode = IACHAR(key(ii:ii))
     IF(65 <= acode .AND. acode <= 90) THEN
        key(ii:ii) = ACHAR(acode + 32)
     END IF
  END DO
END SUBROUTINE key2lower
!
! Read from HPhi/mVMC input files
!
SUBROUTINE read_filename()
  !
  USE fourier_val, ONLY : file_one, file_two, filehead, nsite, nwfc, &
  &                       filetail, calctype, numave, interval
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, lanczos_max, irun, istep, iwfc, idx_start
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
     DO istep = 0, lanczos_max - 1
        DO irun = 0, numave - 1
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
     DO iwfc = 1, nwfc
        WRITE(filetail(iwfc),'(a,i3.3,a)') "_", iwfc - 1 + idx_start, ".dat"
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
  USE fourier_val, ONLY : recipr, box, nsite, phase, irv, rindx, orb, &
  &                       nr, nreq, norb, nnode, knode, nk_line, kname, nkg
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, isite, ii, ir, ipiv(3), irv0(3), i1, i2, i3, inode
  REAL(8) :: phase0(3), work(10), direct(3,3), rrv(3), lenrv, lenrv0
  CHARACTER(256) :: filename
  INTEGER,ALLOCATABLE :: irv1(:,:)
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
  READ(fi,*) phase0(1:3)  
  WRITE(*,*) "    Boundary phase[degree] : "
  WRITE(*,'(4x3f15.10)') phase0(1:3)
  phase0(1:3) = phase0(1:3) * ACOS(-1.0d0) / 180.0d0
  !
  ! Supercell index (a0w, a0l, a1w, a1l)
  !
  DO ii = 1, 3
     READ(fi,*) box(1:3,ii)
  END DO
  WRITE(*,*) "    Supercell Index :"
  WRITE(*,'(3i8)') box(1:3, 1:3)
  !
  ! R-vector and orbital index
  !
  ALLOCATE(irv1(3,nsite), orb(nsite), rindx(nsite))
  !
  nr = 0
  DO isite = 1, nsite
     READ(fi,*) irv0(1:3), orb(isite)
     DO ir = 1, nr
        IF(ALL(irv1(1:3,ir) == irv0(1:3))) THEN
           rindx(isite) = ir
           GOTO 10
        END IF
     END DO
     nr = nr + 1
     irv1(1:3, nr) = irv0(1:3)
     rindx(isite) = nr
     !
10   CONTINUE
     !
  END DO
  orb(1:nsite) = orb(1:nsite) + 1
  norb = MAXVAL(orb)
  WRITE(*,*) "    Number of orbitals :", norb
  !
  ! k-point
  !
  READ(fi,*) nnode, nk_line
  ALLOCATE(knode(3,nnode), kname(nnode))
  WRITE(*,*) "    Number of k-node, and k-points along lines :", nnode, nk_line
  WRITE(*,*) "      k-node :"
  DO inode = 1, nnode
     READ(fi,*) kname(inode), knode(1:3,inode)
     WRITE(*,'(a,a,3f10.5)') "      ", TRIM(kname(inode)), knode(1:3,inode)
  END DO
  READ(fi,*) nkg(1:3)
  WRITE(*,'(a,3i3)') "k-grid for momentum distribution :", nkg(1:3)
  !
  CLOSE(fi)
  !
  ! Compute Reciprocal Lattice Vector
  !
  recipr(1:3,1:3) = transpose(direct(1:3,1:3))
  CALL dgetrf(3, 3, recipr, 3, Ipiv, ii)
  CALL dgetri(3, recipr, 3, ipiv, work, 10, ii)
  WRITE(*,*) "    Reciplocal lattice vector :"
  WRITE(*,'(4x3f15.10)') recipr(1:3, 1:3)
  !
  ! Move original R-vector to the nearest one with periodic boundary cond.
  !
  ALLOCATE(nreq(nsite), irv(3,125,nsite), phase(125,nsite))
  !
  WRITE(*,*) "    Number of R-vector :", nr
  DO ir = 1, nr
     lenrv0 = 1.0d10
     DO i1 = -2, 2
        DO i2 = -2, 2
           DO i3 = -2, 2
              !
              irv0(1:3) = irv1(1:3,ir) + MATMUL(box(1:3,1:3), (/i1,i2,i3/))
              rrv(1:3) = MATMUL(direct(1:3,1:3), DBLE(irv0(1:3)))
              lenrv = SQRT(DOT_PRODUCT(rrv, rrv))
              IF(lenrv < lenrv0 - 1.0d-6) THEN
                 lenrv0 = lenrv
                 nreq(ir) = 1
              ELSE IF(ABS(lenrv - lenrv0) < 1.0d-6) THEN
                 nreq(ir) = nreq(ir) + 1
              ELSE
                 CYCLE
              END IF
              !
              irv(1:3, nreq(ir), ir) = irv0(1:3)
              phase(nreq(ir), ir) = DOT_PRODUCT(DBLE((/i1,i2,i3/)), phase0(1:3))
              !
           END DO ! i3 = -2, 2
        END DO ! i2 = -2, 2
     END DO ! i1 = -2, 2
     !
     DO i1 = 1, nreq(ir)
        WRITE(*,'(3i5,f7.2,a)',advance="no") irv(1:3, i1, ir), phase(i1, ir), ", "
     END DO
     WRITE(*,*)
     !
  END DO ! ir = 1, nr
  !
  DEALLOCATE(irv1)
  !
END SUBROUTINE read_geometry
!
! Set k points
!
SUBROUTINE set_kpoints()
  !
  USE fourier_val, ONLY : nk, kvec, nnode, nk_line, knode, nkg
  !
  IMPLICIT NONE
  !
  INTEGER :: inode, ik, i1, i2, i3
  REAL(8) :: xx
  !
  nk = nk_line * (nnode - 1) + 1 + PRODUCT(nkg(1:3))
  WRITE(*,*) "     Number of k : ", nk
  ALLOCATE(kvec(3,nk))
  !
  kvec(1:3,1) = knode(1:3,1)
  nk = 1
  DO inode = 1, nnode - 1
     DO ik = 1, nk_line
        xx = DBLE(ik) / DBLE(nk_line)
        nk = nk + 1
        kvec(1:3,nk) = (1d0 - xx)*knode(1:3,inode) + xx*knode(1:3,inode+1)
     END DO
  END DO
  !
  ! k-grid for momentum distribution
  !
  DO i1 = 1, nkg(1)
     DO i2 = 1, nkg(2)
        DO i3 = 1, nkg(3)
           nk = nk + 1
           kvec(1:3,nk) = DBLE((/i1,i2,i3/)-1) / DBLE(nkg(1:3))
        END DO
     END DO
  END DO
  !
END SUBROUTINE set_kpoints
!
! Read Correlation Function
!
SUBROUTINE read_corrindx()
  !
  USE fourier_val, ONLY : file_one, file_two, ncor1, ncor2, ncor, indx, &
  &                       calctype, nr, rindx, orb, norb, irv
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, itmp(8), icor, ir, iorb, jorb
  CHARACTER(100) :: ctmp
  !
  WRITE(*,*) 
  WRITE(*,*) "#####  Read Correlation Index File  #####" 
  WRITE(*,*) 
  !
  ALLOCATE(indx(nr,8,norb,norb))
  indx(1:nr,1:8,1:norb,1:norb) = 0
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
     !
     READ(fi,*) itmp(1:4)
     !
     IF(ANY(irv(1:3,1,rindx(itmp(1)+1)) /= 0)) THEN
        WRITE(*,'(a,i0,a)') "       REMARK : The left operator at # ", icor, " is not at R=0."
        CYCLE
     END IF
     !
     IF(itmp(2) == 0 .AND. itmp(4) == 0) THEN
        !
        ! Up-Up correlation
        !
        indx(rindx(itmp(3) + 1), 1, orb(itmp(1) + 1), orb(itmp(3) + 1)) = icor 
     ELSE IF (itmp(2) == 1 .AND. itmp(4) == 1) THEN
        !
        ! Down-Down correlation
        !
        indx(rindx(itmp(3) + 1), 2, orb(itmp(1) + 1), orb(itmp(3) + 1)) = icor 
     END IF
  END DO
  !
  WRITE(*,*) "    Number of Up-Up Index : ",     COUNT(indx(1:nr, 1, 1:norb, 1:norb) /= 0)
  WRITE(*,*) "    Number of Down-Down Index : ", COUNT(indx(1:nr, 2, 1:norb, 1:norb) /= 0)
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
        IF(ANY(irv(1:3,1,rindx(itmp(1)+1)) /= 0)) THEN
           WRITE(*,'(a,i0,a)') "       REMARK : The left operator at # ", icor, " is not at R=0."
           CYCLE
        END IF
        !
        IF(itmp(2) == 0 .AND. itmp(4) == 0) THEN
           !
           IF(itmp(6) == 0 .AND. itmp(8) == 0) THEN
              !
              ! UpUpUpUp
              !
              indx(rindx(itmp(5) + 1), 3, orb(itmp(1) + 1), orb(itmp(5) + 1)) = icor 
           ELSE IF(itmp(6) == 1 .AND. itmp(8) == 1) THEN
              !
              ! UpUpDownDown
              !
              indx(rindx(itmp(5) + 1), 4, orb(itmp(1) + 1), orb(itmp(5) + 1)) = icor 
           END IF
           !
        ELSE IF(itmp(2) == 1 .AND. itmp(4) == 1) THEN
           !
           IF(itmp(6) == 0 .AND. itmp(8) == 0) THEN
              !
              ! DownDownUpUp
              !
              indx(rindx(itmp(5) + 1), 5, orb(itmp(1) + 1), orb(itmp(5) + 1)) = icor 
           ELSE IF(itmp(6) == 1 .AND. itmp(8) == 1) THEN
              !
              ! DownDownDownDown
              !
              indx(rindx(itmp(5) + 1), 6, orb(itmp(1) + 1), orb(itmp(5) + 1)) = icor 
           END IF
           !
        ELSE IF(calctype /= 4) THEN
           !
           ! S+S- & S-S+ for HPhi
           !
           IF((itmp(2) == 0 .AND. itmp(4) == 1) .AND. (itmp(6) == 1 .AND. itmp(8) == 0)) THEN
              !
              ! Up-Down-Down-Up = S+S-
              !
              indx(rindx(itmp(5) + 1), 7, orb(itmp(1) + 1), orb(itmp(5) + 1)) = icor 
           ELSE IF((itmp(2) == 1 .AND. itmp(4) == 0) .AND. (itmp(6) == 0 .AND. itmp(8) == 1)) THEN
              !
              ! Down-Up-Up-Down = S-S+
              !
              indx(rindx(itmp(5) + 1), 8, orb(itmp(1) + 1), orb(itmp(5) + 1)) = icor 
           END IF
           !
        END IF ! (calctype /= 4)
        !
     END IF
     !
     IF(calctype == 4 .AND. (itmp(1) == itmp(7) .AND. itmp(3) == itmp(5))) THEN
        !
        IF(ANY(irv(1:3,1,rindx(itmp(1)+1)) /= 0)) THEN
           WRITE(*,'(a,i0,a)') "       REMARK : The left operator at # ", icor, " is not at R=0."
           CYCLE
        END IF
        !
        ! S+S- & S-S+ for mVMC
        !
        IF((itmp(2) == 0 .AND. itmp(4) == 0) .AND. (itmp(6) == 1 .AND. itmp(8) == 1)) THEN
           !
           ! Up-Down-Down-Up = S+S-
           !
           indx(rindx(itmp(5) + 1), 7, orb(itmp(1) + 1), orb(itmp(5) + 1)) = icor 
        ELSE IF((itmp(2) == 1 .AND. itmp(4) == 1) .AND. (itmp(6) == 0 .AND. itmp(8) == 0)) THEN
           !
           ! Down-Up-Up-Down = S-S+
           !
           indx(rindx(itmp(5) + 1), 8, orb(itmp(1) + 1), orb(itmp(5) + 1)) = icor 
        END IF
        !
     END IF ! (calctype == 4 .AND. (itmp(1) == itmp(7) .AND. itmp(3) == itmp(5)))
     !
  END DO
  !
  IF(COUNT(indx(1:nr,3:8,1:norb,1:norb) == 0) /= 0) THEN
     WRITE(*,*) "ERROR! The following correlation function is missed:"
     WRITE(*,*) "R,   kind,   orb1,   orb2"
     DO icor = 3, 8
        DO ir = 1, nr
           DO iorb = 1, norb
              DO jorb = 1, norb
                 IF(indx(ir,icor,iorb,jorb) == 0) THEN
                    IF(icor == 3) THEN
                       WRITE(*,'(i0,a,i0,2x,i0)') ir, " Up-Up-Up-Up         ", iorb, jorb
                    ELSE IF(icor == 4) THEN
                       WRITE(*,'(i0,a,i0,2x,i0)') ir, " Up-Up-Down-Down     ", iorb, jorb
                    ELSE IF(icor == 5) THEN
                       WRITE(*,'(i0,a,i0,2x,i0)') ir, " Down-Down-Up-Up     ", iorb, jorb
                    ELSE IF(icor == 6) THEN
                       WRITE(*,'(i0,a,i0,2x,i0)') ir, " Down-Down-Down-Down ", iorb, jorb
                    ELSE IF(icor == 7) THEN
                       WRITE(*,'(i0,a,i0,2x,i0)') ir, " Up-Down-Down-Up     ", iorb, jorb
                    ELSE
                       WRITE(*,'(i0,a,i0,2x,i0)') ir, " Down-Up-Up-Down     ", iorb, jorb
                    END IF
                 END IF
              END DO
           END DO
        END DO
     END DO
     STOP "Missing indices for the Green function."
  END IF
  !
  WRITE(*,*) "    Number of UpUpUpUp         Index : ", COUNT(indx(1:nr, 3, 1:norb, 1:norb) /= 0)
  WRITE(*,*) "    Number of UpUpDownDown     Index : ", COUNT(indx(1:nr, 4, 1:norb, 1:norb) /= 0)
  WRITE(*,*) "    Number of DownDownUpUp     Index : ", COUNT(indx(1:nr, 5, 1:norb, 1:norb) /= 0)
  WRITE(*,*) "    Number of DownDownDownDown Index : ", COUNT(indx(1:nr, 6, 1:norb, 1:norb) /= 0)
  WRITE(*,*) "    Number of Plus-Minus       Index : ", COUNT(indx(1:nr, 7, 1:norb, 1:norb) /= 0)
  WRITE(*,*) "    Number of Minus-Plus       Index : ", COUNT(indx(1:nr, 8, 1:norb, 1:norb) /= 0)
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
  &                       ncor1, ncor2, indx, cor, norb, nr, irv
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, icor, itmp(8), iwfc, iorb, jorb, ir, ir0
  COMPLEX(8),ALLOCATABLE :: cor0(:)
  REAL(8) :: cor0_r(2)
  CHARACTER(256) :: filename
  !
  ALLOCATE(cor(nr,6,norb,norb,nwfc))
  ALLOCATE(cor0(0:MAX(ncor1,ncor2)))
  cor(1:nr,1:6,1:norb,1:norb,1:nwfc) = CMPLX(0d0, 0d0, KIND(1d0))
  cor0(0) = CMPLX(0d0, 0d0, KIND(1d0))
  !
  DO ir = 1, nr
     IF(all(irv(1:3, 1, ir) == 0)) THEN
        ir0 = ir
        EXIT
     END IF
  END DO
  !
  DO iwfc = 1, nwfc
     !
     ! One Body Correlation Function
     !
     filename = TRIM(filehead) // "_cisajs" // TRIM(filetail(iwfc))
     OPEN(fi, file = TRIM(filename))
     !
     DO icor = 1, ncor1
        READ(fi,*) itmp(1:4), cor0_r(1:2)
        cor0(icor) = CMPLX(cor0_r(1), cor0_r(2), KIND(1d0))
     END DO
     !
     CLOSE(fi)
     !
     ! Map it into Up-Up(1) and Down-Down(2) Correlation
     !
     DO iorb = 1, norb
        DO jorb = 1, norb
           DO ir = 1, nr
              cor(ir, 1:2, jorb, iorb, iwfc) = cor0(indx(ir, 1:2, jorb, iorb))
           END DO
        END DO
     END DO
     !
     ! Two Body Correlation function
     !
     filename = TRIM(filehead) // "_cisajscktalt" // TRIM(filetail(iwfc))
     OPEN(fi, file = TRIM(filename))
     !
     DO icor = 1, ncor2
        READ(fi,*) itmp(1:8), cor0_r(1:2)
        cor0(icor) = CMPLX(cor0_r(1), cor0_r(2), KIND(1d0))
     END DO
     !
     CLOSE(fi)
     !
     ! Map it into Density-Density(3), Sz-Sz(4), S+S-(5), S-S+(6) Correlation
     !
     ! Up-Up-Up-Up and Down-Down-Down-Down into Density-Density & Sz-Sz
     !
     DO iorb = 1, norb
        DO jorb = 1, norb
           DO ir = 1, nr
              !
              cor(ir, 3, jorb, iorb, iwfc) = cor0(indx(ir, 3, jorb, iorb)) &
              &                            + cor0(indx(ir, 4, jorb, iorb)) &
              &                            + cor0(indx(ir, 5, jorb, iorb)) &
              &                            + cor0(indx(ir, 6, jorb, iorb))
              !
              cor(ir, 3, jorb, iorb, iwfc) = cor(ir, 3,   jorb, iorb, iwfc) &
              &                       - SUM(cor(ir0, 1:2, iorb, iorb, iwfc)) &
              &                       * SUM(cor(ir0, 1:2, jorb, jorb, iwfc))
              !
              cor(ir, 4, jorb, iorb, iwfc) = cor0(indx(ir, 3, jorb, iorb)) &
              &                            - cor0(indx(ir, 4, jorb, iorb)) &
              &                            - cor0(indx(ir, 5, jorb, iorb)) &
              &                            + cor0(indx(ir, 6, jorb, iorb))
              !
              cor(ir, 4, jorb, iorb, iwfc) = cor(ir, 4, jorb, iorb, iwfc) * 0.25d0
              !
              ! Up-Down-Down-Up(S+S-) and Down-Up-Up-Down(S-S+)
              !
              cor(ir, 5:6, jorb, iorb, iwfc) = cor0(indx(ir, 7:8, jorb, iorb))
              !
           END DO ! ir = 1, nr
        END DO ! jorb = 1, norb
     END DO ! iorb = 1, norb
     !
     ! For mVMC
     !
     !   Ciu+ Cid Cjd+ Cju = delta_{ij} Ciu+ Ciu - Ciu+ Cju Cjd+ Cid
     !   Cid+ Ciu Cju+ Cjd = delta_{ij} Cid+ Cid - Cid+ Cjd Cju+ Ciu
     !
     IF (calctype == 4) THEN
        cor(1:nr, 5:6, 1:norb, 1:norb, iwfc) = - cor(1:nr, 5:6, 1:norb, 1:norb, iwfc)
        DO iorb = 1, norb
           cor(1, 5:6, iorb, iorb, iwfc) = cor(1, 5:6, iorb, iorb, iwfc) &
           &                             + cor(1, 1:2, iorb, iorb, iwfc)
        END DO
     END IF
     !
     ! S.S = Sz Sz + 0.5 * (S+S- + S-S+)
     !
     cor(1:nr, 6, 1:norb, 1:norb, iwfc) = cor(1:nr, 4, 1:norb, 1:norb, iwfc) &
     &                        + 0.5d0 * ( cor(1:nr, 5, 1:norb, 1:norb, iwfc) &
     &                                  + cor(1:nr, 6, 1:norb, 1:norb, iwfc) )
     !
  END DO ! iwfc = 1, nwfc
  !
  DEALLOCATE(cor0, indx)
  !
END SUBROUTINE read_corrfile
!
! Fourier transformation
!
SUBROUTINE fourier_cor()
  !
  USE fourier_val, ONLY : cor, cor_k, kvec, nwfc, nk, nr, nreq, norb, irv, phase
  IMPLICIT NONE
  !
  INTEGER :: ik, ir, ireq
  REAL(8) :: tpi = 2.0 * ACOS(-1d0), theta
  COMPLEX(8),ALLOCATABLE :: fmat(:,:)
  !
  ALLOCATE(fmat(nk,nr), cor_k(nk,6,norb,norb,nwfc))
  !
  ! Matirx for Fourier trans. exp(-i k R)
  !
  DO ik = 1, nk
     DO ir = 1, nr
        fmat(ik,ir) = CMPLX(0d0, 0d0, KIND(1d0))
        DO ireq = 1, nreq(ir)
           theta = - tpi * DOT_PRODUCT(kvec(1:3,ik), DBLE(irv(1:3,ireq,ir))) &
           &     + phase(ireq,ir)
           fmat(ik,ir) = fmat(ik,ir) + CMPLX(COS(theta), SIN(theta), KIND(1d0))
        END DO
        fmat(ik,ir) = fmat(ik,ir) / DBLE(nreq(ir))
     END DO ! ir = 1, nr
  END DO ! ik = 1, nk
  !
  CALL zgemm('N', 'N', nk, 6*norb*norb*nwfc, nr, CMPLX(1d0, 0d0, KIND(1d0)), fmat, nk, &
  &          cor, nr, CMPLX(0d0,0d0,KIND(1d0)), cor_k, nk)
  !
  cor_k(1:nk,1:2,1:norb,1:norb,1:nwfc) = cor_k(1:nk,1:2,1:norb,1:norb,1:nwfc)
  cor_k(1:nk,3:6,1:norb,1:norb,1:nwfc) = cor_k(1:nk,3:6,1:norb,1:norb,1:nwfc) / dble(nr)
  !
  DEALLOCATE(fmat, cor)
  !
END SUBROUTINE fourier_cor
!
! Output Fourier component of Correlation function
!
SUBROUTINE output_cor()
  !
  USE fourier_val, ONLY : cor_k, nk, nnode, knode, nk_line, kname, norb, interval, &
  &                       nwfc, recipr, filehead, filetail, calctype, nkg, numave
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, ik, iwfc, inode, iorb, jorb, ii, ikk, iwfc1, iwfc2, istep
  REAL(8) :: dk(3), dk_cart(3), xk(nk), &
  &          xk_label(nnode), klength
  CHARACTER(256) :: filename
  COMPLEX(8),ALLOCATABLE :: cor_ave(:,:,:,:), cor_err(:,:,:,:)
  !
  ! Compute x-position for plotting band
  !
  xk(1) = 0.0
  ikk = 1
  DO inode = 1, nnode - 1
     dk(1:3) = knode(1:3, inode+1) - knode(1:3, inode)
     dk_cart(1:3) = MATMUL(recipr(1:3,1:3), dk(1:3))
     klength = SQRT(DOT_PRODUCT(dk_cart, dk_cart)) / DBLE(nk_line)
     xk_label(inode) = xk(ikk)
     DO ik = 1, nk_line
        xk(ikk+1) = xk(ikk) + klength
        ikk = ikk + 1
     END DO
  END DO
  xk_label(nnode) = xk(ikk)
  !
  ! Output Correlation function in the 1st BZ
  !
  WRITE(*,*) 
  WRITE(*,*) "#####  Output Files  #####" 
  WRITE(*,*) 
  !
  WRITE(*,*) "  Correlation in k-space : ", TRIM(filehead) // "_corr", "*.dat"
  !
  IF(calctype == 1 .OR. calctype == 4) THEN
     !
     ! TPQ/mVMC
     !
     ALLOCATE(cor_ave(ikk,6,norb,norb), cor_err(ikk,6,norb,norb))
     !
     DO istep = 1, nwfc / numave
        !
        iwfc1 = numave*(istep-1) + 1
        iwfc2 = numave*istep
        !
        ! Average
        !
        cor_ave(1:ikk,1:6,1:norb,1:norb) = SUM(cor_k(1:ikk,1:6,1:norb,1:norb,&
        &                                      iwfc1:iwfc2), 5) / DBLE(numave)
        !
        ! Variance
        !
        cor_err(1:ikk,1:6,1:norb,1:norb) = 0d0
        DO iwfc = iwfc1, iwfc2
           cor_err(1:ikk,1:6,1:norb,1:norb) = cor_err(1:ikk,1:6,1:norb,1:norb) &
           & + CMPLX( DBLE(cor_k(1:ikk,1:6,1:norb,1:norb,iwfc) - cor_ave(1:ikk,1:6,1:norb,1:norb))**2, &
           &         AIMAG(cor_k(1:ikk,1:6,1:norb,1:norb,iwfc) - cor_ave(1:ikk,1:6,1:norb,1:norb))**2, &
           &         KIND(0d0))
        END DO
        !
        ! Standard Error
        !
        IF(numave == 1) THEN
           cor_err(1:ikk,1:6,1:norb,1:norb) = CMPLX(0d0, 0d0, KIND(0d0))
        ELSE
           cor_err(1:ikk,1:6,1:norb,1:norb) = CMPLX(SQRT( DBLE(cor_err(1:ikk,1:6,1:norb,1:norb))), &
           &                                       SQRT(AIMAG(cor_err(1:ikk,1:6,1:norb,1:norb))), KIND(0d0)) &
           &                               / SQRT(DBLE(numave * (numave - 1)))
        END IF
        !
        IF(calctype == 1)THEN
           WRITE(filename,'(a,a,i0,a)') TRIM(filehead), "_corr_step", interval*(istep-1), ".dat"
        ELSE 
           filename = TRIM(filehead) // "_corr.dat"
        END IF
        OPEN(fo, file = TRIM(filename))
        !
        WRITE(fo,*) "# k-length[1]"
        ii = 1
        DO iorb = 1, norb
           DO jorb = 1, norb
              WRITE(fo,'(a,i3,a,i3)') "# Orbital", iorb, " to Orbital", jorb
              WRITE(fo,'(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)') &
              & "#  UpUp[", ii+1, ",", ii+2, ",", ii+13, ",", ii+14, &
              & "] (Re. Im. Err.) DownDown[", ii+3, ",", ii+4, ",", ii+15, ",", ii+16, "]"
              WRITE(fo,'(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)') &
              & "#  Density[", ii+5, ",", ii+6, ",", ii+17, ",", ii+18, &
              & "] SzSz[", ii+7, ",", ii+8, ",", ii+19, ",", ii+20, &
              & "] S+S-[", ii+9, ",", ii+10, ",", ii+21, ",", ii+22, &
              & "] S.S[", ii+11, ",", ii+12, ",", ii+23, ",", ii+24, "]"
              ii = ii+24
           END DO
        END DO
        !
        DO ik = 1, ikk
           WRITE(fo,'(e15.5)',advance="no") xk(ik)
           DO iorb = 1, norb
              DO jorb = 1, norb
                 WRITE(fo,'(24e15.5)',advance="no") cor_ave(ik,1:6, jorb, iorb), cor_err(ik,1:6, jorb, iorb)
              END DO
           END DO
           WRITE(fo,*)
        END DO
        !
        CLOSE(fo)
        !
     END DO ! istep = 1, nwfc / numave
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
        WRITE(fo,*) "# k-length[1]"
        ii = 1
        DO iorb = 1, norb
           DO jorb = 1, norb
              WRITE(fo,'(a,i3,a,i3)') "# Orbital", iorb, " to Orbital", jorb
              WRITE(fo,'(a,i4,a,i4,a,i4,a,i4,a)') &
              & "#  UpUp[", ii+1, ",", ii+2, "] (Re. Im.) DownDown[", ii+3, ",", ii+4, "]"
              WRITE(fo,'(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)') &
              & "#  Density[", ii+5, ",", ii+6, "] SzSz[", ii+7, ",", ii+8, &
              & "] S+S-[", ii+9, ",", ii+10, "] S.S[", ii+11, ",", ii+12, "]"
              ii = ii+12
           END DO
        END DO
        !
        DO ik = 1, ikk
           WRITE(fo,'(1000e15.5)') xk(ik), cor_k(ik, 1:6, 1:norb, 1:norb, iwfc)
        END DO
        !
        CLOSE(fo)
        !
     END DO
     !
  END IF ! IF(calctype == 4)
  !
  OPEN(fo, file = "kpath.gp")
  !
  WRITE(fo,'(a)',advance="no") "set xtics ("
  DO inode = 1, nnode - 1
     WRITE(fo,'(a,a,a,f10.5,a)',advance="no") "'", TRIM(kname(inode)), "'  ", xk_label(inode), ", "
  END DO
  WRITE(fo,'(a,a,a,f10.5,a)') "'", TRIM(kname(nnode)), "' ", xk_label(nnode), ")"
  WRITE(fo,'(a)') "set ylabel 'Correlation function'"
  WRITE(fo,'(a)') "set grid xtics lt 1 lc 0"
  !
  CLOSE(fo)
  !
  ! FermiSuerfer file
  !
  DO iwfc = 1, nwfc
     !
     filename = TRIM(filehead) // "_corr" // TRIM(filetail(iwfc)) // ".frmsf"
     OPEN(fo, file = TRIM(filename))
     !
     WRITE(fo,*) nkg(1:3)
     WRITE(fo,*) 1
     WRITE(fo,*) norb
     WRITE(fo,*) REAL(recipr(1:3,1))
     WRITE(fo,*) REAL(recipr(1:3,2))
     WRITE(fo,*) REAL(recipr(1:3,3))
     DO iorb = 1, norb
        DO ik = ikk+1, nk
           WRITE(fo,*) SUM(REAL(cor_k(ik, 1:2, iorb, iorb, iwfc)))
        END DO
     END DO
     DO iorb = 1, norb
        DO ik = ikk+1, nk
           WRITE(fo,*) REAL(iorb) 
        END DO
     END DO
     CLOSE(fo)
     !
  END DO ! iwfc = 1, nwfc
  !
  DEALLOCATE(cor_k)
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
