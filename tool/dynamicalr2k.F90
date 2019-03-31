MODULE fourier_val
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & nomega, &
  & nk_line, & ! Numberof along each k line
  & nnode, & ! Number of node of k-path
  & nr,   & ! Number of R-vector
  & norb, & ! Number of orbitals per unit cell
  & box(3,3), & ! Supercell index
  & nsite,   &  ! Number of sites
  & ncor,   &  ! Nomber of Correlation function
  & nk          ! Number of k to be computed
  !
  REAL(8),SAVE :: &
  & omegamin, &
  & omegamax, &
  & recipr(3,3)    ! Reciprocal lattice vector
  !
  CHARACTER(256),SAVE :: &
  & filehead, & ! Filename header for correlation functions
  & file_gindx    ! Filename for index of Correlation
  !
  INTEGER,ALLOCATABLE,SAVE :: &
  & nreq(:),    & ! (nr) Number of equivalent R-vector for each R
  & irv(:,:,:), & ! (3,125,nr) R-vector
  & rindx(:),   & ! (nsite) Index of R
  & orb(:),     & ! (nsite) Index of orbital
  & indx(:,:) ! (nr,norb) Mapping index for each Correlation function
  !
  REAL(8),ALLOCATABLE,SAVE :: &
  & knode(:,:), & ! (3,nnode) Nodes of k path
  & phase(:,:), & ! (125,nr) Boundary phase 
  & kvec(:,:)     ! (3,nk) k-vector in the 1st BZ
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & cor(:,:,:), & ! (nr,norb,nomega) Correlation function in real space
  & cor_k(:,:,:)  ! (nk,norb,nomega) Correlation function in the k-space
  !
  CHARACTER(256),ALLOCATABLE :: &
  & kname(:) ! (nnode) Label of k-point node
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
  USE fourier_val, ONLY : file_gindx, filehead, nsite, omegamin, omegamax, nomega
  IMPLICIT NONE
  !
  INTEGER :: fi = 10
  CHARACTER(256) :: modpara, keyname, namelist
  REAL(8) :: eig0
  !
  WRITE(*,*) 
  WRITE(*,*) "#####  Read HPhi Input Files  #####" 
  WRITE(*,*) 
  !
  namelist = ""
  CALL GETARG(1, namelist)
  !
  ! Read from NameList file
  !
  OPEN(fi,file = TRIM(namelist))
  !
  DO
     READ(fi,*,END=10) keyname
     BACKSPACE(fi)
     CALL key2lower(keyname)
     !
     IF(TRIM(ADJUSTL(keyname)) == "singleexcitation") THEN
        READ(fi,*) keyname, file_gindx
     ELSE IF(TRIM(ADJUSTL(keyname)) == "pairexcitation") THEN
        READ(fi,*) keyname, file_gindx
     ELSE IF(TRIM(ADJUSTL(keyname)) == "modpara") THEN
        READ(fi,*) keyname, modpara
     ELSE
        READ(fi,*) keyname
     END IF
  END DO
  !
10 CONTINUE
  WRITE(*,*) "  Read from ", TRIM(namelist)
  CLOSE(FI)
  !
  WRITE(*,*) "    Excitation Index file : ", TRIM(ADJUSTL(file_gindx))
  WRITE(*,*) "    ModPara file : ", TRIM(ADJUSTL(modpara)) 
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
     ELSE IF(TRIM(ADJUSTL(keyname)) == "nomega") THEN
        READ(fi,*) keyname, nomega
     ELSE IF(TRIM(ADJUSTL(keyname)) == "omegamax") THEN
        READ(fi,*) keyname, omegamax
     ELSE IF(TRIM(ADJUSTL(keyname)) == "omegamin") THEN
        READ(fi,*) keyname, omegamin
     ELSE
        READ(fi,*) keyname
     END IF
  END DO
  !
20 CONTINUE
  WRITE(*,*) "  Read from ", TRIM(modpara)
  WRITE(*,*) "    FileHead : ", TRIM(ADJUSTL(filehead))
  WRITE(*,*) "    Number of site : ", nsite
  WRITE(*,*) "    Number of omega : ", nomega
  WRITE(*,*) "    Minimum Omega : ", omegamin
  WRITE(*,*) "    Maximum Omega : ", omegamax
  CLOSE(FI)
  !
  filehead = "output/" // TRIM(ADJUSTL(filehead))
  !
  OPEN(fi,file = TRIM(filehead)//"_energy.dat")
  READ(fi,*) keyname
  READ(fi,*) keyname, eig0
  CLOSE(fi)
  WRITE(*,*) "    Minimum energy : ", eig0
  omegamin = omegamin - eig0
  omegamax = omegamax - eig0
  !
END SUBROUTINE read_filename
!
! Read geometry from file
!
SUBROUTINE read_geometry()
  !
  USE fourier_val, ONLY : recipr, box, nsite, phase, irv, rindx, orb, &
  &                       nr, nreq, norb, nnode, knode, nk_line, kname
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
  USE fourier_val, ONLY : nk, kvec, nnode, nk_line, knode
  !
  IMPLICIT NONE
  !
  INTEGER :: inode, ik
  REAL(8) :: xx
  !
  nk = nk_line * (nnode - 1) + 1
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
END SUBROUTINE set_kpoints
!
! Read Correlation Function
!
SUBROUTINE read_corrindx()
  !
  USE fourier_val, ONLY : file_gindx, ncor, indx, nr, rindx, orb, norb
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, itmp, icor, nops, iops
  CHARACTER(100) :: ctmp
  !
  WRITE(*,*) 
  WRITE(*,*) "#####  Read Correlation Index File  #####" 
  WRITE(*,*) 
  !
  ALLOCATE(indx(nr,norb))
  indx(1:nr,1:norb) = 0
  !
  ! Read index for the One-Body Correlation
  !
  OPEN(fi, file = TRIM(file_gindx))
  WRITE(*,*) "  Read from ", TRIM(file_gindx)
  READ(fi,*) ctmp
  READ(fi,*) ctmp, ncor
  ncor = ncor - 1
  IF(ncor /= nr*norb) STOP "Number of correlation and NR*Norb is different."
  READ(fi,*) ctmp
  READ(fi,*) ctmp
  READ(fi,*) ctmp
  WRITE(*,*) "    Number of Correlation Function : ", ncor
  !
  ! Search mapping index for up-up and down-down correlation
  !
  READ(fi,*) nops
  DO iops = 1, nops
     READ(fi,*) itmp
  END DO
  !
  DO icor = 1, ncor
     READ(fi,*) nops
     DO iops = 1, nops
        READ(fi,*) itmp
     END DO
     indx(rindx(itmp + 1), orb(itmp + 1)) = icor 
  END DO
  !
  WRITE(*,*) "    Number of Index : ",     COUNT(indx(1:nr, 1:norb) /= 0)
  !
  CLOSE(fi)
  !
END SUBROUTINE read_corrindx
!
! Read Correlation Function
!
SUBROUTINE read_corrfile()
  !
  USE fourier_val, ONLY : filehead, ncor, indx, cor, norb, nr, nomega
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, icor, iorb, ir, iomega
  COMPLEX(8),ALLOCATABLE :: cor0(:,:)
  REAL(8) :: dtmp(4)
  !
  ALLOCATE(cor(nr,norb,nomega))
  ALLOCATE(cor0(nomega,ncor))
  cor(1:nr,1:norb,1:nomega) = CMPLX(0d0, 0d0, KIND(1d0))
  !
  OPEN(fi, file = TRIM(filehead) // "_DynamicalGreen.dat")
  !
  DO icor = 1, ncor
     DO iomega = 1, nomega
        READ(fi,*) dtmp(1:4)
        cor0(iomega,icor) = CMPLX(dtmp(3), dtmp(4), KIND(1d0))
     END DO
  END DO
  !
  CLOSE(fi)
  !
  ! Map it into Up-Up(1) and Down-Down(2) Correlation
  !
  DO iorb = 1, norb
     DO ir = 1, nr
        cor(ir, iorb, 1:nomega) = cor0(1:nomega, indx(ir, iorb))
     END DO
  END DO
  !
  DEALLOCATE(cor0, indx)
  !
END SUBROUTINE read_corrfile
!
! Fourier transformation
!
SUBROUTINE fourier_cor()
  !
  USE fourier_val, ONLY : cor, cor_k, kvec, nk, nr, nreq, norb, irv, phase, nomega
  IMPLICIT NONE
  !
  INTEGER :: ik, ir, ireq
  REAL(8) :: tpi = 2.0 * ACOS(-1d0), theta
  COMPLEX(8),ALLOCATABLE :: fmat(:,:)
  !
  ALLOCATE(fmat(nk,nr), cor_k(nk,norb,nomega))
  !
  ! Matirx for Fourier trans. exp(-i k R)
  !
  DO ik = 1, nk
     DO ir = 1, nr
        fmat(ik,ir) = CMPLX(0d0, 0d0, KIND(1d0))
        DO ireq = 1, nreq(ir)
           theta = - tpi * DOT_PRODUCT(kvec(1:3,ik), DBLE(irv(1:3,ireq,ir))) &
           &     + tpi * phase(ireq,ir)
           fmat(ik,ir) = fmat(ik,ir) + CMPLX(COS(theta), SIN(theta), KIND(1d0))
        END DO
        fmat(ik,ir) = fmat(ik,ir) / DBLE(nreq(ir))
     END DO ! ir = 1, nr
  END DO ! ik = 1, nk
  !
  CALL zgemm('N', 'N', nk, norb*nomega, nr, CMPLX(1d0, 0d0, KIND(1d0)), fmat, nk, &
  &          cor, nr, CMPLX(0d0,0d0,KIND(1d0)), cor_k, nk)
  !
  DEALLOCATE(fmat, cor)
  !
END SUBROUTINE fourier_cor
!
! Output Fourier component of Correlation function
!
SUBROUTINE output_cor()
  !
  USE fourier_val, ONLY : cor_k, nk, nnode, knode, nk_line, kname, norb, &
  &                       recipr, filehead, nomega, omegamin, omegamax
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, ik, inode, ikk, iomega
  REAL(8) :: dk(3), dk_cart(3), xk(nk), &
  &          xk_label(nnode), klength, omega
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
  WRITE(*,*) "  Correlation in k-space : ", TRIM(filehead) // "_dyn.dat"
  !
  OPEN(fo, file = TRIM(filehead) // "_dyn.dat")
  !
  DO ik = 1, ikk
     DO iomega = 1, nomega
        omega = (omegamax - omegamin) * DBLE(iomega - 1) / DBLE(nomega) + omegamin
        WRITE(fo,'(1000e15.5)') xk(ik), omega, cor_k(ik, 1:norb, iomega)
     END DO
     WRITE(fo,*)
  END DO
  !
  CLOSE(fo)
  !
  OPEN(fo, file = "kpath.gp")
  !
  WRITE(fo,'(a)',advance="no") "set xtics ("
  DO inode = 1, nnode - 1
     WRITE(fo,'(a,a,a,f10.5,a)',advance="no") "'", TRIM(kname(inode)), "'  ", xk_label(inode), ", "
  END DO
  WRITE(fo,'(a,a,a,f10.5,a)') "'", TRIM(kname(nnode)), "' ", xk_label(nnode), ")"
  WRITE(fo,'(a)') "set ylabel 'Energy from E_0'"
  WRITE(fo,'(a)') "set zlabel 'Spectrum' rotate"
  WRITE(fo,'(a)') "set ticslevel 0"
  WRITE(fo,'(a)') "set xzeroaxis"
  WRITE(fo,'(a)') "set grid xtics lt 1 lc 0"
  !
  CLOSE(fo)
  !
  DEALLOCATE(cor_k)
  !
END SUBROUTINE output_cor
!
END MODULE fourier_routine
!
! Main routine
!
PROGRAM dynamicalr2k
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
END PROGRAM dynamicalr2k
