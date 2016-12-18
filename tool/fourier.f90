MODULE fourier_val
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
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
  & file_geom, & ! Filename for geometry
  & file_one, &  ! Filename for One-body Correlation
  & file_two     ! Filename for Two-body Correlation
  !
  INTEGER,ALLOCATABLE,SAVE :: &
  & indx(:,:,:) ! (3,nsite*nsite,8) Mapping index for each Correlation function
  & equivk(:)   ! (nktot)
  !
  REAL(8),ALLOCATABLE,SAVE :: &
  & kvec(:,:), & ! (2,nk) k-vector in the !st BZ
  & kvec_tot(:,:), & ! (2,nktot)
  & site(:,:)    ! (2,nsite) Site geometry in the fractional coordinate
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & cor(:,:,:), & ! (nsite,nsite,6) Correlation function in real space (See below)
  & cor_k(:,:)    ! (nk,6) Correlation function in the k-space (See below)
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
!
!
SUBROUTINE read_filename()
  !
  USE fourier_val, ONLY : file_one, file_two, nsite
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, itmp
  CHARACTER(256) :: ctmp, modpara, calcmod
  !
  CALL getarg(1, namelist)
  !
  ! Read from NameList file
  !
  OPEN(fi,file = TRIM(namelist))
  !
  DO
     READ(fi,*,EOF=10) keyname, ctmp
     CALL key2lower(keyname)
     !
     IF(TRIM(ADJUSTL(keyname)) == "onebodyg") THEN
        file_oneindx = TRIM(ADJUSTL(ctmp))
     ELSE IF(TRIM(ADJUSTL(keyname)) == "twobodyg") THEN
        file_twoindx = TRIM(ADJUSTL(ctmp))
     ELSE IF(TRIM(ADJUSTL(keyname)) == "modpara") THEN
        modpara = TRIM(ADJUSTL(ctmp))
     ELSE IF(TRIM(ADJUSTL(keyname)) == "calcmod") THEN
        calcmod = TRIM(ADJUSTL(ctmp))
     END IF
  END DO
  !
10 WRITE(*,*) "  Read from ", TRIM(namelist)
  CLOSE(FI)
  !
  WRITE(*,*) "    OneBodyG file : ", TRIM(ADJUSTL(file_oneindx))
  WRITE(*,*) "    TwoBodyG file : ", TRIM(ADJUSTL(file_twoindx))
  WRITE(*,*) "    ModPara file : ", TRIM(ADJUSTL(modpara)) 
  WRITE(*,*) "    CalcMod file : ", TRIM(ADJUSTL(calcmod))
  !
  ! Read from Modpara file
  !
  OPEN(fi,file = TRIM(modpara))
  !
  READ(fi,*) ctmp
  READ(fi,*) ctmp
  READ(fi,*) ctmp
  READ(fi,*) ctmp
  READ(fi,*) ctmp
  READ(fi,*) keyname, filehead
  READ(fi,*) ctmp
  READ(fi,*) ctmp
  DO
     READ(fi,*,EOF=20) keyname, itmp
     CALL key2lower(keyname)
     !
     IF(TRIM(ADJUSTL(keyname)) == "nsite") THEN
        nsite = itmp
     ELSE IF(TRIM(ADJUSTL(keyname)) == "numave") THEN
        numave = itmp
     ELSE IF(TRIM(ADJUSTL(keyname)) == "lanczos_max") THEN
        lanczos_max = itmp
     ELSE IF(TRIM(ADJUSTL(keyname)) == "expecinterval") THEN
        interval = itmp
     ELSE IF(TRIM(ADJUSTL(keyname)) == "exct") THEN
        nwfc = itmp
     END IF
  END DO
  !
20 WRITE(*,*) "  Read from ", TRIM(modpara)
  WRITE(*,*) "    FileHead : ", TRIM(ADJUSTL(filehead))
  WRITE(*,*) "    Number of site : ", nsite
  WRITE(*,*) "    Number of run : ", numave
  WRITE(*,*) "    Maximum iteration : ", lanczos_max
  WRITE(*,*) "    Expectation interval : ", interval
  CLOSE(FI)
  !
  ! Read from CalcMod file
  !
  OPEN(fi,file = TRIM(calcmod))
  !
  DO
     READ(fi,*,EOF=30) keyname, itmp
     CALL key2lower(keyname)
     !
     IF(TRIM(ADJUSTL(keyname)) == "calctype") THEN
        calctype = itmp
     END IF
  END DO
  !
30 WRITE(*,*) "  Read from ", TRIM(calcmod)
  WRITE(*,*) "    CalcType : ", calctype
  CLOSE(FI)
  !
  ! Spefify the tail of file name
  !
  IF(calctype == 0) THEN ! Lanczos
     !
     nwfc = 1
     ALLOCATE(filetail(nwfc))
     filetail = ".dat"
     !
  ELSE IF (calctype == 1) THEN ! TPQ
     !
     nwfc = numave * (1 + (lanczos_max - 1) / interval)
     ALLOCATE(filetail(nwfc))
     !
     nwfc = 0
     DO irun = 0, numave - 1
        DO istep = 0, lanczos_max - 1
           IF(istep % interval == 0) THEN
              nwfc = nwfc + 1
              WRITE(cirun,*) irun
              WRITE(cistep,*) istep
              WRITE(filetail(nwfc,'(a,a,a,a,a)') "_set", TRIM(ADJUSTL(cirun)), "step", TRIM(ADJUSTL(cistep)), ".dat"
           END IF
        END DO
     END DO
     !
  ELSE IF (calctype == 2) THEN ! FullDiag
     !
     OPEN(fi, file = "output/CHECK_Memory.dat")
     READ(fi, '("  MAX DIMENSION idim_max=1", i)') nwfc
     CLOSE(fi)
     ALLOCATE(filehead(nwfc))
     !
     DO iwfc = 1, nwfc
        WRITE(cirun,*) iwfc
        WRITE(filetail(iwfc),'(a,a,a)') "_eigen", TRIM(ADJUSTL(cirun)), ".dat"
     END DO
     !     
  ELSE ! LOBCG
     !
     ALLOCATE(filehead(nwfc))
     DO iwfc = 1, nwfc
        WRITE(cirun,*) iwfc
        WRITE(filetail(iwfc),'(a,a,a)') "_eigen", TRIM(ADJUSTL(cirun)), ".dat"
     END DO
  END IF
  !
  WRITE(*,*) "  Number of states : ", nwfc
  !
END SUBROUTINE read_filename
!
! Read geometry from file
!
SUBROUTINE read_geometry()
  !
  USE fourier_val, ONLY : file_geom, file_one, file_two, nsample, &
  &                       direct, recipr, avec, bvec, nsite, site, nk
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, isite
  REAL(8) :: det
  NAMELIST /fourier/ nsite, direct, avec, file_one, file_two
  !
  CALL getarg(1, file_geom)
  !
  OPEN(fi, file = TRIM(file_geom))
  !  
  READ(*,control,err=10)
  !  
  WRITE(*,*) "  Number of samples : ", nsample
  ALLOCATE(site(2,nsite))
  !
  WRITE(*,*) "  Direct LATTICE VECTOR :"
  WRITE(*,'(2f15.10)') direct(1:2, 1:2)
  !
  WRITE(*,*) "  Supercell Index :"
  WRITE(*,'(2i8)') avec(1:2, 1:2)
  !
  DO isite = 1, nsite
     READ(fi,*) site(1:2,isite)
  END DO
  !
  WRITE(*,*) "  One Body Correlation File: ", TRIM(file_one)
  WRITE(*,*) "  Two Body Correlation File: ", TRIM(file_two)
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
  WRITE(*,*) "  Number of k point : ", nk
  WRITE(*,*) "  Reciplocal superlattice vector :"
  WRITE(*,'(2i8)') bvec(1:2, 1:2)
  !
  ! Compute Reciprocal Lattice Vector
  !
  det = recipr(1,1) * recipr(2,2) - recipr(1,2) * recipr(2,1)
  recipr(1:2,1) = (/  direct(2,2), - direct(1,2)/) / det
  recipr(1:2,2) = (/- direct(2,1),   direct(1,1)/) / det
  WRITE(*,*) "  Reciplocal lattice vector :"
  WRITE(*,'(2f15.10)') recipr(1:2, 1:2)
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
  INTEGER :: imax(2), imin(2), i1, i2, ii, edge(2,4), nk0, ik, jk
  INTEGER,ALLOCATABLE :: ikvec(:,:), ikvec_tot(:,:)
  !
  ii = 0
  DO i2 = -1, 1, 2
     DO i1 = -1, 1, 2
        ii = ii + 1
        edge(1:2,ii) = MATMUL((/i1, i2/) * avec(1:2,1:2))
     END DO
  END DO
  !
  DO idim = 1, 2
     imin(idim) = MINVAL(edge(idim,1:4))
     imax(idim) = MAXVAL(edge(idim,1:4))
  END DO
  !
  nktot = PRODUCT(imax(1:2) - imin(1:2) + 1)
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
        IF(ALL(0 <= ikvec_tot(1:2,nktot)) .AND. ALL(ikvec_tot(1:2,nktot) < nk)) THEN
           nk0 = nk0 + 1
           ikvec(1:2,nk0) = ikvec_tot(1:2,nktot)
           kvec( 1:2,nk0) = kvec_tot( 1:2,nktot)
        END IF
        !
     END DO
  END DO
  WRITE(*,*) "  Number of k : ", nk0
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
  USE fourier_val, ONLY : file_oneindx, file_twoindx, ncor1, ncor2, ncor, indx
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, itmp(8)
  CHARACTER(100) :: ctmp
  !
  ALLOCATE(indx(3,nsite*nsite,8))
  !
  ! Read index for the One-Body Correlation
  !
  OPEN(fi, file = TRIM(file_oneindx))
  READ(fi,*) ctmp
  READ(fi,*) ctmp, ncor1
  READ(fi,*) ctmp
  READ(fi,*) ctmp
  READ(fi,*) ctmp
  WRITE(*,*) "  Number of Correlation Function (One body) : ", ncor1
  !
  ncor(1:2) = 0
  DO icor = 1, ncor1
     READ(fi,*) itmp(1), itmp(2), itmp(3), itmp(4)
     IF(itmp(2) == 0 .AND. itmp(4) == 0) THEN
        ncor(1) = ncor(1) + 1
        indx(1:3,ncor(1),1) = (/itmp(1), itmp(3), icor/)
     ELSE IF (itmp(2) == 1 .AND. itmp(4) == 1) THEN
        ncor(2) = ncor(2) + 1
        indx(1:3,ncor(2),2) = (/itmp(1), itmp(3), icor/)
     END IF
  END DO
  !
  CLOSE(fi)
  !
  ! Read index for the One-Body Correlation
  !
  OPEN(fi, file = TRIM(file_oneindx))
  READ(fi,*) ctmp
  READ(fi,*) ctmp, ncor2
  READ(fi,*) ctmp
  READ(fi,*) ctmp
  READ(fi,*) ctmp
  WRITE(*,*) "  Number of Correlation Function (Two body) : ", ncor2
  !
  ncor(3:8) = 0
  DO icor = 1, ncor2
     !
     READ(fi,*) itmp(1:8)
     IF(itmp(1) /= itmp(3) .OR. itmp(5) /= itmp(7)) CYCLE
     !
     IF(itmp(2) == 0 .AND. itmp(4) == 0) THEN
        IF(itmp(6) == 0 .AND. itmp(8) == 0) THEN
           ncor(3) = ncor(3) + 1
           indx(1:3,ncor(3),3) = (/isite, jsite, icor/)
        ELSE IF(itmp(6) == 1 .AND. itmp(8) == 1) THEN
           ncor(4) = ncor(4) + 1
           indx(1:3,ncor(4),4) = (/isite, jsite, icor/)
        END IF
     ELSE IF(itmp(2) == 1 .AND. itmp(4) == 1) THEN
        IF(itmp(6) == 0 .AND. itmp(8) == 0) THEN
           ncor(5) = ncor(5) + 1
           indx(1:3,ncor(5),5) = (/isite, jsite, icor/)
        ELSE IF(itmp(6) == 1 .AND. itmp(8) == 1) THEN
           ncor(6) = ncor(6) + 1
           indx(1:3,ncor(6),6) = (/isite, jsite, icor/)
        END IF
     ELSE IF((itmp(2) == 0 .AND. itmp(4) == 1) .AND. (itmp(6) == 1 .AND. itmp(4) == 0)) THEN
        ncor(7) = ncor(7) + 1
        indx(1:3,ncor(7),7) = (/isite, jsite, icor/)
     ELSE IF((itmp(2) == 1 .AND. itmp(4) == 0) .AND. (itmp(6) == 0 .AND. itmp(4) == 1)) THEN
        ncor(8) = ncor(8) + 1
        indx(1:3,ncor(8),8) = (/isite, jsite, icor/)
     END IF
  END DO
  !
  CLOSE(fi)
  !
END SUBROUTINE read_corrindx
!
! Read Correlation Function
!
SUBROUTINE read_corrfile()
  !
  USE fourier_val, ONLY : file_one, file_two, ncor1, ncor2, ncor, indx, cor, nsite
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, icor, itmp(8), ii
  COMPLEX(8),ALLOCATABLE :: cor0(:)
  CHARACTER(256) :: filename, set, step
  !
  ALLOCATE(cor(nsite,nsite,6,nwfc), cor0(MAX(ncor1,ncor2))
  cor(1:nsite,1:nsite,1:6,1:nwfc) = CMPLX(0d0, 0d0)
  !
  ! One Body Correlation Function
  !
  DO iwfc = 1, nwfc
     !
     filename = TRIM(filehead) // "_cisajs_" // TRIM(filetail(iwfc))
     OPEN(fi, file = TRIM(filename))
     !
     DO icor = 1, ncor1
        READ(fi,*) itmp(1:4), cor1(icor)
     END DO
     CLOSE(fi)
     !
     ! Map it into Up-Up(1) and Down-Down(2) Correlation
     !
     DO ii = 1, 2
        DO icor = 1, ncor(ii)
           cor(indx(1,icor,ii), indx(2,icor,ii), ii, iwfc) = cor1(indx(3,icor,ii))
        END DO
     END DO
     !
     ! Two Body Correlation function
     !
     filename = TRIM(filehead) // "_cisajscktalt_" // TRIM(filetail(iwfc))
     OPEN(fi, file = TRIM(filename))
     !
     DO icor = 1, ncor2
        READ(fi,*) itmp(1:8), cor2(icor)
     END DO
     CLOSE(fi)
     !
     ! Map it into Density-Density(3), Sz-Sz(4), S+S-(5), S-S+(6) Correlation
     !
     DO ii = 3, 6, 3
        DO icor = 1, ncor(3)
           cor(indx(1,icor,ii), indx(2,icor,ii), 3, iwfc) &
           &  = cor( indx(1,icor,ii), indx(2,icor,ii), 3, iwfc) &
           &  + cor2(indx(3,icor,ii))
           cor(indx(1,icor,ii), indx(2,icor,ii), 4, iwfc) &
           &  = cor( indx(1,icor,ii), indx(2,icor,ii), 4, iwfc) &
           &  + cor2(indx(3,icor,ii))
        END DO
     END DO
     DO ii = 4, 5
        DO icor = 1, ncor(3)
           cor(indx(1,icor,ii), indx(2,icor,ii), 3, iwfc) &
           &  = cor( indx(1,icor,ii), indx(2,icor,ii), 3) &
           &  + cor2(indx(3,icor,ii))
           cor(indx(1,icor,ii), indx(2,icor,ii), 4, iwfc) &
           &  = cor( indx(1,icor,ii), indx(2,icor,ii), 4) &
           &  - cor2(indx(3,icor,ii))
        END DO
     END DO
     cor(1:nsite,1:nsite, 4, iwfc) = cor(1:nsite,1:nsite, 4, iwfc) * 0.25d0
     DO ii = 7, 8
        DO icor = 1, ncor(ii)
           cor(indx(1,icor,ii), indx(2,icor,ii), ii-2, iwfc) = cor2(indx(3,icor,ii))
        END DO
     END DO
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
  USE fourier_val, ONLY : nsite, cor, cor_k, site
  IMPLICIT NONE
  !
  INTEGER :: isite, jsite
  REAL(8) :: tpi = 2.0 * ACOS(-1d0), theta
  COMPLEX(8),ALLOCATABLE :: fmat(:,:)
  !
  ALLOCATE(fmat(nk,nsite,nsite), cor_k(nk,6,nwfc))
  !
  DO jsite = 1, nsite
     DO isite = 1, nsite
        DO ik = 1, nk
           theta = - tpi * DOT_PRODUCT(kvec(1:2,ik), (site(1:2,isite) - site(1:2,jsite)))
           fmat(ik,isite) = CMPLX(COS(theta), SIN(theta))
        END DO
     END DO
  END DO
  !
  CALL zgemm('N', 'N', nk, 6*nwfc, nsite*nsite, CMPLX(1d0, 0d0), fmat, nk, &
  &          cor, nsite*nsite, CMPLX(0d0,0d0), cor_k, nk)
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
  USE fourier_val, ONLY : cor_k, nk, nk_tot, nk_row, kvec, kvec_tot, equiv
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, ik
  !
  ! Output Correlation function in the 1st BZ
  !
  DO iwfc = 1, nwfc
     !
     filename = TRIM(filehead) // "_corr_k_" // TRIM(filetail(iwfc))
     OPEN(fo, file = TRIM(filename))
     !
     WRITE(fo,*) "# k1[3] k2[4](Cart.) UpUp[5,6] DownDown[7,8] Density[9:10] SzSz[11,12] S+S-[13,14] S-S+[15,16]"
     !
     DO ik = 1, nk
        WRITE(fo,'(20e15.5)') MATMUL(recipr(1:2,1:2), kvec(1:2,ik)), cor_k(ik,1:6,iwfc)
     END DO
     !
     CLOSE(fo)
     !
  END DO
  !
  ! k-points in the larger area
  !
  OPEN(fo, file = "kpoints.dat")
  !
  WRITE(fo,*) nktot, nk_row
  !
  DO ik = 1, nk_tot
     WRITE(fo,'(2e15.5,i7)') MATMUL(recipr(1:2,1:2), kvec_tot(1:2,ik)), equiv(ik)
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
PROGRAM fourier
  !
  USE fourier_routine, ONLY : read_filename, read_geometry, set_kpints, &
  &                           read_corrindx, read_corrfile, fourier_cor, output_cor
  IMPLICIT NONE
  !
  CALL read_filename()
  CALL read_geometry()
  CALL set_kpints()
  CALL read_corrindx()
  CALL read_corrfile()
  CALL fourier_cor()
  CALL output_cor()
  !
END PROGRAM fourier
