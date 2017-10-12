MODULE cTPQ_mod
  !
  IMPLICIT NONE
  !
CONTAINS
  !
SUBROUTINE read_modpara(nLan, nSite, nAve, Large, ss_rand, norm_rand)
  !
#if defined(FUJITSU)
  USE service_routines, ONLY : IARGC
#endif
  IMPLICIT NONE
  !
  INTEGER,INTENT(OUT) :: nLan, nSite, nAve
  REAL(8),INTENT(OUT) :: Large
  CHARACTER(256),ALLOCATABLE,INTENT(OUT) :: ss_rand(:), norm_rand(:)
  !
  INTEGER :: fi = 10, iAve
  CHARACTER(256) :: modpara, command
#if defined(SR)
  INTEGER,INTRINSIC :: IARGC
#endif
  !
  WRITE(*,*)
  WRITE(*,*) "#################################################"
  WRITE(*,*) "  CAUTION ! This program is very EXPERIMENTAL !  "
  WRITE(*,*) "        Please use it at your own risk !         "
  WRITE(*,*) "#################################################"
  WRITE(*,*)
  !
  IF(IARGC() /= 1) THEN
     WRITE(*,*)
     WRITE(*,*) "Usage :"
     WRITE(*,*) "$ cTPQ {Modpara file}"
     WRITE(*,*) 
     STOP
  END IF
  !
  CALL getarg(1, modpara)
  !
  WRITE(command, *) "grep -i largevalue ", TRIM(modpara), "| awk '{print $2}' > __HPhi_temp__"
  WRITE(*,*) TRIM(command)
  CALL system(TRIM(command))
  !
  WRITE(command, *) "grep -i lanczos_max ", TRIM(modpara), "| awk '{print $2}' >> __HPhi_temp__"
  WRITE(*,*) TRIM(command)
  CALL system(TRIM(command))
  !
  WRITE(command, *) "grep -i nsite ", TRIM(modpara), "| awk '{print $2}' >> __HPhi_temp__"
  WRITE(*,*) TRIM(command)
  CALL system(TRIM(command))
  !
  CALL system("find ./ -name 'SS_rand*.dat' | wc -l >> __HPhi_temp__")
  CALL system("find ./ -name 'SS_rand*.dat' | sort >> __HPhi_temp__")
  CALL system("find ./ -name 'Norm_rand*.dat' | sort >> __HPhi_temp__")
  !
  OPEN(fi, file = "__HPhi_temp__")
  !
  READ(fi,*) large
  READ(fi,*) nLan
  nLan = nLan - 1
  READ(fi,*) nSite
  READ(fi,*) nAve
  !
  WRITE(*,*) "  Lanczos_max : ", nLan + 1
  WRITE(*,*) "        Nsite : ", nSite
  WRITE(*,*) "   LargeValue : ", Large
  WRITE(*,*) "       NumAve : ", nAve
  !
  ALLOCATE(ss_rand(nAve), norm_rand(nAve))
  !
  DO iAve = 1, nAve
     READ(fi,'(a)') ss_rand(iAve)
  END DO
  DO iAve = 1, nAve
     READ(fi,'(a)') norm_rand(iAve)
  END DO
  CLOSE(fi)
  CALL system("rm __HPhi_temp__")
  !
END SUBROUTINE read_modpara
!
SUBROUTINE read_ssrand(nAve,nLan,nSite,ss_rand,norm_rand,norm_m,ene_m,ene2_m,beta_m)
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: nAve, nLan, nSite
  CHARACTER(256),INTENT(IN) :: ss_rand(nAve), norm_rand(nAve)
  REAL(8),INTENT(OUT) :: norm_m(0:nLan-1, nave), ene_m(0:nLan-1, nave), &
  &                      ene2_m(0:nLan-1, nave), beta_m(0:nLan-1, nave)
  !
  INTEGER :: fi = 10, iave, iLan
  REAL(8) :: coeff, beta
  CHARACTER(256) :: ctemp
  !
  DO iave = 1, nave
     !
     ! Read SS_rand
     !
     OPEN(fi, file = TRIM(ss_rand(iAve)))
     WRITE(*,*) "  Read from ", TRIM(ss_rand(iAve))
     READ(fi,*) ctemp
     !
     DO iLan = 0, nLan - 1
        READ(fi,*) beta_m(iLan, iave), ene_m(iLan, iave), ene2_m(iLan, iave)
     END DO
     !
     CLOSE(fi)
     !
     ! Read Norm_rand
     !
     OPEN(fi, file = TRIM(norm_rand(iAve)))
     WRITE(*,*) "  Read from ", TRIM(norm_rand(iAve))
     READ(fi,*) ctemp
     !
     coeff = 1.0d0
     !
     DO iLan = 0, nLan - 1
        READ(fi,*) beta, norm_m(iLan, iave)
        norm_m(iLan, iave) = norm_m(iLan, iave)**2 * coeff
        coeff = DBLE(nSite * nSite) / DBLE(2 * (iLan+1) * (2*iLan + 1))
     END DO
     !
     CLOSE(fi)
     !
  ENDdo
  !
END SUBROUTINE read_ssrand
!
!
!
SUBROUTINE compute_cTPQ(nLan,nAve,nSite,norm_m,ene_m,ene2_m,beta_m,Large,norm_c,ene_c,ene2_c)
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: nLan, nAve, nSite
  REAL(8),INTENT(IN) :: norm_m(0:nLan-1, nave), ene_m(0:nLan-1, nave), &
  &                     ene2_m(0:nLan-1, nave), beta_m(0:nLan-1, nave), Large
  REAL(8),INTENT(OUT) :: norm_c(0:nLan-1, nave), ene_c(0:nLan-1, nave), ene2_c(0:nLan-1, nave)
  !
  INTEGER :: iave, iLan, jLan
  REAL(8) :: beta, coeff, bb, bHb, bHHb, norm, E, E2
  LOGICAL :: lscale(0:nLan-1, 0:nLan-1)
  !
  DO iave = 1, nave
     !
     DO iLan = 0, nLan - 1
        !
        beta = beta_m(iLan, 1)
        coeff = 1.0d0
        !
        bb   = 0.0d0
        bHb  = 0.0d0
        bHHb = 0.0d0
        !
        DO jLan = nLan - 1, 0, -1
           !
           norm  = norm_m(jLan, iave)
           E  =  ene_m(jLan, iave)
           E2 = ene2_m(jLan, iave)
           !
           bb = beta * beta * norm * &
           &  (bb + coeff * (1.0d0 + beta * (DBLE(nSite) * large- E) / DBLE(2*jLan + 1)))
           bHb = beta * beta * norm * &
           &  (bHb + coeff * (E + beta * (DBLE(nSite) * large * E - E2) / DBLE(2*jLan + 1)))
           bHHb = beta * beta * norm * &
           &  (bHHb + coeff * ( E2 - DBLE(2*jLan) / beta * E &
           &                  - beta * DBLE(nSite) * large * (E2 - nsite*Large*E) / DBLE(2*jLan + 1) &
           &                  ) &
           &  )
           !
           IF(iave == 1) THEN
              IF(bb > 1.0d+10) THEN
                 bb    =    bb * (1.0d-10)
                 bHb   =   bHb * (1.0d-10)
                 bHHb  =  bHHb * (1.0d-10)
                 coeff = coeff * (1.0d-10)
                 lscale(jLan,iLan) = .TRUE.
              ELSE
                 lscale(jLan,iLan) = .FALSE.
              END IF
           ELSE
              IF(lscale(jLan,iLan)) THEN
                 bb    =    bb * (1.0d-10)
                 bHb   =   bHb * (1.0d-10)
                 bHHb  =  bHHb * (1.0d-10)
                 coeff = coeff * (1.0d-10)
              END IF
           END IF
           !
        END DO ! jLan = nlan-1, 0, -1
        !
        norm_c(iLan, iave) = bb
        ene_c( iLan, iave) = bHb
        ene2_c(iLan, iave) = bHHb
        !
     END DO ! iLan = 0, nLan - 1
     !
  END DO ! iave = 1, nave
  !
END SUBROUTINE compute_cTPQ
!
END MODULE cTPQ_mod
!
PROGRAM cTPQ
  !
  USE cTPQ_mod
  IMPLICIT NONE
  !
  INTEGER :: iLan, fo = 20, nAve, nLan, nSite
  REAL(8) :: Large, beta, bb, bHb, bHHb, dbb, dbHb, dbHHb
  REAL(8),ALLOCATABLE :: norm_m(:,:), ene_m(:,:), ene2_m(:,:), beta_m(:,:), &
  &                      norm_c(:,:), ene_c(:,:), ene2_c(:,:)
  CHARACTER(256),ALLOCATABLE :: ss_rand(:), norm_rand(:)
  !
  CALL read_modpara(nLan, nSite, nAve, Large, ss_rand, norm_rand)
  !
  ALLOCATE(norm_m(0:nLan-1, nave), ene_m(0:nLan-1, nave), ene2_m(0:nLan-1, nave), beta_m(0:nLan-1, nave), &
  &        norm_c(0:nLan-1, nave), ene_c(0:nLan-1, nave), ene2_c(0:nLan-1, nave))
  !
  CALL read_ssrand(nAve, nLan,nSite,ss_rand, norm_rand, norm_m, ene_m, ene2_m, beta_m)
  !
  CALL compute_cTPQ(nLan,nAve,nSite,norm_m,ene_m,ene2_m,beta_m,Large,norm_c,ene_c,ene2_c)
  !
  OPEN(fo, file = "sh.dat")
  !
  WRITE(fo,'(a)') "# Temperature, Energy, Energy-error, Variance, Variance-error, Specific heat, its error"
  !
  DO iLan = 0, nLan - 1
     !
     beta = beta_m(iLan, 1)
     !
     write(50,'(100e15.5)') norm_c(iLan, 1:nAve)
     bb   = SUM(norm_c(iLan, 1:nAve)) / DBLE(nAve)
     bHb  = SUM( ene_c(iLan, 1:nAve)) / DBLE(nAve)
     bHHb = SUM(ene2_c(iLan, 1:nAve)) / DBLE(nAve)
     !
     IF(nAve == 1) THEN
        dbb   = 0d0
        dbHb   = 0d0
        dbHHb  = 0d0
     ELSE
        dbb    = SQRT(SUM((norm_c(iLan, 1:nAve) -   bb)**2)) / SQRT(DBLE(nAve*(nAve-1)))
        dbHb   = SQRT(SUM(( ene_c(iLan, 1:nAve) -  bHb)**2)) / SQRT(DBLE(nAve*(nAve-1)))
        dbHHb  = SQRT(SUM((ene2_c(iLan, 1:nAve) - bHHb)**2)) / SQRT(DBLE(nAve*(nAve-1)))
     END IF
     !
     dbb   = dbb   / bb
     dbHb  = dbHb  / bHb
     dbHHb = dbHHb / bHHb
     !
     bHb  = bHb  / bb
     bHHb = bHHb / bb
     !
     WRITE(fo,'(10e25.15)') beta, &
     &  bHb,  ABS(bHb )*(ABS(dbHb )+ABS(dbb)), &
     &  bHHb, ABS(bHHb)*(ABS(dbHHb)+ABS(dbb)), &
     &  beta**2 * (bHHb - bHb**2), &
     &  beta**2 * ( ABS(bHHb) * ABS(dbHHb) + 2 * bHb**2 * ABS(dbHb ) &
     &            + ABS((bHHb) - 2 * bHb**2) * ABS(dbb) )
     !
  END DO
  !
  CLOSE(fo)
  !
  WRITE(*,*)
  WRITE(*,*) "DONE."
  WRITE(*,*)
END PROGRAM cTPQ
