PROGRAM cTPQ
  !
  IMPLICIT NONE
  !
  INTEGER :: iave, iLan, jLan, fi = 10, fo = 20, nAve, nLan, nSite
  REAL(8) :: Large, coeff, beta, bb, bHb, bHHb, vecN, vecH, vecH2
  REAL(8),ALLOCATABLE :: norm_m(:,:), ene_m(:,:), ene2_m(:,:), beta_m(:,:), &
  &                      norm_c(:,:), ene_c(:,:), ene2_c(:,:)
  LOGICAL,ALLOCATABLE :: lscale(:,:)
  CHARACTER(256) :: ctemp, command
  CHARACTER(256),ALLOCATABLE :: ss_rand(:), norm_rand(:)
  !
  IF(iargc() /= 1) THEN
     WRITE(*,*)
     WRITE(*,*) "Usage :"
     WRITE(*,*) "$ cTPQ.x {Modpara file}"
     WRITE(*,*) 
     STOP
  END IF
  !
  CALL getarg(1, ctemp)
  !
  WRITE(command, *) "grep -i largevalue ", TRIM(ctemp), "| awk '{print $2}' > __HPhi_temp__"
  WRITE(*,*) TRIM(command)
  CALL system(TRIM(command))
  !
  WRITE(command, *) "grep -i lanczos_max ", TRIM(ctemp), "| awk '{print $2}' >> __HPhi_temp__"
  WRITE(*,*) TRIM(command)
  CALL system(TRIM(command))
  !
  WRITE(command, *) "grep -i nsite ", TRIM(ctemp), "| awk '{print $2}' >> __HPhi_temp__"
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
  ALLOCATE(norm_m(0:nLan-1, nave), ene_m(0:nLan-1, nave), ene2_m(0:nLan-1, nave), beta_m(0:nLan-1, nave), &
  &        norm_c(0:nLan-1, nave), ene_c(0:nLan-1, nave), ene2_c(0:nLan-1, nave), lscale(nLan, 0:nLan-1), &
  &        ss_rand(nAve), norm_rand(nAve)  )
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
  DO iave = 1, nave
     !
     DO iLan = 0, nLan - 1
        !
        beta = beta_m(iLan, 1)!minval(beta_m(iLan, 1:nave))
        coeff = 1.0d0
        !
        bb   = 0.0d0
        bHb  = 0.0d0
        bHHb = 0.0d0
        !
        DO jLan = 1, nLan
           !
           vecN  = norm_m(nLan-jLan, iave)
           vecH  =  ene_m(nLan-jLan, iave)
           vecH2 = ene2_m(nLan-jLan, iave)
           !
           bb = beta * beta * vecN * &
           &  (bb + coeff * ((1.0d0 + beta * DBLE(nSite) * large / DBLE(2 * (nLan-jLan) + 1)) &
           &                - beta * vecH / DBLE(2*(nLan-jLan) + 1) &
           &                ) &
           &  )
           bHb = beta * beta * vecN * &
           &  (bHb + coeff * ((1.0d0 + beta * DBLE(nSite) * large / DBLE(2 * (nLan-jLan) + 1)) * vecH &
           &                 - beta*vecH2 / DBLE(2 * (nLan-jLan) + 1) &
           &                 ) &
           &  )
           bHHb = beta * beta * vecN * &
           &  (bHHb + coeff * ((1.0d0 - beta * DBLE(nSite) * large / DBLE(2 * (nLan-jLan) + 1)) * vecH2 &
           &                  + (beta * DBLE(nSite*nSite) * large * large / DBLE(2 * (nLan-jLan) + 1) &
           &                    - DBLE(2 * (nLan-jLan)) / beta &
           &                    ) * vecH &
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
        ENDdo
        !
        norm_c(iLan, iave) = bb
        ene_c( iLan, iave) = bHb
        ene2_c(iLan, iave) = bHHb
        !
     END DO
     !
  END DO
  !
  OPEN(fo, file = "sh.dat")
  !
  DO iLan = 0, nLan - 1
     !
     beta = beta_m(iLan, 1)!minval(beta_m(iLan, 1:nave))
     !
     bb   = SUM(norm_c(iLan, 1:nAve)) / DBLE(nAve)
     bHb  = SUM( ene_c(iLan, 1:nAve)) / DBLE(nAve)
     bHHb = SUM(ene2_c(iLan, 1:nAve)) / DBLE(nAve)
     !
     vecN   = SQRT(SUM((norm_c(iLan, 1:nAve) -   bb)**2)) / DBLE(nAve)
     vecH   = SQRT(SUM(( ene_c(iLan, 1:nAve) -  bHb)**2)) / DBLE(nAve)
     vecH2  = SQRT(SUM((ene2_c(iLan, 1:nAve) - bHHb)**2)) / DBLE(nAve)
     !
     WRITE(fo,'(10e25.15)') beta, &
     &  bHb / bb, ABS(vecH / bb) + ABS(bHb * vecN / bb**2), &
     &  bHHb / bb, ABS(vecH2 / bb) + ABS(bHHb * vecN / bb**2), &
     &  beta**2 * (bHHb / bb - (bHb / bb)**2), &
     &  beta**2 * ( ABS(vecH2 / bb) + ABS(2 * bHHb * vecH / bb**2) &
     &            + ABS(-bHHb + 2*bHb**2/bb) / bb**2) * vecN
     !
  END DO
  !
  CLOSE(fo)
  !
END PROGRAM cTPQ
