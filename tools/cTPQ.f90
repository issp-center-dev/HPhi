program cTPQ
  !
  implicit none
  !
  integer :: iave, iLan, jLan, fi = 10, fo = 20, nAve, nLan, nSite
  real(8) :: Large, coeff, beta, bb, bHb, bHHb, vecN, vecH, vecH2
  real(8),allocatable :: norm_m(:,:), ene_m(:,:), ene2_m(:,:), beta_m(:,:), &
  &                      norm_c(:,:), ene_c(:,:), ene2_c(:,:)
  logical,allocatable :: lscale(:,:)
  character(256) :: ctemp, command
  character(256),allocatable :: ss_rand(:), norm_rand(:)
  !
  if(iargc() /= 1) then
     write(*,*)
     write(*,*) "Usage :"
     write(*,*) "$ cTPQ.x {Modpara file}"
     write(*,*) 
     stop
  end if
  !
  call getarg(1, ctemp)
  !
  write(command, *) "grep -i largevalue ", trim(ctemp), "| awk '{print $2}' > __HPhi_temp__"
  write(*,*) trim(command)
  call system(trim(command))
  !
  write(command, *) "grep -i lanczos_max ", trim(ctemp), "| awk '{print $2}' >> __HPhi_temp__"
  write(*,*) trim(command)
  call system(trim(command))
  !
  write(command, *) "grep -i nsite ", trim(ctemp), "| awk '{print $2}' >> __HPhi_temp__"
  write(*,*) trim(command)
  call system(trim(command))
  !
  call system("find ./ -name 'SS_rand*.dat' | wc -l >> __HPhi_temp__")
  call system("find ./ -name 'SS_rand*.dat' | sort >> __HPhi_temp__")
  call system("find ./ -name 'Norm_rand*.dat' | sort >> __HPhi_temp__")
  !
  open(fi, file = "__HPhi_temp__")
  !
  read(fi,*) large
  read(fi,*) nLan
  nLan = nLan - 1
  read(fi,*) nSite
  read(fi,*) nAve
  !
  write(*,*) "  Lanczos_max : ", nLan + 1
  write(*,*) "        Nsite : ", nSite
  write(*,*) "   LargeValue : ", Large
  write(*,*) "       NumAve : ", nAve
  !
  allocate(norm_m(0:nLan-1, nave), ene_m(0:nLan-1, nave), ene2_m(0:nLan-1, nave), beta_m(0:nLan-1, nave), &
  &        norm_c(0:nLan-1, nave), ene_c(0:nLan-1, nave), ene2_c(0:nLan-1, nave), lscale(nLan, 0:nLan-1), &
  &        ss_rand(nAve), norm_rand(nAve)  )
  !
  do iAve = 1, nAve
     read(fi,'(a)') ss_rand(iAve)
  end do
  do iAve = 1, nAve
     read(fi,'(a)') norm_rand(iAve)
  end do
  close(fi)
  call system("rm __HPhi_temp__")
  !
  do iave = 1, nave
     !
     ! Read SS_rand
     !
     open(fi, file = trim(ss_rand(iAve)))
     write(*,*) "  Read from ", trim(ss_rand(iAve))
     read(fi,*) ctemp
     !
     do iLan = 0, nLan - 1
        read(fi,*) beta_m(iLan, iave), ene_m(iLan, iave), ene2_m(iLan, iave)
     end do
     !
     close(fi)
     !
     ! Read Norm_rand
     !
     open(fi, file = trim(norm_rand(iAve)))
     write(*,*) "  Read from ", trim(norm_rand(iAve))
     read(fi,*) ctemp
     !
     coeff = 1.0d0
     !
     do iLan = 0, nLan - 1
        read(fi,*) beta, norm_m(iLan, iave)
        norm_m(iLan, iave) = norm_m(iLan, iave)**2 * coeff
        coeff = dble(nSite * nSite) / dble(2 * (iLan+1) * (2*iLan + 1))
     end do
     !
     close(fi)
     !
  enddo
  !
  do iave = 1, nave
     !
     do iLan = 0, nLan - 1
        !
        beta = beta_m(iLan, 1)!minval(beta_m(iLan, 1:nave))
        coeff = 1.0d0
        !
        bb   = 0.0d0
        bHb  = 0.0d0
        bHHb = 0.0d0
        !
        do jLan = 1, nLan
           !
           vecN  = norm_m(nLan-jLan, iave)
           vecH  =  ene_m(nLan-jLan, iave)
           vecH2 = ene2_m(nLan-jLan, iave)
           !
           bb = beta * beta * vecN * &
           &  (bb + coeff * ((1.0d0 + beta * dble(nSite) * large / dble(2 * (nLan-jLan) + 1)) &
           &                - beta * vecH / dble(2*(nLan-jLan) + 1) &
           &                ) &
           &  )
           bHb = beta * beta * vecN * &
           &  (bHb + coeff * ((1.0d0 + beta * dble(nSite) * large / dble(2 * (nLan-jLan) + 1)) * vecH &
           &                 - beta*vecH2 / dble(2 * (nLan-jLan) + 1) &
           &                 ) &
           &  )
           bHHb = beta * beta * vecN * &
           &  (bHHb + coeff * ((1.0d0 - beta * dble(nSite) * large / dble(2 * (nLan-jLan) + 1)) * vecH2 &
           &                  + (beta * dble(nSite*nSite) * large * large / dble(2 * (nLan-jLan) + 1) &
           &                    - dble(2 * (nLan-jLan)) / beta &
           &                    ) * vecH &
           &                  ) &
           &  )
           !
           if(iave == 1) then
              if(bb > 1.0d+10) then
                 bb    =    bb * (1.0d-10)
                 bHb   =   bHb * (1.0d-10)
                 bHHb  =  bHHb * (1.0d-10)
                 coeff = coeff * (1.0d-10)
                 lscale(jLan,iLan) = .true.
              else
                 lscale(jLan,iLan) = .false.
              endif
           else
              if(lscale(jLan,iLan)) then
                 bb    =    bb * (1.0d-10)
                 bHb   =   bHb * (1.0d-10)
                 bHHb  =  bHHb * (1.0d-10)
                 coeff = coeff * (1.0d-10)
              endif
           end if
           !
        enddo
        !
        norm_c(iLan, iave) = bb
        ene_c( iLan, iave) = bHb
        ene2_c(iLan, iave) = bHHb
        !
     end do
     !
  enddo
  !
  open(fo, file = "sh.dat")
  !
  do iLan = 0, nLan - 1
     !
     beta = beta_m(iLan, 1)!minval(beta_m(iLan, 1:nave))
     !
     bb   = sum(norm_c(iLan, 1:nAve)) / dble(nAve)
     bHb  = sum( ene_c(iLan, 1:nAve)) / dble(nAve)
     bHHb = sum(ene2_c(iLan, 1:nAve)) / dble(nAve)
     !
     vecN   = sqrt(sum((norm_c(iLan, 1:nAve) -   bb)**2)) / dble(nAve)
     vecH   = sqrt(sum(( ene_c(iLan, 1:nAve) -  bHb)**2)) / dble(nAve)
     vecH2  = sqrt(sum((ene2_c(iLan, 1:nAve) - bHHb)**2)) / dble(nAve)
     !
     write(fo,'(10e25.15)') beta, &
     &  bHb / bb, abs(vecH / bb) + abs(bHb * vecN / bb**2), &
     &  bHHb / bb, abs(vecH2 / bb) + abs(bHHb * vecN / bb**2), &
     &  beta**2 * (bHHb / bb - (bHb / bb)**2), &
     &  beta**2 * ( abs(vecH2 / bb) + abs(2 * bHHb * vecH / bb**2) &
     &            + abs(-bHHb + 2*bHb**2/bb) / bb**2) * vecN
     !
  end do
  !
  close(fo)
  !
end program cTPQ
