program proc

  implicit none

  integer :: j,i,ell,L,N
  real(8) :: large
  real(8) :: coeff
  real(8) :: tmpb,tmpn1,tmpn2
  real(8) :: tmpnew,tmpold,beta
  real(8) :: bb,bHb,bHHb
  real(8),dimension(:),allocatable :: vecH,vecH2,vecN
  !
  ! input: fort.10
  ! use L Lanczos steps for the N site system
  ! large is the large value in TPQ
  !
  read(10,*)L,N,large
  allocate(vecH(0:L-1))
  allocate(vecH2(0:L-1))
  allocate(vecN(0:L-1))
  !
  tmpold = 1.0d0
  coeff = 1.0d0
  !
  do j = 0, L - 1
     !
     ! input: fort.11 SS_rand*.dat (please remove header starting with #) 
     !
     read(11,*)tmpb,vecH(j),vecH2(j),tmpn1,tmpn2,i
     !
     ! input: fort.12 Norm_rand*.dat (please remove header starting with #) 
     !
     read(12,*)tmpb,tmpnew,tmpn1,i
     tmpold=tmpnew*tmpnew*coeff
     vecN(j)=tmpold
     coeff=dble(N*N)/dble(2*(j+1)*(2*j+1))
     write(*,*)"#vecN",j,vecN(j)
     !
  enddo
  !
  write(*,*)"#",vecH(0),vecN(L-1)
  !
  ! loop for beta (inverse temperature)
  !
  do ell=0,999
     !
     beta=0.01d0*dexp(0.009d0*dble(ell))
     coeff=1.0d0
     bb  =0.0d0
     bHb =0.0d0
     bHHb=0.0d0
     !
     do j=1,L
        bb = beta * beta * vecN(L-j) * &
        &  (bb + coeff * ((1.0d0 + beta * dble(N) * large / dble(2 * (L-j) + 1)) &
        &                - beta * vecH(L-j) / dble(2*(L-j) + 1) &
        &                ) &
        &  )
        bHb = beta * beta * vecN(L-j) * &
        &  (bHb + coeff * ((1.0d0 + beta * dble(N) * large / dble(2 * (L-j) + 1)) * vecH(L-j) &
        &                 - beta*vecH2(L-j) / dble(2*(L-j)+1) &
        &                 ) &
        &  )
        bHHb = beta * beta * vecN(L-j) * &
        &  (bHHb + coeff * ((1.0d0 - beta * dble(N) * large / dble(2 * (L-j) + 1)) * vecH2(L-j) &
        &                  + (beta * dble(N*N) * large * large / dble(2 * (L-j) + 1) &
        &                    - dble(2 * (L-j)) / beta &
        &                    ) * vecH(L-j) &
        &                  ) &
        &  )
        !
        if(bb > 1.0d+10) then
           bb = bb * (1.0d-10)
           bHb = bHb * (1.0d-10)
           bHHb = bHHb * (1.0d-10)
           coeff = coeff * (1.0d-10)
        endif
        !
     enddo
     !
     write(*,*)beta,bHb/bb,bHHb/bb
     !
  enddo
  !
end program proc
