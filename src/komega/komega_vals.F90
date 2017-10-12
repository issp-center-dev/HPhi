!
! ISSP Math Library - A library for solving linear systems in materials science
! Copyright (C) 2016 Mitsuaki Kawamura
! 
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
! 
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
! 
! For more details, See `COPYING.LESSER' in the root directory of this library.
!
MODULE komega_parameter
  !
  IMPLICIT NONE
  !
  REAL(8),PARAMETER :: &
  & almost0 = 1d-50
  !
  INTEGER,SAVE :: &
  & comm,    & !< Communicator for MPI
  & iz_seed, & !< Index of frequency of seed
  & ndim,    & !< Dimension of Hamiltonian
  & nl,      & !< Dimension of projection
  & nz,      & !< Number of Frequency (Shift)
  & itermax, & !< Maximum number of iteration
  & iter       !< Counter of iteration
  !
  LOGICAL,SAVE :: &
  & lmpi !< Use MPI or not
  !
  REAL(8),SAVE :: &
  & threshold, & !< Convergence threshold
  & resnorm      !< Residual norm
  !
  LOGICAL,ALLOCATABLE,SAVE :: &
  & lz_conv(:) !< If converged at this frequency -> .TRUE.
  !
#if defined(__KOMEGA_THREAD)
  !$OMP THREADPRIVATE(comm, iz_seed, ndim, nl, nz, itermax, &
  !$OMP &             iter, threshold, resnorm, lz_conv)
#endif
  !
END MODULE komega_parameter
!
!> Variables for CG
!! \f$\alpha\f$
!
MODULE komega_vals_r
  !
  IMPLICIT NONE
  !
  REAL(8),SAVE :: z_seed    !< Seed frequency
  REAL(8),SAVE :: rho       !< \f$\rho\f$ of BiCG
  REAL(8),SAVE :: alpha     !< \f$\alpha\f$ of BiCG
  REAL(8),SAVE :: alpha_old !< \f$\alpha\f$ at the previous step
  REAL(8),SAVE :: beta      !< \f$\beta\f$ of BiCG
  !
  REAL(8),ALLOCATABLE,SAVE :: z(:)          !< Frequencies
  REAL(8),ALLOCATABLE,SAVE :: pi(:)         !< \f$\pi\f$ of BiCG
  REAL(8),ALLOCATABLE,SAVE :: pi_old(:)     !< \f$\pi\f$ at the previous step
  REAL(8),ALLOCATABLE,SAVE :: pi_save(:,:)  !< \f$\pi\f$ saved for the restart
  REAL(8),ALLOCATABLE,SAVE :: alpha_save(:) !< \f$\alpha\f$ saved for restart
  REAL(8),ALLOCATABLE,SAVE :: beta_save(:)  !< \f$\beta\f$ saved for restart
  !
#if defined(__KOMEGA_THREAD)
  !$OMP THREADPRIVATE(z_seed, rho, alpha, alpha_old, beta, z, pi, &
  !$OMP &             pi_old, pi_save, alpha_save, beta_save)
#endif
  !
END MODULE komega_vals_r
!
!
!
MODULE komega_vecs_r
  !
  IMPLICIT NONE
  !
  REAL(8),ALLOCATABLE,SAVE :: &
  & v3(:),        &
  & p(:,:),       &
  & r_l_save(:,:)
  !
#if defined(__KOMEGA_THREAD)
  !$OMP THREADPRIVATE(v3, p, r_l_save)
#endif
  !
END MODULE komega_vecs_r
!
!
!
MODULE komega_vals_c
  !
  IMPLICIT NONE
  !
  COMPLEX(8),SAVE :: &
  & z_seed, &
  & rho, &
  & alpha, &
  & alpha_old, &
  & beta
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & z(:), &
  & pi(:), &
  & pi_old(:), &
  & pi_save(:,:), &
  & alpha_save(:), &
  & beta_save(:)
  !
#if defined(__KOMEGA_THREAD)
  !$OMP THREADPRIVATE(z_seed, rho, alpha, alpha_old, beta, z, pi, pi_old, &
  !$OMP &             pi_save, alpha_save, beta_save)
#endif
  !
END MODULE komega_vals_c
!
!
!
MODULE komega_vecs_c
  !
  IMPLICIT NONE
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & v3(:), &
  & v5(:), &
  & p(:,:), &
  & r_l_save(:,:)
  !
#if defined(__KOMEGA_THREAD)
  !$OMP THREADPRIVATE(v3, v5, p, r_l_save)
#endif
  !
END MODULE komega_vecs_c
