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
MODULE komega_math
  !
  IMPLICIT NONE
  !
  INTERFACE
     !
     REAL(8) FUNCTION ddot(n,dx,incx,dy,incy)
       REAL(8) dx(*),dy(*)
       INTEGER          incx,incy,n
     END FUNCTION ddot
     !
     COMPLEX(8) FUNCTION zdotc(n,zx,incx,zy,incy)
       COMPLEX(8) zx(*),zy(*)
       INTEGER        incx,incy,n
     END FUNCTION zdotc
     !
     COMPLEX(8) FUNCTION zdotu(n,zx,incx,zy,incy)
       COMPLEX(8) zx(*),zy(*)
       INTEGER        incx,incy,n
     END FUNCTION zdotu
     !
     SUBROUTINE  dscal(n,da,dx,incx)
       REAL(8) da,dx(*)
       INTEGER          incx,n
     END SUBROUTINE dscal
     !
     SUBROUTINE  zscal(n,za,zx,incx)
       COMPLEX(8) za,zx(*)
       INTEGER        incx,n
     END SUBROUTINE zscal
     !
     SUBROUTINE  dcopy(n,dx,incx,dy,incy)
       REAL(8) dx(*),dy(*)
       INTEGER          incx,incy,n
     END SUBROUTINE dcopy
     !
     SUBROUTINE  zcopy(n,zx,incx,zy,incy)
       COMPLEX(8) zx(*),zy(*)
       INTEGER        incx,incy,n
     END SUBROUTINE zcopy
     !
     SUBROUTINE daxpy(n,da,dx,incx,dy,incy)
       REAL(8) dx(*),dy(*),da
       INTEGER          incx,incy,n
     END SUBROUTINE daxpy
     !
     SUBROUTINE zaxpy(n,za,zx,incx,zy,incy)
       COMPLEX(8) zx(*),zy(*),za
       INTEGER        incx,incy,n
     END SUBROUTINE zaxpy
     !
  END INTERFACE
  !
CONTAINS
!
! ddot with MPI allreduce
!
FUNCTION ddotMPI(n,dx,dy) RESULT(prod)
  !
#if defined(MPI)
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
  USE komega_parameter, ONLY : comm, lmpi
#endif
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: n
  REAL(8),INTENT(IN) :: dx(n), dy(n)
  REAL(8) prod
  !
#if defined(MPI)
  INTEGER :: ierr
#endif
  !
  prod = ddot(n,dx,1,dy,1)
  !
#if defined(MPI)
  IF(lmpi) &
  &  call MPI_allREDUCE(MPI_IN_PLACE, prod, 1, &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
#endif
  !
END FUNCTION ddotMPI
!
! zdotc with MPI allreduce
!
FUNCTION zdotcMPI(n,zx,zy) RESULT(prod)
  !
#if defined(MPI)
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_COMPLEX, MPI_SUM
  USE komega_parameter, ONLY : comm, lmpi
#endif
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: n
  COMPLEX(8),INTENT(IN) :: zx(n), zy(n)
  COMPLEX(8) prod
  !
#if defined(MPI)
  INTEGER :: ierr
#endif
  !
#if defined(__NO_ZDOT)
  prod = DOT_PRODUCT(zx,zy)
#else
  prod = zdotc(n,zx,1,zy,1)
#endif
  !
#if defined(MPI)
  IF(lmpi) &
  &  call MPI_allREDUCE(MPI_IN_PLACE, prod, 1, &
  &                  MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
#endif
  !
END FUNCTION zdotcMPI
!
! zdotu with MPI allreduce
!
FUNCTION zdotuMPI(n,zx,zy) RESULT(prod)
  !
#if defined(MPI)
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_COMPLEX, MPI_SUM
  USE komega_parameter, ONLY : comm, lmpi
#endif
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: n
  COMPLEX(8),INTENT(IN) :: zx(n), zy(n)
  COMPLEX(8) prod
  !
#if defined(MPI)
  INTEGER :: ierr
#endif
  !
#if defined(__NO_ZDOT)
  prod = SUM(zx(1:n) * zy(1:n))
#else
  prod = zdotu(n,zx,1,zy,1)
#endif
  !
#if defined(MPI)
  IF(lmpi) &
  &  call MPI_allREDUCE(MPI_IN_PLACE, prod, 1, &
  &                  MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
#endif
  !
END FUNCTION zdotuMPI
!
! MAXVAL with MPI allreduce (for real(8))
!
FUNCTION dabsmax(array, n) RESULT(maxarray)
  !
#if defined(MPI)
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_MAX
  USE komega_parameter, ONLY : comm, lmpi
#endif
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: n
  REAL(8),INTENT(IN) :: array(n)
  REAL(8) maxarray
  !
#if defined(MPI)
  INTEGER :: ierr
#endif
  !
  maxarray = MAXVAL(ABS(array))
  !
#if defined(MPI)
  IF(lmpi) &
  &  call MPI_allREDUCE(MPI_IN_PLACE, maxarray, 1, &
  &                  MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
#endif
  !
END FUNCTION dabsmax
!
end MODULE komega_math
