! SUBROUTINE fortran_bessel_jn(n1, n2, x, y)
!   implicit none
!   integer(4) :: n1,n2
!   real(8) :: x
!   real(8) :: y(n2-n1+1)
!
!   y = bessel_jn(n1,n2,x)
!
! END SUBROUTINE

SUBROUTINE fortran_bessel_jn(n1, n2, x, y)
  USE Complex_Bessel
  IMPLICIT NONE

  INTEGER       :: np, ierrj, i
  INTEGER       :: n1, n2
  INTEGER       :: kode = 1, nzj = 1
  COMPLEX (dp)  :: z, cj(5), cj_1(4)
  REAL (dp)     :: fnu, x
  REAL (dp)     :: y(1:5)

  fnu = n1
  np = 5

  z = CMPLX(x, 0.0, KIND=dp)

IF (fnu >= 0) THEN
  CALL cbesj(z, fnu, kode, np, cj, nzj, ierrj)
ELSE
  fnu = 0
  np = 4
  CALL cbesj(z, fnu, kode, np, cj_1, nzj, ierrj)
  cj(1) = COMPLEX( -REAL(cj_1(2)), AIMAG(cj_1(2)) )

  DO i = 2, 5
   cj(i) = cj_1(i-1)
  END DO

END IF

DO i = 1, 5
  y(i) = REAL(cj(i))
END DO

END SUBROUTINE

!program test_besjn
!
!  real(8) :: y(3)
!  real(8) :: ff
!  integer(4) :: i,j
!
!  i = 1000000
!  j = 1000002
!  ff = 1d6
!
!  call fortran_bessel_jn(i,j,ff*0.9999,y)
!  print*,y
!
!end program
