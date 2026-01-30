program test_amos_bes

  implicit none
  external :: ZBESJ

  integer(4) :: err
  integer(4) :: NZ
  integer(4) :: KODE = 1
  integer(4) :: N2 = 1
  ! integer(4) :: n1 = 1
  ! integer(4) :: n2 = 5

  real(8) :: n1 = 1.0D0
  real(8) :: Re = 0.0, Im = 0.0

  real(8) :: r1 = 0.0D0
  real(8) :: r2 = 0.0D0

  real(8) :: y = 0.0D0


  call ZBESJ(Re,Im,n1,KODE,N2,r1,r2,NZ,err)
  print*,r1, r2
	print*,"err = ",err


  y = bessel_jn(1,Re)

  print*, y

end program test_amos_bes
