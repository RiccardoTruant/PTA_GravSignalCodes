program test_besjn
 real(8) :: x = 1.0_8
 x = bessel_jn(100000,100000*0.999999)
 print*,x
end program test_besjn
