program check
  character(200) :: a
  integer :: n, istatus
  complex, dimension (:), allocatable :: amh, amo
  real :: amax, amax1

  open (unit=1, file="tmp.out.hydra", status="old")
  open (unit=2, file="tmp.out", status="old")

  read (1,*) n
  read (2,*) n
  print *, n

  allocate (amh(n), amo(n))

  amax = 0.0

  do
     istatus = 0
     read (1,*,iostat=istatus) a
     if (istatus /= 0) exit
     read (2,"(a200)") a
     read (1,"(6e12.6)") amh
     read (2,"(6e12.6)") amo
     amax1 = maxval(abs(2.0*(amh-amo))/(abs(amh+amo)+epsilon(0.0)))
     print *, trim(a), " max: ", amax1
     print "(6(x,e12.6))", abs(2.0*(amh-amo))/(abs(amh+amo)+epsilon(0.0))
     amax = max(amax,amax1)
  end do

  print *, "max difference: ", amax

end program check
