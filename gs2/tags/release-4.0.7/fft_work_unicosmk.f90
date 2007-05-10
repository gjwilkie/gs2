module fft_work
  ! ccfft
  !   UNICOS:
  !    table(100 + 8*n)
  !    work(8*n)
  !   UNICOS/mk:
  !    table(2*n)
  !    work(4*n)
  !   Origin:
  !    table(30 + 2*n)
  !    work(2*n)

  integer, parameter :: ccfft_table0 = 0
  integer, parameter :: ccfft_table1 = 2
  integer, parameter :: ccfft_work0 = 0
  integer, parameter :: ccfft_work1 = 4

  ! csfft/scfft
  !   UNICOS:
  !    table(100 + 4*n)
  !    work(4 + 4*n)
  !   UNICOS/mk:
  !    table(2*n)
  !    work(2*n)
  !   Origin:
  !    table(15 + n)
  !    work(n)

  integer, parameter :: csfft_table0 = 0
  integer, parameter :: csfft_table1 = 2
  integer, parameter :: csfft_work0 = 0
  integer, parameter :: csfft_work1 = 2

end module fft_work
