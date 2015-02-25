!> Helper functions for running benchmarks
module benchmarks
  implicit none

  private

  public :: benchmark_identifier
contains
  !> A string which is used as an extension for timing files
  !! and which identifies the time when and system on which the
  !! benchmark was carried out
  function benchmark_identifier()
    use runtime_tests, only: build_identifier
    !use job_manage, only: timer_local
    character(len=61) :: benchmark_identifier
    character(len=11) :: timenowstr
    character(len=8) :: date
    integer :: timenow
    call system_clock(timenow)
    timenow = timenow !/ 1000 ! Convert to seconds
    write(timenowstr, "(I11)") timenow
    call date_and_time(date=date)
    benchmark_identifier = &
      '.timing.'//date//'.'//trim(build_identifier())
  end function benchmark_identifier
end module benchmarks
