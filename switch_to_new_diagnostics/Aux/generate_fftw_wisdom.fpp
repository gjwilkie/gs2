
!> A program that calls init_transforms, to be used for generating fftw wisdom 
!! for a given input file.
!!
!! This is free software released under the MIT licence
!!   Written by: Edmund Highcock (edmundhighcock@users.sourceforge.net)
program generate_fftw_wisdom
  !use unit_tests
  use mp, only: init_mp, finish_mp, proc0, broadcast, nproc
  use file_utils, only: init_file_utils, run_name
  !use file_utils, only: append_output_file, close_output_file
  !use species, only: init_species, nspec, spec
  !use constants, only: pi
  !use kt_grids, only: naky, ntheta0, init_kt_grids
  !use theta_grid, only: ntgrid, init_theta_grid
  !use gs2_layouts, only: init_gs2_layouts, g_lo, ie_idx
  use nonlinear_terms, only: init_nonlinear_terms, finish_nonlinear_terms
  !use dist_fn_arrays, only: g
  !use nonlinear_terms, only: nonlinear_terms_unit_test_time_add_nl
  !use kt_grids, only: ntheta0, naky
  use job_manage, only: time_message
  !use runtime_tests, only: get_svn_rev, get_compiler_name
  use benchmarks, only: benchmark_identifier
  implicit none
    character (500), target :: cbuff
  real :: tstart
  logical :: dummy=.false.
  real :: time_taken(2) = 0.0
  !integer :: i
  !integer :: timing_unit


  ! Set up depenencies
  call init_mp
  if (proc0) call init_file_utils(dummy, name="gs")
       if (proc0) then
          cbuff = trim(run_name)
       end if
       
       call broadcast (cbuff)
       if (.not. proc0) run_name => cbuff





  if (proc0) call time_message(.false., time_taken, "Init transforms")

  call init_nonlinear_terms

  if (proc0) then
    call time_message(.false., time_taken, "Init transforms")
    write(*, '(" Time taken to calculate fftw plans ",I6," procs: ",F3.1," s")') nproc, time_taken(1)
    write(*,*)
!    call append_output_file(timing_unit, &
!      benchmark_identifier())
!    write(timing_unit, '(I6,"   ",F9.3)') nproc, time_taken(1)
!    call close_output_file(timing_unit)
  end if

  !call finish_nonlinear_terms


  call finish_mp

end program generate_fftw_wisdom
