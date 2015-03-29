
!> This unit tests the gs2 main interface,
!! and by repeatedly initializing and finalizing
!! gs2, tests that gs2 is being properly tidied up,
!! variables deallocated etc.
!
program test_gs2_gryfx_zonal
  use gs2_main
  use gs2_gryfx_zonal
  use unit_tests
  use mp, only: init_mp, finish_mp, mp_comm
  use theta_grid, only: ntgrid, theta
  use kt_grids, only: ntheta0, naky
  use species, only: nspec
  use mp, only: proc0
  use file_utils, only: run_name
  implicit none
  real :: eps
  integer :: i, gs2_counter
  type(gs2_program_state_type) :: gs2_state
  type(gryfx_parameters_type) :: gryfx_parameters
  logical :: first_half_step
  character (1000) :: file_name
  real*8, dimension (:), allocatable :: gryfx_theta
  complex*8, dimension (:), allocatable :: dens_ky0, upar_ky0, tpar_ky0, &
                                         tprp_ky0, qpar_ky0, qprp_ky0, phi_ky0


  eps = 1.0e-7
  if (precision(eps).lt. 11) eps = eps * 1000.0
  
  call init_mp
  
  call announce_module_test("gs2_gryfx_zonal")


  gs2_state%mp_comm_external = .true.
  gs2_state%mp_comm = mp_comm  !this will come from GryfX
  test_flag = .true.

  call announce_test('initializing and finalizing gs2')

  call initialize_gs2(gs2_state)
  call finalize_gs2(gs2_state)

  call initialize_gs2(gs2_state)
  call finalize_gs2(gs2_state)

  call initialize_gs2(gs2_state)
  call finalize_gs2(gs2_state)
  call process_test(.true.,'initializing and finalizing gs2')

  call announce_test('initializing equations and getting run name')

  call initialize_gs2(gs2_state)
  call initialize_equations(gs2_state)
  call finalize_equations(gs2_state)
  call finalize_gs2(gs2_state)

  call initialize_gs2(gs2_state)
  call initialize_equations(gs2_state)
  file_name = trim(run_name) // '.in'  !this will come from GryfX
  if(proc0) write (*,*) 'run_name is ', file_name
  allocate(gryfx_theta(2*ntgrid))
  gryfx_theta(1:2*ntgrid) = theta(-ntgrid:ntgrid-1)
  call finalize_equations(gs2_state)
  call finalize_gs2(gs2_state)
  call process_test(.true.,'initializing equations and getting run name')

  call announce_test('initializing gs2_gryfx_zonal')
  call init_gs2_gryfx(len_trim(file_name), file_name, gs2_state%mp_comm, &
                                gryfx_theta, gryfx_parameters)
  call finish_gs2_gryfx
  call process_test(.true.,'initializing gs2_gryfx_zonal')

  !!program gs2
    !!type(gs2_program_gs2_state_type) :: gs2_state
  call announce_test('initialize and evolve equations')
    call initialize_gs2(gs2_state)
    call initialize_equations(gs2_state)
    call initialize_diagnostics(gs2_state)
    !if (gs2_state%eigsolve) then 
      !call solve_eigenproblem
    !else
    call evolve_equations(gs2_state, gs2_state%nstep/2)
    call evolve_equations(gs2_state, gs2_state%nstep/2)
    !! This call should do nothing and print a warning
    call evolve_equations(gs2_state, gs2_state%nstep/2)
    call evolve_equations(gs2_state, gs2_state%nstep/2)

    call calculate_outputs(gs2_state)


    call finalize_diagnostics(gs2_state)
    call finalize_equations(gs2_state)
    call finalize_gs2(gs2_state)
  !!end program gs2
  call process_test(.true., 'initialize and evolve equations')


  call announce_test('initialize gs2 gryfx II')
  !begin hybrid gs2_gryfx_zonal algorithm
  call init_gs2_gryfx(len_trim(file_name), file_name, gs2_state%mp_comm, &
                                gryfx_theta, gryfx_parameters)
  if(proc0) write (*,*) 'naky = ', naky
  ! dens_ky0, upar_ky0, etc will come from gryfx. 
  ! in this test, we need to allocate and initialize them
  ! let's use dummy values for now, eventually these can
  ! be calculated from an actual hybrid simulation
    allocate(dens_ky0(ntheta0*2*ntgrid*nspec))
    allocate(upar_ky0(ntheta0*2*ntgrid*nspec))
    allocate(tpar_ky0(ntheta0*2*ntgrid*nspec))
    allocate(tprp_ky0(ntheta0*2*ntgrid*nspec))
    allocate(qpar_ky0(ntheta0*2*ntgrid*nspec))
    allocate(qprp_ky0(ntheta0*2*ntgrid*nspec))
    allocate(phi_ky0(ntheta0*2*ntgrid))
     
    ! in GryfX, only proc0 will know values of arrays
    if(proc0) then
      dens_ky0 = 1.e-10
      upar_ky0 = 1.e-10
      tpar_ky0 = 1.e-10
      tprp_ky0 = 1.e-10
      qpar_ky0 = 1.e-10
      qprp_ky0 = 1.e-10
      phi_ky0 = 1.e-10
    endif
  call process_test(.true., 'initialize gs2 gryfx II')

  call announce_test('first half step')
  gs2_counter = 1
  !first_half_step will be set on all procs in GryfX
  first_half_step = .true.
  call advance_gs2_gryfx(gs2_counter, dens_ky0, upar_ky0, tpar_ky0, tprp_ky0, qpar_ky0, &
                        qprp_ky0, phi_ky0, first_half_step)
  call process_test(.true., 'first half step')
  call announce_test('second half step')
  gs2_counter = gs2_counter + 1
  first_half_step = .false.
  call advance_gs2_gryfx(gs2_counter, dens_ky0, upar_ky0, tpar_ky0, tprp_ky0, qpar_ky0, &
                        qprp_ky0, phi_ky0, first_half_step)
  gs2_counter = gs2_counter + 1
  call process_test(.true., 'second half step')
  call announce_test('finish_gs2_gryfx II')
  call finish_gs2_gryfx
  call process_test(.true., 'finish_gs2_gryfx II')


  call init_gs2_gryfx(len_trim(file_name), file_name, gs2_state%mp_comm, &
                                gryfx_theta, gryfx_parameters)
    if(proc0) then
      dens_ky0 = 1.e-10
      upar_ky0 = 1.e-10
      tpar_ky0 = 1.e-10
      tprp_ky0 = 1.e-10
      qpar_ky0 = 1.e-10
      qprp_ky0 = 1.e-10
      phi_ky0 = 1.e-10
    endif
  gs2_counter = 1

  ! if nsteps is set by input file, this will produce a warning for last two
  ! calls to advance. need to overwrite nsteps -> 2*nsteps for use with gryfx
  do i=1,gs2_state%nstep/2 + 1 
  first_half_step = .true.
  call advance_gs2_gryfx(gs2_counter, dens_ky0, upar_ky0, tpar_ky0, tprp_ky0, qpar_ky0, &
                        qprp_ky0, phi_ky0, first_half_step)
  gs2_counter = gs2_counter + 1
    if(proc0) then
      dens_ky0 = 1.e-10
      upar_ky0 = 1.e-10
      tpar_ky0 = 1.e-10
      tprp_ky0 = 1.e-10
      qpar_ky0 = 1.e-10
      qprp_ky0 = 1.e-10
      phi_ky0 = 1.e-10
    endif
  first_half_step = .false.
  call advance_gs2_gryfx(gs2_counter, dens_ky0, upar_ky0, tpar_ky0, tprp_ky0, qpar_ky0, &
                        qprp_ky0, phi_ky0, first_half_step)
  gs2_counter = gs2_counter + 1
    if(proc0) then
      dens_ky0 = 1.e-10
      upar_ky0 = 1.e-10
      tpar_ky0 = 1.e-10
      tprp_ky0 = 1.e-10
      qpar_ky0 = 1.e-10
      qprp_ky0 = 1.e-10
      phi_ky0 = 1.e-10
    endif
  end do

  !write(*,*) gs2_2_gryfx_grid(:,:)

  ! interpolation leaves unchanged because here
  ! the 'gryfx grid' is the same as the gs2 one
  call announce_test(' interpolation leaves unchanged')
  call interpolate_theta(gs2_2_gryfx_grid, dens_ky0, .false.)
  call process_test(agrees_with(real(real(dens_ky0)), real(real(upar_ky0)), eps), &
                     ' interpolation leaves unchanged')

  


  call finish_gs2_gryfx



    deallocate(dens_ky0)
    deallocate(upar_ky0)
    deallocate(tpar_ky0)
    deallocate(tprp_ky0)
    deallocate(qpar_ky0)
    deallocate(qprp_ky0)
    deallocate(phi_ky0)
    !!call init_mp

    !!test_driver_flag = .true.
    !!functional_test_flag = .true.



  call finalize_overrides(gs2_state)
  call close_module_test("gs2_gryfx_zonal")

  !call finish_gs2


  call finish_mp



contains
  


end program test_gs2_gryfx_zonal
