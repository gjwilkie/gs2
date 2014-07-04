
!> A program that runs unit tests on the fields_local module.
!! Currently success is defined as giving the same results as
!! the old implicit module...could be improved by an absolute test
!!
!! This is free software released under GPLv3
!!   Written by: Edmund Highcock (edmundhighcock@sourceforge.net)
program test_fields_local
  use unit_tests
  use fields_local, only: fields_local_unit_test_init_fields_matrixlocal, advance_local, fields_local_functional, minnrow, do_smart_update
  use fields_implicit, only: fields_implicit_unit_test_init_fields_implicit, advance_implicit

  use fields, only: fields_pre_init
  use egrid
  !use general_f0, only: init_general_f0
  use mp, only: init_mp, finish_mp, proc0, broadcast
  use file_utils, only: init_file_utils
  use species, only: init_species, nspec, spec
  use constants, only: pi
  !use fields, only: init_fields
  !use fields_arrays, only: phi
  use dist_fn, only: init_dist_fn
  !use dist_fn_arrays, only: g
  use kt_grids, only: naky, ntheta0, init_kt_grids
  use theta_grid, only: ntgrid, init_theta_grid
  use gs2_layouts, only: init_gs2_layouts, g_lo, ie_idx
  implicit none
  real :: eps
  logical :: dummy


  ! General config
  eps = 1.0e-7

  if (precision(eps).lt. 11) eps = eps * 100.0

  ! Set up depenencies
  call init_mp
  if (proc0) call init_file_utils(dummy)
  call init_species
  !call init_general_f0



  call announce_module_test('fields_local')

  !call init_kt_grids
  !call init_theta_grid
  !call init_dist_fn
  call fields_pre_init
  call broadcast(MinNrow)
  call broadcast(do_smart_update)
  !call init_fields_implicit

  call announce_test('init_fields_implicit')
  call process_test(fields_implicit_unit_test_init_fields_implicit(), 'init_fields_implicit')

  if (fields_local_functional()) then

    call announce_test('init_fields_matrixlocal')
    call process_test(fields_local_unit_test_init_fields_matrixlocal(), 'init_fields_matrixlocal')

    call announce_test('advance')
    call process_test(test_advance(eps), 'advance')

  else 

    write (*,*) "WARNING: fields_local is non-functional in your build. &
      & Skipping the fields_local unit test. &
      & If you are using the PGI compilers this is to be expected. "
  end if 


  call close_module_test('fields_local')
  call finish_mp
contains

  function test_advance(eps)
    use init_g, only: ginit
    use dist_fn_arrays, only: g, gnew
    use dist_fn, only: get_init_field
    use fields_arrays
    use run_parameters, only: nstep
    use fields, only: remove_zonal_flows_switch
    complex, dimension (:,:,:), allocatable :: gbak 
    complex, dimension (:,:,:), allocatable :: phi_imp, apar_imp, bpar_imp
    complex, dimension (:,:,:), allocatable :: phi_loc, apar_loc, bpar_loc
    character(len=29) :: message
    real, intent(in) :: eps
    logical :: test_advance
    logical :: check_result
    logical :: restarted
    integer :: istep, ik, it

    allocate(gbak(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))
    allocate(phi_imp(-ntgrid:ntgrid,ntheta0,naky))
    allocate(apar_imp(-ntgrid:ntgrid,ntheta0,naky))
    allocate(bpar_imp(-ntgrid:ntgrid,ntheta0,naky))
    allocate(phi_loc(-ntgrid:ntgrid,ntheta0,naky))
    allocate(apar_loc(-ntgrid:ntgrid,ntheta0,naky))
    allocate(bpar_loc(-ntgrid:ntgrid,ntheta0,naky))

    test_advance = .true.
    !Now we want to fill g with some data, for now just use initialisation routines
    call ginit(restarted)

    !Backup this initial g
    gbak=g

    !//USE 'IMPLICIT'

    !Now we setup the initial fields to be consistent with g
    !gnew =0.0
    call get_init_field(phinew,aparnew,bparnew)
    phi = phinew; apar = aparnew; bpar = bparnew; g = gnew

    !Now we can do a timestep (or lots)
    do istep=1,nstep
        call advance_implicit(istep, remove_zonal_flows_switch)
    enddo

    !!Now we store the results
    phi_imp=phinew
    apar_imp=aparnew
    bpar_imp=bparnew

    !//USE 'LOCAL'
    !Restore original g
    g=gbak
    gnew=g

    !Now we setup the initial fields to be consistent with g
    call get_init_field(phinew,aparnew,bparnew)
    phi = phinew; apar = aparnew; bpar = bparnew; g=gnew

    !Now we can do a timestep (or lots)
    do istep=1,nstep
        call advance_local(istep, remove_zonal_flows_switch)
    enddo

    !Now we store the results
    phi_loc=phinew
    apar_loc=aparnew
    bpar_loc=bparnew

    !ig = maxexponent(eps)/4

    !write (*,*) 'max exponent', maxexponent(eps), 10.0**(-real(maxexponent(eps)-100)), ig, 10.0**(-ig), 1.0e-300 .lt. 10.0**(-ig)
    do ik = 1,naky
      do it = 1,ntheta0
        if (it==1 .and. ik==1) cycle
        write(message, fmt="(A19, I2, A6, I2)") 'value of phi,  it =', it, ' ik = ', ik
        call announce_check(message)
        call process_check(test_advance, agrees_with(phi_imp(:, it, ik), phi_loc(:, it, ik), eps), message)
        write(message, fmt="(A19, I2, A6, I2)") 'value of apar, it =', it, ' ik = ', ik
        call announce_check(message)
        call process_check(test_advance, agrees_with(apar_imp(:, it, ik), apar_loc(:, it, ik), eps), message)
        write(message, fmt="(A19, I2, A6, I2)") 'value of bpar, it =', it, ' ik = ', ik
        call announce_check(message)
        call process_check(test_advance, agrees_with(bpar_imp(:, it, ik), bpar_loc(:, it, ik), eps), message)
        !check_result =  agrees_with(phi_imp(ik, it, :), phi_loc(ik, it, :), eps) .and. check_result
        !if (check_result) write (*,*) it,ik
      end do
    end do

  end function test_advance

end program test_fields_local
