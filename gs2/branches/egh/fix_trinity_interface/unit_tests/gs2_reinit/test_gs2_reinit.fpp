

module checks_mod
  use unit_tests
  public checks
  contains
    function checks()
      logical :: checks
      checks = .true.
    end function checks
end module checks_mod

!> A program that tests the gs2_reinit module. It  runs 
!! a  linear cyclone test case and then forces a reinit
!! and checks that the distribution function, fields and
!! response matrix are the same
!!
!! This is free software released under the MIT license
!!   Written by: Edmund Highcock (edmundhighcock@users.sourceforge.net)

program test_gs2_reinit
  !use functional_tests
  !use checks_mod
  !call test_gs2('Linear CBC (unit test) to test new diagnostics', checks)
    use gs2_reinit, only: reset_time_step
    use gs2_init, only: save_fields_and_dist_fn, init
    use gs2_reinit, only: gs2_reinit_unit_test_set_in_memory
    use gs2_main, only: run_gs2, finish_gs2
    use gs2_init, only:  init, init_level_list
    use gs2_main, only: gs2_main_unit_test_reset_gs2
    use gs2_main, only: finalize_diagnostics, initialize_diagnostics
    use gs2_time, only: code_dt_cfl
    use gs2_main, only: old_iface_state
    use fields, only: finish_fields, init_fields
    use fields_local, only: init_fields_local, finish_fields_local
    use unit_tests
    use unit_tests, only: should_print
    use mp, only: init_mp, mp_comm, proc0, test_driver_flag, finish_mp
    use mp, only: broadcast, barrier, iproc
    use fields_arrays, only: phinew, aparnew, bparnew
    use dist_fn_arrays, only: gnew
    use gs2_diagnostics, only: finish_gs2_diagnostics
  use kt_grids, only: naky, ntheta0, init_kt_grids
  use theta_grid, only: ntgrid, init_theta_grid
  use gs2_layouts, only: init_gs2_layouts, g_lo, ie_idx
  use fields_local, only: fields_local_functional, fieldmat
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: finish_gs2_diagnostics_new
    use gs2_diagnostics_new, only: gnostics
#endif
    implicit none
    integer :: n_vars, n_file_names
    integer :: i, it, ik
    real :: eps
    complex, dimension (:,:,:), allocatable :: gbak 
    complex, dimension (:,:,:), allocatable :: phi_bak, apar_bak, bpar_bak
    complex, dimension (:,:), allocatable :: rowbloc
    character(len=29) :: message
    logical :: test_result
    logical :: dummy
    integer :: supercell_idx



    allocate(gbak(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate(phi_bak(-ntgrid:ntgrid,ntheta0,naky))
    allocate(apar_bak(-ntgrid:ntgrid,ntheta0,naky))
    allocate(bpar_bak(-ntgrid:ntgrid,ntheta0,naky))

      !variables = (/'lambda', 'phi'/), &
      !n_lines = (/'3', '30'/)
  ! General config


  eps = 1.0e-7

  if (precision(eps).lt. 11) eps = eps * 100.0


    call init_mp


    test_driver_flag = .true.
    functional_test_flag = .true.

   call announce_module_test("gs2_reinit")

    call run_gs2(mp_comm)

    call announce_test('init down and up')
    call finalize_diagnostics(old_iface_state)
    !call finish_fields_local
    !call init_fields_local
    call finish_fields
    call init_fields
    call init(old_iface_state%init, init_level_list%collisions)
    call init(old_iface_state%init, init_level_list%full)
    call init(old_iface_state%init, init_level_list%gs2)
    call init(old_iface_state%init, init_level_list%full)
    call initialize_diagnostics(old_iface_state)
    call process_test(.true., 'init down and up')

    supercell_idx = min(iproc+2, size(fieldmat%kyb(2)%supercells))
    allocate(rowbloc( &
      size(fieldmat%kyb(2)%supercells(supercell_idx)%cells(1)%rb(1)%data, 1), &
      size(fieldmat%kyb(2)%supercells(supercell_idx)%cells(1)%rb(1)%data, 2)))


    phi_bak = phinew
    apar_bak = aparnew
    bpar_bak = bparnew
    gbak = gnew
    rowbloc = fieldmat%kyb(2)%supercells(supercell_idx)%cells(1)%rb(1)%data

    call save_fields_and_dist_fn
    !call reinit_gk_and_field_equations(reset_antenna=.true.)
    call init(old_iface_state%init, init_level_list%override_timestep)
    call init(old_iface_state%init, init_level_list%full)

    call announce_test('Values of fields and dist fn after reinitialising')
    call process_test(test_fields_and_dist(), 'Values of fields and dist fn after reinitialising')

    call announce_test('Value of response matrix after reinitialising')
    call process_test(response_unchanged(), &
      'Value of response matrix after reinitialising')


    call gs2_reinit_unit_test_set_in_memory(.true.)

    call save_fields_and_dist_fn
    call init(old_iface_state%init, init_level_list%override_timestep)
    call init(old_iface_state%init, init_level_list%full)
    !call reinit_gk_and_field_equations(reset_antenna=.true.) 

    call announce_test('Values of fields and dist fn after reinitialising in memory')
    call process_test(test_fields_and_dist(), 'Values of fields and dist fn after reinitialising in memory')

    call announce_test('Value of response matrix after reinitialising in memory')
    call process_test(response_unchanged(),&
      'Value of response matrix after reinitialising in memory')

    call reset_time_step(old_iface_state%init, 0, dummy)

    call announce_test('Values of fields and dist fn after calling reset_time_step')
    call process_test(test_fields_and_dist(), &
      'Values of fields and dist fn after calling reset_time_step')

    call announce_test('Value of response matrix after calling reset_time_step')
    call process_test(response_unchanged(),&
      'Value of response matrix after calling reset_time_step')

    code_dt_cfl = 0.001
    call reset_time_step(old_iface_state%init, 0, dummy)

    call announce_test('Values of fields and dist fn after changing timestep')
    call process_test(test_fields_and_dist(), 'Values of fields and dist fn after changing timestep')

    call announce_test('Value of response matrix should have changed after changing timestep')
    call process_test(.not. response_unchanged() ,&
      'Value of response matrix should have changed after changing timestep')

    code_dt_cfl = 1.0
    call reset_time_step(old_iface_state%init, 0, dummy)
    call announce_test('Values of fields and dist fn after restoring timestep')
    call process_test(test_fields_and_dist(), 'Values of fields and dist fn after restoring timestep')

    call announce_test('Value of response matrix after restoring timestep')
    call process_test(response_unchanged() ,&
      'Value of response matrix after restoring timestep')

    dummy = gs2_main_unit_test_reset_gs2(1.0)
    call announce_test('Values of fields and dist fn after calling gs2_reset')
    call process_test(test_fields_and_dist(), 'Values of fields and dist fn after calling gs2_reset')

    call announce_test('Value of response matrix  after calling gs2_reset')
    call process_test(response_unchanged() ,&
      'Value of response matrix  after calling gs2_reset')


    dummy = gs2_main_unit_test_reset_gs2(1.5)
    call announce_test('Values of fields and dist fn after calling gs2_reset with a factor')
    call process_test(test_fields_and_dist(), 'Values of fields and dist fn after calling gs2_reset with a factor')

    call announce_test('Value of response matrix should have changed after calling gs2_reset with a factor')
    call process_test(.not. response_unchanged() ,&
      'Value of response matrix should have changed after calling gs2_reset with a factor')



    call barrier


    !call announce_test('results')
    !call process_test(test_function(), 'results')

      !call announce_test("average heat flux")
      !call process_test( &
        !"average heat flux")


    call finish_gs2_diagnostics(ilast_step)

#ifdef NEW_DIAG
    call finish_gs2_diagnostics_new
#endif
    call finish_gs2

   call close_module_test("gs2_reinit")

    call finish_mp

    deallocate(gbak)
    deallocate(phi_bak)
    deallocate(apar_bak)
    deallocate(bpar_bak)

contains
  
  function response_unchanged()
    use mp, only: sum_allreduce
    integer :: result_int
    logical :: response_unchanged

    response_unchanged = agrees_with( &
      fieldmat%kyb(2)%supercells(supercell_idx)%cells(1)%rb(1)%data(4,:), rowbloc(4,:), eps)

    ! Response is only unchanged if it is unchanged
    ! on all procs.
    ! Of course it's possible to achieve the following 
    ! with the right mpi OR-reduce command, but I 
    ! am lazy
    if (response_unchanged) then 
      result_int = 0
    else
      result_int = 1
    end if
    call sum_allreduce(result_int)
    if (result_int.eq.0) then
      response_unchanged = .true.
    else 
      response_unchanged = .false.
    end if
  end function response_unchanged
  function test_fields_and_dist()
    use theta_grid, only: ntheta
    logical :: test_fields_and_dist
    logical :: test_result
    integer :: ig, isgn
    test_result = .true.

    do ik = 1,naky
      do it = 1,ntheta0
        if (it==1 .and. ik==1) cycle
        write(message, fmt="(A19, I2, A6, I2)") 'value of phi,  it =', it, ' ik = ', ik
        call announce_check(message)
        call process_check(test_result, agrees_with(phinew(:, it, ik), phi_bak(:, it, ik), eps), message)
        write(message, fmt="(A19, I2, A6, I2)") 'value of apar, it =', it, ' ik = ', ik
        call announce_check(message)
        call process_check(test_result, agrees_with(aparnew(:, it, ik), apar_bak(:, it, ik), eps), message)
        write(message, fmt="(A19, I2, A6, I2)") 'value of bpar, it =', it, ' ik = ', ik
        call announce_check(message)
        call process_check(test_result, agrees_with(bparnew(:, it, ik), bpar_bak(:, it, ik), eps), message)
        !check_result =  agrees_with(phi_imp(ik, it, :), phi_loc(ik, it, :), eps) .and. check_result
        !if (check_result) write (*,*) it,ik
      end do
    end do
    do ig = -ntgrid,ntgrid
      do isgn = 1,2
        write(message, fmt="(A19, I2, A6, I2)") 'value of gnew, ig =', ig, ' isgn=', isgn
        call announce_check(message)
        call process_check(test_result, agrees_with(gnew(ig, isgn, :), gbak(ig, isgn,:), eps), message)
      end do
    end do

    test_fields_and_dist = test_result
  end function test_fields_and_dist


end program test_gs2_reinit
