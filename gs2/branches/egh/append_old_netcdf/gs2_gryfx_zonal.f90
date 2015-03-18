module gs2_gryfx_zonal

  use gs2_main, only: gs2_program_state_type

  implicit none
  private
  type(gs2_program_state_type) :: state

  public :: init_gs2_gryfx, advance_gs2_gryfx, finish_gs2_gryfx
  public :: allocate_gryfx_zonal_arrays
  public :: deallocate_gryfx_zonal_arrays

contains
  subroutine allocate_gryfx_zonal_arrays
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use nonlinear_terms, only: gryfx_zonal
    
    implicit none

    allocate(gryfx_zonal%NLdens_ky0(naky*ntheta0*2*ntgrid*nspec))
    allocate(gryfx_zonal%NLupar_ky0(naky*ntheta0*2*ntgrid*nspec))
    allocate(gryfx_zonal%NLtpar_ky0(naky*ntheta0*2*ntgrid*nspec))
    allocate(gryfx_zonal%NLtprp_ky0(naky*ntheta0*2*ntgrid*nspec))
    allocate(gryfx_zonal%NLqpar_ky0(naky*ntheta0*2*ntgrid*nspec))
    allocate(gryfx_zonal%NLqprp_ky0(naky*ntheta0*2*ntgrid*nspec))
  end subroutine allocate_gryfx_zonal_arrays    

  subroutine deallocate_gryfx_zonal_arrays
    use nonlinear_terms, only: gryfx_zonal
    
    implicit none

    deallocate(gryfx_zonal%NLdens_ky0)
    deallocate(gryfx_zonal%NLupar_ky0)
    deallocate(gryfx_zonal%NLtpar_ky0)
    deallocate(gryfx_zonal%NLtprp_ky0)
    deallocate(gryfx_zonal%NLqpar_ky0)
    deallocate(gryfx_zonal%NLqprp_ky0)
  end subroutine deallocate_gryfx_zonal_arrays    

  subroutine init_gs2_gryfx_c(strlen, run_name, mp_comm) &
                            bind(c, name='init_gs2')
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: strlen
    character(kind=c_char), intent(in) :: run_name
    integer(c_int), intent(in) :: mp_comm
    call init_gs2_gryfx(strlen, run_name, mp_comm)

  end subroutine init_gs2_gryfx_c

  subroutine init_gs2_gryfx(strlen, file_name, mp_comm)
    use gs2_main, only: initialize_gs2
    use gs2_main, only: initialize_equations
    use gs2_main, only: initialize_diagnostics
    use nonlinear_terms, only: gryfx_zonal
    use file_utils, only: run_name
    use mp, only: proc0

    implicit none
    integer, intent(in) :: strlen
    character (len=strlen), intent (in) :: file_name
    integer, intent(in) :: mp_comm

    !gryfx_zonal%on = .true.
    state%mp_comm_external = .true.
    state%mp_comm = mp_comm  

    state%run_name_external = .true.
    state%run_name = file_name
 

    call initialize_gs2(state)
    if(proc0) write(*,*) 'initialize_gs2 complete. gs2 thinks run_name is ', run_name
    !set overrides of naky, x0, y0, dt, vnewk here. 
    !currently this is hard-coded in kt_grids.f90 and run_parameters.f90
    !just search for "GRYFX"
    call initialize_equations(state)
    if(proc0) write(*,*) 'initialize_equations complete.'
    call initialize_diagnostics(state)
    if(proc0) write(*,*) 'initialize_diagnostics complete.'
    
    call allocate_gryfx_zonal_arrays

  end subroutine init_gs2_gryfx

  subroutine advance_gs2_gryfx_c(istep, dens_ky0, upar_ky0, tpar_ky0, tprp_ky0, &
                               qpar_ky0, qprp_ky0, phi_ky0, first_half_step) &
                                    bind(c, name='advance_gs2')
    use iso_c_binding
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    implicit none
    complex (c_float_complex), intent (inout) :: & 
                                dens_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                upar_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                tpar_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                tprp_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                qpar_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                qprp_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (out) :: phi_ky0(naky*ntheta0*2*ntgrid)
    logical(c_int), intent (in) :: first_half_step
    integer(c_int), intent (in) :: istep
    call advance_gs2_gryfx(istep, dens_ky0, upar_ky0, tpar_ky0, tprp_ky0, &
                               qpar_ky0, qprp_ky0, phi_ky0, first_half_step)


  end subroutine advance_gs2_gryfx_c

  subroutine advance_gs2_gryfx(istep, dens_ky0, upar_ky0, tpar_ky0, tprp_ky0, &
                               qpar_ky0, qprp_ky0, phi_ky0, first_half_step)

    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use nonlinear_terms, only: gryfx_zonal
    use gs2_main, only: evolve_equations
    use dist_fn, only: getmoms_gryfx
    use mp, only: proc0, broadcast

    implicit none
    complex*8, dimension (naky*ntheta0*2*ntgrid*nspec), intent (inout) :: &
                    dens_ky0, upar_ky0, tpar_ky0, tprp_ky0, qpar_ky0, qprp_ky0
    complex*8, dimension (naky*ntheta0*2*ntgrid), intent (out) :: phi_ky0
    logical, intent (in) :: first_half_step
    integer, intent (in) :: istep

    state%istep_end = istep
    state%print_times = .false.

    ! only proc0 knows values of arrays, so broadcast
    call broadcast(dens_ky0)
    call broadcast(upar_ky0)
    call broadcast(tpar_ky0)
    call broadcast(tprp_ky0)
    call broadcast(qpar_ky0)
    call broadcast(qprp_ky0)
    ! all procs know value of first_half_step already

    ! now that all procs know values, copy into gryfx_zonal structure
    gryfx_zonal%first_half_step = first_half_step 
    gryfx_zonal%NLdens_ky0 = dens_ky0
    gryfx_zonal%NLupar_ky0 = upar_ky0
    gryfx_zonal%NLtpar_ky0 = tpar_ky0
    gryfx_zonal%NLtprp_ky0 = tprp_ky0
    gryfx_zonal%NLqpar_ky0 = qpar_ky0
    gryfx_zonal%NLqprp_ky0 = qprp_ky0

    ! only advance one timestep
    call evolve_equations(state, 1)
    ! getmoms_gryfx takes moments of g and puts results into arguments.
    ! since these are the same arguments that are passed into 
    ! advance_gs2_gryfx, we don't need to copy them like on the way in
    call getmoms_gryfx(dens_ky0, upar_ky0, tpar_ky0, &
                 tprp_ky0, qpar_ky0, qprp_ky0, phi_ky0)

  end subroutine advance_gs2_gryfx

  subroutine finish_gs2_gryfx_c bind(c, name='finish_gs2')

    use iso_c_binding
    implicit none
    call finish_gs2_gryfx

  end subroutine finish_gs2_gryfx_c

  subroutine finish_gs2_gryfx
    use gs2_main, only: finalize_diagnostics
    use gs2_main, only: finalize_equations
    use gs2_main, only: finalize_gs2
    implicit none

    call deallocate_gryfx_zonal_arrays
    call finalize_diagnostics(state)
    call finalize_equations(state)
    call finalize_gs2(state)
  end subroutine finish_gs2_gryfx

  subroutine getmoms_gryfx_c (dens_ky0, upar_ky0, tpar_ky0, tprp_ky0, &
                              qpar_ky0, qprp_ky0, phi_ky0) &
                                   bind(c, name='getmoms_gryfx')
    use dist_fn, only: getmoms_gryfx
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use iso_c_binding
    implicit none
    complex (c_float_complex), intent (inout) :: & 
                                dens_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                upar_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                tpar_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                tprp_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                qpar_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                qprp_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (out) :: phi_ky0(naky*ntheta0*2*ntgrid)
    
    call getmoms_gryfx(dens_ky0, upar_ky0, tpar_ky0, &
            tprp_ky0, qpar_ky0, qprp_ky0, phi_ky0)
  end subroutine getmoms_gryfx_c

  subroutine broadcast_integer_c(a) bind(c, name='broadcast_integer')
    use iso_c_binding
    use mp, only : broadcast
    implicit none
    integer(c_int), intent(in out) :: a

    call broadcast(a) 
  end subroutine broadcast_integer_c
  
end module gs2_gryfx_zonal


