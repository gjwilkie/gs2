! lowflow terms include higher-order corrections to GK equation
! such as parallel nonlinearity that require derivatives in v-space.
! most efficient way to take these derivatives is to go from g_lo to le_lo,
! i.e., bring all energies and lambdas onto each processor
# ifdef LOWFLOW
# ifndef USE_LE_LAYOUT
# define USE_LE_LAYOUT on
# endif
# endif

module le_derivatives

  implicit none

  public :: vspace_derivatives

  private

contains

  subroutine vspace_derivatives (g, gold, g1, phi, apar, bpar, phinew, aparnew, bparnew, diagnostics)

    use redistribute, only: gather, scatter
    use dist_fn_arrays, only: c_rate, g_adjust
    use gs2_layouts, only: g_lo, le_lo
    use gs2_time, only: code_dt
    use run_parameters, only: fphi, fbpar
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment, nxi, negrid, g2le
    use collisions, only: solfp1, colls, adjust, heating, hyper_colls
    ! TMP FOR TESTING -- MAB
    use mp, only: proc0

    implicit none
    
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, gold, g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar, phinew, aparnew, bparnew
    integer, optional, intent (in) :: diagnostics

    complex, dimension (:,:,:), allocatable :: gle
    complex, dimension (:,:,:), allocatable :: gc1, gc2, gc3
    logical :: heating_flag

    heating_flag = heating .and. present(diagnostics)

# ifdef USE_LE_LAYOUT

# ifndef LOWFLOW
    if (colls) then
# endif

       if (adjust) then
          call g_adjust (g, phinew, bparnew, fphi, fbpar)
          if (heating_flag) call g_adjust (gold, phi, bpar, fphi, fbpar)
       end if
       if (heating_flag) then
          allocate (gc3(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          gc3 = g
       end if
       
       allocate (gle(nxi+1,negrid+1,le_lo%llim_proc:le_lo%ulim_alloc)) ; gle = 0.
       
       ! map data from g_layout to le_layout
       call gather (g2le, g, gle)
       
# ifdef LOWFLOW
       ! add lowflow stuff here

       if (colls) then
# endif
       ! update distribution function to take into account collisions
       if (present(diagnostics)) then
          call solfp1 (gle, diagnostics)
       else
          call solfp1 (gle)
       end if

# ifdef LOWFLOW
       end if
# endif

       ! remap from le_layout to g_layout
       call scatter (g2le, gle, g)

       deallocate (gle)

       if (heating_flag) then
          ! form (h_i+1 + h_i)/2 * C(h_i+1) and integrate.  
          gc3 = 0.5*conjg(g+gold)*(g-gc3)/code_dt
          
          call integrate_moment (gc3, c_rate(:,:,:,:,3))
          
          deallocate (gc3)
       end if

       if (adjust) then
          call g_adjust (g, phinew, bparnew, -fphi, -fbpar)
          if (heating_flag) call g_adjust (gold, phi, bpar, -fphi, -fbpar)
       end if

# ifndef LOWFLOW
    end if
# endif LOWFLOW

# else

    if (colls) then

       if (heating_flag) then
          allocate (gc1(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          if (hyper_colls) then
             allocate (gc2(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          else
             allocate (gc2(1,1,1))
          end if
          allocate (gc3(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       else
          allocate (gc1(1,1,1))
          allocate (gc2(1,1,1))
          allocate (gc3(1,1,1))
       end if
       gc1 = 0. ; gc2 = 0. ; gc3 = 0.

       if (adjust) then
          call g_adjust (g, phinew, bparnew, fphi, fbpar)
          if (heating_flag) call g_adjust (gold, phi, bpar, fphi, fbpar)
       end if

       if (heating_flag) gc3 = g

       ! update distribution function to take into account collisions
       if (present(diagnostics)) then
          call solfp1 (g, g1, gc1, gc2, diagnostics)
       else
          call solfp1 (g, g1, gc1, gc2)
       end if

       if (heating_flag) then
          call integrate_moment (gc1, c_rate(:,:,:,:,1))
          deallocate (gc1)
          
          if (hyper_colls) call integrate_moment (gc2, c_rate(:,:,:,:,2))
          deallocate (gc2)
          
          ! form (h_i+1 + h_i)/2 * C(h_i+1) and integrate.  
          gc3 = 0.5*conjg(g+gold)*(g-gc3)/code_dt
          
          call integrate_moment (gc3, c_rate(:,:,:,:,3))
          
          deallocate (gc3)
       end if

       if (adjust) then
          call g_adjust (g, phinew, bparnew, -fphi, -fbpar)
          if (heating_flag) call g_adjust (gold, phi, bpar, -fphi, -fbpar)
       end if
    end if

# endif

  end subroutine vspace_derivatives

end module le_derivatives
