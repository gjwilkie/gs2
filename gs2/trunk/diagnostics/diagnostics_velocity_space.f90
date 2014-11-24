!> A module which writes out quantities which writes out 
!! diagnostic quantities which assess whether the velocity 
!! space resolution is sufficient.
module diagnostics_velocity_space
  implicit none
  private

  public :: init_diagnostics_velocity_space, write_velocity_space_checks

contains
  
  subroutine init_diagnostics_velocity_space(gnostics)
    use le_grids, only: init_weights
    use mp, only: proc0
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics 
    
    if (.not. gnostics%write_verr) return
    
    ! initialize weights for less accurate integrals used
    ! to provide an error estimate for v-space integrals (energy and untrapped)
    if (proc0) call init_weights
  end subroutine init_diagnostics_velocity_space
  
  
  subroutine write_velocity_space_checks(gnostics)
    use dist_fn, only: get_verr, get_gtran
    use mp, only: proc0
    use le_grids, only: nlambda, ng2
    use fields_arrays, only: phinew, bparnew
    use gs2_time, only: user_time
    use collisions, only: vnmult
    use species, only: spec
    use diagnostics_config, only: diagnostics_type
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_dimensions, only: dim_string
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real, dimension (:,:), allocatable :: errest
    integer, dimension (:,:), allocatable :: erridx
    real :: geavg, glavg, gtavg
    
    allocate(errest(5,2), erridx(5,3))
    errest = 0.0; erridx = 0
    
    ! error estimate obtained by comparing standard integral with less-accurate integral
    call get_verr (errest, erridx, phinew, bparnew)
    
    ! error estimate based on monitoring amplitudes of legendre polynomial coefficients
    call get_gtran (geavg, glavg, gtavg, phinew, bparnew)
    
    !errest(5,:) = -4
    
    if (.not. gnostics%vary_vnew_only) then
       call create_and_write_variable(gnostics, gnostics%rtype, "vspace_lpcfrac", &
            dim_string([gnostics%dims%generic_3,gnostics%dims%time]), &
            "Fraction of free energy contained in the high order coefficients of &
            & the Legendre polynomial transform of (1) energy space, (2) untrapped &
            & pitch angles and (3) trapped pitch angles (each should ideally be < 0.1).  &
            & Note that there are no trapped pitch angles for certain geometries", &
            "1", (/geavg, glavg, gtavg/))
       call create_and_write_variable(gnostics, gnostics%rtype, "vspace_err", &
            dim_string([gnostics%dims%generic_5,gnostics%dims%generic_2,gnostics%dims%time]), &
            "Estimate of the (1) absolute and (2) relative errors resulting from &
            & velocity space integrals in the calculation of the following quantities &
            & in the given dimensions: (1) k phi, energy (2) k phi, untrapped pitch angles &
            & (3) k phi, trapped pitch angles, (4) k apar, energy, (5) k apar, untrapped &
            & angles. Relative errors should be < 0.1. ", &
            "absolute error measures have units T_r/(e rho_r)", errest)
       call create_and_write_variable(gnostics, gnostics%rtype, "vspace_vnewk", &
            dim_string([gnostics%dims%generic_2,gnostics%dims%time]), &
            "If the simulation is set to vary the collisionality in order to keep &
            & error in velocity integrals to acceptable levels, contains species 1 &
            & collisionality in (1) pitch angle and (2) energy  ", &
            "v_thr/a", (/vnmult(1)*spec(1)%vnewk, vnmult(2)*spec(1)%vnewk/))
       ! This next statement causes annoying printout because of an error in netcdf 4.1
       ! The error is fixed in 4.2
       ! See
       ! http://www.unidata.ucar.edu/software/netcdf/docs/known_problems.html#f90-debug-segfault
      if (gnostics%write_max_verr) &
           call create_and_write_variable(gnostics, gnostics%itype, "vspace_err_maxindex", &
            dim_string([gnostics%dims%generic_5,gnostics%dims%generic_3,gnostics%dims%time]), &
           "Gives the (1) theta index, (2) ky index and (3) kx index of the maximum &
           & error resulting from the &
           & velocity space integrals in the calculation of the following quantities &
           & in the given dimensions: (1) k phi, energy (2) k phi, untrapped pitch angles &
           & (3) k phi, trapped pitch angles, (4) k apar, energy, (5) k apar, untrapped &
           & angles. Relative errors should be < 0.1. ", &
           "1", erridx)
   end if
   
   if (proc0 .and. gnostics%write_ascii) call write_ascii
   deallocate(errest,erridx)
   
 contains
   subroutine write_ascii
     if (nlambda - ng2 > 1) then
        write(gnostics%ascii_files%lpc,"(4(1x,e13.6))") user_time, geavg, glavg, gtavg
     else
        write(gnostics%ascii_files%lpc,"(3(1x,e13.6))") user_time, geavg, glavg
     end if
     write(gnostics%ascii_files%vres,"(8(1x,e13.6))") user_time, errest(1,2), errest(2,2), errest(3,2), &
          errest(4,2), errest(5,2), vnmult(1)*spec(1)%vnewk, vnmult(2)*spec(1)%vnewk
     if (gnostics%write_max_verr) then
        write(gnostics%ascii_files%vres2,"(3(i8),(1x,e13.6),3(i8),(1x,e13.6),3(i8),(1x,e13.6),3(i8),(1x,e13.6),3(i8),(1x,e13.6))") &
             erridx(1,1), erridx(1,2), erridx(1,3), errest(1,1), &
             erridx(2,1), erridx(2,2), erridx(2,3), errest(2,1), &
             erridx(3,1), erridx(3,2), erridx(3,3), errest(3,1), &
             erridx(4,1), erridx(4,2), erridx(4,3), errest(4,1), &
             erridx(5,1), erridx(5,2), erridx(5,3), errest(5,1)
     end if
   end subroutine write_ascii
 end subroutine write_velocity_space_checks
end module diagnostics_velocity_space
