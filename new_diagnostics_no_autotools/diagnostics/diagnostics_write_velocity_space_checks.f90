!> A module which writes out quantities which writes out 
!! diagnostic quantities which assess whether the velocity 
!! space resolution is sufficient.
module diagnostics_write_velocity_space_checks
  
  use diagnostics_create_and_write
  use simpledataio
  use diagnostics_config, only: diagnostics_type
  public :: init_diagnostics_write_velocity_space_checks

contains

  subroutine init_diagnostics_write_velocity_space_checks(gnostics)
    use le_grids, only: init_weights
    use mp, only: proc0
    type(diagnostics_type), intent(in) :: gnostics 
    

    if (.not. gnostics%write_verr) return

    ! initialize weights for less accurate integrals used
    ! to provide an error estimate for v-space integrals (energy and untrapped)
    if (proc0) call init_weights
  end subroutine init_diagnostics_write_velocity_space_checks

  
  subroutine write_velocity_space_checks(gnostics)
    use dist_fn, only: get_verr, get_gtran
    use mp, only: proc0
    use le_grids, only: nlambda, ng2
    use fields_arrays, only: phinew, bparnew
    use gs2_time, only: user_time
    use collisions, only: vnmult
    use species, only: spec
    use diagnostics_config, only: diagnostics_type
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

    call create_and_write_variable(gnostics, gnostics%rtype, "vspace_lpcfrac", "3t", &
      "Fraction of free energy contained in the high order coefficients of &
        & the Legendre polynomial transform of (1) energy space, (2) untrapped &
        & pitch angles and (3) trapped pitch angles (each should ideally be < 0.1).  &
        & Note that there are no trapped pitch angles for certain geometries", &
        "1", (/geavg, glavg, gtavg/))
    call create_and_write_variable(gnostics, gnostics%rtype, "vspace_err", "52t", &
      "Estimate of the (1) absolute and (2) relative errors resulting from &
      & velocity space integrals in the calculation of the following quantities &
      & in the given dimensions: (1) k phi, energy (2) k phi, untrapped pitch angles &
      & (3) k phi, trapped pitch angles, (4) k apar, energy, (5) k apar, untrapped &
      & angles. Relative errors should be < 0.1. ", &
      "absolute error measures have units T_r/(e rho_r)", errest)
    call create_and_write_variable(gnostics, gnostics%rtype, "vspace_vnewk", "2t", &
      "If the simulation is set to vary the collisionality in order to keep &
      & error in velocity integrals to acceptable levels, contains species 1 &
      & collisionality in (1) pitch angle and (2) energy  ", &
      "v_thr/a", (/vnmult(1)*spec(1)%vnewk, vnmult(2)*spec(1)%vnewk/))
    ! This next statement causes annoying printout because of an error in netcdf 4.1
    ! The error is fixed in 4.2
    ! See
    ! http://www.unidata.ucar.edu/software/netcdf/docs/known_problems.html#f90-debug-segfault
    if (gnostics%write_max_verr) &
      call create_and_write_variable(gnostics, gnostics%rtype, "vspace_err_maxindex", "52t", &
        "Gives the (1) theta index, (2) ky index and (3) kx index of the maximum &
        & error resulting from the &
        & velocity space integrals in the calculation of the following quantities &
        & in the given dimensions: (1) k phi, energy (2) k phi, untrapped pitch angles &
        & (3) k phi, trapped pitch angles, (4) k apar, energy, (5) k apar, untrapped &
        & angles. Relative errors should be < 0.1. ", &
        "1", erridx)
    !call create_and_write_variable(gnostics, gnostics%rtype, "vspace_lpcfrac_energy", "t", &
      !"Fraction of free energy contained in the high order coefficients of &
        !& the Legendre polynomial transform of energy space (should ideally be < 0.1). ", "1", geavg)
    !call create_and_write_variable(gnostics, gnostics%rtype, "vspace_lpcfrac_passing", "t", &
      !"Fraction of free energy contained in the high order coefficients of &
        !& the Legendre polynomial transform of untrapped pitch angles (should ideally be < 0.1). ", "1", glavg)
    !if (nlambda - ng2 > 1) &
      !call create_and_write_variable(gnostics, gnostics%rtype, "vspace_lpcfrac_trapped", "t", &
        !"Fraction of free energy contained in the high order coefficients of &
          !& the Legendre polynomial transform of trapped pitch angles (should ideally be < 0.1). ", "1", gtavg)
    !call create_and_write_variable(gnostics, gnostics%rtype, "vspace_errabs_energy", "t", &
      !"Estimate of the absolute error in phi resulting from the integral over energy space. ", "T_r/(e rho_r)", errest(1,1))
    !call create_and_write_variable(gnostics, gnostics%rtype, "vspace_err_energy", "t", &
      !"Estimate of the relative error in phi resulting from the integral over energy &
      !& space (should ideally be < 0.1). ", "1", errest(1,2))
    !call create_and_write_variable(gnostics, gnostics%rtype, "vspace_errabs_passing", "t", &
      !"Estimate of the absolute error in phi resulting from the integral over passing pitch &
      !& angles.", "1", errest(2,1))
    !call create_and_write_variable(gnostics, gnostics%rtype, "vspace_err_passing", "t", &
      !"Estimate of the relative error in phi resulting from the integral over passing pitch &
      !& angles (should ideally be < 0.1). ", "1", errest(2,2))
    !if (nlambda - ng2 > 1) then
      !call create_and_write_variable(gnostics, gnostics%rtype, "vspace_errabs_trapped", "t", &
        !"Estimate of the absolute error in phi resulting from the integral over passing pitch &
        !& angles.", "1", errest(2,1))
      !call create_and_write_variable(gnostics, gnostics%rtype, "vspace_err_trapped", "t", &
        !"Estimate of the relative error in phi resulting from the integral over passing pitch &
        !& angles (should ideally be < 0.1). ", "1", errest(2,2))
    !endif
    !if (proc0) then
       !! write error estimates to .nc file          
       !!          call nc_loop_vres (nout, errest_by_mode, lpcoef_by_mode)

       !! write error estimates for ion dist. fn. at outboard midplane with ik=it=1 to ascii files
       !if (write_ascii) then
          !if (nlambda - ng2 > 1) then
             !write(lpc_unit,"(4(1x,e13.6))") user_time, geavg, glavg, gtavg
          !else
             !write(lpc_unit,"(3(1x,e13.6))") user_time, geavg, glavg
          !end if
          !write(res_unit,"(8(1x,e13.6))") user_time, errest(1,2), errest(2,2), errest(3,2), &
               !errest(4,2), errest(5,2), vnmult(1)*spec(1)%vnewk, vnmult(2)*spec(1)%vnewk
          !if (write_max_verr) then
             !write(res_unit2,"(3(i8),(1x,e13.6),3(i8),(1x,e13.6),3(i8),(1x,e13.6),3(i8),(1x,e13.6),3(i8),(1x,e13.6))") &
                  !erridx(1,1), erridx(1,2), erridx(1,3), errest(1,1), &
                  !erridx(2,1), erridx(2,2), erridx(2,3), errest(2,1), &
                  !erridx(3,1), erridx(3,2), erridx(3,3), errest(3,1), &
                  !erridx(4,1), erridx(4,2), erridx(4,3), errest(4,1), &
                  !erridx(5,1), erridx(5,2), erridx(5,3), errest(5,1)
          !end if
       !end if
    !end if
    deallocate(errest,erridx)
  end subroutine write_velocity_space_checks
end module diagnostics_write_velocity_space_checks
