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
      call create_and_write_variable(gnostics, SDATIO_INT, "vspace_err_maxindex", "53t", &
        "Gives the (1) theta index, (2) ky index and (3) kx index of the maximum &
        & error resulting from the &
        & velocity space integrals in the calculation of the following quantities &
        & in the given dimensions: (1) k phi, energy (2) k phi, untrapped pitch angles &
        & (3) k phi, trapped pitch angles, (4) k apar, energy, (5) k apar, untrapped &
        & angles. Relative errors should be < 0.1. ", &
        "1", erridx)
    deallocate(errest,erridx)
  end subroutine write_velocity_space_checks
end module diagnostics_write_velocity_space_checks
