
!> DO NOT EDIT THIS FILE
!! This file has been automatically generated using generate_diagnostics_config.rb

!> A module for handling the configuration of the diagnostics
!! module via the namelist diagnostics_config.
module diagnostics_config

  use simpledataio, only: sdatio_file
  use diagnostics_ascii, only: diagnostics_ascii_type



  public :: init_diagnostics_config
  public :: finish_diagnostics_config
  public ::diagnostics_type
  public :: results_summary_type
  public :: override_screen_printout_options


  !> A type for storing the current results of the simulation
  type results_summary_type
    real :: phi2
    real :: apar2
    real :: bpar2
    real :: total_heat_flux
    real :: total_momentum_flux
    real :: total_particle_flux
    real :: max_growth_rate

    ! Individual heat fluxes
    real, dimension(:), allocatable :: species_es_heat_flux
    real, dimension(:), allocatable :: species_apar_heat_flux
    real, dimension(:), allocatable :: species_bpar_heat_flux

    ! Total fluxes
    real, dimension(:), allocatable :: species_heat_flux
    real, dimension(:), allocatable :: species_momentum_flux
    real, dimension(:), allocatable :: species_particle_flux
    real, dimension(:), allocatable :: species_energy_exchange

    ! Average total fluxes
    real, dimension(:), allocatable :: species_heat_flux_avg
    real, dimension(:), allocatable :: species_momentum_flux_avg
    real, dimension(:), allocatable :: species_particle_flux_avg

    ! Heating
    real, dimension(:), allocatable :: species_heating
    real, dimension(:), allocatable :: species_heating_avg
  end type results_summary_type

  !> A type for storing the diagnostics configuration,
  !! a reference to the output file, and current 
  !! results of the simulation
  type diagnostics_type
   type(sdatio_file) :: sfile
   !type(sdatio_file) :: sfilemovie
   type(diagnostics_ascii_type) :: ascii_files
   type(results_summary_type) :: current_results
   !> Integer below gives the sdatio type 
   !! which corresponds to a gs2 real
   integer :: rtype
   integer :: istep
   logical :: create
   logical :: wryte
   logical :: distributed
   logical :: parallel
   logical :: exit
   logical :: vary_vnew_only
   logical :: calculate_fluxes
   real :: user_time
   real :: user_time_old
   real, dimension(:), allocatable :: fluxfac
   integer :: nwrite
   integer :: nwrite_large
   logical :: write_any
   logical :: enable_parallel
   logical :: serial_netcdf4
   integer :: igomega
   logical :: print_line
   logical :: print_flux_line
   logical :: write_line
   logical :: write_flux_line
   logical :: write_fields
   logical :: write_phi_over_time
   logical :: write_apar_over_time
   logical :: write_bpar_over_time
   logical :: write_movie
   logical :: dump_fields_periodically
   logical :: write_moments
   logical :: write_full_moments_notgc
   logical :: write_ntot_over_time
   logical :: write_density_over_time
   logical :: write_upar_over_time
   logical :: write_tperp_over_time
   logical :: write_fluxes
   logical :: write_fluxes_by_mode
   logical :: write_symmetry
   logical :: write_parity
   logical :: write_omega
   integer :: navg
   real :: omegatinst
   real :: omegatol
   logical :: exit_when_converged
   logical :: write_verr
   logical :: write_max_verr
   integer :: ncheck
   logical :: write_heating
   logical :: write_ascii
   logical :: write_gyx
   logical :: write_g
   logical :: write_lpoly
   logical :: write_cerr
   integer :: conv_nstep_av
   real :: conv_test_multiplier
   integer :: conv_min_step
   integer :: conv_max_step
   integer :: conv_nsteps_converged
   logical :: use_nonlin_convergence
   logical :: write_cross_phase
   logical :: write_correlation
   logical :: write_correlation_extend
   logical :: write_jext
   logical :: write_lorentzian
   logical :: write_eigenfunc
   logical :: write_final_fields
   logical :: write_kpar
   logical :: write_final_epar
   logical :: write_final_db
   logical :: write_final_moments
   logical :: write_final_antot
   logical :: write_gs
   logical :: save_for_restart
   logical :: save_distfn
  end type diagnostics_type


  private

  !> Used for testing... causes screen printout to be 
  !! generated regardless of the values of print_line 
  !! and print_flux_line if set to true
  logical :: override_screen_printout_options = .false.

contains
  subroutine init_diagnostics_config(gnostics)
    use file_utils, only: open_output_file
    type(diagnostics_type), intent(out) :: gnostics
    call read_parameters(gnostics)
    call allocate_current_results(gnostics)
  end subroutine init_diagnostics_config

  subroutine finish_diagnostics_config(gnostics)
    type(diagnostics_type), intent(out) :: gnostics
    call deallocate_current_results(gnostics)
  end subroutine finish_diagnostics_config

  subroutine allocate_current_results(gnostics)
    use species, only: nspec
    type(diagnostics_type), intent(inout) :: gnostics
    allocate(gnostics%current_results%species_es_heat_flux(nspec))
    allocate(gnostics%current_results%species_apar_heat_flux(nspec))
    allocate(gnostics%current_results%species_bpar_heat_flux(nspec))

    allocate(gnostics%current_results%species_heat_flux(nspec))
    allocate(gnostics%current_results%species_momentum_flux(nspec))
    allocate(gnostics%current_results%species_particle_flux(nspec))
    allocate(gnostics%current_results%species_energy_exchange(nspec))
    allocate(gnostics%current_results%species_heat_flux_avg(nspec))
    allocate(gnostics%current_results%species_momentum_flux_avg(nspec))
    allocate(gnostics%current_results%species_particle_flux_avg(nspec))
    allocate(gnostics%current_results%species_heating(nspec))
    allocate(gnostics%current_results%species_heating_avg(nspec))
  end subroutine allocate_current_results

  subroutine deallocate_current_results(gnostics)
    type(diagnostics_type), intent(inout) :: gnostics
    ! This routine needs to be fixed: I don't know 
    ! how to correctly deallocate the derived type
    return
    
    ! One call deallocates gnostics and all allocatable arrays 
    ! within it
    !deallocate(gnostics%current_results)
    write (*,*) "DEALLOCATING"
    deallocate(gnostics%current_results%species_heat_flux)
    deallocate(gnostics%current_results%species_momentum_flux)
    deallocate(gnostics%current_results%species_particle_flux)
    deallocate(gnostics%current_results%species_heat_flux_avg)
    deallocate(gnostics%current_results%species_momentum_flux_avg)
    deallocate(gnostics%current_results%species_particle_flux_avg)
    deallocate(gnostics%current_results%species_heating)
    deallocate(gnostics%current_results%species_heating_avg)
  end subroutine deallocate_current_results


  subroutine read_parameters(gnostics)
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    implicit none
    type(diagnostics_type), intent(out) :: gnostics
    integer :: nwrite
    integer :: nwrite_large
    logical :: write_any
    logical :: enable_parallel
    logical :: serial_netcdf4
    integer :: igomega
    logical :: print_line
    logical :: print_flux_line
    logical :: write_line
    logical :: write_flux_line
    logical :: write_fields
    logical :: write_phi_over_time
    logical :: write_apar_over_time
    logical :: write_bpar_over_time
    logical :: write_movie
    logical :: dump_fields_periodically
    logical :: write_moments
    logical :: write_full_moments_notgc
    logical :: write_ntot_over_time
    logical :: write_density_over_time
    logical :: write_upar_over_time
    logical :: write_tperp_over_time
    logical :: write_fluxes
    logical :: write_fluxes_by_mode
    logical :: write_symmetry
    logical :: write_parity
    logical :: write_omega
    integer :: navg
    real :: omegatinst
    real :: omegatol
    logical :: exit_when_converged
    logical :: write_verr
    logical :: write_max_verr
    integer :: ncheck
    logical :: write_heating
    logical :: write_ascii
    logical :: write_gyx
    logical :: write_g
    logical :: write_lpoly
    logical :: write_cerr
    integer :: conv_nstep_av
    real :: conv_test_multiplier
    integer :: conv_min_step
    integer :: conv_max_step
    integer :: conv_nsteps_converged
    logical :: use_nonlin_convergence
    logical :: write_cross_phase
    logical :: write_correlation
    logical :: write_correlation_extend
    logical :: write_jext
    logical :: write_lorentzian
    logical :: write_eigenfunc
    logical :: write_final_fields
    logical :: write_kpar
    logical :: write_final_epar
    logical :: write_final_db
    logical :: write_final_moments
    logical :: write_final_antot
    logical :: write_gs
    logical :: save_for_restart
    logical :: save_distfn
    namelist /diagnostics_config/ &
      nwrite, &
      nwrite_large, &
      write_any, &
      enable_parallel, &
      serial_netcdf4, &
      igomega, &
      print_line, &
      print_flux_line, &
      write_line, &
      write_flux_line, &
      write_fields, &
      write_phi_over_time, &
      write_apar_over_time, &
      write_bpar_over_time, &
      write_movie, &
      dump_fields_periodically, &
      write_moments, &
      write_full_moments_notgc, &
      write_ntot_over_time, &
      write_density_over_time, &
      write_upar_over_time, &
      write_tperp_over_time, &
      write_fluxes, &
      write_fluxes_by_mode, &
      write_symmetry, &
      write_parity, &
      write_omega, &
      navg, &
      omegatinst, &
      omegatol, &
      exit_when_converged, &
      write_verr, &
      write_max_verr, &
      ncheck, &
      write_heating, &
      write_ascii, &
      write_gyx, &
      write_g, &
      write_lpoly, &
      write_cerr, &
      conv_nstep_av, &
      conv_test_multiplier, &
      conv_min_step, &
      conv_max_step, &
      conv_nsteps_converged, &
      use_nonlin_convergence, &
      write_cross_phase, &
      write_correlation, &
      write_correlation_extend, &
      write_jext, &
      write_lorentzian, &
      write_eigenfunc, &
      write_final_fields, &
      write_kpar, &
      write_final_epar, &
      write_final_db, &
      write_final_moments, &
      write_final_antot, &
      write_gs, &
      save_for_restart, &
      save_distfn

    integer :: in_file
    logical :: exist

    if (proc0) then
      nwrite = 10
      nwrite_large = 100
      write_any = .true.
      enable_parallel = .false.
      serial_netcdf4 = .true.
      igomega = 0
      print_line = .false.
      print_flux_line = .false.
      write_line = .true.
      write_flux_line = .true.
      write_fields = .true.
      write_phi_over_time = .false.
      write_apar_over_time = .false.
      write_bpar_over_time = .false.
      write_movie = .false.
      dump_fields_periodically = .false.
      write_moments = .true.
      write_full_moments_notgc = .false.
      write_ntot_over_time = .false.
      write_density_over_time = .false.
      write_upar_over_time = .false.
      write_tperp_over_time = .false.
      write_fluxes = .true.
      write_fluxes_by_mode = .false.
      write_symmetry = .false.
      write_parity = .false.
      write_omega = .true.
      navg = 10
      omegatinst = 1.0e6
      omegatol = -0.001
      exit_when_converged = .true.
      write_verr = .true.
      write_max_verr = .false.
      ncheck = 10
      write_heating = .false.
      write_ascii = .true.
      write_gyx = .false.
      write_g = .false.
      write_lpoly = .false.
      write_cerr = .false.
      conv_nstep_av = 4000
      conv_test_multiplier = 4e-1
      conv_min_step = 4000
      conv_max_step = 80000
      conv_nsteps_converged = 10000
      use_nonlin_convergence = .false.
      write_cross_phase = .false.
      write_correlation = .true.
      write_correlation_extend = .false.
      write_jext = .false.
      write_lorentzian = .false.
      write_eigenfunc = .false.
      write_final_fields = .false.
      write_kpar = .false.
      write_final_epar = .false.
      write_final_db = .false.
      write_final_moments = .false.
      write_final_antot = .false.
      write_gs = .false.
      save_for_restart = .true.
      save_distfn = .false.

      in_file = input_unit_exist ("diagnostics_config", exist)
      if (exist) read (unit=in_file, nml=diagnostics_config)

      gnostics%nwrite = nwrite
      gnostics%nwrite_large = nwrite_large
      gnostics%write_any = write_any
      gnostics%enable_parallel = enable_parallel
      gnostics%serial_netcdf4 = serial_netcdf4
      gnostics%igomega = igomega
      gnostics%print_line = print_line
      gnostics%print_flux_line = print_flux_line
      gnostics%write_line = write_line
      gnostics%write_flux_line = write_flux_line
      gnostics%write_fields = write_fields
      gnostics%write_phi_over_time = write_phi_over_time
      gnostics%write_apar_over_time = write_apar_over_time
      gnostics%write_bpar_over_time = write_bpar_over_time
      gnostics%write_movie = write_movie
      gnostics%dump_fields_periodically = dump_fields_periodically
      gnostics%write_moments = write_moments
      gnostics%write_full_moments_notgc = write_full_moments_notgc
      gnostics%write_ntot_over_time = write_ntot_over_time
      gnostics%write_density_over_time = write_density_over_time
      gnostics%write_upar_over_time = write_upar_over_time
      gnostics%write_tperp_over_time = write_tperp_over_time
      gnostics%write_fluxes = write_fluxes
      gnostics%write_fluxes_by_mode = write_fluxes_by_mode
      gnostics%write_symmetry = write_symmetry
      gnostics%write_parity = write_parity
      gnostics%write_omega = write_omega
      gnostics%navg = navg
      gnostics%omegatinst = omegatinst
      gnostics%omegatol = omegatol
      gnostics%exit_when_converged = exit_when_converged
      gnostics%write_verr = write_verr
      gnostics%write_max_verr = write_max_verr
      gnostics%ncheck = ncheck
      gnostics%write_heating = write_heating
      gnostics%write_ascii = write_ascii
      gnostics%write_gyx = write_gyx
      gnostics%write_g = write_g
      gnostics%write_lpoly = write_lpoly
      gnostics%write_cerr = write_cerr
      gnostics%conv_nstep_av = conv_nstep_av
      gnostics%conv_test_multiplier = conv_test_multiplier
      gnostics%conv_min_step = conv_min_step
      gnostics%conv_max_step = conv_max_step
      gnostics%conv_nsteps_converged = conv_nsteps_converged
      gnostics%use_nonlin_convergence = use_nonlin_convergence
      gnostics%write_cross_phase = write_cross_phase
      gnostics%write_correlation = write_correlation
      gnostics%write_correlation_extend = write_correlation_extend
      gnostics%write_jext = write_jext
      gnostics%write_lorentzian = write_lorentzian
      gnostics%write_eigenfunc = write_eigenfunc
      gnostics%write_final_fields = write_final_fields
      gnostics%write_kpar = write_kpar
      gnostics%write_final_epar = write_final_epar
      gnostics%write_final_db = write_final_db
      gnostics%write_final_moments = write_final_moments
      gnostics%write_final_antot = write_final_antot
      gnostics%write_gs = write_gs
      gnostics%save_for_restart = save_for_restart
      gnostics%save_distfn = save_distfn

    end if

    call broadcast (gnostics%nwrite)
    call broadcast (gnostics%nwrite_large)
    call broadcast (gnostics%write_any)
    call broadcast (gnostics%enable_parallel)
    call broadcast (gnostics%serial_netcdf4)
    call broadcast (gnostics%igomega)
    call broadcast (gnostics%print_line)
    call broadcast (gnostics%print_flux_line)
    call broadcast (gnostics%write_line)
    call broadcast (gnostics%write_flux_line)
    call broadcast (gnostics%write_fields)
    call broadcast (gnostics%write_phi_over_time)
    call broadcast (gnostics%write_apar_over_time)
    call broadcast (gnostics%write_bpar_over_time)
    call broadcast (gnostics%write_movie)
    call broadcast (gnostics%dump_fields_periodically)
    call broadcast (gnostics%write_moments)
    call broadcast (gnostics%write_full_moments_notgc)
    call broadcast (gnostics%write_ntot_over_time)
    call broadcast (gnostics%write_density_over_time)
    call broadcast (gnostics%write_upar_over_time)
    call broadcast (gnostics%write_tperp_over_time)
    call broadcast (gnostics%write_fluxes)
    call broadcast (gnostics%write_fluxes_by_mode)
    call broadcast (gnostics%write_symmetry)
    call broadcast (gnostics%write_parity)
    call broadcast (gnostics%write_omega)
    call broadcast (gnostics%navg)
    call broadcast (gnostics%omegatinst)
    call broadcast (gnostics%omegatol)
    call broadcast (gnostics%exit_when_converged)
    call broadcast (gnostics%write_verr)
    call broadcast (gnostics%write_max_verr)
    call broadcast (gnostics%ncheck)
    call broadcast (gnostics%write_heating)
    call broadcast (gnostics%write_ascii)
    call broadcast (gnostics%write_gyx)
    call broadcast (gnostics%write_g)
    call broadcast (gnostics%write_lpoly)
    call broadcast (gnostics%write_cerr)
    call broadcast (gnostics%conv_nstep_av)
    call broadcast (gnostics%conv_test_multiplier)
    call broadcast (gnostics%conv_min_step)
    call broadcast (gnostics%conv_max_step)
    call broadcast (gnostics%conv_nsteps_converged)
    call broadcast (gnostics%use_nonlin_convergence)
    call broadcast (gnostics%write_cross_phase)
    call broadcast (gnostics%write_correlation)
    call broadcast (gnostics%write_correlation_extend)
    call broadcast (gnostics%write_jext)
    call broadcast (gnostics%write_lorentzian)
    call broadcast (gnostics%write_eigenfunc)
    call broadcast (gnostics%write_final_fields)
    call broadcast (gnostics%write_kpar)
    call broadcast (gnostics%write_final_epar)
    call broadcast (gnostics%write_final_db)
    call broadcast (gnostics%write_final_moments)
    call broadcast (gnostics%write_final_antot)
    call broadcast (gnostics%write_gs)
    call broadcast (gnostics%save_for_restart)
    call broadcast (gnostics%save_distfn)
    
    if (override_screen_printout_options) then 
      gnostics%print_line = .true.
      gnostics%print_flux_line = .true.
    end if



  end subroutine read_parameters
end module diagnostics_config



