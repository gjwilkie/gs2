! what about boundary condition contributions?  In the presence 
! of magnetic shear, there are not always zero amplitudes at the ends
! of the supercells.

module gs2_diagnostics
  use gs2_heating, only: heating_diagnostics

  implicit none

  private

  public :: read_parameters
  public :: init_gs2_diagnostics
  public :: finish_gs2_diagnostics
  public :: check_gs2_diagnostics
  public :: wnml_gs2_diagnostics
  public :: loop_diagnostics
  public :: ensemble_average
  public :: reset_init
  public :: pflux_avg, qflux_avg, heat_avg, vflux_avg, start_time

  public :: get_omegaavg

  interface get_vol_average
     module procedure get_vol_average_one, get_vol_average_all
  end interface

  interface get_vol_int
     module procedure get_vol_int_one, get_vol_int_all
  end interface

  interface get_fldline_avg
     module procedure get_fldline_avg_r, get_fldline_avg_c
  end interface

!CMR, 17/11/2009:   read_parameters now public so ingen can USE instead of copy
!
! Why are these variables public?  This is not good.
  real, public :: omegatol, omegatinst
  logical, public :: print_line, print_flux_line
  logical, public :: print_summary, write_line, write_flux_line
  logical, public :: write_omega, write_omavg, write_ascii
  logical, public :: write_gs, write_lpoly
  logical, public :: write_g, write_gg, write_gyx
  logical, public :: write_eigenfunc, write_fields, write_final_fields, write_final_antot
  logical, public :: write_final_moments, write_avg_moments, write_parity
  logical, public :: write_moments, write_final_db
  logical, public :: write_full_moments_notgc, write_cross_phase = .false.
  logical, public :: write_final_epar, write_kpar
  logical, public :: write_hrate, write_lorentzian, write_hrate_e
  logical, public :: write_nl_flux
  logical, public :: write_verr, write_cerr, write_max_verr
  logical, public :: exit_when_converged
  logical, public :: dump_check1, dump_check2
  logical, public :: dump_fields_periodically, make_movie
  logical, public :: save_for_restart
  logical, public :: save_distfn !<DD> Added for saving distribution function
  logical, public :: write_symmetry, write_correlation_extend, write_correlation
  integer, public :: nwrite, igomega, nmovie
  integer, public :: navg, nsave, nwrite_mult
  logical, public :: write_flux_emu, write_flux_e, write_eflux

  logical, public :: write_phi_over_time, write_apar_over_time, write_bpar_over_time !EGH
!>GGH
  logical :: write_jext=.false.
!<GGH

  ! internal
  logical :: write_any, write_any_fluxes, dump_any
  logical, private :: initialized = .false.

  namelist /gs2_diagnostics_knobs/ print_line, print_flux_line, &
         write_line, write_flux_line, &
         write_omega, write_omavg, write_ascii, write_kpar, &
         write_gs, write_gyx, write_g, write_gg, write_hrate, write_hrate_e, write_lpoly, &
         write_eigenfunc, write_fields, write_final_fields, write_final_antot, &
         write_final_epar, write_moments, write_final_moments, write_cerr, &
         write_verr, write_max_verr, write_nl_flux, write_final_db, &
         nwrite, nmovie, nsave, navg, omegatol, omegatinst, igomega, write_lorentzian, &
         exit_when_converged, write_avg_moments, &
         write_full_moments_notgc, write_cross_phase, &
         dump_check1, dump_check2, &
         dump_fields_periodically, make_movie, &
         save_for_restart, write_parity, write_symmetry, save_distfn, & !<DD> Added for saving distribution function
         write_correlation_extend, nwrite_mult, write_correlation, &
         write_phi_over_time, write_apar_over_time, write_bpar_over_time, write_flux_emu, write_flux_e, write_eflux

  integer :: out_unit, kp_unit, heat_unit, polar_raw_unit, polar_avg_unit, heat_unit2, lpc_unit
  integer :: jext_unit   !GGH Additions
  integer :: phase_unit
  integer :: dump_check1_unit, dump_check2_unit
  integer :: res_unit, res_unit2, parity_unit

  complex, dimension (:,:,:), allocatable :: omegahist
  ! (navg,ntheta0,naky)
  type (heating_diagnostics), save :: h
  type (heating_diagnostics), dimension(:), save, allocatable :: h_hist
  type (heating_diagnostics), dimension(:,:), save, allocatable :: hk
  type (heating_diagnostics), dimension(:), save, allocatable :: he
  type (heating_diagnostics), dimension(:,:,:), save, allocatable :: hk_hist
  type (heating_diagnostics), dimension(:,:), save, allocatable :: he_hist
  !GGH J_external
  real, dimension(:,:,:), allocatable ::  j_ext_hist

  real, dimension (:,:,:,:), allocatable ::  qheat, qmheat, qbheat
  ! (ntheta0,naky,nspec,3)

  real, dimension (:,:,:), allocatable ::  pflux,  vflux, vflux_par, vflux_perp
  real, dimension (:,:,:), allocatable :: vflux0, vflux1  ! low flow correction to turbulent momentum flux
  real, dimension (:,:,:), allocatable :: pmflux, vmflux
  real, dimension (:,:,:), allocatable :: pbflux, vbflux
  real, dimension (:,:,:), allocatable :: exchange
  real,dimension(:,:,:,:), allocatable :: pflux_emu, mflux_emu, bflux_emu
  real,dimension(:,:), allocatable :: pflux_e, mflux_e, bflux_e, trapped_flux, passing_flux
  real,dimension(:,:), allocatable :: peflux, meflux, beflux

  ! (ntheta0,naky,nspec)

  real :: start_time = 0.0
  real, dimension (:), allocatable :: pflux_avg, qflux_avg, heat_avg, vflux_avg

  integer :: ntg_out, ntg_extend, nth0_extend
  integer :: nout = 1
  integer :: nout_movie = 1
  integer :: nout_big = 1
  complex :: wtmp_old = 0.
  logical :: exist

contains
  !> Define NetCDF vars, call real_init, which calls read_parameters; broadcast all the different write flags. 
   subroutine wnml_gs2_diagnostics(unit)
   implicit none
   integer :: unit
       if (.not.exist) return
       write (unit, *)
       write (unit, fmt="(' &',a)") "gs2_diagnostics_knobs"
       write (unit, fmt="(' save_for_restart = ',L1)") save_for_restart
       write (unit, fmt="(' print_line = ',L1)") print_line 
       write (unit, fmt="(' write_line = ',L1)") write_line
       write (unit, fmt="(' print_flux_line = ',L1)") print_flux_line
       write (unit, fmt="(' write_flux_line = ',L1)") write_flux_line
       write (unit, fmt="(' nmovie = ',i6)") nmovie
       write (unit, fmt="(' nwrite_mult = ',i6)") nwrite_mult
       write (unit, fmt="(' nwrite = ',i6)") nwrite
       write (unit, fmt="(' nsave = ',i6)") nsave
       write (unit, fmt="(' navg = ',i6)") navg
       write (unit, fmt="(' omegatol = ',e16.10)") omegatol
       write (unit, fmt="(' omegatinst = ',e16.10)") omegatinst
! should be legal -- not checked yet
       if (igomega /= 0) write (unit, fmt="(' igomega = ',i6)") igomega  
       
       if (write_ascii) then
          write (unit, fmt="(' write_ascii = ',L1)") write_ascii
          write (unit, fmt="(' write_omega = ',L1)") write_omega
          write (unit, fmt="(' write_omavg = ',L1)") write_omavg
       end if
       write (unit, fmt="(' write_hrate = ',L1)") write_hrate
       write (unit, fmt="(' write_hrate_e = ',L1)") write_hrate_e
       write (unit, fmt="(' write_lorentzian = ',L1)") write_lorentzian
       write (unit, fmt="(' write_eigenfunc = ',L1)") write_eigenfunc
       write (unit, fmt="(' write_final_fields = ',L1)") write_final_fields
       write (unit, fmt="(' write_final_epar = ',L1)") write_final_epar
       write (unit, fmt="(' write_final_db = ',L1)") write_final_db
       write (unit, fmt="(' write_final_moments = ',L1)") write_final_moments
       write (unit, fmt="(' write_final_antot = ',L1)") write_final_antot
       write (unit, fmt="(' write_nl_flux = ',L1)") write_nl_flux
       write (unit, fmt="(' write_flux_emu = ',L1)") write_flux_emu
       write (unit, fmt="(' write_flux_e = ',L1)") write_flux_e
       write (unit, fmt="(' write_eflux = ',L1)") write_eflux
       write (unit, fmt="(' exit_when_converged = ',L1)") exit_when_converged
       if (write_avg_moments) write (unit, fmt="(' write_avg_moments = ',L1)") write_avg_moments
       if (dump_check1) write (unit, fmt="(' dump_check1 = ',L1)") dump_check1
       if (dump_check2) write (unit, fmt="(' dump_check2 = ',L1)") dump_check2
       if (dump_fields_periodically) &
            write (unit, fmt="(' dump_fields_periodically = ',L1)") dump_fields_periodically
       if (make_movie) &
            write (unit, fmt="(' make_movie = ',L1)") make_movie

       write (unit, fmt="(' /')")       
   end subroutine wnml_gs2_diagnostics

   subroutine check_gs2_diagnostics(report_unit)
   use file_utils, only: run_name
   use nonlinear_terms, only: nonlinear_mode_switch, nonlinear_mode_on
   use dist_fn, only : def_parity, even 
   use kt_grids, only : gridopt_switch, gridopt_box
   use init_g, only : restart_file
   implicit none
   integer :: report_unit
    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 
    write (report_unit, fmt="('Diagnostic control section.')")

    if (print_line) then
       write (report_unit, fmt="('print_line = T:            Estimated frequencies &
          & output to the screen every ',i4,' steps.')") nwrite
    else
       ! nothing
    end if

    if (write_line) then
       if (write_ascii) then
          write (report_unit, fmt="('write_line = T:            Estimated frequencies output to ',a,' every ',i4,' steps.')") &
               & trim(run_name)//'.out',  nwrite
       end if
       write (report_unit, fmt="('write_line = T:            Estimated frequencies output to ',a,' every ',i4,' steps.')") &
            & trim(run_name)//'.out.nc',  nwrite
    else
       ! nothing
    end if

    if (print_flux_line) then
       write (report_unit, fmt="('print_flux_line = T:       Instantaneous fluxes output to screen every ', &
             & i4,' steps.')") nwrite
    else
       ! nothing
    end if

    if (write_flux_line) then
       if (write_ascii) then
          write (report_unit, fmt="('write_flux_line = T:       Instantaneous fluxes output to ',a,' every ',i4,' steps.')") &
               & trim(run_name)//'.out',  nwrite
       end if
       write (report_unit, fmt="('write_flux_line = T:       Instantaneous fluxes output to ',a,' every ',i4,' steps.')") &
            & trim(run_name)//'.out.nc',  nwrite
    else
       ! nothing
    end if

    if (write_omega) then
       if (write_ascii) then
          write (report_unit, fmt="('write_omega = T:           Instantaneous frequencies written to ',a)") trim(run_name)//'.out'
       else
          write (report_unit, fmt="('write_omega = T:           No effect.')")
       end if
       write (report_unit, fmt="('                           Frequencies calculated at igomega = ',i4)") igomega
       if (def_parity .and. .not. even) then
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('   You probably want igomega /= 0 for odd parity modes.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if
    end if

    if (write_omavg) then
       if (write_ascii) then
          write (report_unit, fmt="('write_omavg = T:           Time-averaged frequencies written to ',a)") trim(run_name)//'.out'
          write (report_unit, fmt="('                           Averages taken over ',i4,' timesteps.')") navg
       else
          write (report_unit, fmt="('write_omavg = T:           No effect.')")
       end if
    end if

    if (write_ascii) then
       write (report_unit, fmt="('write_ascii = T:           Write some data to ',a)") trim(run_name)//'.out'
    end if

    if (write_eigenfunc) then
       if (write_ascii) then
          write (report_unit, fmt="('write_eigenfunc = T:       Normalized Phi(theta) written to ',a)") trim(run_name)//'.eigenfunc'
       end if
       write (report_unit, fmt="('write_eigenfunc = T:       Normalized Phi(theta) written to ',a)") trim(run_name)//'.out.nc'
    end if

    if (write_final_fields) then
       if (write_ascii) then
          write (report_unit, fmt="('write_final_fields = T:    Phi(theta), etc. written to ',a)") trim(run_name)//'.fields'
       end if
       write (report_unit, fmt="('write_final_fields = T:    Phi(theta), etc. written to ',a)") trim(run_name)//'.out.nc'
    end if

    if (write_final_antot) then
       if (write_ascii) then
          write (report_unit, fmt="('write_final_antot = T:          Sources for Maxwell eqns. written to ',a)") &
          	& trim(run_name)//'.antot'
       end if
       write (report_unit, fmt="('write_final_antot = T:          Sources for Maxwell eqns. written to ',a)") &
	& trim(run_name)//'.out.nc'
    end if

    if (write_final_moments) then
       if (write_ascii) then
          write (report_unit, fmt="('write_final_moments = T:   Low-order moments of g written to ',a)") &
               & trim(run_name)//'.moments'
          write (report_unit, fmt="('write_final_moments = T:   int dl/B average of low-order moments of g written to ',a)") &
               & trim(run_name)//'.amoments'
       end if
       write (report_unit, fmt="('write_final_moments = T:   Low-order moments of g written to ',a)") &
            & trim(run_name)//'.out.nc'
       write (report_unit, fmt="('write_final_moments = T:   int dl/B average of low-order moments of g written to ',a)") &
            & trim(run_name)//'.out.nc'
    end if

    if (write_avg_moments) then
       if (gridopt_switch /= gridopt_box) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('write_avg_moments = T:          Ignored unless grid_option=box')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
       else
          if (write_ascii) then
             write (report_unit, fmt="('write_avg_moments = T:     Flux surface averaged low-order moments of g written to ',a)") &
                  & trim(run_name)//'.moments'
          end if
          write (report_unit, fmt="('write_avg_moments = T:     Flux surface averaged low-order moments of g written to ',a)") &
               & trim(run_name)//'.out.nc'
       end if
    end if

    if (write_final_epar) then
       if (write_ascii) then
          write (report_unit, fmt="('write_final_epar = T:      E_parallel(theta) written to ',a)") trim(run_name)//'.epar'
       end if
       write (report_unit, fmt="('write_final_epar = T:      E_parallel(theta) written to ',a)") trim(run_name)//'.out.nc'
    end if

    if (write_nl_flux) then
       if (write_ascii) then
          write (report_unit, fmt="('write_nl_flux = T:         Phi**2(kx, ky) written to ',a)") trim(run_name)//'.out'
       end if
    else
       write (report_unit, fmt="('write_nl_flux = F:         Phi**2(kx, ky) NOT written to ',a)") trim(run_name)//'.out'
    end if

    if (dump_check1) then
       write (report_unit, fmt="('dump_check1 = T:          Field-line avg of Phi written to ',a)") 'dump.check1'
       write (report_unit, fmt="('This option is usually used for Rosenbluth-Hinton calculations.')") 
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
    end if

    if (dump_check2) then
       write (report_unit, fmt="('dump_check2 = T:           Apar(kx, ky, igomega) written to ',a)") trim(run_name)//'.dc2'
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
    end if

    if (dump_fields_periodically) then
       write (report_unit, fmt="('dump_fields_periodically = T:          Phi, Apar, Bpar written to ',a)") 'dump.fields.t=(time)'
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.  IT IS EXPENSIVE.')") 
    end if

    if (save_for_restart) then
       write (report_unit, fmt="('save_for_restart = T:      Restart files written to ',a)") trim(restart_file)//'.(PE)'
    else
       if (nonlinear_mode_switch == nonlinear_mode_on) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('save_for_restart = F:              This run cannot be continued.')")
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if
    end if

  end subroutine check_gs2_diagnostics


  subroutine init_gs2_diagnostics (list, nstep)
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids, ntheta0, naky
    use run_parameters, only: init_run_parameters
    use species, only: init_species, nspec
    use dist_fn, only: init_dist_fn
    use init_g, only: init_init_g
    use gs2_io, only: init_gs2_io
    use gs2_heating, only: init_htype, init_hetype
    use collisions, only: collision_model_switch, init_lorentz_error
    use mp, only: broadcast, proc0
    use le_grids, only: init_weights, init_le_grids, negrid

    implicit none
    logical, intent (in) :: list
    integer, intent (in) :: nstep
    integer :: nmovie_tot, nwrite_big_tot
    logical :: dum1,dum2

    if (initialized) return
    initialized = .true.

    call init_theta_grid
    call init_kt_grids
    call init_run_parameters
    call init_species
    call init_le_grids(dum1,dum2)
    call init_init_g
    call init_dist_fn

    call real_init (list)
    call broadcast (navg)
    call broadcast (nwrite)
    call broadcast (nmovie)
    call broadcast (nwrite_mult)
    call broadcast (nsave)
    call broadcast (write_any)
    call broadcast (write_any_fluxes)
    call broadcast (write_cross_phase)
    call broadcast (write_nl_flux)
    call broadcast (write_omega)
    call broadcast (dump_any)
    call broadcast (dump_check1)
    call broadcast (dump_check2)
    call broadcast (write_fields)
    call broadcast (write_moments)
    call broadcast (dump_fields_periodically)
    call broadcast (make_movie)
    call broadcast (save_for_restart)
    call broadcast (save_distfn) !<DD> Added for saving distribution function
    call broadcast (write_gs)
    call broadcast (write_g)
    call broadcast (write_gyx)
    call broadcast (write_gg)
    call broadcast (write_final_antot)
    call broadcast (write_verr)
    call broadcast (write_max_verr)
    call broadcast (write_lpoly)
    call broadcast (write_cerr)

    call broadcast (ntg_out)
    call broadcast (write_hrate)
    call broadcast (write_hrate_e)
    call broadcast (write_lorentzian)
    call broadcast (write_eigenfunc)
    call broadcast (write_flux_emu)
    call broadcast (write_flux_e)
    call broadcast (write_eflux)

    call broadcast (write_full_moments_notgc)
    call broadcast (write_phi_over_time)
    call broadcast (write_apar_over_time)
    call broadcast (write_bpar_over_time)

    nmovie_tot = nstep/nmovie
    nwrite_big_tot = nstep/(nwrite*nwrite_mult)-nstep/4/(nwrite*nwrite_mult)

! initialize weights for less accurate integrals used
! to provide an error estimate for v-space integrals (energy and untrapped)
    if (write_verr .and. proc0) call init_weights

! allocate heating diagnostic data structures
    if (write_hrate) then
       allocate (h_hist(0:navg-1))
       call init_htype (h_hist,  nspec)

       allocate (hk_hist(ntheta0,naky,0:navg-1))
       call init_htype (hk_hist, nspec)

       call init_htype (h,  nspec)

       allocate (hk(ntheta0, naky))
       call init_htype (hk, nspec)

    else
       allocate (h_hist(0))
       allocate (hk(1,1))
       allocate (hk_hist(1,1,0))
    end if

    if (write_hrate_e) then
       allocate (he_hist(negrid,0:navg-1))
       call init_hetype (he_hist,  nspec)

       allocate (he(negrid))
       call init_hetype (he, nspec)
    else
       allocate (he_hist(1,0))
       allocate (he(1))
    end if
       
!GGH Allocate density and velocity perturbation diagnostic structures
    if (write_jext) allocate (j_ext_hist(ntheta0, naky,0:navg-1)) 

    call init_gs2_io (write_nl_flux, write_omega, &
         write_hrate, write_hrate_e, write_final_antot, &
         write_eigenfunc, make_movie, nmovie_tot, write_verr, &
         write_fields, write_moments, write_full_moments_notgc, &
         write_symmetry, &
         write_correlation, nwrite_big_tot, write_correlation_extend, &
         write_phi_over_time, write_apar_over_time, write_bpar_over_time,write_flux_emu,write_flux_e, write_eflux)
    
    if (write_cerr) then
       if (collision_model_switch == 1 .or. collision_model_switch == 5) then
          call init_lorentz_error
       else
          write_cerr = .false.
       end if
    end if

    allocate (pflux_avg(nspec), qflux_avg(nspec), heat_avg(nspec), vflux_avg(nspec))
    pflux_avg = 0.0 ; qflux_avg = 0.0 ; heat_avg = 0.0 ; vflux_avg = 0.0

  end subroutine init_gs2_diagnostics
 

  subroutine real_init (list)
    use le_grids, only: negrid, nlambda
    use run_parameters, only: fapar
    use file_utils, only: open_output_file, get_unused_unit
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, aky, akx
    use gs2_layouts, only: yxf_lo
    use species, only: nspec
    use mp, only: proc0, iproc, broadcast
    use constants
    implicit none
    logical, intent (in) :: list
    character(20) :: datestamp, timestamp, zone

    !<doc> Call read_parameters </doc>
    call read_parameters (list)
    !<doc> Open the various ascii output files (depending on the write flags) </doc>
    if (proc0) then
       if (write_ascii) then
          call open_output_file (out_unit, ".out")          
!          if (write_kpar) call open_output_file (kp_unit, ".kp")
       end if             

       if (write_cross_phase .and. write_ascii) then
          call open_output_file (phase_unit, ".phase")
       end if

       if (write_hrate .and. write_ascii) then
          call open_output_file (heat_unit, ".heat")
          call open_output_file (heat_unit2, ".heat2")
       end if

       if (write_verr .and. write_ascii) then
          call open_output_file (res_unit, ".vres")
          call open_output_file (lpc_unit, ".lpc")
          if (write_max_verr) call open_output_file (res_unit2, ".vres2")
       end if

       if (write_parity .and. write_ascii) then
          call open_output_file (parity_unit, ".parity")
       end if

       !GGH J_external, only if A_parallel is being calculated.
       if (write_jext .and. fapar > epsilon(0.0)) then
          if (write_ascii) then
             call open_output_file (jext_unit, ".jext")
          end if
       else
          write_jext = .false.
       end if

       if (dump_check1) then
          call get_unused_unit (dump_check1_unit)
          open (unit=dump_check1_unit, file="dump.check1", status="unknown")
       end if
       
       if (dump_check2) then
          call get_unused_unit (dump_check2_unit)
          call open_output_file (dump_check2_unit, ".dc2")
       end if
       
       if (write_ascii) then
          write (unit=out_unit, fmt="('gs2')")
          datestamp(:) = ' '
          timestamp(:) = ' '
          zone(:) = ' '
          call date_and_time (datestamp, timestamp, zone)
          write (unit=out_unit, fmt="('Date: ',a,' Time: ',a,1x,a)") &
               trim(datestamp), trim(timestamp), trim(zone)
       end if
       
       allocate (omegahist(0:navg-1,ntheta0,naky))
       omegahist = 0.0
    end if

    call broadcast(write_flux_emu)
    call broadcast(write_flux_e)
    call broadcast(write_eflux)

    !<doc> Allocate arrays for storing the various fluxes which the diagnostics will output </doc>
    allocate (pflux (ntheta0,naky,nspec)) ; pflux = 0.
    allocate (qheat (ntheta0,naky,nspec,3)) ; qheat = 0.
    allocate (vflux (ntheta0,naky,nspec)) ; vflux = 0.
    allocate (exchange (ntheta0,naky,nspec)) ; exchange = 0.

    allocate (vflux_par (ntheta0,naky,nspec)) ; vflux_par = 0.
    allocate (vflux_perp (ntheta0,naky,nspec)) ; vflux_perp = 0.

    allocate (vflux0 (ntheta0,naky,nspec)) ; vflux0 = 0.
    allocate (vflux1 (ntheta0,naky,nspec)) ; vflux1 = 0.

    allocate (pmflux(ntheta0,naky,nspec)) ; pmflux = 0.
    allocate (qmheat(ntheta0,naky,nspec,3)) ; qmheat = 0.
    allocate (vmflux(ntheta0,naky,nspec)) ; vmflux = 0.

    allocate (pbflux(ntheta0,naky,nspec)) ; pbflux = 0.
    allocate (qbheat(ntheta0,naky,nspec,3)) ; qbheat = 0.
    allocate (vbflux(ntheta0,naky,nspec)) ; vbflux = 0.

    if (write_flux_emu) then
       allocate (pflux_emu(negrid,nlambda,2,nspec)) ; pflux_emu = 0.
       allocate (mflux_emu(negrid,nlambda,2,nspec)) ; mflux_emu = 0.
       allocate (bflux_emu(negrid,nlambda,2,nspec)) ; bflux_emu = 0.
    end if
      
    if (write_flux_e) then
       allocate (pflux_e(negrid,nspec)) ; pflux_e = 0.
       allocate (mflux_e(negrid,nspec)) ; mflux_e = 0.
       allocate (bflux_e(negrid,nspec)) ; bflux_e = 0.
       allocate (trapped_flux(negrid,nspec)) ; trapped_flux = 0.
       allocate (passing_flux(negrid,nspec)) ; trapped_flux = 0.
    end if
      
    if (write_eflux) then
       allocate (peflux(negrid,nspec)) ; peflux = 0.
       allocate (meflux(negrid,nspec)) ; meflux = 0.
       allocate (beflux(negrid,nspec)) ; beflux = 0.
    end if

  end subroutine real_init

  subroutine read_parameters (list)
!CMR, 17/11/2009:   namelist gs2_diagnostics_knobs made public
!                   so that ingen can just USE it instead of copying it!
!
    use file_utils, only: input_unit, input_unit_exist
    use theta_grid, only: nperiod, ntheta
    use kt_grids, only: box, nx, ny
    use mp, only: proc0
    implicit none
    integer :: in_file
    logical, intent (in) :: list
    !<doc> Set defaults for the gs2_diagnostics_knobs</doc>		
    if (proc0) then
	!<wkdoc> Set defaults for the gs2_diagnostics_knobs</wkdoc>		
       print_line = .true.
       print_flux_line = .false.
       write_line = .true.
       write_flux_line = .true.
       write_kpar = .false.
       write_hrate = .false.
       write_hrate_e = .false.
       write_gs = .false.
       write_g = .false.
       write_gyx = .false.
       write_lpoly = .false.
       write_gg = .false.
       write_lorentzian = .false.
       write_omega = .false.
       write_ascii = .true.
       write_omavg = .false.
       write_nl_flux = .false.
       write_eigenfunc = .false.
       write_moments = .false.
       write_final_moments = .false.
       write_avg_moments = .false.
       write_parity = .false.
       write_symmetry = .false.
       write_correlation_extend = .false.
       write_correlation = .false.
       write_fields = .false.
       write_full_moments_notgc = .false.
       write_final_fields = .false.
       write_final_antot = .false.
       write_final_epar = .false.
       write_final_db = .false.
       write_verr = .false.
       write_max_verr = .false.
       write_cerr = .false.
       nwrite = 100
       nmovie = 1000
       nwrite_mult = 10
       navg = 100
       nsave = -1
       omegatol = 1e-3
       omegatinst = 1.0
       igomega = 0
       exit_when_converged = .true.
       dump_check1 = .false.
       dump_check2 = .false.
       dump_fields_periodically = .false.
       make_movie = .false.
       save_for_restart = .false.
       save_distfn = .false. !<DD> Added for saving distribution function
       write_phi_over_time = .false.
       write_bpar_over_time = .false.
       write_apar_over_time = .false.
       in_file = input_unit_exist ("gs2_diagnostics_knobs", exist)
       write_flux_emu = .false.
       write_flux_e = .false.
       write_eflux = .false.

	!<doc> Read in parameters from the namelist gs2_diagnostics_knobs, if the namelist exists </doc>
!       if (exist) read (unit=input_unit("gs2_diagnostics_knobs"), nml=gs2_diagnostics_knobs)
       if (exist) read (unit=in_file, nml=gs2_diagnostics_knobs)

       print_summary = (list .and. (print_line .or. print_flux_line)) 

       if (list) then
          print_line = .false.
          print_flux_line = .false.
       end if

       if (.not. save_for_restart) nsave = -1
! changed temporarily for testing -- MAB
!       write_avg_moments = write_avg_moments .and. box
       write_avg_moments = write_avg_moments

       write_any = write_line .or. write_omega     .or. write_omavg &
            .or. write_flux_line                   .or. write_nl_flux  &
            .or. write_kpar   .or. write_hrate     .or. write_lorentzian &
            .or. write_gs .or. write_flux_emu .or. write_flux_e .or. write_eflux
       write_any_fluxes =  write_flux_line .or. print_flux_line  &
            .or. write_nl_flux  .or. write_flux_emu .or. write_flux_e .or. write_eflux
       dump_any = dump_check1  .or. dump_fields_periodically &
            .or.  dump_check2 .or. make_movie .or. print_summary &
            .or.  write_full_moments_notgc .or. write_hrate_e

       ntg_out = ntheta/2 + (nperiod-1)*ntheta
    end if 
 end subroutine read_parameters

  subroutine finish_gs2_diagnostics (istep)
    use file_utils, only: open_output_file, close_output_file, get_unused_unit
    use mp, only: proc0, broadcast, nproc, iproc, sum_reduce
    use species, only: nspec, spec
    use run_parameters, only: fphi, fapar, fbpar
    use theta_grid, only: ntgrid, theta, delthet, jacob, gradpar, nperiod, bmag
    use theta_grid, only: Rplot, Zplot, aplot, Rprime, Zprime, aprime
    use theta_grid, only: drhodpsi, qval, shape
    use kt_grids, only: naky, ntheta0, theta0, nx, ny, akx, aky
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use dist_fn, only: getan, get_epar, getmoms, par_spectrum
    use dist_fn, only: write_f, write_fyx, def_parity, even
    use dist_fn, only: get_verr, get_gtran, write_poly, collision_error
    use dist_fn_arrays, only: g_adjust
    use collisions, only: vnmult
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: xxf_lo
    use gs2_transforms, only: transform2, inverse2
    use gs2_save, only: gs2_save_for_restart
    use constants
    use gs2_time, only: user_time, user_dt, code_dt
    use gs2_io, only: nc_eigenfunc, nc_final_fields, nc_final_epar, nc_final_an
    use gs2_io, only: nc_final_moments, nc_finish
    use antenna, only: dump_ant_amp
    use splines, only: fitp_surf1, fitp_surf2
    use gs2_heating, only: del_htype, del_hetype
!    use gs2_dist_io, only: write_dist
    implicit none
    integer, intent (in) :: istep
    integer :: ig, ik, it, il, ie, is, unit, ierr
    complex, dimension (:,:,:,:), allocatable :: fcheck_f
!    complex, dimension (:,:,:), allocatable :: xphi, xapar, xbpar
    real, dimension (:), allocatable :: total
    real, dimension (:,:,:), allocatable :: xphi
    real, dimension (:,:,:), allocatable :: bxf, byf, vxf, vyf, bxfsavg, byfsavg
    real, dimension (:,:,:), allocatable :: bxfs, byfs, vxfs, vyfs, rvx, rvy, rx, ry
    complex, dimension (:,:,:), allocatable :: bx, by, vx, vy, vx2, vy2
    complex, dimension (:,:,:), allocatable :: phi2, apar2, bpar2, antot, antota, antotp, epar
    complex, dimension (:,:,:,:), allocatable :: ntot, density, upar, tpar, tperp
    complex, dimension (:,:,:,:), allocatable :: qparflux, pperpj1, qpperpj1
    real, dimension (:), allocatable :: dl_over_b
    complex, dimension (ntheta0, naky) :: phi0, dbfac
    complex, dimension (-ntgrid:ntgrid, ntheta0, naky) :: db
    real, dimension (ntheta0, naky) :: phi02
    real, dimension (2*ntgrid) :: kpar
    real, dimension (:), allocatable :: xx4, yy4, dz
    real, dimension (:,:), allocatable :: bxs, bys, vxs, vys
    real, dimension (:,:), allocatable :: bxsavg, bysavg
    real, dimension (:), allocatable :: stemp, zx1, zxm, zy1, zyn, xx, yy
    real :: zxy11, zxym1, zxy1n, zxymn, L_x, L_y, rxt, ryt, bxt, byt
    integer :: istatus, nnx, nny, nnx4, nny4, ulim, llim, iblock, i, g_unit
    logical :: last = .true.

    if (write_gyx) call write_fyx (phinew,bparnew,last)
    if (write_g) call write_f (last)
    if (write_lpoly) call write_poly (phinew, bparnew, last, istep)
    if (write_cerr) call collision_error (phinew, bparnew, last)
!    if (write_gg) call write_dist (g)

    phi0 = 1.

    if (proc0) then
       if (write_ascii .and. write_parity) then
          call close_output_file (parity_unit)
       end if
       if (write_ascii .and. write_verr) then
          call close_output_file (res_unit)
          call close_output_file (lpc_unit)
          if (write_max_verr) call close_output_file (res_unit2)
       end if
       if (write_ascii) call close_output_file (out_unit)
       if (write_ascii .and. write_cross_phase) call close_output_file (phase_unit)
       if (write_ascii .and. write_hrate) call close_output_file (heat_unit)
       if (write_ascii .and. write_hrate) call close_output_file (heat_unit2)
       if (write_ascii .and. write_jext) call close_output_file (jext_unit)

       if (dump_check1) close (dump_check1_unit)
       if (dump_check2) call close_output_file (dump_check2_unit)

       if (write_eigenfunc) then
          if (write_ascii) call open_output_file (unit, ".eigenfunc")
          phi0 = phi(0,:,:)

          if (def_parity .and. fapar > 0 .and. (.not. even)) phi0 = apar(0, :, :)

          where (abs(phi0) < 10.0*epsilon(0.0)) 
             phi0 = phi(1,:,:)/(theta(1)-theta(0))
          end where

          where (abs(phi0) < 10.0*epsilon(0.0)) 
             phi0 = 1.0
          end where

          if (write_ascii) then
             do ik = 1, naky
                do it = 1, ntheta0
                   do ig = -ntg_out, ntg_out
                      write (unit, "(9(1x,e12.5))") &
                           theta(ig), theta0(it,ik), aky(ik), &
                             phi(ig,it,ik)/phi0(it,ik), &
                            apar(ig,it,ik)/phi0(it,ik), &
                           bpar(ig,it,ik)/phi0(it,ik)
                   end do
                   write (unit, "()")
                end do
             end do
             call close_output_file (unit)
          end if
          call nc_eigenfunc (phi0)
       end if

       if (write_final_fields) then
          if (write_ascii) then
             call open_output_file (unit, ".fields")
             do ik = 1, naky
                do it = 1, ntheta0
                   do ig = -ntg_out, ntg_out
                      write (unit, "(15(1x,e12.5))") &
!                           theta(ig), theta0(it,ik), aky(ik), &
                           theta(ig), aky(ik), akx(it), &
                           phi(ig,it,ik), &
                           apar(ig,it,ik), &
                           bpar(ig,it,ik), &
                           theta(ig) - theta0(it,ik), &
                           cabs(phi(ig,it,ik))
                   end do
                   write (unit, "()")
                end do
             end do
             call close_output_file (unit)
          end if
          call nc_final_fields
       end if
       
       if (write_kpar) then

          allocate (phi2(-ntgrid:ntgrid,ntheta0,naky))   ; phi2 = 0.
          allocate (apar2(-ntgrid:ntgrid,ntheta0,naky))  ; apar2 = 0.
          allocate (bpar2(-ntgrid:ntgrid,ntheta0,naky)) ; bpar2 = 0.

          if (fphi > epsilon(0.0)) then
             call par_spectrum(phi, phi2)
          end if
          if (fapar > epsilon(0.0)) then
             call par_spectrum(apar, apar2)
          endif
          if (fbpar > epsilon(0.0)) then
             call par_spectrum(bpar, bpar2)
          endif

          call open_output_file (unit, ".kpar")
          do ig = 1, ntgrid
             kpar(ig) = (ig-1)*gradpar(ig)/real(2*nperiod-1)
             kpar(2*ntgrid-ig+1)=-(ig)*gradpar(ig)/real(2*nperiod-1)
          end do
          do ik = 1, naky
             do it = 1, ntheta0
                do ig = ntgrid+1,2*ntgrid
                   write (unit, "(9(1x,e12.5))") &
                        kpar(ig), aky(ik), akx(it), &
                        phi2(ig-ntgrid-1,it,ik), &
                        apar2(ig-ntgrid-1,it,ik), &
                        bpar2(ig-ntgrid-1,it,ik)                        
                end do
                do ig = 1, ntgrid
                   write (unit, "(9(1x,e12.5))") &
                        kpar(ig), aky(ik), akx(it), &
                        phi2(ig-ntgrid-1,it,ik), &
                        apar2(ig-ntgrid-1,it,ik), &
                        bpar2(ig-ntgrid-1,it,ik)
                end do
                write (unit, "()")
             end do
          end do
          call close_output_file (unit)
          deallocate (phi2, apar2, bpar2)
       end if

       if (write_final_epar) then
          allocate (epar(-ntgrid:ntgrid,ntheta0,naky))   ; epar = 0.

          call get_epar (phi, apar, phinew, aparnew, epar)
          if (write_ascii) then
             call open_output_file (unit, ".epar")
             do ik = 1, naky
                do it = 1, ntheta0
                   do ig = -ntg_out, ntg_out-1
                      write (unit, "(6(1x,e12.5))") &
!                           theta(ig), theta0(it,ik), aky(ik), &
                           theta(ig), aky(ik), akx(it), &
                           epar(ig,it,ik), &
                           theta(ig) - theta0(it,ik)
                   end do
                   write (unit, "()")
                end do
             end do
             call close_output_file (unit)
          end if
          call nc_final_epar (epar  )
          deallocate (epar)
       end if
   
       if (write_final_db) then  ! definition here assumes we are not using wstar_units
          db = 0.
          do ik = 1, naky
             do it = 1, ntheta0
                dbfac(it,ik) = 1./sum(delthet/bmag/gradpar)/maxval(cabs(phinew(:,it,ik)),1) &
                     * cabs(log(aparnew(1,it,ik)/apar(1,it,ik)))/code_dt
                ig = -ntg_out
                db(ig, it, ik) = aparnew(ig,it,ik)*delthet(ig)/bmag(ig)/gradpar(ig)*dbfac(it,ik)
                do ig = -ntg_out+1, ntg_out-1
                   db(ig, it, ik) = db(ig-1, it, ik) + aparnew(ig,it,ik)*delthet(ig)/bmag(ig)/gradpar(ig)*dbfac(it,ik)
                end do
             end do
          end do

!          db = db * cabs(omega)

          if (write_ascii) then
             call open_output_file (unit, ".db")
             do ik = 1, naky
                do it = 1, ntheta0
                   do ig = -ntg_out, ntg_out-1
                      write (unit, "(5(1x,e12.5))") &
                           theta(ig), aky(ik), akx(it), real(db(ig, it,ik)), aimag(db(ig, it, ik))
                   end do
                   write (unit, "()")
                end do
             end do
             call close_output_file (unit)
          end if
       end if
    end if

    call broadcast (write_final_moments)
    if (write_final_moments) then

       allocate (ntot(-ntgrid:ntgrid,ntheta0,naky,nspec))
       allocate (density(-ntgrid:ntgrid,ntheta0,naky,nspec))
       allocate (upar(-ntgrid:ntgrid,ntheta0,naky,nspec))
       allocate (tpar(-ntgrid:ntgrid,ntheta0,naky,nspec))
       allocate (tperp(-ntgrid:ntgrid,ntheta0,naky,nspec))
       allocate (qparflux(-ntgrid:ntgrid,ntheta0,naky,nspec))
       allocate (pperpj1(-ntgrid:ntgrid,ntheta0,naky,nspec))
       allocate (qpperpj1(-ntgrid:ntgrid,ntheta0,naky,nspec))
       call getmoms (ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)

       if (proc0) then
          if (write_ascii) then
             call open_output_file (unit, ".moments")
             if (.not. write_eigenfunc) phi0 = 1.
             do is  = 1, nspec
                do ik = 1, naky
                   do it = 1, ntheta0
                      do ig = -ntg_out, ntg_out
                         ! TEMP FOR TESTING -- MAB
!                         write (unit, *) &
!                              real(ntot(ig,it,ik,is)/phi0(it,ik)), &
!                              real(upar(ig,it,ik,is)/phi0(it,ik)), &
!                              aimag(upar(ig,it,ik,is)/phi0(it,ik)), &
!                              real(tperp(ig,it,ik,is)/phi0(it,ik))
                         write (unit, "(15(1x,e12.5))") &
                              theta(ig), aky(ik), akx(it), &
                              ntot(ig,it,ik,is)/phi0(it,ik), &
                              density(ig,it,ik,is)/phi0(it,ik), &
                              upar(ig,it,ik,is)/phi0(it,ik), &
                              tpar(ig,it,ik,is)/phi0(it,ik), &
                              tperp(ig,it,ik,is)/phi0(it,ik), &
                              theta(ig) - theta0(it,ik), &
                              real(is)
                      end do
                      write (unit, "()")
                   end do
                end do
             end do
             call close_output_file (unit)          
          end if
          call nc_final_moments (ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)

          if (write_ascii) then
             call open_output_file (unit, ".mom2")
             if (.not. write_eigenfunc) phi0 = 1.
             phi02=real(phi0*conjg(phi0))
             do is  = 1, nspec
                do ik = 1, naky
                   do it = 1, ntheta0
                      do ig = -ntg_out, ntg_out
                         write (unit, "(15(1x,e12.5))") &
                              theta(ig), aky(ik), akx(it), &
                              real(ntot(ig,it,ik,is)*conjg(ntot(ig,it,ik,is)))/phi02(it,ik), &
                              real(density(ig,it,ik,is)*conjg(density(ig,it,ik,is)))/phi02(it,ik), &
                              real(upar(ig,it,ik,is)*conjg(upar(ig,it,ik,is)))/phi02(it,ik), &
                              real(tpar(ig,it,ik,is)*conjg(tpar(ig,it,ik,is)))/phi02(it,ik), &
                              real(tperp(ig,it,ik,is)*conjg(tperp(ig,it,ik,is)))/phi02(it,ik), &
                              theta(ig) - theta0(it,ik), &
                              real(is)
                      end do
                      write (unit, "()")
                   end do
                end do
             end do

             call close_output_file (unit)          
             call open_output_file (unit, ".amoments")
             write (unit,*) 'type    kx     re(phi)    im(phi)    re(ntot)   im(ntot)   ',&
                  &'re(dens)   im(dens)   re(upar)   im(upar)   re(tpar)',&
                  &'   im(tpar)   re(tperp)  im(tperp)'
             
             allocate (dl_over_b(-ntg_out:ntg_out))
             
             dl_over_b = delthet(-ntg_out:ntg_out)*jacob(-ntg_out:ntg_out)
             dl_over_b = dl_over_b / sum(dl_over_b)
             
             do is  = 1, nspec
                do it = 2, ntheta0/2+1
                   write (unit, "(i2,14(1x,e10.3))") spec(is)%type, akx(it), &
                        sum(phinew(-ntg_out:ntg_out,it,1)*dl_over_b), &
                        sum(ntot(-ntg_out:ntg_out,it,1,is)*dl_over_b), &
                        sum(density(-ntg_out:ntg_out,it,1,is)*dl_over_b), &
                        sum(upar(-ntg_out:ntg_out,it,1,is)*dl_over_b), &
                        sum(tpar(-ntg_out:ntg_out,it,1,is)*dl_over_b), &
                        sum(tperp(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                end do
             end do
             deallocate (dl_over_b)
             call close_output_file (unit)          
          end if
       end if
       deallocate (ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)
    end if

    if (write_final_antot) then
       
       allocate ( antot(-ntgrid:ntgrid,ntheta0,naky)) ; antot = 0.
       allocate (antota(-ntgrid:ntgrid,ntheta0,naky)) ; antota = 0. 
       allocate (antotp(-ntgrid:ntgrid,ntheta0,naky)) ; antotp = 0.
       call getan (antot, antota, antotp)
       if (proc0) then
          if (write_ascii) then
             call open_output_file (unit, ".antot")
             do ik = 1, naky
                do it = 1, ntheta0
                   do ig = -ntg_out, ntg_out
                      write (unit, "(10(1x,e12.5))") &
                           theta(ig), theta0(it,ik), aky(ik), &
                           antot(ig,it,ik), &
                           antota(ig,it,ik), &
                           antotp(ig,it,ik), &
                           theta(ig) - theta0(it,ik)
                   end do
                   write (unit, "()")
                end do
             end do
             call close_output_file (unit)
          end if
          call nc_final_an (antot, antota, antotp)
       end if
       deallocate (antot, antota, antotp)
    end if

    if (save_for_restart) then
       call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, &
            fphi, fapar, fbpar, .true.)
    end if

    !<DD> Added for saving distribution function
    if (save_distfn) then
    	!Convert h to distribution function
    	call g_adjust(gnew,phinew,bparnew,fphi,fbpar)
    	
    	!Save dfn, fields and velocity grids to file
       	call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, &
          	fphi, fapar, fbpar, .true.,.true.)
    	
        !Convert distribution function back to h
        call g_adjust(gnew,phinew,bparnew,-fphi,-fbpar)
    end if

    !</DD> Added for saving distribution function
    
    call nc_finish

    if (proc0) call dump_ant_amp

    if (write_gs) then
       nny = 2*ny
       nnx = 2*nx
       nnx4 = nnx+4
       nny4 = nny+4

       allocate (bx(-ntgrid:ntgrid,ntheta0,naky))
       allocate (by(-ntgrid:ntgrid,ntheta0,naky))
       allocate (vx(-ntgrid:ntgrid,ntheta0,naky))
       allocate (vy(-ntgrid:ntgrid,ntheta0,naky))

       do ik=1,naky
          do it=1,ntheta0
             do ig=-ntgrid, ntgrid
                bx(ig,it,ik) =  zi*aky(ik)*apar(ig,it,ik)
                by(ig,it,ik) = -zi*akx(it)*apar(ig,it,ik)
                vx(ig,it,ik) = -zi*aky(ik)*phi(ig,it,ik)
                vy(ig,it,ik) =  zi*akx(it)*phi(ig,it,ik)
             end do
          end do
       end do

       allocate (bxf(nnx,nny,-ntgrid:ntgrid))
       allocate (byf(nnx,nny,-ntgrid:ntgrid))
       allocate (vxf(nnx,nny,-ntgrid:ntgrid))
       allocate (vyf(nnx,nny,-ntgrid:ntgrid))

       call transform2 (bx, bxf, nny, nnx)
       call transform2 (by, byf, nny, nnx)
       call transform2 (vx, vxf, nny, nnx)
       call transform2 (vy, vyf, nny, nnx)
       
       ! fields come out as (x, y, z)

       deallocate (bx, by)

       allocate (   bxfs(nnx4, nny4, -ntgrid:ntgrid))
       allocate (   byfs(nnx4, nny4, -ntgrid:ntgrid))
       allocate (bxfsavg(nnx4, nny4, -ntgrid:ntgrid))
       allocate (byfsavg(nnx4, nny4, -ntgrid:ntgrid))
       allocate (   vxfs(nnx4, nny4, -ntgrid:ntgrid))
       allocate (   vyfs(nnx4, nny4, -ntgrid:ntgrid))

       do ig=-ntgrid,ntgrid
          do ik=1,2
             do it=3,nnx4-2
                bxfs(it,ik,ig) = bxf(it-2,nny-2+ik,ig)
                byfs(it,ik,ig) = byf(it-2,nny-2+ik,ig)
                vxfs(it,ik,ig) = vxf(it-2,nny-2+ik,ig)
                vyfs(it,ik,ig) = vyf(it-2,nny-2+ik,ig)

                bxfs(it,nny4-2+ik,ig) = bxf(it-2,ik,ig)
                byfs(it,nny4-2+ik,ig) = byf(it-2,ik,ig)
                vxfs(it,nny4-2+ik,ig) = vxf(it-2,ik,ig)
                vyfs(it,nny4-2+ik,ig) = vyf(it-2,ik,ig)
             end do
          end do
          do ik=3,nny4-2
             do it=3,nnx4-2
                bxfs(it,ik,ig) = bxf(it-2,ik-2,ig)
                byfs(it,ik,ig) = byf(it-2,ik-2,ig)
                vxfs(it,ik,ig) = vxf(it-2,ik-2,ig)
                vyfs(it,ik,ig) = vyf(it-2,ik-2,ig)
             end do
          end do
          do ik=1,nny4
             do it=1,2
                bxfs(it,ik,ig) = bxfs(nnx4-4+it,ik,ig)
                byfs(it,ik,ig) = byfs(nnx4-4+it,ik,ig)
                vxfs(it,ik,ig) = vxfs(nnx4-4+it,ik,ig)
                vyfs(it,ik,ig) = vyfs(nnx4-4+it,ik,ig)

                bxfs(nnx4-2+it,ik,ig) = bxfs(it+2,ik,ig)
                byfs(nnx4-2+it,ik,ig) = byfs(it+2,ik,ig)
                vxfs(nnx4-2+it,ik,ig) = vxfs(it+2,ik,ig)
                vyfs(nnx4-2+it,ik,ig) = vyfs(it+2,ik,ig)
             end do
          end do
       end do

       deallocate (vxf, vyf)

       allocate (xx4(nnx4), xx(nnx))
       allocate (yy4(nny4), yy(nny))
       
       L_x = 2.0*pi/akx(2)
       L_y = 2.0*pi/aky(2)

       do it = 1, nnx
          xx4(it+2) = real(it-1)*L_x/real(nnx)
          xx(it) = real(it-1)*L_x/real(nnx)
       end do
       do it=1,2
          xx4(it) = xx4(nnx4-4+it)-L_x
          xx4(nnx4-2+it) = xx4(it+2)+L_x
       end do

       do ik = 1, nny
          yy4(ik+2) = real(ik-1)*L_y/real(nny)
          yy(ik)    = real(ik-1)*L_y/real(nny)
       end do
       do ik=1,2
          yy4(ik) = yy4(nny4-4+ik)-L_y
          yy4(nny4-2+ik) = yy4(ik+2)+L_y
       end do

       allocate (dz(-ntgrid:ntgrid))
       dz = delthet*jacob

       allocate (   bxs(3*nnx4*nny4,-ntgrid:ntgrid)) ; bxs = 0.
       allocate (   bys(3*nnx4*nny4,-ntgrid:ntgrid)) ; bys = 0.
       allocate (   vxs(3*nnx4*nny4,-ntgrid:ntgrid)) ; vxs = 0.
       allocate (   vys(3*nnx4*nny4,-ntgrid:ntgrid)) ; vys = 0.

       allocate (bxsavg(3*nnx4*nny4,-ntgrid:ntgrid))
       allocate (bysavg(3*nnx4*nny4,-ntgrid:ntgrid))

       allocate (stemp(nnx4+2*nny4))
       allocate (zx1(nny4), zxm(nny4), zy1(nnx4), zyn(nnx4))

       do ig=-ntgrid, ntgrid
          call fitp_surf1(nnx4, nny4, xx, yy, bxfs(:,:,ig), &
               nnx4, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
               255, bxs(:,ig), stemp, 1., ierr)

          call fitp_surf1(nnx4, nny4, xx, yy, byfs(:,:,ig), &
               nnx4, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
               255, bys(:,ig), stemp, 1., ierr)

          call fitp_surf1(nnx4, nny4, xx, yy, vxfs(:,:,ig), &
               nnx4, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
               255, vxs(:,ig), stemp, 1., ierr)

          call fitp_surf1(nnx4, nny4, xx, yy, vyfs(:,:,ig), &
               nnx4, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
               255, vys(:,ig), stemp, 1., ierr)
       end do

       deallocate (zx1, zxm, zy1, zyn, stemp)

       do ig=-ntgrid, ntgrid-1
          bxsavg(:,ig) = 0.5*(bxs(:,ig)+bxs(:,ig+1))
          bysavg(:,ig) = 0.5*(bys(:,ig)+bys(:,ig+1))

          bxfsavg(:,:,ig) = 0.5*(bxfs(:,:,ig)+bxfs(:,:,ig+1))
          byfsavg(:,:,ig) = 0.5*(byfs(:,:,ig)+byfs(:,:,ig+1))
       end do

       ! now, integrate to find a field line

       allocate ( rx(nnx,nny,-ntgrid:ntgrid))
       allocate ( ry(nnx,nny,-ntgrid:ntgrid))
       allocate (rvx(nnx,nny,-ntgrid:ntgrid)) ; rvx = 0.
       allocate (rvy(nnx,nny,-ntgrid:ntgrid)) ; rvy = 0.

       do ik=1,nny
          do it=1,nnx
             rx(it,ik,-ntgrid) = xx(it)
             ry(it,ik,-ntgrid) = yy(ik)
          end do
       end do

       iblock = (nnx*nny-1)/nproc + 1
       llim = 1 + iblock * iproc
       ulim = min(nnx*nny, llim+iblock-1)

       do i=llim, ulim
          it = 1 + mod(i-1, nnx)
          ik = 1 + mod((i-1)/nnx, nny)
          
          ig = -ntgrid
          
          rxt = rx(it,ik,ig)
          ryt = ry(it,ik,ig)
          
          rvx(it,ik,ig) = fitp_surf2(rxt, ryt, nnx4, nny4, xx4, yy4, vxfs(:,:,ig), nnx4, vxs(:,ig), 1.)
          rvy(it,ik,ig) = fitp_surf2(rxt, ryt, nnx4, nny4, xx4, yy4, vyfs(:,:,ig), nnx4, vys(:,ig), 1.)
          
          do ig=-ntgrid,ntgrid-1
             
             bxt = fitp_surf2(rxt, ryt, nnx4, nny4, xx4, yy4, bxfs(:,:,ig), nnx4, bxs(:,ig), 1.)
             byt = fitp_surf2(rxt, ryt, nnx4, nny4, xx4, yy4, byfs(:,:,ig), nnx4, bys(:,ig), 1.)
             
             rxt = rx(it,ik,ig) + 0.5*dz(ig)*bxt  
             ryt = ry(it,ik,ig) + 0.5*dz(ig)*byt  
             
             if (rxt > L_x) rxt = rxt - L_x
             if (ryt > L_y) ryt = ryt - L_y
             
             if (rxt < 0.) rxt = rxt + L_x
             if (ryt < 0.) ryt = ryt + L_y
             
             bxt = fitp_surf2(rxt, ryt, nnx4, nny4, xx4, yy4, bxfsavg(:,:,ig), nnx4, bxsavg(:,ig), 1.)
             byt = fitp_surf2(rxt, ryt, nnx4, nny4, xx4, yy4, byfsavg(:,:,ig), nnx4, bysavg(:,ig), 1.)
             
             rxt = rx(it,ik,ig) + dz(ig)*bxt  
             ryt = ry(it,ik,ig) + dz(ig)*byt  
             
             if (rxt > L_x) rxt = rxt - L_x
             if (ryt > L_y) ryt = ryt - L_y
             
             if (rxt < 0.) rxt = rxt + L_x
             if (ryt < 0.) ryt = ryt + L_y
             
             rx(it,ik,ig+1) = rxt
             ry(it,ik,ig+1) = ryt
             
             rvx(it,ik,ig+1) = fitp_surf2(rxt, ryt, nnx4, nny4, xx4, yy4, vxfs(:,:,ig+1), nnx4, vxs(:,ig+1), 1.)
             rvy(it,ik,ig+1) = fitp_surf2(rxt, ryt, nnx4, nny4, xx4, yy4, vyfs(:,:,ig+1), nnx4, vys(:,ig+1), 1.)
          end do
       end do

       deallocate (bxfs, byfs, bxfsavg, byfsavg, vxfs, vyfs)
       deallocate (rx, ry, bxs, bys, vxs, vys, bxsavg, bysavg)

       allocate (total(2*nnx*nny*(2*ntgrid+1)))

       i=1
       do ig=-ntgrid,ntgrid
          do ik=1,nny
             do it=1,nnx
                total(i) = rvx(it,ik,ig)
                total(i+1) = rvy(it,ik,ig)
                i = i + 2
             end do
          end do
       end do
       
       call sum_reduce(total, 0)

       i=1
       do ig=-ntgrid,ntgrid
          do ik=1,nny
             do it=1,nnx
                rvx(it,ik,ig) = total(i)
                rvy(it,ik,ig) = total(i+1)
                i = i + 2
             end do
          end do
       end do
       
       if (proc0) then
          call inverse2 (rvx, vx, nny, nnx)
          call inverse2 (rvy, vy, nny, nnx)
       
          allocate (vx2(-ntgrid:ntgrid,ntheta0,naky))
          allocate (vy2(-ntgrid:ntgrid,ntheta0,naky))

          call par_spectrum (vx, vx2)
          call par_spectrum (vy, vy2)

          call open_output_file (unit, ".gs")
          do ig = 1, ntgrid
             kpar(ig) = (ig-1)*gradpar(ig)/real(2*nperiod-1)
             kpar(2*ntgrid-ig+1)=-(ig)*gradpar(ig)/real(2*nperiod-1)
          end do
          do ik = 1, naky
             do it = 1, ntheta0
                do ig = ntgrid+1,2*ntgrid
                   write (unit, "(9(1x,e12.5))") &
                        kpar(ig), aky(ik), akx(it), &
                        real(vx2(ig-ntgrid-1,it,ik)), &
                        real(vy2(ig-ntgrid-1,it,ik))
                end do
                do ig = 1, ntgrid
                   write (unit, "(9(1x,e12.5))") &
                        kpar(ig), aky(ik), akx(it), &
                        real(vx2(ig-ntgrid-1,it,ik)), &
                        real(vy2(ig-ntgrid-1,it,ik))
                end do
                write (unit, "()")
             end do
          end do
          call close_output_file (unit)
          deallocate (vx2, vy2)
       end if

       deallocate (vx, vy, rvx, rvy)
    
    end if
    
    if (proc0) then
       call open_output_file (g_unit, ".g")
       write (g_unit,fmt="('# shape: ',a)") trim(shape)
       write (g_unit,fmt="('# q = ',e10.4,' drhodpsi = ',e10.4)") qval, drhodpsi
       write (g_unit,fmt="('# theta1             R2                  Z3               alpha4      ', &
            &   '       Rprime5              Zprime6           alpha_prime7 ')")
       do i=-ntgrid,ntgrid
          write (g_unit,1000) theta(i),Rplot(i),Zplot(i),aplot(i), &
               Rprime(i),Zprime(i),aprime(i)
       enddo
       call close_output_file (g_unit)
    end if

    if (write_hrate) then
       call del_htype (h)
       call del_htype (h_hist)
       call del_htype (hk_hist)
       call del_htype (hk)
    end if
    if (write_hrate_e) then
       call del_hetype (he)
       call del_hetype (he_hist)
    end if
    if (allocated(h_hist)) deallocate (h_hist, hk_hist, hk)
    if (allocated(he_hist)) deallocate (he_hist, he)
    if (allocated(j_ext_hist)) deallocate (j_ext_hist)
    if (allocated(omegahist)) deallocate (omegahist)
    if (allocated(pflux)) deallocate (pflux, qheat, vflux, vflux_par, vflux_perp, pmflux, qmheat, vmflux, &
         pbflux, qbheat, vbflux, vflux0, vflux1, exchange)
    if (allocated(pflux_emu)) deallocate(pflux_emu,mflux_emu,bflux_emu)
    if (allocated(pflux_e)) deallocate(pflux_e,mflux_e,bflux_e,trapped_flux,passing_flux)
    if (allocated(peflux)) deallocate(peflux,meflux,beflux)
    if (allocated(bxf)) deallocate (bxf, byf, xx4, xx, yy4, yy, dz, total)
    if (allocated(pflux_avg)) deallocate (pflux_avg, qflux_avg, heat_avg, vflux_avg)

    wtmp_old = 0. ; nout = 1 ; nout_movie = 1 ; nout_big = 1
    initialized = .false.

1000  format(20(1x,1pg18.11))
  end subroutine finish_gs2_diagnostics

  subroutine loop_diagnostics (istep, exit, debopt)
    use species, only: nspec, spec, has_electron_species
    use theta_grid, only: theta, ntgrid, delthet, jacob
    use theta_grid, only: gradpar, nperiod
    use kt_grids, only: naky, ntheta0, theta0, aky, akx
    use kt_grids, only: nkpolar, jtwist_out !, akpolar
    use run_parameters, only: woutunits, tunits, fapar, fphi, fbpar, eqzip
!    use run_parameters, only: nstep, include_lowflow
    use run_parameters, only: nstep
    use fields, only: phinew, aparnew, bparnew, phi,apar,bpar
    use fields, only: kperp, fieldlineavgphi, phinorm
    use dist_fn, only: flux, write_f, write_fyx, flux_emu, flux_e, eflux
    use dist_fn, only: omega0, gamma0, getmoms, par_spectrum
    use dist_fn, only: get_verr, get_gtran, write_poly, collision_error
    use dist_fn, only: getmoms_notgc, lf_flux, eexchange
    use dist_fn, only: flux_vs_theta_vs_vpa
    use dist_fn_arrays, only: g, gnew, aj0, vpa, g_adjust
    use collisions, only: ncheck, vnmult, vary_vnew
    use mp, only: proc0, broadcast, iproc, send, receive
    use file_utils, only: get_unused_unit, flush_output_file
    use prof, only: prof_entering, prof_leaving
    use gs2_time, only: user_time
    use gs2_io, only: nc_qflux, nc_vflux, nc_pflux, nc_loop, nc_loop_moments, nc_flux_emu, nc_flux_e, nc_eflux
    use gs2_io, only: nc_loop_fullmom, nc_loop_sym, nc_loop_corr, nc_loop_corr_extend
    use gs2_io, only: nc_loop_vres
    use gs2_io, only: nc_loop_movie, nc_write_fields, nc_write_moments
    use gs2_layouts, only: yxf_lo, g_lo
! MAB>
    use gs2_layouts, only: init_parity_layouts, idx, idx_local, proc_id
    use gs2_layouts, only: p_lo, ie_idx, is_idx, il_idx, it_idx, ik_idx
! <MAB
    use gs2_transforms, only: init_transforms, transform2
    use le_grids, only: nlambda, ng2, integrate_moment, negrid, integrate_kysum
    use nonlinear_terms, only: nonlin
    use antenna, only: antenna_w
    use gs2_heating, only: heating_diagnostics, del_htype
    use constants
    use parameter_scan_arrays, only: scan_hflux => hflux_tot 
    use parameter_scan_arrays, only: scan_momflux => momflux_tot 
    use parameter_scan_arrays, only: scan_phi2_tot => phi2_tot 
    use parameter_scan_arrays, only: scan_nout => nout

    implicit none
!    integer :: nout = 1
!    integer :: nout_movie = 1
    integer, intent (in) :: istep
    logical, intent (out) :: exit
    real, dimension(:,:,:), allocatable :: yxphi, yxapar, yxbpar
    complex, dimension (ntheta0, naky) :: omega, omegaavg

    !GGH J_external
    real, dimension(:,:), allocatable ::  j_ext

    real, dimension (ntheta0, naky) :: phitot, akperp
    complex, dimension (ntheta0, naky, nspec) :: pfluxneo,qfluxneo
    real :: phi2, apar2, bpar2
    real, dimension (ntheta0, naky) :: phi2_by_mode, apar2_by_mode, bpar2_by_mode
    real, dimension (ntheta0, naky, nspec) :: ntot2_by_mode, ntot20_by_mode
    real, dimension (ntheta0, naky, nspec) :: tpar2_by_mode, tperp2_by_mode
!    real, dimension (:,:,:,:), allocatable :: errest_by_mode
    integer, dimension (:,:), allocatable :: erridx
    real, dimension (:,:), allocatable :: errest
!MAB> arrays needed for parity diagnostic
    integer :: iplo, iglo, sgn2, isgn, il, ie
    real, dimension (:,:,:), allocatable :: vflx_sym
    complex, dimension (:,:,:,:), allocatable :: gparity, gmx, gpx
    complex, dimension (:,:,:), allocatable :: g0, gm, gp
    complex, dimension (:,:,:), allocatable :: g_kint, gm_kint, gp_kint
    complex, dimension (:,:), allocatable :: g_avg, gnorm_avg, phim
    complex, dimension (:), allocatable :: g_all_tot, g_nokx_tot, g_nosig_tot, gtmp
    complex, dimension (:), allocatable :: gnorm_all_tot, gnorm_nokx_tot, gnorm_nosig_tot
    real, dimension (:,:,:), allocatable :: gmnorm, gmint, gpint, gmnormint, gmavg, gpavg, gmnormavg
    real, dimension (:), allocatable :: gmtot, gptot, gtot, gmnormtot
    real, dimension (:), allocatable :: gm_nokx_tot, gp_nokx_tot, gmnorm_nokx_tot
    real, dimension (:), allocatable :: gm_nosig_tot, gp_nosig_tot, gmnorm_nosig_tot
    logical :: first = .true.
    ! arrays needed for parallel autocorrelation diagnostic
    real, save :: tcorr0 = 0.0
    real, dimension (:,:,:), allocatable :: phi2_extend
    complex, dimension (:,:,:), allocatable :: phi_corr
    complex, dimension (:,:), allocatable :: phi_corr_2pi
    real, dimension (:,:,:), allocatable, save :: phiextend_sum
    complex, dimension (:,:,:), allocatable, save :: phicorr_sum
    complex, dimension (:,:), allocatable, save :: phicorr_2pi_sum
!<MAB
    real :: geavg, glavg, gtavg
    real :: t, denom
    integer :: ig, ik, it, is, unit, i, j, nnx, nny, ifield, write_mod
    complex :: sourcefac
!    complex :: phiavg
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky,nspec) :: ntot, density, &
         upar, tpar, tperp, qparflux, pperpj1, qpperpj1, upartot, uperptot, ttot
    complex, dimension (ntheta0, nspec) :: ntot00, density00, upar00, tpar00, tperp00
    complex, dimension (ntheta0) :: phi00
    complex, save :: wtmp_new
!    complex :: wtmp_old = 0.
    real, dimension (:), allocatable :: dl_over_b
    real, dimension (ntheta0, nspec) :: x_qmflux
    real, dimension (nspec) :: ntot2, ntot20, tpar2, tperp2
    real, dimension (nspec) ::  heat_fluxes,  part_fluxes, mom_fluxes, parmom_fluxes, perpmom_fluxes
    real, dimension (nspec) :: lfmom_fluxes, vflux1_avg  ! low-flow correction to turbulent momentum fluxes
    real, dimension (nspec) :: mheat_fluxes, mpart_fluxes, mmom_fluxes
    real, dimension (nspec) :: bheat_fluxes, bpart_fluxes, bmom_fluxes
    real, dimension (nspec) :: energy_exchange
    real, dimension (nspec) ::  heat_par,  heat_perp
    real, dimension (nspec) :: mheat_par, mheat_perp
    real, dimension (nspec) :: bheat_par, bheat_perp
    real, dimension (naky) :: fluxfac
    real :: phase_tot, phase_theta
!    real, dimension (:), allocatable :: phi_by_k, apar_by_k, bpar_by_k
    real :: hflux_tot, zflux_tot, vflux_tot
    real, dimension(nspec) :: tprim_tot, fprim_tot
    real, save :: t_old = 0.
    character(200) :: filename
    logical :: last = .false.
    logical,optional:: debopt
    logical:: debug=.false.

    if (present(debopt)) debug=debopt

    part_fluxes = 0.0 ; mpart_fluxes = 0.0 ; bpart_fluxes = 0.0
    heat_fluxes = 0.0 ; mheat_fluxes = 0.0 ; bheat_fluxes = 0.0
    mom_fluxes = 0.0 ; mmom_fluxes = 0.0 ; bmom_fluxes = 0.0
    energy_exchange = 0.0

    phase_tot = 0.0 ;  phase_theta = 0.0

    call prof_entering ("loop_diagnostics")

    exit = .false.

!    if (eqzip .or. .not. nonlin) then
! MR, 10/3/2009: avoid calling get_omegaavg in nonlinear calculations
       if (proc0) then
          if (debug) write(6,*) "loop_diagnostics: proc0 call get_omegaavg"
          call get_omegaavg (istep, exit, omegaavg, debug)
          if (debug) write(6,*) "loop_diagnostics: proc0 done called get_omegaavg"
       endif
       call broadcast (exit)
!    endif

    if (write_hrate) call heating (istep, h, hk)

    if (write_hrate_e) call heating_e (istep, he)

!>GGH
    !Write Jexternal vs. time
    if (write_jext) then
       allocate (j_ext(ntheta0, naky)); j_ext=0.
       call calc_jext(istep,j_ext)
    endif
!<GGH

    call prof_leaving ("loop_diagnostics")
if (debug) write(6,*) "loop_diagnostics: call update_time"

    if (make_movie .and. mod(istep,nmovie)==0) then
       t = user_time
       ! EAB 09/17/03 -- modify dump_fields_periodically to print out inverse fft of fields in x,y
       ! write(*,*) "iproc now in dump_fields_periodically case", iproc 
       nnx = yxf_lo%nx
       nny = yxf_lo%ny
       if (fphi > epsilon(0.0)) then
          allocate (yxphi(nnx,nny,-ntgrid:ntgrid))
          call getmoms (ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)
!          call transform2 (ntot, yxphi, nny, nnx)
          call transform2 (phinew, yxphi, nny, nnx)
       end if
!       if (fapar > epsilon(0.0)) then
          allocate (yxapar(nnx,nny,-ntgrid:ntgrid))
          call transform2 (ntot, yxapar, nny, nnx)
!          call transform2 (aparnew, yxapar, nny, nnx)
!       end if
       if (fbpar > epsilon(0.0)) then 
          allocate (yxbpar(nnx,nny,-ntgrid:ntgrid))
          call transform2 (bparnew, yxbpar, nny, nnx)
       end if

       if (proc0) then
          call nc_loop_movie(nout_movie, t, yxphi, yxapar, yxbpar)
       end if

       if (fphi > epsilon(0.0)) deallocate (yxphi)
       if (fapar > epsilon(0.0)) deallocate (yxapar)
       if (fbpar > epsilon(0.0)) deallocate (yxbpar)
       nout_movie = nout_movie + 1

    end if

    if (write_gyx .and. mod(istep,nmovie) == 0) call write_fyx (phinew,bparnew,last)

    if (vary_vnew) then
       write_mod = mod(istep,ncheck)
    else
       write_mod = mod(istep,nwrite)
    end if

    if (write_max_verr) write_verr = .true.
    if (write_verr .and. write_mod == 0) then

       allocate(errest(5,2), erridx(5,3))

       errest = 0.0; erridx = 0

! error estimate obtained by comparing standard integral with less-accurate integral
       call get_verr (errest, erridx, phinew, bparnew)

! error estimate based on monitoring amplitudes of legendre polynomial coefficients
       call get_gtran (geavg, glavg, gtavg, phinew, bparnew, istep)

       if (proc0) then

! write error estimates to .nc file          
!          call nc_loop_vres (nout, errest_by_mode, lpcoef_by_mode)

! write error estimates for ion dist. fn. at outboard midplane with ik=it=1 to ascii files
          if (write_ascii) then
             t = user_time
             if (nlambda - ng2 > 1) then
                write(lpc_unit,"(4(1x,e12.6))") t, geavg, glavg, gtavg
             else
                write(lpc_unit,"(3(1x,e12.6))") t, geavg, glavg
             end if
             write(res_unit,"(8(1x,e12.6))") t, errest(1,2), errest(2,2), errest(3,2), &
                  errest(4,2), errest(5,2), vnmult(1)*spec(1)%vnewk, vnmult(2)*spec(1)%vnewk
             if (write_max_verr) then
                write(res_unit2,"(3(i8),(1x,e12.6),3(i8),(1x,e12.6),3(i8),(1x,e12.6),3(i8),(1x,e12.6),3(i8),(1x,e12.6))") &
                     erridx(1,1), erridx(1,2), erridx(1,3), errest(1,1), &
                     erridx(2,1), erridx(2,2), erridx(2,3), errest(2,1), &
                     erridx(3,1), erridx(3,2), erridx(3,3), errest(3,1), &
                     erridx(4,1), erridx(4,2), erridx(4,3), errest(4,1), &
                     erridx(5,1), erridx(5,2), erridx(5,3), errest(5,1)
             end if
          end if
       end if
       deallocate(errest,erridx)
    end if

    if (mod(istep,nwrite) /= 0 .and. .not. exit) return
    t = user_time

    if (write_g) call write_f (last)
    if (write_lpoly) call write_poly (phinew, bparnew, last, istep)

    call prof_entering ("loop_diagnostics-1")
if (debug) write(6,*) "loop_diagnostics: -1"

    if (proc0) then
       omega = omegahist(mod(istep,navg),:,:)
       sourcefac = exp(-zi*omega0*t+gamma0*t)
       call phinorm (phitot)
       call get_vol_average (phinew, phinew, phi2, phi2_by_mode)
       if (fapar > epsilon(0.0)) then
          call get_vol_average (aparnew, aparnew, apar2, apar2_by_mode)
          apar2 = apar2
          apar2_by_mode = apar2_by_mode
       end if
       if (fbpar > epsilon(0.0)) then
          call get_vol_average (bparnew, bparnew, bpar2, bpar2_by_mode)
          bpar2 = bpar2
          bpar2_by_mode = bpar2_by_mode
       end if
    end if

    if (write_any_fluxes) then
       call g_adjust (gnew, phinew, bparnew, fphi, fbpar)
       if (write_flux_emu) then
          call flux_emu (phinew, aparnew, bparnew, pflux_emu, mflux_emu, bflux_emu)
          do is = 1,nspec
             pflux_emu(:,:,:,is) = pflux_emu(:,:,:,is) * spec(is)%dens
             mflux_emu(:,:,:,is) = mflux_emu(:,:,:,is) * spec(is)%dens
             bflux_emu(:,:,:,is) = bflux_emu(:,:,:,is) * spec(is)%dens
          end do
       end if
       if (write_flux_e) then
          call flux_e (phinew, aparnew, bparnew, pflux_e, mflux_e, bflux_e,trapped_flux,passing_flux)
          do is = 1,nspec
             pflux_e(:,is) = pflux_e(:,is) * spec(is)%dens
             mflux_e(:,is) = mflux_e(:,is) * spec(is)%dens
             bflux_e(:,is) = bflux_e(:,is) * spec(is)%dens
             trapped_flux(:,is) = trapped_flux(:,is) * spec(is)%dens
             passing_flux(:,is) = passing_flux(:,is) * spec(is)%dens
          end do
       end if
       if (write_eflux) then
          call eflux (phi,apar,bpar,phinew, aparnew, bparnew, peflux, meflux, beflux)
          do is = 1,nspec
             peflux(:,is) = peflux(:,is) * spec(is)%dens
             meflux(:,is) = meflux(:,is) * spec(is)%dens
             beflux(:,is) = beflux(:,is) * spec(is)%dens
          end do
       end if

       call flux (phinew, aparnew, bparnew, &
            pflux,  qheat,  vflux, vflux_par, vflux_perp, &
            pmflux, qmheat, vmflux, pbflux, qbheat, vbflux)
#ifdef LOWFLOW
       ! lowflow terms only implemented in electrostatic limit at present
       call lf_flux (phinew, vflux0, vflux1)
#endif
       call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar)
       call eexchange (phinew, exchange)

       if (proc0) then
          if (fphi > epsilon(0.0)) then
             do is = 1, nspec
                qheat(:,:,is,1) = qheat(:,:,is,1) * spec(is)%dens*spec(is)%temp
                call get_volume_average (qheat(:,:,is,1), heat_fluxes(is))

                qheat(:,:,is,2) = qheat(:,:,is,2) * spec(is)%dens*spec(is)%temp
                call get_volume_average (qheat(:,:,is,2), heat_par(is))

                qheat(:,:,is,3) = qheat(:,:,is,3) * spec(is)%dens*spec(is)%temp
                call get_volume_average (qheat(:,:,is,3), heat_perp(is))
                
                pflux(:,:,is) = pflux(:,:,is) * spec(is)%dens
                call get_volume_average (pflux(:,:,is), part_fluxes(is))

                vflux(:,:,is) = vflux(:,:,is) * spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
                call get_volume_average (vflux(:,:,is), mom_fluxes(is))

                vflux_par(:,:,is) = vflux_par(:,:,is) * spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
                call get_volume_average (vflux_par(:,:,is), parmom_fluxes(is))

                vflux_perp(:,:,is) = vflux_perp(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%stm
                call get_volume_average (vflux_perp(:,:,is), perpmom_fluxes(is))

                exchange(:,:,is) = exchange(:,:,is) * spec(is)%dens*spec(is)%z
                call get_volume_average (exchange(:,:,is), energy_exchange(is))

#ifdef LOWFLOW
                vflux0(:,:,is) = vflux0(:,:,is) * spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
                call get_volume_average (vflux0(:,:,is), lfmom_fluxes(is))

                vflux1(:,:,is) = vflux1(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%temp/spec(is)%z
                call get_volume_average (vflux1(:,:,is), vflux1_avg(is))

! TMP UNTIL VFLUX0 IS TESTED
!                   mom_fluxes = mom_fluxes + lfmom_fluxes
#endif
             end do
          end if
          if (fapar > epsilon(0.0)) then
             do is = 1, nspec
                qmheat(:,:,is,1)=qmheat(:,:,is,1) * spec(is)%dens*spec(is)%temp
                call get_volume_average (qmheat(:,:,is,1), mheat_fluxes(is))

                call get_surf_average (qmheat(:,:,is,1), x_qmflux(:,is))

                qmheat(:,:,is,2)=qmheat(:,:,is,2) * spec(is)%dens*spec(is)%temp
                call get_volume_average (qmheat(:,:,is,2), mheat_par(is))

                qmheat(:,:,is,3)=qmheat(:,:,is,3) * spec(is)%dens*spec(is)%temp
                call get_volume_average (qmheat(:,:,is,3), mheat_perp(is))
                
                pmflux(:,:,is)=pmflux(:,:,is) * spec(is)%dens
                call get_volume_average (pmflux(:,:,is), mpart_fluxes(is))

                vmflux(:,:,is)=vmflux(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%stm
                call get_volume_average (vmflux(:,:,is), mmom_fluxes(is))
             end do
          end if
          if (fbpar > epsilon(0.0)) then
             do is = 1, nspec
                qbheat(:,:,is,1)=qbheat(:,:,is,1) * spec(is)%dens*spec(is)%temp
                call get_volume_average (qbheat(:,:,is,1), bheat_fluxes(is))

                qbheat(:,:,is,2)=qbheat(:,:,is,2) * spec(is)%dens*spec(is)%temp
                call get_volume_average (qbheat(:,:,is,2), bheat_par(is))

                qbheat(:,:,is,3)=qbheat(:,:,is,3) * spec(is)%dens*spec(is)%temp
                call get_volume_average (qbheat(:,:,is,3), bheat_perp(is))
                
                pbflux(:,:,is)=pbflux(:,:,is) * spec(is)%dens
                call get_volume_average (pbflux(:,:,is), bpart_fluxes(is))

                vbflux(:,:,is)=vbflux(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%stm
                call get_volume_average (vbflux(:,:,is), bmom_fluxes(is))
             end do
          end if
          pflux_avg = pflux_avg + (part_fluxes + mpart_fluxes + bpart_fluxes)*(t-t_old)
          qflux_avg = qflux_avg + (heat_fluxes + mheat_fluxes + bheat_fluxes)*(t-t_old)
          vflux_avg = vflux_avg + (mom_fluxes + mmom_fluxes + bmom_fluxes)*(t-t_old)
          if (write_hrate) heat_avg = heat_avg + h%imp_colls*(t-t_old)
!          t_old = t
       end if
    end if

    call broadcast (pflux_avg)
    call broadcast (qflux_avg)
    call broadcast (vflux_avg)
    if (write_hrate) call broadcast (heat_avg)

    fluxfac = 0.5
    fluxfac(1) = 1.0

    if (proc0) then
       if (print_flux_line) then
          write (unit=*, fmt="('t= ',e16.10,' <phi**2>= ',e12.6, &
               & ' heat fluxes: ', 5(1x,e12.6))") &
               t, phi2, heat_fluxes(1:min(nspec,5))
          write (unit=*, fmt="('t= ',e16.10,' <phi**2>= ',e12.6, &
               & ' energy exchange: ', 5(1x,e12.6))") &
               t, phi2, energy_exchange(1:min(nspec,5))
          if (fapar > epsilon(0.0)) then
             write (unit=*, fmt="('t= ',e16.10,' <apar**2>= ',e10.4, &
                  & ' heat flux m: ', 5(1x,e10.4))") &
                  t, apar2, mheat_fluxes(1:min(nspec,5))
          end if
          if (fbpar > epsilon(0.0)) then
             write (unit=*, fmt="('t= ',e16.10,' <bpar**2>= ',e10.4, &
                  & ' heat flux b: ', 5(1x,e10.4))") &
                  t, bpar2, bheat_fluxes(1:min(nspec,5))
          end if
       end if
       if (print_line) then
          do ik = 1, naky
             do it = 1, ntheta0
!                   write (unit=*, fmt="('aky=',f5.2, ' th0=',f7.2, &
!                          &' om=',2f8.3,' omav=', 2f8.3,' phtot=',e10.4)") &
!                        aky(ik), theta0(it,ik), &
!                        real( omega(it,ik)*woutunits(ik)), &
!                        aimag(omega(it,ik)*woutunits(ik)), &
!                        real( omegaavg(it,ik)*woutunits(ik)), &
!                        aimag(omegaavg(it,ik)*woutunits(ik)), &
!                        phitot(it,ik)
                write (unit=*, fmt="('ky=',f7.4, ' kx=',f7.4, &
                     &' om=',2f8.3,' omav=', 2f8.3,' phtot=',e8.2)") &
                     aky(ik), akx(it), &
                     real( omega(it,ik)*woutunits(ik)), &
                     aimag(omega(it,ik)*woutunits(ik)), &
                     real( omegaavg(it,ik)*woutunits(ik)), &
                     aimag(omegaavg(it,ik)*woutunits(ik)), &
                     phitot(it,ik)
             end do
          end do
          write (*,*) 
       end if
    end if

    i=istep/nwrite

    if (write_fields) call nc_write_fields (nout, phinew, aparnew, bparnew)  !MR
    if (write_moments) then !CMR
       call getmoms (ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)
       call nc_write_moments(nout, ntot)
    endif

    if (write_cross_phase .and. has_electron_species(spec)) then
       call get_cross_phase (phase_tot, phase_theta)
       if (proc0) write (unit=phase_unit, fmt="('t= ',e16.10,' phase_tot= ',e10.4,' phase_theta= ',e10.4)") &
            & t, phase_tot, phase_theta
    end if

    call prof_leaving ("loop_diagnostics-1")

    if (.not. (write_any .or. dump_any)) return

if (debug) write(6,*) "loop_diagnostics: -2"
    call prof_entering ("loop_diagnostics-2")

    call kperp (ntg_out, akperp)

    if (proc0 .and. write_any) then
       if (write_ascii) write (unit=out_unit, fmt=*) 'time=', t
       if (write_ascii .and. write_hrate) then
!
! For case with two species:
!
! Column     Item               
!   1        time              
!   2        Energy              
!   3        dEnergy/dt            
!   4        J_ant.E             
!   5        [h_(i+1)*h_*]/2 * C[h_(i+1)] * T_0 for species 1
!   6        [h_(i+1)*h_*]/2 * C[h_(i+1)] * T_0 for species 2
!   7       -[h H(h) * T_0]_1
!   8       -[h H(h) * T_0]_2
!   9       -[h C(h) * T_0]_1 
!  10       -[h C(h) * T_0]_2
!  11        [h w_* h]_1
!  12        [h w_* h]_2
!  13        [h * (q dchi/dt - dh/dt * T0)]_1
!  14        [h * (q dchi/dt - dh/dt * T0)]_2
!  15      sum (h C(h) * T_0)  in total, as in 5, 6      
!  16     -sum (h H(h) * T_0)      
!  17     -sum (h C(h) * T_0)   
!  18      sum (h w_* h)  
!  19      sum [h (q dchi/dt - dh/dt * T0)]
!  20      3 + 4 + 18 + 19
!  21      (k_perp A)**2
!  22      B_par**2
!  23      df_1 ** 2
!  24      df_2 ** 2
!  25      h_1 ** 2
!  26      h_2 ** 2
!  27      Phi_bar_1 ** 2
!  28      Phi_bar_2 ** 2
!
!
! For case with one species:
!
! Column     Item               
!   1        time              
!   2        Energy              
!   3        dEnergy/dt            
!   4        J_ant.E             
!   5        [h_(i+1)*h_*]/2 * C[h_(i+1)] * T_0 
!   6       -[h H(h) * T_0]
!   7       -[h C(h) * T_0]
!   8        [h w_* h]
!   9        [h * (q dchi/dt - dh/dt * T0)]_1
!  10      sum (h C(h) * T_0)  in total, as in 5, 6      
!  11     -sum (h H(h) * T_0)      
!  12     -sum (h C(h) * T_0)   
!  13      sum (h w_* h)  
!  14      sum [h (q dchi/dt - dh/dt * T0)]
!  15      3 + 4 + 9 + 10
!  16      (k_perp A)**2
!  17      B_par**2
!  18      df ** 2
!  19      h ** 2
!  20      Phi_bar ** 2

          write (unit=heat_unit, fmt="(28es12.4)") t,h % energy,  &
               h % energy_dot, h % antenna, h % imp_colls, h % hypercoll, h % collisions, &
               h % gradients, h % heating, sum(h % imp_colls), sum(h % hypercoll), sum(h % collisions), &
               sum(h % gradients), sum(h % heating),sum(h%heating)+h%antenna+sum(h%gradients)+h%energy_dot, &
               h % eapar, h % ebpar, h % delfs2(:),  h % hs2(:), h % phis2(:)
          
          do is=1,nspec
             write (unit=heat_unit2, fmt="(15es12.4)") t,h % energy,  &
                  h % energy_dot, h % antenna, h % imp_colls(is), h % hypercoll(is), h % collisions(is), &
                  h % gradients(is), h % heating(is), &
                  h % eapar, h % ebpar, h % delfs2(is),  h % hs2(is), h % phis2(is), real(is)
             write (unit=heat_unit2, fmt=*)
          end do
          write (unit=heat_unit2, fmt=*)

          call flush_output_file (heat_unit, ".heat")
          call flush_output_file (heat_unit2, ".heat2")

!GGH          write (unit=heat_unit, fmt="('t= ',e12.6,' energy= ',e12.6)") t, h % energy
!GGH          write (unit=heat_unit, fmt="('t= ',e12.6,' energy_dot= ',e12.6)") t, h % energy_dot
!GGH          write (unit=heat_unit, fmt="('t= ',e12.6,' J_ant.E= ',e12.6)") t, h % antenna

!GGH          write (unit=heat_unit, fmt="('t= ',e12.6,' hyperC= ',12(1x,e12.6))") t, h % hypercoll
!GGH          write (unit=heat_unit, fmt="('t= ',e12.6,' hCh= ',12(1x,e12.6))") t, h % collisions
!GGH          write (unit=heat_unit, fmt="('t= ',e12.6,' hw*= ',12(1x,e12.6))") t, h % gradients
!GGH!         write (unit=heat_unit, fmt="('t= ',e12.6,' hwd= ',12(1x,e12.6))") t, h % curvature
!GGH          write (unit=heat_unit, fmt="('t= ',e12.6,' heating= ',12(1x,e12.6))") t, h % heating

!GGH          write (unit=heat_unit, fmt="('t= ',e12.6,' total_hvisc= ',e12.6)") t, sum(h % hypervisc)
!GGH          write (unit=heat_unit, fmt="('t= ',e12.6,' total_hyperC= ',e12.6)") t, sum(h % hypercoll)
!GGH          write (unit=heat_unit, fmt="('t= ',e12.6,' total_hCh= ',e12.6)") t, sum(h % collisions)
!GGH          write (unit=heat_unit, fmt="('t= ',e12.6,' total_hw*= ',e12.6)") t, sum(h % gradients)
!GGH          write (unit=heat_unit, fmt="('t= ',e12.6,' total_heating= ',e12.6)") t, sum(h % heating)

!GGH          write (unit=heat_unit, fmt="('t= ',e12.6,' total_power= ',e12.6)") t, &
!GGH               sum(h%heating)+h%antenna+sum(h%gradients)+h%energy_dot
          !GGH TEST try adding sqrt(2.) to the edot
!GGH          write (unit=heat_unit, fmt="('t= ',e12.6,' total_power= ',e12.6)") t, &
!GGH               sum(h%heating)+h%antenna+sum(h%gradients)+h%energy_dot*sqrt(2.)
!GGH          write (unit=heat_unit, fmt='(a)') ''
       end if

       !Write out data for j_external
       if (write_ascii .and. write_jext) then
          do ik=1,naky
             do it = 1, ntheta0
                if (j_ext(ik,it) .ne. 0.) then
                   write (unit=jext_unit, fmt="(es12.4,i4,i4,es12.4)")  &
                        t,it,ik,j_ext(ik,it)
                endif
             enddo
          enddo
       end if
!>GGH

       if (write_flux_line) then
          hflux_tot = 0.
          zflux_tot = 0.
          if (fphi > epsilon(0.0)) then
             if (write_ascii) then
                write (unit=out_unit, fmt="('t= ',e16.10,' <phi**2>= ',e10.4, &
                     & ' heat fluxes: ', 5(1x,e10.4),' qflux_avg: ', 5(1x,e10.4))") &
!                     t, phi2, heat_fluxes(1:min(nspec,5)), qflux_avg(1:min(nspec,5))/t
                     t, phi2, heat_fluxes(1:min(nspec,5)), energy_exchange(1:min(nspec,5))
                write (unit=out_unit, fmt="('t= ',e16.10,' <phi**2>= ',e10.4, &
                     & ' part fluxes: ', 5(1x,e10.4),' pflux_avg: ', 5(1x,e10.4))") &
                     t, phi2, part_fluxes(1:min(nspec,5)), pflux_avg(1:min(nspec,5))/t
                write (unit=out_unit, fmt="('t= ',e16.10,' <phi**2>= ',e10.4, &
                     & ' mom fluxes: ', 5(1x,e10.4),' vflux_avg: ', 5(1x,e10.4))") &
                     t, phi2, mom_fluxes(1:min(nspec,5)), vflux_avg(1:min(nspec,5))/t

#ifdef LOWFLOW
                   write (unit=out_unit, fmt="('t= ',e16.10,' <phi**2>= ',e10.4, &
                        & ' lfmom fluxes: ', 5(1x,e10.4),' lfvflx1: ', 5(1x,e10.4))") &
                        t, phi2, lfmom_fluxes(1:min(nspec,5)), vflux1_avg(1:min(nspec,5))
#endif
             end if
             hflux_tot = sum(heat_fluxes)
             vflux_tot = sum(mom_fluxes)
             zflux_tot = sum(part_fluxes*spec%z)
          end if
          if (fapar > epsilon(0.0)) then
             if (write_lorentzian .and. write_ascii) then
                wtmp_new = antenna_w()
                if (real(wtmp_old) /= 0. .and. wtmp_new /= wtmp_old) &
                     write (unit=out_unit, fmt="('w= ',e16.10, &
                     &  ' amp= ',e16.10)") real(wtmp_new), sqrt(2.*apar2)
                wtmp_old = wtmp_new                
             end if
             if (write_ascii) then
                write (unit=out_unit, fmt="('t= ',e16.10,' <apar**2>= ',e10.4, &
                     & ' heat mluxes: ', 5(1x,e10.4))") &
                     t, apar2, mheat_fluxes(1:min(nspec,5))
                write (unit=out_unit, fmt="('t= ',e16.10,' <apar**2>= ',e10.4, &
                     & ' part mluxes: ', 5(1x,e10.4))") &
                     t, apar2, mpart_fluxes(1:min(nspec,5))
             end if
             hflux_tot = hflux_tot + sum(mheat_fluxes)
             vflux_tot = vflux_tot + sum(mmom_fluxes)
             zflux_tot = zflux_tot + sum(mpart_fluxes*spec%z)
          end if
          if (fbpar > epsilon(0.0)) then
             if (write_ascii) then
                write (unit=out_unit, fmt="('t= ',e16.10,' <bpar**2>= ',e10.4, &
                     & ' heat bluxes: ', 5(1x,e10.4))") &
                     t, bpar2, bheat_fluxes(1:min(nspec,5))
                write (unit=out_unit, fmt="('t= ',e16.10,' <bpar**2>= ',e10.4, &
                     & ' part bluxes: ', 5(1x,e10.4))") &
                     t, bpar2, bpart_fluxes(1:min(nspec,5))
             end if
             hflux_tot = hflux_tot + sum(bheat_fluxes)
             vflux_tot = vflux_tot + sum(bmom_fluxes)
             zflux_tot = zflux_tot + sum(bpart_fluxes*spec%z)
          end if
          if (write_ascii) write (unit=out_unit, fmt="('t= ',e16.10,' h_tot= ',e10.4, &
               & ' z_tot= ',e10.4)") t, hflux_tot, zflux_tot
          if (write_nl_flux) then
             call nc_qflux (nout, qheat(:,:,:,1), qmheat(:,:,:,1), qbheat(:,:,:,1), &
                  heat_par, mheat_par, bheat_par, &
                  heat_perp, mheat_perp, bheat_perp, &
                  heat_fluxes, mheat_fluxes, bheat_fluxes, x_qmflux, hflux_tot, &
                  energy_exchange)
                  ! Update the target array in parameter_scan_arrays
! below line gives out-of-bounds array for runs inside trinity
!                  scan_hflux(nout) = hflux_tot
                  scan_hflux(mod(nout-1,nstep/nwrite+1)+1) = hflux_tot
!          call nc_qflux (nout, qheat, qmheat, qbheat, &
!               heat_par, mheat_par, bheat_par, &
!               heat_perp, mheat_perp, bheat_perp, &
!               heat_fluxes, mheat_fluxes, bheat_fluxes, hflux_tot)
             call nc_vflux (nout, vflux, vmflux, vbflux, &
                  mom_fluxes, mmom_fluxes, bmom_fluxes, vflux_tot, &
                  vflux_par, vflux_perp, vflux0, vflux1)
                  ! Update the target array in parameter_scan_arrays
! below line gives out-of-bounds array for runs inside trinity
!                  scan_momflux(nout) = vflux_tot
                  scan_momflux(mod(nout-1,nstep/nwrite+1)+1) = vflux_tot
             call nc_pflux (nout, pflux, pmflux, pbflux, &
                  part_fluxes, mpart_fluxes, bpart_fluxes, zflux_tot)
             if (write_flux_emu) call nc_flux_emu(nout,pflux_emu,mflux_emu,bflux_emu)
             if (write_flux_e) call nc_flux_e(nout,pflux_e,mflux_e,bflux_e,trapped_flux,passing_flux)
             if (write_eflux) call nc_eflux(nout,peflux,meflux,beflux)
          end if
          call nc_loop (nout, t, fluxfac, &
               phinew(igomega,:,:), phi2, phi2_by_mode, &
               aparnew(igomega,:,:), apar2, apar2_by_mode, &
               bparnew(igomega,:,:), bpar2, bpar2_by_mode, &
               h, hk, he, omega, omegaavg, woutunits, phitot, write_omega, write_hrate, write_hrate_e)
! below line gives out-of-bounds array for runs inside trinity
!               scan_phi2_tot(nout) = phi2
               scan_phi2_tot(mod(nout-1,nstep/nwrite+1)+1) = phi2
       end if
       if (write_ascii) then
          do ik = 1, naky
             do it = 1, ntheta0
!                write (out_unit,*) "aky=",aky(ik)," theta0=",theta0(it,ik)

                if (write_line) then
                   write (out_unit, "('t= ',e16.10,' aky= ',1pe12.4, ' akx= ',1pe12.4, &
                        &' om= ',1p2e12.4,' omav= ', 1p2e12.4,' phtot= ',1pe12.4)") &
                        t, aky(ik), akx(it), &
                        real( omega(it,ik)*woutunits(ik)), &
                        aimag(omega(it,ik)*woutunits(ik)), &
                        real( omegaavg(it,ik)*woutunits(ik)), &
                        aimag(omegaavg(it,ik)*woutunits(ik)), &
                        phitot(it,ik)
                end if
                
!                if (write_omega) write (out_unit, *) &
!                     'omega=', omega(it,ik), &
!                     ' omega/(vt/a)=', real(omega(it,ik)*woutunits(ik)), &
!                     aimag(omega(it,ik)*woutunits(ik))
                if (write_omega) write (out_unit,&
                         fmt='(" omega= (",1pe12.4,",",1pe12.4,")",t45,"omega/(vt/a)= (",1pe12.4,",",1pe12.4,")")') &
                     omega(it,ik)/tunits(ik), omega(it,ik)*woutunits(ik)
                if (write_omavg) write (out_unit,&
                         fmt='(" omavg= (",1pe12.4,",",1pe12.4,")",t45,"omavg/(vt/a)= (",1pe12.4,",",1pe12.4,")")') &
                         omegaavg(it,ik)/tunits(ik), omegaavg(it,ik)*woutunits(ik)
!                if (write_omavg) write (out_unit, *) &
!                     'omavg=', omegaavg(it,ik), &
!                     ' omavg/(vt/a)=', real(omegaavg(it,ik)*woutunits(ik)), &
!                     aimag(omegaavg(it,ik)*woutunits(ik))
                
             end do
          end do
       end if
       do ik = 1, naky
          do it = 1, ntheta0

             if (write_nl_flux) then
                if (write_ascii) then
!                   write (out_unit,"('ik,it,aky,akx,<phi**2>,t: ', &
!                        & 2i5,4(1x,e12.6))") &
!                        ik, it, aky(ik), akx(it), phi2_by_mode(it,ik), t
!                   write (out_unit,"('ik,it,aky,akx,<apar**2>,t: ', &
!                        & 2i5,4(1x,e12.6))") &
!                        ik, it, aky(ik), akx(it), apar2_by_mode(it,ik), t
                end if
!                if (ntheta0 > 1 .and. ik == 1) then
!                   write (out_unit, "('it,akx,<phi**2>,t: ',i5,3(1x,e12.6))") &
!                        it, akx(it), sum(phi2_by_mode(it,:)*fluxfac(:)), t
!                end if
!                if (naky > 1 .and. it == 1) then
!                   write (out_unit, "('ik,aky,<phi**2>,t: ',i5,3(1x,e12.6))") &
!                        ik, aky(ik), sum(phi2_by_mode(:,ik))*fluxfac(ik), t
!                end if
             end if
          end do
       end do
    end if

    call broadcast(scan_hflux)
    call broadcast(scan_momflux)
    call broadcast(scan_phi2_tot)

    if (write_cerr) call collision_error(phinew,bparnew,last)

    call broadcast (write_symmetry)
    if (write_symmetry) then
       allocate (vflx_sym(-ntgrid:ntgrid,nlambda*negrid,nspec))
       call flux_vs_theta_vs_vpa (vflx_sym)
       if (proc0) call nc_loop_sym (nout, vflx_sym)
       deallocate (vflx_sym)
    end if

    call broadcast (write_correlation)
    if (write_correlation) then
       allocate (phi_corr_2pi(-ntgrid:ntgrid,naky))
       call correlation (phi_corr_2pi)
       if (proc0) call nc_loop_corr (nout,phi_corr_2pi)
       deallocate (phi_corr_2pi)
    end if

    call broadcast (write_correlation_extend)
    if (write_correlation_extend .and. istep > nstep/4) then
       if (.not. allocated(phicorr_sum)) then
          ntg_extend = (2*ntgrid+1)*((ntheta0-1)/jtwist_out+1)
          nth0_extend = min(ntheta0,jtwist_out*(naky-1))
          allocate (phicorr_sum(ntg_extend,ntheta0,naky)) ; phicorr_sum = 0.0
          allocate (phiextend_sum(ntg_extend,ntheta0,naky)) ; phiextend_sum = 0.0
          tcorr0 = t_old
       end if
       allocate (phi_corr(ntg_extend,ntheta0,naky))
       allocate (phi2_extend(ntg_extend,ntheta0,naky))
       call correlation_extend (phi_corr, phi2_extend)
       phicorr_sum = phicorr_sum + phi_corr*(t-t_old)
       phiextend_sum = phiextend_sum + phi2_extend*(t-t_old)
       if (proc0 .and. mod(istep,nwrite_mult*nwrite)==0) then
          call nc_loop_corr_extend (nout_big, phicorr_sum/(t-tcorr0), phiextend_sum/(t-tcorr0))
          nout_big = nout_big + 1
       end if
       deallocate (phi_corr, phi2_extend)
    end if

    call broadcast (write_parity)
    if (write_parity) then

       ! initialize layouts for parity diagnostic
       if (first) then
          call init_parity_layouts (naky, nlambda, negrid, nspec)
          first = .false.
       end if

       allocate (gparity(-ntgrid:ntgrid,ntheta0,2,p_lo%llim_proc:p_lo%ulim_alloc))
       allocate (g0(-ntgrid:ntgrid,2,p_lo%llim_proc:p_lo%ulim_alloc))
       allocate (gm(-ntgrid:ntgrid,2,p_lo%llim_proc:p_lo%ulim_alloc))
       allocate (gp(-ntgrid:ntgrid,2,p_lo%llim_proc:p_lo%ulim_alloc))
       allocate (gmnorm(-ntgrid:ntgrid,2,p_lo%llim_proc:p_lo%ulim_alloc))
       allocate (gmint(-ntgrid:ntgrid,naky,nspec))
       allocate (gpint(-ntgrid:ntgrid,naky,nspec))
       allocate (gmnormint(-ntgrid:ntgrid,naky,nspec))
       allocate (gmavg(ntheta0,naky,nspec))
       allocate (gpavg(ntheta0,naky,nspec))
       allocate (gmnormavg(ntheta0,naky,nspec))
       allocate (gmtot(nspec), gm_nokx_tot(nspec), gm_nosig_tot(nspec))
       allocate (gptot(nspec), gp_nokx_tot(nspec), gp_nosig_tot(nspec))
       allocate (gtot(nspec), gmnormtot(nspec), gmnorm_nokx_tot(nspec), gmnorm_nosig_tot(nspec))
       allocate (phim(-ntgrid:ntgrid,naky))
       allocate (g_kint(-ntgrid:ntgrid,2,nspec), gm_kint(-ntgrid:ntgrid,2,nspec), gp_kint(-ntgrid:ntgrid,2,nspec))
       allocate (g_avg(ntheta0,nspec), gnorm_avg(ntheta0,nspec))
       allocate (g_all_tot(nspec), g_nokx_tot(nspec), g_nosig_tot(nspec), gtmp(nspec))
       allocate (gnorm_all_tot(nspec), gnorm_nokx_tot(nspec), gnorm_nosig_tot(nspec))

       ! convert from g to h
       call g_adjust (gnew, phinew, bparnew, fphi, fbpar)

       ! below we define gparity = J0 * h, where delta f = h - (q phi / T) F0
       ! because we're ultimately interested in int J0 h * phi (i.e. particle flux)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             do iglo = g_lo%llim_world, g_lo%ulim_world
                ik = ik_idx(g_lo,iglo)
                ie = ie_idx(g_lo,iglo)
                il = il_idx(g_lo,iglo)
                is = is_idx(g_lo,iglo)
                it = it_idx(g_lo,iglo)
                ! count total index for given (ik,il,ie,is) combo (dependent on layout)
                iplo = idx(p_lo,ik,il,ie,is)
                ! check to see if value of g corresponding to iglo is on this processor
                if (idx_local(g_lo,iglo)) then
                   if (idx_local(p_lo,iplo)) then
                      ! if g value corresponding to iplo should be on this processor, then get it
                      gparity(ig,it,isgn,iplo) = gnew(ig,isgn,iglo)*aj0(ig,iglo)
                   else
                      ! otherwise, send g for this iplo value to correct processor
                      call send (gnew(ig,isgn,iglo)*aj0(ig,iglo),proc_id(p_lo,iplo))
                   end if
                else if (idx_local(p_lo,iplo)) then
                   ! if g for this iplo not on processor, receive it
                   call receive (gparity(ig,it,isgn,iplo),proc_id(g_lo,iglo))
                end if
             end do
          end do
       end do

       ! convert from h back to g
       call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar)
       
       ! first diagnostic is phase space average of 
       ! |J0*(h(z,vpa,kx) +/- h(-z,-vpa,-kx))|**2 / ( |J0*h(z,vpa,kx)|**2 + |J0*h(-z,-vpa,-kx)|**2 )
       do it = 1, ntheta0
          ! have to treat kx=0 specially
          if (it == 1) then
             do isgn = 1, 2
                sgn2 = mod(isgn,2)+1
                do ig = -ntgrid, ntgrid
                   g0(ig,isgn,:) = gparity(-ig,1,sgn2,:)
                end do
             end do
          else
             do isgn = 1, 2
                sgn2 = mod(isgn,2)+1
                do ig = -ntgrid, ntgrid
                   g0(ig,isgn,:) = gparity(-ig,ntheta0-it+2,sgn2,:)
                end do
             end do
          end if
          gm = gparity(:,it,:,:)-g0
          gp = gparity(:,it,:,:)+g0
          gmnorm = real(g0*conjg(g0))
          ! integrate out velocity dependence
          call integrate_moment (real(gm*conjg(gm)),gmint,1)
          call integrate_moment (real(gp*conjg(gp)),gpint,1)
          call integrate_moment (gmnorm,gmnormint,1)
          ! average along field line
          do is = 1, nspec
             do ik = 1, naky
                call get_fldline_avg (gmint(:,ik,is),gmavg(it,ik,is))
                call get_fldline_avg (gpint(:,ik,is),gpavg(it,ik,is))
                call get_fldline_avg (gmnormint(:,ik,is),gmnormavg(it,ik,is))
             end do
          end do

          ! phim(theta,kx) = phi(-theta,-kx)
          ! have to treat kx=0 specially
          if (it == 1) then
             do ig = -ntgrid, ntgrid
                phim(ig,:) = phinew(-ig,1,:)
             end do
          else
             do ig = -ntgrid, ntgrid
                phim(ig,:) = phinew(-ig,ntheta0-it+2,:)
             end do
          end if

          ! want int dtheta sum_{kx} int d3v | sum_{ky} [ J0*(h+ * conjg(phi-) + h- * conjg(phi+)) ky ] |**2
          ! J0*(h+ * conjg(phi-) + h- * conjg(phi+)) = h*conjg(phi) - h(-theta,-vpa,-kx)*conjg(phi(-theta,-kx))
          do iplo = p_lo%llim_proc, p_lo%ulim_proc
             ik = ik_idx(p_lo,iplo)
             do isgn = 1, 2
                gm(:,isgn,iplo) = gparity(:,it,isgn,iplo)*conjg(phinew(:,it,ik)) - g0(:,isgn,iplo)*conjg(phim(:,ik))
             end do
          end do
          do isgn = 1, 2
             do ig = -ntgrid, ntgrid
                call integrate_kysum (gm(ig,isgn,p_lo%llim_proc:p_lo%ulim_proc),ig,g_kint(ig,isgn,:),1)
             end do
          end do
          do is = 1, nspec
             call get_fldline_avg (g_kint(:,1,is)+g_kint(:,2,is),g_avg(it,is))
          end do

          ! get normalizing term for above diagnostic
          do iplo = p_lo%llim_proc, p_lo%ulim_proc
             ik = ik_idx(p_lo,iplo)
             do isgn = 1, 2
                gm(:,isgn,iplo) = gparity(:,it,isgn,iplo)*conjg(phinew(:,it,ik))
                gp(:,isgn,iplo) = g0(:,isgn,iplo)*conjg(phim(:,ik))
             end do
          end do
          do isgn = 1, 2
             do ig = -ntgrid, ntgrid
                call integrate_kysum (gm(ig,isgn,:),ig,gm_kint(ig,isgn,:),1)
                call integrate_kysum (gp(ig,isgn,:),ig,gp_kint(ig,isgn,:),1)
             end do
          end do
          do is = 1, nspec
             call get_fldline_avg (gm_kint(:,1,is)+gm_kint(:,2,is) &
                  + gp_kint(:,1,is)+gp_kint(:,2,is),gnorm_avg(it,is))
          end do
       end do

       ! average over perp plane
       do is = 1, nspec
          call get_volume_average (gmavg(:,:,is), gmtot(is))
          call get_volume_average (gpavg(:,:,is), gptot(is))
          call get_volume_average (gmnormavg(:,:,is), gmnormtot(is))    
          g_all_tot(is) = sum(g_avg(:,is))
          gnorm_all_tot(is) = sum(gnorm_avg(:,is))
       end do

       allocate (gmx(-ntgrid:ntgrid,ntheta0,2,p_lo%llim_proc:p_lo%ulim_alloc))
       allocate (gpx(-ntgrid:ntgrid,ntheta0,2,p_lo%llim_proc:p_lo%ulim_alloc))

       ! now we want diagnostic of phase space average of
       ! |J0*(h(z,vpa) +/- h(-z,-vpa))|**2 / ( |J0*h(z,vpa)|**2 + |J0*h(-z,-vpa)|**2 )
       do it = 1, ntheta0
          do isgn = 1, 2
             sgn2 = mod(isgn,2)+1
             do ig = -ntgrid, ntgrid
                g0(ig,isgn,:) = gparity(-ig,it,sgn2,:)
             end do
          end do
          gm = gparity(:,it,:,:)-g0
          gp = gparity(:,it,:,:)+g0
          gmnorm = real(g0*conjg(g0))
          ! integrate out velocity dependence
          call integrate_moment (real(gm*conjg(gm)),gmint,1)
          call integrate_moment (real(gp*conjg(gp)),gpint,1)
          call integrate_moment (gmnorm,gmnormint,1)
          ! average along field line
          do is = 1, nspec
             do ik = 1, naky
                call get_fldline_avg (gmint(:,ik,is),gmavg(it,ik,is))
                call get_fldline_avg (gpint(:,ik,is),gpavg(it,ik,is))
                call get_fldline_avg (gmnormint(:,ik,is),gmnormavg(it,ik,is))
             end do
          end do

          ! phim(theta) = phi(-theta)
          ! have to treat kx=0 specially
          do ig = -ntgrid, ntgrid
             phim(ig,:) = phinew(-ig,it,:)
          end do
          
          ! want int dtheta int d3v | sum_{kx} sum_{ky} [ J0*(h+ * conjg(phi-) + h- * conjg(phi+)) ky ] |**2
          ! J0*(h+ * conjg(phi-) + h- * conjg(phi+)) = h*conjg(phi) - h(-theta,-vpa)*conjg(phi(-theta))
          do iplo = p_lo%llim_proc, p_lo%ulim_proc
             ik = ik_idx(p_lo,iplo)
             do isgn = 1, 2
                gmx(:,it,isgn,iplo) = g0(:,isgn,iplo)*conjg(phim(:,ik))
                gpx(:,it,isgn,iplo) = gparity(:,it,isgn,iplo)*conjg(phinew(:,it,ik)) - gmx(:,it,isgn,iplo) 
             end do
          end do
       end do

       ! sum over kx
       gp = sum(gpx,2)
       deallocate (gpx)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             call integrate_kysum (gp(ig,isgn,:),ig,g_kint(ig,isgn,:),1)
          end do
       end do
       do is = 1, nspec
          call get_fldline_avg (g_kint(:,1,is)+g_kint(:,2,is), g_nokx_tot(is))
       end do

       ! get normalizing terms for above diagnostic
       gm = sum(gmx,2)
       deallocate (gmx)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             call integrate_kysum (gm(ig,isgn,:),ig,gm_kint(ig,isgn,:),1)
             call integrate_kysum (gm(ig,isgn,:)+gp(ig,isgn,:),ig,gp_kint(ig,isgn,:),1)
          end do
       end do
       do is = 1, nspec
          call get_fldline_avg (gm_kint(:,1,is)+gm_kint(:,2,is) &
               + gp_kint(:,1,is)+gp_kint(:,2,is), gnorm_nokx_tot(is))
       end do

       ! average over perp plane
       do is = 1, nspec
          call get_volume_average (gmavg(:,:,is), gm_nokx_tot(is))
          call get_volume_average (gpavg(:,:,is), gp_nokx_tot(is))
          call get_volume_average (gmnormavg(:,:,is), gmnorm_nokx_tot(is))    
       end do

       ! final diagnostic is phase space average of 
       ! |J0*(h(z,kx) +/- h(-z,-kx))|**2 / ( |J0*h(z,kx)|**2 + |J0*h(-z,-kx)|**2 )
       do it = 1, ntheta0
          ! have to treat kx=0 specially
          if (it == 1) then
             do ig = -ntgrid, ntgrid
                g0(ig,:,:) = gparity(-ig,1,:,:)
             end do
          else
             do ig = -ntgrid, ntgrid
                g0(ig,:,:) = gparity(-ig,ntheta0-it+2,:,:)
             end do
          end if
          gm = gparity(:,it,:,:)-g0
          gp = gparity(:,it,:,:)+g0
          gmnorm = real(g0*conjg(g0))
          ! integrate out velocity dependence
          call integrate_moment (real(gm*conjg(gm)),gmint,1)
          call integrate_moment (real(gp*conjg(gp)),gpint,1)
          call integrate_moment (gmnorm,gmnormint,1)
          ! average along field line
          do is = 1, nspec
             do ik = 1, naky
                call get_fldline_avg (gmint(:,ik,is),gmavg(it,ik,is))
                call get_fldline_avg (gpint(:,ik,is),gpavg(it,ik,is))
                call get_fldline_avg (gmnormint(:,ik,is),gmnormavg(it,ik,is))
             end do
          end do

          ! phim(theta,kx) = phi(-theta,-kx)
          ! have to treat kx=0 specially
          if (it == 1) then
             do ig = -ntgrid, ntgrid
                phim(ig,:) = phinew(-ig,1,:)
             end do
          else
             do ig = -ntgrid, ntgrid
                phim(ig,:) = phinew(-ig,ntheta0-it+2,:)
             end do
          end if

          ! want int dtheta sum_{kx} int dE dmu | sum_{sigma} sum_{ky} [ J0*(h+ * conjg(phi-) + h- * conjg(phi+)) ky ] |**2
          ! J0*(h+ * conjg(phi-) + h- * conjg(phi+)) = h*conjg(phi) - h(-theta,-kx)*conjg(phi(-theta,-kx))
          do iplo = p_lo%llim_proc, p_lo%ulim_proc
             ik = ik_idx(p_lo,iplo)
             do isgn = 1, 2
                gm(:,isgn,iplo) = gparity(:,it,isgn,iplo)*conjg(phinew(:,it,ik)) - g0(:,isgn,iplo)*conjg(phim(:,ik))
             end do
          end do
          do ig = -ntgrid, ntgrid
             call integrate_kysum (gm(ig,1,:)+gm(ig,2,:),ig,g_kint(ig,1,:),1)
          end do
          do is = 1, nspec
             call get_fldline_avg (g_kint(:,1,is),g_avg(it,is))
          end do

          ! get normalizing term for above diagnostic
          do iplo = p_lo%llim_proc, p_lo%ulim_proc
             ik = ik_idx(p_lo,iplo)
             do isgn = 1, 2
                gm(:,isgn,iplo) = gparity(:,it,isgn,iplo)*conjg(phinew(:,it,ik))
                gp(:,isgn,iplo) = g0(:,isgn,iplo)*conjg(phim(:,ik))
             end do
          end do
          do ig = -ntgrid, ntgrid
             call integrate_kysum (gm(ig,1,:)+gm(ig,2,:),ig,gm_kint(ig,1,:),1)
             call integrate_kysum (gp(ig,1,:)+gp(ig,2,:),ig,gp_kint(ig,1,:),1)
          end do
          do is = 1, nspec
             call get_fldline_avg (gm_kint(:,1,is)+gp_kint(:,1,is),gnorm_avg(it,is))
          end do
       end do

       ! average over perp plane
       do is = 1, nspec
          call get_volume_average (gmavg(:,:,is), gm_nosig_tot(is))
          call get_volume_average (gpavg(:,:,is), gp_nosig_tot(is))
          call get_volume_average (gmnormavg(:,:,is), gmnorm_nosig_tot(is))    
          g_nosig_tot(is) = sum(g_avg(:,is))
          gnorm_nosig_tot(is) = sum(gnorm_avg(:,is))
       end do

       deallocate (gparity) ; allocate (gparity(-ntgrid:ntgrid,ntheta0,naky,nspec))
       ! obtain normalization factor = int over phase space of |g|**2
       call g_adjust (gnew, phinew, bparnew, fphi, fbpar)
       call integrate_moment (spread(aj0,2,2)*spread(aj0,2,2)*gnew*conjg(gnew), gparity, 1)
       call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar)
       do is = 1, nspec
          do ik = 1, naky
             do it = 1, ntheta0
                call get_fldline_avg (real(gparity(:,it,ik,is)),gmavg(it,ik,is))
             end do
          end do
          call get_volume_average (gmavg(:,:,is), gtot(is))
       end do

       ! normalize g(theta,vpa,kx) - g(-theta,-vpa,-kx) terms
       where (gtot+gmnormtot > epsilon(0.0))
          gmtot = sqrt(gmtot/(gtot+gmnormtot)) ; gptot = sqrt(gptot/(gtot+gmnormtot))
       elsewhere
          gmtot = sqrt(gmtot) ; gptot = sqrt(gptot)
       end where

       where (real(gnorm_all_tot) > epsilon(0.0))
          gtmp = sqrt(real(g_all_tot)/real(gnorm_all_tot))
       elsewhere
          gtmp = sqrt(real(g_all_tot))
       end where
       where (aimag(gnorm_all_tot) > epsilon(0.0))
          g_all_tot = gtmp + zi*sqrt(aimag(g_all_tot)/aimag(gnorm_all_tot))
       elsewhere
          g_all_tot = gtmp + zi*sqrt(aimag(g_all_tot))
       end where

       ! normalize g(theta,vpa) +/- g(-theta,-vpa) terms
       where (gtot+gmnorm_nokx_tot > epsilon(0.0))
          gm_nokx_tot = sqrt(gm_nokx_tot/(gtot+gmnorm_nokx_tot)) ; gp_nokx_tot = sqrt(gp_nokx_tot/(gtot+gmnorm_nokx_tot))
       elsewhere
          gm_nokx_tot = sqrt(gm_nokx_tot) ; gp_nokx_tot = sqrt(gp_nokx_tot)
       end where
       
       where (real(gnorm_nokx_tot) > epsilon(0.0))
          gtmp = sqrt(real(g_nokx_tot)/real(gnorm_nokx_tot))
       elsewhere
          gtmp = sqrt(real(g_nokx_tot))
       end where
       where (aimag(gnorm_nokx_tot) > epsilon(0.0))
          g_nokx_tot = gtmp + zi*sqrt(aimag(g_nokx_tot)/aimag(gnorm_nokx_tot))
       elsewhere
          g_nokx_tot = gtmp + zi*sqrt(aimag(g_nokx_tot))
       end where

       ! normalize g(theta,kx) +/ g(-theta,-kx) terms
       where (gtot+gmnorm_nosig_tot > epsilon(0.0))
          gm_nosig_tot = sqrt(gm_nosig_tot/(gtot+gmnorm_nosig_tot)) ; gp_nosig_tot = sqrt(gp_nosig_tot/(gtot+gmnorm_nosig_tot))
       elsewhere
          gm_nosig_tot = sqrt(gm_nosig_tot) ; gp_nosig_tot = sqrt(gp_nosig_tot)
       end where

       where (real(gnorm_nosig_tot) > epsilon(0.0))
          gtmp = sqrt(real(g_nosig_tot)/real(gnorm_nosig_tot))
       elsewhere
          gtmp = sqrt(real(g_nosig_tot))
       end where
       where (aimag(gnorm_nosig_tot) > epsilon(0.0))
          g_nosig_tot = gtmp + zi*sqrt(aimag(g_nosig_tot)/aimag(gnorm_nosig_tot))
       elsewhere
          g_nosig_tot = gtmp + zi*sqrt(aimag(g_nosig_tot))
       end where

       if (proc0) write (parity_unit,"(19(1x,e12.5))") t, gmtot, gptot, real(g_all_tot), aimag(g_all_tot), &
            real(gnorm_all_tot), aimag(gnorm_all_tot), gm_nokx_tot, gp_nokx_tot, real(g_nokx_tot), aimag(g_nokx_tot), &
            real(gnorm_nokx_tot), aimag(gnorm_nokx_tot), gm_nosig_tot, gp_nosig_tot, &
            real(g_nosig_tot), aimag(g_nosig_tot), real(gnorm_nosig_tot), aimag(gnorm_nosig_tot)

       deallocate (gparity, g0)
       deallocate (gm, gp, gmnorm, gmint, gpint, gmnormint)
       deallocate (gmavg, gpavg, gmnormavg)
       deallocate (gmtot, gm_nokx_tot, gm_nosig_tot)
       deallocate (gptot, gp_nokx_tot, gp_nosig_tot)
       deallocate (gtot, gmnormtot, gmnorm_nokx_tot, gmnorm_nosig_tot)
       deallocate (g_avg, gnorm_avg)
       deallocate (g_kint, gm_kint, gp_kint)
       deallocate (g_all_tot, g_nokx_tot, g_nosig_tot, gtmp)
       deallocate (gnorm_all_tot, gnorm_nokx_tot, gnorm_nosig_tot)
       deallocate (phim)
    end if

    call broadcast (write_avg_moments)
    if (write_avg_moments) then

       call getmoms (ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)

       if (proc0) then
          allocate (dl_over_b(-ntg_out:ntg_out))

          dl_over_b = delthet(-ntg_out:ntg_out)*jacob(-ntg_out:ntg_out)
          dl_over_b = dl_over_b / sum(dl_over_b)

          do is = 1, nspec
             do it = 1, ntheta0
                ntot00(it, is)   = sum( ntot(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                density00(it,is) = sum(density(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                upar00(it, is)   = sum( upar(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                tpar00(it, is)   = sum( tpar(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                tperp00(it, is)  = sum(tperp(-ntg_out:ntg_out,it,1,is)*dl_over_b)
             end do
          end do

          do it = 1, ntheta0
             phi00(it)   = sum( phinew(-ntg_out:ntg_out,it,1)*dl_over_b)
          end do

          deallocate (dl_over_b)

          do is = 1, nspec
             call get_vol_average (ntot(:,:,:,is), ntot(:,:,:,is), &
                  ntot2(is), ntot2_by_mode(:,:,is))
             call get_vol_average (tpar(:,:,:,is), tpar(:,:,:,is), &
                  tpar2(is), tpar2_by_mode(:,:,is))
             call get_vol_average (tperp(:,:,:,is), tperp(:,:,:,is), &
                  tperp2(is), tperp2_by_mode(:,:,is))
          end do

          do is = 1, nspec
             call get_vol_average (ntot(0,:,:,is), ntot(0,:,:,is), &
                  ntot20(is), ntot20_by_mode(:,:,is))
          end do

!          write (out_unit, "('t,u0 ',2(1x,e12.6))") t, aimag(upar00 (1,1))

          call nc_loop_moments (nout, ntot2, ntot2_by_mode, ntot20, &
               ntot20_by_mode, phi00, ntot00, density00, upar00, &
               tpar00, tperp00, tpar2_by_mode, tperp2_by_mode)

       end if
    end if

    ! RN> output not guiding center moments in x-y plane
    if (write_full_moments_notgc) then !RN
       call getmoms_notgc(density,upar,tpar,tperp,ntot)
       if(proc0) then
          call nc_loop_fullmom(nout,t, &
               & ntot(igomega,:,:,:),density(igomega,:,:,:), &
               & upar(igomega,:,:,:), &
               & tpar(igomega,:,:,:),tperp(igomega,:,:,:) )
       endif
    endif

    if (proc0 .and. dump_any) then
!
! I have not checked the units in this section. BD
!       
       do ik = 1, naky
          do it = 1, ntheta0
!             it = 2
             if (dump_check1) then
                denom=sum(delthet(-ntg_out:ntg_out-1)*jacob(-ntg_out:ntg_out-1))
                write (dump_check1_unit, "(20(1x,e12.5))") &
                     t, aky(ik), akx(it), &
                     sum(phinew(-ntg_out:ntg_out-1,it,ik) &
                     *delthet(-ntg_out:ntg_out-1) &
                         *jacob(-ntg_out:ntg_out-1))/denom, &
                     sum(phinew(-ntg_out:ntg_out-1,it,ik) &
                         *delthet(-ntg_out:ntg_out-1) &
                         *jacob(-ntg_out:ntg_out-1)) &
                         /denom/sourcefac
             end if
             if (dump_check2) then
                write (dump_check2_unit, "(5(1x,e12.5))") &
                     t, aky(ik), akx(it), aparnew(igomega,it,ik)
             end if
          end do
       end do
    end if

    if (dump_fields_periodically .and. mod(istep,10*nwrite) == 0) then
       call get_unused_unit (unit)
       write (filename, "('dump.fields.t=',e12.6)") t
       open (unit=unit, file=filename, status="unknown")
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntg_out, ntg_out
                write (unit=unit, fmt="(20(1x,e12.5))") &
                     theta(ig), aky(ik), theta0(it,ik), &
                     phinew(ig,it,ik), aparnew(ig,it,ik), &
                     bparnew(ig,it,ik), t, akx(it)
             end do
             write (unit, "()")
          end do
       end do
       close (unit=unit)
    end if
    
    ! Update the counter in parameter_scan_arrays
    scan_nout = nout

    nout = nout + 1


    if (write_ascii .and. mod(nout, 10) == 0 .and. proc0) then
       call flush_output_file (out_unit, ".out")
       if (write_verr) then
          call flush_output_file (res_unit, ".vres")
          call flush_output_file (lpc_unit, ".lpc")
       end if
       if (write_parity) call flush_output_file (parity_unit, ".parity")
    end if

! Deallocate variable for Jexternal
    if (write_jext) deallocate(j_ext)
!<GGH

    t_old = t
    call prof_leaving ("loop_diagnostics-2")

if (debug) write(6,*) "loop_diagnostics: done"
  end subroutine loop_diagnostics

  subroutine heating (istep, h, hk)

    use mp, only: proc0, iproc
    use dist_fn, only: get_heat
    use fields, only: phi, apar, bpar, phinew, aparnew, bparnew
    use species, only: nspec, spec
    use kt_grids, only: naky, ntheta0, aky, akx
    use theta_grid, only: ntgrid, delthet, jacob
    use nonlinear_terms, only: nonlin
    use dist_fn_arrays, only: c_rate
    use le_grids, only: negrid, w
    use gs2_heating, only: heating_diagnostics, avg_h, avg_hk, avg_he, zero_htype
    use general_f0, only: df0dE
    implicit none
    integer, intent (in) :: istep
    type (heating_diagnostics) :: h
    type (heating_diagnostics), dimension(:,:) :: hk

    real, dimension(-ntgrid:ntgrid) :: wgt
    real :: fac
    integer :: is, ik, it, ig, ie
    
    !Zero out variables for heating diagnostics
    call zero_htype(h)
    call zero_htype(hk)

    if (proc0) then
       
       !GGH NOTE: Here wgt is 1/(2*ntgrid+1)
       wgt = delthet*jacob
       wgt = wgt/sum(wgt)
          
       do is = 1, nspec
          do ik = 1, naky
             fac = 0.5
             if (aky(ik) < epsilon(0.)) fac = 1.0
             do it = 1, ntheta0
                if (aky(ik) < epsilon(0.0) .and. abs(akx(it)) < epsilon(0.0)) cycle
                do ig = -ntgrid, ntgrid
                   
                   !Sum heating by k over all z points (ig)
                   hk(it, ik) % collisions(is) = hk(it, ik) % collisions(is) &
                        + real(c_rate(ig,it,ik,is,1))*fac*wgt(ig)*spec(is)%temp*spec(is)%dens

                   hk(it, ik) % hypercoll(is) = hk(it, ik) % hypercoll(is) &
                        + real(c_rate(ig,it,ik,is,2))*fac*wgt(ig)*spec(is)%temp*spec(is)%dens

                   hk(it, ik) % imp_colls(is) = hk(it, ik) % imp_colls(is) &
                        + real(c_rate(ig,it,ik,is,3))*fac*wgt(ig)*spec(is)%temp*spec(is)%dens

                end do
                h % collisions(is) = h % collisions(is) + hk(it, ik) % collisions(is)
                h % hypercoll(is)  = h % hypercoll(is)  + hk(it, ik) % hypercoll(is)
                h % imp_colls(is)  = h % imp_colls(is)  + hk(it, ik) % imp_colls(is)
             end do
          end do

       end do
    end if

    call get_heat (h, hk, phi, apar, bpar, phinew, aparnew, bparnew)    

    call avg_h(h, h_hist, istep, navg)
    call avg_hk(hk, hk_hist, istep, navg)

  end subroutine heating

  subroutine heating_e(istep,he)
    use dist_fn, only: get_heat_e
    use fields, only: phi, apar, bpar, phinew, aparnew, bparnew
    use gs2_heating, only: heating_diagnostics, avg_he, zero_hetype
    implicit none
    integer, intent (in) :: istep
    type (heating_diagnostics), dimension(:) :: he

    !Zero out variables for heating diagnostics
    call zero_hetype(he)

    call get_heat_e (he, phi, apar, bpar, phinew, aparnew, bparnew)    
    call avg_he(he, he_hist, istep, navg)

  end subroutine heating_e
!=============================================================================
! Density: Calculate Density perturbations
!=============================================================================
 subroutine calc_jext (istep, j_ext)
    use mp, only: proc0
    use dist_fn, only: get_jext
    implicit none
    !Passed
    integer, intent (in) :: istep
    real, dimension(:,:) ::  j_ext
    !Local 
    integer :: i

    !Call routine to calculate density and velocity perturbations
    call get_jext(j_ext)    
    
    !Do averages with a history variable
    if (proc0) then
       !Save variable to history
       if (navg > 1) then
          if (istep > 1) &
               j_ext_hist(:,:,mod(istep,navg))= j_ext(:,:)

          !Use average of history
          if (istep >= navg) then
             j_ext=0.
             do i=0,navg-1
                j_ext(:,:) = j_ext(:,:) + j_ext_hist(:,:,i)/ real(navg)
             end do
          end if
       end if
    end if

  end subroutine calc_jext
!=============================================================================
!<GGH

  subroutine get_omegaavg (istep, exit, omegaavg, debopt)
    use kt_grids, only: naky, ntheta0
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use gs2_time, only: code_dt
    use constants, only: zi
    implicit none
    integer, intent (in) :: istep
    logical, intent (in out) :: exit
    complex, dimension (:,:), intent (out) :: omegaavg
!CMR, 7/11/2008: save here crucial to avoid crippling memory leak with mpixl!
    complex, allocatable, save, dimension (:,:,:) :: domega
    integer :: j
    logical, optional :: debopt
    logical :: debug=.false.

    !write (*,*) 'size of omegahist', size(omegahist)
if (.not. allocated(domega)) allocate(domega(navg,ntheta0,naky))
    if (present(debopt)) debug=debopt
if (debug) write(6,*) "get_omeaavg: start"
    j = igomega
    where (abs(phinew(j,:,:)+aparnew(j,:,:)+bparnew(j,:,:)) < epsilon(0.0) &
           .or. abs(phi(j,:,:)+apar(j,:,:)+bpar(j,:,:)) < epsilon(0.0))
       omegahist(mod(istep,navg),:,:) = 0.0
    elsewhere
       omegahist(mod(istep,navg),:,:) &
            = log((phinew(j,:,:) + aparnew(j,:,:) + bparnew(j,:,:)) &
                  /(phi(j,:,:)   + apar(j,:,:)    + bpar(j,:,:)))*zi/code_dt
    end where
    omegaavg = sum(omegahist/real(navg),dim=1)
if (debug) write(6,*) "get_omegaavg: omegaavg=",omegaavg

    if (istep > navg) then
       domega = spread(omegaavg,1,navg) - omegahist
       if (all(sqrt(sum(abs(domega)**2/real(navg),dim=1)) &
            < min(abs(omegaavg),1.0)*omegatol)) &
       then
          if (write_ascii) write (out_unit, "('*** omega converged')")
          exit = .true. .and. exit_when_converged
       end if

       if (any(abs(omegaavg)*code_dt > omegatinst)) then
          if (write_ascii) write (out_unit, "('*** numerical instability detected')") 
          exit = .true.
       end if
    end if
if (debug) write(6,*) "get_omegaavg: done"
  end subroutine get_omegaavg

!NOTE: Here we calculate the correction factors for each possible kperp
  subroutine get_vol_average_all (a, b, axb, axb_by_mode)
    use theta_grid, only: ntgrid, delthet, jacob
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: a, b
    real, intent (out) :: axb
    real, dimension (:,:), intent (out) :: axb_by_mode

    integer :: ik, it
    integer :: ng
    real, dimension (-ntg_out:ntg_out) :: wgt
    real :: anorm

    ng = ntg_out
    wgt = delthet(-ng:ng)*jacob(-ng:ng)
    anorm = sum(wgt)

    do ik = 1, naky
       do it = 1, ntheta0
          axb_by_mode(it,ik) &
               = sum(real(conjg(a(-ng:ng,it,ik))*b(-ng:ng,it,ik))*wgt)/anorm
       end do
    end do

    call get_volume_average (axb_by_mode, axb)
  end subroutine get_vol_average_all

  subroutine get_vol_average_one (a, b, axb, axb_by_mode)

    implicit none
    complex, dimension (:,:), intent (in) :: a, b
    real, intent (out) :: axb
    real, dimension (:,:), intent (out) :: axb_by_mode

    axb_by_mode = real(conjg(a)*b)

    call get_volume_average (axb_by_mode, axb)

  end subroutine get_vol_average_one

  subroutine get_volume_average (f, favg)
    use mp, only: iproc
    use kt_grids, only: naky, ntheta0, aky
    implicit none
    real, dimension (:,:), intent (in) :: f
    real, intent (out) :: favg
    real :: fac
    integer :: ik, it

! ky=0 modes have correct amplitudes; rest must be scaled
! note contrast with scaling factors in FFT routines.

    favg = 0.
    do ik = 1, naky
       fac = 0.5
       if (aky(ik) == 0.) fac = 1.0
       do it = 1, ntheta0
          favg = favg + f(it, ik) * fac
       end do
    end do

  end subroutine get_volume_average

  subroutine get_surf_average (f, favg)
    use mp, only: iproc
    use kt_grids, only: naky, ntheta0, aky
    implicit none
    real, dimension (:,:), intent (in) :: f
    real, dimension (:), intent (out) :: favg
    real :: fac
    integer :: ik, it

! ky=0 modes have correct amplitudes; rest must be scaled
! note contrast with scaling factors in FFT routines.

    favg = 0.
    do ik = 1, naky
       fac = 0.5
       if (aky(ik) == 0.) fac = 1.0
       do it = 1, ntheta0
          favg(it) = favg(it) + f(it, ik) * fac
       end do
    end do

  end subroutine get_surf_average

  subroutine reset_init

    use gs2_time, only: user_time

    implicit none

    start_time = user_time
    pflux_avg = 0.0 ; qflux_avg = 0.0 ; heat_avg = 0.0 ; vflux_avg = 0.0

  end subroutine reset_init

  subroutine ensemble_average (nensembles, time_int)
    use mp, only: scope, allprocs, group_to_all, broadcast
    use gs2_time, only: user_time
    use species, only: nspec
    implicit none
    integer, intent (in) :: nensembles
    real, intent (out) :: time_int
    integer :: is
    real, dimension (nensembles) :: dt_global
    real, dimension (nensembles,nspec) :: pflx_global, qflx_global, heat_global, vflx_global
    time_int=user_time-start_time
    call scope (allprocs)
    call group_to_all (time_int, dt_global, nensembles)
    call broadcast (dt_global)
    time_int = sum(dt_global)
    call group_to_all (pflux_avg, pflx_global, nensembles)
    call group_to_all (qflux_avg, qflx_global, nensembles)
    call group_to_all (vflux_avg, vflx_global, nensembles)
    do is = 1, nspec
       call broadcast (pflx_global(:,is))
       call broadcast (qflx_global(:,is))
       call broadcast (vflx_global(:,is))
       pflux_avg = sum(pflx_global(:,is))
       qflux_avg = sum(qflx_global(:,is))
       vflux_avg = sum(vflx_global(:,is))
    end do
    if (write_hrate) then
       call group_to_all (heat_avg, heat_global, nensembles)
       do is = 1, nspec
          call broadcast (heat_global(:,is))
          heat_avg = sum(heat_global(:,is))
       end do
    end if
  end subroutine ensemble_average

  subroutine get_fldline_avg_r (fld_in, fld_out)

    use theta_grid, only: delthet, jacob

    implicit none

    real, dimension (:), allocatable :: dl_over_b

    real, dimension (-ntg_out:), intent (in) :: fld_in
    real, intent (out) :: fld_out

    allocate (dl_over_b(-ntg_out:ntg_out))

    dl_over_b = delthet(-ntg_out:ntg_out)*jacob(-ntg_out:ntg_out)
    dl_over_b = dl_over_b / sum(dl_over_b)

    fld_out = sum(fld_in*dl_over_b)

    deallocate (dl_over_b)

  end subroutine get_fldline_avg_r

  subroutine get_fldline_avg_c (fld_in, fld_out)

    use theta_grid, only: delthet, jacob

    implicit none

    real, dimension (:), allocatable :: dl_over_b

    complex, dimension (-ntg_out:), intent (in) :: fld_in
    complex, intent (out) :: fld_out

    allocate (dl_over_b(-ntg_out:ntg_out))

    dl_over_b = delthet(-ntg_out:ntg_out)*jacob(-ntg_out:ntg_out)
    dl_over_b = dl_over_b / sum(dl_over_b)

    fld_out = sum(fld_in*dl_over_b)

    deallocate (dl_over_b)

  end subroutine get_fldline_avg_c

  subroutine get_cross_phase (phase_tot, phase_theta)

! <doc> This is a highly simplified synthetic diagnostic which 
! calculates the cross phase between the electron density and the 
! perpendicular electron temperature for comparisons with DIII-D.  
! Returns the value of the cross-phase at the outboard midplane and 
! integrated over all v and x. We can generalize this routine to 
! other fields at some point, but for now this is just a skeleton for 
! a more realistic synthetic diagnostic. </doc>

    use species, only: nspec, spec, electron_species
    use kt_grids, only: ntheta0, naky
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: gnew
    use dist_fn, only: getemoms
    use mp, only: proc0

    implicit none
    real, intent (out) :: phase_tot, phase_theta
    complex, dimension (:,:,:,:), allocatable :: ntot, tperp
    complex, dimension (ntheta0, naky) :: nTp_by_mode
    complex :: nTp
    real, dimension (ntheta0, naky) :: n2_by_mode, T2_by_mode
    real :: n2, T2

    integer :: it, ik, is, isgn, ig
    integer :: iglo

    allocate ( ntot(-ntgrid:ntgrid,ntheta0,naky,nspec))
    allocate (tperp(-ntgrid:ntgrid,ntheta0,naky,nspec))

    call getemoms (ntot, tperp)

    do is = 1,nspec
       if (spec(is)%type == electron_species) then
          ! get cross_phase at outboard midplane
          call get_vol_int (ntot(0,:,:,is), tperp(0,:,:,is), nTp, nTp_by_mode)
!          call get_vol_average (ntot(0,:,:,is), ntot(0,:,:,is), n2, n2_by_mode)
!          call get_vol_average (tperp(0,:,:,is), tperp(0,:,:,is), T2, T2_by_mode)
          phase_theta = atan2(aimag(nTp),real(nTp))!/sqrt(n2*T2)
          ! get integrated cross_phase 
          call get_vol_int (ntot(:,:,:,is), tperp(:,:,:,is), nTp, nTp_by_mode)
!          call get_vol_average (ntot(:,:,:,is), ntot(:,:,:,is), n2, n2_by_mode)
!          call get_vol_average (tperp(:,:,:,is), tperp(:,:,:,is), T2, T2_by_mode)
          phase_tot = atan2(aimag(nTp),real(nTp))!/sqrt(n2*T2)
       end if
    end do

    deallocate (ntot, tperp)

  end subroutine get_cross_phase

  subroutine get_vol_int_all (a, b, axb, axb_by_mode)
    use theta_grid, only: ntgrid, delthet, jacob
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: a, b
    complex, intent (out) :: axb
    complex, dimension (:,:), intent (out) :: axb_by_mode

    integer :: ik, it
    integer :: ng
    real, dimension (-ntg_out:ntg_out) :: wgt
    real :: anorm

    ng = ntg_out
    wgt = delthet(-ng:ng)*jacob(-ng:ng)
    anorm = sum(wgt)

    do ik = 1, naky
       do it = 1, ntheta0
          axb_by_mode(it,ik) &
               = sum((conjg(a(-ng:ng,it,ik))*b(-ng:ng,it,ik))*wgt)/anorm
       end do
    end do

    call get_volume_int (axb_by_mode, axb)
  end subroutine get_vol_int_all

  subroutine get_vol_int_one (a, b, axb, axb_by_mode)

    implicit none
    complex, dimension (:,:), intent (in) :: a, b
    complex, intent (out) :: axb
    complex, dimension (:,:), intent (out) :: axb_by_mode

    axb_by_mode = conjg(a)*b
    call get_volume_int (axb_by_mode, axb)

  end subroutine get_vol_int_one

  subroutine get_volume_int (f, favg)
    use mp, only: iproc
    use kt_grids, only: naky, ntheta0, aky
    implicit none
    complex, dimension (:,:), intent (in) :: f
    complex, intent (out) :: favg
    real :: fac
    integer :: ik, it

! ky=0 modes have correct amplitudes; rest must be scaled
! note contrast with scaling factors in FFT routines.

    favg = 0.
    do ik = 1, naky
       fac = 0.5
       if (aky(ik) == 0.) fac = 1.0
       do it = 1, ntheta0
          favg = favg + f(it, ik) * fac
       end do
    end do

  end subroutine get_volume_int

!  subroutine autocorrelation (cfnc,phi2extend,cfnc_2pi)
  subroutine correlation_extend (cfnc,phi2extend)

    use constants, only: pi
    use fields_arrays, only: phinew
    use theta_grid, only: ntgrid, theta, jacob, delthet
    use kt_grids, only: theta0, ntheta0, naky, jtwist_out

    implicit none

    real, dimension (:,:,:), intent (out) :: phi2extend
    complex, dimension (:,:,:), intent (out) :: cfnc

    integer :: ig, it, ik, im, igmod
    integer :: itshift, nconnect, offset
    real :: fac

    real, dimension (:), allocatable :: dl_over_b
    complex, dimension (:,:,:), allocatable :: phiextend
    complex, dimension (:,:), allocatable :: phir
!    real, dimension (:), allocatable :: phisum

    allocate (dl_over_b(2*ntgrid+1))
    dl_over_b = delthet*jacob

!    allocate (theta_extend(ntg_extend)) ; theta_extend = 0.0
!    do ig = 1, (ntheta0-1)/jtwist_out+1
!       theta_extend((ig-1)*(2*ntgrid+1)+1:ig*(2*ntgrid+1)) = pi+theta+(ig-1)*2.*pi
!    end do
!    theta_extend = theta_extend - theta_extend(size(theta_extend))/2.

    allocate (phiextend(ntg_extend,ntheta0,naky)) ; phiextend = 0.0
    allocate (phir(-ntgrid:ntgrid,ntheta0))
!    allocate (phisum(naky))

!    phisum = 0.0 ; cfnc_2pi = 0.0
!     do ik = 1, naky
!        if (ik==1) then
!           fac = 1.0
!        else
!           fac = 0.5
!        end if
!        ! integrate over x
!        do it = 1, ntheta0
!           do ig = -ntgrid, ntgrid
!              cfnc_2pi(ig,ik) = cfnc_2pi(ig,ik) + phinew(0,it,ik)*conjg(phinew(ig,it,ik))*fac
!           end do
! !          phisum(ik) = phisum(ik) + phinew(0,it,ik)*conjg(0,it,ik)*fac
!        end do
! !       cfnc_2pi(:,ik) = cfnc_2pi(:,ik)/phisum(ik)
!     end do

    cfnc = 0.0 ; phiextend = 0.0
    offset = ((ntheta0-1)/jtwist_out)*(2*ntgrid+1)/2
    do ig = -ntgrid, ntgrid
       call reorder_kx (phinew(ig,:,1), phir(ig,:))
    end do
    phiextend(offset+1:offset+(2*ntgrid+1),:,1) = phir
    do it = 1, ntheta0
       do im = 1, 2*ntgrid+1
          do ig = im+offset, offset+(2*ntgrid+1)
             igmod = mod(ig-offset-1,2*ntgrid+1)+1
             cfnc(im,it,1) = cfnc(im,it,1) &
                  + phiextend(ig,it,1)*conjg(phiextend(ig-im+1,it,1)) &
                  * dl_over_b(igmod)
          end do
       end do
       cfnc(:,it,1) = cfnc(:,it,1) / cfnc(1,it,1)
    end do

    do ik = 2, naky
       do ig = -ntgrid, ntgrid
          call reorder_kx (phinew(ig,:,ik), phir(ig,:))
       end do
       ! shift in kx due to parallel boundary condition
       ! also the number of independent theta0s
       itshift = jtwist_out*(ik-1)
       do it = 1, min(itshift,ntheta0)
          ! number of connections between kx's
          nconnect = (ntheta0-it)/itshift
          ! shift of theta index to account for fact that not all ky's
          ! have same length in extended theta
          offset = (2*ntgrid+1)*((ntheta0-1)/jtwist_out-nconnect)/2
          do ig = 1, nconnect+1
             phiextend(offset+(ig-1)*(2*ntgrid+1)+1:offset+ig*(2*ntgrid+1),it,ik) &
                     = phir(:,ntheta0-it-(ig-1)*itshift+1)
          end do
          do im = 1, (2*ntgrid+1)*(nconnect+1)
             do ig = im+offset, offset+(2*ntgrid+1)*(nconnect+1)
                igmod = mod(ig-offset-1,2*ntgrid+1)+1
                cfnc(im,it,ik) = cfnc(im,it,ik) &
                     + phiextend(ig,it,ik)*conjg(phiextend(ig-im+1,it,ik)) &
                     * dl_over_b(igmod)
             end do
          end do
          cfnc(:,it,ik) = cfnc(:,it,ik) / cfnc(1,it,ik)
       end do
    end do
    
    phi2extend = phiextend*conjg(phiextend)

!    deallocate (dl_over_b, phir, phiextend, phisum)
    deallocate (dl_over_b, phir, phiextend)

  end subroutine correlation_extend

  subroutine correlation (cfnc_2pi)

    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use fields_arrays, only: phinew

    implicit none

    complex, dimension (-ntgrid:,:), intent (out) :: cfnc_2pi

    integer :: ik, it, ig
    real :: fac

    cfnc_2pi = 0.0

    do ik = 1, naky
       if (ik==1) then
          fac = 1.0
       else
          fac = 0.5
       end if
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid
             cfnc_2pi(ig,ik) = cfnc_2pi(ig,ik) + phinew(0,it,ik)*conjg(phinew(ig,it,ik))*fac
          end do
       end do
    end do

  end subroutine correlation

  subroutine reorder_kx (unord, ord)

    use kt_grids, only: ntheta0

    implicit none

    complex, dimension (:), intent (in) :: unord
    complex, dimension (:), intent (out) :: ord

    ord(:ntheta0/2) = unord(ntheta0/2+2:)
    ord(ntheta0/2+1:) = unord(:ntheta0/2+1)

  end subroutine reorder_kx

end module gs2_diagnostics
