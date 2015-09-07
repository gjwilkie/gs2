! what about boundary condition contributions?  In the presence 
! of magnetic shear, there are not always zero amplitudes at the ends
! of the supercells.

module gs2_diagnostics

  use gs2_save, only: save_many

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

  interface get_vol_average
     module procedure get_vol_average_one, get_vol_average_all
  end interface

!  interface get_vol_int
!     module procedure get_vol_int_one, get_vol_int_all
!  end interface

!  interface get_fldline_avg
!     module procedure get_fldline_avg_r, get_fldline_avg_c
!  end interface

!CMR, 17/11/2009:   read_parameters now public so ingen can USE instead of copy
!
! Why are these variables public?  This is not good.
  real, public :: omegatol, omegatinst
  logical, public :: print_line, print_flux_line
  logical, public :: print_summary, write_line, write_flux_line
  logical, public :: write_omega, write_omavg, write_ascii
  logical, public :: write_gs
  logical, public :: write_g, write_gg, write_gyx
  logical, public :: write_eigenfunc, write_fields, write_final_fields, write_final_antot
  logical, public :: write_final_moments, write_avg_moments
  logical, public :: write_moments, write_final_db
  logical, public :: write_kpar
  logical, public :: write_lorentzian
  logical, public :: write_nl_flux
  logical, public :: exit_when_converged
  logical, public :: dump_check1, dump_check2
  logical, public :: dump_fields_periodically, make_movie
  logical, public :: save_for_restart
  logical, public :: write_symmetry
  logical, public :: write_correlation_extend, write_correlation
  integer, public :: nwrite, igomega, nmovie
  integer, public :: navg, nsave, nwrite_mult

  logical, public :: write_phi_over_time, write_apar_over_time, write_bpar_over_time !EGH

  ! internal
  logical :: write_any, write_any_fluxes, dump_any
  logical, private :: initialized = .false.

  namelist /gs2_diagnostics_knobs/ print_line, print_flux_line, &
         write_line, write_flux_line, &
         write_omega, write_omavg, write_ascii, write_kpar, &
         write_gs, write_gyx, write_g, write_gg, &
         write_eigenfunc, write_fields, write_final_fields, write_final_antot, &
         write_moments, write_final_moments, &
         write_nl_flux, write_final_db, &
         nwrite, nmovie, nsave, navg, omegatol, omegatinst, igomega, write_lorentzian, &
         exit_when_converged, write_avg_moments, &
         dump_check1, dump_check2, &
         dump_fields_periodically, make_movie, &
         save_for_restart, save_many, &
         write_symmetry, &
         write_correlation_extend, nwrite_mult, write_correlation, &
         write_phi_over_time, write_apar_over_time, write_bpar_over_time

  integer :: out_unit, kp_unit, heat_unit, polar_raw_unit, polar_avg_unit, heat_unit2, funit
  integer :: phase_unit
  integer :: dump_check1_unit, dump_check2_unit

  complex, dimension (:,:,:), allocatable :: omegahist
  ! (navg,ntheta0,naky)

  real, dimension (:,:,:,:), allocatable ::  qheat, qmheat, qbheat, qheath
  ! (ntheta0,naky,nspec,3)

  real, dimension (:,:,:), allocatable ::  pflux,  vflux, vflux_par, vflux_perp
  real, dimension (:,:,:), allocatable :: pfluxh, vfluxh
  real, dimension (:,:,:), allocatable :: vflux0, vflux1  ! low flow correction to turbulent momentum flux
  real, dimension (:,:,:), allocatable :: pmflux, vmflux
  real, dimension (:,:,:), allocatable :: pbflux, vbflux
  real, dimension (:,:,:), allocatable :: exchange

  ! (ntheta0,naky,nspec)

  real :: start_time = 0.0
  real, dimension (:), allocatable :: pflux_avg, qflux_avg, heat_avg, vflux_avg

  integer :: ntg_out, ntg_extend, nth0_extend
  integer :: nout = 1
  integer :: nout_movie = 1
  integer :: nout_big = 1
  complex :: wtmp_old = 0.
  logical :: exist

  character (3) :: nspec_str
  character (100) :: fluxes_str
  
contains
  !> Define NetCDF vars, call real_init, which calls read_parameters; broadcast all the different write flags. 
   subroutine wnml_gs2_diagnostics(unit)
   implicit none
   integer :: unit
       if (.not.exist) return
       write (unit, *)
       write (unit, fmt="(' &',a)") "gs2_diagnostics_knobs"
       write (unit, fmt="(' save_for_restart = ',L1)") save_for_restart
       write (unit, fmt="(' save_many = ',L1)") save_many
       write (unit, fmt="(' print_line = ',L1)") print_line 
       write (unit, fmt="(' write_line = ',L1)") write_line
       write (unit, fmt="(' print_flux_line = ',L1)") print_flux_line
       write (unit, fmt="(' write_flux_line = ',L1)") write_flux_line
       write (unit, fmt="(' nmovie = ',i6)") nmovie
       write (unit, fmt="(' nwrite_mult = ',i6)") nwrite_mult
       write (unit, fmt="(' nwrite = ',i6)") nwrite
       write (unit, fmt="(' nsave = ',i6)") nsave
       write (unit, fmt="(' navg = ',i6)") navg
       write (unit, fmt="(' omegatol = ',e17.10)") omegatol
       write (unit, fmt="(' omegatinst = ',e17.10)") omegatinst
! should be legal -- not checked yet
       if (igomega /= 0) write (unit, fmt="(' igomega = ',i6)") igomega  
       
       if (write_ascii) then
          write (unit, fmt="(' write_ascii = ',L1)") write_ascii
          write (unit, fmt="(' write_omega = ',L1)") write_omega
          write (unit, fmt="(' write_omavg = ',L1)") write_omavg
       end if
       write (unit, fmt="(' write_lorentzian = ',L1)") write_lorentzian
       write (unit, fmt="(' write_eigenfunc = ',L1)") write_eigenfunc
       write (unit, fmt="(' write_final_fields = ',L1)") write_final_fields
       write (unit, fmt="(' write_final_db = ',L1)") write_final_db
       write (unit, fmt="(' write_final_moments = ',L1)") write_final_moments
       write (unit, fmt="(' write_final_antot = ',L1)") write_final_antot
       write (unit, fmt="(' write_nl_flux = ',L1)") write_nl_flux
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
       write (report_unit, fmt="('print_line = T: Estimated frequencies output to the screen every ',i4,' steps.')") nwrite
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
       write (report_unit, fmt='(a22,a45,i4,a7)') 'print_flux_line = T: ', &
            'instantaneous fluxes output to screen every ', &
            nwrite, ' steps.'
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
    use kt_grids, only: init_kt_grids
    use run_parameters, only: init_run_parameters
    use species, only: init_species, nspec
    use dist_fn, only: init_dist_fn
    use init_g, only: init_init_g
    use gs2_io, only: init_gs2_io
    use mp, only: broadcast

    implicit none
    logical, intent (in) :: list
    integer, intent (in) :: nstep
    integer :: nmovie_tot, nwrite_big_tot

    if (initialized) return
    initialized = .true.

    call init_theta_grid
    call init_kt_grids
    call init_run_parameters
    call init_species
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
!    call broadcast (write_cross_phase)
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
    call broadcast (save_many)
    call broadcast (write_gs)
    call broadcast (write_g)
    call broadcast (write_gyx)
    call broadcast (write_gg)
    call broadcast (write_final_antot)

    call broadcast (ntg_out)
    call broadcast (write_lorentzian)
    call broadcast (write_eigenfunc)

    call broadcast (write_phi_over_time)
    call broadcast (write_apar_over_time)
    call broadcast (write_bpar_over_time)

    nmovie_tot = nstep/nmovie
    nwrite_big_tot = nstep/(nwrite*nwrite_mult)-nstep/4/(nwrite*nwrite_mult)

    call init_gs2_io (write_nl_flux, write_omega, &
         write_final_antot, &
         write_eigenfunc, make_movie, nmovie_tot, &
         write_fields, write_moments, &
         write_symmetry, &
         write_correlation, nwrite_big_tot, write_correlation_extend, &
         write_phi_over_time, write_apar_over_time, write_bpar_over_time)
    
    allocate (pflux_avg(nspec), qflux_avg(nspec), heat_avg(nspec), vflux_avg(nspec))
    pflux_avg = 0.0 ; qflux_avg = 0.0 ; heat_avg = 0.0 ; vflux_avg = 0.0

    write (nspec_str,'(i3)') nspec
    fluxes_str = trim('(a3,e17.9,a11,e12.4,a13,'//trim(nspec_str)//'(1x,e12.4),a11,'//trim(nspec_str)//'(1x,e12.4))')

  end subroutine init_gs2_diagnostics
 

  subroutine real_init (list)
    use run_parameters, only: fapar
    use file_utils, only: open_output_file, get_unused_unit
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use mp, only: proc0
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

!       if (write_cross_phase .and. write_ascii) then
!          call open_output_file (phase_unit, ".phase")
!       end if

       if (dump_check1) then
          call get_unused_unit (dump_check1_unit)
          open (unit=dump_check1_unit, file="dump.check1", status="unknown")
       end if
       
       ! TMP FOR TESTING -- MAB
!       call open_output_file (funit,".funit")

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

    !<doc> Allocate arrays for storing the various fluxes which the diagnostics will output </doc>
    allocate (pflux (ntheta0,naky,nspec)) ; pflux = 0.
    allocate (qheat (ntheta0,naky,nspec,3)) ; qheat = 0.
    allocate (vflux (ntheta0,naky,nspec)) ; vflux = 0.
    allocate (pfluxh (ntheta0,naky,nspec)) ; pfluxh = 0.
    allocate (qheath (ntheta0,naky,nspec,3)) ; qheath = 0.
    allocate (vfluxh (ntheta0,naky,nspec)) ; vfluxh = 0.
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
      
  end subroutine real_init

  subroutine read_parameters (list)
!CMR, 17/11/2009:   namelist gs2_diagnostics_knobs made public
!                   so that ingen can just USE it instead of copying it!
!
    use file_utils, only: input_unit, input_unit_exist
    use theta_grid, only: nperiod, ntheta
    use kt_grids, only: box
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
       write_gs = .false.
       write_g = .false.
       write_gyx = .false.
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
       write_symmetry = .false.
       write_correlation_extend = .false.
       write_correlation = .false.
       write_fields = .false.
       write_final_fields = .false.
       write_final_antot = .false.
       write_final_db = .false.
       nwrite = 100
       nmovie = 100000
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
       save_many = .false.
       write_phi_over_time = .false.
       write_bpar_over_time = .false.
       write_apar_over_time = .false.
       in_file = input_unit_exist ("gs2_diagnostics_knobs", exist)

	!<doc> Read in parameters from the namelist gs2_diagnostics_knobs, if the namelist exists </doc>
!       if (exist) read (unit=input_unit("gs2_diagnostics_knobs"), nml=gs2_diagnostics_knobs)
       if (exist) read (unit=in_file, nml=gs2_diagnostics_knobs)

       print_summary = (list .and. (print_line .or. print_flux_line)) 

       if (list) then
          print_line = .false.
          print_flux_line = .false.
       end if

       if (.not. save_for_restart) nsave = -1
       write_avg_moments = write_avg_moments .and. box

       write_any = write_line .or. write_omega     .or. write_omavg &
            .or. write_flux_line                   .or. write_nl_flux  &
            .or. write_kpar   .or. write_lorentzian  .or. write_gs
       write_any_fluxes =  write_flux_line .or. print_flux_line .or. write_nl_flux 
       dump_any = dump_check1  .or. dump_fields_periodically &
            .or.  dump_check2 .or. make_movie .or. print_summary

       ntg_out = ntheta/2 + (nperiod-1)*ntheta
    end if 
 end subroutine read_parameters

  subroutine finish_gs2_diagnostics
    use file_utils, only: open_output_file, close_output_file, get_unused_unit
    use mp, only: proc0, broadcast, nproc, iproc, sum_reduce
    use species, only: nspec, spec
    use run_parameters, only: fphi, fapar, fbpar
    use theta_grid, only: ntgrid, theta, delthet, jacob, gradpar, nperiod, bmag
    use theta_grid, only: Rplot, Zplot, aplot, Rprime, Zprime, aprime
    use theta_grid, only: drhodpsi, qval, shape
    use kt_grids, only: naky, ntheta0, theta0, nx, ny, akx, aky
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew
    use fields_arrays, only: phip
    use dist_fn, only: getan, getmoms, par_spectrum
!    use dist_fn, only: write_f, write_fyx, def_parity, even
    use dist_fn, only: write_f, def_parity, even
    use dist_fn_arrays, only: gnew
    use gs2_transforms, only: transform2, inverse2
    use gs2_save, only: gs2_save_for_restart
    use constants
    use gs2_time, only: user_time, user_dt, code_dt
    use gs2_io, only: nc_eigenfunc, nc_final_fields, nc_final_an
    use gs2_io, only: nc_final_moments, nc_finish
    use antenna, only: dump_ant_amp
    use splines, only: fitp_surf1, fitp_surf2
    implicit none
    integer :: ig, ik, it, is, unit, ierr
!    complex, dimension (:,:,:), allocatable :: xphi, xapar, xbpar
    real, dimension (:), allocatable :: total
    real, dimension (:,:,:), allocatable :: bxf, byf, vxf, vyf, bxfsavg, byfsavg
    real, dimension (:,:,:), allocatable :: bxfs, byfs, vxfs, vyfs, rvx, rvy, rx, ry
    complex, dimension (:,:,:), allocatable :: bx, by, vx, vy, vx2, vy2
    complex, dimension (:,:,:), allocatable :: phi2, apar2, bpar2, antot, antota, antotp
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

!    if (write_gyx) call write_fyx (phinew,bparnew,last)
    if (write_g) call write_f (last)
!    if (write_gg) call write_dist (g)

    phi0 = 1.

    if (proc0) then
       if (write_ascii) call close_output_file (out_unit)
!       if (write_ascii .and. write_cross_phase) call close_output_file (phase_unit)

       ! TMP FOR TESTING -- MAB
!       call close_output_file (funit)

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
             ! TMP FOR TESTING -- MAB
             call open_output_file (unit, ".fields2")
             do ik = 1, naky
                do it = 1, ntheta0
                   do ig = -ntg_out, ntg_out
                      write (unit, "(7(1x,e12.5))") &
                           theta(ig), aky(ik), akx(it), &
                           phip(ig,it,ik), &
                           theta(ig) - theta0(it,ik), &
                           cabs(phip(ig,it,ik))
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
       call getan (gnew, antot, antota, antotp)
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
       call gs2_save_for_restart (gnew, user_time, user_dt, istatus, &
            fphi, fapar, fbpar, exit_in=.true.)
    end if

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
       write (g_unit,fmt="('# q = ',e11.4,' drhodpsi = ',e11.4)") qval, drhodpsi
       write (g_unit,fmt='(7a18)') '# theta', 'R', 'Z', 'alpha', 'Rprim', 'Zprim', 'alpha_prim'
       do i=-ntgrid,ntgrid
          write (g_unit,fmt='(7(1x,1pg18.11))') theta(i),Rplot(i),Zplot(i),aplot(i), &
               Rprime(i),Zprime(i),aprime(i)
       enddo
       call close_output_file (g_unit)
    end if

    if (allocated(omegahist)) deallocate (omegahist)
    if (allocated(pflux)) deallocate (pflux, qheat, vflux, vflux_par, vflux_perp, pmflux, qmheat, vmflux, &
         pbflux, qbheat, vbflux, vflux0, vflux1, exchange, pfluxh, qheath, vfluxh)
    if (allocated(bxf)) deallocate (bxf, byf, xx4, xx, yy4, yy, dz, total)
    if (allocated(pflux_avg)) deallocate (pflux_avg, qflux_avg, heat_avg, vflux_avg)

    wtmp_old = 0. ; nout = 1 ; nout_movie = 1 ; nout_big = 1
    initialized = .false.

  end subroutine finish_gs2_diagnostics

  subroutine loop_diagnostics (istep, exit, debopt)
    use species, only: nspec, spec, has_electron_species
    use theta_grid, only: theta, ntgrid, delthet, jacob
    use kt_grids, only: naky, ntheta0, theta0, aky, akx
    use kt_grids, only: jtwist_out
    use run_parameters, only: woutunits, tunits, fapar, fphi, fbpar, eqzip
!    use run_parameters, only: nstep, include_lowflow
    use run_parameters, only: nstep
    use fields, only: phinew, aparnew, bparnew
!    use fields, only: phihnew, aparhnew!, bparhnew
    use fields, only: phipnew!, aparpnew, bparpnew
    use fields, only: kperp, phinorm
    use dist_fn, only: flux, write_f!, write_fyx
    use dist_fn, only: omega0, gamma0, getmoms, par_spectrum
    use dist_fn, only: eexchange
    use dist_fn, only: flux_vs_theta_vs_vpa
!# ifdef LOWFLOW
!    use dist_fn, only: lf_flux
!# endif
    use dist_fn_arrays, only: gnew, g_adjust!, ghnew
    use vpamu_grids, only: nvgrid
    use mp, only: proc0, broadcast, send, receive
    use file_utils, only: get_unused_unit, flush_output_file
    use gs2_time, only: user_time, code_dt
    use gs2_io, only: nc_qflux, nc_vflux, nc_pflux, nc_loop, nc_loop_moments
    use gs2_io, only: nc_loop_fullmom, nc_loop_corr, nc_loop_corr_extend
    use gs2_io, only: nc_loop_movie, nc_write_fields, nc_write_moments
    use gs2_io, only: nc_loop_sym
    use gs2_layouts, only: yxf_lo
    use gs2_layouts, only: idx, idx_local, proc_id, imu_idx, is_idx
    use gs2_transforms, only: init_transforms, transform2
    use nonlinear_terms, only: nonlin
    use antenna, only: antenna_w
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
!    complex, dimension (ntheta0, naky) :: omega, omegaavg
    complex, dimension (:, :), allocatable :: omega, omegaavg

    real, dimension (ntheta0, naky) :: phitot, akperp
    real :: phi2, apar2, bpar2, phip2
    real, dimension (ntheta0, naky) :: phi2_by_mode, apar2_by_mode, bpar2_by_mode, phip2_by_mode
    real, dimension (ntheta0, naky, nspec) :: ntot2_by_mode, ntot20_by_mode
    real, dimension (ntheta0, naky, nspec) :: tpar2_by_mode, tperp2_by_mode
    ! arrays needed for symmetry diagnostic
    real, dimension (:,:,:), allocatable :: pflx_sym, vflx_sym, qflx_sym
    real, dimension (:,:,:), allocatable :: pflxh_sym, vflxh_sym, qflxh_sym
    ! arrays needed for parallel autocorrelation diagnostic
    real, save :: tcorr0 = 0.0
    real, dimension (:,:,:), allocatable :: phi2_extend
    complex, dimension (:,:,:), allocatable :: phi_corr
    complex, dimension (:,:), allocatable :: phi_corr_2pi
    real, dimension (:,:,:), allocatable, save :: phiextend_sum
    complex, dimension (:,:,:), allocatable, save :: phicorr_sum
    real :: t, denom
    integer :: ig, ik, it, is, unit, i, nnx, nny
    complex :: sourcefac
!    complex :: phiavg
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky,nspec) :: ntot, density, &
         upar, tpar, tperp, qparflux, pperpj1, qpperpj1
    complex, dimension (ntheta0, nspec) :: ntot00, density00, upar00, tpar00, tperp00
    complex, dimension (ntheta0) :: phi00
    complex, save :: wtmp_new
!    complex :: wtmp_old = 0.
    real, dimension (:), allocatable :: dl_over_b
    real, dimension (ntheta0, nspec) :: x_qmflux
    real, dimension (nspec) :: ntot2, ntot20, tpar2, tperp2
    real, dimension (nspec) ::  heat_fluxes,  part_fluxes, mom_fluxes, parmom_fluxes, perpmom_fluxes
    real, dimension (nspec) ::  heath_fluxes,  parth_fluxes, momh_fluxes
# ifdef USE_LOWFLOW
!    real, dimension (nspec) :: lfmom_fluxes, vflux1_avg  ! low-flow correction to turbulent momentum fluxes
# endif
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
    real, save :: t_old = 0.
    character(200) :: filename
    logical :: last = .false.
    logical,optional:: debopt
    logical:: debug=.false.

    if (present(debopt)) debug=debopt

    part_fluxes = 0.0 ; mpart_fluxes = 0.0 ; bpart_fluxes = 0.0
    heat_fluxes = 0.0 ; mheat_fluxes = 0.0 ; bheat_fluxes = 0.0
    mom_fluxes = 0.0 ; mmom_fluxes = 0.0 ; bmom_fluxes = 0.0
    heath_fluxes = 0.0 ; parth_fluxes = 0.0 ; momh_fluxes = 0.0
    energy_exchange = 0.0

    phase_tot = 0.0 ;  phase_theta = 0.0

    allocate (omega(ntheta0,naky)) ; omega = 0.
    allocate (omegaavg(ntheta0,naky)) ; omegaavg = 0.

    exit = .false.

    if (eqzip .or. .not. nonlin) then
! MR, 10/3/2009: avoid calling get_omegaavg in nonlinear calculations
       if (proc0) then
          if (debug) write(6,*) "loop_diagnostics: proc0 call get_omegaavg"
          call get_omegaavg (istep, exit, omegaavg, debug)
          if (debug) write(6,*) "loop_diagnostics: proc0 done called get_omegaavg"
       endif
       call broadcast (exit)
    endif

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

!    if (write_gyx .and. mod(istep,nmovie) == 0) call write_fyx (phinew,bparnew,last)

    if (mod(istep,nwrite) /= 0 .and. .not. exit) return
    t = user_time

    if (write_g .and. mod(istep,nmovie) == 0) call write_f (last)

if (debug) write(6,*) "loop_diagnostics: -1"

    if (proc0) then
       omega = omegahist(mod(istep,navg),:,:)
       sourcefac = exp(-zi*omega0*t+gamma0*t)
       call phinorm (phinew, aparnew, bparnew, phitot)
       call get_vol_average (phinew, phinew, phi2, phi2_by_mode)
       call get_vol_average (phipnew, phipnew, phip2, phip2_by_mode)
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
       ! switch from g=<f> to h = nonadiabatic part of f
       call g_adjust (gnew, phinew, bparnew, fphi, fbpar)
!       call g_adjust (ghnew, phihnew, bparhnew, fphi, fbpar)

       ! get higher order correction to fluxes first
       ! call flux (ghnew, phinew, aparnew, bparnew, &
       !      pflux,  qheat,  vflux, vflux_par, vflux_perp, &
       !      pmflux, qmheat, vmflux, pbflux, qbheat, vbflux)
       ! call flux (gnew, phihnew, aparhnew, bparhnew, &
       !      pfluxh,  qheath,  vfluxh, vflux_par, vflux_perp, &
       !      pmflux, qmheat, vmflux, pbflux, qbheat, vbflux)
!       pfluxh = pflux + pfluxh
!       vfluxh = vflux + vfluxh
!       qheath = qheat + qheath
       pfluxh = 0. ; vfluxh = 0. ; qheath = 0.

       ! now get lowest order fluxes
       call flux (gnew, phinew, aparnew, bparnew, &
            pflux,  qheat,  vflux, vflux_par, vflux_perp, &
            pmflux, qmheat, vmflux, pbflux, qbheat, vbflux)
!#ifdef LOWFLOW
!       ! lowflow terms only implemented in electrostatic limit at present
!       call lf_flux (phinew, vflux0, vflux1)
!#endif

       ! convert back from h to g
       call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar)
!       call g_adjust (ghnew, phihnew, bparhnew, -fphi, -fbpar)

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

                qheath(:,:,is,1) = qheath(:,:,is,1) * spec(is)%dens*spec(is)%temp
                call get_volume_average (qheath(:,:,is,1), heath_fluxes(is))

                pfluxh(:,:,is) = pfluxh(:,:,is) * spec(is)%dens
                call get_volume_average (pfluxh(:,:,is), parth_fluxes(is))

                vfluxh(:,:,is) = vfluxh(:,:,is) * spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
                call get_volume_average (vfluxh(:,:,is), momh_fluxes(is))

                vflux_par(:,:,is) = vflux_par(:,:,is) * spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
                call get_volume_average (vflux_par(:,:,is), parmom_fluxes(is))

                vflux_perp(:,:,is) = vflux_perp(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%stm
                call get_volume_average (vflux_perp(:,:,is), perpmom_fluxes(is))

                exchange(:,:,is) = exchange(:,:,is) * spec(is)%dens*spec(is)%z
                call get_volume_average (exchange(:,:,is), energy_exchange(is))

#ifdef LOWFLOW
!                vflux0(:,:,is) = vflux0(:,:,is) * spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
!                call get_volume_average (vflux0(:,:,is), lfmom_fluxes(is))
!
!                vflux1(:,:,is) = vflux1(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%temp/spec(is)%z
!                call get_volume_average (vflux1(:,:,is), vflux1_avg(is))
!
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
!          t_old = t
       end if
    end if

    call broadcast (pflux_avg)
    call broadcast (qflux_avg)
    call broadcast (vflux_avg)

    fluxfac = 0.5
    fluxfac(1) = 1.0

    ! TMP FOR TESTING -- MAB
!     do iglo=g_lo%llim_world, g_lo%ulim_world
!        imu = imu_idx(g_lo,iglo) ; if (imu /= 4) cycle
!        is = is_idx(g_lo,iglo) ; if (is /= 1) cycle
!        do iv = -nvgrid, nvgrid
!           do ig = -ntgrid, ntgrid
!              write (funit,'(a9,6e12.4)') 'distfnc', t, theta(ig), vpa(iv), &
!                   real(gnew(ig,iv,1,iglo)), aimag(gnew(ig,iv,1,iglo)), bmag(ig)
!           end do
!           write (funit,*)
!        end do
!     end do

    if (proc0) then
       if (print_flux_line) then
          if (fphi > epsilon(0.0)) then
             write (unit=*, fmt=107) &
                  't=', t, '<phi**2>=', phi2, 'heat flux:', heat_fluxes(1:min(nspec,5))
             write (unit=*, fmt=107) &
                  't=', t, '<phip**2>=', phip2, 'eexchange:', energy_exchange(1:min(nspec,5))
          end if
          if (fapar > epsilon(0.0)) then
             write (unit=*, fmt=107) &
                  't=', t, '<apar**2>=', apar2, 'heat flux m:', mheat_fluxes(1:min(nspec,5))
          end if
          if (fbpar > epsilon(0.0)) then
             write (unit=*, fmt=107) &
                  't=', t, '<bpar**2>=', bpar2, 'heat flux b:', bheat_fluxes(1:min(nspec,5))
          end if
       end if
       if (print_line) then
          write (*,*) 'time: ', t
          do ik = 1, naky
             do it = 1, ntheta0
                write (unit=*, fmt=106) &
                     'ky=', aky(ik), 'kx=', akx(it), &
                     'om=', real( omega(it,ik)*woutunits(ik)), &
                     aimag(omega(it,ik)*woutunits(ik)), &
                     'omav=', real( omegaavg(it,ik)*woutunits(ik)), &
                     aimag(omegaavg(it,ik)*woutunits(ik)), &
                     'phtot=', phitot(it,ik)
             end do
          end do
          write (*,*) 
       end if
    end if

    i=istep/nwrite

    if (write_fields) call nc_write_fields (nout, phinew, aparnew, bparnew)  !MR
    if (write_moments) then !CMR
       call getmoms (ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)
       if (proc0) write (*,'(a7,3e12.4)') 'phiphi', t, real(phinew(0,1,1)), real(aparnew(0,1,1))
       call nc_write_moments(nout, ntot)
    endif

!     if (write_cross_phase .and. has_electron_species(spec)) then
!        call get_cross_phase (phase_tot, phase_theta)
!        if (proc0) write (unit=phase_unit, fmt="('t= ',e17.10,' phase_tot= ',e11.4,' phase_theta= ',e11.4)") &
!             & t, phase_tot, phase_theta
!     end if

    if (.not. (write_any .or. dump_any)) return

    if (debug) write(6,*) "loop_diagnostics: -2"

    call kperp (ntg_out, akperp)

    if (proc0 .and. write_any) then
       if (write_ascii) write (unit=out_unit, fmt=*) 'time=', t

       if (write_flux_line) then
          hflux_tot = 0.
          zflux_tot = 0.
          if (fphi > epsilon(0.0)) then
             if (write_ascii) then
                write (unit=out_unit, fmt=fluxes_str) &
                     't=', t, '<phi**2>=', phi2, 'heat_fluxes:', heat_fluxes(1:min(nspec,5)), &
                     'qflux_avg:', qflux_avg(1:min(nspec,5))/(t+code_dt)
                write (unit=out_unit, fmt=fluxes_str) &
                     't=', t, '<phi**2>=', phi2, 'part_fluxes:', part_fluxes(1:min(nspec,5)), &
                     'pflux_avg:', pflux_avg(1:min(nspec,5))/(t+code_dt)
                write (unit=out_unit, fmt=fluxes_str) &
                     't=', t, '<phi**2>=', phi2, 'mom_fluxes:', mom_fluxes(1:min(nspec,5)), &
                     'vflux_avg:', vflux_avg(1:min(nspec,5))/(t+code_dt)
                write (unit=out_unit, fmt=105) &
                     't=', t, '<phi**2>=', phi2, '<phip**2>=', phip2, &
                     'heating:', energy_exchange(1:min(nspec,5))
                
#ifdef LOWFLOW
!                   write (unit=out_unit, fmt="('t= ',e17.10,' <phi**2>= ',e11.4, &
!                        & ' lfmom fluxes: ', 5(1x,e11.4),' lfvflx1: ', 5(1x,e11.4))") &
!                        t, phi2, lfmom_fluxes(1:min(nspec,5)), vflux1_avg(1:min(nspec,5))
#endif
             end if
             hflux_tot = sum(heat_fluxes)
             vflux_tot = sum(mom_fluxes)
             zflux_tot = sum(part_fluxes*spec%z)
          end if
          if (fapar > epsilon(0.0)) then
             if (write_lorentzian .and. write_ascii) then
                wtmp_new = antenna_w()
                if (abs(real(wtmp_old)) > epsilon(0.) &
                     .and. abs(wtmp_new-wtmp_old) > epsilon(0.)) then
                   write (unit=out_unit,fmt='(a3,e17.10,a5,e17.10)') &
                        'w=', real(wtmp_new), 'amp=', sqrt(2.*apar2)
                end if
                wtmp_old = wtmp_new
             end if
             if (write_ascii) then
                write (unit=out_unit, fmt=103) &
                     't=', t, '<apar**2>=', apar2, &
                     'heat mluxes:', mheat_fluxes(1:min(nspec,5))
                write (unit=out_unit, fmt=103) &
                     't=', t, '<apar**2>=', apar2, &
                     'part mluxes:', mpart_fluxes(1:min(nspec,5))
             end if
             hflux_tot = hflux_tot + sum(mheat_fluxes)
             vflux_tot = vflux_tot + sum(mmom_fluxes)
             zflux_tot = zflux_tot + sum(mpart_fluxes*spec%z)
          end if
          if (fbpar > epsilon(0.0)) then
             if (write_ascii) then
                write (unit=out_unit, fmt=103) &
                     't=', t, '<bpar**2>=', bpar2, &
                     'heat bluxes:', bheat_fluxes(1:min(nspec,5))
                write (unit=out_unit, fmt=103) &
                     't=', t, '<bpar**2>=', bpar2, &
                     'part bluxes:', bpart_fluxes(1:min(nspec,5))
             end if
             hflux_tot = hflux_tot + sum(bheat_fluxes)
             vflux_tot = vflux_tot + sum(bmom_fluxes)
             zflux_tot = zflux_tot + sum(bpart_fluxes*spec%z)
          end if
          if (write_ascii) write (unit=out_unit,fmt=102) 't=', t, 'h_tot=', hflux_tot, 'z_tot=', zflux_tot

          if (write_nl_flux) then
             call nc_qflux (nout, qheat(:,:,:,1), qmheat(:,:,:,1), qbheat(:,:,:,1), &
                  heat_par, mheat_par, bheat_par, &
                  heat_perp, mheat_perp, bheat_perp, &
                  heat_fluxes, mheat_fluxes, bheat_fluxes, hflux_tot, &
                  energy_exchange, heath_fluxes)
                  ! Update the target array in parameter_scan_arrays
! below line gives out-of-bounds array for runs inside trinity
!                  scan_hflux(nout) = hflux_tot
                  scan_hflux(mod(nout-1,nstep/nwrite+1)+1) = hflux_tot
             call nc_vflux (nout, vflux, vmflux, vbflux, &
                  mom_fluxes, mmom_fluxes, bmom_fluxes, vflux_tot, &
                  vflux_par, vflux_perp, vflux0, vflux1, momh_fluxes)
                  ! Update the target array in parameter_scan_arrays
! below line gives out-of-bounds array for runs inside trinity
!                  scan_momflux(nout) = vflux_tot
                  scan_momflux(mod(nout-1,nstep/nwrite+1)+1) = vflux_tot
             call nc_pflux (nout, pflux, pmflux, pbflux, &
                  part_fluxes, mpart_fluxes, bpart_fluxes, zflux_tot, parth_fluxes)
          end if
          call nc_loop (nout, t, fluxfac, &
               phinew(igomega,:,:), phi2, phi2_by_mode, &
               aparnew(igomega,:,:), apar2, apar2_by_mode, &
               bparnew(igomega,:,:), bpar2, bpar2_by_mode, &
               omega, omegaavg, woutunits, write_omega)
! below line gives out-of-bounds array for runs inside trinity
!               scan_phi2_tot(nout) = phi2
               scan_phi2_tot(mod(nout-1,nstep/nwrite+1)+1) = phi2
       end if
       if (write_ascii) then
          do ik = 1, naky
             do it = 1, ntheta0
!                write (out_unit,*) "aky=",aky(ik)," theta0=",theta0(it,ik)

                if (write_line) then
                   write (out_unit,101) &
                        't=', t, 'aky=', aky(ik), 'akx=', akx(it), &
                        'om=', real( omega(it,ik)*woutunits(ik)), &
                        aimag(omega(it,ik)*woutunits(ik)), &
                        'omav=', real( omegaavg(it,ik)*woutunits(ik)), &
                        aimag(omegaavg(it,ik)*woutunits(ik)), &
                        'phtot=', phitot(it,ik)
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

    end if

    call broadcast(scan_hflux)
    call broadcast(scan_momflux)
    call broadcast(scan_phi2_tot)

    call broadcast (write_symmetry)
    if (write_symmetry) then
       allocate (pflx_sym(-ntgrid:ntgrid,-nvgrid:nvgrid,nspec))
       allocate (vflx_sym(-ntgrid:ntgrid,-nvgrid:nvgrid,nspec))
       allocate (qflx_sym(-ntgrid:ntgrid,-nvgrid:nvgrid,nspec))
       allocate (pflxh_sym(-ntgrid:ntgrid,-nvgrid:nvgrid,nspec)) ; pflxh_sym = 0.
       allocate (vflxh_sym(-ntgrid:ntgrid,-nvgrid:nvgrid,nspec)) ; vflxh_sym = 0.
       allocate (qflxh_sym(-ntgrid:ntgrid,-nvgrid:nvgrid,nspec)) ; qflxh_sym = 0.

!       call flux_vs_theta_vs_vpa (ghnew, phinew, pflx_sym, vflx_sym, qflx_sym)
!       call flux_vs_theta_vs_vpa (gnew, phihnew, pflxh_sym, vflxh_sym, qflxh_sym)
!       pflxh_sym = pflx_sym + pflxh_sym
!       vflxh_sym = vflx_sym + vflxh_sym
!       qflxh_sym = qflx_sym + qflxh_sym
       
       call flux_vs_theta_vs_vpa (gnew, phinew, pflx_sym, vflx_sym, qflx_sym)
       if (proc0) call nc_loop_sym (nout, pflx_sym, vflx_sym, qflx_sym, &
            pflxh_sym, vflxh_sym, qflxh_sym)
       deallocate (pflx_sym, vflx_sym, qflx_sym, pflxh_sym, vflxh_sym, qflxh_sym)
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

!          write (out_unit, "('t,u0 ',2(1x,e13.6))") t, aimag(upar00 (1,1))

          call nc_loop_moments (nout, ntot2, ntot2_by_mode, ntot20, &
               ntot20_by_mode, phi00, ntot00, density00, upar00, &
               tpar00, tperp00, tpar2_by_mode, tperp2_by_mode)

       end if
    end if

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
       write (filename, "('dump.fields.t=',e13.6)") t
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
    end if

    deallocate (omega, omegaavg)

    t_old = t

101 format (a3,e17.9,a5,1pe12.4,a5,1pe12.4,a4,1p2e12.4,a6,1p2e12.4,a7,1pe12.4)
102 format (a3,e17.10,a7,e11.4,a7,e11.4)
103 format (a3,e17.10,a11,e11.4,a13,5(1x,e11.4))
105 format (a3,e17.10,a10,e11.4,a11,e11.4,a9,5(1x,e11.4))
106 format (a4,f11.4,a4,f11.4,a4,2f10.3,a6,2f10.3,a7,e9.2)
107 format (a3,e17.10,a11,e11.4,a13,5(1x,e11.4))
    
if (debug) write(6,*) "loop_diagnostics: done"
  end subroutine loop_diagnostics

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
       if (aky(ik) < epsilon(0.)) fac = 1.0
       do it = 1, ntheta0
          favg = favg + f(it, ik) * fac
       end do
    end do

  end subroutine get_volume_average

  subroutine get_surf_average (f, favg)
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
       if (aky(ik) < epsilon(0.)) fac = 1.0
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
    real, dimension (nensembles,nspec) :: pflx_global, qflx_global, vflx_global
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
  end subroutine ensemble_average

!   subroutine get_fldline_avg_r (fld_in, fld_out)

!     use theta_grid, only: delthet, jacob

!     implicit none

!     real, dimension (:), allocatable :: dl_over_b

!     real, dimension (-ntg_out:), intent (in) :: fld_in
!     real, intent (out) :: fld_out

!     allocate (dl_over_b(-ntg_out:ntg_out))

!     dl_over_b = delthet(-ntg_out:ntg_out)*jacob(-ntg_out:ntg_out)
!     dl_over_b = dl_over_b / sum(dl_over_b)

!     fld_out = sum(fld_in*dl_over_b)

!     deallocate (dl_over_b)

!   end subroutine get_fldline_avg_r

!   subroutine get_fldline_avg_c (fld_in, fld_out)

!     use theta_grid, only: delthet, jacob

!     implicit none

!     real, dimension (:), allocatable :: dl_over_b

!     complex, dimension (-ntg_out:), intent (in) :: fld_in
!     complex, intent (out) :: fld_out

!     allocate (dl_over_b(-ntg_out:ntg_out))

!     dl_over_b = delthet(-ntg_out:ntg_out)*jacob(-ntg_out:ntg_out)
!     dl_over_b = dl_over_b / sum(dl_over_b)

!     fld_out = sum(fld_in*dl_over_b)

!     deallocate (dl_over_b)

!   end subroutine get_fldline_avg_c

!   subroutine get_cross_phase (phase_tot, phase_theta)

! ! <doc> This is a highly simplified synthetic diagnostic which 
! ! calculates the cross phase between the electron density and the 
! ! perpendicular electron temperature for comparisons with DIII-D.  
! ! Returns the value of the cross-phase at the outboard midplane and 
! ! integrated over all v and x. We can generalize this routine to 
! ! other fields at some point, but for now this is just a skeleton for 
! ! a more realistic synthetic diagnostic. </doc>

!     use species, only: nspec, spec, electron_species
!     use kt_grids, only: ntheta0, naky
!     use theta_grid, only: ntgrid
!     use dist_fn_arrays, only: gnew
!     use dist_fn, only: getemoms
!     use mp, only: proc0

!     implicit none
!     real, intent (out) :: phase_tot, phase_theta
!     complex, dimension (:,:,:,:), allocatable :: ntot, tperp
!     complex, dimension (ntheta0, naky) :: nTp_by_mode
!     complex :: nTp
!     real, dimension (ntheta0, naky) :: n2_by_mode, T2_by_mode
!     real :: n2, T2

!     integer :: it, ik, is, isgn, ig
!     integer :: iglo

!     allocate ( ntot(-ntgrid:ntgrid,ntheta0,naky,nspec))
!     allocate (tperp(-ntgrid:ntgrid,ntheta0,naky,nspec))

!     call getemoms (ntot, tperp)

!     do is = 1,nspec
!        if (spec(is)%type == electron_species) then
!           ! get cross_phase at outboard midplane
!           call get_vol_int (ntot(0,:,:,is), tperp(0,:,:,is), nTp, nTp_by_mode)
! !          call get_vol_average (ntot(0,:,:,is), ntot(0,:,:,is), n2, n2_by_mode)
! !          call get_vol_average (tperp(0,:,:,is), tperp(0,:,:,is), T2, T2_by_mode)
!           phase_theta = atan2(aimag(nTp),real(nTp))!/sqrt(n2*T2)
!           ! get integrated cross_phase 
!           call get_vol_int (ntot(:,:,:,is), tperp(:,:,:,is), nTp, nTp_by_mode)
! !          call get_vol_average (ntot(:,:,:,is), ntot(:,:,:,is), n2, n2_by_mode)
! !          call get_vol_average (tperp(:,:,:,is), tperp(:,:,:,is), T2, T2_by_mode)
!           phase_tot = atan2(aimag(nTp),real(nTp))!/sqrt(n2*T2)
!        end if
!     end do

!     deallocate (ntot, tperp)

!   end subroutine get_cross_phase

!   subroutine get_vol_int_all (a, b, axb, axb_by_mode)
!     use theta_grid, only: ntgrid, delthet, jacob
!     use kt_grids, only: naky, ntheta0
!     implicit none
!     complex, dimension (-ntgrid:,:,:), intent (in) :: a, b
!     complex, intent (out) :: axb
!     complex, dimension (:,:), intent (out) :: axb_by_mode

!     integer :: ik, it
!     integer :: ng
!     real, dimension (-ntg_out:ntg_out) :: wgt
!     real :: anorm

!     ng = ntg_out
!     wgt = delthet(-ng:ng)*jacob(-ng:ng)
!     anorm = sum(wgt)

!     do ik = 1, naky
!        do it = 1, ntheta0
!           axb_by_mode(it,ik) &
!                = sum((conjg(a(-ng:ng,it,ik))*b(-ng:ng,it,ik))*wgt)/anorm
!        end do
!     end do

!     call get_volume_int (axb_by_mode, axb)
!   end subroutine get_vol_int_all

!   subroutine get_vol_int_one (a, b, axb, axb_by_mode)

!     implicit none
!     complex, dimension (:,:), intent (in) :: a, b
!     complex, intent (out) :: axb
!     complex, dimension (:,:), intent (out) :: axb_by_mode

!     axb_by_mode = conjg(a)*b
!     call get_volume_int (axb_by_mode, axb)

!   end subroutine get_vol_int_one

!   subroutine get_volume_int (f, favg)
!     use mp, only: iproc
!     use kt_grids, only: naky, ntheta0, aky
!     implicit none
!     complex, dimension (:,:), intent (in) :: f
!     complex, intent (out) :: favg
!     real :: fac
!     integer :: ik, it

! ! ky=0 modes have correct amplitudes; rest must be scaled
! ! note contrast with scaling factors in FFT routines.

!     favg = 0.
!     do ik = 1, naky
!        fac = 0.5
!        if (aky(ik) == 0.) fac = 1.0
!        do it = 1, ntheta0
!           favg = favg + f(it, ik) * fac
!        end do
!     end do

!   end subroutine get_volume_int

!  subroutine autocorrelation (cfnc,phi2extend,cfnc_2pi)
  subroutine correlation_extend (cfnc,phi2extend)

    use fields_arrays, only: phinew
    use theta_grid, only: ntgrid, jacob, delthet
    use kt_grids, only: ntheta0, naky, jtwist_out

    implicit none

    real, dimension (:,:,:), intent (out) :: phi2extend
    complex, dimension (:,:,:), intent (out) :: cfnc

    integer :: ig, it, ik, im, igmod
    integer :: itshift, nconnect, offset

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
          ! number of connections between kx values
          nconnect = (ntheta0-it)/itshift
          ! shift of theta index to account for fact that not all ky values
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
    
    phi2extend = real(phiextend*conjg(phiextend))

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
