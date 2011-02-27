! what about boundary condition contributions?  In the presence 
! of magnetic shear, there are not always zero amplitudes at the ends
! of the supercells.

module gs2_diagnostics
  use gs2_heating, only: heating_diagnostics,dens_vel_diagnostics

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

  interface get_vol_int
     module procedure get_vol_int_one, get_vol_int_all
  end interface

  interface get_fldline_avg
     module procedure get_fldline_avg_r, get_fldline_avg_c
  end interface

!CMR, 17/11/2009:   read_parameters now public so ingen can USE instead of copy
!

! Why are these variables public?  This is not good.
  real,public :: omegatol, omegatinst
  logical,public :: print_line, print_old_units, print_flux_line
  logical,public :: print_summary, write_line, write_flux_line, write_phi
  logical,public :: write_apar, write_bpar, write_aperp
  logical,public :: write_omega, write_omavg, write_ascii, write_lamavg
  logical,public :: write_eavg, write_qheat, write_pflux, write_vflux
  logical,public :: write_tavg, write_qmheat, write_pmflux, write_vmflux
  logical,public :: write_gs, write_lpoly, write_qbheat
  logical,public :: write_pbflux, write_vbflux, write_g, write_gg, write_gyx
  logical,public ::write_dmix, write_kperpnorm, write_phitot, write_epartot
  logical,public :: write_eigenfunc, write_fields, write_final_fields, write_final_antot
  logical,public :: write_final_moments, write_avg_moments, write_stress, write_parity
  logical,public :: write_full_moments_notgc, write_cross_phase = .false.
  logical,public :: write_fcheck, write_final_epar, write_kpar
  logical,public :: write_vortcheck, write_fieldcheck
  logical,public :: write_fieldline_avg_phi, write_hrate, write_lorentzian
  logical,public :: write_neoclassical_flux, write_nl_flux, write_Epolar
  logical,public :: write_verr, write_cerr, write_max_verr
  logical,public :: exit_when_converged
  logical,public :: dump_neoclassical_flux, dump_check1, dump_check2
  logical,public :: dump_fields_periodically, make_movie
  logical,public :: dump_final_xfields
  logical,public :: use_shmem_for_xfields
  logical,public :: save_for_restart
  logical,public :: save_distfn !<DD> Added for saving distribution function
  logical,public :: test_conserve
  logical,public :: write_symmetry
  integer,public :: nwrite, igomega, nmovie
  integer,public :: navg, nsave
  integer,public :: nperiod_output

!>GGH
  logical, parameter :: write_density_velocity=.false.
  logical :: write_jext=.false.
!<GGH

  ! internal
  logical :: write_any, write_any_fluxes, dump_any
  logical, private :: initialized = .false.

  namelist /gs2_diagnostics_knobs/ print_line, print_old_units, print_flux_line, &
         write_line, write_flux_line, write_phi, write_apar, write_bpar, write_aperp, &
         write_omega, write_omavg, write_ascii, write_kpar, write_lamavg, &
         write_qheat, write_pflux, write_vflux, write_eavg, write_gs, write_gyx, &
         write_qmheat, write_pmflux, write_vmflux, write_tavg, write_g, write_gg, &
         write_qbheat, write_pbflux, write_vbflux, write_hrate, write_lpoly, &
         write_dmix, write_kperpnorm, write_phitot, write_epartot, &
         write_eigenfunc, write_fields, write_final_fields, write_final_antot, &
         write_fcheck, write_final_epar, write_final_moments, write_cerr, &
         write_vortcheck, write_fieldcheck, write_Epolar, write_verr, write_max_verr, &
         write_fieldline_avg_phi, write_neoclassical_flux, write_nl_flux, &
         nwrite, nmovie, nsave, navg, omegatol, omegatinst, igomega, write_lorentzian, &
         exit_when_converged, write_avg_moments, write_stress, &
         write_full_moments_notgc, write_cross_phase, &
         dump_neoclassical_flux, dump_check1, dump_check2, &
         dump_fields_periodically, make_movie, &
         dump_final_xfields, use_shmem_for_xfields, &
         nperiod_output, test_conserve, &
         save_for_restart, write_parity, write_symmetry, save_distfn !<DD> Added for saving distribution function

  integer :: out_unit, kp_unit, heat_unit, polar_raw_unit, polar_avg_unit, heat_unit2, lpc_unit
  integer :: dv_unit, jext_unit   !GGH Additions
  integer :: phase_unit
  integer :: dump_neoclassical_flux_unit, dump_check1_unit, dump_check2_unit
  integer :: res_unit, res_unit2, parity_unit

  complex, dimension (:,:,:), allocatable :: omegahist
  ! (navg,ntheta0,naky)
  type (heating_diagnostics), save :: h
  type (heating_diagnostics), dimension(:), save, allocatable :: h_hist
  type (heating_diagnostics), dimension(:,:), save, allocatable :: hk
  type (heating_diagnostics), dimension(:,:,:), save, allocatable :: hk_hist
  !GGH Density/velocity pertubration diagnostics
  type (dens_vel_diagnostics), dimension(:), allocatable :: dv_hist
  type (dens_vel_diagnostics), dimension(:,:,:), allocatable :: dvk_hist
  !GGH J_external
  real, dimension(:,:,:), allocatable ::  j_ext_hist
  !Polar spectrum variables
  integer :: nbx                                   !Number of polar bins
  real, dimension(:), allocatable :: kpbin,ebincorr
  real, dimension(:), allocatable :: numavg,kpavg
  real, dimension(:,:), pointer :: ebinarray => null ()
  real, dimension(:,:), pointer :: eavgarray => null ()
  real, dimension(:,:), allocatable :: etmp
  integer, dimension(:,:), allocatable :: polar_index
  integer, dimension(:), allocatable :: polar_avg_index
  integer, parameter :: iefluc = 1        !Total fluctuating energy
  integer, parameter :: ieapar = 2        !Energy in A_parallel
  integer, parameter :: iebpar = 3        !Energy in B_parallel
  integer, parameter :: iephis2= 4        !Energy in q_s^ n_s Phi^2/T_s
  integer, parameter :: iehs2   = 5       !Energy in hs^2
  integer, parameter :: iedelfs2= 6       !Energy in del f_s^2

  real, dimension (:,:,:,:), allocatable ::  qheat !,  qheat_par,  qheat_perp
  real, dimension (:,:,:,:), allocatable :: qmheat !, qmheat_par, qmheat_perp
  real, dimension (:,:,:,:), allocatable :: qbheat !, qbheat_par, qbheat_perp
  ! (ntheta0,naky,nspec,3)

  real, dimension (:,:,:), allocatable ::  pflux,  vflux, vflux_par, vflux_perp
  real, dimension (:,:,:), allocatable :: vflux0, vflux1  ! low flow correction to turbulent momentum flux
  real, dimension (:,:,:), allocatable :: pmflux, vmflux
  real, dimension (:,:,:), allocatable :: pbflux, vbflux
  real, dimension (:,:),   allocatable :: theta_pflux, theta_pmflux, theta_pbflux
  real, dimension (:,:),   allocatable :: theta_vflux, theta_vmflux, theta_vbflux
  real, dimension (:,:),   allocatable :: theta_vflux_par, theta_vflux_perp
  real, dimension (:,:,:), allocatable :: theta_qflux, theta_qmflux, theta_qbflux

  ! (ntheta0,naky,nspec)

  real :: start_time = 0.0
  real, dimension (:), allocatable :: pflux_avg, qflux_avg, heat_avg, vflux_avg

  integer :: ntg_out
  integer :: nout = 1
  integer :: nout_movie = 1
  complex :: wtmp_old = 0.
  logical :: exist

contains

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
       write (unit, fmt="(' nwrite = ',i6)") nwrite
       write (unit, fmt="(' nsave = ',i6)") nsave
       write (unit, fmt="(' navg = ',i6)") navg
       write (unit, fmt="(' omegatol = ',e16.10)") omegatol
       write (unit, fmt="(' omegatinst = ',e16.10)") omegatinst
! should be legal -- not checked yet
       if (igomega /= 0) write (unit, fmt="(' igomega = ',i6)") igomega  
!       if (nperiod_output /= nperiod) &
!            write (unit, fmt="(' nperiod_output = ',i3)") nperiod_output
       
       write (unit, fmt="(' print_old_units = ',L1)") print_old_units
       if (write_ascii) then
          write (unit, fmt="(' write_ascii = ',L1)") write_ascii
          write (unit, fmt="(' write_omega = ',L1)") write_omega
          write (unit, fmt="(' write_omavg = ',L1)") write_omavg
          write (unit, fmt="(' write_dmix = ',L1)") write_dmix
          write (unit, fmt="(' write_kperpnorm = ',L1)") write_kperpnorm
       end if
       if (write_Epolar) &
            write (unit, fmt="(' write_Epolar = ',L1)") write_Epolar
       write (unit, fmt="(' write_hrate = ',L1)") write_hrate
       write (unit, fmt="(' write_lorentzian = ',L1)") write_lorentzian
       write (unit, fmt="(' write_eigenfunc = ',L1)") write_eigenfunc
       write (unit, fmt="(' write_final_fields = ',L1)") write_final_fields
       write (unit, fmt="(' write_final_epar = ',L1)") write_final_epar
       write (unit, fmt="(' write_final_moments = ',L1)") write_final_moments
       write (unit, fmt="(' write_final_antot = ',L1)") write_final_antot
       write (unit, fmt="(' write_tavg = ',L1)") write_tavg
       write (unit, fmt="(' write_lamavg = ',L1)") write_lamavg
       write (unit, fmt="(' write_eavg = ',L1)") write_eavg
       if (write_fcheck) write (unit, fmt="(' write_fcheck = ',L1)") write_fcheck
       if (write_vortcheck) write (unit, fmt="(' write_vortcheck = ',L1)") write_vortcheck
       if (write_fieldcheck) write (unit, fmt="(' write_fieldcheck = ',L1)") write_fieldcheck
       if (write_fieldline_avg_phi) &
            write (unit, fmt="(' write_fieldline_avg_phi = ',L1)") write_fieldline_avg_phi
       if (write_neoclassical_flux) &
            write (unit, fmt="(' write_neoclassical_flux = ',L1)") write_neoclassical_flux
       write (unit, fmt="(' write_nl_flux = ',L1)") write_nl_flux
       write (unit, fmt="(' exit_when_converged = ',L1)") exit_when_converged
       if (write_avg_moments) write (unit, fmt="(' write_avg_moments = ',L1)") write_avg_moments
       if (dump_neoclassical_flux) &
            write (unit, fmt="(' dump_neoclassical_flux = ',L1)") dump_neoclassical_flux
       if (dump_check1) write (unit, fmt="(' dump_check1 = ',L1)") dump_check1
       if (dump_check2) write (unit, fmt="(' dump_check2 = ',L1)") dump_check2
       if (dump_fields_periodically) &
            write (unit, fmt="(' dump_fields_periodically = ',L1)") dump_fields_periodically
       if (make_movie) &
            write (unit, fmt="(' make_movie = ',L1)") make_movie
       if (dump_final_xfields) &
            write (unit, fmt="(' dump_final_xfields = ',L1)") dump_final_xfields

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

    if (print_old_units) then
       write (report_unit, fmt="('print_old_units = T:       Frequencies on screen in 1/omega_* units, omega_*=(cT/eB)*ky/L_ref.')")
    end if

    if (.not. write_phi) then
       write (report_unit, fmt="('write_phi = F:             Ignored.')")
    end if

    if (.not. write_apar) then
       write (report_unit, fmt="('write_apar = F:            Ignored.')")
    end if

    if (.not. write_phi) then
       write (report_unit, fmt="('write_aperp = F:           Ignored.')")
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

    if (write_lamavg) then
       write (report_unit, fmt="('write_lamavg = T:          Write particle flux vs. lambda to ',a)") trim(run_name)//'.lam'
       write (report_unit, fmt="('write_lamavg = T:          Write energy flux vs. lambda to ',a)") trim(run_name)//'.lame'
    end if

    if (write_tavg) then
       write (report_unit, fmt="('write_tavg = T:            Write particle flux vs. theta to ',a)") trim(run_name)//'.theta'
       write (report_unit, fmt="('write_tavg = T:          Write energy flux vs. thetaa to ',a)") trim(run_name)//'.thetae'
    end if

    if (write_eavg) then
       write (report_unit, &
         & fmt="('write_eavg = T:            Write particle flux vs. energy to ',a)") trim(run_name)//'.energy'
       write (report_unit, &
         & fmt="('write_eavg = T:            Write energy flux vs. energy to ',a)") trim(run_name)//'.energye'
    end if

    if (write_dmix) then
       if (write_ascii) then
          write (report_unit, fmt="('write_dmix = T:            Write D_ML ',a)") trim(run_name)//'.out'
       else
          write (report_unit, fmt="('write_dmix = T:            Ignored if write_ascii = F')")
       end if
    end if

    if (write_kperpnorm) then
       write (report_unit, fmt="('write_kperpnorm = T:       Ignored.')")
    end if

    if (write_phitot) then
       write (report_unit, fmt="('write_phitot = T:          Ignored.')")
    end if
       
    if (write_epartot) then
       write (report_unit, fmt="('write_epartot = T:         Ignored.')")
    end if

    if (write_fieldline_avg_phi) then
       write (report_unit, fmt="('write_fieldline_avg_phi = T: Ignored.')")
       write (report_unit, fmt="('    Perhaps you want write_avg_moments = T')")
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

    if (write_fcheck) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('write_fcheck = T:               Turns on obscure diagnostics.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (write_vortcheck) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('write_vortcheck = T:              Turns on obscure diagnostics.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (write_fieldcheck) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('write_fieldcheck = T:              Turns on obscure diagnostics.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (write_neoclassical_flux) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('write_neoclassical_flux = T:               Turns on neoclassical flux calc, &
           & but result not written.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('    Perhaps you want dump_neoclassical_flux = T.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (write_nl_flux) then
       if (write_ascii) then
          write (report_unit, fmt="('write_nl_flux = T:         Phi**2(kx, ky) written to ',a)") trim(run_name)//'.out'
       end if
    else
       write (report_unit, fmt="('write_nl_flux = F:         Phi**2(kx, ky) NOT written to ',a)") trim(run_name)//'.out'
    end if

    if (dump_neoclassical_flux) then
       write (report_unit, fmt="('dump_neoclassical_flux = T: Neoclassical fluxes written to ',a)") 'dump.neoflux'
       write (report_unit, fmt="('This option requires an expert user.')") 
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
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

    if (dump_final_xfields) then
       write (report_unit, fmt="('dump_final_xfields is not longer maintained')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
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

    if (write_pflux) write (report_unit, fmt="('write_pflux = T:           Ignored.')")
    if (write_vflux) write (report_unit, fmt="('write_vflux = T:           Ignored.')")
    if (write_qheat) write (report_unit, fmt="('write_qheat = T:           Ignored.')")
    if (write_pmflux) write (report_unit, fmt="('write_pmflux = T:          Ignored.')")
    if (write_vmflux) write (report_unit, fmt="('write_vmflux = T:          Ignored.')")
    if (write_qmheat) write (report_unit, fmt="('write_qmheat = T:          Ignored.')")
    if (write_pbflux) write (report_unit, fmt="('write_pbflux = T:          Ignored.')")
    if (write_vbflux) write (report_unit, fmt="('write_vbflux = T:          Ignored.')")
    if (write_qbheat) write (report_unit, fmt="('write_qbheat = T:          Ignored.')")    

  end subroutine check_gs2_diagnostics


  subroutine init_gs2_diagnostics (list, nstep)
   !<doc> Define NetCDF vars, call real_init, which calls read_parameters; broadcast all the different write flags. </doc>

    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids, ntheta0, naky
    use run_parameters, only: init_run_parameters
    use species, only: init_species, nspec
    use dist_fn, only: init_dist_fn
    use init_g, only: init_init_g
!    use gs2_flux, only: init_gs2_flux
    use gs2_io, only: init_gs2_io
    use gs2_heating, only: init_htype,init_dvtype
    use collisions, only: collision_model_switch, init_lorentz_error
    use mp, only: broadcast, proc0
    use le_grids, only: init_weights

    implicit none
    logical, intent (in) :: list
    integer, intent (in) :: nstep
    integer :: nmovie_tot

    if (initialized) return
    initialized = .true.

    call init_theta_grid
    call init_kt_grids
    call init_run_parameters
    call init_species
    call init_init_g
    call init_dist_fn
!    call init_gs2_flux

    call real_init (list)
    call broadcast (navg)
    call broadcast (nwrite)
    call broadcast (nmovie)
    call broadcast (nsave)
    call broadcast (write_any)
    call broadcast (write_any_fluxes)
    call broadcast (write_cross_phase)
    call broadcast (write_neoclassical_flux)
    call broadcast (write_nl_flux)
    call broadcast (write_omega)
    call broadcast (write_fieldline_avg_phi)
    call broadcast (write_fcheck)
    call broadcast (dump_any)
    call broadcast (dump_neoclassical_flux)
    call broadcast (dump_check1)
    call broadcast (dump_check2)
    call broadcast (write_fields)
    call broadcast (dump_fields_periodically)
    call broadcast (make_movie)
    call broadcast (dump_final_xfields)
    call broadcast (use_shmem_for_xfields)
    call broadcast (save_for_restart)
    call broadcast (save_distfn) !<DD> Added for saving distribution function
    call broadcast (write_gs)
    call broadcast (write_g)
    call broadcast (write_gyx)
    call broadcast (write_gg)
    call broadcast (write_stress)
    call broadcast (write_final_antot)
    call broadcast (write_verr)
    call broadcast (write_max_verr)
    call broadcast (write_lpoly)
    call broadcast (write_cerr)

    call broadcast (write_vortcheck)
    call broadcast (write_fieldcheck)
    call broadcast (ntg_out)
    call broadcast (write_lamavg)
    call broadcast (write_eavg)
    call broadcast (write_tavg)
    call broadcast (write_hrate)
    call broadcast (write_Epolar)
    call broadcast (write_lorentzian)
    call broadcast (write_eigenfunc)

    call broadcast (write_full_moments_notgc)

    nmovie_tot = nstep/nmovie

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
       
!GGH Allocate density and velocity perturbation diagnostic structures
    if (write_density_velocity) then
       allocate (dv_hist(0:navg-1))
       call init_dvtype (dv_hist,  nspec)

       allocate (dvk_hist(ntheta0,naky,0:navg-1))
       call init_dvtype (dvk_hist, nspec)
    end if
       
!GGH Allocate density and velocity perturbation diagnostic structures
    if (write_jext) allocate (j_ext_hist(ntheta0, naky,0:navg-1)) 

!Initialize polar spectrum diagnostic
    if (write_Epolar .and. proc0) call init_polar_spectrum

!       allocate (hratehist(nspec, 7, 0:navg-1));  hratehist = 0.0
!       allocate (hkratehist(ntheta0, naky, nspec, 7, 0:navg-1));  hkratehist = 0.0

    call init_gs2_io (write_nl_flux, write_omega, write_stress, &
         write_fieldline_avg_phi, write_hrate, write_final_antot, &
         write_eigenfunc, make_movie, nmovie_tot, write_verr, &
         write_fields, write_full_moments_notgc, write_symmetry)
    
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
    use run_parameters, only: fapar
    use file_utils, only: open_output_file, get_unused_unit
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, aky_out, akx_out
    use gs2_layouts, only: yxf_lo
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

       !GGH Density and velocity perturbations
       if (write_density_velocity .and. write_ascii) then
          call open_output_file (dv_unit, ".dv")
       end if

       !GGH J_external, only if A_parallel is being calculated.
       if (write_jext .and. fapar > epsilon(0.0)) then
          if (write_ascii) then
             call open_output_file (jext_unit, ".jext")
          end if
       else
          write_jext = .false.
       end if

       if (write_Epolar .and. write_ascii) then
          call open_output_file (polar_raw_unit, ".kspec_raw")
          call open_output_file (polar_avg_unit, ".kspec_avg")
       end if

       if (dump_neoclassical_flux) then
          call get_unused_unit (dump_neoclassical_flux_unit)
          open (unit=dump_neoclassical_flux_unit, file="dump.neoflux", &
               status="unknown")
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

    !<doc> Allocate arrays for storing the various fluxes which the diagnostics will output </doc>
    allocate (pflux (ntheta0,naky,nspec)) ; pflux = 0.
    allocate (qheat (ntheta0,naky,nspec,3)) ; qheat = 0.
!       allocate (qheat_par  (ntheta0,naky,nspec))
!       allocate (qheat_perp (ntheta0,naky,nspec))
    allocate (vflux (ntheta0,naky,nspec)) ; vflux = 0.
    allocate (vflux_par (ntheta0,naky,nspec)) ; vflux_par = 0.
    allocate (vflux_perp (ntheta0,naky,nspec)) ; vflux_perp = 0.
    allocate (vflux0 (ntheta0,naky,nspec)) ; vflux0 = 0.
    allocate (vflux1 (ntheta0,naky,nspec)) ; vflux1 = 0.
    allocate (pmflux(ntheta0,naky,nspec)) ; pmflux = 0.
    allocate (qmheat(ntheta0,naky,nspec,3)) ; qmheat = 0.
!       allocate (qmheat_par (ntheta0,naky,nspec))
!       allocate (qmheat_perp(ntheta0,naky,nspec))
    allocate (vmflux(ntheta0,naky,nspec)) ; vmflux = 0.
    allocate (pbflux(ntheta0,naky,nspec)) ; pbflux = 0.
    allocate (qbheat(ntheta0,naky,nspec,3)) ; qbheat = 0.
!       allocate (qbheat_par (ntheta0,naky,nspec))
!       allocate (qbheat_perp(ntheta0,naky,nspec))
    allocate (vbflux(ntheta0,naky,nspec)) ; vbflux = 0.
       
    allocate (theta_pflux (-ntgrid:ntgrid, nspec)) ; theta_pflux = 0.
    allocate (theta_vflux (-ntgrid:ntgrid, nspec)) ; theta_vflux = 0.
    allocate (theta_vflux_par (-ntgrid:ntgrid, nspec)) ; theta_vflux_par = 0.
    allocate (theta_vflux_perp (-ntgrid:ntgrid, nspec)) ; theta_vflux_perp = 0.
    allocate (theta_qflux (-ntgrid:ntgrid, nspec, 3)) ; theta_qflux = 0.
    
    allocate (theta_pmflux (-ntgrid:ntgrid, nspec)) ; theta_pmflux = 0.
    allocate (theta_vmflux (-ntgrid:ntgrid, nspec)) ; theta_vmflux = 0.
    allocate (theta_qmflux (-ntgrid:ntgrid, nspec, 3)) ; theta_qmflux = 0.
    
    allocate (theta_pbflux (-ntgrid:ntgrid, nspec)) ; theta_pbflux = 0.
    allocate (theta_vbflux (-ntgrid:ntgrid, nspec)) ; theta_vbflux = 0.
    allocate (theta_qbflux (-ntgrid:ntgrid, nspec, 3)) ; theta_qbflux = 0.

  end subroutine real_init

  subroutine read_parameters (list)
!CMR, 17/11/2009:   namelist gs2_diagnostics_knobs made public
!                   so that ingen can just USE it instead of copying it!
!
    use file_utils, only: input_unit, input_unit_exist
    use theta_grid, only: nperiod, ntheta
!    use gs2_flux, only: gs2_flux_adjust
    use dist_fn, only: nperiod_guard
    use kt_grids, only: box, nx, ny
    use mp, only: proc0
    implicit none
    integer :: in_file
    logical, intent (in) :: list

    !<doc> Set defaults for the gs2_diagnostics_knobs</doc>		
    if (proc0) then
       print_line = .true.
       print_old_units = .false.
       print_flux_line = .false.
       write_line = .true.
       write_flux_line = .true.
       write_kpar = .false.
       write_hrate = .false.
       write_Epolar = .false.
       write_gs = .false.
       write_g = .false.
       write_gyx = .false.
       write_lpoly = .false.
       write_gg = .false.
       write_lorentzian = .false.
       write_phi = .true.
       write_apar = .true.
       write_aperp = .true.
       write_bpar = .true.
       write_omega = .false.
       write_ascii = .true.
       write_eavg = .false.
       write_tavg = .false.
       write_lamavg = .false.
       write_omavg = .false.
       write_dmix = .false.
       write_kperpnorm = .false.
       write_phitot = .true.
       write_epartot = .false.
       write_fieldline_avg_phi = .false.
       write_neoclassical_flux = .false.
       write_nl_flux = .false.
       write_eigenfunc = .false.
       write_final_moments = .false.
       write_stress = .false.
       write_avg_moments = .false.
       write_parity = .false.
       write_symmetry = .false.
       write_fields = .false.
       write_full_moments_notgc = .false.
       write_final_fields = .false.
       write_final_antot = .false.
       write_final_epar = .false.
       write_fcheck = .false.
       write_vortcheck = .false.
       write_fieldcheck = .false.
       write_verr = .false.
       write_max_verr = .false.
       write_cerr = .false.
       test_conserve = .false.
       nwrite = 100
       nmovie = 1000
       navg = 100
       nsave = -1
       omegatol = 1e-3
       omegatinst = 1.0
       igomega = 0
       exit_when_converged = .true.
       dump_neoclassical_flux = .false.
       dump_check1 = .false.
       dump_check2 = .false.
       dump_fields_periodically = .false.
       make_movie = .false.
       dump_final_xfields = .false.
       use_shmem_for_xfields = .true.
       nperiod_output = nperiod - nperiod_guard
       save_for_restart = .false.
       save_distfn = .false. !<DD> Added for saving distribution function
       in_file = input_unit_exist ("gs2_diagnostics_knobs", exist)

	!<doc> Read in parameters from the namelist gs2_diagnostics_knobs, if the namelist exists </doc>

!       if (exist) read (unit=input_unit("gs2_diagnostics_knobs"), nml=gs2_diagnostics_knobs)
       if (exist) read (unit=in_file, nml=gs2_diagnostics_knobs)

! If either the old or new variable was set to .false., choose false.
       write_bpar = write_bpar .and. write_aperp

       print_summary = (list .and. (print_line .or. print_flux_line)) 

       if (list) then
          print_line = .false.
          print_flux_line = .false.
       end if

       if (.not. save_for_restart) nsave = -1
! changed temporarily for testing -- MAB
!       write_avg_moments = write_avg_moments .and. box
       write_avg_moments = write_avg_moments
       write_stress = write_stress .and. box

! Only calculate polar integrals in box layout
       write_Epolar = write_Epolar .and. box

! Disable polar integrals if nx /= ny
       if (nx /= ny) write_Epolar = .false.

!       write_avg_moments = write_avg_moments .or. write_neoclassical_flux

       write_any = write_line .or. write_phi       .or. write_apar &
            .or. write_bpar  .or. write_omega     .or. write_omavg &
            .or. write_dmix   .or. write_kperpnorm .or. write_phitot &
            .or. write_fieldline_avg_phi           .or. write_neoclassical_flux &
            .or. write_flux_line                   .or. write_nl_flux .or. write_Epolar &
            .or. write_kpar   .or. write_hrate     .or. write_lorentzian  .or. write_gs
       write_any_fluxes =  write_flux_line .or. print_flux_line .or. write_nl_flux !&
!            .or. gs2_flux_adjust
       dump_any = dump_neoclassical_flux .or. dump_check1  .or. dump_fields_periodically &
            .or. dump_check2 .or. make_movie .or. print_summary &
            .or. write_full_moments_notgc

       nperiod_output = min(nperiod,nperiod_output)
       ntg_out = ntheta/2 + (nperiod_output-1)*ntheta
    end if 
 end subroutine read_parameters

  subroutine finish_gs2_diagnostics (istep)
    use file_utils, only: open_output_file, close_output_file, get_unused_unit
    use mp, only: proc0, broadcast, nproc, iproc, sum_reduce
    use species, only: nspec, spec
    use run_parameters, only: fphi, fapar, fbpar, funits
    use theta_grid, only: ntgrid, theta, delthet, jacob, gradpar, nperiod
    use theta_grid, only: Rplot, Zplot, aplot, Rprime, Zprime, aprime
    use theta_grid, only: drhodpsi, qval, shape
    use kt_grids, only: naky, ntheta0, theta0, nx, ny, aky_out, akx_out, aky, akx
    use le_grids, only: nlambda, negrid, fcheck, al, delal
    use le_grids, only: e, dele
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use dist_fn, only: getan, get_epar, getmoms, par_spectrum, lambda_flux
    use dist_fn, only: e_flux
    use dist_fn, only: write_f, write_fyx
    use dist_fn, only: get_verr, get_gtran, write_poly, collision_error
    use dist_fn, only: g_adjust
    use collisions, only: vnmult
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: xxf_lo
    use gs2_transforms, only: transform2, inverse2
    use gs2_save, only: gs2_save_for_restart
    use constants
    use gs2_time, only: user_time, user_dt
    use gs2_io, only: nc_eigenfunc, nc_final_fields, nc_final_epar, nc_final_an
    use gs2_io, only: nc_final_moments, nc_finish
    use antenna, only: dump_ant_amp
    use splines, only: fitp_surf1, fitp_surf2
    use gs2_heating, only: del_htype, del_dvtype
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
    real, dimension (:), allocatable :: dl_over_b
    real, dimension (:,:,:), allocatable :: lamflux, enflux
    complex, dimension (ntheta0, naky) :: phi0
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
       if (write_ascii .and. write_density_velocity) call close_output_file (dv_unit)
       if (write_ascii .and. write_jext) call close_output_file (jext_unit)
       if (write_ascii .and. write_Epolar) then
          call close_output_file (polar_raw_unit)
          call close_output_file (polar_avg_unit)
       endif
       if (dump_neoclassical_flux) close (dump_neoclassical_flux_unit)

       if (dump_check1) close (dump_check1_unit)
       if (dump_check2) call close_output_file (dump_check2_unit)

       if (write_eigenfunc) then
          if (write_ascii) call open_output_file (unit, ".eigenfunc")
          phi0 = phi(0,:,:)

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
                           theta(ig), theta0(it,ik), aky_out(ik), &
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

       !Finish polar spectrum diagnostic (deallocate variables)
       if (write_Epolar) call finish_polar_spectrum

       if (write_final_fields) then
          if (write_ascii) then
             call open_output_file (unit, ".fields")
             do ik = 1, naky
                do it = 1, ntheta0
                   do ig = -ntg_out, ntg_out
                      write (unit, "(15(1x,e12.5))") &
!                           theta(ig), theta0(it,ik), aky_out(ik), &
                           theta(ig), aky_out(ik), akx_out(it), &
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
                        kpar(ig), aky_out(ik), akx_out(it), &
                        phi2(ig-ntgrid-1,it,ik), &
                        apar2(ig-ntgrid-1,it,ik), &
                        bpar2(ig-ntgrid-1,it,ik)                        
                end do
                do ig = 1, ntgrid
                   write (unit, "(9(1x,e12.5))") &
                        kpar(ig), aky_out(ik), akx_out(it), &
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
!                           theta(ig), theta0(it,ik), aky_out(ik), &
                           theta(ig), aky_out(ik), akx_out(it), &
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
    end if

    if (write_lamavg) then
       allocate (lamflux(nlambda, nspec, 4)) ; lamflux = 0.
       call lambda_flux (phinew, lamflux)
       if (proc0) then

          call open_output_file (unit, ".lam")
          do is = 1,nspec
             lamflux(:,is,1) = lamflux(:,is,1)*funits*spec(is)%dens
             do il=1,nlambda
                write (unit=unit, fmt=*) al(il), lamflux(il,is,1), is, &
                     sum(lamflux(1:il,is,1)*delal(1:il))
             end do
             write (unit=unit, fmt=*) 
          end do
          call close_output_file (unit)

          call open_output_file (unit, ".lame")
          do is = 1,nspec
             lamflux(:,is,2) = lamflux(:,is,2)*funits*spec(is)%dens*spec(is)%temp
             lamflux(:,is,3) = lamflux(:,is,3)*funits*spec(is)%dens*spec(is)%temp
             lamflux(:,is,4) = lamflux(:,is,4)*funits*spec(is)%dens*spec(is)%temp
             do il=1,nlambda
                write (unit=unit, fmt=*) al(il), lamflux(il,is,2), is, &
                     sum(lamflux(1:il,is,2)*delal(1:il)), &
                     lamflux(il,is,3), sum(lamflux(1:il,is,3)*delal(1:il)), &
                     lamflux(il,is,4), sum(lamflux(1:il,is,4)*delal(1:il))
             end do
             write (unit=unit, fmt=*) 
          end do
          call close_output_file (unit)
       end if
       deallocate (lamflux)
    end if
    
    if (write_eavg) then
       allocate (enflux(negrid, nspec, 4))
       call e_flux (phinew, enflux)
       if (proc0) then
          call open_output_file (unit, ".energy")
          do is = 1,nspec
             enflux(:,is,1) = enflux(:,is,1)*funits*spec(is)%dens
             do ie=1,negrid
                write (unit=unit, fmt=*) e(ie,is), enflux(ie,is,1),is, &
                     sum(enflux(1:ie,is,1)*dele(1:ie,is))
             end do
             write (unit=unit, fmt=*) 
          end do
          call close_output_file (unit)

          call open_output_file (unit, ".energye")
          do is = 1,nspec
             enflux(:,is,2) = enflux(:,is,2)*funits*spec(is)%dens*spec(is)%temp
             enflux(:,is,3) = enflux(:,is,3)*funits*spec(is)%dens*spec(is)%temp
             enflux(:,is,4) = enflux(:,is,4)*funits*spec(is)%dens*spec(is)%temp
             do ie=1,negrid
                write (unit=unit, &
                     fmt="(2(1x,e10.4),i3,5(1x,e10.4))") e(ie,is), enflux(ie,is,2), is, &
                     sum(enflux(1:ie,is,2)*dele(1:ie,is)), &
                     enflux(ie,is,3), sum(enflux(1:ie,is,3)*dele(1:ie,is)), &
                     enflux(ie,is,4), sum(enflux(1:ie,is,4)*dele(1:ie,is))
             end do
             write (unit=unit, fmt=*) 
          end do
          call close_output_file (unit)
       end if
       deallocate (enflux)
    end if

    if (write_tavg) then
       if (proc0) then
          call open_output_file (unit, ".theta")
          do is = 1,nspec
             do ig = -ntgrid, ntgrid
                write (unit=unit, fmt=*) theta(ig), theta_pflux(ig,is), is
             end do
             write (unit=unit, fmt=*) 
          end do
          call close_output_file (unit)

          call open_output_file (unit, ".thetae")
          do is=1,nspec
             do ig = -ntgrid,ntgrid
                write (unit=unit, &
                     fmt="(2(1x,e10.4),i3,5(1x,e10.4))") theta(ig), theta_qflux(ig,is,1), is, &
                     theta_qflux(ig,is,2), theta_qflux(ig,is,3)
             end do
             write (unit=unit, fmt=*) 
          end do
          call close_output_file (unit)
       end if
    end if
    
    call broadcast (write_final_moments)
    if (write_final_moments) then

       allocate (ntot(-ntgrid:ntgrid,ntheta0,naky,nspec))
       allocate (density(-ntgrid:ntgrid,ntheta0,naky,nspec))
       allocate (upar(-ntgrid:ntgrid,ntheta0,naky,nspec))
       allocate (tpar(-ntgrid:ntgrid,ntheta0,naky,nspec))
       allocate (tperp(-ntgrid:ntgrid,ntheta0,naky,nspec))
       call getmoms (ntot, density, upar, tpar, tperp)

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
                              theta(ig), aky_out(ik), akx_out(it), &
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
          call nc_final_moments (ntot, density, upar, tpar, tperp)

          if (write_ascii) then
             call open_output_file (unit, ".mom2")
             if (.not. write_eigenfunc) phi0 = 1.
             phi02=real(phi0*conjg(phi0))
             do is  = 1, nspec
                do ik = 1, naky
                   do it = 1, ntheta0
                      do ig = -ntg_out, ntg_out
                         write (unit, "(15(1x,e12.5))") &
                              theta(ig), aky_out(ik), akx_out(it), &
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
                   write (unit, "(i2,14(1x,e10.3))") spec(is)%type, akx_out(it), &
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
       deallocate (ntot, density, upar, tpar, tperp)
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
                           theta(ig), theta0(it,ik), aky_out(ik), &
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

    if (dump_final_xfields) then
       nny = 2*ny
       nnx = 2*nx
       allocate (xphi(nnx,nny,-ntgrid:ntgrid))
       call transform2 (phi, xphi, nny, nnx)

       if (proc0 .and. size(aky) > 1) then
          call get_unused_unit (unit)
          open (unit=unit, file="dump.xfields", status="unknown")
          do ig = -ntgrid, ntgrid
             do ik = 1, nny
                do it = 1, nnx
                   write (unit, "(15(1x,e12.5))") &
                        2.0*pi/akx_out(2)*real(it-1-nnx/2)/real(nnx), &
                        2.0*pi/aky(2)*real(ik-1)/real(nny), &
                        theta(ig), &
                        xphi(it,ik,ig)
                end do
                write (unit, "()")
             end do
          end do
          close (unit=unit)
       end if
       deallocate (xphi)
    end if

    if (write_fcheck) then
       allocate (fcheck_f(nlambda,ntheta0,naky,nspec))
       call fcheck (g, fcheck_f)
       if (proc0) then
          call open_output_file (unit, ".fcheck")
          do il = 2, nlambda
             do is = 1, nspec
                do ik = 1, naky
                   do it = 1, ntheta0
                      write (unit, "(20(1x,e12.6))") 0.5*(al(il) + al(il-1)), &
                           aky_out(ik), akx_out(it), fcheck_f(il,it,ik,is), &
                           sum((al(2:il)-al(1:il-1))*fcheck_f(2:il,it,ik,is))
                   end do
                end do
             end do
          end do
          call close_output_file (unit)
       end if
       deallocate (fcheck_f)
    end if

    if (save_for_restart) then
       call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, &
            fphi, fapar, fbpar, .true.)
    end if

    !<DD> Added for saving distribution function
    IF (save_distfn) THEN
    	!Convert h to distribution function
    	CALL g_adjust(gnew,phinew,bparnew,fphi,fbpar)
    	
    	!Save dfn, fields and velocity grids to file
       	CALL gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, &
          	fphi, fapar, fbpar, .true.,.true.)
    	
        !Convert distribution function back to h
        CALL g_adjust(gnew,phinew,bparnew,-fphi,-fbpar)
    END IF

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
                        kpar(ig), aky_out(ik), akx_out(it), &
                        real(vx2(ig-ntgrid-1,it,ik)), &
                        real(vy2(ig-ntgrid-1,it,ik))
                end do
                do ig = 1, ntgrid
                   write (unit, "(9(1x,e12.5))") &
                        kpar(ig), aky_out(ik), akx_out(it), &
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
    if (write_density_velocity) then
       call del_dvtype (dv_hist)
       call del_dvtype (dvk_hist)
    end if
    if (allocated(h_hist)) deallocate (h_hist, hk_hist, hk)
    if (allocated(dv_hist)) deallocate (dv_hist, dvk_hist)
    if (allocated(j_ext_hist)) deallocate (j_ext_hist)
    if (allocated(omegahist)) deallocate (omegahist)
    if (allocated(pflux)) deallocate (pflux, qheat, vflux, vflux_par, vflux_perp, pmflux, qmheat, vmflux, &
         pbflux, qbheat, vbflux, theta_pflux, theta_vflux, theta_vflux_par, theta_vflux_perp, &
         theta_qflux, theta_pmflux, vflux0, vflux1, &
         theta_vmflux, theta_qmflux, theta_pbflux, theta_vbflux, theta_qbflux)
    if (allocated(bxf)) deallocate (bxf, byf, xx4, xx, yy4, yy, dz, total)
    if (allocated(pflux_avg)) deallocate (pflux_avg, qflux_avg, heat_avg, vflux_avg)

    wtmp_old = 0. ; nout = 1 ; nout_movie = 1
    initialized = .false.

1000  format(20(1x,1pg18.11))
  end subroutine finish_gs2_diagnostics

  subroutine loop_diagnostics (istep, exit, debopt)
    use species, only: nspec, spec, has_electron_species
    use theta_grid, only: theta, ntgrid, delthet, jacob
    use theta_grid, only: gradpar, nperiod
    use kt_grids, only: naky, ntheta0, aky_out, theta0, akx_out, aky, akx
    use kt_grids, only: nkpolar !, akpolar, akpolar_out
    use run_parameters, only: woutunits, funits, tunits, fapar, fphi, fbpar, eqzip
    use fields, only: phinew, aparnew, bparnew
    use fields, only: kperp, fieldlineavgphi, phinorm
    use dist_fn, only: flux, vortcheck, fieldcheck, get_stress, write_f, write_fyx
    use dist_fn, only: neoclassical_flux, omega0, gamma0, getmoms, par_spectrum, gettotmoms
    use dist_fn, only: get_verr, get_gtran, write_poly, collision_error, neoflux
    use dist_fn, only: getmoms_notgc, g_adjust, include_lowflow, lf_flux
    use dist_fn, only: flux_vs_theta_vs_vpa
    use dist_fn_arrays, only: g, gnew, aj0, vpa
    use collisions, only: ncheck, vnmult, vary_vnew
    use mp, only: proc0, broadcast, iproc, send, receive
    use file_utils, only: get_unused_unit, flush_output_file
    use prof, only: prof_entering, prof_leaving
    use gs2_time, only: user_time
    use gs2_io, only: nc_qflux, nc_vflux, nc_pflux, nc_loop, nc_loop_moments
    use gs2_io, only: nc_loop_fullmom, nc_loop_sym
    use gs2_io, only: nc_loop_stress, nc_loop_vres
    use gs2_io, only: nc_loop_movie, nc_write_fields
    use gs2_layouts, only: yxf_lo, g_lo
! MAB>
    use gs2_layouts, only: init_parity_layouts, idx, idx_local, proc_id
    use gs2_layouts, only: p_lo, ie_idx, is_idx, il_idx, it_idx, ik_idx
! <MAB
    use gs2_transforms, only: init_transforms, transform2
    use le_grids, only: nlambda, ng2, integrate_moment, negrid, integrate_kysum
    use nonlinear_terms, only: nonlin
    use antenna, only: antenna_w
!    use gs2_flux, only: check_flux
    use gs2_heating, only: heating_diagnostics, del_htype, &
         dens_vel_diagnostics,init_dvtype, del_dvtype
    use constants

    implicit none
!    integer :: nout = 1
!    integer :: nout_movie = 1
    integer, intent (in) :: istep
    logical, intent (out) :: exit
    real, dimension(:,:,:), allocatable :: yxphi, yxapar, yxbpar
    complex, dimension (ntheta0, naky) :: omega, omegaavg

    type (dens_vel_diagnostics) :: dv
    type (dens_vel_diagnostics), dimension(:,:), allocatable :: dvk
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
!<MAB
    real :: geavg, glavg, gtavg
    real :: dmix, dmix4, dmixx
    real :: t, denom
    integer :: ig, ik, it, is, unit, i, j, nnx, nny, ifield, write_mod
    complex :: sourcefac
!    complex :: phiavg
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky,nspec) :: ntot, density, &
         upar, tpar, tperp, upartot, uperptot, ttot
    complex, dimension (ntheta0, nspec) :: ntot00, density00, upar00, tpar00, tperp00
    complex, dimension (ntheta0, nspec) :: rstress, ustress
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

    phase_tot = 0.0 ;  phase_theta = 0.0

    call prof_entering ("loop_diagnostics")

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

    if (write_hrate) call heating (istep, h, hk)

!>GGH
    !Write density and velocity perturbations
    if (write_density_velocity) then
       call init_dvtype (dv,  nspec)
       allocate (dvk(ntheta0, naky))
       call init_dvtype (dvk, nspec)

       call dens_vel(istep,dv,dvk)
    endif
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
          call getmoms (ntot, density, upar, tpar, tperp)
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
       call flux (phinew, aparnew, bparnew, &
            pflux,  qheat,  vflux, vflux_par, vflux_perp, &
            pmflux, qmheat, vmflux,  &
            pbflux, qbheat, vbflux, &
            theta_pflux, theta_vflux, theta_vflux_par, theta_vflux_perp, theta_qflux, &
            theta_pmflux, theta_vmflux, theta_qmflux, & 
            theta_pbflux, theta_vbflux, theta_qbflux)
       ! lowflow terms only implemented in electrostatic limit at present
       if (include_lowflow) call lf_flux (phinew, vflux0, vflux1)
       call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar)
!       call flux (phinew, aparnew, bparnew, &
!            pflux, qheat, qheat_par, qheat_perp, vflux, &
!            pmflux, qmheat, qmheat_par, qmheat_perp, vmflux, &
!            pbflux, qbheat, qbheat_par, qbheat_perp, vbflux)

       if (proc0) then
          if (fphi > epsilon(0.0)) then
             do is = 1, nspec
                qheat(:,:,is,1) = qheat(:,:,is,1)*funits &
                     *spec(is)%dens*spec(is)%temp
                call get_volume_average (qheat(:,:,is,1), heat_fluxes(is))

                qheat(:,:,is,2) = qheat(:,:,is,2)*funits &
                     *spec(is)%dens*spec(is)%temp
                call get_volume_average (qheat(:,:,is,2), heat_par(is))

                qheat(:,:,is,3) = qheat(:,:,is,3)*funits &
                     *spec(is)%dens*spec(is)%temp
                call get_volume_average (qheat(:,:,is,3), heat_perp(is))
                
                theta_qflux(:,is,1) = theta_qflux(:,is,1)*funits &
                     *spec(is)%dens*spec(is)%temp
                theta_qflux(:,is,2) = theta_qflux(:,is,2)*funits &
                     *spec(is)%dens*spec(is)%temp
                theta_qflux(:,is,3) = theta_qflux(:,is,3)*funits &
                     *spec(is)%dens*spec(is)%temp

                pflux(:,:,is) = pflux(:,:,is)*funits &
                     *spec(is)%dens
                call get_volume_average (pflux(:,:,is), part_fluxes(is))

                theta_pflux(:,is) = theta_pflux(:,is)*funits &
                     *spec(is)%dens
                
                vflux(:,:,is) = vflux(:,:,is)*funits**2 &
                     *spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
                call get_volume_average (vflux(:,:,is), mom_fluxes(is))
                vflux_par(:,:,is) = vflux_par(:,:,is)*funits**2 &
                     *spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
                call get_volume_average (vflux_par(:,:,is), parmom_fluxes(is))
                vflux_perp(:,:,is) = vflux_perp(:,:,is)*funits**2 &
                     *spec(is)%dens*spec(is)%mass*spec(is)%stm
                call get_volume_average (vflux_perp(:,:,is), perpmom_fluxes(is))

                if (include_lowflow) then
                   vflux0(:,:,is) = vflux0(:,:,is)*funits**2 &
                        *spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
                   call get_volume_average (vflux0(:,:,is), lfmom_fluxes(is))
                   vflux1(:,:,is) = vflux1(:,:,is)*funits**2 &
                        *spec(is)%dens*spec(is)%mass*spec(is)%temp/spec(is)%z
                   call get_volume_average (vflux1(:,:,is), vflux1_avg(is))
! TMP UNTIL VFLUX0 IS TESTED
!                   mom_fluxes = mom_fluxes + lfmom_fluxes
                end if

                theta_vflux(:,is) = theta_vflux(:,is)*funits**2 &
                     *spec(is)%dens*spec(is)%mass*spec(is)%stm
                theta_vflux_par(:,is) = theta_vflux_par(:,is)*funits**2 &
                     *spec(is)%dens*spec(is)%mass*spec(is)%stm
                theta_vflux_perp(:,is) = theta_vflux_perp(:,is)*funits**2 &
                     *spec(is)%dens*spec(is)%mass*spec(is)%stm

             end do
          end if
          if (fapar > epsilon(0.0)) then
             do is = 1, nspec
                qmheat(:,:,is,1)=qmheat(:,:,is,1)*funits &
                     *spec(is)%dens*spec(is)%temp
                call get_volume_average (qmheat(:,:,is,1), mheat_fluxes(is))

                call get_surf_average (qmheat(:,:,is,1), x_qmflux(:,is))

                qmheat(:,:,is,2)=qmheat(:,:,is,2)*funits &
                     *spec(is)%dens*spec(is)%temp
                call get_volume_average (qmheat(:,:,is,2), mheat_par(is))

                qmheat(:,:,is,3)=qmheat(:,:,is,3)*funits &
                     *spec(is)%dens*spec(is)%temp
                call get_volume_average (qmheat(:,:,is,3), mheat_perp(is))
                
                pmflux(:,:,is)=pmflux(:,:,is)*funits &
                     *spec(is)%dens
                call get_volume_average (pmflux(:,:,is), mpart_fluxes(is))

                vmflux(:,:,is)=vmflux(:,:,is)*funits &
                     *spec(is)%dens*spec(is)%mass*spec(is)%stm
                call get_volume_average (vmflux(:,:,is), mmom_fluxes(is))
             end do
          end if
          if (fbpar > epsilon(0.0)) then
             do is = 1, nspec
                qbheat(:,:,is,1)=qbheat(:,:,is,1)*funits &
                     *spec(is)%dens*spec(is)%temp
                call get_volume_average (qbheat(:,:,is,1), bheat_fluxes(is))

                qbheat(:,:,is,2)=qbheat(:,:,is,2)*funits &
                     *spec(is)%dens*spec(is)%temp
                call get_volume_average (qbheat(:,:,is,2), bheat_par(is))

                qbheat(:,:,is,3)=qbheat(:,:,is,3)*funits &
                     *spec(is)%dens*spec(is)%temp
                call get_volume_average (qbheat(:,:,is,3), bheat_perp(is))
                
                pbflux(:,:,is)=pbflux(:,:,is)*funits &
                     *spec(is)%dens
                call get_volume_average (pbflux(:,:,is), bpart_fluxes(is))

                vbflux(:,:,is)=vbflux(:,:,is)*funits &
                     *spec(is)%dens*spec(is)%mass*spec(is)%stm
                call get_volume_average (vbflux(:,:,is), bmom_fluxes(is))
             end do
          end if
          pflux_avg = pflux_avg + (part_fluxes + mpart_fluxes + bpart_fluxes)*(t-t_old)
          qflux_avg = qflux_avg + (heat_fluxes + mheat_fluxes + bheat_fluxes)*(t-t_old)
          vflux_avg = vflux_avg + (mom_fluxes + mmom_fluxes + bmom_fluxes)*(t-t_old)
          if (write_hrate) heat_avg = heat_avg + h%imp_colls*(t-t_old)
          t_old = t
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
          if (print_old_units) then
             do ik = 1, naky
                do it = 1, ntheta0
                   write (unit=*, fmt="('aky=',f5.2, ' th0=',f5.2, &
                          &' om=',2f8.3,' omav=', 2f8.3,' phtot=',e10.4)") &
                        aky_out(ik), theta0(it,ik), &
                        real (omega(it,ik)), &
                        aimag(omega(it,ik)), &
                        real( omegaavg(it,ik)), &
                        aimag(omegaavg(it,ik)), &
                        phitot(it,ik)
                end do
             end do
          else
             do ik = 1, naky
                do it = 1, ntheta0
!                   write (unit=*, fmt="('aky=',f5.2, ' th0=',f7.2, &
!                          &' om=',2f8.3,' omav=', 2f8.3,' phtot=',e10.4)") &
!                        aky_out(ik), theta0(it,ik), &
!                        real( omega(it,ik)*woutunits(ik)), &
!                        aimag(omega(it,ik)*woutunits(ik)), &
!                        real( omegaavg(it,ik)*woutunits(ik)), &
!                        aimag(omegaavg(it,ik)*woutunits(ik)), &
!                        phitot(it,ik)
                   write (unit=*, fmt="('ky=',f7.4, ' kx=',f7.4, &
                          &' om=',2f8.3,' omav=', 2f8.3,' phtot=',e8.2)") &
                        aky_out(ik), akx_out(it), &
                        real( omega(it,ik)*woutunits(ik)), &
                        aimag(omega(it,ik)*woutunits(ik)), &
                        real( omegaavg(it,ik)*woutunits(ik)), &
                        aimag(omegaavg(it,ik)*woutunits(ik)), &
                        phitot(it,ik)
                end do
             end do
          end if
          write (*,*) 
       end if
    end if

    i=istep/nwrite
!    call check_flux (i, t, heat_fluxes)
! Polar spectrum calculation----------------------------------------------
    if (write_Epolar .and. proc0) then
       ebinarray(:,:)=0.
       
       !Calculate polar spectrum of energies 
       call get_polar_spectrum (hk%energy, ebinarray(:,iefluc))
       call get_polar_spectrum (hk%eapar,  ebinarray(:,ieapar))
       call get_polar_spectrum (hk%ebpar,  ebinarray(:,iebpar))

       do is=1,nspec
          do ik= 1,naky
             do it=1,ntheta0
                etmp(it,ik)=hk(it,ik)%phis2(is)
             enddo
          enddo
          call get_polar_spectrum(etmp(:,:), ebinarray(:,iephis2  + (is-1)*3))
          do ik= 1,naky
             do it=1,ntheta0
                etmp(it,ik)=hk(it,ik)%hs2(is)
             enddo
          enddo
          call get_polar_spectrum(etmp(:,:), ebinarray(:,iehs2    + (is-1)*3))
          do ik= 1,naky
             do it=1,ntheta0
                etmp(it,ik)=hk(it,ik)%delfs2(is)
             enddo
          enddo
          call get_polar_spectrum(etmp(:,:), ebinarray(:,iedelfs2 + (is-1)*3))
       enddo
       
!
! BD:
! --- must have write_hrate = T for Epolar to work b/c hk array is used

       !Output raw kspectrum to file
       if (nspec == 1) then
          do i=1,nbx
             write (unit=polar_raw_unit, fmt="('t= ',e16.10,' kperp= ',e10.4, &
                  &' E= '     ,e10.4,' Eapar= ',e10.4,' Ebpar= ',e10.4, &
                  &' Ephis2= ',e10.4,' Ehss2= ',e10.4,' Edfs2= ',e10.4)") &
                  & t, kpbin(i),ebinarray(i,1:6)
          end do
       else  !Only writing this data for first two species for now
! Labels assume first species is ion species.
          do i=1,nbx
             write (unit=polar_raw_unit, fmt="('t= ',e16.10,' kperp= ',e10.4, &
                  &' E= '     ,e10.4,' Eapar= ',e10.4,' Ebpar= ',e10.4, &
                  &' Ephii2= ',e10.4,' Ehsi2= ',e10.4,' Edfi2= ',e10.4, &
                  &' Ephie2= ',e10.4,' Ehse2= ',e10.4,' Edfe2= ',e10.4)") &
                  & t, kpbin(i),ebinarray(i,1:9)
          end do
       end if
       write (unit=polar_raw_unit, fmt='(a)') ''      

       !Compute log-averaged polar spectrum
       eavgarray(:,:)=0.
       do ifield = 1, 3+nspec*3

!ERROR- Below is an attempt to be more efficient, but it is not right
!          if (ifield == 2 .and. fapar > epsilon(0.0)) then
!             continue
!          else
!             cycle
!          end if

!          if (ifield == 3 .and. fbpar > epsilon(0.0)) then
!             continue
!          else
!             cycle
!          end if

         do i = 1,nbx
             j=polar_avg_index(i)
             eavgarray(j,ifield)=eavgarray(j,ifield)+log(ebinarray(i,ifield))
          enddo
! Finish log-averaging
          do j=1,nkpolar
             eavgarray(j,ifield)=eavgarray(j,ifield)/numavg(j)
          enddo
          eavgarray(:,ifield)=exp(eavgarray(:,ifield))
       end do

! Output log-averaged kspectrum to file
       if (nspec == 1) then
          do i=1,nkpolar
             write (unit=polar_avg_unit, fmt="('t= ',e16.10,' kperp= ',e10.4, &
                  &' E= '     ,e10.4,' Eapar= ',e10.4,' Ebpar= ',e10.4, &
                  &' Ephis2= ',e10.4,' Ehss2= ',e10.4,' Edfs2= ',e10.4)") &
                  & t, kpavg(i),eavgarray(i,1:6)
          end do
       else ! Only writing this data for first two species right now
! Labels assume first species is ion species.
          do i=1,nkpolar
             write (unit=polar_avg_unit, fmt="('t= ',e16.10,' kperp= ',e10.4, &
                  &' E= '     ,e10.4,' Eapar= ',e10.4,' Ebpar= ',e10.4, &
                  &' Ephii2= ',e10.4,' Ehsi2= ',e10.4,' Edfi2= ',e10.4, &
                  &' Ephie2= ',e10.4,' Ehse2= ',e10.4,' Edfe2= ',e10.4)") &
                  & t, kpavg(i),eavgarray(i,1:9)
          end do
       end if
       write (unit=polar_avg_unit, fmt='(a)') ''      

    end if !END Polar spectrum calculation----------------------------------

    if (write_vortcheck) call vortcheck (phinew, bparnew)
    if (write_fieldcheck) call fieldcheck (phinew, aparnew, bparnew)
    if (write_fields) call nc_write_fields (nout, phinew, aparnew, bparnew)  !MR

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

    if (write_neoclassical_flux .or. dump_neoclassical_flux .and. neoflux) then
!       call neoclassical_flux (pfluxneo, qfluxneo, phinew, bparnew)
       call neoclassical_flux (pfluxneo, qfluxneo)
    end if

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

!<GGH
       !Write out data for density and velocity peturbations
       if (write_ascii .and. write_density_velocity) then
          write (unit=dv_unit, fmt="('t= ',e12.6,' dvpar= ',e12.6,' ',e12.6,' dvperp= ',e12.6,' ',e12.6,' dn= ',e12.6,' ',e12.6)")  &
               t, dv % dvpar(:), dv % dvperp(:), dv % dn(:)
!          write (unit=dv_unit, fmt="('t= ',e12.6,' dvperp= ',e12.6)") t, dv % dvperp
!          write (unit=dv_unit, fmt="('t= ',e12.6,' dn= ',e12.6)") t, dv % dn
!          write (unit=heat_unit, fmt='(a)') ''
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
                     t, phi2, heat_fluxes(1:min(nspec,5)), qflux_avg(1:min(nspec,5))/t
                write (unit=out_unit, fmt="('t= ',e16.10,' <phi**2>= ',e10.4, &
                     & ' part fluxes: ', 5(1x,e10.4),' pflux_avg: ', 5(1x,e10.4))") &
                     t, phi2, part_fluxes(1:min(nspec,5)), pflux_avg(1:min(nspec,5))/t
                write (unit=out_unit, fmt="('t= ',e16.10,' <phi**2>= ',e10.4, &
                     & ' mom fluxes: ', 5(1x,e10.4),' vflux_avg: ', 5(1x,e10.4))") &
                     t, phi2, mom_fluxes(1:min(nspec,5)), vflux_avg(1:min(nspec,5))/t
                write (unit=out_unit, fmt="('t= ',e16.10,' <phi**2>= ',e10.4, &
                     & ' parmom fluxes: ', 5(1x,e10.4))") &
                     t, phi2, parmom_fluxes(1:min(nspec,5))
                write (unit=out_unit, fmt="('t= ',e16.10,' <phi**2>= ',e10.4, &
                     & ' perpmom fluxes: ', 5(1x,e10.4))") &
                     t, phi2, perpmom_fluxes(1:min(nspec,5))
                if (include_lowflow) then
                   write (unit=out_unit, fmt="('t= ',e16.10,' <phi**2>= ',e10.4, &
                        & ' lfmom fluxes: ', 5(1x,e10.4),' lfvflx1: ', 5(1x,e10.4))") &
                        t, phi2, lfmom_fluxes(1:min(nspec,5)), vflux1_avg(1:min(nspec,5))
                end if
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
                  heat_fluxes, mheat_fluxes, bheat_fluxes, x_qmflux, hflux_tot)
!          call nc_qflux (nout, qheat, qmheat, qbheat, &
!               heat_par, mheat_par, bheat_par, &
!               heat_perp, mheat_perp, bheat_perp, &
!               heat_fluxes, mheat_fluxes, bheat_fluxes, hflux_tot)
             call nc_vflux (nout, vflux, vmflux, vbflux, &
                  mom_fluxes, mmom_fluxes, bmom_fluxes, vflux_tot, &
                  vflux_par, vflux_perp, vflux0, vflux1)
             call nc_pflux (nout, pflux, pmflux, pbflux, &
                  part_fluxes, mpart_fluxes, bpart_fluxes, zflux_tot)
          end if
          call nc_loop (nout, t, fluxfac, &
               phinew(igomega,:,:), phi2, phi2_by_mode, &
               aparnew(igomega,:,:), apar2, apar2_by_mode, &
               bparnew(igomega,:,:), bpar2, bpar2_by_mode, &
               h, hk, omega, omegaavg, woutunits, phitot, write_omega, write_hrate)
       end if
       if (write_ascii) then
          do ik = 1, naky
             do it = 1, ntheta0
!                write (out_unit,*) "aky=",aky_out(ik)," theta0=",theta0(it,ik)

                if (write_line) then
                   write (out_unit, "('t= ',e16.10,' aky= ',1pe12.4, ' akx= ',1pe12.4, &
                        &' om= ',1p2e12.4,' omav= ', 1p2e12.4,' phtot= ',1pe12.4)") &
                        t, aky_out(ik), akx_out(it), &
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
                
                if (write_dmix) then
                   if (abs(akperp(it,ik)) < 2.0*epsilon(0.0)) then
                      write (out_unit, *) 'dmix indeterminate'
                   else
                      dmix = 2.0*woutunits(ik)*aimag(omega(it,ik))/akperp(it,ik)**2
                      dmix4 = dmix*dmix/(woutunits(ik)*aimag(omega(it,ik)))
                      dmixx = 2.0*woutunits(ik)*aimag(omega(it,ik))/(akperp(it,ik)**2-aky(ik)**2)
!                      write (out_unit, *) 'dmix=', dmix,' dmix4= ',dmix4,' dmixx= ',dmixx
                      write (out_unit,fmt='(" dmix=",1pe12.4," dmix4=",1pe12.4," dmixx=",1pe12.4)') dmix, dmix4, dmixx
                   end if
                end if
             end do
          end do
       end if
       do ik = 1, naky
          do it = 1, ntheta0
             if (write_kperpnorm) then
                if (aky_out(ik) /= 0.0) then
                   write (out_unit, *) 'kperpnorm=',akperp(it,ik)/aky_out(ik)
                else
! huh? 
                   write (out_unit, *) 'kperpnorm*aky=',akperp(it,ik)
                end if
             end if
!             if (write_phitot) write (out_unit, *) 'phitot=', phitot(it,ik)

!             if (write_fieldline_avg_phi) then
!                call fieldlineavgphi (ntg_out, it, ik, phiavg)
!                write (out_unit, *) '<phi>=',phiavg
!             end if

             if (write_nl_flux) then
                if (write_ascii) then
!                   write (out_unit,"('ik,it,aky,akx,<phi**2>,t: ', &
!                        & 2i5,4(1x,e12.6))") &
!                        ik, it, aky_out(ik), akx_out(it), phi2_by_mode(it,ik), t
!                   write (out_unit,"('ik,it,aky,akx,<apar**2>,t: ', &
!                        & 2i5,4(1x,e12.6))") &
!                        ik, it, aky_out(ik), akx_out(it), apar2_by_mode(it,ik), t
                end if
!                if (ntheta0 > 1 .and. ik == 1) then
!                   write (out_unit, "('it,akx,<phi**2>,t: ',i5,3(1x,e12.6))") &
!                        it, akx_out(it), sum(phi2_by_mode(it,:)*fluxfac(:)), t
!                end if
!                if (naky > 1 .and. it == 1) then
!                   write (out_unit, "('ik,aky,<phi**2>,t: ',i5,3(1x,e12.6))") &
!                        ik, aky_out(ik), sum(phi2_by_mode(:,ik))*fluxfac(ik), t
!                end if
             end if
          end do
       end do
    end if

    if (write_cerr) call collision_error(phinew,bparnew,last)

    if (write_stress) then
       call get_stress (rstress, ustress)
       if (proc0) call nc_loop_stress(nout, rstress, ustress)
    end if

     call broadcast (write_symmetry)
     if (write_symmetry) then
        allocate (vflx_sym(-ntgrid:ntgrid,nlambda*negrid,nspec))
        call flux_vs_theta_vs_vpa (vflx_sym)
        if (proc0) call nc_loop_sym (nout, vflx_sym)
        deallocate (vflx_sym)
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

        ! TMP FOR TESTING -- MAB
!        do ik = 1, naky
!           do it = 1, ntheta0
!              phinew(:,it,ik) = sin(theta)**2
!           end do
!        end do
!         do iglo = g_lo%llim_proc, g_lo%ulim_proc
!            ik = ik_idx(g_lo,iglo)
!            it = it_idx(g_lo,iglo)
!            do isgn = 1, 2
! !              gnew(:,isgn,iglo) = vpa(:,isgn,iglo)**2*(sin(theta)*akx(it)+cos(theta)*akx(it)**2)
! !              gnew(:,isgn,iglo) = sin(theta)*vpa(:,isgn,iglo)
! ! !             gnew(:,isgn,iglo) = vpa(:,isgn,iglo)*cos(theta)
! ! !             gnew(:,isgn,iglo) = vpa(:,isgn,iglo)*cos(theta)*akx(it)
! ! !             gnew(:,isgn,iglo) = akx(it)
!            end do
!         end do

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

    call broadcast (test_conserve)
    call broadcast (write_avg_moments)
    if (write_avg_moments) then

       call getmoms (ntot, density, upar, tpar, tperp)

       if (test_conserve) call gettotmoms (ntot, upartot, uperptot, ttot)

       if (proc0) then
          allocate (dl_over_b(-ntg_out:ntg_out))

          dl_over_b = delthet(-ntg_out:ntg_out)*jacob(-ntg_out:ntg_out)
          dl_over_b = dl_over_b / sum(dl_over_b)

          if (test_conserve) then
             do is = 1, nspec
                do it = 1, ntheta0
                   ntot00(it, is)   = sum( ntot(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                   density00(it,is) = sum(density(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                   upar00(it, is)   = sum( upartot(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                   tpar00(it, is)   = sum( ttot(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                   tperp00(it, is)  = sum(uperptot(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                end do
             end do
             if (nspec == 1) then
                write (out_unit, *) 'moms', t, real(ntot00(1,1)), aimag(upar00(1,1)), real(tperp00(1,1)), &
                     real(tpar00(1,1))
             else
                write (out_unit, *) 'moms', t, real(ntot00(1,1)), real(ntot00(1,2)), aimag(upar00(1,1)), &
                     aimag(upar00(1,2)), aimag(tperp00(1,1)), aimag(tperp00(1,2)), real(tpar00(1,1)), &
                     real(tpar00(1,2))
             end if
          else
             do is = 1, nspec
                do it = 1, ntheta0
                   ntot00(it, is)   = sum( ntot(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                   density00(it,is) = sum(density(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                   upar00(it, is)   = sum( upar(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                   tpar00(it, is)   = sum( tpar(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                   tperp00(it, is)  = sum(tperp(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                end do
             end do
          end if

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
    ! <RN

    if (write_neoclassical_flux .and. neoflux .and. proc0) then
       do is=1,nspec
          tprim_tot(is) = spec(is)%tprim-0.333*aimag(tpar00(1,is))-0.6667*aimag(tperp00(1,is))
          fprim_tot(is) = spec(is)%fprim-aimag(density00(1,is))
          write (out_unit, "('t= ',e12.6,' is= ',i2,' pfluxneo= ',e12.6,' nprim_tot= ',e12.6,' tprim_tot= ',e12.6)") &
               t, is, real(pfluxneo(1,1,is)), fprim_tot(is), tprim_tot(is)
          write (out_unit, "('t= ',e12.6,' is= ',i2,' qfluxneo= ',e12.6,' tprim_tot= ',3(1x,e12.6))") &
               t, is, real(qfluxneo(1,1,is)), tprim_tot(is), -0.333*aimag(tpar00(1,is)), -0.6667*aimag(tperp00(1,is))
       end do
    end if


    if (proc0 .and. dump_any) then
!
! I have not checked the units in this section. BD
!       
       do ik = 1, naky
          do it = 1, ntheta0
!             it = 2
             if (dump_neoclassical_flux .and. neoflux) then
                do is = 1, nspec
                   write (dump_neoclassical_flux_unit, "(20(1x,e12.5))") &
                        t, aky_out(ik), akx_out(it), real(is), pfluxneo(it,ik,is), &
                        qfluxneo(it,ik,is), &
                        pfluxneo(it,ik,is)/sourcefac, &
                        qfluxneo(it,ik,is)/sourcefac
                end do
             end if
             if (dump_check1) then
                denom=sum(delthet(-ntg_out:ntg_out-1)*jacob(-ntg_out:ntg_out-1))
                write (dump_check1_unit, "(20(1x,e12.5))") &
                     t, aky_out(ik), akx_out(it), &
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
                     t, aky_out(ik), akx_out(it), aparnew(igomega,it,ik)
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
                     theta(ig), aky_out(ik), theta0(it,ik), &
                     phinew(ig,it,ik), aparnew(ig,it,ik), &
                     bparnew(ig,it,ik), t, akx_out(it)
             end do
             write (unit, "()")
          end do
       end do
       close (unit=unit)
    end if
    
    nout = nout + 1
    if (write_ascii .and. mod(nout, 10) == 0 .and. proc0) then
       call flush_output_file (out_unit, ".out")
       if (write_verr) then
          call flush_output_file (res_unit, ".vres")
          call flush_output_file (lpc_unit, ".lpc")
       end if
       if (write_parity) call flush_output_file (parity_unit, ".parity")
    end if

!>GGH
    !Deallocate variables for density and velocity perturbations
    if (write_density_velocity) then
       call del_dvtype (dv)
       call del_dvtype (dvk)
       deallocate(dvk)
    endif
! Deallocate variable for Jexternal
    if (write_jext) deallocate(j_ext)
!<GGH

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
    use run_parameters, only: funits
    use nonlinear_terms, only: nonlin
    use dist_fn_arrays, only: c_rate
    use gs2_heating, only: heating_diagnostics, avg_h, avg_hk, htimesx, zero_htype
    implicit none
    integer, intent (in) :: istep
    type (heating_diagnostics) :: h
    type (heating_diagnostics), dimension(:,:) :: hk

    real, dimension(-ntgrid:ntgrid) :: wgt
    real :: fac
    integer :: is, ik, it, ig
    
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

    call htimesx (h, funits)
    call htimesx (hk, funits)

    call avg_h(h, h_hist, istep, navg)
    call avg_hk(hk, hk_hist, istep, navg)

  end subroutine heating
!>GGH
!=============================================================================
! Density: Calculate Density perturbations
!=============================================================================
 subroutine dens_vel (istep, dv, dvk)
    use dist_fn, only: get_dens_vel
    use fields, only: phi, apar, bpar, phinew, aparnew, bparnew
    use gs2_heating, only: dens_vel_diagnostics, avg_dv, avg_dvk
    implicit none
    !Passed
    integer, intent (in) :: istep
    type (dens_vel_diagnostics) :: dv
    type (dens_vel_diagnostics), dimension(:,:) :: dvk

    !Call routine to calculate density and velocity perturbations
!    call get_dens_vel(dv, dvk, phi, apar, bpar, phinew, aparnew, bparnew)    
    call get_dens_vel(dv, dvk, phi, bpar, phinew, bparnew)
    
    !Do averages with a history variable
    call avg_dv(dv, dv_hist, istep, navg)
    call avg_dvk(dvk, dvk_hist, istep, navg)

  end subroutine dens_vel
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
!================================================================================
! Set up corrections for polar energy spectrum
!================================================================================
!NOTE: Here we calculate the correction factors for each possible kperp
  subroutine init_polar_spectrum
    use dist_fn_arrays, only: kperp2
    use kt_grids, only: naky, ntheta0, aky, akx, nkpolar, ikx, iky
    use constants, only: pi
    use nonlinear_terms, only: nonlin
    use species, only: nspec
    implicit none
!    real, dimension(:,:), allocatable :: kp_by_mode
    integer, dimension(:), allocatable :: num, nbin
    real, dimension(:), allocatable :: kp
    real, dimension(:), allocatable :: kpavg_lim
    real :: kpmax,dkp
    integer :: nkperp                           !Total possible number of kperps
    integer :: ik, it, i,j,inbx !,nkx,nky
    integer :: ig = 0

! Not yet properly generalized for kperp2 = kperp2 (theta).  Currently 
! set up using kperp**2 values calculated at theta = 0 (by setting ig = 0 above).

! NOTE: In this routine, a square domain is assumed! (nx=ny)
! In read_parameters, write_Epolar => .false. if nx /= ny

    !Determine total number of possible kperps and allocate array
    nkperp = ntheta0**2 + naky**2
    allocate (num(1:nkperp)); num=0
    allocate (kp(1:nkperp)); kp(:)=huge(kp)
!    allocate(kp_by_mode(ntheta0,naky)) ; kp_by_mode(:,:)=0.
    allocate (polar_index(ntheta0,naky)) ; polar_index(:,:)=0

! Loop through all modes and sum number at each kperp
    do it = 1,ntheta0
       do ik= 1,naky
          if (nonlin .and. it == 1 .and. ik == 1) cycle

! Add to number and calculate magnitude of this kperp
          i = ikx(it)**2 + iky(ik)**2
          num(i)=num(i)+1
! In unsheared slab, this is kp = kperp.  Generalize to kp = sqrt(kperp2)
!         kp(i)=sqrt(akx(it)**2+aky(ik)**2)
          kp(i) = sqrt(kperp2(ig, it, ik))
          polar_index(it,ik)=i
       enddo
    enddo
   
    !Collapse bins to only existing values of kperp
    !Find total number of existing kperps
    nbx=0
    do i=1,nkperp
       if (num(i) > 0) nbx=nbx+1
    enddo
    
    !Allocate bin variables
    allocate (nbin(1:nbx)); nbin=0
    allocate (kpbin(1:nbx)); kpbin=0.
    allocate (ebincorr(1:nbx)); ebincorr=0.
    allocate (ebinarray(1:nbx,3+3*nspec)); ebinarray=0. !Used in loop diagnostics
    allocate (etmp(1:ntheta0,1:naky)); etmp=0.

    !Copy data
    inbx=0
    do i=1,nkperp
       if (num(i) > 0) then
          inbx=inbx+1
          nbin(inbx)=num(i)
          kpbin(inbx)=kp(i)
          !Correct polar_index
          do it = 1,ntheta0
             do ik= 1,naky
                if (polar_index(it,ik) == i) polar_index(it,ik)=inbx
             enddo
          enddo
       endif
    enddo

    !Calculate the correction factor for the discrete values
    ebincorr=pi*kpbin/real(nbin)

    !Deallocate variables
    deallocate(num,kp,nbin)

    !Allocate variables for log-averaged polar spectrum
    allocate(kpavg_lim(1:nkpolar+1)); kpavg_lim(:)=0.
    allocate(numavg(1:nkpolar)) ; numavg(:)=0.
    allocate(kpavg(1:nkpolar)) ; kpavg(:)=0.
    allocate(eavgarray(1:nkpolar,9)); eavgarray(:,:)=0. !Used in loop diagnostics
    allocate(polar_avg_index(1:nbx)) ; polar_avg_index(:)=0

    !Determine limits of log-averaged kperp spectra bins
    kpmax=real(int(kpbin(nbx)/kpbin(1))+1)*kpbin(1)
    dkp=(kpmax-kpbin(1))/real(nkpolar)
    do i=1,nkpolar+1
       kpavg_lim(i)=dkp*real(i)
    enddo

    !Do log-average of kperp in each bin and build index
    do i =1,nbx
       do j=1,nkpolar
          if (kpbin(i) .ge. kpavg_lim(j) .and. kpbin(i) .lt. kpavg_lim(j+1) ) then
             numavg(j)=numavg(j)+1.
             kpavg(j)=kpavg(j)+log(kpbin(i))
             polar_avg_index(i)=j
          endif
       enddo
    enddo
    !Finish log-averaging
    kpavg(:)=kpavg(:)/numavg(:)
    kpavg(:)=exp(kpavg(:))

    deallocate(kpavg_lim)

    !Debug
!    if (0. .eq. 0.) then
!       do i=1,nbx
!          write(*,'(i8,es12.4,i8,es12.4)')i, kpbin(i), polar_avg_index(i)
!       enddo
!       write(*,*)'Total kperps= ',nbx,' Total binned= ',sum(numavg)
!       do i=1,nkpolar
!          write(*,'(i8,es12.4,i8)')i, kpavg(i), int(numavg(i))
!       enddo
!    endif

    !Debug
!    if (0. .eq. 1.) then
!       do i=1,nbx
!          write(*,'(i8,es12.4,i8,es12.4)')i, kpbin(i), int(nbin(i)), ebincorr(i)
!       enddo
!       do it = 1,ntheta0
!          do ik= 1,naky
!             nkx=-(mod((it-1)+int(ntheta0/2),ntheta0)-int(ntheta0/2))
!             nky=ik-1
!            write(*,'(5i8,2es12.4)') it,ik,nkx,nky,polar_index(it,ik),akx(it),aky(ik)
!          enddo
!       enddo
!    endif

  end subroutine init_polar_spectrum
!================================================================================
! Calculate corrected polar spectrum 
!================================================================================
! NOTE: polar_index connects a given mode to the correct binned kperp value
! NOTE: It is assumed here that weighting factors for ky=0 are already 
! incorporated in ee (already done in heating calculations)
  subroutine get_polar_spectrum(ee,ebin)
    use kt_grids, only: naky, ntheta0
    use nonlinear_terms, only: nonlin
    implicit none
    real, dimension (:,:), intent (in) :: ee   !variable ee by (kx,ky) mode
    real, dimension (:), intent (out) :: ebin  !variable ee in kperp bins
    integer :: ik, it

    !Initialize ebin
    ebin=0.

    !Loop through all modes and sum number at each kperp
    do it = 1,ntheta0
       do ik= 1,naky
          if (nonlin .and. it == 1 .and. ik == 1) cycle
          ebin(polar_index(it,ik))= ebin(polar_index(it,ik))+ee(it,ik)
       enddo
    enddo
    
    !Make corrections for discrete spectrum
    ebin=ebin*ebincorr

  end subroutine get_polar_spectrum
!================================================================================
! Calculate corrected polar spectrum with an average in z
!================================================================================
! NOTE: polar_index connects a given mode to the correct binned kperp value
! NOTE: At the moment, this is unused in the code
  subroutine get_polar_spectrum_zavg(a,b,ebin)
    use theta_grid, only: ntgrid, delthet, jacob
    use kt_grids, only: naky, ntheta0,aky
    use nonlinear_terms, only: nonlin
    implicit none
    real, dimension (-ntgrid:,:,:), intent (in) :: a,b  !data by (ig,kx,ky) mode
    real, dimension (:), intent (out) :: ebin  !variable ee in kperp bins
    integer :: ik, it
    integer :: ng
    real, dimension (-ntg_out:ntg_out) :: wgt
    real :: anorm
    real :: fac                                !Factor for ky=0 modes

    !Get weighting for z-average
    ng = ntg_out
    wgt = delthet(-ng:ng)*jacob(-ng:ng)
    anorm = sum(wgt)

    !Initialize ebin
!    write (*,*) size(ebin),' is size of ebin'
    ebin=0.

    !Loop through all modes and sum number at each kperp
    do it = 1,ntheta0
       do ik= 1,naky
          fac = 0.5
          if (aky(ik) < epsilon(0.)) fac = 1.0

          if (nonlin .and. it == 1 .and. ik == 1) cycle

!          write (*,*) polar_index(it,ik),' should be < nbx?'
          ebin(polar_index(it,ik))= ebin(polar_index(it,ik)) + &
               sum(real(a(-ng:ng,it,ik)*b(-ng:ng,it,ik)*wgt(-ng:ng)))/anorm*fac
       enddo
    enddo
    
    !Make corrections for discrete spectrum
    ebin=ebin*ebincorr

  end subroutine get_polar_spectrum_zavg
!================================================================================
!
!================================================================================
! Deallocate variables used for polar spectrum
  subroutine finish_polar_spectrum
    implicit none

    !Deallocate variables
    deallocate(polar_index,kpbin,ebincorr,ebinarray,etmp)
    deallocate(numavg,kpavg,eavgarray,polar_avg_index)

  end subroutine finish_polar_spectrum

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
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension (:,:), intent (in) :: a, b
    real, intent (out) :: axb
    real, dimension (:,:), intent (out) :: axb_by_mode

    integer :: ik, it

    do ik = 1, naky
       do it = 1, ntheta0
          axb_by_mode(it,ik) = real(conjg(a(it,ik))*b(it,ik))
       end do
    end do

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

! bug fixed 3.17.00
!    fac = 1.0
!    if (naky > 1) fac(1) = 0.5
! old field array order: 
!    favg = sum(f*spread(fac,2,ntheta0))

!    fac = 0.5
!    fac(1) = 1.0

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
    real, dimension (nensembles,nspec) :: pflx_global, qflx_global, heat_global
    time_int=user_time-start_time
    call scope (allprocs)
    call group_to_all (time_int, dt_global, nensembles)
    call broadcast (dt_global)
    time_int = sum(dt_global)
    call group_to_all (pflux_avg, pflx_global, nensembles)
    call group_to_all (qflux_avg, qflx_global, nensembles)
    do is = 1, nspec
       call broadcast (pflx_global(:,is))
       call broadcast (qflx_global(:,is))
       pflux_avg = sum(pflx_global(:,is))
       qflux_avg = sum(qflx_global(:,is))
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
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension (:,:), intent (in) :: a, b
    complex, intent (out) :: axb
    complex, dimension (:,:), intent (out) :: axb_by_mode

    integer :: ik, it

    do ik = 1, naky
       do it = 1, ntheta0
          axb_by_mode(it,ik) = conjg(a(it,ik))*b(it,ik)
       end do
    end do

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

end module gs2_diagnostics
