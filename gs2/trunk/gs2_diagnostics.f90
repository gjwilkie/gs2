! what about boundary condition contributions?  In the presence 
! of magnetic shear, there are not always zero amplitudes at the ends
! of the supercells.

module gs2_diagnostics

  use gs2_heating, only: heating_diagnostics,dens_vel_diagnostics

  implicit none

  public :: init_gs2_diagnostics
  public :: finish_gs2_diagnostics
  public :: loop_diagnostics

  interface get_vol_average
     module procedure get_vol_average_one, get_vol_average_all
  end interface

  ! knobs
  real :: omegatol, omegatinst
  logical :: print_line, print_old_units, print_flux_line, print_summary
  logical :: write_line, write_flux_line, write_phi, write_apar, write_bpar, write_aperp
  logical :: write_omega, write_omavg, write_ascii, write_lamavg, write_eavg
  logical :: write_qheat, write_pflux, write_vflux, write_tavg
  logical :: write_qmheat, write_pmflux, write_vmflux, write_gs, write_lpoly
  logical :: write_qbheat, write_pbflux, write_vbflux, write_g, write_gg, write_gyx
  logical :: write_dmix, write_kperpnorm, write_phitot, write_epartot
  logical :: write_eigenfunc, write_final_fields, write_final_antot
  logical :: write_final_moments, write_avg_moments, write_stress
  logical :: write_fcheck, write_final_epar, write_kpar
  logical :: write_vortcheck, write_fieldcheck
  logical :: write_fieldline_avg_phi, write_hrate, write_lorentzian
  logical :: write_neoclassical_flux, write_nl_flux, write_Epolar
  logical :: write_verr, write_cerr, write_max_verr
  logical :: exit_when_converged
  logical :: dump_neoclassical_flux, dump_check1, dump_check2
  logical :: dump_fields_periodically, make_movie
  logical :: dump_final_xfields
  logical :: use_shmem_for_xfields
  logical :: save_for_restart
  logical :: test_conserve
!>GGH
  logical, parameter :: write_density_velocity=.false.
  logical :: write_jext=.false.
!<GGH

  integer :: nperiod_output
  integer :: nwrite, igomega, nmovie
  integer :: navg, nsave

  ! internal
  logical :: write_any, write_any_fluxes, dump_any
  logical, private :: initialized = .false.

  integer :: out_unit, kp_unit, heat_unit, polar_raw_unit, polar_avg_unit, heat_unit2, lpc_unit
  integer :: dv_unit, jext_unit   !GGH Additions
  integer :: dump_neoclassical_flux_unit, dump_check1_unit, dump_check2_unit
  integer :: res_unit, res_unit2

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

  real, dimension (:,:,:), allocatable ::  pflux,  vflux
  real, dimension (:,:,:), allocatable :: pmflux, vmflux
  real, dimension (:,:,:), allocatable :: pbflux, vbflux
  real, dimension (:,:),   allocatable :: theta_pflux, theta_pmflux, theta_pbflux
  real, dimension (:,:),   allocatable :: theta_vflux, theta_vmflux, theta_vbflux
  real, dimension (:,:,:), allocatable :: theta_qflux, theta_qmflux, theta_qbflux

  ! (ntheta0,naky,nspec)

  integer :: ntg_out

contains

  subroutine init_gs2_diagnostics (list, nstep)
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids, ntheta0, naky
    use run_parameters, only: init_run_parameters
    use species, only: init_species, nspec
    use dist_fn, only: init_dist_fn
    use gs2_flux, only: init_gs2_flux
    use gs2_io, only: init_gs2_io
    use gs2_heating, only: init_htype,init_dvtype
    use collisions, only: collision_model_switch, init_lorentz_error
    use mp, only: broadcast, proc0
    use le_grids, only: init_weights

    implicit none
    logical, intent (in) :: list
    integer, intent (in) :: nstep
    real :: denom
    integer :: ik, it, nmovie_tot

    if (initialized) return
    initialized = .true.

    call init_theta_grid
    call init_kt_grids
    call init_run_parameters
    call init_species
    call init_dist_fn
    call init_gs2_flux

    call real_init (list)
    call broadcast (navg)
    call broadcast (nwrite)
    call broadcast (nmovie)
    call broadcast (nsave)
    call broadcast (write_any)
    call broadcast (write_any_fluxes)
    call broadcast (write_neoclassical_flux)
    call broadcast (write_nl_flux)
    call broadcast (write_omega)
    call broadcast (write_fieldline_avg_phi)
    call broadcast (write_fcheck)
    call broadcast (dump_any)
    call broadcast (dump_neoclassical_flux)
    call broadcast (dump_check1)
    call broadcast (dump_check2)
    call broadcast (dump_fields_periodically)
    call broadcast (make_movie)
    call broadcast (dump_final_xfields)
    call broadcast (use_shmem_for_xfields)
    call broadcast (save_for_restart)
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
         write_eigenfunc, make_movie, nmovie_tot, write_verr)
    
    if (write_cerr) then
       if (collision_model_switch == 1 .or. collision_model_switch == 5) then
          call init_lorentz_error
       else
          write_cerr = .false.
       end if
    end if

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
    integer :: ig, ik, it

    call read_parameters (list)
    if (proc0) then
       if (write_ascii) then
          call open_output_file (out_unit, ".out")
!          if (write_kpar) call open_output_file (kp_unit, ".kp")
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

    allocate (pflux (ntheta0,naky,nspec)) ; pflux = 0.
    allocate (qheat (ntheta0,naky,nspec,3)) ; qheat = 0.
!       allocate (qheat_par  (ntheta0,naky,nspec))
!       allocate (qheat_perp (ntheta0,naky,nspec))
    allocate (vflux (ntheta0,naky,nspec)) ; vflux = 0.
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
    allocate (theta_qflux (-ntgrid:ntgrid, nspec, 3)) ; theta_qflux = 0.
    
    allocate (theta_pmflux (-ntgrid:ntgrid, nspec)) ; theta_pmflux = 0.
    allocate (theta_vmflux (-ntgrid:ntgrid, nspec)) ; theta_vmflux = 0.
    allocate (theta_qmflux (-ntgrid:ntgrid, nspec, 3)) ; theta_qmflux = 0.
    
    allocate (theta_pbflux (-ntgrid:ntgrid, nspec)) ; theta_pbflux = 0.
    allocate (theta_vbflux (-ntgrid:ntgrid, nspec)) ; theta_vbflux = 0.
    allocate (theta_qbflux (-ntgrid:ntgrid, nspec, 3)) ; theta_qbflux = 0.

  end subroutine real_init

  subroutine read_parameters (list)
    use file_utils, only: input_unit, run_name, input_unit_exist
    use theta_grid, only: nperiod, ntheta
    use gs2_flux, only: gs2_flux_adjust
    use dist_fn, only: nperiod_guard
    use kt_grids, only: box, nx, ny
    use mp, only: proc0
    implicit none
    integer :: in_file
    logical, intent (in) :: list
    logical :: exist
    namelist /gs2_diagnostics_knobs/ &
         print_line, print_old_units, print_flux_line, &
         write_line, write_flux_line, write_phi, write_apar, write_bpar, write_aperp, &
         write_omega, write_omavg, write_ascii, write_kpar, write_lamavg, &
         write_qheat, write_pflux, write_vflux, write_eavg, write_gs, write_gyx, &
         write_qmheat, write_pmflux, write_vmflux, write_tavg, write_g, write_gg, &
         write_qbheat, write_pbflux, write_vbflux, write_hrate, write_lpoly, &
         write_dmix, write_kperpnorm, write_phitot, write_epartot, &
         write_eigenfunc, write_final_fields, write_final_antot, &
         write_fcheck, write_final_epar, write_final_moments, write_cerr, &
         write_vortcheck, write_fieldcheck, write_Epolar, write_verr, write_max_verr, &
         write_fieldline_avg_phi, write_neoclassical_flux, write_nl_flux, &
         nwrite, nmovie, nsave, navg, omegatol, omegatinst, igomega, write_lorentzian, &
         exit_when_converged, write_avg_moments, write_stress, &
         dump_neoclassical_flux, dump_check1, dump_check2, &
         dump_fields_periodically, make_movie, &
         dump_final_xfields, use_shmem_for_xfields, &
         nperiod_output, test_conserve, &
         save_for_restart

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
       in_file = input_unit_exist ("gs2_diagnostics_knobs", exist)
       if (exist) read (unit=input_unit("gs2_diagnostics_knobs"), nml=gs2_diagnostics_knobs)

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
       write_any_fluxes =  write_flux_line .or. print_flux_line .or. write_nl_flux &
            .or. gs2_flux_adjust
       dump_any = dump_neoclassical_flux .or. dump_check1  .or. dump_fields_periodically &
            .or. dump_check2 .or. make_movie .or. print_summary

       nperiod_output = min(nperiod,nperiod_output)
       ntg_out = ntheta/2 + (nperiod_output-1)*ntheta
    end if 
 end subroutine read_parameters

  subroutine finish_gs2_diagnostics (istep)
    use file_utils, only: open_output_file, close_output_file, get_unused_unit
    use mp, only: proc0, broadcast, nproc, iproc, sum_reduce, barrier
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
!    use gs2_dist_io, only: write_dist
    implicit none
    integer, intent (in) :: istep
    integer :: ig, ik, it, il, ie, is, unit, ierr
    complex, dimension (:,:,:,:), allocatable :: fcheck_f
!    complex, dimension (:,:,:), allocatable :: xphi, xapar, xbpar
    real, dimension (:), allocatable :: total
    real, dimension (:,:,:), allocatable :: xphi, xapar, xbpar
    real, dimension (:,:,:), allocatable :: bxf, byf, vxf, vyf, bxfsavg, byfsavg
    real, dimension (:,:,:), allocatable :: bxfs, byfs, vxfs, vyfs, rvx, rvy, rx, ry
    complex, dimension (:,:,:), allocatable :: bx, by, vx, vy, vx2, vy2
    complex, dimension (:,:,:), allocatable :: phi2, apar2, bpar2, antot, antota, antotp, epar
    complex, dimension (:,:,:,:), allocatable :: ntot, density, upar, tpar, tperp
    real, dimension (:), allocatable :: dl_over_b
    real, dimension (:,:,:), allocatable :: lamflux, enflux, tflux
!    real, dimension (-ntgrid:ntgrid, nspec, 4) :: tflux
    complex, dimension (ntheta0, naky) :: phi0
    real, dimension (ntheta0, naky) :: phi02
    real, dimension (nspec) :: weights
    real, dimension (2*ntgrid) :: kpar
    real, dimension (:), allocatable :: xx4, yy4, dz
    real, dimension (:,:), allocatable :: bxs, bys, vxs, vys
    real, dimension (:,:), allocatable :: bxsavg, bysavg
    real, dimension (:), allocatable :: stemp, zx1, zxm, zy1, zyn, xx, yy
    real :: zxy11, zxym1, zxy1n, zxymn, L_x, L_y, rxt, ryt, bxt, byt
    real :: geavg, glavg
    integer :: istatus, nnx, nny, nnx4, nny4, ulim, llim, iblock, i, g_unit
    logical :: last = .true.

    real, dimension (:,:,:,:), allocatable :: errest_by_mode
    integer, dimension (:,:), allocatable :: erridx
    real, dimension (:,:), allocatable :: errest

    if (write_gyx) call write_fyx (phinew,bparnew,last)
    if (write_g) call write_f (last)
    if (write_lpoly) call write_poly (phinew, bparnew, last, istep)
    if (write_cerr) call collision_error (phinew, bparnew, last, istep)
!    if (write_gg) call write_dist (g)

    phi0 = 1.

    if (proc0) then
       if (write_ascii .and. write_verr) then
          call close_output_file (res_unit)
          call close_output_file (lpc_unit)
          if (write_max_verr) call close_output_file (res_unit2)
       end if
       if (write_ascii) call close_output_file (out_unit)
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
       call getmoms (phinew, ntot, density, upar, tpar, tperp)

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
!                              real(density(ig,it,ik,is)/phi0(it,ik)), &
!                              real(tpar(ig,it,ik,is)/phi0(it,ik))
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
    end if

1000  format(20(1x,1pg18.11))
  end subroutine finish_gs2_diagnostics

  subroutine loop_diagnostics (istep, exit)
    use species, only: nspec, spec
    use theta_grid, only: theta, ntgrid, delthet, jacob
    use theta_grid, only: gradpar, nperiod
    use kt_grids, only: naky, ntheta0, aky_out, theta0, akx_out, aky
    use kt_grids, only: nkpolar !, akpolar, akpolar_out
    use run_parameters, only: woutunits, funits, fapar, fphi, fbpar
    use fields, only: phinew, aparnew, bparnew
    use fields, only: kperp, fieldlineavgphi, phinorm
    use dist_fn, only: flux, vortcheck, fieldcheck, get_stress, write_f, write_fyx
    use dist_fn, only: neoclassical_flux, omega0, gamma0, getmoms, par_spectrum, gettotmoms
    use dist_fn, only: get_verr, get_gtran, write_poly, collision_error, neoflux
    use collisions, only: ncheck, vnmult, vary_vnew
!    use collisions, only: ntot_diff, upar_diff, uperp_diff, ttot_diff
    use mp, only: proc0, broadcast, iproc
    use file_utils, only: get_unused_unit, flush_output_file
    use prof, only: prof_entering, prof_leaving
    use gs2_time, only: update_time, user_time
    use gs2_io, only: nc_qflux, nc_vflux, nc_pflux, nc_loop, nc_loop_moments
    use gs2_io, only: nc_loop_stress, nc_loop_vres
    use gs2_io, only: nc_loop_movie
    use gs2_layouts, only: yxf_lo
    use gs2_transforms, only: init_transforms, transform2
    use le_grids, only: nlambda
    use nonlinear_terms, only: nonlin
    use antenna, only: antenna_w
    use gs2_flux, only: check_flux
    use gs2_heating, only: heating_diagnostics, init_htype, del_htype, &
         dens_vel_diagnostics,init_dvtype, del_dvtype
    use constants
    implicit none
    integer :: nout = 1
    integer :: nout_movie = 1
    integer, intent (in) :: istep
    logical, intent (out) :: exit
    logical :: accelerated
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
    real, dimension (:,:,:,:), allocatable :: errest_by_mode
    integer, dimension (:,:), allocatable :: erridx
    real, dimension (:,:), allocatable :: errest
    real :: geavg, glavg
    real :: dmix, dmix4, dmixx
    real :: t, denom
    integer :: ig, ik, it, is, unit, il, i, j, nnx, nny, ifield, write_mod
    complex :: phiavg, sourcefac
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky,nspec) :: ntot, density, &
         upar, tpar, tperp, upartot, uperptot, ttot
    complex, dimension (ntheta0, nspec) :: ntot00, density00, upar00, tpar00, tperp00
    complex, dimension (ntheta0, nspec) :: upartot00, ttot00
    complex, dimension (ntheta0, nspec) :: rstress, ustress
    complex, dimension (ntheta0) :: phi00
    complex, allocatable, dimension (:,:,:) :: phik2
    complex, save :: wtmp_new
    complex :: wtmp_old = 0.
    real, dimension (:), allocatable :: dl_over_b
    real, dimension (2*ntgrid) :: kpar
    real, dimension (ntheta0, nspec) :: x_qmflux
    real, dimension (nspec) ::  heat_fluxes,  part_fluxes, mom_fluxes,  ntot2, ntot20
    real, dimension (nspec) :: mheat_fluxes, mpart_fluxes, mmom_fluxes
    real, dimension (nspec) :: bheat_fluxes, bpart_fluxes, bmom_fluxes
    real, dimension (nspec) ::  heat_par,  heat_perp
    real, dimension (nspec) :: mheat_par, mheat_perp
    real, dimension (nspec) :: bheat_par, bheat_perp
    real, dimension (naky) :: fluxfac
!    real, dimension (:), allocatable :: phi_by_k, apar_by_k, bpar_by_k
    real :: hflux_tot, zflux_tot, vflux_tot
    real, dimension(nspec) :: tprim_tot, fprim_tot
    character(200) :: filename
    logical :: last = .false.

    complex, dimension (nspec) :: ntdiff, upadiff, upediff, ttdiff

    call prof_entering ("loop_diagnostics")

    exit = .false.

    if (proc0) call get_omegaavg (istep, exit, omegaavg)
    call broadcast (exit)

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

    call update_time

    if (make_movie .and. mod(istep,nmovie)==0) then
       t = user_time
       ! EAB 09/17/03 -- modify dump_fields_periodically to print out inverse fft of fields in x,y
       ! write(*,*) "iproc now in dump_fields_periodically case", iproc 
       nnx = yxf_lo%nx
       nny = yxf_lo%ny
       if (fphi > epsilon(0.0)) then
          allocate (yxphi(nnx,nny,-ntgrid:ntgrid))
          call getmoms (phinew, ntot, density, upar, tpar, tperp)
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
       call get_gtran (geavg, glavg, phinew, bparnew, istep)
       if (proc0) then

! write error estimates to .nc file          
!          call nc_loop_vres (nout, errest_by_mode, lpcoef_by_mode)

! write error estimates for ion dist. fn. at outboard midplane with ik=it=1 to ascii files
          if (write_ascii) then
             t = user_time
             write(lpc_unit,"(3(1x,e12.6))") t, geavg, glavg             
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
       call flux (phinew, aparnew, bparnew, &
             pflux,  qheat,  vflux, &
            pmflux, qmheat, vmflux, &
            pbflux, qbheat, vbflux, &
       theta_pflux, theta_vflux, theta_qflux, &
       theta_pmflux, theta_vmflux, theta_qmflux, & 
       theta_pbflux, theta_vbflux, theta_qbflux)
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
                
                vflux(:,:,is) = vflux(:,:,is)*funits &
                     *spec(is)%dens*spec(is)%mass*spec(is)%stm
                call get_volume_average (vflux(:,:,is), mom_fluxes(is))

                theta_vflux(:,is) = theta_vflux(:,is)*funits &
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
       end if
    end if

    fluxfac = 0.5
    fluxfac(1) = 1.0

    if (proc0) then
       if (print_flux_line) then
          write (unit=*, fmt="('t= ',e16.10,' <phi**2>= ',e10.4, &
               & ' heat fluxes: ', 5(1x,e10.4))") &
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
    call check_flux (i, t, heat_fluxes)
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

    call prof_leaving ("loop_diagnostics-1")

    if (.not. (write_any .or. dump_any)) return

    call prof_entering ("loop_diagnostics-2")

    call kperp (ntg_out, akperp)

    if (write_neoclassical_flux .or. dump_neoclassical_flux .and. neoflux) then
       call neoclassical_flux (pfluxneo, qfluxneo, phinew, bparnew)
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
                     & ' heat fluxes: ', 5(1x,e10.4))") &
                     t, phi2, heat_fluxes(1:min(nspec,5))
                write (unit=out_unit, fmt="('t= ',e16.10,' <phi**2>= ',e10.4, &
                     & ' part fluxes: ', 5(1x,e10.4))") &
                     t, phi2, part_fluxes(1:min(nspec,5))
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
                  mom_fluxes, mmom_fluxes, bmom_fluxes, vflux_tot)
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
                   write (out_unit, "('t= ',e16.10,' aky= ',f5.2, ' akx= ',f5.2, &
                        &' om= ',2f8.3,' omav= ', 2f8.3,' phtot= ',e8.2)") &
                        t, aky_out(ik), akx_out(it), &
                        real( omega(it,ik)*woutunits(ik)), &
                        aimag(omega(it,ik)*woutunits(ik)), &
                        real( omegaavg(it,ik)*woutunits(ik)), &
                        aimag(omegaavg(it,ik)*woutunits(ik)), &
                        phitot(it,ik)
                end if
                
                if (write_omega) write (out_unit, *) &
                     'omega=', omega(it,ik), &
                     ' omega/(vt/a)=', real(omega(it,ik)*woutunits(ik)), &
                     aimag(omega(it,ik)*woutunits(ik))
                if (write_omavg) write (out_unit, *) &
                     'omavg=', omegaavg(it,ik), &
                     ' omavg/(vt/a)=', real(omegaavg(it,ik)*woutunits(ik)), &
                     aimag(omegaavg(it,ik)*woutunits(ik))
                
                if (write_dmix) then
                   if (abs(akperp(it,ik)) < 2.0*epsilon(0.0)) then
                      write (out_unit, *) 'dmix indeterminate'
                   else
                      dmix = 2.0*woutunits(ik)*aimag(omega(it,ik))/akperp(it,ik)**2
                      dmix4 = dmix*dmix/(woutunits(ik)*aimag(omega(it,ik)))
                      dmixx = 2.0*woutunits(ik)*aimag(omega(it,ik))/(akperp(it,ik)**2-aky(ik)**2)
                      write (out_unit, *) 'dmix=', dmix,' dmix4= ',dmix4,' dmixx= ',dmixx
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

    if (write_cerr) call collision_error(phinew,bparnew,last,istep)

    if (write_stress) then
       call get_stress (rstress, ustress)
       if (proc0) call nc_loop_stress(nout, rstress, ustress)
    end if

    call broadcast (test_conserve)
    call broadcast (write_avg_moments)
    if (write_avg_moments) then

       call getmoms (phinew, ntot, density, upar, tpar, tperp)

       if (test_conserve) call gettotmoms (phinew, bparnew, ntot, upartot, uperptot, ttot)

       if (proc0) then
          allocate (dl_over_b(-ntg_out:ntg_out))

          dl_over_b = delthet(-ntg_out:ntg_out)*jacob(-ntg_out:ntg_out)
          dl_over_b = dl_over_b / sum(dl_over_b)

          do is = 1, nspec
             do it = 1, ntheta0
                ntot00(it, is)   = sum( ntot(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                density00(it,is) = sum(density(-ntg_out:ntg_out,it,1,is)*dl_over_b)
!                upar00(it, is)   = sum( upar(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                upar00(it, is)   = sum( upartot(-ntg_out:ntg_out,it,1,is)*dl_over_b)
!                tpar00(it, is)   = sum( tpar(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                tpar00(it, is)   = sum( ttot(-ntg_out:ntg_out,it,1,is)*dl_over_b)
!                tperp00(it, is)  = sum(tperp(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                tperp00(it, is)  = sum(uperptot(-ntg_out:ntg_out,it,1,is)*dl_over_b)

!                if (test_conserve .and. allocated(ntot_diff)) then
!                   ntdiff(is) = sum(ntot_diff(-ntg_out:ntg_out,1,1,is)*dl_over_b)
!                   upadiff(is) = sum(upar_diff(-ntg_out:ntg_out,1,1,is)*dl_over_b)
!                   upediff(is) = sum(uperp_diff(-ntg_out:ntg_out,1,1,is)*dl_over_b)
!                   ttdiff(is) = sum(ttot_diff(-ntg_out:ntg_out,1,1,is)*dl_over_b)
!                end if
             end do
          end do
          
!          if (test_conserve .and. allocated(ntot_diff)) then
          if (test_conserve) then
             if (nspec == 1) then
                write (out_unit, *) 'moms', real(ntot00(1,1)), aimag(upar00(1,1)), real(tperp00(1,1)), real(tpar00(1,1))
!                write (out_unit, *) 'moms', real(ntdiff(1)), aimag(ntdiff(1)), &
!                     real(upadiff(1)), aimag(upadiff(1)), real(upediff(1)), aimag(upediff(1)), &
!                     real(ttdiff(1)), aimag(ttdiff(1))
             else
!                write (out_unit, *) 'moms', real(ntdiff(1)), real(ntdiff(2)), &
!                     aimag(upadiff(1)), aimag(upadiff(2)), real(ttdiff(1)), real(ttdiff(2))
             end if
          end if

          do it = 1, ntheta0
             phi00(it)   = sum( phinew(-ntg_out:ntg_out,it,1)*dl_over_b)
          end do

          deallocate (dl_over_b)

          do is = 1, nspec
             call get_vol_average (ntot(:,:,:,is), ntot(:,:,:,is), &
                  ntot2(is), ntot2_by_mode(:,:,is))
          end do

          do is = 1, nspec
             call get_vol_average (ntot(0,:,:,is), ntot(0,:,:,is), &
                  ntot20(is), ntot20_by_mode(:,:,is))
          end do

!          write (out_unit, "('t,u0 ',2(1x,e12.6))") t, aimag(upar00 (1,1))

          call nc_loop_moments (nout, ntot2, ntot2_by_mode, ntot20, &
               ntot20_by_mode, phi00, ntot00, density00, upar00, &
               tpar00, tperp00)

       end if
    end if

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
    if (write_ascii .and. mod(nout, 10) == 0 .and. proc0) &
         call flush_output_file (out_unit, ".out")

    if (write_ascii .and. write_verr .and. mod(nout, 10) == 0 .and. proc0) &
         call flush_output_file (res_unit, ".vres")

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
    call get_dens_vel(dv, dvk, phi, apar, bpar, phinew, aparnew, bparnew)    
    
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

  subroutine get_omegaavg (istep, exit, omegaavg)
    use kt_grids, only: naky, ntheta0
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use gs2_time, only: code_dt
    use constants
    implicit none
    integer, intent (in) :: istep
    logical, intent (in out) :: exit
    complex, dimension (:,:), intent (out) :: omegaavg
    complex, dimension (navg,ntheta0,naky) :: domega
    integer :: j

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
1    real, dimension(:), allocatable :: kpavg_lim
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

    integer :: ig, ik, it
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

end module gs2_diagnostics
