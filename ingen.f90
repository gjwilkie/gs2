program ingen

! 
! Reads in namelists, checks for consistency of inputs, writes a report
! to runname.report, checks for requested variations, 
! generates new input files, and exits.
! 
! Consistency checks/reports:
!
! wstar_units incompatible with nonlinear
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!        Declarations              !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  use mp, only: init_mp, finish_mp
  use constants 
  use file_utils
  use text_options, only: text_option, get_option_value
  implicit none
  
  character (100) :: pythonin
  integer :: interactive_record, interactive_input

  integer :: nlambda
  integer :: in_file, i, ierr, unit, is, report_unit, iunit, ncut
  logical :: exist, shat_is_known = .true., scan, stdin
  logical :: has_electrons = .false.

  integer, dimension (18), parameter :: nesub_ok = (/ &
       1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20, 24, &
       32, 48, 64 /)
    
  integer, dimension (12), parameter :: nesuper_ok = (/ &
       1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15 /)

  integer, dimension (15), parameter :: ngauss_ok = (/ &
       1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48 /)

! gs2_layouts:
  character (len=5) :: layout
  logical :: local_field_solve

! antenna: 
  complex :: w_antenna, a, b
  real :: amplitude, t0_driver, w_dot
  integer :: nk_stir, kx, ky, kz
  logical :: write_antenna, no_driver = .false., ant_off = .false., travel
  complex, dimension (:), allocatable :: a_ant, b_ant
  integer, dimension (:), allocatable :: kx_stir, ky_stir, kz_stir
  logical, dimension (:), allocatable :: trav

  real :: cfac

! collisions: 
  integer :: collision_model_switch
  real :: vncoef, absom, cfac_nu
  integer :: ivnew
  character (20) :: collision_model
  logical :: conserve_number, conserve_momentum, use_shmem, hypercoll, const_v
  logical :: heating
  integer, parameter :: collision_model_lorentz = 1
  integer, parameter :: collision_model_krook = 2
  integer, parameter :: collision_model_none = 3
  integer, parameter :: collision_model_krook_test = 4
  integer, parameter :: collision_model_lorentz_test = 5
  type (text_option), dimension (7), parameter :: coll_modelopts = &
       (/ text_option('default', collision_model_lorentz), &
       text_option('lorentz', collision_model_lorentz), &
       text_option('krook', collision_model_krook), &
       text_option('krook-test', collision_model_krook_test), &
       text_option('lorentz-test', collision_model_lorentz_test), &
       text_option('none', collision_model_none), &
       text_option('collisionless', collision_model_none) /)

! init_g:
  integer :: ginitopt_switch
  integer, parameter :: ginitopt_default = 1, ginitopt_test1 = 2, &
       ginitopt_xi = 3, ginitopt_xi2 = 4, ginitopt_rh = 5, ginitopt_zero = 6, &
       ginitopt_test3 = 7, ginitopt_convect = 8, ginitopt_restart_file = 9, &
       ginitopt_noise = 10, ginitopt_restart_many = 11, ginitopt_continue = 12, &
       ginitopt_nl = 13, ginitopt_kz0 = 14, ginitopt_restart_small = 15, &
       ginitopt_nl2 = 16, ginitopt_nl3 = 17, ginitopt_nl4 = 18, &
       ginitopt_nl5 = 19, ginitopt_alf = 20, ginitopt_kpar = 21, &
       ginitopt_nl6 = 22, ginitopt_nl7 = 23, ginitopt_gs = 24, ginitopt_recon = 25, &
       ginitopt_nl3r = 26
  character (20) :: ginit_option
  type (text_option), dimension (26), parameter :: ginitopts = &
       (/ text_option('default', ginitopt_default), &
       text_option('noise', ginitopt_noise), &
       text_option('test1', ginitopt_test1), &
       text_option('xi', ginitopt_xi), &
       text_option('xi2', ginitopt_xi2), &
       text_option('zero', ginitopt_zero), &
       text_option('test3', ginitopt_test3), &
       text_option('convect', ginitopt_convect), &
       text_option('rh', ginitopt_rh), &
       text_option('many', ginitopt_restart_many), &
       text_option('small', ginitopt_restart_small), &
       text_option('file', ginitopt_restart_file), &
       text_option('cont', ginitopt_continue), &
       text_option('kz0', ginitopt_kz0), &
       text_option('nl', ginitopt_nl), &
       text_option('nl2', ginitopt_nl2), &
       text_option('nl3', ginitopt_nl3), &
       text_option('nl3r', ginitopt_nl3r), &
       text_option('nl4', ginitopt_nl4), &
       text_option('nl5', ginitopt_nl5), &
       text_option('nl6', ginitopt_nl6), &
       text_option('nl7', ginitopt_nl7), &
       text_option('alf', ginitopt_alf), &
       text_option('gs', ginitopt_gs), &
       text_option('recon', ginitopt_recon), &
       text_option('kpar', ginitopt_kpar) /)
  real :: initk0
  real :: width0, phiinit, k0, imfac, refac, zf_init
  real :: den0, upar0, tpar0, tperp0
  real :: den1, upar1, tpar1, tperp1
  real :: den2, upar2, tpar2, tperp2
  real :: tstart, scale, apar0, dphiinit
  logical :: chop_side, left
  character(300) :: restart_file
  integer, dimension(2) :: ikk, itt

! dist_fn:
  complex, dimension (:), allocatable :: fexp ! (nspec)
  real, dimension (:), allocatable :: bkdiff  ! (nspec)
  integer, dimension (:), allocatable :: bd_exp ! nspec
  real :: gridfac, apfac, driftknob, tpdriftknob, poisfac
  real :: kfilter, afilter, D_kill, noise, g_exb, omprimfac
  real :: t0, omega0, gamma0, source0, thetas, phi_ext, a_ext
  real :: akx_star, aky_star, cfac_df
  integer :: nperiod_guard 
  logical :: mult_imp, test, def_parity, even, save_n, save_u, save_Tpar
  logical :: save_Tperp, test_df
  character (20) :: source_option
  character (20) :: boundary_option
  character (30) :: adiabatic_option
  character (20) :: heating_option
  integer :: adiabatic_option_switch, heating_option_switch
  integer, parameter :: adiabatic_option_default = 1, &
       adiabatic_option_zero = 2, &
       adiabatic_option_fieldlineavg = 3, &
       adiabatic_option_yavg = 4, &
       adiabatic_option_noJ = 5
  integer :: source_option_switch
  integer, parameter :: source_option_full = 1, &
       source_option_zero = 2, source_option_sine = 3, &
       source_option_test1 = 4, source_option_phiext_full = 5, &
       source_option_test2_full = 6, source_option_cosine = 7, &
       source_option_convect_full = 8
  integer :: boundary_option_switch
  integer, parameter :: boundary_option_zero = 1, &
       boundary_option_self_periodic = 2, &
       boundary_option_alternate_zero = 3, &
       boundary_option_linked = 4
  integer, parameter :: heating_option_hammett = 1, &
       heating_option_cowley = 2

  type (text_option), dimension (9), parameter :: sourceopts = &
       (/ text_option('default', source_option_full), &
       text_option('full', source_option_full), &
       text_option('zero', source_option_zero), &
       text_option('sine', source_option_sine), &
       text_option('cosine', source_option_cosine), &
       text_option('test1', source_option_test1), &
       text_option('phiext_full', source_option_phiext_full), &
       text_option('test2_full', source_option_test2_full), &
       text_option('convect_full', source_option_convect_full) /)

  type (text_option), dimension (8), parameter :: boundaryopts = &
       (/ text_option('default', boundary_option_zero), &
       text_option('zero', boundary_option_zero), &
       text_option('unconnected', boundary_option_zero), &
       text_option('self-periodic', boundary_option_self_periodic), &
       text_option('periodic', boundary_option_self_periodic), &
       text_option('kperiod=1', boundary_option_self_periodic), &
       text_option('linked', boundary_option_linked), &
       text_option('alternate-zero', boundary_option_alternate_zero) /)

  type (text_option), dimension (8), parameter :: adiabaticopts = &
       (/ text_option('default', adiabatic_option_default), &
       text_option('no-field-line-average-term', adiabatic_option_default), &
       text_option('field-line-average-term', adiabatic_option_fieldlineavg), &
       ! eventually add in iphi00 = 0 option:
       text_option('iphi00=0', adiabatic_option_default), &
       text_option('iphi00=1', adiabatic_option_default), &
       text_option('iphi00=2', adiabatic_option_fieldlineavg), &
       text_option('iphi00=3', adiabatic_option_yavg), &
       text_option('dimits', adiabatic_option_noJ) /)

    type (text_option), dimension (3), parameter :: heatingopts = &
         (/ text_option('default', heating_option_cowley), &
            text_option('Hammett', heating_option_hammett), &
            text_option('Cowley',  heating_option_cowley) /)

! fields: 
  character (20) :: field_option
  integer :: fieldopt_switch
  integer, parameter :: fieldopt_implicit = 1, fieldopt_test = 2, fieldopt_explicit = 3
  type (text_option), dimension (4), parameter :: fieldopts = &
       (/ text_option('default', fieldopt_implicit), &
       text_option('implicit', fieldopt_implicit), &
       text_option('explicit', fieldopt_explicit), &
       text_option('test', fieldopt_test) /)

! gs2_diagnostics: 
  logical :: print_line, print_old_units, print_flux_line
  logical :: write_line, write_flux_line, write_phi, write_apar, write_aperp
  logical :: write_omega, write_omavg, write_ascii, write_lamavg
  logical :: write_qheat, write_pflux, write_vflux, write_eavg
  logical :: write_qmheat, write_pmflux, write_vmflux, write_tavg
  logical :: write_qbheat, write_pbflux, write_vbflux, write_gs
  logical :: write_dmix, write_kperpnorm, write_phitot, write_epartot
  logical :: write_eigenfunc, write_final_fields, write_final_antot
  logical :: write_final_moments, write_avg_moments
  logical :: write_fcheck, write_final_epar, write_kpar
  logical :: write_intcheck, write_vortcheck, write_fieldcheck
  logical :: write_fieldline_avg_phi, write_hrate, write_lorentzian
  logical :: write_neoclassical_flux, write_nl_flux, write_Epolar
  logical :: exit_when_converged
  logical :: dump_neoclassical_flux, dump_check1, dump_check2
  logical :: dump_fields_periodically, make_movie
  logical :: dump_final_xfields
  logical :: use_shmem_for_xfields
  logical :: save_for_restart

  integer :: nwrite, igomega, nmovie
  integer :: navg, nperiod_output, nsave

  real :: omegatol, omegatinst

! gs2_reinit: 
  real :: delt_adj, delt_minimum

! hyper: 
  character (9) :: hyper_option
  logical :: const_amp, include_kpar, isotropic_shear
  real :: D_hypervisc, D_hyperres, omega_osc, D_hyper
  integer :: hyper_option_switch
  integer, parameter :: hyper_option_none = 1, &
       hyper_option_visc = 2, &
       hyper_option_res  = 3, &
       hyper_option_both = 4
  type (text_option), dimension(5), parameter :: hyperopts = &
       (/ text_option('default', hyper_option_none), &
       text_option('none', hyper_option_none), &
       text_option('visc_only', hyper_option_visc), &
       text_option('res_only', hyper_option_res), &
       text_option('both', hyper_option_both) /)
  logical :: hyper_on = .false.
  logical :: gridnorm

! kt_grids:
  real, dimension (:), allocatable :: aky_tmp, theta0_tmp, akx_tmp
  real :: aky, theta0, akx
  integer :: naky, ntheta0, nx, ny, jtwist
  integer :: gridopt_switch, normopt_switch
  real :: aky_min, aky_max, theta0_min, theta0_max, lx, ly, x0, y0, rtwist
  real :: y0_internal
  integer, parameter :: gridopt_single = 1, gridopt_range = 2, &
       gridopt_specified = 3, gridopt_box = 4, gridopt_xbox = 5
  type (text_option), dimension (7), parameter :: gridopts = &
       (/ text_option('default', gridopt_single), &
       text_option('single', gridopt_single), &
       text_option('range', gridopt_range), &
       text_option('specified', gridopt_specified), &
       text_option('box', gridopt_box), &
       text_option('nonlinear', gridopt_box), &
       text_option('xbox', gridopt_xbox) /)
  character (20) :: grid_option

  integer, parameter :: normopt_mtk = 1, normopt_bd = 2
  type (text_option), dimension(6), parameter :: normopts = &
       (/ text_option('default', normopt_mtk), &
       text_option('with_root_2', normopt_mtk), &
       text_option('mtk', normopt_mtk), &
       text_option('no_root_2', normopt_bd), &
       text_option('bd', normopt_bd), &
       text_option('t_over_m', normopt_bd) /)
  character (20) :: norm_option
  
  integer :: ngauss, negrid, nesuper, nesub
  real :: ecut, bouncefuzz
  logical :: trapped_particles = .true.
  logical :: advanced_egrid = .true.

! nonlinear_terms: 
  integer :: nonlinear_mode_switch
  integer :: flow_mode_switch
  integer, parameter :: nonlinear_mode_none = 1, nonlinear_mode_on = 2
  integer, parameter :: flow_mode_off = 1, flow_mode_on = 2
  type (text_option), dimension (4), parameter :: nonlinearopts = &
       (/ text_option('default', nonlinear_mode_none), &
       text_option('none', nonlinear_mode_none), &
       text_option('off', nonlinear_mode_none), &
       text_option('on', nonlinear_mode_on) /)
  character (20) :: nonlinear_mode
  type (text_option), dimension (3), parameter :: flowopts = &
       (/ text_option('default', flow_mode_off), &
       text_option('off', flow_mode_off), &
       text_option('on', flow_mode_on) /)
  character (20) :: flow_mode
  real :: cfl, c_par, C_perp, p_x, p_y, p_z
  logical :: zip

! run_parameters:
  real :: beta, zeff, tite, rhostar, teti
  real :: fphi, fapar, fbpar, faperp
  real :: delt, margin
  integer :: nstep
  logical :: wstar_units, eqzip, secondary, tertiary
  integer :: delt_option_switch
  integer, parameter :: delt_option_hand = 1, delt_option_auto = 2
  type (text_option), dimension (3), parameter :: deltopts = &
       (/ text_option('default', delt_option_hand), &
       text_option('set_by_hand', delt_option_hand), &
       text_option('check_restart', delt_option_auto) /)
  character (20) :: delt_option
  
! species:
  integer :: nspec
  real :: z, mass, dens, dens0, temp, tprim, fprim, uprim, uprim2, vnewk, vnewk4
  character (20) :: type

  type :: specie
     real :: z
     real :: mass
     real :: dens, dens0
     real :: temp
     real :: tprim
     real :: fprim
     real :: uprim, uprim2
     real :: vnewk, vnewk4
     real :: stm, zstm, tz, smz, zt
     integer :: type
  end type specie
  type (specie), dimension (:), allocatable :: spec

  integer, parameter :: ion_species = 1
  integer, parameter :: electron_species = 2 ! for collision operator
  integer, parameter :: slowing_down_species = 3 ! slowing-down distn

  type (text_option), dimension (8), parameter :: typeopts = &
       (/ text_option('default', ion_species), &
       text_option('ion', ion_species), &
       text_option('electron', electron_species), &
       text_option('e', electron_species), &
       text_option('beam', slowing_down_species), &
       text_option('fast', slowing_down_species), &
       text_option('alpha', slowing_down_species), &
       text_option('slowing-down', slowing_down_species) /)

! theta_grid: 
  real :: rhoc, rmaj, r_geo, eps, epsl, qinp, shat, alpmhd, pk, shift
  real :: akappa, akappri, tri, tripri, kp
  integer :: ntheta, nperiod, npadd
  real :: alknob, epsknob, bpknob, extrknob, tension, thetamax, deltaw, widthw
  real :: alpmhdfac, alpha1
  character (20) :: model_option
  real :: s_hat_input, alpha_input, invLp_input, beta_prime_input, dp_mult
  real :: delrho, rmin, rmax, ak0, k1, k2
  integer :: itor, iflux, irho, bishop, ismooth, isym
  logical :: ppl_eq, gen_eq, vmom_eq, efit_eq, dfit_eq, equal_arc, idfit_eq
  logical :: local_eq, writelots, gs2d_eq, test_le, transp_eq
  logical :: gb_to_cv
  character (80) :: eqfile
  character (200) :: gridout_file
  character (20) :: equilibrium_option
  integer :: model_switch
  integer, parameter :: model_salpha = 1, model_alpha1 = 2, &
       model_nocurve = 3, model_ccurv = 4, model_b2 = 5, &
       model_eps = 6, model_normal_only = 7
  type (text_option), dimension (8), parameter :: sa_modelopts = &
       (/ text_option('default', model_salpha), &
       text_option('s-alpha', model_salpha), &
       text_option('alpha1', model_alpha1), &
       text_option('rogers', model_eps), &
       text_option('b2', model_b2), &
       text_option('normal_only', model_normal_only), &
       text_option('const-curv', model_ccurv), &
       text_option('no-curvature', model_nocurve) /)

  integer :: eqopt_switch
  integer, parameter :: eqopt_eik = 1, eqopt_salpha = 2, eqopt_file = 3

  type (text_option), dimension (5), parameter :: eqopts = &
       (/ text_option('default', eqopt_eik), &
       text_option('eik', eqopt_eik), &
       text_option('s-alpha', eqopt_salpha), &
       text_option('grid.out', eqopt_file), &
       text_option('file', eqopt_file) /)

  logical :: layouts_write = .false.
  logical :: driver_write = .false.
  logical :: stir_write = .false.
  logical :: collisions_write = .false.
  logical :: init_g_write = .false.
  logical :: dist_fn_write = .false.
  logical :: source_write = .false.
  logical :: fields_write = .false.
  logical :: diagnostics_write = .false.
  logical :: reinit_write = .false.
  logical :: hyper_write = .false.
  logical :: kt_range_write = .false.
  logical :: kt_single_write = .false.
  logical :: kt_specified_write = .false.
  logical :: kt_box_write = .false.
  logical :: kt_xbox_write = .false.
  logical :: kt_write = .false.
  logical :: le_write = .false.
  logical :: nonlinear_write = .false.
  logical :: parameters_write = .false.
  logical :: species_write = .false.
  logical :: species_parameters_write = .false.
  logical :: theta_parameters_write = .false.
  logical :: theta_gridgen_write = .false.
  logical :: theta_salpha_write = .false.
  logical :: theta_eik_write = .false.
  logical :: theta_file_write = .false.
  logical :: theta_write = .false.
  logical :: knobs_write = .false.

  integer :: iostat
  real :: tmpfac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!          Namelists               !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! gs2_layouts:
  namelist /layouts_knobs/ layout, local_field_solve

! antenna: 
  namelist /driver/ amplitude, w_antenna, nk_stir, write_antenna, ant_off, w_dot, t0
  namelist /stir/ kx, ky, kz, travel, a, b

! collisions: 
  namelist /collisions_knobs/ collision_model, vncoef, absom, ivnew, &
       conserve_number, conserve_momentum, use_shmem, hypercoll, cfac, heating

! init_g:
  namelist /init_g_knobs/ ginit_option, width0, phiinit, k0, chop_side, &
       restart_file, left, ikk, itt, scale, tstart, zf_init, &
       den0, upar0, tpar0, tperp0, imfac, refac, even, &
       den1, upar1, tpar1, tperp1, &
       den2, upar2, tpar2, tperp2, dphiinit, apar0

! dist_fn:
  namelist /dist_fn_knobs/ boundary_option, gridfac, apfac, &
       driftknob, tpdriftknob, g_exb, omprimfac, &
       nperiod_guard, poisfac, adiabatic_option, cfac, &
       kfilter, afilter, mult_imp, test, def_parity, even, &
       save_n, save_u, save_Tpar, save_Tperp, D_kill, noise, heating_option
  
  namelist /source_knobs/ t0, omega0, gamma0, source0, &
       thetas, k0, phi_ext, source_option, a_ext, aky_star, akx_star

!  namelist /dist_fn_species_knobs/ fexpr, fexpi, bakdif, bd_exp

! fields: 
  namelist /fields_knobs/ field_option

! gs2_diagnostics: 
  namelist /gs2_diagnostics_knobs/ &
       print_line, print_old_units, print_flux_line, &
       write_line, write_flux_line, write_phi, write_apar, write_aperp, &
       write_omega, write_omavg, write_ascii, write_lamavg, write_tavg, &
       write_qheat, write_pflux, write_vflux, write_kpar, &
       write_qmheat, write_pmflux, write_vmflux, write_eavg, &
       write_qbheat, write_pbflux, write_vbflux, write_hrate, write_lorentzian, &
       write_dmix, write_kperpnorm, write_phitot, write_epartot, &
       write_eigenfunc, write_final_fields, write_final_antot, &
       write_fcheck, write_final_epar, write_final_moments, &
       write_intcheck, write_vortcheck, write_fieldcheck, write_Epolar, &
       write_fieldline_avg_phi, write_neoclassical_flux, write_nl_flux, &
       nwrite, nmovie, nsave, navg, omegatol, omegatinst, igomega, &
       exit_when_converged, write_avg_moments, &
       dump_neoclassical_flux, dump_check1, dump_check2, &
       dump_fields_periodically, make_movie, &
       dump_final_xfields, use_shmem_for_xfields, &
       nperiod_output, &
       save_for_restart, write_gs

! gs2_reinit: 
  namelist /reinit_knobs/ delt_adj, delt_minimum

! hyper: 
  namelist /hyper_knobs/ hyper_option, const_amp, include_kpar, &
       isotropic_shear, D_hyperres, D_hypervisc, omega_osc, D_hyper, gridnorm

! kt_grids:
  namelist /kt_grids_single_parameters/ aky, theta0, akx

  namelist /kt_grids_range_parameters/ naky, ntheta0, &
       aky_min, aky_max, theta0_min, theta0_max

  namelist /kt_grids_specified_parameters/ naky, ntheta0, nx, ny

  namelist /kt_grids_specified_element/ aky, theta0, akx

  namelist /kt_grids_box_parameters/ naky, ntheta0, ly, nx, ny, jtwist, &
       y0, x0, rtwist

  namelist /kt_grids_xbox_parameters/ ntheta0, lx, aky, nx

  namelist /kt_grids_knobs/ grid_option, norm_option

! le_grids:

  namelist /le_grids_knobs/ ngauss, negrid, ecut, bouncefuzz, &
       nesuper, nesub, test, trapped_particles, advanced_egrid

! nonlinear_terms: 
  namelist /nonlinear_terms_knobs/ nonlinear_mode, flow_mode, cfl, &
       C_par, C_perp, p_x, p_y, p_z, zip

! run_parameters:
  namelist /parameters/ beta, zeff, tite, rhostar, teti

  namelist /knobs/ fphi, fapar, fbpar, faperp, delt, nstep, wstar_units, eqzip, &
       delt_option, margin, secondary, tertiary

! species 
  namelist /species_knobs/ nspec
  namelist /species_parameters/ &
       z, mass, dens, dens0, temp, tprim, fprim, uprim, uprim2, vnewk, vnewk4, type

! theta_grid: 
  namelist /theta_grid_parameters/ rhoc, rmaj, r_geo, eps, epsl, &
       qinp, shat, alpmhd, pk, shift, akappa, akappri, tri, tripri, &
       ntheta, nperiod, kp

  namelist /theta_grid_gridgen_knobs/ &
       npadd, alknob, epsknob, bpknob, extrknob, tension, thetamax, deltaw, widthw

  namelist /theta_grid_salpha_knobs/ alpmhdfac, alpha1, model_option

  namelist /theta_grid_eik_knobs/ itor, iflux, irho, &
       ppl_eq, gen_eq, vmom_eq, efit_eq, eqfile, dfit_eq, &
       equal_arc, bishop, local_eq, idfit_eq, gs2d_eq, transp_eq, &
       s_hat_input, alpha_input, invLp_input, beta_prime_input, dp_mult, &
       delrho, rmin, rmax, ismooth, ak0, k1, k2, isym, writelots

  namelist /theta_grid_file_knobs/ gridout_file

  namelist /theta_grid_knobs/ equilibrium_option, gb_to_cv

  namelist /ingen_knobs/ ncut, scan, stdin

!CMR
  logical:: debug=.false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!      Main code starts here       !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call init_mp
 
  if (debug) write(6,*) 'ingen: call get_namelists'
  
  call get_namelists
  if (debug) write(6,*) 'ingen: call report'
  call report
  if (debug) write(6,*) 'ingen: call write_namelists'
  call write_namelists

  if (debug) write(6,*) 'ingen: call interactive, scan=',scan
  if (scan) call interactive

  call finish_mp

contains

  subroutine interactive

    integer :: sel, nbeta, j, ise
    real :: beta_low, beta_high, dbeta, beta_save
    real :: fapar_save, faperp_save, pri, pe, alpi, tpe_save, ptot, alp, dbdr
    real :: alt, aln, fac, beta_prime_save, bishop_save, beta_prime
    real, dimension (:), allocatable :: tp_save, fp_save
    character (500) :: tag, tag1, tag2
    logical :: first = .true.

    if (first) then
       call get_unused_unit (interactive_record)
       open (unit=interactive_record, file='.'//trim(run_name)//".record")
       first = .false.

       if (.not. stdin) then
          call get_unused_unit (interactive_input)
          open (unit=interactive_input, file=trim(pythonin))
       else
          interactive_input = 5
       end if
    end if

    call tell ('Interactive specification of a parameter scan')

100 continue
    
    call text
    call text ('Choose a parameter that you would like to vary (1-6):')
    call text ('(1) beta            (4) temperature gradient')
    call text ('(2) beta_prime      (5) density gradient')
    call text ('(3) collisionality  (6) Z_effective')
    call text
    call get_choice (6, sel)
    
    select case (sel)
       
    case default
       call tell ('Try again.  Choose an integer between 1 and 6, inclusively.')
       goto 100

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  1.0  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    case (1) 
       call tell ('You have chosen to vary beta.')

101    continue

       call text
       call text ('Choose from the following:')
       call text ('(1) Vary beta self-consistently')
       call text ('(2) Vary beta with all other parameters held fixed (not self-consistently).')
       call text
       call get_choice (2, sel)

       select case (sel)

       case default
          call tell ('Try again.  Choose an integer between 1 and 2, inclusively.')
          goto 101

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  1.1  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       case (1) 
          call tell ('You have chosen to vary beta self-consistently.')

102       continue
          call text
          call text ('Choose from the following:')
          call text ('(1) Hold beta_prime fixed, vary electron temperature gradient scale length')
          call text ('(2) Hold beta_prime fixed, vary all temperature gradient scale lengths by same factor')
          call text ('(3) Hold beta_prime fixed, vary all density gradient scale lengths by same factor')
          call text

          call get_choice (3, sel)

          select case (sel)
             
          case default
             call tell ('Try again.  Choose an integer between 1 and 2, inclusively.')
             goto 102
             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  1.1.1  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          case (1)  
             call tell ('You have chosen to vary beta and electron temperature gradient at fixed beta_prime.')

             call beta_range_low (beta_low, 'le', 0.)
             call beta_range_high (beta_high, 'le', beta_low)
             call num_runs (nbeta)

             call tell ('Preparing a self-consistent beta scan at fixed beta_prime.', &
                  'The electron temperature gradient scale length will be varied', &
                  'to maintain consistency.')

             call run_number (sel, nbeta)

             write (tag1, fmt='(i3," runs prepared with beta_min = ",e16.10,&
                  &" and beta_max = ",e16.10)') nbeta, beta_low, beta_high
             write (tag2, fmt='("Files are numbered from ",i3," to ",i3)') sel, sel+nbeta-1

             call tell (tag1, tag2)

             ptot = 0.
             alp = 0.
             pri = 0.
             pe = 0.
             alpi = 0.
             ise = 0
             do is=1,nspec
                if (spec(is)%type == 2) then
                   pe = spec(is)%dens * spec(is)%temp
                   ise = is
                else
                   pri = pri + spec(is)%dens * spec(is)%temp
                   alpi = alpi + spec(is)%dens * spec(is)%temp *(spec(is)%fprim + spec(is)%tprim)
                endif
                ptot = ptot + spec(is)%dens * spec(is)%temp
                alp = alp + spec(is)%dens * spec(is)%temp *(spec(is)%fprim + spec(is)%tprim)
             end do
             
             if (.not. has_electrons) call tell ('You really should use electrons for electromagnetic runs.')

             alp = alp/ptot
             dbdr = - beta*ptot*alp

             dbeta = (beta_high-beta_low)/(nbeta-1)
             do j = sel, sel+nbeta-1
                
                beta_save = beta
                beta = beta_low + (j - sel)*dbeta
                
                tpe_save = spec(ise)%tprim
                spec(ise)%tprim = - (spec(ise)%fprim + alpi/pe + dbdr/beta/pe)

                write (tag1, fmt='("Varying beta and L_Te self-consistently with& 
                     & beta_prime fixed")') 

                write (tag2, fmt='("beta = ",e16.10," and electron tprim = ",e16.10)') beta, spec(ise)%tprim 

                call write_namelists (j, tag1, tag2)
                spec(ise)%tprim = tpe_save
             end do
             beta = beta_save

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  1.1.2  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          case (2)

             call tell ('You have chosen to vary beta and all temperature &
                  &gradient scale lengths together at fixed beta_prime.')

             call beta_range_low (beta_low, 'le', 0.)
             call beta_range_high (beta_high, 'le', beta_low)
             call num_runs (nbeta)

             call tell ('Preparing a self-consistent beta scan at fixed beta_prime.', &
                  'All temperature gradient scale lengths will be varied', &
                  'by the same factor to maintain consistency.')

             call run_number (sel, nbeta)

             write (tag1, fmt='(i3," runs prepared with beta_min = ",e16.10,&
                  &" and beta_max = ",e16.10)') nbeta, beta_low, beta_high
             write (tag2, fmt='("Files are numbered from ",i3," to ",i3)') sel, sel+nbeta-1

             call tell (tag1, tag2)

             allocate (tp_save (nspec))

             ptot = 0.
             alt = 0.
             aln = 0.
             do is=1,nspec
                ptot = ptot + spec(is)%dens * spec(is)%temp
                alt = alt + spec(is)%dens * spec(is)%temp *(spec(is)%tprim)
                aln = aln + spec(is)%dens * spec(is)%temp *(spec(is)%fprim)
             end do
             
             if (.not. has_electrons) call tell ('You really should use electrons for electromagnetic runs.')

             alp = (alt+aln)/ptot
             dbdr = - beta*ptot*alp

             dbeta = (beta_high-beta_low)/(nbeta-1)
             do j = sel, sel+nbeta-1
                
                beta_save = beta
                beta = beta_low + (j - sel)*dbeta
                
                fac = -(dbdr/beta+aln)/alt
                tp_save = spec%tprim
                spec%tprim = fac*spec%tprim

                write (tag1, fmt='("Varying beta and all L_T values self-consistently")')
                write (tag2, fmt='("beta = ",e16.10," and tprim values scaled by ",e16.10)') beta, fac

                call write_namelists (j, tag1, tag2) 
                spec%tprim = tp_save
             end do
             beta = beta_save

             deallocate (tp_save)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  1.1.3  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          case (3)

             call tell ('You have chosen to vary beta and all density &
                  &gradient scale lengths together at fixed beta_prime.')

             call beta_range_low (beta_low, 'le', 0.)
             call beta_range_high (beta_high, 'le', beta_low)
             call num_runs (nbeta)

             call tell ('Preparing a self-consistent beta scan at fixed beta_prime.', &
                  'All density gradient scale lengths will be varied', &
                  'by the same factor to maintain consistency.')

             call run_number (sel, nbeta)

             write (tag1, fmt='(i3," runs prepared with beta_min = ",e16.10,&
                  &" and beta_max = ",e16.10)') nbeta, beta_low, beta_high
             write (tag2, fmt='("Files are numbered from ",i3," to ",i3)') sel, sel+nbeta-1

             call tell (tag1, tag2)

             allocate (fp_save (nspec))

             ptot = 0.
             alt = 0.
             aln = 0.
             do is=1,nspec
                ptot = ptot + spec(is)%dens * spec(is)%temp
                alt = alt + spec(is)%dens * spec(is)%temp *(spec(is)%tprim)
                aln = aln + spec(is)%dens * spec(is)%temp *(spec(is)%fprim)
             end do
             
             if (.not. has_electrons) call tell ('You really should use electrons for electromagnetic runs.')

             alp = (alt+aln)/ptot
             dbdr = - beta*ptot*alp

             dbeta = (beta_high-beta_low)/(nbeta-1)
             do j = sel, sel+nbeta-1
                
                beta_save = beta
                beta = beta_low + (j - sel)*dbeta
                 
                fac = -(dbdr/beta+alt)/aln
                fp_save = spec%fprim
                spec%fprim = fac*spec%fprim

                write (tag1, fmt='("Varying beta and all L_n values self-consistently")')
                write (tag2, fmt='("beta = ",e16.10," and tprim values scaled by ",e16.10)') beta, fac

                call write_namelists (j, tag1, tag2) 
                spec%fprim = fp_save
             end do
             beta = beta_save

             deallocate (fp_save)
          end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  1.2  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       case (2)
          call tell ('You have selected to vary beta non-self-consistently.')
          
          call beta_range_low (beta_low, 'lt', 0.)
          call beta_range_high (beta_high, 'le', beta_low)
          call num_runs (nbeta)

          call tell ('Preparing a non-self-consistent beta scan.')

          call run_number (sel, nbeta)

          write (tag1, fmt='(i3," runs prepared with beta_min = ",e16.10,&
               &" and beta_max = ",e16.10)') nbeta, beta_low, beta_high
          write (tag2, fmt='("Files are numbered from ",i3," to ",i3)') sel, sel+nbeta-1
          
          call tell (tag1, tag2)

          dbeta = (beta_high-beta_low)/(nbeta-1)
          do j = sel, sel+nbeta-1
             
             beta_save = beta
             fapar_save = fapar 
             faperp_save = faperp

             beta = beta_low + (j - sel)*dbeta
             if (beta == 0.) then 
                fapar = 0.
                faperp = 0.
             else
                if (fapar == 0. .and. faperp == 0.) then
                   fapar = 1.0 ;  faperp = 1.0
                end if
             end if
             
             write (tag1, fmt='("Varying beta, all else fixed")')
             write (tag2, fmt='("beta = ",e16.10)') beta

             call write_namelists (j, tag1, tag2)

             fapar = fapar_save 
             faperp = faperp_save
             beta = beta_save
          end do
          
       end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  2.0  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    case (2) ! beta_prime
       
       call tell ('You have chosen to vary beta_prime.')

115    continue
       call text
       call text ('Choose from the following:')
       call text ('(1) Vary beta_prime self-consistently')
       call text ('(2) Vary beta_prime with ALL other parameters held fixed (non-self-consistently).')
       call text
       call get_choice (2, sel)

       select case (sel)

       case default
          call tell ('Try again.  Choose an integer between 1 and 2, inclusively.')
          goto 115

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  2.1  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       case (1)
          call tell ('You have chosen to vary beta_prime self-consistently.')

116       continue
          call text
          call text ('Choose from the following:')
          call text ('(1) Hold gradient scale lengths fixed, vary beta')
          call text

          call get_choice (1, sel)

          select case (sel)

          case default
             call tell ('Try again.  Choose an integer between 1 and 1, inclusively.')
             goto 116

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  2.1.1  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          case (1)  
             call tell ('You have chosen to vary beta_prime while holding gradient scale lengths fixed.')

             call beta_prime_range_low (beta_low)
             call beta_prime_range_high (beta_high, beta_low)
             call num_runs (nbeta)

             call tell ('Preparing a self-consistent beta_prime scan.', &
                  'Beta will be varied to maintain consistency.')

             call run_number (sel, nbeta)

             write (tag1, fmt='(i3," runs prepared with beta_prime_min = ",e16.10,&
                  &" and beta_prime_max = ",e16.10)') nbeta, beta_low, beta_high
             write (tag2, fmt='("Files are numbered from ",i3," to ",i3)') sel, sel+nbeta-1
             
             call tell (tag1, tag2)

             ptot = 0.
             alp = 0.
             do is=1,nspec
                ptot = ptot + spec(is)%dens * spec(is)%temp
                alp = alp + spec(is)%dens * spec(is)%temp *(spec(is)%fprim + spec(is)%tprim)
             end do
            
             alp = alp/ptot
             dbdr = - beta*ptot*alp

             if (alp == 0.) then
                call tell ('Cannot proceed, because Lp = infinity', &
                     'No input files for scan written')
                return
             end if

             beta_save = beta
             beta_prime_save = dbdr

             fac = -1./(ptot*alp)

             dbeta = (beta_high-beta_low)/(nbeta-1)   ! actually, this is dbeta_prime
             do j = sel, sel+nbeta-1
                
                beta_prime_save = beta_prime_input
                beta_prime_input = beta_low + (j - sel)*dbeta
                
                beta_save = beta
                beta = beta_prime_input*fac

                fapar_save = fapar ; faperp_save = faperp
                if (beta == 0.) then
                   fapar = 0.      ; faperp = 0.
                else
                   if (fapar == 0. .and. faperp == 0.) then
                      fapar = 1.0 ;  faperp = 1.0
                   end if
                end if

                select case (bishop)
                case default
                   bishop_save = bishop
                   bishop = 6
                case (4)
                   ! nothing, continue to use bishop = 4
                end select

             
                write (tag1, fmt='("Varying beta_prime and beta self-consistently")')
                write (tag2, fmt='("beta_prime = ",e16.10," and beta = ",e16.10)') beta_prime_input, beta

                call write_namelists (j, tag1, tag2)

                fapar = fapar_save 
                faperp = faperp_save
                beta = beta_save
                beta_prime_input = beta_prime_save

                select case (bishop)
                case default
                   bishop = bishop_save
                case (4)
                   ! nothing, continue to use bishop = 4
                end select

             end do
          end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  2.2   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       case (2)

          call tell ('You have selected to vary beta_prime non-self-consistently.')          
          
          call beta_prime_range_low (beta_low)
          call beta_prime_range_high (beta_high, beta_low)
          call num_runs (nbeta)

          call tell ('Preparing a non-self-consistent beta_prime scan.')

          call run_number (sel, nbeta)
          
          write (tag1, fmt='(i3," runs prepared with beta_prime_min = ",e16.10,& 
               &" and beta_prime_max = ",e16.10)') nbeta, beta_low, beta_high 
          write (tag2, fmt='("Files are numbered from ",i3," to ",i3)') sel, sel+nbeta-1 
          
          call tell (tag1, tag2)

          dbeta = (beta_high-beta_low)/(nbeta-1) 
          do j = sel, sel+nbeta-1
             
             beta_prime_save = beta_prime_input
             beta_prime_input = beta_low + (j - sel)*dbeta
             
             select case (bishop)
             case default
                bishop_save = bishop
                bishop = 6
             case (4)
                ! nothing, continue to use bishop = 4
             end select
             
             write (tag1, fmt='("Varying beta_prime only (non-self-consistently)")')
             write (tag2, fmt='("beta_prime = ",e16.10)') beta_prime_input

             call write_namelists (j, tag1, tag2)
             
             beta_prime_input = beta_prime_save
             
             select case (bishop)
             case default
                bishop = bishop_save
             case (4)
                ! nothing, continue to use bishop = 4
             end select
             
          end do
          
       end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   1.3  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    case (3) ! collisionality
       call tell ('You have chosen to vary collisionality.')
       call text ('Not yet implemented (sorry).')
       call text

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   1.4  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    case (4) ! temperature gradient
       call tell ('You have chosen to vary temperature gradient.')
       call text ('Not yet implemented (sorry).')
       call text

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   1.5  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    case (5) ! density gradient
       call tell ('You have chosen to vary density gradient.')
       call text ('Not yet implemented (sorry).')
       call text

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   1.6  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    case (6) ! Z_effective
       call tell ('You have chosen to vary Z_effective.')
       call text ('Not yet implemented (sorry).')
       call text
       
    end select

10  format (a) 

  end subroutine interactive

  subroutine get_namelists

    use theta_grid, only: init_theta_grid, nbset, shat_real => shat
    logical :: list
!CMR
    logical:: debug=.false.
    call init_file_utils (list, name="template")
if (debug) write(6,*) 'get_namelists: called init_file_utils'
    ncut= 100000
    scan = .false.
    stdin = .true.
    pythonin = "."//trim(run_name)//".pythonin"
    in_file=input_unit_exist("ingen_knobs", exist)
    if (exist) read(unit=input_unit("ingen_knobs"), nml=ingen_knobs)

if (debug) write(6,*) 'get_namelists: if (scan), scan=',scan
    if (scan) then
       if (.not. stdin) then
          if (pythonin == "") then
             write (*,*) 'Need to specify pythonin in ingen_knobs.'
          end if 
       end if 
    end if 

if (debug) write(6,*) 'get_namelists: layouts'
    local_field_solve = .false.
!    layout = 'lexys'
    layout = 'lxyes'
    in_file=input_unit_exist("layouts_knobs", exist)
    if (exist) then
       read (unit=input_unit("layouts_knobs"), nml=layouts_knobs)
       layouts_write = .true.
    end if

    ! antenna: 
    w_antenna = (1., 0.0)
    t0 = -1.
    w_dot = 0.
    amplitude = 0.
    nk_stir = 1
    ant_off = .false.

    write_antenna = .false.
    in_file=input_unit_exist("driver",exist)
    if (exist .and. .not. ant_off) then
       read (unit=input_unit("driver"), nml=driver)

       allocate (kx_stir(nk_stir))
       allocate (ky_stir(nk_stir))
       allocate (kz_stir(nk_stir))
       allocate (a_ant(nk_stir))
       allocate (b_ant(nk_stir))
       allocate (trav(nk_stir))

       do i=1,nk_stir
          call get_indexed_namelist_unit (in_file, "stir", i)
          kx=1
          ky=1
          kz=1
          a=0.
          b=0.
          travel = .true.
          read (unit=in_file, nml=stir)
          close(unit=in_file)
          kx_stir(i) = kx 
          ky_stir(i) = ky
          kz_stir(i) = kz
          trav(i) = travel
          if (a == 0.0 .and. b == 0.0) then
             a_ant(i) = amplitude*cmplx(1.,1.)/2. 
             b_ant(i) = amplitude*cmplx(1.,1.)/2.
          else 
             a_ant(i) = a
             b_ant(i) = b
          end if
       end do
       driver_write = .true.       
    else
       no_driver = .true.
    end if

    t0_driver = t0

if (debug) write(6,*) 'get_namelists: collisions'
    ! collisions: 
    collision_model = 'default'
    tmpfac = cfac
    cfac = 0.
    vncoef = 0.6
    absom = 0.5
    ivnew = 0
    conserve_number = .true.
    conserve_momentum = .true.
    hypercoll = .false.
    heating = .false.
    in_file=input_unit_exist("collisions_knobs",exist)
    if (exist) then
       read (unit=input_unit("collisions_knobs"), nml=collisions_knobs)
       collisions_write = .true.
    end if

    ierr = error_unit()
    call get_option_value &
         (collision_model, coll_modelopts, collision_model_switch, &
         ierr, "collision_model in collisions_knobs")

    cfac_nu = cfac
    cfac = tmpfac

if (debug) write(6,*) 'get_namelists: init_g'
    ! init_g:
    tstart = 0.
    scale = 1.0
    ginit_option = "default"
    width0 = -3.5
    refac = 1.
    imfac = 0.
    den0 = 1.
    upar0 = 0.
    tpar0 = 0.
    tperp0 = 0.
    den1 = 0.
    upar1 = 0.
    tpar1 = 0.
    tperp1 = 0.
    den2 = 0.
    upar2 = 0.
    tpar2 = 0.
    tperp2 = 0.
    phiinit = 1.0
    dphiinit = 1.0
    zf_init = 1.0
    k0 = 1.0
    chop_side = .true.
    apar0 = 0.
    left = .true.
    even = .true.
    ikk(1) = 1
    ikk(2) = 2
    itt(1) = 1
    itt(2) = 2
    restart_file = trim(run_name)//".nc"
    in_file=input_unit_exist("init_g_knobs",exist)
    if (exist) then
       read (unit=input_unit("init_g_knobs"), nml=init_g_knobs)
       init_g_write = .true.
       initk0 = k0
    end if

    ierr = error_unit()
    call get_option_value &
         (ginit_option, ginitopts, ginitopt_switch, &
         ierr, "ginit_option in ginit_knobs")

    ! fields
    field_option="default"
    in_file = input_unit_exist("fields_knobs", exist)
    if (exist) then
       read (unit=input_unit("fields_knobs"), nml=fields_knobs)
       fields_write = .true.
    end if

    ierr = error_unit()
    call get_option_value &
         (field_option, fieldopts, fieldopt_switch, &
         ierr, "field_option in fields_knobs")
    
if (debug) write(6,*) 'get_namelists: gs2_reinit'
    ! gs2_reinit:
    delt_adj = 2.0
    delt_minimum = 1.e-5
    in_file = input_unit_exist("reinit_knobs",exist)
    if(exist) then
       read (unit=in_file, nml=reinit_knobs)
       reinit_write = .true.
    end if

if (debug) write(6,*) 'get_namelists: hyper'
    ! hyper:
    const_amp = .false.
    include_kpar = .false.
    isotropic_shear = .true.
    D_hyperres = -10.
    D_hypervisc = -10.
    D_hyper = -10.
    hyper_option = 'default'
    omega_osc = 0.4
    gridnorm = .true.
    in_file=input_unit_exist("hyper_knobs",exist)
    if (exist) then
       read (unit=input_unit("hyper_knobs"), nml=hyper_knobs)
       hyper_write = .true.
    endif

    ierr = error_unit()
    call get_option_value &
         (hyper_option, hyperopts, hyper_option_switch, &
         ierr, "hyper_option in hyper_knobs")

    select case (hyper_option_switch)

       case (hyper_option_none)
          if (D_hyperres > 0.)  D_hyperres = -10.
          if (D_hypervisc > 0.) D_hypervisc = -10.

       case (hyper_option_visc)
          if (D_hyperres > 0.) D_hyperres = -10.

       case (hyper_option_res)
          if (D_hypervisc > 0.) D_hypervisc = -10.

       case (hyper_option_both)
          
          if (D_hyper >= 0.) then
             D_hyperres  = D_hyper
             D_hypervisc = D_hyper
          end if

    end select

if (debug) write(6,*) 'get_namelists: kt_grids'
    ! kt_grids:
    norm_option = 'default'
    grid_option = 'default'
    in_file=input_unit_exist("kt_grids_knobs",exist)
    if (exist) then
       read (unit=input_unit("kt_grids_knobs"), nml=kt_grids_knobs)
       kt_write = .true.
    end if
    
    ierr = error_unit()
    call get_option_value &
         (grid_option, gridopts, gridopt_switch, &
         ierr, "grid_option in kt_grids_knobs")

    ierr = error_unit()
    call get_option_value &
         (norm_option, normopts, normopt_switch, &
         ierr, "norm_option in kt_grids_knobs")

    select case (gridopt_switch) 
    case (gridopt_single) 
       aky = 0.4
       theta0 = 0.0
       akx = 0.0
       in_file=input_unit_exist("kt_grids_single_parameters",exist)
       if (exist) then
          read (unit=input_unit("kt_grids_single_parameters"), &
               nml=kt_grids_single_parameters)
          kt_single_write = .true.
       end if
       naky = 1
       ntheta0 = 1
       
    case (gridopt_range)
       
       naky = 1
       ntheta0 = 1
       aky_min = 0.0
       aky_max = 0.0
       theta0_min = 0.0
       theta0_max = 0.0
       in_file=input_unit_exist("kt_grids_range_parameters",exist)
       if (exist) then
          read (unit=input_unit("kt_grids_range_parameters"), &
               nml=kt_grids_range_parameters)
          kt_range_write = .true.
       end if

    case (gridopt_specified) 
       naky = 1
       ntheta0 = 1
       nx = 0
       ny = 0
       in_file=input_unit_exist("kt_grids_specified_parameters",exist)
       if (exist) then
          read (unit=input_unit("kt_grids_specified_parameters"), &
               nml=kt_grids_specified_parameters)
          kt_specified_write = .true.
       end if

       allocate (aky_tmp(max(naky,ntheta0)))
       allocate (theta0_tmp(max(naky,ntheta0)))
       allocate (akx_tmp(max(naky,ntheta0)))

       do i = 1, max(naky,ntheta0)
          aky = 0.4
          theta0 = 0.0
          akx = 0.0
          call get_indexed_namelist_unit (unit, "kt_grids_specified_element", i)
          read (unit=unit, nml=kt_grids_specified_element)
          close (unit)
          aky_tmp(i) = aky
          theta0_tmp(i) = theta0
          akx_tmp(i) = akx
       end do
          
    case (gridopt_box)
       naky = 0
       ntheta0 = 0
       ly = 0.0
       x0 = 0.
       y0 = 2.0
       nx = 0
       ny = 0
       jtwist = 1
       rtwist = 0.0
       in_file=input_unit_exist("kt_grids_box_parameters",exist)
       if (exist) then
          read (unit=input_unit("kt_grids_box_parameters"), &
               nml=kt_grids_box_parameters)
          kt_box_write = .true.
       end if
       y0_internal = y0
       if (y0 < 0) then
          y0_internal = -1./y0_internal
       else
          y0_internal = y0_internal
       end if
       if (ly == 0.) ly = 2.0*pi*y0_internal
       if (naky == 0) naky = (ny-1)/3 + 1
       if (ntheta0 == 0) ntheta0 = 2*((nx-1)/3) + 1
       if (rtwist == 0.) rtwist = real(jtwist)

    case (gridopt_xbox)
       ntheta0 = 1
       lx = 1.0
       aky = 0.2
       nx = 0
       in_file=input_unit_exist("kt_grids_xbox_parameters",exist)
       if (exist) then
          read (unit=input_unit("kt_grids_xbox_parameters"), &
               nml=kt_grids_xbox_parameters)
          kt_xbox_write = .true.
       end if

    end select

if (debug) write(6,*) 'get_namelists: le_grids'
    ! le_grids:
    nesub = 8
    nesuper = 2
    ngauss = 5
    negrid = -10
    ecut = 6.0
    bouncefuzz = 1e-5
    test = .false.
    in_file=input_unit_exist("le_grids_knobs", exist)
    if (exist) then
       read (unit=input_unit("le_grids_knobs"), nml=le_grids_knobs)
       le_write = .true.
    end if
    test_le = test

    ! user can choose not to set negrid (preferred for old algorithm)
    if (negrid == -10) then
       negrid = nesub + nesuper
    else  ! If user chose negrid, assume nesuper makes sense and check nesub
       if (.not. advanced_egrid) then
          if (negrid - nesuper /= nesub) then
             ! Report problem to error file, and continue, using nesuper and negrid
             ! (Note that nesub is not used anywhere else.)
             nesub = negrid - nesuper
             ierr = error_unit()
             write (unit=ierr, fmt='("Forcing nesub = ",i5)') nesub
          endif
       endif
    endif

if (debug) write(6,*) 'get_namelists: nonlinear terms'
    ! nonlinear_terms:
    nonlinear_mode = 'default'
    flow_mode = 'default'
    cfl = 0.1
    C_par = 0.1
    C_perp = 0.1
    p_x = 6.0
    p_y = 6.0
    p_z = 6.0

    in_file=input_unit_exist("nonlinear_terms_knobs",exist)
    if(exist) then
       read (unit=in_file,nml=nonlinear_terms_knobs)
       nonlinear_write = .true.
    end if

    ierr = error_unit()
    call get_option_value &
         (nonlinear_mode, nonlinearopts, nonlinear_mode_switch, &
         ierr, "nonlinear_mode in nonlinear_terms_knobs")
    call get_option_value &
         (flow_mode, flowopts, flow_mode_switch, &
         ierr, "flow_mode in nonlinear_terms_knobs")

    ! run_parameters:
    fbpar = -1.0
    faperp = 0.0
    beta = 0.0
    zeff = 1.0
    tite = 1.0
    teti = -100.0
    rhostar = 0.1
    wstar_units = .false.
    eqzip = .false.
    secondary = .true.
    tertiary = .false.
    delt_option = 'default'
    margin = 0.05

    in_file=input_unit_exist("parameters",exist)
    if (exist) then
       read (unit=input_unit("parameters"), nml=parameters)
       parameters_write = .true.
    end if

    in_file = input_unit_exist("knobs",exist)
    if (exist) then
       read (unit=input_unit("knobs"), nml=knobs)
       knobs_write = .true.
    end if

    if (secondary .and. tertiary) secondary = .false.

    if (teti /= -100.0) tite = teti

! Allow faperp-style initialization for backwards compatibility.
! Only fbpar is used outside of this subroutine.
       if (fbpar == -1.) then
          fbpar = faperp
       end if


    ierr = error_unit()
    call get_option_value &
         (delt_option, deltopts, delt_option_switch, ierr, &
         "delt_option in knobs")

if (debug) write(6,*) 'get_namelists: species'
    ! species: 
    nspec = 2
    in_file = input_unit_exist("species_knobs", exist)
    if (exist) then
       read (unit=input_unit("species_knobs"), nml=species_knobs)
       species_write = .true.
if (debug) write(6,*) 'get_namelists: species_knobs exist, nspec=',nspec
    end if

if (debug) write(6,*) 'get_namelists: nspec=',nspec
    allocate (spec(nspec))
    do is = 1, nspec
       call get_indexed_namelist_unit (unit, "species_parameters", is)
       dens0 = 1.0
       uprim = 0.0
       uprim2 = 0.0
       vnewk = 0.0
       vnewk4 = 0.0
       type = "default"
       read (unit=unit, nml=species_parameters, iostat=iostat)
       close (unit=unit)

       spec(is)%z = z
       spec(is)%mass = mass
       spec(is)%dens = dens
       spec(is)%dens0 = dens0
       spec(is)%temp = temp
       spec(is)%tprim = tprim
       spec(is)%fprim = fprim
       spec(is)%uprim = uprim
       spec(is)%uprim2 = uprim2
       spec(is)%vnewk = vnewk
       spec(is)%vnewk4 = vnewk4

       ierr = error_unit()
       call get_option_value (type, typeopts, spec(is)%type, &
            ierr, "type in species_parameters_x")
    end do

    gb_to_cv = .false.
if (debug) write(6,*) 'get_namelists: theta_grid_knobs'
    equilibrium_option = 'default'
    in_file= input_unit_exist("theta_grid_knobs", exist)
    if (exist) then
       read (unit=input_unit("theta_grid_knobs"), nml=theta_grid_knobs)
       theta_write = .true.
    end if
       
    ierr = error_unit()
    call get_option_value &
         (equilibrium_option, eqopts, eqopt_switch, &
         ierr, "equilibrium_option in theta_grid_knobs")

    ! theta_grid: 
    rhoc = 0.5
    rmaj = 3.0
    r_geo = 3.0
    eps = 0.3
    epsl = 0.3
    qinp = 1.5
    shat = 0.75
    pk = 0.3
    kp = -1.
    shift = 0.0
    akappa = 1.0
    akappri = 0.0
    tri = 0.0
    tripri = 0.0
    ntheta = 24
    nperiod = 2

if (debug) write(6,*) 'get_namelists: case (eqopt_switch), eqopt_switch=',eqopt_switch
    select case (eqopt_switch)
    case (eqopt_eik)
       itor = 1
       iflux = 0
       irho = 2
       equal_arc = .true.
       bishop = 5
       dp_mult = 1.0
       delrho = 1e-3
       rmin = 1e-3
       rmax = 1.0
       ismooth = 0
       isym = 0
       writelots = .false.
       local_eq = .true.
       transp_eq = .false.
       idfit_eq=.false.
       dfit_eq=.false.
       gs2d_eq=.false.
       beta_prime_input = 0.
       s_hat_input = 0.
       
       in_file = input_unit_exist("theta_grid_parameters", exist)
       if (exist) then
          read (unit=input_unit("theta_grid_parameters"), nml=theta_grid_parameters)
          theta_parameters_write = .true.
       end if

       in_file= input_unit_exist("theta_grid_eik_knobs", exist)
       if (exist) then
if (debug) write(6,*) "get_namelists: reading theta_grid_eik_knobs"
          read (unit=input_unit("theta_grid_eik_knobs"), nml=theta_grid_eik_knobs)
if (debug) write(6,*) "get_namelists: done theta_grid_eik_knobs"
          theta_eik_write = .true.
       end if

     case (eqopt_salpha)
       in_file = input_unit_exist("theta_grid_parameters", exist)
       if (exist) then
          read (unit=input_unit("theta_grid_parameters"), nml=theta_grid_parameters)
          theta_parameters_write = .true.
       end if

       alpmhdfac = 0.0
       alpha1 = 0.0
       model_option = 'default'
       in_file = input_unit_exist("theta_grid_salpha_knobs", exist)
       if (exist) then
          read (unit=input_unit("theta_grid_salpha_knobs"), nml=theta_grid_salpha_knobs)
          theta_salpha_write = .true.
       end if

       ierr = error_unit()
       call get_option_value &
            (model_option, sa_modelopts, model_switch, &
            ierr, "model_option in theta_grid_salpha_knobs")

    case (eqopt_file)

       gridout_file = "grid.out"
       in_file= input_unit_exist("theta_grid_file_knobs", exist)
       if (exist) then
          read (unit=input_unit("theta_grid_file_knobs"), nml=theta_grid_file_knobs)
          theta_file_write = .true.
       end if

    end select


    if (kp > 0.) pk = 2.*kp

    npadd = 2
    alknob = 0.0
    epsknob = 1e-5
    bpknob = 1.e-8
    extrknob = 0.0
    tension = 1.0
    thetamax = 0.0
    deltaw = 0.0
    widthw = 1.0
    in_file = input_unit_exist("theta_grid_gridgen_knobs", exist)
    if (exist) then
       read (unit=input_unit("theta_grid_gridgen_knobs"), nml=theta_grid_gridgen_knobs)
       theta_gridgen_write = .true.
    end if

    ! dist_fn
    save_n = .true.
    save_u = .true.
    save_Tpar = .true.
    save_Tperp = .true.
    boundary_option = 'default'
    adiabatic_option = 'default'
    heating_option = 'default'
    poisfac = 0.0
    gridfac = 5e4
    apfac = 1.0
    driftknob = 1.0
    tpdriftknob = -9.9e9
    t0 = 100.0
    source0 = 1.0
    omega0 = 0.0
    gamma0 = 0.0
    thetas = 1.0
    aky_star = 0.0
    akx_star = 0.0
    phi_ext = 0.0
    a_ext = 0.0
    afilter = 0.0
    kfilter = 0.0
    g_exb = 0.0
    D_kill = -10.0
    noise = -1.
    mult_imp = .false.
    test = .false.
    def_parity = .false.
    even = .true.
    source_option = 'default'
    nperiod_guard = 0
    tmpfac = cfac
    cfac = 1.0

if (debug) write(6,*) 'get_namelists: dist_fn_knobs'
    in_file= input_unit_exist("dist_fn_knobs", exist)
    if (exist) then
       read (unit=input_unit("dist_fn_knobs"), nml=dist_fn_knobs)
       if (tpdriftknob == -9.9e9) tpdriftknob=driftknob
       dist_fn_write = .true.
    end if

    cfac_df = cfac
    cfac = tmpfac
    test_df = test

    in_file= input_unit_exist("source_knobs", exist)
    if (exist) then
       read (unit=input_unit("source_knobs"), nml=source_knobs)
       source_write = .true.
    end if

    allocate (fexp(nspec), bkdiff(nspec), bd_exp(nspec))

    do is = 1, nspec
       fexp(is) = (0.4, 0.)
       bkdiff(is) = 0.
       bd_exp(is) = 0
       call get_indexed_namelist_unit (unit, "dist_fn_species_knobs", is)
       call fill_species_knobs (unit, fexp(is), bkdiff(is), bd_exp(is))
       close (unit=unit)
    end do

if (debug) write(6,*) 'get_namelists: call init_theta_grid'
    call init_theta_grid 
if (debug) write(6,*) 'get_namelists: called init_theta_grid'
    shat = shat_real

    if(abs(shat) <=  1.e-5) boundary_option = 'periodic'

    ierr = error_unit()
    call get_option_value &
         (boundary_option, boundaryopts, boundary_option_switch, &
         ierr, "boundary_option in dist_fn_knobs")
    
    call get_option_value &
         (heating_option, heatingopts, heating_option_switch, &
         ierr, "heating_option in dist_fn_knobs")
    
    call get_option_value &
         (source_option, sourceopts, source_option_switch, &
         ierr, "source_option in source_knobs")
    call get_option_value &
         (adiabatic_option, adiabaticopts, adiabatic_option_switch, &
         ierr, "adiabatic_option in dist_fn_knobs")
    

    ! gs2_diagnostics
    print_line = .true.
    print_old_units = .false.
    print_flux_line = .false.
    write_line = .true.
    write_flux_line = .true.
    write_phi = .true.
    write_kpar = .false.
    write_apar = .true.
    write_aperp = .true.
    write_omega = .false.
    write_ascii = .true.
    write_hrate = .false.
    write_Epolar = .false.
    write_gs = .false.
    write_lorentzian = .false.
    write_tavg = .false.
    write_lamavg = .false.
    write_eavg = .false.
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
    write_avg_moments = .false.
    write_final_fields = .false.
    write_final_antot = .false.
    write_final_epar = .false.
    write_fcheck = .false.
    write_intcheck = .false.
    write_vortcheck = .false.
    write_fieldcheck = .false.
    nsave = -1
    nwrite = 100
    nmovie = 1000
    navg = 100
    nperiod_output = 1
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
    save_for_restart = .false.
    in_file = input_unit_exist("gs2_diagnostics_knobs", exist)
    if (exist) then
       read (unit=input_unit("gs2_diagnostics_knobs"), nml=gs2_diagnostics_knobs)
       diagnostics_write = .true.
    end if

    if (.not. save_for_restart) nsave = -1
    nperiod_output = max(nperiod, 1)
if (debug) write(6,*) 'get_namelists: returning'

  end subroutine get_namelists

  subroutine write_namelists (jr, tag1, tag2)

    integer, intent (in), optional :: jr
    character (*), intent (in), optional :: tag1, tag2

    character (100) :: line
    integer :: th, h, t, u
    integer :: i, j
    character (4) :: suffix
    character(20) :: datestamp, timestamp, zone
    
    call get_unused_unit (unit)

    if (present(jr)) then

       h = jr / 100
       t = (jr - h * 100) / 10
       u = (jr - h * 100 - t * 10)
       suffix = '_'//achar(48+h)//achar(48+t)//achar(48+u)
       open (unit=unit, file=trim(run_name)//suffix//".in")
    else
       open (unit=unit, file=trim(run_name)//".inp")
    endif

    write (unit, *)

    write (unit, fmt="('gs2')")
    datestamp(:) = ' '
    timestamp(:) = ' '
    zone(:) = ' '
    call date_and_time (datestamp, timestamp, zone)
    write (unit=unit, fmt="('Date: ',a,' Time: ',a,1x,a)") &
         trim(datestamp), trim(timestamp), trim(zone)

    if (present(tag1)) then
       write (unit, *) '*****************************************************'
       write (unit, *) trim(tag1)
       if (present(tag2)) write (unit, *) trim(tag2)
       write (unit, *) '*****************************************************'
    end if

    if (layouts_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "layouts_knobs"
       write (unit, fmt="(' layout = ',a)") '"'//trim(layout)//'"'
       write (unit, fmt="(' local_field_solve = ',L1)") local_field_solve
       write (unit, fmt="(' /')")
    end if

    if (collisions_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "collisions_knobs"
       select case (collision_model_switch)
       case (collision_model_lorentz)
          write (unit, fmt="(' collision_model = ',a)") '"lorentz"'
          write (unit, fmt="(' conserve_momentum = ',L1)") conserve_momentum
          if (hypercoll) write (unit, fmt="(' hypercoll = ',L1)") hypercoll
       case (collision_model_lorentz_test)
          write (unit, fmt="(' collision_model = ',a)") '"lorentz-test"'
          write (unit, fmt="(' conserve_momentum = ',L1)") conserve_momentum
       case (collision_model_krook)
          write (unit, fmt="(' collision_model = ',a)") '"krook"'
          write (unit, fmt="(' conserve_number = ',L1)") conserve_number
          write (unit, fmt="(' conserve_momentum = ',L1)") conserve_momentum
          write (unit, fmt="(' vncoef = ',f5.3)") vncoef
       case (collision_model_krook_test)
          write (unit, fmt="(' collision_model = ',a)") '"krook-test"'
          write (unit, fmt="(' conserve_number = ',L1)") conserve_number
          write (unit, fmt="(' conserve_momentum = ',L1)") conserve_momentum
       case (collision_model_none)
          write (unit, fmt="(' collision_model = ',a)") '"collisionless"'
       end select
       write (unit, fmt="(' cfac = ',f5.3)") cfac_nu
       write (unit, fmt="(' heating = ',L1)") heating
       write (unit, fmt="(' cfac = ',f5.3)") cfac
       write (unit, fmt="(' /')")
    end if

    if (init_g_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "init_g_knobs"
       select case (ginitopt_switch)

       case (ginitopt_default)
          write (unit, fmt="(' ginit_option = ',a)") '"default"'
          write (unit, fmt="(' phiinit = ',e16.10)") phiinit
          write (unit, fmt="(' width0 = ',e16.10)") width0
          write (unit, fmt="(' chop_side = ',L1)") chop_side
          if (chop_side) write (unit, fmt="(' left = ',L1)") left

       case (ginitopt_noise)
          write (unit, fmt="(' ginit_option = ',a)") '"noise"'
          write (unit, fmt="(' phiinit = ',e16.10)") phiinit
          write (unit, fmt="(' zf_init = ',e16.10)") zf_init
          write (unit, fmt="(' chop_side = ',L1)") chop_side
          if (chop_side) write (unit, fmt="(' left = ',L1)") left

       case (ginitopt_test1)
          write (unit, fmt="(' ginit_option = ',a)") '"test1"'

       case (ginitopt_xi)
          write (unit, fmt="(' ginit_option = ',a)") '"xi"'
          write (unit, fmt="(' width0 = ',e16.10)") width0

       case (ginitopt_xi2)
          write (unit, fmt="(' ginit_option = ',a)") '"xi2"'
          write (unit, fmt="(' width0 = ',e16.10)") width0

       case (ginitopt_zero)
          write (unit, fmt="(' ginit_option = ',a)") '"zero"'

       case (ginitopt_test3)
          write (unit, fmt="(' ginit_option = ',a)") '"test3"'

       case (ginitopt_convect)
          write (unit, fmt="(' ginit_option = ',a)") '"convect"'
          write (unit, fmt="(' k0 = ',e16.10)") initk0

       case (ginitopt_rh)
          write (unit, fmt="(' ginit_option = ',a)") '"rh"'

       case (ginitopt_restart_many)
          write (unit, fmt="(' ginit_option = ',a)") '"many"'
          write (unit, fmt="(' restart_file = ',a)") '"'//trim(restart_file)//'"'
          write (unit, fmt="(' scale = ',e16.10)") scale

       case (ginitopt_restart_small)
          write (unit, fmt="(' ginit_option = ',a)") '"small"'
          write (unit, fmt="(' phiinit = ',e16.10)") phiinit
          write (unit, fmt="(' zf_init = ',e16.10)") zf_init
          write (unit, fmt="(' chop_side = ',L1)") chop_side
          if (chop_side) write (unit, fmt="(' left = ',L1)") left
          write (unit, fmt="(' restart_file = ',a)") '"'//trim(restart_file)//'"'
          write (unit, fmt="(' scale = ',e16.10)") scale

       case (ginitopt_restart_file)
          write (unit, fmt="(' ginit_option = ',a)") '"file"'
          write (unit, fmt="(' restart_file = ',a)") '"'//trim(restart_file)//'"'
          write (unit, fmt="(' scale = ',e16.10)") scale

       case (ginitopt_continue)
          write (unit, fmt="(' ginit_option = ',a)") '"cont"'

       case (ginitopt_kz0)
          write (unit, fmt="(' ginit_option = ',a)") '"kz0"'
          write (unit, fmt="(' phiinit = ',e16.10)") phiinit
          write (unit, fmt="(' chop_side = ',L1)") chop_side
          if (chop_side) write (unit, fmt="(' left = ',L1)") left

       case (ginitopt_nl)
          write (unit, fmt="(' ginit_option = ',a)") '"nl"'
          write (unit, fmt="(' phiinit = ',e16.10)") phiinit
          write (unit, fmt="(' ikk(1) = ',i3,' itt(1) = ',i3)") ikk(1),itt(1)
          write (unit, fmt="(' ikk(2) = ',i3,' itt(2) = ',i3)") ikk(2), itt(2)
          write (unit, fmt="(' chop_side = ',L1)") chop_side
          if (chop_side) write (unit, fmt="(' left = ',L1)") left

       case (ginitopt_nl2)
          write (unit, fmt="(' ginit_option = ',a)") '"nl2"'
          write (unit, fmt="(' phiinit = ',e16.10)") phiinit
          write (unit, fmt="(' ikk(1) = ',i3,' itt(1) = ',i3)") ikk(1),itt(1)
          write (unit, fmt="(' ikk(2) = ',i3,' itt(2) = ',i3)") ikk(2), itt(2)
          write (unit, fmt="(' chop_side = ',L1)") chop_side
          if (chop_side) write (unit, fmt="(' left = ',L1)") left

       case (ginitopt_nl3)
          write (unit, fmt="(' ginit_option = ',a)") '"nl3"'
          write (unit, fmt="(' phiinit = ',e16.10)") phiinit
          write (unit, fmt="(' width0 = ',e16.10)") width0
          write (unit, fmt="(' refac = ',e16.10)") refac
          write (unit, fmt="(' imfac = ',e16.10)") imfac
          write (unit, fmt="(' ikk(1) = ',i3,' itt(1) = ',i3)") ikk(1),itt(1)
          write (unit, fmt="(' ikk(2) = ',i3,' itt(2) = ',i3)") ikk(2), itt(2)
          write (unit, fmt="(' chop_side = ',L1)") chop_side
          if (chop_side) write (unit, fmt="(' left = ',L1)") left
          write (unit, fmt="(' den0 = ',e16.10)") den0
          write (unit, fmt="(' den1 = ',e16.10)") den1
          write (unit, fmt="(' den2 = ',e16.10)") den2
          write (unit, fmt="(' upar0 = ',e16.10)") upar0
          write (unit, fmt="(' upar1 = ',e16.10)") upar1
          write (unit, fmt="(' upar2 = ',e16.10)") upar2
          write (unit, fmt="(' tpar0 = ',e16.10)") tpar0
          write (unit, fmt="(' tpar1 = ',e16.10)") tpar1
          write (unit, fmt="(' tperp0 = ',e16.10)") tperp0
          write (unit, fmt="(' tperp1 = ',e16.10)") tperp1
          write (unit, fmt="(' tperp2 = ',e16.10)") tperp2

       case (ginitopt_nl4)
          write (unit, fmt="(' ginit_option = ',a)") '"nl4"'
          write (unit, fmt="(' phiinit = ',e16.10)") phiinit
          write (unit, fmt="(' restart_file = ',a)") '"'//trim(restart_file)//'"'
          write (unit, fmt="(' scale = ',e16.10)") scale
          write (unit, fmt="(' ikk(1) = ',i3,' itt(1) = ',i3)") ikk(1),itt(1)
          write (unit, fmt="(' ikk(2) = ',i3,' itt(2) = ',i3)") ikk(2), itt(2)
          write (unit, fmt="(' chop_side = ',L1)") chop_side
          if (chop_side) write (unit, fmt="(' left = ',L1)") left

       case (ginitopt_nl5)
          write (unit, fmt="(' ginit_option = ',a)") '"nl5"'
          write (unit, fmt="(' phiinit = ',e16.10)") phiinit
          write (unit, fmt="(' restart_file = ',a)") '"'//trim(restart_file)//'"'
          write (unit, fmt="(' scale = ',e16.10)") scale
          write (unit, fmt="(' chop_side = ',L1)") chop_side
          if (chop_side) write (unit, fmt="(' left = ',L1)") left

       case (ginitopt_nl6)
          write (unit, fmt="(' ginit_option = ',a)") '"nl6"'
          write (unit, fmt="(' phiinit = ',e16.10)") phiinit
          write (unit, fmt="(' restart_file = ',a)") '"'//trim(restart_file)//'"'
          write (unit, fmt="(' scale = ',e16.10)") scale

       case (ginitopt_alf)
          write (unit, fmt="(' ginit_option = ',a)") '"alf"'
          write (unit, fmt="(' phiinit = ',e16.10)") phiinit

       case (ginitopt_gs)
          write (unit, fmt="(' ginit_option = ',a)") '"gs"'
          write (unit, fmt="(' phiinit = ',e16.10)") phiinit
          write (unit, fmt="(' refac = ',e16.10)") refac
          write (unit, fmt="(' imfac = ',e16.10)") imfac
          write (unit, fmt="(' den1 = ',e16.10)") den1
          write (unit, fmt="(' upar1 = ',e16.10)") upar1
          write (unit, fmt="(' tpar1 = ',e16.10)") tpar1
          write (unit, fmt="(' tperp1 = ',e16.10)") tperp1


       case (ginitopt_kpar)
          write (unit, fmt="(' ginit_option = ',a)") '"kpar"'
          write (unit, fmt="(' phiinit = ',e16.10)") phiinit
          write (unit, fmt="(' width0 = ',e16.10)") width0
          write (unit, fmt="(' refac = ',e16.10)") refac
          write (unit, fmt="(' imfac = ',e16.10)") imfac
          write (unit, fmt="(' den0 = ',e16.10)") den0
          write (unit, fmt="(' den1 = ',e16.10)") den1
          write (unit, fmt="(' den2 = ',e16.10)") den2
          write (unit, fmt="(' upar0 = ',e16.10)") upar0
          write (unit, fmt="(' upar1 = ',e16.10)") upar1
          write (unit, fmt="(' upar2 = ',e16.10)") upar2
          write (unit, fmt="(' tpar0 = ',e16.10)") tpar0
          write (unit, fmt="(' tpar1 = ',e16.10)") tpar1
          write (unit, fmt="(' tperp0 = ',e16.10)") tperp0
          write (unit, fmt="(' tperp1 = ',e16.10)") tperp1
          write (unit, fmt="(' tperp2 = ',e16.10)") tperp2

       end select
       write (unit, fmt="(' /')")
       
    end if

    if (dist_fn_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "dist_fn_knobs"

       select case (boundary_option_switch)

       case (boundary_option_zero)
          write (unit, fmt="(' boundary_option = ',a)") '"default"'

       case (boundary_option_self_periodic)
          write (unit, fmt="(' boundary_option = ',a)") '"periodic"'

       case (boundary_option_linked)
          write (unit, fmt="(' boundary_option = ',a)") '"linked"'

       case (boundary_option_alternate_zero)
          write (unit, fmt="(' boundary_option = ',a)") '"alternate-zero"'

       end select

       write (unit, fmt="(' gridfac = ',e16.10)") gridfac

       if (.not. has_electrons) then
          select case (adiabatic_option_switch)
             
          case (adiabatic_option_default)
             write (unit, *)
             write (unit, fmt="(' adiabatic_option = ',a)") &
                  & '"no-field-line-average-term"'
             
          case (adiabatic_option_fieldlineavg)
             write (unit, fmt="(' adiabatic_option = ',a)") '"field-line-average-term"'
             
          case (adiabatic_option_yavg)
             write (unit, fmt="(' adiabatic_option = ',a)") '"iphi00=3"'
             
          case (adiabatic_option_noJ)
             write (unit, fmt="(' adiabatic_option = ',a)") '"dimits"'
             
          end select
       end if

       if (apfac /= 1.) write (unit, fmt="(' apfac = ',e16.10)") apfac
       if (driftknob /= 1.) write (unit, fmt="(' driftknob = ',e16.10)") driftknob
       if (tpdriftknob /= 1.) write (unit, fmt="(' tpdriftknob = ',e16.10)") tpdriftknob
       if (poisfac /= 0.) write (unit, fmt="(' poisfac = ',e16.10)") poisfac
       if (kfilter /= 0.) write (unit, fmt="(' kfilter = ',e16.10)") kfilter
       if (afilter /= 0.) write (unit, fmt="(' afilter = ',e16.10)") afilter
       if (nperiod_guard /= 0) &
            write (unit, fmt="(' nperiod_guard = ',i2)") nperiod_guard
       if (mult_imp) write (unit, fmt="(' mult_imp = ',L1)") mult_imp
       if (test) write (unit, fmt="(' test = ',L1)") test
       if (def_parity) then
          write (unit, fmt="(' def_parity = ',L1)") def_parity
          if (even) write (unit, fmt="(' even = ',L1)") even
       end if
       if (D_kill > 0.) then
          write (unit, fmt="(' D_kill = ',e16.10)") D_kill
          write (unit, fmt="(' save_n = ',L1)") save_n
          write (unit, fmt="(' save_u = ',L1)") save_u
          write (unit, fmt="(' save_Tpar = ',L1)") save_Tpar
          write (unit, fmt="(' save_Tperp = ',L1)") save_Tperp
       end if
       if (noise > 0.) write (unit, fmt="(' noise = ',e16.10)") noise
       if (cfac_df < 1.) write (unit, fmt="(' cfac = ',e16.10)") cfac_df
       write (unit, fmt="(' /')")
    end if

    if (source_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "source_knobs"
       select case (source_option_switch)

       case (source_option_full)
          write (unit, fmt="(' source_option = ',a)") '"full"'

       case(source_option_phiext_full)
          write (unit, fmt="(' course_option = ',a)") '"phiext_full"'
          write (unit, fmt="(' source0 = ',e16.10)") source0
          write (unit, fmt="(' omega0 = ',e16.10)") omega0
          write (unit, fmt="(' gamma0 = ',e16.10)") gamma0
          write (unit, fmt="(' t0 = ',e16.10)") t0
          write (unit, fmt="(' phi_ext = ',e16.10)") phi_ext
       
       case(source_option_test2_full)
          write (unit, fmt="(' course_option = ',a)") '"test2_full"'
          write (unit, fmt="(' source0 = ',e16.10)") source0
          write (unit, fmt="(' omega0 = ',e16.10)") omega0
          write (unit, fmt="(' gamma0 = ',e16.10)") gamma0
          write (unit, fmt="(' t0 = ',e16.10)") t0
          write (unit, fmt="(' thetas = ',e16.10)") thetas

       case(source_option_convect_full)
          write (unit, fmt="(' course_option = ',a)") '"convect_full"'
          write (unit, fmt="(' source0 = ',e16.10)") source0
          write (unit, fmt="(' omega0 = ',e16.10)") omega0
          write (unit, fmt="(' gamma0 = ',e16.10)") gamma0
          write (unit, fmt="(' t0 = ',e16.10)") t0
          write (unit, fmt="(' k0 = ',e16.10)") k0

       case (source_option_zero)
          write (unit, fmt="(' course_option = ',a)") '"zero"'

       case (source_option_cosine)
          write (unit, fmt="(' course_option = ',a)") '"cosine"'
          write (unit, fmt="(' source0 = ',e16.10)") source0
          write (unit, fmt="(' omega0 = ',e16.10)") omega0
          write (unit, fmt="(' gamma0 = ',e16.10)") gamma0
          write (unit, fmt="(' t0 = ',e16.10)") t0

       case (source_option_sine)
          write (unit, fmt="(' course_option = ',a)") '"sine"'
          write (unit, fmt="(' source0 = ',e16.10)") source0
          write (unit, fmt="(' omega0 = ',e16.10)") omega0
          write (unit, fmt="(' gamma0 = ',e16.10)") gamma0
          write (unit, fmt="(' t0 = ',e16.10)") t0

       case (source_option_test1)
          write (unit, fmt="(' course_option = ',a)") '"test1"'
          write (unit, fmt="(' source0 = ',e16.10)") source0
          write (unit, fmt="(' omega0 = ',e16.10)") omega0
          write (unit, fmt="(' gamma0 = ',e16.10)") gamma0
          write (unit, fmt="(' t0 = ',e16.10)") t0

       end select
       write (unit, fmt="(' /')")
    end if

    if (fields_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "fields_knobs"
       select case (fieldopt_switch)
       case (fieldopt_implicit)
          write (unit, fmt="(' field_option = ',a)") '"implicit"'
       case (fieldopt_explicit)
          write (unit, fmt="(' field_option = ',a)") '"explicit"'
       case (fieldopt_test)
          write (unit, fmt="(' field_option = ',a)") '"test"'
       end select
       write (unit, fmt="(' /')")
    end if

    if (diagnostics_write) then
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
       if (write_intcheck) write (unit, fmt="(' write_intcheck = ',L1)") write_intcheck
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
    end if

    if (nonlinear_mode_switch == nonlinear_mode_on) then
       if (reinit_write) then
          write (unit, *)
          write (unit, fmt="(' &',a)") "reinit_knobs"
          write (unit, fmt="(' delt_adj = ',e16.10)") delt_adj
          write (unit, fmt="(' delt_minimum = ',e16.10)") delt_minimum
          write (unit, fmt="(' /')")       
       end if
    end if


    if (hyper_option_switch /= hyper_option_none) then
       if (hyper_write) then
          
          write (unit, *)
          write (unit, fmt="(' &',a)") "hyper_knobs"

          select case (hyper_option_switch)
             
          case (hyper_option_visc) 
             write (unit, fmt="(' hyper_option = ',a)") '"visc_only"'
             write (unit, fmt="(' D_hypervisc = ',e16.10)") D_hypervisc
             
          case (hyper_option_res) 
             write (unit, fmt="(' hyper_option = ',a)") '"res_only"'
             write (unit, fmt="(' D_hyperres = ',e16.10)") D_hyperres
             
          case (hyper_option_both) 
             write (unit, fmt="(' hyper_option = ',a)") '"both"'
             if (D_hyperres == D_hypervisc) then
                write (unit, fmt="(' D_hyper = ',e16.10)") D_hyper
             else
                write (unit, fmt="(' D_hypervisc = ',e16.10)") D_hypervisc
                write (unit, fmt="(' D_hyperres = ',e16.10)") D_hyperres
             end if
          end select

!          write (unit, fmt="(' include_kpar = ',L1)") include_kpar

          write (unit, fmt="(' const_amp = ',L1)") const_amp
          write (unit, fmt="(' isotropic_shear = ',L1)") isotropic_shear
          if (.not. isotropic_shear) &
               write (unit, fmt="(' omega_osc = ',e16.10)") omega_osc

          write (unit, fmt="(' gridnorm = ',L1)") gridnorm
          write (unit, fmt="(' /')")
       end if
    end if

    if (kt_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "kt_grids_knobs"

       select case (gridopt_switch)
          
       case (gridopt_single)
          write (unit, fmt="(' grid_option = ',a)") '"single"'

       case (gridopt_range)
          write (unit, fmt="(' grid_option = ',a)") '"range"'

       case (gridopt_specified)
          write (unit, fmt="(' grid_option = ',a)") '"specified"'

       case (gridopt_box)
          write (unit, fmt="(' grid_option = ',a)") '"box"'

       case (gridopt_xbox)
          write (unit, fmt="(' grid_option = ',a)") '"xbox"'

       end select

       select case (normopt_switch)
          
       case (normopt_mtk)
          write (unit, fmt="(' norm_option = ',a)") '"with_root_2"'

       case (normopt_bd)
          write (unit, fmt="(' norm_option = ',a)") '"no_root_2"'

       end select

       write (unit, fmt="(' /')")
    end if
    
    if (le_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "le_grids_knobs"
       write (unit, fmt="(' advanced_egrid = ',L1)") advanced_egrid
       if (advanced_egrid) then
          write (unit, fmt="(' negrid = ',i4)") negrid
       else
          write (unit, fmt="(' nesub = ',i4)") nesub
          write (unit, fmt="(' nesuper = ',i4)") nesuper
       end if
       write (unit, fmt="(' ngauss = ',i4)") ngauss
       write (unit, fmt="(' ecut = ',e16.10)") ecut
       write (unit, fmt="(' /')")
    end if


    if (kt_single_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "kt_grids_single_parameters"
       write (unit, fmt="(' aky = ',e16.10)") aky
       write (unit, fmt="(' theta0 = ',e16.10)") theta0
       write (unit, fmt="(' /')")
    end if

    if (kt_range_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "kt_grids_range_parameters"
       write (unit, fmt="(' naky = ',i3)") naky
       write (unit, fmt="(' aky_min = ',e16.10)") aky_min
       write (unit, fmt="(' aky_max = ',e16.10)") aky_max
       write (unit, fmt="(' ntheta0 = ',i3)") ntheta0
       write (unit, fmt="(' theta0_min = ',e16.10)") theta0_min
       write (unit, fmt="(' theta0_max = ',e16.10)") theta0_max
       write (unit, fmt="(' /')")
    end if

    if (kt_specified_write) then
       write(unit, kt_grids_specified_parameters)
       do i=1,max(naky,ntheta0)
          write (unit, *)
          write (line, *) i
          write (unit, fmt="(' &',a)") &
               & trim("kt_grids_specified_element_"//trim(adjustl(line)))
          write (unit, fmt="(' aky = ',e13.6,' theta0 = ',e13.6,'  /')") aky_tmp(i), theta0_tmp(i)
          write (unit, fmt="(' /')")
       end do
    end if

    if (kt_box_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "kt_grids_box_parameters"
       write (unit, fmt="(' nx = ',i4)") nx
       write (unit, fmt="(' ny = ',i4)") ny
       write (unit, fmt="(' Ly = ',e16.10)") ly
       if (jtwist /= 1) then
          write (unit, fmt="(' rtwist = ',e16.10)") rtwist
       else
          write (unit, fmt="(' jtwist = ',i4)") jtwist
       end if
       write (unit, fmt="(' /')")
    end if


    if (kt_xbox_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "kt_grids_xbox_parameters"
       write (unit, fmt="(' nx = ',i4)") nx
       write (unit, fmt="(' ny = ',i4)") ny
       write (unit, fmt="(' Lx = ',e16.10)") lx
       write (unit, fmt="(' /')")
    end if

    if (nonlinear_write .and. nonlinear_mode_switch == nonlinear_mode_on) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "nonlinear_terms_knobs"
       write (unit, fmt="(' nonlinear_mode = ',a)") '"on"'
       write (unit, fmt="(' cfl = ',e16.10)") cfl
       if (zip) write (unit, fmt="(' zip = ',L1)") zip
       write (unit, fmt="(' /')")
    end if

    if (parameters_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "parameters"
       write (unit, fmt="(' beta = ',e16.10)") beta       ! if zero, fapar, fbpar should be zero
       if (collision_model_switch /= collision_model_none) &
            write (unit, fmt="(' zeff = ',e16.10)") zeff
       if (.not. has_electrons)  write (unit, fmt="(' tite = ',e16.10)") tite
       if (zip) write (unit, fmt="(' zip = ',L1)") zip
       write (unit, fmt="(' /')")
    end if

    if (knobs_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "knobs"
       write (unit, fmt="(' fphi   = ',f6.3)") fphi
       write (unit, fmt="(' fapar  = ',f6.3)") fapar
       write (unit, fmt="(' fbpar = ',f6.3)") fbpar
       write (unit, fmt="(' delt = ',e16.10)") delt
       write (unit, fmt="(' nstep = ',i8)") nstep
       write (unit, fmt="(' wstar_units = ',L1)") wstar_units
       if (eqzip) then
          write (unit, fmt="(' eqzip = ',L1)") eqzip
          write (unit, fmt="(' secondary = ',L1)") secondary
          write (unit, fmt="(' tertiary = ',L1)") tertiary
       end if
       write (unit, fmt="(' margin = ',e16.10)") margin
       select case (delt_option_switch)
       case (delt_option_auto)
          write (unit, fmt="(' delt_option = ',a)") '"check_restart"'
       case (delt_option_hand)
          ! nothing
       end select
       write (unit, fmt="(' /')")
    end if

    if (species_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "species_knobs"
       write (unit, fmt="(' nspec = ',i2)") nspec
       write (unit, fmt="(' /')")

       do i=1,nspec
          write (unit, *)
          write (line, *) i
          write (unit, fmt="(' &',a)") &
               & trim("species_parameters_"//trim(adjustl(line)))
          write (unit, fmt="(' z = ',e13.6)") spec(i)%z
          write (unit, fmt="(' mass = ',e13.6)") spec(i)%mass
          write (unit, fmt="(' dens = ',e13.6)") spec(i)%dens
          write (unit, fmt="(' temp = ',e13.6)") spec(i)%temp
          write (unit, fmt="(' tprim = ',e13.6)") spec(i)%tprim
          write (unit, fmt="(' fprim = ',e13.6)") spec(i)%fprim
          write (unit, fmt="(' uprim = ',e13.6)") spec(i)%uprim
          if (spec(i)%uprim2 /= 0.) write (unit, fmt="(' uprim2 = ',e13.6)") spec(i)%uprim2
          write (unit, fmt="(' vnewk = ',e13.6)") spec(i)%vnewk
          if (spec(i)%type == ion_species) &
               write (unit, fmt="(a)") ' type = "ion" /'
          if (spec(i)%type == electron_species) &
               write (unit, fmt="(a)") ' type = "electron"  /'
          if (spec(i)%type == slowing_down_species) &
               write (unit, fmt="(a)") ' type = "fast"  /'
          write (unit, fmt="(' dens0 = ',e13.6)") spec(i)%dens0
       end do
    end if
    do i=1,nspec
       write (unit, *)
       write (line, *) i
       write (unit, fmt="(' &',a)") &
            & trim("dist_fn_species_knobs_"//trim(adjustl(line)))
       write (unit, fmt="(' fexpr = ',e13.6)") real(fexp(i))
       write (unit, fmt="(' bakdif = ',e13.6)") bkdiff(i)
       write (unit, fmt="(' bd_exp = ',i6,'  /')") bd_exp(i)
    end do

    if (theta_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "theta_grid_knobs"

       select case (eqopt_switch)

       case (eqopt_eik)
          write (unit, fmt="(a)") ' equilibrium_option = "eik"'
          
       case (eqopt_salpha)
          write (unit, fmt="(a)") ' equilibrium_option = "s-alpha"'

       case (eqopt_file)
          write (unit, fmt="(a)") ' equilibrium_option = "file"'

       end select
       write (unit, fmt="(' gb_to_cv = ',L1)") gb_to_cv
       write (unit, fmt="(' /')")
    end if

    if (theta_parameters_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "theta_grid_parameters"
       write (unit, fmt="(' ntheta =  ',i4)") ntheta
       write (unit, fmt="(' nperiod = ',i4)") nperiod
       write (unit, fmt="(' rhoc =    ',f7.4)") rhoc
       write (unit, fmt="(' Rmaj =    ',f7.4)") rmaj
       write (unit, fmt="(' R_geo =   ',f7.4)") r_geo
       write (unit, fmt="(' eps =     ',f7.4)") eps
       write (unit, fmt="(' epsl =    ',f7.4)") epsl
       write (unit, fmt="(' qinp =    ',f7.4)") qinp
       write (unit, fmt="(' shat =    ',f7.4)") shat
       write (unit, fmt="(' alpmhd =  ',f7.4)") alpmhd
       write (unit, fmt="(' pk =      ',f7.4)") pk
       write (unit, fmt="(' kp =      ',f7.4)") kp
       write (unit, fmt="(' shift =   ',f7.4)") shift
       write (unit, fmt="(' akappa =  ',f7.4)") akappa
       write (unit, fmt="(' akappri = ',f7.4)") akappri
       write (unit, fmt="(' tri =     ',f7.4)") tri
       write (unit, fmt="(' tripri =  ',f7.4)") tripri
       write (unit, fmt="(' /')")
    end if


    if (theta_gridgen_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "theta_grid_gridgen_knobs"
       write (unit, fmt="(' npadd =    ',i4)") npadd
       write (unit, fmt="(' alknob =   ',e16.10)") alknob
       write (unit, fmt="(' epsknob =  ',e16.10)") epsknob
       write (unit, fmt="(' bpknob =   ',e16.10)") bpknob
       write (unit, fmt="(' extrknob = ',e16.10)") extrknob
       write (unit, fmt="(' tension =  ',e16.10)") tension
       write (unit, fmt="(' thetamax = ',e16.10)") thetamax
       write (unit, fmt="(' deltaw =   ',e16.10)") deltaw
       write (unit, fmt="(' widthw =   ',e16.10)") widthw
       write (unit, fmt="(' /')")
    end if

    if (theta_salpha_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "theta_grid_salpha_knobs"
       write (unit, fmt="(' alpmhdfac = ',e16.10)") alpmhdfac
       write (unit, fmt="(' alpha1 =    ',e16.10)") alpha1

       select case (model_switch)

       case (model_salpha)
          write (unit, fmt="(a)") ' model_option = "s-alpha"'
          
       case (model_alpha1)
          write (unit, fmt="(a)") ' model_option = "alpha1"'

       case (model_eps)
          write (unit, fmt="(a)") ' model_option = "rogers"'
          
       case (model_b2)
          write (unit, fmt="(a)") ' model_option = "b2"'
          
       case (model_normal_only)
          write (unit, fmt="(a)") ' model_option = "normal_only"'

       case (model_ccurv)
          write (unit, fmt="(a)") ' model_option = "const-curv"'

       case (model_nocurve)
          write (unit, fmt="(a)") ' model_option = "no-curvature"'
          
       end select
       write (unit, fmt="(' /')")
    end if

    if (theta_eik_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "theta_grid_eik_knobs"
       write (unit, fmt="(' itor =  ',i2)") itor
       write (unit, fmt="(' iflux =  ',i2)") iflux
       write (unit, fmt="(' irho =  ',i2)") irho
       write (unit, fmt="(' ppl_eq =   ',L1)") ppl_eq
       write (unit, fmt="(' efit_eq =  ',L1)") efit_eq
       write (unit, fmt="(' gen_eq =   ',L1)") gen_eq
       write (unit, fmt="(' vmom_eq =  ',L1)") vmom_eq
       write (unit, fmt="(' dfit_eq =  ',L1)") dfit_eq
!       write (unit, fmt="(' idfit_eq = ',L1)") idfit_eq
       write (unit, fmt="(' local_eq =  ',L1)") local_eq
       write (unit, fmt="(' transp_eq =  ',L1)") transp_eq
       write (unit, fmt="(' gs2d_eq =  ',L1)") gs2d_eq
       write (unit, fmt="(' equal_arc =  ',L1)") equal_arc
       write (unit, fmt="(' bishop =  ',i2)") bishop
       write (unit, fmt="(' s_hat_input =  ',e13.6)") s_hat_input
       write (unit, fmt="(' alpha_input =  ',e13.6)") alpha_input
       write (unit, fmt="(' invLp_input =  ',e13.6)") invLp_input
       write (unit, fmt="(' beta_prime_input =  ',e13.6)") beta_prime_input
       write (unit, fmt="(' dp_mult =  ',e13.6)") dp_mult
       write (unit, fmt="(' delrho =  ',e13.6)") delrho
       write (unit, fmt="(' rmin =  ',e13.6)") rmin
       write (unit, fmt="(' rmax =  ',e13.6)") rmax
       write (unit, fmt="(' ismooth =  ',i1)") ismooth
       write (unit, fmt="(' isym =  ',i1)") isym
       write (unit, fmt="(' writelots =  ',L1)") writelots
       write (unit, fmt="(' ak0 =  ',e13.5)") ak0
!       write (unit, fmt="(' k1 =  ',e13.5)") k1
!       write (unit, fmt="(' k2 =  ',e13.5)") k2
       write (unit, fmt="(' eqfile = ',a)") '"'//trim(eqfile)//'"'
       write (unit, fmt="(' /')")
    end if

    if (theta_file_write) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "theta_grid_file_knobs"
       write (unit, fmt="(' gridout_file = ',a)") '"'//trim(gridout_file)//'"'
       write (unit, fmt="(' /')")
    end if

    if (driver_write) then 
       write (unit, *)
       write (unit, fmt="(' &',a)") "driver"
       write (unit, fmt="(' ant_off = ',L1)") ant_off
       write (unit, fmt="(' write_antenna = ',L1)") write_antenna
       write (unit, fmt="(' amplitude = ',e16.10)") amplitude
       write (unit, fmt="(' w_antenna = (',e16.10,', ',e16.10,')')") w_antenna
       write (unit, fmt="(' w_dot = ',e16.10)") w_dot
       write (unit, fmt="(' t0 = ',e16.10)") t0_driver
       write (unit, fmt="(' nk_stir = ',i3)") nk_stir
       write (unit, fmt="(' /')")

        do i=1,nk_stir
           write (unit, *)
           write (line, *) i
           write(unit, fmt="(' &',a)") trim("stir_"//trim(adjustl(line)))
           write(unit, fmt="(' kx = ',i2,' ky = ',i2,' kz = ',i2)") &
                kx_stir(i), ky_stir(i), kz_stir(i)
           write(unit, fmt="(' travel = ',L1)") trav(i)
           write(unit, fmt="(' a = (',e19.13,',',e19.13,')')") a_ant(i)
           write(unit, fmt="(' b = (',e19.13,',',e19.13,') /')") b_ant(i)
        end do
     end if

     write(unit, fmt=*)
     close (unit)

   end subroutine write_namelists

   subroutine fill_species_knobs (unit, fexp_out, bakdif_out, bd_exp_out)
     implicit none
     integer, intent (in) :: unit
     complex, intent (in out) :: fexp_out
     real, intent (in out) :: bakdif_out
     integer, intent (in out) :: bd_exp_out
     integer :: bd_exp, iostat
     real :: fexpr, fexpi, bakdif
     namelist /dist_fn_species_knobs/ fexpr, fexpi, bakdif, bd_exp

     fexpr = real(fexp_out)
     fexpi = aimag(fexp_out)
     bakdif = bakdif_out
     bd_exp = bd_exp_out
     read (unit=unit, nml=dist_fn_species_knobs, iostat=iostat)
     fexp_out = cmplx(fexpr,fexpi)
     bd_exp_out = bd_exp
     bakdif_out = bakdif
   end subroutine fill_species_knobs

   subroutine check(veq, geq, eeq, peq, leq, deq, ideq, report_unit)
     logical, intent(in) :: veq, geq, eeq, peq, leq, deq, ideq
     integer, intent (in) :: report_unit

     if(veq .and. geq) then
        write (report_unit, *) 
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing vmom_eq = .true. AND gen_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *) 
     endif                      

     if(veq .and. deq) then     
        write (report_unit, *) 
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing vmom_eq = .true. AND dfit_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *) 
     endif                      

     if(veq .and. eeq) then     
        write (report_unit, *) 
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing vmom_eq = .true. AND efit_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *) 
     endif                      

     if(veq .and. leq) then     
        write (report_unit, *) 
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing vmom_eq = .true. AND local_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *) 
     endif                      

     if(veq .and. peq) then     
        write (report_unit, *) 
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing vmom_eq = .true. AND ppl_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *) 
     endif                      

     if(geq .and. deq) then     
        write (report_unit, *) 
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing gen_eq = .true. AND dfit_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *) 
     endif                      

     if(geq .and. eeq) then     
        write (report_unit, *) 
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing gen_eq = .true. AND efit_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *) 
     endif                      

     if(geq .and. peq) then     
        write (report_unit, *) 
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing gen_eq = .true. AND ppl_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *) 
     endif                      

     if(geq .and. leq) then     
        write (report_unit, *) 
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing gen_eq = .true. AND local_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *) 
     endif                      

     if(eeq .and. deq) then     
        write (report_unit, *) 
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing efit_eq = .true. AND dfit_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *) 
     endif                      

     if(eeq .and. leq) then     
        write (report_unit, *) 
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing efit_eq = .true. AND local_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *) 
     endif                      

     if(eeq .and. peq) then     
        write (report_unit, *) 
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing efit_eq = .true. AND ppl_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *) 
     endif                      

     if(deq .and. leq) then     
        write (report_unit, *) 
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing dfit_eq = .true. AND local_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *) 
     endif                      

     if(deq .and. peq) then     
        write (report_unit, *) 
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing dfit_eq = .true. AND ppl_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *) 
     endif                      

     if(peq .and. leq) then     
        write (report_unit, *) 
        write (report_unit, fmt="('################# WARNING #######################')")
        write(report_unit,fmt="('Choosing ppl_eq = .true. AND local_eq = .true. is not permitted.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *) 
     endif

   end subroutine check

   subroutine factors (n, j, div)
     integer, intent (in) :: n
     integer, intent (out) :: j
     integer, dimension (:), intent (out) :: div
     integer :: i, imax

     do i=2,n
        if (mod(n,i)==0) exit
     end do
     imax = n/i
     j=1
     do i=1,imax
        if (mod(n,i)==0) then
           div(j) = i
           j=j+1
        end if
     end do
     div(j) = n
   end subroutine factors

   subroutine pfactors (n, div)
     integer, intent (in) :: n
     integer, dimension (:), intent (out) :: div
     integer, dimension (50), parameter :: primes = (/ &
          2, 3, 5, 7, 11, &
          13, 17, 19, 23, 29, &
          31, 37, 41, 43, 47, &
          53, 59, 61, 67, 71, &
          73, 79, 83, 89, 97, &
          101, 103, 107, 109, 113, &
          127, 131, 137, 139, 149, &
          151, 157, 163, 167, 173, &
          179, 181, 191, 193, 197, &
          199, 211, 223, 227, 229 /)

     integer :: i, ntmp

     ntmp = n 
     i=1
     do while (ntmp > 1 .and. i < 51)
        do while (mod(ntmp, primes(i)) == 0)
           if (i < 4) div(i) = div(i) + 1
           if (i > 3) div(4) = primes(i)
           ntmp = ntmp / primes(i)
        end do
        i=i+1
     end do
   end subroutine pfactors

   subroutine report

     use theta_grid, only: nbset, ntgrid_real => ntgrid
     implicit none
     real :: zeff_calc, charge, aln, alne, ne, ee, alp, ptot, qsf, dbdr, arat, daky, dtheta0
     real :: kxfac, drhodpsi
     character (20) :: datestamp, timestamp, zone
     character (200) :: line
     logical :: coll_on = .false., le_ok = .true.
     integer :: ntgrid, j, nmesh, npe
     integer, dimension(4) :: pfacs

     call get_unused_unit (report_unit)
     call open_output_file (report_unit, ".report")

     write (report_unit, *) 

     write (report_unit, fmt="('GS2')")
     datestamp(:) = ' '
     timestamp(:) = ' '
     zone(:) = ' '
     call date_and_time (datestamp, timestamp, zone)
     write (report_unit, fmt="('Date: ',a,'/',a,'/',a,&
          &  ' Time: ',a,':',a,1x,a)") datestamp(5:6), &
          &  datestamp(7:8), datestamp(1:4), timestamp(1:2), timestamp(3:4), trim(zone)

     write (report_unit, *) 
     write (report_unit, fmt="('------------------------------------------------------------')")
     write (report_unit, *) 

     if (.not. advanced_egrid) then
        if (.not. any(nesub_ok == nesub)) then
           write (report_unit, *) 
           write (report_unit, fmt="('################# WARNING #######################')")
           write (report_unit, fmt="('You have selected nesub = ',i3)") nesub
           write (report_unit, fmt="('This value is not allowed.')")
           write (report_unit, fmt="('THIS IS AN ERROR.')")
           write (report_unit, fmt="&
                &('The allowed values are: 1-10, 12, 14, 16, 20, 24, 32, 48, or 64.')")
           write (report_unit, fmt="('################# WARNING #######################')")
           write (report_unit, *) 
           le_ok = .false.
        end if
        if (.not. any(nesuper_ok == nesuper)) then
           write (report_unit, *) 
           write (report_unit, fmt="('################# WARNING #######################')")
           write (report_unit, fmt="('You have selected nesuper = ',i3)") nesuper
           write (report_unit, fmt="('This value is not allowed.')")
           write (report_unit, fmt="('THIS IS AN ERROR.')")
           write (report_unit, fmt="('The allowed values are: 1-10, 12, or 15.')")
           write (report_unit, fmt="('################# WARNING #######################')")
           write (report_unit, *) 
           le_ok = .false.
        end if
     else
        if (ecut <= 4.0) then
           write (report_unit, *) 
           write (report_unit, fmt="('################# WARNING #######################')")
           write (report_unit, fmt="('You have selected necut = ',f7.3)") ecut
           write (report_unit, fmt="('With the advanced energy grid, this is small.')")
           write (report_unit, fmt="('THIS IS A PROBABLY AN ERROR.')")
           write (report_unit, fmt="('Recommended values are negrid = 16, ecut = 6.0')")
           write (report_unit, fmt="('################# WARNING #######################')")
           write (report_unit, *) 
        end if
     end if

     if (.not. any(ngauss_ok == ngauss)) then
        write (report_unit, *) 
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, fmt="('You have selected ngauss = ',i3)") ngauss
        write (report_unit, fmt="('This value is not allowed.')")
        write (report_unit, fmt="('THIS IS AN ERROR.')")
        write (report_unit, fmt="&
             &('The allowed values are: 1-6, 8, 10, 12, 16, 20, 24, 32, 40, or 48.')")
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *) 
        le_ok = .false.
     end if

     if (le_ok) then
        if (.not. advanced_egrid) negrid = nesub + nesuper
        if (eps < epsilon(0.0)) trapped_particles = .false.
        if (trapped_particles) then
           nlambda = 2*ngauss + nbset
        else
           nlambda = 2*ngauss
        end if

        write (report_unit, *) 
        write (report_unit, fmt="('Number of lambdas: ',i3)") nlambda
        write (report_unit, *) 

        ntgrid = ntgrid_real
        if (nonlinear_mode_switch == nonlinear_mode_on) then
           if (gridopt_switch /= gridopt_box) then
              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('Nonlinear runs must be carried out in a box.')") 
              write (report_unit, fmt="('Set grid_option to box in the kt_grids_knobs namelist')") 
              write (report_unit, fmt="('or set nonlinear_mode to off in the nonlinear_knobs namelist.')") 
              write (report_unit, fmt="('THIS IS AN ERROR.')") 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, *) 
           end if

           write (report_unit, fmt="('nmesh=(2*ntgrid+1)*2*nlambda*negrid*nx*ny*nspec')")
           nmesh = (2*ntgrid+1)*2*nlambda*negrid*nx*ny*nspec

 !
 ! check that nx, ny have no large prime factors
 !
         pfacs = 0
           call pfactors (ny, pfacs)
           if (pfacs(4) /= 0) then
              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('ny is a multiple of ',i4)") pfacs(4)
              write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
              write (report_unit, fmt="('ny should have only 2, 3, and/or 5 as prime factors.')")
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, *) 
           end if
           i = 1
           if (pfacs(1) > 0) i=2**pfacs(1)
           if (pfacs(2) > 0) i=3**pfacs(2)*i
           if (pfacs(3) > 0) i=5**pfacs(3)*i
           if (i /= ny) then
              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('ny = ',i3)") ny
              write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
              write (report_unit, fmt="('ny should have only 2, 3, and/or 5 as prime factors.')")
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, *) 
           end if

           pfacs = 0
           call pfactors (nx, pfacs)
           if (pfacs(4) /= 0) then
              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('nx is a multiple of ',i3)") pfacs(4)
              write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
              write (report_unit, fmt="('nx should have only 2, 3, and/or 5 as prime factors.')")
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, *) 
           end if
           i = 1
           if (pfacs(1) > 0) i=2**pfacs(1)
           if (pfacs(2) > 0) i=3**pfacs(2)*i
           if (pfacs(3) > 0) i=5**pfacs(3)*i
           if (i /= nx) then
              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('nx = ',i3)") nx
              write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
              write (report_unit, fmt="('nx should have only 2, 3, and/or 5 as prime factors.')")
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, *) 
           end if
        else
           write (report_unit, fmt="('nmesh=(2*ntgrid+1)*2*nlambda*negrid*ntheta0*naky*nspec')")
           nmesh = (2*ntgrid+1)*2*nlambda*negrid*ntheta0*naky*nspec
        end if

        write (report_unit, fmt="('ntgrid :    ',i12)") ntgrid
        write (report_unit, fmt="('nlambda:    ',i12)") nlambda
        write (report_unit, fmt="('negrid :    ',i12)") negrid
        write (report_unit, fmt="('ntheta0:    ',i12)") ntheta0
        write (report_unit, fmt="('naky   :    ',i12)") naky
        write (report_unit, fmt="('ny     :    ',i12)") ny
        write (report_unit, fmt="('nx     :    ',i12)") nx
        write (report_unit, fmt="('nspec  :    ',i12)") nspec
        write (report_unit, fmt="('Number of meshpoints:    ',i12)") nmesh

        call nprocs (nmesh)

     end if

     write (report_unit, *) 
     write (report_unit, fmt="('------------------------------------------------------------')")
     write (report_unit, *) 
     write (report_unit, fmt="('Number of species: ',i3)") nspec
     zeff_calc = 0.
     charge = 0.
     aln = 0.
     alne = 0.
     alp = 0.
     ptot = 0.
     do is=1, nspec
        write (report_unit, *) 
        write (report_unit, fmt="('  Species ',i3)") is
        if (spec(is)%type == 1) write (report_unit, fmt="('    Type:             Ion')")
        if (spec(is)%type == 2) write (report_unit, fmt="('    Type:             Electron')")
        if (spec(is)%type == 3) write (report_unit, fmt="('    Type:             Slowing-down')")
        write (report_unit, fmt="('    Charge:         ',f7.3)") spec(is)%z
        write (report_unit, fmt="('    Mass:             ',es10.4)") spec(is)%mass
        write (report_unit, fmt="('    Density:        ',f7.3)") spec(is)%dens
        write (report_unit, fmt="('    Temperature:    ',f7.3)") spec(is)%temp
        write (report_unit, fmt="('    Collisionality:   ',es10.4)") spec(is)%vnewk
        write (report_unit, fmt="('    Normalized Inverse Gradient Scale Lengths:')")
        write (report_unit, fmt="('      Temperature:  ',f7.3)") spec(is)%tprim
        write (report_unit, fmt="('      Density:      ',f7.3)") spec(is)%fprim
        write (report_unit, fmt="('      Parallel v:   ',f7.3)") spec(is)%uprim
!        write (report_unit, fmt="('    Ignore this:')")
!        write (report_unit, fmt="('    D_0: ',es10.4)") spec(is)%dens0
        if (spec(is)%type /= 2) then
           zeff_calc = zeff_calc + spec(is)%dens*spec(is)%z**2
           charge = charge + spec(is)%dens*spec(is)%z
           aln = aln + spec(is)%dens*spec(is)%z*spec(is)%fprim
        else
           alne = alne + spec(is)%dens*spec(is)%z*spec(is)%fprim
           ne = spec(is)%dens
           ee = spec(is)%z
           has_electrons = .true.
        end if
        alp = alp + spec(is)%dens * spec(is)%temp *(spec(is)%fprim + spec(is)%tprim)
        ptot = ptot + spec(is)%dens * spec(is)%temp
     end do

     if (.not. has_electrons) then
        ptot = ptot + 1./tite   ! electron contribution to pressure
        alp = alp + aln/tite    ! assuming charge neutrality, electron contribution to alp
     end if

     alp = alp / ptot

     write (report_unit, *) 
     write (report_unit, fmt="('------------------------------------------------------------')")

     write (report_unit, fmt="('Calculated Z_eff: ',f7.3)") zeff_calc

     if (has_electrons) then
        if (abs(charge+ne*ee) > 1.e-2) then
           if (charge+ne*ee < 0.) then
              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('You are neglecting an ion species.')")
              write (report_unit, fmt="('This species has a charge fraction of ',f7.3)") abs(charge+ne*ee)
              write (report_unit, &
                   & fmt="('and a normalized inverse density gradient scale length of ',f7.3)") &
                   (aln+alne)/(charge+ne*ee)
              write (report_unit, fmt="('################# WARNING #######################')")
           else
              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('There is an excess ion charge fraction of ',f7.3)") abs(charge+ne*ee)
              write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
              write (report_unit, fmt="('################# WARNING #######################')")
           end if
        else
           if (abs(aln+alne) > 1.e-2) then
              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('The density gradients are inconsistent.')")
              write (report_unit, fmt="('################# WARNING #######################')")
           end if
        end if
     else
        if (charge > 1.01) then
           write (report_unit, *) 
           write (report_unit, fmt="('################# WARNING #######################')")
           write (report_unit, fmt="('There is an excess ion charge fraction of ',f7.3)") charge-1.
           write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
           write (report_unit, fmt="('################# WARNING #######################')")
        end if
        if (charge < 0.99) then
           write (report_unit, *) 
           write (report_unit, fmt="('################# WARNING #######################')")
           write (report_unit, fmt="('You are neglecting an ion species.')")
           write (report_unit, fmt="('This species has a charge fraction of ',f7.3)") abs(charge-1.)
           write (report_unit, fmt="('################# WARNING #######################')")
        end if
     end if

     write (report_unit, *) 
     write (report_unit, fmt="('------------------------------------------------------------')")

     write (report_unit, *) 
     write (report_unit, fmt="('GS2 beta parameter = ',f9.4)") beta
     write (report_unit, fmt="('Total beta = ',f9.4)") beta*ptot
     write (report_unit, *) 
     write (report_unit, fmt="('The total normalized inverse pressure gradient scale length is ',f10.4)") alp
     dbdr = -beta*ptot*alp
     write (report_unit, fmt="('corresponding to d beta / d rho = ',f10.4)") dbdr


     write (report_unit, *) 
     write (report_unit, fmt="('------------------------------------------------------------')")
     write (report_unit, *) 

     select case (eqopt_switch)
     case (eqopt_salpha)
 !
 ! Find q, r/R, R/a
 !
        if (epsl > 0.) then
           arat = 2. / epsl

           if (epsl == 2.0) then
              write (report_unit, &
                   & fmt="('Scale lengths are normalized to the major radius, R')")
           else
              write (report_unit, fmt="('The aspect ratio R/a = ',f7.4)") arat
              if (alne == 1.0) then
                 write (report_unit, &
                      & fmt="('Scale lengths are normalized to the density scale length, Ln')")
              end if
           end if
           qsf = epsl/pk
           write (report_unit, fmt="('The safety factor q =      ',f7.4)") qsf
           write (report_unit, fmt="('The magnetic shear s_hat = ',f7.4)") shat
           if (abs(shat) <= 1.e-5) then
              write (report_unit, fmt="('This is effectively zero; periodic boundary conditions are assumed.')")
           end if
           write (report_unit, fmt="('and epsilon == r/R = ',f7.4)") eps
           write (report_unit, *) 
           if (eps > epsilon(0.0)) then
              write (report_unit, fmt="('Trapped particles are included.')")
           else
              write (report_unit, fmt="('Trapped particles are neglected.')")
           end if
           write (report_unit, *) 

           if (shift > -epsilon(0.0)) then
              write (report_unit, fmt="('The s-alpha alpha parameter is ',f7.4)") shift
!CMR 10/11/06: correct sign of dbeta/drho in s-alpha
              write (report_unit, fmt="('corresponding to d beta / d rho = ',f10.4)") -shift/arat/qsf**2
!CMR 10/11/06: correct sign of dbeta/drho in s-alpha in this check
              if (abs(dbdr + shift/arat/qsf**2) > 1.e-2) then
                 write (report_unit, *) 
                 write (report_unit, fmt="('################# WARNING #######################')")
                 write (report_unit, fmt="('This is inconsistent with beta and the pressure gradient.')") 
                 write (report_unit, fmt="('################# WARNING #######################')")
              end if
           else
              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('The s-alpha alpha parameter is less that zero.')") 
              write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
              write (report_unit, fmt="('################# WARNING #######################')")
           end if

        else
           arat = 1.
           write (report_unit, &
                & fmt="('The radius of curvature is infinite.  This is a slab calculation.')")
        end if

        write (report_unit, *) 
        select case (model_switch)

        case (model_salpha,model_b2,model_eps)
           if (epsl > 0.) then
              write (report_unit, fmt="('An s-alpha model equilibrium has been selected.')")
              write (report_unit, fmt="('The curvature and grad-B drifts are equal.')")
              write (report_unit, *) 
              if (model_switch /= model_eps) then
                 write (report_unit, fmt="('For theta0 = 0, each is of the form')")
                 write (report_unit, *) 
                 write (report_unit, fmt="('  epsl*(cos(theta) + (shat*theta-shift*sin(theta))*sin(theta))')")
                 write (report_unit, *) 
              else
                 write (report_unit, fmt="('For theta0 = 0, each is of the form')")
                 write (report_unit, *) 
                 write (report_unit, fmt="('  epsl*(cos(theta) - eps + (shat*theta-shift*sin(theta))*sin(theta))')")
                 write (report_unit, *) 
              end if
              write (report_unit, fmt="('For finite theta0, there is also a term')")
              write (report_unit, *) 
              write (report_unit, fmt="('  -epsl*shat*sin(theta)*theta0')")
              write (report_unit, *)
           end if
           write (report_unit, *) 
           write (report_unit, fmt="('For theta0 = 0, |(grad S)**2| is of the form')")
           write (report_unit, *) 
           write (report_unit, fmt="('  1.0 + (shat*theta-shift*sin(theta))**2')")
           write (report_unit, *) 
           write (report_unit, fmt="('For finite theta0, there is also a term')")
           write (report_unit, *) 
           write (report_unit, fmt="('  -shat*(shat*theta - shift*sin(theta))*theta0')")
           write (report_unit, *) 
           write (report_unit, fmt="('and finally, the term')")
           write (report_unit, *) 
           write (report_unit, fmt="('  shat**2 * theta0**2')")
           write (report_unit, *) 
           if (model_switch == model_eps) then
              write (report_unit, *) 
              write (report_unit, fmt="(' This model differs from the normal s-alpha model')") 
              write (report_unit, fmt="(' only in the curv and grad_B drifts.')")
           end if
           if (model_switch == model_b2) then
              write (report_unit, *) 
              write (report_unit, fmt="(' This model differs from the normal s-alpha model')") 
              write (report_unit, fmt="(' by an additional factor of 1/B(theta)**2 (not shown above)')")
              write (report_unit, fmt="(' in the curv and grad_B drifts.')")
           end if
        case (model_ccurv)
           write (report_unit, fmt="('Constant curvature is assumed.')")
           write (report_unit, fmt="('The grad-B and curvature drifts are each = ',f10.4)") epsl
           write (report_unit, *) 
           write (report_unit, fmt="('For theta0 = 0, |(grad S)**2| is of the form')")
           write (report_unit, *) 
           write (report_unit, fmt="('  1.0 + (shat*theta-shift*sin(theta))**2')")
           write (report_unit, *) 
           write (report_unit, fmt="('For finite theta0, there is also a term')")
           write (report_unit, *) 
           write (report_unit, fmt="('  -shat*shat*theta*theta0')")
           write (report_unit, *) 
           write (report_unit, fmt="('and finally, the term')")
           write (report_unit, *) 
           write (report_unit, fmt="('  shat**2 * theta0**2')")
           write (report_unit, *) 
        case (model_nocurve)
           write (report_unit, fmt="('Zero curvature is assumed.')")
           write (report_unit, *) 
           write (report_unit, fmt="('For theta0 = 0, |(grad S)**2| is of the form')")
           write (report_unit, *) 
           write (report_unit, fmt="('  1.0 + (shat*theta)**2')")
           write (report_unit, *) 
           write (report_unit, fmt="('For finite theta0, there is also a term')")
           write (report_unit, *) 
           write (report_unit, fmt="('  -shat*shat*theta*theta0')")
           write (report_unit, *) 
           write (report_unit, fmt="('and finally, the term')")
           write (report_unit, *) 
           write (report_unit, fmt="('  shat**2 * theta0**2')")
           write (report_unit, *) 
        end select

     case (eqopt_eik)

        call check (vmom_eq, gen_eq, efit_eq, ppl_eq, local_eq, dfit_eq, idfit_eq, report_unit)
        write (report_unit, *)
        if (local_eq .and. iflux == 0) then
           write (report_unit, fmt="('A local equilibrium model has been selected.')")
           if (Rmaj == 1.0) then
              write (report_unit, &
                   & fmt="('Scale lengths are normalized to the major radius, R')")
           else
              write (report_unit, fmt="('The aspect ratio R/a = ',f7.4)") Rmaj
           end if
           if (Rmaj /= R_geo) then
              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('R_geo is not equal to Rmaj.')")
              write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
              write (report_unit, fmt="('################# WARNING #######################')")
           end if
           if (irho /= 2) then
              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('You have selected irho = ',i2)") irho
              write (report_unit, fmt="('For local equilibria, irho=2 is required.')")
              write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
              write (report_unit, fmt="('################# WARNING #######################')")
           end if
           write (report_unit, *) 
           write (report_unit, fmt="('The safety factor q =      ',f7.4)") qinp
           eps = rhoc/R_geo
           write (report_unit, fmt="('and epsilon == r/R =       ',f7.4)") eps
           write (report_unit, *) 
           if (eps > epsilon(0.0)) then
              write (report_unit, fmt="('Trapped particles are included.')")
           else
              write (report_unit, fmt="('Trapped particles are neglected.')")
           end if
           write (report_unit, *) 
           write (report_unit, fmt="('B_poloidal is determined by:')")
           write (report_unit, *) 
           write (report_unit, fmt="('    triangularity, tri =       ',f7.4)") tri
           write (report_unit, fmt="('  & gradient: d tri /d rho =   ',f7.4)") tripri
           write (report_unit, *) 
           write (report_unit, fmt="('    elongation, kappa =        ',f7.4)") akappa
           write (report_unit, fmt="('  & gradient: d kappa /d rho = ',f7.4)") akappri

           write (report_unit, *) 
           write (report_unit, fmt="('The magnetic shear s_hat = ',f7.4)") shat
           write (report_unit, fmt="('This value is set by s_hat_input in the theta_grid_eik_knobs namelist.')") 
           if (abs(shat) <= 1.e-5) then
              write (report_unit, fmt="('This is effectively zero; periodic boundary conditions are assumed.')")
           end if
           select case (bishop)
           case (3) 
              write (report_unit, fmt="('The normalized inverse pressure gradient scale length = ',f8.4)") invLp_input
           case (4) 
              write (report_unit, fmt="('The beta gradient d beta / d rho = ',f8.4)") beta_prime_input
              if (beta_prime_input > epsilon(0.0)) then
                 write (report_unit, *) 
                 write (report_unit, fmt="('################# WARNING #######################')")
                 write (report_unit, fmt="('beta_prime > 0.')")
                 write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
                 write (report_unit, fmt="('################# WARNING #######################')")
                 write (report_unit, *) 
              end if
              if (abs(beta_prime_input - dbdr) > 1.e-2) then
                 write (report_unit, *) 
                 write (report_unit, fmt="('################# WARNING #######################')")
                 write (report_unit, fmt="('beta_prime_input is not consistent with beta and Lp.')")
                 write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
                 write (report_unit, fmt="('################# WARNING #######################')")
                 write (report_unit, *) 
              end if
           case (5) 
              write (report_unit, fmt="('The alpha parameter (R beta_prime q**2) = ',f8.4)") alpha_input
!              write (*,*) alpha_input, dbdr, qinp, Rmaj
              if (abs(alpha_input + dbdr*qinp**2*Rmaj) > 1.e-2) then
                 write (report_unit, *) 
                 write (report_unit, fmt="('################# WARNING #######################')")
                 write (report_unit, fmt="('alpha is not consistent with beta, q, and Lp.')")
                 write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
                 write (report_unit, fmt="('################# WARNING #######################')")
                 write (report_unit, *) 
              end if
           case default
              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('You have selected bishop = ',i2)") bishop
              write (report_unit, fmt="('For local equilibria, bishop = 4 is recommended.')")
              if (bishop == 1) then
                 write (report_unit, fmt="('For d beta / d rho = 0, bishop = 1 is ok.')")
                 write (report_unit, fmt="('Otherwise, ')")
              end if
              write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, *) 
           end select
        end if
        if (local_eq .and. .not. (iflux == 0)) then
           write (report_unit, *) 
           write (report_unit, fmt="('################# WARNING #######################')")
           write (report_unit, fmt="('You have selected a local equilibrium and iflux = ',i2)") iflux
           write (report_unit, fmt="('For local equilibria, iflux=0 is required.')")
           write (report_unit, fmt="('THIS IS AN ERROR.')") 
           write (report_unit, fmt="('################# WARNING #######################')")
           write (report_unit, *) 
        end if
        if (.not. local_eq) then
           if (gen_eq) then
              write (report_unit, *) 
              write (report_unit, fmt="('Equilibrium information obtained from NetCDF file:')")
              write (report_unit, fmt="(a)") trim(eqfile)
           end if
           if (vmom_eq) then
              write (report_unit, *) 
              write (report_unit, fmt="('Equilibrium information obtained from file:')")
              write (report_unit, fmt="(a)") trim(eqfile)
           end if
           if (ppl_eq) then
              write (report_unit, *) 
              write (report_unit, fmt="('Equilibrium information obtained from NetCDF file:')")
              write (report_unit, fmt="(a)") trim(eqfile)
           end if
           if (dfit_eq) then
              write (report_unit, *) 
              write (report_unit, fmt="('Dipole equilibrium information obtained from file:')")
              write (report_unit, fmt="(a)") trim(eqfile)
           end if
           if (idfit_eq) then
              write (report_unit, *) 
              write (report_unit, fmt="('Dipole equilibrium information obtained from file:')")
              write (report_unit, fmt="(a)") trim(eqfile)
           end if
 !          if (mds) then
 !          end if
           if (efit_eq) then
              write (report_unit, *) 
              write (report_unit, fmt="('Equilibrium information obtained from eqdsk:')")
              write (report_unit, fmt="(a)") trim(eqfile)
           end if
           select case (bishop)
           case (1) 
              write (report_unit, *) 
              write (report_unit, fmt="('You have set bishop=1, so dp/drho and s_hat will be found from the equilibrium file.')")
              write (report_unit, *) 
           case (3) 
              write (report_unit, *) 
              write (report_unit, fmt="('You have set bishop=3.')")
              write (report_unit, *) 
              write (report_unit, fmt="('The magnetic shear s_hat = ',f7.4)") s_hat_input
              write (report_unit, fmt="('This value is set by s_hat_input in the theta_grid_eik_knobs namelist.')") 
              if (abs(shat) <= 1.e-5) then
                 write (report_unit, fmt="('This is effectively zero; periodic boundary conditions are assumed.')")
              end if
              write (report_unit, fmt="('The normalized inverse pressure gradient scale length = ',f8.4)") invLp_input
           case (4) 
              write (report_unit, *) 
              write (report_unit, fmt="('You have set bishop=4.')")
              write (report_unit, *) 
              write (report_unit, fmt="('The magnetic shear s_hat = ',f7.4)") s_hat_input
              write (report_unit, fmt="('This value is set by s_hat_input in the theta_grid_eik_knobs namelist.')") 
              if (abs(shat) <= 1.e-5) then
                 write (report_unit, fmt="('This is effectively zero; periodic boundary conditions are assumed.')")
              end if
              write (report_unit, fmt="('The beta gradient d beta / d rho = ',f8.4)") beta_prime_input
              if (beta_prime_input > epsilon(0.0)) then
                 write (report_unit, *) 
                 write (report_unit, fmt="('################# WARNING #######################')")
                 write (report_unit, fmt="('beta_prime > 0.')")
                 write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
                 write (report_unit, fmt="('################# WARNING #######################')")
                 write (report_unit, *) 
              end if
              if (abs(beta_prime_input - dbdr) > 1.e-2) then
                 write (report_unit, *) 
                 write (report_unit, fmt="('################# WARNING #######################')")
                 write (report_unit, fmt="('beta_prime_input is not consistent with beta and Lp.')")
                 write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
                 write (report_unit, fmt="('################# WARNING #######################')")
                 write (report_unit, *) 
              end if
           case (5) 
              write (report_unit, *) 
              write (report_unit, fmt="('You have set bishop=5.')")
              write (report_unit, *) 
              write (report_unit, fmt="('The magnetic shear s_hat = ',f7.4)") s_hat_input
              write (report_unit, fmt="('This value is set by s_hat_input in the theta_grid_eik_knobs namelist.')") 
              if (abs(shat) <= 1.e-5) then
                 write (report_unit, fmt="('This is effectively zero; periodic boundary conditions are assumed.')")
              end if
              write (report_unit, fmt="('The alpha parameter (R beta_prime q**2) = ',f8.4)") alpha_input
              write (*,*) alpha_input, dbdr, qinp, Rmaj
              if (abs(alpha_input + dbdr*qinp**2*Rmaj) > 1.e-2) then
                 write (report_unit, *) 
                 write (report_unit, fmt="('################# WARNING #######################')")
                 write (report_unit, fmt="('alpha is not consistent with beta, q, and Lp.')")
                 write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
                 write (report_unit, fmt="('################# WARNING #######################')")
                 write (report_unit, *) 
              end if
           case (6) 
              write (report_unit, *) 
              write (report_unit, fmt="('You have set bishop=6.')")
              write (report_unit, *) 
              write (report_unit, fmt="('The magnetic shear s_hat = ',f7.4)") s_hat_input
              write (report_unit, fmt="('This value is set by s_hat_input in the theta_grid_eik_knobs namelist.')") 
              if (abs(shat) <= 1.e-5) then
                 write (report_unit, fmt="('This is effectively zero; periodic boundary conditions are assumed.')")
              end if
              write (report_unit, fmt="('The value of dp/drho will be found from the equilibrium file.')") 
           case (7) 
              write (report_unit, *) 
              write (report_unit, fmt="('You have set bishop=7.')")
              write (report_unit, fmt="('The value of s_hat will be found from the equilibrium file.')") 
              write (report_unit, fmt="('The value of dp/drho found from the equilibrium file will be multiplied by',f10.4)") &
	           dp_mult
           case default

              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('You have selected a value for bishop that is not recommended.')")
              write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, *) 

           end select
        end if
     case (eqopt_file)
        write (report_unit, *) 
        write (report_unit, fmt="('Equilibrium information obtained from gridgen output file:')")
        write (report_unit, fmt="(a)") trim(gridout_file)

        call get_unused_unit (iunit)
        open (unit=iunit, file=gridout_file, status="old", err=100)
        read (unit=iunit, fmt="(a)") line
        read (unit=iunit, fmt=*) nbset
        read (unit=iunit, fmt="(a)") line
        do i = 1, nbset
           read (unit=iunit, fmt="(a)") line
        end do

        read (unit=iunit, fmt="(a)") line
        read (unit=iunit, fmt=*) ntgrid, nperiod, ntheta, &
             drhodpsi, rmaj, shat, kxfac

        close (unit=iunit)

        write (report_unit, *) 
        write (report_unit, fmt="('Limited information available:')")
        write (report_unit, *) 
        write (report_unit, fmt="('nbset =     ',i5)") nbset
        write (report_unit, fmt="('ntgrid =    ',i5)") ntgrid
        write (report_unit, fmt="('ntheta =    ',i5)") ntheta
        write (report_unit, fmt="('nperiod =   ',i2)") nperiod
        write (report_unit, fmt="('drhodpsi =  ',f8.4)") drhodpsi
        write (report_unit, fmt="('R =         ',f8.4)") Rmaj
        write (report_unit, fmt="('s_hat =     ',f8.4)") shat
        write (report_unit, fmt="('kxfac =     ',f8.4)") kxfac

        write (report_unit, *)
        write (report_unit, *) 'NOTE: Regardless of the values of ntheta and nperiod'
        write (report_unit, *) '      found in the theta_grid_parameters namelist,'
        write (report_unit, *) '      this calculation will use the values listed here:'
        write (report_unit, fmt="('ntgrid =    ',i5)") ntgrid
        write (report_unit, fmt="('ntheta =    ',i5)") ntheta
        write (report_unit, *) '      These were obtained from the gridgen output file.'
        write (report_unit, *)
 100    continue

     end select

     if (gb_to_cv) then
        write (report_unit, *) 'The grad B drift coefficients have been set equal to the'
        write (report_unit, *) 'values for the curvature drift coefficients.  Do not use'
        write (report_unit, *) 'fbpar = 1.0 in this case.'
        write (report_unit, *)
        write (report_unit, *) 'You got this option by setting gb_to_cv = .true.'
        write (report_unit, *) 
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, fmt="('You have chosen to set the grad B drift equal to the curvature drift.')")
        write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
        write (report_unit, fmt="('################# WARNING #######################')")
        write (report_unit, *) 
     end if

     write (report_unit, *) 
     write (report_unit, fmt="('------------------------------------------------------------')")
     write (report_unit, *) 

     do is=1, nspec
        coll_on = spec(is)%vnewk > epsilon(0.0) .or. coll_on
     end do

     if (coll_on) then
        select case (collision_model_switch)
        case (collision_model_lorentz,collision_model_lorentz_test)
           write (report_unit, fmt="('A Lorentz collision operator has been selected.')")
           if (cfac > 0) write (report_unit, fmt="('This has both terms of Lorentz collision operator (cfac=1.0)')")
           if (cfac == 0) write (report_unit, fmt="('This is only a partial Lorentz collision operator (cfac=0.0)')")
           if (const_v) write (report_unit, fmt="('This is an energy independent Lorentz collision operator (const_v=true)')")  
 !          if (hypercoll) call init_hyper_lorentz
       case (collision_model_krook,collision_model_krook_test)
          write (report_unit, fmt="('A Krook collision operator has been selected.')")
       end select
    else
       write (report_unit, fmt="('All collisionality parameters (vnewk) are zero.')")
       write (report_unit, fmt="('No collision operator will be used.')")
    end if
       
    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 

    select case (fieldopt_switch)
    case (fieldopt_implicit)
       write (report_unit, fmt="('The field equations will be advanced in time implicitly.')")
    case (fieldopt_explicit)
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The field equations will be advanced in time explicitly.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    case (fieldopt_test)
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The field equations will only be tested.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end select

! 
! implicitness parameters
!
    do is = 1, nspec
       if (aimag(fexp(is)) /= 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Species ',i2,' has fexpi = ',e10.4)") is, aimag(fexp(is))
          write (report_unit, fmt="('THIS IS AN ERROR')")
          write (report_unit, fmt="('fexpi should be zero for all species.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if

       write (report_unit, fmt="('Species ',i2,' has fexpr = ', e10.4)") is, real(fexp(is))
    end do
    
    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 

    select case (ginitopt_switch)
    case (ginitopt_default)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, fmt="('  Amplitude:        ',f10.4)") phiinit
       write (report_unit, fmt="('  Width in theta:   ',f10.4)") width0
       if (chop_side) then
          write (report_unit, fmt="('  Parity:   none')") 
       else
          write (report_unit, fmt="('  Parity:   even')") 
       end if

    case (ginitopt_kz0)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, fmt="('  Amplitude:        ',f10.4)") phiinit
       write (report_unit, fmt="('  Constant along field line',f10.4)") width0
       if (chop_side) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('  Parity:   none')") 
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
          write (report_unit, fmt="('Remedy: set chop_side = .false. in init_g_knobs.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if

    case (ginitopt_noise)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, fmt="('  Amplitude:        ',f10.4)") phiinit
       write (report_unit, fmt="('  Noise along field line.')") 
       if (zf_init /= 1.) then
          write (report_unit, fmt="('  Zonal flows adjusted by factor of zf_init = ',f10.4)") zf_init
       end if

    case (ginitopt_kpar)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, fmt="('  Amplitude:             ',f10.4)") phiinit
       write (report_unit, fmt="('  Real part multiplier:  ',f10.4)") refac
       write (report_unit, fmt="('  Imag part multiplier:  ',f10.4)") imfac
       if (width0 > 0.) then
          write (report_unit, fmt="('  Gaussian envelope in theta with width:  ',f10.4)") width0
       end if
       if (chop_side) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('  Parity:   none')") 
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
          write (report_unit, fmt="('Remedy: set chop_side = .false. in init_g_knobs.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if
       if (den0 > epsilon(0.0) .or. den1 > epsilon(0.0) .or. den2 > epsilon(0.0)) then
          write (report_unit, *) 
          write (report_unit, fmt="('Initial density perturbation of the form:')")
          write (report_unit, fmt="('den0   + den1 * cos(theta) + den2 * cos(2.*theta)')")
          write (report_unit, *) 
          write (report_unit, fmt="('with den0 =',f7.4,' den1 = ',f7.4,' den2 = ',f7.4)") den0, den1, den2
       end if
       if (upar0 > epsilon(0.0) .or. upar1 > epsilon(0.0) .or. upar2 > epsilon(0.0)) then
          write (report_unit, *) 
          write (report_unit, fmt="('Initial parallel velocity perturbation of the form:')")
          write (report_unit, fmt="('upar0   + upar1 * cos(theta) + upar2 * cos(2.*theta)')")
          write (report_unit, fmt="('90 degrees out of phase with other perturbations.')")
          write (report_unit, *) 
          write (report_unit, fmt="('with upar0 =',f7.4,' upar1 = ',f7.4,' upar2 = ',f7.4)") upar0, upar1, upar2
       end if
       if (tpar0 > epsilon(0.0) .or. tpar1 > epsilon(0.0) .or. tpar2 > epsilon(0.0)) then
          write (report_unit, *) 
          write (report_unit, fmt="('Initial Tpar perturbation of the form:')")
          write (report_unit, fmt="('tpar0   + tpar1 * cos(theta) + tpar2 * cos(2.*theta)')")
          write (report_unit, *) 
          write (report_unit, fmt="('with tpar0 =',f7.4,' tpar1 = ',f7.4,' tpar2 = ',f7.4)") tpar0, tpar1, tpar2
       end if
       if (tperp0 > epsilon(0.0) .or. tperp1 > epsilon(0.0) .or. tperp2 > epsilon(0.0)) then
          write (report_unit, *) 
          write (report_unit, fmt="('Initial Tperp perturbation of the form:')")
          write (report_unit, fmt="('tperp0   + tperp1 * cos(theta) + tperp2 * cos(2.*theta)')")
          write (report_unit, *) 
          write (report_unit, fmt="('with tperp0 =',f7.4,' tperp1 = ',f7.4,' tperp2 = ',f7.4)") tperp0, tperp1, tperp2
       end if
       if (has_electrons) then
          write (report_unit, *) 
          write (report_unit, fmt="('Field line average of g_electron subtracted off.')")
       end if

    case (ginitopt_gs)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, fmt="('  Randomly phased kpar=1 sines and cosines')") 
       write (report_unit, fmt="('  in density, upar, tpar, or tperp.')") 
       write (report_unit, fmt="('  Real part amplitude:  ',f10.4)") refac*phiinit
       write (report_unit, fmt="('  Imag part amplitude:  ',f10.4)") imfac*phiinit
       if (abs( den1)  > epsilon(0.0)) write (report_unit, fmt="('  Density amplitude:  ',f10.4)") den1
       if (abs( upar1) > epsilon(0.0)) write (report_unit, fmt="('  Upar amplitude:  ',f10.4)") upar1
       if (abs( tpar1) > epsilon(0.0)) write (report_unit, fmt="('  Tpar amplitude:  ',f10.4)") tpar1
       if (abs(tperp1) > epsilon(0.0)) write (report_unit, fmt="('  Tperp amplitude:  ',f10.4)") tperp1

    case (ginitopt_nl)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, fmt="('At most two k_perps excited, with amplitude = ',f10.4)") phiinit
       write (report_unit, fmt="(' First k_perp has ik = ',i3,' it = ',i3)") ikk(1), itt(1)
       write (report_unit, fmt="('Second k_perp has ik = ',i3,' it = ',i3)") ikk(2), itt(2)
       if (chop_side) then
          write (report_unit, fmt="('  Parity:   none')") 
       else
          write (report_unit, fmt="('  Parity:   even')") 
       end if
       write (report_unit, fmt="('Reality condition is enforced.')")
       
    case (ginitopt_nl2)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, fmt="('At most two k_perps excited, with amplitude = ',f10.4)") phiinit
       write (report_unit, fmt="(' First k_perp has ik = ',i3,' it = ',i3)") ikk(1), itt(1)
       write (report_unit, fmt="('Second k_perp has ik = ',i3,' it = ',i3)") ikk(2), itt(2)
       if (chop_side) then
          write (report_unit, fmt="('  Parity:   none')") 
       else
          write (report_unit, fmt="('  Parity:   even')") 
       end if
       write (report_unit, fmt="('Reality condition is enforced.')")
       write (report_unit, fmt="('g perturbation proportional to (1+v_parallel)*sin(theta)')")

    case (ginitopt_nl3)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, fmt="('At most two k_perps excited, with amplitude = ',f10.4)") phiinit
       write (report_unit, fmt="(' First k_perp has ik = ',i3,' it = ',i3)") ikk(1), itt(1)
       write (report_unit, fmt="('Second k_perp has ik = ',i3,' it = ',i3)") ikk(2), itt(2)
       write (report_unit, fmt="('  Real part multiplied by:  ',f10.4)") refac
       write (report_unit, fmt="('  Imag part multiplied by:  ',f10.4)") imfac
       if (width0 > 0.) then
          write (report_unit, fmt="('  Gaussian envelope in theta with width:  ',f10.4)") width0
       end if
       if (chop_side) then
          write (report_unit, fmt="('  Parity:   none')") 
       else
          write (report_unit, fmt="('  Parity:   even')") 
       end if
       write (report_unit, fmt="('Reality condition is enforced.')")
       if (den0 > epsilon(0.0) .or. den1 > epsilon(0.0) .or. den2 > epsilon(0.0)) then
          write (report_unit, *) 
          write (report_unit, fmt="('Initial density perturbation of the form:')")
          write (report_unit, fmt="('den0   + den1 * cos(theta) + den2 * cos(2.*theta)')")
          write (report_unit, *) 
          write (report_unit, fmt="('with den0 =',f7.4,' den1 = ',f7.4,' den2 = ',f7.4)") den0, den1, den2
       end if
       if (upar0 > epsilon(0.0) .or. upar1 > epsilon(0.0) .or. upar2 > epsilon(0.0)) then
          write (report_unit, *) 
          write (report_unit, fmt="('Initial parallel velocity perturbation of the form:')")
          write (report_unit, fmt="('upar0   + upar1 * cos(theta) + upar2 * cos(2.*theta)')")
          write (report_unit, fmt="('90 degrees out of phase with other perturbations.')")
          write (report_unit, *) 
          write (report_unit, fmt="('with upar0 =',f7.4,' upar1 = ',f7.4,' upar2 = ',f7.4)") upar0, upar1, upar2
       end if
       if (tpar0 > epsilon(0.0) .or. tpar1 > epsilon(0.0) .or. tpar2 > epsilon(0.0)) then
          write (report_unit, *) 
          write (report_unit, fmt="('Initial Tpar perturbation of the form:')")
          write (report_unit, fmt="('tpar0   + tpar1 * cos(theta) + tpar2 * cos(2.*theta)')")
          write (report_unit, *) 
          write (report_unit, fmt="('with tpar0 =',f7.4,' tpar1 = ',f7.4,' tpar2 = ',f7.4)") tpar0, tpar1, tpar2
       end if
       if (tperp0 > epsilon(0.0) .or. tperp1 > epsilon(0.0) .or. tperp2 > epsilon(0.0)) then
          write (report_unit, *) 
          write (report_unit, fmt="('Initial Tperp perturbation of the form:')")
          write (report_unit, fmt="('tperp0   + tperp1 * cos(theta) + tperp2 * cos(2.*theta)')")
          write (report_unit, *) 
          write (report_unit, fmt="('with tperp0 =',f7.4,' tperp1 = ',f7.4,' tperp2 = ',f7.4)") tperp0, tperp1, tperp2
       end if
       
    case (ginitopt_nl4)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('Under development for study of secondary instabilities.')")
       write (report_unit, fmt="('Scale factor:   ',f10.4)") scale
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    case (ginitopt_nl5)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('Under development for study of secondary instabilities.')")
       write (report_unit, fmt="('Scale factor:   ',f10.4)") scale
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    case (ginitopt_nl6)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('Change amplitude of a particular mode.')")
       write (report_unit, fmt="('Scale factor:   ',f10.4)") scale
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    case (ginitopt_test1)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('Maxwellian with sin(kr * theta)/(i*kr), amplitude = ',f10.4)") phiinit
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    case (ginitopt_xi)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('Perturbation proportional to pitch angle.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    case (ginitopt_xi2)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('Perturbation proportional to function of pitch angle.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    case (ginitopt_rh)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('Maxwellian perturbation in ik=1 mode.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    case (ginitopt_alf)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('Ion dist fn proportional to v_parallel * sin(theta).')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    case (ginitopt_zero)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('Distribution function = 0.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    case (ginitopt_test3)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    case (ginitopt_convect)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    case (ginitopt_restart_file)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('Restart from a single NetCDF restart file.')") 
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    case (ginitopt_restart_many)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, fmt="('Each PE restarts from its own NetCDF restart file.')") 

    case (ginitopt_restart_small)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, fmt="('Each PE restarts from its own NetCDF restart file.')") 
       write (report_unit, fmt="('with amplitudes scaled by factor of scale = ',f10.4)") scale
       write (report_unit, fmt="('Noise added with amplitude = ',f10.4)") phiinit

    case (ginitopt_continue)
       write (report_unit, fmt="('Initial conditions:')")
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    end select

    if (ginitopt_switch == ginitopt_restart_many) then
       if (delt_option_switch == delt_option_auto) then
          write (report_unit, *) 
          write (report_unit, fmt="('This run is a continuation of a previous run.')") 
          write (report_unit, fmt="('The time step at the beginning of this run')") 
          write (report_unit, fmt="('will be taken from the end of the previous run.')") 
       else
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('This run is a continuation of a previous run.')") 
          write (report_unit, fmt="('The time step is being set by hand.')") 
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
          write (report_unit, fmt="('You probably want to set delt_option to be check_restart in the knobs namelist.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if
    end if

    if (delt_option_switch == delt_option_auto) then
       if (ginitopt_switch /= ginitopt_restart_many) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('This is not a normal continuation run.')") 
          write (report_unit, fmt="('You probably want to set delt_option to be default in the knobs namelist.')") 
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if
    end if
       
    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 

    if (nonlinear_mode_switch == nonlinear_mode_on) then
       write (report_unit, *) 
       write (report_unit, fmt="('This is a nonlinear simulation.')")
       write (report_unit, *) 
       if (wstar_units) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Nonlinear runs require wstar_units = .false. in the knobs namelist.')") 
          write (report_unit, fmt="('THIS IS AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if
       if (gridopt_switch /= gridopt_box) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Nonlinear runs must be carried out in a box.')") 
          write (report_unit, fmt="('Set grid_option to box in the kt_grids_knobs namelist.')") 
          write (report_unit, fmt="('THIS IS AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if

       write (report_unit, fmt="('The minimum delt = ',e10.4)") delt_minimum
       write (report_unit, fmt="('The maximum delt = ',e10.4)") delt
       write (report_unit, fmt="('The maximum delt < ',f10.4,' * min(Delta_perp/v_perp.')") cfl
       write (report_unit, fmt="('When the time step needs to be changed, it is adjusted by a factor of ',f10.4)") delt_adj
       write (report_unit, fmt="('The number of time steps nstep = ',i7)") nstep
       write (report_unit, fmt="('If running in batch mode on the NERSC T3E, the run will stop when ', &
            & f6.4,' % of the time remains.')") 100.*margin

       if (nperiod > 1) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Nonlinear runs usually have nperiod = 1.')") 
          write (report_unit, fmt="('THIS MAY BE AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if

    else
       write (report_unit, *) 
       write (report_unit, fmt="('This is a linear calculation.')")
       write (report_unit, *)        
       write (report_unit, fmt="('The time step (delta t) = ',e10.4)") delt
       write (report_unit, fmt="('The maximum number of time steps is ',i7)") nstep
       if (exit_when_converged) then
          write (report_unit, fmt="('When the frequencies for each k have converged, the run will stop.')")
          write (report_unit, fmt="('The convergence has to be better than one part in ',e10.4)") 1./omegatol
       end if

       if (wstar_units) then
          write (report_unit, *) 
          write (report_unit, fmt="('The timestep for each ky is scaled by a factor of 1/ky.')")
       end if
    end if

    write (report_unit, *) 
    write (report_unit, fmt="('Data will be written to ',a,' every ',i4,' timesteps.')") trim(run_name)//'.out.nc', nwrite
    write (report_unit, *) 

    if (dump_fields_periodically) then
       write (report_unit, *) 
       write (report_unit, fmt="('Data will be written to dump.fields.t=* every ',i4,' timesteps.')") 10*nmovie
       write (report_unit, *) 
    end if

    if (make_movie) then
       write (report_unit, *) 
       write (report_unit, fmt="('Movie data will be written to runname.movie.nc every ',i4,' timesteps.')") nmovie
       write (report_unit, *) 
    end if

    if (save_for_restart .and. nsave > 0) then
       write (report_unit, *) 
       write (report_unit, fmt="('Restart data will be written every ',i4,' timesteps.')") nsave
       write (report_unit, *) 
    end if
    
    if (flow_mode_switch == flow_mode_on) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('flow_mode=on is not allowed.  Flow mode is buggy.')") 
       write (report_unit, fmt="('THIS IS AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 

    select case (gridopt_switch) 
    case (gridopt_single) 
       write (report_unit, *) 
       write (report_unit, fmt="('A single k_perp will be evolved, with: ')")
       write (report_unit, *) 
       write (report_unit, fmt="('ky rho = ',f10.4)") aky
       write (report_unit, fmt="('theta0 = ',f10.4)") theta0
       if (akx /= 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('The value of akx in the kt_grids_single_parameters namelist is ignored.')") 
          write (report_unit, fmt="('You have set akx to a non-zero value.')") 
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if

    case (gridopt_range)
       write (report_unit, *) 
       write (report_unit, fmt="('A range of k_perps will be evolved.')")
       write (report_unit, *) 
       write (report_unit, fmt="('There are ',i3,' values of ky rho and ',i3,' values of theta0:')") naky, ntheta0
       write (report_unit, *) 
          
       daky = 0.0
       if (naky > 1) daky = (aky_max - aky_min)/real(naky - 1)
       dtheta0 = 0.0
       if (ntheta0 > 1) dtheta0 = (theta0_max - theta0_min)/real(ntheta0 - 1)

       do j = 0, naky-1
          do i = 0, ntheta0-1
             write (report_unit, fmt="('ky rho = ',e10.4,' theta0 = ',e10.4)") &
                  aky_min + daky*real(j), theta0_min + dtheta0*real(i)
          end do
       end do

    case (gridopt_specified) 

       write (report_unit, *) 
       i = max (naky, ntheta0)
       write (report_unit, fmt="('A set of ',i3,' k_perps will be evolved.')") i
       write (report_unit, *) 
       do i=1, max(naky,ntheta0)
          write (report_unit, fmt="('ky rho = ',e10.4,' theta0 = ',e10.4)") &
               aky_tmp(i), theta0_tmp(i)
       end do

    case (gridopt_box)

       if (y0_internal /= 2.) then
          if (abs(2.*pi*y0_internal - ly) > epsilon(0.)) then
             write (report_unit, *) 
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('You cannot specify both ly and y0.')")
             write (report_unit, fmt="('THIS IS AN ERROR.')") 
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 
          end if
       end if

       write (report_unit, *) 
       write (report_unit, fmt="('A rectangular simulation domain has been selected.')")
       write (report_unit, *) 
       write (report_unit, fmt="('The domain is ',f10.4,' rho in the y direction.')") ly

       if (abs(shat) <= 1.e-5) then
          if (x0 == 0.) then
             if (rtwist > 0) then
                write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)") abs(rtwist)*ly
             else
                write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)") ly/abs(rtwist)
             end if
          else
             if (x0 > 0) then
                write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)") 2.*pi*x0
             else
                write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)") -2.*pi/x0
             end if
          end if
       else
          lx = ly * rtwist / (2.*pi*shat)
          write (report_unit, fmt="('At theta=0, the domain is ',f10.4,' rho in the x direction.')") lx
       end if
       
       write (report_unit, *) 
       write (report_unit, fmt="('The nonlinear terms will be evaluated on a grid with ',&
            & i4,' points in x and ',i4,' points in y.')") nx, ny
       write (report_unit, *) 
       naky = (ny-1)/3+1
       ntheta0 = 2*((nx-1)/3)+1
       write (report_unit, fmt="('After de-aliasing, there will be ',i4,'  ky >= 0 modes and ',i4,' kx modes.')") naky, ntheta0
       write (report_unit, fmt="('The modes with ky < 0 are determined by the reality condition.')")

    case (gridopt_xbox)
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('You selected grid_option=xbox in kt_grids_knobs')")
          write (report_unit, fmt="('The xbox option is not working.')")
          write (report_unit, fmt="('THIS IS AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 

    end select

    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 

    if (gridfac /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('You selected gridfac = ',e10.4,' in dist_fn_knobs.')") gridfac
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('The normal choice is gridfac = 1.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (apfac /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('You selected apfac = ',e10.4,' in dist_fn_knobs.')") apfac
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('The normal choice is apfac = 1.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (driftknob /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('You selected driftknob = ',e10.4,' in dist_fn_knobs.')") driftknob
       write (report_unit, fmt="('THIS IS EITHER AN ERROR, or you are DELIBERATELY SCALING THE DRIFTS.')") 
       write (report_unit, fmt="('The normal choice is driftknob = 1.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (tpdriftknob /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('You selected tpdriftknob = ',e10.4,' in dist_fn_knobs.')") tpdriftknob
       write (report_unit, fmt="('THIS IS EITHER AN ERROR, or you are DELIBERATELY SCALING THE TRAPPED PARTICLE DRIFTS (either via driftknob or via tpdriftknob).')") 
       write (report_unit, fmt="('The normal choice is tpdriftknob = 1.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    select case (boundary_option_switch)
    case (boundary_option_linked)
       write (report_unit, *) 
       if (gridopt_switch /= gridopt_box) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Linked boundary conditions require a box for a simulation domain.')")
          write (report_unit, fmt="('You have grid_option = ',a)") trim(grid_option)
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       else
          write (report_unit, *) 
          write (report_unit, fmt="('Linked (twist and shift) boundary conditions will be used.')")
          write (report_unit, *) 
       end if
    case (boundary_option_self_periodic)
       write (report_unit, *) 
       write (report_unit, fmt="('Periodic boundary conditions will be used.')")
       write (report_unit, fmt="('(No twist and shift.)')")
       write (report_unit, *) 
    case default
       write (report_unit, *) 
       write (report_unit, fmt="('Outgoing boundary conditions will be used.')")
       write (report_unit, *) 
    end select

    if (.not. has_electrons) then
       select case (adiabatic_option_switch)
          case (adiabatic_option_default)
             write (report_unit, *) 
             write (report_unit, fmt="('The adiabatic electron response is of the form:')")
             write (report_unit, *) 
             write (report_unit, fmt="('             ne = Phi')")
             write (report_unit, *) 
             write (report_unit, fmt="('This is appropriate for an ETG simulation,')") 
             write (report_unit, fmt="('where the role of ions and electrons in GS2 is switched.')")
             write (report_unit, *) 
    
          case (adiabatic_option_fieldlineavg)
             write (report_unit, *) 
             write (report_unit, fmt="('The adiabatic electron response is of the form:')")
             write (report_unit, *) 
             write (report_unit, fmt="('             ne = Phi - <Phi>')")
             write (report_unit, *) 
             write (report_unit, fmt="('The angle brackets denote a proper field line average.')") 
             write (report_unit, fmt="('This is appropriate for an ITG simulation.')") 
             write (report_unit, *) 
             
          case (adiabatic_option_yavg)
             write (report_unit, *) 
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('The adiabatic electron response is of the form:')")
             write (report_unit, *) 
             write (report_unit, fmt="('             ne = Phi - <Phi>_y')")
             write (report_unit, *) 
             write (report_unit, fmt="('The angle brackets denote an average over y only.')") 
             write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
             write (report_unit, fmt="('Perhaps you want field-line-average-term for adiabatic_option.')") 
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 

          case (adiabatic_option_noJ)
             write (report_unit, *) 
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('The adiabatic electron response is of the form:')")
             write (report_unit, *) 
             write (report_unit, fmt="('             ne = Phi - <Phi>')")
             write (report_unit, *) 
             write (report_unit, fmt="('The angle brackets denote a field-line average, but without the proper Jacobian factors.')") 
             write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
             write (report_unit, fmt="('Perhaps you want field-line-average-term for adiabatic_option.')") 
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 
          end select
       end if

       if (poisfac /= 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('Quasineutrality is not enforced.  The ratio (lambda_Debye/rho)**2 = ',e10.4)") poisfac
          write (report_unit, *) 
       end if
          
       if (mult_imp .and. nonlinear_mode_switch == nonlinear_mode_on) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('For nonlinear runs, all species must use the same values of fexpr and bakdif')")
          write (report_unit, fmt="('in the dist_fn_species_knobs_x namelists.')")
          write (report_unit, fmt="('THIS IS AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if

       if (test_df) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Test = T in the dist_fn_knobs namelist will stop the run before ')")
          write (report_unit, fmt="('any significant calculation is done, and will result in several ')")
          write (report_unit, fmt="('variables that determine array sizes to be written to the screen.')")
          write (report_unit, fmt="('THIS MAY BE AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if

       if (def_parity .and. nonlinear_mode_switch == nonlinear_mode_on) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Choosing a definite parity for a nonlinear run has never been tested.')")
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if

       if (def_parity) then
          if (even) then
             write (report_unit, fmt="('Only eigenmodes of even parity will be included.')")
          else
             write (report_unit, fmt="('Only eigenmodes of odd parity will be included.')")
          end if
       end if

       if (D_kill > 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('D_kill > 0 probably does not work correctly.')")
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if

    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 
    write (report_unit, fmt="('The ExB parameter is ',f7.4)") g_exb

    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 

    select case (source_option_switch)
       
    case (source_option_full)
       write (report_unit, *) 
       write (report_unit, fmt="('The standard GK equation will be solved.')")
       write (report_unit, *) 
       
    case(source_option_phiext_full)
       write (report_unit, *) 
       write (report_unit, fmt="('The standard GK equation will be solved,')")
       write (report_unit, fmt="('with an additional source proportional to Phi*F_0')")
       write (report_unit, fmt="('Together with phi_ext = -1., this is the usual way to &
             & calculate the Rosenbluth-Hinton response.')")
       write (report_unit, *) 

    case(source_option_test2_full)
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The standard GK equation will be solved,')")
       write (report_unit, fmt="('with additional developmental sources, determined by ')")
       write (report_unit, fmt="('source_option=test2_full in the source_knobs namelist.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    case(source_option_convect_full)
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The standard GK equation will be solved,')")
       write (report_unit, fmt="('with additional developmental sources, determined by ')")
       write (report_unit, fmt="('source_option=convect_full in the source_knobs namelist.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    case (source_option_zero)
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The GK distribution function will be advanced non-self-consistently.')")
       write (report_unit, fmt="('source_option=zero in the source_knobs namelist.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    case (source_option_sine)
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The GK distribution function will be advanced non-self-consistently.')")
       write (report_unit, fmt="('source_option=sine in the source_knobs namelist.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    case (source_option_test1)
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The GK distribution function will be advanced non-self-consistently.')")
       write (report_unit, fmt="('source_option=test1 in the source_knobs namelist.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    case (source_option_cosine)
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The GK distribution function will be advanced non-self-consistently.')")
       write (report_unit, fmt="('source_option=cosine in the source_knobs namelist.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 

    end select

    if (a_ext /= 0.0) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The variable a_ext in the source_knobs namelist does nothing.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (akx_star /= 0.0) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The variable akx_star in the source_knobs namelist does nothing.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (aky_star /= 0.0) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The variable aky_star in the source_knobs namelist does nothing.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 

    if (no_driver) then
       ! nothing
    else
       if (amplitude == 0.0) then 
          write (report_unit, *) 
          write (report_unit, fmt="('No Langevin antenna included.')")
          write (report_unit, *) 
       else
          write (report_unit, *) 
          write (report_unit, fmt="('A Langevin antenna is included, with characteristics:')")
          write (report_unit, *) 
          write (report_unit, fmt="('Frequency:  (',e10.4,', ',e10.4,')')") w_antenna
          write (report_unit, fmt="('Number of independent k values: ',i3)") nk_stir
          if (write_antenna) then
             write (report_unit, *) 
             write (report_unit, fmt="('Antenna data will be written to ',a)") trim(run_name)//'.antenna'
             write (report_unit, *) 
          end if
          write (report_unit, fmt="('k values:')")
          do i=1,nk_stir
             if (trav(i)) then
                write (report_unit, fmt="('Travelling wave:')")
                write (report_unit, fmt="('   kx = ',i2,'    ky = ',i2,'    kz = ',i2)") &
                     & kx_stir(i), ky_stir(i), kz_stir(i)
             else
                write (report_unit, fmt="('Standing wave:')")
                write (report_unit, fmt="('   kx = ',i2,'    ky = ',i2,'    kz = ',i2)") &
                     & kx_stir(i), ky_stir(i), kz_stir(i)
             end if
          end do
       end if
    end if

    select case (hyper_option_switch)
    case (hyper_option_none)
       if (D_hyperres > 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses no hyperresistivity.  &
	       &D_hyperres ignored.')") trim(hyper_option)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          D_hyperres = -10.
       end if
       if (D_hypervisc > 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses no hyperviscosity.  &
               &D_hypervisc ignored.')") trim(hyper_option)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          D_hypervisc = -10.
       endif

    case (hyper_option_visc)
       hyper_on = .true.
       if (D_hyperres > 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses no hyperresistivity.  &
            &D_hyperres ignored.')") trim(hyper_option)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          D_hyperres = -10.
       end if
       if (D_hypervisc < 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses hyperviscosity but &
               &D_hypervisc < 0.')") trim(hyper_option)
          write (report_unit, fmt="('No hyperviscosity used.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          hyper_on = .false.
       endif

    case (hyper_option_res)
       hyper_on = .true.
       if (D_hyperres < 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses hyperresistivity but D_hyperres < 0.')") trim(hyper_option)
          write(report_unit, fmt="('No hyperresistivity used.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          hyper_on = .false.
       end if
       if (D_hypervisc > 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses no hyperviscosity.  D_hypervisc ignored.')") trim(hyper_option)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          D_hypervisc = -10.
       endif

    case (hyper_option_both)
       hyper_on = .true.
       if (D_hyperres < 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses hyperresistivity but D_hyperres < 0.')") trim(hyper_option)
          write (report_unit, fmt="('No hyperresistivity used.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if
       if (D_hypervisc < 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses hyperviscosity but D_hypervisc < 0.')") trim(hyper_option)
          write (report_unit, fmt="('No hyperviscosity used.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       endif
       if (D_hypervisc < 0. .and. D_hyperres < 0.) hyper_on = .false.

    end select

    if (hyper_on) then
       write (report_unit, *) 
       write (report_unit, fmt="('------------------------------------------------------------')")
       write (report_unit, *) 

       select case (hyper_option_switch)

          case (hyper_option_visc)

          write (report_unit, *) 
          write (report_unit, fmt="('Hyperviscosity included without hyperresistivity.')")
          if (const_amp) then
             write (report_unit, fmt="('Damping rate is ',e10.4,' at highest k_perp.')") D_hypervisc
          else
             write (report_unit, fmt="('The damping coefficent is ',e10.4)") D_hypervisc
             write (report_unit, fmt="('The damping rate is proportional to the RMS amplitude of the turbulence.')")
          end if
          if (isotropic_shear) then
             write (report_unit, fmt="('The hyperviscosity is isotropic in the perpendicular plane.')")
             write (report_unit, fmt="('This is appropriate for MHD-like calculations.')")
          else
             write (report_unit, fmt="('The hyperviscosity is anisotropic in the perpendicular plane.')")
             write (report_unit, fmt="('This is appropriate for drift-type calculations.')")
             write (report_unit, fmt="('omega_osc = ',e10.4)") omega_osc
          end if
          
       case (hyper_option_res)

          write (report_unit, *) 
          write (report_unit, fmt="('Hyperresistivity included without hyperviscosity.')")
          if (const_amp) then
             write (report_unit, fmt="('Damping rate is ',e10.4,' at highest k_perp.')") D_hyperres
          else
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('const_amp = .false. is not implemented for hyperresistivity.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 
          end if
          if (isotropic_shear) then
             write (report_unit, fmt="('The hyperresistivity is isotropic in the perpendicular plane.')")
             write (report_unit, fmt="('This is appropriate for MHD-like calculations.')")
          else
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('isotropic_shear = .false. is not implemented for hyperresistivity.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 
          end if

       case (hyper_option_both)
          
          write (report_unit, *) 
          write (report_unit, fmt="('Hyperresistivity and hyperviscosity included.')")
          if (const_amp) then
             write (report_unit, fmt="('Damping rate is ',e10.4,' at highest k_perp.')") D_hyperres
          else
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('const_amp = .false. is not implemented for hyperresistivity.')")
             write (report_unit, fmt="('THIS IS AN ERROR.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 
          end if
          if (isotropic_shear) then
             write (report_unit, fmt="('The damping is isotropic in the perpendicular plane.')")
             write (report_unit, fmt="('This is appropriate for MHD-like calculations.')")
          else
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('isotropic_shear = .false. is not implemented for hyperresistivity.')")
             write (report_unit, fmt="('THIS IS AN ERROR.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 
          end if
       end select
    end if

    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 
    write (report_unit, fmt="&
         &('Integrals over energy are broken into two regions, 0:Ecut and Ecut:infinity.')")
    write (report_unit, fmt="('Ecut = ',f8.4)") ecut
    write (report_unit, fmt="('There are ',i3,' total energy grid points.')") negrid
    if (test_le) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('Test = T in the le_grids_knobs namelist will stop the run before ')")
       write (report_unit, fmt="('any significant calculation is done.')")
       write (report_unit, fmt="('The lambda and energy grids will be written to the screen, and the run will stop.')")
       write (report_unit, fmt="('THIS MAY BE AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if
    
    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 
    if (fphi /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('fphi in the knobs namelist = ',e10.4)") fphi
       write (report_unit, fmt="('fphi is a scale factor of all instances of Phi (the electrostatic potential).')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (fapar == 0.) then
       write (report_unit, fmt="('A_parallel will not be included in the calculation.')")
    end if
    if (fapar == 1.) then
       write (report_unit, fmt="('A_parallel will be included in the calculation.')")
    end if
    if (fapar /= 0. .and. fapar /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('fapar in the knobs namelist = ',e10.4)") fapar
       write (report_unit, fmt="('fapar is a scale factor of all instances of A_parallel (the parallel vector potential).')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (fbpar == 0.) then
       write (report_unit, fmt="('B_parallel will not be included in the calculation.')")
    end if
    if (fbpar == 1.) then
       write (report_unit, fmt="('B_parallel will be included in the calculation.')")
    end if
    if (fbpar /= 0. .and. fbpar /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('fbpar in the knobs namelist = ',e10.4)") fbpar
       write (report_unit, fmt="('fbpar is a scale factor of all instances of B_parallel &
           & (the perturbed parallel magnetic field).')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (eqzip) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('eqzip = T in the knobs namelist.')")
       write (report_unit, fmt="('This freezes some modes in time for a secondary stability analysis.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
       if (secondary) write (report_unit, fmt="('Mode with kx = 0, ky = ky_min fixed in time')")
       if (tertiary)  write (report_unit, fmt="('Mode with ky = 0, kx = kx_min fixed in time')")
    end if

! diagnostic controls:

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
          write (report_unit, fmt="('write_eigenfunc = T:       Normalized Phi(theta) written to ',a)") &
	& 	trim(run_name)//'.eigenfunc'
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

    if (write_intcheck) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('write_intcheck = T:              Turns on obscure diagnostics.')")
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
    

    call close_output_file (report_unit)

  end subroutine report

  subroutine nprocs (nmesh)

    implicit none
    real :: fac
    integer, intent (in) :: nmesh
    integer :: nefacs, nlfacs, nkyfacs, nkxfacs, nspfacs
    integer, dimension(:,:), allocatable :: facs
    integer :: npe

    if (nonlinear_mode_switch == nonlinear_mode_on) then
       select case (layout)
       case ('lexys')

          write (report_unit, *) 
          write (report_unit, fmt="('Recommended numbers of processors, time on T3E')") 
          allocate (facs(max(nspec,naky,ntheta0)/2+1,3))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (naky, nkyfacs, facs(:,2))
          call factors (ntheta0, nkxfacs, facs(:,3))
          do i=1,nspfacs-1
             npe = facs(i,1)
             if (nmesh/npe > ncut) write &
                  & (report_unit, fmt="('  npe = ',i4,'    time = ',f8.2,'  seconds/time step')") &
                  & npe, 40.*nmesh/1.e7/npe**0.85
          end do
          do i=1,nkyfacs-1
             npe = facs(i,2)*nspec
             if (nmesh/npe > ncut) write &
                  & (report_unit, fmt="('  npe = ',i4,'    time = ',f8.2,'  seconds/time step')") &
                  & npe, 40.*nmesh/1.e7/npe**0.85
          end do
          do i=1,nkxfacs
             npe = facs(i,3)*naky*nspec
             if (nmesh/npe > ncut) write &
                  & (report_unit, fmt="('  npe = ',i4,'    time = ',f8.2,'  seconds/time step')") &
                  & npe, 40.*nmesh/1.e7/npe**0.85
          end do
          deallocate (facs)

       case ('lxyes')

          write (report_unit, *) 
          write (report_unit, fmt="('Recommended numbers of processors, time on SP2')") 
          allocate (facs(max(nspec,negrid,naky,ntheta0)/2+1,4))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (negrid, nefacs, facs(:,2))
          call factors (naky, nkyfacs, facs(:,3))
          call factors (ntheta0, nkxfacs, facs(:,4))
!          faclin = 3.5*(real(nmesh))**1.1/1.e7/5.
          fac = 3.5*(real(nmesh))**1.1/1.e7
          do i=1,nspfacs-1
             npe = facs(i,1)
             if (nmesh/npe > ncut) write &
                  & (report_unit, fmt="('  npe = ',i4,'    time = ',f8.2,'  seconds/time step')") &
                  & npe, fac/npe**0.95
          end do
          do i=1,nefacs-1
             npe = facs(i,2)*nspec
             if (nmesh/npe > ncut) write &
                  & (report_unit, fmt="('  npe = ',i4,'    time = ',f8.2,'  seconds/time step')") &
                  & npe, fac/npe**0.95
          end do
          do i=1,nkyfacs-1
             npe = facs(i,3)*negrid*nspec
             if (nmesh/npe > ncut) write &
                  & (report_unit, fmt="('  npe = ',i4,'    time = ',f8.2,'  seconds/time step')") &
                  & npe, fac/npe**0.95
          end do
          do i=1,nkxfacs
             npe = facs(i,4)*naky*negrid*nspec
             if (nmesh/npe > ncut) write &
                  & (report_unit, fmt="('  npe = ',i4,'    time = ',f8.2,'  seconds/time step')") &
                  & npe, fac/npe**0.95
          end do
          deallocate (facs)

       case ('yxels')

          write (report_unit, *) 
          write (report_unit, fmt="('Recommended numbers of processors:')") 

          allocate (facs(max(nspec,negrid,nlambda)/2+1,3))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (nlambda, nlfacs, facs(:,2))
          call factors (negrid, nefacs, facs(:,3))
          fac = 3.5*(real(nmesh))**1.1/1.e7
          do i=1,nspfacs-1
             npe = facs(i,1)
             if (nmesh/npe > ncut) write &
                  & (report_unit, fmt="('  npe = ',i4,'    time = ',f8.2,'  seconds/time step')") &
                  & npe, fac/npe**0.95
          end do
          do i=1,nlfacs-1
             npe = facs(i,2)*nspec
             if (nmesh/npe > ncut) write &
                  & (report_unit, fmt="('  npe = ',i4,'    time = ',f8.2,'  seconds/time step')") &
                  & npe, fac/npe**0.95
          end do
          do i=1,nefacs
             npe = facs(i,3)*nlambda*nspec
             if (nmesh/npe > ncut) write &
                  & (report_unit, fmt="('  npe = ',i4,'    time = ',f8.2,'  seconds/time step')") &
                  & npe, fac/npe**0.95
          end do
          deallocate (facs)

       case ('yxles')

          write (report_unit, *) 
          write (report_unit, fmt="('Recommended numbers of processors:')") 

          allocate (facs(max(nspec,negrid,nlambda)/2+1,3))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (negrid, nefacs, facs(:,2))
          call factors (nlambda, nlfacs, facs(:,3))
          fac = 3.5*(real(nmesh))**1.1/1.e7
          do i=1,nspfacs-1
             npe = facs(i,1)
             if (nmesh/npe > ncut) write &
                  & (report_unit, fmt="('  npe = ',i4,'    time = ',f8.2,'  seconds/time step')") &
                  & npe, fac/npe**0.95
          end do
          do i=1,nefacs-1
             npe = facs(i,2)*nspec
             if (nmesh/npe > ncut) write &
                  & (report_unit, fmt="('  npe = ',i4,'    time = ',f8.2,'  seconds/time step')") &
                  & npe, fac/npe**0.95
          end do
          do i=1,nlfacs
             npe = facs(i,3)*negrid*nspec
             if (nmesh/npe > ncut) write &
                  & (report_unit, fmt="('  npe = ',i4,'    time = ',f8.2,'  seconds/time step')") &
                  & npe, fac/npe**0.95
          end do
          deallocate (facs)

       end select
    else
       write (report_unit, *) 
       write (report_unit, fmt="('Recommended numbers of processors:')")
       select case (gridopt_switch)

       case (gridopt_single)

          select case (layout)
          case ('lexys', 'yxles', 'lxyes')
             allocate (facs(max(negrid,nspec,nlambda)/2+1,3))
             call factors (nspec, nspfacs, facs(:,1))
             call factors (negrid, nefacs, facs(:,2))
             call factors (nlambda, nlfacs, facs(:,3))
             do i=1,nspfacs-1
                npe = facs(i,1)
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nefacs-1
                npe = facs(i,2)*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nlfacs
                npe = facs(i,3)*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             deallocate (facs)

          case ('yxels')
             allocate (facs(max(negrid,nspec,nlambda)/2+1,3))
             call factors (nspec, nspfacs, facs(:,1))
             call factors (nlambda, nlfacs, facs(:,2))
             call factors (negrid, nefacs, facs(:,3))
             do i=1,nspfacs-1
                npe = facs(i,1)
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nlfacs-1
                npe = facs(i,2)*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nefacs
                npe = facs(i,3)*nlambda*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             deallocate (facs)
          end select

       case (gridopt_range,gridopt_specified,gridopt_box)

          select case (layout)
          case ('lexys')

             allocate (facs(max(nspec,naky,ntheta0,negrid,nlambda)/2+1,5))
             call factors (nspec, nspfacs, facs(:,1))
             call factors (naky, nkyfacs, facs(:,2))
             call factors (ntheta0, nkxfacs, facs(:,3))
             call factors (negrid, nefacs, facs(:,4))
             call factors (nlambda, nlfacs, facs(:,5))
             do i=1,nspfacs-1
                npe = facs(i,1)
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nkyfacs-1
                npe = facs(i,2)*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nkxfacs-1
                npe = facs(i,3)*naky*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nefacs-1
                npe = facs(i,4)*ntheta0*naky*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nlfacs
                npe = facs(i,5)*negrid*ntheta0*naky*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             deallocate (facs)

          case ('lxyes')

             write (report_unit, *) 
             allocate (facs(max(nspec,negrid,naky,ntheta0)/2+1,4))
             call factors (nspec, nspfacs, facs(:,1))
             call factors (negrid, nefacs, facs(:,2))
             call factors (naky, nkyfacs, facs(:,3))
             call factors (ntheta0, nkxfacs, facs(:,4))
!             fac = 3.5*(real(nmesh))**1.1/1.e7/5.  ! okay for large runs, not small
             do i=1,nspfacs-1
                npe = facs(i,1)
!                if (nmesh/npe > ncut) write &
!                     & (report_unit, fmt="('  npe = ',i4,'    time = ',f8.2,'  seconds/time step')") &
!                     & npe, fac/npe**0.95
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nefacs-1
                npe = facs(i,2)*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nkyfacs-1
                npe = facs(i,3)*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nkxfacs
                npe = facs(i,4)*naky*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             deallocate (facs)

          case ('yxles')

             allocate (facs(max(nspec,naky,ntheta0,negrid,nlambda)/2+1,5))
             call factors (nspec, nspfacs, facs(:,1))
             call factors (negrid, nefacs, facs(:,2))
             call factors (nlambda, nlfacs, facs(:,3))
             call factors (ntheta0, nkxfacs, facs(:,4))
             call factors (naky, nkyfacs, facs(:,5))
             do i=1,nspfacs-1
                npe = facs(i,1)
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nefacs-1
                npe = facs(i,2)*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nlfacs-1
                npe = facs(i,3)*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nkxfacs-1
                npe = facs(i,4)*nlambda*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nkyfacs
                npe = facs(i,5)*ntheta0*nlambda*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             deallocate (facs)

          case ('yxels')

             allocate (facs(max(nspec,naky,ntheta0,negrid,nlambda)/2+1,5))
             call factors (nspec, nspfacs, facs(:,1))
             call factors (nlambda, nlfacs, facs(:,2))
             call factors (negrid, nefacs, facs(:,3))
             call factors (ntheta0, nkxfacs, facs(:,4))
             call factors (naky, nkyfacs, facs(:,5))
             do i=1,nspfacs-1
                npe = facs(i,1)
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nlfacs-1
                npe = facs(i,2)*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nefacs-1
                npe = facs(i,3)*nlambda*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nkxfacs-1
                npe = facs(i,4)*negrid*nlambda*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nkyfacs
                npe = facs(i,5)*ntheta0*negrid*nlambda*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             deallocate (facs)

          end select
       end select
    end if

  end subroutine nprocs

  subroutine tell (a, b, c)
    
    integer :: j
    character (*), intent (in) :: a
    character (*), intent (in), optional :: b, c
    character (500) :: hack

    j = len_trim(a)

    if (present(b)) then
       j = max(j, len_trim(b))
       if (present(c)) j = max(j, len_trim(c))
    end if

    write (6, *)
    write (6, *)
    hack = repeat('*', j+2)
    write (6, 10) '  '//trim(hack)
    write (6, 10) '   '//trim(a)
    write (interactive_record, 10) trim(a)
    if (present(b)) then
       write (6, 10) '   '//trim(b)
       write (interactive_record, 10) trim(b)
       if (present(c)) then
          write (6, 10) '   '//trim(c)
          write (interactive_record, 10) trim(c)
       end if
    end if
    write (6, 10) '  '//trim(hack)
    write (6, *)
    write (6, *)

10 format (a)

  end subroutine tell

  subroutine text (a)
    
    character (*), intent (in), optional :: a

    if (present(a)) then
       write (6, 10) trim(a)
    else
       write (*, *)
    end if

10 format (a)

  end subroutine text

  subroutine run_number (sel, nbeta)
    
    integer, intent (out) :: sel
    integer, intent (in) :: nbeta
    
100 continue
    do 
       call text
       call text ('Lowest run number for this scan (0-999):')
       read (interactive_input, *, err = 100) sel
       if (sel + nbeta - 1 > 999) then
          write (6,*) 'The lowest run number for this scan must be less than ',1000-nbeta
          call try_again
       else if (sel < 0) then
          write (6,*) 'The lowest run number for this scan must be greater than zero.'
          call try_again
       else
          return
       end if
    end do
    
  end subroutine run_number
  
  subroutine try_again 
    
    write (6, *) '*****************'
    write (6, *) 'Please try again.'
    write (6, *) '*****************'

  end subroutine try_again

  subroutine get_choice (in, out)

    integer, intent (in) :: in
    integer, intent (out) :: out

    read (interactive_input, *, err=999) out
    if (out <= 0) out = 0
    if (out > in) out = 0

    return
999 out = 0
    
  end subroutine get_choice
  
  subroutine beta_range_low (x, a, x0)
    
    real, intent (out) :: x
    character (*), intent (in) :: a
    real, intent (in) :: x0
    
    do
100    continue
       call text
       if (trim(a) == 'le') then
          call text ('Lower limit of beta for this scan (must be > 0):')
       else if(trim(a) == 'lt') then
          call text ('Lower limit of beta for this scan (must be >= 0):')
       end if
       read (interactive_input, *, err=100) x
       if (trim(a) == 'le') then
          if (x <= x0) then
             call try_again
          else
             return
          end if
       else if (trim(a) == 'lt') then
          if (x < x0) then
             call try_again
          else
             return
          end if
       end if
    end do
    
  end subroutine beta_range_low

  subroutine beta_range_high (x, a, x0)
    
    real, intent (out) :: x
    character (*), intent (in) :: a
    real, intent (in) :: x0
    
    do
100    continue
       call text
       call text ('Upper limit of beta for this scan (must be > lower limit):')
       read (interactive_input, *, err=100) x
       if (trim(a) == 'le') then
          if (x <= x0) then
             call try_again
          else
             return
          end if
       end if
    end do
    
  end subroutine beta_range_high

  subroutine beta_prime_range_low (x)

    real, intent (out) :: x
    
    do
100    continue
       call tell ('Note: You will be asked to enter:         - d beta/d rho', &
            'Since d beta /d rho <= 0, you should enter a number >= 0')
       
       call text ('Weakest beta gradient in scan (zero or positive number):')
       read (interactive_input, *, err=100) x
       x = -x
       if (x > 0.) then
          call try_again
       else
          return
       end if
    end do

  end subroutine beta_prime_range_low
  
  subroutine beta_prime_range_high (x, x0)
    
    real, intent (out) :: x
    real, intent (in) :: x0
    
    do
100    continue
       call text
       call text ('Strongest gradient in scan: (positive number)')
       read (interactive_input, *, err=100) x
       x = -x
       if (x >= x0) then
          call try_again
       else
          return
       end if
    end do

  end subroutine beta_prime_range_high

  subroutine num_runs (n)

    integer, intent (out) :: n
    
    do 
100    continue
       call text
       call text ('How many runs should be done? (n > 1)')
       read (interactive_input, *, err=100) n
       if (n < 2) then
          call try_again
       else
          return
       end if
    end do
  end subroutine num_runs

end program ingen

