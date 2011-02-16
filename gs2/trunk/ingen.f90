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
  integer :: in_file, i, ierr, unit, is, report_unit, iunit, ncut, npmax
  logical :: exist, scan, stdin

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
  logical :: conserve_number, conserve_momentum, hypercoll, const_v
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

! dist_fn:
  complex, dimension (:), allocatable :: fexp ! (nspec)
  real, dimension (:), allocatable :: bkdiff  ! (nspec)
  integer, dimension (:), allocatable :: bd_exp ! nspec
  real :: gridfac, apfac, driftknob, tpdriftknob, poisfac
  real :: kfilter, afilter, D_kill, noise, g_exb, g_exbfac, omprimfac, btor_slab
  real :: t0, omega0, gamma0, source0, thetas, k0, phi_ext, a_ext
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

! le_grids
  integer :: ngauss, negrid, nesuper, nesub
  real :: ecut, vcut, bouncefuzz
  logical :: trapped_particles = .true.
  logical :: advanced_egrid = .true.
  logical :: vgrid = .false.

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


  logical :: test_le
  logical :: layouts_write = .false.
  logical :: driver_write = .false.
  logical :: collisions_write = .false.
  logical :: init_g_write = .false.
  logical :: dist_fn_write = .false.
  logical :: source_write = .false.
  logical :: fields_write = .false.
  logical :: diagnostics_write = .false.
  logical :: reinit_write = .false.
  logical :: hyper_write = .false.
  logical :: le_write = .false.
  logical :: nonlinear_write = .false.
  logical :: run_parameters_write = .false.
  logical :: species_write = .false.
  logical :: theta_parameters_write = .false.
  logical :: theta_gridgen_write = .false.
  logical :: theta_salpha_write = .false.
  logical :: theta_eik_write = .false.
  logical :: theta_file_write = .false.
  logical :: theta_write = .false.

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
       conserve_number, conserve_momentum, hypercoll, cfac, heating

! dist_fn:
  namelist /dist_fn_knobs/ boundary_option, gridfac, apfac, &
       driftknob, tpdriftknob, g_exb, g_exbfac, omprimfac, btor_slab,&
       nperiod_guard, poisfac, adiabatic_option, cfac, &
       kfilter, afilter, mult_imp, test, def_parity, even, &
       save_n, save_u, save_Tpar, save_Tperp, D_kill, noise, heating_option
  
  namelist /source_knobs/ t0, omega0, gamma0, source0, &
       thetas, k0, phi_ext, source_option, a_ext, aky_star, akx_star

!  namelist /dist_fn_species_knobs/ fexpr, fexpi, bakdif, bd_exp

! fields: 
  namelist /fields_knobs/ field_option

! gs2_reinit: 
  namelist /reinit_knobs/ delt_adj, delt_minimum

! hyper: 
  namelist /hyper_knobs/ hyper_option, const_amp, include_kpar, &
       isotropic_shear, D_hyperres, D_hypervisc, omega_osc, D_hyper, gridnorm

! le_grids:

  namelist /le_grids_knobs/ ngauss, negrid, ecut, bouncefuzz, &
       nesuper, nesub, test, trapped_particles, advanced_egrid, &
       vgrid, vcut

! nonlinear_terms: 
  namelist /nonlinear_terms_knobs/ nonlinear_mode, flow_mode, cfl, &
       C_par, C_perp, p_x, p_y, p_z, zip

  namelist /ingen_knobs/ ncut, scan, stdin, npmax

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
    use species, only: spec, nspec, has_electron_species
    use geometry, only: beta_prime_input, bishop
    use run_parameters, only: beta, fapar, fbpar
    integer :: sel, nbeta, j, ise
    real :: beta_low, beta_high, dbeta, beta_save
    real :: fapar_save, fbpar_save, pri, pe, alpi, tpe_save, ptot, alp, dbdr
    real :: alt, aln, fac, beta_prime_save, bishop_save
    real, dimension (:), allocatable :: tp_save, fp_save
    character (500) :: tag1, tag2
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
             
             if (.not. has_electron_species(spec)) call tell ('You really should use electrons for electromagnetic runs.')

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
             
             if (.not. has_electron_species(spec)) call tell ('You really should use electrons for electromagnetic runs.')

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
             
             if (.not. has_electron_species(spec)) call tell ('You really should use electrons for electromagnetic runs.')

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
             fbpar_save = fbpar

             beta = beta_low + (j - sel)*dbeta
             if (beta == 0.) then 
                fapar = 0.
                fbpar = 0.
             else
                if (fapar == 0. .and. fbpar == 0.) then
                   fapar = 1.0 ;  fbpar = 1.0
                end if
             end if
             
             write (tag1, fmt='("Varying beta, all else fixed")')
             write (tag2, fmt='("beta = ",e16.10)') beta

             call write_namelists (j, tag1, tag2)

             fapar = fapar_save 
             fbpar = fbpar_save
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

                fapar_save = fapar ; fbpar_save = fbpar
                if (beta == 0.) then
                   fapar = 0.      ; fbpar = 0.
                else
                   if (fapar == 0. .and. fbpar == 0.) then
                      fapar = 1.0 ;  fbpar = 1.0
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
                fbpar = fbpar_save
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

  end subroutine interactive

  subroutine get_namelists
!CMR, 17/11/2009: use gs2_diagnostics module to pick up public variables
!                 and routines to reduces maintenance of ingen.
!                 Should replicate this for other modules in future.
!CMR, 2/2/2011:   Have extended this to include use of theta_grid module:
!                 Strategy is simply to add two types of routines to modules:
!                      wnml_xxxxx   to write the modules namelists
!                      check_xxxxx  to perform the ingen checking inside the module
!                 More object oriented strategy, easier maintenance of ingen.f90 
!                 which is gradually shrinking.!                  
!
    use run_parameters, only: init_run_parameters
    use init_g, only: init_init_g
    use species, only: init_species, nspec
    use gs2_diagnostics,only: gs2diag_read_parameters=>read_parameters
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids
    use theta_grid, only: init_theta_grid, nbset, shat_real => shat
    use theta_grid_params, only: shat
    use constants, only: pi
    logical :: list
!CMR
    logical:: debug=.false.
    call init_file_utils (list, name="template")
if (debug) write(6,*) 'get_namelists: called init_file_utils'
    ncut= 100000
    npmax=10000
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
    call init_init_g

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
    call init_kt_grids

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
    call init_run_parameters
    run_parameters_write = .true.


if (debug) write(6,*) 'get_namelists: init_species'
    call init_species

!CMR, 2/2/2011:  reduce much duplication by calling init_theta_grid
if (debug) write(6,*) 'get_namelists: call init_theta_grid'
    call init_theta_grid

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
    omprimfac=1.0
    g_exbfac = 1.0
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

! CMR 18/11/2009:  reduce duplication by calling gs2diag_read_parameters 
!
   call gs2diag_read_parameters(.false.)

if (debug) write(6,*) 'get_namelists: returning'

  end subroutine get_namelists

  subroutine write_namelists (jr, tag1, tag2)
    use gs2_diagnostics, only: wnml_gs2_diagnostics
    use run_parameters, only: wnml_run_parameters
    use species, only: wnml_species, nspec, spec, has_electron_species
    use theta_grid_params, only: wnml_theta_grid_params
    use theta_grid_gridgen, only: wnml_theta_grid_gridgen
    use theta_grid_salpha, only : wnml_theta_grid_salpha
    use theta_grid_eik, only : wnml_theta_grid_eik
    use theta_grid_file, only : wnml_theta_grid_file
    use init_g, only : wnml_init_g
    use kt_grids, only: wnml_kt
    integer, intent (in), optional :: jr
    character (*), intent (in), optional :: tag1, tag2

    character (100) :: line
    integer :: h, t, u
    integer :: i
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

    if (init_g_write) call wnml_init_g(unit)

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

       if (.not. has_electron_species(spec)) then
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

    if (diagnostics_write) call wnml_gs2_diagnostics(unit)


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

    call wnml_kt(unit)
    
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

    if (nonlinear_write .and. nonlinear_mode_switch == nonlinear_mode_on) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "nonlinear_terms_knobs"
       write (unit, fmt="(' nonlinear_mode = ',a)") '"on"'
       write (unit, fmt="(' cfl = ',e16.10)") cfl
       if (zip) write (unit, fmt="(' zip = ',L1)") zip
       write (unit, fmt="(' /')")
    end if

    if (run_parameters_write) call wnml_run_parameters(unit)

    if (species_write) call wnml_species(unit)
    do i=1,nspec
       write (unit, *)
       write (line, *) i
       write (unit, fmt="(' &',a)") &
            & trim("dist_fn_species_knobs_"//trim(adjustl(line)))
       write (unit, fmt="(' fexpr = ',e13.6)") real(fexp(i))
       write (unit, fmt="(' bakdif = ',e13.6)") bkdiff(i)
       write (unit, fmt="(' bd_exp = ',i6,'  /')") bd_exp(i)
    end do

    if (theta_write) call wnml_theta_grid(unit)

    if (theta_parameters_write) call wnml_theta_grid_params(unit)

    if (theta_gridgen_write) call wnml_theta_grid_gridgen(unit)

    if (theta_salpha_write) call wnml_theta_grid_salpha(unit)

    if (theta_eik_write) call wnml_theta_grid_eik(unit)

    if (theta_file_write) call wnml_theta_grid_file(unit)

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
     use run_parameters, only: check_run_parameters
     use run_parameters, only: beta, tite, delt, fapar, fbpar, fphi, margin, nstep, wstar_units
     use gs2_diagnostics, only: check_gs2_diagnostics
     use gs2_diagnostics, only: dump_fields_periodically, save_for_restart, nsave, make_movie, nmovie, exit_when_converged, nwrite, omegatol
      use species, only: check_species, spec, nspec, has_electron_species
      use init_g, only: check_init_g
      use theta_grid, only: check_theta_grid
     use theta_grid, only: gb_to_cv, nbset, ntgrid_real => ntgrid
     use theta_grid_params, only: nperiod, ntheta, eps, epsl, rmaj, r_geo
     use theta_grid_params, only: pk, qinp, rhoc, shift, shat
     use theta_grid_params, only: akappa, akappri, tri, tripri
     use kt_grids, only: check_kt_grids, grid_option, gridopt_switch
     use kt_grids, only: gridopt_box, naky, ntheta0, nx, ny

     implicit none
     real :: alne, dbetadrho_spec
     real :: kxfac, drhodpsi
     character (20) :: datestamp, timestamp, zone
     character (200) :: line
     logical :: coll_on = .false., le_ok = .true.
     integer :: ntgrid, j, nmesh
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

     call check_species(report_unit,beta,tite,alne,dbetadrho_spec)

     write (report_unit, *) 
     write (report_unit, fmt="('------------------------------------------------------------')")
     write (report_unit, *) 
     call check_theta_grid(report_unit,alne,dbetadrho_spec)

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

    call check_init_g(report_unit)
       
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

    call check_kt_grids(report_unit)

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

    if (.not. has_electron_species(spec)) then
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
    write (report_unit, fmt="('Perp shear terms will be multiplied by factor',f7.4)") g_exbfac
    write (report_unit, fmt="('Parallel shear term will be multiplied by factor',f7.4)") omprimfac

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
    call check_run_parameters(report_unit)

! diagnostic controls:
    call check_gs2_diagnostics(report_unit)

    call close_output_file (report_unit)

  end subroutine report

  subroutine nprocs (nmesh)

    use species, only : nspec
    use kt_grids, only: gridopt_switch, gridopt_single, gridopt_range, gridopt_specified, gridopt_box, gridopt_xbox
    use kt_grids, only: naky, ntheta0
    implicit none
    real :: fac
    integer, intent (in) :: nmesh
    integer :: nefacs, nlfacs, nkyfacs, nkxfacs, nspfacs
    integer, dimension(:,:), allocatable :: facs
    integer :: npe
    real :: time

    write (report_unit, fmt="('Layout = ',a5)") layout 
    if (nonlinear_mode_switch == nonlinear_mode_on) then
       select case (layout)
       case ('lexys')

          write (report_unit, *) 
          write (report_unit, fmt="('Recommended numbers of processors, time on T3E')") 
          allocate (facs(max(nspec,naky,ntheta0)/2+1,5))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (naky, nkyfacs, facs(:,2))
          call factors (ntheta0, nkxfacs, facs(:,3))
          call factors (negrid, nefacs, facs(:,4))
          call factors (nlambda, nlfacs, facs(:,5))
          do i=1,nspfacs
             npe = facs(i,1)
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= 40.*nmesh/1.e7/npe**0.85
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'s'
          end do
          do i=2,nkyfacs
             npe = facs(i,2)*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= 40.*nmesh/1.e7/npe**0.85
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'y'
          end do
          do i=2,nkxfacs
             npe = facs(i,3)*naky*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= 40.*nmesh/1.e7/npe**0.85
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step  (',a,')')") npe, time,'x'
          end do
          do i=2,nefacs
             npe = facs(i,4)*ntheta0*naky*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= 40.*nmesh/1.e7/npe**0.85
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step  (',a,')')") npe, time,'e'
          end do
          do i=2,nlfacs
             npe = facs(i,5)*negrid*ntheta0*naky*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= 40.*nmesh/1.e7/npe**0.85
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step  (',a,')')") npe, time,'l'
          end do
          deallocate (facs)

       case ('lxyes')

          write (report_unit, *) 
          write (report_unit, fmt="('Recommended numbers of processors, time on SP2')") 
          allocate (facs(max(nspec,negrid,naky,ntheta0)/2+1,5))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (negrid, nefacs, facs(:,2))
          call factors (naky, nkyfacs, facs(:,3))
          call factors (ntheta0, nkxfacs, facs(:,4))
          call factors (nlambda, nlfacs, facs(:,5))
!          faclin = 3.5*(real(nmesh))**1.1/1.e7/5.
          fac = 3.5*(real(nmesh))**1.1/1.e7
          do i=1,nspfacs
             npe = facs(i,1)
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'s'
          end do
          do i=2,nefacs
             npe = facs(i,2)*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'e'
          end do
          do i=2,nkyfacs
             npe = facs(i,3)*negrid*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'y'
          end do
          do i=2,nkxfacs
             npe = facs(i,4)*naky*negrid*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'x'
          end do
          do i=2,nlfacs
             npe = facs(i,5)*naky*ntheta0*negrid*nspec 
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'l'
          end do
          deallocate (facs)

       case ('yxels')

          write (report_unit, *) 
          write (report_unit, fmt="('Recommended numbers of processors:')") 

          allocate (facs(max(nspec,negrid,nlambda)/2+1,5))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (nlambda, nlfacs, facs(:,2))
          call factors (negrid, nefacs, facs(:,3))
          call factors (ntheta0, nkxfacs, facs(:,4))
          call factors (naky, nkyfacs, facs(:,5))
          fac = 3.5*(real(nmesh))**1.1/1.e7
          do i=1,nspfacs
             npe = facs(i,1)
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'s'
          end do
          do i=2,nlfacs
             npe = facs(i,2)*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'l'
          end do
          do i=2,nefacs
             npe = facs(i,3)*nlambda*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'e'
          end do
          do i=2,nkxfacs
             npe = facs(i,4)*negrid*nlambda*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'x'
          end do
          do i=2,nkyfacs
             npe = facs(i,5)*nkxfacs*negrid*nlambda*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'y'
          end do
          deallocate (facs)

       case ('yxles')

          write (report_unit, *) 
          write (report_unit, fmt="('Recommended numbers of processors:')") 

          allocate (facs(max(nspec,negrid,nlambda)/2+1,5))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (negrid, nefacs, facs(:,2))
          call factors (nlambda, nlfacs, facs(:,3))
          call factors (ntheta0, nkxfacs, facs(:,4))
          call factors (naky, nkyfacs, facs(:,5))
          fac = 3.5*(real(nmesh))**1.1/1.e7
          do i=1,nspfacs
             npe = facs(i,1)
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'s'
          end do
          do i=2,nefacs
             npe = facs(i,2)*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'e'
          end do
          do i=2,nlfacs
             npe = facs(i,3)*negrid*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'l'
          end do
          do i=2,nkxfacs
             npe = facs(i,4)*nlambda*negrid*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'x'
          end do
          do i=2,nkyfacs
             npe = facs(i,5)*nkxfacs*nlambda*negrid*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'y'
          end do
          deallocate (facs)

       case ('xyles')
!CMR, 11/11/2009: add processor recommendations for xyles layout
!            NB added recommendations that also parallelise in y and x, 
!                          which may be unwise!
          write (report_unit, *) 
          write (report_unit, fmt="('Recommended numbers of processors:')") 

          allocate (facs(max(nspec,negrid,nlambda,naky,ntheta0)/2+1,5))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (negrid, nefacs, facs(:,2))
          call factors (nlambda, nlfacs, facs(:,3))
          call factors (naky, nkyfacs, facs(:,4))
          call factors (ntheta0, nkxfacs, facs(:,5))
          fac = 3.5*(real(nmesh))**1.1/1.e7
          do i=1,nspfacs
             npe = facs(i,1)
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'s'
          end do
          do i=2,nefacs
             npe = facs(i,2)*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'e'
          end do
          do i=2,nlfacs
             npe = facs(i,3)*negrid*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'l'
          end do
          do i=2,nkyfacs
             npe = facs(i,4)*nlambda*negrid*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'y'
          end do
          do i=2,nkxfacs
             npe = facs(i,5)*naky*nlambda*negrid*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'x'
          end do
          deallocate (facs)

       end select
    else
       write (report_unit, *) 
       write (report_unit, fmt="('Recommended numbers of processors:')")
       select case (gridopt_switch)

       case (gridopt_single)

          select case (layout)
          case ('lexys', 'yxles', 'lxyes', 'xyles')
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

          case ('xyles')

             allocate (facs(max(nspec,naky,ntheta0,negrid,nlambda)/2+1,5))
             call factors (nspec, nspfacs, facs(:,1))
             call factors (negrid, nefacs, facs(:,2))
             call factors (nlambda, nlfacs, facs(:,3))
             call factors (naky, nkyfacs, facs(:,4))
             call factors (ntheta0, nkxfacs, facs(:,5))
             do i=1,nspfacs-1
                npe = facs(i,1)
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i8)") npe
             end do
             do i=1,nefacs-1
                npe = facs(i,2)*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i8)") npe
             end do
             do i=1,nlfacs-1
                npe = facs(i,3)*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i8)") npe
             end do
             do i=1,nkyfacs-1
                npe = facs(i,4)*nlambda*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i8)") npe
             end do
             do i=1,nkxfacs
                npe = facs(i,5)*ntheta0*nlambda*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i8)") npe
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
