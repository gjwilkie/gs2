! Notes from BD, 7/2011:
!
! Need to extend the verr tools to include delta B_parallel integrals
! There are new factors of 1/B here and there which I do not understand.  
!
!>The principal function of this module is to evolve the distribution function, 
!!that is, to advance the discrete gyrokinetic equation. 
!!This involves calculating the source and dealing with the complexities 
!!of the parallel boundary conditions. 
!!In addition it contains a routine for implementing perpendicular
!!velocity shear and calculating the right-hand side of the field 
!!equation, as well as a host of other functions.

module dist_fn
  use redistribute, only: redist_type
  implicit none
  public :: init_dist_fn, finish_dist_fn
  public :: read_parameters, wnml_dist_fn, wnml_dist_fn_species, check_dist_fn
  public :: timeadv, exb_shear, g_exb
  public :: getfieldeq, getan, getmoms, getemoms
  public :: flux, lf_flux, eexchange
  public :: get_epar, get_heat
  public :: t0, omega0, gamma0, source0
  public :: reset_init, write_f, reset_physics, write_poly
  public :: M_class, N_class, i_class, par_spectrum
  public :: l_links, r_links, itright, itleft, boundary
  public :: init_kperp2
  public :: get_jext !GGH
  public :: get_verr, get_gtran, write_fyx, collision_error
  public :: get_init_field
  public :: flux_vs_theta_vs_vpa
  public :: flux_emu, flux_e

  public :: gamtot,gamtot1,gamtot2
  public :: getmoms_notgc
  public :: mom_coeff, ncnt_mom_coeff
  public :: mom_coeff_npara, mom_coeff_nperp
  public :: mom_coeff_tpara, mom_coeff_tperp
  public :: mom_shift_para, mom_shift_perp
!CMR, 25/1/13: 
!  add public variables below for init_g
  public :: boundary_option_switch, boundary_option_linked
  public :: boundary_option_self_periodic, boundary_option_zero
  public :: init_pass_right, pass_right

  private

  ! knobs
  complex, dimension (:), allocatable :: fexp ! (nspec)
  real, dimension (:), allocatable :: bkdiff  ! (nspec)
  integer, dimension (:), allocatable :: bd_exp ! nspec
  real :: gridfac, apfac, driftknob, tpdriftknob, poisfac
  real :: t0, omega0, gamma0, source0
  real :: phi_ext, afilter, kfilter
  real :: wfb, g_exb, g_exbfac, omprimfac, btor_slab, mach
  logical :: dfexist, skexist, nonad_zero, lf_default, lf_decompose

  integer :: adiabatic_option_switch
  integer, parameter :: adiabatic_option_default = 1, &
       adiabatic_option_zero = 2, &
       adiabatic_option_fieldlineavg = 3, &
       adiabatic_option_yavg = 4

  integer :: source_option_switch
  integer, parameter :: source_option_full = 1, &
       source_option_phiext_full = 5
  
  integer :: boundary_option_switch
  integer, parameter :: boundary_option_zero = 1, &
       boundary_option_self_periodic = 2, &
       boundary_option_alternate_zero = 3, &
       boundary_option_linked = 4
  logical, public :: def_parity, even
  logical :: mult_imp, test
  logical :: accelerated_x = .false.
  logical :: accelerated_v = .false.
  logical :: increase = .true., decrease = .false.
  
!! k_parallel filter items
!  real, dimension(:), allocatable :: work, tablekp
!  real :: scale
!  integer :: nwork, ntablekp

  ! internal arrays

!  real, dimension (:,:), allocatable :: wdrift, wcurv
!  ! (-ntgrid:ntgrid, -g-layout-)

#ifdef LOWFLOW
  real, dimension (:,:), allocatable :: wcurv
  real, dimension (:,:,:), allocatable :: wstar_neo
  ! (-ntgrid:ntgrid,ntheta0,naky)
#endif

!  real, dimension (:,:,:,:,:), allocatable :: wdriftttp
  real, dimension (:,:,:,:,:,:), allocatable :: wdriftttp
  ! (-ntgrid:ntgrid,ntheta0,naky,negrid,nspec) replicated

  !> This is omega_drift multiplied by T_s/Tstar_s, 
  !! for totally trapped particles
  !! It multiples phi
  real, dimension (:,:,:,:,:,:), allocatable :: wdriftttpgen

  ! (-ntgrid:ntgrid, 2, -g-layout-)
  real, dimension (:,:,:), allocatable :: wdrift
  ! (-ntgrid:ntgrid, 2, -g-layout-)

  !> This is omega_drift multiplied by T_s/Tstar_s
  !! It multiples phi
  real, dimension (:,:,:), allocatable :: wdriftgen
  ! (-ntgrid:ntgrid, 2, -g-layout-)

  
  real, dimension (:,:,:), allocatable :: wstar
  ! (naky,negrid,nspec) replicated

  ! fieldeq
  real, dimension (:,:,:), allocatable :: gamtot, gamtot1, gamtot2, gamtot3
  ! (-ntgrid:ntgrid,ntheta0,naky) replicated

  complex, dimension (:,:,:), allocatable :: a, b, r, ainv
  ! (-ntgrid:ntgrid, 2, -g-layout-)

  real, dimension (:,:,:), allocatable :: gridfac1
  ! (-ntgrid:ntgrid,ntheta0,naky)

  complex, dimension (:,:,:), allocatable :: g0, g_h
  ! (-ntgrid:ntgrid,2, -g-layout-)

  complex, dimension (:,:,:), allocatable :: g_adj
  ! (N(links), 2, -g-layout-)

!  complex, dimension (:,:,:), allocatable, save :: gnl_1, gnl_2, gnl_3
  complex, dimension (:,:,:), allocatable, save :: gexp_1, gexp_2, gexp_3
  ! (-ntgrid:ntgrid,2, -g-layout-)

  ! momentum conservation
!  complex, dimension (:,:), allocatable :: g3int
!  real, dimension (:,:,:), allocatable :: sq

  ! exb shear
  integer, dimension(:), allocatable :: jump, ikx_indexed

  ! set_source
  real, dimension(:,:), allocatable :: ufac

  ! getfieldeq1
  real, allocatable, dimension(:,:) :: awgt
  complex, allocatable, dimension(:,:) :: fl_avg

  ! get_verr
  real, dimension (:,:), allocatable :: kmax

  ! connected bc

  integer, dimension (:,:), allocatable :: itleft, itright
  ! (naky,ntheta0)

  type :: connections_type
     integer :: iproc_left,  iglo_left
     integer :: iproc_right, iglo_right
     logical :: neighbor
  end type connections_type

  type (connections_type), dimension (:), allocatable :: connections
  ! (-g-layout-)

  ! linked only
  type (redist_type), save :: gc_from_left, gc_from_right
  type (redist_type), save :: links_p, links_h
  type (redist_type), save :: wfb_p, wfb_h
  type (redist_type), save :: pass_right

  integer, dimension (:,:), allocatable :: l_links, r_links
  integer, dimension (:,:,:), allocatable :: n_links
  logical, dimension (:,:), allocatable :: save_h
  logical :: no_comm = .false.
  integer, dimension(:), allocatable :: M_class, N_class
  integer :: i_class

  logical :: initialized = .false.
  logical :: initializing = .true.
  logical :: readinit = .false.
  logical :: bessinit = .false.
  logical :: kp2init = .false.
  logical :: connectinit = .false.
  logical :: feqinit = .false.
  logical :: lpolinit = .false.
  logical :: fyxinit = .false.
  logical :: cerrinit = .false.
  logical :: mominit = .false.

  real, allocatable :: mom_coeff(:,:,:,:)
  real, allocatable :: mom_coeff_npara(:,:,:), mom_coeff_nperp(:,:,:)
  real, allocatable :: mom_coeff_tpara(:,:,:), mom_coeff_tperp(:,:,:)
  real, allocatable :: mom_shift_para(:,:,:), mom_shift_perp(:,:,:)
  integer, parameter :: ncnt_mom_coeff=8

contains

subroutine check_dist_fn(report_unit)
  use kt_grids, only: grid_option, gridopt_box, gridopt_switch
  use nonlinear_terms, only: nonlinear_mode_switch, nonlinear_mode_on
  use species, only: spec, nspec, has_electron_species
  implicit none
  integer :: report_unit
  integer :: is 
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
    end select

    write (report_unit, fmt="('Parallel bc for passing particles at ends of the domain is:')")
    if (nonad_zero) then
       write (report_unit, fmt="(T20,'g_wesson = g_krt = 0')")
       write (report_unit, fmt="('ie NO incoming particles in the nonadiabatic piece of delta(f)')") 
    else
       write (report_unit, fmt="(T20,'g_gs2 = 0')")
       write (report_unit, fmt="('NB this ONLY gives NO incoming particles in the nonadiabatic piece of delta(f)')")
       write (report_unit, fmt="(T20,'if phi and bpar are zero at the ends of the domain')") 
    endif
    write (report_unit, *) 

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

       if (test) then
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

    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 
    write (report_unit, fmt="('The ExB parameter is ',f7.4)") g_exb
    if (abs(g_exb) .gt. epsilon(0.0)) then
       write (report_unit, fmt="('Perp shear terms will be multiplied by factor',f7.4)") g_exbfac
       write (report_unit, fmt="('Parallel shear term will be multiplied by factor',f7.4)") omprimfac
    endif

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

  end subroutine check_dist_fn

  subroutine wnml_dist_fn(unit)
  use species, only: spec, has_electron_species
  implicit none
  integer :: unit
    if (dfexist) then
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

       write (unit, fmt="(' nonad_zero = ',L1)") nonad_zero
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
             
          end select
       end if

       if (apfac /= 1.) write (unit, fmt="(' apfac = ',e16.10)") apfac
       if (driftknob /= 1.) write (unit, fmt="(' driftknob = ',e16.10)") driftknob
       if (tpdriftknob /= 1.) write (unit, fmt="(' tpdriftknob = ',e16.10)") tpdriftknob
       if (poisfac /= 0.) write (unit, fmt="(' poisfac = ',e16.10)") poisfac
       if (kfilter /= 0.) write (unit, fmt="(' kfilter = ',e16.10)") kfilter
       if (afilter /= 0.) write (unit, fmt="(' afilter = ',e16.10)") afilter
       if (mult_imp) write (unit, fmt="(' mult_imp = ',L1)") mult_imp
       if (test) write (unit, fmt="(' test = ',L1)") test
       if (def_parity) then
          write (unit, fmt="(' def_parity = ',L1)") def_parity
          if (even) write (unit, fmt="(' even = ',L1)") even
       end if
       write (unit, fmt="(' /')")
    endif
    if (skexist) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "source_knobs"
       select case (source_option_switch)

       case (source_option_full)
          write (unit, fmt="(' source_option = ',a)") '"full"'

       case(source_option_phiext_full)
          write (unit, fmt="(' source_option = ',a)") '"phiext_full"'
          write (unit, fmt="(' source0 = ',e16.10)") source0
          write (unit, fmt="(' omega0 = ',e16.10)") omega0
          write (unit, fmt="(' gamma0 = ',e16.10)") gamma0
          write (unit, fmt="(' t0 = ',e16.10)") t0
          write (unit, fmt="(' phi_ext = ',e16.10)") phi_ext
       
       end select
       write (unit, fmt="(' /')")
    endif

  end subroutine wnml_dist_fn


  subroutine wnml_dist_fn_species(unit)
     use species, only: nspec
     implicit none
     integer :: unit, i
     character (100) :: line
     do i=1,nspec
         write (unit, *)
         write (line, *) i
         write (unit, fmt="(' &',a)") &
            & trim("dist_fn_species_knobs_"//trim(adjustl(line)))
         write (unit, fmt="(' fexpr = ',e13.6)") real(fexp(i))
         write (unit, fmt="(' bakdif = ',e13.6)") bkdiff(i)
         write (unit, fmt="(' bd_exp = ',i6,'  /')") bd_exp(i)
      end do
  end subroutine wnml_dist_fn_species

  subroutine init_dist_fn
    use mp, only: proc0, finish_mp
    use species, only: init_species, nspec
    use theta_grid, only: init_theta_grid, ntgrid
    use kt_grids, only: init_kt_grids, naky, ntheta0, akx, aky
    use le_grids, only: init_le_grids, nlambda, negrid
    use run_parameters, only: init_run_parameters
    use collisions, only: init_collisions!, vnmult
    use gs2_layouts, only: init_dist_fn_layouts, init_gs2_layouts
    use nonlinear_terms, only: init_nonlinear_terms
    use hyper, only: init_hyper
    implicit none
    logical:: debug=.false.

!    write (*,*) 'entering init_dist_fn'
    if (initialized) return
    initialized = .true.

    if (debug) write(6,*) "init_dist_fn: init_gs2_layouts"
    call init_gs2_layouts

    if (debug) write(6,*) "init_dist_fn: init_species"
    call init_species

    if (debug) write(6,*) "init_dist_fn: init_le_grids"
    call init_le_grids (accelerated_x, accelerated_v)

    if (debug) write(6,*) "init_dist_fn: init_theta_grid"
    call init_theta_grid

    if (debug) write(6,*) "init_dist_fn: init_kt_grids"
    call init_kt_grids

    if (debug) write(6,*) "init_dist_fn: read_parameters"
    call read_parameters

    if (test) then
       if (proc0) then
          write (*,*) 'nspecies = ',nspec
          write (*,*) 'nlambda = ', nlambda
          write (*,*) 'negrid = ',negrid
          write (*,*) 'ntheta0 = ',ntheta0
          write (*,*) 'naky = ',naky
       end if
       call finish_mp
       stop
    end if

    if (debug) write(6,*) "init_dist_fn: run_parameters"
    call init_run_parameters

    if (debug) write(6,*) "init_dist_fn: kperp2"
    call init_kperp2

    if (debug) write(6,*) "init_dist_fn: dist_fn_layouts"
    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec)

    if (debug) write(6,*) "init_dist_fn: nonlinear_terms"
    call init_nonlinear_terms 

    if (debug) write(6,*) "init_dist_fn: allocate_arrays"
    call allocate_arrays

    if (debug) write(6,*) "init_dist_fn: init_vpar"
    call init_vpar

    if (debug) write(6,*) "init_dist_fn: init_wdrift"
    call init_wdrift

    if (debug) write(6,*) "init_dist_fn: init_wstar"
    call init_wstar

    if (debug) write(6,*) "init_dist_fn: init_bessel"
    call init_bessel

    if (debug) write(6,*) "init_dist_fn: init_par_filter"
    call init_par_filter

!    call init_vnmult (vnmult)
    if (debug) write(6,*) "init_dist_fn: init_collisions"
    call init_collisions ! needs to be after init_run_parameters

#ifdef LOWFLOW
    if (debug) write(6,*) "init_dist_fn: init_lowflow"
    call init_lowflow ! needs to be before init_invert_rhs
#endif

    if (debug) write(6,*) "init_dist_fn: init_invert_rhs"
    call init_invert_rhs

    if (debug) write(6,*) "init_dist_fn: init_fieldeq"
    call init_fieldeq

    if (debug) write(6,*) "init_dist_fn: init_hyper"
    call init_hyper

    if (debug) write(6,*) "init_dist_fn: init_mom_coeff"
    call init_mom_coeff

  end subroutine init_dist_fn

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use theta_grid, only: nperiod, shat
    use text_options, only: text_option, get_option_value
    use species, only: nspec
    use mp, only: proc0, broadcast
    use theta_grid, only: itor_over_B
    implicit none
    type (text_option), dimension (3), parameter :: sourceopts = &
         (/ text_option('default', source_option_full), &
            text_option('full', source_option_full), &
            text_option('phiext_full', source_option_phiext_full) /)
    character(20) :: source_option

    type (text_option), dimension (8), parameter :: boundaryopts = &
         (/ text_option('default', boundary_option_zero), &
            text_option('zero', boundary_option_zero), &
            text_option('unconnected', boundary_option_zero), &
            text_option('self-periodic', boundary_option_self_periodic), &
            text_option('periodic', boundary_option_self_periodic), &
            text_option('kperiod=1', boundary_option_self_periodic), &
            text_option('linked', boundary_option_linked), &
            text_option('alternate-zero', boundary_option_alternate_zero) /)
    character(20) :: boundary_option

    type (text_option), dimension (7), parameter :: adiabaticopts = &
         (/ text_option('default', adiabatic_option_default), &
            text_option('no-field-line-average-term', adiabatic_option_default), &
            text_option('field-line-average-term', adiabatic_option_fieldlineavg), &
! eventually add in iphi00 = 0 option:
            text_option('iphi00=0', adiabatic_option_default), &
            text_option('iphi00=1', adiabatic_option_default), &
            text_option('iphi00=2', adiabatic_option_fieldlineavg), &
            text_option('iphi00=3', adiabatic_option_yavg)/)
    character(30) :: adiabatic_option
            
    namelist /dist_fn_knobs/ boundary_option, nonad_zero, gridfac, apfac, &
         driftknob, tpdriftknob, poisfac, adiabatic_option, &
         kfilter, afilter, mult_imp, test, def_parity, even, wfb, &
         g_exb, g_exbfac, omprimfac, btor_slab, mach, lf_default, lf_decompose
    
    namelist /source_knobs/ t0, omega0, gamma0, source0, phi_ext, source_option
    integer :: ierr, is, in_file
    real :: bd

    if (readinit) return
    readinit = .true.

    if (proc0) then
       boundary_option = 'default'
       nonad_zero = .false.
       adiabatic_option = 'default'
       poisfac = 0.0
       gridfac = 1.0  ! used to be 5.e4
       apfac = 1.0
       driftknob = 1.0
       tpdriftknob = -9.9e9
       t0 = 100.0
       source0 = 1.0
       omega0 = 0.0
       gamma0 = 0.0
       phi_ext = 0.0
       afilter = 0.0
       kfilter = 0.0
       g_exb = 0.0
       g_exbfac = 1.0
       mach = 0.0
       omprimfac = 1.0
       btor_slab = 0.0
       wfb = 1.
       mult_imp = .false.
       test = .false.
       def_parity = .false.
       even = .true.
       lf_default = .true.
       lf_decompose = .false.
       source_option = 'default'
       in_file = input_unit_exist("dist_fn_knobs", dfexist)
!       if (dfexist) read (unit=input_unit("dist_fn_knobs"), nml=dist_fn_knobs)
       if (dfexist) read (unit=in_file, nml=dist_fn_knobs)
       if (tpdriftknob == -9.9e9) tpdriftknob=driftknob

       in_file = input_unit_exist("source_knobs", skexist)
!       if (skexist) read (unit=input_unit("source_knobs"), nml=source_knobs)
       if (skexist) read (unit=in_file, nml=source_knobs)

       if(abs(shat) <=  1.e-5) boundary_option = 'periodic'

       ierr = error_unit()
       call get_option_value &
            (boundary_option, boundaryopts, boundary_option_switch, &
            ierr, "boundary_option in dist_fn_knobs")

       call get_option_value &
            (source_option, sourceopts, source_option_switch, &
            ierr, "source_option in source_knobs")
       call get_option_value &
            (adiabatic_option, adiabaticopts, adiabatic_option_switch, &
            ierr, "adiabatic_option in dist_fn_knobs")

    end if
    if (.not.allocated(fexp)) allocate (fexp(nspec), bkdiff(nspec), bd_exp(nspec))
    if (proc0) call read_species_knobs

    call broadcast (boundary_option_switch)
    call broadcast (nonad_zero)
    call broadcast (adiabatic_option_switch)
    call broadcast (gridfac)
    call broadcast (poisfac)
    call broadcast (apfac)
    call broadcast (driftknob)
    call broadcast (tpdriftknob)
    call broadcast (t0)
    call broadcast (source0)
    call broadcast (omega0)
    call broadcast (gamma0)
    call broadcast (phi_ext)
    call broadcast (g_exb)
    call broadcast (g_exbfac)
    call broadcast (mach)
    call broadcast (omprimfac)
    call broadcast (btor_slab)
    call broadcast (afilter)
    call broadcast (kfilter)
    call broadcast (source_option_switch)
    call broadcast (fexp)
    call broadcast (bkdiff)
    call broadcast (bd_exp)
    call broadcast (mult_imp)
    call broadcast (test)
    call broadcast (def_parity)
    call broadcast (lf_default)
    call broadcast (lf_decompose)
    call broadcast (even)
    call broadcast (wfb)

    if (mult_imp) then
       ! nothing -- fine for linear runs, but not implemented nonlinearly
    else
! consistency check for bkdiff
       bd = bkdiff(1)
       do is = 1, nspec
          if (bkdiff(is) /= bd) then
             if (proc0) write(*,*) 'Forcing bkdiff for species ',is,' equal to ',bd
             if (proc0) write(*,*) 'If this is a linear run, and you want unequal bkdiff'
             if (proc0) write(*,*) 'for different species, specify mult_imp = .true.'
             if (proc0) write(*,*) 'in the dist_fn_knobs namelist.'
             bkdiff(is) = bd
          endif
       end do
! consistency check for fexp
!       fe = fexp(1)
!       do is = 1, nspec
!          if (fexp(is) /= fe) then
!             if (proc0) write(*,*) 'Forcing fexp for species ',is,' equal to ',fe
!             if (proc0) write(*,*) 'If this is a linear run, and you want unequal fexp'
!             if (proc0) write(*,*) 'for different species, specify mult_imp = .true.'
!             if (proc0) write(*,*) 'in the dist_fn_knobs namelist.'
!             fexp(is) = fe
!          endif
!       end do
    end if

! consistency check for afilter
!    if (afilter /= 0.0) then
!       if (proc0) write(*,*) 'Forcing afilter = 0.0'
!       afilter = 0.0
!    end if


!CMR, 15/2/11:  Move following lines to here from init_dist_fn, so read_parameters 
!               sets up itor_over_B
!!CMR, 19/10/10:
!! Override itor_over_B, if "dist_fn_knobs" parameter btor_slab ne 0
!! Not ideal to set geometry quantity here, but its historical! 
       if (abs(btor_slab) > epsilon(0.0)) itor_over_B = btor_slab
!! Done for slab, where itor_over_B is determined by angle between B-field 
!! and toroidal flow: itor_over_B = (d(u_z)/dx) / (d(u_y)/dx) = Btor / Bpol
!! u = u0 (phihat) = x d(u0)/dx (phihat) = x d(uy)/dx (yhat + Btor/Bpol zhat)
!! g_exb = d(uy)/dx => d(uz)/dx = g_exb * Btor/Bpol = g_exb * itor_over_B

  end subroutine read_parameters 

  subroutine read_species_knobs
    use species, only: nspec
    use file_utils, only: get_indexed_namelist_unit
    implicit none
    integer :: is, unit
    do is = 1, nspec
       fexp(is) = (0.4,0.0)
       bkdiff(is) = 0.0
       bd_exp(is) = 0
       
       call get_indexed_namelist_unit (unit, "dist_fn_species_knobs", is)
       call fill_species_knobs (unit, fexp(is), bkdiff(is), bd_exp(is))
       close (unit=unit)
    end do
  end subroutine read_species_knobs

  subroutine fill_species_knobs (unit, fexp_out, bakdif_out, bd_exp_out)
    implicit none
    integer, intent (in) :: unit
    complex, intent (in out) :: fexp_out
    real, intent (in out) :: bakdif_out
    integer, intent (in out) :: bd_exp_out
    integer :: bd_exp
    real :: fexpr, fexpi, bakdif
    namelist /dist_fn_species_knobs/ fexpr, fexpi, bakdif, bd_exp

    fexpr = real(fexp_out)
    fexpi = aimag(fexp_out)
    bakdif = bakdif_out
    bd_exp = bd_exp_out
    read (unit=unit, nml=dist_fn_species_knobs)
    fexp_out = cmplx(fexpr,fexpi)
    bd_exp_out = bd_exp
    bakdif_out = bakdif
  end subroutine fill_species_knobs

  subroutine init_wdrift
    use species, only: nspec, spec
    use general_f0, only: df0dE
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, ntheta0
    use le_grids, only: negrid, ng2, nlambda, al, jend, forbid
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use dist_fn_arrays, only: ittp
    implicit none
    integer :: ig, ik, it, il, ie, is
    integer :: iglo
!    logical :: alloc = .true.
!CMR
    logical :: debug = .false.
!CMR

! find totally trapped particles 
!    if (alloc) allocate (ittp(-ntgrid:ntgrid))
    if (.not. allocated(ittp)) allocate (ittp(-ntgrid:ntgrid))
!    ittp = 0
     ittp = nlambda+1
    do ig = -ntgrid+1, ntgrid-1
!       if (jend(ig) > 0 .and. jend(ig) <= nlambda) then
!          if (1.0-al(jend(ig))*bmag(ig+1) < 2.0*epsilon(0.0) &
!               .and. 1.0-al(jend(ig))*bmag(ig-1) < 2.0*epsilon(0.0)) &
!          then
!             ittp(ig) = jend(ig)
!          end if
!       end if

! all pitch angles greater than or equal to ittp are totally trapped or forbidden

       if (nlambda > ng2) then
          do il = ng2+1, nlambda
             if (forbid(ig-1,il) .and. forbid(ig+1, il) .and. .not. forbid(ig, il)) then
                ittp(ig) = il
                exit
             end if
          end do
       end if
    end do

    if (.not. allocated(wdrift)) then
       ! allocate wdrift with sign(vpa) dependence because will contain 
       ! Coriolis as well as magnetic drifts
       allocate (wdrift(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (wdriftgen(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (wdriftttp(-ntgrid:ntgrid,ntheta0,naky,negrid,nspec,2))
       allocate (wdriftttpgen(-ntgrid:ntgrid,ntheta0,naky,negrid,nspec,2))
    end if
    wdrift = 0.  ; wdriftgen = 0. ; wdriftttp = 0. ; wdriftttpgen = 0.
#ifdef LOWFLOW
    if (.not. allocated(wcurv)) allocate (wcurv(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    wcurv = 0.
#endif

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do ig = -ntgrid, ntgrid
!CMR, 13/3/2007:     isolate passing and trapped particle drift knobs
          il=il_idx(g_lo,iglo)
          ie=ie_idx(g_lo,iglo)
          it=it_idx(g_lo,iglo)
          ik=ik_idx(g_lo,iglo)
          is=is_idx(g_lo,iglo)
          if ( jend(ig) > 0 .and. jend(ig) <= nlambda .and. il >= ng2+1 .and. il <= jend(ig)) then
             wdrift(ig,1,iglo) = wdrift_func(ig, il, ie, it, ik, is)*tpdriftknob
#ifdef LOWFLOW
             wcurv(ig,iglo) = wcurv_func(ig, it, ik)*tpdriftknob
#endif
!CMR:  multiply trapped particle drift by tpdriftknob 
!CMR:               (tpdriftknob defaults to driftknob if not supplied)
          else
             wdrift(ig,1,iglo) = wdrift_func(ig, il, ie, it, ik, is)*driftknob
#ifdef LOWFLOW
             wcurv(ig,iglo) = wcurv_func(ig, it, ik)*driftknob
#endif
!CMR:  multiply passing particle drift by driftknob
          endif
!CMRend
          ! add Coriolis drift to magnetic drifts
          wdrift(ig,2,iglo) = wdrift(ig,1,iglo) - wcoriolis_func(ig,il,ie,it,ik,is)
          wdrift(ig,1,iglo) = wdrift(ig,1,iglo) + wcoriolis_func(ig,il,ie,it,ik,is)
       end do
    end do
    wdriftttp = 0.0
    wdriftttpgen = 0.0
    do is = 1, nspec
       do ie = 1, negrid
          do ik = 1, naky
             do it = 1, ntheta0
! moved this loop inside. 4.10.99
                do ig = -ntgrid, ntgrid
                   if (ittp(ig) == nlambda+1) cycle
                   do il = ittp(ig), nlambda
                      if (forbid(ig, il)) exit
!GWH+JAB: should this be calculated only at ittp? or for each totally trapped pitch angle? (Orig logic: there was only one totally trapped pitch angle; now multiple ttp are allowed.
                      wdriftttp(ig,it,ik,ie,is,1) &
                           = wdrift_func(ig,ittp(ig),ie,it,ik,is)*tpdriftknob
!CMR:  totally trapped particle drifts also scaled by tpdriftknob 
                      ! add Coriolis drift to magnetic drifts
                      wdriftttp(ig,it,ik,ie,is,2) = wdriftttp(ig,it,ik,ie,is,1) &
                           - wcoriolis_func(ig,il,ie,it,ik,is)
                      wdriftttp(ig,it,ik,ie,is,1) = wdriftttp(ig,it,ik,ie,is,1) &
                           + wcoriolis_func(ig,il,ie,it,ik,is)
                   end do
                end do
             end do
          end do
       end do
    end do

! This should be weighted by bakdif to be completely consistent
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do ig = -ntgrid, ntgrid-1
          wdrift(ig,:,iglo) = 0.5*(wdrift(ig,:,iglo) + wdrift(ig+1,:,iglo))
#ifdef LOWFLOW
          wcurv(ig,iglo) =  0.5*(wcurv(ig,iglo) + wcurv(ig+1,iglo))
#endif
       end do
    end do

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       wdriftgen(:,:,iglo) = &
               -wdrift(:,:,iglo)*spec(is)%temp*df0dE(ie,is)
       wdriftttpgen(:,it,ik,ie,is,:) = &
        -wdriftttp(:,it,ik,ie,is,:)*spec(is)%temp*df0dE(ie,is)
    end do
!    alloc = .false.
!CMR
    if (debug) write(6,*) 'init_wdrift: driftknob, tpdriftknob=',driftknob,tpdriftknob

  end subroutine init_wdrift

  function wdrift_func (ig, il, ie, it, ik, is)
    use theta_grid, only: bmag, gbdrift, gbdrift0, cvdrift, cvdrift0
    use theta_grid, only: shat
    use kt_grids, only: aky, theta0, akx
    use le_grids, only: energy, al
    use run_parameters, only: wunits
    use gs2_time, only: code_dt
    implicit none
    real :: wdrift_func
    integer, intent (in) :: ig, ik, it, il, ie,is

    ! note that wunits=aky/2 (for wstar_units=F)
    if (aky(ik) == 0.0) then
       wdrift_func = akx(it)/shat &
                    *(cvdrift0(ig)*energy(ie,is)*(1.0 - al(il)*bmag(ig)) &
                      + gbdrift0(ig)*0.5*energy(ie,is)*al(il)*bmag(ig)) &
                     *code_dt/2.0
    else
       wdrift_func = ((cvdrift(ig) + theta0(it,ik)*cvdrift0(ig)) &
                        *energy(ie,is)*(1.0 - al(il)*bmag(ig)) &
                      + (gbdrift(ig) + theta0(it,ik)*gbdrift0(ig)) &
                        *0.5*energy(ie,is)*al(il)*bmag(ig)) &
                     *code_dt*wunits(ik)
    end if
  end function wdrift_func

#ifdef LOWFLOW
  function wcurv_func (ig, it, ik)
    use theta_grid, only: cvdrift, cvdrift0, shat
    use kt_grids, only: aky, theta0, akx
    use run_parameters, only: wunits
    use gs2_time, only: code_dt
    implicit none
    real :: wcurv_func
    integer, intent (in) :: ig, ik, it

    if (aky(ik) == 0.0) then
       wcurv_func = akx(it)/shat &
            * cvdrift0(ig) * code_dt/2.0
    else
       wcurv_func = (cvdrift(ig) + theta0(it,ik)*cvdrift0(ig)) &
            *code_dt*wunits(ik)
    end if
  end function wcurv_func
#endif

  function wcoriolis_func (ig, il, ie, it, ik, is)

    use theta_grid, only: bmag, cdrift, cdrift0, shat
    use kt_grids, only: aky, theta0, akx
    use le_grids, only: energy, al
    use run_parameters, only: wunits
    use gs2_time, only: code_dt
    use species, only: spec

    implicit none

    real :: wcoriolis_func
    integer, intent (in) :: ig, ik, it, il, ie, is

    if (aky(ik) == 0.0) then
       wcoriolis_func = mach*sqrt(max(energy(ie,is)*(1.0-al(il)*bmag(ig)),0.0)) &
            * cdrift0(ig) * code_dt * akx(it)/(2.*shat*spec(is)%stm)
    else
       wcoriolis_func = mach*sqrt(max(energy(ie,is)*(1.0-al(il)*bmag(ig)),0.0)) &
            * (cdrift(ig) + theta0(it,ik)*cdrift0(ig))*code_dt*wunits(ik)/spec(is)%stm
    end if

  end function wcoriolis_func

  subroutine init_vpar
    use general_f0, only: zogtemp, df0dE
    use dist_fn_arrays, only: vpa, vpar, vpargen, vpac,vpacgen, vperp2
    use species, only: spec
    use theta_grid, only: ntgrid, delthet, bmag, gradpar
    use le_grids, only: energy, al
    use run_parameters, only: tunits
    use gs2_time, only: code_dt
    use gs2_layouts, only: g_lo, ik_idx, il_idx, ie_idx, is_idx
    implicit none
    integer :: iglo, ik, is, ie
    real :: al1, e1
    
    if (.not.allocated(vpa)) then
       allocate (vpa    (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (vpac   (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (vpacgen(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (vperp2 (-ntgrid:ntgrid,  g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (vpar   (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (vpargen(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    endif
    vpa = 0. ; vpac = 0. ; vpacgen = 0.; vperp2 = 0. ; vpar = 0. ; vpargen = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       al1 = al(il_idx(g_lo,iglo))
       e1 = energy(ie_idx(g_lo,iglo),is_idx(g_lo,iglo))

       vpa(:,1,iglo) = sqrt(e1*max(0.0, 1.0 - al1*bmag))
       vpa(:,2,iglo) = - vpa(:,1,iglo)
       vperp2(:,iglo) = bmag*al1*e1

       where (1.0 - al1*bmag < 100.0*epsilon(0.0))
          vpa(:,1,iglo) = 0.0
          vpa(:,2,iglo) = 0.0
       end where

! Where vpac /= 1, it could be weighted by bakdif for better consistency??
!CMR, 4/8/2011:
!CMR : vpa is parallel velocity at grid points (normalised to v_ts)
!CMR : vpac is grid centered parallel velocity
!CMR : vpar = q_s/sqrt{T_s m_s}*DELT/DTHETA * vpac |\gradpar(theta)| 
!                                     where gradpar(theta) is centered
!  ie  vpar = q_s/sqrt{T_s m_s} (v_||^GS2). \gradpar(theta)/DTHETA . DELT
! 
!   comments on vpac, vpar
!     (i) surely vpac=0 at or beyond bounce points, so WHY was it set to +-1?
!                                seems unphysical!
!    (ii) should some be weighted by bakdif?
!CMR 

       where (1.0 - al1*0.5*(bmag(-ntgrid:ntgrid-1)+bmag(-ntgrid+1:ntgrid)) &
              < 0.0)
          vpac(-ntgrid:ntgrid-1,1,iglo) = 1.0
          vpac(-ntgrid:ntgrid-1,2,iglo) = -1.0
       elsewhere
          vpac(-ntgrid:ntgrid-1,1,iglo) = &
              0.5*(vpa(-ntgrid:ntgrid-1,1,iglo) + vpa(-ntgrid+1:ntgrid,1,iglo))
          vpac(-ntgrid:ntgrid-1,2,iglo) = &
              0.5*(vpa(-ntgrid:ntgrid-1,2,iglo) + vpa(-ntgrid+1:ntgrid,2,iglo))
       end where
       vpac(ntgrid,:,iglo) = 0.0

       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       
       vpacgen(:,:,iglo) = -vpac(:,:,iglo)*spec(is)%temp *df0dE(ie,is)

       vpar(-ntgrid:ntgrid-1,1,iglo) = &
            spec(is)%zstm*tunits(ik)*code_dt &
            *0.5/delthet(-ntgrid:ntgrid-1) &
            *(abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid)))&
            *vpac(-ntgrid:ntgrid-1,1,iglo)
       vpar(-ntgrid:ntgrid-1,2,iglo) = &
            spec(is)%zstm*tunits(ik)*code_dt &
            *0.5/delthet(-ntgrid:ntgrid-1) &
            *(abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid)))&
            *vpac(-ntgrid:ntgrid-1,2,iglo)
       vpar(ntgrid,:,iglo) = 0.0

       vpargen(-ntgrid:ntgrid-1,1,iglo) = &
            zogtemp(ie,is)* spec(is)%stm*tunits(ik)*code_dt &
            *0.5/delthet(-ntgrid:ntgrid-1) &
            *(abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid)))&
            *vpac(-ntgrid:ntgrid-1,1,iglo)
       vpargen(-ntgrid:ntgrid-1,2,iglo) = &
            zogtemp(ie,is)*spec(is)%stm*tunits(ik)*code_dt &
            *0.5/delthet(-ntgrid:ntgrid-1) &
            *(abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid)))&
            *vpac(-ntgrid:ntgrid-1,2,iglo)
       vpargen(ntgrid,:,iglo) = 0.0
 

    end do

  end subroutine init_vpar

  subroutine init_wstar
    use species, only: nspec, spec
    use kt_grids, only: naky
    use le_grids, only: negrid, energy
    use run_parameters, only: wunits
    use gs2_time, only: code_dt
    use general_f0, only: f0prim

    implicit none
    integer :: ik, ie, is

    if(.not.allocated(wstar)) allocate (wstar(naky,negrid,nspec))

    do is = 1, nspec
       do ie = 1, negrid
          do ik = 1, naky
             wstar(ik,ie,is) = - code_dt*wunits(ik) &
                  * f0prim(ie,is)
          end do
       end do
    end do
  end subroutine init_wstar

  subroutine init_bessel
    use dist_fn_arrays, only: aj0, aj1, aj2, kperp2
    use species, only: spec
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, aky, akx
    use le_grids, only: energy, al
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use spfunc, only: j0, j1

    implicit none

    integer :: ig, ik, it, il, ie, is
    integer :: iglo
    real :: arg

    if (bessinit) return
    bessinit = .true.

    call init_kperp2

    allocate (aj0(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (aj1(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (aj2(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    aj0 = 0. ; aj1 = 0. ; aj2=0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do ig = -ntgrid, ntgrid
          arg = spec(is)%smz*sqrt(energy(ie,is)*al(il)/bmag(ig)*kperp2(ig,it,ik))
          aj0(ig,iglo) = j0(arg)
! CMR 17/1/06: BEWARE, j1 returns and aj1 stores J_1(x)/x (NOT J_1(x)), 
          aj1(ig,iglo) = j1(arg)
          aj2(ig,iglo) = 2.0*aj1(ig,iglo)-aj0(ig,iglo)
       end do
    end do

  end subroutine init_bessel

  subroutine init_kperp2
    use dist_fn_arrays, only: kperp2
    use species, only: spec
    use theta_grid, only: ntgrid, gds2, gds21, gds22, shat
    use kt_grids, only: naky, ntheta0, aky, theta0, akx
    implicit none
    integer :: ik, it

    if (kp2init) return
    kp2init = .true.

    allocate (kperp2(-ntgrid:ntgrid,ntheta0,naky))
    do ik = 1, naky
       if (aky(ik) == 0.0) then
         do it = 1, ntheta0
             kperp2(:,it,ik) = akx(it)*akx(it)*gds22/(shat*shat)
          end do
       else
          do it = 1, ntheta0
             kperp2(:,it,ik) = aky(ik)*aky(ik) &
                  *(gds2 + 2.0*theta0(it,ik)*gds21 &
                  + theta0(it,ik)*theta0(it,ik)*gds22)
          end do
       end if
    end do

  end subroutine init_kperp2

  subroutine init_par_filter
    use theta_grid, only: ntgrid, nperiod
    use gs2_transforms, only: init_zf
    use kt_grids, only: naky, ntheta0

    if ( naky*ntheta0 .eq. 0 ) then
       print *,"WARNING: kt_grids used in init_par_filter before initialised?"
    endif

    call init_zf (ntgrid, nperiod, ntheta0*naky)

  end subroutine init_par_filter

  subroutine par_spectrum(an, an2)

    use gs2_transforms, only: kz_spectrum
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0

    complex, dimension(:,:,:) :: an, an2    
    real :: scale

    call kz_spectrum (an, an2, ntgrid, ntheta0, naky)
    scale = 1./real(4*ntgrid**2)
    an2 = an2*scale

  end subroutine par_spectrum

  subroutine init_invert_rhs
    use mp, only: proc0
    use dist_fn_arrays, only: vpa, vpar, vpac, ittp
    use species, only: spec
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, aky
    use le_grids, only: negrid, nlambda, ng2, forbid
    use constants
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    implicit none
    integer :: iglo
    integer :: ig, ik, it, il, ie, is, isgn
    real :: wd, wdttp, vp, bd, wc

    if (.not.allocated(a)) then
       allocate (a(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (b(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (r(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (ainv(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    endif
    a = 0. ; b = 0. ; r = 0. ; ainv = 0.
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid-1
# ifdef LOWFLOW
             wd = wdrift(ig,isgn,iglo) + wstar_neo(ig,it,ik)*spec(is)%zt
             wdttp = wdriftttp(ig,it,ik,ie,is,isgn) + wstar_neo(ig,it,ik)*spec(is)%zt
# else
             wd = wdrift(ig,isgn,iglo)
             wdttp = wdriftttp(ig,it,ik,ie,is,isgn)
# endif

             ! use positive vpar because we will be flipping sign of d/dz
             ! when doing parallel field solve for -vpar
             ! also, note that vpar is (vpa bhat + v_{magnetic}) . grad theta / delthet
             ! if compile with LOWFLOW defined
             vp = vpar(ig,1,iglo)
             bd = bkdiff(is)

             ainv(ig,isgn,iglo) &
                  = 1.0/(1.0 + bd &
                  + (1.0-fexp(is))*spec(is)%tz*(zi*wd*(1.0+bd) + 2.0*vp))
             r(ig,isgn,iglo) &
                  = (1.0 - bd &
                  + (1.0-fexp(is))*spec(is)%tz*(zi*wd*(1.0-bd) - 2.0*vp)) &
                  *ainv(ig,isgn,iglo)
             a(ig,isgn,iglo) &
                  = 1.0 + bd &
                  + fexp(is)*spec(is)%tz*(-zi*wd*(1.0+bd) - 2.0*vp)
             b(ig,isgn,iglo) &
                  = 1.0 - bd &
                  + fexp(is)*spec(is)%tz*(-zi*wd*(1.0-bd) + 2.0*vp)
          
             if (nlambda > ng2) then
                ! zero out forbidden regions
                if (forbid(ig,il) .or. forbid(ig+1,il)) then
                   r(ig,isgn,iglo) = 0.0
                   ainv(ig,isgn,iglo) = 0.0
                end if
             
! CMR: set ainv=1 at lower bounce points for trapped particles
!      part of multiple trapped particle algorithm
!      NB not applicable to ttp or wfb!
                if (forbid(ig,il) .and. .not. forbid(ig+1,il)) then
                   ainv(ig,isgn,iglo) = 1.0 + ainv(ig,isgn,iglo)
                end if
             
                ! ???? mysterious mucking around with totally trapped particles
                ! part of multiple trapped particle algorithm
                if (il >= ittp(ig) .and. .not. forbid(ig, il)) then
                   ainv(ig,isgn,iglo) = 1.0/(1.0 + zi*(1.0-fexp(is))*spec(is)%tz*wdttp)
                   a(ig,isgn,iglo) = 1.0 - zi*fexp(is)*spec(is)%tz*wdttp
                   r(ig,isgn,iglo) = 0.0
                end if
             end if
          end do
       end do
    end do

    select case (boundary_option_switch)
    case (boundary_option_linked)
       call init_connected_bc
    case default
       if (.not. allocated(l_links)) then
          allocate (l_links(naky, ntheta0))
          allocate (r_links(naky, ntheta0))
          allocate (n_links(2, naky, ntheta0))
       end if
       l_links = 0;   r_links = 0;  n_links = 0
       
       i_class = 1
       if (.not. allocated(M_class)) then
          allocate (M_class(i_class))
          allocate (N_class(i_class))
       end if
       M_class = naky*ntheta0 ; N_class = 1
       
    end select

    initializing = .false.

  end subroutine init_invert_rhs

  subroutine init_connected_bc
    use theta_grid, only: ntgrid, nperiod, ntheta, theta
    use kt_grids, only: naky, ntheta0, aky, theta0
    use le_grids, only: ng2, nlambda
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use gs2_layouts, only: idx, proc_id
    use mp, only: iproc, nproc, max_allreduce, proc0
    use constants
    use redistribute, only: index_list_type, init_fill, delete_list
    implicit none
    type (index_list_type), dimension(0:nproc-1) :: to, from
    integer, dimension (0:nproc-1) :: nn_from, nn_to
    integer, dimension (3) :: to_low, from_low, to_high, from_high
    integer :: ik, it, il, ie, is, iglo, it0, itl, itr, jshift0
    integer :: ip, ipleft, ipright
    integer :: iglo_left, iglo_right, i, j, k
    integer :: iglo_star, it_star, ncell
    integer :: n, n_links_max, nn_max
    integer :: ng
    integer, dimension(naky*ntheta0) :: n_k

    if (connectinit) return
    connectinit = .true.

    ng = ntheta/2 + (nperiod-1)*ntheta

    ! jshift0 corresponds to J (not delta j) from Beer thesis (unsure about +0.01) -- MAB
    ! delta j is the number of akxs (i.e. theta0s) 'skipped' when connecting one 
    ! parallel segment to the next to satisfy the twist and shift
    ! boundary conditions. delta j = (ik - 1) * jshift0 . EGH

    if (naky > 1 .and. ntheta0 > 1) then
       jshift0 = int((theta(ng)-theta(-ng))/(theta0(2,2)-theta0(1,2)) + 0.01)
!       if (proc0) write (*,*) 'init_connected_bc: ',theta0(2,2), theta0(1,2), jshift0
    else if (naky == 1 .and. ntheta0 > 1 .and. aky(1) /= 0.0) then
       jshift0 = int((theta(ng)-theta(-ng))/(theta0(2,1)-theta0(1,1)) + 0.01)
    else
       jshift0 = 1
    end if
    
!    if (proc0) write (*,*) 'init_connected_bc: jshift0 = ',jshift0

    allocate (itleft(naky,ntheta0), itright(naky,ntheta0))
    itleft(1,:) = -1 ! No connected segments when ky = 0. EGH
    itright(1,:) = -1
    do ik = 1, naky
       do it = 1, ntheta0
          if (it > (ntheta0+1)/2) then
             ! akx is negative (-1 shift because akx(1) = 0) -- MAB
             it0 = it - ntheta0 - 1
          else
             ! akx is positive (-1 shift because akx(1) = 0) -- MAB
             it0 = it - 1
          end if

          if (ik == 1) then
             if (aky(ik) /= 0.0 .and. naky == 1) then
                ! for this case, jshift0 is delta j from Beer thesis -- MAB
                itl = it0 + jshift0
                itr = it0 - jshift0
             else
                itl = it0
                itr = it0
             end if
          else
             ! ik = 1 corresponds to aky=0, so shift index by 1
             ! (ik-1)*jshift0 is delta j from Beer thesis -- MAB
             itl = it0 + (ik-1)*jshift0
             itr = it0 - (ik-1)*jshift0
          end if

!          if (ik>1 .and. proc0) write(*,*) 'init_connected_bc: itl, itr = ',itl, itr

          ! remap to usual GS2 indices -- MAB
          if (itl >= 0 .and. itl < (ntheta0+1)/2) then
             itleft(ik,it) = itl + 1
          else if (itl + ntheta0 + 1 > (ntheta0+1)/2 &
             .and. itl + ntheta0 + 1 <= ntheta0) then
             itleft(ik,it) = itl + ntheta0 + 1
          else
             ! for periodic, need to change below -- MAB/BD
             ! j' = j + delta j not included in simulation, so can't connect -- MAB
             itleft(ik,it) = -1
          end if

          ! same stuff for j' = j - delta j -- MAB
          if (itr >= 0 .and. itr < (ntheta0+1)/2) then
             itright(ik,it) = itr + 1
          else if (itr + ntheta0 + 1 > (ntheta0+1)/2 &
             .and. itr + ntheta0 + 1 <= ntheta0) then
             itright(ik,it) = itr + ntheta0 + 1
          else
             ! for periodic, need to change below -- MAB/BD
             itright(ik,it) = -1
          end if
       end do
    end do

    allocate (connections(g_lo%llim_proc:g_lo%ulim_alloc))

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)

       ! get processors and indices for j' (kx') modes connecting
       ! to mode j (kx) so we can set up communication -- MAB

! if non-wfb trapped particle, no connections
       if (nlambda > ng2 .and. il >= ng2+2) then
          connections(iglo)%iproc_left = -1
          connections(iglo)%iglo_left = -1
          connections(iglo)%iproc_right = -1
          connections(iglo)%iglo_right = -1
       else
          if (itleft(ik,it) < 0) then
             connections(iglo)%iproc_left = -1
             connections(iglo)%iglo_left = -1
          else
             connections(iglo)%iproc_left &
                  = proc_id(g_lo,idx(g_lo,ik,itleft(ik,it),il,ie,is))
             connections(iglo)%iglo_left &
                  = idx(g_lo,ik,itleft(ik,it),il,ie,is)
          end if
          if (itright(ik,it) < 0) then
             connections(iglo)%iproc_right = -1
             connections(iglo)%iglo_right = -1
          else
             connections(iglo)%iproc_right &
                  = proc_id(g_lo,idx(g_lo,ik,itright(ik,it),il,ie,is))
             connections(iglo)%iglo_right &
                  = idx(g_lo,ik,itright(ik,it),il,ie,is)
          end if
       end if
       if (connections(iglo)%iproc_left >= 0 .or. &
            connections(iglo)%iproc_right >= 0) then
          connections(iglo)%neighbor = .true.
       else
          connections(iglo)%neighbor = .false.
       end if
    end do

    if (boundary_option_switch == boundary_option_linked) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc          
          ik = ik_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          if (connections(iglo)%iglo_left >= 0 .and. aky(ik) /= 0.0) &
               save_h (1,iglo) = .true.
          if (connections(iglo)%iglo_right >= 0 .and. aky(ik) /= 0.0) &
               save_h (2,iglo) = .true.
! wfb (linked)
          if (nlambda > ng2               .and. &
               il == ng2+1                .and. &
               connections(iglo)%neighbor .and. &
               aky(ik) /= 0.0) &
               save_h (:,iglo) = .true.
       end do

       allocate (l_links(naky, ntheta0))
       allocate (r_links(naky, ntheta0))
       allocate (n_links(2, naky, ntheta0))

       n_links_max = 0
       do it = 1, ntheta0
          do ik = 1, naky
! count the links for each region
! l_links = number of links to the left
             l_links(ik, it) = 0
             it_star = it
             do 
                if (it_star == itleft(ik, it_star)) exit
                if (itleft(ik, it_star) >= 0) then
                   l_links(ik, it) = l_links(ik, it) + 1
                   it_star = itleft(ik, it_star)
                   if (l_links(ik, it) > 5000) then
! abort by deallocating twice
                      write(*,*) 'l_links error'
                      deallocate (l_links)
                      deallocate (l_links)
                   endif
                else
                   exit
                end if
             end do
! r_links = number of links to the right
             r_links(ik, it) = 0
             it_star = it
             do 
                if (it_star == itright(ik, it_star)) exit
                if (itright(ik, it_star) >= 0) then
                   r_links(ik, it) = r_links(ik, it) + 1
                   it_star = itright(ik, it_star)
                   if (r_links(ik, it) > 5000) then
! abort by deallocating twice
                      write(*,*) 'r_links error'
                      deallocate (r_links)
                      deallocate (r_links)
                   endif
                else
                   exit
                end if
             end do
! 'n_links' complex numbers are needed to specify bc for (ik, it) region
! ignoring wfb
! n_links(1,:,:) is for v_par > 0, etc.
             if (l_links(ik, it) == 0) then
                n_links(1, ik, it) = 0
             else 
                n_links(1, ik, it) = 2*l_links(ik, it) - 1
             end if

             if (r_links(ik, it) == 0) then
                n_links(2, ik, it) = 0
             else 
                n_links(2, ik, it) = 2*r_links(ik, it) - 1
             end if
             n_links_max = max(n_links_max, n_links(1,ik,it), n_links(2,ik,it))
          end do
       end do
! wfb
       if (n_links_max > 0) n_links_max = n_links_max + 3
       
! now set up communication pattern:
! excluding wfb

       nn_to = 0
       nn_from = 0 

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          ! Exclude trapped particles (inc wfb)
          if (nlambda > ng2 .and. il >= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

          iglo_star = iglo
          do j = 1, r_links(ik, it)
             call get_right_connection (iglo_star, iglo_right, ipright)
! sender
             if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
! receiver
             if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1

             iglo_star = iglo_right
          end do
             
          iglo_star = iglo
          do j = 1, l_links(ik, it)
             call get_left_connection (iglo_star, iglo_left, ipleft)
! sender
             if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
! receiver
             if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
             iglo_star = iglo_left
          end do

       end do
       
       nn_max = maxval(nn_to)
       call max_allreduce (nn_max)
       if (nn_max == 0) then
          no_comm = .true.
          goto 200
       end if

       do ip = 0, nproc-1
          if (nn_from(ip) > 0) then
             allocate (from(ip)%first(nn_from(ip)))
             allocate (from(ip)%second(nn_from(ip)))
             allocate (from(ip)%third(nn_from(ip)))
          endif
          if (nn_to(ip) > 0) then
             allocate (to(ip)%first(nn_to(ip)))
             allocate (to(ip)%second(nn_to(ip)))
             allocate (to(ip)%third(nn_to(ip)))
          endif
       end do
       
       nn_from = 0
       nn_to = 0          

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (nlambda > ng2 .and. il >= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle
          
          iglo_star = iglo
          do j = 1, r_links(ik, it)
             call get_right_connection (iglo_star, iglo_right, ipright)
! sender
             if (ip == iproc) then
                n = nn_from(ipright) + 1
                nn_from(ipright) = n
                from(ipright)%first(n) = ntgrid
                from(ipright)%second(n) = 1
                from(ipright)%third(n) = iglo
             end if
! receiver
             if (ipright == iproc) then
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to(ip)%first(n) = j 
                to(ip)%second(n) = 1
                to(ip)%third(n) = iglo_right
             end if
             iglo_star = iglo_right
          end do
             
          iglo_star = iglo
          do j = 1, l_links(ik, it)
             call get_left_connection (iglo_star, iglo_left, ipleft)
! sender
             if (ip == iproc) then
                n = nn_from(ipleft) + 1
                nn_from(ipleft) = n
                from(ipleft)%first(n) = -ntgrid
                from(ipleft)%second(n) = 2
                from(ipleft)%third(n) = iglo
             end if
! receiver
             if (ipleft == iproc) then
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to(ip)%first(n) = j 
                to(ip)%second(n) = 2
                to(ip)%third(n) = iglo_left
             end if
             iglo_star = iglo_left
          end do
       end do

       from_low (1) = -ntgrid
       from_low (2) = 1
       from_low (3) = g_lo%llim_proc
       
       from_high(1) = ntgrid
       from_high(2) = 2
       from_high(3) = g_lo%ulim_alloc

       to_low (1) = 1
       to_low (2) = 1 
       to_low (3) = g_lo%llim_proc
       
       to_high(1) = n_links_max
       to_high(2) = 2
       to_high(3) = g_lo%ulim_alloc

       call init_fill (links_p, 'c', to_low, to_high, to, &
            from_low, from_high, from)
       
       call delete_list (from)
       call delete_list (to)
       
! take care of wfb

       nn_to = 0
       nn_from = 0 

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (il /= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle
             
          iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

! v_par > 0:
          call find_leftmost_link (iglo, iglo_right, ipright)
          do j = 1, ncell
! sender
             if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
! receiver
             if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1
             call get_right_connection (iglo_right, iglo_right, ipright)
          end do
             
! v_par < 0:
          call find_rightmost_link (iglo, iglo_left, ipleft)
          do j = 1, ncell
! sender
             if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
! receiver
             if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
             call get_left_connection (iglo_left, iglo_left, ipleft)
          end do
       end do
       
       do ip = 0, nproc-1
          if (nn_from(ip) > 0) then
             allocate (from(ip)%first(nn_from(ip)))
             allocate (from(ip)%second(nn_from(ip)))
             allocate (from(ip)%third(nn_from(ip)))
          endif
          if (nn_to(ip) > 0) then
             allocate (to(ip)%first(nn_to(ip)))
             allocate (to(ip)%second(nn_to(ip)))
             allocate (to(ip)%third(nn_to(ip)))
          endif
       end do
       
       nn_from = 0
       nn_to = 0          

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (il /= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

          iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

! v_par > 0: 
          call find_leftmost_link (iglo, iglo_right, ipright)
          do j = 1, ncell
! sender
             if (ip == iproc) then
                n = nn_from(ipright) + 1
                nn_from(ipright) = n
                from(ipright)%first(n) = ntgrid
                from(ipright)%second(n) = 1
                from(ipright)%third(n) = iglo
             end if
! receiver
             if (ipright == iproc) then
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to(ip)%first(n) = r_links(ik, it) + 1
                to(ip)%second(n) = 1
                to(ip)%third(n) = iglo_right
             end if
             call get_right_connection (iglo_right, iglo_right, ipright)
          end do
             
! v_par < 0: 
          call find_rightmost_link (iglo, iglo_left, ipleft)
          do j = 1, ncell
! sender
             if (ip == iproc) then
                n = nn_from(ipleft) + 1
                nn_from(ipleft) = n
                from(ipleft)%first(n) = -ntgrid
                from(ipleft)%second(n) = 2
                from(ipleft)%third(n) = iglo
             end if
! receiver
             if (ipleft == iproc) then
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to(ip)%first(n) = l_links(ik, it) + 1
                to(ip)%second(n) = 2
                to(ip)%third(n) = iglo_left
             end if
             call get_left_connection (iglo_left, iglo_left, ipleft)
          end do
       end do

       from_low (1) = -ntgrid
       from_low (2) = 1
       from_low (3) = g_lo%llim_proc
       
       from_high(1) = ntgrid
       from_high(2) = 2
       from_high(3) = g_lo%ulim_alloc

       to_low (1) = 1
       to_low (2) = 1 
       to_low (3) = g_lo%llim_proc
       
       to_high(1) = n_links_max
       to_high(2) = 2
       to_high(3) = g_lo%ulim_alloc

       call init_fill (wfb_p, 'c', to_low, to_high, to, from_low, from_high, from)
       
       call delete_list (from)
       call delete_list (to)
       
! n_links_max is typically 2 * number of cells in largest supercell
       allocate (g_adj (n_links_max, 2, g_lo%llim_proc:g_lo%ulim_alloc))

! now set up links_h:
! excluding wfb

       nn_to = 0
       nn_from = 0 

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (nlambda > ng2 .and. il >= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

! If this is not the first link in the chain, continue
          if (l_links(ik, it) > 0) then
! For each link to the right, do:
             iglo_star = iglo
             do j = 1, r_links(ik, it)
                call get_right_connection (iglo_star, iglo_right, ipright)
! sender
                if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
! receiver
                if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1
                iglo_star = iglo_right
             end do
          end if

          if (r_links(ik, it) > 0) then
! For each link to the left, do:
             iglo_star = iglo
             do j = 1, l_links(ik, it)
                call get_left_connection (iglo_star, iglo_left, ipleft)
! sender
                if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
! receiver
                if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
                iglo_star = iglo_left
             end do
          end if
       end do
       
       do ip = 0, nproc-1
          if (nn_from(ip) > 0) then
             allocate (from(ip)%first(nn_from(ip)))
             allocate (from(ip)%second(nn_from(ip)))
             allocate (from(ip)%third(nn_from(ip)))
          endif
          if (nn_to(ip) > 0) then
             allocate (to(ip)%first(nn_to(ip)))
             allocate (to(ip)%second(nn_to(ip)))
             allocate (to(ip)%third(nn_to(ip)))
          endif
       end do
       
       nn_from = 0
       nn_to = 0          

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (nlambda > ng2 .and. il >= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

          if (l_links(ik, it) > 0) then
! For each link to the right, do:
             iglo_star = iglo
             do j = 1, r_links(ik, it)
! get address of link
                call get_right_connection (iglo_star, iglo_right, ipright)
! sender
                if (ip == iproc) then
                   n = nn_from(ipright) + 1
                   nn_from(ipright) = n
                   from(ipright)%first(n) = ntgrid
                   from(ipright)%second(n) = 1
                   from(ipright)%third(n) = iglo
                end if
! receiver
                if (ipright == iproc) then
                   n = nn_to(ip) + 1
                   nn_to(ip) = n
                   to(ip)%first(n) = 2*l_links(ik, it) + j
                   to(ip)%second(n) = 1
                   to(ip)%third(n) = iglo_right
                end if
                iglo_star = iglo_right
             end do
          end if

          if (r_links(ik, it) > 0) then
! For each link to the left, do:
             iglo_star = iglo
             do j = 1, l_links(ik, it)
! get address of link
                call get_left_connection (iglo_star, iglo_left, ipleft)
! sender
                if (ip == iproc) then
                   n = nn_from(ipleft) + 1
                   nn_from(ipleft) = n
                   from(ipleft)%first(n) = -ntgrid
                   from(ipleft)%second(n) = 2   
                   from(ipleft)%third(n) = iglo
                end if
! receiver
                if (ipleft == iproc) then
                   n = nn_to(ip) + 1
                   nn_to(ip) = n
                   to(ip)%first(n) = 2*r_links(ik, it) + j
                   to(ip)%second(n) = 2
                   to(ip)%third(n) = iglo_left
                end if
                iglo_star = iglo_left
             end do
          end if
       end do

       from_low (1) = -ntgrid
       from_low (2) = 1
       from_low (3) = g_lo%llim_proc
       
       from_high(1) = ntgrid
       from_high(2) = 2
       from_high(3) = g_lo%ulim_alloc

       to_low (1) = 1
       to_low (2) = 1 
       to_low (3) = g_lo%llim_proc
       
       to_high(1) = n_links_max
       to_high(2) = 2
       to_high(3) = g_lo%ulim_alloc

       call init_fill (links_h, 'c', to_low, to_high, to, &
            from_low, from_high, from)
       
       call delete_list (from)
       call delete_list (to)

! now take care of wfb (homogeneous part)

       nn_to = 0
       nn_from = 0 

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (il /= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

          iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

! v_par > 0:
          call find_leftmost_link (iglo, iglo_right, ipright)
          do j = 1, ncell
! sender
             if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
! receiver
             if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1
             call get_right_connection (iglo_right, iglo_right, ipright)
          end do

! v_par < 0:
          call find_rightmost_link (iglo, iglo_left, ipleft)
          do j = 1, ncell
! sender
             if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
! receiver
             if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
             call get_left_connection (iglo_left, iglo_left, ipleft)
          end do
       end do
       
       do ip = 0, nproc-1
          if (nn_from(ip) > 0) then
             allocate (from(ip)%first(nn_from(ip)))
             allocate (from(ip)%second(nn_from(ip)))
             allocate (from(ip)%third(nn_from(ip)))
          endif
          if (nn_to(ip) > 0) then
             allocate (to(ip)%first(nn_to(ip)))
             allocate (to(ip)%second(nn_to(ip)))
             allocate (to(ip)%third(nn_to(ip)))
          endif
       end do
       
       nn_from = 0
       nn_to = 0          

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (il /= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

          iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

! v_par > 0:
          call find_leftmost_link (iglo, iglo_right, ipright)
          do j = 1, ncell
! sender
             if (ip == iproc) then
                n = nn_from(ipright) + 1
                nn_from(ipright) = n
                from(ipright)%first(n) = ntgrid
                from(ipright)%second(n) = 1
                from(ipright)%third(n) = iglo
             end if
! receiver
             if (ipright == iproc) then
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to(ip)%first(n) = 2*ncell - r_links(ik, it)
                to(ip)%second(n) = 1
                to(ip)%third(n) = iglo_right
             end if
             call get_right_connection (iglo_right, iglo_right, ipright)
          end do
 
! v_par < 0:
          call find_rightmost_link (iglo, iglo_left, ipleft)
          do j = 1, ncell
! sender
             if (ip == iproc) then
                n = nn_from(ipleft) + 1
                nn_from(ipleft) = n
                from(ipleft)%first(n) = -ntgrid
                from(ipleft)%second(n) = 2   
                from(ipleft)%third(n) = iglo
             end if
                
! receiver
             if (ipleft == iproc) then
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to(ip)%first(n) = 2*ncell - l_links(ik, it)
                to(ip)%second(n) = 2
                to(ip)%third(n) = iglo_left
             end if
             call get_left_connection (iglo_left, iglo_left, ipleft)
          end do
       end do

       from_low (1) = -ntgrid
       from_low (2) = 1
       from_low (3) = g_lo%llim_proc
       
       from_high(1) = ntgrid
       from_high(2) = 2
       from_high(3) = g_lo%ulim_alloc

       to_low (1) = 1
       to_low (2) = 1 
       to_low (3) = g_lo%llim_proc
       
       to_high(1) = n_links_max
       to_high(2) = 2
       to_high(3) = g_lo%ulim_alloc

       call init_fill (wfb_h, 'c', to_low, to_high, to, from_low, from_high, from)
       
       call delete_list (from)
       call delete_list (to)

200    continue

! Now set up class arrays for the implicit fields 
! i_class classes
! N_class(i) = number of linked cells for i_th class 
! M_class(i) = number of members in i_th class

! First count number of linked cells for each (kx, ky) 
       k = 1
       do it = 1, ntheta0
          do ik = 1, naky
             n_k(k) = 1 + l_links(ik, it) + r_links(ik, it)
             k = k + 1
!             if (proc0) write (*,*) 'init_connected_bc: ',ik, it, l_links(ik, it), r_links(ik, it)
          end do
       end do

! Count how many unique values of n_k there are.  This is the number 
! of classes.

! Sort: 
       do j = 1, naky*ntheta0-1
          do k = 1, naky*ntheta0-1                       
             if (n_k(k+1) < n_k(k)) then
                i = n_k(k)
                n_k(k) = n_k(k+1)
                n_k(k+1) = i
             end if
          end do
       end do

! Then count:
       i_class = 1
       do k = 1, naky*ntheta0-1
          if (n_k(k) == n_k(k+1)) cycle
          i_class = i_class + 1
       end do

! Allocate M, N:
       allocate (M_class(i_class))
       allocate (N_class(i_class))

! Initial values
       M_class = 1 ; N_class = 0

! Fill M, N arrays: 
       j = 1
       do k = 2, naky*ntheta0
          if (n_k(k) == n_k(k-1)) then
             M_class(j) = M_class(j) + 1
          else 
             N_class(j) = n_k(k-1)
             M_class(j) = M_class(j)/N_class(j)
             j = j + 1
          end if
       end do       
       j = i_class
       N_class(j) = n_k(naky*ntheta0)
       M_class(j) = M_class(j)/N_class(j)

! Check for consistency:
       
! j is number of linked cells in class structure
       j = 0
       do i = 1, i_class
          j = j + N_class(i)*M_class(i)
       end do

       if (j /= naky*ntheta0) then
          write(*,*) 'PE ',iproc,'has j= ',j,' k= ',naky*ntheta0,' : Stopping'
          stop
       end if

    end if

  end subroutine init_connected_bc



  subroutine init_pass_right
! CMR, 29/1/2013: simple new routine sets up "pass_right" redist_type 
! to pass g(ntgrid,1,iglo) to g(-ntgrid,1,iglo) in right neighbouring flux tube.
! To be used in init_g.
    use gs2_layouts, only: g_lo, il_idx, idx, proc_id
    use le_grids, only: ng2
    use mp, only: iproc, nproc, max_allreduce
    use redistribute, only: redist_type
    use redistribute, only: index_list_type, init_fill, delete_list
    use theta_grid, only:ntgrid 
    implicit none
    type (index_list_type), dimension(0:nproc-1) :: to, from
    integer, dimension (0:nproc-1) :: nn_from, nn_to
    integer, dimension(3) :: from_low, from_high, to_low, to_high
    integer :: il, iglo, ip, iglo_right, ipright, n, nn_max
    logical :: haslinks=.true.
    logical :: debug=.false.

      if (boundary_option_switch .eq. boundary_option_linked) then
! Need communications to satisfy || boundary conditions
! First find required blocksizes 
         nn_to = 0   ! # communicates from ip TO HERE (iproc)
         nn_from = 0 ! # communicates to ip FROM HERE (iproc)
         do iglo = g_lo%llim_world, g_lo%ulim_world
            il = il_idx(g_lo,iglo)          
            if (il > ng2+1) cycle     ! Exclude trapped as not connected
            ip = proc_id(g_lo,iglo)
! ONLY interested in right connections here  (for use in init_g)
            call get_right_connection (iglo, iglo_right, ipright)
            if (ipright .eq. iproc ) then
               nn_to(ip)=nn_to(ip)+1
               if (ip .eq.1 .and. ipright.eq.0) then
                  write(6,*) "TO 0 from 1:",ip,iglo,ipright,iglo_right
               endif
            endif
            if (ip .eq. iproc .and. ipright .ge. 0 ) then
               nn_from(ipright)=nn_from(ipright)+1
               if ( ip .eq. 1 .and. ipright.eq.0 ) then
                  write(6,*) "FROM 1 to 0:",ip,iglo,ipright,iglo_right
               endif
            endif
         end do
         nn_max = maxval(nn_to)
         call max_allreduce (nn_max)
! skip addressing if no links are needed
         if (nn_max == 0) then 
            haslinks=.false.
         endif
      endif
      if (debug) then
         write(6,*) 'init_pass_right processor, nn_to:',iproc,nn_to
         write(6,*) 'init_pass_right processor, nn_from:',iproc,nn_from
      endif

      if (boundary_option_switch.eq.boundary_option_linked .and. haslinks) then
! 
! CMR, 25/1/2013: 
!      communication required to satisfy linked BC
!      allocate indirect addresses for sends/receives 
!     
!      NB communications use "c_fill_3" as g has 3 indices
!      but 1 index sufficient as only communicating g(ntgrid,1,*)! 
!      if "c_fill_1" in redistribute we'd only need allocate: from|to(ip)%first 
!                             could be more efficient
!  
         do ip = 0, nproc-1
            if (nn_from(ip) > 0) then
               allocate (from(ip)%first(nn_from(ip)))
               allocate (from(ip)%second(nn_from(ip)))
               allocate (from(ip)%third(nn_from(ip)))
            endif
            if (nn_to(ip) > 0) then
               allocate (to(ip)%first(nn_to(ip)))
               allocate (to(ip)%second(nn_to(ip)))
               allocate (to(ip)%third(nn_to(ip)))
            endif
         end do
! Now fill the indirect addresses...
         nn_from = 0 ; nn_to = 0          
         do iglo = g_lo%llim_world, g_lo%ulim_world
            il = il_idx(g_lo,iglo)          
            if (il > ng2+1) cycle     ! Exclude disconnected trapped particles
            ip = proc_id(g_lo,iglo)
            call get_right_connection (iglo, iglo_right, ipright)
            if (ip .eq. iproc .and. ipright .ge. 0 ) then 
               n=nn_from(ipright)+1 ; nn_from(ipright)=n
               from(ipright)%first(n)=ntgrid
               from(ipright)%second(n)=1
               from(ipright)%third(n)=iglo
            endif
            if (ipright .eq. iproc ) then 
               n=nn_to(ip)+1 ; nn_to(ip)=n
               to(ip)%first(n)=-ntgrid
               to(ip)%second(n)=1
               to(ip)%third(n)=iglo_right
            endif
         end do
         if (debug) then
            write(6,*) 'init_pass_right processor, nn_to:',iproc,nn_to
            write(6,*) 'init_pass_right processor, nn_from:',iproc,nn_from
         endif

         from_low(1)=-ntgrid ; from_low(2)=1 ; from_low(3)=g_lo%llim_proc       
         from_high(1)=ntgrid; from_high(2)=2; from_high(3)=g_lo%ulim_alloc
         to_low(1)=-ntgrid  ; to_low(2)=1   ; to_low(3)=g_lo%llim_proc       
         to_high(1)=ntgrid ; to_high(2)=2  ; to_high(3)=g_lo%ulim_alloc
         call init_fill (pass_right, 'c', to_low, to_high, to, &
               from_low, from_high, from)

         call delete_list (from)
         call delete_list (to)
      endif
  end subroutine init_pass_right



  subroutine get_left_connection (iglo, iglo_left, iproc_left)
    use gs2_layouts, only: g_lo, proc_id, idx
    use gs2_layouts, only: ik_idx, it_idx, il_idx, ie_idx, is_idx
    implicit none
    integer, intent (in) :: iglo
    integer, intent (out) :: iglo_left, iproc_left
    integer :: ik, it, il, ie, is

    ik = ik_idx(g_lo,iglo)
    it = it_idx(g_lo,iglo)
    
    if (itleft(ik,it) < 0) then
       iglo_left = -1
       iproc_left = -1
       return
    end if

    il = il_idx(g_lo,iglo)
    ie = ie_idx(g_lo,iglo)
    is = is_idx(g_lo,iglo)

    iglo_left = idx(g_lo,ik,itleft(ik,it),il,ie,is)
    iproc_left = proc_id(g_lo,iglo_left)
  end subroutine get_left_connection

  subroutine get_right_connection (iglo, iglo_right, iproc_right)
    use gs2_layouts, only: g_lo, proc_id, idx
    use gs2_layouts, only: ik_idx, it_idx, il_idx, ie_idx, is_idx
    implicit none
    integer, intent (in) :: iglo
    integer, intent (out) :: iglo_right, iproc_right
    integer :: ik, it, il, ie, is

    ik = ik_idx(g_lo,iglo)
    it = it_idx(g_lo,iglo)
    
    if (itright(ik,it) < 0) then
       iglo_right = -1
       iproc_right = -1
       return
    end if

    il = il_idx(g_lo,iglo)
    ie = ie_idx(g_lo,iglo)
    is = is_idx(g_lo,iglo)

    iglo_right = idx(g_lo,ik,itright(ik,it),il,ie,is)
    iproc_right = proc_id(g_lo,iglo_right)
  end subroutine get_right_connection

  subroutine allocate_arrays
    use kt_grids, only: naky, ntheta0, box
    use theta_grid, only: ntgrid, shat
    use dist_fn_arrays, only: g, gnew, gold
    use dist_fn_arrays, only: kx_shift, theta0_shift   ! MR
    use gs2_layouts, only: g_lo
    use nonlinear_terms, only: nonlin
    implicit none
!    logical :: alloc = .true.

!    if (alloc) then
    if (.not. allocated(g)) then
       allocate (g    (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (gnew (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (g0   (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
#ifdef LOWFLOW
       allocate (gexp_1(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (gexp_2(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (gexp_3(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       gexp_1 = 0. ; gexp_2 = 0. ; gexp_3 = 0.
#else
       if (nonlin) then
          allocate (gexp_1(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          allocate (gexp_2(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          allocate (gexp_3(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          gexp_1 = 0. ; gexp_2 = 0. ; gexp_3 = 0.
       else
          allocate (gexp_1(1,2,1), gexp_2(1,2,1), gexp_3(1,2,1))
       end if
#endif
       if (boundary_option_switch == boundary_option_linked) then
          allocate (g_h(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          g_h = 0.
          allocate (save_h(2,g_lo%llim_proc:g_lo%ulim_alloc))
          save_h = .false.
       endif
       if (abs(g_exb*g_exbfac) > epsilon(0.)) then           ! MR 
          if (box .or. shat .eq. 0.0) then
             allocate (kx_shift(naky))
             kx_shift = 0.
          else
             allocate (theta0_shift(naky))
             theta0_shift = 0.
          endif
       endif                           ! MR end
    endif

    g = 0. ; gnew = 0. ; g0 = 0. !; gold = 0.

!    alloc = .false.
  end subroutine allocate_arrays

  !> This function calculates the distribution function at the next timestep. 
  !! It calculates both the inhomogenous part, gnew, due to the sources
  !! (principly the drive terms and the nonlinear term)
  !! and the homogenous part, g1. The actual evolution of the dist func
  !! is done in the subroutine invert_rhs. 
  !!
  !! After solving for the new dist funcs, this routine calls hyper_diff, which
  !! adds hyper diffusion if present, and solfp1, from the collisions module,
  !! which adds collisions if present.

  subroutine timeadv (phi, apar, bpar, phinew, aparnew, bparnew, istep, mode)

    use theta_grid, only: ntgrid
!    use collisions, only: solfp1
    use le_derivatives, only: vspace_derivatives
    use dist_fn_arrays, only: gnew, g, gold
!    use nonlinear_terms, only: add_nonlinear_terms
    use nonlinear_terms, only: add_explicit_terms
    use hyper, only: hyper_diff
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    integer, optional, intent (in) :: mode
    integer :: modep

    modep = 0
    if (present(mode)) modep = mode

    call add_explicit_terms (gexp_1, gexp_2, gexp_3, &
         phi, apar, bpar, istep, bkdiff(1), fexp(1))
    call invert_rhs (phi, apar, bpar, phinew, aparnew, bparnew, istep)
    call hyper_diff (gnew, phinew, bparnew)
    call vspace_derivatives (gnew, g, g0, phi, apar, bpar, phinew, aparnew, bparnew, modep)
!    if (istep == nstep) call write_mpdist (gnew, '.gtmp', last=.true.)

    if (def_parity) then
       if (even) then
          gnew(-ntgrid:-1, 1,:) = gnew( ntgrid: 1:-1,2,:)
          gnew( 1: ntgrid, 1,:) = gnew(-1:-ntgrid:-1,2,:)
       else
          gnew( 1: ntgrid, 1,:) = -gnew(-1:-ntgrid:-1,2,:)
          gnew(-ntgrid:-1, 1,:) = -gnew( ntgrid: 1:-1,2,:)
       end if
    end if
       
  end subroutine timeadv

! communication initializations for exb_shear should be done once and 
! redistribute routines should be used  BD
!
!  subroutine init_exb_shear
!
!    use redistribute, only: index_list_type, delete_list
!    implicit none
!    type (index_list_type), dimension(0:nproc-1) :: to, from
!    integer, dimension (0:nproc-1) :: nn_from, nn_to
!    integer, dimension (3) :: to_low, from_low, to_high, from_high
!
!  end subroutine init_exb_shear

  subroutine exb_shear (g0, phi, apar, bpar)
! MR, 2007: modified Bill Dorland's version to include grids where kx grid
!           is split over different processors
! MR, March 2009: ExB shear now available on extended theta grid (ballooning)
! CMR, May 2009: 2pishat correction factor on extended theta grid (ballooning)
!                so GEXB is same physical quantity in box and ballooning
! CMR, Oct 2010: multiply timestep by tunits(iky) for runs with wstar_units=.t.
! CMR, Oct 2010: add save statements to prevent potentially large and memory 
!                killing array allocations!
    
    use mp, only: iproc, proc0, send, receive, mp_abort
    use gs2_layouts, only: ik_idx, it_idx, g_lo, idx_local, idx, proc_id
    use run_parameters, only: tunits
    use theta_grid, only: ntgrid, ntheta, shat
    use file_utils, only: error_unit
    use kt_grids, only: akx, aky, naky, ikx, ntheta0, box, theta0
    use le_grids, only: negrid, nlambda
    use species, only: nspec
    use run_parameters, only: fphi, fapar, fbpar
    use dist_fn_arrays, only: kx_shift, theta0_shift
    use gs2_time, only: code_dt, code_dt_old
    use constants, only: twopi    

    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g0
    complex, dimension(:,:,:), allocatable :: temp 
    complex, dimension(:,:), allocatable :: temp2
    integer, dimension(1), save :: itmin
    integer :: ierr, j 
    integer :: ik, it, ie, is, il, isgn, to_iglo, from_iglo, to_iproc, from_iproc
    integer:: iib, iit, ileft, iright, i

    real, save :: dkx, dtheta0
    real :: gdt
    logical, save :: exb_first = .true.
    complex , dimension(-ntgrid:ntgrid) :: z
    character(130) :: str

    ierr = error_unit()

! MR, March 2009: remove ExB restriction to box grids
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MR, 2007: Works for ALL layouts (some layouts need no communication)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize kx_shift, jump, idx_indexed
    if (exb_first) then
       exb_first = .false.
       allocate (jump(naky)) 
       jump = 0
       if (box .or. shat .eq. 0.0) then
          allocate (ikx_indexed(ntheta0))
          itmin = minloc (ikx)
          
          do it=itmin(1), ntheta0
             ikx_indexed (it+1-itmin(1)) = it
          end do
          
          do it=1,itmin(1)-1
             ikx_indexed (ntheta0 - itmin(1) + 1 + it)= it
          end do

          if (ntheta0 .gt. 1) then 
             dkx = akx(2)-akx(1)
          else
             write(ierr,*) "exb_shear: ERROR, need ntheta0>1 for sheared flow"
          endif
       else
! MR, March 2009: on extended theta grid theta0_shift tracks ExB shear
! CMR, 25 May 2009: fix 2pi error so that: dtheta0/dt = -GEXB/shat
          if (ntheta0 .gt. 1) then 
             dtheta0 = theta0(2,1)-theta0(1,1) 
          else 
             dtheta0=twopi
          endif
       end if
    end if
    
    
    ! BD: To do: Put the right timestep in here.
    ! For now, approximate Greg's dt == 1/2 (t_(n+1) - t_(n-1))
    ! with code_dt.  
    !
    ! Note: at first time step, there is a difference of a factor of 2.
    !
 
    ! necessary to get factor of 2 right in first time step and
    ! also to get things right after changing time step size
    ! added May 18, 2009 -- MAB
    gdt = 0.5*(code_dt + code_dt_old)
    
! kx_shift is a function of time.   Update it here:  
! MR, 2007: kx_shift array gets saved in restart file
! CMR, 5/10/2010: multiply timestep by tunits(ik) for wstar_units=.true. 
    if ( box .or. shat .eq. 0.0 ) then
       do ik=1, naky
          kx_shift(ik) = kx_shift(ik) - aky(ik)*g_exb*g_exbfac*gdt*tunits(ik)
          jump(ik) = nint(kx_shift(ik)/dkx)
          kx_shift(ik) = kx_shift(ik) - jump(ik)*dkx
       end do
    else
       do ik=1, naky
          theta0_shift(ik) = theta0_shift(ik) - g_exb*g_exbfac*gdt/shat*tunits(ik)
          jump(ik) = nint(theta0_shift(ik)/dtheta0)
          theta0_shift(ik) = theta0_shift(ik) - dtheta0*jump(ik)
       enddo 
    end if

    
    if (.not. box .and. shat .ne. 0.0 ) then
! MR, March 2009: impact of ExB shear on extended theta grid computed here
!                 for finite shat
       do ik =1,naky
          j=jump(ik)
          if (j .eq. 0) cycle     
          if (abs(j) .ge. ntheta0) then
              write(str,fmt='("in exb_shear: jump(ik)=",i4," > ntheta0 =",i4," for ik=",i4". => reduce timestep or increase ntheta0")') j,ik,ntheta0
              write(ierr,*) str
              call mp_abort(str)
          endif 
          allocate(temp2(-ntgrid:ntgrid,abs(j)),temp(-ntgrid:ntgrid,2,abs(j)))
          iit=ntheta0+1-abs(j) ; iib=abs(j)
          ileft = -ntgrid+ntheta ; iright=ntgrid-ntheta

          if (fphi > epsilon(0.0)) then
             if (j < 0) then
                temp2 = phi(:,:iib,ik)
                do i=1,iit-1
                   phi(:,i,ik) = phi(:,i-j,ik)
                enddo
                phi(ileft:,iit:,ik) = temp2(:iright,:)
                phi(:ileft-1,iit:,ik) = 0.0
             else 
                temp2 = phi(:,iit:,ik)
                do i=ntheta0,iib+1,-1 
                   phi(:,i,ik) = phi(:,i-j,ik)
                enddo
                phi(:iright,:iib ,ik) = temp2(ileft:,:)
                phi(iright+1:ntgrid,:iib,:) = 0.0
             endif
          endif
          if (fapar > epsilon(0.0)) then
             if (j < 0) then
                temp2 = apar(:,:iib,ik)
                do i=1,iit-1
                   apar(:,i,ik) = apar(:,i-j,ik)
                enddo
                apar(ileft:,iit:,ik) = temp2(:iright,:)
                apar(:ileft-1,iit:,ik) = 0.0
             else 
                temp2 = apar(:,iit:,ik)
                do i=ntheta0,iib+1,-1 
                   apar(:,i,ik) = apar(:,i-j,ik)
                enddo
                apar(:iright,:iib ,ik) = temp2(ileft:,:)
                apar(iright+1:ntgrid,:iib,:) = 0.0
             endif
          endif
          if (fbpar > epsilon(0.0)) then
             if (j < 0) then
                temp2 = bpar(:,:iib,ik)
                do i=1,iit-1
                   bpar(:,i,ik) = bpar(:,i-j,ik)
                enddo
                bpar(ileft:,iit:,ik) = temp2(:iright,:)
                bpar(:ileft-1,iit:,ik) = 0.0
             else 
                temp2 = bpar(:,iit:,ik)
                do i=ntheta0,iib+1,-1 
                   bpar(:,i,ik) = bpar(:,i-j,ik)
                enddo
                bpar(:iright,:iib ,ik) = temp2(ileft:,:)
                bpar(iright+1:ntgrid,:iib,:) = 0.0
             endif
          end if

! now the distribution functions

          do is=1,nspec
             do ie=1,negrid
                do il=1,nlambda

                   if (j < 0) then
                      do it = 1, iib
                         from_iglo = idx(g_lo, ik, it, il, ie, is)
                         if (idx_local (g_lo, from_iglo)) temp(:,:,it) = g0(:,:,from_iglo)
                      end do

                      do it = 1, iit-1                        
                           to_iglo = idx(g_lo, ik, it,   il, ie, is)
                         from_iglo = idx(g_lo, ik, it-j, il, ie, is)

                         if (idx_local(g_lo, to_iglo).and. idx_local(g_lo, from_iglo)) then
                            g0(:,:,to_iglo) = g0(:,:,from_iglo)
                         else if (idx_local(g_lo, from_iglo)) then
                            do isgn = 1, 2
                               call send(g0(:, isgn, from_iglo), proc_id (g_lo, to_iglo))
                            enddo
                         else if (idx_local(g_lo, to_iglo)) then
                            do isgn = 1, 2
                               call receive(g0(:, isgn, to_iglo), proc_id (g_lo, from_iglo))
                            enddo
                         endif
                      enddo

                      do it = iit, ntheta0                     
                         to_iglo = idx(g_lo, ik, it, il, ie, is)
                         from_iglo = idx(g_lo, ik, it-j-ntheta0, il, ie, is)

                         if (idx_local(g_lo, to_iglo) .and. idx_local(g_lo, from_iglo)) then
                            g0(ileft:,:,to_iglo) = temp(:iright,:,it-iit+1)
                            g0(:ileft-1,:,to_iglo) = 0.0
                         else if (idx_local(g_lo, from_iglo)) then
                            do isgn = 1,2
                               call send(temp(:, isgn, it-iit), proc_id (g_lo, to_iglo))
                            enddo
                         else if (idx_local(g_lo, to_iglo)) then
                            do isgn=1, 2
                               call receive(z, proc_id (g_lo, from_iglo))
                               g0(ileft:,isgn,to_iglo) = z(:iright)
                               g0(:ileft-1,isgn,to_iglo) = 0.0
                            enddo
                         endif
                      enddo

                   else ! j>0

                      do it = 1, j
                         from_iglo = idx(g_lo, ik, iit+it-1, il, ie, is)
                         if (idx_local (g_lo, from_iglo)) temp(:,:,it) = g0(:,:,from_iglo)
                      end do

                      do it = ntheta0, iib+1, -1
                           to_iglo = idx(g_lo, ik, it,   il, ie, is)
                         from_iglo = idx(g_lo, ik, it-j, il, ie, is)

                         if (idx_local(g_lo, to_iglo) .and. idx_local(g_lo, from_iglo)) then
                            g0(:,:,to_iglo) = g0(:,:,from_iglo)
                         else if (idx_local(g_lo, from_iglo)) then
                            do isgn = 1, 2
                               call send(g0(:, isgn, from_iglo), proc_id (g_lo, to_iglo))
                            enddo
                         else if (idx_local(g_lo, to_iglo)) then
                            do isgn = 1, 2
                               call receive(g0(:,isgn,to_iglo), proc_id (g_lo, from_iglo))
                            enddo
                         endif
                      enddo

                      do it = 1, iib
                           to_iglo = idx(g_lo, ik, it,           il, ie, is)
                         from_iglo = idx(g_lo, ik, iit+it-1, il, ie, is)

                         if (idx_local(g_lo, to_iglo).and. idx_local(g_lo, from_iglo)) then
                            g0(:iright,:,to_iglo) = temp(ileft:,:,it)
                            g0(iright+1:,:,to_iglo) = 0.0
                         else if (idx_local(g_lo, from_iglo)) then
                            do isgn = 1, 2
                               call send(temp(:, isgn, it), proc_id (g_lo, to_iglo))
                            enddo
                         else if (idx_local(g_lo, to_iglo)) then
                            do isgn = 1, 2
                               call receive(z, proc_id (g_lo, from_iglo))
                               g0(:iright,isgn,to_iglo) = z(ileft:)
                               g0(iright+1:,isgn,to_iglo) = 0.0
                            enddo
                         endif
                      enddo
                   endif
                enddo
             enddo
          enddo
          deallocate (temp,temp2)
       enddo
    end if
    
    if (box .or. shat .eq. 0.0) then
       do ik = naky, 2, -1
          if (jump(ik) < 0) then
             if (fphi > epsilon(0.0)) then
                do it = 1, ntheta0 + jump(ik)
                   phi(:,ikx_indexed(it),ik) = phi(:,ikx_indexed(it-jump(ik)),ik)
                end do
                do it = ntheta0 + jump(ik) + 1, ntheta0
                   phi(:,ikx_indexed(it),ik) = 0.
                end do
             end if
             if (fapar > epsilon(0.0)) then
                do it = 1, ntheta0 + jump(ik)
                   apar(:,ikx_indexed(it),ik) = apar(:,ikx_indexed(it-jump(ik)),ik)
                end do
                do it = ntheta0 + jump(ik) + 1, ntheta0
                   apar (:,ikx_indexed(it),ik) = 0.
                end do
             end if
             if (fbpar > epsilon(0.0)) then 
                do it = 1, ntheta0 + jump(ik)
                   bpar(:,ikx_indexed(it),ik) = bpar(:,ikx_indexed(it-jump(ik)),ik)
                end do
                do it = ntheta0 + jump(ik) + 1, ntheta0
                   bpar (:,ikx_indexed(it),ik) = 0.
                end do
             end if
             do is=1,nspec
                do ie=1,negrid
                   do il=1,nlambda

                      do it = 1, ntheta0 + jump(ik)                        

                           to_iglo = idx(g_lo, ik, ikx_indexed(it),          il, ie, is)
                         from_iglo = idx(g_lo, ik, ikx_indexed(it-jump(ik)), il, ie, is)

                         if (idx_local(g_lo, to_iglo) .and. idx_local(g_lo, from_iglo)) then
                            g0(:,:,to_iglo) = g0(:,:,from_iglo)
                         else if (idx_local(g_lo, from_iglo)) then
                            do isgn=1, 2
                               call send (g0(:, isgn, from_iglo), proc_id (g_lo, to_iglo))
                            enddo
                         else if (idx_local(g_lo, to_iglo)) then
                            do isgn=1, 2
                               call receive (g0(:, isgn, to_iglo), proc_id (g_lo, from_iglo))
                            enddo
                         endif
                      enddo

                      do it = ntheta0 + jump(ik) + 1, ntheta0                     
                         to_iglo = idx(g_lo, ik, ikx_indexed(it), il, ie, is)
                         if (idx_local (g_lo, to_iglo)) g0(:,:,to_iglo) = 0.
                      enddo

                   enddo
                enddo
             enddo
          endif

          if (jump(ik) > 0) then 
             if (fphi > epsilon(0.0)) then
                do it = ntheta0, 1+jump(ik), -1
                   phi(:,ikx_indexed(it),ik) = phi(:,ikx_indexed(it-jump(ik)),ik)
                end do
                do it = jump(ik), 1, -1
                   phi(:,ikx_indexed(it),ik) = 0.
                end do
             end if
             if (fapar > epsilon(0.0)) then
                do it = ntheta0, 1+jump(ik), -1
                   apar(:,ikx_indexed(it),ik) = apar(:,ikx_indexed(it-jump(ik)),ik)
                end do
                do it = jump(ik), 1, -1
                   apar(:,ikx_indexed(it),ik) = 0.
                end do
             end if
             if (fbpar > epsilon(0.0)) then
                do it = ntheta0, 1+jump(ik), -1
                   bpar(:,ikx_indexed(it),ik) = bpar(:,ikx_indexed(it-jump(ik)),ik)
                end do
                do it = jump(ik), 1, -1
                   bpar(:,ikx_indexed(it),ik) = 0.
                end do
             end if
             do is=1,nspec
                do ie=1,negrid
                   do il=1,nlambda

                      do it = ntheta0, 1+jump(ik), -1

                           to_iglo = idx(g_lo, ik, ikx_indexed(it),          il, ie, is)
                         from_iglo = idx(g_lo, ik, ikx_indexed(it-jump(ik)), il, ie, is)

                         if (idx_local(g_lo, to_iglo) .and. idx_local(g_lo, from_iglo)) then
                            g0(:,:,to_iglo) = g0(:,:,from_iglo)
                         else if (idx_local(g_lo, from_iglo)) then
                            do isgn=1, 2
                               call send(g0(:, isgn, from_iglo), proc_id(g_lo, to_iglo))
                            enddo
                         else if (idx_local(g_lo, to_iglo)) then
                            do isgn=1, 2
                               call receive(g0(:, isgn, to_iglo), proc_id (g_lo, from_iglo))
                            enddo
                         endif
                      enddo

                      do it = jump(ik), 1, -1
                         to_iglo = idx(g_lo, ik, ikx_indexed(it), il, ie, is)
                         if (idx_local (g_lo, to_iglo)) g0(:,:,to_iglo) = 0.
                      enddo

                   enddo
                enddo
             enddo
          endif
       enddo
    end if
  end subroutine exb_shear
     
  subroutine get_source_term &
       (phi, apar, bpar, phinew, aparnew, bparnew, istep, &
        isgn, iglo, sourcefac, source)
#ifdef LOWFLOW
    use dist_fn_arrays, only: hneoc, vparterm, wdfac, wstarfac, wdttpfac
#endif
    use dist_fn_arrays, only: aj0, aj1, vperp2, vpar, vpargen, vpac, vpacgen, g, ittp
    use theta_grid, only: ntgrid, theta, bmag
    use kt_grids, only: aky, akx
    use le_grids, only: nlambda, ng2, lmax, anon, energy, negrid, forbid
    use species, only: spec, nspec
    use run_parameters, only: fphi, fapar, fbpar, wunits, tunits
    use gs2_time, only: code_dt
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use nonlinear_terms, only: nonlin
    use hyper, only: D_res
    use constants
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    integer, intent (in) :: isgn, iglo
    complex, intent (in) :: sourcefac
    complex, dimension (-ntgrid:), intent (out) :: source

    integer :: ig, ik, it, il, ie, is
    complex, dimension (-ntgrid:ntgrid) :: phigavg, apargavg

!    call timer (0, 'get_source_term')

    ik = ik_idx(g_lo,iglo)
    it = it_idx(g_lo,iglo)
    il = il_idx(g_lo,iglo)
    ie = ie_idx(g_lo,iglo)
    is = is_idx(g_lo,iglo)
!CMR, 4/8/2011
! apargavg and phigavg combine to give the GK EM potential chi. 
!          chi = phigavg - apargavg*vpa(:,isgn,iglo)*spec(is)%stm
! phigavg  = phi J0 + 2 T_s/q_s . vperp^2 bpar/bmag J1/Z
! apargavg = apar J0 
! Both quantities are decentred in time and evaluated on || grid points
!
    phigavg  = (fexp(is)*phi(:,it,ik)   + (1.0-fexp(is))*phinew(:,it,ik)) &
                *aj0(:,iglo)*fphi &
             + (fexp(is)*bpar(:,it,ik) + (1.0-fexp(is))*bparnew(:,it,ik))&
                *aj1(:,iglo)*fbpar*2.0*vperp2(:,iglo)*spec(is)%tz
    apargavg = (fexp(is)*apar(:,it,ik)  + (1.0-fexp(is))*aparnew(:,it,ik)) &
                *aj0(:,iglo)*fapar

! source term in finite difference equations
    select case (source_option_switch)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Default choice: solve self-consistent equations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (source_option_full)
       if (il <= lmax) then
          call set_source
       else
          source = 0.0
       end if       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve self-consistent terms + include external i omega_d * F_0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case(source_option_phiext_full)
       if (il <= lmax) then
          call set_source
          if (aky(ik) < epsilon(0.0)) then
             source(:ntgrid-1) = source(:ntgrid-1) &
                  - zi*anon(ie)*wdriftgen(:ntgrid-1,isgn,iglo) &
                  *2.0*phi_ext*sourcefac*aj0(:ntgrid-1,iglo)
          end if
       else
          source = 0.0
       end if

    end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do matrix multiplications
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (il <= lmax) then
       if (isgn == 1) then
          do ig = -ntgrid, ntgrid-1
             source(ig) = source(ig) &
                  + b(ig,1,iglo)*g(ig,1,iglo) + a(ig,1,iglo)*g(ig+1,1,iglo)
          end do
       else
          do ig = -ntgrid, ntgrid-1
             source(ig) = source(ig) &
                  + a(ig,2,iglo)*g(ig,2,iglo) + b(ig,2,iglo)*g(ig+1,2,iglo)
          end do
       end if
    end if

    source(ntgrid) = source(-ntgrid)  ! is this necessary?  BD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! special source term for totally trapped particles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (source_option_switch == source_option_full .or. &
        source_option_switch == source_option_phiext_full) then
       if (nlambda > ng2 .and. isgn == 2) then
          do ig = -ntgrid, ntgrid
             if (il < ittp(ig)) cycle
             source(ig) &
                  = g(ig,2,iglo)*a(ig,2,iglo) &
#ifdef LOWFLOW
                  - anon(ie)*zi*(wdttpfac(ig,it,ik,ie,is,2)*hneoc(ig,2,iglo))*phigavg(ig) &
                  + zi*wstar(ik,ie,is)*hneoc(ig,2,iglo)*phigavg(ig)
#else
                  - anon(ie)*zi*(wdriftttpgen(ig,it,ik,ie,is,2))*phigavg(ig) &
                  ! Where is the PVG term? 
                  + zi*wstar(ik,ie,is)*phigavg(ig)
#endif             
          end do

          if (source_option_switch == source_option_phiext_full .and. &
               aky(ik) < epsilon(0.0)) then
             do ig = -ntgrid, ntgrid
                if (il < ittp(ig)) cycle             
                source(ig) = source(ig) - zi*anon(ie)* &
#ifdef LOWFLOW
                     wdttpfac(ig,it,ik,ie,is,isgn)*2.0*phi_ext*sourcefac*aj0(ig,iglo)
#else
                     wdriftttpgen(ig,it,ik,ie,is,isgn)*2.0*phi_ext*sourcefac*aj0(ig,iglo)
#endif
             end do
          endif

#ifdef LOWFLOW
          select case (istep)
          case (0)
             ! nothing
          case (1)
             do ig = -ntgrid, ntgrid
                if (il < ittp(ig)) cycle
                source(ig) = source(ig) + 0.5*code_dt*gexp_1(ig,isgn,iglo)
             end do
          case (2) 
             do ig = -ntgrid, ntgrid
                if (il < ittp(ig)) cycle
                source(ig) = source(ig) + 0.5*code_dt*( &
                     1.5*gexp_1(ig,isgn,iglo) - 0.5*gexp_2(ig,isgn,iglo))
             end do
          case default
             do ig = -ntgrid, ntgrid
                if (il < ittp(ig)) cycle
                source(ig) = source(ig) + 0.5*code_dt*( &
                     (23./12.)*gexp_1(ig,isgn,iglo) &
                     - (4./3.)  *gexp_2(ig,isgn,iglo) &
                     + (5./12.) *gexp_3(ig,isgn,iglo))
             end do
          end select
#else
! add in nonlinear terms 
          if (nonlin) then         
             select case (istep)
             case (0)
                ! nothing
             case (1)
                do ig = -ntgrid, ntgrid
                   if (il < ittp(ig)) cycle
                   source(ig) = source(ig) + 0.5*code_dt*gexp_1(ig,isgn,iglo)
                end do
             case (2) 
                do ig = -ntgrid, ntgrid
                   if (il < ittp(ig)) cycle
                   source(ig) = source(ig) + 0.5*code_dt*( &
                        1.5*gexp_1(ig,isgn,iglo) - 0.5*gexp_2(ig,isgn,iglo))
                end do                   
             case default
                do ig = -ntgrid, ntgrid
                   if (il < ittp(ig)) cycle
                   source(ig) = source(ig) + 0.5*code_dt*( &
                          (23./12.)*gexp_1(ig,isgn,iglo) &
                        - (4./3.)  *gexp_2(ig,isgn,iglo) &
                        + (5./12.) *gexp_3(ig,isgn,iglo))
                end do
             end select
          end if
#endif

       end if
    end if

!    call timer (1, 'get_source_term')
  contains

    subroutine set_source

      use species, only: spec
      use theta_grid, only: itor_over_B
      use mp, only: proc0

      complex :: apar_p, apar_m, phi_p, phi_m!, bpar_p !GGH added bpar_p
!      real, dimension(:,:), allocatable, save :: ufac
      real :: bd, bdfac_p, bdfac_m
      integer :: i_e, i_s
!      logical :: first = .true.

!      if (first) then
      if (.not. allocated(ufac)) then
!         first = .false.
         allocate (ufac(negrid, nspec))
         do i_e = 1, negrid
            do i_s = 1, nspec
               ufac(i_e, i_s) = (2.0*spec(i_s)%uprim &
                    + spec(i_s)%uprim2*energy(i_e,i_s)**(1.5)*sqrt(pi)/4.0)
            end do
         end do
      endif

! try fixing bkdiff dependence
      bd = bkdiff(1)

      bdfac_p = 1.+bd*(3.-2.*real(isgn))
      bdfac_m = 1.-bd*(3.-2.*real(isgn))


      do ig = -ntgrid, ntgrid-1
!CMR, 4/8/2011:
! Some concerns, may be red herrings !
! (1) no bakdif factors in phi_m, apar_p, apar_m, vpar !!! 
!                        (RN also spotted this for apar_p)
! (2) can interpolations of products be improved? 
!
!  Attempt at variable documentation:
! phigavg  = phi J0 + 2 T_s/q_s . vperp^2 bpar/bmag J1/Z
! apargavg = apar J0                        (decentered in t) 
! NB apargavg and phigavg combine to give the GK EM potential chi
! chi = phigavg - apargavg*vpa(:,isgn,iglo)*spec(is)%stm
! phi_p = 2 phigavg                      .... (roughly!)
! phi_m = d/dtheta (phigavg)*DTHETA 
! apar_p = 2 apargavg  
! apar_m = 2 d/dt (apar)*DELT  (gets multiplied later by J0 and vpa when included in source)
! => phi_p - apar_p*vpa(:,isgn,iglo)*spec(is)%stm = 2 chi  .... (roughly!)  
! vparterm = -2.0*vpar (IN ABSENCE OF LOWFLOW TERMS)
! wdfac = wdrift + wcoriolis/spec(is)%stm (IN ABSENCE OF LOWFLOW TERMS)
! wstarfac = wstar  (IN ABSENCE OF LOWFLOW TERMS)
! Line below is not true for alphas. vpar = q_s/Tstar * v_|| * the rest. EGH/GW
! vpar = q_s/sqrt{T_s m_s} (v_||^GS2). \gradpar(theta)/DTHETA . DELT (centred)
! wdrift =    q_s/T_s  v_d.\grad_perp . DELT 
! wcoriolis = q_s/T_s  v_C.\grad_perp . DELT 
!
! Definition of source:= 2*code_dt*RHS of GKE
! source     appears to contain following physical terms
!   -2q_s/T_s v||.grad(J0 phi + 2 vperp^2 bpar/bmag J1/Z T_s/q_s).delt 
!   -2d/dt(q v|| J0 apar / T).delt
!   +hyperviscosity
!   -2 v_d.\grad_perp (q J0 phi/T + 2 vperp^2 bpar/bmag J1/Z).delt 
!   -coriolis terms
!   2{\chi,f_{0s}}.delt  (allowing for sheared flow)
!CMRend

         phi_p = bdfac_p*phigavg(ig+1)+bdfac_m*phigavg(ig)
         phi_m = phigavg(ig+1)-phigavg(ig)
         ! RN> bdfac factors seem missing for apar_p
         apar_p = apargavg(ig+1)+apargavg(ig)
         apar_m = aparnew(ig+1,it,ik)+aparnew(ig,it,ik) & 
              -apar(ig+1,it,ik)-apar(ig,it,ik)

!MAB, 6/5/2009:
! added the omprimfac source term arising with equilibrium flow shear  
!CMR, 26/11/2010:
! Useful to understand that all source terms propto aky are specified here 
! using Tref=mref vtref^2. See 
! [note by BD and MK on "Microinstabilities in Axisymmetric Configurations"].
! This is converted to  the standard internal gs2 normalisation, 
! Tref=(1/2) mref vtref^2, by wunits, which contains a crucial factor 1/2.
! (Would be less confusing if always used same Tref!)
!
#ifdef LOWFLOW
         source(ig) = anon(ie)*(vparterm(ig,isgn,iglo)*phi_m &
              -spec(is)%zstm*vpac(ig,isgn,iglo) &
              *((aj0(ig+1,iglo) + aj0(ig,iglo))*0.5*apar_m  &
              + D_res(it,ik)*apar_p) &
              -zi*wdfac(ig,isgn,iglo)*phi_p) &
              + zi*(wstarfac(ig,isgn,iglo) &
              + vpac(ig,isgn,iglo)*code_dt*wunits(ik)*ufac(ie,is) &
              -2.0*omprimfac*vpac(ig,isgn,iglo)*code_dt*wunits(ik)*g_exb*itor_over_B(ig)/spec(is)%stm) &
              *(phi_p - apar_p*spec(is)%stm*vpac(ig,isgn,iglo))
#else
         source(ig) = anon(ie)*(-2.0*vpargen(ig,isgn,iglo)*phi_m &
               ! Note line below still needs to be changed for alphas
              -spec(is)%zstm*vpacgen(ig,isgn,iglo) &
              *((aj0(ig+1,iglo) + aj0(ig,iglo))*0.5*apar_m  &
              + D_res(it,ik)*apar_p) &
              -zi*wdriftgen(ig,isgn,iglo)*phi_p) &
              + zi*(wstar(ik,ie,is) &
              + vpac(ig,isgn,iglo)*code_dt*wunits(ik)*ufac(ie,is) &
              -2.0*omprimfac*vpacgen(ig,isgn,iglo)*code_dt*wunits(ik)*g_exb*itor_over_B(ig)/spec(is)%stm) &
              *(phi_p - apar_p*spec(is)%stm*vpac(ig,isgn,iglo))
#endif
      end do
        
#ifdef LOWFLOW
      select case (istep)
      case (0)
         ! nothing
      case (1)
         do ig = -ntgrid, ntgrid-1
            source(ig) = source(ig) + 0.5*code_dt*gexp_1(ig,isgn,iglo)
         end do
      case (2) 
         do ig = -ntgrid, ntgrid-1
            source(ig) = source(ig) + 0.5*code_dt*( &
                 1.5*gexp_1(ig,isgn,iglo) - 0.5*gexp_2(ig,isgn,iglo))
         end do
      case default
         do ig = -ntgrid, ntgrid-1
            source(ig) = source(ig) + 0.5*code_dt*( &
                 (23./12.)*gexp_1(ig,isgn,iglo) &
                 - (4./3.)  *gexp_2(ig,isgn,iglo) &
                 + (5./12.) *gexp_3(ig,isgn,iglo))
         end do
      end select
#else
! add in nonlinear terms 
      if (nonlin) then         
         select case (istep)
         case (0)
            ! nothing
         case (1)
            do ig = -ntgrid, ntgrid-1
               source(ig) = source(ig) + 0.5*code_dt*gexp_1(ig,isgn,iglo)
            end do
         case (2) 
            do ig = -ntgrid, ntgrid-1
               source(ig) = source(ig) + 0.5*code_dt*( &
                    1.5*gexp_1(ig,isgn,iglo) - 0.5*gexp_2(ig,isgn,iglo))
            end do
         case default
            do ig = -ntgrid, ntgrid-1
               source(ig) = source(ig) + 0.5*code_dt*( &
                      (23./12.)*gexp_1(ig,isgn,iglo) &
                    - (4./3.)  *gexp_2(ig,isgn,iglo) &
                    + (5./12.) *gexp_3(ig,isgn,iglo))
            end do
         end select
      end if
#endif

    end subroutine set_source

  end subroutine get_source_term

  subroutine invert_rhs_1 &
       (phi, apar, bpar, phinew, aparnew, bparnew, istep, &
        iglo, sourcefac)
    use dist_fn_arrays, only: gnew, ittp, vperp2, aj1, aj0
    use run_parameters, only: eqzip, secondary, tertiary, harris
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, ng2, lmax, forbid, anon
    use kt_grids, only: aky, ntheta0
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use prof, only: prof_entering, prof_leaving
    use run_parameters, only: fapar, fbpar, fphi, ieqzip
    use species, only: spec
    use general_f0, only: df0dE
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    integer, intent (in) :: iglo
    complex, intent (in) :: sourcefac

    integer :: ig, ik, it, il, ie, isgn, is
    integer :: ilmin
    complex :: beta1
    complex, dimension (-ntgrid:ntgrid,2) :: source, g1, g2
    complex :: adjleft, adjright
    logical :: kperiod_flag, speriod_flag
    integer :: ntgl, ntgr

    call prof_entering ("invert_rhs_1", "dist_fn")

    ik = ik_idx(g_lo,iglo)
    it = it_idx(g_lo,iglo)
    il = il_idx(g_lo,iglo)
    ie = ie_idx(g_lo,iglo)
    is = is_idx(g_lo,iglo)

    if(ieqzip(it_idx(g_lo,iglo),ik_idx(g_lo,iglo))==0) return
    if (eqzip) then
       if (secondary .and. ik == 2 .and. it == 1) return ! do not evolve primary mode
       if (tertiary .and. ik == 1) then
          if (it == 2 .or. it == ntheta0) return ! do not evolve periodic equilibrium
       end if
       if (harris .and. ik == 1) return ! do not evolve primary mode       
    end if

    do isgn = 1, 2
       call get_source_term (phi, apar, bpar, phinew, aparnew, bparnew, &
            istep, isgn, iglo, sourcefac, source(:,isgn))
    end do

    ! gnew is the inhomogeneous solution
    gnew(:,:,iglo) = 0.0

!CMR, 18/4/2012:
! What follows is a selectable improved parallel bc for passing particles.
!                                            (prompted by Greg Hammett)
! Original bc is: g_gs2 = gnew = 0 at ends of domain:
!   ONLY implies zero incoming particles in nonadiabatic delta f if phi=bpar=0
! Here ensure ZERO incoming particles in nonadiabatic delta f at domain ends
!  (exploits code used in subroutine g_adjust to transform g_wesson to g_gs2)
    if ( nonad_zero ) then 
       if (il <= ng2+1) then
          !This ensures that we only apply the new boundary condition to the leftmost
          !cell for sign going from left to right
          if (l_links(ik,it) .eq. 0) then
             adjleft = anon(ie)*2.0*vperp2(-ntgrid,iglo)*aj1(-ntgrid,iglo) &
                  *bparnew(-ntgrid,it,ik)*fbpar &
                  - spec(is)%z*anon(ie)*phinew(-ntgrid,it,ik)*aj0(-ntgrid,iglo) &
                  *df0dE(ie,is)*fphi
             gnew(-ntgrid,1,iglo) = gnew(-ntgrid,1,iglo) - adjleft
          end if
          !This ensures that we only apply the new boundary condition to the rightmost
          !cell for sign going from right to left
          if (r_links(ik,it) .eq. 0) then
             adjright = anon(ie)*2.0*vperp2(ntgrid,iglo)*aj1(ntgrid,iglo) &
                  *bparnew(ntgrid,it,ik)*fbpar &
                  - spec(is)%z*anon(ie)*phinew(ntgrid,it,ik)*aj0(ntgrid,iglo) &
                  *df0dE(ie,is)*fphi
             gnew(ntgrid,2,iglo) = gnew(ntgrid,2,iglo) - adjright
         end if
       endif
    endif

    ! g1 is the homogeneous solution
    g1 = 0.0

!CMR, 17/4/2012:
!  kperiod_flag = T  iff one of following applies:   
!                 boundary_option_self_periodic
!     OR          boundary_option_linked
!     OR          aky=0
!  if kperiod_flag = T, compute homogeneous solution (g1) for passing.
!        otherwise ONLY compute inhomogenous solution for passing particles.
!
!  speriod_flag = T  iff boundary_option_linked AND aky=0
!  if speriod_flag = T, apply self-periodic bcs for passing and wfb.

    select case (boundary_option_switch)
    case (boundary_option_self_periodic)
       kperiod_flag = .true.
    case (boundary_option_linked)
       kperiod_flag = .true.
       speriod_flag = aky(ik) == 0.0
    case default
       kperiod_flag = .false.
    end select

    kperiod_flag = kperiod_flag .or. aky(ik) == 0.0

    ntgl = -ntgrid
    ntgr =  ntgrid

! ng2+1 is WFB

    if (kperiod_flag) then
       if (il <= ng2+1) then
          g1(-ntgrid,1) = 1.0
          g1( ntgrid,2) = 1.0
       end if
    end if

    if (il == ng2+1) then
       g1(-ntgrid,1) = wfb  ! wfb should be unity here; variable is for testing
       g1( ntgrid,2) = wfb  ! wfb should be unity here; variable is for testing
    end if

!CMR
! g2 simply stores trapped homogeneous boundary conditions at bounce points
!     g2(iub,2) = 1.0        iub is UPPER bounce point, trapped bc for vpar < 0
!     g2(ilb,1) = g1(ilb,2)  ilb is LOWER bounce point, trapped bc for vpar > 0
! otherwise g2 is zero
!      -- g2 = 0 for all passing particles
!      -- g2 = 0 for wfb as forbid always false
!      -- g2 = 0 for ttp too
! g2 NOT used for totally trapped particles at il=nlambda 

! NB with trapped particles lmax=nlambda-1, and ttp is at il=nlambda

    g2 = 0.0
    if (nlambda > ng2 .and. il >= ng2+1 .and. il <= lmax) then
!CMR, surely il >= ng2+2 works, as forbid gives no bounce point for wfb?
!      ... but does not matter
       do ig=-ntgrid,ntgrid-1
!CMR: set g2=1 at UPPER bounce point for trapped (not ttp or wfb) at vpar<0
          if (forbid(ig+1,il).and..not.forbid(ig,il)) g2(ig,2) = 1.0
       end do
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! time advance vpar < 0  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

! inhomogeneous part: gnew
! r=ainv=0 if forbid(ig,il) or forbid(ig+1,il), so gnew=0 in forbidden
! region and at upper bounce point
    do ig = ntgrid-1, -ntgrid, -1
       gnew(ig,2,iglo) = -gnew(ig+1,2,iglo)*r(ig,2,iglo) + ainv(ig,2,iglo)*source(ig,2)
    end do

    if (kperiod_flag) then
       ilmin = 1
    else
       ilmin = ng2 + 1              !!! ilmin = ng2 + 2
    end if

! time advance vpar < 0 homogeneous part: g1
!CMR, 17/4/2012: computes homogeneous solutions for il > ilmin
!                il > ilmin includes trapped particles
!                AND passing particles IF kperiod_flag = T

    if (il >= ilmin) then
       do ig = ntgrid-1, -ntgrid, -1
          g1(ig,2) = -g1(ig+1,2)*r(ig,2,iglo) + g2(ig,2)
       end do
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! time advance vpar > 0   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! First set BCs for trapped particles at lower bounce point
    ! (excluding wfb and ttp)

! GWH & JAB: see if this fixes a numerical instability, related to close-to-wfb:
!CMR, 17/4/2012:  il<=lmax excludes ttp 
!CMR, 1/8/2011:  Not sure why wfb is included in following loop
    if (nlambda > ng2 .and. il >= ng2+1 .and. il <= lmax) then
!    if (nlambda > ng2 .and. il >= ng2+2 .and. il <= lmax) then
       ! match boundary conditions at lower bounce point
       do ig = -ntgrid, ntgrid-1
          if (forbid(ig,il) .and. .not. forbid(ig+1,il)) then
!CMR, 17/4/2012: set g2=(ig+1,1) = g1(ig+1,2) where ig+1 is LOWER bounce point
!     (previously g2 was set to 1 at ig just LEFT of lower bounce point
!      but this was handled consistently in integration of g1)
!  
             g2(ig+1,1) = g1(ig+1,2)
!CMR: init_invert_rhs  sets ainv=1 at lower bounce point for trapped 
!     source below ensures gnew(lb,1,iglo)=gnew(lb,2,iglo)
!     where lb is the lower bounce point.
             source(ig,1) = gnew(ig+1,2,iglo)  
          end if
       end do
    end if

    ! time advance vpar > 0 inhomogeneous part
    if (il <= lmax) then
       do ig = -ntgrid, ntgrid-1
          gnew(ig+1,1,iglo) = -gnew(ig,1,iglo)*r(ig,1,iglo) + ainv(ig,1,iglo)*source(ig,1)
       end do
    end if

    ! BD: Are there ever totally trapped particles at ig = ntgrid? 
    ! balancing totally trapped particles
    do ig = -ntgrid, ntgrid
       if (il >= ittp(ig)) then
          if (forbid(ig,il)) then
             gnew(ig,1,iglo) = 0.0
          else
             gnew(ig,1,iglo) = gnew(ig,2,iglo)
          end if
       end if
    end do

    ! time advance vpar > 0 homogeneous part
    if (il >= ilmin) then
       do ig = -ntgrid, ntgrid-1
!CMR, 17/4/2012:  use consistent homogeneous trapped solution (g2) at lbp
          g1(ig+1,1) = -g1(ig,1)*r(ig,1,iglo) + g2(ig+1,1)
       end do
    end if

!CMR, 17/4/2012:  
! self_periodic bc applied as follows:
! boundary_option_linked = T
!      passing + wfb if aky=0
!      isolated wfb (ie no linked domains) 
! boundary_option_linked = F
!      passing and wfb particles if boundary_option_self_periodic = T
!      wfb 
    if (boundary_option_switch == boundary_option_linked) then
       if (speriod_flag .and. il <= ng2+1) then
          call self_periodic
       else 
       ! save homogeneous solution as necessary
          if (save_h (1, iglo)) g_h(:,1,iglo) = g1(:,1)
          if (save_h (2, iglo)) g_h(:,2,iglo) = g1(:,2)
       end if

! wfb (isolated)
       if (il == ng2+1 .and. .not. connections(iglo)%neighbor) &
            call self_periodic

    else       
       ! add correct amount of homogeneous solution now
       if (kperiod_flag .and. il <= ng2+1) then
          call self_periodic
       else if (il == ng2 + 1) then
          call self_periodic
       end if

    end if

    ! add correct amount of homogeneous solution for trapped particles to satisfy boundary conditions
!JAB include wfb (ng2+1)
 !orig   if (il >= ng2+2 .and. il <= lmax) then
    if (il >= ng2+1 .and. il <= lmax) then
       beta1 = 0.0
       do ig = ntgr-1, ntgl, -1
          if (ittp(ig) <= il) cycle
          if (forbid(ig,il)) then
             beta1 = 0.0
             cycle !CMR: to avoid pointless arithmetic later in loop
          else if (forbid(ig+1,il)) then
             beta1 = (gnew(ig,1,iglo) - gnew(ig,2,iglo))/(1.0 - g1(ig,1))
          end if
          gnew(ig,1,iglo) = gnew(ig,1,iglo) + beta1*g1(ig,1)
          gnew(ig,2,iglo) = gnew(ig,2,iglo) + beta1*g1(ig,2)
       end do
    end if

    if (def_parity) then
       if (even) then
          gnew(ntgl:-1,1,iglo) = gnew( ntgr:1:-1,2,iglo)
          gnew(1:ntgr, 1,iglo) = gnew(-1:ntgl:-1,2,iglo)
       else
          gnew(1:ntgr, 1,iglo) = -gnew(-1:ntgl:-1,2,iglo)
          gnew(ntgl:-1,1,iglo) = -gnew( ntgr:1:-1,2,iglo)
       end if
    end if

    ! zero out spurious gnew outside trapped boundary
    where (forbid(:,il))
       gnew(:,1,iglo) = 0.0
       gnew(:,2,iglo) = 0.0
    end where

    call prof_leaving ("invert_rhs_1", "dist_fn")

  contains

    subroutine self_periodic
!CMR: sets sum of homogeneous and inhomogeneous solutions to give a result
!     gnew(ntgr,2) = gnew(ntgl,2)
!     gnew(ntgr,1) = gnew(ntgl,1) 
! ie periodic bcs at the ends of the flux tube.
      if (g1(ntgr,1) /= 1.) then
         beta1 = (gnew(ntgr,1,iglo) - gnew(ntgl,1,iglo))/(1.0 - g1(ntgr,1))
         gnew(:,1,iglo) = gnew(:,1,iglo) + beta1*g1(:,1)
      end if

      if (g1(ntgl,2) /= 1.) then
         beta1 = (gnew(ntgl,2,iglo) - gnew(ntgr,2,iglo))/(1.0 - g1(ntgl,2))
         gnew(:,2,iglo) = gnew(:,2,iglo) + beta1*g1(:,2)
      end if
      
    end subroutine self_periodic

  end subroutine invert_rhs_1

  subroutine invert_rhs_linked &
       (phi, apar, bpar, phinew, aparnew, bparnew, istep, sourcefac)
    use dist_fn_arrays, only: gnew
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, ng2
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx
    use redistribute, only: fill
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    complex, intent (in) :: sourcefac

    complex :: b0, fac, facd
    integer :: il, ik, it, n, i, j
    integer :: iglo, ncell

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       call invert_rhs_1 (phi, apar, bpar, phinew, aparnew, bparnew, &
            istep, iglo, sourcefac)
    end do

    if (no_comm) then
       ! nothing
    else       
       g_adj = 0.
       call fill (links_p, gnew, g_adj)
       call fill (links_h, g_h, g_adj)
       call fill (wfb_p, gnew, g_adj)
       call fill (wfb_h, g_h, g_adj)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)

          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle
! wfb
          if (nlambda > ng2 .and. il == ng2+1) then
             if (save_h(1, iglo)) then

                facd = 1.0
                do j = 1, ncell
                   facd = facd * g_adj(ncell+j,1,iglo)
                end do
                facd = 1./(1.-facd)
                
                b0 = 0.
                do i = 1, ncell-1
                   fac = 1.0
                   do j = i+1, ncell
                      fac = fac * g_adj(ncell+j,1,iglo)
                   end do
                   b0 = b0 + fac * g_adj(ncell+1-i,1,iglo)
                end do
                b0 = (b0 + g_adj(1,1,iglo))*facd

                do i = 1, l_links(ik, it)
                   b0 = b0 * g_adj(ncell+i,1,iglo) + g_adj(ncell+1-i,1,iglo)
                end do
                
                gnew(:,1,iglo) = gnew(:,1,iglo) + b0*g_h(:,1,iglo)
             endif

             if (save_h(2, iglo)) then

                facd = 1.0
                do j = 1, ncell
                   facd = facd * g_adj(ncell+j,2,iglo)
                end do
                facd = 1./(1.-facd)
                
                b0 = 0.
                do i = 1, ncell-1
                   fac = 1.0
                   do j = i+1, ncell
                      fac = fac * g_adj(ncell+j,2,iglo)
                   end do
                   b0 = b0 + fac * g_adj(ncell+1-i,2,iglo)
                end do
                b0 = (b0 + g_adj(1,2,iglo))*facd

                do i = 1, r_links(ik, it)
                   b0 = b0 * g_adj(ncell+i,2,iglo) + g_adj(ncell+1-i,2,iglo)
                end do
                
                gnew(:,2,iglo) = gnew(:,2,iglo) + b0*g_h(:,2,iglo)
             end if
          else
!
! n_links is the number of complex numbers required to fix the boundary 
! conditions in each cell that is a member of a supercell with at least two
! cells, and for which the bounce point is not at theta=pi.
!
             if (save_h(1, iglo)) then
                n = n_links(1, ik, it)
                b0 = 0.0
                do i = 1, l_links(ik, it)
                   fac = 1.0
                   do j = 1, i-1
                      fac = fac * g_adj(n+1-j, 1, iglo)
                   end do
                   b0 = b0 + g_adj(i,1,iglo) * fac
                end do
                
                gnew(:,1,iglo) = gnew(:,1,iglo) + b0*g_h(:,1,iglo)
             end if

             if (save_h(2, iglo)) then
                n = n_links(2, ik, it)
                b0 = 0.0
                do i = 1, r_links(ik, it)
                   fac = 1.0
                   do j = 1, i-1
                      fac = fac * g_adj(n+1-j, 2, iglo)
                   end do
                   b0 = b0 + g_adj(i,2,iglo) * fac
                end do
                
                gnew(:,2,iglo) = gnew(:,2,iglo) + b0*g_h(:,2,iglo)
             end if
          end if

       end do
       
    end if

  end subroutine invert_rhs_linked

  subroutine invert_rhs (phi, apar, bpar, phinew, aparnew, bparnew, istep)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use gs2_time, only: code_time
    use constants
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep

    integer :: iglo

    real :: time
    complex :: sourcefac

    call prof_entering ("invert_rhs", "dist_fn")

    time = code_time
    if (time > t0) then
       sourcefac = source0*exp(-zi*omega0*time+gamma0*time)
    else
       sourcefac = (0.5 - 0.5*cos(pi*time/t0))*exp(-zi*omega0*time+gamma0*time)
    end if

    select case (boundary_option_switch)
    case (boundary_option_linked)
       call invert_rhs_linked &
            (phi, apar, bpar, phinew, aparnew, bparnew, istep, sourcefac) 
    case default
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          call invert_rhs_1 (phi, apar, bpar, phinew, aparnew, bparnew, &
               istep, iglo, sourcefac)
       end do
    end select

    call prof_leaving ("invert_rhs", "dist_fn")
  end subroutine invert_rhs

  subroutine getan (antot, antota, antotp)
    use dist_fn_arrays, only: vpa, vperp2, aj0, aj1, gnew
    use dist_fn_arrays, only: kperp2
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_species
    use run_parameters, only: beta, fphi, fapar, fbpar
    use prof, only: prof_entering, prof_leaving
    use gs2_layouts, only: g_lo
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (out) :: antot, antota, antotp
    real, dimension (nspec) :: wgt

    integer :: isgn, iglo, ig

    call prof_entering ("getan", "dist_fn")

    antot=0. ; antota=0. ; antotp=0.

    if (fphi > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             do ig=-ntgrid, ntgrid
                g0(ig,isgn,iglo) = aj0(ig,iglo)*gnew(ig,isgn,iglo)
             end do
          end do
       end do

       wgt = spec%z*spec%dens
       call integrate_species (g0, wgt, antot)

!    if (kfilter > epsilon(0.0)) call par_filter(antot)
       if (afilter > epsilon(0.0)) antot = antot * exp(-afilter**4*kperp2**2/4.)

    end if

    if (fapar > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             do ig=-ntgrid, ntgrid
                g0(ig,isgn,iglo) = aj0(ig,iglo)*vpa(ig,isgn,iglo)*gnew(ig,isgn,iglo)
             end do
          end do
       end do
       
       wgt = 2.0*beta*spec%z*spec%dens*sqrt(spec%temp/spec%mass)
       call integrate_species (g0, wgt, antota)
!    if (kfilter > epsilon(0.0)) call par_filter(antota)

    end if

    if (fbpar > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             do ig=-ntgrid, ntgrid
                g0(ig,isgn,iglo) = aj1(ig,iglo)*vperp2(ig,iglo)*gnew(ig,isgn,iglo)
             end do
          end do
       end do
       wgt = spec%temp*spec%dens
       call integrate_species (g0, wgt, antotp)
!    if (kfilter > epsilon(0.0)) call par_filter(antotp)

    end if
    call prof_leaving ("getan", "dist_fn")
  end subroutine getan

  subroutine getmoms (ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)
    use general_f0, only: zogtemp
    use dist_fn_arrays, only: vpa, vperp2, aj0, aj1, gnew, g_adjust
    use gs2_layouts, only: is_idx, ie_idx, g_lo, ik_idx, it_idx
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment, anon, energy
    use prof, only: prof_entering, prof_leaving
    use run_parameters, only: fphi, fbpar
    use fields_arrays, only: phinew, bparnew

    implicit none
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: density, &
         upar, tpar, tperp, ntot, qparflux, pperpj1, qpperpj1

    integer :: ik, it, isgn, ie, is, iglo

! returns moment integrals to PE 0
    call prof_entering ("getmoms", "dist_fn")

! DJA+CMR: 17/1/06, use g_adjust routine to extract g_wesson
!                   from gnew, phinew and bparnew.
!           nb  <delta_f> = g_wesson J0 - q phi/T F_m  where <> = gyroaverage
!           ie  <delta_f>/F_m = g_wesson J0 - q phi/T
!
! use g0 as dist_fn dimensioned working space for all moments
! (avoid making more copies of gnew to save memory!)
!
! set gnew = g_wesson, but return gnew to entry state prior to exit 
    call g_adjust(gnew,phinew,bparnew,fphi,fbpar)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! integrate moments over the nonadiabatic part of <delta_f>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set g0= J0 g_wesson = nonadiabatic piece of <delta_f>/F_m 
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo) ; is = is_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo) ; it = it_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo)
       end do
    end do

! CMR: density is the nonadiabatic piece of perturbed density
! NB normalised wrt equ'm density for species s: n_s n_ref  
!    ie multiply by (n_s n_ref) to get abs density pert'n
    call integrate_moment (g0, density)

! DJA/CMR: upar and tpar moments 
! (nb adiabatic part of <delta f> does not contribute to upar, tpar or tperp)
! NB UPAR is normalised to vt_s = sqrt(T_s/m_s) vt_ref
!    ie multiply by spec(is)%stm * vt_ref to get abs upar
    g0 = vpa*g0
    call integrate_moment (g0, upar)

    g0 = 2.*vpa*g0
    call integrate_moment (g0, tpar)
! tpar transiently stores ppar, nonadiabatic perturbed par pressure 
!      vpa normalised to: sqrt(2 T_s T_ref/m_s m_ref)
!  including factor 2 in g0 product ensures 
!     ppar normalised to: n_s T_s n_ref T_ref 
!                         ppar = tpar + density, and so:
    tpar = tpar - density
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = vperp2(:,iglo)*gnew(:,isgn,iglo)*aj0(:,iglo)
       end do
    end do
    call integrate_moment (g0, tperp)
! tperp transiently stores pperp, nonadiabatic perturbed perp pressure
!                          pperp = tperp + density, and so:
    tperp = tperp - density
! NB TPAR, and TPERP are normalised by T_s T_ref
!    ie multiply by T_s T_ref to get abs TPAR, TPERP

! Now compute QPARFLUX
! NB QPARFLUX is normalised to n_s n_ref T_s T_ref v_ts
!    ie multiply by (n_s T_s spec(is)%stm) n_ref T_ref vt_ref to get abs qparflux
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = vpa(:,isgn,iglo)*gnew(:,isgn,iglo)*aj0(:,iglo)*energy(ie_idx(g_lo,iglo),is_idx(g_lo,iglo))
       end do
    end do 
    call integrate_moment (g0, qparflux)
   
! Now compute PPERPJ1, a modified p_perp which gives particle flux from Bpar
! NB PPERPJ1 is normalised to (n_s T_s/q_s)  n_ref T_ref/q_ref 
!    ie multiply by (n_s spec(is)%tz) n_ref T_ref/q_ref to get abs PPERPJ1
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) &
               = gnew(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo)*spec(is)%tz
       end do
    end do
    call integrate_moment (g0, pperpj1)

! Now compute QPPERPJ1, a modified p_perp*energy which gives heat flux from Bpar
! NB QPPERPJ1 is normalised to (n_s T_s^2/q_s)  n_ref  T_ref^2/q_ref
!    ie multiply by (n_s T_s spec(is)%tz) n_ref T_ref^2/q_ref to get abs QPPERPJ1
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       g0(:,:,iglo) = g0(:,:,iglo)*energy(ie_idx(g_lo,iglo),is_idx(g_lo,iglo))
    end do
    call integrate_moment (g0, qpperpj1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now include the adiabatic part of <delta f>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set g0 = <delta_f>/F_m, including the adiabatic term
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo) ; is = is_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo) ; it = it_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo) - anon(ie)*phinew(:,it,ik)*zogtemp(ie,is)
       end do
    end do

! total perturbed density
    call integrate_moment (g0, ntot)

!CMR, now multiply by species dependent normalisations to leave only reference normalisations
    do is=1,nspec
       ntot(:,:,:,is)=ntot(:,:,:,is)*spec(is)%dens
       density(:,:,:,is)=density(:,:,:,is)*spec(is)%dens
       upar(:,:,:,is)=upar(:,:,:,is)*spec(is)%stm
       tpar(:,:,:,is)=tpar(:,:,:,is)*spec(is)%temp
       tperp(:,:,:,is)=tperp(:,:,:,is)*spec(is)%temp
       qparflux(:,:,:,is)=qparflux(:,:,:,is)*spec(is)%dens*spec(is)%temp*spec(is)%stm
       pperpj1(:,:,:,is)=pperpj1(:,:,:,is)*spec(is)%dens*spec(is)%tz
       qpperpj1(:,:,:,is)=qpperpj1(:,:,:,is)*spec(is)%dens*spec(is)%temp*spec(is)%tz
    end do

! return gnew to its initial state, the variable evolved in GS2
    call g_adjust(gnew,phinew,bparnew,-fphi,-fbpar)

    call prof_leaving ("getmoms", "dist_fn")
  end subroutine getmoms

  subroutine getemoms (ntot, tperp)
    use general_f0, only: zogtemp
    use dist_fn_arrays, only: vperp2, aj0, gnew, g_adjust
    use gs2_layouts, only: is_idx, ie_idx, g_lo, ik_idx, it_idx
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment, anon
    use prof, only: prof_entering, prof_leaving
    use run_parameters, only: fphi, fbpar
    use fields_arrays, only: phinew, bparnew

    implicit none
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: tperp, ntot

    integer :: ik, it, isgn, ie, is, iglo

! returns electron density and Tperp moment integrals to PE 0
    call prof_entering ("getemoms", "dist_fn")
!
! What are <delta_f> and g_wesson in the note below?
! g_wesson = <delta_f> + q phi/T    [ignore F_m factors for simplicity]
!
! Electrostatically (for simplicity), g_adjust produces:
!
! h = g_gs2 + q <phi> / T
! 
! then in the subsequent code they calculate for ntot:
!
! ntot = integral[   J0 h - q phi / T  ]
!
! so g_wesson == h.  What is odd in our notation is the LHS.  
! We typically indicate the perturbed distribution at fixed spatial position 
! by delta_f.  In the note below, they must mean delta_f = delta_f (R), so that 
! they are gyro-averaging to get the distribution at fixed spatial position.
!
! In summary, DJA and CMR are calculating the moments at fixed spatial position
! rather than at fixed guiding centers, and using different notation than appears
! in most of our papers.
!
! DJA+CMR: 17/1/06, use g_adjust routine to extract g_wesson
!                   from gnew, phinew and bparnew.
!           nb  <delta_f> = g_wesson J0 - q phi/T F_m  where <> = gyroaverage
!           ie  <delta_f>/F_m = g_wesson J0 - q phi/T
!
! use g0 as dist_fn dimensioned working space for all moments
! (avoid making more copies of gnew to save memory!)
!
! set gnew = g_wesson, but return gnew to entry state prior to exit 
    call g_adjust(gnew, phinew, bparnew, fphi, fbpar)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc       
       ie = ie_idx(g_lo,iglo) ; is = is_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo) ; it = it_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo) - anon(ie)*phinew(:,it,ik)*zogtemp(ie,is)
       end do
    end do

! total perturbed density
    call integrate_moment (g0, ntot)

! vperp**2 moment:
    do iglo = g_lo%llim_proc, g_lo%ulim_proc       
       ie = ie_idx(g_lo,iglo) ; is = is_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo) ; it = it_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo)*vperp2(:,iglo) - anon(ie)*phinew(:,it,ik)*zogtemp(ie,is)*vperp2(:,iglo)
       end do
    end do

! total perturbed perp pressure
    call integrate_moment (g0, tperp)

! tperp transiently stores pperp, 
!       pperp = tperp + density, and so:
    tperp = tperp - ntot

    do is=1,nspec
       ntot(:,:,:,is)=ntot(:,:,:,is)*spec(is)%dens
       tperp(:,:,:,is)=tperp(:,:,:,is)*spec(is)%temp
    end do

! return gnew to its initial state, the variable evolved in GS2
    call g_adjust(gnew,phinew,bparnew,-fphi,-fbpar)

    call prof_leaving ("getemoms", "dist_fn")
  end subroutine getemoms
  
  ! moment at not guiding center coordinate
  subroutine getmoms_notgc (dens, upar, tpar, tper, ntot, jpar)
    use general_f0, only: zogtemp, df0dE
    use dist_fn_arrays, only: vpa, vperp2, aj0, aj1, gnew
    use gs2_layouts, only: g_lo, is_idx, ik_idx, it_idx, ie_idx
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use kt_grids, only: nakx => ntheta0, naky
    use le_grids, only: integrate_moment
    use mp, only: iproc
    use fields_arrays, only: phinew, bparnew
    implicit none
    complex, intent (out) :: dens(-ntgrid:,:,:,:), upar(-ntgrid:,:,:,:), &
         & tpar(-ntgrid:,:,:,:), tper(-ntgrid:,:,:,:)
    complex, intent (out), optional :: ntot(-ntgrid:,:,:,:)
    complex, intent (out), optional :: jpar(-ntgrid:,:,:)
    integer :: isgn, iglo, is, ie

    real :: a, b, tpar2, tper2
    integer :: it, ik, ig

! returns moment integrals to PE 0

! not guiding center n_total
    if(present(ntot)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)

          do isgn = 1, 2
             g0(:,isgn,iglo) = aj0(:,iglo) * gnew(:,isgn,iglo)
          end do
          do isgn = 1, 2
             g0(:,isgn,iglo) = g0(:,isgn,iglo) + phinew(:,it,ik) &
                  & *(aj0(:,iglo)**2-1.0) * zogtemp(ie,is)
          end do
          do isgn = 1, 2
             g0(:,isgn,iglo) = g0(:,isgn,iglo) &
                  & - 2.*vperp2(:,iglo)*aj1(:,iglo)*aj0(:,iglo) &
                  & * bparnew(:,it,ik) * spec(is)%temp &
                  & * df0dE(ie,is)
          end do
       end do
       call integrate_moment (g0, ntot, 1)
    endif

! not guiding center density
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g0, dens, 1)

! not guiding center upar
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = aj0(:,iglo)*vpa(:,isgn,iglo)*gnew(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g0, upar, 1)

! not guiding center tpar
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = 2.*aj0(:,iglo)*vpa(:,isgn,iglo)**2*gnew(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g0, tpar, 1)
    
! not guiding center tperp
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) = 2.*vperp2(:,iglo)*aj1(:,iglo)*gnew(:,isgn,iglo) 
       end do
    end do

    call integrate_moment (g0, tper, 1)

    do ig=-ntgrid,ntgrid
       do it=1,nakx
          do ik=1,naky
             do is=1,nspec
                tpar2=tpar(ig,it,ik,is) &
                     & - dens(ig,it,ik,is)*mom_coeff_npara(it,ik,is)
                tper2=tper(ig,it,ik,is) &
                     & - dens(ig,it,ik,is)*mom_coeff_nperp(it,ik,is)

                a=mom_coeff_tperp(it,ik,is)
                b=mom_coeff_tpara(it,ik,is)

                tpar(ig,it,ik,is)=(   tpar2-a*tper2)/(1.-a*b)
                tper(ig,it,ik,is)=(-b*tpar2+  tper2)/(1.-a*b)
             end do
          end do
       end do
    end do

    do is=1,nspec
       dens(:,:,:,is)=dens(:,:,:,is)*spec(is)%dens
       upar(:,:,:,is)=upar(:,:,:,is)*spec(is)%stm
       tpar(:,:,:,is)=tpar(:,:,:,is)*spec(is)%temp
       tper(:,:,:,is)=tper(:,:,:,is)*spec(is)%temp
    end do

    if(present(jpar)) then
       jpar(:,:,:)=cmplx(0.,0.)
       do is=1,nspec
          jpar(:,:,:)=jpar(:,:,:)+spec(is)%z*spec(is)%dens*upar(:,:,:,is)
       end do
    endif
  end subroutine getmoms_notgc

  subroutine init_fieldeq
    use dist_fn_arrays, only: aj0, aj1, vperp2, kperp2
    use species, only: nspec, spec, has_electron_species
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0, aky
    use le_grids, only: anon, integrate_species
    use gs2_layouts, only: g_lo, ie_idx, is_idx
    use run_parameters, only: tite
    use general_f0, only: df0dE
    implicit none
    integer :: iglo, isgn
    integer :: ik, it, ie, is
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: tot
    real, dimension (nspec) :: wgt

    if (feqinit) return
    feqinit = .true.

    allocate (gridfac1(-ntgrid:ntgrid,ntheta0,naky))
    gridfac1 = 1.0
    select case (boundary_option_switch)
    case (boundary_option_self_periodic)
       ! nothing
    case (boundary_option_linked)
       do it = 1, ntheta0
          do ik = 1, naky
             if (aky(ik) == 0.0) cycle
             if (itleft(ik,it) < 0) gridfac1(-ntgrid,it,ik) = gridfac
             if (itright(ik,it) < 0) gridfac1(ntgrid,it,ik) = gridfac
          end do
       end do
    case default
       do ik = 1, naky
          if (aky(ik) == 0.0) cycle
          gridfac1(-ntgrid,:,ik) = gridfac
          gridfac1(ntgrid,:,ik) = gridfac
       end do
    end select

    allocate (gamtot(-ntgrid:ntgrid,ntheta0,naky))
    allocate (gamtot1(-ntgrid:ntgrid,ntheta0,naky))
    allocate (gamtot2(-ntgrid:ntgrid,ntheta0,naky))
    if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
       allocate (gamtot3(-ntgrid:ntgrid,ntheta0,naky))
    endif
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) = -(1.0 - aj0(:,iglo)**2)*df0dE(ie,is)
       end do
    end do
    wgt = spec%z*spec%z*spec%dens
    call integrate_species (g0, wgt, tot)
    gamtot = real(tot) + kperp2*poisfac
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) = -aj0(:,iglo)*aj1(:,iglo) &
               *2.0*vperp2(:,iglo)*df0dE(ie,is)
       end do
    end do
    wgt = spec%z*spec%dens
    call integrate_species (g0, wgt, tot)
    gamtot1 = real(tot)
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) = -aj1(:,iglo)**2*2.0*vperp2(:,iglo)**2*df0dE(ie,is)
       end do
    end do
    wgt = spec%temp*spec%dens
    call integrate_species (g0, wgt, tot)
    gamtot2 = real(tot)

! adiabatic electrons 
    if (.not. has_electron_species(spec)) then
       if (adiabatic_option_switch == adiabatic_option_yavg) then
          do ik = 1, naky
             if (aky(ik) > epsilon(0.0)) gamtot(:,:,ik) = gamtot(:,:,ik) + tite
          end do
       elseif (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
          gamtot  = gamtot + tite
          gamtot3 = (gamtot-tite) / gamtot
          where (gamtot3 < 2.*epsilon(0.0)) gamtot3 = 1.0
       else
          gamtot = gamtot + tite 
       endif
    endif
  end subroutine init_fieldeq

  subroutine getfieldeq1 (phi, apar, bpar, antot, antota, antotp, &
       fieldeq, fieldeqa, fieldeqp)
    use dist_fn_arrays, only: kperp2
    use theta_grid, only: ntgrid, bmag, delthet, jacob
    use kt_grids, only: naky, ntheta0, aky
    use run_parameters, only: fphi, fapar, fbpar
    use run_parameters, only: beta, tite
    use species, only: spec, has_electron_species
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: antot, antota, antotp
    complex, dimension (-ntgrid:,:,:), intent (out) ::fieldeq,fieldeqa,fieldeqp

    integer :: ik, it
!    logical :: first = .true.
    
!    if (first) allocate (fl_avg(ntheta0, naky))
    if (.not. allocated(fl_avg)) allocate (fl_avg(ntheta0, naky))
    fl_avg = 0.

    if (.not. has_electron_species(spec)) then
       if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
          
!          if (first) then 
          if (.not. allocated(awgt)) then
             allocate (awgt(ntheta0, naky))
             awgt = 0.
             do ik = 1, naky
                do it = 1, ntheta0
                   if (aky(ik) > epsilon(0.0)) cycle
                   awgt(it,ik) = 1.0/sum(delthet*jacob*gamtot3(:,it,ik))
                end do
             end do
          endif
           
          do ik = 1, naky
             do it = 1, ntheta0
                fl_avg(it,ik) = tite*sum(delthet*jacob*antot(:,it,ik)/gamtot(:,it,ik))*awgt(it,ik)
             end do
          end do

!          do it=2,ntheta0
!             fl_avg(it, 1) = fl_avg(it, 1) - 1.0
!          end do
       end if

    end if

    if (fphi > epsilon(0.0)) then
       fieldeq = antot + bpar*gamtot1 - gamtot*gridfac1*phi 

       if (.not. has_electron_species(spec)) then
          do ik = 1, naky
             do it = 1, ntheta0
                fieldeq(:,it,ik) = fieldeq(:,it,ik) + fl_avg(it,ik)
             end do
          end do
       end if
    end if

    if (fapar > epsilon(0.0)) then
       fieldeqa = antota - kperp2*gridfac1*apar
    end if
! bpar == delta B_parallel / B_0(theta) b/c of the factor of 1/bmag(theta)**2
! in the following
    if (fbpar > epsilon(0.0)) then
       fieldeqp = (antotp+bpar*gamtot2+0.5*phi*gamtot1)*beta*apfac
       do ik = 1, naky
          do it = 1, ntheta0
             fieldeqp(:,it,ik) = fieldeqp(:,it,ik)/bmag(:)**2
          end do
       end do
       fieldeqp = fieldeqp + bpar*gridfac1
    end if

!    first = .false.

  end subroutine getfieldeq1

  subroutine getfieldeq (phi, apar, bpar, fieldeq, fieldeqa, fieldeqp)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (out) ::fieldeq,fieldeqa,fieldeqp
    complex, dimension (:,:,:), allocatable :: antot, antota, antotp

    allocate (antot (-ntgrid:ntgrid,ntheta0,naky))
    allocate (antota(-ntgrid:ntgrid,ntheta0,naky))
    allocate (antotp(-ntgrid:ntgrid,ntheta0,naky))

    call getan (antot, antota, antotp)
    call getfieldeq1 (phi, apar, bpar, antot, antota, antotp, &
         fieldeq, fieldeqa, fieldeqp)

    deallocate (antot, antota, antotp)
  end subroutine getfieldeq
  
! MAB> ported from agk
! TT> Given initial distribution function this obtains consistent fields
!CMR, 1/8/2011> corrections below for inhomogeneous bmag
  subroutine get_init_field (phi, apar, bpar)
    ! inverts the field equations:
    !   gamtot * phi - gamtot1 * bpar = antot
    !   kperp2 * apar = antota
    !   beta/2 * gamtot1 * phi + (beta * gamtot2 + 1) * bpar = - beta * antotp
    ! I haven't made any check for use_Bpar=T case.
    use run_parameters, only: beta, fphi, fapar, fbpar
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use dist_fn_arrays, only: aj0, vpa, kperp2

    complex, dimension (-ntgrid:,:,:), intent (out) :: phi, apar, bpar
    real, dimension (-ntgrid:ntgrid,ntheta0,naky) :: denominator
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: antot, antota, antotp
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: numerator
    real, dimension (-ntgrid:ntgrid,ntheta0,naky) :: bmagsp

    phi=0. ; apar=0. ; bpar=0.
    antot=0.0 ; antota=0.0 ; antotp=0.0
!CMR, 1/8/2011:  bmagsp is 3D array containing bmag
    bmagsp=spread(spread(bmag,2,ntheta0),3,naky)
    call getan (antot, antota, antotp)

    ! get phi
    if (fphi > epsilon(0.0)) then

!CMR, 1/8/2011:  bmag corrections here: 
       numerator = (beta * gamtot2 + bmagsp**2) * antot - (beta * gamtot1) * antotp
       denominator = (beta * gamtot2 + bmagsp**2) * gamtot + (beta/2.0) * gamtot1 * gamtot1

       where (abs(denominator) < epsilon(0.0)) ! it == ik == 1 only
          phi = 0.0
       elsewhere
          phi = numerator / denominator
       end where
    end if

    ! get apar
    if (fapar > epsilon(0.0)) then
       denominator = kperp2
       where (abs(denominator) < epsilon(0.0)) ! it == ik == 1 only
          apar = 0.0
       elsewhere
          apar = antota / denominator
       end where
    end if

    ! get bpar
    if (fbpar > epsilon(0.0)) then
!CMR>  bmag corrections here
       numerator = - (beta * gamtot) * antotp - (beta/2.0) * gamtot1 * antot
       denominator = gamtot * (beta * gamtot2 + bmagsp**2) + (beta/2.0) * gamtot1 * gamtot1

       where (abs(denominator) < epsilon(0.0)) ! it == ik == 1 only
          bpar = 0.0
       elsewhere
          bpar = numerator / denominator
       end where
    end if

  end subroutine get_init_field

  subroutine flux (phi, apar, bpar, &
       pflux,  qflux,  vflux, vflux_par, vflux_perp, &
       pmflux, qmflux, vmflux, &
       pbflux, qbflux, vbflux)

!CMR, 15/1/08: 
!  Implemented Clemente Angioni's fix for fluxes by replacing g with gnew 
!  so fields and distribution function are evaluated self-consistently in time.
!  This fixed unphysical oscillations in non-ambipolar particle fluxes 
!
    use species, only: spec
    use theta_grid, only: ntgrid, bmag, gradpar, grho, delthet, drhodpsi
    use theta_grid, only: qval, shat, gds21, gds22
    use kt_grids, only: naky, ntheta0, akx, theta0, aky
    use le_grids, only: energy
    use dist_fn_arrays, only: gnew, aj0, vpac, vpa, aj1, vperp2
    use gs2_layouts, only: g_lo, ie_idx, is_idx, it_idx, ik_idx
    use mp, only: proc0
    use run_parameters, only: woutunits, fphi, fapar, fbpar
    use constants, only: zi
    use geometry, only: rhoc
    use theta_grid, only: Rplot, Bpol
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    real, dimension (:,:,:), intent (out) :: pflux, pmflux, pbflux
    real, dimension (:,:,:), intent (out) :: vflux, vmflux, vbflux, vflux_par, vflux_perp
    real, dimension (:,:,:,:), intent (out) :: qflux, qmflux, qbflux
    real, dimension (:,:,:), allocatable :: dnorm
    integer :: it, ik, is, isgn, ig
    integer :: iglo

    allocate (dnorm (-ntgrid:ntgrid,ntheta0,naky))

    if (proc0) then
       pflux = 0.0;   qflux = 0.0;   vflux = 0.0 ; vflux_par = 0.0 ; vflux_perp = 0.0
       pmflux = 0.0;  qmflux = 0.0;  vmflux = 0.0
       pbflux = 0.0;  qbflux = 0.0;  vbflux = 0.0
    end if

    do ik = 1, naky
       do it = 1, ntheta0
          dnorm(:,it,ik) = delthet/bmag/gradpar*woutunits(ik)
       end do
    end do

    if (fphi > epsilon(0.0)) then
       do isgn = 1, 2
          g0(:,isgn,:) = gnew(:,isgn,:)*aj0
       end do
       call get_flux (phi, pflux, dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          g0(:,:,iglo) = g0(:,:,iglo)*energy(ie_idx(g_lo,iglo),is_idx(g_lo,iglo))
       end do
       call get_flux (phi, qflux(:,:,:,1), dnorm)

       do isgn = 1, 2
          g0(:,isgn,:) = gnew(:,isgn,:)*2.*vpa(:,isgn,:)**2*aj0
       end do
       call get_flux (phi, qflux(:,:,:,2), dnorm)

       do isgn = 1, 2
          g0(:,isgn,:) = gnew(:,isgn,:)*vperp2*aj0
       end do
       call get_flux (phi, qflux(:,:,:,3), dnorm)

       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             g0(ig,isgn,:) = gnew(ig,isgn,:)*aj0(ig,:)*vpac(ig,isgn,:)*Rplot(ig)*sqrt(1.0-Bpol(ig)**2/bmag(ig)**2)
          end do
       end do
       call get_flux (phi, vflux_par, dnorm)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) = -zi*aky(ik)*gnew(:,isgn,iglo)*aj1(:,iglo) &
                  *rhoc*(gds21+theta0(it,ik)*gds22)*vperp2(:,iglo)*spec(is)%smz/(qval*shat*bmag**2)
!             g0(:,isgn,iglo) = zi*akx(it)*grho*gnew(:,isgn,iglo)*aj1(:,iglo) &
!                  *2.0*vperp2(:,iglo)*spec(is)%smz/(bmag**2*drhodpsi)
          end do
       end do
       call get_flux (phi, vflux_perp, dnorm)
       vflux = vflux_par + vflux_perp

    else
       pflux = 0.
       qflux = 0.
!       qflux_par = 0.
!       qflux_perp = 0.
       vflux = 0.
    end if

    if (fapar > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = -gnew(:,isgn,iglo)*aj0(:,iglo)*spec(is)%stm*vpa(:,isgn,iglo)
          end do
       end do
       call get_flux (apar, pmflux, dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          g0(:,:,iglo) = g0(:,:,iglo)*energy(ie_idx(g_lo,iglo),is_idx(g_lo,iglo))
       end do
       call get_flux (apar, qmflux(:,:,:,1), dnorm)
       
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = -gnew(:,isgn,iglo)*aj0(:,iglo)*spec(is)%stm*vpa(:,isgn,iglo) &
                  *2.*vpa(:,isgn,iglo)**2
          end do
       end do
       call get_flux (apar, qmflux(:,:,:,2), dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = -gnew(:,isgn,iglo)*aj0(:,iglo)*spec(is)%stm*vpa(:,isgn,iglo) &
                  *vperp2(:,iglo)
          end do
       end do
       call get_flux (apar, qmflux(:,:,:,3), dnorm)
       
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = -gnew(:,isgn,iglo)*aj0(:,iglo)*spec(is)%stm &
                  *vpa(:,isgn,iglo)*vpac(:,isgn,iglo)
          end do
       end do
       call get_flux (apar, vmflux, dnorm)
    else
       pmflux = 0.
       qmflux = 0.
!       qmflux_par = 0.
!       qmflux_perp = 0.
       vmflux = 0.
    end if

    if (fbpar > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = gnew(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo)*spec(is)%tz
          end do
       end do
       call get_flux (bpar, pbflux, dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          g0(:,:,iglo) = g0(:,:,iglo)*energy(ie_idx(g_lo,iglo),is_idx(g_lo,iglo))
       end do
       call get_flux (bpar, qbflux(:,:,:,1), dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = gnew(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo)*spec(is)%tz &
                    *2.*vpa(:,isgn,iglo)**2
          end do
       end do
       call get_flux (bpar, qbflux(:,:,:,2), dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = gnew(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo)*spec(is)%tz &
                    *vperp2(:,iglo)
          end do
       end do
       call get_flux (bpar, qbflux(:,:,:,3), dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = gnew(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo) &
                  *spec(is)%tz*vpac(:,isgn,iglo)
          end do
       end do
       call get_flux (bpar, vbflux, dnorm)
    else
       pbflux = 0.
       qbflux = 0.
!       qbflux_par = 0.
!       qbflux_perp = 0.
       vbflux = 0.
    end if

    deallocate (dnorm)
  end subroutine flux

  subroutine flux_emu (phi, apar, bpar, pflux, pmflux, pbflux)
    use species, only: spec
    use theta_grid, only: ntgrid, bmag, gradpar, grho, delthet, drhodpsi
    use theta_grid, only: qval, shat, gds21, gds22
    use kt_grids, only: naky, ntheta0, akx, theta0, aky
    use le_grids, only: energy
    use dist_fn_arrays, only: gnew, aj0, vpac, vpa, aj1, vperp2
    use gs2_layouts, only: g_lo, ie_idx, is_idx, it_idx, ik_idx
    use mp, only: proc0
    use run_parameters, only: woutunits, fphi, fapar, fbpar
    use constants, only: zi
    use geometry, only: rhoc
    use theta_grid, only: Rplot, Bpol
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    real, dimension (:,:,:,:), intent (inout) :: pflux, pmflux, pbflux
    real, dimension (:,:,:), allocatable :: dnorm
    integer :: it, ik, is, isgn, ig
    integer :: iglo

    allocate (dnorm (-ntgrid:ntgrid,ntheta0,naky))

    if (proc0) then
        pflux = 0.0;   pmflux = 0.0; pbflux = 0.0 ; 
    end if

    do ik = 1, naky
       do it = 1, ntheta0
          dnorm(:,it,ik) = delthet/bmag/gradpar*woutunits(ik)
       end do
    end do

    if (fphi > epsilon(0.0)) then
       do isgn = 1, 2
          g0(:,isgn,:) = gnew(:,isgn,:)*aj0
       end do
       call get_flux_emu (phi, pflux, dnorm)

    else
       pflux = 0.
    end if

    if (fapar > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = -gnew(:,isgn,iglo)*aj0(:,iglo)*spec(is)%stm*vpa(:,isgn,iglo)
          end do
       end do
       call get_flux_emu (apar, pmflux, dnorm)

    else
       pmflux = 0.
    end if

    if (fbpar > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = gnew(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo)*spec(is)%tz
          end do
       end do
       call get_flux_emu (bpar, pbflux, dnorm)
    else
       pbflux = 0.
    end if

    deallocate (dnorm)
  end subroutine flux_emu

  subroutine flux_e (phi, apar, bpar, pflux, pmflux, pbflux)
    use species, only: spec, nspec
    use theta_grid, only: ntgrid, bmag, gradpar, grho, delthet, drhodpsi
    use theta_grid, only: qval, shat, gds21, gds22
    use kt_grids, only: naky, ntheta0, akx, theta0, aky
    use le_grids, only: energy
    use dist_fn_arrays, only: gnew, aj0, vpac, vpa, aj1, vperp2
    use gs2_layouts, only: g_lo, ie_idx, is_idx, it_idx, ik_idx
    use mp, only: proc0
    use run_parameters, only: woutunits, fphi, fapar, fbpar
    use constants, only: zi
    use geometry, only: rhoc
    use theta_grid, only: Rplot, Bpol
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    real, dimension (:,:), intent (inout) :: pflux, pmflux, pbflux
    real, dimension (:,:,:), allocatable :: dnorm
    integer :: it, ik, is, isgn, ig
    integer :: iglo

    allocate (dnorm (-ntgrid:ntgrid,ntheta0,naky))

    if (proc0) then
        pflux = 0.0;   pmflux = 0.0; pbflux = 0.0 ; 
    end if

    do ik = 1, naky
       do it = 1, ntheta0
          dnorm(:,it,ik) = delthet/bmag/gradpar*woutunits(ik)
       end do
    end do

    if (fphi > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) = aimag(gnew(:,isgn,iglo)*aj0(:,iglo)*aky(ik)*conjg(phi(:,it,ik)))
          end do
       end do
       call get_flux_e (phi, pflux, dnorm)
!       do is = 1,nspec
!          call fsavg_keep_e(g0,pflux(:,is))
!       end do
      

    else
       pflux = 0.
    end if

    if (fapar > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = -gnew(:,isgn,iglo)*aj0(:,iglo)*spec(is)%stm*vpa(:,isgn,iglo)
          end do
       end do
       call get_flux_e (apar, pmflux, dnorm)

    else
       pmflux = 0.
    end if

    if (fbpar > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = gnew(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo)*spec(is)%tz
          end do
       end do
       call get_flux_e (bpar, pbflux, dnorm)
    else
       pbflux = 0.
    end if

    deallocate (dnorm)
  end subroutine flux_e


  subroutine get_flux (fld, flx, dnorm)
    use theta_grid, only: ntgrid, delthet, grho
    use kt_grids, only: ntheta0, aky, naky
    use le_grids, only: integrate_moment
    use species, only: nspec
    use mp, only: proc0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (:,:,:), intent (in out) :: flx
    real, dimension (-ntgrid:,:,:) :: dnorm
    complex, dimension (:,:,:,:), allocatable :: total
    real :: wgt
    integer :: ik, it, is, ig

    allocate (total(-ntgrid:ntgrid,ntheta0,naky,nspec))
    call integrate_moment (g0, total)

    if (proc0) then
       do is = 1, nspec
          do ik = 1, naky
             do it = 1, ntheta0
                wgt = sum(dnorm(:,it,ik)*grho)
                flx(it,ik,is) = sum(aimag(total(:,it,ik,is)*conjg(fld(:,it,ik))) &
                     *dnorm(:,it,ik)*aky(ik))/wgt
             end do
          end do
       end do
       flx = flx*0.5
    end if

    deallocate (total)

  end subroutine get_flux

  subroutine fsavg_keep_e (g,flx)
  !< Integrates a function g over lambda, flux-surface-averages, returning to result to flx on proc0
  !! Similar to get_flux_e
    use theta_grid, only: ntgrid, delthet, grho, bmag, gradpar, jacob
    use run_parameters, only: woutunits
    use kt_grids, only: ntheta0, aky, naky
    use le_grids, only: negrid,nlambda,wl
    use species, only: nspec
    use mp, only: sum_reduce, proc0
    use gs2_layouts, only: g_lo,ie_idx,il_idx,is_idx,it_idx,ik_idx,isign_idx,yxf_lo
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:), intent (inout) :: flx
    real, dimension (-ntgrid:ntgrid,ntheta0,naky) :: dnorm
    complex, dimension (:,:), allocatable :: total
    real, dimension (-ntgrid:ntgrid) :: wgt
    real:: fac
    integer :: ik, it, is, ig, iglo, ie, il, isgn

    allocate(total(negrid,nspec))
    total = 0.
    
    wgt(-ntgrid:) = 0.0
    do ig=-ntgrid,ntgrid-1
       wgt(ig) = delthet(ig)*(jacob(ig)+jacob(ig+1))*0.5
    end do
    wgt(-ntgrid:) = wgt(-ntgrid:) / sum(wgt)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       fac =0.5
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       if (aky(ik) == 0.) fac = 1.0

       total(ie,is) = total(ie,is) + fac*sum((g0(-ntgrid:,1,iglo) + g0(-ntgrid:,2,iglo)) &
                  * dnorm(-ntgrid:,it,ik) * wl(-ntgrid:,il)*wgt(-ntgrid:))

    end do

    call sum_reduce(total,0)
    if (proc0) flx = total

    deallocate(total)

  end subroutine fsavg_keep_e


  subroutine get_flux_e (fld, flx, dnorm)
    use theta_grid, only: ntgrid, delthet, grho
    use kt_grids, only: ntheta0, aky, naky
    use le_grids, only: negrid,nlambda,wl
    use species, only: nspec
    use mp, only: sum_reduce, proc0
    use gs2_layouts, only: g_lo,ie_idx,il_idx,is_idx,it_idx,ik_idx,isign_idx,yxf_lo
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (:,:), intent (inout) :: flx
    real, dimension (-ntgrid:,:,:) :: dnorm
    real, dimension (:,:), allocatable :: total
    real:: wgt,fac
    integer :: ik, it, is, ig, iglo, ie, il, isgn

    allocate(total(negrid,nspec))
    total = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       fac =0.5
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       if (aky(ik) == 0.) fac = 1.0

       wgt = sum(dnorm(-ntgrid:,it,ik)*grho(-ntgrid:))
       total(ie,is) = total(ie,is) + fac*& 
          sum(aimag(g0(-ntgrid:,1,iglo)*conjg(fld(-ntgrid:,it,ik)) * &
          dnorm(-ntgrid:,it,ik) * aky(ik)) * wl(-ntgrid:,il)/wgt)
       total(ie,is) = total(ie,is) + fac*& 
          sum(aimag(g0(-ntgrid:,2,iglo)*conjg(fld(-ntgrid:,it,ik)) * &
          dnorm(-ntgrid:,it,ik) * aky(ik)) * wl(-ntgrid:,il)/wgt)
    end do

    call sum_reduce(total,0)

    if (proc0) flx = 0.5*total  

    deallocate(total)

  end subroutine get_flux_e

  subroutine get_flux_emu (fld, flx, dnorm)
    use theta_grid, only: ntgrid, delthet, grho
    use kt_grids, only: ntheta0, aky, naky
    use le_grids, only: negrid,nlambda,Bovervpar
    use species, only: nspec
    use mp, only: sum_reduce, proc0
    use gs2_layouts, only: g_lo,ie_idx,il_idx,is_idx,it_idx,ik_idx,isign_idx,yxf_lo
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (:,:,:,:), intent (in out) :: flx
    real, dimension (-ntgrid:,:,:) :: dnorm
    real, dimension (:,:,:,:), allocatable :: total
    real:: wgt,fac
    integer :: ik, it, is, ig, iglo, ie, il, isgn

    allocate(total(negrid,nlambda,2,nspec))
    total = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       fac = 0.5
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       if (aky(ik) == 0.) fac = 1.0

       wgt = sum(dnorm(-ntgrid:,it,ik)*grho(-ntgrid:))
       total(ie,il,1,is) = total(ie,il,1,is) + fac*& 
          sum(aimag(g0(-ntgrid:,1,iglo)*conjg(fld(-ntgrid:,it,ik)) * &
          dnorm(-ntgrid:,it,ik) * aky(ik)) * Bovervpar(-ntgrid:,il)/wgt)
       total(ie,il,2,is) = total(ie,il,2,is) + fac*& 
          sum(aimag(g0(-ntgrid:,2,iglo)*conjg(fld(-ntgrid:,it,ik)) * &
          dnorm(-ntgrid:,it,ik) * aky(ik)) * Bovervpar(-ntgrid:,il)/wgt)
    end do

    call sum_reduce(total,0)

    if (proc0) flx = 0.5*total    !< Reality condition from summing over ky

    deallocate(total)

  end subroutine get_flux_emu



  subroutine eexchange (phi, exchange)

    use mp, only: proc0
    use constants, only: zi
    use gs2_layouts, only: g_lo, il_idx, ie_idx, it_idx, ik_idx, is_idx
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: gnew, aj0, vpac
    use theta_grid, only: ntgrid, gradpar, delthet, bmag, jacob
    use kt_grids, only: ntheta0, naky
    use le_grids, only: integrate_moment
    use run_parameters, only: woutunits, fphi
    use species, only: spec, nspec
    use nonlinear_terms, only: nonlin

    implicit none

    complex, dimension (-ntgrid:,:,:), intent (in) :: phi
    real, dimension (:,:,:), intent (out) :: exchange

    integer :: ig, il, ie, it, ik, is, iglo, isgn
    real :: wgt
    real, dimension (:,:,:), allocatable :: dnorm
    complex, dimension (:,:,:,:), allocatable :: total

    allocate (dnorm(-ntgrid:ntgrid, ntheta0, naky)) ; dnorm = 0.0
    allocate (total(-ntgrid:ntgrid, ntheta0, naky, nspec)) ; total = 0.0

    if (proc0) exchange = 0.0

    do ik = 1, naky
       do it = 1, ntheta0
          dnorm(:ntgrid-1,it,ik) = delthet(:ntgrid-1)*(jacob(-ntgrid+1:)+jacob(:ntgrid-1))*0.5
       end do
    end do

    if (fphi > epsilon(0.0)) then
       g0 = 0.
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)          
	  if (nonlin .and. it==1 .and. ik==1) cycle
          do isgn = 1, 2
             ! get v_magnetic piece of g0 at grid points instead of cell centers
             do ig = -ntgrid, ntgrid
                g0(ig,isgn,iglo) = aj0(ig,iglo)*(zi*wdrift_func(ig,il,ie,it,ik,is)/code_dt &
                     * gnew(ig,isgn,iglo)*spec(is)%tz)
             end do
             ! get v_magnetic piece of g0 at cell centers and add in vpar piece at cell centers
             do ig = -ntgrid, ntgrid-1
                g0(ig,isgn,iglo) = 0.5*(g0(ig,isgn,iglo)+g0(ig+1,isgn,iglo)) &
                     + 0.5*vpac(ig,isgn,iglo)*(gradpar(ig)+gradpar(ig+1))/delthet(ig) &
                     * (gnew(ig+1,isgn,iglo)-gnew(ig,isgn,iglo))*spec(is)%stm
             end do

!              do ig = -ntgrid+1, ntgrid-1
!                 g0(ig,isgn,iglo) = aj0(ig,iglo)*(zi*wdrift_func(ig, il, ie, it, ik)/code_dt &
!                      * gnew(ig,isgn,iglo)*spec(is)%tz + vpa(ig,isgn,iglo)*gradpar(ig)/(delthet(ig)+delthet(ig+1)) &
!                      * (gnew(ig+1,isgn,iglo)-gnew(ig-1,isgn,iglo))*spec(is)%stm)
!              end do
!              if (isgn == 1) then
!                 ! g = 0 at ig = -ntgrid for vpa > 0
!                 ig = -ntgrid
!                 g0(ig,isgn,iglo) = aj0(ig,iglo)*(zi*wdrift_func(ig, il, ie, it, ik)/code_dt &
!                      * gnew(ig,isgn,iglo)*spec(is)%tz + vpa(ig,isgn,iglo)*gradpar(ig)/delthet(ig) &
!                      * gnew(ig+1,isgn,iglo)*spec(is)%stm)
!                 ig = ntgrid
!                 g0(ig,isgn,iglo) = aj0(ig,iglo)*(zi*wdrift_func(ig, il, ie, it, ik)/code_dt &
!                      * gnew(ig,isgn,iglo)*spec(is)%tz + vpa(ig,isgn,iglo)*gradpar(ig)/delthet(ig-1) &
!                      * (gnew(ig,isgn,iglo)-gnew(ig-1,isgn,iglo))*spec(is)%stm)
!              else
!                 ! g = 0 at ig = ntgrid for vpa < 0
!                 ig = ntgrid
!                 g0(ig,isgn,iglo) = aj0(ig,iglo)*(zi*wdrift_func(ig, il, ie, it, ik)/code_dt &
!                      * gnew(ig,isgn,iglo)*spec(is)%tz - vpa(ig,isgn,iglo)*gradpar(ig)/delthet(ig-1) &
!                      * gnew(ig-1,isgn,iglo)*spec(is)%stm)
!                 ig = -ntgrid
!                 g0(ig,isgn,iglo) = aj0(ig,iglo)*(zi*wdrift_func(ig, il, ie, it, ik)/code_dt &
!                      * gnew(ig,isgn,iglo)*spec(is)%tz + vpa(ig,isgn,iglo)*gradpar(ig)/delthet(ig) &
!                      * (gnew(ig+1,isgn,iglo)-gnew(ig,isgn,iglo))*spec(is)%stm)
!              end if
          end do
       end do

       call integrate_moment (g0, total)

       if (proc0) then
          do is = 1, nspec
             do ik = 1, naky
                do it = 1, ntheta0
                   if (nonlin .and. it==1 .and. ik==1) cycle
                   wgt = sum(dnorm(:,it,ik))
                   exchange(it,ik,is) = sum(real(total(:,it,ik,is)*conjg(phi(:,it,ik))) &
                        *dnorm(:,it,ik))/wgt
                end do
             end do
          end do
       end if

    end if

    deallocate (dnorm, total)

  end subroutine eexchange

  subroutine lf_flux (phi, vflx0, vflx1)

    use species, only: spec, nspec
    use theta_grid, only: ntgrid, bmag, gradpar, grho, delthet
    use theta_grid, only: qval, shat, gds21, gds22, drhodpsi, IoB
    use kt_grids, only: naky, ntheta0, akx, aky
    use dist_fn_arrays, only: gnew, aj0, vpac, vpa, aj1, vperp2
    use gs2_layouts, only: g_lo, ie_idx, is_idx, it_idx, ik_idx
    use mp, only: proc0
    use run_parameters, only: woutunits, fphi, fapar, fbpar, rhostar
    use constants, only: zi
    use geometry, only: rhoc
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi
    real, dimension (:,:,:), intent (out) :: vflx0, vflx1
    real, dimension (:,:), allocatable :: dum
    real, dimension (:,:,:), allocatable :: dnorm
    complex, dimension (:,:,:), allocatable :: dphi
    integer :: it, ik, is, isgn, ig
    integer :: iglo

    allocate (dnorm (-ntgrid:ntgrid,ntheta0,naky))
    allocate (dum (-ntgrid:ntgrid,nspec))
    allocate (dphi (-ntgrid:ntgrid,ntheta0,naky))

    if (proc0) then
       vflx0 = 0.0 ; vflx1 = 0.0 ; dum = 0.0
    end if

    do ik = 1, naky
       do it = 1, ntheta0
          dnorm(:,it,ik) = delthet/bmag/gradpar*woutunits(ik)
       end do
    end do

    do ig = -ntgrid, ntgrid-1
       dphi(ig,:,:) = (phi(ig+1,:,:)-phi(ig,:,:))/delthet(ig)
    end do
    ! not sure if this is the right way to handle ntgrid point -- MAB
    dphi(ntgrid,:,:) = (phi(ntgrid,:,:)-phi(ntgrid-1,:,:))/delthet(-ntgrid)

    if (fphi > epsilon(0.0)) then
       ! this is the second term in Pi_0^{tb} in toroidal_flux.pdf notes
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             g0(:,isgn,iglo) = -zi*gnew(:,isgn,iglo)*aj0(:,iglo)*vpa(:,isgn,iglo) &
                  *drhodpsi*IoB**2*gradpar*rhostar
          end do
       end do
       call get_lfflux (dphi, vflx0, dnorm)

       ! this is the bracketed part of the first term in Pi_0^{tb} in toroidal_flux.pdf notes
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             g0(:,isgn,iglo) = 0.5*gnew(:,isgn,iglo)*aj0(:,iglo)*vpa(:,isgn,iglo)**2 &
                  *drhodpsi*IoB**2*rhostar
          end do
       end do
       call get_flux (phi, vflx1, dnorm)

!        do isgn = 1, 2
!           do ig = -ntgrid, ntgrid
!              if (allocated(rmajor_geo)) then
!                 g0(ig,isgn,:) = gnew(ig,isgn,:)*aj0(ig,:)*vpac(ig,isgn,:)*rmajor_geo(ig)*sqrt(1.0-bpol_geo(ig)**2/bmag(ig)**2)
!              else
!                 g0(ig,isgn,:) = gnew(ig,isgn,:)*aj0(ig,:)*vpac(ig,isgn,:)
!              end if
!           end do
!        end do
!        call get_flux (phi, vflux_par, dnorm)
    else
       vflx0 = 0. ; vflx1 = 0.
    end if

    deallocate (dnorm,dum,dphi)
  end subroutine lf_flux

  subroutine get_lfflux (fld, flx, dnorm)
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, aky, naky
    use le_grids, only: integrate_moment
    use theta_grid, only: grho
    use species, only: nspec
    use mp, only: proc0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (:,:,:), intent (in out) :: flx
    real, dimension (-ntgrid:,:,:) :: dnorm
    complex, dimension (:,:,:,:), allocatable :: total
    real :: wgt
    integer :: ik, it, is, ig

    allocate (total(-ntgrid:ntgrid,ntheta0,naky,nspec))
    call integrate_moment (g0, total)

    if (proc0) then
       do is = 1, nspec
          do ik = 1, naky
             do it = 1, ntheta0
                wgt = sum(dnorm(:,it,ik)*grho)
                flx(it,ik,is) = sum(aimag(total(:,it,ik,is)*conjg(fld(:,it,ik))) &
                     *dnorm(:,it,ik))/wgt
             end do
          end do
       end do

       flx = flx*0.5

    end if

    deallocate (total)

  end subroutine get_lfflux

  subroutine flux_vs_theta_vs_vpa (vflx)

    use constants, only: zi
    use dist_fn_arrays, only: gnew, vperp2, aj1, aj0, vpac
    use fields_arrays, only: phinew
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: it_idx, ik_idx, is_idx
    use geometry, only: rhoc
    use theta_grid, only: ntgrid, bmag, gds21, gds22, qval, shat
    use theta_grid, only: Rplot, Bpol
    use kt_grids, only: aky, theta0
    use le_grids, only: integrate_volume, nlambda, negrid
    use le_grids, only: get_flux_vs_theta_vs_vpa
    use species, only: spec, nspec

    implicit none

    real, dimension (-ntgrid:,:,:), intent (out) :: vflx
    
    integer :: all = 1
    integer :: iglo, isgn, ig, it, ik, is

    real, dimension (:,:,:), allocatable :: g0r
    real, dimension (:,:,:,:,:), allocatable :: gavg

    allocate (g0r(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (gavg(-ntgrid:ntgrid,nlambda,negrid,2,nspec))

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             g0(ig,isgn,iglo) = gnew(ig,isgn,iglo)*aj0(ig,iglo)*vpac(ig,isgn,iglo) &
                  *Rplot(ig)*sqrt(1.0-Bpol(ig)**2/bmag(ig)**2) &
                  -zi*aky(ik)*gnew(ig,isgn,iglo)*aj1(ig,iglo) &
                  *rhoc*(gds21(ig)+theta0(it,ik)*gds22(ig))*vperp2(ig,iglo)*spec(is)%smz/(qval*shat*bmag(ig)**2)
             g0r(ig,isgn,iglo) = aimag(g0(ig,isgn,iglo)*conjg(phinew(ig,it,ik)))*aky(ik)
          end do
       end do
    end do
    
    call integrate_volume (g0r, gavg, all)
    call get_flux_vs_theta_vs_vpa (gavg, vflx)

    deallocate (gavg)
    deallocate (g0r)

  end subroutine flux_vs_theta_vs_vpa

!=============================================================================
! Density: Calculate Density perturbations
!=============================================================================
  subroutine get_jext(j_ext)
    use dist_fn_arrays, only: kperp2
    use kt_grids, only: ntheta0, naky,aky
    use mp, only: proc0
    use theta_grid, only: jacob, delthet, ntgrid
    use antenna, only: antenna_apar
    implicit none
    !Passed
    real, dimension(:,:) ::  j_ext
    !Local 
    complex, dimension(:,:,:), allocatable :: j_extz
    integer :: ig,ik, it                             !Indices
    real :: fac2                                     !Factor
    real, dimension (:), allocatable :: wgt
    

    !Get z-centered j_ext at current time
    allocate (j_extz(-ntgrid:ntgrid, ntheta0, naky)) ; j_extz = 0.
    call antenna_apar (kperp2, j_extz)       

    !Set weighting factor for z-averages
!    if (proc0) then
       allocate (wgt(-ntgrid:ntgrid))
       !GGH NOTE: Here wgt is 1/(2*ntgrid)
       wgt = 0.
       do ig=-ntgrid,ntgrid-1
          wgt(ig) = delthet(ig)*jacob(ig)
       end do
       wgt = wgt/sum(wgt)         
!    endif

    !Take average over z
    do ig=-ntgrid, ntgrid-1
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, ntheta0
             j_ext(it,ik)=j_ext(it,ik)+real(j_extz(ig,it,ik))*wgt(ig)*fac2
          end do
       end do
    enddo

    if (proc0) deallocate (wgt)
    deallocate (j_extz)

  end subroutine get_jext
!<GGH
!=============================================================================
  subroutine get_heat (h, hk, he, phi, apar, bpar, phinew, aparnew, bparnew)
    use mp, only: proc0, iproc
    use constants, only: pi, zi
    use kt_grids, only: ntheta0, naky, aky, akx
#ifdef LOWFLOW
    use dist_fn_arrays, only: hneoc
#endif
    use dist_fn_arrays, only: vpa, vpac, aj0, aj1, vperp2, g, gnew, kperp2, g_adjust
    use gs2_heating, only: heating_diagnostics
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx, ie_idx
    use le_grids, only: integrate_moment, negrid
    use general_f0, only: df0dE
    use species, only: spec, nspec,has_electron_species
    use theta_grid, only: jacob, delthet, ntgrid
    use run_parameters, only: fphi, fapar, fbpar, tunits, beta, tite
    use gs2_time, only: code_dt
    use nonlinear_terms, only: nonlin
    use antenna, only: antenna_apar, a_ext_data
    use hyper, only: D_v, D_eta, nexp, hypervisc_filter

    implicit none
    type (heating_diagnostics) :: h
    type (heating_diagnostics), dimension(:,:) :: hk
    type (heating_diagnostics), dimension(:) :: he
!    complex, dimension (-ntgrid:,:,:), pointer :: hh, hnew
    complex, dimension (-ntgrid:,:,:) :: phi, apar, bpar, phinew, aparnew, bparnew
    complex, dimension(:,:,:,:), allocatable :: tot
    complex, dimension(:,:), allocatable :: tote
!    complex, dimension (:,:,:), allocatable :: epar
    complex, dimension(:,:,:), allocatable :: bpardot, apardot, phidot, j_ext
    complex :: chi, havg
    complex :: chidot, j0phiavg, j1bparavg, j0aparavg
!    complex :: pstar, pstardot, gdot
    complex :: phi_m, apar_m, bpar_m, hdot
    complex :: phi_avg, bpar_avg, bperp_m, bperp_avg
    complex :: de, denew
    complex :: dgdt_hypervisc
    real, dimension (:), allocatable :: wgt
    real :: fac2, dtinv, akperp4
    integer :: isgn, iglo, ig, is, ik, it, ie

    g0(ntgrid,:,:) = 0.

! ==========================================================================
! Ion/Electron heating------------------------------------------------------
! ==========================================================================

    allocate ( phidot(-ntgrid:ntgrid, ntheta0, naky))
    allocate (apardot(-ntgrid:ntgrid, ntheta0, naky))
    allocate (bpardot(-ntgrid:ntgrid, ntheta0, naky))

    call dot ( phi,  phinew,  phidot, fphi)
    call dot (apar, aparnew, apardot, fapar)
    call dot (bpar, bparnew, bpardot, fbpar)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Next two calls make the variables g, gnew = h, hnew 
!!! until the end of this procedure!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call g_adjust (g,    phi,    bpar,    fphi, fbpar)
    call g_adjust (gnew, phinew, bparnew, fphi, fbpar)

    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       ie = ie_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       dtinv = 1./(code_dt*tunits(ik))
       do isgn=1,2
          
          do ig=-ntgrid, ntgrid-1
             
             chidot = aj0(ig,iglo)*(phidot(ig,it,ik) &
                  - vpac(ig,isgn,iglo) * spec(is)%stm * apardot(ig,it,ik)) &
                  + aj1(ig,iglo)*2.0*vperp2(ig,iglo)*bpardot(ig,it,ik)*spec(is)%tz
             
             hdot = fdot (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo), dtinv)
             
             havg = favg (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo))
             
! First term on RHS and LHS of Eq B-10 of H1:

             g0(ig,isgn,iglo) = spec(is)%dens*conjg(havg)* &
                  (chidot*spec(is)%z+hdot*(1.0/df0dE(ie,is)))

          end do
       end do
    end do

    deallocate (phidot, apardot, bpardot)

    allocate (tot(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (tote(negrid, nspec))

    call integrate_moment (g0, tot)
    call fsavg_keep_e (g0, tote)

    if (proc0) then
       allocate (wgt(-ntgrid:ntgrid))
       wgt = 0.
       do ig=-ntgrid,ntgrid-1
!          wgt(ig) = delthet(ig)*jacob(ig)
! delthet is cell-centered, but jacob is given on grid
          wgt(ig) = delthet(ig)*(jacob(ig)+jacob(ig+1))*0.5
       end do
       wgt = wgt/sum(wgt)         

       do is = 1, nspec
          do ik = 1, naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             do it = 1, ntheta0
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                do ig = -ntgrid, ntgrid-1
                    hk(it,ik) % heating(is) = hk(it,ik) % heating(is) &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2 
                end do
                h % heating(is) = h % heating(is) + hk(it,ik) % heating(is)
             end do
          end do
          do ie = 1,negrid
             he(ie) % heating(is) = he(ie) % heating(is) &
                   + real(tote(ie,is))
          end do
       end do
    end if

! ==========================================================================
! Antenna Power and B-field contribution to E and E_dot---------------------
! ==========================================================================
    if (proc0) then
!       allocate (epar (-ntgrid:ntgrid, ntheta0, naky)) ; epar = 0.
       allocate (j_ext(-ntgrid:ntgrid, ntheta0, naky)) ; j_ext = 0.
       call antenna_apar (kperp2, j_ext)       
       
       if (beta > epsilon(0.)) then
          do ik=1,naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             dtinv = 1./(code_dt*tunits(ik))
             do it = 1,ntheta0
                
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                
                do ig=-ntgrid, ntgrid-1
                   
                   !GGH Time and space averaged estimate of d/dt(apar)
                   apar_m = fdot (apar   (ig  ,it,ik), &
                                  apar   (ig+1,it,ik), &
                                  aparnew(ig  ,it,ik), &
                                  aparnew(ig+1,it,ik), dtinv)*fapar

! J_ext.E when driving antenna only includes A_parallel:

                   hk(it,ik) % antenna = hk(it, ik) % antenna + real(conjg(j_ext(ig,it,ik))*apar_m)*wgt(ig)*fac2

                   !GGH Time and space averaged estimate of d/dt(bperp)
                   bperp_m = fdot (apar   (ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
                                   apar   (ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik)), &
                                   aparnew(ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
                                   aparnew(ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik)), dtinv) * fapar

                   !GGH Time and space averaged estimate of d/dt(bpar)
                   bpar_m = fdot (bpar   (ig  ,it,ik), &
                                  bpar   (ig+1,it,ik), &
                                  bparnew(ig  ,it,ik), &
                                  bparnew(ig+1,it,ik), dtinv)*fbpar

                   !GGH Time and space averaged estimate of bperp
                   bperp_avg = favg (apar   (ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
                                     apar   (ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik)), &
                                     aparnew(ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
                                     aparnew(ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik))) * fapar

                   !GGH Time and space averaged estimate of bpar
                   bpar_avg = favg (bpar   (ig  ,it,ik), &
                                    bpar   (ig+1,it,ik), &
                                    bparnew(ig  ,it,ik), &
                                    bparnew(ig+1,it,ik)) * fbpar

! 1/2 * d/dt B**2
!! GGH: Bug fixed on 2/06; error was in relative weight of B_par**2 and B_perp**2   
                   hk(it,ik) % energy_dot = hk(it,ik) % energy_dot + &
                        real(0.25 * conjg(bperp_m)*bperp_avg + conjg(bpar_m)*bpar_avg) &
                        * wgt(ig)*fac2*(2.0/beta)

! B**2/2
!! GGH: Bug fixed on 2/06; error was in relative weight of B_par**2 and B_perp**2   
                   hk(it,ik) % energy = hk(it,ik) % energy &
                        + 0.5*real((0.25*conjg(bperp_avg)*bperp_avg + conjg(bpar_avg)*bpar_avg)) &
                        * wgt(ig)*fac2*(2.0/beta)

                   !Eapar = int k_perp^2 A_par^2/(8 pi)                   
                   hk(it,ik) % eapar = hk(it,ik) % eapar &
                        + 0.5*real(0.25*conjg(bperp_avg)*bperp_avg) * wgt(ig)*fac2*(2.0/beta)

                   !Ebpar = int B_par^2/(8 pi)
                   hk(it,ik) % ebpar = hk(it,ik) % ebpar &
                        + 0.5*real(conjg(bpar_avg)*bpar_avg) * wgt(ig)*fac2*(2.0/beta)
                   
                end do
                h % antenna = h % antenna + hk(it, ik) % antenna
                h % eapar = h % eapar + hk(it, ik) % eapar
                h % ebpar = h % ebpar + hk(it, ik) % ebpar
             end do
          end do
       else
          hk % antenna = 0.
          he % antenna = 0.
          h  % antenna = 0.
          hk % energy_dot = 0.
          he % energy_dot = 0.
          hk % energy = 0.
          he % energy = 0.
          hk % eapar = 0.
          he % eapar = 0.
          h  % eapar = 0.
          hk % ebpar = 0.
          he % ebpar = 0.
          h  % ebpar = 0.
       end if
       deallocate (j_ext)
    end if

! ==========================================================================
! Finish E_dot--------------------------------------------------------------
! ==========================================================================

!GGH Include response of Boltzmann species for single-species runs

    if (.not. has_electron_species(spec)) then
       if (proc0) then
          !NOTE: It is assumed here that n0i=n0e and zi=-ze
          do ik=1,naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             dtinv = 1./(code_dt*tunits(ik))
             do it = 1,ntheta0

                if (nonlin .and. it == 1 .and. ik == 1) cycle

                do ig=-ntgrid, ntgrid-1

                   phi_avg = favg (phi   (ig  ,it,ik), &
                        phi   (ig+1,it,ik), &
                        phinew(ig  ,it,ik), &
                        phinew(ig+1,it,ik))

                   phi_m   = fdot (phi   (ig  ,it,ik), &
                        phi   (ig+1,it,ik), &
                        phinew(ig  ,it,ik), &
                        phinew(ig+1,it,ik), dtinv)

                   !NOTE: Adiabatic (Boltzmann) species has temperature
                   !       T = spec(1)%temp/tite
                   hk(it,ik) % energy_dot = hk(it,ik) % energy_dot + &
                        fphi * real(conjg(phi_avg)*phi_m) &
                        * spec(1)%dens * spec(1)%z * spec(1)%z * (tite/spec(1)%temp) &
                        * wgt(ig)*fac2

                end do
             end do
          end do
       endif
    endif !END Correction to E_dot for single species runs---------------------
 
!GGH New E_dot calc
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       dtinv = 1./(code_dt*tunits(ik))
       do isgn=1,2

          do ig=-ntgrid, ntgrid-1
             
             !Calculate old fluctuating energy de
             havg = favg_x (g(ig  ,isgn,iglo), &
                            g(ig+1,isgn,iglo))

             j0phiavg = favg_x (aj0(ig  ,iglo) * phi(ig  ,it,ik), &
                                aj0(ig+1,iglo) * phi(ig+1,it,ik)) * fphi * spec(is)%zt

             phi_avg = favg_x (phi(ig  ,it,ik), &
                               phi(ig+1,it,ik)) * fphi * spec(is)%zt

             de=0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg &
                  + conjg(phi_avg)*phi_avg &
                  - conjg(j0phiavg)*havg &
                  - conjg(havg)*j0phiavg)

            !Calculate new fluctuating energy denew
             havg = favg_x (gnew(ig  ,isgn,iglo), &
                            gnew(ig+1,isgn,iglo))

             j0phiavg = favg_x (aj0(ig  ,iglo) * phinew(ig  ,it,ik), &
                                aj0(ig+1,iglo) * phinew(ig+1,it,ik)) * fphi * spec(is)%zt

             phi_avg = favg_x (phinew(ig  ,it,ik), &
                               phinew(ig+1,it,ik)) * fphi * spec(is)%zt

             denew=0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg &
                  + conjg(phi_avg)*phi_avg &
                  - conjg(j0phiavg)*havg &
                  - conjg(havg)*j0phiavg) 

             !Set g0 as the change of energy (denew-de)/dt
             g0(ig,isgn,iglo) = fdot_t(de,denew,dtinv)

          end do
       end do
    end do
    !GGH -END new e_dot calc

   call integrate_moment (g0, tot)

    if (proc0) then
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, ntheta0
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do is = 1, nspec
                do ig = -ntgrid, ntgrid-1
                   hk(it,ik) % energy_dot = hk(it,ik) % energy_dot  &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                end do
             end do
             h % energy_dot = h % energy_dot + hk(it,ik) % energy_dot
          end do
       end do
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ==========================================================================
! Gradient Contributions to Heating-----------------------------------------
! ==========================================================================

    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       ie = ie_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle

       do isgn=1,2
          do ig=-ntgrid, ntgrid-1

             chi = favg (aj0(ig  ,iglo)*phi   (ig  ,it,ik),  &
                         aj0(ig+1,iglo)*phi   (ig+1,it,ik),  &
                         aj0(ig  ,iglo)*phinew(ig  ,it,ik),  &
                         aj0(ig+1,iglo)*phinew(ig+1,it,ik)) &
                         *fphi

!!GGH Bug fix: The apar part should be subtracted (because chi= phi - v|| A|| + B||)
             chi = chi - &
                  favg (aj0(ig  ,iglo)*apar   (ig  ,it,ik)*vpac(ig  ,isgn,iglo),  &
                        aj0(ig+1,iglo)*apar   (ig+1,it,ik)*vpac(ig+1,isgn,iglo),  &
                        aj0(ig  ,iglo)*aparnew(ig  ,it,ik)*vpac(ig  ,isgn,iglo),  &
                        aj0(ig+1,iglo)*aparnew(ig+1,it,ik)*vpac(ig+1,isgn,iglo)) &
                        *spec(is)%stm*fapar
                
             chi = chi + &
                  favg (aj1(ig  ,iglo)*2.0*bpar   (ig  ,it,ik)*vperp2(ig  ,iglo),  &
                        aj1(ig+1,iglo)*2.0*bpar   (ig+1,it,ik)*vperp2(ig+1,iglo),  &
                        aj1(ig  ,iglo)*2.0*bparnew(ig  ,it,ik)*vperp2(ig  ,iglo),  &
                        aj1(ig+1,iglo)*2.0*bparnew(ig+1,it,ik)*vperp2(ig+1,iglo)) &
                        *spec(is)%tz*fbpar

             havg = favg (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo))

#ifdef LOWFLOW
             g0(ig,isgn,iglo) = zi * wstar(ik,ie,is)*hneoc(ig,isgn,iglo)/code_dt * conjg(havg)*chi &
#else
             g0(ig,isgn,iglo) = zi * wstar(ik,ie,is)/code_dt * conjg(havg)*chi &
#endif
                  * spec(is)%dens
            
          end do
       end do
    end do

    call integrate_moment (g0, tot)
   call fsavg_keep_e (g0, tote)

    if (proc0) then
       do is = 1, nspec
          do ik = 1, naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             do it = 1, ntheta0
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                do ig = -ntgrid, ntgrid-1
                   hk(it,ik) % gradients(is) = hk(it,ik) % gradients(is) &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                end do
                h % gradients(is) = h % gradients(is) + hk(it,ik) % gradients(is)
             end do
          end do
          do ie = 1,negrid
             he(ie) % gradients(is) = he(ie) % gradients(is) &
                + real(tote(ie,is))
          end do
       end do
    end if
! ==========================================================================
! Hyperviscosity------------------------------------------------------------
! ==========================================================================

    if (D_v > epsilon(0.)) then

       do iglo=g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo, iglo)
          it = it_idx(g_lo, iglo)
          ik = ik_idx(g_lo, iglo)
          if (nonlin .and. it == 1 .and. ik == 1) cycle
!          akperp4 = (aky(ik)**2 + akx(it)**2)**nexp
           do isgn=1,2
             do ig=-ntgrid, ntgrid-1
                
                havg = favg (g   (ig  ,isgn,iglo), &
                             g   (ig+1,isgn,iglo), &
                             gnew(ig  ,isgn,iglo), &
                             gnew(ig+1,isgn,iglo)) 

                j0phiavg = favg (aj0(ig  ,iglo) * phi(ig  ,it,ik), &
                                 aj0(ig+1,iglo) * phi(ig+1,it,ik), &
                                 aj0(ig  ,iglo) * phinew(ig  ,it,ik), &
                                 aj0(ig+1,iglo) * phinew(ig+1,it,ik)) * fphi * spec(is)%zt

                j1bparavg= favg (aj1(ig  ,iglo)*2.0*bpar   (ig  ,it,ik)*vperp2(ig  ,iglo),  &
                                 aj1(ig+1,iglo)*2.0*bpar   (ig+1,it,ik)*vperp2(ig+1,iglo),  &
                                 aj1(ig  ,iglo)*2.0*bparnew(ig  ,it,ik)*vperp2(ig  ,iglo),  &
                                 aj1(ig+1,iglo)*2.0*bparnew(ig+1,it,ik)*vperp2(ig+1,iglo)) &
                                 *fbpar

                dgdt_hypervisc = 0.5*((1.0-1./hypervisc_filter(ig,it,ik))*gnew(ig,isgn,iglo) &
                     + (1.0-1./hypervisc_filter(ig+1,it,ik))*gnew(ig+1,isgn,iglo))/code_dt

!Set g0 for hyperviscous heating
!                g0(ig,isgn,iglo) = spec(is)%dens*spec(is)%temp*D_v*akperp4* &
!                     ( conjg(havg)*havg - conjg(havg)*j0phiavg - &
!                     conjg(havg)*j1bparavg)
                g0(ig,isgn,iglo) = spec(is)%dens*spec(is)%temp*conjg(havg)*dgdt_hypervisc

             end do
          end do
       end do

       call integrate_moment (g0, tot)
       if (proc0) then
          do ik = 1, naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             do it = 1, ntheta0
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                do is = 1, nspec
                   do ig = -ntgrid, ntgrid-1
                      hk(it,ik) % hypervisc(is) = hk(it,ik) % hypervisc(is) &
                           + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                   end do
                   h % hypervisc(is) = h % hypervisc(is) + hk(it,ik) % hypervisc(is)
                end do
             end do
          end do
       end if

    end if !End Hyperviscous Heating Calculation


! ==========================================================================
! Hyperresistivity------------------------------------------------------------
! ==========================================================================
 
    if (D_eta > epsilon(0.)) then

       do iglo=g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo, iglo)
          it = it_idx(g_lo, iglo)
          ik = ik_idx(g_lo, iglo)
          if (nonlin .and. it == 1 .and. ik == 1) cycle
          akperp4 = (aky(ik)**2 + akx(it)**2)**nexp
           do isgn=1,2
             do ig=-ntgrid, ntgrid-1
                
                havg = favg (g   (ig  ,isgn,iglo), &
                             g   (ig+1,isgn,iglo), &
                             gnew(ig  ,isgn,iglo), &
                             gnew(ig+1,isgn,iglo)) 

                j0aparavg = favg (aj0(ig  ,iglo) * apar(ig  ,it,ik), &
                                 aj0(ig+1,iglo)  * apar(ig+1,it,ik), &
                                 aj0(ig  ,iglo)  * aparnew(ig  ,it,ik), &
                                 aj0(ig+1,iglo)  * aparnew(ig+1,it,ik)) & 
                                 * fapar * spec(is)%zstm * vpac(ig,isgn,iglo)

!Set g0 for hyperresistive heating
                g0(ig,isgn,iglo) = spec(is)%dens*spec(is)%temp*D_eta*akperp4* &
                     conjg(havg)*j0aparavg

             end do
          end do
       end do

       call integrate_moment (g0, tot)
       if (proc0) then
          do ik = 1, naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             do it = 1, ntheta0
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                do is = 1, nspec
                   do ig = -ntgrid, ntgrid-1
                      hk(it,ik) % hyperres(is) = hk(it,ik) % hyperres(is) &
                           + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                   end do
                   h % hyperres(is) = h % hyperres(is) + hk(it,ik) % hyperres(is)
                end do
             end do
          end do
       end if

    end if !End Hyperresistivity Heating Calculation

!==========================================================================
!Finish Energy-------------------------------------------------------------
!==========================================================================

!GGH Calculate hs2-------------------------------------------------------------
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       ie = ie_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       do isgn=1,2

          do ig=-ntgrid, ntgrid-1
             
             havg = favg (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo))

             g0(ig,isgn,iglo) = -0.5*(1.0/df0dE(ie,is))*spec(is)%dens*(conjg(havg)*havg)
          end do
       end do
    end do

    call integrate_moment (g0, tot)
    call fsavg_keep_e (g0, tote)

    if (proc0) then
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, ntheta0
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do is = 1, nspec             
                do ig = -ntgrid, ntgrid-1

                   !hs2 = int_r int_v T/F0 hs^2/2
                   hk(it,ik) % hs2(is) = hk(it,ik) % hs2(is)  &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2

                end do
             end do
             h % hs2(:) = h % hs2(:) + hk(it,ik) % hs2(:)
          end do
       end do
       do ie = 1,negrid
          he(ie) % hs2(:) = he(ie) % hs2(:) + real(tote(ie,:))
       end do
    end if

!Calculate phis2-------------------------------------------------------------
    if (proc0) then
       do ik=1,naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1,ntheta0
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do ig=-ntgrid, ntgrid-1
                do is = 1, nspec
                   phi_avg = favg (phi   (ig  ,it,ik), &
                        phi   (ig+1,it,ik), &
                        phinew(ig  ,it,ik), &
                        phinew(ig+1,it,ik)) * fphi * spec(is)%zt

                   !hs2 = int_r int_v T/F0 hs^2/2
                   hk(it,ik) % phis2(is) = hk(it,ik) % phis2(is)  &
                        +0.5*spec(is)%temp*spec(is)%dens*real(conjg(phi_avg)*phi_avg) &
                        * wgt(ig) * fac2
                enddo
             end do
             h % phis2(:) = h % phis2(:) + hk(it,ik) % phis2(:)
          end do
       end do
    endif

! Calculate delfs2 (rest of energy)-----------------------------------------------

!GGH  Include response of Boltzmann species for single species runs
    if (.not. has_electron_species(spec)) then
       if (proc0) then
          !NOTE: It is assumed here that n0i=n0e and zi=-ze
          do ik=1,naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             dtinv = 1./(code_dt*tunits(ik))
             do it = 1,ntheta0

                if (nonlin .and. it == 1 .and. ik == 1) cycle

                do ig=-ntgrid, ntgrid-1

                   phi_avg = favg (phi   (ig  ,it,ik), &
                        phi   (ig+1,it,ik), &
                        phinew(ig  ,it,ik), &
                        phinew(ig+1,it,ik))

                   !NOTE: Adiabatic (Boltzmann) species has temperature
                   !       T = spec(1)%temp/tite
                   hk(it,ik) % energy = hk(it,ik) % energy + &
                        fphi * real(conjg(phi_avg)*phi_avg) &
                        * 0.5 * spec(1)%dens * spec(1)%z * spec(1)%z * (tite/spec(1)%temp) &
                        * wgt(ig)*fac2

                end do
             end do
          end do
       endif
    endif !END Correction to energy for single species runs---------------------

    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       ie = ie_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       do isgn=1,2

          do ig=-ntgrid, ntgrid-1
             
             havg = favg (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo))

             j0phiavg = favg (aj0(ig  ,iglo)*phi   (ig  ,it,ik), &
                              aj0(ig+1,iglo)*phi   (ig+1,it,ik), &  
                              aj0(ig  ,iglo)*phinew(ig  ,it,ik), &  
                              aj0(ig+1,iglo)*phinew(ig+1,it,ik)) * fphi * spec(is)%zt

             phi_avg = favg (phi   (ig  ,it,ik), &
                             phi   (ig+1,it,ik), &
                             phinew(ig  ,it,ik), &
                             phinew(ig+1,it,ik)) * fphi * spec(is)%zt

             g0(ig,isgn,iglo) = -0.5*(1.0/df0dE(ie,is))*spec(is)%dens*(conjg(havg)*havg &
                  + conjg(phi_avg)*phi_avg &
                  - conjg(j0phiavg)*havg &
                  - conjg(havg)*j0phiavg)
          end do
       end do
    end do

    call integrate_moment (g0, tot)
    call fsavg_keep_e (g0, tote)

    if (proc0) then
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, ntheta0
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do is = 1, nspec             
                do ig = -ntgrid, ntgrid-1
                   hk(it,ik) % energy = hk(it,ik) % energy  &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2

                   !Delfs2 = int_r int_v T/F0 dfs^2/2
                   hk(it,ik) % delfs2(is) = hk(it,ik) % delfs2(is)  &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                end do
             end do
             h % energy = h % energy + hk(it,ik) % energy
             h % delfs2(:) = h % delfs2(:) + hk(it,ik) % delfs2(:)
          end do
       end do
       do ie = 1,negrid
          he(ie) % delfs2(:) = he(ie) % delfs2(:) + real(tote(ie,:))
       end do
       deallocate (wgt)
    end if

    deallocate (tot)
    deallocate (tote)

!!
!! Put g, gnew back to their usual meanings
!!
    call g_adjust (g,    phi,    bpar,    -fphi, -fbpar)
    call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar)

  end subroutine get_heat

  subroutine reset_init

    use dist_fn_arrays, only: gnew, g
    initializing  = .true.
    initialized = .false.
    
    wdrift = 0.
    wdriftgen = 0.
    wdriftttp = 0.
    wdriftttpgen = 0.
    a = 0.
    b = 0.
    r = 0.
    ainv = 0.
    gnew = 0.
    g0 = 0.
    g = 0.

  end subroutine reset_init

  subroutine reset_physics

    call init_wstar

  end subroutine reset_physics

  subroutine get_verr (errest, erridx, phi, bpar)
     
    use mp, only: proc0, broadcast
    use le_grids, only: integrate_species
    use le_grids, only: eint_error, trap_error, lint_error, wdim
    use le_grids, only: ng2, nlambda, negrid, new_trap_int, jend
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky, aky, akx
    use species, only: nspec, spec
    use dist_fn_arrays, only: gnew, aj0, vpa, g_adjust
    use run_parameters, only: fphi, fapar, fbpar, beta
    use gs2_layouts, only: g_lo
    use collisions, only: init_lorentz, init_ediffuse, init_lorentz_conserve, init_diffuse_conserve
    use collisions, only: etol, ewindow, etola, ewindowa
    use collisions, only: vnmult, vary_vnew

    ! TEMP FOR TESTING -- MAB
!    use gs2_layouts, only: il_idx
!    use theta_grid, only: bmag, bmax
!    use le_grids, only: al

    implicit none

    integer :: ig, it, ik, iglo, isgn, ntrap

    integer, dimension (:,:), intent (out) :: erridx
    real, dimension (:,:), intent (out) :: errest
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar

    complex, dimension (:,:,:), allocatable :: phi_app, apar_app
    complex, dimension (:,:,:,:), allocatable :: phi_e, phi_l, phi_t
    complex, dimension (:,:,:,:), allocatable :: apar_e, apar_l, apar_t
    real, dimension (:), allocatable :: wgt, errtmp
    integer, dimension (:), allocatable :: idxtmp

    real :: vnmult_target

!    logical :: increase = .true., decrease = .false., first = .true.
    logical :: trap_flag

    allocate(wgt(nspec))
    allocate(errtmp(2))
    allocate(idxtmp(3))

    if (.not. allocated(kmax)) call broadcast (wdim)  ! Odd-looking statement.  BD

    if (fphi > epsilon(0.0)) then
       allocate(phi_app(-ntgrid:ntgrid,ntheta0,naky))
       allocate(phi_e(-ntgrid:ntgrid,ntheta0,naky,wdim))
       allocate(phi_l(-ntgrid:ntgrid,ntheta0,naky,ng2))
    end if

    if (fapar > epsilon(0.0)) then
       allocate(apar_app(-ntgrid:ntgrid,ntheta0,naky))
       allocate(apar_e(-ntgrid:ntgrid,ntheta0,naky,wdim))
       allocate(apar_l(-ntgrid:ntgrid,ntheta0,naky,ng2))
    end if

! first call to g_adjust converts gyro-averaged dist. fn. (g)
! into nonadiabatic part of dist. fn. (h)

    call g_adjust (gnew, phi, bpar, fphi, fbpar)

! take gyro-average of h at fixed total position (not g.c. position)
    if (fphi > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ! TEMP FOR TESTING -- MAB
!          il = il_idx(g_lo,iglo)
          do isgn = 1, 2
             do ig=-ntgrid, ntgrid
                ! TEMP FOR TESTING -- MAB
                g0(ig,isgn,iglo) = aj0(ig,iglo)*gnew(ig,isgn,iglo)
!                g0(ig,isgn,iglo) = sqrt(max(0.0,1.0-bmag(ig)*al(il)))
!                g0 = 1.
             end do
          end do
       end do

       wgt = spec%z*spec%dens

       call integrate_species (g0, wgt, phi_app)

! integrates dist fn of each species over v-space
! after dropping an energy grid point and returns
! phi_e, which contains the integral approximations
! to phi for each point dropped

       call eint_error (g0, wgt, phi_e)

! integrates dist fn of each species over v-space
! after dropping an untrapped lambda grid point and returns phi_l.
! phi_l contains ng2 approximations for the integral over lambda that
! come from dropping different pts from the gaussian quadrature grid

       call lint_error (g0, wgt, phi_l)

! next loop gets error estimate for trapped particles, if there are any

       if (nlambda > ng2 .and. new_trap_int) then
          ntrap = nlambda - ng2
          allocate(phi_t(-ntgrid:ntgrid,ntheta0,naky,ntrap))       
          phi_t = 0.0
          call trap_error (g0, wgt, phi_t)
       end if

    end if

    if (fapar > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             do ig=-ntgrid, ntgrid
                g0(ig,isgn,iglo) = aj0(ig,iglo)*vpa(ig,isgn,iglo)*gnew(ig,isgn,iglo)
             end do
          end do
       end do
       

       wgt = 2.0*beta*spec%z*spec%dens*sqrt(spec%temp/spec%mass)
       call integrate_species (g0, wgt, apar_app)

       call eint_error (g0, wgt, apar_e)
       call lint_error (g0, wgt, apar_l)
       if (nlambda > ng2 .and. new_trap_int) then
          ntrap = nlambda - ng2
          allocate(apar_t(-ntgrid:ntgrid,ntheta0,naky,ntrap))       
          apar_t = 0.0
          call trap_error (g0, wgt, apar_t)
       end if
       
    end if

! second call to g_adjust converts from h back to g

    call g_adjust (gnew, phi, bpar, -fphi, -fbpar)

    if (.not. allocated(kmax)) then
       allocate (kmax(ntheta0, naky))
       do ik = 1, naky
          do it = 1, ntheta0
             kmax(it,ik) = max(akx(it),aky(ik))
          end do
       end do
    end if
    
    errest = 0.0
    erridx = 0
    
    if (fphi > epsilon(0.0)) then

       call estimate_error (phi_app, phi_e, kmax, errtmp, idxtmp)
       errest(1,:) = errtmp
       erridx(1,:) = idxtmp
       
       call estimate_error (phi_app, phi_l, kmax, errtmp, idxtmp)
       errest(2,:) = errtmp
       erridx(2,:) = idxtmp

       if (nlambda > ng2 .and. new_trap_int) then       
          call estimate_error (phi_app, phi_t, kmax, errtmp, idxtmp, trap_flag)
          errest(3,:) = errtmp
          erridx(3,:) = idxtmp
          deallocate (phi_t)
       end if

    end if
    
    if (fapar > epsilon(0.0)) then

       call estimate_error (apar_app, apar_e, kmax, errtmp, idxtmp)
       errest(4,:) = errtmp
       erridx(4,:) = idxtmp

       call estimate_error (apar_app, apar_l, kmax, errtmp, idxtmp)
       errest(5,:) = errtmp
       erridx(5,:) = idxtmp
    end if

    if (vary_vnew) then
       vnmult_target = vnmult(2)
       
       if (errest(1,2) > etol + ewindow .or. errest(4,2) > etola + ewindowa) then
          call get_vnewk (vnmult(2), vnmult_target, increase)
       else if (errest(1,2) < etol - ewindow .and. errest(4,2) < etola - ewindowa) then
          call get_vnewk (vnmult(2), vnmult_target, decrease)
       end if

       call init_ediffuse (vnmult_target)
       call init_diffuse_conserve
    end if
       
    if (vary_vnew) then
       vnmult_target = vnmult(1)
       
       if (errest(2,2) > etol + ewindow .or. errest(3,2) > etol + ewindow &
            .or. errest(5,2) > etola + ewindowa) then
          call get_vnewk (vnmult(1), vnmult_target, increase)
       else if (errest(2,2) < etol - ewindow .and. errest(5,2) < etola - ewindowa .and. &
            (errest(3,2) < etol - ewindow .or. .not.trap_flag)) then
          call get_vnewk (vnmult(1), vnmult_target, decrease)
       end if
       
       call init_lorentz (vnmult_target)
       call init_lorentz_conserve
    end if
    
    deallocate (wgt, errtmp, idxtmp)
    if (fphi > epsilon(0.0)) deallocate(phi_app, phi_e, phi_l)
    if (fapar > epsilon(0.0)) deallocate(apar_app, apar_e, apar_l)

  end subroutine get_verr

  subroutine get_vnewk (vnm, vnm_target, incr)

    use collisions, only: vnfac, vnslow
    use mp, only: proc0

    implicit none

    logical, intent (in) :: incr
    real, intent (in) :: vnm
    real, intent (out) :: vnm_target

    if (incr) then
       vnm_target = vnm * vnfac
    else 
!       vnmult_target =  vnmult * ((1.0-1.0/vnslow)/vnfac + 1.0/vnslow)
       vnm_target =  vnm * vnslow
    end if

  end subroutine get_vnewk

  subroutine estimate_error (app1, app2, kmax_local, errest, erridx, trap)

    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use le_grids, only: ng2, jend

    implicit none

    complex, dimension (-ntgrid:,:,:), intent (in) :: app1
    complex, dimension (-ntgrid:,:,:,:), intent (in) :: app2
    real, dimension (:,:), intent (in) :: kmax_local
    logical, intent (in), optional :: trap
    real, dimension (:), intent (out) :: errest
    integer, dimension (:), intent (out) :: erridx

    integer :: ik, it, ig, ipt, end_idx
    integer :: igmax, ikmax, itmax, gpcnt
    real :: gdsum, gdmax, gpavg, gnsum, gsmax, gpsum, gptmp

    igmax = 0; ikmax = 0; itmax = 0
    gdsum = 0.0; gdmax = 0.0; gpavg = 0.0; gnsum = 0.0; gsmax = 1.0
    do ik = 1, naky
       do it = 1, ntheta0
          do ig=-ntgrid,ntgrid
             if (.not.present(trap) .or. jend(ig) > ng2+1) then
                gpcnt = 0; gpsum = 0.0

                   if (present(trap)) then
                      end_idx = jend(ig)-ng2
                   else
                      end_idx = size(app2,4)
                   end if

                   do ipt=1,end_idx
                      
                      gptmp = kmax_local(it,ik)*cabs(app1(ig,it,ik) - app2(ig,it,ik,ipt))
                      gpsum = gpsum + gptmp
                      gpcnt = gpcnt + 1
                      
                   end do
                   
                   gpavg = gpsum/gpcnt
                   
                   if (gpavg > gdmax) then
                      igmax = ig
                      ikmax = ik
                      itmax = it
                      gdmax = gpavg
                      gsmax = kmax_local(it,ik)*cabs(app1(ig,it,ik))
                   end if
                   
                   gnsum = gnsum + gpavg
                   gdsum = gdsum + kmax_local(it,ik)*cabs(app1(ig,it,ik))

!                end if
             end if
          end do
       end do
    end do
       
    gdmax = gdmax/gsmax
       
    erridx(1) = igmax
    erridx(2) = ikmax
    erridx(3) = itmax
    errest(1) = gdmax
    errest(2) = gnsum/gdsum

  end subroutine estimate_error

  subroutine get_gtran (geavg, glavg, gtavg, phi, bpar, istep)

    use le_grids, only: legendre_transform, negrid, nlambda, ng2, nesub, jend
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use dist_fn_arrays, only: gnew, aj0, g_adjust
    use run_parameters, only: fphi, fbpar
    use gs2_layouts, only: g_lo
    use mp, only: proc0, broadcast

    implicit none

    integer, intent (in) :: istep
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    real, intent (out) :: geavg, glavg, gtavg

    real, dimension (:), allocatable :: gne2, gnl2, gnt2
    complex, dimension (:,:,:,:,:), allocatable :: getran, gltran, gttran

    real :: genorm, gemax, genum, gedenom
    real :: glnorm, glmax, glnum, gldenom
    real :: gtnorm, gtmax, gtnum, gtdenom
    integer :: ig, it, ik, is, ie, il, iglo, isgn

    allocate(gne2(0:nesub-1))
    allocate(gnl2(0:ng2-1))

    genorm = 0.0 ; glnorm = 0.0
    gne2  = 0.0 ; gnl2 = 0.0
    gemax = 0.0 ; glmax = 0.0
    genum = 0.0 ; gedenom = 0.0
    glnum = 0.0 ; gldenom = 0.0
    gtnum = 0.0 ; gtdenom = 0.0

    allocate(getran(0:nesub-1,-ntgrid:ntgrid,ntheta0,naky,nspec))
    allocate(gltran(0:ng2-1,-ntgrid:ntgrid,ntheta0,naky,nspec))
 
    getran = 0.0; gltran = 0.0

    if (nlambda-ng2 > 0) then
       allocate(gnt2(0:2*(nlambda-ng2-1)))
       allocate(gttran(0:2*(nlambda-ng2-1),-ntgrid:ntgrid,ntheta0,naky,nspec))      
       gtnorm = 0.0 ; gnt2 = 0.0 ; gtmax = 0.0 ; gttran = 0.0
    end if

! transform from g to h
    call g_adjust (gnew, phi, bpar, fphi, fbpar)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
             g0(ig,isgn,iglo) = aj0(ig,iglo)*gnew(ig,isgn,iglo)
          end do
       end do
    end do

! perform legendre transform on dist. fn. to obtain coefficients
! used when expanding dist. fn. in legendre polynomials 
    if (allocated(gttran)) then
       call legendre_transform (g0, getran, gltran, istep, gttran)
    else
       call legendre_transform (g0, getran, gltran, istep)
    end if

! transform from h back to g
    call g_adjust (gnew, phi, bpar, -fphi, -fbpar)

    if (proc0) then
       do is=1,nspec
          do ik=1,naky
             do it=1,ntheta0
                do ig=-ntgrid,ntgrid
                   do ie=0,nesub-1
                      gne2(ie) = real(getran(ie,ig,it,ik,is)*conjg(getran(ie,ig,it,ik,is)))
                   end do
                   do il=0,ng2-1
                      gnl2(il) = real(gltran(il,ig,it,ik,is)*conjg(gltran(il,ig,it,ik,is)))
                   end do
                   genorm = maxval(gne2)
                   if (nesub < 3) then
                      gemax = gne2(size(gne2)-1)
                   else
                      gemax = maxval(gne2(nesub-3:nesub-1))
                   end if
                   glnorm = maxval(gnl2)
                   glmax = maxval(gnl2(ng2-3:ng2-1))

                   genum = genum + gemax
                   gedenom = gedenom + genorm
                   glnum = glnum + glmax
                   gldenom = gldenom + glnorm

                   if (nlambda > ng2) then
                      do il=0,2*(jend(ig)-ng2-1)
                         gnt2(il) = real(gttran(il,ig,it,ik,is)*conjg(gttran(il,ig,it,ik,is)))
                      end do
                      gtnorm = maxval(gnt2(0:2*(jend(ig)-ng2-1)))
                      if (jend(ig) > ng2+1) then
                         gtmax = maxval(gnt2(2*(jend(ig)-ng2-1)-2:2*(jend(ig)-ng2-1)))
                      else
                         gtmax = gnt2(0)
                      end if
                      
                      gtnum = gtnum + gtmax
                      gtdenom = gtdenom + gtnorm
                   end if

                end do
             end do
          end do
       end do
       geavg = genum/gedenom
       glavg = glnum/gldenom
       if (nlambda > ng2) gtavg = gtnum/gtdenom
    end if

    call broadcast (geavg)
    call broadcast (glavg)
    if (nlambda > ng2) then
       call broadcast (gtavg)
       deallocate (gnt2, gttran)
    end if

    deallocate(gne2, gnl2)    
    deallocate(getran, gltran)

  end subroutine get_gtran

  subroutine write_poly (phi, bpar, last, istep)

    use file_utils, only: open_output_file, close_output_file
    use dist_fn_arrays, only: gnew, aj0, g_adjust
    use run_parameters, only: fphi, fbpar
    use gs2_layouts, only: g_lo
    use mp, only: proc0
    use le_grids, only: nterp, negrid, lagrange_interp, xloc
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky

    implicit none

    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    logical, intent (in) :: last
    integer, intent (in) :: istep

    complex, dimension (:,:,:,:,:), allocatable :: lpoly
    integer :: iglo, isgn, ix, ig, tsize
!    logical, save :: first = .true.
    integer, save :: interp_unit

    allocate(lpoly(-ntgrid:ntgrid,ntheta0,naky,negrid,2*nterp-1))

    if (proc0) then
!       if (first) then 
       if (.not. lpolinit) then
          call open_output_file (interp_unit, ".interp")
          lpolinit = .true.
!          first = .false.
       end if
    end if

! why not phinew, bparnew?   MAB
    call g_adjust (gnew, phi, bpar, fphi, fbpar)

! take gyro-average of h at fixed total position (not g.c. position)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
             g0(ig,isgn,iglo) = aj0(ig,iglo)*gnew(ig,isgn,iglo)
          end do
       end do
    end do

! computes and returns lagrange interpolating polynomial for g0
    call lagrange_interp (g0, lpoly, istep)

    call g_adjust (gnew, phi, bpar, -fphi, -fbpar)

    tsize = 2*nterp-1

! write to file lagrange interpolating function value at interpolating pts specified by xloc
    if (proc0) then
       do ix = 1, tsize
          write(interp_unit,*) xloc(0,ix), real(lpoly(0,1,1,1,ix)), cos(0.1*istep*xloc(0,ix)+1.0), ix
       end do
       write(interp_unit,*)
       if (last) call close_output_file (interp_unit)
    end if

    deallocate(lpoly)

  end subroutine write_poly

  subroutine write_f (last)

    use mp, only: proc0, send, receive
    use file_utils, only: open_output_file, close_output_file, get_unused_unit
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx, il_idx, ie_idx
    use gs2_layouts, only: idx_local, proc_id
    use le_grids, only: al, energy, forbid, negrid, nlambda
    use egrid, only: zeroes, x0
    use theta_grid, only: bmax, bmag
    use gs2_time, only: user_time
    use dist_fn_arrays, only: g, gnew
    use species, only: nspec

    integer :: iglo, ik, it, is
    integer :: ie, il, ig
    integer, save :: unit
    real :: vpa, vpe
    complex, dimension(2) :: gtmp
!    logical :: first = .true.
    logical, intent(in)  :: last 

    real, dimension (:,:), allocatable, save :: xpts
    real, dimension (:), allocatable, save :: ypts

    if (.not. allocated(xpts)) then
       allocate(xpts(negrid,nspec))
       allocate(ypts(nlambda))

    ! should really just get rid of xpts and ypts
!    if (first) then
       xpts(1:negrid-1,:) = zeroes(1:negrid-1,:)
       xpts(negrid,:) = x0(:)
!       xpts = 0.0
       ypts = 0.0
       ! change argument of bmag depending on which theta you want to write out
       do il=1,nlambda
          if (1.0-al(il)*bmag(0) .gt. 0.0) ypts(il) = sqrt(1.0-al(il)*bmag(0))
       end do
!    end if

       if (proc0) then
!       if (first) then 
          call open_output_file (unit, ".dist")
          write(unit, *) negrid*nlambda
!          first = .false.
       end if
    endif

    do iglo = g_lo%llim_world, g_lo%ulim_world
       ! writing out g(vpar,vperp) at ik=it=is=1, ig=0
       ik = ik_idx(g_lo, iglo) ; if (ik /= 1) cycle
       it = it_idx(g_lo, iglo) ; if (it /= 1) cycle
       is = is_idx(g_lo, iglo) ; if (is /= 1) cycle
       ie = ie_idx(g_lo, iglo) 
       ig = 0
       il = il_idx(g_lo, iglo)
       if (idx_local (g_lo, ik, it, il, ie, is)) then
          if (proc0) then 
             gtmp = gnew(ig,:,iglo)
          else
             call send (gnew(ig,:,iglo), 0)
          endif
       else if (proc0) then
          call receive (gtmp, proc_id(g_lo, iglo))
       endif
       if (proc0) then
          if (.not. forbid(ig,il)) then
             vpa = sqrt(energy(ie,is)*max(0.0, 1.0-al(il)*bmag(ig)))
             vpe = sqrt(energy(ie,is)*al(il)*bmag(ig))
             write (unit, "(8(1x,e12.6))") vpa, vpe, energy(ie,is), al(il), &
                  xpts(ie,is), ypts(il), real(gtmp(1)), real(gtmp(2))
          end if
       end if
    end do
    if (proc0) write (unit, *)
    if (last .and. proc0) call close_output_file (unit)
    
  end subroutine write_f

  subroutine write_fyx (phi,bpar,last)

    use mp, only: proc0, send, receive, barrier
    use file_utils, only: open_output_file, close_output_file, flush_output_file
    use gs2_layouts, only: il_idx, ig_idx, ik_idx, it_idx, is_idx, isign_idx, ie_idx
    use gs2_layouts, only: idx_local, proc_id, yxf_lo, accelx_lo, g_lo
    use le_grids, only: al, energy, forbid, negrid, nlambda
    use kt_grids, only: naky, ntheta0, nx, ny
    use theta_grid, only: bmag, ntgrid
    use species, only: nspec
    use dist_fn_arrays, only: gnew, g_adjust
    use nonlinear_terms, only: accelerated
    use gs2_transforms, only: transform2, init_transforms
    use run_parameters, only: fphi, fbpar

    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    logical, intent (in) :: last

    real, dimension (:,:), allocatable :: grs, gzf
    real, dimension (:,:,:), allocatable :: agrs, agzf
    real, dimension (:), allocatable :: agp0, agp0zf
    real :: gp0, gp0zf
    integer :: ig, it, ik, il, ie, is, iyxlo, isign, ia, iaclo, iglo, acc
    integer, save :: unit
!    logical :: first = .true.

    if (accelerated) then
       allocate (agrs(2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc))          
       allocate (agzf(2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc))          
       allocate (agp0(2), agp0zf(2))
       agrs = 0.0; agzf = 0.0; agp0 = 0.0; agp0zf = 0.0
    else
       allocate (grs(yxf_lo%ny, yxf_lo%llim_proc:yxf_lo%ulim_alloc))
       allocate (gzf(yxf_lo%ny, yxf_lo%llim_proc:yxf_lo%ulim_alloc))
       grs = 0.0; gzf = 0.0; gp0 = 0.0; gp0zf = 0.0
    end if

!    if (first) then
    if (.not. fyxinit) then
       if (proc0) then
          acc = 0
          call open_output_file (unit,".yxdist")
          if (accelerated) acc = 1
          write(unit,*) 2*negrid*nlambda, bmag(0), acc
       end if

!       first = .false.
       fyxinit = .true.
    end if

    call g_adjust (gnew, phi, bpar, fphi, fbpar)

    g0 = gnew

    if (accelerated) then
       call transform2 (g0, agrs, ia)
    else
       call transform2 (g0, grs)
    end if

    g0 = 0.0
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       if (ik == 1) g0(:,:,iglo) = gnew(:,:,iglo)
    end do

    call g_adjust (gnew, phi, bpar, -fphi, -fbpar)

    if (accelerated) then
       call transform2 (g0, agzf, ia)
    else
       call transform2 (g0, gzf)
    end if

    if (accelerated) then
       do iaclo=accelx_lo%llim_world, accelx_lo%ulim_world
          it = it_idx(accelx_lo, iaclo)
          ik = ik_idx(accelx_lo, iaclo)
          if (it == 1 .and. ik == 1) then
             il = il_idx(accelx_lo, iaclo)
             ig = 0
             if (.not. forbid(ig,il)) then
                ie = ie_idx(accelx_lo, iaclo)
                is = is_idx(accelx_lo, iaclo)

                if (proc0) then
                   if (.not. idx_local(accelx_lo, ik, it, il, ie, is)) then
                      call receive (agp0, proc_id(accelx_lo, iaclo))
                      call receive (agp0zf, proc_id(accelx_lo, iaclo))
                   else
                      agp0 = agrs(ig+ntgrid+1, :, iaclo)
                      agp0zf = agzf(ig+ntgrid+1, :, iaclo)
                   end if
                else if (idx_local(accelx_lo, ik, it, il, ie, is)) then
                   call send (agrs(ig+ntgrid+1, :, iaclo), 0)
                   call send (agzf(ig+ntgrid+1, :, iaclo), 0)
                end if
                
                if (proc0) then
                   write (unit, "(6(1x,e12.6))") energy(ie,is), al(il), &
                        agp0(1), agp0(2), agp0zf(1), agp0zf(2)
                end if
             end if
          end if
          call barrier
       end do
       deallocate(agrs, agzf, agp0, agp0zf)
    else
       do iyxlo=yxf_lo%llim_world, yxf_lo%ulim_world
          
          ig = ig_idx(yxf_lo, iyxlo)
          it = it_idx(yxf_lo, iyxlo)
          if (ig == 0 .and. it == 1) then
             il = il_idx(yxf_lo, iyxlo)
             if (.not. forbid(ig,il)) then
                ik = 1
                ie = ie_idx(yxf_lo, iyxlo)
                is = is_idx(yxf_lo, iyxlo)
                isign = isign_idx(yxf_lo, iyxlo)
                
                if (proc0) then
                   if (.not. idx_local(yxf_lo, ig, isign, it, il, ie, is)) then
                      call receive (gp0, proc_id(yxf_lo, iyxlo))
                      call receive (gp0zf, proc_id(yxf_lo, iyxlo))
                   else
                      gp0 = grs(ik, iyxlo)
                      gp0zf = gzf(ik, iyxlo)
                   end if
                else if (idx_local(yxf_lo, ig, isign, it, il, ie, is)) then
                   call send (grs(ik, iyxlo), 0)
                   call send (gzf(ik, iyxlo), 0)                   
                end if
                
                if (proc0) then
                   write (unit, "(4(1x,e12.6),i8)") energy(ie,is), al(il), &
                        gp0, gp0zf, isign
                end if
             end if
          end if
          call barrier
       end do
       deallocate (grs, gzf)
    end if

    if (proc0) call flush_output_file (unit, ".yxdist")

    if (proc0) then
       write(unit,*)
       if (last) call close_output_file (unit)
    end if

  end subroutine write_fyx

  subroutine collision_error (phi, bpar, last)
    
    use mp, only: proc0, send, receive, barrier
    use le_grids, only: ng2, jend, nlambda, al, forbid, lambda_map
    use theta_grid, only: ntgrid, bmag
    use dist_fn_arrays, only: gnew, aj0, g_adjust
    use run_parameters, only: fphi, fbpar
    use gs2_layouts, only: g_lo, lz_lo, ig_idx, idx_local, proc_id
    use gs2_layouts, only: ik_idx, ie_idx, is_idx, it_idx, il_idx
    use collisions, only: dtot, fdf, fdb
    use redistribute, only: gather, scatter
    use file_utils, only: open_output_file, close_output_file
    use gs2_time, only: user_time
    implicit none

    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    logical, intent (in) :: last

    integer :: je, te, ig, il, ip, ilz, ie, is, ik, it
    integer :: igmax, ikmax, itmax, iemax, ilmax, ismax
    integer, save :: unit
    complex, dimension (:), allocatable :: ltmp, ftmp
    complex, dimension (:,:), allocatable :: lcoll, fdcoll, glze
!    logical :: first = .true.
    real :: etmp, emax, etot, eavg, edenom, ltmax
    real :: time

    allocate (ltmp(2*nlambda), ftmp(2*nlambda))
    allocate (lcoll(2*nlambda,lz_lo%llim_proc:lz_lo%ulim_alloc))
    allocate (fdcoll(2*nlambda,lz_lo%llim_proc:lz_lo%ulim_alloc))
    allocate (glze(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))

    lcoll = 0.0; fdcoll = 0.0; glze = 0.0; ltmp = 0.0; ftmp = 0.0
    etmp = 0.0; emax = 0.0; etot = 0.0; eavg = 0.0; edenom = 1.0; ltmax = 1.0

!    if (first .and. proc0) then
    if (.not.cerrinit .and. proc0) then
       call open_output_file (unit,".cres")
!       first = .false.
       cerrinit = .true.
    end if

! convert gnew from g to h
    call g_adjust(gnew,phi,bpar,fphi,fbpar)

    g0 = gnew

! convert gnew from h back to g
    call g_adjust(gnew,phi,bpar,-fphi,-fbpar)

! map from g0(ig,isgn,iglo) to glze(il,ilz)
    call gather (lambda_map, g0, glze)

! loop over ig, isign, ik, it, ie, is
    do ilz = lz_lo%llim_proc, lz_lo%ulim_proc
       ig = ig_idx(lz_lo,ilz)
       
       if (jend(ig) == 0) then      ! no trapped particles
! je = number of + vpa grid pts
! te = number of + and - vpa grid pts
          je = ng2
          te = 2*ng2
          
! find d/d(xi) ((1+xi**2)( d g(xi)/ d(xi) )) at each xi
! using lagrange (lcoll) and finite difference (fdcoll)
          il = 1
          do ip = il, il+2
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip)*glze(ip,ilz)
          end do

          il = 2
          do ip = il-1, il+1
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip-il+2)*glze(ip,ilz)
          end do

          do il=3,ng2
             do ip=il-2,il+2
                lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip-il+3)*glze(ip,ilz)
             end do
          end do

          do il=ng2+1, 2*ng2-2
             do ip = il-2,il+2
                lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,2*ng2-il+1,il-ip+3)*glze(ip,ilz)
             end do
          end do

          il = 2*ng2-1
          do ip = il-1, il+1
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,2,il-ip+2)*glze(ip,ilz)
          end do

          il = 2*ng2
          do ip = il-2, il
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,1,il-ip+1)*glze(ip,ilz)
          end do

! deal with xi from 1-eps -> eps
          do il=2,ng2
             fdcoll(il,ilz) = fdf(ig,il)*(glze(il+1,ilz) - glze(il,ilz)) - fdb(ig,il)*(glze(il,ilz) - glze(il-1,ilz))
          end do
! deal with xi from -eps -> -1+eps
          do il=ng2+1, 2*ng2-1
             fdcoll(il,ilz) = fdb(ig,2*ng2-il+1)*(glze(il+1,ilz) - glze(il,ilz)) - fdf(ig,2*ng2-il+1)*(glze(il,ilz) - glze(il-1,ilz))
          end do

          fdcoll(1,ilz) = fdf(ig,1)*(glze(2,ilz) - glze(1,ilz))
          fdcoll(2*ng2,ilz) = -fdf(ig,1)*(glze(2*ng2,ilz) - glze(2*ng2-1,ilz))

       else       ! trapped particle runs          
          je = jend(ig)
          te = 2*je - 1

          il = 1
          do ip = il, il+2
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip)*glze(ip,ilz)
          end do

          il = 2
          do ip = il-1, il+1
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip-il+2)*glze(ip,ilz)
          end do

          do il=3,je
             do ip=il-2,il+2
                lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip-il+3)*glze(ip,ilz)
             end do
          end do

          do il=je+1, te-2
             do ip = il-2,il+2
                lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,te-il+1,il-ip+3)*glze(ip,ilz)
             end do
          end do

          il = te-1
          do ip = il-1, il+1
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,2,il-ip+2)*glze(ip,ilz)
          end do

          il = te
          do ip = il-2, il
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,1,il-ip+1)*glze(ip,ilz)
          end do

! is il=je handled correctly here?
          do il=2,je
             fdcoll(il,ilz) = fdf(ig,il)*(glze(il+1,ilz) - glze(il,ilz)) - fdb(ig,il)*(glze(il,ilz) - glze(il-1,ilz))
          end do
          do il=je+1,te-1
             fdcoll(il,ilz) = fdb(ig,te-il+1)*(glze(il+1,ilz) - glze(il,ilz)) - fdf(ig,te-il+1)*(glze(il,ilz) - glze(il-1,ilz))
          end do
          
          fdcoll(1,ilz) = fdf(ig,1)*(glze(2,ilz) - glze(1,ilz))
          fdcoll(te,ilz) = -fdf(ig,1)*(glze(te,ilz) - glze(te-1,ilz))
          
       end if
    end do

    time = user_time

    do ilz=lz_lo%llim_world, lz_lo%ulim_world
       ig = ig_idx(lz_lo, ilz)
       ik = ik_idx(lz_lo, ilz)
       it = it_idx(lz_lo, ilz)
       ie = ie_idx(lz_lo, ilz)
       is = is_idx(lz_lo, ilz)
       je = jend(ig)

       if (je == 0) then
          te = 2*ng2
       else
          te = 2*je-1
       end if

       if (idx_local (lz_lo, ilz)) then
          if (proc0) then 
             ltmp = lcoll(:,ilz)
             ftmp = fdcoll(:,ilz)
          else
             call send (lcoll(:,ilz), 0)
             call send (fdcoll(:,ilz), 0)
          endif
       else if (proc0) then
          call receive (ltmp, proc_id(lz_lo, ilz))
          call receive (ftmp, proc_id(lz_lo, ilz))
       endif
       call barrier

       do il=1,te
          etmp = cabs(ltmp(il) - ftmp(il))
          
          if (etmp > emax) then
             emax = etmp
             ltmax = cabs(ltmp(il))
             ikmax = ik
             itmax = it
             iemax = ie
             ismax = is
             ilmax = il
             igmax = ig
          end if
          
          etot = etot + etmp
          edenom = edenom + cabs(ltmp(il))
       end do
    end do

    eavg = etot/edenom
    emax = emax/ltmax

    if (proc0) then
       write(unit,"((1x,e12.6),6(i8),2(1x,e12.6))") time, &
            igmax, ikmax, itmax, iemax, ilmax, ismax, emax, eavg
       if (last) then
          call close_output_file (unit)
       end if
    end if

    deallocate (lcoll, fdcoll, glze, ltmp, ftmp)
    
  end subroutine collision_error

  subroutine boundary(linked)

    logical :: linked

    call init_dist_fn
    linked = boundary_option_switch == boundary_option_linked

  end subroutine boundary

! This subroutine only returns epar correctly for linear runs.
  subroutine get_epar (phi, apar, phinew, aparnew, epar)

    use theta_grid, only: ntgrid, delthet, gradpar
    use run_parameters, only: tunits, fphi, fapar
    use gs2_time, only: code_dt
    use kt_grids, only: naky, ntheta0
    complex, dimension(-ntgrid:,:,:) :: phi, apar, phinew, aparnew, epar
    complex :: phi_m, apar_m

    integer :: ig, ik, it

    do ik = 1, naky
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid-1
             ! ignoring decentering in time and space for now
             phi_m = 0.5*(phi(ig+1,it,ik)-phi(ig,it,ik) + &
                  phinew(ig+1,it,ik)-phinew(ig,it,ik))*fphi
             apar_m = 0.5*(aparnew(ig+1,it,ik)+aparnew(ig,it,ik) & 
                  -apar(ig+1,it,ik)-apar(ig,it,ik))*fapar
             
             epar(ig,it,ik) = -phi_m/delthet(ig)* &
                  0.5*(abs(gradpar(ig)) + abs(gradpar(ig+1))) &
                  -apar_m/tunits(ik)/code_dt
          end do
       end do
    end do    

  end subroutine get_epar

  subroutine find_leftmost_link (iglo, iglo_left, ipleft)

    integer, intent (in) :: iglo
    integer, intent (in out) :: iglo_left, ipleft
    integer :: iglo_star, iglo_left_star, ipleft_star

    iglo_star = iglo
    do 
       call get_left_connection (iglo_star, iglo_left_star, ipleft_star)
    
       if (ipleft_star == -1) exit
       iglo_star = iglo_left_star
       iglo_left = iglo_left_star
       ipleft = ipleft_star
    end do

  end subroutine find_leftmost_link

  subroutine find_rightmost_link (iglo, iglo_right, ipright)

    integer, intent (in) :: iglo
    integer, intent (in out) :: iglo_right, ipright
    integer :: iglo_star, iglo_right_star, ipright_star

    iglo_star = iglo
    do 
       call get_right_connection (iglo_star, iglo_right_star, ipright_star)
    
       if (ipright_star == -1) exit
       iglo_star = iglo_right_star
       iglo_right = iglo_right_star
       ipright = ipright_star
    end do

  end subroutine find_rightmost_link

  subroutine timer (i, place)
    
    character (len=10) :: zdate, ztime, zzone
    character (*) :: place
    integer, intent (in) :: i
    integer, dimension(8) :: ival
    real, save :: told=0., tnew=0.
    
    call date_and_time (zdate, ztime, zzone, ival)
    tnew = ival(5)*3600.+ival(6)*60.+ival(7)+ival(8)/1000.
    if (i == 0) told = -1.
    if (told > 0.) then
       print *, ': Elapsed time = ',tnew-told,' seconds in '//trim(place)
    end if
    told = tnew
  end subroutine timer

  subroutine dot (a, anew, adot, fac)

! Get a theta-centered and time-centered estimate of the time derivative 
! of a field.
! 
! tunits(ky) == 1. unless the "omega_*" units are chosen.
! omega_* units normalize time by an additional factor of ky.
!    

    use run_parameters, only: tunits
    use gs2_time, only: code_dt
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid

    implicit none
    complex, intent (in), dimension (-ntgrid:,:,:) :: a, anew
    complex, intent (out), dimension (-ntgrid:,:,:) :: adot
    real, intent (in) :: fac
    real :: dtinv
    integer :: ig, it, ik

    do ik=1,naky
       dtinv = 1./(code_dt*tunits(ik))
       do it=1,ntheta0
          do ig=-ntgrid,ntgrid-1
             adot(ig,it,ik) = 0.5*fac*(anew(ig+1,it,ik)+anew(ig,it,ik) - &
                  (a(ig+1,it,ik)+a(ig,it,ik)))*dtinv
          end do
       end do
    end do
    
  end subroutine dot

  complex function fdot (fl, fr, fnewl, fnewr, dtinv)

    complex, intent (in) :: fl, fr, fnewl, fnewr
    real, intent (in) :: dtinv
    
    fdot = 0.5*(fnewl+fnewr-fl-fr)*dtinv

  end function fdot

! construct time and space-centered quantities
! (should use actual bakdif and fexpr values?)
!
  complex function favg (fl, fr, fnewl, fnewr)

    complex, intent (in) :: fl, fr, fnewl, fnewr
    
    favg = 0.25*(fnewl+fnewr+fl+fr)

  end function favg
 
   complex function fdot_t (f,fnew, dtinv)

    complex, intent (in) :: f, fnew
    real, intent (in) :: dtinv
    
    fdot_t = (fnew-f)*dtinv

  end function fdot_t

 complex function favg_x (fl, fr)

    complex, intent (in) :: fl, fr
    
    favg_x = 0.5*(fl+fr)

  end function favg_x

  subroutine init_mom_coeff
    use gs2_layouts, only: g_lo
    use species, only: nspec, spec
    use kt_grids, only: nakx => ntheta0, naky
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: aj0,aj1,vperp2,kperp2
    use dist_fn_arrays, only: vpa, vperp2
    use le_grids, only: integrate_moment
    integer :: i
    integer :: isgn, iglo
    integer :: it,ik,is
    complex, allocatable :: coeff0(:,:,:,:)
    complex, allocatable :: gtmp(:,:,:)
    real, allocatable :: wgt(:)
!    logical, parameter :: analytical = .false.
    logical, parameter :: analytical = .true.
!    logical, save :: initialized=.false.
    real :: bsq
    
!    if (initialized) return
    if (mominit) return
    mominit = .true.

    if(.not.allocated(mom_coeff)) &
         & allocate(mom_coeff(nakx,naky,nspec,ncnt_mom_coeff))
    if(.not.allocated(mom_coeff_npara)) &
         & allocate(mom_coeff_npara(nakx,naky,nspec))
    if(.not.allocated(mom_coeff_nperp)) &
         & allocate(mom_coeff_nperp(nakx,naky,nspec))
    if(.not.allocated(mom_coeff_tpara)) &
         & allocate(mom_coeff_tpara(nakx,naky,nspec))
    if(.not.allocated(mom_coeff_tperp)) &
         & allocate(mom_coeff_tperp(nakx,naky,nspec))
    if(.not.allocated(mom_shift_para)) &
         & allocate(mom_shift_para(nakx,naky,nspec))
    if(.not.allocated(mom_shift_perp)) &
         & allocate(mom_shift_perp(nakx,naky,nspec))

    mom_coeff(:,:,:,:)=0.
    mom_coeff_npara(:,:,:)=0. ; mom_coeff_nperp(:,:,:)=0.
    mom_coeff_tpara(:,:,:)=0. ; mom_coeff_tperp(:,:,:)=0.
    mom_shift_para(:,:,:)=0.  ; mom_shift_perp(:,:,:)=0.

    allocate(wgt(-ntgrid:ntgrid))
    allocate(coeff0(-ntgrid:ntgrid,nakx,naky,nspec))
    allocate(gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    wgt(:)=0.
    coeff0(:,:,:,:)=cmplx(0.,0.)
    gtmp(:,:,:)=cmplx(0.,0.)

    if (analytical) then
       do it=1,nakx
          do ik=1,naky
             do is=1,nspec
                bsq=.25*spec(is)%smz**2*kperp2(0,it,ik)
                mom_coeff(it,ik,is,1) = exp(-bsq)
                mom_coeff(it,ik,is,2) = exp(-bsq) *.5
                mom_coeff(it,ik,is,3) = exp(-bsq) *(1.-bsq)
                mom_coeff(it,ik,is,4) = exp(-bsq) *.75
                mom_coeff(it,ik,is,5) = exp(-bsq) *(1.-bsq)*.5
                mom_coeff(it,ik,is,6) = exp(-bsq) *.5
                mom_coeff(it,ik,is,7) = exp(-bsq) *.25
                mom_coeff(it,ik,is,8) = exp(-bsq) *(1.-.5*bsq)
             end do
          end do
       end do
    else
       do i = 1, ncnt_mom_coeff
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             do isgn = 1,2
                if(i==1) wgt(:)=aj0(:,iglo)
                if(i==2) wgt(:)=aj0(:,iglo)*vpa(:,isgn,iglo)**2
                if(i==3) wgt(:)=aj0(:,iglo)*vperp2(:,iglo)
                if(i==4) wgt(:)=aj0(:,iglo)*vpa(:,isgn,iglo)**4
                if(i==5) wgt(:)=aj0(:,iglo)*vpa(:,isgn,iglo)**2*vperp2(:,iglo)
                if(i==6) wgt(:)=vperp2(:,iglo)*aj1(:,iglo)
                if(i==7) wgt(:)=vperp2(:,iglo)*aj1(:,iglo)*vpa(:,isgn,iglo)**2
                if(i==8) wgt(:)=vperp2(:,iglo)*aj1(:,iglo)*vperp2(:,iglo)
                gtmp(-ntgrid:ntgrid,isgn,iglo) = wgt(-ntgrid:ntgrid)*cmplx(1.,0.)
             end do
          end do
          call integrate_moment(gtmp,coeff0,1)
          where(real(coeff0(0,1:nakx,1:naky,1:nspec)) == 0.)
             mom_coeff(1:nakx,1:naky,1:nspec,i)=1.
          elsewhere
             mom_coeff(1:nakx,1:naky,1:nspec,i)= &
                  & coeff0(0,1:nakx,1:naky,1:nspec)
          end where
       end do
    endif

    mom_shift_para(:,:,:)=mom_coeff(:,:,:,2)/mom_coeff(:,:,:,1)
    mom_shift_perp(:,:,:)=mom_coeff(:,:,:,3)/mom_coeff(:,:,:,1)

    mom_coeff_npara(:,:,:)=2.*mom_coeff(:,:,:,2)/mom_coeff(:,:,:,1)
    mom_coeff_nperp(:,:,:)=2.*mom_coeff(:,:,:,6)/mom_coeff(:,:,:,1)

    mom_coeff_tperp(:,:,:)= &
         & (mom_coeff(:,:,:,5)-mom_shift_perp(:,:,:)*mom_coeff(:,:,:,2)) / &
         & (mom_coeff(:,:,:,8)-mom_shift_perp(:,:,:)*mom_coeff(:,:,:,6))
    mom_coeff_tpara(:,:,:)= &
         & (mom_coeff(:,:,:,7)-mom_shift_para(:,:,:)*mom_coeff(:,:,:,6)) / &
         & (mom_coeff(:,:,:,4)-mom_shift_para(:,:,:)*mom_coeff(:,:,:,2))

    deallocate(gtmp,coeff0,wgt)
    
!    initialized=.true.
  end subroutine init_mom_coeff

  subroutine finish_dist_fn

    use dist_fn_arrays, only: ittp, vpa, vpac, vpacgen, vperp2, vpar, vpargen
    use dist_fn_arrays, only: aj0, aj1, kperp2
    use dist_fn_arrays, only: g, gnew, kx_shift

    implicit none

    accelerated_x = .false. ; accelerated_v = .false.
    no_comm = .false.
    readinit = .false. ; bessinit = .false. ; kp2init = .false. ; connectinit = .false.
    feqinit = .false. ; lpolinit = .false. ; fyxinit = .false. ; cerrinit = .false. ; mominit = .false.
    increase = .true. ; decrease = .false.

    call reset_init

    if (allocated(fexp)) deallocate (fexp, bkdiff, bd_exp)
    if (allocated(ittp)) deallocate (ittp, wdrift, wdriftgen, wdriftttp)
    if (allocated(ittp)) deallocate (wdriftttpgen)
    if (allocated(vpa)) deallocate (vpa, vpac, vpacgen, vperp2, vpar, vpargen)
    if (allocated(wstar)) deallocate (wstar)
    if (allocated(aj0)) deallocate (aj0, aj1)
    if (allocated(kperp2)) deallocate (kperp2)
    if (allocated(a)) deallocate (a, b, r, ainv)
    if (allocated(l_links)) deallocate (l_links, r_links, n_links)
    if (allocated(M_class)) deallocate (M_class, N_class)
    if (allocated(itleft)) deallocate (itleft, itright)
    if (allocated(connections)) deallocate (connections)
    if (allocated(g_adj)) deallocate (g_adj)
    if (allocated(g)) deallocate (g, gnew, g0)
    if (allocated(gexp_1)) deallocate (gexp_1, gexp_2, gexp_3)
    if (allocated(g_h)) deallocate (g_h, save_h)
    if (allocated(kx_shift)) deallocate (kx_shift)
    if (allocated(jump)) deallocate (jump)
    if (allocated(ikx_indexed)) deallocate (ikx_indexed)
    if (allocated(ufac)) deallocate (ufac)
    if (allocated(gridfac1)) deallocate (gridfac1, gamtot, gamtot1, gamtot2)
    if (allocated(gamtot3)) deallocate (gamtot3)
    if (allocated(fl_avg)) deallocate (fl_avg)
    if (allocated(awgt)) deallocate (awgt)
    if (allocated(kmax)) deallocate (kmax)
    if (allocated(mom_coeff)) deallocate (mom_coeff, mom_coeff_npara, mom_coeff_nperp, &
         mom_coeff_tpara, mom_coeff_tperp, mom_shift_para, mom_shift_perp)

    ! gc_from_left, gc_from_right, links_p, links_h, wfb_p, wfb_h

  end subroutine finish_dist_fn

#ifdef LOWFLOW
  subroutine init_lowflow

    use constants, only: zi
    use dist_fn_arrays, only: vparterm, wdfac, vpac, wdttpfac
    use dist_fn_arrays, only: wstarfac, hneoc, vpar
    use species, only: spec, nspec
    use geometry, only: rhoc
    use theta_grid, only: theta, ntgrid, delthet, gradpar, bmag
    use theta_grid, only: gds23, gds24, gds24_noq, cvdrift_th, gbdrift_th
    use theta_grid, only: drhodpsi, qval, shat
    use le_grids, only: energy => energy_maxwell, al, negrid, nlambda, forbid, init_map
    use le_grids, only: get_flux_vs_theta_vs_vpa
    use kt_grids, only: theta0, ntheta0, naky, aky, akx
    use gs2_time, only: code_dt
    use gs2_layouts, only: g_lo, ik_idx, il_idx, ie_idx, is_idx, it_idx
    use run_parameters, only: tunits, wunits, rhostar, neo_test
    use lowflow, only: get_lowflow_terms
    use file_utils, only: open_output_file, close_output_file
    use mp, only: proc0

    implicit none

    integer, save :: neo_unit, neophi_unit, neovpth_unit
    integer :: it, ik, il, ie, is, isgn, iglo, ig
    real, dimension (:,:,:,:,:), allocatable :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
    real, dimension (:,:,:), allocatable :: vpadhdec, dhdec, dhdxic, cdfac, hneovpth
    real, dimension (:), allocatable :: tmp7, tmp8, tmp9

    allocate (vpadhdec (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (dhdec (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (dhdxic (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (cdfac (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    vpadhdec = 0. ; dhdec = 0. ; dhdxic = 0. ; cdfac = 0.

    if (.not. allocated(vparterm)) then
       allocate (vparterm(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (wdfac(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (hneoc(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (wstarfac(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (wdttpfac(-ntgrid:ntgrid,ntheta0,naky,negrid,nspec,2))
       allocate (wstar_neo(-ntgrid:ntgrid,ntheta0,naky))
    end if
    vparterm = 0. ; wdfac = 0. ; hneoc = 1. ; wstarfac = 0. ; wdttpfac = 0. ; wstar_neo = 0.

    allocate (tmp1(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp1 = 0.
    allocate (tmp2(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp2 = 0.
    allocate (tmp3(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp3 = 0.
    allocate (tmp4(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp4 = 0.
    allocate (tmp5(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp5 = 0.
    allocate (tmp6(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp6 = 0.
    allocate (tmp7(-ntgrid:ntgrid), tmp8(-ntgrid:ntgrid), tmp9(-ntgrid:ntgrid))
    tmp7 = 0. ; tmp8 = 0. ; tmp9 = 0.

    allocate (hneovpth(-ntgrid:ntgrid,negrid*nlambda,nspec)) ; hneovpth = 0.

    ! tmp1 is dH^{neo}/dE, tmp2 is dH^{neo}/dxi, tmp3 is vpa*dH^{neo}/dE,
    ! tmp4 is dH^{neo}/dr, tmp5 is dH^{neo}/dtheta, tmp6 is H^{neo}
    ! tmp7 is phi^{neo}/dr, tmp8 is dphi^{neo}/dtheta, and tmp9 phi^{neo}
    call get_lowflow_terms (theta, al, energy, bmag, tmp1, tmp2, tmp3, tmp4, &
         tmp5, tmp6, tmp7, tmp8, tmp9, lf_default, lf_decompose)
    
    if (proc0) then
       call open_output_file (neo_unit,".neodist")
       write (neo_unit,*) "# all quantities given at theta=0 for species 1"
       write (neo_unit,fmt='(10a12)') "# 1) vpa", "2) vpe", "3) energy", "4) vpa/v", &
            "5) dH/dE", "6) dH/dxi", "7) vpa*dH/dE", "8) dH/dr", "9) dH/dtheta", "10) H"
       do isgn = 1, 2
          do il = 1, nlambda
             do ie = 1, negrid
                if (.not. forbid(0,il)) then
                   write (neo_unit,'(10e12.4)') sign(sqrt(energy(ie)*(1.-al(il)*bmag(0))),1.5-real(isgn)), &
                        sqrt(energy(ie)*al(il)*bmag(0)), energy(ie), &
                        sign(sqrt(1.-al(il)*bmag(0)),1.5-real(isgn)), &
                        tmp1(0,il,ie,isgn,1), tmp2(0,il,ie,isgn,1), tmp3(0,il,ie,isgn,1), &
                        tmp4(0,il,ie,isgn,1), tmp5(0,il,ie,isgn,1), tmp6(0,il,ie,isgn,1)
                end if
             end do
             write (neo_unit,*)
          end do
       end do
       call close_output_file (neo_unit)
       
       call open_output_file (neovpth_unit,".neothvp")

       ! Get Fneo(theta,vpa)
       call get_flux_vs_theta_vs_vpa (tmp6, hneovpth)

       write (neovpth_unit,'(3a12)') '1) theta', '2) vpa', '3) Fneo'
       do ie = 1, negrid*nlambda
          do ig = -ntgrid, ntgrid
             write (neovpth_unit,'(3e12.4)'), theta(ig), sqrt(energy(negrid))*(1.-2.*(ie-1)/real(negrid*nlambda-1)), &
                  hneovpth(ig,ie,1)
          end do
          write (neovpth_unit,*)
       end do

       call close_output_file (neovpth_unit)

       call open_output_file (neophi_unit,".neophi")
       write (neophi_unit,*) "# 1) theta, 2) dphi/dr, 3) dphi/dtheta, 4) phi"
       do ig = -ntgrid, ntgrid
          write (neophi_unit,'(4e14.5)') theta(ig), tmp7(ig), tmp8(ig), tmp9(ig)
       end do
       call close_output_file (neophi_unit)
    end if
    
    ! if set neo_test flag to .true. in input file, GS2 exits after writing out
    ! neoclassical quantities of interest
    if (neo_test) stop
    
    ! intialize mappings from g_lo to e_lo and lz_lo or to le_lo to facilitate
    ! energy and lambda derivatives in parallel nonlinearity, etc.
# ifdef USE_LE_LAYOUT
    call init_map (use_lz_layout=.false., use_e_layout=.false., use_le_layout=.true., test=.false.)
# else
    call init_map (use_lz_layout=.true., use_e_layout=.true., use_le_layout=.false., test=.false.)
# endif

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1,2
          dhdec(:,isgn,iglo) = tmp1(:,il,ie,isgn,is)
          dhdxic(:,isgn,iglo) = tmp2(:,il,ie,isgn,is)
          vpadhdec(:,isgn,iglo) = tmp3(:,il,ie,isgn,is)
          hneoc(:,isgn,iglo) = 1.0+tmp6(:,il,ie,isgn,is)
       end do
       
       ! get cell-centered (in theta) values
       
       ! this is the contribution from dH^{neo}/dtheta (part of v_E dot grad F^{neo})
       ! takes care of part of Eq. 60 in MAB's GS2 notes
       if (aky(ik) == 0.0) then
          wstarfac(-ntgrid:ntgrid-1,1,iglo) = 0.25*akx(it)/shat &
               *(gds24(-ntgrid:ntgrid-1)+gds24(-ntgrid+1:ntgrid)) &
               *tmp5(-ntgrid:ntgrid-1,il,ie,1,is)*code_dt
          wstarfac(-ntgrid:ntgrid-1,2,iglo) = 0.25*akx(it)/shat &
               *(gds24(-ntgrid:ntgrid-1)+gds24(-ntgrid+1:ntgrid)) &
               *tmp5(-ntgrid:ntgrid-1,il,ie,2,is)*code_dt
       else
          wstarfac(-ntgrid:ntgrid-1,1,iglo) = 0.5*wunits(ik) &
               *(gds23(-ntgrid:ntgrid-1)+gds23(-ntgrid+1:ntgrid)+theta0(it,ik)*(gds24(-ntgrid:ntgrid-1)+gds24(-ntgrid+1:ntgrid))) &
               *tmp5(-ntgrid:ntgrid-1,il,ie,1,is)*code_dt
          wstarfac(-ntgrid:ntgrid-1,2,iglo) = 0.5*wunits(ik) &
               *(gds23(-ntgrid:ntgrid-1)+gds23(-ntgrid+1:ntgrid)+theta0(it,ik)*(gds24(-ntgrid:ntgrid-1)+gds24(-ntgrid+1:ntgrid))) &
               *tmp5(-ntgrid:ntgrid-1,il,ie,2,is)*code_dt
       end if

       ! this is the contribution from dH^{neo}/dr (part of v_E dot grad F^{neo})
       ! takes care of part of Eq. 60 in MAB's GS2 notes
       wstarfac(:,:,iglo) = wstarfac(:,:,iglo) + tmp4(:,il,ie,:,is)*code_dt*wunits(ik)
       
       ! this is the contribution from v_E^par . grad F0
       ! TMP FOR TESTING -- MAB
!        wstarfac(-ntgrid:ntgrid-1,1,iglo) = wstarfac(-ntgrid:ntgrid-1,1,iglo) &
!             - 0.5*zi*rhostar*(gds24_noq(-ntgrid:ntgrid-1)+gds24_noq(-ntgrid+1:ntgrid)) &
!             *drhodpsi*rhoc/qval*code_dt*(spec(is)%fprim+spec(is)%tprim*(energy(ie)-1.5))
!        wstarfac(-ntgrid:ntgrid-1,2,iglo) = wstarfac(-ntgrid:ntgrid-1,2,iglo) &
!             - 0.5*zi*rhostar*(gds24_noq(-ntgrid:ntgrid-1)+gds24_noq(-ntgrid+1:ntgrid)) &
!             *drhodpsi*rhoc/qval*code_dt*(spec(is)%fprim+spec(is)%tprim*(energy(ie)-1.5))
       
       if (.not. lf_default) then
          ! this is the contribution from the last term of the 2nd line of Eq. 43 in 
          ! MAB's GS2 notes (arises because the NEO dist. fn. is given at v/vt, and
          ! the vt varies radially, so v_E . grad F1 must take this into account)
          ! If lf_default is true this is taken care of by multiplying the vt normalization
          ! of Chebyshev polynomial argument by appropriate temperature factor
          ! when constructing neoclassical distribution function (see lowflow.f90).
          ! Otherwise, the below line of code takes care of it.
          wstarfac(:,:,iglo) = wstarfac(:,:,iglo) + code_dt*wunits(ik) &
               *spec(is)%tprim*energy(ie)*dhdec(:,:,iglo)
       end if
       
       wstarfac(ntgrid,:,iglo) = 0.0
       
       ! wdfac takes care of the gbdrift part of Eq. 60 of MAB's GS2 notes, as well
       ! as part of the curvature drift term of Eq. 54.  the other part of 
       ! the curvature drift term of Eq. 54 is dealt with by cdfac below.
       ! no code_dt in wdfac because it multiplies wdrift, which has code_dt in it
       wdfac(-ntgrid:ntgrid-1,1,iglo) = 0.5*dhdxic(-ntgrid:ntgrid-1,1,iglo)*vpac(-ntgrid:ntgrid-1,1,iglo)/energy(ie)**1.5 &
            - dhdec(-ntgrid:ntgrid-1,1,iglo)
       wdfac(-ntgrid:ntgrid-1,2,iglo) = 0.5*dhdxic(-ntgrid:ntgrid-1,2,iglo)*vpac(-ntgrid:ntgrid-1,2,iglo)/energy(ie)**1.5 &
            - dhdec(-ntgrid:ntgrid-1,2,iglo)
       wdfac(ntgrid,:,iglo) = 0.0

       ! takes care of part of curvature drift term in Eq. 54 of MAB's GS2 notes.
       ! no code_dt in cdfac because it multiples wcurv, which has code_dt in it.
       cdfac(-ntgrid:ntgrid-1,1,iglo) = -0.5*dhdxic(-ntgrid:ntgrid-1,1,iglo)*vpac(-ntgrid:ntgrid-1,1,iglo)/sqrt(energy(ie))
       cdfac(-ntgrid:ntgrid-1,2,iglo) = -0.5*dhdxic(-ntgrid:ntgrid-1,2,iglo)*vpac(-ntgrid:ntgrid-1,2,iglo)/sqrt(energy(ie))
       cdfac(ntgrid,:,iglo) = 0.0
       
       ! this is the first term multiplying dF_1/dE in Eq. 42 of MAB's GS2 notes
       ! i.e. Z_s * e * vpa . grad phi^{tb} d(F1/F0)/dE
       ! note there will be a correction to this added below
       ! because actual term appearing in GKE ~ (1/F0) * d(F1)/dE
       ! also note that vpadhdec is vpa*d(F1/F0)/dE at fixed mu (not xi)
       vparterm(-ntgrid:ntgrid-1,1,iglo) = spec(is)%zstm*tunits(ik)*code_dt &
            /delthet(-ntgrid:ntgrid-1) &
            * (abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid))) &
            * vpadhdec(-ntgrid:ntgrid-1,1,iglo)
       vparterm(-ntgrid:ntgrid-1,2,iglo) = spec(is)%zstm*tunits(ik)*code_dt &
            /delthet(-ntgrid:ntgrid-1) &
            * (abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid))) &
            * vpadhdec(-ntgrid:ntgrid-1,2,iglo)
       vparterm(ntgrid,:,iglo) = 0.0
       
       ! redefine vpar from vpa bhat dot grad theta to
       ! vpa bhat dot grad theta + v_Magnetic dot grad theta
       ! this accounts for the terms in Sec. 5.1 of MAB's GS2 notes
       ! TMP FOR TESTING -- MAB
!         vpar(-ntgrid:ntgrid-1,1,iglo) = vpar(-ntgrid:ntgrid-1,1,iglo) + &
!              0.5*rhostar*tunits(ik)*code_dt/delthet(-ntgrid:ntgrid-1) &
!              *(cvdrift_th(-ntgrid:ntgrid-1)*vpac(-ntgrid:ntgrid-1,1,iglo)**2 &
!              + gbdrift_th(-ntgrid:ntgrid-1)*0.5*(energy(ie)-vpac(-ntgrid:ntgrid-1,1,iglo)**2))
!         vpar(-ntgrid:ntgrid-1,2,iglo) = vpar(-ntgrid:ntgrid-1,2,iglo) + &
!              0.5*rhostar*tunits(ik)*code_dt/delthet(-ntgrid:ntgrid-1) &
!              *(cvdrift_th(-ntgrid:ntgrid-1)*vpac(-ntgrid:ntgrid-1,2,iglo)**2 &
!              + gbdrift_th(-ntgrid:ntgrid-1)*0.5*(energy(ie)-vpac(-ntgrid:ntgrid-1,2,iglo)**2))

       if (aky(ik) == 0) then
          wstar_neo(-ntgrid:ntgrid-1,it,ik) = -0.25*code_dt*akx(it) &
               * tmp8(-ntgrid:ntgrid-1)*(gds24_noq(-ntgrid:ntgrid-1)+gds24_noq(-ntgrid+1:ntgrid))
       else
          wstar_neo(-ntgrid:ntgrid-1,it,ik) = -code_dt*wunits(ik)*(tmp7(-ntgrid:ntgrid-1) &
               + 0.5*tmp8(-ntgrid:ntgrid-1)*(gds23(-ntgrid:ntgrid-1)+gds23(-ntgrid+1:ntgrid) &
               + theta0(it,ik)*shat*(gds24_noq(-ntgrid:ntgrid-1)+gds24_noq(-ntgrid+1:ntgrid))))
       end if
       wstar_neo(ntgrid,it,ik) = 0.0
    end do
    
    ! TMP FOR TESTING -- MAB
!    wdfac = 0. ; cdfac = 0. ; hneoc = 1. ; wstarfac = 0. ; wdttpfac = 0. ; vparterm = 0.
 
    ! vparterm is -2*vpar*(1+H^{neo}) - Ze*(vpa . grad phi^{tb})*(dH^{neo}/dE)
    ! note that vpar has contribution from v_{magnetic} . grad theta in it
    ! hneoc = 1 + H^{neo} below accounts for usual parallel streaming source term,
    ! as well as first of three terms multiplying F_1 in Eq. 42 of MAB's GS2 notes
    vparterm = -2.0*vpar*hneoc + vparterm
    
    ! hneoc below accounts for usual v_magnetic source term, 
    ! as well as second of three terms multiplying F_1 in Eq. 42 of MAB's GS2 notes
    ! wdfac on RHS deals with grad-B drift part of Eq. 60, as well as part of 
    ! curvature drift term in Eq. 54
    wdfac = wdfac + hneoc
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       wdfac(:,1,iglo) = wdfac(:,1,iglo)*wdrift(:,1,iglo) &
            + cdfac(:,1,iglo)*wcurv(:,iglo)
       wdfac(:,2,iglo) = wdfac(:,2,iglo)*wdrift(:,2,iglo) &
            + cdfac(:,2,iglo)*wcurv(:,iglo)
       wdttpfac(:,it,ik,ie,is,1) = wdfac(:,1,iglo)*wdriftttp(:,it,ik,ie,is,1) &
            + cdfac(:,1,iglo)*wcurv(:,iglo)
       wdttpfac(:,it,ik,ie,is,2) = wdfac(:,2,iglo)*wdriftttp(:,it,ik,ie,is,2) &
            + cdfac(:,2,iglo)*wcurv(:,iglo)

       ! add Z^2 * sqrt(m*T^3) * vpa * (bhat . grad(theta)) * d(phineo)/dth term to wdrift
       ! this modification of wdrift takes care of -g term in Eq. 67 of MB's gs2_notes
       ! arises from pulling 1/F0 inside dg/dE term
        wdrift(-ntgrid:ntgrid-1,1,iglo) = wdrift(-ntgrid:ntgrid-1,1,iglo) &
             - zi*code_dt*tunits(ik)*vpac(-ntgrid:ntgrid-1,1,iglo) &
             * spec(is)%zt*spec(is)%zstm*tmp8(-ntgrid:ntgrid-1) &
             * 0.5*(abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid)))
        wdrift(-ntgrid:ntgrid-1,2,iglo) = wdrift(-ntgrid:ntgrid-1,2,iglo) &
             - zi*code_dt*tunits(ik)*vpac(-ntgrid:ntgrid-1,2,iglo) &
             * spec(is)%zt*spec(is)%zstm*tmp8(-ntgrid:ntgrid-1) &
             * 0.5*(abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid)))
       
       ! hneoc below accounts for usual wstar term, as well as last of three terms
       ! multiplying F_1 in Eq. 42 of MAB'S GS2 notes
       ! note that hneo is necessary here because NEO's dist fn is
       ! normalized by F0(r) instead of F0 at center of simulation domain as in GS2
       wstarfac(:,:,iglo) = wstar(ik,ie,is)*hneoc(:,:,iglo) - wstarfac(:,:,iglo)
    end do

    deallocate (tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9)
    deallocate (vpadhdec,dhdec,dhdxic,cdfac,hneovpth)

  end subroutine init_lowflow
#endif

end module dist_fn
