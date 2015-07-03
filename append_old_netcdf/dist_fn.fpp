! Notes from BD, 7/2011:
!
! Need to extend the verr tools to include delta B_parallel integrals
! There are new factors of 1/B here and there which I do not understand.  

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

  private

  public :: init_dist_fn, finish_dist_fn
  
  !> Initializes a limited selection of arrays,
  !! for example g, gnew, vperp2, typically those which 
  !! are needed by other modules that don't need 
  !! dist_fn to be fully initialized (e.g. nonlinear_terms)
  !! This initialization level depends on grid sizes.
  public :: init_dist_fn_arrays
 
  !> Deallocates arrays allocated in init_dist_fn_arrays
  public :: finish_dist_fn_arrays

  !> Reads the dist_fn_knobs, source knobs, and 
  !! dist_fn_species_knobs namelists. This has be done independently
  !! of initializing the distribution function because
  !! it needs to be possible to override parameters such
  !! as g_exb.
  public :: init_dist_fn_parameters
  public :: finish_dist_fn_parameters

  !> Initializes parallel boundary conditions. This level
  !! depends on geometry.
  public :: init_dist_fn_level_1
  public :: finish_dist_fn_level_1

  !> Initializes bessel functions and field_eq. Note 
  !! that this level depends on species paramters.
  public :: init_dist_fn_level_2
  public :: finish_dist_fn_level_2

  !> Fully initialize the dist_fn module. Note that
  !! this level depends on the size of the timestep.
  public :: init_dist_fn_level_3
  public :: finish_dist_fn_level_3

  public :: set_overrides


  !!> init_vpar is called by init_dist_fn should only be called separately 
  !!! for testing purposes
  !public :: init_vpar
  public :: read_parameters, wnml_dist_fn, wnml_dist_fn_species, check_dist_fn
  public :: timeadv, exb_shear, g_exb, g_exbfac
  public :: g_exb_error_limit
  public :: g_exb_start_timestep, g_exb_start_time
  public :: init_bessel, init_fieldeq
  public :: getfieldeq, getan, getmoms, getemoms, getmoms_gryfx
  public :: getfieldeq_nogath
  public :: flux, lf_flux, eexchange
  public :: get_epar, get_heat
  public :: t0, omega0, gamma0, source0
  public :: reset_init, write_f, reset_physics, write_poly
  public :: M_class, N_class, i_class
  public :: l_links, r_links, itright, itleft, boundary
  public :: get_jext !GGH
  public :: get_verr, get_gtran, write_fyx, collision_error
  public :: get_init_field, getan_nogath
  public :: flux_vs_theta_vs_vpa
  public :: pflux_vs_theta_vs_vpa
  
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
  public :: pass_right, pass_left, init_pass_ends
  public :: init_enforce_parity, get_leftmost_it
  public :: gridfac1, awgt, gamtot3, fl_avg, apfac
  public :: adiabatic_option_switch, adiabatic_option_fieldlineavg

  ! knobs
  complex, dimension (:), allocatable :: fexp ! (nspec)
  real, dimension (:), allocatable :: bkdiff  ! (nspec)
  integer, dimension (:), allocatable :: bd_exp ! nspec
  real :: gridfac, apfac, driftknob, tpdriftknob, poisfac
  real :: t0, omega0, gamma0, source0
  real :: phi_ext, afilter, kfilter
  real :: wfb, g_exb, g_exbfac, omprimfac, btor_slab, mach
  real :: g_exb_start_time, g_exb_error_limit
  integer :: g_exb_start_timestep
  !logical :: dfexist, skexist, nonad_zero
  logical :: dfexist, skexist, nonad_zero, lf_default, lf_decompose, esv, opt_init_bc, opt_source
!CMR, 12/9/13: New logical cllc to modify order of operator in timeadv
!CMR, 21/5/14: New logical wfb_cmr to enforce trapping conditions on wfb
  logical :: cllc, wfb_cmr

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
  logical :: zero_forbid
  logical :: mult_imp, test
  logical :: accelerated_x = .false.
  logical :: accelerated_v = .false.
!Not sure why these are module level as only used in one routine
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

  complex, dimension (:,:,:), allocatable :: wdrift
  ! (-ntgrid:ntgrid, 2, -g-layout-)
#else
  real, dimension (:,:,:), allocatable :: wdrift
  ! (-ntgrid:ntgrid, 2, -g-layout-)
#endif

!  real, dimension (:,:,:,:,:), allocatable :: wdriftttp
  real, dimension (:,:,:,:,:,:), allocatable :: wdriftttp
  ! (-ntgrid:ntgrid,ntheta0,naky,negrid,nspec) replicated

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
  complex, dimension (:,:,:), allocatable :: gexp_1, gexp_2, gexp_3
  ! (-ntgrid:ntgrid,2, -g-layout-)

  ! momentum conservation
!  complex, dimension (:,:), allocatable :: g3int
!  real, dimension (:,:,:), allocatable :: sq

  ! exb shear
  integer, dimension(:), allocatable :: jump, ikx_indexed

  ! set_source
  real, dimension(:,:), allocatable :: ufac

  ! set_source_opt
  complex, dimension (:,:,:,:), allocatable :: source_coeffs

  ! getfieldeq1
  real, allocatable, dimension(:,:) :: awgt
  complex, allocatable, dimension(:,:) :: fl_avg !Changed

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

  type (connections_type), dimension (:), allocatable, save :: connections
  ! (-g-layout-)

  ! linked only
  type (redist_type), save :: gc_from_left, gc_from_right
  type (redist_type), save :: links_p, links_h
  type (redist_type), save :: wfb_p, wfb_h
  type (redist_type), save :: pass_right
  type (redist_type), save :: pass_left
  type (redist_type), save :: pass_wfb
  type (redist_type), save :: parity_redist

  integer, dimension (:,:), allocatable :: l_links, r_links
  integer, dimension (:,:,:), allocatable :: n_links
  logical, dimension (:,:), allocatable :: save_h
  logical :: no_comm = .false.
  integer, dimension(:), allocatable :: M_class, N_class
  integer :: i_class

  logical :: initialized = .false.

  logical :: exb_first = .true.

  logical :: initialized_dist_fn_parameters = .false.
  logical :: initialized_dist_fn_arrays = .false.
  logical :: initialized_dist_fn_level_1 = .false.
  logical :: initialized_dist_fn_level_2 = .false.
  logical :: initialized_dist_fn_level_3 = .false.
  !logical :: initializing = .true.
  logical :: readinit = .false.
  logical :: bessinit = .false.
  logical :: connectinit = .false.
  logical :: feqinit = .false.
  logical :: lpolinit = .false.
  logical :: fyxinit = .false.
  logical :: cerrinit = .false.
  logical :: mominit = .false.

  !The following variables are used if write_full_moments_notgc=T
  !or ginit_options='recon3'. These arrays are currently always initialised
  !even if they are not to be used. Can we improve this?
  real, allocatable :: mom_coeff(:,:,:,:)
  real, allocatable :: mom_coeff_npara(:,:,:), mom_coeff_nperp(:,:,:)
  real, allocatable :: mom_coeff_tpara(:,:,:), mom_coeff_tperp(:,:,:)
  real, allocatable :: mom_shift_para(:,:,:), mom_shift_perp(:,:,:)
  integer, parameter :: ncnt_mom_coeff=8


#ifdef NETCDF_PARALLEL
  logical, parameter :: moment_to_allprocs = .true.
#else
  logical, parameter :: moment_to_allprocs = .false.
#endif

contains

  subroutine check_dist_fn(report_unit)
    use kt_grids, only: grid_option, gridopt_box, gridopt_switch
    use nonlinear_terms, only: nonlinear_mode_switch, nonlinear_mode_on
    use species, only: spec, nspec, has_electron_species
    implicit none
    integer, intent(in) :: report_unit
    integer :: is 
    if (gridfac /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('You selected gridfac = ',e11.4,' in dist_fn_knobs.')") gridfac
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('The normal choice is gridfac = 1.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (apfac /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('You selected apfac = ',e11.4,' in dist_fn_knobs.')") apfac
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('The normal choice is apfac = 1.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (driftknob /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('You selected driftknob = ',e11.4,' in dist_fn_knobs.')") driftknob
       write (report_unit, fmt="('THIS IS EITHER AN ERROR, or you are DELIBERATELY SCALING THE DRIFTS.')") 
       write (report_unit, fmt="('The normal choice is driftknob = 1.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (tpdriftknob /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('You selected tpdriftknob = ',e11.4,' in dist_fn_knobs.')") tpdriftknob
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
          if (esv) then 
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('Single valued antot arrays will be enforced.')")
             write (report_unit, fmt="('This can significantly increase the cost of the run.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 
          endif
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
       write (report_unit, fmt="('Quasineutrality is not enforced.  The ratio (lambda_Debye/rho)**2 = ',e11.4)") poisfac
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
          write (report_unit, fmt="('Species ',i2,' has fexpi = ',e11.4)") is, aimag(fexp(is))
          write (report_unit, fmt="('THIS IS AN ERROR')")
          write (report_unit, fmt="('fexpi should be zero for all species.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if

       write (report_unit, fmt="('Species ',i2,' has fexpr = ', e11.4)") is, real(fexp(is))
    end do
    
  end subroutine check_dist_fn

  subroutine wnml_dist_fn(unit)
    use species, only: spec, has_electron_species
    implicit none
    integer, intent(in) :: unit
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
       write (unit, fmt="(' gridfac = ',e17.10)") gridfac
       write (unit, fmt="(' esv = ',L1)") esv
       write (unit, fmt="(' cllc = ',L1)") cllc
       write (unit, fmt="(' wfb_cmr = ',L1)") wfb_cmr
       write (unit, fmt="(' opt_init_bc = ',L1)") opt_init_bc
       write (unit, fmt="(' opt_source = ',L1)") opt_source

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

       if (apfac /= 1.) write (unit, fmt="(' apfac = ',e17.10)") apfac
       if (driftknob /= 1.) write (unit, fmt="(' driftknob = ',e17.10)") driftknob
       if (tpdriftknob /= 1.) write (unit, fmt="(' tpdriftknob = ',e17.10)") tpdriftknob
       if (poisfac /= 0.) write (unit, fmt="(' poisfac = ',e17.10)") poisfac
       if (kfilter /= 0.) write (unit, fmt="(' kfilter = ',e17.10)") kfilter
       if (afilter /= 0.) write (unit, fmt="(' afilter = ',e17.10)") afilter
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
          write (unit, fmt="(' source0 = ',e17.10)") source0
          write (unit, fmt="(' omega0 = ',e17.10)") omega0
          write (unit, fmt="(' gamma0 = ',e17.10)") gamma0
          write (unit, fmt="(' t0 = ',e17.10)") t0
          write (unit, fmt="(' phi_ext = ',e17.10)") phi_ext

       end select
       write (unit, fmt="(' /')")
    endif
  end subroutine wnml_dist_fn


  subroutine wnml_dist_fn_species(unit)
    use species, only: nspec
    implicit none
    integer, intent(in) :: unit
    integer :: i
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

  subroutine init_dist_fn_parameters
    use gs2_layouts, only: init_gs2_layouts
    use kt_grids, only: init_kt_grids
    use le_grids, only: init_le_grids
    use species, only: init_species
    use theta_grid, only: init_theta_grid
    logical,parameter:: debug=.false.

    if (initialized_dist_fn_parameters) return
    initialized_dist_fn_parameters  = .true.

    if (debug) write(6,*) "init_dist_fn: init_gs2_layouts"
    call init_gs2_layouts

    if (debug) write(6,*) "init_dist_fn: init_species"
    call init_species

    if (debug) write(6,*) "init_dist_fn: init_theta_grid"
    call init_theta_grid

    if (debug) write(6,*) "init_dist_fn: init_kt_grids"
    call init_kt_grids

    if (debug) write(6,*) "init_dist_fn: init_le_grids"
    call init_le_grids (accelerated_x, accelerated_v)

    if (debug) write(6,*) "init_dist_fn: read_parameters"
    call read_parameters

  end subroutine init_dist_fn_parameters

  subroutine init_dist_fn_arrays
    use mp, only: proc0, finish_mp, mp_abort
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid
    use gs2_layouts, only: init_dist_fn_layouts
    use nonlinear_terms, only: init_nonlinear_terms
    use run_parameters, only: init_run_parameters
    logical,parameter:: debug=.false.

    if (initialized_dist_fn_arrays) return
    initialized_dist_fn_arrays  = .true.

    call init_dist_fn_parameters

    if (test) then
       if (proc0) then
          write (*,*) 'nspecies = ',nspec
          write (*,*) 'nlambda = ', nlambda
          write (*,*) 'negrid = ',negrid
          write (*,*) 'ntheta0 = ',ntheta0
          write (*,*) 'naky = ',naky
       end if
       call finish_mp
       call mp_abort('only testing')
    end if

    if (debug) write(6,*) "init_dist_fn: run_parameters"
    call init_run_parameters

    if (debug) write(6,*) "init_dist_fn: dist_fn_layouts"
    call init_dist_fn_layouts (naky, ntheta0, nlambda, negrid, nspec)

    if (debug) write(6,*) "init_dist_fn: nonlinear_terms"
    call init_nonlinear_terms 

    if (debug) write(6,*) "init_dist_fn: allocate_arrays"
    call allocate_arrays


  end subroutine init_dist_fn_arrays

  subroutine init_dist_fn_level_1
    if (initialized_dist_fn_level_1) return
    initialized_dist_fn_level_1 = .true.

    call init_dist_fn_arrays

    !if (debug) write(6,*) "init_dist_fn: init_vperp2"
    call init_vperp2

    call init_dist_fn_arrays
    call init_bc
  
  end subroutine init_dist_fn_level_1

  subroutine init_dist_fn_level_2
    logical,parameter:: debug=.false.

    if (initialized_dist_fn_level_2) return
    initialized_dist_fn_level_2 = .true.

    call init_dist_fn_level_1

    if (debug) write(6,*) "init_dist_fn: init_bessel"
    call init_bessel

    if (debug) write(6,*) "init_dist_fn: init_fieldeq"
    call init_fieldeq
  end subroutine init_dist_fn_level_2

  subroutine init_dist_fn
    call init_dist_fn_level_3
  end subroutine init_dist_fn

  subroutine init_dist_fn_level_3
    use mp, only: proc0, finish_mp, mp_abort
    use species, only: init_species, nspec
    use collisions, only: init_collisions!, vnmult
    use hyper, only: init_hyper
    implicit none
    logical,parameter:: debug=.false.

    if (initialized_dist_fn_level_3) return
    initialized_dist_fn_level_3 = .true.

    if (initialized) return
    initialized = .true.

    call init_dist_fn_level_2

    if (debug) write(6,*) "init_dist_fn: init_vpar"
    call init_vpar

    if (debug) write(6,*) "init_dist_fn: init_wdrift"
    call init_wdrift

    if (debug) write(6,*) "init_dist_fn: init_wstar"
    call init_wstar

!    call init_vnmult (vnmult)
    if (debug) write(6,*) "init_dist_fn: init_collisions"
    call init_collisions ! needs to be after init_run_parameters

#ifdef LOWFLOW
    if (debug) write(6,*) "init_dist_fn: init_lowflow"
    call init_lowflow ! needs to be before init_invert_rhs but after
                      ! init_collisions as use_le_layouts is set in 
                      ! collisions
#endif

    if (debug) write(6,*) "init_dist_fn: init_invert_rhs"
    call init_invert_rhs

    if (debug) write(6,*) "init_dist_fn: init_hyper"
    call init_hyper

    if (debug) write(6,*) "init_dist_fn: init_mom_coeff"
    call init_mom_coeff

    if (debug) write(6,*) "init_dist_fn: init_source_term"
    call init_source_term

  end subroutine init_dist_fn_level_3

  subroutine finish_dist_fn
    call finish_dist_fn_parameters
  end subroutine finish_dist_fn

  subroutine finish_dist_fn_parameters
    if (.not. initialized_dist_fn_parameters) return
    initialized_dist_fn_parameters = .false.
    call finish_dist_fn_arrays
    readinit = .false.
    if (allocated(fexp)) deallocate (fexp, bkdiff, bd_exp)
  end subroutine finish_dist_fn_parameters

  subroutine finish_dist_fn_arrays
    use dist_fn_arrays, only: g, gnew, kx_shift, theta0_shift, vperp2

    if (.not. initialized_dist_fn_arrays) return
    initialized_dist_fn_arrays = .false.

    call finish_dist_fn_level_1
    if (allocated(g)) deallocate (g, gnew, g0)
    if (allocated(source_coeffs)) deallocate(source_coeffs)
    if (allocated(gexp_1)) deallocate (gexp_1, gexp_2, gexp_3)
    if (allocated(g_h)) deallocate (g_h, save_h)
    if (allocated(kx_shift)) deallocate (kx_shift)
    if (allocated(theta0_shift)) deallocate (theta0_shift)
    if (allocated(vperp2)) deallocate (vperp2)
    initialized_dist_fn_arrays = .false.
  end subroutine finish_dist_fn_arrays

  subroutine finish_dist_fn_level_1
    use redistribute, only: delete_redist
#ifdef LOWFLOW
    use lowflow, only: finish_lowflow_terms
    use dist_fn_arrays, only: hneoc, vparterm, wdfac, wstarfac, wdttpfac
#endif
    implicit none
    if (.not. initialized_dist_fn_level_1) return
    initialized_dist_fn_level_1 = .false.

    accelerated_x = .false. ; accelerated_v = .false.
    no_comm = .false.
    connectinit = .false.
    lpolinit = .false. ; fyxinit = .false. ; cerrinit = .false. ; mominit = .false.
    increase = .true. ; decrease = .false.
    exb_first = .true.

    call finish_dist_fn_level_2

    if (allocated(l_links)) deallocate (l_links, r_links, n_links)
    if (allocated(M_class)) deallocate (M_class, N_class)
    if (allocated(itleft)) deallocate (itleft, itright)
    if (allocated(connections)) deallocate (connections)
    if (allocated(g_adj)) deallocate (g_adj)
    if (allocated(jump)) deallocate (jump)
    if (allocated(ikx_indexed)) deallocate (ikx_indexed)
    if (allocated(ufac)) deallocate (ufac)
    if (allocated(fl_avg)) deallocate (fl_avg)
    if (allocated(awgt)) deallocate (awgt)
    if (allocated(kmax)) deallocate (kmax)
    if (allocated(mom_coeff)) deallocate (mom_coeff, mom_coeff_npara, mom_coeff_nperp, &
         mom_coeff_tpara, mom_coeff_tperp, mom_shift_para, mom_shift_perp)

    call delete_redist(gc_from_left)
    call delete_redist(gc_from_right)
    call delete_redist(links_p)
    call delete_redist(links_h)
    call delete_redist(wfb_p)
    call delete_redist(wfb_h)
    call delete_redist(pass_right)
    call delete_redist(pass_left)
    call delete_redist(pass_wfb)
    call delete_redist(parity_redist)

    ! gc_from_left, gc_from_right, links_p, links_h, wfb_p, wfb_h
#ifdef LOWFLOW
    call finish_lowflow_terms
    if(allocated(vparterm)) deallocate(vparterm)
    if(allocated(hneoc)) deallocate(hneoc)
    if(allocated(wdfac)) deallocate(wdfac)
    if(allocated(wstarfac)) deallocate(wstarfac)
    if(allocated(wdttpfac)) deallocate(wdttpfac)
    if(allocated(wstar_neo)) deallocate(wstar_neo)
    if(allocated(wcurv)) deallocate(wcurv)
#endif

  end subroutine finish_dist_fn_level_1

  subroutine finish_dist_fn_level_2
    use dist_fn_arrays, only: aj0, aj1

    if (.not. initialized_dist_fn_level_2) return
    initialized_dist_fn_level_2 = .false.

    call finish_dist_fn_level_3
    bessinit = .false. 
    feqinit = .false.  
    if (allocated(aj0)) deallocate (aj0, aj1)
    if (allocated(gridfac1)) deallocate (gridfac1, gamtot, gamtot1, gamtot2)
    if (allocated(gamtot3)) deallocate (gamtot3)
    if (allocated(a)) deallocate (a, b, r, ainv)
  end subroutine finish_dist_fn_level_2
  
  subroutine finish_dist_fn_level_3
    use dist_fn_arrays, only: gnew, g
    use dist_fn_arrays, only: vpa, vpac, vpar
    use dist_fn_arrays, only: ittp
    !initializing  = .true.
    if (.not. initialized_dist_fn_level_3) return
    initialized_dist_fn_level_3 = .false.

    initialized = .false.
    
    wdrift = 0.
    wdriftttp = 0.
    a = 0.
    b = 0.
    r = 0.
    ainv = 0.
    gnew = 0.
    g0 = 0.
    g = 0.

    if (allocated(vpa)) deallocate (vpa, vpac, vpar)
    if (allocated(ittp)) deallocate (ittp, wdrift, wdriftttp)
    if (allocated(wstar)) deallocate (wstar)
  end subroutine finish_dist_fn_level_3

  subroutine set_overrides(prof_ov)
    use overrides, only: profiles_overrides_type
    type(profiles_overrides_type), intent(in) :: prof_ov
    if (prof_ov%override_g_exb) g_exb = prof_ov%g_exb
    if (prof_ov%override_mach) mach = prof_ov%mach
  end subroutine set_overrides


  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use theta_grid, only: shat
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
         g_exb, g_exbfac, omprimfac, btor_slab, mach, cllc, lf_default, &
         g_exb_start_time, g_exb_start_timestep, g_exb_error_limit, &
         lf_decompose, esv, wfb_cmr, opt_init_bc, opt_source, zero_forbid
    
    namelist /source_knobs/ t0, omega0, gamma0, source0, phi_ext, source_option
    integer :: ierr, is, in_file
    real :: bd

    if (readinit) return
    readinit = .true.

    if (proc0) then
       boundary_option = 'default'
       nonad_zero = .true.  ! BD: Default value changed to TRUE  8.15.13
       wfb_cmr= .false.
       cllc = .false.
       esv = .false.
       opt_init_bc = .false.
       opt_source = .false.
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
       g_exb_error_limit = 0.0
       g_exb_start_time = -1
       g_exb_start_timestep = -1
       mach = 0.0
       omprimfac = 1.0
       btor_slab = 0.0
       wfb = 1.
       mult_imp = .false.
       test = .false.
       def_parity = .false.
       zero_forbid = .false.
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


       ierr = error_unit()
       call get_option_value &
            (boundary_option, boundaryopts, boundary_option_switch, &
            ierr, "boundary_option in dist_fn_knobs",.true.)

       if(abs(shat) <=  1.e-5) boundary_option = 'periodic'

       call get_option_value &
            (source_option, sourceopts, source_option_switch, &
            ierr, "source_option in source_knobs",.true.)
       call get_option_value &
            (adiabatic_option, adiabaticopts, adiabatic_option_switch, &
            ierr, "adiabatic_option in dist_fn_knobs",.true.)

    end if
    if (.not.allocated(fexp)) allocate (fexp(nspec), bkdiff(nspec), bd_exp(nspec))
    if (proc0) call read_species_knobs

    call broadcast (boundary_option_switch)
    call broadcast (nonad_zero)
    call broadcast (wfb_cmr)
    call broadcast (cllc)
    call broadcast (esv)
    call broadcast (opt_init_bc)
    call broadcast (opt_source)
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
    call broadcast (g_exb_start_timestep)
    call broadcast (g_exb_start_time)
    call broadcast (g_exb_error_limit)
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
    call broadcast (zero_forbid)
    call broadcast (lf_default)
    call broadcast (lf_decompose)
    call broadcast (even)
    call broadcast (wfb)

    !<DD>Turn off esv if not using linked boundaries
    esv=esv.and.(boundary_option_switch.eq.boundary_option_linked)

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
    integer :: bd_exp,iostat
    real :: fexpr, fexpi, bakdif
    namelist /dist_fn_species_knobs/ fexpr, fexpi, bakdif, bd_exp

    fexpr = real(fexp_out)
    fexpi = aimag(fexp_out)
    bakdif = bakdif_out
    bd_exp = bd_exp_out
    read (unit=unit, nml=dist_fn_species_knobs,iostat=iostat)
    if(iostat /= 0) write(6,*) 'Error ',iostat,'dist_fn species knobs'
    fexp_out = cmplx(fexpr,fexpi)
    bd_exp_out = bd_exp
    bakdif_out = bakdif
  end subroutine fill_species_knobs

  subroutine init_wdrift
    use species, only: nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use le_grids, only: negrid, ng2, nlambda, jend, forbid
!    use le_grids, only: al
!    use theta_grids, only: bmag
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use dist_fn_arrays, only: ittp
    implicit none
    integer :: ig, ik, it, il, ie, is
    integer :: iglo
!    logical :: alloc = .true.
!CMR
    logical,parameter :: debug = .false.
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
       allocate (wdriftttp(-ntgrid:ntgrid,ntheta0,naky,negrid,nspec,2))
    end if
    wdrift = 0.  ; wdriftttp = 0.
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
             wdrift(ig,1,iglo) = wdrift_func(ig, il, ie, it, ik)*tpdriftknob
#ifdef LOWFLOW
             wcurv(ig,iglo) = wcurv_func(ig, it, ik)*tpdriftknob
#endif
!CMR:  multiply trapped particle drift by tpdriftknob 
!CMR:               (tpdriftknob defaults to driftknob if not supplied)
          else
             wdrift(ig,1,iglo) = wdrift_func(ig, il, ie, it, ik)*driftknob
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
                           = wdrift_func(ig,ittp(ig),ie,it,ik)*tpdriftknob
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

!    alloc = .false.
!CMR
    if (debug) write(6,*) 'init_wdrift: driftknob, tpdriftknob=',driftknob,tpdriftknob

  end subroutine init_wdrift

  function wdrift_func (ig, il, ie, it, ik)
    use theta_grid, only: bmag, gbdrift, gbdrift0, cvdrift, cvdrift0
    use theta_grid, only: shat
    use kt_grids, only: aky, theta0, akx
    use le_grids, only: energy, al
    use run_parameters, only: wunits
    use gs2_time, only: code_dt
    implicit none
    real :: wdrift_func
    integer, intent (in) :: ig, ik, it, il, ie

    ! note that wunits=aky/2 (for wstar_units=F)
    if (aky(ik) == 0.0) then
       wdrift_func = akx(it)/shat &
                    *(cvdrift0(ig)*energy(ie)*(1.0 - al(il)*bmag(ig)) &
                      + gbdrift0(ig)*0.5*energy(ie)*al(il)*bmag(ig)) &
                     *code_dt/2.0
    else
       wdrift_func = ((cvdrift(ig) + theta0(it,ik)*cvdrift0(ig)) &
                        *energy(ie)*(1.0 - al(il)*bmag(ig)) &
                      + (gbdrift(ig) + theta0(it,ik)*gbdrift0(ig)) &
                        *0.5*energy(ie)*al(il)*bmag(ig)) &
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
       wcoriolis_func = mach*sqrt(max(energy(ie)*(1.0-al(il)*bmag(ig)),0.0)) &
            * cdrift0(ig) * code_dt * akx(it)/(2.*shat*spec(is)%stm)
    else
       wcoriolis_func = mach*sqrt(max(energy(ie)*(1.0-al(il)*bmag(ig)),0.0)) &
            * (cdrift(ig) + theta0(it,ik)*cdrift0(ig))*code_dt*wunits(ik)/spec(is)%stm
    end if

  end function wcoriolis_func

  subroutine init_vperp2
    use theta_grid, only: ntgrid, bmag
    use gs2_layouts, only: g_lo, ik_idx, il_idx, ie_idx, is_idx
    use dist_fn_arrays, only: vperp2
    use le_grids, only: energy, al
    implicit none
    real :: al1, e1
    integer :: iglo
    if (.not.allocated(vperp2)) then
       allocate (vperp2 (-ntgrid:ntgrid,  g_lo%llim_proc:g_lo%ulim_alloc))
    endif
    vperp2 = 0.
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       al1 = al(il_idx(g_lo,iglo))
       e1 = energy(ie_idx(g_lo,iglo))

       vperp2(:,iglo) = bmag*al1*e1
    end do


  end subroutine init_vperp2

  subroutine init_vpar
    use dist_fn_arrays, only: vpa, vpar, vpac
    use species, only: spec
    use theta_grid, only: ntgrid, delthet, bmag, gradpar
    use le_grids, only: energy, al
    use run_parameters, only: tunits
    use gs2_time, only: code_dt
    use gs2_layouts, only: g_lo, ik_idx, il_idx, ie_idx, is_idx
    implicit none
    integer :: iglo, ik, is
    real :: al1, e1
    
    if (.not.allocated(vpa)) then
       allocate (vpa    (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (vpac   (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       ! EGH Put vperp2 in a separate function
       !allocate (vperp2 (-ntgrid:ntgrid,  g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (vpar   (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    endif
    vpa = 0. ; vpac = 0.; vpar = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       al1 = al(il_idx(g_lo,iglo))
       e1 = energy(ie_idx(g_lo,iglo))

       vpa(:,1,iglo) = sqrt(e1*max(0.0, 1.0 - al1*bmag))
       vpa(:,2,iglo) = - vpa(:,1,iglo)
       !vperp2(:,iglo) = bmag*al1*e1

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
    end do

  end subroutine init_vpar

  subroutine init_wstar
    use species, only: nspec, spec
    use kt_grids, only: naky
    use le_grids, only: negrid, energy
    use run_parameters, only: wunits
    use gs2_time, only: code_dt

    implicit none
    integer :: ik, ie, is

    if(.not.allocated(wstar)) allocate (wstar(naky,negrid,nspec))

    do is = 1, nspec
       do ie = 1, negrid
          do ik = 1, naky
             wstar(ik,ie,is) = code_dt*wunits(ik) &
                  *(spec(is)%fprim+spec(is)%tprim*(energy(ie)-1.5))
          end do
       end do
    end do
  end subroutine init_wstar

  subroutine init_bessel
    use dist_fn_arrays, only: aj0, aj1
    use kt_grids, only: kperp2
    use species, only: spec
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: energy, al
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use spfunc, only: j0, j1

    implicit none

    integer :: ig, ik, it, il, ie, is
    integer :: iglo
    real :: arg

    if (bessinit) return
    bessinit = .true.

    !write(*,*) 'Allocating aj0', ntgrid, g_lo%llim_proc, g_lo%ulim_alloc

    allocate (aj0(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (aj1(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do ig = -ntgrid, ntgrid
          arg = spec(is)%bess_fac*spec(is)%smz*sqrt(energy(ie)*al(il)/bmag(ig)*kperp2(ig,it,ik))
          aj0(ig,iglo) = j0(arg)
! CMR 17/1/06: BEWARE, j1 returns and aj1 stores J_1(x)/x (NOT J_1(x)), 
          aj1(ig,iglo) = j1(arg)
       end do
    end do
  end subroutine init_bessel

  subroutine init_invert_rhs
    use dist_fn_arrays, only: vpar, ittp
    use species, only: spec, nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, ng2, forbid, negrid
    use constants, only: zi
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    implicit none
    integer :: iglo
    integer :: ig, ik, it, il, ie, is, isgn
    real :: wdttp, vp, bd
#ifdef LOWFLOW
    complex :: wd
#else
    real :: wd
#endif
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
             
! CMR, DD, 25/7/2014: 
!  set ainv=1 just left of lower bounce points ONLY for RIGHTWARDS travelling 
!  trapped particles. Part of multiple trapped particle algorithm
!  NB not applicable to ttp or wfb!
if( isgn.eq.1 )then
                if (forbid(ig,il) .and. .not. forbid(ig+1,il)) then
                   ainv(ig,isgn,iglo) = 1.0 + ainv(ig,isgn,iglo)
                end if
endif
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


  end subroutine init_invert_rhs

  subroutine init_bc
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    select case (boundary_option_switch)
    case (boundary_option_linked)
       if(opt_init_bc)then
          call init_connected_bc_opt
       else
          call init_connected_bc
       endif
       if(def_parity)call init_enforce_parity(parity_redist)
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
  end subroutine init_bc

  subroutine init_source_term
    use dist_fn_arrays, only: vpar, vpac, aj0
    use run_parameters, only: wunits, fapar
    use gs2_time, only: code_dt
    use species, only: spec,nspec
    use hyper, only: D_res
    use theta_grid, only: ntgrid, itor_over_b
    use le_grids, only: anon, negrid, energy
    use constants, only: zi,pi
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
#ifdef LOWFLOW
    use dist_fn_arrays, only: vparterm, wdfac, wstarfac
#endif
    implicit none
    integer :: iglo
    integer :: ig, ik, it, il, ie, is, isgn

    !Initialise ufac for use in set_source
    if (.not. allocated(ufac)) then
       allocate (ufac(negrid, nspec))
       do ie = 1, negrid
          do is = 1, nspec
             ufac(ie, is) = (2.0*spec(is)%uprim &
                  + spec(is)%uprim2*energy(ie)**(1.5)*sqrt(pi)/4.0)
          end do
       end do
    endif

    if(.not.opt_source) return

    !Setup the source coefficients
    !See comments in get_source_term and set_source
    !for information on these terms.
    do iglo=g_lo%llim_proc,g_lo%ulim_proc
       ie=ie_idx(g_lo,iglo)
       il=il_idx(g_lo,iglo)
       ik=ik_idx(g_lo,iglo)
       it=it_idx(g_lo,iglo)
       is=is_idx(g_lo,iglo)
       do isgn=1,2
          do ig=-ntgrid,ntgrid-1
#ifdef LOWFLOW
             !Phi m
             source_coeffs(1,ig,isgn,iglo)=-2*anon(ie)*vparterm(ig,isgn,iglo)
             !Phi p
             source_coeffs(2,ig,isgn,iglo)=-zi*anon(ie)*wdfac(ig,isgn,iglo) &
              + zi*(wstarfac(ig,isgn,iglo) &
              + vpac(ig,isgn,iglo)*code_dt*wunits(ik)*ufac(ie,is) &
              -2.0*omprimfac*vpac(ig,isgn,iglo)*code_dt*wunits(ik)*g_exb*itor_over_B(ig)/spec(is)%stm)
             if(fapar.gt.0)then
                !Apar m
                source_coeffs(3,ig,isgn,iglo)=-spec(is)%zstm*anon(ie)*&
                     vpac(ig,isgn,iglo)*(aj0(ig+1,iglo)+aj0(ig,iglo))*0.5
                !Apar p
                source_coeffs(4,ig,isgn,iglo)=anon(ie)*D_res(it,ik)*spec(is)%zstm*&
                     vpac(ig,isgn,iglo)-spec(is)%stm*vpac(ig,isgn,iglo)*&
                     zi*(wstarfac(ig,isgn,iglo) &
                     + vpac(ig,isgn,iglo)*code_dt*wunits(ik)*ufac(ie,is) &
                     -2.0*omprimfac*vpac(ig,isgn,iglo)*code_dt*wunits(ik)*g_exb*itor_over_B(ig)/spec(is)%stm) 
             endif
#else

             !Phi m
             source_coeffs(1,ig,isgn,iglo)=-2*anon(ie)*vpar(ig,isgn,iglo)
             !Phi p
             source_coeffs(2,ig,isgn,iglo)=-zi*anon(ie)*wdrift(ig,isgn,iglo) &
              + zi*(wstar(ik,ie,is) &
              + vpac(ig,isgn,iglo)*code_dt*wunits(ik)*ufac(ie,is) &
              -2.0*omprimfac*vpac(ig,isgn,iglo)*code_dt*wunits(ik)*g_exb*itor_over_B(ig)/spec(is)%stm)
             if(fapar.gt.0)then
                !Apar m
                source_coeffs(3,ig,isgn,iglo)=-spec(is)%zstm*anon(ie)*&
                     vpac(ig,isgn,iglo)*(aj0(ig+1,iglo)+aj0(ig,iglo))*0.5
                !Apar p
                source_coeffs(4,ig,isgn,iglo)=anon(ie)*D_res(it,ik)*spec(is)%zstm*&
                     vpac(ig,isgn,iglo)-spec(is)%stm*vpac(ig,isgn,iglo)*&
                     zi*(wstar(ik,ie,is) &
                     + vpac(ig,isgn,iglo)*code_dt*wunits(ik)*ufac(ie,is) &
                     -2.0*omprimfac*vpac(ig,isgn,iglo)*code_dt*wunits(ik)*g_exb*itor_over_B(ig)/spec(is)%stm) 
             endif
#endif
          enddo
       enddo
    enddo
  end subroutine init_source_term

  subroutine init_connected_bc
    use theta_grid, only: ntgrid, nperiod, ntheta, theta
    use kt_grids, only: naky, ntheta0, aky, theta0
    use le_grids, only: ng2, nlambda
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use gs2_layouts, only: idx, proc_id
    use mp, only: iproc, nproc, max_allreduce, mp_abort
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
          if ( wfb_cmr ) then
! if wfb_cmr = .t.  do NOT save homogeneous solutions for wfb
             if (nlambda > ng2 .and. il == ng2+1) save_h(:,iglo) = .false.
          else
             if (nlambda > ng2               .and. &
                  il == ng2+1                .and. &
                  connections(iglo)%neighbor .and. &
                  aky(ik) /= 0.0) &
                  save_h (:,iglo) = .true.
          endif
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
                      call mp_abort('l_links error')
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
                      call mp_abort('r_links error',.true.)
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
!CMR, 20/10/2013:
!     Communicate pattern sends passing particles to downwind linked cells
!      (ntgrid,1,iglo)  => each RH connected cell (j,1,iglo_right)
!      (-ntgrid,2,iglo) => each LH connected cell (j,2,iglo_left)
!          where j index in receive buffer = #hops in connection
!                         (nb j=1 is nearest connection)

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
    if ( wfb_cmr ) then
       call init_pass_wfb
    else

       nn_to = 0
       nn_from = 0 

       do iglo = g_lo%llim_world, g_lo%ulim_world
!CMR, 20/10/2013:
!     Communicate pattern sends wfb endpoint to ALL linked cells
!      (ntgrid,1,iglo)  => ALL connected cells (j,1,iglo_conn)
!            where j index in receive buffer = r_links(iglo)+1
!      (-ntgrid,2,iglo) => ALL connected cells (j,2,iglo_conn)
!         where j index in receive buffer = l_links(iglo)+1
!

          il = il_idx(g_lo,iglo)
          if (il /= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle
             
          iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

! v_par > 0:
!CMR: introduced iglo_star to make notation below less confusing
!
          call find_leftmost_link (iglo, iglo_star, ipright)
          do j = 1, ncell
! sender
             if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
! receiver
             if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1
             call get_right_connection (iglo_star, iglo_right, ipright)
             iglo_star=iglo_right
          end do
             
! v_par < 0:
          call find_rightmost_link (iglo, iglo_star, ipleft)
          do j = 1, ncell
! sender
             if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
! receiver
             if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
             call get_left_connection (iglo_star, iglo_left, ipleft)
             iglo_star=iglo_left
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
          call find_leftmost_link (iglo, iglo_star, ipright)
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
                to(ip)%third(n) = iglo_star
             end if
             call get_right_connection (iglo_star, iglo_right, ipright)
             iglo_star=iglo_right
          end do
             
! v_par < 0: 
          call find_rightmost_link (iglo, iglo_star, ipleft)
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
                to(ip)%third(n) = iglo_star
             end if
             call get_left_connection (iglo_star, iglo_left, ipleft)
             iglo_star=iglo_left
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
    endif
       
! n_links_max is typically 2 * number of cells in largest supercell
       allocate (g_adj (n_links_max, 2, g_lo%llim_proc:g_lo%ulim_alloc))

! now set up links_h:
! excluding wfb

       nn_to = 0
       nn_from = 0 

       do iglo = g_lo%llim_world, g_lo%ulim_world
!CMR, 20/10/2013:
!     Communicate pattern sends passing homogeneous endpoints to linked cells
!      (ntgrid,1,iglo) (cell > 2) => RH connected cells (j,1,iglo_conn)
!            where j index in receive buffer = 2.l_links(iglo)+ #RH hops
!      (-ntgrid,2,iglo) (cell < ncell-1) => LH connected cells (j,2,iglo_conn)
!         where j index in receive buffer = 2.r_links(iglo)+#LH hops
! NB: homogenous(H) follow inhomogeneous (IH), but orders in j differs!
!     IH#hop=1, IH#hop=2, ... H#hop=2, H#hop=1  
!          

          il = il_idx(g_lo,iglo)
          if (nlambda > ng2 .and. il >= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

! If this is not the first link in the chain, continue
          if (l_links(ik, it) > 0) then
!CMR, 20/10/2013:  
! this sends homogeneous solutions for rightwards passing particles 
! rightwards from all cells APART from the INCOMING cell 
!
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
!CMR, 20/10/2013:  
! this sends homogeneous solutions for leftwards passing particles 
! leftwards from all cells APART from the INCOMING cell 
!
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
    if ( .not. wfb_cmr ) then
       nn_to = 0
       nn_from = 0 

       do iglo = g_lo%llim_world, g_lo%ulim_world
!CMR, 20/10/2013:
!     Communicate pattern sends wfb homogeneous endpoint to ALL linked cells
!      (ntgrid,1,iglo)  => ALL connected cells (j,1,iglo_conn)
!            where j index in receive buffer = 2.ncell - r_links(iglo)
!            1 <= j <= 2.ncell  accommodates inhomog, homog bcs from cells
!               ncell, ncell-1, ncell-2, ..., 1, 1, 2,... ncell-1, ncell  
!      (-ntgrid,2,iglo) => ALL connected cells (j,2,iglo_conn)
!         where j index in receive buffer = 2.ncell - l_links(iglo)
!         1 <= j <= 2.ncell  accommodates inhomog, homog bcs from cells
!                1, 2, .., ncell, ncell, ncell-1, .., 2, 1  
! NB: homogenous follow inhomogeneous bcs, but order of source cell differs!

          il = il_idx(g_lo,iglo)
          if (il /= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

          iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

! v_par > 0:
          call find_leftmost_link (iglo, iglo_star, ipright)
          do j = 1, ncell
! sender
             if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
! receiver
             if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1
             call get_right_connection (iglo_star, iglo_right, ipright)
             iglo_star=iglo_right
          end do

! v_par < 0:
          call find_rightmost_link (iglo, iglo_star, ipleft)
          do j = 1, ncell
! sender
             if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
! receiver
             if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
             call get_left_connection (iglo_star, iglo_left, ipleft)
             iglo_star=iglo_left
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
          call find_leftmost_link (iglo, iglo_star, ipright)
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
                to(ip)%third(n) = iglo_star
             end if
             call get_right_connection (iglo_star, iglo_right, ipright)
             iglo_star=iglo_right
          end do
 
! v_par < 0:
          call find_rightmost_link (iglo, iglo_star, ipleft)
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
                to(ip)%third(n) = iglo_star
             end if
             call get_left_connection (iglo_star, iglo_left, ipleft)
             iglo_star=iglo_left
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
    endif

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
          call mp_abort('problem with connnected bc')
       end if

    end if

  end subroutine init_connected_bc

  subroutine init_connected_bc_opt
    use theta_grid, only: ntgrid, nperiod, ntheta, theta
    use kt_grids, only: naky, ntheta0, aky, theta0
    use le_grids, only: ng2, nlambda
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use gs2_layouts, only: idx, proc_id
    use mp, only: iproc, nproc, max_allreduce, mp_abort
    use redistribute, only: index_list_type, init_fill, delete_list
    implicit none
    type (index_list_type), dimension(0:nproc-1) :: from, from_p, from_h, to_p, to_h
    integer, dimension (0:nproc-1) :: nn_from, nn_to, nn_from_h, nn_to_h
    integer, dimension (3) :: to_low, from_low, to_high, from_high
    integer :: ik, it, il, ie, is, iglo, it0, itl, itr, jshift0
    integer :: ip, ipleft, ipright
    integer :: iglo_left, iglo_right, i, j, k
    integer :: iglo_star, it_star, ncell
    integer :: n, n_h, n_links_max, nn_max
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
    else if (naky == 1 .and. ntheta0 > 1 .and. aky(1) /= 0.0) then
       jshift0 = int((theta(ng)-theta(-ng))/(theta0(2,1)-theta0(1,1)) + 0.01)
    else
       jshift0 = 1
    end if
    
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
       il = il_idx(g_lo,iglo)

       ! get processors and indices for j' (kx') modes connecting
       ! to mode j (kx) so we can set up communication -- MAB

! if non-wfb trapped particle, no connections
       if (nlambda > ng2 .and. il >= ng2+2) then
          connections(iglo)%iproc_left = -1
          connections(iglo)%iglo_left = -1
          connections(iglo)%iproc_right = -1
          connections(iglo)%iglo_right = -1
       else
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
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
                      call mp_abort('l_links error',.true.)
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
                      call mp_abort('r_links error',.true.)
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

!################################
!Setup the linked communications
!################################
       !<DD>Note the communications setup here are often equivalent to an all-to-all type
       !communication, i.e. when nproc --> 2 nproc, t_fill --> 4 t_fill
       !See comments in invert_rhs_linked for more details.

       !Note: This setup currently involves several loops over the entire domain
       !and hence does not scale well (effectively serial code). This can come to
       !dominate initialisation at large core count.

       nn_to = 0
       nn_from = 0 
       nn_to_h = 0
       nn_from_h = 0 

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          ! Exclude trapped particles (inc wfb)
          if (nlambda > ng2 .and. il >= ng2+1) cycle

          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

          ip = proc_id(g_lo,iglo)

          iglo_star = iglo
          do j = 1, r_links(ik, it)
             call get_right_connection (iglo_star, iglo_right, ipright)
! sender
             if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
! receiver
             if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1

             !Special counting for links_h
             if(l_links(ik,it)>0) then
                if (ip == iproc) nn_from_h(ipright) = nn_from_h(ipright) + 1
                if (ipright == iproc) nn_to_h(ip) = nn_to_h(ip) + 1
             endif

             iglo_star = iglo_right
          end do
             
          iglo_star = iglo
          do j = 1, l_links(ik, it)
             call get_left_connection (iglo_star, iglo_left, ipleft)
! sender
             if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
! receiver
             if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1

             !Special counting for links_h
             if(r_links(ik,it)>0) then
                if (ip == iproc) nn_from_h(ipleft) = nn_from_h(ipleft) + 1
                if (ipleft == iproc) nn_to_h(ip) = nn_to_h(ip) + 1
             endif

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
             allocate (from_p(ip)%first(nn_from(ip)))
             allocate (from_p(ip)%second(nn_from(ip)))
             allocate (from_p(ip)%third(nn_from(ip)))
          endif
          if (nn_from_h(ip) > 0) then
             allocate (from_h(ip)%first(nn_from_h(ip)))
             allocate (from_h(ip)%second(nn_from_h(ip)))
             allocate (from_h(ip)%third(nn_from_h(ip)))
          endif
          if (nn_to(ip) > 0) then
             allocate (to_p(ip)%first(nn_to(ip)))
             allocate (to_p(ip)%second(nn_to(ip)))
             allocate (to_p(ip)%third(nn_to(ip)))
          endif
          if (nn_to_h(ip)>0) then
             allocate (to_h(ip)%first(nn_to_h(ip)))
             allocate (to_h(ip)%second(nn_to_h(ip)))
             allocate (to_h(ip)%third(nn_to_h(ip)))
          endif
       end do
       
       nn_from = 0
       nn_from_h=0
       nn_to = 0          
       nn_to_h=0

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (nlambda > ng2 .and. il >= ng2+1) cycle

          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

          ip = proc_id(g_lo,iglo)          

          iglo_star = iglo
          do j = 1, r_links(ik, it)
             call get_right_connection (iglo_star, iglo_right, ipright)
! sender
             if (ip == iproc) then
                n = nn_from(ipright) + 1
                nn_from(ipright) = n
                from_p(ipright)%first(n) = ntgrid
                from_p(ipright)%second(n) = 1
                from_p(ipright)%third(n) = iglo
                !Special restriction for links_h
                if(l_links(ik,it)>0)then
                   n_h=nn_from_h(ipright)+1
                   nn_from_h(ipright)=n_h
                   from_h(ipright)%first(n_h) = ntgrid
                   from_h(ipright)%second(n_h) = 1
                   from_h(ipright)%third(n_h) = iglo
                endif
             end if
! receiver
             if (ipright == iproc) then
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to_p(ip)%first(n) = j 
                to_p(ip)%second(n) = 1
                to_p(ip)%third(n) = iglo_right
                !Special restriction for links_h
                if(l_links(ik,it)>0)then
                   n_h=nn_to_h(ip)+1
                   nn_to_h(ip)=n_h
                   to_h(ip)%first(n_h) = 2*l_links(ik,it)+j
                   to_h(ip)%second(n_h) = 1
                   to_h(ip)%third(n_h) = iglo_right
                endif
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
                from_p(ipleft)%first(n) = -ntgrid
                from_p(ipleft)%second(n) = 2
                from_p(ipleft)%third(n) = iglo
                !Special restriction for links_h
                if(r_links(ik,it)>0)then
                   n_h=nn_from_h(ipleft)+1
                   nn_from_h(ipleft)=n_h
                   from_h(ipleft)%first(n_h) = -ntgrid
                   from_h(ipleft)%second(n_h) = 2
                   from_h(ipleft)%third(n_h) = iglo
                endif
             end if
! receiver
             if (ipleft == iproc) then
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to_p(ip)%first(n) = j 
                to_p(ip)%second(n) = 2
                to_p(ip)%third(n) = iglo_left
                !Special restriction for links_h
                if(r_links(ik,it)>0)then
                   n_h = nn_to_h(ip) + 1
                   nn_to_h(ip) = n_h
                   to_h(ip)%first(n_h) = 2*r_links(ik,it)+j 
                   to_h(ip)%second(n_h) = 2
                   to_h(ip)%third(n_h) = iglo_left
                endif
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

       call init_fill (links_p, 'c', to_low, to_high, to_p, &
            from_low, from_high, from_p)
       call init_fill (links_h, 'c', to_low, to_high, to_h, &
            from_low, from_high, from_h)
       
       call delete_list (from_p)
       call delete_list (from_h)
       call delete_list (to_p)
       call delete_list (to_h)

! take care of wfb

       nn_to = 0
       nn_from = 0 

       !NOTE: No special restriction/counting for wfb_h unlike links_h
       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (il /= ng2+1) cycle

          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

          ip = proc_id(g_lo,iglo)
             
          iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

! v_par > 0:
          call find_leftmost_link (iglo, iglo_star, ipright)
          do j = 1, ncell
! sender
             if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
! receiver
             if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1
             call get_right_connection (iglo_star, iglo_right, ipright)
             iglo_star=iglo_right
          end do
             
! v_par < 0:
          call find_rightmost_link (iglo, iglo_star, ipleft)
          do j = 1, ncell
! sender
             if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
! receiver
             if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
             call get_left_connection (iglo_star, iglo_left, ipleft)
             iglo_star=iglo_left
          end do
       end do
       
       do ip = 0, nproc-1
          if (nn_from(ip) > 0) then
             allocate (from(ip)%first(nn_from(ip)))
             allocate (from(ip)%second(nn_from(ip)))
             allocate (from(ip)%third(nn_from(ip)))
          endif
          if (nn_to(ip) > 0) then
             allocate (to_p(ip)%first(nn_to(ip)))
             allocate (to_p(ip)%second(nn_to(ip)))
             allocate (to_p(ip)%third(nn_to(ip)))
             allocate (to_h(ip)%first(nn_to(ip)))
             allocate (to_h(ip)%second(nn_to(ip)))
             allocate (to_h(ip)%third(nn_to(ip)))
          endif
       end do
       
       nn_from = 0
       nn_to = 0          

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (il /= ng2+1) cycle

          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

          ip = proc_id(g_lo,iglo)

          iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

! v_par > 0: 
          call find_leftmost_link (iglo, iglo_star, ipright)
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
                to_p(ip)%first(n) = r_links(ik, it) + 1
                to_p(ip)%second(n) = 1
                to_p(ip)%third(n) = iglo_star
                to_h(ip)%first(n) = 2*ncell-r_links(ik,it)
                to_h(ip)%second(n) = 1
                to_h(ip)%third(n) = iglo_star
             end if
             call get_right_connection (iglo_star, iglo_right, ipright)
             iglo_star=iglo_right
          end do
             
! v_par < 0: 
          call find_rightmost_link (iglo, iglo_star, ipleft)
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
                to_p(ip)%first(n) = l_links(ik, it) + 1
                to_p(ip)%second(n) = 2
                to_p(ip)%third(n) = iglo_star
                to_h(ip)%first(n) = 2*ncell-l_links(ik, it)
                to_h(ip)%second(n) = 2
                to_h(ip)%third(n) = iglo_star
             end if
             call get_left_connection (iglo_star, iglo_left, ipleft)
             iglo_star=iglo_left
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

       call init_fill (wfb_p, 'c', to_low, to_high, to_p, from_low, from_high, from)
       call init_fill (wfb_h, 'c', to_low, to_high, to_h, from_low, from_high, from)
       
       call delete_list (from)
       call delete_list (to_p)
       call delete_list (to_h)
      
! n_links_max is typically 2 * number of cells in largest supercell
       allocate (g_adj (n_links_max, 2, g_lo%llim_proc:g_lo%ulim_alloc))
!################################
!End of linked comms setup
!################################

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
          call mp_abort('problem with connnected bc')
       end if

    end if

  end subroutine init_connected_bc_opt

  subroutine init_pass_right
    !Create a pass_right object to pass boundary values to the left connection
    !for sigma=1 (i.e. rightwards moving particles).
    call init_pass_ends(pass_right,'r',1,'c')
  end subroutine init_pass_right

  subroutine init_pass_left
    !Create a pass_left object to pass boundary values to the left connection
    !for sigma=2 (i.e. leftwards moving particles).
    call init_pass_ends(pass_left,'l',2,'c')
  end subroutine init_pass_left

  subroutine init_pass_ends(pass_obj,dir,sigma,typestr,ngsend)
    !<DD> 01/02/2013: Copy of CMR's init_pass_right (revision 2085) 
    !but in customisable direction to aid code reuse and with optional
    !range in theta (i.e. so can pass more than just boundary points)
    
    !An example use is:
    ! call init_pass_ends(pass_right,'r',1,'c',3)
    ! call fill(pass_right,gnew,gnew)
    !This will pass gnew(ntgrid-2:ntgrid,1,iglo) to 
    !gnew(-ntgrid-2:-ntgrid,1,iglo_linked)
    
    !TODO:
    !      1) May be helpful to be able to pass both sigmas at once (e.g. for explicit scheme)
    use gs2_layouts, only: g_lo, il_idx, idx, proc_id
    use le_grids, only: ng2
    use mp, only: iproc, nproc, max_allreduce, proc0
    use redistribute, only: index_list_type, init_fill, delete_list, redist_type
    use theta_grid, only:ntgrid

    implicit none

    type (redist_type), intent(out) :: pass_obj !Redist type object to hold communication logic
    character(1),intent(in)::dir                  !Character string for direction of communication, should be 'l' for left and 'r' for right
    character(1),intent(in)::typestr              !Character string for type of data to be communicated. Should be 'c','r','i' or 'l'
    integer,intent(in) :: sigma                   !Which sigma index to send
    integer,intent(in),optional::ngsend           !How many theta grid points are we sending

    !Internal variables
    type (index_list_type), dimension(0:nproc-1) :: to, from
    integer, dimension (0:nproc-1) :: nn_from, nn_to
    integer, dimension(3) :: from_low, from_high, to_low, to_high
    integer :: il, iglo, ip, iglo_con, ipcon, n, nn_max, j
    logical,parameter :: debug=.false.
    integer :: bound_sign
    integer :: local_ngsend

    !Only applies to linked boundary option so exit if not linked
    if (boundary_option_switch .ne. boundary_option_linked) return

    !Handle the direction sign, basically we're either doing
    !     ntgrid --> -ntgrid (passing to right)
    !or 
    !    -ntgrid --> ntgrid  (passing to left)
    if (dir=='r') then
       bound_sign=1
    elseif (dir=='l') then
       bound_sign=-1
    else
       if (proc0) write(6,*) "Invalid direction string passed to init_pass_ends, defaulting to 'r'"
       bound_sign=1
    endif

    !Set the default number of theta grid points to send, if required
    if (present(ngsend)) then
       local_ngsend=ngsend
    else
       local_ngsend=1
    endif

    if (proc0.and.debug) write (6,*) "Initialising redist_type with settings Direction : ",dir," sigma",sigma,"local_ngsend",local_ngsend
    
    ! Need communications to satisfy || boundary conditions
    ! First find required blocksizes 
    
    !Initialise variables used to count how many entries to send and receive
    !for each processor
    nn_to = 0   ! # nn_to(ip) = communicates from ip TO HERE (iproc)
    nn_from = 0 ! # nn_from(ip) = communicates to ip FROM HERE (iproc)
    
    !Now loop over >all< iglo indices and work out how much data needs to be sent and received by each processor
    !Hence this routine does not scale with number of procs, see updated redist object creation in ccfe_opt_test (~r2173)
    !for examples which do scale
    do iglo = g_lo%llim_world, g_lo%ulim_world
       !Get the lambda index so we can skip trapped particles
       il = il_idx(g_lo,iglo)
       
       !Exclude disconnected trapped particles
       if (il > ng2+1) cycle     
       
       !Get the processor id for the proc holding the current iglo index
       ip = proc_id(g_lo,iglo)
       
       !What iglo connects to the current one in the direction of interest (and what proc is it on)
       !Note ipcon is <0 if no connection in direction of interest
       if (bound_sign==1) then
          call get_right_connection (iglo, iglo_con, ipcon)
       else
          call get_left_connection (iglo, iglo_con, ipcon)
       endif
       
       !Is the connected tube's proc the current processor?
       if (ipcon .eq. iproc ) then
          !If so add an entry recording an extra piece of information is to be sent
          !to proc ip (the proc holding iglo) from this proc (ipcon)
          !Note: Here we assume theta grid points are all on same proc 
          nn_to(ip)=nn_to(ip)+local_ngsend
       endif
       
       !Is the proc holding iglo this proc and is there a connection in the direction of interest?
       if (ip .eq. iproc .and. ipcon .ge. 0 ) then
          !If so add an entry recording an extra piece of information is to be sent
          !from this proc (ipcon) to ip
          !Note: Here we assume theta grid points are all on same proc 
          nn_from(ipcon)=nn_from(ipcon)+local_ngsend
       endif
    end do
    
    !Find the maximum amount of data to be received by a given processor
    !(first do it locally)
    nn_max = maxval(nn_to)
    !(now do it globally)
    call max_allreduce (nn_max)
    
    !Bit of debug printing
    if (proc0.and.debug) then
       write(6,*) 'init_pass_ends (1) processor, nn_to, nn_from:',iproc,nn_to,nn_from
    endif

    !Now that we've worked out how much data needs to be sent and received, define what specific
    !data needs to be sent to where
    if (nn_max.gt.0) then
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
       !<DD>, 06/01/2013: This redist object consists of a buffer of length n to hold the 
       !data during transfer and (currently) 3*2 integer arrays each of length n to hold
       !the indices of sent and received data.
       !By using c_fill_1 2*2 of these integer arrays would be removed. Assuming a double complex
       !buffer and long integer indices a 4n long array saving would be equivalent to the buffer
       !size and as such should represent a good memory saving but would not effect the amount
       !of data communicated (obviously).

       !Create to and from list objects for each processor and 
       !create storage to hold information about each specific from/to
       !communication
       do ip = 0, nproc-1
          !If proc ip is sending anything to this processor (iproc)
          if (nn_from(ip) > 0) then
             allocate (from(ip)%first(nn_from(ip)))
             allocate (from(ip)%second(nn_from(ip)))
             allocate (from(ip)%third(nn_from(ip)))
          endif
          !If proc ip is receiving anything from this processor (iproc)
          if (nn_to(ip) > 0) then
             allocate (to(ip)%first(nn_to(ip)))
             allocate (to(ip)%second(nn_to(ip)))
             allocate (to(ip)%third(nn_to(ip)))
          endif
       end do

       !Now fill the indirect addresses...
       
       !Initialise counters used to record how many pieces of data to expect
       nn_from = 0 ; nn_to = 0
       
       !Loop over >all< iglo indices
       do iglo = g_lo%llim_world, g_lo%ulim_world
          !Get the lambda index so we can skip trapped particles
          il = il_idx(g_lo,iglo)
          
          !Exclude disconnected trapped particles
          if (il > ng2+1) cycle     
          
          !What's the processor for the current iglo
          ip = proc_id(g_lo,iglo)
          
          !What iglo connects to the current one in the direction of interest (and what proc is it on)?
          !Note ipcon is <0 if no connection in direction of interest
          if (bound_sign==1) then
             call get_right_connection (iglo, iglo_con, ipcon)
          else
             call get_left_connection (iglo, iglo_con, ipcon)
          endif
          
          !For current proc for current iglo if there's connections in direction
          !then add an entry to the connected procs list of data to expect
          if (ip .eq. iproc .and. ipcon .ge. 0 ) then
             !Loop over theta grid indices
             !Note: Here we assume theta grid points are all on same proc 
             !Note: By looping over theta inside iglo loop we should optimise
             !      memory access compared to looping over theta outside.
             do j=0,local_ngsend-1
                n=nn_from(ipcon)+1 ; nn_from(ipcon)=n
                from(ipcon)%first(n)=bound_sign*(ntgrid-j) !Which theta point to send
                from(ipcon)%second(n)=sigma     !Sigma grid index to send
                from(ipcon)%third(n)=iglo   !iglo index to send
             enddo
          endif
          
          !If target iglo (iglo_con) is on this processor then add an entry recording where
          !we need to put the data when we receive it.
          if (ipcon .eq. iproc ) then 
             !Loop over theta grid indices
             !Note: Here we assume theta grid points are all on same proc 
             do j=0,local_ngsend-1
                n=nn_to(ip)+1 ; nn_to(ip)=n
                to(ip)%first(n)=-bound_sign*(ntgrid-j) !Which theta to store received data
                to(ip)%second(n)=sigma       !Sigma grid index to store received data
                to(ip)%third(n)=iglo_con !iglo index to store received data
             enddo
          endif
       end do
       
       !Bit of debug printing,
       if (debug.and.proc0) then
          write(6,*) 'init_pass_ends (2) processor, nn_to, nn_from:',iproc,nn_to,nn_from
       endif

       !Set data ranges for arrays to be passed, not this just effects how
       !arrays are indexed, not how big the buffer is.
       from_low(1)=-ntgrid ; from_low(2)=1  ; from_low(3)=g_lo%llim_proc       
       from_high(1)=ntgrid ; from_high(2)=2 ; from_high(3)=g_lo%ulim_alloc
       to_low(1)=-ntgrid   ; to_low(2)=1    ; to_low(3)=g_lo%llim_proc       
       to_high(1)=ntgrid   ; to_high(2)=2   ; to_high(3)=g_lo%ulim_alloc
       
       !Initialise fill object
       call init_fill (pass_obj, typestr, to_low, to_high, to, &
            from_low, from_high, from)
       
       !Clean up lists
       call delete_list (from)
       call delete_list (to)
    endif
  end subroutine init_pass_ends

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

  !>Helper function for finding the leftmost it of supercell
  function get_leftmost_it(it,ik)
    implicit none
    integer, intent(in) :: ik, it
    integer :: get_leftmost_it
    integer :: it_cur
    !If not linked then no connections so only one cell in supercell
    if (boundary_option_switch.eq.boundary_option_linked) then
       it_cur=it
       do while(itleft(ik,it_cur).ge.0.and.itleft(ik,it_cur).ne.it_cur)
          !Keep updating it_cur with left connected it until there are no
          !connections to left
          it_cur=itleft(ik,it_cur)
       enddo
    else
       it_cur=it
    endif
    get_leftmost_it=it_cur
  end function get_leftmost_it

  !>Helper function for finding the rightmost it of supercell
  function get_rightmost_it(it,ik)
    implicit none
    integer, intent(in) :: ik, it
    integer :: get_rightmost_it
    integer :: it_cur
    !If not linked then no connections so only one cell in supercell
    if (boundary_option_switch.eq.boundary_option_linked) then
       it_cur=it
       do while(itright(ik,it_cur).ge.0.and.itright(ik,it_cur).ne.it_cur)
          !Keep updating it_cur with right connected it until there are no
          !connections to right
          it_cur=itright(ik,it_cur)
       enddo
    else
       it_cur=it
    endif
    get_rightmost_it=it_cur
  end function get_rightmost_it

  subroutine allocate_arrays
    use kt_grids, only: naky,  box
    use theta_grid, only: ntgrid, shat
    use dist_fn_arrays, only: g, gnew
    use dist_fn_arrays, only: kx_shift, theta0_shift   ! MR
    use gs2_layouts, only: g_lo
    use nonlinear_terms, only: nonlin
    use run_parameters, only: fapar
    implicit none
!    logical :: alloc = .true.

!    if (alloc) then
    if (.not. allocated(g)) then
       allocate (g    (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (gnew (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (g0   (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       if(opt_source)then
          if(fapar.gt.0)then
             allocate (source_coeffs(4,-ntgrid:ntgrid-1,2,g_lo%llim_proc:g_lo%ulim_alloc))
          else
             allocate (source_coeffs(2,-ntgrid:ntgrid-1,2,g_lo%llim_proc:g_lo%ulim_alloc))
          endif
       endif
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

    g = 0. ; gnew = 0. ; g0 = 0.

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
    use le_derivatives, only: vspace_derivatives
    use dist_fn_arrays, only: gnew, g
    use nonlinear_terms, only: add_explicit_terms
    use hyper, only: hyper_diff
    use run_parameters, only: reset
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    integer, optional, intent (in) :: mode
    integer :: modep
    integer,save :: istep_last=-1

    modep = 0
    if (present(mode)) modep = mode

    if (istep.eq.0 .or. .not.cllc) then
!CMR do the usual LC when computing response matrix
       !Calculate the explicit nonlinear terms
       call add_explicit_terms (gexp_1, gexp_2, gexp_3, &
         phi, apar, bpar, istep, bkdiff(1))
       if(reset) return !Return if resetting
       !Solve for gnew
       call invert_rhs (phi, apar, bpar, phinew, aparnew, bparnew, istep)

       !Add hyper terms (damping)
       call hyper_diff (gnew, phinew)
       !Add collisions
       call vspace_derivatives (gnew, g, g0, phi, bpar, phinew, bparnew, modep)
    else if (istep.eq.1 .and. istep.ne.istep_last) then
!CMR on first half step at istep=1 do CL with all redists
       call vspace_derivatives (gnew, g, g0, phi, bpar, phinew, bparnew, modep)
       !Calculate the explicit nonlinear terms
       call add_explicit_terms (gexp_1, gexp_2, gexp_3, &
         phi, apar, bpar, istep, bkdiff(1))
       if(reset) return !Return if resetting
       !Solve for gnew
       call invert_rhs (phi, apar, bpar, phinew, aparnew, bparnew, istep)
       !Add hyper terms (damping)
       call hyper_diff (gnew, phinew)
    else if (istep.ne.istep_last) then
!CMR on first half step do CL with all redists without gtoc redist
       call vspace_derivatives (gnew, g, g0, phi, bpar, phinew, bparnew, modep,gtoc=.false.)
       !Calculate the explicit nonlinear terms
       call add_explicit_terms (gexp_1, gexp_2, gexp_3, &
         phi, apar, bpar, istep, bkdiff(1))
       if(reset) return !Return if resetting
       !Solve for gnew
       call invert_rhs (phi, apar, bpar, phinew, aparnew, bparnew, istep)
       !Add hyper terms (damping)
       call hyper_diff (gnew, phinew)
    else if (istep.eq.istep_last) then
!CMR on second half of all timesteps do LC without ctog redist
       !Calculate the explicit nonlinear terms 
       !NB following call should be unnecessary as NL terms not added on second
       !    half of istep: keeping for now as may be needed by some code
       call add_explicit_terms (gexp_1, gexp_2, gexp_3, &
         phi, apar, bpar, istep, bkdiff(1))
       if(reset) return !Return if resetting

       !Solve for gnew
       call invert_rhs (phi, apar, bpar, phinew, aparnew, bparnew, istep)
       !Add hyper terms (damping)
       call hyper_diff (gnew, phinew)
       call vspace_derivatives (gnew, g, g0, phi, bpar, phinew, bparnew, modep, ctog=.false.)
    endif

    !Enforce parity if desired (also performed in invert_rhs, but this is required
    !as collisions etc. may break parity?)
    if (def_parity) call enforce_parity(parity_redist)
      
    istep_last=istep

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

  subroutine exb_shear (g0, phi, apar, bpar, istep, field_local)
! MR, 2007: modified Bill Dorland's version to include grids where kx grid
!           is split over different processors
! MR, March 2009: ExB shear now available on extended theta grid (ballooning)
! CMR, May 2009: 2pishat correction factor on extended theta grid (ballooning)
!                so GEXB is same physical quantity in box and ballooning
! CMR, Oct 2010: multiply timestep by tunits(iky) for runs with wstar_units=.t.
! CMR, Oct 2010: add save statements to prevent potentially large and memory 
!                killing array allocations!
   
    use mp, only: send, receive, mp_abort, broadcast
    use gs2_layouts, only: ik_idx, it_idx, g_lo, idx_local, idx, proc_id
    use run_parameters, only: tunits
    use theta_grid, only: ntgrid, ntheta, shat
    use file_utils, only: error_unit
    use kt_grids, only: akx, aky, naky, ikx, ntheta0, box, theta0
    use le_grids, only: negrid, nlambda
    use species, only: nspec
    use run_parameters, only: fphi, fapar, fbpar
    use dist_fn_arrays, only: kx_shift, theta0_shift
    use gs2_time, only: code_dt, code_dt_old, code_time
    use constants, only: twopi    

    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g0
    integer, intent(in) :: istep
    complex, dimension(:,:,:), allocatable :: temp 
    complex, dimension(:,:), allocatable :: temp2
    integer, dimension(1), save :: itmin
    integer :: ierr, j 
    integer :: ik, it, ie, is, il, isgn, to_iglo, from_iglo
    integer:: iib, iit, ileft, iright, i
    integer, save :: istep_last=0
    !integer, intent(in) :: istep EGH commented... faulty merge?
    logical, intent(in), optional :: field_local
    logical :: field_local_loc
    real, save :: dkx, dtheta0
    real :: gdt
    !logical, save :: exb_first = .true.
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
             !Warning, dkx has not been set => Should really halt run or find
             !a suitable default definition for dkx.
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

    !Check if we want to exit without applying flow shear
    if (istep.eq.istep_last) return !Don't allow flow shear to happen more than once per step    
    if (g_exb_start_timestep > istep) return !Flow shear not active yet, set in timesteps
    if (g_exb_start_time >= 0 .and. code_time < g_exb_start_time) return !Flow shear not active yet, set in time
    
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

    !Update istep_last
    istep_last = istep

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

    !If using field_option='local' and x/it is not entirely local
    !then we need to make sure that all procs know the full field
    !THIS IS A TEMPORARY FIX AND WE SHOULD PROBABLY DO SOMETHING BETTER
    field_local_loc=.false.
    if(present(field_local)) field_local_loc=field_local
    if(any(jump.ne.0).and.(.not.g_lo%x_local).and.field_local_loc) then
       if(fphi>epsilon(0.0)) call broadcast(phi)
       if(fapar>epsilon(0.0)) call broadcast(apar)
       if(fbpar>epsilon(0.0)) call broadcast(bpar)
    endif

    if (.not. box .and. shat .ne. 0.0 ) then
! MR, March 2009: impact of ExB shear on extended theta grid computed here
!                 for finite shat
       do ik =1,naky
          j=jump(ik)
          if (j .eq. 0) cycle     
          if (abs(j) .ge. ntheta0) then
              write(str,fmt='("in exb_shear: jump(ik)=",i4," > ntheta0 =",i4," for ik=",i4,". => reduce timestep or increase ntheta0")') j,ik,ntheta0
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

  !<DD>Subroutine to setup a redistribute object to be used in enforcing parity
  subroutine init_enforce_parity(redist_obj,ik_ind)
    use theta_grid, only : ntgrid
    use gs2_layouts, only : g_lo,proc_id, ik_idx
    use redistribute, only: index_list_type, init_fill, delete_list, redist_type
    use mp, only: iproc, nproc, max_allreduce
    implicit none
    type(redist_type), intent(out) :: redist_obj
    integer, intent(in), optional :: ik_ind !If present then resulting redistribute object only applies to ik=ik_ind
    integer :: iglo,iglo_conn,ip_conn,ip
    integer :: ig, ndata, nn_max, n,ik_ind_local
    type (index_list_type), dimension(0:nproc-1) :: to, from
    integer, dimension (0:nproc-1) :: nn_from, nn_to
    integer, dimension(3) :: from_low, from_high, to_low, to_high

    !If enforced parity not requested then exit
    if(.not.(def_parity))return
    
    !If not linked then don't need any setup so exit
    if(boundary_option_switch.ne.boundary_option_linked)return
    
    !Deal with optional input
    if(present(ik_ind)) then
       ik_ind_local=ik_ind
    else
       ik_ind_local=-1
    endif

    !NOTE: If we know x is entirely local (g_lo%x_local=.true.) then we
    !know that we don't need to do any comms so our iglo range can be
    !limited to what is on this proc. At the moment we don't take 
    !advantage of this.

    !Use a shorthand for how much data to send per iglo
    ndata=2*ntgrid+1

    !Count how much data is to be sent/recv to/from different procs
    !/First initialise counters
    nn_to =0
    nn_from=0
    
    !/Now loop over iglo to workout what data is to be sent/received
    do iglo=g_lo%llim_world,g_lo%ulim_world
       !Check if we want to skip this
       if(ik_ind_local.gt.0)then
          if(ik_idx(g_lo,iglo).ne.ik_ind_local) cycle
       endif

       !Get proc id of current iglo
       ip=proc_id(g_lo,iglo)

       !Find connected iglo to which we want to send the data
       call get_parity_conn(iglo,iglo_conn,ip_conn)

       !Do we have the data to send?
       if(ip_conn.eq.iproc) nn_to(ip)=nn_to(ip)+ndata

       !Are we going to receive the data?
       if(ip.eq.iproc) nn_from(ip_conn)=nn_from(ip_conn)+ndata

    enddo

    !Now find the maxmimum amount of data to be sent
    nn_max=MAXVAL(nn_to) !Local max
    call max_allreduce(nn_max) !Global max

    !Now define indices of send/receive data
    if(nn_max.gt.0) then
       !Create to/from list objects
       do ip=0,nproc-1
          if(nn_from(ip)>0)then
             allocate(from(ip)%first(nn_from(ip)))
             allocate(from(ip)%second(nn_from(ip)))
             allocate(from(ip)%third(nn_from(ip)))
          endif
          if(nn_to(ip)>0)then
             allocate(to(ip)%first(nn_to(ip)))
             allocate(to(ip)%second(nn_to(ip)))
             allocate(to(ip)%third(nn_to(ip)))
          endif
       enddo

       !Now fill in the lists
       nn_from=0
       nn_to=0

       !Loop over all iglo again (this doesn't scale)
       do iglo=g_lo%llim_world, g_lo%ulim_world
          !Check if we want to skip this
          if(ik_ind_local.gt.0)then
             if(ik_idx(g_lo,iglo).ne.ik_ind_local) cycle
          endif
          
          !Get proc for this iglo
          ip=proc_id(g_lo,iglo)

          !Get connected point
          call get_parity_conn(iglo,iglo_conn,ip_conn)

          !If we're receiving data where do we put it?
          if(ip.eq.iproc)then
             do ig=-ntgrid,ntgrid !Optimised for data send
                n=nn_from(ip_conn)+1
                nn_from(ip_conn)=n
                from(ip_conn)%first(n)=0-ig
                from(ip_conn)%second(n)=1
                from(ip_conn)%third(n)=iglo
             enddo
          endif

          !If we're sending data where do we get it from?
          if(ip_conn.eq.iproc)then
             do ig=-ntgrid,ntgrid !Optimised for data send
                n=nn_to(ip)+1
                nn_to(ip)=n
                to(ip)%first(n)=ig
                to(ip)%second(n)=2
                to(ip)%third(n)=iglo_conn
             enddo
          endif
       enddo

       !Now setup the redistribute object
       from_low(1)=-ntgrid ; from_low(2)=1  ; from_low(3)=g_lo%llim_proc
       from_high(1)=ntgrid ; from_high(2)=2 ; from_high(3)=g_lo%ulim_proc
       to_low(1)=-ntgrid ; to_low(2)=1  ; to_low(3)=g_lo%llim_proc
       to_high(1)=ntgrid ; to_high(2)=2 ; to_high(3)=g_lo%ulim_proc

       !Initialise the fill object
       call init_fill(redist_obj,'c',to_low,to_high,to,&
            from_low,from_high, from)

       !Delete lists
       call delete_list(from)
       call delete_list(to)
    endif
  end subroutine init_enforce_parity
  !</DD>

  !<DD>Subroutine to return the iglo corresponding to the part of the domain given
  !by iglo reflected in theta=theta0
  subroutine get_parity_conn(iglo,iglo_conn,iproc_conn)
    use gs2_layouts, only: ik_idx,it_idx,proc_id,g_lo
    implicit none
    integer, intent(in) :: iglo
    integer, intent(out) :: iglo_conn, iproc_conn
    integer :: it, ik, it_conn, link, tmp

    !Get indices
    it=it_idx(g_lo,iglo)
    ik=ik_idx(g_lo,iglo)

    !Initialise to this cell
    tmp=iglo
    
    !Now check number of links
    if(l_links(ik,it).eq.r_links(ik,it))then
       !If in the middle of the domain then iglo doesn't change
       !Just get the proc id
       iproc_conn=proc_id(g_lo,iglo)
       iglo_conn=tmp
    elseif(l_links(ik,it).gt.r_links(ik,it))then
       !If we're on the right then look left
       do link=1,l_links(ik,it)
          !Get the next connected domain
          call get_left_connection(tmp,iglo_conn,iproc_conn)

          !What is the it index here?
          it_conn=it_idx(g_lo,iglo_conn)

          !Update current iglo
          tmp=iglo_conn

          !If the number of right links now matches the left then we've got the match
          if(r_links(ik,it_conn).eq.l_links(ik,it)) exit
       enddo
    else
       !If we're on the left then look right
       do link=1,r_links(ik,it)
          !Get the next connected domain
          call get_right_connection(tmp,iglo_conn,iproc_conn)

          !What is the it index here?
          it_conn=it_idx(g_lo,iglo_conn)

          !Update current iglo
          tmp=iglo_conn

          !If the number of left links now matches the right then we've got the match
          if(l_links(ik,it_conn).eq.r_links(ik,it)) exit
       enddo
    endif
  end subroutine get_parity_conn
  !</DD>

  !<DD>Subroutine to enforce requested parity
  subroutine enforce_parity(redist_obj,ik_ind)
    use theta_grid, only:ntgrid
    use dist_fn_arrays, only: gnew
    use redistribute, only: scatter,redist_type
    use gs2_layouts, only: g_lo,ik_idx
    implicit none
    type(redist_type),intent(in),optional :: redist_obj
    integer, intent(in),optional :: ik_ind
    type(redist_type) :: redist_local
    integer :: ik_local
    integer :: iglo,mult

    !If enforced parity not requested then exit
    if(.not.(def_parity))return
    
    !Set multiplier
    if(even) then
       mult=1
    else
       mult=-1
    endif
    
    !Behaviour depends upon if we're in flux tube or ballooning space
    if(boundary_option_switch.eq.boundary_option_linked) then !Flux-tube
       !Ensure a redist object is present, if not default to parity_redist
       if(present(redist_obj))then
          redist_local=redist_obj
       else
          redist_local=parity_redist
       endif

       !Redistribute the data
       call scatter(redist_local,gnew,gnew)

       !Multiply by factor if required
       if(mult.ne.1) gnew(:,1,:)=mult*gnew(:,1,:)
    else !Ballooning/extended space
       !Ensure ik_local is specified
       if(present(ik_ind))then
          ik_local=ik_ind
       else
          ik_local=-1
       endif

       !Loop over all local iglo
       do iglo=g_lo%llim_proc,g_lo%ulim_alloc
          !Skip if needed
          if(ik_local.gt.0) then
             if(ik_idx(g_lo,iglo).ne.ik_local) cycle
          endif

          !Apply parity filter
          gnew(-ntgrid:-1,1,iglo)=mult*gnew(ntgrid:1:-1,2,iglo)
          gnew(1:ntgrid,1,iglo)=mult*gnew(-1:-ntgrid:-1,2,iglo)
       enddo
    endif
    
  end subroutine enforce_parity
  !</DD>

  subroutine get_source_term &
       (phi, apar, bpar, phinew, aparnew, bparnew, istep, &
        isgn, iglo,ik,it,il,ie,is, sourcefac, source)
#ifdef LOWFLOW
    use dist_fn_arrays, only: hneoc, vparterm, wdfac, wstarfac, wdttpfac
#endif
    use dist_fn_arrays, only: aj0, aj1, vperp2, vpar, vpac, g, ittp
    use theta_grid, only: ntgrid
    use kt_grids, only: aky
    use le_grids, only: nlambda, ng2, lmax, anon, negrid
    use species, only: spec, nspec
    use run_parameters, only: fphi, fapar, fbpar, wunits
    use gs2_time, only: code_dt
    use gs2_time, only: dt0 => code_dt, dt1 => code_dt_prev1, dt2 => code_dt_prev2
    use nonlinear_terms, only: nonlin
    use hyper, only: D_res
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    integer, intent (in) :: isgn, iglo, ik, it, il, ie, is
    complex, intent (in) :: sourcefac
    complex, dimension (-ntgrid:), intent (out) :: source

    integer :: ig
    complex, dimension (-ntgrid:ntgrid) :: phigavg, apargavg
    real :: c0, c1, c2


    ! NRM, 3/18/2015 
    ! calculate AB3 coefficients generalized to changing timestep
    ! see Tatsuno & Dorland, Physics of Plasmas 13, 092107 (2006)
    ! http://www.cscamm.umd.edu/publications/TatsunoDorland-PoP06_CS-06-11.pdf
    ! if dt0 = dt1 = dt2, this reduces to usual AB3 coefficients
    c0 = (1. / (dt1 + dt2)) * ( &
         ((dt0 + dt1)**2./dt1)*( (dt0+dt1)/3. + dt2/2. ) &
         - dt1*(dt1/3. + dt2/2.) )
    c1 = - dt0**2. * ( dt0/3. + (dt1+dt2)/2. ) / (dt1*dt2)
    c2 = dt0**2. * ( dt0/3. + dt1/2. ) / ( dt2*(dt1+dt2) )

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
                  - zi*anon(ie)*wdrift(:ntgrid-1,isgn,iglo) &
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

!CMR, 21/7/2014: removed redundant line here:    source(ntgrid)=source(-ntgrid) 
!                as source(ntgrid) should never be used.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! special source term for totally trapped particles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CMR, 13/10/2014: 
! Upper limit of following loops setting source changed from ntgrid to ntgrid-1
! Source is allocated as: source(-ntgrid:ntgrid-1), so ntgrid is out of bounds.

    if (source_option_switch == source_option_full .or. &
        source_option_switch == source_option_phiext_full) then
       if (nlambda > ng2 .and. isgn == 2) then
          do ig = -ntgrid, ntgrid-1
             if (il < ittp(ig)) cycle
             source(ig) &
                  = g(ig,2,iglo)*a(ig,2,iglo) &
#ifdef LOWFLOW
                  - anon(ie)*zi*(wdttpfac(ig,it,ik,ie,is,2)*hneoc(ig,2,iglo))*phigavg(ig) &
                  + zi*wstar(ik,ie,is)*hneoc(ig,2,iglo)*phigavg(ig)
#else
                  - anon(ie)*zi*(wdriftttp(ig,it,ik,ie,is,2))*phigavg(ig) &
                  + zi*wstar(ik,ie,is)*phigavg(ig)
#endif             
          end do

          if (source_option_switch == source_option_phiext_full .and. &
               aky(ik) < epsilon(0.0)) then
             do ig = -ntgrid, ntgrid-1
                if (il < ittp(ig)) cycle             
                source(ig) = source(ig) - zi*anon(ie)* &
#ifdef LOWFLOW
                     wdttpfac(ig,it,ik,ie,is,isgn)*2.0*phi_ext*sourcefac*aj0(ig,iglo)
#else
                     wdriftttp(ig,it,ik,ie,is,isgn)*2.0*phi_ext*sourcefac*aj0(ig,iglo)
#endif
             end do
          endif

#ifdef LOWFLOW
          select case (istep)
          case (0)
             ! nothing
          case (1)
             do ig = -ntgrid, ntgrid-1
                if (il < ittp(ig)) cycle
                source(ig) = source(ig) + 0.5*code_dt*gexp_1(ig,isgn,iglo)
             end do
          case (2) 
             do ig = -ntgrid, ntgrid-1
                if (il < ittp(ig)) cycle
                source(ig) = source(ig) + 0.5*code_dt*( &
                     1.5*gexp_1(ig,isgn,iglo) - 0.5*gexp_2(ig,isgn,iglo))
             end do
          case default
             do ig = -ntgrid, ntgrid-1
                if (il < ittp(ig)) cycle
             !   source(ig) = source(ig) + 0.5*code_dt*( &
             !        (23./12.)*gexp_1(ig,isgn,iglo) &
             !        - (4./3.)  *gexp_2(ig,isgn,iglo) &
             !        + (5./12.) *gexp_3(ig,isgn,iglo))
    ! NRM, 3/18/2015 
    ! use AB3 generalized to changing timestep
                 source(ig) = source(ig) + 0.5*( &
                          c0*gexp_1(ig,isgn,iglo) &
                        + c1*gexp_2(ig,isgn,iglo) &
                        + c2*gexp_3(ig,isgn,iglo))
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
                   if (il < ittp(ig)) cycle
                   source(ig) = source(ig) + 0.5*code_dt*gexp_1(ig,isgn,iglo)
                end do
             case (2) 
                do ig = -ntgrid, ntgrid-1
                   if (il < ittp(ig)) cycle
                   source(ig) = source(ig) + 0.5*code_dt*( &
                        1.5*gexp_1(ig,isgn,iglo) - 0.5*gexp_2(ig,isgn,iglo))
                end do                   
             case default
                do ig = -ntgrid, ntgrid-1
                   if (il < ittp(ig)) cycle
              !    source(ig) = source(ig) + 0.5*code_dt*( &
              !           (23./12.)*gexp_1(ig,isgn,iglo) &
              !         - (4./3.)  *gexp_2(ig,isgn,iglo) &
              !         + (5./12.) *gexp_3(ig,isgn,iglo))
    ! NRM, 3/18/2015 
    ! use AB3 generalized to changing timestep
                 source(ig) = source(ig) + 0.5*( &
                          c0*gexp_1(ig,isgn,iglo) &
                        + c1*gexp_2(ig,isgn,iglo) &
                        + c2*gexp_3(ig,isgn,iglo))
                end do
             end select
          end if
#endif

       end if
    end if

  contains

    subroutine set_source

      use species, only: spec
      use theta_grid, only: itor_over_B
      implicit none
      complex :: apar_p, apar_m, phi_p, phi_m!, bpar_p !GGH added bpar_p
      real :: bd, bdfac_p, bdfac_m

! try fixing bkdiff dependence
      bd = bkdiff(1)

      bdfac_p = 1.+bd*(3.-2.*real(isgn))
      bdfac_m = 1.-bd*(3.-2.*real(isgn))



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

      do ig = -ntgrid, ntgrid-1
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
         source(ig) = anon(ie)*(-2.0*vpar(ig,isgn,iglo)*phi_m &
              -spec(is)%zstm*vpac(ig,isgn,iglo) &
              *((aj0(ig+1,iglo) + aj0(ig,iglo))*0.5*apar_m  &
              + D_res(it,ik)*apar_p) &
              -zi*wdrift(ig,isgn,iglo)*phi_p) &
              + zi*(wstar(ik,ie,is) &
              + vpac(ig,isgn,iglo)*code_dt*wunits(ik)*ufac(ie,is) &
              -2.0*omprimfac*vpac(ig,isgn,iglo)*code_dt*wunits(ik)*g_exb*itor_over_B(ig)/spec(is)%stm) &
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
          !  source(ig) = source(ig) + 0.5*code_dt*( &
          !       (23./12.)*gexp_1(ig,isgn,iglo) &
          !       - (4./3.)  *gexp_2(ig,isgn,iglo) &
          !       + (5./12.) *gexp_3(ig,isgn,iglo))
    ! NRM, 3/18/2015 
    ! use AB3 generalized to changing timestep
                 source(ig) = source(ig) + 0.5*( &
                          c0*gexp_1(ig,isgn,iglo) &
                        + c1*gexp_2(ig,isgn,iglo) &
                        + c2*gexp_3(ig,isgn,iglo))
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
             !  source(ig) = source(ig) + 0.5*code_dt*( &
             !         (23./12.)*gexp_1(ig,isgn,iglo) &
             !       - (4./3.)  *gexp_2(ig,isgn,iglo) &
             !       + (5./12.) *gexp_3(ig,isgn,iglo))
    ! NRM, 3/18/2015 
    ! use AB3 generalized to changing timestep
                 source(ig) = source(ig) + 0.5*( &
                          c0*gexp_1(ig,isgn,iglo) &
                        + c1*gexp_2(ig,isgn,iglo) &
                        + c2*gexp_3(ig,isgn,iglo))
            end do
         end select
      end if
#endif

    end subroutine set_source

  end subroutine get_source_term

  !This is a version of get_source_term which does both sign (sigma) together
  !and uses precalculated constant terms. Leads to more memory usage than 
  !original version but can be significantly faster (~50%)
  subroutine get_source_term_opt &
       (phi, apar, bpar, phinew, aparnew, bparnew, istep, &
        iglo,ik,it,il,ie,is, sourcefac, source)
#ifdef LOWFLOW
    use dist_fn_arrays, only: hneoc, vparterm, wdfac, wstarfac, wdttpfac
#endif
    use dist_fn_arrays, only: aj0, aj1, vperp2, g, ittp
    use theta_grid, only: ntgrid
    use kt_grids, only: aky
    use le_grids, only: nlambda, ng2, lmax, anon
    use species, only: spec
    use run_parameters, only: fphi, fapar, fbpar
    use gs2_time, only: code_dt
    use gs2_time, only: dt0 => code_dt, dt1 => code_dt_prev1, dt2 => code_dt_prev2
    use nonlinear_terms, only: nonlin
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    integer, intent (in) :: iglo, ik, it, il, ie, is
    complex, intent (in) :: sourcefac
    complex, dimension (-ntgrid:,:), intent (out) :: source
    integer :: ig, isgn
    complex, dimension (-ntgrid:ntgrid) :: phigavg, apargavg

    real :: c0, c1, c2

    ! NRM, 3/18/2015 
    ! calculate AB3 coefficients generalized to changing timestep
    ! see Tatsuno & Dorland, Physics of Plasmas 13, 092107 (2006)
    ! http://www.cscamm.umd.edu/publications/TatsunoDorland-PoP06_CS-06-11.pdf
    ! if dt0 = dt1 = dt2, this reduces to usual AB3 coefficients
    c0 = (1. / (dt1 + dt2)) * ( &
         ((dt0 + dt1)**2./dt1)*( (dt0+dt1)/3. + dt2/2. ) &
         - dt1*(dt1/3. + dt2/2.) )
    c1 = - dt0**2. * ( dt0/3. + (dt1+dt2)/2. ) / (dt1*dt2)
    c2 = dt0**2. * ( dt0/3. + dt1/2. ) / ( dt2*(dt1+dt2) )


    !Temporally weighted fields
    phigavg  = (fexp(is)*phi(:,it,ik)   + (1.0-fexp(is))*phinew(:,it,ik)) &
                *aj0(:,iglo)*fphi

    if(fbpar.gt.0)then
       phigavg=phigavg+&
            (fexp(is)*bpar(:,it,ik) + (1.0-fexp(is))*bparnew(:,it,ik))&
            *aj1(:,iglo)*fbpar*2.0*vperp2(:,iglo)*spec(is)%tz
    endif

    if(fapar.gt.0)then
       apargavg = (fexp(is)*apar(:,it,ik)  + (1.0-fexp(is))*aparnew(:,it,ik)) &
            *aj0(:,iglo)*fapar
    endif

    !Initialise to zero in case where we don't call set_source_opt
    if(il>lmax)then
       source=0.0
    endif

    !Calculate the part of the source related to EM potentials
    !Do both signs at once to improve memory access
    do isgn=1,2
       ! source term in finite difference equations
       select case (source_option_switch)
       case (source_option_full)
          if (il <= lmax) then
             call set_source_opt
          end if
       case(source_option_phiext_full)
          if (il <= lmax) then
             call set_source_opt
             if (aky(ik) < epsilon(0.0)) then
                source(:ntgrid-1,isgn) = source(:ntgrid-1,isgn) &
                  - zi*anon(ie)*wdrift(:ntgrid-1,isgn,iglo) &
                  *2.0*phi_ext*sourcefac*aj0(:ntgrid-1,iglo)
             end if
          end if
       end select
    enddo
       
    ! Do matrix multiplications
    if (il <= lmax) then
       do ig = -ntgrid, ntgrid-1
          source(ig,1) = source(ig,1) &
               + b(ig,1,iglo)*g(ig,1,iglo) + a(ig,1,iglo)*g(ig+1,1,iglo)
       end do
!CMR, 21/7/2014: removed redundant line here:    source(ntgrid,1)=source(-ntgrid,1) 
!                as source(ntgrid,1) should never be used.
       do ig = -ntgrid, ntgrid-1
          source(ig,2) = source(ig,2) &
               + a(ig,2,iglo)*g(ig,2,iglo) + b(ig,2,iglo)*g(ig+1,2,iglo)
       end do
!CMR, 21/7/2014: removed redundant line here:    source(ntgrid,2)=source(-ntgrid,2) 
!                as source(ntgrid,2) should never be used.
    end if

    ! special source term for totally trapped particles (isgn=2 only)

!CMR, 13/10/2014: 
! Upper limit of following loops setting source changed from ntgrid to ntgrid-1
! Source allocated as: source(-ntgrid:ntgrid-1,2), so ntgrid is out of bounds.

    isgn=2
    if (source_option_switch == source_option_full .or. &
         source_option_switch == source_option_phiext_full) then
       if (nlambda > ng2) then
          do ig = -ntgrid, ntgrid-1
             if (il < ittp(ig)) cycle
             source(ig,isgn) &
                  = g(ig,2,iglo)*a(ig,2,iglo) &
#ifdef LOWFLOW
                  - anon(ie)*zi*(wdttpfac(ig,it,ik,ie,is,2)*hneoc(ig,2,iglo))*phigavg(ig) &
                  + zi*wstar(ik,ie,is)*hneoc(ig,2,iglo)*phigavg(ig)
#else
                  - anon(ie)*zi*(wdriftttp(ig,it,ik,ie,is,2))*phigavg(ig) &
                  + zi*wstar(ik,ie,is)*phigavg(ig)
#endif             
          end do
             
          if (source_option_switch == source_option_phiext_full .and. &
               aky(ik) < epsilon(0.0)) then
             do ig = -ntgrid, ntgrid-1
                if (il < ittp(ig)) cycle             
                source(ig,isgn) = source(ig,isgn) - zi*anon(ie)* &
#ifdef LOWFLOW
                     wdttpfac(ig,it,ik,ie,is,isgn)*2.0*phi_ext*sourcefac*aj0(ig,iglo)
#else
                     wdriftttp(ig,it,ik,ie,is,isgn)*2.0*phi_ext*sourcefac*aj0(ig,iglo)
#endif
             end do
          endif

          ! add in nonlinear terms 
          if (nonlin) then         
             select case (istep)
             case (0)
                ! nothing
             case (1)
                do ig = -ntgrid, ntgrid-1
                   if (il < ittp(ig)) cycle
                   source(ig,isgn) = source(ig,isgn) + 0.5*code_dt*gexp_1(ig,isgn,iglo)
                end do
             case (2) 
                do ig = -ntgrid, ntgrid-1
                   if (il < ittp(ig)) cycle
                   source(ig,isgn) = source(ig,isgn) + 0.5*code_dt*( &
                        1.5*gexp_1(ig,isgn,iglo) - 0.5*gexp_2(ig,isgn,iglo))
                end do
             case default
                do ig = -ntgrid, ntgrid-1
                   if (il < ittp(ig)) cycle
                !   source(ig,isgn) = source(ig,isgn) + 0.5*code_dt*( &
                !        (23./12.)*gexp_1(ig,isgn,iglo) &
                !        - (4./3.)  *gexp_2(ig,isgn,iglo) &
                !        + (5./12.) *gexp_3(ig,isgn,iglo))
    ! NRM, 3/18/2015 
    ! use AB3 generalized to changing timestep
                 source(ig,isgn) = source(ig,isgn) + 0.5*( &
                          c0*gexp_1(ig,isgn,iglo) &
                        + c1*gexp_2(ig,isgn,iglo) &
                        + c2*gexp_3(ig,isgn,iglo))
                end do
             end select
          end if
       end if
    end if

  contains

    subroutine set_source_opt
      implicit none
      complex :: apar_p, apar_m, phi_p, phi_m
      real :: bd, bdfac_p, bdfac_m

! try fixing bkdiff dependence
      bd = bkdiff(1)

      bdfac_p = 1.+bd*(3.-2.*real(isgn))
      bdfac_m = 1.-bd*(3.-2.*real(isgn))

      if(fapar.gt.0)then
         do ig = -ntgrid, ntgrid-1
            phi_p = bdfac_p*phigavg(ig+1)+bdfac_m*phigavg(ig)
            phi_m = phigavg(ig+1)-phigavg(ig)
            ! RN> bdfac factors seem missing for apar_p
            apar_p = apargavg(ig+1)+apargavg(ig)
            apar_m = aparnew(ig+1,it,ik)+aparnew(ig,it,ik) & 
                 -apar(ig+1,it,ik)-apar(ig,it,ik)
            
            source(ig,isgn)=phi_m*source_coeffs(1,ig,isgn,iglo)+&
                 phi_p*source_coeffs(2,ig,isgn,iglo)+&
                 apar_m*source_coeffs(3,ig,isgn,iglo)+&
                 apar_p*source_coeffs(4,ig,isgn,iglo)
         end do
      else
         do ig = -ntgrid, ntgrid-1
            phi_p = bdfac_p*phigavg(ig+1)+bdfac_m*phigavg(ig)
            phi_m = phigavg(ig+1)-phigavg(ig)
            
            source(ig,isgn)=phi_m*source_coeffs(1,ig,isgn,iglo)+&
                 phi_p*source_coeffs(2,ig,isgn,iglo)
         end do
      endif

! add in nonlinear terms 
      if (nonlin) then         
         select case (istep)
         case (0)
            ! nothing
         case (1)
            do ig = -ntgrid, ntgrid-1
               source(ig,isgn) = source(ig,isgn) + 0.5*code_dt*gexp_1(ig,isgn,iglo)
            end do
         case (2) 
            do ig = -ntgrid, ntgrid-1
               source(ig,isgn) = source(ig,isgn) + 0.5*code_dt*( &
                    1.5*gexp_1(ig,isgn,iglo) - 0.5*gexp_2(ig,isgn,iglo))
            end do
         case default
            do ig = -ntgrid, ntgrid-1
            !   source(ig,isgn) = source(ig,isgn) + 0.5*code_dt*( &
            !          (23./12.)*gexp_1(ig,isgn,iglo) &
            !        - (4./3.)  *gexp_2(ig,isgn,iglo) &
            !        + (5./12.) *gexp_3(ig,isgn,iglo))
    ! NRM, 3/18/2015 
    ! use AB3 generalized to changing timestep
                 source(ig,isgn) = source(ig,isgn) + 0.5*( &
                          c0*gexp_1(ig,isgn,iglo) &
                        + c1*gexp_2(ig,isgn,iglo) &
                        + c2*gexp_3(ig,isgn,iglo))
            end do
         end select
      end if

    end subroutine set_source_opt

  end subroutine get_source_term_opt

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
    use run_parameters, only: fbpar, fphi, ieqzip
    use kt_grids, only: kwork_filter
    use species, only: spec
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    integer, intent (in) :: iglo
    complex, intent (in) :: sourcefac

    integer :: ig, ik, it, il, ie, isgn, is
    integer :: ilmin
    complex :: beta1
    complex, dimension (-ntgrid:ntgrid,2) :: g1, g2
    complex, dimension (-ntgrid:ntgrid-1,2) :: source
    complex :: adjleft, adjright
    logical :: use_pass_homog, speriod_flag
    integer :: ntgl, ntgr

!    call prof_entering ("invert_rhs_1", "dist_fn")

    ik = ik_idx(g_lo,iglo)
    it = it_idx(g_lo,iglo)

    !Skip work if we're not interested in this ik and it
    if(kwork_filter(it,ik)) return

    if(ieqzip(it,ik)==0) return
    if (eqzip) then
       if (secondary .and. ik == 2 .and. it == 1) return ! do not evolve primary mode
       if (tertiary .and. ik == 1) then
          if (it == 2 .or. it == ntheta0) return ! do not evolve periodic equilibrium
       end if
       if (harris .and. ik == 1) return ! do not evolve primary mode       
    end if

    il = il_idx(g_lo,iglo)
    ie = ie_idx(g_lo,iglo)
    is = is_idx(g_lo,iglo)

    if(opt_source)then
       call get_source_term_opt (phi, apar, bpar, phinew, aparnew, bparnew, &
            istep, iglo,ik,it,il,ie,is, sourcefac, source)
    else
       do isgn = 1, 2
          call get_source_term (phi, apar, bpar, phinew, aparnew, bparnew, &
               istep, isgn, iglo,ik,it,il,ie,is, sourcefac, source(:,isgn))
       end do
    endif

!CMR, 17/4/2012:
!  use_pass_homog = T  iff one of following applies:   
!                 boundary_option_self_periodic
!     OR          boundary_option_linked
!     OR          aky=0
!  if use_pass_homog = T, compute homogeneous solution (g1) for passing.
!        otherwise ONLY compute inhomogenous solution for passing particles.
!
!  speriod_flag = T  iff boundary_option_linked AND aky=0
!  if speriod_flag = T, apply self-periodic bcs for passing and wfb.

    select case (boundary_option_switch)
    case (boundary_option_self_periodic)
       use_pass_homog = .true.
    case (boundary_option_linked)
       use_pass_homog = .true.
       speriod_flag = aky(ik) == 0.0
    case default
       use_pass_homog = .false.
    end select

    use_pass_homog = use_pass_homog .or. aky(ik) == 0.0

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
       if (il <= ng2) then
          !Only apply the new boundary condition to the leftmost
          !cell for sign going from left to right
          if (l_links(ik,it) .eq. 0) then
             adjleft = anon(ie)*2.0*vperp2(-ntgrid,iglo)*aj1(-ntgrid,iglo) &
                  *bparnew(-ntgrid,it,ik)*fbpar &
                  + spec(is)%z*anon(ie)*phinew(-ntgrid,it,ik)*aj0(-ntgrid,iglo) &
                  /spec(is)%temp*fphi
             gnew(-ntgrid,1,iglo) = gnew(-ntgrid,1,iglo) - adjleft
          end if
          !Only apply the new boundary condition to the rightmost
          !cell for sign going from right to left
          if (r_links(ik,it) .eq. 0) then
             adjright = anon(ie)*2.0*vperp2(ntgrid,iglo)*aj1(ntgrid,iglo) &
                  *bparnew(ntgrid,it,ik)*fbpar &
                  + spec(is)%z*anon(ie)*phinew(ntgrid,it,ik)*aj0(ntgrid,iglo) &
                  /spec(is)%temp*fphi
             gnew(ntgrid,2,iglo) = gnew(ntgrid,2,iglo) - adjright
         end if
       endif
    endif

    ! g1 is the homogeneous solution
    g1 = 0.0

    ntgl = -ntgrid
    ntgr =  ntgrid

! ng2+1 is WFB

    if (use_pass_homog) then
       if (il < ng2+1) then
          g1(-ntgrid,1) = 1.0
          g1( ntgrid,2) = 1.0
       end if
    end if

    if (il == ng2+1) then ! ng2+1 is WFB
       g1(-ntgrid,1) = wfb  ! wfb should be unity here; variable is for testing
       g1( ntgrid,2) = wfb  ! wfb should be unity here; variable is for testing
    endif

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

    if (use_pass_homog) then
       ilmin = 1
    else
       ilmin = ng2 + 1              !!! ilmin = ng2 + 2
    end if

! time advance vpar < 0 homogeneous part: g1
!CMR, 17/4/2012: computes homogeneous solutions for il >= ilmin
!                il >= ilmin includes trapped particles, wfb
!                AND passing particles IF use_pass_homog = T

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

!CMR, 17/4/2012: ng2+1< il <=lmax excludes wfb,ttp 
    if (nlambda > ng2 .and. il > ng2+1 .and. il <= lmax) then
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
! self_periodic bc is applied as follows:
! if boundary_option_linked = .T.
!      isolated wfb (ie no linked domains) 
!      passing + wfb if aky=0
! else (ie boundary_option_linked = .F.)
!      wfb 
!      passing and wfb particles if boundary_option_self_periodic = T
!
!CMR, 6/8/2014:
! Parallel BC for wfb is as follows:
!    ballooning space 
!           wfb ALWAYS self-periodic (-ntgrid:ntgrid)
!    linked fluxtube  
!           wfb self-periodic (-ntgrid:ntgrid) ONLY if ISOLATED, 
!                                               OR if iky=0
!           OTHERWISE wfb treated as passing for now in (-ntgrid:ntgrid)
!                     store homogeneous and inhomog solutions 
!                     AND set amplitude of homog later in invert_rhs_linked 
!                     to satisfy self-periodic bc in fully linked SUPER-CELL

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
       if (use_pass_homog .and. il <= ng2+1) then
          call self_periodic
       else if (il == ng2 + 1) then
          call self_periodic
       end if

    end if

    ! add correct amount of homogeneous solution for trapped particles to satisfy boundary conditions
!CMR, 24/7/2014:  
!Revert to looping from il>= ng2+2, i.e. exclude wfb as: 
!          (1) wfb bc is already handled above
!          (2) forbid never true for wfb, so including ng2+1 in loop is pointless.
    if (il >= ng2+2 .and. il <= lmax) then
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

    if (def_parity) call enforce_parity(parity_redist)

!CMR,DD, 25/7/2014: 
! Not keen on following kludge zeroing forbidden region
! Introduced new flag: zero_forbid defaults to .false.
!  Tested and default is fine linearly, expect should work nonlinearly, 
!  (Can easily restore old behaviour by setting: zero_forbid=.true.)
    ! zero out spurious gnew outside trapped boundary
if(zero_forbid)then
    where (forbid(:,il))
       gnew(:,1,iglo) = 0.0
       gnew(:,2,iglo) = 0.0
    end where
endif
!    call prof_leaving ("invert_rhs_1", "dist_fn")

  contains

    subroutine self_periodic
!CMR: sets sum of homogeneous and inhomogeneous solutions to give a result
!     gnew(ntgr,2) = gnew(ntgl,2)
!     gnew(ntgr,1) = gnew(ntgl,1) 
! ie periodic bcs at the ends of the flux tube.

!CMR, 25/7/2014:
! self-periodic applied to zonal modes (makes sense)
!                       and wfb        (seems strange)
! adding adjr, adjl to cope with application of self-periodic to WFB
!   dadj=adjl-adjr will be ZERO for ky=0 modes, but NOT for WFB.
! This change sets g_wesson (or h) to be self-periodic for wfb, not g !!!
! NB this code change will implement this only in ballooning space, 
! and not in a linked fluxtube.
      implicit none
      complex :: adjl, adjr, dadj

      if (il .eq. ng2+1) then 
         adjl = anon(ie)*2.0*vperp2(ntgl,iglo)*aj1(ntgl,iglo) &
              *bparnew(ntgl,it,ik)*fbpar &
              + spec(is)%z*anon(ie)*phinew(ntgl,it,ik)*aj0(ntgl,iglo) &
              /spec(is)%temp*fphi
         adjr = anon(ie)*2.0*vperp2(ntgr,iglo)*aj1(ntgr,iglo) &
              *bparnew(ntgr,it,ik)*fbpar &
              + spec(is)%z*anon(ie)*phinew(ntgr,it,ik)*aj0(ntgr,iglo) &
              /spec(is)%temp*fphi
         dadj = adjl-adjr
      else 
         dadj = cmplx(0.0,0.0)
      endif

      if (g1(ntgr,1) /= 1.) then
         beta1 = (gnew(ntgr,1,iglo) - gnew(ntgl,1,iglo) - dadj)/(1.0 - g1(ntgr,1))
         gnew(:,1,iglo) = gnew(:,1,iglo) + beta1*g1(:,1)
      end if

      if (g1(ntgl,2) /= 1.) then
         beta1 = (gnew(ntgl,2,iglo) - gnew(ntgr,2,iglo) + dadj)/(1.0 - g1(ntgl,2))
         gnew(:,2,iglo) = gnew(:,2,iglo) + beta1*g1(:,2)
      end if
      
    end subroutine self_periodic

  end subroutine invert_rhs_1

  subroutine invert_rhs_linked &
       (phi, apar, bpar, phinew, aparnew, bparnew, istep, sourcefac)
    use dist_fn_arrays, only: gnew
    use theta_grid, only: bmag, ntgrid
    use le_grids, only: energy, al, nlambda, ng2, anon
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx, idx
    use redistribute, only: fill
    use run_parameters, only: fbpar, fphi
    use species, only: spec
    use spfunc, only: j0, j1
    use kt_grids, only: kperp2

    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    complex, intent (in) :: sourcefac

    complex :: b0, fac, facd
    integer :: il, ik, it, n, i, j
    integer :: iglo, ncell
!
! adding adjr, adjl to cope with application of self-periodic to WFB
!   dadj=adjl-adjr will be ZERO for ky=0 modes, but NOT for WFB.
! This change sets g_wesson (or h) to be self-periodic for wfb, not g !!!
! NB this code change implement this in a linked fluxtube.
!
    integer :: ie, is, itl, itr
    complex :: adjl, adjr, dadj
    real :: vperp2left, vperp2right, argl, argr, aj0l, aj0r, aj1l, aj1r


    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       call invert_rhs_1 (phi, apar, bpar, phinew, aparnew, bparnew, &
            istep, iglo, sourcefac)
    end do

    if (no_comm) then
       ! nothing
    else       
       g_adj = 0. !This shouldn't be needed
       !<DD>Note these fill routines are often equivalent to an all-to-all type
       !communication, i.e. when nproc --> 2 nproc, t_fill --> 4 t_fill
       !By only communicating with our direct neighbours we would significantly
       !reduce the amount of data to communicate and we should improve the communication
       !scaling. However, if we do this then we lose the ability to perform the linked
       !update (i.e. what we do below) in a parallel manner, so the below code would
       !become partially serial.
       call fill (links_p, gnew, g_adj)
       call fill (links_h, g_h, g_adj)
       call fill (wfb_p, gnew, g_adj)
       call fill (wfb_h, g_h, g_adj)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)

          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

          il = il_idx(g_lo,iglo)

! wfb
          if (nlambda > ng2 .and. il == ng2+1) then
             is = is_idx(g_lo,iglo)
             ie = ie_idx(g_lo,iglo)
!CMR, 29/8/2014:
!(1) compute adjl, adjr: the corrections in mapping from g_gs2 to g_wesson
!                     at the extreme left and right of the supercell
!(2) dadj=adjl-adjr then used to apply the self-periodic bc to g_wesson
!    (was previously applied to g)
!      
             vperp2left = bmag(-ntgrid)*al(il)*energy(ie)
             vperp2right = bmag(ntgrid)*al(il)*energy(ie)
             itl=get_leftmost_it(it,ik)
             itr=get_rightmost_it(it,ik)
             argl = spec(is)%bess_fac*spec(is)%smz*sqrt(energy(ie)*al(il)/bmag(-ntgrid)*kperp2(-ntgrid,itl,ik))
             argr = spec(is)%bess_fac*spec(is)%smz*sqrt(energy(ie)*al(il)/bmag(ntgrid)*kperp2(ntgrid,itr,ik))
             aj0l = j0(argl)
             aj0r = j0(argr)
             aj1l = j1(argl)
             aj1r = j1(argr)

             adjl = anon(ie)*2.0*vperp2left*aj1l &
              *bparnew(-ntgrid,itl,ik)*fbpar &
              + spec(is)%z*anon(ie)*phinew(-ntgrid,itl,ik)*aj0l &
              /spec(is)%temp*fphi
             adjr = anon(ie)*2.0*vperp2right*aj1r &
              *bparnew(ntgrid,itr,ik)*fbpar &
              + spec(is)%z*anon(ie)*phinew(ntgrid,itr,ik)*aj0r &
              /spec(is)%temp*fphi
             dadj = adjl-adjr

             if (save_h(1, iglo)) then
!CMR, 7/8/2014:
! This code implements "self-periodic" parallel BC for wfb
! g_adj contains exit point inhomog and homog sol'ns from ALL wfb cells
! g_adj(j,1,iglo):    bcs for rightwards travelling particles (RP)
!  j=1:ncell : inhom sol'n at ntgrid in cell j from RIGHT for RP
!  j=ncell+1:2ncell : hom sol'n at ntgrid in cell (2ncell+1-j) from RIGHT for RP
! g_adj(j,2,iglo):    bcs for leftwards travelling particles (LP)
!  j=1:ncell : inhom sol'n at -ntgrid in cell j from LEFT for LP
!  j=ncell+1:2ncell : hom sol'n at -ntgrid in cell (2ncell+1-j) from LEFT for LP
!
!  facd= 1/(1-\Prod_{j=1,ncell} h^r_j)  (see CMR Parallel BC note)    
!    where h^r_j is the homogeneous exit solution from cell j for RP
!                      
                facd = 1.0
                do j = 1, ncell
                   facd = facd * g_adj(ncell+j,1,iglo)
                end do
                facd = 1./(1.-facd)
                
                b0 = 0.

! i labels cell counting from LEFTmost cell.
                do i = 1, ncell-1
                   fac = 1.0
                   do j = i+1, ncell
! fac is product of all homog sol's from cells to RIGHT on cell i
! g_adj(ncell+j,1,iglo) accesses these homog solutions
                      fac = fac * g_adj(ncell+j,1,iglo)
                   end do
! g_adj(ncell+1-i,1,iglo) accesses inhom solution from cell i
                   b0 = b0 + fac * g_adj(ncell+1-i,1,iglo)
                end do

! b0 computed next line is homog amplitude in leftmost cell  (see CMR note)
                b0 = (b0 + g_adj(1,1,iglo)-dadj)*facd

! BUT we NEED homog amplitude in THIS cell.
! Solve matrix BC equation by cascading homog amplitude, b0, rightwards 
!        from leftmost cell to THIS cell.

                do i = 1, l_links(ik, it)
!  Loop implements cascade to right, using CMR note equ'n:  
!           a^r_{j+1} = a^r_j h^r_j + i^r_j 

                   b0 = b0 * g_adj(ncell+i,1,iglo) + g_adj(ncell+1-i,1,iglo)
                end do

! Resultant b0 is homog amplitude for THIS cell.
! Finally add correct homog amplitude to gnew to match the parallel BC.
                gnew(:,1,iglo) = gnew(:,1,iglo) + b0*g_h(:,1,iglo)
             endif

! Following code repeats same calculation for LEFTWARD travelling wfb.
!CMRend
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
                b0 = (b0 + g_adj(1,2,iglo)+dadj)*facd

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
    use constants, only: zi, pi
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
    !Sourcefac ends up being passed all the way through to get_source_term
    !where it is multiplied by phi_ext (default 0.0) and added to ky<epsilon(0.0)
    !modes source term. Should probably just be calculated in get_source_term and
    !only if min(ky)<epsilon(0.0) & phi_ext/=0 & source_option_switch==source_option_phiext_full
    !
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

  !>Ensure that linked boundary values of passed complex field are single valued (e.g. kperp2(ntgrid,ikx,iky) is
  !!equal to gnew(-ntgrid,ikx_link,iky) as these correspond to the same location).
  subroutine ensure_single_val_fields_pass(arr)
    !Added by <DD>
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use mp, only: broadcast
    implicit none
    integer :: it,ik, link_it
    complex, dimension (-ntgrid:,:,:), intent (in out) :: arr
    if (boundary_option_switch.ne.boundary_option_linked) return
    do ik=1,naky
       do it=1,ntheta0
          link_it=itright(ik,it)
          if (link_it.lt.0) cycle
          arr(-ntgrid,link_it,ik)=arr(ntgrid,it,ik)
       enddo
    enddo
    !This broadcast was used originally but doesn't appear to be required during recent testing
    !This could depend on the mpi library in use
    !    call broadcast(arr)
  end subroutine ensure_single_val_fields_pass

  !>Ensure that linked boundary values of passed complex field are single valued (e.g. kperp2(ntgrid,ikx,iky) is
  !!equal to gnew(-ntgrid,ikx_link,iky) as these correspond to the same location).
  subroutine ensure_single_val_fields_pass_r(arr)
    !Added by <DD>
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use mp, only: broadcast
    implicit none
    integer :: it,ik, link_it
    real, dimension (-ntgrid:,:,:), intent (in out) :: arr
    if (boundary_option_switch.ne.boundary_option_linked) return
    do ik=1,naky
       do it=1,ntheta0
          link_it=itright(ik,it)
          if (link_it.lt.0) cycle
          arr(-ntgrid,link_it,ik)=arr(ntgrid,it,ik)
       enddo
    enddo
    !This broadcast was used originally but doesn't appear to be required during recent testing
    !This could depend on the mpi library in use
    !    call broadcast(arr)
  end subroutine ensure_single_val_fields_pass_r

  subroutine getan (antot, antota, antotp)
    use dist_fn_arrays, only: vpa, vperp2, aj0, aj1, gnew
    use kt_grids, only: kperp2
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

!<DD>
!Don't do this as integrate_species will fill in all values
!    antot=0. ; antota=0. ; antotp=0.
!Need to set individual arrays to zero if not using integrate_species for
!that field. (NOTE this probably isn't actually needed as we shouldn't
!use the various antots if we're not calculating/using the related field).

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
       if(esv) call ensure_single_val_fields_pass(antot) !<DD>
    else
       antot=0.
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
       if(esv) call ensure_single_val_fields_pass(antota) !<DD>
    else
       antota=0.
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
       if(esv) call ensure_single_val_fields_pass(antotp) !<DD>
    else
       antotp=0.
    end if

    call prof_leaving ("getan", "dist_fn")
  end subroutine getan

  subroutine getmoms (phinew,bparnew,ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)
    use dist_fn_arrays, only: vpa, vperp2, aj0, aj1, gnew, g_adjust
    use gs2_layouts, only: is_idx, ie_idx, g_lo, ik_idx, it_idx
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment, anon, energy
    use prof, only: prof_entering, prof_leaving
    use run_parameters, only: fphi, fbpar

    implicit none
    logical, parameter :: full_arr=moment_to_allprocs
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: density, &
         upar, tpar, tperp, ntot, qparflux, pperpj1, qpperpj1
    complex, dimension (-ntgrid:,:,:), intent(in) :: phinew, bparnew

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
    call integrate_moment (g0, density,  moment_to_allprocs, full_arr)

! DJA/CMR: upar and tpar moments 
! (nb adiabatic part of <delta f> does not contribute to upar, tpar or tperp)
! NB UPAR is normalised to vt_s = sqrt(T_s/m_s) vt_ref
!    ie multiply by spec(is)%stm * vt_ref to get abs upar
    g0 = vpa*g0
    call integrate_moment (g0, upar,  moment_to_allprocs, full_arr)

    g0 = 2.*vpa*g0
    call integrate_moment (g0, tpar,  moment_to_allprocs, full_arr)
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
    call integrate_moment (g0, tperp,  moment_to_allprocs, full_arr)
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
          g0(:,isgn,iglo) = vpa(:,isgn,iglo)*gnew(:,isgn,iglo)*aj0(:,iglo)*energy(ie_idx(g_lo,iglo))
       end do
    end do 
    call integrate_moment (g0, qparflux,  moment_to_allprocs, full_arr)
   
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
    call integrate_moment (g0, pperpj1,  moment_to_allprocs, full_arr)

! Now compute QPPERPJ1, a modified p_perp*energy which gives heat flux from Bpar
! NB QPPERPJ1 is normalised to (n_s T_s^2/q_s)  n_ref  T_ref^2/q_ref
!    ie multiply by (n_s T_s spec(is)%tz) n_ref T_ref^2/q_ref to get abs QPPERPJ1
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       g0(:,:,iglo) = g0(:,:,iglo)*energy(ie_idx(g_lo,iglo))
    end do
    call integrate_moment (g0, qpperpj1, moment_to_allprocs, full_arr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now include the adiabatic part of <delta f>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set g0 = <delta_f>/F_m, including the adiabatic term
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo) ; is = is_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo) ; it = it_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo) - anon(ie)*phinew(:,it,ik)*spec(is)%zt
       end do
    end do

! total perturbed density
    call integrate_moment (g0, ntot,  moment_to_allprocs, full_arr)

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

  subroutine getmoms_gryfx (density_gryfx, upar_gryfx, tpar_gryfx, tperp_gryfx, qpar_gryfx, qperp_gryfx, phi_gryfx)
    use dist_fn_arrays, only: vpa, vperp2, aj0, aj1, gnew, g_adjust
    use gs2_layouts, only: is_idx, ie_idx, g_lo, ik_idx, it_idx
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment, anon, energy
    use prof, only: prof_entering, prof_leaving
    use run_parameters, only: fphi, fbpar
    use fields_arrays, only: phinew, bparnew, phi
    use kt_grids, only: ntheta0, naky
    use mp, only: proc0, iproc

    implicit none
    integer :: ik, it, isgn, ie, is, iglo, ig, iz, index_gryfx
    complex*8, dimension(naky*ntheta0*2*ntgrid*nspec), intent(out) :: density_gryfx, upar_gryfx, &
         tpar_gryfx, tperp_gryfx, qpar_gryfx, qperp_gryfx
    complex*8, dimension(naky*ntheta0*2*ntgrid), intent(out) :: phi_gryfx       
    complex, dimension(:,:,:,:), allocatable :: total

    real :: densfac_lin, uparfac_lin, tparfac_lin, tprpfac_lin, qparfac_lin, qprpfac_lin, phifac_lin

    logical :: higher_order_moments

    allocate(total(-ntgrid:ntgrid,ntheta0,naky,nspec))
   

    higher_order_moments = .false.

    ! dens = < f >
    g0 = gnew 
    if(higher_order_moments) g0 = vpa**4.*gnew !rparpar = < vpar^4 f >
    call integrate_moment (g0, total)
    if(proc0) then
    do ig = -ntgrid, ntgrid-1
      iz = ig + ntgrid + 1
      do it = 1,g_lo%ntheta0
        do ik = 1, g_lo%naky
          do is = 1,g_lo%nspec
            index_gryfx = 1 + (ik-1) + g_lo%naky*((it-1)) + &
                            g_lo%naky*g_lo%ntheta0*(iz-1) + &
                            2*ntgrid*g_lo%naky*g_lo%ntheta0*(is-1)
            density_gryfx(index_gryfx) = total(ig, it, ik, is)
          end do
        end do
      end do
    end do
    end if

    ! upar = < (vpar/vti) f > = a/rho_i upar1/vti
    g0 = vpa*g0  ! vpa is normalized to vti, i.e. vpa = vpar/vti
    if(higher_order_moments) then ! rparprp = < vpar^2 vprp^2/2 f >
      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         do isgn = 1, 2
            g0(:,isgn,iglo) = .5*vperp2(:,iglo)*vpa(:,isgn,iglo)*vpa(:,isgn,iglo)*gnew(:,isgn,iglo) 
         end do
      end do
    endif
    call integrate_moment (g0, total)
    if(proc0) then
    do ig = -ntgrid, ntgrid-1
      iz = ig + ntgrid + 1
      do it = 1,g_lo%ntheta0
        do ik = 1, g_lo%naky
          do is = 1,g_lo%nspec
            index_gryfx = 1 + (ik-1) + g_lo%naky*((it-1)) + &
                            g_lo%naky*g_lo%ntheta0*(iz-1) + &
                            2*ntgrid*g_lo%naky*g_lo%ntheta0*(is-1)
            upar_gryfx(index_gryfx) = total(ig, it, ik, is)
          end do
        end do
      end do
    end do
    end if

    ! ppar = < (vpar/vti)^2 f > = a/rho_i ppar1/( n0 m vti^2 )
    g0 = vpa*g0
    if(higher_order_moments) then ! rprpprp = < vprp^4/4 f >
      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         do isgn = 1, 2
            g0(:,isgn,iglo) = .25*vperp2(:,iglo)*vperp2(:,iglo)*gnew(:,isgn,iglo) 
         end do
      end do
    endif
    call integrate_moment (g0, total)
    ! this gives ppar = tpar + dens
    if(proc0) then
    do ig = -ntgrid, ntgrid-1
      iz = ig + ntgrid + 1
      do it = 1,g_lo%ntheta0
        do ik = 1, g_lo%naky
          do is = 1,g_lo%nspec
            index_gryfx = 1 + (ik-1) + g_lo%naky*((it-1)) + &
                            g_lo%naky*g_lo%ntheta0*(iz-1) + &
                            2*ntgrid*g_lo%naky*g_lo%ntheta0*(is-1)
            tpar_gryfx(index_gryfx) = total(ig, it, ik, is)
          end do
        end do
      end do
    end do
    end if

    ! Qpar = < (vpar/vti)^3 f > = a/rho_i Qpar1/( n0 m vti^3 )
    g0 = vpa*g0
    if(higher_order_moments) then ! sparprp = < vpar^3 vprp^2/2 f >
      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         do isgn = 1, 2
            g0(:,isgn,iglo) = .5*vperp2(:,iglo)*vpa(:,isgn,iglo)*vpa(:,isgn,iglo)*vpa(:,isgn,iglo)*gnew(:,isgn,iglo) 
         end do
      end do
    endif
    call integrate_moment (g0, total)
    ! this gives Qpar = qpar + 3*upar
    if(proc0) then
    do ig = -ntgrid, ntgrid-1
      iz = ig + ntgrid + 1
      do it = 1,g_lo%ntheta0
        do ik = 1, g_lo%naky
          do is = 1,g_lo%nspec
            index_gryfx = 1 + (ik-1) + g_lo%naky*((it-1)) + &
                            g_lo%naky*g_lo%ntheta0*(iz-1) + &
                            2*ntgrid*g_lo%naky*g_lo%ntheta0*(is-1)
            qpar_gryfx(index_gryfx) = total(ig, it, ik, is)
          end do
        end do
      end do
    end do
    end if

    ! pprp = < 1/2 (vprp/vti)^2 f > = a/rho_i pprp1/( n0 m vti^2 )
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = .5*vperp2(:,iglo)*gnew(:,isgn,iglo) 
          ! vperp2 is normalized to vti, i.e. vperp2 = (vprp/vti)^2
       end do
    end do
    if(higher_order_moments) g0 = vpa**5.*gnew !sparpar = < vpar^5 f >
    call integrate_moment (g0, total)
    ! this gives pprp = tprp + dens
    if(proc0) then
    do ig = -ntgrid, ntgrid-1
      iz = ig + ntgrid + 1
      do it = 1,g_lo%ntheta0
        do ik = 1, g_lo%naky
          do is = 1,g_lo%nspec
            index_gryfx = 1 + (ik-1) + g_lo%naky*((it-1)) + &
                            g_lo%naky*g_lo%ntheta0*(iz-1) + &
                            2*ntgrid*g_lo%naky*g_lo%ntheta0*(is-1)
            tperp_gryfx(index_gryfx) = total(ig, it, ik, is)
          end do
        end do
      end do
    end do
    end if

    ! Qprp = < 1/2 (vpar/vti) (vprp/vti)^2 f > = a/rho_i Qprp1/( n0 m vti^3 )
    g0 = vpa*g0
    if(higher_order_moments) then ! sprpprp = < vpar vprp^4/4 f >
      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         do isgn = 1, 2
            g0(:,isgn,iglo) = .25*vperp2(:,iglo)*vperp2(:,iglo)*vpa(:,isgn,iglo)*gnew(:,isgn,iglo) 
         end do
      end do
    endif
    call integrate_moment (g0, total)
    ! this gives Qprp = qprp + upar
    if(proc0) then
    do ig = -ntgrid, ntgrid-1
      iz = ig + ntgrid + 1
      do it = 1,g_lo%ntheta0
        do ik = 1, g_lo%naky
          do is = 1,g_lo%nspec
            index_gryfx = 1 + (ik-1) + g_lo%naky*((it-1)) + &
                            g_lo%naky*g_lo%ntheta0*(iz-1) + &
                            2*ntgrid*g_lo%naky*g_lo%ntheta0*(is-1)
            qperp_gryfx(index_gryfx) = total(ig, it, ik, is)
          end do
        end do
      end do
    end do
    end if
 
       densfac_lin = sqrt(2.)           ! ~ rho_i
       uparfac_lin = 2.                 ! ~ rho_i vti
       tparfac_lin = 2.*sqrt(2.)        ! ~ rho_i vti^2
       tprpfac_lin = 2.*sqrt(2.)        ! ~ rho_i vti^2
       qparfac_lin = 4.                 ! ~ rho_i vti^3
       qprpfac_lin = 4.                 ! ~ rho_i vti^3
       phifac_lin = densfac_lin         ! ~ rho_i

    if(higher_order_moments) then 
       densfac_lin = -4.*sqrt(2.)        ! rparpar
       uparfac_lin = -4.*sqrt(2.)        ! rparprp
       tparfac_lin = -4.*sqrt(2.)        ! rprpprp
       tprpfac_lin = 8.                 ! sparpar
       qparfac_lin = 8.                 ! sparprp
       qprpfac_lin = 8.                 ! sprpprp
    endif

    if(proc0 ) then 
      ! normalize to gryfx units
      density_gryfx = densfac_lin*density_gryfx
      upar_gryfx = uparfac_lin*upar_gryfx
      tpar_gryfx = tparfac_lin*tpar_gryfx
      tperp_gryfx = tprpfac_lin*tperp_gryfx
      qpar_gryfx = qparfac_lin*qpar_gryfx
      qperp_gryfx = qprpfac_lin*qperp_gryfx
      if(.not. higher_order_moments) then 
        ! make tpar, tprp, qpar, qprp from ppar, pprp, Qpar, Qprp
        tpar_gryfx = tpar_gryfx - density_gryfx
        tperp_gryfx = tperp_gryfx - density_gryfx
        qpar_gryfx = qpar_gryfx - 3.*upar_gryfx
        qperp_gryfx = qperp_gryfx - upar_gryfx
      endif
    endif
      

    if(proc0) then
    do ig = -ntgrid, ntgrid-1
      iz = ig + ntgrid + 1
      do it = 1, g_lo%ntheta0
        do ik = 1, g_lo%naky
          index_gryfx = 1 + (ik-1) + g_lo%naky*((it-1)) + g_lo%naky*g_lo%ntheta0*(iz-1)
          phi_gryfx(index_gryfx) = phinew(ig, it, ik)*phifac_lin
        end do
      end do
    end do
    end if

    deallocate(total)


  end subroutine getmoms_gryfx


  subroutine getemoms (phinew, bparnew, ntot, tperp)
    use dist_fn_arrays, only: vperp2, aj0, gnew, g_adjust
    use gs2_layouts, only: is_idx, ie_idx, g_lo, ik_idx, it_idx
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment, anon
    use prof, only: prof_entering, prof_leaving
    use run_parameters, only: fphi, fbpar

    implicit none
    logical, parameter :: full_arr=moment_to_allprocs
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: tperp, ntot
    complex, dimension (-ntgrid:,:,:), intent(in) :: phinew, bparnew

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
          g0(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo) - anon(ie)*phinew(:,it,ik)*spec(is)%zt
       end do
    end do

! total perturbed density
    call integrate_moment (g0, ntot,  moment_to_allprocs, full_arr)

! vperp**2 moment:
    do iglo = g_lo%llim_proc, g_lo%ulim_proc       
       ie = ie_idx(g_lo,iglo) ; is = is_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo) ; it = it_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo)*vperp2(:,iglo) - anon(ie)*phinew(:,it,ik)*spec(is)%zt*vperp2(:,iglo)
       end do
    end do

! total perturbed perp pressure
    call integrate_moment (g0, tperp,  moment_to_allprocs, full_arr)

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
  subroutine getmoms_notgc (phinew, bparnew, dens, upar, tpar, tper, ntot, jpar)
    use dist_fn_arrays, only: vpa, vperp2, aj0, aj1, gnew
    use gs2_layouts, only: g_lo, is_idx, ik_idx, it_idx
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use kt_grids, only: nakx => ntheta0, naky
    use le_grids, only: integrate_moment

    implicit none
    logical, parameter :: full_arr=moment_to_allprocs
    complex, intent (out) :: &
         & dens(-ntgrid:,:,:,:), upar(-ntgrid:,:,:,:), &
         & tpar(-ntgrid:,:,:,:), tper(-ntgrid:,:,:,:)
    complex, intent (out), optional :: ntot(-ntgrid:,:,:,:)
    complex, intent (out), optional :: jpar(-ntgrid:,:,:)
    complex, dimension(-ntgrid:,:,:), intent(in) :: phinew, bparnew

    integer :: isgn, iglo, is

    real :: a, b, tpar2, tper2
    integer :: it, ik, ig

! returns moment integrals to PE 0

! not guiding center n_total
    if(present(ntot)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)

          do isgn = 1, 2
             g0(:,isgn,iglo) = aj0(:,iglo) * gnew(:,isgn,iglo)
          end do
          do isgn = 1, 2
             g0(:,isgn,iglo) = g0(:,isgn,iglo) + phinew(:,it,ik) &
                  & *(aj0(:,iglo)**2-1.0) * spec(is)%zt
          end do
          do isgn = 1, 2
             g0(:,isgn,iglo) = g0(:,isgn,iglo) &
                  & + 2.*vperp2(:,iglo)*aj1(:,iglo)*aj0(:,iglo) &
                  & * bparnew(:,it,ik)
          end do
       end do
       call integrate_moment (g0, ntot,  moment_to_allprocs,full_arr)
    endif

! not guiding center density
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g0, dens,  moment_to_allprocs,full_arr)

! not guiding center upar
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = aj0(:,iglo)*vpa(:,isgn,iglo)*gnew(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g0, upar,  moment_to_allprocs,full_arr)

! not guiding center tpar
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = 2.*aj0(:,iglo)*vpa(:,isgn,iglo)**2*gnew(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g0, tpar,  moment_to_allprocs,full_arr)
    
! not guiding center tperp
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = 2.*vperp2(:,iglo)*aj1(:,iglo)*gnew(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g0, tper,  moment_to_allprocs,full_arr)

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
    use dist_fn_arrays, only: aj0, aj1, vperp2
    use species, only: nspec, spec, has_electron_species
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0, aky, kperp2
    use le_grids, only: anon, integrate_species
    use gs2_layouts, only: g_lo, ie_idx, is_idx
    use run_parameters, only: tite
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
          g0(:,isgn,iglo) = (1.0 - aj0(:,iglo)**2)*anon(ie)
       end do
    end do
    wgt = spec%z*spec%z*spec%dens/spec%temp
    call integrate_species (g0, wgt, tot)
    gamtot = real(tot) + kperp2*poisfac
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) = aj0(:,iglo)*aj1(:,iglo) &
               *2.0*vperp2(:,iglo)*anon(ie)
       end do
    end do
    wgt = spec%z*spec%dens
    call integrate_species (g0, wgt, tot)
    gamtot1 = real(tot)
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) = aj1(:,iglo)**2*2.0*vperp2(:,iglo)**2*anon(ie)
       end do
    end do
    wgt = spec%temp*spec%dens
    call integrate_species (g0, wgt, tot)
    gamtot2 = real(tot)

    !<DD>Make sure gamtots are single valued
    if(esv)then
       call ensure_single_val_fields_pass_r(gamtot)
       call ensure_single_val_fields_pass_r(gamtot1)
       call ensure_single_val_fields_pass_r(gamtot2)
    endif

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
    use theta_grid, only: ntgrid, bmag, delthet, jacob
    use kt_grids, only: naky, ntheta0, aky, kperp2
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


!///////////////////////////////////////
!// SOME NO GATHER TEST ROUTINES
!///////////////////////////////////////
  subroutine getfieldeq_nogath (phi, apar, bpar, fieldeq, fieldeqa, fieldeqp)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (out) ::fieldeq,fieldeqa,fieldeqp
    complex, dimension (:,:,:), allocatable :: antot, antota, antotp

    allocate (antot (-ntgrid:ntgrid,ntheta0,naky))
    allocate (antota(-ntgrid:ntgrid,ntheta0,naky))
    allocate (antotp(-ntgrid:ntgrid,ntheta0,naky))

    call getan_nogath (antot, antota, antotp)
    call getfieldeq1_nogath (phi, apar, bpar, antot, antota, antotp, &
         fieldeq, fieldeqa, fieldeqp)

    deallocate (antot, antota, antotp)
  end subroutine getfieldeq_nogath

  subroutine getfieldeq1_nogath (phi, apar, bpar, antot, antota, antotp, &
       fieldeq, fieldeqa, fieldeqp)
    use theta_grid, only: ntgrid, bmag, delthet, jacob
    use kt_grids, only: naky, ntheta0, aky, kperp2
    use run_parameters, only: fphi, fapar, fbpar
    use run_parameters, only: beta, tite
    use kt_grids, only: kwork_filter
    use species, only: spec, has_electron_species
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: antot, antota, antotp
    complex, dimension (-ntgrid:,:,:), intent (out) ::fieldeq,fieldeqa,fieldeqp

    integer :: ik, it
    if (.not. allocated(fl_avg)) allocate (fl_avg(ntheta0, naky))
    if (.not. has_electron_species(spec)) fl_avg = 0.

    if (.not. has_electron_species(spec)) then
       if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
          if (.not. allocated(awgt)) then
             allocate (awgt(ntheta0, naky))
             awgt = 0.
             do ik = 1, naky
                if (aky(ik) > epsilon(0.0)) cycle
                do it = 1, ntheta0
                   if(kwork_filter(it,ik)) cycle
                   awgt(it,ik) = 1.0/sum(delthet*jacob*gamtot3(:,it,ik))
                end do
             end do
          endif
           
          do ik = 1, naky
             do it = 1, ntheta0
                if(kwork_filter(it,ik)) cycle
                fl_avg(it,ik) = tite*sum(delthet*jacob*antot(:,it,ik)/gamtot(:,it,ik))*awgt(it,ik)
             end do
          end do
       end if
    end if

    if (fphi > epsilon(0.0)) then
!       fieldeq = antot + bpar*gamtot1 - gamtot*gridfac1*phi 
       fieldeq = antot - gamtot*gridfac1*phi 
       if(fbpar.gt.epsilon(0.0)) fieldeq=fieldeq + bpar*gamtot1

       if (.not. has_electron_species(spec)) then
          do ik = 1, naky
             do it = 1, ntheta0
                if(kwork_filter(it,ik)) cycle
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
             if(kwork_filter(it,ik)) cycle
             fieldeqp(:,it,ik) = fieldeqp(:,it,ik)/bmag(:)**2
          end do
       end do
       fieldeqp = fieldeqp + bpar*gridfac1
    end if
  end subroutine getfieldeq1_nogath

  subroutine getan_nogath (antot, antota, antotp)
    use dist_fn_arrays, only: vpa, vperp2, aj0, aj1, gnew
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_species
    use run_parameters, only: beta, fphi, fapar, fbpar
    use gs2_layouts, only: g_lo, it_idx,ik_idx
    use kt_grids, only: kwork_filter, kperp2

    implicit none

    complex, dimension (-ntgrid:,:,:), intent (out) :: antot, antota, antotp
    real, dimension (nspec) :: wgt
    integer :: isgn, iglo, it,ik

    if (fphi > epsilon(0.0)) then
       !<DD>NOTE: It's possible to rewrite this loop as simply
       !g0=gnew*spread(aj0,2,2)
       !but this seems to be slower than an explicit loop and 
       !the ability to skip certain it/ik values is lost.
       if(any(kwork_filter))then
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             it=it_idx(g_lo,iglo)
             ik=ik_idx(g_lo,iglo)
             if(kwork_filter(it,ik))cycle
             do isgn = 1, 2
                g0(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo)
             end do
          end do
       else
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             do isgn = 1, 2
                g0(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo)
             end do
          end do
       endif

       wgt = spec%z*spec%dens
       call integrate_species (g0, wgt, antot, nogath=.true.)

       if (afilter > epsilon(0.0)) antot = antot * exp(-afilter**4*kperp2**2/4.)
       !NOTE: We don't do ensure_single_val_fields here as we're not certain we
       !have the full data
    end if

    if (fapar > epsilon(0.0)) then
       if(any(kwork_filter))then
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             it=it_idx(g_lo,iglo)
             ik=ik_idx(g_lo,iglo)
             if(kwork_filter(it,ik))cycle
             do isgn = 1, 2
                   g0(:,isgn,iglo) = aj0(:,iglo)*vpa(:,isgn,iglo)*gnew(:,isgn,iglo)
             end do
          end do
       else
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             do isgn = 1, 2
                   g0(:,isgn,iglo) = aj0(:,iglo)*vpa(:,isgn,iglo)*gnew(:,isgn,iglo)
             end do
          end do
       endif

       wgt = 2.0*beta*spec%z*spec%dens*spec%stm
       call integrate_species (g0, wgt, antota, nogath=.true.)
       !NOTE: We don't do ensure_single_val_fields here as we're not certain we
       !have the full data
    end if

    if (fbpar > epsilon(0.0)) then
       if(any(kwork_filter))then
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             it=it_idx(g_lo,iglo)
             ik=ik_idx(g_lo,iglo)
             if(kwork_filter(it,ik))cycle
             do isgn = 1, 2
                   g0(:,isgn,iglo) = aj1(:,iglo)*vperp2(:,iglo)*gnew(:,isgn,iglo)
             end do
          end do
       else
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             do isgn = 1, 2
                   g0(:,isgn,iglo) = aj1(:,iglo)*vperp2(:,iglo)*gnew(:,isgn,iglo)
             end do
          end do
       endif

       wgt = spec%temp*spec%dens
       call integrate_species (g0, wgt, antotp, nogath=.true.)
       !NOTE: We don't do ensure_single_val_fields here as we're not certain we
       !have the full data
    end if

  end subroutine getan_nogath

!///////////////////////////////////////
!///////////////////////////////////////
!///////////////////////////////////////

  
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
    use kt_grids, only: ntheta0, naky, kperp2
    
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
          !NOTE: denominator=0 for the it==ik==1 only in certain circumstances
          !for example a simulation with beta=0.0 and an adiabatic species does
          !not have a zero denominator.
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
       pbflux, qbflux, vbflux, pflux_tormom)

!CMR, 15/1/08: 
!  Implemented Clemente Angioni's fix for fluxes by replacing g with gnew 
!  so fields and distribution function are evaluated self-consistently in time.
!  This fixed unphysical oscillations in non-ambipolar particle fluxes 
!
    use species, only: spec
    use theta_grid, only: ntgrid, bmag, gradpar, delthet
    use theta_grid, only: qval, shat, gds21, gds22
    use kt_grids, only: naky, ntheta0, theta0, aky
!   use theta_grid, only: drhodpsi, grho
!    use kt_grids, only: akx
    use le_grids, only: energy
    use dist_fn_arrays, only: gnew, aj0, vpac, vpa, aj1, vperp2
    use gs2_layouts, only: g_lo, ie_idx, is_idx, it_idx, ik_idx
    use run_parameters, only: woutunits, fphi, fapar, fbpar
    use constants, only: zi
    use geometry, only: rhoc!Should this not be value from theta_grid_params?
    use theta_grid, only: Rplot, Bpol
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    real, dimension (:,:,:), intent (out) :: pflux, pmflux, pbflux, pflux_tormom
    real, dimension (:,:,:), intent (out) :: vflux, vmflux, vbflux, vflux_par, vflux_perp
    real, dimension (:,:,:,:), intent (out) :: qflux, qmflux, qbflux
    real, dimension (:,:,:), allocatable :: dnorm
    integer :: it, ik, is, isgn, ig
    integer :: iglo

    allocate (dnorm (-ntgrid:ntgrid,ntheta0,naky))

    ! EGH for new parallel I/O everyone
    ! calculates fluxes
    !if (proc0) then
       pflux = 0.0;   qflux = 0.0;   vflux = 0.0 ; vflux_par = 0.0 ; vflux_perp = 0.0
       pmflux = 0.0;  qmflux = 0.0;  vmflux = 0.0
       pbflux = 0.0;  qbflux = 0.0;  vbflux = 0.0
       pflux_tormom = 0.0 
    !end if

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
       call get_flux_tormom (phi, pflux_tormom, dnorm)  !JPL

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          g0(:,:,iglo) = g0(:,:,iglo)*energy(ie_idx(g_lo,iglo))
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
          g0(:,:,iglo) = g0(:,:,iglo)*energy(ie_idx(g_lo,iglo))
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
          g0(:,:,iglo) = g0(:,:,iglo)*energy(ie_idx(g_lo,iglo))
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
       vbflux = 0.
    end if

    deallocate (dnorm)
  end subroutine flux

  subroutine get_flux (fld, flx, dnorm)
    use theta_grid, only: ntgrid, grho
    use kt_grids, only: ntheta0, aky, naky
    use le_grids, only: integrate_moment
    use species, only: nspec
    implicit none
    logical, parameter :: full_arr=moment_to_allprocs
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (:,:,:), intent (out) :: flx
    real, dimension (-ntgrid:,:,:) :: dnorm
    complex, dimension (:,:,:,:), allocatable :: total
    real :: wgt
    integer :: ik, it, is

    allocate (total(-ntgrid:ntgrid,ntheta0,naky,nspec))
    ! EGH added final parameter 'all'
    ! to the call below for new parallel output
    ! This is temporary until distributed fields
    ! are implemented 1/2014
    ! DD added full_arr=.true. to ensure all procs get the 
    ! full array
    call integrate_moment (g0, total, moment_to_allprocs, full_arr)

    ! EGH for new parallel I/O everyone
    ! calculates fluxes
    !if (proc0) then
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
    ! end if

    deallocate (total)

  end subroutine get_flux

  subroutine get_flux_tormom (fld, flx, dnorm) !JPL
    use theta_grid, only: ntgrid
#ifdef LOWFLOW
    use kt_grids, only: ntheta0, naky
    use le_grids, only: integrate_moment
    use species, only: nspec
    use mp, only: proc0
    use theta_grid, only: rplot, grho
    use kt_grids, only: aky
    use lowflow, only: mach_lab
#endif
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (:,:,:), intent (out) :: flx
    real, dimension (-ntgrid:,:,:) :: dnorm
!<DD> Moved directive to here (was previously just around flx= line below)
!     as if we haven't compiled with LOWFLOW then this routine shouldn't
!     do anything (integrate_moment is not that cheap).
#ifdef LOWFLOW 
    complex, dimension (:,:,:,:), allocatable :: total
    real :: wgt
    integer :: ik, it, is
    allocate (total(-ntgrid:ntgrid,ntheta0,naky,nspec))
    call integrate_moment (g0, total)

    if (proc0) then
       do is = 1, nspec
          do ik = 1, naky
             do it = 1, ntheta0
                wgt = sum(dnorm(:,it,ik)*grho)
                flx(it,ik,is) = sum(aimag(total(:,it,ik,is)*conjg(fld(:,it,ik))) &
                     *dnorm(:,it,ik)*aky(ik)*Rplot(:)**2)/wgt*mach_lab(is)
             end do
          end do
       end do
       flx = flx*0.5
    end if

    deallocate (total)
#endif
  end subroutine get_flux_tormom !JPL

  subroutine eexchange (phinew, phi, exchange1, exchange2)

    use mp, only: proc0
    use constants, only: zi
    use gs2_layouts, only: g_lo, il_idx, ie_idx, it_idx, ik_idx, is_idx
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: gnew, aj0, vpac, g
    use theta_grid, only: ntgrid, gradpar, delthet, jacob
    use kt_grids, only: ntheta0, naky
    use le_grids, only: integrate_moment
    use run_parameters, only: fphi
    use species, only: spec, nspec
    use nonlinear_terms, only: nonlin
    use hyper, only: hypervisc_filter

    implicit none

    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, phi
    real, dimension (:,:,:), intent (out) :: exchange1, exchange2

    integer :: ig, il, ie, it, ik, is, iglo, isgn
    real :: wgt
    complex :: dgdt_hypervisc
    real, dimension (:,:,:), allocatable :: dnorm
    complex, dimension (:,:,:,:), allocatable :: total, total2

    allocate (dnorm(-ntgrid:ntgrid, ntheta0, naky)) ; dnorm = 0.0
    allocate (total(-ntgrid:ntgrid, ntheta0, naky, nspec)) ; total = 0.0
    allocate (total2(-ntgrid:ntgrid, ntheta0, naky, nspec)) ; total2 = 0.0

    if (proc0) then
       exchange1 = 0.0 ; exchange2 = 0.0
    end if

    do ik = 1, naky
       do it = 1, ntheta0
          dnorm(:ntgrid-1,it,ik) = delthet(:ntgrid-1)*(jacob(-ntgrid+1:)+jacob(:ntgrid-1))*0.5
       end do
    end do

    if (fphi > epsilon(0.0)) then
       g0 = 0.
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)          
          if (nonlin .and. it==1 .and. ik==1) cycle
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          do isgn = 1, 2
             ! get v_magnetic piece of g0 at grid points instead of cell centers
             do ig = -ntgrid, ntgrid
                g0(ig,isgn,iglo) = aj0(ig,iglo)*(zi*wdrift_func(ig,il,ie,it,ik)/code_dt &
                     * gnew(ig,isgn,iglo)*spec(is)%tz)
             end do

             ! add contribution to g0 from hyperviscosity at grid points
             ! this is -(dg/dt)_hypervisc, equivalent to collisions term in eqn. 5 of PRL 109, 185003
             do ig = -ntgrid, ntgrid
                if (abs(hypervisc_filter(ig,it,ik)-1.0) > epsilon(0.0)) then
                   dgdt_hypervisc = (1.0-1./hypervisc_filter(ig,it,ik))*gnew(ig,isgn,iglo)/code_dt
                   ! should gnew be (gnew+gold)/2?
                   g0(ig,isgn,iglo) = g0(ig,isgn,iglo) - dgdt_hypervisc
                end if
             end do

             ! get v_magnetic piece of g0 at cell centers and add in vpar piece at cell centers
             do ig = -ntgrid, ntgrid-1
                g0(ig,isgn,iglo) = 0.5*(g0(ig,isgn,iglo)+g0(ig+1,isgn,iglo)) &
                     + 0.5*vpac(ig,isgn,iglo)*(gradpar(ig)+gradpar(ig+1))/delthet(ig) &
                     * (gnew(ig+1,isgn,iglo)-gnew(ig,isgn,iglo))*spec(is)%stm
             end do

          end do
       end do

       call integrate_moment (g0, total)

       if (proc0) then
          do is = 1, nspec
             do ik = 1, naky
                do it = 1, ntheta0
                   if (nonlin .and. it==1 .and. ik==1) cycle
                   wgt = sum(dnorm(:,it,ik))
                   exchange1(it,ik,is) = sum(real(total(:,it,ik,is)*conjg(phinew(:,it,ik))) &
                        *dnorm(:,it,ik))/wgt
                end do
             end do
          end do
       end if

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             g0(:,isgn,iglo) = aj0(:,iglo)*0.25*(gnew(:,isgn,iglo)+g(:,isgn,iglo))
          end do
       end do
       call integrate_moment (g0, total)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             g0(:,isgn,iglo) = aj0(:,iglo)*0.25*(gnew(:,isgn,iglo)-g(:,isgn,iglo))
          end do
       end do
       call integrate_moment (g0, total2)

       ! exchange2 is a symmetrized form of energy exchange,
       ! which guarantees species-summed energy exchange is zero
       if (proc0) then
          do is = 1, nspec
             do ik = 1, naky
                do it = 1, ntheta0
                   if (nonlin .and. it==1 .and. ik==1) cycle
                   wgt = sum(dnorm(:,it,ik))*code_dt
                   exchange2(it,ik,is) = sum(real(total(:,it,ik,is) &
                        *conjg(phinew(:,it,ik)-phi(:,it,ik)) &
                        - (phinew(:,it,ik)+phi(:,it,ik))*conjg(total2(:,it,ik,is))) &
                        *dnorm(:,it,ik))/wgt
                end do
             end do
          end do
       end if

    end if

    deallocate (dnorm, total, total2)

  end subroutine eexchange

  subroutine lf_flux (phi, vflx0, vflx1)

    use species, only: nspec
    use theta_grid, only: ntgrid, bmag, gradpar, delthet
    use theta_grid, only: drhodpsi, IoB
    use kt_grids, only: naky, ntheta0
    use dist_fn_arrays, only: gnew, aj0, vpa
!    use dist_fn_arrays, only: vpac
    use gs2_layouts, only: g_lo, ie_idx, is_idx, it_idx, ik_idx
    use mp, only: proc0
    use run_parameters, only: woutunits, fphi, rhostar
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi
    real, dimension (:,:,:), intent (out) :: vflx0, vflx1
    real, dimension (:,:), allocatable :: dum
    real, dimension (:,:,:), allocatable :: dnorm
    complex, dimension (:,:,:), allocatable :: dphi
    integer :: it, ik, isgn, ig
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
    use theta_grid, only: ntgrid, grho
    use kt_grids, only: ntheta0, naky
    use le_grids, only: integrate_moment
    use species, only: nspec
    use mp, only: proc0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (:,:,:), intent (out) :: flx
    real, dimension (-ntgrid:,:,:) :: dnorm
    complex, dimension (:,:,:,:), allocatable :: total
    real :: wgt
    integer :: ik, it, is

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

  subroutine flux_vs_theta_vs_vpa (phinew,vflx)
    use constants, only: zi
    use dist_fn_arrays, only: gnew, vperp2, aj1, aj0, vpac
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: it_idx, ik_idx, is_idx
    use geometry, only: rhoc!Should this not be value from theta_grid_params?
    use theta_grid, only: ntgrid, bmag, gds21, gds22, qval, shat
    use theta_grid, only: Rplot, Bpol
    use kt_grids, only: aky, theta0
    use le_grids, only: integrate_volume, nlambda, negrid
    use le_grids, only: get_flux_vs_theta_vs_vpa
    use species, only: spec, nspec

    implicit none
    complex, dimension(-ntgrid:,:,:), intent(in) :: phinew
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
! JPL: Diagnose particle flux contribution to toroidal momentum flux 
!      in the lab frame in terms of vpar and theta. 
!=============================================================================

  subroutine pflux_vs_theta_vs_vpa (vflx)
#ifdef LOWFLOW   
    use dist_fn_arrays, only: gnew, aj0
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: it_idx, ik_idx, is_idx
    use theta_grid, only: Rplot !JPL
    use fields_arrays, only: phinew
    use kt_grids, only: aky
    use le_grids, only: integrate_volume, nlambda, negrid
    use le_grids, only: get_flux_vs_theta_vs_vpa
    use species, only:  nspec
    use lowflow, only: mach_lab    
#endif
    use theta_grid, only: ntgrid
    implicit none

    real, dimension (-ntgrid:,:,:), intent (out) :: vflx
    !This whole routine will return zero when not using lowflow so move directive to here <DD>
#ifdef LOWFLOW  
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
             g0(ig,isgn,iglo) = gnew(ig,isgn,iglo)*aj0(ig,iglo) 

             g0r(ig,isgn,iglo) = aimag(g0(ig,isgn,iglo)*conjg(phinew(ig,it,ik)))*aky(ik)*Rplot(ig)**2*mach_lab(is)

          end do
       end do
    end do
    
    call integrate_volume (g0r, gavg, all)
    call get_flux_vs_theta_vs_vpa (gavg, vflx)

    deallocate (gavg)
    deallocate (g0r)
#else
    vflx=0.
#endif
  end subroutine pflux_vs_theta_vs_vpa


!=============================================================================
! Density: Calculate Density perturbations
!=============================================================================
  !NOTE... this routine is deprecated... its contents have been moved to the 
  ! new diagnostics module and this routine is only needed by the old
  ! diagnostics module. It is thus scheduled for deletion at some point, 
  ! unless someone pleads for its life :-) . EGH.
  subroutine get_jext(j_ext)
    use kt_grids, only: ntheta0, naky,aky, kperp2
!    use mp, only: proc0
    use theta_grid, only: jacob, delthet, ntgrid
    use antenna, only: antenna_apar
    implicit none
    !Passed
    !Intent(in) only needed as initialised to zero at top level
    real, dimension(:,:), intent(in out) ::  j_ext 
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

!    if (proc0) then
       deallocate (wgt)
!    endif
    deallocate (j_extz)

  end subroutine get_jext
!<GGH
!=============================================================================
  subroutine get_heat (h, hk, phi, apar, bpar, phinew, aparnew, bparnew)
    use mp, only: proc0
    use constants, only: zi
    use kt_grids, only: ntheta0, naky, aky, akx, kperp2
#ifdef LOWFLOW
    use dist_fn_arrays, only: hneoc
#endif
    use dist_fn_arrays, only: vpac, aj0, aj1, vperp2, g, gnew, g_adjust
    use gs2_heating, only: heating_diagnostics
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx, ie_idx
    use le_grids, only: integrate_moment
    use species, only: spec, nspec,has_electron_species
    use theta_grid, only: jacob, delthet, ntgrid
    use run_parameters, only: fphi, fapar, fbpar, tunits, beta, tite
    use gs2_time, only: code_dt
    use nonlinear_terms, only: nonlin
    use antenna, only: antenna_apar, a_ext_data
    use hyper, only: D_v, D_eta, nexp, hypervisc_filter

    implicit none
    type (heating_diagnostics), intent(in out) :: h
    type (heating_diagnostics), dimension(:,:), intent(in out) :: hk
!    complex, dimension (-ntgrid:,:,:), pointer :: hh, hnew
    complex, dimension (-ntgrid:,:,:), intent(in) :: phi, apar, bpar, phinew, aparnew, bparnew
    complex, dimension(:,:,:,:), allocatable :: tot
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
                  (chidot*spec(is)%z-hdot*spec(is)%temp)

          end do
       end do
    end do

    deallocate (phidot, apardot, bpardot)

    allocate (tot(-ntgrid:ntgrid, ntheta0, naky, nspec))

    call integrate_moment (g0, tot)

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
          h  % antenna = 0.
          hk % energy_dot = 0.
          hk % energy = 0.
          hk % eapar = 0.
          h  % eapar = 0.
          hk % ebpar = 0.
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
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       do isgn=1,2

          do ig=-ntgrid, ntgrid-1
             
             havg = favg (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo))

             g0(ig,isgn,iglo) = 0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg)
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

                   !hs2 = int_r int_v T/F0 hs^2/2
                   hk(it,ik) % hs2(is) = hk(it,ik) % hs2(is)  &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2

                end do
             end do
             h % hs2(:) = h % hs2(:) + hk(it,ik) % hs2(:)
          end do
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

             g0(ig,isgn,iglo) = 0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg &
                  + conjg(phi_avg)*phi_avg &
                  - conjg(j0phiavg)*havg &
                  - conjg(havg)*j0phiavg)
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
       deallocate (wgt)
    end if

    deallocate (tot)

!!
!! Put g, gnew back to their usual meanings
!!
    call g_adjust (g,    phi,    bpar,    -fphi, -fbpar)
    call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar)

  end subroutine get_heat

  subroutine reset_init
    call finish_dist_fn_level_3

    !use dist_fn_arrays, only: gnew, g
    !initializing  = .true.
    !initialized = .false.
   ! 
    !wdrift = 0.
    !wdriftttp = 0.
    !a = 0.
    !b = 0.
    !r = 0.
    !ainv = 0.
    !gnew = 0.
    !g0 = 0.
    !g = 0.

  end subroutine reset_init

  subroutine reset_physics
    call init_wstar
  end subroutine reset_physics

  subroutine get_verr (errest, erridx, phi, bpar)
     
    use mp, only: broadcast
    use le_grids, only: integrate_species
    use le_grids, only: eint_error, trap_error, lint_error, wdim
    use le_grids, only: ng2, nlambda, new_trap_int
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
    logical :: trap_flag=.true. !Value doesn't really matter as just used in present() checks

    allocate(wgt(nspec))
    allocate(errtmp(2))
    allocate(idxtmp(3))

    !This should probably be something like if(first) then ; call broadcast ; first=.false. ; endif
    !as we only deallocate kmax during finish_dist_fn.
    !Really wdim should just be broadcast when we first calculate it -- currently slightly complicated
    !as calculation triggered in gs2_diagnostics and only called by proc0.
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

  subroutine get_gtran (geavg, glavg, gtavg, phi, bpar)

    use le_grids, only: legendre_transform, nlambda, ng2, nesub, jend
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use dist_fn_arrays, only: gnew, aj0, g_adjust
    use run_parameters, only: fphi, fbpar
    use gs2_layouts, only: g_lo
    use mp, only: proc0, broadcast

    implicit none

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
       call legendre_transform (g0, getran, gltran, gttran)
    else
       call legendre_transform (g0, getran, gltran)
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
!<DD>WARNING: THIS ROUTINE DOESN'T USE g0 FOR ANYTHING IS THIS CORRECT?
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
    use theta_grid, only: bmag
    use dist_fn_arrays, only: gnew

    integer :: iglo, ik, it, is
    integer :: ie, il, ig
    integer, save :: unit
    real :: vpa, vpe
    complex, dimension(2) :: gtmp
!    logical :: first = .true.
    logical, intent(in)  :: last 

    real, dimension (:), allocatable, save :: xpts, ypts

    if (.not. allocated(xpts)) then
       allocate(xpts(negrid))
       allocate(ypts(nlambda))

    ! should really just get rid of xpts and ypts
!    if (first) then
       xpts(1:negrid-1) = zeroes
       xpts(negrid) = x0
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
             vpa = sqrt(energy(ie)*max(0.0, 1.0-al(il)*bmag(ig)))
             vpe = sqrt(energy(ie)*al(il)*bmag(ig))
             write (unit, "(8(1x,e13.6))") vpa, vpe, energy(ie), al(il), &
                  xpts(ie), ypts(il), real(gtmp(1)), real(gtmp(2))
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
    use theta_grid, only: bmag, ntgrid
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
                   write (unit, "(6(1x,e13.6))") energy(ie), al(il), &
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
                   write (unit, "(4(1x,e13.6),i8)") energy(ie), al(il), &
                        gp0, gp0zf, isign
                end if
             end if
          end if
          call barrier
       end do
       deallocate (grs, gzf)
    end if

    if (proc0) call flush_output_file (unit)

    if (proc0) then
       write(unit,*)
       if (last) call close_output_file (unit)
    end if

  end subroutine write_fyx

  !> This routine is now only needed by the old diagnostics
  !! module and is scheduled for deletion. Its functionality
  !! is now in diagnostics/diagnostics_velocity_space.f90
  subroutine collision_error (phi, bpar, last)
    
    use mp, only: proc0, send, receive, barrier
    use le_grids, only: ng2, jend, nlambda, lambda_map
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: gnew, g_adjust
    use run_parameters, only: fphi, fbpar
    use gs2_layouts, only: lz_lo, ig_idx, idx_local, proc_id
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
       write(unit,"((1x,e13.6),6(i8),2(1x,e13.6))") time, &
            igmax, ikmax, itmax, iemax, ilmax, ismax, emax, eavg
       if (last) then
          call close_output_file (unit)
       end if
    end if

    deallocate (lcoll, fdcoll, glze, ltmp, ftmp)
    
  end subroutine collision_error

  subroutine boundary(linked)
    implicit none
    logical, intent(out) :: linked
    call init_dist_fn
    linked = boundary_option_switch == boundary_option_linked
  end subroutine boundary

! This subroutine only returns epar correctly for linear runs.
  subroutine get_epar (phi, apar, phinew, aparnew, epar)

    use theta_grid, only: ntgrid, delthet, gradpar
    use run_parameters, only: tunits, fphi, fapar
    use gs2_time, only: code_dt
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension(-ntgrid:,:,:), intent(in) :: phi, apar, phinew, aparnew
    complex, dimension(-ntgrid:,:,:), intent(out) :: epar
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
    use gs2_layouts, only: it_idx,ik_idx,g_lo,il_idx,ie_idx,is_idx,idx,proc_id
    implicit none
    integer, intent (in) :: iglo
    integer, intent (out) :: iglo_left, ipleft
    integer :: iglo_star
    integer :: it_cur,ik,it,il,ie,is
    iglo_star = iglo
    it_cur=it_idx(g_lo,iglo)
    it=it_cur
    ik=ik_idx(g_lo,iglo)

    !Now get the leftmost it
    it_cur=get_leftmost_it(it,ik)

    !If we're at the same it then don't need to do much
    if(it.eq.it_cur)then
       iglo_left=iglo
       ipleft=proc_id(g_lo,iglo)
       return
    endif

    !If not then we need to calculate iglo_left and ipleft
    il=il_idx(g_lo,iglo)
    ie=ie_idx(g_lo,iglo)
    is=is_idx(g_lo,iglo)
    iglo_left=idx(g_lo,ik,it_cur,il,ie,is)
    ipleft=proc_id(g_lo,iglo_left)

  end subroutine find_leftmost_link

  subroutine find_rightmost_link (iglo, iglo_right, ipright)
    use gs2_layouts, only: it_idx,ik_idx,g_lo,il_idx,ie_idx,is_idx,idx,proc_id
    implicit none
    integer, intent (in) :: iglo
    integer, intent (out) :: iglo_right, ipright
    integer :: iglo_star
    integer :: it_cur,ik,it,il,ie,is
    iglo_star = iglo
    it_cur=it_idx(g_lo,iglo)
    ik=ik_idx(g_lo,iglo)
    it=it_cur

    !Now get the rightmost it
    it_cur=get_rightmost_it(it,ik)

    !If we're at the same it then don't need to do much
    if(it.eq.it_cur)then
       iglo_right=iglo
       ipright=proc_id(g_lo,iglo)
       return
    endif

    !If not then we need to calculate iglo_left and ipleft
    il=il_idx(g_lo,iglo)
    ie=ie_idx(g_lo,iglo)
    is=is_idx(g_lo,iglo)
    iglo_right=idx(g_lo,ik,it_cur,il,ie,is)
    ipright=proc_id(g_lo,iglo_right)

  end subroutine find_rightmost_link

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
    use kt_grids, only: nakx => ntheta0, naky, kperp2
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: aj0,aj1,vperp2
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
          call integrate_moment(gtmp,coeff0,.true.,full_arr=.true.)
          where(real(coeff0(0,1:nakx,1:naky,1:nspec)) == 0.)
             mom_coeff(1:nakx,1:naky,1:nspec,i)=1.
          elsewhere
             mom_coeff(1:nakx,1:naky,1:nspec,i)= &
                  & coeff0(0,1:nakx,1:naky,1:nspec)
          end where
       end do
    endif

    !<DD>Currently below could include divide by zero if analytical=.true.
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


#ifdef LOWFLOW
  ! This rouinte must be called after init_collisions as it relies on the 
  ! variable use_le_layout which is set in init_collisions based on input 
  ! from the user provided input file.
  subroutine init_lowflow

    use constants, only: zi
    use dist_fn_arrays, only: vparterm, wdfac, vpac, wdttpfac
    use dist_fn_arrays, only: wstarfac, hneoc, vpar
    use species, only: spec, nspec
    use geometry, only: rhoc!Should this not be value from theta_grid_params?
    use theta_grid, only: theta, ntgrid, delthet, gradpar, bmag
    use theta_grid, only: gds23, gds24, gds24_noq, cvdrift_th, gbdrift_th
    use theta_grid, only: drhodpsi, qval, shat
    use le_grids, only: energy, al, negrid, nlambda, forbid, init_map
    use le_grids, only: get_flux_vs_theta_vs_vpa
    use kt_grids, only: theta0, ntheta0, naky, aky, akx
    use gs2_time, only: code_dt
    use gs2_layouts, only: g_lo, ik_idx, il_idx, ie_idx, is_idx, it_idx
    use run_parameters, only: tunits, wunits, rhostar, neo_test
    use lowflow, only: get_lowflow_terms
    use file_utils, only: open_output_file, close_output_file
    use collisions, only: use_le_layout
    use mp, only: proc0, mp_abort

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
    if (neo_test) call mp_abort('stopping as neo_test=.true.')
    
    ! intialize mappings from g_lo to e_lo and lz_lo or to le_lo to facilitate
    ! energy and lambda derivatives in parallel nonlinearity, etc.
    if(use_le_layout) then
      call init_map (use_lz_layout=.false., use_e_layout=.false., use_le_layout=.true., test=.false.)
    else
      call init_map (use_lz_layout=.true., use_e_layout=.true., use_le_layout=.false., test=.false.)
    end if

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

  subroutine init_pass_wfb
! CMR, 21/5/2014: 
!       simple experimental new routine sets up "pass_wfb" redist_type to pass 
!       g(ntgrid,1,iglo)  => g(-ntgrid,1,iglo_right) in right connected cell.
!       g(-ntgrid,2,iglo) => g(ntgrid,1,iglo_right) in left connected cell.
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
    integer :: il, iglo, ip, iglo_left, iglo_right, ipleft, ipright, n, nn_max
    logical :: haslinks=.true.
    logical :: debug=.false.

      if (boundary_option_switch .eq. boundary_option_linked) then
! Need communications to satisfy || boundary conditions
! First find required blocksizes 
         nn_to = 0   ! # communicates from ip TO HERE (iproc)
         nn_from = 0 ! # communicates to ip FROM HERE (iproc)
         do iglo = g_lo%llim_world, g_lo%ulim_world
            il = il_idx(g_lo,iglo)          
            if (il /= ng2+1) cycle     ! ONLY consider wfb 
            ip = proc_id(g_lo,iglo)
! First consider rightward connections 
            call get_right_connection (iglo, iglo_right, ipright)
            if (ipright .eq. iproc ) then
               nn_to(ip)=nn_to(ip)+1
            endif
            if (ip .eq. iproc .and. ipright .ge. 0 ) then
               nn_from(ipright)=nn_from(ipright)+1
            endif
! Then leftward connections 
            call get_left_connection (iglo, iglo_left, ipleft)
            if (ipleft .eq. iproc ) then
               nn_to(ip)=nn_to(ip)+1
            endif
            if (ip .eq. iproc .and. ipleft .ge. 0 ) then
               nn_from(ipleft)=nn_from(ipleft)+1
            endif
         end do
         nn_max = maxval(nn_to+nn_from)
         call max_allreduce (nn_max)
! skip addressing if no links are needed
         if (nn_max == 0) then 
            haslinks=.false.
         endif
      endif
      if (debug) then
         write(6,*) 'init_pass_wfb processor, nn_to:',iproc,nn_to
         write(6,*) 'init_pass_wfb processor, nn_from:',iproc,nn_from
      endif

      if (boundary_option_switch.eq.boundary_option_linked .and. haslinks) then
! 
! CMR, 12/11/2013: 
!      communication required to satisfy linked BC for wfb
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
            if (il /= ng2+1) cycle     ! ONLY consider wfb 
            ip = proc_id(g_lo,iglo)
! First consider rightward connections 
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
! Then leftward connections 
            call get_left_connection (iglo, iglo_left, ipleft)
            if (ip .eq. iproc .and. ipleft .ge. 0 ) then
               n=nn_from(ipleft)+1 ; nn_from(ipleft)=n
               from(ipleft)%first(n)=-ntgrid
               from(ipleft)%second(n)=2
               from(ipleft)%third(n)=iglo
            endif
            if (ipleft .eq. iproc ) then
               n=nn_to(ip)+1 ; nn_to(ip)=n
               to(ip)%first(n)=ntgrid
               to(ip)%second(n)=2
               to(ip)%third(n)=iglo_left
            endif
         end do
         if (debug) then
            write(6,*) 'init_pass_wfb processor, nn_to:',iproc,nn_to
            write(6,*) 'init_pass_wfb processor, nn_from:',iproc,nn_from
         endif

         from_low(1)=-ntgrid ; from_low(2)=1 ; from_low(3)=g_lo%llim_proc       
         from_high(1)=ntgrid; from_high(2)=2; from_high(3)=g_lo%ulim_alloc
         to_low(1)=-ntgrid  ; to_low(2)=1   ; to_low(3)=g_lo%llim_proc       
         to_high(1)=ntgrid ; to_high(2)=2  ; to_high(3)=g_lo%ulim_alloc
         call init_fill (pass_wfb, 'c', to_low, to_high, to, &
               from_low, from_high, from)

         call delete_list (from)
         call delete_list (to)
      endif
  end subroutine init_pass_wfb

end module dist_fn
