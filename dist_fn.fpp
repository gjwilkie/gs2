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
  public :: read_parameters, wnml_dist_fn, check_dist_fn
  public :: timeadv, g_exb, exb_shear
  public :: getan, getfieldeq, getmoms
  public :: flux, eexchange
!  public :: getemoms
!  public :: lf_flux
!  public :: get_heat
  public :: t0, omega0, gamma0
  public :: reset_init, write_f, reset_physics
  public :: M_class, N_class, i_class, par_spectrum
  public :: l_links, r_links, itright, itleft, boundary
  public :: init_kperp2
!  public :: write_fyx
  public :: get_init_field, write_mpdist

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

  private

  real :: apfac, poisfac, driftknob
  real :: t0, omega0, gamma0
  real :: afilter, kfilter
  real :: wfb, g_exb, g_exbfac, omprimfac, btor_slab, mach
  logical :: dfexist, skexist, nonad_zero, lf_default, lf_decompose

  integer :: adiabatic_option_switch
  integer, parameter :: adiabatic_option_default = 1, &
       adiabatic_option_zero = 2, &
       adiabatic_option_fieldlineavg = 3, &
       adiabatic_option_yavg = 4

  integer :: boundary_option_switch
  integer, parameter :: boundary_option_zero = 1, &
       boundary_option_self_periodic = 2, &
       boundary_option_alternate_zero = 3, &
       boundary_option_linked = 4
  logical, public :: def_parity, even
  logical :: test
  logical :: increase = .true., decrease = .false.
  
  ! internal arrays

#ifdef LOWFLOW
  real, dimension (:,:), allocatable :: wcurv
  ! (-ntgrid:ntgrid, -g-layout-)
  real, dimension (:,:,:), allocatable :: wstar_neo
  ! (-ntgrid:ntgrid,ntheta0,naky)
#endif

  real, dimension (:,:,:), allocatable :: wdrift
  ! (-ntgrid:ntgrid, -nvgrid:nvgrid, -g-layout-)

  real, dimension (:,:,:), allocatable :: wstar
  ! (-ntgrid:ntgrid, -nvgrid:nvgrid, -g-layout-)

  ! fieldeq
  real, dimension (:,:,:), allocatable :: gamtot, gamtot1, gamtot2, gamtot3
  ! (-ntgrid:ntgrid,ntheta0,naky) replicated

  !> following arrays needed for implicit solve

  ! coefficients multiplying g_{i+1,j+1}^{n+1}, g_{i+1,j+1}^{n}, g_{i+1,j}^{n+1}, etc.
  real, dimension (:,:,:), allocatable :: pppfac, ppmfac, pmpfac, pmmfac
  real, dimension (:,:,:), allocatable :: mppfac, mpmfac, mmpfac, mmmfac
  ! (-ntgrid:ntgird, -nvgrid:nvgrid, -g-layout-)

  ! response matrix in g
  real, dimension (:,:,:), allocatable :: gresponse
  ! special source term for point where dvpa/dt = 0
  complex, dimension (:,:), allocatable :: source0

  ! matrix needed for reponse matrix approach
  real, dimension (:,:,:), allocatable :: m_mat

  ! these arrays needed to keep track of connections between different
  ! 2pi segments
  integer :: nseg_max
  integer, dimension (:), allocatable :: it_shift, neigen
  integer, dimension (:), allocatable :: ig_low, ig_mid, ig_up, it_shift_left
  integer, dimension (:,:), allocatable :: nsegments, ir_up, iglo_shift

  !< end arrays for implicit solve

  complex, dimension (:,:,:), allocatable :: g0, g_h
  ! (-ntgrid:ntgrid,2, -g-layout-)

  complex, dimension (:,:,:), allocatable :: g_adj
  ! (N(links), 2, -g-layout-)

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
  real, allocatable, dimension(:,:) :: fl_avg, awgt

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
  type (redist_type), save :: pass_left

  integer, dimension (:,:), allocatable :: l_links, r_links
  integer, dimension (:,:,:), allocatable :: n_links
  logical, dimension (:,:), allocatable :: save_h
  logical :: no_comm = .false.
  integer, dimension(:), allocatable :: M_class, N_class
  integer :: i_class

  logical :: initialized = .false.
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

    write (report_unit, *) 
    write (report_unit, fmt="('The standard GK equation will be solved.')")
    write (report_unit, *) 

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
       if (poisfac /= 0.) write (unit, fmt="(' poisfac = ',e16.10)") poisfac
       if (kfilter /= 0.) write (unit, fmt="(' kfilter = ',e16.10)") kfilter
       if (afilter /= 0.) write (unit, fmt="(' afilter = ',e16.10)") afilter
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
       write (unit, fmt="(' source_option = ',a)") '"full"'
       write (unit, fmt="(' /')")
    endif

  end subroutine wnml_dist_fn

  subroutine init_dist_fn

    use mp, only: proc0, finish_mp
    use species, only: init_species, nspec
    use theta_grid, only: init_theta_grid, ntgrid
    use kt_grids, only: init_kt_grids, naky, ntheta0, akx, aky
    use vpamu_grids, only: init_vpamu_grids, nmu, nvgrid
    use run_parameters, only: init_run_parameters
    use gs2_layouts, only: init_dist_fn_layouts, init_gs2_layouts
    use nonlinear_terms, only: init_nonlinear_terms
    use hyper, only: init_hyper

    use gs2_layouts, only: g_lo
    use dist_fn_arrays, only: gnew

    implicit none

    integer :: iglo

    logical:: debug=.false.

    if (initialized) return
    initialized = .true.

    if (debug) write(6,*) "init_dist_fn: init_gs2_layouts"
    call init_gs2_layouts

    if (debug) write(6,*) "init_dist_fn: init_species"
    call init_species

    if (debug) write(6,*) "init_dist_fn: init_theta_grid"
    call init_theta_grid

    if (debug) write(6,*) "init_dist_fn: init_kt_grids"
    call init_kt_grids

    if (debug) write(6,*) "init_dist_fn: init_vpamu_grids"
    call init_vpamu_grids

    if (debug) write(6,*) "init_dist_fn: read_parameters"
    call read_parameters

    if (test) then
       if (proc0) then
          write (*,*) 'nspecies = ',nspec
          write (*,*) 'nvpa = ', 2*nvgrid+1
          write (*,*) 'nmu = ',nmu
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
    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nvgrid, nmu, nspec)

    if (debug) write(6,*) "init_dist_fn: nonlinear_terms"
    call init_nonlinear_terms 

    if (debug) write(6,*) "init_dist_fn: allocate_arrays"
    call allocate_arrays

    if (debug) write(6,*) "init_dist_fn: init_energy"
    call init_vperp2

    if (debug) write(6,*) "init_dist_fn: init_wdrift"
    call init_wdrift

    if (debug) write(6,*) "init_dist_fn: init_wstar"
    call init_wstar

    if (debug) write(6,*) "init_dist_fn: init_bessel"
    call init_bessel

    if (debug) write(6,*) "init_dist_fn: init_par_filter"
    call init_par_filter

#ifdef LOWFLOW
    if (debug) write(6,*) "init_dist_fn: init_lowflow"
    call init_lowflow ! needs to be before init_implicit_solve
#endif

    if (debug) write(6,*) "init_dist_fn: init_implicit_solve"
    call init_implicit_solve

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
            
    namelist /dist_fn_knobs/ boundary_option, nonad_zero, apfac, &
         driftknob, poisfac, adiabatic_option, &
         kfilter, afilter, test, def_parity, even, wfb, &
         g_exb, g_exbfac, omprimfac, btor_slab, mach, lf_default, lf_decompose
    
    namelist /source_knobs/ t0, omega0, gamma0
    integer :: ierr, is, in_file

    if (readinit) return
    readinit = .true.

    if (proc0) then
       boundary_option = 'default'
       nonad_zero = .false.
       adiabatic_option = 'default'
       poisfac = 0.0
       apfac = 1.0
       driftknob = 1.0
       t0 = 100.0
       omega0 = 0.0
       gamma0 = 0.0
       afilter = 0.0
       kfilter = 0.0
       g_exb = 0.0
       g_exbfac = 1.0
       mach = 0.0
       omprimfac = 1.0
       btor_slab = 0.0
       wfb = 1.
       test = .false.
       def_parity = .false.
       even = .true.
       lf_default = .true.
       lf_decompose = .false.
       in_file = input_unit_exist("dist_fn_knobs", dfexist)
       if (dfexist) read (unit=in_file, nml=dist_fn_knobs)

       in_file = input_unit_exist("source_knobs", skexist)
       if (skexist) read (unit=in_file, nml=source_knobs)

       if(abs(shat) <=  1.e-5) boundary_option = 'periodic'

       ierr = error_unit()
       call get_option_value &
            (boundary_option, boundaryopts, boundary_option_switch, &
            ierr, "boundary_option in dist_fn_knobs")

       call get_option_value &
            (adiabatic_option, adiabaticopts, adiabatic_option_switch, &
            ierr, "adiabatic_option in dist_fn_knobs")

    end if

    call broadcast (boundary_option_switch)
    call broadcast (nonad_zero)
    call broadcast (adiabatic_option_switch)
    call broadcast (poisfac)
    call broadcast (apfac)
    call broadcast (driftknob)
    call broadcast (t0)
    call broadcast (omega0)
    call broadcast (gamma0)
    call broadcast (g_exb)
    call broadcast (g_exbfac)
    call broadcast (mach)
    call broadcast (omprimfac)
    call broadcast (btor_slab)
    call broadcast (afilter)
    call broadcast (kfilter)
    call broadcast (test)
    call broadcast (def_parity)
    call broadcast (lf_default)
    call broadcast (lf_decompose)
    call broadcast (even)
    call broadcast (wfb)

!! Override itor_over_B, if "dist_fn_knobs" parameter btor_slab ne 0
       if (abs(btor_slab) > epsilon(0.0)) itor_over_B = btor_slab
!! Done for slab, where itor_over_B is determined by angle between B-field 
!! and toroidal flow: itor_over_B = (d(u_z)/dx) / (d(u_y)/dx) = Btor / Bpol
!! u = u0 (phihat) = x d(u0)/dx (phihat) = x d(uy)/dx (yhat + Btor/Bpol zhat)
!! g_exb = d(uy)/dx => d(uz)/dx = g_exb * Btor/Bpol = g_exb * itor_over_B

  end subroutine read_parameters 

  subroutine init_wdrift

    use centering, only: get_cell_value
    use species, only: nspec
    use theta_grid, only: ntgrid, bmag, thet_imp
    use kt_grids, only: naky, ntheta0
    use vpamu_grids, only: nvgrid, vpa_imp
    use gs2_layouts, only: g_lo, ik_idx, it_idx, imu_idx, is_idx

    implicit none

    integer :: iglo, ig, ik, it, iv, imu, is
    logical :: debug = .false.

    if (.not. allocated(wdrift)) then
       ! allocate wdrift with sign(vpa) dependence because will contain 
       ! Coriolis as well as magnetic drifts
       allocate (wdrift(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    end if
    wdrift = 0.
#ifdef LOWFLOW
    if (.not. allocated(wcurv)) allocate (wcurv(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    wcurv = 0.
#endif

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do iv = -nvgrid, nvgrid
          do ig = -ntgrid, ntgrid
             imu=imu_idx(g_lo,iglo)
             it=it_idx(g_lo,iglo)
             ik=ik_idx(g_lo,iglo)
             is=is_idx(g_lo,iglo)
             ! get grad-B and curvature drifts
             wdrift(ig,iv,iglo) = wdrift_func(ig,iv,imu,it,ik)*driftknob
#ifdef LOWFLOW
             ! get curvature drift without vpa dependence
             wcurv(ig,iglo) = wcurv_func(ig, it, ik)*driftknob
#endif
             ! add Coriolis drift to magnetic drifts
             wdrift(ig,iv,iglo) = wdrift(ig,iv,iglo) + wcoriolis_func(ig,iv,it,ik,is)
          end do
       end do
    end do

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ! get wdrift at center of theta-vpa cells
       call get_cell_value (thet_imp, vpa_imp, &
            wdrift(:,:,iglo), wdrift(:,:,iglo), -ntgrid, -nvgrid)
#ifdef LOWFLOW
       ! get wcurv at center of theta cell
       call get_cell_value (thet_imp, wcurv(:,iglo), wcurv(:,iglo), -ntgrid)
#endif
    end do

    if (debug) write(6,*) 'init_wdrift: driftknob',driftknob

  end subroutine init_wdrift

  function wdrift_func (ig, iv, imu, it, ik)

    use theta_grid, only: bmag, gbdrift, gbdrift0, cvdrift, cvdrift0
    use theta_grid, only: shat
    use kt_grids, only: aky, theta0, akx
    use vpamu_grids, only: vpa, vperp2
    use run_parameters, only: wunits
    use gs2_time, only: code_dt

    implicit none

    real :: wdrift_func
    integer, intent (in) :: ig, ik, it, iv, imu

    ! note that wunits=aky/2 (for wstar_units=F)
    if (aky(ik) == 0.0) then
       wdrift_func = akx(it)/shat &
            * (cvdrift0(ig)*vpa(iv)**2 + gbdrift0(ig)*0.5*vperp2(ig,imu)) &
            * code_dt/2.0
    else
       wdrift_func = ((cvdrift(ig) + theta0(it,ik)*cvdrift0(ig)) * vpa(iv)**2 &
            + (gbdrift(ig) + theta0(it,ik)*gbdrift0(ig)) * 0.5*vperp2(ig,imu)) &
            * code_dt*wunits(ik)
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

  function wcoriolis_func (ig, iv, it, ik, is)

    use theta_grid, only: bmag, cdrift, cdrift0, shat
    use kt_grids, only: aky, theta0, akx
    use vpamu_grids, only: vpa
    use run_parameters, only: wunits
    use gs2_time, only: code_dt
    use species, only: spec

    implicit none

    real :: wcoriolis_func
    integer, intent (in) :: ig, ik, it, iv, is

    if (aky(ik) == 0.0) then
       wcoriolis_func = mach * vpa(iv) &
            * cdrift0(ig) * code_dt * akx(it)/(2.*shat*spec(is)%stm)
    else
       wcoriolis_func = mach * vpa(iv) &
            * (cdrift(ig) + theta0(it,ik)*cdrift0(ig))*code_dt*wunits(ik)/spec(is)%stm
    end if

  end function wcoriolis_func

  subroutine init_vperp2

    use theta_grid, only: ntgrid, bmag
    use vpamu_grids, only: vperp2, nmu, mu, energy, anon, nvgrid, vpa

    implicit none

    integer :: iv
    
    if (.not.allocated(vperp2)) then
       allocate (vperp2 (-ntgrid:ntgrid,nmu))
       allocate (energy (-ntgrid:ntgrid,-nvgrid:nvgrid,nmu))
       allocate (anon (-ntgrid:ntgrid,-nvgrid:nvgrid,nmu))
!       allocate (vpar   (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    endif
    vperp2 = 0. ; energy = 0. ; anon = 0.

    vperp2 = 2.0*spread(mu,1,2*ntgrid+1)*spread(bmag,2,nmu)

    do iv = -nvgrid, nvgrid
       energy(:,iv,:) = vpa(iv)**2 + 2.0*spread(mu,1,2*ntgrid+1)*spread(bmag,2,nmu)
       anon(:,iv,:) = exp(-energy(:,iv,:))
    end do

!    do iglo = g_lo%llim_proc, g_lo%ulim_proc
!CMR : vpar = q_s/sqrt{T_s m_s}*DELT/DTHETA * vpac |\gradpar(theta)| 
!                                     where gradpar(theta) is centered
!  ie  vpar = q_s/sqrt{T_s m_s} (v_||^GS2). \gradpar(theta)/DTHETA . DELT
!CMR 

!        ik = ik_idx(g_lo,iglo)
!        is = is_idx(g_lo,iglo)
!        vpar(-ntgrid:ntgrid-1,1,iglo) = &
!             spec(is)%zstm*tunits(ik)*code_dt &
!             *0.5/delthet(-ntgrid:ntgrid-1) &
!             *(abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid)))&
!             *vpac(-ntgrid:ntgrid-1,1,iglo)
!        vpar(-ntgrid:ntgrid-1,2,iglo) = &
!             spec(is)%zstm*tunits(ik)*code_dt &
!             *0.5/delthet(-ntgrid:ntgrid-1) &
!             *(abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid)))&
!             *vpac(-ntgrid:ntgrid-1,2,iglo)
!        vpar(ntgrid,:,iglo) = 0.0
!    end do

  end subroutine init_vperp2

  subroutine init_wstar

    use species, only: spec
    use theta_grid, only: ntgrid
    use vpamu_grids, only: nvgrid, vpa, vperp2
    use run_parameters, only: wunits
    use gs2_time, only: code_dt
    use gs2_layouts, only: is_idx, ik_idx, imu_idx, g_lo

    implicit none

    integer :: ik, is, imu, ig, iv, iglo

    if (.not. allocated(wstar)) then
       allocate (wstar(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
       wstar = 0.0
    end if

    do iglo = g_lo%llim_proc, g_lo%ulim_alloc
       is = is_idx(g_lo,iglo)
       imu = imu_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       do iv = -nvgrid, nvgrid
          do ig = -ntgrid, ntgrid
             wstar(ig,iv,iglo) = code_dt*wunits(ik) &
                  *(spec(is)%fprim+spec(is)%tprim*(vpa(iv)**2+vperp2(ig,imu)-1.5))
          end do
       end do
    end do

  end subroutine init_wstar

  subroutine init_bessel

    use dist_fn_arrays, only: aj0, aj1, aj2, kperp2
    use species, only: spec
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, aky, akx
    use vpamu_grids, only: vperp2
    use gs2_layouts, only: g_lo, ik_idx, it_idx, imu_idx, is_idx
    use spfunc, only: j0, j1

    implicit none

    integer :: ig, ik, it, imu, is
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
       imu = imu_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do ig = -ntgrid, ntgrid
          arg = spec(is)%smz*sqrt(vperp2(ig,imu)*kperp2(ig,it,ik))/bmag(ig)
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

  subroutine init_implicit_solve

    use gs2_layouts, only: g_lo, imu_idx
    use gs2_time, only: code_dt
    use run_parameters, only: t_imp
    use theta_grid, only: ntgrid, dbdthetc, thet_imp, delthet, ntheta
    use vpamu_grids, only: nvgrid, vpa_imp, vpac, dvpa, vpa, mu
    use kt_grids, only: naky, ntheta0

    implicit none

    integer :: iglo, imu, ig, iv, ntg
    real :: thm_fac, vpm_fac, decay_fac, vpp_fac, thp_fac
    complex :: wd1, wd2

    real, dimension (:,:), allocatable :: dum1, dum2

    call init_connections

    ntg = ntheta/2

    if (.not. allocated(gresponse)) then
       allocate (gresponse(ntg*nseg_max+1,ntg*nseg_max+1,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (m_mat(ntg*nseg_max+1,ntg*nseg_max+1,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (source0(nseg_max,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (pppfac(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; pppfac = 0.0
       allocate (ppmfac(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; ppmfac = 0.0
       allocate (pmpfac(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; pmpfac = 0.0
       allocate (pmmfac(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; pmmfac = 0.0
       allocate (mppfac(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; mppfac = 0.0
       allocate (mpmfac(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; mpmfac = 0.0
       allocate (mmpfac(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; mmpfac = 0.0
       allocate (mmmfac(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; mmmfac = 0.0
    end if

    if (.not. allocated (dum1)) then
       allocate (dum1(-ntgrid:ntgrid-1,-nvgrid:nvgrid-1)) ; dum1 = 0.0
       allocate (dum2(-ntgrid:ntgrid-1,-nvgrid:nvgrid-1)) ; dum2 = 0.0
    end if

    ! need to add in multiplication by b . grad theta
    dum1(-ntgrid:-1,:) = code_dt*spread(vpac(:nvgrid-1,1),1,ntgrid)/spread(delthet(-ntgrid:-1),2,2*nvgrid)
    dum1(0:ntgrid-1,:) = code_dt*spread(vpac(:nvgrid-1,2),1,ntgrid)/spread(delthet(0:ntgrid-1),2,2*nvgrid)

    do iv = -nvgrid, -1
       where (dbdthetc(:ntgrid-1,1) < 0.0)
          dum1(:ntgrid-1,iv) = code_dt*vpac(iv,1)/delthet(:ntgrid-1)
       elsewhere
          dum1(:ntgrid-1,iv) = code_dt*vpac(iv,2)/delthet(:ntgrid-1)
       end where
    end do
    do iv = 0, nvgrid-1
       where (dbdthetc(:ntgrid-1,2) < 0.0)
          dum1(:ntgrid-1,iv) = code_dt*vpac(iv,1)/delthet(:ntgrid-1)
       elsewhere
          dum1(:ntgrid-1,iv) = code_dt*vpac(iv,2)/delthet(:ntgrid-1)
       end where
    end do

    ! TMP FOR TESTING -- MAB
!     do iv = -nvgrid, nvgrid-1
!        do ig = -ntgrid, ntgrid-1
!           write (*,'(a5,5e12.4)') 'dum1', dum1(ig,iv), code_dt, delthet(ig), vpac(iv,1), vpac(iv,2)
!        end do
!     end do

    decay_fac = exp(dvpa(-nvgrid)*(dvpa(-nvgrid)-2.*abs(vpa(-nvgrid))))

    do iglo = g_lo%llim_proc, g_lo%ulim_proc

       ! get im corresponding to iglo
       imu = imu_idx(g_lo,iglo)

       ! MAB FLAG - need to add in b . grad theta !!!
       dum2(:,-nvgrid:-1) = -code_dt*mu(imu) &
            / spread(dvpa(-nvgrid:-1),1,2*ntgrid)*spread(dbdthetc(:ntgrid-1,1),2,nvgrid)
       dum2(:,0:nvgrid-1) = -code_dt*mu(imu) &
            / spread(dvpa(0:nvgrid-1),1,2*ntgrid)*spread(dbdthetc(:ntgrid-1,2),2,nvgrid)

       ! first treat theta<0, vpa>0 quadrant
       ! to upwind here requires using info from ig-1,iv-1
       ! i.e., ig+1 gets (1-thet_imp) and ig gets thet_imp,
       ! iv+1 gets (1-vpa_imp) and iv gets vpa_imp

       do iv = -nvgrid, nvgrid-1
          do ig = -ntgrid, ntgrid-1

             if (dum1(ig,iv) > 0.0) then
                thp_fac = thet_imp
             else
                thp_fac = 1.0-thet_imp
             end if
             if (dum2(ig,iv) > 0.0) then
                vpp_fac = vpa_imp
             else
                vpp_fac = 1.0-vpa_imp
             end if
             thm_fac = 1.0-thp_fac ; vpm_fac = 1.0-vpp_fac

             ! wd1 and wd2 account for contributions from wdrift
!             wd1 = 1.0 + t_imp*zi*wdrift(ig,iv,iglo)
!             wd2 = 1.0 + (t_imp-1.0)*zi*wdrift(ig,iv,iglo)
             ! TMP FOR TESTING -- MAB
             wd1 = 1.0 ; wd2 = 1.0

             pppfac(ig,iv,iglo) = thp_fac*vpp_fac*wd1 &
                  + t_imp*vpp_fac*dum1(ig,iv) &
                  + t_imp*thp_fac*dum2(ig,iv)
             ppmfac(ig,iv,iglo) = thp_fac*vpm_fac*wd1 &
                  + t_imp*vpm_fac*dum1(ig,iv) &
                  - t_imp*thp_fac*dum2(ig,iv)
             pmpfac(ig,iv,iglo) = thm_fac*vpp_fac*wd1 &
                  - t_imp*vpp_fac*dum1(ig,iv) &
                  + t_imp*thm_fac*dum2(ig,iv)
             pmmfac(ig,iv,iglo) = thm_fac*vpm_fac*wd1 &
                  - t_imp*vpm_fac*dum1(ig,iv) &
                  - t_imp*thm_fac*dum2(ig,iv)
             mppfac(ig,iv,iglo) = -thp_fac*vpp_fac*wd2 &
                  + (1.0-t_imp)*vpp_fac*dum1(ig,iv) &
                  + (1.0-t_imp)*thp_fac*dum2(ig,iv)
             mpmfac(ig,iv,iglo) = -thp_fac*vpm_fac*wd2 &
                  + (1.0-t_imp)*vpm_fac*dum1(ig,iv) &
                  - (1.0-t_imp)*thp_fac*dum2(ig,iv)
             mmpfac(ig,iv,iglo) = -thm_fac*vpp_fac*wd2 &
                  - (1.0-t_imp)*vpp_fac*dum1(ig,iv) &
                  + (1.0-t_imp)*thm_fac*dum2(ig,iv)
             mmmfac(ig,iv,iglo) = -thm_fac*vpm_fac*wd2 &
                  - (1.0-t_imp)*vpm_fac*dum1(ig,iv) &
                  - (1.0-t_imp)*thm_fac*dum2(ig,iv)

! TMP FOR TESTING -- MAB
!              write (*,*) '!-------------------- xxxfac -----------------------!'
!              write (*,'(a7,3i4)') 'xxxfac', ig, iv, iglo
!              write (*,'(a7,4e12.4)') 'xxxfac', pppfac(ig,iv,iglo), ppmfac(ig,iv,iglo), pmpfac(ig,iv,iglo), pmmfac(ig,iv,iglo)
!              write (*,'(a7,4e12.4)') 'xxxfac', mppfac(ig,iv,iglo), mpmfac(ig,iv,iglo), mmpfac(ig,iv,iglo), mmmfac(ig,iv,iglo)
!              write (*,'(a7,5e12.4)') 'xxxfac', thm_fac, vpm_fac, t_imp, dum1(ig,iv), dum2(ig,iv)
!              write (*,'(a7,4e12.4)') 'xxxfac', mu(imu), dvpa(iv), dbdthetc(ig,1), dbdthetc(ig,2)
!              write (*,*) '!---------------------------------------------------!'

          end do
       end do

       ! treat theta<0, vpa=-vpa_max and theta>0, vpa=vpa_max specially
       ! in particular, assume g decays like a Maxwellian in vpa at this point
       where (dbdthetc(:ntgrid-1,1) < 0.0)
          pppfac(:ntgrid-1,-nvgrid,iglo) = pppfac(:ntgrid-1,-nvgrid,iglo) &
               + ppmfac(:ntgrid-1,-nvgrid,iglo)*decay_fac
          pmpfac(:ntgrid-1,-nvgrid,iglo) = pmpfac(:ntgrid-1,-nvgrid,iglo) &
               + pmmfac(:ntgrid-1,-nvgrid,iglo)*decay_fac
          mppfac(:ntgrid-1,-nvgrid,iglo) = mppfac(:ntgrid-1,-nvgrid,iglo) &
               + mpmfac(:ntgrid-1,-nvgrid,iglo)*decay_fac
          mmpfac(:ntgrid-1,-nvgrid,iglo) = mmpfac(:ntgrid-1,-nvgrid,iglo) &
               + mmmfac(:ntgrid-1,-nvgrid,iglo)*decay_fac
          ppmfac(:ntgrid-1,-nvgrid,iglo) = 0.
          pmmfac(:ntgrid-1,-nvgrid,iglo) = 0.
          mpmfac(:ntgrid-1,-nvgrid,iglo) = 0.
          mmmfac(:ntgrid-1,-nvgrid,iglo) = 0.
       elsewhere
          ppmfac(:ntgrid-1,nvgrid-1,iglo) = ppmfac(:ntgrid-1,nvgrid-1,iglo) &
               + pppfac(:ntgrid-1,nvgrid-1,iglo)*decay_fac
          pmmfac(:ntgrid-1,nvgrid-1,iglo) = pmmfac(:ntgrid-1,nvgrid-1,iglo) &
               + pmpfac(:ntgrid-1,nvgrid-1,iglo)*decay_fac
          mpmfac(:ntgrid-1,nvgrid-1,iglo) = mpmfac(:ntgrid-1,nvgrid-1,iglo) &
               + mppfac(:ntgrid-1,nvgrid-1,iglo)*decay_fac
          mmmfac(:ntgrid-1,nvgrid-1,iglo) = mmmfac(:ntgrid-1,nvgrid-1,iglo) &
               + mmpfac(:ntgrid-1,nvgrid-1,iglo)*decay_fac
          pppfac(:ntgrid-1,nvgrid-1,iglo) = 0.
          pmpfac(:ntgrid-1,nvgrid-1,iglo) = 0.
          mppfac(:ntgrid-1,nvgrid-1,iglo) = 0.
          mmpfac(:ntgrid-1,nvgrid-1,iglo) = 0.
       end where

!        do iv = -nvgrid, nvgrid
!           do ig = -ntgrid, ntgrid
!              write (*,*) '!-------------------- xxxfac -----------------------!'
!              write (*,'(a7,4e12.4)') 'xxxfac', pppfac(ig,iv,iglo), ppmfac(ig,iv,iglo), pmpfac(ig,iv,iglo), pmmfac(ig,iv,iglo)
!              write (*,'(a7,4e12.4)') 'xxxfac', mppfac(ig,iv,iglo), mpmfac(ig,iv,iglo), mmpfac(ig,iv,iglo), mmmfac(ig,iv,iglo)
!              write (*,*) '!---------------------------------------------------!'
!           end do
!           write (*,*)
!        end do

    end do

    ! get response of g at theta<=0, vpa=-dvpa
    ! to unit impulses in g at theta<0, vpa=0
    call get_gresponse_matrix

! TMP FOR TESTING -- MAB
!     do iglo = g_lo%llim_proc, g_lo%ulim_proc
!        do iv = 1, ntg*nseg_max+1
!           do ig = 1, ntg*nseg_max+1
!              write (*,'(a10,2e12.4)'), 'gresponse', gresponse(ig,iv,iglo), m_mat(ig,iv,iglo)
!           end do
!        end do
!     end do

    ! MAB FLAG -- need to figure out how to deal with i_class, M_class, N_class

    i_class = 1
    if (.not. allocated(M_class)) then
       allocate (M_class(i_class))
       allocate (N_class(i_class))
    end if
    M_class = naky*ntheta0 ; N_class = 1

!     select case (boundary_option_switch)
!     case (boundary_option_linked)
!        call init_connected_bc
!     case default
!        if (.not. allocated(l_links)) then
!           allocate (l_links(naky, ntheta0))
!           allocate (r_links(naky, ntheta0))
!           allocate (n_links(2, naky, ntheta0))
!        end if
!        l_links = 0;   r_links = 0;  n_links = 0
       
!        i_class = 1
!        if (.not. allocated(M_class)) then
!           allocate (M_class(i_class))
!           allocate (N_class(i_class))
!        end if
!        M_class = naky*ntheta0 ; N_class = 1
       
!     end select

    if (allocated(dum1)) deallocate (dum1, dum2)

  end subroutine init_implicit_solve

  subroutine init_connections

    use mp, only: nproc, mp_abort
    use theta_grid, only: nperiod, ntgrid, ntheta
    use kt_grids, only: ntheta0, jtwist_out, naky
    use species, only: nspec
    use vpamu_grids, only: nmu

    implicit none

    integer :: iseg, ik, ie, ntg, it, neigen_max

    ntg = ntheta/2

    if (.not. allocated(it_shift)) then
       allocate (it_shift(naky))
       allocate (neigen(naky))
    end if

    select case (boundary_option_switch)
    case (boundary_option_zero)
       
       it_shift = ntheta0
       neigen = 1 ; neigen_max = 1
       
       if (.not. allocated(it_shift_left)) then
          allocate (it_shift_left(ntheta0))
          allocate (iglo_shift(ntheta0,naky))
       end if
       iglo_shift = 0 ; it_shift_left = 0
       
       if (.not. allocated(nsegments)) then
          allocate (nsegments(neigen_max,naky))
          allocate (ir_up(neigen_max,naky))
       end if
       
       ! this is the number of 2pi poloidal segments in the extended theta domain,
       ! which is needed in initializing the reponse matrix and doing the implicit sweep
       nsegments = 2*(nperiod-1) + 1
       
       nseg_max = maxval(nsegments)
       
       ir_up = ntg*nseg_max+1
       
       if (.not. allocated(ig_low)) then
          allocate (ig_low(nseg_max))
          allocate (ig_mid(nseg_max))
          allocate (ig_up(nseg_max))
       end if

       ! ig_low(j) is the ig index corresponding to the inboard midplane from below (theta=-pi) within the jth segment
       ! ig_mid(j) is the ig index corresponding to the outboard midplane (theta=0) within the jth segment
       do iseg = 1, nseg_max
          ig_low(iseg) = -ntgrid + (iseg-1)*ntheta
          ig_mid(iseg) = ig_low(iseg) + ntheta/2
          ig_up(iseg) = ig_low(iseg) + ntheta
       end do

    case (boundary_option_linked)

       if (nproc > naky*nmu*nspec) then
          write (*,*) 'Parallelization over kx not currently supported'
          write (*,*) 'with twist and shift boundary condition.  Aborting.'
          call mp_abort ('Parallelization over kx not currently supported with twist-and-shift BC.')
       end if
       
       ik = 1
       it_shift(ik) = ntheta0 ; neigen(ik) = it_shift(ik)
       if (naky > 1) then
          do ik = 2, naky
             ! must link different kx's at theta = +/- pi
             it_shift(ik) = (ik-1)*jtwist_out
             ! neigen is the number of independent eigenfunctions along the field line
             neigen(ik) = it_shift(ik)
          end do
       end if

       neigen_max = maxval(neigen)

       if (.not. allocated(iglo_shift)) then
          allocate (it_shift_left(ntheta0))
          allocate (iglo_shift(ntheta0,naky)) ; iglo_shift = 0
       end if

       ! figure out how much to shift iglo by to get to
       ! the left-most (theta-theta0) in each set of connected 2pi segments
       do it = 1, ntheta0
          if (it <= ntheta0/2+1) then
             it_shift_left(it) = ntheta0/2-2*it+2
          else
             it_shift_left(it) = ntheta0/2-2*(it-ntheta0/2-1)+1
          end if
       end do

       do ik = 1, naky
          ! iglo_shift is how much to shift each iglo by to connect
          ! to the next theta0 (from most positive to most negative)
          do it = 1, ntheta0
             if (it > ntheta0/2+1) then
                if (it-it_shift(ik) >= ntheta0/2+1) then
                   iglo_shift(it,ik) = -it_shift(ik)
                end if
             else
                if (it-it_shift(ik) > 0) then
                   iglo_shift(it,ik) = -it_shift(ik)
                else if (it-it_shift(ik)+ntheta0 >= 2+ntheta0/2) then
                   iglo_shift(it,ik) = ntheta0 - it_shift(ik)
                end if
             end if
          end do
       end do

       if (.not. allocated(nsegments)) then
          allocate (nsegments(neigen_max,naky))
          allocate (ir_up(neigen_max,naky))
       end if

       do ik = 1, naky
          if (it_shift(ik) == 0) then
             nsegments(:,ik) = 1
          else
             nsegments(:,ik) = (ntheta0-1)/it_shift(ik)

             do ie = 1, mod(ntheta0-1,it_shift(ik))+1
                nsegments(ie,ik) = nsegments(ie,ik) + 1
             end do
          end if
       end do

       ir_up = ntg*nsegments+1

       nseg_max = maxval(nsegments)

       if (.not. allocated(ig_low)) then
          allocate (ig_low(nseg_max)) ; ig_low = -ntgrid
          allocate (ig_mid(nseg_max)) ; ig_mid = 0
          allocate (ig_up(nseg_max)) ; ig_up = ntgrid
       end if

    end select

  end subroutine init_connections

  subroutine get_gresponse_matrix

    use theta_grid, only: ntgrid, ntheta
    use vpamu_grids, only: nvgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx
    use centering, only: invert_matrix
    use dist_fn_arrays, only: gnew, source

    implicit none

    integer :: ig, iv, iglo, i, j, k, idx, ntg, iseg, it, itmod, ik, kmax
    integer, dimension (:), allocatable :: iglomod
    real, dimension (:,:), allocatable :: p_mat

    ntg = ntheta/2
    
    ! ensure that gnew and source are initialized to zero
    gnew = 0.0 ; source = 0.0 ; gresponse = 0.0 ; source0 = 0.0

    m_mat = 0.0
    allocate (p_mat(ntg*nseg_max+1,ntg*nseg_max+1)) ; p_mat = 0.0

    ! get response of particles with vpa=-dvpa below and including the midplane to a unit impulse
    ! in particles with vpa=0 below (not including) the midplane
    ! this will be used to calculate the distribution of particles with vpa=0 below the midplane, which
    ! starts the implicit sweep in theta and vpa;
    do iglo = g_lo%llim_proc, g_lo%ulim_proc

       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)

       ! only do the following once for each independent set of theta0s
       ! the assumption here is that all kx are on processor and sequential
       if (it > it_shift(ik)) cycle

       allocate (iglomod(nsegments(it,ik)))

       ! remap to start at theta0 = -theta0_max for this set of connected theta0s
       iseg = 1
       iglomod(iseg) = iglo + it_shift_left(it)
       itmod = it + it_shift_left(it)
       if (nsegments(it,ik) > 1) then
          do iseg = 2, nsegments(it,ik)
             iglomod(iseg) = iglomod(iseg-1) + iglo_shift(itmod,ik)
             itmod = itmod + iglo_shift(itmod,ik)
          end do
       end if

       i = 1

       ! for iseg=1 (the left-most segment), g(theta=-pi,vpa=0) is not known
       ! for subsequent segments, g(theta=-pi,vpa=0) will have been solved for
       ! thus iseg=1 must be treated specially
       iseg = 1
       ! give a unit impulse to each theta from -pi up to but not including theta=0 (for vpa=0)
       ! and obtain the corresponding columns in the response matrix
       do ig = ig_low(iseg), ig_mid(iseg)-1
          ! give a unit impulse to g at this theta
          gnew(ig,0,iglomod(iseg)) = 1.0

          ! sweep over entire (theta,vpa) plane and get values of gnew 
          ! for particles with vpa=-dvpa below the midplane
          call implicit_sweep (iglo)

          ! fill in cotribution from jth segment to ith column of initial 
          ! (ntheta/2*nsegment+1) x (ntheta/2*nsegment+1) x (nmu) response matrix          
          ! again, must treat left-most segment specially
          j = 1 ; k = 1 ; kmax = k+ntheta/2
          gresponse(k:kmax,i,iglo) = real(gnew(ig_low(iseg):ig_mid(iseg),-1,iglomod(iseg)))          
          k = kmax + 1
          if (nsegments(it,ik) > 1) then
             do j = 2, nsegments(it,ik)
                kmax = k+ntheta/2-1
                gresponse(k:kmax,i,iglo) = real(gnew(ig_low(iseg)+1:ig_mid(iseg),-1,iglomod(iseg)))
                k = kmax + 1
             end do
          end if

          ! reset g to zero everywhere
          gnew(:,:,iglomod(iseg)) = 0.0

          i = i + 1
       end do

       ! put in blank column for theta=0 since it can be determined without
       ! coupling to other theta-vpa points
       i = i + 1

       if (nsegments(it,ik) > 1) then
          do iseg = 2, nsegments(it,ik)
             do ig = ig_low(iseg)+1, ig_mid(iseg)-1
                ! give a unit impulse to g at this theta
                gnew(ig,0,iglomod(iseg)) = 1.0

                ! sweep over (theta,vpa) plane and get values of gnew for theta<0, vpa=-dvpa
                call implicit_sweep (iglo)
                
                ! fill in cotribution from jth segment to ith column of initial 
                ! (ntheta/2*nsegment+1) x (ntheta/2*nsegment+1) x (nmu) response matrix          
                ! again, must treat left-most segment specially
                j = 1 ; k = 1 ; kmax = k+ntheta/2
                gresponse(k:kmax,i,iglo) = real(gnew(ig_low(iseg):ig_mid(iseg),-1,iglomod(iseg)))
                k = kmax + 1
                if (nsegments(it,ik) > 1) then
                   do j = 2, nsegments(it,ik)
                      kmax = k+ntheta/2-1
                      gresponse(k:kmax,i,iglo) = real(gnew(ig_low(iseg)+1:ig_mid(iseg),-1,iglomod(iseg)))
                      k = kmax + 1
                   end do
                end if

                ! reset g to zero everywhere
                gnew(:,:,iglomod(iseg)) = 0.0

                i = i + 1
             end do
             ! put in blank column for theta=0 since it can be determined without
             ! coupling to other theta-vpa points
             i = i + 1
          end do
       end if

       ! get the matrices multiplying g(theta<=0,vpa=0) and g(theta<=0,dvpa=-dvpa) 
       ! in the implicit solve. these are necessary to obtain the final response matrix

       ! p_mat is matrix multiplying g_{0} in equations for g(theta<0,vpa=0) (M1 from notes)
       ! m_mat is matrix multiplying g_{-1} in same equations (M2 from notes)

       idx = 1

       iseg = 1
       do ig = ig_low(iseg), ig_mid(iseg)-1
          p_mat(idx,idx) = pmpfac(ig,-1,iglomod(iseg))
          p_mat(idx,idx+1) = pppfac(ig,-1,iglomod(iseg))
          m_mat(idx,idx,iglo) = pmmfac(ig,-1,iglomod(iseg))
          m_mat(idx,idx+1,iglo) = ppmfac(ig,-1,iglomod(iseg))
          idx = idx + 1
       end do
       ! account for totally trapped particles, which decouple from the other theta-vpa points
       p_mat(idx,idx) = 1.0
       idx = idx + 1

       if (nsegments(it,ik) > 1) then
          do iseg = 2, nsegments(it,ik)
             do ig = ig_low(iseg)+1, ig_mid(iseg)-1
                p_mat(idx,idx) = pmpfac(ig,-1,iglomod(iseg))
                p_mat(idx,idx+1) = pppfac(ig,-1,iglomod(iseg))
                m_mat(idx,idx,iglo) = pmmfac(ig,-1,iglomod(iseg))
                m_mat(idx,idx+1,iglo) = ppmfac(ig,-1,iglomod(iseg))
                idx = idx + 1
             end do
             ! account for totally trapped particles, which decouple from the other theta-vpa points
             p_mat(idx,idx) = 1.0
             idx = idx + 1
          end do
       end if

       gresponse(:ir_up(it,ik),:ir_up(it,ik),iglo) = p_mat(:ir_up(it,ik),:ir_up(it,ik)) &
            + matmul(m_mat(:ir_up(it,ik),:ir_up(it,ik),iglo),gresponse(:ir_up(it,ik),:ir_up(it,ik),iglo))

       ! take the inverse to get the final response matrix
       call invert_matrix (gresponse(:ir_up(it,ik),:ir_up(it,ik),iglo))

       p_mat = 0.0

       deallocate (iglomod)
    end do

    deallocate (p_mat)

  end subroutine get_gresponse_matrix

  ! implicit_sweep starts with a boundary condition for g^{n+1} along the line (theta<0,vpa=0)
  ! as well as zero BCs at (theta=-theta_max,vpa>0), (theta=theta_max,vpa<0), 
  ! (vpa=-vpa_max,theta<0), and (vpa=vpa_max,theta>0), and solves for g^{n+1}.  
  ! note that dvpa/dt<0 for theta>0 and dvpa/dt>0 for theta<0.
  subroutine implicit_sweep (iglo, source_mod)

    use theta_grid, only: ntgrid, ntheta
    use vpamu_grids, only: nvgrid
    use gs2_layouts, only: g_lo, it_idx, ik_idx
    use dist_fn_arrays, only: gnew, source

    implicit none

    integer, intent (in) :: iglo
    real, dimension (:), intent (in out), optional :: source_mod

    integer :: ig, iv, ntg, k, iseg, kmax, it, itmod, ik
    integer, dimension (:), allocatable :: iglomod
    complex, dimension (:), allocatable :: dgdgn

    it = it_idx(g_lo,iglo)
    ik = ik_idx(g_lo,iglo)

    allocate (iglomod(nsegments(it,ik)))

    iseg = 1
    iglomod(iseg) = iglo + it_shift_left(it)
    itmod = it + it_shift_left(it)
    if (nsegments(it,ik) > 1) then
       do iseg = 2, nsegments(it,ik)
          iglomod(iseg) = iglomod(iseg-1) + iglo_shift(itmod,ik)
          itmod = itmod + iglo_shift(itmod,ik)
       end do
    end if

    ntg = ntheta/2

    if (present(source_mod)) then
       allocate (dgdgn(ntg*nseg_max+1)) ; dgdgn = 0.
    else
       allocate (dgdgn(1)) ; dgdgn = 0.
    end if

    iseg = 1
    ! initialize gnew to zero for all (theta,vpa) in this segment except (theta<0,vpa=0)
    ! which was set earlier as initial condition
    gnew(ig_low(iseg):ig_up(iseg),-nvgrid:-1,iglomod(iseg)) = 0.0
    gnew(ig_low(iseg):ig_up(iseg),1:,iglomod(iseg)) = 0.0
    gnew(ig_mid(iseg):ig_up(iseg),0,iglomod(iseg)) = 0.0

    call implicit_sweep_right (iglomod(iseg), iseg)
    if (nsegments(it,ik) > 1) then
       do iseg = 2, nsegments(it,ik)
          ! initialize gnew to zero for all (theta,vpa) except (theta<0,vpa=0)
          ! which was set earlier as initial condition
          gnew(ig_low(iseg):ig_up(iseg),-nvgrid:-1,iglomod(iseg)) = 0.0
          ! note that gnew(ig_low(iseg),1:) will be set via connection to the previous segment
          gnew(ig_low(iseg)+1:ig_up(iseg),1:,iglomod(iseg)) = 0.0
          gnew(ig_mid(iseg):ig_up(iseg),0,iglomod(iseg)) = 0.0

          ! make connection with previous segment by setting gnew(-pi) for this segment
          ! equal to gnew(pi) from previous segment
          gnew(ig_low(iseg),0:,iglomod(iseg)) = gnew(ig_up(iseg-1),0:,iglomod(iseg-1))

          call implicit_sweep_right (iglomod(iseg), iseg)
       end do
    end if

    iseg = nsegments(it,ik)
    call implicit_sweep_left(iglomod(iseg), iseg)
    if (nsegments(it,ik) > 1) then
       do iseg = nsegments(it,ik)-1, 1, -1
          
          ! make connection with previous segment by setting gnew(pi) for this segment
          ! equal to gnew(-pi) from previous segment
          gnew(ig_up(iseg),-nvgrid:-1,iglomod(iseg)) &
               = gnew(ig_low(iseg+1),-nvgrid:-1,iglomod(iseg+1))

          call implicit_sweep_left (iglomod(iseg), iseg)
       end do
    end if

    if (present(source_mod)) then

       k = 1

       iseg = 1 ; kmax = k+ntg
       dgdgn(k:kmax) = gnew(ig_low(iseg):ig_mid(iseg),-1,iglomod(iseg))
       k = kmax+1

       if (nsegments(it,ik) > 1) then
          do iseg = 2, nsegments(it,ik)
             kmax = k+ntg-1
             dgdgn(k:kmax) = gnew(ig_low(iseg)+1:ig_mid(iseg),-1,iglomod(iseg))
             k = kmax+1
          end do
       end if

       ! will need to remove 'real' operator below later -- MAB
       ! this will require tackling the inversion of a complex matrix gresponse as well
       dgdgn(:ir_up(it,ik)) = matmul(m_mat(:ir_up(it,ik),:ir_up(it,ik),iglo),real(dgdgn(:ir_up(it,ik))))

       k=1
       iseg=1 ; kmax=k+ntg
       source_mod(k:kmax-1) = real(source(ig_low(iseg):ig_mid(iseg)-1,-1,iglomod(iseg))) &
            - real(dgdgn(k:kmax-1))
       source_mod(kmax) = real(source0(iseg,iglomod(iseg)))
       k = kmax+1

       if (nsegments(it,ik) > 1) then
          do iseg = 2, nsegments(it,ik)
             kmax = k+ntg-1
             source_mod(k:kmax-1) = real(source(ig_low(iseg)+1:ig_mid(iseg)-1,-1,iglomod(iseg))) &
                  - real(dgdgn(k:kmax-1))
             source_mod(kmax) = real(source0(iseg,iglomod(iseg)))
             k = kmax+1
          end do
       end if
    end if

    deallocate (dgdgn, iglomod)

  end subroutine implicit_sweep

  ! implicit_sweep_right starts with a boundary condition for g^{n+1} along the line (theta<0,vpa=0)
  ! as well as a BC at (theta=-pi,vpa>0) and at (theta>0,vpa=vpa_max), and solves for g^{n+1}
  ! note that dvpa/dt<0 for theta>0 and dvpa/dt>0 for theta<0.
  subroutine implicit_sweep_right (iglo, iseg)

    use dist_fn_arrays, only: gnew, source
    use theta_grid, only: ntgrid
    use vpamu_grids, only: nvgrid

    implicit none

    integer, intent (in) :: iglo, iseg

    integer :: ig, iv

    ! first get gnew(theta=vpa=0), which is decoupled from other points
    gnew(ig_mid(iseg),0,iglo) = source0(iseg,iglo)

    ! boundary condition for vpa > 0 is g(theta=-infinity) = 0
    ! boundary condition for theta < 0 (dvpa/dt > 0) is g(vpa=-infinity)=0
    do iv = 0, nvgrid-1
       do ig = ig_low(iseg), ig_mid(iseg)-1
          gnew(ig+1,iv+1,iglo) = (source(ig,iv,iglo) - (gnew(ig,iv+1,iglo)*pmpfac(ig,iv,iglo) &
               + gnew(ig+1,iv,iglo)*ppmfac(ig,iv,iglo) &
               + gnew(ig,iv,iglo)*pmmfac(ig,iv,iglo))) / pppfac(ig,iv,iglo)
       end do
    end do

    ! boundary condition for vpa > 0 particles is g(theta=-infinity) = 0
    ! boundary condition for theta > 0 (corresponds to dvpa/dt < 0) is 
    ! g(vpa=infinity) = 0
    do iv = nvgrid-1, 0, -1
       do ig = ig_mid(iseg), ig_up(iseg)-1
          gnew(ig+1,iv,iglo) = (source(ig,iv,iglo) - (gnew(ig+1,iv+1,iglo)*pppfac(ig,iv,iglo) &
               + gnew(ig,iv+1,iglo)*pmpfac(ig,iv,iglo) &
               + gnew(ig,iv,iglo)*pmmfac(ig,iv,iglo))) / ppmfac(ig,iv,iglo)
       end do
    end do
       
  end subroutine implicit_sweep_right

  ! implicit_sweep starts with a boundary condition for g^{n+1} along the line (theta>0,vpa=0)
  ! as well as BCs at (theta=theta_max,vpa<0) and (theta<0,vpa=-vpa_max), and solves for g^{n+1}
  ! note that dvpa/dt<0 for theta>0 and dvpa/dt>0 for theta<0.
  subroutine implicit_sweep_left (iglo, iseg)

    use dist_fn_arrays, only: gnew, source
    use theta_grid, only: ntgrid
    use vpamu_grids, only: nvgrid

    implicit none

    integer, intent (in) :: iglo, iseg

    integer :: ig, iv

    ! boundary condition for vpa < 0 particles is g(theta=infinity) = 0
    ! boundary condition for theta > 0 (dvpa/dt < 0) is g(vpa=infinity)=0
    do iv = -1, -nvgrid, -1
       do ig = ig_up(iseg)-1, ig_mid(iseg), -1
          gnew(ig,iv,iglo) = (source(ig,iv,iglo) - (gnew(ig+1,iv+1,iglo)*pppfac(ig,iv,iglo) &
               + gnew(ig+1,iv,iglo)*ppmfac(ig,iv,iglo) &
               + gnew(ig,iv+1,iglo)*pmpfac(ig,iv,iglo))) / pmmfac(ig,iv,iglo)
       end do
    end do
    
    ! boundary condition for vpa < 0 is g(theta=infinity) = 0
    ! boundary condition for theta < 0 (dvpa/dt > 0) is g(vpa=-infinity)=0
    do iv = -nvgrid, -2
       do ig = ig_mid(iseg)-1, ig_low(iseg), -1
          gnew(ig,iv+1,iglo) = (source(ig,iv,iglo) - (gnew(ig+1,iv+1,iglo)*pppfac(ig,iv,iglo) &
               + gnew(ig+1,iv,iglo)*ppmfac(ig,iv,iglo) &
               + gnew(ig,iv,iglo)*pmmfac(ig,iv,iglo))) / pmpfac(ig,iv,iglo)
       end do
    end do
    
  end subroutine implicit_sweep_left

  subroutine implicit_solve

    use dist_fn_arrays, only: gnew
    use gs2_layouts, only: g_lo, it_idx, ik_idx
    use theta_grid, only: dbdthetc, theta, ntgrid, ntheta
    use vpamu_grids, only: nvgrid

    implicit none

    integer :: iglo, iseg, k, kmax, ntg, it, itmod, ig, ik

    integer, dimension (:), allocatable :: iglomod
    real, dimension (:), allocatable :: source_mod, tmp

    ntg = ntheta/2

    if (.not. allocated(source_mod)) then
       allocate (tmp(ntg*nseg_max+1)) ; tmp = 0.0
       allocate (source_mod(ntg*nseg_max+1)) ; source_mod = 0.0
    end if

    call get_source

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ! initially set g^{n+1}(theta<0,vpa=0) to zero
       where (dbdthetc(:ntgrid-1,1) < -epsilon(0.))
          gnew(:ntgrid-1,0,iglo) = 0.0
       end where
    end do

    do iglo = g_lo%llim_proc, g_lo%ulim_proc

       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)

       if (it > it_shift(ik)) cycle

       allocate (iglomod(nsegments(it,ik)))

       ! sweep through vpa and theta once to obtain the response of 
       ! particles with vpa=-dvpa below and including the midplane
       ! to g^{n} when g^{n+1} for particles with vpa=0 below and including the midplane is zero
       call implicit_sweep (iglo, source_mod=source_mod)
       
       ! get g^{n+1}(theta<=0,vpa=0) using response matrix
       tmp = matmul(gresponse(:ir_up(it,ik),:ir_up(it,ik),iglo), source_mod(:ir_up(it,ik)))

       k=1
       iseg=1 ; kmax=k+ntg
       iglomod(iseg) = iglo + it_shift_left(it)
       itmod = it + it_shift_left(it)
       gnew(ig_low(iseg):ig_mid(iseg),0,iglomod(iseg)) = tmp(k:kmax)

       k = kmax+1
       if (nsegments(it,ik) > 1) then
          do iseg = 2, nsegments(it,ik)
             iglomod(iseg) = iglomod(iseg-1) + iglo_shift(itmod,ik)
             itmod = itmod + iglo_shift(itmod,ik)
             kmax = k+ntg-1
             ! set g(-pi) for this segment equal to g(pi) at the previous segment
             gnew(ig_low(iseg),0,iglomod(iseg)) = gnew(ig_up(iseg-1),0,iglomod(iseg-1))
             ! fill in the remaining theta values in this segment
             gnew(ig_low(iseg)+1:ig_mid(iseg),0,iglomod(iseg)) = tmp(k:kmax)
             k = kmax+1
          end do
       end if

       ! with g^{n+1}(theta=0,vpa>0) specified, sweep once more to get g^{n+1} everywhere else
       call implicit_sweep (iglo)

       deallocate (iglomod)
       
    end do

    if (allocated (tmp)) deallocate (source_mod, tmp)

  end subroutine implicit_solve

  subroutine get_source

    use dist_fn_arrays, only: gold, source
    use gs2_time, only: code_dt
    use theta_grid, only: ntgrid
    use vpamu_grids, only: nvgrid
    use gs2_layouts, only: g_lo, it_idx, ik_idx

    implicit none

    integer :: ig, iv, iglo, iseg, it, iglomod, itmod, ik

    ! note that GKE is evaluated at cell values
    ! ig=-ntgrid is first cell value, ig=ntgrid-1 is last cell value
    ! iv=-nvgrid is first cell value, iv=nvgrid-1 is last cell value

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do iv = -nvgrid, nvgrid-1
          do ig = -ntgrid, ntgrid-1
             source(ig,iv,iglo) = -(gold(ig+1,iv+1,iglo)*mppfac(ig,iv,iglo) &
                  + gold(ig+1,iv,iglo)*mpmfac(ig,iv,iglo) &
                  + gold(ig,iv+1,iglo)*mmpfac(ig,iv,iglo) &
                  + gold(ig,iv,iglo)*mmmfac(ig,iv,iglo))
          end do
       end do
    end do

    ! treat particles with vpa=0 at outboard midplane (totally trapped particles) specially, 
    ! as they are decoupled from other theta-vpa points in the implicit solve
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       if (it > it_shift(ik)) cycle
       iglomod = iglo + it_shift_left(it)
       itmod = it + it_shift_left(it)
       do iseg = 1, nsegments(it,ik)
          ! this is the source at the (theta=0,vpa=0) grid point (not a cell value)
          ! it is needed because that point does not take info from other grid points
          ! except indirectly through the fields
          source0(iseg,iglomod) = gold(ig_mid(iseg),0,iglomod)
          iglomod = iglomod + iglo_shift(itmod,ik)
          itmod = itmod + iglo_shift(itmod,ik)
       end do
    end do

  end subroutine get_source

!   subroutine init_connected_bc

!     use theta_grid, only: ntgrid, nperiod, ntheta, theta
!     use kt_grids, only: naky, ntheta0, aky, theta0
!     use le_grids, only: ng2, nlambda
!     use gs2_layouts, only: g_lo, ik_idx, it_idx, imu_idx, is_idx
!     use gs2_layouts, only: idx, proc_id
!     use mp, only: iproc, nproc, max_allreduce, proc0
!     use constants
!     use redistribute, only: index_list_type, init_fill, delete_list

!     implicit none

!     type (index_list_type), dimension(0:nproc-1) :: to, from
!     integer, dimension (0:nproc-1) :: nn_from, nn_to
!     integer, dimension (3) :: to_low, from_low, to_high, from_high
!     integer :: ik, it, imu, is, iglo, it0, itl, itr, jshift0
!     integer :: ip, ipleft, ipright
!     integer :: iglo_left, iglo_right, i, j, k
!     integer :: iglo_star, it_star, ncell
!     integer :: n, n_links_max, nn_max
!     integer :: ng
!     integer, dimension(naky*ntheta0) :: n_k

!     if (connectinit) return
!     connectinit = .true.

!     ng = ntheta/2 + (nperiod-1)*ntheta

!     ! jshift0 corresponds to J (not delta j) from Beer thesis (unsure about +0.01) -- MAB
!     ! delta j is the number of akxs (i.e. theta0s) 'skipped' when connecting one 
!     ! parallel segment to the next to satisfy the twist and shift
!     ! boundary conditions. delta j = (ik - 1) * jshift0 . EGH

!     if (naky > 1 .and. ntheta0 > 1) then
!        jshift0 = int((theta(ng)-theta(-ng))/(theta0(2,2)-theta0(1,2)) + 0.01)
!     else if (naky == 1 .and. ntheta0 > 1 .and. aky(1) /= 0.0) then
!        jshift0 = int((theta(ng)-theta(-ng))/(theta0(2,1)-theta0(1,1)) + 0.01)
!     else
!        jshift0 = 1
!     end if
    
!     allocate (itleft(naky,ntheta0), itright(naky,ntheta0))
!     itleft(1,:) = -1 ! No connected segments when ky = 0. EGH
!     itright(1,:) = -1
!     do ik = 1, naky
!        do it = 1, ntheta0
!           if (it > (ntheta0+1)/2) then
!              ! akx is negative (-1 shift because akx(1) = 0) -- MAB
!              it0 = it - ntheta0 - 1
!           else
!              ! akx is positive (-1 shift because akx(1) = 0) -- MAB
!              it0 = it - 1
!           end if

!           if (ik == 1) then
!              if (aky(ik) /= 0.0 .and. naky == 1) then
!                 ! for this case, jshift0 is delta j from Beer thesis -- MAB
!                 itl = it0 + jshift0
!                 itr = it0 - jshift0
!              else
!                 itl = it0
!                 itr = it0
!              end if
!           else
!              ! ik = 1 corresponds to aky=0, so shift index by 1
!              ! (ik-1)*jshift0 is delta j from Beer thesis -- MAB
!              itl = it0 + (ik-1)*jshift0
!              itr = it0 - (ik-1)*jshift0
!           end if

!           ! itleft and itright are the kx indices to which the
!           ! local kx connects to on left and right - MAB

!           ! remap to usual GS2 indices -- MAB
!           if (itl >= 0 .and. itl < (ntheta0+1)/2) then
!              itleft(ik,it) = itl + 1
!           else if (itl + ntheta0 + 1 > (ntheta0+1)/2 &
!              .and. itl + ntheta0 + 1 <= ntheta0) then
!              itleft(ik,it) = itl + ntheta0 + 1
!           else
!              ! for periodic, need to change below -- MAB/BD
!              ! j' = j + delta j not included in simulation, so can't connect -- MAB
!              itleft(ik,it) = -1
!           end if

!           ! same stuff for j' = j - delta j -- MAB
!           if (itr >= 0 .and. itr < (ntheta0+1)/2) then
!              itright(ik,it) = itr + 1
!           else if (itr + ntheta0 + 1 > (ntheta0+1)/2 &
!              .and. itr + ntheta0 + 1 <= ntheta0) then
!              itright(ik,it) = itr + ntheta0 + 1
!           else
!              ! for periodic, need to change below -- MAB/BD
!              itright(ik,it) = -1
!           end if
!        end do
!     end do

!     allocate (connections(g_lo%llim_proc:g_lo%ulim_alloc))

!     do iglo = g_lo%llim_proc, g_lo%ulim_proc
!        ik = ik_idx(g_lo,iglo)
!        it = it_idx(g_lo,iglo)
!        imu = imu_idx(g_lo,iglo)
!        is = is_idx(g_lo,iglo)

!        ! get processors and indices for j' (kx') modes connecting
!        ! to mode j (kx) so we can set up communication -- MAB

!        if (itleft(ik,it) < 0) then
!           connections(iglo)%iproc_left = -1
!           connections(iglo)%iglo_left = -1
!        else
!           connections(iglo)%iproc_left &
!                = proc_id(g_lo,idx(g_lo,ik,itleft(ik,it),imu,is))
!           connections(iglo)%iglo_left &
!                = idx(g_lo,ik,itleft(ik,it),imu,is)
!        end if
!        if (itright(ik,it) < 0) then
!           connections(iglo)%iproc_right = -1
!           connections(iglo)%iglo_right = -1
!        else
!           connections(iglo)%iproc_right &
!                = proc_id(g_lo,idx(g_lo,ik,itright(ik,it),imu,is))
!           connections(iglo)%iglo_right &
!                = idx(g_lo,ik,itright(ik,it),imu,is)
!        end if
!        if (connections(iglo)%iproc_left >= 0 .or. &
!             connections(iglo)%iproc_right >= 0) then
!           connections(iglo)%neighbor = .true.
!        else
!           connections(iglo)%neighbor = .false.
!        end if
!     end do

!     if (boundary_option_switch == boundary_option_linked) then
!        do iglo = g_lo%llim_proc, g_lo%ulim_proc          
!           ik = ik_idx(g_lo,iglo)
!           il = il_idx(g_lo,iglo)
!           ! if connecting to the left, save vpa>0 part of g
!           if (connections(iglo)%iglo_left >= 0 .and. aky(ik) /= 0.0) &
!                save_h (1,iglo) = .true.
!           ! if connecting to the right, save vpa<0 part of g
!           if (connections(iglo)%iglo_right >= 0 .and. aky(ik) /= 0.0) &
!                save_h (2,iglo) = .true.
!        end do

!        allocate (l_links(naky, ntheta0))
!        allocate (r_links(naky, ntheta0))
!        allocate (n_links(2, naky, ntheta0))

!        n_links_max = 0
!        do it = 1, ntheta0
!           do ik = 1, naky
! ! count the links for each region
! ! l_links = number of links to the left
!              l_links(ik, it) = 0
!              it_star = it
!              do 
!                 ! if no more connections on left, exit
!                 if (it_star == itleft(ik, it_star)) exit
!                 ! if more connections on left, increment l_links
!                 if (itleft(ik, it_star) >= 0) then
!                    l_links(ik, it) = l_links(ik, it) + 1
!                    it_star = itleft(ik, it_star)
!                    if (l_links(ik, it) > 5000) then
! ! abort by deallocating twice
!                       write(*,*) 'l_links error'
!                       deallocate (l_links)
!                       deallocate (l_links)
!                    endif
!                 else
!                    exit
!                 end if
!              end do
! ! r_links = number of links to the right
!              r_links(ik, it) = 0
!              it_star = it
!              do 
!                 ! if no more connections to right, exit
!                 if (it_star == itright(ik, it_star)) exit
!                 ! if more connections to right, increment r_links
!                 if (itright(ik, it_star) >= 0) then
!                    r_links(ik, it) = r_links(ik, it) + 1
!                    it_star = itright(ik, it_star)
!                    if (r_links(ik, it) > 5000) then
! ! abort by deallocating twice
!                       write(*,*) 'r_links error'
!                       deallocate (r_links)
!                       deallocate (r_links)
!                    endif
!                 else
!                    exit
!                 end if
!              end do
! ! 'n_links' complex numbers are needed to specify bc for (ik, it) region
! ! ignoring wfb
! ! n_links(1,:,:) is for v_par > 0, etc.
!              if (l_links(ik, it) == 0) then
!                 n_links(1, ik, it) = 0
!              else 
!                 n_links(1, ik, it) = 2*l_links(ik, it) - 1
!              end if

!              if (r_links(ik, it) == 0) then
!                 n_links(2, ik, it) = 0
!              else 
!                 n_links(2, ik, it) = 2*r_links(ik, it) - 1
!              end if
!              n_links_max = max(n_links_max, n_links(1,ik,it), n_links(2,ik,it))
!           end do
!        end do
! ! MAB FLAG - don't think wfb coding necessary anymore as won't treat vpa (theta = +/- pi) = 0 specially
! ! wfb
! !       if (n_links_max > 0) n_links_max = n_links_max + 3
       
! ! now set up communication pattern:
! ! excluding wfb

!        nn_to = 0
!        nn_from = 0 

!        do iglo = g_lo%llim_world, g_lo%ulim_world

!           ip = proc_id(g_lo,iglo)
!           ik = ik_idx(g_lo,iglo)
!           it = it_idx(g_lo,iglo)
          
!           ncell = r_links(ik, it) + l_links(ik, it) + 1
!           if (ncell == 1) cycle

!           iglo_star = iglo
!           do j = 1, r_links(ik, it)
!              call get_right_connection (iglo_star, iglo_right, ipright)
! ! sender
!              if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
! ! receiver
!              if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1

!              iglo_star = iglo_right
!           end do
             
!           iglo_star = iglo
!           do j = 1, l_links(ik, it)
!              call get_left_connection (iglo_star, iglo_left, ipleft)
! ! sender
!              if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
! ! receiver
!              if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
!              iglo_star = iglo_left
!           end do

!        end do
       
!        nn_max = maxval(nn_to)
!        call max_allreduce (nn_max)
!        if (nn_max == 0) then
!           no_comm = .true.
!           goto 200
!        end if

!        do ip = 0, nproc-1
!           if (nn_from(ip) > 0) then
!              allocate (from(ip)%first(nn_from(ip)))
!              allocate (from(ip)%second(nn_from(ip)))
!              allocate (from(ip)%third(nn_from(ip)))
!           endif
!           if (nn_to(ip) > 0) then
!              allocate (to(ip)%first(nn_to(ip)))
!              allocate (to(ip)%second(nn_to(ip)))
!              allocate (to(ip)%third(nn_to(ip)))
!           endif
!        end do
       
!        nn_from = 0
!        nn_to = 0          

!        do iglo = g_lo%llim_world, g_lo%ulim_world

!           ip = proc_id(g_lo,iglo)
!           ik = ik_idx(g_lo,iglo)
!           it = it_idx(g_lo,iglo)
          
!           ncell = r_links(ik, it) + l_links(ik, it) + 1
!           if (ncell == 1) cycle
          
!           iglo_star = iglo
!           do j = 1, r_links(ik, it)
!              call get_right_connection (iglo_star, iglo_right, ipright)
! ! sender
!              if (ip == iproc) then
!                 n = nn_from(ipright) + 1
!                 nn_from(ipright) = n
!                 from(ipright)%first(n) = ntgrid
!                 from(ipright)%second(n) = 1
!                 from(ipright)%third(n) = iglo
!              end if
! ! receiver
!              if (ipright == iproc) then
!                 n = nn_to(ip) + 1
!                 nn_to(ip) = n
!                 to(ip)%first(n) = j 
!                 to(ip)%second(n) = 1
!                 to(ip)%third(n) = iglo_right
!              end if
!              iglo_star = iglo_right
!           end do
             
!           iglo_star = iglo
!           do j = 1, l_links(ik, it)
!              call get_left_connection (iglo_star, iglo_left, ipleft)
! ! sender
!              if (ip == iproc) then
!                 n = nn_from(ipleft) + 1
!                 nn_from(ipleft) = n
!                 from(ipleft)%first(n) = -ntgrid
!                 from(ipleft)%second(n) = 2
!                 from(ipleft)%third(n) = iglo
!              end if
! ! receiver
!              if (ipleft == iproc) then
!                 n = nn_to(ip) + 1
!                 nn_to(ip) = n
!                 to(ip)%first(n) = j 
!                 to(ip)%second(n) = 2
!                 to(ip)%third(n) = iglo_left
!              end if
!              iglo_star = iglo_left
!           end do
!        end do

!        from_low (1) = -ntgrid
!        from_low (2) = 1
!        from_low (3) = g_lo%llim_proc
       
!        from_high(1) = ntgrid
!        from_high(2) = 2
!        from_high(3) = g_lo%ulim_alloc

!        to_low (1) = 1
!        to_low (2) = 1 
!        to_low (3) = g_lo%llim_proc
       
!        to_high(1) = n_links_max
!        to_high(2) = 2
!        to_high(3) = g_lo%ulim_alloc

!        call init_fill (links_p, 'c', to_low, to_high, to, &
!             from_low, from_high, from)
       
!        call delete_list (from)
!        call delete_list (to)
       
! ! take care of wfb

!        nn_to = 0
!        nn_from = 0 

!        do iglo = g_lo%llim_world, g_lo%ulim_world

!           il = il_idx(g_lo,iglo)
!           if (il /= ng2+1) cycle

!           ip = proc_id(g_lo,iglo)
!           ik = ik_idx(g_lo,iglo)
!           it = it_idx(g_lo,iglo)
          
!           ncell = r_links(ik, it) + l_links(ik, it) + 1
!           if (ncell == 1) cycle
             
!           iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

! ! v_par > 0:
!           call find_leftmost_link (iglo, iglo_right, ipright)
!           do j = 1, ncell
! ! sender
!              if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
! ! receiver
!              if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1
!              call get_right_connection (iglo_right, iglo_right, ipright)
!           end do
             
! ! v_par < 0:
!           call find_rightmost_link (iglo, iglo_left, ipleft)
!           do j = 1, ncell
! ! sender
!              if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
! ! receiver
!              if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
!              call get_left_connection (iglo_left, iglo_left, ipleft)
!           end do
!        end do
       
!        do ip = 0, nproc-1
!           if (nn_from(ip) > 0) then
!              allocate (from(ip)%first(nn_from(ip)))
!              allocate (from(ip)%second(nn_from(ip)))
!              allocate (from(ip)%third(nn_from(ip)))
!           endif
!           if (nn_to(ip) > 0) then
!              allocate (to(ip)%first(nn_to(ip)))
!              allocate (to(ip)%second(nn_to(ip)))
!              allocate (to(ip)%third(nn_to(ip)))
!           endif
!        end do
       
!        nn_from = 0
!        nn_to = 0          

!        do iglo = g_lo%llim_world, g_lo%ulim_world

!           il = il_idx(g_lo,iglo)
!           if (il /= ng2+1) cycle

!           ip = proc_id(g_lo,iglo)
!           ik = ik_idx(g_lo,iglo)
!           it = it_idx(g_lo,iglo)
          
!           ncell = r_links(ik, it) + l_links(ik, it) + 1
!           if (ncell == 1) cycle

!           iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

! ! v_par > 0: 
!           call find_leftmost_link (iglo, iglo_right, ipright)
!           do j = 1, ncell
! ! sender
!              if (ip == iproc) then
!                 n = nn_from(ipright) + 1
!                 nn_from(ipright) = n
!                 from(ipright)%first(n) = ntgrid
!                 from(ipright)%second(n) = 1
!                 from(ipright)%third(n) = iglo
!              end if
! ! receiver
!              if (ipright == iproc) then
!                 n = nn_to(ip) + 1
!                 nn_to(ip) = n
!                 to(ip)%first(n) = r_links(ik, it) + 1
!                 to(ip)%second(n) = 1
!                 to(ip)%third(n) = iglo_right
!              end if
!              call get_right_connection (iglo_right, iglo_right, ipright)
!           end do
             
! ! v_par < 0: 
!           call find_rightmost_link (iglo, iglo_left, ipleft)
!           do j = 1, ncell
! ! sender
!              if (ip == iproc) then
!                 n = nn_from(ipleft) + 1
!                 nn_from(ipleft) = n
!                 from(ipleft)%first(n) = -ntgrid
!                 from(ipleft)%second(n) = 2
!                 from(ipleft)%third(n) = iglo
!              end if
! ! receiver
!              if (ipleft == iproc) then
!                 n = nn_to(ip) + 1
!                 nn_to(ip) = n
!                 to(ip)%first(n) = l_links(ik, it) + 1
!                 to(ip)%second(n) = 2
!                 to(ip)%third(n) = iglo_left
!              end if
!              call get_left_connection (iglo_left, iglo_left, ipleft)
!           end do
!        end do

!        from_low (1) = -ntgrid
!        from_low (2) = 1
!        from_low (3) = g_lo%llim_proc
       
!        from_high(1) = ntgrid
!        from_high(2) = 2
!        from_high(3) = g_lo%ulim_alloc

!        to_low (1) = 1
!        to_low (2) = 1 
!        to_low (3) = g_lo%llim_proc
       
!        to_high(1) = n_links_max
!        to_high(2) = 2
!        to_high(3) = g_lo%ulim_alloc

!        call init_fill (wfb_p, 'c', to_low, to_high, to, from_low, from_high, from)
       
!        call delete_list (from)
!        call delete_list (to)
       
! ! n_links_max is typically 2 * number of cells in largest supercell
!        allocate (g_adj (n_links_max, 2, g_lo%llim_proc:g_lo%ulim_alloc))

! ! now set up links_h:
! ! excluding wfb

!        nn_to = 0
!        nn_from = 0 

!        do iglo = g_lo%llim_world, g_lo%ulim_world

!           il = il_idx(g_lo,iglo)
!           if (nlambda > ng2 .and. il >= ng2+1) cycle

!           ip = proc_id(g_lo,iglo)
!           ik = ik_idx(g_lo,iglo)
!           it = it_idx(g_lo,iglo)
          
!           ncell = r_links(ik, it) + l_links(ik, it) + 1
!           if (ncell == 1) cycle

! ! If this is not the first link in the chain, continue
!           if (l_links(ik, it) > 0) then
! ! For each link to the right, do:
!              iglo_star = iglo
!              do j = 1, r_links(ik, it)
!                 call get_right_connection (iglo_star, iglo_right, ipright)
! ! sender
!                 if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
! ! receiver
!                 if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1
!                 iglo_star = iglo_right
!              end do
!           end if

!           if (r_links(ik, it) > 0) then
! ! For each link to the left, do:
!              iglo_star = iglo
!              do j = 1, l_links(ik, it)
!                 call get_left_connection (iglo_star, iglo_left, ipleft)
! ! sender
!                 if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
! ! receiver
!                 if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
!                 iglo_star = iglo_left
!              end do
!           end if
!        end do
       
!        do ip = 0, nproc-1
!           if (nn_from(ip) > 0) then
!              allocate (from(ip)%first(nn_from(ip)))
!              allocate (from(ip)%second(nn_from(ip)))
!              allocate (from(ip)%third(nn_from(ip)))
!           endif
!           if (nn_to(ip) > 0) then
!              allocate (to(ip)%first(nn_to(ip)))
!              allocate (to(ip)%second(nn_to(ip)))
!              allocate (to(ip)%third(nn_to(ip)))
!           endif
!        end do
       
!        nn_from = 0
!        nn_to = 0          

!        do iglo = g_lo%llim_world, g_lo%ulim_world

!           il = il_idx(g_lo,iglo)
!           if (nlambda > ng2 .and. il >= ng2+1) cycle

!           ip = proc_id(g_lo,iglo)
!           ik = ik_idx(g_lo,iglo)
!           it = it_idx(g_lo,iglo)
          
!           ncell = r_links(ik, it) + l_links(ik, it) + 1
!           if (ncell == 1) cycle

!           if (l_links(ik, it) > 0) then
! ! For each link to the right, do:
!              iglo_star = iglo
!              do j = 1, r_links(ik, it)
! ! get address of link
!                 call get_right_connection (iglo_star, iglo_right, ipright)
! ! sender
!                 if (ip == iproc) then
!                    n = nn_from(ipright) + 1
!                    nn_from(ipright) = n
!                    from(ipright)%first(n) = ntgrid
!                    from(ipright)%second(n) = 1
!                    from(ipright)%third(n) = iglo
!                 end if
! ! receiver
!                 if (ipright == iproc) then
!                    n = nn_to(ip) + 1
!                    nn_to(ip) = n
!                    to(ip)%first(n) = 2*l_links(ik, it) + j
!                    to(ip)%second(n) = 1
!                    to(ip)%third(n) = iglo_right
!                 end if
!                 iglo_star = iglo_right
!              end do
!           end if

!           if (r_links(ik, it) > 0) then
! ! For each link to the left, do:
!              iglo_star = iglo
!              do j = 1, l_links(ik, it)
! ! get address of link
!                 call get_left_connection (iglo_star, iglo_left, ipleft)
! ! sender
!                 if (ip == iproc) then
!                    n = nn_from(ipleft) + 1
!                    nn_from(ipleft) = n
!                    from(ipleft)%first(n) = -ntgrid
!                    from(ipleft)%second(n) = 2   
!                    from(ipleft)%third(n) = iglo
!                 end if
! ! receiver
!                 if (ipleft == iproc) then
!                    n = nn_to(ip) + 1
!                    nn_to(ip) = n
!                    to(ip)%first(n) = 2*r_links(ik, it) + j
!                    to(ip)%second(n) = 2
!                    to(ip)%third(n) = iglo_left
!                 end if
!                 iglo_star = iglo_left
!              end do
!           end if
!        end do

!        from_low (1) = -ntgrid
!        from_low (2) = 1
!        from_low (3) = g_lo%llim_proc
       
!        from_high(1) = ntgrid
!        from_high(2) = 2
!        from_high(3) = g_lo%ulim_alloc

!        to_low (1) = 1
!        to_low (2) = 1 
!        to_low (3) = g_lo%llim_proc
       
!        to_high(1) = n_links_max
!        to_high(2) = 2
!        to_high(3) = g_lo%ulim_alloc

!        call init_fill (links_h, 'c', to_low, to_high, to, &
!             from_low, from_high, from)
       
!        call delete_list (from)
!        call delete_list (to)

! ! now take care of wfb (homogeneous part)

!        nn_to = 0
!        nn_from = 0 

!        do iglo = g_lo%llim_world, g_lo%ulim_world

!           il = il_idx(g_lo,iglo)
!           if (il /= ng2+1) cycle

!           ip = proc_id(g_lo,iglo)
!           ik = ik_idx(g_lo,iglo)
!           it = it_idx(g_lo,iglo)
          
!           ncell = r_links(ik, it) + l_links(ik, it) + 1
!           if (ncell == 1) cycle

!           iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

! ! v_par > 0:
!           call find_leftmost_link (iglo, iglo_right, ipright)
!           do j = 1, ncell
! ! sender
!              if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
! ! receiver
!              if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1
!              call get_right_connection (iglo_right, iglo_right, ipright)
!           end do

! ! v_par < 0:
!           call find_rightmost_link (iglo, iglo_left, ipleft)
!           do j = 1, ncell
! ! sender
!              if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
! ! receiver
!              if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
!              call get_left_connection (iglo_left, iglo_left, ipleft)
!           end do
!        end do
       
!        do ip = 0, nproc-1
!           if (nn_from(ip) > 0) then
!              allocate (from(ip)%first(nn_from(ip)))
!              allocate (from(ip)%second(nn_from(ip)))
!              allocate (from(ip)%third(nn_from(ip)))
!           endif
!           if (nn_to(ip) > 0) then
!              allocate (to(ip)%first(nn_to(ip)))
!              allocate (to(ip)%second(nn_to(ip)))
!              allocate (to(ip)%third(nn_to(ip)))
!           endif
!        end do
       
!        nn_from = 0
!        nn_to = 0          

!        do iglo = g_lo%llim_world, g_lo%ulim_world

!           il = il_idx(g_lo,iglo)
!           if (il /= ng2+1) cycle

!           ip = proc_id(g_lo,iglo)
!           ik = ik_idx(g_lo,iglo)
!           it = it_idx(g_lo,iglo)
          
!           ncell = r_links(ik, it) + l_links(ik, it) + 1
!           if (ncell == 1) cycle

!           iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

! ! v_par > 0:
!           call find_leftmost_link (iglo, iglo_right, ipright)
!           do j = 1, ncell
! ! sender
!              if (ip == iproc) then
!                 n = nn_from(ipright) + 1
!                 nn_from(ipright) = n
!                 from(ipright)%first(n) = ntgrid
!                 from(ipright)%second(n) = 1
!                 from(ipright)%third(n) = iglo
!              end if
! ! receiver
!              if (ipright == iproc) then
!                 n = nn_to(ip) + 1
!                 nn_to(ip) = n
!                 to(ip)%first(n) = 2*ncell - r_links(ik, it)
!                 to(ip)%second(n) = 1
!                 to(ip)%third(n) = iglo_right
!              end if
!              call get_right_connection (iglo_right, iglo_right, ipright)
!           end do
 
! ! v_par < 0:
!           call find_rightmost_link (iglo, iglo_left, ipleft)
!           do j = 1, ncell
! ! sender
!              if (ip == iproc) then
!                 n = nn_from(ipleft) + 1
!                 nn_from(ipleft) = n
!                 from(ipleft)%first(n) = -ntgrid
!                 from(ipleft)%second(n) = 2   
!                 from(ipleft)%third(n) = iglo
!              end if
                
! ! receiver
!              if (ipleft == iproc) then
!                 n = nn_to(ip) + 1
!                 nn_to(ip) = n
!                 to(ip)%first(n) = 2*ncell - l_links(ik, it)
!                 to(ip)%second(n) = 2
!                 to(ip)%third(n) = iglo_left
!              end if
!              call get_left_connection (iglo_left, iglo_left, ipleft)
!           end do
!        end do

!        from_low (1) = -ntgrid
!        from_low (2) = 1
!        from_low (3) = g_lo%llim_proc
       
!        from_high(1) = ntgrid
!        from_high(2) = 2
!        from_high(3) = g_lo%ulim_alloc

!        to_low (1) = 1
!        to_low (2) = 1 
!        to_low (3) = g_lo%llim_proc
       
!        to_high(1) = n_links_max
!        to_high(2) = 2
!        to_high(3) = g_lo%ulim_alloc

!        call init_fill (wfb_h, 'c', to_low, to_high, to, from_low, from_high, from)
       
!        call delete_list (from)
!        call delete_list (to)

! 200    continue

! ! Now set up class arrays for the implicit fields 
! ! i_class classes
! ! N_class(i) = number of linked cells for i_th class 
! ! M_class(i) = number of members in i_th class

! ! First count number of linked cells for each (kx, ky) 
!        k = 1
!        do it = 1, ntheta0
!           do ik = 1, naky
!              n_k(k) = 1 + l_links(ik, it) + r_links(ik, it)
!              k = k + 1
! !             if (proc0) write (*,*) 'init_connected_bc: ',ik, it, l_links(ik, it), r_links(ik, it)
!           end do
!        end do

! ! Count how many unique values of n_k there are.  This is the number 
! ! of classes.

! ! Sort: 
!        do j = 1, naky*ntheta0-1
!           do k = 1, naky*ntheta0-1                       
!              if (n_k(k+1) < n_k(k)) then
!                 i = n_k(k)
!                 n_k(k) = n_k(k+1)
!                 n_k(k+1) = i
!              end if
!           end do
!        end do

! ! Then count:
!        i_class = 1
!        do k = 1, naky*ntheta0-1
!           if (n_k(k) == n_k(k+1)) cycle
!           i_class = i_class + 1
!        end do

! ! Allocate M, N:
!        allocate (M_class(i_class))
!        allocate (N_class(i_class))

! ! Initial values
!        M_class = 1 ; N_class = 0

! ! Fill M, N arrays: 
!        j = 1
!        do k = 2, naky*ntheta0
!           if (n_k(k) == n_k(k-1)) then
!              M_class(j) = M_class(j) + 1
!           else 
!              N_class(j) = n_k(k-1)
!              M_class(j) = M_class(j)/N_class(j)
!              j = j + 1
!           end if
!        end do       
!        j = i_class
!        N_class(j) = n_k(naky*ntheta0)
!        M_class(j) = M_class(j)/N_class(j)

! ! Check for consistency:
       
! ! j is number of linked cells in class structure
!        j = 0
!        do i = 1, i_class
!           j = j + N_class(i)*M_class(i)
!        end do

!        if (j /= naky*ntheta0) then
!           write(*,*) 'PE ',iproc,'has j= ',j,' k= ',naky*ntheta0,' : Stopping'
!           stop
!        end if

!     end if

!   end subroutine init_connected_bc

!   subroutine init_pass_right
!     !Create a pass_right object to pass boundary values to the left connection
!     !for sigma=1 (i.e. rightwards moving particles).
!     call init_pass_ends(pass_right,'r',1,'c')
!   end subroutine init_pass_right

!   subroutine init_pass_left
!     !Create a pass_left object to pass boundary values to the left connection
!     !for sigma=2 (i.e. leftwards moving particles).
!     call init_pass_ends(pass_left,'l',2,'c')
!   end subroutine init_pass_left

!   subroutine init_pass_ends(pass_obj,dir,sigma,typestr,ngsend)
!     !<DD> 01/02/2013: Copy of CMR's init_pass_right (revision 2085) 
!     !but in customisable direction to aid code reuse and with optional
!     !range in theta (i.e. so can pass more than just boundary points)
    
!     !An example use is:
!     ! call init_pass_ends(pass_right,'r',1,'c',3)
!     ! call fill(pass_right,gnew,gnew)
!     !This will pass gnew(ntgrid-2:ntgrid,1,iglo) to 
!     !gnew(-ntgrid-2:-ntgrid,1,iglo_linked)
    
!     !TODO:
!     !      1) May be helpful to be able to pass both sigmas at once (e.g. for explicit scheme)
!     use gs2_layouts, only: g_lo, il_idx, idx, proc_id
!     use le_grids, only: ng2
!     use mp, only: iproc, nproc, max_allreduce, proc0
!     use redistribute, only: index_list_type, init_fill, delete_list, redist_type
!     use theta_grid, only:ntgrid

!     implicit none

!     type (redist_type), intent(inout) :: pass_obj !Redist type object to hold communication logic
!     character(1),intent(in)::dir                  !Character string for direction of communication, should be 'l' for left and 'r' for right
!     character(1),intent(in)::typestr              !Character string for type of data to be communicated. Should be 'c','r','i' or 'l'
!     integer,intent(in) :: sigma                   !Which sigma index to send
!     integer,intent(in),optional::ngsend           !How many theta grid points are we sending

!     !Internal variables
!     type (index_list_type), dimension(0:nproc-1) :: to, from
!     integer, dimension (0:nproc-1) :: nn_from, nn_to
!     integer, dimension(3) :: from_low, from_high, to_low, to_high
!     integer :: il, iglo, ip, iglo_con, ipcon, n, nn_max, j
!     logical :: debug=.false.
!     integer :: bound_sign
!     integer :: local_ngsend

!     !Only applies to linked boundary option so exit if not linked
!     if (boundary_option_switch .ne. boundary_option_linked) return

!     !Handle the direction sign, basically we're either doing
!     !     ntgrid --> -ntgrid (passing to right)
!     !or 
!     !    -ntgrid --> ntgrid  (passing to left)
!     if (dir=='r') then
!        bound_sign=1
!     elseif (dir=='l') then
!        bound_sign=-1
!     else
!        if (proc0) write(6,*) "Invalid direction string passed to init_pass_ends, defaulting to 'r'"
!        bound_sign=1
!     endif

!     !Set the default number of theta grid points to send, if required
!     if (present(ngsend)) then
!        local_ngsend=ngsend
!     else
!        local_ngsend=1
!     endif

!     if (proc0.and.debug) write (6,*) "Initialising redist_type with settings Direction : ",dir," sigma",sigma,"local_ngsend",local_ngsend
    
!     ! Need communications to satisfy || boundary conditions
!     ! First find required blocksizes 
    
!     !Initialise variables used to count how many entries to send and receive
!     !for each processor
!     nn_to = 0   ! # nn_to(ip) = communicates from ip TO HERE (iproc)
!     nn_from = 0 ! # nn_from(ip) = communicates to ip FROM HERE (iproc)
    
!     !Now loop over >all< iglo indices and work out how much data needs to be sent and received by each processor
!     !Hence this routine does not scale with number of procs, see updated redist object creation in ccfe_opt_test (~r2173)
!     !for examples which do scale
!     do iglo = g_lo%llim_world, g_lo%ulim_world
!        !Get the lambda index so we can skip trapped particles
!        il = il_idx(g_lo,iglo)
       
!        !Exclude disconnected trapped particles
!        if (il > ng2+1) cycle     
       
!        !Get the processor id for the proc holding the current iglo index
!        ip = proc_id(g_lo,iglo)
       
!        !What iglo connects to the current one in the direction of interest (and what proc is it on)
!        !Note ipcon is <0 if no connection in direction of interest
!        if (bound_sign==1) then
!           call get_right_connection (iglo, iglo_con, ipcon)
!        else
!           call get_left_connection (iglo, iglo_con, ipcon)
!        endif
       
!        !Is the connected tube's proc the current processor?
!        if (ipcon .eq. iproc ) then
!           !If so add an entry recording an extra piece of information is to be sent
!           !to proc ip (the proc holding iglo) from this proc (ipcon)
!           !Note: Here we assume theta grid points are all on same proc 
!           nn_to(ip)=nn_to(ip)+local_ngsend
!        endif
       
!        !Is the proc holding iglo this proc and is there a connection in the direction of interest?
!        if (ip .eq. iproc .and. ipcon .ge. 0 ) then
!           !If so add an entry recording an extra piece of information is to be sent
!           !from this proc (ipcon) to ip
!           !Note: Here we assume theta grid points are all on same proc 
!           nn_from(ipcon)=nn_from(ipcon)+local_ngsend
!        endif
!     end do
    
!     !Find the maximum amount of data to be received by a given processor
!     !(first do it locally)
!     nn_max = maxval(nn_to)
!     !(now do it globally)
!     call max_allreduce (nn_max)
    
!     !Bit of debug printing
!     if (proc0.and.debug) then
!        write(6,*) 'init_pass_ends (1) processor, nn_to, nn_from:',iproc,nn_to,nn_from
!     endif

!     !Now that we've worked out how much data needs to be sent and received, define what specific
!     !data needs to be sent to where
!     if (nn_max.gt.0) then
!        ! 
!        ! CMR, 25/1/2013: 
!        !      communication required to satisfy linked BC
!        !      allocate indirect addresses for sends/receives 
!        !     
!        !      NB communications use "c_fill_3" as g has 3 indices
!        !      but 1 index sufficient as only communicating g(ntgrid,1,*)! 
!        !      if "c_fill_1" in redistribute we'd only need allocate: from|to(ip)%first 
!        !                             could be more efficient
!        !  
!        !<DD>, 06/01/2013: This redist object consists of a buffer of length n to hold the 
!        !data during transfer and (currently) 3*2 integer arrays each of length n to hold
!        !the indices of sent and received data.
!        !By using c_fill_1 2*2 of these integer arrays would be removed. Assuming a double complex
!        !buffer and long integer indices a 4n long array saving would be equivalent to the buffer
!        !size and as such should represent a good memory saving but would not effect the amount
!        !of data communicated (obviously).

!        !Create to and from list objects for each processor and 
!        !create storage to hold information about each specific from/to
!        !communication
!        do ip = 0, nproc-1
!           !If proc ip is sending anything to this processor (iproc)
!           if (nn_from(ip) > 0) then
!              allocate (from(ip)%first(nn_from(ip)))
!              allocate (from(ip)%second(nn_from(ip)))
!              allocate (from(ip)%third(nn_from(ip)))
!           endif
!           !If proc ip is receiving anything from this processor (iproc)
!           if (nn_to(ip) > 0) then
!              allocate (to(ip)%first(nn_to(ip)))
!              allocate (to(ip)%second(nn_to(ip)))
!              allocate (to(ip)%third(nn_to(ip)))
!           endif
!        end do

!        !Now fill the indirect addresses...
       
!        !Initialise counters used to record how many pieces of data to expect
!        nn_from = 0 ; nn_to = 0
       
!        !Loop over >all< iglo indices
!        do iglo = g_lo%llim_world, g_lo%ulim_world
!           !Get the lambda index so we can skip trapped particles
!           il = il_idx(g_lo,iglo)
          
!           !Exclude disconnected trapped particles
!           if (il > ng2+1) cycle     
          
!           !What's the processor for the current iglo
!           ip = proc_id(g_lo,iglo)
          
!           !What iglo connects to the current one in the direction of interest (and what proc is it on)?
!           !Note ipcon is <0 if no connection in direction of interest
!           if (bound_sign==1) then
!              call get_right_connection (iglo, iglo_con, ipcon)
!           else
!              call get_left_connection (iglo, iglo_con, ipcon)
!           endif
          
!           !For current proc for current iglo if there's connections in direction
!           !then add an entry to the connected procs list of data to expect
!           if (ip .eq. iproc .and. ipcon .ge. 0 ) then
!              !Loop over theta grid indices
!              !Note: Here we assume theta grid points are all on same proc 
!              !Note: By looping over theta inside iglo loop we should optimise
!              !      memory access compared to looping over theta outside.
!              do j=0,local_ngsend-1
!                 n=nn_from(ipcon)+1 ; nn_from(ipcon)=n
!                 from(ipcon)%first(n)=bound_sign*(ntgrid-j) !Which theta point to send
!                 from(ipcon)%second(n)=sigma     !Sigma grid index to send
!                 from(ipcon)%third(n)=iglo   !iglo index to send
!              enddo
!           endif
          
!           !If target iglo (iglo_con) is on this processor then add an entry recording where
!           !we need to put the data when we receive it.
!           if (ipcon .eq. iproc ) then 
!              !Loop over theta grid indices
!              !Note: Here we assume theta grid points are all on same proc 
!              do j=0,local_ngsend-1
!                 n=nn_to(ip)+1 ; nn_to(ip)=n
!                 to(ip)%first(n)=-bound_sign*(ntgrid-j) !Which theta to store received data
!                 to(ip)%second(n)=sigma       !Sigma grid index to store received data
!                 to(ip)%third(n)=iglo_con !iglo index to store received data
!              enddo
!           endif
!        end do
       
!        !Bit of debug printing,
!        if (debug.and.proc0) then
!           write(6,*) 'init_pass_ends (2) processor, nn_to, nn_from:',iproc,nn_to,nn_from
!        endif

!        !Set data ranges for arrays to be passed, not this just effects how
!        !arrays are indexed, not how big the buffer is.
!        from_low(1)=-ntgrid ; from_low(2)=1  ; from_low(3)=g_lo%llim_proc       
!        from_high(1)=ntgrid ; from_high(2)=2 ; from_high(3)=g_lo%ulim_alloc
!        to_low(1)=-ntgrid   ; to_low(2)=1    ; to_low(3)=g_lo%llim_proc       
!        to_high(1)=ntgrid   ; to_high(2)=2   ; to_high(3)=g_lo%ulim_alloc
       
!        !Initialise fill object
!        call init_fill (pass_obj, typestr, to_low, to_high, to, &
!             from_low, from_high, from)
       
!        !Clean up lists
!        call delete_list (from)
!        call delete_list (to)
!     endif
!   end subroutine init_pass_ends

!   subroutine get_left_connection (iglo, iglo_left, iproc_left)
!     use gs2_layouts, only: g_lo, proc_id, idx
!     use gs2_layouts, only: ik_idx, it_idx, il_idx, ie_idx, is_idx
!     implicit none
!     integer, intent (in) :: iglo
!     integer, intent (out) :: iglo_left, iproc_left
!     integer :: ik, it, il, ie, is

!     ik = ik_idx(g_lo,iglo)
!     it = it_idx(g_lo,iglo)
    
!     if (itleft(ik,it) < 0) then
!        iglo_left = -1
!        iproc_left = -1
!        return
!     end if

!     il = il_idx(g_lo,iglo)
!     ie = ie_idx(g_lo,iglo)
!     is = is_idx(g_lo,iglo)

!     iglo_left = idx(g_lo,ik,itleft(ik,it),il,ie,is)
!     iproc_left = proc_id(g_lo,iglo_left)
!   end subroutine get_left_connection

!   subroutine get_right_connection (iglo, iglo_right, iproc_right)
!     use gs2_layouts, only: g_lo, proc_id, idx
!     use gs2_layouts, only: ik_idx, it_idx, il_idx, ie_idx, is_idx
!     implicit none
!     integer, intent (in) :: iglo
!     integer, intent (out) :: iglo_right, iproc_right
!     integer :: ik, it, il, ie, is

!     ik = ik_idx(g_lo,iglo)
!     it = it_idx(g_lo,iglo)
    
!     if (itright(ik,it) < 0) then
!        iglo_right = -1
!        iproc_right = -1
!        return
!     end if

!     il = il_idx(g_lo,iglo)
!     ie = ie_idx(g_lo,iglo)
!     is = is_idx(g_lo,iglo)

!     iglo_right = idx(g_lo,ik,itright(ik,it),il,ie,is)
!     iproc_right = proc_id(g_lo,iglo_right)
!   end subroutine get_right_connection

  subroutine allocate_arrays

    use kt_grids, only: naky, ntheta0, box
    use theta_grid, only: ntgrid, shat
    use dist_fn_arrays, only: g, gnew, gold, source
    use dist_fn_arrays, only: kx_shift, theta0_shift   ! MR
    use gs2_layouts, only: g_lo
    use nonlinear_terms, only: nonlin
    use vpamu_grids, only: nvgrid

    implicit none

    if (.not. allocated(g)) then
       allocate (g    (-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (gold (-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (gnew (-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (g0   (-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (source(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
#ifdef LOWFLOW
       allocate (gexp_1(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (gexp_2(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (gexp_3(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
       gexp_1 = 0. ; gexp_2 = 0. ; gexp_3 = 0.
#else
       if (nonlin) then
          allocate (gexp_1(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
          allocate (gexp_2(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
          allocate (gexp_3(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
          gexp_1 = 0. ; gexp_2 = 0. ; gexp_3 = 0.
       else
          allocate (gexp_1(1,1,1), gexp_2(1,1,1), gexp_3(1,1,1))
       end if
#endif
       if (boundary_option_switch == boundary_option_linked) then
          allocate (g_h(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
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

  end subroutine allocate_arrays

  !> This function calculates the distribution function at the next timestep. 
  !! It calculates both the inhomogenous part, gnew, due to the sources
  !! (principly the drive terms and the nonlinear term)
  !! and the homogenous part, g1. The actual evolution of the dist func
  !! is done in the subroutine implicit_solve
  !!
  !! After solving for the new dist funcs, this routine calls hyper_diff, which
  !! adds hyper diffusion if present, and solfp1, from the collisions module,
  !! which adds collisions if present.

  subroutine timeadv (phi, apar, bpar, phinew, aparnew, bparnew, istep, mode)

    use theta_grid, only: ntgrid
!    use le_derivatives, only: vspace_derivatives
    use dist_fn_arrays, only: gnew, g, gold
    use nonlinear_terms, only: add_explicit_terms
    use hyper, only: hyper_diff
    use run_parameters, only: nstep

    use vpamu_grids, only: nvgrid
    use gs2_layouts, only: g_lo

    implicit none

    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    integer, optional, intent (in) :: mode
    integer :: modep

    ! TMP FOR TESTING -- MAB
    integer :: ig, iglo, iv

    modep = 0
    if (present(mode)) modep = mode

    if (istep == nstep) call write_mpdist (gnew, '.gtmp', last=.true.)

    ! Calculate the explicit nonlinear terms
    call add_explicit_terms (gexp_1, gexp_2, gexp_3, &
         phi, apar, bpar, istep)
    ! TMP FOR TESTING -- MAB
!     if (istep > 0) then
!        do iglo = g_lo%llim_proc, g_lo%ulim_proc
!           do ig = 1, nseg_max
!              write (*,'(a8,2e12.4)') 'source0', real(source0(ig,iglo)), aimag(source0(ig,iglo))
!           end do
!        end do
!     end if
    ! Implicit Solve for gnew
    call implicit_solve
!     if (istep > 0) then
!        do iglo = g_lo%llim_proc, g_lo%ulim_proc
!           do iv = -nvgrid, nvgrid
!              do ig = -ntgrid, ntgrid
!                 write (*,'(a11,2e12.4)') 'post-solve', real(gnew(ig,iv,iglo)), aimag(gnew(ig,iv,iglo))
!              end do
!           end do
!        end do
!     end if
    ! Add hyper terms (damping)
!    call hyper_diff (gnew, phinew, bparnew)
    ! Add collisions
!    call vspace_derivatives (gnew, g, g0, phi, apar, bpar, phinew, aparnew, bparnew, modep)

  end subroutine timeadv

  ! communication initializations for exb_shear should be done once and 
  ! redistribute routines should be used  BD

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
    use vpamu_grids, only: nmu, nvgrid
    use species, only: nspec
    use run_parameters, only: fphi, fapar, fbpar
    use dist_fn_arrays, only: kx_shift, theta0_shift
    use gs2_time, only: code_dt, code_dt_old
    use constants, only: twopi    

    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,-nvgrid:,g_lo%llim_proc:), intent (in out) :: g0
    complex, dimension(:,:,:), allocatable :: temp 
    complex, dimension(:,:), allocatable :: temp2
    integer, dimension(1), save :: itmin
    integer :: ierr, j 
    integer :: ik, it, is, imu, iv, to_iglo, from_iglo, to_iproc, from_iproc
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
          allocate(temp2(-ntgrid:ntgrid,abs(j)),temp(-ntgrid:ntgrid,-nvgrid:nvgrid,abs(j)))
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
             do imu=1,nmu

                if (j < 0) then
                   do it = 1, iib
                      from_iglo = idx(g_lo, ik, it, imu, is)
                      if (idx_local (g_lo, from_iglo)) temp(:,:,it) = g0(:,:,from_iglo)
                   end do

                   do it = 1, iit-1                        
                      to_iglo = idx(g_lo, ik, it, imu, is)
                      from_iglo = idx(g_lo, ik, it-j, imu, is)

                      if (idx_local(g_lo, to_iglo).and. idx_local(g_lo, from_iglo)) then
                         g0(:,:,to_iglo) = g0(:,:,from_iglo)
                      else if (idx_local(g_lo, from_iglo)) then
                         do iv = -nvgrid, nvgrid
                            call send(g0(:, iv, from_iglo), proc_id (g_lo, to_iglo))
                         enddo
                      else if (idx_local(g_lo, to_iglo)) then
                         do iv = -nvgrid, nvgrid
                            call receive(g0(:, iv, to_iglo), proc_id (g_lo, from_iglo))
                         enddo
                      endif
                   enddo

                   do it = iit, ntheta0                     
                      to_iglo = idx(g_lo, ik, it, imu, is)
                      from_iglo = idx(g_lo, ik, it-j-ntheta0, imu, is)
                      
                      if (idx_local(g_lo, to_iglo) .and. idx_local(g_lo, from_iglo)) then
                         g0(ileft:,:,to_iglo) = temp(:iright,:,it-iit+1)
                         g0(:ileft-1,:,to_iglo) = 0.0
                      else if (idx_local(g_lo, from_iglo)) then
                         do iv = -nvgrid, nvgrid
                            call send(temp(:, iv, it-iit), proc_id (g_lo, to_iglo))
                         enddo
                      else if (idx_local(g_lo, to_iglo)) then
                         do iv = -nvgrid, nvgrid
                            call receive(z, proc_id (g_lo, from_iglo))
                            g0(ileft:,iv,to_iglo) = z(:iright)
                            g0(:ileft-1,iv,to_iglo) = 0.0
                         enddo
                      endif
                   enddo

                else ! j>0

                   do it = 1, j
                      from_iglo = idx(g_lo, ik, iit+it-1, imu, is)
                      if (idx_local (g_lo, from_iglo)) temp(:,:,it) = g0(:,:,from_iglo)
                   end do

                   do it = ntheta0, iib+1, -1
                      to_iglo = idx(g_lo, ik, it, imu, is)
                      from_iglo = idx(g_lo, ik, it-j, imu, is)

                      if (idx_local(g_lo, to_iglo) .and. idx_local(g_lo, from_iglo)) then
                         g0(:,:,to_iglo) = g0(:,:,from_iglo)
                      else if (idx_local(g_lo, from_iglo)) then
                         do iv = -nvgrid, nvgrid
                            call send(g0(:, iv, from_iglo), proc_id (g_lo, to_iglo))
                         enddo
                      else if (idx_local(g_lo, to_iglo)) then
                         do iv = -nvgrid, nvgrid
                            call receive(g0(:,iv,to_iglo), proc_id (g_lo, from_iglo))
                         enddo
                      endif
                   enddo

                   do it = 1, iib
                      to_iglo = idx(g_lo, ik, it, imu, is)
                      from_iglo = idx(g_lo, ik, iit+it-1, imu, is)

                      if (idx_local(g_lo, to_iglo).and. idx_local(g_lo, from_iglo)) then
                         g0(:iright,:,to_iglo) = temp(ileft:,:,it)
                         g0(iright+1:,:,to_iglo) = 0.0
                      else if (idx_local(g_lo, from_iglo)) then
                         do iv = -nvgrid, nvgrid
                            call send(temp(:, iv, it), proc_id (g_lo, to_iglo))
                         enddo
                      else if (idx_local(g_lo, to_iglo)) then
                         do iv = -nvgrid, nvgrid
                            call receive(z, proc_id (g_lo, from_iglo))
                            g0(:iright,iv,to_iglo) = z(ileft:)
                            g0(iright+1:,iv,to_iglo) = 0.0
                         enddo
                      endif
                   enddo
                endif
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
                do imu=1,nmu

                   do it = 1, ntheta0 + jump(ik)                        

                      to_iglo = idx(g_lo, ik, ikx_indexed(it), imu, is)
                      from_iglo = idx(g_lo, ik, ikx_indexed(it-jump(ik)), imu, is)

                      if (idx_local(g_lo, to_iglo) .and. idx_local(g_lo, from_iglo)) then
                         g0(:,:,to_iglo) = g0(:,:,from_iglo)
                      else if (idx_local(g_lo, from_iglo)) then
                         do iv=-nvgrid,nvgrid
                            call send (g0(:, iv, from_iglo), proc_id (g_lo, to_iglo))
                         enddo
                      else if (idx_local(g_lo, to_iglo)) then
                         do iv=-nvgrid,nvgrid
                            call receive (g0(:, iv, to_iglo), proc_id (g_lo, from_iglo))
                         enddo
                      endif
                   enddo

                   do it = ntheta0 + jump(ik) + 1, ntheta0                     
                      to_iglo = idx(g_lo, ik, ikx_indexed(it), imu, is)
                      if (idx_local (g_lo, to_iglo)) g0(:,:,to_iglo) = 0.
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
                do imu=1,nmu

                   do it = ntheta0, 1+jump(ik), -1
                      
                      to_iglo = idx(g_lo, ik, ikx_indexed(it), imu, is)
                      from_iglo = idx(g_lo, ik, ikx_indexed(it-jump(ik)), imu, is)

                      if (idx_local(g_lo, to_iglo) .and. idx_local(g_lo, from_iglo)) then
                         g0(:,:,to_iglo) = g0(:,:,from_iglo)
                      else if (idx_local(g_lo, from_iglo)) then
                         do iv=-nvgrid,nvgrid
                            call send(g0(:, iv, from_iglo), proc_id(g_lo, to_iglo))
                         enddo
                      else if (idx_local(g_lo, to_iglo)) then
                         do iv=-nvgrid,nvgrid
                            call receive(g0(:, iv, to_iglo), proc_id (g_lo, from_iglo))
                         enddo
                      endif
                   enddo

                   do it = jump(ik), 1, -1
                      to_iglo = idx(g_lo, ik, ikx_indexed(it), imu, is)
                      if (idx_local (g_lo, to_iglo)) g0(:,:,to_iglo) = 0.
                   enddo

                enddo
             enddo
          endif
       enddo
    end if
  end subroutine exb_shear

!   subroutine get_source_term &
!        (phi, apar, bpar, phinew, aparnew, bparnew, istep, &
!         isgn, iglo, sourcefac, source)
! #ifdef LOWFLOW
!     use dist_fn_arrays, only: hneoc, vparterm, wdfac, wstarfac
! #endif
!     use dist_fn_arrays, only: aj0, aj1, vpar, g
!     use theta_grid, only: ntgrid, theta, bmag
!     use kt_grids, only: aky, akx
!     use le_grids, only: nlambda, ng2, lmax, anon, energy, negrid, forbid
!     use species, only: spec, nspec
!     use run_parameters, only: fphi, fapar, fbpar, wunits, tunits, t_imp
!     use gs2_time, only: code_dt
!     use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
!     use nonlinear_terms, only: nonlin
!     use hyper, only: D_res
!     use constants
!     implicit none
!     complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
!     complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
!     integer, intent (in) :: istep
!     integer, intent (in) :: isgn, iglo
!     complex, intent (in) :: sourcefac
!     complex, dimension (-ntgrid:), intent (out) :: source

!     integer :: ig, ik, it, il, ie, is
!     complex, dimension (-ntgrid:ntgrid) :: phigavg, apargavg

! !    call timer (0, 'get_source_term')

!     ik = ik_idx(g_lo,iglo)
!     it = it_idx(g_lo,iglo)
!     il = il_idx(g_lo,iglo)
!     ie = ie_idx(g_lo,iglo)
!     is = is_idx(g_lo,iglo)
! !CMR, 4/8/2011
! ! apargavg and phigavg combine to give the GK EM potential chi. 
! !          chi = phigavg - apargavg*vpa(:,isgn,iglo)*spec(is)%stm
! ! phigavg  = phi J0 + 2 T_s/q_s . vperp^2 bpar/bmag J1/Z
! ! apargavg = apar J0 
! ! Both quantities are decentred in time and evaluated on || grid points
! !
!     phigavg  = ((1.0-t_imp)*phi(:,it,ik) + t_imp*phinew(:,it,ik)) &
!                 *aj0(:,iglo)*fphi &
!              + ((1.0-t_imp)*bpar(:,it,ik) + t_imp*bparnew(:,it,ik))&
!                 *aj1(:,iglo)*fbpar*2.0*vperp2(:,iglo)*spec(is)%tz
!     apargavg = ((1.0-t_imp)*apar(:,it,ik) + t_imp*aparnew(:,it,ik)) &
!                 *aj0(:,iglo)*fapar

! ! source term in finite difference equations
!     call set_source

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! Do matrix multiplications
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     if (isgn == 1) then
!        do ig = -ntgrid, ntgrid-1
!           source(ig) = source(ig) &
!                + b(ig,1,iglo)*g(ig,1,iglo) + a(ig,1,iglo)*g(ig+1,1,iglo)
!        end do
!     else
!        do ig = -ntgrid, ntgrid-1
!           source(ig) = source(ig) &
!                + a(ig,2,iglo)*g(ig,2,iglo) + b(ig,2,iglo)*g(ig+1,2,iglo)
!        end do
!     end if

!     source(ntgrid) = source(-ntgrid)  ! is this necessary?  BD

!   contains

!     subroutine set_source

!       use species, only: spec
!       use theta_grid, only: itor_over_B
!       use mp, only: proc0

!       complex :: apar_p, apar_m, phi_p, phi_m!, bpar_p !GGH added bpar_p
! !      real, dimension(:,:), allocatable, save :: ufac
!       real :: bd, bdfac_p, bdfac_m
!       integer :: i_e, i_s
! !      logical :: first = .true.

! !      if (first) then
!       if (.not. allocated(ufac)) then
! !         first = .false.
!          allocate (ufac(negrid, nspec))
!          do i_e = 1, negrid
!             do i_s = 1, nspec
!                ufac(i_e, i_s) = (2.0*spec(i_s)%uprim &
!                     + spec(i_s)%uprim2*energy(i_e)**(1.5)*sqrt(pi)/4.0)
!             end do
!          end do
!       endif

! ! try fixing bkdiff dependence
!       bd = bkdiff(1)

!       bdfac_p = 1.+bd*(3.-2.*real(isgn))
!       bdfac_m = 1.-bd*(3.-2.*real(isgn))


!       do ig = -ntgrid, ntgrid-1
! !CMR, 4/8/2011:
! ! Some concerns, may be red herrings !
! ! (1) no bakdif factors in phi_m, apar_p, apar_m, vpar !!! 
! !                        (RN also spotted this for apar_p)
! ! (2) can interpolations of products be improved? 
! !
! !  Attempt at variable documentation:
! ! phigavg  = phi J0 + 2 T_s/q_s . vperp^2 bpar/bmag J1/Z
! ! apargavg = apar J0                        (decentered in t) 
! ! NB apargavg and phigavg combine to give the GK EM potential chi
! ! chi = phigavg - apargavg*vpa(:,isgn,iglo)*spec(is)%stm
! ! phi_p = 2 phigavg                      .... (roughly!)
! ! phi_m = d/dtheta (phigavg)*DTHETA 
! ! apar_p = 2 apargavg  
! ! apar_m = 2 d/dt (apar)*DELT  (gets multiplied later by J0 and vpa when included in source)
! ! => phi_p - apar_p*vpa(:,isgn,iglo)*spec(is)%stm = 2 chi  .... (roughly!)  
! ! vparterm = -2.0*vpar (IN ABSENCE OF LOWFLOW TERMS)
! ! wdfac = wdrift + wcoriolis/spec(is)%stm (IN ABSENCE OF LOWFLOW TERMS)
! ! wstarfac = wstar  (IN ABSENCE OF LOWFLOW TERMS)
! ! vpar = q_s/sqrt{T_s m_s} (v_||^GS2). \gradpar(theta)/DTHETA . DELT (centred)
! ! wdrift =    q_s/T_s  v_d.\grad_perp . DELT 
! ! wcoriolis = q_s/T_s  v_C.\grad_perp . DELT 
! !
! ! Definition of source:= 2*code_dt*RHS of GKE
! ! source     appears to contain following physical terms
! !   -2q_s/T_s v||.grad(J0 phi + 2 vperp^2 bpar/bmag J1/Z T_s/q_s).delt 
! !   -2d/dt(q v|| J0 apar / T).delt
! !   +hyperviscosity
! !   -2 v_d.\grad_perp (q J0 phi/T + 2 vperp^2 bpar/bmag J1/Z).delt 
! !   -coriolis terms
! !   2{\chi,f_{0s}}.delt  (allowing for sheared flow)
! !CMRend

!          phi_p = bdfac_p*phigavg(ig+1)+bdfac_m*phigavg(ig)
!          phi_m = phigavg(ig+1)-phigavg(ig)
!          ! RN> bdfac factors seem missing for apar_p
!          apar_p = apargavg(ig+1)+apargavg(ig)
!          apar_m = aparnew(ig+1,it,ik)+aparnew(ig,it,ik) & 
!               -apar(ig+1,it,ik)-apar(ig,it,ik)

! !MAB, 6/5/2009:
! ! added the omprimfac source term arising with equilibrium flow shear  
! !CMR, 26/11/2010:
! ! Useful to understand that all source terms propto aky are specified here 
! ! using Tref=mref vtref^2. See 
! ! [note by BD and MK on "Microinstabilities in Axisymmetric Configurations"].
! ! This is converted to  the standard internal gs2 normalisation, 
! ! Tref=(1/2) mref vtref^2, by wunits, which contains a crucial factor 1/2.
! ! (Would be less confusing if always used same Tref!)
! !
! #ifdef LOWFLOW
!          source(ig) = anon(ie)*(vparterm(ig,isgn,iglo)*phi_m &
!               -spec(is)%zstm*vpac(ig,isgn,iglo) &
!               *((aj0(ig+1,iglo) + aj0(ig,iglo))*0.5*apar_m  &
!               + D_res(it,ik)*apar_p) &
!               -zi*wdfac(ig,isgn,iglo)*phi_p) &
!               + zi*(wstarfac(ig,isgn,iglo) &
!               + vpac(ig,isgn,iglo)*code_dt*wunits(ik)*ufac(ie,is) &
!               -2.0*omprimfac*vpac(ig,isgn,iglo)*code_dt*wunits(ik)*g_exb*itor_over_B(ig)/spec(is)%stm) &
!               *(phi_p - apar_p*spec(is)%stm*vpac(ig,isgn,iglo))
! #else
!          source(ig) = anon(ie)*(-2.0*vpar(ig,isgn,iglo)*phi_m &
!               -spec(is)%zstm*vpac(ig,isgn,iglo) &
!               *((aj0(ig+1,iglo) + aj0(ig,iglo))*0.5*apar_m  &
!               + D_res(it,ik)*apar_p) &
!               -zi*wdrift(ig,isgn,iglo)*phi_p) &
!               + zi*(wstar(ig,iv,iglo) &
!               + vpac(ig,isgn,iglo)*code_dt*wunits(ik)*ufac(ie,is) &
!               -2.0*omprimfac*vpac(ig,isgn,iglo)*code_dt*wunits(ik)*g_exb*itor_over_B(ig)/spec(is)%stm) &
!               *(phi_p - apar_p*spec(is)%stm*vpac(ig,isgn,iglo))
! #endif
!       end do
        
! #ifdef LOWFLOW
!       select case (istep)
!       case (0)
!          ! nothing
!       case (1)
!          do ig = -ntgrid, ntgrid-1
!             source(ig) = source(ig) + 0.5*code_dt*gexp_1(ig,isgn,iglo)
!          end do
!       case (2) 
!          do ig = -ntgrid, ntgrid-1
!             source(ig) = source(ig) + 0.5*code_dt*( &
!                  1.5*gexp_1(ig,isgn,iglo) - 0.5*gexp_2(ig,isgn,iglo))
!          end do
!       case default
!          do ig = -ntgrid, ntgrid-1
!             source(ig) = source(ig) + 0.5*code_dt*( &
!                  (23./12.)*gexp_1(ig,isgn,iglo) &
!                  - (4./3.)  *gexp_2(ig,isgn,iglo) &
!                  + (5./12.) *gexp_3(ig,isgn,iglo))
!          end do
!       end select
! #else
! ! add in nonlinear terms 
!       if (nonlin) then         
!          select case (istep)
!          case (0)
!             ! nothing
!          case (1)
!             do ig = -ntgrid, ntgrid-1
!                source(ig) = source(ig) + 0.5*code_dt*gexp_1(ig,isgn,iglo)
!             end do
!          case (2) 
!             do ig = -ntgrid, ntgrid-1
!                source(ig) = source(ig) + 0.5*code_dt*( &
!                     1.5*gexp_1(ig,isgn,iglo) - 0.5*gexp_2(ig,isgn,iglo))
!             end do
!          case default
!             do ig = -ntgrid, ntgrid-1
!                source(ig) = source(ig) + 0.5*code_dt*( &
!                       (23./12.)*gexp_1(ig,isgn,iglo) &
!                     - (4./3.)  *gexp_2(ig,isgn,iglo) &
!                     + (5./12.) *gexp_3(ig,isgn,iglo))
!             end do
!          end select
!       end if
! #endif

!     end subroutine set_source

!   end subroutine get_source_term

  ! MAB FLAG -- need to uncomment all the way down to reset_init
  subroutine getan (antot, antota, antotp)

    use dist_fn_arrays, only: aj0, aj1, gnew, kperp2
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use vpamu_grids, only: vperp2, vpa, integrate_species, nvgrid
    use run_parameters, only: beta, fphi, fapar, fbpar
    use prof, only: prof_entering, prof_leaving
    use gs2_layouts, only: g_lo, imu_idx

    implicit none

    complex, dimension (-ntgrid:,:,:), intent (out) :: antot, antota, antotp
    real, dimension (nspec) :: wgt

    integer :: iv, iglo, ig, imu

    call prof_entering ("getan", "dist_fn")

    if (fphi > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do iv = -nvgrid, nvgrid
             do ig=-ntgrid, ntgrid
                g0(ig,iv,iglo) = aj0(ig,iglo)*gnew(ig,iv,iglo)
             end do
          end do
       end do

       wgt = spec%z*spec%dens
       call integrate_species (g0, wgt, antot)

       if (afilter > epsilon(0.0)) antot = antot * exp(-afilter**4*kperp2**2/4.)
    else
       antot=0.
    end if

    if (fapar > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do iv = -nvgrid, nvgrid
             do ig=-ntgrid, ntgrid
                g0(ig,iv,iglo) = aj0(ig,iglo)*vpa(iv)*gnew(ig,iv,iglo)
             end do
          end do
       end do
       
       wgt = 2.0*beta*spec%z*spec%dens*sqrt(spec%temp/spec%mass)
       call integrate_species (g0, wgt, antota)
    else
       antota=0.
    end if

    if (fbpar > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          imu = imu_idx(g_lo,iglo)
          do iv = -nvgrid, nvgrid
             do ig=-ntgrid, ntgrid
                g0(ig,iv,iglo) = aj1(ig,iglo)*vperp2(ig,imu)*gnew(ig,iv,iglo)
             end do
          end do
       end do
       wgt = spec%temp*spec%dens
       call integrate_species (g0, wgt, antotp)
    else
       antotp=0.
    end if

    call prof_leaving ("getan", "dist_fn")
  end subroutine getan

  subroutine getmoms (ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)

    use dist_fn_arrays, only: aj0, aj1, gnew, g_adjust
    use gs2_layouts, only: is_idx, imu_idx, g_lo, ik_idx, it_idx
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use vpamu_grids, only: vpa, vperp2, integrate_moment, anon, energy, nvgrid
    use prof, only: prof_entering, prof_leaving
    use run_parameters, only: fphi, fbpar
    use fields_arrays, only: phinew, bparnew

    implicit none
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: density, &
         upar, tpar, tperp, ntot, qparflux, pperpj1, qpperpj1

    integer :: ik, it, iv, imu, is, iglo, ig

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
       do iv = -nvgrid, nvgrid
          g0(:,iv,iglo) = aj0(:,iglo)*gnew(:,iv,iglo)
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
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do ig = -ntgrid, ntgrid
          g0(ig,:,iglo) = vpa*g0(ig,:,iglo)
       end do
    end do
    call integrate_moment (g0, upar)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do ig = -ntgrid, ntgrid
          g0(ig,:,iglo) = 2.0*vpa*g0(ig,:,iglo)
       end do
    end do
    call integrate_moment (g0, tpar)
! tpar transiently stores ppar, nonadiabatic perturbed par pressure 
!      vpa normalised to: sqrt(2 T_s T_ref/m_s m_ref)
!  including factor 2 in g0 product ensures 
!     ppar normalised to: n_s T_s n_ref T_ref 
!                         ppar = tpar + density, and so:
    tpar = tpar - density
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       imu = imu_idx(g_lo,iglo)
       do iv = -nvgrid, nvgrid
          g0(:,iv,iglo) = vperp2(:,imu)*gnew(:,iv,iglo)*aj0(:,iglo)
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
       imu = imu_idx(g_lo,iglo)
       do iv = -nvgrid, nvgrid
          do ig = -ntgrid, ntgrid
             g0(ig,iv,iglo) = vpa(iv)*gnew(ig,iv,iglo)*aj0(ig,iglo)*energy(ig,iv,imu)
          end do
       end do
    end do 
    call integrate_moment (g0, qparflux)
   
! Now compute PPERPJ1, a modified p_perp which gives particle flux from Bpar
! NB PPERPJ1 is normalised to (n_s T_s/q_s)  n_ref T_ref/q_ref 
!    ie multiply by (n_s spec(is)%tz) n_ref T_ref/q_ref to get abs PPERPJ1
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       imu = imu_idx(g_lo,iglo)
       do iv = -nvgrid, nvgrid
          g0(:,iv,iglo) &
               = gnew(:,iv,iglo)*aj1(:,iglo)*2.0*vperp2(:,imu)*spec(is)%tz
       end do
    end do
    call integrate_moment (g0, pperpj1)

! Now compute QPPERPJ1, a modified p_perp*energy which gives heat flux from Bpar
! NB QPPERPJ1 is normalised to (n_s T_s^2/q_s)  n_ref  T_ref^2/q_ref
!    ie multiply by (n_s T_s spec(is)%tz) n_ref T_ref^2/q_ref to get abs QPPERPJ1
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       imu = imu_idx(g_lo,iglo)
       do iv = -nvgrid, nvgrid
          g0(:,iv,iglo) = g0(:,iv,iglo)*energy(:,iv,imu)
       end do
    end do
    call integrate_moment (g0, qpperpj1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now include the adiabatic part of <delta f>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! set g0 = <delta_f>/F_m, including the adiabatic term
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       imu = imu_idx(g_lo,iglo) ; is = is_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo) ; it = it_idx(g_lo,iglo)
       do iv = -nvgrid, nvgrid
          g0(:,iv,iglo) = aj0(:,iglo)*gnew(:,iv,iglo) - anon(:,iv,imu)*phinew(:,it,ik)*spec(is)%zt
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

!   subroutine getemoms (ntot, tperp)
!     use dist_fn_arrays, only: vperp2, aj0, gnew, g_adjust
!     use gs2_layouts, only: is_idx, ie_idx, g_lo, ik_idx, it_idx
!     use species, only: nspec, spec
!     use theta_grid, only: ntgrid
!     use le_grids, only: integrate_moment, anon
!     use prof, only: prof_entering, prof_leaving
!     use run_parameters, only: fphi, fbpar
!     use fields_arrays, only: phinew, bparnew

!     implicit none
!     complex, dimension (-ntgrid:,:,:,:), intent (out) :: tperp, ntot

!     integer :: ik, it, isgn, ie, is, iglo

! ! returns electron density and Tperp moment integrals to PE 0
!     call prof_entering ("getemoms", "dist_fn")
! !
! ! What are <delta_f> and g_wesson in the note below?
! ! g_wesson = <delta_f> + q phi/T    [ignore F_m factors for simplicity]
! !
! ! Electrostatically (for simplicity), g_adjust produces:
! !
! ! h = g_gs2 + q <phi> / T
! ! 
! ! then in the subsequent code they calculate for ntot:
! !
! ! ntot = integral[   J0 h - q phi / T  ]
! !
! ! so g_wesson == h.  What is odd in our notation is the LHS.  
! ! We typically indicate the perturbed distribution at fixed spatial position 
! ! by delta_f.  In the note below, they must mean delta_f = delta_f (R), so that 
! ! they are gyro-averaging to get the distribution at fixed spatial position.
! !
! ! In summary, DJA and CMR are calculating the moments at fixed spatial position
! ! rather than at fixed guiding centers, and using different notation than appears
! ! in most of our papers.
! !
! ! DJA+CMR: 17/1/06, use g_adjust routine to extract g_wesson
! !                   from gnew, phinew and bparnew.
! !           nb  <delta_f> = g_wesson J0 - q phi/T F_m  where <> = gyroaverage
! !           ie  <delta_f>/F_m = g_wesson J0 - q phi/T
! !
! ! use g0 as dist_fn dimensioned working space for all moments
! ! (avoid making more copies of gnew to save memory!)
! !
! ! set gnew = g_wesson, but return gnew to entry state prior to exit 
!     call g_adjust(gnew, phinew, bparnew, fphi, fbpar)

!     do iglo = g_lo%llim_proc, g_lo%ulim_proc       
!        ie = ie_idx(g_lo,iglo) ; is = is_idx(g_lo,iglo)
!        ik = ik_idx(g_lo,iglo) ; it = it_idx(g_lo,iglo)
!        do isgn = 1, 2
!           g0(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo) - anon(ie)*phinew(:,it,ik)*spec(is)%zt
!        end do
!     end do

! ! total perturbed density
!     call integrate_moment (g0, ntot)

! ! vperp**2 moment:
!     do iglo = g_lo%llim_proc, g_lo%ulim_proc       
!        ie = ie_idx(g_lo,iglo) ; is = is_idx(g_lo,iglo)
!        ik = ik_idx(g_lo,iglo) ; it = it_idx(g_lo,iglo)
!        do isgn = 1, 2
!           g0(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo)*vperp2(:,iglo) - anon(ie)*phinew(:,it,ik)*spec(is)%zt*vperp2(:,iglo)
!        end do
!     end do

! ! total perturbed perp pressure
!     call integrate_moment (g0, tperp)

! ! tperp transiently stores pperp, 
! !       pperp = tperp + density, and so:
!     tperp = tperp - ntot

!     do is=1,nspec
!        ntot(:,:,:,is)=ntot(:,:,:,is)*spec(is)%dens
!        tperp(:,:,:,is)=tperp(:,:,:,is)*spec(is)%temp
!     end do

! ! return gnew to its initial state, the variable evolved in GS2
!     call g_adjust(gnew,phinew,bparnew,-fphi,-fbpar)

!     call prof_leaving ("getemoms", "dist_fn")
!   end subroutine getemoms
  
  ! moment at not guiding center coordinate
  subroutine getmoms_notgc (dens, upar, tpar, tper, ntot, jpar)

    use dist_fn_arrays, only: aj0, aj1, gnew
    use gs2_layouts, only: g_lo, is_idx, ik_idx, it_idx, imu_idx
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use kt_grids, only: nakx => ntheta0, naky
    use vpamu_grids, only: vpa, vperp2, integrate_moment, nvgrid
    use mp, only: iproc
    use fields_arrays, only: phinew, bparnew

    implicit none

    complex, intent (out) :: &
         & dens(-ntgrid:,:,:,:), upar(-ntgrid:,:,:,:), &
         & tpar(-ntgrid:,:,:,:), tper(-ntgrid:,:,:,:)
    complex, intent (out), optional :: ntot(-ntgrid:,:,:,:)
    complex, intent (out), optional :: jpar(-ntgrid:,:,:)
    integer :: isgn, iglo, is

    real :: a, b, tpar2, tper2
    integer :: it, ik, ig, iv, imu

! returns moment integrals to PE 0

! not guiding center n_total
    if(present(ntot)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          imu = imu_idx(g_lo,iglo)

          do iv = -nvgrid, nvgrid
             g0(:,iv,iglo) = aj0(:,iglo) * gnew(:,iv,iglo)
          end do
          do iv = -nvgrid, nvgrid
             g0(:,iv,iglo) = g0(:,iv,iglo) + phinew(:,it,ik) &
                  & *(aj0(:,iglo)**2-1.0) * spec(is)%zt
          end do
          do iv = -nvgrid, nvgrid
             g0(:,iv,iglo) = g0(:,iv,iglo) &
                  & + 2.*vperp2(:,imu)*aj1(:,iglo)*aj0(:,iglo) &
                  & * bparnew(:,it,ik)
          end do
       end do
       call integrate_moment (g0, ntot, 1)
    endif

! not guiding center density
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do iv = -nvgrid, nvgrid
          g0(:,iv,iglo) = aj0(:,iglo)*gnew(:,iv,iglo)
       end do
    end do

    call integrate_moment (g0, dens, 1)

! not guiding center upar
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do iv = -nvgrid, nvgrid
          g0(:,iv,iglo) = aj0(:,iglo)*vpa(iv)*gnew(:,iv,iglo)
       end do
    end do

    call integrate_moment (g0, upar, 1)

! not guiding center tpar
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do iv = -nvgrid, nvgrid
          g0(:,iv,iglo) = 2.*aj0(:,iglo)*vpa(iv)**2*gnew(:,iv,iglo)
       end do
    end do

    call integrate_moment (g0, tpar, 1)
    
! not guiding center tperp
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do iv = -nvgrid, nvgrid
          imu = imu_idx(g_lo,iglo)
          g0(:,iv,iglo) = 2.*vperp2(:,imu)*aj1(:,iglo)*gnew(:,iv,iglo)
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

    use dist_fn_arrays, only: aj0, aj1, kperp2
    use species, only: nspec, spec, has_electron_species
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0, aky
    use vpamu_grids, only: anon, integrate_species, vperp2, nvgrid
    use gs2_layouts, only: g_lo, is_idx, imu_idx
    use run_parameters, only: tite

    implicit none

    integer :: iglo, iv
    integer :: ik, it, imu, is
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: tot
    real, dimension (nspec) :: wgt

    if (feqinit) return
    feqinit = .true.

    allocate (gamtot(-ntgrid:ntgrid,ntheta0,naky))
    allocate (gamtot1(-ntgrid:ntgrid,ntheta0,naky))
    allocate (gamtot2(-ntgrid:ntgrid,ntheta0,naky))
    if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then	
       allocate (gamtot3(-ntgrid:ntgrid,ntheta0,naky))
    endif
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       imu = imu_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do iv = -nvgrid, nvgrid
          g0(:,iv,iglo) = (1.0 - aj0(:,iglo)**2)*anon(:,iv,imu)
       end do
    end do
    wgt = spec%z*spec%z*spec%dens/spec%temp
    call integrate_species (g0, wgt, tot)
    gamtot = real(tot) + kperp2*poisfac
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       imu = imu_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do iv = -nvgrid, nvgrid
          g0(:,iv,iglo) = aj0(:,iglo)*aj1(:,iglo) &
               *2.0*vperp2(:,imu)*anon(:,iv,imu)
       end do
    end do
    wgt = spec%z*spec%dens
    call integrate_species (g0, wgt, tot)
    gamtot1 = real(tot)
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       imu = imu_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do iv = -nvgrid, nvgrid
          g0(:,iv,iglo) = aj1(:,iglo)**2*2.0*vperp2(:,imu)**2*anon(:,iv,imu)
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
    
    if (.not. allocated(fl_avg)) allocate (fl_avg(ntheta0, naky))
    fl_avg = 0.

    if (.not. has_electron_species(spec)) then
       if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
          
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

       end if

    end if

    if (fphi > epsilon(0.0)) then
       fieldeq = antot + bpar*gamtot1 - gamtot*phi 

       if (.not. has_electron_species(spec)) then
          do ik = 1, naky
             do it = 1, ntheta0
                fieldeq(:,it,ik) = fieldeq(:,it,ik) + fl_avg(it,ik)
             end do
          end do
       end if
    end if

    if (fapar > epsilon(0.0)) then
       fieldeqa = antota - kperp2*apar
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
       fieldeqp = fieldeqp + bpar
    end if

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
  
  ! Given initial distribution function this obtains consistent fields
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
    use dist_fn_arrays, only: aj0, kperp2

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
    use vpamu_grids, only: energy, vpa, vperp2, nvgrid
    use dist_fn_arrays, only: gnew, aj0, aj1
    use gs2_layouts, only: g_lo, imu_idx, is_idx, it_idx, ik_idx
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
    integer :: it, ik, is, iv, ig, imu, iglo

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
       do iv=-nvgrid,nvgrid
          g0(:,iv,:) = gnew(:,iv,:)*aj0
       end do
       call get_flux (phi, pflux, dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          imu = imu_idx(g_lo,iglo)
          g0(:,:,iglo) = g0(:,:,iglo)*energy(:,:,imu)
       end do
       call get_flux (phi, qflux(:,:,:,1), dnorm)

       do iv = -nvgrid, nvgrid
          g0(:,iv,:) = gnew(:,iv,:)*2.*vpa(iv)**2*aj0
       end do
       call get_flux (phi, qflux(:,:,:,2), dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          imu = imu_idx(g_lo,iglo)
          do iv = -nvgrid, nvgrid
             g0(:,iv,iglo) = gnew(:,iv,iglo)*vperp2(:,imu)*aj0(:,iglo)
          end do
       end do
       call get_flux (phi, qflux(:,:,:,3), dnorm)

       do iv = -nvgrid, nvgrid
          do ig = -ntgrid, ntgrid
             g0(ig,iv,:) = gnew(ig,iv,:)*aj0(ig,:)*vpa(iv)*Rplot(ig)*sqrt(1.0-Bpol(ig)**2/bmag(ig)**2)
          end do
       end do
       call get_flux (phi, vflux_par, dnorm)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          imu = imu_idx(g_lo,iglo)
          do iv = -nvgrid, nvgrid
             g0(:,iv,iglo) = -zi*aky(ik)*gnew(:,iv,iglo)*aj1(:,iglo) &
                  *rhoc*(gds21+theta0(it,ik)*gds22)*vperp2(:,imu)*spec(is)%smz/(qval*shat*bmag**2)
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
          do iv = -nvgrid, nvgrid
             g0(:,iv,iglo) &
                  = -gnew(:,iv,iglo)*aj0(:,iglo)*spec(is)%stm*vpa(iv)
          end do
       end do
       call get_flux (apar, pmflux, dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          imu = imu_idx(g_lo,iglo)
          g0(:,:,iglo) = g0(:,:,iglo)*energy(:,:,imu)
       end do
       call get_flux (apar, qmflux(:,:,:,1), dnorm)
       
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do iv = -nvgrid, nvgrid
             g0(:,iv,iglo) &
                  = -gnew(:,iv,iglo)*aj0(:,iglo)*spec(is)%stm*vpa(iv) &
                  *2.*vpa(iv)**2
          end do
       end do
       call get_flux (apar, qmflux(:,:,:,2), dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          imu = imu_idx(g_lo,iglo)
          do iv = -nvgrid, nvgrid
             g0(:,iv,iglo) &
                  = -gnew(:,iv,iglo)*aj0(:,iglo)*spec(is)%stm*vpa(iv) &
                  *vperp2(:,imu)
          end do
       end do
       call get_flux (apar, qmflux(:,:,:,3), dnorm)
       
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do iv = -nvgrid, nvgrid
             g0(:,iv,iglo) &
                  = -gnew(:,iv,iglo)*aj0(:,iglo)*spec(is)%stm &
                  *vpa(iv)**2
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
          imu = imu_idx(g_lo,iglo)
          do iv = -nvgrid, nvgrid
             g0(:,iv,iglo) &
                  = gnew(:,iv,iglo)*aj1(:,iglo)*2.0*vperp2(:,imu)*spec(is)%tz
          end do
       end do
       call get_flux (bpar, pbflux, dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          imu = imu_idx(g_lo,iglo)
          g0(:,:,iglo) = g0(:,:,iglo)*energy(:,:,imu)
       end do
       call get_flux (bpar, qbflux(:,:,:,1), dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          imu = imu_idx(g_lo,iglo)
          do iv = -nvgrid, nvgrid
             g0(:,iv,iglo) &
                  = gnew(:,iv,iglo)*aj1(:,iglo)*2.0*vperp2(:,imu)*spec(is)%tz &
                    *2.*vpa(iv)**2
          end do
       end do
       call get_flux (bpar, qbflux(:,:,:,2), dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          imu = imu_idx(g_lo,iglo)
          do iv = -nvgrid, nvgrid
             g0(:,iv,iglo) &
                  = gnew(:,iv,iglo)*aj1(:,iglo)*2.0*vperp2(:,imu)*spec(is)%tz &
                    *vperp2(:,imu)
          end do
       end do
       call get_flux (bpar, qbflux(:,:,:,3), dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          imu = imu_idx(g_lo,iglo)
          do iv = -nvgrid, nvgrid
             g0(:,iv,iglo) &
                  = gnew(:,iv,iglo)*aj1(:,iglo)*2.0*vperp2(:,imu) &
                  *spec(is)%tz*vpa(iv)
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

    use theta_grid, only: ntgrid, delthet, grho
    use kt_grids, only: ntheta0, aky, naky
    use vpamu_grids, only: integrate_moment
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

  subroutine eexchange (phi, exchange)

    use mp, only: proc0
    use constants, only: zi
    use gs2_layouts, only: g_lo, imu_idx, it_idx, ik_idx, is_idx
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: gnew, aj0
    use theta_grid, only: ntgrid, gradpar, delthet, bmag, jacob, thet_imp
    use kt_grids, only: ntheta0, naky
    use vpamu_grids, only: integrate_moment, nvgrid, vpac, vpa_imp
    use run_parameters, only: woutunits, fphi
    use species, only: spec, nspec
    use nonlinear_terms, only: nonlin
    use centering, only: get_cell_value

    implicit none

    complex, dimension (-ntgrid:,:,:), intent (in) :: phi
    real, dimension (:,:,:), intent (out) :: exchange

    integer :: ig, imu, it, ik, is, iglo, iv
    real :: wgt
    real, dimension (:), allocatable :: gradparc
    real, dimension (:,:,:), allocatable :: dnorm
    complex, dimension (:,:,:,:), allocatable :: total

    allocate (gradparc(-ntgrid:ntgrid)) ; gradparc = 0.0
    allocate (dnorm(-ntgrid:ntgrid, ntheta0, naky)) ; dnorm = 0.0
    allocate (total(-ntgrid:ntgrid, ntheta0, naky, nspec)) ; total = 0.0

    if (proc0) exchange = 0.0

    do ik = 1, naky
       do it = 1, ntheta0
          dnorm(:ntgrid-1,it,ik) = delthet(:ntgrid-1)*(jacob(-ntgrid+1:)+jacob(:ntgrid-1))*0.5
       end do
    end do

    call get_cell_value (thet_imp, gradpar, gradparc, -ntgrid)
    if (fphi > epsilon(0.0)) then
       g0 = 0.
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          imu = imu_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          if (nonlin .and. it==0 .and. ik==1) cycle
          do iv = -nvgrid, nvgrid
             ! get v_magnetic piece of g0 at grid points instead of cell centers
             do ig = -ntgrid, ntgrid
                g0(ig,iv,iglo) = aj0(ig,iglo)*(zi*wdrift_func(ig,iv,imu,it,ik)/code_dt &
                     * gnew(ig,iv,iglo)*spec(is)%tz)
             end do
          end do
          ! get cell centered value in theta and vpa
          call get_cell_value (thet_imp, vpa_imp, g0(:,:,iglo), g0(:,:,iglo), -ntgrid, -nvgrid)
          ! note that gnew below is not cell centered in vpa, which is should be!
          ! MAB FLAG
          do iv = -nvgrid, nvgrid
             ! get v_magnetic piece of g0 at cell centers and add in vpar piece at cell centers
             do ig = -ntgrid, -1
                g0(ig,iv,iglo) = g0(ig,iv,iglo) &
                     + vpac(iv,1)*gradparc(ig)/delthet(ig) &
                     * (gnew(ig+1,iv,iglo)-gnew(ig,iv,iglo))*spec(is)%stm
             end do
             do ig = 0, ntgrid-1
                g0(ig,iv,iglo) = g0(ig,iv,iglo) &
                     + vpac(iv,2)*gradparc(ig)/delthet(ig) &
                     * (gnew(ig+1,iv,iglo)-gnew(ig,iv,iglo))*spec(is)%stm
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
                   exchange(it,ik,is) = sum(real(total(:,it,ik,is)*conjg(phi(:,it,ik))) &
                        *dnorm(:,it,ik))/wgt
                end do
             end do
          end do
       end if

    end if

    deallocate (dnorm, total, gradparc)

  end subroutine eexchange

!   subroutine lf_flux (phi, vflx0, vflx1)

!     use species, only: spec, nspec
!     use theta_grid, only: ntgrid, bmag, gradpar, grho, delthet
!     use theta_grid, only: qval, shat, gds21, gds22, drhodpsi, IoB
!     use kt_grids, only: naky, ntheta0, akx, aky
!     use dist_fn_arrays, only: gnew, aj0, vpac, vpa, aj1, vperp2
!     use gs2_layouts, only: g_lo, ie_idx, is_idx, it_idx, ik_idx
!     use mp, only: proc0
!     use run_parameters, only: woutunits, fphi, fapar, fbpar, rhostar
!     use constants, only: zi
!     use geometry, only: rhoc
!     implicit none
!     complex, dimension (-ntgrid:,:,:), intent (in) :: phi
!     real, dimension (:,:,:), intent (out) :: vflx0, vflx1
!     real, dimension (:,:), allocatable :: dum
!     real, dimension (:,:,:), allocatable :: dnorm
!     complex, dimension (:,:,:), allocatable :: dphi
!     integer :: it, ik, is, isgn, ig
!     integer :: iglo

!     allocate (dnorm (-ntgrid:ntgrid,ntheta0,naky))
!     allocate (dum (-ntgrid:ntgrid,nspec))
!     allocate (dphi (-ntgrid:ntgrid,ntheta0,naky))

!     if (proc0) then
!        vflx0 = 0.0 ; vflx1 = 0.0 ; dum = 0.0
!     end if

!     do ik = 1, naky
!        do it = 1, ntheta0
!           dnorm(:,it,ik) = delthet/bmag/gradpar*woutunits(ik)
!        end do
!     end do

!     do ig = -ntgrid, ntgrid-1
!        dphi(ig,:,:) = (phi(ig+1,:,:)-phi(ig,:,:))/delthet(ig)
!     end do
!     ! not sure if this is the right way to handle ntgrid point -- MAB
!     dphi(ntgrid,:,:) = (phi(ntgrid,:,:)-phi(ntgrid-1,:,:))/delthet(-ntgrid)

!     if (fphi > epsilon(0.0)) then
!        ! this is the second term in Pi_0^{tb} in toroidal_flux.pdf notes
!        do iglo = g_lo%llim_proc, g_lo%ulim_proc
!           do isgn = 1, 2
!              g0(:,isgn,iglo) = -zi*gnew(:,isgn,iglo)*aj0(:,iglo)*vpa(:,isgn,iglo) &
!                   *drhodpsi*IoB**2*gradpar*rhostar
!           end do
!        end do
!        call get_lfflux (dphi, vflx0, dnorm)

!        ! this is the bracketed part of the first term in Pi_0^{tb} in toroidal_flux.pdf notes
!        do iglo = g_lo%llim_proc, g_lo%ulim_proc
!           do isgn = 1, 2
!              g0(:,isgn,iglo) = 0.5*gnew(:,isgn,iglo)*aj0(:,iglo)*vpa(:,isgn,iglo)**2 &
!                   *drhodpsi*IoB**2*rhostar
!           end do
!        end do
!        call get_flux (phi, vflx1, dnorm)

! !        do isgn = 1, 2
! !           do ig = -ntgrid, ntgrid
! !              if (allocated(rmajor_geo)) then
! !                 g0(ig,isgn,:) = gnew(ig,isgn,:)*aj0(ig,:)*vpac(ig,isgn,:)*rmajor_geo(ig)*sqrt(1.0-bpol_geo(ig)**2/bmag(ig)**2)
! !              else
! !                 g0(ig,isgn,:) = gnew(ig,isgn,:)*aj0(ig,:)*vpac(ig,isgn,:)
! !              end if
! !           end do
! !        end do
! !        call get_flux (phi, vflux_par, dnorm)
!     else
!        vflx0 = 0. ; vflx1 = 0.
!     end if

!     deallocate (dnorm,dum,dphi)
!   end subroutine lf_flux

!   subroutine get_lfflux (fld, flx, dnorm)
!     use theta_grid, only: ntgrid
!     use kt_grids, only: ntheta0, aky, naky
!     use le_grids, only: integrate_moment
!     use theta_grid, only: grho
!     use species, only: nspec
!     use mp, only: proc0
!     implicit none
!     complex, dimension (-ntgrid:,:,:), intent (in) :: fld
!     real, dimension (:,:,:), intent (in out) :: flx
!     real, dimension (-ntgrid:,:,:) :: dnorm
!     complex, dimension (:,:,:,:), allocatable :: total
!     real :: wgt
!     integer :: ik, it, is, ig

!     allocate (total(-ntgrid:ntgrid,ntheta0,naky,nspec))
!     call integrate_moment (g0, total)

!     if (proc0) then
!        do is = 1, nspec
!           do ik = 1, naky
!              do it = 1, ntheta0
!                 wgt = sum(dnorm(:,it,ik)*grho)
!                 flx(it,ik,is) = sum(aimag(total(:,it,ik,is)*conjg(fld(:,it,ik))) &
!                      *dnorm(:,it,ik))/wgt
!              end do
!           end do
!        end do

!        flx = flx*0.5

!     end if

!     deallocate (total)

!   end subroutine get_lfflux

! !=============================================================================
!   subroutine get_heat (h, hk, phi, apar, bpar, phinew, aparnew, bparnew)
!     use mp, only: proc0, iproc
!     use constants, only: pi, zi
!     use kt_grids, only: ntheta0, naky, aky, akx
! #ifdef LOWFLOW
!     use dist_fn_arrays, only: hneoc
! #endif
!     use dist_fn_arrays, only: vpa, vpac, aj0, aj1, vperp2, g, gnew, kperp2, g_adjust
!     use gs2_heating, only: heating_diagnostics
!     use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx, ie_idx
!     use le_grids, only: integrate_moment
!     use species, only: spec, nspec,has_electron_species
!     use theta_grid, only: jacob, delthet, ntgrid
!     use run_parameters, only: fphi, fapar, fbpar, tunits, beta, tite
!     use gs2_time, only: code_dt
!     use nonlinear_terms, only: nonlin
!     use antenna, only: antenna_apar, a_ext_data
!     use hyper, only: D_v, D_eta, nexp, hypervisc_filter

!     implicit none
!     type (heating_diagnostics) :: h
!     type (heating_diagnostics), dimension(:,:) :: hk
! !    complex, dimension (-ntgrid:,:,:), pointer :: hh, hnew
!     complex, dimension (-ntgrid:,:,:) :: phi, apar, bpar, phinew, aparnew, bparnew
!     complex, dimension(:,:,:,:), allocatable :: tot
!     complex, dimension(:,:,:), allocatable :: bpardot, apardot, phidot, j_ext
!     complex :: chi, havg
!     complex :: chidot, j0phiavg, j1bparavg, j0aparavg
! !    complex :: pstar, pstardot, gdot
!     complex :: phi_m, apar_m, bpar_m, hdot
!     complex :: phi_avg, bpar_avg, bperp_m, bperp_avg
!     complex :: de, denew
!     complex :: dgdt_hypervisc
!     real, dimension (:), allocatable :: wgt
!     real :: fac2, dtinv, akperp4
!     integer :: isgn, iglo, ig, is, ik, it, ie

!     g0(ntgrid,:,:) = 0.

! ! ==========================================================================
! ! Ion/Electron heating------------------------------------------------------
! ! ==========================================================================

!     allocate ( phidot(-ntgrid:ntgrid, ntheta0, naky))
!     allocate (apardot(-ntgrid:ntgrid, ntheta0, naky))
!     allocate (bpardot(-ntgrid:ntgrid, ntheta0, naky))

!     call dot ( phi,  phinew,  phidot, fphi)
!     call dot (apar, aparnew, apardot, fapar)
!     call dot (bpar, bparnew, bpardot, fbpar)

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!! Next two calls make the variables g, gnew = h, hnew 
! !!! until the end of this procedure!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     call g_adjust (g,    phi,    bpar,    fphi, fbpar)
!     call g_adjust (gnew, phinew, bparnew, fphi, fbpar)

!     do iglo=g_lo%llim_proc, g_lo%ulim_proc
!        is = is_idx(g_lo, iglo)
!        it = it_idx(g_lo, iglo)
!        ik = ik_idx(g_lo, iglo)
!        if (nonlin .and. it == 1 .and. ik == 1) cycle
!        dtinv = 1./(code_dt*tunits(ik))
!        do isgn=1,2
          
!           do ig=-ntgrid, ntgrid-1
             
!              chidot = aj0(ig,iglo)*(phidot(ig,it,ik) &
!                   - vpac(ig,isgn,iglo) * spec(is)%stm * apardot(ig,it,ik)) &
!                   + aj1(ig,iglo)*2.0*vperp2(ig,iglo)*bpardot(ig,it,ik)*spec(is)%tz
             
!              hdot = fdot (g   (ig  ,isgn,iglo), &
!                           g   (ig+1,isgn,iglo), &
!                           gnew(ig  ,isgn,iglo), &
!                           gnew(ig+1,isgn,iglo), dtinv)
             
!              havg = favg (g   (ig  ,isgn,iglo), &
!                           g   (ig+1,isgn,iglo), &
!                           gnew(ig  ,isgn,iglo), &
!                           gnew(ig+1,isgn,iglo))
             
! ! First term on RHS and LHS of Eq B-10 of H1:

!              g0(ig,isgn,iglo) = spec(is)%dens*conjg(havg)* &
!                   (chidot*spec(is)%z-hdot*spec(is)%temp)

!           end do
!        end do
!     end do

!     deallocate (phidot, apardot, bpardot)

!     allocate (tot(-ntgrid:ntgrid, ntheta0, naky, nspec))

!     call integrate_moment (g0, tot)

!     if (proc0) then
!        allocate (wgt(-ntgrid:ntgrid))
!        wgt = 0.
!        do ig=-ntgrid,ntgrid-1
! !          wgt(ig) = delthet(ig)*jacob(ig)
! ! delthet is cell-centered, but jacob is given on grid
!           wgt(ig) = delthet(ig)*(jacob(ig)+jacob(ig+1))*0.5
!        end do
!        wgt = wgt/sum(wgt)         

!        do is = 1, nspec
!           do ik = 1, naky
!              fac2 = 0.5
!              if (aky(ik) < epsilon(0.0)) fac2 = 1.0
!              do it = 1, ntheta0
!                 if (nonlin .and. it == 1 .and. ik == 1) cycle
!                 do ig = -ntgrid, ntgrid-1
!                     hk(it,ik) % heating(is) = hk(it,ik) % heating(is) &
!                         + real(tot(ig,it,ik,is))*wgt(ig)*fac2 
!                 end do
!                 h % heating(is) = h % heating(is) + hk(it,ik) % heating(is)
!              end do
!           end do
!        end do
!     end if

! ! ==========================================================================
! ! Antenna Power and B-field contribution to E and E_dot---------------------
! ! ==========================================================================
!     if (proc0) then
!        allocate (j_ext(-ntgrid:ntgrid, ntheta0, naky)) ; j_ext = 0.
!        call antenna_apar (kperp2, j_ext)       
       
!        if (beta > epsilon(0.)) then
!           do ik=1,naky
!              fac2 = 0.5
!              if (aky(ik) < epsilon(0.0)) fac2 = 1.0
!              dtinv = 1./(code_dt*tunits(ik))
!              do it = 1,ntheta0
                
!                 if (nonlin .and. it == 1 .and. ik == 1) cycle
                
!                 do ig=-ntgrid, ntgrid-1
                   
!                    !GGH Time and space averaged estimate of d/dt(apar)
!                    apar_m = fdot (apar   (ig  ,it,ik), &
!                                   apar   (ig+1,it,ik), &
!                                   aparnew(ig  ,it,ik), &
!                                   aparnew(ig+1,it,ik), dtinv)*fapar

! ! J_ext.E when driving antenna only includes A_parallel:

!                    hk(it,ik) % antenna = hk(it, ik) % antenna + real(conjg(j_ext(ig,it,ik))*apar_m)*wgt(ig)*fac2

!                    !GGH Time and space averaged estimate of d/dt(bperp)
!                    bperp_m = fdot (apar   (ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
!                                    apar   (ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik)), &
!                                    aparnew(ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
!                                    aparnew(ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik)), dtinv) * fapar

!                    !GGH Time and space averaged estimate of d/dt(bpar)
!                    bpar_m = fdot (bpar   (ig  ,it,ik), &
!                                   bpar   (ig+1,it,ik), &
!                                   bparnew(ig  ,it,ik), &
!                                   bparnew(ig+1,it,ik), dtinv)*fbpar

!                    !GGH Time and space averaged estimate of bperp
!                    bperp_avg = favg (apar   (ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
!                                      apar   (ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik)), &
!                                      aparnew(ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
!                                      aparnew(ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik))) * fapar

!                    !GGH Time and space averaged estimate of bpar
!                    bpar_avg = favg (bpar   (ig  ,it,ik), &
!                                     bpar   (ig+1,it,ik), &
!                                     bparnew(ig  ,it,ik), &
!                                     bparnew(ig+1,it,ik)) * fbpar

! ! 1/2 * d/dt B**2
! !! GGH: Bug fixed on 2/06; error was in relative weight of B_par**2 and B_perp**2   
!                    hk(it,ik) % energy_dot = hk(it,ik) % energy_dot + &
!                         real(0.25 * conjg(bperp_m)*bperp_avg + conjg(bpar_m)*bpar_avg) &
!                         * wgt(ig)*fac2*(2.0/beta)

! ! B**2/2
! !! GGH: Bug fixed on 2/06; error was in relative weight of B_par**2 and B_perp**2   
!                    hk(it,ik) % energy = hk(it,ik) % energy &
!                         + 0.5*real((0.25*conjg(bperp_avg)*bperp_avg + conjg(bpar_avg)*bpar_avg)) &
!                         * wgt(ig)*fac2*(2.0/beta)

!                    !Eapar = int k_perp^2 A_par^2/(8 pi)                   
!                    hk(it,ik) % eapar = hk(it,ik) % eapar &
!                         + 0.5*real(0.25*conjg(bperp_avg)*bperp_avg) * wgt(ig)*fac2*(2.0/beta)

!                    !Ebpar = int B_par^2/(8 pi)
!                    hk(it,ik) % ebpar = hk(it,ik) % ebpar &
!                         + 0.5*real(conjg(bpar_avg)*bpar_avg) * wgt(ig)*fac2*(2.0/beta)
                   
!                 end do
!                 h % antenna = h % antenna + hk(it, ik) % antenna
!                 h % eapar = h % eapar + hk(it, ik) % eapar
!                 h % ebpar = h % ebpar + hk(it, ik) % ebpar
!              end do
!           end do
!        else
!           hk % antenna = 0.
!           h  % antenna = 0.
!           hk % energy_dot = 0.
!           hk % energy = 0.
!           hk % eapar = 0.
!           h  % eapar = 0.
!           hk % ebpar = 0.
!           h  % ebpar = 0.
!        end if
!        deallocate (j_ext)
!     end if

! ! ==========================================================================
! ! Finish E_dot--------------------------------------------------------------
! ! ==========================================================================

! !GGH Include response of Boltzmann species for single-species runs

!     if (.not. has_electron_species(spec)) then
!        if (proc0) then
!           !NOTE: It is assumed here that n0i=n0e and zi=-ze
!           do ik=1,naky
!              fac2 = 0.5
!              if (aky(ik) < epsilon(0.0)) fac2 = 1.0
!              dtinv = 1./(code_dt*tunits(ik))
!              do it = 1,ntheta0

!                 if (nonlin .and. it == 1 .and. ik == 1) cycle

!                 do ig=-ntgrid, ntgrid-1

!                    phi_avg = favg (phi   (ig  ,it,ik), &
!                         phi   (ig+1,it,ik), &
!                         phinew(ig  ,it,ik), &
!                         phinew(ig+1,it,ik))

!                    phi_m   = fdot (phi   (ig  ,it,ik), &
!                         phi   (ig+1,it,ik), &
!                         phinew(ig  ,it,ik), &
!                         phinew(ig+1,it,ik), dtinv)

!                    !NOTE: Adiabatic (Boltzmann) species has temperature
!                    !       T = spec(1)%temp/tite
!                    hk(it,ik) % energy_dot = hk(it,ik) % energy_dot + &
!                         fphi * real(conjg(phi_avg)*phi_m) &
!                         * spec(1)%dens * spec(1)%z * spec(1)%z * (tite/spec(1)%temp) &
!                         * wgt(ig)*fac2

!                 end do
!              end do
!           end do
!        endif
!     endif !END Correction to E_dot for single species runs---------------------
 
! !GGH New E_dot calc
!     do iglo=g_lo%llim_proc, g_lo%ulim_proc
!        is = is_idx(g_lo, iglo)
!        it = it_idx(g_lo, iglo)
!        ik = ik_idx(g_lo, iglo)
!        if (nonlin .and. it == 1 .and. ik == 1) cycle
!        dtinv = 1./(code_dt*tunits(ik))
!        do isgn=1,2

!           do ig=-ntgrid, ntgrid-1
             
!              !Calculate old fluctuating energy de
!              havg = favg_x (g(ig  ,isgn,iglo), &
!                             g(ig+1,isgn,iglo))

!              j0phiavg = favg_x (aj0(ig  ,iglo) * phi(ig  ,it,ik), &
!                                 aj0(ig+1,iglo) * phi(ig+1,it,ik)) * fphi * spec(is)%zt

!              phi_avg = favg_x (phi(ig  ,it,ik), &
!                                phi(ig+1,it,ik)) * fphi * spec(is)%zt

!              de=0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg &
!                   + conjg(phi_avg)*phi_avg &
!                   - conjg(j0phiavg)*havg &
!                   - conjg(havg)*j0phiavg)

!             !Calculate new fluctuating energy denew
!              havg = favg_x (gnew(ig  ,isgn,iglo), &
!                             gnew(ig+1,isgn,iglo))

!              j0phiavg = favg_x (aj0(ig  ,iglo) * phinew(ig  ,it,ik), &
!                                 aj0(ig+1,iglo) * phinew(ig+1,it,ik)) * fphi * spec(is)%zt

!              phi_avg = favg_x (phinew(ig  ,it,ik), &
!                                phinew(ig+1,it,ik)) * fphi * spec(is)%zt

!              denew=0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg &
!                   + conjg(phi_avg)*phi_avg &
!                   - conjg(j0phiavg)*havg &
!                   - conjg(havg)*j0phiavg) 

!              !Set g0 as the change of energy (denew-de)/dt
!              g0(ig,isgn,iglo) = fdot_t(de,denew,dtinv)

!           end do
!        end do
!     end do
!     !GGH -END new e_dot calc

!    call integrate_moment (g0, tot)

!     if (proc0) then
!        do ik = 1, naky
!           fac2 = 0.5
!           if (aky(ik) < epsilon(0.0)) fac2 = 1.0
!           do it = 1, ntheta0
!              if (nonlin .and. it == 1 .and. ik == 1) cycle
!              do is = 1, nspec
!                 do ig = -ntgrid, ntgrid-1
!                    hk(it,ik) % energy_dot = hk(it,ik) % energy_dot  &
!                         + real(tot(ig,it,ik,is))*wgt(ig)*fac2
!                 end do
!              end do
!              h % energy_dot = h % energy_dot + hk(it,ik) % energy_dot
!           end do
!        end do
!     end if


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! ==========================================================================
! ! Gradient Contributions to Heating-----------------------------------------
! ! ==========================================================================

!     do iglo=g_lo%llim_proc, g_lo%ulim_proc
!        is = is_idx(g_lo, iglo)
!        it = it_idx(g_lo, iglo)
!        ik = ik_idx(g_lo, iglo)
!        ie = ie_idx(g_lo, iglo)
!        if (nonlin .and. it == 1 .and. ik == 1) cycle

!        do isgn=1,2
!           do ig=-ntgrid, ntgrid-1

!              chi = favg (aj0(ig  ,iglo)*phi   (ig  ,it,ik),  &
!                          aj0(ig+1,iglo)*phi   (ig+1,it,ik),  &
!                          aj0(ig  ,iglo)*phinew(ig  ,it,ik),  &
!                          aj0(ig+1,iglo)*phinew(ig+1,it,ik)) &
!                          *fphi

! !!GGH Bug fix: The apar part should be subtracted (because chi= phi - v|| A|| + B||)
!              chi = chi - &
!                   favg (aj0(ig  ,iglo)*apar   (ig  ,it,ik)*vpac(ig  ,isgn,iglo),  &
!                         aj0(ig+1,iglo)*apar   (ig+1,it,ik)*vpac(ig+1,isgn,iglo),  &
!                         aj0(ig  ,iglo)*aparnew(ig  ,it,ik)*vpac(ig  ,isgn,iglo),  &
!                         aj0(ig+1,iglo)*aparnew(ig+1,it,ik)*vpac(ig+1,isgn,iglo)) &
!                         *spec(is)%stm*fapar
                
!              chi = chi + &
!                   favg (aj1(ig  ,iglo)*2.0*bpar   (ig  ,it,ik)*vperp2(ig  ,iglo),  &
!                         aj1(ig+1,iglo)*2.0*bpar   (ig+1,it,ik)*vperp2(ig+1,iglo),  &
!                         aj1(ig  ,iglo)*2.0*bparnew(ig  ,it,ik)*vperp2(ig  ,iglo),  &
!                         aj1(ig+1,iglo)*2.0*bparnew(ig+1,it,ik)*vperp2(ig+1,iglo)) &
!                         *spec(is)%tz*fbpar

!              havg = favg (g   (ig  ,isgn,iglo), &
!                           g   (ig+1,isgn,iglo), &
!                           gnew(ig  ,isgn,iglo), &
!                           gnew(ig+1,isgn,iglo))

! #ifdef LOWFLOW
!              g0(ig,isgn,iglo) = zi * wstar(ig,iv,iglo)*hneoc(ig,isgn,iglo)/code_dt * conjg(havg)*chi &
! #else
!              g0(ig,isgn,iglo) = zi * wstar(ig,iv,iglo)/code_dt * conjg(havg)*chi &
! #endif
!                   * spec(is)%dens
            
!           end do
!        end do
!     end do

!     call integrate_moment (g0, tot)

!     if (proc0) then
!        do is = 1, nspec
!           do ik = 1, naky
!              fac2 = 0.5
!              if (aky(ik) < epsilon(0.0)) fac2 = 1.0
!              do it = 1, ntheta0
!                 if (nonlin .and. it == 1 .and. ik == 1) cycle
!                 do ig = -ntgrid, ntgrid-1
!                    hk(it,ik) % gradients(is) = hk(it,ik) % gradients(is) &
!                         + real(tot(ig,it,ik,is))*wgt(ig)*fac2
!                 end do
!                 h % gradients(is) = h % gradients(is) + hk(it,ik) % gradients(is)
!              end do
!           end do
!        end do
!     end if
! ! ==========================================================================
! ! Hyperviscosity------------------------------------------------------------
! ! ==========================================================================

!     if (D_v > epsilon(0.)) then

!        do iglo=g_lo%llim_proc, g_lo%ulim_proc
!           is = is_idx(g_lo, iglo)
!           it = it_idx(g_lo, iglo)
!           ik = ik_idx(g_lo, iglo)
!           if (nonlin .and. it == 1 .and. ik == 1) cycle
! !          akperp4 = (aky(ik)**2 + akx(it)**2)**nexp
!            do isgn=1,2
!              do ig=-ntgrid, ntgrid-1
                
!                 havg = favg (g   (ig  ,isgn,iglo), &
!                              g   (ig+1,isgn,iglo), &
!                              gnew(ig  ,isgn,iglo), &
!                              gnew(ig+1,isgn,iglo)) 

!                 j0phiavg = favg (aj0(ig  ,iglo) * phi(ig  ,it,ik), &
!                                  aj0(ig+1,iglo) * phi(ig+1,it,ik), &
!                                  aj0(ig  ,iglo) * phinew(ig  ,it,ik), &
!                                  aj0(ig+1,iglo) * phinew(ig+1,it,ik)) * fphi * spec(is)%zt

!                 j1bparavg= favg (aj1(ig  ,iglo)*2.0*bpar   (ig  ,it,ik)*vperp2(ig  ,iglo),  &
!                                  aj1(ig+1,iglo)*2.0*bpar   (ig+1,it,ik)*vperp2(ig+1,iglo),  &
!                                  aj1(ig  ,iglo)*2.0*bparnew(ig  ,it,ik)*vperp2(ig  ,iglo),  &
!                                  aj1(ig+1,iglo)*2.0*bparnew(ig+1,it,ik)*vperp2(ig+1,iglo)) &
!                                  *fbpar

!                 dgdt_hypervisc = 0.5*((1.0-1./hypervisc_filter(ig,it,ik))*gnew(ig,isgn,iglo) &
!                      + (1.0-1./hypervisc_filter(ig+1,it,ik))*gnew(ig+1,isgn,iglo))/code_dt

! !Set g0 for hyperviscous heating
! !                g0(ig,isgn,iglo) = spec(is)%dens*spec(is)%temp*D_v*akperp4* &
! !                     ( conjg(havg)*havg - conjg(havg)*j0phiavg - &
! !                     conjg(havg)*j1bparavg)
!                 g0(ig,isgn,iglo) = spec(is)%dens*spec(is)%temp*conjg(havg)*dgdt_hypervisc

!              end do
!           end do
!        end do

!        call integrate_moment (g0, tot)
!        if (proc0) then
!           do ik = 1, naky
!              fac2 = 0.5
!              if (aky(ik) < epsilon(0.0)) fac2 = 1.0
!              do it = 1, ntheta0
!                 if (nonlin .and. it == 1 .and. ik == 1) cycle
!                 do is = 1, nspec
!                    do ig = -ntgrid, ntgrid-1
!                       hk(it,ik) % hypervisc(is) = hk(it,ik) % hypervisc(is) &
!                            + real(tot(ig,it,ik,is))*wgt(ig)*fac2
!                    end do
!                    h % hypervisc(is) = h % hypervisc(is) + hk(it,ik) % hypervisc(is)
!                 end do
!              end do
!           end do
!        end if

!     end if !End Hyperviscous Heating Calculation


! ! ==========================================================================
! ! Hyperresistivity------------------------------------------------------------
! ! ==========================================================================
 
!     if (D_eta > epsilon(0.)) then

!        do iglo=g_lo%llim_proc, g_lo%ulim_proc
!           is = is_idx(g_lo, iglo)
!           it = it_idx(g_lo, iglo)
!           ik = ik_idx(g_lo, iglo)
!           if (nonlin .and. it == 1 .and. ik == 1) cycle
!           akperp4 = (aky(ik)**2 + akx(it)**2)**nexp
!            do isgn=1,2
!              do ig=-ntgrid, ntgrid-1
                
!                 havg = favg (g   (ig  ,isgn,iglo), &
!                              g   (ig+1,isgn,iglo), &
!                              gnew(ig  ,isgn,iglo), &
!                              gnew(ig+1,isgn,iglo)) 

!                 j0aparavg = favg (aj0(ig  ,iglo) * apar(ig  ,it,ik), &
!                                  aj0(ig+1,iglo)  * apar(ig+1,it,ik), &
!                                  aj0(ig  ,iglo)  * aparnew(ig  ,it,ik), &
!                                  aj0(ig+1,iglo)  * aparnew(ig+1,it,ik)) & 
!                                  * fapar * spec(is)%zstm * vpac(ig,isgn,iglo)

! !Set g0 for hyperresistive heating
!                 g0(ig,isgn,iglo) = spec(is)%dens*spec(is)%temp*D_eta*akperp4* &
!                      conjg(havg)*j0aparavg

!              end do
!           end do
!        end do

!        call integrate_moment (g0, tot)
!        if (proc0) then
!           do ik = 1, naky
!              fac2 = 0.5
!              if (aky(ik) < epsilon(0.0)) fac2 = 1.0
!              do it = 1, ntheta0
!                 if (nonlin .and. it == 1 .and. ik == 1) cycle
!                 do is = 1, nspec
!                    do ig = -ntgrid, ntgrid-1
!                       hk(it,ik) % hyperres(is) = hk(it,ik) % hyperres(is) &
!                            + real(tot(ig,it,ik,is))*wgt(ig)*fac2
!                    end do
!                    h % hyperres(is) = h % hyperres(is) + hk(it,ik) % hyperres(is)
!                 end do
!              end do
!           end do
!        end if

!     end if !End Hyperresistivity Heating Calculation

! !==========================================================================
! !Finish Energy-------------------------------------------------------------
! !==========================================================================

! !GGH Calculate hs2-------------------------------------------------------------
!     do iglo=g_lo%llim_proc, g_lo%ulim_proc
!        is = is_idx(g_lo, iglo)
!        it = it_idx(g_lo, iglo)
!        ik = ik_idx(g_lo, iglo)
!        if (nonlin .and. it == 1 .and. ik == 1) cycle
!        do isgn=1,2

!           do ig=-ntgrid, ntgrid-1
             
!              havg = favg (g   (ig  ,isgn,iglo), &
!                           g   (ig+1,isgn,iglo), &
!                           gnew(ig  ,isgn,iglo), &
!                           gnew(ig+1,isgn,iglo))

!              g0(ig,isgn,iglo) = 0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg)
!           end do
!        end do
!     end do

!     call integrate_moment (g0, tot)

!     if (proc0) then
!        do ik = 1, naky
!           fac2 = 0.5
!           if (aky(ik) < epsilon(0.0)) fac2 = 1.0
!           do it = 1, ntheta0
!              if (nonlin .and. it == 1 .and. ik == 1) cycle
!              do is = 1, nspec             
!                 do ig = -ntgrid, ntgrid-1

!                    !hs2 = int_r int_v T/F0 hs^2/2
!                    hk(it,ik) % hs2(is) = hk(it,ik) % hs2(is)  &
!                         + real(tot(ig,it,ik,is))*wgt(ig)*fac2

!                 end do
!              end do
!              h % hs2(:) = h % hs2(:) + hk(it,ik) % hs2(:)
!           end do
!        end do
!     end if

! !Calculate phis2-------------------------------------------------------------
!     if (proc0) then
!        do ik=1,naky
!           fac2 = 0.5
!           if (aky(ik) < epsilon(0.0)) fac2 = 1.0
!           do it = 1,ntheta0
!              if (nonlin .and. it == 1 .and. ik == 1) cycle
!              do ig=-ntgrid, ntgrid-1
!                 do is = 1, nspec
!                    phi_avg = favg (phi   (ig  ,it,ik), &
!                         phi   (ig+1,it,ik), &
!                         phinew(ig  ,it,ik), &
!                         phinew(ig+1,it,ik)) * fphi * spec(is)%zt

!                    !hs2 = int_r int_v T/F0 hs^2/2
!                    hk(it,ik) % phis2(is) = hk(it,ik) % phis2(is)  &
!                         +0.5*spec(is)%temp*spec(is)%dens*real(conjg(phi_avg)*phi_avg) &
!                         * wgt(ig) * fac2
!                 enddo
!              end do
!              h % phis2(:) = h % phis2(:) + hk(it,ik) % phis2(:)
!           end do
!        end do
!     endif

! ! Calculate delfs2 (rest of energy)-----------------------------------------------

! !GGH  Include response of Boltzmann species for single species runs
!     if (.not. has_electron_species(spec)) then
!        if (proc0) then
!           !NOTE: It is assumed here that n0i=n0e and zi=-ze
!           do ik=1,naky
!              fac2 = 0.5
!              if (aky(ik) < epsilon(0.0)) fac2 = 1.0
!              dtinv = 1./(code_dt*tunits(ik))
!              do it = 1,ntheta0

!                 if (nonlin .and. it == 1 .and. ik == 1) cycle

!                 do ig=-ntgrid, ntgrid-1

!                    phi_avg = favg (phi   (ig  ,it,ik), &
!                         phi   (ig+1,it,ik), &
!                         phinew(ig  ,it,ik), &
!                         phinew(ig+1,it,ik))

!                    !NOTE: Adiabatic (Boltzmann) species has temperature
!                    !       T = spec(1)%temp/tite
!                    hk(it,ik) % energy = hk(it,ik) % energy + &
!                         fphi * real(conjg(phi_avg)*phi_avg) &
!                         * 0.5 * spec(1)%dens * spec(1)%z * spec(1)%z * (tite/spec(1)%temp) &
!                         * wgt(ig)*fac2

!                 end do
!              end do
!           end do
!        endif
!     endif !END Correction to energy for single species runs---------------------

!     do iglo=g_lo%llim_proc, g_lo%ulim_proc
!        is = is_idx(g_lo, iglo)
!        it = it_idx(g_lo, iglo)
!        ik = ik_idx(g_lo, iglo)
!        if (nonlin .and. it == 1 .and. ik == 1) cycle
!        do isgn=1,2

!           do ig=-ntgrid, ntgrid-1
             
!              havg = favg (g   (ig  ,isgn,iglo), &
!                           g   (ig+1,isgn,iglo), &
!                           gnew(ig  ,isgn,iglo), &
!                           gnew(ig+1,isgn,iglo))

!              j0phiavg = favg (aj0(ig  ,iglo)*phi   (ig  ,it,ik), &
!                               aj0(ig+1,iglo)*phi   (ig+1,it,ik), &  
!                               aj0(ig  ,iglo)*phinew(ig  ,it,ik), &  
!                               aj0(ig+1,iglo)*phinew(ig+1,it,ik)) * fphi * spec(is)%zt

!              phi_avg = favg (phi   (ig  ,it,ik), &
!                              phi   (ig+1,it,ik), &
!                              phinew(ig  ,it,ik), &
!                              phinew(ig+1,it,ik)) * fphi * spec(is)%zt

!              g0(ig,isgn,iglo) = 0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg &
!                   + conjg(phi_avg)*phi_avg &
!                   - conjg(j0phiavg)*havg &
!                   - conjg(havg)*j0phiavg)
!           end do
!        end do
!     end do

!     call integrate_moment (g0, tot)

!     if (proc0) then
!        do ik = 1, naky
!           fac2 = 0.5
!           if (aky(ik) < epsilon(0.0)) fac2 = 1.0
!           do it = 1, ntheta0
!              if (nonlin .and. it == 1 .and. ik == 1) cycle
!              do is = 1, nspec             
!                 do ig = -ntgrid, ntgrid-1
!                    hk(it,ik) % energy = hk(it,ik) % energy  &
!                         + real(tot(ig,it,ik,is))*wgt(ig)*fac2

!                    !Delfs2 = int_r int_v T/F0 dfs^2/2
!                    hk(it,ik) % delfs2(is) = hk(it,ik) % delfs2(is)  &
!                         + real(tot(ig,it,ik,is))*wgt(ig)*fac2
!                 end do
!              end do
!              h % energy = h % energy + hk(it,ik) % energy
!              h % delfs2(:) = h % delfs2(:) + hk(it,ik) % delfs2(:)
!           end do
!        end do
!        deallocate (wgt)
!     end if

!     deallocate (tot)

! !!
! !! Put g, gnew back to their usual meanings
! !!
!     call g_adjust (g,    phi,    bpar,    -fphi, -fbpar)
!     call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar)

!   end subroutine get_heat

  subroutine reset_init

    use dist_fn_arrays, only: gnew, g, source
    initialized = .false.
    
    wdrift = 0.
    source = 0.
    gnew = 0.
    g0 = 0.
    g = 0.

  end subroutine reset_init

  subroutine reset_physics

    call init_wstar

  end subroutine reset_physics

  subroutine write_f (last)

    use mp, only: proc0, send, receive
    use file_utils, only: open_output_file, close_output_file, get_unused_unit
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx, imu_idx
    use gs2_layouts, only: idx_local, proc_id
    use vpamu_grids, only: nvgrid, nmu, vpa, vperp2, mu
    use gs2_time, only: user_time
    use dist_fn_arrays, only: g, gnew

    integer :: iglo, ik, it, is, iv, ig, imu
    integer, save :: unit
    complex, dimension (-nvgrid:nvgrid) :: gtmp
    logical, save :: first = .true.
    logical, intent(in)  :: last 

    if (first .and. proc0) then
       call open_output_file (unit, ".dist")
       write(unit, *) (2*nvgrid+1)*nmu
    endif
    first = .false.

    do iglo = g_lo%llim_world, g_lo%ulim_world
       ! writing out g(vpar,vperp) at ik=it=is=1, ig=0
       ik = ik_idx(g_lo, iglo) ; if (ik /= 1) cycle
       it = it_idx(g_lo, iglo) ; if (it /= 1) cycle
       is = is_idx(g_lo, iglo) ; if (is /= 1) cycle
       imu = imu_idx(g_lo, iglo) 
       ig = 0
       if (idx_local (g_lo, ik, it, imu, is)) then
          if (proc0) then 
             gtmp = gnew(ig,:,iglo)
          else
             call send (gnew(ig,:,iglo), 0)
          endif
       else if (proc0) then
          call receive (gtmp, proc_id(g_lo, iglo))
       endif
       if (proc0) then
          do iv = -nvgrid, nvgrid
             write (unit, "(8(1x,e12.6))") vpa(iv), sqrt(vperp2(ig,imu)), mu(imu), &
                  real(gtmp(1)), real(gtmp(2))
          end do
       end if
    end do
    if (proc0) write (unit, *)
    if (last .and. proc0) call close_output_file (unit)
    
  end subroutine write_f

!   subroutine write_fyx (phi,bpar,last)

!     use mp, only: proc0, send, receive, barrier
!     use file_utils, only: open_output_file, close_output_file, flush_output_file
!     use gs2_layouts, only: iv_idx, ig_idx, ik_idx, it_idx, is_idx, imu_idx
!     use gs2_layouts, only: idx_local, proc_id, yxf_lo, g_lo
!     use vpamu_grids, only: nvgrid, nmu
!     use kt_grids, only: naky, ntheta0, nx, ny
!     use theta_grid, only: bmag, ntgrid
!     use species, only: nspec
!     use dist_fn_arrays, only: gnew, g_adjust
!     use gs2_transforms, only: transform2, init_transforms
!     use run_parameters, only: fphi, fbpar

!     complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
!     logical, intent (in) :: last

!     real, dimension (:,:), allocatable :: grs, gzf
!     real :: gp0, gp0zf
!     integer :: ig, it, ik, is, iyxlo, iv, imu, iglo
!     integer, save :: unit
! !    logical :: first = .true.

!     allocate (grs(yxf_lo%ny, yxf_lo%llim_proc:yxf_lo%ulim_alloc))
!     allocate (gzf(yxf_lo%ny, yxf_lo%llim_proc:yxf_lo%ulim_alloc))
!     grs = 0.0; gzf = 0.0; gp0 = 0.0; gp0zf = 0.0

! !    if (first) then
!     if (.not. fyxinit) then
!        if (proc0) then
!           call open_output_file (unit,".yxdist")
!           write(unit,*) nvgrid*nmu, bmag(0)
!        end if

! !       first = .false.
!        fyxinit = .true.
!     end if

!     call g_adjust (gnew, phi, bpar, fphi, fbpar)

!     g0 = gnew

!     call transform2 (g0, grs)

!     g0 = 0.0
!     do iglo=g_lo%llim_proc, g_lo%ulim_proc
!        ik = ik_idx(g_lo,iglo)
!        if (ik == 1) g0(:,:,iglo) = gnew(:,:,iglo)
!     end do

!     call g_adjust (gnew, phi, bpar, -fphi, -fbpar)

!     call transform2 (g0, gzf)

!     do iyxlo=yxf_lo%llim_world, yxf_lo%ulim_world
          
!        ig = ig_idx(yxf_lo, iyxlo)
!        it = it_idx(yxf_lo, iyxlo)
!        if (ig == 0 .and. it == 1) then
!           ik = 1
!           iv = iv_idx(yxf_lo, iyxlo)
!           imu = imu_idx(yxf_lo, iyxlo)
!           is = is_idx(yxf_lo, iyxlo)
             
!           if (proc0) then
!              if (.not. idx_local(yxf_lo, ig, iv, it, imu, is)) then
!                    call receive (gp0, proc_id(yxf_lo, iyxlo))
!                    call receive (gp0zf, proc_id(yxf_lo, iyxlo))
!                 else
!                    gp0 = grs(ik, iyxlo)
!                    gp0zf = gzf(ik, iyxlo)
!                 end if
!              else if (idx_local(yxf_lo, ig, iv, it, imu, is)) then
!                 call send (grs(ik, iyxlo), 0)
!                 call send (gzf(ik, iyxlo), 0)                   
!              end if
             
!              if (proc0) then
!                 write (unit, "(4(1x,e12.6))") vpa(iv), mu(imu), &
!                      gp0, gp0zf
!              end if
!           end if
!        end if
!        call barrier
!     end do
!     deallocate (grs, gzf)

!     if (proc0) call flush_output_file (unit, ".yxdist")

!     if (proc0) then
!        write(unit,*)
!        if (last) call close_output_file (unit)
!     end if

!   end subroutine write_fyx

  subroutine boundary(linked)

    logical :: linked

    call init_dist_fn
    linked = boundary_option_switch == boundary_option_linked

  end subroutine boundary

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

  ! need to change this to have actual t_imp and thet_imp values

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

    use gs2_layouts, only: g_lo, imu_idx
    use species, only: nspec, spec
    use kt_grids, only: nakx => ntheta0, naky
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: aj0,aj1,kperp2
    use vpamu_grids, only: integrate_moment, vpa, vperp2, nvgrid

    integer :: i, imu
    integer :: iv, iglo
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
    allocate(gtmp(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
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
             imu = imu_idx(g_lo,iglo)
             do iv = -nvgrid, nvgrid
                if(i==1) wgt(:)=aj0(:,iglo)
                if(i==2) wgt(:)=aj0(:,iglo)*vpa(iv)**2
                if(i==3) wgt(:)=aj0(:,iglo)*vperp2(:,imu)
                if(i==4) wgt(:)=aj0(:,iglo)*vpa(iv)**4
                if(i==5) wgt(:)=aj0(:,iglo)*vpa(iv)**2*vperp2(:,imu)
                if(i==6) wgt(:)=vperp2(:,imu)*aj1(:,iglo)
                if(i==7) wgt(:)=vperp2(:,imu)*aj1(:,iglo)*vpa(iv)**2
                if(i==8) wgt(:)=vperp2(:,imu)*aj1(:,iglo)*vperp2(:,imu)
                gtmp(-ntgrid:ntgrid,iv,iglo) = wgt(-ntgrid:ntgrid)*cmplx(1.,0.)
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

  ! subroutine used for testing
  ! takes as input an array using g_lo and
  ! writes it to a .distmp output file
  subroutine write_mpdist (dist, extension, last)

    use mp, only: proc0, send, receive
    use file_utils, only: open_output_file, close_output_file
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use gs2_layouts, only: imu_idx, idx_local, proc_id
    use gs2_time, only: code_time
    use theta_grid, only: ntgrid, bmag, theta
    use vpamu_grids, only: vpa, nvgrid, mu
    use kt_grids, only: theta0

    implicit none
    
    complex, dimension (-ntgrid:,-nvgrid:,g_lo%llim_proc:), intent (in) :: dist
!    real, dimension (-ntgrid:,-nvgrid:,g_lo%llim_proc:), intent (in) :: dist
    character (*), intent (in) :: extension
    logical, intent (in), optional :: last
    
    integer :: iglo, ik, it, is, imu, ig, iv
    integer, save :: unit
    logical, save :: done = .false.
    complex :: gtmp
!    real, dimension (2) :: gtmp
    
    if (.not. done) then
       !        if (proc0) call open_output_file (unit, ".distmp")
       if (proc0) call open_output_file (unit, trim(extension))
       do iglo=g_lo%llim_world, g_lo%ulim_world
          ik = ik_idx(g_lo, iglo)
          it = it_idx(g_lo, iglo)
          is = is_idx(g_lo, iglo) ; if (is /= 1) cycle
          imu = imu_idx(g_lo, iglo)
          do iv = -nvgrid, nvgrid
             do ig = -ntgrid, ntgrid
                if (idx_local (g_lo, ik, it, imu, is)) then
                   if (proc0) then
                      gtmp = dist(ig,iv,iglo)
                   else
                      call send (dist(ig,iv,iglo), 0)
                   end if
                else if (proc0) then
                   call receive (gtmp, proc_id(g_lo, iglo))
                end if
                if (proc0) then
                   write (unit,'(a1,8e14.4)') "", code_time, theta(ig), vpa(iv), mu(imu), bmag(ig), &
                        real(gtmp), aimag(gtmp), theta(ig)-theta0(it,ik)
                end if
             end do
          end do
          if (proc0) then
             write (unit,*)
             write (unit,*)
          end if
       end do
       if (proc0) call close_output_file (unit)
       if (present(last)) done = .true.
    end if
    
  end subroutine write_mpdist
  
  subroutine finish_dist_fn

    use vpamu_grids, only: vperp2, energy, anon
    use dist_fn_arrays, only: aj0, aj1, kperp2
    use dist_fn_arrays, only: g, gnew, kx_shift, source

    implicit none

    no_comm = .false.
    readinit = .false. ; bessinit = .false. ; kp2init = .false. ; connectinit = .false.
    feqinit = .false. ; lpolinit = .false. ; fyxinit = .false. ; cerrinit = .false. ; mominit = .false.
    increase = .true. ; decrease = .false.

    call reset_init

    if (allocated(wdrift)) deallocate (wdrift)
    if (allocated(vperp2)) deallocate (vperp2, energy, anon)
    if (allocated(wstar)) deallocate (wstar)
    if (allocated(aj0)) deallocate (aj0, aj1)
    if (allocated(kperp2)) deallocate (kperp2)
    if (allocated(l_links)) deallocate (l_links, r_links, n_links)
    if (allocated(M_class)) deallocate (M_class, N_class)
    if (allocated(itleft)) deallocate (itleft, itright)
    if (allocated(connections)) deallocate (connections)
    if (allocated(g_adj)) deallocate (g_adj)
    if (allocated(g)) deallocate (g, gnew, g0, source)
    if (allocated(gexp_1)) deallocate (gexp_1, gexp_2, gexp_3)
    if (allocated(g_h)) deallocate (g_h, save_h)
    if (allocated(kx_shift)) deallocate (kx_shift)
    if (allocated(jump)) deallocate (jump)
    if (allocated(ikx_indexed)) deallocate (ikx_indexed)
    if (allocated(ufac)) deallocate (ufac)
    if (allocated(gamtot)) deallocate (gamtot, gamtot1, gamtot2)
    if (allocated(gamtot3)) deallocate (gamtot3)
    if (allocated(fl_avg)) deallocate (fl_avg)
    if (allocated(awgt)) deallocate (awgt)
    if (allocated(mom_coeff)) deallocate (mom_coeff, mom_coeff_npara, mom_coeff_nperp, &
         mom_coeff_tpara, mom_coeff_tperp, mom_shift_para, mom_shift_perp)

    if (allocated(source0)) deallocate (source0)
    if (allocated(gresponse)) deallocate (gresponse)
    if (allocated(m_mat)) deallocate (m_mat)

    if (allocated(pppfac)) then
       deallocate (pppfac)
       deallocate (ppmfac)
       deallocate (pmpfac)
       deallocate (pmmfac)
       deallocate (mppfac)
       deallocate (mpmfac)
       deallocate (mmpfac)
       deallocate (mmmfac)
    end if

    if (allocated(ig_low)) then
       deallocate (ig_low)
       deallocate (ig_mid)
       deallocate (ig_up)
    end if

  end subroutine finish_dist_fn

#ifdef LOWFLOW
  subroutine init_lowflow

    use constants, only: zi
    use dist_fn_arrays, only: vparterm, wdfac, vpac
    use dist_fn_arrays, only: wstarfac, hneoc, vpar
    use species, only: spec, nspec
    use geometry, only: rhoc
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
       allocate (wstar_neo(-ntgrid:ntgrid,ntheta0,naky))
    end if
    vparterm = 0. ; wdfac = 0. ; hneoc = 1. ; wstarfac = 0. ; wstar_neo = 0.

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
                if (.not. forbid(ig0,il)) then
                   write (neo_unit,'(10e12.4)') sign(sqrt(energy(ie)*(1.-al(il)*bmag(ig0))),1.5-real(isgn)), &
                        sqrt(energy(ie)*al(il)*bmag(ig0)), energy(ie), &
                        sign(sqrt(1.-al(il)*bmag(ig0)),1.5-real(isgn)), &
                        tmp1(ig0,il,ie,isgn,1), tmp2(ig0,il,ie,isgn,1), tmp3(ig0,il,ie,isgn,1), &
                        tmp4(ig0,il,ie,isgn,1), tmp5(ig0,il,ie,isgn,1), tmp6(ig0,il,ie,isgn,1)
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
       
       ! this is the contribution from v_E^par . grad F0 (Sec. 5.7 of gs2_notes.pdf)
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
!    wdfac = 0. ; cdfac = 0. ; hneoc = 1. ; wstarfac = 0. ; vparterm = 0.
 
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
       wstarfac(:,:,iglo) = wstar(:,:,iglo)*hneoc(:,:,iglo) - wstarfac(:,:,iglo)
    end do

    deallocate (tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9)
    deallocate (vpadhdec,dhdec,dhdxic,cdfac,hneovpth)

  end subroutine init_lowflow
#endif

end module dist_fn
