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
  public :: read_parameters
  public :: timeadv, g_exb
  public :: getan, getfieldeq, getmoms
!  public :: get_total_g
  public :: flux, eexchange, get_fieldcorrection
  public :: t0, omega0, gamma0
  public :: reset_init, write_f, reset_physics
  public :: M_class, N_class, i_class, par_spectrum
  public :: l_links, r_links, itright, itleft, boundary
  public :: init_kperp2
  public :: get_init_field, write_mpdist
  public :: flux_vs_theta_vs_vpa
!  public :: get_omega_prime, omprim_converge
  public :: profile_variation

  public :: gamtot,gamtot1,gamtot2
  public :: gamtotp, gamtota
  public :: omega_prime, omega_correction
  public :: boundary_option_switch, boundary_option_linked
  public :: boundary_option_self_periodic, boundary_option_zero

  private

  ! TMP FOR TESTING -- MAB
  real, parameter :: globalfac1 = 1.0
  real, parameter :: globalfac2 = 1.0

  real :: apfac, poisfac, driftknob
  real :: t0, omega0, gamma0
  real :: afilter, kfilter
  real :: wfb, g_exb, g_exbfac, omprimfac, btor_slab, mach
  logical :: dfexist, skexist, nonad_zero, lf_default, lf_decompose
  logical :: vpa_bc_zero, theta_bc_zero, profile_variation

  integer :: ntg_out

  integer :: adiabatic_option_switch
  integer, parameter :: adiabatic_option_default = 1, &
       adiabatic_option_zero = 2, &
       adiabatic_option_fieldlineavg = 3, &
       adiabatic_option_yavg = 4

  integer :: boundary_option_switch
  integer, parameter :: boundary_option_zero = 1, &
       boundary_option_self_periodic = 2, &
       boundary_option_linked = 3
  logical, public :: def_parity, even
  logical :: test
  logical :: increase = .true., decrease = .false.
!  logical :: omprim_converge = .false., omcorr_converge = .false.

  ! internal arrays

! #ifdef LOWFLOW
!   real, dimension (:,:), allocatable :: wcurv
!   ! (-ntgrid:ntgrid, -g-layout-)
!   real, dimension (:,:,:), allocatable :: wstar_neo
!   ! (-ntgrid:ntgrid,ntheta0,naky)
! #endif

  ! can probably do away with a full array for wdrift at grid points,
  ! as it is only needed in limited region of phase space
  real, dimension (:,:,:,:), allocatable :: wdrift, wdriftc, wdriftp, wdriftpc
  ! (-ntgrid:ntgrid, -nvgrid:nvgrid, ntheta0, -g-layout-)
  real, dimension (:,:,:), allocatable :: wdriftmod, wdriftmodc
  ! (-ntgrid:ntgrid, -nvgrid:nvgrid, nmu)

  ! can probably do away with a full array for wstar at grid points,
  ! as it is only needed in limited region of phase space
  real, dimension (:,:,:), allocatable :: wstar, wstarc, wstarp, wstarpc
  ! (-ntgrid:ntgrid, -nvgrid:nvgrid, -g-layout-)

  ! fieldeq
  real, dimension (:,:,:), allocatable :: gamtot, gamtot1, gamtot2, gamtot3
  real, dimension (:,:,:), allocatable :: gamtota
  real, dimension (:,:,:), allocatable :: gamtotp ! needed for radial profile variation
  ! (-ntgrid:ntgrid,ntheta0,naky) replicated

  !> following arrays needed for implicit solve

  ! coefficients multiplying g_{i+1,j+1}^{n+1}, g_{i+1,j+1}^{n}, g_{i+1,j}^{n+1}, etc.
  complex, dimension (:,:,:,:), allocatable :: pppfac, ppmfac, pmpfac, pmmfac
  complex, dimension (:,:,:,:), allocatable :: mppfac, mpmfac, mmpfac, mmmfac
  ! (-ntgrid:ntgird, -nvgrid:nvgrid, ntheta0, -g-layout-)

  ! response matrix in g
  complex, dimension (:,:,:,:), allocatable :: gresponse
  ! matrix needed for reponse matrix approach
  complex, dimension (:,:,:,:), allocatable :: m_mat
  ! ((ntheta/2)*nseg_max+1, (ntheta/2)*nseg_max+1, neigen_max, -g-layout-)

  ! special source term for points where dvpa/dt = 0
  complex, dimension (:,:,:), allocatable :: source0
  ! (nseg_max, ntheta0, -g-layout-)
  complex, dimension (:,:,:,:), allocatable :: mu0_source
  real :: decay_fac

  ! these arrays needed to keep track of connections between different
  ! 2pi segments
  integer :: nseg_max, neigen_max
  integer, dimension (:), allocatable :: neigen
  integer, dimension (:), allocatable :: ig_low, ig_mid, ig_up, it_shift_left
  integer, dimension (:,:), allocatable :: nsegments, ir_up, it_shift
  integer, dimension (:,:,:), allocatable :: itmod
  logical, dimension (:), allocatable :: periodic

  ! max number of elements in response matrix
  ! for periodic, this is ntheta/2 * nsegments + 1 + 2*nvgrid
!  ! for periodic, this is ntheta/2 * nsegments + 1 + 2*nvgrid+1
  ! (the number of grid points for vpa=0, theta below midplane,
  ! plus grid points at theta=-ntgrid, vpa>0 and theta=ntgrid, vpa<0).
!  ! plus grid points at theta=-ntgrid, vpa>=0 and theta=ntgrid, vpa<0).
  ! for non-periodic, nresponse = ntheta/2 * nsegments + 1.
  ! no need for grid points at vpa /= 0 since do not need to enfore periodicity
  integer :: nresponse

  !< end arrays for implicit solve

!  complex, dimension (:,:,:), allocatable :: g_h
  ! (-ntgrid:ntgrid, 2, -g-layout-)

  complex, dimension (:,:,:,:), allocatable :: g0
  ! (-ntgrid:ntgrid, -nvgrid:nvgrid, ntheta0, -g-layout-)

!  complex, dimension (:,:,:,:), allocatable :: g_adj
  ! (N(links), -nvgrid:nvgrid, ntheta0, -g-layout-)

  complex, dimension (:,:,:,:), allocatable, save :: gexp_1, gexp_2, gexp_3
  ! (-ntgrid:ntgrid, -nvgrid:nvgrid, ntheta0, -g-layout-)

  !> factors needed for inclusion of profile variation

  real, dimension (:,:,:,:), allocatable :: varfac, varfacc
  ! (-ntgrid:ntgrid,-nvgrid:nvgrid,nmu,nspec)

  real, dimension (:,:,:,:), allocatable :: mirror
  ! (-ntgrid:ntgrid,-nvgrid:nvgrid,nmu,nspec)

  real, dimension (:,:), allocatable :: streamfac
  ! (-ntgrid:ntgrid,2)

  complex, dimension (:,:,:), allocatable :: omega_prime, omega_correction
  ! (naky, ntheta0, navg)

  !< end profile variation factors

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
  type (redist_type), save :: pass_right
  type (redist_type), save :: pass_left

  integer, dimension (:,:), allocatable :: l_links, r_links
!  integer, dimension (:,:,:), allocatable :: n_links
!  logical, dimension (:,:), allocatable :: save_h
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

  interface write_mpdist
     module procedure write_mpdist_real, write_mpdist_complex
  end interface

contains

  subroutine init_dist_fn

    use mp, only: proc0, finish_mp, mp_abort
    use species, only: init_species, nspec
    use theta_grid, only: init_theta_grid, ntgrid, ntheta, nperiod
    use kt_grids, only: init_kt_grids, naky, ntheta0
    use vpamu_grids, only: init_vpamu_grids, nmu, nvgrid
    use run_parameters, only: init_run_parameters
    use gs2_layouts, only: init_dist_fn_layouts, init_gs2_layouts
    use nonlinear_terms, only: init_nonlinear_terms
    use hyper, only: init_hyper

    implicit none

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

    if (debug) write(6,*) "init_dist_fn: dist_fn_layouts"
    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nvgrid, nmu, nspec)

    if (debug) write(6,*) "init_dist_fn: nonlinear_terms"
    call init_nonlinear_terms 

    if (debug) write(6,*) "init_dist_fn: allocate_arrays"
    call allocate_arrays

    if (debug) write(6,*) "init_dist_fn: connections"
    call init_connections

    if (debug) write(6,*) "init_dist_fn: kperp2"
    call init_kperp2

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

    ntg_out = ntheta/2 + (nperiod-1)*ntheta

  end subroutine init_dist_fn

  subroutine read_parameters

    use file_utils, only: input_unit, error_unit, input_unit_exist
    use theta_grid, only: shat
    use text_options, only: text_option, get_option_value
    use species, only: nspec
    use mp, only: proc0, broadcast
    use theta_grid, only: itor_over_B

    implicit none

    type (text_option), dimension (7), parameter :: boundaryopts = &
         (/ text_option('default', boundary_option_zero), &
            text_option('zero', boundary_option_zero), &
            text_option('unconnected', boundary_option_zero), &
            text_option('self-periodic', boundary_option_self_periodic), &
            text_option('periodic', boundary_option_self_periodic), &
            text_option('kperiod=1', boundary_option_self_periodic), &
            text_option('linked', boundary_option_linked) /)
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
         driftknob, poisfac, adiabatic_option, vpa_bc_zero, &
         kfilter, afilter, test, def_parity, even, wfb, theta_bc_zero, &
         g_exb, g_exbfac, omprimfac, btor_slab, mach, lf_default, lf_decompose, &
         profile_variation
    
    namelist /source_knobs/ t0, omega0, gamma0
    integer :: ierr, in_file

    if (readinit) return
    readinit = .true.

    if (proc0) then
       boundary_option = 'default'
       nonad_zero = .true.
       vpa_bc_zero = .false.
       theta_bc_zero = .true.
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
       profile_variation = .false.
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
    call broadcast (vpa_bc_zero)
    call broadcast (theta_bc_zero)
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
    call broadcast (profile_variation)
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
    use theta_grid, only: ntgrid, thet_imp
    use kt_grids, only: ntheta0
    use vpamu_grids, only: nvgrid, vpa_imp, nmu
    use gs2_layouts, only: g_lo, ik_idx, imu_idx, is_idx

    implicit none

    integer :: iglo, ig, ik, it, iv, imu, is
    logical :: debug = .false.

    if (.not. allocated(wdrift)) then
       ! allocate wdrift with sign(vpa) dependence because will contain 
       ! Coriolis as well as magnetic drifts
       allocate (wdrift(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (wdriftc(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (wdriftp(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (wdriftpc(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (wdriftmod(-ntgrid:ntgrid,-nvgrid:nvgrid,nmu))
       allocate (wdriftmodc(-ntgrid:ntgrid,-nvgrid:nvgrid,nmu))
    end if
    wdrift = 0. ; wdriftc = 0. ; wdriftp = 0. ; wdriftpc = 0. ; wdriftmod = 0. ; wdriftmodc = 0.
! #ifdef LOWFLOW
!     if (.not. allocated(wcurv)) allocate (wcurv(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
!     wcurv = 0.
! #endif

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do it = 1, ntheta0
          do iv = -nvgrid, nvgrid
             do ig = -ntgrid, ntgrid
                imu=imu_idx(g_lo,iglo)
                ik=ik_idx(g_lo,iglo)
                is=is_idx(g_lo,iglo)
                ! get grad-B and curvature drifts
                wdrift(ig,iv,it,iglo) = wdrift_func(ig,iv,imu,it,ik)*driftknob
                wdriftp(ig,iv,it,iglo) = wdriftp_func(ig,iv,imu,it,ik)*driftknob
! #ifdef LOWFLOW
!              ! get curvature drift without vpa dependence
!              wcurv(ig,iglo) = wcurv_func(ig, it, ik)*driftknob
! #endif
             ! add Coriolis drift to magnetic drifts
                wdrift(ig,iv,it,iglo) = wdrift(ig,iv,it,iglo) + wcoriolis_func(ig,iv,it,ik,is)
             end do
          end do
       end do
    end do

!    call write_mpdist (wdrift,'.wdrift',last=.true.)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do it = 1, ntheta0
          ! get wdrift at center of theta-vpa cells
          call get_cell_value (thet_imp, vpa_imp, &
               wdrift(:,:,it,iglo), wdriftc(:,:,it,iglo), -ntgrid, -nvgrid)
          call get_cell_value (thet_imp, vpa_imp, &
               wdriftp(:,:,it,iglo), wdriftpc(:,:,it,iglo), -ntgrid, -nvgrid)
! #ifdef LOWFLOW
!        ! get wcurv at center of theta cell
!        call get_cell_value (thet_imp, wcurv(:,iglo), wcurv(:,iglo), -ntgrid)
! #endif
       end do
    end do

    do imu = 1, nmu
       do iv = -nvgrid, nvgrid
          do ig = -ntgrid, ntgrid
             wdriftmod(ig,iv,imu) = wdriftmod_func(ig,iv,imu)*driftknob
          end do
       end do
    end do
    do imu = 1, nmu
       call get_cell_value (thet_imp, vpa_imp, &
            wdriftmod(:,:,imu), wdriftmodc(:,:,imu), -ntgrid, -nvgrid)
    end do

    if (debug) write(6,*) 'init_wdrift: driftknob',driftknob

  end subroutine init_wdrift

  function wdrift_func (ig, iv, imu, it, ik)

    use theta_grid, only: gbdrift, gbdrift0, cvdrift, cvdrift0
    use theta_grid, only: shat
    use kt_grids, only: aky, theta0, akx
    use vpamu_grids, only: vpa, vperp2
    use run_parameters, only: wunits
    use gs2_time, only: code_dt

    implicit none

    real :: wdrift_func
    integer, intent (in) :: ig, ik, it, iv, imu

    ! note that wunits=aky/2 (for wstar_units=F)
    if (aky(ik) < epsilon(0.0)) then
       wdrift_func = akx(it)/shat &
            * (cvdrift0(ig)*vpa(iv)**2 + gbdrift0(ig)*0.5*vperp2(ig,imu)) &
            * code_dt/2.0
    else
       wdrift_func = ((cvdrift(ig) + theta0(it,ik)*cvdrift0(ig)) * vpa(iv)**2 &
            + (gbdrift(ig) + theta0(it,ik)*gbdrift0(ig)) * 0.5*vperp2(ig,imu)) &
            * code_dt*wunits(ik)
    end if
  end function wdrift_func

  function wdriftp_func (ig, iv, imu, it, ik)

    use theta_grid, only: shat, dcvdrift0drho, dgbdrift0drho
    use theta_grid, only: dcvdriftdrho, dgbdriftdrho
    use kt_grids, only: aky, theta0, akx
    use vpamu_grids, only: vpa, mu, nmu
    use gs2_time, only: code_dt
    use run_parameters, only: wunits

    implicit none

    real :: wdriftp_func
    integer, intent (in) :: ig, ik, it, iv, imu

    ! note that wunits=aky/2 (for wstar_units=F)
    if (aky(ik) < epsilon(0.0)) then
       wdriftp_func = akx(it)/shat &
            * (dcvdrift0drho(ig)*vpa(iv)**2 + dgbdrift0drho(ig)*mu(imu)) &
            * code_dt*0.5
    else
       wdriftp_func = ((dcvdriftdrho(ig) + theta0(it,ik)*dcvdrift0drho(ig)) * vpa(iv)**2 &
            + (dgbdriftdrho(ig) + theta0(it,ik)*dgbdrift0drho(ig)) * mu(imu)) &
            * code_dt*wunits(ik)
    end if

  end function wdriftp_func

  function wdriftmod_func (ig, iv, imu)

    use geometry, only: rhoc
    use theta_grid, only: gbdrift0, cvdrift0
    use theta_grid, only: shat, drhodpsi, qval
    use vpamu_grids, only: vpa, vperp2
    use gs2_time, only: code_dt
    use run_parameters, only: rhostar

    implicit none

    real :: wdriftmod_func
    integer, intent (in) :: ig, iv, imu

    ! this is vM . grad psi without the k_psi that usually accompanies it
    wdriftmod_func = rhostar*drhodpsi*(rhoc/qval/shat) &
         * (cvdrift0(ig)*vpa(iv)**2 + gbdrift0(ig)*0.5*vperp2(ig,imu)) &
         * code_dt*0.5

  end function wdriftmod_func

! #ifdef LOWFLOW
!   function wcurv_func (ig, it, ik)

!     use theta_grid, only: cvdrift, cvdrift0, shat
!     use kt_grids, only: aky, theta0, akx
!     use run_parameters, only: wunits
!     use gs2_time, only: code_dt

!     implicit none

!     real :: wcurv_func
!     integer, intent (in) :: ig, ik, it

!     if (aky(ik) == 0.0) then
!        wcurv_func = akx(it)/shat &
!             * cvdrift0(ig) * code_dt/2.0
!     else
!        wcurv_func = (cvdrift(ig) + theta0(it,ik)*cvdrift0(ig)) &
!             *code_dt*wunits(ik)
!     end if
!   end function wcurv_func
! #endif

  function wcoriolis_func (ig, iv, it, ik, is)

    use theta_grid, only: cdrift, cdrift0, shat
    use kt_grids, only: aky, theta0, akx
    use vpamu_grids, only: vpa
    use run_parameters, only: wunits
    use gs2_time, only: code_dt
    use species, only: spec

    implicit none

    real :: wcoriolis_func
    integer, intent (in) :: ig, ik, it, iv, is

    if (aky(ik) < epsilon(0.0)) then
       wcoriolis_func = mach * vpa(iv) &
            * cdrift0(ig) * code_dt * akx(it)/(2.*shat*spec(is)%stm)
    else
       wcoriolis_func = mach * vpa(iv) &
            * (cdrift(ig) + theta0(it,ik)*cdrift0(ig))*code_dt*wunits(ik)/spec(is)%stm
    end if

  end function wcoriolis_func

  subroutine init_vperp2

    use centering, only: get_cell_value
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: vpar, vparp
    use species, only: spec, nspec
    use theta_grid, only: ntgrid, bmag, thet_imp, dbdthetc, gradparc, delthet
    use theta_grid, only: dgradpardrhoc, dgradparbdrhoc
    use vpamu_grids, only: vperp2, nmu, mu, energy, anon, anonc
    use vpamu_grids, only: vpa, vpac, vpa_imp, nvgrid, dvpa

    implicit none

    integer :: iv, imu, is
    
    if (.not.allocated(vperp2)) then
       allocate (vperp2(-ntgrid:ntgrid,nmu)) ; vperp2 = 0.
       allocate (energy(-ntgrid:ntgrid,-nvgrid:nvgrid,nmu)) ; energy = 0.
       allocate (anon(-ntgrid:ntgrid,-nvgrid:nvgrid,nmu)) ; anon = 0.
       allocate (anonc(-ntgrid:ntgrid,-nvgrid:nvgrid,nmu)) ; anonc = 0.
       allocate (vpar(-ntgrid:ntgrid,-nvgrid:nvgrid,nspec)) ; vpar = 0.
       allocate (vparp(-ntgrid:ntgrid,-nvgrid:nvgrid,nspec)) ; vparp = 0.
       allocate (mirror(-ntgrid:ntgrid,-nvgrid:nvgrid,nmu,nspec)) ; mirror = 0.
       allocate (streamfac(-ntgrid:ntgrid,2)) ; streamfac = 0.
    endif

    vperp2 = 2.0*spread(mu,1,2*ntgrid+1)*spread(bmag,2,nmu)

    do iv = -nvgrid, nvgrid
       energy(:,iv,:) = vpa(iv)**2 + vperp2
       anon(:,iv,:) = exp(-energy(:,iv,:))
    end do

    do imu = 1, nmu
       call get_cell_value (thet_imp, vpa_imp, anon(:,:,imu), anonc(:,:,imu), &
            -ntgrid, -nvgrid)
    end do
    
    ! get the proper upwinded cell value for the parallel velocity
    do is = 1, nspec
       do iv = -nvgrid, -1
          where (dbdthetc(:ntgrid-1,1) < 0.0)
             vpar(:ntgrid-1,iv,is) = vpac(iv,1)*code_dt*spec(is)%zstm*gradparc(:ntgrid-1,1)/delthet(:ntgrid-1)
             vparp(:ntgrid-1,iv,is) = vpac(iv,1)*code_dt*spec(is)%stm*dgradpardrhoc(:ntgrid-1,1)/delthet(:ntgrid-1)
          elsewhere
             vpar(:ntgrid-1,iv,is) = vpac(iv,2)*code_dt*spec(is)%zstm*gradparc(:ntgrid-1,1)/delthet(:ntgrid-1)
             vparp(:ntgrid-1,iv,is) = vpac(iv,2)*code_dt*spec(is)%stm*dgradpardrhoc(:ntgrid-1,1)/delthet(:ntgrid-1)
          end where
       end do
       do iv = 0, nvgrid-1
          where (dbdthetc(:ntgrid-1,2) < 0.0)
             vpar(:ntgrid-1,iv,is) = vpac(iv,1)*code_dt*spec(is)%zstm*gradparc(:ntgrid-1,2)/delthet(:ntgrid-1)
             vparp(:ntgrid-1,iv,is) = vpac(iv,1)*code_dt*spec(is)%stm*dgradpardrhoc(:ntgrid-1,2)/delthet(:ntgrid-1)
          elsewhere
             vpar(:ntgrid-1,iv,is) = vpac(iv,2)*code_dt*spec(is)%zstm*gradparc(:ntgrid-1,2)/delthet(:ntgrid-1)
             vparp(:ntgrid-1,iv,is) = vpac(iv,2)*code_dt*spec(is)%stm*dgradpardrhoc(:ntgrid-1,2)/delthet(:ntgrid-1)
          end where
       end do

       do imu = 1, nmu
          ! this factor will be needed for the mirror term appearing as a source in the gprim equation
          mirror(:ntgrid-1,-nvgrid:-1,imu,is) = -code_dt*spec(is)%stm*mu(imu) &
               / spread(dvpa(-nvgrid:-1),1,2*ntgrid)*spread(dgradparbdrhoc(:ntgrid-1,1),2,nvgrid)
          mirror(:ntgrid-1,0:nvgrid-1,imu,is) = -code_dt*spec(is)%stm*mu(imu) &
               / spread(dvpa(0:nvgrid-1),1,2*ntgrid)*spread(dgradparbdrhoc(:ntgrid-1,2),2,nvgrid)
       end do
    end do

    ! this factor will be needed for the streaming term appearing as a source in the gprim equation
    streamfac(:ntgrid-1,:) = dgradpardrhoc(:ntgrid-1,:)/gradparc(:ntgrid-1,:)

  end subroutine init_vperp2

  subroutine init_wstar

    use centering, only: get_cell_value
    use geometry, only: d2psidr2
    use species, only: spec, nspec
    use theta_grid, only: ntgrid, thet_imp, itor_over_b, dBdrho, drhodpsi
    use vpamu_grids, only: nvgrid, vpa, vpa_imp, nmu, energy, mu
    use run_parameters, only: wunits
    use gs2_time, only: code_dt
    use gs2_layouts, only: is_idx, ik_idx, imu_idx, g_lo

    implicit none

    integer :: ik, is, imu, ig, iv, iglo

    if (.not. allocated(wstar)) then
       allocate (wstar(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; wstar = 0.0
       allocate (wstarc(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; wstarc = 0.0
       allocate (wstarp(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; wstarp = 0.0
       allocate (wstarpc(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; wstarpc = 0.0
       allocate (varfac(-ntgrid:ntgrid,-nvgrid:nvgrid,nmu,nspec)) ; varfac = 0.0
       allocate (varfacc(-ntgrid:ntgrid,-nvgrid:nvgrid,nmu,nspec)) ; varfacc = 0.0
    end if

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       imu = imu_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       do iv = -nvgrid, nvgrid
          do ig = -ntgrid, ntgrid
             wstar(ig,iv,iglo) = code_dt*wunits(ik) &
                  ! contribution from background rotation gradient source term
                  *( -2.0*vpa(iv)*g_exb*itor_over_b(ig)/spec(is)%stm &
                  ! contribution from background density and temperature gradients
                  + spec(is)%fprim+spec(is)%tprim*(energy(ig,iv,imu)-1.5))
             ! fdbprim = (1/n)*d2n/dr2, tdbprim = (1/T)*d2T/dr2
             wstarp(ig,iv,iglo) = code_dt*wunits(ik) &
                  * ( -2.*spec(is)%fprim*spec(is)%tprim*(energy(ig,iv,imu)-1.5) &
                  - 2.*mu(imu)*dBdrho(ig)*(spec(is)%fprim+spec(is)%tprim*(energy(ig,iv,imu)-2.5)) &
                  - spec(is)%fdbprim - spec(is)%tdbprim*(energy(ig,iv,imu)-1.5) &
                  - spec(is)%tprim**2*((energy(ig,iv,imu)-1.5)*(energy(ig,iv,imu)-2.5)-energy(ig,iv,imu)) &
                  - (spec(is)%fprim + spec(is)%tprim*(energy(ig,iv,imu)-1.5))*drhodpsi*d2psidr2 )
          end do
       end do

       ! get (theta-vpa) cell values for wstar from grid values
       call get_cell_value (thet_imp, vpa_imp, &
            wstar(:,:,iglo), wstarc(:,:,iglo), -ntgrid, -nvgrid)
       call get_cell_value (thet_imp, vpa_imp, &
            wstarp(:,:,iglo), wstarpc(:,:,iglo), -ntgrid, -nvgrid)
    end do

    ! open (unit=1007,file='wstarp.out',status='unknown')
    ! do ig = -ntgrid, ntgrid
    !    is = 1 ; iv = 48 ; imu = 6 ; ik = 1
    !    write (1007,'(3e12.4)'), theta(ig), wunits(ik)*(spec(is)%fprim+spec(is)%tprim*(energy(ig,iv,imu)-1.5)), &
    !         wunits(ik) * ( -2.*spec(is)%fprim*spec(is)%tprim*(energy(ig,iv,imu)-1.5) &
    !         - 2.*mu(imu)*dBdrho(ig)*(spec(is)%fprim+spec(is)%tprim*(energy(ig,iv,imu)-2.5)) &
    !         - spec(is)%fdbprim - spec(is)%tdbprim*(energy(ig,iv,imu)-1.5) &
    !         - spec(is)%tprim**2*((energy(ig,iv,imu)-1.5)*(energy(ig,iv,imu)-2.5)-energy(ig,iv,imu)) &
    !         - (spec(is)%fprim + spec(is)%tprim*(energy(ig,iv,imu)-1.5))*drhodpsi*d2psidr2 ) &
    !         + wunits(ik)*(spec(is)%fprim+spec(is)%tprim*(energy(ig,iv,imu)-1.5))*(spec(is)%fprim &
    !         + spec(is)%tprim*(energy(ig,iv,imu)-1.5)+2.*mu(imu)*dBdrho(ig))
    ! end do
    ! close (1007)

    ! also calculate factor needed when including profile variation
    ! this is (T/F_M) * d/drho (F_M/T)
    do is = 1, nspec
       do imu = 1, nmu
          do iv = -nvgrid, nvgrid
             varfac(:,iv,imu,is) = -( spec(is)%fprim + spec(is)%tprim*(energy(:,iv,imu)-2.5) &
                  + 2.*mu(imu)*dBdrho )
          end do
          call get_cell_value (thet_imp, vpa_imp, &
               varfac(:,:,imu,is), varfacc(:,:,imu,is), -ntgrid, -nvgrid)
       end do
    end do

  end subroutine init_wstar

  subroutine init_bessel

    use dist_fn_arrays, only: aj0, aj1, kperp2, aj0p, dkperp2dr
    use species, only: spec
    use theta_grid, only: ntgrid, bmag, dBdrho
    use kt_grids, only: ntheta0
    use vpamu_grids, only: vperp2, mu
    use gs2_layouts, only: g_lo, ik_idx, imu_idx, is_idx
    use spfunc, only: j0, j1

    implicit none

    integer :: ig, ik, it, imu, is
    integer :: iglo
    real :: arg

    if (bessinit) return
    bessinit = .true.

    call init_kperp2

    allocate (aj0(-ntgrid:ntgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (aj1(-ntgrid:ntgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (aj0p(-ntgrid:ntgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
    aj0 = 0. ; aj1 = 0. ; aj0p = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       imu = imu_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid
             arg = spec(is)%smz*sqrt(vperp2(ig,imu)*kperp2(ig,it,ik))/bmag(ig)
             aj0(ig,it,iglo) = j0(arg)
             ! note that j1 returns and aj1 stores J_1(x)/x (NOT J_1(x)), 
             aj1(ig,it,iglo) = j1(arg)
             ! aj0p = dJ0/dr = -(J1/arg)*arg*darg/dr = -(1/2)*(J1/arg)*d(arg**2)/dr
             aj0p(ig,it,iglo) = -aj1(ig,it,iglo)*mu(imu)/(bmag(ig)*spec(is)%zstm**2) &
                  * (dkperp2dr(ig,it,ik) - kperp2(ig,it,ik)*dBdrho(ig)/bmag(ig))
          end do
       end do
    end do

!     if (proc0) then
!        open (unit=1008,file='aj0p.out',status='unknown')
!        is = 1 ; imu = 6 ; ik = 1 ; it = 1
!        do ig = -ntgrid, ntgrid
!           arg = spec(is)%smz*sqrt(vperp2(ig,imu)*kperp2(ig,it,ik))/bmag(ig)          
!           write (1008,'(8e12.4)') theta(ig), j0(arg), -j1(arg)*mu(imu)/(bmag(ig)*spec(is)%zstm**2) &
!                * (dkperp2dr(ig,it,ik) - kperp2(ig,it,ik)*dBdrho(ig)/bmag(ig)), &
! !               + j1(arg)*mu(imu)/(bmag(ig)*spec(is)%zstm**2)*kperp2(ig,it,ik)*spec(is)%tprim, &
!                kperp2(ig,it,ik), dkperp2dr(ig,it,ik), bmag(ig), dBdrho(ig), mu(imu)*spec(is)%temp
!        end do
!        close (1008)
!     end if

  end subroutine init_bessel

  subroutine init_kperp2

    use dist_fn_arrays, only: kperp2, dkperp2dr
    use theta_grid, only: ntgrid, gds2, gds21, gds22, shat
    use theta_grid, only: dgds2dr, dgds21dr, dgds22dr
    use kt_grids, only: naky, ntheta0, aky, theta0, akx

    implicit none

    integer :: ik, it

    if (kp2init) return
    kp2init = .true.

    allocate (kperp2(-ntgrid:ntgrid,ntheta0,naky))
    allocate (dkperp2dr(-ntgrid:ntgrid,ntheta0,naky))
    do ik = 1, naky
       if (aky(ik) < epsilon(0.0)) then
         do it = 1, ntheta0
             kperp2(:,it,ik) = akx(it)*akx(it)*gds22/(shat*shat)
             dkperp2dr(:,it,ik) = akx(it)*akx(it)*dgds22dr/shat**2
          end do
       else
          do it = 1, ntheta0
             kperp2(:,it,ik) = aky(ik)*aky(ik) &
                  *(gds2 + 2.0*theta0(it,ik)*gds21 &
                  + theta0(it,ik)*theta0(it,ik)*gds22)
             ! this is dkperp^2/drho
             dkperp2dr(:,it,ik) = aky(ik)**2 &
                  *(dgds2dr + 2.0*theta0(it,ik)*dgds21dr &
                  + theta0(it,ik)**2*dgds22dr)
          end do
       end if
    end do

  end subroutine init_kperp2

  subroutine init_par_filter
    use theta_grid, only: ntgrid
    use gs2_transforms, only: init_zf
    use kt_grids, only: naky, ntheta0

    if ( naky*ntheta0 .eq. 0 ) then
       print *,"WARNING: kt_grids used in init_par_filter before initialised?"
    endif

# if FFT == _FFTW_
    call init_zf (ntgrid)
# elif FFT == _FFTW3_
    call init_zf (ntgrid, ntheta0*naky)
# endif
    
  end subroutine init_par_filter

  subroutine par_spectrum(an, an2)

    use gs2_transforms, only: kz_spectrum
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0

    complex, dimension(:,:,:) :: an, an2    
    real :: scale

    call kz_spectrum (an, an2, ntheta0, naky)
    scale = 1./real(4*ntgrid**2)
    an2 = an2*scale

  end subroutine par_spectrum

  subroutine init_implicit_solve

    use constants, only: zi
    use gs2_layouts, only: g_lo, imu_idx, is_idx, ik_idx
    use dist_fn_arrays, only: vpar
    use gs2_time, only: code_dt
    use run_parameters, only: t_imp, rhostar
    use species, only: spec, nspec
    use theta_grid, only: ntgrid, dbdthetc, gradparc, thet_imp, ntheta
    use vpamu_grids, only: nvgrid, vpa_imp, dvpa, vpa, mu
    use kt_grids, only: naky, ntheta0
    use centering, only: get_cell_value

    implicit none

    integer :: iglo, imu, ig, iv, ntg, is, it, ik, igl, igm, igu, iseg, ie
    real :: thm_fac, vpm_fac, vpp_fac, thp_fac, stm

    real, dimension (:,:), allocatable :: dum1, dum2
    complex, dimension (:,:), allocatable :: wd1, wd2

#ifdef LOWFLOW
    real, dimension (:,:), allocatable :: vp
    real, dimension (:,:), allocatable :: cvdrift_thc, gbdrift_thc
#endif

!    call init_connections

    ntg = ntheta/2

    if (.not. allocated(gresponse)) then
       allocate (gresponse(nresponse,nresponse,neigen_max,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (m_mat(nresponse,nresponse,neigen_max,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (source0(nseg_max,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (mu0_source(-ntgrid:ntgrid,ntheta0,naky,nspec)) ; mu0_source = 0.0
       allocate (pppfac(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc)) ; pppfac = 0.0
       allocate (ppmfac(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc)) ; ppmfac = 0.0
       allocate (pmpfac(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc)) ; pmpfac = 0.0
       allocate (pmmfac(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc)) ; pmmfac = 0.0
       allocate (mppfac(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc)) ; mppfac = 0.0
       allocate (mpmfac(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc)) ; mpmfac = 0.0
       allocate (mmpfac(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc)) ; mmpfac = 0.0
       allocate (mmmfac(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc)) ; mmmfac = 0.0
    end if

    if (.not. allocated(dum1)) then
       allocate (dum1(-ntgrid:ntgrid-1,-nvgrid:nvgrid-1)) ; dum1 = 0.0
       allocate (dum2(-ntgrid:ntgrid-1,-nvgrid:nvgrid-1)) ; dum2 = 0.0
    end if
    if (.not. allocated(wd1)) then
       allocate (wd1(-ntgrid:ntgrid-1,-nvgrid:nvgrid)) ; wd1 = 0.0
       allocate (wd2(-ntgrid:ntgrid-1,-nvgrid:nvgrid)) ; wd2 = 0.0
    end if
#ifdef LOWFLOW
    if (.not. allocated(vp)) then
       allocate (vp(-ntgrid:ntgrid-1,-nvgrid:nvgrid-1)) ; vp = 0.0
    end if
    if (.not. allocated(cvdrift_thc)) then
       allocate (cvdrift_thc(-ntgrid:ntgrid-1,2)) ; cvdrift_thc = 0.0
       allocate (gbdrift_thc(-ntgrid:ntgrid-1,2)) ; gbdrift_thc = 0.0
    end if
#endif

    thm_fac = 0. ; vpm_fac = 0.

    if (vpa_bc_zero) then
       decay_fac = 0.0
    else
       decay_fac = exp(dvpa(-nvgrid)*(dvpa(-nvgrid)-2.*abs(vpa(-nvgrid))))
    end if

    do iglo = g_lo%llim_proc, g_lo%ulim_proc

       ! get im corresponding to iglo
       imu = imu_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)

       dum1 = vpar(:ntgrid-1,:nvgrid-1,is)/spec(is)%zstm

       ! this is sqrt( vt_s / vt_ref )
       stm = spec(is)%stm
       ! TMP FOR TESTING -- MAB
!       stm = spec(2)%stm
         
       if (imu == 1) then
          where (dum1 > 0.0)
             dum2 = thet_imp
          elsewhere
             dum2 = 1.0-thet_imp
          end where
          do ie = 1, neigen(ik)
             do iseg = 1, nsegments(ie,ik)
                
                it = itmod(iseg,ie,ik)
                igl = ig_low(iseg) ; igm = ig_mid(iseg) ; igu = ig_up(iseg)

                ! wd1 and wd2 account for contributions from wdrift
                wd1 = 1.0 + t_imp*zi*wdriftc(:ntgrid-1,:,it,iglo)*spec(is)%tz
                wd2 = 1.0 + (t_imp-1.0)*zi*wdriftc(:ntgrid-1,:,it,iglo)*spec(is)%tz

                ! this is for quadrants 1 and 4
                pppfac(igl:igm-1,-nvgrid:nvgrid-1,it,iglo) = &
                     dum2(igl:igm-1,-nvgrid:nvgrid-1)*wd1(igl:igm-1,-nvgrid:nvgrid-1) &
                     + stm*(t_imp*dum1(igl:igm-1,-nvgrid:nvgrid-1))
                pmpfac(igl:igm-1,-nvgrid:nvgrid-1,it,iglo) = &
                     (1.0-dum2(igl:igm-1,-nvgrid:nvgrid-1))*wd1(igl:igm-1,-nvgrid:nvgrid-1) &
                     - stm*t_imp*dum1(igl:igm-1,-nvgrid:nvgrid-1)
                mppfac(igl:igm-1,-nvgrid:nvgrid-1,it,iglo) = &
                     -dum2(igl:igm-1,-nvgrid:nvgrid-1)*wd2(igl:igm-1,-nvgrid:nvgrid-1) &
                     + stm*((1.0-t_imp)*dum1(igl:igm-1,-nvgrid:nvgrid-1))
                mmpfac(igl:igm-1,-nvgrid:nvgrid-1,it,iglo) = &
                     -(1.0-dum2(igl:igm-1,-nvgrid:nvgrid-1))*wd2(igl:igm-1,-nvgrid:nvgrid-1) &
                     + stm*(-(1.0-t_imp)*dum1(igl:igm-1,-nvgrid:nvgrid-1))
                ppmfac(igl:igm-1,-nvgrid:nvgrid-1,it,iglo) = 0.0
                pmmfac(igl:igm-1,-nvgrid:nvgrid-1,it,iglo) = 0.0
                mpmfac(igl:igm-1,-nvgrid:nvgrid-1,it,iglo) = 0.0
                mmmfac(igl:igm-1,-nvgrid:nvgrid-1,it,iglo) = 0.0
                
                ! this is for quadrants 2 and 3
                pppfac(igm:igu-1,-nvgrid:nvgrid-1,it,iglo) = 0.0
                pmpfac(igm:igu-1,-nvgrid:nvgrid-1,it,iglo) = 0.0
                mppfac(igm:igu-1,-nvgrid:nvgrid-1,it,iglo) = 0.0
                mmpfac(igm:igu-1,-nvgrid:nvgrid-1,it,iglo) = 0.0
                ppmfac(igm:igu-1,-nvgrid:nvgrid-1,it,iglo) = &
                     dum2(igm:igu-1,-nvgrid:nvgrid-1)*wd1(igm:igu-1,-nvgrid:nvgrid-1) &
                     + stm*t_imp*dum1(igm:igu-1,-nvgrid:nvgrid-1)
                pmmfac(igm:igu-1,-nvgrid:nvgrid-1,it,iglo) = &
                     (1.0-dum2(igm:igu-1,-nvgrid:nvgrid-1))*wd1(igm:igu-1,-nvgrid:nvgrid-1) &
                     - stm*t_imp*dum1(igm:igu-1,-nvgrid:nvgrid-1)
                mpmfac(igm:igu-1,-nvgrid:nvgrid-1,it,iglo) = &
                     -dum2(igm:igu-1,-nvgrid:nvgrid-1)*wd2(igm:igu-1,-nvgrid:nvgrid-1) &
                     + stm*(1.0-t_imp)*dum1(igm:igu-1,-nvgrid:nvgrid-1)
                mmmfac(igm:igu-1,-nvgrid:nvgrid-1,it,iglo) = &
                     -(1.0-dum2(igm:igu-1,-nvgrid:nvgrid-1))*wd2(igm:igu-1,-nvgrid:nvgrid-1) &
                     - stm*(1.0-t_imp)*dum1(igm:igu-1,-nvgrid:nvgrid-1)
                
             end do
          end do
       else

          ! get mirror force term multiplying dg/dvpa
          ! will multiply by sqrt(vt_s / vt_r) later
          ! sign here is consistent with appearance on LHS of GK equation
#ifdef LOWFLOW
          ! if next order correction to GK equation is included then the mirror term
          ! that involves dg/dvpa is modified to include a curvature drift contribution 
          ! (see Sec. 6.3 of gs3_notes.pdf)

          ! get the proper upwinded cell value for the parallel velocity
          do iv = -nvgrid, -1
             where (dbdthetc(:ntgrid-1,1) < 0.0)
                vp(:,iv) = vpac(iv,1)
             elsewhere
                vp(:,iv) = vpac(iv,2)
             end where
          end do
          do iv = 0, nvgrid-1
             where (dbdthetc(:ntgrid-1,2) < 0.0)
                vp(:,iv) = vpac(iv,1)
             elsewhere
                vp(:,iv) = vpac(iv,2)
             end where
          end do

          ! get cell-centered values for some geometric coefficients
          call get_cell_value (1.0-thet_imp, cvdrift_th, cvdrift_thc(:,1), -ntgrid)
          call get_cell_value (thet_imp, cvdrift_th, cvdrift_thc(:,2), -ntgrid)
          call get_cell_value (1.0-thet_imp, gbdrift_th, gbdrift_thc(:,1), -ntgrid)
          call get_cell_value (thet_imp, gbdrift_th, gbdrift_thc(:,2), -ntgrid)

          dum2(:,-nvgrid:-1) = -code_dt*mu(imu)*spread(dbdthetc(:ntgrid-1,1),2,nvgrid) &
               / spread(dvpa(-nvgrid:-1),1,2*ntgrid) &
               * (stm*spread(gradparc(:ntgrid-1,1),2,nvgrid) &
               + rhostar*vp(:,-nvgrid:-1)*spread(cvdrift_thc(:,1)-gbdrift_thc(:,1),2,nvgrid)*spec(is)%tz)
          dum2(:,0:nvgrid-1) = -code_dt*mu(imu)*spread(dbdthetc(:ntgrid-1,2),2,nvgrid) &
               / spread(dvpa(0:nvgrid-1),1,2*ntgrid) &
               * (stm*spread(gradparc(:ntgrid-1,2),2,nvgrid) &
               + rhostar*vp(:,0:nvgrid-1)*spread(cvdrift_thc(:,2)-gbdrift_thc(:,2),2,nvgrid)*spec(is)%tz)
#else
          dum2(:,-nvgrid:-1) = -code_dt*mu(imu)*stm &
               / spread(dvpa(-nvgrid:-1),1,2*ntgrid)*spread(dbdthetc(:ntgrid-1,1)*gradparc(:ntgrid-1,1),2,nvgrid)
          dum2(:,0:nvgrid-1) = -code_dt*mu(imu)*stm &
               / spread(dvpa(0:nvgrid-1),1,2*ntgrid)*spread(dbdthetc(:ntgrid-1,2)*gradparc(:ntgrid-1,2),2,nvgrid)
#endif     
                  
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

                do it = 1, ntheta0

                   ! wd1 and wd2 account for contributions from wdrift
                   wd1(ig,iv) = 1.0 + t_imp*zi*wdriftc(ig,iv,it,iglo)*spec(is)%tz
                   wd2(ig,iv) = 1.0 + (t_imp-1.0)*zi*wdriftc(ig,iv,it,iglo)*spec(is)%tz
                   
                   ! note that wd1, dum1, and dum2 have code_dt hidden inside
                   
                   pppfac(ig,iv,it,iglo) = thp_fac*vpp_fac*wd1(ig,iv) &
                        + (stm*t_imp*vpp_fac*dum1(ig,iv) &
                        + t_imp*thp_fac*dum2(ig,iv))
                   ppmfac(ig,iv,it,iglo) = thp_fac*vpm_fac*wd1(ig,iv) &
                        + (stm*t_imp*vpm_fac*dum1(ig,iv) &
                        - t_imp*thp_fac*dum2(ig,iv))
                   pmpfac(ig,iv,it,iglo) = thm_fac*vpp_fac*wd1(ig,iv) &
                        + (-stm*t_imp*vpp_fac*dum1(ig,iv) &
                        + t_imp*thm_fac*dum2(ig,iv))
                   pmmfac(ig,iv,it,iglo) = thm_fac*vpm_fac*wd1(ig,iv) &
                        - (stm*t_imp*vpm_fac*dum1(ig,iv) &
                        + t_imp*thm_fac*dum2(ig,iv))
                   mppfac(ig,iv,it,iglo) = -thp_fac*vpp_fac*wd2(ig,iv) &
                        + (stm*(1.0-t_imp)*vpp_fac*dum1(ig,iv) &
                        + (1.0-t_imp)*thp_fac*dum2(ig,iv))
                   mpmfac(ig,iv,it,iglo) = -thp_fac*vpm_fac*wd2(ig,iv) &
                        + (stm*(1.0-t_imp)*vpm_fac*dum1(ig,iv) &
                        - (1.0-t_imp)*thp_fac*dum2(ig,iv))
                   mmpfac(ig,iv,it,iglo) = -thm_fac*vpp_fac*wd2(ig,iv) &
                        + (-stm*(1.0-t_imp)*vpp_fac*dum1(ig,iv) &
                        + (1.0-t_imp)*thm_fac*dum2(ig,iv))
                   mmmfac(ig,iv,it,iglo) = -thm_fac*vpm_fac*wd2(ig,iv) &
                        - (stm*(1.0-t_imp)*vpm_fac*dum1(ig,iv) &
                        + (1.0-t_imp)*thm_fac*dum2(ig,iv))
                end do
             end do
          end do
       end if

       ! MAB - TMP FOR TESTING
       do it = 1, ntheta0
          ! treat theta<0, vpa=-vpa_max and theta>0, vpa=vpa_max specially
          ! in particular, assume g decays like a Maxwellian in vpa at this point
          where (dbdthetc(:ntgrid-1,1) < 0.0)
             pmpfac(:ntgrid-1,-nvgrid,it,iglo) = pmpfac(:ntgrid-1,-nvgrid,it,iglo) &
                  + pmmfac(:ntgrid-1,-nvgrid,it,iglo)*decay_fac
             pmmfac(:ntgrid-1,-nvgrid,it,iglo) = 0.
             ! should not need to modify anything other than pmpfac and pmmfac for these points
             pppfac(:ntgrid-1,-nvgrid,it,iglo) = pppfac(:ntgrid-1,-nvgrid,it,iglo) &
                  + ppmfac(:ntgrid-1,-nvgrid,it,iglo)*decay_fac
             mppfac(:ntgrid-1,-nvgrid,it,iglo) = mppfac(:ntgrid-1,-nvgrid,it,iglo) &
                  + mpmfac(:ntgrid-1,-nvgrid,it,iglo)*decay_fac
             mmpfac(:ntgrid-1,-nvgrid,it,iglo) = mmpfac(:ntgrid-1,-nvgrid,it,iglo) &
                  + mmmfac(:ntgrid-1,-nvgrid,it,iglo)*decay_fac
             ppmfac(:ntgrid-1,-nvgrid,it,iglo) = 0.
             mpmfac(:ntgrid-1,-nvgrid,it,iglo) = 0.
             mmmfac(:ntgrid-1,-nvgrid,it,iglo) = 0.
          elsewhere
             ppmfac(:ntgrid-1,nvgrid-1,it,iglo) = ppmfac(:ntgrid-1,nvgrid-1,it,iglo) &
                  + pppfac(:ntgrid-1,nvgrid-1,it,iglo)*decay_fac
             pppfac(:ntgrid-1,nvgrid-1,it,iglo) = 0.0
             ! should not need to modify anything other than ppmfac and pppfac for these points
             pmmfac(:ntgrid-1,nvgrid-1,it,iglo) = pmmfac(:ntgrid-1,nvgrid-1,it,iglo) &
                  + pmpfac(:ntgrid-1,nvgrid-1,it,iglo)*decay_fac
             mpmfac(:ntgrid-1,nvgrid-1,it,iglo) = mpmfac(:ntgrid-1,nvgrid-1,it,iglo) &
                  + mppfac(:ntgrid-1,nvgrid-1,it,iglo)*decay_fac
             mmmfac(:ntgrid-1,nvgrid-1,it,iglo) = mmmfac(:ntgrid-1,nvgrid-1,it,iglo) &
                  + mmpfac(:ntgrid-1,nvgrid-1,it,iglo)*decay_fac
             pmpfac(:ntgrid-1,nvgrid-1,it,iglo) = 0.
             mppfac(:ntgrid-1,nvgrid-1,it,iglo) = 0.
             mmpfac(:ntgrid-1,nvgrid-1,it,iglo) = 0.
          end where
       end do

    end do

    ! if non-periodic BC along theta, then
    ! get response of g at theta<=0, vpa=-dvpa
    ! to unit impulses in g at theta<0, vpa=0
    ! note that this is not needed for mu=0 where different vpa points are not connected
    ! response matrix gets more complicated if enforcing periodic BCs in theta,
    ! in which case need to additionally get response of g at theta=pi-dtheta, vpa>0
    ! to unit impulses in g at theta=-pi, vpa>0, etc.
    call get_gresponse_matrix

!    call write_response (gresponse, '.gresponse')
!    call write_response (m_mat, '.mmat')

    if (allocated(dum1)) deallocate (dum1, dum2)
    if (allocated(wd1)) deallocate (wd1, wd2)
#ifdef LOWFLOW
    if (allocated(vp)) deallocate (vp)
    if (allocated(cvdrift_thc)) deallocate (cvdrift_thc, gbdrift_thc)
#endif

  end subroutine init_implicit_solve

  subroutine init_connections

    use mp, only: nproc, mp_abort
    use centering, only: init_centering
    use theta_grid, only: nperiod, ntgrid, ntheta
    use kt_grids, only: ntheta0, jtwist_out, naky, aky
    use species, only: nspec
    use vpamu_grids, only: nmu, nvgrid

    implicit none

    integer :: iseg, ik, ie, ntg, it

    ntg = ntheta/2

    if (.not. allocated(neigen)) allocate (neigen(naky))
    if (.not. allocated(periodic)) allocate (periodic(naky)) ; periodic = .false.

    if (boundary_option_switch==boundary_option_self_periodic) then
       periodic = .true.
    else
       do ik = 1, naky
          if (aky(ik) < epsilon(0.0)) periodic(ik) = .true.
       end do
    end if

    select case (boundary_option_switch)
    case (boundary_option_linked)

       if (nproc > naky*nmu*nspec) then
          write (*,*) 'Parallelization over kx not currently supported'
          write (*,*) 'with twist and shift boundary condition.  Aborting.'
          call mp_abort ('Parallelization over kx not currently supported with twist-and-shift BC.')
       end if

       ik = 1
       neigen(ik) = ntheta0
       if (naky > 1) then
          do ik = 2, naky
             ! must link different kx values at theta = +/- pi
             ! neigen is the number of independent eigenfunctions along the field line
             neigen(ik) = min((ik-1)*jtwist_out,ntheta0)
          end do
       end if

       neigen_max = maxval(neigen)

       if (.not. allocated(it_shift_left)) then
          allocate (it_shift_left(neigen_max))
          allocate (it_shift(ntheta0,naky)) ; it_shift = 0
       end if

       ! figure out how much to shift it by to get to
       ! the left-most (theta-theta0) in each set of connected 2pi segments
       ! note that theta0 goes from 0 to theta0_max and then from theta0_min back
       ! to -dtheta0
       do it = 1, neigen_max
          ! first ntheta0/2+1 theta0s are 0 and all positive theta0 values
          ! remainder are negative theta0s
          ! note that ntheta0 is always positive for box
          if (it <= ntheta0/2+1) then
             it_shift_left(it) = ntheta0/2-2*it+2
          else
!             it_shift_left(it) = 0
             it_shift_left(it) = 3*(ntheta0/2)-2*it+3
          end if
       end do

       do ik = 1, naky
          ! it_shift is how much to shift each it by to connect
          ! to the next theta0 (from most positive to most negative)
          do it = 1, ntheta0
             ! if theta0 is negative, then shifting to more negative
             ! theta0 corresponds to decreasing it
             if (it > ntheta0/2+1) then
                ! if theta0 is sufficiently negative, it has no
                ! more negative theta0 with which it can connect
                if (it-neigen(ik) >= ntheta0/2+1) then
                   it_shift(it,ik) = -neigen(ik)
                end if
             ! theta0 is positive
             else
                ! if theta0 is sufficiently positive, shifting to more
                ! negative theta0 corresponds to decreasing it
                if (it-neigen(ik) > 0) then
                   it_shift(it,ik) = -neigen(ik)
                ! if a positive theta0 connects to a negative theta0
                ! must do more complicated mapping of it
                else if (it-neigen(ik)+ntheta0 >= 2+ntheta0/2) then
                   it_shift(it,ik) = ntheta0 - neigen(ik)
                end if
             end if
          end do
       end do

       if (.not. allocated(nsegments)) then
          allocate (nsegments(neigen_max,naky))
!          allocate (ir_up(neigen_max,naky))
       end if

       do ik = 1, naky
          if (neigen(ik) == 0) then
             nsegments(:,ik) = 1
          else
             nsegments(:,ik) = (ntheta0-1)/neigen(ik)

             do ie = 1, mod(ntheta0-1,neigen(ik))+1
                nsegments(ie,ik) = nsegments(ie,ik) + 1
             end do
          end if
       end do

!       ir_up = ntg*nsegments+1

       nseg_max = maxval(nsegments)

       if (.not. allocated(ig_low)) then
          allocate (ig_low(nseg_max)) ; ig_low = -ntgrid
          allocate (ig_mid(nseg_max)) ; ig_mid = 0
          allocate (ig_up(nseg_max)) ; ig_up = ntgrid
       end if

       call init_connected_bc

    case default
       
!       neigen = 1 ; neigen_max = 1
       neigen = ntheta0 ; neigen_max = ntheta0
       
       if (.not. allocated(it_shift_left)) then
          allocate (it_shift_left(neigen_max))
          allocate (it_shift(ntheta0,naky))
       end if
       it_shift = 0 ; it_shift_left = 0
       
       if (.not. allocated(nsegments)) then
          allocate (nsegments(neigen_max,naky))
!          allocate (ir_up(neigen_max,naky))
       end if
       
       ! this is the number of 2pi poloidal segments in the extended theta domain,
       ! which is needed in initializing the reponse matrix and doing the implicit sweep
       nsegments = 2*(nperiod-1) + 1
       
       nseg_max = maxval(nsegments)
       
!       ir_up = ntg*nsegments+1
       
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

       i_class = 1
       if (.not. allocated(M_class)) then
          allocate (M_class(i_class))
          allocate (N_class(i_class))
       end if
       M_class = naky*ntheta0 ; N_class = 1

    end select

    if (.not. allocated(ir_up)) allocate (ir_up(neigen_max,naky))
    ir_up = ntg*nsegments+1

    if (.not. allocated(itmod)) allocate (itmod(nseg_max,neigen_max,naky))
    do ik = 1, naky
       if (periodic(ik)) ir_up(:,ik) = ir_up(:,ik) + 2*nvgrid
!       if (periodic(ik)) ir_up(:,ik) = ir_up(:,ik) + 2*nvgrid+1

       ! only do the following once for each independent set of theta0s
       ! the assumption here is that all kx are on processor and sequential
       do it = 1, neigen(ik)
          
          ! remap to start at theta0 = -theta0_max for this set of connected theta0s
          iseg = 1
          itmod(iseg,it,ik) = it + it_shift_left(it)
          if (nsegments(it,ik) > 1) then
             do iseg = 2, nsegments(it,ik)
                itmod(iseg,it,ik) = itmod(iseg-1,it,ik) + it_shift(itmod(iseg-1,it,ik),ik)
             end do
          end if
       end do
    end do

    call init_centering (nperiod, ig_low, ig_mid, ig_up)
!    nresponse = nseg_max*ntg + 1

    nresponse = nseg_max*ntg + 2*nvgrid + 1
!    nresponse = nseg_max*ntg+1 + 2*(nvgrid+1)
    
!    call write_connectinfo

  end subroutine init_connections

  subroutine init_connected_bc

    use theta_grid, only: nperiod, ntheta, theta
    use kt_grids, only: naky, ntheta0, aky, theta0
    use gs2_layouts, only: ik_idx, imu_idx
    use gs2_layouts, only: idx, proc_id
    use mp, only: iproc, max_allreduce
    use constants
    use redistribute, only: index_list_type, delete_list
    use vpamu_grids, only: nvgrid

    implicit none

    integer :: ik, it, it0, itl, itr, jshift0
    integer :: i, j, k
    integer :: it_star
    integer :: ng
    integer, dimension(naky*ntheta0) :: n_k

    if (connectinit) return
    connectinit = .true.

    ng = ntheta/2 + (nperiod-1)*ntheta

    ! jshift0 corresponds to J (not delta j) from Beer thesis
    ! delta j is the number of akxs (i.e. theta0s) 'skipped' when connecting one 
    ! parallel segment to the next to satisfy the twist and shift
    ! boundary conditions. delta j = (ik - 1) * jshift0

    if (naky > 1 .and. ntheta0 > 1) then
       jshift0 = int((theta(ng)-theta(-ng))/(theta0(2,2)-theta0(1,2)) + 0.01)
    else if (naky == 1 .and. ntheta0 > 1 .and. aky(1) > epsilon(0.0)) then
       jshift0 = int((theta(ng)-theta(-ng))/(theta0(2,1)-theta0(1,1)) + 0.01)
    else
       jshift0 = 1
    end if
    
    allocate (itleft(naky,ntheta0), itright(naky,ntheta0))
    itleft(1,:) = -1 ! No connected segments when ky = 0.
    itright(1,:) = -1
    do ik = 1, naky
       do it = 1, ntheta0
          if (it > (ntheta0+1)/2) then
             ! akx is negative (-1 shift because akx(1) = 0)
             it0 = it - ntheta0 - 1
          else
             ! akx is positive (-1 shift because akx(1) = 0)
             it0 = it - 1
          end if

          if (ik == 1) then
             if (aky(ik) > epsilon(0.0) .and. naky == 1) then
                ! for this case, jshift0 is delta j from Beer thesis
                itl = it0 + jshift0
                itr = it0 - jshift0
             else
                itl = it0
                itr = it0
             end if
          else
             ! ik = 1 corresponds to aky=0, so shift index by 1
             ! (ik-1)*jshift0 is delta j from Beer thesis
             itl = it0 + (ik-1)*jshift0
             itr = it0 - (ik-1)*jshift0
          end if

          ! remap to usual GS2 indices
          if (itl >= 0 .and. itl < (ntheta0+1)/2) then
             itleft(ik,it) = itl + 1
          else if (itl + ntheta0 + 1 > (ntheta0+1)/2 &
             .and. itl + ntheta0 + 1 <= ntheta0) then
             itleft(ik,it) = itl + ntheta0 + 1
          else
             ! for periodic, need to change below -- MAB/BD
             ! jhat = j + delta j not included in simulation, so cannot connect
             itleft(ik,it) = -1
          end if

          ! same stuff for jhat = j - delta j
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

    ! allocate (connections(g_lo%llim_proc:g_lo%ulim_alloc))

    ! do iglo = g_lo%llim_proc, g_lo%ulim_proc
    !    ik = ik_idx(g_lo,iglo)
    !    is = is_idx(g_lo,iglo)
    !    imu = imu_idx(g_lo,iglo)

    !    do it = 1, ntheta0

    !       ! get processors and indices for jhat (kxhat) modes connecting
    !       ! to mode j (kx) so we can set up communication
          
    !       if (itleft(ik,it) < 0) then
    !          connections(iglo)%iproc_left = -1
    !          connections(iglo)%iglo_left = -1
    !       else
    !          connections(iglo)%iproc_left &
    !               = proc_id(g_lo,idx(g_lo,ik,itleft(ik,it),imu,is))
    !          connections(iglo)%iglo_left &
    !               = idx(g_lo,ik,itleft(ik,it),imu,is)
    !       end if
    !       if (itright(ik,it) < 0) then
    !          connections(iglo)%iproc_right = -1
    !          connections(iglo)%iglo_right = -1
    !       else
    !          connections(iglo)%iproc_right &
    !               = proc_id(g_lo,idx(g_lo,ik,itright(ik,it),imu,is))
    !          connections(iglo)%iglo_right &
    !               = idx(g_lo,ik,itright(ik,it),imu,is)
    !       end if
    !       if (connections(iglo)%iproc_left >= 0 .or. &
    !            connections(iglo)%iproc_right >= 0) then
    !          connections(iglo)%neighbor = .true.
    !       else
    !          connections(iglo)%neighbor = .false.
    !       end if
    !    end do
    ! end do

    if (boundary_option_switch == boundary_option_linked) then
       ! do iglo = g_lo%llim_proc, g_lo%ulim_proc          
       !    ik = ik_idx(g_lo,iglo)
       !    if (connections(iglo)%iglo_left >= 0 .and. aky(ik) /= 0.0) &
       !         save_h (1,iglo) = .true.
       !    if (connections(iglo)%iglo_right >= 0 .and. aky(ik) /= 0.0) &
       !         save_h (2,iglo) = .true.
       ! end do
          
       allocate (l_links(naky, ntheta0))
       allocate (r_links(naky, ntheta0))
!       allocate (n_links(2, naky, ntheta0))

!       n_links_max = 0
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
             ! ! 'n_links' complex numbers are needed to specify bc for (ik, it) region
             ! ! ignoring wfb
             ! ! n_links(1,:,:) is for v_par > 0, etc.
             ! if (l_links(ik, it) == 0) then
             !    n_links(1, ik, it) = 0
             ! else 
             !    n_links(1, ik, it) = 2*l_links(ik, it) - 1
             ! end if
             
             ! if (r_links(ik, it) == 0) then
             !    n_links(2, ik, it) = 0
             ! else 
             !    n_links(2, ik, it) = 2*r_links(ik, it) - 1
             ! end if
             ! n_links_max = max(n_links_max, n_links(1,ik,it), n_links(2,ik,it))
          end do
       end do
! wfb
!       if (n_links_max > 0) n_links_max = n_links_max + 3
       
       ! now set up communication pattern:
       ! excluding wfb
       
!        nn_to = 0
!        nn_from = 0 
       
!        do iglo = g_lo%llim_world, g_lo%ulim_world
          
!           ! Exclude trapped particles (inc wfb)
! !          if (nlambda > ng2 .and. il >= ng2+1) cycle

!           ip = proc_id(g_lo,iglo)
!           ik = ik_idx(g_lo,iglo)

!           do it = 1, ntheta0
!              ncell = r_links(ik, it) + l_links(ik, it) + 1
!              if (ncell == 1) cycle
             
!              iglo_star = iglo
!              do j = 1, r_links(ik, it)
!                 call get_right_connection (iglo_star, iglo_right, ipright)
                
!                 if (ip /= ipright) then
!                    ! sender
!                    if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
!                    ! receiver
!                    if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1
!                 end if
                
!                 iglo_star = iglo_right
!              end do
             
!              iglo_star = iglo
!              do j = 1, l_links(ik, it)
!                 call get_left_connection (iglo_star, iglo_left, ipleft)
!                 if (ip /= ipleft) then
!                    ! sender
!                    if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
!                    ! receiver
!                    if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
!                 end if
!                 iglo_star = iglo_left
!              end do

!           end do
!        end do

!        nn_max = maxval(nn_to)
!        call max_allreduce (nn_max)
!        if (nn_max == 0) then
!           no_comm = .true.
! !          write (*,*) 'no communication in init_connected_bc'
!           goto 200
!        end if
! !       write (*,*) 'communication in init_connected_bc'

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

! !          if (nlambda > ng2 .and. il >= ng2+1) cycle

!           ip = proc_id(g_lo,iglo)
!           ik = ik_idx(g_lo,iglo)

!           do it = 1, ntheta0
          
!              ncell = r_links(ik, it) + l_links(ik, it) + 1
!              if (ncell == 1) cycle
             
!              iglo_star = iglo
!              do j = 1, r_links(ik, it)
!                 call get_right_connection (iglo_star, iglo_right, ipright)
!                 ! sender
!                 if (ip == iproc) then
!                    n = nn_from(ipright) + 1
!                    nn_from(ipright) = n
!                    from(ipright)%first(n) = ntgrid
!                    from(ipright)%second(n) = 1
!                    from(ipright)%third(n) = iglo
!                 end if
!                 ! receiver
!                 if (ipright == iproc) then
!                    n = nn_to(ip) + 1
!                    nn_to(ip) = n
!                    to(ip)%first(n) = j 
!                    to(ip)%second(n) = 1
!                    to(ip)%third(n) = iglo_right
!                 end if
!                 iglo_star = iglo_right
!              end do
             
!              iglo_star = iglo
!              do j = 1, l_links(ik, it)
!                 call get_left_connection (iglo_star, iglo_left, ipleft)
!                 ! sender
!                 if (ip == iproc) then
!                    n = nn_from(ipleft) + 1
!                    nn_from(ipleft) = n
!                    from(ipleft)%first(n) = -ntgrid
!                    from(ipleft)%second(n) = 2
!                    from(ipleft)%third(n) = iglo
!                 end if
!                 ! receiver
!                 if (ipleft == iproc) then
!                    n = nn_to(ip) + 1
!                    nn_to(ip) = n
!                    to(ip)%first(n) = j 
!                    to(ip)%second(n) = 2
!                    to(ip)%third(n) = iglo_left
!                 end if
!                 iglo_star = iglo_left
!              end do
!           end do
!        end do

!        from_low (1) = -ntgrid
!        ! perhaps from_low(2), from_high(2), to_low(2), to_high(2)
!        ! below should be 1 and 2 instead of -nvgrid and nvgrid?
!        from_low (2) = -nvgrid
!        from_low (3) = g_lo%llim_proc
       
!        from_high(1) = ntgrid
!        from_high(2) = nvgrid
!        from_high(3) = g_lo%ulim_alloc

!        to_low (1) = 1
!        to_low (2) = -nvgrid
!        to_low (3) = g_lo%llim_proc
       
!        to_high(1) = n_links_max
!        to_high(2) = nvgrid
!        to_high(3) = g_lo%ulim_alloc

!        call init_fill (links_p, 'c', to_low, to_high, to, &
!             from_low, from_high, from)
       
!        call delete_list (from)
!        call delete_list (to)
       
! ! n_links_max is typically 2 * number of cells in largest supercell
!        allocate (g_adj (n_links_max, 2*nvgrid+1, g_lo%llim_proc:g_lo%ulim_alloc))

! ! now set up links_h:
! ! excluding wfb

!        nn_to = 0
!        nn_from = 0 

!        do iglo = g_lo%llim_world, g_lo%ulim_world

! !          if (nlambda > ng2 .and. il >= ng2+1) cycle

!           ip = proc_id(g_lo,iglo)
!           ik = ik_idx(g_lo,iglo)

!           do it = 1, ntheta0
          
!              ncell = r_links(ik, it) + l_links(ik, it) + 1
!              if (ncell == 1) cycle
             
!              ! If this is not the first link in the chain, continue
!              if (l_links(ik, it) > 0) then
!                 ! For each link to the right, do:
!                 iglo_star = iglo
!                 do j = 1, r_links(ik, it)
!                    call get_right_connection (iglo_star, iglo_right, ipright)
!                    ! sender
!                    if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
!                    ! receiver
!                    if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1
!                    iglo_star = iglo_right
!                 end do
!              end if
             
!              if (r_links(ik, it) > 0) then
!                 ! For each link to the left, do:
!                 iglo_star = iglo
!                 do j = 1, l_links(ik, it)
!                    call get_left_connection (iglo_star, iglo_left, ipleft)
!                    ! sender
!                    if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
!                    ! receiver
!                    if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
!                    iglo_star = iglo_left
!                 end do
!              end if
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

! !          if (nlambda > ng2 .and. il >= ng2+1) cycle

!           ip = proc_id(g_lo,iglo)
!           ik = ik_idx(g_lo,iglo)

!           do it = 1, ntheta0
          
!              ncell = r_links(ik, it) + l_links(ik, it) + 1
!              if (ncell == 1) cycle
             
!              if (l_links(ik, it) > 0) then
!                 ! For each link to the right, do:
!                 iglo_star = iglo
!                 do j = 1, r_links(ik, it)
!                    ! get address of link
!                    call get_right_connection (iglo_star, iglo_right, ipright)
!                    ! sender
!                    if (ip == iproc) then
!                       n = nn_from(ipright) + 1
!                       nn_from(ipright) = n
!                       from(ipright)%first(n) = ntgrid
!                       from(ipright)%second(n) = 1
!                       from(ipright)%third(n) = iglo
!                    end if
!                    ! receiver
!                    if (ipright == iproc) then
!                       n = nn_to(ip) + 1
!                       nn_to(ip) = n
!                       to(ip)%first(n) = 2*l_links(ik, it) + j
!                       to(ip)%second(n) = 1
!                       to(ip)%third(n) = iglo_right
!                    end if
!                    iglo_star = iglo_right
!                 end do
!              end if
             
!              if (r_links(ik, it) > 0) then
!                 ! For each link to the left, do:
!                 iglo_star = iglo
!                 do j = 1, l_links(ik, it)
!                    ! get address of link
!                    call get_left_connection (iglo_star, iglo_left, ipleft)
!                    ! sender
!                    if (ip == iproc) then
!                       n = nn_from(ipleft) + 1
!                       nn_from(ipleft) = n
!                       from(ipleft)%first(n) = -ntgrid
!                       from(ipleft)%second(n) = 2   
!                       from(ipleft)%third(n) = iglo
!                    end if
!                    ! receiver
!                    if (ipleft == iproc) then
!                       n = nn_to(ip) + 1
!                       nn_to(ip) = n
!                       to(ip)%first(n) = 2*r_links(ik, it) + j
!                       to(ip)%second(n) = 2
!                       to(ip)%third(n) = iglo_left
!                    end if
!                    iglo_star = iglo_left
!                 end do
!              end if
!           end do
!        end do

!        ! second index of arrays below perhaps should be -nvgrid, nvgrid?
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

! 200    continue

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

  ! subroutine get_left_connection (iglo, iglo_left, iproc_left)
  !   use gs2_layouts, only: g_lo, proc_id, idx
  !   use gs2_layouts, only: ik_idx, imu_idx, is_idx
  !   use kt_grids, only: ntheta0
  !   implicit none
  !   integer, intent (in) :: iglo
  !   integer, intent (out) :: iglo_left, iproc_left
  !   integer :: ik, is, imu

  !   ik = ik_idx(g_lo,iglo)
  !   imu = imu_idx(g_lo,iglo)
  !   is = is_idx(g_lo,iglo)
       
  !   iglo_left = idx(g_lo,ik,imu,is)
  !   iproc_left = proc_id(g_lo,iglo_left)
    
  ! end subroutine get_left_connection

  ! subroutine get_right_connection (iglo, iglo_right, iproc_right)
  !   use gs2_layouts, only: g_lo, proc_id, idx
  !   use gs2_layouts, only: ik_idx, imu_idx, is_idx
  !   implicit none
  !   integer, intent (in) :: iglo
  !   integer, intent (out) :: iglo_right, iproc_right
  !   integer :: ik, imu, is

  !   ik = ik_idx(g_lo,iglo)
  !   imu = imu_idx(g_lo,iglo)
  !   is = is_idx(g_lo,iglo)

  !   iglo_right = idx(g_lo,ik,imu,is)
  !   iproc_right = proc_id(g_lo,iglo_right)
  ! end subroutine get_right_connection

  subroutine get_gresponse_matrix

    use theta_grid, only: ntgrid, ntheta
    use vpamu_grids, only: nvgrid, energy
    use gs2_layouts, only: g_lo, ik_idx, imu_idx
    use matrix_inversion, only: invert_matrix
    use dist_fn_arrays, only: gnew, source

    implicit none

    integer :: ig, iglo, idx, ntg, iseg, it, ik, iup, imu, iv
    complex, dimension (:,:), allocatable :: p_mat

    ntg = ntheta/2
    
    ! ensure that gnew and source are initialized to zero
    gnew = 0.0 ; source = 0.0 ; gresponse = 0.0 ; source0 = 0.0 ; mu0_source = 0.0

    m_mat = 0.0
    allocate (p_mat(nresponse,nresponse)) ; p_mat = 0.0

    ! get response of particles with vpa=-dvpa below and including the midplane to a unit impulse
    ! in particles with vpa=0 below (not including) the midplane
    ! this will be used to calculate the distribution of particles with vpa=0 below the midplane, 
    ! which starts the implicit sweep in theta and vpa;
    do iglo = g_lo%llim_proc, g_lo%ulim_proc

       ik = ik_idx(g_lo,iglo)
       imu = imu_idx(g_lo,iglo)

! TMP FOR TESTING -- MAB       
!       if (imu==1) cycle ! treat mu=0 specially

       do it = 1, neigen(ik)

          call fill_response (iglo, it, ik, impulse=cmplx(1.0,0.0), response=gresponse(:,:,it,iglo))

          ! get the matrices multiplying g(theta<=0,vpa=0) and g(theta<=0,dvpa=-dvpa) 
          ! in the implicit solve. these are necessary to obtain the final response matrix

          ! p_mat is matrix multiplying g_{0} in equations for g(theta<0,vpa=0) (M1 from notes)
          ! m_mat is matrix multiplying g_{-1} in same equations (M2 from notes)

          idx = 1

          iseg = 1          
          do ig = ig_low(iseg), ig_mid(iseg)-1
             p_mat(idx,idx) = pmpfac(ig,-1,itmod(iseg,it,ik),iglo)
             p_mat(idx,idx+1) = pppfac(ig,-1,itmod(iseg,it,ik),iglo)
             m_mat(idx,idx,it,iglo) = pmmfac(ig,-1,itmod(iseg,it,ik),iglo)
             m_mat(idx,idx+1,it,iglo) = ppmfac(ig,-1,itmod(iseg,it,ik),iglo)
             idx = idx + 1
          end do

          ! RESPONSE CHANGE
          if (periodic(ik)) then
             p_mat(1,1) = globalfac1*p_mat(1,1)+globalfac2*ppmfac(ntgrid-1,0,itmod(nsegments(it,ik),it,ik),iglo)
             p_mat(1,2) = globalfac1*p_mat(1,2)
             m_mat(1,1:2,it,iglo) = globalfac1*m_mat(1,1:2,it,iglo)
          end if

          ! account for totally trapped particles, which decouple from the other theta-vpa points
          p_mat(idx,idx) = 1.0
          idx = idx + 1

          if (nsegments(it,ik) > 1) then
             do iseg = 2, nsegments(it,ik)
                do ig = ig_low(iseg)+1, ig_mid(iseg)-1
                   p_mat(idx,idx) = pmpfac(ig,-1,itmod(iseg,it,ik),iglo)
                   p_mat(idx,idx+1) = pppfac(ig,-1,itmod(iseg,it,ik),iglo)
                   m_mat(idx,idx,it,iglo) = pmmfac(ig,-1,itmod(iseg,it,ik),iglo)
                   m_mat(idx,idx+1,it,iglo) = ppmfac(ig,-1,itmod(iseg,it,ik),iglo)
                   idx = idx + 1
                end do
                ! account for totally trapped particles, which decouple from the other theta-vpa points
                p_mat(idx,idx) = 1.0
                
                idx = idx + 1
             end do
          end if
          
          ! account for vpa=vpa_min, which is determined by vpa BC
          p_mat(idx,idx) = 1.0
          p_mat(idx,idx+1) = -exp(energy(-ntgrid,-nvgrid+1,imu)-energy(-ntgrid,-nvgrid,imu))
          idx = idx + 1
          do iv = -nvgrid, -2
             if (iv /= -2) m_mat(idx,idx,it,iglo) = pppfac(-ntgrid,iv,itmod(1,it,ik),iglo)
             m_mat(idx,idx-1,it,iglo) = ppmfac(-ntgrid,iv,itmod(1,it,ik),iglo)
             p_mat(idx,idx) = pmpfac(-ntgrid,iv,itmod(1,it,ik),iglo)
             p_mat(idx,idx-1) = pmmfac(-ntgrid,iv,itmod(1,it,ik),iglo)
             idx = idx + 1
          end do

          ! RESPONSE CHANGE
          if (periodic(ik)) then
             ! this is entry from M_{21} in GS3 notes
             m_mat(idx-1,2,it,iglo) = pppfac(-ntgrid,-2,itmod(1,it,ik),iglo)
             ! this is entry from P_{13} in GS3 notes
             p_mat(1,idx) = globalfac2*pppfac(ntgrid-1,0,itmod(nsegments(it,ik),it,ik),iglo)
             ! these are entries from M_{13} in GS3 notes
             m_mat(1,idx-1,it,iglo) = globalfac2*pmmfac(ntgrid-1,0,itmod(nsegments(it,ik),it,ik),iglo)
             m_mat(1,idx,it,iglo) = globalfac2*pmpfac(ntgrid-1,0,itmod(nsegments(it,ik),it,ik),iglo)
          end if

!          do iv = 0, nvgrid-2
          do iv = 1, nvgrid-1
             p_mat(idx,idx) = ppmfac(ntgrid-1,iv,itmod(nsegments(it,ik),it,ik),iglo)
             p_mat(idx,idx+1) = pppfac(ntgrid-1,iv,itmod(nsegments(it,ik),it,ik),iglo)
             m_mat(idx,idx,it,iglo) = pmmfac(ntgrid-1,iv,itmod(nsegments(it,ik),it,ik),iglo)
             m_mat(idx,idx+1,it,iglo) = pmpfac(ntgrid-1,iv,itmod(nsegments(it,ik),it,ik),iglo)
             idx = idx + 1
          end do
          ! account for vpa=vpa_max, which is determined by vpa BC
          p_mat(idx,idx) = 1.0
          p_mat(idx,idx-1) = -exp(energy(ntgrid,nvgrid-1,imu)-energy(ntgrid,nvgrid,imu))
          idx = idx + 1

          iup = ir_up(it,ik)

          ! these correspond to (M1+M2*R1) and (M2*R2) from MAB notes, respectively
          gresponse(:iup,:iup,it,iglo) = p_mat(:iup,:iup) &
               + matmul(m_mat(:iup,:iup,it,iglo),gresponse(:iup,:iup,it,iglo))

!          write (*,*) 'entering invert_matrix'
          call invert_matrix (gresponse(:iup,:iup,it,iglo))

          p_mat = 0.0

       end do
    end do
    
    deallocate (p_mat)

  end subroutine get_gresponse_matrix

  subroutine fill_response (iglo, it, ik, impulse, response)

    use theta_grid, only: ntheta, ntgrid
    use dist_fn_arrays, only: gnew
    use vpamu_grids, only: nvgrid

    implicit none

    integer, intent (in) :: iglo, it, ik
    complex, intent (in) :: impulse
    complex, dimension (:,:), intent (out) :: response

    integer :: i, j, k, kmax, iseg, ig, iv

    ! initialize entries of response matrix to zero
    response = 0.0

    i = 1

    ! for iseg=1 (the left-most segment), g(theta=-pi,vpa=0) is not known
    ! for subsequent segments, g(theta=-pi,vpa=0) will have been solved for
    ! thus iseg=1 must be treated specially
    iseg = 1
    ! give a unit impulse to each theta from -pi up to but not including theta=0 (for vpa=0)
    ! and obtain the corresponding columns in the response matrix

! below should not be necessary because I am now enforcing periodicity in 
! implicit_sweep_right and implicit_sweep_left
!     ! RESPONSE CHANGE
!     if (periodic(ik)) then
!        ! must treat vpa=0, theta=theta_max specially in periodic case
!        ! in particular, it must be set equal to vpa=0, theta=theta_min
!        ! which is the first point to get an impulse below
!        gnew(ig_up(nsegments(it,ik)),0,itmod(nsegments(it,ik),it,ik),iglo) = impulse
!     end if

    do ig = ig_low(iseg), ig_mid(iseg)-1
       ! give a unit impulse to g at this theta
       gnew(ig,0,itmod(iseg,it,ik),iglo) = impulse

       ! sweep over entire (theta,vpa) plane and get values of gnew 
       ! for particles with vpa=-dvpa below the midplane
       call implicit_sweep (iglo, it, gnew)

       ! fill in cotribution from jth segment to ith column of initial 
       ! (ntheta/2*nsegment+1) x (ntheta/2*nsegment+1) x (nmu) response matrix          
       ! again, must treat left-most segment specially
!       ! RESPONSE CHANGE
!!       j = 1 ; k = 1 ; kmax = k+ntheta/2
!       j = 1 ; k = 1 ; kmax = k+ntheta/2-1
       j = 1 ; k = 1 ; kmax = k+ntheta/2

!       ! RESPONSE CHANGE
!       ! for leftmost theta of leftmost segment, use periodicity
!       iglow = ig_low(j)
!       if (periodic(ik)) iglow = ig_low(j)+1
!       ! RESPONSE CHANGE
!!       response(k:kmax,i) = gnew(ig_low(j):ig_mid(j),-1,itmod(j,it,ik),iglo)
!       response(k:kmax,i) = gnew(iglow:ig_mid(j),-1,itmod(j,it,ik),iglo)
       ! this is R_{11} from GS3 notes
       response(k:kmax,i) = gnew(ig_low(j):ig_mid(j),-1,itmod(j,it,ik),iglo)
       k = kmax + 1
       if (nsegments(it,ik) > 1) then
          do j = 2, nsegments(it,ik)
             ! this is R_{11} from GS3 notes
             kmax = k+ntheta/2-1
             response(k:kmax,i) = gnew(ig_low(j)+1:ig_mid(j),-1,itmod(j,it,ik),iglo)
             k = kmax + 1
          end do
       end if

!       ! if this ky has a periodic BC in theta, get response of g for vpa>0 at theta=theta_max-delthet
       ! if this ky has a periodic BC in theta, get response of g for vpa>=0 at theta=theta_max-delthet
       ! and for vpa<0 at theta=theta_min+delthet
       if (periodic(ik)) then

          ! RESPONSE CHANGE
!          kmax = k+nvgrid-1
!          response(k:kmax,i) = gnew(-ntgrid+1,-nvgrid:-1,itmod(1,it,ik),iglo)
          ! this is R_{31} from GS3 notes
          kmax = k+nvgrid-2
          response(k:kmax,i) = gnew(-ntgrid+1,-nvgrid:-2,itmod(1,it,ik),iglo)
          k = kmax + 1
          ! RESPONSE CHANGE
!          kmax = k+nvgrid-1
!          response(k:kmax,i) = gnew(ntgrid-1,1:nvgrid,itmod(nsegments(it,ik),it,ik),iglo)
          ! this is R_{21} from GS3 notes
          kmax = k+nvgrid
          response(k:kmax,i) = gnew(ntgrid-1,0:nvgrid,itmod(nsegments(it,ik),it,ik),iglo)

       end if

       do j = 1, nsegments(it,ik)
          ! reset g to zero everywhere
          gnew(:,:,itmod(j,it,ik),iglo) = 0.0
       end do

       i = i + 1
    end do
    
    ! put in blank column for (theta=0,vpa=0) since it can be determined without
    ! coupling to other theta-vpa points
    i = i + 1
    
    if (nsegments(it,ik) > 1) then
       do iseg = 2, nsegments(it,ik)
          do ig = ig_low(iseg)+1, ig_mid(iseg)-1
             ! give a unit impulse to g at this theta
             gnew(ig,0,itmod(iseg,it,ik),iglo) = impulse
             
             ! sweep over (theta,vpa) plane and get values of gnew for theta<0, vpa=-dvpa
             call implicit_sweep (iglo, it, gnew)
             
             ! fill in cotribution from jth segment to ith column of initial 
             ! (ntheta/2*nsegment+1) x (ntheta/2*nsegment+1) x (nmu) response matrix          
             ! again, must treat left-most segment specially
             j = 1 ; k = 1 ; kmax = k+ntheta/2
             response(k:kmax,i) = gnew(ig_low(j):ig_mid(j),-1,itmod(j,it,ik),iglo)
             k = kmax + 1
             if (nsegments(it,ik) > 1) then
                do j = 2, nsegments(it,ik)
                   kmax = k+ntheta/2-1
                   response(k:kmax,i) = gnew(ig_low(j)+1:ig_mid(j),-1,itmod(j,it,ik),iglo)
                   k = kmax + 1
                end do
             end if

             ! if this ky has a periodic BC in theta, get response of g for vpa>0 at theta=theta_max-delthet
             ! and for vpa<0 at theta=theta_min+delthet
             if (periodic(ik)) then
!                kmax = k+nvgrid-1
!                response(k:kmax,i,it,iglo) = gnew(ntgrid-1,1:nvgrid,itmod(nsegments(it,ik),it,ik),iglo)
!                k = kmax + 1
!                kmax = k+nvgrid-1
!                response(k:kmax,i,it,iglo) = gnew(-ntgrid+1,-1:-nvgrid:-1,itmod(1,it,ik),iglo)
                ! TMP FOR TESTING - MAB
!                response(k:k+1,i,it,iglo) = gnew(ntgrid-1,-nvgrid,itmod(nsegments(it,ik),it,ik),iglo)
!                k = kmax + 1

                ! RESPONSE CHANGE
!                kmax = k+nvgrid-1
!                response(k:kmax,i) = gnew(-ntgrid+1,-nvgrid:-1,itmod(1,it,ik),iglo)
                kmax = k+nvgrid-2
                response(k:kmax,i) = gnew(-ntgrid+1,-nvgrid:-2,itmod(1,it,ik),iglo)
                k = kmax + 1
                ! RESPONSE CHANGE
!                kmax = k+nvgrid-1
!                response(k:kmax,i) = gnew(ntgrid-1,1:nvgrid,itmod(nsegments(it,ik),it,ik),iglo)
                kmax = k+nvgrid
                response(k:kmax,i) = gnew(ntgrid-1,0:nvgrid,itmod(nsegments(it,ik),it,ik),iglo)
!                kmax = k+nvgrid
!                response(k:kmax,i) = gnew(ntgrid-1,0:nvgrid,itmod(nsegments(it,ik),it,ik),iglo)

                ! TMP FOR TESTING - MAB
!                k = kmax + 1
!                response(k:k+1,i,it,iglo) = gnew(-ntgrid+1,nvgrid,itmod(1,it,ik),iglo)
                
             end if
             
             do j = 1, nsegments(it,ik)
                ! reset g to zero everywhere
                gnew(:,:,itmod(j,it,ik),iglo) = 0.0
             end do
                
             i = i + 1
          end do
          ! put in blank column for theta=0 since it can be determined without
          ! coupling to other theta-vpa points
          i = i + 1
       end do
    end if

    ! if this ky is periodic, need to expand response matrix
    ! to get g(ntgrid-1,vpa>0) for unit impulses at g(-ntgrid,vpa>0)
    ! and g(-ntgrid+1,vpa<0) for unit impulses at g(ntgrid,vpa<0)
    if (periodic(ik)) then
!       ! for iseg=1 (the left-most segment), g(theta=-pi,vpa=0) is not known
!       ! for subsequent segments, g(theta=-pi,vpa=0) will have been solved for
!       ! thus iseg=1 must be treated specially
!       iseg = 1

       ! put blank column in for vpa=vpa_min, where vpa BC determines g
       i = i + 1

       ! for negative vpa, BC is at rightmost theta
       ! vpa BC at vpa_min determines g there
       ! TMP FOR TESTING -- MAB
!       do iv = -nvgrid, -1
       do iv = -nvgrid+1, -1
!       do iv = -1, -nvgrid+1, -1
          ! give a unit impulse to g at this vpa and leftmost theta
          gnew(ntgrid,iv,itmod(nsegments(it,ik),it,ik),iglo) = impulse
          ! TMP FOR TESTING -- MAB
          ! note that we are imposing vpa BC at rightmost theta and vpa_min
!          if (iv==-nvgrid+1 .and. (.not.vpa_bc_zero)) &
!               gnew(ntgrid,-nvgrid,itmod(nsegments(it,ik),it,ik),iglo) = impulse*decay_fac

          ! sweep over entire (theta,vpa) plane and get values of gnew
          ! for particles with vpa<0 at the leftmost theta
          call implicit_sweep (iglo, it, gnew)

          ! fill in cotribution from jth segment to ith column of initial 
          ! (ntheta/2*nsegment+1) x (ntheta/2*nsegment+1) x (nmu) response matrix          
          ! again, must treat left-most segment specially
          j = 1 ; k = 1 ; kmax = k+ntheta/2
          response(k:kmax,i) = gnew(ig_low(j):ig_mid(j),-1,itmod(j,it,ik),iglo)
          k = kmax + 1
          if (nsegments(it,ik) > 1) then
             do j = 2, nsegments(it,ik)
                kmax = k+ntheta/2-1
                response(k:kmax,i) = gnew(ig_low(j)+1:ig_mid(j),-1,itmod(j,it,ik),iglo)
                k = kmax + 1
             end do
          end if

          ! this ky has a periodic BC in theta, so get response of g for vpa>0 at theta=theta_max-delthet
          ! and for vpa<0 at theta=theta_min+delthet
          ! TMP FOR TESTING -- MAB
!          response(k:k+1,i,it,iglo) = gnew(ntgrid-1,-nvgrid,itmod(nsegments(it,ik),it,ik),iglo)
!          k = kmax + 1

          ! RESPONSE CHANGE
!          kmax = k+nvgrid-1
!          response(k:kmax,i) = gnew(-ntgrid+1,-nvgrid:-1,itmod(1,it,ik),iglo) 
          kmax = k+nvgrid-2
          response(k:kmax,i) = gnew(-ntgrid+1,-nvgrid:-2,itmod(1,it,ik),iglo)
          k = kmax + 1
          ! RESPONSE CHANGE
!          kmax = k+nvgrid-1
!          response(k:kmax,i) = gnew(ntgrid-1,1:nvgrid,itmod(nsegments(it,ik),it,ik),iglo)
          kmax = k+nvgrid
          response(k:kmax,i) = gnew(ntgrid-1,0:nvgrid,itmod(nsegments(it,ik),it,ik),iglo)

          ! TMP FOR TESTING -- MAB
!          k = kmax + 1
!          response(k:k+1,i,it,iglo) = gnew(-ntgrid+1,nvgrid,itmod(1,it,ik),iglo)

          do j = 1, nsegments(it,ik)
             ! reset g to zero everywhere
             gnew(:,:,itmod(j,it,ik),iglo) = 0.0
          end do

          i = i + 1
       end do
       
       ! for positive vpa, BC is at leftmost theta
       ! we know gnew(theta>0,vpa_max) via vpa BC
       ! TMP FOR TESTING -- MAB
       do iv = 1, nvgrid-1
!       do iv = 1, nvgrid
          ! give a unit impulse to g at this vpa and leftmost theta
          gnew(-ntgrid,iv,itmod(1,it,ik),iglo) = impulse
          ! TMP FOR TESTING -- MAB
          ! note that we are imposing vpa BC at leftmost theta and vpa_max
!          if (iv==nvgrid-1 .and. (.not.vpa_bc_zero)) gnew(-ntgrid,nvgrid,itmod(1,it,ik),iglo) = impulse*decay_fac

          ! sweep over entire (theta,vpa) plane and get values of gnew
          ! for particles with vpa>0 at the rightmost theta
          call implicit_sweep (iglo, it, gnew)

          ! fill in cotribution from jth segment to ith column of initial
          ! (ntheta/2*nsegment+1) x (ntheta/2*nsegment+1) x (nmu) response matrix          
          ! again, must treat left-most segment specially
          j = 1 ; k = 1 ; kmax = k+ntheta/2
          response(k:kmax,i) = gnew(ig_low(j):ig_mid(j),-1,itmod(j,it,ik),iglo)
          k = kmax + 1
          if (nsegments(it,ik) > 1) then
             do j = 2, nsegments(it,ik)
                kmax = k+ntheta/2-1
                response(k:kmax,i) = gnew(ig_low(j)+1:ig_mid(j),-1,itmod(j,it,ik),iglo)
                k = kmax + 1
             end do
          end if

          ! this ky has a periodic BC in theta, so get response of g for vpa>0 at theta=theta_max-delthet
          ! and for vpa<0 at theta=theta_min+delthet
          ! TMP FOR TESTING -- MAB
!          response(k:k+1,i,it,iglo) = gnew(ntgrid-1,-nvgrid,itmod(nsegments(it,ik),it,ik),iglo)
!          k = kmax + 1

          ! RESPONSE CHANGE
!          kmax = k+nvgrid-1
!          response(k:kmax,i) = gnew(-ntgrid+1,-nvgrid:-1,itmod(1,it,ik),iglo)
          kmax = k+nvgrid-2
          response(k:kmax,i) = gnew(-ntgrid+1,-nvgrid:-2,itmod(1,it,ik),iglo)
          k = kmax + 1
          ! RESPONSE CHANGE
!          kmax = k+nvgrid-1
!          response(k:kmax,i) = gnew(ntgrid-1,1:nvgrid,itmod(nsegments(it,ik),it,ik),iglo)
          kmax = k+nvgrid
          response(k:kmax,i) = gnew(ntgrid-1,0:nvgrid,itmod(nsegments(it,ik),it,ik),iglo)

          ! TMP FOR TESTING -- MAB
!          k = kmax + 1
!          response(k:k+1,i,it,iglo) = gnew(-ntgrid+1,nvgrid,itmod(1,it,ik),iglo)

          do j = 1, nsegments(it,ik)
             ! reset g to zero everywhere
             gnew(:,:,itmod(j,it,ik),iglo) = 0.0
          end do

          i = i + 1
       end do
       ! put blank column in for vpa=vpa_max, where vpa BC determines g
       i = i + 1


!        if (nsegments(it,ik) > 1) then
!           do iseg = 2, nsegments(it,ik)
!              do iv = 1, nvgrid
!                 ! give a unit impulse to g at this vpa
!                 gnew(-ntgrid,iv,itmod(1,it,ik),iglo) = impulse
                
!                 ! sweep over (theta,vpa) plane and get values of gnew for theta<0, vpa=-dvpa
!                 call implicit_sweep (iglo, it, gnew)
                
!                 ! fill in cotribution from jth segment to ith column of initial 
!                 ! (ntheta/2*nsegment+1) x (ntheta/2*nsegment+1) x (nmu) response matrix          
!                 ! again, must treat left-most segment specially
!                 j = 1 ; k = 1 ; kmax = k+ntheta/2
!                 response(k:kmax,i,it,iglo) = gnew(ig_low(j):ig_mid(j),-1,itmod(j,it,ik),iglo)
!                 k = kmax + 1
!                 if (nsegments(it,ik) > 1) then
!                    do j = 2, nsegments(it,ik)
!                       kmax = k+ntheta/2-1
!                       response(k:kmax,i,it,iglo) = gnew(ig_low(j)+1:ig_mid(j),-1,itmod(j,it,ik),iglo)
!                       k = kmax + 1
!                    end do
!                 end if
                
!                 ! if this ky has a periodic BC in theta, get response of g for vpa>0 at theta=theta_max-delthet
!                 ! and for vpa<0 at theta=theta_min+delthet
!                 kmax = k+nvgrid-1
!                 response(k:kmax,i,it,iglo) = gnew(ntgrid-1,1:nvgrid,itmod(nsegments(it,ik),it,ik),iglo)
!                 k = kmax + 1
!                 kmax = k+nvgrid-1
!                 response(k:kmax,i,it,iglo) = gnew(-ntgrid+1,-1:-nvgrid:-1,itmod(1,it,ik),iglo)
                
!                 do j = 1, nsegments(it,ik)
!                    ! reset g to zero everywhere
!                    gnew(:,:,itmod(j,it,ik),iglo) = 0.0
!                 end do
                
!                 i = i + 1
!              end do
! !             do iv = -nvgrid, -1
!              do iv = -1, -nvgrid, -1
!                 ! give a unit impulse to g at this vpa
!                 gnew(ntgrid,iv,itmod(nsegments(it,ik),it,ik),iglo) = impulse
                
!                 ! sweep over (theta,vpa) plane and get values of gnew for theta<0, vpa=-dvpa
!                 call implicit_sweep (iglo, it, gnew)
                
!                 ! fill in cotribution from jth segment to ith column of initial 
!                 ! (ntheta/2*nsegment+1) x (ntheta/2*nsegment+1) x (nmu) response matrix          
!                 ! again, must treat left-most segment specially
!                 j = 1 ; k = 1 ; kmax = k+ntheta/2
!                 response(k:kmax,i,it,iglo) = gnew(ig_low(j):ig_mid(j),-1,itmod(j,it,ik),iglo)
!                 k = kmax + 1
!                 if (nsegments(it,ik) > 1) then
!                    do j = 2, nsegments(it,ik)
!                       kmax = k+ntheta/2-1
!                       response(k:kmax,i,it,iglo) = gnew(ig_low(j)+1:ig_mid(j),-1,itmod(j,it,ik),iglo)
!                       k = kmax + 1
!                    end do
!                 end if
                
!                 ! if this ky has a periodic BC in theta, get response of g for vpa>0 at theta=theta_max-delthet
!                 ! and for vpa<0 at theta=theta_min+delthet
!                 kmax = k+nvgrid-1
!                 response(k:kmax,i,it,iglo) = gnew(ntgrid-1,1:nvgrid,itmod(nsegments(it,ik),it,ik),iglo)
!                 k = kmax + 1
!                 kmax = k+nvgrid-1
!                 response(k:kmax,i,it,iglo) = gnew(-ntgrid+1,-1:-nvgrid:-1,itmod(1,it,ik),iglo)
                
!                 do j = 1, nsegments(it,ik)
!                    ! reset g to zero everywhere
!                    gnew(:,:,itmod(j,it,ik),iglo) = 0.0
!                 end do
                
!                 i = i + 1
!              end do
!           end do
!        end if
    end if

  end subroutine fill_response
    
  ! implicit_sweep starts with a boundary condition for g^{n+1} along the line (theta<0,vpa=0)
  ! as well as zero BCs at (theta=-theta_max,vpa>0), (theta=theta_max,vpa<0), 
  ! (vpa=-vpa_max,theta<0), and (vpa=vpa_max,theta>0), and solves for g^{n+1}.  
  ! note that dvpa/dt<0 for theta>0 and dvpa/dt>0 for theta<0.
  subroutine implicit_sweep (iglo, it, gfnc, source_mod)

    use theta_grid, only: ntgrid, ntheta, theta, delthet
    use vpamu_grids, only: nvgrid
    use gs2_layouts, only: g_lo, ik_idx, imu_idx
    use dist_fn_arrays, only: source

    implicit none

    integer, intent (in) :: iglo, it
    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in out) :: gfnc
    complex, dimension (:), intent (in out), optional :: source_mod

    integer :: ntg, k, ik, iseg, kmax, iup, imu
    complex, dimension (:), allocatable :: dgdgn

    ik = ik_idx(g_lo,iglo)
    imu = imu_idx(g_lo,iglo)

    ntg = ntheta/2

    if (present(source_mod)) then
       allocate (dgdgn(nresponse)) ; dgdgn = 0.
    else
       allocate (dgdgn(1)) ; dgdgn = 0.
    end if

    ! reset gfnc to zero for all phase space points
    ! except those which are potentially determined already
    ! via parallel BC or response matrix
!    call reset_gfnc (gfnc, iglo, it)

    iseg = 1

    call implicit_sweep_right (iglo, itmod(iseg,it,ik), ik, iseg, imu, gfnc)

    if (nsegments(it,ik) > 1) then
       do iseg = 2, nsegments(it,ik)
          ! make connection with previous segment by setting gfnc(-pi) for this segment
          ! equal to gfnc(pi) from previous segment
          gfnc(ig_low(iseg),0:,itmod(iseg,it,ik),iglo) = gfnc(ig_up(iseg-1),0:,itmod(iseg-1,it,ik),iglo)
          call implicit_sweep_right (iglo, itmod(iseg,it,ik), ik, iseg, imu, gfnc)
       end do
    end if

    iseg = nsegments(it,ik)
    call implicit_sweep_left(iglo, itmod(iseg,it,ik), ik, iseg, gfnc)
    if (nsegments(it,ik) > 1) then

       do iseg = nsegments(it,ik)-1, 1, -1
          
          ! make connection with previous segment by setting gfnc(pi) for this segment
          ! equal to gfnc(-pi) from previous segment
          gfnc(ig_up(iseg),-nvgrid:-1,itmod(iseg,it,ik),iglo) &
               = gfnc(ig_low(iseg+1),-nvgrid:-1,itmod(iseg+1,it,ik),iglo)

          call implicit_sweep_left (iglo, itmod(iseg,it,ik), ik, iseg, gfnc)
       end do
    end if

    if (present(source_mod)) then

       k = 1

       iseg = 1 ; kmax = k+ntg
       dgdgn(k:kmax) = gfnc(ig_low(iseg):ig_mid(iseg),-1,itmod(iseg,it,ik),iglo)
       k = kmax+1

       if (nsegments(it,ik) > 1) then
          do iseg = 2, nsegments(it,ik)
             kmax = k+ntg-1
             dgdgn(k:kmax) = gfnc(ig_low(iseg)+1:ig_mid(iseg),-1,itmod(iseg,it,ik),iglo)
             k = kmax+1
          end do
       end if

       if (periodic(ik)) then
          ! RESPONSE CHANGE
!          kmax = k+nvgrid-1
!          dgdgn(k:kmax) = gfnc(-ntgrid+1,-nvgrid:-1,itmod(1,it,ik),iglo)
          kmax = k+nvgrid-2
          dgdgn(k:kmax) = gfnc(-ntgrid+1,-nvgrid:-2,itmod(1,it,ik),iglo)
          k = kmax+1
          ! RESPONSE CHANGE
!          kmax = k+nvgrid-1
!          dgdgn(k:kmax) = gfnc(ntgrid-1,1:nvgrid,itmod(nsegments(it,ik),it,ik),iglo)
          kmax = k+nvgrid
          dgdgn(k:kmax) = gfnc(ntgrid-1,0:nvgrid,itmod(nsegments(it,ik),it,ik),iglo)
       end if

       iup = ir_up(it,ik)

       ! this is M1*r_{-1} from MAB notes
       dgdgn(:iup) = matmul(m_mat(:iup,:iup,it,iglo),dgdgn(:iup))

       k=1
       iseg=1 ; kmax=k+ntg
       source_mod(k:kmax-1) = source(ig_low(iseg):ig_mid(iseg)-1,-1,itmod(iseg,it,ik),iglo) &
            - dgdgn(k:kmax-1)
       source_mod(kmax) = source0(iseg,itmod(iseg,it,ik),iglo)
       k = kmax+1

       if (nsegments(it,ik) > 1) then
          do iseg = 2, nsegments(it,ik)
             kmax = k+ntg-1
             source_mod(k:kmax-1) = source(ig_low(iseg)+1:ig_mid(iseg)-1,-1,itmod(iseg,it,ik),iglo) &
                  - dgdgn(k:kmax-1)
             source_mod(kmax) = source0(iseg,itmod(iseg,it,ik),iglo)
             k = kmax+1
          end do
       end if

       if (periodic(ik)) then
          kmax = k+nvgrid-1
!          source_mod(k:kmax) = source(-ntgrid+1,-nvgrid:-1,itmod(1,it,ik),iglo) &
          source_mod(k) = 0.
          source_mod(k+1:kmax) = source(-ntgrid,-nvgrid:-2,itmod(1,it,ik),iglo) &
               - dgdgn(k+1:kmax)
          k = kmax+1
          kmax = k+nvgrid-1
!          source_mod(k:kmax) = source(ntgrid-1,1:nvgrid,itmod(nsegments(it,ik),it,ik),iglo) &
          source_mod(k:kmax-1) = source(ntgrid-1,1:nvgrid-1,itmod(nsegments(it,ik),it,ik),iglo) &
               - dgdgn(k:kmax-1)
!          kmax = k+nvgrid
!!          source_mod(k:kmax) = source(ntgrid-1,1:nvgrid,itmod(nsegments(it,ik),it,ik),iglo) &
!          source_mod(k:kmax-1) = source(ntgrid-1,0:nvgrid-1,itmod(nsegments(it,ik),it,ik),iglo) &
!               - dgdgn(k:kmax-1)
          source_mod(kmax) = 0.
       end if

    end if

    ! set BCs at +/- theta_max
    if (.not. theta_bc_zero) then
       gfnc(-ntgrid,1:nvgrid,it,iglo) = &
            gfnc(-ntgrid+1,1:nvgrid,it,iglo)*exp(-delthet(-ntgrid)*(delthet(-ntgrid)-2.*theta(-ntgrid+1)))
       gfnc(ntgrid,-nvgrid:-1,it,iglo) &
            = gfnc(ntgrid-1,-nvgrid:-1,it,iglo)*exp(-delthet(ntgrid-1)*(delthet(ntgrid-1)+2.*theta(ntgrid-1)))
    end if

    ! TMP FOR TESTING -- MAB
!    if (imu==1) gfnc(:,:,:,iglo) = 0.0

    deallocate (dgdgn)

  end subroutine implicit_sweep

  ! implicit_sweep_right starts with a boundary condition for g^{n+1} along the line (theta<0,vpa=0)
  ! as well as a BC at (theta=-pi,vpa>0) and at (theta>0,vpa=vpa_max), and solves for g^{n+1}
  ! note that dvpa/dt<0 for theta>0 and dvpa/dt>0 for theta<0.
  subroutine implicit_sweep_right (iglo, it, ik, iseg, imu, gfnc)

    use gs2_layouts, only: ik_idx, is_idx, g_lo
    use theta_grid, only: ntgrid
    use vpamu_grids, only: nvgrid
    use dist_fn_arrays, only: source

    implicit none

    integer, intent (in) :: iglo, it, ik, iseg, imu
    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in out) :: gfnc

    integer :: ig, iv, iv_low, igup

    ! RESPONSE CHANGE
    ! do not solve for rightmost theta of rightmost segment if periodic
    ! instead is set equal to leftmost theta of leftmost segment earlier
    igup = ig_up(iseg)
    if (periodic(ik) .and. iseg==nsegments(it,ik)) then
       ! if periodic and rightmost segment, set g(theta_max,vpa>=0)=g(theta_min,vpa>=0)
       gfnc(igup,0:nvgrid,it,iglo) = gfnc(-ntgrid,0:nvgrid,it,iglo)
       ! since gfnc(igup) is fixed by periodicity, no need to solve for it later
       igup = igup-1
    end if

    ! treat mu=0 specially since it means different vpars are uncoupled
    if (imu == 1) then
       iv_low = 1

       ! first get gfnc(theta,vpa=0,mu=0), which is decoupled from other points
       gfnc(ig_low(iseg):ig_up(iseg),0,it,iglo) = mu0_source(ig_low(iseg):ig_up(iseg), &
            it,ik_idx(g_lo,iglo),is_idx(g_lo,iglo))
    else
       iv_low = 0

       ! first get gfnc(theta=vpa=0), which is decoupled from other points
       gfnc(ig_mid(iseg),0,it,iglo) = source0(iseg,it,iglo)
    end if

    ! boundary condition for vpa > 0 is g(theta=-infinity) = 0 or periodic
    ! either way, it has been obtained already
    ! boundary condition for theta < 0 (dvpa/dt > 0) is g(vpa=-infinity)=0 or decaying Maxwellian
    do iv = 0, nvgrid-1
       do ig = ig_low(iseg), ig_mid(iseg)-1
          gfnc(ig+1,iv+1,it,iglo) = (source(ig,iv,it,iglo) - (gfnc(ig,iv+1,it,iglo)*pmpfac(ig,iv,it,iglo) &
               + gfnc(ig+1,iv,it,iglo)*ppmfac(ig,iv,it,iglo) &
               + gfnc(ig,iv,it,iglo)*pmmfac(ig,iv,it,iglo))) / pppfac(ig,iv,it,iglo)
       end do
    end do

    ! note that for periodic BC do not need to solve for ig = ig_up at all

    ! boundary condition for vpa > 0 particles is g(theta=-infinity) = 0
    ! boundary condition for theta > 0 (corresponds to dvpa/dt < 0) is 
    ! g(vpa=infinity) = 0 or decaying Maxwellian
    iv = nvgrid-1
    do ig = ig_mid(iseg), igup-1
       gfnc(ig+1,iv,it,iglo) = (source(ig,iv,it,iglo) - (gfnc(ig+1,iv+1,it,iglo)*pppfac(ig,iv,it,iglo) &
            + gfnc(ig,iv+1,it,iglo)*pmpfac(ig,iv,it,iglo) &
            + gfnc(ig,iv,it,iglo)*pmmfac(ig,iv,it,iglo))) / ppmfac(ig,iv,it,iglo)
       ! might need to worry about periodicity for ig_up(nsegments(it,ik)),nvgrid,it,iglo
       ! deal with vpa = vpa_max BC for theta > 0
       gfnc(ig+1,nvgrid,it,iglo) = gfnc(ig+1,nvgrid-1,it,iglo)*decay_fac
    end do

    do iv = nvgrid-2, iv_low, -1
       do ig = ig_mid(iseg), igup-1
          gfnc(ig+1,iv,it,iglo) = (source(ig,iv,it,iglo) - (gfnc(ig+1,iv+1,it,iglo)*pppfac(ig,iv,it,iglo) &
               + gfnc(ig,iv+1,it,iglo)*pmpfac(ig,iv,it,iglo) &
               + gfnc(ig,iv,it,iglo)*pmmfac(ig,iv,it,iglo))) / ppmfac(ig,iv,it,iglo)
       end do
    end do

    ! TMP FOR TESTING -- MAB
    ! might need to worry about periodicity for ig_up(nsegments(it,ik)),nvgrid,it,iglo
    ! deal with vpa = vpa_max BC for theta > 0
!    gfnc(ig_mid(iseg)+1:ig_up(iseg),nvgrid,it,iglo) &
!         = gfnc(ig_mid(iseg)+1:ig_up(iseg),nvgrid-1,it,iglo)*decay_fac

  end subroutine implicit_sweep_right

  ! implicit_sweep starts with a boundary condition for g^{n+1} along the line (theta>0,vpa=0)
  ! as well as BCs at (theta=theta_max,vpa<0) and (theta<0,vpa=-vpa_max), and solves for g^{n+1}
  ! note that dvpa/dt<0 for theta>0 and dvpa/dt>0 for theta<0.
  subroutine implicit_sweep_left (iglo, it, ik, iseg, gfnc)

    use gs2_layouts, only: g_lo
    use theta_grid, only: ntgrid
    use vpamu_grids, only: nvgrid
    use dist_fn_arrays, only: source

    implicit none

    integer, intent (in) :: iglo, it, ik, iseg
    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in out) :: gfnc

    integer :: ig, iv, iglow

    iglow = ig_low(iseg)
    if (periodic(ik) .and. iseg==1) then
       ! if periodic and leftmost segment, set g(theta_min,vpa<0) = g(theta_max,vpa<0)
       gfnc(iglow,-nvgrid:-1,it,iglo) = gfnc(ntgrid,-nvgrid:-1,it,iglo)
       ! since g(theta_min,vpa<0) fixed by periodicity, no need to solve for it later
       iglow = iglow+1
    end if

    ! boundary condition for vpa < 0 particles is g(theta=infinity) = 0
    ! boundary condition for theta > 0 (dvpa/dt < 0) is g(vpa=infinity)=0
    do iv = -1, -nvgrid, -1
       do ig = ig_up(iseg)-1, ig_mid(iseg), -1
          gfnc(ig,iv,it,iglo) = (source(ig,iv,it,iglo) - (gfnc(ig+1,iv+1,it,iglo)*pppfac(ig,iv,it,iglo) &
               + gfnc(ig+1,iv,it,iglo)*ppmfac(ig,iv,it,iglo) &
               + gfnc(ig,iv+1,it,iglo)*pmpfac(ig,iv,it,iglo))) / pmmfac(ig,iv,it,iglo)
       end do
    end do
    
    ! boundary condition for vpa < 0 is g(theta=infinity) = 0
    ! boundary condition for theta < 0 (dvpa/dt > 0) is g(vpa=-infinity)=0
    iv = -nvgrid
    ! RESPONSE CHANGE
!    do ig = ig_mid(iseg)-1, ig_low(iseg), -1
    do ig = ig_mid(iseg)-1, iglow, -1
       gfnc(ig,iv+1,it,iglo) = (source(ig,iv,it,iglo) - (gfnc(ig+1,iv+1,it,iglo)*pppfac(ig,iv,it,iglo) &
            + gfnc(ig+1,iv,it,iglo)*ppmfac(ig,iv,it,iglo) &
            + gfnc(ig,iv,it,iglo)*pmmfac(ig,iv,it,iglo))) / pmpfac(ig,iv,it,iglo)
       ! may need to deal with ig_low(1),-nvgrid,it,iglo
       ! deal with vpa = -vpa_max BC for theta < 0
       gfnc(ig,-nvgrid,it,iglo) = gfnc(ig,-nvgrid+1,it,iglo)*decay_fac
   end do

    do iv = -nvgrid+1, -2
       ! RESPONSE CHANGE
!       do ig = ig_mid(iseg)-1, ig_low(iseg), -1
       do ig = ig_mid(iseg)-1, iglow, -1
          gfnc(ig,iv+1,it,iglo) = (source(ig,iv,it,iglo) - (gfnc(ig+1,iv+1,it,iglo)*pppfac(ig,iv,it,iglo) &
               + gfnc(ig+1,iv,it,iglo)*ppmfac(ig,iv,it,iglo) &
               + gfnc(ig,iv,it,iglo)*pmmfac(ig,iv,it,iglo))) / pmpfac(ig,iv,it,iglo)
       end do
    end do

    ! may need to deal with ig_low(1),-nvgrid,it,iglo
    ! deal with vpa = -vpa_max BC for theta < 0
!    gfnc(ig_low(iseg):ig_mid(iseg)-1,-nvgrid,it,iglo) &
!         = gfnc(ig_low(iseg):ig_mid(iseg)-1,-nvgrid+1,it,iglo)*decay_fac
    
  end subroutine implicit_sweep_left

  ! subroutine reset_gfnc (gfnc, iglo, it)
    
  !   use gs2_layouts, only: g_lo, ik_idx
  !   use theta_grid, only: ntgrid
  !   use vpamu_grids, only: nvgrid
    
  !   implicit none
    
  !   complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in out) :: gfnc
  !   integer, intent (in) :: iglo, it
    
  !   integer :: iseg, ik
  !   integer :: iglow, igup
    
  !   ik = ik_idx(g_lo,iglo)
    
  !   do iseg = 1, nsegments(it,ik)

  !      iglow = ig_low(iseg)
  !      igup = ig_up(iseg)
  !      if (periodic(ik)) then
  !         if (iseg==1) then
  !            ! do not zero out leftmost theta for vpa < 0 if periodic,
  !            ! as these points have been determined via response matrix
  !            iglow = ig_low(iseg)+1
  !         else if (iseg==nsegments(it,ik)) then
  !            ! do not zero out rightmost theta for vpa >= 0 if periodic,
  !            ! as these points have been determined via response matrix
  !            igup = ig_up(iseg)-1
  !         end if
  !      end if

  !      ! zero out vpa=0, theta above and including midplane
  !      gfnc(ig_mid(iseg):igup,0,itmod(iseg,it,ik),iglo) = 0.0
  !      ! zero out vpa > 0, all theta except leftmost theta in each segment
  !      ! note that leftmost theta in all segments except first will be 
  !      ! immediately set by linking to rightmost theta of previous segment
  !      ! if periodic, rightmost theta of last segment also not zeroed out
  !      ! as it will have been fixed already by periodicity
  !      gfnc(ig_low(iseg)+1:igup,1:,itmod(iseg,it,ik),iglo) = 0.0
  !      ! zero out vpa < 0, all tehta except rightmost theta in each segment
  !      ! note that rightmost theta in all segments except first will be
  !      ! immediately set by linking to leftmost theat of previous segment
  !      ! if periodic, leftmost theta of final segment also not zeroed out
  !      ! as it will have been fixed already by periodicity
  !      gfnc(iglow:ig_up(iseg)-1,-nvgrid:-1,itmod(iseg,it,ik),iglo) = 0.0
  !   end do

  ! end subroutine reset_gfnc

  subroutine implicit_solve (gfnc, gfncold, phi, phinew, &
       apar, aparnew, istep, equation)

    use gs2_layouts, only: g_lo, ik_idx, imu_idx
    use theta_grid, only: dbdthetc, ntgrid, ntheta
    use vpamu_grids, only: nvgrid
    use kt_grids, only: ntheta0

    implicit none

    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in out) :: gfnc, gfncold
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, phinew
    complex, dimension (-ntgrid:,:,:), intent (in) :: apar, aparnew
    integer, intent (in) :: istep, equation

    integer :: iglo, iseg, k, kmax, ntg, it, ik, iup, imu

    complex, dimension (:), allocatable :: source_mod, tmp

    ntg = ntheta/2

    if (.not. allocated(source_mod)) then
!       allocate (tmp(ntg*nseg_max+1)) ; tmp = 0.0
!       allocate (source_mod(ntg*nseg_max+1)) ; source_mod = 0.0
       allocate (tmp(nresponse)) ; tmp = 0.0
       allocate (source_mod(nresponse)) ; source_mod = 0.0
    end if
    
    call get_source (gfncold, phi, phinew, apar, aparnew, istep)
    if (equation==1) then
       call add_gprim_source
    else if (equation==2) then
       call add_higher_order_source
    end if
    
    call set_parallel_bc (gfnc)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do it = 1, ntheta0
          ! initially set g^{n+1}(theta<0,vpa=0) to zero
          ! also set g^{n+1}(theta_min,vpa>0)
          ! and g^{n+1}(theta_max,vpa<0) to zero
          where (dbdthetc(:ntgrid-1,1) < -epsilon(0.))
             gfnc(:ntgrid-1,0,it,iglo) = 0.0
          end where
       end do
    end do

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       
       ik = ik_idx(g_lo,iglo)
       imu = imu_idx(g_lo,iglo)
       
       do it = 1, neigen(ik)
          
          iup = ir_up(it,ik)
          
          ! treat mu=0 specially since different vpa points are not connected
! TMP FOR TESTING
!          if (imu /= 1) then
             ! sweep through vpa and theta once to obtain the response of 
             ! particles with vpa=-dvpa below and including the midplane
             ! as well as particles with (theta_max-dtheta,vpa>0) and (theta_min+dtheta,vpa<0),
             ! to g^{n} when g^{n+1} for particles with vpa=0 below and including the midplane is zero
             ! as is g^{n+1} for particles with (theta_min,vpa>0) and (theta_max,vpa<0)
             call implicit_sweep (iglo, it, gfnc, source_mod=source_mod)

             ! get g^{n+1}(theta<=0,vpa=0), g^{n+1}(theta_min,vpa>0), and g^{n+1}(theta_max,vpa<0) using response matrix
             tmp(:iup) = matmul(gresponse(:iup,:iup,it,iglo), source_mod(:iup))

             k=1
             iseg=1 ; kmax=k+ntg
             gfnc(ig_low(iseg):ig_mid(iseg),0,itmod(iseg,it,ik),iglo) = tmp(k:kmax)
             
             k = kmax+1
             if (nsegments(it,ik) > 1) then
                do iseg = 2, nsegments(it,ik)
                   kmax = k+ntg-1
                   ! set g(-pi) for this segment equal to g(pi) at the previous segment
                   gfnc(ig_low(iseg),0,itmod(iseg,it,ik),iglo) &
                        = gfnc(ig_up(iseg-1),0,itmod(iseg-1,it,ik),iglo)
                   ! fill in the remaining theta values in this segment
                   gfnc(ig_low(iseg)+1:ig_mid(iseg),0,itmod(iseg,it,ik),iglo) = tmp(k:kmax)
                   k = kmax+1
                end do
             end if

             if (periodic(ik)) then
                kmax = k+nvgrid-1
                gfnc(ntgrid,-nvgrid:-1,itmod(nsegments(it,ik),it,ik),iglo) = tmp(k:kmax)
                gfnc(ntgrid,-nvgrid,itmod(nsegments(it,ik),it,ik),iglo) &
                     = gfnc(ntgrid,-nvgrid+1,itmod(nsegments(it,ik),it,ik),iglo)*decay_fac
                k = kmax+1
                kmax = k+nvgrid-1
                gfnc(-ntgrid,1:nvgrid,itmod(1,it,ik),iglo) = tmp(k:kmax)
!                kmax = k+nvgrid
!                gfnc(-ntgrid,0:nvgrid,itmod(1,it,ik),iglo) = tmp(k:kmax)
                gfnc(-ntgrid,nvgrid,itmod(1,it,ik),iglo) = gfnc(-ntgrid,nvgrid-1,itmod(1,it,ik),iglo)*decay_fac
             end if

!          end if

          ! with g^{n+1}(theta=0,vpa>0) specified, sweep once more to get g^{n+1} everywhere else
          call implicit_sweep (iglo, it, gfnc)

       end do
    end do

    if (allocated (tmp)) deallocate (source_mod, tmp)

  end subroutine implicit_solve

  subroutine get_source (gfnc, phifnc, phinewfnc, aparfnc, aparnewfnc, istep)

    use constants, only: zi
    use centering, only: get_cell_value
    use dist_fn_arrays, only: aj0, vpar, source
    use gs2_time, only: code_dt
    use species, only: spec
    use run_parameters, only: t_imp, fapar
    use theta_grid, only: ntgrid, thet_imp
    use vpamu_grids, only: nvgrid, anon, vpac, anonc, vpa_imp
    use gs2_layouts, only: g_lo, ik_idx, imu_idx, is_idx
    use nonlinear_terms, only: nonlin
    use kt_grids, only: ntheta0

    implicit none

    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in) :: gfnc
    complex, dimension (-ntgrid:,:,:), intent (in) :: phifnc, phinewfnc
    complex, dimension (-ntgrid:,:,:), intent (in) :: aparfnc, aparnewfnc
    integer, intent (in) :: istep

    integer :: ig, iv, iglo, iseg, it, ik, idx, imu, is
    complex, dimension (:), allocatable :: phi_m
    complex, dimension (:,:), allocatable :: phic
    complex, dimension (:,:), allocatable :: apar_m
    complex, dimension (:,:), allocatable :: aparc

    allocate (phic(-ntgrid:ntgrid,3))
    allocate (phi_m(-ntgrid:ntgrid))

    allocate (aparc(-ntgrid:ntgrid,2))
    allocate (apar_m(-ntgrid:ntgrid,2))

    if (fapar < epsilon(0.0)) then
       aparc = 0.
       apar_m = 0.
    end if

    ! note that GKE is evaluated at cell values
    ! ig=-ntgrid is first cell value, ig=ntgrid-1 is last cell value
    ! iv=-nvgrid is first cell value, iv=nvgrid-1 is last cell value

    do iglo = g_lo%llim_proc, g_lo%ulim_proc

       ik = ik_idx(g_lo,iglo)
       imu = imu_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)

!       do ie = 1, neigen(ik)

!          do iseg = 1, nsegments(ie,ik)
!             it = itmod(iseg,ie,ik)
       do it = 1, ntheta0
          ! get space-time cell values for J0*phi.  where to evaluate in 
          ! spatial cell depends on vpa in order to get upwinding right
          phic(:,3) = aj0(:,it,iglo)*(t_imp*phinewfnc(:,it,ik) + (1.0-t_imp)*phifnc(:,it,ik))
          call get_cell_value (1.0-thet_imp, phic(:,3), phic(:,1), -ntgrid)
          call get_cell_value (thet_imp, phic(:,3), phic(:,2), -ntgrid)

          ! d(J0*phi)/dtheta = del(J0*phi) / delthet...(1/delthet) is already present in vpar,
          ! which multiplies phi_m below
          phi_m(:ntgrid-1) = phic(-ntgrid+1:,3) - phic(:ntgrid-1,3)

          if (fapar > epsilon(0.0)) then
             aparc(:,2) = spec(is)%stm*aj0(:,it,iglo)*(t_imp*aparnewfnc(:,it,ik) + (1.0-t_imp)*aparfnc(:,it,ik))
             call get_cell_value (1.0-thet_imp, aparc(:,2), aparc(:,1), -ntgrid)
             call get_cell_value (thet_imp, aparc(:,2), aparc(:,2), -ntgrid)

             apar_m(:,2) = spec(is)%zstm*aj0(:,it,iglo)*(aparnewfnc(:,it,ik)-aparfnc(:,it,ik))
             call get_cell_value (1.0-thet_imp, apar_m(:,2), apar_m(:,1), -ntgrid)
             call get_cell_value (thet_imp, apar_m(:,2), apar_m(:,2), -ntgrid)

             ! nvgrid index never used later so serves as convenient dummy index
!             aparc(:,nvgrid) = aj0(:,it,iglo)*(t_imp*aparnewfnc(:,it,ik) + (1.0-t_imp)*aparfnc(:,it,ik))
!             call get_cell_value (1.0-thet_imp, aparc(:,nvgrid), aparc(:,-nvgrid), -ntgrid)
!             call get_cell_value (thet_imp, aparc(:,nvgrid), aparc(:,0), -ntgrid)

!             aparc(:,-nvgrid+1:-1) = spread(aparc(:,-nvgrid),2,nvgrid-1)
!             aparc(:,1:nvgrid-1) = spread(aparc(:,0),2,nvgrid-1)
!
!             apar_m(:ntgrid-1,:) = vpc*spread(aparc(-ntgrid+1:,nvgrid)-aparc(:ntgrid-1,nvgrid),2,nvpa)

!             aparc(-ntgrid:ntgrid-1,:) = vpc*aparc(-ntgrid:ntgrid-1,:)
             
             ! TMP FOR TESTING -- MAB
!             apar_m = 0.
!             aparc = 0.
          end if

          do iv = -nvgrid, nvgrid-1
             ! idx needed to distinguish between +/- vpar for upwinding
             if (iv < 0) then
                idx = 1
             else
                idx = 2
             end if
             
!             do ig = ig_low(iseg), ig_up(iseg)-1
             do ig = -ntgrid, -1
                source(ig,iv,it,iglo) = -(gfnc(ig+1,iv+1,it,iglo)*mppfac(ig,iv,it,iglo) &
                     + gfnc(ig+1,iv,it,iglo)*mpmfac(ig,iv,it,iglo) &
                     + gfnc(ig,iv+1,it,iglo)*mmpfac(ig,iv,it,iglo) &
                     + gfnc(ig,iv,it,iglo)*mmmfac(ig,iv,it,iglo)) &
                     ! from here on are actual source terms (not contributions from g at old time level)
                     - anonc(ig,iv,imu)*(vpar(ig,iv,is)*phi_m(ig) &
                     + zi*(wdriftc(ig,iv,it,iglo)*phic(ig,idx) &
                     - wstarc(ig,iv,iglo)*(phic(ig,idx)-vpac(iv,1)*aparc(ig,1))) &
                     + vpac(iv,1)*apar_m(ig,1))
             end do
             do ig = 0, ntgrid-1
                source(ig,iv,it,iglo) = -(gfnc(ig+1,iv+1,it,iglo)*mppfac(ig,iv,it,iglo) &
                     + gfnc(ig+1,iv,it,iglo)*mpmfac(ig,iv,it,iglo) &
                     + gfnc(ig,iv+1,it,iglo)*mmpfac(ig,iv,it,iglo) &
                     + gfnc(ig,iv,it,iglo)*mmmfac(ig,iv,it,iglo)) &
                     ! from here on are actual source terms (not contributions from g at old time level)
                     - anonc(ig,iv,imu)*(vpar(ig,iv,is)*phi_m(ig) &
                     + zi*(wdriftc(ig,iv,it,iglo)*phic(ig,idx) &
                     - wstarc(ig,iv,iglo)*(phic(ig,idx)-vpac(iv,2)*aparc(ig,2))) &
                     + vpac(iv,2)*apar_m(ig,2))
             end do
          end do

          if (nonlin) then
             select case (istep)
             case (0)
                ! do nothing
             case (1)
                do ig = -ntgrid, ntgrid-1
!                   source(ig,:,it,iglo) = source(ig,:,it,iglo) + 0.5*code_dt*gexp_1(ig,:,it,iglo)
                   source(ig,:,it,iglo) = source(ig,:,it,iglo) + code_dt*gexp_1(ig,:,it,iglo)
                end do
             case (2)
                do ig = -ntgrid, ntgrid-1
                   source(ig,:,it,iglo) = source(ig,:,it,iglo) + code_dt*( &
!                   source(ig,:,it,iglo) = source(ig,:,it,iglo) + 0.5*code_dt*( &
                        1.5*gexp_1(ig,:,it,iglo) - 0.5*gexp_2(ig,:,it,iglo))
                end do
             case default
                do ig = -ntgrid, ntgrid-1
                   source(ig,:,it,iglo) = source(ig,:,it,iglo) + code_dt*( &
!                   source(ig,:,it,iglo) = source(ig,:,it,iglo) + 0.5*code_dt*( &
                        (23./12.)*gexp_1(ig,:,it,iglo) &
                        - (4./3.)*gexp_2(ig,:,it,iglo) &
                        + (5./12.)*gexp_3(ig,:,it,iglo))
                end do
             end select
          end if

          if (periodic(ik)) then
             source(-ntgrid,-1,it,iglo) = globalfac1*source(-ntgrid,-1,it,iglo) &
                  + globalfac2*source(ntgrid-1,0,it,iglo)
          end if

          ! treat mu=0,vpa=0 points specially, as they are decoupled from other points
          if (imu==1) then
             
             mu0_source(:,it,ik,is) = gfnc(:,0,it,iglo) + zi*anon(:,0,imu) &
                  * (wstar(:,0,iglo)-wdrift(:,0,it,iglo))*phic(:,3)

             if (nonlin) then
                select case (istep)
                case (0)
                   ! do nothing
                case (1)
!                   do ig = -ntgrid, ntgrid-1
                   mu0_source(:,it,ik,is) = mu0_source(:,it,ik,is) + code_dt*gexp_1(:,0,it,iglo)
!                   mu0_source(:,it,ik,is) = mu0_source(:,it,ik,is) + 0.5*code_dt*gexp_1(:,0,it,iglo)
!                   end do
                case (2)
!                   do ig = -ntgrid, ntgrid-1
                   mu0_source(:,it,ik,is) = mu0_source(:,it,ik,is) + code_dt*( &
!                   mu0_source(:,it,ik,is) = mu0_source(:,it,ik,is) + 0.5*code_dt*( &
                        1.5*gexp_1(:,0,it,iglo) - 0.5*gexp_2(:,0,it,iglo))
!                   end do
                case default
!                   do ig = -ntgrid, ntgrid-1
                   mu0_source(:,it,ik,is) = mu0_source(:,it,ik,is) + code_dt*( &
!                   mu0_source(:,it,ik,is) = mu0_source(:,it,ik,is) + 0.5*code_dt*( &
                        (23./12.)*gexp_1(:,0,it,iglo) &
                        - (4./3.)*gexp_2(:,0,it,iglo) &
                        + (5./12.)*gexp_3(:,0,it,iglo))
!                   end do
                end select
             end if
                
          end if

          if (it > neigen(ik)) cycle

          do iseg = 1, nsegments(it,ik)
             ! this is the source at the (theta=0,vpa=0) grid point (not a cell value)
             ! it is needed because that point does not take info from other grid points
             ! except indirectly through the fields
             source0(iseg,itmod(iseg,it,ik),iglo) = gfnc(ig_mid(iseg),0,itmod(iseg,it,ik),iglo) + zi*anon(ig_mid(iseg),0,imu) &
                  * (wstar(ig_mid(iseg),0,iglo)-wdrift(ig_mid(iseg),0,itmod(iseg,it,ik),iglo)) &
                  * aj0(ig_mid(iseg),itmod(iseg,it,ik),iglo)*(t_imp*phinewfnc(ig_mid(iseg),itmod(iseg,it,ik),ik) &
                  + (1.0-t_imp)*phifnc(ig_mid(iseg),itmod(iseg,it,ik),ik))
             
             if (nonlin) then
                select case (istep)
                case (0)
                   ! do nothing
                case (1)
                   source0(iseg,itmod(iseg,it,ik),iglo) = source0(iseg,itmod(iseg,it,ik),iglo) + code_dt*gexp_1(ig_mid(iseg),0,itmod(iseg,it,ik),iglo)
!                   source0(iseg,itmod(iseg,it,ik),iglo) = source0(iseg,itmod(iseg,it,ik),iglo) + 0.5*code_dt*gexp_1(ig_mid(iseg),0,itmod(iseg,it,ik),iglo)
                case (2)
                   source0(iseg,itmod(iseg,it,ik),iglo) = source0(iseg,itmod(iseg,it,ik),iglo) + code_dt*( &
!                   source0(iseg,itmod(iseg,it,ik),iglo) = source0(iseg,itmod(iseg,it,ik),iglo) + 0.5*code_dt*( &
                        1.5*gexp_1(ig_mid(iseg),0,itmod(iseg,it,ik),iglo) - 0.5*gexp_2(ig_mid(iseg),0,itmod(iseg,it,ik),iglo))
                case default
                   source0(iseg,itmod(iseg,it,ik),iglo) = source0(iseg,itmod(iseg,it,ik),iglo) + code_dt*( &
!                   source0(iseg,itmod(iseg,it,ik),iglo) = source0(iseg,itmod(iseg,it,ik),iglo) + 0.5*code_dt*( &
                        (23./12.)*gexp_1(ig_mid(iseg),0,itmod(iseg,it,ik),iglo) &
                        - (4./3.)*gexp_2(ig_mid(iseg),0,itmod(iseg,it,ik),iglo) &
                        + (5./12.)*gexp_3(ig_mid(iseg),0,itmod(iseg,it,ik),iglo))
                end select
             end if
          end do
          
       end do
    end do
    
    deallocate (phic, phi_m, aparc, apar_m)

  end subroutine get_source

  subroutine add_gprim_source

    use constants, only: zi
    use gs2_layouts, only: g_lo, is_idx, imu_idx, ik_idx
    use fields_arrays, only: phinew, phi
    use dist_fn_arrays, only: source, gnew, gold, vparp, aj0, vpar, aj0p
    use species, only: spec
    use run_parameters, only: t_imp
    use centering, only: get_cell_value
    use theta_grid, only: ntgrid, dbdthet, thet_imp
    use vpamu_grids, only: nvgrid, vpa_imp, anonc, anon
    use kt_grids, only: ntheta0

    implicit none

    integer :: ig, iglo, is, iv, imu, it, ik
    integer :: iseg
    complex :: phitc
    complex, dimension (:,:), allocatable :: gtc, gc, gtvc, g_m, dgdv, dgdvc
    complex, dimension (:,:), allocatable :: phic, jpphic
    complex, dimension (:), allocatable :: phi_m, jpphi_m

    allocate (gtc(-ntgrid:ntgrid,-nvgrid:nvgrid))
    allocate (gc(-ntgrid:ntgrid,-nvgrid:nvgrid))
    allocate (gtvc(-ntgrid:ntgrid,-nvgrid:nvgrid))
    allocate (g_m(-ntgrid:ntgrid,-nvgrid:nvgrid))
    allocate (dgdv(-ntgrid:ntgrid,-nvgrid:nvgrid))
    allocate (dgdvc(-ntgrid:ntgrid,-nvgrid:nvgrid))

    allocate (phic(-ntgrid:ntgrid,3), jpphic(-ntgrid:ntgrid,3))
    allocate (phi_m(-ntgrid:ntgrid), jpphi_m(-ntgrid:ntgrid))

    ! add in rhs of gprim equation here to usual source
    do iglo = g_lo%llim_proc, g_lo%ulim_proc

       ik = ik_idx(g_lo,iglo)
       imu = imu_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)

       do it = 1, ntheta0

          ! get space-time cell values for J0*phi and dJ0/dr*phi.  where to evaluate in 
          ! spatial cell depends on vpa in order to get upwinding right
          phic(:,3) = aj0(:,it,iglo)*(t_imp*phinew(:,it,ik) + (1.0-t_imp)*phi(:,it,ik))
          call get_cell_value (1.0-thet_imp, phic(:,3), phic(:,1), -ntgrid)
          call get_cell_value (thet_imp, phic(:,3), phic(:,2), -ntgrid)
          jpphic(:,3) = aj0p(:,it,iglo)*(t_imp*phinew(:,it,ik) + (1.0-t_imp)*phi(:,it,ik))
          call get_cell_value (1.0-thet_imp, jpphic(:,3), jpphic(:,1), -ntgrid)
          call get_cell_value (thet_imp, jpphic(:,3), jpphic(:,2), -ntgrid)

          ! d(J0*phi)/dtheta = del(J0*phi) / delthet...(1/delthet) is already present in vpar,
          ! which multiplies phi_m below
          phi_m(:ntgrid-1) = phic(-ntgrid+1:,3) - phic(:ntgrid-1,3)
          jpphi_m(:ntgrid-1) = jpphic(-ntgrid+1:,3) - jpphic(:ntgrid-1,3)

          ! get time-centered values for g
          gtc = t_imp*gnew(:,:,it,iglo) + (1.0-t_imp)*gold(:,:,it,iglo)

          ! get time, vpa, and theta cell values for g
          call get_cell_value (thet_imp, vpa_imp, gtc(:,:), gc(:,:), -ntgrid, -nvgrid)

          ! get time and v-space cell values for g
          do ig = -ntgrid, ntgrid
             if (dbdthet(ig) < 0.0) then
                call get_cell_value (vpa_imp, gtc(ig,:), gtvc(ig,:), -nvgrid)
             else
                call get_cell_value (1.0-vpa_imp, gtc(ig,:), gtvc(ig,:), -nvgrid)
             end if
          end do

          ! get theta derivative of gtvc with respect to theta
          g_m(:ntgrid-1,:nvgrid-1) = gtvc(-ntgrid+1:,:nvgrid-1) - gtvc(:ntgrid-1,:nvgrid-1)

          ! get vpa derivative of gtc
          dgdv(:,:nvgrid-1) = gtc(:,-nvgrid+1:)-gtc(:,:nvgrid-1)
       
          ! get time and theta cell values for g
          do iv = -nvgrid, -1
             call get_cell_value (1.0-thet_imp, dgdv(:,iv), dgdvc(:,iv), -ntgrid)
          end do
          do iv = 0, nvgrid-1
             call get_cell_value (thet_imp, dgdv(:,iv), dgdvc(:,iv), -ntgrid)
          end do

          ! these are the terms appearing on the RHS of the d/dt(dg/dr) equation
          source(:ntgrid-1,-nvgrid:-1,it,iglo) = source(:ntgrid-1,-nvgrid:-1,it,iglo) &
               ! vparp is parallel streaming term with b.grad(theta) -> d(b.grad(theta))/dr
               ! mirror is the mirror term (dg/dvpa) with b.grad(B) -> d(b.grad(B))/dr
               ! varfacc term is the parallel electric field term with F_M/T -> d[F_M/T]/dr
               ! streamfac term is the parallel electric field term with b.grad[theta] -> d[b.grad[theta]]/dr
               ! wdriftpc is the wdrift term with v_M.k -> d[v_M.k]/dr
               ! wdriftc term is the usual wdrift with F_M/T -> d[F_M/T]/dr
               ! wstarpc is the wstar term with wstar -> d(wstar)/dr
               - vparp(:ntgrid-1,-nvgrid:-1,is)*g_m(:ntgrid-1,-nvgrid:-1) &
               - mirror(:ntgrid-1,-nvgrid:-1,imu,is)*dgdvc(:ntgrid-1,-nvgrid:-1) &
               - zi*wdriftpc(:ntgrid-1,-nvgrid:-1,it,iglo)*spec(is)%tz*gc(:ntgrid-1,-nvgrid:-1) &
               - anonc(:ntgrid-1,-nvgrid:-1,imu) & ! <-- begin anonc multiplication
               * ( vpar(:ntgrid-1,-nvgrid:-1,is) & ! <-- begin vpar multiplication
               * (spread(phi_m(:ntgrid-1),2,nvgrid)*(varfacc(:ntgrid-1,-nvgrid:-1,imu,is) &
               + spread(streamfac(:ntgrid-1,1),2,nvgrid)) &
               + spread(jpphi_m(:ntgrid-1),2,nvgrid) ) & ! <-- end vpar multiplication
               + zi*(spread(phic(:ntgrid-1,1),2,nvgrid)*( wdriftpc(:ntgrid-1,-nvgrid:-1,it,iglo) &
               + wdriftc(:ntgrid-1,-nvgrid:-1,it,iglo)*varfacc(:ntgrid-1,-nvgrid:-1,imu,is) &
               - wstarpc(:ntgrid-1,-nvgrid:-1,iglo)) &
               + (wdriftc(:ntgrid-1,-nvgrid:-1,it,iglo)-wstarc(:ntgrid-1,-nvgrid:-1,iglo)) &
               * spread(jpphic(:ntgrid-1,1),2,nvgrid)) ) ! <-- end anonc multiplication
          source(:ntgrid-1,0:nvgrid-1,it,iglo) = source(:ntgrid-1,0:nvgrid-1,it,iglo) &
               - vparp(:ntgrid-1,0:nvgrid-1,is)*g_m(:ntgrid-1,0:nvgrid-1) &
               - mirror(:ntgrid-1,0:nvgrid-1,imu,is)*dgdvc(:ntgrid-1,0:nvgrid-1) &
               - zi*wdriftpc(:ntgrid-1,0:nvgrid-1,it,iglo)*spec(is)%tz*gc(:ntgrid-1,0:nvgrid-1) &
               - anonc(:ntgrid-1,0:nvgrid-1,imu) & ! <-- begin anonc multiplication
               * ( vpar(:ntgrid-1,0:nvgrid-1,is) & ! <-- begin vpar multiplication
               * (spread(phi_m(:ntgrid-1),2,nvgrid)*(varfacc(:ntgrid-1,0:nvgrid-1,imu,is) &
               + spread(streamfac(:ntgrid-1,2),2,nvgrid)) &
               + spread(jpphi_m(:ntgrid-1),2,nvgrid) ) & ! <-- end vpar multiplication
               + zi*(spread(phic(:ntgrid-1,2),2,nvgrid)*( wdriftpc(:ntgrid-1,0:nvgrid-1,it,iglo) &
               + wdriftc(:ntgrid-1,0:nvgrid-1,it,iglo)*varfacc(:ntgrid-1,0:nvgrid-1,imu,is) &
               - wstarpc(:ntgrid-1,0:nvgrid-1,iglo)) &
               + (wdriftc(:ntgrid-1,0:nvgrid-1,it,iglo)-wstarc(:ntgrid-1,0:nvgrid-1,iglo)) &
               * spread(jpphic(:ntgrid-1,2),2,nvgrid)) ) ! <-- end anonc multiplication
          
          ! treat mu=0,vpa=0 points specially, as they are decoupled from other points
          if (imu==1) then
             mu0_source(:,it,ik,is) = mu0_source(:,it,ik,is) &
                  + zi &  ! <-- begin zi multiplication
                  * ( -wdriftp(:,0,it,iglo)*spec(is)%tz*gtc(:,0) & ! <--remove )
                  + anon(:,0,imu) & ! <-- begin anon multiplication
                  * ( phic(:,3)*(wstarp(:,0,iglo)-wdriftp(:,0,it,iglo) &
                  - wdrift(:,0,it,iglo)*varfac(:,0,imu,is)) &
                  - jpphic(:,3)*(wdrift(:,0,it,iglo)-wstar(:,0,iglo)) )) ! <-- end zi,anon multiplications
          end if
          
          if (it > neigen(ik)) cycle
          do iseg = 1, nsegments(it,ik)
             phitc = t_imp*phinew(ig_mid(iseg),itmod(iseg,it,ik),ik) + (1.0-t_imp)*phi(ig_mid(iseg),itmod(iseg,it,ik),ik)
             ! this is the source at the (theta=0,vpa=0) grid point (not a cell value)
             ! it is needed because that point does not take info from other grid points
             ! except indirectly through the fields
             source0(iseg,itmod(iseg,it,ik),iglo) = source0(iseg,itmod(iseg,it,ik),iglo) &
                  - zi*wdriftp(ig_mid(iseg),0,itmod(iseg,it,ik),iglo)*spec(is)%tz*gtc(ig_mid(iseg),0) &
                  + zi*anon(ig_mid(iseg),0,imu) &
                  * aj0(ig_mid(iseg),itmod(iseg,it,ik),iglo)*phitc &
                  * ( -wdriftp(ig_mid(iseg),0,itmod(iseg,it,ik),iglo) + wstarp(ig_mid(iseg),0,iglo) &
                  - wdrift(ig_mid(iseg),0,itmod(iseg,it,ik),iglo)*varfac(ig_mid(iseg),0,imu,is)) &
                  - zi*anon(ig_mid(iseg),0,imu)*(wdrift(ig_mid(iseg),0,itmod(iseg,it,ik),iglo)-wstar(ig_mid(iseg),0,iglo)) &
                  * aj0p(ig_mid(iseg),itmod(iseg,it,ik),iglo)*phitc
          end do
          
       end do
    end do

    deallocate (gtc, gc, gtvc, g_m, dgdv, dgdvc, phic, phi_m, jpphic, jpphi_m)

  end subroutine add_gprim_source

  subroutine add_higher_order_source

    use centering, only: get_cell_value
    use fields_arrays, only: phi, phinew, phip, phipnew
    use dist_fn_arrays, only: aj0, source, gpnew, aj0p, gpold
    use species, only: spec
    use run_parameters, only: t_imp
    use theta_grid, only: ntgrid, thet_imp
    use vpamu_grids, only: nvgrid, anon, anonc, vpa_imp
    use gs2_layouts, only: g_lo, ik_idx, imu_idx, is_idx
    use kt_grids, only: ntheta0

    implicit none

    integer :: ig, iv, iglo, iseg, it, ik, idx, imu, is
    complex, dimension (:,:), allocatable :: phipc, gpc, gptc

    allocate (phipc(-ntgrid:ntgrid,3))
    allocate (gptc(-ntgrid:ntgrid,-nvgrid:nvgrid), gpc(-ntgrid:ntgrid,-nvgrid:nvgrid))

    ! note that GKE is evaluated at cell values
    ! ig=-ntgrid is first cell value, ig=ntgrid-1 is last cell value
    ! iv=-nvgrid is first cell value, iv=nvgrid-1 is last cell value

    do iglo = g_lo%llim_proc, g_lo%ulim_proc

       ik = ik_idx(g_lo,iglo)
       imu = imu_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)

       do it = 1, ntheta0

          ! this is time-cell value for d(J0*phi)/drho and dg/drho
          phipc(:,3) = aj0(:,it,iglo)*(t_imp*phipnew(:,it,ik) + (1.0-t_imp)*phip(:,it,ik)) &
               + aj0p(:,it,iglo)*(t_imp*phinew(:,it,ik) + (1.0-t_imp)*phi(:,it,ik))
          gptc = t_imp*gpnew(:,:,it,iglo) + (1.0-t_imp)*gpold(:,:,it,iglo)

          ! get theta-cell value for d(J0*phi)/drho
          call get_cell_value (1.0-thet_imp, phipc(:,3), phipc(:,1), -ntgrid)
          call get_cell_value (thet_imp, phipc(:,3), phipc(:,2), -ntgrid)

          ! get time, vpa, and theta cell values for g
          call get_cell_value (thet_imp, vpa_imp, gptc, gpc, -ntgrid, -nvgrid)

          do iv = -nvgrid, nvgrid-1
             ! idx needed to distinguish between +/- vpar for upwinding
             if (iv < 0) then
                idx = 1
             else
                idx = 2
             end if

             do iseg = 1, nsegments(it,ik)
                do ig = ig_low(iseg), ig_up(iseg)-1
                   source(ig,iv,it,iglo) = source(ig,iv,it,iglo) &
                        - wdriftmodc(ig,iv,imu)*(gpc(ig,iv)*spec(is)%tz &
                        + anonc(ig,iv,imu)*phipc(ig,idx))
                end do
             end do
          end do
          
          ! treat mu=0,vpa=0 points specially, as they are decoupled from other points
          if (imu==1) then
             mu0_source(:,it,ik,is) = mu0_source(:,it,ik,is) &
                  - wdriftmod(:,0,imu)*(gptc(:,0)*spec(is)%tz &
                  + anon(:,0,imu)*phipc(:,3))
          end if
          
          if (it > neigen(ik)) cycle
          do iseg = 1, nsegments(it,ik)
             ! this is the source at the (theta=0,vpa=0) grid point (not a cell value)
             ! it is needed because that point does not take info from other grid points
             ! except indirectly through the fields
             source0(iseg,itmod(iseg,it,ik),iglo) = source0(iseg,itmod(iseg,it,ik),iglo) &
                  - wdriftmod(ig_mid(iseg),0,imu)*(gptc(ig_mid(iseg),0)*spec(is)%tz &
                  + anon(ig_mid(iseg),0,imu)*phipc(ig_mid(iseg),3))
          end do

       end do
    end do

    deallocate (phipc, gpc, gptc)

  end subroutine add_higher_order_source

  subroutine set_parallel_bc (gfnc)

    use theta_grid, only: ntgrid
    use vpamu_grids, only: nvgrid
    use gs2_layouts, only: g_lo, ik_idx

    implicit none

    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in out) :: gfnc

    integer :: iglo, it, ik
    integer :: itleft, itright

    do iglo = g_lo%llim_proc, g_lo%ulim_proc

       ik = ik_idx(g_lo,iglo)
       do it = 1, neigen(ik)
          itleft = itmod(1,it,ik)
          itright = itmod(nsegments(it,ik),it,ik)

          gfnc(-ntgrid,1:nvgrid,itleft,iglo) = 0.0
          gfnc(ntgrid,-nvgrid:-1,itright,iglo) = 0.0
       end do

    end do

  end subroutine set_parallel_bc

!   subroutine get_total_g

!     use fields_arrays, only: phinew, phipnew
!     use dist_fn_arrays, only: gnew, gpnew, aj0, aj0p
!     use species, only: spec
!     use theta_grid, only: ntgrid
!     use vpamu_grids, only: nvgrid, anon
!     use gs2_layouts, only: g_lo, ik_idx, imu_idx, is_idx

!     implicit none

!     integer :: iv, iglo, it, ik, imu, is
!     complex, dimension (:), allocatable :: jphip

!     allocate (jphip(-ntgrid:ntgrid))

!     do iglo = g_lo%llim_proc, g_lo%ulim_proc

!        ik = ik_idx(g_lo,iglo)
!        imu = imu_idx(g_lo,iglo)
!        is = is_idx(g_lo,iglo)

!        jphip = aj0(:,it,iglo)*phipnew(:,it,ik) + aj0p(:,it,iglo)*phinew(:,it,ik)

! ! TMP FOR TESTING -- MAB
!        do iv = -nvgrid, nvgrid
!           gnew(:,iv,iglo) = gnew(:,iv,iglo) &
!                - wdriftmod(:,iv,imu) &
!                * (gpnew(:,iv,iglo)*spec(is)%tz &
!                + anon(:,iv,imu)*jphip)
!        end do
       
!     end do

!     deallocate (jphip)

!   end subroutine get_total_g

  subroutine allocate_arrays

    use kt_grids, only: naky, ntheta0, box
    use theta_grid, only: ntgrid, shat
    use dist_fn_arrays, only: g, gnew, gold, source
    use dist_fn_arrays, only: gpnew, gpold!, ghnew, ghold
    use dist_fn_arrays, only: kx_shift, theta0_shift   ! MR
    use gs2_layouts, only: g_lo
    use nonlinear_terms, only: nonlin
    use vpamu_grids, only: nvgrid

    implicit none

    if (.not. allocated(g)) then
       allocate (g    (-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (gold (-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (gnew (-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (g0   (-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (source(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (gpold(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (gpnew(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
!       allocate (ghold(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
!       allocate (ghnew(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
! #ifdef LOWFLOW
!        allocate (gexp_1(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
!        allocate (gexp_2(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
!        allocate (gexp_3(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc))
!        gexp_1 = 0. ; gexp_2 = 0. ; gexp_3 = 0.
! #else
       if (nonlin) then
          allocate (gexp_1(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
          allocate (gexp_2(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
          allocate (gexp_3(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
          gexp_1 = 0. ; gexp_2 = 0. ; gexp_3 = 0.
       else
          allocate (gexp_1(1,1,1,1), gexp_2(1,1,1,1), gexp_3(1,1,1,1))
       end if
! #endif
!       if (boundary_option_switch == boundary_option_linked) then
!          allocate (g_h(-ntgrid:ntgrid,-nvgrid:nvgrid,ntheta0,g_lo%llim_proc:g_lo%ulim_alloc))
!          g_h = 0.
!          allocate (save_h(2,g_lo%llim_proc:g_lo%ulim_alloc))
!          save_h = .false.
!       endif
       if (abs(g_exb*g_exbfac) > epsilon(0.)) then           ! MR 
          if (box .or. abs(shat) < epsilon(0.0)) then
             allocate (kx_shift(naky))
             kx_shift = 0.
          else
             allocate (theta0_shift(naky))
             theta0_shift = 0.
          endif
       endif                           ! MR end
    endif

    g = 0. ; gnew = 0. ; g0 = 0. ; gold = 0.
    gpnew = 0. ; gpold = 0. !; ghnew = 0. ; ghold = 0.

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

  subroutine timeadv (gfnc, gfncold, phi, apar, bpar, phinew, aparnew, istep, equation, mode)

    use theta_grid, only: ntgrid
!    use le_derivatives, only: vspace_derivatives
    use nonlinear_terms, only: add_explicit_terms
    use hyper, only: hyper_diff
    use run_parameters, only: nstep
    use vpamu_grids, only: nvgrid
    use gs2_layouts, only: g_lo

    implicit none

    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in out) :: gfnc, gfncold
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phinew, aparnew
    integer, intent (in) :: istep, equation
    integer, optional, intent (in) :: mode
    integer :: modep

    modep = 0
    if (present(mode)) modep = mode

    ! Calculate the explicit nonlinear terms
    call add_explicit_terms (gexp_1, gexp_2, gexp_3, &
         phi, apar, bpar, istep)
    ! Implicit Solve for gfnc
    call implicit_solve (gfnc, gfncold, phi, phinew, &
         apar, aparnew, istep, equation)
    ! Add hyper terms (damping)
    call hyper_diff (gfnc, phinew)
    ! Add collisions
!    call vspace_derivatives (gnew, g, g0, phi, apar, bpar, phinew, aparnew, bparnew, modep)

  end subroutine timeadv

  subroutine get_fieldcorrection (gfnc, phifnc, phipfnc)

    use dist_fn_arrays, only: aj0p, aj0
    use gs2_layouts, only: g_lo
    use species, only: spec, nspec
    use theta_grid, only: ntgrid, bmag, dBdrho
    use vpamu_grids, only: nvgrid, integrate_species
    use kt_grids, only: ntheta0, naky

    implicit none

    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in) :: gfnc
    complex, dimension (-ntgrid:,:,:), intent (in) :: phifnc
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phipfnc

    integer :: iglo, iv, it
    real, dimension (nspec) :: wgt
    complex, dimension (:,:,:), allocatable :: tot

    allocate (tot(-ntgrid:ntgrid,ntheta0,naky))

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do it = 1, ntheta0
          do iv = -nvgrid, nvgrid
             g0(:,iv,it,iglo) = (aj0p(:,it,iglo)+dBdrho*aj0(:,it,iglo)/bmag)*gfnc(:,iv,it,iglo)
          end do
       end do
    end do
    wgt = spec%z
    call integrate_species (g0, wgt, tot)
    where (abs(gamtot) > epsilon(0.))
       phipfnc = phipfnc + tot/gamtot + gamtotp*phifnc
    end where

    deallocate (tot)

  end subroutine get_fieldcorrection


  ! communication initializations for exb_shear should be done once and 
  ! redistribute routines should be used  BD

!   subroutine exb_shear (g0, phi, apar, bpar)
! ! MR, 2007: modified Bill Dorland version to include grids where kx grid
! !           is split over different processors
! ! MR, March 2009: ExB shear now available on extended theta grid (ballooning)
! ! CMR, May 2009: 2pishat correction factor on extended theta grid (ballooning)
! !                so GEXB is same physical quantity in box and ballooning
! ! CMR, Oct 2010: multiply timestep by tunits(iky) for runs with wstar_units=.t.
! ! CMR, Oct 2010: add save statements to prevent potentially large and memory 
! !                killing array allocations!
    
!     use mp, only: iproc, proc0, send, receive, mp_abort
!     use gs2_layouts, only: ik_idx, g_lo, idx_local, idx, proc_id
!     use run_parameters, only: tunits
!     use theta_grid, only: ntgrid, ntheta, shat
!     use file_utils, only: error_unit
!     use kt_grids, only: akx, aky, naky, ikx, ntheta0, box, theta0
!     use vpamu_grids, only: nmu, nvgrid
!     use species, only: nspec
!     use run_parameters, only: fphi, fapar, fbpar
!     use dist_fn_arrays, only: kx_shift, theta0_shift
!     use gs2_time, only: code_dt, code_dt_old
!     use constants, only: twopi    

!     complex, dimension (-ntgrid:,:,:), intent (in out) :: phi,    apar,    bpar
!     complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in out) :: g0
!     complex, dimension(:,:,:), allocatable :: temp 
!     complex, dimension(:,:), allocatable :: temp2
!     integer, dimension(1), save :: itmin
!     integer :: ierr, j 
!     integer :: ik, it, is, imu, iv, to_iglo, from_iglo
!     integer:: iib, iit, ileft, iright, i

!     real, save :: dkx, dtheta0
!     real :: gdt
!     logical, save :: exb_first = .true.
!     complex , dimension(-ntgrid:ntgrid) :: z
!     character(130) :: str

!     ierr = error_unit()

! ! MR, March 2009: remove ExB restriction to box grids
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! MR, 2007: Works for ALL layouts (some layouts need no communication)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ! Initialize kx_shift, jump, idx_indexed
!     if (exb_first) then
!        exb_first = .false.
!        allocate (jump(naky)) 
!        jump = 0
!        if (box .or. shat .eq. 0.0) then
!           allocate (ikx_indexed(ntheta0))
!           itmin = minloc (ikx)
          
!           do it=itmin(1), ntheta0
!              ikx_indexed (it+1-itmin(1)) = it
!           end do
          
!           do it=1,itmin(1)-1
!              ikx_indexed (ntheta0 - itmin(1) + 1 + it)= it
!           end do

!           if (ntheta0 .gt. 1) then 
!              dkx = akx(2)-akx(1)
!           else
!              write(ierr,*) "exb_shear: ERROR, need ntheta0>1 for sheared flow"
!           endif
!        else
! ! MR, March 2009: on extended theta grid theta0_shift tracks ExB shear
! ! CMR, 25 May 2009: fix 2pi error so that: dtheta0/dt = -GEXB/shat
!           if (ntheta0 .gt. 1) then 
!              dtheta0 = theta0(2,1)-theta0(1,1) 
!           else 
!              dtheta0=twopi
!           endif
!        end if
!     end if
    
    
!     ! BD: To do: Put the right timestep in here.
!     ! For now, approximate dt == 1/2 (t_(n+1) - t_(n-1))
!     ! with code_dt.  
!     !
!     ! Note: at first time step, there is a difference of a factor of 2.
!     !
 
!     ! necessary to get factor of 2 right in first time step and
!     ! also to get things right after changing time step size
!     ! added May 18, 2009 -- MAB
!     gdt = 0.5*(code_dt + code_dt_old)
    
! ! kx_shift is a function of time.   Update it here:  
! ! MR, 2007: kx_shift array gets saved in restart file
! ! CMR, 5/10/2010: multiply timestep by tunits(ik) for wstar_units=.true. 
!     if ( box .or. shat .eq. 0.0 ) then
!        do ik=1, naky
!           kx_shift(ik) = kx_shift(ik) - aky(ik)*g_exb*g_exbfac*gdt*tunits(ik)
!           jump(ik) = nint(kx_shift(ik)/dkx)
!           kx_shift(ik) = kx_shift(ik) - jump(ik)*dkx
!        end do
!     else
!        do ik=1, naky
!           theta0_shift(ik) = theta0_shift(ik) - g_exb*g_exbfac*gdt/shat*tunits(ik)
!           jump(ik) = nint(theta0_shift(ik)/dtheta0)
!           theta0_shift(ik) = theta0_shift(ik) - dtheta0*jump(ik)
!        enddo 
!     end if

    
!     if (.not. box .and. shat .ne. 0.0 ) then
! ! MR, March 2009: impact of ExB shear on extended theta grid computed here
! !                 for finite shat
!        do ik =1,naky
!           j=jump(ik)
!           if (j .eq. 0) cycle     
!           if (abs(j) .ge. ntheta0) then
!               write(str,fmt='("in exb_shear: jump(ik)=",i4," > ntheta0 =",i4," for ik=",i4". => reduce timestep or increase ntheta0")') j,ik,ntheta0
!               write(ierr,*) str
!               call mp_abort(str)
!           endif 
!           allocate(temp2(-ntgrid:ntgrid,abs(j)),temp(-ntgrid:ntgrid,-nvgrid:nvgrid,abs(j)))
!           iit=ntheta0+1-abs(j) ; iib=abs(j)
!           ileft = -ntgrid+ntheta ; iright=ntgrid-ntheta

!           if (fphi > epsilon(0.0)) then
!              if (j < 0) then
!                 temp2 = phi(:,:iib,ik)
!                 do i=1,iit-1
!                    phi(:,i,ik) = phi(:,i-j,ik)
!                 enddo
!                 phi(ileft:,iit:,ik) = temp2(:iright,:)
!                 phi(:ileft-1,iit:,ik) = 0.0
!              else 
!                 temp2 = phi(:,iit:,ik)
!                 do i=ntheta0,iib+1,-1 
!                    phi(:,i,ik) = phi(:,i-j,ik)
!                 enddo
!                 phi(:iright,:iib ,ik) = temp2(ileft:,:)
!                 phi(iright+1:ntgrid,:iib,:) = 0.0
!              endif
!           endif
!           if (fapar > epsilon(0.0)) then
!              if (j < 0) then
!                 temp2 = apar(:,:iib,ik)
!                 do i=1,iit-1
!                    apar(:,i,ik) = apar(:,i-j,ik)
!                 enddo
!                 apar(ileft:,iit:,ik) = temp2(:iright,:)
!                 apar(:ileft-1,iit:,ik) = 0.0
!              else 
!                 temp2 = apar(:,iit:,ik)
!                 do i=ntheta0,iib+1,-1 
!                    apar(:,i,ik) = apar(:,i-j,ik)
!                 enddo
!                 apar(:iright,:iib ,ik) = temp2(ileft:,:)
!                 apar(iright+1:ntgrid,:iib,:) = 0.0
!              endif
!           endif
!           if (fbpar > epsilon(0.0)) then
!              if (j < 0) then
!                 temp2 = bpar(:,:iib,ik)
!                 do i=1,iit-1
!                    bpar(:,i,ik) = bpar(:,i-j,ik)
!                 enddo
!                 bpar(ileft:,iit:,ik) = temp2(:iright,:)
!                 bpar(:ileft-1,iit:,ik) = 0.0
!              else 
!                 temp2 = bpar(:,iit:,ik)
!                 do i=ntheta0,iib+1,-1 
!                    bpar(:,i,ik) = bpar(:,i-j,ik)
!                 enddo
!                 bpar(:iright,:iib ,ik) = temp2(ileft:,:)
!                 bpar(iright+1:ntgrid,:iib,:) = 0.0
!              endif
!           end if

! ! now the distribution functions

!           do is=1,nspec
!              do imu=1,nmu

!                 if (j < 0) then
!                    do it = 1, iib
!                       from_iglo = idx(g_lo, ik, it, imu, is)
!                       if (idx_local (g_lo, from_iglo)) temp(:,:,it) = g0(:,:,from_iglo)
!                    end do

!                    do it = 1, iit-1                        
!                       to_iglo = idx(g_lo, ik, it, imu, is)
!                       from_iglo = idx(g_lo, ik, it-j, imu, is)

!                       if (idx_local(g_lo, to_iglo).and. idx_local(g_lo, from_iglo)) then
!                          g0(:,:,to_iglo) = g0(:,:,from_iglo)
!                       else if (idx_local(g_lo, from_iglo)) then
!                          do iv = -nvgrid, nvgrid
!                             call send(g0(:, iv, from_iglo), proc_id (g_lo, to_iglo))
!                          enddo
!                       else if (idx_local(g_lo, to_iglo)) then
!                          do iv = -nvgrid, nvgrid
!                             call receive(g0(:, iv, to_iglo), proc_id (g_lo, from_iglo))
!                          enddo
!                       endif
!                    enddo

!                    do it = iit, ntheta0                     
!                       to_iglo = idx(g_lo, ik, it, imu, is)
!                       from_iglo = idx(g_lo, ik, it-j-ntheta0, imu, is)
                      
!                       if (idx_local(g_lo, to_iglo) .and. idx_local(g_lo, from_iglo)) then
!                          g0(ileft:,:,to_iglo) = temp(:iright,:,it-iit+1)
!                          g0(:ileft-1,:,to_iglo) = 0.0
!                       else if (idx_local(g_lo, from_iglo)) then
!                          do iv = -nvgrid, nvgrid
!                             call send(temp(:, iv, it-iit), proc_id (g_lo, to_iglo))
!                          enddo
!                       else if (idx_local(g_lo, to_iglo)) then
!                          do iv = -nvgrid, nvgrid
!                             call receive(z, proc_id (g_lo, from_iglo))
!                             g0(ileft:,iv,to_iglo) = z(:iright)
!                             g0(:ileft-1,iv,to_iglo) = 0.0
!                          enddo
!                       endif
!                    enddo

!                 else ! j>0

!                    do it = 1, j
!                       from_iglo = idx(g_lo, ik, iit+it-1, imu, is)
!                       if (idx_local (g_lo, from_iglo)) temp(:,:,it) = g0(:,:,from_iglo)
!                    end do

!                    do it = ntheta0, iib+1, -1
!                       to_iglo = idx(g_lo, ik, it, imu, is)
!                       from_iglo = idx(g_lo, ik, it-j, imu, is)

!                       if (idx_local(g_lo, to_iglo) .and. idx_local(g_lo, from_iglo)) then
!                          g0(:,:,to_iglo) = g0(:,:,from_iglo)
!                       else if (idx_local(g_lo, from_iglo)) then
!                          do iv = -nvgrid, nvgrid
!                             call send(g0(:, iv, from_iglo), proc_id (g_lo, to_iglo))
!                          enddo
!                       else if (idx_local(g_lo, to_iglo)) then
!                          do iv = -nvgrid, nvgrid
!                             call receive(g0(:,iv,to_iglo), proc_id (g_lo, from_iglo))
!                          enddo
!                       endif
!                    enddo

!                    do it = 1, iib
!                       to_iglo = idx(g_lo, ik, it, imu, is)
!                       from_iglo = idx(g_lo, ik, iit+it-1, imu, is)

!                       if (idx_local(g_lo, to_iglo).and. idx_local(g_lo, from_iglo)) then
!                          g0(:iright,:,to_iglo) = temp(ileft:,:,it)
!                          g0(iright+1:,:,to_iglo) = 0.0
!                       else if (idx_local(g_lo, from_iglo)) then
!                          do iv = -nvgrid, nvgrid
!                             call send(temp(:, iv, it), proc_id (g_lo, to_iglo))
!                          enddo
!                       else if (idx_local(g_lo, to_iglo)) then
!                          do iv = -nvgrid, nvgrid
!                             call receive(z, proc_id (g_lo, from_iglo))
!                             g0(:iright,iv,to_iglo) = z(ileft:)
!                             g0(iright+1:,iv,to_iglo) = 0.0
!                          enddo
!                       endif
!                    enddo
!                 endif
!              enddo
!           enddo
!           deallocate (temp,temp2)
!        enddo
!     end if
    
!     if (box .or. shat .eq. 0.0) then
!        do ik = naky, 2, -1
!           if (jump(ik) < 0) then
!              if (fphi > epsilon(0.0)) then
!                 do it = 1, ntheta0 + jump(ik)
!                    phi(:,ikx_indexed(it),ik) = phi(:,ikx_indexed(it-jump(ik)),ik)
!                 end do
!                 do it = ntheta0 + jump(ik) + 1, ntheta0
!                    phi(:,ikx_indexed(it),ik) = 0.
!                 end do
!              end if
!              if (fapar > epsilon(0.0)) then
!                 do it = 1, ntheta0 + jump(ik)
!                    apar(:,ikx_indexed(it),ik) = apar(:,ikx_indexed(it-jump(ik)),ik)
!                 end do
!                 do it = ntheta0 + jump(ik) + 1, ntheta0
!                    apar (:,ikx_indexed(it),ik) = 0.
!                 end do
!              end if
!              if (fbpar > epsilon(0.0)) then 
!                 do it = 1, ntheta0 + jump(ik)
!                    bpar(:,ikx_indexed(it),ik) = bpar(:,ikx_indexed(it-jump(ik)),ik)
!                 end do
!                 do it = ntheta0 + jump(ik) + 1, ntheta0
!                    bpar (:,ikx_indexed(it),ik) = 0.
!                 end do
!              end if
!              do is=1,nspec
!                 do imu=1,nmu

!                    do it = 1, ntheta0 + jump(ik)                        

!                       to_iglo = idx(g_lo, ik, ikx_indexed(it), imu, is)
!                       from_iglo = idx(g_lo, ik, ikx_indexed(it-jump(ik)), imu, is)

!                       if (idx_local(g_lo, to_iglo) .and. idx_local(g_lo, from_iglo)) then
!                          g0(:,:,to_iglo) = g0(:,:,from_iglo)
!                       else if (idx_local(g_lo, from_iglo)) then
!                          do iv=-nvgrid,nvgrid
!                             call send (g0(:, iv, from_iglo), proc_id (g_lo, to_iglo))
!                          enddo
!                       else if (idx_local(g_lo, to_iglo)) then
!                          do iv=-nvgrid,nvgrid
!                             call receive (g0(:, iv, to_iglo), proc_id (g_lo, from_iglo))
!                          enddo
!                       endif
!                    enddo

!                    do it = ntheta0 + jump(ik) + 1, ntheta0                     
!                       to_iglo = idx(g_lo, ik, ikx_indexed(it), imu, is)
!                       if (idx_local (g_lo, to_iglo)) g0(:,:,to_iglo) = 0.
!                    enddo

!                 enddo
!              enddo
!           endif

!           if (jump(ik) > 0) then 
!              if (fphi > epsilon(0.0)) then
!                 do it = ntheta0, 1+jump(ik), -1
!                    phi(:,ikx_indexed(it),ik) = phi(:,ikx_indexed(it-jump(ik)),ik)
!                 end do
!                 do it = jump(ik), 1, -1
!                    phi(:,ikx_indexed(it),ik) = 0.
!                 end do
!              end if
!              if (fapar > epsilon(0.0)) then
!                 do it = ntheta0, 1+jump(ik), -1
!                    apar(:,ikx_indexed(it),ik) = apar(:,ikx_indexed(it-jump(ik)),ik)
!                 end do
!                 do it = jump(ik), 1, -1
!                    apar(:,ikx_indexed(it),ik) = 0.
!                 end do
!              end if
!              if (fbpar > epsilon(0.0)) then
!                 do it = ntheta0, 1+jump(ik), -1
!                    bpar(:,ikx_indexed(it),ik) = bpar(:,ikx_indexed(it-jump(ik)),ik)
!                 end do
!                 do it = jump(ik), 1, -1
!                    bpar(:,ikx_indexed(it),ik) = 0.
!                 end do
!              end if
!              do is=1,nspec
!                 do imu=1,nmu

!                    do it = ntheta0, 1+jump(ik), -1
                      
!                       to_iglo = idx(g_lo, ik, ikx_indexed(it), imu, is)
!                       from_iglo = idx(g_lo, ik, ikx_indexed(it-jump(ik)), imu, is)

!                       if (idx_local(g_lo, to_iglo) .and. idx_local(g_lo, from_iglo)) then
!                          g0(:,:,to_iglo) = g0(:,:,from_iglo)
!                       else if (idx_local(g_lo, from_iglo)) then
!                          do iv=-nvgrid,nvgrid
!                             call send(g0(:, iv, from_iglo), proc_id(g_lo, to_iglo))
!                          enddo
!                       else if (idx_local(g_lo, to_iglo)) then
!                          do iv=-nvgrid,nvgrid
!                             call receive(g0(:, iv, to_iglo), proc_id (g_lo, from_iglo))
!                          enddo
!                       endif
!                    enddo

!                    do it = jump(ik), 1, -1
!                       to_iglo = idx(g_lo, ik, ikx_indexed(it), imu, is)
!                       if (idx_local (g_lo, to_iglo)) g0(:,:,to_iglo) = 0.
!                    enddo

!                 enddo
!              enddo
!           endif
!        enddo
!     end if
!   end subroutine exb_shear

  subroutine getan (gfnc, antot, antota, antotp)

    use dist_fn_arrays, only: aj0, aj1, kperp2
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use vpamu_grids, only: vperp2, vpa, integrate_species, nvgrid
    use run_parameters, only: beta, fphi, fapar, fbpar
    use gs2_layouts, only: g_lo, imu_idx, ik_idx, is_idx
    use kt_grids, only: ntheta0

    implicit none

    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in) :: gfnc
    complex, dimension (-ntgrid:,:,:), intent (out) :: antot, antota, antotp
    real, dimension (nspec) :: wgt

    integer :: iv, iglo, ig, imu, it

    if (fphi > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do it = 1, ntheta0
             do iv = -nvgrid, nvgrid
                do ig=-ntgrid, ntgrid
                   g0(ig,iv,it,iglo) = aj0(ig,it,iglo)*gfnc(ig,iv,it,iglo)
                end do
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
          do it = 1, ntheta0
             do iv = -nvgrid, nvgrid
                do ig=-ntgrid, ntgrid
                   g0(ig,iv,it,iglo) = aj0(ig,it,iglo)*vpa(iv)*gfnc(ig,iv,it,iglo)
                end do
             end do
          end do
       end do
       
       wgt = 2.0*beta*spec%z*spec%dens*spec%stm
       call integrate_species (g0, wgt, antota)
    else
       antota=0.
    end if

    if (fbpar > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          imu = imu_idx(g_lo,iglo)
          do iv = -nvgrid, nvgrid
             do ig=-ntgrid, ntgrid
                g0(ig,iv,:,iglo) = aj1(ig,:,iglo)*vperp2(ig,imu)*gfnc(ig,iv,:,iglo)
             end do
          end do
       end do
       wgt = spec%temp*spec%dens
       call integrate_species (g0, wgt, antotp)
    else
       antotp=0.
    end if

  end subroutine getan

  subroutine getmoms (ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)

    use dist_fn_arrays, only: aj0, aj1, gnew, g_adjust
    use gs2_layouts, only: is_idx, imu_idx, g_lo, ik_idx
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use vpamu_grids, only: vpa, vperp2, integrate_moment, anon, energy, nvgrid
    use run_parameters, only: fphi, fbpar
    use fields_arrays, only: phinew, bparnew
    use kt_grids, only: ntheta0

    implicit none
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: density, &
         upar, tpar, tperp, ntot, qparflux, pperpj1, qpperpj1

    integer :: ik, it, iv, imu, is, iglo, ig

    ! returns moment integrals to PE 0

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
          g0(:,iv,:,iglo) = aj0(:,:,iglo)*gnew(:,iv,:,iglo)
       end do
    end do

! CMR: density is the nonadiabatic piece of perturbed density
! NB normalised wrt equilibrium density for species s: n_s n_ref  
!    ie multiply by (n_s n_ref) to get abs density perturbation
    call integrate_moment (g0, density)

! DJA/CMR: upar and tpar moments 
! (nb adiabatic part of <delta f> does not contribute to upar, tpar or tperp)
! NB UPAR is normalised to vt_s = sqrt(T_s/m_s) vt_ref
!    ie multiply by spec(is)%stm * vt_ref to get abs upar
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid
             g0(ig,:,it,iglo) = vpa*g0(ig,:,it,iglo)
          end do
       end do
    end do
    call integrate_moment (g0, upar)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid
             g0(ig,:,it,iglo) = 2.0*vpa*g0(ig,:,it,iglo)
          end do
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
       do it = 1, ntheta0
          do iv = -nvgrid, nvgrid
             g0(:,iv,it,iglo) = vperp2(:,imu)*gnew(:,iv,it,iglo)*aj0(:,it,iglo)
          end do
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
             g0(ig,iv,:,iglo) = vpa(iv)*gnew(ig,iv,:,iglo)*aj0(ig,:,iglo)*energy(ig,iv,imu)
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
       do it = 1, ntheta0
          do iv = -nvgrid, nvgrid
             g0(:,iv,it,iglo) &
                  = gnew(:,iv,it,iglo)*aj1(:,it,iglo)*2.0*vperp2(:,imu)*spec(is)%tz
          end do
       end do
    end do
    call integrate_moment (g0, pperpj1)

! Now compute QPPERPJ1, a modified p_perp*energy which gives heat flux from Bpar
! NB QPPERPJ1 is normalised to (n_s T_s^2/q_s)  n_ref  T_ref^2/q_ref
!    ie multiply by (n_s T_s spec(is)%tz) n_ref T_ref^2/q_ref to get abs QPPERPJ1
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       imu = imu_idx(g_lo,iglo)
       do it = 1, ntheta0
          do iv = -nvgrid, nvgrid
             g0(:,iv,it,iglo) = g0(:,iv,it,iglo)*energy(:,iv,imu)
          end do
       end do
    end do
    call integrate_moment (g0, qpperpj1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now include the adiabatic part of <delta f>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! set g0 = <delta_f>/F_m, including the adiabatic term
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       imu = imu_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       do it = 1, ntheta0
          do iv = -nvgrid, nvgrid
             g0(:,iv,it,iglo) = aj0(:,it,iglo)*gnew(:,iv,it,iglo) - anon(:,iv,imu)*phinew(:,it,ik)*spec(is)%zt
          end do
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

  end subroutine getmoms

  subroutine init_fieldeq

    use dist_fn_arrays, only: aj0, aj1, kperp2, aj0p
    use species, only: nspec, spec, has_electron_species
    use theta_grid, only: ntgrid, dBdrho, bmag
    use kt_grids, only: naky, ntheta0, aky
    use vpamu_grids, only: anon, integrate_species, vperp2, nvgrid, mu, energy, vpa
    use gs2_layouts, only: g_lo, is_idx, imu_idx, ik_idx
    use run_parameters, only: tite, beta, fapar

    implicit none

    integer :: iglo, iv
    integer :: ik, imu, is, it
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: tot
    real, dimension (nspec) :: wgt

    if (feqinit) return
    feqinit = .true.

    allocate (gamtot(-ntgrid:ntgrid,ntheta0,naky))
    allocate (gamtot1(-ntgrid:ntgrid,ntheta0,naky))
    allocate (gamtot2(-ntgrid:ntgrid,ntheta0,naky))
    if (fapar > epsilon(0.)) then
       allocate (gamtota(-ntgrid:ntgrid,ntheta0,naky)) ; gamtota = 0.
    else
       allocate (gamtota(1,1,1))
    end if
    if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
       allocate (gamtot3(-ntgrid:ntgrid,ntheta0,naky))
    endif
    allocate (gamtotp(-ntgrid:ntgrid,ntheta0,naky)) ! needed for radial profile variation

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       imu = imu_idx(g_lo,iglo)
       do it = 1, ntheta0
          do iv = -nvgrid, nvgrid
             g0(:,iv,it,iglo) = (1.0 - aj0(:,it,iglo)**2)*anon(:,iv,imu)
          end do
       end do
    end do
    wgt = spec%z*spec%z*spec%dens/spec%temp
    call integrate_species (g0, wgt, tot)
    gamtot = real(tot) + kperp2*poisfac

    if (fapar > epsilon(0.)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          imu = imu_idx(g_lo,iglo)
          do it = 1, ntheta0
             do iv = -nvgrid, nvgrid
                g0(:,iv,it,iglo) = (vpa(iv)*aj0(:,it,iglo))**2*anon(:,iv,imu)
             end do
          end do
       end do
       wgt = 2.0*beta*spec%dens*spec%z**2/spec%mass
       call integrate_species (g0, wgt, tot)
       gamtota = real(tot)
    end if
    
    ! for now, ignoring possibility of poisfac =/= 0
    ! also may need to be modified for adiabatic electrons
!    wgt = wgt*(spec%tprim-spec%fprim)
!    call integrate_species (g0, wgt, tot)
!    gamtotp = real(tot)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       imu = imu_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do it = 1, ntheta0
          do iv = -nvgrid, nvgrid
             g0(:,iv,it,iglo) = anon(:,iv,imu)*(2.*aj0(:,it,iglo)*aj0p(:,it,iglo) &
                  + (1.0-aj0(:,it,iglo)**2)*(spec(is)%fprim+spec(is)%tprim*(energy(:,iv,imu)-2.5) &
                  + dBdrho*(2.*mu(imu)*bmag-1./bmag)))
          end do
       end do
    end do
    call integrate_species (g0, wgt, tot)
    where (abs(gamtot) < epsilon(0.))
       gamtotp = 0.
    elsewhere
       gamtotp = real(tot)/gamtot
    end where

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       imu = imu_idx(g_lo,iglo)
       do it = 1, ntheta0
          do iv = -nvgrid, nvgrid
             g0(:,iv,it,iglo) = aj0(:,it,iglo)*aj1(:,it,iglo) &
                  *2.0*vperp2(:,imu)*anon(:,iv,imu)
          end do
       end do
    end do
    wgt = spec%z*spec%dens
    call integrate_species (g0, wgt, tot)
    gamtot1 = real(tot)
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       imu = imu_idx(g_lo,iglo)
       do it = 1, ntheta0
          do iv = -nvgrid, nvgrid
             g0(:,iv,it,iglo) = aj1(:,it,iglo)**2*2.0*vperp2(:,imu)**2*anon(:,iv,imu)
          end do
       end do
    end do
    wgt = spec%temp*spec%dens
    call integrate_species (g0, wgt, tot)
    gamtot2 = real(tot)

! adiabatic electrons 
! assumes that reference density is electron density
! and tite is really Tref/Te
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

  subroutine getfieldeq1 (phifnc, aparfnc, bparfnc, antot, antota, antotp, &
       fieldeq, fieldeqa, fieldeqp)

    use dist_fn_arrays, only: kperp2
    use theta_grid, only: ntgrid, bmag, delthet, jacob
    use kt_grids, only: naky, ntheta0, aky
    use run_parameters, only: fphi, fapar, fbpar
    use run_parameters, only: beta, tite
    use species, only: spec, has_electron_species

    implicit none

    complex, dimension (-ntgrid:,:,:), intent (in) :: phifnc, aparfnc, bparfnc
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
       ! fieldeq is Poisson equation grouped so that RHS = 0; i.e., fieldeq = 0
       fieldeq = antot + bparfnc*gamtot1 - gamtot*phifnc

       if (.not. has_electron_species(spec)) then
          do ik = 1, naky
             do it = 1, ntheta0
                fieldeq(:,it,ik) = fieldeq(:,it,ik) + fl_avg(it,ik)
             end do
          end do
       end if
    end if

    if (fapar > epsilon(0.0)) then
!       fieldeqa = antota - (kperp2+gamtota)*apar
       fieldeqa = antota - kperp2*aparfnc
    end if
    ! bpar == delta B_parallel / B_0(theta) b/c of the factor of 1/bmag(theta)**2
    ! in the following
    if (fbpar > epsilon(0.0)) then
       fieldeqp = (antotp+bparfnc*gamtot2+0.5*phifnc*gamtot1)*beta*apfac
       do ik = 1, naky
          do it = 1, ntheta0
             fieldeqp(:,it,ik) = fieldeqp(:,it,ik)/bmag(:)**2
          end do
       end do
       fieldeqp = fieldeqp + bparfnc
    end if

  end subroutine getfieldeq1

  subroutine getfieldeq (gfnc, phifnc, aparfnc, bparfnc, &
       fieldeq, fieldeqa, fieldeqp)

    use gs2_layouts, only: g_lo
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use vpamu_grids, only: nvgrid

    implicit none

    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in) :: gfnc
    complex, dimension (-ntgrid:,:,:), intent (in) :: phifnc, aparfnc, bparfnc
    complex, dimension (-ntgrid:,:,:), intent (out) ::fieldeq,fieldeqa,fieldeqp

    complex, dimension (:,:,:), allocatable :: antot, antota, antotp

    allocate (antot (-ntgrid:ntgrid,ntheta0,naky))
    allocate (antota(-ntgrid:ntgrid,ntheta0,naky))
    allocate (antotp(-ntgrid:ntgrid,ntheta0,naky))

    ! getan returns velocity space integrals of g needed to 
    ! solve Maxwell equations; e.g., antot = sum_s Z_s int d3v J0_s * g_s / n_ref
    call getan (gfnc, antot, antota, antotp)

!    write (*,*) 'antot', antot(0,1,naky)

    ! getfieldeq1 returns field equations;
    ! e.g., fieldeq = antot + sum_s (Gam0_s - 1)*Z_s^2*e*phi/T_s * n_s / n_ref = 0
    call getfieldeq1 (phifnc, aparfnc, bparfnc, antot, antota, antotp, &
         fieldeq, fieldeqa, fieldeqp)

    deallocate (antot, antota, antotp)

  end subroutine getfieldeq
  
  ! Given initial distribution function this obtains consistent fields
  subroutine get_init_field (gfnc, phi, apar, bpar)
    ! inverts the field equations:
    !   gamtot * phi - gamtot1 * bpar = antot
    !   (kperp2 + gamtota) * apar = antota
    !   beta/2 * gamtot1 * phi + (beta * gamtot2 + 1) * bpar = - beta * antotp
    ! I have not made any check for use_Bpar=T case.
    use gs2_layouts, only: g_lo
    use run_parameters, only: beta, fphi, fapar, fbpar
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use dist_fn_arrays, only: kperp2
    use vpamu_grids, only: nvgrid

    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in) :: gfnc
    complex, dimension (-ntgrid:,:,:), intent (out) :: phi, apar, bpar

    real, dimension (-ntgrid:ntgrid,ntheta0,naky) :: denominator
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: antot, antota, antotp
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: numerator
    real, dimension (-ntgrid:ntgrid,ntheta0,naky) :: bmagsp

    phi=0. ; apar=0. ; bpar=0.
    antot=0.0 ; antota=0.0 ; antotp=0.0
    !CMR, 1/8/2011:  bmagsp is 3D array containing bmag
    bmagsp=spread(spread(bmag,2,ntheta0),3,naky)
    call getan (gfnc, antot, antota, antotp)

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
       denominator = kperp2 + gamtota
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
 
  subroutine flux (gfnc, phi, apar, bpar, &
       pflux,  qflux,  vflux, vflux_par, vflux_perp, &
       pmflux, qmflux, vmflux, &
       pbflux, qbflux, vbflux)

    use species, only: spec
    use theta_grid, only: ntgrid, bmag, gradpar, delthet
    use theta_grid, only: qval, shat, gds21, gds22
    use kt_grids, only: naky, ntheta0, theta0, aky
    use vpamu_grids, only: energy, vpa, vperp2, nvgrid
    use dist_fn_arrays, only: aj0, aj1
    use gs2_layouts, only: g_lo, imu_idx, is_idx, ik_idx
    use mp, only: proc0
    use run_parameters, only: woutunits, fphi, fapar, fbpar
    use constants, only: zi
    use geometry, only: rhoc
    use theta_grid, only: Rplot, Bpol

    implicit none

    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in) :: gfnc
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
          g0(:,iv,:,:) = gfnc(:,iv,:,:)*aj0
       end do
       call get_flux (phi, pflux, dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          imu = imu_idx(g_lo,iglo)
          do it = 1, ntheta0
             g0(:,:,it,iglo) = g0(:,:,it,iglo)*energy(:,:,imu)
          end do
       end do
       call get_flux (phi, qflux(:,:,:,1), dnorm)

       do iv = -nvgrid, nvgrid
          g0(:,iv,:,:) = gfnc(:,iv,:,:)*2.*vpa(iv)**2*aj0
       end do
       call get_flux (phi, qflux(:,:,:,2), dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          imu = imu_idx(g_lo,iglo)
          do it = 1, ntheta0
             do iv = -nvgrid, nvgrid
                g0(:,iv,it,iglo) = gfnc(:,iv,it,iglo)*vperp2(:,imu)*aj0(:,it,iglo)
             end do
          end do
       end do
       call get_flux (phi, qflux(:,:,:,3), dnorm)

       do iv = -nvgrid, nvgrid
          do ig = -ntgrid, ntgrid
             g0(ig,iv,:,:) = gfnc(ig,iv,:,:)*aj0(ig,:,:)*vpa(iv)*Rplot(ig)*sqrt(1.0-Bpol(ig)**2/bmag(ig)**2)
          end do
       end do
       call get_flux (phi, vflux_par, dnorm)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          imu = imu_idx(g_lo,iglo)
          do it = 1, ntheta0
             do iv = -nvgrid, nvgrid
                g0(:,iv,it,iglo) = -zi*aky(ik)*gfnc(:,iv,it,iglo)*aj1(:,it,iglo) &
                     *rhoc*(gds21+theta0(it,ik)*gds22)*vperp2(:,imu)*spec(is)%smz/(qval*shat*bmag**2)
             end do
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
             g0(:,iv,:,iglo) &
                  = -gfnc(:,iv,:,iglo)*aj0(:,:,iglo)*spec(is)%stm*vpa(iv)
          end do
       end do
       call get_flux (apar, pmflux, dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          imu = imu_idx(g_lo,iglo)
          do it = 1, ntheta0
             g0(:,:,it,iglo) = g0(:,:,it,iglo)*energy(:,:,imu)
          end do
       end do
       call get_flux (apar, qmflux(:,:,:,1), dnorm)
       
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do iv = -nvgrid, nvgrid
             g0(:,iv,:,iglo) &
                  = -gfnc(:,iv,:,iglo)*aj0(:,:,iglo)*spec(is)%stm*vpa(iv) &
                  *2.*vpa(iv)**2
          end do
       end do
       call get_flux (apar, qmflux(:,:,:,2), dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          imu = imu_idx(g_lo,iglo)
          do it = 1, ntheta0
             do iv = -nvgrid, nvgrid
                g0(:,iv,it,iglo) &
                     = -gfnc(:,iv,it,iglo)*aj0(:,it,iglo)*spec(is)%stm*vpa(iv) &
                     *vperp2(:,imu)
             end do
          end do
       end do
       call get_flux (apar, qmflux(:,:,:,3), dnorm)
       
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do iv = -nvgrid, nvgrid
             g0(:,iv,:,iglo) &
                  = -gfnc(:,iv,:,iglo)*aj0(:,:,iglo)*spec(is)%stm &
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
          do it = 1, ntheta0
             do iv = -nvgrid, nvgrid
                g0(:,iv,it,iglo) &
                     = gfnc(:,iv,it,iglo)*aj1(:,it,iglo)*2.0*vperp2(:,imu)*spec(is)%tz
             end do
          end do
       end do
       call get_flux (bpar, pbflux, dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          imu = imu_idx(g_lo,iglo)
          do it = 1, ntheta0
             g0(:,:,it,iglo) = g0(:,:,it,iglo)*energy(:,:,imu)
          end do
       end do
       call get_flux (bpar, qbflux(:,:,:,1), dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          imu = imu_idx(g_lo,iglo)
          do it = 1, ntheta0
             do iv = -nvgrid, nvgrid
                g0(:,iv,it,iglo) &
                     = gfnc(:,iv,it,iglo)*aj1(:,it,iglo)*2.0*vperp2(:,imu)*spec(is)%tz &
                     *2.*vpa(iv)**2
             end do
          end do
       end do
       call get_flux (bpar, qbflux(:,:,:,2), dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          imu = imu_idx(g_lo,iglo)
          do it = 1, ntheta0
             do iv = -nvgrid, nvgrid
                g0(:,iv,it,iglo) &
                     = gfnc(:,iv,it,iglo)*aj1(:,it,iglo)*2.0*vperp2(:,imu)*spec(is)%tz &
                     *vperp2(:,imu)
             end do
          end do
       end do
       call get_flux (bpar, qbflux(:,:,:,3), dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          imu = imu_idx(g_lo,iglo)
          do it = 1, ntheta0
             do iv = -nvgrid, nvgrid
                g0(:,iv,it,iglo) &
                     = gfnc(:,iv,it,iglo)*aj1(:,it,iglo)*2.0*vperp2(:,imu) &
                     *spec(is)%tz*vpa(iv)
             end do
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
    use vpamu_grids, only: integrate_moment
    use species, only: nspec
    use mp, only: proc0

    implicit none

    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (:,:,:), intent (in out) :: flx
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
    use gs2_layouts, only: g_lo, imu_idx, ik_idx, is_idx
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: gnew, aj0
    use theta_grid, only: ntgrid, gradpar, delthet, jacob, thet_imp
    use kt_grids, only: ntheta0, naky
    use vpamu_grids, only: integrate_moment, nvgrid, vpac, vpa_imp
    use run_parameters, only: fphi
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
          ik = ik_idx(g_lo,iglo)
          do it = 1, ntheta0
             if (nonlin .and. it==0 .and. ik==1) cycle
             do iv = -nvgrid, nvgrid
                ! get v_magnetic piece of g0 at grid points instead of cell centers
                do ig = -ntgrid, ntgrid
                   g0(ig,iv,it,iglo) = aj0(ig,it,iglo)*(zi*wdrift_func(ig,iv,imu,it,ik)/code_dt &
                        * gnew(ig,iv,it,iglo)*spec(is)%tz)
                end do
             end do
             ! get cell centered value in theta and vpa
             call get_cell_value (thet_imp, vpa_imp, &
                  g0(:,:,it,iglo), g0(:,:,it,iglo), -ntgrid, -nvgrid)
             ! note that gnew below is not cell centered in vpa, which is should be!
             ! MAB FLAG
             do iv = -nvgrid, nvgrid
                ! get v_magnetic piece of g0 at cell centers and add in vpar piece at cell centers
                do ig = -ntgrid, -1
                   g0(ig,iv,it,iglo) = g0(ig,iv,it,iglo) &
                        + vpac(iv,1)*gradparc(ig)/delthet(ig) &
                        * (gnew(ig+1,iv,it,iglo)-gnew(ig,iv,it,iglo))*spec(is)%stm
                end do
                do ig = 0, ntgrid-1
                   g0(ig,iv,it,iglo) = g0(ig,iv,it,iglo) &
                        + vpac(iv,2)*gradparc(ig)/delthet(ig) &
                        * (gnew(ig+1,iv,it,iglo)-gnew(ig,iv,it,iglo))*spec(is)%stm
                end do
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

  subroutine flux_vs_theta_vs_vpa (gfnc, phifnc, pflx, vflx, qflx)

    use constants, only: zi
    use dist_fn_arrays, only: aj1, aj0
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: ik_idx, is_idx, imu_idx
    use geometry, only: rhoc
    use theta_grid, only: ntgrid, bmag, gds21, gds22, qval, shat
    use theta_grid, only: Rplot, Bpol
    use kt_grids, only: aky, theta0, ntheta0
    use vpamu_grids, only: nmu, nvgrid, vpa, vperp2, integrate_mu, energy
    use species, only: spec, nspec

    implicit none

    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in) :: gfnc
    complex, dimension (-ntgrid:,:,:), intent (in) :: phifnc
    real, dimension (-ntgrid:,-nvgrid:,:), intent (out) :: pflx, vflx, qflx
    
    integer :: all = 1
    integer :: iglo, ig, it, ik, is, imu, iv

    real, dimension (:,:), allocatable :: gavg

    allocate (gavg(nmu,nspec))

    ! first compute particle flux contribution
    ! from particles with a given theta and vpa
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       do it = 1, ntheta0
          do iv = -nvgrid, nvgrid
             do ig = -ntgrid, ntgrid
                g0(ig,iv,it,iglo) = gfnc(ig,iv,it,iglo)*aj0(ig,it,iglo)
                g0(ig,iv,it,iglo) = aimag(g0(ig,iv,it,iglo)*conjg(phifnc(ig,it,ik))) &
                     * aky(ik)
             end do
          end do
       end do
    end do

    do iv = -nvgrid, nvgrid
       do ig = -ntgrid, ntgrid
          ! integrate_volume averages over x and y,
          ! leaving gavg a function of mu and species
          call integrate_volume (real(g0(ig,iv,:,:)), gavg, all)
          ! integrate over mu, returning a function of species
          call integrate_mu (ig, gavg, pflx(ig,iv,:))
       end do
    end do

    ! next compute heat flux contribution
    ! from particles as position theta and speed vpa
    ! first compute particle flux contribution
    ! from particles with a given theta and vpa
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       imu = imu_idx(g_lo,iglo)
       do iv = -nvgrid, nvgrid
          do ig = -ntgrid, ntgrid
             g0(ig,iv,:,iglo) = gfnc(ig,iv,:,iglo)*aj0(ig,:,iglo)*energy(ig,iv,imu)
             g0(ig,iv,:,iglo) = aimag(g0(ig,iv,:,iglo)*conjg(phifnc(ig,:,ik))) &
                  * aky(ik)
          end do
       end do
    end do

    do iv = -nvgrid, nvgrid
       do ig = -ntgrid, ntgrid
          ! integrate_volume averages over x and y,
          ! leaving gavg a function of mu and species
          call integrate_volume (real(g0(ig,iv,:,:)), gavg, all)
          ! integrate over mu, returning a function of species
          call integrate_mu (ig, gavg, qflx(ig,iv,:))
       end do
    end do

    ! finally, get contribution to momentum flux
    ! arising from particles with given theta, vpa
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       imu = imu_idx(g_lo,iglo)
       do iv = -nvgrid, nvgrid
          do ig = -ntgrid, ntgrid
             g0(ig,iv,:,iglo) = gfnc(ig,iv,:,iglo)*aj0(ig,:,iglo)*vpa(iv) &
                  *Rplot(ig)*sqrt(1.0-Bpol(ig)**2/bmag(ig)**2) &
                  -zi*aky(ik)*gfnc(ig,iv,:,iglo)*aj1(ig,:,iglo) &
                  *rhoc*(gds21(ig)+theta0(:,ik)*gds22(ig)) &
                  * vperp2(ig,imu)*spec(is)%smz/(qval*shat*bmag(ig)**2)
             g0(ig,iv,:,iglo) = aimag(g0(ig,iv,:,iglo)*conjg(phifnc(ig,:,ik))) &
                  * aky(ik)
          end do
       end do
    end do

    do iv = -nvgrid, nvgrid
       do ig = -ntgrid, ntgrid
          ! integrate_volume averages over x and y,
          ! leaving gavg a function of mu and species
          call integrate_volume (real(g0(ig,iv,:,:)), gavg, all)
          ! integrate over mu, returning a function of species
          call integrate_mu (ig, gavg, vflx(ig,iv,:))
       end do
    end do

    deallocate (gavg)

  end subroutine flux_vs_theta_vs_vpa

  subroutine reset_init

    use dist_fn_arrays, only: gnew, g, source
    use dist_fn_arrays, only: gpnew!, ghnew
    initialized = .false.
    
    wdrift = 0.
    source = 0.
    gnew = 0.
    g0 = 0.
    g = 0.
    gpnew = 0.
!    ghnew = 0.

  end subroutine reset_init

  subroutine reset_physics

    call init_wstar

  end subroutine reset_physics

  subroutine write_f (last)

    use mp, only: proc0, send, receive
    use file_utils, only: open_output_file, close_output_file, get_unused_unit
    use gs2_layouts, only: g_lo, ik_idx, is_idx, imu_idx
    use gs2_layouts, only: idx_local, proc_id
    use vpamu_grids, only: nvgrid, nmu, vpa, vperp2, mu
    use dist_fn_arrays, only: gnew

    integer :: iglo, ik, it, is, iv, ig, imu
    integer, save :: unit
    complex, dimension (-nvgrid:nvgrid) :: gtmp
    logical, save :: first = .true.
    logical, intent(in)  :: last 

    if (first .and. proc0) then
       call open_output_file (unit, ".dist")
!       write(unit, *) (2*nvgrid+1)*nmu
    endif
    first = .false.

    do iglo = g_lo%llim_world, g_lo%ulim_world
       ! writing out g(vpar,vperp) at ik=it=is=1, ig=0
       ik = ik_idx(g_lo, iglo) ; if (ik /= 1) cycle
       is = is_idx(g_lo, iglo) !; if (is /= 1) cycle
       imu = imu_idx(g_lo, iglo) 
       ig = 1 ; it = 1
       if (idx_local (g_lo, ik, imu, is)) then
          if (proc0) then 
             gtmp = gnew(ig,:,it,iglo)
          else
             call send (gnew(ig,:,it,iglo), 0)
          endif
       else if (proc0) then
          call receive (gtmp, proc_id(g_lo, iglo))
       endif
       if (proc0) then
          do iv = -nvgrid, nvgrid
             write (unit, "(5(1x,e12.5),i4)") vpa(iv), sqrt(vperp2(ig,imu)), mu(imu), &
                  real(gtmp(iv)), aimag(gtmp(iv)), is
          end do
          write (unit, *)
       end if
    end do
    if (proc0) write (unit, *)
    if (last .and. proc0) call close_output_file (unit)
    
  end subroutine write_f

  subroutine boundary(linked)

    logical :: linked

    call init_dist_fn
    linked = boundary_option_switch == boundary_option_linked

  end subroutine boundary

  ! subroutine timer (i, place)
    
  !   character (len=10) :: zdate, ztime, zzone
  !   character (*) :: place
  !   integer, intent (in) :: i
  !   integer, dimension(8) :: ival
  !   real, save :: told=0., tnew=0.
    
  !   call date_and_time (zdate, ztime, zzone, ival)
  !   tnew = ival(5)*3600.+ival(6)*60.+ival(7)+ival(8)/1000.
  !   if (i == 0) told = -1.
  !   if (told > 0.) then
  !      print *, ': Elapsed time = ',tnew-told,' seconds in '//trim(place)
  !   end if
  !   told = tnew
  ! end subroutine timer

!   subroutine dot (a, anew, adot, fac)

! ! Get a theta-centered and time-centered estimate of the time derivative 
! ! of a field.
! ! 
! ! tunits(ky) == 1. unless the "omega_*" units are chosen.
! ! omega_* units normalize time by an additional factor of ky.
! !    

!     use run_parameters, only: tunits
!     use gs2_time, only: code_dt
!     use kt_grids, only: naky, ntheta0
!     use theta_grid, only: ntgrid

!     implicit none
!     complex, intent (in), dimension (-ntgrid:,:,:) :: a, anew
!     complex, intent (out), dimension (-ntgrid:,:,:) :: adot
!     real, intent (in) :: fac
!     real :: dtinv
!     integer :: ig, it, ik

!     do ik=1,naky
!        dtinv = 1./(code_dt*tunits(ik))
!        do it=1,ntheta0
!           do ig=-ntgrid,ntgrid-1
!              adot(ig,it,ik) = 0.5*fac*(anew(ig+1,it,ik)+anew(ig,it,ik) - &
!                   (a(ig+1,it,ik)+a(ig,it,ik)))*dtinv
!           end do
!        end do
!     end do
    
!   end subroutine dot

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

!   subroutine get_omega_prime (phipnew_ratio, phipold_ratio, &
!        phihnew_ratio, phihold_ratio)

!     use mp, only: proc0
!     use constants, only: zi
!     use gs2_time, only: code_dt
!     use theta_grid, only: ntgrid, theta
!     use kt_grids, only: naky, ntheta0

!     implicit none

!     complex, dimension (-ntgrid:,:,:), intent (in) :: phipnew_ratio, phipold_ratio
!     complex, dimension (-ntgrid:,:,:), intent (in) :: phihnew_ratio, phihold_ratio

!     integer, save :: idx = 0
!     logical, save :: avg_flag = .false.
!     logical, save :: prim_save_flag = .true., corr_save_flag = .true.
!     integer :: it, ik
!     integer :: navg = 500
!     real :: tol = 0.0005
!     complex, save :: omcorr_save = 0.0, omprim_save = 0.0
!     complex, dimension (-ntgrid:ntgrid) :: dum
!     complex, dimension (ntheta0,naky) :: omega_prime_avg1, omega_prime_avg2
!     complex, dimension (ntheta0,naky) :: omega_correction_avg1, omega_correction_avg2

!     if (.not. allocated(omega_prime)) allocate (omega_prime(ntheta0,naky,navg))
!     if (.not. allocated(omega_correction)) allocate (omega_correction(ntheta0,naky,navg))

!     idx = mod(idx,navg)+1
!     if (idx==navg .and. .not.avg_flag) avg_flag = .true.
!     do ik = 1, naky
!        do it = 1, ntheta0
!           dum = zi*(phipnew_ratio(:,it,ik) - phipold_ratio(:,it,ik))/code_dt
!           call get_fldline_avg (dum, omega_prime(it,ik,idx))

!           dum = zi*(phihnew_ratio(:,it,ik) - phihold_ratio(:,it,ik))/code_dt
!           call get_fldline_avg (dum, omega_correction(it,ik,idx))
!        end do
!     end do

! !    if (avg_flag .and. .not. omcorr_converge) then
!     if (avg_flag) then
!        if (idx < navg/2) then
!           do ik = 1, naky
!              do it = 1, ntheta0
!                 omega_prime_avg1(it,ik) = sum(omega_prime(it,ik,idx+1:idx+navg/2))/(navg/2)
!                 omega_prime_avg2(it,ik) = (sum(omega_prime(it,ik,:idx)) &
!                      + sum(omega_prime(it,ik,navg/2+idx+1:navg)))/(navg/2)
!                 omega_correction_avg1(it,ik) = sum(omega_correction(it,ik,idx+1:idx+navg/2))/(navg/2)
!                 omega_correction_avg2(it,ik) = (sum(omega_correction(it,ik,:idx)) &
!                      + sum(omega_correction(it,ik,navg/2+idx+1:navg)))/(navg/2)
!              end do
!           end do
!        else if (idx > navg/2) then
!           do ik = 1, naky
!              do it = 1, ntheta0
!                 omega_prime_avg1(it,ik) = (sum(omega_prime(it,ik,:idx-navg/2)) &
!                      + sum(omega_prime(it,ik,idx+1:navg)))/(navg/2)
!                 omega_prime_avg2(it,ik) = sum(omega_prime(it,ik,idx-navg/2+1:idx))/(navg/2)
!                 omega_correction_avg1(it,ik) = (sum(omega_correction(it,ik,:idx-navg/2)) &
!                      + sum(omega_correction(it,ik,idx+1:navg)))/(navg/2)
!                 omega_correction_avg2(it,ik) = sum(omega_correction(it,ik,idx-navg/2+1:idx))/(navg/2)

!              end do
!           end do
!        else
!           do ik = 1, naky
!              do it = 1, ntheta0
!                 omega_prime_avg1(it,ik) = sum(omega_prime(it,ik,idx+1:))/(navg/2)
!                 omega_prime_avg2(it,ik) = sum(omega_prime(it,ik,:idx))/(navg/2)
!                 omega_correction_avg1(it,ik) = sum(omega_correction(it,ik,idx+1:))/(navg/2)
!                 omega_correction_avg2(it,ik) = sum(omega_correction(it,ik,:idx))/(navg/2)
!              end do
!           end do
!        end if
       
! !       omprim_converge = all(0.5*abs(omega_prime_avg2-omega_prime_avg1) &
! !            / abs(omega_prime_avg2+omega_prime_avg1) < tol)
! !       omcorr_converge = all(0.5*abs(omega_correction_avg2-omega_correction_avg1) &
! !            / abs(omega_correction_avg2+omega_correction_avg1) < tol)
!        omprim_converge = all(0.5*abs(real(omega_prime_avg2-omega_prime_avg1) &
!             / real(omega_prime_avg2+omega_prime_avg1)) < tol) &
!             .and. all(0.5*abs(aimag(omega_prime_avg2-omega_prime_avg1) &
!             / aimag(omega_prime_avg2+omega_prime_avg1)) < tol)
!        omcorr_converge = all(0.5*abs(real(omega_correction_avg2-omega_correction_avg1) &
!             / real(omega_correction_avg2+omega_correction_avg1)) < tol) &
!             .and. all(0.5*abs(aimag(omega_correction_avg2-omega_correction_avg1) &
!             / aimag(omega_correction_avg2+omega_correction_avg1)) < tol)
!        if (omprim_converge .and. prim_save_flag) then
!           omprim_save = omega_prime_avg2(1,1)
!           prim_save_flag = .false.
!        end if
!        if (omcorr_converge .and. corr_save_flag) then
!           omcorr_save = omega_correction_avg2(1,1)
!           corr_save_flag = .false.
!        end if
!        if (proc0 .and. mod(idx,navg)==0) then
!           do ik = 1, naky
!              do it = 1, ntheta0
!                 write (*,'(a12,4e12.4,l4)') 'omega_prime', real(omega_prime(it,ik,idx)), &
!                      aimag(omega_prime(it,ik,idx)), &
! !                     real(omega_correction(it,ik,idx)), aimag(omega_correction(it,ik,idx)), &
!                      real(omprim_save), aimag(omprim_save), &!real(omcorr_save), aimag(omcorr_save), &
!                      omprim_converge!, omcorr_converge
!              end do
!           end do
!        end if

!     end if

!   end subroutine get_omega_prime

  ! subroutine get_fldline_avg (fld_in, fld_out)

  !   use theta_grid, only: delthet, jacob

  !   implicit none

  !   real, dimension (:), allocatable :: dl_over_b

  !   complex, dimension (-ntg_out:), intent (in) :: fld_in
  !   complex, intent (out) :: fld_out

  !   allocate (dl_over_b(-ntg_out:ntg_out))

  !   dl_over_b = delthet(-ntg_out:ntg_out)*jacob(-ntg_out:ntg_out)
  !   dl_over_b = dl_over_b / sum(dl_over_b)

  !   fld_out = sum(fld_in*dl_over_b)

  !   deallocate (dl_over_b)

  ! end subroutine get_fldline_avg

  ! subroutine takes g(x,y,mu,spec) in g_lo
  ! and averages over x,y, leaving g(mu,spec)
  subroutine integrate_volume (gxy, gavg, all)

    use mp, only: nproc, sum_reduce, sum_allreduce
    use gs2_layouts, only: g_lo, ik_idx, is_idx, imu_idx
    use kt_grids, only: aky, ntheta0

    implicit none

    real, dimension (:,g_lo%llim_proc:), intent (in) :: gxy
    real, dimension (:,:), intent (out) :: gavg
    integer, optional, intent (in) :: all

    integer :: iglo, ik, is, imu, it
    real :: fac

    ! initialize to zero
    gavg = 0.

    ! do integral over local x-y space
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       imu = imu_idx(g_lo,iglo)

       if (aky(ik) < epsilon(0.)) then
          fac = 1.0
       else
          fac = 0.5
       end if

       do it = 1, ntheta0
          gavg(imu,is) = gavg(imu,is) + fac*gxy(it,iglo)
       end do
    end do

    if (nproc > 1) then
       if (present(all)) then
          call sum_allreduce (gavg)
       else
          call sum_reduce (gavg, 0)
       end if
    end if

  end subroutine integrate_volume

  ! subroutine used for testing
  ! takes as input an array using g_lo and
  ! writes it to a .distmp output file
  subroutine write_mpdist_real (dist, extension, last)

    use mp, only: proc0, send, receive
    use file_utils, only: open_output_file, close_output_file
    use gs2_layouts, only: g_lo, ik_idx, is_idx
    use gs2_layouts, only: imu_idx, idx_local, proc_id
    use gs2_time, only: code_time
    use theta_grid, only: ntgrid, bmag, theta
    use vpamu_grids, only: vpa, nvgrid, mu
    use kt_grids, only: theta0, ntheta0

    implicit none
    
    real, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in) :: dist
    character (*), intent (in) :: extension
    logical, intent (in), optional :: last
    
    integer :: iglo, ik, it, is, imu, ig, iv
    integer, save :: unit
    logical, save :: done = .false.
    real :: gtmp
    
    if (.not. done) then
       !        if (proc0) call open_output_file (unit, ".distmp")
       if (proc0) call open_output_file (unit, trim(extension))
       do iglo=g_lo%llim_world, g_lo%ulim_world
          ik = ik_idx(g_lo, iglo)
          is = is_idx(g_lo, iglo) !; if (is /= 1) cycle
          imu = imu_idx(g_lo, iglo)
          do it = 1, ntheta0
             do iv = -nvgrid, nvgrid
                do ig = -ntgrid, ntgrid
                   if (idx_local (g_lo, ik, imu, is)) then
                      if (proc0) then
                         gtmp = dist(ig,iv,it,iglo)
                      else
                         call send (dist(ig,iv,it,iglo), 0)
                      end if
                   else if (proc0) then
                      call receive (gtmp, proc_id(g_lo, iglo))
                   end if
                   if (proc0) then
                      write (unit,'(a1,8e14.4,3i4)') "", code_time, theta(ig), vpa(iv), mu(imu), bmag(ig), &
                           gtmp, theta(ig)-theta0(it,ik), theta0(it,ik), is, ik, it
                   end if
                end do
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
    
  end subroutine write_mpdist_real

  ! subroutine used for testing
  ! takes as input an array using g_lo and
  ! writes it to a .distmp output file
  subroutine write_mpdist_complex (dist, extension, last)

    use mp, only: proc0, send, receive
    use file_utils, only: open_output_file, close_output_file
    use gs2_layouts, only: g_lo, ik_idx, is_idx
    use gs2_layouts, only: imu_idx, idx_local, proc_id
    use gs2_time, only: code_time
    use theta_grid, only: ntgrid, bmag, theta
    use vpamu_grids, only: vpa, nvgrid, mu
    use kt_grids, only: theta0, ntheta0

    implicit none
    
    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in) :: dist
    character (*), intent (in) :: extension
    logical, intent (in), optional :: last
    
    integer :: iglo, ik, it, is, imu, ig, iv
    integer, save :: unit
    logical, save :: done = .false.
    complex :: gtmp
    
    if (.not. done) then
       !        if (proc0) call open_output_file (unit, ".distmp")
       if (proc0) call open_output_file (unit, trim(extension))
       do iglo=g_lo%llim_world, g_lo%ulim_world
          ik = ik_idx(g_lo, iglo)
          is = is_idx(g_lo, iglo) !; if (is /= 1) cycle
          imu = imu_idx(g_lo, iglo)
          do it = 1, ntheta0
             do iv = -nvgrid, nvgrid
                do ig = -ntgrid, ntgrid
                   if (idx_local (g_lo, ik, imu, is)) then
                      if (proc0) then
                         gtmp = dist(ig,iv,it,iglo)
                      else
                         call send (dist(ig,iv,it,iglo), 0)
                      end if
                   else if (proc0) then
                      call receive (gtmp, proc_id(g_lo, iglo))
                   end if
                   if (proc0) then
                      write (unit,'(a1,9e14.4,3i4)') "", code_time, theta(ig), vpa(iv), mu(imu), bmag(ig), &
                           real(gtmp), aimag(gtmp), theta(ig)-theta0(it,ik), theta0(it,ik), is, ik, it
                   end if
                end do
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
    
  end subroutine write_mpdist_complex

  ! subroutine write_response (dist, extension)

  !   use mp, only: proc0, send, receive
  !   use file_utils, only: open_output_file, close_output_file
  !   use gs2_layouts, only: g_lo, ik_idx, is_idx
  !   use gs2_layouts, only: imu_idx, idx_local, proc_id
  !   use gs2_time, only: code_time
  !   use theta_grid, only: ntgrid, bmag, theta, ntheta
  !   use vpamu_grids, only: vpa, nvgrid, mu
  !   use kt_grids, only: theta0, ntheta0, naky

  !   implicit none
    
  !   complex, dimension (:,:,:,g_lo%llim_proc:), intent (in) :: dist
  !   character (*), intent (in) :: extension
    
  !   integer :: iglo, ik, it, is, imu, ig, iv
  !   integer, save :: unit
  !   complex :: gtmp
    
  !   if (proc0) call open_output_file (unit, trim(extension))
  !   do iglo=g_lo%llim_world, g_lo%ulim_world
  !      ik = ik_idx(g_lo, iglo) ; if (ik /= naky) cycle
  !      is = is_idx(g_lo, iglo) !; if (is /= 1) cycle
  !      imu = imu_idx(g_lo, iglo)
  !      do it = 1, neigen(ik)
  !         do iv = 1, (ntheta/2)*nsegments(it,ik)+1
  !            do ig = 1, (ntheta/2)*nsegments(it,ik)+1
  !               if (idx_local (g_lo, ik, imu, is)) then
  !                  if (proc0) then
  !                     gtmp = dist(ig,iv,it,iglo)
  !                  else
  !                     call send (dist(ig,iv,it,iglo), 0)
  !                  end if
  !               else if (proc0) then
  !                  call receive (gtmp, proc_id(g_lo, iglo))
  !               end if
  !               if (proc0) then
  !                  write (unit,'(a1,6i4,2e14.6)') "", ig, iv, it, ik, is, imu, &
  !                       real(gtmp), aimag(gtmp)
  !               end if
  !            end do
  !         end do
  !      end do
  !      if (proc0) then
  !         write (unit,*)
  !         write (unit,*)
  !      end if
  !   end do
  !   if (proc0) call close_output_file (unit)
    
  ! end subroutine write_response
  
  subroutine finish_dist_fn

    use vpamu_grids, only: vperp2, energy, anon, anonc
    use dist_fn_arrays, only: aj0, aj1, aj0p, kperp2, dkperp2dr
    use dist_fn_arrays, only: g, gnew, kx_shift, source, vpar
    use dist_fn_arrays, only: gpnew, vparp, gpold, gold !ghnew, ghold

    implicit none

    no_comm = .false.
    readinit = .false. ; bessinit = .false. ; kp2init = .false. ; connectinit = .false.
    feqinit = .false. ; lpolinit = .false. ; fyxinit = .false. ; cerrinit = .false. ; mominit = .false.
    increase = .true. ; decrease = .false.

    call reset_init

    if (allocated(wdrift)) deallocate (wdrift, wdriftc)
    if (allocated(wdriftp)) deallocate (wdriftp, wdriftpc)
    if (allocated(wdriftmod)) deallocate (wdriftmod, wdriftmodc)
    if (allocated(wstarp)) deallocate (wstarp, wstarpc)
    if (allocated(vperp2)) deallocate (vperp2, energy, anon, anonc, streamfac)
    if (allocated(wstar)) deallocate (wstar, wstarc, varfac, varfacc)
    if (allocated(aj0)) deallocate (aj0, aj1)
    if (allocated(aj0p)) deallocate (aj0p)
    if (allocated(kperp2)) deallocate (kperp2)
    if (allocated(dkperp2dr)) deallocate (dkperp2dr)
    if (allocated(l_links)) deallocate (l_links, r_links)!, n_links)
    if (allocated(M_class)) deallocate (M_class, N_class)
    if (allocated(itleft)) deallocate (itleft, itright)
    if (allocated(connections)) deallocate (connections)
!    if (allocated(g_adj)) deallocate (g_adj)
    if (allocated(g)) deallocate (g, gnew, gold, g0, source, gpnew, gpold)!, ghnew, ghold)
    if (allocated(gexp_1)) deallocate (gexp_1, gexp_2, gexp_3)
!    if (allocated(g_h)) deallocate (g_h, save_h)
    if (allocated(kx_shift)) deallocate (kx_shift)
    if (allocated(jump)) deallocate (jump)
    if (allocated(ikx_indexed)) deallocate (ikx_indexed)
    if (allocated(ufac)) deallocate (ufac)
    if (allocated(gamtot)) deallocate (gamtot, gamtot1, gamtot2, gamtotp)
    if (allocated(gamtota)) deallocate (gamtota)
    if (allocated(gamtot3)) deallocate (gamtot3)
    if (allocated(fl_avg)) deallocate (fl_avg)
    if (allocated(awgt)) deallocate (awgt)

    if (allocated(source0)) deallocate (source0)
    if (allocated(mu0_source)) deallocate (mu0_source)
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

    if (allocated(vpar)) deallocate (vpar)
    if (allocated(mirror)) deallocate (mirror)
    if (allocated(vparp)) deallocate (vparp)
    if (allocated(omega_prime)) deallocate (omega_prime)
    if (allocated(omega_correction)) deallocate (omega_correction)

  end subroutine finish_dist_fn

#ifdef LOWFLOW
  subroutine init_lowflow

    use constants, only: zi
    use centering, only: get_cell_value
!    use dist_fn_arrays, only: vparterm, wdfac, vpac
!    use dist_fn_arrays, only: wstarfac, hneoc, vpar
    use dist_fn_arrays, only: vpar
    use species, only: spec, nspec
    use geometry, only: rhoc
    use theta_grid, only: theta, ntgrid, delthet, gradpar, bmag
    use theta_grid, only: gds23, gds24, gds24_noq, cvdrift_th, gbdrift_th
    use theta_grid, only: drhodpsi, qval, shat, dbdthetc, bmagc, thet_imp
    use kt_grids, only: theta0, ntheta0, naky, aky, akx
    use gs2_time, only: code_dt
    use gs2_layouts, only: g_lo, ik_idx, imu_idx, is_idx
    use run_parameters, only: tunits, wunits, rhostar
    use vpamu_grids, only: nvgrid, mu, vpac, vpa, vperp2, vpa_imp
!    use lowflow, only: get_lowflow_terms
    use file_utils, only: open_output_file, close_output_file
    use mp, only: proc0

    implicit none

!    integer, save :: neo_unit, neophi_unit, neovpth_unit
    integer :: it, ik, imu, is, iv, iglo, ig
!    real, dimension (:,:,:,:,:), allocatable :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
!    real, dimension (:,:,:), allocatable :: vpadhdec, dhdec, dhdxic, hneovpth
!    real, dimension (:), allocatable :: tmp7, tmp8, tmp9
!    real, dimension (:,:,:), allocatable :: cdfac
    real, dimension (:,:), allocatable :: cvdrift_thc, gbdrift_thc
    real, dimension (:,:), allocatable :: vp

    allocate (cvdrift_thc(-ntgrid:ntgrid,2)) ; cvdrift_thc = 0.
    allocate (gbdrift_thc(-ntgrid:ntgrid,2)) ; gbdrift_thc = 0.
    allocate (vp(-ntgrid:ntgrid,-nvgrid:nvgrid)) ; vp = 0.

!    allocate (vpadhdec (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc)) ; vpadhdec = 0.
!    allocate (dhdec (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc)) ; dhdec = 0.
!    allocate (dhdxic (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc)) ; dhdxic = 0.
!    allocate (cdfac (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc)) ; cdfac = 0.

!    if (.not. allocated(vparterm)) then
!       allocate (vparterm(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc)) ; vparterm = 0.
!       allocate (wdfac(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc)) ; wdfac = 0.
!       allocate (hneoc(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc)) ; hneoc = 1.
!       allocate (wstarfac(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc)) ; wstarfac = 0.
!       allocate (wstar_neo(-ntgrid:ntgrid,ntheta0,naky)) ; wstar_neo = 0.
!    end if

    ! allocate (tmp1(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp1 = 0.
    ! allocate (tmp2(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp2 = 0.
    ! allocate (tmp3(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp3 = 0.
    ! allocate (tmp4(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp4 = 0.
    ! allocate (tmp5(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp5 = 0.
    ! allocate (tmp6(-ntgrid:ntgrid,nlambda,negrid,2,nspec)) ; tmp6 = 0.
    ! allocate (tmp7(-ntgrid:ntgrid), tmp8(-ntgrid:ntgrid), tmp9(-ntgrid:ntgrid))
    ! tmp7 = 0. ; tmp8 = 0. ; tmp9 = 0.

!    allocate (hneovpth(-ntgrid:ntgrid,negrid*nlambda,nspec)) ; hneovpth = 0.

    ! tmp1 is dH^{neo}/dE, tmp2 is dH^{neo}/dxi, tmp3 is vpa*dH^{neo}/dE,
    ! tmp4 is dH^{neo}/dr, tmp5 is dH^{neo}/dtheta, tmp6 is H^{neo}
    ! tmp7 is phi^{neo}/dr, tmp8 is dphi^{neo}/dtheta, and tmp9 phi^{neo}
!    call get_lowflow_terms (theta, al, energy, bmag, tmp1, tmp2, tmp3, tmp4, &
!         tmp5, tmp6, tmp7, tmp8, tmp9, lf_default, lf_decompose)
    
    ! if (proc0) then
    !    call open_output_file (neo_unit,".neodist")
    !    write (neo_unit,*) "# all quantities given at theta=0 for species 1"
    !    write (neo_unit,fmt='(10a12)') "# 1) vpa", "2) vpe", "3) energy", "4) vpa/v", &
    !         "5) dH/dE", "6) dH/dxi", "7) vpa*dH/dE", "8) dH/dr", "9) dH/dtheta", "10) H"
    !    do isgn = 1, 2
    !       do il = 1, nlambda
    !          do ie = 1, negrid
    !             if (.not. forbid(ig0,il)) then
    !                write (neo_unit,'(10e12.4)') sign(sqrt(energy(ie)*(1.-al(il)*bmag(ig0))),1.5-real(isgn)), &
    !                     sqrt(energy(ie)*al(il)*bmag(ig0)), energy(ie), &
    !                     sign(sqrt(1.-al(il)*bmag(ig0)),1.5-real(isgn)), &
    !                     tmp1(ig0,il,ie,isgn,1), tmp2(ig0,il,ie,isgn,1), tmp3(ig0,il,ie,isgn,1), &
    !                     tmp4(ig0,il,ie,isgn,1), tmp5(ig0,il,ie,isgn,1), tmp6(ig0,il,ie,isgn,1)
    !             end if
    !          end do
    !          write (neo_unit,*)
    !       end do
    !    end do
    !    call close_output_file (neo_unit)
       
    !    call open_output_file (neovpth_unit,".neothvp")

    !    ! Get Fneo(theta,vpa)
    !    call get_flux_vs_theta_vs_vpa (tmp6, hneovpth)

    !    write (neovpth_unit,'(3a12)') '1) theta', '2) vpa', '3) Fneo'
    !    do ie = 1, negrid*nlambda
    !       do ig = -ntgrid, ntgrid
    !          write (neovpth_unit,'(3e12.4)'), theta(ig), sqrt(energy(negrid))*(1.-2.*(ie-1)/real(negrid*nlambda-1)), &
    !               hneovpth(ig,ie,1)
    !       end do
    !       write (neovpth_unit,*)
    !    end do

    !    call close_output_file (neovpth_unit)

    !    call open_output_file (neophi_unit,".neophi")
    !    write (neophi_unit,*) "# 1) theta, 2) dphi/dr, 3) dphi/dtheta, 4) phi"
    !    do ig = -ntgrid, ntgrid
    !       write (neophi_unit,'(4e14.5)') theta(ig), tmp7(ig), tmp8(ig), tmp9(ig)
    !    end do
    !    call close_output_file (neophi_unit)
    ! end if
    
    ! if set neo_test flag to .true. in input file, GS2 exits after writing out
    ! neoclassical quantities of interest
!    if (neo_test) stop
    
    ! intialize mappings from g_lo to e_lo and lz_lo or to le_lo to facilitate
    ! energy and lambda derivatives in parallel nonlinearity, etc.
!# ifdef USE_LE_LAYOUT
!    call init_map (use_lz_layout=.false., use_e_layout=.false., use_le_layout=.true., test=.false.)
!# else
!    call init_map (use_lz_layout=.true., use_e_layout=.true., use_le_layout=.false., test=.false.)
!# endif

    ! get cell-centered values for some geometric coefficients
    call get_cell_value (1.0-thet_imp, cvdrift_th, cvdrift_thc(:,1), -ntgrid)
    call get_cell_value (thet_imp, cvdrift_th, cvdrift_thc(:,2), -ntgrid)
    call get_cell_value (1.0-thet_imp, gbdrift_th, gbdrift_thc(:,1), -ntgrid)
    call get_cell_value (thet_imp, gbdrift_th, gbdrift_thc(:,2), -ntgrid)
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       imu = imu_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do it = 1, ntheta0
       ! do isgn = 1,2
       !    dhdec(:,isgn,iglo) = tmp1(:,il,ie,isgn,is)
       !    dhdxic(:,isgn,iglo) = tmp2(:,il,ie,isgn,is)
       !    vpadhdec(:,isgn,iglo) = tmp3(:,il,ie,isgn,is)
       !    hneoc(:,isgn,iglo) = 1.0+tmp6(:,il,ie,isgn,is)
       ! end do
       
       ! get cell-centered (in theta) values
       
       ! this is the contribution from dH^{neo}/dtheta (part of v_E dot grad F^{neo})
       ! takes care of part of Eq. 60 in MAB's GS2 notes
       ! if (aky(ik) == 0.0) then
       !    wstarfac(-ntgrid:ntgrid-1,1,iglo) = 0.25*akx(it)/shat &
       !         *(gds24(-ntgrid:ntgrid-1)+gds24(-ntgrid+1:ntgrid)) &
       !         *tmp5(-ntgrid:ntgrid-1,il,ie,1,is)*code_dt
       !    wstarfac(-ntgrid:ntgrid-1,2,iglo) = 0.25*akx(it)/shat &
       !         *(gds24(-ntgrid:ntgrid-1)+gds24(-ntgrid+1:ntgrid)) &
       !         *tmp5(-ntgrid:ntgrid-1,il,ie,2,is)*code_dt
       ! else
       !    wstarfac(-ntgrid:ntgrid-1,1,iglo) = 0.5*wunits(ik) &
       !         *(gds23(-ntgrid:ntgrid-1)+gds23(-ntgrid+1:ntgrid)+theta0(it,ik)*(gds24(-ntgrid:ntgrid-1)+gds24(-ntgrid+1:ntgrid))) &
       !         *tmp5(-ntgrid:ntgrid-1,il,ie,1,is)*code_dt
       !    wstarfac(-ntgrid:ntgrid-1,2,iglo) = 0.5*wunits(ik) &
       !         *(gds23(-ntgrid:ntgrid-1)+gds23(-ntgrid+1:ntgrid)+theta0(it,ik)*(gds24(-ntgrid:ntgrid-1)+gds24(-ntgrid+1:ntgrid))) &
       !         *tmp5(-ntgrid:ntgrid-1,il,ie,2,is)*code_dt
       ! end if

       ! ! this is the contribution from dH^{neo}/dr (part of v_E dot grad F^{neo})
       ! ! takes care of part of Eq. 60 in MAB's GS2 notes
       ! wstarfac(:,:,iglo) = wstarfac(:,:,iglo) + tmp4(:,il,ie,:,is)*code_dt*wunits(ik)
       
       ! this is the contribution from v_E^par . grad F0 (Sec. 5.2 of gs2_notes.pdf)
       ! at grid points
          do iv = -nvgrid, nvgrid
             wstar(:,iv,iglo) = wstar(:,iv,iglo) &
                  - zi*rhostar*gds24_noq*drhodpsi*rhoc/qval*code_dt &
                  * (spec(is)%fprim+spec(is)%tprim*(vpa(iv)**2+vperp2(:,imu)-1.5))
          end do
          ! get (theta-vpa) cell values for wstar from grid values
          call get_cell_value (thet_imp, vpa_imp, &
               wstar(:,:,iglo), wstarc(:,:,iglo), -ntgrid, -nvgrid)


       ! if (.not. lf_default) then
       !    ! this is the contribution from the last term of the 2nd line of Eq. 43 in 
       !    ! MAB's GS2 notes (arises because the NEO dist. fn. is given at v/vt, and
       !    ! the vt varies radially, so v_E . grad F1 must take this into account)
       !    ! If lf_default is true this is taken care of by multiplying the vt normalization
       !    ! of Chebyshev polynomial argument by appropriate temperature factor
       !    ! when constructing neoclassical distribution function (see lowflow.f90).
       !    ! Otherwise, the below line of code takes care of it.
       !    wstarfac(:,:,iglo) = wstarfac(:,:,iglo) + code_dt*wunits(ik) &
       !         *spec(is)%tprim*energy(ie)*dhdec(:,:,iglo)
       ! end if
       
       ! wstarfac(ntgrid,:,iglo) = 0.0
       
       ! wdfac takes care of the gbdrift part of Eq. 60 of MAB's GS2 notes, as well
       ! as part of the curvature drift term of Eq. 54.  the other part of 
       ! the curvature drift term of Eq. 54 is dealt with by cdfac below.
       ! no code_dt in wdfac because it multiplies wdrift, which has code_dt in it
       ! wdfac(-ntgrid:ntgrid-1,1,iglo) = 0.5*dhdxic(-ntgrid:ntgrid-1,1,iglo)*vpac(-ntgrid:ntgrid-1,1,iglo)/energy(ie)**1.5 &
       !      - dhdec(-ntgrid:ntgrid-1,1,iglo)
       ! wdfac(-ntgrid:ntgrid-1,2,iglo) = 0.5*dhdxic(-ntgrid:ntgrid-1,2,iglo)*vpac(-ntgrid:ntgrid-1,2,iglo)/energy(ie)**1.5 &
       !      - dhdec(-ntgrid:ntgrid-1,2,iglo)
       ! wdfac(ntgrid,:,iglo) = 0.0

       ! takes care of part of curvature drift term in Eq. 54 of MAB's GS2 notes.
       ! no code_dt in cdfac because it multiples wcurv, which has code_dt in it.
       ! cdfac(-ntgrid:ntgrid-1,1,iglo) = -0.5*dhdxic(-ntgrid:ntgrid-1,1,iglo)*vpac(-ntgrid:ntgrid-1,1,iglo)/sqrt(energy(ie))
       ! cdfac(-ntgrid:ntgrid-1,2,iglo) = -0.5*dhdxic(-ntgrid:ntgrid-1,2,iglo)*vpac(-ntgrid:ntgrid-1,2,iglo)/sqrt(energy(ie))
       ! cdfac(ntgrid,:,iglo) = 0.0
       
       ! this is the first term multiplying dF_1/dE in Eq. 42 of MAB's GS2 notes
       ! i.e. Z_s * e * vpa . grad phi^{tb} d(F1/F0)/dE
       ! note there will be a correction to this added below
       ! because actual term appearing in GKE ~ (1/F0) * d(F1)/dE
       ! also note that vpadhdec is vpa*d(F1/F0)/dE at fixed mu (not xi)
       ! vparterm(-ntgrid:ntgrid-1,1,iglo) = spec(is)%zstm*tunits(ik)*code_dt &
       !      /delthet(-ntgrid:ntgrid-1) &
       !      * (abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid))) &
       !      * vpadhdec(-ntgrid:ntgrid-1,1,iglo)
       ! vparterm(-ntgrid:ntgrid-1,2,iglo) = spec(is)%zstm*tunits(ik)*code_dt &
       !      /delthet(-ntgrid:ntgrid-1) &
       !      * (abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid))) &
       !      * vpadhdec(-ntgrid:ntgrid-1,2,iglo)
       ! vparterm(ntgrid,:,iglo) = 0.0
       
          ! get the proper upwinded cell value for the parallel velocity
          do iv = -nvgrid, -1
             where (dbdthetc(:ntgrid-1,1) < 0.0)
                vp(:ntgrid-1,iv) = vpac(iv,1)
             elsewhere
                vp(:ntgrid-1,iv) = vpac(iv,2)
             end where
          end do
          do iv = 0, nvgrid-1
             where (dbdthetc(:ntgrid-1,2) < 0.0)
                vp(:ntgrid-1,iv) = vpac(iv,1)
             elsewhere
                vp(:ntgrid-1,iv) = vpac(iv,2)
             end where
          end do

          ! redefine vpar from vpa bhat dot grad theta to
          ! vpa bhat dot grad theta + v_Magnetic dot grad theta
          ! this accounts for the terms in Sec. 6.1 of MAB's GS3 notes
          do iv = -nvgrid, -1
             vpar(-ntgrid:ntgrid-1,iv,is) = vpar(-ntgrid:ntgrid-1,iv,is) + &
                  rhostar*tunits(ik)*code_dt/delthet(-ntgrid:ntgrid-1) &
                  *(cvdrift_thc(-ntgrid:ntgrid-1,1)*vp(-ntgrid:ntgrid-1,iv)**2 &
                  + gbdrift_thc(-ntgrid:ntgrid-1,1)*0.5*mu(imu)*bmagc(-ntgrid:ntgrid-1,1))
          end do
          do iv = 0, nvgrid-1
             vpar(-ntgrid:ntgrid-1,iv,is) = vpar(-ntgrid:ntgrid-1,iv,is) + &
                  rhostar*tunits(ik)*code_dt/delthet(-ntgrid:ntgrid-1) &
                  *(cvdrift_thc(-ntgrid:ntgrid-1,2)*vp(-ntgrid:ntgrid-1,iv)**2 &
                  + gbdrift_thc(-ntgrid:ntgrid-1,2)*0.5*mu(imu)*bmagc(-ntgrid:ntgrid-1,2))
          end do

       ! if (aky(ik) == 0) then
       !    wstar_neo(-ntgrid:ntgrid-1,it,ik) = -0.25*code_dt*akx(it) &
       !         * tmp8(-ntgrid:ntgrid-1)*(gds24_noq(-ntgrid:ntgrid-1)+gds24_noq(-ntgrid+1:ntgrid))
       ! else
       !    wstar_neo(-ntgrid:ntgrid-1,it,ik) = -code_dt*wunits(ik)*(tmp7(-ntgrid:ntgrid-1) &
       !         + 0.5*tmp8(-ntgrid:ntgrid-1)*(gds23(-ntgrid:ntgrid-1)+gds23(-ntgrid+1:ntgrid) &
       !         + theta0(it,ik)*shat*(gds24_noq(-ntgrid:ntgrid-1)+gds24_noq(-ntgrid+1:ntgrid))))
       ! end if
       ! wstar_neo(ntgrid,it,ik) = 0.0
       end do
    end do

    ! TMP FOR TESTING -- MAB
!    wdfac = 0. ; cdfac = 0. ; hneoc = 1. ; wstarfac = 0. ; vparterm = 0.
 
    ! vparterm is -2*vpar*(1+H^{neo}) - Ze*(vpa . grad phi^{tb})*(dH^{neo}/dE)
    ! note that vpar has contribution from v_{magnetic} . grad theta in it
    ! hneoc = 1 + H^{neo} below accounts for usual parallel streaming source term,
    ! as well as first of three terms multiplying F_1 in Eq. 42 of MAB's GS2 notes
!    vparterm = -2.0*vpar*hneoc + vparterm
    
    ! hneoc below accounts for usual v_magnetic source term, 
    ! as well as second of three terms multiplying F_1 in Eq. 42 of MAB's GS2 notes
    ! wdfac on RHS deals with grad-B drift part of Eq. 60, as well as part of 
    ! curvature drift term in Eq. 54
!    wdfac = wdfac + hneoc
    
    ! do iglo = g_lo%llim_proc, g_lo%ulim_proc
    !    ! need to put in do it = 1, ntheta0 loop here !!!
    !    ik = ik_idx(g_lo,iglo)
    !    ie = ie_idx(g_lo,iglo)
    !    is = is_idx(g_lo,iglo)
    !    wdfac(:,1,iglo) = wdfac(:,1,iglo)*wdriftc(:,1,iglo) &
    !         + cdfac(:,1,iglo)*wcurv(:,iglo)
    !    wdfac(:,2,iglo) = wdfac(:,2,iglo)*wdriftc(:,2,iglo) &
    !         + cdfac(:,2,iglo)*wcurv(:,iglo)

    !    ! add Z^2 * sqrt(m*T^3) * vpa * (bhat . grad(theta)) * d(phineo)/dth term to wdrift
    !    ! this modification of wdrift takes care of -g term in Eq. 67 of MB's gs2_notes
    !    ! arises from pulling 1/F0 inside dg/dE term
    !     wdriftc(-ntgrid:ntgrid-1,1,iglo) = wdriftc(-ntgrid:ntgrid-1,1,iglo) &
    !          - zi*code_dt*tunits(ik)*vpac(-ntgrid:ntgrid-1,1,iglo) &
    !          * spec(is)%zt*spec(is)%zstm*tmp8(-ntgrid:ntgrid-1) &
    !          * 0.5*(abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid)))
    !     wdriftc(-ntgrid:ntgrid-1,2,iglo) = wdriftc(-ntgrid:ntgrid-1,2,iglo) &
    !          - zi*code_dt*tunits(ik)*vpac(-ntgrid:ntgrid-1,2,iglo) &
    !          * spec(is)%zt*spec(is)%zstm*tmp8(-ntgrid:ntgrid-1) &
    !          * 0.5*(abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid)))
       
    !    ! hneoc below accounts for usual wstar term, as well as last of three terms
    !    ! multiplying F_1 in Eq. 42 of MAB'S GS2 notes
    !    ! note that hneo is necessary here because NEO's dist fn is
    !    ! normalized by F0(r) instead of F0 at center of simulation domain as in GS2
    !    wstarfac(:,:,iglo) = wstar(:,:,iglo)*hneoc(:,:,iglo) - wstarfac(:,:,iglo)
    ! end do

    deallocate (cvdrift_thc, gbdrift_thc, vp)
!    deallocate (tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9)
!    deallocate (vpadhdec,dhdec,dhdxic,cdfac,hneovpth)

  end subroutine init_lowflow
#endif

  ! subroutine write_connectinfo

  !   use mp, only: proc0, send, receive, nproc, iproc, mp_gather, barrier
  !   use file_utils, only: open_output_file, close_output_file, get_unused_unit
  !   use gs2_layouts, only: g_lo, ik_idx, is_idx, imu_idx
  !   use gs2_layouts, only: idx_local, proc_id
  !   use vpamu_grids, only: nvgrid, nmu, vpa, vperp2, mu
  !   use gs2_time, only: user_time
  !   use dist_fn_arrays, only: g, gnew

  !   integer :: iglo, ik, it, ip, iseg, is, imu
  !   integer, save :: unit
  !   integer, dimension (0:nproc-1) :: ulim, llim, ulim_alloc

  !   integer, dimension (:), allocatable :: igmod, itmod

  !   if (proc0) call open_output_file (unit, ".connectinfo")

  !   allocate (itmod(nseg_max)) ; itmod = -1
  !   allocate (igmod(nseg_max)) ; igmod = -1

  !   do iglo = g_lo%llim_proc, g_lo%ulim_proc

  !      ik = ik_idx(g_lo,iglo)

  !      do it = 1, neigen(ik)

  !         ! remap to start at theta0 = -theta0_max for this set of connected theta0s
  !         iseg = 1
  !         itmod(iseg) = it + it_shift_left(it)
  !         if (nsegments(it,ik) > 1) then
  !            do iseg = 2, nsegments(it,ik)
  !               itmod(iseg) = itmod(iseg-1) + iglo_shift(itmod(iseg-1),ik)
  !            end do
  !         end if
  !      end do
  !   end do

  !   call mp_gather (g_lo%llim_proc, llim)
  !   call mp_gather (g_lo%ulim_proc, ulim)
  !   call mp_gather (g_lo%ulim_alloc, ulim_alloc)

  !   do iglo = g_lo%llim_world, g_lo%ulim_world

  !      ip = proc_id(g_lo,iglo)

  !      ik = ik_idx(g_lo, iglo)
  !      is = is_idx(g_lo, iglo)
  !      imu = imu_idx(g_lo, iglo)

  !      if (idx_local (g_lo, ik, imu, is)) then
  !         if (proc0) then 
  !            igmod = iglomod(:,iglo)
  !         else
  !            call send (iglomod(:,iglo), 0)
  !         endif
  !      else if (proc0) then
  !         call receive (igmod, ip)
  !      endif

  !      if (proc0) then
  !         do iseg = 1, nsegments(it,ik)
  !            write (unit,'(9i6)') ip, llim(ip), ulim(ip), ulim_alloc(ip), iglo, iseg, it, ik, igmod(iseg)
  !         end do
  !      end if

  !   end do



  !   ! do ip = 0, nproc-1
  !   !    if (ip==iproc) then
  !   !       if (proc0) then
  !   !          llim(ip) = g_lo%llim_proc
  !   !          ulim(ip) = g_lo%ulim_proc
  !   !       else
  !   !          call send (g_lo%llim_proc, 0)
  !   !          call send (g_lo%ulim_proc, 0)
  !   !       endif
  !   !    else if (proc0) then
  !   !       call receive (llim(ip), ip)
  !   !       call receive (ulim(ip), ip)
  !   !    endif
  !   ! end do

  ! end subroutine write_connectinfo

end module dist_fn
