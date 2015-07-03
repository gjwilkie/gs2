module nonlinear_terms
!
! not correct for slowing down species
! terms like df/dU should involve anon, etc.
!
! missing factors of B_a/B(theta) in A_perp terms??
!
  implicit none

  private

  public :: init_nonlinear_terms, finish_nonlinear_terms
  public :: read_parameters, wnml_nonlinear_terms, check_nonlinear_terms
  public :: add_explicit_terms, nonlinear_mode_switch
  public :: finish_init, reset_init, algorithm, nonlin, accelerated
  public :: nonlinear_terms_unit_test_time_add_nl, cfl
  public :: nonlinear_mode_none, nonlinear_mode_on
  public :: gryfx_zonal

  type gs2_gryfx_zonal_type

    logical :: on = .false.
    complex*8, dimension (:), pointer :: NLdens_ky0, NLupar_ky0, NLtpar_ky0, &
                                         NLtprp_ky0, NLqpar_ky0, NLqprp_ky0
    logical :: first_half_step = .false.

  end type gs2_gryfx_zonal_type

  type(gs2_gryfx_zonal_type) gryfx_zonal

  integer :: istep_last = 0

  ! knobs
  integer :: nonlinear_mode_switch
  integer :: flow_mode_switch !THIS IS NOT SUPPORTED, SHOULD BE REMOVED?

  integer, parameter :: nonlinear_mode_none = 1, nonlinear_mode_on = 2
  integer, parameter :: flow_mode_off = 1, flow_mode_on = 2

  real, dimension (:,:), allocatable :: ba, gb, bracket
  ! yxf_lo%ny, yxf_lo%llim_proc:yxf_lo%ulim_alloc

  real, dimension (:,:,:), allocatable :: aba, agb, abracket
  ! 2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc

! CFL coefficients
  real :: cfl, cflx, cfly

! hyperviscosity coefficients
  real :: C_par, C_perp, p_x, p_y, p_z

  integer :: algorithm = 1
  logical :: nonlin = .false.
  logical :: initialized = .false.
  logical :: initializing = .true.
  logical :: alloc = .true.
  logical :: zip = .false.
  logical :: nl_forbid_force_zero = .true.
  logical :: accelerated = .false.

  logical :: exist
  
contains
  
  subroutine check_nonlinear_terms(report_unit,delt_adj)
    use gs2_time, only: code_dt_min
    use kt_grids, only: box
    use run_parameters, only: margin, code_delt_max, nstep, wstar_units
    use theta_grid, only: nperiod
    implicit none
    integer, intent(in) :: report_unit
    real, intent(in) :: delt_adj
    if (nonlin) then
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
       if ( .not. box) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Nonlinear runs must be carried out in a box.')") 
          write (report_unit, fmt="('Set grid_option to box in the kt_grids_knobs namelist.')") 
          write (report_unit, fmt="('THIS IS AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if

       if (nperiod > 1) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Nonlinear runs usually have nperiod = 1.')") 
          write (report_unit, fmt="('THIS MAY BE AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
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
       write (report_unit, fmt="('The minimum delt ( code_dt_min ) = ',e11.4)") code_dt_min
       write (report_unit, fmt="('The maximum delt (code_delt_max) = ',e11.4)") code_delt_max
       write (report_unit, fmt="('The maximum delt < ',f10.4,' * min(Delta_perp/v_perp). (cfl)')") cfl
       write (report_unit, fmt="('When the time step needs to be changed, it is adjusted by a factor of ',f10.4)") delt_adj
       write (report_unit, fmt="('The number of time steps nstep = ',i7)") nstep
       write (report_unit, fmt="('If running in batch mode on the NERSC T3E, the run will stop when ', &
            & f7.4,' % of the time remains.')") 100.*margin
    endif
  end subroutine check_nonlinear_terms

  subroutine wnml_nonlinear_terms(unit)
    implicit none
    integer, intent(in) :: unit
    if (.not. exist) return
    if (nonlinear_mode_switch == nonlinear_mode_on) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "nonlinear_terms_knobs"
       write (unit, fmt="(' nonlinear_mode = ',a)") '"on"'
       write (unit, fmt="(' cfl = ',e17.10)") cfl
       write (unit, fmt="(' nl_forbid_force_zero = ',L1)") nl_forbid_force_zero
       if (zip) write (unit, fmt="(' zip = ',L1)") zip
       write (unit, fmt="(' /')")
    endif
  end subroutine wnml_nonlinear_terms

  subroutine init_nonlinear_terms 
    use theta_grid, only: init_theta_grid, ntgrid
    use kt_grids, only: init_kt_grids, naky, ntheta0, nx, ny, akx, aky
    use le_grids, only: init_le_grids, nlambda, negrid
    use species, only: init_species, nspec
    use gs2_layouts, only: init_dist_fn_layouts, yxf_lo, accelx_lo
    use gs2_layouts, only: init_gs2_layouts
    use gs2_transforms, only: init_transforms
    implicit none
    logical :: dum1, dum2
    logical, parameter :: debug=.false.

    if (initialized) return
    initialized = .true.
    
    if (debug) write(6,*) "init_nonlinear_terms: init_gs2_layouts"
    call init_gs2_layouts
    if (debug) write(6,*) "init_nonlinear_terms: init_theta_grid"
    call init_theta_grid
    if (debug) write(6,*) "init_nonlinear_terms: init_kt_grids"
    call init_kt_grids
    if (debug) write(6,*) "init_nonlinear_terms: init_le_grids"
    call init_le_grids (dum1, dum2)
    if (debug) write(6,*) "init_nonlinear_terms: init_species"
    call init_species
    if (debug) write(6,*) "init_nonlinear_terms: init_dist_fn_layouts"
    call init_dist_fn_layouts (naky, ntheta0, nlambda, negrid, nspec)
    
    call read_parameters

    if (debug) write(6,*) "init_nonlinear_terms: init_transforms"
    if (nonlinear_mode_switch /= nonlinear_mode_none) then
       call init_transforms (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)

       if (debug) write(6,*) "init_nonlinear_terms: allocations"
       if (alloc) then
          if (accelerated) then
             allocate (     aba(2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc))
             allocate (     agb(2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc))
             allocate (abracket(2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc))
             aba = 0. ; agb = 0. ; abracket = 0.
          else
             allocate (     ba(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc))
             allocate (     gb(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc))
             allocate (bracket(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc))
             ba = 0. ; gb = 0. ; bracket = 0.
          end if
          alloc = .false.
       end if

       cfly = aky(naky)/cfl*0.5
       cflx = akx((ntheta0+1)/2)/cfl*0.5
    end if

  end subroutine init_nonlinear_terms

  subroutine read_parameters
    use file_utils, only: input_unit_exist, error_unit
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    implicit none
    type (text_option), dimension (4), parameter :: nonlinearopts = &
         (/ text_option('default', nonlinear_mode_none), &
            text_option('none', nonlinear_mode_none), &
            text_option('off', nonlinear_mode_none), &
            text_option('on', nonlinear_mode_on) /)
    character(20) :: nonlinear_mode
    type (text_option), dimension (3), parameter :: flowopts = &
         (/ text_option('default', flow_mode_off), &
            text_option('off', flow_mode_off), &
            text_option('on', flow_mode_on) /)
    character(20) :: flow_mode
    namelist /nonlinear_terms_knobs/ nonlinear_mode, flow_mode, cfl, &
         C_par, C_perp, p_x, p_y, p_z, zip, nl_forbid_force_zero
    integer :: ierr, in_file

    if (proc0) then
       nonlinear_mode = 'default'
       flow_mode = 'default'
       cfl = 0.1
       C_par = 0.1
       C_perp = 0.1
       p_x = 6.0
       p_y = 6.0
       p_z = 6.0

       in_file=input_unit_exist("nonlinear_terms_knobs",exist)
       if(exist) read (unit=in_file,nml=nonlinear_terms_knobs)

       ierr = error_unit()
       call get_option_value &
            (nonlinear_mode, nonlinearopts, nonlinear_mode_switch, &
            ierr, "nonlinear_mode in nonlinear_terms_knobs",.true.)
       call get_option_value &
            (flow_mode, flowopts, flow_mode_switch, &
            ierr, "flow_mode in nonlinear_terms_knobs",.true.)
    end if

    call broadcast (nonlinear_mode_switch)
    call broadcast (flow_mode_switch)
    call broadcast (cfl)
    call broadcast (C_par) 
    call broadcast (C_perp) 
    call broadcast (p_x)
    call broadcast (p_y)
    call broadcast (p_z)
    call broadcast (nl_forbid_force_zero)
    call broadcast (zip)

    istep_last = 0

    if (flow_mode_switch == flow_mode_on) then
       if (proc0) write(*,*) 'Forcing flow_mode = off'
       flow_mode_switch = flow_mode_off
    endif

    if (nonlinear_mode_switch == nonlinear_mode_on)  then
       algorithm = 1 
       nonlin = .true.
    end if

  end subroutine read_parameters

  subroutine add_explicit_terms (g1, g2, g3, phi, apar, bpar, istep, bd)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use gs2_time, only: save_dt_cfl
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1, g2, g3
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    integer, intent (in) :: istep
    real, intent (in) :: bd
    real :: dt_cfl
    logical, parameter :: nl = .true.

    select case (nonlinear_mode_switch)
    case (nonlinear_mode_none)
!!! NEED TO DO SOMETHING HERE...  BD GGH
       dt_cfl = 1.e8
       call save_dt_cfl (dt_cfl)
#ifdef LOWFLOW
       if (istep /=0) &
            call add_explicit (g1, g2, g3, phi, apar, bpar, istep, bd)
#endif
    case (nonlinear_mode_on)
!       if (istep /= 0) call add_nl (g1, g2, g3, phi, apar, bpar, istep, bd, fexp)
       if (istep /= 0) then
         call add_explicit (g1, g2, g3, phi, apar, bpar, istep, bd, nl)
       endif
    end select
  end subroutine add_explicit_terms

  subroutine add_explicit (g1, g2, g3, phi, apar, bpar, istep, bd,  nl)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use dist_fn_arrays, only: g
    use gs2_time, only: save_dt_cfl
    use run_parameters, only: reset
    use unit_tests, only: debug_message
    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1, g2, g3
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    integer, intent (in) :: istep
    real, intent (in) :: bd
    logical, intent (in), optional :: nl

    integer, parameter :: verb = 4
    integer :: iglo, ik, it
    real :: zero
    real :: dt_cfl

    call debug_message(verb, 'nonlinear_terms::add_explicit starting')

    if (initializing) then
       if (present(nl)) then
          dt_cfl = 1.e8
          call save_dt_cfl (dt_cfl)
       end if
       return
    endif


!
! Currently not self-starting.  Need to fix this.
!

    call debug_message(verb, 'nonlinear_terms::add_explicit copying old terms')

    if (istep /= istep_last) then

       zero = epsilon(0.0)
       g3 = g2
       g2 = g1

       call debug_message(verb, 'nonlinear_terms::add_explicit copied old terms')

       ! if running nonlinearly, then compute the nonlinear term at grid points
       ! and store it in g1
       if (present(nl)) then
          call debug_message(verb, 'nonlinear_terms::add_explicit calling add_nl')
          if(gryfx_zonal%on) then
            call add_nl_gryfx (g1) 
          else 
            call add_nl (g1, phi, apar, bpar)
          endif
          if(reset) return !Return if resetting
          ! takes g1 at grid points and returns 2*g1 at cell centers
          call center (g1)
       else
          g1 = 0.
       end if

#ifdef LOWFLOW
       ! do something
#endif

    end if

    if (zip) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
!          if (it == 3 .or. it == ntheta0-1) then
!          if (it /= 1) then
          if (ik /= 1) then
!          if (ik == 2 .and. it == 1) then
             g (:,1,iglo) = 0.
             g (:,2,iglo) = 0.
             g1(:,1,iglo) = 0.
             g1(:,2,iglo) = 0.
          end if
       end do
    end if
     
    istep_last = istep

  contains

    ! subroutine 'center' takes input array evaluated at theta grid points
    ! and overwrites it with array evaluated at cell centers
    ! note that there is an extra factor of 2 in output array
    subroutine center (gtmp)
!
!CMR, 13/10/2014
! Fixing some (cosmetic) issues from prior to R3021.
! (1) forbid(-ntgrid:ntgrid) refers to grid points at cell boundaries
!     gtmp(-ntgrid:ntgrid-1) refers to cell centers
!     => applying forbid to output gtmp did not zero the correct things!!!
! (2) totally trapped particles are special, so return source without upwinding is fine
!     BUT source for other trapped particles should NOT enter forbidden region.
!     if ig or ig+1 forbidden (i.e. ig+1 or ig is bounce point) now return gtmp(ig,:,iglo)=0
!
      use dist_fn_arrays, only: ittp
      use gs2_layouts, only: g_lo, il_idx
      use theta_grid, only: ntgrid
      use le_grids, only: forbid, ng2, jend
      use mp, only: mp_abort

      implicit none

      integer :: iglo, il, ig
      complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: gtmp
      logical :: trapped = .false.

      if (minval(jend) .gt. ng2) trapped=.true.

      ! factor of one-half appears elsewhere
      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         il = il_idx(g_lo, iglo)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CMR, 7/10/2014: 
! Incoming gtmp SHOULD vanish in forbidden region, 
! New logical in nonlinear_terms_knobs namelist: nl_forbid_force_zero
!     nl_forbid_force_zero =.t. : force zeroing    (default)
!     nl_forbid_force_zero =.f. : NO forced zeroing
!                                 ie assume forbidden gtmp is zero on entry 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if ( nl_forbid_force_zero ) then
         ! force spurious gtmp outside trapped boundary to be zero
            where (forbid(:,il))
               gtmp(:,1,iglo) = 0.0
               gtmp(:,2,iglo) = 0.0
            end where
         endif

         do ig = -ntgrid, ntgrid-1
!
!CMR, 7/10/2014: 
! loop sets gtmp to value at upwinded cell center RIGHT of theta(ig)
!           except for ttp where upwinding makes no sense!
!
            if (il >= ittp(ig)) cycle
            if ( trapped .and. ( il > jend(ig) .or. il > jend(ig+1)) ) then
               !
!CMR, 7/10/2014: 
!   if either ig or ig+1 is forbidden, no source possible in a cell RIGHT of theta(ig) 
!   => gtmp(ig,1:2,iglo)=0
!
               gtmp(ig,1:2,iglo) = 0.0
            else 
!
!CMR, 7/10/2014: 
!    otherwise ig and ig+1 BOTH allowed, and upwinding in cell RIGHT of theta(ig) is fine
!
               gtmp(ig,1,iglo) = (1.+bd)*gtmp(ig+1,1,iglo) + (1.-bd)*gtmp(ig,1,iglo)
               gtmp(ig,2,iglo) = (1.-bd)*gtmp(ig+1,2,iglo) + (1.+bd)*gtmp(ig,2,iglo)

            endif
         end do
      end do
    end subroutine center
  end subroutine add_explicit

  subroutine nonlinear_terms_unit_test_time_add_nl(g1, phi, apar, bpar)
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use theta_grid, only: ntgrid
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    call add_nl(g1, phi, apar, bpar)
  end subroutine nonlinear_terms_unit_test_time_add_nl

  subroutine add_nl (g1, phi, apar, bpar)
    use mp, only: max_allreduce
    use theta_grid, only: ntgrid, kxfac
    use gs2_layouts, only: g_lo, ik_idx, it_idx
    use gs2_layouts, only: accelx_lo, yxf_lo
    use dist_fn_arrays, only: g, g_adjust
    use species, only: spec
    use gs2_transforms, only: transform2, inverse2
    use run_parameters, only: fapar, fbpar, fphi, reset, immediate_reset
    use kt_grids, only: aky, akx
    use gs2_time, only: save_dt_cfl, check_time_step_too_large
    use constants, only: zi
    use unit_tests, only: debug_message
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    integer :: i, j, k
    real :: max_vel, zero
    real :: dt_cfl
    integer, parameter :: verb = 4

    integer :: iglo, ik, it, is, ig, ia
    
    call debug_message(verb, 'nonlinear_terms::add_nl starting')

    !Initialise zero so we can be sure tests are sensible
    zero = epsilon(0.0)

    !Form g1=i*kx*chi
    if (fphi > zero) then
       call load_kx_phi
    else
       g1 = 0.
    end if

    if (fbpar > zero) call load_kx_bpar
    if (fapar  > zero) call load_kx_apar

    call debug_message(verb, 'nonlinear_terms::add_nl calling transform2')
    !Transform to real space
    if (accelerated) then
       call transform2 (g1, aba, ia)
    else
       call transform2 (g1, ba)
    end if

    call debug_message(verb, 'nonlinear_terms::add_nl calculated dphi/dx')
    
    !Form g1=i*ky*g_wesson
    g1=g
    call g_adjust(g1,phi,bpar,fphi,fbpar)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       g1(:,:,iglo)=g1(:,:,iglo)*zi*aky(ik)
    enddo

    !Transform to real space    
    if (accelerated) then
       call transform2 (g1, agb, ia)
    else
       call transform2 (g1, gb)
    end if

    call debug_message(verb, 'nonlinear_terms::add_nl calculated dg/dy')
    
    !It should be possible to write the following with vector notation
    !To find max_vel we'd then use MAXVAL rather than MAX
    !Calculate (d Chi /dx).(d g_wesson/dy)
    if (accelerated) then
       max_vel = 0.
       do k = accelx_lo%llim_proc, accelx_lo%ulim_proc
          do j = 1, 2
             do i = 1, 2*ntgrid+1
                abracket(i,j,k) = aba(i,j,k)*agb(i,j,k)*kxfac
                max_vel = max(max_vel, abs(aba(i,j,k)))
             end do
          end do
       end do
       max_vel = max_vel * cfly
!       max_vel = maxval(abs(aba)*cfly)
    else
       max_vel = 0.
       do j = yxf_lo%llim_proc, yxf_lo%ulim_proc
          do i = 1, yxf_lo%ny
             bracket(i,j) = ba(i,j)*gb(i,j)*kxfac
             max_vel = max(max_vel,abs(ba(i,j)))
          end do
       end do
       max_vel = max_vel*cfly
!       max_vel = maxval(abs(ba)*cfly)
    endif

    call debug_message(verb, 'nonlinear_terms::add_nl calculated dphi/dx.dg/dy')

    !Form g1=i*ky*chi
    if (fphi > zero) then
       call load_ky_phi
    else
       g1 = 0.
    end if
    
    if (fbpar > zero) call load_ky_bpar
    if (fapar  > zero) call load_ky_apar

    !Transform to real space    
    if (accelerated) then
       call transform2 (g1, aba, ia)
    else
       call transform2 (g1, ba)
    end if

    !Form g1=i*kx*g_wesson
    g1=g
    call g_adjust(g1,phi,bpar,fphi,fbpar)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       g1(:,:,iglo)=g1(:,:,iglo)*zi*akx(it)
    enddo
    
    !Transform to real space
    if (accelerated) then
       call transform2 (g1, agb, ia)
    else
       call transform2 (g1, gb)
    end if

    !It should be possible to write the following with vector notation
    !To find max_vel we'd then use MAXVAL rather than MAX   
    !Calculate (d Chi /dy).(d g_wesson/dx) and subtract from (d Chi /dx).(d g_wesson/dy)
    if (accelerated) then
       do k = accelx_lo%llim_proc, accelx_lo%ulim_proc
          do j = 1, 2
             do i = 1, 2*ntgrid+1
                abracket(i,j,k) = abracket(i,j,k) - aba(i,j,k)*agb(i,j,k)*kxfac
                max_vel = max(max_vel, abs(aba(i,j,k))*cflx)
             end do
          end do
       end do
    else
       do j = yxf_lo%llim_proc, yxf_lo%ulim_proc
          do i = 1, yxf_lo%ny
             bracket(i,j) = bracket(i,j) - ba(i,j)*gb(i,j)*kxfac
             max_vel = max(max_vel,abs(ba(i,j))*cflx)
          end do
       end do
    end if

    !Estimate the global cfl limit based on max_vel
    call max_allreduce(max_vel)    
    dt_cfl = 1./max_vel
    call save_dt_cfl (dt_cfl)

    !Now check to see if we've violated the 
    !cfl condition if requested
    if(immediate_reset)then
       call check_time_step_too_large(reset)

       !If we have violated cfl then return immediately
       if(reset)return
    endif

    !Transform NL source back to spectral space
    if (accelerated) then
       call inverse2 (abracket, g1, ia)
    else
       call inverse2 (bracket, g1)
    end if
    
  contains
    !NOTE: These routines don't contain anon(ie) factor which may be desired to
    !      allow more general cases in the future.
    subroutine load_kx_phi

      use dist_fn_arrays, only: aj0
      complex :: fac

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         it = it_idx(g_lo,iglo)
         ik = ik_idx(g_lo,iglo)
         do ig = -ntgrid, ntgrid
            fac = zi*akx(it)*aj0(ig,iglo)*phi(ig,it,ik)*fphi
            g1(ig,1,iglo) = fac
            g1(ig,2,iglo) = fac
         end do
      end do

    end subroutine load_kx_phi

    subroutine load_ky_phi

      use dist_fn_arrays, only: aj0
      complex :: fac

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         it = it_idx(g_lo,iglo)
         ik = ik_idx(g_lo,iglo)
         do ig = -ntgrid, ntgrid
            fac = zi*aky(ik)*aj0(ig,iglo)*phi(ig,it,ik)*fphi
            g1(ig,1,iglo) = fac
            g1(ig,2,iglo) = fac
         end do
      end do

    end subroutine load_ky_phi

! should I use vpa or vpac in next two routines??

    subroutine load_kx_apar

      use dist_fn_arrays, only: vpa, aj0
      use gs2_layouts, only: is_idx

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         it = it_idx(g_lo,iglo)
         ik = ik_idx(g_lo,iglo)
         is = is_idx(g_lo,iglo)
         do ig = -ntgrid, ntgrid
            g1(ig,1,iglo) = g1(ig,1,iglo) - zi*akx(it)*aj0(ig,iglo)*spec(is)%stm &
                 *vpa(ig,1,iglo)*apar(ig,it,ik)*fapar 
         end do
         do ig = -ntgrid, ntgrid
            g1(ig,2,iglo) = g1(ig,2,iglo) - zi*akx(it)*aj0(ig,iglo)*spec(is)%stm &
                 *vpa(ig,2,iglo)*apar(ig,it,ik)*fapar 
         end do
      end do

    end subroutine load_kx_apar

    subroutine load_ky_apar

      use dist_fn_arrays, only: vpa, aj0
      use gs2_layouts, only: is_idx

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         it = it_idx(g_lo,iglo)
         ik = ik_idx(g_lo,iglo)
         is = is_idx(g_lo,iglo)
         do ig = -ntgrid, ntgrid
            g1(ig,1,iglo) = g1(ig,1,iglo) - zi*aky(ik)*aj0(ig,iglo)*spec(is)%stm &
                 *vpa(ig,1,iglo)*apar(ig,it,ik)*fapar 
         end do
         do ig = -ntgrid, ntgrid
            g1(ig,2,iglo) = g1(ig,2,iglo) - zi*aky(ik)*aj0(ig,iglo)*spec(is)%stm &
                 *vpa(ig,2,iglo)*apar(ig,it,ik)*fapar 
         end do
      end do

    end subroutine load_ky_apar

    subroutine load_kx_bpar

      use dist_fn_arrays, only: vperp2, aj1
      use gs2_layouts, only: is_idx
      complex :: fac

! Is this factor of two from the old normalization?

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         it = it_idx(g_lo,iglo)
         ik = ik_idx(g_lo,iglo)
         is = is_idx(g_lo,iglo)
         do ig = -ntgrid, ntgrid
            fac = g1(ig,1,iglo) + zi*akx(it)*aj1(ig,iglo) &
                 *2.0*vperp2(ig,iglo)*spec(is)%tz*bpar(ig,it,ik)*fbpar
            g1(ig,1,iglo) = fac
            g1(ig,2,iglo) = fac
         end do
      end do

    end subroutine load_kx_bpar

    subroutine load_ky_bpar

      use dist_fn_arrays, only: vperp2, aj1
      use gs2_layouts, only: is_idx
      complex :: fac

! Is this factor of two from the old normalization?

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         it = it_idx(g_lo,iglo)
         ik = ik_idx(g_lo,iglo)
         is = is_idx(g_lo,iglo)
         do ig = -ntgrid, ntgrid
            fac = g1(ig,1,iglo) + zi*aky(ik)*aj1(ig,iglo) &
                 *2.0*vperp2(ig,iglo)*spec(is)%tz*bpar(ig,it,ik)*fbpar
            g1(ig,1,iglo) = fac 
            g1(ig,2,iglo) = fac
         end do
      end do

    end subroutine load_ky_bpar

  end subroutine add_nl

  subroutine add_nl_gryfx (g1)
    use mp, only: max_allreduce
    use theta_grid, only: ntgrid, kxfac
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use gs2_layouts, only: accelx_lo, yxf_lo
    use dist_fn_arrays, only: g, g_adjust, vpa, vperp2
    use species, only: spec
    use gs2_transforms, only: transform2, inverse2
    use run_parameters, only: fapar, fbpar, fphi, reset, immediate_reset
    use kt_grids, only: aky, akx
    use gs2_time, only: save_dt_cfl, check_time_step_too_large
    use constants, only: zi
    use mp, only: broadcast
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g1
    integer :: i, j, k
    real :: max_vel, zero
    real :: dt_cfl

    integer :: iglo, ik, it, ig, iz, is, isgn, index_gryfx

    real :: densfac, uparfac, tparfac, tprpfac, qparfac, qprpfac
    densfac = 1.
    uparfac = 1./sqrt(2.)
    tparfac = 1./2.
    tprpfac = 1./2.
    qparfac = 1./sqrt(8.)
    qprpfac = 1./sqrt(8.)


    do iglo = g_lo%llim_proc, g_lo%ulim_proc
    ik = ik_idx(g_lo,iglo)
    it = it_idx(g_lo,iglo)
    is = is_idx(g_lo,iglo)
    do isgn = 1, 2
      do ig = -ntgrid, ntgrid
      iz = ig + ntgrid + 1
      !write (*,*) 'ik', ik, 'it', it, 'is', is, 'isgn', isgn, 'ig',ig
      if(ig==ntgrid) iz = 1 ! periodic point not included in gryfx arrays
        ! max is 1 + naky-1 + naky*(ntheta0-1) + naky*ntheta0*(ntheta-1),
        ! * 
        index_gryfx = 1 + (ik-1) + g_lo%naky*((it-1)) + &
                      g_lo%naky*g_lo%ntheta0*(iz-1) + &
                      (2*ntgrid)*g_lo%naky*g_lo%ntheta0*(is-1)
        !index_gryfx = 1
!        g1(ig,isgn,iglo) =  densfac*gryfx_zonal%NLdens_ky0(index_gryfx) + &
!                 (vpa(ig,isgn,iglo)*vpa(ig,isgn,iglo) - 0.5)* &
!                     tparfac*gryfx_zonal%NLtpar_ky0(index_gryfx) + &
!                 (vperp2(ig,iglo) - 1.0)* &
!                     tprpfac*gryfx_zonal%NLtprp_ky0(index_gryfx) + &
!                 2.*vpa(ig,isgn,iglo)* &
!                     uparfac*gryfx_zonal%NLupar_ky0(index_gryfx) + &
!                 (2./3.*vpa(ig,isgn,iglo)**3. - vpa(ig,isgn,iglo))* &
!                     (qparfac*gryfx_zonal%NLqpar_ky0(index_gryfx) - &
!                     3.*uparfac*gryfx_zonal%NLupar_ky0(index_gryfx) ) + &
!                 2*vpa(ig,isgn,iglo)*(vperp2(ig,iglo) - 1.)* &
!                     (qprpfac*gryfx_zonal%NLqprp_ky0(index_gryfx) - &
!                     uparfac*gryfx_zonal%NLupar_ky0(index_gryfx) )

        ! leave vpa, vperp2 in GS2 units. use momfacs to get moments in GS2 units
        ! this expression is consistent with moment definitions, such that
        ! dn = < dg >, n0 du = < vpar dg >, etc.
        ! this has been checked with a mathematica script
        g1(ig,isgn,iglo) =  densfac*gryfx_zonal%NLdens_ky0(index_gryfx) + &
                 0.5*(vpa(ig,isgn,iglo)*vpa(ig,isgn,iglo) - 1)* &
                     tparfac*gryfx_zonal%NLtpar_ky0(index_gryfx) + &
                 0.5*(vperp2(ig,iglo) - 2.0)* &
                     tprpfac*gryfx_zonal%NLtprp_ky0(index_gryfx) + &
                 vpa(ig,isgn,iglo)* &
                     uparfac*gryfx_zonal%NLupar_ky0(index_gryfx) + &
                 0.5*(vpa(ig,isgn,iglo)**3./3. - vpa(ig,isgn,iglo))* &
                     qparfac*gryfx_zonal%NLqpar_ky0(index_gryfx) + &
                 vpa(ig,isgn,iglo)*(vperp2(ig,iglo)/2. - 1.)* &
                     qprpfac*gryfx_zonal%NLqprp_ky0(index_gryfx) 
         end do
       end do
    end do
    g1 = -g1 !left-handed / right-handed conversion


    
  end subroutine add_nl_gryfx

  

  subroutine reset_init
    
    initialized = .false.
    initializing = .true.

  end subroutine reset_init

  subroutine finish_init

    initializing = .false.

  end subroutine finish_init

  subroutine finish_nonlinear_terms
    use gs2_transforms, only: finish_transforms

    implicit none

    if (allocated(aba)) deallocate (aba, agb, abracket)
    if (allocated(ba)) deallocate (ba, gb, bracket)

    nonlin = .false. ; alloc = .true. ; zip = .false. ; accelerated = .false.
    initialized = .false. ; initializing = .true.

    call finish_transforms
  end subroutine finish_nonlinear_terms
end module nonlinear_terms


