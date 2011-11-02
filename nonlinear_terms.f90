module nonlinear_terms
!
! not correct for slowing down species
! terms like df/dU should involve anon, etc.
!
! missing factors of B_a/B(theta) in A_perp terms??
!
  implicit none

  public :: init_nonlinear_terms, finish_nonlinear_terms
  public :: read_parameters, wnml_nonlinear_terms, check_nonlinear_terms
!  public :: add_nonlinear_terms, finish_nl_terms
  public :: add_explicit_terms, finish_nl_terms
  public :: finish_init, reset_init, algorithm, nonlin, accelerated
  public :: cfl

  private

  ! knobs
  integer, public :: nonlinear_mode_switch
  integer :: flow_mode_switch

  integer, public, parameter :: nonlinear_mode_none = 1, nonlinear_mode_on = 2
  integer, public, parameter :: flow_mode_off = 1, flow_mode_on = 2

  !complex, dimension(:,:), allocatable :: phi_avg, apar_avg, bpar_avg  

  real, save, dimension (:,:), allocatable :: ba, gb, bracket
  ! yxf_lo%ny, yxf_lo%llim_proc:yxf_lo%ulim_alloc

  real, save, dimension (:,:,:), allocatable :: aba, agb, abracket
  ! 2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc

  !complex, dimension (:,:), allocatable :: xax, xbx, g_xf
  ! xxf_lo%nx, xxf_lo%llim_proc:xxf_lo%ulim_alloc

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
  logical :: accelerated = .false.

  logical :: exist
  
contains
  
  subroutine check_nonlinear_terms(report_unit,delt_adj)
    use gs2_time, only: code_dt_min
    use kt_grids, only: box
    use run_parameters, only: margin, code_delt_max, nstep, wstar_units
    use theta_grid, only: nperiod
    implicit none
    integer :: report_unit
    real :: delt_adj
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
       write (report_unit, fmt="('The minimum delt ( code_dt_min ) = ',e10.4)") code_dt_min
       write (report_unit, fmt="('The maximum delt (code_delt_max) = ',e10.4)") code_delt_max
       write (report_unit, fmt="('The maximum delt < ',f10.4,' * min(Delta_perp/v_perp). (cfl)')") cfl
       write (report_unit, fmt="('When the time step needs to be changed, it is adjusted by a factor of ',f10.4)") delt_adj
       write (report_unit, fmt="('The number of time steps nstep = ',i7)") nstep
       write (report_unit, fmt="('If running in batch mode on the NERSC T3E, the run will stop when ', &
            & f6.4,' % of the time remains.')") 100.*margin
    endif
  end subroutine check_nonlinear_terms


  subroutine wnml_nonlinear_terms(unit)
  implicit none
  integer :: unit
    if (.not. exist) return
    if (nonlinear_mode_switch == nonlinear_mode_on) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "nonlinear_terms_knobs"
       write (unit, fmt="(' nonlinear_mode = ',a)") '"on"'
       write (unit, fmt="(' cfl = ',e16.10)") cfl
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
    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec)

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
         C_par, C_perp, p_x, p_y, p_z, zip
    integer :: ierr, in_file
!    logical :: done = .false.

!    if (done) return
!    done = .true.

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
            ierr, "nonlinear_mode in nonlinear_terms_knobs")
       call get_option_value &
            (flow_mode, flowopts, flow_mode_switch, &
            ierr, "flow_mode in nonlinear_terms_knobs")
    end if

    call broadcast (nonlinear_mode_switch)
    call broadcast (flow_mode_switch)
    call broadcast (cfl)
    call broadcast (C_par) 
    call broadcast (C_perp) 
    call broadcast (p_x)
    call broadcast (p_y)
    call broadcast (p_z)
    call broadcast (zip)

    if (flow_mode_switch == flow_mode_on) then
       if (proc0) write(*,*) 'Forcing flow_mode = off'
       flow_mode_switch = flow_mode_off
    endif

    if (nonlinear_mode_switch == nonlinear_mode_on)  then
       algorithm = 1 
       nonlin = .true.
    end if

  end subroutine read_parameters

!  subroutine add_nonlinear_terms (g1, g2, g3, phi, apar, bpar, istep, bd, fexp)
  subroutine add_explicit_terms (g1, g2, g3, phi, apar, bpar, istep, bd, fexp)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use gs2_time, only: save_dt_cfl
    use run_parameters, only: include_lowflow
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1, g2, g3
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    integer, intent (in) :: istep
    real, intent (in) :: bd
    complex, intent (in) :: fexp
    real :: dt_cfl
    logical :: nl = .true.

    select case (nonlinear_mode_switch)
    case (nonlinear_mode_none)
!!! NEED TO DO SOMETHING HERE...  BD GGH
       dt_cfl = 1.e8
       call save_dt_cfl (dt_cfl)
       if (include_lowflow .and. istep /= 0) &
            call add_explicit (g1, g2, g3, phi, apar, bpar, istep, bd, fexp)
    case (nonlinear_mode_on)
!       if (istep /= 0) call add_nl (g1, g2, g3, phi, apar, bpar, istep, bd, fexp)
       if (istep /= 0) call add_explicit (g1, g2, g3, phi, apar, bpar, istep, bd, fexp, nl)
    end select
  end subroutine add_explicit_terms
!  end subroutine add_nonlinear_terms

!  subroutine add_nl (g1, g2, g3, phi, apar, bpar, istep, bd, fexp)
  subroutine add_explicit (g1, g2, g3, phi, apar, bpar, istep, bd, fexp, nl)
    use mp, only: max_allreduce, proc0
    use theta_grid, only: ntgrid, kxfac
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use gs2_layouts, only: accelx_lo, yxf_lo
    use dist_fn_arrays, only: g, ittp
    use species, only: spec
    use gs2_transforms, only: transform2, inverse2
    use run_parameters, only: fapar, fbpar, fphi, include_lowflow
    use kt_grids, only: aky, akx
    use le_grids, only: forbid
    use gs2_time, only: save_dt_cfl
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1, g2, g3
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    integer, intent (in) :: istep
    real, intent (in) :: bd
    complex, intent (in) :: fexp
    logical, intent (in), optional :: nl
    integer :: istep_last = 0
    integer :: i, j, k
    real :: max_vel, zero
    real :: dt_cfl

    integer :: iglo, ik, it, is, ig, il, ia, isgn
    
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

    if (istep /= istep_last) then

       zero = epsilon(0.0)
       g3 = g2
       g2 = g1

       ! if running nonlinearly, then compute the nonlinear term at grid points
       ! and store it in g1
       if (present(nl)) then
          call add_nl (g1, phi, apar, bpar)
          
          ! takes g1 at grid points and returns 2*g1 at cell centers
          call center (g1)
       else
          g1 = 0.
       end if

       ! if including lowflow terms, add in the explicit ones at grid points
       if (include_lowflow) call add_lf (g1)

    endif

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

      use dist_fn_arrays, only: ittp
      use gs2_layouts, only: g_lo, il_idx
      use theta_grid, only: ntgrid
      use le_grids, only: forbid

      implicit none

      integer :: iglo, il, ig
      complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: gtmp

      ! factor of one-half appears elsewhere
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          il = il_idx(g_lo, iglo)
          ! Totally trapped particles get no bakdif 
          do ig = -ntgrid, ntgrid-1
             !JAB & GWH: orig if (il == ittp(ig)) cycle 
             if (il >= ittp(ig)) cycle
             gtmp(ig,1,iglo) = (1.+bd)*gtmp(ig+1,1,iglo) + (1.-bd)*gtmp(ig,1,iglo)
             gtmp(ig,2,iglo) = (1.-bd)*gtmp(ig+1,2,iglo) + (1.+bd)*gtmp(ig,2,iglo)
          end do
          ! zero out spurious gtmp outside trapped boundary
          where (forbid(:,il))
             gtmp(:,1,iglo) = 0.0
             gtmp(:,2,iglo) = 0.0
          end where
       end do

    end subroutine center

  end subroutine add_explicit
!  end subroutine add_nl

  subroutine add_nl (g1, phi, apar, bpar)
  
    use mp, only: max_allreduce, proc0
    use theta_grid, only: ntgrid, kxfac
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use gs2_layouts, only: accelx_lo, yxf_lo
    use dist_fn_arrays, only: g, ittp
    use species, only: spec
    use gs2_transforms, only: transform2, inverse2
    use run_parameters, only: fapar, fbpar, fphi
    use kt_grids, only: aky, akx
    use le_grids, only: forbid
    use gs2_time, only: save_dt_cfl
    use constants, only: zi
    implicit none
    integer :: i, j, k
    real :: max_vel, zero
    real :: dt_cfl
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    
    integer :: iglo, ik, it, is, ig, il, ia, isgn
    
    if (fphi > zero) then
       call load_kx_phi
    else
       g1 = 0.
    end if
    
    if (fbpar > zero) call load_kx_bpar
    if (fapar  > zero) call load_kx_apar
    
    if (accelerated) then
       call transform2 (g1, aba, ia)
    else
       call transform2 (g1, ba)
    end if
    
    if (fphi > zero) then
       call load_ky_phi
    else
       g1 = 0.
    end if
    if (fbpar > zero) call load_ky_bpar
    
    ! more generally, there should probably be a factor of anon...
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             g1(ig,isgn,iglo) = g1(ig,isgn,iglo)*spec(is)%zt + zi*aky(ik)*g(ig,isgn,iglo)
          end do
       end do
    end do
    
    if (accelerated) then
       call transform2 (g1, agb, ia)
    else
       call transform2 (g1, gb)
    end if
    
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
!          max_vel = maxval(abs(aba)*cfly)
    else
       max_vel = 0.
       do j = yxf_lo%llim_proc, yxf_lo%ulim_proc
          do i = 1, yxf_lo%ny
             bracket(i,j) = ba(i,j)*gb(i,j)*kxfac
             max_vel = max(max_vel,abs(ba(i,j)))
          end do
       end do
       max_vel = max_vel*cfly
!          max_vel = maxval(abs(ba)*cfly)
    endif

    if (fphi > zero) then
       call load_ky_phi
    else
       g1 = 0.
    end if
    
    if (fbpar > zero) call load_ky_bpar
    if (fapar  > zero) call load_ky_apar
    
    if (accelerated) then
       call transform2 (g1, aba, ia)
    else
       call transform2 (g1, ba)
    end if
    
    if (fphi > zero) then
       call load_kx_phi
    else
       g1 = 0.
    end if
    
    if (fbpar > zero) call load_kx_bpar
    
    ! more generally, there should probably be a factor of anon...
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             g1(ig,isgn,iglo) = g1(ig,isgn,iglo)*spec(is)%zt + zi*akx(it)*g(ig,isgn,iglo)
          end do
       end do
    end do
    
    if (accelerated) then
       call transform2 (g1, agb, ia)
    else
       call transform2 (g1, gb)
    end if
    
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
    
    call max_allreduce(max_vel)
    
    dt_cfl = 1./max_vel
    call save_dt_cfl (dt_cfl)
    
    if (accelerated) then
       call inverse2 (abracket, g1, ia)
    else
       call inverse2 (bracket, g1)
    end if

  contains

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

  ! add_lf subroutine computes contributions from lowflow terms that are
  ! treated explicitly and adds them to the explicit nonlinear terms
  subroutine add_lf (g1)
  
    use dist_fn_arrays, only: vpac
    use gs2_layouts, only: g_lo, is_idx, ie_idx
    use lowflow, only: dphidth
    use species, only: spec
    use theta_grid, only: ntgrid, gradpar
    use le_grids, only: energy

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g1

    integer :: ig, is, isgn, iglo, ie
    complex, dimension (:,:,:), allocatable :: dgde

    allocate (dgde(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))

    ! compute dg/dE_{\mu}, where E is energy and \mu is magnetic moment
    call get_dgde (dgde)

    ! this is the contribution from the vpar . grad (phi^{nc}) dg/dE term at cell centers
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2               
          dgde(-ntgrid:ntgrid-1,isgn,iglo) = spec(is)%zstm &
               *0.125*(gradpar(-ntgrid+1:)+gradpar(:ntgrid-1)) &
               *(dphidth(2:)+dphidth(:2*ntgrid)) &
               *(dgde(-ntgrid+1:,isgn,iglo)+dgde(:ntgrid-1,isgn,iglo))
          dgde(ntgrid,isgn,iglo) = 0.
       end do
    end do

    call write_mpdist (dgde)

    ! add explicit lowflow terms to explicit nonlinear terms
    g1 = g1 + dgde

    deallocate (dgde)

  end subroutine add_lf

  subroutine get_dgde (dgde)

    use redistribute, only: gather, scatter
    use dist_fn_arrays, only: g
    use gs2_layouts, only: g_lo, le_lo, ig_idx
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: negrid, nlambda, ng2, energy, le_map, jend
    use le_grids, only: al

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: dgde

    integer :: ilelo, ie, nxi, ixi, je, te, ig
    real, dimension (:), allocatable :: xi, de, dxi
    complex, dimension (:,:,:), allocatable :: gle, dgle, dgle2

    nxi = max(2*nlambda-1,2*ng2)

    allocate (de(negrid-1))

    de = energy(2:)-energy(:negrid-1)

    allocate (gle(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
    allocate (dgle(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
    allocate (dgle2(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))

    ! set dgde = g
    dgde = g
    ! rearrange data so that all energies and pitch-angles are on same processor
    call gather (le_map, dgde, gle)

    ! compute dg/de and store it in ged
    do ilelo = le_lo%llim_proc, le_lo%ulim_proc

       ig = ig_idx(le_lo, ilelo)
       ! at this theta value, je is the number of grid points with vpar >= 0
       je = jend(ig)
       ! te is total number of velocity grid points at this theta value
       ! minus one avoids double counting of vpar = 0
       te = 2*jend(ig)-1

       allocate (xi(te+1))
       allocate (dxi(te-1))

       xi = 0.

       ! fill in negative pitch-angles for this theta value
       do ixi = 1, je
          xi(ixi) = sqrt(1.-al(ixi)*bmag(ig))
       end do
       ! fill in posative pitch-angles for this theta value
       xi(je+1:te) = -xi(je-1:1:-1)
       xi(te+1)=0.

       dxi = xi(2:te)-xi(:te-1)

       ! 3-pt finite difference at interior points (2nd order accurate)
       ! this is dg/dE_{\xi}
       dgle(:te+1,2:negrid-1,ilelo) = ((gle(:te+1,3:negrid,ilelo)-gle(:te+1,2:negrid-1,ilelo)) &
            * spread(de(:negrid-2)/de(2:), 1, te+1) &
            + (gle(:te+1,2:negrid-1,ilelo)-gle(:te+1,:negrid-2,ilelo)) &
            * spread(de(2:)/de(:negrid-2), 1, te+1)) &
            / spread(de(:negrid-2)+de(2:), 1, te+1)

       ! 2-pt compact finite difference at boundaries (2nd order accurate)
       dgle(:te+1,negrid,ilelo) = 2.*(gle(:te+1,negrid,ilelo)-gle(:te+1,negrid-1,ilelo)) &
            / spread(de(negrid-1),1,te+1) - dgle(:te+1,negrid-1,ilelo)
       dgle(:te+1,1,ilelo) = 2.*(gle(:te+1,2,ilelo)-gle(:te+1,1,ilelo)) &
            / spread(de(1),1,te+1) - dgle(:te+1,2,ilelo)

       ! need vpa*((dg/dE)_{\xi} - g) because g normalized by exp(-E)
       ! chain rules gives additional -g
       dgle(:te+1,:negrid,ilelo) = spread(xi, 2, negrid)*spread(sqrt(energy), 1, te+1) &
            * (dgle(:te+1,:negrid,ilelo)-gle(:te+1,:negrid,ilelo))

       ! 3-pt finite differences at interior points (2nd order accurate)
       ! this is dg/dxi
       dgle2(2:te-1,:negrid,ilelo) = ((gle(3:te,:negrid,ilelo)-gle(2:te-1,:negrid,ilelo)) &
            * spread(dxi(:te-2)/dxi(2:te-1), 2, negrid) &
            + (gle(2:te-1,:negrid,ilelo)-gle(:te-2,:negrid,ilelo)) &
            * spread(dxi(2:te-1)/dxi(:te-2), 2, negrid)) &
            / spread(dxi(:te-2)+dxi(2:te-1), 2, negrid)

       ! 2-pt compact finite difference at boundaries (2nd order accurate)
       dgle2(te,:negrid,ilelo) = 2.*(gle(te,:negrid,ilelo)-gle(te-1,:negrid,ilelo)) &
            / spread(dxi(te-1),1,negrid) - dgle2(te-1,:negrid,ilelo)
       dgle2(1,:negrid,ilelo) = 2.*(gle(2,:negrid,ilelo)-gle(1,:negrid,ilelo)) &
            / spread(dxi(1),1,negrid) - dgle2(2,:negrid,ilelo)
       ! deal with vpa = 0 point (ixi=je), which is duplicated at ixi=te+1 for mapping
       dgle2(te+1,:negrid,ilelo) = dgle2(je,:negrid,ilelo)

       ! vpa*dg/dE_{\mu} = vpa*dg/dE_{\xi} + (1-xi^2)/(2 sqrt(E)) * dg/dxi
       dgle(:te+1,:negrid,ilelo) = dgle(:te+1,:negrid,ilelo) &
            + (1.-spread(xi**2,2,negrid))*0.5/spread(sqrt(energy),1,te+1)*dgle2(:te+1,:negrid,ilelo)

       deallocate (xi, dxi)

    end do

    ! revert back to original arrangement of data over processors
    call scatter (le_map, dgle, dgde)

!    call write_mpdist (dgde)

    deallocate (gle, dgle, dgle2, de)

  end subroutine get_dgde

  subroutine write_mpdist (dist)

    use mp, only: proc0, send, receive
    use file_utils, only: open_output_file, close_output_file
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use gs2_layouts, only: ie_idx, il_idx, idx_local, proc_id
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: forbid, energy, al

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:) :: dist

    integer :: iglo, ik, it, is, ie, il, ig
    integer, save :: unit
    logical :: first = .true.
    real :: vpa1
    complex, dimension (2) :: gtmp

     if (first) then
        if (proc0) call open_output_file (unit, ".distmp")
        do iglo=g_lo%llim_world, g_lo%ulim_world
           ik = ik_idx(g_lo, iglo)
           it = it_idx(g_lo, iglo)
           is = is_idx(g_lo, iglo)
           ie = ie_idx(g_lo, iglo)
           il = il_idx(g_lo, iglo)
           do ig = -ntgrid, ntgrid
              if (idx_local (g_lo, ik, it, il, ie, is)) then
                 if (proc0) then
                    gtmp = dist(ig,:,iglo)
                 else
                    call send (dist(ig,:,iglo), 0)
                 end if
              else if (proc0) then
                 call receive (gtmp, proc_id(g_lo, iglo))
              end if
              if (proc0) then
                 if (.not. forbid(ig,il)) then
                    vpa1 = sqrt(energy(ie)*max(0.0,1.0-al(il)*bmag(ig)))
                    write (unit,'(6e14.5)') vpa1, energy(ie), &
                         real(gtmp(1)), aimag(gtmp(1)), real(gtmp(2)), aimag(gtmp(2))
                 end if
              end if
           end do
        end do
        if (proc0) call close_output_file (unit)
        first = .false.
     end if

   end subroutine write_mpdist

  subroutine finish_nl_terms

    if (nonlinear_mode_switch == nonlinear_mode_none) return
!    deallocate (ba)
!    deallocate (gb)
!    deallocate (bracket)
!    alloc = .true.

  end subroutine finish_nl_terms

  subroutine reset_init
    
    initialized = .false.
    initializing = .true.
    call finish_nl_terms

  end subroutine reset_init

  subroutine finish_init

    initializing = .false.

  end subroutine finish_init

  subroutine finish_nonlinear_terms

    implicit none

    if (allocated(aba)) deallocate (aba, agb, abracket)
    if (allocated(ba)) deallocate (ba, gb, bracket)

    nonlin = .false. ; alloc = .true. ; zip = .false. ; accelerated = .false.
    initialized = .false. ; initializing = .true.

  end subroutine finish_nonlinear_terms

end module nonlinear_terms


