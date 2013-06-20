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
  public :: finish_init, reset_init, algorithm, nonlin
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
    use vpamu_grids, only: init_vpamu_grids, nvgrid, nmu
    use species, only: init_species, nspec
    use gs2_layouts, only: init_dist_fn_layouts, yxf_lo
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
    if (debug) write(6,*) "init_nonlinear_terms: init_vpamu_grids"
    call init_vpamu_grids
    if (debug) write(6,*) "init_nonlinear_terms: init_species"
    call init_species
    if (debug) write(6,*) "init_nonlinear_terms: init_dist_fn_layouts"
    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nvgrid, nmu, nspec)

    call read_parameters

    if (debug) write(6,*) "init_nonlinear_terms: init_transforms"
    if (nonlinear_mode_switch /= nonlinear_mode_none) then
       call init_transforms (ntgrid, naky, ntheta0, nvgrid, nmu, nspec, nx, ny)

       if (debug) write(6,*) "init_nonlinear_terms: allocations"
       if (alloc) then
          allocate (     ba(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc))
          allocate (     gb(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc))
          allocate (bracket(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc))
          ba = 0. ; gb = 0. ; bracket = 0.
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

  subroutine add_explicit_terms (g1, g2, g3, phi, apar, bpar, istep)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use gs2_time, only: save_dt_cfl
    use vpamu_grids, only: nvgrid
    implicit none
    complex, dimension (-ntgrid:,-nvgrid:,g_lo%llim_proc:), intent (in out) :: g1, g2, g3
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    integer, intent (in) :: istep
    real :: dt_cfl
    logical, save :: nl = .true.

    select case (nonlinear_mode_switch)
    case (nonlinear_mode_none)
!!! NEED TO DO SOMETHING HERE...  BD GGH
       dt_cfl = 1.e8
       call save_dt_cfl (dt_cfl)
#ifdef LOWFLOW
       if (istep /=0) &
            call add_explicit (g1, g2, g3, phi, apar, bpar, istep)
#endif
    case (nonlinear_mode_on)
       if (istep /= 0) call add_explicit (g1, g2, g3, phi, apar, bpar, istep, nl)
    end select

  end subroutine add_explicit_terms

  subroutine add_explicit (g1, g2, g3, phi, apar, bpar, istep, nl)

    use theta_grid, only: ntgrid, thet_imp
    use vpamu_grids, only: nvgrid, vpa_imp
    use gs2_layouts, only: g_lo, ik_idx, it_idx
    use dist_fn_arrays, only: g
    use gs2_time, only: save_dt_cfl
    use centering, only: get_cell_value

    implicit none

    complex, dimension (-ntgrid:,-nvgrid:,g_lo%llim_proc:), intent (in out) :: g1, g2, g3
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    integer, intent (in) :: istep
    logical, intent (in), optional :: nl

    integer :: istep_last = 0
    integer :: iglo, ik, it
    real :: dt_cfl, zero

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
          
          ! take g1 at grid points and return 2*g1 at cell centers
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             call get_cell_value (thet_imp, vpa_imp, g1(:,:,iglo), g1(:,:,iglo), -ntgrid, -nvgrid)
          end do
          g1 = 2.*g1
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
             g (:,:,iglo) = 0.
             g1(:,:,iglo) = 0.
          end if
       end do
    end if
     
    istep_last = istep

  end subroutine add_explicit

  subroutine add_nl (g1, phi, apar, bpar)

    use mp, only: max_allreduce
    use theta_grid, only: ntgrid, kxfac
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use gs2_layouts, only: yxf_lo
    use dist_fn_arrays, only: g
    use species, only: spec
    use gs2_transforms, only: transform2, inverse2
    use run_parameters, only: fapar, fbpar, fphi
    use kt_grids, only: aky, akx
    use vpamu_grids, only: nvgrid
    use gs2_time, only: save_dt_cfl
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,-nvgrid:,g_lo%llim_proc:), intent (in out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    integer :: i, j, k
    real :: max_vel, zero
    real :: dt_cfl

    integer :: iglo, ik, it, is, ig, iv
    
    !Initialise zero so we can be sure tests are sensible
    zero = epsilon(0.0)

    if (fphi > zero) then
       call load_kx_phi
    else
       g1 = 0.
    end if

    if (fbpar > zero) call load_kx_bpar
    if (fapar  > zero) call load_kx_apar

    call transform2 (g1, ba)
    
    if (fphi > zero) then
       call load_ky_phi
    else
       g1 = 0.
    end if
    if (fbpar > zero) call load_ky_bpar
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do iv = -nvgrid, nvgrid
          do ig = -ntgrid, ntgrid
             g1(ig,iv,iglo) = g1(ig,iv,iglo)*spec(is)%zt + zi*aky(ik)*g(ig,iv,iglo)
          end do
       end do
    end do
    
    call transform2 (g1, gb)
    
    max_vel = 0.
    do j = yxf_lo%llim_proc, yxf_lo%ulim_proc
       do i = 1, yxf_lo%ny
          bracket(i,j) = ba(i,j)*gb(i,j)*kxfac
          max_vel = max(max_vel,abs(ba(i,j)))
       end do
    end do
    max_vel = max_vel*cfly

    if (fphi > zero) then
       call load_ky_phi
    else
       g1 = 0.
    end if
    
    if (fbpar > zero) call load_ky_bpar
    if (fapar  > zero) call load_ky_apar
    
    call transform2 (g1, ba)
    
    if (fphi > zero) then
       call load_kx_phi
    else
       g1 = 0.
    end if
    
    if (fbpar > zero) call load_kx_bpar
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do iv = -nvgrid, nvgrid
          do ig = -ntgrid, ntgrid
             g1(ig,iv,iglo) = g1(ig,iv,iglo)*spec(is)%zt + zi*akx(it)*g(ig,iv,iglo)
          end do
       end do
    end do
    
    call transform2 (g1, gb)
    
    do j = yxf_lo%llim_proc, yxf_lo%ulim_proc
       do i = 1, yxf_lo%ny
          bracket(i,j) = bracket(i,j) - ba(i,j)*gb(i,j)*kxfac
          max_vel = max(max_vel,abs(ba(i,j))*cflx)
       end do
    end do
    
    call max_allreduce(max_vel)
    
    dt_cfl = 1./max_vel
    call save_dt_cfl (dt_cfl)
    
    call inverse2 (bracket, g1)
    
  contains

    subroutine load_kx_phi

      use dist_fn_arrays, only: aj0
      complex :: fac

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         it = it_idx(g_lo,iglo)
         ik = ik_idx(g_lo,iglo)
         do ig = -ntgrid, ntgrid
            fac = zi*akx(it)*aj0(ig,iglo)*phi(ig,it,ik)*fphi
            g1(ig,:,iglo) = fac
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
            g1(ig,:,iglo) = fac
         end do
      end do

    end subroutine load_ky_phi

    subroutine load_kx_apar

      use dist_fn_arrays, only: aj0
      use gs2_layouts, only: is_idx
      use vpamu_grids, only: vpa

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         it = it_idx(g_lo,iglo)
         ik = ik_idx(g_lo,iglo)
         is = is_idx(g_lo,iglo)
         do ig = -ntgrid, ntgrid
            g1(ig,:,iglo) = g1(ig,:,iglo) - zi*akx(it)*aj0(ig,iglo)*spec(is)%stm &
                 *vpa*apar(ig,it,ik)*fapar 
         end do
      end do

    end subroutine load_kx_apar

    subroutine load_ky_apar

      use dist_fn_arrays, only: aj0
      use vpamu_grids, only: vpa
      use gs2_layouts, only: is_idx

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         it = it_idx(g_lo,iglo)
         ik = ik_idx(g_lo,iglo)
         is = is_idx(g_lo,iglo)
         do ig = -ntgrid, ntgrid
            g1(ig,:,iglo) = g1(ig,:,iglo) - zi*aky(ik)*aj0(ig,iglo)*spec(is)%stm &
                 *vpa*apar(ig,it,ik)*fapar 
         end do
      end do

    end subroutine load_ky_apar

    subroutine load_kx_bpar

      use dist_fn_arrays, only: aj1
      use gs2_layouts, only: is_idx, ik_idx, it_idx, imu_idx
      use vpamu_grids, only: mu, vperp2

      integer :: it, ik, imu, is
      complex :: fac

! Is this factor of two from the old normalization?

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         it = it_idx(g_lo,iglo)
         ik = ik_idx(g_lo,iglo)
         is = is_idx(g_lo,iglo)
         imu = imu_idx(g_lo,iglo)
         do ig = -ntgrid, ntgrid
            fac = g1(ig,1,iglo) + zi*akx(it)*aj1(ig,iglo) &
                 *2.0*vperp2(ig,imu)*spec(is)%tz*bpar(ig,it,ik)*fbpar
            g1(ig,1,iglo) = fac
            g1(ig,2,iglo) = fac
         end do
      end do

    end subroutine load_kx_bpar

    subroutine load_ky_bpar

      use dist_fn_arrays, only: aj1
      use gs2_layouts, only: is_idx, ik_idx, it_idx, imu_idx
      use vpamu_grids, only: vperp2
      use kt_grids, only: aky

      integer :: imu, it, ik, is, ig
      complex :: fac

! Is this factor of two from the old normalization?

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         it = it_idx(g_lo,iglo)
         ik = ik_idx(g_lo,iglo)
         is = is_idx(g_lo,iglo)
         do ig = -ntgrid, ntgrid
            fac = g1(ig,1,iglo) + zi*aky(ik)*aj1(ig,iglo) &
                 *2.0*vperp2(ig,imu)*spec(is)%tz*bpar(ig,it,ik)*fbpar
            g1(ig,1,iglo) = fac 
            g1(ig,2,iglo) = fac
         end do
      end do

    end subroutine load_ky_bpar

  end subroutine add_nl

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

    if (allocated(ba)) deallocate (ba, gb, bracket)

    nonlin = .false. ; alloc = .true. ; zip = .false.
    initialized = .false. ; initializing = .true.

  end subroutine finish_nonlinear_terms

end module nonlinear_terms


