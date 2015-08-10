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
  logical :: test_nonlinear

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
       write (report_unit, fmt="('The minimum delt ( code_dt_min ) = ',e11.4)") code_dt_min
       write (report_unit, fmt="('The maximum delt (code_delt_max) = ',e11.4)") code_delt_max
       write (report_unit, fmt="('The maximum delt < ',f10.4,' * min(Delta_perp/v_perp). (cfl)')") cfl
       write (report_unit, fmt="('When the time step needs to be changed, it is adjusted by a factor of ',f10.4)") delt_adj
       write (report_unit, fmt="('The number of time steps nstep = ',i7)") nstep
       write (report_unit, fmt="('If running in batch mode on the NERSC T3E, the run will stop when ', &
            & f11.4,' % of the time remains.')") 100.*margin
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
       write (unit, fmt="(' cfl = ',e17.10)") cfl
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
         C_par, C_perp, p_x, p_y, p_z, zip, test_nonlinear
    integer :: ierr, in_file
!    logical :: done = .false.

!    if (done) return
!    done = .true.

    if (proc0) then
       nonlinear_mode = 'default'
       flow_mode = 'default'
       test_nonlinear = .false.
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
    call broadcast (test_nonlinear)
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
    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in out) :: g1, g2, g3
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi, apar, bpar
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

    use mp, only: mp_abort
    use theta_grid, only: ntgrid, thet_imp
    use vpamu_grids, only: nvgrid, vpa_imp
    use gs2_layouts, only: g_lo, ik_idx
    use dist_fn_arrays, only: g
    use gs2_time, only: save_dt_cfl
    use centering, only: get_cell_value
    use kt_grids, only: ntheta0

    implicit none

    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in out) :: g1, g2, g3
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi, apar, bpar
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

    ! Implicit terms must be calculated twice per time step
    ! due to nature of Beam-Warming algorithm.
    ! However, explicit terms such as nonlinearity need only 
    ! be computed once.    
    if (istep /= istep_last) then

       zero = epsilon(0.0)
       g3 = g2
       g2 = g1

       ! if running nonlinearly, then compute the nonlinear term at grid points
       ! and store it in g1
       if (present(nl)) then
          call add_nl (g1, phi, apar, bpar)
          
          ! take g1 at grid points and return g1 at cell centers
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             do it = 1, ntheta0
                call get_cell_value (thet_imp, vpa_imp, g1(:,:,it,iglo), g1(:,:,it,iglo), &
                     -ntgrid, -nvgrid)
             end do
          end do

          if (test_nonlinear) then
             call write_mpdist (g1, '.g1')
             call mp_abort ('testing calculation of nonlinear term. results in .g1 file')
          end if

       else
          g1 = 0.
       end if

#ifdef LOWFLOW
       ! do something
#endif

    end if

    if (zip) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          if (ik /= 1) then
             g (:,:,:,iglo) = 0.
             g1(:,:,:,iglo) = 0.
          end if
       end do
    end if
     
    istep_last = istep

  end subroutine add_explicit

  subroutine add_nl (g1, phi, apar, bpar)

    use mp, only: max_allreduce, mp_abort
    use theta_grid, only: ntgrid, kxfac, theta, shat
    use theta_grid, only: delthet, gds23, gds24_noq
    use gs2_layouts, only: g_lo, ik_idx, is_idx, imu_idx, iv_idx, ig_idx
    use gs2_layouts, only: yxf_lo
    use dist_fn_arrays, only: g
    use species, only: spec
    use gs2_transforms, only: transform2, inverse2
    use run_parameters, only: fapar, fbpar, fphi
    use kt_grids, only: aky, akx, ntheta0
    use vpamu_grids, only: nvgrid, vpa, anon
    use gs2_time, only: save_dt_cfl
    use constants, only: zi, pi

    implicit none

    integer :: i, j, iglo, ik, it, is, ig, iv, imu
    real :: max_vel, zero
    real :: dt_cfl

    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in out) :: g1
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi, apar, bpar

    !Initialise zero so we can be sure tests are sensible
    zero = epsilon(0.0)

    if (test_nonlinear) call init_test (phi, g)

    if (fphi > zero) then
       call load_kx_phi (phi, g1)
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
             g1(ig,iv,:,iglo) = g1(ig,iv,:,iglo)*spec(is)%zt + zi*aky(ik)*g(ig,iv,:,iglo)
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
       call load_kx_phi (phi, g1)
    else
       g1 = 0.
    end if
    
    if (fbpar > zero) call load_kx_bpar
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       do it = 1, ntheta0
          do iv = -nvgrid, nvgrid
             do ig = -ntgrid, ntgrid
                g1(ig,iv,it,iglo) = g1(ig,iv,it,iglo)*spec(is)%zt + zi*akx(it)*g(ig,iv,it,iglo)
             end do
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

#ifdef LOWFLOW
    if (fphi > zero) then
       call load_dthet_phi (phi, g1)
    else
       g1 = 0.
    end if

    call transform2 (g1, ba)
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       do it = 1, ntheta0
          do iv = -nvgrid, nvgrid
             g1(:,iv,it,iglo) = zi*(akx(it)*gds24_noq-aky(ik)*gds23)*g(:,iv,it,iglo)
          end do
       end do
    end do
    
    call transform2 (g1, gb)
    
    do j = yxf_lo%llim_proc, yxf_lo%ulim_proc
       do i = 1, yxf_lo%ny
          bracket(i,j) = bracket(i,j) + ba(i,j)*gb(i,j)
! assume rhostar terms do not set max_vel
!          max_vel = max(max_vel,abs(ba(i,j)))
       end do
    end do

    if (fphi > zero) then
       call load_kxky_phi (phi, g1)
    else
       g1 = 0.
    end if

    call transform2 (g1, ba)
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do it = 1, ntheta0
          do iv = -nvgrid, nvgrid
             g1(:ntgrid-1,iv,it,iglo) = (g(-ntgrid+1:,iv,it,iglo)-g(:ntgrid-1,iv,it,iglo))/delthet(:ntgrid-1)
          end do
       end do
    end do
    
    call transform2 (g1, gb)

    do j = yxf_lo%llim_proc, yxf_lo%ulim_proc
       do i = 1, yxf_lo%ny
          bracket(i,j) = bracket(i,j) - ba(i,j)*gb(i,j)
! assume max_vel is not set by terms small in rhostar
!          max_vel = max(max_vel,abs(ba(i,j))*cflx)
       end do
    end do    
#endif

    call max_allreduce(max_vel)
    
    dt_cfl = 1./max_vel
    call save_dt_cfl (dt_cfl)

    call inverse2 (bracket, g1)

  contains

    subroutine init_test (phitest, gtest)

      use theta_grid, only: ntgrid
      use vpamu_grids, only: nvgrid
      use gs2_layouts, only: g_lo, ik_idx

      implicit none

      complex, dimension (-ntgrid:,:,:), intent (out) :: phitest
      complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (out) :: gtest

      integer :: iglo, ik, it

      ! hardwire g to be cos(dkx * x + dky * y)
      gtest = 0.
      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         ik = ik_idx(g_lo,iglo) ; if (ik /= 2) cycle 
         it = 2
         gtest(:,:,it,iglo) = 1.0
      end do
      
      ! hardwire phi to be cos(2dkx * x + dky * y)
      phitest = 0.
      ik = 2 ; it = 3
      phitest(:,it,ik) = 1.0
      
    end subroutine init_test

    subroutine load_kx_phi (phi, g1)

      use constants, only: zi
      use dist_fn_arrays, only: aj0
      use gs2_layouts, only: g_lo, ik_idx
      use kt_grids, only: akx, ntheta0
      use theta_grid, only: ntgrid
      use run_parameters, only: fphi
      use vpamu_grids, only: nvgrid

      implicit none

      complex, dimension (-ntgrid:,:,:), intent (in) :: phi
      complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (out) :: g1

      integer :: iglo, ik, it, ig

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         ik = ik_idx(g_lo,iglo)
         do iv = -nvgrid, nvgrid
            do ig = -ntgrid, ntgrid
               g1(ig,iv,:,iglo) = zi*akx*aj0(ig,:,iglo)*phi(ig,:,ik)*fphi
            end do
         end do
      end do

    end subroutine load_kx_phi

    subroutine load_ky_phi

      use dist_fn_arrays, only: aj0
      complex :: fac

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         ik = ik_idx(g_lo,iglo)
         do it = 1, ntheta0
            do ig = -ntgrid, ntgrid
               fac = zi*aky(ik)*aj0(ig,it,iglo)*phi(ig,it,ik)*fphi
               g1(ig,:,it,iglo) = fac
            end do
         end do
      end do

    end subroutine load_ky_phi

    subroutine load_kx_apar

      use dist_fn_arrays, only: aj0
      use gs2_layouts, only: is_idx
      use vpamu_grids, only: vpa

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         ik = ik_idx(g_lo,iglo)
         is = is_idx(g_lo,iglo)
         do it = 1, ntheta0
            do ig = -ntgrid, ntgrid
               g1(ig,:,it,iglo) = g1(ig,:,it,iglo) - zi*akx(it)*aj0(ig,it,iglo)*spec(is)%stm &
                    *vpa*apar(ig,it,ik)*fapar 
            end do
         end do
      end do

    end subroutine load_kx_apar

    subroutine load_ky_apar

      use dist_fn_arrays, only: aj0
      use vpamu_grids, only: vpa
      use gs2_layouts, only: is_idx

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         ik = ik_idx(g_lo,iglo)
         is = is_idx(g_lo,iglo)
         do it = 1, ntheta0
            do ig = -ntgrid, ntgrid
               g1(ig,:,it,iglo) = g1(ig,:,it,iglo) - zi*aky(ik)*aj0(ig,it,iglo)*spec(is)%stm &
                    *vpa*apar(ig,it,ik)*fapar 
            end do
         end do
      end do

    end subroutine load_ky_apar

    subroutine load_kx_bpar

      use dist_fn_arrays, only: aj1
      use gs2_layouts, only: is_idx, ik_idx, imu_idx
      use vpamu_grids, only: mu, vperp2

      integer :: it, ik, imu, is
      complex :: fac

! Is this factor of two from the old normalization?

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         ik = ik_idx(g_lo,iglo)
         is = is_idx(g_lo,iglo)
         imu = imu_idx(g_lo,iglo)
         do it = 1, ntheta0
            do ig = -ntgrid, ntgrid
               g1(ig,:,it,iglo) = g1(ig,:,it,iglo) + zi*akx(it)*aj1(ig,it,iglo) &
                    *2.0*vperp2(ig,imu)*spec(is)%tz*bpar(ig,it,ik)*fbpar
            end do
         end do
      end do

    end subroutine load_kx_bpar

    subroutine load_ky_bpar

      use dist_fn_arrays, only: aj1
      use gs2_layouts, only: is_idx, ik_idx, imu_idx
      use vpamu_grids, only: vperp2
      use kt_grids, only: aky

      integer :: imu, it, ik, is, ig
      complex :: fac

! Is this factor of two from the old normalization?

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         ik = ik_idx(g_lo,iglo)
         is = is_idx(g_lo,iglo)
         do it = 1, ntheta0
            do ig = -ntgrid, ntgrid
               g1(ig,:,it,iglo) = g1(ig,:,it,iglo) + zi*aky(ik)*aj1(ig,it,iglo) &
                    *2.0*vperp2(ig,imu)*spec(is)%tz*bpar(ig,it,ik)*fbpar
            end do
         end do
      end do

    end subroutine load_ky_bpar

#ifdef LOWFLOW
    subroutine load_dthet_phi (phi, g1)

      use dist_fn_arrays, only: aj0
      use theta_grid, only: delthet, ntgrid
      use gs2_layouts, only: g_lo, ik_idx
      use kt_grids, only: ntheta0
      use run_parameters, only: fphi, rhostar
      use vpamu_grids, only: nvgrid

      implicit none

      integer :: iglo, ik, it, ig
      complex :: fac

      complex, dimension (-ntgrid:,:,:), intent (in) :: phi
      complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (out) :: g1

      ! obtain d<phi>/dtheta

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         ik = ik_idx(g_lo,iglo)
         do it = 1, ntheta0
            do ig = -ntgrid, ntgrid-1
               fac = (aj0(ig+1,it,iglo)*phi(ig+1,it,ik)-aj0(ig,it,iglo)*phi(ig,it,ik))*fphi/delthet(ig)
               g1(ig,:,it,iglo) = fac*0.5*rhostar
            end do
         end do
      end do

    end subroutine load_dthet_phi

    subroutine load_kxky_phi (phi, g1)

      use dist_fn_arrays, only: aj0
      use theta_grid, only: ntgrid, gds23, gds24_noq, thet_imp
      use gs2_layouts, only: g_lo, ik_idx
      use kt_grids, only: ntheta0, akx, aky
      use run_parameters, only: fphi, rhostar
      use vpamu_grids, only: nvgrid, vpa_imp
      use constants, only: zi
      use centering, only: get_cell_value

      implicit none

      integer :: iglo, ik, it, ig
      complex :: fac

      complex, dimension (-ntgrid:,:,:), intent (in) :: phi
      complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (out) :: g1

      ! obtain v_E . grad theta in k-space at cell centers

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         ik = ik_idx(g_lo,iglo)
         do it = 1, ntheta0
            do ig = -ntgrid, ntgrid-1
               fac = zi*(akx(it)*gds24_noq(ig)-aky(ik)*gds23(ig))*aj0(ig,it,iglo)*phi(ig,it,ik)*fphi
               g1(ig,:,it,iglo) = fac*rhostar*0.5
            end do
            ! convert from grid values to cell values
            call get_cell_value (thet_imp, vpa_imp, g1(:,:,it,iglo), g1(:,:,it,iglo), &
                 -ntgrid, -nvgrid)
         end do
      end do

    end subroutine load_kxky_phi
#endif

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

  subroutine write_mpdistyxf (dist, extension)

    use mp, only: proc0, send, receive
    use file_utils, only: open_output_file, close_output_file
    use gs2_layouts, only: yxf_lo, ig_idx, iv_idx, is_idx, it_idx
    use gs2_layouts, only: imu_idx, idx_local, proc_id
    use gs2_time, only: code_time
    use theta_grid, only: ntgrid, bmag, theta, shat
    use vpamu_grids, only: vpa, nvgrid, mu
    use kt_grids, only: theta0, lx, ly
    use constants, only: pi

    implicit none
    
    real, dimension (:,yxf_lo%llim_proc:), intent (in) :: dist
    character (*), intent (in) :: extension
    
    integer :: i, j, it, is, imu, ig, iv
    integer, save :: unit
    logical, save :: done = .false.
    real :: gtmp
    
    if (.not. done) then
       if (proc0) call open_output_file (unit, trim(extension))
       do j=yxf_lo%llim_world, yxf_lo%ulim_world
          ig = ig_idx(yxf_lo, j)
          iv = iv_idx(yxf_lo, j)
          it = it_idx(yxf_lo, j)
          is = is_idx(yxf_lo, j) ; if (is /= 1) cycle
          imu = imu_idx(yxf_lo, j)
          do i = 1, yxf_lo%ny
             if (idx_local (yxf_lo, ig, iv, it, imu, is)) then
                if (proc0) then
                   gtmp = dist(i,j)
                else
                   call send (dist(i,j), 0)
                end if
             else if (proc0) then
                call receive (gtmp, proc_id(yxf_lo, j))
             end if
             if (proc0) then
                write (unit,'(a1,8e14.4)') "", code_time, theta(ig), vpa(iv), mu(imu), bmag(ig), &
                     (lx/yxf_lo%nx)*(it-1), (ly/yxf_lo%ny)*(i-1), gtmp
             end if
          end do
          if (proc0) then
             write (unit,*)
             write (unit,*)
          end if
       end do
       if (proc0) call close_output_file (unit)
    end if
    
  end subroutine write_mpdistyxf

  ! subroutine used for testing
  ! takes as input an array using g_lo and
  ! writes it to a .distmp output file
  subroutine write_mpdist (dist, extension)

    use mp, only: proc0, send, receive
    use file_utils, only: open_output_file, close_output_file
    use gs2_layouts, only: g_lo, ik_idx, is_idx
    use gs2_layouts, only: imu_idx, idx_local, proc_id
    use gs2_time, only: code_time
    use theta_grid, only: ntgrid, bmag, theta
    use vpamu_grids, only: vpa, nvgrid, mu
    use kt_grids, only: theta0, ntheta0, akx, aky

    implicit none
    
    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in) :: dist
    character (*), intent (in) :: extension
    
    integer :: iglo, ik, it, is, imu, ig, iv
    integer, save :: unit
    complex :: gtmp
    
    if (proc0) call open_output_file (unit, trim(extension))
    do iglo=g_lo%llim_world, g_lo%ulim_world
       ik = ik_idx(g_lo, iglo)
       is = is_idx(g_lo, iglo)
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
                   write (unit,'(a1,11e14.4,6i4)') "", code_time, theta(ig), vpa(iv), mu(imu), bmag(ig), &
                        real(gtmp), aimag(gtmp), theta(ig)-theta0(it,ik), theta0(it,ik), aky(ik), akx(it), ig, imu, iv, ik, it, is
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
    
  end subroutine write_mpdist

end module nonlinear_terms


