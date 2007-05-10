module additional_linear_terms
  implicit none

  public :: init_additional_linear_terms
  public :: add_additional_linear_terms
  public :: reset_init

  private

  logical :: use_shmem

  logical :: any_additional_term

  logical :: phi0_term
  ! phi0_term variables
  real :: fexp_phi0
  complex, dimension (:,:), allocatable :: fac_phi0
  ! nx, naky

  logical :: wstar_term
  ! wstar_term variables
  real :: fexp_wstar
  real, dimension (:,:,:,:), allocatable :: dwstar
  ! nx, naky, negrid, nspec

  logical :: initialized = .false.
  logical :: accelerated_x = .false.
  logical :: accelerated_v = .false.

contains

  subroutine init_additional_linear_terms
    use theta_grid, only: init_theta_grid
    use gs2_layouts, only: init_gs2_layouts
    implicit none

    if (initialized) return
    initialized = .true.

    call init_gs2_layouts
    call init_theta_grid
    call read_parameters

    if (phi0_term) call init_phi0_term
    if (wstar_term) call init_wstar_term

    any_additional_term = phi0_term .or. wstar_term
  end subroutine init_additional_linear_terms

  subroutine read_parameters
    use file_utils, only: input_unit, input_unit_exist
    use mp, only: proc0, broadcast
    implicit none
    integer in_file
    logical exist
    namelist /additional_linear_terms_knobs/ phi0_term, wstar_term, use_shmem

    if (proc0) then
       phi0_term = .false.
       wstar_term = .false.
       use_shmem = .true.
       in_file = input_unit_exist("additional_linear_terms_knobs",exist)
       if(exist) read (unit=in_file, nml=additional_linear_terms_knobs)
    end if

    call broadcast (phi0_term)
    call broadcast (wstar_term)
    call broadcast (use_shmem)
  end subroutine read_parameters

  subroutine add_additional_linear_terms &
       (g, g1, phi, apar, aperp, phinew, aparnew, aperpnew)
    use theta_grid, only: ntgrid
    use gs2_transforms, only: transform_x, inverse_x
    use gs2_layouts, only: g_lo, xxf_lo, ik_idx, it_idx, ie_idx, is_idx
    use prof, only: prof_entering, prof_leaving
    use dist_fn_arrays, only: vpa, vperp2, aj0, aj1
    use species, only: spec
    use run_parameters, only: fphi, fapar, faperp, wunits
    use constants
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    aperp
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, aperpnew

    complex, dimension (:,:), allocatable:: xxf, phixxf
    complex, dimension (-ntgrid:ntgrid) :: phigavg, apargavg
    integer :: iglo, ixxf, ik, it, ie, is

    if (.not. any_additional_term) return

    call prof_entering ("add_additional_linear_terms", &
         "additional_linear_terms")

    allocate (xxf(xxf_lo%nx,xxf_lo%llim_proc:xxf_lo%ulim_alloc))

    call transform_x (g, xxf)

    if (phi0_term) then
       do ixxf = xxf_lo%llim_proc, xxf_lo%ulim_proc
          ik = ik_idx(xxf_lo,ixxf)
          xxf(:,ixxf) = xxf(:,ixxf)*fac_phi0(:,ik)
       end do
    end if

    if (wstar_term) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          phigavg = (fexp_wstar*phi(:,it,ik) &
                     + (1.0-fexp_wstar)*phinew(:,it,ik)) &
                    *aj0(:,iglo)*fphi &
                  + (fexp_wstar*aperp(:,it,ik) &
                     + (1.0-fexp_wstar)*aperpnew(:,it,ik)) &
                    *aj1(:,iglo)*faperp*2.0*vperp2(:,iglo)*spec(is)%tz
          apargavg = (fexp_wstar*apar(:,it,ik) &
                      + (1.0-fexp_wstar)*aparnew(:,it,ik)) &
                     *aj0(:,iglo)*fapar
          g1(:,1,iglo) &
               = phigavg - spec(is)%stm*vpa(:,1,iglo)*wunits(ik)*apargavg
          g1(:,2,iglo) &
               = phigavg - spec(is)%stm*vpa(:,2,iglo)*wunits(ik)*apargavg
       end do
       allocate (phixxf(xxf_lo%nx,xxf_lo%llim_proc:xxf_lo%ulim_alloc))
       call transform_x (g1, phixxf)

       do ixxf = xxf_lo%llim_proc, xxf_lo%ulim_proc
          ik = ik_idx(xxf_lo,ixxf)
          ie = ie_idx(xxf_lo,ixxf)
          is = is_idx(xxf_lo,ixxf)
          xxf(:,ixxf) = xxf(:,ixxf) + zi*phixxf(:,ixxf)*dwstar(:,ik,ie,is)
       end do
       deallocate (phixxf)
    end if

    call inverse_x (xxf, g)

    deallocate (xxf)

    call prof_leaving ("add_additional_linear_terms", &
         "additional_linear_terms")
  end subroutine add_additional_linear_terms

  subroutine init_phi0_term
    use theta_grid, only: init_theta_grid, ntgrid
    use kt_grids, only: init_kt_grids, naky, ntheta0, aky, akx, nx
    use le_grids, only: init_le_grids, nlambda, negrid
    use species, only: init_species, nspec
    use gs2_layouts, only: init_x_transform_layouts 
    use gs2_transforms, only: init_x_transform
    use run_parameters, only: init_run_parameters, rhostar, delt, tunits
    use file_utils, only: input_unit, error_unit
    use mp, only: proc0, broadcast
    use text_options
    use constants
    implicit none
    integer, parameter :: phi0_option_none = 1, phi0_option_constant_ve = 2, &
         phi0_option_linear_ve = 3, phi0_option_sine = 4
    type (text_option), dimension (5), parameter :: phi0opts = &
         (/ text_option('default', phi0_option_none), &
            text_option('none', phi0_option_none), &
            text_option('constant-ve', phi0_option_constant_ve), &
            text_option('linear-ve', phi0_option_linear_ve), &
            text_option('sine', phi0_option_sine) /)
    character(20) :: phi0_option
    integer :: phi0_option_switch
    real :: phi0
    namelist /phi0_term_knobs/ phi0, phi0_option, fexp_phi0

    real, dimension (:), allocatable :: x, dphi0
    real :: lx
    integer :: i

    call init_run_parameters
    call init_theta_grid
    call init_kt_grids
    call init_le_grids (accelerated_x, accelerated_v)
    call init_species
    call init_x_transform_layouts &
         (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)
    call init_x_transform (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)

    if (proc0) then
       phi0_option = 'default'
       phi0 = 1.0
       fexp_phi0 = 0.5
       read (unit=input_unit("phi0_term_knobs"), nml=phi0_term_knobs)

       i = error_unit()
       call get_option_value &
            (phi0_option, phi0opts, phi0_option_switch, &
             i, "phi0_option in phi0_term_knobs")
    end if
    call broadcast (phi0_option_switch)
    call broadcast (phi0)
    call broadcast (fexp_phi0)

    lx = 2.0*pi/akx(2)

    allocate (x(nx), dphi0(nx))
    x = (/ (lx*real(i)/real(nx), i=-nx/2+1,-nx/2+nx) /)

    select case (phi0_option_switch)
    case (phi0_option_none)
       dphi0 = 0.0
    case (phi0_option_constant_ve)
       dphi0 = phi0/rhostar
    case (phi0_option_linear_ve)
       dphi0 = phi0*x/rhostar
    case (phi0_option_sine)
       dphi0 = phi0*2.0*pi/lx*cos(2.0*pi*x/lx)/rhostar
    end select

    if(.not.allocated(fac_phi0)) allocate (fac_phi0(nx,naky))
    do i = 1, naky
       fac_phi0(:,i) = cmplx(1.0,-aky(i)*tunits(i)*delt*fexp_phi0*dphi0) &
                       /cmplx(1.0,aky(i)*tunits(i)*delt*(1.0-fexp_phi0)*dphi0)
    end do
    deallocate (x, dphi0)
  end subroutine init_phi0_term

  subroutine init_wstar_term
    use file_utils, only: input_unit, error_unit
    use text_options
    use theta_grid, only: init_theta_grid, ntgrid
    use kt_grids, only: init_kt_grids, naky, ntheta0, akx, nx
    use le_grids, only: init_le_grids, nlambda, negrid, e
    use species, only: init_species, nspec
    use gs2_layouts, only: init_x_transform_layouts 
    use gs2_transforms, only: init_x_transform
    use run_parameters, only: init_run_parameters, delt, wunits
    use mp, only: proc0, broadcast
    use constants
    implicit none
    integer, parameter :: wstar_option_none = 1, wstar_option_constant = 2, &
         wstar_option_localized = 3
    type (text_option), dimension (4), parameter :: wstaropts = &
         (/ text_option("default", wstar_option_none), &
            text_option("none", wstar_option_none), &
            text_option("constant", wstar_option_constant), &
            text_option("localized", wstar_option_localized) /)
    character(20) :: wstar_option
    real :: dfprim0, dtprim0, dtp_width, dfp_width
    namelist /wstar_term_knobs/ wstar_option, fexp_wstar, dfprim0, dtprim0, &
         dtp_width, dfp_width

    integer :: wstar_option_switch
    integer :: i, ik, ie, is
    real, dimension (:), allocatable :: x
    real :: lx

    call init_run_parameters
    call init_theta_grid
    call init_kt_grids
    call init_le_grids (accelerated_x, accelerated_v)
    call init_species
    call init_x_transform_layouts &
         (ntgrid, naky, ntheta0, nlambda, negrid, nx, nspec)
    call init_x_transform (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)

    if (proc0) then
       wstar_option = 'default'
       dfprim0 = 0.0
       dtprim0 = 0.0
       dfp_width = 1e10
       dtp_width = 1e10
       fexp_wstar = 0.4
       read (unit=input_unit("wstar_term_knobs"), nml=wstar_term_knobs)

       i = error_unit()
       call get_option_value &
            (wstar_option, wstaropts, wstar_option_switch, &
             i, "wstar_option in wstar_term_knobs")
    end if
    call broadcast (wstar_option_switch)
    call broadcast (dfprim0)
    call broadcast (dtprim0)
    call broadcast (dfp_width)
    call broadcast (dtp_width)
    call broadcast (fexp_wstar)

    lx = 2.0*pi/akx(2)
    allocate (x(nx))
    x = (/ (lx*real(i)/real(nx), i=-nx/2+1,-nx/2+nx) /)

    if(.not.allocated(dwstar)) allocate (dwstar(nx,naky,negrid,nspec))

    select case (wstar_option_switch)
    case (wstar_option_none)
       dwstar = 0.0
    case (wstar_option_constant)
       do is = 1, nspec
          do ie = 1, negrid
             dwstar(:,:,ie,is) = dfprim0 + dtprim0*(e(ie,is)-1.5)
          end do
       end do
    case (wstar_option_localized)
       do is = 1, nspec
          do ie = 1, negrid
             do ik = 1, naky
                dwstar(:,ik,ie,is) &
                     = dfprim0*exp(-(x/dfp_width)**2) &
                       + dtprim0*exp(-(x/dtp_width)**2)*(e(ie,is)-1.5)
             end do
          end do
       end do
    end select

    do ik = 1, naky
       dwstar(:,ik,:,:) = dwstar(:,ik,:,:)*delt*wunits(ik)
    end do

    deallocate (x)
  end subroutine init_wstar_term

  subroutine reset_init

    initialized = .false.

  end subroutine reset_init

end module additional_linear_terms
