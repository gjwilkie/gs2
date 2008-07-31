module init_g
  implicit none

  public :: ginit
  public :: init_init_g
  public :: width0, k0
  public :: tstart
  public :: reset_init
  public :: init_vnmult
  private

  ! knobs
  integer :: ginitopt_switch
  integer, parameter :: ginitopt_default = 1, ginitopt_test1 = 2, &
       ginitopt_xi = 3, ginitopt_xi2 = 4, ginitopt_rh = 5, ginitopt_zero = 6, &
       ginitopt_test3 = 7, ginitopt_convect = 8, ginitopt_restart_file = 9, &
       ginitopt_noise = 10, ginitopt_restart_many = 11, ginitopt_continue = 12, &
       ginitopt_nl = 13, ginitopt_kz0 = 14, ginitopt_restart_small = 15, &
       ginitopt_nl2 = 16, ginitopt_nl3 = 17, ginitopt_nl4 = 18, & 
       ginitopt_nl5 = 19, ginitopt_alf = 20, ginitopt_kpar = 21, &
       ginitopt_nl6 = 22, ginitopt_nl7 = 23, ginitopt_gs = 24, ginitopt_recon = 25, &
       ginitopt_nl3r = 26, ginitopt_smallflat = 27, ginitopt_harris = 28
  real :: width0, dphiinit, phiinit, k0, imfac, refac, zf_init
  real :: den0, upar0, tpar0, tperp0
  real :: den1, upar1, tpar1, tperp1
  real :: den2, upar2, tpar2, tperp2
  real :: tstart, scale, apar0
  logical :: chop_side, left, even
  character(300) :: restart_file
  integer, dimension(2) :: ikk, itt
  
contains

  subroutine init_init_g
    use gs2_save, only: init_save
    use gs2_layouts, only: init_gs2_layouts
    use mp, only: proc0, broadcast
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_gs2_layouts

    if (proc0) call read_parameters

    call broadcast (ginitopt_switch)
    call broadcast (width0)
    call broadcast (refac)
    call broadcast (imfac)
    call broadcast (den0)
    call broadcast (upar0)
    call broadcast (tpar0)
    call broadcast (tperp0)
    call broadcast (den1)
    call broadcast (upar1)
    call broadcast (tpar1)
    call broadcast (tperp1)
    call broadcast (den2)
    call broadcast (upar2)
    call broadcast (tpar2)
    call broadcast (tperp2)
    call broadcast (phiinit)
    call broadcast (dphiinit)
    call broadcast (zf_init)
    call broadcast (k0)
    call broadcast (apar0)
    call broadcast (tstart)
    call broadcast (chop_side)
    call broadcast (even)
    call broadcast (left)
    call broadcast (restart_file)
    call broadcast (ikk)
    call broadcast (itt) 
    call broadcast (scale)

    call init_save (restart_file)

  end subroutine init_init_g

  subroutine ginit (restarted)

    use gs2_save, only: init_tstart
    logical, intent (out) :: restarted
    real :: t0
    integer :: istatus

    restarted = .false.
    select case (ginitopt_switch)
    case (ginitopt_default)
       call ginit_default
    case (ginitopt_kz0)
       call ginit_kz0
    case (ginitopt_noise)
       call ginit_noise
    case (ginitopt_kpar)
       call ginit_kpar
    case (ginitopt_gs)
       call ginit_gs
    case (ginitopt_nl)
       call ginit_nl
    case (ginitopt_nl2)
       call ginit_nl2
    case (ginitopt_nl3)
       call ginit_nl3
    case (ginitopt_nl3r)
       call ginit_nl3r
    case (ginitopt_harris)
       call ginit_harris
    case (ginitopt_nl4)
! in an old version, this line was commented out.  Thus, to recover some old
! results, you might need to comment out this line and...
!       t0 = tstart
       call ginit_nl4
       call init_tstart (tstart, istatus)
! this line:
!       tstart = t0
       restarted = .true.
       scale = 1.
    case (ginitopt_nl5)
       t0 = tstart
       call init_tstart (tstart, istatus)
       call ginit_nl5
       tstart = t0
       restarted = .true.
       scale = 1.
    case (ginitopt_recon)
       t0 = tstart
       call init_tstart (tstart, istatus)
       call ginit_recon
       tstart = t0
       restarted = .true.
       scale = 1.
    case (ginitopt_nl6)
       t0 = tstart
       call init_tstart (tstart, istatus)
       call ginit_nl6
       tstart = t0
       restarted = .true.
       scale = 1.
    case (ginitopt_nl7)
       call ginit_nl7
    case (ginitopt_test1)
       call ginit_test1
    case (ginitopt_xi)
       call ginit_xi
    case (ginitopt_xi2)
       call ginit_xi2
    case (ginitopt_rh)
       call ginit_rh
    case (ginitopt_alf)
       call ginit_alf
    case (ginitopt_zero)
       call ginit_zero
    case (ginitopt_test3)
       call ginit_test3
    case (ginitopt_convect)
       call ginit_convect
    case (ginitopt_restart_file)
       call ginit_restart_file 
       call init_tstart (tstart, istatus)
       restarted = .true.
       scale = 1.
    case (ginitopt_restart_many)
       call ginit_restart_many 
       call init_tstart (tstart, istatus)
       restarted = .true.
       scale = 1.
    case (ginitopt_restart_small)
       call ginit_restart_small
       call init_tstart (tstart, istatus)
       restarted = .true.
       scale = 1.
    case (ginitopt_smallflat)
       call ginit_restart_smallflat
       call init_tstart (tstart, istatus)
       restarted = .true.
       scale = 1.
    case (ginitopt_continue)
       restarted = .true.
       scale = 1.
    end select
  end subroutine ginit

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, run_name, input_unit_exist
    use text_options, only: text_option, get_option_value
    implicit none

    type (text_option), dimension (28), parameter :: ginitopts = &
         (/ text_option('default', ginitopt_default), &
            text_option('noise', ginitopt_noise), &
            text_option('test1', ginitopt_test1), &
            text_option('xi', ginitopt_xi), &
            text_option('xi2', ginitopt_xi2), &
            text_option('zero', ginitopt_zero), &
            text_option('test3', ginitopt_test3), &
            text_option('convect', ginitopt_convect), &
            text_option('rh', ginitopt_rh), &
            text_option('many', ginitopt_restart_many), &
            text_option('small', ginitopt_restart_small), &
            text_option('file', ginitopt_restart_file), &
            text_option('cont', ginitopt_continue), &
            text_option('kz0', ginitopt_kz0), &
            text_option('nl', ginitopt_nl), &
            text_option('nl2', ginitopt_nl2), &
            text_option('nl3', ginitopt_nl3), &
            text_option('nl3r', ginitopt_nl3r), &
            text_option('nl4', ginitopt_nl4), &
            text_option('nl5', ginitopt_nl5), &
            text_option('nl6', ginitopt_nl6), &
            text_option('nl7', ginitopt_nl7), &
            text_option('alf', ginitopt_alf), &
            text_option('gs', ginitopt_gs), &
            text_option('kpar', ginitopt_kpar), &
            text_option('smallflat', ginitopt_smallflat), &
            text_option('harris', ginitopt_harris), &
            text_option('recon', ginitopt_recon) /)
    character(20) :: ginit_option
    namelist /init_g_knobs/ ginit_option, width0, phiinit, k0, chop_side, &
         restart_file, left, ikk, itt, scale, tstart, zf_init, &
         den0, upar0, tpar0, tperp0, imfac, refac, even, &
         den1, upar1, tpar1, tperp1, &
         den2, upar2, tpar2, tperp2, dphiinit, apar0

    integer :: ierr, in_file
    logical :: exist

    tstart = 0.
    scale = 1.0
    ginit_option = "default"
    width0 = -3.5
    refac = 1.
    imfac = 0.
    den0 = 1.
    upar0 = 0.
    tpar0 = 0.
    tperp0 = 0.
    den1 = 0.
    upar1 = 0.
    tpar1 = 0.
    tperp1 = 0.
    den2 = 0.
    upar2 = 0.
    tpar2 = 0.
    tperp2 = 0.
    phiinit = 1.0
    dphiinit = 1.0
    zf_init = 1.0
    k0 = 1.0
    apar0 = 0.
    chop_side = .true.
    left = .true.
    even = .true.
    ikk(1) = 1
    ikk(2) = 2
    itt(1) = 1
    itt(2) = 2
    restart_file = trim(run_name)//".nc"
    in_file = input_unit_exist ("init_g_knobs", exist)
    if (exist) read (unit=input_unit("init_g_knobs"), nml=init_g_knobs)

    ierr = error_unit()
    call get_option_value &
         (ginit_option, ginitopts, ginitopt_switch, &
         ierr, "ginit_option in ginit_knobs")
  end subroutine read_parameters

  subroutine ginit_default
    use species, only: spec
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, theta0, aky, reality
    use le_grids, only: forbid
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    implicit none
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phi
    logical :: right
    integer :: iglo
    integer :: ig, ik, it, il, is

    right = .not. left

    do ig = -ntgrid, ntgrid
       phi(ig,:,:) = exp(-((theta(ig)-theta0(:,:))/width0)**2)*cmplx(1.0,1.0)            
    end do
    if (chop_side .and. left) phi(:-1,:,:) = 0.0
    if (chop_side .and. right) phi(1:,:,:) = 0.0
    
    if (reality) then
       phi(:,1,1) = 0.0

       if (naky > 1 .and. aky(1) == 0.0) then
          phi(:,:,1) = 0.0
       end if

! not used:
! reality condition for k_theta = 0 component:
       do it = 1, ntheta0/2
          phi(:,it+(ntheta0+1)/2,1) = conjg(phi(:,(ntheta0+1)/2+1-it,1))
       enddo
    end if

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       g(:,1,iglo) = phi(:,it,ik)*spec(is)%z*phiinit
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       g(:,2,iglo) = g(:,1,iglo)
    end do
    gnew = g
  end subroutine ginit_default

  subroutine ginit_kz0
    use species, only: spec
    use theta_grid, only: ntgrid 
    use kt_grids, only: naky, ntheta0, aky, reality
    use le_grids, only: forbid
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    implicit none
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phi
    logical :: right
    integer :: iglo
    integer :: ik, it, il, is

    right = .not. left

    phi = cmplx(1.0,1.0)
    if (chop_side .and. left) phi(:-1,:,:) = 0.0
    if (chop_side .and. right) phi(1:,:,:) = 0.0
    
    if (reality) then
       phi(:,1,1) = 0.0

       if (naky > 1 .and. aky(1) == 0.0) then
          phi(:,:,1) = 0.0
       end if

! not used:
! reality condition for k_theta = 0 component:
       do it = 1, ntheta0/2
          phi(:,it+(ntheta0+1)/2,1) = conjg(phi(:,(ntheta0+1)/2+1-it,1))
       enddo
    end if

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       g(:,1,iglo) = -phi(:,it,ik)*spec(is)%z*phiinit
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       g(:,2,iglo) = g(:,1,iglo)
    end do
    gnew = g
  end subroutine ginit_kz0

!  subroutine ginit_noise
!    use species, only: spec
!    use theta_grid, only: ntgrid 
!    use kt_grids, only: naky, ntheta0, aky, reality
!    use le_grids, only: forbid
!    use dist_fn_arrays, only: g, gnew
!    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
!    use ran
!    implicit none
!    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phi
!    real :: a, b
!    integer :: iglo
!    integer :: ig, ik, it, il, is
!
!! keep old (it, ik) loop order to get old results exactly: 
!    do it = 1, ntheta0
!       do ik = 1, naky
!          do ig = -ntgrid, ntgrid
!             a = ranf()-0.5
!             b = ranf()-0.5
!!             phi(:,it,ik) = cmplx(a,b)
!             phi(ig,it,ik) = cmplx(a,b)
!          end do
!          if (chop_side) then
!             if (left) then
!                phi(:-1,it,ik) = 0.0
!             else
!                phi(1:,it,ik) = 0.0
!             end if
!          end if
!       end do
!    end do
!
!    if (naky > 1 .and. aky(1) == 0.0) then
!       phi(:,:,1) = phi(:,:,1)*zf_init
!    end if
!! reality condition for k_theta = 0 component:
!    if (reality) then
!       do it = 1, ntheta0/2
!          phi(:,it+(ntheta0+1)/2,1) = conjg(phi(:,(ntheta0+1)/2+1-it,1))
!       enddo
!    end if
!       
!
!    do iglo = g_lo%llim_proc, g_lo%ulim_proc
!       ik = ik_idx(g_lo,iglo)
!       it = it_idx(g_lo,iglo)
!       il = il_idx(g_lo,iglo)
!       is = is_idx(g_lo,iglo)
!       g(:,1,iglo) = -phi(:,it,ik)*spec(is)%z*phiinit
!       where (forbid(:,il)) g(:,1,iglo) = 0.0
!       g(:,2,iglo) = g(:,1,iglo)
!    end do
!    gnew = g
!  end subroutine ginit_noise
  
  subroutine ginit_noise
    use species, only: spec, tracer_species
    use theta_grid, only: ntgrid 
    use kt_grids, only: naky, ntheta0, aky, reality
    use le_grids, only: forbid
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use ran
    implicit none
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phi, phit
    real :: a, b
    integer :: iglo
    integer :: ig, ik, it, il, is, nn

    phit = 0.
    do it=2,ntheta0/2+1
       nn = it-1
! extra factor of 4 to reduce the high k wiggles for now
       phit (:, it, 1) = (-1)**nn*exp(-8.*(real(nn)/ntheta0)**2)
    end do

! keep old (it, ik) loop order to get old results exactly: 
    do it = 1, ntheta0
       do ik = 1, naky
          do ig = -ntgrid, ntgrid
             a = ranf()-0.5
             b = ranf()-0.5
!             phi(:,it,ik) = cmplx(a,b)
             phi(ig,it,ik) = cmplx(a,b)
           end do
          if (chop_side) then
             if (left) then
                phi(:-1,it,ik) = 0.0
             else
                phi(1:,it,ik) = 0.0
             end if
          end if
       end do
    end do

    if (naky > 1 .and. aky(1) == 0.0) then
       phi(:,:,1) = phi(:,:,1)*zf_init
    end if
! reality condition for k_theta = 0 component:
    if (reality) then
       do it = 1, ntheta0/2
          phi(:,it+(ntheta0+1)/2,1) = conjg(phi(:,(ntheta0+1)/2+1-it,1))
          phit(:,it+(ntheta0+1)/2,1) = conjg(phit(:,(ntheta0+1)/2+1-it,1))
       enddo
    end if
       
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       if (spec(is)%type == tracer_species) then          
          g(:,1,iglo) =-phit(:,it,ik)*spec(is)%z*phiinit
       else
          g(:,1,iglo) = -phi(:,it,ik)*spec(is)%z*phiinit
       end if
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       g(:,2,iglo) = g(:,1,iglo)
    end do
    gnew = g

  end subroutine ginit_noise
  
  subroutine ginit_nl
    use species, only: spec
    use theta_grid, only: ntgrid 
    use kt_grids, only: naky, ntheta0 
    use le_grids, only: forbid, ng2
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use ran
    implicit none
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phi
    integer :: iglo
    integer :: ig, ik, it, il, is, j
    
    phi = 0.
    do j = 1, 2
       ik = ikk(j)
       it = itt(j)
       do ig = -ntgrid, ntgrid
!          phi(ig,it,ik) = cmplx(1.0, 0.0)*sin(theta(ig))
          phi(ig,it,ik) = cmplx(1.0, 0.0)!*sin(theta(ig))
       end do
!       do ig = -ntgrid/2, ntgrid/2
!          phi(ig,it,ik) = cmplx(1.0, 1.0)
!       end do
       if (chop_side) then
          if (left) then
             phi(:-1,it,ik) = 0.0
          else
             phi(1:,it,ik) = 0.0
          end if
       end if
    end do
    
! reality condition for k_theta = 0 component:
    do it = 1, ntheta0/2
       phi(:,it+(ntheta0+1)/2,1) = conjg(phi(:,(ntheta0+1)/2+1-it,1))
    enddo

! charge dependence keeps initial Phi from being too small
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)       
       g(:,1,iglo) = -phi(:,it,ik)*phiinit*spec(is)%z
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       g(:,2,iglo) = g(:,1,iglo)
       if (il == ng2+1) g(:,:,iglo) = 0.0
    end do
    gnew = g
  end subroutine ginit_nl

  subroutine ginit_nl2

    use species, only: spec
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0 
    use le_grids, only: forbid, ng2
    use dist_fn_arrays, only: g, gnew, vpa
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use ran
    implicit none
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phi
    integer :: iglo
    integer :: ig, ik, it, il, is, j
    
    phi = 0.
    do j = 1, 2
       ik = ikk(j)
       it = itt(j)
       do ig = -ntgrid, ntgrid
          phi(ig,it,ik) = cmplx(1.0, 0.0)!*sin(theta(ig))
       end do
!       do ig = -ntgrid/2, ntgrid/2
!          phi(ig,it,ik) = cmplx(1.0, 1.0)
!       end do

       if (chop_side) then
          if (left) then
             phi(:-1,it,ik) = 0.0
          else
             phi(1:,it,ik) = 0.0
          end if
       end if

    end do
    
! reality condition for k_theta = 0 component:
    do it = 1, ntheta0/2
       phi(:,it+(ntheta0+1)/2,1) = conjg(phi(:,(ntheta0+1)/2+1-it,1))
    enddo
    
! charge dependence keeps initial Phi from being too small
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)       
       g(:,1,iglo) = -phi(:,it,ik)*phiinit*(1.+vpa(:,1,iglo)*sin(theta))*spec(is)%z
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       g(:,2,iglo) = -phi(:,it,ik)*phiinit*(1.+vpa(:,2,iglo)*sin(theta))*spec(is)%z
!       g(:,1,iglo)*vpa(:,2,iglo)
       if (il == ng2+1) g(:,:,iglo) = 0.0
    end do

    gnew = g
  end subroutine ginit_nl2

  subroutine ginit_nl3
    use species, only: spec, has_electron_species, electron_species
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, theta0, reality
    use le_grids, only: forbid
!    use le_grids, only: ng2
    use dist_fn_arrays, only: g, gnew, vpa, vperp2
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use constants
    use ran
    implicit none
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phi, odd
    real, dimension (-ntgrid:ntgrid) :: dfac, ufac, tparfac, tperpfac, ct, st, c2t, s2t
    integer :: iglo
    integer :: ig, ik, it, il, is, j
    
    phi = 0.
    odd = 0.
    do j = 1, 2
       ik = ikk(j)
       it = itt(j)
       if (width0 > 0.) then
          do ig = -ntgrid, ntgrid
             phi(ig,it,ik) = exp(-((theta(ig)-theta0(it,ik))/width0)**2)*cmplx(refac, imfac)
          end do
       else
          do ig = -ntgrid, ntgrid
             phi(ig,it,ik) = cmplx(refac, imfac)
          end do
       end if
       if (chop_side) then
          if (left) then
             phi(:-1,it,ik) = 0.0
          else
             phi(1:,it,ik) = 0.0
          end if
       end if
    end do

    odd = zi * phi
    
! reality condition for k_theta = 0 component:
    if (reality) then
       do it = 1, ntheta0/2
          phi(:,it+(ntheta0+1)/2,1) = conjg(phi(:,(ntheta0+1)/2+1-it,1))
          odd(:,it+(ntheta0+1)/2,1) = conjg(odd(:,(ntheta0+1)/2+1-it,1))
       enddo
    end if

    if (even) then
       ct = cos(theta)
       st = sin(theta)

       c2t = cos(2.*theta)
       s2t = sin(2.*theta)
    else
       ct = sin(theta)
       st = cos(theta)

       c2t = sin(2.*theta)
       s2t = cos(2.*theta)
    end if

    dfac     = den0   + den1 * ct + den2 * c2t
    ufac     = upar0  + upar1* st + upar2* s2t
    tparfac  = tpar0  + tpar1* ct + tpar2* c2t
    tperpfac = tperp0 + tperp1*ct + tperp2*c2t


! charge dependence keeps initial Phi from being too small
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)       
!       if (spec(is)%type /= electron_species) cycle
       g(:,1,iglo) = phiinit* &!spec(is)%z* &
            ( dfac*spec(is)%dens0            * phi(:,it,ik) &
            + 2.*ufac* vpa(:,1,iglo)         * odd(:,it,ik) &
            + tparfac*(vpa(:,1,iglo)**2-0.5) * phi(:,it,ik) &
            +tperpfac*(vperp2(:,iglo)-1.)    * phi(:,it,ik))
       where (forbid(:,il)) g(:,1,iglo) = 0.0

       g(:,2,iglo) = phiinit* &!spec(is)%z* &
            ( dfac*spec(is)%dens0            * phi(:,it,ik) &
            + 2.*ufac* vpa(:,2,iglo)         * odd(:,it,ik) &
            + tparfac*(vpa(:,2,iglo)**2-0.5) * phi(:,it,ik) &
            +tperpfac*(vperp2(:,iglo)-1.)    * phi(:,it,ik))
       where (forbid(:,il)) g(:,2,iglo) = 0.0

!       if (il == ng2+1) g(:,:,iglo) = 0.0
    end do

!    if (has_electron_species(spec)) then
!       call flae (g, gnew)
!       g = g - gnew
!    end if

    gnew = g
  end subroutine ginit_nl3

  subroutine ginit_nl3r
    use species, only: spec, has_electron_species, electron_species
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, theta0, reality
    use le_grids, only: forbid
!    use le_grids, only: ng2
    use fields_arrays, only: apar, bpar
    use fields_arrays, only: aparnew, bparnew
    use dist_fn_arrays, only: g, gnew, vpa, vperp2
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use constants
    use ran
    implicit none
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phi, odd
    real, dimension (-ntgrid:ntgrid) :: dfac, ufac, tparfac, tperpfac, ct, st, c2t, s2t
    integer :: iglo
    integer :: ig, ik, it, il, is, j
    
    phi = 0.
    odd = 0.
    do j = 1, 2
       ik = ikk(j)
       it = itt(j)
       if (width0 > 0.) then
          do ig = -ntgrid, ntgrid
             phi(ig,it,ik) = exp(-((theta(ig)-theta0(it,ik))/width0)**2)*cmplx(refac, imfac)
          end do
       else
          do ig = -ntgrid, ntgrid
             phi(ig,it,ik) = cmplx(refac, imfac)
             apar(ig,it,ik) = apar0*cmplx(refac, imfac)
          end do
       end if
       if (chop_side) then
          if (left) then
             phi(:-1,it,ik) = 0.0
          else
             phi(1:,it,ik) = 0.0
          end if
       end if
    end do

    odd = zi * phi
    
! reality condition for k_theta = 0 component:
    if (reality) then
       do it = 1, ntheta0/2
          apar(:,it+(ntheta0+1)/2,1) = conjg(apar(:,(ntheta0+1)/2+1-it,1))
          phi(:,it+(ntheta0+1)/2,1) = conjg(phi(:,(ntheta0+1)/2+1-it,1))
          odd(:,it+(ntheta0+1)/2,1) = conjg(odd(:,(ntheta0+1)/2+1-it,1))
       enddo
    end if

    aparnew = apar

    if (even) then
       ct = cos(theta)
       st = sin(theta)

       c2t = cos(2.*theta)
       s2t = sin(2.*theta)
    else
       ct = sin(theta)
       st = cos(theta)

       c2t = sin(2.*theta)
       s2t = cos(2.*theta)
    end if

    dfac     = den0   + den1 * ct + den2 * c2t
    ufac     = upar0  + upar1* st + upar2* s2t
    tparfac  = tpar0  + tpar1* ct + tpar2* c2t
    tperpfac = tperp0 + tperp1*ct + tperp2*c2t


! charge dependence keeps initial Phi from being too small
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)       
       g(:,1,iglo) = phiinit* &!spec(is)%z* &
            ( dfac*spec(is)%dens0                * phi(:,it,ik) &
            + 2.*ufac* vpa(:,1,iglo)*spec(is)%u0 * odd(:,it,ik) &
            + tparfac*(vpa(:,1,iglo)**2-0.5)     * phi(:,it,ik) &
            +tperpfac*(vperp2(:,iglo)-1.)        * phi(:,it,ik))
       where (forbid(:,il)) g(:,1,iglo) = 0.0

       g(:,2,iglo) = phiinit* &!spec(is)%z* &
            ( dfac*spec(is)%dens0                * phi(:,it,ik) &
            + 2.*ufac* vpa(:,2,iglo)*spec(is)%u0 * odd(:,it,ik) &
            + tparfac*(vpa(:,2,iglo)**2-0.5)     * phi(:,it,ik) &
            +tperpfac*(vperp2(:,iglo)-1.)        * phi(:,it,ik))
       where (forbid(:,il)) g(:,2,iglo) = 0.0

!       if (il == ng2+1) g(:,:,iglo) = 0.0
    end do

!    if (has_electron_species(spec)) then
!       call flae (g, gnew)
!       g = g - gnew
!    end if

    gnew = g
  end subroutine ginit_nl3r

  subroutine ginit_harris
    use species, only: spec, has_electron_species, electron_species
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, theta0, reality, nx
    use le_grids, only: forbid
    use dist_fn_arrays, only: g, gnew, vpa, vperp2, aj0
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx, il_idx
    use constants
    use ran
    implicit none
    complex, dimension (ntheta0,naky) :: phi
    integer :: iglo
    integer :: ig, ik, it, il, is, j, j1
    real, dimension(nx) :: lx_pr, a_pr
    complex,dimension(nx) :: ff_pr
    real:: L, dx_pr
    
! 
! Specifying function on x grid, with nx points.  Transforming to kx grid, with 
! nkx < nx points.  Result is function that will be used for initial condition.
! But Fourier transform is the same if you just use nkx points in the first 
! place, right?
!    

! Can specify x0 along with y0; no need to use this construction, although it is okay.
    L=2.*pi*k0
    
    
! nx is available from kt_grids:
    dx_pr=L/real(nx)
    
    do j = 1, nx
       lx_pr(j)=dx_pr*real(j-1)
    end do
    
    do j = 1,nx
       a_pr(j)=1./cosh((lx_pr(j)-L/4.)/width0)- &
         1./cosh((lx_pr(j)-3.*L/4.)/width0)
    end do
    
    a_pr=a_pr/real(nx)
        
    do j = 1,nx
       ff_pr(j) = 0.
       do j1 = 1,nx
    
          ff_pr(j)=ff_pr(j)+a_pr(j1)*exp(-zi*2.*pi*real(j-1)*real(j1-1)/real(nx))
    
       end do
    end do
    
    phi = 0.

! Dealiasing here:
    do j = 1, ntheta0/2-1
       phi(j,1) = ff_pr(j)
    end do

!
! Note: if the j=1 coefficient is non-zero, it might be a problem, because this 
! coefficient is never used.
!
    
    do j = ntheta0/2, ntheta0
       phi(j,1) = ff_pr(nx-ntheta0+j)
    end do

!
! presumably this is not needed, but leave it in for now:
!
! reality condition for k_theta = 0 component:
    if (reality) then
       do it = 1, ntheta0/2
          phi(it+(ntheta0+1)/2,1) = conjg(phi((ntheta0+1)/2+1-it,1))
       enddo
    end if
    
    g = 0.
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo) 
       il = il_idx(g_lo,iglo) 
       is = is_idx(g_lo,iglo)       
    
! ions, kx/=0:
       if ((is==1) .and.(.not.(it==1))) then
       
          g(:,1,iglo) =  2.* vpa(:,1,iglo)*spec(is)%u0 * phi(it,ik)/aj0(:,iglo)  
          g(:,2,iglo) =  2.* vpa(:,2,iglo)*spec(is)%u0 * phi(it,ik)/aj0(:,iglo)

          where (forbid(:,il)) g(:,1,iglo) = 0.0
          where (forbid(:,il)) g(:,2,iglo) = 0.0
          
       end if
    end do

    gnew = g

  end subroutine ginit_harris

  subroutine ginit_recon
    use mp, only: proc0
    use species, only: spec, has_electron_species
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, theta0, reality
    use le_grids, only: forbid
!    use le_grids, only: ng2
    use gs2_save, only: gs2_restore
    use dist_fn_arrays, only: g, gnew, vpa, vperp2
    use fields_arrays, only: phi, apar, bpar
    use fields_arrays, only: phinew, aparnew, bparnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use file_utils, only: error_unit
    use run_parameters, only: fphi, fapar, fbpar
    use constants
    use ran
    implicit none
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phiz, odd
    real, dimension (-ntgrid:ntgrid) :: dfac, ufac, tparfac, tperpfac, ct, st, c2t, s2t
    integer :: iglo, istatus, ierr
    integer :: ig, ik, it, il, is, j
    logical :: many = .true.
    
    call gs2_restore (g, scale, istatus, fphi, fapar, fbpar, many)
    if (istatus /= 0) then
       ierr = error_unit()
       if (proc0) write(ierr,*) "Error reading file: ", trim(restart_file)
       g = 0.
    end if

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       if (ik == 1) cycle

       g (:,1,iglo) = 0.
       g (:,2,iglo) = 0.
    end do
	
    phinew(:,:,2:naky) = 0.
    aparnew(:,:,2:naky) = 0.
    bparnew(:,:,2:naky) = 0.
    phi(:,:,2:naky) = 0.
    apar(:,:,2:naky) = 0.
    bpar(:,:,2:naky) = 0.

    phiz = 0.
    odd = 0.
    do j = 1, 2
       ik = ikk(j)
       it = itt(j)
       do ig = -ntgrid, ntgrid
          phiz(ig,it,ik) = cmplx(refac, imfac)
       end do
    end do

    odd = zi * phiz
    
! reality condition for k_theta = 0 component:
    if (reality) then
       do it = 1, ntheta0/2
          phiz(:,it+(ntheta0+1)/2,1) = conjg(phiz(:,(ntheta0+1)/2+1-it,1))
          odd (:,it+(ntheta0+1)/2,1) = conjg(odd (:,(ntheta0+1)/2+1-it,1))
       enddo
    end if

    if (even) then
       ct = cos(theta)
       st = sin(theta)

       c2t = cos(2.*theta)
       s2t = sin(2.*theta)
    else
       ct = sin(theta)
       st = cos(theta)

       c2t = sin(2.*theta)
       s2t = cos(2.*theta)
    end if

    dfac     = den0   + den1 * ct + den2 * c2t
    ufac     = upar0  + upar1* st + upar2* s2t
    tparfac  = tpar0  + tpar1* ct + tpar2* c2t
    tperpfac = tperp0 + tperp1*ct + tperp2*c2t


! charge dependence keeps initial Phi from being too small
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)       

       g(:,1,iglo) = phiinit* &!spec(is)%z* &
            ( dfac*spec(is)%dens0            * phiz(:,it,ik) &
            + 2.*ufac* vpa(:,1,iglo)         * odd (:,it,ik) &
            + tparfac*(vpa(:,1,iglo)**2-0.5) * phiz(:,it,ik) &
            +tperpfac*(vperp2(:,iglo)-1.)    * phiz(:,it,ik))
       where (forbid(:,il)) g(:,1,iglo) = 0.0

       g(:,2,iglo) = phiinit* &!spec(is)%z* &
            ( dfac*spec(is)%dens0            * phiz(:,it,ik) &
            + 2.*ufac* vpa(:,2,iglo)         * odd (:,it,ik) &
            + tparfac*(vpa(:,2,iglo)**2-0.5) * phiz(:,it,ik) &
            +tperpfac*(vperp2(:,iglo)-1.)    * phiz(:,it,ik))
       where (forbid(:,il)) g(:,2,iglo) = 0.0

!       if (il == ng2+1) g(:,:,iglo) = 0.0
    end do

!    if (has_electron_species(spec)) then
!       call flae (g, gnew)
!       g = g - gnew
!    end if

    gnew = g
  end subroutine ginit_recon

! nl4 has been monkeyed with over time.  Do not expect to recover old results
! that use this startup routine.
  subroutine ginit_nl4
    use mp, only: proc0
    use species, only: spec
    use gs2_save, only: gs2_restore
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use dist_fn_arrays, only: g, gnew
    use le_grids, only: forbid
    use fields_arrays, only: phi, apar, bpar
    use fields_arrays, only: phinew, aparnew, bparnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx, il_idx
    use file_utils, only: error_unit
    use run_parameters, only: fphi, fapar, fbpar
    use ran
    implicit none
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phiz
    integer :: iglo, istatus
    integer :: ig, ik, it, is, il, ierr
    logical :: many = .true.
    
    call gs2_restore (g, scale, istatus, fphi, fapar, fbpar, many)
    if (istatus /= 0) then
       ierr = error_unit()
       if (proc0) write(ierr,*) "Error reading file: ", trim(restart_file)
       g = 0.
    end if

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
!       if ((it == 2 .or. it == ntheta0) .and. ik == 1) cycle
!       if (ik == 1) cycle
       if (it == 1 .and. ik == 2) cycle

       g (:,1,iglo) = 0.
       g (:,2,iglo) = 0.
    end do

!    do ik = 1, naky
!       if (ik /= 2) then
!          phinew(:,:,ik) = 0.
!          aparnew(:,:,ik) = 0.
!          bparnew(:,:,ik) = 0.
!          phi(:,:,ik) = 0.
!          apar(:,:,ik) = 0.
!          bpar(:,:,ik) = 0.
!       else
!          phinew(:,2:ntheta0,ik) = 0.
!          aparnew(:,2:ntheta0,ik) = 0.
!          bparnew(:,2:ntheta0,ik) = 0.
!          phi(:,2:ntheta0,ik) = 0.
!          apar(:,2:ntheta0,ik) = 0.
!          bpar(:,2:ntheta0,ik) = 0.
!       end if
!    end do

    do ik = 1, naky
       do it=1,ntheta0
          if (it == 1 .and. ik == 2) cycle
          phinew(:,it,ik) = 0.
          aparnew(:,it,ik) = 0.
          bparnew(:,it,ik) = 0.
          phi(:,it,ik) = 0.
          apar(:,it,ik) = 0.
          bpar(:,it,ik) = 0.          
       end do
    end do
    
    do ik = 1, naky
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid
             phiz(ig,it,ik) = cmplx(ranf()-0.5,ranf()-0.5)
          end do
          if (chop_side) then
             if (left) then
                phiz(:-1,it,ik) = 0.0
             else
                phiz(1:,it,ik) = 0.0
             end if
          end if
       end do
    end do

! reality condition for k_theta = 0 component:
    do it = 1, ntheta0/2
       phiz(:,it+(ntheta0+1)/2,1) = conjg(phiz(:,(ntheta0+1)/2+1-it,1))
    enddo

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       g(:,1,iglo) = g(:,1,iglo)-phiz(:,it,ik)*spec(is)%z*phiinit
       g(:,2,iglo) = g(:,2,iglo)-phiz(:,it,ik)*spec(is)%z*phiinit
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       where (forbid(:,il)) g(:,2,iglo) = 0.0
    end do
    gnew = g

  end subroutine ginit_nl4

  subroutine ginit_nl5
    use mp, only: proc0
    use species, only: spec
    use gs2_save, only: gs2_restore
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use dist_fn_arrays, only: g, gnew
    use le_grids, only: forbid
    use fields_arrays, only: phi, apar, bpar
    use fields_arrays, only: phinew, aparnew, bparnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx, il_idx
    use file_utils, only: error_unit
    use run_parameters, only: fphi, fapar, fbpar
    use ran
    implicit none
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phiz
    integer :: iglo, istatus
    integer :: ig, ik, it, is, il, ierr
    logical :: many = .true.
    
    call gs2_restore (g, scale, istatus, fphi, fapar, fbpar, many)
    if (istatus /= 0) then
       ierr = error_unit()
       if (proc0) write(ierr,*) "Error reading file: ", trim(restart_file)
       g = 0.
    end if

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       if (ik == 1) cycle

       g (:,1,iglo) = 0.
       g (:,2,iglo) = 0.
    end do
	
    phinew(:,:,2:naky) = 0.
    aparnew(:,:,2:naky) = 0.
    bparnew(:,:,2:naky) = 0.
    phi(:,:,2:naky) = 0.
    apar(:,:,2:naky) = 0.
    bpar(:,:,2:naky) = 0.

    do ik = 1, naky
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid
             phiz(ig,it,ik) = cmplx(ranf()-0.5,ranf()-0.5)
          end do
          if (chop_side) then
             if (left) then
                phiz(:-1,it,ik) = 0.0
             else
                phiz(1:,it,ik) = 0.0
             end if
          end if
       end do
    end do

    do ik = 1, naky
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid
             phiz(ig,it,ik) = cmplx(ranf()-0.5,ranf()-0.5)
          end do
          if (chop_side) then
             if (left) then
                phiz(:-1,it,ik) = 0.0
             else
                phiz(1:,it,ik) = 0.0
             end if
          end if
       end do
    end do

! reality condition for k_theta = 0 component:
    do it = 1, ntheta0/2
       phiz(:,it+(ntheta0+1)/2,1) = conjg(phiz(:,(ntheta0+1)/2+1-it,1))
    enddo

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       g(:,1,iglo) = g(:,1,iglo)-phiz(:,it,ik)*phiinit*spec(is)%z
       g(:,2,iglo) = g(:,2,iglo)-phiz(:,it,ik)*phiinit*spec(is)%z
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       where (forbid(:,il)) g(:,2,iglo) = 0.0
    end do
    gnew = g

  end subroutine ginit_nl5

  subroutine ginit_nl6
    use mp, only: proc0
    use species, only: spec
    use gs2_save, only: gs2_restore
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use dist_fn_arrays, only: g, gnew
    use le_grids, only: forbid
    use fields_arrays, only: phi, apar, bpar
    use fields_arrays, only: phinew, aparnew, bparnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx, il_idx
    use file_utils, only: error_unit
    use run_parameters, only: fphi, fapar, fbpar
    use ran
    implicit none
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phiz
    integer :: iglo, istatus
    integer :: ig, ik, it, is, il, ierr
    logical :: many = .true.
    
    call gs2_restore (g, scale, istatus, fphi, fapar, fbpar, many)
    if (istatus /= 0) then
       ierr = error_unit()
       if (proc0) write(ierr,*) "Error reading file: ", trim(restart_file)
       g = 0.
    end if

    gnew = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       if (ik == 1 .and. it == 2) cycle
       if (ik == 1 .and. it == ntheta0) cycle

       g (:,:,iglo) = gnew(:,:,iglo)
    end do
	
    phinew(:,2,1) = 0.   
    aparnew(:,2,1) = 0.
    bparnew(:,2,1) = 0.

    phi(:,2,1) = 0.
    apar(:,2,1) = 0.
    bpar(:,2,1) = 0.

    phinew(:,ntheta0,1) = 0.
    aparnew(:,ntheta0,1) = 0.
    bparnew(:,ntheta0,1) = 0.

    phi(:,ntheta0,1) = 0.
    apar(:,ntheta0,1) = 0.
    bpar(:,ntheta0,1) = 0.

    gnew = g

  end subroutine ginit_nl6

  subroutine ginit_nl7
    use species, only: spec, has_electron_species
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, theta0, reality
    use le_grids, only: forbid
    use dist_fn_arrays, only: g, gnew, vpa, vperp2
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use constants
    use ran
    implicit none
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phi, odd
    real, dimension (-ntgrid:ntgrid) :: dfac, ufac, tparfac, tperpfac, ct, st, c2t, s2t
    integer :: iglo
    integer :: ig, ik, it, il, is, j
    
    phi = 0.
    odd = 0.
    do it=1,ntheta0
       do ik=1,naky
          phi(:,it,ik) = cmplx(refac,imfac)*dphiinit
       end do
    end do

    do j = 1, 2
       ik = ikk(j)
       it = itt(j)
       do ig = -ntgrid, ntgrid
          phi(ig,it,ik) = cmplx(refac, imfac)
       end do
    end do

    odd = zi * phi
    
! reality condition for k_theta = 0 component:
    if (reality) then
       do it = 1, ntheta0/2
          phi(:,it+(ntheta0+1)/2,1) = conjg(phi(:,(ntheta0+1)/2+1-it,1))
          odd(:,it+(ntheta0+1)/2,1) = conjg(odd(:,(ntheta0+1)/2+1-it,1))
       enddo
    end if

    if (even) then
       ct = cos(theta)
       st = sin(theta)

       c2t = cos(2.*theta)
       s2t = sin(2.*theta)
    else
       ct = sin(theta)
       st = cos(theta)

       c2t = sin(2.*theta)
       s2t = cos(2.*theta)
    end if

    dfac     = den0   + den1 * ct + den2 * c2t
    ufac     = upar0  + upar1* st + upar2* s2t
    tparfac  = tpar0  + tpar1* ct + tpar2* c2t
    tperpfac = tperp0 + tperp1*ct + tperp2*c2t

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)       

       g(:,1,iglo) = phiinit* &
            ( dfac*spec(is)%dens0            * phi(:,it,ik) &
            + 2.*ufac* vpa(:,1,iglo)         * odd(:,it,ik) &
            + tparfac*(vpa(:,1,iglo)**2-0.5) * phi(:,it,ik) &
            +tperpfac*(vperp2(:,iglo)-1.)    * phi(:,it,ik))
       where (forbid(:,il)) g(:,1,iglo) = 0.0

       g(:,2,iglo) = phiinit* &
            ( dfac*spec(is)%dens0            * phi(:,it,ik) &
            + 2.*ufac* vpa(:,2,iglo)         * odd(:,it,ik) &
            + tparfac*(vpa(:,2,iglo)**2-0.5) * phi(:,it,ik) &
            +tperpfac*(vperp2(:,iglo)-1.)    * phi(:,it,ik))
       where (forbid(:,il)) g(:,2,iglo) = 0.0

    end do

    gnew = g

  end subroutine ginit_nl7

  subroutine ginit_kpar
    use species, only: spec, has_electron_species
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, theta0
    use le_grids, only: forbid
    use dist_fn_arrays, only: g, gnew, vpa, vperp2
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use constants
    use ran
    implicit none
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phi, odd
    real, dimension (-ntgrid:ntgrid) :: dfac, ufac, tparfac, tperpfac
    integer :: iglo
    integer :: ig, ik, it, il, is, j
    
    phi = 0.
    odd = 0.
    if (width0 > 0.) then
       do ig = -ntgrid, ntgrid
          phi(ig,:,:) = exp(-((theta(ig)-theta0(:,:))/width0)**2)*cmplx(refac, imfac)
       end do
    else
       do ig = -ntgrid, ntgrid
          phi(ig,:,:) = cmplx(refac, imfac)
       end do
    end if
    if (chop_side) then
       if (left) then
          phi(:-1,:,:) = 0.0
       else
          phi(1:,:,:) = 0.0
       end if
    end if

    odd = zi * phi
        
    dfac     = den0   + den1 * cos(theta) + den2 * cos(2.*theta) 
    ufac     = upar0  + upar1* sin(theta) + upar2* sin(2.*theta) 
    tparfac  = tpar0  + tpar1* cos(theta) + tpar2* cos(2.*theta) 
    tperpfac = tperp0 + tperp1*cos(theta) + tperp2*cos(2.*theta) 

! charge dependence keeps initial Phi from being too small
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)       

       ! TEMP FOR TESTING -- MAB
!       if (is == 2) then
!          ufac = 2.*upar0
!       else
!          ufac = upar0
!       end if

       g(:,1,iglo) = phiinit* &!spec(is)%z* &
            ( dfac                           * phi(:,it,ik) &
            + 2.*ufac* vpa(:,1,iglo)         * odd(:,it,ik) &
            + tparfac*(vpa(:,1,iglo)**2-0.5) * phi(:,it,ik) &
            +tperpfac*(vperp2(:,iglo)-1.)    * phi(:,it,ik))
       where (forbid(:,il)) g(:,1,iglo) = 0.0

       g(:,2,iglo) = phiinit* &!spec(is)%z* &
            ( dfac                           * phi(:,it,ik) &
            + 2.*ufac* vpa(:,2,iglo)         * odd(:,it,ik) &
            + tparfac*(vpa(:,2,iglo)**2-0.5) * phi(:,it,ik) &
            +tperpfac*(vperp2(:,iglo)-1.)    * phi(:,it,ik))
       where (forbid(:,il)) g(:,2,iglo) = 0.0

    end do

    if (has_electron_species(spec)) then
       call flae (g, gnew)
       g = g - gnew
    end if

    gnew = g
  end subroutine ginit_kpar

  subroutine ginit_gs
    use species, only: spec, has_electron_species
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, theta0, reality
    use le_grids, only: forbid
    use dist_fn_arrays, only: g, gnew, vpa, vperp2
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use constants
    use ran
    implicit none
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phi, odd
    integer :: iglo
    integer :: ig, ik, it, il, is
    real :: phase
    
    phi = 0.
    odd = 0.
    do ik=1,naky
       do it=1,ntheta0
          phase = 2.*pi*ranf()
          phi(:,it,ik) = cos(theta+phase)*cmplx(refac,imfac)
          odd(:,it,ik) = sin(theta+phase)*cmplx(refac,imfac) * zi
       end do
    end do
    
! reality condition for k_theta = 0 component:
    do it = 1, ntheta0/2
       phi(:,it+(ntheta0+1)/2,1) = conjg(phi(:,(ntheta0+1)/2+1-it,1))
       odd(:,it+(ntheta0+1)/2,1) = conjg(odd(:,(ntheta0+1)/2+1-it,1))
    enddo

! charge dependence keeps initial Phi from being too small
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)       

       g(:,1,iglo) = phiinit* &!spec(is)%z* &
            ( den1                         * phi(:,it,ik) &
            + 2.*upar1* vpa(:,1,iglo)      * odd(:,it,ik) &
            + tpar1*(vpa(:,1,iglo)**2-0.5) * phi(:,it,ik) &
            + tperp1*(vperp2(:,iglo)-1.)   * phi(:,it,ik))
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       
       g(:,2,iglo) = phiinit* &!spec(is)%z* &
            ( den1                         * phi(:,it,ik) &
            + 2.*upar1* vpa(:,2,iglo)      * odd(:,it,ik) &
            + tpar1*(vpa(:,2,iglo)**2-0.5) * phi(:,it,ik) &
            + tperp1*(vperp2(:,iglo)-1.)   * phi(:,it,ik))
       where (forbid(:,il)) g(:,2,iglo) = 0.0

    end do

    if (has_electron_species(spec)) then
       call flae (g, gnew)
       g = g - gnew
    end if

    gnew = g
  end subroutine ginit_gs

  subroutine ginit_test1
    use species, only: spec
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, akr
    use le_grids, only: e, forbid
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use constants
    implicit none
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phi
    integer :: iglo
    integer :: ig, ik, it, il, ie, is

    do ik = 1, naky
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid
             phi(ig,it,ik) = sin(akr(ig,it)*theta(ig))/(akr(ig,it))*zi
          end do
       end do
    end do

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       g(:,1,iglo) = -phi(:,it,ik)*spec(is)%z*phiinit*exp(-e(ie,is))
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       g(:,2,iglo) = g(:,1,iglo)
    end do
    gnew = g
  end subroutine ginit_test1

  subroutine ginit_xi
    use theta_grid, only: ntgrid, theta, bmag
    use le_grids, only: forbid, al
    use kt_grids, only: theta0
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use constants
    implicit none
    integer :: iglo
    integer :: ik, it, il, ie, is
    real, dimension(-ntgrid:ntgrid) :: xi

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       xi = sqrt(max(1.0-bmag*al(il),0.0))
       g(:,1,iglo) = xi*exp(-((theta-theta0(it,ik))/width0)**2)
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       g(:,2,iglo) = -g(:,1,iglo)
    end do
    gnew = g
  end subroutine ginit_xi

  subroutine ginit_xi2
    use theta_grid, only: ntgrid, theta, bmag
    use le_grids, only: forbid, al
    use kt_grids, only: theta0
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use constants
    implicit none
    integer :: iglo
    integer :: ik, it, il, ie, is
    real, dimension(-ntgrid:ntgrid) :: xi

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       xi = sqrt(max(1.0-bmag*al(il),0.0))
       g(:,1,iglo) = (1.0 - 3.0*xi*xi)*exp(-((theta-theta0(it,ik))/width0)**2)
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       g(:,2,iglo) = g(:,1,iglo)
    end do
    gnew = g
  end subroutine ginit_xi2

  subroutine ginit_rh
    use le_grids, only: forbid, e
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use constants
    implicit none
    integer :: iglo
    integer :: ik, it, il, ie, is

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       g(:,1,iglo) = exp(-e(ie,is))
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       g(:,2,iglo) = g(:,1,iglo)
    end do
    gnew = g
  end subroutine ginit_rh

  subroutine ginit_alf
    use theta_grid, only: theta
    use le_grids, only: forbid
    use dist_fn_arrays, only: g, gnew, vpa
    use gs2_layouts, only: g_lo, il_idx, is_idx
    use species, only: spec, electron_species

    implicit none
    integer :: iglo
    integer :: il, is

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       if (spec(is_idx(g_lo, iglo))%type == electron_species) cycle
       il = il_idx(g_lo,iglo)
       g(:,1,iglo) = sin(theta)*vpa(:,1,iglo)*spec(is)%z
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       g(:,2,iglo) = sin(theta)*vpa(:,2,iglo)*spec(is)%z
       where (forbid(:,il)) g(:,2,iglo) = 0.0
    end do
    g = phiinit * g 
    gnew = g
  end subroutine ginit_alf

  subroutine ginit_zero
    use dist_fn_arrays, only: g, gnew
    implicit none
    g = 0.0
    gnew = 0.0
  end subroutine ginit_zero

  subroutine ginit_test3
    use dist_fn_arrays, only: g, gnew, vpa
    use theta_grid, only: ntgrid, delthet, bmag
    use kt_grids, only: akx
    use theta_grid_params, only: eps, epsl, pk
    use gs2_layouts, only: g_lo, ik_idx, il_idx
    use mp, only: broadcast
    use constants
    implicit none
    integer :: iglo, ik, il
    real :: c1, c2

    call broadcast (epsl)
    call broadcast (eps)
    call broadcast (pk)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       if (any(vpa(-ntgrid:ntgrid-1,:,iglo) == 0.0)) then
          c1 = 0.0
          c2 = 0.0
       else
          c2 = -akx(ik)*epsl/eps/pk &
               *sum(delthet(-ntgrid:ntgrid-1)/bmag(-ntgrid:ntgrid-1))
          c1 = c2/sum(delthet(-ntgrid:ntgrid-1)/vpa(-ntgrid:ntgrid-1,1,iglo))
          c2 = c2/sum(delthet(-ntgrid:ntgrid-1)/vpa(-ntgrid:ntgrid-1,2,iglo))
       end if
       g(:,1,iglo) = -zi*akx(ik)*epsl/eps/pk*vpa(:,1,iglo)/bmag - zi*c1
       g(:,2,iglo) = -zi*akx(ik)*epsl/eps/pk*vpa(:,2,iglo)/bmag - zi*c2
    end do
    gnew = g
  end subroutine ginit_test3

  subroutine ginit_convect
    use dist_fn_arrays, only: g, gnew
    use theta_grid, only: theta
    use kt_grids, only: theta0
    use gs2_layouts, only: g_lo, it_idx, ik_idx
    use constants
    implicit none
    integer :: it, ik, iglo

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       g(:,1,iglo) = exp(cmplx(-(theta-theta0(it,ik))**2,k0*theta))
       g(:,2,iglo) = g(:,1,iglo)
    end do
    gnew = g
  end subroutine ginit_convect

  subroutine ginit_restart_file
    use dist_fn_arrays, only: g, gnew
    use gs2_save, only: gs2_restore
    use mp, only: proc0
    use file_utils, only: error_unit
    use run_parameters, only: fphi, fapar, fbpar
    implicit none
    integer :: istatus, ierr

    call gs2_restore (g, scale, istatus, fphi, fapar, fbpar)
    if (istatus /= 0) then
       ierr = error_unit()
       if (proc0) write(ierr,*) "Error reading file: ", trim(restart_file)
       g = 0.
    end if
    gnew = g

  end subroutine ginit_restart_file

  subroutine ginit_restart_many
    use dist_fn_arrays, only: g, gnew
    use gs2_save, only: gs2_restore
    use mp, only: proc0
    use file_utils, only: error_unit
    use run_parameters, only: fphi, fapar, fbpar
    implicit none
    integer :: istatus, ierr
    logical :: many = .true.

    call gs2_restore (g, scale, istatus, fphi, fapar, fbpar, many)

    if (istatus /= 0) then
       ierr = error_unit()
       if (proc0) write(ierr,*) "Error reading file: ", trim(restart_file)
       g = 0.
    end if
    gnew = g

  end subroutine ginit_restart_many

  subroutine ginit_restart_small
    use dist_fn_arrays, only: g, gnew
    use gs2_save, only: gs2_restore
    use mp, only: proc0
    use file_utils, only: error_unit
    use run_parameters, only: fphi, fapar, fbpar
    implicit none
    integer :: istatus, ierr
    logical :: many = .true.

    call ginit_noise

    call gs2_restore (g, scale, istatus, fphi, fapar, fbpar, many)
    if (istatus /= 0) then
       ierr = error_unit()
       if (proc0) write(ierr,*) "Error reading file: ", trim(restart_file)
       g = 0.
    end if
    g = g + gnew
    gnew = g 

  end subroutine ginit_restart_small

  subroutine ginit_restart_smallflat

    use gs2_save, only: gs2_restore
    use mp, only: proc0
    use file_utils, only: error_unit
    use species, only: spec
    use theta_grid, only: ntgrid 
    use kt_grids, only: naky, ntheta0, aky, reality
    use le_grids, only: forbid
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use run_parameters, only: fphi, fapar, fbpar
    use ran

    implicit none
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phi
    integer :: istatus, ierr
    logical :: many = .true.
    real :: a, b
    integer :: iglo
    integer :: ig, ik, it, il, is

    do it = 1, ntheta0
       do ik = 1, naky
          a = ranf()-0.5
          b = ranf()-0.5
          phi(:,it,ik) = cmplx(a,b)
       end do
    end do

    if (naky > 1 .and. aky(1) == 0.0) then
       phi(:,:,1) = phi(:,:,1)*zf_init
    end if
! reality condition for k_theta = 0 component:
    if (reality) then
       do it = 1, ntheta0/2
          phi(:,it+(ntheta0+1)/2,1) = conjg(phi(:,(ntheta0+1)/2+1-it,1))
       enddo
    end if

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       g(:,1,iglo) = -phi(:,it,ik)*spec(is)%z*phiinit
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       g(:,2,iglo) = g(:,1,iglo)
    end do
    gnew = g

    call gs2_restore (g, scale, istatus, fphi, fapar, fbpar, many)
    if (istatus /= 0) then
       ierr = error_unit()
       if (proc0) write(ierr,*) "Error reading file: ", trim(restart_file)
       g = 0.
    end if
    g = g + gnew
    gnew = g 

  end subroutine ginit_restart_smallflat

  subroutine reset_init

    ginitopt_switch = ginitopt_restart_many

  end subroutine reset_init

  subroutine flae (g, gavg)

    use species, only: spec, electron_species 
    use theta_grid, only: ntgrid, delthet, jacob
    use kt_grids, only: aky
    use gs2_layouts, only: g_lo, is_idx, ik_idx
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: gavg

    real :: wgt
    integer :: iglo
    
    gavg = 0.
    wgt = 1./sum(delthet*jacob)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc       
       if (spec(is_idx(g_lo, iglo))%type /= electron_species) cycle
       if (aky(ik_idx(g_lo, iglo)) /= 0.) cycle
       gavg(:,1,iglo) = sum(g(:,1,iglo)*delthet*jacob)*wgt
       gavg(:,2,iglo) = sum(g(:,2,iglo)*delthet*jacob)*wgt
    end do

  end subroutine flae

  subroutine init_vnmult (vnm)

    use gs2_save, only: init_vnm

    real, dimension (2) :: vnm
    integer :: istatus

    call init_vnm (vnm, istatus)

  end subroutine init_vnmult

end module init_g

