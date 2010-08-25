module init_g
  implicit none

  public :: ginit
  public :: init_init_g, finish_init_g
  public :: width0
  public :: tstart
  public :: reset_init
  public :: init_vnmult
  public :: new_field_init
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
       ginitopt_nl3r = 26, ginitopt_smallflat = 27, ginitopt_harris = 28, &
       ginitopt_recon3 = 29, ginitopt_zonal_only = 30, ginitopt_single_parallel_mode = 31
  real :: width0, dphiinit, phiinit, imfac, refac, zf_init, phifrac
  real :: den0, upar0, tpar0, tperp0
  real :: den1, upar1, tpar1, tperp1
  real :: den2, upar2, tpar2, tperp2
  real :: tstart, scale, apar0
  logical :: chop_side, left, even, new_field_init
  character(300) :: restart_file
  character (len=150) :: restart_dir
  integer, dimension(2) :: ikk, itt
  integer, dimension(3) :: ikkk,ittt

  ! <doc> These are used for the function ginit_single_parallel_mode, and specify the kparallel to initialize. In the case of  zero magnetic shear, of course, the box is periodic in the parallel direction, and so only harmonics of the box size (i.e. ikpar_init) are allowed  EGH</doc>

  integer :: ikpar_init
  real :: kpar_init

  ! RN> for recon3
  real :: phiinit0 ! amplitude of equilibrium
  real :: a0,b0 ! amplitude of equilibrium Apara or By 
  ! if b0 /= 0, u0 is rescaled to give this By amplitude by overriding a0
  ! if a0 /= 0, u0 is rescaled to give this Apara amplitude
  logical :: null_phi, null_bpar, null_apar ! nullify fields
  ! if use_{Phi,Bpar,Apar} = F, override corresponding flags
  integer :: adj_spec ! adjust input parameter of this spec
  ! if = 0, just warn if given conditions are not satisfied
  character (len=20) :: eq_type = ''
  ! equilibrium type = 'sinusoidal', 'porcelli', 'doubleharris'
  real :: prof_width=-.05
  ! width of porcelli, doubleharris profile
  ! if negative, this gives the ratio to the box size
  integer :: eq_mode_n=0, eq_mode_u=1
  ! mode number in x for sinusoidal equilibrium
  logical :: input_check_recon=.false.
  complex :: nkxy_pt(3), ukxy_pt(3)
  ! <RN
  
  logical :: debug = .false.
  logical :: initialized = .false.

contains

  subroutine init_init_g
    use gs2_save, only: init_save
    use gs2_layouts, only: init_gs2_layouts
    use mp, only: proc0, broadcast, job
    implicit none
!    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_gs2_layouts

    if (proc0) call read_parameters

    ! prepend restart_dir to restart_file
    ! append trailing slash if not exists
    if(restart_dir(len_trim(restart_dir):) /= "/") &
         restart_dir=trim(restart_dir)//"/"
    restart_file=trim(restart_dir)//trim(restart_file)

    ! MAB - allows for ensemble averaging of multiple flux tube calculations
    ! job=0 if not doing multiple flux tube calculations, so phiinit unaffected
    phiinit = phiinit * (job*phifrac+1.0)

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
    call broadcast (phifrac)
    call broadcast (dphiinit)
    call broadcast (zf_init)
    call broadcast (apar0)
    call broadcast (tstart)
    call broadcast (chop_side)
    call broadcast (even)
    call broadcast (left)
    call broadcast (restart_file)
    call broadcast (ikk)
    call broadcast (itt) 
    call broadcast (ikkk)
    call broadcast (ittt) 
    call broadcast (scale)
    call broadcast (new_field_init)

    ! RN>
    call broadcast (phiinit0)
    call broadcast (a0)
    call broadcast (b0)
    call broadcast (null_phi)
    call broadcast (null_bpar)
    call broadcast (null_apar)
    call broadcast (adj_spec)
    call broadcast (eq_type)
    call broadcast (prof_width)
    call broadcast (eq_mode_n)
    call broadcast (eq_mode_u)
    call broadcast (input_check_recon)
    call broadcast (nkxy_pt)
    call broadcast (ukxy_pt)
    ! <RN
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
    case (ginitopt_single_parallel_mode)
       call ginit_single_parallel_mode
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
    case (ginitopt_zonal_only)
       call ginit_restart_zonal_only
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
    case (ginitopt_recon3)
       call ginit_recon3
    end select
  end subroutine ginit

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, run_name, input_unit_exist
    use text_options, only: text_option, get_option_value
    implicit none

    type (text_option), dimension (31), parameter :: ginitopts = &
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
            text_option('recon', ginitopt_recon), &
            text_option('recon3', ginitopt_recon3), &
            text_option('zonal_only', ginitopt_zonal_only), &
            text_option('single_parallel_mode', ginitopt_single_parallel_mode) &
            /)
    character(20) :: ginit_option
    namelist /init_g_knobs/ ginit_option, width0, phiinit, chop_side, &
         restart_file, restart_dir, left, ikk, itt, scale, tstart, zf_init, &
         den0, upar0, tpar0, tperp0, imfac, refac, even, &
         den1, upar1, tpar1, tperp1, &
         den2, upar2, tpar2, tperp2, dphiinit, apar0, &
         new_field_init, &
         phiinit0, a0, b0, null_phi, null_bpar, null_apar, adj_spec, &
         eq_type, prof_width, eq_mode_u, eq_mode_n, &
         input_check_recon, nkxy_pt, ukxy_pt, &
         ikkk, ittt, phifrac, ikpar_init, kpar_init

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
    apar0 = 0.
    chop_side = .true.
    left = .true.
    even = .true.
    new_field_init = .true.
    ikk(1) = 1
    ikk(2) = 2
    itt(1) = 1
    itt(2) = 2
    phifrac = 0.1

    ! >RN
    ikkk(1) = 1 ; ikkk(2) = 2 ; ikkk(3) = 2
    ittt(1) = 1 ; ittt(2) = 2 ; ittt(3) = 2
    phiinit0 = 0.
    a0 = 0.
    b0 = 0.
    null_phi = .false.
    null_bpar = .false.
    null_apar = .false.
    adj_spec = 0
    eq_type = 'sinusoidal'
    prof_width=-0.1
    eq_mode_n=0
    eq_mode_u=1
    input_check_recon=.false.
    nkxy_pt=cmplx(0.,0.)
    ukxy_pt=cmplx(0.,0.)

    ikpar_init = 0
    kpar_init = 0.0

    ! <RN
    restart_file = trim(run_name)//".nc"
    restart_dir = "./"
    in_file = input_unit_exist ("init_g_knobs", exist)
!    if (exist) read (unit=input_unit("init_g_knobs"), nml=init_g_knobs)
    if (exist) read (unit=in_file,nml=init_g_knobs)

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

	!<doc> Initialize with a single parallel mode. Only makes sense in a linear calculation. EGH </doc> 

  subroutine ginit_single_parallel_mode
    use species, only: spec, tracer_species
    use theta_grid, only: ntgrid, shat, theta 
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

    if (abs(shat) < 1.0e-05) then 
      if (ikpar_init+1 > ntgrid) then 
        write (*,*) "Error: this value of k_parallel is too large. Increase ntheta or decrease ikpar_init."
        exit
      end if
      kpar_init = ikpar_init
    end if


    do it = 1, ntheta0
       do ik = 1, naky
          do ig = -ntgrid, ntgrid
              !<doc> Set the field to cos(kpar_init*theta), where we remember that the gridpoints are not necessarily evenly spaced in the parallel direction, so we use theta(ig)</doc>
            a = cos(kpar_init * theta(ig)) 
            b = cos(kpar_init * theta(ig))
              

!              a = ranf()-0.5
!              b = ranf()-0.5
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
    !<doc> reality condition for ky = 0 component: </doc>
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

  end subroutine ginit_single_parallel_mode

  
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
    use run_parameters, only: k0
    implicit none
    complex, dimension (ntheta0,naky) :: phi
    integer :: iglo
!    integer :: ig, ik, it, il, is, j, j1
    integer :: ik, it, il, is, j, j1
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
!    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phiz
    integer :: iglo, istatus
!    integer :: ig, ik, it, is, il, ierr
    integer :: ik, it, ierr
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
!    integer :: ig, ik, it, il, is, j
    integer :: ig, ik, it, il, is
    
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
!    integer :: ig, ik, it, il, is
    integer :: ik, it, il, is
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
    use gs2_layouts, only: g_lo, it_idx, il_idx, ie_idx, is_idx
    use constants
    implicit none
    integer :: iglo
!    integer :: it, il, ie, is
    integer :: il, ie, is

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
    use run_parameters, only: k0
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

  subroutine ginit_recon3
    use mp, only: proc0
    use mp, only: iproc
    use species, only: nspec, spec, has_electron_species
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, nakx => ntheta0, akx, aky, reality
    use kt_grids, only: nx,ny
    use dist_fn_arrays, only: g, gnew, vpa, vperp2
    use dist_fn_arrays, only: aj0,aj1,vperp2
    use dist_fn_arrays, only: kperp2
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use gs2_transforms, only: inverse2,transform2
    use run_parameters, only: beta
    use le_grids, only: integrate_moment
    use dist_fn, only: getmoms_notgc, get_init_field
    use dist_fn, only: mom_coeff, mom_shift_para, mom_shift_perp
    use dist_fn, only: gamtot, gamtot1, gamtot2
    use constants, only: pi, zi
    use file_utils, only: error_unit, get_unused_unit
    implicit none
    integer :: iglo
!    integer :: ig, ik, it, is, is2
    integer :: ig, ik, it, is
    integer :: i,j

    ! nkxyz and ukxyz determine profiles in kx-ky plane
    ! nkxyz is for N, T, and ukxyz is for U
    ! [do not control amplitude by these variables]
    complex :: nkxyz(-ntgrid:ntgrid,nakx,naky)
    complex :: ukxyz(-ntgrid:ntgrid,nakx,naky)
    complex :: nkxyz_eq(-ntgrid:ntgrid,nakx,naky)
    complex :: ukxyz_eq(-ntgrid:ntgrid,nakx,naky)
    complex, allocatable :: g_eq(:,:,:)

    ! equilibrium and perturbation
    complex :: nkxy_eq(nakx,naky), ukxy_eq(nakx,naky)

    ! *fac determine z profile
    ! [do not control amplitude by these variables]
    real :: dfac(-ntgrid:ntgrid), ufac(-ntgrid:ntgrid)
    real :: tparfac(-ntgrid:ntgrid), tperpfac(-ntgrid:ntgrid)
    real :: ct(-ntgrid:ntgrid), st(-ntgrid:ntgrid)
    real :: c2t(-ntgrid:ntgrid), s2t(-ntgrid:ntgrid)

    ! normalization
    real :: fac
    real, allocatable :: nnrm(:,:,:),unrm(:,:,:)
    real, allocatable :: tparanrm(:,:,:),tperpnrm(:,:,:)

    real :: check(3)=0.,current=0.
    character (len=2) :: cfit='A0'
    complex, allocatable :: dens(:,:,:,:),upar(:,:,:,:)
    complex, allocatable :: tpar(:,:,:,:),tper(:,:,:,:)
    complex, allocatable :: phi(:,:,:),apar(:,:,:),bpar(:,:,:)
    complex, allocatable :: jpar(:,:,:)

    real :: save_dens0, save_tperp0, save_u0
    real :: ratio

    ! real space profile to be Fourier transformed
    real :: xx,dx,lx,ly
    integer, parameter :: nfx=2**10
    real, allocatable :: nxy(:,:),uxy(:,:)
    real, allocatable :: phixy(:,:,:),aparxy(:,:,:),bparxy(:,:,:)
    real, allocatable :: jparxy(:,:,:)
    real, allocatable :: densxy(:,:,:,:),uparxy(:,:,:,:)
    real, allocatable :: tparxy(:,:,:,:),tperxy(:,:,:,:)

    complex, allocatable :: by(:,:,:)
    real, allocatable :: byxy(:,:,:)

    integer :: nnx,nny
    integer :: unit
    real :: xmax, bymax, xa, xl1, xl2
    real :: fmin, fmax

    if(debug.and.proc0) write(6,*) 'Initialization recon3'

!!! adjust input parameters to kill initial field if wanted
    if(debug.and.proc0) write(6,*) 'check parameters'

    do is=1,nspec
       if(adj_spec == is) cycle
       check(1)=check(1)+spec(is)%z*spec(is)%dens*spec(is)%dens0
       check(2)=check(2)+spec(is)%dens*spec(is)%temp* &
            & (spec(is)%dens0+spec(is)%tperp0)
       check(3)=check(3)+spec(is)%z*spec(is)%dens*spec(is)%stm*spec(is)%u0
    end do

    if(adj_spec == 0) then
       ! just warn if fields don't satisfy given conditions
       if(null_phi.and.null_bpar) then
          if(sum(check(1:2)) /= 0.) then
             if(proc0) write(6,'(a)') &
                  'WARNING: finite Phi or Bpar in initial condition'
          endif
       else if(null_bpar.and..not.null_phi) then
          ratio=gamtot1(0,eq_mode_n+1,1)/gamtot(0,eq_mode_n+1,1)
          if(check(1)/check(2) /= ratio) then
             if(proc0) write(6,'(a)') &
                  'WARNING: finite Bpar in initial condition'
          endif
       else if(null_phi.and..not.null_bpar) then
          ratio=-(2.*gamtot2(0,eq_mode_n+1,1)+2./beta) &
               & /gamtot1(0,eq_mode_n+1,1)
          if(check(1)/check(2) /= ratio) then
             if(proc0) write(6,'(a)') &
                  'WARNING: finite Bpar in initial condition'
          endif
       endif
       if(null_apar) then
          if(check(3) /= 0.) then
             if(proc0) write(6,'(a)') &
                  'WARNING: finite Apar in initial condition'
          endif
       endif
    else
       ! adjust input parameter to satisfy given conditions
       if(null_phi.and.null_bpar) then
          save_dens0=spec(adj_spec)%dens0
          save_tperp0=spec(adj_spec)%tperp0
          spec(adj_spec)%dens0=-check(1)/(spec(adj_spec)%z*spec(adj_spec)%dens)
          spec(adj_spec)%tperp0=-spec(adj_spec)%dens0 &
               & -check(2)/(spec(adj_spec)%dens*spec(adj_spec)%temp)

          if(spec(adj_spec)%dens0 /= save_dens0) then
             if(proc0) write(6,'(a,i0,a,f10.2,a)') &
                  & 'WARNING: Initial density of spec=', spec(adj_spec)%type, &
                  & ' is adjusted to ', spec(adj_spec)%dens0, &
                  & ' to kill Phi and Bpar'
          endif
          if(spec(adj_spec)%tperp0 /= save_tperp0) then
             if(proc0) write(6,'(a,i0,a,f10.2,a)') &
                  & 'WARNING: Initial Tperp of spec=', spec(adj_spec)%type, &
                  & ' is adjusted to ', spec(adj_spec)%tperp0, &
                  & ' to kill Phi and Bpar'
          endif
       else if(null_bpar.and..not.null_phi.and.eq_type.eq.'sinusoidal') then
          save_tperp0=spec(adj_spec)%tperp0
          check(1)=check(1)+ &
               & spec(adj_spec)%z*spec(adj_spec)%dens*spec(adj_spec)%dens0
          ratio=gamtot1(0,eq_mode_n+1,1)/gamtot(0,eq_mode_n+1,1)

          spec(adj_spec)%tperp0=(-check(1)*ratio-check(2)) &
               & /(spec(adj_spec)%dens*spec(adj_spec)%temp) &
               & -spec(adj_spec)%dens0
          if(spec(adj_spec)%tperp0 /= save_tperp0) then
             if(proc0) write(6,'(a,i0,a,f10.2,a)') &
                  & 'WARNING: Initial Tperp of spec=', spec(adj_spec)%type, &
                  & ' is adjusted to ', spec(adj_spec)%tperp0, &
                  & ' to kill Bpar'
          endif
       else if(null_phi.and..not.null_bpar.and.eq_type.eq.'sinusoidal') then
          save_tperp0=spec(adj_spec)%tperp0
          check(1)=check(1)+ &
               & spec(adj_spec)%z*spec(adj_spec)%dens*spec(adj_spec)%dens0
          ratio=-(2.*gamtot2(0,eq_mode_n+1,1)+2./beta) &
               & /gamtot1(0,eq_mode_n+1,1)

          spec(adj_spec)%tperp0=(-check(1)*ratio-check(2)) &
               & /(spec(adj_spec)%dens*spec(adj_spec)%temp) &
               & -spec(adj_spec)%dens0

          if(spec(adj_spec)%tperp0 /= save_tperp0) then
             if(proc0) write(6,'(a,i0,a,f10.2,a)') &
                  & 'WARNING: Initial Tperp of spec=', spec(adj_spec)%type, &
                  & ' is adjusted to ', spec(adj_spec)%tperp0, &
                  & ' to kill Phi'
          endif
       endif
    
       if (null_apar) then
          save_u0=spec(adj_spec)%u0
          spec(adj_spec)%u0=-check(3)/ &
               & (spec(adj_spec)%z*spec(adj_spec)%dens*spec(adj_spec)%stm)
          if(spec(adj_spec)%u0 /= save_u0) then
             if(proc0) write(6,'(a,i0,a,f10.2,a)') &
                  & 'WARNING: Initial U of spec=', spec(adj_spec)%type, &
                  & ' is adjusted to ', spec(adj_spec)%u0, &
                  & ' to kill Apar'
          endif
       endif
    endif

!!! initialize
    if(debug.and.proc0) write(6,*) 'initialize variable'
    nkxyz(-ntgrid:ntgrid,1:nakx,1:naky)=cmplx(0.,0.)
    ukxyz(-ntgrid:ntgrid,1:nakx,1:naky)=cmplx(0.,0.)
    nkxyz_eq(-ntgrid:ntgrid,1:nakx,1:naky)=cmplx(0.,0.)
    ukxyz_eq(-ntgrid:ntgrid,1:nakx,1:naky)=cmplx(0.,0.)

    nkxy_eq(1:nakx,1:naky)=cmplx(0.,0.)
    ukxy_eq(1:nakx,1:naky)=cmplx(0.,0.)

    dfac(-ntgrid:ntgrid)=1.
    ufac(-ntgrid:ntgrid)=1.
    tparfac(-ntgrid:ntgrid)=1.
    tperpfac(-ntgrid:ntgrid)=1.

!!! equilibrium
    if(phiinit0 /= 0.) then
       if(nakx==1 .and. naky==1) then ! grid_option = single case
          nkxy_eq(1,1)=cmplx(.5,0.)
          ukxy_eq(1,1)=cmplx(.5,0.)
       else
          if(debug.and.proc0) write(6,*) 'set equilibrium profile'
          allocate(nxy(nfx,ny),uxy(nfx,ny))
!          lx=2.*pi*x0; ly=2.*pi*y0
          lx=2.*pi/(akx(2)-akx(1)); ly=2.*pi/(aky(2)-aky(1))
          dx=lx/nfx
          ! if width is negative, it gives the ratio to the box size
          if(prof_width < 0. ) prof_width=-prof_width*lx
          select case (eq_type)
          case ('sinusoidal')
             ! this gives n,Tpara,Tperp \propto cos^2 (2 pi/Lx)
             ! nxy(eq_mode1(1),eq_mode1(2))=cmplx(.25, 0.)
             ! this gives Apara \propto cos(2 pi/Lx), By \propto sin(2 pi x/Lx)
             ! uxy(eq_mode2(1),eq_mode2(2))=cmplx(.5, 0.)
             do it=1,nfx
                xx=dx*(it-1)
                nxy(it,1:ny)=0.
                uxy(it,1:ny)=cos(2.*pi/lx*xx*eq_mode_u)
             end do
          case ('porcelli')
             do it=1,nfx
                xx=dx*(it-1)
                nxy(it,1:ny)=0.
                uxy(it,1:ny)=1./cosh((xx-.5*lx)/prof_width)**2 &
                     & * (tanh(2.*pi/lx*xx)**2+tanh(2.*pi/lx*(xx-lx))**2 &
                     & - tanh(2.*pi)**2) / (2.*tanh(pi)**2-tanh(2.*pi)**2)
             end do
          case ('doubleharris')
             do it=1,nfx
                fmin=tanh(.75*lx/prof_width)-tanh(.25*lx/prof_width)
                fmax=2.*tanh(.25*lx/prof_width)
                xx=dx*(it-1)
                nxy(it,1:ny)=0.
                uxy(it,1:ny)= 2.*( &
                     & tanh((xx-.25*lx)/prof_width) - &
                     & tanh((xx-.75*lx)/prof_width) - fmin)/(fmax-fmin)-1.
             end do
          case default
             if(proc0) write(error_unit(),'(2a)') &
                  & 'ERROR: Invalid equilibrium type', eq_type
             stop
          end select
          ! subtract (0,0) mode
          ! since it (constant part of potential) does nothing
          do ik=1,ny
             uxy(1:nfx,ik)=uxy(1:nfx,ik) &
                  - sum(uxy(1:nfx,ik))/nfx
          end do
          call inverse2(nxy, nkxy_eq, ny, nfx)
          call inverse2(uxy, ukxy_eq, ny, nfx)
          deallocate(nxy,uxy)
       endif
    endif

    do ig = -ntgrid, ntgrid
       nkxyz(ig,1:nakx,1:naky) = phiinit0*nkxy_eq(1:nakx,1:naky)
       ukxyz(ig,1:nakx,1:naky) = phiinit0*ukxy_eq(1:nakx,1:naky)
    end do
    if(.not.(nakx==1 .and. naky==1)) then ! grid_option /= single case
       nkxyz(-ntgrid:ntgrid,1,1)=cmplx(0.,0.)
       ukxyz(-ntgrid:ntgrid,1,1)=cmplx(0.,0.)
    endif

!!! save equilibrium profile
    nkxyz_eq(-ntgrid:ntgrid,1:nakx,1:naky)= &
         & nkxyz(-ntgrid:ntgrid,1:nakx,1:naky)
    ukxyz_eq(-ntgrid:ntgrid,1:nakx,1:naky)= &
         & ukxyz(-ntgrid:ntgrid,1:nakx,1:naky)

!!! perturbation
    if(phiinit /= 0.) then
       if(debug.and.proc0) write(6,*) 'set perturbation profile'
       do j = 1, 3
          if(ikkk(j).eq.0) ukxy_pt(j)=.5*ukxy_pt(j) !reality
          ik = ikkk(j)+1
          if(ittt(j) >= 0) then
             it = ittt(j)+1
          else
             it = nakx + ittt(j) + 1
          endif
          do ig = -ntgrid, ntgrid
             ukxyz(ig,it,ik) = ukxyz(ig,it,ik) + phiinit*ukxy_pt(j)
          end do
       end do
    endif

!!! reality condition for k_theta = 0 component:
    if(debug.and.proc0) write(6,*) 'set reality condition'
    if (reality) then
       do it = 1, nakx/2
          nkxyz(-ntgrid:ntgrid,it+(nakx+1)/2,1) = &
               & conjg(nkxyz(-ntgrid:ntgrid,(nakx+1)/2+1-it,1))
          ukxyz(-ntgrid:ntgrid,it+(nakx+1)/2,1) = &
               & conjg(ukxyz(-ntgrid:ntgrid,(nakx+1)/2+1-it,1))
          nkxyz_eq(-ntgrid:ntgrid,it+(nakx+1)/2,1) = &
               & conjg(nkxyz_eq(-ntgrid:ntgrid,(nakx+1)/2+1-it,1))
          ukxyz_eq(-ntgrid:ntgrid,it+(nakx+1)/2,1) = &
               & conjg(ukxyz_eq(-ntgrid:ntgrid,(nakx+1)/2+1-it,1))
       enddo
    end if

!!! parallel profile
    if(debug.and.proc0) write(6,*) 'set parallel profile'
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

    dfac     = dfac     + den1   * ct + den2   * c2t
    ufac     = ufac     + upar1  * st + upar2  * s2t
    tparfac  = tparfac  + tpar1  * ct + tpar2  * c2t
    tperpfac = tperpfac + tperp1 * ct + tperp2 * c2t

!!! normalization
    if(debug.and.proc0) write(6,*) 'normalization'
    if(b0 /= 0.) then
       if(.not.(nakx==1 .and. naky==1)) then ! grid_option /= single case
          if(eq_type == 'porcelli') then
             xmax=.5*prof_width*log(2.+sqrt(3.))+.5*lx
             xmax=get_xmax(xmax)
             xa=(xmax-.5*lx)/prof_width
             xl1=2.*pi/lx*xmax; xl2=xl1-2.*pi
             bymax=(-2./prof_width*sinh(xa)/cosh(xa)**3* &
                  & (tanh(xl1)**2+tanh(xl2)**2-tanh(2.*pi)**2) &
                  & + 1./cosh(xa)**2*4.*pi/lx* &
                  & (sinh(xl1)/cosh(xl1)**3+sinh(xl2)/cosh(xl2)**3) ) &
                  & / (2.*tanh(pi)**2-tanh(2.*pi)**2)
          else if(eq_type == 'doubleharris') then
             bymax=2.*tanh(.5*lx/prof_width)**2/(prof_width*(fmax-fmin))
          else
             bymax=akx(eq_mode_u+1)
          endif
       else
          bymax=sqrt(kperp2(0,1,1))
       endif
       a0 = b0/bymax
       cfit='B0'
    endif

    allocate(nnrm(nakx,naky,nspec),unrm(nakx,naky,nspec))
    allocate(tparanrm(nakx,naky,nspec),tperpnrm(nakx,naky,nspec))
    do is=1,nspec
       nnrm(:,:,is)=mom_coeff(:,:,is,1)
       unrm(:,:,is)=2.*mom_coeff(:,:,is,2)
       tparanrm(:,:,is)=2.*(mom_coeff(:,:,is,4)- &
            & 2.*mom_coeff(:,:,is,2)*mom_shift_para(:,:,is))
       tperpnrm(:,:,is)=2.*(mom_coeff(:,:,is,8)- &
             & mom_coeff(:,:,is,6)*mom_shift_perp(:,:,is))

       if(a0 /= 0.) then
          ! rescale parallel momentum to obtain a given amplitude of Apar
          current=sum(spec(1:nspec)%z*spec(1:nspec)%dens &
               & *spec(1:nspec)%stm*spec(1:nspec)%u0)
          if(current==0.) then
             if(proc0) write(error_unit(),'(a)') &
                  & 'ERROR in init_g: invalid input a0, u0'
             stop
          endif
          do it=1,nakx
             do ik=1,naky
                if(cabs(ukxyz(0,it,ik)) /= 0. .and. spec(is)%u0 /= 0.) then
                   fac=2.*beta*current/kperp2(0,it,ik)/a0
                   unrm(it,ik,is)=unrm(it,ik,is)*fac
                   if(proc0) write(6,'(a,a2,a,3(i3,1x),f20.12)') &
                        & 'INFO: u0 is rescaled to fit ',cfit,' input', &
                        & it,ik,is,spec(is)%u0/fac
                end if
             end do
          end do
       end if
    end do

    if(debug.and.proc0) write(6,*) 'calculate g'
    allocate(g_eq(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)

       g(:,1,iglo) = ( &
            & ( dfac*spec(is)%dens0 / nnrm(it,ik,is) &
            & + tparfac*spec(is)%tpar0 / tparanrm(it,ik,is) &
            & * (vpa(:,1,iglo)**2-mom_shift_para(it,ik,is)) &
            & + tperpfac*spec(is)%tperp0 / tperpnrm(it,ik,is) &
            & * (vperp2(:,iglo)-mom_shift_perp(it,ik,is)) &
            & ) * nkxyz(:,it,ik) &
            & + ufac*2.*vpa(:,1,iglo)*spec(is)%u0 / unrm(it,ik,is) &
            & * ukxyz(:,it,ik) &
            )
       g(:,2,iglo) = ( &
            & ( dfac*spec(is)%dens0 / nnrm(it,ik,is) &
            & + tparfac*spec(is)%tpar0 / tparanrm(it,ik,is) &
            & * (vpa(:,2,iglo)**2-mom_shift_para(it,ik,is)) &
            & + tperpfac*spec(is)%tperp0 / tperpnrm(it,ik,is) &
            & * (vperp2(:,iglo)-mom_shift_perp(it,ik,is)) &
            & ) * nkxyz(:,it,ik) &
            & + ufac*2.*vpa(:,2,iglo)*spec(is)%u0 / unrm(it,ik,is) &
            & * ukxyz(:,it,ik) &
            )

       g_eq(:,1,iglo) = ( &
            & ( dfac*spec(is)%dens0 / nnrm(it,ik,is) &
            & + tparfac*spec(is)%tpar0 / tparanrm(it,ik,is) &
            & * (vpa(:,1,iglo)**2-mom_shift_para(it,ik,is)) &
            & + tperpfac*spec(is)%tperp0 / tperpnrm(it,ik,is) &
            & * (vperp2(:,iglo)-mom_shift_perp(it,ik,is)) &
            & ) * nkxyz_eq(:,it,ik) &
            & + ufac*2.*vpa(:,1,iglo)*spec(is)%u0 / unrm(it,ik,is) &
            & * ukxyz_eq(:,it,ik) &
            )
       g_eq(:,2,iglo) = ( &
            & ( dfac*spec(is)%dens0 / nnrm(it,ik,is) &
            & + tparfac*spec(is)%tpar0 / tparanrm(it,ik,is) &
            & * (vpa(:,2,iglo)**2-mom_shift_para(it,ik,is)) &
            & + tperpfac*spec(is)%tperp0 / tperpnrm(it,ik,is) &
            & * (vperp2(:,iglo)-mom_shift_perp(it,ik,is)) &
            & ) * nkxyz_eq(:,it,ik) &
            & + ufac*2.*vpa(:,2,iglo)*spec(is)%u0 / unrm(it,ik,is) &
            & * ukxyz_eq(:,it,ik) &
            )
    end do

    deallocate(nnrm,unrm,tparanrm,tperpnrm)

!!! store equilibrium fields
    allocate(phi(-ntgrid:ntgrid,1:nakx,1:naky))
    allocate(apar(-ntgrid:ntgrid,1:nakx,1:naky))
    allocate(bpar(-ntgrid:ntgrid,1:nakx,1:naky))
    allocate(jpar(-ntgrid:ntgrid,1:nakx,1:naky))
    allocate(dens(-ntgrid:ntgrid,nakx,naky,nspec))
    allocate(upar(-ntgrid:ntgrid,nakx,naky,nspec))
    allocate(tpar(-ntgrid:ntgrid,nakx,naky,nspec))
    allocate(tper(-ntgrid:ntgrid,nakx,naky,nspec))

    phi(:,:,:)=cmplx(0.,0.); apar(:,:,:)=cmplx(0.,0.);
    bpar(:,:,:)=cmplx(0.,0.); jpar(:,:,:)=cmplx(0.,0.)
    dens(:,:,:,:)=cmplx(0.,0.); upar(:,:,:,:)=cmplx(0.,0.)
    tpar(:,:,:,:)=cmplx(0.,0.); tper(:,:,:,:)=cmplx(0.,0.)

    if(.not.(nakx==1 .and. naky==1)) then ! grid_option /= single case
       if(debug.and.proc0) write(6,*) 'save equilibrium profile'
       if(nx>nakx) then
          nnx=nx
       else ! nx=nakx=1
          nnx=(3*nakx+1)/2
       endif
       if(ny>naky) then
          nny=ny
       else ! ny=naky=1
          nny=3*naky
       endif

       allocate(phixy(nnx,nny,-ntgrid:ntgrid))
       allocate(aparxy(nnx,nny,-ntgrid:ntgrid))
       allocate(bparxy(nnx,nny,-ntgrid:ntgrid))
       allocate(jparxy(nnx,nny,-ntgrid:ntgrid))
       allocate(densxy(nnx,nny,-ntgrid:ntgrid,nspec))
       allocate(uparxy(nnx,nny,-ntgrid:ntgrid,nspec))
       allocate(tparxy(nnx,nny,-ntgrid:ntgrid,nspec))
       allocate(tperxy(nnx,nny,-ntgrid:ntgrid,nspec))
       phixy (:,:,:)=0.; aparxy(:,:,:)=0.
       bparxy(:,:,:)=0.; jparxy(:,:,:)=0.
       densxy(:,:,:,:)=0.; uparxy(:,:,:,:)=0.
       tparxy(:,:,:,:)=0.; tperxy(:,:,:,:)=0.

       allocate(by(-ntgrid:ntgrid,1:nakx,1:naky))
       allocate(byxy(nnx,nny,-ntgrid:ntgrid))
       by(:,:,:)=cmplx(0.,0.)
       byxy(:,:,:)=0.

       ! get equilibrium fields
       gnew(:,:,:)=g_eq(:,:,:)
       call get_init_field(phi,apar,bpar)
       call getmoms_notgc(dens,upar,tpar,tper,jpar=jpar)
       call transform2(phi,phixy,nny,nnx)
       call transform2(apar,aparxy,nny,nnx)
       call transform2(bpar,bparxy,nny,nnx)
       call transform2(jpar,jparxy,nny,nnx)
       do is=1,nspec
          call transform2(dens(:,:,:,is),densxy(:,:,:,is),nny,nnx)
          call transform2(upar(:,:,:,is),uparxy(:,:,:,is),nny,nnx)
          call transform2(tpar(:,:,:,is),tparxy(:,:,:,is),nny,nnx)
          call transform2(tper(:,:,:,is),tperxy(:,:,:,is),nny,nnx)
       end do

       ! store equilibrium fields
       if(proc0) then
          call get_unused_unit(unit)
          open(unit,file='equilibrium.dat',form='unformatted')
          write(unit) phixy(:,:,0),aparxy(:,:,0),bparxy(:,:,0),jparxy(:,:,0)
          write(unit) densxy(:,:,0,:),uparxy(:,:,0,:), &
               & tparxy(:,:,0,:),tperxy(:,:,0,:)
          close(unit)
       endif
    endif

    deallocate(g_eq)

    ! now store the actual g in gnew
    gnew(:,:,:) = g(:,:,:)
    call get_init_field(phi,apar,bpar)

!!! check
    if(debug.and.proc0) write(6,*) 'check'
    if(input_check_recon) then

       call getmoms_notgc(dens,upar,tpar,tper)
       call get_init_field(phi,apar,bpar)

       if(nakx==1 .and. naky==1) then ! grid_option = single case
          if(proc0) then
             do is=1,nspec
                write(6,'(" spec = ",i0)') is
                write(6,'(" upar=",2e25.17)') &
                     & real(upar(0,1,1,is)),aimag(upar(0,1,1,is))
                write(6,'(" dens=",2e25.17)') &
                     & real(dens(0,1,1,is)),aimag(dens(0,1,1,is))
                write(6,'(" tpar=",2e25.17)') &
                     & real(tpar(0,1,1,is)),aimag(tpar(0,1,1,is))
                write(6,'(" tper=",2e25.17)') &
                     & real(tper(0,1,1,is)),aimag(tper(0,1,1,is))
             end do
             write(6,'(" phi  = ",2e25.17)') &
                  & real(phi (0,1,1)),aimag(phi (0,1,1))
             write(6,'(" apar = ",2e25.17)') &
                  & real(apar(0,1,1)),aimag(apar(0,1,1))
             write(6,'(" bpar = ",2e25.17)') &
                  & real(bpar(0,1,1)),aimag(bpar(0,1,1))
          endif
       else
          do is=1,nspec
             call transform2(dens(:,:,:,is),densxy(:,:,:,is),nny,nnx)
             call transform2(upar(:,:,:,is),uparxy(:,:,:,is),nny,nnx)
             call transform2(tpar(:,:,:,is),tparxy(:,:,:,is),nny,nnx)
             call transform2(tper(:,:,:,is),tperxy(:,:,:,is),nny,nnx)
          end do

          if(proc0) then
             do is=1,nspec
                write(6,'(" spec = ",i0)') is
                write(6,'(" upar: mode=",2i4," amp=",e25.17)') &
                     & eq_mode_u,0, &
                     & real(upar(0,eq_mode_u+1,1,is))
                write(6,'(" dens: mode=",2i4," amp=",e25.17)') &
                     & eq_mode_n,0, &
                     & real(dens(0,eq_mode_n+1,1,is))
                write(6,'(" tpar: mode=",2i4," amp=",e25.17)') &
                     & eq_mode_n,0, &
                     & real(tpar(0,eq_mode_n+1,1,is))
                write(6,'(" tper: mode=",2i4," amp=",e25.17)') &
                     & eq_mode_n,0, &
                     & real(tper(0,eq_mode_n+1,1,is))
             end do

             call get_unused_unit(unit)
             open(unit,file='check_moment.dat')
             write(unit,'(a1,1x,2a20,4(4a20))') '#','x','y', &
                  & ('dens','upar','tpar','tperp',i=1,nspec)
             do it=1,nx+1
                i=mod(it-1,nx)+1
                do ik=1,ny+1
                   j=mod(ik-1,ny)+1
                   write(unit,'(2x,22d20.12)') lx/nx*(it-1),ly/ny*(ik-1), &
                        & (densxy(i,j,0,is),uparxy(i,j,0,is), &
                        &  tparxy(i,j,0,is),tperxy(i,j,0,is),is=1,nspec)
                end do
                write(unit,*)
             end do
             close(unit)
          endif

          write(6,*)

          call transform2(phi,phixy,nny,nnx)
          call transform2(apar,aparxy,nny,nnx)
          call transform2(bpar,bparxy,nny,nnx)
          do it=1,nakx
             by(:,it,:)=-zi*akx(it)*apar(:,it,:)
          end do
          call transform2(by,byxy,nny,nnx)

          if(proc0) then
             write(6,'(" phi : mode=",2i4," amp=",e25.17)') &
                  & eq_mode_n,0, &
                  & real(phi (0,eq_mode_n+1,1))
             write(6,'(" apar: mode=",2i4," amp=",e25.17)') &
                  & eq_mode_u,0, &
                  & real(apar(0,eq_mode_u+1,1))
             write(6,'(" bpar: mode=",2i4," amp=",e25.17)') &
                  & eq_mode_n,0, &
                  & real(bpar(0,eq_mode_n+1,1))

             call get_unused_unit(unit)
             open(unit,file='check_field.dat')
             write(unit,'(a1,1x,5a20)') '#','x','y','phi','apar','bpar','by'
             do it=1,nx+1
                i=mod(it-1,nx)+1
                do ik=1,ny+1
                   j=mod(ik-1,ny)+1
                   write(unit,'(2x,6d20.12)') lx/nx*(it-1),ly/ny*(ik-1), &
                        & phixy(i,j,0),aparxy(i,j,0),bparxy(i,j,0),byxy(i,j,0)
                end do
                write(unit,*)
             end do
             close(unit)
          endif

       endif
    endif

    deallocate(dens,upar,tpar,tper)
    deallocate(phi,apar,bpar)
    deallocate(jpar)
    if(.not.(nakx==1 .and. naky==1)) then ! grid_option /= single case
       deallocate(densxy,uparxy,tparxy,tperxy)
       deallocate(phixy,aparxy,bparxy)
       deallocate(jparxy)
       deallocate(by,byxy)
    endif

  contains
    function get_xmax(xguess) result(xsol)
      real :: xsol
      real, intent(in) :: xguess
      real, parameter :: tol=1.e-20
      real :: xold
      integer :: i
      integer :: itmax=1000
      xold=xguess
      do i=1,itmax
         xsol=xold-f(xold)/fprime(xold)
         if(abs(xsol-xold) > tol) then
            xold=xsol
         else
            return
         endif
      end do
    end function get_xmax

    function f(x)
      real :: f
      real, intent(in) :: x
      real :: xa, xl1, xl2
      xa=(x-.5*lx)/prof_width
      xl1=2.*pi/lx*x
      xl2=xl1-2.*pi

      f= &
           & 2./prof_width**2*(cosh(2.*xa)-2.)/cosh(xa)**4* & ! f''
           & (tanh(xl1)**2+tanh(xl2)**2-tanh(2.*pi)**2) & ! g
           & + 2. * &
           & (-2./prof_width)*sinh(xa)/cosh(xa)**3* & ! f'
           & 4.*pi/lx*(sinh(xl1)/cosh(xl1)**3+sinh(xl2)/cosh(xl2)**3) & ! g'
           & + &
           & 1./cosh(xa)**2* & ! f
           & 8.*(pi/lx)**2*((2.-cosh(2.*xl1))/cosh(xl1)**4 & ! g''
           &               +(2.-cosh(2.*xl2))/cosh(xl2)**4)

    end function f

    function fprime(x)
      real :: fprime
      real, intent(in) :: x
      real :: xa, xl1, xl2
      xa=(x-.5*lx)/prof_width
      xl1=2.*pi/lx*x
      xl2=xl1-2.*pi

      fprime = &
           & 8./prof_width**3*(2.*sinh(xa)-sinh(xa)**3)/cosh(xa)**5* & ! f''' 
           & (tanh(xl1)**2+tanh(xl2)**2-tanh(2.*pi)**2) & ! g
           & + 3. * &
           & 2./prof_width**2*(cosh(2.*xa)-2.)/cosh(xa)**4* & ! f''
           & 4.*pi/lx*(sinh(xl1)/cosh(xl1)**3+sinh(xl2)/cosh(xl2)**3) & ! g'
           & + 3. * &
           & (-2./prof_width)*sinh(xa)/cosh(xa)**3* & ! f'
           & 8.*(pi/lx)**2*((2.-cosh(2.*xl1))/cosh(xl1)**4 & ! g''
           &               +(2.-cosh(2.*xl2))/cosh(xl2)**4) &
           & + &
           & 1./cosh(xa)**2* & ! f
           & (-64.)*(pi/lx)**3 * ( & ! g'''
           & (2.*sinh(xl1)-sinh(xl1)**3)/cosh(xl1)**5 + &
           & (2.*sinh(xl2)-sinh(xl2)**3)/cosh(xl2)**5 )
    end function fprime

  end subroutine ginit_recon3

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
!    integer :: ig, ik, it, il, is
    integer :: ik, it, il, is

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


	!<doc> This subroutine removes all turbulence except the zonal flow (ky = 0) component upon 
	!restarting. It can be selected by setting the input parameter ginit to "zonal_only". The size of the zonal flows can be adjusted using the input parameter zf_init. Author EGH</doc> 

  subroutine ginit_restart_zonal_only

    use gs2_save, only: gs2_restore
    use mp, only: proc0
    use file_utils, only: error_unit
    use species, only: spec
    use theta_grid, only: ntgrid 
    use kt_grids, only: naky, ntheta0, aky, reality
    use le_grids, only: forbid
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    use fields_arrays, only: phinew, aparnew, bparnew
    use fields_arrays, only: phi, apar, bpar
    use run_parameters, only: fphi, fapar, fbpar
    use ran

    implicit none
!     complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: phi
    integer :: istatus, ierr
    logical :: many = .true.
    real :: a, b
    integer :: iglo
!    integer :: ig, ik, it, il, is
    integer :: ik, it, il, is


! 		if (phiinit > epsilon(0.0)) then 
! 
!     end if 
       

		! <doc> Load phi and g from the restart file
    call gs2_restore (g, scale, istatus, fphi, fapar, fbpar, many)
    if (istatus /= 0) then
       ierr = error_unit()
       if (proc0) write(ierr,*) "Error reading file: ", trim(restart_file)
       g = 0.
    end if
		write (*,*) "Initialised g"

		!<doc> Set all non-zonal components of phi to 0</doc>
    do it = 1, ntheta0
       do ik = 2, naky ! Starting at 2 is the crucial bit!!
          phinew(:,it,ik) = cmplx(0.0,0.0)
       end do
    end do

	! <doc>Allow adjustment of the size of the zonal flows via the input parameter zf_init</doc>
     phinew(:,:,1) = phinew(:,:,1)*zf_init


		!<doc> Apply reality condition for k_theta = 0 component</doc>
!     if (reality) then
!        do it = 1, ntheta0/2
!           phinew(:,it+(ntheta0+1)/2,1) = conjg(phinew(:,(ntheta0+1)/2+1-it,1))
!        enddo
!     end if

		!<doc> Set non-zonal components of g to zero using phi</doc>

		write(*,*) "Zeroing g"

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
				if (ik > 1) then
					g(:,1,iglo) = 0.0 !-phi(:,it,ik)*spec(is)%z*phiinit
					where (forbid(:,il)) g(:,1,iglo) = 0.0
					g(:,2,iglo) = g(:,1,iglo)
				end if
    end do

    ! <doc> If phiinit > 0, add some noise</doc>
    if (phiinit >  epsilon(0.0)) then 
     call ginit_noise
     g = g + gnew
    end if

    gnew = g



  end subroutine ginit_restart_zonal_only





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

  subroutine finish_init_g

    implicit none

    initialized = .false.

  end subroutine finish_init_g

end module init_g
