module legendre

! Taken from Numerical Recipes, incorporated here by Tomo Tatsuno
! Aug, 2005

  implicit none

  public :: nrgauleg

  private

contains

  subroutine nrgauleg (x1, x2, x, w, eps)

    real, intent(in) :: x1, x2
    real, dimension(:), intent(out) :: x, w
    real, intent(in) :: eps

    integer :: its, j, m, n
    integer, parameter :: maxit=100
    double precision :: xl, xm, pi
    double precision, dimension((size(x)+1)/2) :: p1, p2, p3, pp, z, z1
    logical, dimension((size(x)+1)/2) :: unfinished

    n = size(x)
    pi = asin(real(1.0,kind(pi)))*2.0
    m = (n+1)/2

    xm = real(0.5,kind(xm)) * (x1+x2)   ! middle of the section
    xl = real(0.5,kind(xl)) * (x2-x1)   ! signed half length of the section
    z = (/ (cos(pi*(j-0.25)/(n+0.5)), j=1,m) /)
    unfinished = .true.

    do its=1, maxit
       where(unfinished)
          p1 = real(1.0,kind(p1(1)))
          p2 = real(0.0,kind(p2(1)))
       end where
       do j=1, n
          where (unfinished)
             p3 = p2
             p2 = p1
             p1 = ((2*j-1) * z * p2 - (j-1) * p3) / j
          end where
       end do
! p1 now contains the desired legendre polynomials.
       where (unfinished)
          pp = n * (z * p1 - p2) / (z**2 - 1.0)
          z1 = z
          z = z1 - p1 / pp
          unfinished = (abs(z-z1) > eps)
       end where
       if (.not. any(unfinished)) exit
    end do

    if (its == maxit+1) then
       print*, 'too many iterations in nrgauleg'
       stop
    end if
    x(1:m) = xm - xl * z
    x(n:n-m+1:-1) = xm + xl * z
    w(1:m) = 2.0 * abs(xl) / ((1.0 - z**2) * pp**2)
    w(n:n-m+1:-1) = w(1:m)

  end subroutine nrgauleg
  
end module legendre

module egrid

! By Tomo Tatsuno, Aug 2005
! Improved accuracy and speed and maximum number of energy grid points
!

  implicit none

  public :: energy, xgrid, setegrid

  private

  interface xgrid
     module procedure xgrid_s, xgrid_v
  end interface

contains

  subroutine setegrid (Ecut, negrid, epts, wgts)

    use legendre, only: nrgauleg
    implicit none
    
    integer, intent (in) :: negrid
    real, intent (in) :: ecut
    real, dimension(:), intent (out) :: epts, wgts
    integer :: ie, np
    
    real :: eps=1.e-15, x0
    real, dimension(:), allocatable :: zeroes
    
    real :: x
    
    if (Ecut > 20.0) then
! should go to runname.error?
       write (*,*) 'E -> x transformation is not numerically correct for such a large Ecut'
       write (*,*) 'x(E=20) =', xgrid(20.), ' is too close to 1.0'
    end if
    
    np = negrid-1
    allocate (zeroes(np))
    
    x0 = xgrid(ecut)      ! function xgrid_s (single element)
    
    call nrgauleg(0., x0, zeroes, wgts(1:np), eps**1.5)
    
    do ie=1,np
       epts(ie) = energy(zeroes(ie), Ecut)
    end do

    epts(np+1) = ecut
    wgts(np+1) = 1.-x0
    
    if (wgts(negrid) > wgts(np)) then
       write (*,*) 'WARNING: weight at ecut: ', wgts(negrid)
       write (*,*) 'WARNING: is larger than the one at lower grid: ', wgts(np)
       write (*,*) 'WARNING: Recommend using fewer energy grid points'
       if (ecut < 20.0) &
            & write (*,*) 'WARNING: or you should increase ecut (<= 20)'
    end if

    deallocate (zeroes)
    
  end subroutine setegrid

  function energy (xeval, ecut)
    real, intent (in) :: xeval, ecut
    real :: xerrbi, xerrsec, a, b, energy
    integer :: ier

    xerrbi = 1.e-5
    xerrsec = 1.e-13

    a = 0.0
    b = ecut
!    b = 1.0
    call roote (xeval, a, b, xerrbi, xerrsec, 1, ier, energy)

  end function energy

! called as roote (xeval, a, b, xerrbi, xerrsec, 1, ier, energy)

  subroutine roote(fval,a,b,xerrbi,xerrsec,nsolv,ier,soln)

! TT
! solve xgrid(E) == fval for E
!   xerrbi: error max in x for bisection routine
!   xerrsec: error max in x for secant (Newton-Raphson) method
! TT

    integer, intent(in) :: nsolv
    integer, intent(out) :: ier
    real, intent(in) :: fval, a, b, xerrbi, xerrsec
    real, intent(out) :: soln
    integer, parameter :: maxit=30
    integer :: i, niter, isolv
    real :: a1, b1, f1, f2, f3, trial, aold

    ier=0
    a1=a
    b1=b
    f1=xgrid(a1)-fval
    f2=xgrid(b1)-fval

    if (xerrbi > 0.) then
       if (f1*f2 < 0.) then
          niter = int(log(abs(b-a)/xerrbi)/log(2.))+1
          do i=1, niter
             trial = 0.5*(a1+b1)
             f3 = xgrid(trial) - fval
             if (f3*f1 > 0.) then
                a1 = trial
                f1 = f3
             else
                b1 = trial
                f2 = f3
             endif
             !      write(11,*) 'i,a1,f1,b1,f2 ',i,a1,f1,b1,f2
          enddo
       else
          ! TT: This should be connected to '$(run_name).error'?
          write(11,*) 'f1 and f2 have same sign in bisection routine'
          write(11,*) 'a1,f1,b1,f2=',a1,f1,b1,f2
          write(11,*) 'fval=',fval
          ier=1
       endif
    end if

    ! to make (a1,f1) the closest
    if ( abs(f1) > abs(f2) ) then
       f1=f2
       aold=a1
       a1=b1
       b1=aold
    endif

! Newton-Raphson method (formerly it was secant method)
    isolv = 0
    do i=1, maxit
       b1 = a1
       a1 = a1 - f1 / xgrid_prime(a1)
       if (abs(a1-b1) < xerrsec) isolv = isolv+1
       if (isolv >= nsolv) exit
       f1 = xgrid(a1) - fval
    end do

    if (i > maxit) then
! TT: This should be connected to '$(run_name).error'?
       write (11,*) 'bad convergence in roote'
       ier = 1
    end if

    soln = a1

  end subroutine roote

  function xgrid_s (e)
! TT
! this is a function
!              2           E
!   x(E) = ---------- * int   exp(-e) sqrt(e) de
!           sqrt(pi)       0
! which gives energy integral with a Maxwellian weight
! x is a monotonic function of E with
!    E | 0 -> infinity
!   -------------------
!    x | 0 -> 1
! The integral is an error function and numerical evaluation is
! obtained from the formula of incomplete gamma function
! TT

    double precision :: xg, denom, pi
    real :: e, xgrid_s
    integer :: kmax, k, j

    pi = asin(real(1.0,kind(pi))) * 2.0
    kmax = 100
    xg = 0.0

    denom = 1.0
    do k = 0, kmax
       denom = denom * (1.5+k)
       xg = xg + e**(1.5+k) / denom
    end do

    xgrid_s = xg * exp(-e) * 2. / sqrt(pi)

  end function xgrid_s

  function xgrid_prime (e)

! TT
! this is a function
!               2
!   x'(E) = ---------- * exp(-e) sqrt(e)
!            sqrt(pi)
! TT

    real :: e, xgrid_prime, pi

    pi = asin(1.0) * 2.0
    xgrid_prime = exp(-e) * sqrt(e) * 2. / sqrt(pi)

  end function xgrid_prime

  function xgrid_v (e) result (xg)

    real, dimension (:) :: e
    real, dimension (size(e)) :: xg
    real :: denom, pi
    integer :: kmax, k

    pi = asin(1.0) * 2.0

    xg = 0.
    kmax = 100
    denom=1.
    do k = 0, kmax
       denom = denom * (1.5 + k)
       xg = xg + e**(1.5+k) / denom
    end do
    xg = xg * exp(-e) * 2. / sqrt(pi)

  end function xgrid_v

end module egrid

module le_grids
  
  use redistribute, only: redist_type

  implicit none

  public :: init_le_grids
  public :: integrate, lintegrate, integrate_species
  public :: pintegrate, pe_integrate, integrate_stress
  public :: e, anon, al, delal, jend, forbid, dele
  public :: negrid, nlambda, ng2, lmax, integrate_moment
  public :: geint2g, gint2g, orbit_avg
  public :: intcheck
  public :: fcheck

  private

  real, dimension (:,:), allocatable :: e, w, anon, dele ! (negrid,nspec)
  real, dimension (:,:), allocatable :: orbit_avg
  real, dimension (:), allocatable :: al, delal ! (nlambda)
  real, dimension (:,:), allocatable :: wl ! (nlambda,-ntgrid:ntgrid)
  integer, dimension (:), allocatable :: jend ! (-ntgrid:ntgrid)
  logical, dimension (:,:), allocatable :: forbid ! (-ntgrid:ntgrid,nlambda)

  integer :: lint_lo, lint_hi, eint_lo, eint_hi
  integer :: geint2g_lo, geint2g_hi
  complex, dimension (:,:), allocatable :: integration_work
  ! (-ntgrid:ntgrid, -*- processor-dependent -*-)

 ! knobs
  integer :: ngauss, negrid, nesuper, nesub
  real :: ecut, bouncefuzz

  integer :: nlambda, ng2, lmax
  logical :: accel_x = .false.
  logical :: accel_v = .false.
  logical :: test = .false.
  logical :: trapped_particles = .true.
  logical :: advanced_egrid = .true.

  type (redist_type), save :: lambda_map, gint_map, eint_map

contains

  subroutine init_le_grids (accelerated_x, accelerated_v)
    use mp, only: proc0, finish_mp
    use species, only: init_species
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids
    use gs2_layouts, only: init_gs2_layouts
    implicit none
    logical, intent (out) :: accelerated_x, accelerated_v
    logical, save :: initialized = .false.
    integer :: il, ie

    if (initialized) return
    initialized = .true.

    call init_gs2_layouts
    call init_species
    call init_theta_grid
    call init_kt_grids

    if (proc0) then
       call read_parameters
       call set_grids
    end if
    call broadcast_results
    call init_integrations

    accelerated_x = accel_x
    accelerated_v = accel_v

    if (test) then
       if (proc0) then
          do il = 1, nlambda
             write(*,*) al(il)
          end do
          write(*,*) 
          do ie = 1, negrid
             write(*,*) e(ie,1)
          end do
       end if
       call finish_mp
       stop
    endif
    
  end subroutine init_le_grids

  subroutine broadcast_results
    use mp, only: proc0, broadcast
    use species, only: nspec
    use theta_grid, only: ntgrid
    implicit none
    integer :: il, is, ie

    call broadcast (ngauss)
    call broadcast (negrid)
    call broadcast (nesuper)
    call broadcast (nesub)
    call broadcast (ecut)
    call broadcast (bouncefuzz)
    call broadcast (nlambda)
    call broadcast (ng2)
    call broadcast (lmax)
    call broadcast (test)
    call broadcast (trapped_particles)
    call broadcast (advanced_egrid)

    if (.not. proc0) then
       allocate (e(negrid,nspec), w(negrid,nspec), anon(negrid,nspec))
       allocate (dele(negrid,nspec))
       allocate (al(nlambda), delal(nlambda))
       allocate (wl(-ntgrid:ntgrid,nlambda))
       allocate (jend(-ntgrid:ntgrid))
       allocate (forbid(-ntgrid:ntgrid,nlambda))
       allocate (orbit_avg(ng2,negrid))
    end if

    call broadcast (al)
    call broadcast (delal)
    call broadcast (jend)

    do ie=1,negrid
       call broadcast (orbit_avg(:,ie))
    end do

    do is = 1, nspec
       call broadcast (e(:,is))
       call broadcast (w(:,is))
       call broadcast (anon(:,is))
       call broadcast (dele(:,is))
    end do
    do il = 1, nlambda
       call broadcast (wl(:,il))
       call broadcast (forbid(:,il))
    end do
  end subroutine broadcast_results

  subroutine read_parameters
    use species, only: spec, has_slowing_down_species
    use file_utils, only: input_unit, error_unit, input_unit_exist
    implicit none
    integer :: ierr, in_file
    logical :: exist
    namelist /le_grids_knobs/ ngauss, negrid, ecut, bouncefuzz, &
         nesuper, nesub, test, trapped_particles, advanced_egrid

! New default scheme: advanced_egrid is default, and user should set 
! negrid (as in the original code)

    if (advanced_egrid .and. has_slowing_down_species(spec)) then
       ierr = error_unit()
       write (unit=ierr, fmt='("Slowing down species approximated by Maxwellian")')
       write (unit=ierr, fmt='("because advanced_egrid = .true.")')
    end if

    nesub = 8
    nesuper = 2
    ngauss = 5
    negrid = -10
    ecut = 6.0  ! new default value for advanced scheme
    bouncefuzz = 1e-5
    in_file=input_unit_exist("le_grids_knobs", exist)
    if (exist) read (unit=input_unit("le_grids_knobs"), nml=le_grids_knobs)

! user can choose not to set negrid (preferred for old scheme)
    if (negrid == -10) then
       negrid = nesub + nesuper

! New scheme for energy grid, based on GYRO algorithm (Candy and Waltz)
! requires negrid to be set; nesub, nesuper ignored.

       if (advanced_egrid) then
          ierr = error_unit()
          write (unit=ierr, fmt='("Using advanced energy grid with negrid = ",i5)') negrid
       end if

! If user chose negrid, assume nesuper makes sense and check nesub if necessary
    else
       if (.not. advanced_egrid) then
          if (negrid - nesuper /= nesub) then
! Report problem to error file, and continue, using nesuper and negrid
! (Note that nesub is not used anywhere else.)
             nesub = negrid - nesuper
             ierr = error_unit()
             write (unit=ierr, fmt='("Forcing nesub = ",i5)') nesub
          endif
       end if
    endif

  end subroutine read_parameters

  subroutine init_integrations
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use gs2_layouts, only: init_dist_fn_layouts, pe_layout
    use gs2_layouts, only: g_lo, gint_lo, geint_lo 
    use gs2_layouts, only: ik_idx, it_idx, ie_idx, is_idx, idx
    use mp, only: nproc
    implicit none
    integer :: ig, igint, igeint
    integer :: lint_ulim, geint2g_ulim, eint_ulim, ulim
    character (1) :: char
    logical :: first = .true.

    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec)

    if (first) then
       first = .false.
       call pe_layout (char)
       if (char == 'x') then
          accel_x = mod(ntheta0*naky*nspec, nproc) == 0
          accel_v = .false.
       end if
       if (char == 'v') then
          accel_x = .false.
          accel_v = mod(negrid*nlambda*nspec, nproc) == 0
       end if
    end if
          

  end subroutine init_integrations

  subroutine init_slow_integrals
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, gint_lo, geint_lo 
    use gs2_layouts, only: ik_idx, it_idx, ie_idx, is_idx, idx
    implicit none
    integer :: ig, igint, igeint
    integer :: lint_ulim, geint2g_ulim, eint_ulim, ulim
    logical :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_lintegrate 
    call init_eintegrate

    if (accel_x) return

    lint_lo = gint_lo%ulim_world
    lint_hi = gint_lo%llim_world
    geint2g_lo = geint_lo%ulim_world
    geint2g_hi = geint_lo%llim_world
    do ig = g_lo%llim_proc, g_lo%ulim_proc
       igint = idx(gint_lo,ik_idx(g_lo,ig),it_idx(g_lo,ig), &
                   ie_idx(g_lo,ig),is_idx(g_lo,ig))
       lint_lo = min(igint,lint_lo)
       lint_hi = max(igint,lint_hi)
       igeint = idx(geint_lo,ik_idx(gint_lo,igint),it_idx(gint_lo,igint), &
                    is_idx(gint_lo,igint))
       geint2g_lo = min(igeint,geint2g_lo)
       geint2g_hi = max(igeint,geint2g_lo)
    end do
    lint_ulim = lint_hi - lint_lo
    geint2g_ulim = geint2g_hi - geint2g_lo

    eint_lo = geint_lo%ulim_world
    eint_hi = geint_lo%llim_world
    do igint = gint_lo%llim_proc, gint_lo%ulim_proc
       igeint = idx(geint_lo,ik_idx(gint_lo,igint),it_idx(gint_lo,igint), &
                    is_idx(gint_lo,igint))
       eint_lo = min(igeint,eint_lo)
       eint_hi = max(igeint,eint_hi)
    end do
    eint_ulim = eint_hi - eint_lo

    ulim = max(eint_ulim,geint2g_ulim,lint_ulim)
    if (ulim >= 0) then
       allocate (integration_work(-ntgrid:ntgrid,0:ulim))
    end if

  end subroutine init_slow_integrals

  subroutine geint2g (geint, g)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use gs2_layouts, only: g_lo, geint_lo
    use gs2_layouts, only: proc_id, idx_local
    use gs2_layouts, only: ik_idx, it_idx, is_idx, idx
    use mp, only: broadcast, nproc
    implicit none
    complex, dimension (-ntgrid:,geint_lo%llim_proc:), intent (in) :: geint
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g

    complex, dimension (:), allocatable :: collector
    integer :: iglo, igeint, ik, it, is, ig, j
    integer :: igeint_lo, igeint_hi

    ! essentially, spread geint_lo laid out data into g_lo layout,
    ! replicating along the virtual sign,negrid,nlambda dimensions

    ! Often no communication is needed, so do that case separately.

    if (accel_x) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          igeint = iglo/(nlambda*negrid)
          g(:,1,iglo) = geint(:,igeint)
          g(:,2,iglo) = geint(:,igeint)
       end do
    else
       allocate (collector((2*ntgrid+1)*geint_lo%blocksize))
       do igeint_lo = geint_lo%llim_world, geint_lo%ulim_world, geint_lo%blocksize
          igeint_hi = min(igeint_lo + geint_lo%blocksize - 1, geint_lo%ulim_world)
          if (idx_local(geint_lo, igeint_lo)) then
             do igeint = igeint_lo, igeint_hi
                do ig = -ntgrid, ntgrid
                   j = 1 + ig + ntgrid + (2*ntgrid+1)*(igeint-igeint_lo)
                   collector(j) = geint(ig,igeint)
                end do
             end do
          end if
          call broadcast (collector, proc_id(geint_lo, igeint_lo))
          
          do igeint = igeint_lo, igeint_hi
             if (igeint >= geint2g_lo .and. igeint <= geint2g_hi) then
                do ig = -ntgrid, ntgrid
                   j = 1 + ig + ntgrid + (2*ntgrid+1)*(igeint-igeint_lo)
                   integration_work(ig,igeint-geint2g_lo) = collector(j)
                end do
             end if
          end do
       end do
       deallocate (collector)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          igeint = idx(geint_lo,ik,it,is)
          g(:,1,iglo) = integration_work(:,igeint-geint2g_lo)
          g(:,2,iglo) = integration_work(:,igeint-geint2g_lo)
       end do
    end if

  end subroutine geint2g

  subroutine gint2g (gint, g)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use gs2_layouts, only: g_lo, gint_lo
    use gs2_layouts, only: proc_id, idx_local
    use gs2_layouts, only: ik_idx, it_idx, ie_idx, is_idx, idx
    use mp, only: broadcast, nproc
    implicit none
    complex, dimension (-ntgrid:,gint_lo%llim_proc:), intent (in) :: gint
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g

    complex, dimension (:), allocatable :: collector
    integer :: iglo, igint, ik, it, ie, is, ig, j
    integer :: igint_lo, igint_hi

    ! essentially, spread gint_lo laid out data into g_lo layout,
    ! replicating along the virtual sign, nlambda dimensions

    ! Often no communication is needed, so do that case separately.

    if (accel_x) then 
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          igint = iglo/nlambda
          g(:,1,iglo) = gint(:,igint)
          g(:,2,iglo) = gint(:,igint)
       end do       
    else
       allocate (collector((2*ntgrid+1)*gint_lo%blocksize))
       do igint_lo = gint_lo%llim_world, gint_lo%ulim_world, gint_lo%blocksize
          igint_hi = min(igint_lo + gint_lo%blocksize - 1, gint_lo%ulim_world)
          if (idx_local (gint_lo, igint_lo)) then          
             do igint = igint_lo, igint_hi
                do ig = -ntgrid, ntgrid
                   j = 1 + ig + ntgrid + (2*ntgrid+1)*(igint-igint_lo)
                   collector(j) = gint(ig,igint)
                end do
             end do
          end if
          call broadcast (collector, proc_id(gint_lo,igint_lo))
          
          do igint = igint_lo, igint_hi
             if (igint >= lint_lo .and. igint <= lint_hi) then
                do ig = -ntgrid, ntgrid
                   j = 1 + ig + ntgrid + (2*ntgrid+1)*(igint-igint_lo)
                   integration_work(ig,igint-lint_lo) = collector(j)
                end do
             end if
          end do
       end do
       deallocate (collector)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          igint = idx(gint_lo,ik,it,ie,is)
          g(:,1,iglo) = integration_work(:,igint-lint_lo)
          g(:,2,iglo) = integration_work(:,igint-lint_lo)
       end do
    end if

  end subroutine gint2g

  subroutine integrate (g1, geint)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, gint_lo, geint_lo
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g1
    complex, dimension (-ntgrid:,geint_lo%llim_proc:), intent (out) :: geint
    complex, dimension (:,:), allocatable :: gint

    call init_slow_integrals

    allocate (gint (-ntgrid:ntgrid, gint_lo%llim_proc:gint_lo%ulim_alloc))
    call lintegrate (g1, gint)
    call eintegrate (gint, geint)

    deallocate (gint)

  end subroutine integrate

  subroutine sum_species (geint, weights, total)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: geint_lo, ik_idx, it_idx, is_idx
    use mp, only: sum_allreduce
    implicit none
    complex, dimension (-ntgrid:,geint_lo%llim_proc:), intent (in) :: geint
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:), intent (out) :: total
    complex, dimension (:), allocatable :: work

    integer :: i, ig, ik, it, is
    integer :: igeint

    total = 0.0
    do igeint = geint_lo%llim_proc, geint_lo%ulim_proc
       ik = ik_idx(geint_lo,igeint)
       it = it_idx(geint_lo,igeint)
       is = is_idx(geint_lo,igeint)
       total(:,it,ik) = total(:,it,ik) + weights(is)*geint(:,igeint)
    end do

    allocate (work((2*ntgrid+1)*naky*ntheta0)) ; work = 0.
    i = 0
    do ik = 1, naky
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid
             i = i + 1
             work(i) = total(ig, it, ik)
          end do
       end do
    end do
    
    call sum_allreduce (work) 

    i = 0
    do ik = 1, naky
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid
             i = i + 1
             total(ig, it, ik) = work(i)
          end do
       end do
    end do
    deallocate (work)

  end subroutine sum_species

  subroutine geint20 (geint, total)
    use theta_grid, only: ntgrid
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: geint_lo, ik_idx, it_idx, is_idx
    use mp, only: sum_reduce, proc0
    implicit none
    complex, dimension (-ntgrid:,geint_lo%llim_proc:), intent (in) :: geint
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
! work could be smaller on PEs other than 0.  I will do this later.
    complex, dimension (:), allocatable :: work  

    integer :: i, ig, ik, it, is
    integer :: igeint

    total = 0.0
    do igeint = geint_lo%llim_proc, geint_lo%ulim_proc
       ik = ik_idx(geint_lo,igeint)
       it = it_idx(geint_lo,igeint)
       is = is_idx(geint_lo,igeint)
       total(:,it,ik,is) = total(:,it,ik,is) + geint(:,igeint)
    end do

! reduce number of calls to sum_reduce 

    allocate (work((2*ntgrid+1)*naky*ntheta0*nspec)); work = 0.
    i = 0
    do is = 1, nspec
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntgrid, ntgrid
                i = i + 1
                work(i) = total(ig, it, ik, is)
             end do
          end do
       end do
    end do

    call sum_reduce (work, 0)

    if (proc0) then
       i = 0
       do is = 1, nspec
          do ik = 1, naky
             do it = 1, ntheta0
                do ig = -ntgrid, ntgrid
                   i = i + 1
                   total(ig, it, ik, is) = work(i)
                end do
             end do
          end do
       end do
    end if
    deallocate (work)

  end subroutine geint20

  subroutine integrate_species (g, weights, total)
    use theta_grid, only: ntgrid
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, idx, idx_local
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_allreduce
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:), intent (out) :: total
! total = total(theta, kx, ky)
    complex, dimension (:,:), allocatable :: geint
    complex, dimension (:), allocatable :: work
    real :: fac
    integer :: is, il, ie, ik, it, iglo, ig, i

!    total = 0.
!    do is = 1, nspec
!       do ie = 1, negrid
!          fac = weights(is)*w(ie,is)
!          do ik = 1, naky
!             do it = 1, ntheta0
!                do il = 1, nlambda
!                   iglo = idx (g_lo, ik, it, il, ie, is)
!                   if (.not. idx_local (g_lo, iglo)) cycle
!                   do ig = -ntgrid, ntgrid
!!                      total(ig, it, ik) = total(ig, it, ik) + &
!!                              weights(is)*w(ie,is)*wl(ig,il)*(g(ig,1,iglo)+g(ig,2,iglo))
!                      total(ig, it, ik) = total(ig, it, ik) + &
!                              fac*wl(ig,il)*(g(ig,1,iglo)+g(ig,2,iglo))
!                   end do
!                end do
!             end do
!          end do
!       end do
!    end do

    total = 0.
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       fac = weights(is)*w(ie,is)

       total(:, it, ik) = total(:, it, ik) + fac*wl(:,il)*(g(:,1,iglo)+g(:,2,iglo))
    end do

    allocate (work((2*ntgrid+1)*naky*ntheta0)) ; work = 0.
    i = 0
    do ik = 1, naky
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid
             i = i + 1
             work(i) = total(ig, it, ik)
          end do
       end do
    end do
    
    call sum_allreduce (work) 

    i = 0
    do ik = 1, naky
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid
             i = i + 1
             total(ig, it, ik) = work(i)
          end do
       end do
    end do
    deallocate (work)

  end subroutine integrate_species

  subroutine integrate_stress (g, stress)

    use theta_grid, only: ntgrid, delthet, grho, bmag, gradpar
    use species, only: nspec
    use kt_grids, only: naky, ntheta0, aky
    use gs2_layouts, only: g_lo, idx, idx_local
    use mp, only: sum_reduce, proc0
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:), intent (out) :: stress
    complex, dimension (:), allocatable :: work
    real, dimension (-ntgrid:ntgrid) :: delnorm
    real :: fac, wgt
    integer :: is, ie, ik, it, il, iglo, ig, i

    delnorm=delthet*grho/bmag/gradpar
    delnorm=delnorm/sum(delnorm)
    
    wgt=0.
    do ig=-ntgrid,ntgrid
       wgt=wgt+delnorm(ig)
    end do
    delnorm=delnorm/wgt

    stress = 0.
    do is = 1, nspec
       do ie = 1, negrid
          fac = w(ie,is)
          do ik = 1, naky
             if (aky(ik) > epsilon(0.0)) cycle
             do it = 1, ntheta0
                do il = 1, nlambda
                   iglo = idx (g_lo, ik, it, il, ie, is)
                   if (idx_local (g_lo, iglo)) then
                      stress(it, is) = stress(it, is) + &
                           sum(fac*wl(:,il)*(g(:,1,iglo)+g(:,2,iglo))*delnorm(:))
                   end if
                end do
             end do
          end do
       end do
    end do

    allocate (work(ntheta0*nspec)) ; work = 0.
    i = 0
    do is=1, nspec
       do it = 1, ntheta0
          i = i + 1
          work(i) = stress(it, is)
       end do
    end do

    call sum_reduce (work, 0)

    if (proc0) then
       i = 0
       do is = 1, nspec
          do it = 1, ntheta0
             i = i + 1
             stress(it, is) = work(i)
          end do
       end do
    end if
    deallocate (work)

  end subroutine integrate_stress

  subroutine integrate_moment (g, total)
! returns moments to PE 0
    use mp, only: nproc
    use theta_grid, only: ntgrid
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, idx, idx_local
    use mp, only: sum_reduce, proc0
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total

    complex, dimension (:), allocatable :: work
    real :: fac
    integer :: is, il, ie, ik, it, iglo, ig, i

    total = 0.
    do is = 1, nspec
       do ie = 1, negrid
          fac = w(ie,is)
          do il = 1, nlambda
             do it = 1, ntheta0
                do ik = 1, naky
                   iglo = idx (g_lo, ik, it, il, ie, is)
                   if (idx_local (g_lo, iglo)) then
                      do ig = -ntgrid, ntgrid
                         total(ig, it, ik, is) = total(ig, it, ik, is) + &
                              fac*wl(ig,il)*(g(ig,1,iglo)+g(ig,2,iglo))
                      end do
                   end if
                end do
             end do
          end do
       end do
    end do

    if (nproc > 1) then
       allocate (work((2*ntgrid+1)*naky*ntheta0*nspec)) ; work = 0.
       i = 0
       do is = 1, nspec
          do ik = 1, naky
             do it = 1, ntheta0
                do ig = -ntgrid, ntgrid
                   i = i + 1
                   work(i) = total(ig, it, ik, is)
                end do
             end do
          end do
       end do
       
       call sum_reduce (work, 0)
       
       if (proc0) then
          i = 0
          do is = 1, nspec
             do ik = 1, naky
                do it = 1, ntheta0
                   do ig = -ntgrid, ntgrid
                      i = i + 1
                      total(ig, it, ik, is) = work(i)
                   end do
                end do
             end do
          end do
       end if
       deallocate (work)
    end if

  end subroutine integrate_moment

  subroutine set_grids
    use species, only: init_species, nspec
    use theta_grid, only: init_theta_grid, ntgrid, nbset, bset, eps
    implicit none

    call init_theta_grid
    call init_species

    allocate (e(negrid,nspec), w(negrid,nspec), anon(negrid,nspec))
    allocate (dele(negrid,nspec))
    call egridset

    dele(1,:) = e(1,:)
    dele(2:,:) = e(2:,:)-e(:negrid-1,:)

    ng2 = 2*ngauss
!    if (eps > epsilon(0.0)) then
    if (trapped_particles .and. eps > epsilon(0.0)) then
       nlambda = ng2+nbset
       lmax = nlambda-1
    else
       nlambda = ng2
       lmax = nlambda
    end if
    allocate (al(nlambda), delal(nlambda))
    allocate (wl(-ntgrid:ntgrid,nlambda))
    allocate (jend(-ntgrid:ntgrid))
    allocate (forbid(-ntgrid:ntgrid,nlambda))
    allocate (orbit_avg(ng2,negrid))
    if (nlambda-ng2 > 0) then
       al(ng2+1:nlambda) = 1.0/bset
    end if
    call lgridset
    delal(1) = al(1)
    delal(2:) = al(2:) - al(:nlambda-1)
  end subroutine set_grids

  subroutine egridset
    use species, only: nspec, spec, slowing_down_species
    use legendre, only: nrgauleg
    use constants
    use egrid, only: setegrid
    implicit none

    real, dimension (1), parameter :: esub1 = (/ &
         0.5000000000E+00 /)
    real, dimension (2), parameter :: esub2 = (/ &
         0.2113248654E+00, &
         0.7886751346E+00 /)
    real, dimension (3), parameter :: esub3 = (/ &
         0.1127016654E+00, &
         0.5000000000E+00, &
         0.8872983346E+00 /)
    real, dimension (4), parameter :: esub4 = (/ &
         0.6943184420E-01, &
         0.3300094782E+00, &
         0.6699905218E+00, &
         0.9305681558E+00 /)
    real, dimension (5), parameter :: esub5 = (/ &
         0.4691007703E-01, &
         0.2307653449E+00, &
         0.5000000000E+00, &
         0.7692346551E+00, &
         0.9530899230E+00 /)
    real, dimension (6), parameter :: esub6 = (/ &
         0.3376524290E-01, &
         0.1693953068E+00, &
         0.3806904070E+00, &
         0.6193095930E+00, &
         0.8306046932E+00, &
         0.9662347571E+00 /)
! A&S
    real, dimension (7), parameter :: esub7 = (/ &
         0.0254460438, &
         0.1292344072, &
         0.2970774243, &
         0.5000000000, & 
         0.7029225757, &
         0.8707655928, &
         0.9745539562 /)
    real, dimension (8), parameter :: esub8 = (/ &
         .1985507175E-01, &
         .1016667613E+00, &
         .2372337950E+00, &
         .4082826788E+00, &
         .5917173212E+00, &
         .7627662050E+00, &
         .8983332387E+00, &
         .9801449282E+00 /)
! CRC
    real, dimension (9), parameter :: esub9 = (/ &
         0.01592, &
         0.08199, &
         0.19332, &
         0.33788, & 
         0.50000, &
         0.66213, &
         0.80669, &
         0.91802, &
         0.98408 /)
    real, dimension (10), parameter :: esub10 = (/ &
         0.1304673574E-01, &
         0.6746831666E-01, &
         0.1602952159E+00, &
         0.2833023029E+00, &
         0.4255628305E+00, &
         0.5744371695E+00, &
         0.7166976971E+00, &
         0.8397047841E+00, &
         0.9325316833E+00, &
         0.9869532643E+00 /)
    real, dimension (12), parameter :: esub12 = (/ &
         0.9219682877E-02, &
         0.4794137181E-01, &
         0.1150486629E+00, &
         0.2063410229E+00, &
         0.3160842505E+00, &
         0.4373832957E+00, &
         0.5626167043E+00, &
         0.6839157495E+00, &
         0.7936589771E+00, &
         0.8849513371E+00, &
         0.9520586282E+00, &
         0.9907803171E+00 /)
    real, dimension (14), parameter :: esub14 = (/ &
         0.6858095652E-02, &
         0.3578255817E-01, &
         0.8639934247E-01, &
         0.1563535476E+00, &
         0.2423756818E+00, &
         0.3404438155E+00, &
         0.4459725256E+00, &
         0.5540274744E+00, &
         0.6595561845E+00, &
         0.7576243182E+00, &
         0.8436464524E+00, &
         0.9136006575E+00, &
         0.9642174418E+00, &
         0.9931419043E+00 /)
    real, dimension (16), parameter :: esub16 = (/ &
         0.5299532504E-02, &
         0.2771248846E-01, &
         0.6718439881E-01, &
         0.1222977958E+00, &
         0.1910618778E+00, &
         0.2709916112E+00, &
         0.3591982246E+00, &
         0.4524937451E+00, &
         0.5475062549E+00, &
         0.6408017754E+00, &
         0.7290083888E+00, &
         0.8089381222E+00, &
         0.8777022042E+00, &
         0.9328156012E+00, &
         0.9722875115E+00, &
         0.9947004675E+00 /)
    real, dimension (20), parameter :: esub20 = (/ &
         0.3435700407E-02, &
         0.1801403636E-01, &
         0.4388278587E-01, &
         0.8044151409E-01, &
         0.1268340468E+00, &
         0.1819731596E+00, &
         0.2445664990E+00, &
         0.3131469556E+00, &
         0.3861070744E+00, &
         0.4617367394E+00, &
         0.5382632606E+00, &
         0.6138929256E+00, &
         0.6868530444E+00, &
         0.7554335010E+00, &
         0.8180268404E+00, &
         0.8731659532E+00, &
         0.9195584859E+00, &
         0.9561172141E+00, &
         0.9819859636E+00, &
         0.9965642996E+00 /)
    real, dimension (24), parameter :: esub24 = (/ &
         0.2406390001E-02, &
         0.1263572201E-01, &
         0.3086272400E-01, &
         0.5679223650E-01, &
         0.8999900701E-01, &
         0.1299379042E+00, &
         0.1759531740E+00, &
         0.2272892643E+00, &
         0.2831032462E+00, &
         0.3424786602E+00, &
         0.4044405663E+00, &
         0.4679715536E+00, &
         0.5320284464E+00, &
         0.5955594337E+00, &
         0.6575213398E+00, &
         0.7168967538E+00, &
         0.7727107357E+00, &
         0.8240468260E+00, &
         0.8700620958E+00, &
         0.9100009930E+00, &
         0.9432077635E+00, &
         0.9691372760E+00, &
         0.9873642780E+00, &
         0.9975936100E+00 /)
    real, dimension (32), parameter :: esub32 = (/ &
         0.1368069075E-02, 0.7194244227E-02, 0.1761887221E-01, 0.3254696203E-01, &
         0.5183942212E-01, 0.7531619313E-01, 0.1027581020E+00, 0.1339089406E+00, &
         0.1684778665E+00, 0.2061421214E+00, 0.2465500455E+00, 0.2893243619E+00, &
         0.3340656989E+00, 0.3803563189E+00, 0.4277640192E+00, 0.4758461672E+00, &
         0.5241538328E+00, 0.5722359808E+00, 0.6196436811E+00, 0.6659343011E+00, &
         0.7106756381E+00, 0.7534499545E+00, 0.7938578786E+00, 0.8315221335E+00, &
         0.8660910594E+00, 0.8972418980E+00, 0.9246838069E+00, 0.9481605779E+00, &
         0.9674530380E+00, 0.9823811278E+00, 0.9928057558E+00, 0.9986319309E+00 /)
    real, dimension (48), parameter :: esub48 = (/ &
         0.6144963738E-03, 0.3234913867E-02, 0.7937708139E-02, 0.1470420373E-01, &
         0.2350614842E-01, 0.3430665465E-01, 0.4706043164E-01, 0.6171398986E-01, &
         0.7820586919E-01, 0.9646689799E-01, 0.1164204837E+00, 0.1379829345E+00, &
         0.1610638102E+00, 0.1855663016E+00, 0.2113876370E+00, 0.2384195126E+00, &
         0.2665485476E+00, 0.2956567590E+00, 0.3256220569E+00, 0.3563187563E+00, &
         0.3876181048E+00, 0.4193888220E+00, 0.4514976504E+00, 0.4838099145E+00, &
         0.5161900855E+00, 0.5485023496E+00, 0.5806111780E+00, 0.6123818952E+00, &
         0.6436812437E+00, 0.6743779431E+00, 0.7043432410E+00, 0.7334514524E+00, &
         0.7615804874E+00, 0.7886123630E+00, 0.8144336984E+00, 0.8389361898E+00, &
         0.8620170655E+00, 0.8835795163E+00, 0.9035331020E+00, 0.9217941308E+00, &
         0.9382860101E+00, 0.9529395684E+00, 0.9656933454E+00, 0.9764938516E+00, &
         0.9852957963E+00, 0.9920622919E+00, 0.9967650861E+00, 0.9993855036E+00 /)
    real, dimension (64), parameter :: esub64 = (/ &
         0.3474791321E-03, 0.1829941614E-02, 0.4493314262E-02, 0.8331873058E-02, 0.1333658611E-01, &
         0.1949560017E-01, 0.2679431257E-01, 0.3521541393E-01, 0.4473893146E-01, 0.5534227700E-01, &
         0.6700030092E-01, 0.7968535187E-01, 0.9336734244E-01, 0.1080138205E+00, 0.1235900464E+00, &
         0.1400590749E+00, 0.1573818435E+00, 0.1755172644E+00, 0.1944223224E+00, 0.2140521769E+00, &
         0.2343602680E+00, 0.2552984271E+00, 0.2768169914E+00, 0.2988649210E+00, 0.3213899208E+00, &
         0.3443385640E+00, 0.3676564189E+00, 0.3912881781E+00, 0.4151777898E+00, 0.4392685904E+00, &
         0.4635034391E+00, 0.4878248537E+00, 0.5121751463E+00, 0.5364965609E+00, 0.5607314096E+00, &
         0.5848222102E+00, 0.6087118219E+00, 0.6323435811E+00, 0.6556614360E+00, 0.6786100792E+00, &
         0.7011350790E+00, 0.7231830086E+00, 0.7447015729E+00, 0.7656397320E+00, 0.7859478231E+00, &
         0.8055776776E+00, 0.8244827356E+00, 0.8426181565E+00, 0.8599409251E+00, 0.8764099536E+00, &
         0.8919861795E+00, 0.9066326576E+00, 0.9203146481E+00, 0.9329996991E+00, 0.9446577230E+00, &
         0.9552610685E+00, 0.9647845861E+00, 0.9732056874E+00, 0.9805043998E+00, 0.9866634139E+00, &
         0.9916681269E+00, 0.9955066857E+00, 0.9981700584E+00, 0.9996525209E+00 /)
    real, dimension (1), parameter :: wgt1 = (/ &
         0.1000000000E+01 /)
    real, dimension (2), parameter :: wgt2 = (/ &
         0.5000000000E+00, &
         0.5000000000E+00 /)
    real, dimension (3), parameter :: wgt3 = (/ &
         0.2777777778E+00, &
         0.4444444444E+00, &
         0.2777777778E+00 /)
    real, dimension (4), parameter :: wgt4 = (/ &
         0.1739274226E+00, &
         0.3260725774E+00, &
         0.3260725774E+00, &
         0.1739274226E+00 /)
    real, dimension (5), parameter :: wgt5 = (/ &
         0.1184634425E+00, &
         0.2393143352E+00, &
         0.2844444444E+00, &
         0.2393143352E+00, &
         0.1184634425E+00 /)
    real, dimension (6), parameter :: wgt6 = (/ &
         0.8566224619E-01, &
         0.1803807865E+00, &
         0.2339569673E+00, &
         0.2339569673E+00, &
         0.1803807865E+00, &
         0.8566224619E-01 /)
    real, dimension (7), parameter :: wgt7 = (/ &
         0.0647424831, &
         0.1398526957, &
         0.1909150253, &
         0.2089795918, &
         0.1909150253, &
         0.1398526957, &
         0.0647424831 /)
    real, dimension (8), parameter :: wgt8 = (/ &
         0.5061426815E-01, &
         0.1111905172E+00, &
         0.1568533229E+00, &
         0.1813418917E+00, &
         0.1813418917E+00, &
         0.1568533229E+00, &
         0.1111905172E+00, &
         0.5061426815E-01 /)
    real, dimension (9), parameter :: wgt9 = (/ &
         0.04064, &
         0.09033, &
         0.13031, &
         0.15618, &
         0.16512, &
         0.15618, &
         0.13031, &
         0.09033, &
         0.04064 /)
    real, dimension (10), parameter :: wgt10 = (/ &
         0.3333567215E-01, &
         0.7472567458E-01, &
         0.1095431813E+00, &
         0.1346333597E+00, &
         0.1477621124E+00, &
         0.1477621124E+00, &
         0.1346333597E+00, &
         0.1095431813E+00, &
         0.7472567458E-01, &
         0.3333567215E-01 /)
    real, dimension (12), parameter :: wgt12 = (/ &
         0.2358766819E-01, &
         0.5346966300E-01, &
         0.8003916427E-01, &
         0.1015837134E+00, &
         0.1167462683E+00, &
         0.1245735229E+00, &
         0.1245735229E+00, &
         0.1167462683E+00, &
         0.1015837134E+00, &
         0.8003916427E-01, &
         0.5346966300E-01, &
         0.2358766819E-01 /)
    real, dimension (14), parameter :: wgt14 = (/ &
         0.1755973017E-01, &
         0.4007904358E-01, &
         0.6075928534E-01, &
         0.7860158358E-01, &
         0.9276919874E-01, &
         0.1025992319E+00, &
         0.1076319267E+00, &
         0.1076319267E+00, &
         0.1025992319E+00, &
         0.9276919874E-01, &
         0.7860158358E-01, &
         0.6075928534E-01, &
         0.4007904358E-01, &
         0.1755973017E-01 /)
    real, dimension (16), parameter :: wgt16 = (/ &
         0.1357622971E-01, &
         0.3112676197E-01, &
         0.4757925584E-01, &
         0.6231448563E-01, &
         0.7479799441E-01, &
         0.8457825970E-01, &
         0.9130170752E-01, &
         0.9472530523E-01, &
         0.9472530523E-01, &
         0.9130170752E-01, &
         0.8457825970E-01, &
         0.7479799441E-01, &
         0.6231448563E-01, &
         0.4757925584E-01, &
         0.3112676197E-01, &
         0.1357622971E-01 /)
    real, dimension (20), parameter :: wgt20 = (/ &
         0.8807003570E-02, &
         0.2030071490E-01, &
         0.3133602417E-01, &
         0.4163837079E-01, &
         0.5096505991E-01, &
         0.5909726598E-01, &
         0.6584431922E-01, &
         0.7104805466E-01, &
         0.7458649324E-01, &
         0.7637669357E-01, &
         0.7637669357E-01, &
         0.7458649324E-01, &
         0.7104805466E-01, &
         0.6584431922E-01, &
         0.5909726598E-01, &
         0.5096505991E-01, &
         0.4163837079E-01, &
         0.3133602417E-01, &
         0.2030071490E-01, &
         0.8807003570E-02 /)
    real, dimension (24), parameter :: wgt24 = (/ &
         0.6170614900E-02, &
         0.1426569431E-01, &
         0.2213871941E-01, &
         0.2964929246E-01, &
         0.3667324071E-01, &
         0.4309508077E-01, &
         0.4880932605E-01, &
         0.5372213506E-01, &
         0.5775283403E-01, &
         0.6083523646E-01, &
         0.6291872817E-01, &
         0.6396909767E-01, &
         0.6396909767E-01, &
         0.6291872817E-01, &
         0.6083523646E-01, &
         0.5775283403E-01, &
         0.5372213506E-01, &
         0.4880932605E-01, &
         0.4309508077E-01, &
         0.3667324071E-01, &
         0.2964929246E-01, &
         0.2213871941E-01, &
         0.1426569431E-01, &
         0.6170614900E-02 /)
    real, dimension (32), parameter :: wgt32 = (/ &
         0.3509305005E-02, 0.8137197365E-02, 0.1269603265E-01, 0.1713693146E-01, 0.2141794901E-01, &
         0.2549902963E-01, 0.2934204674E-01, 0.3291111139E-01, 0.3617289705E-01, 0.3909694789E-01, &
         0.4165596211E-01, 0.4382604650E-01, 0.4558693935E-01, 0.4692219954E-01, 0.4781936004E-01, &
         0.4827004426E-01, 0.4827004426E-01, 0.4781936004E-01, 0.4692219954E-01, 0.4558693935E-01, &
         0.4382604650E-01, 0.4165596211E-01, 0.3909694789E-01, 0.3617289705E-01, 0.3291111139E-01, &
         0.2934204674E-01, 0.2549902963E-01, 0.2141794901E-01, 0.1713693146E-01, 0.1269603265E-01, &
         0.8137197365E-02, 0.3509305005E-02 /)
    real, dimension (48), parameter :: wgt48 = (/ &
         0.1576673026E-02, 0.3663776951E-02, 0.5738617290E-02, 0.7789657861E-02, 0.9808080229E-02, &
         0.1178538042E-01, 0.1371325485E-01, 0.1558361392E-01, 0.1738861128E-01, 0.1912067553E-01, &
         0.2077254147E-01, 0.2233728043E-01, 0.2380832925E-01, 0.2517951778E-01, 0.2644509474E-01, &
         0.2759975185E-01, 0.2863864605E-01, 0.2955741985E-01, 0.3035221958E-01, 0.3101971158E-01, &
         0.3155709614E-01, 0.3196211929E-01, 0.3223308222E-01, 0.3236884841E-01, 0.3236884841E-01, &
         0.3223308222E-01, 0.3196211929E-01, 0.3155709614E-01, 0.3101971158E-01, 0.3035221958E-01, &
         0.2955741985E-01, 0.2863864605E-01, 0.2759975185E-01, 0.2644509474E-01, 0.2517951778E-01, &
         0.2380832925E-01, 0.2233728043E-01, 0.2077254147E-01, 0.1912067553E-01, 0.1738861128E-01, &
         0.1558361392E-01, 0.1371325485E-01, 0.1178538042E-01, 0.9808080229E-02, 0.7789657861E-02, &
         0.5738617290E-02, 0.3663776951E-02, 0.1576673026E-02 /)
    real, dimension (64), parameter :: wgt64 = (/ &
         0.8916403608E-03, 0.2073516630E-02, 0.3252228984E-02, 0.4423379913E-02, 0.5584069730E-02, &
         0.6731523948E-02, 0.7863015238E-02, 0.8975857888E-02, 0.1006741158E-01, 0.1113508690E-01, &
         0.1217635128E-01, 0.1318873486E-01, 0.1416983631E-01, 0.1511732854E-01, 0.1602896418E-01, &
         0.1690258092E-01, 0.1773610663E-01, 0.1852756427E-01, 0.1927507659E-01, 0.1997687057E-01, &
         0.2063128162E-01, 0.2123675756E-01, 0.2179186226E-01, 0.2229527908E-01, 0.2274581396E-01, &
         0.2314239829E-01, 0.2348409141E-01, 0.2377008286E-01, 0.2399969430E-01, 0.2417238112E-01, &
         0.2428773372E-01, 0.2434547850E-01, 0.2434547850E-01, 0.2428773372E-01, 0.2417238112E-01, &
         0.2399969430E-01, 0.2377008286E-01, 0.2348409141E-01, 0.2314239829E-01, 0.2274581396E-01, &
         0.2229527908E-01, 0.2179186226E-01, 0.2123675756E-01, 0.2063128162E-01, 0.1997687057E-01, &
         0.1927507659E-01, 0.1852756427E-01, 0.1773610663E-01, 0.1690258092E-01, 0.1602896418E-01, &
         0.1511732854E-01, 0.1416983631E-01, 0.1318873486E-01, 0.1217635128E-01, 0.1113508690E-01, &
         0.1006741158E-01, 0.8975857888E-02, 0.7863015238E-02, 0.6731523948E-02, 0.5584069730E-02, &
         0.4423379913E-02, 0.3252228984E-02, 0.2073516630E-02, 0.8916403608E-03 /)
    real, dimension (2), parameter :: xsup2 = (/ &
         0.5857864376, &
         3.4142135624 /)
    real, dimension (3), parameter :: xsup3 = (/ &
         0.4157745568, &
         2.2942803603, &
         6.2899450829 /)
    real, dimension (4), parameter :: xsup4 = (/ &
         0.3225476896, &
         1.7457611012, &
         4.5366202969, &
         9.3950709123 /)
    real, dimension (5), parameter :: xsup5 = (/ &
         0.2635603197, &
         1.4134030591, &
         3.5964257710, &
         7.0858100059, &
         12.6408008443 /)
    real, dimension (6), parameter :: xsup6 = (/ &
	 0.2228466042, &
         1.1889321017, &
	 2.9927363261, &
	 5.7751435691, &
	 9.8374674184, &
         15.9828739806 /)
    real, dimension (7), parameter :: xsup7 = (/ &
	 0.1930436766, &
	 1.0266648953, &
	 2.5678767450, &
	 4.9003530845, &
	 8.1821534446, &
	 12.7341802918, &
	 19.3957278623 /)
    real, dimension (8), parameter :: xsup8 = (/ &
	 0.1702796323, &
	 0.9037017768, &
	 2.2510866299, &
	 4.2667001703, &
	 7.0459054023, &
	 10.7585160102, &
	 15.7406786413, &
	 22.8631317369 /)
    real, dimension (9), parameter :: xsup9 = (/ &
	 0.1523222277, & 
	 0.8072200227, & 
	 2.0051351556, & 
	 3.7834739733, & 
	 6.2049567779, & 
	 9.3729852517, & 
	 13.4662369111, & 
	 18.8335977890, &
	 26.3740718909 /)
    real, dimension (10), parameter :: xsup10 = (/ &
	 0.1377934705, &
	 0.7294545495, &
	 1.8083429017, &
	 3.4014336979, &
	 5.5524961401, &
	 8.3301527468, &
	 11.8437858379, &
	 16.2792578314, &
	 21.9965858120, &
	 29.9206970123 /)
    real, dimension (12), parameter :: xsup12 = (/ &
	 0.1157221174, &
	 0.6117574845, &
	 1.5126102698, &
	 2.8337513377, &
	 4.5992276394, &
	 6.8445254531, &
	 9.6213168425, &
	 13.0060549933, &
	 17.1168551875, &
	 22.1510903794, &
	 28.4879672510, &
	 37.0991210445 /)
    real, dimension (15), parameter :: xsup15 = (/ &
	 0.0933078120, &
	 0.4926917403, &
	 1.2155954121, &
	 2.2699495262, &
	 3.6676227218, &
	 5.4253366274, &
	 7.5659162266, &
	 10.1202285680, &
	 13.1302824822, &
	 16.6544077083, &
	 20.7764788994, &
	 25.6238942267, &
	 31.4075191698, &
	 38.5306833065, &
	 48.0260855727 /) 

    real, dimension (2), parameter :: wsup2 = (/ &
         8.53553390593e-1, &
         1.46446609407e-1 /)
    real, dimension (3), parameter :: wsup3 = (/ &
         7.11093009929e-1, &
         2.78517733569e-1, &
         1.03892565016e-2 /)
    real, dimension (4), parameter :: wsup4 = (/ &
         6.03154104342e-1, &
         3.57418692438e-1, &
         3.88879085150e-2, &
         5.39294705561e-4 /)
    real, dimension (5), parameter :: wsup5 = (/ &
         5.21755610583e-1, &
         3.98666811083e-1, &
         7.59424496817e-2, &
	 3.61175867992e-3, &
	 2.33699723858e-5 /)
    real, dimension (6), parameter :: wsup6 = (/ &
         4.58964673950e-1, & 
	 4.17000830772e-1, & 
 	 1.13373382074e-1, &
	 1.03991974531e-2, &
	 2.61017202815e-4, &
	 8.98547906430e-7 /)		
    real, dimension (7), parameter :: wsup7 = (/ &
	 4.09318951701e-1, &
	 4.21831277862e-1, &
	 1.47126348658e-1, &
	 2.06335144687e-2, &
	 1.07401014328e-3, &
	 1.58654643486e-5, &
	 3.17031547900e-8 /)
    real, dimension (8), parameter :: wsup8 = (/ &
	 3.69188589342e-1, &
	 4.18786780814e-1, &
 	 1.75794986637e-1, & 
	 3.33434922612e-2, &
	 2.79453623523e-3, &
	 9.07650877336e-5, &
	 8.48574671627e-7, &
	 1.04800117487e-9 /)
    real, dimension (9), parameter :: wsup9 = (/ &
	 3.36126421798e-1, &
	 4.11213980424e-1, &
	 1.99287525371e-1, &
	 4.74605627657e-2, &
	 5.59962661079e-3, &
	 3.05249767093e-4, &
	 6.59212302608e-6, &
	 4.11076933035e-8, &
	 3.29087403035e-11 /)
    real, dimension (10), parameter :: wsup10 = (/ &
	 3.08441115765e-1, &
	 4.01119929155e-1, &
	 2.18068287612e-1, &
	 6.20874560987e-2, &
	 9.50151697518e-3, &
	 7.53008388588e-4, &
	 2.82592334960e-5, &
	 4.24931398496e-7, &
	 1.83956482398e-9, &
	 9.91182721961e-13 /)
    real, dimension (12), parameter :: wsup12 = (/ &
	 2.64731371055e-1, &
	 3.77759275873e-1, &
	 2.44082011320e-1, &
	 9.04492222117e-2, &
	 2.01023811546e-2, &
	 2.66397354187e-3, &
	 2.03231592663e-4, &
	 8.36505585682e-6, &
	 1.66849387654e-7, &
	 1.34239103052e-9, &
	 3.06160163504e-12, &
	 8.14807746743e-16 /)
    real, dimension (15), parameter :: wsup15 = (/ &
	 2.18234885940e-1, &
	 3.42210177923e-1, &
	 2.63027577942e-1, &
	 1.26425818106e-1, &
	 4.02068649210e-2, &
	 8.56387780361e-3, &
	 1.21243614721e-3, &
	 1.11674392344e-4, &
	 6.45992676202e-6, &
	 2.22631690710e-7, &
	 4.22743038498e-9, &
	 3.92189726704e-11, &
	 1.45651526407e-13, &
	 1.48302705111e-16, &
	 1.60059490621e-20 /)	

    real, dimension (nesub) :: esub, wsub
    real, dimension (nesuper) :: xsup, wsup
    integer :: is, isup
    integer :: ng1
    real :: cut
    real :: eps=1.e-15

    if (.not. advanced_egrid) then

!       call nrgauleg (0., 1.0, esub, wsub, eps**1.5)

       select case (nesub)
       case (1)  
          esub = esub1
          wsub = wgt2
       case (2)
          esub = esub2
          wsub = wgt2
       case (3)
          esub = esub3
          wsub = wgt3
       case (4)
          esub = esub4
          wsub = wgt4
       case (5)
          esub = esub5
          wsub = wgt5
       case (6)
          esub = esub6
          wsub = wgt6
       case (7)
          esub = esub7
          wsub = wgt7
       case (8)
          esub = esub8
          wsub = wgt8
       case (9)
          esub = esub9
          wsub = wgt9
       case (10)
          esub = esub10
          wsub = wgt10
       case (12)
          esub = esub12
          wsub = wgt12
       case (14)
          esub = esub14
          wsub = wgt14
       case (16)
          esub = esub16
          wsub = wgt16
       case (20)
          esub = esub20
          wsub = wgt20
       case (24)
          esub = esub24
          wsub = wgt24
       case (32)
          esub = esub32
          wsub = wgt32
       case (48)
          esub = esub48
          wsub = wgt48
       case (64)
          esub = esub64
          wsub = wgt64
       case default
          call stop_invalid ("nesub", nesub)
       end select
       
       select case (nesuper)
       case (1)  ! only for debugging; not correct
          xsup = xsup2(1)
          wsup = wsup2(1)
       case (2)
          xsup = xsup2
          wsup = wsup2
       case (3)
          xsup = xsup3
          wsup = wsup3
       case (4)
          xsup = xsup4
          wsup = wsup4
       case (5)
          xsup = xsup5
          wsup = wsup5
       case (6)
          xsup = xsup6
          wsup = wsup6
       case (7)
          xsup = xsup7
          wsup = wsup7
       case (8)
          xsup = xsup8
          wsup = wsup8
       case (9)
          xsup = xsup9
          wsup = wsup9
       case (10)
          xsup = xsup10
          wsup = wsup10
       case (12)
          xsup = xsup12
          wsup = wsup12
       case (15)
          xsup = xsup15
          wsup = wsup15
       case default
          call stop_invalid ("nesuper", nesuper)
       end select

       if (negrid < 1) call stop_invalid ("negrid", negrid)

       do is = 1, nspec
          if (spec(is)%type == slowing_down_species) then
             e(:negrid-1,is) = esub
             w(:negrid-1,is) = wsub*esub**2
             anon(:negrid-1,is) = -1.5/e(:negrid-1,is)
             e(negrid,is) = 1.0
             w(negrid,is) = 1e-6
             anon(negrid,is) = 1e6
             w(:,is) = w(:,is)*0.75
          else
             ng1 = max(negrid-nesuper,0)
             select case (ng1)
             case (0)  
                cut = 0.0
! these values not correct if nesuper = 1; included for debugging only
                do isup = 1, nesuper
                   e(isup,is) = cut + xsup(isup)
                   w(isup,is) = wsup(isup)*exp(-cut)*sqrt(e(isup,is))
                end do
             case default
                cut = ecut
                do isup = 1, nesuper
                   e(ng1+isup,is) = cut + xsup(isup)
                   w(ng1+isup,is) = wsup(isup)*exp(-cut)*sqrt(e(ng1+isup,is))
                end do
                e(:ng1,is) = cut*esub**(2.0/3.0)
                w(:ng1,is) = (2.0/3.0)*exp(-e(:ng1,is))*wsub*cut**1.5             
             end select
             w(:,is) = w(:,is)*0.5/sqrt(pi)
             anon(:,is) = 1.0
          end if
       end do
       
    else

       call setegrid (ecut, negrid, e(:,1), w(:,1))
       do is = 2, nspec
          e(:,is) = e(:,1)
          w(:,is) = w(:,1)
       end do
       w = 0.25*w
       anon = 1.0

    end if

  end subroutine egridset

  subroutine lgridset

! Modified to use nrgauleg routine, Tomo Tatsuno, Aug 2005

    use theta_grid, only: ntgrid, bmag, bmax, eps
    use legendre, only: nrgauleg
    use constants
    implicit none

! note that xgauss and wgauss are transposed wrt original code

    real :: tiny = 1.e-15
    real, dimension (2*ngauss) :: wx, xx
    real :: ww
    integer :: ig, il

    call nrgauleg(1., 0., xx, wx, tiny**1.5)

    wl = 0.0
    
    al(:ng2) = (1.0 - xx(:ng2)**2)/bmax

    do il = 1, ng2
       do ig = -ntgrid, ntgrid
          wl(ig,il) = wx(il)*2.0*sqrt((bmag(ig)/bmax) &
               *((1.0/bmax-al(il))/(1.0/bmag(ig)-al(il))))
       end do
    end do

    jend = 0
    forbid = .false.

!    if (eps <= epsilon(0.0)) return
    if (.not. trapped_particles .or. eps <= epsilon(0.0)) return

    jend = ng2 + 1
    do il = ng2+1, nlambda-1
       do ig = -ntgrid, ntgrid
          if (1.0-al(il)*bmag(ig) > -bouncefuzz &
               .and. 1.0-al(il+1)*bmag(ig) > -bouncefuzz) &
          then
             jend(ig) = jend(ig) + 1
             ww = sqrt(max(1.0 -   al(il)*bmag(ig),0.0)) - &
                  sqrt(max(1.0 - al(il+1)*bmag(ig),0.0))
             wl(ig,il)   = wl(ig,il)   + ww
             wl(ig,il+1) = wl(ig,il+1) + ww
          end if
       end do
    end do

    do il = ng2+1, nlambda
       do ig = -ntgrid, ntgrid
          forbid(ig,il) = 1.0 - al(il)*bmag(ig) < -bouncefuzz
       end do
    end do

    call init_orbit_average

  end subroutine lgridset

  subroutine init_orbit_average

! not sure about kxfac in the following

    use theta_grid, only: drhodpsi, delthet, ntgrid, bmag
    use constants, only: pi
    implicit none
    real :: vp
    integer :: il, ig, ie

! presently assuming all species are Maxwellian

    orbit_avg = 0.
    do ie=1,negrid
       do il=1,ng2
          do ig=-ntgrid,ntgrid
             vp = sqrt(e(ie,1)*max(0.0,1.0-al(il)*bmag(ig)))
             orbit_avg (il,ie) = orbit_avg(il,ie)+delthet(ig)/vp
          end do
       end do
    end do
    
    orbit_avg = 2.*pi*drhodpsi/orbit_avg

!    write (*,*) '# drhorpsi= ',drhodpsi
!    do ie=1,negrid
!       do il=1,ng2
!          write (*,*) il, al(il), ie, e(ie,1), orbit_avg(il, ie)
!       end do
!    end do

  end subroutine init_orbit_average

  subroutine stop_invalid (name, val)
    use file_utils, only: error_unit
    use mp, only: proc0, finish_mp
    implicit none
    character(*), intent (in) :: name
    integer, intent (in) :: val
    integer :: ierr

    if (proc0) then
       ierr = error_unit()
       write (unit=ierr, fmt='("Invalid value for ",a,": ",i5)') name, val
    end if
    call finish_mp
    stop
  end subroutine stop_invalid

  subroutine intcheck (g0)
    use file_utils, only: open_output_file, close_output_file
    use species, only: nspec 
    use theta_grid, only: ntgrid, bmag, bmax
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, gint_lo, idx_local, proc_id
    use gs2_layouts, only: idx, il_idx, ie_idx 
    use mp, only: proc0, send, receive
    use constants
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g0
    integer :: unit
    integer :: iglo, igint, ig, ik, it, il, ie, is
    complex, dimension (:,:), allocatable :: gint
    complex, dimension (negrid) :: dumout

    call open_output_file (unit, ".intcheck")

    write (unit,*) "bmax ", bmax

    write (unit,*) "untrapped check:"
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       if (il <= 2*ngauss) then
          g0(:,1,iglo) = bmax/bmag*(1.0 + zi)*sqrt(1.0 - al(il)*bmag) &
               *(1.0 - bmax*al(il))**(0.5*real(ie-2))
       else
          g0(:,1,iglo) = 0.0
       end if
       g0(:,2,iglo) = g0(:,1,iglo)
    end do
    allocate (gint(-ntgrid:ntgrid,gint_lo%llim_proc:gint_lo%ulim_alloc))
    call lintegrate (g0, gint)
    do is = 1, nspec
       do it = 1, ntheta0
          do ik = 1, naky
             do ig = -ntgrid, ntgrid
                do ie = 1, negrid
                   igint = idx(gint_lo,ik,it,ie,is)
                   if (proc0) then
                      if (idx_local(gint_lo,igint)) then
                         dumout(ie) = gint(ig,igint)
                      else
                         call receive (dumout(ie), proc_id(gint_lo,igint))
                      end if
                   else if (idx_local(gint_lo,igint)) then
                      call send (gint(ig,igint), 0)
                   end if
                end do
                write (unit,*) "is,j,fac ", &
                     is, ig, max(0.0,1.0-bmag(ig)/bmag)**0.5
                write (unit,"(20(1x,1pe12.5))") 0.25*real(ie)*dumout
             end do
          end do
       end do
    end do

    write (unit,*)
    write (unit,*)
    write (unit,*) "trapped check:"
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       if (il <= 2*ngauss) then
          g0(:,1,iglo) = 0.0
       else
          g0(:,1,iglo) = (max(0.0,1.0-al(il)*bmag)**0.5)**real(ie-1)
       end if
       g0(:,2,iglo) = g0(:,1,iglo)
    end do
    call lintegrate (g0, gint)
    do is = 1, nspec
       do it = 1, ntheta0
          do ik = 1, naky
             do ig = -ntgrid, ntgrid
                do ie = 1, negrid
                   igint = idx(gint_lo,ik,it,ie,is)
                   if (proc0) then
                      if (idx_local(gint_lo,igint)) then
                         dumout(ie) = gint(ig,igint)
                      else
                         call receive (dumout(ie), proc_id(gint_lo,igint))
                      end if
                   else if (idx_local(gint_lo,igint)) then
                      call send (gint(ig,igint), 0)
                   end if
                end do
                write (unit,*) "is,j,fac ", &
                     is, ig, max(0.0,1.0-bmag(ig)/bmag)**0.5
                write (unit,"(20(1x,1pe12.5))") &
                     0.5*real(ie)*dumout &
                     /(max(0.0,1.0-bmag(ig)/bmag)**0.5)**real(ie)
             end do
          end do
       end do
    end do
    deallocate (gint)
    
    do ig = -ntgrid, ntgrid
       write (unit,*) "j ", ig
       do il = 1, nlambda
          write (unit,*) "il,wl,sq ",il,wl(ig,il), &
               sqrt(max(0.0,1.0-al(il)*bmag(ig)))
       end do
    end do

    call close_output_file (unit)
  end subroutine intcheck

  subroutine fcheck (g, f)
    use species, only: nspec
    use theta_grid, only: ntgrid, gradpar, bmag, delthet
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use mp, only: sum_allreduce
    use constants
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,:,:), intent (out) :: f
    integer :: iglo, ik, it, il, ie, is

    f = 0.0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       if (il > 1) then
          f(il,it,ik,is) = f(il,it,ik,is) &
               + pi/(al(il)-al(il-1))*w(ie,is) &
                 *sum(delthet(:ntgrid-1)/gradpar(:ntgrid-1)/bmag(:ntgrid-1) &
                      *(g(:ntgrid-1,1,iglo) + g(:ntgrid-1,2,iglo)) &
                      *sqrt(max(0.0,1.0-al(il)*bmag(:ntgrid-1)))) &
                 /sum(delthet(:ntgrid-1)/gradpar(:ntgrid-1)/bmag(:ntgrid-1))
       end if
       if (il < nlambda) then
          f(il+1,it,ik,is) = f(il+1,it,ik,is) &
               - pi/(al(il+1)-al(il))*w(ie,is) &
                 *sum(delthet(:ntgrid-1)/gradpar(:ntgrid-1)/bmag(:ntgrid-1) &
                      *(g(:ntgrid-1,1,iglo) + g(:ntgrid-1,2,iglo)) &
                      *sqrt(max(0.0,1.0-al(il)*bmag(:ntgrid-1)))) &
                 /sum(delthet(:ntgrid-1)/gradpar(:ntgrid-1)/bmag(:ntgrid-1))
       end if
    end do

    do is = 1, nspec 
       do ik = 1, naky
          do it = 1, ntheta0
             call sum_allreduce (f(:,it,ik,is))
          end do
       end do
    end do
  end subroutine fcheck

  subroutine pintegrate (g, fld, result)
    
    use theta_grid, only: ntgrid, grho, bmag, gradpar, kxfac
    use species, only: nspec
    use kt_grids, only: naky, ntheta0, aky
    use gs2_layouts, only: g_lo, idx, ik_idx, it_idx, is_idx, ie_idx, il_idx
    use mp, only: sum_reduce

    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (nlambda, nspec) :: result
    real :: wgt, fac
    integer :: iglo, is, il, ie, ik, it, ig 

    result = 0.
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       wgt = 1./sum(grho/bmag/gradpar)
       fac = 0.5
       if (aky(ik) == 0.) fac = 1.0
       do ig = -ntgrid, ntgrid
          result(il,is) = result(il,is) + aimag( &  
               w(ie,is)*wl(ig,il)* &  
               (g(ig,1,iglo)+g(ig,2,iglo))*conjg(fld(ig,it,ik))*aky(ik)* &
               grho(ig)/bmag(ig)/gradpar(ig)*wgt*0.5*kxfac*fac)
       end do
    end do
    
    do is = 1,nspec
       result(:,is) = result(:,is)/delal(:)
       call sum_reduce (result(:,is), 0)
    end do
    
  end subroutine pintegrate

  subroutine pe_integrate (g, fld, result)
    
    use theta_grid, only: ntgrid, grho, bmag, gradpar, kxfac
    use species, only: nspec
    use kt_grids, only: naky, ntheta0, aky
    use gs2_layouts, only: g_lo, idx, ik_idx, it_idx, is_idx, ie_idx, il_idx
    use mp, only: sum_reduce, proc0

    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (negrid, nspec) :: result
    real, dimension (:), allocatable :: work
    real :: wgt, fac
    integer :: iglo, is, il, ie, ik, it, ig, i

    result = 0.
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       wgt = 1./sum(grho/bmag/gradpar)
       fac = 0.5
       if (aky(ik) == 0.) fac = 1.0
       do ig = -ntgrid, ntgrid
          result(ie,is) = result(ie,is) + aimag( &  
               w(ie,is)*wl(ig,il)* &  
               (g(ig,1,iglo)+g(ig,2,iglo))*conjg(fld(ig,it,ik))*aky(ik)* &
               grho(ig)/bmag(ig)/gradpar(ig)*wgt*0.5*kxfac*fac)
       end do
    end do
    
    result = result/dele
    
    allocate (work(negrid*nspec)) ; work = 0.
    i = 0
    do is = 1, nspec
       do ie = 1, negrid
          i = i + 1
          work(i) = result(ie,is)
       end do
    end do

    call sum_reduce (work, 0)
    
    if (proc0) then
       i = 0
       do is = 1, nspec
          do ie = 1, negrid
             i = i + 1
             result(ie,is) = work(i) 
          end do
       end do

    end if
    deallocate (work)
    
  end subroutine pe_integrate

  subroutine init_lintegrate
    use species, only: nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use mp, only: nproc 
!    use gs2_layouts, only: init_lambda_layouts, lambda_lo
!    use gs2_layouts, only: gidx2lamidx, lamidx2gintidx
    use gs2_layouts, only: g_lo, gint_lo, gidx2gintidx
    use gs2_layouts, only: idx_local, proc_id, init_dist_fn_layouts
    use redistribute, only: index_list_type, init_fill, delete_list
    implicit none
    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension (3) :: from_low, from_high
    integer, dimension (4) :: to_low, to_high
!    integer, dimension (2) :: to_gint_low, from_lambda_low, to_lhigh, from_lhigh
    integer :: iglo, isign, ig, ip, n
    integer :: ilam, il, igint
    logical :: done = .false.

    if (done) return
    done = .true.

    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec)
         
    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       call gidx2gintidx (g_lo, iglo, gint_lo, il, ilam)
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
             if (idx_local(g_lo,iglo)) &
                  nn_from(proc_id(gint_lo,ilam)) = nn_from(proc_id(gint_lo,ilam)) + 1
             if (idx_local(gint_lo,ilam)) &
                  nn_to(proc_id(g_lo,iglo)) = nn_to(proc_id(g_lo,iglo)) + 1
          end do
       end do
    end do

    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first(nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
          allocate (from_list(ip)%third(nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first(nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
          allocate (to_list(ip)%third(nn_to(ip)))
          allocate (to_list(ip)%fourth(nn_to(ip)))
       end if
    end do

    ! get local indices of elements distributed to/from other processors
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       call gidx2gintidx (g_lo, iglo, gint_lo, il, ilam)
       if (idx_local(g_lo,iglo)) then
          ip = proc_id(gint_lo,ilam)
          do isign = 1, 2
             do ig = -ntgrid, ntgrid
                n = nn_from(ip) + 1
                nn_from(ip) = n
                from_list(ip)%first(n) = ig
                from_list(ip)%second(n) = isign
                from_list(ip)%third(n) = iglo
             end do
          end do
       end if
       if (idx_local(gint_lo,ilam)) then
          ip = proc_id(g_lo,iglo)
          do isign = 1, 2
             do ig = -ntgrid, ntgrid
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to_list(ip)%first(n) = ig
                to_list(ip)%second(n) = isign
                to_list(ip)%third(n) = il
                to_list(ip)%fourth(n) = ilam
             end do
          end do
       end if
    end do

    from_low (1) = -ntgrid
    from_low (2) = 1
    from_low (3) = g_lo%llim_proc

    from_high (1) = ntgrid
    from_high (2) = 2
    from_high (3) = g_lo%ulim_alloc

    to_low (1) = -ntgrid
    to_low (2) = 1
    to_low (3) = 1
    to_low (4) = gint_lo%llim_proc

    to_high (1) = ntgrid
    to_high (2) = 2
    to_high (3) = max(2*nlambda, 2*ng2+1)
    to_high (4) = gint_lo%ulim_alloc

    call init_fill (lambda_map, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)

! There is old coding that allows a lambda layout different from the gint layout.
! This was never used much, but the coding appears after the end of the module

  end subroutine init_lintegrate

  subroutine lintegrate (g1, gint)
    use theta_grid, only: ntgrid
    use redistribute, only: gather, fill
!    use gs2_layouts, only: lambda_lo
    use gs2_layouts, only: g_lo, gint_lo 
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g1
    complex, dimension (-ntgrid:,gint_lo%llim_proc:), intent (out) :: gint

    complex, dimension (:,:,:,:), allocatable :: glam
    complex, dimension (:,:), allocatable :: work

    integer :: igint, il, ig
    logical :: first = .true.

    if (first) call init_slow_integrals
    first = .false.

    allocate (glam (-ntgrid:ntgrid, 2, max(2*nlambda,2*ng2+1), &
         gint_lo%llim_proc:gint_lo%ulim_alloc))
    allocate (work (-ntgrid:ntgrid, gint_lo%llim_proc:gint_lo%ulim_alloc))

    call gather (lambda_map, g1, glam)

    work = 0.
    do igint = gint_lo%llim_proc, gint_lo%ulim_proc
       do il = 1, nlambda
          do ig = -ntgrid, ntgrid
             work(ig,igint) = work(ig,igint) &
                  + wl(ig,il)*(glam(ig,1,il,igint)+ glam(ig,2,il,igint))
          end do
       end do
    end do

    gint = work    
    deallocate (glam, work)
!
! fill statement only needed if lambda_lo and gint_lo diverge in future.  
! (see init_lintegrate; coding moved after end of module)
!
!    if (.false.) then
!       call fill (gint_map, work, gint)
!    endif

  end subroutine lintegrate

  subroutine init_eintegrate
    use theta_grid, only: ntgrid
    use mp, only: nproc 
    use gs2_layouts, only: gint_lo, geint_lo, gintidx2geidx
    use gs2_layouts, only: idx_local, proc_id
    use redistribute, only: index_list_type, init_fill, delete_list
    implicit none
    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension (2) :: from_low, from_high
    integer, dimension (3) :: to_low, to_high
    integer :: ig, ip, n, ie
    integer :: igint, igeint
    logical :: done = .false.

    if (done) return
    done = .true.

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do igint = gint_lo%llim_world, gint_lo%ulim_world
       call gintidx2geidx(gint_lo, igint, ie, geint_lo, igeint)
       if (idx_local(gint_lo,igint)) &
            nn_from(proc_id(geint_lo,igeint)) = nn_from(proc_id(geint_lo,igeint)) + 2*ntgrid + 1
       if (idx_local(geint_lo,igeint)) &
            nn_to(proc_id(gint_lo,igint)) = nn_to(proc_id(gint_lo,igint)) + 2*ntgrid + 1
    end do

    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first(nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first(nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
          allocate (to_list(ip)%third(nn_to(ip)))
       end if
    end do

    ! get local indices of elements distributed to/from other processors
    nn_to = 0
    nn_from = 0
    do igint = gint_lo%llim_world, gint_lo%ulim_world
       call gintidx2geidx(gint_lo, igint, ie, geint_lo, igeint)
       if (idx_local(gint_lo,igint)) then
          ip = proc_id(geint_lo,igeint)
          do ig = -ntgrid, ntgrid
             n = nn_from(ip) + 1
             nn_from(ip) = n
             from_list(ip)%first(n) = ig
             from_list(ip)%second(n) = igint
          end do
       end if
       if (idx_local(geint_lo,igeint)) then
          ip = proc_id(gint_lo,igint)
          do ig = -ntgrid, ntgrid
             n = nn_to(ip) + 1
             nn_to(ip) = n
             to_list(ip)%first(n) = ig
             to_list(ip)%second(n) = ie
             to_list(ip)%third(n) = igeint
          end do
       end if
    end do

    from_low (1) = -ntgrid
    from_low (2) = gint_lo%llim_proc

    to_low (1) = -ntgrid
    to_low (2) = 1
    to_low (3) = geint_lo%llim_proc

    to_high (1) = ntgrid
    to_high (2) = negrid
    to_high (3) = geint_lo%ulim_alloc
 
    from_high (1) = ntgrid
    from_high (2) = gint_lo%ulim_alloc

    call init_fill (eint_map, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)

  end subroutine init_eintegrate

  subroutine eintegrate (gint, geint)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: gint_lo, geint_lo, is_idx
    use redistribute, only: gather
    implicit none
    complex, dimension (-ntgrid:,gint_lo%llim_proc:), intent (in) :: gint
    complex, dimension (-ntgrid:,geint_lo%llim_proc:), intent (out) :: geint

    complex, dimension (:,:,:), allocatable :: work

    integer :: igeint, ie, is

    allocate (work(-ntgrid:ntgrid, negrid, geint_lo%llim_proc:geint_lo%ulim_alloc))
    work = 0. ; geint = 0.

    call gather (eint_map, gint, work)

    do igeint = geint_lo%llim_proc, geint_lo%ulim_proc
       is = is_idx(geint_lo, igeint)
       do ie = 1, negrid
          geint(:,igeint) = geint(:,igeint) + w(ie,is)*work(:,ie,igeint)
       end do
    end do

    deallocate (work)

  end subroutine eintegrate

end module le_grids


 
