module egrid

! By Tomo Tatsuno, Aug 2005
! Improved accuracy and speed and maximum number of energy grid points
!

  implicit none

  public :: setvgrid, init_egrid
  public :: zeroes, x0

  private

  real :: x0
  real, dimension(:), allocatable, save :: zeroes

contains

  subroutine init_egrid (negrid)
    
    integer, intent (in) :: negrid

    if (.not. allocated(zeroes)) then
       allocate (zeroes(negrid-1)) ; zeroes = 0.
    end if

  end subroutine init_egrid

  subroutine setvgrid (vcut, negrid, epts, wgts, nesub)

    use constants, only: pi => dpi
    use gauss_quad, only: get_legendre_grids_from_cheb, get_laguerre_grids

    implicit none
    
    integer, intent (in) :: negrid
    real, intent (in) :: vcut
    integer, intent (in) :: nesub
    real, dimension(:), intent (out) :: epts, wgts

    call init_egrid (negrid)
    
    ! get grid points in v up to vcut (epts is not E yet)
    call get_legendre_grids_from_cheb (0., vcut, epts(:nesub), wgts(:nesub))

    ! change from v to E
    epts(:nesub) = epts(:nesub)**2

    ! absorb exponential and volume element in weights
    wgts(:nesub) = wgts(:nesub)*epts(:nesub)*exp(-epts(:nesub))/sqrt(pi)

    if (negrid > nesub) then

       ! get grid points in y = E - vcut**2 (epts not E yet)
       call get_laguerre_grids (epts(nesub+1:), wgts(nesub+1:))

       ! change from y to E
       epts(nesub+1:) = epts(nesub+1:) + vcut**2

       ! absort exponential and volume element in weights
       wgts(nesub+1:) = wgts(nesub+1:)*0.5*sqrt(epts(nesub+1:)/pi)*exp(-vcut**2)

    end if

    zeroes = sqrt(epts(:negrid-1))
    x0 = sqrt(epts(negrid))

  end subroutine setvgrid

end module egrid

module le_grids
  
  use redistribute, only: redist_type

  implicit none

  public :: init_le_grids, finish_le_grids
  public :: read_parameters, wnml_le_grids
  public :: integrate_species
  public :: energy, anon, al, delal, jend, forbid, dele, wl, w
  public :: negrid, nlambda, ng2, lmax, integrate_moment, nesub
  public :: xloc
  public :: xx, nterp, testfac, new_trap_int, vcut
  public :: init_weights, legendre_transform, lagrange_interp, lagrange_coefs
  public :: eint_error, lint_error, trap_error, integrate_test, wdim
  public :: integrate_kysum, integrate_volume ! MAB
  public :: get_flux_vs_theta_vs_vpa

  private

  interface integrate_moment
     module procedure integrate_moment_c34
     module procedure integrate_moment_lec
     module procedure integrate_moment_r33
  end interface

  interface integrate_volume
     module procedure integrate_volume_c
     module procedure integrate_volume_r
  end interface

  interface get_hermite_polynomials
     module procedure get_hermite_polynomials_1d
     module procedure get_hermite_polynomials_4d
  end interface

  real, dimension (:), allocatable :: xx ! (ng2)
  real, dimension (:,:), allocatable, save :: werr, wlerr, xloc ! mbmark
  real, dimension (:,:,:), allocatable, save :: wlterr
  real, dimension (:,:,:,:), allocatable, save :: lgrnge
  real, dimension (:,:,:), allocatable, save :: wlmod
  real, dimension (:,:,:), allocatable, save :: wtmod
  real, dimension (:,:), allocatable, save :: wmod
  real, dimension (:,:), allocatable, save :: lpe, lpl
  real, dimension (:,:,:,:), allocatable, save :: lpt

  real, dimension (:), allocatable :: energy, w, anon, dele ! (negrid,nspec)
  real, dimension (:), allocatable :: al, delal ! (nlambda)
  real, dimension (:,:), allocatable :: wl ! (nlambda,-ntgrid:ntgrid)
  integer, dimension (:), allocatable :: jend ! (-ntgrid:ntgrid)
  logical, dimension (:,:), allocatable :: forbid ! (-ntgrid:ntgrid,nlambda)

  integer :: wdim
  integer :: lint_lo, lint_hi, eint_lo, eint_hi
  integer :: geint2g_lo, geint2g_hi
  complex, dimension (:,:), allocatable :: integration_work
  ! (-ntgrid:ntgrid, -*- processor-dependent -*-)
  logical :: exist

 ! knobs
  integer :: ngauss, negrid, nesuper, nesub
  real :: bouncefuzz, vcut

  integer :: nlambda, ng2, lmax
  logical :: accel_x = .false.
  logical :: accel_v = .false.
  logical :: test = .false.
  logical :: trapped_particles = .true.
  logical :: new_trap_int = .false.

  logical :: intinit = .false.
  logical :: slintinit = .false.
  logical :: lintinit = .false.
  logical :: eintinit = .false.
  logical :: initialized = .false.

  integer :: testfac = 1
  integer :: nmax = 500
  integer :: nterp = 100

  real :: wgt_fac = 10.0

contains

  subroutine wnml_le_grids(unit)
    implicit none
    integer :: unit
    if (.not. exist) return
       write (unit, *)
       write (unit, fmt="(' &',a)") "le_grids_knobs"
       write (unit, fmt="(' nesub = ',i4)") nesub
       write (unit, fmt="(' nesuper = ',i4)") nesuper
       write (unit, fmt="(' ngauss = ',i4)") ngauss
       write (unit, fmt="(' vcut = ',e16.10)") vcut
       write (unit, fmt="(' /')")
  end subroutine wnml_le_grids

  subroutine init_le_grids (accelerated_x, accelerated_v)
    use mp, only: proc0, finish_mp
    use species, only: init_species
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids
    use gs2_layouts, only: init_gs2_layouts
    implicit none
    logical, intent (out) :: accelerated_x, accelerated_v
!    logical, save :: initialized = .false.
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
             write(*,*) energy(ie)
          end do
       end if
       call finish_mp
       stop
    endif
    
  end subroutine init_le_grids

  subroutine broadcast_results
    use mp, only: proc0, broadcast
    use species, only: nspec
    use egrid, only: zeroes, x0, init_egrid
    use theta_grid, only: ntgrid
    implicit none
    integer :: il, is, ie, ipt, isgn, tsize

    tsize = 2*nterp-1

    call broadcast (ngauss)
    call broadcast (negrid)
    call broadcast (nesuper)
    call broadcast (nesub)
    call broadcast (vcut)
    call broadcast (bouncefuzz)
    call broadcast (nlambda)
    call broadcast (ng2)
    call broadcast (lmax)
    call broadcast (test)
    call broadcast (testfac)
    call broadcast (nmax)
    call broadcast (trapped_particles)
    call broadcast (wgt_fac)
    call broadcast (new_trap_int)
    call broadcast (nterp)

    if (.not. proc0) then
       allocate (energy(negrid), w(negrid), anon(negrid))
       allocate (dele(negrid))
       allocate (al(nlambda), delal(nlambda))
       allocate (wl(-ntgrid:ntgrid,nlambda))
       allocate (jend(-ntgrid:ntgrid))
       allocate (forbid(-ntgrid:ntgrid,nlambda))
       allocate (xx(ng2))
       allocate (lgrnge(-ntgrid:ntgrid,nlambda,tsize,2))
       allocate (xloc(-ntgrid:ntgrid,tsize))
    end if

    call init_egrid (negrid)
    call broadcast (xx)
    call broadcast (x0)
    call broadcast (zeroes)

    call broadcast (al)
    call broadcast (delal)
    call broadcast (jend)
 
    call broadcast (energy)
    call broadcast (dele)
    call broadcast (w)
    call broadcast (anon)

    do il = 1, nlambda
       call broadcast (wl(:,il))
       call broadcast (forbid(:,il))
    end do

    do ipt = 1, tsize
       call broadcast (xloc(:,ipt))
       do isgn=1,2
          do il=1, nlambda
             call broadcast (lgrnge(:,il,ipt,isgn))
          end do
       end do
    end do

  end subroutine broadcast_results

  subroutine read_parameters
    use species, only: spec, has_slowing_down_species
    use file_utils, only: input_unit, error_unit, input_unit_exist
    implicit none
    integer :: ierr, in_file
    namelist /le_grids_knobs/ ngauss, negrid, bouncefuzz, &
         nesuper, nesub, test, trapped_particles, &
         testfac, nmax, wgt_fac, new_trap_int, nterp, vcut

    nesub = 8
    nesuper = 2
    ngauss = 5
    negrid = -10
    vcut = 4.0
    bouncefuzz = 1e-5
    in_file=input_unit_exist("le_grids_knobs", exist)
!    if (exist) read (unit=input_unit("le_grids_knobs"), nml=le_grids_knobs)
    if (exist) read (unit=in_file, nml=le_grids_knobs)

! user can choose not to set negrid (preferred for old scheme)
    if (negrid == -10) then
       negrid = nesub + nesuper

! If user chose negrid, assume nesuper makes sense and check nesub if necessary
    else 
       nesuper = min(negrid/16+1, 3)
       nesub = negrid - nesuper
    endif

  end subroutine read_parameters

  subroutine init_integrations
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use gs2_layouts, only: init_dist_fn_layouts, pe_layout
    use mp, only: nproc
    implicit none
    character (1) :: char

    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec)

    if (.not. intinit) then
       intinit = .true.
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

  subroutine init_weights

    use file_utils, only: open_output_file, close_output_file
    use egrid, only: zeroes
    use constants, only: pi => dpi

    implicit none

    real, dimension (:), allocatable :: modzeroes, werrtmp  ! (negrid-2)
    real, dimension (:), allocatable :: lmodzeroes, wlerrtmp ! (ng2-1)
    integer :: ipt, ndiv, divmax
    logical :: eflag = .false.

    integer :: ie, il

    allocate(lmodzeroes(ng2-1), wlerrtmp(ng2-1))
    allocate(wlerr(ng2,ng2))

    wlerr = 0.0; lmodzeroes = 0.0; wlerrtmp = 0.0

    wdim = nesub
    allocate(modzeroes(nesub-1), werrtmp(nesub-1))
    allocate(werr(negrid-1,nesub))
    
    werr = 0.0 ; modzeroes = 0.0 ; werrtmp = 0.0
    
    ! loop to obtain weights for energy grid points.  negrid-1 sets
    ! of weights are needed because we want to compute integrals
    ! for negrid-1 sets of energy points (corresponding to negrid-1
    ! points that we can choose to drop from the guassian quadrature)
    
    do ipt=1,nesub
       
       ! drops the point corresponding to ipt from the energy grid
       
       if (ipt /= 1) modzeroes(:ipt-1) = zeroes(:ipt-1)
       if (ipt /= nesub) modzeroes(ipt:nesub-1) = zeroes(ipt+1:nesub)
       
       ! get weights for energy grid points
       
       call get_weights (nmax,0.0,vcut,modzeroes,werrtmp,ndiv,divmax,eflag)
       
       ! a zero is left in the position corresponding to the dropped point
       
       if (ipt /= 1) werr(:ipt-1,ipt) = werrtmp(:ipt-1)
       if (ipt /= nesub) werr(ipt+1:nesub,ipt) = werrtmp(ipt:nesub-1)
       werr(nesub+1:,ipt) = w(nesub+1:negrid-1)
       
       ! absorbing volume element into weights
       werr(:nesub,ipt) = werr(:nesub,ipt)*energy(:nesub)*exp(-energy(:nesub))/sqrt(pi)
       
    end do

    ! same thing done here for lamdba as was
    ! done earlier for energy space

    do ipt=1,ng2

       if (ipt /= 1) lmodzeroes(:ipt-1) = xx(:ipt-1)
       if (ipt /= ng2) lmodzeroes(ipt:ng2-1) = xx(ipt+1:)

       call get_weights (nmax,1.0,0.0,lmodzeroes,wlerrtmp,ndiv,divmax,eflag)

       if (ipt /= 1) wlerr(:ipt-1,ipt) = wlerrtmp(:ipt-1)
       if (ipt /= ng2) wlerr(ipt+1:,ipt) = wlerrtmp(ipt:)

!       do il = 1, ng2
          ! TEMP FOR TESTING -- MAB
!          write (*,*) 'wlerr', ipt, il, xx(il), wlerr(il,ipt), ndiv, divmax
!       end do

    end do

    deallocate(modzeroes,werrtmp,lmodzeroes,wlerrtmp)
    eflag = .false.

  end subroutine init_weights

! the get_weights subroutine determines how to divide up the integral into 
! subintervals and how many grid points should be in each subinterval

  subroutine get_weights (maxpts_in, llim, ulim, nodes, wgts, ndiv, divmax, err_flag)

    implicit none

    integer, intent (in) :: maxpts_in
    real, intent (in) :: llim, ulim
    real, dimension (:), intent (in) :: nodes
    real, dimension (:), intent (out) :: wgts
    logical, intent (out) :: err_flag
    integer, intent (out) :: ndiv, divmax

    integer :: npts, rmndr, basepts, divrmndr, base_idx, idiv, epts, im, maxpts
    integer, dimension (:), allocatable :: divpts

    real :: wgt_max

! npts is the number of grid points in the integration interval
    npts = size(nodes)

    wgts = 0.0; epts = npts; basepts = nmax; divrmndr = 0; ndiv = 1; divmax = npts

! maxpts is the max number of pts in an integration subinterval
    maxpts = min(maxpts_in,npts)

    do

!       wgt_max = wgt_fac/maxpts
       wgt_max = abs(ulim-llim)*wgt_fac/npts

! only need to subdivide integration interval if maxpts < npts
       if (maxpts .ge. npts) then
          call get_intrvl_weights (llim, ulim, nodes, wgts)
       else
          rmndr = mod(npts-maxpts,maxpts-1)
          
! if rmndr is 0, then each subinterval contains maxpts pts
          if (rmndr == 0) then
! ndiv is the number of subintervals
             ndiv = (npts-maxpts)/(maxpts-1) + 1
             allocate (divpts(ndiv))
! divpts is an array containing the # of pts for each subinterval
             divpts = maxpts
          else
             ndiv = (npts-maxpts)/(maxpts-1) + 2
             allocate (divpts(ndiv))
! epts is the effective number of pts after taking into account double
! counting of some grid points (those that are boundaries of subintervals
! are used twice)
             epts = npts + ndiv - 1
             basepts = epts/ndiv
             divrmndr = mod(epts,ndiv)
             
! determines if all intervals have same # of pts
             if (divrmndr == 0) then
                divpts = basepts
             else
                divpts(:divrmndr) = basepts + 1
                divpts(divrmndr+1:) = basepts
             end if
          end if
          
          base_idx = 0
          
! loop calls subroutine to get weights for each subinterval
          do idiv=1,ndiv
             if (idiv == 1) then
                call get_intrvl_weights (llim, nodes(base_idx+divpts(idiv)), &
                     nodes(base_idx+1:base_idx+divpts(idiv)),wgts(base_idx+1:base_idx+divpts(idiv)))
             else if (idiv == ndiv) then
                call get_intrvl_weights (nodes(base_idx+1), ulim, &
                     nodes(base_idx+1:base_idx+divpts(idiv)),wgts(base_idx+1:base_idx+divpts(idiv)))
             else
                call get_intrvl_weights (nodes(base_idx+1), nodes(base_idx+divpts(idiv)), &
                     nodes(base_idx+1:base_idx+divpts(idiv)),wgts(base_idx+1:base_idx+divpts(idiv)))
             end if
             base_idx = base_idx + divpts(idiv) - 1
          end do
          
          divmax = maxval(divpts)

          deallocate (divpts)
       end if

! check to make sure the weights do not get too large
       if (abs(maxval(wgts)) .gt. wgt_max) then
          if (maxpts .lt. 3) then
             err_flag = .true.
             exit
          end if
          maxpts = divmax - 1
       else
          exit
       end if

       wgts = 0.0; epts = npts; divrmndr = 0; basepts = nmax
    end do

  end subroutine get_weights

  subroutine get_intrvl_weights (llim, ulim, nodes, wgts)
    use gauss_quad, only: get_legendre_grids_from_cheb
    
    implicit none
    
    ! llim (ulim) is lower (upper) limit of integration
    real, intent (in) :: llim, ulim
    real, dimension (:), intent (in) :: nodes
    real, dimension (:), intent (in out) :: wgts
    
    ! stuff needed to do guassian quadrature 
    real, dimension (:), allocatable :: gnodes, gwgts, omprod
    integer :: ix, iw

    allocate (gnodes(size(nodes)/2+1), gwgts(size(wgts)/2+1), omprod(size(nodes)/2+1))
    
    call get_legendre_grids_from_cheb (llim, ulim, gnodes, gwgts)

    do iw=1,size(wgts)
       omprod = 1.0
       
       do ix=1,size(nodes)
          if (ix /= iw) omprod = omprod*(gnodes - nodes(ix))/(nodes(iw) - nodes(ix))
       end do
       
       do ix=1,size(gwgts)
          wgts(iw) = wgts(iw) + omprod(ix)*gwgts(ix)
       end do
    end do

    deallocate (gnodes, gwgts, omprod)
       
  end subroutine get_intrvl_weights

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
!    complex, dimension (:,:), allocatable :: geint
    complex, dimension (:), allocatable :: work
    real :: fac
    integer :: is, il, ie, ik, it, iglo, ig, i

    total = 0.
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       fac = weights(is)*w(ie)

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

  subroutine integrate_test (g, weights, total, istep)

    use theta_grid, only: ntgrid, bmag
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, idx, idx_local
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_allreduce
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:), intent (out) :: total
    integer, intent (in) :: istep


    complex, dimension (:), allocatable :: work
    real :: fac
    integer :: is, il, ie, ik, it, iglo, ig, i

    real, dimension (:,:), allocatable :: ypt

    allocate(ypt(-ntgrid:ntgrid,nlambda))
    ypt = 0.0

    total = 0.
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       fac = weights(is)*w(ie)

       total(:, it, ik) = total(:, it, ik) + fac*wl(:,il)
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

    deallocate(ypt)

  end subroutine integrate_test

  subroutine legendre_transform (g, tote, totl, istep, tott)
    
    use egrid, only: zeroes
    use mp, only: nproc, broadcast
    use theta_grid, only: ntgrid, bmag, bmax
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, idx, idx_local
    use mp, only: sum_reduce, proc0
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (0:,-ntgrid:,:,:,:), intent (out) :: tote, totl
    complex, dimension (0:,-ntgrid:,:,:,:), intent (out), optional :: tott
    integer, intent (in) :: istep

    complex :: totfac
    complex, dimension (:), allocatable :: worke, workl, workt
    real :: fac, ulim
    integer :: is, il, ie, ik, it, iglo, ig, i, j, im, ntrap, k
    integer, save :: lpesize
!    logical :: first = .true.

    real, dimension (:), allocatable :: nodes
    real, dimension (:,:), allocatable :: lpltmp, lpttmp
!    real, dimension (:,:), allocatable, save :: lpe, lpl
!    real, dimension (:,:,:,:), allocatable, save :: lpt

!    if (first) then

    if (.not. allocated(lpl)) then
       allocate(lpltmp(ng2,0:ng2-1))
       allocate(lpl(nlambda,0:ng2-1))

       lpesize = nesub
       allocate(lpe(negrid,0:lpesize-1)) ; lpe = 0.0
       
       ! get value of first nesub legendre polynomials
       ! at each of the grid points on (0,vcut)
       call legendre_polynomials (0.0,vcut,zeroes(:lpesize),lpe(:lpesize,:))
       ! TEMP FOR TESTING -- MAB
       !          lpe = 2.*lpe/vcut
       
       ! get value of first ng2 legendre polynomials
       ! at each of the grid points on (0,1)
       call legendre_polynomials (0.0,1.0,xx,lpltmp)

       lpl = 0.0
       lpl(1:ng2,:) = lpltmp

       if (present(tott)) then
          allocate (lpt(nlambda,0:2*(nlambda-ng2-1),-ntgrid:ntgrid,2))
          lpt = 0.0
          do ig = -ntgrid, ntgrid
             ntrap = 1
             if (jend(ig) > ng2+1) then
                ntrap = jend(ig)-ng2
                allocate (nodes(2*ntrap-1))
                allocate (lpttmp(2*ntrap-1,0:2*(ntrap-1)))
                do il = 1, ntrap
                   nodes(il) = -sqrt(max(0.0,1.0-al(ng2+il)*bmag(ig)))
                end do
                nodes(ntrap+1:) = -nodes(ntrap-1:1:-1)
! TEMP FOR TESTING -- MAB
!                nodes = nodes + sqrt(1.0-bmag(ig)/bmax)
!                ulim = 2.*sqrt(1.0-bmag(ig)/bmax)
                ulim = sqrt(1.0-bmag(ig)/bmax)
                call legendre_polynomials (-ulim,ulim,nodes,lpttmp)
                lpt(ng2+1:jend(ig),0:2*(ntrap-1),ig,2) = lpttmp(1:ntrap,:)
                lpt(ng2+1:jend(ig)-1,0:2*(ntrap-1),ig,1) = lpttmp(2*ntrap-1:ntrap+1:-1,:)
!                lpt(ng2+1:jend(ig),0:2*(ntrap-1),ig,1) = lpttmp(2*ntrap-1:ntrap+1:-1,:)
!                do ie = 0, 2*(ntrap-1)
!                   do il = 1, 2*ntrap-1
!                      if (proc0) write (*,*) 'lptrap', ig, ntrap, ulim, ie, il, nodes(il), lpttmp(il,ie)
!                   end do
!                end do
!                do ie = 0, 2*(ntrap-1)
!                   do il = ng2+1,jend(ig)
!                      write (*,*) 'lpt', ig, ie, il, lpt(il,ie,ig,1), lpt(il,ie,ig,2)
!                   end do
!                end do
                deallocate (nodes, lpttmp)
             end if
          end do
       end if

       deallocate (lpltmp)
!       first = .false.
    end if

    ! carry out legendre transform to get coefficients of
    ! legendre polynomial expansion of g
    totfac = 0. ; tote = 0. ; totl = 0.
    if (present(tott)) tott = 0.
    do is = 1, nspec
       do ie = 1, negrid
          do il = 1, nlambda
             do it = 1, ntheta0
                do ik = 1, naky
                   iglo = idx (g_lo, ik, it, il, ie, is)
                   if (idx_local (g_lo, iglo)) then
                      do ig=-ntgrid,ntgrid
                         totfac = w(ie)*wl(ig,il)*(g(ig,1,iglo)+g(ig,2,iglo))
                         do im=0,lpesize-1
                            tote(im, ig, it, ik, is) = tote(im, ig, it, ik, is) + totfac*lpe(ie,im)*(2*im+1)
                         end do
                         do im=0,ng2-1
                            totl(im, ig, it, ik, is) = totl(im, ig, it, ik, is) + totfac*lpl(il,im)*(2*im+1)
                         end do
                         if (present(tott)) then
                            do im=0,2*(jend(ig)-ng2-1)
                               tott(im, ig, it, ik, is) = tott(im, ig, it, ik, is) + &
                                    w(ie)*wl(ig,il)*(lpt(il,im,ig,1)*g(ig,1,iglo)+lpt(il,im,ig,2)*g(ig,2,iglo))*(2*im+1)
                            end do
                         end if
                      end do
                   end if
                end do
             end do
          end do
       end do
    end do

    if (nproc > 1) then
       allocate (worke((2*ntgrid+1)*naky*ntheta0*nspec*lpesize)) ; worke = 0.
       allocate (workl((2*ntgrid+1)*naky*ntheta0*nspec*ng2)) ; workl = 0.

       if (present(tott)) then
          allocate (workt((2*ntgrid+1)*naky*ntheta0*nspec*(2*(nlambda-ng2)-1)))
          workt = 0.
       end if

       i = 0 ; j = 0 ; k = 0

       do is = 1, nspec
          do ik = 1, naky
             do it = 1, ntheta0
                do ig = -ntgrid, ntgrid
                   do im = 0, lpesize-1
                      i = i + 1
                      worke(i) = tote(im, ig, it, ik, is)
                   end do
                   do im = 0, ng2-1
                      j = j + 1
                      workl(j) = totl(im, ig, it, ik, is)
                   end do
                   if (present(tott)) then
                      do im = 0, 2*(nlambda-ng2-1)
                         k = k + 1
                         workt(k) = tott(im, ig, it, ik, is)
                      end do
                   end if
                end do
             end do
          end do
       end do

       call sum_reduce (worke, 0)
       call sum_reduce (workl, 0)
       if (present(tott)) call sum_reduce (workt, 0)

       if (proc0) then
          i = 0 ; j = 0 ; k = 0
          do is = 1, nspec
             do ik = 1, naky
                do it = 1, ntheta0
                   do ig = -ntgrid, ntgrid
                      do im = 0, lpesize-1
                         i = i + 1
                         tote(im, ig, it, ik, is) = worke(i)
                      end do
                      do im = 0, ng2-1
                         j = j + 1
                         totl(im, ig, it, ik, is) = workl(j)
                      end do
                      if (present(tott)) then
                         do im = 0, 2*(nlambda-ng2-1)
                            k = k + 1
                            tott(im, ig, it, ik, is) = workt(k)
                         end do
                      end if
                   end do
                end do
             end do
          end do
       end if

       deallocate (worke,workl)
       if (present(tott)) deallocate (workt)
    end if

  end subroutine legendre_transform

! TEMP FOR TESTING -- MAB
!  subroutine legendre_polynomials (ulim, xptsdum, lpdum)
  subroutine legendre_polynomials (llim, ulim, xptsdum, lpdum)

    double precision, dimension (:), allocatable :: lp1, lp2, lp3, zshift

    real, intent (in) :: ulim, llim
    real, dimension (:), intent (in)   :: xptsdum
    real, dimension (:,0:), intent(out) :: lpdum

    integer :: j, im, mmax

    lpdum = 0.0

!    nmax = size(lpdum(1,:))
    mmax = size(xptsdum)

    allocate(lp1(mmax),lp2(mmax),lp3(mmax),zshift(mmax))

    lp1 = real(1.0,kind(lp1(1)))
    lp2 = real(0.0,kind(lp2(1)))

    lpdum(:,0) = real(1.0,kind(lpdum))

! TEMP FOR TESTING -- MAB
!    zshift = real(2.0,kind(zshift))*xptsdum/ulim - real(1.0,kind(zshift))
    zshift = real(2.0,kind(zshift))*(xptsdum-llim)/(ulim-llim) - real(1.0,kind(zshift))

    do j=1, size(lpdum(1,:))-1
       lp3 = lp2
       lp2 = lp1
       lp1 = ((2*j-1) * zshift * lp2 - (j-1) * lp3) / j
       lpdum(:,j) = lp1
    end do

    deallocate(lp1,lp2,lp3,zshift)

  end subroutine legendre_polynomials

  subroutine lagrange_interp (g, poly, istep, all)

    use mp, only: nproc, iproc
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, idx, idx_local
    use mp, only: sum_reduce, proc0, sum_allreduce
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:,:,:), intent (out) :: poly
    integer, optional, intent(in) :: all
    integer, intent(in) :: istep

    complex, dimension (:), allocatable :: work
    real, dimension (:,:,:), allocatable :: ypt
    integer :: il, ie, ik, it, iglo, ig, i, ix, tsize

    allocate(ypt(-ntgrid:ntgrid,nlambda,2))

    ypt = 0.0

    poly = 0.
    do ie = 1, negrid
       do il = ng2+1, nlambda
          do it = 1, ntheta0
             do ik = 1, naky
                iglo = idx (g_lo, ik, it, il, ie, 1)
                if (idx_local (g_lo, iglo)) then
                   do ig = -ntgrid, ntgrid
                      if (.not. forbid(ig,il)) then
                         ypt(ig,il,1) = sqrt(max(1.0 - bmag(ig)*al(il),0.0))                      
                         ypt(ig,il,2) = -ypt(ig,il,1)
                      end if
                      poly(ig, it, ik, ie, :) = poly(ig, it, ik, ie, :) + &
                           lgrnge(ig,il,:,1)*cos(0.1*istep*ypt(ig,il,1)+1.0) + &
                           lgrnge(ig,il,:,2)*cos(0.1*istep*ypt(ig,il,2)+1.0)  !(g(ig,1,iglo)+g(ig,2,iglo))
                   end do
                end if
             end do
          end do
       end do
    end do

    tsize = 2*nterp-1

    if (nproc > 1) then

       allocate (work((2*ntgrid+1)*naky*ntheta0*negrid*tsize)) ; work = 0.

       i = 0
       do ix = 1, tsize
          do ie = 1, negrid
             do ik = 1, naky
                do it = 1, ntheta0
                   do ig = -ntgrid, ntgrid
                      i = i + 1
                      work(i) = poly(ig, it, ik, ie, ix)
                   end do
                end do
             end do
          end do
       end do
       
       if (present(all)) then
          call sum_allreduce (work)
       else
          call sum_reduce (work, 0)
       end if

       if (proc0 .or. present(all)) then
          i = 0
          do ix = 1, tsize
             do ie = 1, negrid
                do ik = 1, naky
                   do it = 1, ntheta0
                      do ig = -ntgrid, ntgrid
                         i = i + 1
                         poly(ig, it, ik, ie, ix) = work(i)
                      end do
                   end do
                end do
             end do
          end do
       end if
       deallocate (work)
    end if

    deallocate(ypt)

  end subroutine lagrange_interp

!  subroutine integrate_moment (g, total, all)
  subroutine integrate_moment_c34 (g, total, all)
! returns results to PE 0 [or to all processors if 'all' is present in input arg list]
! NOTE: Takes f = f(x, y, z, sigma, lambda, E, species) and returns int f, where the integral
! is over all velocity space
    use mp, only: nproc, iproc
    use theta_grid, only: ntgrid
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
! TT>
!    use gs2_layouts, only: g_lo, idx, idx_local
    use gs2_layouts, only: g_lo, is_idx, ik_idx, it_idx, ie_idx, il_idx
! <TT
    use mp, only: sum_reduce, proc0, sum_allreduce
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    integer, optional, intent(in) :: all

    complex, dimension (:), allocatable :: work
    real :: fac
    integer :: is, il, ie, ik, it, iglo, ig, i
!    logical :: only = .true.

    total = 0.
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)

       do ig = -ntgrid, ntgrid
          total(ig, it, ik, is) = total(ig, it, ik, is) + &
               w(ie)*wl(ig,il)*(g(ig,1,iglo)+g(ig,2,iglo))
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
       
       if (present(all)) then
          call sum_allreduce (work)
       else
          call sum_reduce (work, 0)
       end if

       if (proc0 .or. present(all)) then
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

  end subroutine integrate_moment_c34

!  subroutine integrate_moment (g, total, all)
  subroutine integrate_moment_r33 (g, total, all)
! returns results to PE 0 [or to all processors if 'all' is present in input arg list]
! NOTE: Takes f = f(y, z, sigma, lambda, E, species) and returns int f, where the integral
! is over all velocity space
    use mp, only: nproc, iproc
    use theta_grid, only: ntgrid
    use species, only: nspec
    use kt_grids, only: naky
    use gs2_layouts, only: p_lo, is_idx, ik_idx, ie_idx, il_idx
    use mp, only: sum_reduce, proc0, sum_allreduce
    implicit none
    real, dimension (-ntgrid:,:,p_lo%llim_proc:), intent (in) :: g
    real, dimension (-ntgrid:,:,:), intent (out) :: total
    integer, optional, intent(in) :: all

    real, dimension (:), allocatable :: work
    real :: fac
    integer :: is, il, ie, ik, iplo, ig, i

    total = 0.
    do iplo = p_lo%llim_proc, p_lo%ulim_proc
       ik = ik_idx(p_lo,iplo)
       ie = ie_idx(p_lo,iplo)
       is = is_idx(p_lo,iplo)
       il = il_idx(p_lo,iplo)

       do ig = -ntgrid, ntgrid
          total(ig, ik, is) = total(ig, ik, is) + &
               w(ie)*wl(ig,il)*(g(ig,1,iplo)+g(ig,2,iplo))
       end do
    end do

    if (nproc > 1) then
       allocate (work((2*ntgrid+1)*naky*nspec)) ; work = 0.
       i = 0
       do is = 1, nspec
          do ik = 1, naky
             do ig = -ntgrid, ntgrid
                i = i + 1
                work(i) = total(ig, ik, is)
             end do
          end do
       end do

       if (present(all)) then
          call sum_allreduce (work)
       else
          call sum_reduce (work, 0)
       end if

       if (proc0 .or. present(all)) then
          i = 0
          do is = 1, nspec
             do ik = 1, naky
                do ig = -ntgrid, ntgrid
                   i = i + 1
                   total(ig, ik, is) = work(i)
                end do
             end do
          end do
       end if
       deallocate (work)
    end if

  end subroutine integrate_moment_r33

  subroutine integrate_moment_lec (lo, g, total)

    use theta_grid, only: ntgrid
    use layouts_type, only: le_layout_type
    use gs2_layouts, only: ig_idx, it_idx, ik_idx, is_idx
    type (le_layout_type), intent (in) :: lo
    complex, dimension (:,:,lo%llim_proc:), intent (in) :: g
!    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    complex, dimension (lo%llim_proc:), intent (out) :: total
    integer :: ixi, nxi, ie, il, ile, is, it, ik, ig
    real :: fac

    nxi = max(2*nlambda-1,2*ng2)
    total = 0.0
    do ile = lo%llim_proc, lo%ulim_proc
       ig = ig_idx (lo,ile)
       it = it_idx (lo,ile)
       ik = ik_idx (lo,ile)
       is = is_idx (lo,ile)
       do ie=1, negrid
          do ixi=1, nxi
             il = min(ixi, nxi+1-ixi)
             fac = w(ie) * wl(ig,il)
!             total(ig,it,ik,is) = total(ig,it,ik,is) + fac * g(ixi,ie,ile)
             total(ile) = total(ile) + fac * g(ixi,ie,ile)
          end do
       end do
    end do

    ! No need for communication since all velocity grid points are together
    ! and each prcessor does not touch the unset place
    ! They actually don't need to keep all 4D array
    ! Do we stay in le_layout for total?
    ! --- ile contains necessary and sufficient information for (ig,it,ik,is)

  end subroutine integrate_moment_lec

  subroutine integrate_kysum (g, ig, total, all)
! returns results to PE 0 [or to all processors if 'all' is present in input arg list]
! NOTE: Takes f = f(y, lambda, E, species) and returns int sum_{ky} f, where the integral
! is over energy and lambda (not sigma)
    use constants, only: zi
    use mp, only: nproc, iproc
    use theta_grid, only: ntgrid
    use species, only: nspec
    use kt_grids, only: naky, ntheta0, aky
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx, p_lo
    use mp, only: sum_reduce, proc0, sum_allreduce
    implicit none
    complex, dimension (p_lo%llim_proc:), intent (in) :: g
    integer, intent (in) :: ig
    complex, dimension (:), intent (out) :: total
    integer, optional, intent(in) :: all

    complex, dimension (negrid,nlambda,nspec) :: gksum
    complex, dimension (:), allocatable :: work
    real :: fac
    integer :: is, il, ie, ik, it, iplo, i

    total = 0. ; gksum = 0.
    do iplo = p_lo%llim_proc, p_lo%ulim_proc
       ik = ik_idx(p_lo,iplo)
       ie = ie_idx(p_lo,iplo)
       is = is_idx(p_lo,iplo)
       il = il_idx(p_lo,iplo)
       gksum(ie,il,is) = gksum(ie,il,is) + real(aky(ik)*g(iplo)) + zi*aimag(aky(ik)*g(iplo))
    end do
    ! real part of gksum is | sum_{ky} ky * J0 * real[ ky*(conjg(phi+)*h- + conjg(phi-)*h+ ] |**2
    ! imag part of gksum is | sum_{ky} ky * J0 * aimag[ ky*(conjg(phi+)*h- + conjg(phi-)*h+ ] |**2
    gksum = real(gksum)**2 + zi*aimag(gksum)**2

    do iplo = p_lo%llim_proc, p_lo%ulim_proc
       ie = ie_idx(p_lo,iplo)
       is = is_idx(p_lo,iplo)
       il = il_idx(p_lo,iplo)

       total(is) = total(is) + w(ie)*wl(ig,il)*gksum(ie,il,is)
    end do

    if (nproc > 1) then
       allocate (work(nspec)) ; work = 0.
       i = 0
       do is = 1, nspec
          i = i + 1
          work(i) = total(is)
       end do
       
       if (present(all)) then
          call sum_allreduce (work)
       else
          call sum_reduce (work, 0)
       end if

       if (proc0 .or. present(all)) then
          i = 0
          do is = 1, nspec
             i = i + 1
             total(is) = work(i)
          end do
       end if
       deallocate (work)
    end if

  end subroutine integrate_kysum

  subroutine lint_error (g, weights, total)
    use theta_grid, only: ntgrid, bmag, bmax
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, idx, idx_local
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_allreduce, proc0, broadcast
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    complex, dimension (:), allocatable :: work
!    real, dimension (:,:,:), allocatable, save :: wlmod
    real :: fac
    integer :: is, il, ie, ik, it, iglo, ig, i, ipt
!    logical, save :: first = .true.

!    if (first) then
    if (.not. allocated (wlmod)) then
       if (proc0) then
          allocate (wlmod(-ntgrid:ntgrid,nlambda,ng2))
          wlmod = 0.0
          do ipt = 1, ng2
             do il = 1, ng2
                do ig = -ntgrid, ntgrid
                   wlmod(ig,il,ipt) = wlerr(il,ipt)*2.0*sqrt((bmag(ig)/bmax) &
                        *((1.0/bmax-al(il))/(1.0/bmag(ig)-al(il))))
                end do
             end do
             if (nlambda > ng2) wlmod(:,ng2+1:,ipt) = wl(:,ng2+1:)
          end do
       end if
       
       if (.not. proc0) then
          allocate(wlmod(-ntgrid:ntgrid,nlambda,ng2))
          wlmod = 0.0
       end if
       
       do ipt=1,ng2
          do il=1,nlambda
             call broadcast (wlmod(:,il,ipt))
          end do
       end do
       
!       first = .false.
    end if

    do ipt=1,ng2
       total(:,:,:,ipt) = 0.
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          fac = weights(is)*w(ie)

          total(:, it, ik, ipt) = total(:, it, ik, ipt) + fac*wlmod(:,il,ipt)*(g(:,1,iglo)+g(:,2,iglo))
       end do

       allocate (work((2*ntgrid+1)*naky*ntheta0)) ; work = 0.
       i = 0
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntgrid, ntgrid
                i = i + 1
                work(i) = total(ig, it, ik, ipt)
             end do
          end do
       end do
       
       call sum_allreduce (work) 
       
       i = 0
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntgrid, ntgrid
                i = i + 1
                total(ig, it, ik, ipt) = work(i)
             end do
          end do
       end do
       deallocate (work)
    end do

  end subroutine lint_error

  subroutine trap_error (g, weights, total)
    use theta_grid, only: ntgrid, bmag, bmax
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, idx, idx_local
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_allreduce, proc0, broadcast
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    complex, dimension (:), allocatable :: work
!    real, dimension (:,:,:), allocatable, save :: wtmod
!    real, dimension (:,:), allocatable, save :: ypts2
    real :: fac
    integer :: is, il, ie, ik, it, iglo, ig, i, ipt, ntrap
!    logical, save :: first = .true.

    ntrap = nlambda - ng2

!    if (first) then
    if (.not. allocated(wtmod)) then
       if (proc0) then
          allocate (wtmod(-ntgrid:ntgrid,nlambda,ntrap))
          do ipt=1,ntrap
             wtmod(:,:ng2,ipt) = wl(:,:ng2)
          end do
! next line only to be used when testing!!!!
!          wtmod(:,:ng2,:) = 0.
          wtmod(:,ng2+1:,:) = wlterr(:,ng2+1:,:)
       else
          allocate (wtmod(-ntgrid:ntgrid,nlambda,ntrap))
          wtmod = 0.0
       end if

       do ipt=1,ntrap
          do il=1,nlambda
             call broadcast (wtmod(:,il,ipt))
          end do
       end do

!       allocate(ypts2(-ntgrid:ntgrid,nlambda))
!       ypts2 = 0.0

!       first = .false.
    end if

    total = 0.
    do ipt=1,ntrap
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          fac = weights(is)*w(ie)

!          do ig=-ntgrid,ntgrid
!             if (.not. forbid(ig,il)) then
!                ypts2(ig,il) = sqrt(max(1.0-bmag(ig)*al(il),0.0))
!             else
!                ypts2(ig,il) = 0
!             end if
!          end do

          total(:, it, ik, ipt) = total(:, it, ik, ipt) + fac*wtmod(:,il,ipt)*(g(:,1,iglo)+g(:,2,iglo))
       end do

       allocate (work((2*ntgrid+1)*naky*ntheta0)) ; work = 0.
       i = 0
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntgrid, ntgrid
                i = i + 1
                work(i) = total(ig, it, ik, ipt)
             end do
          end do
       end do
       
       call sum_allreduce (work) 
       
       i = 0
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntgrid, ntgrid
                i = i + 1
                total(ig, it, ik, ipt) = work(i)
             end do
          end do
       end do
       deallocate (work)
    end do

  end subroutine trap_error

  subroutine eint_error (g, weights, total)

    use theta_grid, only: ntgrid
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, idx, idx_local
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_allreduce, proc0, broadcast
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    complex, dimension (:), allocatable :: work
!    real, dimension (:,:), allocatable, save :: wmod
    real :: fac
    integer :: is, il, ie, ik, it, iglo, ig, i, ipt
!    logical, save :: first = .true.

!    if (first) then
    if (.not. allocated(wmod)) then
       if (proc0) then
          allocate (wmod(negrid,wdim))
          wmod = 0.0
          wmod(:negrid-1,:) = werr(:,:)
          wmod(negrid,:) = w(negrid)  
       end if

       if (.not. proc0) then
          allocate (wmod(negrid,wdim))
       end if

       do ie = 1, wdim
          call broadcast (wmod(:,ie))
       end do

!       first = .false.
    end if

    do ipt=1,wdim
       total(:,:,:,ipt) = 0.
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          fac = weights(is)*wmod(ie,ipt)

          total(:, it, ik, ipt) = total(:, it, ik, ipt) + fac*wl(:,il)*(g(:,1,iglo)+g(:,2,iglo))
       end do

       allocate (work((2*ntgrid+1)*naky*ntheta0)) ; work = 0.
       i = 0
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntgrid, ntgrid
                i = i + 1
                work(i) = total(ig, it, ik, ipt)
             end do
          end do
       end do
       
       call sum_allreduce (work) 
       
       i = 0
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntgrid, ntgrid
                i = i + 1
                total(ig, it, ik, ipt) = work(i)
             end do
          end do
       end do
       deallocate (work)
    end do

  end subroutine eint_error

  subroutine set_grids
    use species, only: init_species, nspec
    use egrid, only: setvgrid
    use theta_grid, only: init_theta_grid, ntgrid, nbset, bset, eps
    implicit none

    integer :: tsize

    call init_theta_grid
    call init_species

    allocate (energy(negrid), w(negrid), anon(negrid), dele(negrid))

    call setvgrid (vcut, negrid, energy, w, nesub)

    anon = 1.0

    tsize = 2*nterp-1

    dele(1) = energy(1)
    dele(2:) = energy(2:)-energy(:negrid-1)

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
    allocate (lgrnge(-ntgrid:ntgrid,nlambda,tsize,2))
    allocate (xloc(-ntgrid:ntgrid,tsize))
    if (nlambda-ng2 > 0) then
       al(ng2+1:nlambda) = 1.0/bset
    end if
    call lgridset
    delal(1) = al(1)
    delal(2:) = al(2:) - al(:nlambda-1)
  end subroutine set_grids

  subroutine lgridset

    use theta_grid, only: ntgrid, bmag, bmax, eps, ntheta
    use gauss_quad, only: get_legendre_grids_from_cheb
    use constants
    use file_utils, only: open_output_file, close_output_file

    use species, only: nspec

    implicit none

! note that xgauss and wgauss are transposed wrt original code

    real, dimension (2*ngauss) :: wx
    real, dimension (:), allocatable :: ytmp, yb, yberr, wb, wberrtmp
    real, dimension (:,:), allocatable :: wberr
    integer :: icnt, npts, ix, ntrap
    real :: wwo, llim, ulim
    logical :: eflag = .false.

    integer :: ig, il, ndiv, divmax, divmaxerr, ndiverr

    integer :: ie, is

    allocate (xx(2*ngauss))

    call get_legendre_grids_from_cheb (1., 0., xx, wx)

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

! pick integration method (new=high-order interp, old=finite difference)
    if (new_trap_int) then

! wlterr contains weights for less accurate trapped integrals (for error estimation)

       ntrap = nlambda - ng2                ! max number of trapped particles (occurs at outboard midplane)
       allocate(wlterr(-ntgrid:ntgrid,nlambda,ntrap))
       
       wlterr = 0.0
       
       do ig = -ntgrid, ntgrid
          npts = 0
          
! npts is the number of al values in the trapped integral (varies with theta)
          do il = ng2+1, nlambda
             if (1.0 - al(il)*bmag(ig) > -bouncefuzz) then
                npts = npts + 1
             end if
          end do

          if (npts > 0) then

! jend(ig) is total # of al grid points at each theta value
             jend(ig) = ng2 + npts
          
! ytmp is an array containing pitch angle grid points (for vpa >= 0) 
             allocate(ytmp(npts), yb(2*npts-1), wb(2*npts-1))
!             allocate(wberr(npts,npts))
             allocate(wberr(2*npts-1,npts))

             ytmp = 0.0; yb = 0.0; wb = 0.0; wberr = 0.0
          
!             icnt = 1

! loop computes transformed variable of integration
!             do il = ng2+1, nlambda
             do il = ng2+1, ng2+npts
                ytmp(il-ng2) = sqrt(max(1 - al(il)*bmag(ig), 0.0))
!                icnt = icnt + 1
!                write (*,*) 'ytmp', il-ng2, ytmp(il-ng2)
             end do

! define array (yb) with pitch-angle gridpts corresponding to both positive and negative vpa
             if (npts > 1) yb(:npts-1) = -ytmp(:npts-1)
             do ix=1,npts
                yb(ix+npts-1) = ytmp(npts-ix+1)
             end do

! get lagrange coefficients for construction of lagrange interpolating polynomial
! this is not necessary for integration, but can be useful diagnostic for integral accuracy
!          call lagrange_coefs (ig, yb, lgrnge(ig,:,:,:), xloc(ig,:))

! gets grid point weights for trapped particle integral
             ulim = sqrt(max(1.0-bmag(ig)/bmax,0.0))
             llim = -ulim

             if (ulim > 0.0) call get_weights (nmax, llim, ulim, yb, wb, ndiv, divmax, eflag) 

!             do il = 1, 2*npts-1
!                write (*,*) ig, il, 2*npts-1, ulim, yb(il), wb(il), ndiv, divmax
!             end do

! gets weights for less accurate trapped particle integrals (for purposes of error estimation)

! can't get error estimate for npts = 0 or 1
             if (npts > 1) then
                do ix=1,npts

                   if (ix == 1) then
! drop the first and last grid points from the integral
                      allocate (yberr(2*npts-3),wberrtmp(2*npts-3))
                      yberr = 0.0; wberrtmp = 0.0
                      yberr = yb(2:2*npts-2)
                   else if (ix == npts) then
! drop the vpa=0 grid point from the integral
                      allocate (yberr(2*npts-2),wberrtmp(2*npts-2))
                      yberr = 0.0; wberrtmp = 0.0
                      yberr(:npts-1) = yb(:npts-1)
                      yberr(npts:) = yb(npts+1:)
                   else
! drop the grid points corresponding to ix and its negative from the integral
                      allocate (yberr(2*npts-3),wberrtmp(2*npts-3))
                      yberr = 0.0; wberrtmp = 0.0
                      yberr(:ix-1) = yb(:ix-1)
                      yberr(ix:2*npts-ix-2) = yb(ix+1:2*npts-ix-1)
                      yberr(2*npts-ix-1:) = yb(2*npts-ix+1:)
                   end if

                   call get_weights (nmax, llim, ulim, yberr, wberrtmp, ndiverr, divmaxerr, eflag) 

! insert a weight of zero into indices corresponding to ix and its conjugate               
!                   if (ix /= 1) wberr(:ix-1,ix) = wberrtmp(:ix-1)
!                   if (ix /= npts) wberr(ix+1:,ix) = wberrtmp(ix:npts-1)
                   if (ix == 1) then
                      wberr(2:2*npts-2,1) = wberrtmp
                   else if (ix == npts) then
                      wberr(:npts-1,npts) = wberrtmp(:npts-1)
                      wberr(npts+1:,npts) = wberrtmp(npts:)
                   else
                      wberr(:ix-1,ix) = wberrtmp(:ix-1)
                      wberr(ix+1:2*npts-ix-1,ix) = wberrtmp(ix:2*npts-ix-2)
                      wberr(2*npts-ix+1,ix) = wberrtmp(2*npts-ix-1)
                   end if

!                wberr(npts,ix) = wberr(npts,ix)*0.5

! TEMP FOR TESTING -- MAB
!                   do il = 1, size(yb)
!                      write (*,*) 'wberr', ig, ix, yb(il), wberr(il,ix), npts
!                   end do
                   
                   deallocate (yberr,wberrtmp)

                end do
             end if

!             icnt = 1

!             do il = ng2+1, nlambda
!                if (1.0 - al(il)*bmag(ig) > -bouncefuzz) then
! avoid double counting of gridpoint at vpa=0               
!                if (icnt .eq. npts) then
!                   wl(ig,il) = wb(icnt)
!                   wlterr(ig,il,:npts) = wberr(icnt,:)
!                else
!                   wl(ig,il) = 2.0*wb(icnt)
!                   wlterr(ig,il,:npts) = 2.0*wberr(icnt,:)
!                end if
!                icnt = icnt + 1
!             end if
!             end do

             if (npts > 1) then
                do il = ng2+1, ng2+npts-1
!                   wl(ig,il) = 2.0*wb(il-ng2)
!                   wlterr(ig,il,:npts) = 2.0*wberr(il-ng2,:)
! take into account possible asymmetry of weights about xi = 0
! due to unequal # of grid points per integration interval
!                   wl(ig,il) = wb(il-ng2) + wb(il-ng2+npts)
                   wl(ig,il) = wb(il-ng2) + wb(2*npts-il+ng2)
!                   write (*,*) ig, il, wl(ig,il), wb(il-ng2), wb(2*npts-il+ng2), ng2, npts
                   wlterr(ig,il,:npts) = wberr(il-ng2,:) + wberr(2*npts-il+ng2,:)
                end do
             end if

             ! avoid double counting of gridpoint at vpa=0
             wl(ig,ng2+npts) = wb(npts)
             wlterr(ig,ng2+npts,:npts) = wberr(npts,:)

             deallocate(ytmp, yb, wb, wberr)

!             do il = ng2+1, nlambda
!                write (*,*) 'wgts', ig, il, sqrt(max(0.0,1.0-al(il)*bmag(ig))), wl(ig,il), sum(wl(ig,ng2+1:jend(ig))), sqrt(1.0-bmag(ig)/bmax), sum(wl(ig,ng2+1:jend(ig))*(1.0-bmag(ig)*al(ng2+1:jend(ig))))
!             end do
          end if
       end do

! old (finite-difference) integration scheme
    else
       jend = ng2 + 1
!       wlterr = 0.0
       do il = ng2+1, nlambda-1
          do ig = -ntgrid, ntgrid
             if (1.0-al(il)*bmag(ig) > -bouncefuzz &
                  .and. 1.0-al(il+1)*bmag(ig) > -bouncefuzz) &
                  then
                jend(ig) = jend(ig) + 1
                wwo = sqrt(max(1.0 -   al(il)*bmag(ig),0.0)) - &
                     sqrt(max(1.0 - al(il+1)*bmag(ig),0.0))
                wl(ig,il)   = wl(ig,il)   + wwo
                wl(ig,il+1) = wl(ig,il+1) + wwo
             end if
          end do
       end do

    end if

    do il = 1, ng2
       do ig = -ntgrid, ntgrid
          forbid(ig,il) = 1.0 - al(il)*bmag(ig) < -bouncefuzz
          if (forbid(ig,il)) then 
            call stop_message("Fatal error: supposedly passing particle was trapped, in legridset, in le_grids.f90")
         end if
       end do
    end do

    do il = ng2+1, nlambda
       do ig = -ntgrid, ntgrid
          forbid(ig,il) = 1.0 - al(il)*bmag(ig) < -bouncefuzz
       end do
    end do

    eflag = .false.

  end subroutine lgridset

  subroutine lagrange_coefs (ig, nodes, lfac, xloc)
    
    use theta_grid, only: bmag, bmax

    implicit none

    integer, intent (in) :: ig
    real, dimension (:), intent (in) :: nodes
    real, dimension (:,:,:), intent (out) :: lfac
    real, dimension (:), intent (out) :: xloc

    integer :: ix, il, j, icnt, tsize

    xloc(1) = -sqrt(max(1.0-bmag(ig)/bmax,0.0))
    xloc(2*nterp-1) = sqrt(max(1.0-bmag(ig)/bmax,0.0))
    do ix=2,2*nterp-1
       xloc(ix) = xloc(1) + (ix-1)*(xloc(2*nterp-1)-xloc(1))/(2*nterp-2)
    end do

    lfac = 0.
    icnt = 0

    tsize = size(nodes)

    do il=ng2+1,nlambda
       if (.not. forbid(ig,il)) then
          lfac(il,:,:) = 1.0
          icnt = icnt + 1
          if (icnt .eq. tsize-icnt+1) then
             lfac(il,:,1) = 0.
             do j=1,tsize
                if (j /= icnt) lfac(il,:,2) = lfac(il,:,2)*(xloc - nodes(j))/(nodes(icnt)-nodes(j))
             end do
          else
             do j=1,tsize
                if (j /= icnt) lfac(il,:,2) = lfac(il,:,2)*(xloc - nodes(j))/(nodes(icnt)-nodes(j))
                if (j /= tsize-icnt+1) lfac(il,:,1) = lfac(il,:,1)*(xloc - nodes(j))/(nodes(tsize-icnt+1)-nodes(j))
             end do
          end if
       end if
    end do

  end subroutine lagrange_coefs

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

  subroutine stop_message (message)
!JAB: print an error message and end program
    use file_utils, only: error_unit
    use mp, only: proc0, finish_mp
    implicit none
    character(*), intent (in) :: message
    integer :: ierr

    if (proc0) then
       ierr = error_unit()
       write (unit=ierr, fmt='(a)') message
    end if
    call finish_mp
    stop
  end subroutine stop_message

  subroutine integrate_volume_c (g, total, all)
! returns results to PE 0 [or to all processors if 'all' is present in input arg list]
! NOTE: Takes f = f(x, y, z, sigma, lambda, E, species) and returns int f, where the integral
! is over x-y space
    use mp, only: nproc, iproc
    use theta_grid, only: ntgrid
    use kt_grids, only: aky
    use species, only: nspec
    use gs2_layouts, only: g_lo, is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_reduce, proc0, sum_allreduce
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:,:,:), intent (out) :: total
    integer, optional, intent(in) :: all

    complex, dimension (:), allocatable :: work
    real :: fac
    integer :: is, il, ie, ik, it, iglo, ig, i, isgn

    total = 0.
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       if (aky(ik) == 0.) then
          fac = 1.0
       else
          fac = 0.5
       end if

       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             total(ig, il, ie, isgn, is) = total(ig, il, ie, isgn, is) + &
                  fac*g(ig,isgn,iglo)
          end do
       end do
    end do

    if (nproc > 1) then
       allocate (work((2*ntgrid+1)*nlambda*negrid*nspec*2)) ; work = 0.
       i = 0
       do is = 1, nspec
          do isgn = 1, 2
             do ie = 1, negrid
                do il = 1, nlambda
                   do ig = -ntgrid, ntgrid
                      i = i + 1
                      work(i) = total(ig, il, ie, isgn, is)
                   end do
                end do
             end do
          end do
       end do

       if (present(all)) then
          call sum_allreduce (work)
       else
          call sum_reduce (work, 0)
       end if

       if (proc0 .or. present(all)) then
          i = 0
          do is = 1, nspec
             do isgn = 1, 2
                do ie = 1, negrid
                   do il = 1, nlambda
                      do ig = -ntgrid, ntgrid
                         i = i + 1
                         total(ig, il, ie, isgn, is) = work(i)
                      end do
                   end do
                end do
             end do
          end do
       end if
       deallocate (work)
    end if

  end subroutine integrate_volume_c

  subroutine integrate_volume_r (g, total, all)
! returns results to PE 0 [or to all processors if 'all' is present in input arg list]
! NOTE: Takes f = f(x, y, z, sigma, lambda, E, species) and returns int f, where the integral
! is over x-y space
    use mp, only: nproc, iproc
    use theta_grid, only: ntgrid
    use kt_grids, only: aky
    use species, only: nspec
    use gs2_layouts, only: g_lo, is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_reduce, proc0, sum_allreduce
    implicit none
    real, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (-ntgrid:,:,:,:,:), intent (out) :: total
    integer, optional, intent(in) :: all

    real, dimension (:), allocatable :: work
    real :: fac
    integer :: is, il, ie, ik, it, iglo, ig, i, isgn

    total = 0.
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       if (aky(ik) == 0.) then
          fac = 1.0
       else
          fac = 0.5
       end if

       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             total(ig, il, ie, isgn, is) = total(ig, il, ie, isgn, is) + &
                  fac*g(ig,isgn,iglo)
          end do
       end do
    end do

    if (nproc > 1) then
       allocate (work((2*ntgrid+1)*nlambda*negrid*nspec*2)) ; work = 0.
       i = 0
       do is = 1, nspec
          do isgn = 1, 2
             do ie = 1, negrid
                do il = 1, nlambda
                   do ig = -ntgrid, ntgrid
                      i = i + 1
                      work(i) = total(ig, il, ie, isgn, is)
                   end do
                end do
             end do
          end do
       end do

       if (present(all)) then
          call sum_allreduce (work)
       else
          call sum_reduce (work, 0)
       end if

       if (proc0 .or. present(all)) then
          i = 0
          do is = 1, nspec
             do isgn = 1, 2
                do ie = 1, negrid
                   do il = 1, nlambda
                      do ig = -ntgrid, ntgrid
                         i = i + 1
                         total(ig, il, ie, isgn, is) = work(i)
                      end do
                   end do
                end do
             end do
          end do
       end if
       deallocate (work)
    end if

  end subroutine integrate_volume_r

  ! calculates and returns toroidal momentum flux as a function
  ! of vpar and theta
  subroutine get_flux_vs_theta_vs_vpa (f, vflx)

    use constants, only: pi
    use theta_grid, only: ntgrid, bmag
    use species, only: nspec

    implicit none

    real, dimension (-ntgrid:,:,:,:,:), intent (in) :: f
    real, dimension (-ntgrid:,:,:), intent (out) :: vflx

    real, dimension (:,:,:), allocatable :: favg
    real, dimension (:), allocatable, save :: vpa1d
    real, dimension (:,:), allocatable, save :: hermp1d
    real, dimension (:,:,:,:), allocatable, save :: vpapts
    real, dimension (:,:,:,:,:), allocatable, save :: hermp

    real :: fac
    integer :: is, il, ie, ig, isgn, iv
    integer :: norder

    norder = min(negrid, nlambda)/2

    allocate (favg(-ntgrid:ntgrid,nspec,0:norder-1))

    if (.not. allocated(vpapts)) then
       allocate (vpa1d(negrid*nlambda))
       allocate (hermp1d(negrid*nlambda,0:norder-1))
       allocate (vpapts(-ntgrid:ntgrid,nlambda,negrid,2))
       allocate (hermp(-ntgrid:ntgrid,nlambda,negrid,2,0:norder-1))
       vpapts = 0.0 ; hermp = 0.0 ; vpa1d = 0.0 ; hermp1d = 0.0

       do ie = 1, negrid
          do il = 1, nlambda
             do ig = -ntgrid, ntgrid
                vpapts(ig,il,ie,1) = sqrt(energy(ie)*max(0.0, 1.0-al(il)*bmag(ig)))
                vpapts(ig,il,ie,2) = -vpapts(ig,il,ie,1)
             end do
          end do
       end do

       do iv = 1, negrid*nlambda
          vpa1d(iv) = sqrt(energy(negrid))*(1. - 2.*(iv-1)/real(negrid*nlambda-1))
       end do

       call get_hermite_polynomials (vpa1d, hermp1d)
       call get_hermite_polynomials (vpapts, hermp)
    end if

    favg = 0.
    do is = 1, nspec
       do ie = 1, negrid
          do il = 1, nlambda
             do ig = -ntgrid, ntgrid
                favg(ig,is,:) = favg(ig,is,:) &
                     +w(ie)*wl(ig,il)*(hermp(ig,il,ie,1,:)*f(ig,il,ie,1,is) &
                     +hermp(ig,il,ie,2,:)*f(ig,il,ie,2,is))
             end do
          end do
       end do
    end do

    do is = 1, nspec
       do iv = 1, negrid*nlambda
          do ig = -ntgrid, ntgrid
             vflx(ig,iv,is) = sum(favg(ig,is,:)*hermp1d(iv,:))*exp(-vpa1d(iv)**2)
          end do
       end do
    end do

    deallocate (favg)
    
  end subroutine get_flux_vs_theta_vs_vpa

  subroutine get_hermite_polynomials_4d (xptsdum, hpdum)

    ! this subroutine returns Gn = Hn / sqrt(2^n n!) / pi^(1/4),
    ! where Hn are the hermite polynomials
    ! i.e. int dx Gm * Gn exp(-x^2) = 1

    use constants, only: pi
    use theta_grid, only: ntgrid

    implicit none

    real, dimension (-ntgrid:,:,:,:), intent (in)   :: xptsdum
    real, dimension (-ntgrid:,:,:,:,0:), intent (out) :: hpdum

    integer :: j
    double precision, dimension (:,:,:,:), allocatable :: hp1, hp2, hp3

    hpdum = 0.0

    allocate (hp1(-ntgrid:ntgrid,nlambda,negrid,2))
    allocate (hp2(-ntgrid:ntgrid,nlambda,negrid,2))
    allocate (hp3(-ntgrid:ntgrid,nlambda,negrid,2))

    hp1 = real(1.0,kind(hp1(0,1,1,1)))
    hp2 = real(0.0,kind(hp2(0,1,1,1)))

    hpdum(:,:,:,:,0) = 1.0

    do j=1, size(hpdum,5)-1
       hp3 = hp2
       hp2 = hp1
       hp1 = sqrt(2./j)*xptsdum*hp2 - sqrt(real(j-1)/j)*hp3
       hpdum(:,:,:,:,j) = hp1
    end do

    hpdum = hpdum/pi**(0.25)

    deallocate (hp1,hp2,hp3)

  end subroutine get_hermite_polynomials_4d

  subroutine get_hermite_polynomials_1d (xptsdum, hpdum)

    ! this subroutine returns Gn = Hn / sqrt(2^n n!) / pi^(1/4),
    ! where Hn are the hermite polynomials
    ! i.e. int dx Gm * Gn exp(-x^2) = 1

    use constants, only: pi

    implicit none

    real, dimension (:), intent (in)   :: xptsdum
    real, dimension (:,0:), intent (out) :: hpdum

    integer :: j
    double precision, dimension (:), allocatable :: hp1, hp2, hp3

    hpdum = 0.0

    allocate (hp1(size(xptsdum)))
    allocate (hp2(size(xptsdum)))
    allocate (hp3(size(xptsdum)))

    hp1 = real(1.0,kind(hp1(1)))
    hp2 = real(0.0,kind(hp2(1)))

    hpdum(:,0) = 1.0

    do j=1, size(hpdum,2)-1
       hp3 = hp2
       hp2 = hp1
       hp1 = sqrt(2./j)*xptsdum*hp2 - sqrt(real(j-1)/j)*hp3
       hpdum(:,j) = hp1
    end do

    hpdum = hpdum/pi**(0.25)

    deallocate (hp1,hp2,hp3)

  end subroutine get_hermite_polynomials_1d

  subroutine finish_le_grids

    use egrid, only: zeroes

    implicit none

    if (allocated(zeroes)) deallocate (zeroes)
    if (allocated(energy)) deallocate (energy, dele, al, wl, jend, forbid, xx, lgrnge, xloc)
    if (allocated(integration_work)) deallocate (integration_work)
    if (allocated(wlerr)) deallocate (wlerr)
    if (allocated(werr)) deallocate (werr)
    if (allocated(wlterr)) deallocate (wlterr)
    if (allocated(w)) deallocate (w)
    if (allocated(anon)) deallocate (anon)
    if (allocated(delal)) deallocate (delal)
    if (allocated(lpl)) deallocate (lpl, lpe)
    if (allocated(lpt)) deallocate (lpt)
    if (allocated(wlmod)) deallocate (wlmod)
    if (allocated(wmod)) deallocate (wmod)
    if (allocated(wtmod)) deallocate (wtmod)

    accel_x = .false. ; accel_v = .false.
    test = .false. ; trapped_particles = .true.
    new_trap_int = .false. 

    intinit = .false. ; slintinit = .false. ; lintinit = .false. ; eintinit = .false.

    initialized = .false.

  end subroutine finish_le_grids

end module le_grids


 
