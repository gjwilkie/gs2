module collisions

  use redistribute, only: redist_type

  implicit none

  public :: init_collisions
  public :: solfp1
  public :: reset_init

  private

  ! knobs
  real :: vncoef, absom
  integer :: ivnew
  logical :: conserve_number, conserve_momentum
  integer :: collision_model_switch
  logical :: use_shmem

  integer, parameter :: collision_model_lorentz = 1
  integer, parameter :: collision_model_krook = 2
  integer, parameter :: collision_model_none = 3
  integer, parameter :: collision_model_krook_test = 4
  integer, parameter :: collision_model_lorentz_test = 5

  real, dimension (:,:,:), allocatable :: vnew
  ! (naky,negrid,nspec) replicated

  ! only for krook
  real, dimension (:), allocatable :: vnewfe
  ! (-*- g_layout -*-)

  ! only for krook
  real, dimension (:,:), allocatable :: aintnorm
  ! (-ntgrid:ntgrid, -*- geint_layout -*-)

  ! only for lorentz
  real, dimension (:,:,:), allocatable :: sq
  ! (-ntgrid:ntgrid,nlambda,2) replicated

  ! only for lorentz
  real, dimension (:,:), allocatable :: c1, betaa, ql
  complex, dimension (:,:), allocatable :: glz
  ! (max(2*nlambda,2*ng2+1), -*- lz_layout -*-)
  ! (-ntgrid:ntgrid,naky,max(2*nlambda,2*ng2+1),negrid,nspec)

  ! momentum conservation
  complex, dimension (:,:), allocatable :: g3int
  ! (-ntgrid:ntgrid, -*- gint_layout -*-)

  type (redist_type), save :: lorentz_map

  logical :: initialized = .false.

contains

  subroutine init_collisions
    use species, only: init_species, nspec
    use theta_grid, only: init_theta_grid, ntgrid
    use kt_grids, only: init_kt_grids, naky, ntheta0
    use le_grids, only: init_le_grids, nlambda, negrid 
    use run_parameters, only: init_run_parameters
    use gs2_layouts, only: init_dist_fn_layouts
    implicit none

    if (initialized) return
    initialized = .true.

    call init_species
    call init_theta_grid
    call init_kt_grids
    call init_le_grids
    call init_run_parameters
    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec)

    call read_parameters
    call init_arrays
  end subroutine init_collisions

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit
    use text_options
    use mp, only: proc0, broadcast
    use shmem
    implicit none
    type (text_option), dimension (7), parameter :: modelopts = &
         (/ text_option('default', collision_model_lorentz), &
            text_option('lorentz', collision_model_lorentz), &
            text_option('krook', collision_model_krook), &
            text_option('krook-test', collision_model_krook_test), &
            text_option('lorentz-test', collision_model_lorentz_test), &
            text_option('none', collision_model_none), &
            text_option('collisionless', collision_model_none) /)
    character(20) :: collision_model
    namelist /collisions_knobs/ collision_model, vncoef, absom, ivnew, &
         conserve_number, conserve_momentum, use_shmem
    integer :: ierr

    if (proc0) then
       collision_model = 'default'
       vncoef = 0.6
       absom = 0.5
       ivnew = 0
       conserve_number = .true.
       conserve_momentum = .true.
       use_shmem = .true.
       read (unit=input_unit("collisions_knobs"), nml=collisions_knobs)

       ierr = error_unit()
       call get_option_value &
            (collision_model, modelopts, collision_model_switch, &
            ierr, "collision_model in collisions_knobs")
       use_shmem = use_shmem .and. shmem_available
    end if
    call broadcast (vncoef)
    call broadcast (absom)
    call broadcast (ivnew)
    call broadcast (conserve_number)
    call broadcast (conserve_momentum)
    call broadcast (collision_model_switch)
    call broadcast (use_shmem)
  end subroutine read_parameters

  subroutine init_arrays
    use species, only: nspec
    use le_grids, only: negrid
    implicit none
    real, dimension (negrid,nspec) :: hee

    if (collision_model_switch == collision_model_none) return

    call init_vnew (hee)
    if (all(abs(vnew(:,1,:)) <= 2.0*epsilon(0.0))) then
       collision_model_switch = collision_model_none
       return
    end if
    call init_g3int

    select case (collision_model_switch)
    case (collision_model_lorentz,collision_model_lorentz_test)
       call init_lorentz
    case (collision_model_krook,collision_model_krook_test)
       call init_krook (hee)
    end select
  end subroutine init_arrays

  subroutine init_vnew (hee)
    use species, only: nspec, spec, electron_species
    use le_grids, only: negrid, e
    use kt_grids, only: naky 
    use run_parameters, only: zeff, tunits
    use constants
    real, dimension (:,:), intent (out) :: hee
    integer :: ik, ie, is
    real :: v

    do is = 1, nspec
       do ie = 1, negrid
          v = sqrt(e(ie,is))
          hee(ie,is) = 1.0/sqrt(pi)/v*exp(-e(ie,is)) &
               + (1.0 - 0.5/e(ie,is)) &
                  *(1.0 - 1.0/(1.0          + v &
                             *(0.0705230784 + v &
                             *(0.0422820123 + v &
                             *(0.0092705272 + v &
                             *(0.0001520143 + v &
                             *(0.0002765672 + v &
                             *(0.0000430638)))))))**16)
       end do
    end do

    if(.not.allocated(vnew)) allocate (vnew(naky,negrid,nspec))
    do is = 1, nspec
       if (spec(is)%type == electron_species) then
          do ie = 1, negrid
             do ik = 1, naky
                vnew(ik,ie,is) = spec(is)%vnewk/e(ie,is)**1.5 &
                     *(zeff + hee(ie,is))*0.5*tunits(ik)
             end do
          end do
       else
          do ie = 1, negrid
             do ik = 1, naky
                vnew(ik,ie,is) = spec(is)%vnewk/e(ie,is)**1.5 &
                     *hee(ie,is)*0.5*tunits(ik)
             end do
          end do
       end if
    end do
  end subroutine init_vnew

  subroutine init_g3int
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: nlambda, al, lintegrate
    use gs2_layouts, only: g_lo, gint_lo, il_idx
    implicit none
    complex, dimension (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc) :: g
    real :: x
    integer :: il, ig, iglo

    if (.not. allocated(sq)) allocate (sq(-ntgrid:ntgrid,nlambda,2))
    do il = 1, nlambda
       do ig = -ntgrid, ntgrid
          x = sqrt(max(0.0, 1.0 - al(il)*bmag(ig)))
          sq(ig,il,1) =  x
          sq(ig,il,2) = -x
       end do
    end do

    if (.not.allocated(g3int)) &
         allocate (g3int(-ntgrid:ntgrid,gint_lo%llim_proc:gint_lo%ulim_alloc))
    g3int = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       x = al(il_idx(g_lo,iglo))
       g(:,1,iglo) = max(0.0, 1.0 - x*bmag(:))
       g(:,2,iglo) = max(0.0, 1.0 - x*bmag(:))
    end do
    call lintegrate (g, g3int)
  end subroutine init_g3int

  subroutine init_krook (hee)
    use species, only: spec, electron_species, ion_species
    use theta_grid, only: ntgrid, eps
    use le_grids, only: integrate, ng2, al, e
    use gs2_layouts, only: g_lo, geint_lo, ik_idx, il_idx, ie_idx, is_idx
    use run_parameters, only: delt, zeff, tunits
    implicit none
    real, dimension (:,:), intent (in) :: hee
    complex, dimension (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc) :: g
    complex, dimension (-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_alloc) &
         :: geint
    integer :: iglo, ik, il, ie, is, ige
    real :: zeff1, vep, vhat, delta00

    if (.not.allocated(aintnorm)) &
         allocate (aintnorm(-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_alloc))
    aintnorm = 0. ; g = 1.0
    call integrate (g, geint)

    do ige = geint_lo%llim_proc, geint_lo%ulim_proc
       aintnorm(:,ige) = 1.0/real(geint(:,ige))
    end do

    if (.not.allocated(vnewfe)) allocate (vnewfe(g_lo%llim_proc:g_lo%ulim_alloc))
    vnewfe = 0.

    if (collision_model_switch == collision_model_krook_test) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          vnewfe(iglo) = abs(spec(is)%vnewk)*tunits(ik)*delt
       end do
    else
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)

          if (spec(is)%type == electron_species &
               .or. spec(is)%type == ion_species) &
          then
             zeff1 = zeff
          else
             zeff1 = 0.0
          end if

          vep = abs(spec(is)%vnewk)*tunits(ik)
          vhat = sqrt(e(ie,is))
          if (ivnew > 0) then
             delta00 = absom/((vep + 2.0*epsilon(0.0))*zeff)
             if (ivnew >= 2) delta00 = (delta00*eps/37.2)**.333333333
             vnewfe(iglo) = delt*vep*eps*(zeff1 + hee(ie,is)) &
                  /(vhat**3*((1.0-eps-al(il))**2 + 1e-8)) &
                  *(.111*delta00+1.31)/(11.79*delta00+1.0)
          else
             vnewfe(iglo) = 0.00941/((1.0 - eps - al(il))**2 + 1e-8)
             if (il > ng2) vnewfe(iglo) = vnewfe(iglo) + vncoef/eps**2
             vnewfe(iglo) = vnewfe(iglo)*delt*vep*eps &
                  *(zeff1 + hee(ie,is))/vhat**3
          end if
       end do
    end if
  end subroutine init_krook

  subroutine init_lorentz
    use species, only: nspec, spec
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid, al, jend, ng2
    use run_parameters, only: delt, tunits
    use gs2_layouts, only: init_lorentz_layouts
    use gs2_layouts, only: lz_lo
    use gs2_layouts, only: ig_idx, ik_idx, ie_idx, is_idx
    implicit none
    integer :: ig, il, ilz, ik, ie, is, je
    real, dimension (nlambda+1) :: aa, bb, cc
    real, dimension (max(2*nlambda,2*ng2+1)) :: a1, b1
    real :: slb0, slb1, slb2, slbl, slbr, vn

    call init_lorentz_layouts &
         (ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2)
    call init_lorentz_redistribute

    if (.not.allocated(glz)) &
         allocate (glz(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
    glz = 0.0

    if (.not.allocated(c1)) then
       allocate (c1(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
       allocate (betaa(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
       allocate (ql(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
    endif

    c1 = 0.0 ; betaa = 0.0 ; ql = 0.0

    do ilz = lz_lo%llim_proc, lz_lo%ulim_proc
       is = is_idx(lz_lo,ilz)
       ik = ik_idx(lz_lo,ilz)
       ie = ie_idx(lz_lo,ilz)
       ig = ig_idx(lz_lo,ilz)
       je = jend(ig)
       if (collision_model_switch == collision_model_lorentz_test) then
          vn = abs(spec(is)%vnewk)*tunits(ik)
       else
          vn = vnew(ik,ie,is)
       end if
       if (je == 0) then
          do il = 2, ng2-1
             slb0 = sqrt(abs(1.0 - bmag(ig)*al(il-1)))
             slb1 = sqrt(abs(1.0 - bmag(ig)*al(il)))
             slb2 = sqrt(abs(1.0 - bmag(ig)*al(il+1)))

             slbl = (slb1 + slb0)/2.0
             slbr = (slb1 + slb2)/2.0

             cc(il) = -vn*delt*(1.0 - slbr*slbr)/(slbr - slbl)/(slb2 - slb1)
             aa(il) = -vn*delt*(1.0 - slbl*slbl)/(slbr - slbl)/(slb1 - slb0)
             bb(il) = 1.0 - (aa(il) + cc(il))
          end do

          slb1 = sqrt(abs(1.0-bmag(ig)*al(1)))
          slb2 = sqrt(abs(1.0-bmag(ig)*al(2)))

          slbr = (slb1 + slb2)/2.0

          cc(1) = -vn*delt*(-1.0 - slbr)/(slb2-slb1)
          aa(1) = 0.0
          bb(1) = 1.0 - (aa(1) + cc(1))

          il = ng2
          slb0 = sqrt(abs(1.0 - bmag(ig)*al(il-1)))
          slb1 = sqrt(abs(1.0 - bmag(ig)*al(il)))
          slb2 = -slb1

          slbl = (slb1 + slb0)/2.0
          slbr = (slb1 + slb2)/2.0

          cc(il) = -vn*delt*(1.0 - slbr*slbr)/(slbr - slbl)/(slb2 - slb1)
          aa(il) = -vn*delt*(1.0 - slbl*slbl)/(slbr - slbl)/(slb1 - slb0)
          bb(il) = 1.0 - (aa(il) + cc(il))

          a1(:ng2) = aa(:ng2)
          b1(:ng2) = bb(:ng2)
          c1(:ng2,ilz) = cc(:ng2)

          a1(ng2+1:2*ng2) = cc(ng2:1:-1)
          b1(ng2+1:2*ng2) = bb(ng2:1:-1)
          c1(ng2+1:2*ng2,ilz) =aa(ng2:1:-1)

          betaa(1,ilz) = 1.0/b1(1)
          do il = 1, 2*ng2-1
             ql(il+1,ilz) = a1(il+1)*betaa(il,ilz)
             betaa(il+1,ilz) = 1.0/(b1(il+1)-ql(il+1,ilz)*c1(il,ilz))
          end do

          ql(1,ilz) = 0.0
          ql(2*ng2+1:,ilz) = 0.0
          c1(2*ng2+1:,ilz) = 0.0
          betaa(2*ng2+1:,ilz) = 0.0

       else
          do il = 2, je-1
             slb0 = sqrt(abs(1.0 - bmag(ig)*al(il-1)))
             slb1 = sqrt(abs(1.0 - bmag(ig)*al(il)))
             slb2 = sqrt(abs(1.0 - bmag(ig)*al(il+1)))

             slbl = (slb1 + slb0)/2.0
             slbr = (slb1 + slb2)/2.0

             cc(il) = -vn*delt*(1.0 - slbr*slbr)/(slbr - slbl)/(slb2 - slb1)
             aa(il) = -vn*delt*(1.0 - slbl*slbl)/(slbr - slbl)/(slb1 - slb0)
             bb(il) = 1.0 - (aa(il) + cc(il))
          end do

          slb1 = sqrt(abs(1.0-bmag(ig)*al(1)))
          slb2 = sqrt(abs(1.0-bmag(ig)*al(2)))

          slbr = (slb1 + slb2)/2.0

          cc(1) = -vn*delt*(-1.0 - slbr)/(slb2-slb1)
          aa(1) = 0.0
          bb(1) = 1.0 - (aa(1) + cc(1))

          il = je
          slb0 = sqrt(abs(1.0-bmag(ig)*al(il-1)))
          slbl = slb0/2.0

          cc(il) = -0.5*vn*delt*(1.0-slbl*slbl)/slbl/slb0
          aa(il) = cc(il)
          bb(il) = 1.0 - (aa(il) + cc(il))

          a1(:je) = aa(:je)
          b1(:je) = bb(:je)
          c1(:je,ilz) = cc(:je)

          a1(je+1:2*je-1) = cc(je-1:1:-1)
          b1(je+1:2*je-1) = bb(je-1:1:-1)
          c1(je+1:2*je-1,ilz) =aa(je-1:1:-1)

          betaa(1,ilz) = 1.0/b1(1)
          do il = 1, 2*je-2
             ql(il+1,ilz) = a1(il+1)*betaa(il,ilz)
             betaa(il+1,ilz) = 1.0/(b1(il+1)-ql(il+1,ilz)*c1(il,ilz))
          end do

          ql(1,ilz) = 0.0
          c1(2*je:,ilz) = 0.0
          betaa(2*je:,ilz) = 0.0
       end if
    end do
  end subroutine init_lorentz

  subroutine init_lorentz_redistribute
    use mp, only: nproc
    use species, only: nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid, jend, ng2
    use gs2_layouts, only: init_lorentz_layouts
    use gs2_layouts, only: g_lo, lz_lo
    use gs2_layouts, only: idx_local, proc_id
    use gs2_layouts, only: gidx2lzidx 
    use redistribute, only: index_list_type, init_redist, delete_list
    implicit none
    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension(3) :: from_low, from_high
    integer, dimension(2) :: to_high
    integer :: to_low
    integer :: ig, isign, iglo, il, ilz
    integer :: n, ip
    logical :: done = .false.

    if (done) return

    call init_lorentz_layouts &
         (ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2)

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
             call gidx2lzidx (ig, isign, g_lo, iglo, lz_lo, ntgrid, jend, il, ilz)
             if (idx_local(g_lo,iglo)) &
                  nn_from(proc_id(lz_lo,ilz)) = nn_from(proc_id(lz_lo,ilz)) + 1
             if (idx_local(lz_lo,ilz)) &
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
       end if
    end do

    ! get local indices of elements distributed to/from other processors
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
             call gidx2lzidx (ig, isign, g_lo, iglo, lz_lo, ntgrid, jend, il, ilz)
!             write(*,*) ig,':',isign,':',iglo,':',il,':',ilz
             if (idx_local(g_lo,iglo)) then
                ip = proc_id(lz_lo,ilz)
                n = nn_from(ip) + 1
                nn_from(ip) = n
                from_list(ip)%first(n) = ig
                from_list(ip)%second(n) = isign
                from_list(ip)%third(n) = iglo
             end if
             if (idx_local(lz_lo,ilz)) then
                ip = proc_id(g_lo,iglo)
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to_list(ip)%first(n) = il
                to_list(ip)%second(n) = ilz
             end if
          end do
       end do
    end do

    from_low (1) = -ntgrid
    from_low (2) = 1
    from_low (3) = g_lo%llim_proc

    to_low = lz_lo%llim_proc

    to_high(1) = max(2*nlambda, 2*ng2+1)
    to_high(2) = lz_lo%ulim_alloc

    from_high(1) = ntgrid
    from_high(2) = 2
    from_high(3) = g_lo%ulim_alloc

    call init_redist (lorentz_map, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)

    done = .true.

  end subroutine init_lorentz_redistribute

  subroutine solfp1 (g, g1, phi, apar, aperp)
    use theta_grid, only: ntgrid
    use run_parameters, only: fphi, faperp
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: g, g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, aperp

    if (collision_model_switch == collision_model_none) return

    call g_adjust (g, phi, aperp, fphi, faperp)

    select case (collision_model_switch)
    case (collision_model_lorentz,collision_model_lorentz_test)
       call solfp_lorentz (g, g1)
    case (collision_model_krook,collision_model_krook_test)
       call solfp_krook (g, g1)
    end select

    call g_adjust (g, phi, aperp, -fphi, -faperp)
  end subroutine solfp1

  subroutine g_adjust (g, phi, aperp, facphi, facaperp)
    use species, only: spec
    use theta_grid, only: ntgrid
    use le_grids, only: anon
    use dist_fn_arrays, only: vperp2, aj0, aj1
    use gs2_layouts, only: g_lo, ik_idx, it_idx, ie_idx, is_idx
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, aperp
    real, intent (in) :: facphi, facaperp

    integer :: iglo, ig, ik, it, ie, is
    complex :: adj

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do ig = -ntgrid, ntgrid
          adj = anon(ie,is)*2.0*vperp2(ig,iglo)*aj1(ig,iglo) &
                  *aperp(ig,it,ik)*facaperp &
               + spec(is)%z*anon(ie,is)*phi(ig,it,ik)*aj0(ig,iglo) &
                  /spec(is)%temp*facphi
          g(ig,1,iglo) = g(ig,1,iglo) + adj
          g(ig,2,iglo) = g(ig,2,iglo) + adj
       end do
    end do
  end subroutine g_adjust

  subroutine solfp_krook (g, g1)
    use species, only: spec, electron_species
    use theta_grid, only: ntgrid
    use le_grids, only: integrate, lintegrate, geint2g, gint2g
    use gs2_layouts, only: g_lo, gint_lo, geint_lo
    use gs2_layouts, only: ik_idx, il_idx, it_idx, is_idx
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, g1

    complex, dimension (-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_alloc) &
         :: g0eint, g1eint
    complex, dimension (-ntgrid:ntgrid,gint_lo%llim_proc:gint_lo%ulim_alloc) &
         :: g1int, g2int
    integer :: iglo, igint, igeint, ig, ik, it, il, is

    call prof_entering ("solfp_krook", "collisions")

    if (conserve_momentum) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) then
             g1(:,:,iglo) = 0.0
          else
             g1(:,:,iglo) = g(:,:,iglo)*sq(:,il,:)
          end if
       end do
       call lintegrate (g1, g1int)
    end if

    if (conserve_number) call integrate (g, g0eint)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       do ig = -ntgrid, ntgrid
          g1(ig,1,iglo) = g(ig,1,iglo)/(1.0 + vnewfe(iglo))
          g1(ig,2,iglo) = g(ig,2,iglo)/(1.0 + vnewfe(iglo))
       end do
    end do

    if (conserve_number) then
       call integrate (g1, g1eint)
       do igeint = geint_lo%llim_proc, geint_lo%ulim_proc
          g0eint(:,igeint) = (g0eint(:,igeint) - g1eint(:,igeint)) &
               *aintnorm(:,igeint)
       end do
       call geint2g (g0eint, g)
       g = g1 + g
    else
       g = g1
    end if

    if (conserve_momentum) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) then
             g1(:,:,iglo) = 0.0
          else
             g1(:,:,iglo) = g(:,:,iglo)*sq(:,il,:)
          end if
       end do
       call lintegrate (g1, g2int)

       do igint = gint_lo%llim_proc, gint_lo%ulim_proc
          g1int(:,igint) = (g1int(:,igint) - g2int(:,igint))/g3int(:,igint)
       end do
       call gint2g (g1int, g1)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) then
          else if (spec(is)%type == electron_species) then
          else
             g(:,:,iglo) = g(:,:,iglo) + sq(:,il,:)*g1(:,:,iglo)
          end if
       end do
    end if

    call prof_leaving ("solfp_krook", "collisions")
  end subroutine solfp_krook

  subroutine solfp_lorentz (g, g1)
    use species, only: spec, electron_species
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, jend, lintegrate, gint2g, ng2
    use gs2_layouts, only: g_lo, gint_lo, lz_lo
    use gs2_layouts, only: ig_idx, ik_idx, il_idx, is_idx
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: gather, scatter
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, g1

    complex, dimension (max(2*nlambda,2*ng2+1)) :: delta
    complex, dimension (-ntgrid:ntgrid,gint_lo%llim_proc:gint_lo%ulim_alloc) &
         :: g1int, g2int
    integer :: iglo, igint, ilz, ig, ik, il, is, je

    call prof_entering ("solfp_lorentz", "collisions")

!    call check_g ('beg', g)

    if (conserve_momentum) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) then
             g1(:,:,iglo) = 0.0
          else
             g1(:,:,iglo) = g(:,:,iglo)*sq(:,il,:)
          end if
       end do
       call lintegrate (g1, g1int)
    end if

    call gather (lorentz_map, g, glz)

!    call check_glz ('beg', glz)
    ! solve for glz row by row
    do ilz = lz_lo%llim_proc, lz_lo%ulim_proc
       ig = ig_idx(lz_lo,ilz)
       ik = ik_idx(lz_lo,ilz)
       is = is_idx(lz_lo,ilz)
       if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) cycle
       je = 2*jend(ig)

       if (je == 0) then
          je = 2*ng2+1
       end if

       glz(je:,ilz) = 0.0

       delta(1) = glz(1,ilz)
       do il = 1, je-1
          delta(il+1) = glz(il+1,ilz) - ql(il+1,ilz)*delta(il)
       end do
       
       glz(je,ilz) = delta(je)*betaa(je,ilz)
       do il = je-1, 1, -1
          glz(il,ilz) = (delta(il) - c1(il,ilz)*glz(il+1,ilz))*betaa(il,ilz)
       end do
       glz(je,ilz) = glz(jend(ig),ilz)
    end do

!    call check_glz ('end', glz)
    call scatter (lorentz_map, glz, g)

    if (conserve_momentum) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) then
             g1(:,:,iglo) = 0.0
          else
             g1(:,:,iglo) = g(:,:,iglo)*sq(:,il,:)
          end if
       end do
       call lintegrate (g1, g2int)

       do igint = gint_lo%llim_proc, gint_lo%ulim_proc
          g1int(:,igint) = (g1int(:,igint) - g2int(:,igint))/g3int(:,igint)
       end do
       call gint2g (g1int, g1)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) then
          else if (spec(is)%type == electron_species) then
          else
             g(:,:,iglo) = g(:,:,iglo) + sq(:,il,:)*g1(:,:,iglo)
          end if
       end do
    end if

!    call check_g ('end', g)

    call prof_leaving ("solfp_lorentz", "collisions")
  end subroutine solfp_lorentz

  subroutine check_g (str, g)

    use gs2_layouts, only: g_lo
    use theta_grid, only: ntgrid
    
    character (3) :: str 
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    integer :: ig, iglo

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do ig = -ntgrid, 0
          if (g(ig,1,iglo) /= g(-ig,2,iglo)) then
             if (abs((g(ig,1,iglo)-g(-ig,2,iglo))/g(ig,1,iglo)) > 1.e-8) then
                write(*,*) str,ig,iglo,g(ig,1,iglo)-g(-ig,2,iglo)
             endif
          endif
       end do
    end do

  end subroutine check_g

  subroutine check_glz (str, glz)

    use gs2_layouts, only: lz_lo, ig_idx, ie_idx
    use le_grids, only: jend
    use theta_grid, only: ntgrid
    
    character (3) :: str 
    complex, dimension (:,lz_lo%llim_proc:), intent (in) :: glz
    integer :: il, ilz, ig, ie, je, ilzp, ilp

    do ilz = lz_lo%llim_proc, lz_lo%ulim_proc
       ig = ig_idx (lz_lo, ilz)
       ie = ie_idx (lz_lo, ilz)
       je = jend(ig)
       ilzp = -ig + ntgrid + (ie-1)*(2*ntgrid+1)
       do il = 1, je
          ilp = 2*je+1-il
          if (glz(il,ilz) /= glz(ilp,ilzp)) then
             if (abs((glz(il,ilz)-glz(ilp,ilzp))/glz(il,ilz)) > 1.e-8) then
                write(*,*) str,il,ilz,glz(il,ilz),glz(ilp,ilzp)
             endif
          endif
       end do
    end do

  end subroutine check_glz

  subroutine reset_init
    
    initialized = .false.

  end subroutine reset_init

end module collisions




