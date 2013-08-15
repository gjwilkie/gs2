module egrid

! By Tomo Tatsuno, Aug 2005
! Improved accuracy and speed and maximum number of energy grid points
!

  implicit none

  public :: setvgrid, init_egrid
  public :: setvgrid_genquad
  public :: zeroes, x0, zeroes_maxwell, x0_maxwell
  !> Expose this for unit testing
  public :: energy_grid

  

  private

  real :: x0_maxwell
  real,dimension(:), allocatable :: x0
  real, dimension(:,:), allocatable, save :: zeroes
  real, dimension(:,:), allocatable, save :: energy_grid
  real, dimension(:), allocatable, save :: zeroes_maxwell

contains

  subroutine init_egrid (negrid)
    use species, only: nspec
    
    integer, intent (in) :: negrid


    if (.not. allocated(zeroes)) then
       allocate (energy_grid(negrid, nspec))
       allocate (zeroes(negrid-1, nspec)) ; zeroes = 0.
       allocate (zeroes_maxwell(negrid-1)) ; zeroes = 0.
       allocate (x0(nspec))
    end if

  end subroutine init_egrid

  subroutine setvgrid (vcut, negrid, epts, wgts, nesub)

    use general_f0, only: calculate_f0_arrays, f0_values
    use constants, only: pi => dpi
    use general_f0, only: energy_0_general_f0 => energy_0
    use gauss_quad, only: get_legendre_grids_from_cheb, get_laguerre_grids
    use species, only: nspec, spec, alpha_species

    implicit none
    
    integer, intent (in) :: negrid
    real, intent (in) :: vcut
    integer, intent (in) :: nesub
    real, dimension(:,:), intent (out) :: epts
    real, dimension(:,:), intent (out) :: wgts
    real :: energy_0, vcut_local
    integer :: is,ie

    call init_egrid (negrid)

    do is = 1,nspec

      if (spec(is)%type .eq. alpha_species) then
        vcut_local = 1.0
        energy_0 = energy_0_general_f0
        !energy_0 = 0.2
      else
        vcut_local = vcut
        energy_0 = 0.0
      end if
      
      !write (*,*) 'For species ', is, ' energy_0 is ', energy_0
      !write (*,*) 'nesub is ', nesub

    
      ! get grid points in v up to vcut (epts is not E yet)
      call get_legendre_grids_from_cheb (energy_0**0.5, vcut_local, epts(:nesub,is), wgts(:nesub, is))

      ! change from v to E
      epts(:nesub,is) = epts(:nesub,is)**2

      ! absorb exponential and volume element in weights
      !wgts(:nesub) = wgts(:nesub)*epts(:nesub)*exp(-epts(:nesub))/sqrt(pi)
      ! No longer absorb maxwellian... allow for arbitrary f0. EGH/GW
      ! See eq. 4.12 of M. Barnes's thesis
      !wgts(:nesub, is) = wgts(:nesub, is)*epts(:nesub,is)*2.0*pi
      wgts(:nesub, is) = wgts(:nesub, is)*epts(:nesub,is)*pi

      if (negrid > nesub) then

         ! get grid points in y = E - vcut**2 (epts not E yet)
         call get_laguerre_grids (epts(nesub+1:,is), wgts(nesub+1:, is))

         ! change from y to E
         epts(nesub+1:,is) = epts(nesub+1:,is) + vcut_local**2

         ! absort exponential and volume element in weights
         ! See eq. 4.13 of M. Barnes's thesis
      !   wgts(nesub+1:) = wgts(nesub+1:)*0.5*sqrt(epts(nesub+1:)/pi)*exp(-vcut**2)

         ! No longer absorb maxwellian... allow for arbitrary f0. EGH/GW
         ! Note that here this means adding an exponential e^y
         ! See eq. 4.13 of M. Barnes's thesis
         wgts(nesub+1:, is) = wgts(nesub+1:, is)*exp(epts(nesub+1:,is))*pi*0.5*sqrt(epts(nesub+1:,is))*exp(-vcut_local**2)

      end if

    end do


    !write (*,*) 'Weights in le_grids are', wgts
    !write (*,*) 'epts', epts
    call calculate_f0_arrays(epts, wgts, vcut)
    !do ie = 1,negrid
    !   write(*,*) ie, sqrt(epts(ie,3)), wgts(ie,3), f0_values(ie,3)
    !end do
    !do is = 1,nspec
      !wgts(:,is) = wgts(:,1)
    !end do 
    !do is = 1,nspec
      !wgts(:, is) = wgts(:, is) * f0_values(:, is)
    !end do

    energy_grid = epts
    zeroes(:,:) = sqrt(epts(:negrid-1,:))
    x0(:) = sqrt(epts(negrid,:))

  end subroutine setvgrid

  subroutine setvgrid_genquad (ne_int,vcut, negrid, epts, wgts, nesub)

    use general_f0, only: calculate_f0_arrays, f0_values
    use constants, only: pi => dpi
    use general_f0, only: energy_0_general_f0 => energy_0
    use general_quad, only: get_general_weights_from_grid
    use gauss_quad, only: get_legendre_grids_from_cheb, get_laguerre_grids
    use species, only: nspec, spec, alpha_species

    implicit none
    
    integer, intent (in) :: negrid, ne_int
    real, intent (in) :: vcut
    integer, intent (in) :: nesub
    real, dimension(:,:), intent (out) :: epts
    real, dimension(:,:), intent (out) :: wgts
    real :: energy_0, vcut_local
    real, dimension(1:ne_int):: vpts_temp, wgts_temp, F0_temp
    integer :: is

    call init_egrid (negrid)

    do is = 1,nspec

      if (spec(is)%type .eq. alpha_species) then
        vcut_local = 1.0
        energy_0 = energy_0_general_f0
        !energy_0 = 0.2
      else
        vcut_local = vcut
        energy_0 = 0.0
      end if

      call get_legendre_grids_from_cheb (energy_0**0.5, vcut_local, vpts_temp(1:ne_int), wgts_temp(1:ne_int))

      F0_temp(1:ne_int) = exp(-vpts_temp(1:ne_int)**2) 

      call get_general_weights_from_grid (ne_int, vpts_temp, wgts_temp, F0_temp, energy_0**0.5, vcut_local, &
                                          nesub,epts(1:nesub,is),wgts(1:nesub,is))

!      ! get grid points in v up to vcut (epts is not E yet)
!      call get_legendre_grids_from_cheb (energy_0**0.5, vcut_local, epts(1:nesub,is), wgts(1:nesub, is))

      ! change from v to E
      epts(:nesub,is) = epts(:nesub,is)**2

      ! absorb exponential and volume element in weights
      !wgts(:nesub) = wgts(:nesub)*epts(:nesub)*exp(-epts(:nesub))/sqrt(pi)
      ! No longer absorb maxwellian... allow for arbitrary f0. EGH/GW
      ! See eq. 4.12 of M. Barnes's thesis
      !wgts(:nesub, is) = wgts(:nesub, is)*epts(:nesub,is)*2.0*pi
      wgts(:nesub, is) = wgts(:nesub, is)*epts(:nesub,is)*pi*exp(epts(:nesub,is))

      if (negrid > nesub) then

         ! get grid points in y = E - vcut**2 (epts not E yet)
         call get_laguerre_grids (epts(nesub+1:,is), wgts(nesub+1:, is))

         ! change from y to E
         epts(nesub+1:,is) = epts(nesub+1:,is) + vcut_local**2

         ! absort exponential and volume element in weights
         ! See eq. 4.13 of M. Barnes's thesis
      !   wgts(nesub+1:) = wgts(nesub+1:)*0.5*sqrt(epts(nesub+1:)/pi)*exp(-vcut**2)

         ! No longer absorb maxwellian... allow for arbitrary f0. EGH/GW
         ! Note that here this means adding an exponential e^y
         ! See eq. 4.13 of M. Barnes's thesis
         wgts(nesub+1:, is) = wgts(nesub+1:, is)*exp(epts(nesub+1:,is))*pi*0.5*sqrt(epts(nesub+1:,is))*exp(-vcut_local**2)

      end if

    end do


!    write (*,*) 'Weights in le_grids are', wgts
!    write (*,*) 'epts', epts
    call calculate_f0_arrays(epts, wgts, vcut)
    !do is = 1,nspec
      !wgts(:,is) = wgts(:,1)
    !end do 
    !do is = 1,nspec
      !wgts(:, is) = wgts(:, is) * f0_values(:, is)
    !end do

    energy_grid = epts
    zeroes(:,:) = sqrt(epts(:negrid-1,:))
    x0(:) = sqrt(epts(negrid,:))

  end subroutine setvgrid_genquad

end module egrid

module le_grids
  
  use redistribute, only: redist_type

  implicit none

  public :: init_le_grids, finish_le_grids
  public :: read_parameters, wnml_le_grids
  public :: integrate_species, write_mpdist, write_mpdist_le
  public :: energy, energy_maxwell, speed_maxwell
  public :: anon, al, delal, jend, forbid, dele, wl, w, w_maxwell
  public :: negrid, nlambda, ng2, nxi, lmax, integrate_moment, nesub
  public :: xloc, sgn, ixi_to_il, ixi_to_isgn, speed, xi
  public :: xx, nterp, new_trap_int, vcut
  public :: init_weights, legendre_transform, lagrange_interp, lagrange_coefs
  public :: eint_error, lint_error, trap_error, wdim
  public :: integrate_kysum, integrate_volume
  public :: get_flux_vs_theta_vs_vpa
  public :: lambda_map, energy_map, g2le, init_map
  public :: genquad

  !> Unit tests
  public :: le_grids_unit_test_init_le_grids
  public :: le_grids_unit_test_integrate_species

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
  real, dimension (:,:), allocatable, save :: wlerr, xloc ! mbmark
  real, dimension (:,:,:), allocatable, save :: werr ! mbmark
  real, dimension (:,:,:), allocatable, save :: wlterr
  real, dimension (:,:,:,:), allocatable, save :: lgrnge
  real, dimension (:,:,:), allocatable, save :: wlmod
  real, dimension (:,:,:), allocatable, save :: wtmod
  real, dimension (:,:,:), allocatable, save :: wmod
  real, dimension (:,:,:), allocatable, save :: lpe ! now depends on species
  real, dimension (:,:), allocatable, save :: lpl
  real, dimension (:,:,:,:), allocatable, save :: lpt

  real, dimension (:), allocatable :: anon, energy_maxwell, speed_maxwell ! (negrid)
  real, dimension (:,:), allocatable :: energy, speed, dele ! (negrid,nspec)
  real, dimension (:,:), allocatable :: w ! (negrid,nspec)
  ! Maxwellian weights, for the collision operator...
  ! This is a hack because collisions only work for Maxwellian
  ! distributions currently. 
  real, dimension (:), allocatable :: w_maxwell ! (negrid)
  real, dimension (:), allocatable :: al, delal ! (nlambda)
  real, dimension (:,:), allocatable :: wl ! (nlambda,-ntgrid:ntgrid)
  integer, dimension (:), allocatable :: jend ! (-ntgrid:ntgrid)
  logical, dimension (:,:), allocatable :: forbid ! (-ntgrid:ntgrid,nlambda)

  real, dimension (:,:), allocatable :: xi
  integer, dimension (:,:), allocatable :: ixi_to_il, ixi_to_isgn
  integer, dimension (2) :: sgn

  integer :: wdim
  integer :: lint_lo, lint_hi, eint_lo, eint_hi
  integer :: geint2g_lo, geint2g_hi
  complex, dimension (:,:), allocatable :: integration_work
  ! (-ntgrid:ntgrid, -*- processor-dependent -*-)
  logical :: exist

 ! knobs
  integer :: ngauss, negrid, nesuper, nesub
  real :: bouncefuzz, vcut
  logical :: genquad = .false.
  integer :: ne_int_genquad = 200

  integer :: nlambda, ng2, lmax, nxi
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
  logical :: leinit = .false.
  logical :: lzinit = .false.
  logical :: einit = .false.

  integer :: nmax = 500
  integer :: nterp = 100

  real :: wgt_fac = 10.0

  type (redist_type), save :: lambda_map
  type (redist_type), save :: energy_map
  type (redist_type), save :: g2le

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
    end if
    !write (*,*) 'negrid', negrid, 'proc0', proc0
    call broadcast_parameters
    !write (*,*) 'negrid', negrid, 'proc0', proc0
    call set_vgrid
    if (proc0) then
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
             write(*,*) energy(ie,:)
          end do
       end if
       call finish_mp
       stop
    endif
    
  end subroutine init_le_grids

  function le_grids_unit_test_init_le_grids(sizes, energy_results, err)
    use unit_tests
    use egrid
    use species, only: nspec
    use mp, only: proc0
    integer, dimension(:), intent(in) :: sizes
    real, dimension(:,:,:), intent(in) :: energy_results
    real, intent(in) :: err
    logical :: le_grids_unit_test_init_le_grids
    logical :: tr ! Test result
    logical :: accelerated_x, accelerated_v
    integer :: i
    character(2) :: istr

    tr = .true. 
    
    call init_le_grids(accelerated_x, accelerated_v)
    !call broadcast_results
     !write (*,*)  'zeros', zeroes

    call announce_check('Size of energy arrays')
    tr = tr .and. agrees_with(size(energy), sizes(1)) 
    tr = tr .and. agrees_with(size(w), sizes(1)) 
    call process_check(tr, 'Size of energy arrays')

    ! Energy grids only get calulcated on proc0
    if (proc0) then
      do i = 1,nspec
        write(istr, '(I2)') i
        call announce_check('values of the energy grid for species '//istr)
        tr = tr .and. agrees_with(energy_grid(:,i), energy_results(:,i,1), err)
        call process_check(tr, 'values of the energy grid for species '//istr)
      end do
    end if
    do i = 1,nspec
      write(istr, '(I2)') i
      call announce_check('values of energy weights for species '//istr)
      tr = tr .and. agrees_with(w(:,i), energy_results(:,i,2), err)
      call process_check(tr, 'values of energy weights for species '//istr)
    end do
    
    le_grids_unit_test_init_le_grids = tr
  end function le_grids_unit_test_init_le_grids

  subroutine broadcast_parameters
    use mp, only: proc0, broadcast
    !write (*,*) 'Broadcasting ', negrid
    call broadcast (ngauss)
    call broadcast (negrid)
    call broadcast (nesuper)
    call broadcast (nesub)
    call broadcast (vcut)
    call broadcast (bouncefuzz)
    call broadcast (test)
    call broadcast (nmax)
    call broadcast (trapped_particles)
    call broadcast (wgt_fac)
    call broadcast (new_trap_int)
    call broadcast (nterp)
    call broadcast (genquad)
    call broadcast (ne_int_genquad)
    !write (*,*) 'Broadcasting ', negrid
  end subroutine broadcast_parameters

  subroutine broadcast_results
    use mp, only: proc0, broadcast
    use species, only: nspec
    use egrid, only: zeroes, x0, init_egrid, zeroes_maxwell, x0_maxwell
    use theta_grid, only: ntgrid
    use general_f0, only: broadcast_arrays_gen_f0 => broadcast_arrays

    implicit none
    integer :: ig, il, is, ie, ipt, isgn, tsize

    tsize = 2*nterp-1



    call broadcast (lmax)
    call broadcast (ng2)
    call broadcast (nlambda)
    call broadcast (nxi)

    if (.not. proc0) then
       !allocate (energy(negrid,nspec), w(negrid, nspec), anon(negrid), speed(negrid,nspec))
       allocate (anon(negrid), speed(negrid,nspec))
       !if (.not. allocated(zeroes)) &
         !allocate (zeroes(negrid-1,nspec), zeroes_maxwell(negrid-1))
       !allocate (x0(nspec))
       allocate (w_maxwell(negrid))
       allocate (energy_maxwell(negrid))
       allocate (speed_maxwell(negrid))
       allocate (dele(negrid, nspec))
       allocate (al(nlambda), delal(nlambda))
       allocate (wl(-ntgrid:ntgrid,nlambda))
       allocate (jend(-ntgrid:ntgrid))
       allocate (forbid(-ntgrid:ntgrid,nlambda))
       allocate (xx(ng2))
       allocate (lgrnge(-ntgrid:ntgrid,nlambda,tsize,2))
       allocate (xloc(-ntgrid:ntgrid,tsize))
       allocate (xi(-ntgrid:ntgrid, 2*nlambda))
       allocate (ixi_to_il(-ntgrid:ntgrid, 2*nlambda))
       allocate (ixi_to_isgn(-ntgrid:ntgrid, 2*nlambda))
    end if

    !call init_egrid (negrid)
    call broadcast (xx)
    !call broadcast (x0)
    !call broadcast (x0_maxwell)
    !call broadcast (zeroes)
    !call broadcast (zeroes_maxwell)

    call broadcast (al)
    call broadcast (delal)
    call broadcast (jend)
 
    !call broadcast (energy)
    call broadcast (energy_maxwell)
    call broadcast (speed)
    call broadcast (speed_maxwell)
    call broadcast (dele)
    !call broadcast (w)
    call broadcast (w_maxwell)
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

    do ig = -ntgrid, ntgrid
       call broadcast (xi(ig,:))
       call broadcast (ixi_to_il(ig,:))
       call broadcast (ixi_to_isgn(ig,:))
    end do
    call broadcast (sgn)

    !call broadcast_arrays_gen_f0

  end subroutine broadcast_results

  subroutine read_parameters
    use species, only: spec, has_slowing_down_species
    use file_utils, only: input_unit, error_unit, input_unit_exist
    implicit none
    integer :: ierr, in_file
    namelist /le_grids_knobs/ ngauss, negrid, bouncefuzz, &
         nesuper, nesub, test, trapped_particles, &
         nmax, wgt_fac, new_trap_int, nterp, vcut, &
         genquad, ne_int_genquad

    nesub = 8
    nesuper = 2
    ngauss = 5
    negrid = -10
    vcut = 2.5
    bouncefuzz = 1e-5
    genquad = .false.
    ne_int_genquad = 200

    in_file=input_unit_exist("le_grids_knobs", exist)
    if (exist) read (unit=in_file, nml=le_grids_knobs)

! user can choose not to set negrid (preferred for old scheme)
    if (negrid == -10) then
       negrid = nesub + nesuper

! If user chose negrid, assume nesuper makes sense and check nesub if necessary
    else 
       nesuper = min(negrid/10+1, 4)
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
    use general_f0, only: f0_values
    use egrid, only: zeroes
    use constants, only: pi => dpi
    use species, only: nspec

    implicit none

    real, dimension (:), allocatable :: modzeroes  ! (negrid-2)
    real, dimension (:), allocatable :: werrtmp  ! (negrid-2)
    real, dimension (:), allocatable :: lmodzeroes, wlerrtmp ! (ng2-1)
    integer :: ipt, ndiv, divmax
    logical :: eflag = .false.

    integer :: ie, il, is

    allocate(lmodzeroes(ng2-1), wlerrtmp(ng2-1))
    allocate(wlerr(ng2,ng2))

    wlerr = 0.0; lmodzeroes = 0.0; wlerrtmp = 0.0

    wdim = nesub
    allocate(modzeroes(nesub-1), werrtmp(nesub-1))
    allocate(werr(negrid-1,nesub,nspec))
    
    werr = 0.0 ; modzeroes = 0.0 ; werrtmp = 0.0
    
    ! loop to obtain weights for energy grid points.  negrid-1 sets
    ! of weights are needed because we want to compute integrals
    ! for negrid-1 sets of energy points (corresponding to negrid-1
    ! points that we can choose to drop from the guassian quadrature)
    
    do ipt=1,nesub
       do is = 1,nspec
       
         ! drops the point corresponding to ipt from the energy grid
         
         if (ipt /= 1) modzeroes(:ipt-1) = zeroes(:ipt-1,is)
         if (ipt /= nesub) modzeroes(ipt:nesub-1) = zeroes(ipt+1:nesub,is)
         
         ! get weights for energy grid points
         
         call get_weights (nmax,0.0,vcut,modzeroes,werrtmp,ndiv,divmax,eflag)
         
         ! a zero is left in the position corresponding to the dropped point
         
         if (ipt /= 1) werr(:ipt-1,ipt,is) = werrtmp(:ipt-1)
         if (ipt /= nesub) werr(ipt+1:nesub,ipt,is) = werrtmp(ipt:nesub-1)
         werr(nesub+1:,ipt,is) = w(nesub+1:negrid-1,is)
         
         ! absorbing volume element into weights
         !werr(:nesub,ipt,is) = werr(:nesub,ipt)*energy(:nesub)*exp(-energy(:nesub))/sqrt(pi)
         werr(:nesub,ipt,is) = &
              werr(:nesub,ipt,is)*energy(:nesub,is)*2.0*pi*f0_values(:nesub,is)
       end do
       
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
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_allreduce

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:), intent (out) :: total
    integer :: is, il, ie, ik, it, iglo

    !Ensure array is zero to begin
    total = 0.

    !Performed integral (weighted sum) over local velocity space and species
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       !Convert from iglo to the separate indices
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)

       !Sum up weighted g
       total(:, it, ik) = total(:, it, ik) + weights(is)*w(ie,is)*wl(:,il)*(g(:,1,iglo)+g(:,2,iglo))
    end do

    !Reduce sum across all procs to make integral over all velocity space and species
    call sum_allreduce (total) 

  end subroutine integrate_species

  function le_grids_unit_test_integrate_species(g, weights, sizes, rslt, err)
    use unit_tests
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    integer, dimension (:), intent (in) :: sizes
    complex, dimension (:,:,:), allocatable :: total
    complex, dimension (-ntgrid:,:,:) :: rslt
    real, intent(in) :: err
    logical :: le_grids_unit_test_integrate_species
    logical :: tr

    tr = .true.

    call announce_check('Size of naky')
    tr = tr .and. agrees_with(naky, sizes(1))
    call process_check(tr, 'Size of naky')
    call announce_check('Size of ntheta0')
    tr = tr .and. agrees_with(ntheta0, sizes(2))
    call process_check(tr, 'Size of ntheta0')
    call announce_check('Size of ntgrid')
    tr = tr .and. agrees_with(ntgrid, sizes(3))
    call process_check(tr, 'Size of ntgrid')

    allocate(total(-ntgrid:ntgrid,naky,ntheta0))

    total = cmplx(0.0,0.0)


    call integrate_species(g, weights, total)

    call announce_check('total')
    tr = tr .and. agrees_with(total(:,1,1), rslt(:,1,1), err)
    call process_check(tr, 'total ')


    deallocate(total)

    le_grids_unit_test_integrate_species = tr



  end function le_grids_unit_test_integrate_species

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
    real :: fac, ulim
    integer :: is, il, ie, ik, it, iglo, ig, i, j, im, ntrap, k
    integer, save :: lpesize

    real, dimension (:), allocatable :: nodes
    real, dimension (:,:), allocatable :: lpltmp, lpttmp
!    real, dimension (:,:), allocatable, save :: lpe, lpl
!    real, dimension (:,:,:,:), allocatable, save :: lpt

    if (.not. allocated(lpl)) then
       allocate(lpltmp(ng2,0:ng2-1))
       allocate(lpl(nlambda,0:ng2-1))

       lpesize = nesub
       allocate(lpe(negrid,0:lpesize-1,nspec)) ; lpe = 0.0
       
       ! get value of first nesub legendre polynomials
       ! at each of the grid points on (0,vcut)
       do is = 1,nspec
        call legendre_polynomials(0.0,vcut,zeroes(:lpesize,is),lpe(:lpesize,:,is))
       end do
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
!Can we remove this?
! TEMP FOR TESTING -- MAB
!                nodes = nodes + sqrt(1.0-bmag(ig)/bmax)
!                ulim = 2.*sqrt(1.0-bmag(ig)/bmax)
                ulim = sqrt(1.0-bmag(ig)/bmax)
                call legendre_polynomials (-ulim,ulim,nodes,lpttmp)
                lpt(ng2+1:jend(ig),0:2*(ntrap-1),ig,2) = lpttmp(1:ntrap,:)
                lpt(ng2+1:jend(ig)-1,0:2*(ntrap-1),ig,1) = lpttmp(2*ntrap-1:ntrap+1:-1,:)
!Can we remove this?
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
    end if

    ! carry out legendre transform to get coefficients of
    ! legendre polynomial expansion of g
    totfac = 0. ; tote = 0. ; totl = 0.
    if (present(tott)) tott = 0.

    !Loop over all indices, note this loop is optimal only for layout 'xyles' (at least in terms of
    !g memory access)
    do is = 1, nspec
       do ie = 1, negrid
          do il = 1, nlambda
             do ik = 1, naky       !Swapped ik and it loop order.
                do it = 1, ntheta0
                   iglo = idx (g_lo, ik, it, il, ie, is)
                   if (idx_local (g_lo, iglo)) then
                      do ig=-ntgrid,ntgrid
                         totfac = w(ie,is)*wl(ig,il)*(g(ig,1,iglo)+g(ig,2,iglo))
                         do im=0,lpesize-1
                            tote(im, ig, it, ik, is) = tote(im, ig, it, ik, is) + totfac*lpe(ie,im,is)*(2*im+1)
                         end do
                         do im=0,ng2-1
                            totl(im, ig, it, ik, is) = totl(im, ig, it, ik, is) + totfac*lpl(il,im)*(2*im+1)
                         end do
                         if (present(tott)) then
                            do im=0,2*(jend(ig)-ng2-1)
                               tott(im, ig, it, ik, is) = tott(im, ig, it, ik, is) + &
                                    w(ie,is)*wl(ig,il)*(lpt(il,im,ig,1)*g(ig,1,iglo)+lpt(il,im,ig,2)*g(ig,2,iglo))*(2*im+1)
                            end do
                         end if
                      end do
                   end if
                end do
             end do
          end do
       end do
    end do

    !Do we really need this if?
    if (nproc > 1) then
       !Now complete velocity integral, bringing back results to proc0
       call sum_reduce (tote, 0)
       call sum_reduce (totl, 0)
       if (present(tott)) call sum_reduce (tott, 0)
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
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, idx, idx_local
    use mp, only: sum_reduce, proc0, sum_allreduce, nproc

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:,:,:), intent (out) :: poly
    integer, optional, intent(in) :: all
    integer, intent(in) :: istep
    real :: ypt
    integer :: il, ie, ik, it, iglo, ig

    !Initialise to zero
    poly = 0.
    do ie = 1, negrid
       do il = ng2+1, nlambda
          do it = 1, ntheta0
             do ik = 1, naky
                iglo = idx (g_lo, ik, it, il, ie, 1)
                if (idx_local (g_lo, iglo)) then
                   do ig = -ntgrid, ntgrid
                      !Do we need this if?
                      if (.not. forbid(ig,il)) then
                         ypt=sqrt(max(1.0 - bmag(ig)*al(il),0.0))                      
                      else
                         ypt=0.0
                      end if
                      poly(ig, it, ik, ie, :) = poly(ig, it, ik, ie, :) + &
                           lgrnge(ig,il,:,1)*cos(0.1*istep*ypt+1.0) + &
                           lgrnge(ig,il,:,2)*cos(-0.1*istep*ypt+1.0)
                   end do
                end if
             end do
          end do
       end do
    end do

    !Do we really need this if?
    if (nproc > 1) then
       !Reduce integral across procs
       if (present(all)) then
          call sum_allreduce (poly)
       else
          call sum_reduce (poly, 0)
       end if
    end if

  end subroutine lagrange_interp

!  subroutine integrate_moment (g, total, all)
  subroutine integrate_moment_c34 (g, total, all)
! returns results to PE 0 [or to all processors if 'all' is present in input arg list]
! NOTE: Takes f = f(x, y, z, sigma, lambda, E, species) and returns int f, where the integral
! is over all velocity space
! TT>
    use gs2_layouts, only: g_lo, is_idx, ik_idx, it_idx, ie_idx, il_idx
! <TT
    use theta_grid, only: ntgrid
    use mp, only: sum_reduce, proc0, sum_allreduce, nproc

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    integer, optional, intent(in) :: all
    integer :: is, il, ie, ik, it, iglo

    !Ensure array is zero
    total = 0.

    !Integrate over local velocity space
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)

       !Perform local sum
       total(:, it, ik, is) = total(:, it, ik, is) + &
            w(ie,is)*wl(:,il)*(g(:,1,iglo)+g(:,2,iglo))

    end do

    !Not sure that we really need to limit this to nproc>1 as if
    !we run with 1 proc MPI calls should still work ok
    if (nproc > 1) then     
       if (present(all)) then
          !Complete integral over distributed velocity space and ensure all procs know the result
          call sum_allreduce (total)
       else
          !Complete integral over distributed velocity space but only proc0 knows the answer
          call sum_reduce (total, 0)
       end if
    end if

  end subroutine integrate_moment_c34

!  subroutine integrate_moment (g, total, all)
  subroutine integrate_moment_r33 (g, total, all)
! returns results to PE 0 [or to all processors if 'all' is present in input arg list]
! NOTE: Takes f = f(y, z, sigma, lambda, E, species) and returns int f, where the integral
! is over all velocity space
    use mp, only: nproc
    use gs2_layouts, only: p_lo, is_idx, ik_idx, ie_idx, il_idx
    use mp, only: sum_reduce, proc0, sum_allreduce
    use theta_grid, only: ntgrid

    implicit none

    real, dimension (-ntgrid:,:,p_lo%llim_proc:), intent (in) :: g
    real, dimension (-ntgrid:,:,:), intent (out) :: total
    integer, optional, intent(in) :: all
    integer :: is, il, ie, ik, iplo

    !Ensure zero to start
    total = 0.

    !Do local velocity space integral
    do iplo = p_lo%llim_proc, p_lo%ulim_proc
       ik = ik_idx(p_lo,iplo)
       ie = ie_idx(p_lo,iplo)
       is = is_idx(p_lo,iplo)
       il = il_idx(p_lo,iplo)

       total(:, ik, is) = total(:, ik, is) + &
            w(ie,is)*wl(:,il)*(g(:,1,iplo)+g(:,2,iplo))

    end do

    !Do we really need this if?
    if (nproc > 1) then
       !Complete distributed integral
       if (present(all)) then
          !Return result to all procs
          call sum_allreduce (total)
       else
          !Only proc0 knows the result
          call sum_reduce (total, 0)
       end if
    end if

  end subroutine integrate_moment_r33

  subroutine integrate_moment_lec (lo, g, total)
!Perform an integral over velocity space whilst in the LE_LAYOUT in 
!which we have ensured that all of velocity space is local. As such
!we don't need any calls to MPI reduction routines. Note that this means
!the processors for different distributed spatial points (x,y) don't know
!the results at other points.
    use layouts_type, only: le_layout_type
    use gs2_layouts, only: ig_idx, it_idx, ik_idx, is_idx

    implicit none

    type (le_layout_type), intent (in) :: lo
    complex, dimension (:,:,lo%llim_proc:), intent (in) :: g
    complex, dimension (lo%llim_proc:), intent (out) :: total
    integer :: ixi, ie, il, ile, ig,is

    total = 0.0
    do ile = lo%llim_proc, lo%ulim_proc
       ig = ig_idx (lo,ile)
       is = is_idx (lo,ile)
       do ie=1, negrid
          do ixi=1, nxi
             il = ixi_to_il(ig,ixi)
             total(ile) = total(ile) + w(ie,is) * wl(ig,il) * g(ixi,ie,ile)
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
    use species, only: nspec
    use kt_grids, only: aky
    use constants, only: zi
    use gs2_layouts, only: is_idx, ik_idx, ie_idx, il_idx, p_lo
    use mp, only: sum_reduce, proc0, sum_allreduce, nproc

    implicit none

    complex, dimension (p_lo%llim_proc:), intent (in) :: g
    integer, intent (in) :: ig
    complex, dimension (:), intent (out) :: total
    integer, optional, intent(in) :: all

    complex, dimension (negrid,nlambda,nspec) :: gksum
    integer :: is, il, ie, ik, it, iplo

    !Initialise both arrays to zero
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

       total(is) = total(is) + w(ie,is)*wl(ig,il)*gksum(ie,il,is)
    end do

    !Do we really need this if?
    if (nproc > 1) then
       if (present(all)) then
          call sum_allreduce (total)
       else
          call sum_reduce (total, 0)
       end if
    end if

  end subroutine integrate_kysum

  subroutine lint_error (g, weights, total)
    use theta_grid, only: ntgrid, bmag, bmax
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_allreduce, proc0, broadcast

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    integer :: is, il, ie, ik, it, iglo, ipt

    !If the weights array hasn't been filled in then do so now
    if (.not. allocated (wlmod)) then
       !Allocate array, don't need to initialise as below loop ensures
       !all elements are assigned a value
       allocate (wlmod(-ntgrid:ntgrid,nlambda,ng2))

       if (proc0) then
          do ipt = 1, ng2
             do il = 1, ng2
                wlmod(:,il,ipt) = wlerr(il,ipt)*2.0*sqrt((bmag(:)/bmax) &
                     *((1.0/bmax-al(il))/(1.0/bmag(:)-al(il))))
             end do
             !If we have trapped particles use the precalculated weights
             !in wlmod as above is only for passing particles
             if (nlambda > ng2) wlmod(:,ng2+1:,ipt) = wl(:,ng2+1:)
          end do
       end if

       !Now send the calculated value from proc0 to all other procs
       !We could just do the above calculations on all procs?
       call broadcast (wlmod)
    end if

    !Initialise to zero
    total = 0.

    !For each (passing) lambda point do velocity space integral
    do ipt=1,ng2
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)

          total(:, it, ik, ipt) = total(:, it, ik, ipt) + weights(is)*w(ie,is)*wlmod(:,il,ipt)*(g(:,1,iglo)+g(:,2,iglo))
       end do
    end do

    !Moved this outside of the ipt loop above
    call sum_allreduce (total) 

  end subroutine lint_error

  subroutine trap_error (g, weights, total)
    use theta_grid, only: ntgrid, bmag, bmax
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_allreduce, proc0, broadcast

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    integer :: is, il, ie, ik, it, iglo, ipt, ntrap

    !How many trapped pitch angles are there?
    ntrap = nlambda - ng2

    !If weights not calculated yet do so now
    if (.not. allocated(wtmod)) then
       !Allocate array, don't need to initialise as below loops
       !ensure every element is assigned a value
       allocate (wtmod(-ntgrid:ntgrid,nlambda,ntrap))
          
       if (proc0) then
          do ipt=1,ntrap
             wtmod(:,:ng2,ipt) = wl(:,:ng2)
          end do

!Left below comments, but are we done testing this now?
! next line only to be used when testing!!!!
!          wtmod(:,:ng2,:) = 0.

          wtmod(:,ng2+1:,:) = wlterr(:,ng2+1:,:)
       endif

       !Send from proc0 to all others | We could just do the above calculations on all procs?
       call broadcast (wtmod)
    end if


    !Initialise to zero
    total = 0.

    !Loop over number of trapped points
    do ipt=1,ntrap
       !Do local velocity integral
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)

          total(:, it, ik, ipt) = total(:, it, ik, ipt) + weights(is)*w(ie,is)*wtmod(:,il,ipt)*(g(:,1,iglo)+g(:,2,iglo))
       end do
    end do

    !Moved this out of ipt loop above
    call sum_allreduce (total) 

  end subroutine trap_error

  subroutine eint_error (g, weights, total)
    use theta_grid, only: ntgrid
    use species, only: nspec
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: is_idx, ik_idx, it_idx, ie_idx, il_idx
    use mp, only: sum_allreduce, proc0, broadcast

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    integer :: is, il, ie, ik, it, iglo, ipt

    !If we don't have the weights then calculate them now
    if (.not. allocated(wmod)) then
       !Allocate array, don't initialise as we fill in all values below
       allocate (wmod(negrid,wdim,nspec))

       if (proc0) then
          wmod(:negrid-1,:,:) = werr(:,:,:)
          do is = 1,nspec
            wmod(negrid,:,is) = w(negrid,is)
          end do   
       end if

       !send from proc0 to everywhere else
       call broadcast(wmod)
    end if

    !Initialise to zero
    total=0.

    !Do velocity space integral for each ipt (for all energy grid points)
    do ipt=1,wdim
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)

          total(:, it, ik, ipt) = total(:, it, ik, ipt) + weights(is)*wmod(ie,ipt,is)*wl(:,il)*(g(:,1,iglo)+g(:,2,iglo))
       end do
    end do

    !Moved this out of the above loop over ipt
    call sum_allreduce (total) 

  end subroutine eint_error

  subroutine set_vgrid
    use species, only: nspec, init_species
    use egrid, only: setvgrid, setvgrid_genquad
    call init_species

    !write (*,*) 'negrid is', negrid, 'nspec is ', nspec

    allocate (energy(negrid,nspec), w(negrid,nspec))
    if (genquad) then
       call setvgrid_genquad (ne_int_genquad, vcut, negrid, energy, w, nesub)
    else
       call setvgrid (vcut, negrid, energy, w, nesub)
    end if
  end subroutine set_vgrid

  subroutine set_grids
    use species, only: init_species, nspec, spec
    use species, only: ion_species, electron_species
    use egrid, only: setvgrid, zeroes, zeroes_maxwell, x0, x0_maxwell
    use theta_grid, only: init_theta_grid, ntgrid, nbset, bset, eps
    implicit none

    integer :: tsize,is
    logical :: has_maxwellian_species

    call init_theta_grid
    !call init_species

    allocate (anon(negrid), dele(negrid,nspec), speed(negrid,nspec))
    allocate (w_maxwell(negrid), energy_maxwell(negrid), speed_maxwell(negrid))


    w_maxwell = 0.0
    energy_maxwell = 0.0
    speed_maxwell = 0.0
    has_maxwellian_species = .false. 
    do is = 1,nspec
      !if (spec(is)%type .eq. ion_species .or. &
           !spec(is)%type .eq. electron_species) then
      if (spec(is)%is_maxwellian) then
          has_maxwellian_species = .true. 
          w_maxwell(:) = w(:,is)
          energy_maxwell(:) = energy(:,is)
          speed_maxwell(:) = speed(:,is)
          zeroes_maxwell(:) = zeroes(:,is)
          x0_maxwell = x0(is)
          exit
       end if
    end do
    if (.not. has_maxwellian_species) write(*,*) &
        'Warning; no maxwellian species; collisions will fail'


         
    speed = sqrt(energy)
    
    anon = 1.0

    tsize = 2*nterp-1

    dele(1,:) = energy(1,:)
    dele(2:,:) = energy(2:,:)-energy(:negrid-1,:)

    ng2 = 2*ngauss
!    if (eps > epsilon(0.0)) then
    if (trapped_particles .and. eps > epsilon(0.0)) then
       nlambda = ng2+nbset
       lmax = nlambda-1
    else
       nlambda = ng2
       lmax = nlambda
    end if
    nxi = max(2*nlambda-1, 2*ng2)
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
    if (.not. trapped_particles .or. eps <= epsilon(0.0)) then
       call xigridset
       return
    end if

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

    call xigridset

    eflag = .false.

  end subroutine lgridset

  subroutine xigridset

    use theta_grid, only: ntgrid, bmag

    implicit none

    integer :: ig, je, ixi

    if (.not. allocated(xi)) allocate (xi(-ntgrid:ntgrid,2*nlambda))
    if (.not. allocated(ixi_to_il)) allocate (ixi_to_il(-ntgrid:ntgrid,2*nlambda))
    if (.not. allocated(ixi_to_isgn)) allocate (ixi_to_isgn(-ntgrid:ntgrid,2*nlambda))

    ! define array 'sgn' that returns sign of vpa associated with isgn=1,2
    sgn(1) = 1
    sgn(2) = -1

    ! xi is vpa / v and goes from 1 --> -1
    xi = 0.0
    do ig = -ntgrid, ntgrid
       xi(ig,:jend(ig)) = sqrt(max(1.0 - al(:jend(ig))*spread(bmag(ig),1,jend(ig)),0.0))
       xi(ig,jend(ig)+1:2*jend(ig)-1) = -xi(ig,jend(ig)-1:1:-1)
!       do ixi = 1, nxi
!          write (*,*) 'xi', xi(ig,ixi)
!       end do
    end do

!    xi = max(xi,0.0)
!    xi = sqrt(xi)
    
    ! get mapping from ixi (which runs between 1 and 2*nlambda) and il (runs between 1 and nlambda)
    do ig = -ntgrid, ntgrid
       je = jend(ig)
       ! if no trapped particles
       if (je == 0) then
          do ixi = 1, 2*nlambda
             if (ixi > nlambda) then
                ixi_to_isgn(ig,ixi) = 2
                ixi_to_il(ig,ixi) = 2*nlambda - ixi + 1
             else
                ixi_to_isgn(ig,ixi) = 1
                ixi_to_il(ig,ixi) = ixi
             end if
          end do
       else
          do ixi = 1, 2*nlambda
             if (ixi >= nlambda + je + 1) then
                ixi_to_isgn(ig,ixi) = 2
                ixi_to_il(ig,ixi) = 2*nlambda + je + 1 - ixi
             else if (ixi >= 2*je) then
                ixi_to_isgn(ig,ixi) = 1
                ixi_to_il(ig,ixi) = ixi - je
             else if (ixi >= je) then
                ixi_to_isgn(ig,ixi) = 2
                ixi_to_il(ig,ixi) = 2*je - ixi
             else
                ixi_to_isgn(ig,ixi) = 1
                ixi_to_il(ig,ixi) = ixi
             end if
          end do
       end if
    end do

  end subroutine xigridset

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
    use theta_grid, only: ntgrid
    use kt_grids, only: aky
    use gs2_layouts, only: g_lo, is_idx, ik_idx, ie_idx, il_idx
    use mp, only: nproc, sum_reduce, proc0, sum_allreduce

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:,:,:), intent (out) :: total
    integer, optional, intent(in) :: all
    real :: fac
    integer :: is, il, ie, ik, it, iglo, isgn

    !Initialise to zero
    total = 0.

    !Do integral over local x-y space
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)

       !Pick the weighting factor
       if (aky(ik) == 0.) then
          fac = 1.0
       else
          fac = 0.5
       end if

       !For both signs of vpar do sum
       !May be more efficient to move ign loop above iglo loop (good for total but bad for g memory access)
       do isgn = 1, 2
          total(:, il, ie, isgn, is) = total(:, il, ie, isgn, is) + &
               fac*g(:,isgn,iglo)
       end do
    end do

    !Do we really need this if statement?
    if (nproc > 1) then
       if (present(all)) then
          call sum_allreduce (total)
       else
          call sum_reduce (total, 0)
       end if
    end if

  end subroutine integrate_volume_c

  subroutine integrate_volume_r (g, total, all)
! returns results to PE 0 [or to all processors if 'all' is present in input arg list]
! NOTE: Takes f = f(x, y, z, sigma, lambda, E, species) and returns int f, where the integral
! is over x-y space
    use theta_grid, only: ntgrid
    use kt_grids, only: aky
    use gs2_layouts, only: g_lo, is_idx, ik_idx, ie_idx, il_idx
    use mp, only: nproc,sum_reduce, proc0, sum_allreduce

    implicit none

    real, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (-ntgrid:,:,:,:,:), intent (out) :: total
    integer, optional, intent(in) :: all
    real :: fac
    integer :: is, il, ie, ik, it, iglo, isgn

    !Initialise to zero
    total = 0.

    !Do integral over local x-y space
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)

       !Pick the weighting factor
       if (aky(ik) == 0.) then
          fac = 1.0
       else
          fac = 0.5
       end if

       !For both signs of vpar do sum
       !May be more efficient to move ign loop above iglo loop (good for total but bad for g memory access)
       do isgn = 1, 2
          total(:, il, ie, isgn, is) = total(:, il, ie, isgn, is) + &
               fac*g(:,isgn,iglo)
       end do
    end do

    !Do we really need this if statement?
    if (nproc > 1) then
       if (present(all)) then
          call sum_allreduce (total)
       else
          call sum_reduce (total, 0)
       end if
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
    real, dimension (:,:), allocatable, save :: vpa1d
    real, dimension (:,:,:), allocatable, save :: hermp1d
    real, dimension (:,:,:,:,:), allocatable, save :: vpapts
    real, dimension (:,:,:,:,:,:), allocatable, save :: hermp

    real :: fac
    integer :: is, il, ie, ig, isgn, iv
    integer :: norder

    norder = min(negrid, nlambda)/2

    allocate (favg(-ntgrid:ntgrid,nspec,0:norder-1))

    if (.not. allocated(vpapts)) then
       allocate (vpa1d(negrid*nlambda, nspec))
       allocate (hermp1d(negrid*nlambda,0:norder-1, nspec))
       allocate (vpapts(-ntgrid:ntgrid,nlambda,negrid,2, nspec))
       allocate (hermp(-ntgrid:ntgrid,nlambda,negrid,2,0:norder-1, nspec))
       vpapts = 0.0 ; hermp = 0.0 ; vpa1d = 0.0 ; hermp1d = 0.0

       do ie = 1, negrid
          do il = 1, nlambda
             do ig = -ntgrid, ntgrid
                vpapts(ig,il,ie,1,:) = sqrt(energy(ie,:)*max(0.0, 1.0-al(il)*bmag(ig)))
                vpapts(ig,il,ie,2,:) = -vpapts(ig,il,ie,1,:)
             end do
          end do
       end do

       do iv = 1, negrid*nlambda
          vpa1d(iv,:) = sqrt(energy(negrid,:))*(1. - 2.*(iv-1)/real(negrid*nlambda-1))
       end do

       do is = 1,nspec
         call get_hermite_polynomials (vpa1d(:,is), hermp1d(:,:,is))
         call get_hermite_polynomials (vpapts(:,:,:,:,is), hermp(:,:,:,:,:,is))
       end do
    end if

    favg = 0.
    do is = 1, nspec
       do ie = 1, negrid
          do il = 1, nlambda
             do ig = -ntgrid, ntgrid
                favg(ig,is,:) = favg(ig,is,:) &
                     +w(ie,is)*wl(ig,il)*(hermp(ig,il,ie,1,:,is)*f(ig,il,ie,1,is) &
                     +hermp(ig,il,ie,2,:,is)*f(ig,il,ie,2,is))
             end do
          end do
       end do
    end do

    do is = 1, nspec
       do iv = 1, negrid*nlambda
          do ig = -ntgrid, ntgrid
             vflx(ig,iv,is) = sum(favg(ig,is,:)*hermp1d(iv,:,is))*exp(-vpa1d(iv,is)**2)
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

  subroutine init_map (use_lz_layout, use_e_layout, use_le_layout, test)

    use mp, only: finish_mp, proc0
    use redistribute, only: report_map_property

    implicit none

    logical, intent (in) :: use_lz_layout, use_e_layout, use_le_layout, test

    ! initialize maps from g_lo to lz_lo, e_lo, and/or le_lo

    if (use_lz_layout) then
       ! init_lambda_layout is called in redistribute
       call init_lambda_redistribute
       if (test) then
          if (proc0) print *, '=== Lambda map property ==='
          call report_map_property (lambda_map)
       end if
    end if

    if (use_e_layout) then
       ! init_energy_layout is called in redistribute
       call init_energy_redistribute
       if (test) then
          if (proc0) print *, '=== Energy map property ==='
          call report_map_property (energy_map)
       end if
    end if

    if (use_le_layout) then
       call init_g2le_redistribute
       if (test) call check_g2le
    end if
    
    if (test) then
       if (proc0) print *, 'init_map done'
!!$    call finish_mp
!!$    stop
    end if

  end subroutine init_map

  subroutine init_g2le_redistribute

    use mp, only: nproc
    use species, only: nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: init_le_layouts
    use gs2_layouts, only: g_lo, le_lo
    use gs2_layouts, only: idx_local, proc_id
    use gs2_layouts, only: ik_idx, it_idx, ie_idx, is_idx, il_idx, idx
    use redistribute, only: index_list_type, init_redist, delete_list

    implicit none

    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension (0:nproc-1) :: nn_to, nn_from
    integer, dimension (3) :: from_low, from_high
    integer, dimension (3) :: to_high
    integer :: to_low
    integer :: ig, isign, iglo, il, ile
    integer :: ik, it, ie, is, ixi
    integer :: n, ip, je

    if (leinit) return

    call init_le_layouts (ntgrid, naky, ntheta0, nspec)

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do ig = -ntgrid, ntgrid
          ile = idx (le_lo, ig, ik, it, is)
          if (idx_local(g_lo,iglo)) nn_from(proc_id(le_lo,ile)) = nn_from(proc_id(le_lo,ile)) + 2
          if (idx_local(le_lo,ile)) nn_to(proc_id(g_lo,iglo)) = nn_to(proc_id(g_lo,iglo)) + 2
       end do
    end do

    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first (nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
          allocate (from_list(ip)%third (nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first (nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
          allocate (to_list(ip)%third (nn_to(ip)))
       end if
    end do

    ! get local indices of elements distributed to/from other processors
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
             il = il_idx(g_lo,iglo)
             je = jend(ig)
             if (je == 0) then
                if (isign == 2) then
                   il = 2*g_lo%nlambda+1 - il
                end if
             else
                if (il == je) then
                   if (isign == 1) il = 2*je  ! throw this info away
                else if (il > je) then 
                   if (isign == 1) il = il + je
                   if (isign == 2) il = 2*g_lo%nlambda + 1 - il + je
                else
                   if (isign == 2) il = 2*je - il !+ 1
                end if
             end if
             ile = idx (le_lo, ig, ik, it, is)
             if (idx_local(g_lo,iglo)) then
                ip = proc_id(le_lo,ile)
                n = nn_from(ip) + 1
                nn_from(ip) = n
                from_list(ip)%first(n)  = ig
                from_list(ip)%second(n) = isign
                from_list(ip)%third(n)  = iglo
             end if
             if (idx_local(le_lo,ile)) then
                ip = proc_id(g_lo,iglo)
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to_list(ip)%first(n)  = il
                to_list(ip)%second(n) = ie
                to_list(ip)%third(n)  = ile
             end if
          end do
       end do
    end do

    from_low (1) = -ntgrid
    from_low (2) = 1
    from_low (3) = g_lo%llim_proc

    from_high(1) = ntgrid
    from_high(2) = 2
    from_high(3) = g_lo%ulim_alloc

    to_low = le_lo%llim_proc

    to_high(1) = max(2*nlambda, 2*ng2+1)
    to_high(2) = negrid + 1  ! TT: just followed convention with +1.
    ! TT: It may be good to avoid bank conflict.
    to_high(3) = le_lo%ulim_alloc

    call init_redist (g2le, 'c', to_low, to_high, to_list, from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)

    leinit = .true.

  end subroutine init_g2le_redistribute

  subroutine check_g2le

    use file_utils, only: error_unit
    use mp, only: finish_mp, iproc, nproc, proc0
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, le_lo
    use gs2_layouts, only: ig_idx, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use redistribute, only: gather, scatter, report_map_property

    implicit none

    integer :: iglo, ile, ig, isgn, ik, it, il, ie, is, ierr
    integer :: ip, ixi, je
    complex, dimension (:,:,:), allocatable :: gtmp, letmp

    if (proc0) then
       ierr = error_unit()
    else
       ierr = 6
    end if

    ! report the map property
    if (proc0) print *, '=== g2le map property ==='
    call report_map_property (g2le)

    allocate (gtmp(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (letmp(nlambda*2+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))

    ! gather check
    gtmp = 0.0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             gtmp(ig,isgn,iglo) = rule(ig,isgn,ik,it,il,ie,is)
          end do
       end do
    end do
    call gather (g2le, gtmp, letmp)
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       ig = ig_idx(le_lo,ile)
       je = jend(ig)
       ik = ik_idx(le_lo,ile)
       it = it_idx(le_lo,ile)
       is = is_idx(le_lo,ile)
       do ixi = 1, 2*nlambda
          isgn = ixi_to_isgn(ig,ixi)
          il = ixi_to_il(ig,ixi)
          if (int(real(letmp(ixi,ie,ile))) /= rule(ig,isgn,ik,it,il,ie,is)) &
               write (ierr,'(a,8i6)') 'ERROR: gather by g2le broken!', iproc
       end do

    end do
    if (proc0) write (ierr,'(a)') 'g2le gather check done'

    ! scatter check
    gtmp = 0.0
    call scatter (g2le, letmp, gtmp)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             if (gtmp(ig,isgn,iglo) /= rule(ig,isgn,ik,it,il,ie,is)) &
                  write (ierr,'(a,i6)') 'ERROR: scatter by g2le broken!', iproc
          end do
       end do
    end do
    if (proc0) write (ierr,'(a)') 'g2le scatter check done'

    deallocate (gtmp,letmp)

!    call finish_mp
!    stop

  contains

    function rule (ig, isgn, ik, it, il, ie, is)
      integer, intent (in) :: ig, isgn, ik, it, il, ie, is
      integer :: rule
      rule = ig + isgn + ik + it + il + ie + is  ! make whatever you want
    end function rule

  end subroutine check_g2le

  subroutine init_lambda_redistribute

    use mp, only: nproc
    use species, only: nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: init_lambda_layouts
    use gs2_layouts, only: g_lo, lz_lo
    use gs2_layouts, only: idx_local, proc_id
    use gs2_layouts, only: ik_idx, it_idx, ie_idx, is_idx, il_idx, idx
    use redistribute, only: index_list_type, init_redist, delete_list

    implicit none

    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension(3) :: from_low, from_high
    integer, dimension(2) :: to_high
    integer :: to_low
    integer :: ig, isign, iglo, il, ilz
    integer :: ik, it, ie, is, je
    integer :: n, ip

    if (lzinit) return

    call init_lambda_layouts &
         (ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2)

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
             ilz = idx(lz_lo, ig, ik, it, ie, is)
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
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
             je = jend(ig)
             il = il_idx(g_lo,iglo)
             if (je == 0) then
                if (isign == 2) then
                   il = 2*g_lo%nlambda+1 - il
                end if
             else
                if (il == je) then
                   if (isign == 1) il = 2*je  ! throw this info away
                else if (il > je) then 
                   if (isign == 1) il = il + je
                   if (isign == 2) il = 2*g_lo%nlambda + 1 - il + je
                else
                   if (isign == 2) il = 2*je - il !+ 1
                end if
             end if
             ilz = idx(lz_lo, ig, ik, it, ie, is)
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
    call init_redist (lambda_map, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)
    call delete_list (to_list)
    call delete_list (from_list)

    lzinit = .true.

  end subroutine init_lambda_redistribute

  subroutine init_energy_redistribute

    use mp, only: nproc
    use species, only: nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: init_energy_layouts
    use gs2_layouts, only: g_lo, e_lo, ie_idx, ik_idx, it_idx, il_idx, is_idx
    use gs2_layouts, only: idx_local, proc_id, idx
    use redistribute, only: index_list_type, init_redist, delete_list

    implicit none

    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension(3) :: from_low, from_high
    integer, dimension(2) :: to_high
    integer :: to_low
    integer :: ig, isign, iglo, ik, it, il, ie, is, ielo
    integer :: n, ip

    if (einit) return
    
!DD, March 2009: Nullify pointers on initialisation, so do not pass association
!                test when not allocated. 
!  Problem arose on Pascali compiler (York) with  to_list(ip)%third and fourth
    do ip=0, (nproc-1)
       nullify(to_list(ip)%first,from_list(ip)%first,to_list(ip)%second,from_list(ip)%second,to_list(ip)%third,from_list(ip)%third,to_list(ip)%fourth,from_list(ip)%fourth)
    end do
!<DD
	
    call init_energy_layouts &
         (ntgrid, naky, ntheta0, nlambda, nspec)

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
             ielo = idx(e_lo,ig,isign,ik,it,il,is)
             if (idx_local(g_lo,iglo)) &
                  nn_from(proc_id(e_lo,ielo)) = nn_from(proc_id(e_lo,ielo)) + 1
             if (idx_local(e_lo,ielo)) &
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
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
             ielo = idx(e_lo,ig,isign,ik,it,il,is)

             if (idx_local(g_lo,iglo)) then
                ip = proc_id(e_lo,ielo)
                n = nn_from(ip) + 1
                nn_from(ip) = n
                from_list(ip)%first(n) = ig
                from_list(ip)%second(n) = isign
                from_list(ip)%third(n) = iglo
             end if
             if (idx_local(e_lo,ielo)) then
                ip = proc_id(g_lo,iglo)
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to_list(ip)%first(n) = ie
                to_list(ip)%second(n) = ielo
             end if
          end do
       end do
    end do

    from_low (1) = -ntgrid
    from_low (2) = 1
    from_low (3) = g_lo%llim_proc

    to_low = e_lo%llim_proc

    to_high(1) = negrid+1
    to_high(2) = e_lo%ulim_alloc

    from_high(1) = ntgrid
    from_high(2) = 2
    from_high(3) = g_lo%ulim_alloc
    
    call init_redist (energy_map, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)

    einit = .true.

  end subroutine init_energy_redistribute


  ! subroutine used for testing
  ! takes as input an array using g_lo and
  ! writes it to a .distmp output file
  subroutine write_mpdist (dist, extension, last)

    use mp, only: proc0, send, receive
    use file_utils, only: open_output_file, close_output_file
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use gs2_layouts, only: ie_idx, il_idx, idx_local, proc_id
    use theta_grid, only: ntgrid, bmag, theta

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: dist
!    real, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: dist
    character (*), intent (in) :: extension
    logical, intent (in), optional :: last

    integer :: iglo, ik, it, is, ie, il, ig
    integer, save :: unit
    logical :: done = .false.
    real :: vpa1
    complex, dimension (2) :: gtmp
!    real, dimension (2) :: gtmp

     if (.not. done) then
!        if (proc0) call open_output_file (unit, ".distmp")
        if (proc0) call open_output_file (unit, trim(extension))
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
                    vpa1 = sqrt(energy_maxwell(ie)*max(0.0,1.0-al(il)*bmag(ig)))
                    write (unit,'(7e14.5)') theta(ig), vpa1, energy_maxwell(ie), &
                         real(gtmp(1)), aimag(gtmp(1)), real(gtmp(2)), aimag(gtmp(2))
!                    write (unit,'(5e14.5)') theta(ig), vpa1, energy(ie), gtmp(1), gtmp(2)
                 end if
              end if
           end do
           if (proc0) then
              write (unit,*)
              write (unit,*)
           end if
        end do
        if (proc0) call close_output_file (unit)
        if (present(last)) done = .true.
     end if

   end subroutine write_mpdist

  ! subroutine used for testing
  ! takes as input an array using le_lo and
  ! writes it to a .distmp output file
  subroutine write_mpdist_le (dist, extension, last)

    use mp, only: proc0, send, receive
    use file_utils, only: open_output_file, close_output_file
    use gs2_layouts, only: le_lo, ig_idx, ik_idx, it_idx, is_idx
    use gs2_layouts, only: idx_local, proc_id
    use theta_grid, only: ntgrid, bmag, theta

    implicit none

    complex, dimension (:,:,le_lo%llim_proc:), intent (in) :: dist
!    real, dimension (:,:,le_lo%llim_proc:), intent (in) :: dist
    character (*), intent (in) :: extension
    logical, intent (in), optional :: last

    integer :: ile, ixi, nxi, ik, it, is, ie, il, ig, isgn
    integer, save :: unit
    logical :: done = .false.
    real :: vpa1
    complex :: gtmp
!    real :: gtmp

     if (.not. done) then
        if (proc0) call open_output_file (unit, trim(extension))
        nxi = max(2*nlambda-1, 2*ng2)
        do ile=le_lo%llim_world, le_lo%ulim_world
           ig = ig_idx(le_lo, ile)
           ik = ik_idx(le_lo, ile)
           it = it_idx(le_lo, ile)
           is = is_idx(le_lo, ile)
           do ie = 1, negrid
              do ixi = 1, nxi
                 il = ixi_to_il(ig,ixi)
                 isgn = ixi_to_isgn(ig,ixi)
                 if (idx_local (le_lo, ig, ik, it, is)) then
                    if (proc0) then
                       gtmp = dist(ixi,ie,ile)
                    else
                       call send (dist(ixi,ie,ile), 0)
                    end if
                 else if (proc0) then
                    call receive (gtmp, proc_id(le_lo, ile))
                 end if
                 if (proc0) then
                    if (.not. forbid(ig,il)) then
                       vpa1 = sqrt(energy_maxwell(ie)*max(0.0,1.0-al(il)*bmag(ig)))
                       write (unit,'(5e14.5)') theta(ig), sgn(isgn)*vpa1, energy_maxwell(ie), &
                            real(gtmp), aimag(gtmp)
!                    write (unit,'(5e14.5)') theta(ig), vpa1, energy(ie), gtmp(1), gtmp(2)
                    end if
                 end if
              end do
           end do
        end do
        if (proc0) call close_output_file (unit)
        if (present(last)) done = .true.
     end if

   end subroutine write_mpdist_le

  subroutine finish_le_grids

    use egrid, only: zeroes, x0, zeroes_maxwell

    implicit none

    if (allocated(zeroes)) deallocate (zeroes,x0,zeroes_maxwell)
    if (allocated(energy)) deallocate (energy, energy_maxwell, dele, al, wl, jend, forbid, xx, lgrnge, xloc, speed_maxwell, speed)
    if (allocated(integration_work)) deallocate (integration_work)
    if (allocated(wlerr)) deallocate (wlerr)
    if (allocated(werr)) deallocate (werr)
    if (allocated(wlterr)) deallocate (wlterr)
    if (allocated(w)) deallocate (w, w_maxwell)
    if (allocated(anon)) deallocate (anon)
    if (allocated(delal)) deallocate (delal)
    if (allocated(lpl)) deallocate (lpl, lpe)
    if (allocated(lpt)) deallocate (lpt)
    if (allocated(wlmod)) deallocate (wlmod)
    if (allocated(wmod)) deallocate (wmod)
    if (allocated(wtmod)) deallocate (wtmod)
    if (allocated(xi)) deallocate (xi)
    if (allocated(ixi_to_il)) deallocate (ixi_to_il)
    if (allocated(ixi_to_isgn)) deallocate (ixi_to_isgn)

    accel_x = .false. ; accel_v = .false.
    test = .false. ; trapped_particles = .true.
    new_trap_int = .false. 

    intinit = .false. ; slintinit = .false. ; lintinit = .false. ; eintinit = .false.

    initialized = .false.

  end subroutine finish_le_grids

end module le_grids


 