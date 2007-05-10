module fields_implicit
  implicit none

  public :: init_fields_implicit
  public :: advance_implicit
  public :: init_phi_implicit

  private

  complex, dimension (:,:), allocatable :: aminv
  ! (nidx, -*- f_lo -*-)

  integer, save :: nidx, nfield

contains

  subroutine init_fields_implicit
    use theta_grid, only: init_theta_grid, ntgrid
    use kt_grids, only: init_kt_grids, aky, naky, ntheta0
    use mp, only: proc0
    implicit none
    logical, save :: initialized = .false.
    integer :: ik

    if (initialized) return
    initialized = .true.

    call init_theta_grid
    call init_kt_grids

    if (proc0) call read_parameters
    call broadcast_parameters
    call init_response_matrix
  end subroutine init_fields_implicit

  subroutine read_parameters
    use file_utils, only: input_unit
    implicit none
  end subroutine read_parameters

  subroutine broadcast_parameters
    use mp, only: broadcast
    implicit none
  end subroutine broadcast_parameters

  subroutine init_phi_implicit
    use fields_arrays, only: phi, apar, aperp, phinew, aparnew, aperpnew
    implicit none

    call getfield
    phi = phinew; apar = aparnew; aperp = aperpnew
  end subroutine init_phi_implicit

  subroutine init_response_matrix
    use fields_arrays, only: phi, apar, aperp, phinew, aparnew, aperpnew
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use dist_fn_arrays, only: g
    use run_parameters, only: fphi, fapar, faperp
    use gs2_layouts, only: init_fields_layouts, f_lo
    use mp, only: proc0
    use prof, only: prof_entering, prof_leaving
    implicit none
    integer :: ig, ifield
    complex, dimension (:,:), allocatable :: am
    integer :: iflo

    call prof_entering ("init_response_matrix", "fields_implicit")

    nfield = 0
    if (fphi > epsilon(0.0)) nfield = nfield + 1
    if (fapar > epsilon(0.0)) nfield = nfield + 1
    if (faperp > epsilon(0.0)) nfield = nfield + 1
    nidx = (2*ntgrid+1)*nfield

    call init_fields_layouts (nidx, naky, ntheta0)
    allocate (am(nidx,f_lo%llim_proc:f_lo%ulim_proc))

    g = 0.0

    phi = 0.0
    apar = 0.0
    aperp = 0.0
    phinew = 0.0
    aparnew = 0.0
    aperpnew = 0.0

    do ig = -ntgrid, ntgrid
       ifield = 0
       if (fphi > epsilon(0.0)) then
          ifield = ifield + 1
          phinew(ig,:,:) = 1.0
          call init_response_row (ig, ifield, am)
          phinew(ig,:,:) = 0.0
       end if

       if (fapar > epsilon(0.0)) then
          ifield = ifield + 1
          aparnew(ig,:,:) = 1.0
          call init_response_row (ig, ifield, am)
          aparnew(ig,:,:) = 0.0
       end if

       if (faperp > epsilon(0.0)) then
          ifield = ifield + 1
          aperpnew(ig,:,:) = 1.0
          call init_response_row (ig, ifield, am)
          aperpnew(ig,:,:) = 0.0
       end if
    end do

    call init_inverse_matrix (am)
    deallocate (am)

    call prof_leaving ("init_response_matrix", "fields_implicit")

  end subroutine init_response_matrix

  subroutine init_response_row (ig, ifield, am)
    use fields_arrays, only: phi, apar, aperp, phinew, aparnew, aperpnew
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use dist_fn, only: getfieldeq, timeadv
    use run_parameters, only: fphi, fapar, faperp
    use gs2_layouts, only: f_lo, idx, idx_local
    use prof, only: prof_entering, prof_leaving
    implicit none
    integer, intent (in) :: ig, ifield
    complex, dimension (:,f_lo%llim_proc:), intent (in out) :: am
    complex, dimension (-ntgrid:ntgrid,naky,ntheta0) :: &
         fieldeq, fieldeqa, fieldeqp
    integer :: irow, istart, ik, it, iflo

    call prof_entering ("init_response_row", "fields_implicit")

    call timeadv (phi, apar, aperp, phinew, aparnew, aperpnew, 0)
    call getfieldeq (phinew, aparnew, aperpnew, fieldeq, fieldeqa, fieldeqp)

    irow = (ig+ntgrid)*nfield + ifield
    do it = 1, ntheta0
       do ik = 1, naky
          iflo = idx(f_lo, irow, ik, it)
          if (idx_local(f_lo, iflo)) then
             istart = 0
             if (fphi > epsilon(0.0)) then
                istart = istart + 1
                am(istart:nidx:nfield,iflo) = fieldeq(:,ik,it)
             end if
             if (fapar > epsilon(0.0)) then
                istart = istart + 1
                am(istart:nidx:nfield,iflo) = fieldeqa(:,ik,it)
             end if
             if (faperp > epsilon(0.0)) then
                istart = istart + 1
                am(istart:nidx:nfield,iflo) = fieldeqp(:,ik,it)
             end if
          end if
       end do
    end do

    call prof_leaving ("init_response_row", "fields_implicit")
  end subroutine init_response_row

  subroutine init_inverse_matrix (am)
    use file_utils, only: error_unit
    use kt_grids, only: naky, ntheta0
    use mp, only: broadcast
    use gs2_layouts, only: f_lo, idx, idx_local, proc_id
    use gs2_layouts, only: if_idx, ik_idx, it_idx
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (:,f_lo%llim_proc:), intent (in out) :: am
    complex, dimension (nidx,naky,ntheta0) :: lhscol, rhsrow
    complex :: fac
    integer :: i, ik, it
    integer :: irow, ilo, jlo

    call prof_entering ("init_inverse_matrix", "fields_implicit")

    allocate (aminv(nidx,f_lo%llim_proc:f_lo%ulim_proc))

    aminv = 0.0
    do ilo = f_lo%llim_proc, f_lo%ulim_proc
       aminv(if_idx(f_lo,ilo),ilo) = 1.0
    end do

    ! Gauss-Jordan elimination
    do i = 1, nidx
       do it = 1, ntheta0
          do ik = 1, naky
             ilo = idx(f_lo,i,ik,it)
             if (idx_local(f_lo,ilo)) then
                lhscol(:,ik,it) = am(:,ilo)
                rhsrow(:,ik,it) = aminv(:,ilo)
             end if
             call broadcast (lhscol(:,ik,it), proc_id(f_lo,ilo))
             call broadcast (rhsrow(:,ik,it), proc_id(f_lo,ilo))
          end do
       end do

       do jlo = f_lo%llim_proc, f_lo%ulim_proc
          irow = if_idx(f_lo,jlo)
          ik = ik_idx(f_lo,jlo)
          it = it_idx(f_lo,jlo)

          fac = am(i,jlo)/lhscol(i,ik,it)
          am(i,jlo) = am(i,jlo)/lhscol(i,ik,it)
          am(:i-1,jlo) = am(:i-1,jlo) - lhscol(:i-1,ik,it)*fac
          am(i+1:,jlo) = am(i+1:,jlo) - lhscol(i+1:,ik,it)*fac

          if (irow == i) then
             aminv(:,jlo) = aminv(:,jlo)/lhscol(i,ik,it)
          else
             aminv(:,jlo) = aminv(:,jlo) &
                  - rhsrow(:,ik,it)*lhscol(irow,ik,it)/lhscol(i,ik,it)
          end if
       end do
    end do

    call prof_leaving ("init_inverse_matrix", "fields_implicit")
  end subroutine init_inverse_matrix

  subroutine get_field_vector (fl)
    use fields_arrays, only: phi, apar, aperp, phinew, aparnew, aperpnew
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use dist_fn, only: getfieldeq
    use run_parameters, only: fphi, fapar, faperp
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (:,:,:), intent (out) :: fl
    complex, dimension (-ntgrid:ntgrid,naky,ntheta0) :: &
         fieldeq, fieldeqa, fieldeqp
    integer :: istart

    call prof_entering ("get_field_vector", "fields_implicit")

    call getfieldeq (phinew, aparnew, aperpnew, fieldeq, fieldeqa, fieldeqp)

    istart = 0

    if (fphi > epsilon(0.0)) then
       istart = istart + 1
       fl(istart:nidx:nfield,:,:) = fieldeq
    end if

    if (fapar > epsilon(0.0)) then
       istart = istart + 1
       fl(istart:nidx:nfield,:,:) = fieldeqa
    end if

    if (faperp > epsilon(0.0)) then
       istart = istart + 1
       fl(istart:nidx:nfield,:,:) = fieldeqp
    end if

    call prof_leaving ("get_field_vector", "fields_implicit")
  end subroutine get_field_vector

  subroutine get_field_solution (u)
    use fields_arrays, only: phi, apar, aperp, phinew, aparnew, aperpnew
    use theta_grid, only: ntgrid
    use run_parameters, only: fphi, fapar, faperp
    use gs2_layouts, only: f_lo, if_idx, ik_idx, it_idx
    use gs2_layouts, only: idx_local, proc_id
    use mp, only: broadcast
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (f_lo%llim_proc:), intent (in) :: u
    integer :: istart, iflo, if, ik, it, ig

    call prof_entering ("get_field_solution", "fields_implicit")

    istart = 0

    if (fphi > epsilon(0.0)) then
       istart = istart + 1
       do iflo = f_lo%llim_world, f_lo%ulim_world
          if = if_idx(f_lo,iflo)
          if (mod((if - istart),nfield) /= 0) cycle
          ik = ik_idx(f_lo,iflo)
          it = it_idx(f_lo,iflo)
          ig = -ntgrid + (if-istart)/nfield
          if (idx_local(f_lo,iflo)) phinew(ig,ik,it) = -u(iflo)
          call broadcast (phinew(ig,ik,it), proc_id(f_lo,iflo))
       end do
    end if

    if (fapar > epsilon(0.0)) then
       istart = istart + 1
       do iflo = f_lo%llim_world, f_lo%ulim_world
          if = if_idx(f_lo,iflo)
          if (mod((if - istart),nfield) /= 0) cycle
          ik = ik_idx(f_lo,iflo)
          it = it_idx(f_lo,iflo)
          ig = -ntgrid + (if-istart)/nfield
          if (idx_local(f_lo,iflo)) aparnew(ig,ik,it) = -u(iflo)
          call broadcast (aparnew(ig,ik,it), proc_id(f_lo,iflo))
       end do
    end if

    if (faperp > epsilon(0.0)) then
       istart = istart + 1
       do iflo = f_lo%llim_world, f_lo%ulim_world
          if = if_idx(f_lo,iflo)
          if (mod((if - istart),nfield) /= 0) cycle
          ik = ik_idx(f_lo,iflo)
          it = it_idx(f_lo,iflo)
          ig = -ntgrid + (if-istart)/nfield
          if (idx_local(f_lo,iflo)) aperpnew(ig,ik,it) = -u(iflo)
          call broadcast (aperpnew(ig,ik,it), proc_id(f_lo,iflo))
       end do
    end if

    call prof_leaving ("get_field_solution", "fields_implicit")
  end subroutine get_field_solution

  subroutine getfield
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: f_lo, ik_idx, it_idx
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (nidx, naky, ntheta0) :: fl
    complex, dimension (f_lo%llim_proc:f_lo%ulim_proc) :: u
    integer :: iflo

    call prof_entering ("getfield", "fields_implicit")

    ! am*u = fl, Poisson's and Ampere's law, u is phi, apar, aperp interleaved
    ! u = aminv*fl

    call get_field_vector (fl)

    do iflo = f_lo%llim_proc, f_lo%ulim_proc
       u(iflo) = sum(aminv(:,iflo)*fl(:,ik_idx(f_lo,iflo),it_idx(f_lo,iflo)))
    end do

    call get_field_solution (u)

    call prof_leaving ("getfield", "fields_implicit")
  end subroutine getfield

  subroutine advance_implicit (istep)
    use fields_arrays, only: phi, apar, aperp, phinew, aparnew, aperpnew
    use dist_fn, only: timeadv
    use dist_fn_arrays, only: g, gnew
    implicit none
    integer, intent (in) :: istep

    phi = phinew; apar = aparnew; aperp = aperpnew
    g = gnew
    call timeadv (phi, apar, aperp, phinew, aparnew, aperpnew, istep)
    call getfield
    phinew   = phinew   + phi
    aparnew  = aparnew  + apar
    aperpnew = aperpnew + aperp
    call timeadv (phi, apar, aperp, phinew, aparnew, aperpnew, istep)
  end subroutine advance_implicit

end module fields_implicit
