module fields_implicit
  use fields_arrays, only: nidx
  implicit none

  public :: init_fields_implicit
  public :: advance_implicit
  public :: init_phi_implicit
  public :: nidx
  public :: reset_init

  private

  integer, save :: nfield
  logical :: initialized = .false.
  logical :: linked = .false.

contains

  subroutine init_fields_implicit
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids
    implicit none

    if (initialized) return
    initialized = .true.

    call init_theta_grid
    call init_kt_grids
    call read_parameters
    call init_response_matrix

  end subroutine init_fields_implicit

  subroutine read_parameters
  end subroutine read_parameters

  subroutine init_phi_implicit
    use fields_arrays, only: phi, apar, aperp, phinew, aparnew, aperpnew
    implicit none

    call init_fields_implicit
    call getfield (phinew, aparnew, aperpnew)
    phi = phinew; apar = aparnew; aperp = aperpnew
  end subroutine init_phi_implicit

  subroutine get_field_vector (fl, phi, apar, aperp)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use dist_fn, only: getfieldeq
    use run_parameters, only: fphi, fapar, faperp
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, aperp
    complex, dimension (:,:,:), intent (out) :: fl
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: fieldeq, fieldeqa, fieldeqp
    integer :: istart, ifin

    call prof_entering ("get_field_vector", "fields_implicit")

    call getfieldeq (phi, apar, aperp, fieldeq, fieldeqa, fieldeqp)

    ifin = 0

    if (fphi > epsilon(0.0)) then
       istart = ifin + 1
       ifin = (istart-1) + 2*ntgrid+1
       fl(istart:ifin,:,:) = fieldeq
    end if

    if (fapar > epsilon(0.0)) then
       istart = ifin + 1
       ifin = (istart-1) + 2*ntgrid+1
       fl(istart:ifin,:,:) = fieldeqa
    end if

    if (faperp > epsilon(0.0)) then
       istart = ifin + 1
       ifin = (istart-1) + 2*ntgrid+1
       fl(istart:ifin,:,:) = fieldeqp
    end if

    call prof_leaving ("get_field_vector", "fields_implicit")
  end subroutine get_field_vector

  subroutine get_field_solution (u)
    use fields_arrays, only: phinew, aparnew, aperpnew
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: fphi, fapar, faperp
    use gs2_layouts, only: jf_lo, ij_idx
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (0:), intent (in) :: u
    integer :: ik, it, ifield, ll, lr

    call prof_entering ("get_field_solution", "fields_implicit")

    ifield = 0

    if (fphi > epsilon(0.0)) then
       ifield = ifield + 1
       do ik = 1, naky
          do it = 1, ntheta0
             ll = ij_idx (jf_lo, -ntgrid, ifield, ik, it)
             lr = ll + 2*ntgrid
             phinew(:,it,ik) = u(ll:lr)
          end do
       end do
    endif

    if (fapar > epsilon(0.0)) then
       ifield = ifield + 1
       do ik = 1, naky
          do it = 1, ntheta0
             ll = ij_idx (jf_lo, -ntgrid, ifield, ik, it)
             lr = ll + 2*ntgrid
             aparnew(:,it,ik) = u(ll:lr)
          end do
       end do
    endif

    if (faperp > epsilon(0.0)) then
       ifield = ifield + 1
       do ik = 1, naky
          do it = 1, ntheta0
             ll = ij_idx (jf_lo, -ntgrid, ifield, ik, it)
             lr = ll + 2*ntgrid
             aperpnew(:,it,ik) = u(ll:lr)
          end do
       end do
    endif

    call prof_leaving ("get_field_solution", "fields_implicit")
  end subroutine get_field_solution

  subroutine getfield (phi, apar, aperp)
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: f_lo, jf_lo, ij, mj, dj
    use prof, only: prof_entering, prof_leaving
    use fields_arrays, only: aminv
    use theta_grid, only: ntgrid
    use dist_fn, only: N_class
    use mp, only: sum_allreduce
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, aperp
    complex, dimension (nidx, ntheta0, naky) :: fl
    complex, dimension (0:nidx*ntheta0*naky-1) :: u
    integer :: jflo, ik, it, nl, nr, i, m, n, dc

    call prof_entering ("getfield", "fields_implicit")

    ! am*u = fl, Poisson's and Ampere's law, u is phi, apar, aperp 
    ! u = aminv*fl

    call get_field_vector (fl, phi, apar, aperp)
    
    u = 0.
    do jflo = jf_lo%llim_proc, jf_lo%ulim_proc

       i = ij(jflo)
       m = mj(jflo)
       ik = f_lo(i)%ik(m,1)  ! For fixed i and m, ik does not change as n varies 
       dc = dj(i,jflo)
       
       do n = 1, N_class(i)
                    
          it = f_lo(i)%it(m,n)
          
          nl = 1 + nidx*(n-1)
          nr = nl + nidx - 1

          u(jflo)=u(jflo)-sum(aminv(i)%dcell(dc)%supercell(nl:nr)*fl(:, it, ik)) 

       end do
    end do

    call sum_allreduce (u)

    call get_field_solution (u)

    call prof_leaving ("getfield", "fields_implicit")
  end subroutine getfield

  subroutine advance_implicit (istep, dt_cfl)
    use fields_arrays, only: phi, apar, aperp, phinew, aparnew, aperpnew
    use fields_arrays, only: apar_ext, phi_ext
    use dist_fn, only: timeadv, get_apar_ext, get_phi_ext
    use dist_fn_arrays, only: g, gnew
    use nonlinear_terms, only: algorithm !, nonlin
    implicit none
    integer, intent (in) :: istep
    real :: dt_cfl

    if (algorithm == 1) then

       g = gnew
       phi = phinew 
       apar = aparnew
       aperp = aperpnew       
       
       call timeadv (phi, apar, aperp, phinew, aparnew, aperpnew, istep, dt_cfl)
       call getfield (phinew, aparnew, aperpnew)

       call get_apar_ext (apar_ext)
       call get_phi_ext (phi_ext)

       phinew   = phinew   + phi + phi_ext
       aparnew  = aparnew  + apar + apar_ext
       aperpnew = aperpnew + aperp
                 
       call timeadv (phi, apar, aperp, phinew, aparnew, aperpnew, istep, dt_cfl)
    end if

  end subroutine advance_implicit

  subroutine reset_init

    use fields_arrays, only: aminv
    integer :: i, j
    initialized = .false.

    if (.not. allocated (aminv)) return
    do i = 1, size(aminv)
       if (.not. associated (aminv(i)%dcell)) cycle
       do j = 1, size(aminv(i)%dcell)
          if (associated (aminv(i)%dcell(j)%supercell)) &
               deallocate(aminv(i)%dcell(j)%supercell)
       end do
       if (associated (aminv(i)%dcell)) deallocate (aminv(i)%dcell)
    end do
    deallocate (aminv)

  end subroutine reset_init

  subroutine init_response_matrix
    use mp, only: barrier
    use fields_arrays, only: phi, apar, aperp, phinew, aparnew, aperpnew
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use dist_fn_arrays, only: g
    use dist_fn, only: M_class, N_class, i_class
    use run_parameters, only: fphi, fapar, faperp
    use gs2_layouts, only: init_fields_layouts, f_lo
    use gs2_layouts, only: init_jfields_layouts
    use prof, only: prof_entering, prof_leaving
    implicit none
    integer :: ig, ifield, it, ik, i, m, n
    complex, dimension(:,:), allocatable :: am
    logical :: endpoint

    call prof_entering ("init_response_matrix", "fields_implicit")

    nfield = 0
    if (fphi > epsilon(0.0)) nfield = nfield + 1
    if (fapar > epsilon(0.0)) nfield = nfield + 1
    if (faperp > epsilon(0.0)) nfield = nfield + 1
    nidx = (2*ntgrid+1)*nfield

    call init_fields_layouts (nfield, nidx, naky, ntheta0, M_class, N_class, i_class)
    call init_jfields_layouts (nfield, nidx, naky, ntheta0, i_class)
    call finish_fields_layouts

! keep storage cost down by doing one class at a time
    do i = i_class, 1, -1

       call barrier
       allocate (am(nidx*N_class(i), f_lo(i)%llim_proc:f_lo(i)%ulim_alloc))

       am = 0.0
       g = 0.0
       
       phi = 0.0
       apar = 0.0
       aperp = 0.0
       phinew = 0.0
       aparnew = 0.0
       aperpnew = 0.0

       do n = 1, N_class(i)
          do ig = -ntgrid, ntgrid
             endpoint = n > 1
             endpoint = ig == -ntgrid .and. endpoint
             ifield = 0
             if (fphi > epsilon(0.0)) then
                ifield = ifield + 1
                if (endpoint) then
                   do m = 1, M_class(i)
                      ik = f_lo(i)%ik(m,n-1)
                      it = f_lo(i)%it(m,n-1)
                      phinew(ntgrid,it,ik) = 1.0
                   end do
                endif
                do m = 1, M_class(i)
                   ik = f_lo(i)%ik(m,n)
                   it = f_lo(i)%it(m,n)
                   phinew(ig,it,ik) = 1.0
                end do
                call init_response_row (ig, ifield, am, i, n)
                phinew = 0.0
             end if
             
             if (fapar > epsilon(0.0)) then
                ifield = ifield + 1
                if (endpoint) then
                   do m = 1, M_class(i)
                      ik = f_lo(i)%ik(m,n-1)
                      it = f_lo(i)%it(m,n-1)
                      aparnew(ntgrid,it,ik) = 1.0
                   end do
                endif
                do m = 1, M_class(i)
                   ik = f_lo(i)%ik(m,n)
                   it = f_lo(i)%it(m,n)
                   aparnew(ig,it,ik) = 1.0
                end do
                call init_response_row (ig, ifield, am, i, n)
                aparnew = 0.0
             end if
             
             if (faperp > epsilon(0.0)) then
                ifield = ifield + 1
                if (endpoint) then
                   do m = 1, M_class(i)
                      ik = f_lo(i)%ik(m,n-1)
                      it = f_lo(i)%it(m,n-1)
                      aperpnew(ntgrid,it,ik) = 1.0
                   end do
                endif
                do m = 1, M_class(i)
                   ik = f_lo(i)%ik(m,n)
                   it = f_lo(i)%it(m,n)
                   aperpnew(ig,it,ik) = 1.0
                end do
                call init_response_row (ig, ifield, am, i, n)
                aperpnew = 0.0
             end if
          end do
       end do

       call init_inverse_matrix (am, i)
       deallocate (am)
    end do

    call prof_leaving ("init_response_matrix", "fields_implicit")

  end subroutine init_response_matrix

  subroutine init_response_row (ig, ifield, am, ic, n)
    use fields_arrays, only: phi, apar, aperp, phinew, aparnew, aperpnew
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use dist_fn, only: getfieldeq, timeadv, M_class, N_class
    use run_parameters, only: fphi, fapar, faperp
    use gs2_layouts, only: f_lo, idx, idx_local
    use prof, only: prof_entering, prof_leaving
    implicit none
    integer, intent (in) :: ig, ifield, ic, n
    complex, dimension(:,f_lo(ic)%llim_proc:), intent (in out) :: am
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: fieldeq, fieldeqa, fieldeqp
    integer :: irow, istart, iflo, ik, it, ifin, m, nn
    real :: dt_cfl

    call prof_entering ("init_response_row", "fields_implicit")

    call timeadv (phi, apar, aperp, phinew, aparnew, aperpnew, 0, dt_cfl)
    call getfieldeq (phinew, aparnew, aperpnew, fieldeq, fieldeqa, fieldeqp)

    do nn = 1, N_class(ic)
       do m = 1, M_class(ic)

          it = f_lo(ic)%it(m,nn)
          ik = f_lo(ic)%ik(m,nn)
       
          irow = ifield + nfield*((ig+ntgrid) + (2*ntgrid+1)*(n-1))
          
          iflo = idx (f_lo(ic), irow, m)
          
          if (idx_local(f_lo(ic), iflo)) then
             
             istart = 0 + nidx*(nn-1)
             
             if (fphi > epsilon(0.0)) then
                ifin = istart + nidx
                istart = istart + 1
                am(istart:ifin:nfield,iflo) = fieldeq(:,it,ik) 
             end if
             
             if (fapar > epsilon(0.0)) then
                ifin = istart + nidx
                istart = istart + 1
                am(istart:ifin:nfield,iflo) = fieldeqa(:,it,ik)
             end if
             
             if (faperp > epsilon(0.0)) then
                ifin = istart + nidx
                istart = istart + 1
                am(istart:ifin:nfield,iflo) = fieldeqp(:,it,ik)
             end if
             
          end if
                    
       end do
    end do

    call prof_leaving ("init_response_row", "fields_implicit")
  end subroutine init_response_row

  subroutine init_inverse_matrix (am, ic)
    use mp, only: iproc
    use file_utils, only: error_unit
    use kt_grids, only: aky, akx
    use theta_grid, only: ntgrid
    use mp, only: broadcast, send, receive, barrier
    use gs2_layouts, only: f_lo, idx, idx_local, proc_id, jf_lo
    use gs2_layouts, only: if_idx, im_idx, in_idx
    use gs2_layouts, only: ig_idx, ifield_idx, ij_idx, mj, dj
    use prof, only: prof_entering, prof_leaving
    use fields_arrays, only: aminv
    use dist_fn, only: i_class, M_class, N_class
    implicit none
    integer, intent (in) :: ic
    complex, dimension(:,f_lo(ic)%llim_proc:), intent (in out) :: am
    complex, dimension(:,:), allocatable :: a_inv
    complex, dimension (nidx*N_class(ic),M_class(ic)) :: lhscol, rhsrow
    complex, dimension (:), allocatable :: am_tmp
    complex :: fac
    integer :: i, j, k, ik, it, m, n, nn, if, ig, jsc, jf, jg, jc
    integer :: irow, ilo, jlo, dc, iflo, ierr
    logical :: iskip, jskip

    call prof_entering ("init_inverse_matrix", "fields_implicit")
    
    call barrier
   
    j = nidx*N_class(ic)
    allocate (a_inv(j,f_lo(ic)%llim_proc:f_lo(ic)%ulim_alloc))
    a_inv = 0.0
    
    do ilo = f_lo(ic)%llim_proc, f_lo(ic)%ulim_proc
       a_inv(if_idx(f_lo(ic),ilo),ilo) = 1.0
    end do

    ! Gauss-Jordan elimination, leaving out internal points at multiples of ntgrid 
    ! for each supercell
    do i = 1, nidx*N_class(ic)
       iskip = N_class(ic) > 1
       iskip = i <= nidx*N_class(ic) - nfield .and. iskip
       iskip = mod((i+nfield-1)/nfield, 2*ntgrid+1) == 0 .and. iskip
       iskip = i > nfield .and. iskip
       if (iskip) cycle
 
       do m = 1, M_class(ic)
          ilo = idx(f_lo(ic),i,m)
          if (idx_local(f_lo(ic),ilo)) then
             lhscol(:,m) = am(:,ilo)
             rhsrow(:,m) = a_inv(:,ilo)
          end if
          call broadcast (lhscol(:,m), proc_id(f_lo(ic),ilo))
          call broadcast (rhsrow(:,m), proc_id(f_lo(ic),ilo))
       end do

       do jlo = f_lo(ic)%llim_proc, f_lo(ic)%ulim_proc

          jskip = N_class(ic) > 1
          jskip = ig_idx(f_lo(ic), jlo) == ntgrid .and. jskip

          n = in_idx(f_lo(ic),jlo)

          jskip = n < N_class(ic) .and. jskip
          if (jskip) cycle          

          m = im_idx(f_lo(ic),jlo)

          ik = f_lo(ic)%ik(m,n)
          it = f_lo(ic)%it(m,n)
          
          irow = if_idx(f_lo(ic),jlo)

          if (aky(ik) /= 0.0 .or. akx(it) /= 0.0) then
             fac = am(i,jlo)/lhscol(i,m)
             am(i,jlo) = fac
             am(:i-1,jlo) = am(:i-1,jlo) - lhscol(:i-1,m)*fac
             am(i+1:,jlo) = am(i+1:,jlo) - lhscol(i+1:,m)*fac

             if (irow == i) then
                a_inv(:,jlo) = a_inv(:,jlo)/lhscol(i,m)
             else
                a_inv(:,jlo) = a_inv(:,jlo) &
                     - rhsrow(:,m)*lhscol(irow,m)/lhscol(i,m)
             end if
          else
             a_inv(:,jlo) = 0.0
          end if
   
       end do
    end do

! fill in skipped points for each field and supercell:
! Do not include internal ntgrid points in sum over supercell
    do i = 1, nidx*N_class(ic)
       iskip = N_class(ic) > 1
       iskip = i <= nidx*N_class(ic) - nfield .and. iskip
       iskip = mod((i+nfield-1)/nfield, 2*ntgrid+1) == 0 .and. iskip
       iskip = i > nfield .and. iskip
       if (iskip) then
          a_inv(i,:) = 0
          cycle
       end if
    end do
! Make response at internal ntgrid points identical to response
! at internal -ntgrid points:
    do jlo = f_lo(ic)%llim_world, f_lo(ic)%ulim_world
       jskip = N_class(ic) > 1
       jskip = ig_idx(f_lo(ic), jlo) == ntgrid .and. jskip
       jskip = in_idx(f_lo(ic), jlo) < N_class(ic) .and. jskip
       if (jskip) then
          ilo = jlo + nfield
          if (idx_local(f_lo(ic), ilo)) then
             if (idx_local(f_lo(ic), jlo)) then
                a_inv(:,jlo) = a_inv(:,ilo)
             else
                call send(a_inv(:,ilo), proc_id(f_lo(ic), jlo))
             endif
          else
             if (idx_local(f_lo(ic), jlo)) then
                call receive(a_inv(:,jlo), proc_id(f_lo(ic), ilo))
             end if
          end if
       end if
       call barrier
    end do

    am = a_inv
    deallocate (a_inv)

! Re-sort this class of aminv for runtime application.  

    if (.not.allocated(aminv)) allocate (aminv(i_class))

! only need this large array for particular values of jlo.
! To save space, count how many this is and only allocate
! required space:

    dc = 0
! check all members of this class
    do ilo = f_lo(ic)%llim_world, f_lo(ic)%ulim_world

! find supercell coordinates
       m = im_idx(f_lo(ic), ilo)
       n = in_idx(f_lo(ic), ilo)

! find standard coordinates
       ig = ig_idx(f_lo(ic), ilo)
       if = ifield_idx(f_lo(ic), ilo)
       ik = f_lo(ic)%ik(m,n)
       it = f_lo(ic)%it(m,n)

! translate to fast field coordinates
       jlo = ij_idx(jf_lo, ig, if, ik, it)
          
! Locate this jlo, count it, and save address
       if (idx_local(jf_lo,jlo)) then
! count it
          dc = dc + 1
! save dcell address
          dj(ic,jlo) = dc
! save supercell address
          mj(jlo) = m
       endif
          
    end do

! allocate dcells and supercells in this class on this PE:

    do jlo = jf_lo%llim_proc, jf_lo%ulim_proc
          
       if (.not.associated(aminv(ic)%dcell)) then
          allocate (aminv(ic)%dcell(dc))
       else
          j = size(aminv(ic)%dcell)
          if (j /= dc) then
             ierr = error_unit()
             write(ierr,*) 'Error (1) in init_inverse_matrix: ',&
                  iproc,':',jlo,':',dc,':',j
          endif
       endif
       
       k = dj(ic,jlo)
       if (k > 0) then
          jc = nidx*N_class(ic)
          if (.not.associated(aminv(ic)%dcell(k)%supercell)) then
             allocate (aminv(ic)%dcell(k)%supercell(jc))
          else
             j = size(aminv(ic)%dcell(k)%supercell)
             if (j /= jc) then
                ierr = error_unit()
                write(ierr,*) 'Error (2) in init_inverse_matrix: ', &
                     iproc,':',jlo,':',jc,':',j
             end if
          end if
       end if
    end do

! Now fill aminv for this class:

    allocate (am_tmp(nidx*N_class(ic)))

    do ilo = f_lo(ic)%llim_world, f_lo(ic)%ulim_world

       m = im_idx(f_lo(ic), ilo)
       n = in_idx(f_lo(ic), ilo)
       
       ig = ig_idx(f_lo(ic), ilo)
       if = ifield_idx(f_lo(ic), ilo)
       ik = f_lo(ic)%ik(m,n)
       it = f_lo(ic)%it(m,n)
       
       iflo = ij_idx(jf_lo, ig, if, ik, it)
 
       if (idx_local(f_lo(ic),ilo)) then
          if (idx_local(jf_lo,iflo)) then
             am_tmp = am(:,ilo)
          else
             call send(am(:,ilo), proc_id(jf_lo,iflo))
          endif
       else
          if (idx_local(jf_lo,iflo)) then
             call receive(am_tmp, proc_id(f_lo(ic),ilo))
          end if
       end if

       if (idx_local(jf_lo, iflo)) then

          dc = dj(ic,iflo)

          do jlo = 0, nidx*N_class(ic)-1
             
             nn = in_idx(f_lo(ic), jlo)
             
             jg = ig_idx(f_lo(ic), jlo)
             jf = ifield_idx(f_lo(ic), jlo)
             
             jsc = ij_idx(f_lo(ic), jg, jf, nn) + 1

             aminv(ic)%dcell(dc)%supercell(jsc) = am_tmp(jlo+1)
             
          end do
       end if
       call barrier
    end do

    deallocate (am_tmp)

    call prof_leaving ("init_inverse_matrix", "fields_implicit")
  end subroutine init_inverse_matrix

  subroutine finish_fields_layouts

    use dist_fn, only: N_class, i_class, itright, boundary
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: f_lo, jf_lo, ij, ik_idx, it_idx
    implicit none
    integer :: i, m, n, ii, ik, it, itr, jflo

    call boundary(linked)
    if (linked) then

! Complication comes from having to order the supercells in each class
       do ii = 1, i_class
          m = 1
          do it = 1, ntheta0
             do ik = 1, naky
                call kt2ki (i, n, ik, it)
                ! If (ik, it) is in this class, continue:
                if (i == ii) then
                   ! Find left end of links
                   if (n == 1) then
                      f_lo(i)%ik(m,n) = ik
                      f_lo(i)%it(m,n) = it
                      itr = it
                      ! Follow links to the right
                      do n = 2, N_class(i)
                         itr = itright (ik, itr)
                         f_lo(i)%ik(m,n) = ik
                         f_lo(i)%it(m,n) = itr
                      end do
                      m = m + 1
                   end if
                end if
             end do
          end do
       end do
       
    ! initialize ij matrix
       
       do jflo = jf_lo%llim_proc, jf_lo%ulim_proc
          ik = ik_idx(jf_lo, jflo)
          it = it_idx(jf_lo, jflo)
          
          call kt2ki (ij(jflo), n, ik, it)
          
       end do

    else
       m = 0
       do it = 1, ntheta0
          do ik = 1, naky
             m = m + 1
             f_lo(1)%ik(m,1) = ik
             f_lo(1)%it(m,1) = it
          end do
       end do
       
       ij = 1
    end if

  end subroutine finish_fields_layouts

  subroutine kt2ki (i, n, ik, it)

    use file_utils, only: error_unit
    use dist_fn, only: l_links, r_links, N_class, i_class

    integer, intent (in) :: ik, it
    integer, intent (out) :: i, n

    integer :: nn, ierr
!
! Get size of this supercell
!
    nn = 1 + l_links(ik,it) + r_links(ik,it)
!
! Find i = N_class**-1(nn)
!
    do i = 1, i_class
       if (N_class(i) == nn) exit
    end do
!
! Consistency check:
!
    if (N_class(i) /= nn) then
       ierr = error_unit()
       write(ierr,*) 'Error in kt2ki:'
       write(ierr,*) 'i = ',i,' ik = ',ik,' it = ',it,&
            ' N(i) = ',N_class(i),' nn = ',nn
       stop
    end if
! 
! Get position in this supercell, counting from the left
!
    n = 1 + l_links(ik, it)

  end subroutine kt2ki

end module fields_implicit

