module fields_implicit
  use fields_arrays, only: nidx
  implicit none

  public :: init_fields_implicit
  public :: advance_implicit
  public :: init_phi_implicit
  public :: nidx
  public :: reset_init
  public :: time_field
  public :: set_scan_parameter

  private

  integer, save :: nfield
  logical :: initialized = .false.
  logical :: linked = .false.
  real, save :: time_field(2)=0.

contains

  subroutine init_fields_implicit
    use antenna, only: init_antenna
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids
    use gs2_layouts, only: init_gs2_layouts
    use parameter_scan_arrays, only: run_scan
    implicit none
    logical:: debug=.false.
    logical :: dummy

    if (initialized) return
    initialized = .true.

    if (debug) write(6,*) "init_fields_implicit: gs2_layouts"
    call init_gs2_layouts
    if (debug) write(6,*) "init_fields_implicit: theta_grid"
    call init_theta_grid
    if (debug) write(6,*) "init_fields_implicit: kt_grids"
    call init_kt_grids
    if (debug) write(6,*) "init_fields_implicit: read_parameters"
    call read_parameters
    if (debug .and. run_scan) &
        write(6,*) "init_fields_implicit: set_scan_parameter"
        ! Must be done before resp. m.
        if (run_scan) call set_scan_parameter(dummy)
    if (debug) write(6,*) "init_fields_implicit: response_matrix"
    call init_response_matrix
    if (debug) write(6,*) "init_fields_implicit: antenna"
    call init_antenna

  end subroutine init_fields_implicit

  
  subroutine set_scan_parameter(reset)
    !use parameter_scan_arrays, only: current_scan_parameter_value
    !use parameter_scan_arrays, only: scan_parameter_switch
    !use parameter_scan_arrays, only: scan_parameter_tprim
    !use parameter_scan_arrays, only: scan_parameter_g_exb
    use parameter_scan_arrays
    use species, only: spec 
    use dist_fn, only: g_exb
    use mp, only: proc0
    logical, intent (inout) :: reset
     
    select case (scan_parameter_switch)
    case (scan_parameter_tprim)
       spec(scan_spec)%tprim = current_scan_parameter_value
       if (proc0) write (*,*) &
         "Set scan parameter tprim_1 to ", spec(scan_spec)%tprim
       reset = .true.
    case (scan_parameter_g_exb)
       g_exb = current_scan_parameter_value
       if (proc0) write (*,*) &
         "Set scan parameter g_exb to ", g_exb
       reset = .false.
    end select
  end subroutine set_scan_parameter

  subroutine read_parameters
  end subroutine read_parameters

  subroutine init_phi_implicit

    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use fields_arrays, only: phip, aparp, bparp, phipnew, aparpnew, bparpnew
    use fields_arrays, only: phih, aparh, bparh, phihnew, aparhnew, bparhnew
    use dist_fn_arrays, only: g, gnew, gpnew, gpold !ghnew, ghold
    use dist_fn, only: get_init_field
    use init_g, only: new_field_init

    implicit none

    if (new_field_init) then
       ! get iniitial phik corresponding to lowest order gk
       call get_init_field (gnew, phinew, aparnew, bparnew)
       phi = phinew; apar = aparnew; bpar = bparnew; g = gnew

       ! get initial dphi/drho corresponding to dgk/drho
       call get_init_field (gpnew, phipnew, aparpnew, bparpnew)
       phip = phipnew ; aparp = aparpnew ; bparp = bparpnew ; gpold = gpnew

       ! get iniitial phik corresponding to lowest order gk
! TMP UNTIL phih needed
!       call get_init_field (ghnew, phihnew, aparhnew, bparhnew)
       phihnew = 0. ; aparhnew = 0. ; bparhnew = 0.
       phih = phihnew; aparh = aparhnew; bparh = bparhnew! ; ghold = ghnew
    else
       ! get initial phik corresponding to gk
       call getfield (gnew, phinew, aparnew, bparnew)
       phi = phinew; apar = aparnew; bpar = bparnew

       ! get initial dphi/drho corresponding to dgk/drho
       call getfield (gpnew, phipnew, aparpnew, bparpnew)
       phip = phipnew; aparp = aparpnew; bparp = bparpnew

       ! get iniitial phik corresponding to lowest order gk
!       call getfield (ghnew, phihnew, aparhnew, bparhnew)
       phihnew = 0. ; aparhnew = 0. ; bparhnew = 0.
       phih = phihnew; aparh = aparhnew; bparh = bparhnew
    end if

  end subroutine init_phi_implicit

  subroutine get_field_vector (fl, gfnc, phi, apar, bpar)

    use gs2_layouts, only: g_lo
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use dist_fn, only: getfieldeq
    use run_parameters, only: fphi, fapar, fbpar
    use prof, only: prof_entering, prof_leaving
    use vpamu_grids, only: nvgrid

    implicit none

    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in) :: gfnc
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (:,:,:), intent (out) :: fl
    complex, dimension (:,:,:), allocatable :: fieldeq, fieldeqa, fieldeqp
    integer :: istart, ifin

    call prof_entering ("get_field_vector", "fields_implicit")

    allocate (fieldeq (-ntgrid:ntgrid,ntheta0,naky))
    allocate (fieldeqa(-ntgrid:ntgrid,ntheta0,naky))
    allocate (fieldeqp(-ntgrid:ntgrid,ntheta0,naky))

    call getfieldeq (gfnc, phi, apar, bpar, fieldeq, fieldeqa, fieldeqp)

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

    if (fbpar > epsilon(0.0)) then
       istart = ifin + 1
       ifin = (istart-1) + 2*ntgrid+1
       fl(istart:ifin,:,:) = fieldeqp
    end if

    deallocate (fieldeq, fieldeqa, fieldeqp)

    call prof_leaving ("get_field_vector", "fields_implicit")

  end subroutine get_field_vector

  subroutine get_field_solution (u, phifnc, aparfnc, bparfnc)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: fphi, fapar, fbpar
    use gs2_layouts, only: jf_lo, ij_idx
    use prof, only: prof_entering, prof_leaving

    implicit none

    complex, dimension (0:), intent (in) :: u
    complex, dimension (-ntgrid:,:,:), intent (out) :: phifnc, aparfnc, bparfnc

    integer :: ik, it, ifield, ll, lr

    call prof_entering ("get_field_solution", "fields_implicit")

    ifield = 0

    if (fphi > epsilon(0.0)) then
       ifield = ifield + 1
       do ik = 1, naky
          do it = 1, ntheta0
             ll = ij_idx (jf_lo, -ntgrid, ifield, ik, it)
             lr = ll + 2*ntgrid
             phifnc(:,it,ik) = u(ll:lr)
          end do
       end do
    endif

    if (fapar > epsilon(0.0)) then
       ifield = ifield + 1
       do ik = 1, naky
          do it = 1, ntheta0
             ll = ij_idx (jf_lo, -ntgrid, ifield, ik, it)
             lr = ll + 2*ntgrid
             aparfnc(:,it,ik) = u(ll:lr)
          end do
       end do
    endif

    if (fbpar > epsilon(0.0)) then
       ifield = ifield + 1
       do ik = 1, naky
          do it = 1, ntheta0
             ll = ij_idx (jf_lo, -ntgrid, ifield, ik, it)
             lr = ll + 2*ntgrid
             bparfnc(:,it,ik) = u(ll:lr)
          end do
       end do
    endif

    call prof_leaving ("get_field_solution", "fields_implicit")

  end subroutine get_field_solution

  subroutine getfield (gfnc, phi, apar, bpar)

    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: f_lo, jf_lo, ij, mj, dj, g_lo
    use prof, only: prof_entering, prof_leaving
    use fields_arrays, only: aminv
    use theta_grid, only: ntgrid
    use dist_fn, only: N_class
    use mp, only: sum_allreduce, proc0
    use job_manage, only: time_message
    use vpamu_grids, only: nvgrid

    implicit none

    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in) :: gfnc
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi, apar, bpar
    complex, dimension (:,:,:), allocatable :: fl
    complex, dimension (:), allocatable :: u
    integer :: jflo, ik, it, nl, nr, i, m, n, dc

    if (proc0) call time_message(.false.,time_field,' Field Solver')

    call prof_entering ("getfield", "fields_implicit")
    allocate (fl(nidx, ntheta0, naky))
    allocate (u (0:nidx*ntheta0*naky-1))

    ! am*u = fl, Poisson's and Ampere's law, u is phi, apar, bpar 
    ! u = aminv*fl

    ! setup fl to contain the field equations, which should be zero
    ! for the right choice of phi, apar, bpar.  They will be nonzero
    ! though because they have been evaluated using old values of phi, apar, bpar.
    call get_field_vector (fl, gfnc, phi, apar, bpar)

!    write (*,*) 'fl', fl(:,1,naky)

    u = 0.
    do jflo = jf_lo%llim_proc, jf_lo%ulim_proc

       ! i here identifies the 2*pi segment
       i = ij(jflo)
       m = mj(jflo)
       ik = f_lo(i)%ik(m,1)  ! For fixed i and m, ik does not change as n varies 
       dc = dj(i,jflo)
       
       ! N_class is the number of (2*ntgrid+1) segments for this ky
       ! needed to form extended grid in theta
       do n = 1, N_class(i)
                    
          it = f_lo(i)%it(m,n)
          
          nl = 1 + nidx*(n-1)
          nr = nl + nidx - 1

          u(jflo)=u(jflo)-sum(aminv(i)%dcell(dc)%supercell(nl:nr)*fl(:, it, ik)) 

       end do
    end do

    deallocate (fl)
    call sum_allreduce (u)

    call get_field_solution (u, phi, apar, bpar)

    deallocate (u)

    call prof_leaving ("getfield", "fields_implicit")

    if (proc0) call time_message(.false.,time_field,' Field Solver')

  end subroutine getfield

  subroutine advance_implicit (istep, remove_zonal_flows_switch)

    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use fields_arrays, only: phip, aparp, bparp, phipnew, aparpnew, bparpnew
!    use fields_arrays, only: phih, aparh, bparh, phihnew, aparhnew, bparhnew
    use fields_arrays, only: apar_ext !, phi_ext
    use antenna, only: antenna_amplitudes
    use dist_fn, only: timeadv!, get_omega_prime, omprim_converge, get_total_g
    use dist_fn, only: get_fieldcorrection, profile_variation
    use dist_fn_arrays, only: g, gold, gnew
    use dist_fn_arrays, only: gpold, gpnew!, ghold, ghnew
!    use mp, only: proc0
!    use kt_grids, only: naky

    implicit none

    integer :: diagnostics = 1
    integer, intent (in) :: istep
    logical, intent (in) :: remove_zonal_flows_switch

    !GGH NOTE: apar_ext is initialized in this call
    call antenna_amplitudes (apar_ext)

    ! ! need to worry about how this should be done if splitting g1 and g2
    ! if (allocated(kx_shift) .or. allocated(theta0_shift)) then
    !    call exb_shear (gnew, phinew, aparnew, bparnew)
    !    call exb_shear (gpnew, phipnew, aparpnew, bparpnew)
    !    call exb_shear (ghnew, phihnew, aparhnew, bparhnew)
    ! end if

    ! TMP FOR TESTING -- MAB
!    phinew = 0.
    g = gnew ; gold = gnew ; gpold = gpnew !; ghold = ghnew
    phi = phinew ; phip = phipnew !; phih = phihnew
    apar = aparnew ; aparp = aparpnew !; aparh = aparhnew
    bpar = bparnew ; bparp = bparpnew !; bparh = bparhnew

    ! first solve for g^{n+1} and phi^{n+1} using g'^{n} and phi'^{n}

    ! get lowest order gk^{*}. the ^{*} indicates evaluation of g^{n+1} obtained with phi=phi^{n}
    if (profile_variation) then
       call timeadv (gnew, gold, phi, apar, bpar, phinew, aparnew, istep, equation=2)
    else
       call timeadv (gnew, gold, phi, apar, bpar, phinew, aparnew, istep, equation=0)
    end if

!    if (proc0) write (*,*) 'phi1', istep, real(phinew(0,1,naky)), aimag(phinew(0,1,naky))

    aparnew = aparnew + apar_ext

    ! get phik, apark, bpark
    call getfield (gnew, phinew, aparnew, bparnew)
! TMP FOR TESTING -- MAB
!    phinew = 0.
    phinew = phinew + phi
    aparnew = aparnew + apar
    bparnew = bparnew + bpar

!    if (proc0) write (*,*) 'phi2', istep, real(phinew(0,1,naky)), aimag(phinew(0,1,naky))

    if (remove_zonal_flows_switch) call remove_zonal_flows
    
    ! get lowest order gk^{n+1}
    if (profile_variation) then
       call timeadv (gnew, gold, phi, apar, bpar, phinew, &
            aparnew, istep, equation=2, mode=diagnostics)
    else
       call timeadv (gnew, gold, phi, apar, bpar, phinew, &
            aparnew, istep, equation=0, mode=diagnostics)
    end if

!    if (proc0) write (*,*) 'phi3', istep, real(phinew(0,1,naky)), aimag(phinew(0,1,naky))

    if (profile_variation) then
       ! next solve for g'^{n+1} and phi'^{n+1} using g^{n+1} and phi^{n+1}

       ! get (dgk/drho)^{*}, which is is (dgk/drho)^{n+1} obtained with (dphik/drho)=(dphik/drho)^{n}
       call timeadv (gpnew, gpold, phip, aparp, bparp, phipnew, &
            aparpnew, istep, equation=1)
       ! get next order gk^{*}. the ^{*} indicates evaluation of g^{n+1} obtained with phi=phi^{n}
       !    call timeadv (ghnew, ghold, phih, aparh, bparh, phihnew, aparhnew, bparhnew, istep, equation=2)

       ! get phipk, aparpk, bparpk
       call getfield (gpnew, phipnew, aparpnew, bparpnew)
       ! get phihk, aparhk, bparhk
       !   call getfield (ghnew, phihnew, aparhnew, bparhnew)
       
       phipnew = phipnew + phip
       aparpnew = aparpnew + aparp
       bparpnew = bparpnew + bparp
       
       call get_fieldcorrection (gnew, phinew, phipnew)
       
       ! get (dgk/drho)^{n+1}
       call timeadv (gpnew, gpold, phip, aparp, bparp, phipnew, &
            aparpnew, istep, equation=1)
       
       ! correct lowest order g with profile variation piece
!       call get_total_g

       ! get next order gk^{n+1}
       !    call timeadv (ghnew, ghold, phih, aparh, bparh, phihnew, aparhnew, bparhnew, istep, equation=2)
       
       !    call get_omega_prime (phipnew/phinew, phip/phi, phihnew/phinew, phih/phi)
    end if


  end subroutine advance_implicit

  subroutine remove_zonal_flows
    use fields_arrays, only: phinew
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    
    complex, dimension(:,:,:), allocatable :: phi_avg

    allocate(phi_avg(-ntgrid:ntgrid,ntheta0,naky)) 
    phi_avg = 0.
    ! fieldline_average_phi will calculate the field line average of phinew and 
    ! put it into phi_avg, but only for ik = 1 (the last parameter of the call)
    call fieldline_average_phi(phinew, phi_avg, 1)
    phinew = phinew - phi_avg
    deallocate(phi_avg)
  end subroutine remove_zonal_flows

  !> This generates a field line average of phi_in and writes it to 
  !! phi_average. If ik_only is supplied, it will only calculate the
  !! field line average for that ky, leaving the rest of phi_avg unchanged. EGH
  
  ! It replaces the routines fieldlineavgphi_loc and fieldlineavgphi_tot,
  ! in fields.f90, which I  think are defunct, as phi is always on every processor.

  subroutine fieldline_average_phi (phi_in, phi_average, ik_only)
    use theta_grid, only: ntgrid, drhodpsi, gradpar, bmag, delthet
    use kt_grids, only: ntheta0, naky

    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi_in
    complex, dimension (-ntgrid:,:,:), intent (out) :: phi_average
    integer, intent (in), optional :: ik_only
    real, dimension (-ntgrid:ntgrid) :: jac
    !complex, dimension (-ntgrid:ntgrid) :: phi_avg_line
    complex :: phi_avg_line
    integer it, ik, ik_only_actual
    ik_only_actual = -1
    if (present(ik_only)) ik_only_actual = ik_only

    jac = 1.0/abs(drhodpsi*gradpar*bmag)
    if (ik_only_actual .gt. 0) then
      do it = 1,ntheta0
         phi_avg_line = sum(phi_in(-ntgrid:ntgrid,it,ik_only_actual)* &
            jac(-ntgrid:ntgrid)*delthet(-ntgrid:ntgrid))/ &
            sum(delthet(-ntgrid:ntgrid)*jac(-ntgrid:ntgrid))
           phi_average(:, it, ik_only_actual) = phi_avg_line
      end do
    else
      do it = 1,ntheta0
        do ik = 1,naky
          phi_average(:, it, ik) = sum(phi_in(-ntgrid:ntgrid,it,ik)*jac*delthet)/sum(delthet*jac)
        end do
      end do
    end if

  end subroutine fieldline_average_phi

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
    use mp, only: barrier!, proc0
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use dist_fn_arrays, only: g
    use dist_fn, only: M_class, N_class, i_class
    use run_parameters, only: fphi, fapar, fbpar
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
    if (fbpar > epsilon(0.0)) nfield = nfield + 1
    nidx = (2*ntgrid+1)*nfield

    call init_fields_layouts (nfield, nidx, naky, ntheta0, M_class, N_class, i_class)
    call init_jfields_layouts (nfield, nidx, naky, ntheta0, i_class)
    call finish_fields_layouts
!
! keep storage cost down by doing one class at a time
! Note: could define a superclass (of all classes), a structure containing all am, 
! then do this all at once.  This would be faster, especially for large runs in a 
! sheared domain, and could be triggered by local_field_solve 
! 
    do i = i_class, 1, -1

       call barrier
!       if (proc0) write(*,*) 'beginning class ',i,' with size ',nidx*N_class(i)
       allocate (am(nidx*N_class(i), f_lo(i)%llim_proc:f_lo(i)%ulim_alloc))

       am = 0.0
       g = 0.0
       
       phi = 0.0
       apar = 0.0
       bpar = 0.0
       phinew = 0.0
       aparnew = 0.0
       bparnew = 0.0

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
             
             if (fbpar > epsilon(0.0)) then
                ifield = ifield + 1
                if (endpoint) then
                   do m = 1, M_class(i)
                      ik = f_lo(i)%ik(m,n-1)
                      it = f_lo(i)%it(m,n-1)
                      bparnew(ntgrid,it,ik) = 1.0
                   end do
                endif
                do m = 1, M_class(i)
                   ik = f_lo(i)%ik(m,n)
                   it = f_lo(i)%it(m,n)
                   bparnew(ig,it,ik) = 1.0
                end do
                call init_response_row (ig, ifield, am, i, n)
                bparnew = 0.0
             end if
          end do
       end do

       call init_inverse_matrix (am, i)
       deallocate (am)
!       if (proc0) write(*,*) 'finished class ',i

    end do

    call prof_leaving ("init_response_matrix", "fields_implicit")

  end subroutine init_response_matrix

  subroutine init_response_row (ig, ifield, am, ic, n)

    use dist_fn_arrays, only: gnew, gold
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use dist_fn, only: getfieldeq, timeadv, M_class, N_class
    use run_parameters, only: fphi, fapar, fbpar
    use gs2_layouts, only: f_lo, idx, idx_local
    use prof, only: prof_entering, prof_leaving

    implicit none

    integer, intent (in) :: ig, ifield, ic, n
    complex, dimension(:,f_lo(ic)%llim_proc:), intent (in out) :: am
    complex, dimension (:,:,:), allocatable :: fieldeq, fieldeqa, fieldeqp
    integer :: irow, istart, iflo, ik, it, ifin, m, nn

    call prof_entering ("init_response_row", "fields_implicit")

    allocate (fieldeq (-ntgrid:ntgrid, ntheta0, naky))
    allocate (fieldeqa(-ntgrid:ntgrid, ntheta0, naky))
    allocate (fieldeqp(-ntgrid:ntgrid, ntheta0, naky))

    ! get response matrix using lowest order equation for g
    ! no need to get different response matrix for higher order equations
    ! since the form of the implicit solve is the same for all of them
    call timeadv (gnew, gold, phi, apar, bpar, phinew, &
         aparnew, 0, equation=0)

    call getfieldeq (gnew, phinew, aparnew, bparnew, fieldeq, fieldeqa, fieldeqp)

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
             
             if (fbpar > epsilon(0.0)) then
                ifin = istart + nidx
                istart = istart + 1
                am(istart:ifin:nfield,iflo) = fieldeqp(:,it,ik)
             end if
             
          end if
                    
       end do
    end do
    deallocate (fieldeq, fieldeqa, fieldeqp)
    call prof_leaving ("init_response_row", "fields_implicit")
  end subroutine init_response_row

  subroutine init_inverse_matrix (am, ic)

    use file_utils, only: error_unit
    use kt_grids, only: aky, akx
    use theta_grid, only: ntgrid
    use mp, only: broadcast, send, receive, barrier, iproc
    use gs2_layouts, only: f_lo, idx, idx_local, proc_id, jf_lo
    use gs2_layouts, only: if_idx, im_idx, in_idx, local_field_solve
    use gs2_layouts, only: ig_idx, ifield_idx, ij_idx, mj, dj
    use prof, only: prof_entering, prof_leaving
    use fields_arrays, only: aminv
    use dist_fn, only: i_class, M_class, N_class

    implicit none

    integer, intent (in) :: ic
    complex, dimension(:,f_lo(ic)%llim_proc:), intent (in out) :: am
    complex, dimension(:,:), allocatable :: a_inv, lhscol, rhsrow
    complex, dimension (:), allocatable :: am_tmp
    complex :: fac
    integer :: i, j, k, ik, it, m, n, nn, if, ig, jsc, jf, jg, jc
    integer :: irow, ilo, jlo, dc, iflo, ierr
    logical :: iskip, jskip

    call prof_entering ("init_inverse_matrix", "fields_implicit")
    
    allocate (lhscol (nidx*N_class(ic),M_class(ic)))
    allocate (rhsrow (nidx*N_class(ic),M_class(ic)))
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
 
       if (local_field_solve) then
          do m = 1, M_class(ic)
             ilo = idx(f_lo(ic),i,m)
             if (idx_local(f_lo(ic),ilo)) then
                lhscol(:,m) = am(:,ilo)
                rhsrow(:,m) = a_inv(:,ilo)
             end if
          end do
       else
          do m = 1, M_class(ic)
             ilo = idx(f_lo(ic),i,m)
             if (idx_local(f_lo(ic),ilo)) then
                lhscol(:,m) = am(:,ilo)
                rhsrow(:,m) = a_inv(:,ilo)
             end if
             call broadcast (lhscol(:,m), proc_id(f_lo(ic),ilo))
             call broadcast (rhsrow(:,m), proc_id(f_lo(ic),ilo))
          end do
       end if

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

    deallocate (lhscol, rhsrow)

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

  subroutine timer
    
    character (len=10) :: zdate, ztime, zzone
    integer, dimension(8) :: ival
    real, save :: told=0., tnew=0.
    
    call date_and_time (zdate, ztime, zzone, ival)
    tnew = ival(5)*3600.+ival(6)*60.+ival(7)+ival(8)/1000.
    if (told > 0.) then
       print *, 'Fields_implicit: Time since last called: ',tnew-told,' seconds'
    end if
    told = tnew
  end subroutine timer

end module fields_implicit

