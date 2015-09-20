module fields_implicit
  use fields_arrays, only: nidx
  implicit none

  public :: init_fields_implicit
  public :: advance_implicit
  public :: remove_zonal_flows
  public :: init_allfields_implicit
  public :: nidx
  public :: reset_init
  public :: field_subgath, dump_response, read_response
  public :: dump_response_to_file_imp

  !> The number of steps between recalculating
  !! the response matrix for the alternative flow shear implementation. EGH
  public :: n_recalc_response     

  !> Unit tests
  public :: fields_implicit_unit_test_init_fields_implicit

  private

  !> A variable to help with running benchmarks... do not set true
  !! unless you know what you are doing. If true, the response matrix
  !! will not be initialised and set to zero. The results of any 
  !! simulation will be garbage
  logical, public :: skip_initialisation = .false.

  integer, save :: nfield
  logical :: initialized = .false.
  logical :: linked = .false.
  logical :: field_subgath
  logical :: dump_response=.false., read_response=.false.

  !> For running with the explicit flow shear implementation
  !! when you want to run using some implicitness 
  integer, save :: g_exb_error_check_cycle = 0
  integer, save:: n_recalc_response = 1
contains

  subroutine init_fields_implicit
    use antenna, only: init_antenna
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids
    use gs2_layouts, only: init_gs2_layouts
!    use parameter_scan_arrays, only: run_scan
    implicit none
    logical, parameter :: debug=.false.
!    logical :: dummy

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
 !   if (debug .and. run_scan) &
 !       write(6,*) "init_fields_implicit: set_scan_parameter"
        ! Must be done before resp. m.
        !if (run_scan) call set_scan_parameter(dummy)
    if (debug) write(6,*) "init_fields_implicit: response_matrix"
    call init_response_matrix
    if (debug) write(6,*) "init_fields_implicit: antenna"
    call init_antenna
  end subroutine init_fields_implicit

  function fields_implicit_unit_test_init_fields_implicit()
    logical :: fields_implicit_unit_test_init_fields_implicit

    call init_fields_implicit

    fields_implicit_unit_test_init_fields_implicit = .true.

  end function fields_implicit_unit_test_init_fields_implicit

  subroutine read_parameters
  end subroutine read_parameters

  subroutine init_allfields_implicit
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use dist_fn_arrays, only: g, gnew
    use dist_fn, only: get_init_field
    use init_g, only: new_field_init

    implicit none

    ! MAB> new field init option ported from agk
    if (new_field_init) then
       call get_init_field (phinew, aparnew, bparnew)
       phi = phinew; apar = aparnew; bpar = bparnew; g = gnew
    else
       call getfield (phinew, aparnew, bparnew)
       phi = phinew; apar = aparnew; bpar = bparnew
    end if
    ! <MAB

  end subroutine init_allfields_implicit

  subroutine get_field_vector (fl, phi, apar, bpar)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use dist_fn, only: getfieldeq
    use run_parameters, only: fphi, fapar, fbpar
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (:,:,:), intent (out) :: fl
    complex, dimension (:,:,:), allocatable :: fieldeq, fieldeqa, fieldeqp
    integer :: istart, ifin

    call prof_entering ("get_field_vector", "fields_implicit")

    allocate (fieldeq (-ntgrid:ntgrid,ntheta0,naky))
    allocate (fieldeqa(-ntgrid:ntgrid,ntheta0,naky))
    allocate (fieldeqp(-ntgrid:ntgrid,ntheta0,naky))

    call getfieldeq (phi, apar, bpar, fieldeq, fieldeqa, fieldeqp)

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

  subroutine get_field_solution (u)
    use fields_arrays, only: phinew, aparnew, bparnew
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: fphi, fapar, fbpar
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

    if (fbpar > epsilon(0.0)) then
       ifield = ifield + 1
       do ik = 1, naky
          do it = 1, ntheta0
             ll = ij_idx (jf_lo, -ntgrid, ifield, ik, it)
             lr = ll + 2*ntgrid
             bparnew(:,it,ik) = u(ll:lr)
          end do
       end do
    endif

    call prof_leaving ("get_field_solution", "fields_implicit")
  end subroutine get_field_solution

  subroutine getfield (phi, apar, bpar)
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: f_lo, jf_lo, ij, mj, dj
    use prof, only: prof_entering, prof_leaving
    use fields_arrays, only: aminv, time_field
    use theta_grid, only: ntgrid
    use dist_fn, only: N_class
    use mp, only: sum_allreduce, allgatherv, iproc,nproc, proc0
    use job_manage, only: time_message
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (:,:,:), allocatable :: fl
    complex, dimension (:), allocatable :: u
    complex, dimension (:), allocatable :: u_small
    integer :: jflo, ik, it, nl, nr, i, m, n, dc
    integer, dimension(:), allocatable,save :: recvcnts, displs

    if (proc0) call time_message(.false.,time_field,' Field Solver')

    call prof_entering ("getfield", "fields_implicit")
    allocate (fl(nidx, ntheta0, naky))

    !On first call to this routine setup the receive counts (recvcnts)
    !and displacement arrays (displs)
    if ((.not.allocated(recvcnts)).and.field_subgath) then
       allocate(recvcnts(nproc),displs(nproc)) !Note there's no matching deallocate
       do i=0,nproc-1
          displs(i+1)=MIN(i*jf_lo%blocksize,jf_lo%ulim_world+1) !This will assign a displacement outside the array for procs with no data
          recvcnts(i+1)=MIN(jf_lo%blocksize,jf_lo%ulim_world-displs(i+1)+1) !This ensures that we expect no data from procs without any
       enddo
    endif

    ! am*u = fl, Poisson's and Ampere's law, u is phi, apar, bpar 
    ! u = aminv*fl

    call get_field_vector (fl, phi, apar, bpar)

    !Initialise array, if not gathering then have to zero entire array
    if(field_subgath) then
       allocate(u_small(jf_lo%llim_proc:jf_lo%ulim_proc))
    else
       allocate(u_small(0:nidx*ntheta0*naky-1))
    endif
    u_small=0.

    !Should this really be to ulim_alloc instead?
    do jflo = jf_lo%llim_proc, jf_lo%ulim_proc
       
       !Class index
       i = ij(jflo)
       
       !Class member index (i.e. which member of the class)
       m = mj(jflo)

       !Get ik index
       ik = f_lo(i)%ik(m,1)  ! For fixed i and m, ik does not change as n varies 

       !Get d(istributed) cell index
       dc = dj(i,jflo)
       
       !Loop over cells in class (these are the 2pi domains in flux tube/box mode)
       do n = 1, N_class(i)
          
          !Get it index
          it = f_lo(i)%it(m,n)
          
          !Get extent of current cell in extended/ballooning space domain
          nl = 1 + nidx*(n-1)
          nr = nl + nidx - 1
          
          !Perform section of matrix vector multiplication
          u_small(jflo)=u_small(jflo)-sum(aminv(i)%dcell(dc)%supercell(nl:nr)*fl(:, it, ik)) 
          
       end do
    end do

    !Free memory
    deallocate (fl)

    !Gather/reduce the remaining data
    if(field_subgath) then
       allocate (u (0:nidx*ntheta0*naky-1))
       call allgatherv(u_small,recvcnts(iproc+1),u,recvcnts,displs)
       deallocate(u_small)
    else
       call sum_allreduce(u_small)
    endif

    !Reshape data into field arrays and free memory
    if(field_subgath)then
       call get_field_solution (u)
       deallocate(u)
    else
       call get_field_solution (u_small)
       deallocate(u_small)
    endif

    !For profiling
    call prof_leaving ("getfield", "fields_implicit")

    !For timing
    if (proc0) call time_message(.false.,time_field,' Field Solver')

  end subroutine getfield

  subroutine exb_shear (istep)
    use dist_fn, only:  exb_shear_d => exb_shear, g_exb, g_exb_error_limit
    use dist_fn, only:  g_exb_start_timestep, g_exb_start_time
    use dist_fn, only:  init_bessel, init_fieldeq
    use dist_fn_arrays, only: gnew, g_store
    use dist_fn_arrays, only: theta0_shift, kx_shift
    use fields_arrays, only: phinew, aparnew, bparnew
    use fields_arrays, only: phi, apar, bpar
    use fields_arrays, only: phi_store, apar_store, bpar_store
    use kt_grids, only: single, naky, ntheta0, akx
    use kt_grids, only: theta0, aky 
    use kt_grids, only: calculate_kt_grids
    use theta_grid, only: ntgrid, delthet, jacob
    use gs2_time, only: code_time, code_dt, code_dt_old
    use mp, only: proc0, iproc

    integer, intent (in) :: istep
    logical :: recalc_response
    logical,save :: check_g_exb_error = .false.
    complex, dimension(:, :, :), allocatable :: phi_temp, apar_temp, bpar_temp
    real, dimension (-ntgrid:ntgrid-1) :: wgt
    real, save :: wait_time = 0.
    logical, save :: check_error = .false.
    real :: phi_error, anorm, phitot, gdt
    integer :: ik, it
    !integer, parameter :: max_n = 1000 
    integer :: max_n  
    if (.not. single) then
      if (allocated(kx_shift) .or. allocated(theta0_shift)) call exb_shear_d (gnew, phinew, aparnew, bparnew, istep)                             ! See Hammett & Loureiro, APS 2006
      return
    end if

    if (abs(g_exb)<epsilon(0.0)) return

    gdt = 0.5*(code_dt + code_dt_old)

    if (single) then
      allocate(phi_temp(-ntgrid:ntgrid,ntheta0,naky))
      allocate(apar_temp(-ntgrid:ntgrid,ntheta0,naky))
      allocate(bpar_temp(-ntgrid:ntgrid,ntheta0,naky))
      if (g_exb_start_timestep > 0) then 
        if (istep < g_exb_start_timestep) then 
          return
        else if (g_exb_start_timestep == istep) then
          wait_time = code_time
        end if 
      end if 
      if (g_exb_start_time >= 0) then 
        if (code_time < g_exb_start_time) then 
          return
        else if (wait_time .eq. 0.0) then
          wait_time = code_time
        end if 
      end if

      ! if g_exb_error_limit < 0, always recalc response matrices
      if (g_exb_error_limit .lt. 0.0) g_exb_error_check_cycle = 2

      if (g_exb_error_check_cycle .eq. 3) then ! test error
        ! phi_store was calculated without updating response matrix
        phi_temp = phinew - phi_store
        !phi_temp = phinew*conjg(phinew) - phi_store*conjg(phi_store)
        phi_error = 0
        wgt = delthet(-ntgrid:ntgrid-1)*jacob(-ntgrid:ntgrid-1)  
        anorm = sum(wgt)
        do ik = 1, naky
           do it = 1, ntheta0
              phitot =  sum( real(conjg(phinew(-ntgrid:ntgrid-1,it,ik)) &
               * phinew(-ntgrid:ntgrid-1,it,ik)) * wgt ) / anorm   
              phi_error = phi_error +  (sum( real(conjg(phi_temp(-ntgrid:ntgrid-1,it,ik)) &
               * phi_temp(-ntgrid:ntgrid-1,it,ik)) * wgt ) / anorm) / phitot     
              !phi_error = phi_error +  (sum( real(phi_temp(-ntgrid:ntgrid-1,it,ik)) &
                !* wgt ) / anorm) / phitot     
           end do
        end do
        if (proc0) write (*,*) "phi error with no recalc: ", phi_error

        max_n = max(floor(1.0/abs(g_exb)/gdt) / 10, 1)

        if (abs(phi_error) .lt. 1.0e-100) then
          n_recalc_response = max_n 
          ! Maximum number of timesteps without rechecking
        else
          n_recalc_response = floor((g_exb_error_limit / phi_error)**(1.0/3.0))
        end if
        if (n_recalc_response .eq. 0) n_recalc_response = 1
        if (n_recalc_response .gt. max_n) n_recalc_response = max_n
        g_exb_error_check_cycle = 0 ! Normal


      end if

      if (mod(istep, n_recalc_response) .eq. 0) then
        recalc_response = .true.
      else
        recalc_response = .false.
      end if
      
      if (recalc_response .and. g_exb_error_check_cycle .eq. 0) then
        g_exb_error_check_cycle = 1 ! Calculate without recalculating response 
        g_store = gnew
        phi_store = phinew; apar_store = aparnew; bpar_store = bparnew
        call advance_implicit(istep, .false.)
        g_exb_error_check_cycle = 2 ! Calculate with response
        gnew = g_store 
        phi = phinew; apar = aparnew; bpar = bparnew
        phinew = phi_store; aparnew = apar_store; bparnew = bpar_store
        phi_store = phi; apar_store = apar; bpar_store = bpar
      end if


      if (recalc_response .and. g_exb_error_check_cycle == 2 ) then
        g_store = gnew; phi_temp = phinew; apar_temp =  aparnew;
        bpar_temp =  bparnew
        call reset_init
      end if
      call calculate_kt_grids(g_exb, code_time - wait_time)
      ! set akx in kt_grids equal to akx in kt_grids_single
      !akx = single_akx
      
      !call init_kperp2
      call init_bessel
      call init_fieldeq
      if (recalc_response .and. g_exb_error_check_cycle .eq. 2  ) then
        call init_response_matrix
        gnew = g_store; phinew = phi_temp; aparnew = apar_temp; bparnew =  bpar_temp
        g_exb_error_check_cycle = 3 ! Test error 
      end if

      deallocate(phi_temp)
      deallocate(apar_temp)
      deallocate(bpar_temp)
      return
    end if  ! if (single) 



  end subroutine exb_shear


  subroutine advance_implicit (istep, remove_zonal_flows_switch)
    use run_parameters, only: reset
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use fields_arrays, only: apar_ext !, phi_ext
    use antenna, only: antenna_amplitudes, no_driver
    use dist_fn, only: timeadv !, exb_shear
    use dist_fn_arrays, only: g, gnew, kx_shift, theta0_shift
    !use init_g, only: single_kpar, force_single_kpar
    implicit none
    integer :: diagnostics = 1
    integer, intent (in) :: istep
    logical, intent (in) :: remove_zonal_flows_switch


    !GGH NOTE: apar_ext is initialized in this call
    if(.not.no_driver) call antenna_amplitudes (apar_ext)
       
    call exb_shear  (istep) !EGH

    
    g = gnew
    phi = phinew
    apar = aparnew 
    bpar = bparnew       
    
    call timeadv (phi, apar, bpar, phinew, aparnew, bparnew, istep)
    if(reset) return !Return is resetting

    if(.not.no_driver) aparnew = aparnew + apar_ext 
    
    call getfield (phinew, aparnew, bparnew)
    
    phinew   = phinew  + phi
    aparnew  = aparnew + apar
    bparnew  = bparnew + bpar

    if (remove_zonal_flows_switch) call remove_zonal_flows
    
    call timeadv (phi, apar, bpar, phinew, aparnew, bparnew, istep, diagnostics)
    
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
    use gs2_layouts, only: finish_fields_layouts, finish_jfields_layouts
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

    call finish_fields_layouts
    call finish_jfields_layouts
  end subroutine reset_init

  subroutine init_response_matrix
    use mp, only: barrier
!   use mp, only: proc0
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

    !Either read the reponse
    if(read_response) then
        call read_response_from_file_imp
      !elseif(skip_initialisation) then
       !do i = i_class, 1, -1
          !!Pretty sure this barrier is not needed
          !call barrier
          !!       if (proc0) write(*,*) 'beginning class ',i,' with size ',nidx*N_class(i)
          !!Allocate matrix am. First dimension is basically theta along the entire
          !!connected domain for each field. Second dimension is the local section
          !!of the M_class(i)*N_Class(i)*(2*ntgrid+1)*nfield compound domain.
          !!Clearly this will 
          !allocate (am(nidx*N_class(i), f_lo(i)%llim_proc:f_lo(i)%ulim_alloc))


          !!Do we need to zero all 8 arrays on every loop? This can be more expensive than might think.
          !am = 0.0
          !call init_inverse_matrix (am, i)

          !!Free memory
          !deallocate (am)
       !end do
    else
    !or calculate it

!
! keep storage cost down by doing one class at a time
! Note: could define a superclass (of all classes), a structure containing all am, 
! then do this all at once.  This would be faster, especially for large runs in a 
! sheared domain, and could be triggered by local_field_solve 
! 

!<DD> Comments
!A class refers to a class of connected domain.
!These classes are defined by the extent of the connected domain, there can be 
!many members of each class.
!There are i_class classes in total.
!N_class(ic) is a count of how many 2pi domains there are in members of class ic
!M_class(ic) is how many members of class ic there are.
!Sum N_class(ic)*M_class(ic) for ic=1,i_class is naky*ntheta0
!In comments cell refers to a 2pi domain whilst supercell is the connected domain,
!i.e. we have classes of supercells based on the number of cells they contain.

       do i = i_class, 1, -1
          !Pretty sure this barrier is not needed
          call barrier
          !       if (proc0) write(*,*) 'beginning class ',i,' with size ',nidx*N_class(i)
          !Allocate matrix am. First dimension is basically theta along the entire
          !connected domain for each field. Second dimension is the local section
          !of the M_class(i)*N_Class(i)*(2*ntgrid+1)*nfield compound domain.
          !Clearly this will 
          allocate (am(nidx*N_class(i), f_lo(i)%llim_proc:f_lo(i)%ulim_alloc))


          !Do we need to zero all 8 arrays on every loop? This can be more expensive than might think.
          am = 0.0
          g = 0.0

          phi = 0.0
          apar = 0.0
          bpar = 0.0
          phinew = 0.0
          aparnew = 0.0
          bparnew = 0.0

          !Loop over individual 2pi domains / cells
          do n = 1, N_class(i)
             !Loop over theta grid points in cell
             !This is like a loop over nidx as we also handle all the fields in this loop
             do ig = -ntgrid, ntgrid
                !Are we at a connected boundary point on the lower side (i.e. left hand end of a
                !tube/cell connected to the left)
                endpoint = n > 1
                endpoint = ig == -ntgrid .and. endpoint

                !Start counting fields
                ifield = 0

                !Find response to phi
                if (fphi > epsilon(0.0)) then
                   ifield = ifield + 1
                   if (endpoint) then
                      !Do all members of supercell together
                      do m = 1, M_class(i)
                         ik = f_lo(i)%ik(m,n-1)
                         it = f_lo(i)%it(m,n-1)
                         phinew(ntgrid,it,ik) = 1.0
                      end do
                   endif
                   !Do all members of supercell together
                   do m = 1, M_class(i)
                      ik = f_lo(i)%ik(m,n)
                      it = f_lo(i)%it(m,n)
                      phinew(ig,it,ik) = 1.0
                   end do
                   if (.not. skip_initialisation) call init_response_row (ig, ifield, am, i, n)
                   phinew = 0.0
                end if

                !Find response to apar
                if (fapar > epsilon(0.0)) then
                   ifield = ifield + 1
                   if (endpoint) then
                      !Do all members of supercell together
                      do m = 1, M_class(i)
                         ik = f_lo(i)%ik(m,n-1)
                         it = f_lo(i)%it(m,n-1)
                         aparnew(ntgrid,it,ik) = 1.0
                      end do
                   endif
                   !Do all members of supercell together
                   do m = 1, M_class(i)
                      ik = f_lo(i)%ik(m,n)
                      it = f_lo(i)%it(m,n)
                      aparnew(ig,it,ik) = 1.0
                   end do
                   call init_response_row (ig, ifield, am, i, n)
                   aparnew = 0.0
                end if

                !Find response to bpar
                if (fbpar > epsilon(0.0)) then
                   ifield = ifield + 1
                   if (endpoint) then
                      !Do all members of supercell together
                      do m = 1, M_class(i)
                         ik = f_lo(i)%ik(m,n-1)
                         it = f_lo(i)%it(m,n-1)
                         bparnew(ntgrid,it,ik) = 1.0
                      end do
                   endif
                   !Do all members of supercell together
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

          !Invert the matrix
          call init_inverse_matrix (am, i)

          !Free memory
          deallocate (am)

       end do
    endif 

    if(dump_response) call dump_response_to_file_imp
    call prof_leaving ("init_response_matrix", "fields_implicit")

  end subroutine init_response_matrix

  subroutine init_response_row (ig, ifield, am, ic, n)
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

    !For profiling
    call prof_entering ("init_response_row", "fields_implicit")

    !Always the same size so why bother doing this each time?
    allocate (fieldeq (-ntgrid:ntgrid, ntheta0, naky))
    allocate (fieldeqa(-ntgrid:ntgrid, ntheta0, naky))
    allocate (fieldeqp(-ntgrid:ntgrid, ntheta0, naky))

    !Find response to delta function fields
    !NOTE:Timeadv will loop over all iglo even though only one ik
    !has any amplitude, this is quite a waste. Should ideally do all
    !ik at once
    !NOTE:We currently do each independent supercell of the same length
    !together, this may not be so easy if we do all the ik together but it should
    !be possible.
    call timeadv (phi, apar, bpar, phinew, aparnew, bparnew, 0)
    call getfieldeq (phinew, aparnew, bparnew, fieldeq, fieldeqa, fieldeqp)

    !Loop over 2pi domains / cells
    do nn = 1, N_class(ic)

       !Loop over members of the current class (separate supercells/connected domains)
       do m = 1, M_class(ic)

          !Get corresponding it and ik indices
          it = f_lo(ic)%it(m,nn)
          ik = f_lo(ic)%ik(m,nn)
       
          !Work out which row of the matrix we're looking at
          !corresponds to iindex, i.e. which of the nindex points in the
          !supercell we're looking at.
          irow = ifield + nfield*((ig+ntgrid) + (2*ntgrid+1)*(n-1))
          
          !Convert iindex and m to iflo index
          iflo = idx (f_lo(ic), irow, m)
          
          !If this is part of our local iflo range then store
          !the response data
          if (idx_local(f_lo(ic), iflo)) then
             !Where abouts in the supercell does this 2pi*nfield section start
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

    !Free memory
    deallocate (fieldeq, fieldeqa, fieldeqp)

    !For profiling
    call prof_leaving ("init_response_row", "fields_implicit")
  end subroutine init_response_row

  subroutine init_inverse_matrix (am, ic)
    use file_utils, only: error_unit
    use kt_grids, only: aky, akx
    use theta_grid, only: ntgrid
    use mp, only: broadcast, send, receive, iproc
    use gs2_layouts, only: f_lo, idx, idx_local, proc_id, jf_lo
    use gs2_layouts, only: if_idx, im_idx, in_idx, local_field_solve
    use gs2_layouts, only: ig_idx, ifield_idx, ij_idx, mj, dj
    use prof, only: prof_entering, prof_leaving
    use fields_arrays, only: aminv
    use dist_fn, only: i_class, M_class, N_class
    implicit none
    integer, intent (in) :: ic
    complex, dimension(:,f_lo(ic)%llim_proc:), intent (in out) :: am
    complex, dimension(:,:), allocatable :: a_inv, lhscol, rhsrow, col_row_tmp
    complex, dimension (:), allocatable :: am_tmp
    complex :: fac
    integer :: i, j, k, ik, it, m, n, nn, if, ig, jsc, jf, jg, jc
    integer :: irow, ilo, jlo, dc, iflo, ierr
    logical :: iskip, jskip

    call prof_entering ("init_inverse_matrix", "fields_implicit")
    
    allocate (lhscol (nidx*N_class(ic),M_class(ic)))
    allocate (rhsrow (nidx*N_class(ic),M_class(ic)))
   
    !This is the length of a supercell
    j = nidx*N_class(ic)

    !Create storage space
    allocate (a_inv(j,f_lo(ic)%llim_proc:f_lo(ic)%ulim_alloc))
    a_inv = 0.0
    
    if (.not. skip_initialisation) then
      !Set (ifield*ig,ilo) "diagonal" to 1?
      do ilo = f_lo(ic)%llim_proc, f_lo(ic)%ulim_proc
         a_inv(if_idx(f_lo(ic),ilo),ilo) = 1.0
      end do

      ! Gauss-Jordan elimination, leaving out internal points at multiples of ntgrid 
      ! for each supercell
      !Loop over parallel gridpoints in supercell
      do i = 1, nidx*N_class(ic)
         !iskip is true iff the theta grid point(ig) corresponding to i
         !is at the upper end of a 2pi domain/cell and is not the rightmost gridpoint
         iskip = N_class(ic) > 1 !Are the multiple cells => are there connections/boundaries
         iskip = i <= nidx*N_class(ic) - nfield .and. iskip !Are we not near the upper boundary of the supercell
         iskip = mod((i+nfield-1)/nfield, 2*ntgrid+1) == 0 .and. iskip !Are we at a theta grid point corresponding to the rightmost point of a 2pi domain
         iskip = i > nfield .and. iskip !Are we not at the lower boundary of the supercell
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
            allocate(col_row_tmp(nidx*N_class(ic),2)) ; col_row_tmp = 0.
            !Loop over classes (supercell lengths)
            do m = 1, M_class(ic)
               !Convert to f_lo index
               ilo = idx(f_lo(ic),i,m)
               !Is ilo on this proc?
               if (idx_local(f_lo(ic),ilo)) then
                  !If so store column/row
                  !lhscol(:,m) = am(:,ilo)
                  !rhsrow(:,m) = a_inv(:,ilo)
                  col_row_tmp(:,1) = am(:,ilo)
                  col_row_tmp(:,2) = a_inv(:,ilo)
               end if
               !Here we send lhscol and rhscol sections to all procs
               !from the one on which it is currently known
               !Can't do this outside m loop as proc_id depends on m
               !These broadcasts can be relatively expensive so local_field_solve
               !may be preferable
               !call broadcast (lhscol(:,m), proc_id(f_lo(ic),ilo))
               !call broadcast (rhsrow(:,m), proc_id(f_lo(ic),ilo))
               call broadcast (col_row_tmp, proc_id(f_lo(ic),ilo))
               lhscol(:,m) = col_row_tmp(:,1)
               rhsrow(:,m) = col_row_tmp(:,2)
            end do
            !All procs will have the same lhscol and rhsrow after this loop+broadcast
            deallocate(col_row_tmp)
         end if

         !Loop over field compound dimension
         do jlo = f_lo(ic)%llim_proc, f_lo(ic)%ulim_proc
            !jskip is true similarly to iskip
            jskip = N_class(ic) > 1 !Are there any connections?
            jskip = ig_idx(f_lo(ic), jlo) == ntgrid .and. jskip !Are we at a theta grid point corresponding to the upper boundary?
            !Get 2pi domain/cell number out of total for this supercell
            n = in_idx(f_lo(ic),jlo)
            jskip = n < N_class(ic) .and. jskip !Are we not in the last cell (i.e. not at the rightmost grid point/upper end of supercell)?
            if (jskip) cycle  !Skip this point if appropriate

            !Now get m (class number)
            m = im_idx(f_lo(ic),jlo)

            !Convert class number and cell number to ik and it
            ik = f_lo(ic)%ik(m,n)
            it = f_lo(ic)%it(m,n)
            
            !Work out what the compound theta*field index is.
            irow = if_idx(f_lo(ic),jlo)

            !If ky or kx are not 0 (i.e. skip zonal 0,0 mode) then workout the array
            if (aky(ik) /= 0.0 .or. akx(it) /= 0.0) then
               !Get factor
               fac = am(i,jlo)/lhscol(i,m)

               !Store array element
               am(i,jlo) = fac

               !Store other elements
               am(:i-1,jlo) = am(:i-1,jlo) - lhscol(:i-1,m)*fac
               am(i+1:,jlo) = am(i+1:,jlo) - lhscol(i+1:,m)*fac
               !WOULD the above three commands be better written as
               !am(:,jlo)=am(:,jlo)-lhscol(:,m)*fac
               !am(i,jlo)=fac

               !Fill in a_inv
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

      !Free memory
      deallocate (lhscol, rhsrow)

  ! fill in skipped points for each field and supercell:
  ! Do not include internal ntgrid points in sum over supercell

      do i = 1, nidx*N_class(ic)
         !iskip is true iff the theta grid point(ig) corresponding to i
         !is at the upper end of a 2pi domain/cell and is not the rightmost gridpoint
         iskip = N_class(ic) > 1 !Are the multiple cells => are there connections/boundaries
         iskip = i <= nidx*N_class(ic) - nfield .and. iskip  !Are we not near the upper boundary of the supercell
         iskip = mod((i+nfield-1)/nfield, 2*ntgrid+1) == 0 .and. iskip !Are we at a theta grid point corresponding to the rightmost point of a 2pi domain
         iskip = i > nfield .and. iskip !Are we not at the lower boundary of the supercell
         !Zero out skipped points
         if (iskip) then
            a_inv(i,:) = 0
            cycle !Seems unnexessary
         end if
      end do
  ! Make response at internal ntgrid points identical to response
  ! at internal -ntgrid points:
      do jlo = f_lo(ic)%llim_world, f_lo(ic)%ulim_world
         !jskip is true similarly to iskip
         jskip = N_class(ic) > 1 !Are there any connections?
         jskip = ig_idx(f_lo(ic), jlo) == ntgrid .and. jskip  !Are we at a theta grid point corresponding to the upper boundary?
         jskip = in_idx(f_lo(ic), jlo) < N_class(ic) .and. jskip  !Are we not in the last cell (i.e. not at the rightmost grid point/upper end of supercell)?
         !If we previously skipped this point then we want to fill it in from the matched/connected point
         if (jskip) then
            !What is the index of the matched point?
            ilo = jlo + nfield
            !If we have ilo on this proc send it to...
            if (idx_local(f_lo(ic), ilo)) then
               !jlo on this proc
               if (idx_local(f_lo(ic), jlo)) then
                  a_inv(:,jlo) = a_inv(:,ilo)
               !jlo on proc which has jlo
               else
                  call send(a_inv(:,ilo), proc_id(f_lo(ic), jlo))
               endif
            else
               !If this proc has jlo then get ready to receive
               if (idx_local(f_lo(ic), jlo)) then
                  call receive(a_inv(:,jlo), proc_id(f_lo(ic), ilo))
               end if
            end if
         end if
      end do
      !The send receives in the above loop should be able to function in a
      !non-blocking manner fairly easily, but probably don't cost that much
      !Would require WAITALL before doing am=a_inv line below

      !Update am
      am = a_inv
    end if ! .not. skip_initialisation

    !Free memory
    deallocate (a_inv)

! Re-sort this class of aminv for runtime application.  

    !Now allocate array to store matrices for each class
    if (.not.allocated(aminv)) allocate (aminv(i_class))

! only need this large array for particular values of jlo.
! To save space, count how many this is and only allocate
! required space:

    !Initialise counter
    dc = 0
! check all members of this class
    do ilo = f_lo(ic)%llim_world, f_lo(ic)%ulim_world

! find supercell coordinates
       !i.e. what is my class of supercell and which cell am I looking at
       m = im_idx(f_lo(ic), ilo)
       n = in_idx(f_lo(ic), ilo)

! find standard coordinates
       !Get theta, field, kx and ky indexes for current point
       ig = ig_idx(f_lo(ic), ilo)
       if = ifield_idx(f_lo(ic), ilo)
       ik = f_lo(ic)%ik(m,n)
       it = f_lo(ic)%it(m,n)

! translate to fast field coordinates
       jlo = ij_idx(jf_lo, ig, if, ik, it)
          
! Locate this jlo, count it, and save address
       !Is this point on this proc, if so increment counter
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
    !Loop over "fast field" index
    do jlo = jf_lo%llim_proc, jf_lo%ulim_proc
          
       !Allocate store in this class, on this proc to store the jlo points
       if (.not.associated(aminv(ic)%dcell)) then
          allocate (aminv(ic)%dcell(dc))
       else
          !Just check the array is the correct size
          j = size(aminv(ic)%dcell)
          if (j /= dc) then
             ierr = error_unit()
             write(ierr,*) 'Error (1) in init_inverse_matrix: ',&
                  iproc,':',jlo,':',dc,':',j
          endif
       endif
       
       !Get the current "dcell" adress
       k = dj(ic,jlo)

       !No dcell should be 0 but this is a guard
       if (k > 0) then
          !How long is the supercell for this class?
          jc = nidx*N_class(ic)

          !Allocate storage for the supercell if required
          if (.not.associated(aminv(ic)%dcell(k)%supercell)) then
             allocate (aminv(ic)%dcell(k)%supercell(jc))
          else
             !Just check the array is the correct size
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

    !Allocate temporary supercell storage
    allocate (am_tmp(nidx*N_class(ic)))

    !Loop over all grid points
    do ilo = f_lo(ic)%llim_world, f_lo(ic)%ulim_world

       !Get supercell type (class) and cell index
       m = im_idx(f_lo(ic), ilo)
       n = in_idx(f_lo(ic), ilo)
       
       !Convert to theta,field,kx and ky indexes
       ig = ig_idx(f_lo(ic), ilo)
       if = ifield_idx(f_lo(ic), ilo)
       ik = f_lo(ic)%ik(m,n)
       it = f_lo(ic)%it(m,n)
       
       !Get fast field index
       iflo = ij_idx(jf_lo, ig, if, ik, it)
 
       !If this ilo is local then...
       if (idx_local(f_lo(ic),ilo)) then
          ! send the am data to...
          if (idx_local(jf_lo,iflo)) then
             !the local proc
             am_tmp = am(:,ilo)
          else
             !the remote proc
             call send(am(:,ilo), proc_id(jf_lo,iflo))
          endif
       else
          !Get ready to receive the data
          if (idx_local(jf_lo,iflo)) then
             call receive(am_tmp, proc_id(f_lo(ic),ilo))
          end if
       end if

       !If the fast field index is on this processor
       if (idx_local(jf_lo, iflo)) then
          !Get "dcell" adress
          dc = dj(ic,iflo)

          !Loop over supercell size
          do jlo = 0, nidx*N_class(ic)-1
             !Convert to cell/2pi domain index
             nn = in_idx(f_lo(ic), jlo)
             
             !Get theta grid point
             jg = ig_idx(f_lo(ic), jlo)
             !Get field index
             jf = ifield_idx(f_lo(ic), jlo)
             
             !Convert index
             jsc = ij_idx(f_lo(ic), jg, jf, nn) + 1

             !Store inverse matrix data in appropriate supercell position
             aminv(ic)%dcell(dc)%supercell(jsc) = am_tmp(jlo+1)
             
          end do
       end if
    end do

    !Free memory
    deallocate (am_tmp)

    !For profiling
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
    use mp, only: mp_abort
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
       call mp_abort('Error in kt2ki')
    end if
! 
! Get position in this supercell, counting from the left
!
    n = 1 + l_links(ik, it)

  end subroutine kt2ki

  !>A routine to dump the current response matrix to file
  subroutine dump_response_to_file_imp(suffix)
    use file_utils, only: run_name
    use fields_arrays, only: aminv
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use dist_fn, only: i_class, N_class, M_class, get_leftmost_it, itright
    use gs2_layouts, only: f_lo, jf_lo, ij_idx, idx_local, idx, proc_id,idx_local,dj
    use mp, only: proc0, send, receive
    use gs2_save, only: gs2_save_response
    implicit none
    character(len=*), optional, intent(in) :: suffix !If passed then use as part of file suffix
    character(len=64) :: suffix_local, suffix_default='.response'
    character(len=256) :: file_name
    complex, dimension(:,:), allocatable :: tmp_arr, tmp_arr_full
    complex, dimension(:), allocatable :: tmp_vec_full, tmp_vec
    integer :: ic, im, ik, it, itmin, supercell_length, supercell_length_bound, in, ifld, ig, is_tmp
    integer :: jflo, dc, nn, in_tmp, icount, it_tmp, nl, nr, ifld_tmp, ext_dom_length, ig_tmp, cur_idx
    integer, dimension(:,:), allocatable :: it_to_is, leftmost_it
    integer, dimension(:), allocatable :: tmp_ints
    logical :: is_local
    !Set file suffix
    suffix_local=suffix_default
    if(present(suffix)) suffix_local=suffix

    !Make a lookup array to convert itmin (the leftmost it in a connected domain)
    !to the supercell index "is" used in the local fields. This will be used to 
    !ensure equivalent files can be given the same name.
    allocate(it_to_is(ntheta0,naky),leftmost_it(ntheta0,naky),tmp_ints(ntheta0))
    it_to_is=0
    !//Note the following code is mostly borrowed from fm_init in the local fields
    
    !First find all the leftmost it
    do ik=1,naky
       do it=1,ntheta0
          leftmost_it(it,ik)=get_leftmost_it(it,ik)
       enddo
    enddo

    !Now find supercell ids for each ky at a time
    do ik=1,naky
       tmp_ints=leftmost_it(:,ik)
       it_tmp=0
       is_tmp=0
       do while(sum(tmp_ints).ne.-1*ntheta0)
          it_tmp=it_tmp+1
          cur_idx=tmp_ints(it_tmp)

          !If we've seen this domain skip
          if(cur_idx.eq.-1)cycle

          !Increment counter
          is_tmp=is_tmp+1

          !Here we store the value
          it_to_is(it_tmp,ik)=is_tmp

          !Now we set all other connected locations to -1
          !and store the appropriate is value
          do it=1,ntheta0
             if(tmp_ints(it).eq.cur_idx) then
                tmp_ints(it)=-1
                it_to_is(it,ik)=is_tmp
             endif
          enddo
       enddo
    enddo

    !Cleanup
    deallocate(tmp_ints)

    !/End of borrowed code

    !Notation recap:
    ! A class refers to all independent domains with the same length
    ! i_class is how many classes we have
    ! N_class(i_class) is how many 2Pi domains are in each member of i_class
    ! M_class(i_class) is how many independent domains are in i_class

    allocate(tmp_vec(nfield*(2*ntgrid+1)))
    allocate(tmp_arr(1+(2*ntgrid),nfield))

    !Loop over classes (supercell length)
    do ic=1,i_class
       !Work out how long the supercell is
       supercell_length=1+(2*ntgrid)*nfield*N_class(ic) !Without boundary points
       supercell_length_bound=(1+2*ntgrid)*nfield*N_class(ic) !With boundary points
       !Extended domain length
       ext_dom_length=1+(2*ntgrid)*N_class(ic)

       !Make storage
       allocate(tmp_arr_full(supercell_length,supercell_length))
       allocate(tmp_vec_full(supercell_length))

       !Now loop over all members of this class
       do im=1,M_class(ic)
          !Now we are thinking about a single supercell
          !we can get certain properties before looping
          !over the individual elements
          
          !Get the ik index
          ik=f_lo(ic)%ik(im,1)

          !Get the leftmost it index (named itmin to match local field routines)
          !This is currently used to identify the supercell like "is" is used in
          !the local field routines. It would be nice to also use "is" here (or
          !"itmin" there).
          itmin=leftmost_it(f_lo(ic)%it(im,1),ik)
          
          !Now we have the basic properties we want to loop over the elements
          !First initialise "it"
          it=itmin

          !Initialise counter
          icount=1

          !Loop over the different it (2Pi domains)
          do in=1,N_class(ic)
             !Loop over the fields
             do ifld=1,nfield
                !Loop over theta
                do ig=-ntgrid,ntgrid
                   !Skip the duplicate boundary points
                   if((ig.eq.ntgrid).and.(in.ne.N_class(ic))) cycle

                   !Convert to jf_lo index
                   jflo=ij_idx(jf_lo,ig,ifld,ik,it)

                   !See if it's local
                   is_local=idx_local(jf_lo,jflo)

                   !If it's not local then we have nothing to do
                   !unless we're the proc who writes (proc0).
                   if(.not.(is_local.or.proc0)) cycle

                   !Now pack tmp_vec and do communications if needed
                   if(is_local)then
                      !Get dcell index
                      dc=dj(ic,jflo)

                      !Now we pack the tmp_vec in the correct order
                      !whilst ignoring the repeated boundary points
                      !We need to pick the value of "n" in the right order
                      it_tmp=itmin
                      do in_tmp=1,N_class(ic)
                         !Pick the correct n
                         do nn=1,N_class(ic)
                            if(f_lo(ic)%it(im,nn).eq.it_tmp) exit
                         enddo

                         !Now we can get supercell range (including boundaries)
                         nl=1+nidx*(nn-1)
                         nr=nl+nidx-1
                         
                         !Extract section
                         tmp_vec=aminv(ic)%dcell(dc)%supercell(nl:nr)

                         !All that remains now is to ignore the boundary points
                         !To do this we just split on the field so we can ignore
                         !boundary if we want
                         do ifld_tmp=1,nfield
                            nl=1+(ifld_tmp-1)*(2*ntgrid+1)
                            nr=nl+2*ntgrid
                            tmp_arr(:,ifld_tmp)=tmp_vec(nl:nr)
                         enddo

                         !Now we need to work out where to put things in tmp_vec_full
                         !In doing this we try to match the local fields data layout
                         !to aid comparisons
                         do ifld_tmp=1,nfield
                            do ig_tmp=1,2*ntgrid+1
                               !Skip boundary points
                               if((ig_tmp.eq.(2*ntgrid+1)).and.(in_tmp.ne.N_class(ic))) cycle

                               !Get index
                               cur_idx=ig_tmp+(2*ntgrid)*(in_tmp-1)+(ifld_tmp-1)*ext_dom_length

                               !Store data
                               tmp_vec_full(cur_idx)=tmp_arr(ig_tmp,ifld_tmp)
                            enddo
                         enddo

                         !Increment it
                         it_tmp=itright(ik,it_tmp)
                      enddo

                      !No comms needed if on proc0
                      if(.not.proc0) call send(tmp_vec_full,0)
                   else
                      !Only proc0 should get here but test anyway
                      if(proc0) call receive(tmp_vec_full,proc_id(jf_lo,jflo))
                   endif

                   !Now we need to store in the full array
                   !May need to check index order matches local case.
                   if(proc0) then
                      tmp_arr_full(:,icount)=tmp_vec_full
                   endif

                   !Increment counter
                   icount=icount+1
                enddo
             enddo

             !Increment it
             it=itright(ik,it)
          enddo

          !Now make file name
          if(proc0)then
             write(file_name,'(A,"_ik_",I0,"_is_",I0,A)') trim(run_name),ik,it_to_is(itmin,ik),trim(suffix_local)
             call gs2_save_response(tmp_arr_full,file_name)
          endif
       end do
          
       deallocate(tmp_arr_full,tmp_vec_full)
    end do

    !Tidy
    deallocate(tmp_vec,tmp_arr,leftmost_it,it_to_is)

  end subroutine dump_response_to_file_imp

  !>A routine to read the response matrix from file and populate the implicit
  !response storage, note we also allocate the response storage objects
  subroutine read_response_from_file_imp(suffix)
    use file_utils, only: run_name
    use fields_arrays, only: aminv
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use dist_fn, only: i_class, N_class, M_class, get_leftmost_it, itright
    use gs2_layouts, only: f_lo, jf_lo, ij_idx, idx_local, idx, proc_id,idx_local,dj
    use mp, only: proc0, send, receive
    use gs2_save, only: gs2_restore_response
    implicit none
    character(len=*), optional, intent(in) :: suffix !If passed then use as part of file suffix
    character(len=64) :: suffix_local, suffix_default='.response'
    character(len=256) :: file_name
    complex, dimension(:,:), allocatable :: tmp_arr, tmp_arr_full
    complex, dimension(:), allocatable :: tmp_vec_full, tmp_vec
    integer :: ic, im, ik, it, itmin, supercell_length, supercell_length_bound, in, ifld, ig, is_tmp
    integer :: jflo, dc, nn, in_tmp, icount, it_tmp, nl, nr, ifld_tmp, ext_dom_length, ig_tmp, cur_idx
    integer :: jflo_dup, dc_dup
    integer, dimension(:,:), allocatable :: it_to_is, leftmost_it
    integer, dimension(:), allocatable :: tmp_ints
    logical :: is_local, is_local_dup
    !Set file suffix
    suffix_local=suffix_default
    if(present(suffix)) suffix_local=suffix

    !First allocate the matrix storage
    call alloc_response_objects

    !Make a lookup array to convert itmin (the leftmost it in a connected domain)
    !to the supercell index "is" used in the local fields. This will be used to 
    !ensure equivalent files can be given the same name.
    allocate(it_to_is(ntheta0,naky),leftmost_it(ntheta0,naky),tmp_ints(ntheta0))
    it_to_is=0
    !//Note the following code is mostly borrowed from fm_init in the local fields
    
    !First find all the leftmost it
    do ik=1,naky
       do it=1,ntheta0
          leftmost_it(it,ik)=get_leftmost_it(it,ik)
       enddo
    enddo

    !Now find supercell ids for each ky at a time
    do ik=1,naky
       tmp_ints=leftmost_it(:,ik)
       it_tmp=0
       is_tmp=0
       do while(sum(tmp_ints).ne.-1*ntheta0)
          it_tmp=it_tmp+1
          cur_idx=tmp_ints(it_tmp)

          !If we've seen this domain skip
          if(cur_idx.eq.-1)cycle

          !Increment counter
          is_tmp=is_tmp+1

          !Here we store the value
          it_to_is(it_tmp,ik)=is_tmp

          !Now we set all other connected locations to -1
          !and store the appropriate is value
          do it=1,ntheta0
             if(tmp_ints(it).eq.cur_idx) then
                tmp_ints(it)=-1
                it_to_is(it,ik)=is_tmp
             endif
          enddo
       enddo
    enddo

    !Cleanup
    deallocate(tmp_ints)

    !/End of borrowed code

    !Notation recap:
    ! A class refers to all independent domains with the same length
    ! i_class is how many classes we have
    ! N_class(i_class) is how many 2Pi domains are in each member of i_class
    ! M_class(i_class) is how many independent domains are in i_class

    allocate(tmp_vec(nfield*(2*ntgrid+1)))
    allocate(tmp_arr(1+(2*ntgrid),nfield))

    !Loop over classes (supercell length)
    do ic=1,i_class
       !Work out how long the supercell is
       supercell_length=1+(2*ntgrid)*nfield*N_class(ic) !Without boundary points
       supercell_length_bound=(1+2*ntgrid)*nfield*N_class(ic) !With boundary points
       !Extended domain length
       ext_dom_length=1+(2*ntgrid)*N_class(ic)

       !Make storage
       allocate(tmp_arr_full(supercell_length,supercell_length))
       allocate(tmp_vec_full(supercell_length))

       !Now loop over all members of this class
       do im=1,M_class(ic)
          tmp_arr_full=0.
          tmp_vec_full=0.

          !Now we are thinking about a single supercell
          !we can get certain properties before looping
          !over the individual elements
          
          !Get the ik index
          ik=f_lo(ic)%ik(im,1)

          !Get the leftmost it index (named itmin to match local field routines)
          !This is currently used to identify the supercell like "is" is used in
          !the local field routines. It would be nice to also use "is" here (or
          !"itmin" there).
          itmin=leftmost_it(f_lo(ic)%it(im,1),ik)
          
          !Now we have the basic properties we want to loop over the elements
          !First initialise "it"
          it=itmin

          !Now make file name
          if(proc0)then
             write(file_name,'(A,"_ik_",I0,"_is_",I0,A)') trim(run_name),ik,it_to_is(itmin,ik),trim(suffix_local)
             call gs2_restore_response(tmp_arr_full,file_name)
          endif

          !Initialise counter
          icount=1

          !Loop over the different it (2Pi domains)
          do in=1,N_class(ic)
             !Loop over the fields
             do ifld=1,nfield
                !Loop over theta
                do ig=-ntgrid,ntgrid
                   !Skip the duplicate boundary points -- This is no good here. !<DD>
!                   if((ig.eq.ntgrid).and.(in.ne.N_class(ic))) cycle

                   !Convert to jf_lo index
                   jflo=ij_idx(jf_lo,ig,ifld,ik,it)

                   !See if it's local
                   is_local=idx_local(jf_lo,jflo)

                   !If it's not local then we have nothing to do
                   !unless we're the proc who writes (proc0).
                   if(.not.(is_local.or.proc0)) cycle

                   !Get row
                   if(proc0)then
                      tmp_vec_full=tmp_arr_full(:,icount)
                      
                      !Increment counter
                      if(.not.(ig.eq.ntgrid.and.in.ne.N_Class(ic))) icount=icount+1
                   endif

                   !Now unpack tmp_vec_full and do communications if needed
                   if(is_local)then
                      !No comms needed if local
                      if(.not.proc0) call receive(tmp_vec_full,0)

                      !Get dcell index
                      dc=dj(ic,jflo)

                      !Now we pack the tmp_vec in the correct order
                      !We must fill in the boundary points
                      !We need to pick the value of "n" in the right order
                      it_tmp=itmin
                      do in_tmp=1,N_class(ic)
                         tmp_arr=0
                         tmp_vec=0

                         !Now we need to work out where to put things in tmp_vec_full
                         !In doing this we try to match the local fields data layout
                         !to aid comparisons
                         do ifld_tmp=1,nfield
                            do ig_tmp=1,2*ntgrid+1
                               !Skip boundary points
                               if((ig_tmp.eq.(2*ntgrid+1)).and.(in_tmp.ne.N_class(ic))) cycle

                               !Get index
                               cur_idx=ig_tmp+(2*ntgrid)*(in_tmp-1)+(ifld_tmp-1)*ext_dom_length

                               !Store data
                               tmp_arr(ig_tmp,ifld_tmp)=tmp_vec_full(cur_idx)
                            enddo
                         enddo

                         !<DD>It may be anticipated that we need to fix the boundary points
                         !here but we don't actually need to do anything.
                         !Because we sum over the entire supercell in getfield we only want
                         !the repeated boundary point to be included once.
                         !We still need to calculate the field at the repeated point but the
                         !fix for that is handled at the bottom of the routine
                         !In other words we don't need something of the form:
                         ! !Fix boundary points
                         ! if(in_tmp.ne.N_class(ic))then
                         !    do ifld_tmp=1,nfield
                         !       cur_idx=1+(2*ntgrid)*(in_tmp)+(ifld_tmp-1)*ext_dom_length
                         !       tmp_arr(2*ntgrid+1,ifld_tmp)=tmp_vec_full(cur_idx)
                         !    enddo
                         ! endif

                         !Store in correct order
                         do ifld_tmp=1,nfield
                            nl=1+(ifld_tmp-1)*(2*ntgrid+1)
                            nr=nl+2*ntgrid
                            tmp_vec(nl:nr)=tmp_arr(:,ifld_tmp)
                         enddo

                         !Pick the correct n
                         do nn=1,N_class(ic)
                            if(f_lo(ic)%it(im,nn).eq.it_tmp) exit
                         enddo

                         !Now we can get supercell range (including boundaries)
                         nl=1+nidx*(nn-1)
                         nr=nl+nidx-1

                         !Store section
                         aminv(ic)%dcell(dc)%supercell(nl:nr)=tmp_vec

                         !Increment it
                         it_tmp=itright(ik,it_tmp)
                      enddo
                   else
                      !Only proc0 should get here but test anyway
                      if(proc0) call send(tmp_vec_full,proc_id(jf_lo,jflo))
                   endif
                enddo
             enddo

             !Increment it
             it=itright(ik,it)
          enddo

          !Now we need to fill in the repeated boundary points

          !If there are no boundary points then advance
          if(N_class(ic).eq.1) cycle
          it=itmin
          do in=1,N_class(ic)-1
             do ifld=1,nfield
                !First get the index of the point we want to fill
                jflo=ij_idx(jf_lo,ntgrid,ifld,ik,it)

                !Now we get the index of the point which has this data
                jflo_dup=ij_idx(jf_lo,-ntgrid,ifld,ik,itright(ik,it))

                !Now get locality
                is_local=idx_local(jf_lo,jflo)
                is_local_dup=idx_local(jf_lo,jflo_dup)

                !Get dcell values
                if(is_local) dc=dj(ic,jflo)
                if(is_local_dup) dc_dup=dj(ic,jflo_dup)

                !Now copy/communicate
                if(is_local)then
                   if(is_local_dup)then
                      aminv(ic)%dcell(dc)%supercell=aminv(ic)%dcell(dc_dup)%supercell
                   else
                      call receive(aminv(ic)%dcell(dc)%supercell,proc_id(jf_lo,jflo_dup))
                   endif
                elseif(is_local_dup)then
                   call send(aminv(ic)%dcell(dc_dup)%supercell,proc_id(jf_lo,jflo))
                endif
             enddo

             !Increment it
             it=itright(ik,it)
          enddo
       end do
       
       !Free
       deallocate(tmp_arr_full,tmp_vec_full)
    end do

    !Tidy
    deallocate(tmp_vec,tmp_arr,leftmost_it,it_to_is)
  end subroutine read_response_from_file_imp

  !>A subroutine to allocate the response matrix storage objects
  subroutine alloc_response_objects
    use dist_fn, only: i_class, N_class
    use fields_arrays, only: aminv
    use gs2_layouts, only: jf_lo, f_lo, im_idx, in_idx, ig_idx, ifield_idx, ij_idx,dj,mj, idx_local
    use theta_grid, only: ntgrid
    implicit none
    integer :: ic, idc, sc_len, ilo, dc, im, in, ig, ifld, ik, it, jlo

    !Top level, one object for each class (length of supercell)
    if(.not.allocated(aminv)) allocate(aminv(i_class))

    !Loop over each class
    do ic=1,i_class
       !Get the supercell length
       sc_len=(2*ntgrid+1)*nfield*N_class(ic)

       !Count how many dcell we have locally and fill related data
       dc=0
       do ilo=f_lo(ic)%llim_world,f_lo(ic)%ulim_world
          !i.e. what is my class of supercell and which cell am I looking at
          im = im_idx(f_lo(ic), ilo)
          in = in_idx(f_lo(ic), ilo)

          ! find standard coordinates
          !Get theta, field, kx and ky indexes for current point
          ig = ig_idx(f_lo(ic), ilo)
          ifld = ifield_idx(f_lo(ic), ilo)
          ik = f_lo(ic)%ik(im,in)
          it = f_lo(ic)%it(im,in)
          
          ! translate to fast field coordinates
          jlo = ij_idx(jf_lo, ig, ifld, ik, it)
          
          ! Locate this jlo, count it, and save address
          !Is this point on this proc, if so increment counter
          if (idx_local(jf_lo,jlo)) then
             ! count it
             dc = dc + 1
             ! save dcell address
             dj(ic,jlo) = dc
             ! save supercell address
             mj(jlo) = im
          endif
       enddo

       !Next level, one object for each point in the class
       if(.not.associated(aminv(ic)%dcell))then
          allocate(aminv(ic)%dcell(dc))
       endif

       !Now loop over each point and allocate storage for the response data
       do idc=1,dc
          !Bottom level, this is actually where data is stored
          if(.not.associated(aminv(ic)%dcell(idc)%supercell)) then
             allocate(aminv(ic)%dcell(idc)%supercell(sc_len))
          endif
       enddo
    enddo

  end subroutine alloc_response_objects
end module fields_implicit

