!> A small module used to calculate ideal ballooning stability.
!! Based on the original ball program in geo/ball.f90 of GS2
!! but tweaked to integrate closer to full GS2 runs etc.
module ballstab
  implicit none

  private

  !Routines
  public :: init_ballstab, finish_ballstab
  public :: run_stability_check, is_unstable
  public :: write_stability_ascii_1d, write_stability_ascii_2d

  !Vars
  public :: stability

  logical :: make_salpha
  integer :: n_shat, n_beta
  real :: shat_min, shat_max
  real :: beta_mul, beta_div
  real, dimension(:), allocatable :: shat_arr, beta_arr, dbdrho_arr
  integer, dimension(:,:), allocatable :: stability
  real :: diff

  logical :: initialised = .false. 
  character(len=14) :: namelist_name="ballstab_knobs"

contains
  
  !> Initialise grids, other modules etc.
  subroutine init_ballstab
    use theta_grid, only: init_theta_grid
    use mp, only: proc0
    implicit none

    !Other modules
    call init_theta_grid

    !This module -- Note only proc0 will be allowed to do any work
    if(proc0) call real_init_ballstab
  end subroutine init_ballstab

  !> Finalise this module
  subroutine finish_ballstab
    implicit none
    !>Deallocate
    call dealloc_arrays

    initialised = .false.
  end subroutine finish_ballstab

  !> Initialise this module
  subroutine real_init_ballstab
    implicit none

    if(initialised) return
    initialised = .true.

    !Get settings
    call read_parameters

    !Verify parameters
    call verify_parameters

    !Create arrays
    call setup_arrays

  end subroutine real_init_ballstab

  !> Read parameters
  subroutine read_parameters
    use file_utils, only: input_unit, input_unit_exist
    use theta_grid, only: shat
    implicit none
    namelist/ballstab_knobs/make_salpha, n_shat, shat_min, shat_max, &
         n_beta, beta_mul, beta_div, diff
    integer :: in_file
    logical :: exist

    !Set defaults
    make_salpha = .false.
    n_shat = 1
    shat_min = shat
    shat_max = shat
    n_beta = 1
    beta_mul = 1.
    beta_div = 1.
    diff = 0.

    !Check if namelist is present in input file
    in_file = input_unit_exist (namelist_name, exist)

    !If it is then read in parameters
    if(exist) read(unit = in_file, nml = ballstab_knobs)
  end subroutine read_parameters

  !> Check parameters are consistent etc.
  subroutine verify_parameters
    use geometry, only: bishop, iflux
    use file_utils, only: error_unit
    implicit none
    integer :: ierr

    ierr = error_unit()

    !Disable scan in s_hat when incompatible bishop used
    if(n_shat.gt.1)then
       select case(bishop)
       case(2,3,4,5,8,9)
          !These are valid values
       case default
          if((bishop.ne.0).and.(iflux.ne.1)) then
             !This is also ok
          else
             !This is not ok
             n_shat=1
             write(ierr,'("ERROR : Cannot scan in shat for current combination of bishop and iflux")')
          endif
       end select
    endif

    !Disable scan in beta_prime when incompatible bishop used
    if(n_beta.gt.1)then
       select case(bishop)
       case(2,3,4,5,6,7,8,9)
          !These are valid values
       case default
          !This is not ok
          n_beta=1
          write(ierr,'("ERROR : Cannot scan in beta_prime for current combination of bishop and iflux")')
       end select
    endif

  end subroutine verify_parameters

  !> Allocate+populate arrays
  subroutine setup_arrays
    implicit none
    integer :: is, ib
    real :: ds, db, b_min, b_max, bp
    
    allocate(shat_arr(n_shat))
    allocate(beta_arr(n_beta))
    allocate(dbdrho_arr(n_beta))
    dbdrho_arr = 0. !This will be populated later

    allocate(stability(n_shat,n_beta))
    stability = 0

    ds = 0
    if(n_shat.gt.1) ds = (shat_max-shat_min)/(n_shat-1)
    do is=1,n_shat
       shat_arr(is)=shat_min+ds*(is-1)
    enddo

    db = 0
    call get_beta_prime(bp)
    b_max = bp*beta_mul
    b_min = bp/beta_div
    if(n_beta.gt.1) then
       db = (b_max-b_min)/(n_beta-1)
    end if
    do ib=1,n_beta
       beta_arr(ib)=b_min+db*(ib-1)
    enddo
  end subroutine setup_arrays

  !> Deallocate arrays
  subroutine dealloc_arrays
    implicit none
    if(allocated(shat_arr)) deallocate(shat_arr)
    if(allocated(beta_arr)) deallocate(beta_arr)
    if(allocated(dbdrho_arr)) deallocate(dbdrho_arr)
    if(allocated(stability)) deallocate(stability)
  end subroutine dealloc_arrays

  !> Write a namelist based on current values
  subroutine wnml_ballstab(unit)
    implicit none
    integer, intent(in) :: unit

    !Header
    write(unit,'()')
    write(unit,'("&",A)') trim(adjustl(namelist_name))

    !Values
    write(unit,'(" make_salpha = ",L1)') make_salpha
    write(unit,'(" n_shat      = ",I0)') n_shat
    if(n_shat.gt.1) then
       write(unit,'(" shat_min    = ",F12.5)') shat_min
       write(unit,'(" shat_max    = ",F12.5)') shat_max
    endif
    write(unit,'(" n_beta      = ",I0)') n_beta

    !Footer
    write(unit,'("/")')
  end subroutine wnml_ballstab

  !> Check value consistency etc.
  subroutine check_ballstab(report_unit)
    implicit none
    integer, intent(in) :: report_unit
  end subroutine check_ballstab

  !> Gets the current value of shat
  subroutine get_shat(shat_out)
    use geometry, only: shat
    implicit none
    real, intent(out) :: shat_out
    shat_out = shat
  end subroutine get_shat

  !> Sets the value of shat
  subroutine set_shat(shat_in)
    use geometry, only: s_hat_input, shat
    implicit none
    real, intent(in) :: shat_in
    !Note for bishop = 6,7 or not in 2-9 and iflux==1 cannot change shat
    !Should flag this somewhere (during init)
    s_hat_input = shat_in
    shat = shat_in
  end subroutine set_shat

  !> Gets the current value of beta_prime (or equivalent var)
  !! Really just gets the variable that we can use to control 
  !! the pressure gradient.
  subroutine get_beta_prime(bp_out)
    use geometry, only: bishop, p_prime_input, invLp_input, beta_prime_input, alpha_input, dp_mult
    implicit none
    real, intent(out) :: bp_out

    select case (bishop)
    case(2)
       bp_out = p_prime_input
    case(3)
       bp_out = invLp_input
    case(4)
       bp_out = beta_prime_input
    case(5)
       bp_out = alpha_input
    case(6)
       bp_out = beta_prime_input
    case(7)
       bp_out = dp_mult
    case(8)
       bp_out = beta_prime_input
    case(9)
       bp_out = beta_prime_input
    case default
       bp_out = 0.0
       print*,"ERROR: For bishop not in 2-9 cannot modify the pressure gradient used"
    end select
  end subroutine get_beta_prime

  !> Sets the current value of beta_prime (or equivalent var)
  !! Really just sets the variable that we can use to control 
  !! the pressure gradient.
  subroutine set_beta_prime(bp_in)
    use geometry, only: bishop, p_prime_input, invLp_input, beta_prime_input, alpha_input, dp_mult
    implicit none
    real, intent(in) :: bp_in

    select case (bishop)
    case(2)
       p_prime_input = bp_in
    case(3)
       invLp_input = bp_in
    case(4)
       beta_prime_input = bp_in
    case(5)
       alpha_input = bp_in
    case(6)
       beta_prime_input = bp_in
    case(7)
       dp_mult = bp_in
    case(8)
       beta_prime_input = bp_in
    case(9)
       beta_prime_input = bp_in
    case default
       print*,"ERROR: For bishop not in 2-9 cannot modify the pressure gradient used"
    end select
  end subroutine set_beta_prime

  !> Runs the stability check scan
  subroutine run_stability_check
    use file_utils, only: open_output_file, close_output_file
    implicit none
    integer :: stab
    integer :: is, ib

    !Exit if not initialised
    if(.not.initialised) return

    !Check stability over range
    do ib=1,n_beta
       do is=1,n_shat
          call check_stability(shat_arr(is),beta_arr(ib),ib,stab)
          stability(is,ib)=stab
       enddo
    enddo
  end subroutine run_stability_check

  !> Routine to write out stability data to 1D ascii file
  subroutine write_stability_ascii_1d
    use file_utils, only: open_output_file, close_output_file
    implicit none
    integer :: is, ib, unit

    call open_output_file(unit,".ballstab_1d")
    write(unit,'(4(A12," "))') "shear", "beta_prime", "dbetadrho", "stability"
    do ib=1,n_beta
       do is=1,n_shat
          write(unit,'(3(F12.4," "),I12)') shat_arr(is), beta_arr(ib), dbdrho_arr(ib), stability(is,ib)
       enddo
    enddo
    call close_output_file(unit)
  end subroutine write_stability_ascii_1d

  !> Routine to write out stability data to 2D ascii file + 1d axis data
  subroutine write_stability_ascii_2d
    use file_utils, only: open_output_file, close_output_file
    implicit none
    integer :: is, ib, unit

    !/2D -- data
    call open_output_file(unit,".ballstab_2d")
    do ib=1,n_beta
       do is=1,n_shat
          write(unit,'(I0," ")',advance="no") stability(is,ib)
       enddo
       write(unit,'()') !Move to the next line
    enddo

    !/2D -- axis
    call open_output_file(unit,".ballstab_shat")
    !Commented lines here can be used to give row output
    do is=1,n_shat
       !write(unit,'(F9.2," ")',advance="no") shat_arr(is)
       write(unit,'(F9.2)') shat_arr(is)
    enddo
    !write(unit,'()')
    call close_output_file(unit)

    call open_output_file(unit,".ballstab_bp")
    do ib=1,n_beta
       write(unit,'(F12.4)') beta_arr(ib)
    enddo
    call close_output_file(unit)

    call open_output_file(unit,".ballstab_dbdr")
    do ib=1,n_beta
       write(unit,'(F12.4)') dbdrho_arr(ib)
    enddo
    call close_output_file(unit)
  end subroutine write_stability_ascii_2d

  !> Test if given shat/beta is unstable
  subroutine check_stability(shat_in, beta_prime_in, ib, iunstable, restore)
    !Note iunstable>0 means unstable - Could use a logical but size may be interesting
    use geometry, only: dbetadrho, eikcoefs
    implicit none
    real, intent(in) :: shat_in, beta_prime_in
    integer, intent(in) :: ib
    integer, intent(out) :: iunstable
    logical, intent(in), optional :: restore
    real :: shat_bak, beta_prime_bak
    
    !Store initial state
    call get_shat(shat_bak)
    call get_beta_prime(beta_prime_bak)

    !Setup geometry
    !NOTE: As how we can influence beta_prime (and shat) varies with bishop
    !we actually end up setting different parameters in different cases. The
    !physical parameter dbetadrho should really be the relevant parameter.
    call set_beta_prime(beta_prime_in)
    call set_shat(shat_in)
    call eikcoefs !Recalculate geometry terms

    !Here we store the value of dbetadrho, this should be bishop value independent
    dbdrho_arr(ib)=dbetadrho

    !Test stability
    iunstable=is_unstable()

    !Restore geometry if requested
    if(present(restore))then
       if(restore) then
          call set_beta_prime(beta_prime_bak)
          call set_shat(shat_bak)
          call eikcoefs
       endif
    endif
  end subroutine check_stability

  !> Test stability of current system -- return integer
  integer function is_unstable()
    !Note iunstable>0 means unstable - Could use a logical but size may be interesting
    !Assumes a call has been made to eikcoefs so that geometry terms are setup correctly
!<DD> NEED TO IMPROVE VARIABLE NAMES
    use geometry, only: gds2, gradpar, bmag, cvdrift, dbetadrho, theta
    implicit none
    real, dimension(:), allocatable :: g, c, psi, c1, c2, ch, g1, g2, gh, delthet
    integer :: ig, ntgrid
    real :: cflmax, cflmin, one_m_diff, psi_prime
    logical, parameter :: debug = .false.

    !Initialise
    is_unstable=0
    one_m_diff=1.-diff

    !Note we don't use theta_grid:ntgrid as in the case that
    !gridgen resizes the theta grid used by GS2 we don't end up
    !resizing any of the geometry module arrays (that we use directly here)
    !Really we would like to use everything from theta_grid directly but we'd
    !need some way to force theta_grid to update it's internal arrays when we
    !change parameters
    ntgrid=ubound(theta,1)
 
    !Calculate geometry arrays
    allocate(delthet(-ntgrid:ntgrid))
    delthet(:ntgrid-1) = theta(-ntgrid+1:) - theta(:ntgrid-1)

    allocate(g(-ntgrid:ntgrid),c(-ntgrid:ntgrid))
    g = gds2*gradpar/bmag
    c = -0.5*dbetadrho*cvdrift/(bmag*gradpar)

    allocate(g1(-ntgrid:ntgrid), c1(-ntgrid:ntgrid))
    allocate(g2(-ntgrid:ntgrid), c2(-ntgrid:ntgrid))
    allocate(gh(-ntgrid:ntgrid), ch(-ntgrid:ntgrid))
    allocate(psi(-ntgrid:ntgrid))

    do ig=-ntgrid+1, ntgrid
       ch(ig)=0.5*(c(ig)+c(ig-1))
       gh(ig)=0.5*(g(ig)+g(ig-1))
    enddo

    if(debug) then
       cflmax = delthet(-ntgrid)**2
       cflmax = cflmax*(c(-ntgrid+1)+c(-ntgrid))/(g(-ntgrid+1)+g(-ntgrid))
       cflmin = cflmax
       do ig=-ntgrid+1,ntgrid
          cflmax=max(cflmax,delthet(ig-1)**2*ch(ig)/gh(ig))
          cflmin=min(cflmin,delthet(ig-1)**2*ch(ig)/gh(ig))
       enddo
       write(6,'("|cflmax| = ",F12.5)') abs(cflmax)
       write(6,'("|cflmin| = ",F12.5)') abs(cflmin)
       write(6,'("These should be much less than 1")') !Could actually do a quantitative check
    endif

    !Calculate coefficients
    do ig=-ntgrid+1,ntgrid-1
       c1(ig) = -delthet(ig)*(one_m_diff*c(ig)+0.5*diff*ch(ig+1))&
            -delthet(ig-1)*(one_m_diff*c(ig)+0.5*diff*ch(ig))&
            -delthet(ig-1)*0.5*diff*ch(ig)
       c1(ig)=0.5*c1(ig)
    enddo

    !Calculate coefficients
    do ig=-ntgrid+1,ntgrid
       c2(ig) = -0.25*diff*ch(ig)*delthet(ig-1)
       g1(ig) = gh(ig)/delthet(ig-1)
       g2(ig) = 1.0/(0.25*diff*ch(ig)*delthet(ig-1)+gh(ig)/delthet(ig-1))
    enddo

    !Set boundary values
    psi(-ntgrid)=0.
    psi(-ntgrid+1)=delthet(-ntgrid)
    psi_prime=psi(-ntgrid+1)/g2(-ntgrid+1)!-g1(-ntgrid+1)*psi(-ntgrid) !Second term zero by defn

    !Integrate along in theta
    do ig=-ntgrid+1,ntgrid-1
       psi_prime=psi_prime+c1(ig)*psi(ig)+c2(ig)*psi(ig-1)
       psi(ig+1)=(g1(ig+1)*psi(ig)+psi_prime)*g2(ig+1)
    enddo

    !Find how many times psi crosses axis || Actually how many times it hits zero
    !If we really do just want to count crossings then need to replace le with lt
    do ig=-ntgrid+1,ntgrid-1
       if(psi(ig)*psi(ig+1) .le. 0 ) is_unstable=is_unstable+1
    enddo

    !Probably want to think about writing out psi if debug

    !Clean up
    deallocate(g,c,g1,c1,g2,c2,gh,ch)
  end function is_unstable
end module ballstab
