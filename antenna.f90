module antenna
!!  
!! Originally by Hammett, Dorland and Quataert, Dec. 5, 2001
!!
!! Provides random forcing for Alfven wave problem
!!
!! init_antenna should be called once per run
!! antenna_amplitudes provides the necessary info per time step
!!
!!
!! two kinds of namelists should be added to the input file:
!! 
!! driver:
!!       
!!    amplitude : RMS amplitude of | apar_antenna |
!!    w_antenna : frequency of driving, normalized to kpar v_Alfven
!!                (assuming pk=2. or kp=1. in the theta_grid_parameters)
!!                Note: the imaginary part is the decorrelation rate
!!    nk_stir   : Number of k modes that should be driven
!!    write_antenna: .true. for direct antenna output to runname.antenna
!!                   default is .false.
!!
!! stir: (stir_1, stir_2, stir_3, etc., one for each k we want to drive)
!!
!!    kx, ky, kz : each of size nk_stir (not larger!)
!!               : list of k components that are to be driven; integers
!!    travel     : logical; choose false to make a standing wave
!!               : default value is .true.

  implicit none
  complex :: w_antenna
  complex, dimension(:), allocatable :: a_ant, b_ant, w_stir
  integer, dimension(:), allocatable :: kx_stir, ky_stir, kz_stir
  real :: amplitude
  integer :: nk_stir, out_unit
  logical :: no_driver = .false.
  logical :: write_antenna = .false.
  logical, dimension (:), allocatable :: trav
  logical :: ant_off = .false.

  private
  public :: init_antenna, antenna_amplitudes, dump_ant_amp

contains

  subroutine init_antenna
    use species, only: spec, nspec, init_species
    use run_parameters, only: beta
    use theta_grid, only: gradpar
    implicit none
    logical, save :: initialized = .false.
    real :: beta_s
    integer :: i

    if (initialized) return
    initialized = .true.

    call init_species

    call read_parameters
    if (no_driver) return

    allocate (w_stir(nk_stir))
    
    beta_s = 0.
    if (beta > epsilon(0.0)) then
       do i=1,nspec
          beta_s = beta_s + beta*spec(i)%dens*spec(i)%mass
       end do
    else
       beta_s = 2.
    end if

! assumes unsheared box 
    do i=1,nk_stir
       w_stir(i) = gradpar(0)*kz_stir(i)*sqrt(2./beta_s)*w_antenna
    end do

  end subroutine init_antenna
    
  subroutine read_parameters 
    use file_utils
    use mp, only: proc0, broadcast
    implicit none
    complex :: a, b
    integer :: kx, ky, kz, i, unit, in_file
    logical :: exist, travel

    namelist /driver/ amplitude, w_antenna, nk_stir, write_antenna, ant_off
    namelist /stir/ kx, ky, kz, travel, a, b

    w_antenna = (1., 0.0)
    amplitude = 0.
    nk_stir = 1
    ant_off = .false.

    if (proc0) then
       in_file=input_unit_exist("driver",exist)
       if (.not. exist) then
          no_driver = .true.
       else
	  if (ant_off) no_driver = .true.

          read (unit=input_unit("driver"), nml=driver)
          
          allocate (kx_stir(nk_stir))
          allocate (ky_stir(nk_stir))
          allocate (kz_stir(nk_stir))
          allocate (a_ant(nk_stir))
          allocate (b_ant(nk_stir))
          allocate (trav(nk_stir))
          
          do i=1,nk_stir
             call get_indexed_namelist_unit (in_file, "stir", i)
             kx=1
             ky=1
             kz=1
             a=0.
             b=0.
             travel = .true.
             read (unit=in_file, nml=stir)
             close(unit=in_file)
             kx_stir(i) = kx 
             ky_stir(i) = ky
             kz_stir(i) = kz
             trav(i) = travel
	     if (a == 0.0 .and. b == 0.0) then
                a_ant(i) = amplitude*cmplx(1.,1.)/2. 
                b_ant(i) = amplitude*cmplx(1.,1.)/2.
	     else 
                a_ant(i) = a
                b_ant(i) = b
             end if 
          end do
       end if
    end if

    call broadcast (no_driver)
    if (no_driver) return

    if (proc0) then
       call get_unused_unit (out_unit)
       call open_output_file (out_unit, ".antenna")
    end if

    call broadcast (w_antenna)
    call broadcast (amplitude)
    call broadcast (nk_stir)
    call broadcast (write_antenna)

    if (.not. proc0) then
       allocate (kx_stir(nk_stir))
       allocate (ky_stir(nk_stir))
       allocate (kz_stir(nk_stir))
       allocate (a_ant(nk_stir))
       allocate (b_ant(nk_stir))
       allocate (trav(nk_stir))
    end if
    
    call broadcast (trav)
    call broadcast (kx_stir)
    call broadcast (ky_stir)
    call broadcast (kz_stir)
    call broadcast (a_ant)
    call broadcast (b_ant)
       
  end subroutine read_parameters

  subroutine antenna_amplitudes (apar)

    use mp
    use gs2_time, only: simdt, stime
    use kt_grids, only: naky, ntheta0, reality
    use run_parameters, only: tnorm
    use theta_grid, only: theta, ntgrid
    use ran
    use constants

    complex, dimension (-ntgrid:,:,:) :: apar
    complex :: force_a, force_b
    real :: dt, sigma, time
    integer :: i, it, ik
    
    if (no_driver) return

    dt = simdt()/tnorm
    time = stime()/tnorm

    apar = 0.

    do i=1, nk_stir
       
       it = 1 + mod(kx_stir(i)+ntheta0,ntheta0)
       ik = 1 + mod(ky_stir(i)+naky,naky)

       sigma = sqrt(12.*abs(aimag(w_stir(i)))/dt)*amplitude

       force_a = cmplx(ranf() - 0.5, ranf() - 0.5) * sigma
       if (trav(i)) then
          force_b = cmplx(ranf() - 0.5, ranf() - 0.5) * sigma
       else
          force_b = force_a
       end if

       if (trav(i) .and. abs(aimag(w_stir(i))) < epsilon(0.0)) then
          a_ant(i) = a_ant(i)*exp(-zi*w_stir(i)*dt)+force_a*dt
          b_ant(i) = 0.
       else
          a_ant(i) = a_ant(i)*exp(-zi*w_stir(i)*dt)+force_a*dt
          b_ant(i) = b_ant(i)*exp(zi*conjg(w_stir(i))*dt)+force_b*dt
       end if

       apar(:,it,ik) = apar(:,it,ik) &
            + (a_ant(i)+b_ant(i))/sqrt(2.)*exp(zi*kz_stir(i)*theta) 

       if (write_antenna) then
         if (proc0) write(out_unit, fmt='(7(1x,e13.6))') &
	 & time, real((a_ant(i)+b_ant(i))/sqrt(2.)), &
         &     aimag((a_ant(i)+b_ant(i))/sqrt(2.)),real(a_ant(i)), aimag(a_ant(i)), &
         &     real(b_ant(i)), aimag(b_ant(i))
       end if
    end do

    if (reality) then
       apar(:,1,1) = 0.0
       
       do it = 1, ntheta0/2
          apar(:,it+(ntheta0+1)/2,1) = conjg(apar(:,(ntheta0+1)/2+1-it,1))
       enddo
    end if

  end subroutine antenna_amplitudes

  subroutine dump_ant_amp
    use gs2_time, only: simdt, stime
    use run_parameters, only: tnorm
    use mp, only: proc0
    implicit none
    real :: time
    integer :: i

    if (no_driver) return
    if (.not. write_antenna) then
       time = stime()/tnorm
       do i=1,nk_stir
         if (proc0) write(out_unit, fmt='(7(1x,e13.6))') &
	 & time, real((a_ant(i)+b_ant(i))/sqrt(2.)), &
         &     aimag((a_ant(i)+b_ant(i))/sqrt(2.)),real(a_ant(i)), aimag(a_ant(i)), &
         &     real(b_ant(i)), aimag(b_ant(i))
  
       end do
    end if

  end subroutine dump_ant_amp

end module antenna

