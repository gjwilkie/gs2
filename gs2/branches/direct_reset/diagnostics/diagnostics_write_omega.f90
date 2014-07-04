!> A module which calculates and writes the growth rates and frequencies
module diagnostics_write_omega

  use diagnostics_config, only: diagnostics_type
  use kt_grids, only: ntheta0, naky
  
  implicit none

  public :: write_omega
  public :: init_diagnostics_write_omega
  public :: finish_diagnostics_write_omega
  public :: calculate_omega
  public :: omega_average
  public :: debug

  complex, dimension (:,:), allocatable :: omega_average
  logical :: debug =.false.

  private
  complex, dimension (:,:,:), allocatable :: omegahist
  complex, allocatable, save, dimension (:,:,:) :: domega


contains 
  !> Allocate arrays for storing omega history and averages
  subroutine init_diagnostics_write_omega(gnostics)
    type(diagnostics_type), intent(in) :: gnostics
     allocate (omegahist(0:gnostics%navg-1,ntheta0,naky))
     allocate(domega(gnostics%navg,ntheta0,naky))
     allocate(omega_average(ntheta0, naky))
     omegahist = 0.0
     omega_average = 0.0
  end subroutine init_diagnostics_write_omega

  subroutine finish_diagnostics_write_omega
    deallocate(omegahist)
    deallocate(domega)
    deallocate(omega_average)
  end subroutine finish_diagnostics_write_omega

  !> Write omega, omega_average as functions of kx, ky and time
  !! to the new netcdf file
  subroutine write_omega(gnostics)
    use simpledataio
    use simpledataio_write
    use fields_parallelization, only: field_k_local
    type(diagnostics_type), intent(in) :: gnostics
    integer :: ik, it


    ! Here we have to separate creating and writing because the fields
    ! may be distributed
    if (gnostics%create) then
      call create_variable(gnostics%sfile, gnostics%rtype, "omega", "rXYt", &
        "Complex frequency as a function of kx, ky and time", "a/v_thr")
      call create_variable(gnostics%sfile, gnostics%rtype, "omega_average", "rXYt", &
        "Complex frequency, averaged over navg timesteps, as a function of kx, ky and time", "a/v_thr")
     call set_count(gnostics%sfile, "omega", "X", 1)
     call set_count(gnostics%sfile, "omega", "Y", 1)
     call set_count(gnostics%sfile, "omega_average", "X", 1)
     call set_count(gnostics%sfile, "omega_average", "Y", 1)
    end if

    
    if (gnostics%wryte) then
      ! Dummy writes, necessary for parallel for some reason
      call write_variable_with_offset(gnostics%sfile, "omega", omegahist(mod(gnostics%istep,gnostics%navg), :, :))
      call write_variable_with_offset(gnostics%sfile, "omega_average", omega_average)

      do ik=1,naky
        do it=1,ntheta0
          if (.not. gnostics%distributed .or. field_k_local(it, ik)) then 
             call set_start(gnostics%sfile, "omega", "X", it)
             call set_start(gnostics%sfile, "omega", "Y", ik)
             call write_variable_with_offset(gnostics%sfile, "omega", omegahist(mod(gnostics%istep,gnostics%navg), :, :))

             call set_start(gnostics%sfile, "omega_average", "X", it)
             call set_start(gnostics%sfile, "omega_average", "Y", ik)
             call write_variable_with_offset(gnostics%sfile, "omega_average", omega_average)
          end if
        end do
      end do
    end if

   
  end subroutine write_omega

  !> Calculates omega and stores it in omegahist. Calculates omega_average, 
  !! the average of omega over navg timesteps
  subroutine calculate_omega (gnostics)
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use gs2_time, only: code_dt
    use constants, only: zi
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics

    real :: fac
    integer :: j

    if (debug) write(6,*) "get_omeaavg: start"

    j = gnostics%igomega
    if(gnostics%omegatol.eq.0)then
       fac=1.0
    else
       fac=1000/abs(gnostics%omegatol)
    endif

    !<DD>The logic below was originally
    !    where (abs(phinew(j,:,:)+aparnew(j,:,:)+bparnew(j,:,:)) < epsilon(0.0) &
    !           .or. abs(phi(j,:,:)+apar(j,:,:)+bpar(j,:,:)) < epsilon(0.0))
    !This results in not calculating the frequency whenever the fields drop below epsilon(0.0) [~1d-16 for DP]
    !What we really want to do is prevent dividing by too small a number. 
    ! Calculate omega
    where (abs(phinew(j,:,:)+aparnew(j,:,:)+bparnew(j,:,:)) < tiny(0.0)*fac &
           .or. abs(phi(j,:,:)+apar(j,:,:)+bpar(j,:,:)) < tiny(0.0)*fac)
       omegahist(mod(gnostics%istep,gnostics%navg),:,:) = 0.0
    elsewhere
       omegahist(mod(gnostics%istep,gnostics%navg),:,:) &
            = log((phinew(j,:,:) + aparnew(j,:,:) + bparnew(j,:,:)) &
                  /(phi(j,:,:)   + apar(j,:,:)    + bpar(j,:,:)))*zi/code_dt
    end where

    ! Calculate omega_average
    omega_average = sum(omegahist/real(gnostics%navg),dim=1)

    if (debug) write(6,*) "calculate_omega: omega_average=",omega_average

    ! Check if converged by comparing the difference between
    ! omega_average and all values in omegahist, to the value of omega_average
    if (gnostics%istep > gnostics%navg) then
       domega = spread(omega_average,1,gnostics%navg) - omegahist
       if (all(sqrt(sum(abs(domega)**2/real(gnostics%navg),dim=1)) &
            .le. min(abs(omega_average),1.0)*gnostics%omegatol)) &
       then
          write (*, "('*** omega converged')")
          gnostics%exit = gnostics%exit_when_converged
       end if

       if (any(abs(omega_average)*code_dt > gnostics%omegatinst)) then
          write (*, "('*** numerical instability detected')") 
          gnostics%exit = .true.
       end if
    end if

    if (debug) write(6,*) "calculate_omega: done"
  end subroutine calculate_omega

end module diagnostics_write_omega
