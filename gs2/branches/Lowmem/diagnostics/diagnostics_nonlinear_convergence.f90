!> A module which contains methods for checking whether nonlinear
!! runs have reached a saturated steady state
module diagnostics_nonlinear_convergence
  implicit none

  private

  public :: check_nonlin_convergence
  public :: init_nonlinear_convergence
  public :: finish_nonlinear_convergence

  integer :: trin_istep = 0
  integer :: conv_isteps_converged = 0
  real, save, allocatable, dimension(:) :: conv_heat
  real, save :: heat_sum_av = 0, heat_av = 0, heat_av_test = 0
contains

  subroutine init_nonlinear_convergence(gnostics)
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    if(.not. allocated(conv_heat)) allocate(conv_heat(0:gnostics%conv_nstep_av/gnostics%nwrite-1))
  end subroutine  init_nonlinear_convergence
  
  subroutine finish_nonlinear_convergence(gnostics)
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    deallocate(conv_heat)
  end subroutine  finish_nonlinear_convergence

  !> Trinity convergence condition - simple and experimental
  !! look for the averaged differential of the summed averaged heat flux to drop below a threshold
  subroutine check_nonlin_convergence(gnostics)
    use job_manage, only: trin_job
    use gs2_time, only: user_time
    use mp, only: proc0, broadcast
    use diagnostics_config, only: diagnostics_type
    implicit none
    type (diagnostics_type), intent(inout) :: gnostics
    
    logical :: exit
    integer :: istep, nwrite
    real :: heat_flux

! Variables for convergence condition (HJL)
    real :: heat_av_new, heat_av_diff
    integer :: place, iwrite, nwrite_av
    logical :: debug = .false.

    istep = gnostics%istep
    nwrite = gnostics%nwrite
    heat_flux = gnostics%current_results%total_heat_flux
    exit = gnostics%exit
    
    if(proc0 .and. .not. (trin_istep .ne. 0 .and. istep .eq. 0)) then
       if(istep .gt. 0) trin_istep = trin_istep + nwrite ! Total number of steps including restarted trinity runs
       iwrite = trin_istep/nwrite ! Number of diagnostic write steps written
       nwrite_av = gnostics%conv_nstep_av/nwrite ! Number of diagnostic write steps to average over        
       
       heat_sum_av = (heat_sum_av * iwrite + heat_flux) / (iwrite+1) ! Cumulative average of heat flux
       
       place = mod(trin_istep,gnostics%conv_nstep_av)/gnostics%nwrite
       conv_heat(place) = heat_flux
       
       if (debug) write(6,'(A,I5,A,e11.4,A,I6,I6)') 'Job ',trin_job, &
            ' time = ',user_time, ' step = ',trin_istep
       if (debug) write(6,'(A,I5,A,e11.4,A,e11.4)') 'Job ',trin_job, &
            ' heat = ',heat_flux, ' heatsumav = ',heat_av
       
       if (trin_istep .ge. gnostics%conv_nstep_av) then
          heat_av_new = sum(conv_heat) / nwrite_av
          heat_av_diff = heat_av_new - heat_av
          if(debug) write(6,'(A,I5,A,e11.4,A,e11.4)') 'Job ',trin_job, &
               ' heat_sum_av_diff = ',heat_sum_av
          heat_av = heat_av_new
          ! Convergence test - needs to be met conv_nsteps_converged/nwrite times in succession
          if (abs(heat_av_diff) .lt. heat_av_test) then
             conv_isteps_converged = conv_isteps_converged + 1
          else
             conv_isteps_converged = 0
             heat_av_test = heat_sum_av * gnostics%conv_test_multiplier
          endif
          
          if ((conv_isteps_converged .ge. gnostics%conv_nsteps_converged/nwrite) .and. &
               (trin_istep .ge. gnostics%conv_min_step)) then
             if (debug) write(6,'(A,I5,A,I6,I3)')'Job ',trin_job, &  
                  ' Reached convergence condition after step ',trin_istep
             exit = .true. .and. gnostics%exit_when_converged
          endif
          
          if(trin_istep .gt. gnostics%conv_max_step) then
             write(6,'(A,I5,A,I7)') '*** Warning. Job ',trin_job, &
                  ' did not meet the convergence condition after ',trin_istep
             exit = .true. .and. gnostics%exit_when_converged
          endif
       endif
    endif
    
    call broadcast(exit)
    
    gnostics%exit = exit

  end subroutine check_nonlin_convergence
end module diagnostics_nonlinear_convergence
