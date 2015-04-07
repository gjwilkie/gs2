!> Driver program for ballstab module, used to calculate ideal ballooning stability for a given (GS2) input file.
program ideal_ball
  use ballstab, only: init_ballstab, finish_ballstab, run_stability_check, write_stability_ascii_2d
  use job_manage, only: job_fork
  use mp, only: init_mp, finish_mp, proc0, broadcast
  use file_utils, only: init_file_utils, run_name
  use theta_grid, only: finish_theta_grid
  implicit none
  logical :: list
  character(len=500), target :: cbuff

  !Setup mpi
  call init_mp

  !Setup files
  if(proc0) call init_file_utils(list, name="gs")
  call broadcast(list)
  if(list) call job_fork

  !Set run name
  if(proc0) cbuff = trim(run_name)
  call broadcast(cbuff)
  if(.not.proc0) run_name => cbuff

  !Initialise
  call init_ballstab

  !Run
  if(proc0) then 
     call run_stability_check
     call write_stability_ascii_2d
  endif

  !Finish
  call finish_ballstab
  call finish_theta_grid
  call finish_mp
end program ideal_ball
