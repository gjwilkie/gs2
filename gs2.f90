program gs2
  use mp, only: init_mp, finish_mp, proc0
  use file_utils, only: init_file_utils, finish_file_utils
  use fields, only: init_fields
  use run_parameters, only: init_run_parameters
  use gs2_diagnostics, only: init_gs2_diagnostics, finish_gs2_diagnostics

  use run_parameters, only: nstep
  use fields, only: advance
  use gs2_diagnostics, only: loop_diagnostics

  implicit none
  integer :: istep
  logical :: exit

  call init_mp
  if (proc0) call init_file_utils (name="gs")
  call init_fields
  call init_gs2_diagnostics

  do istep = 1, nstep
     call advance (istep)
     call loop_diagnostics (istep, nstep, exit)
     if (exit) exit
  end do

  call finish_gs2_diagnostics
  if (proc0) call finish_file_utils
  call finish_mp

end program gs2
