!> \mainpage
!! \section intro Introduction
!! 
!! GS2 is an initial value nonlinear code which solves the gyrokinetic equation. This is the source code documentation for GS2, and is aimed at developers. 
!!
!! \subsection udoc User Documentation
!! For a general introduction and user documentation please go to the Gyrokinetics Wiki: http://gyrokinetics.sourceforge.net/wiki/index.php/Main_Page
!! 
!! Some useful pages are: 
!!  - A beginners guide to downloading and installing: http://gyrokinetics.sourceforge.net/wiki/index.php/GS2:_A_Guide_for_Beginners
!!  - A general introduction to GS2: http://gyrokinetics.sourceforge.net/wiki/index.php/An_Introduction_to_GS2
!!  - A comprehensive guide to the GS2 input parameters: http://gyrokinetics.sourceforge.net/wiki/index.php/GS2_Input_Parameters
!!
!! \section doc Documentation Structure
!! Documentation is categorized by namespaces (i.e. modules), class (which in
!! fortran means custom types), and files. Within each namespace the subroutines
!! are documented in various ways, including developer comments and (very useful)
!! a chart showing which other subroutines  call the subroutine and which are
!! called by it.
!!
!! \section starting Starting Out
!!
!!  If you want to start at the beginning and work out what GS2 does, start at
!! gs2.f90, and use the call graphs to follow the algorithm.
!!
!! \section update Updating this documentation.
!! \subsection incode Updating Source Code Comments
!! This documentation is generated from commments added to the source code; to
!! add to it, add more comments to the <i>trunk</i> source. All documenting
!! comments begin with <tt>!></tt> and are continued with <tt>!!</tt>. For more help see the
!! doxygen help: http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html
!!
!! DO NOT add comments in between ifdef/endif preprocessor blocks. Such
!! comments will probably be ignored.
!!
!! DO NOT  add documentation for variables which also happen to be
!! input parameters. Your changes will be ignored in favour of the
!! automatically generated input parameter documentation. To edit 
!! this input parameter documentation, please go to the wiki page: 
!! http://gyrokinetics.sourceforge.net/wiki/index.php/Gs2_Input_Parameters. 
!! \subsection gen Updating this Documentation
!!
!! - Install doxygen: http://www.stack.nl/~dimitri/doxygen/
!! - Go to the trunk folder and type 
!! 
!!  <tt> make doc sync_doc USER=[sourceforge user name]</tt>

!> Main program. Used when running GS2
!! standalone, as opposed as a library for, e.g., Trinity.
!! Essentially this initializes a gs2_program_state_type 
!! object, and then calls the standard sequence of subroutines
!! from gs2_main to run the program. See gs2_main for more 
!! information.

program gs2

  ! make_lib is a compiler flag used if running with 
  ! an old version of trinity (coupled flux tube code)
  ! MAKE_LIB is now deprecated.

# ifndef MAKE_LIB 
  !use optimisation_config, only: optimisation_type
  use gs2_optimisation, only: initialize_gs2_optimisation
  use gs2_optimisation, only: finalize_gs2_optimisation
  use gs2_optimisation, only: optimise_gs2
  use gs2_main, only: gs2_program_state_type
  use gs2_main, only: initialize_wall_clock_timer
  use gs2_main, only: initialize_gs2
  use gs2_main, only: initialize_equations
  use gs2_main, only: initialize_diagnostics
  use gs2_main, only: evolve_equations
  use gs2_main, only: run_eigensolver
  use gs2_main, only: finalize_diagnostics
  use gs2_main, only: finalize_equations
  use gs2_main, only: finalize_gs2

  implicit none
  type(gs2_program_state_type) :: state
  !type(optimisation_type) :: optim
  call initialize_wall_clock_timer
  call initialize_gs2_optimisation(state)

  if (state%optim%on) then
    call optimise_gs2(state)
  end if
  if (state%optim%auto .or. .not. state%optim%on) then
    call initialize_gs2(state)
    call initialize_equations(state)
    call initialize_diagnostics(state)
    state%print_times = .false.
    if (state%do_eigsolve) then 
      call run_eigensolver(state)
    else
      call evolve_equations(state, state%nstep)
    end if
    call finalize_diagnostics(state)
    call finalize_equations(state)
    state%print_times = .true.
    state%print_full_timers = .true.
    call finalize_gs2(state)
  end if

  call finalize_gs2_optimisation(state)

# else

  implicit none
  call run_gs2

# endif



end program gs2
