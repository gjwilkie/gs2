!> This module sits on top of gs2_main and provides
!! a toolkit for measuring and optimising performance
module gs2_optimisation
  use gs2_main, only: gs2_program_state_type


  contains

    subroutine initialize_gs2_optimisation(state)
      use optimisation_config, only: init_optimisation_config
      use gs2_main, only: initialize_gs2, finalize_gs2
      use gs2_main, only: initialize_wall_clock_timer
      use mp, only: init_mp, mp_comm
      type(gs2_program_state_type), intent(inout) :: state
      call init_mp
      state%mp_comm_external = .true.
      state%mp_comm = mp_comm
      allocate(state%optim%sorted_optimisations(0))
      allocate(state%optim%sorted_results(0))
      !write (*,*) 'size', size(state%optim%sorted_results)
      ! We have to initialize_gs2 so that we 
      ! can read the optimisation_config namelist
      !call initialize_wall_clock_timer
      call initialize_gs2(state)
      call init_optimisation_config(state%optim)
      call finalize_gs2(state)
    end subroutine initialize_gs2_optimisation

    subroutine finalize_gs2_optimisation(state)
      use mp, only: finish_mp
      use optimisation_config, only: finish_optimisation_config
      type(gs2_program_state_type), intent(inout) :: state
      call finish_optimisation_config(state%optim)
      call finish_mp
    end subroutine finalize_gs2_optimisation

    subroutine optimise_gs2(state)
      use mp, only: proc0, mp_abort
      use fields_implicit, only: skip_initialisation
      type(gs2_program_state_type), intent(inout) :: state
      integer i,n
      !real :: optcost, opttime


      state%print_times = .false.
      state%print_full_timers = .false.
      state%is_external_job = .true.
      ! Initialize optimisation results
      !call measure_timestep(state)
      state%optim%results%optimal_time = -1.0
      state%optim%results%optimal_cost = -1.0
      state%optim%results%optimal = .false.

      skip_initialisation = .true.

      if (state%optim%warm_up) then
      
        call optimise_layout(state)
      
        deallocate(state%optim%sorted_results)
        deallocate(state%optim%sorted_optimisations)
        allocate(state%optim%sorted_optimisations(0))
        allocate(state%optim%sorted_results(0))
        state%optim%results%optimal_time = -1.0
        state%optim%results%optimal_cost = -1.0
        state%optim%results%optimal = .false.
      end if

      call optimise_layout(state)
      skip_initialisation = .false.
      call output_results(state)

      if (state%optim%auto) then 
      ! Find the optimal configuration which satisfies
      ! constraints. Abort if one can't be found.

        n = size(state%optim%sorted_optimisations)
        do i = 1,n
         ! Find the most optimal configuration that satisfies
         ! max_unused_procs and max_imbalance and min_efficiency
          if (&
               (state%optim%max_unused_procs .lt. 0 .or.&
               (state%optim%nproc_max - state%optim%sorted_results(i)%nproc) .le.&
                state%optim%max_unused_procs)  &
               .and. &
               (state%optim%max_imbalance .lt. 0.0 .or.&
               (state%optim%nproc_max - state%optim%sorted_results(i)%nproc) / &
                state%optim%nproc_max .le.&
                 state%optim%max_imbalance) &
               .and. &
               (state%optim%min_efficiency .lt. 0.0 .or. &
                state%optim%sorted_results(i)%efficiency .gt. &
                  state%optim%min_efficiency) & 
               ) exit
          if (i .eq. n) then
            call mp_abort("Could not satisfy min_efficiency without &
              & violating max_imbalance or max_unused_procs", .true.)
          end if
        end do



        !> This is the line which optimises GS2, by copying 
        !! the optimal set of overrides into the init structure
        state%init%opt_ov = state%optim%sorted_optimisations(i)

      end if


    end subroutine optimise_gs2

    subroutine output_results(state)
      use mp, only: proc0
      implicit none
      type(gs2_program_state_type), intent(inout) :: state
      integer :: i,n
       

      if (proc0) then
        write (state%optim%outunit, '(A10," ",A10," ",A10," ",A6," ",A6)') &
          'wallclocktime', 'efficiency', 'cost', 'nproc', 'layout' 
        n = size(state%optim%sorted_results)
        do i = 1,n
            write(state%optim%outunit, &
              '(E10.4," ",F10.6," ",E10.4," ",I6," ",A6)') &
            state%optim%sorted_results(i)%time, &
            state%optim%sorted_results(i)%efficiency, &
            state%optim%sorted_results(i)%cost, &
            state%optim%sorted_results(i)%nproc, &
            state%optim%sorted_optimisations(i)%layout
        end do
      end if

    end subroutine output_results

    subroutine optimise_nprocs(state)
      use ingen_mod, only: init_ingen, finish_ingen, report
      use ingen_mod, only: sweet_spots, n_sweet_spots
      use gs2_main, only: initialize_gs2, initialize_equations, initialize_diagnostics
      use gs2_main, only: finalize_gs2, finalize_equations, finalize_diagnostics
      type(gs2_program_state_type), intent(inout) :: state
      integer, dimension (4) :: nproc_values = (/8,4,2,1/)
      integer :: i

      state%init%opt_ov%override_nproc = .false.
      ! First measure performance using all procs
      call optimise_flags(state)

      call init_ingen
      call initialize_gs2(state)
      call initialize_equations(state)
      call initialize_diagnostics(state)
      call report
      call finalize_diagnostics(state)
      call finalize_equations(state)
      call finalize_gs2(state)

      ! Loop through all sweet spots and measure performance
      do i = 1,n_sweet_spots
        if (sweet_spots(i)%nproc .gt. state%optim%nproc_max) exit
        ! If asked to check for inefficencies, check all proc numbers
        ! otherwise only check proc numbers that satisfy
        ! max_imbalance and max_unused_procs
        if (.not. (state%optim%min_efficiency .gt. 0)) then
          if ( state%optim%max_unused_procs .gt. 0 .and. &
              (state%optim%nproc_max - sweet_spots(i)%nproc) .gt.&
                state%optim%max_unused_procs ) cycle
          if ( state%optim%max_imbalance .gt. 0.0 .and. &
              (state%optim%nproc_max - sweet_spots(i)%nproc) / &
               state%optim%nproc_max .gt.&
                state%optim%max_imbalance ) cycle
        end if
        state%init%opt_ov%override_nproc = .true.
        state%init%opt_ov%nproc = sweet_spots(i)%nproc
        call optimise_flags(state)
      end do
      call finish_ingen
      !end if


    end subroutine optimise_nprocs

    subroutine optimise_layout(state)
      use gs2_main, only: prepare_optimisations_overrides
      use mp, only: proc0
      type(gs2_program_state_type), intent(inout) :: state
      call prepare_optimisations_overrides(state)
      state%init%opt_ov%override_layout = .true.
      state%init%opt_ov%layout = 'lxyes'
      call optimise_nprocs(state)
      state%init%opt_ov%layout = 'lexys'
      call optimise_nprocs(state)
      state%init%opt_ov%layout = 'xyles'
      call optimise_nprocs(state)
      state%init%opt_ov%layout = 'yxles'
      call optimise_nprocs(state)
    end subroutine optimise_layout

    subroutine optimise_flags(state)
      type(gs2_program_state_type), intent(inout) :: state
      call measure_timestep(state)
    end subroutine optimise_flags

    subroutine measure_timestep(state)
      use gs2_main, only: gs2_program_state_type
      use gs2_main, only: initialize_gs2
      use gs2_main, only: initialize_equations
      use gs2_main, only: initialize_diagnostics
      use gs2_main, only: evolve_equations
      use gs2_main, only: run_eigensolver
      use gs2_main, only: finalize_diagnostics
      use gs2_main, only: finalize_equations
      use gs2_main, only: finalize_gs2
      use overrides, only: optimisations_overrides_type
      use optimisation_config, only: optimisation_results_type
      use mp, only: proc0, broadcast, mp_abort

      implicit none
      type(gs2_program_state_type), intent(inout) :: state
      type(optimisations_overrides_type), &
        dimension(:), allocatable :: sorted_opts_temp
      type(optimisation_results_type), &
        dimension(:), allocatable :: sorted_res_temp
      integer :: i,n
      real :: t, cost
      logical :: completed_steps

      completed_steps = .true.

      call initialize_gs2(state)
      call initialize_equations(state)
      call initialize_diagnostics(state)
      call evolve_equations(state, state%optim%nstep_measure)
      if (state%included .and. state%istep_end .ne. state%optim%nstep_measure) then
        completed_steps = .false.
        write(*,*) 'istep_end', state%istep_end, state%optim%nstep_measure
      end if
      call finalize_diagnostics(state)
      call finalize_equations(state)
      call finalize_gs2(state)
      if (.not. completed_steps) &
        call mp_abort('Optimisation has failed because gs2 is not completing &
          & the required number of steps. It may be hitting a convergence &
          & criterion, or a time limit, or it may be a numerical instability. &
          & Check exit_when_converged, avail_cpu_time, omegatol, omegatinst.',&
          .true.)


      if (state%optim%measure_all) then 
        t = state%timers%advance(1)/real(state%optim%nstep_measure)
      else
        t = state%timers%timestep(1)/real(state%optim%nstep_measure)
        !t = state%timers%timestep(1)
      endif
      cost = t*real(state%nproc_actual)
      call broadcast(t)
      call broadcast(cost)
      state%optim%results%nproc = state%nproc_actual
      call broadcast(state%optim%results%nproc)


      state%optim%results%time = t
      state%optim%results%cost = cost
      !state%optim%results%cost = t
      if (t .lt.  state%optim%results%optimal_time .or. &
        state%optim%results%optimal_time .lt. 0.0) then
        state%optim%results%optimal_time = t
        state%optim%results%optimal = .true.
      end if
      if (cost .lt.  state%optim%results%optimal_cost .or. &
        state%optim%results%optimal_cost .lt. 0.0) then
        state%optim%results%optimal_cost = cost
        !if (proc0) write(*,*) 'optimal_cost', state%optim%results%optimal_cost, cost
      end if

      n = size(state%optim%sorted_results)
      !write (*,*) 'size2', size(state%optim%sorted_results)
      allocate(sorted_opts_temp(n), sorted_res_temp(n))
      if (n>0) then
        do i = 1,n
          sorted_opts_temp(i) = state%optim%sorted_optimisations(i)
          sorted_res_temp(i) = state%optim%sorted_results(i)
          end do
      end if

      deallocate(state%optim%sorted_optimisations)
      deallocate(state%optim%sorted_results)
      allocate(state%optim%sorted_optimisations(n+1))
      allocate(state%optim%sorted_results(n+1))

      i=1
      do 
        if (i>n) exit
        if (sorted_res_temp(i)%time > t) exit
        state%optim%sorted_optimisations(i) = sorted_opts_temp(i)
        state%optim%sorted_results(i) = sorted_res_temp(i)
        i = i+1
      end do
      state%optim%sorted_optimisations(i) = state%init%opt_ov
      state%optim%sorted_results(i) = state%optim%results
      !if (proc0) write(*,*) 'cost', state%optim%sorted_results(i)%cost, cost
      i = i + 1
      do while (i < n+2)
        state%optim%sorted_optimisations(i) = sorted_opts_temp(i-1)
        state%optim%sorted_results(i) = sorted_res_temp(i-1)
        i = i + 1
      end do
      do i = 1,size(state%optim%sorted_results)
        state%optim%sorted_results(i)%optimal_cost = &
          state%optim%results%optimal_cost
        state%optim%sorted_results(i)%optimal_time = &
          state%optim%results%optimal_time
        state%optim%sorted_results(i)%efficiency = &
          state%optim%sorted_results(i)%optimal_cost / &
          state%optim%sorted_results(i)%cost
      end do
        
      deallocate(sorted_opts_temp, sorted_res_temp)
      
    end subroutine measure_timestep

end module gs2_optimisation
