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
      use fields_local, only: fieldmat
      implicit none
      type(gs2_program_state_type), intent(inout) :: state
      integer i,n
      real,dimension(10) :: time_array
      real :: mean,sd
      !real :: optcost, opttime


      !state%report_nprocs = .false.
      state%print_times = .false.
      state%print_full_timers = .false.
      state%is_external_job = .true.
      ! Initialize optimisation results
      !call measure_timestep(state)
      state%optim%results%optimal_time = -1.0
      state%optim%results%optimal_cost = -1.0
      state%optim%results%optimal = .false.

      skip_initialisation = .true.
      fieldmat%no_prepare = .true.
      fieldmat%no_populate = .true.
      state%dont_change_timestep = .true.

      state%init%tstep_ov%init = .true.
      state%init%tstep_ov%override_immediate_reset = .true.
      state%init%tstep_ov%immediate_reset = .false.

      if (state%optim%estimate_timing_error) then
        do i = 1,10
          call measure_timestep(state)
          time_array(i) = state%optim%results%time
        end do
        mean = sum(time_array(1:10)) / real(10)
        sd = sqrt (sum((time_array(1:10)-mean)**2) / real(10.0)) 
        state%optim%timing_rel_error = sd/mean
        state%optim%timing_max_rel_error = &
          (maxval(time_array)-minval(time_array))/mean

        write (*,*) 'Timing', mean, sd, sd/mean
        deallocate(state%optim%sorted_results)
        deallocate(state%optim%sorted_optimisations)
        allocate(state%optim%sorted_optimisations(0))
        allocate(state%optim%sorted_results(0))
        state%optim%results%optimal_time = -1.0
        state%optim%results%optimal_cost = -1.0
        state%optim%results%optimal = .false.
      else
        state%optim%timing_rel_error = -1.0

      end if
      if (state%optim%warm_up) then
      
        call optimise_simple(state)
      
        deallocate(state%optim%sorted_results)
        deallocate(state%optim%sorted_optimisations)
        allocate(state%optim%sorted_optimisations(0))
        allocate(state%optim%sorted_results(0))
        state%optim%results%optimal_time = -1.0
        state%optim%results%optimal_cost = -1.0
        state%optim%results%optimal = .false.
      end if

      call optimise_simple(state)

      skip_initialisation = .false.
      fieldmat%no_prepare = .false.
      fieldmat%no_populate = .false.
      state%dont_change_timestep = .false.

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
      character(len=*), parameter :: formt =  &
        '(A10," ",A10," ",A10," ",A6," ",A6," ",A1," ",A1," ",A1," ",A1," ",A1,&
        &" ",A1," ",A1," ",A1," ",A1," ",A1," ",A1," ",A1," ",A7)'
      character(len=*), parameter :: bk = '               '
      character(len=*), parameter :: ul = '---------------'
      character(len=*), parameter :: h1 = "GS2 Timing"
      character(len=*), parameter :: h2 = "Data      "
      character(len=*), parameter :: h3 = " Est. Rel."
      character(len=*), parameter :: ha = " Est. Max."
      character(len=*), parameter :: h4 = "Timing Err"
      character(len=*), parameter :: h5 = " %        "
      character(len=6) :: er
      character(len=6) :: em
      integer :: i,n, ou

      ou = state%optim%outunit

      write(er, "(F6.2)") state%optim%timing_rel_error * 100.0
      write(em, "(F6.2)") state%optim%timing_max_rel_error * 100.0
       

      if (proc0) then

write(ou,formt)h1,h2,bk,bk,bk,'o', bk, bk, bk, bk, bk, bk, bk,'d','f', bk, bk,bk
write(ou,formt)h3,h4,er,h5,bk,'p', bk, bk,'l', bk, bk, bk,'f','|','l', bk, bk,bk
write(ou,formt)ha,h4,em,h5,bk,'t', bk, bk,'o', bk, bk, bk,'i','s','o', bk,'f',bk
write(ou,formt)bk,bk,bk,bk,bk,'|', bk, bk,'c', bk, bk,'i','e','m','c', bk,'i',bk
write(ou,formt)bk,bk,bk,bk,bk,'r', bk, bk,'a','o','i','n','l','a','|', bk,'e',bk
write(ou,formt)bk,bk,bk,bk,bk,'e', bk, bk,'l','p','n','t','d','r','a', bk,'l',bk
write(ou,formt)bk,bk,bk,bk,bk,'d','|','|','|','t','t','s','|','t','l', bk,'d',bk
write(ou,formt)bk,bk,bk,bk,bk,'i','p','o','f','|','m','p','s','|','l', bk,'|',bk
write(ou,formt)bk,bk,bk,bk,bk,'s','e','v','|','s','o','e','u','u','r', bk,'o',bk
write(ou,formt)bk,bk,bk,bk,bk,'t','r','e','s','o','m','c','b','p','e', bk,'p',bk
write(ou,formt)bk,bk,bk,bk,bk,'|','s','r','o','u','|','|','g','d','d','|','t',bk
write(ou,formt)bk,bk,bk,bk,bk,'n','i','l','l','r','s','s','a','a','u','s','i',bk
write(ou,formt)bk,bk,bk,bk,bk,'b','s','a','v','c','u','u','t','t','c','u','o',bk
      write (ou, formt) &
        'wallclocktime', 'efficiency', 'cost', 'nproc', 'layout', &
                              'k','t','p','e','e','b','b','h','e','e','b','n',&
                                    'minnrow'
      write (ou,formt) ul, ul, ul, ul, ul, ul, ul, ul, ul, ul, ul, ul, ul, ul,&
                       ul, ul, ul, ul
                        
      n = size(state%optim%sorted_results)
      do i = 1,n
        call write_summary(state%optim%outunit,&
          state%optim%sorted_results(i), &
          state%optim%sorted_optimisations(i))
      end do
      end if

    end subroutine output_results

    subroutine write_summary(unt, results, optimisations)
      use optimisation_config, only: optimisation_results_type
      use overrides, only: optimisations_overrides_type
      implicit none
      integer, intent(in) :: unt
      type(optimisation_results_type), intent(in) :: results
      type(optimisations_overrides_type), intent(in) :: optimisations
      write(unt, &
        '(E10.4," ",F10.6," ",E10.4," ",I6," ",A6," ",&
        &L1," ",L1," ",L1," ",L1," ",L1," ",L1," ",L1," ",&
        &L1," ",L1," ",L1," ",L1," ",A1," ",I7)') &
      results%time, &
      results%efficiency, &
      results%cost, &
      results%nproc, &
      optimisations%layout, &
      optimisations%opt_redist_nbk, &
      optimisations%opt_redist_persist, &
      optimisations%opt_redist_persist_overlap,&
      optimisations%local_field_solve, &
      optimisations%opt_source, &
      optimisations%intmom_sub,&
      optimisations%intspec_sub,&
      optimisations%field_subgath,&
      optimisations%do_smart_update,&
      optimisations%field_local_allreduce,&
      optimisations%field_local_allreduce_sub,&
      optimisations%field_option(1:1), &
      optimisations%minnrow
    end subroutine write_summary
    
    subroutine optimise_simple(state)
      use gs2_main, only: prepare_optimisations_overrides
      type(gs2_program_state_type), intent(inout) :: state
      call prepare_optimisations_overrides(state)


      call optimise_layout(state)
    end subroutine  optimise_simple

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
      call optimise_fields(state)

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
          if ( state%optim%max_unused_procs .ge. 0 .and. &
              (state%optim%nproc_max - sweet_spots(i)%nproc) .gt.&
                state%optim%max_unused_procs ) cycle
          if ( state%optim%max_imbalance .gt. 0.0 .and. &
              (state%optim%nproc_max - sweet_spots(i)%nproc) / &
               state%optim%nproc_max .gt.&
                state%optim%max_imbalance ) cycle
        end if
        state%init%opt_ov%override_nproc = .true.
        state%init%opt_ov%nproc = sweet_spots(i)%nproc
        call optimise_fields(state)
      end do
      call finish_ingen
      !end if


    end subroutine optimise_nprocs

    subroutine optimise_layout(state)
      use mp, only: proc0
      type(gs2_program_state_type), intent(inout) :: state
      !> Measure default layout
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

      logical :: l1, l2, l3

      l1=.false.
      l2=.false.
      l3=.false.
      state%init%opt_ov%override_opt_redist_nbk = .true.
      state%init%opt_ov%opt_redist_nbk = .false.
      state%init%opt_ov%override_opt_redist_persist = .true.
      state%init%opt_ov%opt_redist_persist = .false.
      state%init%opt_ov%override_opt_redist_persist_overlap = .true.
      state%init%opt_ov%opt_redist_persist_overlap = .false.

      state%init%opt_ov%override_local_field_solve = .true.
      state%init%opt_ov%local_field_solve = .false.

      state%init%opt_ov%override_opt_source = .true.
      state%init%opt_ov%opt_source = .false.

      state%init%opt_ov%override_intmom_sub = .true.
      state%init%opt_ov%intmom_sub = .false.
      state%init%opt_ov%override_intspec_sub = .true.
      state%init%opt_ov%intspec_sub = .false.

      state%init%opt_ov%override_field_subgath = .true.
      state%init%opt_ov%field_subgath = .false.
      state%init%opt_ov%override_do_smart_update = .true.
      state%init%opt_ov%do_smart_update = .false.

      state%init%opt_ov%override_field_local_allreduce = .true.
      state%init%opt_ov%field_local_allreduce = .false.
      state%init%opt_ov%override_field_local_allreduce_sub = .true.
      state%init%opt_ov%field_local_allreduce_sub = .false.

      call measure_timestep(state)
      state%init%opt_ov%opt_redist_nbk = .true.
      call measure_timestep(state)
      l1 = state%optim%results%optimal
      state%init%opt_ov%opt_redist_persist = .true.
      call measure_timestep(state)
      l2 = state%optim%results%optimal
      state%init%opt_ov%opt_redist_persist_overlap = .true.
      call measure_timestep(state)
      l3 = state%optim%results%optimal

      ! Here we pick the optimal solution
      if (.not. l3) then
        state%init%opt_ov%opt_redist_persist_overlap = .false.
        if (.not. l2) then 
          state%init%opt_ov%opt_redist_persist_overlap = .false.
          if (.not. l1) then
            state%init%opt_ov%opt_redist_nbk = .false.
          end if
        end if
      end if

      state%init%opt_ov%opt_source = .true.
      call measure_timestep(state)
      l1 = state%optim%results%optimal
      if (.not. l1) then
        state%init%opt_ov%opt_source = .false.
      end if

      state%init%opt_ov%local_field_solve = .true.
      call measure_timestep(state)
      l1 = state%optim%results%optimal
      if (.not. l1) then
        state%init%opt_ov%local_field_solve = .false.
      end if

      !if (state%init%opt_ov%layout .eq. 'xyles' .or. &
          !state%init%opt_ov%layout .eq. 'yxles') then
      !if (.true.) then
      state%init%opt_ov%intmom_sub = .true.
      call measure_timestep(state)
      if (.not. state%optim%results%optimal) then
        state%init%opt_ov%intmom_sub = .false.
      end if
      state%init%opt_ov%intspec_sub = .true.
      call measure_timestep(state)
      l1 = state%optim%results%optimal
      if (state%init%opt_ov%field_option .eq. "local") then
        state%init%opt_ov%field_local_allreduce = .true.
        call measure_timestep(state)
        l2 = state%optim%results%optimal
        state%init%opt_ov%field_local_allreduce_sub = .true.
        call measure_timestep(state)
        l3 = state%optim%results%optimal
        if (.not. l3) then
          state%init%opt_ov%field_local_allreduce_sub = .false.
          if (.not. l2) then 
            state%init%opt_ov%field_local_allreduce = .false.
            if (.not. l1) then
              state%init%opt_ov%intspec_sub = .false.
            end if
          end if
        end if
      else
        if (.not. l1) then
          state%init%opt_ov%intspec_sub = .false.
        end if
      end if
      !end if


      if (state%init%opt_ov%field_option .eq. "implicit") then
        state%init%opt_ov%field_subgath = .true.
        call measure_timestep(state)
        if (.not. state%optim%results%optimal) then
          state%init%opt_ov%field_subgath = .false.
        end if
      else if (state%init%opt_ov%field_option .eq. "local") then
        state%init%opt_ov%do_smart_update = .true.
        call measure_timestep(state)
        if (.not. state%optim%results%optimal) then
          state%init%opt_ov%do_smart_update = .false.
        end if
      end if
    end subroutine optimise_flags

    subroutine optimise_fields(state)
      type(gs2_program_state_type), intent(inout) :: state
      state%init%opt_ov%override_field_option = .true.
      state%init%opt_ov%field_option = "implicit"
      state%init%opt_ov%override_minnrow = .true.
      state%init%opt_ov%minnrow = 64
      call optimise_flags(state)
      state%init%opt_ov%field_option = "local"
      call optimise_flags(state)
      !state%init%opt_ov%minnrow = 16
      !call measure_timestep(state)
      state%init%opt_ov%minnrow = 32
      call optimise_flags(state)
      state%init%opt_ov%minnrow = 128
      call optimise_flags(state)
      !state%init%opt_ov%minnrow = 256
      !call measure_timestep(state)
      !state%init%opt_ov%minnrow = 512
      !call optimise_flags(state)
    end subroutine optimise_fields

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
      integer :: i,n, iresult
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
      iresult = i
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
      if(proc0) call write_summary(6,&
        state%optim%sorted_results(iresult), &
        state%optim%sorted_optimisations(iresult))
        
      deallocate(sorted_opts_temp, sorted_res_temp)
      
    end subroutine measure_timestep

end module gs2_optimisation
