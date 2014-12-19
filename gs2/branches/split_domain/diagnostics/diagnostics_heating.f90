!> Subroutines for writing out the turbulent heating...
!! the turbulent energy dissipated in collisions and hyperviscosity
module diagnostics_heating
  use gs2_heating, only: heating_diagnostics
  implicit none
  private

  public :: init_diagnostics_heating, finish_diagnostics_heating
  public :: calculate_heating, write_heating

  type (heating_diagnostics) :: h
  type (heating_diagnostics), dimension(:), allocatable :: h_hist
  type (heating_diagnostics), dimension(:,:), allocatable :: hk
  type (heating_diagnostics), dimension(:,:,:), allocatable :: hk_hist
contains
    
  subroutine finish_diagnostics_heating(gnostics)
    use gs2_heating, only: del_htype
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent (in) :: gnostics

    if (gnostics%write_heating) then
       call del_htype (h)
       call del_htype (h_hist)
       call del_htype (hk_hist)
       call del_htype (hk)
    endif
    
    if (allocated(h_hist)) deallocate (h_hist, hk_hist, hk)
  end subroutine finish_diagnostics_heating

  subroutine init_diagnostics_heating(gnostics)
    use gs2_heating, only: init_htype
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent (in) :: gnostics
    integer :: navg
    navg = gnostics%navg
    
    ! allocate heating diagnostic data structures
    if (gnostics%write_heating) then
       allocate (h_hist(0:navg-1))
       call init_htype (h_hist,  nspec)
       
       allocate (hk_hist(ntheta0,naky,0:navg-1))
       call init_htype (hk_hist, nspec)
       
       call init_htype (h,  nspec)

       allocate (hk(ntheta0, naky))
       call init_htype (hk, nspec)
    else
       allocate (h_hist(0))
       allocate (hk(1,1))
       allocate (hk_hist(1,1,0))
    end if
  end subroutine init_diagnostics_heating
  
  subroutine calculate_heating (gnostics)
    use mp, only: proc0
    use dist_fn, only: get_heat
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use species, only: nspec, spec
    use kt_grids, only: naky, ntheta0, aky, akx
    use theta_grid, only: ntgrid, delthet, jacob
    use dist_fn_arrays, only: c_rate
    use gs2_heating, only: heating_diagnostics, avg_h, avg_hk, zero_htype
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent (in) :: gnostics
    real, dimension(-ntgrid:ntgrid) :: wgt
    real :: fac
    integer :: is, ik, it, ig
    
    !Zero out variables for heating diagnostics
    call zero_htype(h)
    call zero_htype(hk)
    
    if (proc0) then       
       !GGH NOTE: Here wgt is 1/(2*ntgrid+1)
       wgt = delthet*jacob
       wgt = wgt/sum(wgt)
       
       do is = 1, nspec
          do ik = 1, naky
             fac = 0.5
             if (aky(ik) < epsilon(0.)) fac = 1.0
             do it = 1, ntheta0
                if (aky(ik) < epsilon(0.0) .and. abs(akx(it)) < epsilon(0.0)) cycle
                do ig = -ntgrid, ntgrid
                   
                   !Sum heating by k over all z points (ig)
                   hk(it, ik) % collisions(is) = hk(it, ik) % collisions(is) &
                        + real(c_rate(ig,it,ik,is,1))*fac*wgt(ig)*spec(is)%temp*spec(is)%dens
                   
                   hk(it, ik) % hypercoll(is) = hk(it, ik) % hypercoll(is) &
                        + real(c_rate(ig,it,ik,is,2))*fac*wgt(ig)*spec(is)%temp*spec(is)%dens
                   
                   hk(it, ik) % imp_colls(is) = hk(it, ik) % imp_colls(is) &
                        + real(c_rate(ig,it,ik,is,3))*fac*wgt(ig)*spec(is)%temp*spec(is)%dens
                   
                end do
                h % collisions(is) = h % collisions(is) + hk(it, ik) % collisions(is)
                h % hypercoll(is)  = h % hypercoll(is)  + hk(it, ik) % hypercoll(is)
                h % imp_colls(is)  = h % imp_colls(is)  + hk(it, ik) % imp_colls(is)
             end do
          end do
       end do
    end if
    
    call get_heat (h, hk, phi, apar, bpar, phinew, aparnew, bparnew)    
    
    call avg_h(h, h_hist, gnostics%istep, gnostics%navg)
    call avg_hk(hk, hk_hist, gnostics%istep, gnostics%navg)
  end subroutine calculate_heating

  subroutine write_heating(gnostics)
    use species, only: nspec
    use file_utils, only: flush_output_file
    use mp, only: proc0
    use diagnostics_config, only: diagnostics_type
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_dimensions, only: dim_string
    implicit none
    type(diagnostics_type), intent (inout) :: gnostics
    real :: t
    integer :: is
    
    t = gnostics%user_time
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_energy", &
         dim_string([gnostics%dims%time]), &
         "Total free (turbulent) energy", "F_r^2?? TBC", h%energy)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_energy_dot", &
         dim_string([gnostics%dims%time]), &
         "d/dt of total free (turbulent) energy", "F_r^2 a/v_thr?? TBC", h%energy_dot)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_antenna", &
         dim_string([gnostics%dims%time]), &
         "Energy injection resulting from antenna J_ant.E", "TBC", h%antenna)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_imp_colls", &
         dim_string([gnostics%dims%species,gnostics%dims%time]), &
         "?? [h_(i+1)*h_*]/2 * C[h_(i+1)] * T_0", "TBC", h%imp_colls)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_imp_colls_sum", &
         dim_string([gnostics%dims%time]), &
         "?? [h_(i+1)*h_*]/2 * C[h_(i+1)] * T_0", "TBC", sum(h%imp_colls))
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_hypercoll", &
         dim_string([gnostics%dims%species,gnostics%dims%time]), &
         "Heating resulting from hyperviscosity -[h H(h) * T_0]", "TBC", h%hypercoll)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_hypercoll_sum", &
         dim_string([gnostics%dims%time]), &
         "Heating resulting from hyperviscosity -[h H(h) * T_0]", "TBC", sum(h%hypercoll))
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_collisions", &
         dim_string([gnostics%dims%species,gnostics%dims%time]), &
         "Heating resulting from collisions -[h C(h) * T_0]", "TBC", h%collisions)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_collisions_sum", &
         dim_string([gnostics%dims%time]), &
         "Heating resulting from collisions -[h C(h) * T_0]", "TBC", sum(h%collisions))
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_gradients", &
         dim_string([gnostics%dims%species,gnostics%dims%time]), &
         "Energy injection resulting from gradients [h omega_* h]", "TBC", h%gradients)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_gradients_sum", &
         dim_string([gnostics%dims%time]), &
         "Energy injection resulting from gradients [h omega_* h]", "TBC", sum(h%gradients))
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_heating", &
         dim_string([gnostics%dims%species,gnostics%dims%time]), &
         "Total heating sum [h (q dchi/dt - dh/dt * T0)]", "TBC", h%heating)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_heating_sum", &
         dim_string([gnostics%dims%time]), &
         "Total heating sum [h (q dchi/dt - dh/dt * T0)]", "TBC", sum(h%heating))
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_energy_balance", &
         dim_string([gnostics%dims%time]), &
         "Sum of total heating, total injection and total change in stored &
         & free energy (should be zero)", "TBC", sum(h%heating)+h%antenna+sum(h%gradients)+h%energy_dot)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_energy_apar", &
         dim_string([gnostics%dims%time]), &
         "Free energy in perp magnetic field: (k_perp A)**2", "TBC", h%eapar)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_energy_bpar", &
         dim_string([gnostics%dims%time]), &
         "Free energy in par magnetic field: B_par**2", "TBC", h%ebpar)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_energy_delfs2", &
         dim_string([gnostics%dims%species,gnostics%dims%time]), &
         "Total free energy as a function of species dfs^2  ", "TBC", h%delfs2)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_energy_hs2", &
         dim_string([gnostics%dims%species,gnostics%dims%time]), &
         "Non-adiabatic free energy as a function of species hs^2  ", "TBC", h%hs2)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_energy_phis2", &
         dim_string([gnostics%dims%species,gnostics%dims%time]), &
         "Adiabatic free energy as a function of species phi_bar^2  ", "TBC", h%phis2)
    
    gnostics%current_results%species_heating = h%imp_colls
    gnostics%current_results%species_heating_avg = gnostics%current_results%species_heating_avg + &
         h%imp_colls*(gnostics%user_time - gnostics%user_time_old)
    
    if (gnostics%write_ascii .and. proc0) call write_ascii
    
  contains
    subroutine write_ascii
      integer :: heat_unit
      integer :: heat_unit2
      
      heat_unit = gnostics%ascii_files%heat
      heat_unit2 = gnostics%ascii_files%heat2
      !
      !
      ! For case with one species:
      !
      ! Column     Item               
      !   1        time              
      !   2        Energy              
      !   3        dEnergy/dt            
      !   4        J_ant.E             
      !   5        [h_(i+1)*h_*]/2 * C[h_(i+1)] * T_0 
      !   6       -[h H(h) * T_0]
      !   7       -[h C(h) * T_0]
      !   8        [h w_* h]
      !   9        [h * (q dchi/dt - dh/dt * T0)]_1
      !  10      sum (h C(h) * T_0)  in total, as in 5, 6      
      !  11     -sum (h H(h) * T_0)      
      !  12     -sum (h C(h) * T_0)   
      !  13      sum (h w_* h)  
      !  14      sum [h (q dchi/dt - dh/dt * T0)]
      !  15      3 + 4 + 9 + 10
      !  16      (k_perp A)**2
      !  17      B_par**2
      !  18      df ** 2
      !  19      h ** 2
      !  20      Phi_bar ** 2
      
      
      write (unit=heat_unit, fmt="(28es12.4)") t,h % energy,  &
           h % energy_dot, h % antenna, h % imp_colls, h % hypercoll, h % collisions, &
           h % gradients, h % heating, sum(h % imp_colls), sum(h % hypercoll), sum(h % collisions), &
           sum(h % gradients), sum(h % heating),sum(h%heating)+h%antenna+sum(h%gradients)+h%energy_dot, &
           h % eapar, h % ebpar, h % delfs2(:),  h % hs2(:), h % phis2(:)
      
      do is=1,nspec
         write (unit=heat_unit2, fmt="(15es12.4)") t,h % energy,  &
              h % energy_dot, h % antenna, h % imp_colls(is), h % hypercoll(is), h % collisions(is), &
              h % gradients(is), h % heating(is), &
              h % eapar, h % ebpar, h % delfs2(is),  h % hs2(is), h % phis2(is), real(is)
         write (unit=heat_unit2, fmt=*)
      end do
      write (unit=heat_unit2, fmt=*)
      
      call flush_output_file (heat_unit)
      call flush_output_file (heat_unit2)
      
      !!GGH          write (unit=heat_unit, fmt="('t= ',e13.6,' energy= ',e13.6)") t, h % energy
      !!GGH          write (unit=heat_unit, fmt="('t= ',e13.6,' energy_dot= ',e13.6)") t, h % energy_dot
      !!GGH          write (unit=heat_unit, fmt="('t= ',e13.6,' J_ant.E= ',e13.6)") t, h % antenna
      
      !!GGH          write (unit=heat_unit, fmt="('t= ',e13.6,' hyperC= ',12(1x,e13.6))") t, h % hypercoll
      !!GGH          write (unit=heat_unit, fmt="('t= ',e13.6,' hCh= ',12(1x,e13.6))") t, h % collisions
      !!GGH          write (unit=heat_unit, fmt="('t= ',e13.6,' hw*= ',12(1x,e13.6))") t, h % gradients
      !!GGH!         write (unit=heat_unit, fmt="('t= ',e13.6,' hwd= ',12(1x,e13.6))") t, h % curvature
      !!GGH          write (unit=heat_unit, fmt="('t= ',e13.6,' heating= ',12(1x,e13.6))") t, h % heating
      
      !!GGH          write (unit=heat_unit, fmt="('t= ',e13.6,' total_hvisc= ',e13.6)") t, sum(h % hypervisc)
      !!GGH          write (unit=heat_unit, fmt="('t= ',e13.6,' total_hyperC= ',e13.6)") t, sum(h % hypercoll)
      !!GGH          write (unit=heat_unit, fmt="('t= ',e13.6,' total_hCh= ',e13.6)") t, sum(h % collisions)
      !!GGH          write (unit=heat_unit, fmt="('t= ',e13.6,' total_hw*= ',e13.6)") t, sum(h % gradients)
      !!GGH          write (unit=heat_unit, fmt="('t= ',e13.6,' total_heating= ',e13.6)") t, sum(h % heating)
      
      !!GGH          write (unit=heat_unit, fmt="('t= ',e13.6,' total_power= ',e13.6)") t, &
      !!GGH               sum(h%heating)+h%antenna+sum(h%gradients)+h%energy_dot
      !!GGH TEST try adding sqrt(2.) to the edot
      !!GGH          write (unit=heat_unit, fmt="('t= ',e13.6,' total_power= ',e13.6)") t, &
      !!GGH               sum(h%heating)+h%antenna+sum(h%gradients)+h%energy_dot*sqrt(2.)
      !!GGH          write (unit=heat_unit, fmt='(a)') ''
      !!end if
    end subroutine write_ascii
  end subroutine write_heating
end module diagnostics_heating
