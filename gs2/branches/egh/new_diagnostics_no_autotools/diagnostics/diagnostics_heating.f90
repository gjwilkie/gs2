!> Subroutines for writing out the turbulent heating...
!! the turbulent energy dissipated in collisions and hyperviscosity
module diagnostics_heating
  use diagnostics_config
  use diagnostics_create_and_write, only: create_and_write_variable
  use gs2_heating, only: heating_diagnostics
  use species, only: nspec
  use kt_grids, only: naky, ntheta0
  type (heating_diagnostics), save :: h
  type (heating_diagnostics), dimension(:), save, allocatable :: h_hist
  type (heating_diagnostics), dimension(:,:), save, allocatable :: hk
  type (heating_diagnostics), dimension(:,:,:), save, allocatable :: hk_hist
  contains

  subroutine finish_diagnostics_heating(gnostics)
    use gs2_heating, only: del_htype
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
    implicit none
    type(diagnostics_type), intent (in) :: gnostics
    !type (heating_diagnostics) :: h
    !type (heating_diagnostics), dimension(:,:) :: hk

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
    implicit none
    type(diagnostics_type), intent (in) :: gnostics
    real :: t
    integer :: is

    t = gnostics%user_time
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_energy", "t", &
      "Total free (turbulent) energy", "F_r^2?? TBC", h%energy)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_energy_dot", "t", &
      "d/dt of total free (turbulent) energy", "F_r^2 a/v_thr?? TBC", h%energy_dot)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_antenna", "t", &
      "Energy injection resulting from antenna J_ant.E", "TBC", h%antenna)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_imp_colls", "st", &
      "?? [h_(i+1)*h_*]/2 * C[h_(i+1)] * T_0", "TBC", h%imp_colls)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_imp_colls_sum", "t", &
      "?? [h_(i+1)*h_*]/2 * C[h_(i+1)] * T_0", "TBC", sum(h%imp_colls))
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_hypercoll", "st", &
      "Heating resulting from hyperviscosity -[h H(h) * T_0]", "TBC", h%hypercoll)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_hypercoll_sum", "t", &
      "Heating resulting from hyperviscosity -[h H(h) * T_0]", "TBC", sum(h%hypercoll))
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_collisions", "st", &
      "Heating resulting from collisions -[h C(h) * T_0]", "TBC", h%collisions)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_collisions_sum", "t", &
      "Heating resulting from collisions -[h C(h) * T_0]", "TBC", sum(h%collisions))
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_gradients", "st", &
      "Energy injection resulting from gradients [h omega_* h]", "TBC", h%gradients)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_gradients_sum", "t", &
      "Energy injection resulting from gradients [h omega_* h]", "TBC", sum(h%gradients))
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_heating", "st", &
      "Total heating sum [h (q dchi/dt - dh/dt * T0)]", "TBC", h%heating)
    call create_and_write_variable(gnostics, gnostics%rtype, "heating_heating_sum", "t", &
      "Total heating sum [h (q dchi/dt - dh/dt * T0)]", "TBC", sum(h%heating))

            !h % energy_dot, h % antenna, h % imp_colls, h % hypercoll, h % collisions, &
            !h % gradients, h % heating, sum(h % imp_colls), sum(h % hypercoll), sum(h % collisions), &
            !sum(h % gradients), sum(h % heating),sum(h%heating)+h%antenna+sum(h%gradients)+h%energy_dot, &
            !h % eapar, h % ebpar, h % delfs2(:),  h % hs2(:), h % phis2(:)

       !  13        [h * (q dchi/dt - dh/dt * T0)]_1
       !  14        [h * (q dchi/dt - dh/dt * T0)]_2
       !  15      sum (h C(h) * T_0)  in total, as in 5, 6      
       !  16     -sum (h H(h) * T_0)      
       !  17     -sum (h C(h) * T_0)   
       !  18      sum (h w_* h)  
       !  19      sum [h (q dchi/dt - dh/dt * T0)]
       !  20      3 + 4 + 18 + 19
       !  21      (k_perp A)**2
       !  22      B_par**2
       !  23      df_1 ** 2
       !  24      df_2 ** 2
       !  25      h_1 ** 2
       !  26      h_2 ** 2
       !  27      Phi_bar_1 ** 2
       !  28      Phi_bar_2 ** 2
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
       

       !write (unit=heat_unit, fmt="(28es12.4)") t,h % energy,  &
            !h % energy_dot, h % antenna, h % imp_colls, h % hypercoll, h % collisions, &
            !h % gradients, h % heating, sum(h % imp_colls), sum(h % hypercoll), sum(h % collisions), &
            !sum(h % gradients), sum(h % heating),sum(h%heating)+h%antenna+sum(h%gradients)+h%energy_dot, &
            !h % eapar, h % ebpar, h % delfs2(:),  h % hs2(:), h % phis2(:)

       !do is=1,nspec
          !write (unit=heat_unit2, fmt="(15es12.4)") t,h % energy,  &
               !h % energy_dot, h % antenna, h % imp_colls(is), h % hypercoll(is), h % collisions(is), &
               !h % gradients(is), h % heating(is), &
               !h % eapar, h % ebpar, h % delfs2(is),  h % hs2(is), h % phis2(is), real(is)
          !write (unit=heat_unit2, fmt=*)
       !end do
       !write (unit=heat_unit2, fmt=*)

       !call flush_output_file (heat_unit)
       !call flush_output_file (heat_unit2)

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
  end subroutine write_heating

end module diagnostics_heating
