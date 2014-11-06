!> A module for printing summary information to stdout during
!! the simulation
module diagnostics_screen_printout
  use diagnostics_config
  use run_parameters, only: fphi, fapar, fbpar
  use species, only: nspec
  contains
    !> Print out growth rates to screen.
    !! Doesn't quite conform to the principle of having all
    !! useful results stored in gnostics... takes stuff 
    !! straight from diagnostics_omega 
    subroutine print_line(gnostics)
      use mp, only: proc0, sum_reduce
      use kt_grids, only: naky, ntheta0, aky, akx, theta0
      use run_parameters, only: woutunits
      use fields_arrays, only: phinew
      use volume_averages, only: average_theta
      use diagnostics_omega, only: omegahist, omega_average
      implicit none
      type(diagnostics_type), intent(in) :: gnostics
      real, dimension (ntheta0, naky) :: phitot
      integer :: ik, it

      call average_theta(phinew, phinew, phitot, gnostics%distributed)

      if (gnostics%distributed) call sum_reduce(phitot, 0)

      if(.not.(proc0 .and. gnostics%wryte)) return
      do ik = 1, naky
         do it = 1, ntheta0
            write (unit=*, fmt="('ky=',1pe9.2, ' kx=',1pe9.2, &
                 & ' om=',e9.2,1x,e9.2,' omav=',e9.2,1x,e9.2, &
                 & ' phtot=',e9.2,' theta0=',1pe9.2)") &
                 aky(ik), akx(it), &
                 real( omegahist(mod(gnostics%istep,gnostics%navg),it,ik)*woutunits(ik)), &
                 aimag(omegahist(mod(gnostics%istep,gnostics%navg),it,ik)*woutunits(ik)), &
                 real( omega_average(it,ik)*woutunits(ik)), &
                 aimag(omega_average(it,ik)*woutunits(ik)), &
                 phitot(it,ik), theta0(it,ik)
         end do
      end do
      write (*,*) 
    end subroutine print_line
    
    !> Print instaneous heat fluxes to screen
    subroutine print_flux_line(gnostics)
      type(diagnostics_type), intent(in) :: gnostics
      if (fphi > epsilon(0.0)) then
         write (unit=*, fmt="('t= ',e17.10,' <phi**2>= ',e13.6, &
              & ' heat fluxes: ', 5(1x,e13.6))") &
              gnostics%user_time, &
              gnostics%current_results%phi2, &
              gnostics%current_results%species_es_heat_flux(1:min(nspec,5))
         write (unit=*, fmt="('t= ',e17.10,' <phi**2>= ',e13.6, &
              & ' energy exchange: ', 5(1x,e13.6))") &
              gnostics%user_time, &
              gnostics%current_results%phi2, &
              gnostics%current_results%species_energy_exchange(1:min(nspec,5))
      end if
      if (fapar > epsilon(0.0)) then
         write (unit=*, fmt="('t= ',e17.10,' <apar**2>= ',e11.4, &
              & ' heat flux m: ', 5(1x,e11.4))") &
              gnostics%user_time, &
              gnostics%current_results%apar2, &
              gnostics%current_results%species_apar_heat_flux(1:min(nspec,5))
      end if
      if (fbpar > epsilon(0.0)) then
         write (unit=*, fmt="('t= ',e17.10,' <bpar**2>= ',e11.4, &
              & ' heat flux b: ', 5(1x,e11.4))") &
              gnostics%user_time, &
              gnostics%current_results%bpar2, & 
              gnostics%current_results%species_bpar_heat_flux(1:min(nspec,5))
      end if
    end subroutine print_flux_line

end module diagnostics_screen_printout
