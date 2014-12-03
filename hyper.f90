! nexp has been changed in a rush.  Only known to be correct for nexp=2
!
!
module hyper
  
  implicit none

  private

  public :: init_hyper, finish_hyper, hyper_diff, D_res
  public :: read_parameters, wnml_hyper, check_hyper
  public :: D_v, D_eta, nexp
  public :: hypervisc_filter

! D_v, D_eta are hyper coefficients, normalized correctly 
! i.e., either by unity or by 1/k_perp**2*nexp

  real :: D_v, D_eta
  real :: D_hypervisc, D_hyperres, omega_osc, D_hyper
  real :: akx4_max, aky4_max, aky_max, akperp4_max

  integer :: hyper_option_switch, nexp
  integer, parameter :: hyper_option_none = 1, &
       hyper_option_visc = 2, &
       hyper_option_res  = 3, &
       hyper_option_both = 4
  character(9) :: hyper_option
  logical :: const_amp, include_kpar, isotropic_shear, damp_zonal_only
  logical :: hyper_on = .false.
  logical :: gridnorm

  real, dimension (:,:), allocatable :: D_res
  ! (it, ik)

  real, dimension (:,:,:), allocatable :: hypervisc_filter
  ! (-ntgrid:ntgrid, ntheta0, naky)

  logical :: initialized = .false.

contains

  subroutine check_hyper(report_unit)
    implicit none
    integer, intent(in) :: report_unit
    select case (hyper_option_switch)
    case (hyper_option_none)
       if (D_hyperres > 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses no hyperresistivity.  &
               &D_hyperres ignored.')") trim(hyper_option)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          D_hyperres = -10.
       end if
       if (D_hypervisc > 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses no hyperviscosity.  &
               &D_hypervisc ignored.')") trim(hyper_option)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          D_hypervisc = -10.
       endif

    case (hyper_option_visc)
       hyper_on = .true.
       if (D_hyperres > 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses no hyperresistivity.  &
               &D_hyperres ignored.')") trim(hyper_option)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          D_hyperres = -10.
       end if
       if (D_hypervisc < 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses hyperviscosity but &
               &D_hypervisc < 0.')") trim(hyper_option)
          write (report_unit, fmt="('No hyperviscosity used.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          hyper_on = .false.
       endif

    case (hyper_option_res)
       hyper_on = .true.
       if (D_hyperres < 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses hyperresistivity but D_hyperres < 0.')") trim(hyper_option)
          write(report_unit, fmt="('No hyperresistivity used.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          hyper_on = .false.
       end if
       if (D_hypervisc > 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses no hyperviscosity.  D_hypervisc ignored.')") trim(hyper_option)
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
          D_hypervisc = -10.
       endif

    case (hyper_option_both)
       hyper_on = .true.
       if (D_hyperres < 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses hyperresistivity but D_hyperres < 0.')") trim(hyper_option)
          write (report_unit, fmt="('No hyperresistivity used.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if
       if (D_hypervisc < 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('hyper_option = ',a,' chooses hyperviscosity but D_hypervisc < 0.')") trim(hyper_option)
          write (report_unit, fmt="('No hyperviscosity used.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       endif
       if (D_hypervisc < 0. .and. D_hyperres < 0.) hyper_on = .false.

    end select

    if (hyper_on) then
       write (report_unit, *) 
       write (report_unit, fmt="('------------------------------------------------------------')")
       write (report_unit, *) 

       select case (hyper_option_switch)

       case (hyper_option_visc)

          write (report_unit, *) 
          write (report_unit, fmt="('Hyperviscosity included without hyperresistivity.')")
          if (const_amp) then
             write (report_unit, fmt="('Damping rate is ',e11.4,' at highest k_perp.')") D_hypervisc
          else
             write (report_unit, fmt="('The damping coefficent is ',e11.4)") D_hypervisc
             write (report_unit, fmt="('The damping rate is proportional to the RMS amplitude of the turbulence.')")
          end if
          if (isotropic_shear) then
             write (report_unit, fmt="('The hyperviscosity is isotropic in the perpendicular plane.')")
             write (report_unit, fmt="('This is appropriate for MHD-like calculations.')")
          else
             write (report_unit, fmt="('The hyperviscosity is anisotropic in the perpendicular plane.')")
             write (report_unit, fmt="('This is appropriate for drift-type calculations.')")
             write (report_unit, fmt="('omega_osc = ',e11.4)") omega_osc
          end if

       case (hyper_option_res)

          write (report_unit, *) 
          write (report_unit, fmt="('Hyperresistivity included without hyperviscosity.')")
          if (const_amp) then
             write (report_unit, fmt="('Damping rate is ',e11.4,' at highest k_perp.')") D_hyperres
          else
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('const_amp = .false. is not implemented for hyperresistivity.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 
          end if
          if (isotropic_shear) then
             write (report_unit, fmt="('The hyperresistivity is isotropic in the perpendicular plane.')")
             write (report_unit, fmt="('This is appropriate for MHD-like calculations.')")
          else
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('isotropic_shear = .false. is not implemented for hyperresistivity.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 
          end if

       case (hyper_option_both)

          write (report_unit, *) 
          write (report_unit, fmt="('Hyperresistivity and hyperviscosity included.')")
          if (const_amp) then
             write (report_unit, fmt="('Damping rate is ',e11.4,' at highest k_perp.')") D_hyperres
          else
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('const_amp = .false. is not implemented for hyperresistivity.')")
             write (report_unit, fmt="('THIS IS AN ERROR.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 
          end if
          if (isotropic_shear) then
             write (report_unit, fmt="('The damping is isotropic in the perpendicular plane.')")
             write (report_unit, fmt="('This is appropriate for MHD-like calculations.')")
          else
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('isotropic_shear = .false. is not implemented for hyperresistivity.')")
             write (report_unit, fmt="('THIS IS AN ERROR.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 
          end if
       end select
    end if
  end subroutine check_hyper

  subroutine wnml_hyper(unit)
    implicit none
    integer, intent(in) :: unit          
    if (.not. hyper_on) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "hyper_knobs"

    select case (hyper_option_switch)

    case (hyper_option_visc) 
       write (unit, fmt="(' hyper_option = ',a)") '"visc_only"'
       write (unit, fmt="(' D_hypervisc = ',e17.10)") D_hypervisc

    case (hyper_option_res) 
       write (unit, fmt="(' hyper_option = ',a)") '"res_only"'
       write (unit, fmt="(' D_hyperres = ',e17.10)") D_hyperres

    case (hyper_option_both) 
       write (unit, fmt="(' hyper_option = ',a)") '"both"'
       if (D_hyperres == D_hypervisc) then
          write (unit, fmt="(' D_hyper = ',e17.10)") D_hyper
       else
          write (unit, fmt="(' D_hypervisc = ',e17.10)") D_hypervisc
          write (unit, fmt="(' D_hyperres = ',e17.10)") D_hyperres
       end if
    end select

!    write (unit, fmt="(' include_kpar = ',L1)") include_kpar
    
    write (unit, fmt="(' const_amp = ',L1)") const_amp
    write (unit, fmt="(' isotropic_shear = ',L1)") isotropic_shear
    if (.not. isotropic_shear) &
         write (unit, fmt="(' omega_osc = ',e17.10)") omega_osc

    write (unit, fmt="(' gridnorm = ',L1)") gridnorm
    write (unit, fmt="(' /')")
  end subroutine wnml_hyper

  subroutine init_hyper
    use kt_grids, only: ntheta0, naky, akx, aky
    use gs2_time, only: code_dt
    use gs2_layouts, only: init_gs2_layouts
    implicit none
    integer :: ik, it

    if (initialized) return
    initialized = .true.
    
    call init_gs2_layouts
    call read_parameters
    call allocate_arrays

! Define variables used in hyperviscosity and hyperresistivity models

    if (gridnorm) then
       akx4_max    = akx(ntheta0/2 + 1) ** (2*nexp)
       aky_max     = aky(naky)
       aky4_max     = aky(naky) ** (2*nexp)
       akperp4_max = ( akx(ntheta0/2 + 1) ** 2  +  aky(naky) ** 2) ** (nexp)
    else
       akx4_max = 1.
       aky_max = 1.
       aky4_max = 1.
       akperp4_max = 1.
    end if

! Get D_res set up if needed

    if (D_hyperres > 0.) then
       do ik = 1, naky
          do it = 1, ntheta0
             D_res(it, ik) = D_hyperres*code_dt &
                  * (aky(ik)**2 + akx(it)**2)**nexp/akperp4_max
          end do
       end do
       D_eta = D_hyperres/akperp4_max
    else
       D_res = 0.
       D_eta = 0.
    end if

    if (D_hypervisc > 0.) then
       D_v = D_hypervisc/akperp4_max
    else
       D_v = 0.
    end if

  end subroutine init_hyper

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    implicit none
    type (text_option), dimension(5), parameter :: hyperopts = &
         (/ text_option('default', hyper_option_none), &
            text_option('none', hyper_option_none), &
            text_option('visc_only', hyper_option_visc), &
            text_option('res_only', hyper_option_res), &
            text_option('both', hyper_option_both) /)
    integer :: ierr, in_file
    logical :: exist
    
    namelist /hyper_knobs/ hyper_option, const_amp, include_kpar, &
         isotropic_shear, D_hyperres, D_hypervisc, D_hyper, omega_osc, gridnorm,&
         nexp, damp_zonal_only
    
    if (proc0) then
       const_amp = .false.
       include_kpar = .false.
       isotropic_shear = .true.
       damp_zonal_only = .false.
       nexp = 2
       D_hyperres = -10.
       D_hypervisc = -10.
       D_hyper = -10.
       hyper_option = 'default'
       omega_osc = 0.4
       gridnorm = .true.
       in_file=input_unit_exist("hyper_knobs",exist)
       if (exist) read (unit=input_unit("hyper_knobs"), nml=hyper_knobs)

       ierr = error_unit()

       call get_option_value &
            (hyper_option, hyperopts, hyper_option_switch, &
            ierr, "hyper_option in hyper_knobs",.true.)

       if (.not. isotropic_shear .and. nexp /=2) then
          write (ierr, *) 'Forcing nexp = 2.  Higher values not implemented for anisotropic shear model.'
          nexp = 2
       end if

       select case (hyper_option_switch)

       case (hyper_option_none)
          if (D_hyperres > 0.) then
             write(ierr, *) 'hyper_option = ',hyper_option, &
                  ' chooses no hyperresistivity.  D_hyperres ignored.'
             D_hyperres = -10.
          end if
          if (D_hypervisc > 0.) then
             write(ierr, *) 'hyper_option = ',hyper_option, &
                  ' chooses no hyperviscosity.  D_hypervisc ignored.'
             D_hypervisc = -10.
          endif

       case (hyper_option_visc)
          hyper_on = .true.
          write (ierr, *) 'WARNING: It is inconsistent to set D_hypervisc different from ', &
               'D_hyperres.  Recommend: Set them equal.'
          if (D_hyperres > 0.) then
             write(ierr, *) 'hyper_option = ',hyper_option, &
                  ' chooses no hyperresistivity.  D_hyperres ignored.'
             D_hyperres = -10.
          end if
          if (D_hypervisc < 0.) then
             write(ierr, *) 'hyper_option = ',hyper_option, &
                  ' chooses hyperviscosity but D_hypervisc < 0 is illegal.'
             write(ierr, *) 'No hyperviscosity used.'
             hyper_on = .false.
          endif

       case (hyper_option_res)
          hyper_on = .true.
          write (ierr, *) 'WARNING: It is inconsistent to set D_hypervisc different from ', &
               'D_hyperres.  Recommend: Set them equal.'
          if (D_hyperres < 0.) then
             write(ierr, *) 'hyper_option = ',hyper_option, &
                  ' chooses hyperresistivity but  D_hyperres < 0 is illegal.'
             write(ierr, *) 'No hyperresistivity used.'
             hyper_on = .false.
          end if
          if (D_hypervisc > 0.) then
             write(ierr, *) 'hyper_option = ',hyper_option, &
                  ' chooses no hyperviscosity.  D_hypervisc ignored.'
             D_hypervisc = -10.
          endif

       case (hyper_option_both)
          hyper_on = .true.
          
          if (D_hyper < 0.) then
             if (D_hyperres /= D_hyperres) then
                write (ierr, *) 'WARNING: It is inconsistent to set D_hypervisc different from ', &
                     'D_hyperres.  Recommend: Set them equal.'
             end if
          else
             write (ierr, *) 'WARNING: Setting D_hypervisc = D_hyperres, each to value of D_hyper'
             D_hyperres  = D_hyper
             D_hypervisc = D_hyper
          end if

          if (D_hyperres < 0.) then
             write(ierr, *) 'hyper_option = ',hyper_option, &
                  ' chooses hyperresistivity but  D_hyperres < 0 is illegal.'
             write(ierr, *) 'No hyperresistivity used.'
          end if
          if (D_hypervisc < 0.) then
             write(ierr, *) 'hyper_option = ',hyper_option, &
                  ' chooses hyperviscosity but D_hypervisc < 0 is illegal.'
             write(ierr, *) 'No hyperviscosity used.'
          endif
          if (D_hypervisc < 0. .and. D_hyperres < 0.) hyper_on = .false.

       end select
    end if
    
    call broadcast (const_amp)
    call broadcast (include_kpar)
    call broadcast (isotropic_shear)
    call broadcast (damp_zonal_only)
    call broadcast (D_hyper)
    call broadcast (D_hyperres)
    call broadcast (D_hypervisc)
    call broadcast (hyper_option_switch)
    call broadcast (hyper_on)
    call broadcast (omega_osc)
    call broadcast (gridnorm)
    call broadcast (nexp)

  end subroutine read_parameters

  subroutine allocate_arrays
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    implicit none

    if (.not. allocated(D_res)) then
       allocate (D_res(ntheta0, naky)) 
    endif
    if (.not. allocated(hypervisc_filter)) then
       allocate (hypervisc_filter(-ntgrid:ntgrid,ntheta0,naky)) ; hypervisc_filter = 1.0
    end if
    D_res = 0.

  end subroutine allocate_arrays

  subroutine hyper_diff (g0, phi)

    use gs2_layouts, only: ik_idx, it_idx, is_idx
    use theta_grid, only: ntgrid
    use gs2_time, only: code_dt
    use gs2_layouts, only: g_lo
    use kt_grids, only: aky, akx, naky, ntheta0

    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g0
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi

    real, dimension (-ntgrid:ntgrid) :: shear_rate_nz, shear_rate_z, shear_rate_z_nz

    integer :: iglo, ik, it
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Includes models by Belli and Hammett
! to calculate the x-y avged shearing rate
! S(theta)^2 = <|grad_perp|^4 |phi|^2> 
!            =  sum_over_k(kperp^4 * |phi|^2)
!               (times crazy fac factor due to FFT conventions.)
!
! and to implement this anisotropically in k_perp, taking into 
! account properties of zonal flows.
!
! Begun December 2001
!
! Number conservation added April 2006
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not. hyper_on) return
    if (D_hypervisc < 0.) return

    if(isotropic_shear) then
       call iso_shear
    else
       call aniso_shear
    end if

  contains

    subroutine aniso_shear 

      real, dimension(naky) :: fac
      
! Do the Belli-Hammett anisotropic calculation 
! which accounts for some zonal/non-zonal differences
       
      fac = 0.5
      fac(1) = 1.0
         
! shearing rate due to non-zonal modes (on nonzonal modes)
      shear_rate_nz = 0.
      do ik = 2, naky
         do it = 1, ntheta0
            shear_rate_nz(:) = shear_rate_nz(:) + real(conjg(phi(:,it,ik)) * &
                 phi(:,it,ik)) * (akx(it)**2 + aky(ik)**2)**2 * fac(ik)
         end do
      end do
      shear_rate_nz = 0.5 * ( -omega_osc + (omega_osc ** 2 + 2 * shear_rate_nz) ** 0.5 )
       
! shearing rate due to zonal modes (on nonzonal modes)
      shear_rate_z = 0.
      do ik = 1, 1
         do it = 1, ntheta0
            shear_rate_z(:) = shear_rate_z(:) + real(conjg(phi(:,it,ik)) * &
                 phi(:,it,ik)) * (akx(it)**2 + aky(ik)**2)**2 * fac(ik)
         end do
      end do
! shear_rate_z = shear_rate_z ** 0.5
      shear_rate_z = 0.5 * ( -omega_osc + (omega_osc ** 2 + 2 * shear_rate_z) ** 0.5 )
      
! shearing rate due to nonzonal modes (on zonal modes)
      shear_rate_z_nz = 0.
      do ik = 2, naky
         do it = 1, ntheta0
            shear_rate_z_nz(:) = shear_rate_z_nz(:) + real(conjg(phi(:,it,ik)) * &
                 phi(:,it,ik)) *  aky(ik)**4 * fac(ik)
         end do
      end do
! shear_rate_z_nz = shear_rate_z_nz ** 0.5
      shear_rate_z_nz = 0.5 * ( -omega_osc + (omega_osc ** 2 + 2 * shear_rate_z_nz) ** 0.5 )
       
! end of anisotropic shearing rate calculations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         ik = ik_idx(g_lo, iglo)
         it = it_idx(g_lo, iglo)
         
         if(aky(ik) == 0.) then
            hypervisc_filter(:,it,ik) = exp(- (D_hypervisc * code_dt &
                 * ( shear_rate_z_nz(:) * akx(it) ** 4 / akx4_max )))
         else
            hypervisc_filter(:,it,ik) = exp(- ( D_hypervisc * code_dt & 
                 * ( shear_rate_nz(:) *  (aky(ik) ** 2 + akx(it) ** 2 )**nexp / akperp4_max & 
                 + shear_rate_z(:) * akx(it) ** 4/ akx4_max * aky(ik) / aky_max )))
         endif
         
         g0(:,1,iglo) = g0(:,1,iglo) * hypervisc_filter(:,it,ik)
         g0(:,2,iglo) = g0(:,2,iglo) * hypervisc_filter(:,it,ik)
      end do
    
    end subroutine aniso_shear

    subroutine iso_shear

      real, dimension(naky) :: fac
      
      if (const_amp) then
         shear_rate_nz = 1.
      else
         fac = 0.5
         fac(1) = 1.0
         shear_rate_nz = 0.
         do ik = 1, naky
            do it = 1, ntheta0
               shear_rate_nz(:) = shear_rate_nz(:) &
                    + real(conjg(phi(:,it,ik))*phi(:,it,ik)) &
                    * (akx(it)**2 + aky(ik)**2)**2 * fac(ik)
            end do
         end do
         shear_rate_nz = shear_rate_nz**0.5
      end if

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         ik = ik_idx(g_lo, iglo)
         it = it_idx(g_lo, iglo)
         if (damp_zonal_only .and. .not. aky(ik)==epsilon(0.0)) cycle
         hypervisc_filter(:,it,ik) = exp(- ( D_hypervisc * code_dt &
              * ( shear_rate_nz(:) *  (aky(ik) ** 2 + akx(it) ** 2 )**nexp / akperp4_max)))
         
         g0(:,1,iglo) = g0(:,1,iglo) * hypervisc_filter(:,it,ik)
         g0(:,2,iglo) = g0(:,2,iglo) * hypervisc_filter(:,it,ik)
      end do
      
    end subroutine iso_shear

  end subroutine hyper_diff

  subroutine finish_hyper

    implicit none

    hyper_on = .false.
    initialized = .false.
    if (allocated(D_res)) deallocate (D_res)
    if (allocated(hypervisc_filter)) deallocate (hypervisc_filter)

  end subroutine finish_hyper

end module hyper
