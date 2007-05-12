! nexp has been changed in a rush.  Only known to be correct for nexp=2
!
!

module hyper
  
  implicit none
  private
  public :: init_hyper, hyper_diff, D_res
  public :: D_v, D_eta, nexp

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
  logical :: const_amp, include_kpar, isotropic_shear
  logical :: hyper_on = .false.
  logical :: gridnorm

  real, dimension (:,:), allocatable, save :: D_res
  ! (it, ik)

  logical :: initialized = .false.

contains

  subroutine init_hyper

    use kt_grids, only: ntheta0, naky, akx, aky
    use run_parameters, only: tnorm
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
             D_res(it, ik) = D_hyperres*code_dt/tnorm &
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
    use text_options
    use mp, only: proc0, broadcast
    implicit none
    type (text_option), dimension(5), parameter :: hyperopts = &
         (/ text_option('default', hyper_option_none), &
            text_option('none', hyper_option_none), &
            text_option('visc_only', hyper_option_visc), &
            text_option('res_only', hyper_option_res), &
            text_option('both', hyper_option_both) /)
    character(9) :: hyper_option
    integer :: ierr, in_file
    logical :: exist
    
    namelist /hyper_knobs/ hyper_option, const_amp, include_kpar, &
         isotropic_shear, D_hyperres, D_hypervisc, D_hyper, omega_osc, gridnorm, nexp
    
    if (proc0) then
       const_amp = .false.
       include_kpar = .false.
       isotropic_shear = .true.
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
            ierr, "hyper_option in hyper_knobs")

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
    use kt_grids, only: ntheta0, naky
    implicit none
    logical :: alloc = .true.

    if (alloc) then
       allocate (D_res(ntheta0, naky)) 
    endif
    D_res = 0.
    alloc = .false.

  end subroutine allocate_arrays

  subroutine hyper_diff (g0, g1, phi, bpar)

    use mp, only: proc0
    use gs2_layouts, only: ik_idx, it_idx, is_idx
    use theta_grid, only: ntgrid
    use run_parameters, only: tnorm, fphi, fbpar
    use gs2_time, only: code_dt
    use gs2_layouts, only: g_lo, geint_lo, gint_lo
    use kt_grids, only: aky, akx, naky, ntheta0
    use le_grids, only: integrate, geint2g
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g0, g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    complex, dimension (:,:), allocatable :: g0eint, g1eint
    real, dimension (:,:), allocatable, save :: aintnorm
    logical :: first = .true.
    logical :: conserve_number = .false.
    logical :: conserve_inth = .true.

    real, dimension (-ntgrid:ntgrid) :: shear_rate_nz, &
         shear_rate_z, shear_rate_z_nz

    integer :: iglo, ik, it, ige
    integer :: ncall = 0 ! variables declared with value are automatically saved
 
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


    if (conserve_number) then ! get initial particle number
       allocate (g0eint(-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_alloc))
       allocate (g1eint(-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_alloc))

       if (first) then
          first = .false.
          allocate (aintnorm(-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_alloc))
          aintnorm = 0. ; g1 = 1.0
          call integrate (g1, g1eint)
          do ige = geint_lo%llim_proc, geint_lo%ulim_proc
             aintnorm(:,ige) = 1.0/real(g1eint(:,ige))
          end do
       end if

       if (conserve_inth) call g_adjust (g0, phi, bpar, fphi, fbpar)
       call integrate (g0, g0eint)
       if (conserve_inth) call g_adjust (g0, phi, bpar, -fphi, -fbpar)

    end if

    if(isotropic_shear) then
       call iso_shear
    else
       call aniso_shear
    end if

    if (conserve_number) then ! restore particle number
       if (conserve_inth) call g_adjust (g1, phi, bpar, fphi, fbpar)
       call integrate (g1, g1eint)
       if (conserve_inth) call g_adjust (g1, phi, bpar, -fphi, -fbpar)
       do ige = geint_lo%llim_proc, geint_lo%ulim_proc
          g0eint(:,ige) = (g0eint(:,ige) - g1eint(:,ige)) &
               *aintnorm(:,ige)
       end do
       call geint2g (g0eint, g0)
       g0 = g0 + g1
       deallocate (g0eint, g1eint)
    else
       g0 = g1
    endif

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
            
            g1(:,1,iglo) = g0(:,1,iglo) * exp(- ( D_hypervisc * code_dt / tnorm &
                 * ( shear_rate_z_nz(:) *  akx(it) ** 4 / akx4_max )))
            
            g1(:,2,iglo) = g0(:,2,iglo) * exp(- ( D_hypervisc * code_dt / tnorm & 
                 * ( shear_rate_z_nz(:) *  akx(it) ** 4 / akx4_max )))
            
         else
            
            g1(:,1,iglo) = g0(:,1,iglo) * exp(- ( D_hypervisc * code_dt / tnorm & 
                 * ( shear_rate_nz(:) *  (aky(ik) ** 2 + akx(it) ** 2 )**nexp / akperp4_max & 
                 + shear_rate_z(:) * akx(it) ** 4/ akx4_max * aky(ik) / aky_max )))
            
            g1(:,2,iglo) = g0(:,2,iglo) * exp(- ( D_hypervisc * code_dt / tnorm &
                 * ( shear_rate_nz(:) *  (aky(ik) ** 2 + akx(it) ** 2 )**nexp / akperp4_max & 
                 + shear_rate_z(:) * akx(it) ** 4 / akx4_max * aky(ik) / aky_max )))
            
         endif
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
         
         g1(:,1,iglo) = g0(:,1,iglo) * exp(- ( D_hypervisc * code_dt / tnorm &
              * ( shear_rate_nz(:) *  (aky(ik) ** 2 + akx(it) ** 2 )**nexp / akperp4_max)))
         
         g1(:,2,iglo) = g0(:,2,iglo) * exp(- ( D_hypervisc * code_dt / tnorm &
              * ( shear_rate_nz(:) *  (aky(ik) ** 2 + akx(it) ** 2 )**nexp / akperp4_max)))
      end do
      
    end subroutine iso_shear
  end subroutine hyper_diff

  subroutine g_adjust (g, phi, bpar, facphi, facbpar)
    use species, only: spec
    use theta_grid, only: ntgrid
    use le_grids, only: anon
    use dist_fn_arrays, only: vperp2, aj0, aj1
    use gs2_layouts, only: g_lo, ik_idx, it_idx, ie_idx, is_idx
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    real, intent (in) :: facphi, facbpar

    integer :: iglo, ig, ik, it, ie, is
    complex :: adj

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do ig = -ntgrid, ntgrid
          adj = anon(ie,is)*2.0*vperp2(ig,iglo)*aj1(ig,iglo) &
                  *bpar(ig,it,ik)*facbpar &
               + spec(is)%z*anon(ie,is)*phi(ig,it,ik)*aj0(ig,iglo) &
                  /spec(is)%temp*facphi
          g(ig,1,iglo) = g(ig,1,iglo) + adj
          g(ig,2,iglo) = g(ig,2,iglo) + adj
       end do
    end do
  end subroutine g_adjust


end module hyper
