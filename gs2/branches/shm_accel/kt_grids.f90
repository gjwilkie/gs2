module kt_grids_single
! <doc> Set up values of kx and ky for linear runs that use a single k_perp mode.
! </doc>

  implicit none

  public :: init_kt_grids_single, single_get_sizes, single_get_grids
  public :: check_kt_grids_single, wnml_kt_grids_single 

  private

  real :: akx, aky, theta0, rhostar_single
  integer :: n0

contains

  subroutine init_kt_grids_single
!CMR, 14/10/2013: 
! New namelist variables n0, rhostar_single to set aky using toroidal mode number.
! Toroidal modenumber used if n0> 0 prescribed in input file. 
    use file_utils, only: input_unit, input_unit_exist
    use theta_grid, only: drhodpsi
    implicit none
    integer :: in_file
    logical :: exist
    namelist /kt_grids_single_parameters/ n0, aky, theta0, akx, rhostar_single

    aky = 0.4   ;  theta0 = 0.0   ;   akx = 0.0
    n0 = 0      ; rhostar_single=1.0e-4

    in_file = input_unit_exist ("kt_grids_single_parameters", exist)
    if (exist) read (in_file, nml=kt_grids_single_parameters)

    if (n0 .gt. 0) then
!CMR if n0>0 then override aky inputs and use n0 to determine aky
       aky=n0*drhodpsi*rhostar_single
    endif

  end subroutine init_kt_grids_single

  subroutine wnml_kt_grids_single(unit)
   implicit none
   integer:: unit
     write (unit, *)
     write (unit, fmt="(' &',a)") "kt_grids_single_parameters"
     write (unit, fmt="(' aky = ',e16.10)") aky
     write (unit, fmt="(' theta0 = ',e16.10)") theta0
     write (unit, fmt="(' /')")
  end subroutine wnml_kt_grids_single

  subroutine single_get_sizes (naky, ntheta0, nx, ny)
    implicit none
    integer, intent (out) :: naky, ntheta0, nx, ny

    naky = 1  ;  ntheta0 = 1
    nx = 0    ;  ny = 0

  end subroutine single_get_sizes

  subroutine single_get_grids (aky_out, theta0_out, akx_out)
    implicit none
    real, dimension (:), intent (out) :: aky_out
    real, dimension (:,:), intent (out) :: theta0_out
    real, dimension (:), intent (out) :: akx_out

    aky_out = aky
    theta0_out = theta0
    akx_out = akx

  end subroutine single_get_grids

  subroutine check_kt_grids_single(report_unit)
    implicit none
    integer:: report_unit

       write (report_unit, *) 
       write (report_unit, fmt="('A single k_perp will be evolved, with: ')")
       if (n0 .gt.0) write (report_unit, fmt="('ky set using toroidal mode number, n0=',i8/T24,'rhostar_single=',1pe12.4)") n0, rhostar_single
       write (report_unit, *) 
       write (report_unit, fmt="('ky rho = ',f10.4)") aky
       write (report_unit, fmt="('theta_0 = ',f10.4)") theta0
       if (akx /= 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('The value of akx in the kt_grids_single_parameters namelist is ignored.')") 
          write (report_unit, fmt="('You have set akx to a non-zero value.')") 
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if
  end subroutine check_kt_grids_single

end module kt_grids_single

module kt_grids_range
! <doc> Set up ranges of kx and ky for linear runs.
! </doc>
  implicit none

  public :: init_kt_grids_range, range_get_sizes, range_get_grids
  public :: check_kt_grids_range
  public :: wnml_kt_grids_range

  private

  integer :: naky, ntheta0, nn0, n0_min, n0_max
  real :: aky_min, aky_max, theta0_min, theta0_max
  real :: akx_min, akx_max, rhostar_range
  namelist /kt_grids_range_parameters/ naky, ntheta0, nn0, n0_min, n0_max, &
       aky_min, aky_max, theta0_min, theta0_max, akx_min, akx_max, rhostar_range

contains

  subroutine init_kt_grids_range
!CMR, 14/10/2013: 
! New namelist variables nn0, n0_min, n0_max, rhostar_range to set ky grid 
!                                             using toroidal mode numbers.
! Toroidal modenumbers are used if n0_min> 0 prescribed in input file. 
    use file_utils, only: input_unit, input_unit_exist
    use theta_grid, only: drhodpsi
    implicit none
    integer :: in_file
    logical :: exist

    naky = 1          ;  ntheta0 = 1  
    aky_min = 0.0     ;  aky_max = 0.0
    theta0_min = 0.0  ;  theta0_max = 0.0
    akx_min = 0.0     ;  akx_max = 0.0
    nn0 = 1           ;  ntheta0 = 1  
    n0_min = 0        ;  n0_max = 0
    rhostar_range=1.0e-4

    in_file = input_unit_exist ("kt_grids_range_parameters", exist)
    if (exist) read (in_file, nml=kt_grids_range_parameters)
    if (n0_min .gt. 0) then
!CMR if n0_min>0 then override aky inputs and use nn0, n0_min, n0_max to determine aky range
       aky_min=n0_min*drhodpsi*rhostar_range
       if (mod(n0_max-n0_min,nn0).ne.0 .or. nn0 .eq. 1 .or. n0_min.ge.n0_max) then
          naky=1
          aky_max=aky_min
       else
          naky=nn0
          aky_max=n0_max*drhodpsi*rhostar_range
       endif 
    endif

  end subroutine init_kt_grids_range

  subroutine wnml_kt_grids_range(unit)
   implicit none
   integer:: unit
     write (unit, *)
     write (unit, fmt="(' &',a)") "kt_grids_range_parameters"
     write (unit, fmt="(' naky = ',i3)") naky
     write (unit, fmt="(' aky_min = ',e16.10)") aky_min
     write (unit, fmt="(' aky_max = ',e16.10)") aky_max
     write (unit, fmt="(' nn0 = ',i3)") nn0
     write (unit, fmt="(' n0_min = ',i10)") n0_min
     write (unit, fmt="(' n0_max = ',i10)") n0_max
     write (unit, fmt="(' rhostar_range = ',e16.10)") rhostar_range
     write (unit, fmt="(' ntheta0 = ',i3)") ntheta0
     write (unit, fmt="(' theta0_min = ',e16.10)") theta0_min
     write (unit, fmt="(' theta0_max = ',e16.10)") theta0_max
     write (unit, fmt="(' akx_min = ',e16.10)") akx_min
     write (unit, fmt="(' akx_max = ',e16.10)") akx_max
     write (unit, fmt="(' /')")
  end subroutine wnml_kt_grids_range

  subroutine range_get_sizes (naky_x, ntheta0_x, nx, ny)
    implicit none
    integer, intent (out) :: naky_x, ntheta0_x, nx, ny
    naky_x = naky  ;  ntheta0_x = ntheta0
    nx = 0         ;  ny = 0

  end subroutine range_get_sizes

! BD: Could add some logic here to set theta0 if akx is given?  When do we need what?

  subroutine range_get_grids (aky, theta0, akx)
    use theta_grid, only: shat
    implicit none
    real, dimension (:), intent (out) :: akx, aky
    real, dimension (:,:), intent (out) :: theta0

    real :: dkx, dky, dtheta0
    integer :: i, j

    if ( size(aky) /= naky) then
       write(6,*) 'range_get_grids: size(aky) /= naky'       ;  stop
    endif

    if ( size(akx) /= ntheta0) then
       write(6,*) 'range_get_grids: size(akx) /= ntheta0'    ;  stop
    endif

    dky = 0.0
    if (naky > 1) dky = (aky_max - aky_min)/real(naky - 1)
    aky = (/ (aky_min + dky*real(i), i = 0,naky-1) /)
 
! set default theta0 to 0
    theta0=0.0

!
! BD: Assumption here differs from convention that abs(shat) <= 1.e-5 triggers periodic bc
!
    if (shat /= 0.0d0) then  ! ie assumes boundary_option .eq. 'linked'
       dtheta0 = 0.0
       if (ntheta0 > 1) dtheta0 = (theta0_max - theta0_min)/real(ntheta0 - 1)

       do j = 1, naky
          theta0(:,j) &
               = (/ (theta0_min + dtheta0*real(i), i=0,ntheta0-1) /)
       end do
       akx = theta0(:,1) * shat * aky(1)

    else

!CMR, 22/9/2010:  ie here assume boundary_option .eq. 'periodic'
!new code for periodic finite kx ballooning space runs with shat=0

       dkx = 0.0
       if (ntheta0 > 1) dkx = (akx_max - akx_min)/real(ntheta0 - 1)
       akx = (/ (akx_min + dkx*real(i), i = 0,ntheta0-1) /)

    endif

  end subroutine range_get_grids

  subroutine check_kt_grids_range(report_unit)
    use constants, only: twopi
    use theta_grid, only: shat
    implicit none
    real :: dky, dtheta0, dkx
    integer :: report_unit, i, j

       write (report_unit, *) 
       write (report_unit, fmt="('A range of k_perps will be evolved.')")
       if (n0_min .gt.0) write (report_unit, fmt="('ky set using toroidal mode numbers with n0_min=',i8/T34,'rhostar_range=',1pe12.4)") n0_min,rhostar_range
       write (report_unit, *) 
       write (report_unit, fmt="('There are ',i3,' values of ky rho and ',i3,' values of theta_0/kx rho:')") naky, ntheta0
       write (report_unit, *) 
          
       dky = 0.0        ;  if (naky > 1) dky = (aky_max - aky_min)/real(naky - 1)
       dkx = 0.0        ;  if (ntheta0 > 1) dkx = (akx_max - akx_min)/real(ntheta0 - 1)
       dtheta0 = 0.0    ;  if (ntheta0 > 1) dtheta0 = (theta0_max - theta0_min)/real(ntheta0 - 1)

       do j = 0, naky-1
          do i = 0, ntheta0-1
             write (report_unit, fmt="('ky rho = ',e10.4,' theta0 = ',e10.4,' kx rho = ',e10.4)") &
                  aky_min + dky*real(j), theta0_min + dtheta0*real(i), akx_min + dkx*real(i)
          end do
       end do

! CMR, add some !!!error checking!!! for ballooning space runs for shat /= 0 
! using flow shear: check that the constraints on theta0 grid are satisfied!

       if (shat /= 0) then
         if (abs(mod(twopi-theta0_max+theta0_min,twopi)-dtheta0) > 1.0e-3*dtheta0) then
             write (report_unit, *) 
             write (report_unit, fmt="('IF using perp ExB flow shear in BALLOONING SPACE there is an ERROR that will corrupt results.')")
             write (report_unit, fmt="('check_kt_grids_range: inappropriate theta0 grid')")
             write (report_unit, fmt="('In ballooning space with sheared flow, 2pi-theta0_max+theta0_min =',e10.4,' must be set equal to dtheta = ',e10.4)") twopi-theta0_max+theta0_min, dtheta0
         endif
       endif

   end subroutine check_kt_grids_range

end module kt_grids_range

module kt_grids_specified
! <doc> Set up sets of (kx, ky) values for linear runs.
! </doc>
  implicit none

  public :: init_kt_grids_specified, specified_get_sizes, specified_get_grids
  public :: check_kt_grids_specified, wnml_kt_grids_specified

  private

  integer :: naky, ntheta0, nx, ny
  namelist /kt_grids_specified_parameters/ naky, ntheta0, nx, ny

contains

  subroutine init_kt_grids_specified
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    integer :: in_file
    logical :: exist

    naky = 1    ;  ntheta0 = 1
    nx = 0      ;  ny = 0

    in_file = input_unit_exist("kt_grids_specified_parameters", exist)
    if (exist) read (in_file, nml=kt_grids_specified_parameters)

  end subroutine init_kt_grids_specified

  subroutine wnml_kt_grids_specified (unit, aky, theta0)
    implicit none
    real, dimension(:) :: aky
    real, dimension(:,:) :: theta0
    integer :: unit, i
    character (200) :: line

     write(unit, kt_grids_specified_parameters)
     do i=1,max(naky,ntheta0)
        write (unit, *)
        write (line, *) i
        write (unit, fmt="(' &',a)") &
             & trim("kt_grids_specified_element_"//trim(adjustl(line)))
        write (unit, fmt="(' aky = ',e13.6,' theta0 = ',e13.6,'  /')") aky(i), theta0(i,1)
        write (unit, fmt="(' /')")
     end do
  end subroutine wnml_kt_grids_specified

  subroutine specified_get_sizes (naky_x, ntheta0_x, nx_x, ny_x)
    implicit none
    integer, intent (out) :: naky_x, ntheta0_x, nx_x, ny_x

    naky_x = naky  ;   ntheta0_x = ntheta0
    nx_x = nx      ;   ny_x = ny

  end subroutine specified_get_sizes

  subroutine specified_get_grids (aky, theta0, akx)
    implicit none
    real, dimension (:), intent (out) :: akx, aky
    real, dimension (:,:), intent (out) :: theta0
    real :: aky_dummy, theta0_dummy, akx_dummy
    integer :: i, naky, ntheta0

    naky = size(aky)  ;  ntheta0 = size(akx)

    do i = 1, max(naky,ntheta0)
       call read_element (i, aky_dummy, theta0_dummy, akx_dummy)
       if (i <= naky) aky(i) = aky_dummy
       if (i <= ntheta0) theta0(i,:) = theta0_dummy
       if (i <= ntheta0) akx(i) = akx_dummy
    end do

  end subroutine specified_get_grids

  subroutine read_element (i, aky_dummy, theta0_dummy, akx_dummy)
    use file_utils, only: get_indexed_namelist_unit
    implicit none
    integer, intent (in) :: i
    real, intent (out) :: aky_dummy, theta0_dummy, akx_dummy
    real :: akx, aky, theta0
    integer :: unit

    namelist /kt_grids_specified_element/ aky, theta0, akx

    aky = 0.4
    theta0 = 0.0
    akx = 0.0
    call get_indexed_namelist_unit (unit, "kt_grids_specified_element", i)
    read (unit=unit, nml=kt_grids_specified_element)
    close (unit)
    aky_dummy = aky
    theta0_dummy = theta0
    akx_dummy = akx
  end subroutine read_element

  subroutine check_kt_grids_specified (report_unit, aky, theta0)
    implicit none
    integer :: report_unit
    real, dimension (:), intent (in) :: aky
    real, dimension (:), intent (in) :: theta0
    integer :: i

    write (report_unit, *) 
    write (report_unit, fmt="('A set of ',i3,' k_perps will be evolved.')") max(naky,ntheta0)
    write (report_unit, *) 
    do i=1, max(naky,ntheta0)
       write (report_unit, fmt="('ky rho = ',e10.4,' theta0 = ',e10.4)") aky(i), theta0(i)
    end do
  end subroutine check_kt_grids_specified

end module kt_grids_specified

module kt_grids_box
! <doc> Set the perpendicular box size and resolution for linear or nonlinear runs.
! </doc>
  implicit none

  public :: init_kt_grids_box, box_get_sizes, box_get_grids
  public :: check_kt_grids_box, wnml_kt_grids_box
  public :: x0, y0, jtwist !RN> Caution: these are not broadcasted!

  private

  real :: ly, y0, x0, rtwist, rhostar_box
  integer :: naky_private, ntheta0_private, nx_private, ny_private
  integer :: nkpolar_private, n0
  integer :: jtwist

contains

  subroutine init_kt_grids_box
!CMR, 14/10/2013: 
! New namelist variables: n0, rhostar_box. 
! If n0 and rhostar_box defined, set ky(1) using toroidal mode number.

    use theta_grid, only: init_theta_grid, shat, drhodpsi
    use file_utils, only: input_unit, input_unit_exist
    use constants
    implicit none
    integer :: naky, ntheta0, nx, ny, nkpolar
    integer :: in_file
    logical :: exist
    namelist /kt_grids_box_parameters/ naky, ntheta0, ly, nx, ny, n0, jtwist, &
         y0, rtwist, x0, nkpolar, rhostar_box

    call init_theta_grid

    nkpolar = 0   ;   naky = 0    ;  ntheta0 = 0
    ly = 0.0      ;   y0 = 2.0    ;  x0 = 0.
    nx = 0        ;   ny = 0      
    n0=0          ;   rhostar_box=0.0

    jtwist = max(int(2.0*pi*shat + 0.5),1)  ! default jtwist -- MAB
    rtwist = 0.0

    in_file = input_unit_exist("kt_grids_box_parameters", exist)
    if (exist) read (in_file, nml=kt_grids_box_parameters)

    if (rhostar_box .gt. 0.0 .and. n0 .gt. 0) y0=1.0/(n0*rhostar_box*drhodpsi)

    if (y0 < 0) y0 = -1./y0

    if (ly == 0.) ly = 2.0*pi*y0
    if (naky == 0) naky = (ny-1)/3 + 1
    if (ntheta0 == 0) ntheta0 = 2*((nx-1)/3) + 1
    if (rtwist == 0.) rtwist = real(jtwist)
    if (nkpolar == 0) nkpolar = int(real(naky-1.)*sqrt(2.))
    
    nkpolar_private = nkpolar
    naky_private = naky
    ntheta0_private = ntheta0
    nx_private = nx
    ny_private = ny
  end subroutine init_kt_grids_box

  subroutine wnml_kt_grids_box (unit)
   implicit none
   integer :: unit

     write (unit, *)
     write (unit, fmt="(' &',a)") "kt_grids_box_parameters"
     write (unit, fmt="(' nx = ',i4)") nx_private
     write (unit, fmt="(' ny = ',i4)") ny_private
     write (unit, fmt="(' Ly = ',e16.10)") ly
     if (rtwist /= 0.) then
        write (unit, fmt="(' rtwist = ',e16.10)") rtwist
     else
        write (unit, fmt="(' jtwist = ',i4)") jtwist
     end if
     write (unit, fmt="(' /')")
  end subroutine wnml_kt_grids_box

  subroutine box_get_sizes (naky, ntheta0, nx, ny, nkpolar)
    implicit none
    integer, intent (out) :: naky, ntheta0, nx, ny, nkpolar
    naky = naky_private
    ntheta0 = ntheta0_private
    nx = nx_private
    ny = ny_private
    nkpolar = nkpolar_private
  end subroutine box_get_sizes

  subroutine box_get_grids (aky, theta0, akx, ikx, iky)
    use theta_grid, only: shat
    use constants
    implicit none
    real, dimension (:), intent (out) :: akx, aky
    real, dimension (:,:), intent (out) :: theta0
    integer, dimension (:), intent (out) :: ikx, iky

    real :: dkx, dky, ratio
    integer :: i, naky, ntheta0

    naky = size(aky)    
    ntheta0 = size(akx)

    dky = 2.0*pi/ly

    if(abs(shat) <=  1.e-5) then   ! non-quantized b/c assumed to be periodic instead linked boundary conditions       

       if (x0 == 0.) then          
          
          if (rtwist > 0) then 
             ratio = rtwist
          else 
             ratio = 1. / abs(rtwist)
          end if
          
          dkx = dky / ratio
          
       else

          if (x0 > 0.) then
             dkx = 1./x0
          else
             dkx = -x0
          end if

       end if

    else
       if (jtwist /= 0) then
          dkx = dky * 2.0*pi*abs(shat)/real(jtwist)
       else
          dkx = dky
       end if
    endif

    do i = 1, naky
       iky(i) = i-1
       aky(i) = real(i-1)*dky
    end do

    do i = 1, (ntheta0+1)/2
       ikx(i) = i-1
       akx(i) = real(i-1)*dkx
    end do

    do i = (ntheta0+1)/2+1, ntheta0
       ikx(i) = i-ntheta0-1
       akx(i) = real(i-ntheta0-1)*dkx
    end do

    if (shat /= 0.) then
       do i = 1, ntheta0
          theta0(i,1) = 0.0
          theta0(i,2:) = akx(i)/(aky(2:)*shat)
       end do
    else
       do i = 1, ntheta0
          theta0(i,1) = 0.0
          theta0(i,2:) = - akx(i)/aky(2:)   ! meaningless, so be careful
       end do
    end if

  end subroutine box_get_grids

  subroutine check_kt_grids_box(report_unit)
    use theta_grid, only: shat_real => shat
    use constants, only: pi
    implicit none
    integer :: report_unit
    real :: lx, shat
    integer :: nx, ny, naky, ntheta0

    naky=naky_private
    ntheta0=ntheta0_private
    nx=nx_private
    ny=ny_private
    shat=shat_real

    if (y0 /= 2.) then
       if (abs(2.0*pi*y0 - ly) > 1.0e-7) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('You cannot specify both ly and y0.')")
          write (report_unit, fmt="(' ly=',e12.4,'  2.0*pi*y0=',2e12.4)") ly
          write (report_unit, fmt="('THIS IS AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if
    end if

    write (report_unit, *) 
    write (report_unit, fmt="('A rectangular simulation domain has been selected.')")
    if (rhostar_box .gt. 0.0 .and. n0 .gt. 0) write (report_unit, fmt="('The flux-tube size corresponds to toroidal mode number, n0=',i8/T44,' at rhostar_box=',1pe12.4)") n0,rhostar_box
    write (report_unit, *) 
    write (report_unit, fmt="('The domain is ',f10.4,' rho in the y direction.')") ly

    
    if (abs(shat) <= 1.e-5) then
       if (x0 == 0.) then
          if (rtwist > 0) then
             write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)")  ly*rtwist
          else
             write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)") -ly/rtwist
          end if
       else
          if (x0 > 0) then
             write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)")  2.*pi*x0
          else
             write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)") -2.*pi/x0
          end if
       end if
    else
       lx = ly * jtwist / (2.*pi*abs(shat))
       write (report_unit, fmt="('At theta=0, the domain is ',f10.4,' rho in the x direction.')") lx
       write (report_unit,*) ly, rtwist, jtwist, pi, shat
    end if
    
    write (report_unit, *) 
    write (report_unit, fmt="('The nonlinear terms will be evaluated on a grid with ',&
         & i4,' points in x and ',i4,' points in y.')") nx, ny
    write (report_unit, *) 
    naky = (ny-1)/3+1
    ntheta0 = 2*((nx-1)/3)+1
    write (report_unit, fmt="('After de-aliasing, there will be ',i4,'  ky >= 0 modes and ',i4,' kx modes.')") naky, ntheta0
    write (report_unit, fmt="('The modes with ky < 0 are determined by the reality condition.')")
  end subroutine check_kt_grids_box

end module kt_grids_box

module kt_grids
!  <doc> Set up the perpendicular wavenumbers by calling the appropriate sub-modules. 
! </doc>
  use kt_grids_box, only: jtwist
  implicit none

  public :: init_kt_grids, box, finish_kt_grids, check_kt_grids, wnml_kt
  public :: aky, theta0, akx
  public :: naky, ntheta0, nx, ny, reality
  public :: nkpolar 
  public :: ikx, iky, jtwist_out
  public :: gridopt_switch, grid_option
  public :: gridopt_single, gridopt_range, gridopt_specified, gridopt_box
  public :: kwork_filter
  private

  logical, dimension(:,:), allocatable :: kwork_filter
  real, dimension (:,:), allocatable :: theta0
  real, dimension (:), allocatable :: aky, akx
  integer, dimension(:), allocatable :: ikx, iky
  integer :: naky, ntheta0, nx, ny, nkpolar, jtwist_out
  character(20) :: grid_option

  namelist /kt_grids_knobs/ grid_option

  ! internal variables
  integer :: gridopt_switch
  integer, parameter :: gridopt_single = 1, gridopt_range = 2, &
       gridopt_specified = 3, gridopt_box = 4
  logical :: reality = .false.
  logical :: box = .false.
  logical :: initialized = .false.
  logical :: nml_exist

contains

  subroutine init_kt_grids
    use theta_grid, only: init_theta_grid
    use mp, only: proc0, broadcast
    implicit none

    integer :: ik

    if (initialized) return
    initialized = .true.

    call init_theta_grid

    if (proc0) then
       nkpolar = 0   ! will be set to non-zero value only in box case; only used for an MHD diagnostic
       call read_parameters
       call get_sizes
       jtwist_out = jtwist
    end if

    call broadcast (reality)
    call broadcast (box)
    call broadcast (naky)
    call broadcast (nkpolar)
    call broadcast (ntheta0)
    call broadcast (ny)
    call broadcast (nx)
    call broadcast (gridopt_switch)
    call allocate_arrays

    if (proc0) call get_grids
    call broadcast (ikx)     ! MR
    call broadcast (aky)
    call broadcast (akx)
    call broadcast (jtwist_out)
    do ik = 1, naky
       call broadcast (theta0(:,ik))
    end do
    allocate(kwork_filter(ntheta0,naky))
    kwork_filter=.false.
  end subroutine init_kt_grids

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    implicit none
    type (text_option), dimension (6), parameter :: gridopts = &
         (/ text_option('default', gridopt_single), &
            text_option('single', gridopt_single), &
            text_option('range', gridopt_range), &
            text_option('specified', gridopt_specified), &
            text_option('box', gridopt_box), &
            text_option('nonlinear', gridopt_box) /)
    integer :: ierr, in_file

    grid_option = 'default'
    in_file = input_unit_exist ("kt_grids_knobs", nml_exist)
    if (nml_exist) read (unit=in_file, nml=kt_grids_knobs)

    ierr = error_unit()
    call get_option_value (grid_option, gridopts, gridopt_switch, &
         ierr, "grid_option in kt_grids_knobs")

  end subroutine read_parameters

  subroutine wnml_kt(unit)
    use kt_grids_single, only: wnml_kt_grids_single
    use kt_grids_range, only: wnml_kt_grids_range
    use kt_grids_specified, only: wnml_kt_grids_specified
    use kt_grids_box, only: wnml_kt_grids_box
    implicit none
    integer :: unit

    if (.not. nml_exist) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "kt_grids_knobs"
    select case (gridopt_switch)          
    case (gridopt_single)
       write (unit, fmt="(' grid_option = ',a)") '"single"'
    case (gridopt_range)
       write (unit, fmt="(' grid_option = ',a)") '"range"'
    case (gridopt_specified)
       write (unit, fmt="(' grid_option = ',a)") '"specified"'
    case (gridopt_box)
       write (unit, fmt="(' grid_option = ',a)") '"box"'
    end select

    select case (gridopt_switch)          
    case (gridopt_single)
       call wnml_kt_grids_single (unit)
    case (gridopt_range)
       call wnml_kt_grids_range (unit)
    case (gridopt_specified)
       call wnml_kt_grids_specified (unit, aky, theta0)
    case (gridopt_box)
       call wnml_kt_grids_box (unit)
    end select
  end subroutine wnml_kt

  subroutine allocate_arrays
    implicit none
    allocate (akx(ntheta0))
    allocate (aky(naky))
    allocate (theta0(ntheta0,naky))
    allocate (ikx(ntheta0))
    allocate (iky(naky))
  end subroutine allocate_arrays

  subroutine get_sizes
    use kt_grids_single, only: init_kt_grids_single, single_get_sizes
    use kt_grids_range, only: init_kt_grids_range, range_get_sizes
    use kt_grids_specified, only: init_kt_grids_specified, specified_get_sizes
    use kt_grids_box, only: init_kt_grids_box, box_get_sizes
    implicit none
    select case (gridopt_switch)
    case (gridopt_single)
       call init_kt_grids_single
       call single_get_sizes (naky, ntheta0, nx, ny)
    case (gridopt_range)
       call init_kt_grids_range
       call range_get_sizes (naky, ntheta0, nx, ny)
    case (gridopt_specified)
       call init_kt_grids_specified
       call specified_get_sizes (naky, ntheta0, nx, ny)
    case (gridopt_box)
       call init_kt_grids_box
       call box_get_sizes (naky, ntheta0, nx, ny, nkpolar)
       reality = .true.
       box = .true.
    end select
  end subroutine get_sizes

  subroutine get_grids
    use kt_grids_single, only: single_get_grids
    use kt_grids_range, only: range_get_grids
    use kt_grids_specified, only: specified_get_grids
    use kt_grids_box, only: box_get_grids
    implicit none
    select case (gridopt_switch)
    case (gridopt_single)
       call single_get_grids (aky, theta0, akx)
    case (gridopt_range)
       call range_get_grids (aky, theta0, akx)
    case (gridopt_specified)
       call specified_get_grids (aky, theta0, akx)
    case (gridopt_box)
       call box_get_grids (aky, theta0, akx, ikx, iky)
    end select
  end subroutine get_grids

  subroutine finish_kt_grids

    implicit none

    if (allocated(aky)) deallocate (akx, aky, theta0, ikx, iky)
    if (allocated(kwork_filter)) deallocate(kwork_filter)
    reality = .false. ; box = .false.
    initialized = .false.

  end subroutine finish_kt_grids

  subroutine check_kt_grids(report_unit)
!CMR, 22/9/2010:
! add routine to GS2 to perform ingen checking and reporting
    use kt_grids_single, only: check_kt_grids_single
    use kt_grids_range,  only: check_kt_grids_range
    use kt_grids_box,    only: check_kt_grids_box
    use kt_grids_specified, only: check_kt_grids_specified
    implicit none
    integer :: report_unit

    select case (gridopt_switch) 
    case (gridopt_single) 
       call check_kt_grids_single (report_unit)
    case (gridopt_range)
       call check_kt_grids_range (report_unit)
    case (gridopt_specified) 
       call check_kt_grids_specified (report_unit, aky, theta0(:,1))
    case (gridopt_box)
       call check_kt_grids_box (report_unit)
    end select

  end subroutine check_kt_grids

end module kt_grids

