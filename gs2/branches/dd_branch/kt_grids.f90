module kt_grids_single
! <doc> Set up values of kx and ky for linear runs that use a single k_perp mode.
! </doc>

  implicit none

  public :: init_kt_grids_single, single_get_sizes, single_get_grids
  public :: check_kt_grids_single, wnml_kt_grids_single 

  private

  real :: aky, theta0, akx

contains

  subroutine init_kt_grids_single
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    integer :: in_file
    logical :: exist
    namelist /kt_grids_single_parameters/ aky, theta0, akx

    aky = 0.4
    theta0 = 0.0
    akx = 0.0
    in_file = input_unit_exist ("kt_grids_single_parameters", exist)
!    if (exist) read (unit=input_unit("kt_grids_single_parameters"), &
!         nml=kt_grids_single_parameters)
    if (exist) read (unit=in_file,nml=kt_grids_single_parameters)
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
    naky = 1
    ntheta0 = 1
    nx = 0
    ny = 0
  end subroutine single_get_sizes

  subroutine single_get_grids (aky, theta0, akx)
    implicit none
    real, dimension (:), intent (out) :: aky
    real, dimension (:,:), intent (out) :: theta0
    real, dimension (:), intent (out) :: akx
    call get_grids (aky, theta0, akx)
  end subroutine single_get_grids

  subroutine get_grids (aky_out, theta0_out, akx_out)
    implicit none
    real, dimension (:), intent (out) :: aky_out
    real, dimension (:,:), intent (out) :: theta0_out
    real, dimension (:), intent (out) :: akx_out
    aky_out(1) = aky
    theta0_out(1,1) = theta0
    akx_out(1) = akx
  end subroutine get_grids

  subroutine check_kt_grids_single(report_unit)
    implicit none
    integer:: report_unit

       write (report_unit, *) 
       write (report_unit, fmt="('A single k_perp will be evolved, with: ')")
       write (report_unit, *) 
       write (report_unit, fmt="('ky rho = ',f10.4)") aky
       write (report_unit, fmt="('theta0 = ',f10.4)") theta0
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

  integer :: naky, ntheta0
  real :: aky_min, aky_max, theta0_min, theta0_max
!CMR: add akx_min,max for periodic finite kx ballooning space runs with shat=0
  real :: akx_min, akx_max
    namelist /kt_grids_range_parameters/ naky, ntheta0, &
         aky_min, aky_max, theta0_min, theta0_max, akx_min, akx_max

contains

  subroutine init_kt_grids_range
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    integer :: in_file
    logical :: exist

    naky = 1
    ntheta0 = 1
    aky_min = 0.0
    aky_max = 0.0
    theta0_min = 0.0
    theta0_max = 0.0
    akx_min = 0.0
    akx_max = 0.0
    in_file = input_unit_exist ("kt_grids_range_parameters", exist)
    if (exist) read (unit=input_unit("kt_grids_range_parameters"), &
         nml=kt_grids_range_parameters)
  end subroutine init_kt_grids_range

  subroutine wnml_kt_grids_range(unit)
   implicit none
   integer:: unit
     write (unit, *)
     write (unit, fmt="(' &',a)") "kt_grids_range_parameters"
     write (unit, fmt="(' naky = ',i3)") naky
     write (unit, fmt="(' aky_min = ',e16.10)") aky_min
     write (unit, fmt="(' aky_max = ',e16.10)") aky_max
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
    naky_x = naky
    ntheta0_x = ntheta0
    nx = 0
    ny = 0
  end subroutine range_get_sizes

  subroutine range_get_grids (aky, theta0, akx)
    use theta_grid, only: shat
    implicit none
    real, dimension (:), intent (out) :: aky
    real, dimension (:,:), intent (out) :: theta0
    real, dimension (:), intent (out) :: akx
    real :: dkx
    real :: daky, dtheta0
    integer :: i, j

    if ( size(aky) .ne. naky) then
       write(6,*) 'range_get_grids: size(aky) /= naky'
       stop
    endif
    if ( size(akx) .ne. ntheta0) then
       write(6,*) 'range_get_grids: size(akx) /= ntheta0'
       stop
    endif

    daky = 0.0
    if (naky > 1) daky = (aky_max - aky_min)/real(naky - 1)
    aky = (/ (aky_min + daky*real(i), i=0,naky-1) /)

!CMR:
! set default theta0 to 0
    theta0=0.0d0
!CMRend

    if (shat .ne. 0.0d0) then  ! ie assumes boundary_option .eq. 'linked'
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
       akx = (/ (akx_min + dkx*real(i), i=0,ntheta0-1) /)
    endif
  end subroutine range_get_grids

  subroutine check_kt_grids_range(report_unit)
    use constants, only: twopi
    use theta_grid, only: shat
    implicit none
    integer:: report_unit
    real :: daky, dtheta0, dakx
    integer:: i,j

       write (report_unit, *) 
       write (report_unit, fmt="('A range of k_perps will be evolved.')")
       write (report_unit, *) 
       write (report_unit, fmt="('There are ',i3,' values of ky rho and ',i3,' values of theta0/kx rho:')") naky, ntheta0
       write (report_unit, *) 
          
       daky = 0.0
       if (naky > 1) daky = (aky_max - aky_min)/real(naky - 1)
       dtheta0 = 0.0
       if (ntheta0 > 1) dtheta0 = (theta0_max - theta0_min)/real(ntheta0 - 1)
       dakx = 0.0
       if (ntheta0 > 1) dakx = (akx_max - akx_min)/real(ntheta0 - 1)

       do j = 0, naky-1
          do i = 0, ntheta0-1
             write (report_unit, fmt="('ky rho = ',e10.4,' theta0 = ',e10.4,' kx rho = ',e10.4)") &
                  aky_min + daky*real(j), theta0_min + dtheta0*real(i), akx_min + dakx*real(i)
          end do
       end do
!CMR, add some !!!error checking!!! for ballooning space runs for shat/=0 
! using flow shear: check that the constraints on theta0 grid are satisfied!
       if (shat .ne. 0) then
         if (abs(mod(twopi-theta0_max+theta0_min,twopi)-dtheta0) .gt. 1.0e-3*dtheta0) then
             write (report_unit, *) 
             write (report_unit, fmt="('IF using perp ExB flow shear in BALLOONING SPACE there is a ERROR that will corrupt results.')")
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
  public :: check_kt_grids_specified
  public :: wnml_kt_grids_specified

  private

  integer :: naky, ntheta0, nx, ny
  namelist /kt_grids_specified_parameters/ naky, ntheta0, nx, ny

contains

  subroutine init_kt_grids_specified
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    integer :: in_file
    logical :: exist

    naky = 1
    ntheta0 = 1
    nx = 0
    ny = 0
    in_file = input_unit_exist("kt_grids_specified_parameters", exist)
    if (exist) read (unit=input_unit("kt_grids_specified_parameters"), &
         nml=kt_grids_specified_parameters)
  end subroutine init_kt_grids_specified

  subroutine wnml_kt_grids_specified(unit,aky,theta0)
    implicit none
    integer:: unit, i
    real, dimension(:):: aky
    real, dimension(:,:):: theta0
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

    naky_x = naky
    ntheta0_x = ntheta0
    nx_x = nx
    ny_x = ny
  end subroutine specified_get_sizes

  subroutine specified_get_grids (aky, theta0, akx)
    implicit none
    real, dimension (:), intent (out) :: aky
    real, dimension (:,:), intent (out) :: theta0
    real, dimension (:), intent (out) :: akx
    integer :: i, naky, ntheta0
    real :: aky_dummy, theta0_dummy, akx_dummy

    naky = size(aky)
    ntheta0 = size(akx)
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
    integer :: unit
    real :: aky, theta0, akx
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

  subroutine check_kt_grids_specified(report_unit,aky,theta0)
    implicit none
    integer :: report_unit
    real, dimension (:), intent (in) :: aky
    real, dimension (:), intent (in) :: theta0
    integer :: i
       write (report_unit, *) 
       write (report_unit, fmt="('A set of ',i3,' k_perps will be evolved.')") max(naky,ntheta0)
       write (report_unit, *) 
       do i=1, max(naky,ntheta0)
          write (report_unit, fmt="('ky rho = ',e10.4,' theta0 = ',e10.4)") &
               aky(i), theta0(i)
       end do
  end subroutine check_kt_grids_specified

end module kt_grids_specified

module kt_grids_box
! <doc> Set the perpendicular box size and resolution for linear or nonlinear runs.
! </doc>
  implicit none

  public :: init_kt_grids_box, box_get_sizes, box_get_grids
  public :: check_kt_grids_box
  public :: wnml_kt_grids_box
  public :: x0, y0 !RN> Caution: these are not broadcasted!

  private

  integer :: naky_private, ntheta0_private, nx_private, ny_private
  integer :: nkpolar_private
  integer :: jtwist
  real :: ly, y0, x0, rtwist

contains

  subroutine init_kt_grids_box
    use theta_grid, only: init_theta_grid, shat
    use file_utils, only: input_unit, input_unit_exist
    use constants
    implicit none
    integer :: naky, ntheta0, nx, ny, nkpolar
    integer :: in_file
    logical :: exist
    namelist /kt_grids_box_parameters/ naky, ntheta0, ly, nx, ny, jtwist, &
	y0, rtwist, x0, nkpolar

    call init_theta_grid

    nkpolar = 0
    naky = 0
    ntheta0 = 0
    ly = 0.0
    y0 = 2.0
    nx = 0
    ny = 0
!    jtwist = 1
    ! default jtwist -- MAB
    jtwist = max(int(2.0*pi*shat + 0.5),1)
    rtwist = 0.0
    x0 = 0.
    in_file = input_unit_exist("kt_grids_box_parameters", exist)
    if (exist) read (unit=input_unit("kt_grids_box_parameters"), nml=kt_grids_box_parameters)
    if (y0 < 0) then
       y0 = -1./y0
    end if
    if (ly == 0.) ly = 2.0*pi*y0
    if (naky == 0) naky = (ny-1)/3 + 1
    if (ntheta0 == 0) ntheta0 = 2*((nx-1)/3) + 1
    if (rtwist == 0.) rtwist = real(jtwist)
!    if (nkpolar == 0) nkpolar = naky  ! should be generalized later.  
                                       ! For now, assuming square domain with nx=ny
    if (nkpolar == 0) nkpolar = int(real(naky-1.)*sqrt(2.))
    
    nkpolar_private = nkpolar
    naky_private = naky
    ntheta0_private = ntheta0
    nx_private = nx
    ny_private = ny
  end subroutine init_kt_grids_box

  subroutine wnml_kt_grids_box(unit)
   implicit none
   integer:: unit, nx, ny
     nx=nx_private
     ny=ny_private
     write (unit, *)
     write (unit, fmt="(' &',a)") "kt_grids_box_parameters"
     write (unit, fmt="(' nx = ',i4)") nx
     write (unit, fmt="(' ny = ',i4)") ny
     write (unit, fmt="(' Ly = ',e16.10)") ly
     if (jtwist /= 1) then
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

!  subroutine box_get_grids (aky, theta0, akx, akpolar)
  subroutine box_get_grids (aky, theta0, akx, ikx, iky)
    use theta_grid, only: shat
    use constants
    implicit none
    real, dimension (:), intent (out) :: aky
    real, dimension (:,:), intent (out) :: theta0
    real, dimension (:), intent (out) :: akx
    integer, dimension (:), intent (out) :: ikx, iky
!    real, dimension (:), intent (out) :: akpolar
    integer :: i, naky, ntheta0
    real :: dkx

    naky = size(aky)
    ntheta0 = size(akx)

    if(abs(shat) <=  1.e-5) then
       if (x0 == 0.) then
          if (rtwist > 0) then
             dkx = 2.0*pi/ly/rtwist
          else
             dkx = 2.0*pi/ly*rtwist
          endif
       else
          if (x0 > 0.) then
             dkx = 1./x0
          else
             dkx = -x0
          end if
       end if
    else
       if (jtwist /= 0) then
          dkx = 2.0*pi/real(jtwist)* 2.0*pi/ly*shat
       else
          dkx = 2.0*pi/ly
       end if
    endif

    do i = 1, naky
       iky(i) = i-1
       aky(i) = real(i-1)*2.*pi/ly
    end do

    do i = 1, (ntheta0+1)/2
       ikx(i) = i-1
       akx(i) = real(i-1)*dkx
    end do

    do i = (ntheta0+1)/2+1, ntheta0
       ikx(i) = i-ntheta0-1
       akx(i) = real(i-ntheta0-1)*dkx
    end do

!!! Is there a negative sign missing here?  For up-down symmetric
!!! equilibria, it would not matter, but otherwise it could??
!!! This is mainly used to define wdrift and kperp2
    do i = 1, ntheta0
       theta0(i,2:) = akx(i)/(aky(2:)*shat)
       theta0(i,1) = 0.0
    end do
    
! currently assuming square domain with nx=ny for MHD energy diagnostics:
!    akpolar = aky
    
  end subroutine box_get_grids

  subroutine check_kt_grids_box(report_unit)
    use theta_grid, only: shat_real => shat
    use constants, only: pi
    implicit none
    integer:: report_unit
    real :: lx, naky, ntheta0, nx, ny, shat

    naky=naky_private
    ntheta0=ntheta0_private
    nx=nx_private
    ny=ny_private
    shat=shat_real

       if (y0 /= 2.) then
          if (abs(2.*pi*y0 - ly) > epsilon(0.)) then
             write (report_unit, *) 
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('You cannot specify both ly and y0.')")
             write (report_unit, fmt="('THIS IS AN ERROR.')") 
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 
          end if
       end if

       write (report_unit, *) 
       write (report_unit, fmt="('A rectangular simulation domain has been selected.')")
       write (report_unit, *) 
       write (report_unit, fmt="('The domain is ',f10.4,' rho in the y direction.')") ly

       if (abs(shat) <= 1.e-5) then
          if (x0 == 0.) then
             if (rtwist > 0) then
                write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)") abs(rtwist)*ly
             else
                write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)") ly/abs(rtwist)
             end if
          else
             if (x0 > 0) then
                write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)") 2.*pi*x0
             else
                write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)") -2.*pi/x0
             end if
          end if
       else
          lx = ly * rtwist / (2.*pi*shat)
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

module kt_grids_xbox
! <doc> Deprecated.
! </doc>
  implicit none

  public :: init_kt_grids_xbox, xbox_get_sizes, xbox_get_grids
  public :: wnml_kt_grids_xbox
  private

  integer :: ntheta0_private, nx_private
  real :: lx, aky_private

contains

  subroutine init_kt_grids_xbox
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    integer :: ntheta0, nx
    integer :: in_file
    logical :: exist
    real :: aky
    namelist /kt_grids_xbox_parameters/ ntheta0, lx, aky, nx

    ntheta0 = 1
    lx = 1.0
    aky = 0.2
    nx = 0
    in_file = input_unit_exist ("kt_grids_xbox_parameters", exist)
    if (exist) read (unit=input_unit("kt_grids_xbox_parameters"), &
         nml=kt_grids_xbox_parameters)
    ntheta0_private = ntheta0
    aky_private = aky
    nx_private = nx
  end subroutine init_kt_grids_xbox

  subroutine wnml_kt_grids_xbox(unit)
    implicit none
    integer:: unit, nx, ntheta0
    nx=nx_private
    ntheta0=ntheta0_private
     write (unit, *)
     write (unit, fmt="(' &',a)") "kt_grids_xbox_parameters"
     write (unit, fmt="(' nx = ',i4)") nx
     write (unit, fmt="(' ntheta0 = ',i4)") ntheta0
     write (unit, fmt="(' Lx = ',e16.10)") lx
     write (unit, fmt="(' /')")
  end subroutine wnml_kt_grids_xbox

  subroutine xbox_get_sizes (naky, ntheta0, nx, ny)
    implicit none
    integer, intent (out) :: naky, ntheta0, nx, ny
    naky = 1
    ntheta0 = ntheta0_private
    nx = nx_private
    ny = 0
  end subroutine xbox_get_sizes

  subroutine xbox_get_grids (aky, theta0, akx)
    use theta_grid, only: shat
    use constants
    implicit none
    real, dimension (:), intent (out) :: aky
    real, dimension (:,:), intent (out) :: theta0
    real, dimension (:), intent (out) :: akx
    integer :: i, ntheta0

    aky(1) = aky_private

    ntheta0 = size(akx)
    akx(:(ntheta0+1)/2) = (/ (real(2*(i-1))*pi/lx, i=1,(ntheta0+1)/2) /)
    akx((ntheta0+1)/2+1:) &
         = (/ (real(2*(i-ntheta0-1))*pi/lx, i=(ntheta0+1)/2+1,ntheta0) /)
    theta0(:,1) = akx(:)/(aky(1)*shat)
  end subroutine xbox_get_grids

end module kt_grids_xbox

module kt_grids
!  <doc> Set up the perpendicular wavenumbers by calling the appropriate sub-modules.  Also, set the normalizing
! length in the perpendicular directions, depending on whether the thermal velocity that appears in the gyroradius
! expression contains a factor of square root of two or not.
! </doc>
  implicit none

  public :: init_kt_grids, box, finish_kt_grids, check_kt_grids, wnml_kt
  public :: aky, theta0, akx, akr
  public :: aky_out, akx_out, akr_out
  public :: naky, ntheta0, nx, ny, reality
  public :: nkpolar 
  public :: ikx, iky
  public :: gridopt_switch, grid_option, norm_option
  public :: gridopt_single, gridopt_range, gridopt_specified, gridopt_box, gridopt_xbox
  private

  integer :: naky, ntheta0, nx, ny, nkpolar
  integer, dimension(:), allocatable :: ikx, iky
  character(20) :: grid_option, norm_option
  namelist /kt_grids_knobs/ grid_option, norm_option
  real, dimension (:), allocatable :: aky, aky_out
  real, dimension (:,:), allocatable :: theta0
  real, dimension (:), allocatable :: akx, akx_out
  real, dimension (:,:), allocatable :: akr, akr_out
!  real, dimension (:), allocatable :: akpolar, akpolar_out

  ! internal variables
  integer :: gridopt_switch, normopt_switch
  integer, parameter :: gridopt_single = 1, gridopt_range = 2, &
       gridopt_specified = 3, gridopt_box = 4, gridopt_xbox = 5
  integer, parameter :: normopt_mtk = 1, normopt_bd = 2
  logical :: reality = .false.
  logical :: box = .false.
  logical :: initialized = .false.
  logical :: nml_exist

contains

  subroutine init_kt_grids (tnorm)
    use theta_grid, only: init_theta_grid, shat, gds22
    use mp, only: proc0, broadcast
    implicit none

    real, optional, intent (out) :: tnorm
    integer :: ik, it
    real :: tfac = 1.0
!    logical, save :: initialized = .false.

    if (present(tnorm)) tnorm = tfac

    if (initialized) return
    initialized = .true.

    call init_theta_grid

    if (proc0) then
       nkpolar = 0   ! will be set to non-zero value only in box case; only used for an MHD diagnostic
       call read_parameters
       call get_sizes

       select case (normopt_switch)
       case (normopt_mtk)
          tfac = 1.
       case (normopt_bd)
          tfac = sqrt(2.)
       end select
    end if

    call broadcast (tfac)
    if (present(tnorm)) tnorm = tfac

    call broadcast (reality)
    call broadcast (box)
    call broadcast (naky)
    call broadcast (nkpolar)
    call broadcast (ntheta0)
    call broadcast (ny)
    call broadcast (nx)
    call allocate_arrays
    if (proc0) call get_grids
    call broadcast (ikx)     ! MR
    call broadcast (aky)
    call broadcast (akx)
    call broadcast (ikx)
    do ik = 1, naky
       call broadcast (theta0(:,ik))
    end do
!    if (nkpolar > 0) call broadcast (akpolar)

    if (abs(shat) > epsilon(0.)) then
       do it = 1, ntheta0
          akr(:,it) = akx(it)*sqrt(abs(gds22))/abs(shat)
       end do
    else
       akr = 1.
    end if

    select case (normopt_switch)
    case (normopt_mtk)
       akr_out = akr
       akx_out = akx
       aky_out = aky
!       if (nkpolar > 0) akpolar_out = akpolar
    case (normopt_bd)
       akr_out = akr / sqrt(2.)
       akx_out = akx / sqrt(2.)
       aky_out = aky / sqrt(2.)
!       if (nkpolar > 0) akpolar_out = akpolar / sqrt(2.)
    end select

  end subroutine init_kt_grids

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    implicit none
    type (text_option), dimension (7), parameter :: gridopts = &
         (/ text_option('default', gridopt_single), &
            text_option('single', gridopt_single), &
            text_option('range', gridopt_range), &
            text_option('specified', gridopt_specified), &
            text_option('box', gridopt_box), &
            text_option('nonlinear', gridopt_box), &
            text_option('xbox', gridopt_xbox) /)
    ! 'default' 'specified': specify grid in namelists
    type (text_option), dimension(6), parameter :: normopts = &
         (/ text_option('default', normopt_mtk), &
            text_option('with_root_2', normopt_mtk), &
            text_option('mtk', normopt_mtk), &
            text_option('no_root_2', normopt_bd), &
            text_option('bd', normopt_bd), &
            text_option('t_over_m', normopt_bd) /)
    integer :: ierr, in_file

    norm_option = 'default'
    grid_option = 'default'
    in_file = input_unit_exist ("kt_grids_knobs", nml_exist)
!    if (exist) read (unit=input_unit("kt_grids_knobs"), nml=kt_grids_knobs)
    if (nml_exist) read (unit=in_file, nml=kt_grids_knobs)

    ierr = error_unit()
    call get_option_value &
         (grid_option, gridopts, gridopt_switch, &
         ierr, "grid_option in kt_grids_knobs")

    ierr = error_unit()
    call get_option_value &
         (norm_option, normopts, normopt_switch, &
         ierr, "norm_option in kt_grids_knobs")

  end subroutine read_parameters

  subroutine wnml_kt(unit)
    use kt_grids_single, only: wnml_kt_grids_single
    use kt_grids_range, only: wnml_kt_grids_range
    use kt_grids_specified, only: wnml_kt_grids_specified
    use kt_grids_box, only: wnml_kt_grids_box
    use kt_grids_xbox, only: wnml_kt_grids_xbox
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
    case (gridopt_xbox)
       write (unit, fmt="(' grid_option = ',a)") '"xbox"'
    end select

    select case (normopt_switch)             
    case (normopt_mtk)
       write (unit, fmt="(' norm_option = ',a)") '"with_root_2"'
    case (normopt_bd)
       write (unit, fmt="(' norm_option = ',a)") '"no_root_2"'
    end select
    write (unit, fmt="(' /')")

    select case (gridopt_switch)          
    case (gridopt_single)
       call wnml_kt_grids_single(unit)
    case (gridopt_range)
       call wnml_kt_grids_range(unit)
    case (gridopt_specified)
       call wnml_kt_grids_specified(unit,aky,theta0)
    case (gridopt_box)
       call wnml_kt_grids_box(unit)
    case (gridopt_xbox)
       call wnml_kt_grids_xbox(unit)
    end select
  end subroutine wnml_kt


  subroutine allocate_arrays
    use theta_grid, only: ntgrid
    implicit none
    allocate (aky(naky))
    allocate (aky_out(naky))
    allocate (theta0(ntheta0,naky))
    allocate (akx(ntheta0))
    allocate (akx_out(ntheta0))
    allocate (akr(-ntgrid:ntgrid,ntheta0))
    allocate (akr_out(-ntgrid:ntgrid,ntheta0))

    allocate (ikx(ntheta0))
    allocate (iky(naky))

!    if (nkpolar > 0) &
!         allocate (akpolar(nkpolar), akpolar_out(nkpolar))
  end subroutine allocate_arrays

  subroutine get_sizes
    use kt_grids_single, only: init_kt_grids_single, single_get_sizes
    use kt_grids_range, only: init_kt_grids_range, range_get_sizes
    use kt_grids_specified, only: init_kt_grids_specified, specified_get_sizes
    use kt_grids_box, only: init_kt_grids_box, box_get_sizes
    use kt_grids_xbox, only: init_kt_grids_xbox, xbox_get_sizes
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
    case (gridopt_xbox)
       call init_kt_grids_xbox
       call xbox_get_sizes (naky, ntheta0, nx, ny)
    end select
  end subroutine get_sizes

  subroutine get_grids
    use kt_grids_single, only: single_get_grids
    use kt_grids_range, only: range_get_grids
    use kt_grids_specified, only: specified_get_grids
    use kt_grids_box, only: box_get_grids
    use kt_grids_xbox, only: xbox_get_grids
    implicit none
    select case (gridopt_switch)
    case (gridopt_single)
       call single_get_grids (aky, theta0, akx)
    case (gridopt_range)
       call range_get_grids (aky, theta0, akx)
    case (gridopt_specified)
       call specified_get_grids (aky, theta0, akx)
    case (gridopt_box)
!       call box_get_grids (aky, theta0, akx, akpolar)
       call box_get_grids (aky, theta0, akx, ikx, iky)
    case (gridopt_xbox)
       call xbox_get_grids (aky, theta0, akx)
    end select

    select case (normopt_switch)
    case (normopt_mtk)
       ! nothing -- this is how the original code is designed
    case (normopt_bd)
       aky = aky * sqrt(2.)
       akx = akx * sqrt(2.)
!       if (nkpolar > 0) akpolar = akpolar * sqrt(2.)
    end select

  end subroutine get_grids

  subroutine finish_kt_grids

    implicit none

    if (allocated(aky)) deallocate (aky, aky_out, theta0, akx, akx_out, akr, akr_out, ikx, iky)

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
    integer :: i

    select case (gridopt_switch) 
    case (gridopt_single) 
       call check_kt_grids_single(report_unit)
    case (gridopt_range)
       call check_kt_grids_range(report_unit)
    case (gridopt_specified) 
       call check_kt_grids_specified(report_unit,aky,theta0(:,1))
    case (gridopt_box)
       call check_kt_grids_box(report_unit)
    case (gridopt_xbox)
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('You selected grid_option=xbox in kt_grids_knobs')")
          write (report_unit, fmt="('The xbox option is not working.')")
          write (report_unit, fmt="('THIS IS AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
    end select
  end subroutine check_kt_grids

end module kt_grids

