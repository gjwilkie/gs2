!> Reads in the geometry using a CHEASE output file.
!! The CHEASE output file is read using the helper module
!! read_chease.
module ceq

  implicit none

  public :: ceq_init, ceqin, gradient, eqitem, bgradient, Hahm_Burrell, ceq_finish
  public :: invR,     initialize_invR
  public :: Rpos
  public :: Zpos
  public :: rcenter,  initialize_rcenter 
  public :: diameter, initialize_diameter
  public :: btori,    initialize_btori
  public :: dbtori,   initialize_dbtori
  public :: qfun,     initialize_q
  public :: pfun,     initialize_pressure
  public :: dpfun,    initialize_dpressure
  public :: betafun,  initialize_beta
  public :: psi,      initialize_psi

  private

  integer :: nr, nt, i_sym
  
  real, allocatable, dimension (:)     :: rho_d, eqpsi, psi_bar, fp, beta
  real, allocatable, dimension (:)     :: pressure, diam, rc, qsf, rho_b
  real, allocatable, dimension (:,:)   :: R_psi, Z_psi !, B_psi
  real, allocatable, dimension (:,:,:) :: drm, dzm, dbtm, dpm, dtm  !, dbm
  real, allocatable, dimension (:,:,:) :: dpcart, dbcart, dtcart, dbtcart
  real, allocatable, dimension (:,:,:) :: dpbish, dbbish, dtbish, dbtbish

  real :: psi_0, psi_a, B_T, beta_0, B_norm_chease
  real :: R_mag, Z_mag, aminor
  
  logical :: init_rcenter = .true.
  logical :: init_diameter = .true.
  logical :: init_btori = .true.
  logical :: init_dbtori = .true.
  logical :: init_q = .true.
  logical :: init_pressure = .true.
  logical :: init_dpressure = .true.
  logical :: init_beta = .true.
  logical :: init_psi = .true.
  logical :: init_invR = .true.
  logical :: transp = .false.

  logical, parameter :: debug = .true.

contains
  function chease_chi_index(nchi,itheta)
    integer, intent(in) :: nchi, itheta
    integer :: chease_chi_index
    ! when ichi = 0, itheta = (ntheta-1)/2 + 1
    ! when itheta = 0, ichi = ntheta / 2
    ! itheta  1 2 3 4 5 6 7 8 
    ! ichi    5 6 7 8 1 2 3 4         
    if (itheta > nchi/2) then 
       chease_chi_index = itheta - nchi/2
    else
       chease_chi_index = itheta + nchi/2
    end if
  end function chease_chi_index
    
  subroutine ceqin(eqfile, psi_0_out, psi_a_out, rmaj, B_T0, &
       avgrmid, initeq, in_nt, nthg) 
    use constants, only: pi
    use read_chease, only: read_infile, npsi_chease, nchi_chease, b0exp_chease
    use read_chease, only: psi_chease, f_chease, q_chease, p_chease
    use read_chease, only: r_chease, z_chease
    implicit none
    character(len=80), intent(in) :: eqfile
    real, intent(out) :: psi_0_out, psi_a_out, rmaj, B_T0, avgrmid
    integer, intent(in) :: initeq
    logical, intent(in) :: in_nt 
    integer, intent(out) :: nthg
    real :: d, R_geo
    integer :: i,j
    integer :: nchar
    character (len=80) :: filename

!
! what is the best way to handle the netcdf single/double problem?
!
    real :: f_N, psi_N


!    Read the data

    if(initeq == 0) return
    if (.not.in_nt) then
       nchar=index(eqfile,' ')-1
       filename=eqfile(1:nchar)
!       filename=trim(eqfile) ?
       write (*,*)  'Reading CHEASE input file: ', eqfile
       call read_infile(filename)
!     netcdf read scalar: nr
!
!     nz2 == number of points in radial array
!     nr == number of actual grid points in radial array

!       id = ncdid (ncid, 'z2', ifail)
!       call ncdinq (ncid, id, fortrancrap, nz2, ifail)
       !istatus = nf90_inq_dimid (ncid, 'z2', id)
       !if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='z2')
       !istatus = nf90_inquire_dimension (ncid, id, len=nz2)
       !if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=id)


!       id = ncdid (ncid, 'z1', ifail)
!       call ncdinq (ncid, id,fortrancrap, nz1, ifail)
!       istatus = nf90_inq_dimid (ncid, 'z1', id)
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='z1')
!       istatus = nf90_inquire_dimension (ncid, id, len=nz1)
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=id)
!
!!       id = ncdid (ncid, 'npsi', ifail)
!!       call ncdinq (ncid, id, fortrancrap, nr, ifail)
!       istatus = nf90_inq_dimid (ncid, 'npsi', id)
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='npsi')
!       istatus = nf90_inquire_dimension (ncid, id, len=nr)
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=id)

       !nz2 = npsi
       nr = npsi_chease
       
!       start(1) = 8
!       id = ncvid (ncid, 'nxy', ifail)
!       call ncvgt1 (ncid, id, start, iwork, ifail)
!       i_sym = iwork
!       istatus = nf90_inq_varid (ncid, 'nxy', id)
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='nxy')
!       istatus = nf90_get_var (ncid, id, i_sym, start=(/ 8 /))
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)

!     netcdf read scalar: nt
!
!     nt == number of theta grid points in theta eq grid
       nt = nchi_chease + 1
       B_norm_chease = b0exp_chease

     
      !nt = nchi_chease

!       id = ncdid (ncid, 'nthe', ifail)
!       call ncdinq (ncid, id, fortrancrap, nt, ifail)
!!       nt = iwork
!       istatus = nf90_inq_dimid (ncid, 'nthe', id)
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='nthe')
!       istatus = nf90_inquire_dimension (ncid, id, len=nt)
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=id)

       call alloc_arrays(nr, nt)

!     netcdf read vectors: rho_d, eqpsi, psi_bar, fp, beta, pressure
!
!     rho_d(1:nr) == half diameters of flux surfaces at elevation 
!                             of Z_mag on the radial grid
!     psi is the poloidal flux
!     eqpsi(1:nr) == values of psi on the radial grid
!     psi_bar(1:nr) == values of psi_bar on the radial grid
!     [psi_bar == (eqpsi - psi_0)/(psi_a - psi_0) if not available]
!     fp(1:nr) == the function that satisfies 
!              B = fp grad zeta + grad zeta x grad psi
!     beta(1:nr) == local beta, with the magnetic field defined 
!     to be vacuum magnetic field on axis
!     pressure(1:nr) == pressure profile on the radial grid,
!     normalized to the value at the magnetic axis.     

!       allocate(workr(nz2))
!       workr = 0.
!       start(1) = 1
!       cnt(1) = nz2

!       id = ncvid (ncid, 'rho', ifail)
!       call ncvgt (ncid, id, start, cnt, workr, ifail)
!       rho_b = workr(1:nr)
! rho_d must be defined by hand

!       id = ncvid (ncid, 'psivec', ifail)
!       call ncvgt (ncid, id, start, cnt, workr, ifail)
!       eqpsi = workr(1:nr)    !*1.e-8
!       istatus = nf90_inq_varid (ncid, 'psivec', id)
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='psivec')
!       istatus = nf90_get_var (ncid, id, eqpsi, count=(/ nr /))
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
       eqpsi = psi_chease

!       id = ncvid (ncid, 'fvec', ifail)
!       call ncvgt (ncid, id, start, cnt, workr, ifail)
!       fp = workr(1:nr)
!       istatus = nf90_inq_varid (ncid, 'fvec', id)
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='fvec')
!       istatus = nf90_get_var (ncid, id, fp, count=(/ nr /))
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
       fp = f_chease

!       id = ncvid (ncid, 'pvec', ifail)
!       call ncvgt (ncid, id, start, cnt, workr, ifail)       
!       pressure = workr(1:nr)
!       istatus = nf90_inq_varid (ncid, 'pvec', id)
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='pvec')
!       istatus = nf90_get_var (ncid, id, eqpsi, count=(/ nr /))
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
!
       qsf = q_chease
       pressure = p_chease
!     assign scalars: psi_0, psi_a, B_T
!
!     psi_0 == value of psi at the magnetic axis
!     psi_a == value of psi at the boundary of the plasma
!     B_T == vacuum toroidal magnetic field at R_center
!
       
       psi_0 = eqpsi(1)
       psi_a = eqpsi(nr)
       
       psi_bar = (eqpsi-psi_0)/(psi_a-psi_0)
       !write (*,*) 'psi_bar', psi_bar

!     netcdf read 2d field: R_psi,Z_psi and B_psi (mod(B))
!     R_psi(1:nr, 1:nt) -- transpose of Menard's storage ...
!     Z_psi(1:nr, 1:nt)
!     B_psi(1:nr, 1:nt)
!         

       !allocate(work(nz1*nz2))

!       start(1) = 1
!       start(2) = 1
!       cnt(1) = nz1
!       cnt(2) = nz2
!       id = ncvid (ncid, 'x', ifail)
!       call ncvgt (ncid, id, start, cnt, work, ifail)
!       istatus = nf90_inq_varid (ncid, 'x', id)
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='x')
!       istatus = nf90_get_var (ncid, id, work, count=(/ nz1*nz2 /))
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
!
!       if(i_sym == 0) then
!          write(*,*) 'Up-down asymmetries not yet handled correctly.'
!       endif

       do j=1,nt
          do i=1,nr
             !R_psi(i,j) = work(3+j-1+nz1*(i-1))
             !if j < nt/2 then 
             R_psi(i,j) = r_chease(i,chease_chi_index(nchi_chease, j))
             !R_psi(i,j) = r_chease(i,j)
            !else
          enddo
       enddo
       
!       cnt(1) = nz1
!       cnt(2) = nz2
!       id = ncvid (ncid, 'z', ifail)
!       call ncvgt (ncid, id, start, cnt, work, ifail)
!       istatus = nf90_inq_varid (ncid, 'z', id)
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='z')
!       istatus = nf90_get_var (ncid, id, work, count=(/ nz1, nz2 /))
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, id)
!
!       if(i_sym == 0) then
!          write(*,*) 'Up-down asymmetries not yet handled correctly.'
!       endif

       do j=1,nt
          do i=1,nr
             !Z_psi(i,j) = z_chease(i,j) * 10.0
             Z_psi(i,j) = z_chease(i,chease_chi_index(nchi_chease, j)) 
!             Z_psi(i,j) = work(3+j-1+nz1*(i-1))
          enddo
       enddo

!       cnt(1) = nz1
!       cnt(2) = nz2
!       id = ncvid (ncid, 'B_psi', ifail)
!       call ncvgt (ncid, id, start, cnt, work, ifail)
!       do j=1,nt
!          do i=1,nr
!             B_psi(i,j) = work(1+i-1+nr*(j-1))
!          enddo
!       enddo


!     assign scalars: R_mag, Z_mag, aminor
!

!     R_mag == R position of magnetic axis
!     Z_mag == Z position of magnetic axis
!     aminor    == half diameter of last flux surface at elevation of Z_mag

       !d = R_psi(nr,1) - R_psi(nr,nt)
       !d = R_psi(nr,(nt - 1)/2+1) - R_psi(nr,1)
       d = R_psi(nr,nt/2+1) - R_psi(nr,1)
       !d = R_psi(nr,1) - R_psi(nr,nt/2)
       !R_geo = (R_psi(nr,1) + R_psi(nr, nt))/2.
       !R_geo = (R_psi(nr,(nt-1)/2+1) + R_psi(nr, 1))/2.
       R_geo = (R_psi(nr,nt/2+1) + R_psi(nr, 1))/2.
       !R_geo = (R_psi(nr,1) + R_psi(nr, nt/2))/2.
       aminor = d/2.
       
       B_T = fp(nr)/R_geo
       !B_T = fp(1)/R_psi(1,1)
       !B_T = b0exp_chease

       R_psi = R_psi / aminor
       Z_psi = Z_psi / aminor
       
       R_mag = R_psi(1,1)
       Z_mag = Z_psi(1,1)

       !beta = 8. * pi * 1.e-7 * pressure / B_T**2  ! not correct definition yet
       !beta = 4. * pi * 1.e-7 * pressure / B_T**2  ! not correct definition yet
       beta = 8. * pi * 1.e-7 * pressure / b0exp_chease**2  ! not correct definition yet
       beta_0 = beta(1)
       pressure = pressure/pressure(1)

!       call ncclos (ncid, ifail)

!       deallocate(work,workr)
!       deallocate(work)

    endif   ! end of external reads
                       
!
!     Normalize, rename quantities 
!

    avgrmid = aminor
    ! rmaj changes depending on iflux.. FIX!!
    ! see geometr... for Miller is the rgeo for the flux surface of 
    ! interest... rgeo is rgeo of LCFS... check!!!
    ! actually when iflux = 1, it is the major radius of the magnetic axis,
    ! so below is OK.
    rmaj = R_mag 
    B_T0 = abs(B_T)

    psi_N = B_T0 * avgrmid**2
    psi_a = psi_a / psi_N
    psi_0 = psi_0 / psi_N
    psi_a_out = psi_a 
    psi_0_out = psi_0 
    eqpsi = eqpsi / psi_N

    f_N = B_T0*avgrmid
    fp=fp/f_N

    nthg=nt

    if (debug) then
       write (*,*) "Finished ceqin... imported CHEASE equilibrium"
       write (*,*) 'Some important quantities:'
       write (*,*) "aminor", aminor
       write (*,*) 'R_mag', R_mag
       write (*,*) 'B_T0', B_T0
       write (*,*) 'f_N', f_N
       write (*,*) 'nthg', nthg
       write (*,*) 'beta', beta_0
    end if
  end subroutine ceqin

  subroutine alloc_arrays(nr, nt)
    implicit none
    integer, intent(in) :: nr, nt

    allocate(rho_d(nr), eqpsi(nr), psi_bar(nr), fp(nr), beta(nr), pressure(nr), &
         rc(nr), diam(nr), qsf(nr), rho_b(nr))
    rho_d = 0. ; eqpsi = 0. ; psi_bar = 0. ; fp = 0. ; beta = 0. ; pressure = 0. 
    rc = 0. ; diam = 0. ; qsf = 0. ; rho_b = 0.
    allocate(R_psi(nr, nt), Z_psi(nr, nt))
    R_psi = 0.  ; Z_psi = 0.
    !    allocate(B_psi(nr, nt))
    allocate(drm(nr, nt, 2), dzm(nr, nt, 2), dbtm(nr, nt, 2), &
         dpm(nr, nt, 2), dtm(nr, nt, 2))
    drm = 0. ; dzm = 0. ; dbtm = 0. ; dpm = 0. ; dtm = 0.
!    allocate(dbm(nr, nt, 2), dbcart(nr, nt, 2), dbbish(nr, nt, 2))
    allocate(dpcart(nr, nt, 2), dtcart(nr, nt, 2), dbtcart(nr, nt, 2))
    allocate(dpbish(nr, nt, 2), dtbish(nr, nt, 2), dbtbish(nr, nt, 2))

    dpcart = 0. ; dtcart = 0. ; dbtcart = 0. 
    dpbish = 0. ; dtbish = 0. ; dbtbish = 0.

  end subroutine alloc_arrays

  subroutine dealloc_arrays
    implicit none
    if(allocated(rho_d)) deallocate(rho_d,eqpsi,psi_bar,fp,beta,pressure,rc,diam,qsf,rho_b)
    if(allocated(R_psi)) deallocate(R_psi,Z_psi)
    if(allocated(drm)) deallocate(drm,dzm,dbtm,dpm,dtm)
    if(allocated(dpcart)) deallocate(dpcart,dtcart,dbtcart)
    if(allocated(dpbish)) deallocate(dpbish,dtbish,dbtbish)
  end subroutine dealloc_arrays
  
  subroutine ceq_finish
    implicit none
    call dealloc_arrays
  end subroutine ceq_finish
  
  subroutine ceq_init
    use constants, only: pi
    implicit none
    real, dimension(nr,nt) :: eqpsi1, eqth, eqbtor
    integer i, j
   
    do j=1,nt
       do i=1,nr
          eqbtor(i,j) = fp(i)/R_psi(i,j)
          eqpsi1(i,j) = eqpsi(i)
       enddo
    enddo
    
    !if (transp) then
       do j=1,nt
          eqth(:,j) = (j-1)*2.*pi/float(nt-1)-pi
       enddo
    !else
       !do j=1,nt
          !eqth(:,j) = (j-1)*pi/float(nt-1)
       !enddo
    !end if

    call derm(eqth,   dtm,  'T')
    call derm(R_psi,  drm,  'E')
    call derm(Z_psi,  dzm,  'O')
!    call derm(B_psi,  dbm,  'E')
    call derm(eqbtor, dbtm, 'E')
    call derm(eqpsi1, dpm,  'E')
    

! diagnostics
!      do j=1,nt
!         do i=1,nr
!            write(*,*) i,j
!            write(*,100) drm(i,j,1),drm(i,j,2),R_psi(i,j)
!            write(*,101) dzm(i,j,1),dzm(i,j,2),Z_psi(i,j)
!            write(*,102) dzm(i,j,1),dtm(i,j,2),eqth(i,j)
!         enddo
!      enddo
! 100  format('(gr R)1 ',g10.4,' (gr R)2 ',g10.4,' R ',g10.4)
! 101  format('(gr Z)1 ',g10.4,' (gr Z)2 ',g10.4,' Z ',g10.4)
! 102  format('(gr t)1 ',g10.4,' (gr t)2 ',g10.4,' t ',g10.4)
!      write(*,*) nr, nt
!      stop

! grad(psi) in cartesian form 
    call eqdcart(dpm, dpcart)
! grad(psi) in Bishop form 
    call eqdbish(dpcart, dpbish)

! grad(B) in cartesian form
!    call eqdcart(dbm, dbcart)
! grad(B) in Bishop form
!    call eqdbish(dbcart, dbbish)

! grad(BT) in cartesian form
    call eqdcart(dbtm, dbtcart)
! grad(BT) in Bishop form
    call eqdbish(dbtcart, dbtbish)

! grad(theta) in cartesian form
    call eqdcart(dtm, dtcart)
! grad(theta) in Bishop form
    call eqdbish(dtcart, dtbish)

! diagnostics
!      call inter_cspl(nr, eqpsi,dpcart(1,1,1),1,rp,f)
!      write(*,*) f
!      call inter_cspl(nr, eqpsi,dpcart(1,1,2),1,rp,f)
!      write(*,*) f

  end subroutine ceq_init

  subroutine derm(f, dfm, char)

    use constants, only: pi
    implicit none
    integer i, j
    character(1), intent(in) :: char
    real, intent(in) :: f(:,:)
    real, intent(out) :: dfm(:,:,:)
    
    i=1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))         
    
    i=nr
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))
   
    !if (.not. transp) then
!! assume up-down symmetry for now:
       
       !select case (char)
       !case ('E') 
          !j=1
          !dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,j+1))
          
          !j=nt      
          !dfm(:,j,2) = -0.5*(f(:,j-1)-f(:,j-1))
       !case ('O')
          !j=1
          !dfm(:,j,2) = 0.5*(f(:,j+1)+f(:,j+1))
          
          !j=nt      
          !dfm(:,j,2) = -0.5*(f(:,j-1)+f(:,j-1))
       !case ('T')
          !j=1
          !dfm(:,j,2) = f(:,j+1)
          
          !j=nt      
          !dfm(:,j,2) = pi - f(:,j-1)
       !end select
       
    !else

       select case (char)
       case ('E') 
          j=1
          dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,nt-1))
          
          j=nt      
          dfm(:,j,2) = 0.5*(f(:,2)-f(:,j-1))
       case ('O')
          j=1
          dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,nt-1))
          
          j=nt      
          dfm(:,j,2) = 0.5*(f(:,2)-f(:,j-1))
       case ('T')
          j=1
          dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,nt-1)+2.*pi)
          
          j=nt
          dfm(:,j,2) = 0.5*(2.*pi + f(:,2) - f(:,j-1))
       end select
       
    !end if
    
    do i=2,nr-1
       dfm(i,:,1)=0.5*(f(i+1,:)-f(i-1,:))
    enddo
    
    do j=2,nt-1
       dfm(:,j,2)=0.5*(f(:,j+1)-f(:,j-1))
    enddo

!    if (char=='T') then
!       do i=1,nr
!          write(*,*) i,'1'
!          write(*,*) dfm(i,:,1)
!          write(*,*) i,'2'
!          write(*,*) dfm(i,:,2)
!          write(*,*)
!       enddo
!    end if

  end subroutine derm

  subroutine gradient(rgrid, theta, grad, char, rp, nth_used, ntm)

    use splines, only: inter_d_cspl
    implicit none
    
    integer, intent(in) :: nth_used, ntm
    character(1), intent(in) :: char
    real, intent(in) :: rgrid(-ntm:), theta(-ntm:)
    real, intent(out) :: grad(-ntm:,:)
    real, intent(in) :: rp
    real :: tmp(2), aa(1), daa(1), rpt(1)
    real, dimension(nr,nt,2) :: dcart
    real, dimension(nr,nt) :: f
    integer :: i
    
    select case(char)
    case('B') 
       dcart = dbcart
       !write(*,*) 'error: bishop = 1 not allowed with ceq.'; stop
    case('D')  ! diagnostic 
       dcart = dbtcart
    case('P') 
       dcart = dpcart
    case('R') 
       dcart = dpcart  ! dpcart is correct for 'R'
    case('T')
       dcart = dtcart
    end select
    
    do i=-nth_used,-1
       f(:,:) = dcart(:,:,1)
       call eqitem(rgrid(i), theta(i), f, tmp(1))
       f(:,:) = dcart(:,:,2)
       call eqitem(rgrid(i), theta(i), f, tmp(2))
       !if(char == 'T' .and. .not. transp) then
          !grad(i,1)=-tmp(1)
          !grad(i,2)=-tmp(2)
       !else
          grad(i,1)=tmp(1)
          grad(i,2)=tmp(2)
       !endif
    enddo

    do i=0,nth_used
       f(:,:) = dcart(:,:,1)
       call eqitem(rgrid(i),theta(i),f,tmp(1))
       f(:,:) = dcart(:,:,2)
       call eqitem(rgrid(i),theta(i),f,tmp(2))
       grad(i,1)=tmp(1)
       grad(i,2)=tmp(2)
    enddo

!     to get grad(pressure), multiply grad(psi) by dpressure/dpsi

    if(char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(nr, eqpsi, pressure, 1, rpt, aa, daa)
       do i=-nth_used, nth_used
          grad(i,1)=grad(i,1)*daa(1) * 0.5*beta_0
          grad(i,2)=grad(i,2)*daa(1) * 0.5*beta_0
       enddo
    endif

  end subroutine gradient

  subroutine bgradient(rgrid, theta, grad, char, rp, nth_used, ntm)
    use mp, only: mp_abort
    use splines, only: inter_d_cspl
    implicit none
    
    integer, intent(in) :: nth_used, ntm
    character(1), intent(in) :: char
    real, intent(in) :: rgrid(-ntm:), theta(-ntm:), rp
    real, intent(out) :: grad(-ntm:,:)
    real :: aa(1), daa(1), rpt(1)
    real, dimension(nr,nt,2) ::  dbish
    integer :: i

    select case(char)
    case('B') 
!       dbish = dbbish
       call mp_abort('error: bishop = 1 not allowed with ceq. (2)',.true.)
    case('D')  ! diagnostic
       dbish = dbtbish
    case('P') 
       dbish = dpbish
    case('R') 
       dbish = dpbish  ! dpcart is correct for 'R'
    case('T')
       dbish = dtbish
    end select

    do i=-nth_used,nth_used
       call eqitem(rgrid(i), theta(i), dbish(:,:,1), grad(i,1))
       call eqitem(rgrid(i), theta(i), dbish(:,:,2), grad(i,2))
    enddo

    !if (char == 'T' .and. .not. transp) then
       !where (theta(-nth_used:nth_used) < 0.0)
          !grad(-nth_used:nth_used,1) = -grad(-nth_used:nth_used,1)
          !grad(-nth_used:nth_used,2) = -grad(-nth_used:nth_used,2)
       !end where
    !end if

!     to get grad(pressure), multiply grad(psi) by dpressure/dpsi

    if(char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(nr, eqpsi, pressure, 1, rpt, aa, daa)
       do i=-nth_used, nth_used
          grad(i,1)=grad(i,1)*daa(1) * 0.5*beta_0
          grad(i,2)=grad(i,2)*daa(1) * 0.5*beta_0
       enddo
    endif

  end subroutine bgradient

  subroutine eqitem(r, theta_in, f, fstar)
    use mp, only: mp_abort
    use constants, only: pi
    real, intent(in) :: r, theta_in
    real, intent(out) :: fstar
    real, dimension(:,:), intent(inout) :: f
    integer :: i, j, istar, jstar
    real :: thet, sign
    real :: st, dt, sr, dr
    real, dimension(size(f,2)) :: mtheta

! check for axis evaluation
      
    if(r == eqpsi(1)) then
       write(*,*) 'no evaluation at axis allowed in eqitem'
       write(*,*) r, theta_in, eqpsi(1)
       call mp_abort('no evaluation at axis allowed in eqitem')
    endif
    
! allow psi(r) to be a decreasing function

    sign=1.
    if(eqpsi(2) < eqpsi(1)) sign=-1.
    
    if(r < sign*eqpsi(1)) then
       write(*,*) 'r < Psi_0 in eqitem'
       write(*,*) r,sign,eqpsi(1)
       call mp_abort('r < Psi_0 in eqitem')
    endif
      
! find r on psi mesh

! disallow evaluations on or outside the surface for now

    if(r >= eqpsi(nr)) then
       write(*,*) 'No evaluation of eqitem allowed on or outside surface'
       write(*,*) '(Could this be relaxed a bit?)'
       write(*,*) r, theta_in, eqpsi(nr), sign
       call mp_abort('No evaluation of eqitem allowed on or outside surface')
    endif
    
    istar=0
    do i=2,nr
       if(r < sign*eqpsi(i)) then
          dr = r - sign*eqpsi(i-1)
          sr = sign*eqpsi(i) - r
          istar=i-1
          exit
       endif
    enddo
    
! no gradients available at axis, so do not get too close

    if(istar == 1) then
       write(*,*) 'Too close to axis in eqitem'
       write(*,*) r, theta_in, eqpsi(1), eqpsi(2)
       call mp_abort('Too close to axis in eqitem')
    endif
  
! Now do theta direction

    thet = mod2pi(theta_in)

!    write(*,*) 'theta_in, thet = ',theta_in, thet

! assume up-down symmetry

    !if (.not. transp) then
       !tp=abs(thet)
       !tps=1.
       !if(char == 'Z' .and. tp >= 1.e-10) tps=thet/tp

!! get thet on theta mesh

       !do j=1,nt
          !mtheta(j)=(j-1)*pi/float(nt-1)
       !enddo
       
!! note that theta(1)=0 for gen_eq theta 

       !jstar=-1
       !do j=1,nt
          !if(jstar /= -1) cycle
          !if(tp < mtheta(j)) then
             !dt = tp - mtheta(j-1)
             !st = mtheta(j) - tp
             !jstar=j-1
          !endif
       !enddo
      
!! treat theta = pi separately
  
       !if(jstar == -1) then
          !jstar=nt-1
          !dt=mtheta(jstar+1)-mtheta(jstar)
          !st=0.
       !endif

!! use opposite area stencil to interpolate

!!    if(char == 'R') i=1
!!    if(char == 'Z') i=2
       !fstar=f(istar    , jstar    ) * sr * st &
            !+f(istar + 1, jstar    ) * dr * st &
            !+f(istar    , jstar + 1) * sr * dt &
            !+f(istar + 1, jstar + 1) * dr * dt
       !fstar=fstar*tps &
            !/abs(eqpsi(istar+1)-eqpsi(istar)) &
            !/(mtheta(jstar+1)-mtheta(jstar))

    !else

! get thet on theta mesh

       do j=1,nt
          mtheta(j)=(j-1)*2.*pi/float(nt-1)-pi
       enddo
       
! note that theta(1)=0 for gen_eq theta 

       jstar=-1
       do j=1,nt
          if(jstar /= -1) cycle
          if(thet < mtheta(j)) then
             dt = thet - mtheta(j-1)
             st = mtheta(j) - thet
             jstar=j-1
          endif
       enddo
      
! treat pi separately
  
       if(jstar == -1) then
          jstar=nt-1
          dt=mtheta(jstar+1)-mtheta(jstar)
          st=0.
       endif

! use opposite area stencil to interpolate

!    if(char == 'R') i=1
!    if(char == 'Z') i=2
       fstar=f(istar    , jstar    ) * sr * st &
            +f(istar + 1, jstar    ) * dr * st &
            +f(istar    , jstar + 1) * sr * dt &
            +f(istar + 1, jstar + 1) * dr * dt
       fstar=fstar &
            /abs(eqpsi(istar+1)-eqpsi(istar)) &
            /(mtheta(jstar+1)-mtheta(jstar))
    !end if
!     write(*,*) i, dr, dt, sr, st, istar, jstar
!     write(*,*) f(istar,jstar+1),f(istar+1,jstar+1)
!     write(*,*) f(istar,jstar),f(istar+1,jstar)
!     write(*,*) eqpsi(istar),eqpsi(istar+1)
!     write(*,*) mtheta(jstar),mtheta(jstar+1)
!     write(*,*) 

  end subroutine eqitem

  subroutine eqdcart(dfm, dfcart)
      
    implicit none

    real, dimension (:,:,:), intent(in)  :: dfm
    real, dimension (:,:,:), intent(out) :: dfcart
    real, dimension (size(dfm,1),size(dfm,2)) :: denom
    integer :: i, j
      
    denom(:,:) = drm(:,:,1)*dzm(:,:,2) - drm(:,:,2)*dzm(:,:,1)

!!!
!    dfcart = 0.
    
    dfcart(:,:,1) =   dfm(:,:,1)*dzm(:,:,2) - dzm(:,:,1)*dfm(:,:,2)
    dfcart(:,:,2) = - dfm(:,:,1)*drm(:,:,2) + drm(:,:,1)*dfm(:,:,2)

    do j=1,nt
       do i=2,nr
          dfcart(i,j,:)=dfcart(i,j,:)/denom(i,j)
       enddo
    enddo    

  end subroutine eqdcart

  subroutine eqdbish(dcart, dbish)
    implicit none
    real, dimension(:, :, :), intent (in) :: dcart
    real, dimension(:, :, :), intent(out) :: dbish
    real, dimension(size(dcart,1),size(dcart,2)) :: denom
    integer :: i, j

    denom(:,:) = sqrt(dpcart(:,:,1)**2 + dpcart(:,:,2)**2)

    dbish(:,:,1) = dcart(:,:,1)*dpcart(:,:,1) + dcart(:,:,2)*dpcart(:,:,2)
    dbish(:,:,2) =-dcart(:,:,1)*dpcart(:,:,2) + dcart(:,:,2)*dpcart(:,:,1)
    
    do j=1,nt
       do i=2,nr
          dbish(i,j,:) = dbish(i,j,:)/denom(i,j)
       enddo
    enddo

  end subroutine eqdbish

  function initialize_invR (init) 

    integer :: init, initialize_invR
    
    init_invR = .false.
    if(init == 1) init_invR = .true.
    initialize_invR = 1

  end function initialize_invR

  function invR (r, theta)
   
    real, intent (in) :: r, theta
    real :: f, invR
    real :: th
    
    th = mod2pi( theta)
    
    call eqitem(r, th, R_psi, f)
    !if (debug) write(*,*) 'invR', 'th', th, 'R', f
    invR=1./f
    
  end function invR

  function Rpos (r, theta)
   
    real, intent (in) :: r, theta
    real :: f, Rpos
    real :: th
    
    th = mod2pi( theta)
    
    call eqitem(r, th, R_psi, f)
    Rpos=f
    
  end function Rpos

  function Zpos (r, theta)
   
    real, intent (in) :: r, theta
    real :: f, Zpos
    real :: th
    
    th = mod2pi( theta)
    
    call eqitem(r, th, Z_psi, f)
    Zpos=f
    
  end function Zpos

  function initialize_psi (init) 

    integer, intent(in) :: init
    integer :: initialize_psi
    
    init_psi = .false.
    if(init == 1) init_psi = .true.
    initialize_psi = 1

  end function initialize_psi

  function psi (r)
    real, intent (in) :: r
    real :: psi

    psi = r
    
  end function psi

  function mod2pi (theta)
    use constants, only: pi
    real, intent(in) :: theta
    real :: th, mod2pi
    logical :: out
    
    if(theta <= pi .and. theta >= -pi) then
       mod2pi = theta
       return
    endif
    
    th=theta
    out=.true.
    do while(out)
       if(th > pi) th = th - 2.*pi
       if(th <-pi) th = th + 2.*pi
       if(th <= pi .and. th >= -pi) out=.false.
    enddo
    mod2pi=th
  end function mod2pi
   
  function initialize_diameter (init) 
    integer, intent(in) :: init
    integer :: initialize_diameter
    
    init_diameter = .false.
    if(init == 1) init_diameter = .true.
    initialize_diameter = 1

  end function initialize_diameter

  function diameter (rp)
    use splines, only: new_spline, splint, spline
    real, intent(in) :: rp
    real :: diameter
    type (spline), save :: spl

    if(init_diameter) then
       diam(:) = abs(R_psi(:,nt/2) - R_psi(:,1))
       call new_spline(nr, eqpsi, diam, spl)
       init_diameter = .false.
    endif

    diameter = splint(rp, spl)

  end function diameter

  function initialize_rcenter (init) 
    integer, intent(in) :: init
    integer :: initialize_rcenter
    
    init_rcenter = .false.
    if(init == 1) init_rcenter = .true.
    initialize_rcenter = 1

  end function initialize_rcenter

  function rcenter (rp)
    use splines, only: new_spline, splint, spline
    real, intent(in) :: rp
    real :: rcenter
    type (spline), save :: spl

    if(init_rcenter) then
       rc(:) = 0.5*(R_psi(:,1) + R_psi(:,nt/2))
       call new_spline(nr, eqpsi, rc, spl)
       init_rcenter = .false.
    endif

    rcenter = splint(rp, spl)

  end function rcenter

  function initialize_dbtori (init) 

    integer, intent(in) :: init
    integer :: initialize_dbtori
    
    init_dbtori = .false.
    if(init == 1) init_dbtori = .true.
    initialize_dbtori = 1

  end function initialize_dbtori

  function dbtori (pbar)
    use splines, only: new_spline, dsplint, spline
    real, intent(in) :: pbar
    real :: dbtori
    type (spline), save :: spl

    if(init_dbtori) call new_spline(nr, psi_bar, fp, spl)
    init_dbtori=.false.

    dbtori = dsplint(pbar, spl)/(psi_a-psi_0)

  end function dbtori

  function initialize_btori (init) 
    integer, intent(in) :: init
    integer :: initialize_btori
    
    init_btori = .false.
    if(init == 1) init_btori = .true.
    initialize_btori = 1
  end function initialize_btori

  function btori (pbar)
    use splines, only: new_spline, splint, spline
    real, intent(in) :: pbar
    real :: btori
    type (spline), save :: spl

    if(init_btori) call new_spline(nr, psi_bar, fp, spl)
    init_btori=.false.

    btori = splint(pbar, spl)
  end function btori

  function initialize_q (init) 

    integer, intent(in) :: init
    integer :: initialize_q
    
    init_q = .false.
    if(init == 1) init_q = .true.
    initialize_q = 1

  end function initialize_q

  function qfun (pbar)
  
    use splines, only: new_spline, splint, spline
    real, intent(in) :: pbar
    real :: qfun
    type (spline), save :: spl

    if(init_q) call new_spline(nr, psi_bar, qsf, spl)
    init_q = .false.

    qfun = splint(pbar, spl)

  end function qfun

  function initialize_pressure (init) 

    integer, intent(in) :: init
    integer :: initialize_pressure
    
    init_pressure = .false.
    if(init == 1) init_pressure = .true.
    initialize_pressure = 1

  end function initialize_pressure

  function pfun (pbar)
  
    use splines, only: new_spline, splint, spline
    real, intent(in) :: pbar
    real :: pfun
    type (spline), save :: spl

    if(init_pressure) call new_spline(nr, psi_bar, pressure, spl)
    init_pressure = .false.
!
! p_N would be B**2/mu_0 => p = beta/2 in our units
!
    pfun = 0.5*beta_0*splint(pbar, spl)

  end function pfun
  
  function initialize_dpressure (init) 

    integer, intent(in) :: init
    integer :: initialize_dpressure
    
    init_dpressure = .false.
    if(init == 1) init_dpressure = .true.
    initialize_dpressure = 1

  end function initialize_dpressure

  function dpfun (pbar)
  
    use splines, only: new_spline, dsplint, spline
    real, intent(in) :: pbar
    real :: dpfun
    type (spline), save :: spl
!
! p_N would be B**2/mu_0 => p = beta/2 in our units
!
    if(init_dpressure) then
       call new_spline(nr, psi_bar, pressure, spl)
       init_dpressure = .false.
    endif

    dpfun = dsplint(pbar, spl)/(psi_a-psi_0) * beta_0/2.

  end function dpfun

  function initialize_beta (init) 

    integer, intent(in) :: init
    integer :: initialize_beta
    
    init_beta = .false.
    if(init == 1) init_beta = .true.
    initialize_beta = 1

  end function initialize_beta

  function betafun (pbar)
  
    use splines, only: new_spline, splint, spline
    real, intent(in) :: pbar
    real :: betafun
    type (spline), save :: spl

    if(pbar == 0.) then
       betafun=beta(1)
       return
    endif

    if(init_beta) call new_spline(nr, psi_bar, beta, spl)
    init_beta = .false.

    betafun = splint(pbar, spl)

  end function betafun

  subroutine Hahm_Burrell(irho, a) 

    real, intent(in) :: a
    integer, intent(in) :: irho
    integer :: i
    real :: gradpsi, mag_B, rho_eq, rp1, rp2, rho1, rho2, drhodpsiq
    real, dimension(nr) :: gamma, pbar, dp, d2p, pres

    gamma = 0.

    pbar = (eqpsi-eqpsi(1))/(eqpsi(nr)-eqpsi(1))

    do i=2, nr-1
       dp(i)  = dpfun(pbar(i))
    enddo

    do i=3,nr-2
       d2p(i) = (dp(i+1)-dp(i-1))/(eqpsi(i+1)-eqpsi(i-1))
    enddo

    pres=0.
    do i=3, nr-2
       rp1=eqpsi(i+1)
       rp2=eqpsi(i-1)
       rho1=0.5*diameter(rp1)
       rho2=0.5*diameter(rp2)
       drhodpsiq=(rho1-rho2)/(rp1-rp2)
       
       call eqitem(eqpsi(i), 0., dpbish(:,:,1), gradpsi)
       
       mag_B=sqrt( (btori(pbar(i))**2 + gradpsi**2))/Rpos(eqpsi(i),0.)

       pres(i) = pfun(pbar(i))

       gamma(i) = 0.01*gradpsi**2*(d2p(i)/pres(i)-a*(dp(i)/pres(i))**2) &
            /mag_B*(2.*pres(i)/beta_0)**((1-a)/2.) &
            *(-pres(i)/(dp(i)/drhodpsiq))            
    enddo
    
    do i=3,nr-2
       if(irho == 1) then
          rho_eq = eqpsi(i)
       else if(irho == 2) then
          rho_eq = 0.5 * diameter(eqpsi(i))
       else if(irho == 3) then
          rho_eq = pbar(i)
       endif
       write(24,1000) i, pbar(i), rho_eq, pres(i), gamma(i)
    enddo

!    write(*,*) 1/(psi_a-psi_0)
!    if(irho == 1) write(* ,*) '# irho = 1 produces psi instead of rho_eq'
    if(irho == 1) write(24,*) '# irho = 1 produces psi instead of rho_eq'

1000 format(i5,11(1x,e16.9))

  end subroutine Hahm_Burrell

end module ceq
