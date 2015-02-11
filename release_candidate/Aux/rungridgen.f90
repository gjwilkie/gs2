program rungridgen
  use file_utils
  use gridgen4mod
  use constants, only: pi, twopi
  implicit none

  ! input parameters
  character(200) :: source, gingrid, gsource
  integer :: nthetaout, nlambdaout, nperiodout
  integer :: npadd, iperiod
  real :: alknob, epsknob, extrknob, bpknob, smoothknob
  integer :: nfinegrid
  real :: thetamax, deltaw, widthw, tension
  logical :: screenout
  logical :: auto_width, three_dim
  real :: cv_fraction, delth_max
  integer :: max_autoiter
  logical :: list = .false.


  ! work variables
  integer :: ntheta, nlambda, ntgrid
  integer :: ntgridin, nperiodin, nthetain
  real, dimension (:), allocatable :: thetain, bmagin, bmagsm
  real, dimension (:), allocatable :: thetaout, bmagout, alambdaout

  real :: drhodpsi, rmaj, shat, kxfac, qval
  real, dimension (:), allocatable :: thetagrid
  real, dimension (:), allocatable :: gbdrift, gradpar, grho
  real, dimension (:), allocatable :: cvdrift, gds2, bmag
  real, dimension (:), allocatable :: gds21, gds22
  real, dimension (:), allocatable :: cvdrift0, gbdrift0
  real, dimension (:), allocatable :: Rplot, Rprime
  real, dimension (:), allocatable :: Zplot, Zprime
  real, dimension (:), allocatable :: aplot, aprime
  
  call init_file_utils (list, name="grid")
  call read_parameters
  if (three_dim) then
     call allocate_arrays_3d
     call get_initial_grids_3d
  else
     call allocate_arrays
     call get_initial_grids
  end if
  call do_smoothing
  call generate_grids
  call write_output
  call finish_file_utils
  stop

contains

  subroutine read_parameters
    implicit none
    namelist /testgridgen/ source, gsource, &
         nthetaout,nlambdaout,nperiodout, &
         npadd,alknob,epsknob,bpknob,extrknob,smoothknob, nfinegrid, &
         thetamax, deltaw, widthw, tension, gingrid, screenout, &
         auto_width, cv_fraction, delth_max, max_autoiter, three_dim, &
         iperiod
    gsource = "eik6.out"
    source = "eik.out"
    nthetaout = 32
    nlambdaout = 20
    nperiodout = 2
    npadd = 2
    alknob = 0.1
    bpknob = 1.e-8
    epsknob = 1e-9
    extrknob = 10.0
    smoothknob = 0.0
    nfinegrid=200
    thetamax = 0.0
    deltaw = 0.0
    widthw = 1.0
    tension = 1.0
    screenout = .false.
    auto_width = .false.
    cv_fraction = 0.6
    delth_max = 0.5
    max_autoiter = 3
    gingrid = "gingrid"
    three_dim = .false.
    iperiod = 1
    read (unit=input_unit("testgridgen"), nml=testgridgen)
    nthetaout = nthetaout + 1
  end subroutine read_parameters

  subroutine allocate_arrays
    implicit none
    integer :: unit
    character(200) :: line

    call get_unused_unit(unit)
    open (unit=unit, file=trim(source), status="old")
    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt=*,err=10) ntgridin, nperiodin, nthetain, drhodpsi, rmaj, shat, kxfac, qval
10  continue
    close (unit=unit)
    
    allocate (thetain(nthetain+1),bmagin(nthetain+1),bmagsm(nthetain+1))
    allocate (thetaout(nthetaout),bmagout(nthetaout),alambdaout(nlambdaout))
    allocate (thetagrid(-ntgridin:ntgridin))
    allocate (gbdrift(-ntgridin:ntgridin))
    allocate (gradpar(-ntgridin:ntgridin))
    allocate (grho(-ntgridin:ntgridin))
    allocate (cvdrift(-ntgridin:ntgridin))
    allocate (gds2(-ntgridin:ntgridin))
    allocate (bmag(-ntgridin:ntgridin))
    allocate (gds21(-ntgridin:ntgridin))
    allocate (gds22(-ntgridin:ntgridin))
    allocate (cvdrift0(-ntgridin:ntgridin))
    allocate (gbdrift0(-ntgridin:ntgridin))
    allocate (Rplot(-ntgridin:ntgridin))
    allocate (Rprime(-ntgridin:ntgridin))
    allocate (Zplot(-ntgridin:ntgridin))
    allocate (Zprime(-ntgridin:ntgridin))
    allocate (aplot(-ntgridin:ntgridin))
    allocate (aprime(-ntgridin:ntgridin))
  end subroutine allocate_arrays

  subroutine allocate_arrays_3d
    implicit none
    integer :: unit, ntor
!    character(200) :: line

    call get_unused_unit(unit)
    open (unit=unit, file=trim(source), status="old")
    read (unit=unit, fmt=*) ntgridin, nthetain, ntor
    close (unit=unit)
    
    drhodpsi = 1.0  ! assumes our radial variable matches VVBAL

    nthetain = nthetain * iperiod

    allocate (thetain(nthetain+1),bmagin(nthetain+1),bmagsm(nthetain+1))
    allocate (thetaout(nthetaout),bmagout(nthetaout),alambdaout(nlambdaout))
    allocate (thetagrid(-ntgridin:ntgridin))
    allocate (gbdrift(-ntgridin:ntgridin))
    allocate (gradpar(-ntgridin:ntgridin))
    allocate (grho(-ntgridin:ntgridin))
    allocate (cvdrift(-ntgridin:ntgridin))
    allocate (gds2(-ntgridin:ntgridin))
    allocate (bmag(-ntgridin:ntgridin))
    allocate (gds21(-ntgridin:ntgridin))
    allocate (gds22(-ntgridin:ntgridin))
    allocate (cvdrift0(-ntgridin:ntgridin))
    allocate (gbdrift0(-ntgridin:ntgridin))
    allocate (Rplot(-ntgridin:ntgridin))
    allocate (Rprime(-ntgridin:ntgridin))
    allocate (Zplot(-ntgridin:ntgridin))
    allocate (Zprime(-ntgridin:ntgridin))
    allocate (aplot(-ntgridin:ntgridin))
    allocate (aprime(-ntgridin:ntgridin))

    cvdrift0 = 0.   ! assumes theta0 = 0.
    gbdrift0 = 0.   ! assumes theta0 = 0.
    gds21 = 0.      ! assumes theta0 = 0.
    gds22 = 0.      ! assumes theta0 = 0.
    grho  = 1.      ! assumes linear calculation
    Rplot = 1.      ! pure fiction
    Rprime = 1.     ! pure fiction
    Zplot = 1.      ! pure fiction
    Zprime = 1.     ! pure fiction
    aplot = 1.      ! pure fiction
    aprime = 1.     ! pure fiction
    
  end subroutine allocate_arrays_3d

  subroutine get_initial_grids

    use file_utils, only: get_unused_unit
    implicit none
    integer :: unit
    integer :: i
!    real :: discard
    character(200) :: line

    call get_unused_unit(unit)
    open (unit=unit, file=trim(source), status="old")
    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt="(a)") line

    ! gbdrift gradpar grho theta
    read (unit=unit, fmt="(a)") line
    do i = -ntgridin,ntgridin
       read (unit=unit, fmt=*) gbdrift(i),gradpar(i),grho(i),thetagrid(i)
    end do

    ! cvdrift gds2 bmag theta
    read (unit=unit, fmt="(a)") line
    do i = -ntgridin,ntgridin
       read (unit=unit, fmt=*) cvdrift(i),gds2(i),bmag(i),thetagrid(i)
    end do

    ! gds21 gds22 theta
    read (unit=unit, fmt="(a)") line
    do i = -ntgridin,ntgridin
       read (unit=unit, fmt=*) gds21(i),gds22(i),thetagrid(i)
    end do

    ! cvdrift0 gbdrift0 theta
    read (unit=unit, fmt="(a)") line
    do i = -ntgridin,ntgridin
       read (unit=unit, fmt=*) cvdrift0(i),gbdrift0(i),thetagrid(i)
    end do

    close (unit=unit)

    thetain = thetagrid(-nthetain/2:nthetain/2)
    bmagin = bmag(-nthetain/2:nthetain/2)
  end subroutine get_initial_grids

  subroutine get_initial_grids_3d
    use file_utils, only: get_unused_unit
    implicit none
    integer :: unit
    integer :: i, ntor
    real :: dpdpsi, pres, bavg, cvavg

    call get_unused_unit(unit)
    open (unit=unit, file=trim(source), status="old", action="read")

    read (unit=unit, fmt=*) ntgridin, nthetain, ntor

    nthetain = nthetain * iperiod
    do i = -ntgridin,ntgridin
       read (unit=unit, fmt=*) thetagrid(i),bmag(i),gradpar(i),gds2(i),cvdrift(i),gbdrift(i),dpdpsi,pres
    end do

    close (unit=unit)

    bavg = sum(bmag(-nthetain/2:nthetain/2))/(nthetain+1)
    cvavg = sum(cvdrift)/size(cvdrift)
    
    write(*,*) 'bavg = ',bavg
    write(*,*) 'cvavg = ',cvavg,' and remember that cvdrift assumed beta_prime=0'

!    bmag = bmag / bavg

    gds2 = gds2 * ntor**2 
    gbdrift = gbdrift * ntor * 2. / bmag**2
!    cvdrift = cvdrift * ntor * 2. / bmag**3    ! not correct for beta_prime term, so use: 
    cvdrift = gbdrift 

!    alp = -1./pres*dpdpsi

!
! these are theta and mod(b) grids that are periodic
!
    thetain = thetagrid(-nthetain/2:nthetain/2)
    bmagin = bmag(-nthetain/2:nthetain/2)
  end subroutine get_initial_grids_3d

  subroutine do_smoothing
    implicit none
    real :: var
    integer :: ifail
    if (smoothknob == 0.0) then
       bmagsm = bmagin
    else
       var = smoothknob
       ifail = 0
       call smooth (nthetain+1,thetain,bmagin,var,bmagsm,ifail)
       if (ifail /= 0) then
          print *, "smooth failed: ",ifail
          select case (ifail)
          case (129)
             print *, "IC IS LESS THAN N-1."
          case (130)
             print *, "N IS lESS THAN 3."
          case (131)
             print *, "INPUT ABSCISSAE ARE NOT ORDERED SO THAT X(I).LT.X(I+1)."
          case (132)
             print *, "DF(I) IS NOT POSITIVE FOR SOME I."
          case (133)
             print *, "JOB IS NOT 0 OR 1."
          case default
          end select
          stop
       end if
    end if
  end subroutine do_smoothing

  subroutine generate_grids
    implicit none
    integer :: i, iter, nin
    real :: d, deltaw_ok

    if (.not. auto_width) then
       ntheta = nthetaout
       nlambda = nlambdaout
       call gridgen4_1 (iperiod,nthetain+1,thetain,bmagsm, npadd, &
            alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw,tension, &
            ntheta,nlambda,thetaout,bmagout,alambdaout)
       return
    end if

    widthw = pi
    do i = 0, nthetain/2
       if (cvdrift(i) < 0.0) then
          widthw = thetagrid(i)
          exit
       end if
    end do

    print *, "widthw: ", widthw
    deltaw_ok = 0.0
    do iter = 1, max_autoiter
       ntheta = nthetaout
       nlambda = nlambdaout
       call gridgen4_1 (iperiod,nthetain+1,thetain,bmagsm, npadd, &
            alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw,tension, &
            ntheta,nlambda,thetaout,bmagout,alambdaout)
       print *, "iter: ", iter
       d = maxval(thetaout(2:ntheta) - thetaout(1:ntheta-1))
       print *, "max deltheta: ", d
       nin = count(abs(thetaout(1:ntheta) - thetamax) < widthw &
                   .or. abs(thetaout(1:ntheta) + thetamax) < widthw)
       print *, "count in +ve cvdrift: ", nin
       print *, "fraction in +ve cvdrift: ", real(nin)/real(ntheta)
       print *, "deltaw: ", deltaw
       if (d > delth_max) then
          if (deltaw_ok /= 0.0) then
             deltaw = sqrt(deltaw*deltaw_ok)
          else
             deltaw = 0.5*deltaw
          end if
          print *, "max deltheta too large: ", d
          if (deltaw < deltaw_ok) exit
       else if (real(nin)/real(ntheta) < cv_fraction) then
          print *, "fraction in +ve cvdrift too small: ", &
               real(nin)/real(ntheta)
          deltaw_ok = deltaw
          deltaw = deltaw*(cv_fraction/(real(nin)/real(ntheta)))**2
       end if
    end do

    if (deltaw_ok /= 0.0 .and. d > delth_max) then
       deltaw = deltaw_ok
       ntheta = nthetaout
       nlambda = nlambdaout
       call gridgen4_1 (iperiod,nthetain+1,thetain,bmagsm, npadd, &
            alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw,tension, &
            ntheta,nlambda,thetaout,bmagout,alambdaout)
       d = maxval(thetaout(2:ntheta) - thetaout(1:ntheta-1))
       print *, "max deltheta: ", d
       nin = count(abs(thetaout(1:ntheta) - thetamax) < widthw &
                   .or. abs(thetaout(1:ntheta) + thetamax) < widthw)
       print *, "count in +ve cvdrift: ", nin
       print *, "fraction in +ve cvdrift: ", real(nin)/real(ntheta)
       print *, "deltaw: ", deltaw
    end if

  end subroutine generate_grids

  subroutine spline (nin, xin, yin, nout, xout, yout)
    use splines, only: fitp_curv1, fitp_curv2
    implicit none
    integer, intent (in) :: nin
    real, dimension (:), intent (in) :: xin, yin
    integer, intent (in) :: nout
    real, dimension (:), intent (in) :: xout
    real, dimension (:), intent (out) :: yout

    real, dimension (:), allocatable :: work, tmp
    integer :: ierr
    integer :: i

    allocate (work(nin), tmp(2*nin))

    ierr = 0
    call fitp_curv1 (nin, xin, yin, 0.0, 0.0, 3, work, tmp, 1.0, ierr)
    do i = 1, nout
       yout(i) = fitp_curv2 (xout(i), nin, xin, yin, work, 1.0)
    end do

    deallocate (work, tmp)
  end subroutine spline

  subroutine write_output
    use splines, only: fitp_curvp2
    implicit none
    integer :: unit
    integer :: i
!    real, allocatable, dimension (:) :: bmaginaux, bmagsmaux, tmp
    real, allocatable, dimension (:) :: thetagridout
    real, allocatable, dimension (:) :: gbdriftout, gradparout, grhoout
    real, allocatable, dimension (:) :: cvdriftout, gds2out, bmaggridout
    real, allocatable, dimension (:) :: gds21out, gds22out
    real, allocatable, dimension (:) :: cvdrift0out, gbdrift0out
    real, allocatable, dimension (:) :: Rplotout, Rprimeout
    real, allocatable, dimension (:) :: Zplotout, Zprimeout
    real, allocatable, dimension (:) :: aplotout, aprimeout
!    real, external :: fitp_curvp2
!    real :: th, bmin
    real :: bmin
!    integer :: ierr
    integer, dimension (1) :: minloca

    call open_output_file (unit, ".input.out")
    write (unit=unit, fmt="('#',i5)") nthetain+1
    do i = 1, nthetain+1
       write (unit=unit, fmt="(3(1x,g19.12))") thetain(i), bmagin(i), bmagsm(i)
    enddo
    call close_output_file (unit)

!    allocate (bmaginaux(nthetain+1), bmagsmaux(nthetain+1), tmp(2*nthetain))
!    call fitp_curvp1 (nthetain,thetain,bmagin,twopi,bmaginaux,tmp,1.0,ierr)
!    call fitp_curvp1 (nthetain,thetain,bmagsm,twopi,bmagsmaux,tmp,1.0,ierr)
!    call open_output_file (unit, ".input.fine")
!    do i = -nfinegrid, nfinegrid
!       th = pi*real(i)/real(nfinegrid)
!       write (unit=unit, fmt="(3(x,g19.12))") th, &
!            fitp_curvp2(th,nthetain,thetain,bmagin,twopi,bmaginaux,1.0), &
!            fitp_curvp2(th,nthetain,thetain,bmagsm,twopi,bmagsmaux,1.0)
!    end do
!    call close_output_file (unit)
!    deallocate (bmaginaux, bmagsmaux, tmp)

    call open_output_file (unit, ".bmag.out")
    write (unit=unit, fmt="('#',i5)") ntheta
    do i = 1, ntheta
       write (unit=unit, fmt="(2(1x,g19.12))") thetaout(i), bmagout(i)
    enddo
    call close_output_file (unit)

    call open_output_file (unit, ".lambda.out")
    write (unit=unit, fmt="('#',i5)") nlambda
    do i = 1, nlambda
       write (unit=unit, fmt=*) alambdaout(i)
    enddo
    call close_output_file (unit)

    ntgrid = ntheta/2 + (nperiodout-1)*ntheta
    allocate (thetagridout(-ntgrid:ntgrid))
    allocate (gbdriftout(-ntgrid:ntgrid))
    allocate (gradparout(-ntgrid:ntgrid))
    allocate (grhoout(-ntgrid:ntgrid))
    allocate (cvdriftout(-ntgrid:ntgrid))
    allocate (gds2out(-ntgrid:ntgrid))
    allocate (bmaggridout(-ntgrid:ntgrid))
    allocate (gds21out(-ntgrid:ntgrid))
    allocate (gds22out(-ntgrid:ntgrid))
    allocate (cvdrift0out(-ntgrid:ntgrid))
    allocate (gbdrift0out(-ntgrid:ntgrid))
    allocate (Rplotout(-ntgrid:ntgrid))
    allocate (Rprimeout(-ntgrid:ntgrid))
    allocate (Zplotout(-ntgrid:ntgrid))
    allocate (Zprimeout(-ntgrid:ntgrid))
    allocate (aplotout(-ntgrid:ntgrid))
    allocate (aprimeout(-ntgrid:ntgrid))

    thetagridout(-ntheta/2:ntheta/2-1) = thetaout(:ntheta)
    bmaggridout(-ntheta/2:ntheta/2-1) = bmagout(:ntheta)
    thetagridout(ntheta/2) = thetagridout(-ntheta/2) + twopi*iperiod
    bmaggridout(ntheta/2) = bmaggridout(-ntheta/2)
    do i = 1, nperiodout - 1
       thetagridout(-ntheta/2-ntheta*i:ntheta/2-ntheta*i-1) &
            = thetagridout(-ntheta/2:ntheta/2-1) - twopi*i
       thetagridout(-ntheta/2+ntheta*i+1:ntheta/2+ntheta*i) &
            = thetagridout(-ntheta/2+1:ntheta/2) + twopi*i
       bmaggridout(-ntheta/2-ntheta*i:ntheta/2-ntheta*i-1) &
            = bmaggridout(-ntheta/2:ntheta/2-1)
       bmaggridout(-ntheta/2+ntheta*i+1:ntheta/2+ntheta*i) &
            = bmaggridout(-ntheta/2+1:ntheta/2)
    end do

    call spline (2*ntgridin+1,thetagrid,gbdrift, &
         2*ntgrid+1,thetagridout,gbdriftout)
    call spline (2*ntgridin+1,thetagrid,gradpar, &
         2*ntgrid+1,thetagridout,gradparout)
    call spline (2*ntgridin+1,thetagrid,grho, &
         2*ntgrid+1,thetagridout,grhoout)
    call spline (2*ntgridin+1,thetagrid,cvdrift, &
         2*ntgrid+1,thetagridout,cvdriftout)
    call spline (2*ntgridin+1,thetagrid,gds2, &
         2*ntgrid+1,thetagridout,gds2out)
    call spline (2*ntgridin+1,thetagrid,bmag, &
         2*ntgrid+1,thetagridout,bmaggridout)
    call spline (2*ntgridin+1,thetagrid,gds21, &
         2*ntgrid+1,thetagridout,gds21out)
    call spline (2*ntgridin+1,thetagrid,gds22, &
         2*ntgrid+1,thetagridout,gds22out)
    call spline (2*ntgridin+1,thetagrid,cvdrift0, &
         2*ntgrid+1,thetagridout,cvdrift0out)
    call spline (2*ntgridin+1,thetagrid,gbdrift0, &
         2*ntgrid+1,thetagridout,gbdrift0out)
    call spline (2*ntgridin+1,thetagrid,Rplot, &
         2*ntgrid+1,thetagridout,Rplotout)
    call spline (2*ntgridin+1,thetagrid,Rprime, &
         2*ntgrid+1,thetagridout,Rprimeout)
    call spline (2*ntgridin+1,thetagrid,Zplot, &
         2*ntgrid+1,thetagridout,Zplotout)
    call spline (2*ntgridin+1,thetagrid,Zprime, &
         2*ntgrid+1,thetagridout,Zprimeout)
    call spline (2*ntgridin+1,thetagrid,aplot, &
         2*ntgrid+1,thetagridout,aplotout)
    call spline (2*ntgridin+1,thetagrid,aprime, &
         2*ntgrid+1,thetagridout,aprimeout)

    call open_output_file (unit, ".out")
    write (unit=unit, fmt=*) 'nlambda'
    write (unit=unit, fmt=*) nlambda
    write (unit=unit, fmt=*) 'lambda'
    do i = 1, nlambda
       write (unit=unit, fmt=*) alambdaout(i)
    enddo

    write (unit=unit, fmt=*) 'ntgrid nperiod ntheta drhodpsi rmaj shat kxfac q'
    write (unit=unit, fmt="(3i4,5(1x,g19.10))") ntgrid, nperiodout, ntheta, drhodpsi, rmaj, shat, kxfac, qval

    write (unit=unit, fmt=*) 'gbdrift gradpar grho tgrid'
    do i = -ntgrid,ntgrid
       write (unit=unit, fmt="(8(1x,g19.10))") &
            gbdriftout(i),gradparout(i),grhoout(i),thetagridout(i)
    end do

    write (unit=unit, fmt=*) 'cvdrift gds2 bmag tgrid'
    do i = -ntgrid,ntgrid
       write (unit=unit, fmt="(8(1x,g19.10))") &
            cvdriftout(i),gds2out(i),bmaggridout(i),thetagridout(i)
    end do

    write (unit=unit, fmt=*) 'gds21 gds22 tgrid'
    do i = -ntgrid,ntgrid
       write (unit=unit, fmt="(8(1x,g19.10))") &
            gds21out(i),gds22out(i),thetagridout(i)
    end do

    write (unit=unit, fmt=*) 'cvdrift0 gbdrift0 tgrid'
    do i = -ntgrid,ntgrid
       write (unit=unit, fmt="(8(1x,g19.10))") &
            cvdrift0out(i),gbdrift0out(i),thetagridout(i)
    end do

    write (unit=unit, fmt=*) 'Rplot Rprime tgrid'
    do i = -ntgrid,ntgrid
       write (unit=unit, fmt="(8(1x,g19.10))") &
            Rplotout(i),Rprimeout(i),thetagridout(i)
    end do

    write (unit=unit, fmt=*) 'Zplot Zprime tgrid'
    do i = -ntgrid,ntgrid
       write (unit=unit, fmt="(8(1x,g19.10))") &
            Zplotout(i),Zprimeout(i),thetagridout(i)
    end do

    write (unit=unit, fmt=*) 'aplot aprime tgrid'
    do i = -ntgrid,ntgrid
       write (unit=unit, fmt="(8(1x,g19.10))") &
            aplotout(i),aprimeout(i),thetagridout(i)
    end do

    call close_output_file (unit)

    if (screenout) then
       write (*, *) 'cvdrift       gds2          bmag          theta'
       do i = -ntheta/2,ntheta/2
          write (unit=*, fmt="(4(1x,g13.5))") &
               cvdriftout(i),gds2out(i),bmaggridout(i),thetagridout(i)
       end do
    end if

    minloca = minloc(bmagout(:ntheta))
    bmin = bmagout(minloca(1))
    print *, "theta,bmag minimum: ", thetaout(minloca(1)), bmin

    minloca = minloc(bmagin)
    print *, "theta,bmag input minimum: ", &
         thetain(minloca(1)), bmagin(minloca(1))

    if (bmin < bmagin(minloca(1))) print *, "WARNING: interpolated new minimum"

  end subroutine write_output

end program rungridgen
