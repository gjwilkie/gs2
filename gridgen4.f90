subroutine gridgen4read (filename,ntheta,nperiod,ntgrid,nlambda, &
     theta,alambda, &
     gbdrift,gradpar,cvdrift,gds2,bmag,gds21,gds22,cvdrift0,gbdrift0)
  implicit none
  character(*), intent (in) :: filename
  integer, intent (in out) :: ntheta, nperiod, ntgrid, nlambda
  real, dimension (-ntgrid:ntgrid) :: theta
  real, dimension (nlambda) :: alambda
  real, dimension (-ntgrid:ntgrid) :: gbdrift, gradpar
  real, dimension (-ntgrid:ntgrid) :: cvdrift, gds2, bmag
  real, dimension (-ntgrid:ntgrid) :: gds21, gds22
  real, dimension (-ntgrid:ntgrid) :: cvdrift0, gbdrift0
  !
  integer :: nthetain, nperiodin, ntgridin, nlambdain
  integer :: i
  integer :: unit
  character(20) :: read, write, readwrite
  character(200) :: line

  unit = 10
  do
     inquire (unit=unit, read=read, write=write, readwrite=readwrite)
     if (read == "UNKNOWN" .and. write == "UNKNOWN" &
          .and. readwrite == "UNKNOWN") &
     then
        exit
     end if
     unit = unit + 1
  end do
  open (unit=unit, file=filename, status="old")

  read (unit=unit, fmt="(a)") line
  read (unit=unit, fmt=*) nlambdain
  if (nlambdain > nlambda) then
     print *, "grid.out:nlambda > nlambda: ", nlambdain, nlambda
     stop
  end if
  nlambda = nlambdain
  read (unit=unit, fmt="(a)") line
  do i = 1, nlambda
     read (unit=unit, fmt=*) alambda(i)
  end do

  read (unit=unit, fmt="(a)") line
  read (unit=unit, fmt=*) ntgridin, nperiodin, nthetain
  if (ntgridin > ntgrid) then
     print *, "grid.out:ntgrid > ntgrid: ", ntgridin, ntgrid
     stop
  end if
  ntgrid = ntgridin
  nperiod = nperiodin
  ntheta = nthetain

  read (unit=unit, fmt="(a)") line
  do i = -ntgrid, ntgrid
     read (unit=unit, fmt=*) gbdrift(i), gradpar(i)
  end do

  read (unit=unit, fmt="(a)") line
  do i = -ntgrid, ntgrid
     read (unit=unit, fmt=*) cvdrift(i), gds2(i), bmag(i), theta(i)
  end do

  read (unit=unit, fmt="(a)") line
  do i = -ntgrid, ntgrid
     read (unit=unit, fmt=*) gds21(i), gds22(i)
  end do

  read (unit=unit, fmt="(a)") line
  do i = -ntgrid, ntgrid
     read (unit=unit, fmt=*) cvdrift0(i), gbdrift0(i)
  end do

  close (unit=unit)
end subroutine gridgen4read

subroutine gridgen4 (nbmag,thetain,bmagin, npadd, &
     alknob,epsknob,extrknob,thetamax,deltaw,widthw, &
     ntheta,nlambda,thetagrid,bmaggrid,alambdagrid)
  implicit none
  integer, intent (in) :: nbmag
  real, dimension (nbmag), intent (in) :: thetain, bmagin
  integer, intent (in) :: npadd
  real, intent (in) :: alknob, epsknob, extrknob
  real, intent (in) :: thetamax, deltaw, widthw
  integer, intent (in out) :: ntheta, nlambda
  real, dimension (ntheta), intent (out) :: thetagrid, bmaggrid
  real, dimension (nlambda), intent (out) :: alambdagrid

  call gridgen4_1 (nbmag,thetain,bmagin, npadd, &
     alknob,epsknob,extrknob,thetamax,deltaw,widthw,1.0, &
     ntheta,nlambda,thetagrid,bmaggrid,alambdagrid)
end subroutine gridgen4

subroutine gridgen4_1 (nbmag,thetain,bmagin, npadd, &
     alknob,epsknob,extrknob,thetamax,deltaw,widthw,tension, &
     ntheta,nlambda,thetagrid,bmaggrid,alambdagrid)
  implicit none
  integer, intent (in) :: nbmag
  real, dimension (nbmag), intent (in) :: thetain, bmagin
  integer, intent (in) :: npadd
  real, intent (in) :: alknob, epsknob, extrknob
  real, intent (in) :: thetamax, deltaw, widthw, tension
  integer, intent (in out) :: ntheta, nlambda
  real, dimension (ntheta), intent (out) :: thetagrid, bmaggrid
  real, dimension (nlambda), intent (out) :: alambdagrid

  call gridgen4_2 (nbmag,thetain,bmagin, npadd, &
     alknob,epsknob,extrknob,thetamax,deltaw,widthw,tension, &
     ntheta,nlambda,thetagrid,bmaggrid,alambdagrid)
  alambdagrid(1:nlambda) = 1.0/alambdagrid(1:nlambda)
end subroutine gridgen4_1

subroutine gridgen4_2 (nbmag,thetain,bmagin, npadd, &
     alknob,epsknob,extrknob,thetamax,deltaw,widthw,tension, &
     ntheta,nbset,thetagrid,bmaggrid,bset)
  implicit none
  integer, intent (in) :: nbmag
  real, dimension (nbmag), intent (in) :: thetain, bmagin
  integer, intent (in) :: npadd
  real, intent (in) :: alknob, epsknob, extrknob
  real, intent (in) :: thetamax, deltaw, widthw, tension
  integer, intent (in out) :: ntheta, nbset
  real, dimension (ntheta), intent (out) :: thetagrid, bmaggrid
  real, dimension (nbset), intent (out) :: bset
  !
  real, parameter :: pi = 3.1415926535897931, twopi = 6.2831853071795862
  integer, parameter :: maxiter = 100
  logical, parameter :: debug_output = .true.
  !
  real, dimension (nbmag) :: bmagspl
  !
  integer :: nstart
  real, allocatable, dimension (:) :: thetastart, bmagstart
  logical, allocatable, dimension (:) :: essentialstart
  real, allocatable, dimension (:) :: bmin, bmax, bextr, thetaextr
  real, allocatable, dimension (:) :: bprime
  !
  integer :: nset, nsetset
  real, allocatable, dimension (:) :: thetasetset, bmagset
  integer, allocatable, dimension (:) :: ibmagsetset, icollsetset
  logical, allocatable, dimension (:) :: essentialset
  !
  integer, allocatable, dimension (:) :: ithetasort, ibmagsort
  !
  real, allocatable, dimension (:) :: thetares, alambdares
  !
  integer :: nthetaout, nlambdaout
  !
  integer :: debug_unit

  call gg4init
  call gg4debug (nbmag, thetain, bmagin, "input grid")

! 1 Set up starting grid.
  call gg4start
  call gg4debug (nstart, thetastart, thetastart, "starting grid")

! 2 Collect all bounce points associated with starting grid.
  call gg4collect
  call gg4debug (nset, bmagset, bmagset, "bmagset")
  call gg4debugi (nsetset,thetasetset,ibmagsetset,icollsetset,"thetasetset")

! 3 Sort collected grids.
  call gg4sort

! 4 Calculate spacing and resolution metrics.
  call gg4metrics

! 5 Remove sets until the number of points is small enough.
  call gg4remove

! 6 Build output grids.
  call gg4results
  call gg4debug (nthetaout, thetagrid, bmaggrid, "output grid")
  call gg4debug (nlambdaout, bset, bset, "output 1/lambda")

  ntheta = 2*(nthetaout/2)
  nbset = nlambdaout

! 7 Clean up.
  call gg4finish

contains
  
  subroutine gg4init
    implicit none
    character(20) :: read, write, readwrite
    real, dimension (2*nbmag) :: tmp
    integer :: ierr, i

    ierr = 0
    call fitp_curvp1 (nbmag-1,thetain,bmagin,twopi,bmagspl,tmp,tension,ierr)
    if (ierr /= 0) then
       print *, "CURVP1: IERR=", ierr
       select case (ierr)
       case (1)
          print *, "N is less than 2"
       case (2)
          print *, "P is less than or equal to X(N)-X(1)"
       case (3)
          print *, "X values are not strictly increasing"
       end select
       stop
    end if
    !
    nstart = 0
    !
    nset = 0
    nsetset = 0
    !
    if (debug_output) then
       debug_unit = 10
       do
          inquire (unit=debug_unit,read=read,write=write,readwrite=readwrite)
          if (read == "UNKNOWN" .and. write == "UNKNOWN" &
               .and. readwrite == "UNKNOWN") &
          then
             exit
          end if
          debug_unit = debug_unit + 1
       end do
       open (unit=debug_unit, file="gridgen.200", status="unknown")
       write (unit=debug_unit, fmt=*) "nbmag=", nbmag
       write (unit=debug_unit, fmt=*) "thetain,bmagin="
       do i=1,nbmag
          write (unit=debug_unit, fmt=*) thetain(i), bmagin(i)
       end do
       write (unit=debug_unit, fmt=*) "alknob=", alknob
       write (unit=debug_unit, fmt=*) "epsknob=", epsknob
       write (unit=debug_unit, fmt=*) "extrknob=", extrknob
       write (unit=debug_unit, fmt=*) "thetamax=", thetamax
       write (unit=debug_unit, fmt=*) "deltaw,widthw=", deltaw, widthw
       write (unit=debug_unit, fmt=*) "tension=", tension
       write (unit=debug_unit, fmt=*) "ntheta,nlambda=", ntheta, nbset
    end if
  end subroutine gg4init

  subroutine gg4finish
    implicit none

    if (debug_output) then
       close (unit=debug_unit)
    end if
    !
    deallocate (thetastart)
    deallocate (bmagstart)
    deallocate (bmin,bmax,bextr,thetaextr,bprime)
    deallocate (essentialstart)
    !
    deallocate (thetasetset, bmagset)
    deallocate (ibmagsetset, icollsetset)
    deallocate (essentialset)
    !
    deallocate (ithetasort, ibmagsort)
    !
    deallocate (thetares, alambdares)

    if (mod(nthetaout,2) /= 1) then
       print *, "gridgen4_1:gg4results:nthetaout=",nthetaout
       print *, "nthetaout is not odd, so there is a problem with gridgen"
       stop
    end if
  end subroutine gg4finish

  subroutine old_gg4finish
    implicit none

    deallocate (thetastart, essentialstart)
    !
    deallocate (thetasetset, bmagset)
    deallocate (ibmagsetset, icollsetset)
    deallocate (essentialset)
    !
    deallocate (ithetasort, ibmagsort)
    !
    deallocate (thetares, alambdares)
    !
    if (debug_output) then
       close (unit=debug_unit)
    end if
  end subroutine old_gg4finish

  subroutine gg4start
    implicit none
    integer :: i, j
    real :: thetal, thetar, theta0
    real :: bprime0
    integer :: nextr

    allocate (bmin(nbmag-1),bmax(nbmag-1),bextr(nbmag-1),thetaextr(nbmag-1))
    allocate (bprime(nbmag))
    bmin = min(bmagin(1:nbmag-1),bmagin(2:nbmag))
    bmax = max(bmagin(1:nbmag-1),bmagin(2:nbmag))
    bextr = 0.0
    thetaextr = 0.0

    do i = 1, nbmag
       bprime(i) = bmagp(thetain(i))
    end do

! Find all extrema.
    nextr = 0
    do i = 2, nbmag-2
       if (bprime(i) == 0.0) then
          bextr(i) = bmagin(i)
          thetaextr(i) = thetain(i)
          nextr = nextr + 1
       else if (bprime(i) < 0.0 .and. bprime(i+1) > 0.0) then
          thetal = thetain(i)
          thetar = thetain(i+1)
          do j = 1, maxiter
             theta0 = 0.5*(thetal+thetar)
             bprime0 = bmagp(theta0)
             if (bprime0 < 0.0) then
                thetal = theta0
             else if (bprime0 > 0.0) then
                thetar = theta0
             else
                exit
             end if
          end do
          if (theta0-thetain(i) < epsknob) then
             theta0 = thetain(i)
             bprime(i) = 0.0
          else if (thetain(i+1)-theta0 < epsknob) then
             theta0 = thetain(i+1)
             bprime(i+1) = 0.0
          end if
          thetaextr(i) = theta0
          bextr(i) = bmagint(theta0)
          nextr = nextr + 1
       else if (bprime(i) > 0.0 .and. bprime(i+1) < 0.0) then
          thetal = thetain(i)
          thetar = thetain(i+1)
          do j = 1, maxiter
             theta0 = 0.5*(thetal+thetar)
             bprime0 = bmagp(theta0)
             if (bprime0 > 0.0) then
                thetal = theta0
             else if (bprime0 < 0.0) then
                thetar = theta0
             else
                exit
             end if
          end do
          if (theta0-thetain(i) < epsknob) then
             theta0 = thetain(i)
             bprime(i) = 0.0
          else if (thetain(i+1)-theta0 < epsknob) then
             theta0 = thetain(i+1)
             bprime(i+1) = 0.0
          end if
          thetaextr(i) = theta0
          bextr(i) = bmagint(theta0)
          nextr = nextr + 1
       end if
    end do

! Collect -pi, all local extrema, original grid points, points between
! original grid points
    nstart = 1 + nextr + (nbmag-1)*(1+npadd)

    allocate (thetastart(nstart), essentialstart(nstart), bmagstart(nstart))
    essentialstart(:1+nextr) = .true.
    essentialstart(2+nextr:) = .false.
    thetastart(1) = -pi
    thetastart(2:nextr+1) = pack(thetaextr,bextr /= 0.0)

    nstart = nextr + 1
    do i = 1, nbmag-1
       do j = 0, npadd
          thetastart(nstart+1) &
               = thetain(i) + real(j)/real(npadd+1)*(thetain(i+1)-thetain(i))
          if (all(abs(thetastart(nstart+1)-thetastart(1:nstart)) > epsknob)) &
          then
             nstart = nstart + 1
          end if
       end do
    end do

    do i = 1, nstart
       bmagstart(i) = bmagint(thetastart(i))
    end do
  end subroutine gg4start

  subroutine old_gg4start
    implicit none
    integer :: i, j
    real :: bprimel, bprimer, bprime0
    real :: thetal, thetar, theta0
    real, parameter :: tol=1e-9
    
    allocate (thetastart(nbmag*(2+npadd)))
    allocate (essentialstart(nbmag*(2+npadd)))

! 1 Collect essential points: -pi + all local extrema.
! 1.1 Collect -pi.
    call add_start (-pi, .true.)
! 1.2 Collect extrema.
    bprimer = bmagp(thetain(2))
    do i = 2, nbmag-2
       bprimel = bprimer
       bprimer = bmagp(thetain(i+1))
       if (bprimel == 0.0) then
          call add_start (thetain(i), .true.)
       else if (bprimel < 0.0 .and. bprimer > 0.0) then
! 1.2.1 Find local minimum.
          thetal = thetain(i)
          thetar = thetain(i+1)
          do j = 1, maxiter
             theta0 = 0.5*(thetal+thetar)
             bprime0 = bmagp(theta0)
             if (bprime0 < 0.0) then
                thetal = theta0
             else if (bprime0 > 0.0) then
                thetar = theta0
             else
                exit
             end if
          end do
          if (theta0-thetain(i) < tol) theta0 = thetain(i)
          if (thetain(i+1)-theta0 < tol) theta0 = thetain(i+1)
          call add_start (theta0, .true.)
       else if (bprimel > 0.0 .and. bprimer < 0.0) then
! 1.2.2 Find local maximum.
          thetal = thetain(i)
          thetar = thetain(i+1)
          do j = 1, maxiter
             theta0 = 0.5*(thetal+thetar)
             bprime0 = bmagp(theta0)
             if (bprime0 > 0.0) then
                thetal = theta0
             else if (bprime0 < 0.0) then
                thetar = theta0
             else
                exit
             end if
          end do
          if (theta0-thetain(i) < tol) then
             theta0 = thetain(i)
             bprime(i) = 0.0
          else if (thetain(i+1)-theta0 < tol) then
             theta0 = thetain(i+1)
             bprime(i+1) = 0.0
          end if
          call add_start (theta0, .true.)
       end if
    end do

! 2 Collect original grid, except for extrema.
    do i = 2, nbmag-1
       if (all(abs(thetain(i)-thetastart(1:nstart)) > tol)) then
          call add_start(thetain(i), .false.)
       end if
    end do

! 3 Collect points between original grid points.
    do i = 1, nbmag-1
       do j = 1, npadd
          theta0 = thetain(i) + real(j)/real(npadd+1)*(thetain(i+1)-thetain(i))
          if (all(abs(theta0-thetastart(1:nstart)) > tol)) then
             call add_start(theta0, .false.)
          end if
       end do
    end do
  end subroutine old_gg4start

  subroutine add_start (theta, essential)
    implicit none
    real, intent (in) :: theta
    logical, intent (in) :: essential
    !
    real, allocatable, dimension (:) :: thetatmp
    logical, allocatable, dimension (:) :: essentialtmp

    if (nstart >= size(thetastart)) then
       allocate (thetatmp(nstart))
       allocate (essentialtmp(nstart))
       thetatmp = thetastart(1:nstart)
       essentialtmp = essentialstart(1:nstart)
       deallocate (thetastart,essentialstart)
       allocate (thetastart(nstart+nbmag))
       allocate (essentialstart(nstart+nbmag))
       thetastart(1:nstart) = thetatmp
       essentialstart(1:nstart) = essentialtmp
       deallocate (thetatmp,essentialtmp)
    end if
    nstart = nstart + 1
    thetastart(nstart) = theta
    essentialstart(nstart) = essential
  end subroutine add_start

  subroutine gg4collect
    implicit none
    integer :: i, iset
    integer :: nsetsetmax
    real :: thetai, bmagi

! Estimate upper bound on number of points.
    nsetsetmax = count( &
     (spread(bmagstart,2,nbmag-1) >= spread(bmin,1,nstart) &
      .and. spread(bmagstart,2,nbmag-1) <= spread(bmax,1,nstart))) &
             + 2*count( &
     (spread(bextr,1,nstart) /= 0 &
      .and. &
       ((spread(bmagstart,2,nbmag-1) >= spread(bmax,1,nstart) &
         .and. spread(bmagstart,2,nbmag-1) <= spread(bextr,1,nstart)) &
        .or. &
        (spread(bmagstart,2,nbmag-1) <= spread(bmin,1,nstart) &
         .and. spread(bmagstart,2,nbmag-1) >= spread(bextr,1,nstart)))))

! Allocate
    allocate (thetasetset(nsetsetmax))
    allocate (ibmagsetset(nsetsetmax))
    allocate (icollsetset(nsetsetmax))

    nsetset = 0
    starting_points: do iset = 1, nstart
       thetai = thetastart(iset)
       bmagi = bmagstart(iset)

! 1 For extrema, check previous sets to attach to.
       if (essentialstart(iset)) then
          do i = 1, iset-1
             if (abs(bmagstart(i)-bmagi) < epsknob) then
! 1.1.1 If the extremum does belong in this set, eliminate points
!       near the extremum from this set.
                where (ibmagsetset(1:nsetset) == i &
                       .and. abs(thetasetset(1:nsetset)-thetai) < epsknob)
                   ibmagsetset(1:nsetset) = 0
                end where
! 1.1.2 Attach the extremum to this set.
                call add_setset (thetai, i, 112)
                bmagstart(iset) = 0.0
                cycle starting_points
             end if
          end do
       end if

! 2 Start a new set.
       call add_setset (thetai, iset, 2)

! 3 Check each original grid interval for matching bmag.
       grid_interval: do i = 1, nbmag-1
! 3.0.1 Stoopid problems near -pi.
          if (iset == 1 .and. i == 1) cycle grid_interval
          if (bmagin(i) > bmagi .and. thetain(i) /= thetai) then
! 3.1 Consider when the left grid point is greater than the target bmag.
!    Then, there are three cases in which there are matching points.
!    (1) The right grid point is equal to the target bmag, and the slope
!        at the right point is positive.
!    (2) The right gridpoint is less than the target bmag.
!    (3) The right grid point is greater than the target bmag,
!        and the interval is concave down, and the minimum in
!        this interval, which is guaranteed to exist, is less than
!        the target bmag.
             if ((bmagin(i+1) == bmagi .or. thetain(i+1) == thetai) &
                  .and. bprime(i+1) > 0.0) &
             then
! 3.1.1 Consider when the right grid point is equal to the target bmag.
                call add_setset_root (thetain(i+1),thetain(i),bmagi,iset,311)
                cycle grid_interval
             end if
             if (bmagin(i+1) < bmagi) then
! 3.1.2 Consider when the right grid point is less than the target bmag.
! 3.1.2.1 If this interval bounds the starting theta point, that is what
!         the target point is, and it is already collected.
                if (thetai >= thetain(i) .and. thetai <= thetain(i+1)) then
                   cycle grid_interval
                end if
! 3.1.2.2 Otherwise, find and collect the target point.
                call add_setset_root (thetain(i+1),thetain(i),bmagi,iset,3122)
                cycle grid_interval
             end if
! 3.1.3 Check if the grid interval is concave down.
! 3.1.3.1 If not, skip to the next interval.
             if (bprime(i) >= 0.0 .or. bprime(i+1) <= 0.0) cycle grid_interval
! 3.1.3.2 Consider the case where the starting theta point is within
!         this interval.
             if (thetai > thetain(i) .and. thetai < thetain(i+1)) then
! 3.1.3.2.1 If the starting point is an extremum, skip to next interval.
                if (essentialstart(iset)) cycle grid_interval
! 3.1.3.2.2 Otherwise, the other target point is right of the starting
!           point is the slope is negative, and left if positive.
                if (bprime(i) < 0.0) then
                   call add_setset_root (thetai,thetain(i+1),bmagi,iset,313221)
                else
                   call add_setset_root (thetai,thetain(i),bmagi,iset,313222)
                end if
                cycle grid_interval
             end if
! 3.1.3.3 If the minimum within this interval is less than the target bmag,
!         then there are two target points in this interval, one on each side
!         of the minimum.
! 3.1.3.3.1 If this interval is on the edge, there will not be any minimum.
             if (i == 1 .or. i == nbmag-1) cycle grid_interval
! 3.1.3.3.2 Find the minimum in this interval.
             if (bextr(i) == 0.0) then
                print *, "gridgen4.f90:gg4collect:3.1.3.3.2:", &
                     " missing extremum"
                print *, "iset,i:",iset,i
                print *, "bmagi:",bmagi
                print *, "bprimel,bprimer:", bprime(i),bprime(i+1)
                print *, "thetain(i):", thetain(i)
                print *, "thetain(i+1):", thetain(i+1)
                print *, "thetai:", thetai
                stop
             end if
! 3.1.3.3.2.1 If the minimum is greater than the target bmag, skip to
!             the next interval.
             if (bextr(i) > bmagi) cycle grid_interval
! 3.1.3.3.2.2 Collect the point left of the minimum.
             call add_setset_root (thetaextr(i),thetain(i),&
                  bmagi,iset,313322)
! 3.1.3.3.2.3 Collect the point right of the minimum.
             call add_setset_root (thetaextr(i),thetain(i+1),&
                  bmagi,iset,313323)
             cycle grid_interval
          else if (bmagin(i) < bmagi .and. thetain(i) /= thetai) then
! 3.2 Consider then the left grid point is less than the target bmag.
!     Then, there are three cases in which there are matching points.
!     (1) The right grid point is equal to the target bmag, and the
!         slope at the right point is negative.
!     (2) The right grid point is greater than the target bmag.
!     (3) The right grid point is less than the target bmag,
!         and the interval is concave up, and the maximum in
!         this interval, which is guaranteed to exist, is greater
!         than the target bmag.
! 3.2.1 Consider when the right grid point is equal to the target bmag.
             if ((bmagin(i+1) == bmagi .or. thetain(i+1) == thetai) &
                  .and. bprime(i+1) < 0.0) &
             then
                call add_setset_root (thetain(i),thetain(i+1),bmagi,iset,321)
                cycle grid_interval
             end if
             if (bmagin(i+1) > bmagi) then
! 3.2.2 Consider when the right grid point is greater than the target bmag.
! 3.2.2.1 If this interval bounds the starting theta point, that is what
!         the target point is, and it is already collected.
                if (thetai >= thetain(i) .and. thetai <= thetain(i+1)) then
                   cycle grid_interval
                end if
! 3.2.2.2 Otherwise, find and collect the target point.
                call add_setset_root (thetain(i),thetain(i+1),bmagi,iset,3222)
                cycle grid_interval
             end if
! 3.2.3 Check if the grid interval is concave up.
! 3.2.3.1 If not, skip to the next interval.
             if (bprime(i) <= 0.0 .or. bprime(i+1) >= 0.0) cycle grid_interval
! 3.2.3.2 Consider the case where the starting theta point is within
!         this interval.
             if (thetai > thetain(i) .and. thetai < thetain(i+1)) then
! 3.2.3.2.1 If the starting point is an extremum, skip to next interval.
                if (essentialstart(iset)) cycle grid_interval
! 3.2.3.2.2 Otherwise, the other target point is right of the starting
!           point if the slope is positive, and left if negative.
                if (bprime(i) > 0.0) then
                   call add_setset_root (thetain(i+1),thetai,bmagi,iset,323221)
                else
                   call add_setset_root (thetain(i),thetai,bmagi,iset,323222)
                end if
                cycle grid_interval
             end if
! 3.2.3.3 If the maximum within this interval is greater than the target bmag,
!         then there are two target points in this interval, one on each side
!         of the maximum.
! 3.2.3.3.1 If this interval is on the edge, there will not be any maximum.
             if (i == 1 .or. i == nbmag-1) cycle grid_interval
! 3.2.3.3.2 Find the maximum in this interval.
             if (bextr(i) == 0.0) then
                print *, "gridgen4.f90:gg4collect:3.2.3.3.2:", &
                     " missing extremum"
                print *, "iset,i:",iset,i
                print *, "bmagi:",bmagi
                print *, "bprimel,bprimer:", bprime(i),bprime(i+1)
                print *, "thetain(i):", thetain(i)
                print *, "thetain(i+1):", thetain(i+1)
                print *, "thetai:", thetai
                stop
             end if
! 3.2.3.3.2.1 If the maximum is less than the target bmag, skip to
!             the next interval.
             if (bextr(i) <= bmagi) cycle grid_interval
! 3.2.3.3.2.2 Collect the point left of the maximum.
             call add_setset_root (thetain(i),thetaextr(i),bmagi,iset,323322)
! 3.2.3.3.2.3 Collect the point right of the maximum.
             call add_setset_root (thetain(i+1),thetaextr(i),bmagi,iset,323323)
             cycle grid_interval
          else if (bmagin(i) == bmagi .or. thetain(i) == thetai) then
! 3.3 Consider when then left grid point is equal to the target bmag.
! 3.3.1 Add the point if it is not the starting grid point.
             if (thetai /= thetain(i)) then
                call add_setset (thetain(i), iset, 331)
             end if
! 3.3.2 Check if there is another matching target bmag in the interval.
             if (bprime(i) > 0.0 .and. bmagin(i+1) < bmagi &
                  .and. i /= 1 .and. i /= nbmag-1) &
             then
                call add_setset_root (thetain(i+1),thetain(i), &
                     bmagi,iset,3321)
                cycle grid_interval
             else if (bprime(i) < 0.0 .and. bmagin(i+1) > bmagi &
                  .and. i /= 1 .and. i /= nbmag-1) &
             then
                call add_setset_root (thetain(i),thetain(i+1), &
                     bmagi,iset,3322)
                cycle grid_interval
             end if
          end if

       end do grid_interval
    end do starting_points

! 4 Compatibility with the old gg4collect.
    allocate (bmagset(nstart))
    allocate (essentialset(nstart))

    nset = nstart
    bmagset = bmagstart(1:nstart)
    essentialset = essentialstart(1:nstart)

  end subroutine gg4collect

  subroutine old_gg4collect
    implicit none
    integer :: iset, i, ii, imatch
    real :: thetai, bmagi
    real :: bprimer, bprimel, bprime0

    allocate (bmagset(nbmag*(2+npadd)))
    allocate (essentialset(nbmag*(2+npadd)))
    allocate (thetasetset(nbmag*(2+npadd)*2))
    allocate (ibmagsetset(nbmag*(2+npadd)*2))
    allocate (icollsetset(nbmag*(2+npadd)*2))

    starting_points: do iset = 1, nstart
       thetai = thetastart(iset)
       bmagi = bmagint(thetai)

! 1 For extrema, check previous sets to attach to.
       if (essentialstart(iset)) then
          do i = 1, nset
! 1.1 Check if the extremum should be attached to each previous set
             if (abs(bmagi-bmagset(i)) < epsknob) then
! 1.1.1 If the extremum does belong in this set, eliminate points
!       near the extremum from this set.
                where (ibmagsetset(1:nsetset) == i &
                       .and. abs(thetasetset(1:nsetset)-thetai) < epsknob)
                   ibmagsetset(1:nsetset) = 0
                end where
! 1.1.2 Attach the extremum to this set.
                call add_setset (thetai, i, 112)
                cycle starting_points
             end if
          end do
       end if

! 2 Start a new set.
       call add_set (bmagi, essentialstart(iset))
       call add_setset (thetai, nset, 2)

! 3 Check each original grid interval for matching bmag.
       grid_interval: do i = 1, nbmag-1
! 3.0.1 Stoopid problems near -pi.
          if (nset == 1 .and. i == 1) cycle grid_interval
          if (bmagin(i) > bmagi .and. thetain(i) /= thetai) then
! 3.1 Consider when the left grid point is greater than the target bmag.
!    Then, there are three cases in which there are matching points.
!    (1) The right grid point is equal to the target bmag, and the slope
!        at the right point is positive.
!    (2) The right gridpoint is less than the target bmag.
!    (3) The right grid point is greater than the target bmag,
!        and the interval is concave down, and the minimum in
!        this interval, which is guaranteed to exist, is less than
!        the target bmag.
             bprimer = bmagp(thetain(i+1))
             if ((bmagin(i+1) == bmagi .or. thetain(i+1) == thetai) &
                  .and. bprimer > 0.0) &
             then
! 3.1.1 Consider when the right grid point is equal to the target bmag.
                call add_setset_root (thetain(i+1),thetain(i),bmagi,nset,311)
                cycle grid_interval
             end if
             if (bmagin(i+1) < bmagi) then
! 3.1.2 Consider when the right grid point is less than the target bmag.
! 3.1.2.1 If this interval bounds the starting theta point, that is what
!         the target point is, and it is already collected.
                if (thetai >= thetain(i) .and. thetai <= thetain(i+1)) then
                   cycle grid_interval
                end if
! 3.1.2.2 Otherwise, find and collect the target point.
                call add_setset_root (thetain(i+1),thetain(i),bmagi,nset,3122)
                cycle grid_interval
             end if
! 3.1.3 Check if the grid interval is concave down.
             bprimel = bmagp (thetain(i))
             bprimer = bmagp (thetain(i+1))
! 3.1.3.1 If not, skip to the next interval.
             if (bprimel >= 0.0 .or. bprimer <= 0.0) cycle grid_interval
! 3.1.3.2 Consider the case where the starting theta point is within
!         this interval.
             if (thetai > thetain(i) .and. thetai < thetain(i+1)) then
! 3.1.3.2.1 If the starting point is an extremum, skip to next interval.
                if (essentialstart(iset)) cycle grid_interval
! 3.1.3.2.2 Otherwise, the other target point is right of the starting
!           point is the slope is negative, and left if positive.
                bprime0 = bmagp(thetai)
                if (bprime0 < 0) then
                   call add_setset_root (thetai,thetain(i+1),bmagi,nset,313221)
                else
                   call add_setset_root (thetai,thetain(i),bmagi,nset,313222)
                end if
                cycle grid_interval
             end if
! 3.1.3.3 If the minimum within this interval is less than the target bmag,
!         then there are two target points in this interval, one on each side
!         of the minimum.
! 3.1.3.3.1 If this interval is on the edge, there will not be any minimum.
             if (i == 1 .or. i == nbmag-1) cycle grid_interval
! 3.1.3.3.2 Find the minimum in this interval.
             imatch = 0
             do ii = 1, nstart
                ! All essential points are at the beginning.
                if (.not. essentialstart(ii)) exit
                if (thetain(i) < thetastart(ii) &
                     .and. thetastart(ii) < thetain(i+1)) &
                then
                   if (imatch /= 0) then
                      print *, "gridgen4.f90:gg4collect:3.1.3.3.2:", &
                           " multiple extrema in interval"
                      print *, "iset,i,ii,imatch:",iset,i,ii,imatch
                      print *, "bmagi:",bmagi
                      print *, "bprimel,bprimer:", bprimel,bprimer
                      print *, "thetain(i):", thetain(i)
                      print *, "thetain(i+1):", thetain(i+1)
                      print *, "thetai:", thetai
                      stop
                   end if
                   imatch = ii
                end if
             end do
             if (imatch == 0) then
                print *, "gridgen4.f90:gg4collect:3.1.3.3.2:", &
                     " missing extremum"
                print *, "iset,i,ii,imatch:",iset,i,ii,imatch
                print *, "bmagi:",bmagi
                print *, "bprimel,bprimer:", bprimel,bprimer
                print *, "thetain(i):", thetain(i)
                print *, "thetain(i+1):", thetain(i+1)
                print *, "thetai:", thetai
                stop
             end if
! 3.1.3.3.2.1 If the minimum is greater than the target bmag, skip to
!             the next interval.
             if (bmagint(thetastart(imatch)) > bmagi) cycle grid_interval
! 3.1.3.3.2.2 Collect the point left of the minimum.
             call add_setset_root (thetastart(imatch),thetain(i),&
                  bmagi,nset,313322)
! 3.1.3.3.2.3 Collect the point right of the minimum.
             call add_setset_root (thetastart(imatch),thetain(i+1),&
                  bmagi,nset,313323)
             cycle grid_interval
          else if (bmagin(i) < bmagi .and. thetain(i) /= thetai) then
! 3.2 Consider then the left grid point is less than the target bmag.
!     Then, there are three cases in which there are matching points.
!     (1) The right grid point is equal to the target bmag, and the
!         slope at the right point is negative.
!     (2) The right grid point is greater than the target bmag.
!     (3) The right grid point is less than the target bmag,
!         and the interval is concave up, and the maximum in
!         this interval, which is guaranteed to exist, is greater
!         than the target bmag.
             bprimer = bmagp(thetain(i+1))
! 3.2.1 Consider when the right grid point is equal to the target bmag.
             if ((bmagin(i+1) == bmagi .or. thetain(i+1) == thetai) &
                  .and. bprimer < 0.0) &
             then
                call add_setset_root (thetain(i),thetain(i+1),bmagi,nset,321)
                cycle grid_interval
             end if
             if (bmagin(i+1) > bmagi) then
! 3.2.2 Consider when the right grid point is greater than the target bmag.
! 3.2.2.1 If this interval bounds the starting theta point, that is what
!         the target point is, and it is already collected.
                if (thetai >= thetain(i) .and. thetai <= thetain(i+1)) then
                   cycle grid_interval
                end if
! 3.2.2.2 Otherwise, find and collect the target point.
                call add_setset_root (thetain(i),thetain(i+1),bmagi,nset,3222)
                cycle grid_interval
             end if
! 3.2.3 Check if the grid interval is concave up.
             bprimel = bmagp(thetain(i))
             bprimer = bmagp(thetain(i+1))
! 3.2.3.1 If not, skip to the next interval.
             if (bprimel <= 0.0 .or. bprimer >= 0.0) cycle grid_interval
! 3.2.3.2 Consider the case where the starting theta point is within
!         this interval.
             if (thetai > thetain(i) .and. thetai < thetain(i+1)) then
! 3.2.3.2.1 If the starting point is an extremum, skip to next interval.
                if (essentialstart(iset)) cycle grid_interval
! 3.2.3.2.2 Otherwise, the other target point is right of the starting
!           point if the slope is positive, and left if negative.
                bprime0 = bmagp(thetai)
                if (bprime0 > 0.0) then
                   call add_setset_root (thetain(i+1),thetai,bmagi,nset,323221)
                else
                   call add_setset_root (thetain(i),thetai,bmagi,nset,323222)
                end if
                cycle grid_interval
             end if
! 3.2.3.3 If the maximum within this interval is greater than the target bmag,
!         then there are two target points in this interval, one on each side
!         of the maximum.
! 3.2.3.3.1 If this interval is on the edge, there will not be any maximum.
             if (i == 1 .or. i == nbmag-1) cycle grid_interval
! 3.2.3.3.2 Find the maximum in this interval.
             imatch = 0
             do ii = 1, nstart
                ! All essential points are at the beginning.
                if (.not. essentialstart(ii)) exit
                if (thetain(i) < thetastart(ii) &
                     .and. thetastart(ii) < thetain(i+1)) &
                then
                   if (imatch /= 0) then
                      print *, "gridgen4.f90:gg4collect:3.2.3.3.2:", &
                           " multiple extrema in interval"
                      print *, "iset,i,ii,imatch:",iset,i,ii,imatch
                      print *, "bmagi:",bmagi
                      print *, "bprimel,bprimer:", bprimel,bprimer
                      print *, "thetain(i):", thetain(i)
                      print *, "thetain(i+1):", thetain(i+1)
                      print *, "thetai:", thetai
                      stop
                   end if
                   imatch = ii
                end if
             end do
             if (imatch == 0) then
                print *, "gridgen4.f90:gg4collect:3.2.3.3.2:", &
                     " missing extremum"
                print *, "iset,i,ii,imatch:",iset,i,ii,imatch
                print *, "bmagi:",bmagi
                print *, "bprimel,bprimer:", bprimel,bprimer
                print *, "thetain(i):", thetain(i)
                print *, "thetain(i+1):", thetain(i+1)
                print *, "thetai:", thetai
                stop
             end if
! 3.2.3.3.2.1 If the maximum is less than the target bmag, skip to
!             the next interval.
             if (bmagint(thetastart(imatch)) <= bmagi) cycle grid_interval
! 3.2.3.3.2.2 Collect the point left of the maximum.
             call add_setset_root (thetain(i),thetastart(imatch),&
                  bmagi,nset,323322)
! 3.2.3.3.2.3 Collect the point right of the maximum.
             call add_setset_root (thetain(i+1),thetastart(imatch),&
                  bmagi,nset,323323)
             cycle grid_interval
          else if (bmagin(i) == bmagi .or. thetain(i) == thetai) then
! 3.3 Consider when then left grid point is equal to the target bmag.
! 3.3.1 Add the point if it is not the starting grid point.
             if (thetai /= thetain(i)) then
                call add_setset (thetain(i), nset, 331)
! 3.3.2 Check if there is another matching target bmag in the interval.
                bprime0 = bmagp(thetain(i))
                if (bprime0 > 0.0 .and. bmagin(i+1) < bmagi &
                     .and. i /= 1 .and. i /= nbmag-1) &
                then
                   call add_setset_root (thetain(i+1),thetain(i), &
                        bmagi,nset,3321)
                   cycle grid_interval
                else if (bprime0 < 0.0 .and. bmagin(i+1) > bmagi &
                     .and. i /= 1 .and. i /= nbmag-1) &
                then
                   call add_setset_root (thetain(i),thetain(i+1), &
                        bmagi,nset,3322)
                   cycle grid_interval
                end if
             end if
          end if
       end do grid_interval
    end do starting_points
  end subroutine old_gg4collect

  subroutine add_set (bmag, essential)
    implicit none
    real, intent (in) :: bmag
    logical, intent (in) :: essential
    !
    real, allocatable, dimension (:) :: bmagtmp
    logical, allocatable, dimension (:) :: essentialtmp

    if (nset >= size(bmagset)) then
       allocate (bmagtmp(nset))
       allocate (essentialtmp(nset))
       bmagtmp = bmagset(1:nset)
       essentialtmp = essentialset(1:nset)
       deallocate (bmagset, essentialset)
       allocate (bmagset(nset+nbmag))
       allocate (essentialset(nset+nbmag))
       bmagset(1:nset) = bmagtmp
       essentialset(1:nset) = essentialtmp
       deallocate (bmagtmp, essentialtmp)
    end if
    nset = nset + 1
    bmagset(nset) = bmag
    essentialset(nset) = essential
  end subroutine add_set

  subroutine add_setset (theta, ibmag, icoll)
    implicit none
    real, intent (in) :: theta
    integer, intent (in) :: ibmag, icoll
    !
    real, allocatable, dimension (:) :: thetatmp
    integer, allocatable, dimension (:) :: ibmagtmp, icolltmp

    if (nsetset >= size(thetasetset)) then
       allocate (thetatmp(nsetset))
       allocate (ibmagtmp(nsetset), icolltmp(nsetset))
       thetatmp = thetasetset(1:nsetset)
       ibmagtmp = ibmagsetset(1:nsetset)
       icolltmp = icollsetset(1:nsetset)
       deallocate (thetasetset, ibmagsetset, icollsetset)
       allocate (thetasetset(nsetset+nbmag*(2+npadd)))
       allocate (ibmagsetset(nsetset+nbmag*(2+npadd)))
       allocate (icollsetset(nsetset+nbmag*(2+npadd)))
       thetasetset(1:nsetset) = thetatmp
       ibmagsetset(1:nsetset) = ibmagtmp
       icollsetset(1:nsetset) = icolltmp
       deallocate (thetatmp, ibmagtmp, icolltmp)
    end if
    nsetset = nsetset + 1
    thetasetset(nsetset) = theta
    ibmagsetset(nsetset) = ibmag
    icollsetset(nsetset) = icoll
    if (any(ibmagsetset(1:nsetset-1) == ibmag .and. &
            abs(thetasetset(1:nsetset-1)-theta) < epsknob)) &
    then
       ibmagsetset(nsetset) = 0
    end if
  end subroutine add_setset

  subroutine add_setset_root (thetal, thetag, bmagi, ibmag, icoll)
    implicit none
    real, intent (in) :: thetal, thetag, bmagi
    integer, intent (in) :: ibmag, icoll
    !
    real :: thl, thg, theta0, bmag0
    integer :: i

    thl = thetal
    thg = thetag
    do i = 1, maxiter
       theta0 = 0.5*(thl+thg)
       bmag0 = bmagint(theta0)
       if (bmag0 > bmagi) then
          thg = theta0
       else if (bmag0 < bmagi) then
          thl = theta0
       else
          exit
       end if
    end do
    call add_setset (theta0, ibmag, icoll)
  end subroutine add_setset_root

  subroutine gg4sort
    implicit none
    allocate (ithetasort(nsetset), ibmagsort(nset))
    call heapsort (nsetset, thetasetset(1:nsetset), ithetasort)
    call heapsort (nset, bmagset(1:nset), ibmagsort)
  end subroutine gg4sort

  subroutine heapsort (n, v, index)
    implicit none
    integer, intent (in) :: n
    real, dimension (n), intent (in) :: v
    integer, dimension (n), intent (out) :: index
    integer :: i, j, l, ir, ira

    index = (/ (i, i=1,n) /)

    l = n/2+1
    ir = n

    do
       if (l > 1) then
          l = l - 1
          ira = index(l)
       else
          ira = index(ir)
          index(ir) = index(1)
          ir = ir - 1
          if (ir == 1) then
             index(1) = ira
             return
          end if
       end if
       i = l
       j = l + l
       do while (j <= ir)
          if (j < ir) then
             if (v(index(j)) < v(index(j+1))) j = j + 1
          end if
          if (v(ira) < v(index(j))) then
             index(i) = index(j)
             i = j
             j = j + j
          else
             j = ir + 1
          end if
       end do
       index(i) = ira
    end do
  end subroutine heapsort

  subroutine gg4metrics
    implicit none
    integer :: i

    allocate (thetares(nset), alambdares(nset))

    do i = 1, nset
       if (bmagset(i) /= 0.0) then
          call get_thetares (i, thetares(i))
          call get_lambdares (i, alambdares(i))
       else
          thetares(i) = 0.0
          alambdares(i) = 0.0
       end if
    end do
  end subroutine gg4metrics

  subroutine get_thetares (iset, thetares)
    implicit none
    integer, intent (in) :: iset
    real, intent (out) :: thetares
    !
    integer :: i, ileft, iright
    real :: dthetal, dthetar
    real :: res
    integer :: npts

! 1 Return large value for essential sets with odd number of points or
!   for set containing -pi.
    if (essentialset(iset)) then
! 1.1 Set containing -pi.
       if (iset == 1) then
          thetares = 2e20
          return
       end if
! 1.2 Sets with odd number of points.
! 1.2.1 Count points.
       npts = count(ibmagsetset == iset)
       if (mod(npts,2) == 1) then
          thetares = 1e20
          return
       end if
       res = extrknob/real(npts*npts)
    else
       res = 0.0
    end if

! 2 Look through all points for points in the set.
    npts = 0
    points_in_set: do i = 1, nsetset
       if (ibmagsetset(ithetasort(i)) /= iset) cycle
! 2.1 Look for the point to the left.
       ileft = i
       do
          ileft = ileft - 1
          if (ileft < 1) cycle points_in_set
          if (ibmagsetset(ithetasort(ileft)) == 0) cycle
          if (bmagset(ibmagsetset(ithetasort(ileft))) /= 0.0) exit
       end do
! 2.2 Look for the point to the right.
       iright = i
       do
          iright = iright + 1
          if (iright > nsetset) cycle points_in_set
          if (ibmagsetset(ithetasort(iright)) == 0) cycle
          if (bmagset(ibmagsetset(ithetasort(iright))) /= 0.0) exit
       end do
! 2.3 Add contribution from this interval.
       dthetal=abs(thetasetset(ithetasort(i))-thetasetset(ithetasort(ileft)))
       dthetar=abs(thetasetset(ithetasort(i))-thetasetset(ithetasort(iright)))
       res = res + tfunc(thetasetset(ithetasort(i)))*dthetal*dthetar &
            /(dthetal+dthetar+1e-20)
       npts = npts + 1
    end do points_in_set
    if (npts > 0) then
       thetares = res/real(npts)
    else
       thetares = res
    end if
  end subroutine get_thetares

  subroutine get_lambdares (iset, alambdares)
    implicit none
    integer, intent (in) :: iset
    real, intent (out) :: alambdares
    !
    integer :: i, iplus, iminus
    real :: al, alplus, alminus, dalplus, dalminus
    real :: res
    integer :: npts

    al = 1.0/bmagset(iset)
    do i = 1, nset
! 1 Look for target lambda
       if (ibmagsort(i) == iset) then
! 2 Look for bordering lambdas.
          alplus = 0.0
          do iplus = i-1, 1, -1
             if (ibmagsort(iplus) == 0) cycle
             if (bmagset(ibmagsort(iplus)) /= 0.0) then
                alplus = 1.0/bmagset(ibmagsort(iplus))
                exit
             end if
          end do
          alminus = 0.0
          do iminus = i+1, nset
             if (ibmagsort(iminus) == 0) cycle
             if (bmagset(ibmagsort(iminus)) /= 0.0) then
                alminus = 1.0/bmagset(ibmagsort(iminus))
                exit
             end if
          end do
          exit
       end if
    end do

! 3 Add up contributions to the result.
    res = 0.0
    npts = 0
    do i = 1, nbmag
       dalplus = abs(sqrt(max(1.0-alplus*bmagin(i),0.0)) &
            -sqrt(max(1.0-al*bmagin(i),0.0)))
       dalminus = abs(sqrt(max(1.0-alminus*bmagin(i),0.0)) &
            -sqrt(max(1.0-al*bmagin(i),0.0)))
       if (dalplus+dalminus /= 0.0) then
          npts = npts + 1
          res = res + dalplus*dalminus/(dalplus+dalminus+1e-20)
       end if
    end do
    if (npts /= 0) then
       alambdares = res/real(npts)
    else
       alambdares = res
    end if
  end subroutine get_lambdares

  subroutine gg4remove
    implicit none
    integer :: idel, i
    integer :: ntheta_left, nlambda_left
    integer, dimension (nset) :: work

    nlambda_left = nset
    ntheta_left = nsetset
    do i = 1, nset
       if (.not. any(ibmagsetset(1:nsetset) == i)) then
          nlambda_left = nlambda_left - 1
       end if
    end do
! 1 Find the set with the minimum resolution metric.
    do while (ntheta_left > ntheta .or. nlambda_left > nbset)
       work(1:1) = minloc(thetares(1:nset)+alknob*alambdares(1:nset), &
            bmagset(1:nset) /= 0.0)
       idel = work(1)

! 2 Delete the set just found.
       if (idel == 0) then
          print *, "gridgen4.f:gg4remove:2: This cannot happen."
          stop
       end if
       call delete_set (idel, work, ntheta_left)
       nlambda_left = nlambda_left - 1
    end do
  end subroutine gg4remove

  subroutine delete_set (idel, work, ntheta_left)
    implicit none
    integer, intent (in) :: idel
    integer, intent (in out), dimension (nset) :: work
    integer, intent (out) :: ntheta_left
    !
    integer :: i, j

    work = 0
! 1 Mark neighboring lambda sets to be recalculated.
    do i = 1, nset
       if (ibmagsort(i) == idel) then
          do j = i-1, 1, -1
             if (bmagset(ibmagsort(j)) /= 0.0) then
                work(ibmagsort(j)) = ibmagsort(j)
                exit
             end if
          end do
          do j = i+1, nset
             if (bmagset(ibmagsort(j)) /= 0.0) then
                work(ibmagsort(j)) = ibmagsort(j)
                exit
             end if
          end do
          exit
       end if
    end do

! 2 Mark lambda sets with neighboring theta points to have their resolution
!   metrics recalculated, counting the remaining theta points.
    ntheta_left = 0
    do i = 1, nsetset
       if (ibmagsetset(ithetasort(i)) == 0) cycle
       if (bmagset(ibmagsetset(ithetasort(i))) == 0.0) cycle
       if (ibmagsetset(ithetasort(i)) /= idel) then
          ntheta_left = ntheta_left + 1
          cycle
       end if
! 2.1 Found point to be deleted.
! 2.1.1 Mark set of neighboring point to the left to be recalculated.
       do j = i-1, 1, -1
          if (ibmagsetset(ithetasort(j)) == idel) cycle
          if (ibmagsetset(ithetasort(j)) == 0) cycle
          if (bmagset(ibmagsetset(ithetasort(j))) == 0.0) cycle
          work(ibmagsetset(ithetasort(j))) = ibmagsetset(ithetasort(j))
          exit
       end do
! 2.1.2 Mark set of neighboring point to the right to be recalculated.
       do j = i+1, nsetset
          if (ibmagsetset(ithetasort(j)) == idel) cycle
          if (ibmagsetset(ithetasort(j)) == 0) cycle
          if (bmagset(ibmagsetset(ithetasort(j))) == 0.0) cycle
          work(ibmagsetset(ithetasort(j))) = ibmagsetset(ithetasort(j))
          exit
       end do
    end do

! 3 Delete this set.
    bmagset(idel) = 0.0

! 4 Recalculate resolution metric for affected sets.
    do i = 1, nset
       if (work(i) /= 0) then
          call get_thetares (work(i), thetares(work(i)))
          call get_lambdares (work(i), alambdares(work(i)))
       end if
    end do
  end subroutine delete_set

  subroutine gg4results
    implicit none
    integer :: i, n

    n = 0
    do i = 1, nsetset
       if (ibmagsetset(ithetasort(i)) == 0) cycle
       if (bmagset(ibmagsetset(ithetasort(i))) == 0.0) cycle
       n = n + 1
       thetagrid(n) = thetasetset(ithetasort(i))
       bmaggrid(n) = bmagset(ibmagsetset(ithetasort(i)))
    end do

! Check point at +pi
    if (bmaggrid(n) /= bmaggrid(1)) then
       n = n + 1
       thetagrid(n) = pi
       bmaggrid(n) = bmaggrid(1)
    end if
    nthetaout = n

    n = 0
    do i = nset, 1, -1
       if (bmagset(ibmagsort(i)) == 0.0) cycle
       n = n + 1
       bset(n) = bmagset(ibmagsort(i))
    end do
    nlambdaout = n
  end subroutine gg4results

  subroutine gg4debug (n, x, y, label)
    implicit none
    integer, intent (in) :: n
    real, dimension (:), intent (in) :: x, y
    character(*), intent (in) :: label
    !
    integer :: i

    if (debug_output) then
       write (unit=debug_unit, fmt="('#',a)") label
       do i = 1, n
          write (unit=debug_unit, fmt="(i5,x,2(g19.12,x))") i, x(i), y(i)
       end do
    end if
  end subroutine gg4debug

  subroutine gg4debugi (n, x, i1, i2, label)
    implicit none
    integer, intent (in) :: n
    real, dimension (:), intent (in) :: x
    integer, dimension (:), intent (in) :: i1, i2
    character(*), intent (in) :: label
    !
    integer :: i

    if (debug_output) then
       write (unit=debug_unit, fmt="('#',a)") label
       do i = 1, n
          write (unit=debug_unit, fmt="(i5,x,g22.15,x,i4,x,i10)") &
               i, x(i), i1(i), i2(i)
       end do
    end if
  end subroutine gg4debugi

  function bmagint (theta)
    implicit none
    real, intent (in) :: theta
    real :: bmagint
    real, external :: fitp_curvp2

    bmagint = fitp_curvp2(theta,nbmag-1,thetain,bmagin,twopi,bmagspl,tension)
  end function bmagint

  function bmagp (theta)
    implicit none
    real, intent (in) :: theta
    real :: bmagp
    real, parameter :: diff = 1e-5
    !
    real :: th1, th2, bm1, bm2
    th1 = theta + 0.5*diff
    th2 = theta - 0.5*diff
    bm1 = bmagint (th1)
    bm2 = bmagint (th2)
    bmagp = (bm1-bm2)/diff
  end function bmagp

  real function tfunc (theta)
    implicit none
    real :: theta
    tfunc = 1.0 + deltaw*(1.0/(1.0+(theta-thetamax)**2/widthw**2) &
         +1.0/(1.0+(theta+thetamax)**2/widthw**2))
  end function tfunc

end subroutine gridgen4_2
