      subroutine gridgen3read
     >     (unit,ntheta,nlambda,thetagrid,bmaggrid,alambdagrid)
      implicit none
cin
      integer unit
cinout
      integer ntheta, nlambda
cout
      real thetagrid(ntheta), bmaggrid(ntheta), alambdagrid(nlambda)
c
      integer nthetain, nlambdain
      integer i
c
      read (unit,*) nthetain
      if (nthetain .gt. ntheta) then
         write (*,*) 'gridgen3read:nthetain.gt.ntheta:',nthetain,ntheta
         stop
      endif
      do i = 1, nthetain
         read (unit,*) thetagrid(i), bmaggrid(i)
      enddo
      ntheta = nthetain

      read (unit,*) nlambdain
      if (nlambdain .gt. nlambda) then
         write (*,*) 'gridgen3read:nlambdain.gt.nlambda:',
     >        nlambdain,nlambda
         stop
      endif
      do i = 1, nlambdain
         read (unit,*) alambdagrid(i)
      enddo
      nlambda = nlambdain
      return
      end

      subroutine gridgen3
     >     (nbmag,theta,bmag, npadd, tfunc,alknob,epsknob,
     >     ntheta,nlambda,thetagrid,bmaggrid,alambdagrid)
      implicit none
cin
      integer nbmag
      real theta(nbmag), bmag(nbmag)
      integer npadd
      real tfunc
      external tfunc
      real alknob,epsknob
cinout
      integer ntheta
      integer nlambda
cout
      real thetagrid(ntheta)
      real bmaggrid(ntheta)
      real alambdagrid(nlambda)
clocal
      integer nmax
      parameter (nmax=10000)
      real pi, twopi
      parameter (pi=3.1415926535897931, twopi=6.2831853071795862)
      real tol
      parameter (tol=1e-12)
c
      real bmagaux(nmax)
c
      integer nstart
      real thetastart(nmax)
      logical essentialstart(nmax)
c
      integer nsetset
      integer nset
      real thetasetset(nmax)
      integer ibmagsetset(nmax)
      integer icollsetset(nmax)
      real bmagset(nmax)
      logical essentialset(nmax)
c
      integer ithetasort(nmax)
      integer ibmagsort(nmax)
      real thetares(nmax)
      real alambdares(nmax)
c
      integer nthetadelete, nthetaout
      integer nlambdadelete, nlambdaout
c
      integer i, idel, iset, isetset
      logical flag
c
      include 'interpolate.if'
      real gridgen_dspl
      external gridgen_dspl
c
      if (nbmag .ge. nmax) then
         write (*,*) 'gridgen3:0:nmax too small, nbmag,nmax=',nbmag,nmax
         stop
      endif
      call inter_getpspl (nbmag-1, theta, twopi, bmag, bmagaux)
      call gg3debug (nbmag,theta,bmag,'input grid')

c 1 Set up starting grid.
      call gg3start (nbmag,theta,twopi,bmag,bmagaux,npadd,
     >     pi,nmax, nstart,thetastart,essentialstart)
      call gg3debug (nstart,thetastart,thetastart,'starting grid')

c 2 Collect all points with same bmag as starting grid.
      call gg3collect (nbmag,theta,twopi,bmag,bmagaux,
     >     pi,nmax,epsknob, nstart,thetastart,essentialstart,
     >     nsetset,nset,thetasetset,ibmagsetset,icollsetset,
     >     bmagset,essentialset)
      call gg3debug (nset,bmagset,bmagset,'bmagset')
      call gg3debugi
     >     (nsetset,thetasetset,ibmagsetset,icollsetset,'thetasetset')

c 3 Sort theta points and bmag sets.
      call heapsort (nsetset, thetasetset, ithetasort)
      call heapsort (nset, bmagset, ibmagsort)

c 4 Calculate spacing and resolution metrics.
      do i = 1, nset
         call gg3thetares
     >        (nset,nsetset,thetasetset,ithetasort,ibmagsetset,bmagset,
     >         essentialset,tfunc,i, thetares(i))
         call gg3lambdares
     >        (nbmag,theta,twopi,bmag,bmagaux,
     >        nset,bmagset,ibmagsort,i, alambdares(i))

      enddo

c 5 Remove sets until the number of points and sets are small enoough.
      nthetadelete = 0
      nlambdadelete = 0
      goto 2
 1    continue
c 5.1 Find the set with the smallest resolution metric.
      idel = 0
      do i = 1, nset
         if (bmagset(i).ne.0.0) then
            flag = .false.
            if (idel .eq. 0) flag = .true.
            if (.not. flag) then
               if (thetares(i) + alknob*alambdares(i)
     >              .le. thetares(idel) + alknob*alambdares(idel))
     >         then
                  flag = .true.
               endif
            endif
            if (flag) idel = i
         endif
      enddo
c 5.2 Delete the set just found.
      if (idel.eq.0) then
         write (*,*) 'gridgen3.f:gridgen:5.2: this cannot happen'
         stop
      endif
      call gg3delete (nbmag,theta,twopi,bmag,bmagaux,
     >     nset,nsetset,thetasetset,ithetasort,ibmagsetset,
     >     ibmagsort,essentialset,tfunc,idel,
     >     bmagset,thetares,alambdares, nthetadelete)
      nlambdadelete = nlambdadelete + 1

 2    continue
      if (nset - nlambdadelete .gt. nlambda
     >     .or. nsetset - nthetadelete .gt. ntheta)
     >then
         goto 1
      endif

c 6 Build theta and lambda grids.
      i = 0
      do isetset = 1, nsetset
         if (bmagset(ibmagsetset(ithetasort(isetset))) .ne. 0.0) then
            i = i + 1
            thetagrid(i) = thetasetset(ithetasort(isetset))
            bmaggrid(i) = bmagset(ibmagsetset(ithetasort(isetset)))
         endif
      enddo
c 6.1 Check point at +pi.
      if (bmaggrid(i) .ne. bmaggrid(1)) then
         i = i + 1
         thetagrid(i) = pi
         bmaggrid(i) = bmaggrid(1)
      endif
      nthetaout = i
      i = 0
      do iset = nset, 1, -1
         if (bmagset(ibmagsort(iset)) .ne. 0.0) then
            i = i + 1
            alambdagrid(i) = 1.0/bmagset(ibmagsort(iset))
         endif
      enddo
      nlambdaout = i

      call gg3debug (nthetaout,thetagrid,bmaggrid,'output grid')
      call gg3debug (nlambdaout,alambdagrid,alambdagrid,'output lambda')
      ntheta = nthetaout
      nlambda = nlambdaout
      return
      end

      subroutine gg3start
     >     (nbmag,theta,twopi,bmag,bmagaux,npadd,
     >     pi,nmax,nstart,thetastart,essentialstart)
      implicit none
cin
      integer nbmag
      real theta(nbmag),twopi,bmag(nbmag),bmagaux(nbmag)
      integer npadd
      real pi
      integer nmax
cout
      integer nstart
      real thetastart(nmax)
      logical essentialstart(nmax)
c
      include 'interpolate.if'
      real gridgen_dspl
      external gridgen_dspl
clocal
      integer i,j,k
      real bprimel, bprimer, bprime0
      real thetal, thetar, theta0
      logical flag
      integer maxiter
      parameter (maxiter=100)
      real tol
      parameter (tol=1.0e-9)

c 1 Collect essential points: -pi + all local extrema
c 1.1 Collect -pi
      nstart = 1
      thetastart(1) = -pi
      essentialstart(1) = .true.
c 1.2 Collect extrema
      bprimer = gridgen_dspl(theta(2),nbmag-1,theta,twopi,bmag,bmagaux)
      do i = 2, nbmag-2
         bprimel = bprimer
         bprimer =
     >        gridgen_dspl(theta(i+1),nbmag-1,theta,twopi,bmag,bmagaux)
         if (bprimel.eq.0.0) then
            nstart = nstart + 1
            thetastart(nstart) = theta(i)
            essentialstart(nstart) = .true.
         elseif (bprimel.lt.0.0 .and. bprimer.gt.0.0) then
c 1.2.1 Find local minumum
            thetal = theta(i)
            thetar = theta(i+1)
            do j = 1, maxiter
               theta0 = 0.5*(thetal+thetar)
               bprime0 = gridgen_dspl
     >              (theta0,nbmag-1,theta,twopi,bmag,bmagaux)
               if (bprime0.lt.0.0) then
                  thetal = theta0
               elseif (bprime0.gt.0.0) then
                  thetar = theta0
               endif
            enddo
            if (theta0-theta(i) .lt. tol) theta0 = theta(i)
            if (theta(i+1)-theta0 .lt. tol) theta0 = theta(i+1)
            nstart = nstart + 1
            if (nstart .gt. nmax) then
               write (*,*) 'gg3start:1.2.1: nmax too small'
               stop
            endif
            thetastart(nstart) = theta0
            essentialstart(nstart) = .true.
         elseif (bprimel.gt.0.0 .and. bprimer.lt.0.0) then
c 1.2.2 Find local maximum
            thetal = theta(i)
            thetar = theta(i+1)
            do j = 1, maxiter
               theta0 = 0.5*(thetal+thetar)
               bprime0 = gridgen_dspl
     >              (theta0,nbmag-1,theta,twopi,bmag,bmagaux)
               if (bprime0.gt.0.0) then
                  thetal = theta0
               elseif (bprime0.lt.0.0) then
                  thetar = theta0
               endif
            enddo
            if (theta0-theta(i) .lt. tol) theta0 = theta(i)
            if (theta(i+1)-theta0 .lt. tol) theta0 = theta(i+1)
            nstart = nstart + 1
            if (nstart .gt. nmax) then
               write (*,*) 'gg3start:1.2.2: nmax too small'
               stop
            endif
            thetastart(nstart) = theta0
            essentialstart(nstart) = .true.
         endif            
      enddo

c 2 Collect original grid, except for extrema
      do i = 2, nbmag-1
         flag = .true.
         do j = 1, nstart
            if (abs(theta(i)-thetastart(j)) .lt. tol) flag = .false.
         enddo
         if (flag) then
            nstart = nstart + 1
            if (nstart .gt. nmax) then
               write (*,*) 'gg3start:2: nmax too small'
               stop
            endif
            thetastart(nstart) = theta(i)
            essentialstart(nstart) = .false.
         endif
      enddo

c 3 Collect points between original grid points
      do i = 1, nbmag-1
         do k = 1, npadd
            theta0 = theta(i)
     >           + real(k)/real(npadd+1)*(theta(i+1)-theta(i))
            flag = .true.
            do j = 1, nstart
               if (abs(theta0-thetastart(j)) .lt. tol) flag = .false.
            enddo
            if (flag) then
               nstart = nstart + 1
               if (nstart .gt. nmax) then
                  write (*,*) 'gg3start:3: nmax too small'
                  stop
               endif
               thetastart(nstart) = theta0
               essentialstart(nstart) = .false.
            endif
         enddo
      enddo
      return
      end

      subroutine gg3collect (nbmag,theta,twopi,bmag,bmagaux,
     >     pi,nmax,epsilon, nstart,thetastart,essentialstart,
     >     nsetset,nset,thetasetset,ibmagsetset,icollsetset,
     >     bmagset,essentialset)
      implicit none
cin
      integer nbmag
      real theta(nbmag),twopi,bmag(nbmag),bmagaux(nbmag)
      real pi
      integer nmax
      real epsilon
      integer nstart
      real thetastart(nstart)
      logical essentialstart(nstart)
cout
      integer nsetset,nset
      real thetasetset(nmax)
      integer ibmagsetset(nmax)
      integer icollsetset(nmax)
      real bmagset(nmax)
      logical essentialset(nmax)
c
      include 'interpolate.if'
      real gridgen_dspl
      external gridgen_dspl
clocal
      integer iset,i,j,ii,imatch
      real thetai,bmagi
      real thetal, thetar, theta0, bmag0
      real bprimel, bprimer, bprime0
      integer maxiter
      parameter (maxiter=100)

      nset = 0
      nsetset = 0
      do iset = 1, nstart
         thetai = thetastart(iset)
         bmagi = inter_psplint(thetai,nbmag-1,theta,twopi,bmag,bmagaux)

c 1 For extrema, check previous sets to attach to.
         if (essentialstart(iset)) then
            do i = 1, nset
c 1.1 Check if the extremum should be attached to each previous set.
               if (abs(bmagi-bmagset(i)) .lt. epsilon) then
c 1.1.1 If the extrema does belong in this set, eliminate points
c       near the extrema from this set.
                  do j = 1, nsetset
                     if (ibmagsetset(j).eq.i
     >                    .and. abs(thetasetset(j)-thetai).lt.epsilon)
     >               then
                        ibmagsetset(j) = 0
                     endif
                  enddo
c 1.1.2 Attach the extremum to that set.
                  nsetset = nsetset + 1
                  if (nsetset .gt. nmax) then
                     write (*,*) 'gg3collect:1.1.2: nmax too small'
                     stop
                  endif
                  thetasetset(nsetset) = thetai
                  ibmagsetset(nsetset) = i
                  icollsetset(nsetset) = 112
c 1.1.3 Then skip to the next thetastart.
                  goto 2
               endif
            enddo
         endif

c 2 Start a new set.
         nset = nset + 1
         bmagset(nset) = bmagi
         essentialset(nset) = essentialstart(iset)
         nsetset = nsetset + 1
         if (nsetset .gt. nmax) then
            write (*,*) 'gg3collect:2: nmax too small'
            stop
         endif
         thetasetset(nsetset) = thetai
         ibmagsetset(nsetset) = nset
         icollsetset(nsetset) = 2

c 3 Check each original grid interval for matching bmag
         do i = 1, nbmag-1
c 3.0.1 This is really annoying.  Stoopid small inaccuracies near -pi
c       often cause another point to be found near -pi.  Skip the
c       first interval if collecting for the set containing -pi.
            if (nset.eq.1 .and. i.eq.1) goto 1
            if (bmag(i).gt.bmagi .and. theta(i).ne.thetai) then
c 3.1 Consider when the left grid point is greater than the target bmag.
c     Then, there are three cases in which there are matching points.
c     (1) The right grid point is equal to the target bmag, and the
c         slope at the right point is positive.
c     (2) The right grid point is less than the target bmag.
c     (3) The right grid point is greater than the target bmag,
c         and the interval is concave down, and the minimum in
c         this interval, which is guaranteed to exist, is less than
c         the target bmag.
               bprimer = gridgen_dspl
     >              (theta(i+1),nbmag-1,theta,twopi,bmag,bmagaux)
               if ((bmag(i+1).eq.bmagi .or. theta(i+1).eq.thetai)
     >              .and. bprimer.gt.0.0)
     >         then
c 3.1.1 Consider when the right grid point is equal to the target bmag.
                  thetal = theta(i)
                  thetar = theta(i+1)
                  do j = 1, maxiter
                     theta0 = 0.5*(thetal+thetar)
                     bmag0 = inter_psplint
     >                    (theta0,nbmag-1,theta,twopi,bmag,bmagaux)
                     if (bmag0 .gt. bmagi) then
                        thetal = theta0
                     elseif (bmag0 .lt. bmagi) then
                        thetar = theta0
                     endif
                  enddo
                  nsetset = nsetset + 1
                  if (nsetset .gt. nmax) then
                     write (*,*) 'gg3collect:3.1.1: nmax too small'
                     stop
                  endif
                  thetasetset(nsetset) = theta0
                  ibmagsetset(nsetset) = nset
                  icollsetset(nsetset) = 311
                  goto 1
               endif
               if (bmag(i+1) .lt. bmagi) then
c 3.1.2 Consider when right grid point is less than the target bmag.
c 3.1.2.1 If this interval bounds the starting theta point, the
c         target point is the starting theta point, so skip to the
c         next interval.
                  if (thetai.ge.theta(i) .and. thetai.le.theta(i+1))
     >            then
                     goto 1
                  endif
c 3.1.2.2 Otherwise, find and collect the target point.
                  thetal = theta(i)
                  thetar = theta(i+1)
                  do j = 1, maxiter
                     theta0 = 0.5*(thetal+thetar)
                     bmag0 = inter_psplint
     >                    (theta0,nbmag-1,theta,twopi,bmag,bmagaux)
                     if (bmag0 .gt. bmagi) then
                        thetal = theta0
                     elseif (bmag0 .lt. bmagi) then
                        thetar = theta0
                     endif
                  enddo
                  nsetset = nsetset + 1
                  if (nsetset .gt. nmax) then
                     write (*,*) 'gg3collect:3.1.2.2: nmax too small'
                     stop
                  endif
                  thetasetset(nsetset) = theta0
                  ibmagsetset(nsetset) = nset
                  icollsetset(nsetset) = 3122
                  goto 1
               endif
c 3.1.3 Check if the grid interval is concave down.
               bprimel = gridgen_dspl
     >              (theta(i),nbmag-1,theta,twopi,bmag,bmagaux)
               bprimer = gridgen_dspl
     >              (theta(i+1),nbmag-1,theta,twopi,bmag,bmagaux)
c 3.1.3.1 If not, skip to the next interval.
               if (bprimel.ge.0.0 .or. bprimer.le.0.0) goto 1
c 3.1.3.2 Consider the case where the starting theta point is within
c         this interval.
               if (thetai.gt.theta(i) .and. thetai.lt.theta(i+1)) then
c 3.1.3.2.1 If the starting point is an extremum, skip to next interval.
                  if (essentialstart(iset)) goto 1
c 3.1.3.2.2 Otherwise, the other target point is right of the starting
c           point if the slope is negative, and left if positive.
                  bprime0 = gridgen_dspl
     >                 (thetai,nbmag-1,theta,twopi,bmag,bmagaux)
                  if (bprime0 .lt. 0.0) then
                     thetal = thetai
                     thetar = theta(i+1)
                  else
                     thetal = theta(i)
                     thetar = thetai
                  endif
                  do j = 1, maxiter
                     theta0 = 0.5*(thetal+thetar)
                     bmag0 = inter_psplint
     >                    (theta0,nbmag-1,theta,twopi,bmag,bmagaux)
                     if (bmag0 .gt. bmagi) then
                        thetal = theta0
                     elseif (bmag0 .lt. bmagi) then
                        thetar = theta0
                     endif
                  enddo
                  nsetset = nsetset + 1
                  if (nsetset .gt. nmax) then
                     write (*,*) 'gg3collect:3.1.3.2.2: nmax too small'
                     stop
                  endif
                  thetasetset(nsetset) = theta0
                  ibmagsetset(nsetset) = nset
                  icollsetset(nsetset) = 31322
                  goto 1
               endif
c 3.1.3.3 If the minimum within this interval is less than the target bmag,
c         then there are two target points in this interval, one on each side
c         of the minimum.
c 3.1.3.3.1 If this interval is on the edge, there won't be any minimum.
               if (i.eq.1 .or. i.eq.nbmag-1) goto 1
c 3.1.3.3.2 Find the minimum in this interval.
               imatch = 0
               do ii = 1, nstart
                  if (essentialstart(ii)
     >                 .and. theta(i) .lt. thetastart(ii)
     >                 .and. thetastart(ii) .le. theta(i+1))
     >            then
                     if (imatch .ne. 0) then
                        write (*,*) 'gridgen3.f:gg3collect:3.1.3.3.2:',
     >                       ' multiple extrema in interval'
                        write (*,*) 'iset,i,ii,imatch:',iset,i,ii,imatch
                        write (*,*) 'bmagi:',bmagi
                        write (*,*) 'bprimel,bprimer:',bprimel,bprimer
                        write (*,*) 'theta(i):',theta(i)
                        write (*,*) 'theta(i+1):',theta(i+1)
                        write (*,*) 'thetai:',thetai
                        stop
                     endif
                     imatch = ii
                  endif
               enddo
               if (imatch .eq. 0) then
                  write (*,*) 'gridgen3.f:gg3collect:3.1.3.3.2:',
     >                 ' missing extremum'
                  write (*,*) 'iset,i,ii,imatch:',iset,i,ii,imatch
                  write (*,*) 'bmagi:',bmagi
                  write (*,*) 'bprimel,bprimer:',bprimel,bprimer
                  write (*,*) 'theta(i):',theta(i)
                  write (*,*) 'theta(i+1):',theta(i+1)
                  write (*,*) 'thetai:',thetai
                  stop
               endif
c 3.1.3.3.2.1 If the minimum is greater than the target bmag, skip to
c             the next interval.
               bmag0 = inter_psplint
     >            (thetastart(imatch),nbmag-1,theta,twopi,bmag,bmagaux)
               if (bmag0 .ge. bmagi) goto 1
c 3.1.3.3.2.2 Collect the point left of the minimum.
               thetal = theta(i)
               thetar = thetastart(imatch)
               do j = 1, maxiter
                  theta0 = 0.5*(thetal+thetar)
                  bmag0 = inter_psplint
     >                 (theta0,nbmag-1,theta,twopi,bmag,bmagaux)
                  if (bmag0 .gt. bmagi) then
                     thetal = theta0
                  elseif (bmag0 .lt. bmagi) then
                     thetar = theta0
                  endif
               enddo
               nsetset = nsetset + 1
               if (nsetset .gt. nmax) then
                  write (*,*) 'gg3collect:3.1.3.3.2.2: nmax too small'
                  stop
               endif
               thetasetset(nsetset) = theta0
               ibmagsetset(nsetset) = nset
               icollsetset(nsetset) = 313322
c 3.1.3.3.2.3 Collect the point right of the minimum.
               thetal = thetastart(imatch)
               thetar = theta(i+1)
               do j = 1, maxiter
                  theta0 = 0.5*(thetal+thetar)
                  bmag0 = inter_psplint
     >                 (theta0,nbmag-1,theta,twopi,bmag,bmagaux)
                  if (bmag0 .lt. bmagi) then
                     thetal = theta0
                  elseif (bmag0 .gt. bmagi) then
                     thetar = theta0
                  endif
               enddo
               nsetset = nsetset + 1
               if (nsetset .gt. nmax) then
                  write (*,*) 'gg3collect:3.1.3.3.2.3: nmax too small'
                  stop
               endif
               thetasetset(nsetset) = theta0
               ibmagsetset(nsetset) = nset
               icollsetset(nsetset) = 313323
               goto 1
            elseif (bmag(i).lt.bmagi .and. theta(i).ne.thetai) then
c 3.2 Consider when the left grid point is less than the target bmag.
c     Then, there are three cases in which there are matching points.
c     (1) The right grid point is equal to the target bmag, and the
c         slope at the right point is negative.
c     (2) The right grid point is greater than the target bmag.
c     (3) The right grid point is less thant the target bmag,
c         and the interval is concave up, and the maximum in
c         this interval, which is guaranteed to exists, is greater
c         thant the target bmag.
               bprimer = gridgen_dspl
     >              (theta(i+1),nbmag-1,theta,twopi,bmag,bmagaux)
c 3.2.1 Consider when the right grid point is equal to the target bmag.
               if ((bmag(i+1).eq.bmagi .or. theta(i+1).eq.thetai)
     >              .and. bprimer.lt.0.0)
     >         then
                  thetal = theta(i)
                  thetar = theta(i+1)
                  do j = 1, maxiter
                     theta0 = 0.5*(thetal+thetar)
                     bmag0 = inter_psplint
     >                    (theta0,nbmag-1,theta,twopi,bmag,bmagaux)
                     if (bmag0 .lt. bmagi) then
                        thetal = theta0
                     elseif (bmag0 .gt. bmagi) then
                        thetar = theta0
                     endif
                  enddo
                  nsetset = nsetset + 1
                  if (nsetset .gt. nmax) then
                     write (*,*) 'gg3collect:3.2.1: nmax too small'
                     stop
                  endif
                  thetasetset(nsetset) = theta0
                  ibmagsetset(nsetset) = nset
                  icollsetset(nsetset) = 321
                  goto 1
               endif
               if (bmag(i+1) .gt. bmagi) then
c 3.2.2 Consider when right grid point is greater than the target bmag.
c 3.2.2.1 If this interval bounds the starting theta point, the
c         target point is the starting theta point, so skip to the
c         next interval.
                  if (thetai.ge.theta(i) .and. thetai.le.theta(i+1))
     >            then
                     goto 1
                  endif
c 3.2.2.2 Otherwise, find and collect the target point.
                  thetal = theta(i)
                  thetar = theta(i+1)
                  do j = 1, maxiter
                     theta0 = 0.5*(thetal+thetar)
                     bmag0 = inter_psplint
     >                    (theta0,nbmag-1,theta,twopi,bmag,bmagaux)
                     if (bmag0 .lt. bmagi) then
                        thetal = theta0
                     elseif (bmag0 .gt. bmagi) then
                        thetar = theta0
                     endif
                  enddo
                  nsetset = nsetset + 1
                  if (nsetset .gt. nmax) then
                     write (*,*) 'gg3collect:3.2.2.2: nmax too small'
                     stop
                  endif
                  thetasetset(nsetset) = theta0
                  ibmagsetset(nsetset) = nset
                  icollsetset(nsetset) = 3222
                  goto 1
               endif
c 3.2.3 Check if the grid interval is concave up.
               bprimel = gridgen_dspl
     >              (theta(i),nbmag-1,theta,twopi,bmag,bmagaux)
               bprimer = gridgen_dspl
     >              (theta(i+1),nbmag-1,theta,twopi,bmag,bmagaux)
c 3.2.3.1 If not, skip to the next interval.
               if (bprimel.le.0.0 .or. bprimer.ge.0.0) goto 1
c 3.2.3.2 Consider the case where the starting theta point is within
c         this interval.
               if (thetai.gt.theta(i) .and. thetai.lt.theta(i+1)) then
c 3.2.3.2.1 If the starting point is an extremum, skip to next interval.
                  if (essentialstart(iset)) goto 1
c 3.2.3.2.2 Otherwise, the other target point is right of the starting
c           point if the slope is positive, and left if negative.
                  bprime0 = gridgen_dspl
     >                 (thetai,nbmag-1,theta,twopi,bmag,bmagaux)
                  if (bprime0 .gt. 0.0) then
                     thetal = thetai
                     thetar = theta(i+1)
                  else
                     thetal = theta(i)
                     thetar = thetai
                  endif
                  do j = 1, maxiter
                     theta0 = 0.5*(thetal+thetar)
                     bmag0 = inter_psplint
     >                    (theta0,nbmag-1,theta,twopi,bmag,bmagaux)
                     if (bmag0 .gt. bmagi) then
                        thetal = theta0
                     elseif (bmag0 .lt. bmagi) then
                        thetar = theta0
                     endif
                  enddo
                  nsetset = nsetset + 1
                  if (nsetset .gt. nmax) then
                     write (*,*) 'gg3collect:3.2.3.2.2: nmax too small'
                     stop
                  endif
                  thetasetset(nsetset) = theta0
                  ibmagsetset(nsetset) = nset
                  icollsetset(nsetset) = 32322
                  goto 1
               endif
c 3.2.3.3 If the maximum within this interval is greater than the target bmag,
c         then there are two target points in this interval, one on each side
c         of the maximum.
c 3.2.3.3.1 If this interval is on the edge, there won't be any maximum.
               if (i.eq.1 .or. i.eq.nbmag-1) goto 1
c 3.2.3.3.2 Find the maximum in this interval.
               imatch = 0
               do ii = 1, nstart
                  if (essentialstart(ii)
     >                 .and. theta(i) .lt. thetastart(ii)
     >                 .and. thetastart(ii) .le. theta(i+1))
     >            then
                     if (imatch .ne. 0) then
                        write (*,*) 'gridgen3.f:gg3collect:3.2.3.3.2:',
     >                       ' multiple extrema in interval'
                        write (*,*) 'iset,i,ii,imatch:',iset,i,ii,imatch
                        write (*,*) 'bmagi:',bmagi
                        write (*,*) 'bprimel,bprimer:',bprimel,bprimer
                        write (*,*) 'theta(i):',theta(i)
                        write (*,*) 'theta(i+1):',theta(i+1)
                        write (*,*) 'thetai:',thetai
                        stop
                     endif
                     imatch = ii
                  endif
               enddo
               if (imatch .eq. 0) then
                  write (*,*) 'gridgen3.f:gg3collect:3.2.3.3.2:',
     >                 ' missing extremum'
                  write (*,*) 'iset,i,ii,imatch:',iset,i,ii,imatch
                  write (*,*) 'bmagi:',bmagi
                  write (*,*) 'bprimel,bprimer:',bprimel,bprimer
                  write (*,*) 'theta(i):',theta(i)
                  write (*,*) 'theta(i+1):',theta(i+1)
                  write (*,*) 'thetai:',thetai
                  stop
               endif
c 3.2.3.3.2.1 If the maximum is less than the target bmag, skip to
c             the next interval.
               bmag0 = inter_psplint
     >            (thetastart(imatch),nbmag-1,theta,twopi,bmag,bmagaux)
               if (bmag0 .le. bmagi) goto 1
c 3.2.3.3.2.2 Collect the point left of the maximum.
               thetal = theta(i)
               thetar = thetastart(imatch)
               do j = 1, maxiter
                  theta0 = 0.5*(thetal+thetar)
                  bmag0 = inter_psplint
     >                 (theta0,nbmag-1,theta,twopi,bmag,bmagaux)
                  if (bmag0 .lt. bmagi) then
                     thetal = theta0
                  elseif (bmag0 .gt. bmagi) then
                     thetar = theta0
                  endif
               enddo
               nsetset = nsetset + 1
               if (nsetset .gt. nmax) then
                  write (*,*) 'gg3collect:3.2.3.3.2.2: nmax too small'
                  stop
               endif
               thetasetset(nsetset) = theta0
               ibmagsetset(nsetset) = nset
               icollsetset(nsetset) = 323322
c 3.2.3.3.2.3 Collect the point right of the maximum.
               thetal = thetastart(imatch)
               thetar = theta(i+1)
               do j = 1, maxiter
                  theta0 = 0.5*(thetal+thetar)
                  bmag0 = inter_psplint
     >                 (theta0,nbmag-1,theta,twopi,bmag,bmagaux)
                  if (bmag0 .gt. bmagi) then
                     thetal = theta0
                  elseif (bmag0 .lt. bmagi) then
                     thetar = theta0
                  endif
               enddo
               nsetset = nsetset + 1
               if (nsetset .gt. nmax) then
                  write (*,*) 'gg3collect:3.2.3.3.2.3: nmax too small'
                  stop
               endif
               thetasetset(nsetset) = theta0
               ibmagsetset(nsetset) = nset
               icollsetset(nsetset) = 323323
               goto 1
            elseif (bmag(i).eq.bmagi .or. theta(i).eq.thetai) then
c 3.3 Consider when the left grid point is equal to the target bmag.
c 3.3.1 Add the point if it isn't the starting theta point.
               if (thetai .ne. theta(i)) then
                  nsetset = nsetset + 1
                  if (nsetset .gt. nmax) then
                     write (*,*) 'gg3collect:3.3.1: nmax too small'
                     stop
                  endif
                  thetasetset(nsetset) = theta(i)
                  ibmagsetset(nsetset) = nset
                  icollsetset(nsetset) = 331
               endif
c 3.3.2 Check if there is another matching target bmag in the interval.
               bprime0 = gridgen_dspl
     >              (theta(i),nbmag-1,theta,twopi,bmag,bmagaux)
               if (bprime0.gt.0.0 .and. bmag(i+1).lt.bmagi
     >              .and. i.ne.1 .and. i.ne.nbmag-1)
     >         then
                  thetal = theta(i)
                  thetar = theta(i+1)
                  do j = 1, maxiter
                     theta0 = 0.5*(thetal+thetar)
                     bmag0 = inter_psplint
     >                    (theta0,nbmag-1,theta,twopi,bmag,bmagaux)
                     if (bmag0 .gt. bmagi) then
                        thetal = theta0
                     elseif (bmag0 .lt. bmagi) then
                        thetar = theta0
                     endif
                  enddo
                  nsetset = nsetset + 1
                  if (nsetset .gt. nmax) then
                     write (*,*) 'gg3collect:3.3.2.1: nmax too small'
                     stop
                  endif
                  thetasetset(nsetset) = theta0
                  ibmagsetset(nsetset) = nset
                  icollsetset(nsetset) = 3321
                  goto 1
               elseif (bprime0.lt.0.0 .and. bmag(i+1).gt.bmagi
     >              .and. i.ne.1 .and. i.ne.nbmag-1)
     >         then
                  thetal = theta(i)
                  thetar = theta(i+1)
                  do j = 1, maxiter
                     theta0 = 0.5*(thetal+thetar)
                     bmag0 = inter_psplint
     >                    (theta0,nbmag-1,theta,twopi,bmag,bmagaux)
                     if (bmag0 .lt. bmagi) then
                        thetal = theta0
                     elseif (bmag0 .gt. bmagi) then
                        thetar = theta0
                     endif
                  enddo
                  nsetset = nsetset + 1
                  if (nsetset .gt. nmax) then
                     write (*,*) 'gg3collect:3.3.2.2: nmax too small'
                     stop
                  endif
                  thetasetset(nsetset) = theta0
                  ibmagsetset(nsetset) = nset
                  icollsetset(nsetset) = 3322
                  goto 1
               endif
            endif
c 4 Ready for next original grid interval.
 1          continue
         enddo
c 5 Ready for next thetastart set.
 2       continue
      enddo
      return
      end

      subroutine gg3thetares
     >     (nset,nsetset,thetasetset,ithetasort,ibmagsetset,bmagset,
     >     essentialset,tfunc,iset, thetares)
      implicit none
cin
      integer nset, nsetset
      real thetasetset(nsetset)
      integer ithetasort(nsetset),ibmagsetset(nsetset)
      real bmagset(nset)
      logical essentialset(nset)
      real tfunc
      external tfunc
      integer iset
cout
      real thetares
clocal
      integer i, ileft, iright
      real dthetal, dthetar
      real res
      integer npts
c
c 1 If this set contains an odd number of points, return large value.
      if (essentialset(iset)) then
c 1.1 The point at -pi is always kept, and always returns a large value.
         if (iset.eq.1) then
            thetares = 2e20
            return
         endif
c 1.2 Count points
         npts = 0
         do i = 1, nsetset
            if (ibmagsetset(i).eq.iset) then
               npts = npts + 1
            endif
         enddo
         if (mod(npts,2) .eq. 1) then
            thetares = 1e20
            return
         endif
      endif

c 2 Look through all points for points in the set.
      res = 0.0
      npts = 0
      do i = 1, nsetset
         if (ibmagsetset(ithetasort(i)) .eq. iset) then
c 2.1 Look for the point to the left.
            do ileft = i-1, 1, -1
               if (bmagset(ibmagsetset(ithetasort(ileft))) .ne. 0.0)
     >         then
                  goto 1
               endif
            enddo
c 2.1.1 If there is none, this is the left boundary point, skip.
            goto 3
 1          continue
c 2.2 Look for the point to the right.
            do iright = i+1, nsetset
               if (bmagset(ibmagsetset(ithetasort(iright))) .ne. 0.0)
     >         then
                  goto 2
               endif
            enddo
c 2.2.1 If there is none, this is the right boundary point, skip.
            goto 3
 2          continue
c 2.3 Add contribution from this interval.
            dthetal = abs(thetasetset(ithetasort(i))
     >           -thetasetset(ithetasort(ileft)))
            dthetar = abs(thetasetset(ithetasort(i))
     >           -thetasetset(ithetasort(ileft)))
            res = res + 1.0/tfunc(thetasetset(ithetasort(i)))
     >           *dthetal*dthetar/(dthetal+dthetar+1e-20)
            npts = npts + 1
         endif
c 2.4 Ready to check next point.
 3       continue
      enddo
      if (npts.gt.0) then
         thetares = res/real(npts)
      else
         thetares = res
      endif
      return
      end

      subroutine gg3lambdares
     >     (nbmag,theta,twopi,bmag,bmagaux,
     >     nset,bmagset,ibmagsort,iset, alambdares)
      implicit none
cin
      integer nbmag
      real theta(nbmag),twopi,bmag(nbmag),bmagaux(nbmag)
      integer nset
      real bmagset(nset)
      integer ibmagsort(nset)
      integer iset
cout
      real alambdares
clocal
      integer i,iplus,iminus
      real al,alplus,alminus,dalplus,dalminus
      real res
      integer npts
c
      al = 1.0/bmagset(iset)
      do i = 1, nset
c 1 Look for target lambda.
         if (ibmagsort(i) .eq. iset) then
c 2 Look for bordering lambdas.
            do iplus = i-1, 1, -1
               if (bmagset(ibmagsort(iplus)) .ne. 0.0) then
                  alplus = 1.0/bmagset(ibmagsort(iplus))
                  goto 1
               endif
            enddo
            alplus = 0.0
 1          continue
            do iminus = i+1, nset
               if (bmagset(ibmagsort(iminus)) .ne. 0.0) then
                  alminus = 1.0/bmagset(ibmagsort(iminus))
                  goto 2
               endif
            enddo
            alminus = 0.0
 2          continue
         endif
      enddo

c 3 Add up contributions to the result.
      res = 0.0
      npts = 0
      do i = 1, nbmag
         dalplus = abs(sqrt(max(1.0-alplus*bmag(i),0.0))
     >        -sqrt(max(1.0-al*bmag(i),0.0)))
         dalminus = abs(sqrt(max(1.0-alminus*bmag(i),0.0))
     >        -sqrt(max(1.0-al*bmag(i),0.0)))
         if (dalplus+dalminus .ne. 0.0) then
            npts = npts + 1
            res = res + dalplus*dalminus/(dalplus+dalminus+1e-20)
         endif
      enddo
      if (npts .ne. 0) then
         alambdares = res/real(npts)
      else
         alambdares = res
      endif
      return
      end

      subroutine gg3delete
     >     (nbmag,theta,twopi,bmag,bmagaux,
     >     nset,nsetset,thetasetset,ithetasort,ibmagsetset,
     >     ibmagsort,essentialset,tfunc,iset,
     >     bmagset,thetares,alambdares, ndeleted)
      implicit none
cin
      integer nbmag
      real theta(nbmag),twopi,bmag(nbmag),bmagaux(nbmag)
      integer nset,nsetset
      real thetasetset(nsetset)
      integer ithetasort(nsetset),ibmagsetset(nsetset)
      integer ibmagsort(nset)
      logical essentialset(nset)
      real tfunc
      external tfunc
      integer iset
cinout
      real bmagset(nset)
      real thetares(nset),alambdares(nset)
cout
      integer ndeleted
clocal
      integer nmax
      parameter (nmax=5000)
      integer nsetrecalc, isetrecalc(nmax)
      integer i,j,k
c
      nsetrecalc = 0
c 1 Mark neighboring lambda sets to have their resolution metrics
c   recalculated.
      do i = 1, nset
         if (ibmagsort(i) .eq. iset) then
            do j = i-1, 1, -1
               if (bmagset(ibmagsort(j)) .ne. 0.0) then
                  nsetrecalc = nsetrecalc + 1
                  isetrecalc(nsetrecalc) = ibmagsort(j)
                  goto 1
               endif
            enddo
 1          continue
            do j = i+1, nset
               if (bmagset(ibmagsort(j)) .ne. 0.0) then
                  nsetrecalc = nsetrecalc + 1
                  isetrecalc(nsetrecalc) = ibmagsort(j)
                  goto 2
               endif
            enddo
 2          continue
         endif
      enddo

c 2 Mark lambda sets with neighboring theta points to have their
c   resolution metrics recalculated.
      ndeleted = 0
      do i = 1, nsetset
c 2.1 Count already deleted points.
         if (ibmagsetset(ithetasort(i)) .eq. 0) then
            ndeleted = ndeleted + 1
         elseif (bmagset(ibmagsetset(ithetasort(i))) .eq. 0.0) then
            ndeleted = ndeleted + 1
c 2.2 Pick out theta points in the set to be deleted.
         elseif (ibmagsetset(ithetasort(i)) .eq. iset) then
            ndeleted = ndeleted + 1
c 2.2.1 Pick out the neighboring theta point to the left.
            do j = i-1, 1, -1
               if (bmagset(ibmagsetset(ithetasort(j))) .ne. 0.0) then
c 2.2.1.1 If the neighboring point is in an already marked set,
c         skip to considering the neighboring theta point to the right.
                  do k = 1, nsetrecalc
                     if (isetrecalc(k) .eq. ibmagsetset(ithetasort(j)))
     >               then
                        goto 3
                     endif
                  enddo
                  if (ibmagsetset(ithetasort(j)) .eq. iset) goto 3
                  nsetrecalc = nsetrecalc + 1
                  isetrecalc(nsetrecalc) = ibmagsetset(ithetasort(j))
                  goto 3
               endif
            enddo
 3          continue
c 2.2.2 Pick out the neighboring theta point to the right.
            do j = i+1, nsetset
               if (bmagset(ibmagsetset(ithetasort(j))) .ne. 0.0) then
c 2.2.2.1 If the neighboring point is in an already marked set,
c         skip to picking out the next theta point.
                  do k = 1, nsetrecalc
                     if (isetrecalc(k) .eq. ibmagsetset(ithetasort(j)))
     >               then
                        goto 4
                     endif
                  enddo
                  if (ibmagsetset(ithetasort(j)) .eq. iset) goto 4
                  nsetrecalc = nsetrecalc + 1
                  isetrecalc(nsetrecalc) = ibmagsetset(ithetasort(j))
                  goto 4
               endif
            enddo
 4          continue
         endif
      enddo

c 3 Delete this set.
      bmagset(iset) = 0.0

c 4 Recalculate resolution metrics for the affected sets.
      do i = 1, nsetrecalc
         call gg3thetares
     >        (nset,nsetset,thetasetset,ithetasort,ibmagsetset,bmagset,
     >        essentialset,tfunc,isetrecalc(i), thetares(isetrecalc(i)))
         call gg3lambdares
     >        (nbmag,theta,twopi,bmag,bmagaux,
     >        nset,bmagset,ibmagsort,isetrecalc(i),
     >        alambdares(isetrecalc(i)))
      enddo
      return
      end

      subroutine gg3debug (n,theta,bmag,caption)
      implicit none
cin
      integer n
      real theta(n), bmag(n)
      character*(*) caption
c
c     call gg3dump (n,theta,bmag,caption)
      return
      end

      subroutine gg3dump (n,theta,bmag,caption)
      implicit none
cin
      integer n
      real theta(n), bmag(n)
      character*(*) caption
clocal
      integer i
c
      write (200, '(''#'',a)') caption
      do i=1,n
         write (200, '(i4,x,2(g19.12,x))') i, theta(i), bmag(i)
      enddo
      return
      end

      subroutine gg3debugi (n,theta,ibmag1,ibmag2,caption)
      implicit none
cin
      integer n
      real theta(n)
      integer ibmag1(n),ibmag2(n)
      character*(*) caption
c
c     call gg3dumpi (n,theta,ibmag1,ibmag2,caption)
      return
      end

      subroutine gg3dumpi (n,theta,ibmag1,ibmag2,caption)
      implicit none
cin
      integer n
      real theta(n)
      integer ibmag1(n),ibmag2(n)
      character*(*) caption
clocal
      integer i
c
      write (200, '(''#'',a)') caption
      do i=1,n
         write (200, '(i4,x,g22.15,x,i4,x,i10)')
     >        i,theta(i),ibmag1(i),ibmag2(i)
      enddo
      return
      end

      real function gridgen_dspl (x0,n,x,p,y,y2)
      implicit none
      real x0
      integer n
      real x(n),p,y(n),y2(n)
c
      real xa,xb,ya,yb
      real diff
      parameter (diff=1e-5)
c
      include 'interpolate.if'
c
      xa = x0 + 0.5*diff
      xb = x0 - 0.5*diff
      ya = inter_psplint(xa,n,x,p,y,y2)
      yb = inter_psplint(xb,n,x,p,y,y2)
      gridgen_dspl = (ya-yb)/diff
      return
      end

      subroutine heapsort (n, v, index)
      implicit none
cin
      integer n
      real v(n)
cout
      integer index(n)
clocal
      integer i, j, l, ir, ira
c
      do i=1,n
         index(i) = i
      enddo

      l = n/2+1
      ir = n

 1    continue
      if (l .gt. 1) then
         l = l - 1
         ira = index(l)
      else
         ira = index(ir)
         index(ir) = index(1)
         ir = ir - 1
         if (ir .eq. 1) then
            index(1) = ira
            return
         endif
      endif
      i = l
      j = l + l
 2    continue
      if (j .le. ir) then
         if (j .lt. ir) then
            if (v(index(j)) .lt. v(index(j+1))) j = j + 1
         endif
         if (v(ira) .lt. v(index(j))) then
            index(i) = index(j)
            i = j
            j = j + j
         else
            j = ir + 1
         endif
         goto 2
      endif
      index(i) = ira
      goto 1
      end
