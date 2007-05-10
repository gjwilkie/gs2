      subroutine u2read (dirname,prefix,suffix,shotnum,
     >     xdata,ydata,zdata,max,nx,ny,ifail)
c
c Read a 2D ufile, in default (real*8) precision.
c
      implicit none
      character dirname*(*)
      character prefix*(*)
      character suffix*(*)
      character shotnum*(*)
      integer max
c
      real xdata(max),ydata(max),zdata(max*max)
      integer nx,ny
      integer ifail
c
      integer lnblank
      external lnblank
c
      integer iunit
      character*500 line
      integer id,ip,is,iline,isx
      integer i,j,k
      logical onretry
c
      ifail = 0
      iunit = 1

      onretry = .false.
 10   continue

      id = lnblank(dirname)
      ip = lnblank(prefix)
      is = lnblank(shotnum)
      isx = lnblank(suffix)

      if (dirname(id:id) .ne. '/' .and. dirname(id:id) .ne. ']') then
         line = dirname(1:id) // '/' // prefix(1:ip)
     >        // shotnum(1:is) // '.' // suffix(1:isx)
      else
         line = dirname(1:id) // prefix(1:ip)
     >        // shotnum(1:is) // '.' // suffix(1:isx)
      endif
      iline = lnblank(line)
      open (unit=iunit,file=line(1:iline),status='old',err=1)

c     header
      read (iunit,'(a)',err=2) line
      read (iunit,'(a)',err=2) line
c     scalars
      read (iunit,'(i5)',err=2) is
      do i=1,is*2
         read (iunit,'(a)',err=2) line
      enddo
c     labels
      read (iunit,'(a)',err=2) line
      read (iunit,'(a)',err=2) line
      read (iunit,'(a)',err=2) line
      read (iunit,'(a)',err=2) line

c     number of points
      read (iunit,'(i11)',err=2) nx
      read (iunit,'(i11)',err=2) ny
      if (ny .gt. max .or. nx .gt. max) goto 3

      do i=1,nx*ny
         zdata(i) = 0.0
      enddo

c     read xdata
      do i=1,nx/6
         read (iunit, *, err=2) (xdata(j),j=1+(i-1)*6,i*6)
      enddo
      if ((nx/6)*6 .ne. nx) then
         read (iunit, *, err=2) (xdata(j),j=(nx/6)*6+1,nx)
      endif

c     read ydata
      do i=1,ny/6
         read (iunit, *, err=2) (ydata(j),j=1+(i-1)*6,i*6)
      enddo
      if ((ny/6)*6 .ne. ny) then
         read (iunit, *, err=2) (ydata(j),j=(ny/6)*6+1,ny)
      endif

c     read zdata 1x6e13.6  
      do i=1,nx*ny/6
         read (iunit, '(1x,6e13.6)', err=2)
     >        (zdata(j),j=1+(i-1)*6,i*6)
      enddo
      if ((nx*ny/6)*6 .ne. nx*ny) then
         read (iunit, '(1x,6e13.6)', err=2)
     >        (zdata(j),j=(nx*ny/6)*6+1,nx*ny)
      endif


      close (unit=iunit)
      return

 1    continue
      ifail = 1
      return

 2    continue
      ifail = 2
      close (iunit)
      return

 3    continue
      ifail = 3
      close (iunit)
      return
      end
