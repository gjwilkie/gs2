      program tix2nt
c
c Convert Experimental Ti+Error_bars ufiles into ascii tables for plotting
c and comparisons with NT.
c
      implicit none
      integer mxgrid
      parameter (mxgrid=100)
      integer mtimes
      parameter (mtimes=1000)

      real Ti_r(mxgrid,mtimes)
      real Ti_err_r(mxgrid,mtimes)
      real Rmajp(mxgrid)
      real times(mtimes)

      real dum1, dum2,twant
      integer k,ngrid,kerr,ifail
      integer icount,iargc
      character runname*80
      integer lrunname
      character line*256
      character cshot*60
      character ctime*60
      character udir*80
      integer lshot,nx,ny
c
c Determine the runname for the input and output files.
c If not passed on the command line, then prompt for it:
c
      icount=iargc()
      if(icount .eq. 3) then
         call getarg(1,udir)
         call getarg(2,cshot)
         call getarg(3,ctime)
      else
         write(6,*) 'Usage: tix2nt ufile-directory shot time'
         write(6,*) 'where ufile-directory is the directory where the '
         write(6,*) '                cort*.ti and fsm*.teb ufiles are.'
         write(6,*) '                /u/hammett/vax/transp2u/sw/'
         write(6,*) '      shot is the shot #'
         write(6,*) '      time is the time at which to extract the Ti data'
         write(6,*) ' '
         write(6,*) 'Example: tix2nt /u/hammett/vax/transp2u/sw/ 88582 4.419'
         stop
      endif
      lshot=index(cshot,' ')-1
c will not compile on cray
c      read(ctime,'(f)') twant

c Read the Ti data:
      call u2read(udir,'fsm','teb',cshot,
     .     Rmajp,times,Ti_err_r,mtimes*mxgrid,nx,ny,ifail)

c Read the Ti error bars:
      call u2read(udir,'cort','ti',cshot,
     .     Rmajp,times,Ti_r,mtimes*mxgrid,nx,ny,ifail)

      call writedat('t'//cshot(1:lshot)//'_tix.dat',
     .     Rmajp,times,Ti_r,Ti_err_r,nx,ny,twant)

      stop
      end

      subroutine writedat(fname, xdata, ydata, fdata, fdata_err, nx, ny, ywant)
      implicit none
      character fname*(*)
      integer nx,ny
      real xdata(nx)
      real ydata(ny)
      real fdata(nx,ny)
      real fdata_err(nx,ny)
      real ywant

c local variables:
      integer j,k,jwant

c find the desired time:

      do j=1,ny
         jwant=j
         if(ydata(j) .ge. ywant) goto 30
      enddo
 30   if(j .gt. 1) then
         if(abs(ydata(j-1)-ywant) .lt. abs(ydata(j)-ywant)) jwant=j-1
      endif

      open(unit=21,file=fname,status='unknown')
      write(21,'(a,f8.3,a)') 
     &     '# Rmaj        Ti_exp     Ti_err   at t=',ywant,' sec.'
      write(21,'(a)') '#  (m)         (keV)      (keV)'
      do k=1,nx
         write(21,55) xdata(k)/100., fdata(k,jwant)/1000., 
     .        fdata_err(k,jwant)/1000.
      enddo
 55   format(10(1x,1pe11.4))
      close(unit=21)

      return
      end

