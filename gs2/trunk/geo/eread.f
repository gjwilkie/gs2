      program eread
      character*20 fname
      integer lnblank
      character*5 char
      real dummy(129),psi(129,129)
      real q(129),p(129)

      write(*,*) 'filename?'
      read(*,*) fname
      fname=fname(1:lnblank(fname,80))

      open(unit=10,file=fname,form='formatted')
      
      read(10,1000)char,char,char,char,char,i,nw,nh
      read(10,2020)rwid,zhei,rcentr,rleft,xdum      
      read(10,2020)R_0,Z_0,psi_0,psi_a,bcentr
      do i=1,2
         read(10,2020)xdum,xdum,xdum,xdum,xdum
      enddo
c fp
      read(10,2020)(dummy(j),j=1,nw)
c p
      read(10,2020)(p(j),j=1,nw)
c ff'
      read(10,2020)(dummy(j),j=1,nw)
c p'
      read(10,2020)(dummy(j),j=1,nw)
c
      read(10,2020)((psi(i,j),i=1,nw),j=1,nh)
c q
      read(10,2020)(q(j),j=1,nw)
      write(*,*) 'q         p'
      do j=1,nw
         write(*,2020) q(j),p(j)
      enddo

c      read(10,2022) nbbbs,limitr
c      read(10,2020)(rbbbs(i),zbbbs(i),i=1,nbbbs)


 1000 format(5(a10),i2,i4,i4)
 2020 format (5e16.9)
 2022 format (2i5)      
      stop
      end




      integer function lnblank(string,max)
c
c Finds the FIRST non-blank character in a string
c
      character string*(*)
      integer i,j,max,max2

      max2=min(max,len(string))  ! don't search past the end of the string.

      i=0
      do j=1,max2
         if( string(j:j) .ne. ' ' ) then
            i = i + 1
         else
            lnblank=i
            return
         endif
      enddo
      lnblank=i
      return 
      end


