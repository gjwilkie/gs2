!  Encapsulation of netCDF API parameters and interface blocks. 
!  Derived from "netcdf.inc" in the netCDF 2.4/HDF 4.0.1  distribution.
!  Russ Rew, Unidata, and Robert Pincus, NASA/GSFC

module netcdf
! Netcdf data types
  integer, parameter :: NCBYTE = 1, NCCHAR = 2, NCSHORT = 3, &
                        NCLONG = 4, NCFLOAT = 5, NCDOUBLE = 6
      
!     masks for the struct NC flag field; passed in as 'mode' arg to
!     nccreate and ncopen.
     
  integer, parameter :: & 
    NCRDWR   = 1,   & ! read/write, 0 => readonly 
    NCCREAT  = 2,   & ! in create phase, cleared by ncendef 
    NCEXCL   = 4,   & ! on create destroy existing file 
    NCINDEF  = 8,   & ! in define mode, cleared by ncendef 
    NCNSYNC  = 16,  & ! synchronise whole header on change (X'20')
    NCHSYNC  = 32,  & ! synchronise numrecs on change (X'10')
    NCNDIRTY = 64,  & ! numrecs has changed (X'40')
    NCHDIRTY = 128, & ! header info has changed (X'80')
    NCFILL   = 0,   & ! prefill vars on endef and increase of record (default)
    NCNOFILL = 256, & ! don't fill vars on endef and increase of record (X'100')
    NCLINK   = 32768  ! isa link (X'8000')
      
! 'mode' arguments for nccreate and ncopen
     
  integer, parameter ::  &
    NCNOWRIT = 0, NCWRITE = NCRDWR, &
    NCCLOB = 11, NCNOCLOB = 15
     
  integer, parameter :: &
    NCUNLIM  = 0 ! 'size' argument to ncdimdef for an unlimited dimension
      
  integer, parameter :: &
    NCGLOBAL = 0 ! attribute id to put/get a global attribute
     
! Advisory Maximums
     
  integer, parameter :: &
    MAXNCOP = 32, MAXNCDIM = 32, MAXNCATT = 512, MAXNCVAR = 512, &
    MAXNCNAM = 128, MAXVDIMS = MAXNCDIM ! Not enforced 
      
     
!     Global netcdf error status variable
!     Initialized in error.c
      
  integer, parameter :: &
    NCNOERR  = 0, &  ! No Error 
    NCEBADID = 1, &  ! Too many netcdfs open 
    NCENFILE = 2, &  ! Not a netcdf id 
    NCEEXIST = 3, &  ! netcdf file exists && NCNOCLOB
    NCEINVAL = 4, &  ! Invalid Argument 
    NCEPERM  = 5, &  ! Write to read only 
    NCENOTIN = 6, &  ! Operation not allowed in data mode 
    NCEINDEF = 7, &  ! Operation not allowed in define mode 
    NCECOORD = 8, &  ! Coordinates out of Domain 
    NCEMAXDS = 9, &  ! String match to name in use 
    NCENAME  = 10, & ! MAXNCDIMS exceeded 
    NCENOATT = 11, & ! Attribute not found 
    NCEMAXAT = 12, & ! MAXNCATTRS exceeded 
    NCEBADTY = 13, & ! Not a netcdf data type 
    NCEBADD  = 14, & ! Invalid dimension id 
    NCEUNLIM = 15, & ! NCUNLIMITED in the wrong index 
    NCEMAXVS = 16, & ! MAXNCVARS exceeded 
    NCENOTVR = 17, & ! Variable not found 
    NCEGLOB  = 18, & ! Action prohibited on NCGLOBAL varid 
    NCENOTNC = 19, & ! Not a netcdf file 
    NCESTS   = 20, &
    NCENTOOL = 21, &    
    NCFOOBAR = 32, &
    NCSYSERR = -1
      
!     Global options variable. Used to determine behavior of error handler.
!     Initialized in lerror.c
     
 integer, parameter :: NCFATAL = 1, NCVERBOS = 2

! Functions in the FORTRAN interface

  interface
    integer function nccre(filename, cmode, rcode)
      character (len = *), intent ( in) :: filename
      integer            , intent ( in) :: cmode
      integer            , intent (out) :: rcode
    end function nccre

    integer function ncopn(filename, rwmode, rcode)
      character (len = *), intent ( in) :: filename
      integer            , intent ( in) :: rwmode
      integer            , intent (out) :: rcode
    end function ncopn

    subroutine ncredf (ncid, rcode)
      integer            , intent ( in) :: ncid
      integer            , intent (out) :: rcode
    end subroutine ncredf

    subroutine ncendf (ncid, rcode)
      integer            , intent ( in) :: ncid
      integer            , intent (out) :: rcode
    end subroutine ncendf

    subroutine ncclos (ncid, rcode)
      integer            , intent ( in) :: ncid
      integer            , intent (out) :: rcode
    end subroutine ncclos

    subroutine ncinq (ncid, ndims, nvars, natts, recdim, rcode)
      integer            , intent ( in) :: ncid
      integer            , intent (out) :: ndims, nvars, natts, recdim
      integer            , intent (out) :: rcode
    end subroutine ncinq

    subroutine ncsnc (ncid, rcode)
      integer            , intent ( in) :: ncid
      integer            , intent (out) :: rcode
    end subroutine ncsnc

    subroutine ncabor (ncid, rcode)
      integer            , intent ( in) :: ncid
      integer            , intent (out) :: rcode
    end subroutine ncabor

    integer function ncddef(ncid, dimname, dimsize, rcode)
      integer            , intent ( in) :: ncid, dimsize
      character (len = *), intent ( in) :: dimname
      integer            , intent (out) :: rcode
    end function ncddef

    integer function ncdid(ncid, dimname, rcode)
      integer            , intent ( in) :: ncid
      character (len = *), intent ( in) :: dimname
      integer            , intent (out) :: rcode
    end function ncdid

    subroutine ncdinq (ncid,dimid, dimname,size,rcode)
      integer            , intent ( in) :: ncid, dimid
      character (len = *), intent (out) :: dimname
      integer            , intent (out) :: size
      integer            , intent (out) :: rcode
    end subroutine ncdinq

    subroutine ncdren (ncid,dimid,dimname, rcode)
      integer            , intent ( in) :: ncid, dimid
      character (len = *), intent ( in) :: dimname
      integer            , intent (out) :: rcode
    end subroutine ncdren

    integer function ncvdef(ncid, varname, vartype, nvdims, vdims, rcode)
      integer              , intent ( in) :: ncid, vartype, nvdims
      character (len = *)  , intent ( in) :: varname
      integer, dimension(:), intent ( in) :: vdims
      integer              , intent (out) :: rcode
    end function ncvdef
      
    integer function ncvid(ncid, varname, rcode)
      integer            , intent ( in) :: ncid
      character (len = *), intent ( in) :: varname
      integer            , intent (out) :: rcode
    end function ncvid

    subroutine ncvinq (ncid,varid, varname,datatype,nvdims,vdims,nvatts,rcode)
      integer            , intent ( in) :: ncid, varid
      character (len = *), intent (out) :: varname
      integer            , intent (out) :: datatype, nvdims, nvatts, rcode
      integer, dimension(:), &
			   intent (out) :: vdims
    end subroutine ncvinq

! In some of the following routines, values may legally have any type and any
! dimension up to 7 (the Fortran limit). We haven't worked out a way to 
! treat this in Fortran-90, so the interface blocks have been commented out. 

!   subroutine ncvpt1 (ncid,varid,indices,value, rcode)
!     integer            , intent ( in) :: ncid, varid
!     integer, dimension(:), &
!                          intent ( in) :: indices
!     real               , intent ( in) :: value
!     integer            , intent (out) :: rcode
!   end subroutine ncvpt1

    subroutine ncvp1c (ncid,varid,indices, chval, rcode)
      integer            , intent ( in) :: ncid, varid
      integer, dimension(:), &
			   intent ( in) :: indices
      character          , intent ( in) :: chval
      integer            , intent (out) :: rcode
    end subroutine ncvp1c

!   subroutine ncvgt1 (ncid,varid,indices, value, rcode)
!     integer            , intent ( in) :: ncid, varid
!     integer, dimension(:), &            
!                          intent ( in) :: indices
!     real               , intent (out) :: value
!     integer            , intent (out) :: rcode
!   end subroutine ncvgt1

    subroutine ncvg1c (ncid,varid,indices, chval, rcode)
      integer            , intent ( in) :: ncid, varid
      integer, dimension(:), &            
			   intent ( in) :: indices
      character          , intent (out) :: chval
      integer            , intent (out) :: rcode
    end subroutine ncvg1c

!   subroutine ncvpt (ncid,varid,start,counts,values, rcode)
!     integer            , intent ( in) :: ncid, varid
!     integer, dimension(:), &
!                          intent ( in) :: start, counts
!     real               , intent ( in) :: values(*)
!     integer            , intent (out) :: rcode
!   end subroutine ncvpt

!   subroutine ncvptc (ncid,varid,start,counts,string,lenstr, rcode)
!     integer            , intent ( in) :: ncid, varid
!     integer, dimension(:), &
!                          intent ( in) :: start, counts
!     character (len = *), intent ( in) :: string
!     integer            , intent ( in) :: lenstr
!     integer            , intent (out) :: rcode
!   end subroutine ncvptc

!   subroutine ncvptg (ncid,varid,start,counts,strides,imap,values, rcode)
!     integer            , intent ( in) :: ncid, varid
!     integer, dimension(:), &
!                          intent ( in) :: start, counts, strides, imap
!     real               , intent ( in) :: values(*)
!     integer            , intent (out) :: rcode
!   end subroutine ncvptg

!   subroutine ncvpgc (ncid,varid,start,counts,strides,imap,string,rcode)
!     integer            , intent ( in) :: ncid, varid
!     integer, dimension(:), &
!                          intent ( in) :: start, counts, strides, imap
!     character (len = *), intent ( in) :: string
!     integer            , intent (out) :: rcode
!   end subroutine ncvpgc

!   subroutine ncvgt (ncid,varid,start,counts, values,rcode)
!     integer            , intent ( in) :: ncid, varid
!     integer, dimension(:), &
!                          intent ( in) :: start, counts
!     real               , intent (out) :: values
!     integer            , intent (out) :: rcode
!   end subroutine ncvgt

!   subroutine ncvgtc (ncid,varid,start,counts, string,lenstr,rcode)
!     integer            , intent ( in) :: ncid, varid
!     integer, dimension(:), &
!                          intent ( in) :: start, counts
!     character (len = *), intent (out) :: string
!     integer            , intent ( in) :: lenstr
!     integer            , intent (out) :: rcode
!   end subroutine ncvgtc

!   subroutine ncvgtg (ncid,varid,start,counts,strides,imap,values,rcode)
!     integer            , intent ( in) :: ncid, varid
!     integer, dimension(:), &
!                          intent ( in) :: start, counts, strides, imap
!     real               , intent (out) :: values
!     integer            , intent (out) :: rcode
!   end subroutine ncvgtg

!   subroutine ncvggc (ncid,varid,start,counts,strides,imap,string,rcode)
!     integer            , intent ( in) :: ncid, varid
!     integer, dimension(:), &
!                          intent ( in) :: start, counts
!     integer            , intent ( in) :: strides(*)
!     integer            , intent ( in) :: imap(*)
!     character (len = *), intent (out) :: string
!     integer            , intent (out) :: rcode
!   end subroutine ncvggc

    subroutine ncvren (ncid,varid,varname, rcode)
      integer            , intent ( in) :: ncid, varid
      character (len = *), intent ( in) :: varname
      integer            , intent (out) :: rcode
    end subroutine ncvren

!   subroutine ncapt (ncid,varid,attname,datatype,attlen,values, rcode)
!     integer            , intent ( in) :: ncid, varid
!     character (len = *), intent ( in) :: attname
!     integer            , intent ( in) :: datatype, attlen
!     real               , intent ( in) :: values(*)
!     integer            , intent (out) :: rcode
!   end subroutine ncapt

    subroutine ncaptc (ncid,varid,attname,datatype,lenstr,string, rcode)
      integer            , intent ( in) :: ncid, varid
      character (len = *), intent ( in) :: attname
      integer            , intent ( in) :: datatype, lenstr
      character (len = *), intent ( in) :: string
      integer            , intent (out) :: rcode
    end subroutine ncaptc

    subroutine ncainq (ncid,varid,attname, datatype,attlen,rcode)
      integer            , intent ( in) :: ncid, varid
      character (len = *), intent ( in) :: attname
      integer            , intent (out) :: datatype, attlen
      integer            , intent (out) :: rcode
    end subroutine ncainq

!   subroutine ncagt (ncid,varid,attname, values,rcode)
!     integer            , intent ( in) :: ncid, varid
!     character (len = *), intent ( in) :: attname
!     real               , intent (out) :: values
!     integer            , intent (out) :: rcode
!   end subroutine ncagt

    subroutine ncagtc (ncid,varid,attname, string,lenstr,rcode)
      integer            , intent ( in) :: ncid, varid
      character (len = *), intent ( in) :: attname
      character (len = *), intent (out) :: string
      integer            , intent ( in) :: lenstr
      integer            , intent (out) :: rcode
    end subroutine ncagtc

    subroutine ncacpy (inncid,invarid,attname,outncid,outvarid, rcode)
      integer            , intent ( in) :: inncid, invarid
      character (len = *), intent ( in) :: attname
      integer            , intent ( in) :: outncid, outvarid
      integer            , intent (out) :: rcode
    end subroutine ncacpy

    subroutine ncanam (ncid,varid,attnum, attname,rcode)
      integer            , intent ( in) :: ncid, varid, attnum
      character (len = *), intent (out) :: attname
      integer            , intent (out) :: rcode
    end subroutine ncanam

    subroutine ncaren (ncid,varid,attname,newname, rcode)
      integer            , intent ( in) :: ncid, varid
      character (len = *), intent ( in) :: attname, newname
      integer            , intent (out) :: rcode
    end subroutine ncaren

    subroutine ncadel (ncid,varid,attname, rcode)
      integer            , intent ( in) :: ncid, varid
      character (len = *), intent ( in) :: attname
      integer            , intent (out) :: rcode
    end subroutine ncadel
    
    integer function nctlen(type, rcode)
      integer            , intent ( in) :: type
      integer            , intent (out) :: rcode
    end function nctlen

    subroutine ncpopt (ncopts)
      integer            , intent ( in) :: ncopts
    end subroutine ncpopt

    subroutine ncgopt (ncopts)
      integer            , intent (out) :: ncopts
    end subroutine ncgopt

    integer function ncsfil(ncid, fillmode, rcode)
      integer            , intent ( in) :: ncid, fillmode
      integer            , intent (out) :: rcode
    end function ncsfil

  end interface
end module netcdf
