      integer function system(command)
c
c The "system" call on Sun workstations is called "ishell" on the cray:
c
      character command*(*)
      istat=ishell(command)
      system=istat
      return
      end

c Convert double precision names on workstations to single precision
c on the Crays:

C LAPACK routines:
      subroutine dgbtrf(m,n,kl,ku,ab,ldab,ipiv,info)
      implicit none
      integer info,kl,ku,ldab,m,n
      integer ipiv(*)
      real ab(ldab,*)
      call sgbtrf(m,n,kl,ku,ab,ldab,ipiv,info)
      return
      end

      subroutine dgbtrs(TRANS, N, KL, KU, NRHS, AB,  LDAB,  IPIV,
     &     B, LDB, INFO )
      implicit none
      CHARACTER      TRANS
      INTEGER        INFO, KL, KU, LDAB, LDB, N, NRHS
      INTEGER        IPIV( * )
      real AB( LDAB, * ), B( LDB, * )
      call sgbtrs(TRANS, N, KL, KU, NRHS, AB,  LDAB,  IPIV,
     &     B, LDB, INFO )
      return
      end
