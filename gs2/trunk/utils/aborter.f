        Subroutine Aborter(iunit,ierrmsg)
C----------------------------------------------------------------------
C  ABORT A PROGRAM AFTER A FATAL ERROR CONDITION IS DETECTED.
C
c
c input: iunit  Unit Number of the file to write error messages to.
c        ierrmsg  An error message to write ilunerr
c
c       The advantage of using this subroutine is that it will
c       generate a traceback showing the chain of subroutines which
c       eventually bombed, and it forces an arithmetic error which
c       the job control system can detect as a fatal error.
C
        character ierrmsg*(*)
        common /abortcmn/ zz0,zz1
c
c zz0 is in a common block to prevent an optimizing compiler
c from evaluating 1.0/zz0 during compilation rather than during
c execution.
c
        write(iunit,1001)
1001    format(//
     *   ' %ABORTER:  ** FATAL ERROR.  ABORT SUBROUTINE CALLED **'
     *    //)

        write(iunit,1002)ierrmsg
1002    format(1x,a,//)

c on CRAY's and VAXes:
c generate a divide-by-zero error:
        zz1=1.0/zz0         !^
        ilunerr=5           !^

c on the DecStation:
C^      call exit(1)

        write(ilunerr,1010) zz0
1010    format(' ?ABORTER-- ZZ0= ',1PE11.4,' AND STILL EXECUTING...')
        STOP
        END


