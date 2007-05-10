#include <sys/types.h>
#include <sys/jtab.h>
#include <errno.h>
int MPPTIME(double *timelimit, double *timeleft, double *timeused) {
  job_t job;
  int iret;

  if (job_cntl(0, J_GET_ALL, (int)&job) != -1) {
    *timelimit=job.j_mpptimelimit;
    *timeused=job.j_mpptimeused;
    *timeleft=job.j_mpptimelimit-job.j_mpptimeused;
    iret=0;
  } else {
    *timelimit=0.0;
    *timeused=0.0;
    *timeleft=0.0;
    iret=errno;
  }
  return iret;
}

