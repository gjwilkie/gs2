      integer nw_max,nh_max,nbbbs,limitr,mb
      integer bigfac,big_nw,big_nh,big_nw_max,big_nh_max
c
c Apparently nh_max, nw_max, etc. should be odd
c
      parameter(nw_max=129,nh_max=129)
      parameter(bigfac=8)
      parameter(big_nw_max=1+bigfac*(nw_max-1))
      parameter(big_nh_max=1+bigfac*(nh_max-1))
      parameter(mb=100)

      real R_0,Z_0,rwid,rleft,zhei
      real big_R(big_nw_max),big_Z(big_nh_max),big_q(big_nw_max)
      real big_psi(big_nw_max,big_nh_max),big_fp(big_nw_max)
      real big_p(big_nw_max),big_ba(big_nw_max)
      real rbbbs(mb),zbbbs(mb),big_pbar(big_nw_max)

      common/efit1/big_nw,big_nh,big_R,big_Z,big_psi,big_fp,nbbbs,limitr
      common/efit2/R_0,Z_0,rwid,rleft,zhei,big_q,big_pbar,big_p,big_ba
      common/efit3/rbbbs,zbbbs      
