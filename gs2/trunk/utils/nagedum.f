      subroutine e01sae(n,x,y,f,ct,cg,ifail)
      real x(n),y(n),f(n),ct(7*n),cg(2,n)
      integer m
      parameter (m=8192)
      double precision dx(m),dy(m),df(m),dct(7*m),dcg(2,m)
      integer ifail,i

      if(n.gt.m) then
         write(*,*) 'm too small in dummy NAG routine e01sae in nage.f'
         stop
      endif

      do i=1,n
         dx(i)=x(i)
         dy(i)=y(i)
         df(i)=f(i)
         dcg(1,i)=cg(1,i)
         dcg(2,i)=cg(2,i)
      enddo
      do i=1,7*n
         dct(i)=ct(i)
      enddo
      call e01saf(n,dx,dy,df,dct,dcg,ifail)
      do i=1,n
         cg(1,i)=dcg(1,i)
         cg(2,i)=dcg(2,i)
      enddo
      do i=1,7*n
         ct(i)=dct(i)
      enddo
      return
      end

      subroutine e01sbe(n,x,y,f,ct,cg,vx,vy,vf,ifail)
      real x(n),y(n),f(n),ct(7*n),cg(2,n),vx,vy,vf
      integer m
      parameter (m=8192)
      double precision dx(m),dy(m),df(m),dct(7*m),dcg(2,m),dvx,dvy,dvf

      if(n.gt.m) then
         write(*,*) 'm too small in dummy NAG routine e01sbe in nage.f'
         stop
      endif

      dvx=vx
      dvy=vy
      do i=1,n
         dx(i)=x(i)
         dy(i)=y(i)
         df(i)=f(i)
         dcg(1,i)=cg(1,i)
         dcg(2,i)=cg(2,i)
      enddo
      do i=1,7*n
         dct(i)=ct(i)
      enddo
      call e01sbf(n,dx,dy,df,dct,dcg,dvx,dvy,dvf,ifail)
      vf=dvf
      return
      end
