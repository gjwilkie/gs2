       allocate (   bxs(3*nnx*nny,-ntgrid:ntgrid-1))
       allocate (   bys(3*nnx*nny,-ntgrid:ntgrid-1))
       allocate (   vxs(3*nnx*nny,-ntgrid:ntgrid-1))
       allocate (   vys(3*nnx*nny,-ntgrid:ntgrid-1))

       allocate (bxsavg(3*nnx*nny,-ntgrid:ntgrid-1))
       allocate (bysavg(3*nnx*nny,-ntgrid:ntgrid-1))
       allocate (vxsavg(3*nnx*nny,-ntgrid:ntgrid-1))
       allocate (vysavg(3*nnx*nny,-ntgrid:ntgrid-1))

       allocate (stemp(nnx+2*nny))
       allocate (zx1(nny), zxm(nny), zy1(nnx), zyn(nnx))

       do ig=-ntgrid, ntgrid-1
          call fitp_surf1(nnx, nny, xx, yy, bxf(:,:,ig), &
               nnx, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
               255, bxs(:,ig), stemp, 1., ierr)

          call fitp_surf1(nnx, nny, xx, yy, byf(:,:,ig), &
               nnx, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
               255, bys(:,ig), stemp, 1., ierr)

          call fitp_surf1(nnx, nny, xx, yy, vxf(:,:,ig), &
               nnx, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
               255, vxs(:,ig), stemp, 1., ierr)

          call fitp_surf1(nnx, nny, xx, yy, vyf(:,:,ig), &
               nnx, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
               255, vys(:,ig), stemp, 1., ierr)
       end do

       deallocate (zx1, zxm, zy1, zyn, stemp)

       do ig=-ntgrid, ntgrid-1
          bxsavg(:,ig) = 0.5*(bxs(:,ig)+bxs(:,ig+1))
          bysavg(:,ig) = 0.5*(bys(:,ig)+bys(:,ig+1))
          vxsavg(:,ig) = 0.5*(vxs(:,ig)+vxs(:,ig+1))  ! unnecessary
          vysavg(:,ig) = 0.5*(vys(:,ig)+vys(:,ig+1))  ! unnecessary
       end do

       ! now, integrate to find a field line

       allocate (rx(nnx,nny,-ntgrid,ntgrid-1))
       allocate (ry(nnx,nny,-ntgrid,ntgrid-1))

       do ik=1,nny
          do it=1,nnx
             rx(it,ik,-ntgrid) = xx(it)
             ry(it,ik,-ntgrid) = yy(ik)
          end do
       end do

       do ig=-ntgrid,ntgrid-1
          rxt = rx + 0.5*dz*bxf(ig)  ! need mod or something like that to stay in the box
          ryt = ry + 0.5*dz*byf(ig)  ! need mod or something like that to stay in the box
          
          bxt = bx(rxt, ryt)
          byt = by(rxt, ryt)
          
          rx(ig+1) = rx(ig) + dz*bxt  ! need mod or something like that to stay in the box
          ry(ig+1) = ry(ig) + dz*byt  ! need mod or something like that to stay in the box
       end do

