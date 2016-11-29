module mg_vert_grids

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_grids
  use mg_mpi_exchange
  use mg_gather
  use mg_namelist
  use mg_netcdf_out

  implicit none

contains

  !----------------------------------------
  subroutine set_vert_grids(z_r,z_w)

    real(kind=rp), dimension(:,:,:), pointer, intent(in)  :: z_r,z_w

    integer(kind=ip) :: lev
    integer(kind=ip) :: nx,ny,nz
    integer(kind=ip) :: nxf,nxc
    integer(kind=ip) :: nyf,nyc
    integer(kind=ip) :: nzf,nzc

    real(kind=rp), dimension(:,:,:), pointer :: zr,zw
 
    real(kind=rp), dimension(:,:,:), pointer :: zrf,zrc 
    real(kind=rp), dimension(:,:,:), pointer :: zwf,zwc

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary
    real(kind=rp), dimension(:,:,:), pointer :: zydx,zxdy
    real(kind=rp), dimension(:,:,:), pointer :: cw

    integer(kind=ip) :: i,j,k

    do lev = 1, nlevs

       !! fill and coarsen zr and zw

       if (myrank==0) write(*,*)'   lev=',lev

       nx=grid(lev)%nx
       ny=grid(lev)%ny
       nz=grid(lev)%nz

       if (lev == 1) then ! zr,zw from croco

          grid(lev)%zr(1:nz  ,0:ny+1,0:nx+1) = z_r !only 1 extra point ??????????
          grid(lev)%zw(1:nz+1,0:ny+1,0:nx+1) = z_w !only 1 extra point ??????????

       else               ! coarsen zr,zw (needed when directly discretizing on coarser grids)

          nxf = grid(lev-1)%nx
          nyf = grid(lev-1)%ny
          nzf = grid(lev-1)%nz

          zrf => grid(lev-1)%zr
          zwf => grid(lev-1)%zw

          if (grid(lev)%gather == 1) then
             nxc = nx/grid(lev)%ngx
             nyc = ny/grid(lev)%ngy
             nzc = nz
             allocate(zrc(1:nzc  ,0:nyc+1,0:nxc+1))
             allocate(zwc(1:nzc+1,0:nyc+1,0:nxc+1))
          else
             nxc = nx
             nyc = ny
             nzc = nz
             zrc => grid(lev)%zr
             zwc => grid(lev)%zw
          endif

          zrc(1:nzc,1:nyc,1:nxc) = eighth * (      & !only interior points ?????????? 
               zrf(1:nzf  :2,1:nyf  :2,1:nxf  :2) + &
               zrf(1:nzf  :2,2:nyf+1:2,1:nxf  :2) + &
               zrf(1:nzf  :2,1:nyf  :2,2:nxf+1:2) + &
               zrf(1:nzf  :2,2:nyf+1:2,2:nxf+1:2) + &
               zrf(2:nzf+1:2,1:nyf  :2,1:nxf  :2) + &
               zrf(2:nzf+1:2,2:nyf+1:2,1:nxf  :2) + &
               zrf(2:nzf+1:2,1:nyf  :2,2:nxf+1:2) + &
               zrf(2:nzf+1:2,2:nyf+1:2,2:nxf+1:2) )

          zwc(1:nzc,1:nyc,1:nxc) = qrt * (       &
               zwf(1:nzf:2,1:nyf  :2,1:nxf  :2) + &
               zwf(1:nzf:2,2:nyf+1:2,1:nxf  :2) + &
               zwf(1:nzf:2,1:nyf  :2,2:nxf+1:2) + &
               zwf(1:nzf:2,2:nyf+1:2,2:nxf+1:2) )      
          zwc(nzc+1,1:nyc,1:nxc) = qrt * (       &
               zwf(nzf+1  ,1:nyf  :2,1:nxf  :2) + &
               zwf(nzf+1  ,2:nyf+1:2,1:nxf  :2) + &
               zwf(nzf+1  ,1:nyf  :2,2:nxf+1:2) + &
               zwf(nzf+1  ,2:nyf+1:2,2:nxf+1:2) ) 

          if (grid(lev)%gather == 1) then
             call gather(lev,zrc,grid(lev)%zr)
             call gather(lev,zwc,grid(lev)%zw)
             deallocate(zrc)
             deallocate(zwc)
          endif

       end if

       call fill_halo(lev,grid(lev)%zr) ! special fill_halo of zr
       call fill_halo(lev,grid(lev)%zw) ! special fill_halo of zw

       if (netcdf_output) then
          call write_netcdf(grid(lev)%zr,vname='zr',netcdf_file_name='zr.nc',rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%zw,vname='zw',netcdf_file_name='zw.nc',rank=myrank,iter=lev)
       endif

       !! compute derived qties

       dx => grid(lev)%dx
       dy => grid(lev)%dy
       zr => grid(lev)%zr
       zw => grid(lev)%zw

       dzw => grid(lev)%dzw
       Arx => grid(lev)%Arx
       Ary => grid(lev)%Ary 
       zxdy => grid(lev)%zxdy
       zydx => grid(lev)%zydx
       cw => grid(lev)%cw

!!$    !! Cell heights
!!$    allocate(dz(nz,0:ny+1,0:nx+1))
!!$    do i = 0,nx+1
!!$       do j = 0,ny+1
!!$          do k = 1,nz
!!$             dz(k,j,i) = zw(k+1,j,i)-zw(k,j,i)
!!$          enddo
!!$       enddo
!!$    enddo

       do i = 0,nx+1
          do j = 0,ny+1
             dzw(1,j,i) = zr(1,j,i)-zw(1,j,i) 
             do k = 2,nz
                dzw(k,j,i) = zr(k,j,i)-zr(k-1,j,i) 
             enddo
             dzw(nz+1,j,i) = zw(nz+1,j,i)-zr(nz,j,i) 
          enddo
       enddo

       !!  Areas
       do i = 1,nx+1
          do j = 0,ny+1
             do k = 1,nz
                !Arx(k,j,i) = qrt*(dz(k,j,i)+dz(k,j,i-1))*(dy(j,i)+dy(j,i-1)) 
                Arx(k,j,i) = qrt * & 
                     ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
                     ( dy(j,i) + dy(j,i-1) ) 
             enddo
          enddo
       enddo
       do i = 0,nx+1
          do j = 1,ny+1
             do k = 1,nz
                !Ary(k,j,i) = qrt*(dz(k,j,i)+dz(k,j-1,i))*(dx(j,i)+dx(j-1,i)) 
                Ary(k,j,i) = qrt * &
                     ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
                     ( dx(j,i) + dx(j-1,i) )
             enddo
          enddo
       enddo
!!$    allocate(Arz(0:ny+1,0:nx+1))
!!$    do i = 0,nx+1
!!$       do j = 0,ny+1
!!$          Arz(j,i) = dx(j,i)*dy(j,i)
!!$       enddo
!!$    enddo

!!$    !! Slopes in x- and y-direction defined at rho-points
!!$    allocate(zx(nz,0:ny+1,0:nx+1))
!!$    allocate(zy(nz,0:ny+1,0:nx+1))
!!$    do i = 1,nx
!!$       do j = 1,ny
!!$          do k = 1,nz
!!$             zx(k,j,i) = hlf*(zr(k,j,i+1)-zr(k,j,i-1))/dx(j,i)
!!$             zy(k,j,i) = hlf*(zr(k,j+1,i)-zr(k,j-1,i))/dy(j,i)
!!$          enddo
!!$       enddo
!!$    enddo
!!$    call fill_halo(1,zx) ! copy interior value in the halo     
!!$    call fill_halo(1,zy) ! copy interior value in the halo    

       do i = 0,nx+1        ! because of the special fill_halo of zr, the slopes zxdy and 
          do j = 0,ny+1     ! zydx in the halo are equal to the slopes at the first interior points 
             do k = 1, nz
                zxdy(k,j,i) = hlf * (zr(k,j  ,i+1)-zr(k,j  ,i-1)) / dx(j,i) * dy(j,i)
                zydx(k,j,i) = hlf * (zr(k,j+1,i  )-zr(k,j-1,i  )) / dy(j,i) * dx(j,i)
             enddo
          enddo
       enddo

       do i = 0,nx+1
          do j = 0,ny+1
             do k = 1,nz+1
                cw(k,j,i) = dx(j,i)*dy(j,i) / dzw(k,j,i) * &
                     ( one + &
                     + ( hlf * (zw(k,j  ,i+1)-zw(k,j  ,i-1)) / dx(j,i) )**2 + &
                     + ( hlf * (zw(k,j+1,i  )-zw(k,j-1,i  )) / dy(j,i) )**2 )
             enddo
          enddo
       enddo

       if (netcdf_output) then
          call write_netcdf(grid(lev)%zxdy,vname='zxdy',netcdf_file_name='zxdy.nc',rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%zydx,vname='zydx',netcdf_file_name='zydx.nc',rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%cw,vname='cw',netcdf_file_name='cw.nc',rank=myrank,iter=lev)
       endif

    enddo

  end subroutine set_vert_grids

end module mg_vert_grids
