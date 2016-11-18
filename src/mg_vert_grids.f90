module mg_vert_grids

  use mg_mpi
  use mg_tictoc
  use mg_grids
  use mg_mpi_exchange
  use mg_gather
  use mg_namelist

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
    real(kind=rp) :: Arz

    integer(kind=ip) :: i,j,k

    real(kind=rp), parameter :: two  = 2._rp
    real(kind=rp), parameter :: one  = 1._rp
    real(kind=rp), parameter :: hlf  = 0.5_rp
    real(kind=rp), parameter :: qrt  = 0.25_rp
    real(kind=rp), parameter :: zero = 0._rp

    !! fill and coarsen zr and zw

    do lev = 1, nlevs

       if (myrank==0) write(*,*)'   lev=',lev

       nx=grid(lev)%nx
       ny=grid(lev)%ny
       nz=grid(lev)%nz

       if (lev == 1) then ! fill zr,zw

          grid(lev)%zr(1:nz  ,0:ny+1,0:nx+1) = z_r !only 1 extra point ??????????
          grid(lev)%zw(1:nz+1,0:ny+1,0:nx+1) = z_w !only 1 extra point ??????????

       else ! coarsen zr,zw (needed when directly discretizing on coarser grids)

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

          zrc(1:nzc,1:nyc,1:nxc) = 0.125_8 * (      & !only interior points ?????????? 
               zrf(1:nzf  :2,1:nyf  :2,1:nxf  :2) + &
               zrf(1:nzf  :2,2:nyf+1:2,1:nxf  :2) + &
               zrf(1:nzf  :2,1:nyf  :2,2:nxf+1:2) + &
               zrf(1:nzf  :2,2:nyf+1:2,2:nxf+1:2) + &
               zrf(2:nzf+1:2,1:nyf  :2,1:nxf  :2) + &
               zrf(2:nzf+1:2,2:nyf+1:2,1:nxf  :2) + &
               zrf(2:nzf+1:2,1:nyf  :2,2:nxf+1:2) + &
               zrf(2:nzf+1:2,2:nyf+1:2,2:nxf+1:2) )

          zwc(1:nzc,1:nyc,1:nxc) = 0.25_8 * (       &
               zwf(1:nzf:2,1:nyf  :2,1:nxf  :2) + &
               zwf(1:nzf:2,2:nyf+1:2,1:nxf  :2) + &
               zwf(1:nzf:2,1:nyf  :2,2:nxf+1:2) + &
               zwf(1:nzf:2,2:nyf+1:2,2:nxf+1:2) )      
          zwc(nzc+1,1:nyc,1:nxc) = 0.25_8 * (       &
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

          call fill_halo(lev,grid(lev)%zr) ! necesary fill_halo ?????????? 
          call fill_halo(lev,grid(lev)%zw)

       end if

    end do

    !! compute derived qties

    do lev = 1, nlevs

       if (myrank==0) write(*,*)'   lev=',lev

       nx=grid(lev)%nx
       ny=grid(lev)%ny
       nz=grid(lev)%nz

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
!!$    !! Cell widths
!!$    allocate(dxu(0:ny+1,nx+1))
!!$    do i = 1,nx+1
!!$       do j = 0,ny+1
!!$          dxu(j,i) = 0.5_8*(dx(j,i)+dx(j,i-1))
!!$       enddo
!!$    enddo
!!$    allocate(dyv(ny+1,0:nx+1))
!!$    do i = 0,nx+1
!!$       do j = 1,ny+1
!!$          dyv(j,i) = 0.5_8*(dy(j,i)+dy(j-1,i))
!!$       enddo
!!$    enddo
       !!  Areas

       do i = 1,nx+1
          do j = 0,ny+1
             do k = 1,nz
                !Arx(k,j,i) = 0.25_8*(dz(k,j,i)+dz(k,j,i-1))*(dy(j,i)+dy(j,i-1)) 
                Arx(k,j,i) = qrt * & 
                     ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
                     ( dy(j,i) + dy(j,i-1) ) 
             enddo
          enddo
       enddo

       do i = 0,nx+1
          do j = 1,ny+1
             do k = 1,nz
                !Ary(k,j,i) = 0.25_8*(dz(k,j,i)+dz(k,j-1,i))*(dx(j,i)+dx(j-1,i)) 
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
!!$             zy(k,j,i) = 0.5_8*(zr(k,j+1,i)-zr(k,j-1,i))/dy(j,i)
!!$             zx(k,j,i) = 0.5_8*(zr(k,j,i+1)-zr(k,j,i-1))/dx(j,i)
!!$          enddo
!!$       enddo
!!$    enddo
!!$    call fill_halo(1,zy) ! copy interior value in the halo
!!$    call fill_halo(1,zx) ! copy interior value in the halo       

       do i = 0,nx+1
          do j = 0,ny+1
             do k = 1, nz
                zydx(k,j,i) = hlf * (( zr(k,j+1,i  ) - zr(k,j-1,i  ) ) / dy(j,i) ) * dx(j,i)
                zxdy(k,j,i) = hlf * (( zr(k,j  ,i+1) - zr(k,j  ,i-1) ) / dx(j,i) ) * dy(j,i)
             enddo
          enddo
       enddo
!!$    allocate(zyw(nz+1,0:ny+1,0:nx+1))
!!$    allocate(zxw(nz+1,0:ny+1,0:nx+1))
!!$    do i = 1,nx
!!$       do j = 1,ny
!!$          do k = 1,nz+1
!!$             zyw(k,j,i) = 0.5_8*(zw(k,j+1,i)-zw(k,j-1,i))/dy(j,i)
!!$             zxw(k,j,i) = 0.5_8*(zw(k,j,i+1)-zw(k,j,i-1))/dx(j,i)
!!$          enddo
!!$       enddo
!!$    enddo
!!$    call fill_halo(1,zyw) ! copy interior value in the halo
!!$    call fill_halo(1,zxw) ! copy interior value in the halo

       do i = 0,nx+1
          do j = 0,ny+1

             !do k = 1,nz+1
             !cw(k,j,i) = Arz(j,i)/dzw(k,j,i) * (1._8 + zxw(k,j,i)*zxw(k,j,i)+zyw(k,j,i)*zyw(k,j,i))

             Arz = dx(j,i)*dy(j,i)

             k=1
             cw(k,j,i) = ( Arz / dzw(k,j,i) ) * &
                  (one + &
                  ( hlf * (zw(k,j  ,i+1)-zw(k,j  ,i-1)) / dx(j,i) )**2 + &
                  ( hlf * (zw(k,j+1,i  )-zw(k,j-1,i  )) / dy(j,i) )**2 )

             do k = 2,nz
                cw(k,j,i) = ( Arz / dzw(k,j,i) ) * &
                     (one + &
                     ( hlf * (zw(k,j  ,i+1)-zw(k,j  ,i-1)) / dx(j,i) )**2 + &
                     ( hlf * (zw(k,j+1,i  )-zw(k,j-1,i  )) / dy(j,i) )**2 )

             enddo

             k=nz+1
             cw(k,j,i) = ( Arz / dzw(k,j,i) ) * &
                  (one + &
                  ( hlf * (zw(k,j  ,i+1)-zw(k,j  ,i-1)) / dx(j,i) )**2 + &
                  ( hlf * (zw(k,j+1,i  )-zw(k,j-1,i  )) / dy(j,i) )**2 )

          enddo
       enddo

    enddo

  end subroutine set_vert_grids

end module mg_vert_grids
