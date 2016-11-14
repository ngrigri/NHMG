module mg_set_bbc

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine set_bbc(zr,zw,u,v,w)

    real(kind=rp), dimension(:,:,:), pointer, intent(in)  :: zr,zw
    real(kind=rp), dimension(:,:,:), pointer, intent(in)  :: u,v
    real(kind=rp), dimension(:,:,:), pointer, intent(inout)  :: w

    integer(kind=ip):: k,j,i
    integer(kind=ip):: nx,ny,nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
!    real(kind=rp), dimension(:,:,:), pointer :: zr,zw
!    real(kind=rp), dimension(:,:,:), pointer :: dzw
!    real(kind=rp), dimension(:,:,:), pointer :: zydx,zxdy
!    real(kind=rp), dimension(:,:,:), pointer :: cw

    real(kind=rp), dimension(:,:)  ,   pointer :: dxu,dyv
    real(kind=rp), dimension(:,:,:),   pointer :: Arx,Ary
    real(kind=rp), dimension(:,:)  ,   pointer :: Arz
    real(kind=rp), dimension(:,:,:),   pointer :: dz,dzw
    real(kind=rp), dimension(:,:,:),   pointer :: zy,zx
    real(kind=rp), dimension(:,:,:),   pointer :: zydx,zxdy
    real(kind=rp), dimension(:,:,:),   pointer :: zxw,zyw
    real(kind=rp), dimension(:,:,:),   pointer :: cw

    real(kind=rp), parameter :: two  = 2._rp
    real(kind=rp), parameter :: one  = 1._rp
    real(kind=rp), parameter :: hlf  = 0.5_rp
    real(kind=rp), parameter :: qrt  = 0.25_rp
    real(kind=rp), parameter :: zero = 0._rp

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx => grid(1)%dx
    dy => grid(1)%dy

    !! Cell heights
    allocate(dz(nz,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz
!             dz(k,j,i) = zw(k+1,j,i)-zw(k,j,i)
             dz(k,j,i) = zw(k,j,i)-zw(k-1,j,i) !because zw indexed as croco from 0 to nz
          enddo
       enddo
    enddo
    allocate(dzw(nz+1,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
!          dzw(1,j,i) = zr(1,j,i)-zw(1,j,i) 
          dzw(1,j,i) = zr(1,j,i)-zw(0,j,i) !because zw indexed as croco from 0 to nz
          do k = 2,nz
             dzw(k,j,i) = zr(k,j,i)-zr(k-1,j,i) 
          enddo
!          dzw(nz+1,j,i) = zw(nz+1,j,i)-zr(nz,j,i) 
         dzw(nz+1,j,i) = zw(nz,j,i)-zr(nz,j,i) !because zw indexed as croco from 0 to nz
       enddo
    enddo
    !! Cell widths
    allocate(dxu(0:ny+1,nx+1))
    do i = 1,nx+1
       do j = 0,ny+1
          dxu(j,i) = 0.5_8*(dx(j,i)+dx(j,i-1))
       enddo
    enddo
    allocate(dyv(ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 1,ny+1
          dyv(j,i) = 0.5_8*(dy(j,i)+dy(j-1,i))
       enddo
    enddo
    !!  Areas
    allocate(Arx(nz,0:ny+1,nx+1))
    do i = 1,nx+1
       do j = 0,ny+1
          do k = 1,nz
             Arx(k,j,i) = 0.25_8*(dz(k,j,i)+dz(k,j,i-1))*(dy(j,i)+dy(j,i-1)) 
          enddo
       enddo
    enddo
    allocate(Ary(nz,ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 1,ny+1
          do k = 1,nz
             Ary(k,j,i) = 0.25_8*(dz(k,j,i)+dz(k,j-1,i))*(dx(j,i)+dx(j-1,i)) 
          enddo
       enddo
    enddo
    allocate(Arz(0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          Arz(j,i) = dx(j,i)*dy(j,i)
       enddo
    enddo
    !! Slopes in x- and y-direction defined at rho-points
    allocate(zx(nz,0:ny+1,0:nx+1))
    allocate(zy(nz,0:ny+1,0:nx+1))
    do i = 1,nx
       do j = 1,ny
          do k = 1,nz
             zy(k,j,i) = 0.5_8*(zr(k,j+1,i)-zr(k,j-1,i))/dy(j,i)
             zx(k,j,i) = 0.5_8*(zr(k,j,i+1)-zr(k,j,i-1))/dx(j,i)
          enddo
       enddo
    enddo
    call fill_halo(1,zy) ! copy interior value in the halo
    call fill_halo(1,zx) ! copy interior value in the halo
    allocate(zxdy(nz,0:ny+1,0:nx+1))
    allocate(zydx(nz,0:ny+1,0:nx+1))
    do k = 1,nz
       zydx(k,:,:) = zy(k,:,:)*dx(:,:)
       zxdy(k,:,:) = zx(k,:,:)*dy(:,:)
    enddo
    allocate(zyw(nz+1,0:ny+1,0:nx+1))
    allocate(zxw(nz+1,0:ny+1,0:nx+1))
    do i = 1,nx
       do j = 1,ny
          do k = 1,nz+1
!             zyw(k,j,i) = 0.5_8*(zw(k,j+1,i)-zw(k,j-1,i))/dy(j,i)
             zyw(k,j,i) = 0.5_8*(zw(k-1,j+1,i)-zw(k-1,j-1,i))/dy(j,i) !because zw indexed as croco from 0 to nz
!             zxw(k,j,i) = 0.5_8*(zw(k,j,i+1)-zw(k,j,i-1))/dx(j,i)
             zxw(k,j,i) = 0.5_8*(zw(k-1,j,i+1)-zw(k-1,j,i-1))/dx(j,i) !because zw indexed as croco from 0 to nz
          enddo
       enddo
    enddo
    call fill_halo(1,zyw) ! copy interior value in the halo
    call fill_halo(1,zxw) ! copy interior value in the halo
    allocate(cw(nz+1,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz+1
             cw(k,j,i) = Arz(j,i)/dzw(k,j,i) * (1._8 + zxw(k,j,i)*zxw(k,j,i)+zyw(k,j,i)*zyw(k,j,i))
          enddo
       enddo
    enddo

    !- bottom vertical momentum -!

    do i = 1,nx  
       do j = 1,ny

          k = 0 ! bottom
          w(k,j,i) = ( &
               + hlf * hlf * ( &
               + zxdy(k+1,j,i) * ( dx(j,i) + dx(j,i-1) ) * u(k+1,j,i) &
               + zxdy(k+1,j,i) * ( dx(j,i+1) + dx(j,i) ) * u(k+1,j,i+1) ) &
               + hlf * hlf * ( &
               + zydx(k+1,j,i) * ( dy(j,i) + dy(j-1,i) ) * v(k+1,j,i) &
               + zydx(k+1,j,i) * ( dy(j+1,i) + dy(j,i) ) * v(k+1,j+1,i) ) &
                      ) / (cw(k+1,j,i)*dzw(k+1,j,i))

       enddo
    enddo

    deallocate(dz)
    deallocate(dzw)
    deallocate(dxu)
    deallocate(dyv)
    deallocate(Arx)
    deallocate(Ary)
    deallocate(Arz)
    deallocate(zx)
    deallocate(zy)
    deallocate(zxdy)
    deallocate(zydx)
    deallocate(zxw)
    deallocate(zyw)
    deallocate(cw)

  end subroutine set_bbc
 
end module mg_set_bbc
