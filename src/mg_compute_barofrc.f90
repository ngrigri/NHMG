module mg_compute_barofrc

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine compute_barofrc(zr,zw,dt,ru,rv)

    real(kind=rp), dimension(:,:,:), pointer, intent(in)  :: zr,zw
    real(kind=rp),                            intent(in)  :: dt
    real(kind=rp), dimension(:,:)  , pointer, intent(out) :: ru,rv

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
!    real(kind=rp), dimension(:,:,:), pointer :: zw
    real(kind=rp), dimension(:,:)  , pointer :: dxu,dyv
    real(kind=rp), dimension(:,:,:), pointer :: dz
    real(kind=rp), dimension(:,:,:),   pointer :: zy,zx
    real(kind=rp), dimension(:,:,:),   pointer :: zydx,zxdy
    real(kind=rp), dimension(:,:,:), pointer :: p

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
!    zw => grid(1)%zw

    if (myrank==0) write(*,*)'- compute barofrc:'

    !! Cell heights
    allocate(dz(nz,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz
             dz(k,j,i) = zw(k+1,j,i)-zw(k,j,i)
          enddo
       enddo
    enddo
    !! Cell widths
    allocate(dxu(0:ny+1,nx+1))
    do i = 1,nx+1
       do j = 0,ny+1
          dxu(j,i) = hlf * (dx(j,i)+dx(j,i-1))
       enddo
    enddo
    allocate(dyv(ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 1,ny+1
          dyv(j,i) = hlf * (dy(j,i)+dy(j-1,i))
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
    call fill_halo(1,zy)
    call fill_halo(1,zx) 
    allocate(zxdy(nz,0:ny+1,0:nx+1))
    allocate(zydx(nz,0:ny+1,0:nx+1))
    do k = 1,nz
       zydx(k,:,:) = zy(k,:,:)*dx(:,:)
       zxdy(k,:,:) = zx(k,:,:)*dy(:,:)
    enddo

    !! Compute
    p => grid(1)%p

    ru(:,:) = 0._8
    rv(:,:) = 0._8

    do i = 1,nx+1
       do j = 0,ny+1 
!       do j = 1,ny

          k = 1

             ru(i,j) = ru(i,j) &

                     - (hlf*(dy(j,i)+dy(j,i-1)))*(hlf*(dz(k,j,i)+dz(k,j,i-1))) &
                       *(p(k,j,i)-p(k,j,i-1)) &

                     + qrt * (  & ! zero in our flat bottom case
                             + dx(j,i  )*zxdy(k,j,i  ) * (p(k+1,j,i)-p(k  ,j,i)) & 
                             +  & ! zero in our flat bottom case
                             + dx(j,i-1)*zxdy(k,j,i-1) * (p(k+1,j,i-1)-p(k  ,j,i-1)) &
                             )

          do k = 2,nz-1

             ru(i,j) = ru(i,j) &

                     - (hlf*(dy(j,i)+dy(j,i-1)))*(hlf*(dz(k,j,i)+dz(k,j,i-1))) &
                       *(p(k,j,i)-p(k,j,i-1)) &

                     + qrt * ( dx(j,i  )*zxdy(k,j,i  ) * (p(k  ,j,i)-p(k-1,j,i)) & 
                             + dx(j,i  )*zxdy(k,j,i  ) * (p(k+1,j,i)-p(k  ,j,i)) & 
                             + dx(j,i-1)*zxdy(k,j,i-1) * (p(k  ,j,i-1)-p(k-1,j,i-1)) & 
                             + dx(j,i-1)*zxdy(k,j,i-1) * (p(k+1,j,i-1)-p(k  ,j,i-1)) &
                             )

          enddo

          k = nz

             ru(i,j) = ru(i,j) &

                     - (hlf*(dy(j,i)+dy(j,i-1)))*(hlf*(dz(k,j,i)+dz(k,j,i-1))) &
                       *(p(k,j,i)-p(k,j,i-1)) &

                     + qrt * ( dx(j,i  )*zxdy(k,j,i  ) * (p(k  ,j,i)-p(k-1,j,i)) & 
                             + dx(j,i  )*zxdy(k,j,i  ) * two * (-p(k  ,j,i)) & 
                             + dx(j,i-1)*zxdy(k,j,i-1) * (p(k  ,j,i-1)-p(k-1,j,i-1)) & 
                             + dx(j,i-1)*zxdy(k,j,i-1) * two * (-p(k  ,j,i-1)) &
                             )

          ru(i,j) = ru(i,j)/dt

       enddo
    enddo
 
    do i = 0,nx+1
!    do i = 1,nx
       do j = 1,ny+1 

          k = 1

             rv(i,j) = rv(i,j) &

                     - (hlf*(dx(j,i)+dx(j-1,i)))*(hlf*(dz(k,j,i)+dz(k,j-1,i))) &
                       *(p(k,j,i)-p(k,j-1,i)) &

                     + qrt * (  & ! zero in our flat bottom case
                             + dy(j  ,i)*zydx(k,j  ,i) * (p(k+1,j,i)-p(k  ,j,i)) & 
                             +  & ! zero in our flat bottom case
                             + dy(j-1,i)*zydx(k,j-1,i) * (p(k+1,j-1,i)-p(k  ,j-1,i)) &
                             )

          do k = 2,nz-1

             rv(i,j) = rv(i,j) &

                     - (hlf*(dx(j,i)+dx(j-1,i)))*(hlf*(dz(k,j,i)+dz(k,j-1,i))) &
                       *(p(k,j,i)-p(k,j-1,i)) &

                     + qrt * ( dy(j  ,i)*zydx(k,j  ,i) * (p(k  ,j,i)-p(k-1,j,i)) & 
                             + dy(j  ,i)*zydx(k,j  ,i) * (p(k+1,j,i)-p(k  ,j,i)) & 
                             + dy(j-1,i)*zydx(k,j-1,i) * (p(k  ,j-1,i)-p(k-1,j-1,i)) & 
                             + dy(j-1,i)*zydx(k,j-1,i) * (p(k+1,j-1,i)-p(k  ,j-1,i)) &
                             )

          enddo

          k = nz

             rv(i,j) = rv(i,j) &

                     - (hlf*(dx(j,i)+dx(j-1,i)))*(hlf*(dz(k,j,i)+dz(k,j-1,i))) &
                       *(p(k,j,i)-p(k,j-1,i)) &

                     + qrt * ( dy(j  ,i)*zydx(k,j  ,i) * (p(k  ,j,i)-p(k-1,j,i)) & 
                             + dy(j  ,i)*zydx(k,j  ,i) * two * (-p(k  ,j,i)) & 
                             + dy(j-1,i)*zydx(k,j-1,i) * (p(k  ,j-1,i)-p(k-1,j-1,i)) & 
                             + dy(j-1,i)*zydx(k,j-1,i) * two * (-p(k  ,j-1,i)) &
                             )

          rv(i,j) = rv(i,j)/dt

       enddo
    enddo

    deallocate(dz)
    deallocate(dxu)
    deallocate(dyv)
    deallocate(zx)
    deallocate(zy)
    deallocate(zxdy)
    deallocate(zydx)

  end subroutine compute_barofrc

end module mg_compute_barofrc
