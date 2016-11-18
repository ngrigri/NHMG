module mg_bbc

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine set_bbc(u,v,w)

    real(kind=rp), dimension(:,:,:), pointer, intent(in)     :: u,v
    real(kind=rp), dimension(:,:,:), pointer, intent(inout)  :: w

    integer(kind=ip):: k,j,i
    integer(kind=ip):: nx,ny,nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: zr,zw
    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:), pointer :: cw

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
    zr => grid(1)%zr
    zw => grid(1)%zw

    dzw => grid(1)%dzw
    cw  => grid(1)%cw
    zxdy => grid(1)%zxdy
    zydx => grid(1)%zydx

    !- bottom vertical momentum -!

    do i = 1,nx  
       do j = 1,ny

          k = 1 ! bottom
          w(k,j,i) = ( &
               + hlf * hlf * ( &
               + zxdy(k,j,i) * ( dx(j,i  ) + dx(j,i-1) ) * u(k,j,i  ) &
               + zxdy(k,j,i) * ( dx(j,i+1) + dx(j,i  ) ) * u(k,j,i+1) ) &
               + hlf * hlf * ( &
               + zydx(k,j,i) * ( dy(j  ,i) + dy(j-1,i) ) * v(k,j  ,i) &
               + zydx(k,j,i) * ( dy(j+1,i) + dy(j  ,i) ) * v(k,j+1,i) ) &
                      ) / (cw(k,j,i)*dzw(k,j,i))

       enddo
    enddo

  end subroutine set_bbc
 
end module mg_bbc
