module mg_bbc

  use mg_cst
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
    real(kind=rp), dimension(:,:,:), pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:), pointer :: alpha

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx    => grid(1)%dx
    dy    => grid(1)%dy
    zxdy  => grid(1)%zxdy
    zydx  => grid(1)%zydx
    alpha => grid(1)%alpha

    !- bottom vertical momentum -!
    do i = 1,nx  
       do j = 1,ny

          k = 1 ! bottom

          w(k,j,i) = ( &
               - zxdy(k,j,i) / dy(j,i) * hlf * ( u(k,j,i) + u(k,j  ,i+1)) &
               - zydx(k,j,i) / dx(j,i) * hlf * ( v(k,j,i) + v(k,j+1,i  )) ) &
               / alpha(k,j,i)

       enddo
    enddo

  end subroutine set_bbc
 
end module mg_bbc
