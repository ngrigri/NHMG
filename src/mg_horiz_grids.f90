module mg_horiz_grids

  use mg_mpi
  use mg_tictoc
  use mg_grids
  use mg_mpi_exchange
  use mg_gather
  use mg_namelist

  implicit none

contains

  !----------------------------------------
  subroutine fill_horiz_grids(dx,dy)

    real(kind=rp), dimension(:,:), intent(in) :: dx
    real(kind=rp), dimension(:,:), intent(in) :: dy

    integer(kind=ip)::  lev

    real(kind=rp), dimension(:,:), pointer :: dxf
    real(kind=rp), dimension(:,:), pointer :: dxc

    real(kind=rp), dimension(:,:), pointer :: dyf
    real(kind=rp), dimension(:,:), pointer :: dyc

    integer(kind=ip) :: ny,nx
    integer(kind=ip) :: nyf,nxf
    integer(kind=ip) :: nyc,nxc

    real(kind=rp), parameter :: hlf  = 0.5_rp

    if (myrank==0) write(*,*)'- define horizontal grids :'

    do lev = 1, nlevs

       if (myrank==0) write(*,*)'   lev=',lev

       nx=grid(lev)%nx
       ny=grid(lev)%ny

       if (lev == 1) then

          !NG: WARNING dx, dy and h model have to be defined 
          !NG: on (ny,nx) and not on (nx,ny) !!
          grid(lev)%dx(0:ny+1,0:nx+1) = dx
          grid(lev)%dy(0:ny+1,0:nx+1) = dy

       else

          nxf =grid(lev-1)%nx
          nyf =grid(lev-1)%ny

          dxf => grid(lev-1)%dx
          dyf => grid(lev-1)%dy

          if (grid(lev)%gather == 1) then
             nxc= nx/grid(lev)%ngx
             nyc= ny/grid(lev)%ngy
             allocate(dxc(0:nyc+1,0:nxc+1))
             allocate(dyc(0:nyc+1,0:nxc+1))
          else
             nxc = nx
             nyc = ny
             dxc => grid(lev)%dx
             dyc => grid(lev)%dy
          endif

          ! Coarsen dx, dy and h
          dxc(1:nyc,1:nxc) = hlf      * ( &
               dxf(1:nyf  :2,1:nxf  :2) + &
               dxf(2:nyf+1:2,1:nxf  :2) + &
               dxf(1:nyf  :2,2:nxf+1:2) + &
               dxf(2:nyf+1:2,2:nxf+1:2) )

          dyc(1:nyc,1:nxc) = hlf      * ( &
               dyf(1:nyf  :2,1:nxf  :2) + &
               dyf(2:nyf+1:2,1:nxf  :2) + &
               dyf(1:nyf  :2,2:nxf+1:2) + &
               dyf(2:nyf+1:2,2:nxf+1:2) )

          if (grid(lev)%gather == 1) then

             call gather(lev,dxc,grid(lev)%dx)
             call gather(lev,dyc,grid(lev)%dy)

             deallocate(dxc)
             deallocate(dyc)

          endif

       endif ! lev == 1

       call fill_halo(lev, grid(lev)%dx)
       call fill_halo(lev, grid(lev)%dy)

    enddo

  end subroutine fill_horiz_grids

end module mg_horiz_grids
