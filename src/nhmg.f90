module nhmg

  use mg_cst
  use mg_mpi
  use mg_grids
  use mg_namelist
  use mg_tictoc
  use mg_mpi_exchange
  use mg_horiz_grids
  use mg_vert_grids
  use mg_projection
  use mg_solvers
  use mg_netcdf_out

  implicit none

contains

  !--------------------------------------------------------------
  subroutine nhmg_init(nx,ny,nz,npxg,npyg)
      
    integer(kind=ip), intent(in) :: nx, ny, nz
    integer(kind=ip), intent(in) :: npxg, npyg

    call mg_mpi_init()

    if (myrank==0) write(*,*)' nhmg_init:'

    call read_nhnamelist(vbrank=myrank)

    call define_grids(npxg,npyg,nx,ny,nz)

    call define_neighbours()

    call print_grids()

  end subroutine nhmg_init

  !--------------------------------------------------------------
  subroutine nhmg_comp_rw(nx,ny,nz,rua,rva,rwa)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz  ), target, intent(in) :: rua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz  ), target, intent(in) :: rva
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(out):: rwa

    real(kind=rp), dimension(:,:,:), pointer :: ru,rv,rw

!!! dirty reshape arrays indexing ijk -> kji !!!
    real(kind=rp), dimension(:,:,:), allocatable, target :: rub,rvb
!!! dirty reshape arrays indexing ijk -> kji !!!

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: dz
    real(kind=rp), dimension(:,:,:), pointer :: zxdy,zydx

    real(kind=rp),dimension(nz) :: dzdhp
    integer(kind=ip) :: i,j,k

    dx    => grid(1)%dx
    dy    => grid(1)%dy
    dz    => grid(1)%dz
    zxdy  => grid(1)%zxdy
    zydx  => grid(1)%zydx

!!! dirty reshape arrays indexing ijk -> kji !!!
    allocate(rub(1:nz,0:ny+1,0:nx+1))
    allocate(rvb(1:nz,0:ny+1,0:nx+1))
    do i = 1,nx+1
      do j = 0,ny+1
        do k = 1,nz
          rub(k,j,i) = rua(i,j,k)
        enddo
      enddo
    enddo
    do i = 0,nx+1
      do j = 1,ny+1
        do k = 1,nz
           rvb(k,j,i) = rva(i,j,k)
        enddo
      enddo
    enddo
    ru => rub
    rv => rvb
    call fill_halo(1,ru)
    call fill_halo(1,rv)
!!! 
    rw => rwa

    !XXXXXXXXXXXXXXXXXXXXXXXXX
    !care with dz. properly updated?

    do i = 1,nx
       do j = 1,ny
          do k = 1,nz
             dzdhp(k) = zxdy(k,j,i)/dy(j,i)*(               &
                  ru(k,j,i  )/(dz(k,j,i) + dz(k,j,i-1))     &
                  + ru(k,j,i+1)/(dz(k,j,i) + dz(k,j,i+1)) ) &
                  + zydx(k,j,i)/dx(j,i)*(                   &
                  rv(k,j  ,i)/(dz(k,j,i) + dz(k,j-1,i))     &
                  + rv(k,j+1,i)/(dz(k,j,i) + dz(k,j+1,i)) )
          enddo
          do k = 2,nz
             rw(i,j,k) =  - 0.5*( dzdhp(k) + dzdhp(k-1) ) &
                          * 0.5*( dz(k,j,i) + dz(k-1,j,i) )
          enddo
          rw(i,j,nz+1) =  - dzdhp(nz) *0.5*dz(nz,j,i)
       enddo
    enddo

!!! dirty reshape arrays indexing kji -> ijk !!!
    deallocate(rub)
    deallocate(rvb)
!!! dirty reshape arrays indexing kji -> ijk !!!

  end subroutine nhmg_comp_rw

  !--------------------------------------------------------------
  subroutine nhmg_matrices(nx,ny,nz,zra,Hza,dxa,dya)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(0:nx+1,0:ny+1,1:nz)          , intent(in) :: zra
    real(kind=rp), dimension(0:nx+1,0:ny+1,1:nz)          , intent(in) :: Hza
    real(kind=rp), dimension(0:nx+1,0:ny+1)     , optional, intent(in) :: dxa
    real(kind=rp), dimension(0:nx+1,0:ny+1)     , optional, intent(in) :: dya

    real(kind=rp), dimension(:,:,:), pointer :: zr,Hz
    real(kind=rp), dimension(:,:)  , pointer :: dx,dy

!!! dirty reshape arrays indexing ijk -> kji !!!
    integer(kind=ip) :: i,j,k
    real(kind=rp), dimension(1:nz,0:ny+1,0:nx+1), target :: zrb
    real(kind=rp), dimension(1:nz,0:ny+1,0:nx+1), target :: Hzb 
    real(kind=rp), dimension(0:ny+1,0:nx+1),      target :: dxb,dyb
!!!

    integer(kind=ip), save :: iter_matrices=0
    iter_matrices = iter_matrices + 1

!    if (myrank==0) write(*,*)' nhmg_matrices: ',iter_matrices

    !--------------------!
    !- Horizontal grids -!
    !--------------------!
    if (present(dxa) .and. present(dya)) then

       dxb = transpose(dxa)
       dyb = transpose(dya)
       dx => dxb
       dy => dyb

       call set_horiz_grids(dx,dy)

!!$       if (check_output) then
!!$          call write_netcdf(grid(1)%dx,vname='dx',netcdf_file_name='mat.nc',rank=myrank,iter=iter_matrices)
!!$          call write_netcdf(grid(1)%dy,vname='dy',netcdf_file_name='mat.nc',rank=myrank,iter=iter_matrices)
!!$       endif

       if (associated(dx)) dx => null()
       if (associated(dy)) dy => null()

    end if

    !------------------!
    !- Vertical grids -!
    !------------------!
!!! dirty reshape arrays indexing ijk -> kji !!!
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz
             zrb(k,j,i) = zra(i,j,k)
             Hzb(k,j,i) = Hza(i,j,k)
          enddo
       enddo
    enddo
    zr => zrb
    Hz => Hzb
!!!

    call set_vert_grids(zr,Hz)

!!$    if (check_output) then
!!$       !if ((iter_matrices .EQ. 1) .OR. (iter_matrices .EQ. 2)) then
!!$          call write_netcdf(grid(1)%zr,vname='zr',netcdf_file_name='mat.nc',rank=myrank,iter=iter_matrices)
!!$          call write_netcdf(grid(1)%dz,vname='dz',netcdf_file_name='mat.nc',rank=myrank,iter=iter_matrices)
!!$       !endif
!!$    endif

    if (associated(zr)) zr => null()
    if (associated(Hz)) Hz => null()

    !------------!
    !- matrices -!
    !------------!
    call set_matrices()

!!$    if (check_output) then
!!$       !if ((iter_matrices .EQ. 1) .OR. (iter_matrices .EQ. 2)) then
!!$          call write_netcdf(grid(1)%cA,vname='cA',netcdf_file_name='mat.nc',rank=myrank,iter=iter_matrices)
!!$       !endif
!!$    endif

  end subroutine nhmg_matrices

  !--------------------------------------------------------------
  subroutine nhmg_solve(ua,va,wa,Hza,fill_hz)

    real(kind=rp), dimension(:,:,:), intent(in) :: ua
    real(kind=rp), dimension(:,:,:), intent(in) :: va
    real(kind=rp), dimension(:,:,:), intent(in) :: wa    
    real(kind=rp), dimension(:,:,:), intent(in) :: Hza
    logical :: fill_hz

    real(kind=rp), dimension(:,:),   pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: u,v,w,dz

    integer(kind=ip) :: i,j,k,is,js,ishift
    integer(kind=ip) :: nx,ny,nz

    integer(kind=ip), save :: iter_solve=0
    iter_solve = iter_solve + 1

    !    if (myrank==0) write(*,*)' nhmg_solve:',iter_solve

    call tic(1,'nhmg_solve')

    !    write(*,*) 'rank',myrank,'lbound(ua)',lbound(ua)
    !    write(*,*) 'rank',myrank,'ubound(ua)',ubound(ua)
    !    write(*,*) 'rank',myrank,'shape(ua)',shape(ua)

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx => grid(1)%dx
    dy => grid(1)%dy
    u  => grid(1)%u
    v  => grid(1)%v
    w  => grid(1)%w

    ! need to update dz because define_matrices may not be called every time step
    dz => grid(1)%dz

    ishift=2

    do k=1,nz
       do j=0,ny+1
          js=j+ishift
          do i=0,nx+1
             is=i+ishift
             dz(k,j,i) = Hza(is,js,k)
          enddo
       enddo
    enddo

    if (fill_hz) then
       call fill_halo(1,dz)
    endif

    ! set fluxes
    do k=1,nz
       do j=1,ny
          js=j+ishift
          do i=1,nx+1
             is=i+ishift
             u(k,j,i) = ua(is,js,k) * &
                  qrt * (dz(k,j,i) + dz(k,j,i-1)) * (dy(j,i)+dy(j,i-1))
          enddo
       enddo
       do j=1,ny+1
          js=j+ishift
          do i=1,nx
             is=i+ishift
             v(k,j,i) = va(is,js,k) * &
                  qrt * (dz(k,j,i) + dz(k,j-1,i)) * (dx(j,i)+dx(j-1,i))
          enddo
       enddo
       do j=1,ny
          js=j+ishift
          do i=1,nx
             is=i+ishift
             w(k+1,j,i) = wa(is,js,k+1) * &
                  dx(j,i) * dy(j,i)
          enddo
       enddo
    enddo
    w(1,:,:) = zero

    !- set rhs and solve for p
    call set_rhs()
    call solve_p()

    if (check_output) then
       !if ((iter_solve .EQ. 1) .OR. (iter_solve .EQ. 2)) then
       call write_netcdf(grid(1)%b,vname='b',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
       call write_netcdf(grid(1)%p,vname='p',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
       call write_netcdf(grid(1)%r,vname='r',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
       !endif
    endif

    !- correct
    call correct_uvw()

    if (check_output) then
       !if ((iter_solve .EQ. 1) .OR. (iter_solve .EQ. 2)) then
       call write_netcdf(u,vname='uout',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
       call write_netcdf(v,vname='vout',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
       call write_netcdf(w,vname='wout',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
       !endif
    endif

    !- check step - non-divergence of the projected u,v,w
    if (check_output) then
       call set_rhs()
       !if ((iter_solve .EQ. 1) .OR. (iter_solve .EQ. 2)) then
       call write_netcdf(grid(1)%b,vname='bout',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
       !endif
    endif

    !- check step - the projected u,v,w do not work

    !    if (associated(u)) u => null()
    !    if (associated(v)) v => null()
    !    if (associated(w)) w => null()

    call toc(1,'nhmg_solve')

  end subroutine nhmg_solve

  !--------------------------------------------------------------
  subroutine nhmg_clean()

    real(kind=rp) :: tstart,tend,perf

    call grids_dealloc()

    call print_tictoc()

  end subroutine nhmg_clean

end module nhmg

