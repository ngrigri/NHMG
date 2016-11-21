program mg_testcuc

  use mg_mpi 
  use mg_tictoc
  use mg_zr_zw
  use mg_mpi_exchange_ijk
  use mg_setup_tests
  use nhmg

  implicit none

  integer(kind=4):: nxg        ! global x dimension
  integer(kind=4):: nyg        ! global y dimension
  integer(kind=4):: nzg        ! z dimension
  integer(kind=4):: npxg       ! number of processes in x
  integer(kind=4):: npyg       ! number of processes in y
  integer(kind=4):: nx, ny, nz ! local dimensions

  integer(kind=4):: ierr, np, rank
  integer(kind=ip) :: it, nit

  real(kind=8), dimension(:,:,:), pointer :: u,v,w

  real(kind=8) :: Lx, Ly, Htot
  real(kind=8) :: hc, theta_b, theta_s
  real(kind=8), dimension(:,:), pointer :: dx, dy, zeta, h
  real(kind=rp), dimension(:,:), pointer :: dxu, dyv

  real(kind=rp), dimension(:,:,:), pointer :: z_r
  real(kind=rp), dimension(:,:,:), pointer :: z_w

  call tic(1,'mg_testcuc')

  nit = 1

  !---------------!
  !- Ocean model -!
  !---------------!
  nxg  = 1024
  nyg  = 1024
  nzg  =   64

  npxg  = 2
  npyg  = 2

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call mpi_comm_size(mpi_comm_world, np, ierr)

  if (np /= (npxg*npyg)) then
     write(*,*) "Error: in number of processes !"
     stop -1
  endif

  nx = nxg / npxg
  ny = nyg / npyg
  nz = nzg

  !---------------------!
  !- Initialise nhmg -!
  !---------------------!
  if (rank == 0) write(*,*)'Initialise nhmg grids'

  call nhmg_init(nx,ny,nz,npxg,npyg)

  !---------------!
  !- Setup CUC   -!
  !---------------!
  if (rank == 0) write(*,*)'Initialise cuc bench'

  Lx   = 200d3
  Ly   = 200d3
  Htot = 4d3

  hc      = 250._rp
  theta_b =   6._rp
  theta_s =   6._rp

  if (myrank==0) then
     write(*,*)''
     write(*,*)'Lx, Ly, Htot:',Lx, Ly, Htot
     write(*,*)'hc, theta_b, theta_s:',hc, theta_b, theta_s
  endif

  !-------------------------------------------------!
  !- dx,dy,h and U,V,W initialisation (model vars) -!
  !-------------------------------------------------!
  allocate(   h(0:ny+1,0:nx+1))
  allocate(zeta(0:ny+1,0:nx+1))

  allocate(  dx(0:ny+1,0:nx+1))
  allocate(  dy(0:ny+1,0:nx+1))
  allocate(  dxu(0:nx+1,0:ny+1))
  allocate(  dyv(0:nx+1,0:ny+1))

  allocate(   z_r(0:nx+1,0:ny+1,1:nz))
  allocate(   z_w(0:nx+1,0:ny+1,1:nz+1))

  call setup_cuc(       &
       nx,ny,npxg,npyg, &
       dx,dy,           &
       dxu,dyv,         &
       zeta,h)

  !TODO -> calculate (or read) correct dxu, dyu
  call nhmg_set_horiz_grids(nx,ny,dx,dy,dxu,dyv)

  call setup_zr_zw(hc,theta_b,theta_s,zeta,h,z_r,z_w,'new_s_coord')
  call nhmg_set_vert_grids(nx,ny,nz,z_r,z_w)

  !-------------------------------------!
  !- U,V,W initialisation (model vars) -!
  !-------------------------------------!
  allocate(u(0:nx+1,0:ny+1,  nz))
  allocate(v(0:nx+1,0:ny+1,  nz))
  allocate(w(0:nx+1,0:ny+1,0:nz))

  do it = 1, nit

     u(:,:,:)      =  0._8
     v(:,:,:)      =  0._8

     w(:,:,0)      =  0._8
     w(:,:,1:nz-1) = -1._8
     w(:,:,nz)     =  0._8

     if (netcdf_output) then
        call write_netcdf(u,vname='u',netcdf_file_name='u.nc',rank=myrank,iter=it)
        call write_netcdf(v,vname='v',netcdf_file_name='v.nc',rank=myrank,iter=it)
        call write_netcdf(w,vname='w',netcdf_file_name='w.nc',rank=myrank,iter=it)
     endif

     !----------------------!
     !- Call nhmg solver -!
     !----------------------!
     if (rank == 0) write(*,*)'Call nhmg solver'
     call  nhmg_solve(nx,ny,nz,u,v,w)

     if (netcdf_output) then
        call write_netcdf(u,vname='uc',netcdf_file_name='uc.nc',rank=myrank,iter=it)
        call write_netcdf(v,vname='vc',netcdf_file_name='vc.nc',rank=myrank,iter=it)
        call write_netcdf(w,vname='wc',netcdf_file_name='wc.nc',rank=myrank,iter=it)
     endif

  enddo

  !---------------------!
  !- Deallocate memory -!
  !---------------------!
  if (rank == 0) write(*,*)'Cleaning memory before to finish the program.'
  call nhmg_clean()

  !------------------!
  !- End test-model -!
  !------------------!
  call mpi_finalize(ierr)

  call toc(1,'mg_testcuc')
  if(myrank == 0) call print_tictoc(myrank)

end program mg_testcuc
