program mg_testcuc

  use mg_mpi 
  use mg_tictoc
  use mg_setup_tests
  use nhydro

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
  real(kind=8), dimension(:,:), pointer :: rmask

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
  !- Initialise nhydro -!
  !---------------------!
  if (rank == 0) write(*,*)'Initialise nhydro grids'

  call nhydro_init(nx,ny,nz,npxg,npyg)

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
  allocate(  dx(0:ny+1,0:nx+1))
  allocate(zeta(0:ny+1,0:nx+1))
  allocate(  dy(0:ny+1,0:nx+1))

  call setup_cuc(nx,ny,nz,npxg,npyg,Lx,Ly,Htot,dx,dy,zeta,h)

  allocate(rmask(0:ny+1,0:nx+1))
  rmask(:,:) = 1._rp
  if (bmask) then
     !- Mask the boundaries
     call fill_halo_2D_bmask(1,rmask)   !- Generic n procs
  endif

  call nhydro_matrices(dx,dy,zeta,h,rmask,hc,theta_b,theta_s)

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
     !- Call nhydro solver -!
     !----------------------!
     if (rank == 0) write(*,*)'Call nhydro solver'
     call  nhydro_solve(nx,ny,nz,rmask,u,v,w)

     if (netcdf_output) then
        call write_netcdf(u,vname='uc',netcdf_file_name='uc.nc',rank=myrank,iter=it)
        call write_netcdf(v,vname='vc',netcdf_file_name='vc.nc',rank=myrank,iter=it)
        call write_netcdf(w,vname='wc',netcdf_file_name='wc.nc',rank=myrank,iter=it)
     endif

     !------------------------------------------------------------!
     !- Call nhydro correct to check if nh correction is correct -!
     !------------------------------------------------------------!
     if (rank == 0) write(*,*)'Check nondivergence'
     call nhydro_check_nondivergence(nx,ny,nz,rmask,u,v,w)

     if (netcdf_output) then
        call write_netcdf(grid(1)%b,vname='bc',netcdf_file_name='bc.nc',rank=myrank,iter=it)
     endif

  enddo

  !---------------------!
  !- Deallocate memory -!
  !---------------------!
  if (rank == 0) write(*,*)'Cleaning memory before to finish the program.'
  call nhydro_clean()

  !------------------!
  !- End test-model -!
  !------------------!
  call mpi_finalize(ierr)

  call toc(1,'mg_testcuc')
  if(myrank == 0) call print_tictoc(myrank)

end program mg_testcuc
