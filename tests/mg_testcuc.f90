program mg_testcuc

  use mg_cst
  use mg_mpi 
  use mg_tictoc
  use mg_zr_zw
  use mg_mpi_exchange_ijk
  use mg_setup_tests
  use nhmg

  implicit none

  integer(kind=ip) :: it
  integer(kind=4):: ierr, np, rank
  integer(kind=4):: nx, ny, nz ! local dimensions

  real(kind=rp), dimension(:,:), pointer :: dx, dy 
  real(kind=rp), dimension(:,:), pointer :: zeta, h
  real(kind=rp), dimension(:,:), pointer :: dxu, dyv
  real(kind=rp), dimension(:,:,:), pointer :: z_r
  real(kind=rp), dimension(:,:,:), pointer :: z_w
  real(kind=rp), dimension(:,:,:), pointer :: u,v,w

  !- Namelist parameters -!
  integer(kind=ip) :: nit=1
  integer(kind=ip):: nxg  = 1024       ! global x dimension
  integer(kind=ip):: nyg  = 1024       ! global y dimension
  integer(kind=ip):: nzg  = 64         ! z dimension
  integer(kind=ip):: npxg  = 2         ! number of processes in x
  integer(kind=ip):: npyg  = 2         ! number of processes in y
  real(kind=rp) :: Lx      = 200.e3_rp ! Domain length in x [meter]
  real(kind=rp) :: Ly      = 200.e3_rp ! Domain length in y [meter]
  real(kind=rp) :: Htot    = 4.e3_rp   ! Depth [meter]
  real(kind=rp) :: hc      = 250._rp   !
  real(kind=rp) :: theta_b =   6._rp   !
  real(kind=rp) :: theta_s =   6._rp   !

  integer(kind=ip)  :: lun_nml = 4 ! Logical Unit Number
  logical :: nml_exist=.false.

  call tic(1,'mg_testcuc')

  namelist/cucparam/ &
       nit        , &
       nxg        , &
       nyg        , &
       nzg        , &
       npxg       , &
       npyg       , &
       Lx         , &
       Ly         , &
       Htot       , &
       hc         , &
       theta_b    , &
       theta_s

  !---------------------!
  !- Namelist (or not) -!
  !---------------------!
  !- Check if a cuc_namelist file exist
  inquire(file='cuc_namelist', exist=nml_exist)

  !- Read namelist file if it is present, else use default values
  if (nml_exist) then
     if (myrank == 0) write(*,*)'- Reading cuc_namelist file'
     open(unit=lun_nml, File='cuc_namelist', ACTION='READ')
     rewind(unit=lun_nml)
     read(unit=lun_nml, nml=cucparam)
  endif

  if (myrank == 0) then
     write(*,*)'test CUC parameters:'
     write(*,*)'  - nit     : ', nit 
     write(*,*)'  - nxg     : ', nxg
     write(*,*)'  - nyg     : ', nyg
     write(*,*)'  - nzg     : ', nzg
     write(*,*)'  - npxg    : ', npxg 
     write(*,*)'  - npyg    : ', npyg
     write(*,*)'  - Lx      : ', Lx
     write(*,*)'  - Ly      : ', Ly
     write(*,*)'  - Htot    : ', Htot
     write(*,*)'  - hc      : ', hc
     write(*,*)'  - theta_b : ', theta_b
     write(*,*)'  - theta_s : ', theta_s
     write(*,*)'  '
  endif

  !---------------!
  !- Ocean model -!
  !---------------!
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
