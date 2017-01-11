program mg_testfillhalo

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_mpi_exchange
  use nhmg

  implicit none

  integer(kind=ip) :: nx, ny, nz  ! local dimensions
  integer(kind=ip) :: np, ierr, rank

  real(kind=rp), dimension(:,:)  , pointer :: a2D
  real(kind=rp), dimension(:,:,:), pointer :: a3D
  real(kind=rp), dimension(:,:,:), pointer :: a3Dp

  integer(kind=ip):: nxg  = 64       ! global x dimension
  integer(kind=ip):: nyg  = 64       ! global y dimension
  integer(kind=ip):: nzg  = 64       ! z dimension
  integer(kind=ip):: npxg  = 1       ! number of processes in x
  integer(kind=ip):: npyg  = 1       ! number of processes in y

  integer(kind=ip)  :: lun_nml = 4
  logical :: nml_exist=.false.

  namelist/fhparam/ &
       nxg        , &
       nyg        , &
       nzg        , &
       npxg       , &
       npyg 

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call mpi_comm_size(mpi_comm_world, np, ierr)

  !---------------------!
  !- Namelist (or not) -!
  !---------------------!
  !- Check if a fh_namelist file exist
  inquire(file='fh_namelist', exist=nml_exist)

  !- Read namelist file if it is present, else use default values
  if (nml_exist) then
     if (rank == 0) write(*,*)' Reading fh_namelist file'
     open(unit=lun_nml, File='fh_namelist', ACTION='READ')
     rewind(unit=lun_nml)
     read(unit=lun_nml, nml=fhparam)
  endif

  if (rank == 0) then
     write(*,*)'  test fill halo parameters:'
     write(*,*)'  - nxg     : ', nxg
     write(*,*)'  - nyg     : ', nyg
     write(*,*)'  - nzg     : ', nzg
     write(*,*)'  - npxg    : ', npxg 
     write(*,*)'  - npyg    : ', npyg
     write(*,*)'  '
  endif

  !---------------------!
  !- Global/local dim  -!
  !---------------------!
  if (np /= (npxg*npyg)) then
     write(*,*) "Error: in number of processes !", np, npxg, npyg
     stop -1
  endif

  nx = nxg / npxg
  ny = nyg / npyg
  nz = nzg

  !-------------------!
  !- Initialise nhmg -!
  !-------------------!
  call nhmg_init(nx,ny,nz,npxg,npyg)

  !---------------------!
  !- Setup fill halo   -!
  !---------------------!
  if (rank == 0) write(*,*)' Allocate fill halo test arrays'

  allocate( a2D (       0:ny+1,0:ny+1))
  allocate( a3D (1:nz  ,0:ny+1,0:nx+1))
  allocate( a3Dp(1:nz+1,0:ny+1,0:nx+1))

  if (rank == 0) write(*,*)' Initialize fill halo test arrays'

  a2D (:,:)   = zero
  a3D (:,:,:) = zero
  a3Dp(:,:,:) = zero
  a2D (  1:ny,1:nx) = rank
  a3D (:,1:ny,1:nx) = rank
  a3Dp(:,1:ny,1:nx) = rank

  call write_netcdf(a2D ,vname='a2D' ,netcdf_file_name='a2D.nc' ,rank=rank)
  call write_netcdf(a3D ,vname='a3D' ,netcdf_file_name='a3D.nc' ,rank=rank)
  call write_netcdf(a3Dp,vname='a3Dp',netcdf_file_name='a3Dp.nc',rank=rank)

  call fill_halo(1,a2D)
  call fill_halo(1,a3D)
  call fill_halo(1,a3Dp)
  if (rank == 0) write(*,*)' Writting arrays after fill_halo'
  call write_netcdf(a2D ,vname='a2D' ,netcdf_file_name='a2Dfh.nc' ,rank=rank)
  call write_netcdf(a3D ,vname='a3D' ,netcdf_file_name='a3Dfh.nc' ,rank=rank)
  call write_netcdf(a3Dp,vname='a3Dp',netcdf_file_name='a3Dpfh.nc',rank=rank)

  if (rank == 0) write(*,*)' Initialize fill halo test arrays (u)'

  a2D (:,:)   = zero
  a3D (:,:,:) = zero
  a2D (  1:ny,1:nx) = rank
  a3D (:,1:ny,1:nx) = rank

  call fill_halo(1,a2D,lbc_null='u')
  call fill_halo(1,a3D,lbc_null='u')
  if (rank == 0) write(*,*)' Writting arrays (u) after fill_halo'
  call write_netcdf(a2D ,vname='a2D' ,netcdf_file_name='a2Dufh.nc' ,rank=rank)
  call write_netcdf(a3D ,vname='a3D' ,netcdf_file_name='a3Dufh.nc' ,rank=rank)

  if (rank == 0) write(*,*)' Initialize fill halo test arrays (v)'

  a2D (:,:)   = zero
  a3D (:,:,:) = zero
  a2D (  1:ny,1:nx) = rank
  a3D (:,1:ny,1:nx) = rank

  call fill_halo(1,a2D,lbc_null='v')
  call fill_halo(1,a3D,lbc_null='v')

  if (rank == 0) write(*,*)' Writting arrays (v) after fill_halo'
  call write_netcdf(a2D ,vname='a2D' ,netcdf_file_name='a2Dvfh.nc' ,rank=rank)
  call write_netcdf(a3D ,vname='a3D' ,netcdf_file_name='a3Dvfh.nc' ,rank=rank)

  !---------------------!
  !- Deallocate memory -!
  !---------------------!
  if (rank == 0) write(*,*)' Clean memory before to finish the program.'
  call nhmg_clean()

  !----------------------!
  !- End Bench-fillhalo -!
  !----------------------!
  call mpi_finalize(ierr)

end program mg_testfillhalo

