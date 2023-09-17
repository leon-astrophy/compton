!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read terminal input and alter default parameters
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE TERMINAL
USE DEFINITION
USE CLA
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize CLA
CALL CLA_INIT

! Register? IDK !
CALL CLA_REGISTER('-o', '--outfile', 'str', cla_char, output_file)
CALL CLA_REGISTER('-s', '--snapshot', 'str', cla_char, grmhd_file)

! Get values !
CALL CLA_GET_CHAR('--outfile', output_file)
CALL CLA_GET_CHAR('--snapshot', grmhd_file)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Open hdf5 space for data input
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HDF5_OPEN
USE DEFINITION
USE HDF5
IMPLICIT NONE

! File identifier
INTEGER(HID_T) :: file_id   

! Dataset identifier
INTEGER(HID_T) :: dset_id      

! Error flag
INTEGER :: error 

! Data dims !
INTEGER(HSIZE_T), DIMENSION(1) :: data_dims

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize FORTRAN interface !
CALL h5open_f(error)

! Open an existing file !
CALL h5fopen_f (grmhd_file, H5F_ACC_RDWR_F, file_id, error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Data length !
data_dims(1) = 1

! Open an existing dataset.
CALL h5dopen_f(file_id, "n1", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, nx, data_dims, error)

! Close the dataset.
CALL h5dclose_f(dset_id, error)

! Open an existing dataset.
CALL h5dopen_f(file_id, "n2", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, ny, data_dims, error)

! Close the dataset.
CALL h5dclose_f(dset_id, error)

! Open an existing dataset.
CALL h5dopen_f(file_id, "n3", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, nz, data_dims, error)

! Close the dataset.
CALL h5dclose_f(dset_id, error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
! Close the file.
CALL h5fclose_f(file_id, error)

! Close FORTRAN interface.
CALL h5close_f(error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Load hdf5 arrays into xraypol
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HDF5_LOAD
USE DEFINITION
USE HDF5
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! File identifier
INTEGER(HID_T) :: file_id   

! Dataset identifier
INTEGER(HID_T) :: dset_id, dspace_id    

! Error flag
INTEGER :: error 

! Data dims !
INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

integer(hsize_t), dimension(3)               :: dims,maxdims

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize FORTRAN interface !
CALL h5open_f(error)

! Open an existing file !
CALL h5fopen_f (grmhd_file, H5F_ACC_RDWR_F, file_id, error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Data length !
data_dims(1) = nx
data_dims(2) = ny
data_dims(3) = nz

! Open an existing dataset.
CALL h5dopen_f(file_id, "r", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, r_grid, data_dims, error)

! Close the dataset.
CALL h5dclose_f(dset_id, error)

! Open an existing dataset.
CALL h5dopen_f(file_id, "th", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, th_grid, data_dims, error)

! Close the dataset.
CALL h5dclose_f(dset_id, error)

! Open an existing dataset.
CALL h5dopen_f(file_id, "phi", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, phi_grid, data_dims, error)

! Close the dataset.
CALL h5dclose_f(dset_id, error)

! Open an existing dataset.
CALL h5dopen_f(file_id, "vol", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vol, data_dims, error)

! Close the dataset.
CALL h5dclose_f(dset_id, error)

! Open an existing dataset.
CALL h5dopen_f(file_id, "gamma", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, gamma, data_dims, error)

! Close the dataset.
CALL h5dclose_f(dset_id, error)

! Open an existing dataset.
CALL h5dopen_f(file_id, "nebar", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, nebar, data_dims, error)

! Close the dataset.
CALL h5dclose_f(dset_id, error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
! Close the file.
CALL h5fclose_f(file_id, error)

! Close FORTRAN interface.
CALL h5close_f(error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Open hdf5 space for data output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HDF5_OUT
USE DEFINITION
USE HDF5
IMPLICIT NONE

! for HDF5 !
character(len=99) :: filename
integer :: error, space_rank
integer(HSIZE_T) :: angular_dims(5), data_dims(3), stokes_dims(2), angle_dims(1)
integer(HID_T) :: file_id, dspace_id, dset_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! assign !
filename = trim(adjustl(output_file)) //'.hdf5'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! create interface !
call h5open_f(error)

! open the file !
call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
angle_dims(1) = n_angle

! open dataspace !
call h5screate_simple_f(space_rank,angle_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"theta0",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,theta_0,angle_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 2
stokes_dims(1) = 3
stokes_dims(2) = n_angle

! open dataspace !
call h5screate_simple_f(space_rank,stokes_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"stokes",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,stokes_total,stokes_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
data_dims(1) = nx
data_dims(2) = ny
data_dims(3) = nz

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"r-dir",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,r_grid,data_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
data_dims(1) = nx
data_dims(2) = ny
data_dims(3) = nz

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"th-dir",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,th_grid,data_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
data_dims(1) = nx
data_dims(2) = ny
data_dims(3) = nz

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"phi-dir",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi_grid,data_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
data_dims(1) = nx
data_dims(2) = ny
data_dims(3) = nz

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"electron",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,nebar,data_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
data_dims(1) = nx
data_dims(2) = ny
data_dims(3) = nz

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"volume",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,vol,data_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
data_dims(1) = nx
data_dims(2) = ny
data_dims(3) = nz

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"gammam1",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,gammam1,data_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 5
angular_dims(1) = 2
angular_dims(2) = 3
angular_dims(3) = nx
angular_dims(4) = ny
angular_dims(5) = nz

! open dataspace !
call h5screate_simple_f(space_rank,angular_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"angular",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,ang_fac_out,angular_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! close the file !
call h5fclose_f(file_id,error)

! close interface !
call h5close_f(error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE
