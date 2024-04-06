!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read terminal input and alter default parameters
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE TERMINAL
USE DEFINITION
USE CLA
IMPLICIT NONE

! String !
CHARACTER (len=99) :: freq, angles
CHARACTER (len=99) :: g_0, s_g, g_j
CHARACTER (len=99) :: p_j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! convert to string !
write(angles, *) theta_observe
write(freq, *) nu_light
write(g_0, *) gamma_0
write(g_j, *) gamma_j
write(s_g, *) s_gamma
write(p_j, *) phi_j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize CLA
CALL CLA_INIT

! Register? IDK !
CALL CLA_REGISTER('-o', '--outfile', 'str', cla_char, output_file)
CALL CLA_REGISTER('-s', '--snapshot', 'str', cla_char, grmhd_file)
CALL CLA_REGISTER('-f', '--frequency', 'real', cla_float, freq)
CALL CLA_REGISTER('-t', '--theta', 'real', cla_float, angles)
CALL CLA_REGISTER('-a', '--gamma_0', 'real', cla_float, g_0)
CALL CLA_REGISTER('-b', '--sgamma', 'real', cla_float, s_g)
CALL CLA_REGISTER('-c', '--gamma_j', 'real', cla_float, g_j)
CALL CLA_REGISTER('-d', '--phi_j', 'real', cla_float, p_j)

! Get values !
CALL CLA_GET_CHAR('--outfile', output_file)
CALL CLA_GET_CHAR('--snapshot', grmhd_file)
CALL CLA_GET_FLOAT_R8('--frequency', nu_light)
CALL CLA_GET_FLOAT_R8('--theta', theta_observe)
CALL CLA_GET_FLOAT_R8('--gamma_0', gamma_0)
CALL CLA_GET_FLOAT_R8('--gamma_j', gamma_j)
CALL CLA_GET_FLOAT_R8('--sgamma', s_gamma)
CALL CLA_GET_FLOAT_R8('--phi_j', phi_j)

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Open an existing dataset.
CALL h5dopen_f(file_id, "n1", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, nx, data_dims, error)

! Close the dataset.
CALL h5dclose_f(dset_id, error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Open an existing dataset.
CALL h5dopen_f(file_id, "n2", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, ny, data_dims, error)

! Close the dataset.
CALL h5dclose_f(dset_id, error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Open an existing dataset.
CALL h5dopen_f(file_id, "r", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, r_grid, data_dims, error)

! Close the dataset.
CALL h5dclose_f(dset_id, error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Open an existing dataset.
CALL h5dopen_f(file_id, "th", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, th_grid, data_dims, error)

! Close the dataset.
CALL h5dclose_f(dset_id, error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Open an existing dataset.
CALL h5dopen_f(file_id, "phi", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, phi_grid, data_dims, error)

! Close the dataset.
CALL h5dclose_f(dset_id, error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Open an existing dataset.
CALL h5dopen_f(file_id, "gamma", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, gamma, data_dims, error)

! Close the dataset.
CALL h5dclose_f(dset_id, error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Open an existing dataset.
CALL h5dopen_f(file_id, "nhat_r", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, nhat_r, data_dims, error)

! Close the dataset.
CALL h5dclose_f(dset_id, error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Open an existing dataset.
CALL h5dopen_f(file_id, "nhat_th", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, nhat_th, data_dims, error)

! Close the dataset.
CALL h5dclose_f(dset_id, error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Open an existing dataset.
CALL h5dopen_f(file_id, "nhat_phi", dset_id, error)

! Read the dataset.
CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, nhat_phi, data_dims, error)

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
integer(HSIZE_T) :: angle_disc_dims(3)
integer(HSIZE_T) :: data_dims(3), stokes_dims(2), angle_dims(1)
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
space_rank = 1
angle_dims(1) = n_angle

! open dataspace !
call h5screate_simple_f(space_rank,angle_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"gamma_j",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,gamma_j,angle_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
angle_dims(1) = n_angle

! open dataspace !
call h5screate_simple_f(space_rank,angle_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"theta_j",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,theta_j,angle_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
angle_dims(1) = n_angle

! open dataspace !
call h5screate_simple_f(space_rank,angle_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"gamma_0",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,gamma_0,angle_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
angle_dims(1) = n_angle

! open dataspace !
call h5screate_simple_f(space_rank,angle_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"s_gamma",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,s_gamma,angle_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
angle_dims(1) = n_angle

! open dataspace !
call h5screate_simple_f(space_rank,angle_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"s_power",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,s_power,angle_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
angle_dims(1) = n_angle

! open dataspace !
call h5screate_simple_f(space_rank,angle_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"c_power",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,c_power,angle_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
angle_dims(1) = n_angle

! open dataspace !
call h5screate_simple_f(space_rank,angle_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"nu_light",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,nu_light,angle_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
angle_dims(1) = n_angle

! open dataspace !
call h5screate_simple_f(space_rank,angle_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"bbody_temp",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bbody_temp,angle_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
angle_dims(1) = n_angle

! open dataspace !
call h5screate_simple_f(space_rank,angle_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"isctaui0",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,isctaui0,angle_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
angle_dims(1) = n_angle

! open dataspace !
call h5screate_simple_f(space_rank,angle_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"qsctaui0",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,qsctaui0,angle_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
angle_dims(1) = n_angle

! open dataspace !
call h5screate_simple_f(space_rank,angle_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"usctaui0",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,usctaui0,angle_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Choose according to mode !

! Output map !
IF(map) THEN

  ! define DIMENSION !
  space_rank = 3
  angle_disc_dims(1) = nx
  angle_disc_dims(2) = ny
  angle_disc_dims(3) = nz

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! open dataspace !
  call h5screate_simple_f(space_rank,angle_disc_dims,dspace_id,error)

  ! create dataset !
  call h5dcreate_f(file_id,"fsc",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

  ! write dataset !
  call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,fsc,angle_disc_dims,error)

  ! close dataset !
  call h5dclose_f(dset_id,error)

  ! close data space !
  call h5sclose_f(dspace_id,error)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! open dataspace !
  call h5screate_simple_f(space_rank,angle_disc_dims,dspace_id,error)

  ! create dataset !
  call h5dcreate_f(file_id,"gsc",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

  ! write dataset !
  call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,gsc,angle_disc_dims,error)

  ! close dataset !
  call h5dclose_f(dset_id,error)

  ! close data space !
  call h5sclose_f(dspace_id,error)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! open dataspace !
  call h5screate_simple_f(space_rank,angle_disc_dims,dspace_id,error)

  ! create dataset !
  call h5dcreate_f(file_id,"hsc",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

  ! write dataset !
  call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,hsc,angle_disc_dims,error)

  ! close dataset !
  call h5dclose_f(dset_id,error)

  ! close data space !
  call h5sclose_f(dspace_id,error)

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! close the file !
call h5fclose_f(file_id,error)

! close interface !
call h5close_f(error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE
