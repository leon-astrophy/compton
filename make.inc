# Define F90 compiler #
F90=/home/leon/HDF5/bin/h5fc

# Choose between CPU or GPU threading #
#MULTICORE = -DCPU -fopenmp 

# Define F90 flags #
F90FLAGS= -static -cpp -Ofast -O4 -fstack-arrays -ffree-line-length-none $(MULTICORE)
LDFLAGS= -static -cpp -Ofast -O4 -fstack-arrays -ffree-line-length-none $(MULTICORE)
