# Define F90 compiler #
F90=h5fc

# Define F90 flags #
F90FLAGS= -cpp -Ofast -Mipa=fast -O4 -Mvect=simd -Mvect=levels:5 -Msmartalloc -Mstack_arrays -Mprefetch -fopenmp -ffree-line-length-none
LDFLAGS= -cpp -Ofast -Mipa=fast -O4 -Mvect=simd -Mvect=levels:5 -Msmartalloc -Mstack_arrays -Mprefetch -fopenmp -ffree-line-length-none
