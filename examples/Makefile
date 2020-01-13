all: core mpi

make: core_make mpi_make
clean: core_clean mpi_clean



core: core_clean core_make

core_make: core_clean
	( cd core; make -j )

core_clean:
	( cd core; make clean; )



mpi: mpi_clean mpi_make

mpi_make: mpi_clean
	( cd mpi; make -j )

mpi_clean:
	( cd mpi; make clean; )
