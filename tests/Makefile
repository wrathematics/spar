all: core mpi

make: core_make mpi_make
run: core_run mpi_run
miniclean: core_miniclean mpi_miniclean
clean: core_clean mpi_clean



core: core_miniclean core_make core_run

core_make: core_miniclean
	( cd core; make -j )

core_run: core_make
	( cd core; ./core )

core_miniclean:
	( cd core; make miniclean; )

core_clean:
	( cd core; make clean; )



mpi: mpi_miniclean mpi_make mpi_run

mpi_make: mpi_miniclean
	( cd mpi; make -j )

mpi_run: mpi_make
	( cd mpi; mpirun -np 2 ./mpi )

mpi_miniclean:
	( cd mpi; make miniclean; )

mpi_clean:
	( cd mpi; make clean; )
