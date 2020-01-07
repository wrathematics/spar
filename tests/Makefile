all: core

make: core_make
run: core_run
miniclean: core_miniclean
clean: core_clean


core: core_miniclean core_make core_run

core_make: core_miniclean
	( cd core; make )

core_run: core_make
	( cd core; ./core )

core_miniclean:
	( cd core; make miniclean; )

core_clean:
	( cd core; make clean; )
