#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <cstdio>

int size;
int rank;


int main(int argc, char *argv[])
{
  int num_failed_tests;
  MPI_Init(NULL, NULL);
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0)
    fclose(stdout);
  
  num_failed_tests = Catch::Session().run(argc, argv);
  
  MPI_Finalize();
  
  return num_failed_tests;
}
