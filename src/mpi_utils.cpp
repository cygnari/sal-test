#include <cassert>
#include <mpi.h>
#include <iostream>
#define assertm(exp, msg) assert(((void)msg, exp))

bool test_is_same(const int x, MPI_Comm mpi_communicator = MPI_COMM_WORLD) { // test if all processes have the same value for a variable
  int p[2];
  p[0] = -x;
  p[1] = x;
  MPI_Allreduce(MPI_IN_PLACE, p, 2, MPI_INT, MPI_MIN, mpi_communicator);
  return (p[0] == -p[1]);
}
