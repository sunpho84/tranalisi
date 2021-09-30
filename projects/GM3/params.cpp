#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <mpi.h>

#include <GM3/phiFunction.hpp>
#include <GM3/params.hpp>

void init(int narg,char **arg)
{
  omp_set_num_threads(1);
  
  hashedPhiAndDerivCalculator.plot("plots/phi.xmg","plots/phiDeriv.xmg");
  
  MPI_Init(&narg,&arg);
  
  int temp;
  MPI_Comm_size(MPI_COMM_WORLD,&temp);
  nMPIranks=temp;
  MPI_Comm_rank(MPI_COMM_WORLD,&temp);
  MPIrank=temp;
  
  console.open((MPIrank==0)?"/dev/stdout":"/dev/null");
  
  set_njacks(30);
}

void close()
{
  MPI_Finalize();
}
