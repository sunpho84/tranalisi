#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#ifdef USE_MPI
# include <mpi.h>
#endif

#include <GM3/phiFunction.hpp>
#include <GM3/params.hpp>

void init(int narg,char **arg)
{
  omp_set_num_threads(1);
  
  hashedPhiAndDerivCalculator.plot("plots/phi.xmg","plots/phiDeriv.xmg");

#ifdef USE_MPI
  MPI_Init(&narg,&arg);
  
  int temp;
  MPI_Comm_size(MPI_COMM_WORLD,&temp);
  nMPIranks=temp;
  MPI_Comm_rank(MPI_COMM_WORLD,&temp);
  MPIrank=temp;
#else
  nMPIranks=1;
  MPIrank=0;
#endif
  console.open((MPIrank==0)?"/dev/stdout":"/dev/null");
  
  set_njacks(30);
}

void close()
{
#ifdef USE_MPI
  MPI_Finalize();
#endif
}
