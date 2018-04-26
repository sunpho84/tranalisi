#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <gevp.hpp>

#include <Eigen/Dense>

using namespace Eigen;

vector<djvec_t> gevp(const vector<djvec_t> &d,const size_t t0)
{
  //Check correlators and get the time extent
  if(d.size()==0) CRASH("Needs at least 1 correlator");
  const size_t T=d[0].size();
  
  //Compute the number of eigenstates
  const size_t n=sqrt(d.size());
  if(n*n!=d.size()) CRASH("%zu correlators, not a square",d.size());
  
  //Work data
  typedef Matrix<double,Dynamic,Dynamic> Matr;
  GeneralizedEigenSolver<Matr> ges;
  
  //Output
  vector<djvec_t> l(n,djvec_t(T));
  
  for(size_t ijack=0;ijack<=njacks;ijack++)
    {
      //Fill the rhs matrix
      Matr b(n,n);
      for(size_t i=0;i<n;i++)
	for(size_t j=0;j<n;j++)
	  b(i,j)=d[i*n+j][t0][ijack];
      
      for(size_t t=0;t<T;t++)
	{
	  //Fill the lhs matrix
	  Matr a(n,n);
	  for(size_t i=0;i<n;i++)
	    for(size_t j=0;j<n;j++)
	      a(i,j)=d[i*n+j][t][ijack];
	  
	  //Compute and store
	  ges.compute(a,b);
	  for(size_t i=0;i<n;i++)
	    l[i][t][ijack]=ges.eigenvalues()(i).real();
	}
    }
  
  return l;
}
