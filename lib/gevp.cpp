#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <gevp.hpp>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;

/// Original eigenvectors are those needed to build the optimal operator
tuple<vector<djvec_t>,vector<djvec_t>,vector<djvec_t>> gevp(const vector<djvec_t> &d,const size_t t0)
{
  //Check correlators and get the time extent
  if(d.size()==0) CRASH("Needs at least 1 correlator");
  const size_t T=d[0].size();
  
  //Compute the number of eigenstates
  const size_t n=sqrt(d.size());
  if(n*n!=d.size()) CRASH("%zu correlators, not a square",d.size());
  
  //Work data
  typedef Matrix<double,Dynamic,Dynamic> Matr;
  EigenSolver<Matr> ge;
  // GeneralizedEigenSolver<Matr> ges;
  
  //Output
  vector<djvec_t> l(n,djvec_t(T));
  vector<djvec_t> recastEigvec(n*n,djvec_t(T));
  vector<djvec_t> origEigvec(n*n,djvec_t(T));
  
  for(size_t ijack=0;ijack<=njacks;ijack++)
    {
      //Fill the rhs matrix
      Matr b(n,n);
      for(size_t i=0;i<n;i++)
	for(size_t j=0;j<n;j++)
	  b(i,j)=d[i*n+j][t0][ijack];
      
      const Matr binv=b.inverse();
      
      for(size_t t=0;t<T;t++)
	{
	  //Fill the lhs matrix
	  Matr a(n,n);
	  for(size_t i=0;i<n;i++)
	    for(size_t j=0;j<n;j++)
	      a(i,j)=d[i*n+j][t][ijack];
	  
	  //Recast problem
	  const Matr c=binv*a;
	  
	  //Compute and store
	  ge.compute(c);
	  for(size_t i=0;i<n;i++)
	    l[i][t][ijack]=ge.eigenvalues()(i).real();
	  
	  auto vEig=(b.inverse().sqrt()*ge.eigenvectors()).eval();
	  
	  for(size_t i=0;i<n;i++)
	    for(size_t j=0;j<n;j++)
	      {
		recastEigvec[i*n+j][t][ijack]=ge.eigenvectors()(i,j).real();
		origEigvec[i*n+j][t][ijack]=vEig(i,j).real();
	      }
	}
    }
  
  return {l,recastEigvec,origEigvec};
}
