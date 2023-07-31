#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#include <iostream>
#include <random>
#include <vector>

#ifndef EXTERN_RANDOM
 #define EXTERN_RANDOM extern
 #define INIT_TO(A)
#else
 #define INIT_TO(A) =A
#endif

#include <ave_err.hpp>

using namespace std;

//! random generator
template <typename BaseGen=mt19937_64>
struct gen_t:
  public BaseGen
{
  //! init with a seed
  gen_t(const int& seed) :
    BaseGen(seed)
  {
  }
  
  //! return a real random number in the range [min,max)
  double get_double(const double& min,
		    const double& max)
  {
    return uniform_real_distribution<double>(min,max)(*this);
  }
  
  //! return an integer random number in the range [min,max)
  double get_int(const int& min,
		 const int& max)
  {
    return uniform_int_distribution<>(min,max-1)(*this);
  }
  
  //! return an integer random number in the range [min,max)
  double get_gauss(const double& ave,
		   const double& sig)
  {
    return normal_distribution<>(ave,sig)(*this);
  }
  
private:
  
  //! init without a seed
  gen_t() :
    BaseGen()
  {
  }
};

///////////////////////////////////////////////////////////// gauss_filler_t /////////////////////////////////////////////////////

//! allows to fill from gauss
class gauss_filler_t : pair<ave_err_t,int>
{
public:
  //! fill from ave_err and seed
  gauss_filler_t(const ave_err_t &ext_ae,int ext_seed)
  {
    ae=ext_ae;
    seed=ext_seed;
  }
  
  //! fill from ave, err and seed
  gauss_filler_t(double ave,double err,int ext_seed) : gauss_filler_t(ave_err_t(ave,err),ext_seed) {}
  
  //! rebind ave_err
  ave_err_t &ae=first;
  
  //! rebind seed
  int &seed=second;
};

template <typename G>
inline vector<double> multivariate(const vector<double> &_mean,const vector<double>& _sigma,gen_t<G>& gen)
{
  const int n=_mean.size();
  
  Eigen::VectorXd mean(n);
  for(int i=0;i<n;i++) mean(i)=_mean[i];
  Eigen::MatrixXd sigma(n,n);
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      sigma(i,j)=_sigma[j+n*i];
  
  Eigen::VectorXd x(n);
  for(int i=0;i<n;i++)
    x(i)=gen.get_gauss(0,1);
  
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(sigma);
  const Eigen::MatrixXd eigenvectors=eigen_solver.eigenvectors().real();
  
  const Eigen::MatrixXd eigenvalues=eigen_solver.eigenvalues().real().asDiagonal();
  
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(eigenvalues);
  const Eigen::MatrixXd sqrt_eigenvalues=es.operatorSqrt();
  const Eigen::MatrixXd Q=eigenvectors*sqrt_eigenvalues;
  
  const Eigen::VectorXd temp=Q*x+mean;
  
  vector<double> out(n);
  for(int i=0;i<n;i++) out[i]=temp(i);
  
  return out;
}

//! if we ever needed a global generator...
//EXTERN_RANDOM gen_t glb_gen;

#undef EXTERN_RANDOM
#undef INIT_TO

#endif
