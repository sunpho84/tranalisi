#include "PrecFloat.hpp"
#include "hpIntegrator.hpp"
#include <tranalisi.hpp>

using Real=double;
using Jack=jack_t<Real>;
using Jvec=vmeas_t<Jack>;
using Vect=Vector<Real,Dynamic>;
using Matr=Matrix<Real,Dynamic,Dynamic>;

const size_t T=26;

Real basis(const size_t& iT,
	     const Real& E)
{
  const Real t=2*iT+1;
  
  return exp(-t*E)+exp(-(T-t)*E);
}

const double& get(const double& x)
{
  return x;
}

double get(const PrecFloat& x)
{
  return x.get();
}

template <typename F>
Real integrate(F&& f,
	       const double& lower=0)
{
  if constexpr(std::is_same_v<Real,double>)
    return gslIntegrateUpToInfinity(f,lower);
  else
    return precIntegrateUpToInfinity(f,lower);
}

int main()
{
  PrecFloat::setDefaultPrecision(64);
  
  const size_t nT=(T/2+1)/2-1;

  const double sigma=1.0/4;
  //for(double sigma=1;sigma>=0.001;sigma/=2)
    {
      /// Target function
      const auto targ=
	[sigma](const Real& E) -> Real
	{
	  if(E==0)
	    return sqr(2/sigma*M_PI)*sigma;
	  else
	    return sqr(2/sigma*M_PI)*E/sinh(E/sigma);
	};
      
      Vect e(nT);
      Matr f(nT,nT);
      
      for(size_t iT=0;iT<nT;iT++)
	{
	  e(iT)=
	    integrate([targ,iT](const Real& E)
	    {
	      return targ(E)*basis(iT,E);
	    });
	  
	  for(size_t iS=iT;iS<nT;iS++)
	    f(iT,iS)=f(iS,iT)=
	      integrate([iT,iS](const Real& E)
	      {
		return basis(iT,E)*basis(iS,E);
	      });
	}
      
      Vect L(nT);
      Matr QSyst(nT,nT);
      
      /// Prepares the functional components
      for(size_t iT=0;iT<nT;iT++)
	{
	  L(iT)=e(iT);
	  for(size_t iS=0;iS<nT;iS++)
	    QSyst(iT,iS)=f(iT,iS);
	}
      
      /// Track the results
      Vect g;
      
      Matr Q=
	QSyst;
      
      g=Q.inverse()*L;
      
      // Plots the reconstruction of the target function
      grace_file_t recoPlot("/tmp/"+to_string(sigma)+".xmg");
      recoPlot.write_line([targ](const double& E){return get(targ(E));},1e-6,10*sigma);
      
      const auto reco=
	[&](const double& E)
      {
	Real R=0;
	for(size_t iT=0;iT<nT;iT++)
	  R+=basis(iT,E)*g[iT];
	return get(R);
      };
      
      recoPlot.write_line(reco,0,10*sigma);
      
      const Real err=
	integrate([&targ, &reco](const Real& E)
	{
	  return sqr(targ(E)-reco(get(E)));
	});
      
      cout<<"Err: "<<get(err)<<endl;
      
    }
  
  return 0;
}
