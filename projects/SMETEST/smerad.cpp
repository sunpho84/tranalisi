#include <tranalisi.hpp>

//compute expected radius
double expected_rad(size_t n,double kappa,double plaq)
{
  kappa*=pow(plaq,0.25);
  return sqrt(n*kappa/(1.0+6.0*kappa));
}

//! integrate
djack_t integrate_corr(const djvec_t &c,const vector<double> &d,double powd=0.0)
{
  djack_t out;
  out=0.0;
  for(size_t iel=0;iel<c.size()-1;iel++)
    {
      double x1=d[iel+1];
      double x0=d[iel];
      double dx=x1-x0;
      djack_t y1=c[iel+1]*pow(x1,powd);
      djack_t y0=c[iel]*pow(x0,powd);
      
      out+=dx*(y1+y0)/2.0;
    }
  
  return 2.0*out;
}

int main(int narg,char **arg)
{
  raw_file_t obs("density","r");
  
  //set njacks
  size_t T=obs.read<size_t>("T");
  set_njacks(T);
  
  double kappa=obs.read<double>("Kappa");
  size_t nlevels=obs.read<int>("NLevels");
  size_t meas_each=obs.read<int>("MeasEach");
  obs.read<double>("TimePlaquette");
  double plaqs=obs.read<double>("SpatPlaquette");
  
  //keep distances
  vector<double> dists;
  
  for(size_t ilev=0;ilev<=nlevels;ilev+=meas_each)
    {
      double ex_rad=expected_rad(ilev,kappa,plaqs);
      
      //read the smearing level, to be ignored
      obs.read<size_t>("Smearlevel");
      
      //read ndists and resize
      size_t ndists=obs.read<size_t>("NDists");
      dists.resize(ndists);
      djvec_t dens(ndists);
      
      //read the whole distances
      for(size_t t=0;t<T;t++)
	{
	  //read over t, to be ignored
	  obs.read<size_t>("t");
	  
	  //read all distances and densities
	  for(size_t idist=0;idist<ndists;idist++)
	    {
	      size_t d2;
	      obs.read(d2);
	      dists[idist]=sqrt(d2);
	      obs.read(dens[idist][t]);
	    }
	}
      dens.clusterize();
      dens/=djack_t(dens[0]);
      
      //normalize the density
      dens/=integrate_corr(dens,dists);
      
      //write the measured distribution
      grace_file_t dens_plot("dens_"+to_string(ilev)+".xmg");
      dens_plot.write_vec_ave_err(dists,dens.ave_err());
      dens_plot.write_line([ex_rad](double x){return gauss_distr(x,0.0,ex_rad);},0.0,dists.back());
      
      //compute average squared radius
      djack_t r2=integrate_corr(dens,dists,2.0);
      djack_t r4=integrate_corr(dens,dists,4.0);
      djack_t r6=integrate_corr(dens,dists,6.0);
      djack_t r=sqrt(r2);
      djack_t m4=r4-3*pow(r,4.0);
      djack_t m6=r6-15*pow(r,6.0);
      cout<<"Radius at "<<ilev<<": "<<r.ave_err()<<" (expected: "<<ex_rad<<")"<<endl;
      cout<<"Deviation of fourth moment from gaussian at "<<ilev<<": "<<m4.ave_err()<<endl;
      cout<<"Deviation of sixth moment from gaussian at "<<ilev<<": "<<m6.ave_err()<<endl;
      
      //perform analysis
      djvec_t eff_rad=sqrt(-2.0*log(dens/dens[0]));
      for(size_t idist=1;idist<ndists;idist++)
	eff_rad[idist]=dists[idist]/eff_rad[idist];
      
      //plot
      grace_file_t plot("effrad_"+to_string(ilev)+".xmg");
      plot.write_vec_ave_err(dists,eff_rad.ave_err());
      plot.write_line([ex_rad](double x){return ex_rad;},0.0,dists.back());
    }
  
  return 0;
}
