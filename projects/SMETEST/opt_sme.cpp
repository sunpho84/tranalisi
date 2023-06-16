#include <tranalisi.hpp>

const int nSme=10;
const int sme[nSme]={0,5,10,20,30,40,80,120,160,240};
const double rho=0.4;
const int nMes=6;
const string mes[nMes]={"P","K","ETAS","D","Ds","ETAC"};
const size_t T=64;
const size_t tStar=5;

int main()
{
  const index_t id({{"",nMes},{"",nSme},{"",nSme}});
  set_njacks(20);
  
  raw_file_t infile("data/data.dat","r");
  vector<djvec_t> data(id.max(),djvec_t(T));
  infile.bin_read(data);
  
  grace_file_t resRad("plots/resRad.xmg");
  grace_file_t resSup("plots/resSup.xmg");
  for(size_t iMes=0;iMes<nMes;iMes++)
    {
      vector<pair<double,djvec_t>> mEff;
      for(size_t iSm1=0;iSm1<nSme;iSm1++)
	for(size_t iSm2=iSm1;iSm2<nSme;iSm2++)
	  {
	    const djvec_t tmp=data[id({iMes,iSm1,iSm2})].symmetrized();
	    
	    const double smeRad=sqrt(sme[iSm1]+sme[iSm2]);
	    
	    mEff.emplace_back(smeRad,effective_mass(tmp));
	  }
      
      sort(mEff.begin(),mEff.end());
      size_t iMin=0;
      double min=1e300;
      for(size_t i=0;i<mEff.size();i++)
	{
	  const double mt=mEff[i].second[tStar].ave();
	  if(mt<min)
	    {
	      min=mt;
	      iMin=i;
	    }
	}
      
      const djack_t mMin=constant_fit(mEff[iMin].second,14,20,"plots/mEff"+mes[iMes]+".xmg");
      
      vector<pair<double,djack_t>> estF;
      for(const auto& m : mEff)
	{
	  const djack_t est=(m.second[tStar]-mMin)/(mEff[0].second[tStar]-mMin);
	  estF.emplace_back(m.first,est);
	}
      
      vector<bool> discr;
      const ave_err_t y0=estF[iMin].second.ave_err();
      for(const auto& est : estF)
	{
	  const ave_err_t y=est.second.ave_err();
	  discr.push_back((y.ave()-y0.ave())/sqrt(sqr(y.err())+sqr(y0.err()))<3);
	}
      
      jack_fit_t fitter;
      size_t nDiscr=0;
      for(size_t i=0;i<estF.size();i++)
	if(discr[i])
	  nDiscr++;
      
      vector<double> x(nDiscr);
      djvec_t y(nDiscr);
      size_t iDiscr=0;
      for(size_t i=0;i<estF.size();i++)
	if(discr[i])
	  {
	    x[iDiscr]=log(estF[i].first);
	    y[iDiscr]=estF[i].second;
	    iDiscr++;
	  }
      
      const djvec_t p=poly_fit(x,y,2);
      
      grace_file_t plot("plots/est"+mes[iMes]+".xmg");
      plot.set_no_line();
      for(const auto& [cond,color] : {make_pair(false,grace::RED),{true,grace::GREEN4}})
	{
	  plot.set_all_colors(color);
	  for(size_t i=0;i<estF.size();i++)
	    if(discr[i]==cond)
	      {
		const ave_err_t y=estF[i].second.ave_err();
		plot.write_ave_err(estF[i].first,y);
	      }
	  plot.new_data_set();
	}
      
      plot.write_polygon([p](const double& x)
      {
	return poly_eval(p,log(x));
      },exp(x.front()),exp(x.back()),grace::BLUE);
      
      const djack_t x0=exp(-p[1]/(2*p[2]));
      const djack_t yVal=poly_eval(p,log(x0));
      plot.write_ave_err(x0.ave(),yVal.ave_err());
      
      resRad.write_ave_err(mMin.ave_err(),(1/x0).ave_err());
      resSup.write_ave_err(mMin.ave_err(),yVal.ave_err());
    }
  
  return 0;
}
