#include <set>
#include <tranalisi.hpp>

const int nSme=10;
const int sme[nSme]={0,5,10,20,30,40,80,120,160,240};
const double rho=0.4;
const int nMes=10;
const string mes[nMes]={"P","K","ETAS","D","Ds","ETAC","Hl","Hs","Hc","HH"};
const size_t T=64;
const size_t tStar=4;

int main()
{
  const index_t id({{"",nMes},{"",nSme},{"",nSme}});
  set_njacks(20);
  
  raw_file_t infile("data/data.dat","r");
  vector<djvec_t> data(id.max(),djvec_t(T));
  infile.bin_read(data);
  for(auto& d : data)
    d=d.symmetrized();
  
  grace_file_t resRad("plots/resRad.xmg");
  grace_file_t resSup("plots/resSup.xmg");
  for(size_t iMes=0;iMes<nMes;iMes++)
    {
      vector<pair<double,djvec_t>> mEff;
      for(size_t iSm1=0;iSm1<nSme;iSm1++)
	for(size_t iSm2=0;iSm2<nSme;iSm2++)
	  {
	    const djvec_t tmp=data[id({iMes,iSm1,iSm2})];
	    
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
      
      constexpr size_t nMinToFit=4;
      vector<size_t> listToFit;
      vector<pair<double,size_t>> potentialFitList;
      for(size_t i=2;i<estF.size();i++)
	potentialFitList.emplace_back(estF[i].second.ave(),i);
      sort(potentialFitList.begin(),potentialFitList.end());
      
      set<double> indepX;
      
      for(size_t iScan=0;iScan<potentialFitList.size();iScan++)
	{
	  const size_t i=potentialFitList[iScan].second;
	  const size_t i0=potentialFitList.front().second;
	  const double t=estF[i].first;
	  const djack_t y=estF[i].second;
	  const djack_t y0=estF[i0].second;
	  bool acc=(y.ave()-y0.ave())<sqrt(sqr(y.err())+sqr(y0.err()))*3.5;
	  
	  acc|=indepX.size()<nMinToFit;
	  
	  if(not acc)
	    {
	      size_t nIndepX=0;
	      acc=true;
	      for(size_t j=0;j<listToFit.size();j++)
		{
		  const double& x=estF[listToFit[j]].first;
		  acc&=fabs(t-x)>fabs(t+x)*0.0001;
		}
	    }
	  
	  if(acc)
	    listToFit.push_back(i);
	}
      
      vector<double> x(listToFit.size());
      djvec_t y(listToFit.size());
      
      for(size_t i=0;i<listToFit.size();i++)
	{
	  x[i]=log(estF[listToFit[i]].first);
	  y[i]=estF[listToFit[i]].second;
	}

      const djvec_t p=poly_fit(x,y,2);
      
      grace_file_t plot("plots/est"+mes[iMes]+".xmg");
      plot.set_no_line();
      for(const auto& [extCond,color] :
	    {make_pair(false,grace::RED),{true,grace::GREEN4}})
	{
	  plot.set_all_colors(color);
	  for(size_t i=0;i<estF.size();i++)
	    {
	      bool cond=false;
	      for(size_t iList=0;iList<listToFit.size();iList++)
		{
		  const size_t j=listToFit[iList];
		  cond|=(iList==j);
		}
	      
	      if(extCond==cond)
		{
		  const ave_err_t y=estF[i].second.ave_err();
		  plot.write_ave_err(estF[i].first,y);
		}
	    }
	  plot.new_data_set();
	}
      
      plot.write_polygon([p](const double& x)
      {
	return poly_eval(p,log(x));
      },exp(*min_element(x.begin(),x.end())),exp(*max_element(x.begin(),x.end())),grace::BLUE);
      
      const djack_t x0=exp(-p[1]/(2*p[2]));
      const djack_t yVal=poly_eval(p,log(x0));
      plot.write_ave_err(x0.ave(),yVal.ave_err());
      
      resRad.write_ave_err(mMin.ave_err(),(1/x0).ave_err());
      resSup.write_ave_err(mMin.ave_err(),yVal.ave_err());
    }
  
  map<size_t,vector<pair<size_t,size_t>>> d;
  for(size_t iSm1=0;iSm1<nSme;iSm1++)
    for(size_t iSm2=0;iSm2<nSme;iSm2++)
      d[sme[iSm1]+sme[iSm2]].emplace_back(iSm1,iSm2);
  
  vector<pair<size_t,size_t>> maxPair;
  for(const auto& [dum,p] : d)
    if(p.size()>maxPair.size())
      maxPair=p;
  
  grace_file_t u("/tmp/u.xmg");
  for(const auto& [iSm1,iSm2] : maxPair)
    u.write_xy(sme[iSm1]/((double)sme[iSm1]+sme[iSm2]),effective_mass(data[id({3,iSm1,iSm2})])[24].err());
  
  return 0;
}
