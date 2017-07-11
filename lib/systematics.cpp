#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <systematics.hpp>

vector<double> syst_analysis_sep_bis(const vector<ave_err_t> &v,const index_t &ind)
{
  if(ind.max()!=v.size()) CRASH("v has size %zu, fact has product %zu",v.size(),ind.max());
  
  vector<double> out(ind.rank(),0.0);
  vector<double> partials=out; //<! ignoring mixed terms
  
  for(int iter=0;iter<2;iter++)
    for(size_t i=0;i<v.size();i++)
      for(size_t j=0;j<v.size();j++)
	{
	  vector<size_t> ci=ind(i);
	  vector<size_t> cj=ind(j);
	  // cout<<ci[0]<<" "<<ci[1]<<" "<<ci[2]<<endl;
	  // cout<<cj[0]<<" "<<cj[1]<<" "<<cj[2]<<endl;
	  double delta=v[i].ave()*(v[i].ave()-v[j].ave());
	  
	  //get the list of different components
	  vector<size_t> mudiff;
	  for(size_t mu=0;mu<ind.rank();mu++)
	    if(ci[mu]!=cj[mu])
	      mudiff.push_back(mu);
	  //cout<<ndiff<<endl;
	  if(iter==0) {if(mudiff.size()==1) partials[mudiff[0]]+=delta;}
	  else
	    if(mudiff.size()!=0)
	      {
		double norm=0;
		for(auto &mu : mudiff) norm+=partials[mu];
		if(norm) for(auto &mu : mudiff) out[mu]+=delta*partials[mu]/norm;
	      }
	}
  
  //normalize
  for(size_t mu=0;mu<ind.rank();mu++) out[mu]=sqrt(out[mu]/(v.size()*(v.size()-1)));
  
  return out;
}
