#include "common.hpp"

size_t nHits{};

int main()
{
  setConfs();
  cout<<"nConfs: "<<nConfs<<endl;
  
  set_njacks(10);
  
  T=128;
  
  map<string,vector<vector<double>>> rawData;
  map<string,djvec_t> data;
  
  for(const char* fName : {"2p_Ds_star","2p_Ds_star_dth","2p_Ss","2p_Ss_DTH","4p","2p_Ds"})
    for(const auto& [tag,data] : getRawData(fName))
      rawData[fName+("__"+tag)]=data;
  
  for(const auto& it : rawData)
    for(size_t ri=0;ri<2;ri++)
      {
	const string& tag=it.first;
	cout<<tag<<endl;
	
	const vector<vector<double>>& v=it.second;
	
	if(nHits==0)
	  {
	    nHits=v.size()/nConfs;
	    cout<<"nHits: "<<nHits<<endl;
	    
	    if(nHits==0)
	      CRASH("Unable to set nHits");
	    if(nHits*nConfs!=v.size())
	      CRASH("Length %zu for raw data set %s is not divisible by nConfs %zu",v.size(),tag.c_str(),nConfs);
	  }
	else
	  if(nHits*nConfs!=v.size())
	    CRASH("Length %zu for raw data set %s is not equal to nHits (=%zu) * nConfs(=%zu) = %zu",v.size(),tag.c_str(),nHits,nConfs,nConfs*nHits);
	
	djvec_t& d=data[tag+"_ri"+to_string(ri)]=djvec_t(T);
	
	jackknivesFill(nConfs,
		       [&d,
			&ri,
			&v](const size_t& iConf,
			    const size_t& iClust,
			    const double& weight)
		       {
			 for(size_t iHit=0;iHit<nHits;iHit++)
			   for(size_t t=0;t<T;t++)
			     d[t][iClust]+=v[iHit+nHits*iConf][t+T*ri]*weight;
		       });
	d.clusterize((double)nConfs/njacks)/=nHits;
	d.ave_err().write("plots/"+tag+"_ri"+to_string(ri)+".xmg");
	
	effective_mass(d.symmetrized()).ave_err().write("plots/eff_"+tag+".xmg");
      }
  
  djvec_t& DS=data.at("2p_Ds_star__C_T0__S_T0,__V1V1_ri0");
  djvec_t& DS_DTH=data.at("2p_Ds_star_dth__C_T0__S_DTH_S_T0,__V1V1_ri1");
  
  (DS_DTH/DS).symmetrized().ave_err().write("plots/dDs.xmg");
  
  const djack_t mDs=constant_fit(effective_mass(data.at("2p_Ds__C_T0__S_T0,__P5P5_ri0").symmetrized()),25,50,"plots/fit_Ds.xmg");
  djvec_t ff=data.at("4p__C_T0__L48_S_DTH_S,__V1V1_ri0");
  for(size_t t=0;t<T;t++)
    ff[t]/=exp(-mDs*t);
  const djack_t dcDs=constant_fit(ff,25,30,"plots/ff.xmg");
  djvec_t ffp=data.at("4p__C_T0__L48_S_DTH_S,__V1V1_ri0");
  for(size_t t=0;t<T;t++)
    ffp[t]-=dcDs*exp(-mDs*t);
  effective_mass(-ffp.symmetrized(0)).ave_err().write("plots/ffp.xmg");
  
  djvec_t test(T);
  jackknivesFill(nConfs,
		 [&test,
		  &Ds=rawData.at("2p_Ds_star__C_T0__S_T0,__V1V1"),
		  &ss=rawData.at("2p_Ss__sm_S_sm_T1__sm_S_sm_T1,__P5P5")
		  // &Ds=rawData.at("2p_Ds_star_dth__C_T0__S_DTH_S_T0,__V1V1"),
		  // &ss=rawData.at("2p_Ss_DTH__sm_S_DTH_S_sm_T1__sm_S_sm_T1,__P5P5")
		  ](const size_t& iConf,
										    const size_t& iClust,
										    const double& weight)
		 {
		   for(size_t iHit=0;iHit<nHits;iHit++)
		     for(size_t i=0;i<T;i++)
		       {
			 const size_t t=(i+T+0)%T;
			 const size_t u=(i+T-1)%T;
			 
			 const complex cDs{Ds[iHit+nHits*iConf][t+T*0],Ds[iHit+nHits*iConf][t+T*1]};
			 const complex css{ss[iHit+nHits*iConf][u+T*0],ss[iHit+nHits*iConf][u+T*1]};
			 test[i][iClust]+=(cDs.real()*css.real())*weight;
			 //test[t][iClust]+=css.real()*weight;
			 //test[t][iClust]+=cDs.imag()*css.imag()*weight;
		       }
		 });
  test.clusterize((double)nConfs/njacks)/=nHits;
  test.ave_err().write("plots/test.xmg");
  effective_mass(test.symmetrized(0)).ave_err().write("plots/eff_test.xmg");
  
  return 0;
}
