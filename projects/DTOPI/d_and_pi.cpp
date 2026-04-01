#include "ave_err.hpp"
#include "common.hpp"

int main()
{
  const size_t locNjacks=10;
  
  setConfs();
  T=128;
  set_njacks(locNjacks);
  cout<<"NConfs: "<<nConfs<<endl;
  
  const auto rawData=getRawData("2p");
  map<string,djvec_t> data;
  for(const auto& it : rawData)
    {
      const string& tag=it.first;
      const vector<vector<double>>& v=it.second;
      djvec_t& d=data[tag];
      d.resize(T);
      jackknivesFill(nConfs,[&v,
			     &d](const size_t& iConf,
				 const size_t& iClust,
				 const double& weight)
      {
	for(size_t t=0;t<T;t++)
	  d[t][iClust]+=v[iConf][t]*weight;
      });
      d.clusterize((double)nConfs/njacks);
    }
  
  const size_t n[]{0,1,2,4,8,16,32,64,128,256};
  grace_file_t plot("plots/eff.xmg");
  map<size_t,map<size_t,djack_t>> mt;
    for(size_t i : n)
      for(size_t j : n)
	{
	  const djvec_t& c=data.at(to_string(i)+"_L_"+to_string(i)+"__"+to_string(j)+"_L_"+to_string(j)+"__P5P5");
	  plot.write_vec_ave_err(effective_mass(c.symmetrized()).ave_err());
	  
	  for(size_t t=4;t<12;t++)
	    mt[t][i+j]=effective_mass(c.symmetrized())[t];
	}
  
  grace_file_t mtplot("plots/m.xmg");
  for(size_t t=4;t<12;t++)
    {
      for(const auto& [n,j] : mt[t])
	mtplot.write_ave_err(n,j.ave_err());
      mtplot.new_data_set();
    }
  
  effective_mass(data.at("128_L_128__128_L_128__P5P5").symmetrized()).ave_err().write("plots/m128.xmg");
  
  return 0;
}
