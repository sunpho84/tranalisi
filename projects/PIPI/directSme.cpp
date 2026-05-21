#include "common.hpp"

void computeDirect()
{
  const std::string directDataPath="out";
  const std::vector<std::string> confs=getConfs("confsDirectList.dat",directDataPath,"finished");
  const size_t nConfs=confs.size();
  cout<<"NConfs: "<<nConfs<<endl;
  
  const auto rawData=
     getRaw("rawDirect.dat",
	    "mes_contr_dir",
	    {""},
	    T,
	    directDataPath,
	    confs);
  
  const size_t nHits=rawData.begin()->second.front().size();
  cout<<"NHits: "<<nHits<<endl;
}

int main()
{
  computeDirect();
  
  return 0;
}
