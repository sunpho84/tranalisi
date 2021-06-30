#include <tranalisi.hpp>

const size_t T=128;
const index_t id({{"r2",2},{"r1",2},{"corr",61},{"T",T},{"ri",2}});

std::pair<vector<double>,vector<double>> read(const char* path)
{
  ifstream input(path);
  
  if(not input.good())
    CRASH("unable to open %s",path);
  
  vector<double> data;
  
  string is;
  while(input>>is)
    {
      istringstream iss(is);
      double d;
      if(iss>>d)
	data.push_back(d);
    }
  
  if(data.size()!=id.max())
    CRASH("mismatch, %zu %zu",data.size(),id.max());
  
  vector<double> VV(T/2+1,0.0),PP(T/2+1,0.0);
  for(size_t r1=0;r1<2;r1++)
    //for(size_t r2=0;r2<2;r2++)
      for(size_t _t=0;_t<T;_t++)
	{
	  const size_t t=std::min(T-_t,_t);
	  const int c=(t==0)?2:1;
	  
	  const size_t r2=r1;
	  for(size_t corr : {35,36,37})
	    VV[t]+=data[id({r2,r1,corr,_t,0})]*c;
	  PP[t]+=data[id({r2,r1,5,_t,0})]*c;
	}
  
  for(auto& pp : PP)
    pp/=2*2;
  
  for(auto& vv : VV)
    vv/=2*2*3;
  
  return {PP,VV};
}

int main(int narg,char** arg)
{
  if(narg<2)
    CRASH("use %s file",arg[0]);
  
  vector<double> PP,VV;
  tie(PP,VV)=read(arg[1]);
  
  for(auto& pp : PP)
    cout<<pp<<endl;
  
  // for(const double& d : data)
  //   cout<<d<<endl;
  
  return 0;
}
