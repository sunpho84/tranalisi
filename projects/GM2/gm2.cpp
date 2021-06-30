#include <tranalisi.hpp>

vector<double> read(const char* path)
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
  
  return data;
}

index_t id({{"r2",2},{"r1",2},{"corr",61}});

int main(int narg,char** arg)
{
  vector<double> data=read("/tmp/test");
  if(data.size()!=id.max())
    CRASH("mismatch, %zu %zu",data.size(),id.max());
  
  // for(const double& d : data)
  //   cout<<d<<endl;
  
  return 0;
}
