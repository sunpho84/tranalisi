#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <unistd.h>

int main(int narg,char **arg)
{
  if(narg<2) CRASH("Use: %s file",arg[0]);
  
  vector<tuple<double,double,double>> tmp;
  
  ifstream file(arg[1]);
  if(not file.good()) CRASH("Unable to open %s",arg[1]);
  
  double x,y,ey;
  while(file>>x>>y>>ey)
    tmp.push_back(make_tuple(x,y,ey));
  
  int d=1;
  
  vector <double> Al(2*d+1,0.0);
  vector<double> c(d+1,0.0);
  
  for(int p=0;p<(int)tmp.size();p++)
    {
      //calculate the weight
      double w=pow(get<2>(tmp[p]),-2);
      
      //compute Al and c
      for(int f=0;f<=2*d;f++)
	{
	  Al[f]+=w;
	  if(f<=d) c[f]+=get<1>(tmp[p])*w;
	  w*=get<0>(tmp[p]);
	}
    }
  
  vector<double> A((d+1)*(d+1));
  for(int i=0;i<=d;i++)
    for(int j=0;j<=d;j++)
      A[i*(d+1)+j]=Al[i+j];
  
  vector<double> res=lin_solve<vector<double>,double>(A,c);
  
  cout<<res[0]<<endl;
  
  return 0;
}
