#include "tranalisi.hpp"

int main()
{
  set_njacks(30);
  
  const int niters=2000;
  
  for(int iter=0;iter<niters;iter++)
    {
      gen_t gen(3443645+iter);
      
      double Amean=1e4;
      double Bmean=1e3;
      vector<double> means={Amean,Bmean};
      double eA=1;
      double eB=1;
      double corr=0.1;
      const int nC=300;
      
      double sAA=eA*eA;
      double sBB=eB*eB;
      double sAB=sqrt(sAA*sBB)*corr;
      
      vector<array<double,2>> x;
      array<double,2> derive{};
      for(int i=0;i<nC;i++)
	{
	  const auto y=multivariate({Amean,Bmean},{sAA,sAB,sAB,sBB},gen);
	  
	  x.push_back({y[0]+derive[0],y[1]+derive[1]});
	  
	  for(int j=0;j<2;j++)
	    derive[j]=derive[j]*0.1+y[j]-means[j];
	}
      
      djack_t ab,a,b;
      
      jackknivesFill(nC,
		     [&a,&b,&ab,&x](const size_t iConf,const size_t iJack,const double w)
		     {
		       a[iJack]+=x[iConf][0]*w;
		       b[iJack]+=x[iConf][1]*w;
		       ab[iJack]+=x[iConf][0]*x[iConf][1]*w;
		     });
      
      const double clustSize=(double)nC/njacks;
      a.clusterize(clustSize);
      b.clusterize(clustSize);
      ab.clusterize(clustSize);
      
      const djack_t c=ab-a*b;
      const djack_t d=ab-Amean*Bmean;

      cout.precision(16);
      cout<<"a "<<(a)<<" "<<Amean<<endl;
      cout<<"b "<<(b)<<" "<<Bmean<<endl;
      cout<<"ab "<<(c)<<" "<<(d)<<" "<<sAB<<endl;
    }
  
  return 0;
}
