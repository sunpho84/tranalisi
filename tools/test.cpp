#include <tranalisi.hpp>

int main()
{
  set_njacks(40);
  
  const int n=100000;
  double se1=0,s2e1=0,se2=0,s2e2=0,see2=0,ssk=0,s2sk=0;
  for(int i=0;i<n;i++)
    {
      djack_t t;
      t.fill_gauss(1.3, 0.64, 23144+i);
      
      const double e1=t.err()*njacks/(njacks-1.0);
      const pair<double,double> e2ee2=err_with_err(t);
      
      const pair<double,double> sk=skewness(t);
      ssk+=sk.first;
      s2sk+=sqr(sk.first);
      
      const double e2=e2ee2.first;
      
      const double ee2=e2ee2.second;
      see2+=ee2;
      
      // cerr<<" "<<e1<<endl;
      se1+=e1;
      se2+=e2;
      
      s2e1+=e1*e1;
      s2e2+=e2*e2;
    }
  se1/=n;
  se2/=n;
  s2e1/=n;
  s2e2/=n;
  see2/=n;
  ssk/=n;
  s2sk/=n;
  
  s2e1=sqrt((s2e1-se1*se1)*n/(n-1));
  s2e2=sqrt((s2e2-se2*se2)*n/(n-1));
  s2sk=sqrt((s2sk-ssk*ssk)*n/(n-1));
  
  cout<<se1<<" "<<s2e1<<endl;
  cout<<se2<<" "<<s2e2<<endl;
  cout<<see2<<endl;
  cout<<ssk<<" "<<s2sk<<endl;
  
  djack_t t;
  t.fill_gauss(1.3, 0.64, 23144);
  const auto ae=err_with_err(t);
  cout<<ae.first<<" "<<ae.second<<endl;
  
  return 0;
}
