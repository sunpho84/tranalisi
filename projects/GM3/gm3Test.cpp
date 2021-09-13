#include <tranalisi.hpp>

#include <aLaLuscherRepresentation.hpp>

using namespace std;

int main()
{
  PhiCalculator pc(20.0);
    
  const double mPi=
    1.3498e-01;
  
  const double L=
    4.1046/mPi;
  
  grace_file_t phi("/tmp/phi.xmg");
  grace_file_t rawPhi("/tmp/rawPhi.xmg");
  phi.set_no_symbol();
  rawPhi.set_no_symbol();
  double pr=0;
  for(double q2=0.1;q2<20;q2+=0.01)
    {
      const double q=sqrt(q2);
      const double p=pc(q);
      const double r=pc.rawPhi(q);
      phi.write_xy(q2,p);
      if(fabs(pr-r)>2)
	{
	  rawPhi.new_data_set();
	  rawPhi.set_no_symbol();
	}
      pr=r;
      rawPhi.write_xy(q2,r);
    }
  // const double q=2.0430e-01;//*L/(2*M_PI);
  
  // cout<<" phi: "<<pc(q)<<endl;
  
  
  // return 0;
  
  // int L=
  
  const double mRho=
    7.7526e-01;
  
  const double g=
    5.435509;
  
  const double g2=
    sqr(g);
  
  cout<<"L: "<<L<<endl;
  
  const vector<double> kNlist{2.0430e-01,2.8392e-01,3.4428e-01,3.7889e-01,4.1793e-01,4.7176e-01,0.1,5.1510e-01};
  
  // cout<<"Verify"<<endl;
  // for(int n=1;n<=8;n++)
  //   {
  //     const double& kN=
  // 	kNlist[n-1];
      
  //     const double d11=
  // 	delta11(mPi,mRho,g2,kN);
      
  //     // cout<<kN<<" "<<d11<<endl;
      
  //     const double phi=
  // 	pc(kN*L/(2*M_PI));
      
  //     cout<<n<<": "<<kN<<" "<<tan(d11+phi)<<endl;
  //   }
  
  momentumFinder(mPi,mRho,g2*0,L,10,"/tmp/free.xmg");
  momentumFinder(mPi,mRho,g2,L,10,"/tmp/interacting.xmg");
  
  // grace_file_t lhs("/tmp/lhs.xmg");
  // lhs.set_no_symbol();
  // for(double kN=0;kN<0.8;kN+=0.001)
  //   {
  //     const double d11=
  // 	delta11(mPi,mRho,g2,kN);
      
  //     const double phi=
  // 	pc(kN*L/(2*M_PI));
      
  //     const double y=
  // 	tan(d11+phi);
      
  //     if(fabs(y)>1)
  // 	{
  // 	  lhs.new_data_set();
  // 	  lhs.set_no_symbol();
  // 	}
  //     else
  // 	lhs.write_xy(kN,y);
  //   }
  // const int nIntervals=200;
  
  // PhiCalculator phiCalculator(0.,15,nIntervals);
  
  // grace_file_t out("/tmp/test.xmg");
  // for(double x=0;x<15;x+=0.1)
  //   out.write_xy(x,phiCalculator(x));
  
  return 0;
}
