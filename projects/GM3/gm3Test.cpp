#include <tranalisi.hpp>

#define EXTERN
#define INIT_EXTERN_TO(ARGS...) ARGS

#include <aLaLuscherRepresentation.hpp>

using namespace std;

// const double mPi=
//   1.3498e-01;

// const double mRho=
//   7.7526e-01;

// const double g=
//   5.435509;

// const double g2=
//   sqr(g);

const double mPi=
  0.1389;

const double mRho=
  0.7678;

const double g=
  5.3589;

const double g2=
  sqr(g);

// void test()
// {
//   TanPhiCalculator pc;
  
//   HashedTanPhiAndDerivFunction hashedTanPhiFunction("test.hash",4.3,1000);
  
//   grace_file_t original("/tmp/original.xmg");
//   stopwatch_t a;
//   a.start();
//   double t=0;
//   for(double q=0.1;q<4.3;q+=0.001)
//     {
//       //original.write_xy(q*q,pc(q));
//       t=pc(q);
//     }
//   a.stop();
//   cout<<a<<endl;
//       original.write_xy(0,t);

//   original.set_no_symbol();
//   original.new_data_set();
  
//   grace_file_t interpo("/tmp/interpo.xmg");
//   stopwatch_t b;
//   b.start();
//   for(double q=0.1;q<4.3;q+=0.001)
//     {
//       //original.write_xy(q*q,hashedTanPhiFunction(q));
//       t=hashedTanPhiFunction(q);
//     }
//   b.stop();
//   cout<<b<<endl;
//       original.write_xy(0,t);
//   original.set_no_symbol();
//   original.new_data_set();
  
//   grace_file_t error("/tmp/error.xmg");
//   for(double q=0.1;q<4.3;q+=0.001)
//     {
//       original.write_xy(q*q,hashedTanPhiFunction(q)/pc(q)-1.0);
//     }
//   original.set_no_symbol();
//   original.new_data_set();
  
//   ALaLuscherRepresentation interacting(pc,mPi,L,g2,mRho);
  
//   ALaLuscherRepresentation approximated(hashedTanPhiFunction,mPi,L,g2,mRho);
//   //ALaLuscherRepresentation free(mPi,L,0.0,mRho);
  
//   //free.findLevels(10,"/tmp/free.xmg");
//   stopwatch_t c;
//   c.start();
//   interacting.findLevels(10,"/tmp/interacting.xmg");
//   c.stop();
//   cout<<c<<endl;
  
//   stopwatch_t d;
//   d.start();
//   approximated.findLevels(10,"/tmp/approximated.xmg");
//   d.stop();
//   cout<<d<<endl;
// }

int main()
{
  // L=
  //   4.1046/mPi;
  
  L=28.04;
  
  // const vector<double> kNlist{2.0430e-01,2.8392e-01,3.4428e-01,3.7889e-01,4.1793e-01,4.7176e-01,0.1,5.1510e-01};
  
  const double qMax=
    4.3;
  
  const int nIntervals=
    1000;
  
  HashedTanPhiAndDerivFunction hashedPhiAndDerivCalculator("lookupTables/phi.hash","lookupTables/phiDeriv.hash",qMax,nIntervals);
  hashedPhiAndDerivCalculator.plot("plots/phi.xmg","plots/phiDeriv.xmg");
  
  for(double f=1;f>=0.3;f/=1.01)
    {
      ALaLuscherRepresentationCalculator interacting(hashedPhiAndDerivCalculator,mPi*f,L/f,g2,mRho*f);
      
      //for(const auto& p : interacting.getLevelsPars(4).coeffs)
      auto p=interacting.getLevelsPars(4,"/tmp/glp"+to_string(f)+".xmg").coeffs.front();
      cout<<f<<" "<<p.energy/f<<endl;
    }
  
  return 0;
}
