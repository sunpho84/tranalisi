#include <tranalisi.hpp>

size_t T,TH;
size_t L,spatVol;

int main(int narg,char **arg)
{
  raw_file_t input("input","r");
  
  L=input.read<size_t>("L");
  T=input.read<size_t>("T");
  njacks=input.read<size_t>("NJacks");
  const size_t tMin=input.read<size_t>("TMin");
  const size_t tMax=input.read<size_t>("TMax");
  const double extraNorm=input.read<double>("ExtraNorm");
  
  TH=T/2;
  spatVol=L*L*L;
  
  const djvec_t P5P5=read_djvec("jacks/P5P5",T).symmetrized();
  P5P5.ave_err().write("plots/P5P5.xmg");
  
  const djvec_t exchange=read_djvec("jacks/exchange",T).symmetrized();
  exchange.ave_err().write("plots/exchange.xmg");
  const djvec_t exchangeRatio=exchange/P5P5;
  exchangeRatio.ave_err().write("plots/exchangeRatio.xmg");
  
  const djvec_t handcuffs=read_djvec("jacks/handcuffs",T).symmetrized()*extraNorm;
  handcuffs.ave_err().write("plots/handcuffs.xmg");
  const djvec_t handcuffsRatio=handcuffs/P5P5;
  handcuffsRatio.ave_err().write("plots/handcuffsRatio.xmg");
  
  const djvec_t effMass=effective_mass(P5P5);
  const djack_t M=constant_fit(effMass,tMin,tMax,"plots/P5P5Eff.xmg");
  cout<<"M: "<<smart_print(M)<<endl;
  
  const djvec_t effSlopeExchange=effective_slope(exchangeRatio,effMass,TH);
  const djack_t dMExchange=constant_fit(effSlopeExchange,tMin,tMax,"plots/exchangeRat.xmg");
  
  cout<<"Exc: "<<// smart_print
    (dMExchange.ave_err())<<endl;
  
  const djvec_t effSlopeHandcuffs=effective_slope(handcuffsRatio,effMass,TH);
  const djack_t dMHandcuffs=constant_fit(effSlopeHandcuffs,tMin,tMax,"plots/handcuffsRat.xmg");
  
  cout<<"Hands: "<<// smart_print
    (dMHandcuffs.ave_err())<<endl;
  
  const djack_t dMRat=dMHandcuffs/dMExchange;
  cout<<smart_print(dMRat)<<endl;
  
  const djack_t dMTot=dMHandcuffs+dMExchange;
  cout<<"Tot: "<<// smart_print
    (dMTot.ave_err())<<endl;
  
  return 0;
}
