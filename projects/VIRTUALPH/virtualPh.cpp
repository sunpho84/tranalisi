#include <tranalisi.hpp>

#include "perens.hpp"

constexpr int q1=2,q2=-1;

int main()
{
  raw_file_t input("analysis.txt","r");
  
  set_njacks(input.read<size_t>("NJacks"));
  
  useSmeared=input.read<size_t>("Smeared");
  const size_t subtract=input.read<size_t>("Subtract");
  const size_t tMinA=input.read<size_t>("TIntsA");
  const size_t tMaxA=input.read<size_t>();
  const size_t tMinV=input.read<size_t>("TIntsV");
  const size_t tMaxV=input.read<size_t>();
  
  perens_t ens;
  
  cout<<"m1: "<<ens.mu1<<" , q1: "<<q1<<endl;
  cout<<"m2: "<<ens.mu2<<" , q2: "<<q2<<endl;
  
  const djack_t aFp=(ens.mu1+ens.mu2)*ens.ZPloc/sqr(ens.mP);
  cout<<"a*fP: "<<aFp.ave_err()<<endl;
  
  auto g=[&ens](const size_t iAV,const size_t iX) -> djvec_t
	 {
	   constexpr int coeff[2]={-1,+1};
	   
	   return
	     (ens.c3[ens.iC3({0,iAV,iX})]*q1+
	      coeff[iAV]*
	      ens.c3[ens.iC3({1,iAV,iX})]*q2)/3;
	 };
  
  djvec_t Rtot[2],fA(ens.nX),fV(ens.nX);
  for(size_t iAV=0;iAV<2;iAV++)
    Rtot[iAV]=djvec_t(ens.nX);
  
  for(size_t iX=0;iX<ens.nX;iX++)
    {
      cout<<"X["<<iX<<"]: "<<ens.xG[iX].ave_err()<<" , Eg: "<<ens.eG[iX]<<endl;
      
      const djvec_t RA=g(HA,iX)/g(HA,0)-1.0;
      const djvec_t RV=(g(HV,iX)-g(HV,0)*subtract)/g(HA,0);
      
      g(HA,iX).ave_err().write(combine("plots/C_A_x%zu.xmg",iX));
      g(HV,iX).ave_err().write(combine("plots/C_V_x%zu.xmg",iX));
      
      // for(size_t iIns=0;iIns<2;iIns++)
      // 	ens.c3[ens.iC3({iIns,HV,iX})].ave_err().write(combine("plots/Vrec_x%zu_ins%zu.xmg",iX,iIns));
      
      djack_t y[2][2];
      for(size_t iIns=0;iIns<2;iIns++)
	for(size_t iAV=0;iAV<2;iAV++)
	  {
	    const djvec_t e=ens.c3[ens.iC3({iIns,iAV,iX})]/ens.norm;
	    y[iIns][iAV]=constant_fit(e,16,20,combine("plots/H_%s_ins%zu_x%zu.xmg",((iAV==HA)?"A":"V"),iIns,iX));
	  }
      
      for(size_t iAV=0;iAV<2;iAV++)
	{
	  const djvec_t e=g(HA,iX)/ens.norm;
	  constant_fit(e,16,20,combine("plots/H_%s_x%zu.xmg",((iAV==HA)?"A":"V"),iX));
	}
      
      Rtot[HA][iX]=constant_fit(RA,tMinA,tMaxA,combine("plots/R_A_x%zu.xmg",iX));
      Rtot[HV][iX]=constant_fit(RV,tMinV,tMaxV,combine("plots/R_V_x%zu.xmg",iX));
      
      fA[iX]=2*aFp*Rtot[HA][iX]/(ens.mP*ens.xG[iX]);
      fV[iX]=aFp*Rtot[HV][iX]/ens.kZ[iX];
    }
  
  grace_file_t fAplot("plots/fAFunX.xmg");
  grace_file_t fVplot("plots/fVFunX.xmg");
  
  fAplot.write_vec_ave_err(ens.xG.ave(),fA.ave_err());
  fVplot.write_vec_ave_err(ens.xG.ave(),fV.ave_err());
  
  return 0;
}
