#include <tranalisi.hpp>

#define EXTERN_GLOBAL
 #include "global.hpp"

#include "perens.hpp"

string ensName;

/// Parse the analysis.txt file
void readAnalysis()
{
  raw_file_t input("analysis.txt","r");
  
  set_njacks(input.read<size_t>("NJacks"));
  
  nX=input.read<size_t>("Nx");
  
  useSmeared=input.read<size_t>("Smeared");
  subtractZeroV=input.read<size_t>("Subtract");
  
  input.expect("EnsembleList");
  ensName=input.read<string>();
}

int main()
{
  readAnalysis();
  
  perens_t ens(ensName);
  
  cout<<"m1: "<<ens.mu1<<" , q1: "<<q1<<endl;
  cout<<"m2: "<<ens.mu2<<" , q2: "<<q2<<endl;
  
  const djack_t aFp=(ens.mu1+ens.mu2)*ens.ZP[0]/sqr(ens.mP);
  cout<<"a*fP: "<<aFp.ave_err()<<endl;
  
  auto g=[&ens](const size_t iAV,const size_t iX) -> djvec_t
	 {
	   constexpr int coeff[2]={-1,+1};
	   
	   return
	     (ens.c3[ens.iC3({0,iAV,iX})]*q1+
	      coeff[iAV]*
	      ens.c3[ens.iC3({1,iAV,iX})]*q2)/3;
	 };
  
  djvec_t Rtot[2],fA(nX),fV(nX);
  for(size_t iAV=0;iAV<2;iAV++)
    Rtot[iAV]=djvec_t(nX);
  
  for(size_t iX=0;iX<nX;iX++)
    {
      cout<<"X["<<iX<<"]: "<<ens.xG[iX].ave_err()<<" , Eg: "<<ens.eG[iX]<<endl;
      
      const djvec_t RA=g(HA,iX)/g(HA,0)-1.0;
      const djvec_t RV=(g(HV,iX)-g(HV,0)*subtractZeroV)/g(HA,0);
      
      g(HA,iX).ave_err().write(ens.plotsPath+combine("C_A_x%zu.xmg",iX));
      g(HV,iX).ave_err().write(ens.plotsPath+combine("C_V_x%zu.xmg",iX));
      
      // for(size_t iIns=0;iIns<2;iIns++)
      // 	ens.c3[ens.iC3({iIns,HV,iX})].ave_err().write(plotsPath+combine("Vrec_x%zu_ins%zu.xmg",iX,iIns));
      
      djack_t y[2][2];
      for(size_t iIns=0;iIns<2;iIns++)
	for(size_t iAV=0;iAV<2;iAV++)
	  {
	    const djvec_t e=ens.c3[ens.iC3({iIns,iAV,iX})]/ens.norm;
	    y[iIns][iAV]=constant_fit(e,16,20,ens.plotsPath+combine("H_%s_ins%zu_x%zu.xmg",((iAV==HA)?"A":"V"),iIns,iX));
	  }
      
      for(size_t iAV=0;iAV<2;iAV++)
	{
	  const djvec_t e=g(HA,iX)/ens.norm;
	  constant_fit(e,16,20,ens.plotsPath+combine("H_%s_x%zu.xmg",((iAV==HA)?"A":"V"),iX));
	}
      
      Rtot[HA][iX]=constant_fit(RA,ens.tMinA[iX],ens.tMaxA[iX],ens.plotsPath+combine("R_A_x%zu.xmg",iX));
      Rtot[HV][iX]=constant_fit(RV,ens.tMinV[iX],ens.tMaxV[iX],ens.plotsPath+combine("R_V_x%zu.xmg",iX));
      
      fA[iX]=2*aFp*Rtot[HA][iX]/(ens.mP*ens.xG[iX]);
      fV[iX]=aFp*Rtot[HV][iX]/ens.kZ[iX];
    }
  
  grace_file_t fAplot(ens.plotsPath+"fAFunX.xmg");
  grace_file_t fVplot(ens.plotsPath+"fVFunX.xmg");
  
  fAplot.write_vec_ave_err(ens.xG.ave(),fA.ave_err());
  fVplot.write_vec_ave_err(ens.xG.ave(),fV.ave_err());
  
  return 0;
}
