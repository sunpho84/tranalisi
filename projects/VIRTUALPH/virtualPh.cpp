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

djack_t excFit(const djvec_t& c,const size_t tSep,const string& plotPath)
{
  jack_fit_t fit;
  
  djvec_t pGuess(5);
  pGuess[0]=c[tSep/2];
  djvec_t y=c-c[tSep/2];
  djvec_t i=y;
  for(size_t t=0;t<y.size();t++)
    i[t]/=y[(t+1)%y.size()];
  i=log(i);
  pGuess[2]=constant_fit(i,1,tSep/2,plotPath+"2");
  i=y;
  for(size_t t=0;t<y.size();t++)
    i[t]/=exp(-pGuess[2]*t);
  pGuess[1]=constant_fit(i,1,tSep/2,plotPath+"1");
  i=y;
  for(size_t t=0;t<y.size();t++)
    i[t]=y[(t-1+y.size())%y.size()]/y[t];
  i=-log(i);
  pGuess[4]=constant_fit(i,tSep/2,tSep-2,plotPath+"4");
  i=y;
  for(size_t t=0;t<y.size();t++)
    i[t]/=exp(-pGuess[4]*(tSep-t));
  pGuess[3]=constant_fit(i,tSep/2,tSep-1,plotPath+"3");
  
  cout<<plotPath<<" guess: "<<endl;
  cout<<pGuess.ave_err()<<endl;
  
  djvec_t pFit(5);
  fit.add_fit_par(pFit[0],"M",pGuess[0].ave_err());
  fit.add_fit_par(pFit[1],"A1",pGuess[1].ave_err());
  fit.add_fit_par(pFit[2],"M1",pGuess[2].ave_err());
  fit.add_fit_par(pFit[3],"A2",pGuess[3].ave_err());
  fit.add_fit_par(pFit[4],"M2",pGuess[4].ave_err());
  
  const size_t dt=1;
  for(size_t it=0+dt;it<tSep-dt;it++)
    {
      fit.add_point(c[it],
		    [=]
		    (const vector<double>& p,int ijack)
		    {
		      return p[0]+p[1]*exp(-p[2]*it)+p[3]*exp(-p[4]*(tSep-it));
		    });
    }
  
  auto status=fit.fit();
  
  cout<<plotPath<<" res: "<<endl;
  cout<<pFit.ave_err()<<endl;
  
  grace_file_t plot(plotPath);
  plot.write_vec_ave_err(c.ave_err());
  for(auto& p : {pGuess,pFit})
    plot.write_polygon([&p,tSep](const double x)->djack_t
		       {
			 return p[0]+p[1]*exp(-p[2]*x)+p[3]*exp(-p[4]*(tSep-x));
		       }
	,0+dt,tSep-dt,grace::RED);
  
  return pFit[0];
}

int main()
{
  readAnalysis();
  
  perens_t ens(ensName);
  
  ens.norm.ave_err().write("plots/norm.xmg");
  
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
  
  djack_t y[nX][2][2];
  for(size_t iX=0;iX<nX;iX++)
    {
      cout<<"X["<<iX<<"]: "<<ens.xG[iX].ave_err()<<" , Eg: "<<ens.eG[iX]<<endl;
      
      const djvec_t RA=g(HA,iX)/g(HA,0)-1.0;
      const djvec_t RV=(g(HV,iX)-g(HV,0)*subtractZeroV)/g(HA,0);
      
      g(HA,iX).ave_err().write(ens.plotsPath+combine("C_A_x%zu.xmg",iX));
      g(HV,iX).ave_err().write(ens.plotsPath+combine("C_V_x%zu.xmg",iX));
      
      // for(size_t iIns=0;iIns<2;iIns++)
      // 	ens.c3[ens.iC3({iIns,HV,iX})].ave_err().write(plotsPath+combine("Vrec_x%zu_ins%zu.xmg",iX,iIns));
      
      for(size_t iIns=0;iIns<2;iIns++)
	for(size_t iAV=0;iAV<2;iAV++)
	  {
	    const djvec_t e=ens.c3[ens.iC3({iIns,iAV,iX})]/ens.norm;
	    y[iX][iIns][iAV]=constant_fit(e,(iIns==1)?12:16,(iIns==1)?16:20,ens.plotsPath+combine("H_%s_ins%zu_x%zu.xmg",((iAV==HA)?"A":"V"),iIns,iX));
	    
	    // const djack_t ye=excFit(e,ens.T/2,ens.plotsPath+combine("H_%s_ins%zu_x%zu_exc.xmg",((iAV==HA)?"A":"V"),iIns,iX));
	    // cout<<"Fit: "<<y[iX][iIns][iAV].ave_err()<<" "<<ye.ave_err()<<endl;
	  }
      
      for(size_t iAV=0;iAV<2;iAV++)
	{
	  const djvec_t e=g(iAV,iX)/ens.norm;
	  constant_fit(e,16,20,ens.plotsPath+combine("H_%s_x%zu.xmg",((iAV==HA)?"A":"V"),iX));
	}
      
      Rtot[HA][iX]=constant_fit(RA,ens.tMinA[iX],ens.tMaxA[iX],ens.plotsPath+combine("R_A_x%zu.xmg",iX));
      Rtot[HV][iX]=constant_fit(RV,ens.tMinV[iX],ens.tMaxV[iX],ens.plotsPath+combine("R_V_x%zu.xmg",iX));
      
      const djack_t RAreco=(y[iX][0][HA]*q1-y[iX][1][HA]*q2)/(y[0][0][HA]*q1-y[0][1][HA]*q2)-1;
      cout<<"RbarA: "<<Rtot[HA][iX].ave_err()<<" , reco: "<<RAreco.ave_err()<<endl;
      
      fA[iX]=2*aFp*Rtot[HA][iX]/(ens.mP*ens.xG[iX]);
      fV[iX]=aFp*Rtot[HV][iX]/ens.kZ[iX];
    }
  
  grace_file_t fAplot(ens.plotsPath+"fAFunX.xmg");
  grace_file_t fVplot(ens.plotsPath+"fVFunX.xmg");
  
  fAplot.write_vec_ave_err(ens.xG.ave(),fA.ave_err());
  fVplot.write_vec_ave_err(ens.xG.ave(),fV.ave_err());
  
  return 0;
}
