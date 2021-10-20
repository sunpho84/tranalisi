#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "fit.hpp"

#include <GM3/params.hpp>
#include <GM3/perens.hpp>

djvec_t perens_t::determineRenoConst() const
{
  djvec_t Z(2),Zsilv(2);
  
  djack_t mP[2],ZA0[2],ZP5[2];
  const djvec_t corrP5P5[2]=
    {(getAve(0,nSources,idP5P5,0)+
      getAve(0,nSources,idP5P5,2))/2.0,
      getAve(0,nSources,idP5P5,1)};
  // const djvec_t corrA0P5[2]=
  //   {(getAve(0,nSources,idA0P5,0)+
  //     getAve(0,nSources,idA0P5,2))/2.0,
  //     getAve(0,nSources,idA0P5,1)};
  const djvec_t corrP5A0[2]=
    {(-getAve(0,nSources,idP5A0,0)+
      -getAve(0,nSources,idP5A0,2))/2.0,
      -getAve(0,nSources,idP5A0,1)};
  
  for(size_t iMes=0;iMes<2;iMes++)
    {
      djack_t Z2P5separated,MP5separated;
      two_pts_fit(Z2P5separated,MP5separated,corrP5P5[iMes],TH,tMinP5P5[iMes],tMaxP5P5[iMes],"plots/fitP5P5forZmes"+to_string(iMes)+".xmg");
      djvec_t testA0=corrP5A0[iMes];
      const djack_t ZP5separated=sqrt(Z2P5separated);
      for(size_t t=0;t<=TH;t++)
	testA0[t]/=two_pts_corr_fun(ZP5separated,MP5separated,TH,t,-1);
      const djack_t ZA0separated=constant_fit(testA0,tMinP5P5[iMes],tMaxP5P5[iMes],"plots/fitP5A0forZmes"+to_string(iMes)+".xmg");
      console<<"from separated fit, mP: "<<MP5separated.ave_err()<<" , ZA0: "<<ZA0separated.ave_err()<<" , ZP5: "<<ZP5separated.ave_err()<<endl;
      
      two_pts_SL_fit(ZP5[iMes],ZA0[iMes],mP[iMes],corrP5A0[iMes],corrP5P5[iMes],TH,tMinP5P5[iMes],tMaxP5P5[iMes],combine("plots/A0P5FitMes%zu.xmg",iMes),-1,+1);
      console<<"from combined fit, mP: "<<mP[iMes].ave_err()<<" , ZA0: "<<ZA0[iMes].ave_err()<<" , ZP5: "<<ZP5[iMes].ave_err()<<endl;
      
      const djack_t fPfromP=2.0*ZP5[0]*amq/sqr(mP[0]);
      const djack_t fPfromA=ZA0[iMes]/mP[iMes];
      Z[iMes]=fPfromP/fPfromA;
      
      console<<"Z"<<((iMes==0)?"V":"A")<<": "<<Z[iMes].ave_err()<<endl;
    }
  
  const djack_t mPiTM=mP[0]/a;
  const djack_t mPiOS=mP[1]/a;
  
  const djack_t a2=a*a;
  cout<<"TM pion in phys units: "<<a2.ave()<<" "<<mPiTM.ave_err()<<endl;
  cout<<"OS pion in phys units: "<<a2.ave()<<" "<<mPiOS.ave_err()<<endl;
  
  const djack_t ZV_fr_ZA=Z[0]/Z[1];
  console<<"Zv/Za: "<<ZV_fr_ZA.ave_err()<<endl;
  
  const djvec_t derP5A0TM=-symmetric_derivative(corrP5A0[0]);
  derP5A0TM.ave_err().write("plots/derP5A0_regoTM.xmg");
  const djvec_t ZvSilvCorr=2*amq*corrP5P5[0]/derP5A0TM;
  const djack_t ZvSilv=constant_fit(ZvSilvCorr,tMinP5P5[0],tMaxP5P5[0],"plots/ZvSilv.xmg");
  console<<"ZvSilv: "<<ZvSilv.ave_err()<<endl;
  Zsilv[0]=ZvSilv;
  
  const djack_t ZaSilvCorrectingFactor=ZP5[0]/ZP5[1]*mP[1]*sinh(mP[1])/(mP[0]*sinh(mP[0]));
  console<<"ZaSilvCorrectingFactor: "<<ZaSilvCorrectingFactor.ave_err()<<endl;
  
  const djvec_t derP5A0OS=-symmetric_derivative(corrP5A0[1]);
  derP5A0OS.ave_err().write("plots/derP5A0_regoOS.xmg");
  const djvec_t ZaSilvUnCorr=2*amq*corrP5P5[1]/derP5A0OS;
  const djvec_t ZaSilvCorr=ZaSilvUnCorr*ZaSilvCorrectingFactor;
  const djack_t ZaSilvUn=constant_fit(ZaSilvUnCorr,tMinP5P5[0],tMaxP5P5[0],"plots/ZaSilvUncorr.xmg");
  const djack_t ZaSilv=constant_fit(ZaSilvCorr,tMinP5P5[0],tMaxP5P5[0],"plots/ZaSilv.xmg");
  console<<"ZaSilv: "<<ZaSilv.ave_err()<<endl;
  Zsilv[1]=ZaSilv;
  
  return Zsilv;
}
