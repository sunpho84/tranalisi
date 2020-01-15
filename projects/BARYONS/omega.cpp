#include <tranalisi.hpp>

int T;

djvec_t read(const char *tag,const int iproj,const int iWick)
{
  string path=combine("jacks/bar_alt_contr_%s_proj_%d_Wick_%d",tag,iproj,iWick);
  
  return -read_djvec(path,T).symmetrized(-1);
}

djvec_t read_phi(const char* tag)
{
  djvec_t phi_ss(T);
  
  phi_ss=0.0;
  
  for(int i=1;i<=3;i++)
    phi_ss+=read_djvec(combine("jacks/mes_contr_%s_V%dV%d",tag,i,i),T);
  
  return phi_ss.symmetrized()/3;
 }

int main()
{
  raw_file_t fin("input.txt","r");
  
  set_njacks(fin.read<int>("NJacks"));
  T=fin.read<size_t>("T");
  const int tMin=fin.read<size_t>("Tmin");
  const int tMax=fin.read<size_t>("Tmax");
  const double a_ave=fin.read<double>("a")/0.197;
  const double a_err=fin.read<double>()/0.197;
  const djack_t a(gauss_filler_t(a_ave,a_err,2354235));
  
  for(int im=2;im<=2;im++)
    {
      std::string tag_cpp=combine("SP%dS_SP%dS",im,im);
      const char* tag=tag_cpp.c_str();
      
      //const char SM_SO[2]="S",SM_SI[]="S";
      const djvec_t dir_55=read(tag,0,0);
      const djvec_t exc_55=read(tag,0,1);
      const djvec_t dir_ii=read(tag,1,0);
      const djvec_t exc_ii=read(tag,1,1);
      const djvec_t dir_ij=read(tag,2,0);
      const djvec_t exc_ij=read(tag,2,1);
      
      const djvec_t sig0_55=dir_55;
      const djvec_t nuc_55=dir_55-exc_55;
      const djvec_t boh_55=dir_55-2.0*exc_55;
      
      const djvec_t ii=dir_ii-2*exc_ii;
      const djvec_t ij=dir_ij-2*exc_ij;
      ii.ave_err().write("plots/ii.xmg");
      ij.ave_err().write("plots/ij.xmg");
      // const djvec_t s05=sA-sB;
      // const djvec_t s15=sA+0.5*sB;
      // const djvec_t n05=nA-nB;
      // const djvec_t n15=nA+0.5*nB;
      const djvec_t nuc_alt=ii-ij;
      const djvec_t omega=ii+0.5*ij;
      sig0_55.ave_err().write("plots/corr_sig_55.xmg");
      nuc_55.ave_err().write("plots/corr_nuc_55.xmg");
      nuc_alt.ave_err().write("plots/corr_nuc_alt.xmg");
      omega.ave_err().write("plots/corr_omega.xmg");
      constant_fit(effective_mass(sig0_55,T/2,-1),tMin,tMax,"plots/sig_55.xmg");
      constant_fit(effective_mass(nuc_55,T/2,-1),tMin,tMax,"plots/nuc_55.xmg");
      constant_fit(effective_mass(nuc_alt,T/2,-1),tMin,tMax,"plots/nuc_alt.xmg");
      
      const djack_t aMomega=constant_fit(effective_mass(omega,T/2,-1),tMin,tMax,"plots/omega.xmg");
      const djack_t mOmega=aMomega/a;
      const double mOmega_exp=1.672;
      cout<<"aMOmega: "<<smart_print(aMomega)<<endl;
      cout<<"MOmega: "<<smart_print(mOmega)<<" GeV"<<endl;
      
      if(file_exists("jacks/bar_alt_contr_LPPPL_proj_4_Wick_0"))
	{
	  const djvec_t dir_00=read(tag,3,0);
	  const djvec_t exc_00=read(tag,3,1);
	  // const djvec_t dir_0i=read(tag,4,0);
	  // const djvec_t exc_0i=read(tag,4,1);
	  // const djvec_t dir_i0=read(tag,5,0);
	  // const djvec_t exc_i0=read(tag,5,1);
	  
	  const djvec_t _00=dir_00-2*exc_00;
	  // const djvec_t _0i=dir_0i-2*exc_0i;
	  // const djvec_t _i0=dir_i0-2*exc_i0;
	  constant_fit(effective_mass(_00,T/2,-1),tMin,tMax,"plots/00.xmg");
	  // constant_fit(effective_mass(_0i,T/2,-1),tMin,tMax,"plots/0i.xmg");
	  // constant_fit(effective_mass(_i0,T/2,-1),tMin,tMax,"plots/i0.xmg");
	  const djack_t aM00=constant_fit(effective_mass(_00,T/2,-1),tMin,tMax,"plots/00.xmg");
	  const djack_t m00=aM00/a;
	  cout<<"M00: "<<smart_print(m00)<<" GeV"<<endl;
	}

      const std::string eta_path=combine("jacks/mes_contr_%s_P5P5",tag);
      if(file_exists(eta_path))
	{
	  const djvec_t eta_ss=read_djvec(eta_path,T).symmetrized();
	  const double mEta_ss_exp=0.672;
	  const djack_t aMeta_ss=constant_fit(effective_mass(eta_ss),tMin,tMax,"plots/eta_ss.xmg");
	  
	  const djack_t rat1=aMomega/aMeta_ss;
	  const double rat1_exp=mOmega_exp/mEta_ss_exp;
	  const djack_t rat1_offset=rat1_exp/rat1-1;
	  cout<<rat1_offset.ave_err()<<endl;
	}
      
      if(file_exists(combine("jacks/mes_contr_%s_V1V1",tag)))
	{
	  const djvec_t phi_ss=read_phi(tag);
	  const double mPhi_ss_exp=1.020;
	  const djack_t aMphi_ss=constant_fit(effective_mass(phi_ss),tMin,tMax,"plots/phi_ss.xmg");
	  
	  const djack_t rat2=aMomega/aMphi_ss;
	  const double rat2_exp=mOmega_exp/mPhi_ss_exp;
	  const djack_t rat2_offset=rat2_exp/rat2-1;
	  cout<<rat2_offset.ave_err()<<endl;
	}
	
      // const djack_t rat3=aMphi_ss/aMeta_ss;
      // const double rat3_exp=mPhi_ss_exp/mEta_ss_exp;
      // const djack_t rat3_offset=rat3_exp/rat3-1;
      // cout<<rat3_offset.ave_err()<<endl;
      
      // effective_mass(s05,T/2,-1).ave_err().write("plots/s05.xmg");
      // effective_mass(s15,T/2,-1).ave_err().write("plots/s15.xmg");
      // effective_mass(n05,T/2,-1).ave_err().write("plots/n05.xmg");
      // effective_mass(n15,T/2,-1).ave_err().write("plots/n15.xmg");
      // effective_mass(o05,T/2,-1).ave_err().write("plots/o05.xmg");
      // effective_mass(o15,T/2,-1).ave_err().write("plots/o15.xmg");
      // effective_mass(cB,T/2,-1).ave_err().write("plots/cB.xmg");
      // const djack_t am05=constant_fit(effective_mass(c05,T/2,-1),tMin,tMax,"plots/c05.xmg");
      // const djack_t am15=constant_fit(effective_mass(c15,T/2,-1),tMin,tMax,"plots/c15.xmg");
      
      // const djack_t m05=am05/a;
      // const djack_t m15=am15/a;
      
      // cout<<"M05: "<<smart_print(m05)<<endl;
    }
  
  return 0;
}
