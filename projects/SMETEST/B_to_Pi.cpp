#include <tranalisi.hpp>

size_t T=48,TH=T/2;
size_t tmin_Pi=7,tmax_Pi=9;
size_t tmin_B=9,tmax_B=11;
size_t tsep_min=8,tsep_max=12;

djvec_t load_part(const string path,bool reim=RE,int tsep=T)
{return read_vec_meas<djvec_t>("data/mes_contr_"+path,T,reim).subset(0,tsep);}

int main()
{
  set_njacks(15);
  
  djvec_t Pi_sme_loc=load_part("Pi_sme_loc_P5P5").symmetrized();
  djvec_t Pi_sme_sme=load_part("Pi_sme_sme_P5P5").symmetrized();
  
  djack_t E_Pi,Z_Pi_sme,Z_Pi_loc;
  two_pts_SL_fit(Z_Pi_sme,Z_Pi_loc,E_Pi,Pi_sme_loc,Pi_sme_sme,TH,tmin_Pi,tmax_Pi,"plots/Pi_eff_mass.xmg");
  
  djvec_t B_sme_loc=load_part("B_sme_loc_P5P5").symmetrized();
  djvec_t B_sme_sme=load_part("B_sme_sme_P5P5").symmetrized();
  
  djack_t M_B,Z_B_sme,Z_B_loc;
  two_pts_SL_fit(Z_B_sme,Z_B_loc,M_B,B_sme_loc,B_sme_sme,TH,tmin_B,tmax_B,"plots/B_eff_mass.xmg");
  

  
  grace_file_t V0P5_mel_plateau_an_out("plots/V0P5_mel_plateaux_an.xmg");
  grace_file_t V0P5_mel_plateau_nu_out("plots/V0P5_mel_plateaux_nu.xmg");
  for(size_t tsep=tsep_min;tsep<=tsep_max;tsep+=2)
    {
      vector<double> x(tsep+1);
      djvec_t yan(tsep+1);
      djvec_t ynu(tsep+1);
      djvec_t P5_V0_P5=load_part(combine("Semi_%02zu_V0P5",tsep),RE,tsep+1);
      djvec_t P5_V0_P5_test=load_part(combine("Demi_%02zu_V0P5",tsep),RE,tsep+1);
      for(size_t t=0;t<=tsep;t++)
	{
	  x[t]=(double)t/tsep;
	  yan[t]=// P5_V0_P5_test[t]/
	    (Z_B_sme*Z_Pi_sme*exp(-M_B*t)*exp(-E_Pi*(tsep-t)));
	  ynu[t]=// P5_V0_P5[t]/
	    (B_sme_loc[t]*Pi_sme_loc[tsep-t]);
	}
      V0P5_mel_plateau_an_out.write_vec_ave_err(x,yan.ave_err());
      V0P5_mel_plateau_nu_out.write_vec_ave_err(x,ynu.ave_err());
      
      P5_V0_P5.ave_err().write(x,combine("plots/V0P5_%02zu.xmg",tsep));
      P5_V0_P5_test.ave_err().write(x,combine("plots/V0P5_%02zu_test.xmg",tsep));
    }
  
  return 0;
}
