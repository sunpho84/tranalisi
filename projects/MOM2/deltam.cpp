#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <MOM2/contractions.hpp>
#include <MOM2/deltam.hpp>
#include <MOM2/perens.hpp>
#include <MOM2/prop.hpp>
#include <MOM2/timings.hpp>

perens_t& perens_t::get_deltam()
{
  vector<tuple<djvec_t*,string,bool>> delta_tasks{
    {&deltam_cr,"cr",pars::use_deltam_cr_ct},
    {&deltam_tm,"tm",pars::use_deltam_tm_ct}};
  
  //if file exists open it, otherwise compute it
  const string deltam_path=dir_path+"/deltam.dat";
  if(file_exists(deltam_path))
    {
      cout<<"File "<<deltam_path<<" found, opening"<<endl;
      raw_file_t deltam_file(deltam_path,"r");
      for(auto &dtu : delta_tasks)
	get<0>(dtu)->bin_read(deltam_file);
    }
  else
    {
      cout<<"File "<<deltam_path<<" not found, computing"<<endl;
      prepare_list_of_confs();
      
      deltam_time.start();
      for(size_t im=0;im<nm;im++)
	for(size_t r=0;r<nr;r++)
	  compute_deltam(im,r);
      deltam_time.stop();
      
      raw_file_t deltam_file(deltam_path,"w");
      for(auto &dtu : delta_tasks)
	get<0>(dtu)->bin_write(deltam_file);
    }
  
  //Prepare plot
  for(auto &dtu : delta_tasks)
    {
      djvec_t &deltam=*get<0>(dtu);
      const string &tag=get<1>(dtu);
      
      grace_file_t plot(dir_path+"/plots/deltam_"+tag+"_vs_am.xmg");
      
      for(size_t r=0;r<nr;r++)
	{
	  plot.new_data_set();
	  for(size_t im=0;im<nm;im++) plot.write_ave_err(am[im],deltam[im_r_ind({im,r})].ave_err());
	}
    }
  
  //Print, putting it to zero if not asked to use
  for(size_t im=0;im<nm;im++)
    for(size_t r=0;r<nr;r++)
      {
	cout<<"Deltam[m="<<im<<",r="<<r<<"]\t";
	for(auto &dtu : delta_tasks)
	  {
	    djvec_t &deltam=*get<0>(dtu);
	    const string &tag=get<1>(dtu);
	    const bool use=get<2>(dtu);
	    
	    cout<<" "<<tag<<": "<<deltam[im_r_ind({im,r})];
	    if(not use)
	      {
		deltam[im_r_ind({im,r})]=0.0;
		cout<<" -> "<<0.0;
	      }
	    cout<<"\t";
	  }
	cout<<endl;
      }
  
  return *this;
}

perens_t& perens_t::get_mPCAC()
{
  mPCAC.resize(nm);
  
  //if file exists open it, otherwise compute it
  const string mPCAC_path=dir_path+"/mPCAC.dat";
  if(file_exists(mPCAC_path))
    {
      cout<<"File "<<mPCAC_path<<" found, opening"<<endl;
      mPCAC.bin_read(mPCAC_path);
    }
  else
    {
      cout<<"File "<<mPCAC_path<<" not found, computing"<<endl;
      
      prepare_list_of_confs();
      for(size_t im=0;im<nm;im++)
	  mPCAC[im]=compute_mPCAC(to_string(im));
      mPCAC.bin_write(mPCAC_path);
    }
  
  grace_file_t mPCAC2_plot(dir_path+"/plots/mPCAC.xmg");
  mPCAC2_plot.new_data_set();
  for(size_t im=0;im<nm;im++)
    mPCAC2_plot.write_ave_err(am[im],mPCAC[im].ave_err());
  
  //sea
  const string mPCAC_sea_path=dir_path+"/mPCAC_sea.dat";
  if(im_sea>=0 and im_sea<(int)nm) mPCAC_sea=mPCAC[im_im_ind({(size_t)im_sea,(size_t)im_sea})];
  else
    if(file_exists(mPCAC_sea_path))
      {
	cout<<"File "<<mPCAC_sea_path<<" found, opening"<<endl;
	mPCAC_sea.bin_read(mPCAC_sea_path);
      }
    else
      {
	cout<<"File "<<mPCAC_sea_path<<" not found, computing"<<endl;
	
	prepare_list_of_confs();
	mPCAC_sea=compute_mPCAC("sea");
	mPCAC_sea.bin_write(mPCAC_sea_path);
      }
  mPCAC2_plot.write_ave_err(0,mPCAC_sea.ave_err());
  cout<<"Sea PCAC: "<<smart_print(mPCAC_sea.ave_err())<<endl;
  
  return *this;
}

void perens_t::compute_deltam(const size_t im,const size_t rfw)
{
  //ext_reim and rpar is relative to non-inserted
  auto get=[im,this]
    (size_t kbw,size_t kfw,const string &ID,const size_t ext_reim,const int tpar,const size_t rfw,const int rdiff)
    {
      const string tag_bw=kind_tag[kbw];
      const string tag_fw=kind_tag[kfw];
      
      const size_t rbw=(rfw+rdiff)%nr;
      
      //Compute the coefficient
      dcompl_t c_coeff=
      conj(m_r_mom_conf_qprops_t::coeff_to_read(kbw,rbw))*
      m_r_mom_conf_qprops_t::coeff_to_read(kfw,rfw);
      
      //Include -i if asking the imaginary part
      if(ext_reim==1) c_coeff*=dcompl_t(0.0,-1.0);
      
      const string name="M"+to_string(im)+"_R"+to_string(rbw)+"_"+tag_bw+"_M"+to_string(im)+"_R"+to_string(rfw)+"_"+tag_fw;
      const djvec_t res=get_contraction(name,ID,c_coeff,tpar);
      res.ave_err().write(dir_path+"/plots/"+ID+"_"+name+".xmg");
      
      return res;
    };
  
  //measure mcrit according to eq.3 of hep-lat/0701012
  {
    const djvec_t P5P5_00=get(_LO,_LO,"P5P5",RE,EVN,rfw,0);
    const djvec_t V0P5_00=get(_LO,_LO,"V0P5",IM,ODD,rfw,0);
    const djvec_t m_cr_corr=2.0*forward_derivative(V0P5_00)/(2.0*P5P5_00);
    const djvec_t m_cr_corr_symm=2.0*symmetric_derivative(V0P5_00)/(2.0*P5P5_00);
    const djack_t m_cr=constant_fit(m_cr_corr,tmin,tmax,dir_path+"/plots/m_cr_"+to_string(im)+".xmg");
    const djack_t m_cr_symm=constant_fit(m_cr_corr_symm,tmin,tmax,dir_path+"/plots/m_cr_symm_"+to_string(im)+".xmg");
    cout<<"m_cr[m="<<im<<",rfw="<<rfw<<"]: "<<smart_print(m_cr.ave_err())<<", symm: "<<smart_print(m_cr_symm.ave_err())<<endl;
  }
  
  if(pars::use_QED)
    {
      int rdiff=0;
      
      //load corrections
      const djvec_t P5P5_00=get(_LO,_LO,"P5P5",RE,UNK,rfw,rdiff);
      const djvec_t P5P5_LL=get(_F,_F,"P5P5",RE,UNK,rfw,rdiff);
      const djvec_t P5P5_0M=get(_LO,_FF,"P5P5",RE,UNK,rfw,rdiff);
      const djvec_t P5P5_M0=get(_FF,_LO,"P5P5",RE,UNK,rfw,rdiff);
      const djvec_t P5P5_0T=get(_LO,_T,"P5P5",RE,UNK,rfw,rdiff);
      const djvec_t P5P5_T0=get(_T,_LO,"P5P5",RE,UNK,rfw,rdiff);
      //load the derivative wrt counterterm
      const djvec_t P5P5_0P=get(_LO,_P,"P5P5",RE,UNK,rfw,rdiff);
      const djvec_t P5P5_P0=get(_P,_LO,"P5P5",RE,UNK,rfw,rdiff);
      //load the derivative wrt mass
      const djvec_t P5P5_0S=get(_LO,_S,"P5P5",RE,UNK,rfw,rdiff);
      const djvec_t P5P5_S0=get(_S,_LO,"P5P5",RE,UNK,rfw,rdiff);
      
      //load corrections
      const djvec_t V0P5_00=get(_LO,_LO,"V0P5",IM,UNK,rfw,rdiff);
      const djvec_t V0P5_LL=get(_F,_F,"V0P5",IM,UNK,rfw,rdiff);
      const djvec_t V0P5_0M=get(_LO,_FF,"V0P5",IM,UNK,rfw,rdiff);
      const djvec_t V0P5_M0=get(_FF,_LO,"V0P5",IM,UNK,rfw,rdiff);
      const djvec_t V0P5_0T=get(_LO,_T,"V0P5",IM,UNK,rfw,rdiff);
      const djvec_t V0P5_T0=get(_T,_LO,"V0P5",IM,UNK,rfw,rdiff);
      //load the derivative wrt counterterm
      const djvec_t V0P5_0P=get(_LO,_P,"V0P5",IM,UNK,rfw,rdiff);
      const djvec_t V0P5_P0=get(_P,_LO,"V0P5",IM,UNK,rfw,rdiff);
      //load the derivative wrt mass
      const djvec_t V0P5_0S=get(_LO,_S,"V0P5",IM,UNK,rfw,rdiff);
      const djvec_t V0P5_S0=get(_S,_LO,"V0P5",IM,UNK,rfw,rdiff);
      
      const djvec_t
      V0P5_LO=
	V0P5_00,
      V0P5_QED=
	V0P5_0T+
	V0P5_0M+
	V0P5_T0+
	V0P5_M0+
	V0P5_LL,
      V0P5_S=
	V0P5_0S+
	V0P5_S0,
      V0P5_P=
	V0P5_0P+
	V0P5_P0;
      const djvec_t
      P5P5_LO=
	P5P5_00,
      P5P5_QED=
	P5P5_0T+
	P5P5_0M+
	P5P5_T0+
	P5P5_M0+
	P5P5_LL,
      P5P5_S=
	P5P5_0S+
	P5P5_S0,
      P5P5_P=
	P5P5_0P+
	P5P5_P0;
      
      const size_t T=P5P5_LO.size();
      const djvec_t eff_P5P5=effective_mass(P5P5_LO.symmetrized());
      
      djvec_t a=djvec_t(symmetric_derivative(V0P5_QED)/P5P5_LO-symmetric_derivative(V0P5_LO)/sqr(P5P5_LO)*P5P5_QED).symmetrized().subset(0,T/2-1);
      djvec_t b=djvec_t(symmetric_derivative(V0P5_S)/P5P5_LO-symmetric_derivative(V0P5_LO)/sqr(P5P5_LO)*P5P5_S).symmetrized().subset(0,T/2-1);
      djvec_t c=djvec_t(symmetric_derivative(V0P5_P)/P5P5_LO-symmetric_derivative(V0P5_LO)/sqr(P5P5_LO)*P5P5_P).symmetrized().subset(0,T/2-1);
      djvec_t d=effective_slope(djvec_t(P5P5_QED/P5P5_LO).symmetrized(),eff_P5P5,T/2);
      djvec_t e=effective_slope(djvec_t(P5P5_S/P5P5_LO).symmetrized(),eff_P5P5,T/2);
      djvec_t f=effective_slope(djvec_t(P5P5_P/P5P5_LO).symmetrized(),eff_P5P5,T/2);
      
      a.ave_err().write(dir_path+"/plots/deltam_a_m"+to_string(im)+"_rfw"+to_string(rfw)+".xmg");
      b.ave_err().write(dir_path+"/plots/deltam_b_m"+to_string(im)+"_rfw"+to_string(rfw)+".xmg");
      c.ave_err().write(dir_path+"/plots/deltam_c_m"+to_string(im)+"_rfw"+to_string(rfw)+".xmg");
      d.ave_err().write(dir_path+"/plots/deltam_d_m"+to_string(im)+"_rfw"+to_string(rfw)+".xmg");
      e.ave_err().write(dir_path+"/plots/deltam_e_m"+to_string(im)+"_rfw"+to_string(rfw)+".xmg");
      f.ave_err().write(dir_path+"/plots/deltam_f_m"+to_string(im)+"_rfw"+to_string(rfw)+".xmg");
      
      enum strategy_t{only_cr,both_cr_and_tm};
      const size_t t_probe=T/8;
      const double tol=1e-13;
      const strategy_t strategy=
	(fabs(b[t_probe].ave())>tol or
	 fabs(e[t_probe].ave())>tol)
	?
	both_cr_and_tm
	:
	only_cr;
      
      if(strategy==both_cr_and_tm)
	{
	  const djvec_t old_corr=djvec_t(symmetric_derivative(V0P5_QED)/symmetric_derivative(V0P5_P)).symmetrized();
	  const djack_t old=constant_fit(old_corr,tmin,tmax,dir_path+"/plots/old_deltam_cr_m"+to_string(im)+"_rfw"+to_string(rfw)+".xmg");
	  
	  const djvec_t den=b*f-c*e;
	  const djvec_t deltam_tm_corr=(-a*f+c*d)/den;
	  const djvec_t deltam_cr_corr=(-b*d+a*e)/den;
	  
	  const size_t imr=im_r_ind({im,rfw});
	  
	  for(auto &oit : vector<tuple<djvec_t*,const djvec_t*,string>>{{&deltam_cr,&deltam_cr_corr,"cr"},{&deltam_tm,&deltam_tm_corr,"tm"}})
	    {
	      djack_t &out=(*::get<0>(oit))[imr];
	      const djvec_t &in=*::get<1>(oit);
	      const string tag=::get<2>(oit);
	      
	      out=constant_fit(in,tmin,tmax,dir_path+"/plots/deltam_"+tag+"_m"+to_string(im)+"_rfw"+to_string(rfw)+".xmg");
	    }
	}
      else
	{
	  cout<<"Switching to determining only deltam_cr"<<endl;
	  
	  djvec_t a=djvec_t(symmetric_derivative(V0P5_QED)/P5P5_LO).symmetrized().subset(0,T/2-1);
	  djvec_t c=djvec_t(symmetric_derivative(V0P5_P)/P5P5_LO).symmetrized().subset(0,T/2-1);
	  djvec_t d=P5P5_QED.symmetrized();
	  djvec_t f=P5P5_P.symmetrized();
	  
	  const djvec_t ma_fr_c=-a/c;
	  const djvec_t md_fr_f=-d/f;
	  
	  deltam_cr=constant_fit(md_fr_f,tmin,tmax,dir_path+"/plots/deltam_cr_m_"+to_string(im)+"_rfw"+to_string(rfw)+"_P5P5.xmg");
	  deltam_cr=constant_fit(ma_fr_c,tmin,tmax,dir_path+"/plots/deltam_cr_m_"+to_string(im)+"_rfw"+to_string(rfw)+"_V0P5.xmg");
	  deltam_tm=0;
	}
    }
}

void perens_t::val_chir_extrap_deltam(perens_t &out) const
{
  //slice m
  vector<double> x(nm);
  djvec_t y(nm);
  for(size_t im=0;im<nm;im++)
    if(pars::chir_extr_method==chir_extr::MQUARK) x[im]=am[im];
    else                                          x[im]=sqr(meson_mass[im_im_ind({im,im})].ave());
  
  for(auto &t : vector<tuple<djvec_t*,const djvec_t*,const string>>{{&out.deltam_cr,&deltam_cr,"deltam_cr"},{&out.deltam_tm,&deltam_tm,"deltam_tm"}})
    {
      djvec_t &res=*get<0>(t);
      const djvec_t &in=*get<1>(t);
      const string &tag=get<2>(t);
      
      //open the plot file if needed
      const string plot_path=dir_path+"/plots/chir_extr_"+tag+".xmg";
      grace_file_t plot(plot_path);
      
      for(size_t r=0;r<nr;r++)
	{
	  //slice m
	  djvec_t y(nm);
	  for(size_t im=0;im<nm;im++)
	    y[im]=in[im_r_ind({im,r})];
	  
	  //fit, store and write the result
	  djvec_t coeffs=poly_fit(x,y,1);
	  res[out.im_r_ind({0,r})]=coeffs[0];
	  
	  auto xminmax=minmax_element(x.begin(),x.end());
	  double xmax=*xminmax.second*1.1;
	  write_fit_plot(plot,0,xmax,bind(poly_eval<djvec_t>,coeffs,_1),x,y);
	  plot.write_ave_err(0,coeffs[0].ave_err());
	}
    }
}
