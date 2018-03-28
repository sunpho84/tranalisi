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
  deltam_tm.resize(im_r_ind.max());
  deltam_cr.resize(im_r_ind.max());
  
  vector<tuple<djvec_t*,string,bool>> delta_tasks{
    {&deltam_cr,"cr",pars::use_deltam_cr_ct},
    {&deltam_tm,"tm",pars::use_deltam_tm_ct}};
  
  for(auto &dtu : delta_tasks)
    {
      djvec_t &deltam=*get<0>(dtu);
      const string &tag=get<1>(dtu);
      
      //if file exists open it, otherwise compute it
      const string deltam_path=dir_path+"/deltam_"+tag+".dat";
      if(file_exists(deltam_path))
	{
	  cout<<"File "<<deltam_path<<" found, opening"<<endl;
	  deltam.bin_read(deltam_path);
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
	  deltam.bin_write(deltam_path);
	}
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
    const djvec_t m_cr_corr=forward_derivative(V0P5_00)/(2.0*P5P5_00);
    const djvec_t m_cr_corr_symm=symmetric_derivative(V0P5_00)/(2.0*P5P5_00);
    const djack_t m_cr=constant_fit(m_cr_corr,tmin,tmax,dir_path+"/plots/m_cr_"+to_string(im)+".xmg");
    const djack_t m_cr_symm=constant_fit(m_cr_corr_symm,tmin,tmax,dir_path+"/plots/m_cr_symm_"+to_string(im)+".xmg");
    cout<<"m_cr[m="<<im<<",rfw="<<rfw<<"]: "<<m_cr.ave_err()<<", symm: "<<m_cr_symm.ave_err()<<endl;
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
	a=
	V0P5_0T+
	V0P5_0M+
	V0P5_T0+
	V0P5_M0+
	V0P5_LL,
	b=
	V0P5_0S+
	V0P5_S0,
	c=
	V0P5_0P+
	V0P5_P0;
      const djvec_t
	d=
	P5P5_0T+
	P5P5_0M+
	P5P5_T0+
	P5P5_M0+
	P5P5_LL,
	e=
	P5P5_0S+
	P5P5_S0,
	f=
	P5P5_0P+
	P5P5_P0;
      
      djvec_t del_t=symmetric_derivative(a)/(2.0*P5P5_00)-symmetric_derivative(V0P5_00)/(2.0*P5P5_00*P5P5_00)*d;
      djack_t del=constant_fit(del_t,tmin,tmax,dir_path+"/plots/test_deltam_cr_m"+to_string(im)+"_rfw"+to_string(rfw)+".xmg");
      
      const djvec_t den=b*f-c*e;
      const djvec_t deltam_tm_corr=djvec_t((-a*f+c*d)/den).symmetrized();
      const djvec_t deltam_cr_corr=djvec_t((-b*d+a*e)/den).symmetrized();
      
      const size_t imr=im_r_ind({im,rfw});
      
      for(auto &oit : vector<tuple<djvec_t*,const djvec_t*,string>>{{&deltam_cr,&deltam_cr_corr,"cr"},{&deltam_tm,&deltam_tm_corr,"tm"}})
	{
	  djack_t &out=(*::get<0>(oit))[imr];
	  const djvec_t &in=*::get<1>(oit);
	  const string tag=::get<2>(oit);
	  
	  out=constant_fit(in,tmin,tmax,dir_path+"/plots/deltam_"+tag+"_m"+to_string(im)+"_rfw"+to_string(rfw)+".xmg");
	}
    }
}
