#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <MOM2/contractions.hpp>
#include <MOM2/deltam_cr.hpp>
#include <MOM2/perens.hpp>
#include <MOM2/timings.hpp>

perens_t& perens_t::get_deltam_cr()
{
  deltam_cr.resize(nm);
  
  //if file exists open it, otherwise compute it
  const string deltam_cr_path=dir_path+"/deltam_cr.dat";
  if(file_exists(deltam_cr_path))
    {
      cout<<"File "<<deltam_cr_path<<" found, opening"<<endl;
      deltam_cr.bin_read(deltam_cr_path);
    }
  else
    {
      cout<<"File "<<deltam_cr_path<<" not found, computing"<<endl;
      prepare_list_of_confs();
      
      deltam_cr_time.start();
      for(size_t im=0;im<nm;im++)
	deltam_cr[im]=compute_deltam_cr(im);
      deltam_cr_time.stop();
      deltam_cr.bin_write(deltam_cr_path);
    }
  
  //Print, putting it to zero if not asked to use
  for(size_t im=0;im<nm;im++)
    {
      cout<<"Deltam cr["<<im<<"]: "<<deltam_cr[im];
      if(not pars::use_deltam_cr_ct)
	{
	  deltam_cr[im]=0.0;
	  cout<<" -> "<<0.0;
	}
      cout<<endl;
    }
  
  return *this;
}

djack_t perens_t::compute_deltam_cr(const size_t im)
{
  auto get=[im,this]
    (string tag_bw,string tag_fw,const string &ID,const size_t reim,const int tpar,const int rpar,const int rdiff)
    {
      djvec_t res(L[0]/2+1);
      res=0.0;
      for(size_t r=0;r<((rdiff==0)?nr:1);r++) //looping on rfw
	{
	  const size_t rfw=r;
	  const size_t rbw=(r+rdiff)%nr;
	  const string name="M"+to_string(im)+"_R"+to_string(rbw)+"_"+tag_bw+"_M"+to_string(im)+"_R"+to_string(rfw)+"_"+tag_fw;
	  const djvec_t contr=get_contraction(name,ID,reim,tpar)*((r==0)?1:rpar)/(1+abs(rpar));
	  contr.ave_err().write(dir_path+"/plots/"+ID+"_"+name+".xmg");
	  res+=contr; //Warning, might be ill normalized if nr==1
	}
      
      return res;
    };
  
  //measure mcrit according to eq.3 of hep-lat/0701012
  {
    const djvec_t P5P5_00=get("0","0","P5P5",RE,EVN,EVN,0);
    const djvec_t V0P5_00=get("0","0","V0P5",IM,ODD,ODD,0);
    const djvec_t m_cr_corr=forward_derivative(V0P5_00)/(2.0*P5P5_00);
    const djvec_t m_cr_corr_symm=(forward_derivative(V0P5_00)+backward_derivative(V0P5_00))/(4.0*P5P5_00);
    const djack_t m_cr=constant_fit(m_cr_corr,tmin,tmax,dir_path+"/plots/m_cr_"+to_string(im)+".xmg");
    const djack_t m_cr_symm=constant_fit(m_cr_corr_symm,tmin,tmax,dir_path+"/plots/m_cr_symm_"+to_string(im)+".xmg");
    cout<<"m_cr["<<im<<"]: "<<m_cr.ave_err()<<", symm: "<<m_cr_symm.ave_err()<<endl;
  }
  
  if(pars::use_QED)
    {
      int rdiff=pars::use_deltam_cr_ropp;
      if(pars::use_deltam_cr_ropp and nr!=2) CRASH("Cannot use opposite r when computing deltam_cr if nr is not 2");
      
      const djvec_t P5P5_00=get("0","0","P5P5",RE,EVN,EVN,rdiff);
      //load corrections
      //const djvec_t V0P5_LL=get("F","F","V0P5",IM,ODD,ODD,rdiff);
    const djvec_t V0P5_00=get("0","0","V0P5",IM,ODD,ODD,rdiff);
      const djvec_t V0P5_0M=get("0","FF","V0P5",IM,ODD,ODD,rdiff);
      //const djvec_t V0P5_M0=get("FF","0","V0P5",IM,ODD,ODD,rdiff);
      const djvec_t V0P5_0T=get("0","T","V0P5",IM,ODD,ODD,rdiff);
      const djvec_t V0P5_0S=get("0","S","V0P5",IM,ODD,ODD,rdiff);
      //const djvec_t V0P5_T0=get("T","0","V0P5",IM,ODD,ODD,rdiff);
      //load the derivative wrt counterterm
      const djvec_t V0P5_0P=get("0","P","V0P5",RE,ODD,EVN,rdiff);
      //const djvec_t V0P5_P0=get("P","0","V0P5",RE,ODD,EVN,rdiff);
      
      //build numerator
      const djvec_t num_deltam_cr_corr=
	//+V0P5_LL
	+V0P5_0M
	//+V0P5_M0
	+V0P5_0T
	//+V0P5_T0
	;
      //build denominator
      const djvec_t den_deltam_cr_corr=
	//+V0P5_P0
	-V0P5_0P;
      const djvec_t deltam_cr_corr=num_deltam_cr_corr/den_deltam_cr_corr;
      const djack_t deltam_cr=constant_fit(deltam_cr_corr,tmin,tmax,dir_path+"/plots/deltam_cr.xmg");
      
      return deltam_cr;
    }
  else return djack_t{};
}
