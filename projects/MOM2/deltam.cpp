#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <MOM2/contractions.hpp>
#include <MOM2/deltam.hpp>
#include <MOM2/perens.hpp>
#include <MOM2/prop.hpp>
#include <MOM2/timings.hpp>

template <class T>
auto sigma_ansatz(const T &p,double p2,double p4_fr_p2=0.0)
  -> remove_reference_t<decltype(p[0])>
{
  return p[0]+p2*p[1]+p2*p2*p[2]+p4_fr_p2*p[3];
}

void perens_t::print_deltam(ostream &out) const
{
  for(size_t im=0;im<nm;im++)
    for(size_t r=0;r<nr;r++)
      {
	//get index
	const size_t imr=im_r_ind({im,r});
	
	//output
	out<<"m: "<<im<<", r: "<<r<<", deltam_tm: "<<deltam_tm[imr].ave_err()<<endl;
	out<<"m: "<<im<<", r: "<<r<<", deltam_cr: "<<deltam_cr[imr].ave_err()<<endl;
      }
}

void perens_t::compute_deltam_from_prop()
{
  const size_t mom_probe=linmoms.size()/8;
  const size_t iprobe=im_r_ilinmom_isigmaproj_isigmains_ind({0,0,mom_probe,sigma::SIGMA2,sigma::TM});
  cout<<"Probe sigma: "<<sigma[iprobe].ave_err()<<endl;
  const bool both=(fabs(sigma[iprobe].ave())>1e-10);
  if(not both) cout<<"Determining only Critical correction"<<endl;
  
  for(size_t im=0;im<nm;im++)
    for(size_t r=0;r<nr;r++)
      {
	//get index
	const size_t imr=im_r_ind({im,r});
	
	//allocate
	vector<double> x(linmoms.size());
	djvec_t deltam_tm_ct_corr(linmoms.size());
	djvec_t deltam_cr_ct_corr(linmoms.size());
	
	//degree and range of fit
	const size_t degree=3;
	const double p2_min=0.01,p2_max=2.0;
	
	for(size_t ilinmom=0;ilinmom<linmoms.size();ilinmom++)
	  {
	    //ascissa
	    x[ilinmom]=all_moms[linmoms[ilinmom][0]].p(L).norm2();
	    
	    auto get_sigma=[&](sigma::proj proj,sigma::ins ins) -> djack_t
	      {
		return sigma[im_r_ilinmom_isigmaproj_isigmains_ind({im,r,ilinmom,proj,ins})];
	      };
	    
	    //elements to solve the system:
	    //
	    // b x + c y - a = 0
	    // e x + f y - d = 0
	    using namespace sigma;
	    const djack_t& a=get_sigma(SIGMA2,PH);
	    const djack_t& b=get_sigma(SIGMA2,TM);
	    const djack_t& c=get_sigma(SIGMA2,CR);
	    const djack_t& d=get_sigma(SIGMA3,PH);
	    const djack_t& e=get_sigma(SIGMA3,TM);
	    const djack_t& f=get_sigma(SIGMA3,CR);
	    
	    //non singular case_of
	    if(both)
	      {
		const djack_t den=b*f-c*e;
		deltam_tm_ct_corr[ilinmom]=(-a*f+c*d)/den;
		deltam_cr_ct_corr[ilinmom]=(-b*d+a*e)/den;
	      }
	    else
	      {
		deltam_tm_ct_corr[ilinmom]=0.0;
		deltam_cr_ct_corr[ilinmom]=-a/c;
	      }
	  }
	
	//polynomial fit of the two deltas
	const djvec_t deltam_tm_ct_pars=poly_fit(x,deltam_tm_ct_corr,degree,p2_min,p2_max,
						 dir_path+"/plots/fit_deltam_tm_ct_m"+to_string(im)+"_r"+to_string(r)+".xmg");
	const djvec_t deltam_cr_ct_pars=poly_fit(x,deltam_cr_ct_corr,degree,p2_min,p2_max,
						 dir_path+"/plots/fit_deltam_cr_ct_m"+to_string(im)+"_r"+to_string(r)+".xmg");
	
	//get extrapolation and overwrite in singular case
	deltam_tm[imr]=deltam_tm_ct_pars[0];
	deltam_cr[imr]=deltam_cr_ct_pars[0];
	if(not both) deltam_tm[imr]=0.0;
      }
  
  print_deltam();
}

void perens_t::recompute_deltam()
{
  using namespace pars;
  
  if(use_QED==2) CRASH("Cannot compute delta without QED classic way");
  
  deltam_time.start();
  switch(deltam_method)
    {
    case FROM_PROP:
      compute_deltam_from_prop();
      break;
    case FROM_CORR:
      compute_deltam_from_corr();
      break;
    }
  deltam_time.stop();
}

void perens_t::bin_read_deltam()
{
  raw_file_t deltam_file(deltam_path(),"r");
  for(auto dtu : get_deltam_tasks())
    dtu.out->bin_read(deltam_file);
}

void perens_t::bin_write_deltam()
{
  raw_file_t deltam_file(deltam_path(),"w");
  for(auto dtu : get_deltam_tasks({this}))
    dtu.in[0]->bin_write(deltam_file);
}

perens_t& perens_t::get_deltam()
{
  //if file exists open it, otherwise compute it
  if(file_exists(deltam_path()))
    {
      cout<<"File "<<deltam_path()<<" found, opening"<<endl;
      bin_read_deltam();
    }
  else
    {
      cout<<"File "<<deltam_path()<<" not found, computing"<<endl;
      prepare_list_of_confs();
      
      recompute_deltam();
      
      bin_write_deltam();
    }
  
  //Prepare plot
  for(auto dtu : get_deltam_tasks())
    {
      djvec_t &deltam=*dtu.out;
      const string &tag=dtu.tag;
      
      grace_file_t plot(dir_path+"/plots/"+tag+"_vs_am.xmg");
      
      for(size_t r=0;r<nr;r++)
	{
	  plot.new_data_set();
	  for(size_t im=0;im<nm;im++) plot.write_ave_err(am[im],deltam[im_r_ind({im,r})].ave_err());
	}
    }
  
  //put to zero if not asked to use
  if(not pars::use_deltam_cr_ct) deltam_cr=0.0;
  if(not pars::use_deltam_tm_ct) deltam_tm=0.0;
  
  //Print
  for(size_t im=0;im<nm;im++)
    for(size_t r=0;r<nr;r++)
      {
	cout<<"Deltam[m="<<im<<",r="<<r<<"]\t";
	for(auto dtu : get_deltam_tasks())
	  {
	    djvec_t &deltam=*dtu.out;
	    const string &tag=dtu.tag;
	    
	    cout<<" "<<tag<<": "<<deltam[im_r_ind({im,r})];
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
	  mPCAC[im]=compute_mPCAC(im);
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
	mPCAC_sea=compute_mPCAC(nm);
	mPCAC_sea.bin_write(mPCAC_sea_path);
      }
  mPCAC2_plot.write_ave_err(0,mPCAC_sea.ave_err());
  cout<<"Sea PCAC: "<<smart_print(mPCAC_sea.ave_err())<<endl;
  
  return *this;
}

void perens_t::compute_deltam_from_corr()
{
  prepare_list_of_confs();
  
  for(size_t im=0;im<nm;im++)
    for(size_t rfw=0;rfw<nr;rfw++)
      {
	//ext_reim and rpar is relative to non-inserted
	
	//measure mcrit according to eq.3 of hep-lat/0701012
	{
	  using namespace qprop;
	  
	  const djvec_t P5P5_00=get_contraction(im,LO,im,LO,"P5P5",RE,EVN,rfw,0);
	  const djvec_t V0P5_00=get_contraction(im,LO,im,LO,"V0P5",IM,ODD,rfw,0);
	  const djvec_t m_cr_corr=2.0*forward_derivative(V0P5_00)/(2.0*P5P5_00);
	  const djvec_t m_cr_corr_symm=2.0*symmetric_derivative(V0P5_00)/(2.0*P5P5_00);
	  const djack_t m_cr=constant_fit(m_cr_corr,tmin,tmax,dir_path+"/plots/m_cr_"+to_string(im)+".xmg");
	  const djack_t m_cr_symm=constant_fit(m_cr_corr_symm,tmin,tmax,dir_path+"/plots/m_cr_symm_"+to_string(im)+".xmg");
	  cout<<"m_cr[m="<<im<<",rfw="<<rfw<<"]: "<<smart_print(m_cr.ave_err())<<", symm: "<<smart_print(m_cr_symm.ave_err())<<endl;
	}
	
	using namespace qprop;
	
	int rdiff=0;
	
	switch(pars::use_QED)
	  {
	  case 0:
	    break;
	  case 1:
	    {
	      //load corrections
	      const djvec_t P5P5_00=get_contraction(im,LO,im,LO,"P5P5",RE,UNK,rfw,rdiff);
	      const djvec_t P5P5_LL=get_contraction(im,F,im,F,"P5P5",RE,UNK,rfw,rdiff);
	      const djvec_t P5P5_0M=get_contraction(im,LO,im,FF,"P5P5",RE,UNK,rfw,rdiff);
	      const djvec_t P5P5_M0=get_contraction(im,FF,im,LO,"P5P5",RE,UNK,rfw,rdiff);
	      const djvec_t P5P5_0T=get_contraction(im,LO,im,T,"P5P5",RE,UNK,rfw,rdiff);
	      const djvec_t P5P5_T0=get_contraction(im,T,im,LO,"P5P5",RE,UNK,rfw,rdiff);
	      //load the derivative wrt counterterm
	      const djvec_t P5P5_0P=get_contraction(im,LO,im,P,"P5P5",RE,UNK,rfw,rdiff);
	      const djvec_t P5P5_P0=get_contraction(im,P,im,LO,"P5P5",RE,UNK,rfw,rdiff);
	      //load the derivative wrt mass
	      const djvec_t P5P5_0S=get_contraction(im,LO,im,S,"P5P5",RE,UNK,rfw,rdiff);
	      const djvec_t P5P5_S0=get_contraction(im,S,im,LO,"P5P5",RE,UNK,rfw,rdiff);
	      
	      //load corrections
	      const djvec_t V0P5_00=get_contraction(im,LO,im,LO,"V0P5",IM,UNK,rfw,rdiff);
	      const djvec_t V0P5_LL=get_contraction(im,F,im,F,"V0P5",IM,UNK,rfw,rdiff);
	      const djvec_t V0P5_0M=get_contraction(im,LO,im,FF,"V0P5",IM,UNK,rfw,rdiff);
	      const djvec_t V0P5_M0=get_contraction(im,FF,im,LO,"V0P5",IM,UNK,rfw,rdiff);
	      const djvec_t V0P5_0T=get_contraction(im,LO,im,T,"V0P5",IM,UNK,rfw,rdiff);
	      const djvec_t V0P5_T0=get_contraction(im,T,im,LO,"V0P5",IM,UNK,rfw,rdiff);
	      //load the derivative wrt counterterm
	      const djvec_t V0P5_0P=get_contraction(im,LO,im,P,"V0P5",IM,UNK,rfw,rdiff);
	      const djvec_t V0P5_P0=get_contraction(im,P,im,LO,"V0P5",IM,UNK,rfw,rdiff);
	      //load the derivative wrt mass
	      const djvec_t V0P5_0S=get_contraction(im,LO,im,S,"V0P5",IM,UNK,rfw,rdiff);
	      const djvec_t V0P5_S0=get_contraction(im,S,im,LO,"V0P5",IM,UNK,rfw,rdiff);
	      
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
	      
	      const size_t imr=im_r_ind({im,rfw});
	      
	      if(strategy==both_cr_and_tm)
		{
		  const djvec_t old_corr=djvec_t(symmetric_derivative(V0P5_QED)/symmetric_derivative(V0P5_P)).symmetrized();
		  const djack_t old=constant_fit(old_corr,tmin,tmax,dir_path+"/plots/old_deltam_cr_m"+to_string(im)+"_rfw"+to_string(rfw)+".xmg");
		  
		  const djvec_t den=b*f-c*e;
		  const djvec_t deltam_tm_corr=(-a*f+c*d)/den;
		  const djvec_t deltam_cr_corr=(-b*d+a*e)/den;
		  
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
		  
		  djvec_t a=V0P5_QED.symmetrized(-1);
		  djvec_t c=V0P5_P.symmetrized(-1);
		  djvec_t d=P5P5_QED.symmetrized();
		  djvec_t f=P5P5_P.symmetrized();
		  
		  a.ave_err().write(dir_path+"/plots/deltam_a_m"+to_string(im)+"_rfw"+to_string(rfw)+".xmg");
		  c.ave_err().write(dir_path+"/plots/deltam_c_m"+to_string(im)+"_rfw"+to_string(rfw)+".xmg");
		  
		  const djvec_t ma_fr_c=-a/c;
		  const djvec_t md_fr_f=-d/f;
		  
		  deltam_cr[imr]=constant_fit(md_fr_f,tmin,tmax,dir_path+"/plots/deltam_cr_m"+to_string(im)+"_rfw"+to_string(rfw)+"_P5P5.xmg");
		  deltam_cr[imr]=constant_fit(ma_fr_c,tmin,tmax,dir_path+"/plots/deltam_cr_m"+to_string(im)+"_rfw"+to_string(rfw)+"_V0P5.xmg");
		  // deltam_cr[imr]=-2*0.16285705871085078618;
		  deltam_tm[imr]=0;
		}
	    }
	    break;
	  case 2:
	    CRASH("Cannot be implemented");
	    break;
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
  
  for(auto t : out.get_deltam_tasks({this}))
    {
      djvec_t &res=*t.out;
      const djvec_t &in=*t.in[0];
      const string &tag=t.tag;;
      
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
	  
	  cout<<tag<<" "<<out.im_r_ind({0,r})<<" "<<&res[out.im_r_ind({0,r})]<<" poss1: "<<&out.deltam_cr[out.im_r_ind({0,r})]<<" poss2: "<<&out.deltam_tm[out.im_r_ind({0,r})]<<endl;
	}
    }
  
  out.print_deltam();
}

vector<perens_t::task_t> perens_t::get_deltam_tasks(const vector<const perens_t*> &ens)
{
  vector<const djvec_t*> in_deltam_tm;
  vector<const djvec_t*> in_deltam_cr;
  
  for(auto &e : ens)
    {
      in_deltam_tm.push_back(&e->deltam_tm);
      in_deltam_cr.push_back(&e->deltam_cr);
    }
  
  vector<task_t> deltam_tasks=
    {{&deltam_tm,in_deltam_tm,im_r_ind,"deltam_tm",QED_task},
     {&deltam_cr,in_deltam_cr,im_r_ind,"deltam_cr",QED_task}};
  
  return deltam_tasks;
}

void perens_t::average_r_deltam(perens_t &out) const
{
  cout<<"Averaging r for deltam"<<endl;
  
  if(nr>1)
    for(size_t im=0;im<nm;im++)
      {
	out.deltam_cr[im]=out.deltam_tm[im]=0.0;
	for(size_t r=0;r<nr;r++)
	  {
	    out.deltam_cr[im]+=deltam_cr[im_r_ind({im,r})];
	    out.deltam_tm[im]+=deltam_tm[im_r_ind({im,r})];
	  }
	out.deltam_cr[im]/=nr;
	out.deltam_tm[im]/=nr;
      }
}
