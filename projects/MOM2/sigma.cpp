#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_SIGMA
 #include <MOM2/sigma.hpp>

#include <analysis.hpp>
#include <corrections.hpp>
#include <perens.hpp>
#include <timings.hpp>

namespace sigma
{
  void set_ins()
  {
    switch(pars::use_QED)
      {
      case 0:
	ins_list={LO};
	break;
      case 1:
	ins_list={LO , QED , CR , TM , PH};
	break;
      case 2:
	ins_list={LO , QED};
	break;
     }
    
    //set sigma proj
    proj_list={SIGMA1,SIGMA2,SIGMA3};
    if(pars::compute_ri)
      proj_list.push_back(SIGMA4);
    nproj=proj_list.size();
    cout<<"Sigma, nproj: "<<nproj<<endl;
    
    //set sigma ins
    iins_of_ins.resize(ins_tag.size());
    for(size_t iins=0;iins<ins_list.size();iins++)
      iins_of_ins[ins_list[iins]]=iins;
    nins=ins_list.size();
    cout<<"Sigma, nins: "<<nins<<endl;
  }
}

void perens_t::plot_sigma(const string &suffix)
{
  cout<<"Plotting all sigma of "<<dir_path<<" for suffix: \""<<suffix<<"\""<<endl;
  
  for(size_t iins=0;iins<sigma::nins;iins++)
    for(size_t iproj=0;iproj<sigma::nproj;iproj++)
      {
	const sigma::ins ins=sigma::ins_list[iins];
	grace_file_t out(dir_path+"/plots/sigma"+to_string(iproj+1)+"_"+sigma::ins_tag[ins]+(suffix!=""?("_"+suffix):string(""))+".xmg");
	
	for(size_t im=0;im<nm;im++)
	  for(size_t r=0;r<nr;r++)
	    {
	      out.new_data_set();
	      for(size_t imom=0;imom<linmoms.size();imom++)
		{
		  const double p2tilde=all_moms[linmoms[imom][0]].p(L).tilde().norm2();
		  out.write_ave_err(p2tilde,sigma[im_r_ilinmom_isigmaproj_isigmains_ind({im,r,imom,iproj,iins})].ave_err());
	      }
	  }
    }
}

void perens_t::subtract_Oa2_sigma()
{
  const size_t iins=0;
  const size_t iproj=sigma::SIGMA1;
  for(size_t im=0;im<nm;im++)
    for(size_t r=0;r<nr;r++)
      for(size_t ilinmom=0;ilinmom<linmoms.size();ilinmom++)
	{
	  const size_t imom=linmoms[ilinmom][0];
	  const imom_t mom=all_moms[imom];
	  const double sub=g2to_correct()*sig1_a2(pars::act,gf::LANDAU,group::SU3,mom,L);
	  
	  sigma[im_r_ilinmom_isigmaproj_isigmains_ind({im,r,imom,iproj,iins})]-=sub;
	}
}

void perens_t::evolve_sigma(perens_t &out) const
{
  CRASH("");
}

vector<perens_t::task_t> perens_t::get_sigma_tasks(const vector<const perens_t*> &ens)
{
  vector<const djvec_t*> in_sigma;
  
  for(auto &e : ens)
    in_sigma.push_back(&e->sigma);
  
  vector<task_t> sigma_tasks={{&sigma,in_sigma,im_r_ilinmom_isigmaproj_isigmains_ind,"sigma",QCD_task}};
  
  return sigma_tasks;
}

perens_t& perens_t::compute_sigmas()
{
  //! list of all combination of transformations to be applied
  vector<pair<size_t,size_t>> map;
  
#define ADD_COMBO(MAP,INS1,INS2) MAP.push_back({sigma::iins_of_ins[sigma::INS1],jqprop::iins_of_ins[jqprop::INS2]})
#define ADD_COMBO_SAME(MAP,INS) ADD_COMBO(MAP,INS,INS)
  
  ADD_COMBO_SAME(map,LO);
  switch(pars::use_QED)
    {
    case 0:
      break;
    case 1:
      ADD_COMBO_SAME(map,PH);
      ADD_COMBO_SAME(map,CR);
      ADD_COMBO_SAME(map,TM);
      break;
    case 2:
      ADD_COMBO_SAME(map,QED);
      break;
    }
#undef ADD_COMBO_SAME
  
  //map for RI
  vector<pair<size_t,size_t>> map_ri;
  if(pars::compute_ri)
    {
      ADD_COMBO(map_ri,LO,RI);
      switch(pars::use_QED)
	{
	case 0:
	  break;
	case 1:
	  CRASH("Not implemented yet");
	  break;
	case 2:
	  ADD_COMBO(map_ri,QED,RI_QED);
	  break;
	}
    }
#undef ADD_COMBO
  
  vector<raw_file_t> files=setup_read_all_qprops_mom(conf_list);
  
  for(size_t ilinmom=0;ilinmom<linmoms.size();ilinmom++)
    {
      const size_t mom=linmoms[ilinmom][0];
      vector<jqprop_t> jprops(im_r_ijqins_ind.max()); //!< jackknived props
      
      for(size_t i_in_clust=0;i_in_clust<clust_size;i_in_clust++)
	for(size_t ihit=0;ihit<nhits_to_use;ihit++)
	  {
	    const size_t i_in_clust_ihit=i_in_clust_ihit_ind({i_in_clust,ihit});
	    const size_t mom=linmoms[ilinmom][0];
	    cout<<"Reading qprop, "
	      "clust_entry "<<i_in_clust+1<<"/"<<clust_size<<", "
	      "hit "<<ihit+1<<"/"<<nhits<<", "
	      "momentum "<<ilinmom+1<<"/"<<linmoms.size()<<", "
	      "mom: "<<mom<<endl;
	    read_time.start();
	    const vector<qprop_t> props=read_all_qprops_mom(files,i_in_clust_ihit,mom);
	    read_time.stop();
	    
	    //build all props
	    build_props_time.start();
	    build_all_mr_jackkniffed_qprops(jprops,props);
	    build_props_time.stop();
	  }
      
      //clusterize
      clust_time.start();
      for(auto &j : jprops) j.clusterize(clust_size);
      clust_time.stop();
      
      const vector<jqprop_t> jprops_inv=get_inverse_propagators(jprops);
      
#pragma omp parallel for reduction(+:sigma_time)
      for(size_t im_r_ijack=0;im_r_ijack<im_r_ijackp1_ind.max();im_r_ijack++)
	{
	  //decript indices
	  const vector<size_t> im_r_ijack_comps=im_r_ijackp1_ind(im_r_ijack);
	  const size_t im=im_r_ijack_comps[0],r=im_r_ijack_comps[1],ijack=im_r_ijack_comps[2];
	  
	  //cout<<" Computing projected propagator for combo im="<<im<<" , r="<<r<<" , ijack="<<ijack<<endl;
	  
	  const p_t ptilde=all_moms[mom].p(L).tilde();
	  
	  //! function to compute sigma 1,2,3
	  auto compute_sigma=[&](const qprop_t prop_inv,sigma::proj proj)
	    {
	      using namespace sigma;
	      
	      double out=0.0;
	      
	      switch(proj)
		{
		case SIGMA1:
		  //trace with pslash
		  out=0.0;
		  for(size_t mu=0;mu<NDIM;mu++)
		    if(fabs(ptilde[mu])>1e-10)
		      out+=
			(prop_inv*quaGamma[igmu[mu]]).trace().imag()/
			(12.0*ptilde[mu]*V*all_moms[mom].Np());
		  break;
		case SIGMA2:
		  //trace with identity
		  out=(prop_inv*quaGamma[0]).trace().real()/(12.0*V);
		  break;
		case SIGMA3:
		  //trace with gamma5
		  out=(prop_inv*quaGamma[5]).trace().imag()/(12.0*V);
		case SIGMA4:
		  //trace with gamma_mu
		  out=0.0;
		  for(size_t mu=0;mu<NDIM;mu++)
		    if(mu==1)
		      out+=
			(prop_inv*quaGamma[igmu[mu]]).trace().real()/(12.0*V);
		  break;
		};
	      
	      return out;
	    };
	  
	  sigma_time.start();
	  for(size_t iproj=0;iproj<sigma::nproj;iproj++)
	    {
	      const sigma::proj proj=sigma::proj_list[iproj];
	      //cout<<"   Computing proj "<<iproj<<"/"<<sigma::nproj<<endl;
	      
	      
	      for(auto m : (proj==sigma::SIGMA4 ? map_ri : map))
		{
		  const size_t isins=get<0>(m);
		  const size_t ijqins=get<1>(m);
		  
		  const size_t im_r_ilinmom_isigmaproj_isigmains=im_r_ilinmom_isigmaproj_isigmains_ind({im,r,ilinmom,iproj,isins});
		  const size_t im_r_ijqins=im_r_ijqins_ind({im,r,ijqins});
		  sigma[im_r_ilinmom_isigmaproj_isigmains][ijack]=compute_sigma(jprops_inv[im_r_ijqins][ijack],proj);
		}
	    }
	  sigma_time.stop();
	  
#undef COMPUTE_SIGMA
	}
    }
  
  return *this;
}

void perens_t::average_r_sigma(perens_t &out) const
{
  cout<<"Averaging r for sigma"<<endl;
  
  djvec_t &sigma_rave=out.sigma;
  
#pragma omp parallel for
  for(size_t out_i=0;out_i<out.im_r_ilinmom_isigmaproj_isigmains_ind.max();out_i++)
    {
      const vector<size_t> out_comps=out.im_r_ilinmom_isigmaproj_isigmains_ind(out_i);
      vector<size_t> in_comps=out_comps;
      const sigma::proj isigmaproj=sigma::proj_list[in_comps[3]];
      
      sigma_rave[out_i]=0.0;
      for(size_t r=0;r<nr;r++)
	{
	  in_comps[1]=r;
	  const size_t in_i=im_r_ilinmom_isigmaproj_isigmains_ind(in_comps);
	  
	  //include a -1 on SIGMA3 and second r
	  const int coeff=(isigmaproj==sigma::SIGMA3 and r==1)?-1:+1;
	  
	  sigma_rave[out_i]+=sigma[in_i]*coeff;
	}
      sigma_rave[out_i]/=nr;
    }
}

void perens_t::average_equiv_momenta_sigma(perens_t &out,const vector<vector<size_t>> &equiv_linmom_combos) const
{
#pragma omp parallel for
  for(size_t iout=0;iout<out.im_r_ilinmom_isigmaproj_isigmains_ind.max();iout++)
    {
      const vector<size_t> out_comps=out.im_r_ilinmom_isigmaproj_isigmains_ind(iout);
      const size_t out_ilinmom_combo=out_comps[2];
      
      djack_t &ave=out.sigma[iout];
      ave=0.0;
      for(const size_t ieq : equiv_linmom_combos[out_ilinmom_combo])
	{
	  vector<size_t> in_comps=out_comps;
	  in_comps[2]=ieq;
	  const size_t iin=im_r_ilinmom_isigmaproj_isigmains_ind(in_comps);
	  
	  ave+=sigma[iin];
	}
      ave/=equiv_linmom_combos[out_ilinmom_combo].size();
    }
}

void perens_t::val_chir_extrap_sigma(perens_t &out) const
{
  index_t r_ilinmom_isigmaproj_ind({{"r",nr},{"linmom",linmoms.size()},{"sigmaproj",sigma::nproj}});
  index_t r_ilinmom_isigmaproj_isigmains_ind({{"r",nr},{"linmom",linmoms.size()},{"sigmaproj",sigma::nproj},{"sigmains",sigma::nins}});
  
  if(pars::use_QED) needs_to_read_assembled_QED_greenfunctions();
  
  for(auto &t : out.get_sigma_tasks({this}))
    {
      const djvec_t &sigma=*t.in.front();
      djvec_t &sigma_chir=*t.out;
      const string &tag=t.tag;
      
#pragma omp parallel for
      for(size_t r_ilinmom_isigmaproj=0;r_ilinmom_isigmaproj<r_ilinmom_isigmaproj_ind.max();r_ilinmom_isigmaproj++)
	{
	  const vector<size_t> r_ilinmom_isigmaproj_comps=r_ilinmom_isigmaproj_ind(r_ilinmom_isigmaproj);
	  const size_t r=r_ilinmom_isigmaproj_comps[0];
	  const size_t ilinmom=r_ilinmom_isigmaproj_comps[1];
	  const size_t isigmaproj=r_ilinmom_isigmaproj_comps[2];
	  
	  vector<djvec_t> coeffs(sigma::nins);
	  
	  for(size_t isigmains=0;isigmains<sigma::nins;isigmains++)
	  {
	    const vector<size_t> r_ilinmom_isigmaproj_isigmains_comps={r,ilinmom,isigmaproj,isigmains};
	    const size_t r_ilinmom_isigmaproj_isigmains=r_ilinmom_isigmaproj_isigmains_ind({r_ilinmom_isigmaproj_isigmains_comps});
	    
	    //open the plot file if needed
	    const string plot_path=dir_path+"/plots/chir_extr_"+tag+"_"+r_ilinmom_isigmaproj_isigmains_ind.descr(r_ilinmom_isigmaproj_isigmains)+".xmg";
	    grace_file_t *plot=nullptr;
	    if(ilinmom%pars::print_each_mom==0) plot=new grace_file_t(plot_path);
	    
	    //slice m and fit iter_on_n3
	    vector<double> x(nm);
	    djvec_t y(nm);
	    for(size_t im=0;im<nm;im++)
	      {
		const size_t imeson=im_im_ind({im,im});
		
		//compute mass sum
		if(pars::chir_extr_method==chir_extr::MQUARK) x[im]=2.0*am[im];
		else                                          x[im]=sqr(meson_mass[imeson].ave());
		
		y[im]=sigma[im_r_ilinmom_isigmaproj_isigmains_ind(concat(im,r_ilinmom_isigmaproj_isigmains_comps))];

		//if QED case take into account variation due to leading pole
		if(pars::use_QED and isigmains>0 and pars::sub_meson_mass_shift_when_no_pole)
		  {
		    const djack_t M=meson_mass[imeson],dM=meson_mass_QED[imeson];

		    const djack_t b0=coeffs[sigma::LO][1];
		    const djack_t varb=2.0*b0*dM*M;
		    y[im]-=varb;
		  }
	      }
	    
	    //fit, store and write the result
	    coeffs[isigmains]=poly_fit(x,y,1);
	    sigma_chir[out.im_r_ilinmom_isigmaproj_isigmains_ind(concat((size_t)0,r_ilinmom_isigmaproj_isigmains_comps))]=coeffs[isigmains][0];
	    if(plot!=nullptr)
	      {
		auto xminmax=minmax_element(x.begin(),x.end());
		double xmax=*xminmax.second*1.1;
		write_fit_plot(*plot,0,xmax,bind(poly_eval<djvec_t>,coeffs[isigmains],_1),x,y);
		plot->write_ave_err(0,coeffs[isigmains][0].ave_err());
	      }
	    
	  if(plot!=nullptr) delete plot;
	  }
	}
    }
}

void perens_t::assemble_sigma_QED_greenfunctions()
{
  cout<<"Assembling sigma QED greenfunctions"<<endl;
  
#pragma omp parallel for
  for(size_t im_r_ilinmom=0;im_r_ilinmom<im_r_ilinmom_ind.max();im_r_ilinmom++)
    {
      const vector<size_t> comps=im_r_ilinmom_ind(im_r_ilinmom);
      
      const size_t im=comps[0];
      const size_t r=comps[1];
      const size_t ilinmom=comps[2];
      const size_t im_r=im_r_ind({im,r});
     
      using namespace sigma;
      auto sigma1=sigma_ins_getter(im,r,ilinmom,SIGMA1);
      
      switch(pars::use_QED)
	{
	case 0:
	  break;
	case 1:
	  sigma1(QED)=
	    sigma1(PH)+
	    sigma1(CR)*deltam_cr[im_r]+
	    sigma1(TM)*deltam_tm[im_r];
	  break;
	case 2:
	  break;
	}
    }
}
