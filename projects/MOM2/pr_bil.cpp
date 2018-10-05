#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_PR_BIL
 #include <MOM2/pr_bil.hpp>

#include <MOM2/corrections.hpp>
#include <MOM2/perens.hpp>
#include <MOM2/prop.hpp>
#include <MOM2/timings.hpp>

namespace pr_bil
{
  void set_ins()
  {
    switch(pars::use_QED)
      {
      case 0:
	ins_list={LO};
	break;
      case 1:
	ins_list={LO , QED , CR_OU , CR_IN , TM_OU , TM_IN , PH_OU , PH_IN , EX };
	break;
      case 2:
	ins_list={LO , QED , QED_OU , QED_IN, EX};
	break;
     }
    
    iins_of_ins.resize(ins_tag.size());
    for(size_t iins=0;iins<ins_list.size();iins++)
      iins_of_ins[ins_list[iins]]=iins;
    nins=ins_list.size();
    cout<<"Pr_bil, nins: "<<nins<<endl;
  }
}

void perens_t::build_all_mr_gbil_jackkniffed_verts(vector<jqprop_t>& jbil,const vector<qprop_t>& props_in,const vector<qprop_t>& props_ou) const
{
  //! help finding the bilinear/clust combo
  index_t ind({{"i",im_r_im_r_igam_ind.max()},{"iclust",njacks}});
  
  //! list of all combination of transformations to be applied
  vector<tuple<size_t,size_t,size_t>> map;
#define ADD_COMBO(V,O,I)			\
  map.push_back({pr_bil::iins_of_ins[pr_bil::V],qprop::iins_of_ins[qprop::O],qprop::iins_of_ins[qprop::I]})
  ADD_COMBO(LO,    LO, LO);
  switch(pars::use_QED)
    {
    case 0:
      break;
    case 1:
      ADD_COMBO(EX,    F,  F);
      ADD_COMBO(PH_IN, FF, LO);
      ADD_COMBO(PH_IN, T,  LO);
      ADD_COMBO(PH_OU, LO, FF);
      ADD_COMBO(PH_OU, LO, T);
      ADD_COMBO(CR_IN, P,  LO);
      ADD_COMBO(CR_OU, LO, P);
      ADD_COMBO(TM_IN, S,  LO);
      ADD_COMBO(TM_OU, LO, S);
      break;
    case 2:
      ADD_COMBO(EX,     F,  F);
      ADD_COMBO(QED_IN, QED,LO);
      ADD_COMBO(QED_OU, LO, QED);
      break;
    }
#undef ADD_COMBO
  //cout<<"pr_bil, needs to compute "<<map.size()<<" kind of vertices"<<endl;
  
#pragma omp parallel for
  for(size_t i=0;i<ind.max();i++)
    {
      //decript bilinar/jack
      vector<size_t> ind_comp=ind(i);
      const size_t im_r_im_r_igam=ind_comp[0],iclust=ind_comp[1];
      
      //decript props to combine and which Gamma to create
      const vector<size_t> im_r_im_r_igam_comp=im_r_im_r_igam_ind(im_r_im_r_igam);
      const size_t im_ou=im_r_im_r_igam_comp[0],r_ou=im_r_im_r_igam_comp[1];
      const size_t im_in=im_r_im_r_igam_comp[2],r_in=im_r_im_r_igam_comp[3];
      const size_t iG=im_r_im_r_igam_comp[4];
      
      for(auto t : map)
	{
	  const size_t ipr_ins=get<0>(t);
	  const size_t iq_ins_in=get<1>(t);
	  const size_t iq_ins_ou=get<2>(t);
	  
	  //cout<<" pr_bil , im_ou="<<im_ou<<" , r_ou="<<r_ou<< ", im_in="<<im_in<<" , r_in="<<r_in<<" iG="<<iG<<" , ipr_ins="<<ipr_ins<<" , iq_ins_in="<<iq_ins_in<<" , iq_ins_ou="<<iq_ins_ou<<endl;
	  //proxy for vector and props
	  const qprop_t &p_in=props_in[im_r_iqins_ijack_ind({im_in,r_in,iq_ins_in,iclust})];
	  const qprop_t &p_ou=props_ou[im_r_iqins_ijack_ind({im_ou,r_ou,iq_ins_ou,iclust})];
	  qprop_t &pr_bil=jbil[im_r_im_r_bilins_igam_ind({im_ou,r_ou,im_in,r_in,ipr_ins,iG})][iclust];
	  
	  pr_bil+=p_ou*quaGamma[iG]*quaGamma[5]*p_in.adjoint()*quaGamma[5];
	}
    }
}

void perens_t::compute_proj_bil(const vector<jqprop_t>& jprop_inv_in,const vector<jqprop_t>& jverts,const vector<jqprop_t>& jprop_inv_ou,const size_t ibilmom)
{
  const index_t ind({{"rest",im_r_im_r_ibil_ind.max()},{"ijack",njacks+1}});
  
  vector<tuple<size_t,size_t,size_t,size_t>> map;
#define ADD_COMBO(A,I,V,O)			\
  map.push_back({pr_bil::iins_of_ins[pr_bil::A],jqprop::iins_of_ins[jqprop::I],pr_bil::iins_of_ins[pr_bil::V],jqprop::iins_of_ins[jqprop::O]})
  ADD_COMBO(LO,  LO, LO, LO);
  switch(pars::use_QED)
    {
    case 0:
      break;
    case 1:
      ADD_COMBO(EX    , LO , EX    , LO);
      ADD_COMBO(CR_IN , LO , CR_IN , LO);
      ADD_COMBO(CR_OU , LO , CR_OU , LO);
      ADD_COMBO(TM_IN , LO , TM_IN , LO);
      ADD_COMBO(TM_OU , LO , TM_OU , LO);
      ADD_COMBO(PH_IN , LO , PH_IN , LO);
      ADD_COMBO(PH_OU , LO , PH_OU , LO);
      ///
      ADD_COMBO(CR_IN , CR , LO, LO);
      ADD_COMBO(CR_OU , LO , LO, CR);
      ADD_COMBO(TM_IN , TM , LO, LO);
      ADD_COMBO(TM_OU , LO , LO, TM);
      ADD_COMBO(PH_IN , PH , LO, LO);
      ADD_COMBO(PH_OU , LO , LO, PH);
      break;
    case 2:
      ADD_COMBO(EX     , LO  , EX     , LO);
      ADD_COMBO(QED_IN , LO  , QED_IN , LO);
      ADD_COMBO(QED_OU , LO  , QED_OU , LO);
      ///
      ADD_COMBO(QED_IN , QED , LO     , LO);
      ADD_COMBO(QED_OU , LO  , LO     , QED);
      break;
    }
#undef ADD_COMBO
  
#pragma omp parallel for
  for(size_t i=0;i<ind.max();i++)
    {
      //split im_r_im_r_ibil and ijack
      vector<size_t> i_comp=ind(i);
      const size_t im_r_im_r_ibil=i_comp[0],ijack=i_comp[1];
      const vector<size_t> im_r_im_r_ibil_comp=im_r_im_r_ibil_ind(im_r_im_r_ibil);
      
      //get im and r for fw and back, and ibil
      const size_t im_ou=im_r_im_r_ibil_comp[0],r_ou=im_r_im_r_ibil_comp[1];
      const size_t im_in=im_r_im_r_ibil_comp[2],r_in=im_r_im_r_ibil_comp[3];
      const size_t ibil=im_r_im_r_ibil_comp[4];
      
      for(auto t : map)
	{
	  const size_t pr_bilins=get<0>(t);
	  const size_t ijqins_in=get<1>(t);
	  const size_t bilins=get<2>(t);
	  const size_t ijqins_ou=get<3>(t);
	  
	  const size_t im_r_ijqins_ou=im_r_ijqins_ind({im_ou,r_ou,ijqins_ou});
	  const size_t im_r_ijqins_in=im_r_ijqins_ind({im_in,r_in,ijqins_in});
	  
	  const qprop_t& pinv_ou=jprop_inv_ou[im_r_ijqins_ou][ijack];
	  const qprop_t& pinv_in=jprop_inv_in[im_r_ijqins_in][ijack];
	  
	  //loop on all gammas
	  for(auto & iG : iG_of_bil[ibil])
	    {
	      vector<size_t> im_r_im_r_iG_comp=im_r_im_r_ibil_comp;
	      im_r_im_r_iG_comp[4]=iG;
	      const size_t im_r_im_r_bilins_igam=im_r_im_r_bilins_igam_ind({im_ou,r_ou,im_in,r_in,bilins,iG});
	      const size_t im_r_im_r_bilins_ibil_ibilmom=im_r_im_r_bilins_ibil_ibilmom_ind({im_ou,r_ou,im_in,r_in,pr_bilins,ibil,ibilmom});
	      const qprop_t &jv=jverts[im_r_im_r_bilins_igam][ijack];
	      
	      pr_bil[im_r_im_r_bilins_ibil_ibilmom][ijack]+=
		(pinv_ou*jv*quaGamma[5]*pinv_in.adjoint()*quaGamma[5]*quaGamma[iG].adjoint()).trace().real()/(12.0*iG_of_bil[ibil].size());
	    }
	}
    }
}

void perens_t::mom_compute_bil()
{
  //reset
  for(auto &t : get_pr_bil_tasks())
    *t.out=0.0;
  
  vector<raw_file_t> files=setup_read_all_qprops_mom(conf_list);
  
  for(size_t ibilmom=0;ibilmom<bilmoms.size();ibilmom++)
    {
      const size_t imom_ou=bilmoms[ibilmom][1];
      const size_t imom_in=bilmoms[ibilmom][2];
      const bool read_in=(imom_ou!=imom_in);
      
      vector<jqprop_t> jprops_ou(im_r_ijqins_ind.max());         //!< jackknived props out
      vector<jqprop_t> jprops_in(im_r_ijqins_ind.max());         //!< jackknived props in
      vector<jqprop_t> jverts(im_r_im_r_bilins_igam_ind.max());  //!< jackknived vertex
      
      for(size_t i_in_clust=0;i_in_clust<clust_size;i_in_clust++)
	for(size_t ihit=0;ihit<nhits_to_use;ihit++)
	  {
	    const size_t mom_ou=linmoms[imom_ou][0];
	    const size_t mom_in=linmoms[imom_in][0];
	    const size_t i_in_clust_ihit=i_in_clust_ihit_ind({i_in_clust,ihit});
	    cout<<"Working on qbil, "
	      "clust_entry "<<i_in_clust+1<<"/"<<clust_size<<", "
	      "hit "<<ihit+1<<"/"<<nhits<<", "
	      "momentum combo "<<ibilmom+1<<"/"<<bilmoms.size()<<", "
	      "moms: "<<mom_ou<<" "<<mom_in<<endl;
	    
	    //read
	    read_time.start();
	    const vector<qprop_t> props_ou=read_all_qprops_mom(files,i_in_clust_ihit,mom_ou);
	    const vector<qprop_t> props_in=(read_in?read_all_qprops_mom(files,i_in_clust_ihit,mom_in):props_ou);
	    read_time.stop();
	    
	    //build all props
	    build_props_time.start();
	    build_all_mr_jackkniffed_qprops(jprops_ou,props_ou);
	    build_all_mr_jackkniffed_qprops(jprops_in,props_in);
	    build_props_time.stop();
	    
	    //build all bilinear verts
	    build_verts_time.start();
	    build_all_mr_gbil_jackkniffed_verts(jverts,props_in,props_ou);
	    build_verts_time.stop();
	  }
      
      //clusterize
      clust_time.start();
      for(auto o : {&jprops_ou,&jprops_in,&jverts})
	for(auto &j : *o)
	  j.clusterize(clust_size);
      clust_time.stop();
      
      vector<jqprop_t> jprops_inv_ou=get_inverse_propagators(jprops_ou); //!< inverse propagator_ou
      vector<jqprop_t> jprops_inv_in=get_inverse_propagators(jprops_in); //!< inverse propagator_in
      
      proj_time.start();
      compute_proj_bil(jprops_inv_in,jverts,jprops_inv_ou,ibilmom);
      proj_time.stop();
    }
}

void perens_t::subtract_Oa2_pr_bil()
{
  const size_t iins=0;
  
  for(size_t im_in=0;im_in<nm;im_in++)
    for(size_t r_in=0;r_in<nr;r_in++)
      for(size_t im_ou=0;im_ou<nm;im_ou++)
	for(size_t r_ou=0;r_ou<nr;r_ou++)
	  for(size_t ibil=0;ibil<nbil;ibil++)
	    for(size_t ibilmom=0;ibilmom<bilmoms.size();ibilmom++)
	      {
		const size_t imom=bilmoms[0][0];
		pr_bil[im_r_im_r_bilins_ibil_ibilmom_ind({im_in,r_in,im_ou,r_ou,iins,ibil,ibilmom})]-=
		  g2to_correct()*pr_bil_a2(pars::act,gf::LANDAU,group::SU3,all_moms[bilmoms[imom][0]],L,ibil);
	      }
}

void perens_t::evolve_pr_bil(perens_t &out) const
{
  CRASH("");
}

vector<perens_t::task_t> perens_t::get_pr_bil_tasks(const vector<const perens_t*> &ens)
{
  vector<const djvec_t*> in_pr_bil;
  
  for(auto &e : ens)
    in_pr_bil.push_back(&e->pr_bil);
  
  vector<task_t> pr_bil_tasks;
  
  if(pars::compute_bilinears) pr_bil_tasks.push_back({&pr_bil,in_pr_bil,im_r_im_r_bilins_ibil_ibilmom_ind,"pr_bil",QCD_task});
  
  return pr_bil_tasks;
}

void perens_t::average_r_pr_bil(perens_t &out) const
{
  cout<<"Averaging r for pr_bil"<<endl;
  
  for(auto &t : out.get_pr_bil_tasks({this}))
    {
      cout<<" "<<t.tag<<endl;
      
      const djvec_t &pr=*t.in.front();
      djvec_t &pr_rave=*t.out;
      
#pragma omp parallel for
      for(size_t out_i=0;out_i<out.im_r_im_r_bilins_ibil_ibilmom_ind.max();out_i++)
	{
	  const vector<size_t> out_im_r_im_r_bilins_ibil_ibilmom_comp=out.im_r_im_r_bilins_ibil_ibilmom_ind(out_i);
	  vector<size_t> im_r_im_r_bilins_ibil_ibilmom_comp=out_im_r_im_r_bilins_ibil_ibilmom_comp;
	  
	  pr_rave[out_i]=0.0;
	  for(size_t r=0;r<nr;r++)
	    {
	      im_r_im_r_bilins_ibil_ibilmom_comp[1]=
	        im_r_im_r_bilins_ibil_ibilmom_comp[3]=r;
	      const size_t i=im_r_im_r_bilins_ibil_ibilmom_ind(im_r_im_r_bilins_ibil_ibilmom_comp);
	      pr_rave[out_i]+=pr[i];
	    }
	  pr_rave[out_i]/=nr;
	}
    }
}

void perens_t::average_equiv_momenta_pr_bil(perens_t &out,const vector<vector<size_t>> &equiv_bilmom_combos) const
{
#pragma omp parallel for
  for(size_t i=0;i<out.im_r_im_r_bilins_ibil_ibilmom_ind.max();i++)
    {
      const vector<size_t> out_im_r_im_r_bilins_ibil_ibilmom_comp=out.im_r_im_r_bilins_ibil_ibilmom_ind(i);
      const size_t out_imom_combo=out_im_r_im_r_bilins_ibil_ibilmom_comp[6];
      
      for(const auto &t : out.get_pr_bil_tasks({this}))
  	{
  	  djack_t &ave=(*t.out)[i];
  	  ave=0.0;
  	  for(const size_t ieq_mom : equiv_bilmom_combos[out_imom_combo])
	    {
	      vector<size_t> in_im_r_im_r_bilins_ibil_ibilmom_comp=out_im_r_im_r_bilins_ibil_ibilmom_comp;
	      in_im_r_im_r_bilins_ibil_ibilmom_comp[6]=ieq_mom;
	      const size_t ieq=im_r_im_r_bilins_ibil_ibilmom_ind(in_im_r_im_r_bilins_ibil_ibilmom_comp);
	      ave+=(*t.in.front())[ieq];
	    }
  	  ave/=equiv_bilmom_combos[out_imom_combo].size();
  	}
    }
}

void perens_t::val_chir_extrap_pr_bil(perens_t &out) const
{
  const index_t r_r_ibil_ibilmom_ind({{"r",nr},{"r",nr},{"bil",nbil},{"bilmoms",bilmoms.size()}});
  const index_t r_r_bilins_ibil_ibilmom_ind({{"r",nr},{"r",nr},{"bilins",pr_bil::nins},{"bil",nbil},{"bilmoms",bilmoms.size()}});
  
  for(auto &t : out.get_pr_bil_tasks({this}))
#pragma omp parallel for
    for(size_t r_r_ibil_ibilmom=0;r_r_ibil_ibilmom<r_r_ibil_ibilmom_ind.max();r_r_ibil_ibilmom++)
      {
	const vector<size_t> r_r_ibil_ibilmom_comps=r_r_ibil_ibilmom_ind(r_r_ibil_ibilmom);
	const size_t r_in=r_r_ibil_ibilmom_comps[0];
	const size_t r_ou=r_r_ibil_ibilmom_comps[1];
	const size_t ibil=r_r_ibil_ibilmom_comps[2];
	const size_t ibilmom=r_r_ibil_ibilmom_comps[3];
	
	const djvec_t &pr=*t.in.front();
	djvec_t &pr_chir=*t.out;
	const string &tag=t.tag;
	
	//check if we need to subtract the pole
	const bool sub_pole=(ibil==iS or ibil==iP);
	const size_t x_pow=(sub_pole?1:0);
	
	vector<djvec_t> coeffs(pr_bil::nins);
	for(size_t bilins=0;bilins<pr_bil::nins;bilins++)
	  {
	    const vector<size_t> r_r_bilins_ibil_ibilmom_comps={r_in,r_ou,bilins,ibil,ibilmom};
	    const size_t r_r_bilins_ibil_ibilmom=r_r_bilins_ibil_ibilmom_ind(r_r_bilins_ibil_ibilmom_comps);
	    
	    //open the plot file if needed
	    const string plot_path=dir_path+"/plots/chir_extr_"+tag+"_"+r_r_bilins_ibil_ibilmom_ind.descr(r_r_bilins_ibil_ibilmom)+".xmg";
	    grace_file_t *plot=nullptr;
	    if(ibilmom%pars::print_each_mom==0) plot=new grace_file_t(plot_path);
	    
	    //slice m and fit it
	    djvec_t y(nm*(nm+1)/2),y_plot(nm*(nm+1)/2);
	    vector<double> x(nm*(nm+1)/2);
	    for(size_t im_ou=0;im_ou<nm;im_ou++)
	      for(size_t im_in=im_ou;im_in<nm;im_in++)
		{
		  const size_t imeson=im_im_ind({im_ou,im_in});
		  
		  //compute mass sum
		  if(pars::chir_extr_method==chir_extr::MQUARK) x[imeson]=am[im_ou]+am[im_in];
		  else                                          x[imeson]=sqr(meson_mass[imeson].ave());
		  
		  if(std::isnan(x[imeson])) CRASH("Nanning %d",imeson);
		  
		  //compute y and y_plot
		  y_plot[imeson]=pr[im_r_im_r_bilins_ibil_ibilmom_ind({im_ou,r_ou,im_in,r_in,bilins,ibil,ibilmom})];
		  
		  //if QED case and pole must be subtracted, take into account variation due to leading pole
		  if(pars::use_QED and bilins>0)
		    {
		      const djack_t M=meson_mass[imeson],dM=meson_mass_QED[imeson];
		      if(sub_pole)
			{
			  const djack_t b0=coeffs[pr_bil::LO][2],c0=coeffs[pr_bil::LO][0];
			  const djack_t varb=2.0*b0*dM*M;
			  const djack_t varc=-2.0*c0*dM/(M*M*M);
			  y_plot[imeson]-=varb+varc;
			}
		      else
			if(pars::sub_meson_mass_shift_when_no_pole)
			  {
			    const djack_t b0=coeffs[pr_bil::LO][1];
			    const djack_t varb=2.0*b0*dM*M;
			    y_plot[imeson]-=varb;
			  }
		    }
		  
		  //fit x*y if pole present
		  y[imeson]=pow(x[imeson],x_pow)*y_plot[imeson];
		}
	    
	    //fit, store and write the result
	    coeffs[bilins]=poly_fit(x,y,1+x_pow);
	    const size_t iout=out.im_r_im_r_bilins_ibil_ibilmom_ind({0,r_ou,0,r_in,bilins,ibil,ibilmom});
	    pr_chir[iout]=coeffs[bilins][x_pow];
	    if(plot!=nullptr)
	      {
		auto xminmax=minmax_element(x.begin(),x.end());
		const double xmin=*xminmax.first*(sub_pole?0.5:0.01);
		const double xmax=*xminmax.second*1.1;
		write_fit_plot(*plot,xmin,xmax,[&coeffs,x_pow,bilins](double x)->djack_t{return poly_eval<djvec_t>(coeffs[bilins],x)/pow(x,x_pow);},x,y_plot);
		plot->write_ave_err(0.0,pr_chir[iout].ave_err());
	      }
	    
	    if(plot) delete plot;
	  }
      }
}

void perens_t::assemble_pr_bil_QED_greenfunctions()
{
  cout<<"Assembling pr_bil QED greenfunctions"<<endl;
  
  for(size_t ibilmom=0;ibilmom<bilmoms.size();ibilmom++)
    for(size_t im_r_im_r_ibil=0;im_r_im_r_ibil<im_r_im_r_ibil_ind.max();im_r_im_r_ibil++)
      {
	using namespace pr_bil;
	
	const vector<size_t> im_r_im_r_ibil_comp=im_r_im_r_ibil_ind(im_r_im_r_ibil);
	const size_t im_ou=im_r_im_r_ibil_comp[0];
	const size_t r_ou=im_r_im_r_ibil_comp[1];
	const size_t im_in=im_r_im_r_ibil_comp[2];
	const size_t r_in=im_r_im_r_ibil_comp[3];
	const size_t ibil=im_r_im_r_ibil_comp[4];
	const size_t im_r_ou=im_r_ind({im_ou,r_ou});
	const size_t im_r_in=im_r_ind({im_in,r_in});
	auto pr=pr_bil_ins_getter(im_ou,r_ou,im_in,r_in,ibil,ibilmom);
	
	switch(pars::use_QED)
	  {
	  case 0:
	    break;
	  case 1:
	    //basic part
	    pr(QED)=
	      pr(EX);
	    
	    //add critical counterterm
	    if(pars::use_deltam_cr_ct)
	      pr(QED)+=
		pr(CR_OU)*deltam_cr[im_r_ou]+
		pr(CR_IN)*deltam_cr[im_r_in];
	    
	    //add tm counterterm
	    if(pars::use_deltam_tm_ct)
	      pr(QED)+=
		pr(TM_OU)*deltam_tm[im_r_ou]+
		pr(TM_IN)*deltam_tm[im_r_in];
	    
	    //include self energy if needed
	    if(pars::include_self_energy_in_bilinears)
	      pr(QED)+=
		pr(PH_OU)+
		pr(PH_IN);
	    break;
	  case 2:
	    //basic part
	    pr(QED)=
	      pr(EX);
	    
	    //include self energy if needed
	    if(pars::include_self_energy_in_bilinears)
	      pr(QED)+=
		pr(QED_OU)+
		pr(QED_IN);
	    break;
	  }
      }
}
