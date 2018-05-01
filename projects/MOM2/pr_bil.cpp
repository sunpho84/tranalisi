#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_pr_BIL
 #include <MOM2/pr_bil.hpp>

#include <MOM2/perens.hpp>
#include <MOM2/prop.hpp>
#include <MOM2/timings.hpp>

void build_jackkniffed_vert_Gamma(jqprop_t &jvert,const qprop_t &prop_in,size_t iG,const qprop_t &prop_ou,size_t iclust)
{
  auto c=prop_ou*quaGamma[iG]*quaGamma[5]*prop_in.adjoint()*quaGamma[5];
  jvert[iclust]+=c;
  
  if(0)
  cout
    <<"iG: "<<iG<<
    ", iclust: "<<iclust<<
    ", prop_in: "<<prop_in(0,0)<<
    ", prop_ou: "<<prop_ou(0,0)<<
    ", res: "<<c(0,0)<<
    endl;
}

void perens_t::build_all_mr_gbil_jackkniffed_verts(vector<jbil_vert_t>& jbil,const vector<m_r_mom_conf_qprops_t>& props_in,const vector<m_r_mom_conf_qprops_t>& props_ou) const
{
  //! help finding the bilinear/clust combo
  index_t ind({{"i",im_r_im_r_igam_ind.max()},{"iclust",njacks}});
  
#pragma omp parallel for
  for(size_t i=0;i<ind.max();i++)
    {
      //decript bilinar/jack
      vector<size_t> ind_comp=ind(i);
      const size_t im_r_im_r_igam=ind_comp[0],iclust=ind_comp[1];
      
      //decript props to combine and which Gamma to create
      const vector<size_t> im_r_im_r_igam_comp=im_r_im_r_igam_ind(im_r_im_r_igam);
      const size_t im_fw=im_r_im_r_igam_comp[0],r_fw=im_r_im_r_igam_comp[1];
      const size_t im_bw=im_r_im_r_igam_comp[2],r_bw=im_r_im_r_igam_comp[3];
      const size_t iG=im_r_im_r_igam_comp[4];
      
      //proxy for vector and props
      const m_r_mom_conf_qprops_t &p_in=props_in[im_r_ijack_ind({im_fw,r_fw,iclust})];
      const m_r_mom_conf_qprops_t &p_ou=props_ou[im_r_ijack_ind({im_bw,r_bw,iclust})];
      jbil_vert_t &bil=jbil[im_r_im_r_igam];

      //create list of operations
      vector<tuple<jqprop_t*,const qprop_t*,const qprop_t*>> list={{&bil.LO,&p_in.LO,&p_ou.LO}};
      if(pars::use_QED)
	for(auto &o : vector<tuple<jqprop_t*,const qprop_t*,const qprop_t*>>({
	      {&bil.PH,&p_in.F,&p_ou.F},
	      {&bil.PH,&p_in.FF,&p_ou.LO},
	      {&bil.PH,&p_in.LO,&p_ou.FF},
	      {&bil.PH,&p_in.T,&p_ou.LO},
	      {&bil.PH,&p_in.LO,&p_ou.T},
	      {&bil.CR_CT_in,&p_in.P,&p_ou.LO},
	      {&bil.CR_CT_ou,&p_in.LO,&p_ou.P},
	      {&bil.TM_CT_in,&p_in.S,&p_ou.LO},
	      {&bil.TM_CT_ou,&p_in.LO,&p_ou.S}
	    }))
	  list.push_back(o);
      
      //create the vertex
      for(auto &o : list)
	build_jackkniffed_vert_Gamma(*get<0>(o),*get<1>(o),iG,*get<2>(o),iclust);
    }
}

void perens_t::clusterize_all_mr_jackkniffed_bilverts(vector<jbil_vert_t>& jverts) const
{
#pragma omp parallel for
  for(size_t ivert=0;ivert<jverts.size();ivert++)
    jverts[ivert].clusterize_all(clust_size);
}

void perens_t::compute_proj_bil(const vector<jm_r_mom_qprops_t>& jprop_inv_in,const vector<jbil_vert_t>& jverts,const vector<jm_r_mom_qprops_t>& jprop_inv_ou)
{
  for(auto &t : get_pr_bil_tasks())
    *t.out=0.0;
  
  const index_t ind({{"rest",im_r_im_r_ibil_ind.max()},{"ijack",njacks+1}});
  
#pragma omp parallel for
  for(size_t i=0;i<ind.max();i++)
    {
      //split im_r_im_r_ibil and ijack
      vector<size_t> i_comp=ind(i);
      const size_t im_r_im_r_ibil=i_comp[0],ijack=i_comp[1];
      const vector<size_t> im_r_im_r_ibil_comp=im_r_im_r_ibil_ind(im_r_im_r_ibil);
      
      //get im and r for fw and back, and ibil
      const size_t im_fw=im_r_im_r_ibil_comp[0],r_fw=im_r_im_r_ibil_comp[1];
      const size_t im_bw=im_r_im_r_ibil_comp[2],r_bw=im_r_im_r_ibil_comp[3];
      const size_t ibil=im_r_im_r_ibil_comp[4];
      
      const size_t im_r_in=im_r_ind({im_fw,r_fw});
      const size_t im_r_ou=im_r_ind({im_bw,r_bw});
      
      const jm_r_mom_qprops_t& pinv_ou=jprop_inv_ou[im_r_ou];
      const jm_r_mom_qprops_t& pinv_in=jprop_inv_in[im_r_in];
      
      //loop on all gammas
      for(auto & iG : iG_of_bil[ibil])
	{
	  vector<size_t> im_r_im_r_iG_comp=im_r_im_r_ibil_comp;
	  im_r_im_r_iG_comp[4]=iG;
	  const size_t im_r_im_r_iG=im_r_im_r_igam_ind(im_r_im_r_iG_comp);
	  const jbil_vert_t &jv=jverts[im_r_im_r_iG];
	  
#define PROJ(OUT,SIGN,OU,V,IN)						\
	  pr_bil_ ## OUT[im_r_im_r_ibil][ijack]+=			\
	    SIGN								\
	    (pinv_ou.OU[ijack]*jv.V[ijack]*quaGamma[5]*pinv_in.IN[ijack].adjoint()*quaGamma[5]*quaGamma[iG].adjoint()).trace().real()/(12.0*iG_of_bil[ibil].size())
	  
	  PROJ(LO, +,LO,LO,LO);
	  //
	  PROJ(PH, -,PH,LO,LO);
	  PROJ(PH, +,LO,PH,LO);
	  PROJ(PH, -,LO,LO,PH);
	  //
	  PROJ(CR_CT1, -,CR_CT,LO,LO);
	  PROJ(TM_CT1, -,TM_CT,LO,LO);
	  PROJ(CR_CT1, +,LO,CR_CT_ou,LO);
	  PROJ(TM_CT1, +,LO,TM_CT_ou,LO);
	  //
	  PROJ(CR_CT2, -,LO,LO,CR_CT);
	  PROJ(TM_CT2, -,LO,LO,TM_CT);
	  PROJ(CR_CT2, +,LO,CR_CT_in,LO);
	  PROJ(TM_CT2, +,LO,TM_CT_in,LO);
	  
	  #undef PROJ
	}
    }
}

void perens_t::mom_compute_bil()
{
  vector<raw_file_t> files=setup_read_all_qprops_mom(conf_list);
  
  for(size_t ibilmom=0;ibilmom<bilmoms.size();ibilmom++)
    {
      const size_t imom1=bilmoms[ibilmom][1];
      const size_t imom2=bilmoms[ibilmom][2];
      const bool read2=(imom1!=imom2);
      
      vector<jm_r_mom_qprops_t> jprops1(im_r_ind.max());      //!< jackknived props
      vector<jm_r_mom_qprops_t> jprops2(im_r_ind.max());      //!< jackknived props
      vector<jbil_vert_t> jverts(im_r_im_r_igam_ind.max());   //!< jackknived vertex
      
      for(size_t i_in_clust=0;i_in_clust<clust_size;i_in_clust++)
	for(size_t ihit=0;ihit<nhits_to_use;ihit++)
	  {
	    const size_t mom1=linmoms[imom1][0];
	    const size_t mom2=linmoms[imom2][0];
	    const size_t i_in_clust_ihit=i_in_clust_ihit_ind({i_in_clust,ihit});
	    cout<<"Working on qbil, "
	      "clust_entry "<<i_in_clust+1<<"/"<<clust_size<<", "
	      "hit "<<ihit+1<<"/"<<nhits<<", "
	      "momentum combo "<<ibilmom+1<<"/"<<bilmoms.size()<<", "
	      "moms: "<<mom1<<" "<<mom2<<endl;
	    
	    //read
	    read_time.start();
	    const vector<m_r_mom_conf_qprops_t> props1=read_all_qprops_mom(files,i_in_clust_ihit,mom1);
	    const vector<m_r_mom_conf_qprops_t> props2=(read2?read_all_qprops_mom(files,i_in_clust_ihit,mom2):props1);
	    read_time.stop();
	    
	    //build all props
	    build_props_time.start();
	    build_all_mr_jackkniffed_qprops(jprops1,props1);
	    build_all_mr_jackkniffed_qprops(jprops2,props2);
	    build_props_time.stop();
	    
	    //build all bilinear verts
	    build_verts_time.start();
	    build_all_mr_gbil_jackkniffed_verts(jverts,props1,props2);
	    build_verts_time.stop();
	  }
      
      //clusterize
      clust_time.start();
      clusterize_all_mr_jackkniffed_qprops(jprops1);
      clusterize_all_mr_jackkniffed_qprops(jprops2);
      clusterize_all_mr_jackkniffed_bilverts(jverts);
      clust_time.stop();
      
      vector<jm_r_mom_qprops_t> jprops_inv1; //!< inverse propagator1
      vector<jm_r_mom_qprops_t> jprops_inv2; //!< inverse propagator2
      
      get_inverse_propagators(jprops_inv1,jprops1);
      get_inverse_propagators(jprops_inv2,jprops2);
      
      CRASH("");
      proj_time.start();
      compute_proj_bil(jprops_inv1,jverts,jprops_inv2);
      
      // //QED
      // djvec_t pr_bil_QED_temp;
      // if(pars::use_QED)
      // 	{
      // 	  const djvec_t pr_bil_QED=compute_proj_bil(jprop_inv1,jverts.QED,jprop_inv2);
      // 	  const djvec_t pr_bil_a=compute_proj_bil(jprop_QED_inv1,jverts.LO,jprop_inv2);
      // 	  const djvec_t pr_bil_b=compute_proj_bil(jprop_inv1,jverts.LO,jprop_QED_inv2);
      // 	  pr_bil_QED_temp=
      // 	    (-pr_bil_a
      // 	     -pr_bil_b
      // 	     +pr_bil_QED);
      // 	}
      
      // //! an index running on all packed combo, and momenta
      // const index_t all_ibilmom_ind({{"All",im_r_im_r_ibil_ind.max()},{"bilmom",bilmoms.size()}});
      
      // //store
      // vector<pair<djvec_t*,djvec_t*>> store_task={{&pr_bil,&pr_bil_temp}};
      // if(pars::use_QED) store_task.push_back({&pr_bil_QED,&pr_bil_QED_temp});
      // for(auto &t : store_task)
      // 	for(size_t iall=0;iall<im_r_im_r_ibil_ind.max();iall++)
      // 	  (*t.first)[all_ibilmom_ind({iall,ibilmom})]=(*t.second)[iall];
      
      // proj_time.stop();
    }
}

vector<perens_t::task_t> perens_t::get_pr_bil_tasks(const vector<const perens_t*> &ens)
{
  vector<const djvec_t*>
    in_pr_bil_LO,
    in_pr_bil_CR_CT1,
    in_pr_bil_CR_CT2,
    in_pr_bil_TM_CT1,
    in_pr_bil_TM_CT2,
    in_pr_bil_PH;
  for(auto &e : ens)
    {
      in_pr_bil_LO.push_back(&e->pr_bil_LO);
      if(pars::use_QED)
	{
	  in_pr_bil_CR_CT1.push_back(&e->pr_bil_CR_CT1);
	  in_pr_bil_CR_CT2.push_back(&e->pr_bil_CR_CT2);
	  in_pr_bil_TM_CT1.push_back(&e->pr_bil_TM_CT1);
	  in_pr_bil_TM_CT2.push_back(&e->pr_bil_TM_CT2);
	  in_pr_bil_PH.push_back(&e->pr_bil_PH);
	}
    }
  
  vector<task_t> pr_bil_tasks={{&pr_bil_LO,in_pr_bil_LO,im_r_im_r_ibil_ibilmom_ind,"pr_bil_LO",QCD_task}};
  if(pars::use_QED)
    {
      pr_bil_tasks.push_back({&pr_bil_CR_CT1,in_pr_bil_CR_CT1,im_r_im_r_ibil_ibilmom_ind,"pr_bil_CR_CT1",QED_task});
      pr_bil_tasks.push_back({&pr_bil_CR_CT2,in_pr_bil_CR_CT2,im_r_im_r_ibil_ibilmom_ind,"pr_bil_CR_CT2",QED_task});
      pr_bil_tasks.push_back({&pr_bil_TM_CT1,in_pr_bil_TM_CT1,im_r_im_r_ibil_ibilmom_ind,"pr_bil_TM_CT1",QED_task});
      pr_bil_tasks.push_back({&pr_bil_TM_CT2,in_pr_bil_TM_CT2,im_r_im_r_ibil_ibilmom_ind,"pr_bil_TM_CT2",QED_task});
      pr_bil_tasks.push_back({&pr_bil_PH,in_pr_bil_PH,im_r_im_r_ibil_ibilmom_ind,"pr_bil_PH",QED_task});
    }
  
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
      
      for(size_t out_i=0;out_i<out.im_r_im_r_ibil_ibilmom_ind.max();out_i++)
	{
	  const vector<size_t> out_im_r_im_r_ibil_ibilmom_comp=out.im_r_im_r_ibil_ibilmom_ind(out_i);
	  vector<size_t> im_r_im_r_ibil_ibilmom_comp=out_im_r_im_r_ibil_ibilmom_comp;
	  
	  pr_rave[out_i]=0.0;
	  for(size_t r=0;r<nr;r++)
	    {
	      im_r_im_r_ibil_ibilmom_comp[1]=
		im_r_im_r_ibil_ibilmom_comp[3]=r;
	      const size_t i=im_r_im_r_ibil_ibilmom_ind(im_r_im_r_ibil_ibilmom_comp);
	      pr_rave[out_i]+=pr[i];
	    }
	  pr_rave[out_i]/=nr;
	}
    }
}

void perens_t::average_equiv_momenta_pr_bil(perens_t &out,const vector<vector<size_t>> &equiv_bilmom_combos) const
{
  for(size_t i=0;i<out.im_r_im_r_ibil_ibilmom_ind.max();i++)
    {
      const vector<size_t> out_im_r_im_r_ibil_ibilmom_comp=out.im_r_im_r_ibil_ibilmom_ind(i);
      const size_t out_imom_combo=out_im_r_im_r_ibil_ibilmom_comp[5];
      
      for(const auto &t : out.get_pr_bil_tasks({this}))
  	{
  	  djack_t &ave=(*t.out)[i];
  	  ave=0.0;
  	  for(const size_t ieq_mom : equiv_bilmom_combos[out_imom_combo])
	    {
	      vector<size_t> im_r_im_r_ibil_ibilmom_comp=out_im_r_im_r_ibil_ibilmom_comp;
	      im_r_im_r_ibil_ibilmom_comp[5]=ieq_mom;
	      const size_t ieq=im_r_im_r_ibil_ibilmom_ind(im_r_im_r_ibil_ibilmom_comp);
	      ave+=(*t.in.front())[ieq];
	    }
  	  ave/=equiv_bilmom_combos[out_imom_combo].size();
  	}
    }
}

void perens_t::val_chir_extrap_pr_bil(perens_t &out) const
{
  for(size_t ibilmom=0;ibilmom<bilmoms.size();ibilmom++)
    for(auto &t : out.get_pr_bil_tasks({this}))
      for(size_t ibil=0;ibil<nbil;ibil++)
	{
	  const djvec_t &pr=*t.in.front();
	  djvec_t &pr_chir=*t.out;
	  const string &tag=t.tag;
	  
	  //check if we need to subtract the pole
	  const bool sub_pole=(ibil==iS or ibil==iP);
	  const size_t x_pow=(sub_pole?
			      (t.QCD_QED_task==QED_task?1:1)
			      :0);
	  
	  //open the plot file if needed
	  const string plot_path=dir_path+"/plots/chir_extr_"+tag+"_"+bil_tag[ibil]+"_bilmom_"+to_string(ibilmom)+".xmg";
	  grace_file_t *plot=nullptr;
	  if(ibilmom%pars::print_each_mom==0) plot=new grace_file_t(plot_path);
	  
	  for(size_t r1=0;r1<nr;r1++)
	    for(size_t r2=0;r2<nr;r2++)
	      {
		//slice m and fit it
		djvec_t y(nm*(nm+1)/2),y_plot(nm*(nm+1)/2);
		vector<double> x(nm*(nm+1)/2);
		int i=0;
		for(size_t im1=0;im1<nm;im1++)
		  for(size_t im2=im1;im2<nm;im2++)
		    {
		      //compute mass sum
		      if(pars::chir_extr_method==chir_extr::MQUARK) x[i]=am[im1]+am[im2];
		      else                                          x[i]=sqr(meson_mass[im_im_ind({im1,im2})].ave());
		      
		      if(std::isnan(x[i])) CRASH("Nanning %d",i);
		      
		      //compute y and y_plot
		      y_plot[i]=pr[im_r_im_r_ibil_ibilmom_ind({im1,r1,im2,r2,ibil,ibilmom})];
		      //fit x*y if pole present
		      y[i]=pow(x[i],x_pow)*y_plot[i];
		      //increment the number of mass combos
		      i++;
		    }
		
		//fit, store and write the result
		const djvec_t coeffs=poly_fit(x,y,1+x_pow);
		const size_t iout=out.im_r_im_r_ibil_ibilmom_ind({0,r1,0,r2,ibil,ibilmom});
		pr_chir[iout]=coeffs[x_pow];
		if(plot!=nullptr)
		  {
		    auto xminmax=minmax_element(x.begin(),x.end());
		    double xmin=*xminmax.first*0.5;
		    double xmax=*xminmax.second*1.1;
		    write_fit_plot(*plot,xmin,xmax,[&coeffs,sub_pole,x_pow](double x)->djack_t{return poly_eval<djvec_t>(coeffs,x)/pow(x,x_pow);},x,y_plot);
		    plot->write_ave_err(0.0,pr_chir[iout].ave_err());
		  }
	      }
	  if(plot) delete plot;
	}
}
