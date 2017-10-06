#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <tranalisi.hpp>

#include <contractions.hpp>
#include <corrections.hpp>
#include <geometry.hpp>
#include <types.hpp>

#include <prop.hpp>
#include <sig3.hpp>
#include <Zbil.hpp>
#include <Zq.hpp>
#include <Zq_sig1.hpp>

const bool use_QED=true;

string suff_hit="";

index_t conf_ind; //!< index of a conf given ijack and i_in_clust
index_t im_r_ind; //!< index of im,r
index_t r_imom_ind; //!< index of r,imom combo
index_t im_r_imom_ind; //!< index of im,r,imom combo
index_t r_ind_imom_ind; //!< index of r,imom combo
index_t im_r_ind_imom_ind; //!< index of im,r,imom combo
index_t i_in_clust_ihit_ind; //!< index of i_in_clust,ihit

//! prepare a list of reading task, to be executed in parallel
vector<task_list_t> prepare_read_prop_taks(vector<m_r_mom_conf_props_t> &props,const vector<size_t> &conf_list)
{
  const index_t im_r_ijack_ind=im_r_ind*index_t({{"ijack",njacks}}); //!< index of im,r,ijack combo
  props.resize(im_r_ijack_ind.max());
  
  //! tasks
  vector<task_list_t> read_tasks(i_in_clust_ihit_ind.max());
  
  for(size_t ijack=0;ijack<njacks;ijack++)
    {
      for(size_t im=0;im<nm;im++)
	for(size_t r=0;r<nr;r++)
	  {
	    size_t im_r_ijack=im_r_ijack_ind({im,r,ijack});
	    m_r_mom_conf_props_t &l=props[im_r_ijack];
	    
	    //add EM if asked
	    using tup_in_t=tuple<prop_t*,string,dcompl_t>;
	    vector<tup_in_t> list={{&l.LO,"0",1}};
	    if(use_QED)
	      {
		list.push_back(tup_in_t(&l.P,"P",dcompl_t(0,tau3[r])));
		list.push_back(tup_in_t(&l.S,"S",dcompl_t(-1,0)));
		list.push_back(tup_in_t(&l.T,"T",dcompl_t(1,0)));
		list.push_back(tup_in_t(&l.F,"F",dcompl_t(1,0)));
		list.push_back(tup_in_t(&l.FF,"FF",dcompl_t(1,0)));
	      }
	    
	    for(auto &psc : list)
	      for(size_t i_i_in_clust_ihit=0;i_i_in_clust_ihit<i_in_clust_ihit_ind.max();i_i_in_clust_ihit++)
		{
		  vector<size_t> i_in_clust_ihit=i_in_clust_ihit_ind(i_i_in_clust_ihit);
		  size_t i_in_clust=i_in_clust_ihit[0],ihit=i_in_clust_ihit[1];
		  size_t iconf=conf_ind({ijack,i_in_clust});
		  string path=combine("out/%04zu/fft_",conf_list[iconf])+get_prop_tag(im,r,get<1>(psc))+combine(suff_hit.c_str(),ihit);
		  read_tasks[i_i_in_clust_ihit].push_back(incapsulate_task(read_prop,get<0>(psc),raw_file_t(path,"r"),get<2>(psc)));
		}
	  }
    }
  
  return read_tasks;
}

//! write a given Z
void write_Z(const string &name,const djvec_t &Z,const vector<double> &pt2)
{
  grace_file_t outf("plots/"+name+".xmg");
  outf.write_vec_ave_err(pt2,Z.ave_err());
}

//! linearly fit a given Z
djvec_t linfit_Z(const djvec_t &Z,const string &name,double band_val=0)
{
  vector<double> pt2=get_indep_pt2();
  double p2max=*max_element(pt2.begin(),pt2.end());
  djvec_t pars=poly_fit(pt2,Z,1,1.0,2.0);
  
  grace_file_t outf("plots/"+name+".xmg");
  outf.write_vec_ave_err(pt2,Z.ave_err());
  outf.write_polygon(bind(poly_eval<djvec_t>,pars,_1),0,p2max);
  if(band_val!=0) outf.write_line([band_val](double x){return band_val;},0,p2max);
  
  return pars;
}

int main(int narg,char **arg)
{
  const size_t ntimers=9;
  time_stats_t ts(ntimers);
  //compute the partial times
  stopwatch_t &read_time=ts.add("read propagators");
  stopwatch_t &build_props_time=ts.add("build props");
  stopwatch_t &build_verts_time=ts.add("build verts");
  stopwatch_t &clust_time=ts.add("clusterize");
  stopwatch_t &invert_time=ts.add("invert the props");
  stopwatch_t &finish_EM_time=ts.add("finish EM");
  stopwatch_t &proj_time=ts.add("project bilinears");
  stopwatch_t &Zq_time=ts.add("compute Zq");
  stopwatch_t &deltam_cr_time=ts.add("compute deltam_cr");
  
  ///////////////////////////////////////////////////////////////////////////
  
  //read input file
  string input_path="input.txt";
  if(narg>=2) input_path=arg[1];
  
  //open input
  raw_file_t input(input_path,"r");
  
  size_t Ls=input.read<size_t>("L"); //!< lattice spatial size
  L[0]=input.read<size_t>("T");
  
  const string act_str=input.read<string>("Action"); //!< action name
  auto act_key=gaz_decr.find(act_str); //!< key in the map of act
  if(act_key==gaz_decr.end()) CRASH("Unable to decript %s",act_str.c_str());
  gaz_t act=act_key->second; //!< action
  
  const double beta=input.read<double>("Beta"); //!< beta
  const double plaq=input.read<double>("Plaq"); //!< plaquette
  nm=input.read<double>("Nm");
  am.resize(nm);
  for(size_t im=0;im<nm;im++) am[im]=input.read<double>();
  const double am_max=*max_element(am.begin(),am.end())*1.1;
  const double am_min=*min_element(am.begin(),am.end())/1.1;
  nr=input.read<double>("Nr");
  
  size_t im_sea=input.read<double>("ImSea"); //!< index of sea mass
  
  const string mom_list_path=input.read<string>("MomList"); //!< list of momenta
  const size_t ext_njacks=input.read<size_t>("NJacks"); //!< number of jacks
  
  //! conf range
  range_t conf_range;
  input.expect("ConfRange");
  conf_range.start=input.read<size_t>();
  conf_range.each=input.read<size_t>();
  conf_range.end=input.read<size_t>();
  
  const size_t nhits=input.read<size_t>("NHits"); //!< number of hits
  const size_t nhits_to_use=input.read<size_t>("NHitsToUse"); //!< number of hits to be used
  
  const size_t tmin=input.read<size_t>("Tmin");
  const size_t tmax=input.read<size_t>("Tmax");
  
  //////////////////////////////////////////////////
  
  //set the number of jackknives
  set_njacks(ext_njacks);
  
  //! sufffix if a number of hits is different from 1
  if(nhits>1) suff_hit="_hit_%zu";
  
  double g2=6.0/beta; //!< coupling
  double g2tilde=g2/plaq; //!< boosted coupling
  cout<<"g2tilde: "<<g2tilde<<endl;
  
  //set the coefficients
  set_pr_bil_a2(act);
  
  //set spatial sizes
  V=L[0]*pow(Ls,NDIM-1);
  for(size_t mu=1;mu<NDIM;mu++) L[mu]=Ls;
  
  //initialize momenta
  set_list_of_moms(mom_list_path);
  set_class_of_equiv_moms();
  //list_all_smom_pairs();
  
  string test_path="out/%04zu/fft_S_M0_R0_0";
  if(nhits>1) test_path+="_hit_0";
  vector<size_t> conf_list=get_existing_paths_in_range(test_path,conf_range); //!< list of existing confs
  if(conf_list.size()==0) CRASH("list of configurations is empty! check %s ",test_path.c_str());
  
  //compute deltam_cr
  deltam_cr_time.start();
  djvec_t deltam_cr(nm);
  for(size_t im=0;im<nm;im++)
    {
      deltam_cr[im]=compute_deltam_cr(conf_list,tmin,tmax,im);
      cout<<"Deltam cr["<<im<<"]: "<<smart_print(deltam_cr[im])<<endl;
    }
  deltam_cr_time.stop();
  
  size_t clust_size=trim_to_njacks_multiple(conf_list,true); //!< cluster size
  
  conf_ind.set_ranges({{"ijack",njacks},{"i_in_clust",clust_size}});
  im_r_ind.set_ranges({{"m",nm},{"r",nr}});
  r_imom_ind.set_ranges({{"r",nr},{"imom",imoms.size()}});
  im_r_imom_ind.set_ranges({{"m",nm},{"r",nr},{"imom",imoms.size()}});
  r_ind_imom_ind.set_ranges({{"r",nr},{"ind_mom",equiv_imoms.size()}});
  im_r_ind_imom_ind.set_ranges({{"m",nm},{"r",nr},{"ind_mom",equiv_imoms.size()}});
  i_in_clust_ihit_ind.set_ranges({{"i_in_clust",clust_size},{"ihit",nhits_to_use}});
  const index_t im_r_im_r_igam_ind=im_r_ind*im_r_ind*index_t({{"igamma",nGamma}});
  const index_t r_r_iZbil_ind({{"r",nr},{"r",nr},{"iZbil",nZbil}});
  const index_t im_r_im_r_iZbil_ind=im_r_ind*im_r_ind*index_t({{"iZbil",nZbil}});
  const index_t r_r_iZbil_imom_ind=r_r_iZbil_ind*index_t({{"imom",imoms.size()}});
  const index_t im_r_im_r_iZbil_imom_ind=im_r_im_r_iZbil_ind*index_t({{"imom",imoms.size()}});
  const index_t im_r_ijack_ind=im_r_ind*index_t({{"ijack",njacks}});
  const index_t im_r_ijackp1_ind=im_r_ind*index_t({{"ijack",njacks+1}});
  
  //Zq for all moms, with and without EM
  djvec_t Zq_allmoms(im_r_imom_ind.max());
  djvec_t Zq_sig1_allmoms(im_r_imom_ind.max());
  djvec_t Zq_sig1_EM_allmoms;
  if(use_QED) Zq_sig1_EM_allmoms.resize(im_r_imom_ind.max());
  
  //Subtracted Zq, with and without EM, all moms
  djvec_t Zq_chir_allmoms(r_imom_ind.max());
  djvec_t Zq_sig1_chir_allmoms(r_imom_ind.max());
  djvec_t Zq_sig1_EM_chir_allmoms;
  if(use_QED) Zq_sig1_EM_chir_allmoms.resize(r_imom_ind.max());
  
  //! task to chiral extrapolate
  using Z_task_t=tuple<djvec_t*,djvec_t*,string>;
  
  //! list of task to chirally extrapolate Zq
  vector<Z_task_t> Zq_tasks{
    {&Zq_allmoms,&Zq_chir_allmoms,string("Zq")},
    {&Zq_sig1_allmoms,&Zq_sig1_chir_allmoms,"Zq_sig1"}};
  if(use_QED) Zq_tasks.push_back(make_tuple(&Zq_sig1_EM_allmoms,&Zq_sig1_EM_chir_allmoms,"Zq_sig1_EM"));
  
  djvec_t Zbil_allmoms(im_r_im_r_iZbil_imom_ind.max());
  djvec_t Zbil_chir_allmoms(r_r_iZbil_imom_ind.max());
  djvec_t Zbil_QED_chir_allmoms(im_r_im_r_iZbil_imom_ind.max());
  djvec_t Zbil_QED_allmoms(r_r_iZbil_imom_ind.max());
  
  vector<m_r_mom_conf_props_t> props; //!< store props for individual conf
  
  vector<task_list_t> read_tasks=prepare_read_prop_taks(props,conf_list);
  vector<jm_r_mom_props_t> jprops(im_r_ind.max()); //!< jackknived props
  jbil_vert_t jverts(im_r_im_r_igam_ind.max(),use_QED); //!< jackknived vertex
  
  for(size_t imom=0;imom<imoms.size();imom++)
    {
      for(size_t i_in_clust=0;i_in_clust<clust_size;i_in_clust++)
	for(size_t ihit=0;ihit<nhits_to_use;ihit++)
	  {
	    size_t i_in_clust_hit=i_in_clust_ihit_ind({i_in_clust,ihit});
	    cout<<"Working on clust_entry "<<i_in_clust+1<<"/"<<clust_size<<", hit "<<ihit+1<<"/"<<nhits<<", momentum "<<imom+1<<"/"<<imoms.size()<<endl;
	    read_time.start();
	    read_tasks[i_in_clust_hit].assolve_all(RECYCLE);
	    read_time.stop();
	    
	    //build all props
	    build_props_time.start();
	    build_all_mr_jackknifed_props(jprops,props,use_QED,im_r_ind);
	    build_props_time.stop();
	    
	    build_verts_time.start();
	    build_all_mr_gbil_jackknifed_verts(jverts,props,im_r_im_r_igam_ind,im_r_ijack_ind,use_QED);
	    build_verts_time.stop();
	  }
      
      //clusterize
      clust_time.start();
      clusterize_all_mr_jackknifed_props(jprops,use_QED,clust_size);
      jverts.clusterize_all(use_QED,clust_size);
      clust_time.stop();
      
      //finish EM
      if(use_QED)
	{
	  finish_EM_time.start();
	  finish_jverts_EM(jverts,deltam_cr[im_sea]);
	  finish_jprops_EM(jprops,deltam_cr[im_sea]);
	  finish_EM_time.stop();
	}
      
      vector<jprop_t> jprop_inv(im_r_ind.max()); //!< inverse propagator
      vector<jprop_t> jprop_EM_inv(im_r_ind.max()); //!< inverse propagator with em insertion
#pragma omp parallel for reduction(+:invert_time,Zq_time)
      for(size_t im_r_ijack=0;im_r_ijack<im_r_ijackp1_ind.max();im_r_ijack++)
	{
	  //decript indices
	  vector<size_t> im_r_ijack_comps=im_r_ijackp1_ind(im_r_ijack);
	  size_t im=im_r_ijack_comps[0],r=im_r_ijack_comps[1],ijack=im_r_ijack_comps[2];
	  size_t im_r=im_r_ind({im,r});
	  size_t im_r_imom=im_r_imom_ind({im,r,imom});
	  
	  //compute inverse
	  invert_time.start();
	  prop_t prop_inv=jprops[im_r].LO[ijack].inverse();
	  jprop_inv[im_r][ijack]=prop_inv;
	  invert_time.stop();
	  
	  //compute Zq
	  Zq_time.start();
	  Zq_allmoms[im_r_imom][ijack]=compute_Zq(prop_inv,imom);
	  Zq_sig1_allmoms[im_r_imom][ijack]=compute_Zq_sig1(prop_inv,imom);
	  Zq_time.stop();
	  
	  //do the same with QED
	  if(use_QED)
	    {
	      invert_time.start();
	      prop_t prop_EM_inv=prop_inv*jprops[im_r].EM[ijack]*prop_inv;
	      jprop_EM_inv[im_r][ijack]=prop_EM_inv;
	      invert_time.stop();
	      
	      Zq_time.start();
	      Zq_sig1_EM_allmoms[im_r_imom][ijack]=-compute_Zq_sig1(prop_EM_inv,imom);
	      Zq_time.stop();
	    }
	}
      
      //extrapolate to chiral limit Zq
      for(size_t r=0;r<nr;r++)
	for(auto & p : Zq_tasks)
	  {
	    djvec_t y(nm);
	    const djvec_t &Z=(*get<0>(p));
	    djvec_t &Z_chir=(*get<1>(p));
	    const string &tag=get<2>(p);
	    
	    //slice m
	    double am_max=*max_element(am.begin(),am.end())*1.1;
	    for(size_t im=0;im<nm;im++) y[im]=Z[im_r_imom_ind({im,r,imom})];
	    //fit and write the result
	    djvec_t coeffs=poly_fit(am,y,1,am_min,am_max);
	    if(imom%20==0) write_poly_fit_plot("plots/chir_extr_"+tag+"_r_"+to_string(r)+"_mom_"+to_string(imom)+".xmg",0,am_max,coeffs,am,y);
	    //extrapolated value
	    Z_chir[r_imom_ind({r,imom})]=coeffs[0];
	  }
      
      proj_time.start();
      djvec_t pr_bil_allmoms=compute_proj_bil(jprop_inv,jverts.LO,jprop_inv,im_r_ind);
      djvec_t pr_bil_chir_allmoms(r_r_iZbil_ind.max());
      vector<Z_task_t> pr_bil_tasks{{&pr_bil_allmoms,&pr_bil_chir_allmoms,string("pr_bil")}};
      
      //QED
      djvec_t pr_bil_EM_allmoms,pr_bil_a_allmoms,pr_bil_b_allmoms;
      djvec_t pr_bil_EM_chir_allmoms(r_r_iZbil_ind.max()),pr_bil_a_chir_allmoms(r_r_iZbil_ind.max()),pr_bil_b_chir_allmoms(r_r_iZbil_ind.max());
      if(use_QED)
	{
	  pr_bil_EM_allmoms=compute_proj_bil(jprop_inv,jverts.EM,jprop_inv,im_r_ind);
	  pr_bil_a_allmoms=compute_proj_bil(jprop_EM_inv,jverts.LO,jprop_inv,im_r_ind);
	  pr_bil_b_allmoms=compute_proj_bil(jprop_inv,jverts.LO,jprop_EM_inv,im_r_ind);
	  
	  pr_bil_tasks.push_back(make_tuple(&pr_bil_EM_allmoms,&pr_bil_EM_chir_allmoms,string("pr_bil_EM")));
	  pr_bil_tasks.push_back(make_tuple(&pr_bil_a_allmoms,&pr_bil_a_chir_allmoms,string("pr_bil_a")));
	  pr_bil_tasks.push_back(make_tuple(&pr_bil_b_allmoms,&pr_bil_b_chir_allmoms,string("pr_bil_b")));
	}
      proj_time.stop();
      
      //extrapolate to chiral limit the projected bilinears
      for(size_t r1=0;r1<nr;r1++)
	for(size_t r2=0;r2<nr;r2++)
	  for(size_t iZbil=0;iZbil<nZbil;iZbil++)
	    for(auto & p : pr_bil_tasks)
	      {
		const djvec_t &pr=(*get<0>(p));
		djvec_t &pr_chir=(*get<1>(p));
		const string &tag=get<2>(p);
		
		//check if we need to subtract the pole
		const bool sub_pole=(iZbil==iZS or iZbil==iZP);
		
		//slice m and fit it
		djvec_t y(nm*(nm+1)/2),y_plot(nm*(nm+1)/2);
		vector<double> x(nm*(nm+1)/2);
		int i=0;
		for(size_t im1=0;im1<nm;im1++)
		  for(size_t im2=im1;im2<nm;im2++)
		    {
		      //compute mass sum
		      x[i]=am[im1]+am[im2];
		      //compute y and y_plot
		      y_plot[i]=pr[im_r_im_r_iZbil_ind({im1,r1,im2,r2,iZbil})];
		      if(sub_pole) y[i]=x[i]*y_plot[i];
		      else         y[i]=y_plot[i];
		      //increment
		      i++;
		    }
		
		//fit and write the result
		djvec_t coeffs=poly_fit(x,y,(sub_pole?2:1),2.0*am_min,2.0*am_max);
		if(imom%20==0)
		  {
		    const string path="plots/chir_extr_"+tag+"_"+Zbil_tag[iZbil]+"_r1_"+to_string(r1)+"_r2_"+to_string(r2)+"_mom_"+to_string(imom)+".xmg";
		    write_fit_plot(path,2*am_min,2*am_max,[&coeffs,sub_pole](double x)->djack_t{return poly_eval<djvec_t>(coeffs,x)/(sub_pole?x:1);},x,y_plot);
		  }
		//extrapolated value
		pr_chir[r_r_iZbil_ind({r1,r2,iZbil})]=coeffs[sub_pole?1:0];
	      }
      
      //build Z
      for(size_t im_r_im_r_iZbil=0;im_r_im_r_iZbil<im_r_im_r_iZbil_ind.max();im_r_im_r_iZbil++)
	{
	  const vector<size_t> im_r_im_r_iZbil_comp=im_r_im_r_iZbil_ind(im_r_im_r_iZbil);
	  const vector<size_t> im_r1_comp=subset(im_r_im_r_iZbil_comp,0,2);
	  const vector<size_t> im_r2_comp=subset(im_r_im_r_iZbil_comp,2,4);
	  const size_t iZbil=im_r_im_r_iZbil_comp[4];
	  const size_t im_r1_imom=im_r_imom_ind(concat(im_r1_comp,imom));
	  const size_t im_r2_imom=im_r_imom_ind(concat(im_r2_comp,imom));
	  
	  const size_t im_r_im_r_iZbil_imom=im_r_im_r_iZbil_imom_ind(concat(im_r1_comp,im_r2_comp,vector<size_t>({iZbil,imom})));
	  
	  Zbil_allmoms[im_r_im_r_iZbil_imom]=
	    sqrt(Zq_allmoms[im_r1_imom]*Zq_allmoms[im_r2_imom])/pr_bil_allmoms[im_r_im_r_iZbil];
	  
	  if(use_QED)
	    {
	      Zbil_QED_allmoms[im_r_im_r_iZbil_imom]=
		(pr_bil_a_allmoms[im_r_im_r_iZbil]+pr_bil_b_allmoms[im_r_im_r_iZbil]-pr_bil_EM_allmoms[im_r_im_r_iZbil])/pr_bil_allmoms[im_r_im_r_iZbil]+
		(Zq_sig1_EM_allmoms[im_r1_imom]/Zq_sig1_allmoms[im_r1_imom]+Zq_sig1_EM_allmoms[im_r2_imom]/Zq_sig1_allmoms[im_r2_imom])/2.0;
	    }
	}
      
      //build Z
      for(size_t r_r_iZbil=0;r_r_iZbil<r_r_iZbil_ind.max();r_r_iZbil++)
	{
	  const vector<size_t> r_r_iZbil_comp=r_r_iZbil_ind(r_r_iZbil);
	  const size_t r1=r_r_iZbil_comp[0],r2=r_r_iZbil_comp[1];
	  const size_t r_r_iZbil_imom=r_r_iZbil_imom_ind(concat(r_r_iZbil_comp,imom));
	  const size_t r1_imom=im_r_imom_ind({r1,imom});
	  const size_t r2_imom=im_r_imom_ind({r2,imom});
	  
	  Zbil_chir_allmoms[r_r_iZbil_imom]=
	    sqrt(Zq_chir_allmoms[r1_imom]*Zq_allmoms[r2_imom])/pr_bil_chir_allmoms[r_r_iZbil];
	  
	  if(use_QED)
	    {
	      Zbil_QED_chir_allmoms[r_r_iZbil_imom]=
		(pr_bil_a_chir_allmoms[r_r_iZbil]+pr_bil_b_chir_allmoms[r_r_iZbil]-pr_bil_EM_chir_allmoms[r_r_iZbil])/pr_bil_chir_allmoms[r_r_iZbil]+
		(Zq_sig1_EM_chir_allmoms[r1_imom]/Zq_sig1_chir_allmoms[r1_imom]+Zq_sig1_EM_chir_allmoms[r2_imom]/Zq_sig1_chir_allmoms[r2_imom])/2.0;
	    }
	}
      
    }
  
  // //compute Zq, Zq_sig1 and Zbil for all moms
  // djvec_t sig3_allmoms=compute_sig3(jprop_inv);
  // djvec_t sig3_em_allmoms=compute_sig3(jprop_em_inv);
  
  
  // //correct Z LO
  // djvec_t Zq_allmoms_sub=Zq_allmoms;
  // djvec_t Zq_sig1_allmoms_sub=Zq_sig1_allmoms;
  // vector<djvec_t> pr_bil_allmoms_sub=pr_bil_allmoms;
  // for(size_t imom=0;imom<imoms.size();imom++)
  //   {
  //     imom_t mom=imoms[imom];
  //     Zq_allmoms_sub[imom]-=g2tilde*sig1_a2(act,mom,L);
  //     Zq_sig1_allmoms_sub[imom]-=g2tilde*sig1_a2(act,mom,L);
  //     for(size_t iZbil=0;iZbil<nZbil;iZbil++)
  // 	pr_bil_allmoms_sub[iZbil][imom]-=g2tilde*pr_bil_a2(act,mom,L,iZbil);
  //   }
  
  const index_t r_r_iZbil_ind_imom_ind=r_r_iZbil_ind*index_t({{"ind_mom",equiv_imoms.size()}});
  const index_t im_r_im_r_iZbil_ind_imom_ind=im_r_im_r_iZbil_ind*index_t({{"ind_mom",equiv_imoms.size()}});
  
  //average equiv moms
  const djvec_t Zq=average_equiv_moms(Zq_allmoms,im_r_ind_imom_ind,im_r_imom_ind);
  const djvec_t Zq_sig1=average_equiv_moms(Zq_sig1_allmoms,im_r_ind_imom_ind,im_r_imom_ind);
  djvec_t Zbil=average_equiv_moms(Zbil_allmoms,im_r_im_r_iZbil_ind_imom_ind,im_r_im_r_iZbil_imom_ind);
  djvec_t Zbil_QED=use_QED?average_equiv_moms(Zbil_QED_allmoms,im_r_im_r_iZbil_ind_imom_ind,im_r_im_r_iZbil_imom_ind):djvec_t();
  
  //chirally extrapolated ones
  const djvec_t Zq_chir=average_equiv_moms(Zq_chir_allmoms,r_ind_imom_ind,r_imom_ind);
  const djvec_t Zq_sig1_chir=average_equiv_moms(Zq_sig1_chir_allmoms,r_ind_imom_ind,r_imom_ind);
  djvec_t Zbil_chir=average_equiv_moms(Zbil_chir_allmoms,r_r_iZbil_ind_imom_ind,r_r_iZbil_imom_ind);
  djvec_t Zbil_QED_chir=use_QED?average_equiv_moms(Zbil_QED_chir_allmoms,r_r_iZbil_ind_imom_ind,r_r_iZbil_imom_ind):djvec_t();
  
  // djvec_t Zq_sub=average_equiv_moms(Zq_allmoms_sub);
  // djvec_t Zq_sig1=average_equiv_moms(Zq_sig1_allmoms);
  // djvec_t sig3=average_equiv_moms(sig3_allmoms);
  // djvec_t sig3_em=average_equiv_moms(sig3_em_allmoms);
  // djvec_t Zq_sig1_em=average_equiv_moms(Zq_sig1_em_allmoms);
  // djvec_t Zq_sig1_sub=average_equiv_moms(Zq_sig1_allmoms_sub);
  // vector<djvec_t> Zbil(nZbil);
  // vector<djvec_t> Zbil_QED_allmoms(nZbil),Zbil_QED(nZbil);
  // vector<djvec_t> Zbil_sub(nZbil);
  // for(size_t iZbil=0;iZbil<nZbil;iZbil++)
  //   {
  //     Zbil[iZbil]=average_equiv_moms(Zq_allmoms/pr_bil_allmoms[iZbil]);
  //     Zbil_QED_allmoms[iZbil]=(pr_bil_a_allmoms[iZbil]+pr_bil_b_allmoms[iZbil]-pr_bil_em_allmoms[iZbil])/pr_bil_allmoms[iZbil]+Zq_sig1_em_allmoms/Zq_sig1_allmoms;
  //     Zbil_sub[iZbil]=average_equiv_moms(Zq_allmoms_sub/pr_bil_allmoms_sub[iZbil]);
  //     Zbil_QED[iZbil]=average_equiv_moms(Zbil_QED_allmoms[iZbil]);
  //   }
  
  // write_Z("Zq_sig1_allmoms",Zq_sig1_allmoms,get_pt2());
  // write_Z("sig3_allmoms",sig3_allmoms,get_pt2());
  // write_Z("sig3_em_allmoms",sig3_em_allmoms,get_pt2());
  
  // write_Z("Zq",Zq,get_indep_pt2());
  // write_Z("Zq_sub",Zq_sub,get_indep_pt2());
  // write_Z("Zq_sig1",Zq_sig1,get_indep_pt2());
  // write_Z("sig3",sig3,get_indep_pt2());
  // write_Z("sig3_em",sig3_em,get_indep_pt2());
  // write_Z("Zq_sig1_em",Zq_sig1_em,get_indep_pt2());
  // write_Z("Zq_sig1_sub",Zq_sig1_sub,get_indep_pt2());
  // for(size_t iZbil=0;iZbil<nZbil;iZbil++)
  //   {
  //     write_Z(combine("Z%c",Zbil_tag[iZbil]),Zbil[iZbil],get_indep_pt2());
  //     write_Z(combine("Z%c_QED_allmoms",Zbil_tag[iZbil]),Zbil_QED_allmoms[iZbil],get_pt2());
  //     write_Z(combine("Z%c_QED",Zbil_tag[iZbil]),Zbil_QED[iZbil],get_indep_pt2());
  //     write_Z(combine("Z%c_sub",Zbil_tag[iZbil]),Zbil_sub[iZbil],get_indep_pt2());
  //   }
  
  // write_Z("Zq_sub",Zq_sub,get_indep_pt2());
  // write_Z("Zq_sig1_sub",Zq_sig1_sub,get_indep_pt2());
  
  // linfit_Z(Zbil_QED[iZV]*sqr(4*M_PI),"ZV_QED",-20.6178);
  // linfit_Z(Zbil_QED[iZA]*sqr(4*M_PI),"ZA_QED",-15.7963);

  // //! list of p2tile
  // vector<double> p2tilde(equiv_imoms.size());
  // for(size_t ind_imom=0;ind_imom<equiv_imoms.size();ind_imom++)
  //   {
  //     size_t imom=equiv_imoms[ind_imom].first;
  //     p2tilde[imom]=imoms[imom].p(L).tilde().norm2();
  //   }

  for(auto &p : Zq_tasks)
    {
      const djvec_t &Zq=(*get<0>(p));
      djvec_t &Zq_chir=(*get<1>(p));
      const string &tag=get<2>(p);
      
      grace_file_t out("plots/"+tag+".xmg");
      out.new_data_set();
      
      //m and r
      for(size_t im_r=0;im_r<im_r_ind.max();im_r++)
	{
	  vector<size_t> im_r_comps=im_r_ind(im_r);
	  size_t im=im_r_comps[0],r=im_r_comps[1];
	  for(size_t ind_imom=0;ind_imom<equiv_imoms.size();ind_imom++)
	    {
	      size_t im_r_imom=im_r_ind_imom_ind({im,r,ind_imom});
	      size_t imom=equiv_imoms[ind_imom].first;
	      out.write_ave_err(imoms[imom].p(L).tilde().norm2(),Zq[im_r_imom].ave_err());
	    }
	  out.new_data_set();
	}
      
      //chir extr
      for(size_t r=0;r<nr;r++)
	for(size_t ind_imom=0;ind_imom<equiv_imoms.size();ind_imom++)
	  {
	    size_t r_imom=r_ind_imom_ind({r,ind_imom});
	    size_t imom=equiv_imoms[ind_imom].first;
	    out.write_ave_err(imoms[imom].p(L).tilde().norm2(),Zq_chir[r_imom].ave_err());
	  }
      out.new_data_set();
    }
  
  //! list of task to chiral extrapolate bilinears
  vector<Z_task_t> Zbil_tasks{{&Zbil,&Zbil_chir,string("Zbil")}};
  if(use_QED) Zbil_tasks.push_back(make_tuple(&Zbil_QED,&Zbil_QED_chir,"Zbil_EM"));
  
  for(auto &p : Zq_tasks)
    {
      const djvec_t &Z=(*get<0>(p));
      const djvec_t &Z_chir=(*get<1>(p));
      const string &tag=get<2>(p);
      
      for(size_t r_r_iZbil=0;r_r_iZbil<r_r_iZbil_ind.max();r_r_iZbil++)
	{
	  const vector<size_t> r_r_iZbil_comp=r_r_iZbil_ind(r_r_iZbil);
	  const size_t r1=r_r_iZbil_comp[0],r2=r_r_iZbil_comp[1];
	  const size_t iZbil=r_r_iZbil_comp[2];
	  grace_file_t out("plots/"+tag+"_Z"+Zbil_tag[iZbil]+"_r1_"+to_string(r1)+"_r2_"+to_string(r2)+".xmg");
	  out.new_data_set();
	  
	  for(size_t im1=0;im1<nm;im1++)
	    for(size_t im2=im1;im2<nm;im2++)
	      {
		for(size_t ind_imom=0;ind_imom<equiv_imoms.size();ind_imom++)
		  out.write_ave_err(imoms[equiv_imoms[ind_imom].first].p(L).tilde().norm2(),Z[im_r_im_r_iZbil_ind_imom_ind({im1,r1,im2,r2,iZbil,ind_imom})].ave_err());
		out.new_data_set();
	      }
	  for(size_t ind_imom=0;ind_imom<equiv_imoms.size();ind_imom++)
	    out.write_ave_err(imoms[equiv_imoms[ind_imom].first].p(L).tilde().norm2(),Z_chir[r_r_iZbil_ind_imom_ind({r1,r2,iZbil,ind_imom})].ave_err());
	  out.new_data_set();
	}
    }
  
  cout<<ts<<endl;
  
  return 0;
}
