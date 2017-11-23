#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <tranalisi.hpp>

#include <contractions.hpp>
#include <corrections.hpp>
#include <evolutions.hpp>
#include <geometry.hpp>
#include <types.hpp>

#include <prop.hpp>
#include <sig3.hpp>
#include <Zbil.hpp>
#include <Zq.hpp>
#include <Zq_sig1.hpp>

bool use_QED;
size_t print_each_mom;

string suff_hit="";

index_t conf_ind; //!< index of a conf given ijack and i_in_clust
index_t im_r_ind; //!< index of im,r
index_t imom_ind; //!< index of imom (trivial)
index_t r_imom_ind; //!< index of r,imom combo
index_t im_r_imom_ind; //!< index of im,r,imom combo
index_t indep_imom_ind; //!< index of indep imom (trivial)
index_t r_indep_imom_ind; //!< index of r,indep imom combo
index_t im_r_indep_imom_ind; //!< index of im,r,indep imom combo
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
		  //cout<<"Binding "<<im_r_ijack<<"=(im="<<im<<"r"<<r<<"ijack"<<ijack<<") to file "<<path<<endl;
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
  auto act_key=gaz::decr.find(act_str); //!< key in the map of act
  if(act_key==gaz::decr.end()) CRASH("Unable to decript %s",act_str.c_str());
  gaz::type_t act=act_key->second; //!< action
  
  const double beta=input.read<double>("Beta"); //!< beta
  const double plaq=input.read<double>("Plaq"); //!< plaquette
  nm=input.read<double>("Nm");
  am.resize(nm);
  for(size_t im=0;im<nm;im++) am[im]=input.read<double>();
  const double am_max=*max_element(am.begin(),am.end())*1.1;
  const double am_min=*min_element(am.begin(),am.end())/1.1;
  nr=input.read<double>("Nr");
  
  size_t im_sea=input.read<double>("ImSea"); //!< index of sea mass
  cout<<"ImSea: "<<im_sea<<endl;
  
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
  
  use_QED=input.read<bool>("UseQED");
  print_each_mom=input.read<bool>("PrintEachMom");
  
  //////////////////////////////////////////////////
  
  //set the number of jackknives
  set_njacks(ext_njacks);
  
  //! sufffix if a number of hits is different from 1
  if(nhits>1) suff_hit="_hit_%zu";
  
  double g2=6.0/beta; //!< coupling
  double g2tilde=g2/plaq; //!< boosted coupling
  cout<<"g2tilde: "<<g2tilde<<endl;
  
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
      deltam_cr[im]=compute_deltam_cr(conf_list,tmin,tmax,im,nr);
      cout<<"Deltam cr["<<im<<"]: "<<deltam_cr[im]<<endl;
    }
  deltam_cr_time.stop();
  
  size_t clust_size=trim_to_njacks_multiple(conf_list,true); //!< cluster size
  
  conf_ind.set_ranges({{"ijack",njacks},{"i_in_clust",clust_size}});
  im_r_ind.set_ranges({{"m",nm},{"r",nr}});
  imom_ind.set_ranges({{"imom",imoms.size()}});
  r_imom_ind.set_ranges({{"r",nr},{"imom",imoms.size()}});
  im_r_imom_ind.set_ranges({{"m",nm},{"r",nr},{"imom",imoms.size()}});
  indep_imom_ind.set_ranges({{"indep_mom",equiv_imoms.size()}});
  r_indep_imom_ind.set_ranges({{"r",nr},{"indep_mom",equiv_imoms.size()}});
  im_r_indep_imom_ind.set_ranges({{"m",nm},{"r",nr},{"indep_mom",equiv_imoms.size()}});
  i_in_clust_ihit_ind.set_ranges({{"i_in_clust",clust_size},{"ihit",nhits_to_use}});
  const index_t im_r_im_r_igam_ind=im_r_ind*im_r_ind*index_t({{"igamma",nGamma}});
  const index_t r_r_iZbil_ind({{"r",nr},{"r",nr},{"iZbil",nZbil}});
  const index_t im_r_im_r_iZbil_ind=im_r_ind*im_r_ind*index_t({{"iZbil",nZbil}});
  const index_t iZbil_imom_ind({{"iZbil",nZbil},{"imom",imoms.size()}});
  const index_t im_r_im_r_iZbil_imom_ind=im_r_im_r_iZbil_ind*index_t({{"imom",imoms.size()}});
  const index_t im_r_ijack_ind=im_r_ind*index_t({{"ijack",njacks}});
  const index_t im_r_ijackp1_ind=im_r_ind*index_t({{"ijack",njacks+1}});
  
  //Zq for all moms, with and without EM
  djvec_t Zq_allmoms(im_r_imom_ind.max());
  djvec_t Zq_sig1_allmoms(im_r_imom_ind.max());
  djvec_t Zq_sig1_EM_allmoms;
  if(use_QED) Zq_sig1_EM_allmoms.resize(im_r_imom_ind.max());
  
  //Subtracted Zq, with and without EM, all moms, averaged r
  djvec_t Zq_chir_allmoms(imoms.size());
  djvec_t Zq_chir_allmoms_sub(imoms.size());
  djvec_t Zq_sig1_chir_allmoms(imoms.size());
  djvec_t Zq_sig1_chir_allmoms_sub(imoms.size());
  djvec_t Zq_sig1_EM_chir_allmoms;
  djvec_t Zq_sig1_EM_chir_allmoms_sub;
  if(use_QED)
    {
      Zq_sig1_EM_chir_allmoms.resize(imoms.size());
      Zq_sig1_EM_chir_allmoms_sub.resize(imoms.size());
    }
  
  //! list of task to chirally extrapolate Zq
  vector<tuple<djvec_t*,djvec_t*,string>> Zq_chirextr_tasks{
    {&Zq_allmoms,&Zq_chir_allmoms,string("Zq")},
    {&Zq_sig1_allmoms,&Zq_sig1_chir_allmoms,"Zq_sig1"}};
  if(use_QED) Zq_chirextr_tasks.push_back(make_tuple(&Zq_sig1_EM_allmoms,&Zq_sig1_EM_chir_allmoms,"Zq_sig1_EM"));
  
  djvec_t Zbil_allmoms(im_r_im_r_iZbil_imom_ind.max());
  djvec_t Zbil_chir_allmoms_sub(iZbil_imom_ind.max());
  djvec_t Zbil_chir_allmoms(iZbil_imom_ind.max());
  djvec_t Zbil_QED_allmoms(im_r_im_r_iZbil_imom_ind.max());
  djvec_t Zbil_QED_chir_allmoms(iZbil_imom_ind.max());
  djvec_t Zbil_QED_chir_allmoms_sub(iZbil_imom_ind.max());
  
  vector<m_r_mom_conf_props_t> props; //!< store props for individual conf
  
  vector<task_list_t> read_tasks=prepare_read_prop_taks(props,conf_list);
  
  for(size_t imom=0;imom<imoms.size();imom++)
    {
      vector<jm_r_mom_props_t> jprops(im_r_ind.max()); //!< jackknived props
      jbil_vert_t jverts(im_r_im_r_igam_ind.max(),use_QED); //!< jackknived vertex
      
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
	    build_all_mr_jackknifed_props(jprops,props,use_QED,im_r_ind,deltam_cr);
	    build_props_time.stop();
	    
	    build_verts_time.start();
	    build_all_mr_gbil_jackknifed_verts(jverts,props,im_r_im_r_igam_ind,im_r_ijack_ind,use_QED,deltam_cr);
	    build_verts_time.stop();
	  }
      
      //clusterize
      clust_time.start();
      clusterize_all_mr_jackknifed_props(jprops,use_QED,clust_size);
      jverts.clusterize_all(use_QED,clust_size);
      clust_time.stop();
      
      vector<jprop_t> jprop_inv(im_r_ind.max()); //!< inverse propagator
      vector<jprop_t> jprop_EM_inv(im_r_ind.max()); //!< inverse propagator with em insertion
#pragma omp parallel for reduction(+:invert_time,Zq_time)
      for(size_t im_r_ijack=0;im_r_ijack<im_r_ijackp1_ind.max();im_r_ijack++)
	{
	  //decript indices
	  const vector<size_t> im_r_ijack_comps=im_r_ijackp1_ind(im_r_ijack);
	  const size_t im=im_r_ijack_comps[0],r=im_r_ijack_comps[1],ijack=im_r_ijack_comps[2];
	  const size_t im_r=im_r_ind({im,r});
	  const size_t im_r_imom=im_r_imom_ind({im,r,imom});
	  
	  //compute inverse
	  invert_time.start();
	  prop_t prop_inv=jprops[im_r].LO[ijack].inverse();
	  //if(im_r_ijack==0) cout<<jprops[im_r].LO[ijack](0,0)<<endl;
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
	  
	  //if(ijack==0) cout<<"im_r_imom: "<<im_r_imom<<"=(im"<<im<<"r"<<r<<"imom"<<imom<<"): "<<Zq_sig1_allmoms[im_r_imom][0]<<endl;
	}
      
      //extrapolate to chiral limit Zq
      for(auto & p : Zq_chirextr_tasks)
	{
	  djvec_t y(nm);
	  const djvec_t &Zq=(*get<0>(p));
	  djvec_t &Zq_chir=(*get<1>(p));
	  const string &tag=get<2>(p);
	  
	  //slice m
	  double am_max=*max_element(am.begin(),am.end())*1.1;
	  for(size_t im=0;im<nm;im++)
	    {
	      //averages r if both asked
	      y[im]=0.0;
	      for(size_t r=0;r<nr;r++)
		{
		  size_t i=im_r_imom_ind({im,r,imom});
		  y[im]+=Zq[i]/nr;
		  cout<<tag<<"["<<i<<"=(im"<<im<<"r"<<r<<"imom"<<imom<<")], mom "<<imoms[imom].p(L).norm2()<<": "<<Zq[i].ave_err()<<endl;
		}
	      cout<<tag<<"[im"<<im<<"imom"<<imom<<"]: "<<y[im][0]<<" "<<y[im].err()<<endl;
	    }
	  
	  //fit and write the result
	  djvec_t coeffs=poly_fit(am,y,1,am_min,am_max);
	  if(imom%print_each_mom==0)
	    {
	      grace_file_t plot("plots/chir_extr_"+tag+"_mom_"+to_string(imom)+".xmg");
	      write_fit_plot(plot,0,am_max,bind(poly_eval<djvec_t>,coeffs,_1),am,y);
	      plot.write_ave_err(0,coeffs[0].ave_err());
	    }
	  //extrapolated value
	  Zq_chir[imom]=coeffs[0];
	}
      
      //subtract Zq
      imom_t mom=imoms[imom];
      double sub=g2tilde*sig1_a2(act,gf::LANDAU,group::SU3,mom,L);
      Zq_chir_allmoms_sub[imom]=Zq_chir_allmoms[imom]-sub;
      Zq_sig1_chir_allmoms_sub[imom]=Zq_sig1_chir_allmoms[imom]-sub;
      if(use_QED)
	{
	  double sub_EM=1.0*/*coupling*/sig1_a2(gaz::PLAQ,gf::FEYNMAN,group::U1,mom,L);
	  Zq_sig1_EM_chir_allmoms_sub[imom]=Zq_sig1_EM_chir_allmoms[imom]-sub_EM;
	}
      
      proj_time.start();
      djvec_t pr_bil_mom=compute_proj_bil(jprop_inv,jverts.LO,jprop_inv,im_r_ind);
      djvec_t pr_bil_chir_mom(nZbil);
      
      //QED
      djvec_t pr_bil_QED_mom;
      djvec_t pr_bil_QED_chir_mom(nZbil);
      if(use_QED)
	{
	  djvec_t pr_bil_EM_mom=compute_proj_bil(jprop_inv,jverts.EM,jprop_inv,im_r_ind);
	  djvec_t pr_bil_a_mom=compute_proj_bil(jprop_EM_inv,jverts.LO,jprop_inv,im_r_ind);
	  djvec_t pr_bil_b_mom=compute_proj_bil(jprop_inv,jverts.LO,jprop_EM_inv,im_r_ind);
	  pr_bil_QED_mom=pr_bil_a_mom+pr_bil_b_mom-pr_bil_EM_mom;
	  
	}
      proj_time.stop();
      
      //compute subtractions
      djvec_t pr_bil_mom_correction(nZbil);
      djvec_t pr_bil_QED_mom_correction(nZbil);
      for(size_t iZbil=0;iZbil<nZbil;iZbil++)
	{
	  pr_bil_mom_correction[iZbil]=g2tilde*pr_bil_a2(act,gf::LANDAU,group::SU3,mom,L,iZbil);
	  if(use_QED) pr_bil_QED_mom_correction[iZbil]=1.0*pr_bil_a2(gaz::PLAQ,gf::FEYNMAN,group::U1,mom,L,iZbil);
	}
      
      //extrapolate to chiral limit the projected bilinears
      djvec_t pr_bil_chir_mom_sub(nZbil);
      djvec_t pr_bil_QED_chir_mom_sub(nZbil);
      vector<tuple<djvec_t*,djvec_t*,djvec_t*,djvec_t*,string>>
	pr_bil_chirextr_tasks{{&pr_bil_mom,&pr_bil_chir_mom,&pr_bil_chir_mom_sub,&pr_bil_mom_correction,string("pr_bil")}};
      if(use_QED) pr_bil_chirextr_tasks.push_back(make_tuple(&pr_bil_QED_mom,&pr_bil_QED_chir_mom,&pr_bil_QED_chir_mom_sub,&pr_bil_QED_mom_correction,string("pr_bil_QED")));
      
      for(size_t iZbil=0;iZbil<nZbil;iZbil++)
	for(auto & p : pr_bil_chirextr_tasks)
	  {
	    const djvec_t &pr=(*get<0>(p));
	    djvec_t &pr_chir=(*get<1>(p));
	    djvec_t &pr_chir_sub=(*get<2>(p));
	    djvec_t &pr_corr=(*get<3>(p));
	    const string &tag=get<4>(p);
	    
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
		  y_plot[i]=0.0;
		  for(size_t r=0;r<nr;r++) y_plot[i]+=pr[im_r_im_r_iZbil_ind({im1,r,im2,r,iZbil})]/nr;
		  
		  if(sub_pole) y[i]=x[i]*y_plot[i];
		  else         y[i]=y_plot[i];
		  //increment
		  i++;
		}
	    
	    //fit and store extrapolated value
	    djvec_t coeffs=poly_fit(x,y,(sub_pole?2:1),2.0*am_min,2.0*am_max);
	    pr_chir[iZbil]=coeffs[sub_pole?1:0];
	    
	    //subtract from bilinear
	    pr_chir_sub[iZbil]=pr_chir[iZbil]-pr_corr[iZbil];
	    
	    //plot
	    if(imom%print_each_mom==0)
	      {
		grace_file_t plot("plots/chir_extr_"+tag+"_"+Zbil_tag[iZbil]+"_mom_"+to_string(imom)+".xmg");
		write_fit_plot(plot,2*am_min,2*am_max,[&coeffs,sub_pole](double x)->djack_t{return poly_eval<djvec_t>(coeffs,x)/(sub_pole?x:1);},x,y_plot);
		plot.write_ave_err(0.0,pr_chir[iZbil].ave_err());
		plot.new_data_set();
		plot.write_ave_err(0.0,pr_chir_sub[iZbil].ave_err());
	      }
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
	    sqrt(Zq_sig1_allmoms[im_r1_imom]*Zq_sig1_allmoms[im_r2_imom])/pr_bil_mom[im_r_im_r_iZbil];
	  
	  if(use_QED)
	    {
	      Zbil_QED_allmoms[im_r_im_r_iZbil_imom]=
		pr_bil_QED_mom[im_r_im_r_iZbil]/pr_bil_mom[im_r_im_r_iZbil]+
		(Zq_sig1_EM_allmoms[im_r1_imom]/Zq_sig1_allmoms[im_r1_imom]+Zq_sig1_EM_allmoms[im_r2_imom]/Zq_sig1_allmoms[im_r2_imom])/2.0;
	    }
	}
      
      //build Z in the chiral limit
      for(size_t iZbil=0;iZbil<nZbil;iZbil++)
	{
	  const size_t iZbil_imom=iZbil_imom_ind({iZbil,imom});
	  Zbil_chir_allmoms[iZbil_imom]=sqrt(Zq_sig1_chir_allmoms[imom]*Zq_sig1_chir_allmoms[imom])/pr_bil_chir_mom[iZbil];
	  Zbil_chir_allmoms_sub[iZbil_imom]=sqrt(Zq_sig1_chir_allmoms_sub[imom]*Zq_sig1_chir_allmoms_sub[imom])/pr_bil_chir_mom_sub[iZbil];
	  
	  if(use_QED)
	    {
	      Zbil_QED_chir_allmoms[iZbil_imom]=
		pr_bil_QED_chir_mom[iZbil]/pr_bil_chir_mom[iZbil]+
		    (Zq_sig1_EM_chir_allmoms[imom]/Zq_sig1_chir_allmoms[imom]+Zq_sig1_EM_chir_allmoms[imom]/Zq_sig1_chir_allmoms[imom])/2.0;
	      Zbil_QED_chir_allmoms_sub[iZbil_imom]=
		pr_bil_QED_chir_mom_sub[iZbil]/pr_bil_chir_mom_sub[iZbil]+
		    (Zq_sig1_EM_chir_allmoms_sub[imom]/Zq_sig1_chir_allmoms_sub[imom]+Zq_sig1_EM_chir_allmoms_sub[imom]/Zq_sig1_chir_allmoms_sub[imom])/2.0;
	    }
	}
    }
  
  const index_t iZbil_indep_imom_ind({{"Zbil",nZbil},{"indep_mom",equiv_imoms.size()}});
  const index_t im_r_im_r_iZbil_indep_imom_ind=im_r_im_r_iZbil_ind*index_t({{"indep_mom",equiv_imoms.size()}});
  
  //average equiv moms
  djvec_t Zq=average_equiv_moms(Zq_allmoms,im_r_indep_imom_ind,im_r_imom_ind);
  djvec_t Zq_sig1=average_equiv_moms(Zq_sig1_allmoms,im_r_indep_imom_ind,im_r_imom_ind);
  djvec_t Zq_sig1_EM=average_equiv_moms(Zq_sig1_EM_allmoms,im_r_indep_imom_ind,im_r_imom_ind);
  djvec_t Zbil=average_equiv_moms(Zbil_allmoms,im_r_im_r_iZbil_indep_imom_ind,im_r_im_r_iZbil_imom_ind);
  djvec_t Zbil_QED=use_QED?average_equiv_moms(Zbil_QED_allmoms,im_r_im_r_iZbil_indep_imom_ind,im_r_im_r_iZbil_imom_ind):djvec_t();
  
  //chirally extrapolated ones
  djvec_t Zq_chir=average_equiv_moms(Zq_chir_allmoms,indep_imom_ind,imom_ind);
  djvec_t Zq_chir_sub=average_equiv_moms(Zq_chir_allmoms_sub,indep_imom_ind,imom_ind);
  djvec_t Zq_sig1_chir=average_equiv_moms(Zq_sig1_chir_allmoms,indep_imom_ind,imom_ind);
  djvec_t Zq_sig1_EM_chir=average_equiv_moms(Zq_sig1_EM_chir_allmoms,indep_imom_ind,imom_ind);
  djvec_t Zq_sig1_chir_sub=average_equiv_moms(Zq_sig1_chir_allmoms_sub,indep_imom_ind,imom_ind);
  djvec_t Zq_sig1_EM_chir_sub=average_equiv_moms(Zq_sig1_EM_chir_allmoms_sub,indep_imom_ind,imom_ind);
  
  djvec_t Zbil_chir=average_equiv_moms(Zbil_chir_allmoms,iZbil_indep_imom_ind,iZbil_imom_ind);
  djvec_t Zbil_chir_sub=average_equiv_moms(Zbil_chir_allmoms_sub,iZbil_indep_imom_ind,iZbil_imom_ind);
  djvec_t Zbil_QED_chir=use_QED?average_equiv_moms(Zbil_QED_chir_allmoms,iZbil_indep_imom_ind,iZbil_imom_ind):djvec_t();
  djvec_t Zbil_QED_chir_sub=use_QED?average_equiv_moms(Zbil_QED_chir_allmoms_sub,iZbil_indep_imom_ind,iZbil_imom_ind):djvec_t();
  
  using Z_plot_task_t=tuple<djvec_t*,djvec_t*,djvec_t*,string>;
  
  //! list of task to plot chiral extrapolation Zq
  vector<Z_plot_task_t> Zq_plot_tasks{
    {&Zq,&Zq_chir,&Zq_chir_sub,string("Zq")},
    {&Zq_sig1,&Zq_sig1_chir,&Zq_sig1_chir_sub,"Zq_sig1"}};
  if(use_QED) Zq_plot_tasks.push_back(make_tuple(&Zq_sig1_EM,&Zq_sig1_EM_chir,&Zq_sig1_EM_chir_sub,"Zq_sig1_EM"));
  
  for(auto &p : Zq_plot_tasks)
    {
      //decript the tuple
      const djvec_t &Zq=(*get<0>(p));
      djvec_t &Zq_chir=(*get<1>(p));
      djvec_t &Zq_chir_sub=(*get<2>(p));
      const string &tag=get<3>(p);
      
      //loop over all r
      for(size_t r=0;r<nr;r++)
	{
	  //open the file
	  grace_file_t out("plots/"+tag+"_r_"+to_string(r)+".xmg");
	  out.new_data_set();
	  
	  //m
	  for(size_t im=0;im<nm+2;im++)
	    {
	      for(size_t indep_imom=0;indep_imom<equiv_imoms.size();indep_imom++)
		{
		  djack_t y;
		  if(im<nm) y=Zq[im_r_indep_imom_ind({im,r,indep_imom})];
		  if(im==nm+0) y=Zq_chir[indep_imom];
		  if(im==nm+1) y=Zq_chir_sub[indep_imom];
		  
		  out.write_ave_err(imoms[equiv_imoms[indep_imom].first].p(L).tilde().norm2(),y.ave_err());
		}
	      out.new_data_set();
	    }
	}
    }
  
  //! list of task to print the chiral extrapolate bilinears
  vector<Z_plot_task_t> Zbil_tasks{{&Zbil,&Zbil_chir,&Zbil_chir_sub,string("Zbil")}};
  if(use_QED) Zbil_tasks.push_back(make_tuple(&Zbil_QED,&Zbil_QED_chir,&Zbil_QED_chir_sub,"Zbil_EM"));
  
  for(auto &p : Zbil_tasks)
    {
      //decript tuple
      const djvec_t &Z=(*get<0>(p));
      const djvec_t &Z_chir=(*get<1>(p));
      const djvec_t &Z_chir_sub=(*get<2>(p));
      const string &tag=get<3>(p);
      
      for(size_t iZbil=0;iZbil<nZbil;iZbil++)
	{
	  //open the file
	  grace_file_t out("plots/"+tag+"_Z"+Zbil_tag[iZbil]+".xmg");
	  out.new_data_set();
	  
	  //write mass by mass, only half of the combos
	  for(size_t im1=0;im1<nm;im1++)
	    for(size_t im2=im1;im2<nm;im2++)
	      {
		for(size_t indep_imom=0;indep_imom<equiv_imoms.size();indep_imom++)
		  {
		    djack_t y;
		    y=0.0;
		    for(size_t r=0;r<nr;r++) y+=Z[im_r_im_r_iZbil_indep_imom_ind({im1,r,im2,r,iZbil,indep_imom})]/nr;
		    
		    out.write_ave_err(imoms[equiv_imoms[indep_imom].first].p(L).tilde().norm2(),y.ave_err());
		  }
		out.new_data_set();
	      }
	  
	  //write chiral extrap and subtracted
	  for(auto &Ztag : vector<tuple<const djvec_t*,string>>{{&Z_chir,"chir"},{&Z_chir_sub,"sub"}})
	    {
	      out.set_legend(get<1>(Ztag));
	      for(size_t indep_imom=0;indep_imom<equiv_imoms.size();indep_imom++)
		out.write_ave_err(imoms[equiv_imoms[indep_imom].first].p(L).tilde().norm2(),(*get<0>(Ztag))[iZbil_indep_imom_ind({iZbil,indep_imom})].ave_err());
	      out.new_data_set();
	    }
	}
    }
  
  
  //print time statistics
  cout<<ts<<endl;
  
  return 0;
}
