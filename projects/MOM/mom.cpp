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

string suff_hit="";

index_t conf_ind; //!< index of a conf given ijack and i_in_clust
index_t im_r_ind; //!< index of im,r
index_t im_r_ijack_ind; //!< index of im,r,ijack combo
index_t im_r_imom_ind; //!< index of im,r,imom combo
index_t im_r_ind_imom_ind; //!< index of im,r,imom combo
index_t i_in_clust_ihit_ind; //!< index of i_in_clust,ihit
index_t im_r_im_r_igam_ind; //!< index for gamma
index_t im_r_im_r_igam_ijack_ind; //!< index for jgamma
index_t im_r_im_r_ibil_ijack_ind; //!< index for jbil

//! prepare a list of reading task, to be executed in parallel
vector<task_list_t> prepare_read_prop_taks(vector<m_r_mom_conf_props_t> &props,const vector<size_t> &conf_list,bool use_QED)
{
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
		list.push_back(tup_in_t(&l.P,"P",dcompl_t(0,-1)));
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
  djack_t deltam_cr=compute_deltam_cr(conf_list,tmin,tmax,im_sea);
  cout<<"Deltam cr: "<<deltam_cr<<endl;
  
  size_t clust_size=trim_to_njacks_multiple(conf_list,true); //!< cluster size
  
  bool use_QED=true;
  
  conf_ind.set_ranges({{"ijack",njacks},{"i_in_clust",clust_size}});
  im_r_ind.set_ranges({{"m",nm},{"r",nr}});
  im_r_ijack_ind.set_ranges({{"m",nm},{"r",nr},{"ijack",njacks}});
  im_r_imom_ind.set_ranges({{"m",nm},{"r",nr},{"imom",imoms.size()}});
  im_r_ind_imom_ind.set_ranges({{"m",nm},{"r",nr},{"ind_mom",equiv_imoms.size()}});
  i_in_clust_ihit_ind.set_ranges({{"i_in_clust",clust_size},{"ihit",nhits_to_use}});
  im_r_im_r_igam_ind.set_ranges({{"im",nm},{"r",nr},{"im",nm},{"r",nr},{"igamma",nGamma}});
  im_r_im_r_igam_ijack_ind.set_ranges({{"im",nm},{"r",nr},{"im",nm},{"r",nr},{"igamma",nGamma},{"ijack",njacks}});
  im_r_im_r_ibil_ijack_ind.set_ranges({{"im",nm},{"r",nr},{"im",nm},{"r",nr},{"ibil",nZbil},{"ijack",njacks}});
  
  //Zq for all moms, with and without em
  djvec_t Zq_allmoms(im_r_imom_ind.max());
  djvec_t Zq_sig1_allmoms(im_r_imom_ind.max());
  djvec_t Zq_sig1_EM_allmoms(im_r_imom_ind.max());
  
  vector<m_r_mom_conf_props_t> props(im_r_ijack_ind.max()); //!< store props for individual conf
  
  vector<task_list_t> read_tasks=prepare_read_prop_taks(props,conf_list,use_QED);
  for(size_t imom=0;imom<imoms.size();imom++)
    {
      vector<jm_r_mom_props_t> jprops(im_r_ind.max()); //!< jackknived props
      jbil_vert_t jverts(im_r_im_r_igam_ind.max(),use_QED); //!< jackknived vertex
      
      for(size_t i_in_clust=0;i_in_clust<clust_size;i_in_clust++)
	for(size_t ihit=0;ihit<nhits_to_use;ihit++)
	  {
	    size_t i_in_clust_hit=i_in_clust_ihit_ind({i_in_clust,ihit});
	    cout<<"Working on clust_entry "<<i_in_clust<<"/"<<clust_size<<", hit "<<ihit<<"/"<<nhits<<", momentum "<<imom+1<<"/"<<imoms.size()<<endl;
	    read_tasks[i_in_clust_hit].assolve_all(RECYCLE);
	    
	    build_all_mr_jackknifed_props(jprops,props,use_QED,im_r_ijack_ind,im_r_ind);
	    build_all_mr_gbil_jackknifed_verts(jverts,props,im_r_im_r_igam_ind,im_r_ijack_ind,use_QED);
	  }
      
      clusterize_all_mr_jackknifed_props(jprops,use_QED,clust_size);
      jverts.clusterize_all(use_QED,clust_size);
      
      finish_jverts_EM(jverts,deltam_cr);
      finish_jprops_EM(jprops,deltam_cr);
      
      vector<jprop_t> jprop_inv(im_r_ind.max()); //!< inverse propagator
      vector<jprop_t> jprop_EM_inv(im_r_ind.max()); //!< inverse propagator with em insertion
#pragma omp parallel for
      for(size_t im_r_ijack=0;im_r_ijack<im_r_ijack_ind.max();im_r_ijack++)
	{
	  //decript indices
	  vector<size_t> im_r_ijack_comps=im_r_ijack_ind(im_r_ijack);
	  size_t im=im_r_ijack_comps[0],r=im_r_ijack_comps[1],ijack=im_r_ijack_comps[2];
	  size_t im_r=im_r_ind({im,r});
	  size_t im_r_imom=im_r_imom_ind({im,r,imom});
	  
	  //compute inverse
	  prop_t prop_inv=get_from_jackknife(jprops[im_r].LO,ijack).inverse();
	  prop_t prop_EM_inv=prop_inv*get_from_jackknife(jprops[im_r].EM,ijack)*prop_inv;
	  
	  //store inverses
	  put_into_jackknife(jprop_inv[im_r],prop_inv,ijack);
	  put_into_jackknife(jprop_EM_inv[im_r],prop_EM_inv,ijack);
	  
	  //compute Zq
	  Zq_allmoms[im_r_imom][ijack]=compute_Zq(prop_inv,imom);
	  Zq_sig1_allmoms[im_r_imom][ijack]=compute_Zq_sig1(prop_inv,imom);
	  
	  Zq_sig1_EM_allmoms[im_r_imom][ijack]=-compute_Zq_sig1(prop_EM_inv,imom);
	}
      
      djvec_t pr_bil_allmoms=compute_proj_bil(jprop_inv,jverts.LO,jprop_inv,im_r_ind);
      //QED
      djvec_t pr_bil_EM_allmoms,pr_bil_a_allmoms,pr_bil_b_allmoms;
      if(use_QED)
	{
	  pr_bil_EM_allmoms=compute_proj_bil(jprop_inv,jverts.EM,jprop_inv,im_r_ind);
	  pr_bil_a_allmoms=compute_proj_bil(jprop_EM_inv,jverts.LO,jprop_inv,im_r_ind);
	  pr_bil_b_allmoms=compute_proj_bil(jprop_inv,jverts.LO,jprop_EM_inv,im_r_ind);
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
  
  //average equiv moms
  djvec_t Zq=average_equiv_moms(Zq_allmoms,im_r_ind_imom_ind,im_r_imom_ind);
  djvec_t Zq_sig1=average_equiv_moms(Zq_sig1_allmoms,im_r_ind_imom_ind,im_r_imom_ind);
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
  
  grace_file_t out("plots/Zq_sig1_new.xmg");
  out.set_settype(grace::XYDY);
  for(size_t im_r=0;im_r<im_r_ind.max();im_r++)
    {
      vector<size_t> im_r_comps=im_r_ind(im_r);
      size_t im=im_r_comps[0],r=im_r_comps[1];
      for(size_t ind_imom=0;ind_imom<equiv_imoms.size();ind_imom++)
	{
	  size_t im_r_imom=im_r_ind_imom_ind({im,r,ind_imom});
	  size_t imom=equiv_imoms[ind_imom].first;
	  out<<imoms[imom].p(L).tilde().norm2()<<" "<<Zq_sig1[im_r_imom]<<endl;
	}
      out.new_data_set();
    }
  
  return 0;
}
