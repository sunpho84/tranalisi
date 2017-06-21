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
vector<string> ins_list={"0"};

//! open all files in battery
map<string,vector<raw_file_t>> open_all_prop_files(const vector<size_t> &conf_list,size_t nhits_to_use,bool use_QED)
{
  //add EM if asked
  if(use_QED) for(auto &ins : {"P","S","T","F","FF"}) ins_list.push_back(ins);
  
  //set indices
  m_r_ind.set_ranges({{"m",nm},{"r",nr}});
  conf_hit_ind.set_ranges({{"conf",conf_list.size()},{"hit",nhits_to_use}});
  m_r_conf_hit_ind.set_ranges({{"m",nm},{"r",nr},{"conf",conf_list.size()},{"hit",nhits_to_use}});
  
  //resize the list of prop for each ins
  map<string,vector<raw_file_t>> prop_files;
  for(auto &ins : ins_list) prop_files[ins].resize(m_r_conf_hit_ind.max());
  
  //open for real
  for(auto &ins : ins_list)
    for(size_t im=0;im<nm;im++)
      for(size_t r=0;r<nr;r++)
	for(size_t iconf=0;iconf<conf_list.size();iconf++)
	  for(size_t ihit=0;ihit<nhits_to_use;ihit++)
	    {
	      string path=combine("out/%04zu/fft_",conf_list[iconf])+get_prop_tag(im,r,ins)+combine(suff_hit.c_str(),ihit);
	      prop_files[ins][m_r_conf_hit_ind({im,r,iconf,ihit})].open(path,"r");
	    }
  
  return prop_files;
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
  
  //////////////////////////////////////////////////
  
  //set the number of jackknives
  set_njacks(ext_njacks);
  
  //! sufffix if a number of hits is different from 1
  if(nhits>1) suff_hit="_hit_%zu";
  
  double g2=6.0/beta; //!< coupling
  double g2tilde=g2/plaq; //!< boosted coupling
  
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
  size_t tmin=12,tmax=23;
  djack_t deltam_cr=compute_deltam_cr(conf_list,tmin,tmax,im_sea);
  cout<<"Deltam cr: "<<deltam_cr<<endl;
  
  size_t clust_size=trim_to_njacks_multiple(conf_list,true); //!< cluster size
  
  bool use_QED=true;
  
  set_mr_ind(nm,nr);
  set_conf_props(use_QED);
  set_jprops(use_QED);
  
  set_mr_gbil_ind(nm,nr);
  set_mr_Zbil_ind(nm,nr);
  set_jbil_verts(use_QED);
  
  map<string,vector<raw_file_t>> prop_files=open_all_prop_files(conf_list,nhits,use_QED);
  
  for(size_t imom=0;imom<imoms.size();imom++)
    {
      cout<<"Reading momentum "<<imom+1<<"/"<<imoms.size()<<endl;
#pragma omp parallel for
      for(size_t ijack=0;ijack<njacks;ijack++)
	for(size_t iconf=ijack*clust_size;iconf<(ijack+1)*clust_size;iconf++)
	  for(size_t ihit=0;ihit<nhits_to_use;ihit++)
	    {
	      printf("Thread %d/%d reading conf %zu/%zu hit %zu/%zu\n",omp_get_thread_num()+1,omp_get_num_threads(),iconf+1,conf_list.size(),ihit+1,nhits_to_use);
	      
	      read_all_mr_props(use_QED,prop_files,iconf,ihit);
	      build_all_mr_jackknifed_props(use_QED,ijack);
	      build_all_mr_gbil_jackknifed_verts(use_QED,ijack);
	    }
      
      clusterize_all_mr_props(use_QED,clust_size);
      clusterize_all_mr_gbil_verts(use_QED,clust_size);
      
      jverts_em=jverts_em-jverts_P*SC(deltam_cr);
      jprop_2=jprop_2-jprop_P*SC(deltam_cr);
      
      vjprop_t jprop_inv=get_all_mr_props_inv(jprop_0); //!< inverse prop
      vjprop_t jprop_em_inv=jprop_inv*jprop_2*jprop_inv;
    }
  
  // //compute Zq, Zq_sig1 and Zbil for all moms
  // djvec_t Zq_allmoms=compute_Zq(jprop_inv);
  // djvec_t Zq_sig1_allmoms=compute_Zq_sig1(jprop_inv);
  // djvec_t sig3_allmoms=compute_sig3(jprop_inv);
  // djvec_t sig3_em_allmoms=compute_sig3(jprop_em_inv);
  // djvec_t Zq_sig1_em_allmoms=-compute_Zq_sig1(jprop_em_inv);
  // vector<djvec_t> pr_bil_allmoms=compute_proj_bil(jprop_inv,jverts,jprop_inv);
  
  // //QED
  // vector<djvec_t> pr_bil_em_allmoms=compute_proj_bil(jprop_inv,jverts_em,jprop_inv);
  // vector<djvec_t> pr_bil_a_allmoms=compute_proj_bil(jprop_em_inv,jverts,jprop_inv);
  // vector<djvec_t> pr_bil_b_allmoms=compute_proj_bil(jprop_inv,jverts,jprop_em_inv);
  
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
  
  // //average equiv moms
  // djvec_t Zq=average_equiv_moms(Zq_allmoms);
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
  
  return 0;
}
