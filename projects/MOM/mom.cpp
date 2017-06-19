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

//! write a given Z
void write_Z(const string &name,const djvec_t &Z,const vector<double> &pt2)
{
  grace_file_t outf("plots/"+name+".xmg");
  outf<<fixed;
  outf.precision(8);
  outf.write_vec_ave_err(pt2,Z.ave_err());
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
  
  const string mom_list_path=input.read<string>("MomList"); //!< list of momenta
  const size_t ext_njacks=input.read<size_t>("NJacks"); //!< number of jacks
  
  //! conf range
  range_t conf_range;
  input.expect("ConfRange");
  conf_range.start=input.read<size_t>();
  conf_range.each=input.read<size_t>();
  conf_range.end=input.read<size_t>();
  
  const size_t nhits=input.read<size_t>("NHits"); //!< number of hits
  
  //! template path
  string template_path=input.read<string>("TemplatePath");
  
  //! include or not tadpole
  const double use_tad=input.read<double>("UseTad");
  
  //////////////////////////////////////////////////
  
  //set the number of jackknives
  set_njacks(ext_njacks);
  
  //! sufffix if a number of hits is different from 1
  string suff_hit="";
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
  
  string test_path=template_path;
  if(nhits>1) test_path+="_hit_0";
  vector<size_t> conf_list=get_existing_paths_in_range(test_path,conf_range); //!< list of existing confs
  if(conf_list.size()==0) CRASH("list of configurations is empty! check %s ",test_path.c_str());
  
  //compute deltam_cr
  size_t tmin=12,tmax=23;
  djack_t deltam_cr=compute_deltam_cr(conf_list,tmin,tmax,use_tad);
  cout<<"Deltam cr: "<<deltam_cr<<endl;
  
  size_t clust_size=trim_to_njacks_multiple(conf_list,true); //!< cluster size
  vjprop_t jprop(imoms.size()); //!< jackkniffed LO propagator
  vector<jverts_t> jverts(imoms.size()); //!< jackkniffed vertex
  
  vector<jverts_t> jverts_em(imoms.size()); //!< jackkniffed em vertex
  vjprop_t jprop_em(imoms.size()); //!< jackkniffed FF+T+P propagator
  
  //scope for 2T and P
  {
    vjprop_t jprop_2T(imoms.size()); //!< jackkniffed FF+T propagator
    vjprop_t jprop_P(imoms.size()); //!< jackkniffed FF+T propagator
    vector<jverts_t> jverts_2T(imoms.size()); //!< jackkniffed FF + T vertex
    vector<jverts_t> jverts_P(imoms.size()); //!< jackkniffed P vertex
    
#pragma omp parallel for
    for(size_t ijack=0;ijack<njacks;ijack++)
      for(size_t iconf=ijack*clust_size;iconf<(ijack+1)*clust_size;iconf++)
	for(size_t ihit=0;ihit<nhits;ihit++)
	  {
	    printf("Thread %d/%d reading conf %zu/%zu hit %zu/%zu\n",omp_get_thread_num()+1,omp_get_num_threads(),iconf+1,conf_list.size(),ihit+1,nhits);
	    
	    vprop_t prop=read_prop(template_path+""+suff_hit,conf_list[iconf],ihit);
	    
	    vprop_t prop_FF=read_prop(template_path+"_FF"+suff_hit,conf_list[iconf],ihit);
	    vprop_t prop_F=read_prop(template_path+"_F"+suff_hit,conf_list[iconf],ihit);
	    vprop_t prop_T=read_prop(template_path+"_T"+suff_hit,conf_list[iconf],ihit);
	    vprop_t prop_P=read_prop(template_path+"_P"+suff_hit,conf_list[iconf],ihit);
	    for(auto &pP : prop_P) pP*=dcompl_t(0.0,-1.0);
	    
	    vprop_t prop_2T=prop_FF+prop_T*use_tad;
	    
	    build_jackknifed_prop(jprop,prop,ijack);
	    build_jackknifed_prop(jprop_2T,prop_2T,ijack);
	    build_jackknifed_prop(jprop_P,prop_P,ijack);
	    
	    build_jackknifed_verts(jverts,prop,prop,ijack);
	    
	    build_jackknifed_verts(jverts_2T,prop_F,prop_F,ijack);
	    if(use_tad)
	      {
		build_jackknifed_verts(jverts_2T,prop_2T,prop,ijack);
		build_jackknifed_verts(jverts_2T,prop,prop_2T,ijack);
	      }
	    
	    build_jackknifed_verts(jverts_P,prop_P,prop,ijack);
	    build_jackknifed_verts(jverts_P,prop,prop_P,ijack);
	  }
    
    //clusterize
    clusterize_prop(jprop,clust_size*nhits);
    clusterize_prop(jprop_2T,clust_size*nhits);
    clusterize_verts(jverts,clust_size*nhits);
    clusterize_verts(jverts_2T,clust_size*nhits);
    clusterize_verts(jverts_P,clust_size*nhits);
    
    for(size_t imom=0;imom<imoms.size();imom++)
      {
	for(size_t iGamma=0;iGamma<nGamma;iGamma++) jverts_em[imom][iGamma]=jverts_2T[imom][iGamma]-deltam_cr*jverts_P[imom][iGamma];
	jprop_em[imom]=jprop_2T[imom]-deltam_cr*jprop_P[imom];
      }
  }
  
  vjprop_t jprop_inv=get_prop_inv(jprop); //!< inverse prop
  vjprop_t jprop_em_inv=jprop_inv*jprop_em*jprop_inv;
  
  //compute Zq, Zq_sig1 and Zbil for all moms
  djvec_t Zq_allmoms=compute_Zq(jprop_inv);
  djvec_t Zq_sig1_allmoms=compute_Zq_sig1(jprop_inv);
  djvec_t Zq_sig1_em_allmoms=compute_Zq_sig1(jprop_em_inv);
  vector<djvec_t> pr_bil_allmoms=compute_proj_bil(jprop_inv,jverts,jprop_inv);
  
  //QED
  vector<djvec_t> pr_bil_em_allmoms=compute_proj_bil(jprop_inv,jverts_em,jprop_inv);
  vector<djvec_t> pr_bil_a_allmoms=compute_proj_bil(jprop_em_inv,jverts,jprop_inv);
  vector<djvec_t> pr_bil_b_allmoms=compute_proj_bil(jprop_inv,jverts,jprop_em_inv);
  
  //correct Z LO
  djvec_t Zq_allmoms_sub=Zq_allmoms;
  djvec_t Zq_sig1_allmoms_sub=Zq_sig1_allmoms;
  vector<djvec_t> pr_bil_allmoms_sub=pr_bil_allmoms;
  for(size_t imom=0;imom<imoms.size();imom++)
    {
      imom_t mom=imoms[imom];
      Zq_allmoms_sub[imom]-=g2tilde*sig1_a2(act,mom,L);
      Zq_sig1_allmoms_sub[imom]-=g2tilde*sig1_a2(act,mom,L);
      for(size_t iZbil=0;iZbil<nZbil;iZbil++)
	pr_bil_allmoms_sub[iZbil][imom]-=g2tilde*pr_bil_a2(act,mom,L,iZbil);
    }
  
  //average equiv moms
  djvec_t Zq=average_equiv_moms(Zq_allmoms);
  djvec_t Zq_sub=average_equiv_moms(Zq_allmoms_sub);
  djvec_t Zq_sig1=average_equiv_moms(Zq_sig1_allmoms);
  djvec_t Zq_sig1_em=average_equiv_moms(Zq_sig1_em_allmoms);
  djvec_t Zq_sig1_sub=average_equiv_moms(Zq_sig1_allmoms_sub);
  vector<djvec_t> Zbil(nZbil);
  vector<djvec_t> Zbil_QED_allmoms(nZbil),Zbil_QED(nZbil);
  vector<djvec_t> Zbil_sub(nZbil);
  for(size_t iZbil=0;iZbil<nZbil;iZbil++)
    {
      Zbil[iZbil]=average_equiv_moms(Zq_allmoms/pr_bil_allmoms[iZbil]);
      Zbil_QED_allmoms[iZbil]=(pr_bil_a_allmoms[iZbil]+pr_bil_b_allmoms[iZbil]-pr_bil_em_allmoms[iZbil])/pr_bil_allmoms[iZbil]+Zq_sig1_em_allmoms/Zq_sig1_allmoms;
      Zbil_sub[iZbil]=average_equiv_moms(Zq_allmoms_sub/pr_bil_allmoms_sub[iZbil]);
      Zbil_QED[iZbil]=average_equiv_moms(Zbil_QED_allmoms[iZbil]);
    }
  
  write_Z("Zq_sig1_allmoms",Zq_sig1_allmoms,get_pt2());
  
  write_Z("Zq",Zq,get_indep_pt2());
  write_Z("Zq_sub",Zq_sub,get_indep_pt2());
  write_Z("Zq_sig1",Zq_sig1,get_indep_pt2());
  write_Z("Zq_sig1_em",Zq_sig1_em,get_indep_pt2());
  write_Z("Zq_sig1_sub",Zq_sig1_sub,get_indep_pt2());
  for(size_t iZbil=0;iZbil<nZbil;iZbil++)
    {
      write_Z(combine("Z%c",Zbil_tag[iZbil]),Zbil[iZbil],get_indep_pt2());
      write_Z(combine("Z%c_QED_allmoms",Zbil_tag[iZbil]),Zbil_QED_allmoms[iZbil],get_pt2());
      write_Z(combine("Z%c_QED",Zbil_tag[iZbil]),Zbil_QED[iZbil],get_indep_pt2());
      write_Z(combine("Z%c_sub",Zbil_tag[iZbil]),Zbil_sub[iZbil],get_indep_pt2());
    }
  
  write_Z("Zq_sub",Zq_sub,get_indep_pt2());
  write_Z("Zq_sig1_sub",Zq_sig1_sub,get_indep_pt2());
  
  return 0;
}
