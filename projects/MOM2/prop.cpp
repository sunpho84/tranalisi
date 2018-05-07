#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_PROP
 #include <MOM2/prop.hpp>

#include <MOM2/pars.hpp>
#include <MOM2/perens.hpp>
#include <MOM2/timings.hpp>

namespace qprop
{
  void set_ins()
  {
    if(pars::use_QED)
      {
	ins_list={LO , FF , F , T , S , P};
	ins_tag={"0" ,"FF","F","T","S","P"};
      }
    else
      {
	ins_list={LO};
	ins_tag={"0"};
      }
    
    nins=ins_list.size();
  }
}

namespace jqprop
{
  void set_ins()
  {
    if(pars::use_QED)
      {
	ins_list={LO , PH , CR , TM};
	ins_tag={"LO","PH","CR","TM"};
      }
    else
      {
	ins_list={LO};
	ins_tag={"LO"};
      }
    
    nins=ins_list.size();
  }
}

namespace lprop
{
  void set_ins()
  {
    if(pars::use_QED)
      {
	ins_list={LO , F };
	ins_tag={"0" ,"F"};
      }
    else
      {
	ins_list={LO};
	ins_tag={"0"};
      }
    
    nins=ins_list.size();
  }
}

dcompl_t coeff_to_read(const qprop::ins ins,const size_t r)
{
  switch(ins)
    {
    case qprop::P:
      return dcompl_t(0,tau3[r]);
      break;
    case qprop::S:
      return -1.0;
      break;
    default:
      return 1.0;
    }
}

string get_qprop_filename(const size_t im,const size_t ir,const qprop::ins ins)
{
  return combine("S_M%zu_R%zu_%s",im,ir,qprop::ins_tag[ins].c_str());
}

string get_lprop_filename(const lprop::ins ins)
{
  return combine("L_%s",lprop::ins_tag[ins].c_str());
}

vector<raw_file_t> perens_t::setup_read_all_qprops_mom(const vector<size_t> &conf_list) const
{
  using namespace pars;
  
  string suff_hit="";
  if(nhits>1) suff_hit="_hit_%zu";
  
  vector<raw_file_t> files(im_r_iconf_ihit_iqins_ind.max());
  
#pragma omp parallel for
  for(size_t i=0;i<im_r_iconf_ihit_iqins_ind.max();i++)
    {
      const vector<size_t> comps=im_r_iconf_ihit_iqins_ind(i);
      const size_t im=comps[0];
      const size_t r=comps[1];
      const size_t iconf=comps[2];
      const size_t ihit=comps[3];
      const qprop::ins ins=qprop::ins_list[comps[4]];
      const string path_base=combine("%s/%s/%04zu/fft_",dir_path.c_str(),prop_hadr_path.c_str(),conf_list[iconf]);
      const string path_suff=combine(suff_hit.c_str(),ihit);
      const string path=path_base+get_qprop_filename(im,r,ins)+path_suff;
      
      files[i].open(path,"r");
    }
  
  return files;
}

vector<raw_file_t> perens_t::setup_read_all_lprops_mom(const vector<size_t> &conf_list) const
{
  vector<raw_file_t> files(iconf_ihit_ilins_ind.max());
  
  string suff_hit="";
  if(nhits>1) suff_hit="_hit_%zu";
  
#pragma omp parallel for
  for(size_t i=0;i<iconf_ihit_ilins_ind.max();i++)
    {
      const vector<size_t> comps=iconf_ihit_ilins_ind(i);
      const size_t iconf=comps[0];
      const size_t ihit=comps[1];
      const lprop::ins ins=lprop::ins_list[comps[2]];
      const string path_base=combine("%s/%s/%04zu/fft_",dir_path.c_str(),prop_lep_path.c_str(),conf_list[iconf]);
      const string path_suff=combine(suff_hit.c_str(),ihit);
      
      files[i].open(path_base+get_lprop_filename(ins)+path_suff,"r");
    }
  
  return files;
}

template <typename T>
T get_rotator(const vector<T> &gamma,const int r)
{
  switch(r)
    {
    case 0:
    case 1:
      return (gamma[0]+dcompl_t(0,1)*tau3[r]*gamma[5])/sqrt(2);
      break;
    default:
      return gamma[0];
    }
}

qprop_t read_qprop(raw_file_t &file,const dcompl_t &fact,const size_t imom,const int r_si,const int r_so)
{
  //cout<<"Seeking file "<<file.get_path()<<" to mom "<<imom<<" position "<<imom*sizeof(dcompl_t)*NSPIN*NSPIN*NCOL*NCOL<<" from "<<file.get_pos()<<endl;
  file.set_pos(imom*sizeof(dcompl_t)*NSPIN*NSPIN*NCOL*NCOL);
  
  qprop_t temp;
  for(size_t is_so=0;is_so<NSPIN;is_so++)
    for(size_t ic_so=0;ic_so<NCOL;ic_so++)
      for(size_t is_si=0;is_si<NSPIN;is_si++)
	for(size_t ic_si=0;ic_si<NCOL;ic_si++)
	  {
	    dcompl_t c;
	    file.bin_read(c);
	    temp(isc(is_si,ic_si),isc(is_so,ic_so))=c*fact;
	  }
  
  if(pars::twisted_run and pars::phys_basis)
    {
      auto rot_si=get_rotator(quaGamma,r_si);
      auto rot_so=get_rotator(quaGamma,r_so);
      
      return rot_si*temp*rot_so;
    }
  else
    return temp;
}

lprop_t read_lprop(raw_file_t &file,const dcompl_t &fact,const size_t imom,const int r_si,const int r_so)
{
  //cout<<"Seeking file "<<file.get_path()<<" to mom "<<imom<<" position "<<imom*sizeof(dcompl_t)*NSPIN*NSPIN*NCOL*NCOL<<" from "<<file.get_pos()<<endl;
  file.set_pos(imom*sizeof(dcompl_t)*NSPIN*NSPIN*NCOL*NCOL);
  
  lprop_t temp;
  for(size_t is_so=0;is_so<NSPIN;is_so++)
    for(size_t ic_so=0;ic_so<NCOL;ic_so++)
      for(size_t is_si=0;is_si<NSPIN;is_si++)
	for(size_t ic_si=0;ic_si<NCOL;ic_si++)
	  {
	    dcompl_t c;
	    file.bin_read(c);
	    temp(is_si,is_so)=c*fact;
	  }
  
  if(pars::twisted_run and pars::phys_basis)
    {
      auto rot_si=get_rotator(lepGamma,r_si);
      auto rot_so=get_rotator(lepGamma,r_so);
      auto out=rot_si*temp*rot_so;
      
      return out;
    }
  else
    return temp;
}

vector<qprop_t> perens_t::read_all_qprops_mom(vector<raw_file_t> &files,const size_t i_in_clust_ihit,const size_t imom)
{
  //! output
  vector<qprop_t> props(im_r_iqins_ijack_ind.max());
  
  //decompose the outer index
  const vector<size_t> i_in_clust_ihit_comp=i_in_clust_ihit_ind(i_in_clust_ihit);
  const size_t i_in_clust=i_in_clust_ihit_comp[0],ihit=i_in_clust_ihit_comp[1];
  
#pragma omp parallel for
  for(size_t im_r_iqins_ijack=0;im_r_iqins_ijack<im_r_iqins_ijack_ind.max();im_r_iqins_ijack++)
    {
      const vector<size_t> comps=im_r_iqins_ijack_ind(im_r_iqins_ijack);
      const size_t im=comps[0];
      const size_t r=comps[1];
      const size_t iqins=comps[2];
      const size_t ijack=comps[3];
      
      //! index of the conf built from ijack and i_in_clust
      const size_t iconf=i_in_clust+clust_size*ijack;
      
      //! index of the file to use
      const size_t im_r_iconf_ihit_iqins=im_r_iconf_ihit_iqins_ind({im,r,iconf,ihit,iqins});
      
      props[im_r_iqins_ijack]=read_qprop(files[im_r_iconf_ihit_iqins],coeff_to_read(qprop::ins_list[iqins],r),imom,r,r);
    }
  
  return props;
}

vector<lprop_t> perens_t::read_all_lprops_mom(vector<raw_file_t> &files,const size_t i_in_clust_ihit,const size_t imom)
{
  //! output
  vector<lprop_t> props(ilins_ijack_ind.max());
  
  //decompose the outer index
  const vector<size_t> i_in_clust_ihit_comp=i_in_clust_ihit_ind(i_in_clust_ihit);
  const size_t i_in_clust=i_in_clust_ihit_comp[0],ihit=i_in_clust_ihit_comp[1];
  
#pragma omp parallel for
  for(size_t ilins_ijack=0;ilins_ijack<ilins_ijack_ind.max();ilins_ijack++)
    {
      const vector<size_t> ilins_ijack_comps=ilins_ijack_ind(ilins_ijack);
      const size_t ilins=ilins_ijack_comps[0];
      const size_t ijack=ilins_ijack_comps[1];
      
      //! index of the conf built from ijack and i_in_clust
      const size_t iconf=i_in_clust+clust_size*ijack;
      
      //! index of the file to use
      const size_t iconf_ihit_ilins=iconf_ihit_ilins_ind({iconf,ihit,ilins});
      
      props[ilins_ijack]=read_lprop(files[iconf_ihit_ilins],1.0,imom,0,0); //r of lprop is always 0
      
      cout<<ilins<<" "<<ijack<<endl;
      cout<<"/////////////////////////////////////////////////////////////////"<<endl;
      cout<<&props[ilins_ijack](0,0).real()<<" "<<props[ilins_ijack]<<endl;
    }
  
  return props;
}

void perens_t::build_all_mr_jackkniffed_qprops(vector<jqprop_t>& jprops,const vector<qprop_t>& props,const size_t imom) const
{
  //! list of all combination of transformations to be applied
  vector<tuple<jqprop::ins,qprop::ins,int>> map;
  
#define ADD_COMBO(JQ,Q,SIGN) map.push_back({jqprop::JQ,qprop::Q,SIGN})
  ADD_COMBO(LO,LO,+1);
  if(pars::use_QED)
    {
      ADD_COMBO(PH,FF,+1);
      ADD_COMBO(PH,T,+1);
      ADD_COMBO(CR,P,+1);
      ADD_COMBO(TM,S,+1);
    }
#undef ADD_COMBO
  
#pragma omp parallel for
  for(size_t im_r_ijack=0;im_r_ijack<im_r_ijack_ind.max();im_r_ijack++)
    {
      const vector<size_t> im_r_ijack_comps=im_r_ijack_ind(im_r_ijack);
      const size_t im=im_r_ijack_comps[0],r=im_r_ijack_comps[1],ijack=im_r_ijack_comps[2];
      
      for(auto t : map)
	{
	  jqprop::ins jins=get<jqprop::ins>(t);
	  qprop::ins ins=get<qprop::ins>(t);
	  int sign=get<int>(t);
	  
	  const size_t im_r_ijqins=im_r_ijqins_ind({im,r,jins});
	  const size_t im_r_iqins_ijack=im_r_iqins_ijack_ind({im,r,ins,ijack});
	  qprop_t &j=jprops[im_r_ijqins][ijack];
	  const qprop_t &p=props[im_r_iqins_ijack];
	  
	  j+=sign*p;
	}
    }
}

vector<jqprop_t> perens_t::get_inverse_propagators(const vector<jqprop_t>& jqprops) const
{
  invert_time.start();
  
  vector<jqprop_t> jqprops_inv(im_r_ijqins_ind.max());
  
#pragma omp parallel for
  for(size_t im_r_ijack=0;im_r_ijack<im_r_ijackp1_ind.max();im_r_ijack++)
    {
      //decript indices
      const vector<size_t> im_r_ijack_comps=im_r_ijackp1_ind(im_r_ijack);
      const size_t im=im_r_ijack_comps[0];
      const size_t r=im_r_ijack_comps[1];
      const size_t ijack=im_r_ijack_comps[2];
      
      //compute inverse
      const size_t im_r_LO=im_r_ijqins_ind({im,r,jqprop::LO});
      qprop_t prop_inv=jqprops_inv[im_r_LO][ijack]=jqprops[im_r_LO][ijack].inverse();
      
      //do the same with QED
      if(pars::use_QED)
	for(jqprop::ins ijqins : {jqprop::CR,jqprop::TM,jqprop::PH})
	  {
	    const size_t im_r_ijqins=im_r_ijqins_ind({im,r,ijqins});
	    jqprops_inv[im_r_ijqins][ijack]=-prop_inv*jqprops[im_r_ijqins][ijack]*prop_inv;
	  }
    }
  
  invert_time.stop();
  
  return jqprops_inv;
}

// void perens_t::clusterize_all_mr_jackkniffed_qprops(vector<jm_r_mom_qprops_t> &jprops) const
// {
// #pragma omp parallel for
//   for(size_t iprop=0;iprop<jprops.size();iprop++)
//     jprops[iprop].clusterize_all_mr_props(pars::use_QED,clust_size);
// }
