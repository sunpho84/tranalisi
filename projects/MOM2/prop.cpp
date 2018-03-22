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

const char m_r_mom_conf_qprops_t::tag[NPROP_WITH_QED][3]={"0","FF","F","T","S","P"};
const char mom_conf_lprops_t::tag[NPROP_WITH_QED][3]={"0","F"};

string get_qprop_tag(const size_t im,const size_t ir,const size_t ikind)
{
  return combine("S_M%zu_R%zu_%s",im,ir,m_r_mom_conf_qprops_t::tag[ikind]);
}

string get_lprop_tag(const size_t ikind)
{
  return combine("L_%s",mom_conf_lprops_t::tag[ikind]);
}

vector<raw_file_t> perens_t::setup_read_all_qprops_mom(const vector<size_t> &conf_list) const
{
  using namespace pars;
  
  string suff_hit="";
  if(nhits>1) suff_hit="_hit_%zu";
  
  vector<raw_file_t> files(im_r_iconf_ihit_iqkind_ind.max());
  
#pragma omp parallel for
  for(size_t i=0;i<im_r_iconf_ihit_iqkind_ind.max();i++)
    {
      const vector<size_t> comps=im_r_iconf_ihit_iqkind_ind(i);
      const size_t im=comps[0];
      const size_t r=comps[1];
      const size_t iconf=comps[2];
      const size_t ihit=comps[3];
      const size_t ikind=comps[4];
      const string path_base=combine("%s/%s/%04zu/fft_",dir_path.c_str(),prop_hadr_path.c_str(),conf_list[iconf]);
      const string path_suff=combine(suff_hit.c_str(),ihit);
      const string path=path_base+get_qprop_tag(im,r,ikind)+path_suff;
      
      files[i].open(path,"r");
    }
  
  return files;
}

vector<raw_file_t> perens_t::setup_read_all_lprops_mom(const vector<size_t> &conf_list) const
{
  vector<raw_file_t> files(iconf_ihit_ilkind_ind.max());
  
  string suff_hit="";
  if(nhits>1) suff_hit="_hit_%zu";
  
#pragma omp parallel for
  for(size_t i=0;i<iconf_ihit_ilkind_ind.max();i++)
    {
      const vector<size_t> comps=iconf_ihit_ilkind_ind(i);
      const size_t iconf=comps[0];
      const size_t ihit=comps[1];
      const size_t ikind=comps[2];
      const string path_base=combine("%s/%s/%04zu/fft_",dir_path.c_str(),prop_lep_path.c_str(),conf_list[iconf]);
      const string path_suff=combine(suff_hit.c_str(),ihit);
      
      files[i].open(path_base+get_lprop_tag(ikind)+path_suff,"r");
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

void read_qprop(qprop_t *prop,raw_file_t &file,const dcompl_t &fact,const size_t imom,const int r_si,const int r_so)
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
  
  auto rot_si=get_rotator(quaGamma,r_si);
  auto rot_so=get_rotator(quaGamma,r_so);
  *prop=rot_si*temp*rot_so;
}

void read_lprop(lprop_t *prop,raw_file_t &file,const dcompl_t &fact,const size_t imom,const int r_si,const int r_so)
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
  
  auto rot_si=get_rotator(lepGamma,r_si);
  auto rot_so=get_rotator(lepGamma,r_so);
  *prop=rot_si*temp*rot_so;
}

void incorporate_charge(vector<m_r_mom_conf_qprops_t> &props,const double ch)
{
  const double ch2=ch*ch;
  
  for(auto &p : props)
    {
      p.F*=ch;
      p.FF*=ch2;
      p.T*=ch2;
      p.P*=ch2;
      p.S*=ch2;
    }
}

void incorporate_charge(vector<mom_conf_lprops_t> &props,const double ch)
{
  for(auto &p : props)
    p.F*=ch;
}

vector<m_r_mom_conf_qprops_t> perens_t::read_all_qprops_mom(vector<raw_file_t> &files,const size_t i_in_clust_ihit,const size_t imom)
{
  //! output
  vector<m_r_mom_conf_qprops_t> props(im_r_ijack_ind.max());
  
  //decompose the outer index
  const vector<size_t> i_in_clust_ihit_comp=i_in_clust_ihit_ind(i_in_clust_ihit);
  const size_t i_in_clust=i_in_clust_ihit_comp[0],ihit=i_in_clust_ihit_comp[1];
  
  //! index of all that must be read
  const index_t im_r_ijack_ikind_ind=im_r_ijack_ind*index_t({{"ikind",m_r_mom_conf_qprops_t::nprop_kind()}});
#pragma omp parallel for
  for(size_t im_r_ijack_ikind=0;im_r_ijack_ikind<im_r_ijack_ikind_ind.max();im_r_ijack_ikind++)
    {
      const vector<size_t> im_r_ijack_ikind_comps=im_r_ijack_ikind_ind(im_r_ijack_ikind);
      const size_t im=im_r_ijack_ikind_comps[0];
      const size_t r=im_r_ijack_ikind_comps[1];
      const size_t ijack=im_r_ijack_ikind_comps[2];
      const size_t ikind=im_r_ijack_ikind_comps[3];
      
      //! index of the conf built from ijack and i_in_clust
      const size_t iconf=i_in_clust+clust_size*ijack;
      
      //! index of the file to use
      //cout<<im<<" "<<r<<" "<<iconf<<" "<<ihit<<" "<<ikind<<endl;
      const size_t im_r_iconf_ihit_ikind=im_r_iconf_ihit_iqkind_ind({im,r,iconf,ihit,ikind});
      
      //! index of the output propagator
      const size_t im_r_ijack=im_r_ijack_ind({im,r,ijack});
      
      read_qprop(props[im_r_ijack].kind[ikind],files[im_r_iconf_ihit_ikind],m_r_mom_conf_qprops_t::coeff_to_read(ikind,r),imom,r,r);
    }
  
  return props;
}

vector<mom_conf_lprops_t> perens_t::read_all_lprops_mom(vector<raw_file_t> &files,const size_t i_in_clust_ihit,const size_t imom)
{
  //! output
  vector<mom_conf_lprops_t> props(njacks);
  
  //decompose the outer index
  const vector<size_t> i_in_clust_ihit_comp=i_in_clust_ihit_ind(i_in_clust_ihit);
  const size_t i_in_clust=i_in_clust_ihit_comp[0],ihit=i_in_clust_ihit_comp[1];
  
  //! index of all that must be read
  const index_t ijack_ikind_ind=index_t({{"ijack",njacks},{"ikind",mom_conf_lprops_t::nprop_kind()}});
#pragma omp parallel for
  for(size_t ijack_ikind=0;ijack_ikind<ijack_ikind_ind.max();ijack_ikind++)
    {
      const vector<size_t> ijack_ikind_comps=ijack_ikind_ind(ijack_ikind);
      const size_t ijack=ijack_ikind_comps[0];
      const size_t ikind=ijack_ikind_comps[1];
      
      //! index of the conf built from ijack and i_in_clust
      const size_t iconf=i_in_clust+clust_size*ijack;
      
      //! index of the file to use
      const size_t iconf_ihit_ikind=iconf_ihit_ilkind_ind({iconf,ihit,ikind});
      
      read_lprop(props[ijack].kind[ikind],files[iconf_ihit_ikind],1.0,imom,0,0); //r of lprop is always 0
    }
  
  return props;
}

void perens_t::build_all_mr_jackkniffed_qprops(vector<jm_r_mom_qprops_t> &jprops,const vector<m_r_mom_conf_qprops_t> &props) const
{
#pragma omp parallel for
  for(size_t i_im_r_ijack=0;i_im_r_ijack<im_r_ijack_ind.max();i_im_r_ijack++)
    {
      vector<size_t> im_r_ijack=im_r_ijack_ind(i_im_r_ijack);
      size_t im=im_r_ijack[0],r=im_r_ijack[1],ijack=im_r_ijack[2];
      size_t im_r=im_r_ind({im,r});
      //cout<<"Building jack prop im="<<im<<", r="<<r<<", ijack="<<ijack<<endl;
      
      jm_r_mom_qprops_t &j=jprops[im_r];
      const m_r_mom_conf_qprops_t &p=props[i_im_r_ijack];
      
      j.LO[ijack]+=p.LO;
      if(pars::use_QED)
	for(auto &jp_p : vector<tuple<jqprop_t*,const qprop_t*>>({
	      {&j.PH,&p.FF},
	      {&j.PH,&p.T},
	      {&j.CT,&p.P},
	      {&j.S,&p.S}}))
	  (*get<0>(jp_p))[ijack]+=*get<1>(jp_p);
    }
}

void perens_t::get_inverse_propagators(vector<jqprop_t> &jprop_inv,vector<jqprop_t> &jprop_QED_inv,
				       const vector<jm_r_mom_qprops_t> &jprops) const
{
  jprop_inv.resize(im_r_ind.max());
  jprop_QED_inv.resize(im_r_ind.max());
  
#pragma omp parallel for reduction(+:invert_time)
  for(size_t im_r_ijack=0;im_r_ijack<im_r_ijackp1_ind.max();im_r_ijack++)
    {
      //decript indices
      const vector<size_t> im_r_ijack_comps=im_r_ijackp1_ind(im_r_ijack);
      const size_t im=im_r_ijack_comps[0],r=im_r_ijack_comps[1],ijack=im_r_ijack_comps[2];
      const size_t im_r=im_r_ind({im,r});
      
      //compute inverse
      invert_time.start();
      qprop_t prop_inv=jprop_inv[im_r][ijack]=jprops[im_r].LO[ijack].inverse();
      invert_time.stop();
      
      //do the same with QED
      if(pars::use_QED)
	{
	  invert_time.start(); //This misses a sign -1 coming from the original inverse
	  jprop_QED_inv[im_r][ijack]=prop_inv*jprops[im_r].QED[ijack]*prop_inv;
	  invert_time.stop();
	}
    }
}

void perens_t::clusterize_all_mr_jackkniffed_qprops(vector<jm_r_mom_qprops_t> &jprops) const
{
#pragma omp parallel for
  for(size_t iprop=0;iprop<jprops.size();iprop++)
    jprops[iprop].clusterize_all_mr_props(pars::use_QED,clust_size);
  
  if(pars::use_QED)
#pragma omp parallel for
    for(size_t i=0;i<im_r_ind.max();i++)
      {
	// const vector<size_t> comps=im_r_ind(i);
	// const size_t r=comps[1];
	for(size_t ijack=0;ijack<=njacks;ijack++)
	  jprops[i].QED[ijack]=
	    jprops[i].PH[ijack]-// +tau3[r]*
	    deltam_cr[i][ijack]*jprops[i].CT[ijack];
      }
}

void perens_t::mom_compute_qprop()
{
  vector<raw_file_t> files=setup_read_all_qprops_mom(conf_list);
  
  for(size_t ilinmom=0;ilinmom<linmoms.size();ilinmom++)
    {
      const size_t mom=linmoms[ilinmom][0];
      vector<jm_r_mom_qprops_t> jprops(im_r_ind.max()); //!< jackknived props
      
      for(size_t i_in_clust=0;i_in_clust<clust_size;i_in_clust++)
	for(size_t ihit=0;ihit<nhits_to_use;ihit++)
	  {
	    const size_t i_in_clust_ihit=i_in_clust_ihit_ind({i_in_clust,ihit});
	    const size_t mom=linmoms[ilinmom][0];
	    cout<<"Working on qprop, "
	      "clust_entry "<<i_in_clust+1<<"/"<<clust_size<<", "
	      "hit "<<ihit+1<<"/"<<nhits<<", "
	      "momentum "<<ilinmom+1<<"/"<<linmoms.size()<<", "
	      "mom: "<<mom<<endl;
	    read_time.start();
	    const vector<m_r_mom_conf_qprops_t> props=read_all_qprops_mom(files,i_in_clust_ihit,mom);
	    read_time.stop();
	    
	    //build all props
	    build_props_time.start();
	    build_all_mr_jackkniffed_qprops(jprops,props);
	    build_props_time.stop();
	  }
      
      //clusterize
      clust_time.start();
      clusterize_all_mr_jackkniffed_qprops(jprops);
      clust_time.stop();
      
      vector<jqprop_t> jprop_inv; //!< inverse propagator
      vector<jqprop_t> jprop_QED_inv; //!< inverse propagator with em insertion
      
      get_inverse_propagators(jprop_inv,jprop_QED_inv,jprops);
      
#pragma omp parallel for reduction(+:invert_time,Zq_time)
      for(size_t im_r_ijack=0;im_r_ijack<im_r_ijackp1_ind.max();im_r_ijack++)
	{
	  //decript indices
	  const vector<size_t> im_r_ijack_comps=im_r_ijackp1_ind(im_r_ijack);
	  const size_t im=im_r_ijack_comps[0],r=im_r_ijack_comps[1],ijack=im_r_ijack_comps[2];
	  const size_t im_r=im_r_ind({im,r});
	  const size_t im_r_ilinmom=im_r_ilinmom_ind({im,r,ilinmom});
	  
	  //compute Zq
	  Zq_time.start();
	  Zq[im_r_ilinmom][ijack]=compute_Zq(jprop_inv[im_r][ijack],mom);
	  Zq_sig1[im_r_ilinmom][ijack]=compute_Zq_sig1(jprop_inv[im_r][ijack],mom);
	  Zq_time.stop();
	  
	  //do the same with QED
	  if(pars::use_QED)
	    {
	      Zq_time.start();
	      Zq_QED[im_r_ilinmom][ijack]=-compute_Zq(jprop_QED_inv[im_r][ijack],mom);
	      Zq_sig1_QED[im_r_ilinmom][ijack]=-compute_Zq_sig1(jprop_QED_inv[im_r][ijack],mom);
	      Zq_time.stop();
	    }
	}
    }
}
