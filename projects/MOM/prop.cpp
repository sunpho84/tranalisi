#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_PROP
 #include <prop.hpp>

#include <oper.hpp>
#include <timings.hpp>

#include <corrections.hpp>
#include <geometry.hpp>
#include <types.hpp>

const char m_r_mom_conf_qprops_t::tag[NPROP_WITH_QED][3]={"0","FF","F","T","S","P"};
const char mom_conf_lprops_t::tag[NPROP_WITH_QED][3]={"0","F"};

void read_qprop(qprop_t *prop,raw_file_t &file,const dcompl_t &fact,const size_t imom)
{
  //cout<<"Seeking file "<<file.get_path()<<" to mom "<<imom<<" position "<<imom*sizeof(dcompl_t)*NSPIN*NSPIN*NCOL*NCOL<<" from "<<file.get_pos()<<endl;
  file.set_pos(imom*sizeof(dcompl_t)*NSPIN*NSPIN*NCOL*NCOL);
  
  for(size_t is_so=0;is_so<NSPIN;is_so++)
    for(size_t ic_so=0;ic_so<NCOL;ic_so++)
      for(size_t is_si=0;is_si<NSPIN;is_si++)
	for(size_t ic_si=0;ic_si<NCOL;ic_si++)
	  {
	    dcompl_t c;
	    file.bin_read(c);
	    (*prop)(isc(is_si,ic_si),isc(is_so,ic_so))=c*fact;
	  }
}

void read_lprop(lprop_t *prop,raw_file_t &file,const dcompl_t &fact,const size_t imom)
{
  //cout<<"Seeking file "<<file.get_path()<<" to mom "<<imom<<" position "<<imom*sizeof(dcompl_t)*NSPIN*NSPIN*NCOL*NCOL<<" from "<<file.get_pos()<<endl;
  file.set_pos(imom*sizeof(dcompl_t)*NSPIN*NSPIN*NCOL*NCOL);
  
  for(size_t is_so=0;is_so<NSPIN;is_so++)
    for(size_t ic_so=0;ic_so<NCOL;ic_so++)
      for(size_t is_si=0;is_si<NSPIN;is_si++)
	for(size_t ic_si=0;ic_si<NCOL;ic_si++)
	  {
	    dcompl_t c;
	    file.bin_read(c);
	    (*prop)(is_si,is_so)=c*fact;
	  }
}

string get_qprop_tag(const size_t im,const size_t ir,const size_t ikind)
{
  return combine("S_M%zu_R%zu_%s",im,ir,m_r_mom_conf_qprops_t::tag[ikind]);
}

string get_lprop_tag(const size_t ikind)
{
  return combine("L_%s",mom_conf_lprops_t::tag[ikind]);
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

void build_all_mr_jackkniffed_qprops(vector<jm_r_mom_qprops_t> &jprops,const vector<m_r_mom_conf_qprops_t> &props,bool use_QED,const index_t &im_r_ind)
{
  const index_t im_r_ijack_ind=im_r_ind*index_t({{"ijack",njacks}});
  
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
      if(use_QED)
	for(auto &jp_p : vector<tuple<jqprop_t*,const qprop_t*>>({
	      {&j.PH,&p.FF},
	      // {&j.PH,&p.T},
	      {&j.CT,&p.P},
	      {&j.S,&p.S}}))
	  (*get<0>(jp_p))[ijack]+=*get<1>(jp_p);
    }
}

void clusterize_all_mr_jackkniffed_qprops(vector<jm_r_mom_qprops_t> &jprops,bool use_QED,size_t clust_size,const index_t &im_r_ind,const djvec_t &deltam_cr)
{
#pragma omp parallel for
  for(size_t iprop=0;iprop<jprops.size();iprop++)
    jprops[iprop].clusterize_all_mr_props(use_QED,clust_size);
  
  if(use_QED)
#pragma omp parallel for
    for(size_t i=0;i<im_r_ind.max();i++)
      for(size_t ijack=0;ijack<=njacks;ijack++)
	jprops[i].QED[ijack]=
	  jprops[i].PH[ijack]-
 #warning excluding counterterm on prop
 	  0.0*
	  deltam_cr[im_r_ind(i)[0]][ijack]*jprops[i].CT[ijack];
}

void get_inverse_propagators(vector<jqprop_t> &jprop_inv,vector<jqprop_t> &jprop_QED_inv,
			     const vector<jm_r_mom_qprops_t> &jprops,
			     const index_t &im_r_ijackp1_ind)
{
  jprop_inv.resize(glb::im_r_ind.max());
  jprop_QED_inv.resize(glb::im_r_ind.max());
  
#pragma omp parallel for reduction(+:invert_time)
  for(size_t im_r_ijack=0;im_r_ijack<im_r_ijackp1_ind.max();im_r_ijack++)
    {
      //decript indices
      const vector<size_t> im_r_ijack_comps=im_r_ijackp1_ind(im_r_ijack);
      const size_t im=im_r_ijack_comps[0],r=im_r_ijack_comps[1],ijack=im_r_ijack_comps[2];
      const size_t im_r=glb::im_r_ind({im,r});
      
      //compute inverse
      invert_time.start();
      qprop_t prop_inv=jprop_inv[im_r][ijack]=jprops[im_r].LO[ijack].inverse();
      invert_time.stop();
      
      //do the same with QED
      if(use_QED)
	{
	  invert_time.start(); //This misses a sign -1 coming from the original inverse
	  jprop_QED_inv[im_r][ijack]=prop_inv*jprops[im_r].QED[ijack]*prop_inv;
	  invert_time.stop();
	}
    }
}

double m0_of_kappa(double kappa)
{
  return 1.0/(2.0*kappa)-4;
}

double M_of_p(const p_t &p,double kappa)
{
  return m0_of_kappa(kappa)+p.hat().norm2()/2.0;
}

lprop_t free_prop(const imom_t &pi,double mu,double kappa,size_t r)
{
  const p_t p=pi.p(L);
  const p_t ptilde=p.tilde();
  
  double M=M_of_p(p,kappa);
  dcompl_t I=dcompl_t(0.0,1.0);
  
  lprop_t out=-I*lep_slash(ptilde)+mu*lepGamma[0]+I*M*lepGamma[5]*tau3[r];
  out/=sqr(M)+sqr(mu)+ptilde.norm2();
  out/=V;
  
  return out;
}
