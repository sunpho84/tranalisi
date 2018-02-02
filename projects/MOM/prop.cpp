#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_PROP
 #include <prop.hpp>

#include <oper.hpp>

#include <corrections.hpp>
#include <geometry.hpp>
#include <types.hpp>

const char m_r_mom_conf_props_t::tag[NPROP_WITH_QED][3]={"0","FF","F","T","S","P"};

void read_prop(qprop_t *prop,raw_file_t &file,const dcompl_t &fact,const size_t imom)
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

string get_prop_tag(const size_t im,const size_t ir,const size_t ikind)
{
  return combine("S_M%zu_R%zu_%s",im,ir,m_r_mom_conf_props_t::tag[ikind]);
}

vjqprop_t get_all_mr_props_inv(const vjqprop_t &jprop)
{
  cout<<"Inverting all props"<<endl;
  
  vjqprop_t jprop_inv(jprop.size());
  
#pragma omp parallel for
  for(size_t imr=0;imr<glb::nmr;imr++)
    for(size_t ijack=0;ijack<=njacks;ijack++)
      jprop_inv[imr][ijack]=jprop[imr][ijack].inverse();
  
  return jprop_inv;
}

void build_all_mr_jackknifed_props(vector<jm_r_mom_props_t> &jprops,const vector<m_r_mom_conf_props_t> &props,bool set_QED,const index_t &im_r_ind,const djvec_t &deltam_cr)
{
  const index_t im_r_ijack_ind=im_r_ind*index_t({{"ijack",njacks}});
  
#pragma omp parallel for
  for(size_t i_im_r_ijack=0;i_im_r_ijack<im_r_ijack_ind.max();i_im_r_ijack++)
    {
      vector<size_t> im_r_ijack=im_r_ijack_ind(i_im_r_ijack);
      size_t im=im_r_ijack[0],r=im_r_ijack[1],ijack=im_r_ijack[2];
      size_t im_r=im_r_ind({im,r});
      //cout<<"Building jack prop im="<<im<<", r="<<r<<", ijack="<<ijack<<endl;
      
      jm_r_mom_props_t &j=jprops[im_r];
      const m_r_mom_conf_props_t &p=props[i_im_r_ijack];
      
      j.LO[ijack]+=p.LO;
      if(set_QED)
	for(auto &jp_p : vector<tuple<jqprop_t*,const qprop_t*,double>>({{&j.EM,&p.FF,1.0},{&j.EM,&p.T,1.0},{&j.EM,&p.P,-deltam_cr[im][ijack]},{&j.S,&p.S,1.0}}))
	  (*get<0>(jp_p))[ijack]+=*get<1>(jp_p)*get<2>(jp_p);
    }
}

void clusterize_all_mr_jackknifed_props(vector<jm_r_mom_props_t> &jprops,bool use_QED,size_t clust_size)
{
#pragma omp parallel for
  for(size_t iprop=0;iprop<jprops.size();iprop++)
    jprops[iprop].clusterize_all_mr_props(use_QED,clust_size);
}

double m0_of_kappa(double kappa)
{
  return 1.0/(2.0*kappa)-4;
}

double M_of_p(const p_t &p,double kappa)
{
  return m0_of_kappa(kappa)+p.hat().norm2()/2.0;
}

qprop_t free_prop(const imom_t &pi,double mu,double kappa,size_t r)
{
  const p_t p=pi.p(L);
  const p_t ptilde=p.tilde();
  
  double M=M_of_p(p,kappa);
  dcompl_t I=dcompl_t(0.0,1.0);
  
  qprop_t out=-I*slash(ptilde)+mu*Gamma[0]+I*M*Gamma[5]*tau3[r];
  out/=sqr(M)+sqr(mu)+ptilde.norm2();
  out/=V;
  
  return out;
}
