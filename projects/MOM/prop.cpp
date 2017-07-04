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

void read_prop(prop_t *prop,raw_file_t &file,const dcompl_t &fact)
{
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

string get_prop_tag(size_t im,size_t ir,const string &ins)
{return combine("S_M%zu_R%zu_%s",im,ir,ins.c_str());}

vjprop_t get_all_mr_props_inv(const vjprop_t &jprop)
{
  cout<<"Inverting all props"<<endl;
  
  vjprop_t jprop_inv(jprop.size());
  
#pragma omp parallel for
  for(size_t imr=0;imr<nmr;imr++)
    for(size_t ijack=0;ijack<=njacks;ijack++)
      jprop_inv[imr][ijack]=jprop[imr][ijack].inverse();
  
  return jprop_inv;
}

void build_all_mr_jackknifed_props(vector<jm_r_mom_props_t> &jprops,const vector<m_r_mom_conf_props_t> &props,bool set_QED,const index_t &im_r_ind)
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
	for(auto &jp_p : vector<pair<jprop_t*,const prop_t*>>({{&j.EM,&p.FF},{&j.EM,&p.T},{&j.P,&p.P},{&j.S,&p.S}}))
	  (*jp_p.first)[ijack]+=*jp_p.second;
    }
}

void clusterize_all_mr_jackknifed_props(vector<jm_r_mom_props_t> &jprops,bool use_QED,size_t clust_size)
{
#pragma omp parallel for
  for(size_t iprop=0;iprop<jprops.size();iprop++)
    jprops[iprop].clusterize_all_mr_props(use_QED,clust_size);
}

void finish_jprops_EM(vector<jm_r_mom_props_t> &jprops,const djack_t &deltam_cr)
{
  index_t ind({{"iprop",jprops.size()},{"ijack",njacks+1}});
#pragma omp parallel for
  for(size_t i=0;i<ind.max();i++)
    {
      const vector<size_t> comp=ind(i);
      const size_t iprop=comp[0],ijack=comp[1];
      jprops[iprop].EM[ijack]-=jprops[iprop].P[ijack]*deltam_cr[ijack];
    }
}

double m0_of_kappa(double kappa)
{return 1.0/(2.0*kappa)-4;}

double M_of_p(const p_t &p,double kappa)
{return m0_of_kappa(kappa)+p.hat().norm2()/2.0;}

prop_t free_prop(const imom_t &pi,double mu,double kappa,size_t r)
{
  double tau3[2]={-1.0,+1.0};
  
  const p_t p=pi.p(L);
  const p_t ptilde=p.tilde();
  
  double M=M_of_p(p,kappa);
  dcompl_t I=dcompl_t(0.0,1.0);
  
  prop_t out=-I*slash(ptilde)+mu*Gamma[0]+I*M*Gamma[5]*tau3[r];
  out/=sqr(M)+sqr(mu)+ptilde.norm2();
  out/=V;
  
  return out;
}
