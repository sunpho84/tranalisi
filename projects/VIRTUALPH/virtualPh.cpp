#include <tranalisi.hpp>

#include "nazarioStructures.hpp"

const int L=32;
int T;
// const char havTag[2][2]={"A","V"};

vector<djvec_t> threePts;
vector<double> eG;

// constexpr int eps[4][4]={{0,0,0,0},
// 			 {0,-1,-1,0},
// 			 {0,+1,-1,0},
// 			 {0,0,0,0}};

double Eg(const double th0,const double tht)
{
  const double mom0=2*M_PI*th0/L;
  const double momt=2*M_PI*tht/L;
  
  //cout<<"Th0: "<<th0<<", tht: "<<tht<<endl;
  
  const double kDec=mom0-momt;
  const double kHatDec=2*sin(kDec/2);
  
  return 2*asinh(fabs(kHatDec)/2);
}

double Eg(int icomb)
{
  return Eg(f.comb[icomb].th0[2],f.comb[icomb].tht[2]);
}

djvec_t load3pts(const HAV hAV,const size_t icomb,const size_t isl,const size_t mu,const size_t alpha)
{
  const size_t ire=(hAV==HAV::HA)?0:1;
  
  const size_t ic=index3pts({hAV,icomb,isl,mu,alpha,ire});
  cout<<"ic: "<<ic<<" , "<<index3pts.descr(ic)<<endl;
  
  djvec_t out=threePts[ic];
  
  for(int it=0;it<T;it++)
    {
      const double dt=fabs(T/2.0-it);
      const double arg=-eG[icomb]*dt;
      out[it]/=exp(arg);
    }
  
  return out.symmetrized((hAV==0)?+1:-1);
}

// /// Load a given polarization r, and weak current alpha
// djvec_t load3ptsPol(const HAV hAV,const size_t icomb,const size_t isl,const size_t alpha,const size_t r) // r and alpha go from 1 to 2
// {
//   auto t=bind(load3pts,hAV,icomb,isl,_1,alpha);
  
//   return t(1)*eps[r][1]+t(2)*eps[r][2];
// }

djvec_t load3ptsHA(const size_t icomb,const size_t isl)
{
  const djvec_t a11=load3pts(HA,icomb,isl,1,1);
  const djvec_t a22=load3pts(HA,icomb,isl,2,2);
  
  return (a11+a22);
}

djvec_t load3ptsHV(const size_t icomb,const size_t isl)
{
  const djvec_t a12=load3pts(HV,icomb,isl,1,2);
  const djvec_t a21=load3pts(HV,icomb,isl,2,1);
  
  return (a21-a12);
}

int main()
{
  set_njacks(3);
  
  threePts=readData();
  const auto twoPts=readData2();
  
  T=f.tmax;
  // fills eG
  eG.resize(f.ncomb);
  for(int icomb=0;icomb<f.ncomb;++icomb)
    {
      eG[icomb]=Eg(icomb);
      cout<<"Eg["<<icomb<<"]: "<<eG[icomb]<<endl;
    }
  
  cout<<"f.nqsml: "<<f.nqsml<<endl;
  
  enum{PP,PA0,PA1,PA2,PA3};
  const int is=0,il=1; //index to be fetched from inv list
  const int tmin=15,tmax=28;
  djvec_t mP(f.nqsml),ZP(f.nqsml);
  for(size_t isl=0;isl<(size_t)f.nqsml;isl++)
    {
      djack_t Z;
      const djvec_t c=twoPts[index2pts(PP,is,il,isl,0)].symmetrized();
      c.ave_err().write(combine("plots/raw_PP_ll_sm%zu.xmg",isl));
      
      two_pts_fit(Z,mP[isl],c,T/2,tmin,tmax,combine("plots/PP_ll_sm%zu.xmg",isl));
      if(isl==0) ZP[isl]=sqrt(Z);
      else ZP[isl]=Z/ZP[0];
      cout<<"mP["<<isl<<"]: "<<mP[isl].ave_err()<<endl;
      cout<<"Z["<<isl<<"]: "<<ZP[isl].ave_err()<<endl;
    }
  
  djvec_t norm[2];
  for(size_t isl=0;isl<(size_t)f.nqsml;isl++)
    {
      norm[isl]=djvec_t(T/2+1);
      for(int it=0;it<=T/2;it++)
	norm[isl][it]=ZP[isl]*exp(-mP[isl]*it)/(4*mP[isl]*mP[isl]);//the second M is in some
    }
  
  // for(size_t i=0;i<index3pts.max();i++)
  //   {
  //     const size_t isl=index3pts(i)[2];
      
  //     threePts[i].ave_err().write(combine("plots/threeptsRaw/%s.xmg",index3pts.descr(i).c_str()));
  //     const djvec_t tsymm=threePts[i].symmetrized(+1)/norm[isl];
  //     const djvec_t tasymm=threePts[i].symmetrized(-1)/norm[isl];
  //     tasymm.ave_err().write(combine("plots/threeptsRaw/antisymm_%s.xmg",index3pts.descr(i).c_str()));
  //     tsymm.ave_err().write(combine("plots/threeptsRaw/symm_%s.xmg",index3pts.descr(i).c_str()));
  //   }
  
  // for(size_t hAV=0;hAV<2;hAV++)
  //   for(size_t icomb=0;icomb<(size_t)f.ncomb;icomb++)
  //     for(size_t isl=0;isl<(size_t)f.nqsml;isl++)
  // 	{
  // 	  const djvec_t out=((hAV==0)?load3ptsHA:load3ptsHV)(icomb,isl)/norm[isl];
  // 	  out.ave_err().write(combine("plots/threepts/H%s_comb%zu_isl%zu.xmg",havTag[hAV],icomb,isl));
	  
  // 	  const djvec_t outImpr=out-((hAV==0)?load3ptsHA:load3ptsHV)(0,isl)/norm[isl];
  // 	  outImpr.ave_err().write(combine("plots/threeptsImpr/H%s_comb%zu_isl%zu.xmg",havTag[hAV],icomb,isl));
	  
  // 	  for(size_t alpha=1; alpha<=2;alpha++)
  // 	    for(size_t r=1;r<=2;r++)
  // 	      {
  // 		const djvec_t out=load3ptsPol((HAV)hAV,icomb,isl,alpha,r)/norm[isl];
  // 		out.ave_err().write(combine("plots/threeptsPol/H%s_comb%zu_isl%zu_al%zu_r%zu.xmg",havTag[hAV],icomb,isl,alpha,r));
  // 	      }
  // 	}
  
  // const djvec_t with_without_sme_ratio=
  //   twoPts[index2pts(PP,is,il,1/*with*/,0)].symmetrized()/twoPts[index2pts(PP,is,il,0/*without*/,0)].symmetrized();
  // with_without_sme_ratio.ave_err().write("plots/with_without_PP.xmg");
  
  grace_file_t a0_plot("plots/a0.xmg");
  a0_plot.set_title("at rest, combination of both insertions");
  grace_file_t a0_insS_plot("plots/a0_insS.xmg");
  a0_insS_plot.set_title("at rest, insertion on S");
  grace_file_t a0_insT_plot("plots/a0_insT.xmg");
  a0_insT_plot.set_title("at rest, insertion on T");
  
  grace_file_t v0_plot("plots/v0.xmg");
  v0_plot.set_title("at rest, combination of both insertions");
  grace_file_t v0_insS_plot("plots/v0_insS.xmg");
  v0_insS_plot.set_title("at rest, insertion on S");
  grace_file_t v0_insT_plot("plots/v0_insT.xmg");
  v0_insT_plot.set_title("at rest, insertion on T");
  
  grace_file_t a1_plot("plots/a1.xmg");
  a1_plot.set_title("in motion, combination of both insertions");
  grace_file_t a1_insS_plot("plots/a1_insS.xmg");
  a1_insS_plot.set_title("in motion, insertion on S");
  grace_file_t a1_insT_plot("plots/a1_insT.xmg");
  a1_insT_plot.set_title("in motion, insertion on T");
  grace_file_t a10_insS_plot("plots/a10_insS.xmg");
  a10_insS_plot.set_title("in motion, insertion on S");
  grace_file_t a10_insT_plot("plots/a10_insT.xmg");
  a10_insT_plot.set_title("in motion, insertion on T");
  
  grace_file_t v1_plot("plots/v1.xmg");
  v1_plot.set_title("in motion, combination of both insertions");
  grace_file_t v1_insS_plot("plots/v1_insS.xmg");
  v1_insS_plot.set_title("in motion, insertion on S");
  grace_file_t v1_insT_plot("plots/v1_insT.xmg");
  v1_insT_plot.set_title("in motion, insertion on T");
  grace_file_t v10_insS_plot("plots/v10_insS.xmg");
  v10_insS_plot.set_title("in motion, insertion on S");
  grace_file_t v10_insT_plot("plots/v10_insT.xmg");
  v10_insT_plot.set_title("in motion, insertion on T");
  
  grace_file_t v10_plot("plots/v10.xmg");
  v10_plot.set_title("in motion diff");
  grace_file_t a10_plot("plots/a10.xmg");
  a10_plot.set_title("in motion diff");
  
  grace_file_t rA_plot("plots/rA.xmg");
  grace_file_t rV_plot("plots/rV.xmg");
  grace_file_t rVbis_plot("plots/rVbis.xmg");
  
  for(size_t isl=0;isl<(size_t)f.nqsml;isl++)
    {
      const double eT=+2.0/3;
      const double eS=-1.0/3;
      
      const double coeffT[2]={-1.0,+1.0};
      
      const size_t AT_REST_S=0;
      const size_t AT_REST_T=f.iOppComb(AT_REST_S);
      
      cout<<"---"<<endl;
      const djvec_t v0_insS=load3ptsHV(AT_REST_S,isl)/norm[isl];
      const djvec_t a0_insS=load3ptsHA(AT_REST_S,isl)/norm[isl];
      cout<<"---"<<endl;
      const djvec_t v0_insT=load3ptsHV(AT_REST_T,isl)*coeffT[HV]/norm[isl];
      const djvec_t a0_insT=load3ptsHA(AT_REST_T,isl)*coeffT[HA]/norm[isl];
      cout<<"---"<<endl;
      
      a0_insS_plot.write_vec_ave_err(a0_insS.ave_err());
      v0_insS_plot.write_vec_ave_err(v0_insS.ave_err());
      a0_insT_plot.write_vec_ave_err(a0_insT.ave_err());
      v0_insT_plot.write_vec_ave_err(v0_insT.ave_err());
      
      const djvec_t a0=
	a0_insS*eS+
	a0_insT*eT;
      
      a0_plot.write_vec_ave_err(a0.ave_err());
      
      const djvec_t v0=
	v0_insS*eS+
	v0_insT*eT;
      
      v0_plot.write_vec_ave_err(v0.ave_err());
      
      /////////////////////////////////////////////////////////////////
      
      const size_t IN_MOTION_S=1;
      const size_t IN_MOTION_T=f.iOppComb(IN_MOTION_S);
      
      const djvec_t v1_insS=load3ptsHV(IN_MOTION_S,isl)/norm[isl];
      const djvec_t a1_insS=load3ptsHA(IN_MOTION_S,isl)/norm[isl];
      const djvec_t v1_insT=load3ptsHV(IN_MOTION_T,isl)*coeffT[HV]/norm[isl];
      const djvec_t a1_insT=load3ptsHA(IN_MOTION_T,isl)*coeffT[HA]/norm[isl];
      
      v1_insS_plot.write_vec_ave_err(v1_insS.ave_err());
      a1_insS_plot.write_vec_ave_err(a1_insS.ave_err());
      v1_insT_plot.write_vec_ave_err(v1_insT.ave_err());
      a1_insT_plot.write_vec_ave_err(a1_insT.ave_err());
      
      const djvec_t a1=
	a1_insS*eS+
	a1_insT*eT;
      
      a1_plot.write_vec_ave_err(a1.ave_err());
      
      const djvec_t v1=
	v1_insS*eS+
	v1_insT*eT;
      
      v1_plot.write_vec_ave_err(v1.ave_err());
      
      /////////////////////////////////////////////////////////////////
      
      const djvec_t v10_insS=v1_insS-v0_insS;
      const djvec_t a10_insS=a1_insS-a0_insS;
      const djvec_t v10_insT=v1_insT-v0_insT;
      const djvec_t a10_insT=a1_insT-a0_insT;
      
      const djvec_t v10=v1-v0;
      const djvec_t a10=a1-a0;
      
      v10_insS_plot.write_vec_ave_err(v10_insS.ave_err());
      a10_insS_plot.write_vec_ave_err(a10_insS.ave_err());
      v10_insT_plot.write_vec_ave_err(v10_insT.ave_err());
      a10_insT_plot.write_vec_ave_err(a10_insT.ave_err());
      
      v10_plot.write_vec_ave_err(v10.ave_err());
      a10_plot.write_vec_ave_err(a10.ave_err());
      
      // // effective_mass(a0,T/2,-1).ave_err().write("plots/a0_sml"+to_string(isl)+".xmg");
      // // effective_mass(a1,T/2,-1).ave_err().write("plots/a1_sml"+to_string(isl)+".xmg");
      
      const djvec_t rV=v10/a0;
      rV_plot.write_vec_ave_err(rV.ave_err());
      
      const djvec_t rVbis=v1/a0;
      rVbis_plot.write_vec_ave_err(rVbis.ave_err());
      
      const djvec_t rA=a10/a0;
      rA_plot.write_vec_ave_err(rA.ave_err());
      
      // // const djack_t xG=2*eG/mP[isl];
      // // cout<<"Xg: "<<smart_print(xG)<<endl;
      // r.ave_err().write(combine("plots/threePts_sml%zu.xmg",isl));
      
      for(auto& p : {  &a0_plot,&a0_insS_plot,&a0_insT_plot,
      		       &a1_plot,&a1_insS_plot,&a1_insT_plot,
      		       &a10_plot,&a10_insS_plot,&a10_insT_plot,
      		       &v0_plot,&v0_insS_plot,&v0_insT_plot,
      		       &v1_plot,&v1_insS_plot,&v1_insT_plot,
      		       &v10_plot,&v10_insS_plot,&v10_insT_plot,
      		       &rA_plot})
      	{
      	  p->set_legend((isl==0)?"NOT SMEARED":"SMEARED");
      	  p->set_all_colors((isl==0)?grace::RED:grace::BLUE);
      	}
    }
  
  // for(size_t icombo=0;icombo<(size_t)f.ncomb;icombo++)
  //   for(size_t isl=0;isl<1;isl++)
  //     for(size_t mu=0;mu<ndim;mu++)
  // 	for(size_t alpha=0;alpha<ndim;alpha++)
  // 	for(size_t ri=0;ri<2;ri++)
  // 	  {
  // 	    const size_t i=index3pts({HA,icombo,isl,mu,alpha,ri});
  // 	    auto b=threePts[i];
  // 	    b.ave_err().write(combine("plots/naz/HA_%s.xmg",index3pts.escaped_descr(i).c_str()));
  // 	  }
  // dt.ave_err().write("plots/naz/dt.xmg");
  
  // for(double th=1.9;th<2.1;th+=0.01)
  //   cout<<th<<" "<<2*Eg(0.0,th)/mP[0].ave()<<endl;
  
  // for(int isl=0;isl<2;isl++)
  // for(int icomb=0;icomb<4;icomb++)
  //   {
  //     const djvec_t a=load3pts(HV,icomb,isl,1,1)/norm[isl];
  //     const djvec_t b=load3pts(HV,icomb,isl,1,2)/norm[isl];
  //     const djvec_t b0=load3pts(HV,icomb-icomb%2,isl,1,2)/norm[isl];
  //     const djvec_t c=load3pts(HV,icomb,isl,2,1)/norm[isl];
  //     const djvec_t c0=load3pts(HV,icomb-icomb%2,isl,2,1)/norm[isl];
  //     const djvec_t bc=(b-c)/2;
  //     const djvec_t bc0=(b0-c0)/2;
  //     const djvec_t bc_sub=bc-bc0;
  //     const djvec_t d=load3pts(HV,icomb,isl,2,2)/norm[isl];
      
  //     a.ave_err().write(combine("/tmp/Va_isl%d_icomb%d.xmg",isl,icomb));
  //     b.ave_err().write(combine("/tmp/Vb_isl%d_icomb%d.xmg",isl,icomb));
  //     bc.ave_err().write(combine("/tmp/Vbc_isl%d_icomb%d.xmg",isl,icomb));
  //     bc_sub.ave_err().write(combine("/tmp/Vbc_sub_isl%d_icomb%d.xmg",isl,icomb));
  //     c.ave_err().write(combine("/tmp/Vc_isl%d_icomb%d.xmg",isl,icomb));
  //     d.ave_err().write(combine("/tmp/Vd_isl%d_icomb%d.xmg",isl,icomb));
  //   }
  
  // for(int isl=0;isl<2;isl++)
  // for(int icomb=0;icomb<4;icomb++)
  //   {
  //     const djvec_t a=load3pts(HA,icomb,isl,1,1)/norm[isl];
  //     const djvec_t a0=load3pts(HA,icomb-icomb%2,isl,1,1)/norm[isl];
  //     const djvec_t b=load3pts(HA,icomb,isl,1,2)/norm[isl];
  //     const djvec_t c=load3pts(HA,icomb,isl,2,1)/norm[isl];
  //     const djvec_t d=load3pts(HA,icomb,isl,2,2)/norm[isl];
  //     const djvec_t d0=load3pts(HA,icomb-icomb%2,isl,2,2)/norm[isl];
  //     const djvec_t ad=(a+d)/2;
  //     const djvec_t ad0=(a0+d0)/2;
  //     const djvec_t ad_sub=ad-ad0;
      
  //     a.ave_err().write(combine("/tmp/Aa_isl%d_icomb%d.xmg",isl,icomb));
  //     b.ave_err().write(combine("/tmp/Ab_isl%d_icomb%d.xmg",isl,icomb));
  //     ad.ave_err().write(combine("/tmp/Aad_isl%d_icomb%d.xmg",isl,icomb));
  //     ad_sub.ave_err().write(combine("/tmp/Aad_sub_isl%d_icomb%d.xmg",isl,icomb));
  //     c.ave_err().write(combine("/tmp/Ac_isl%d_icomb%d.xmg",isl,icomb));
  //     d.ave_err().write(combine("/tmp/Ad_isl%d_icomb%d.xmg",isl,icomb));
  //   }
  
  return 0;
}
