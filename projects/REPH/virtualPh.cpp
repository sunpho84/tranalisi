#include <tranalisi.hpp>

const int L=32;
int T;
const int ndim=4;

enum HAV{HA=0,HV=1};
const char havTag[2][2]={"A","V"};

vector<djvec_t> threePts;
vector<double> eG;

constexpr int eps[4][4]={{0,0,0,0},
			 {0,-1,-1,0},
			 {0,+1,-1,0},
			 {0,0,0,0}};

index_t index3pts;

struct flavour_t
{
  int qhat;
  double kappa;
  double mu;
  double su3csw;
  double u1csw;
  double cF;
  double cF_prime;
  double th1;
  double th2;
  double th3;
};

struct inv_t
{
  int isolv;
  double mu,th[3];
};

struct einv_t
{
  int isolv;
  double mu,th0[3],tht[3],off;
};

struct combination_t
{
  int i0,it,is;
  double mu1,mu2,off;
  double th0[3],tht[3],ths[3];
};

ostream& operator<<(ostream& os,const combination_t& c)
{
#define PRINT(A)				\
  os<<" "<<#A<<" "<<c.A<<endl
  
  PRINT(i0);
  PRINT(it);
  PRINT(is);
  PRINT(mu1);
  PRINT(mu2);
  PRINT(off);
  PRINT(th0[0]);
  PRINT(th0[1]);
  PRINT(th0[2]);
  PRINT(tht[0]);
  PRINT(tht[1]);
  PRINT(tht[2]);
  PRINT(ths[0]);
  PRINT(ths[1]);
  PRINT(ths[2]);
  os<<endl;
  
#undef PRINT
  
  return os;
}

struct file_head_t
{
  int tmax;
  int x0;
  int stype;
  int phptype;
  int z0;
  int ninv;
  int neinv;
  int nsolv;
  int nhits;
  int ncomb;
  int ngsm;
  double epsgsm;
  int nqsml,nqsm0,nqsm;
  double epsqsm;
  flavour_t gflv;
  vector<combination_t> comb;
  vector<inv_t> inv;
  vector<einv_t> einv;
  
  size_t iOppComb(const size_t icomb) const
  {
    size_t iout=0;
    
    const combination_t& c=comb[icomb];
    
    bool found=false;
    
    do
      {
	const combination_t& d=comb[iout];
	found=(c.mu1==d.mu2 and
	       c.mu2==d.mu1 and
	       c.th0[2]==d.th0[2] and
	       c.tht[2]==d.tht[2] and
	       c.ths[2]==d.ths[2]);
	
	if(not found) iout++;
      }
    while(iout<comb.size() and not found);
    
    return iout;
  }
};

file_head_t f;

struct data_t
{
  int nc;   ///< Conf id
  int size;
  vector<double> HA,HV;
};

struct data2_t
{
  int nc,size;
};

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

void readHeader(raw_file_t& fin)
{
  for(auto& i : {&f.tmax,&f.x0,&f.stype,&f.phptype,&f.z0,&f.ninv,&f.neinv,&f.nsolv,&f.nhits,&f.ncomb,&f.ngsm,&f.nqsml,&f.nqsm0,&f.nqsm})
    fin.bin_read(*i);
  
  cout<<f.ncomb<<endl;
  
  T=f.tmax;
  f.inv.resize(f.ninv);
  f.comb.resize(f.ncomb);
  f.einv.resize(f.neinv);
  
  for(auto& d : {&f.epsgsm,&f.epsqsm})
    fin.bin_read(*d);
   
  fin.bin_read(f.gflv.qhat);
  
  for(auto& d : {&f.gflv.kappa,&f.gflv.mu,&f.gflv.su3csw,&f.gflv.u1csw,&f.gflv.cF,&f.gflv.cF_prime,&f.gflv.th1,&f.gflv.th2,&f.gflv.th3})
    fin.bin_read(*d);
  
  for(int icomb=0;icomb<f.ncomb;++icomb)
    {
      auto& c=f.comb[icomb];
      
      for(auto& i : {&c.i0,&c.it,&c.is})
	fin.bin_read(*i);
      
      for(auto& d : {&c.mu1,&c.mu2,&c.off,&c.th0[0],&c.th0[1],&c.th0[2],&c.tht[0],&c.tht[1],&c.tht[2],&c.ths[0],&c.ths[1],&c.ths[2]})
	fin.bin_read(*d);
      
      cout<<c<<endl;
    }
  
  // fills eG
  eG.resize(f.ncomb);
  for(int icomb=0;icomb<f.ncomb;++icomb)
    {
      eG[icomb]=Eg(icomb);
      cout<<"Eg["<<icomb<<"]: "<<eG[icomb]<<endl;
    }
  
  // prints opposite of all combo
  for(int icomb=0;icomb<f.ncomb;++icomb)
    cout<<"Opposite of comb "<<icomb<<": "<<f.iOppComb(icomb)<<endl;
  
  for(int inv=0;inv<f.ninv;++inv)
    {
      auto& v=f.inv[inv];
      
      for(auto& d : {&v.mu,&v.th[0],&v.th[1],&v.th[2]})
	fin.bin_read(*d);
    }
}

/// three pts
vector<djvec_t> readData()
{
  raw_file_t fin("data/conf.virtualph.dat","r");
  readHeader(fin);
  
  const int ncorrs=2;
  index3pts.set_ranges({{"corr",ncorrs},{"comb",f.ncomb},{"sl",f.nqsml},{"mu",ndim},{"alpha",ndim},{"ri",2}});
  
  vector<djvec_t> out(index3pts.max(),djvec_t(f.tmax));
  
  data_t data;
  data.size=2*ndim*ndim*f.tmax*f.ncomb*f.nqsml;
  
  int nTotConfs=(fin.size()-fin.get_pos())/(2*data.size*sizeof(double)+sizeof(int));
  cout<<"nTotConfs: "<<nTotConfs<<endl;
  
  const int clustSize=nTotConfs/njacks;
  nTotConfs=clustSize*njacks;
  cout<<"nTotConfs after rounding: "<<nTotConfs<<endl;
  
  vector<double> corr(data.size);
  
  for(int ic=0;ic<nTotConfs;ic++)
    {
      fin.bin_read(data.nc);
      //cout<<"nc: "<<data.nc<<endl;
      
      const int iclust=ic/clustSize;
      
      //HA,HV
      for(size_t icorr=0;icorr<ncorrs;icorr++)
	{
	  fin.bin_read(corr);
	  
	  for(size_t icomb=0;icomb<(size_t)f.ncomb;icomb++)
	    for(size_t isl=0;isl<(size_t)f.nqsml;isl++)
	      for(size_t mu=0;mu<ndim;mu++)
		for(size_t alpha=0;alpha<ndim;alpha++)
		  for(size_t t=0;t<(size_t)f.tmax;t++)
		    for(size_t ire=0;ire<2;ire++)
		      {
			const int iin=ire+2*(t+f.tmax*(alpha+ndim*(mu+ndim*(isl+f.nqsml*icomb))));
			const int iout=index3pts({icorr,icomb,isl,mu,alpha,ire});
			
			// if(icorr==0 and icomb==0 and isl==0 and mu==0 and alpha==0 and ire==0)
			//   cout<<"t "<<t<<" "<<corr[iin]<<endl;
			
			out[iout][t][iclust]+=corr[iin];
		      }
	}
    }
  
  for(auto& j : out)
    j.clusterize(clustSize);
  
  return out;
}

int index2pts(int icorr,int iinv2,int iinv1,int isl,int ire)
{
  return ire+2*(isl+f.nqsml*(iinv1+f.ninv*(iinv2+f.ninv*icorr)));
}

vector<djvec_t> readData2()
{
  raw_file_t fin("data/conf.virtualph.dat2","r");
  readHeader(fin);
  
  const int ncorrs=5;
  vector<djvec_t> out(index2pts(ncorrs,0,0,0,0),djvec_t(f.tmax));
  
  data2_t data2;
  data2.size=2*f.tmax*f.ninv*f.ninv*f.nqsml;
  
  int nTotConfs=(fin.size()-fin.get_pos())/(ncorrs*data2.size*sizeof(double)+sizeof(int));
  cout<<"nTotConfs: "<<nTotConfs<<endl;
  
  const int clustSize=nTotConfs/njacks;
  nTotConfs=clustSize*njacks;
  cout<<"nTotConfs after rounding: "<<nTotConfs<<endl;
  
  vector<double> corr(data2.size);
  
  for(int ic=0;ic<nTotConfs;ic++)
    {
      fin.bin_read(data2.nc);
      
      const int iclust=ic/clustSize;
      
      //PP,PA0,PA1,PA2,PA3
      for(int icorr=0;icorr<5;icorr++)
	{
	  fin.bin_read(corr);
	  
	  for(int iinv2=0;iinv2<f.ninv;iinv2++)
	    for(int iinv1=0;iinv1<f.ninv;iinv1++)
	      for(int isl=0;isl<f.nqsml;isl++)
	      for(int t=0;t<f.tmax;t++)
		for(int ire=0;ire<2;ire++)
		  {
		    const int iin=ire+2*(t+f.tmax*(isl+f.nqsml*(iinv1+f.ninv*iinv2)));
		    const int iout=index2pts(icorr,iinv2,iinv1,isl,ire);
		    
		    out[iout][t][iclust]+=corr[iin];
		  }
	}
    }
  
  for(auto& j : out)
    j.clusterize(clustSize);
  
  return out;
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
  
  return out.symmetrized(0*((hAV==0)?+1:-1));
}

/// Load a given polarization r, and weak current alpha
djvec_t load3ptsPol(const HAV hAV,const size_t icomb,const size_t isl,const size_t alpha,const size_t r) // r and alpha go from 1 to 2
{
  auto t=bind(load3pts,hAV,icomb,isl,_1,alpha);
  
  return t(1)*eps[r][1]+t(2)*eps[r][2];
}

djvec_t load3ptsHA(const size_t icomb,const size_t isl)
{
  djvec_t res(T/2+1);
  for(int alpha=1;alpha<=2;alpha++)
    for(int r=1;r<=2;r++)
      res+=load3ptsPol(HA,icomb,isl,alpha,r)/eps[r][alpha];
  
  return res;
}

djvec_t load3ptsHV(const size_t icomb,const size_t isl)
{
  const double c[2][2]={{-1,+1},{-1,-1}};
  djvec_t res(T/2+1);
  for(int alpha=1;alpha<=2;alpha++)
    for(int r=1;r<=2;r++)
      res+=load3ptsPol(HV,icomb,isl,alpha,r)*c[r-1][alpha-1];
  
  return res;
}

int main()
{
  set_njacks(10);
  
  threePts=readData();
  const auto twoPts=readData2();
  
  cout<<"f.nqsml: "<<f.nqsml<<endl;
  
  for(size_t i=0;i<index3pts.max();i++)
    {
      threePts[i].ave_err().write(combine("plots/threeptsRaw/%s.xmg",index3pts.descr(i).c_str()));
      threePts[i].symmetrized(-1).ave_err().write(combine("plots/threeptsRaw/antisymm_%s.xmg",index3pts.descr(i).c_str()));
      threePts[i].symmetrized(+1).ave_err().write(combine("plots/threeptsRaw/symm_%s.xmg",index3pts.descr(i).c_str()));
    }
  
  enum{PP,PA0,PA1,PA2,PA3};
  const int is=0,il=1; //index to be fetched from inv list
  const int tmin=20,tmax=25;
  djvec_t mP(f.nqsml),ZP(f.nqsm);
  for(size_t isl=0;isl<(size_t)f.nqsml;isl++)
    {
      djack_t Z;
      two_pts_fit(Z,mP[isl],twoPts[index2pts(PP,is,il,isl,0)].symmetrized(),T/2,tmin,tmax,combine("plots/PP_ll_sm%zu.xmg",isl));
      if(isl==0) ZP[isl]=sqrt(Z);
      else ZP[isl]=Z/ZP[0];
      cout<<"mP["<<isl<<"]: "<<mP[isl].ave_err()<<endl;
    }
  
  djvec_t norm[2];
  for(size_t isl=0;isl<(size_t)f.nqsml;isl++)
    {
      norm[isl]=djvec_t(T/2+1);
      for(int it=0;it<T/2;it++)
	norm[isl][it]=ZP[isl]*exp(-mP[isl]*it)/(2*mP[isl]);
    }
  
  for(size_t hAV=0;hAV<2;hAV++)
    for(size_t icomb=0;icomb<(size_t)f.ncomb;icomb++)
      for(size_t isl=0;isl<(size_t)f.nqsml;isl++)
	{
	  const djvec_t out=((hAV==0)?load3ptsHA:load3ptsHV)(icomb,isl)/norm[isl];
	  out.ave_err().write(combine("plots/threepts/H%s_comb%zu_isl%zu.xmg",havTag[hAV],icomb,isl));
	  
	  const djvec_t outImpr=out-((hAV==0)?load3ptsHA:load3ptsHV)(0,isl)/norm[isl];
	  outImpr.ave_err().write(combine("plots/threeptsImpr/H%s_comb%zu_isl%zu.xmg",havTag[hAV],icomb,isl));
	  
	  for(size_t alpha=1; alpha<=2;alpha++)
	    for(size_t r=1;r<=2;r++)
	      {
		const djvec_t out=load3ptsPol((HAV)hAV,icomb,isl,alpha,r)/norm[isl];
		out.ave_err().write(combine("plots/threeptsPol/H%s_comb%zu_isl%zu_al%zu_r%zu.xmg",havTag[hAV],icomb,isl,alpha,r));
	      }
	}
  
  // const djvec_t with_without_sme_ratio=
  //   twoPts[index2pts(PP,is,il,1/*with*/,0)].symmetrized()/twoPts[index2pts(PP,is,il,0/*without*/,0)].symmetrized();
  // with_without_sme_ratio.ave_err().write("plots/with_without_PP.xmg");
  
  grace_file_t t0_plot("plots/t0.xmg");
  t0_plot.set_title("at rest, combination of both insertions");
  grace_file_t t0_insS_plot("plots/t0_insS.xmg");
  t0_insS_plot.set_title("at rest, insertion on S");
  grace_file_t t0_insT_plot("plots/t0_insT.xmg");
  t0_insT_plot.set_title("at rest, insertion on T");
  
  grace_file_t t1_plot("plots/t1.xmg");
  t1_plot.set_title("in motion, combination of both insertions");
  grace_file_t t1_insS_plot("plots/t1_insS.xmg");
  t1_insS_plot.set_title("in motion, insertion on S");
  grace_file_t t1_insT_plot("plots/t1_insT.xmg");
  t1_insT_plot.set_title("in motion, insertion on T");
  
  grace_file_t r_plot("plots/r.xmg");
  
  for(size_t isl=0;isl<(size_t)f.nqsml;isl++)
    {
      const double eT=+2.0/3;
      const double eS=-1.0/3;
      
      const double coeffT[2]={-1.0,+1.0};
      
      const size_t AT_REST_S=0;
      const size_t AT_REST_T=f.iOppComb(AT_REST_S);
      
      cout<<"---"<<endl;
      const djvec_t t0_insS=load3ptsHV(AT_REST_S,isl)/norm[isl];
      cout<<"---"<<endl;
      const djvec_t t0_insT=load3ptsHV(AT_REST_T,isl)*coeffT[HV]/norm[isl];
      cout<<"---"<<endl;
      
      t0_insS_plot.write_vec_ave_err(t0_insS.ave_err());
      t0_insT_plot.write_vec_ave_err(t0_insT.ave_err());
      
      const djvec_t t0=
	t0_insS*eS+
	t0_insT*eT;
      
      t0_plot.write_vec_ave_err(t0.ave_err());
      
      /////////////////////////////////////////////////////////////////
      
      const size_t IN_MOTION_S=1;
      const size_t IN_MOTION_T=f.iOppComb(IN_MOTION_S);
      
      const djvec_t t1_insS=load3ptsHV(IN_MOTION_S,isl)/norm[isl];
      const djvec_t t1_insT=load3ptsHV(IN_MOTION_T,isl)*coeffT[HV]/norm[isl];
      
      t1_insS_plot.write_vec_ave_err(t1_insS.ave_err());
      t1_insT_plot.write_vec_ave_err(t1_insT.ave_err());
      
      const djvec_t t1=
	t1_insS*eS+
	t1_insT*eT;
      
      t1_plot.write_vec_ave_err(t1.ave_err());
      
      // // effective_mass(t0,T/2,-1).ave_err().write("plots/t0_sml"+to_string(isl)+".xmg");
      // // effective_mass(t1,T/2,-1).ave_err().write("plots/t1_sml"+to_string(isl)+".xmg");
      
      const djvec_t r=t1/t0;//-1.0;
      r_plot.write_vec_ave_err(r.ave_err());
      
      // // const djack_t xG=2*eG/mP[isl];
      // // cout<<"Xg: "<<smart_print(xG)<<endl;
      // r.ave_err().write(combine("plots/threePts_sml%zu.xmg",isl));
      
      for(auto& p :  {&t0_plot,&t0_insS_plot,&t0_insT_plot,
			&t1_plot,&t1_insS_plot,&t1_insT_plot,
			&r_plot})
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
  
  for(int isl=0;isl<2;isl++)
  for(int icomb=0;icomb<2;icomb++)
    {
      const djvec_t a=load3pts(HV,icomb,isl,1,1)/norm[isl];
      const djvec_t b=load3pts(HV,icomb,isl,1,2)/norm[isl];
      const djvec_t c=load3pts(HV,icomb,isl,2,1)/norm[isl];
      const djvec_t d=load3pts(HV,icomb,isl,2,2)/norm[isl];
      
      a.ave_err().write(combine("/tmp/a_isl%d_icomb%d.xmg",isl,icomb));
      b.ave_err().write(combine("/tmp/b_isl%d_icomb%d.xmg",isl,icomb));
      c.ave_err().write(combine("/tmp/c_isl%d_icomb%d.xmg",isl,icomb));
      d.ave_err().write(combine("/tmp/d_isl%d_icomb%d.xmg",isl,icomb));
    }
  return 0;
}
