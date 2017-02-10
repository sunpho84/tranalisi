#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <common.hpp>
#include <gm2_IB_integrators.hpp>

const size_t ilight=0,istrange=1,icharm=2;

size_t nm,nr;
index_t<4> ind;

double Za_ave[nbeta]={0.731,0.737,0.762}; //M1
double Za_err[nbeta]={0.008,0.005,0.004};
// double Za_ave[nbeta]={0.703,0.714,0.752}; //M2
// double Za_err[nbeta]={0.002,0.002,0.002};
int Za_seed[nbeta]={13124,862464,76753};
dbvec_t Za(nbeta);
boot_init_t bi;

const int im=2;
const double eq[3]={eu,es,ec};

//! hold the data for a single ensemble
class ens_data_t
{
public:
  size_t iult; //< input in the ultimate file
  size_t ib,T,L;
  double aml;
  string path;
  
  size_t tmin[3],tmax[3];
  djack_t deltam_cr;
};
vector<ens_data_t> ens_data;

//! read a single vector, for a specific mass and r, real or imaginary
dbvec_t read(const char *what,const ens_data_t &ens,size_t im,size_t r,size_t reim)
{return dbvec_t(bi,read_djvec(combine("%s/data/corr%s",ens.path.c_str(),what),ens.T,ind({im,im,r,reim})));}

//! read a combination of r and return appropriately simmetrized
dbvec_t read(const char *what,const ens_data_t &ens,int tpar,size_t im,int rpar,size_t reim)
{
  dbvec_t o(ens.T);
  o=0.0;
  
  for(size_t r=0;r<nr;r++) o+=read(what,ens,im,r,reim)*(r==0?1:rpar);
  return o.symmetrized(tpar)/(1+abs(rpar));
}

//! read VV averaging the three channels
dbvec_t read_VV_ren(const char *what,const ens_data_t &ens,int tpar,size_t im,int rpar,size_t reim)
{
  return
    dbvec_t(read(combine("%s_V1V1",what).c_str(),ens,tpar,im,rpar,reim)+
	    read(combine("%s_V2V2",what).c_str(),ens,tpar,im,rpar,reim)+
	    read(combine("%s_V3V3",what).c_str(),ens,tpar,im,rpar,reim))*
    dboot_t(sqr(Za[ens.ib])/3.0);
}

//! integrate
dboot_t integrate(const dbvec_t &corr,size_t T,const dboot_t &a)
{
  dboot_t out;
  out=0.0;
  for(size_t t=1;t<T/2-4;t++) out+=corr[t]*ftilde_t(t,a);
  out*=4*sqr(alpha_em)*sqr(eq[im]);
  
  return out;
}

//! compute the critical deltam
dboot_t compute_deltam_cr(const ens_data_t &ens,size_t iq)
{
  dbvec_t V0P5_LL=read("LL_V0P5",ens,-1,iq,-1,IM);
  dbvec_t V0P5_0M=read("0M_V0P5",ens,-1,iq,-1,IM);
  dbvec_t V0P5_0T=read("0T_V0P5",ens,-1,iq,-1,IM);
  dbvec_t num_deltam_cr=forward_derivative(dbvec_t(V0P5_LL+2.0*dbvec_t(V0P5_0M+V0P5_0T)));
  num_deltam_cr.ave_err().write(combine("%s/plots/num_deltam_cr.xmg",ens.path.c_str()));

  dbvec_t V0P5_0P=read("0P_V0P5",ens,-1,iq,+1,RE);
  dbvec_t den_deltam_cr=forward_derivative(V0P5_0P);
  den_deltam_cr.ave_err().write(combine("%s/plots/den_deltam_cr.xmg",ens.path.c_str()));
  
  dboot_t deltam_cr=constant_fit(dbvec_t(-num_deltam_cr/(2.0*den_deltam_cr)),ens.tmin[iq],ens.tmax[iq],combine("%s/plots/deltam_cr_t.xmg",ens.path.c_str()));
  
  return deltam_cr;
}

//! load VV for the leading order
dbvec_t load_LO(const ens_data_t &ens,const int im)
{return read_VV_ren("00",ens,1,im,1,RE);}

//! load QED corrections
dbvec_t load_QED(const ens_data_t &ens,const int im,const dboot_t &deltam_cr,const dbvec_t &VV_LO)
{
  dbvec_t VV_0T=read_VV_ren("0T",ens,1,im,1,RE);
  dbvec_t VV_0M=read_VV_ren("0M",ens,1,im,1,RE);
  dbvec_t VV_LL=read_VV_ren("LL",ens,1,im,1,RE);
  dbvec_t(VV_0T/VV_LO).ave_err().write(combine("%s/plots/VV_0T.xmg",ens.path.c_str()));
  dbvec_t(VV_0M/VV_LO).ave_err().write(combine("%s/plots/VV_0M.xmg",ens.path.c_str()));
  dbvec_t(VV_LL/VV_LO).ave_err().write(combine("%s/plots/VV_LL.xmg",ens.path.c_str()));
  dbvec_t c=dbvec_t(VV_LL+2.0*dbvec_t(VV_0T+VV_0M));
  
  dbvec_t VV_0P=read_VV_ren("0P",ens,1,im,-1,IM);
  dbvec_t(VV_0P/VV_LO).ave_err().write(combine("%s/plots/VV_0P.xmg",ens.path.c_str()));
  dbvec_t d=-(deltam_cr*VV_0P);
  return dbvec_t(c+2.0*d)*e2*sqr(eq[im]);
}

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz(const Tpars &C,const Tpars &Kpi,const Tm &ml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L2dep,const size_t an_flag)
{return C+Kpi*ml+a*a*(adep+ml*adep_ml)+L2dep/sqr(L);}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit(const dbvec_t &a,const dbvec_t &z,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  
  ave_err_t adep_guess({6e-4,1e-4});
  ave_err_t adep_ml_guess({0.0,0.0});
  pars.add_az_pars(a,z,boot_fit);
  
  pars.add_adep_pars(adep_guess,adep_ml_guess,boot_fit);
  pars.iC=boot_fit.add_fit_par(pars.C,"C",1e-9,1e-9);
  pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",0,1e-9);
  pars.iL2dep=boot_fit.add_fit_par(pars.L2dep,"L2dep",0.0,0.07);
  
  //boot_fit.fix_par(pars.iadep);
  boot_fit.fix_par(pars.iL2dep);
  //boot_fit.fix_par(pars.iKPi);
  boot_fit.fix_par(pars.iadep_ml);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_ansatz(p[pars.iC],p[pars.iKPi],ml,ac,p[pars.iadep],p[pars.iadep_ml],L,p[pars.iL2dep],an_flag);},cov_flag);
  
  double a_cont=0;
  dboot_t phys_res=cont_chir_ansatz(pars.C,pars.KPi,ml_phys,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L2dep,an_flag);
  cout<<"result: "<<phys_res.ave_err()<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz<double,double,double>
		    (pars.C.ave(),pars.KPi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),inf_vol,pars.L2dep.ave(),an_flag);},
		bind(cont_chir_ansatz<dboot_t,double,double>,pars.C,pars.KPi,_1,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L2dep,an_flag),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse-without_with_fse*pars.L2dep/sqr(ext_data[idata].L));},
		ml_phys,phys_res,"$$a_\\mu");
  
  return phys_res;
}


int main(int narg,char **arg)
{
  //open input file
  string name="input_global.txt";
  if(narg>=2) name=arg[1];
  raw_file_t input(name,"r");
  
  cout.precision(16);
  
  //read where to read input and how many ensemble
  string ens_pars=input.read<string>("UltimatePath");
  nm=input.read<size_t>("NMass");
  nr=input.read<size_t>("NR");
  ind.set_ranges({nm,nm,nr,2});
  init_common_IB(ens_pars);
  size_t nens_used=input.read<int>("NEnsemble");
  
  //fill Za
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    Za[ibeta].fill_gauss(Za_ave[ibeta],Za_err[ibeta],Za_seed[ibeta]);
  
  input.expect({"Ens","beta","L","T","aml","tint_cr","tint_ss","tint_cc","path"});
  ens_data.resize(nens_used);
  int input_an_id=0;
  dbvec_t LO(ens_data.size());
  dbvec_t QED(ens_data.size());
  dbvec_t ratio(ens_data.size());
  dbvec_t alist(nbeta),zlist(nbeta);
  for(size_t iens=0;iens<nens_used;iens++)
    {
      cout<<"----------------------- "<<iens<<" ---------------------"<<endl;
      ens_data_t &ens=ens_data[iens];
      input.read(ens.iult);
      input.read(ens.ib);
      input.read(ens.L);
      input.read(ens.T);
      input.read(ens.aml);
      for(size_t iq=0;iq<3;iq++)
	{
	  input.read(ens.tmin[iq]);
	  input.read(ens.tmax[iq]);
	}
      input.read(ens.path);
      
      size_t ib=ens.ib;
      bi=jack_index[input_an_id][ens.iult];
      
      dboot_t a=1/lat_par[input_an_id].ainv[ib];
      
      dboot_t deltam_cr=compute_deltam_cr(ens,ilight);
      
      dbvec_t VV_LO=load_LO(ens,im);
      effective_mass(VV_LO).ave_err().write(combine("%s/plots/VV_LO.xmg",ens.path.c_str()));
      dbvec_t VV_QED=load_QED(ens,im,deltam_cr,VV_LO);
      
      // dbvec_t c=2.0*dbvec_t(VV_0T+VV_0M)+VV_LL;
      // c=c.subset(0,T/2-1);
      // dbvec_t d=deltam_cr*VV_0P.subset(0,T/2-1);
      // c.ave_err().write("plots/c.xmg");
      // d.ave_err().write("plots/d.xmg");
      
      LO[iens]=integrate(VV_LO,ens.T,a);
      cout<<"amu: "<<LO[iens].ave_err()<<endl;
      
      QED[iens]=integrate(VV_QED,ens.T,a);
      cout<<"amu_QED: "<<QED[iens].ave_err()<<endl;
      
      ratio[iens]=QED[iens]/LO[iens];
      cout<<" Ratio: "<<dboot_t(ratio[iens]).ave_err()<<endl;
    }
  
  //prepare the list of a and z
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    {
      alist[ibeta]=1.0/lat_par[input_an_id].ainv[ibeta];
      zlist[ibeta]=lat_par[input_an_id].Z[ibeta];
    }
  
  const size_t an_flag=0;
  vector<cont_chir_fit_data_t> data_LO;
  for(size_t iens=0;iens<ens_data.size();iens++)
    data_LO.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,ens_data[iens].ib,ens_data[iens].L,ratio[iens],ratio[iens]));
  cont_chir_fit(alist,zlist,data_LO,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dM2Pi_flag%zu_an%zu.xmg",an_flag,input_an_id,"%s"),an_flag,false);
  close_integrators();
  
  return 0;
}
