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
// double Za_ave[nbeta]={0.703,0.714,0.774}; //M2
// double Za_err[nbeta]={0.002,0.002,0.004};
double Zt_ave[nbeta]={0.711,0.724,0.762}; //M1
double Zt_err[nbeta]={0.005,0.004,0.004};
// double Zt_ave[nbeta]={0.700,0.711,0.767}; //M2
// double Zt_err[nbeta]={0.003,0.002,0.002};
int Za_seed[nbeta]={13124,862464,76753};
int Zt_seed[nbeta]={5634,917453,324338};
dbvec_t Za(nbeta),Zt(nbeta);
boot_init_t bi;

const int im=0;
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

//! read averaging the three channels
dbvec_t read(const char *what,const char *pat,const ens_data_t &ens,int tpar,size_t im,int rpar,size_t reim)
{
  return
    dbvec_t(read(combine("%s_%c1%c1",what,pat[0],pat[1]).c_str(),ens,tpar,im,rpar,reim)+
	    read(combine("%s_%c2%c2",what,pat[0],pat[1]).c_str(),ens,tpar,im,rpar,reim)+
	    read(combine("%s_%c3%c3",what,pat[0],pat[1]).c_str(),ens,tpar,im,rpar,reim))
    /3.0;
}

//! read PP
dbvec_t read_PP(const char *what,const ens_data_t &ens,size_t im,int rpar,size_t reim)
{return read(combine("%s_%s",what,"P5P5").c_str(),ens,1,im,rpar,reim);}

//! read VV
dbvec_t read_VV(const char *what,const ens_data_t &ens,size_t im,int rpar,size_t reim)
{return sqr(Za[ens.ib])*read(what,"VV",ens,1,im,rpar,reim);}

//! read TV
dbvec_t read_TV(const char *what,const ens_data_t &ens,size_t im,int rpar,size_t reim)
{return dboot_t(-Za[ens.ib]*Zt[ens.ib])*read(what,"TV",ens,-1,im,rpar,reim);}

//! read VT
dbvec_t read_VT(const char *what,const ens_data_t &ens,size_t im,int rpar,size_t reim)
{return dboot_t(Zt[ens.ib]*Za[ens.ib])*read(what,"VT",ens,-1,im,rpar,reim);}

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

//! read QED corrections
dbvec_t read_QED(const char *pat,const ens_data_t &ens,const int tpar,const int im,const dboot_t &deltam_cr,const dbvec_t &c_LO)
{
  dbvec_t c_0T=read("0T",pat,ens,tpar,im,1,RE);
  dbvec_t c_0M=read("0M",pat,ens,tpar,im,1,RE);
  dbvec_t c_LL=read("LL",pat,ens,tpar,im,1,RE);
  dbvec_t(c_0T/c_LO).ave_err().write(combine("%s/plots/%s_0T.xmg",ens.path.c_str(),pat));
  dbvec_t(c_0M/c_LO).ave_err().write(combine("%s/plots/%s_0M.xmg",ens.path.c_str(),pat));
  dbvec_t(c_LL/c_LO).ave_err().write(combine("%s/plots/%s_LL.xmg",ens.path.c_str(),pat));
  dbvec_t c=dbvec_t(c_LL+2.0*dbvec_t(c_0T+c_0M));
  
  dbvec_t c_0P=read("0P",pat,ens,tpar,im,-1,IM);
  dbvec_t(c_0P/c_LO).ave_err().write(combine("%s/plots/%s_0P.xmg",ens.path.c_str(),pat));
  dbvec_t d=-(deltam_cr*c_0P);
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
  //ave_err_t adep_guess({0,0.1});
  ave_err_t adep_ml_guess({0.0,0.0});
  //ave_err_t adep_ml_guess({0.0,0.1});
  pars.add_az_pars(a,z,boot_fit);
  
  pars.add_adep_pars(adep_guess,adep_ml_guess,boot_fit);
  pars.iC=boot_fit.add_fit_par(pars.C,"C",1e-9,1e-9);
  //pars.iC=boot_fit.add_fit_par(pars.C,"C",1,1e-3);
  pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",0,1e-9);
  //pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",0,0.1);
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
  
  //fill Za,Zt
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    {
      Za[ibeta].fill_gauss(Za_ave[ibeta],Za_err[ibeta],Za_seed[ibeta]);
      Zt[ibeta].fill_gauss(Zt_ave[ibeta],Zt_err[ibeta],Zt_seed[ibeta]);
    }
  
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
      size_t TH=ens.T/2;
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
      
      //load LO
      dbvec_t PP_LO=read_PP("00",ens,im,1,RE);
      dboot_t Z_P,M_P;
      two_pts_fit(Z_P,M_P,PP_LO, TH, ens.tmin[im], ens.tmax[im],combine("%s/plots/PP_LO_emass.xmg",ens.path.c_str()));
      dbvec_t VV_LO=read_VV("00",ens,im,1,RE),eff_VV_LO=effective_mass(VV_LO,TH);
      eff_VV_LO.ave_err().write(combine("%s/plots/VV_LO_emass.xmg",ens.path.c_str()));
      dbvec_t TV_LO=read_TV("00",ens,im,1,RE),eff_TV_LO=effective_mass(TV_LO,TH,-1);
      eff_TV_LO.ave_err().write(combine("%s/plots/TV_LO_emass.xmg",ens.path.c_str()));
      //load QED
      dbvec_t VV_QED=read_QED("VV",ens,1,im,deltam_cr,VV_LO),VV_rat=VV_QED/VV_LO;
      dbvec_t TV_QED=-read_QED("TV",ens,-1,im,deltam_cr,TV_LO),TV_rat=TV_QED/TV_LO;
      dbvec_t(VV_rat).ave_err().write(combine("%s/plots/VV_QED_ratio.xmg",ens.path.c_str()));
      dbvec_t(TV_rat).ave_err().write(combine("%s/plots/TV_QED_ratio.xmg",ens.path.c_str()));
      dbvec_t(VV_rat/TV_rat).ave_err().write(combine("%s/plots/VV_TV_ratio.xmg",ens.path.c_str()));
      
      dboot_t Z_V,M_V,A_V,SL_V;
      two_pts_with_ins_ratio_fit(Z_V,M_V,A_V,SL_V,VV_LO,VV_QED,TH,ens.tmin[im],ens.tmax[im],
				 combine("%s/plots/VV_LO.xmg",ens.path.c_str()),
				 combine("%s/plots/VV_QED.xmg",ens.path.c_str()));
      dboot_t sqrt_Z_T_Z_V,M_T,A_T,SL_T;
      two_pts_with_ins_ratio_fit(sqrt_Z_T_Z_V,M_T,A_T,SL_T,TV_LO,TV_QED,TH,ens.tmin[im],ens.tmax[im],
				 combine("%s/plots/TV_LO.xmg",ens.path.c_str()),
				 combine("%s/plots/TV_QED.xmg",ens.path.c_str()));
      
      LO[iens]=integrate(VV_LO,ens.T,a);
      cout<<"amu: "<<LO[iens].ave_err()<<endl;
      
      QED[iens]=integrate(VV_QED,ens.T,a);
      cout<<"amu_QED: "<<QED[iens].ave_err()<<endl;
      
      ratio[iens]=QED[iens]/LO[iens];
      cout<<" Ratio: "<<dboot_t(ratio[iens]).ave_err()<<endl;
      
      //parameters to fit
      {
      	MnUserParameters pars;
      	size_t i_Z_V=0;pars.Add("Z_V",Z_V[0],Z_V.err());
      	size_t i_sqrt_Z_T_Z_V=1;pars.Add("sqrt_Z_T_Z_V",sqrt_Z_T_Z_V[0],sqrt_Z_T_Z_V.err());
      	size_t i_M=2;pars.Add("M_V",M_V[0],M_V.err());
      	size_t i_A_V=3;pars.Add("A_V",A_V[0],A_V.err());
      	size_t i_A_T=4;pars.Add("A_T",A_T[0],A_T.err());
      	size_t i_SL=5;pars.Add("SL",SL_V[0],SL_V.err());
	
      	size_t iel=0;
	const size_t dcT=size_t(2.0*TH/24);
      	auto x=vector_up_to<double>(VV_LO.size());
      	multi_ch2_t<dbvec_t> four_fit_obj({x,x,x,x},
      					  {ens.tmin[im],ens.tmin[im],ens.tmin[im]-dcT,ens.tmin[im]-dcT},
      					  {ens.tmax[im],ens.tmax[im],ens.tmax[im],ens.tmax[im]},
      					  {VV_LO,TV_LO,VV_QED/VV_LO,TV_QED/TV_LO},
      					  {two_pts_corr_fun_t(TH,1),two_pts_corr_fun_t(TH,-1),two_pts_corr_with_ins_fun_t(TH,1),two_pts_corr_with_ins_fun_t(TH,1)},
      					  [i_Z_V,i_sqrt_Z_T_Z_V,i_M,i_A_V,i_A_T,i_SL](const vector<double> &p,size_t icontr)
      					  {
      					    switch(icontr)
      					      {
      					      case 0:return vector<double>({p[i_Z_V],p[i_M]});break;
      					      case 1:return vector<double>({p[i_sqrt_Z_T_Z_V],p[i_M]});break;
      					      case 2:return vector<double>({p[i_M],p[i_A_V],p[i_SL]});break;
      					      case 3:return vector<double>({p[i_M],p[i_A_T],p[i_SL]});break;
      					      default: CRASH("Unknown contr %zu",icontr);return p;
      					      }
      					  },iel);
	
      	MnMigrad migrad(four_fit_obj,pars);
	
      	dboot_t Z_V_bis,Z_T_Z_V_bis,M,A_V_bis,A_T_bis,SL,ch2;
      	for(iel=0;iel<=nboots;iel++)
      	  {
      	    //minimize and print the result
      	    FunctionMinimum min=migrad();
      	    MinimumParameters par_min=min.Parameters();
      	    Z_V_bis[iel]=par_min.Vec()[i_Z_V];
      	    Z_T_Z_V_bis[iel]=par_min.Vec()[i_sqrt_Z_T_Z_V];
      	    M[iel]=par_min.Vec()[i_M];
      	    A_V_bis[iel]=par_min.Vec()[i_A_V];
      	    A_T_bis[iel]=par_min.Vec()[i_A_T];
      	    SL[iel]=par_min.Vec()[i_SL];
      	    ch2[iel]=min.Fval();
      	  }
	
      	cout<<"Z: "<<Z_V.ave_err()<<" "<<Z_V_bis.ave_err()<<endl;
      	cout<<"sqrt_Z_T_Z_V: "<<sqrt_Z_T_Z_V.ave_err()<<" "<<sqrt_Z_T_Z_V.ave_err()<<endl;
      	cout<<"M: "<<M_V.ave_err()<<" "<<M.ave_err()<<endl;
      	cout<<"A_V: "<<A_V.ave_err()<<" "<<A_V_bis.ave_err()<<endl;
      	cout<<"A_T: "<<A_T.ave_err()<<" "<<A_T_bis.ave_err()<<endl;
      	cout<<"SL: "<<SL_V.ave_err()<<" "<<SL.ave_err()<<endl;
      	cout<<"Ch2: "<<ch2.ave_err()<<endl;
	
	write_constant_fit_plot(combine("%s/plots/VV_LO_bis.xmg",ens.path.c_str()),ens.tmin[im],ens.tmax[im],M,eff_VV_LO);
	write_constant_fit_plot(combine("%s/plots/TV_LO_bis.xmg",ens.path.c_str()),ens.tmin[im]-dcT,ens.tmax[im],M,eff_TV_LO);
	write_fit_plot(combine("%s/plots/VV_QED_bis.xmg",ens.path.c_str()),ens.tmin[im],ens.tmax[im],[M,A_V_bis,SL,TH](double x)->dboot_t{return two_pts_corr_with_ins_ratio_fun(M,A_V_bis,SL,TH,x,1);},dbvec_t(VV_QED/VV_LO));
	write_fit_plot(combine("%s/plots/TV_QED_bis.xmg",ens.path.c_str()),ens.tmin[im]-dcT,ens.tmax[im],[M,A_T_bis,SL,TH](double x)->dboot_t{return two_pts_corr_with_ins_ratio_fun(M,A_T_bis,SL,TH,x,1);},dbvec_t(TV_QED/TV_LO));
	
	for(size_t t=ens.tmin[im];t<=TH;t++) VV_LO[t]=two_pts_corr_fun(Z_V_bis,M,TH,t,1);
	LO[iens]=integrate(VV_LO,ens.T,a);
	cout<<"amu: "<<LO[iens].ave_err()<<endl;
	
	for(size_t t=ens.tmin[im];t<=TH;t++) VV_QED[t]=two_pts_corr_with_ins_ratio_fun(M,A_V_bis,SL,TH,t,1)*VV_LO[t];
	QED[iens]=integrate(VV_QED,ens.T,a);
	cout<<"amu_QED: "<<QED[iens].ave_err()<<endl;
	
	ratio[iens]=QED[iens]/LO[iens];
	cout<<" Ratio: "<<dboot_t(ratio[iens]).ave_err()<<endl;
      }
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
