#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <gm2_IB_common.hpp>

#include <omp.h>

int mul[10000];
double Zn20[10000],E2pi0[10000],Zn2[10000],E2pi[10000];


//evaluate multiplicity
void multply()
{
  int Nup=100;
  int Nterm=pow(Nup,2);
  int i2=0;

  for(int ix=-Nup;ix<Nup+1;ix++)
    for(int iy=-Nup;iy<Nup+1;iy++)
      for(int iz=-Nup;iz<Nup+1;iz++)
	{
	  i2=ix*ix+iy*iy+iz*iz;
	  if(i2<=Nterm) mul[i2]=mul[i2]+1;
	}
}

struct paramss_t
{
  int m2;
  double z2;
  paramss_t(int m2,double z2) : m2(m2),z2(z2){}
};

double ff(double t,void *params)
{
  int m2=((paramss_t*)params)->m2;
  double z2=((paramss_t*)params)->z2;

  double t0=0.5*M_PI;
  if(z2>2.0) t0=M_PI/z2;

  double sumt=0.0;
  for(int i2=0;i2<m2+1;i2++)
    sumt=sumt+mul[i2]*exp(-t*i2);

  double hp=0.0;
  for(int i2=1;i2<m2+4;i2++)
    hp=hp+mul[i2]*exp(-pow(M_PI,2)*i2/t);

  double h=0.0;
  for(int i2=m2+1;i2<m2+4;i2++)
    h=h+mul[i2]*exp(-t*i2);

  if(t<=t0)
    {
      return (exp(z2*t)*(1.0+hp)-1.0)*pow(M_PI/t,1.5)-exp(t*z2+log(sumt));
    }
  else
    {
      return exp(t*z2+log(h));
    }
}

//phi function
double phiml(int n2,double z2)
{
  if(z2<=1.0e-08 or abs(z2-n2)<1.0e-08) return 0.0;

  double tt0=0.5*M_PI;
  if(z2>2.0) tt0=M_PI/z2;

  int m2=n2+2;
  if(mul[m2]==0) m2=m2+1;

  double sum0=0.0;
  for(int i2=0;i2<m2+1;i2++)
    sum0=sum0+mul[i2]/(i2-z2);

  // gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
  gsl_integration_glfixed_table *tab=gsl_integration_glfixed_table_alloc(60);

  paramss_t param(m2,z2);

  gsl_function F;
  F.function=ff;
  F.params=&param;
  
  // //integrate
  // double result;
  // double abserr;
  // double start=0.0,epsabs=0.0,epsrel=1.0e-08;
  // gsl_integration_qagiu(&F,start,epsabs,epsrel,1000,w,&result,&abserr);
  
  // gsl_integration_workspace_free(w);

  double result;
  double lb=0.0,ub=10*M_PI;
  result=gsl_integration_glfixed(&F,lb,ub,tab);

  gsl_integration_glfixed_table_free(tab);

  double res;
  res=result-2.0*M_PI*sqrt(M_PI/tt0)+sum0;

  double phi=atan(-2.0*pow(M_PI,2)*sqrt(z2)/res);
  
  return phi;
}

//Gounaris-Sakurai model
double GSmodel(double k,double M_PS,double M_V,double g_Vpipi,int g_Vpipi_free)
{
  double kV=0.25*M_V*M_V-M_PS*M_PS;
  
  if(k<=1.0e-08 or kV<=1.0e-08) return 0.0;

  kV=sqrt(kV);
  double hV=2.0*kV*log((M_V+2.0*kV)/(2.0*M_PS))/(M_PI*M_V);
  double hpV=(1.0/M_PI+hV*M_PS*M_PS/(kV*kV))/M_V;
  double omg=2.0*sqrt(M_PS*M_PS+k*k);
  double homg=2.0*k*log((omg+2.0*k)/(2.0*M_PS))/(M_PI*omg);
  double b;
  //double b=-(4.0*kV*kV*kV/(M_V*M_V*G_V)+hV+2.0*kV*kV*hpV/M_V);
  if(g_Vpipi_free==0)
    {
      b=-(4.0*0.638*g_Vpipi+hV+2.0*kV*kV*hpV/M_V);
    }
  else
    {
      b=-(24.0*M_PI/(g_Vpipi*g_Vpipi)+hV+2.0*kV*kV*hpV/M_V);
    }
  double delta11=omg*(k*k*homg-kV*kV*hV+b*(k*k-kV*kV))/(k*k*k);

  double sign=1.0;
  if(delta11<0.0) sign=-1.0;

  delta11=asin(sign/sqrt(1.0+delta11*delta11));
  
  return delta11;
}

double rootlm(int n2,double zeta2,int L,double M_PS,double M_V,double g_Vpipi,int g_Vpipi_free)
{
  if(zeta2<=1.0e-08) return 0.0;

  double phi=phiml(n2,zeta2);

  double k=2.0*M_PI*sqrt(zeta2)/L;

  double delta11=GSmodel(k,M_PS,M_V,g_Vpipi,g_Vpipi_free);

  double root=delta11+phi;
  
  return root;
}

//pion form factor
double F2piGS(double omg,double M_PS,double M_V,double g_Vpipi,int g_Vpipi_free)
{
  double kV=0.25*M_V*M_V-M_PS*M_PS;

  if(kV<=1.0e-08) return 0.0;

  kV=sqrt(kV);
  double hV=2.0*kV*log((M_V+2.0*kV)/(2.0*M_PS))/(M_PI*M_V);
  double hpV=(1.0/M_PI+hV*M_PS*M_PS/(kV*kV))/M_V;
  double k=0.25*omg*omg-M_PS*M_PS;

  if(k<=1.0e-08) return 0.0;

  k=sqrt(k);
  double homg=2.0*k*log((omg+2.0*k)/(2.0*M_PS))/(M_PI*omg);
  double b;
  //double b=-(4.0*kV*kV*kV/(M_V*M_V*G_V)+hV+2.0*kV*kV*hpV/M_V);
  if(g_Vpipi_free==0)
    {
      b=-(4.0*0.638*g_Vpipi+hV+2.0*kV*kV*hpV/M_V);
    }
  else
    {
      b=-(24.0*M_PI/(g_Vpipi*g_Vpipi)+hV+2.0*kV*kV*hpV/M_V);
    }
  double f0=-(M_PS*M_PS/M_PI+kV*kV*hV+0.25*b*M_V*M_V);
  double delta11=omg*(k*k*homg-kV*kV*hV+b*(k*k-kV*kV))/(k*k*k);
  double F2pi=omg*abs(f0)/(k*k*k*sqrt(1.0+delta11*delta11));

  return F2pi;
}

void states(int Nstati,int L,double M_PS,double M_V,double g_Vpipi,int g_Vpipi_free)
{
  double conv=L/(2.0*M_PI);

  double stepz=0.0001;

  // FILE* pFile;

  // pFile=fopen("energy_2pi_test.dat","w+");
  // fprintf(pFile,"M_PS(GeV) = %f  M_V(GeV) = %f\n\n",M_PS,M_V);
  // fprintf(pFile,"L(fm) = %i\n\n",L);
  // fprintf(pFile,"Nstates = %i\n\n\n\n",Nstati);
  // fprintf(pFile,"n2    mul    k0(GeV)        E2pi0(GeV)        Z0^2/mul        z2        k(GeV)        E2pi(GeV)        Delta(GeV)        Z^2/mul        delta11        \n\n");
  // fclose(pFile);

  for(int n2=1;n2<Nstati+1;n2++)
    {
      Zn20[n2]=0.0;
      E2pi0[n2]=0.0;
      Zn2[n2]=0.0;
      E2pi[n2]=0.0;
      if(mul[n2]==0) continue;
      double k0=sqrt(n2)/conv;
      E2pi0[n2]=2.0*sqrt(M_PS*M_PS+k0*k0);
      double phim=phiml(n2,n2-stepz);
      double phip=phiml(n2,n2+stepz);      
      double dllk0=n2*(phip-phim)/stepz;
      double F2pi0=F2piGS(E2pi0[n2],M_PS,M_V,g_Vpipi,g_Vpipi_free);
      Zn20[n2]=2.0*F2pi0*F2pi0*k0*k0*k0*k0*k0/(3.0*M_PI*E2pi0[n2]*E2pi0[n2]*dllk0);
      double step=0.01;
      double z2root=n2+0.5*step;
      double root=rootlm(n2,z2root,L,M_PS,M_V,g_Vpipi,g_Vpipi_free);
      while(abs(root)>=1.0e-08)
  	{
  	  double root0=root;
  	  z2root=z2root-step;
  	  root=rootlm(n2,z2root,L,M_PS,M_V,g_Vpipi,g_Vpipi_free);
  	  if((root0*root)<0.0)
  	    {
  	      if(abs(root-root0)>1.0) continue;
  	      step=-0.5*step;
  	    }
  	}
      if(z2root<1.0e-08) z2root=n2;
      double kroot=sqrt(z2root)/conv;
      E2pi[n2]=2.0*sqrt(M_PS*M_PS+kroot*kroot);
      phim=phiml(n2,z2root-stepz);
      phip=phiml(n2,z2root+stepz);
      double delta11m=GSmodel(kroot*(1.0-stepz),M_PS,M_V,g_Vpipi,g_Vpipi_free);
      double delta11p=GSmodel(kroot*(1.0+stepz),M_PS,M_V,g_Vpipi,g_Vpipi_free);
      double dllk=z2root*(phip-phim)/stepz+(delta11p-delta11m)/(2.0*stepz);
      double F2pi=F2piGS(E2pi[n2],M_PS,M_V,g_Vpipi,g_Vpipi_free);
      Zn2[n2]=2.0*F2pi*F2pi*kroot*kroot*kroot*kroot*kroot/(3.0*M_PI*E2pi[n2]*E2pi[n2]*dllk);

      // cout<<"n2: "<<n2<<"  "<<"Zn2: "<<Zn2[n2]*pow(0.19731/(0.474/5.31),3)<<"  "<<"E2pi: "<<E2pi[n2]*0.19731/(0.474/5.31)<<endl;

      // double delta11=GSmodel(kroot,M_PS,M_V,G_V);

      // pFile=fopen("energy_2pi_test.dat","a+");
      // fprintf(pFile,"%i    %i    %lg        %lg        %lg        %lg        %lg        %lg        %lg        %lg        %lg        \n",n2,mul[n2],k0,E2pi0[n2],Zn20[n2],z2root,kroot,E2pi[n2],(E2pi[n2]-E2pi0[n2]),Zn2[n2],delta11);
      // fclose(pFile);
    }
}

double V_corr_ll(double t,int Nstati,int L,double M_PS,double M_V,double g_Vpipi,int g_Vpipi_free)
{  
  double sum=0;
  
  states(Nstati,L,M_PS,abs(M_V),abs(g_Vpipi),g_Vpipi_free);
  
  for(int i=1;i<Nstati+1;i++)
    {
      sum=sum+Zn2[i]*exp(-E2pi[i]*t);
    }

  return sum;
}

int main(int narg,char **arg)
{
  //we have to limit the number of threads
  omp_set_num_threads(1);
  int start=time(0);
  
  gm2_initialize(narg,arg);

  const double qu=2.0/3.0;
  const double qd=1.0/3.0;

  index_t ind_tab({{"Input",ninput_an}});
  
  vector<djvec_t> jPP_LO(nens_used),jVV_LO(nens_used);
  djvec_t jZ_P(nens_used),jM_P(nens_used),jZ_V(nens_used),jZ_V_nocharge(nens_used),jM_V(nens_used),jg_Vpipi(nens_used),jM_V_ll(nens_used),jg_Vpipi_ll(nens_used),jM_V_ll_s(nens_used),jg_Vpipi_ll_s(nens_used);
  valarray<djack_t> ja_mu_lat_tot(nens_used),ja_mu_lat_res(nens_used),ja_mu_twopion_fit(nens_used),ja_mu_twopion_fit_s(nens_used);

  multply();

  for(size_t iens=3;iens<4;iens++)
      {
	ens_data_t &ens=ens_data[iens];
	size_t TH=ens.T/2;
      
	size_t ib=ens.ib;

	size_t L=ens.L;

	double Za=Za_ae[0][ib].ave();
	djack_t a;a=1/lat_par[0].ainv[ib].ave();

	//load LO
	jPP_LO[iens]=read_PP("00",ens,im,1,RE);
	jVV_LO[iens]=(sqr(qu)+sqr(qd))*read_VV("00",ens,im,1,RE)*sqr(Za);

	two_pts_fit(jZ_P[iens],jM_P[iens],jPP_LO[iens],TH,ens.tmin[im],ens.tmax[im],combine("%s/plots/PP_LO_emass.xmg",ens.path.c_str()));
	two_pts_fit(jZ_V[iens],jM_V[iens],jVV_LO[iens],TH,ens.tmin_VV,ens.tmax_VV,combine("%s/plots/VV_LO_emass.xmg",ens.path.c_str()));

	jPP_LO[iens].ave_err().write(ens.path+"/plots/jPP_LO_corr.xmg");
	jVV_LO[iens].ave_err().write(ens.path+"/plots/jVV_LO_corr.xmg");

	jZ_V_nocharge[iens]=jZ_V[iens]*9.0/5.0;
	jg_Vpipi[iens]=sqrt(3.0)*jM_V[iens]*jM_V[iens]/sqrt(jZ_V_nocharge[iens]);

	cout<<"---------------------"<<iens<<" "<<ens.path<<"------------------"<<endl<<endl<<endl;
	cout<<"aM_P: "<<jM_P[iens].ave_err()<<endl;
	cout<<"aM_V: "<<jM_V[iens].ave_err()<<endl;
	cout<<"a^4Z_V: "<<jZ_V_nocharge[iens].ave_err()<<endl;
	cout<<"g_Vpipi: "<<jg_Vpipi[iens].ave_err()<<endl;

	double Emax=4.0;
	double conv=L*a.ave()/(2.0*M_PI);
	double k2max=max(1.0e-08,0.25*pow(Emax,2)-jM_P[iens].ave()*jM_P[iens].ave()*a.ave()*a.ave());
	double zmax=sqrt(k2max)*conv;
  
	int Nstati=pow(zmax,2)+1.0e-08;

	if(Nstati<2) Nstati=2;
	if(mul[Nstati]==0) Nstati=Nstati+1;

	//fits
	vector<double> t_x(ens.tmax_VV-ens.tmin_VV+1);

	for(size_t t=ens.tmin_VV;t<ens.tmax_VV+1;t++)
	  t_x[t-ens.tmin_VV]=(double)t;

	// cout<<"---------------------"<<"aM_V_single_fit"<<"------------------"<<endl<<endl;
	
	// jack_fit_t fit_ll_M;
	// fit_ll_M.add_fit_par_limits(jM_V_ll_s[iens],"aM_V",jM_V[iens].ave_err(),1.0e-08,10.0);
	// for(size_t t=ens.tmin_VV;t<ens.tmax_VV+1;t++)
	//   fit_ll_M.add_point(jVV_LO[iens][t],[&t_x,t,ens,Nstati,L,iens,jM_P](const vector<double> &p,int iel){return V_corr_ll(t_x[t-ens.tmin_VV],Nstati,L,jM_P[iens][iel],p[0],1.0,0);});
	
	// fit_ll_M.fit();

	// cout<<"aM_V_fit_s: "<<jM_V_ll_s[iens].ave_err()<<endl;

	// cout<<"---------------------"<<"g_Vpipi_single_fit"<<"------------------"<<endl<<endl;

	// jack_fit_t fit_ll_g;
	// fit_ll_g.add_fit_par_limits(jg_Vpipi_ll_s[iens],"g_Vpipi",jg_Vpipi[iens].ave_err(),1.0,20.0);
	// for(size_t t=ens.tmin_VV;t<ens.tmax_VV+1;t++)
	//   fit_ll_g.add_point(jVV_LO[iens][t],[&t_x,t,ens,Nstati,L,iens,jM_P,jM_V_ll_s](const vector<double> &p,int iel){return V_corr_ll(t_x[t-ens.tmin_VV],Nstati,L,jM_P[iens][iel],jM_V_ll_s[iens][iel],p[0],1);});
	
	// fit_ll_g.fit();

	// cout<<"g_Vpipi_fit_s: "<<jg_Vpipi_ll_s[iens].ave_err()<<endl;

	cout<<"---------------------"<<"GLOBAL FIT"<<"------------------"<<endl<<endl;
		
	jack_fit_t fit_ll;
	fit_ll.add_fit_par_limits(jM_V_ll[iens],"aM_V",jM_V[iens].ave_err(),1.0e-08,10.0);
	fit_ll.add_fit_par_limits(jg_Vpipi_ll[iens],"g_Vpipi",jg_Vpipi[iens].ave_err(),1.0,20.0);
	// fit_ll.add_self_fitted_point(jg_Vpipi_ll[iens],"jg_Vpipi",jg_Vpipi[iens].ave_err());
	for(size_t t=ens.tmin_VV;t<ens.tmax_VV+1;t++)
	  fit_ll.add_point(jVV_LO[iens][t],[&t_x,t,ens,Nstati,L,iens,jM_P](const vector<double> &p,int iel){return V_corr_ll(t_x[t-ens.tmin_VV],Nstati,L,jM_P[iens][iel],p[0],p[1],1);});
	
	fit_ll.fit();

	cout<<"aM_V_fit: "<<jM_V_ll[iens].ave_err()<<endl;
	cout<<"g_Vpipi_fit: "<<jg_Vpipi_ll[iens].ave_err()<<endl;

	// valarray<djack_t> jVV_LO_twopion_ll(TH+1),jVV_LO_twopion_ll_s(TH+1);
	
	// //a_mu
	// for(size_t t=0;t<TH+1;t++)
	//   for(size_t i=0;i<jVV_LO_twopion_ll[t].size();i++)
	//     jVV_LO_twopion_ll[t][i]=V_corr_ll(t,Nstati,L,jM_P[iens][i],jM_V_ll[iens][i],jg_Vpipi_ll[iens][i],1);

	// for(size_t t=0;t<TH+1;t++)
	//   for(size_t i=0;i<jVV_LO_twopion_ll_s[t].size();i++)
	//     jVV_LO_twopion_ll_s[t][i]=V_corr_ll(t,Nstati,L,jM_P[iens][i],jM_V_ll_s[iens][i],jg_Vpipi_ll_s[iens][i],1);

	// ja_mu_lat_tot[iens]=integrate_jcorr_light_times_kern(jVV_LO[iens],ens.T,a,0,TH);
	// ja_mu_lat_res[iens]=integrate_jcorr_light_times_kern(jVV_LO[iens],ens.T,a,3,ens.tmin_VV-1);
	// ja_mu_twopion_fit[iens]=integrate_jcorr_light_times_kern(jVV_LO_twopion_ll,ens.T,a,ens.tmin_VV,ens.tmax_VV);
	// ja_mu_twopion_fit_s[iens]=integrate_jcorr_light_times_kern(jVV_LO_twopion_ll_s,ens.T,a,ens.tmin_VV,ens.tmax_VV);

	// cout<<"a_mu lattice total: "<<ja_mu_lat_tot[iens].ave_err()<<endl;
	// cout<<"a_mu lattice res: "<<ja_mu_lat_res[iens].ave_err()<<endl;
	// cout<<"a_mu two pion global fit: "<<ja_mu_twopion_fit[iens].ave_err()<<endl;
	// cout<<"a_mu two pion single fit: "<<ja_mu_twopion_fit_s[iens].ave_err()<<endl;

	// //print results
	// cout<<"---------------------"<<iens<<" "<<ens.path<<"------------------"<<endl<<endl<<endl;
	
	// cout<<"---------------------"<<"aM_V_single_fit"<<"------------------"<<endl;
	// for(size_t i=0;i<jM_V_ll_s[iens].size();i++)
	//   cout<<jM_V_ll_s[iens][i]<<endl;

	// cout<<"average: "<<jM_V_ll_s[iens].ave_err()<<endl<<endl;

	// cout<<"---------------------"<<"g_Vpipi_single_fit"<<"------------------"<<endl;
	// for(size_t i=0;i<jg_Vpipi_ll_s[iens].size();i++)
	//   cout<<jg_Vpipi_ll_s[iens][i]<<endl;

	// cout<<"average: "<<jg_Vpipi_ll_s[iens].ave_err()<<endl<<endl;

	// cout<<"---------------------"<<"aM_V_global_fit"<<"------------------"<<endl;
	// for(size_t i=0;i<jM_V_ll[iens].size();i++)
	//   cout<<jM_V_ll[iens][i]<<endl;

	// cout<<"average: "<<jM_V_ll[iens].ave_err()<<endl<<endl;

	// cout<<"---------------------"<<"g_Vpipi_global_fit"<<"------------------"<<endl;
	// for(size_t i=0;i<jg_Vpipi_ll[iens].size();i++)
	//   cout<<jg_Vpipi_ll[iens][i]<<endl;

	// cout<<"average: "<<jg_Vpipi_ll[iens].ave_err()<<endl<<endl;

	// cout<<"---------------------"<<"a_mu_lat_tot"<<"------------------"<<endl;
	// cout<<"average: "<<ja_mu_lat_tot[iens].ave_err()<<endl<<endl;

	// cout<<"---------------------"<<"a_mu_lat_res"<<"------------------"<<endl;
	// cout<<"average: "<<ja_mu_lat_res[iens].ave_err()<<endl<<endl;

	// cout<<"---------------------"<<"a_mu_twopion_fit"<<"------------------"<<endl;
	// cout<<"average: "<<ja_mu_twopion_fit[iens].ave_err()<<endl<<endl;

	// cout<<"---------------------"<<"a_mu_twopion_fit_s"<<"------------------"<<endl;
	// cout<<"average: "<<ja_mu_twopion_fit_s[iens].ave_err()<<endl<<endl;


	// for(size_t t=0;t<TH+1;t++)
	//   cout<<V_corr_ll(t,Nstati,L,jM_P[iens].ave(),jM_V[iens].ave(),1.0)<<endl;
	// for(size_t t=0;t<TH+1;t++)
	//   cout<<t<<endl;
	    
	  
	//bootstrap
	for(size_t itab=0;itab<ind_tab.max();itab++)
	{
	  vector<size_t> comps=ind_tab(itab);
	  size_t input_an_id=comps[0];
	  
	  boot_init_t &bi=jack_index[input_an_id][ens.iult];
      
	  dboot_t a=1/lat_par[input_an_id].ainv[ib];

	  dbvec_t PP_LO(bi,jPP_LO[iens]),VV_LO(bi,jVV_LO[iens]);
	}
      }
  
  cout<<endl<<"Total time: "<<time(0)-start<<endl;
  
  return 0;
}

// const double sizeL=23.949;
  
  // double conv=sizeL/(2.0*M_PI*0.19731);
  // double kV=sqrt(max(1.0e-08,0.25*pow(M_V,2)-pow(M_PS,2)));
  
  // int Nstati=200;

  // multply();

  // FILE* pFile;

  // pFile=fopen("phi.dat","w+");
  // fprintf(pFile,"M_PS(GeV) = %f  M_V(GeV) = %f  kV(GeV) = %f  L(fm) = %f\n\n\n\n",M_PS,M_V,kV,sizeL);
  // fclose(pFile);

  // for(int n2=0;n2<Nstati+1;n2++)
  //   {
  //     if(mul[n2]==0) continue;

  //     pFile=fopen("phi.dat","a+");
  //     fprintf(pFile,"n2 = %i\n",n2);
  //     fprintf(pFile,"k        zeta^2        delta11        phi        root        \n\n");
  //     fclose(pFile);

  //     double z2min=n2-0.95;
  //     if(z2min<=0.0) z2min=0.0;
  //     double z2max=n2+0.95;
  //     int nstep=190;
  //     double stepz2=(z2max-z2min)/((float)nstep);
  //     double zeta2=z2min;
  //     for(int i=0;i<nstep+1;i++)
  // 	{
  // 	  double phi=phiml(n2,zeta2);
  // 	  double k=sqrt(zeta2)/conv;
  // 	  double delta11=GSmodel(k);
  // 	  double root=delta11+phi;
  // 	  pFile=fopen("phi.dat","a+");
  // 	  fprintf(pFile,"%lg        %lg        %lg        %lg        %lg        \n",k,zeta2,delta11,phi,root);
  // 	  fclose(pFile);
  // 	  zeta2=zeta2+stepz2;
  // 	}
  //   }

// FILE* pFile;

//   pFile=fopen("energy_2pi.dat","w+");
//   fprintf(pFile,"M_PS(GeV) = %f  M_V(GeV) = %f  kV(GeV) = %f\n\n",M_PS,M_V,kV);
//   fprintf(pFile,"L(fm) = %f  M_PS*L = %f\n\n",sizeL,M_PSL);
//   fprintf(pFile,"Emax(GeV) = %f  Nstates = %i\n\n\n\n",Emax,Nstati);
//   fprintf(pFile,"n2    mul    k0(GeV)        E2pi0(GeV)        Z0^2/mul        z2        k(GeV)        E2pi(GeV)        Delta(GeV)        Z^2/mul        delta11        \n\n");
//   fclose(pFile);

// pFile=fopen("energy_2pi.dat","a+");
//       fprintf(pFile,"%i    %i    %lg        %lg        %lg        %lg        %lg        %lg        %lg        %lg        %lg        \n",n2,mul[n2],k0,E2pi0,Zn20,z2root,kroot,E2pi,(E2pi-E2pi0),Zn2,delta11);
//       fclose(pFile);
