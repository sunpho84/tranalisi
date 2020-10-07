#include <tranalisi.hpp>

#include <phys_point.hpp>

constexpr double kappa=2.837297;
const double M2PiPhys=sqr(0.139);

string beta_of_ibeta[3]={"1.90","1.95","2.10"};
map<char,size_t> ibeta_of_id{{'A',0},{'B',1},{'D',2}};
bool isLoc;

static constexpr int propagateLatErr=0;
static constexpr int includeLog=1;
static constexpr int includeHands=0;

template <size_t B>
auto parseName(const string& name)
{
  string tempPref;
  string tempL;
  
  bool passed=false;
  
  bool discardedFirst=false;
  for(auto& t : name)
    {
      if(discardedFirst)
	{
	  if(t=='.')
	    passed=true;
	  else
	    (passed?tempL:tempPref)+=t;
	}
      discardedFirst=true;
    }
  
  return std::get<B>(make_tuple((size_t)atoi(tempL.c_str()),atoi(tempPref.c_str())*1e-4));
}

template <typename TL,
	  typename TM>
TM FVEuniv(const TL& L,const TM& M)
{
  return alpha_em*kappa/sqr(L)*(2+M*L);
}

struct perens_t
{
  const string name;
  
  const double am;
  
  const size_t T;
  
  const size_t Lfra;
  
  const size_t tmin;
  
  const size_t tmax;
  
  const size_t ibeta;
  
  djack_t aMPiFitted;
  
  djack_t daMPiFitted;
  
  djack_t daMPiFittedHandcuffs;
  
  /////////////////////////////////////////////////////////////////
  
  boot_init_t bi;
  
  dboot_t aM;
  
  dboot_t daM;
  
  dboot_t daMHandcuffs;
  
   bool in{false};
  
  template <typename T>
  struct ToFit
  {
    T Zv;
    
    T Za;
    
    T ainv;
    
    T M;
    
    T L;
    
    T dM;
    
    T FVE;
    
    T dM2;
    
    T da2M2;
    
    T dM2UnivCorrected;
    
    T da2M2UnivCorrected;
  };
  
  // template <typename T,
  // 	    typename...Oth>
  // static const T& g(const T& x,Oth...)
  // {
  //   return x;
  // }
  
  static const double& g(const dboot_t& x,const size_t& iboot)
  {
    return x[iboot];
  }
  
  template <typename T>
  static const T& g(const T& x)
  {
    return x;
  }
  
  static const double& g(const double& x,const size_t& iboot)
  {
    return x;
  }
  
  template <typename T,
	    typename...Oth>
  ToFit<T> getToFit(const T& ainv,const T& Zv,const T& Za,
		    Oth&&...oth) const
  {
    ToFit<T> out;
    
    out.ainv=ainv;
    out.Zv=Zv;
    out.M=g(aM,oth...)*ainv;
    out.L=g((double)Lfra,oth...)/ainv;
    if(isLoc)
      {
	out.dM=g(e2*daM/2,oth...)*sqr(Zv)*ainv;
	if(includeHands)
	  out.dM+=g(e2*daMHandcuffs/2,oth...)*sqr(Za)*ainv;
      }
    else
      out.dM=g(e2*daM/2,oth...)*ainv;
    
    out.FVE=FVEuniv(g(out.L,oth...),g(out.M,oth...));
    out.dM2=out.dM*2*out.M;
    out.da2M2=out.dM2/ainv/ainv;
    out.dM2UnivCorrected=out.dM2+out.FVE;
    out.da2M2UnivCorrected=out.dM2UnivCorrected/ainv/ainv;
    
    return out;
  }
  
  void convertToBoot(const size_t iult)
  {
    if(not in)
      bi.fill(234213423);
    in=true;
    aM=dboot_t(bi,aMPiFitted);
    cout<<"am errs: "<<aM.err()<<" "<<aMPiFitted.err()<<endl;
    daM=dboot_t(bi,daMPiFitted);
    cout<<"dam errs: "<<daM.err()<<" "<<daMPiFitted.err()<<endl;
    daMHandcuffs=dboot_t(bi,daMPiFittedHandcuffs);
  }
  
  perens_t(const string& name,const size_t& T,const size_t& tmin,const size_t& tmax) :
    name(name),am(parseName<1>(name)),T(T),Lfra(parseName<0>(name)),tmin(tmin),tmax(tmax),ibeta(ibeta_of_id[name[0]])
  {
  }
};

//! Reads the list of ensembles to be used
vector<perens_t> readEnsList()
{
  //! List returned
  vector<perens_t> output;
  
  //! Input file
  raw_file_t input("ensemble_list.txt","r");
  
  //! Total number of ensembles
  const size_t nEns=input.read<size_t>("NEns");
  output.reserve(nEns);
  
  set_njacks(input.read<size_t>("NJacks"));
  
  isLoc=input.read<size_t>("Local");
  
  for(size_t iEns=0;iEns<nEns;iEns++)
    {
      const string name=input.read<string>();
      const size_t T=input.read<size_t>();
      const size_t tmin=input.read<size_t>();
      const size_t tmax=T/2-5;
      
      output.emplace_back(name,T,tmin,tmax);
    }
  
  return output;
}

template <typename T,
	  typename TA,
	  typename TL>
T FVEansatz(const T& r2,const T& r2a2,const T& r2mpi,const TA& a2,const double m2,const TL& L)
{
  return 4.0*M_PI*alpha_em/3.0*(r2*(1+r2mpi*(m2-M2PiPhys)+r2a2*a2))*sqrt(m2)/(L*L*L);
}

const vector<vector<ave_err_t>> Zv_ae({{{0.587,0.004},{0.603,0.003},{0.655,0.003}},{{0.608,0.003},{0.614,0.002},{0.657,0.002}}});
const vector<vector<ave_err_t>> Za_ae({{{0.731,0.008},{0.737,0.005},{0.762,0.004}},{{0.703,0.002},{0.714,0.002},{0.752,0.002}}});
size_t iCurv,iCurvA2,iA1,if0,iD,iDm,iOffset,iR2,iR2a2,iR2Mpi,ia[nbeta],izv[nbeta],iza[nbeta];

template <typename T,
	  typename TA,
	  typename TL>
T ansatz(const vector<T>& p,const double& m2,const TA& a,const TL& L)
{
  const TA a2=a*a;
  const T& Curv0=p[iCurv];
  const T& CurvA2=p[iCurvA2];
  const T Curv=Curv0+a2*CurvA2;
  const T& A1=p[iA1];
  const T& f0=p[if0];
  const T& D=p[iD];
  const T& Dm=p[iDm];
  const T& Offset=p[iOffset];
  
  const T Q=4*Offset/pow(f0,4);
  const T W=m2/sqr((T)(4*M_PI*f0));
  
  const T uncorrected=e2*sqr(f0)*(Q-Curv*W*log(W)*includeLog)+A1*m2*alpha_em/(4*M_PI)+(D+Dm*m2)*a2;
  const T correction=FVEansatz(p[iR2],p[iR2a2],p[iR2Mpi],a2,m2,L);
  
  if(fit_debug)
    cout<<" {   "<<a2<<" "<<m2<<" "<<L<<" "<<uncorrected<<" "<<correction<<"   } "<<endl;
  return uncorrected+correction;
}

int main()
{
  auto ensList=readEnsList();
  
  loadUltimateInput("ultimate_input.txt");
  
  for(auto& ens : ensList)
    {
      mkdir(ens.name+"/plots");
      
      const string& dir=ens.name;
      const size_t& T=ens.T;
      
      auto read=
	[&dir,&T]
	(const string& corr)
	{
	  const djvec_t out=read_djvec(dir+"/jacks/"+corr,T).symmetrized();
	  
	  out.ave_err().write(dir+"/plots/"+corr+".xmg");
	  
	  return out;
	};
      
      const djvec_t P5P5_00=read("mes_contr_00");
      const djvec_t P5P5_LL=read("mes_contr_LL");
      const djvec_t handcuffs=includeHands?read("handcuffs"):(P5P5_LL*0.0);
      const djvec_t ratio_LL=P5P5_LL/P5P5_00;
      const djvec_t ratio_handcuffs=handcuffs/P5P5_00;
      const djvec_t eff_mass=effective_mass(P5P5_00);
      const djvec_t eff_slope_LL=effective_slope(ratio_LL,eff_mass,T/2);
      const djvec_t eff_slope_handcuffs=effective_slope(ratio_handcuffs,eff_mass,T/2);
      
      auto fit=
	[&ens,&dir]
	(const djvec_t& vec,const string& corr)
	{
	  return constant_fit(vec,ens.tmin,ens.tmax,dir+"/plots/"+corr+".xmg");
	};
      
      ens.aMPiFitted=fit(eff_mass,"eff_mass");
      ens.daMPiFitted=fit(eff_slope_LL,"eff_slope_LL");
      ens.daMPiFittedHandcuffs=fit(eff_slope_handcuffs,"eff_slope_handcuffs");
      
      cout<<ens.name<<" "<<smart_print(ens.aMPiFitted)
	  <<" "<<smart_print(ens.daMPiFitted)<<" "<<smart_print((djack_t)(2*(ens.aMPiFitted*ens.daMPiFitted)));
      if(includeHands)
	cout<<" "<<smart_print(ens.daMPiFittedHandcuffs);
      cout<<endl;
    }
  
  // for(size_t iens=2;iens<6;iens++)
  // 	ensList[iens].aMPiFitted=ensList[1].aMPiFitted;
  
  const size_t nUlt=2;
  for(size_t iult=0;iult<nUlt;iult++)
    {
      const size_t imethod=iult/4;
      
      dbvec_t Zv(nbeta);
      dbvec_t Za(nbeta);
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	{
	  Zv[ibeta].fill_gauss(Zv_ae[imethod][ibeta],234212+ibeta);
	  Za[ibeta].fill_gauss(Za_ae[imethod][ibeta],342342+ibeta);
	}
      
      for(auto& ens : ensList)
	ens.convertToBoot(iult);
      
      const size_t nFitPars=10+(1+isLoc*(1+includeHands))*nbeta;
      vector<dboot_t> pFit(nFitPars);
      vector<string> pname{"Curv","CurvA2","A1","f0","D","Dm","Offset","R2","R2a2","R2Mpi","ainv0","ainv1","ainv2","Zv0","Zv1","Zv2","Za0","Za1","Za2"};
      boot_fit_t fit;
      
      const double A1Guess=isLoc?-2e-2:-5.7;
      const double R2Guess=sqr(0.672/0.197);
      const double R2a2Guess=isLoc?1.40479:2.5;
      const double R2MpiGuess=isLoc?72.5:0.0;
      
      iCurv=fit.add_fit_par(pFit[0],pname[0],3.7,0.7);
      iCurvA2=fit.add_fit_par(pFit[1],pname[1],0,1e-6);
      iA1=fit.add_fit_par(pFit[2],pname[2],A1Guess,2e-2);
      if0=fit.add_self_fitted_point(pFit[3],pname[3],lat_par[iult].f0,-1);
      iD=fit.add_fit_par(pFit[4],pname[4],2e-3,1e-5);
      iDm=fit.add_fit_par(pFit[5],pname[5],0.0,0.001);
      iOffset=fit.add_fit_par(pFit[6],pname[6],1e-5,0.001);
      iR2=fit.add_fit_par(pFit[7],pname[7],R2Guess,1);
      iR2a2=fit.add_fit_par(pFit[8],pname[8],R2a2Guess,0.1);
      iR2Mpi=fit.add_fit_par(pFit[9],pname[9],R2MpiGuess,0.1);
      
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	ia[ibeta]=fit.add_self_fitted_point(pFit[10+ibeta],pname[10+ibeta],lat_par[iult].ainv[ibeta],-1);

      if(isLoc)
	{
	  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	    izv[ibeta]=fit.add_self_fitted_point(pFit[13+ibeta],pname[13+ibeta],Zv[ibeta],-1);
	  if(includeHands)
	    for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	      iza[ibeta]=fit.add_self_fitted_point(pFit[16+ibeta],pname[16+ibeta],Za[ibeta],-1);
	}
      
      ofstream data_out("data_out.txt");
      data_out.precision(16);
      auto pr=[&data_out](const string descr,const dboot_t a)
	      {
		data_out<<descr<<"\t"<<a[0]<<"\t"<<a.err()<<endl;
	      };
      
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	pr("ainv["+beta_of_ibeta[ibeta]+"]",lat_par[iult].ainv[ibeta]);
      pr("f0",lat_par[iult].f0);
      data_out<<
	"---------------------------\n"
	" Ensemble 	 pion_mass 	 DM^2 	        DM^2_err\n";
      
      fit.fix_par(iCurvA2);
      fit.fix_par(iR2);
      fit.fix_par(iR2Mpi);
      //fit.fix_par(iR2a2);
      
      //fit.fix_par(iOffset);
      if(isLoc)
	{
	  fit.fix_par(iD);
	  fit.fix_par(iDm);
	}
      else
	{
	  //fit.fix_par(iD);
	  fit.fix_par(iDm);
	}
      
      for(auto& ens : ensList)
	{
	  double err;
	  if(propagateLatErr)
	    err=ens.getToFit(lat_par[iult].ainv[ens.ibeta],Zv[ens.ibeta],Za[ens.ibeta]).da2M2.err();
	  else
	    {
	      const auto data=ens.getToFit(lat_par[iult].ainv[ens.ibeta],Zv[ens.ibeta],Za[ens.ibeta]);
	      err=data.da2M2.err();
	      
	      dboot_t y=data.dM/data.ainv*2;
	      pr(ens.name,y);
	      
	      if(isLoc)
		err*=sqr(Zv[ens.ibeta].ave());
	    }
	  
	  fit.add_point([=]
			(const vector<double>& p,int iboot)
			{
			  const auto data=ens.getToFit(p[ia[ens.ibeta]],p[izv[ens.ibeta]],p[iza[ens.ibeta]],iboot);
			  return data.da2M2;
			},
			[=]
			(const vector<double>& p,int iboot)
			{
			  //fit_debug=(iboot==nboots-1);
			  auto data=ens.getToFit(p[ia[ens.ibeta]],p[izv[ens.ibeta]],p[iza[ens.ibeta]],iboot);
			  
			  const double& ainv=data.ainv;
			  const double a=1/ainv;
			  const double& L=data.L;
			  const double x=sqr(data.M);
			  
			  return (ansatz(p,x,a,L)-FVEuniv(L,data.M))*a*a;
			},
			err);
	}
      
      // ens.getToFit(lat_par[iult].ainv[ens.ibeta]).dM2UnivCorrected.err());
      
      // for(size_t ibeta=0;ibeta<nbeta;ibeta++)
      // 	fit.fix_par(ia[ibeta]);
      // auto status=fit.fit();
      
      decltype(fit.fit()) status;
      
      for(int iit=0;iit<2;iit++)
      	{
	  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	    {
	      if(iit==0)
		fit.fix_par(ia[ibeta]);
	      else
		fit.unfix_par(ia[ibeta]);
	      
	      if(isLoc)
		{
		  if(iit==0)
		    {
		      fit.fix_par(izv[ibeta]);
		      fit.fix_par(iza[ibeta]);
		    }
		  else
		    {
		      fit.unfix_par(izv[ibeta]);
		      fit.unfix_par(iza[ibeta]);
		    }
		}
	    }
	  
	  status=fit.fit();
	}
      
      const dboot_t chi2fit=std::get<0>(status);
      
      /////////////////////////////////////////////////////////////////
      
      const size_t iboot_to_print=63;
      double chi2;
      chi2=0.0;
      auto add_pr_ch2=[&chi2](const dboot_t& teo,const dboot_t& num)
		   {
		     const double c=(teo[iboot_to_print]-num[iboot_to_print])/num.err();
		     const double r=sqr(c);
		     chi2+=r;
		     //cout<<"  "<<teo<<" "<<num<<"     "<<r.ave_err()<<"     "<<chi2.ave_err()<<endl;
		   };
      add_pr_ch2(pFit[3],lat_par[iult].f0);
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	add_pr_ch2(pFit[10+ibeta],lat_par[iult].ainv[ibeta]);
      
      ofstream file_out("parsandchi2.txt");
      file_out.precision(16);
      
#define PR(A) file_out<<#A<<": "<<A<<endl
      
      /// Fit pars and external inputs
      PR(pFit[ia[0]][iboot_to_print]);
      PR(pFit[ia[1]][iboot_to_print]);
      PR(pFit[ia[2]][iboot_to_print]);
      const double Curv0=pFit[iCurv][iboot_to_print];PR(Curv0);
      const double CurvA2=pFit[iCurvA2][iboot_to_print];PR(CurvA2);
      const double A1=pFit[iA1][iboot_to_print];PR(A1);
      const double f0=pFit[if0][iboot_to_print];PR(f0);
      const double D=pFit[iD][iboot_to_print];PR(D);
      const double Dm=pFit[iDm][iboot_to_print];PR(Dm);
      const double Offset=pFit[iOffset][iboot_to_print];PR(Offset);
      const double r2=pFit[iR2][iboot_to_print];PR(r2);
      const double r2m2=pFit[iR2Mpi][iboot_to_print];PR(r2m2);
      const double r2a2=pFit[iR2a2][iboot_to_print];PR(r2a2);
      
      file_out<<"# ===="<<endl;
      for(auto& ens : ensList)
	{
	  file_out<<"# Ens "<<ens.name<<" "<<ens.ibeta<<endl;
	  double err;
	  {
	    const dboot_t ainv=lat_par[iult].ainv[ens.ibeta];
	    const dboot_t M=ens.aM*ainv;
	    const dboot_t L=ens.Lfra/ainv;
	    const dboot_t dM=e2*ens.daM/2*ainv;
	    cout<<dM[iboot_to_print]<<endl;
	    //const dboot_t FVE=alpha_em*kappa/sqr(L)*(2+M*L);
	    const dboot_t dM2=dM*2*M;
	    //const dboot_t dM2UnivCorrected=dM2+FVE;
	    const dboot_t da2M2=dM2/ainv/ainv;
	    const djack_t da2M2proerr=e2*ens.daMPiFitted*2*ens.aMPiFitted;
	    err=da2M2proerr.err();
	  }
	  file_out<<"error: "<<err<<endl;
	  
	  const double ainv=pFit[ia[ens.ibeta]][iboot_to_print];
	  const double aM=ens.aM[iboot_to_print];PR(aM);
	  const size_t Lfra=ens.Lfra;file_out<<"Lfra: "<<Lfra<<endl;
	  const double daM=ens.daM[iboot_to_print];PR(daM);
	  
	  ///////////////////////
	  
	  const double M=aM*ainv;
	  const double L=Lfra/ainv;
	  const double dM=e2*daM/2*ainv;
	  const double dM2=dM*2*M;
	  const double da2M2=dM2/ainv/ainv;
	  //const double dM2UnivCorrected=dM2+FVE;
	  
	  const double a=1/ainv;
	  const double a2=a*a;
	  const double Curv=Curv0+a2*CurvA2;
	  const double m2=sqr(M);
	  
	  const double Q=4*Offset/pow(f0,4);
	  const double W=m2/sqr(4*M_PI*f0);
	  
	  const double uncorrected=e2*sqr(f0)*(Q-Curv*W*log(W)*includeLog)+A1*m2*alpha_em/(4*M_PI)+(D+Dm*m2)*a2;
	  const double univFVE=alpha_em*kappa/sqr(L)*(2+M*L);
	  const double nonunivFVE=4.0*M_PI*alpha_em/3.0*(r2*(1+r2m2*(m2-M2PiPhys)+r2a2*a2))*sqrt(m2)/(L*L*L);
	  
	  const double num=da2M2;
	  const double teo=(uncorrected-univFVE+nonunivFVE)*a*a;
	  const double r=(teo-num)/err;
	  const double c=sqr(r);
	  file_out<<"numerical_value: "<<num<<endl;
	  file_out<<"ansatz: "<<teo<<endl;
	  file_out<<"contribution: "<<c<<endl;
	  file_out<<"== =="<<endl;
	  chi2+=c;
	}
      file_out<<"chi2["<<iboot_to_print<<"]: "<<chi2<<" ori fit: "<<chi2fit[iboot_to_print]<<endl;
      
      fit_debug=false;
      /////////////////////////////////////////////////////////////////
      
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	{
	  const dboot_t& ainv_fit=pFit[ia[ibeta]];
	  const dboot_t& ainv_old=lat_par[iult].ainv[ibeta];
	  
	  const double ainv_diff=ainv_fit.ave()-ainv_old.ave();
	  const dboot_t ainv_sum=ainv_fit+ainv_old;
	  
	  cout<<"Significativity ainv beta "<<ibeta<<": "<<
	    ainv_diff/ainv_sum.err()<<endl;
	}
      
      if(isLoc and includeHands)
	for(auto& l : {make_tuple(izv,Zv),make_tuple(iza,Za)})
	  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	    {
	      auto i{get<0>(l)};
	      auto Z{get<1>(l)};
	      const dboot_t& z_fit=pFit[i[ibeta]];
	      
	      const double z_diff=z_fit.ave()-Z[ibeta].ave();
	      const dboot_t z_sum=z_fit+Z[ibeta];
	      
	      cout<<"Significativity z beta "<<ibeta<<": "<<
		z_diff/z_sum.err()<<endl;
	    }
      
      cout<<"Parameters:"<<endl;
      for(size_t ipar=0;ipar<nFitPars;ipar++)
	cout<<pname[ipar]<<" "<<pFit[ipar].ave_err()<<endl;
      
      /////////////////////////////////////////////////////////////////
      
      const string ult_plot_dir="plots/"+to_string(iult);
      mkdir(ult_plot_dir);
      grace_file_t dM2pi_plot(ult_plot_dir+"/dM2pi.xmg");
      grace_file_t da2M2pi_plot(ult_plot_dir+"/da2M2pi.xmg");
      grace_file_t daMpiHandcuffs_plot(ult_plot_dir+"/daMpiHandcuffs.xmg");
      
      for(auto p : {&dM2pi_plot,&da2M2pi_plot,
		      &daMpiHandcuffs_plot})
	{
	  p->set_settype(grace::XYDY);
	  p->set_xaxis_label("M\\s\\x\\p\\0\\N\\S2\\N (GeV\\S-1\\N)");
	  p->set_yaxis_label("2\\xd\\0M\\s\\x\\p\\0\\N M\\s\\x\\p\\0\\N (GeV\\S-2\\N)");
	}
      vector<grace::color_t> colors{grace::color_t::BLACK,grace::color_t::RED,grace::color_t::BLUE};
      
      const grace::line_style_t lineStyle[3]={grace::DASHED_LINE,grace::SHORT_DASHED_LINE,grace::CONTINUOUS_LINE};
      for(size_t FVEswitch=0;FVEswitch<3;FVEswitch++)
	for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	  {
	    for(auto p : {&dM2pi_plot,&da2M2pi_plot,&daMpiHandcuffs_plot})
	      {
		p->new_data_set();
		p->set_all_colors(colors[ibeta]);
		p->set_line_style(lineStyle[FVEswitch]);
		p->set_transparency(0.33+0.33*FVEswitch);
	      }
	    
	    for(auto& ens : ensList)
	      if(ens.ibeta==ibeta)
		{
		  const auto data=ens.getToFit(pFit[ia[ibeta]],pFit[izv[ibeta]],pFit[iza[ibeta]]);
		  const dboot_t x=sqr(data.M);
		  const vector<dboot_t> y=
		    {data.dM2,
		     data.dM2UnivCorrected,
		     data.dM2UnivCorrected-FVEansatz(pFit[iR2],pFit[iR2a2],pFit[iR2Mpi],1/sqr(pFit[ia[ibeta]]),x.ave(),data.L)};
		  
		  if(FVEswitch==0)
		    daMpiHandcuffs_plot.write_ave_err(x.ave(),((djack_t)(ens.daMPiFittedHandcuffs*data.ainv.ave())).ave_err());
		  
		  dM2pi_plot.write_ave_err(x.ave(),y[FVEswitch].ave_err());
		  da2M2pi_plot.write_ave_err(x.ave(),((dboot_t)(y[FVEswitch]/sqr((dboot_t)(pFit[ia[ibeta]]/lat_par[iult].ainv[ibeta].ave())))).ave_err());
		}
	    
	  }
      
      dboot_t Linf;
      Linf=1e30;
      
      dboot_t a0;
      a0=0.0;
      
      const double dM2PiPDG=1261.2e-6;
      const dboot_t dM2PiPhys=ansatz(pFit,M2PiPhys,a0,Linf);
      cout<<"dM2PiPhys: "<<smart_print(dM2PiPhys)<<endl;
      const dboot_t dM2PiDiff=dM2PiPhys-dM2PiPDG;
      cout<<"diff with phys: "<<smart_print(dM2PiDiff)<<endl;
      const ave_err_t dM2PiPublished(1137e-6,63e-6);
      
      vector<grace_file_t*> plots{&dM2pi_plot,&da2M2pi_plot};
      size_t plot_pow[2]={0,2};
      for(size_t iplot=0;iplot<2;iplot++)
	{
	  auto& p=*plots[iplot];
	  
	  if(not isLoc)
	    for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	      p.write_polygon([&pFit,&Linf,ibeta,iplot,plot_pow,iult]
			      (const double& x) -> dboot_t
			      {
				const dboot_t a=1/pFit[ia[ibeta]];
				
				return ansatz(pFit,x,a,Linf)*pow(a*lat_par[iult].ainv[ibeta].ave(),plot_pow[iplot]);
			      },1e-3,0.25,colors[ibeta]);
	  
	  p.write_polygon([&pFit,&a0,&Linf]
			   (const double& x) -> dboot_t
			   {
			     return ansatz(pFit,x,a0,Linf);
			   },1e-3,0.25,grace::YELLOW);
	  
	  p.set_symbol_fill_pattern(grace::FILLED_SYMBOL);
	  p.write_ave_err(M2PiPhys,dM2PiPhys.ave_err());
	  
	  p.new_data_set();
	  p.write_ave_err(M2PiPhys,dM2PiPublished);
	  
	  p.write_line([&dM2PiPDG](const double x){return dM2PiPDG;},1e-3,0.25,grace::GREEN4,2);
	}
      
      grace_file_t A40slice(ult_plot_dir+"/A40slice.xmg");
      grace_file_t A40handcuffsSlice(ult_plot_dir+"/A40handcuffsSlice.xmg");
      grace_file_t A40Mslice(ult_plot_dir+"/A40Mslice.xmg");
      
      for(auto& ens: ensList)
	if(ens.am==0.0040)
	  A40Mslice.write_ave_err(1.0/ens.Lfra,ens.aM.ave_err());
      
      /// handcuffs
      for(auto& ens: ensList)
	if(ens.am==0.0040)
	  A40handcuffsSlice.write_ave_err(1.0/ens.Lfra,ens.daMPiFitted.ave_err());
      A40handcuffsSlice.set_legend("exchange");
      A40handcuffsSlice.new_data_set(grace::GREEN4,grace::DIAMOND);
      for(auto& ens: ensList)
	if(ens.am==0.0040)
	  A40handcuffsSlice.write_ave_err(1.0/ens.Lfra,ens.daMPiFittedHandcuffs.ave_err());
      A40handcuffsSlice.set_legend("handcuffs");
      A40handcuffsSlice.new_data_set();
      
      for(int FSEflag=0;FSEflag<3;FSEflag++)
	{
	  A40slice.new_data_set();
	  A40slice.set_transparency(0.33+0.33*FSEflag);
	  
	  for(auto& ens: ensList)
	    if(ens.am==0.0040)
	      {
		const auto data=ens.getToFit(pFit[ia[0]],pFit[izv[0]],pFit[iza[0]]);
		
		djack_t y=e2*ens.daMPiFitted*2*ens.aMPiFitted/2;
		if(FSEflag>0)
		  y+=FVEuniv(ens.Lfra,ens.aMPiFitted);
		if(FSEflag>1)
		  y-=FVEansatz(pFit[iR2],pFit[iR2a2],pFit[iR2Mpi],1/sqr(pFit[ia[0]].ave()),sqr(data.M.ave()),data.L)/sqr(data.ainv);
	    
	    A40slice.write_ave_err(pow(ens.Lfra,-3),y.ave_err());
	  }
	}
      
      const double slice_mass=0.28;
      std::vector<perens_t*> slice(nbeta);
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	{
	  double closest_mass=0;
	  perens_t* ref_ens=nullptr;
	  
	  for(auto& ens : ensList)
	    if(ens.ibeta==ibeta)
	      {
		const auto data=ens.getToFit(pFit[ia[ibeta]],pFit[izv[ibeta]],pFit[iza[ibeta]]);
		if(fabs(data.M.ave()-slice_mass)<fabs(closest_mass-slice_mass))
		  {
		    closest_mass=data.M.ave();
		    ref_ens=&ens;
		  }
	      }
	  
	  if(ref_ens==nullptr)
	    CRASH("Not sliced");
	  
	  slice[ibeta]=ref_ens;
	  cout<<ref_ens->name<<endl;
	}
      
      grace_file_t daMpi_a2_plot(ult_plot_dir+"/daMpi_a2.xmg");
      daMpi_a2_plot.write_polygon([&pFit,&Linf,slice_mass]
				  (const double& x) -> dboot_t
				  {
				    return ansatz(pFit,slice_mass,x,Linf);
				  },1e-3,0.25,grace::YELLOW);
      
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	{
	  const auto& ainv=pFit[ia[ibeta]];
	  const auto& zv=pFit[izv[ibeta]];
	  const auto& za=pFit[iza[ibeta]];
	  const auto data=slice[ibeta]->getToFit(ainv,zv,za);
	  
	  const dboot_t y=data.dM2UnivCorrected-FVEansatz(pFit[iR2],pFit[iR2a2],pFit[iR2Mpi],1e10,sqr(data.M.ave()),data.L);
	  daMpi_a2_plot.write_ave_err(sqr(1/ainv.ave()),y.ave_err());
	}
    }
  
  return 0;
}
