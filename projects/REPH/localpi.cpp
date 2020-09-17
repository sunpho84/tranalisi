#include <tranalisi.hpp>

#include <phys_point.hpp>

constexpr double kappa=2.837297;
const double M2PiPhys=sqr(0.139);

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
    
    T dM2UnivCorrected;
    
    T daM2UnivCorrected;
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
    out.dM2UnivCorrected=out.dM2+out.FVE;
    out.daM2UnivCorrected=out.dM2UnivCorrected/ainv/ainv;
    
    return out;
  }
  
  void convertToBoot(const size_t iult)
  {
    bi=jack_index[iult][0];
    aM=dboot_t(bi,aMPiFitted);
    daM=dboot_t(bi,daMPiFitted);
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
  return (4.0*M_PI*alpha_em/3.0*(r2+r2mpi*(m2-M2PiPhys))+r2a2*a2)*sqrt(m2)/(L*L*L);
}

const vector<vector<ave_err_t>> Zv_ae({{{0.587,0.004},{0.603,0.003},{0.655,0.003}},{{0.608,0.003},{0.614,0.002},{0.657,0.002}}});
const vector<vector<ave_err_t>> Za_ae({{{0.731,0.008},{0.737,0.005},{0.762,0.004}},{{0.703,0.002},{0.714,0.002},{0.752,0.002}}});
size_t iC,iCa,iA1,if0,iD,iDm,iQshift,iR2,iR2a2,iR2Mpi,ia[nbeta],izv[nbeta],iza[nbeta];

template <typename T,
	  typename TA,
	  typename TL>
T ansatz(const vector<T>& p,const double& m2,const TA& a,const TL& L)
{
  const TA a2=a*a;
  const T& C0=p[iC];
  const T& Ca=p[iCa];
  const T C=C0+a2*Ca;
  const T& A1=p[iA1];
  const T& f0=p[if0];
  const T& D=p[iD];
  const T& Dm=p[iDm];
  const T& Qshift=p[iQshift];
  
  const T Q=4*C/pow(f0,4);
  const T W=m2/sqr((T)(4*M_PI*f0));
  
  T uncorrected=e2*sqr(f0)*(Q+Qshift-(3+4*Q)*W*log(W)*includeLog)+A1*m2*alpha_em/(4*M_PI)+(D+Dm*m2)*a2;
  
  return uncorrected+FVEansatz(p[iR2],p[iR2a2],p[iR2Mpi],a2,m2,L);
}

int main()
{
  set_njacks(15);
  loadUltimateInput("ultimate_input.txt");
  
  auto ensList=readEnsList();
  
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
	  <<" "<<smart_print(ens.daMPiFitted);
      if(includeHands)
	cout<<" "<<smart_print(ens.daMPiFittedHandcuffs);
      cout<<endl;
    }
  
  // for(size_t iens=2;iens<6;iens++)
  // 	ensList[iens].aMPiFitted=ensList[1].aMPiFitted;
  
  //fit_debug=true;
  const size_t nUlt=1;
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
      
      const size_t nFitPars=10+3*nbeta;
      vector<dboot_t> pFit(nFitPars);
      vector<string> pname{"C","Ca","A1","f0","D","Dm","Qshift","R2","R2a2","R2Mpi","ainv0","ainv1","ainv2","Zv0","Zv1","Zv2","Za0","Za1","Za2"};
      boot_fit_t fit;
      
      const double A1Guess=isLoc?-3.5:-5.7;
      const double R2Guess=11;
      const double R2a2Guess=isLoc?1.40479:1.8;
      const double R2MpiGuess=isLoc?72.5:46.0;
      
      iC=fit.add_fit_par(pFit[0],pname[0],4e-05,1e-6);
      iCa=fit.add_fit_par(pFit[1],pname[1],0,1e-6);
      iA1=fit.add_fit_par(pFit[2],pname[2],A1Guess,0.2);
      if0=fit.add_self_fitted_point(pFit[3],pname[3],lat_par[iult].f0,-1);
      iD=fit.add_fit_par(pFit[4],pname[4],0.0,1e-5);
      iDm=fit.add_fit_par(pFit[5],pname[5],0.0,0.001);
      iQshift=fit.add_fit_par(pFit[6],pname[6],0.0,0.001);
      iR2=fit.add_fit_par(pFit[7],pname[7],R2Guess,1);
      iR2a2=fit.add_fit_par(pFit[8],pname[8],R2a2Guess,0.1);
      iR2Mpi=fit.add_fit_par(pFit[9],pname[9],R2MpiGuess,0.1);
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	ia[ibeta]=fit.add_self_fitted_point(pFit[10+ibeta],pname[10+ibeta],lat_par[iult].ainv[ibeta],-1);
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	izv[ibeta]=fit.add_self_fitted_point(pFit[13+ibeta],pname[13+ibeta],Zv[ibeta],-1);
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	iza[ibeta]=fit.add_self_fitted_point(pFit[16+ibeta],pname[16+ibeta],Za[ibeta],-1);
      
      fit.fix_par(iCa);
      fit.fix_par(iR2);
      //fit.fix_par(iFVEnonUniv2);
      
      fit.fix_par(iQshift);
      if(isLoc)
	{
	  fit.fix_par(iD);
	  fit.fix_par(iDm);
	}
      
      for(auto& ens : ensList)
	{
	  double err;
	  if(propagateLatErr)
	    err=ens.getToFit(lat_par[iult].ainv[ens.ibeta],Zv[ens.ibeta],Za[ens.ibeta]).dM2UnivCorrected.err();
	  else
	    {
	      err=ens.getToFit(lat_par[iult].ainv[ens.ibeta],Zv[ens.ibeta],Za[ens.ibeta]).daM2UnivCorrected.err();
	      err*=sqr(lat_par[iult].ainv[ens.ibeta].ave());
	      
	      if(isLoc)
		err*=sqr(Zv[ens.ibeta].ave());
	    }
	  
	  
	  fit.add_point([=]
			(const vector<double>& p,int iboot)
			{
			  const auto data=ens.getToFit(p[ia[ens.ibeta]],p[izv[ens.ibeta]],p[iza[ens.ibeta]],iboot);
			  return data.dM2UnivCorrected;
			},
			[=]
			(const vector<double>& p,int iboot)
			{
			  auto data=ens.getToFit(p[ia[ens.ibeta]],p[izv[ens.ibeta]],p[iza[ens.ibeta]],iboot);
			  
			  const double& ainv=data.ainv;
			  const double a=1/ainv;
			  const double& L=data.L;
			  const double x=sqr(data.M);
			  
			  return ansatz(p,x,a,L);
			},
			err);
	}
      
      // ens.getToFit(lat_par[iult].ainv[ens.ibeta]).dM2UnivCorrected.err());
      
      // for(size_t ibeta=0;ibeta<nbeta;ibeta++)
      // 	fit.fix_par(ia[ibeta]);
      // auto status=fit.fit();
	  
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
	  
	    auto status=fit.fit();
	}
      
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	{
	  dboot_t& ainv_fit=pFit[ia[ibeta]];
	  dboot_t& ainv_old=lat_par[iult].ainv[ibeta];
	  
	  dboot_t ainv_diff=ainv_fit-ainv_old;
	  dboot_t ainv_sum=ainv_fit+ainv_old;
	  
	  cout<<"Significativity ainv beta "<<ibeta<<": "<<
	    ainv_diff.err()/ainv_sum.ave()<<endl;
	}
      
      for(auto& l : {make_tuple(izv,Zv),make_tuple(iza,Za)})
	for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	  {
	    auto i{get<0>(l)};
	    auto Z{get<1>(l)};
	    dboot_t& z_fit=pFit[i[ibeta]];
	  
	    dboot_t z_diff=z_fit-Z[ibeta];
	    dboot_t z_sum=z_fit+Z[ibeta];
	    
	    cout<<"Significativity z beta "<<ibeta<<": "<<
	      z_diff.err()/z_sum.ave()<<endl;
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
