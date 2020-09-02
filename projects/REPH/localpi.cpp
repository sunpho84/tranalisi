#include <tranalisi.hpp>

#include <phys_point.hpp>

constexpr double kappa=2.837297;

map<char,size_t> ibeta_of_id{{'A',0},{'B',1},{'D',2}};

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
  
  /////////////////////////////////////////////////////////////////
  
  boot_init_t bi;
  
  dboot_t aM;
  
  dboot_t daM;
  
  template <typename T>
  struct ToFit
  {
    T ainv;
    
    T M;
    
    T L;
    
    T dM;
    
    T FVE;
    
    T dM2;
    
    T dM2UnivCorrected;
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
  ToFit<T> getToFit(const T& ainv,Oth&&...oth) const
  {
    ToFit<T> out;
    
    out.ainv=ainv;
    out.M=g(aM,oth...)*ainv;
    out.L=g((double)Lfra,oth...)/ainv;
    out.dM=g(e2*daM/2,oth...)*ainv;
    out.FVE=FVEuniv(g(out.L,oth...),g(out.M,oth...));
    out.dM2=out.dM*2*out.M;
    out.dM2UnivCorrected=out.dM2+out.FVE;
    
    return out;
  }
  
  void convertToBoot(const size_t iult)
  {
    bi=jack_index[iult][0];
    aM=dboot_t(bi,aMPiFitted);
    daM=dboot_t(bi,daMPiFitted);
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
	  typename TL>
T FVEansatz(const T& p,const double m2,const TL& L)
{
  return p*sqrt(m2)/(L*L*L);
}

size_t iC,iA1,if0,iD,iFVEnonUniv,ia[nbeta];

template <typename T>
T ansatz(const vector<T>& p,const double& m2,const T& a,const T& L)
{
  const T& C=p[iC];
  const T& A1=p[iA1];
  const T& f0=p[if0];
  const T& D=p[iD];
  
  T Q=4*C/pow(f0,4);
  T W=m2/sqr((T)(4*M_PI*f0));
  
  T uncorrected=e2*sqr(f0)*(Q-(3+4*Q)*W*log(W))+A1*m2+D*sqr(a);
  
  return uncorrected+FVEansatz(p[iFVEnonUniv],m2,L);
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
	(const string& suff)
	{
	  const string corr="mes_contr_"+suff;
	  
	  const djvec_t out=read_djvec(dir+"/jacks/"+corr,T).symmetrized();
	  
	  out.ave_err().write(dir+"/plots/"+corr+".xmg");
	  
	  return out;
	};
      
      const djvec_t P5P5_00=read("00");
      const djvec_t P5P5_LL=read("LL");
      const djvec_t ratio_LL=P5P5_LL/P5P5_00;
      const djvec_t eff_mass=effective_mass(P5P5_00);
      const djvec_t eff_slope=effective_slope(ratio_LL,eff_mass,T/2);
      
      auto fit=
	[&ens,&dir]
	(const djvec_t& vec,const string& corr)
	{
	  return constant_fit(vec,ens.tmin,ens.tmax,dir+"/plots/"+corr+".xmg");
	};
      
      ens.aMPiFitted=fit(eff_mass,"eff_mass");
      ens.daMPiFitted=fit(eff_slope,"eff_slope");
      
      cout<<ens.name<<" "<<smart_print(ens.aMPiFitted)<<" "
	<<smart_print(ens.daMPiFitted)<<endl;
    }
  
  for(size_t iult=0;iult<8;iult++)
    {
      for(auto& ens : ensList)
	ens.convertToBoot(iult);
      
      const size_t nFitPars=8;
      vector<dboot_t> pFit(nFitPars);
      boot_fit_t fit;
      iC=fit.add_fit_par(pFit[0],"C",3.2e-5,0.05);
      iA1=fit.add_fit_par(pFit[1],"A1",-7.1e-4,0.05);
      if0=fit.add_self_fitted_point(pFit[2],"f0",lat_par[iult].f0,-1);
      iD=fit.add_fit_par(pFit[3],"D",2.6e-3,0.05);
      iFVEnonUniv=fit.add_fit_par(pFit[4],"DVEnonUniv",0.94,0.05);
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	ia[ibeta]=fit.add_self_fitted_point(pFit[5+ibeta],"a"+to_string(ibeta),lat_par[iult].ainv[ibeta],-1);
      
      for(auto& ens : ensList)
	fit.add_point([=]
		      (const vector<double>& p,int iboot)
		      {
			const size_t i=ia[ens.ibeta];
			const auto data=ens.getToFit(p[i],iboot);
			
			return data.dM2UnivCorrected;
		      },
		      [=]
		      (const vector<double>& p,int iboot)
		      {
			auto data=ens.getToFit(p[ia[ens.ibeta]],iboot);
			
			const double ainv=data.ainv;
			const double a=1/ainv;
			const double L=data.L;
			const double x=sqr(data.M);
			
			return ansatz(p,x,a,L);
		      },
		      ens.getToFit(lat_par[iult].ainv[ens.ibeta]).dM2UnivCorrected.err());
      
      auto status=fit.fit();
      
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	{
	  dboot_t& afit=pFit[ia[ibeta]];
	  dboot_t& aold=lat_par[iult].ainv[ibeta];
	  
	  dboot_t adiff=afit-aold;
	  
	  cout<<"Significativity beta "<<ibeta<<": "<<adiff.ave_err().significativity()<<endl;
	}
      
      /////////////////////////////////////////////////////////////////
      
      const string ult_plot_dir="plots/"+to_string(iult);
      mkdir(ult_plot_dir);
      grace_file_t daMpi_plot(ult_plot_dir+"/daMpi.xmg");
      
      daMpi_plot.set_settype(grace::XYDY);
      daMpi_plot.set_xaxis_label("M\\s\\x\\p\\0\\N\\S2\\N (GeV\\S-1\\N)");
      daMpi_plot.set_yaxis_label("2\\xd\\0M\\s\\x\\p\\0\\N M\\s\\x\\p\\0\\N (GeV\\S-2\\N)");
      vector<grace::color_t> colors{grace::color_t::BLACK,grace::color_t::RED,grace::color_t::BLUE};
      
      const grace::line_style_t lineStyle[3]={grace::DASHED_LINE,grace::SHORT_DASHED_LINE,grace::CONTINUOUS_LINE};
      for(size_t FVEswitch=0;FVEswitch<3;FVEswitch++)
	for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	  {
	    daMpi_plot.new_data_set();
	    daMpi_plot.set_all_colors(colors[ibeta]);
	    daMpi_plot.set_line_style(lineStyle[FVEswitch]);
	    
	    daMpi_plot.set_transparency(0.33+0.33*FVEswitch);
	    
	    for(auto& ens : ensList)
	      if(ens.ibeta==ibeta)
		{
		  const auto data=ens.getToFit(pFit[ia[ibeta]]);
		  const dboot_t x=sqr(data.M);
		  const vector<dboot_t> y=
		    {data.dM2,
		     data.dM2UnivCorrected,
		     data.dM2UnivCorrected-FVEansatz(pFit[iFVEnonUniv],x.ave(),data.L)};
		  
		  daMpi_plot.write_ave_err(x.ave(),y[FVEswitch].ave_err());
		}
	  }
      
      dboot_t Linf;
      Linf=1e30;
      
      dboot_t a0;
      a0=0.0;
      
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	daMpi_plot.write_polygon([&pFit,&Linf,ibeta,iult]
				 (const double& x) -> dboot_t
				 {
				   const dboot_t a=1/lat_par[iult].ainv[ibeta];
				   
				   return ansatz(pFit,x,a,Linf);
				 },1e-3,0.25,colors[ibeta]);
      
      daMpi_plot.write_polygon([&pFit,&a0,&Linf]
			       (const double& x) -> dboot_t
			       {
				 return ansatz(pFit,x,a0,Linf);
			       },1e-3,0.25,grace::YELLOW);
      
      const double M2PiPhys=sqr(0.139);
      dboot_t dM2PiPhys=ansatz(pFit,M2PiPhys,a0,Linf);
      cout<<"dM2PiPhys: "<<smart_print(dM2PiPhys)<<endl;
      const double dM2PiPDG=1261.2e-6;
      const ave_err_t dM2PiPublished(1137e-6,63e-6);
      
      daMpi_plot.write_ave_err(M2PiPhys,dM2PiPhys.ave_err());
      daMpi_plot.write_ave_err(M2PiPhys,dM2PiPublished);
      daMpi_plot.write_line([&dM2PiPDG](const double x){return dM2PiPDG;},1e-3,0.25,grace::GREEN4,2);
      
      grace_file_t A40slice(ult_plot_dir+"/A40slice.xmg");
      
      for(int FSEflag=0;FSEflag<3;FSEflag++)
	{
	  A40slice.new_data_set();
	  A40slice.set_transparency(0.33+0.33*FSEflag);
	  for(auto& ens: ensList)
	    if(ens.am==0.0040)
	      {
		const auto data=ens.getToFit(pFit[ia[0]]);
		
		djack_t y=e2*ens.daMPiFitted*2*ens.aMPiFitted/2;
		if(FSEflag>0)
		  y+=FVEuniv(ens.Lfra,ens.aMPiFitted);
		if(FSEflag>1)
		  y-=FVEansatz(pFit[3],sqr(data.M.ave()),data.L)/sqr(data.ainv);
	    
	    A40slice.write_ave_err(pow(ens.Lfra,-3),y.ave_err());
	  }
	}
    }
  
  return 0;
}
