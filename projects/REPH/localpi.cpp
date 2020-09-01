#include <tranalisi.hpp>

#include <phys_point.hpp>

map<char,size_t> ibeta_of_id{{'A',0},{'B',1},{'D',2}};

const size_t parseL(const string& name)
{
  string temp;
  
  bool passed=false;
  
  for(auto& t : name)
    {
      if(passed) temp+=t;
      if(t=='.') passed=true;
    }
  
  cout<<temp<<endl;
  
  return atoi(temp.c_str());
}

struct perens_t
{
  const string name;
  
  const size_t T;
  
  const size_t L;
  
  const size_t tmin;
  
  const size_t tmax;
  
  const size_t ibeta;
  
  djack_t aMPi;
  
  djack_t daMPi;
  
  perens_t(const string& name,const size_t& T,const size_t& tmin,const size_t& tmax) :
    name(name),T(T),L(parseL(name)),tmin(tmin),tmax(tmax),ibeta(ibeta_of_id[name[0]])
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
      
      ens.aMPi=fit(eff_mass,"eff_mass");
      ens.daMPi=fit(eff_slope,"eff_slope");
      
      cout<<smart_print(ens.daMPi)<<endl;
    }
  
  mkdir("plots");
  grace_file_t daMpi_plot("plots/daMpi.xmg");
  daMpi_plot.set_settype(grace::XYDY);
  vector<grace::color_t> colors{grace::color_t::RED,grace::color_t::ORANGE,grace::color_t::GREEN4};
  
  for(size_t without_with_FVE=0;without_with_FVE<2;without_with_FVE++)
    for(size_t ibeta=0;ibeta<nbeta;ibeta++)
      {
	const double ainv=lat_par[0].ainv[ibeta].ave();
	
	daMpi_plot.new_data_set();
	daMpi_plot.set_all_colors(colors[ibeta]);
	daMpi_plot.set_line_style(without_with_FVE?grace::CONTINUOUS_LINE:grace::DASHED_LINE);
	
	daMpi_plot.set_transparency(without_with_FVE?1:0.5);
	
	for(auto& ens : ensList)
	  if(ens.ibeta==ibeta)
	    {
	      const djack_t aM=ens.aMPi;
	      const djack_t daM=ens.daMPi;
	      const djack_t M=aM*ainv;
	      const djack_t dM=daM*ainv;
	      const double kappa=2.837297;
	      const djack_t FVE=kappa/sqr(ens.L)*(1+2*M);
	      
	      const djack_t x=sqr(M);
	      const djack_t y=dM*2*M+FVE*without_with_FVE;
	      
	      daMpi_plot.write_ave_err(x.ave(),y.ave_err());
	    }
    }
  
  return 0;
}
