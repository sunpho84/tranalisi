#include <tranalisi.hpp>

const size_t T=48;

int main()
{
  if(1){
    raw_file_t in("tm_tuning","r");
  
  char line[1024];
  
  std::map<double,vector<std::array<double,6>>> raw;
  
  bool goon=true;
  do
    {
      in.get_line(line);
      istringstream is(line);
      string dum;
      
      size_t conf,flav;
      double mass,kappa;
      is>>
  	dum>>dum>>conf>>
  	dum>>dum>>dum>>flav>>
  	dum>>dum>>dum>>mass>>
  	dum>>dum>>dum>>kappa;
      
      goon=(not in.feof());
      
      dcomplex l[6];
      
      if(goon)
  	{
  	  for(size_t t=0;t<T;t++)
  	    {
  	      in.expect(to_string(t).c_str());
	      
  	      for(int i=0;i<6;i++)
  		{
  		  const double re=in.read<double>();
  		  in.expect(",");
  		  const double im=in.read<double>();
		  
  		  l[i]={re,im};
  		}
	      
  	      const double P5P5=l[0].real();
  	      const double V0P5=l[1].imag();
  	      const double P5P5_ins_S=2*l[2].real();
  	      const double V0P5_ins_S=2*l[3].imag();
  	      const double P5P5_ins_P=2*l[4].imag();
  	      const double V0P5_ins_P=-2*l[5].real();
	      
  	      raw[kappa].push_back({P5P5,V0P5,P5P5_ins_S,V0P5_ins_S,P5P5_ins_P,V0P5_ins_P});
  	    }
	  
  	  in.get_line(line);
  	  in.get_line(line);
  	}
    }
  while(goon);
  
  grace_file_t P5P5_plot("P5P5.xmg");
  grace_file_t mcr_plot("mcr.xmg");
  grace_file_t dmcr_plot("dmcr.xmg");
  grace_file_t V0P5_plot("V0P5.xmg");
  grace_file_t P5P5_der_P_plot("P5P5_der_P.xmg");
  grace_file_t V0P5_der_P_plot("V0P5_der_P.xmg");
  
  grace_file_t mcr_fun_kappa_plot("mcr_fun_kappa.xmg");
  grace_file_t mpi_fun_kappa_plot("mpi_fun_kappa.xmg");
  
  for(auto& s : raw)
    {
      const double kappa=get<0>(s);
      const auto d=get<1>(s);
      const size_t n=d.size()/T;
      const size_t dump=0;
      njacks=n-dump;
      
      if(njacks>1)
  	{
  	  cout<<njacks<<" "<<kappa<<endl;
	  
  	  djvec_t P5P5(T);
  	  djvec_t V0P5(T);
  	  djvec_t P5P5_ins_P(T);
  	  djvec_t V0P5_ins_P(T);
	  
  	  for(size_t iconf=dump;iconf<n;iconf++)
  	    for(size_t t=0;t<T;t++)
  	      {
  		const auto& a=d[t+T*iconf];
		
  		P5P5[t][iconf-dump]=a[0];
  		V0P5[t][iconf-dump]=a[1];
  		P5P5_ins_P[t][iconf-dump]=a[4];
  		V0P5_ins_P[t][iconf-dump]=a[5];
  	      }
	  
  	  P5P5.clusterize();
  	  V0P5.clusterize();
  	  P5P5_ins_P.clusterize();
  	  V0P5_ins_P.clusterize();
	  
  	  P5P5.symmetrize();
  	  V0P5.symmetrize(-1);
	  P5P5_ins_P.symmetrize();
  	  V0P5_ins_P.symmetrize(-1);
	  
  	  const djvec_t eff_pi=effective_mass(P5P5);
  	  P5P5_plot.write_vec_ave_err(P5P5.ave_err());
  	  const djvec_t mcr_corr=forward_derivative(V0P5)/P5P5;
	  const djack_t mcr=mcr_corr[T/4];
	  
  	  mcr_plot.write_vec_ave_err(mcr_corr.ave_err());
	  
  	  mcr_fun_kappa_plot.write_ave_err(kappa,mcr.ave_err());
  	  mpi_fun_kappa_plot.write_ave_err(kappa,eff_pi[T*3/8].ave_err());
	  
  	  V0P5_plot.write_vec_ave_err(V0P5.ave_err());
	  
	  const djvec_t P5P5_der_P=P5P5_ins_P/P5P5;
  	  P5P5_der_P_plot.write_vec_ave_err(P5P5_der_P.ave_err());
	  
	  const djvec_t V0P5_der_P=V0P5_ins_P/V0P5;
  	  V0P5_der_P_plot.write_vec_ave_err(V0P5_der_P.ave_err());
	  
	  const djvec_t dmcr_corr=forward_derivative(V0P5_ins_P)/P5P5-forward_derivative(V0P5)*P5P5_ins_P/sqr(P5P5);
	  const djack_t dmcr=dmcr_corr[T/4];
  	  dmcr_plot.write_vec_ave_err(dmcr_corr.ave_err());
	  
  	  for(auto& plot : {&P5P5_plot,&mcr_plot,&dmcr_plot,&V0P5_plot,&P5P5_der_P_plot,&V0P5_der_P_plot})
  	    {
  	      plot->set_settype(grace::XYDY);
  	      plot->set_legend(to_string(kappa));
  	    }
	  
	  const double m0=0.5/kappa;
	  const djack_t dm0=-mcr/dmcr;
	  const djack_t kappa_prime=0.5/(m0+dm0);
	  cout<<"dm0: "<<dm0.ave_err()<<endl;
	  cout<<"kappa_prime: "<<kappa_prime.ave_err()<<endl;
	}
    }
}
  if(0)
  {
    raw_file_t input("analysis.txt","r");
  njacks=input.read<size_t>("NJacks");
  const double kappa=input.read<double>("Kappa");
  
  const djvec_t V0P5=read_djvec("00_V0P5",48,1).symmetrized(0);
  const djvec_t P5P5=read_djvec("00_P5P5",48,0).symmetrized(0);
  const djvec_t V0P5_ins_P=-2*read_djvec("0P_V0P5",48,0).symmetrized(0);
  const djvec_t P5P5_ins_P=2*read_djvec("0P_P5P5",48,1).symmetrized(0);
  const djvec_t mcr_corr=forward_derivative(V0P5)/P5P5;
  const djvec_t dmcr_corr=forward_derivative(V0P5_ins_P)/P5P5-forward_derivative(V0P5)*P5P5_ins_P/sqr(P5P5);
  mcr_corr.ave_err().write("mcr.xmg");
  dmcr_corr.ave_err().write("dmcr.xmg");
  
  const djack_t mcr=mcr_corr[T/4];
  cout<<"mcr: "<<1/(2*kappa)<<" "<<mcr.ave_err()<<endl;
  const djack_t dmcr=dmcr_corr[T/4];
  cout<<"dmcr: "<<dmcr.ave_err()<<endl;
  
  const double m0=0.5/kappa;
  const djack_t dm0=-mcr/dmcr;
  const double kappa_prime=0.5/(m0+dm0.ave());
  cout<<"dm0: "<<dm0.ave_err()<<endl;
  cout<<"kappa_prime: "<<kappa_prime<<endl;
  
  // const djack_t =1/dmcr[T/4];
  // cout<<der.ave_err()<<endl;
  
  // const double m0=1/(2*kappa)-4;
  // const double new_kappa=1/(2*(m0+4+der.ave()));
  // cout<<new_kappa<<endl;
  }
  return 0;
}
