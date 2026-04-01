#include <tranalisi.hpp>

const size_t T=64;

int main()
{
  set_njacks(32);
  djvec_t m[3][3];

  for(size_t i1=0;i1<3;i1++)
    for(size_t i2=0;i2<3;i2++)
      m[i1][i2]=read_djvec("data/mes_contr_2PT_P5P5",T,i1+3*i2).symmetrized();
  
  for(size_t i1=0;i1<3;i1++)
    for(size_t i2=i1+1;i2<3;i2++)
      m[i1][i2]=m[i2][i1]=(m[i1][i2]+m[i2][i1])/2;
  
  const size_t tStar=6;
  djack_t th,ph;
  for(size_t iJack=0;iJack<=njacks;iJack++)
    {
      fun_minuit_wrapper_t f([&m,&iJack](const vector<double>& p)
      {
	const double th=p[0];
	const double ph=p[1];
	double c[2]{};
	const double f[3]={sin(th)*cos(ph),sin(th)*sin(ph),cos(th)};
	for(size_t i1=0;i1<3;i1++)
	  for(size_t i2=0;i2<3;i2++)
	    for(size_t i=0;i<2;i++)
	      c[i]+=m[i1][i2][tStar+i][iJack]*f[i1]*f[i2];
	
	const double m=effective_mass(c[0],c[1],tStar,T/2);
	// if(iJack==0) cout<<th<<" "<<ph<<" "<<m<<endl;
	
	return m;
      });
      
      minimizer_pars_t pars;
      pars.add("th",0.5,0.001);
      pars.add("ph",0.1,0.001);
      pars.setlimits("th",0,M_PI);
      pars.setlimits("ph",-M_PI,M_PI);
      minimizer_t min(f,pars);
      
      const vector<double> p=min.minimize();
      th[iJack]=p[0];
      ph[iJack]=p[1];
    }
  
  djvec_t best(T/2+1);
  const djack_t f[3]={sin(th)*cos(ph),sin(th)*sin(ph),cos(th)};
  for(size_t i1=0;i1<3;i1++)
    for(size_t i2=0;i2<3;i2++)
      best+=m[i1][i2]*f[i1]*f[i2];
  
  cout<<th.ave_err()<<endl;
  cout<<ph.ave_err()<<endl;
  cout<<"Weights:"<<endl;
  for(size_t i=0;i<3;i++)
    cout<<f[i].ave_err()<<endl;
  
  grace_file_t plot("plots/best.xmg");
  plot.write_vec_ave_err(effective_mass(best).ave_err());
  for(size_t i=0;i<3;i++)
    plot.write_vec_ave_err(effective_mass(m[i][i]).ave_err());
  
  vector<djvec_t> corr(4);
  const int r[2]={1,2};
  for(size_t i1=0;i1<2;i1++)
    for(size_t i2=0;i2<2;i2++)
      corr[i2+2*i1]=m[r[i1]][r[i2]];
  
  auto [eig,recastEigvec,ignore]=gevp(corr,2);
  cout<<"gevp: "<<endl;
  for(size_t i1=0;i1<2;i1++)
    for(size_t i2=0;i2<2;i2++)
      cout<<recastEigvec[i1+2*i2][tStar].ave_err()<<endl;
  effective_mass(eig[0]).ave_err().write("plots/gevp.xmg");
  effective_mass(eig[1]).ave_err().write("plots/gevpPrime.xmg");
  //effective_mass(eig[2]).ave_err().write("plots/gevpSecond.xmg");
  
  // const int n[3]={18,36,72};
  // grace_file_t gaussPlot("plots/gauss.xmg");
  // for(size_t i=0;i<3;i++)
  //   gaussPlot.write_line([i,&n](const double& x)
  //   {
  //     return sqrt(-log(gauss_distr(x,0.0,sqrt(0.4*n[i]))/gauss_distr(0,0.0,sqrt(0.4*n[i]))))/x;
  //   },0.1,20);
  
  // gaussPlot.write_line([&f,&n](const double& x)
  // {
  //   double y=0,z=0;
  //   for(size_t i=0;i<3;i++)
  //     {
  // 	y+=gauss_distr(x,0.0,sqrt(0.4*n[i]))*f[i].ave();
  // 	z+=gauss_distr(0.0,0.0,sqrt(0.4*n[i]))*f[i].ave();
  //     }
    
  //   return sqrt(-log(y/z))/x;
  // },0.1,20);
  
  // grace_file_t excPlot("plots/exc.xmg");
  // for(int ist=0;ist<3;ist++)
  //   {
  //     djack_t b[3];
  //     for(size_t i2=0;i2<3;i2++)
  // 	b[i2]=recastEigvec[ist+3*i2][tStar];
      
  //     excPlot.write_line([&b,&n](const double& x)
  //     {
  // 	double y=0,z=0;
  // 	for(size_t i=0;i<3;i++)
  // 	  {
  // 	    y+=gauss_distr(x,0.0,sqrt(0.4*n[i]))*b[i].ave();
  // 	    z+=gauss_distr(0.0,0.0,sqrt(0.4*n[i]))*b[i].ave();
  // 	  }
	
  // 	return y/z;
  //     },0.1,20);
  //   }
  
  return 0;
}
