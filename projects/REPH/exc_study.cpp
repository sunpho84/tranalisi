#include <tranalisi.hpp>

int TH;

struct data_t
{
  djvec_t r;
  djack_t X;
  double p;
  double Eg;
  double k;
  size_t c[3];
  
  data_t() : r(TH+1) {}
};

vector<data_t> d;
size_t n;

void load()
{
  set_njacks(15);
  
  const string path="exc_study.dat";
  
  raw_file_t input(path,"r");
  
  n=input.bin_read<size_t>();
  TH=input.bin_read<int>();
  
  cout<<TH<<endl;
  
  d.resize(n);
  
  ofstream out("kinematics.txt");
  
  for(size_t i=0;i<n;i++)
    {
      d[i].r.bin_read(input);
      d[i].X=input.bin_read<djack_t>();
      d[i].p=input.bin_read<double>();
      d[i].Eg=input.bin_read<double>();
      d[i].k=input.bin_read<double>();
      for(int j=0;j<3;j++) input.bin_read(d[i].c[j]);
      
      out<<d[i].X.ave()<<" "<<d[i].p<<" "<<d[i].Eg<<" "<<d[i].k<<"    "<<i<<"   "<<d[i].c[0]<<" "<<d[i].c[1]<<" "<<d[i].c[2]<<endl;
    };
}

template <typename TV> auto single_ansatz(const TV& p,const double t) -> std::remove_reference_t<decltype(p[0])>
{
  return
    p[0]*exp(-p[1]*t)+
    p[2]*exp(-p[3]*(TH-t))+
    p[4];
}

djvec_t single_fit(const vector<djvec_t>& guesses,const string dir="single")
{
  mkdir(dir);
  
  const size_t nFitPars=5;
  vector<djvec_t> pFit(n,djvec_t(nFitPars));
  vector<tuple<djack_t,int>> status(n);
  vector<bool> skip(n);
  
  int nSkip=0;
  for(size_t i=0;i<n;i++)
    {
      djvec_t& p=pFit[i];
      cout<<"Fitting "<<i<<endl;
      
      // const double min1=-0.05e2;
      jack_fit_t fit;
      for(int j=0;j<5;j++)
	fit.add_fit_par(p[j],combine("p[%d]",j),guesses[i][j].ave(),0.1);
      
      const int tmin=1;
      const int tmax=TH-2;
      
      djack_t e;
      e.fill_gauss(0,20,23423);
      fit.add_point(e,[](const vector<double>& p,int ijack){return 1/p[1];});
      fit.add_point(e,[](const vector<double>& p,int ijack){return 1/p[3];});
      
      for(int t=tmin;t<tmax;t++)
	fit.add_point(d[i].r[t],[=](const vector<double>& p,int ijack)
				{
				  return single_ansatz(p,t);
				});
      
      status[i]=fit.fit();
      
      const string outPath=dir+combine("/%zu.xmg",i);
      grace_file_t plot(outPath);
      
      plot.write_vec_ave_err(d[i].r.ave_err());
      
      plot.write_polygon([=](const double x) -> djack_t
			 {
			   return single_ansatz(p,x);
			 },tmin,tmax);
      
      plot.write_constant_band(tmin,tmax,p[4]);
      
      const djack_t& ch2=get<0>(status[i]);
      const size_t& nDof=get<1>(status[i]);
      skip[i]=(ch2.ave()/nDof>1.5 // or fabs(p[1].ave()-min1)<1e-2
	       );
      if(skip[i])
	nSkip++;
      
      {
      const string outPath=dir+combine("/range_%zu.xmg",i);
      grace_file_t plot(outPath);
      
      djvec_t t=p;
      const djack_t o=t[4];
      t[4]=0.0;
      
      plot.write_polygon([=](const double x) -> djack_t
			 {
			   return abs(single_ansatz(t,x)/o);
			 },tmin,tmax);
	  }
    }
  
  djvec_t out(10);
  for(size_t iPar=0;iPar<nFitPars;iPar++)
    {
      vector<double> x(n-nSkip);
      djvec_t y(n-nSkip);
      
      for(int what=0;what<4;what++)
	{
	  string whatTag[]={"x","Eg","p","k"};
	  
	  const string outPath=dir+combine("/par_%zu_fun_%s.xmg",iPar,whatTag[what].c_str());
	  
	  int j=0;
	  for(size_t i=0;i<n;i++)
	    if(not skip[i])
	      {
		y[j]=pFit[i][iPar];
		
		switch(what)
		  {
		  case 0:
		    x[j]=d[i].X.ave();
		    if(iPar==2)
		      y[j]*=d[i].X;
		    break;
		  case 1:
		    x[j]=d[i].Eg;
		    break;
		  case 2:
		    x[j]=sqr(d[i].p);
		    break;
		  case 3:
		    x[j]=d[i].k;
		    break;
		  }
		
		j++;
	      }
	  
	  const djvec_t p=poly_fit(x,y,1,0.0,2.0,outPath);
	  cout<<p<<endl;
	  
	  if(what==0)
	    for(int k=0;k<2;k++)
	      out[2*iPar+k]=p[k];
	}
    }
  
  return out;
}

template <typename TV,
	  typename T=std::remove_reference_t<decltype((*(TV*)nullptr)[0])>>
TV get_ansatz_par_all(const TV& p,const T& x,const double& Eg,const double& p2,const double& k)
{
  TV pp(5);
  
  pp[0]=p[0]+p[1]*Eg+p[2]*k+p[3]*p2;
  pp[1]=p[4]+p[5]*Eg+p[6]*k+p[7]*p2;
  pp[2]=(p[8]+p[9]*Eg+p[10]*k+p[11]*p2)/x;
  pp[3]=p[12]+p[13]*Eg+p[14]*k+p[15]*p2;
  pp[4]=p[16]+p[17]*Eg+p[18]*k+p[19]*p2;
  
  return pp;
}

template <typename TV,
	  typename T=std::remove_reference_t<decltype((*(TV*)nullptr)[0])>>
T all_ansatz(const TV& p,const T& x,const double& Eg,const double& p2,const double& k,const double t)
{
  TV pp=get_ansatz_par_all(p,x,Eg,p2,k);
  
  return
    single_ansatz(pp,t);
}

djvec_t all_fit(const djvec_t& in)
{
  const string dir="all";
  
  const size_t nFitPars=5*4;
  djvec_t pFit(nFitPars);
  jack_fit_t fit;
  //Z1
  fit.add_fit_par(pFit[0],"p[0]",in[0].ave_err());
  fit.add_fit_par(pFit[1],"p[1]",0,0.1);
  fit.add_fit_par(pFit[2],"p[2]",0,0.1);
  fit.add_fit_par(pFit[3],"p[3]",0,0.1);
  //M1
  fit.add_fit_par_limits(pFit[4],"p[4]",in[2].ave_err(),0.1,4);
  fit.add_fit_par(pFit[5],"p[5]",0,0.1);
  fit.add_fit_par(pFit[6],"p[6]",0,0.1);
  fit.add_fit_par(pFit[7],"p[7]",0,0.1);
  //Z2
  fit.add_fit_par(pFit[8],"p[8]",in[4].ave_err());
  fit.add_fit_par(pFit[9],"p[9]",0,0.1);
  fit.add_fit_par(pFit[10],"p[10]",0,0.1);
  fit.add_fit_par(pFit[11],"p[11]",0,0.1);
  //M2
  fit.add_fit_par_limits(pFit[12],"p[12]",in[6].ave_err(),0.1,4);
  fit.add_fit_par(pFit[13],"p[13]",0,0.1);
  fit.add_fit_par(pFit[14],"p[14]",0,0.1);
  fit.add_fit_par(pFit[15],"p[15]",0,0.1);
  //C
  fit.add_fit_par(pFit[16],"p[16]",in[8].ave_err());
  fit.add_fit_par(pFit[17],"p[17]",0,0.1);
  fit.add_fit_par(pFit[18],"p[18]",0,0.1);
  fit.add_fit_par(pFit[19],"p[19]",0,0.1);
  
  // for(size_t i=0;i<nFitPars;i++)
  //   fit.fix_par(i);
  
  const int tmin=1;
  const int tmax=TH-1;
  
  for(size_t i=0;i<n;i++)
    {
      for(int t=tmin;t<tmax;t++)
	fit.add_point(d[i].r[t],[=](const vector<double>& p,int ijack)->double
			   {
			     return
			       all_ansatz(p,d[i].X[ijack],d[i].Eg,sqr(d[i].p),d[i].k,t);
			   });
    };
  
  fit.fit();
  
  mkdir(dir);
  
  for(size_t i=0;i<n;i++)
    {
      const string outPath=dir+combine("/%zu.xmg",i);
      grace_file_t plot(outPath);
      plot.set_title(combine("x: %lg, p: %lg, k: %lg",d[i].X.ave(),d[i].p,d[i].k));
      
      plot.write_vec_ave_err(d[i].r.ave_err());
      
      plot.write_polygon([=](const double t) -> djack_t
			 {
			   return all_ansatz(pFit,d[i].X,d[i].Eg,sqr(d[i].p),d[i].k,t);
			 },tmin,tmax);
      
      // const string outRange=dir+combine("/range_%zu.xmg",i);
      // grace_file_t plotRange(outRange);
      
      // djvec_t o=pFit;
      // o[8]=o[9]=0;
      // plotRange.write_line([=](const double t) -> double
      // 			 {
      // 			   return ((djack_t)log(sqr(all_ansatz(o,d[i].X,d[i].Eg,sqr(d[i].p),t)))).ave();
      // 			 },tmin,tmax);
      
    }
  
  cout<<pFit.ave_err()<<endl;
  
  return pFit;
}

int main(int narg,char **arg)
{
  
  load();
  
  vector<djvec_t> guesses(n,djvec_t(5));
  for(size_t i=0;i<n;i++)
    {
      guesses[i][0]=1.26e-1;
      guesses[i][1]=7e-1;
      guesses[i][2]=-9e-1;
      guesses[i][3]=3e-1;
      guesses[i][4]=0.0;
    }
  
  const djvec_t guess=single_fit(guesses);
  const djvec_t pFit=all_fit(guess);
  
  for(size_t i=0;i<n;i++)
    guesses[i]=get_ansatz_par_all(pFit,d[i].X,d[i].Eg,sqr(d[i].p),d[i].k);
  
  single_fit(guesses,"single_again");
  
  return 0;
}
