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

void single_fat(const string dir="single")
{
  mkdir(dir);
  
  const size_t nFatPars=3;
  vector<djvec_t> pFat(n,djvec_t(nFatPars));
  vector<tuple<djack_t,int,vector<bool>>> status(n);
  vector<bool> skip(n);
  
  grace_file_t ddd(dir+"/ddd.xmg");
  for(size_t i=0;i<n;i++)
    {
      djvec_t& gp=pFat[i];
      cout<<"Fitting "<<i<<endl;
      
      const djvec_t& r=d[i].r;
      
      const djack_t center=r[TH/2];
      
      const djvec_t check=backward_derivative((djvec_t)log(-forward_derivative(r)));//(djvec_t)log(abs(forward_derivative(r))));
      check.ave_err().write(dir+combine("/check_%zu.xmg",i));
      
      const djvec_t o=r-center;
      o.ave_err().write(dir+combine("/off_%zu.xmg",i));
      
      const int L=4;
      int tmin=TH/2;
      int tmax=TH/2+L;
      bool s=false;
      while(tmax>TH/4 and not s)
	{
	  tmax++;
	  s=true;
	  for(int u=tmax-L;u<tmax;u++)
	    {
	      cout<<u<<" "<<o[u].significativity()<<endl;
	      s&=(o[u].significativity()>1 and o[u].ave()<0);
	    }
	}
      
      // y=A*exp(-B*(TH-t));
      // y=A*exp(-B*TH)*exp(-B*t);
      // ln(-y)=ln(-A)-B*TH+B*t;
      // al=ln(-A)-B*TH;
      // A=-exp(al+B*TH)
      
      const djvec_t rsh=poly_fit((djvec_t)log(abs(-o)),1,tmax-L,tmax,dir+combine("/rshould_%zu.xmg",i));
      
      gp[0]=center;
      gp[1]=-exp(rsh[0]+rsh[1]*TH);
      gp[2]=rsh[1];
      
      cout<<gp.ave_err()<<endl;
      cout<<"."<<endl;
      
      bool goOn;
      djack_t ch2;
      size_t ndof;
      size_t iter=0;
      bool worked;
      djvec_t p(3);
      do
	{
	  jack_fit_t fat;
	  for(int j=0;j<3;j++)
	    // if(j!=0)
	      fat.add_fit_par(p[j],combine("p[%d]",j),gp[j].ave_err());
	  //else fat.add_fit_par_limits(p[j],combine("p[%d]",j),p[j].ave_err(),gp[j].ave()-3*gp[j].err(),gp[j].ave()+3*gp[j].err());
	  
	  for(int u=tmin;u<tmax;u++)
	    fat.add_point(r[u],[=](const vector<double>& p,int ijack)
			       {
				 return p[0]+p[1]*exp(-(TH-u)*p[2]);
			       });
	  
	  const auto status=fat.fit();
	  const vector<bool>& ws=get<2>(status);
	  ch2=get<0>(status);
	  ndof=get<1>(status);
	  worked=accumulate(ws.begin(),ws.end(),1,multiplies<double>());
	  
	  goOn=not worked;
	  
	  cout<<"Ch2: "<<ch2.ave_err()<<"/"<<ndof<<endl;
	  // if(worked)
	  //   {
	  
	  bool solutionProposed=false;
	  
	  
	  if(not solutionProposed)
	    {
	      const double p1s=p[1].significativity();
	      const double p2s=p[2].significativity();
	      
	      cout<<"Significativity of par 1: "<<p1s<<", par 2: "<<p2s<<endl;
	      
	      if(p1s<2 or p2s<2)
		{
		  cout<<"Trying to extend right"<<endl;
		  
		  if(tmax<TH)
		    {
		      solutionProposed=true;
		      goOn=true;
		      tmax++;
		    }
		  else
		    cout<<"Cannot extend further right"<<endl;
		}
	    }
	  
	  if(not solutionProposed)
	    {
	      const djack_t lv=p[1]*exp(-(TH-tmin)*p[2]);
	      const double a0=p[0].err(),lva=abs(lv.ave())*3;
	      
	      cout<<"Error of offset: "<<a0<<", leftermost value of the right expo: "<<lva<<endl;
	      
	      if(a0<lva)
		{
		  cout<<a0<<"<"<<lva<<", trying to extend left"<<endl;
		  
		  if(tmin>1)
		    {
		      solutionProposed=true;
		      goOn=true;
		      tmin--;
		    }
		  else
		    cout<<"Cannot extend further left"<<endl;
		}
	    }
	  
	  if(not solutionProposed)
	    {
	      cout<<"No solution propesed, quitting"<<endl;
	      goOn=false;
	    }
	  
	  cout<<p.ave_err()<<endl;
	  cout<<"--"<<endl;
	  
	  iter++;
	}
      while(goOn);
      
      const string outPath=dir+combine("/r_%zu.xmg",i);
      grace_file_t plot(outPath);
      
      plot.write_vec_ave_err(d[i].r.ave_err());
      
      plot.write_polygon([=](const double x) -> djack_t
			 {
			   return p[0]+p[1]*exp(-(TH-x)*p[2]);
			 },tmin,tmax);
      
      cout<<"FINAL "<<i<<" "<<(tmin+tmax)/2.0<<" "<<(tmax-tmin)/2.0<<endl;
      
      plot.write_constant_band(tmin,tmax,p[0]);
      
      ddd.write_ave_err(d[i].X.ave(),p[0].ave_err());
      
      if(not worked)
	CRASH("NOT WORKED");
    }
}

djvec_t single_fit(const vector<djvec_t>& guesses,const string dir="single")
{
  mkdir(dir);
  
  const size_t nFitPars=5;
  vector<djvec_t> pFit(n,djvec_t(nFitPars));
  vector<tuple<djack_t,int,vector<bool>>> status(n);
  vector<bool> skip(n);
  
  int nSkip=0;
  for(size_t i=0;i<n;i++)
    {
      djvec_t& p=pFit[i];
      cout<<"Fitting "<<i<<endl;
      
      int tmin=2;
      int tmax=TH-2;
      
      const djvec_t& r=d[i].r;
      const djvec_t o=symmetric_derivative((djvec_t)log(abs(r-r[TH/2])));
      //o.ave_err().write(dir+combine("/off_%zu.xmg",i));
      
      const djvec_t f=forward_derivative(d[i].r)/d[i].r;
      while(tmin<TH and f[tmin].ave()>0 and f[tmin].significativity()>1)
	tmin++;
      const djvec_t b=backward_derivative(d[i].r)/d[i].r;
      while(tmax>0 and b[tmax].ave()<0 and b[tmax].significativity()>1)
	tmax++;
      
      cout<<" T: "<<tmin<<" "<<tmax<<endl;
      
      // const double min1=-0.05e2;
      jack_fit_t fit;
      for(int j=0;j<4;j++)
	fit.add_fit_par(p[j],combine("p[%d]",j),guesses[i][j].ave(),0.1);
      
      {
	const djack_t center=r[TH/2];
	const double guess=center.ave();
	const double err=center.err();
	const double min=guess-3*err;
	const double max=guess+3*err;
	fit.add_fit_par_limits(p[4],combine("p[%d]",4),guess,err,min,max);
	
	// const djvec_t o=r-center;
	// o.ave_err().write(dir+combine("/off_%zu.xmg",i));
	
	// const int L=4;
	// int t=TH/2+L;
	// bool s=false;
	// while(t>TH/4 and not s)
	//   {
	//     t++;
	//     s=true;
	//     for(int u=t-L;u<t;u++)
	//       {
	// 	cout<<u<<" "<<o[u].significativity()<<endl;
	// 	s&=(o[u].significativity()>1);
	//       }
	//   }
	
	// poly_fit((djvec_t)log(-o),1,t-L,t,dir+combine("/rshould_%zu.xmg",i));
      }
      
      // djack_t e;
      // e.fill_gauss(0,20,23423);
      // fit.add_point(e,[](const vector<double>& p,int ijack){return 1/p[1];});
      // fit.add_point(e,[](const vector<double>& p,int ijack){return 1/p[3];});
      
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
      skip[i]=(ch2.ave()/nDof>3 // or fabs(p[1].ave()-min1)<1e-2
	       );
      if(skip[i])
	nSkip++;
      
      {
	const string outPath=dir+combine("/range_%zu.xmg",i);
	grace_file_t plot(outPath);
	
	djvec_t t=p;
	const djack_t o=t[4];
	t[4]=0.0;
	
	plot.write_line([=](const double x) -> double
			{
			  return abs(single_ansatz(t,x)).significativity();
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
TV get_ansatz_par_all(const TV& p,const T& x,const double& Eg,const double& pi,const double& k)
{
  TV pp(5);
  
  pp[0]=p[0]+p[1]*Eg+p[2]*x+p[3]*pi*pi; //dipende da Eg e k
  pp[1]=p[4]+p[5]*Eg+p[6]*x+p[7]*pi; // dipende da k e pi
  pp[2]=(p[8]+p[9]*Eg+p[10]*x+p[11]*pi*pi)/x; // non dipende da nulla
  pp[3]=p[12]+p[13]*Eg+p[14]*x+p[15]*pi*pi;
  pp[4]=p[16]+p[17]*Eg+p[18]*x+p[19]*pi*pi;
  
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

djvec_t all_fit(const djvec_t& in,const string dir="all")
{
  int tmin=1;
  int tmax=TH-2;
  
  bool ok;
  const size_t nFitPars=5*4;
  djvec_t pFit(nFitPars);
  do
    {
      cout<<"Tmin: "<<tmin<<", Tmax: "<<tmax<<endl;
      if(tmax-tmin<TH/2)
	CRASH("tmax-tmin too close");
      
      jack_fit_t fit;
      //Z1
      fit.add_fit_par(pFit[0],"p[0]",in[0].ave_err());
      fit.add_fit_par(pFit[1],"p[1]",0.2,0.1);//fit.fix_par(1);
      fit.add_fit_par(pFit[2],"p[2]",0,0.1);fit.fix_par(2);
      fit.add_fit_par(pFit[3],"p[3]",0,0.1);fit.fix_par(3);
      //M1
      fit.add_fit_par_limits(pFit[4],"p[4]",in[2].ave_err(),0.1,4);
      fit.add_fit_par(pFit[5],"p[5]",0.0,0.1);fit.fix_par(5);
      fit.add_fit_par(pFit[6],"p[6]",0,0.1);//fit.fix_par(6);
      fit.add_fit_par(pFit[7],"p[7]",0,0.1);fit.fix_par(7);
      //Z2
      fit.add_fit_par(pFit[8],"p[8]",in[4].ave_err());
      fit.add_fit_par(pFit[9],"p[9]",0,0.1);fit.fix_par(9);
      fit.add_fit_par(pFit[10],"p[10]",0,0.1);//fit.fix_par(10);
      fit.add_fit_par(pFit[11],"p[11]",0,0.1);fit.fix_par(11);
      //M2
      fit.add_fit_par_limits(pFit[12],"p[12]",in[6].ave_err(),0.1,4);
      fit.add_fit_par(pFit[13],"p[13]",-1,0.1);fit.fix_par(13);
      fit.add_fit_par(pFit[14],"p[14]",0,0.1);//fit.fix_par(14);
      fit.add_fit_par(pFit[15],"p[15]",0,0.1);fit.fix_par(15);
      //C
      fit.add_fit_par(pFit[16],"p[16]",in[8].ave_err());
      fit.add_fit_par(pFit[17],"p[17]",0,0.1);fit.fix_par(17);
      fit.add_fit_par(pFit[18],"p[18]",0,0.1);//fit.fix_par(18);
      fit.add_fit_par(pFit[19],"p[19]",0,0.1);fit.fix_par(19);
      
      // for(size_t i=0;i<nFitPars;i++)
      //   fit.fix_par(i);
      
      for(size_t i=0;i<n;i++)
	{
	  for(int t=tmin;t<tmax;t++)
	    fit.add_point(d[i].r[t],[=](const vector<double>& p,int ijack)->double
				    {
				      return
					all_ansatz(p,d[i].X[ijack],d[i].Eg,sqr(d[i].p),d[i].k,t);
				    });
	};
      
      const auto status=fit.fit();
      
      ok=true;
      for(const bool& s : get<2>(status))
	ok&=s;
      
      if(not ok)
	{
	  tmin++;
	  tmax--;
	}
    }
  while(not ok);
  
  mkdir(dir);
  
  djvec_t const_fix_range(n);
  grace_file_t plot_const_fix_range(dir+"/const_fix_range.xmg");
  
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
      
      {
	djvec_t o=pFit;
	for(int i=16;i<20;i++)
	  o[i]=0;
	
	vector<double> c(TH+1);
	
	ofstream signif(dir+combine("/sign_%zu.xmg",i));
	for(int x=0;x<=TH;x++)
	  {
	    c[x]=abs(all_ansatz(o,d[i].X,d[i].Eg,sqr(d[i].p),d[i].k,x)).significativity();
	    signif<<x<<" "<<c[x]<<endl;
	  }
	
	auto it=min_element(c.begin(),c.end());
	int l=distance(c.begin(),it),r=l;
	for(int y=0;y<4;y++)
	  {
	    const int lp=(l-1+TH+1)%(TH+1);
	    const int rp=(r+1+TH+1)%(TH+1);
	    if(lp<l and c[lp]<c[rp])
	      l=lp;
	    else
	      r=rp;
	  }
	
	cout<<l<<" "<<r<<endl;
	
	const_fix_range[i]=constant_fit(d[i].r,l,r,dir+combine("/const_fit_%zu.xmg",i));
	
	plot_const_fix_range.write_ave_err(d[i].X.ave(),const_fix_range[i].ave_err());
      }
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
  
  single_fat();
  const djvec_t guess=single_fit(guesses);
  const djvec_t pFit=all_fit(guess);
  
  for(size_t i=0;i<n;i++)
    guesses[i]=get_ansatz_par_all(pFit,d[i].X,d[i].Eg,sqr(d[i].p),d[i].k);
  
  single_fit(guesses,"single_again");
  
  return 0;
}
