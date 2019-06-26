#ifndef _PERMES_COMBO_HPP
#define _PERMES_COMBO_HPP

#include <REPH/perens.hpp>

#include <set>

//! Hold all data to draw ff
template <typename TV=djvec_t>
struct permes_t
{
  //! Reference ensemble
  const perens_t& ens;
  
  //! Combination relative to the meson combination
  const std::string mesTag;
  
  //! Build milled path
  string milledPath() const
  {
    const string path=ens.dirPath+"/milledData/"+mesTag;
    static set<string> createdDir;
    
    if(createdDir.find(path)==createdDir.end())
      {
	mkdir(path);
	createdDir.insert(path);
      }
    
    return path;
  }
  
  //! Energy
  TV E;
  
  //! Form factor independent variable
  TV X;
  
  //! Default counstructor
  permes_t(const perens_t& ens,const string& mesTag) : ens(ens),mesTag(mesTag){}
  
  template <typename F,typename S>
  auto interpolate(const vector<double>& m1,const vector<double>& m2,const F& f,const vector<S>& data);
  
  //! 2 masses interpolation
  template <typename F,typename S>
  auto interpolate2D(const vector<double>& m1,const vector<double>& m2,const F& f,const vector<S>& data);
  
  //! 1 mass interpolation
  template <typename F,typename S>
  auto interpolate1D(const vector<double>& m,const F& f,const vector<S>& data);
  
  //! 0 mass interpolation
  template <typename F,typename S>
  auto interpolate0D(const F& f,const vector<S>& data);
  
  //! Interpolate from a list
  template <typename S>
  void interpolate(const vector<double>& m1,const vector<double>& m2,const vector<S>& data)
  {
#define INTERP(A) \
    A=interpolate(m1,m2,[](const auto& D){return D.A;},data)
    
    INTERP(E);
    INTERP(X);
    INTERP(ff[0]);
    INTERP(ff[1]);
    
#undef INTERP
  }
  
  //! Form factors for V and A
  array<djvec_t,2> ff;
  
  //! Plot ff
  void plotFf();
};

//! Holds all info for a combination of quarks
template <typename TV=djvec_t>
struct permes_combo_t : public permes_t<TV>
{
  using permes_t<TV>::ens;
  using permes_t<TV>::ff;
  using permes_t<TV>::E;
  using permes_t<TV>::X;
  using permes_t<TV>::milledPath;
  using permes_t<TV>::mesTag;
  
  //! Index of the spectator quark
  const size_t iMs;
  
  //! Index of the forward line quark
  const size_t iMt;
  
  //! Charge of the spectator quark
  const double eS;
  
  //! Charge of the forward line quark
  const double eT;
  
  //! Pseudoscalar coupling
  djvec_t ZP;
  
  //! Axial coupling
  djvec_t ZA;
  
  //! Decay constant taken from Pseudoscalar current
  djvec_t fP;
  
  //! Decay constant taken from Axial current
  djvec_t fPbare;
  
  //! Effective mass
  vector<djvec_t> eEff;
  
  //! Difference of energies between initial and final state
  djvec_t dEdec;
  
  //! Quadrimomentum product
  djvec_t PKdec;
  
  //! Normalization to be used for 3pts
  vector<djvec_t> normaliz;
  
  //! Pseudoscalar correlation function
  vector<djvec_t> corrPP;
  
  //! T component of axial current
  vector<djvec_t> corrA0P;
  
  //! Z component of axial current
  vector<djvec_t> corrA3P;
  
  enum{insOnS,insOnT};
  
  //! Decay correlators for V and A, with the insertion on S or T
  array<array<vector<djvec_t>,2>,2> corrPX;
  
  //! Time interval for 3pts fit
  array<array<vector<Range>,2>,2> tint3pts;
  
  //! Time interval for 2pts fit
  vector<Range> tint2pts;
  
  //! Load the PP correlation function
  djvec_t load2ptsPP(const size_t iMoms,const size_t iMomt);
  
  //! Load the AP correlation function
  djvec_t load2ptsAP(const size_t iMoms,const size_t iMomt,const size_t iGamma);
  
  //! Load the three points correlation functions
  djvec_t load3pts(const size_t iVA,const size_t iMs,const size_t iMt,const size_t iMoms,const size_t iMomt,const size_t iMom0);
  
  //! Computes the axial couplings
  void computeAxialPseudoCouplings(const string& mesPlotsPath);
  
  //! Plot the dispersion relation
  void plotDispRel(const string& mesPlotsPath) const;
  
  //! Load all 2pts
  void load2pts(const bool forceLoad);
  
  //! Load all 3pts
  void load3pts(const bool forceLoad);
  
  //! Load all data
  void load(const bool forceLoad=false)
  {
    
    load2pts(forceLoad);
    load3pts(forceLoad);
  }
  
  //! Plots the form factors
  permes_combo_t& plotFf()
  {
    static_cast<permes_t<TV>*>(this)->plotFf();
    
    return *this;
  }
  
  //! Constructor
  permes_combo_t(const perens_t& ens,const string& mesName,const size_t& iMs,const size_t& iMt,const double& eS,const double& eT) :
    permes_t<TV>(ens,combine("mes_%s/iMs%zu" "_" "iMt%zu",mesName.c_str(),iMs,iMt)),
    iMs(iMs),
    iMt(iMt),
    eS(eS),
    eT(eT)
  {
    resizeListOfContainers({&ZP,&ZA,&fP,&fPbare,&this->E},ens.nMesKin);
    resizeListOfContainers({&eEff},ens.nMesKin);
    resizeListOfContainers({&tint2pts},ens.nMesKin);
    
    resizeListOfContainers({&corrA0P,&corrA3P,&corrPP},ens.nMesKin,djvec_t{(size_t)ens.T/2+1});
    
    resizeListOfContainers({&this->ff[0],&this->ff[1]},ens.nDecKin);
    
    for(size_t iVA=0;iVA<2;iVA++)
      for(size_t iST=0;iST<2;iST++)
	{
	  resizeListOfContainers({&corrPX[iVA][iST]},ens.nDecKin,djvec_t{(size_t)ens.T/2+1});
	  resizeListOfContainers({&tint3pts[iVA][iST]},ens.nDecKin,Range{0,0});
	}
    
    resizeListOfContainers({&dEdec,&PKdec,&this->X},ens.nDecKin);
    
    normaliz.resize(ens.indDecKin.max(),djvec_t{(size_t)ens.T/2+1});
    
    load();
  }
  
  //! Prepare the 3pts normalization
  permes_combo_t& prepare3ptsNormalization();
  
  //! Perform the 2pts fit
  permes_combo_t& fit2pts(const char* fitTag,const bool forceRechoose=false);
  
  //! Chosses the time interval for 3pts
  permes_combo_t& choose3ptsTint(const string& mesPlotsPath,const bool forceRechoose=false);
  
  //! Chosses the time interval for 2pts
  permes_combo_t& choose2ptsTint(const string& mesPlotsPath,const bool forceRechoose=false);
  
  //! Perform the 3pts fit
  permes_combo_t& fit3pts(const char* fitTag,const bool forceRechoose=false);
  
  //! X study
  void babababab()
  {
    grace_file_t test("/tmp/o.xmg");
    const double m=this->E[0].ave();
    
    double thS=0.0;
    double thT=0.0;
    
    for(double x=1e-6;x<=1+1e-6;x+=0.05)
      {
	// test.write_line([=](const double& thS){return x-ens.getX(m,thS,thT,th0);},-100,100,1000);
	
	class thFinder : public minimizer_fun_t
	{
	public:
	  
	  const double x;
	  const perens_t& ens;
	  const double m;
	  
	  double y(const vector<double>& p,const double& th0) const
	  {
	    const double& thS=p[0];
	    const double& thT=p[1];
	    
	    return
	      x-ens.getX(m,thS,thT,th0);
	  };
	  
	  double th0forX(const vector<double>& p) const
	  {
	    double a=0;
	    double b=0;
	    
	    while(y(p,a)*y(p,b)>=0)
	      {
		//cout<<a<<" "<<b<<"    "<<y(p,a)<<" "<<y(p,b)<<endl;
		a=-b;
		b=(b+1e-3)*1.01;
		
		if(b>1e3) CRASH("No solution found");
	      };
	    
	    return
	      Brent_solve([this,&p](const double& x){return y(p,x);},a,b);
	  };
	  
	  //! Constructor
	  thFinder(const double x,const perens_t& ens,const double m) : x(x),ens(ens),m(m)
	  {
	  }
	  
	  double operator()(const double& thS,const double& thT,const double& th0) const
	  {
	    const double y=sqr(thS)+sqr(thT)+sqr(th0);
	    
	    return y;
	  }
	  
	  double operator()(const vector<double> &p) const
	  {
	    const double th0=th0forX(p);
	    //cout<<p[0]<<" "<<p[1]<<" "<<th0<<endl;
	    return (*this)(p[0],p[1],th0);
	  }
	  
	  double Up() const {return 1;}
	};
	
	thFinder thetas(x,this->ens,m);
	
	minimizer_pars_t pars;
	pars.add("thS",thS,0.1);
	pars.add("thT",thT,0.1);
	
	minimizer_t min(thetas,pars);
	
	const vector<double> minPars=min.minimize();
	cout<<x<<endl;
	const double thS=minPars[0];
	const double thT=minPars[1];
	const double th0=thetas.th0forX(minPars);
	cout<<"ThS: "<<thS<<endl;
	cout<<"ThT: "<<thT<<endl;
	cout<<"Th0: "<<th0<<endl;
	cout<<"x: "<<this->ens.getX(m,thS,thT,th0)<<endl;
	cout<<"func: "<<thetas(minPars)<<endl;
	cout<<"x: "<<this->ens.getX(m,0,0,0)<<endl;
	cout<<"func: "<<thetas(0.0,0.0,0.0)<<endl;
      }
  }

};

#endif
