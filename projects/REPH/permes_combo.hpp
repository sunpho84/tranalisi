#ifndef _PERMES_COMBO_HPP
#define _PERMES_COMBO_HPP

#include <REPH/perens.hpp>
#include <REPH/phys_point.hpp>

#include <set>

//! Hold all data to draw ff
template <typename TV=djvec_t>
struct permes_t
{
  //! Scalar type
  using T=remove_reference_t<decltype(TV{}[0])>;
  
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
  
  //! Maximal X
  T xMax() const
  {
    T m;
    bool isFirst=true;
    for(size_t iDecKin=0;iDecKin<ens.nDecKin;iDecKin++)
      if(isFirst or m.ave()<X[iDecKin].ave())
	{
	  isFirst=false;
	  m=X[iDecKin];
	}
    
    return m;
  }
  
  permes_t(const perens_t& ens,const string& mesTag) : ens(ens),mesTag(mesTag){}
  
  //! Form factors for V and A
  array<TV,2> ff;
  
  //! Quality of the point
  array<vector<size_t>,2> quality;
  
  //! Plot ff
  void plotFf(const string& tag="") const;
  
  //! Return a bootstrap version
  static permes_t<dbvec_t> getBoot(const permes_t<djvec_t>& in,const boot_init_t& jack_of_boot)
  {
    permes_t<dbvec_t> out(in.ens,in.mesTag);
    
    out.E=bvec_from_jvec(jack_of_boot,in.E);
    out.X=bvec_from_jvec(jack_of_boot,in.X);
    for(int i=0;i<2;i++)
      out.ff[i]=bvec_from_jvec(jack_of_boot,in.ff[i]);
    
    return out;
  }
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
  
  //! Effective coupling
  vector<djvec_t> zEff;
  
  //! Difference of energies between initial and final state
  djvec_t dEdec;
  
  //! Quadrimomentum product
  djvec_t PKdec;
  
  //! Pseudoscalar correlation function
  vector<djvec_t> corrPP;
  
  //! T component of axial current
  vector<djvec_t> corrA0P;
  
  //! Z component of axial current
  vector<djvec_t> corrA3P;
  
  //! Decay correlators for V and A
  array<vector<djvec_t>,2> corrPX;
  
  //! Time interval for 3pts fit
  array<vector<Range>,2> tint3pts;
  
  //! Common range for all three points
  array<Range,2> commonTint3pts;
  
  //! Time interval for 2pts fit
  vector<Range> tint2pts;
  
  //! Load the PP correlation function
  djvec_t load2ptsPP(const size_t iMoms,const size_t iMomt);
  
  //! Load the AP correlation function
  djvec_t load2ptsAP(const size_t iMoms,const size_t iMomt,const size_t iGamma);
  
  //! Load the three points correlation functions
  djvec_t load3pts(const string& plotDirPath,const size_t iVA,const size_t iMs,const size_t iMt,const size_t iMoms,const size_t iMomt,const size_t iMom0,const char* extraTag);
  
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
  
  //! Print all kinematics info
  void printKin() const
  {
    const string path=ens.dirPath+"/plots/"+mesTag+"/kinematics.txt";
    ofstream out(path);
    if(not out.good())
      CRASH("Unable to create output %s",path.c_str());
	
    for(size_t iDecKin=0;iDecKin<ens.nDecKin;iDecKin++)
      {
	const int iMesKin=ens.iMesKinOfDecKin[iDecKin];
	
	out<<ens.indDecKin.descr(iDecKin)<<endl;
	out<<" X: "<<X[iDecKin].ave_err()<<endl;
	out<<" Pmes_z: "<<ens.pMes[iMesKin]<<endl;
	out<<" Khat_z: "<<ens.kHatDec[iDecKin]<<endl;
	out<<" PK: "<<PKdec[iDecKin].ave_err()<<endl;
	out<<" Eg: "<<ens.Eg[iDecKin]<<endl;
	out<<" EgT: "<<ens.EgT(iDecKin)<<endl;
	out<<endl;
      }
    
    out<<"/////////////////////////////////////////////////////////////////"<<endl;
    
    for(size_t iMesKin=0;iMesKin<ens.nMesKin;iMesKin++)
      {
	out<<ens.indMesKin.descr(iMesKin)<<endl;
	out<<" E: "<<E[iMesKin].ave_err()<<endl;
	out<<" fP: "<<fP[iMesKin].ave_err()<<endl;
	out<<" fPbare: "<<fPbare[iMesKin].ave_err()<<endl;
	out<<" ZP: "<<ZP[iMesKin].ave_err()<<endl;
	out<<" ZA: "<<ZA[iMesKin].ave_err()<<endl;
	out<<endl;
      }
  }
  
  //! Plots the form factors
  permes_combo_t& plotFf(const string& tag="")
  {
    static_cast<permes_t<TV>*>(this)->plotFf(tag);
    
    return *this;
  }
  
  //! Constructor
  permes_combo_t(const perens_t& ens,const string& mesName,const size_t& iMs,const size_t& iMt,const double& eS,const double& eT) :
    permes_t<TV>(ens,combine("%s/iMs%zu" "_" "iMt%zu",mesName.c_str(),iMs,iMt)),
    iMs(iMs),
    iMt(iMt),
    eS(eS),
    eT(eT)
  {
    resizeListOfContainers({&ZP,&ZA,&fP,&fPbare,&this->E},ens.nMesKin);
    resizeListOfContainers({&eEff,&zEff},ens.nMesKin);
    resizeListOfContainers({&tint2pts},ens.nMesKin);
    
    resizeListOfContainers({&corrA0P,&corrA3P,&corrPP},ens.nMesKin,djvec_t{(size_t)ens.T/2+1});
    
    resizeListOfContainers({&this->ff[0],&this->ff[1]},ens.nDecKin);
    resizeListOfContainers({&this->quality[0],&this->quality[1]},ens.nDecKin,true);
    
    for(size_t iVA=0;iVA<2;iVA++)
      {
	resizeListOfContainers({&corrPX[iVA]},ens.nDecKin,djvec_t{(size_t)ens.T/2+1});
	resizeListOfContainers({&tint3pts[iVA]},ens.nDecKin,Range{0,0});
      }
    
    resizeListOfContainers({&dEdec,&PKdec,&this->X},ens.nDecKin);
    
    load();
  }
  
  //! Prepare the 3pts normalization
  permes_combo_t& prepareKinematics(const bool useAnalytic);
  
  //! Perform the 2pts fit
  permes_combo_t& fit2pts(const char* fitTag,const bool forceRechoose=false);
  
  //! Chooses the time interval for 3pts
  permes_combo_t& choose3ptsTint(const string& mesPlotsPath,const char* fitTag,const bool forceRechoose=false);
  
  //! Chooses a common range for three points
  permes_combo_t& choose3ptsTintCommon(const size_t& length=4);
  
  //! Chosses the time interval for 2pts
  permes_combo_t& choose2ptsTint(const string& mesPlotsPath,const bool forceRechoose=false);
  
  //! Gets the correlation ratio
  TV getCorrRat(const int iVA,const size_t iDecKin);
  
  //! Perform the 3pts fit
  permes_combo_t& fit3pts(const bool& useCommonRange,const char* fitTag,const bool forceRechoose=false);
  
  //! X study
  void babababab()
  {
    grace_file_t test("/tmp/o.xmg");
    const double m=this->E[0].ave();
    
    double thS=0.0;
    double thT=0.0;
    
    for(double x=1e-6;x<=1+1e-6;x+=0.05)
      {
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

//! All combinations of a physical meson
using AllMesCombos=vector<permes_combo_t<>>;

#endif
