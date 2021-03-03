#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <tranalisi.hpp>

#include "data.hpp"
#include "node.hpp"

#include "fitter/fitter_parser.hpp"
#include "operations.hpp"

using namespace std;

const size_t nArgsPerFile=5;

///Exit with help message
void exitWithHelp(const char* execName)
{
  CRASH("Use %s fileData xmin xmax formula filePlot [ fileData xmin xmax formula filePlot ... ] [ -- ansatz ]",execName);
}

vector<data_t> dataFit;
vector<size_t> parIds;
map<size_t,size_t> idOfPars;
vector<double> guessPars;

/// Parse all arguments
void parseArgs(int narg,char** arg)
{
  int iDivAnsatz=0;
  while(iDivAnsatz<narg and strcasecmp(arg[iDivAnsatz],"--"))
    iDivAnsatz++;
  cout<<"DivAnsatz: "<<iDivAnsatz<<endl;
  
  /// Number of files to be fit
  const int nFiles=(iDivAnsatz-1)/nArgsPerFile;
  if(iDivAnsatz!=(int)(1+nFiles*nArgsPerFile))
    exitWithHelp(arg[0]);
  cout<<"NFiles: "<<nFiles<<endl;
  if(nFiles<1)
    CRASH("Please provide at least one file to fit");
  
  set<size_t> parIdsSet;
  for(int iFile=0;iFile<nFiles;iFile++)
    {
      enum{PathData=0,Xmin=1,Xmax=2,Formula=3,PathPlot=4};
      char** pars=arg+1+iFile*nArgsPerFile;
      
      const string pathData=pars[PathData];
      const double xMin=strtod(pars[Xmin],nullptr);
      const double xMax=strtod(pars[Xmax],nullptr);
      const string formula=pars[Formula];
      const string pathPlot=pars[PathPlot];
      
      dataFit.emplace_back(pathData,xMin,xMax,formula,pathPlot);
      
      for(const size_t& i : dataFit.back().formula->getParIds())
	parIdsSet.insert(i);
    }
  
  parIds.reserve(parIdsSet.size());
  for(auto& i : parIdsSet)
    {
      idOfPars[i]=parIds.size();
      parIds.push_back(i);
    }
  
  for(size_t iArg=iDivAnsatz+1;iArg<(size_t)narg;iArg++)
    guessPars.push_back(strtod(arg[iArg],nullptr));
}

size_t ijack;

//! perform a simple fit using x, a function and data
class ch2_t : public minimizer_fun_t
{
public:
  //! constructor
  ch2_t()
  {
  }
  
  //! compute the function
  double operator()(const vector<double> &p) const
  {
    double ch2=0;
    
    size_t ntot_contr=0;
    for(auto& d : dataFit)
      {
	size_t ncontr=0;
	for(const auto& c : d.points)
	  {
	    const double& x=c.x;
	    const double& n=c.y[ijack];
	    const double& e=c.e;
	    const double t=d.formula->eval(p,x);
	    double contr=sqr((n-t)/e);
	    // cout<<contr<<" = [("<<n<<"-f("<<x<<")="<<t<<")/"<<e<<"]^2]"<<endl;
	    if(x>=d.xmin and x<=d.xmax)
	      {
		ncontr++;
		ch2+=contr;
	      }
	    // else
	    //   cout<<"Discarded, "<<x<<" not in range: "<<d.xmin<<" "<<d.xmax<<endl;
	  }
	ntot_contr+=ncontr;
	// cout<<"NContr: "<<ncontr<<endl;
      }
    // cout<<"ch2: "<<ch2<<endl;
    
    return ch2;
  }
  
  double Up() const {return 1;}
};

int main(int narg,char** arg)
{
  set_njacks(100);
  
  parseArgs(narg,arg);
  
  const size_t nPars=parIds.size();
  cout<<"NPars: "<<nPars<<endl;
  guessPars.resize(nPars);
  
  minimizer_pars_t fitPars;
  for(size_t iPar=0;iPar<nPars;iPar++)
    fitPars.add(combine("[%zu]",parIds[iPar]),guessPars[iPar],fabs(guessPars[iPar]/10));
  
  ch2_t ch2;
  
  djvec_t pars(nPars);
  djack_t ch2min;
  for(ijack=0;ijack<=njacks;ijack++)
    {
      minimizer_t minimizer(ch2,fitPars);
      vector<double> parsEl=minimizer.minimize();
      for(size_t iPar=0;iPar<nPars;iPar++)
	pars[iPar][ijack]=parsEl[iPar];
      
      ch2min[ijack]=ch2(parsEl);
    }
  cout<<"Ch2: "<<ch2min[njacks]<<endl;
  
  cout<<pars.ave_err()<<endl;
  
  for(auto& data : dataFit)
    {
      grace_file_t plot(data.filePlot);
      
      for(auto& d : data.points)
	plot.write_ave_err(d.x,{d.y.ave(),d.e});
      
      plot.write_polygon([&](const double& x){return data.formula->jackEval(pars,x);},data.xmin,data.xmax,grace::GREEN4);
    }
  
  return 0;
}

