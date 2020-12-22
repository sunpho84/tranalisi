#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <tranalisi.hpp>

#include "node.hpp"

#include "fitter/fitter_parser.hpp"
#include "operations.hpp"

using namespace std;

node_t* parseFormula(const string& formulaString);

const size_t nArgsPerFile=5;

///Exit with help message
void exitWithHelp(const char* execName)
{
  CRASH("Use %s fileData xmin xmax formula filePlot [, fileData xmin xmax formula filePlot, ... ] [ -- ansatz ]",execName);
}

/// Data
struct data_fit_t
{
  const string fileData;
  const double xmin,xmax;
  const node_t* formula;
  const string filePlot;
  
  vector<array<double,3>> data;
  
  data_fit_t(const string& fileData,const double xmin,const double xmax,const string& formulaString,const string& filePlot) :
    fileData(fileData),xmin(xmin),xmax(xmax),formula(parseFormula(formulaString)),filePlot(filePlot)
  {
    /// Open the file
    obs_file_t fin(fileData,3,{0,1,2});
    
    /// Count the lines
    const size_t nLines=fin.length(1)/3;
    data.resize(nLines);
    
    for(size_t iLine=0;iLine<nLines;iLine++)
      {
	/// Single line
	const vector<double> l=fin.read(1);
	for(size_t i=0;i<3;i++) data[iLine][i]=l[i];
      }
  }
};

vector<data_fit_t> dataFit;
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
    
    for(auto& d : dataFit)
      for(const auto& c : d.data)
	{
	  const double& x=c[0];
	  const double& n=c[1];
	  const double& e=c[2];
	  const double t=d.formula->eval(p,x);
	  double contr=sqr((n-t)/e);
	  // cout<<contr<<" = [("<<n<<"-f("<<x<<")="<<t<<")/"<<e<<"]^2]"<<endl;
	  if(x>=d.xmin and x<=d.xmax)
	    ch2+=contr;
	}
    
    // cout<<"ch2: "<<ch2<<endl;
    
    return ch2;
  }
  
  double Up() const {return 1;}
};

int main(int narg,char** arg)
{
  parseArgs(narg,arg);
  
  const size_t nPars=parIds.size();
  cout<<"NPars: "<<nPars<<endl;
  guessPars.resize(nPars);
  
  minimizer_pars_t fitPars;
  for(size_t iPar=0;iPar<nPars;iPar++)
    fitPars.add(combine("[%zu]",parIds[iPar]),guessPars[iPar],fabs(guessPars[iPar]/10));
  
  ch2_t ch2;
  
  minimizer_t minimizer(ch2,fitPars);
  
  vector<double> pars=minimizer.minimize();
  cout<<pars<<endl;
  
  for(auto& data : dataFit)
    {
      grace_file_t plot(data.filePlot);
      
      for(auto& d : data.data)
	plot.write_ave_err(d[0],{d[1],d[2]});
      
      plot.write_line([&data,&pars](const double& x){return data.formula->eval(pars,x);},data.xmin,data.xmax,grace::RED);
    }
  
  return 0;
}

