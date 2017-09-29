#ifndef _SYSTEMATICS_HPP
#define _SYSTEMATICS_HPP

#include <vector>

#include <index.hpp>

using namespace std;

class syst_t
{
public:
  double ave; //!< average
  double wave; //!< weighted average
  double stat; //!< pure statistics
  double wstat; //!< weighted statistics
  double totsy; //!< total systematics
  double tot; //!< all statistics and systematics summed
  vector<double> syst; //!< separated systematics
};

vector<double> syst_analysis_sep_bis(const vector<ave_err_t> &v,const index_t &ind);

//! perform the eq.28-like analysis
template <class TV>
syst_t perform_analysis(const TV &data,const index_t &ind)
{
  syst_t r;
  double &ave=r.ave=0,&err=r.stat=0,&wave=r.wave=0,&werr=r.wstat=0,weight=0;
  for(auto &d : data)
    {
      double a=d.ave(),e=d.err(),w=1/sqr(e);
      ave+=a;
      err+=e;
      wave+=a*w;
      werr+=e*w;
      weight+=w;
    }
  ave/=data.size();
  err/=data.size();
  wave/=weight;
  werr/=weight;
  
  r.syst=syst_analysis_sep_bis(data.ave_err(),ind);
  double &totsy=r.totsy=0,&tot=r.tot=0;
  for(auto &s : r.syst) totsy+=s*s;
  totsy=sqrt(totsy);
  tot=sqrt(sqr(totsy)+sqr(err));
  
  return r;
}

//! write average, errors and systematics
template <class TV>
void perform_analysis(const TV &data,const index_t &ind,const string &name)
{
  cout<<" ================== "<<name<<" ================== "<<endl;
  auto r=perform_analysis(data,ind);
  
  for(auto el : vector<pair<string,double>>{{" Ave",r.ave},{"Weig",r.wave}})
    cout<<el.first.substr(0,6)<<":\t"<<el.second<<endl;
  
  vector<pair<string,double>> top;
  top.push_back({" Stat",r.stat});
  top.push_back({"Werr",r.wstat});
  for(size_t i=0;i<ind.rank();i++) top.push_back({ind.name(i),r.syst[i]});
  top.push_back({"TotSy",r.totsy});
  top.push_back({" Tot",r.tot});
  
  for(auto t : top)
    {
      cout<<t.first.substr(0,6);
      cout<<"\t=\t";
      cout<<t.second;
      cout<<"\t=\t";
      streamsize prec=cout.precision();
      cout.precision(3);
      cout<<fabs(t.second/r.ave)*100<<" %";
      cout.precision(prec);
      cout<<endl;
    }
  cout<<endl;
  cout<<"In sysntesis: "<<smart_print(r.ave,r.syst)<<endl;
}

#endif
