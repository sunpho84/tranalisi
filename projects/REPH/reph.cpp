#include <tranalisi.hpp>

//               ~~~~~~~~~~~~~~~~~~~~~~~~~~ T/2
//               }                           |
//     --- -> Q0 --- X --- -> Qt ---         |
//    /                             \        |
//   /                               \       |
// P ---------- <- QS -------------V-A (t)

#include <reph.hpp>

using AllMesCombos=vector<permes_combo_t<>>;
using AllMes=vector<AllMesCombos>;

void readQuarkList(raw_file_t& input)
{
  const size_t nQuarks=input.read<size_t>("NQuarks");
  for(auto& q : {"Quark","Charge","IMPhys","NMasses","IMasses"})
    input.expect(q);
  
  for(size_t iQuark=0;iQuark<nQuarks;iQuark++)
    {
      const string name=input.read<string>();
      
      const double ch=input.read<double>();
      get<0>(quarkList[name])=ch;
      
      const size_t iMphys=input.read<size_t>();
      get<1>(quarkList[name])=iMphys;
      
      const size_t nMasses=input.read<size_t>();
      for(size_t iMass=0;iMass<nMasses;iMass++)
	{
	  const size_t jMass=input.read<size_t>();
	  get<2>(quarkList[name]).push_back(jMass);
	}
    }
  
  //Report
  for(const auto& q : quarkList)
    {
      cout<<q.first<<" "<<get<0>(q.second)<<"/3 "<<" iMphys: "<<get<1>(q.second);
      for(const auto& i : get<2>(q.second))
	cout<<" "<<i;
      cout<<endl;
    }
}

void readMesonList(raw_file_t& input)
{
  const size_t nMesons=input.read<size_t>("NMesons");
  for(auto& q : {"Name","Bw","Fw"})
    input.expect(q);
  
  for(size_t iMeson=0;iMeson<nMesons;iMeson++)
    {
      const string name=input.read<string>();
      const string quark1=input.read<string>();
      const string quark2=input.read<string>();
      
      for(const auto& q : {quark1,quark2})
	if(quarkList.find(q)==quarkList.end())
	  CRASH("Unable to find quark %s",q.c_str());
      
      mesonList.push_back({name,quark1,quark2});
    }
}

void readInput(const string& path)
{
  raw_file_t input(path,"r");
  
  readQuarkList(input);
  
  readMesonList(input);
}

//! Read and compute ll combination for a given meson
AllMesCombos computeAllMesCombos(vector<double>& mS,vector<double>& mT,const perens_t& ens,const meson_t& mes)
{
  //! Result
  AllMesCombos res;
  
  //! Name of the meson
  const string& mesName=get<0>(mes);
  
  //! Name of spectator quark
  const quark_t& quarkS=quarkList[get<1>(mes)];
  
  //! Name of other quark
  const quark_t& quarkT=quarkList[get<2>(mes)];
  
  //! Charge of the spectator and forward line
  const double eS=get<0>(quarkS)/3.0;
  const double eT=get<0>(quarkT)/3.0;
  
  //! Loop on all meson combos
  const size_t nMs=get<2>(quarkS).size();
  const size_t nMt=get<2>(quarkT).size();
  const index_t indMesCombo({{"iMs",nMs},{"iMt",nMt}});
  
  mS.resize(nMs);
  mT.resize(nMt);
  
  const size_t nMesCombos=indMesCombo.max();
  
  //! Setup all meson combos
  for(size_t iMesCombo=0;iMesCombo<nMesCombos;iMesCombo++)
    {
      const vector<size_t> c=indMesCombo(iMesCombo);
      const size_t iMs=get<2>(quarkS)[c[0]];
      const size_t iMt=get<2>(quarkT)[c[1]];
      
      cout<<endl<<" === "<<mesName<<" "<<indMesCombo.descr(iMesCombo)<<" ==="<<endl;
      
      res.emplace_back(ens,mesName,iMs,iMt,eS,eT);
      
      mS[c[0]]=ens.mass[iMs];
      mT[c[1]]=ens.mass[iMt];
    }
  
  for(size_t iMesCombo=0;iMesCombo<nMesCombos;iMesCombo++)
    res[iMesCombo]
      .fit2pts("selfChosenTint")
      .prepare3ptsNormalization()
      .fit3pts("selfChosenTint")
      .plotFf();
  
  return res;
}

int main(int narg,char **arg)
{
  readInput("input.txt");
  
  loadUltimateInput("ultimate_input.txt");
  
  set_njacks(15);
  def_nboots=nboots;
  
  perens_t ens(".");
  const size_t& iBeta=ens.iBeta;
  cout<<"Beta: "<<iBeta<<endl;
  
  const size_t nMes=mesonList.size();
  
  //! Holds ff energy etc for each meson and combination
  AllMes mesCombos(nMes);
  
  for(size_t iMes=0;iMes<nMes;iMes++)
    {
      const meson_t& mes=mesonList[iMes];
      
      vector<double> amS,amT;
      
      mesCombos[iMes]=computeAllMesCombos(amS,amT,ens,mes);
      
      //! Index of S quark in the physical list
      const size_t iMphysS=get<1>(quarkList[get<1>(mes)]);
      
      //! Index of T quark in the physical list
      const size_t iMphysT=get<1>(quarkList[get<2>(mes)]);
      
      const size_t iult=0;
      
      //! Bare mass of S quark
      const dboot_t amSbare=lat_par[iult].amBare(iMphysS,iBeta);
      
      //! Bare mass of T quark
      const dboot_t amTbare=lat_par[iult].amBare(iMphysT,iBeta);
      
      cout<<"amSbare: "<<smart_print(amSbare);
      for(auto& amSi: amS) cout<<" "<<amSi;
      cout<<endl;
      
      cout<<"amTbare: "<<smart_print(amTbare);
      for(auto& amTi: amT) cout<<" "<<amTi;
      cout<<endl;
      
      //! Name of the meson_
      const string& mesName=get<0>(mes);
      
      const permes_t<dbvec_t> inte=interpolate(ens,mesName,amS,amT,mesCombos[iMes],jack_index[0][0],amSbare,amTbare);
      inte.plotFf();
    }
  
  /*
    - * loop on all physical mesons
    - | loop on all ensembles
    - | | * read the input
    - | | * loop on all meson combo
    - | | | * read 2pts correlators
    - | | |
    - | | | loop on all algorithms to select 2pts fit range
    - | | | | * output to diferent meson combo folder
    - | | | | * do the fit of PP and AP
    - | | | | * check dispersion relation
    - | | | |
    - | | | * read V and A correlation functions
    - | | | | loop on all ratios definitions
    - | | | | | * build the ratios
    - | | | | | ..
    - | | | | |
    - | | | | | loop on all 3pts fit range
    - | | | | | | *fit 3pts
    - | | | |
    - | loop on all systematics so far
    - | | loop on all 8 analysis
    - | | | * interpolate to physical meson_corrs
    - | | | fit all ensembles with a multilinear function
   */
  
  return 0;
}
