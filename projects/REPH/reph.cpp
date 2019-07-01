#include <tranalisi.hpp>

//               ~~~~~~~~~~~~~~~~~~~~~~~~~~ T/2
//               }                           |
//     --- -> Q0 --- X --- -> Qt ---         |
//    /                             \        |
//   /                               \       |
// P ---------- <- QS -------------V-A (t)

#include <reph.hpp>

//! Reads the list of quarks
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

//! Reads the list of mesons
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

//! Reads the input file
void readInput(const string& path)
{
  //! Input file
  raw_file_t input(path,"r");
  
  readQuarkList(input);
  
  readMesonList(input);
}

//! Reads and compute ll combination for a given meson
AllMesCombos computeAllMesCombos(const perens_t& ens,const meson_t& mes)
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
  
  const size_t nMesCombos=indMesCombo.max();
  
  //! Setup all meson combos
  for(size_t iMesCombo=0;iMesCombo<nMesCombos;iMesCombo++)
    {
      const vector<size_t> c=indMesCombo(iMesCombo);
      const size_t iMs=get<2>(quarkS)[c[0]];
      const size_t iMt=get<2>(quarkT)[c[1]];
      
      cout<<endl<<" === "<<mesName<<" "<<indMesCombo.descr(iMesCombo)<<" ==="<<endl;
      
      res.emplace_back(ens,mesName,iMs,iMt,eS,eT);
      
      res[iMesCombo]
	.fit2pts("selfChosenTint")
	.prepare3ptsNormalization()
	.fit3pts("selfChosenTint")
	.plotFf();
      
      cout<<endl;
    }
  
  return res;
}

//! Perform a loop over ensemble, putting a proper header
template <typename F>
void mesonLoop(const F& f)
{
  const size_t nMes=mesonList.size();
  
  for(size_t iMes=0;iMes<nMes;iMes++)
    {
      const meson_t& mesComposition=mesonList[iMes];
      cout<<endl<<"//// Meson: "<<get<0>(mesComposition)<<" ////"<<endl;
      
      f(mesComposition);
    }
}

//! Perform a loop over ensemble, putting a proper header
template <typename F>
void ensembleLoop(const vector<perens_t>& ens,const F& f)
{
  const size_t& nens=ens.size();
  
  for(size_t iens=0;iens<nens;iens++)
    {
      cout<<endl<<"**** "<<ens[iens].dirPath<<" ****"<<endl<<endl;
      f(ens[iens],iens);
    }
}

//! Perform a loop over ultimate analysis
template <typename F>
void ultimateAnalysisLoop(const F& f)
{
  for(size_t inputAn=0;inputAn<ninput_an;inputAn++)
    {
      cout<<endl<<"//// Input analysis: "<<inputAn<<" ////"<<endl<<endl;
      
      f(inputAn);
    }
}

int main(int narg,char **arg)
{
  readInput("input.txt");
  
  loadUltimateInput("ultimate_input.txt");
  
  set_njacks(15);
  def_nboots=nboots;
  
  const vector<perens_t> ens=readEnsList();
  
  // const size_t nMes=mesonList.size();
  
  mesonLoop([&](const meson_t& mesComposition)
	    {
	      //! Holds ff energy etc for each meson and combination, for each ensemble
	      vector<AllMesCombos> mesCombos;
	      ensembleLoop(ens,[&](const perens_t& e,size_t){mesCombos.push_back(computeAllMesCombos(e,mesComposition));});
	      
	      // Loop over analysis
	      ultimateAnalysisLoop([&](const size_t& inputAn)
				   {
				     //! Interpolated data for each ensemble
				     vector<permes_t<dbvec_t>> inte;
				     
				     ensembleLoop(ens,[&](const perens_t& e,const size_t& iens)
						      {
							const map<string,vector<double>> am=getMassList(e);
							
							inte.emplace_back(interpolate(mesCombos[iens],mesComposition,e,am,inputAn));
							inte.back().plotFf();
						      });
				     
				   });
	    });
  
    /*
    - * loop on all physical mesons
    - | * loop on all ensembles
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
