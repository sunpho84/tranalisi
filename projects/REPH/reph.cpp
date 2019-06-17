#include <tranalisi.hpp>

//               ~~~~~~~~~~~~~~~~~~~~~~~~~~ T/2
//               }                           |
//     --- -> Q0 --- X --- -> Qt ---         |
//    /                             \        |
//   /                               \       |
// P ---------- <- QS -------------V-A (t)

#include <reph.hpp>

void readInput(const string& path)
{
  raw_file_t input(path,"r");
  
  const size_t nQuarks=input.read<size_t>("NQuarks");
  for(auto& q : {"Quark","Charge","NMasses","IMasses"})
    input.expect(q);
  
  for(size_t iQuark=0;iQuark<nQuarks;iQuark++)
    {
      const string name=input.read<string>();
      
      const double ch=input.read<double>();
      get<0>(quarks[name])=ch;
      
      const size_t nMasses=input.read<size_t>();
      for(size_t iMass=0;iMass<nMasses;iMass++)
	{
	  const size_t jMass=input.read<size_t>();
	  get<1>(quarks[name]).push_back(jMass);
	}
    }
  
  //Report
  for(const auto& q : quarks)
    {
      cout<<q.first<<" "<<get<0>(q.second)<<"/3 ";
      for(const auto& i : get<1>(q.second))
	cout<<i<<" ";
      cout<<endl;
    }
  
  const size_t nMesons=input.read<size_t>("NMesons");
  for(auto& q : {"Name","Bw","Fw"})
    input.expect(q);
  
  for(size_t iMeson=0;iMeson<nMesons;iMeson++)
    {
      const string name=input.read<string>();
      const string quark1=input.read<string>();
      const string quark2=input.read<string>();
      
      for(const auto& q : {quark1,quark2})
	if(quarks.find(q)==quarks.end())
	  CRASH("Unable to find quark %s",q.c_str());
      
      mesons[name]=make_pair(quark1,quark2);
    }
}

int main(int narg,char **arg)
{
  readInput("input.txt");
  
  set_njacks(15);
  
  perens_t ens(".");
  
  for(const auto& mes : mesons)
    {
      const string& mesName=mes.first;
      cout<<"Meson: "<<mesName<<endl;
      
      const quark_t& quarkS=quarks[get<1>(mes).first];
      const quark_t& quarkT=quarks[get<1>(mes).second];
      
      //! Charge of the spectator and forward line
      const double eS=get<0>(quarkS)/3.0;
      const double eT=get<0>(quarkT)/3.0;
      
      //! Loop on all meson combos
      const size_t nMs=get<1>(quarkS).size();
      const size_t nMt=get<1>(quarkT).size();
      const index_t indMesCombo({{"iMs",nMs},{"iMt",nMt}});
      
      //! Holds the 2pts info for each meson combination
      vector<permes_combo_t> mesCombos;
      const size_t nMesCombos=indMesCombo.max();
      
      //! Setup all meson combos
      for(size_t iMesCombo=0;iMesCombo<nMesCombos;iMesCombo++)
	{
	  const vector<size_t> c=indMesCombo(iMesCombo);
	  const size_t iMs=get<1>(quarkS)[c[0]];
	  const size_t iMt=get<1>(quarkT)[c[1]];
	  mesCombos.emplace_back(ens,mesName,iMs,iMt,eS,eT);
	}
      
      for(size_t iMesCombo=0;iMesCombo<nMesCombos;iMesCombo++)
	{
	  mesCombos[iMesCombo]
	    .fit2pts("fixed_tfit")
	    .prepare3ptsNormalization()
	    .chooseTint()
	    .fit3pts("tfixed_tfit")
	    .plotFf();
	}
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
    - | | | interpolate to physical meson_corrs
    - | | | fit all ensembles with a multilinear function
   */
  
  return 0;
}
