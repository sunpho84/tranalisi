#include <tranalisi.hpp>

//               ~~~~~~~~~~~~~~~~~~~~~~~~~~ T/2
//               }                           |
//     --- -> Q0 --- X --- -> Qt ---         |
//    /                             \        |
//   /                               \       |
// P ---------- <- QS -------------V-A (t)

#include <reph.hpp>

void readQuarkList(raw_file_t& input)
{
  const size_t nQuarks=input.read<size_t>("NQuarks");
  for(auto& q : {"Quark","Charge","NMasses","IMasses"})
    input.expect(q);
  
  for(size_t iQuark=0;iQuark<nQuarks;iQuark++)
    {
      const string name=input.read<string>();
      
      const double ch=input.read<double>();
      get<0>(quarkList[name])=ch;
      
      const size_t nMasses=input.read<size_t>();
      for(size_t iMass=0;iMass<nMasses;iMass++)
	{
	  const size_t jMass=input.read<size_t>();
	  get<1>(quarkList[name]).push_back(jMass);
	}
    }
  
  //Report
  for(const auto& q : quarkList)
    {
      cout<<q.first<<" "<<get<0>(q.second)<<"/3 ";
      for(const auto& i : get<1>(q.second))
	cout<<i<<" ";
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

int main(int narg,char **arg)
{
  readInput("input.txt");
  
  set_njacks(15);
  
  perens_t ens(".");
  
  const size_t nMes=mesonList.size();
  
  //! Holds ff energy etc for each meson and combination
  vector<vector<permes_combo_t<>>> mesCombos(nMes);
  
  for(size_t iMes=0;iMes<nMes;iMes++)
    {
      const meson_t& mes=mesonList[iMes];
      
      const string& mesName=get<0>(mes);
      cout<<"Meson: "<<mesName<<endl;
      
      const quark_t& quarkS=quarkList[get<1>(mes)];
      const quark_t& quarkT=quarkList[get<2>(mes)];
      
      //! Charge of the spectator and forward line
      const double eS=get<0>(quarkS)/3.0;
      const double eT=get<0>(quarkT)/3.0;
      
      //! Loop on all meson combos
      const size_t nMs=get<1>(quarkS).size();
      const size_t nMt=get<1>(quarkT).size();
      const index_t indMesCombo({{"iMs",nMs},{"iMt",nMt}});
      vector<double> mS(nMs),mT(nMt);
      
      const size_t nMesCombos=indMesCombo.max();
      
      //! Setup all meson combos
      for(size_t iMesCombo=0;iMesCombo<nMesCombos;iMesCombo++)
	{
	  const vector<size_t> c=indMesCombo(iMesCombo);
	  const size_t iMs=get<1>(quarkS)[c[0]];
	  const size_t iMt=get<1>(quarkT)[c[1]];
	  
	  cout<<" === "<<mesName<<" "<<indMesCombo.descr(iMesCombo)<<" ==="<<endl;
	  
	  mesCombos[iMes].emplace_back(ens,mesName,iMs,iMt,eS,eT);
	  
	  mS[c[0]]=ens.mass[iMs];
	  mT[c[1]]=ens.mass[iMt];
	}
      
      for(size_t iMesCombo=0;iMesCombo<nMesCombos;iMesCombo++)
	mesCombos[iMes][iMesCombo]
	  .fit2pts("selfChosenTint")
	  .prepare3ptsNormalization()
	  .fit3pts("selfChosenTint")
	  .plotFf();
      
      cout<<"Interpolating"<<endl;
      permes_t<> mesInterpolated(ens,combine("mes_%s",mesName.c_str()));
      mesInterpolated.interpolate(mS,mT,mesCombos[iMes]);
      
      for(size_t i=0;i<mesInterpolated.ff[0].size();i++)
	{
	  cout<<i<<" "<<mesInterpolated.ff[0][i]<<" "<<mesCombos[iMes][0].ff[0][i]<<endl;
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
    - | | | * interpolate to physical meson_corrs
    - | | | fit all ensembles with a multilinear function
   */
  
  return 0;
}
