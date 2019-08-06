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
  for(auto& q : {"Name","Bw","Fw","Mass"})
    input.expect(q);
  
  for(size_t iMeson=0;iMeson<nMesons;iMeson++)
    {
      const string name=input.read<string>();
      const string quark1=input.read<string>();
      const string quark2=input.read<string>();
      const double mass=input.read<double>();
      
      for(const auto& q : {quark1,quark2})
	if(quarkList.find(q)==quarkList.end())
	  CRASH("Unable to find quark %s",q.c_str());
      
      mesonList.push_back({name,quark1,quark2,mass});
    }
}

//! Reads the physics file
void readPhysics(const string& path)
{
  //! Input file
  raw_file_t input(path,"r");
  
  readQuarkList(input);
  
  readMesonList(input);
}

//! Reads and compute ll combination for a given meson
AllMesCombos computeAllMesCombos(const perens_t& ens,const meson_t& mes,const size_t& timeDependentEnergy,const size_t& useCommonRange,const string& totTag)
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
      
      const bool DO_NOT_USE_ANALYTIC=true;
      
      res[iMesCombo]
	.fit2pts("selfChosenTint")
	.prepare3ptsNormalization(DO_NOT_USE_ANALYTIC,timeDependentEnergy,totTag)
	.fit3pts(useCommonRange,totTag.c_str())
	.plotFf(totTag);
      
      res[iMesCombo].printKin();
      
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

//! Perform a loop over decay kinematic
template <typename F>
void decKinLoop(const perens_t& e,const F& f)
{
  for(size_t iDecKin=0;iDecKin<e.nDecKin;iDecKin++)
    if(e.considerDec[iDecKin])
      f(iDecKin);
}

//! Ansatz fit for fA
template <typename TV>
auto ansatzA(const TV& p,const double M,const double a2,const double x) -> std::remove_reference_t<decltype(p[0])>
{
  //using T=std::remove_reference_t<decltype(p[0])>;
  
  return M*(p[0]+p[1]*M*M+p[2]*M*M*x+p[3]*a2);
}

//! Ansatz fit for fV
template <typename TV>
auto ansatzV(const TV& p,const double M,const double a2,const double x) -> std::remove_reference_t<decltype(p[0])>
{
  //using T=std::remove_reference_t<decltype(p[0])>;
  
  return M*(p[0]+p[1]*M*M+p[2]*M*M*x+p[3]*a2);
}

//! Ansatz fit
template <typename TV>
auto ansatz(const size_t iVA,const TV& p,const double M,const double a2,const double x) -> std::remove_reference_t<decltype(p[0])>
{
  if(iVA==0)
    return ansatzV(p,M,a2,x);
  else
    return ansatzA(p,M,a2,x);
}

int main(int narg,char **arg)
{
  readPhysics("physics.txt");
  
  loadUltimateInput("ultimate_input.txt");
  
  //set_printlevel(3);
  set_njacks(15);
  def_nboots=nboots;
  
  // Initialize the renormalization constants
  for(auto& q : vector<tuple<const vector<ave_err_t>*,dbvec_t*,int>>{{&Za_ae,&Za,221533},{&Zv_ae,&Zv,854312},{&Zp_ae,&Zp,235346},{&Zt_ae,&Zt,6342423}})
    {
      const vector<ave_err_t>& zae=*get<0>(q);
      dbvec_t& z=*get<1>(q);
      const size_t& seed=get<2>(q);
      
      for(size_t i=0;i<2*nbeta;i++)
	z[i].fill_gauss(zae[i],seed*2*nbeta+i);
    }
  
  const vector<perens_t> ens=readEnsList();
  
  // const size_t nMes=mesonList.size();
  
  mesonLoop([&](const meson_t& mesComposition)
	    {
	      const double& mPhys=get<3>(mesComposition);
	      
	      
	      
	      for(size_t timeDependentEnergy=0;timeDependentEnergy<2;timeDependentEnergy++)
		for(size_t useCommonRange=0;useCommonRange<2;useCommonRange++)
		  {
		    string totTag=combine("%s_%s",rangeTag[useCommonRange],timeDepEnTag[timeDependentEnergy]);
		    
		    //! Holds ff energy etc for each meson and combination, for each ensemble
		    vector<AllMesCombos> mesCombos;
		    ensembleLoop(ens,[&](const perens_t& e,size_t){mesCombos.push_back(computeAllMesCombos(e,mesComposition,timeDependentEnergy,useCommonRange,totTag));});
		    
		    // Loop over analysis
		    ultimateAnalysisLoop([&](const size_t& inputAn)
					 {
					   const size_t iM12=inputAn/4;
					   
					   //! Interpolated data for each ensemble
					   vector<permes_t<dbvec_t>> inte;
					   
					   ensembleLoop(ens,[&](const perens_t& e,const size_t& iens)
							    {
							      const map<string,vector<double>> am=getMassList(e);
							      
							      inte.emplace_back(interpolate(mesCombos[iens],mesComposition,e,am,inputAn,totTag));
							      inte.back().plotFf();
							    });
					   
					   //fit
					   for(int iVA=0;iVA<2;iVA++)
					     {
					       cout<<"Fitting "<<VA_tag[iVA]<<endl;
					       
					       dbvec_t pFit(4);
					       boot_fit_t fit;
					       for(size_t i=0;i<4;i++)
						 fit.add_fit_par(pFit[i],combine("p[%zu]",i),0.0,0.1);
					       
					       ensembleLoop(ens,[&](const perens_t& e,const size_t& iens)
								{
								  const dboot_t z=((iVA==0)?Za:Zv)[iM12*nbeta+e.iBeta];
								  cout<<e.dirPath<<" "<<z.ave_err()<<endl;
								  cout<<"Mass: "<<((dboot_t)(inte[iens].E[0]*lat_par[inputAn].ainv[e.iBeta])).ave()<<endl;
								  
								  decKinLoop(e,[&,iens](const size_t iDecKin)
									       {
										 fit.add_point(inte[iens].ff[iVA][iDecKin]*z
											       ,[=,&inte](const vector<double>& p,int iboot)
												{
												  const size_t iBeta=e.iBeta;
												  const double aInv=lat_par[inputAn].ainv[iBeta][iboot];
												  const double a2=1/sqr(aInv);
												  const double x=inte[iens].X[iDecKin][iboot];
												  
												  const double M=inte[iens].E[0][iboot]*aInv;
												  
												  return ansatz(iVA,p,M,a2,x);
												});
									       });
								});
					       
					       fit.fit();
					       
					       ensembleLoop(ens,[&](const perens_t& e,const size_t& iens)
								{
								  const string dirPath=e.dirPath+"/plots/"+inte[iens].mesTag+"_"+totTag;
								  mkdir(dirPath);
								  
								  grace_file_t plot(dirPath+"/ff_"+VA_tag[iVA]+"_"+totTag+"_Fit.xmg");
								  plot.set_no_line();
								  
								  decKinLoop(e,[&,iens](const size_t iDecKin)
									       {
										 plot.write_ave_err(inte[iens].X[iDecKin].ave_err(),inte[iens].ff[iVA][iDecKin].ave_err());
									       });
								  
								  const size_t& iBeta=e.iBeta;
								  const dboot_t aInv=lat_par[inputAn].ainv[iBeta];
								  const double a2=((dboot_t)(1/sqr(aInv))).ave();
								  const double M=((dboot_t)(inte[iens].E[0]*aInv)).ave();
								  const double xMax=inte[iens].xMax()[0];
								  
								  plot.write_polygon([&](const double x) -> dboot_t{return ansatz(iVA,pFit,M,a2,x);},0,xMax);
								  plot.write_polygon([&](const double x) -> dboot_t{return ansatz(iVA,pFit,mPhys,0,x);},0,xMax);
								});
					       
					       cout<<"f: "<<ansatz(iVA,pFit,mPhys,0.0,0.0).ave_err()<<endl;
					     }
					 });
		  }
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
    - | | * loop on all 8 analysis
    - | | | * interpolate to physical meson_corrs
    - | | | * fit all ensembles with a multilinear function
   */
  
  return 0;
}
