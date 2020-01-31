#include <tranalisi.hpp>

#include <base.hpp>

void read3ptsTint()
{
  raw_file_t input("tints.txt","r");
  input.expect("Ens");
  
  auto getToken=[](char *token,char *line,int& pos)
		{
		  // Discard blank chars
		  while(line[pos]==' ' or line[pos]=='\t')
		    pos++;
		  
		  size_t rc=sscanf(line+pos,"%s",token);
		  
		  // Advance
		  if(rc==1)
		    pos+=strlen(token);
		  
		  return rc;
		};
  
  char line[1024];
  input.get_line(line);
  
  int rc,pos=0;
  size_t iMes=0;
  do
    {
      // Read
      char token[1024];
      rc=getToken(token,line,pos);
      
      // Advance
      if(rc==1)
	mesMap[token]=iMes++;
    }
  while(rc==1);
  
  cout<<"Known mesons: "<<endl;
  for(auto& it : mesMap)
    cout<<it.first<<endl;
  
  size_t iEns=0;
  while(input.get_line(line))
    {
      int rc;
      pos=0;
      
      char token[1024];
      rc=getToken(token,line,pos);
      if(rc==1)
	{
	  ensMap[token]=iEns++;
	  
	  for(size_t i=0;i<mesMap.size();i++)
	    for(int iVA=0;iVA<2;iVA++)
	      {
		size_t o[2];
		for(int mM=0;mM<2;mM++)
		  {
		    if(getToken(token,line,pos)!=1) CRASH("Parsing %s",line);
		    if(sscanf(token,"%zu",&o[mM])!=1) CRASH("Parsing %s",token);
		  }
		tint3pts.push_back({o[0],o[1]});
	      }
	}
    }
  
  cout<<"Known ensembles: "<<endl;
  for(auto& it : ensMap)
    cout<<it.first<<endl;
  
  tint3ptsIdx.set_ranges({{"Ens",ensMap.size()},{"Mes",mesMap.size()},{"VA",2}});
}

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
      
      if(mesMap.find(name)==mesMap.end())
	CRASH("Unable to find meson %s",name.c_str());
      
      mesonList.push_back({name,quark1,quark2,mass});
    }
}

//! Reads the physics file
void readPhysics(const string& path)
{
  //! Input file
  raw_file_t input(path,"r");
  
  oldNormalization=input.read<int>("OldNormalization");
  
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
      
      const bool USE_ANALYTIC=true;
      
      res[iMesCombo]
	.fit2pts("selfChosenTint")
	.prepareKinematics(USE_ANALYTIC)
	.fit3pts()
	.plotFf();
      
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
void ensembleLoop(const vector<perens_t>& ens,const F& f,const bool& verbose=false)
{
  const size_t& nens=ens.size();
  
  for(size_t iens=0;iens<nens;iens++)
    {
      if(verbose) cout<<endl<<"**** "<<ens[iens].dirPath<<" ****"<<endl<<endl;
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
  //return p[0]+p[1]*M*M+p[2]*M*M*x+p[3]*a2;
  return (p[0]+p[2]*a2)/(1-(p[1]+p[3]*a2)*M*M*(1-x));
}

//! Ansatz fit for fV
template <typename TV>
auto ansatzV(const TV& p,const double M,const double a2,const double x) -> std::remove_reference_t<decltype(p[0])>
{
  return p[0]+p[1]*M*M+p[2]*M*M*x+p[3]*a2;
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
  read3ptsTint();
  
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
  
  mesonLoop([&](const meson_t& mesComposition)
	    {
	      const double& mPhys=get<3>(mesComposition);
	      
	      const index_t indSyst({{"inputAn",ninput_an}});
	      
	      const size_t nFitPars=4;
	      array<vector<dbvec_t>,2> storePars;
	      for(size_t iVA=0;iVA<2;iVA++)
		storePars[iVA]=vector<dbvec_t>(indSyst.max(),dbvec_t(nFitPars));
	      
	      //! Holds ff energy etc for each meson and combination, for each ensemble
	      vector<AllMesCombos> mesCombos;
	      ensembleLoop(ens,[&](const perens_t& e,size_t){mesCombos.push_back(computeAllMesCombos(e,mesComposition));});
	      
	      // Loop over analysis
	      ultimateAnalysisLoop([&](const size_t& inputAn)
				   {
				     /// Method 1 or 2
				     //const size_t iM12=inputAn/4;
				     
				     //! Interpolated data for each ensemble
				     vector<permes_t<dbvec_t>> inte;
				     
				     ensembleLoop(ens,[&](const perens_t& e,const size_t& iens)
						      {
							const map<string,vector<double>> am=getMassList(e);
							
							inte.emplace_back(interpolate(mesCombos[iens],mesComposition,e,am,inputAn));
							
							inte.back().correctFf(mesComposition,inputAn);
							inte.back().plotFf();
						      });
				     
				     //fit
				     for(size_t iVA=0;iVA<2;iVA++)
				       {
					 cout<<"Fitting "<<VA_tag[iVA]<<endl;
					 
					 dbvec_t pFit(nFitPars);
					 boot_fit_t fit;
					 for(size_t i=0;i<nFitPars;i++)
					   fit.add_fit_par(pFit[i],combine("p[%zu]",i),0.0,0.1);
					 
					 fit.fix_par(1);
					 
					 ensembleLoop(ens,[&](const perens_t& e,const size_t& iens)
							  {
							    // const dboot_t z=((iVA==0)?Za:Zv)[iM12*nbeta+e.iBeta];
							    // cout<<e.dirPath<<" "<<z.ave_err()<<endl;
							    // cout<<"Mass: "<<((dboot_t)(inte[iens].E[0]*lat_par[inputAn].ainv[e.iBeta])).ave()<<endl;
							    
							    decKinLoop(e,[&,iens](const size_t iDecKin)
									 {
									   fit.add_point(inte[iens].ff[iVA][iDecKin]// *z
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
					 cout<<"Fit pars\n"<<pFit.ave_err()<<endl;
					 const double xMin=1e-3;
					 double maxXMax=0;
					 
					 const string fitDirPath="plots/"+get<0>(mesComposition)+"/"+to_string(inputAn);
					 grace_file_t fitPlot(fitDirPath+"/ff_"+VA_tag[iVA]+"_Fit.xmg");
					 grace_file_t sliceA2Plot(fitDirPath+"/ff_"+VA_tag[iVA]+"_funA2_Fit.xmg");
					 for(auto& p : {&fitPlot,&sliceA2Plot}) p->set_no_line();
					 
					 vector<grace::color_t> colors{grace::color_t::RED,grace::color_t::BLACK,grace::color_t::GREEN4};
					 
					 const double targetX=0.75;
					 ensembleLoop(ens,[&](const perens_t& e,const size_t& iens)
							  {
							    const size_t& iBeta=e.iBeta;
							    const string ensDirPath=e.dirPath+"/plots/"+inte[iens].mesTag;
							    
							    // const dboot_t z=((iVA==0)?Za:Zv)[iM12*nbeta+e.iBeta];
							    
							    grace_file_t ensPlot(ensDirPath+"/ff_"+VA_tag[iVA]+"_Fit.xmg");
							    for(auto& p : {&fitPlot,&sliceA2Plot})
							      {
								p->new_data_set();
								p->set_comment(ensDirPath);
							      }
							    
							    for(auto& p : {&ensPlot,&fitPlot,&sliceA2Plot})
							      {
								p->set_no_line();
								p->set_all_colors(colors[iBeta]);
							      }
							    
							    double closestX=0;
							    dboot_t closestY;
							    
							    decKinLoop(e,[&,iens](const size_t iDecKin)
									 {
									   const ave_err_t X=inte[iens].X[iDecKin].ave_err();
									   const dboot_t Y=inte[iens].ff[iVA][iDecKin];
									   
									   if(closestX==0 or fabs(closestX-targetX)>fabs(X.ave()-targetX))
									     {
									       closestX=X.ave();
									       closestY=Y;
									     }
									   
									   for(auto& p : {&ensPlot,&fitPlot})
									     p->write_ave_err(X,Y.ave_err());
									 });
							    
							    const dboot_t aInv=lat_par[inputAn].ainv[iBeta];
							    const double a2=((dboot_t)(1/sqr(aInv))).ave();
							    const double M=((dboot_t)(inte[iens].E[0]*aInv)).ave();
							    const double xMax=inte[iens].xMax()[0];
							    maxXMax=max(maxXMax,xMax);
							    
							    ensPlot.write_polygon([&](const double x) -> dboot_t{return ansatz(iVA,pFit,M,a2,x);},xMin,xMax);
							    fitPlot.write_line([&](const double x){return ansatz(iVA,pFit,M,a2,x).ave();},xMin,xMax,colors[iBeta]);
							    
							    sliceA2Plot.write_ave_err(a2+M/get<3>(mesComposition)/100,closestY.ave_err());
							  });
					 
					 fitPlot.write_polygon([&](const double x) -> dboot_t{return ansatz(iVA,pFit,mPhys,0,x);},xMin,maxXMax,grace::color_t::VIOLET);
					 sliceA2Plot.write_polygon([&](const double x) -> dboot_t{return ansatz(iVA,pFit,mPhys,x,targetX);},1e-3,0.3,grace::color_t::VIOLET);
					 
					 const dboot_t cph=ansatz(iVA,pFit,mPhys,0.0,xMin);
					 cout<<"f: "<<cph.ave_err()<<endl;
					 
					 //! Store for future uses
					 storePars[iVA][indSyst(vector<size_t>{inputAn})]=pFit;
				       }
				   });
	      
	      /////////////////////////////////////////////////////////////////
	      
	      // Total systematics analysis
	      
	      //! Directory
	      const string dirPath="plots/"+get<0>(mesComposition);
	      mkdir(dirPath);
	      
	      for(size_t iVA=0;iVA<2;iVA++)
		{
		  grace_file_t plot(dirPath+"/ff_"+VA_tag[iVA]+".xmg");
		  plot.set_no_line();
		  
		  const double xMax=1.0;
		  const vector<double> x=vector_up_to(xMax,0.0,0.05);
		  vec_ave_err_t y(x.size());
		  
		  for(size_t i=0;i<x.size();i++)
		    {
		      dbvec_t temp(indSyst.max());
		      for(size_t iSyst=0;iSyst<indSyst.max();iSyst++)
			temp[iSyst]=ansatz(iVA,storePars[iVA][iSyst],mPhys,0.0,x[i]);
		      
		      const syst_t s=perform_analysis(temp,indSyst);
		      
		      y[i].ave()=s.ave;
		      y[i].err()=s.tot;
		    }
		  
		  plot.write_polygon(x,y);
		}
	      
	      // plot.write_polygon([&](const double x) -> dboot_t{return ansatz(iVA,pFit,mPhys,0,x);},0,xMax);
	      //perform_analysis(fTest,iTest,"ciccio");
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
