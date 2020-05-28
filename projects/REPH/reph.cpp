#include <tranalisi.hpp>

#include <base.hpp>

/// Physical pion mass
const double mPiPhys=0.135;

//               ~~~~~~~~~~~~~~~~~~~~~~~~~~ T/2
//               }                           |
//     --- -> Q0 --- X --- -> Qt ---         |
//    /                             \        |
//   /                               \       |
// P ---------- <- QS -------------V-A (t)

#include <reph.hpp>

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
	    {
	      auto readTints=
		[&](vector<Range>& r)
		{
		  array<size_t,2> o;
		  
		  for(int mM=0;mM<2;mM++)
		    {
		      if(getToken(token,line,pos)!=1) CRASH("Parsing %s",line);
		      if(sscanf(token,"%zu",&o[mM])!=1) CRASH("Parsing %s",token);
		    }
		  
		  r.push_back({o[0],o[1]});
		};
	      
	      readTints(tint2pts);
	      
	      for(int iVA=0;iVA<2;iVA++)
		readTints(tint3pts);
	    }
	}
    }
  
  cout<<"Known ensembles: "<<endl;
  for(auto& it : ensMap)
    cout<<it.first<<endl;
  
  tint2ptsIdx.set_ranges({{"Ens",ensMap.size()},{"Mes",mesMap.size()}});
  tint3ptsIdx.set_ranges({{"Ens",ensMap.size()},{"Mes",mesMap.size()},{"VA",2}});
}

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
  
  const string& firstMeson=get<0>(mesonList[0]);
  if(firstMeson!="Pi+")
    CRASH("First meson in the list is %s, must be Pi+",firstMeson.c_str());
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
      
      //cout<<endl<<" === "<<mesName<<" "<<indMesCombo.descr(iMesCombo)<<" ==="<<endl;
      
      res.emplace_back(ens,mesName,iMs,iMt,eS,eT);
      
      const bool USE_ANALYTIC=true;
      
      res[iMesCombo]
	.fit2pts()
	.prepareKinematics(USE_ANALYTIC)
	.fit3pts()
	.plotFf();
      
      res[iMesCombo].printKin();
      
      //cout<<endl;
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
      //cout<<endl<<"//// Input analysis: "<<inputAn<<" ////"<<endl<<endl;
      
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
auto chirAnsatzA(const TV& p,const double MPi,const double a2,const double x) -> std::remove_reference_t<decltype(p[0])>
{
  return (p[0]+p[2]*a2)/(1-(p[1]+p[3]*a2)*MPi*MPi*(1-x));
}

//! Ansatz fit for fV
template <typename TV>
auto chirAnsatzV(const TV& p,const double MPi,const double a2,const double x) -> std::remove_reference_t<decltype(p[0])>
{
  return (p[0]+p[2]*a2)/(1-(p[1]+p[3]*a2)*MPi*MPi*(1-x));
}

//! Ansatz fit
template <typename TV>
auto chirAnsatz(const size_t iVA,const TV& p,const double M,const double a2,const double x) -> std::remove_reference_t<decltype(p[0])>
{
  if(iVA==0)
    return chirAnsatzV(p,M,a2,x);
  else
    return chirAnsatzA(p,M,a2,x);
}

//! Ansatz fit for fA
template <typename TV>
auto linearAnsatzA(const TV& p,const double MPi,const double a2,const double x) -> std::remove_reference_t<decltype(p[0])>
{
  return p[0]+p[1]*MPi+p[2]*x+p[3]*a2;
}

//! Ansatz fit for fV
template <typename TV>
auto linearAnsatzV(const TV& p,const double MPi,const double a2,const double x) -> std::remove_reference_t<decltype(p[0])>
{
  return p[0]+p[1]*MPi+p[2]*x+p[3]*a2;
}

//! Ansatz fit
template <typename TV>
auto linearAnsatz(const size_t iVA,const TV& p,const double M,const double a2,const double x) -> std::remove_reference_t<decltype(p[0])>
{
  if(iVA==0)
    return linearAnsatzV(p,M,a2,x);
  else
    return linearAnsatzA(p,M,a2,x);
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
  
  const index_t ensInputInd({{"Ens",ens.size()},{"Input",ninput_an}});
  
  dbvec_t MPi(ensInputInd.max());
  
  mesonLoop([&](const meson_t& mesComposition)
	    {
	      const string& mesName=get<0>(mesComposition);
	      // const double& mPhys=get<3>(mesComposition);
	      
	      const index_t indSyst({{"inputAn",ninput_an}});
	      
	      const size_t nFitPars=4;
	      array<dbvec_t,2> storeCh2;
	      array<size_t,2> storeNDof;
	      array<vector<dbvec_t>,2> storePars;
	      for(size_t iVA=0;iVA<2;iVA++)
		{
		  storePars[iVA]=vector<dbvec_t>(indSyst.max(),dbvec_t(nFitPars));
		  storeCh2[iVA]=dbvec_t(indSyst.max());
		}
	      
	      //! Holds ff energy etc for each meson and combination, for each ensemble
	      vector<AllMesCombos> mesCombos;
	      ensembleLoop(ens,[&](const perens_t& e,size_t){mesCombos.push_back(computeAllMesCombos(e,mesComposition));});
	      
	      const bool heavy=(mesName[0]!='P');
	      
	      auto dbvec_ansatz=chirAnsatz<dbvec_t>;
	      auto double_ansatz=chirAnsatz<vector<double>>;
	      
	      if(heavy)
		{
		  dbvec_ansatz=linearAnsatz<dbvec_t>;
		  double_ansatz=linearAnsatz<vector<double>>;
		}
	      
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
							
							if(mesName=="Pi+")
							  {
							    const dboot_t aInv=lat_par[inputAn].ainv[ens[iens].iBeta];
							    MPi[ensInputInd({iens,inputAn})]=inte[iens].E[0]*aInv;
							  }
						      });
				     
				     //fit
				     for(size_t iVA=0;iVA<2;iVA++)
				       {
					 //cout<<"Fitting "<<VA_tag[iVA]<<endl;
					 
					 dbvec_t pFit(nFitPars);
					 const double guess[nFitPars]={0.046,-0.07,-0.1,-0.19};
					 boot_fit_t fit;
					 for(size_t i=0;i<nFitPars;i++)
					   fit.add_fit_par(pFit[i],combine("p[%zu]",i),guess[i],0.05);
					 
					 // fit.fix_par(1);
					 if(not heavy) fit.fix_par(3);
					 
					 ensembleLoop(ens,[&](const perens_t& e,const size_t& iens)
							  {
							    decKinLoop(e,[&,iens](const size_t iDecKin)
									 {
									   fit.add_point(inte[iens].ff[iVA][iDecKin]
											 ,[=,&inte](const vector<double>& p,int iboot)
											  {
											    const size_t iBeta=e.iBeta;
											    const double aInv=lat_par[inputAn].ainv[iBeta][iboot];
											    const double a2=1/sqr(aInv);
											    const double x=inte[iens].X[iDecKin][iboot];
											    
											    const double M=MPi[ensInputInd({iens,inputAn})][iboot];
											    
											    return double_ansatz(iVA,p,M,a2,x);
											  });
									 });
							  });
					 
					 auto status=fit.fit();
					 
					 //cout<<"Fit pars\n"<<pFit.ave_err()<<endl;
					 const double xMin=1e-3;
					 double maxXMax=0;
					 
					 const string fitDirPath="plots/"+mesName+"/"+to_string(inputAn);
					 mkdir(fitDirPath);
					 grace_file_t fitPlot(fitDirPath+"/ff_"+VA_tag[iVA]+"_Fit.xmg");
					 grace_file_t sliceA2Plot(fitDirPath+"/ff_"+VA_tag[iVA]+"_funA2_Fit.xmg");
					 grace_file_t sliceMpiPlot(fitDirPath+"/ff_"+VA_tag[iVA]+"_funMpi_Fit.xmg");
					 for(auto& p : {&fitPlot,&sliceA2Plot}) p->set_no_line();
					 
					 vector<grace::color_t> colors{grace::color_t::RED,grace::color_t::ORANGE,grace::color_t::GREEN4};
					 
					 map<string,double> targetXlist={{"Pi+",0.75},{"K+",0.75},{"D+",0.3},{"Ds",0.3}};
					 const double& targetX=targetXlist[mesName];
					 ensembleLoop(ens,[&](const perens_t& e,const size_t& iens)
							  {
							    const size_t& iBeta=e.iBeta;
							    const string ensDirPath=e.dirPath+"/plots/"+inte[iens].mesTag;
							    
							    // const dboot_t z=((iVA==0)?Za:Zv)[iM12*nbeta+e.iBeta];
							    
							    grace_file_t ensPlot(ensDirPath+"/ff_"+VA_tag[iVA]+"_Fit.xmg");
							    for(auto& p : {&fitPlot,&sliceA2Plot,&sliceMpiPlot})
							      {
								p->new_data_set();
								p->set_comment(ensDirPath);
							      }
							    
							    for(auto& p : {&ensPlot,&fitPlot,&sliceA2Plot,&sliceMpiPlot})
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
							    
							    const double M=MPi[ensInputInd({iens,inputAn})].ave();
							    const double xMax=inte[iens].xMax()[0];
							    maxXMax=max(maxXMax,xMax);
							    
							    ensPlot.write_polygon([&](const double x) -> dboot_t{return dbvec_ansatz(iVA,pFit,M,a2,x);},xMin,xMax);
							    fitPlot.write_line([&](const double x){return dbvec_ansatz(iVA,pFit,M,a2,x).ave();},xMin,xMax,colors[iBeta]);
							    
							    sliceMpiPlot.write_ave_err(M,closestY.ave_err());
							    sliceA2Plot.write_ave_err(a2+M/get<3>(mesComposition)/100,closestY.ave_err());
							  });
					 
					 fitPlot.write_polygon([&](const double x) -> dboot_t{return dbvec_ansatz(iVA,pFit,mPiPhys,0,x);},xMin,maxXMax,grace::color_t::VIOLET);
					 sliceA2Plot.write_polygon([&](const double x) -> dboot_t{return dbvec_ansatz(iVA,pFit,mPiPhys,x,targetX);},1e-3,0.3,grace::color_t::VIOLET);
					 
					 sliceMpiPlot.write_ave_err(mPiPhys,dbvec_ansatz(iVA,pFit,mPiPhys,0,targetX).ave_err());
					 
					 //f as a function of mpi
					 for(size_t ib=0;ib<3;ib++)
					   {
					     const dboot_t aInv=lat_par[inputAn].ainv[ib];
					     const double a2=((dboot_t)(1/sqr(aInv))).ave();
					     sliceMpiPlot.write_line([&](const double x){return dbvec_ansatz(iVA,pFit,x,a2,targetX).ave();},1e-3,0.5,colors[ib]);
					   }
					 //cont
					 sliceMpiPlot.write_polygon([&](const double x) -> dboot_t{return dbvec_ansatz(iVA,pFit,x,0,targetX);},1e-3,0.5,grace::color_t::VIOLET);
					 
					 /// Index of the systematic study
					 const size_t iSyst=indSyst(vector<size_t>{inputAn});
					 
					 //! Store for future uses
					 storePars[iVA][iSyst]=pFit;
					 storeCh2[iVA][iSyst]=get<0>(status);
					 storeNDof[iVA]=get<1>(status);
				       }
				   });
	      
	      /////////////////////////////////////////////////////////////////
	      
	      // Total systematics analysis
	      
	      //! Directory
	      const string dirPath="plots/"+mesName;
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
	      		temp[iSyst]=dbvec_ansatz(iVA,storePars[iVA][iSyst],mPiPhys,0.0,x[i]);
		      
	      	      const syst_t s=perform_analysis(temp,indSyst);
		      
	      	      y[i].ave()=s.ave;
	      	      y[i].err()=s.tot;
	      	    }
		  
	      	  plot.write_polygon(x,y);
		  
		  const syst_t ch2=perform_analysis(storeCh2[iVA],indSyst);
		  
		  cout<<mesName<<", ff"<<VA_tag[iVA]<<"(0): "<<smart_print(y[0])<<" , ch2: "<<smart_print({ch2.ave,ch2.tot})<<"/"<<storeNDof[iVA]<<endl;
		}
	      
	      // // plot.write_polygon([&](const double x) -> dboot_t{return dbvec_ansatz(iVA,pFit,mPhys,0,x);},0,xMax);
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
