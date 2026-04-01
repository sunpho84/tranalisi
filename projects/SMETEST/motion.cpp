#include <tranalisi.hpp>

double kappa;
vector<string> mesons;
vector<double> thetas;
vector<size_t> smeLevs;
const size_t T=64,L=32;

int main()
{
  input_file_t parsFile("pars.txt");
  set_njacks(parsFile.read<size_t>("NJacks"));
  kappa=parsFile.read<double>("Kappa");
  
  const size_t nMeson=parsFile.read<size_t>("Mesons");
  mesons.resize(nMeson);
  for(size_t iMes=0;iMes<nMeson;iMes++)
    mesons[iMes]=parsFile.read<string>();
  
  const size_t nTheta=parsFile.read<size_t>("Thetas");
  thetas.resize(nTheta);
  for(size_t iTheta=0;iTheta<nTheta;iTheta++)
    thetas[iTheta]=parsFile.read<double>();
  
  const size_t nSmeLev=parsFile.read<size_t>("SmeLevs");
  smeLevs.resize(nSmeLev);
  for(size_t iSme=0;iSme<nSmeLev;iSme++)
    smeLevs[iSme]=parsFile.read<size_t>();
  
  const index_t id({{"",nMeson},{"",nTheta},{"",nTheta},{"",nSmeLev},{"",nSmeLev},{"T",T}});
  
  auto getCorr=
    [&id,rawData=read_djvec("data/mes_contr_2PT_P5P5",id.max())](const size_t iMes,
								 const size_t iTheta1,
								 const size_t iTheta2,
								 const size_t iSme1,
								 const size_t iSme2)
    {
      djvec_t res(T);
      for(size_t t=0;t<T;t++)
	res[t]=rawData[id({iMes,iTheta1,iTheta2,iSme1,iSme2,t})];
      
      return res.symmetrized();
    };
  
  // Reorder the smearing combo in order by the squared radius
  vector<array<size_t,3>> rad2Pairs;
  for(size_t iSme1=0;iSme1<nSmeLev;iSme1++)
    for(size_t iSme2=0;iSme2<nSmeLev;iSme2++)
      rad2Pairs.push_back({smeLevs[iSme1]+smeLevs[iSme2],iSme1,iSme2});
  sort(rad2Pairs.begin(),rad2Pairs.end());
  
  for(size_t iMes=0;iMes<nMeson;iMes++)
    for(size_t iTheta1=0;iTheta1<nTheta;iTheta1++)
      for(size_t iTheta2=0;iTheta2<nTheta;iTheta2++)
	{
	  const string suff=mesons[iMes]+"_th_"+to_string(iTheta1)+"_"+to_string(iTheta2);
	  
	  auto m=[&](const size_t iSme1,const size_t iSme2)
	  {
	    return effective_mass(getCorr(iMes,iTheta1,iTheta2,iSme1,iSme2));
	  };
	  
	  // Computes the effective mass at fixed time
	  vector<djack_t> comp(rad2Pairs.size());
	  for(size_t i=0;i<rad2Pairs.size();i++)
	    {
	      const auto& [rad2,iSme1,iSme2]=rad2Pairs[i];
	      
	      const djvec_t mEff=m(iSme1,iSme2);
	      comp[i]=mEff[6];
	    }
	  
	  vector<double> fitMinX;
	  vector<djack_t> _fitMinY;
	  for(size_t i=0;i<rad2Pairs.size();i++)
	    if(not isnan(comp[i].err()))
	      {
		const double x=(i==0)?rad2Pairs[1][0]/2.0:rad2Pairs[i][0];
		fitMinX.push_back(log(x));
		_fitMinY.push_back(comp[i]);
		// cout<<fitMinX.back()<<" "<<_fitMinY.back().ave_err()<<endl;
	      }
	  djvec_t fitMinY(_fitMinY.size());
	  for(size_t i=0;i<_fitMinY.size();i++)
	    fitMinY[i]=_fitMinY[i];
	  
	  const auto fitmMinPars=poly_fit(fitMinX,fitMinY,2,fitMinX[2],fitMinX.back()/1.0,"plots/optRad_"+suff+".xmg");
	  const djack_t xMin=exp(-fitmMinPars[1]/(2*fitmMinPars[2]));
	  
	  int iNewChampion=-1;
	  int significativity=3;
	  do
	    {
	      for(size_t i=1;i<rad2Pairs.size();i++)
		{
		  const djack_t dist=log(rad2Pairs[i][0])-log(xMin);
		  
		  if(dist.significativity()<significativity)
		    if(iNewChampion==-1 or comp[i].err()<comp[iNewChampion].err())
		      iNewChampion=i;
		}
	      significativity++;
	    }
	  while(iNewChampion==-1);

	  // map<size_t,vector<pair<double,djack_t>>> fracPerRad2;
	  // for(size_t i=0;i<rad2Pairs.size();i++)
	  //   {
	  //     const auto& [rad2,iSme1,iSme2]=rad2Pairs[i];
	      
	  //     //mEff.ave_err().write("plots/effMass_"+mesons[iMes]+"_th_"+to_string(iTheta1)+"_"+to_string(iTheta2)+"_iSme1_"+to_string(iSme1)+"_"+to_string(iSme2)+".xmg");
	  //     const double f=smeLevs[iSme1]/double(smeLevs[iSme1]+smeLevs[iSme2]);
	  //     fracPerRad2[rad2].emplace_back(f,comp);
	      
	  //     if(rad2==0 or comp.ave()<champion.ave() or (rad2==rad2Pairs[iChampion][0] and comp.err()<champion.err()))
	  // 	{
	  // 	  // cout<<"Old rad2: "<<rad2Pairs[iChampion][0]<<", "<<champion.ave_err()<<endl;
		  
	  // 	  iChampion=i;
	  // 	  cout<<"New rad2: "<<rad2Pairs[iChampion][0]<<" "<<rad2Pairs[iChampion][1]<<" "<<rad2Pairs[iChampion][2]<<", "<<champion.ave_err()<<endl;
	  // 	}
	  //     plot.write_ave_err(sqrt(rad2),mEff[6].ave_err());
	  //   }
	  
	  grace_file_t file("plots/effMass_"+suff+".xmg");
	  file.write_vec_ave_err(m(0,0).ave_err());
	  const auto [rad2C,iSmeC1,iSmeC2]=rad2Pairs[iNewChampion];
	  file.write_vec_ave_err(m(iSmeC1,iSmeC2).ave_err());
	  
	  cout<<mesons[iMes]<<" th1="<<iTheta1<<" th2="<<iTheta2<<" minfit: "<<xMin.ave_err()<<" closest: "<<rad2C<<" "<<iSmeC1<<" "<<iSmeC2<<endl;
	  
	  map<size_t,vector<pair<double,djack_t>>> fracPerRad2;
	  for(size_t i=0;i<rad2Pairs.size();i++)
	    {
	      const auto& [rad2,iSme1,iSme2]=rad2Pairs[i];
	      
	      //mEff.ave_err().write("plots/effMass_"+mesons[iMes]+"_th_"+to_string(iTheta1)+"_"+to_string(iTheta2)+"_iSme1_"+to_string(iSme1)+"_"+to_string(iSme2)+".xmg");
	      const double f=smeLevs[iSme1]/double(smeLevs[iSme1]+smeLevs[iSme2]);
	      fracPerRad2[rad2].emplace_back(f,comp[i]);
	    }
	  
	  for(auto& f : fracPerRad2)
	    if(f.second.size()>2)
	      {
		grace_file_t plot("plots/fracs_"+suff+"_rad2_"+to_string(f.first)+".xmg");
		sort(f.second.begin(),f.second.end());
		for(const auto& fi : f.second)
		  plot.write_ave_err(fi.first,fi.second.ave_err());
	      }
	}
  
  const double m=0.896;
  
  auto q2=[](const double& mB)
  {
    const double th=3.9254;
    const double m=0.896;
    const double pi=th*3.14159/32;
    const double e=latt_en(m,pi);
    const double q2=sqr(mB-e)-sqr(pi);
    
    return q2;
  };

  cout<<"mB: "<<Brent_solve(q2,m,100)<<endl;
  grace_file_t ooo("/tmp/ooo.xmg");
  ooo.write_line(q2,m,2*m);
  
  return 0;
}
