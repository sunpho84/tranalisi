#include <tranalisi.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int twist;
int nF;
int nHits;
array<int,4> dims;
int T;
int nMass;
int nMoms;

double beta;
double kSea;
double muSea;
double cSW;

vector<double> kappa;
vector<double> mu;
vector<array<double,4>> moms;
int nobs;

bool saveText=false;
bool saveJacks=false;

void convert_file(const char *corr,const bool is2pt,const size_t nGamma)
{
  const string input_path=combine("data/%s_conf.realph.dat",corr);
  const string text_path=combine("converted/%s_conf.realph.txt",corr);
  const string jack_path=combine("jacks/%s",corr);
  
  if(saveText) cout<<input_path<<" -> "<<text_path<<endl;
  if(saveJacks) cout<<input_path<<" -> "<<jack_path<<endl;
  
  raw_file_t fin(input_path,"r");
  
  for(auto &t : {&twist,&nF,&nHits})
    fin.bin_read(*t);
  printf(" - twisted: %s\n",(twist ? "yes" : "no"));
  printf(" - flavors: %d\n",nF);
  printf(" - hits: %d\n",nHits);
  
  fin.bin_read(dims);
  printf(" - volume: %dx%dx%dx%d\n",dims[0],dims[1],dims[2],dims[3]);
  T=dims[0];
  
  for(auto &t : {&nMass,&nMoms})
    fin.bin_read(*t);
  printf(" - masses: %d\n",nMass);
  printf(" - momenta: %d\n",nMoms);
  
  kappa.resize(nMass);
  mu.resize(nMass);
  moms.resize(nMoms);
  
  for(auto &t : {&beta,&kSea,&muSea,&cSW})
    fin.bin_read(*t);
  
  printf(" - beta: %lf\n",beta);
  printf(" - kappa: %lf\n",kSea);
  printf(" - mu: %lf\n",muSea);
  printf(" - csw: %lf\n",cSW);
  
  for(auto &t : {&kappa,&mu})
    fin.bin_read(*t);
  
  for(size_t n=0;n<(size_t)nMass;n++)
    printf(" - mass#%lu: (%f,%f)\n",n,kappa[n],mu[n]);
  
  fin.bin_read(moms);
  
  for(size_t n=0;n<(size_t)nMoms;n++)
    printf(" - momenta#%lu: (%f,%f,%f,%f)\n",n,moms[n][0],moms[n][1],moms[n][2],moms[n][3]);
  
  const size_t nMom0=(is2pt?1:nMoms);
  const size_t nPol=(is2pt?1:2);
  
  index_t indIn({{"iks",nMass},
		 {"ikt",nMass},
		 {"moms",nMoms},
		 {"momt",nMoms},
		 {"mom0",nMom0},
		 {"t",T},
		 {"pol",nPol},
		 {"gamma",nGamma},
		 {"reim",2}});
  
  index_t indOut({{"iks",nMass},
		  {"ikt",nMass},
		  {"moms",nMoms},
		  {"momt",nMoms},
		  {"mom0",nMom0},
		  {"pol",nPol},
		  {"gamma",nGamma},
		  {"reim",2},
		  {"t",T}});
  
  const int nDoublesExpected=indIn.max();
  const int nDoubles=fin.bin_read<int>();
  if(nDoubles!=nDoublesExpected)
    CRASH("Expecting %d doubles, obtained %d",nDoublesExpected,nDoubles);
  
  const size_t corrPerConfSize=nDoubles*sizeof(double);
  
  vector<double> temp(indIn.max());
  
  const size_t restFileSize=fin.size()-fin.get_pos();
  const int nConfs=restFileSize/corrPerConfSize;
  const size_t tagsSize=restFileSize-nConfs*corrPerConfSize;
  const int nConfsBis=tagsSize/sizeof(int);
  
  if(nConfs!=nConfsBis)
    CRASH("nConfs determined from rest of the file size %d and tags %d do not match",nConfs,nConfsBis);
  
  raw_file_t fout_text;
  if(saveText)
    fout_text.open(text_path,"w");
  
  djvec_t *out=nullptr;
  if(saveJacks)
    {
      out=new djvec_t(indOut.max());
      *out=0.0;
    }
  
  const int clustSize=saveJacks?(nConfs/njacks):0;
  for(size_t iConf=0;iConf<clustSize*njacks;iConf++)
    {
      int iClust=0;
      if(saveJacks)
	iClust=iConf/clustSize;
      
      int tag;
      fin.bin_read(tag);
      cout<<tag<<endl;
      
      fin.bin_read(temp);
      
      for(size_t iks=0;iks<(size_t)nMass;iks++)
	for(size_t ikt=0;ikt<(size_t)nMass;ikt++)
	  for(size_t imoms=0;imoms<(size_t)nMoms;imoms++)
	    for(size_t imomt=0;imomt<(size_t)nMoms;imomt++)
	      for(size_t imom0=0;imom0<(size_t)nMom0;imom0++)
		for(size_t iPol=0;iPol<nPol;iPol++)
		  {
		    if(saveText)
		      {
			fout_text.printf("\n # mass_s %lg , mass_t %lg , mom_s ( ",mu[iks],mu[ikt]);
			
			for(int j=1;j<4;j++)
			  fout_text.printf("%lg ",moms[imoms][j]);
			fout_text.printf(") , mom_t ( ");
			for(int j=1;j<4;j++)
			  fout_text.printf("%lg ",moms[imomt][j]);
			
			if(not is2pt)
			  {
			    fout_text.printf(") , mom_0 ( ");
			    for(int j=1;j<4;j++)
			      fout_text.printf("%lg ",moms[imom0][j]);
			    
			    fout_text.printf(") pol %zu",iPol);
			  }
			fout_text.printf("\n");
		      }
		    
		    for(size_t iGamma=0;iGamma<nGamma;iGamma++)
		      {
			if(saveText) fout_text.printf("\n # Gamma %lu\n\n",iGamma);
			
			for(size_t t=0;t<(size_t)T;t++)
			  for(size_t ri=0;ri<2;ri++)
			    {
			      const int iIn=indIn({iks,ikt,imoms,imomt,imom0,t,iPol,iGamma,ri});
			      const int iOut=indOut({iks,ikt,imoms,imomt,imom0,iPol,iGamma,ri,t});
			      
			      if(saveJacks)
				(*out)[iOut][iClust]+=temp[iIn];
			      
			      if(saveText)
				{
				  fout_text.printf("%+.16lg",temp[iIn]);
				  fout_text.printf(ri==RE?"\t":"\n");
				}
			    }
		      }
		  }
    }
  
  if(saveJacks)
    {
      out->clusterize(clustSize);
      out->bin_write(jack_path);
      
      delete out;
    }
  
  fin.close();
}

int main(int narg,char **arg)
{
  //parse opts
  int c;
  while((c=getopt(narg,arg,"j:th"))!= -1)
    switch (c)
      {
      case 't': saveText=true;break;
      case 'j':
	set_njacks(to_int(optarg));
	saveJacks=true;
	break;
      case 'h':
	cout<<"Use "<<arg[0]<<" -jnjacks -t -h"<<endl;
	exit(0);
	break;
      default: CRASH("Unknown option -%c, try -h for help ",optopt);
      }
  
  convert_file("oPPo-ss",true,1);
  convert_file("oAmuPo-ss",true,4);
  convert_file("oPGPo-gs",false,1);
  convert_file("oAmuGPo-gs",false,4);
  convert_file("oVmuGPo-gs",false,4);
  
  if(saveJacks)
    {
      raw_file_t input_out("jacks/input.txt","w");
      
      input_out.printf("L %d\n",dims[1]);
      input_out.printf("T %d\n",dims[0]);
      input_out.printf("NMass %d\n",nMass);
      for(auto &m : mu)
	input_out.printf(" %lg\n",m);
      
      input_out.printf("NMoms %d\n",nMoms);
      for(auto &m : moms)
	{
	  for(int i=1;i<4;i++)
	    input_out.printf("%lg ",m[i]);
	  input_out.printf("\n");
	}
    }
  
  return 0;
}
