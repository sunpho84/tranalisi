#include <tranalisi.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

void check_file(const char *corr,const bool is2pt,const size_t nGamma)
{
  const string input_path=combine("data/%s_conf.realph.dat",corr);
  const string output_path=combine("converted/%s_conf.realph.dat",corr);
  
  cout<<input_path<<" -> "<<output_path<<endl;
  
  if(file_exists(output_path))
     return;
  
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
  
  index_t ind({{"iks",nMass},
		{"ikt",nMass},
		{"moms",nMoms},
		{"momt",nMoms},
		{"mom0",nMom0},
		{"t",T},
		{"pol",nPol},
		{"gamma",nGamma},
		{"reim",2}});
  
  // for(iks=0;iks<nk;iks++)
  //    for(ikt=0;ikt<nk;ikt++)
  //    for(imoms=0;imoms<nmoms;imoms++) // spectator momentum (associated for the iks mass)
  //    for(imomt=0;imomt<nmoms;imomt++) // momentum after photon insertion
  //    for(imom0=0;imom0<nmoms;imom0++) // momentum before photon insertion
  //    for(t=0;t<T;t++)
  //    for(ipol=0;ipol<2;ipol++)
  //    for(igamma=0;igamma<Ngamma;igamma++)
  //    for(ire=0;ire<2;ire++)
  
  const int nDoublesExpected=ind.max();
  const int nDoubles=fin.bin_read<int>();
  if(nDoubles!=nDoublesExpected)
    CRASH("Expecting %d doubles, obtained %d",nDoublesExpected,nDoubles);
  
  const size_t corrPerConfSize=nDoubles*sizeof(double);
  
  vector<double> temp(ind.max());
  
  const size_t restFileSize=fin.size()-fin.get_pos();
  const int nConfs=restFileSize/corrPerConfSize;
  const size_t tagsSize=restFileSize-nConfs*corrPerConfSize;
  const int nConfsBis=tagsSize/sizeof(int);
  
  if(nConfs!=nConfsBis)
    CRASH("nConfs determined from rest of the file size %d and tags %d do not match",nConfs,nConfsBis);
  
  raw_file_t fout(output_path,"w");
  
  for(int iConf=0;iConf<nConfs;iConf++)
    {
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
		    fout.printf("\n # mass_s %lg , mass_t %lg , mom_s ( ",mu[iks],mu[ikt]);
		    for(int j=1;j<4;j++)
		      fout.printf("%lg ",moms[imoms][j]);
		    fout.printf(") , mom_t ( ");
		    for(int j=1;j<4;j++)
		      fout.printf("%lg ",moms[imomt][j]);
		    
		    if(not is2pt)
		      {
			fout.printf(") , mom_0 ( ");
			for(int j=1;j<4;j++)
			  fout.printf("%lg ",moms[imom0][j]);
			
			fout.printf(") pol %zu",iPol);
		      }
		    fout.printf("\n");
		    
		    for(size_t iGamma=0;iGamma<nGamma;iGamma++)
		      {
			fout.printf("\n # Gamma %lu\n",iGamma);
			
			for(size_t t=0;t<(size_t)T;t++)
			  for(size_t ri=0;ri<2;ri++)
			    {
			      const int i=ind({iks,ikt,imoms,imomt,imom0,t,iPol,iGamma,ri});
			      fout.printf("%+.16lg",temp[i]);
			      fout.printf(ri==RE?"\t":"\n");
			    }
		      }
		  }
    }
  
  // size = sizeof(double)*lh.nobs + sizeof(int);
  // cpos = ftell(fp);
  // fseek(fp, 0, SEEK_END);
  // epos = ftell(fp);
  // epos = epos-cpos;
  // nc = epos/size;
  
  // if(epos % size)
  //   {
  //     printf(" - ERROR: file size wrong\n\n");
  //     fclose(fp);
  //     exit(1);
  //   }
  
  // fseek(fp, -size, SEEK_END);
  // fread(&last, sizeof(int), 1, fp);
  // fseek(fp, cpos, SEEK_SET);
  // fread(&first, sizeof(int), 1, fp);
  
  // if(init == 0)
  //   {
  //     step = (last-first+1)/nc;
  //     head_size = cpos;
  //   }
  
  // cpos = (last-first+1)/nc;
  // printf(" - configurations: %d\n", nc);
  // printf(" - range: %d - %d\n", first, last);
  // printf(" - spacing: %d\n", cpos);
  
  // if(cpos != step)
  //   {
  //     printf(" - ERROR: configuration spacing is different\n\n");
  // 		fclose(fp);
  // 		exit(1);
  //   }
  
  // if(init == 0)
  //   {
  //     c0 = first;
  //     c1 = last;
  //   }
  // else
  //   {
  //     if(c1+step != first)
  // 	{
  // 	  printf(" - ERROR: range does not continue previous file\n\n");
  // 	  fclose(fp);
  // 	  exit(1);
  // 	}
  //     c1 = last;
  //   }
  
  // init = 1;
  
  fin.close();
}

int main(int narg,char **arg)
{
  check_file("oPPo-ss",true,1);
  check_file("oAmuPo-ss",true,4);
  check_file("oPGPo-gs",false,1);
  check_file("oAmuGPo-gs",false,4);
  check_file("oVmuGPo-gs",false,4);
  
  return 0;
}
