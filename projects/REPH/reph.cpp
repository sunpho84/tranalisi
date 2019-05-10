#include <tranalisi.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct
{
  int twist;
  int nf;
  int nhits;
  int dim_t;
  int dim_x;
  int dim_y;
  int dim_z;
  int nk;
  int nm;
  double beta;
  double ksea;
  double musea;
  double csw;
} __attribute__((packed)) par_t;

typedef struct
{
  double p0;
  double p1;
  double p2;
  double p3;
} __attribute__((packed)) mom_t;

typedef struct
{
  par_t pars;
  double kappa[64];
  double mu[64];
  mom_t mom[64];
  int nobs;
} file_head;

static int init;
static file_head gh;
static int c0,c1;
static int head_size;
static int step;

void check_file(const char *path)
{
  // int first, last;
  // int cpos, epos;
  // int size, nc;
  // FILE *fp;
  
  // fp = fopen(fn, "rb");
  // if(fp == NULL)
  //   {
  //     perror(fn);
  //   }
  // else
  //   {
  //     printf("\nChecking file: %s\n", fn);
  //   }

  /// File to be read
  FILE *fin=fopen(path,"r");
  if(fin==nullptr) CRASH("Unable to open file %s",path);
  
  file_head lh;
  fread(&lh.pars,sizeof(par_t),1,fin);
  fread(lh.kappa,sizeof(double),lh.pars.nk,fin);
  fread(lh.mu,sizeof(double),lh.pars.nk,fin);
  fread(lh.mom,sizeof(mom_t),lh.pars.nm,fin);
  fread(&lh.nobs,sizeof(int),1,fin);
  
  printf(" - twisted: %s\n", (lh.pars.twist ? "yes" : "no"));
  printf(" - flavors: %d\n", lh.pars.nf);
  printf(" - hits: %d\n", lh.pars.nhits);
  printf(" - volume: %dx%dx%dx%d\n", lh.pars.dim_t, lh.pars.dim_x, lh.pars.dim_y, lh.pars.dim_z);
  printf(" - masses: %d\n", lh.pars.nk);
  printf(" - momenta: %d\n", lh.pars.nm);
  printf(" - beta: %lf\n", lh.pars.beta);
  printf(" - kappa: %lf\n", lh.pars.ksea);
  printf(" - mu: %lf\n", lh.pars.musea);
  printf(" - csw: %lf\n", lh.pars.csw);
  
  for(int n = 0; n < lh.pars.nk; n++)
    printf(" - mass#%d: ( %f , %f )\n", n, lh.kappa[n], lh.mu[n]);
  
  for(int n = 0; n < lh.pars.nm; n++)
    printf(" - momenta#%d: ( %f , %f , %f, %f )\n", n, lh.mom[n].p0, lh.mom[n].p1, lh.mom[n].p2, lh.mom[n].p3);
  
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
  
  fclose(fin);
}

int main(int narg,char **arg)
{
  check_file("oAmuGPo-gs_conf.realph.dat");
  
  return 0;
}
