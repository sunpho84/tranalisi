#ifndef _BASE_HPP
#define _BASE_HPP

#include <array>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <tranalisi.hpp>
#include <compatibility.hpp>

const char VA_tag[2][2]={"V","A"};

const char ST_tag[2]={'S','T'};

//! Holds all info for a given physical quark
using quark_t=tuple<double,vector<size_t>>;

//! List of physical quarks
//!
//! Each entry in the map is connected to a physical quark (Up, Down, etc)
//! Each value is a tuple: the first component is the charge, the second is a vector, listing all masses in the data
map<string,quark_t> quarkList;

//! Holds all info for a given meson: the first entry is the name, the second is the bw quatk, the third is the fw one
using meson_t=tuple<string,string,string>;

//! List of all known mesons
//!
//! Each entry is connected to a physical meson (Pi+, D, etc)
//! Each value is a pair of strings, pointing to the physical quark
vector<meson_t> mesonList;

/////////////////////////////////////////////////////////////////

//! total number of possible ensemble
const size_t nens_total=15;

//! number of beta
const size_t nbeta=3;

//! class to hold results from mass determination
class lat_par_t
{
public:
  dboot_t ml,ms,mc,r0,f0,B0;
  dbvec_t ainv,Z;
 
  lat_par_t() : ainv(nbeta),Z(nbeta) {}
};

//! number of input analysis
const size_t ninput_an=8;

vector<lat_par_t> lat_par(ninput_an);

//! number of boostrap
const size_t nboots=100;

boot_init_t jack_index[ninput_an][nens_total];

dboot_t read_boot(const raw_file_t &file)
{
  dboot_t out;
  for(size_t ib=0;ib<nboots;ib++) file.read(out[ib]);
  out.fill_ave_with_components_ave();
  return out;
}

void a()
{
  raw_file_t file("input_global","r");
  
  double dum;
  file.expect({"ml","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].ml=read_boot(file);
  file.expect({"ms","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].ms=read_boot(file);
  file.expect({"mc(2","GeV)","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].mc=read_boot(file);
  file.expect({"a^-1","(GeV)","(1.90","1.95","2.10)"});
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
    {
      for(size_t iboot=0;iboot<nboots;iboot++)
	for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	  file.read(lat_par[input_an_id].ainv[ibeta][iboot]);
      for(size_t ibeta=0;ibeta<nbeta;ibeta++) lat_par[input_an_id].ainv[ibeta].fill_ave_with_components_ave();
    }
  file.expect({"r0","(GeV^-1)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].r0=read_boot(file);
  file.expect({"Zp","(1.90","1.95","2.10)"});
  for(int ihalf=0;ihalf<2;ihalf++)
    {
      //skip average of the half
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	file.read(dum);
      
      size_t batch=ninput_an/2;
      for(size_t input_an_id=batch*ihalf;input_an_id<batch*(ihalf+1);input_an_id++)
	{
	  for(size_t iboot=0;iboot<nboots;iboot++)
	    for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	      file.read(lat_par[input_an_id].Z[ibeta][iboot]);
	  for(size_t ibeta=0;ibeta<nbeta;ibeta++) lat_par[input_an_id].Z[ibeta].fill_ave_with_components_ave();
	}
    }
  file.expect({"Jackknife","numbers","(","0.0030(32),0.0040(32),","0.0050(32),","0.0040(24)","...)"});
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
    for(size_t iboot=0;iboot<nboots;iboot++)
      for(size_t iens=0;iens<nens_total;iens++)
	{
	  size_t ijack_plus_one;
	  file.read(ijack_plus_one);
	  jack_index[input_an_id][iens][iboot]=ijack_plus_one-1;
	}
  file.expect({"f0","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].f0=read_boot(file);
  file.expect({"2*B0","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].B0=read_boot(file)/2.0;
}

#endif
