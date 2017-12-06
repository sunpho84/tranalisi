#ifndef _INGREDIENTS_HPP
#define _INGREDIENTS_HPP

#include <tranalisi.hpp>

#include <contractions.hpp>
#include <read.hpp>
#include <timings.hpp>

#ifndef EXTERN_INGREDIENTS
 #define EXTERN_INGREDIENTS extern
 #define INIT_TO(...)
#else
 #define INIT_TO(...) = __VA_ARGS__
#endif

EXTERN_INGREDIENTS index_t im_r_im_r_igam_ind;
EXTERN_INGREDIENTS index_t r_r_iZbil_ind;
EXTERN_INGREDIENTS index_t im_r_ijack_ind;
EXTERN_INGREDIENTS index_t im_r_ijackp1_ind;

//! kind of scheme supported
enum scheme_t{RI_MOM};

//! reno scheme used
EXTERN_INGREDIENTS scheme_t scheme;

//! passed to force recomputing Zbil
static const bool FORCE_RECOMPUTE_ZBIL=true;

struct ingredients_t
{
  size_t _nm;
  size_t _nr;
  vector<double> _am;
  
  vector<array<size_t,1>> linmoms; //!< list of momenta used for Z, relative to glb list
  vector<array<size_t,3>> bilmoms; //!< list of momenta used for bilinear, forst relative to glb list, then to linmoms
  
  //Zq, with and without EM
  djvec_t Zq;
  djvec_t Zq_sig1;
  djvec_t Zq_sig1_EM;
  
  //projected bilinears with and without EM
  djvec_t pr_bil;
  djvec_t pr_bil_QED;
  
  //bilinear Z
  bool Zbil_computed{false};
  djvec_t Zbil;
  djvec_t Zbil_QED;
  
  index_t r_ilinmom_ind;
  index_t im_r_ilinmom_ind;
  index_t iZbil_ibilmom_ind;
  index_t im_r_ind;
  index_t im_r_im_r_iZbil_ind;
  index_t im_r_im_r_iZbil_ibilmom_ind;
  
  //! fill all the ingredients according to ri-mom
  void ri_mom();
  
  //! read from binary
  void bin_read(raw_file_t &file);
  
  //! write to binary
  void bin_write(raw_file_t &file) const;
  
  //! read from a path
  inline void bin_read(const string &path)
  {
    raw_file_t fin(path,"r");
    this->bin_read(fin);
  }
  
  //! write to a path
  inline void bin_write(const string &path) const
  {
    raw_file_t fout(path,"w");
    this->bin_write(fout);
  }
  
  //! create using files
  void create_from_scratch(const string ingredients_path="ingredients.dat");
  
  //! set momenta for ri-moms
  void set_ri_mom_moms();
  
  //! set basic pars
  void set_pars_for_scratch();
  
  //! set all indices
  void set_indices();
  
  //! computes all bilinear Z
  void compute_Zbil();
  
  //allocate all vectors
  void allocate();
  
  //! extrapolate to chiral limit
  ingredients_t chir_extrap() const;
  
  //! subtract O(a^2)
  ingredients_t subtract_Oa2(const bool recompute_Zbil=false) const;
  
  //! evolve to a common scale
  ingredients_t evolve() const;
  
  //! helper to fill the momenta index of output
  template <size_t N>
  void fill_output_equivalent_momenta(vector<array<size_t,N>> &out_equiv,const vector<vector<size_t>> &equiv_linmom_combos,
				      const vector<vector<size_t>> &equiv_combos,const vector<array<size_t,N>> &mom_combos) const
  {
    out_equiv.clear();
    out_equiv.reserve(equiv_combos.size());
    
    for(const auto &p : equiv_combos)
      {
	const size_t &ref=p[0];
	const size_t &imom=mom_combos[ref][0];
	
	array<size_t,N> temp;
	temp[0]=imom;
	
	for(size_t i=1;i<N;i++)
	  {
	    //search in which output lin to search for lin
	    const size_t ilin_search=mom_combos[ref][i];
	    size_t &iout_lin=temp[i];
	    size_t nfound=0;
	    for(size_t ilin_combo=0;ilin_combo<equiv_linmom_combos.size();ilin_combo++)
	      {
		const vector<size_t> &equiv_list=equiv_linmom_combos[ilin_combo];
		auto ei=find(equiv_list.begin(),equiv_list.end(),ilin_search);
		if(ei!=equiv_list.end())
		  {
		    nfound++;
		    iout_lin=ilin_combo;
		  }
	      }
	    if(nfound!=1) CRASH("found %zu times %zu",nfound,ilin_search);
	  }
	
	out_equiv.push_back(temp);
      }
  }
  
  void fill_output_equivalent_momenta() const;
  
  //! average equivalent momenta
  ingredients_t average_equiv_momenta(const bool recompute_Zbil=false) const;
  
  //! plots all Z
  void plot_Z(const string &suffix="") const;
  
  //! return the list of momenta equivalent, in the combo passed
  template <size_t N>
  vector<vector<size_t>> get_equiv_list(const vector<array<size_t,N>> &combo,const string &listpath) const
  {
    //periodicity
    enum pe_t{UNK,PER,APE};
    const string pe_tag[3]={"UNK","PE","APE"};
    pe_t pe=UNK;
    double p0=fabs(ph_mom[0]);
    if(fabs(p0)<1e-10) pe=PER;
    if(fabs(p0-0.5)<1e-10) pe=APE;
    cout<<"Periodicity: "<<pe_tag[pe]<<endl;
    
    //prepare the index of the output
    map<array<imom_t,N>,vector<size_t>> equiv_combo_map;
    for(size_t icombo=0;icombo<combo.size();icombo++)
      {
	//get representative
	array<imom_t,N> cr;
	for(size_t ip=0;ip<N;ip++)
	  {
	    size_t imom=combo[icombo][ip];
	    if(ip>0) imom=linmoms[imom][0];
	    
	    //decide time component
	    cr[ip][0]=glb_moms[imom][0];
	    if(cr[ip][0]<0)
	      switch(pe)
		{
		case UNK: CRASH("phase on momentum 0 cannot be %lg",ph_mom[0]);break;
		case PER: cr[ip][0]=-cr[ip][0];break;
		case APE: cr[ip][0]=-cr[ip][0]-1;break;
		}
	    
	    //decide space components
	    for(size_t mu=1;mu<NDIM;mu++) cr[ip][mu]=abs(glb_moms[imom][mu]);
	    sort(&cr[ip][1],cr[ip].end());
	  }
	sort(cr.begin(),cr.end());
	
	//store the index to equvalents
	equiv_combo_map[cr].push_back(icombo);
      }
    
    //store the output vector to which the input contributes
    const size_t nequiv=equiv_combo_map.size();
    vector<vector<size_t>> equiv_of_combo;
    equiv_of_combo.reserve(nequiv);
    
    //fill the output index and write it to file
    ofstream mom_out(listpath);
    mom_out<<"Found "<<nequiv<<" independent momenta combo"<<endl;
    for(auto &combo_class : equiv_combo_map)
      {
	//trasform map into vector
	equiv_of_combo.push_back(combo_class.second);
	const size_t combo_repr=combo_class.second[0]; //take the first real combo as the key could not exist
	const size_t imom=combo[combo_repr][0];
	
	//print details
	mom_out<<imom<<" , p2hat: "<<glb_moms[imom].p(L).tilde().norm2()<<endl;
	mom_out<<" Equivalent to: "<<combo_class.second.size()<<" mom combos: "<<endl;
	for(auto &eq_combo : combo_class.second)
	  {
	    for(size_t ip=0;ip<N;ip++)
	      {
		const size_t imom=combo[eq_combo][ip];
		mom_out<<"  "<<imom<<"=";
		//components
		mom_out<<"{";
		for(size_t mu=0;mu<NDIM;mu++)
		  {
		    if(mu) mom_out<<",";
		    mom_out<<glb_moms[imom][mu];
		  }
		mom_out<<"}";
	      }
	    mom_out<<endl;
	  }
      }
    
  return equiv_of_combo;
  }
};

#undef EXTERN_INGREDIENTS
#undef INIT_TO

#endif
