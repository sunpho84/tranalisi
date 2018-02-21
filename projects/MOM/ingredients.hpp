#ifndef _INGREDIENTS_HPP
#define _INGREDIENTS_HPP

#include <tranalisi.hpp>

#include <contractions.hpp>
#include <read.hpp>
#include <timings.hpp>
#include <types.hpp>

#ifndef EXTERN_INGREDIENTS
 #define EXTERN_INGREDIENTS extern
 #define INIT_TO(...)
#else
 #define INIT_TO(...) = __VA_ARGS__
#endif

EXTERN_INGREDIENTS index_t iGl_ipGl_iclust_ind;    //meslep
EXTERN_INGREDIENTS index_t iop_ipGl_iclust_ind;    //meslep
EXTERN_INGREDIENTS index_t iop_iproj_iclust_ind;   //meslep

//! kind of scheme supported
namespace reno_scheme
{
  enum type_t{RI_MOM,SMOM};
  const map<string,type_t> decr={{"RI_MOM",RI_MOM},{"SMOM",SMOM}};
  const size_t n=2;
  PROVIDE_DECRYPTER;
}

//! chiral extrapolation
namespace chir_extr
{
  enum type_t{MQUARK,MMES};
  const map<string,type_t> decr={{"MQuark",MQUARK},{"MMes",MMES}};
  const size_t n=2;
  PROVIDE_DECRYPTER;
}
EXTERN_INGREDIENTS chir_extr::type_t chir_extr_method;

//! reno scheme used
EXTERN_INGREDIENTS reno_scheme::type_t scheme;

//! passed to force recomputing Zbil
static const bool RECOMPUTE_ZBIL=true;

//! passed to disable recomputing Zbil
static const bool DONT_RECOMPUTE_ZBIL=false;

struct ingredients_t
{
  size_t _nm;
  size_t _nr;
  vector<double> _am;
  
  vector<array<size_t,1>> linmoms; //!< list of momenta used for Z, relative to glb list
  vector<array<size_t,3>> bilmoms; //!< list of momenta used for bilinear, first relative to glb list, then to linmoms
  vector<array<size_t,3>> &meslepmoms=bilmoms; //!< list of momenta used for meslep
  
  //Zq, with and without QED
  djvec_t Zq;
  djvec_t Zq_sig1;
  djvec_t Zq_sig1_QED;
  
  //! task to do someting
  struct task_t
  {
    //! input vector
    const djvec_t *in;
    //! output vector
    djvec_t *out;
    //! name of the task
    const string tag;
    //! contructor
    task_t(const djvec_t *in,djvec_t *out,const string tag) : in(in),out(out),tag(tag) {}
  };
  
  //! return a list of tasks for Zq
  vector<task_t> get_Zq_tasks(ingredients_t &out) const
  {
    vector<task_t> Zq_tasks={{&Zq,&out.Zq,"Zq"},{&Zq_sig1,&out.Zq_sig1,"Zq_sig1"}};
    if(use_QED) Zq_tasks.push_back({&Zq_sig1_QED,&out.Zq_sig1_QED,"Zq_sig1_QED"});
    return Zq_tasks;
  }
  
  //projected bilinears with and without QED
  djvec_t pr_bil;
  djvec_t pr_bil_QED;
  
  //! return a list of tasks for bilinears projected vertex
  vector<task_t> get_pr_bil_tasks(ingredients_t &out) const
  {
    vector<task_t> pr_bil_tasks={{&pr_bil,&out.pr_bil,"pr_bil"}};
    if(use_QED) pr_bil_tasks.push_back({&pr_bil_QED,&out.pr_bil_QED,"pr_bil_QED"});
    return pr_bil_tasks;
  }
  
  //bilinear Z
  djvec_t Zbil;
  djvec_t Zbil_QED;
  
  //! return a list of tasks for bilinears Z
  vector<task_t> get_Zbil_tasks(ingredients_t &out) const
  {
    vector<task_t> Zbil_tasks={{&Zbil,&out.Zbil,"Zbil"}};
    if(use_QED) Zbil_tasks.push_back({&Zbil_QED,&out.Zbil_QED,"Zbil_QED"});
    
    return Zbil_tasks;
  }
  
  index_t im_r_im_r_igam_ind;
  index_t r_r_iZbil_ind;
  index_t im_r_ijack_ind;
  index_t im_r_ijackp1_ind;
  index_t r_ilinmom_ind;
  index_t im_r_ilinmom_ind;
  index_t iZbil_ibilmom_ind;
  index_t im_r_ind;
  index_t im_r_im_r_iZbil_ind;
  index_t im_r_im_r_iZbil_ibilmom_ind;
  
  index_t im_r_im_r_iop_ilistpGl_ind;
  index_t im_r_im_r_iop_iproj_imeslepmom_ind;
  
  //! compute all mom-scheme propagator
  void mom_compute_qprop();
  
  //! compute all mom-scheme vertices
  void mom_compute_bil();
  
  //! projected meslep
  djvec_t pr_meslep;
  
  //! renormalization meslep
  djvec_t Zmeslep;
  index_t im_r_im_r_ilistGl_ipGl_ind;
  index_t im_r_im_r_iop_iproj_ind;
  
  //! compute all mom-scheme mesoleptonic vertices
  void mom_compute_meslep();
  
  //! compute according to ri-mom scheme
  void ri_mom();
  
  //! compute according to smom scheme
  void smom();
  
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
  
  //! set momenta for ri-mom
  void set_ri_mom_moms();
  
  //! set momenta for smom
  void set_smom_moms();
  
  //! set basic pars
  void set_pars_for_scratch();
  
  //! set all indices
  void set_indices();
  
  //! computes all bilinear Z
  void compute_Zbil();
  
  //! average r
  ingredients_t average_r(const bool recompute_Zbil) const;
  
  //allocate all vectors
  void allocate();
  
  //! extrapolate to chiral limit
  ingredients_t chir_extrap() const;
  
  //! subtract O(a^2)
  ingredients_t subtract_Oa2(const bool recompute_Zbil) const;
  
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
  ingredients_t average_equiv_momenta(const bool recompute_Zbil) const;
  
  //! plots all Z
  void plot_Z(const string &suffix="") const;
  
  //! return the list of momenta equivalent, in the combo passed
  template <size_t N>
  vector<vector<size_t>> get_equiv_list(const vector<array<size_t,N>> &combo,const string &listpath) const
  {
    //prepare the index of the output
    map<array<imom_t,N>,vector<size_t>> equiv_combo_map;
    for(size_t icombo=0;icombo<combo.size();icombo++)
      {
	//get representative
	array<imom_t,N> cr;
	for(size_t ip=0;ip<N;ip++)
	  {
	    size_t imom=combo[icombo][ip];
	    
	    //dereference components larger than 0
	    if(ip>0) imom=linmoms[imom][0];
	    
	    //decide time component
	    cr[ip][0]=all_moms[imom][0];
	    if(cr[ip][0]<0)
	      {
		using namespace temporal_bc;
		
		switch(bc)
		  {
		  case PERIODIC:     cr[ip][0]=-cr[ip][0];break;
		  case ANTIPERIODIC: cr[ip][0]=-cr[ip][0]-1;break;
		  }
	      }
	    
	    //decide space components
	    for(size_t mu=1;mu<NDIM;mu++) cr[ip][mu]=abs(all_moms[imom][mu]);
	    
	    //sort the components, treating differently the time component if it has different L[0] or phase
	    size_t ifirst=0;
	    if(temporal_bc::bc!=temporal_bc::PERIODIC or L[0]!=L[1]) ifirst=1;
	    sort(&cr[ip][ifirst],cr[ip].end());
	  }
	
	//sort starting from the defining momenta
	sort(&cr[1],cr.end());
	
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
	mom_out<<imom<<" , p2hat: "<<all_moms[imom].p(L).tilde().norm2()<<endl;
	mom_out<<" Equivalent to: "<<combo_class.second.size()<<" mom combos: "<<endl;
	for(auto &eq_combo : combo_class.second)
	  {
	    for(size_t ip=0;ip<N;ip++)
	      {
		size_t imom=combo[eq_combo][ip];
		
		//dereference
		if(ip>0) imom=linmoms[imom][0];
		
		mom_out<<"  "<<imom<<all_moms[imom];
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
