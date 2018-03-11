#ifndef _PERENS_HPP
#define _PERENS_HPP

#include <tranalisi.hpp>

#include <MOM2/geometry.hpp>
#include <MOM2/pars.hpp>
#include <MOM2/Zmeslep.hpp>
#include <MOM2/Zbil.hpp>

struct perens_t
{
  //! path where to find all data
  string dir_path;
  
  //! beta
  double beta;
  //! plaquette
  double plaq;
  
  //! minimal time for fit
  size_t tmin;
  //! maximal time for fit
  size_t tmax;
  //! inverse lattice spacing
  double ainv;
  //! conf range
  range_t conf_range;
  //! hadron path
  string prop_hadr_path;
  //! lepton path
  string prop_lep_path;
  
  //! list of existing confs
  vector<size_t> conf_list;
  //! cluster size
  size_t clust_size;
  //! number of hits
  size_t nhits;
  //! number of hits to use
  size_t nhits_to_use;
  
  //! lattice sizes
  coords_t L;
  //! lattice volume
  size_t V;
  
  //! number of masses
  size_t nm;
  //! number of r
  size_t nr;
  //! list of masses
  vector<double> am;
  
  //! minimal mass
  const double& am_min() const
  {
    return *min_element(am.begin(),am.end());
  }
  
  //! maximal mass
  const double& am_max() const
  {
    return *max_element(am.begin(),am.end());
  }
  
  /////////////////////////////////////////////////////////////////
  
  double g2() const
  {
    return 6.0/beta;
  }
  
  double g2tilde() const
  {
    return g2()/plaq;
  }
  
  /////////////////////////////////////////////////////////////////
  
  //! list of momenta registered in the system
  vector<imom_t> all_moms;
  //! number of momenta for which the propagator has been computed (first in list)
  size_t ncomp_moms;
  //! store if a momentum passed the filter
  vector<bool> filt_moms;
  
  //! list of momenta used for Z, relative to glb list
  vector<array<size_t,1>> linmoms;
  //! list of momenta used for bilinear, first relative to glb list, then to linmoms
  vector<array<size_t,3>> bilmoms;
  
  //! list of momenta used for meslep
  const vector<array<size_t,3>> &meslepmoms() const
  {
    return bilmoms;
  }
  
  /////////////////////////////////////////////////////////////////
  
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
  
  //Zq, with and without QED
  djvec_t Zq;
  djvec_t Zq_sig1;
  djvec_t Zq_QED;
  djvec_t Zq_sig1_QED;
  
  //! return a list of tasks for Zq
  vector<task_t> get_Zq_tasks(perens_t &out) const
  {
    vector<task_t> Zq_tasks={{&Zq,&out.Zq,"Zq"},{&Zq_sig1,&out.Zq_sig1,"Zq_sig1"}};
    if(pars::use_QED)
      {
	Zq_tasks.push_back({&Zq_QED,&out.Zq_QED,"Zq_QED"});
	Zq_tasks.push_back({&Zq_sig1_QED,&out.Zq_sig1_QED,"Zq_sig1_QED"});
      }
    return Zq_tasks;
  }
  
  //projected bilinears with and without QED
  djvec_t pr_bil;
  djvec_t pr_bil_QED;
  
  //! return a list of tasks for bilinears projected vertex
  vector<task_t> get_pr_bil_tasks(perens_t &out) const
  {
    vector<task_t> pr_bil_tasks={{&pr_bil,&out.pr_bil,"pr_bil"}};
    if(pars::use_QED) pr_bil_tasks.push_back({&pr_bil_QED,&out.pr_bil_QED,"pr_bil_QED"});
    return pr_bil_tasks;
  }
  
  //bilinear Z
  djvec_t Zbil;
  djvec_t Zbil_QED;
  
  //! return a list of tasks for bilinears Z
  vector<task_t> get_Zbil_tasks(perens_t &out) const
  {
    vector<task_t> Zbil_tasks={{&Zbil,&out.Zbil,"Zbil"}};
    if(pars::use_QED) Zbil_tasks.push_back({&Zbil_QED,&out.Zbil_QED,"Zbil_QED"});
    
    return Zbil_tasks;
  }
  
  //! projected meslep
  djvec_t pr_meslep;
  djvec_t pr_meslep_QED;
  
  //! return a list of tasks for meslep projected vertex
  vector<task_t> get_pr_meslep_tasks(perens_t &out) const
  {
    vector<task_t> pr_meslep_tasks;
    if(pars::compute_meslep)
      {
	pr_meslep_tasks.push_back({&pr_meslep,&out.pr_meslep,"pr_meslep"});
	if(pars::use_QED) pr_meslep_tasks.push_back({&pr_meslep_QED,&out.pr_meslep_QED,"pr_meslep_QED"});
      }
    return pr_meslep_tasks;
  }
  
  //! renormalization meslep
  djvec_t Zmeslep;
  djvec_t Zmeslep_QED;
  
  //! return a list of tasks for meslep projected vertex
  vector<task_t> get_Zmeslep_tasks(perens_t &out) const
  {
    vector<task_t> Zmeslep_tasks;
    if(pars::compute_meslep)
      {
	Zmeslep_tasks.push_back({&Zmeslep,&out.Zmeslep,"Zmeslep"});
	if(pars::use_QED) Zmeslep_tasks.push_back({&Zmeslep_QED,&out.Zmeslep_QED,"Zmeslep_QED"});
      }
    return Zmeslep_tasks;
  }
  
  //! allocate all data
  perens_t& allocate();
  
  //! return a list of all tasks
  vector<task_t> get_all_tasks(perens_t &out) const
  {
    return concat(get_Zq_tasks(out),
		  get_pr_bil_tasks(out),
		  get_Zbil_tasks(out),
		  get_pr_meslep_tasks(out),
		  get_Zmeslep_tasks(out));
  }
  
  /////////////////////////////////////////////////////////////////
  
  djvec_t deltam_cr;
  
  djvec_t meson_mass;
  
  //! read or compute and write deltam_cr
  perens_t& get_deltam_cr();
  
  //! read or compute and write meson mass
  perens_t& get_meson_mass();
  
  /////////////////////////////////////////////////////////////////
  
  index_t im_im_ind;
  
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
  index_t im_r_im_r_ilistGl_ipGl_ind;
  index_t im_r_im_r_iop_iproj_ind;
  
  index_t iGl_ipGl_iclust_ind;
  index_t iop_ipGl_iclust_ind;
  index_t iop_iproj_iclust_ind;
  
  index_t im_r_iconf_ihit_iqkind_ind;
  index_t iconf_ihit_ilkind_ind;
  
  index_t i_in_clust_ihit_ind;
  index_t conf_ind;
  
  index_t ilistGl_ilistpGl_iclust_ind;
  
  perens_t& set_indices();
  
  /////////////////////////////////////////////////////////////////
  
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
  
  /////////////////////////////////////////////////////////////////
  
  //! load parameters and create structures
  perens_t& read_pars(const string &path);
  
  //! set all momenta
  void set_comp_list_of_moms(const string &mom_list_path,double filter_thresh);
  
  //! gets the mirrored site
  size_t get_mir_mom(size_t imom,size_t imir);
  
  //! set momenta for ri-mom
  void set_ri_mom_moms();
  
  //! set momenta for s-mom
  void set_smom_moms();
  
  /////////////////////////////////////////////////////////////////
  
  //! set parameters for the scratch
  perens_t& set_pars_for_scratch();
  
  //! open to read all files
  vector<raw_file_t> setup_read_all_qprops_mom(const vector<size_t> &conf_list) const;
  
  //! opem to read all lepton files
  vector<raw_file_t> setup_read_all_lprops_mom(const vector<size_t> &conf_list) const;
  
  //! prepares the list of configurations to use
  void prepare_list_of_confs();
  
  //! try to read, otherwise compute
  perens_t& read_or_compute();
  
  //! computes the basic Z
  perens_t& compute_basic(const string& ingredients_path);
  
  //! compute according to mom scheme
  void mom()
  {
    mom_compute_qprop();
    mom_compute_bil();
    if(pars::compute_meslep) mom_compute_meslep();
  }
  
  double compute_Zq(const qprop_t &prop_inv,const size_t glb_mom);
  
  djack_t compute_Zq(const jqprop_t &jprop_inv,const size_t glb_mom);
  
  double compute_Zq_sig1(const qprop_t &prop_inv,const size_t glb_mom);
  
  djack_t compute_Zq_sig1(const jqprop_t &jprop_inv,const size_t glb_mom);
  
  vector<m_r_mom_conf_qprops_t> read_all_qprops_mom(vector<raw_file_t> &files,const size_t i_in_clust_ihit,const size_t imom);
  
  void get_inverse_propagators(vector<jqprop_t> &jprop_inv,vector<jqprop_t> &jprop_QED_inv,
			       const vector<jm_r_mom_qprops_t> &jprops) const;
  
  //! compute all props
  void mom_compute_qprop();
  
  vector<mom_conf_lprops_t> read_all_lprops_mom(vector<raw_file_t> &files,const size_t i_in_clust_ihit,const size_t imom);
  
  void build_all_mr_jackkniffed_qprops(vector<jm_r_mom_qprops_t> &jprops,const vector<m_r_mom_conf_qprops_t> &props) const;
  
  void clusterize_all_mr_jackkniffed_qprops(vector<jm_r_mom_qprops_t> &jprops) const;
  
  /////////////////////////////////////////////////////////////////
  
  djack_t compute_meson_mass(const size_t im1,const size_t im2);
  
  djvec_t get_contraction(const string &combo,const string &ID,const size_t reim,const int tpar);
  
  djack_t compute_deltam_cr(const size_t im);
  
  /////////////////////////////////////////////////////////////////
  
  //! computes all bilinear Z
  void compute_Zbil();
  
  //! compute all mom-scheme vertices
  void mom_compute_bil();
  
  void build_all_mr_gbil_jackkniffed_verts(jbil_vert_t &jbil,const vector<m_r_mom_conf_qprops_t> &props_in,const vector<m_r_mom_conf_qprops_t> &props_ou) const;
  
  djvec_t compute_proj_bil(const vjqprop_t &jprop_inv_in,const vector<jqprop_t> &jverts,const vjqprop_t &jprop_inv_ou) const;
  
  /////////////////////////////////////////////////////////////////
  
  //! computes all meslep Z
  void compute_Zmeslep();
  
  //! compute all mom-scheme mesoleptonic vertices
  void mom_compute_meslep();
  
  meslep::mesloop_t build_mesloop(const vector<mom_conf_lprops_t> &props_lep) const;
  
  void build_all_mr_gmeslep_jackkniffed_verts(jmeslep_vert_t &j,const vector<m_r_mom_conf_qprops_t> &props_in,const vector<m_r_mom_conf_qprops_t> &props_ou,
					      const vector<mom_conf_lprops_t> &props_lep) const;
  
  djvec_t compute_proj_meslep(const vjqprop_t &jprop_inv_in,const vector<jqprop_t> &jverts,const vjqprop_t &jprop_inv_ou) const;
  
  /////////////////////////////////////////////////////////////////
  
  void plot_Z(const string &suffix)
  {
    plot_Zq(suffix);
    plot_Zbil(suffix);
    plot_Zmeslep(suffix);
  }
  
  void plot_Zq(const string &suffix);
  void plot_Zbil(const string &suffix);
  void plot_Zmeslep(const string &suffix);
  
  /////////////////////////////////////////////////////////////////
  
  perens_t average_r() const;
  
  //! average Zq
  void average_r_Zq(perens_t &out) const;
  //! average bil
  void average_r_Zbil(perens_t &out) const;
  //! average Zmeslep
  void average_r_Zmeslep(perens_t &out) const;
  
  /////////////////////////////////////////////////////////////////
  
  perens_t val_chir_extrap() const;
  
  //extrapolate to chiral limit Zq
  void val_chir_extrap_Zq(perens_t &out) const;
  
  //extrapolate to chiral limit Zbil
  void val_chir_extrap_Zbil(perens_t &out) const;
  
  //extrapolate to chiral limit Zmeslep
  void val_chir_extrap_Zmeslep(perens_t &out) const;
  
  /////////////////////////////////////////////////////////////////
  
  perens_t average_equiv_momenta() const;
  
  //! average Zq
  void average_equiv_momenta_Zq(perens_t &out,const vector<vector<size_t>> &equiv_linmom_combos) const;
  //! average bil
  void average_equiv_momenta_Zbil(perens_t &out,const vector<vector<size_t>> &equiv_bilmom_combos) const;
  //! average Zmeslep
  void average_equiv_momenta_Zmeslep(perens_t &out,const vector<vector<size_t>> &equiv_meslepmom_combos) const;
  
  template <size_t N>
  void fill_output_equivalent_momenta(vector<array<size_t,N>> &out_equiv,const vector<vector<size_t>> &equiv_linmom_combos,
				      const vector<vector<size_t>> &equiv_combos,const vector<array<size_t,N>> &mom_combos) const;
  
  template <size_t N>
  vector<vector<size_t>> get_equiv_list(const vector<array<size_t,N>> &combo,const string &listpath) const;
};

//! return the list of momenta equivalent, in the combo passed
template <size_t N>
vector<vector<size_t>> perens_t::get_equiv_list(const vector<array<size_t,N>> &combo,const string &listpath) const
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

//! helper to fill the momenta index of output
template <size_t N>
void perens_t::fill_output_equivalent_momenta(vector<array<size_t,N>> &out_equiv,const vector<vector<size_t>> &equiv_linmom_combos,
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

#endif
