#ifndef _PERENS_HPP
#define _PERENS_HPP

#include <tranalisi.hpp>

#include <MOM2/geometry.hpp>
#include <MOM2/pars.hpp>
#include <MOM2/sigma.hpp>
#include <MOM2/pr_bil.hpp>
#include <MOM2/pr_meslep.hpp>
#include <MOM2/Zbil.hpp>
#include <MOM2/Zmeslep.hpp>

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
  
  //! index of sea mass
  int im_sea;
  
  /////////////////////////////////////////////////////////////////
  
  enum QCD_QED_task_t{QCD_task,QED_task};
  
  /////////////////////////////////////////////////////////////////
  
  double g2() const
  {
    return 6.0/beta;
  }
  
  double g2to_correct() const
  {
    if(pars::correct_Oa2_using_gtilde)
      return g2()/plaq;
    else
      return g2();
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
    const vector<const djvec_t*> in;
    //! output vector
    djvec_t* out;
    //! index of the output
    const index_t ind;
    //! name of the task
    const string tag;
    //! store if is QED
    const QCD_QED_task_t QCD_QED_task;
    //! contructor
    task_t(djvec_t *out,const vector<const djvec_t*>& in,const index_t& ind,const string tag,const QCD_QED_task_t QCD_QED_task) :
      in(in),out(out),ind(ind),tag(tag),QCD_QED_task(QCD_QED_task) {}
  };
  
  /////////////////////////////////////////////////////////////////
  
  //! Sigmas
  djvec_t sigma;
  
  //! returns a view on a sigma, with fixed im,r,ilinmom and proj
  auto sigma_ins_getter(const size_t im,const size_t r,const size_t ilinmom,sigma::proj proj)
  {
    using namespace sigma;
    
    return [this,im,r,ilinmom,proj](sigma::ins ins)->djack_t&
      {
	return sigma[im_r_ilinmom_isigmaproj_isigmains_ind({im,r,ilinmom,sigma::iproj_of_proj[proj],sigma::iins_of_ins[ins]})];
      };
  }
  
  //! compute all sigmas
  perens_t& compute_sigmas();
  
  //! plot all sigmas
  void plot_sigma(const string &suffix);
  
  //! return a list of tasks for sigma
  vector<task_t> get_sigma_tasks(const vector<const perens_t*> &ens={});
  
  //! put together the sigma needed for QED
  void assemble_sigma_QED_greenfunctions();
  
  /////////////////////////////////////////////////////////////////
  
  //! projected bilinears
  djvec_t pr_bil;
  
  //! returns a view on a pr_bil with a given im and r in and out
  auto pr_bil_ins_getter(const size_t im_ou,const size_t r_ou,const size_t im_in,const size_t r_in,const size_t ibil,const size_t ibilmom)
  {
    using namespace pr_bil;
    
    return [this,im_ou,r_ou,im_in,r_in,ibil,ibilmom](pr_bil::ins ins)->djack_t&
      {
	return pr_bil[im_r_im_r_bilins_ibil_ibilmom_ind({im_ou,r_ou,im_in,r_in,pr_bil::iins_of_ins[ins],ibil,ibilmom})];
      };
  }
  
  //! return a list of tasks for bilinears projected vertex
  vector<task_t> get_pr_bil_tasks(const vector<const perens_t*> &ens={});
  
  void build_all_mr_gbil_jackkniffed_verts(vector<jqprop_t>& jbil,const vector<qprop_t>& props_in,const vector<qprop_t>& props_ou) const;
  
  //! compute all bilinears
  void compute_proj_bil(const vector<jqprop_t>& jprop_inv_in,const vector<jqprop_t>& jverts,const vector<jqprop_t>& jprop_inv_ou,const size_t ibilmom);
  
  //! compute all mom-scheme vertices
  void mom_compute_bil();
  
  //! put together the pr_bil needed for QED
  void assemble_pr_bil_QED_greenfunctions();
  
  /////////////////////////////////////////////////////////////////
  
  //! projected meslep
  djvec_t pr_meslep;
  
  //! returns a view on a pr_meslep
  auto pr_meslep_ins_getter(const size_t im_ou,const size_t r_ou,const size_t im_in,const size_t r_in,const size_t iop,const size_t iproj,const size_t imeslepmom)
  {
    using namespace pr_meslep;
    
    return [this,im_ou,r_ou,im_in,r_in,iop,iproj,imeslepmom](pr_meslep::ins ins)->djack_t&
      {
	return pr_meslep[im_r_im_r_meslepins_iop_iproj_imeslepmom_ind({im_ou,r_ou,im_in,r_in,pr_meslep::iins_of_ins[ins],iop,iproj,imeslepmom})];
      };
  }
  
  //! return a list of tasks for meslep projected vertex
  vector<task_t> get_pr_meslep_tasks(const vector<const perens_t*> &ens={});
  
  //! compute all mom-scheme mesoleptonic vertices
  void mom_compute_meslep();
  
  //! put together the pr_meslep needed for QED
  void assemble_pr_meslep_QED_greenfunctions();
  
  //! make meslep Z of QED absolute
  void make_Zmeslep_QED_absolute();
  
  /////////////////////////////////////////////////////////////////
  
  //! allocate all data
  perens_t& allocate();
  
  //! get all ingredients
  vector<task_t> get_all_ingredients(const vector<const perens_t*> &ens={})
  {
    return concat(get_sigma_tasks(ens),
		  get_pr_bil_tasks(ens),
		  get_pr_meslep_tasks(ens),
		  get_deltam_tasks(ens),
		  get_meson_mass_tasks(ens));
  }
  
  //! return a list of all Ztasks
  vector<task_t> get_all_Ztasks(const vector<const perens_t*> &ens={})
  {
    return concat(get_Zq_tasks(ens),
		  get_Zbil_tasks(ens),
  		  get_Zmeslep_tasks(ens));
  }
  
  /////////////////////////////////////////////////////////////////
  
  //! correction to m_cr
  djvec_t deltam_cr;
  
  //! correction to m_tm
  djvec_t deltam_tm;
  
  //! print the value of deltam
  void print_deltam(ostream &out=cout) const;
  
  //! compose the ingredients_path
  string deltam_path() const
  {
    return dir_path+"/deltam.dat";
  }
  
  //! write the deltam
  void bin_write_deltam();
  
  //! read the deltam
  void bin_read_deltam();
  
  //! compose the meson mass path
  string meson_mass_path()
  {
    return dir_path+"/meson_mass.dat";
  }
  
  //! compose the meson mass QEC path
  string meson_mass_QED_path()
  {
    return dir_path+"/meson_mass_QED.dat";
  }
  
  //! compose the meson mass sea path
  string meson_mass_sea_path()
  {
    return dir_path+"/meson_mass_sea.dat";
  }
  
  //! mass of pseudoscalar meson
  djvec_t meson_mass;
  
  //! mass of pseudoscalar meson QED
  djvec_t meson_mass_QED;
  
  //! mass of sea pseudoscalar meson
  djack_t meson_mass_sea;
  
  //! write the meson_mass
  void bin_write_meson_mass();
  
  //! read the meson_mass
  void bin_read_meson_mass();
  
  //! write the meson_mass in QED
  void bin_write_meson_mass_QED();
  
  //! read the meson_mass in QED
  void bin_read_meson_mass_QED();
  
  //! write the meson_mass_sea
  void bin_write_meson_mass_sea();
  
  //! read the meson_mass_sea
  void bin_read_meson_mass_sea();
  
  //! compute or recompute deltam
  void recompute_deltam();
  
  //! read or compute and write deltam
  perens_t& get_deltam();
  
  //! read or compute and write meson mass
  perens_t& get_meson_mass();
  
  //! mass PCAC
  djvec_t mPCAC;
  
  //! mass PCAC for sea meson
  djack_t mPCAC_sea;
  
  //! read or compute and write mPCAC
  perens_t& get_mPCAC();
  
  //! returns a list of all tasks for deltam
  vector<task_t> get_deltam_tasks(const vector<const perens_t*> &ens={});
  
  //! returns a list of all tasks for meson masses
  vector<task_t> get_meson_mass_tasks(const vector<const perens_t*> &ens={});
  
  /////////////////////////////////////////////////////////////////
  
  index_t im_im_ind;
  
  index_t im_r_im_r_igam_ind;
  index_t im_r_im_r_bilins_igam_ind;
  index_t im_r_im_r_bilins_ibil_ibilmom_ind;
  
  index_t im_r_im_r_meslepins_iop_iproj_ind;
  index_t im_r_im_r_meslepins_iop_iproj_imeslepmom_ind;
  
  index_t r_r_ibil_ind;
  index_t im_r_ijack_ind;
  index_t im_r_ijackp1_ind;
  index_t r_ilinmom_ind;
  index_t im_r_ilinmom_ind;
  index_t ibil_ibilmom_ind;
  index_t im_r_ind;
  index_t im_r_im_r_ibil_ind;
  index_t im_r_im_r_ibil_ibilmom_ind;
  
  index_t im_r_ilinmom_isigmaproj_isigmains_ind;
  
  index_t im_r_im_r_iop_ilistpGl_ind;
  index_t im_r_im_r_iop_iproj_imeslepmom_ind;
  index_t im_r_im_r_ilistGl_ipGl_ind;
  index_t im_r_im_r_iop_iproj_ind;
  
  index_t iGl_ipGl_iclust_ind;
  index_t iop_ipGl_iclust_ind;
  index_t iop_iproj_iclust_ind;
  
  index_t im_r_iconf_ihit_iqins_ind;
  index_t im_r_iqins_ijack_ind;
  index_t im_r_ijqins_ijack_ind;
  index_t im_r_ijqins_ind;
  index_t iconf_ihit_ilins_ind;
  index_t ilins_ijack_ind;
  
  index_t i_in_clust_ihit_ind;
  index_t conf_ind;
  
  index_t ilistGl_ilistpGl_iclust_ind;
  index_t ilistGl_ilistpGl_lins_iclust_ind;
  index_t im_r_im_r_iop_ilistpGl_meslepins_ind;
  
  perens_t& set_indices();
  
  /////////////////////////////////////////////////////////////////
  
  //! read from binary ingredients
  void bin_read_ingredients(raw_file_t &file);
  
  //! write to binary ingredients
  void bin_write_ingredients(raw_file_t &file);
  
  //! read ingredients from a path
  inline void bin_read_ingredients()
  {
    raw_file_t fin(ingredients_path(),"r");
    this->bin_read_ingredients(fin);
  }
  
  //! write ingredients to a path
  inline void bin_write_ingredients()
  {
    raw_file_t fout(ingredients_path(),"w");
    this->bin_write_ingredients(fout);
  }
  
  /////////////////////////////////////////////////////////////////
  
  //! load parameters and create structures
  perens_t& read_pars(const string &path);
  
  //! save parameters
  void write_pars(const string &path) const;
  
  //! set all momenta
  void set_comp_list_of_moms(const string &mom_list_path,double filter_thresh);
  
  //! write all momenta
  void write_comp_list_of_moms(const string &mom_list_path) const;
  
  //! gets the mirrored site
  size_t get_mir_mom(size_t imom,size_t imir);
  
  //! set momenta for ri-mom
  void set_ri_mom_moms();
  
  //! set momenta for s-mom
  void set_smom_moms();
  
  //! print all momenta with discretizations
  void print_discr();
  
  /////////////////////////////////////////////////////////////////
  
  //! set parameters for the scratch
  perens_t& set_pars_for_scratch();
  
  //! open to read all files
  vector<raw_file_t> setup_read_all_qprops_mom(const vector<size_t> &conf_list) const;
  
  //! opem to read all lepton files
  vector<raw_file_t> setup_read_all_lprops_mom(const vector<size_t> &conf_list) const;
  
  //! prepares the list of configurations to use
  void prepare_list_of_confs();
  
  //! compose the ingredients_path
  string ingredients_path() const
  {
    return dir_path+"/ingredients.dat";
  }
  
  //! try to read, otherwise compute the ingredients (sigma, projected bil, etc)
  perens_t& read_or_compute_ingredients();
  
  //! write the current status
  perens_t write_checkpoint();
  
  //! computes the basic quantities
  perens_t& compute_ingredients();
  
  //! compute according to mom scheme
  void mom()
  {
    compute_sigmas();
    if(pars::compute_bilinears) mom_compute_bil();
    if(pars::compute_meslep) mom_compute_meslep();
  }
  
  vector<qprop_t> read_all_qprops_mom(vector<raw_file_t> &files,const size_t i_in_clust_ihit,const size_t imom);
  
  vector<lprop_t> read_all_lprops_mom(vector<raw_file_t> &files,const size_t i_in_clust_ihit,const size_t imom);
  
  void build_all_mr_jackkniffed_qprops(vector<jqprop_t> &jprops,const vector<qprop_t> &props) const;
  
  //! returns the inverse propagators
  vector<jqprop_t> get_inverse_propagators(const vector<jqprop_t>& qprop) const;
  
  /////////////////////////////////////////////////////////////////
  
  djack_t compute_meson_mass(const string& m1_tag,const string& m2_tag);
  
  djack_t compute_meson_mass_QED(const size_t im1,const size_t im2);
  
  djack_t compute_mPCAC(const size_t im);
  
  djvec_t get_contraction_by_name(const string &suffix,const string &bil_name,const dcompl_t &coeff,const int tpar);
  
  //! load a given correlation function, forming the name on the basis of insertion, mass index, etc
  djvec_t get_contraction(const int imbw,qprop::ins kbw,const int imfw,qprop::ins kfw,const string &ID,const size_t ext_reim,const int tpar,const size_t rfw,const int rdiff);
  
  //! compute deltam from correlators
  void compute_deltam_from_corr();
  
  //! compute deltam from propagators
  void compute_deltam_from_prop();
  
  //! compute deltam
  void compute_deltam();
  
  /////////////////////////////////////////////////////////////////
  
  //Zq, with and without QED
  djvec_t Zq;
  djvec_t Zq_QED_rel;
  djvec_t Zq_RI;
  djvec_t Zq_RI_QED_rel;
  
  //! return a list of tasks for Zq
  vector<task_t> get_Zq_tasks(const vector<const perens_t*>& ens={});
  
  //! plot all Zq
  void plot_Zq(const string &suffix);
  
  //! compute all Zq
  void compute_Zq(const bool also_QCD,const bool also_QED);
  
  /////////////////////////////////////////////////////////////////
  
  //bilinear Z
  
  djvec_t Zbil;
  djvec_t Zbil_QED_rel;
  
  //! return a list of tasks for bilinears Z
  vector<task_t> get_Zbil_tasks(const vector<const perens_t*> &ens={});
  
  //! plot all Zbil
  void plot_Zbil(const string &suffix);
  
  //! computes all bilinear Z
  void compute_Zbil(const bool also_QCD,const bool also_QED);
  
  /////////////////////////////////////////////////////////////////
  
  //meslep Z
  
  djvec_t Zmeslep;
  djvec_t Zmeslep_QED_rel;
  
  //! return a list of tasks for meslep projected vertex
  vector<task_t> get_Zmeslep_tasks(const vector<const perens_t*> &ens={});
  
  //! computes all meslep Z
  void compute_Zmeslep(const bool also_QCD,const bool also_QED);
  
  void plot_Zmeslep(const string &suffix);
  
  vector<dcompl_t> build_mesloop(const vector<lprop_t> &props_lep) const;
  
  void build_all_mr_gmeslep_jackkniffed_verts(vector<jqprop_t> &j,const vector<qprop_t> &props_in,const vector<qprop_t> &props_ou,
					      const vector<lprop_t> &props_lep) const;
  
  void compute_proj_meslep(const vector<jqprop_t> &jprop_inv_in,const vector<jqprop_t> &jverts,const vector<jqprop_t> &jprop_inv_ou,const size_t imeslepmom);
  
  /////////////////////////////////////////////////////////////////
  
  //! assemble the green function
  perens_t assemble_QED_greenfunctions() const;
  
  //! compute all Z
  void compute_Z(const bool also_QCD=true,const bool also_QED=true)
  {
    if(also_QED and not pars::can_compute_QED) CRASH("Cannot compute QED now");
    
    compute_Zq(also_QCD,also_QED);
    if(pars::compute_bilinears) compute_Zbil(also_QCD,also_QED);
    if(pars::compute_meslep) compute_Zmeslep(also_QCD,also_QED);
  }
  
  //! compute all Z but not QED
  void compute_Z_QCD()
  {
    compute_Z(true,false);
  }
  
  void make_Z_QED_absolute();
  
  void plot_Z(const string &suffix);
  
  void print_Z(ofstream& file);
  
  /////////////////////////////////////////////////////////////////
  
  //! triggers r-averaging of all quantities
  perens_t average_r(const array<double,2> w={1,1}) const;
  
  //! average sigma
  void average_r_sigma(perens_t &out,const array<double,2> w={1,1}) const;
  
  //! average pr_bil
  void average_r_pr_bil(perens_t &out,const array<double,2> w={1,1}) const;
  
  //! average pr_meslep
  void average_r_pr_meslep(perens_t &out,const array<double,2> w={1,1}) const;
  
  //! average deltam
  void average_r_deltam(perens_t &out,const array<double,2> w={1,1}) const;
  
  /////////////////////////////////////////////////////////////////
  
  //! triggers r-selecting of all quantities
  perens_t select_r(const size_t r) const;
  
  /////////////////////////////////////////////////////////////////
  
  //! triggers cut-off effects remotion  of all quantities
  perens_t subtract_Oa2() const;
  
  //! correct cut-off effects for sigma
  void subtract_Oa2_sigma();
  
  //! correct cut-off effects for pr_bil
  void subtract_Oa2_pr_bil();
  
  //! correct cut-off effects for pr_meslep
  void subtract_Oa2_pr_meslep();
  
  /////////////////////////////////////////////////////////////////
  
  //! triggers evolve  of all quantities
  perens_t evolve() const;
  
  //! evolve sigma
  void evolve_sigma(perens_t& out) const;
  
  //! evolve pr_bil
  void evolve_pr_bil(perens_t& out) const;
  
  //! evolve pr_meslep
  void evolve_pr_meslep(perens_t& out) const;
  
  /////////////////////////////////////////////////////////////////
  
  //! triggers match to W-reg all quantities
  perens_t match_to_W_reg() const;
  
  //! match Zmeslep
  void match_to_W_reg_Zmeslep(perens_t& out) const;
  
  /////////////////////////////////////////////////////////////////
  
  //! returns the range in which x is contained
  template <size_t N>
  pair<double,double> get_a2p2tilde_range_bracketting(const vector<array<size_t,N>> &list,const double a2p2_ref,const size_t n=5) const;
  
  //! interpolate or exptrapolate to a fixed p2
  perens_t interpolate_or_extrapolate_to_single_p2_preamble() const;
  
  //! interpolate all quantities to reference p2
  perens_t interpolate_to_p2ref() const;
  
  //! interpolate Zq to reference p2
  void interpolate_Zq_to_p2ref(perens_t &out) const;
  
  //! interpolate Zbil to reference p2
  void interpolate_Zbil_to_p2ref(perens_t &out) const;
  
  //! interpolate Zmeslep to reference p2
  void interpolate_Zmeslep_to_p2ref(perens_t &out) const;
  
  /////////////////////////////////////////////////////////////////
  
  perens_t val_chir_extrap() const;
  
  //! extrapolate to chiral limit sigma
  void val_chir_extrap_sigma(perens_t &out) const;
  
  //! extrapolate to chiral limit pr_bil
  void val_chir_extrap_pr_bil(perens_t &out) const;
  
  //! extrapolate to chiral limit pr_meslep
  void val_chir_extrap_pr_meslep(perens_t &out) const;
  
  //! extrapolate to chiral limit deltas
  void val_chir_extrap_deltam(perens_t &out) const;
  
  /////////////////////////////////////////////////////////////////
  
  //! extrapolate to null momentum
  perens_t extrapolate_to_0_p2() const;
  
  //! internal part
  template <size_t N>
  void extrapolate_to_0_p2_internal(const vector<array<size_t,N>> &moms,const task_t &t) const;
  
  /////////////////////////////////////////////////////////////////
  
  perens_t average_equiv_momenta() const;
  
  //! average sigma
  void average_equiv_momenta_sigma(perens_t &out,const vector<vector<size_t>> &equiv_linmom_combos) const;
  //! average bil
  void average_equiv_momenta_pr_bil(perens_t &out,const vector<vector<size_t>> &equiv_bilmom_combos) const;
  //! average meslep
  void average_equiv_momenta_pr_meslep(perens_t &out,const vector<vector<size_t>> &equiv_meslepmom_combos) const;
  
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
      mom_out<<imom<<" , p2tilde: "<<all_moms[imom].p(L).tilde().norm2()<<endl;
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

template <size_t N>
pair<double,double> perens_t::get_a2p2tilde_range_bracketting(const vector<array<size_t,N>> &list,const double a2p2_ref,const size_t n) const
{
  //get the list of a2p2 by proximity
  vector<pair<double,size_t>> prox_list;
  vector<double> x(list.size());
  for(size_t imom=0;imom<list.size();imom++)
    {
      x[imom]=all_moms[list[imom][0]].p(L).tilde().norm2();
      const double dist=fabs(x[imom]-a2p2_ref);
      prox_list.push_back({dist,imom});
    }
  sort(prox_list.begin(),prox_list.end());
  
  if(n>prox_list.size()) CRASH("Interpolating with a number of points %zu, but we dispose of %zu",n,prox_list.size());
  
  //get min max
  double a2p2min=1e300,a2p2max=1e-300;
  for(size_t i=0;i<n;i++)
    {
      const size_t imom=prox_list[i].second;
      a2p2min=min(a2p2min,x[imom]);
      a2p2max=max(a2p2max,x[imom]);
    }
  
  return {a2p2min,a2p2max};
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

template <size_t N>
void perens_t::extrapolate_to_0_p2_internal(const vector<array<size_t,N>> &moms,const task_t &t) const
{
  for(size_t iout=0;iout<t.ind.max();iout++)
    {
      vector<double> x(moms.size());
      djvec_t y(moms.size());
      for(size_t imom=0;imom<moms.size();imom++)
	{
	  size_t iin=iout*moms.size()+imom;
	  x[imom]=all_moms[moms[imom][0]].p(L).tilde().norm2();
	  y[imom]=(*t.in[0])[iin];
	}
      
      const double xmin=pars::extrapolate_to_0_p2_range.first;
      const double xmax=pars::extrapolate_to_0_p2_range.second;
      
      const djvec_t coeffs=poly_fit(x,y,1,xmin,xmax);
      (*t.out)[iout]=coeffs[0];
      
      //plot
      grace_file_t plot(dir_path+"/plots/extrap_0_p2_"+t.tag+"_"+t.ind.descr(iout)+".xmg");
      plot.write_polygon(bind(poly_eval<djvec_t>,coeffs,_1),0.0,*max_element(x.begin(),x.end())*1.1);
      write_fit_plot(plot,xmin,xmax,bind(poly_eval<djvec_t>,coeffs,_1),x,y);
      plot.write_ave_err(0.0,coeffs[0].ave_err());
    }
}

#endif
