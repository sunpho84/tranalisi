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
  
  /////////////////////////////////////////////////////////////////
  
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
  
  index_t i_in_clust_ihit_ind;
  index_t conf_ind;
  
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
  
  /////////////////////////////////////////////////////////////////
  
  //! prepares the list of configurations to use
  void prepare_list_of_confs();
  
  //! try to read, otherwise compute
  perens_t& read_or_compute();
  
  //! computes the basic Z
  perens_t& compute_basic(const string& ingredients_path);
};

#endif
