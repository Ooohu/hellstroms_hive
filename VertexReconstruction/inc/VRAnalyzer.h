

#ifndef VRANALYZER_H
#define VRANALYZER_H

#include "Analyzer.h"
#include "Storage.h"
#include "TTree.h"
#include "ParticleAssociations.h"
#include "VertexBuilder.h"
#include "VertexQuality.h"
#include "FillTreeVariables.h"


class VRAnalyzer : public Analyzer {

 public:

  VRAnalyzer(std::string const & name = "vrana", VertexQuality * vq = nullptr);
  ~VRAnalyzer();

  void SetVerbose(bool const verbose = true);
  void SetFlashTimings(double, double);
  void SetProducers(std::string const & track_producer,
		    std::string const & shower_producer,
		    std::string const & hit_producer,
		    std::string const & opflash_producer,
		    std::string const & trigger_product);
  void SetVBVariables(double const start_prox,
		      double const shower_prox,
		      double const max_bp_dist,
		      double const cpoa_vert_prox,
		      double const cpoa_trackend_prox);
  void RunPandora(bool const run_pandora = true);

  void GetPOT(bool const get_pot = true);
  void RunFillTreeVariables(bool const run_fill_tree_variables = true);
  void RunVerticesPerEvent();

  void Initialize();
  bool Run();
  void Finalize();

 private:

  void AddTracks(DetectorObjects & detos,
		 std::string const & producer,
		 bool const track_original_indices = false);
  void AddShowers(DetectorObjects & detos,
		  std::string const & producer,
		  bool const track_original_indices = false);
  geoalgo::Point_t GetSecondaryTrackPandoraVertex(ParticleAssociations const & pas,
						  size_t const track_index,
						  geoalgo::Point_t const & vertex);
  geoalgo::Point_t GetSecondaryTrackPandoraVertex(ParticleAssociations const & pas,
						  std::vector<size_t> const & associated_indices,
						  geoalgo::Point_t const & vertex);
  void FillLoneShowers(ParticleAssociations & pas,
		       std::vector<size_t> const & associated_showers);
  void FillWPandoraTrackRec(ParticleAssociations & pas,
			    size_t const pfp_index,
			    geoalgo::Point_t const & vertex);
  void PrintPandoraRec();
  void PrintPandoraRec(size_t const pfp_index,
		       std::string indent = "");
  void FillWPandora(ParticleAssociations & pas);
  void FillPOTTree();
  void VerticesPerEvent(ParticleAssociations const & pas);

  std::string fname;

  bool fverbose;

  bool fmcordata;
  std::string ftrack_producer;
  std::string fswtrigger_product;
  std::string fshower_producer;
  std::string fhit_producer;
  std::string fopflash_producer;

  double fflash_start;
  double fflash_end;

  double fstart_prox;
  double fshower_prox;
  double fmax_bp_dist;
  double fcpoa_vert_prox;
  double fcpoa_trackend_prox;

  bool frun_pandora;
  bool fget_pot;
  bool frun_fill_tree_variables;

  VertexQuality * fvq;
  FillTreeVariables fftv;

  TTree * fvertices_per_event_tree;
  int frun_number;
  int fsubrun_number;
  int fevent_number;
  int fvertex_number;

};


#endif
