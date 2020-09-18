
#ifndef SCAFFOLD_TEST_SCAFDATA_H
#define SCAFFOLD_TEST_SCAFDATA_H

#include <Eigen/Dense>
#include <map>
#include"Common/Common.h"
#include"Common/TriangleInterface.h"
struct ScafData {
  // dimension for domain of parameterization/deformation
 public:
  ScafData() {
	  dim = 2;
	  mv_num = 0;
	  mf_num = 0;
	  sf_num = 0;
	  sv_num = 0;
	  mesh_measure = 0;
	  w_uv.resize(0, 0);
  };

  void add_new_patch(const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                     const Eigen::RowVectorXd &center);

  void mesh_improve(bool);

  void update_scaffold();

  double scaffold_factor = 10;

// Optional Input

  double exp_factor = 1.0; // used for exponential energies, ignored otherwise

// Output
  double energy; // objective value

// INTERNAL
  long mv_num, mf_num;
  long sv_num, sf_num;
  Eigen::MatrixXd m_V; // input initial mesh V
  Eigen::MatrixXi m_T; // input initial mesh F/T

  Eigen::MatrixXd w_uv; // whole domain uv: mesh + free vertices
  Eigen::MatrixXi s_T; // scaffold domain tets: scaffold tets
  Eigen::MatrixXi w_T;

  Eigen::VectorXd m_M; // mesh area or volume
  Eigen::VectorXd s_M; // scaffold area or volume
  Eigen::VectorXd w_M; // area/volume weights for whole
  double mesh_measure; // area or volume
  double avg_edge_length;
  long v_num;
  long f_num;
  double proximal_p = 1e-8; //unused

  Eigen::VectorXi frame_ids;
  double interval;
 public: // public for ser
  // caching
  Eigen::VectorXi internal_bnd;
  Eigen::MatrixXd rect_frame_V;

  // multi-chart support
  std::vector<int> component_sizes;
  std::vector<int> bnd_sizes;

  //3D
 public:
  Eigen::MatrixXi surface_F;

  int dim=2; // dimension for ambient space. Same for mesh/scaf
  
  };


  #endif //SCAFFOLD_TEST_SCAFDATA_H
