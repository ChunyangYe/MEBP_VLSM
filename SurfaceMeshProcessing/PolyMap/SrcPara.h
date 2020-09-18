#pragma once
#include"Common\Common.h"
#include"SharedData.h"
#include <boost/property_tree/ptree.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/coarsening/rigid_body_modes.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/preconditioner/runtime.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/coarsening/ruge_stuben.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/iluk.hpp>
#include <amgcl/relaxation/ilut.hpp>
#include <amgcl/solver/cg.hpp>
#include<amgcl/solver/bicgstab.hpp>
#include<amgcl/solver/lgmres.hpp>
#include <amgcl/profiler.hpp>
#include<amgcl/io/mm.hpp>
class SrcPara
{
public:
	SrcPara(SharedData& shardd_);
	~SrcPara();

private:
	void prepare_pardiso_with_fixed_boundary();
	void assembly_pardiso_with_fixed_boundary();
	void assembly_pardiso_scaffold();
	double max_step_general(const Eigen::MatrixXd& decrease_);
	void backtracking_general(const Eigen::MatrixXd& neg_grad_, const Eigen::MatrixXd& decrease_, double& step_);
	double calc_energy(const Eigen::MatrixXd& cur_pos,bool is_add_scaf=true);
	void global2local_with_fixed_boundary();
	void init();
	void solve_equation();
	void graph_color_GS(const std::vector<std::set<int>>& VVs);
	void GS_solver();
	void amgcl_solver();
	void amg_write();
	void assemble_hessian_gradient_slim();
public:
	void run();
	void set_method(int meth_type);

private:
	SharedData &sD;

	bool is_amgcl_solver = true;

	int v_w_N;
	int v_free_N;

	//PardisoSolver ps_local;
	PardisoSolver* ps_local = NULL;

	Eigen::SparseMatrix<int> SM_find_id;
	std::vector<std::vector<int>> phase;
	vector<double> result_gs;

	bool is_pardiso_solver = false;
	vector<int> var2src;
	vector<int> src2var;

	int MAX_ITERATION_N = 1000;
	double CONVERGENCE_RATE = 1e-4;

};

