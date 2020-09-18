#pragma once
#include"Common\Common.h"
#include"SharedData.h"
class LocalPara
{
public:
	LocalPara(SharedData& shardd_, const set<int>& fids_p_, vector<int>& f_visit_);
	~LocalPara();

private:
	void prepare_pardiso_with_fixed_boundary();
	void assembly_pardiso_with_fixed_boundary();
	double max_step_general(const Eigen::MatrixXd& decrease_);
	void backtracking_general(const Eigen::MatrixXd& neg_grad_, const Eigen::MatrixXd& decrease_, double& step_);
	double calc_energy(const Eigen::MatrixXd& cur_pos);
	void global2local_with_fixed_boundary();
	void init();
	void solve_equation();
	void update_src_pos();

public:
	void run();

private:
	SharedData &sD;
	vector<int>& f_visit;

	vector<int> fids_p;

	int v_src_N;
	int v_w_N;
	int v_free_N;
	double weight_total;

	Eigen::MatrixXi F_src_p;
	Eigen::MatrixXi F_p;
	Eigen::MatrixXd V_p;

	vector<double> weight_p;
	PardisoSolver ps_local;
	Eigen::SparseMatrix<int> SM_find_id;
	std::vector<Eigen::Vector4d> src_inv_p;
	Eigen::VectorXi v_l2g, v_g2l;

	int MAX_ITERATION_N = 5;
	double CONVERGENCE_RATE = 1e-2;

};