#pragma once
#include "Common\Common.h"
#include "PolyMap\SharedData.h"

#include "BSplineCurve.h"
#include "BSplineSurface.h"

using namespace Eigen;
class Parafun
{
public:
	Parafun(SharedData & data);
	~Parafun();

	void init(int grid_num_=50);
	void Pre_calculate_spline();

	void after_mesh_improve();

	void spline_init();
	void cntrl_vs_index();

	double perform_1iteration_spline(bool use_CM);
	void CM_spline();
	void SLIM_spline();
	void BPE_spline(bool is_slim=true);
	void iteratively_BPE_spline();

	void backtracking_line_search(const Eigen::MatrixXd &x, const Eigen::MatrixXd &d, const double &neg_grad_dot_d, double &alpha);
	void max_step(const Eigen::MatrixXd &x, const Eigen::MatrixXd &d, double &step);


	void calc_gradient_norm(Eigen::MatrixXd& negative_grad_norm);
	void Reconstruct(Eigen::MatrixXd& direction);

	double compute_energy(const Eigen::MatrixXd & x, bool whole = false);


	void leave_this_solver();
	void enter_this_solver();
	void clear_samples();

	//data 
	SharedData& d_;
	int set_ctrl_num = 1200;

	int total_num;
	int V_N;
	int F_N;

	BSplineSurface	*spline_surface;
	vector<double>	u_x_m;
	vector<double>	u_y_m;
	vector<double>	v_x_m;
	vector<double>	v_y_m;
	double			max_energy=25000;

	size_t			Control_N;
	size_t			total_cntrl;
	size_t			Sample_per_face = 4;
	size_t			Sample_per_scaffold = 4;
	vector<double>	samplepointsx, samplepointsy;
	size_t			Sample_N;
	size_t			Sample_N_m;
	vector<int>		sample_at_face;

	VectorXd		pos_of_control;

	vector<int>		var_cntrls;
	vector<int>		cntrl2index;

	PardisoSolver* spline_pardiso;
	SparseMatrix<int> find_cntrl_in_rows;

	//Assist
	int Assist_N;
	int Interval_N;

	//information
	double time_consumption;

	int	iter_general_cout = 0;

	double convgence_con_rate = 1e-4;
	int MAX_ITER_NUM = 500;

};

