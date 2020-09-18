#pragma once
#include"Common\TriangleInterface.h"
#include"Common\Common.h"
#include"FindIndex.h"
#include"InnerScafData.h"
#include"SharedData.h"

class HybridMap
{
public:
	HybridMap(SharedData& sharedd_);
	~HybridMap();

public:
	void prepare_data();

	void tutte_map();
	void run_1_iter(int iter_num=-1);
	void autorun(double conv_rate_=1e-4);
	void reTriangulation();
	void setOrder();

	void update_pos_from_sD();

	void find_index_fine2coarse();
	void gen_s_mesh(const Eigen::MatrixXd& V_, const Eigen::MatrixXi& F_);
	void setNodes();
	void setIdentity();
	void assembly_CM_para();
	void prepare_pardiso();
	void solve_equation();
	void calc_src_decrease_direction(const vector<double>& c_d_, Eigen::MatrixXd& p_d_);
	void calc_src_gradient(const vector<double>& c_g_, Eigen::MatrixXd& p_g_);

	void assembly_SLIM_para();
	void assembly_scaf_part_SLIM_para();

	double max_step(const Eigen::MatrixXd& p_d_);
	void backtracking_line_search(const Eigen::MatrixXd& p_d_, const Eigen::MatrixXd& neg_grad, double & step_);
	double calc_energy(const Eigen::MatrixXd& pos_cur,bool with_scaf=false);
	Eigen::Vector3d calc_img(const Eigen::Vector3d& kesi, const std::vector<std::pair<int, BezierBase>>& bases_);
	void detect_inner_engine_seeds();

	void find_outer_rings(set<int>& part, std::vector<tuple<std::vector<int>, bool, int>>& boundaryloop);
	void set_fix_src();

	bool judge_loop_direction(const vector<int>& loop_);
	void engine_cond_num(set<int>& seeds_,int seeds_num);
	void set_minarea_for_retri(double minarea_);
	void engine_convergence_points(const vector<pair<double, int>>&e_dis_iso, set<int>& seeds_, int loops_num);
	void load_2d_mesh(const char* filename_);
	void leave_this_solver();
	void enter_this_solver();

	void info_src_hessian();
	pair<int,int> info_sim_hessian();
	
	//scaffold
	void add_scaffold();
	void meshid_vs_varid();
	void assembly_scaf_part_para();
	void calc_scaf_jacobian();

	//test

private:
	void calc_bc(Eigen::Vector3d& bc_, const int& f_id, const vector<double>& areas);

	template<typename T>
	void out_vec_file(const vector<T>& vec_, const string& file_);

public:
	pair<double, double> e_order_non1;
	pair<double, double> e_order_1;
	int order_set1 = 0;
	SharedData& sD;
	//Mesh mesh_whole;
	int V_N;
	int F_N;
	//vector<double> src_weight;
	std::vector<tuple<std::vector<int>, bool, int>> boundary_loops_part;//clockwise

	set<int> order1_fs;


	//update_src

	bool is_slim_meth=true;

	//simply-mesh after retriangulation
	MyMesh s_mesh;
	Eigen::MatrixXd s_V;
	Eigen::MatrixXi s_F;
	//Eigen::VectorXd pos_of_Smesh;

	//the fixed part of the simply-mesh
	//fix_ mean boundary 3-rings ,which are consisted with origin mesh
	//Eigen::MatrixXi fix_F;
	int fix_F_N;
	int fix_V_N;
	Eigen::VectorXi fix_v_src2s;// size: V_N, others -1;
	//Eigen::VectorXi fix_b_v_src2s;// size: V_N, others -1;
	Eigen::VectorXi fix_v_s2src;// size: fix_V.rows();
	Eigen::VectorXi fix_f_src2s;
	

	vector<pair<int,Eigen::Vector2d>> fine2coarse_v_id_pos;

	vector<vector<int>> samples_fid_sim;
	vector<vector<Eigen::Vector2d>> samples_bc;

	int sample_num_per_tri=4;
	vector<int> flag_order;

	//solver related
	Eigen::VectorXd variables_;
	int VAR_N;
	Eigen::SparseMatrix<int> find_idx_in_hess;
	PardisoSolver* pardiso_solver=NULL;
	int MAX_ITER_NUM = 500;
	double convergence_rate = 1e-4;
	double min_area_for_retri;

	InnerScafData scafd;
	int order1_fix_F_N;
	int order1_fix_V_N;
	vector<int> mid2varid;//smid2varid
	vector<int> varid2mid;//varid2smid

	int pAs;//pos_of_mesh And scaf
	Eigen::MatrixXi scaf_F;

	vector<pair<double,int>> steps_;

	//======================================
	set<int> engine_seeds_v;//

	pair<int,int> block_isscaf_fid;
	bool is_init_2d=false;

};

template<typename T>
inline void HybridMap::out_vec_file(const vector<T>& vec_, const string & file_)
{
	ofstream ofs(file_, ios::trunc);
	for (const T &var : vec_)
		ofs << var << endl;
	ofs.close();
}
